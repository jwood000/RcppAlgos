#!/usr/bin/env python3

import os
import re
import sys
import fnmatch
import argparse
import subprocess
from pathlib import Path
from openai import OpenAI
from typing import Optional


MODEL = os.getenv("AI_REVIEW_MODEL", "gpt-4o-mini")
MAX_OUTPUT_TOKENS = 800
TEMPERATURE = 0.2

SYSTEM_PROMPT = """\
You are an expert C++ software engineer with deep knowledge of
performance, correctness, and undefined behavior.

You will receive:
1) The full source file (for context only)
2) A git diff showing recent changes

Focus primarily on the diff.
Use the full file only to understand context and invariants.

Priorities:
1. Undefined behavior / correctness issues
2. Performance implications
3. Readability / maintainability
4. API consistency
5. Only use C++17 and below

Be concise and practical.

Output format:

# Summary
# Blockers
# Risks / correctness
# Performance
# Readability
# Tests / edge cases
"""

USE_COLOR = os.getenv("AI_REVIEW_COLOR", "1") != "0" and sys.stdout.isatty()

ANSI = {
    "reset": "\033[0m",
    "bold": "\033[1m",
    "dim": "\033[2m",
    "red": "\033[31m",
    "green": "\033[32m",
    "yellow": "\033[33m",
    "blue": "\033[34m",
    "magenta": "\033[35m",
    "cyan": "\033[36m",
    "underline": "\033[4m"
}


def color(s: str, *styles: str) -> str:
    if not USE_COLOR:
        return s

    prefix = "".join(ANSI.get(t, "") for t in styles)
    return f"{prefix}{s}{ANSI['reset']}"


def md_inline_to_ansi(s: str) -> str:
    # Inline code: `code`
    def repl_code(m):
        inner = m.group(1)
        return color(inner, "cyan")  # code-ish

    # Bold: **text**
    def repl_bold(m):
        inner = m.group(1)
        return color(inner, "bold")

    # Order matters: code first, then bold.
    s = re.sub(r"`([^`]+)`", repl_code, s)
    s = re.sub(r"\*\*([^*]+)\*\*", repl_bold, s)
    return s


def render_markdown(text: str) -> str:
    lines = text.splitlines()
    out = []
    for line in lines:
        # Headings
        if line.startswith("### "):
            out.append(color("─" * 60, "dim"))
            out.append(color(md_inline_to_ansi(line[4:]), "bold", "magenta"))
            continue
        if line.startswith("## "):
            out.append(color("─" * 60, "dim"))
            out.append(color(md_inline_to_ansi(line[3:]), "bold", "blue"))
            continue
        if line.startswith("# "):
            out.append(color("─" * 60, "dim"))
            title = md_inline_to_ansi(line[2:])
            # Give Blockers a stronger color if present
            if "Blockers" in line:
                out.append(color(title, "bold", "red"))
            else:
                out.append(color(title, "bold", "green"))
            continue

        # Bullets
        stripped = line.lstrip()
        indent = len(line) - len(stripped)

        if stripped.startswith("- "):
            bullet = color("•", "dim")
            content = md_inline_to_ansi(stripped[2:])
            out.append((" " * indent) + bullet + " " + content)
            continue

        if stripped.startswith("* "):
            bullet = color("•", "dim")
            content = md_inline_to_ansi(stripped[2:])
            out.append((" " * indent) + bullet + " " + content)
            continue

        # Default: just do inline markdown rendering
        out.append(md_inline_to_ansi(line))

    return "\n".join(out)


def get_repo_root() -> Path:
    try:
        out = subprocess.check_output(
            ["git", "rev-parse", "--show-toplevel"],
            stderr=subprocess.DEVNULL,
            text=True,
        ).strip()
        return Path(out)
    except Exception:
        # Fall back to current working dir if not in a git repo
        return Path.cwd()


def safe_resolve_path(repo_root: Path, rel_path: str) -> Optional[Path]:
    # Normalize diff paths like "a/foo" or "b/foo" already stripped earlier
    candidate = (repo_root / rel_path).resolve()

    # Python 3.9: use commonpath
    root = str(repo_root.resolve())
    cand = str(candidate)
    if os.path.commonpath([root, cand]) != root:
        return None
    return candidate


DENY_GLOBS = [
    ".env", "*.env",
    ".Renviron", ".Rprofile",
    "*.pem", "*.key", "*.p12", "*.crt",
    "id_rsa*", "id_ed25519*",
    "*secret*", "*token*", "*credential*",
    ".git/*", ".git/**",
]

def is_denied(rel_path: str) -> bool:
    # match both basename and full path globs
    base = os.path.basename(rel_path)
    for pat in DENY_GLOBS:
        if fnmatch.fnmatch(rel_path, pat) or fnmatch.fnmatch(base, pat):
            return True
    return False


def read_file_safely(path: str, max_bytes: int = 200_000) -> str:
    try:
        p = Path(path)
        if not p.exists() or not p.is_file():
            return f"[Could not read file: {path}]"

        with p.open("rb") as f:
            data = f.read(max_bytes + 1)  # +1 so we can detect truncation

        truncated = len(data) > max_bytes
        if truncated:
            data = data[:max_bytes]

        text = data.decode("utf-8", errors="replace")
        if truncated:
            text += "\n\n[...file truncated...]\n"
        return text

    except Exception as e:
        return f"[Error reading file: {e}]"


def collect_paths_from_diff(diff: str) -> set[str]:
    files = set()

    for line in diff.splitlines():
        if line.startswith("diff --git "):
            # line format: diff --git a/foo/bar.cpp b/foo/bar.cpp
            parts = line.split()
            if len(parts) >= 4:
                a_path = parts[2]
                b_path = parts[3]
                if a_path.startswith("a/"):
                    a_path = a_path[2:]
                if b_path.startswith("b/"):
                    b_path = b_path[2:]
                # Prefer the "b/" path (new path), but keep "a/" too if different.
                if b_path != "/dev/null":
                    files.add(b_path)
                if a_path != "/dev/null":
                    files.add(a_path)

    return files


def main():
    parser = argparse.ArgumentParser(description="AI-assisted code review")

    parser.description = (
        "AI-assisted code review tool.\n"
        "Default: diff-only. Use --with-context to include full files."
    )

    parser.add_argument(
        "--with-context",
        action="store_true",
        help="Include full file contents for files in the diff (may be sensitive).",
    )

    parser.add_argument(
        "--diff-only",
        action="store_true",
        help="Review only the diff (default).",
    )

    args = parser.parse_args()

    diff = sys.stdin.read()
    print(f">>> diff chars: {len(diff)}", file=sys.stderr)

    if not diff.strip():
        print("No diff provided.", file=sys.stderr)
        return 1

    # Collect file paths from diff
    files = collect_paths_from_diff(diff)

    file_contexts = []
    use_context = args.with_context and not args.diff_only

    # Only include full file context if explicitly requested
    if use_context:
        print(
            ">>> WARNING: --with-context will upload full file contents to the API. "
            "Use with care.",
            file=sys.stderr,
        )

        repo_root = get_repo_root()
        files_for_context = set()

        for rel_path in sorted(files):
            if is_denied(rel_path):
                continue

            abs_path = safe_resolve_path(repo_root, rel_path)
            if abs_path is None:
                # Path traversal / outside-repo path; skip
                continue

            file_contexts.append(
                f"\n===== FILE: {rel_path} =====\n"
                f"{read_file_safely(str(abs_path))}"
            )

            files_for_context.add(rel_path)

        count = len(files_for_context)
        print(
            f">>> Including context for {count} file"
            f"{'' if count == 1 else 's'}:",
            file=sys.stderr,
        )

        for f in sorted(files_for_context):
            print(f">>>   {f}", file=sys.stderr)

    full_context = "\n".join(file_contexts)
    print(f">>> full_context chars: {len(full_context)}", file=sys.stderr)

    api_key = os.getenv("RCPPALGOS_OPENAI_KEY")

    if not api_key:
        print("Missing API key.", file=sys.stderr)
        return 1

    client = OpenAI(api_key=api_key)

    if use_context:
        user_content = f"FULL FILE CONTEXT:\n{full_context}\n\nDIFF:\n{diff}"
    else:
        user_content = f"DIFF:\n{diff}"

    response = client.responses.create(
        model=MODEL,
        input=[
            {"role": "system", "content": SYSTEM_PROMPT},
            {"role": "user", "content": user_content},
        ],
        max_output_tokens=MAX_OUTPUT_TOKENS,
        temperature=TEMPERATURE,
    )

    text = (response.output_text or "").strip()

    if not text:
        # Helpful debug if the model produced no text output
        print(">>> Empty output_text from model.", file=sys.stderr)
        print(">>> Raw response id:", getattr(response, "id", None), file=sys.stderr)
    else:
        print(render_markdown(text))


if __name__ == "__main__":
    raise SystemExit(main())

