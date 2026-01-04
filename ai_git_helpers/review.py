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

MAX_FILE_BYTES = 200_000
MAX_DIFF_BYTES_TOTAL = 400_000
MODEL = os.getenv("AI_REVIEW_MODEL", "gpt-4o-mini")
MAX_OUTPUT_TOKENS = 800
TEMPERATURE = 0.2

SYSTEM_PROMPT = """\
You are an expert C++ software engineer with deep knowledge of
performance, correctness, and undefined behavior. You also have
knowledge of Python and R.

You will receive:
1) The full source file (for context only)
2) A git diff showing recent changes

Focus primarily on the diff.
Use the full file only to understand context and invariants.

When only FULL FILE CONTEXT is provided, treat this as an
intentional holistic code review request. Do not mention
or warn about missing diffs.

Priorities:
1. Undefined behavior / correctness issues
2. Performance implications
3. Readability / maintainability
4. API consistency
5. Only use C++17 and below
6. Python files are personal developer tools, not production libraries

Be concise and practical.

For Python files, treat them as personal developer tools rather than
widely distributed production libraries. Prefer pragmatic, readable
solutions over enterprise-level robustness or defensive overengineering.

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


def run_cmd_limited(cmd: list[str], max_bytes: int, timeout_s: int = 10) -> str:
    p = None
    try:
        p = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
        )
        out, _ = p.communicate(timeout=timeout_s)

        truncated = len(out) > max_bytes
        out = out[:max_bytes]

        text = out.decode("utf-8", errors="replace")
        if truncated:
            text += "\n\n[...output truncated...]\n"
        return text

    except subprocess.TimeoutExpired:
        if p is not None:
            p.kill()
            p.communicate()
        return f"[Command timed out: {' '.join(cmd)}]"

    except Exception as e:
        if p is not None:
            try:
                p.kill()
                p.communicate()
            except Exception:
                pass
        return f"[Error running command: {e}]"


def git_diff_file(path: str, ref: str, max_bytes: int = MAX_FILE_BYTES) -> str:
    return run_cmd_limited(
        ["git", "diff", "--no-color", "--no-ext-diff", "-U0", ref, "--", path],
        max_bytes=max_bytes,
        timeout_s=10,
    )


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


def read_file_safely(path: str, max_bytes: int=200_000) -> str:
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
        "Default: diff-only. Use --with-context to include full files.\n"
        "If no diff is provided, review the file(s) holistically."
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

    parser.add_argument(
        "--file",
        action="append",
        default=[],
        help="Review this file (can be provided multiple times).",
    )

    parser.add_argument(
        "--diff-ref",
        metavar="REF",
        help="With --file, include git diff against this ref (e.g. HEAD, main).",
    )

    args = parser.parse_args()

    if args.diff_ref and not args.file:
        print(
            "Use: ai-review --file <path> --diff-ref REF",
            file=sys.stderr
        )
        return 1

    stdin_text = ""
    if not sys.stdin.isatty():
        stdin_text = sys.stdin.read()

    use_file_mode = bool(args.file)
    use_stdin_mode = bool(stdin_text.strip()) and not use_file_mode

    if not (use_stdin_mode or use_file_mode):
        print("No input provided.", file=sys.stderr)
        return 1

    if use_file_mode:
        print(
            f">>> stdin chars (ignored due to --file): {len(stdin_text)}",
            file=sys.stderr
        )
    else:
        print(f">>> diff chars: {len(stdin_text)}", file=sys.stderr)

    file_contexts = []
    files_for_context = []
    repo_root = get_repo_root()
    user_content = ""

    if use_file_mode:
        for rel_path in dict.fromkeys(args.file):
            if is_denied(rel_path):
                continue

            abs_path = safe_resolve_path(repo_root, rel_path)
            if abs_path is None:
                continue

            rel_path_norm = abs_path.relative_to(repo_root).as_posix()
            file_contexts.append(
                f"\n===== FILE: {rel_path_norm} =====\n{read_file_safely(str(abs_path))}"
            )
            files_for_context.append(rel_path_norm)

        diff_text = ""

        if args.diff_ref:
            # Include diffs for focus
            diffs = []
            remaining = MAX_DIFF_BYTES_TOTAL

            for rel_path_norm in files_for_context:
                if remaining <= 0:
                    diffs.append("\n[...additional diffs omitted due to size limit...]\n")
                    break

                d = git_diff_file(
                    rel_path_norm, args.diff_ref,
                    max_bytes=min(remaining, MAX_FILE_BYTES)
                )

                # Budget is approximate (post-decode + truncation markers),
                # but keeps requests bounded.
                remaining -= len(d.encode("utf-8", errors="replace"))

                if d.strip():
                    diffs.append(d)

            diff_text = "\n".join(diffs)

        if diff_text.strip():
            user_content = (
                "FULL FILE CONTEXT:\n"
                + "\n".join(file_contexts)
                + "\n\nDIFF (for focus):\n"
                + diff_text
            )
        else:
            user_content = "FULL FILE CONTEXT:\n" + "\n".join(file_contexts)

    else:
        # Collect file paths from diff
        files = collect_paths_from_diff(stdin_text)
        use_context = args.with_context and not args.diff_only

        # Only include full file context if explicitly requested
        if use_context:
            print(
                ">>> WARNING: --with-context will upload full file contents to the API. "
                "Use with care.",
                file=sys.stderr,
            )

            for rel_path in sorted(files):
                if is_denied(rel_path):
                    continue

                abs_path = safe_resolve_path(repo_root, rel_path)
                if abs_path is None:
                    continue

                rel_path_norm = abs_path.relative_to(repo_root).as_posix()
                file_contexts.append(
                    f"\n===== FILE: {rel_path_norm} =====\n"
                    f"{read_file_safely(str(abs_path))}"
                )
                files_for_context.append(rel_path_norm)

            unique_files = sorted(set(files_for_context))
            count = len(unique_files)

            print(
                f">>> Including context for {count} file"
                f"{'' if count == 1 else 's'}:",
                file=sys.stderr,
            )

            for f in unique_files:
                print(f">>>   {f}", file=sys.stderr)

        if use_context:
            full_context = "\n".join(file_contexts)
            print(f">>> full_context chars: {len(full_context)}", file=sys.stderr)
            user_content = f"FULL FILE CONTEXT:\n{full_context}\n\nDIFF:\n{stdin_text}"
        else:
            user_content = f"DIFF:\n{stdin_text}"


    # File mode: explicit, precise failure
    if use_file_mode and not file_contexts:
        print(
            "No readable files provided (denied/missing/outside repo).",
            file=sys.stderr,
        )
        return 1

    # Generic final guard (covers all paths)
    if not user_content.strip():
        print(
            "No input provided. Pipe a diff or pass --file <path>.",
            file=sys.stderr,
        )
        return 1

    api_key = os.getenv("RCPPALGOS_OPENAI_KEY")

    if not api_key:
        print("Missing API key.", file=sys.stderr)
        return 1

    client = OpenAI(api_key=api_key)

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

