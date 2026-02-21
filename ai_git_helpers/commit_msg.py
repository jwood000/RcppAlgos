#!/usr/bin/env python3

import os
import sys
import re
import argparse
import textwrap
from openai import OpenAI

_FENCE_RE = re.compile(r"^\s*```")
_SUBJECT_RE = re.compile(
    r"^(feat|fix|docs|style|refactor|perf|test|chore|build|ci)"
    r"(\([^)]+\))?:\s+\S",
    re.IGNORECASE,
)

_BULLET_RE = re.compile(r'^(\s*)([-*]|\d+\.)\s+')

ALLOWED_PREFIXES = {
    "feat",
    "fix",
    "docs",
    "style",
    "refactor",
    "perf",
    "test",
    "chore",
    "build",
    "ci"
}

MODEL = os.getenv("AI_REVIEW_MODEL", "gpt-4o-mini")
MAX_OUTPUT_TOKENS = 200
TEMPERATURE = 0.2
SUBJECT_MAX = 72
BODY_WIDTH = 72


# -----------------------------
# Formatting helpers
# -----------------------------

def _strip_code_fences(raw: str) -> str:
    """Remove surrounding markdown code fences if present."""
    s = (raw or "").strip()
    lines = s.splitlines()
    if lines and _FENCE_RE.match(lines[0].strip()):
        lines = lines[1:]
        while lines and _FENCE_RE.match(lines[-1].strip()):
            lines = lines[:-1]
        s = "\n".join(lines).strip()
    return s


def _first_subject_and_rest(raw: str):
    """
    Return (subject_line, body_text) by choosing the first line that looks like a
    conventional-commit subject. If none match, fall back to first meaningful line.
    """
    lines = [ln.rstrip() for ln in (raw or "").splitlines()]
    cleaned = [ln.strip() for ln in lines if ln.strip() and ln.strip() != "```"]

    if not cleaned:
        return "", ""

    subj_idx = None
    for i, ln in enumerate(cleaned):
        if _SUBJECT_RE.match(ln):
            subj_idx = i
            break

    if subj_idx is None:
        subject = cleaned[0]
        body = "\n".join(cleaned[1:]).strip()
        return subject, body

    subject = cleaned[subj_idx]
    body = "\n".join(cleaned[subj_idx + 1:]).strip()
    return subject, body


def _ensure_prefix(subject: str, default_prefix: str = "chore") -> str:
    """
    Ensure the subject starts with a Conventional Commits prefix.
    """
    if not subject:
        return f"{default_prefix}:"

    head, sep, tail = subject.partition(":")
    if sep and head.strip().lower() in ALLOWED_PREFIXES:
        return subject

    return f"{default_prefix}: {subject}"


def _truncate(text: str, max_len: int) -> str:
    if len(text) <= max_len:
        return text

    truncated = text[:max_len]
    last_space = truncated.rfind(" ")
    if last_space > 0:
        truncated = truncated[:last_space]

    return truncated.rstrip()


def _format_commit_message(raw: str,
                           subject_max: int = SUBJECT_MAX,
                           body_width: int = BODY_WIDTH) -> str:

    raw = _strip_code_fences(raw)
    raw = (raw or "").strip()
    if not raw:
        return ""

    subject, body_text = _first_subject_and_rest(raw)
    if not subject:
        return ""

    subject = re.sub(r"\s+", " ", subject).strip()
    subject = subject.rstrip(".")

    subject = _ensure_prefix(subject)
    subject = _truncate(subject, subject_max)

    if not body_text:
        return subject + "\n"

    parts = []
    paragraphs = re.split(r"\n\s*\n", body_text)

    for para in paragraphs:
        para = para.rstrip()
        if not para.strip():
            continue

        para_lines = [ln.rstrip() for ln in para.splitlines() if ln.strip()]

        if para_lines and all(_BULLET_RE.match(ln) for ln in para_lines):
            wrapped_bullets = []
            for ln in para_lines:
                m = _BULLET_RE.match(ln)
                indent = m.group(1)
                marker = m.group(2)
                prefix = f"{indent}{marker} "
                content = ln[m.end():].strip()

                wrapped_bullets.append(textwrap.fill(
                    content,
                    width=body_width,
                    initial_indent=prefix,
                    subsequent_indent=" " * len(prefix),
                    break_long_words=False,
                    break_on_hyphens=False,
                ))
            parts.append("\n".join(wrapped_bullets))
        else:
            collapsed = re.sub(r"\s+", " ", " ".join(para_lines)).strip()
            parts.append(textwrap.fill(
                collapsed,
                width=body_width,
                break_long_words=False,
                break_on_hyphens=False,
            ))

    return subject + "\n\n" + "\n\n".join(parts).rstrip() + "\n"


# -----------------------------
# Prompts
# -----------------------------

SYSTEM_PROMPT = """\
You are an expert C++ and R developer.
You write concise, meaningful, high-quality git commit messages.

Use the Conventional Commits style.

Allowed types:
- fix: when behavior changes to correct incorrect output
- refactor: when restructuring without behavioral change
- perf: when improving time or memory efficiency
- feat: when adding new externally visible functionality
- docs: documentation only
- test: adding or fixing tests
- chore: maintenance, tooling, formatting, or cleanup
- build: build system or dependency changes
- ci: CI configuration changes

Rules:
- Subject format: <type>: <short description>
- Use imperative mood (e.g., "fix", not "fixed")
- Subject line max 72 characters
- No trailing period in subject
- Choose the most accurate type
- Only include a body if it adds real value
- Do not mention tooling, AI, or ChatGPT
- Do not claim changes that are not supported by the diff.
  - THEME/CONTEXT describe intent or motivation.
  - Use them to guide framing and "why", not to invent "what".
  - If intent mentions work not shown, phrase it as support or motivation.

Format:
<type>: <subject>

<body (optional)>
"""

USER_PROMPT_TEMPLATE = """\
Write a Conventional Commit message for the following staged diff.

Choose the most appropriate type and include a short body if it adds value.

THEME (intent / goal; may be empty):
{theme}

ADDITIONAL CONTEXT (may be empty):
{context}

Diff:
{diff}
"""


# -----------------------------
# Main
# -----------------------------

def _read_text_file(path: str) -> str:
    with open(path, "r", encoding="utf-8") as f:
        return f.read().strip()


def main():
    parser = argparse.ArgumentParser(
        description="Generate a Conventional Commit message from a staged diff."
    )
    parser.add_argument(
        "--theme", "--intent",
        dest="theme",
        default="",
        help="High-level intent or goal for this commit (guides framing only).",
    )
    parser.add_argument(
        "--context-file",
        dest="context_file",
        default="",
        help="Path to a text file with additional context (optional).",
    )
    parser.add_argument(
        "--ticket",
        dest="ticket",
        default="",
        help="Ticket or issue id to include in context (optional).",
    )

    args = parser.parse_args()
    diff = sys.stdin.read()

    if not diff.strip():
        print("No diff provided on stdin.", file=sys.stderr)
        sys.exit(1)

    api_key = os.getenv("RCPPALGOS_OPENAI_KEY")
    if not api_key:
        print("RCPPALGOS_OPENAI_KEY is not set.", file=sys.stderr)
        sys.exit(1)

    theme = (args.theme or "").strip()
    context_chunks = []

    if args.ticket.strip():
        context_chunks.append(f"Ticket: {args.ticket.strip()}")

    if args.context_file.strip():
        try:
            context_chunks.append(_read_text_file(args.context_file.strip()))
        except OSError as e:
            print(f"Could not read --context-file: {e}", file=sys.stderr)
            sys.exit(1)

    context = "\n\n".join([c for c in context_chunks if c]).strip()

    client = OpenAI(api_key=api_key)

    response = client.responses.create(
        model=MODEL,
        input=[
            {"role": "system", "content": SYSTEM_PROMPT},
            {"role": "user", "content": USER_PROMPT_TEMPLATE.format(
                diff=diff,
                theme=theme,
                context=context
            )},
        ],
        max_output_tokens=MAX_OUTPUT_TOKENS,
        temperature=TEMPERATURE,
    )

    message = (response.output_text or "").strip()

    if not message:
        print(">>> Empty output_text from model.", file=sys.stderr)
        print(">>> Raw response id:", getattr(response, "id", None), file=sys.stderr)
    else:
        print(_format_commit_message(message))


if __name__ == "__main__":
    main()
