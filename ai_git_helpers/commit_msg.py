#!/usr/bin/env python3

import os
import sys
import re
import textwrap
from openai import OpenAI

_FENCE_RE = re.compile(r"^\s*```")
_SUBJECT_RE = re.compile(
    r"^(feat|fix|docs|style|refactor|perf|test|chore|build|ci)"
    r"(\([^)]+\))?:\s+\S",
    re.IGNORECASE,
)


def _strip_code_fences(raw: str) -> str:
    """Remove surrounding markdown code fences if present."""
    s = (raw or "").strip()

    # If the whole output is fenced, remove first/last fence lines.
    lines = s.splitlines()
    if lines and _FENCE_RE.match(lines[0].strip()):
        lines = lines[1:]
        # Drop trailing fence if present (allow whitespace before it)
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
    # Remove empty lines and standalone fences
    cleaned = [ln.strip() for ln in lines if ln.strip() and ln.strip() != "```"]

    if not cleaned:
        return "", ""

    # Prefer the first line that *already* looks like a conventional subject
    subj_idx = None
    for i, ln in enumerate(cleaned):
        if _SUBJECT_RE.match(ln):
            subj_idx = i
            break

    if subj_idx is None:
        # Fallback: first meaningful line
        subject = cleaned[0]
        body = "\n".join(cleaned[1:]).strip()
        return subject, body

    subject = cleaned[subj_idx]
    body = "\n".join(cleaned[subj_idx + 1:]).strip()
    return subject, body


MODEL = os.getenv("AI_REVIEW_MODEL", "gpt-4o-mini")
MAX_OUTPUT_TOKENS = 200
TEMPERATURE = 0.2
SUBJECT_MAX = 72
BODY_WIDTH = 72

SYSTEM_PROMPT = """\
You are an expert C++ and R developer.
You write concise, high-quality git commit messages.

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

Format:
<type>: <subject>

<body (optional)>
"""

USER_PROMPT_TEMPLATE = """\
Write a Conventional Commit message for the following staged diff.

Choose the most appropriate type and include a short body if it adds value.

Diff:
{diff}
"""

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

def _ensure_prefix(subject: str, default_prefix: str = "chore") -> str:
    """
    Ensure the subject starts with a Conventional Commits prefix.
    Accepts existing prefixes case-insensitively, but does not rewrite case.
    """
    if not subject:
        return f"{default_prefix}:"

    # Extract candidate prefix before first colon
    head, sep, tail = subject.partition(":")
    if sep:
        if head.strip().lower() in ALLOWED_PREFIXES:
            # Valid prefix already present (Fix:, FIX:, fix:, etc.)
            return subject

    # No valid prefix → prepend default
    return f"{default_prefix}: {subject}"


def _truncate(text: str, max_len: int) -> str:
    """
    Truncate text to at most max_len characters.
    Prefer truncating at a word boundary. If no boundary is found,
    perform a hard truncation.
    """
    if len(text) <= max_len:
        return text

    # Cut to max_len first
    truncated = text[:max_len]

    # If there's a space, truncate back to the last full word
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

    # Normalize subject spacing, remove trailing period, enforce max length
    subject = re.sub(r"\s+", " ", subject).strip()
    subject = subject.rstrip(".")

    if len(subject) > subject_max:
        subject = subject[:subject_max].rstrip()

    subject = _ensure_prefix(subject)
    subject = _truncate(subject, subject_max)

    if not body_text:
        return subject + "\n"

    # Wrap body to body_width, preserving paragraphs and bullets.
    parts = []
    paragraphs = re.split(r"\n\s*\n", body_text)

    for para in paragraphs:
        para = para.rstrip()
        if not para.strip():
            continue

        para_lines = [ln.rstrip() for ln in para.splitlines() if ln.strip()]

        # Bullet paragraph: wrap each bullet with hanging indent
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
            # Normal paragraph: collapse internal whitespace then wrap
            collapsed = re.sub(r"\s+", " ", " ".join(para_lines)).strip()
            parts.append(textwrap.fill(
                collapsed,
                width=body_width,
                break_long_words=False,
                break_on_hyphens=False,
            ))

    return subject + "\n\n" + "\n\n".join(parts).rstrip() + "\n"

# -----------------------------
# Main
# -----------------------------

def main():
    diff = sys.stdin.read()

    if not diff.strip():
        print("No diff provided on stdin.", file=sys.stderr)
        sys.exit(1)

    api_key = os.getenv("RCPPALGOS_OPENAI_KEY")
    if not api_key:
        print("RCPPALGOS_OPENAI_KEY is not set.", file=sys.stderr)
        sys.exit(1)

    client = OpenAI(api_key=api_key)

    response = client.responses.create(
        model=MODEL,
        input=[
            {"role": "system", "content": SYSTEM_PROMPT},
            {"role": "user", "content": USER_PROMPT_TEMPLATE.format(diff=diff)},
        ],
        max_output_tokens=MAX_OUTPUT_TOKENS,
        temperature=TEMPERATURE,
    )

    message = (response.output_text or "").strip()

    if not message:
        # Helpful debug if the model produced no text output
        print(">>> Empty output_text from model.", file=sys.stderr)
        print(">>> Raw response id:", getattr(response, "id", None), file=sys.stderr)
    else:
        print(_format_commit_message(message))


if __name__ == "__main__":
    main()

