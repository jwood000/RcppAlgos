#!/usr/bin/env python3

import os
import sys
import re
import textwrap
from openai import OpenAI

# -----------------------------
# Configuration
# -----------------------------

MODEL = "gpt-5.2"
MAX_OUTPUT_TOKENS = 200

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

ALLOWED_PREFIXES = (
    "fix:", "feat:", "refactor:", "perf:",
    "docs:", "test:", "chore:", "build:", "ci:"
)

def _ensure_prefix(subject: str) -> str:

    if any(subject.startswith(p) for p in ALLOWED_PREFIXES):
        return subject

    return f"chore: {subject}"


_BULLET_RE = re.compile(r'^(\s*)([-*]|\d+\.)\s+')


def _format_commit_message(raw: str,
                           subject_max: int = SUBJECT_MAX,
                           body_width: int = BODY_WIDTH) -> str:
    raw = raw.strip()
    if not raw:
        return ""

    lines = raw.splitlines()

    # Subject = first non-empty line
    i = 0
    while i < len(lines) and not lines[i].strip():
        i += 1
    subject = lines[i].strip() if i < len(lines) else ""
    i += 1

    # Rest = body
    body_text = "\n".join(lines[i:]).strip()

    # Normalize subject spacing, remove trailing period, enforce max length
    subject = re.sub(r"\s+", " ", subject).strip()
    subject = subject.rstrip(".")

    if len(subject) > subject_max:
        subject = subject[:subject_max].rstrip()

    subject = _ensure_prefix(subject)

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

    response = client.chat.completions.create(
        model=MODEL,
        messages=[
            {"role": "system", "content": SYSTEM_PROMPT},
            {"role": "user", "content": USER_PROMPT_TEMPLATE.format(diff=diff)},
        ],
        max_completion_tokens=MAX_OUTPUT_TOKENS,
        temperature=0.2,
    )

    message = response.choices[0].message.content or ""
    print(_format_commit_message(message))


if __name__ == "__main__":
    main()
