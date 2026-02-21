#!/usr/bin/env python3

import os
import sys
from openai import OpenAI

MODEL = "gpt-5.2"

SYSTEM_PROMPT = """
You are helping maintain the R package RcppAlgos. You are an expert C++
developer and understand the R ecosystem very well.

Your task is to generate a NEWS.md entry from a git diff.

CRITICAL REQUIREMENTS:

• Write in CRAN NEWS.md style
• Focus ONLY on user-visible changes
• Ignore internal refactoring unless it affects users
• Ignore formatting-only changes
• Ignore test changes unless they fix bugs

Group into sections:

## NEW FEATURES
## BUG FIXES
## IMPROVEMENTS
## PERFORMANCE
## INTERNAL

Be concise and professional.

Do NOT mention specific file names.

Output markdown only.
"""

def main():

    diff = sys.stdin.read().strip()

    if not diff:
        print("No diff provided.", file=sys.stderr)
        sys.exit(1)

    api_key = os.getenv("RCPPALGOS_OPENAI_KEY")
    client = OpenAI(api_key=api_key)

    resp = client.responses.create(
        model=MODEL,
        input=[
            {"role": "system", "content": SYSTEM_PROMPT},
            {"role": "user", "content": diff},
        ],
        max_output_tokens=2000,
    )

    print(resp.output_text.strip())


if __name__ == "__main__":
    main()
