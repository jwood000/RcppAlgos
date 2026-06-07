#!/usr/bin/env python3

import os
import sys
from openai import OpenAI

MODEL = os.getenv("AI_REVIEW_MODEL", "gpt-5.5")

SYSTEM_PROMPT = """
You are helping maintain the R package RcppAlgos. You are an expert C++
developer, algorithm designer, and understand the R and CRAN ecosystem deeply.

Your task is to generate a professional GitHub Pull Request summary from a git diff.

CRITICAL REQUIREMENTS:

• This summary will be used as the PR description.
• Focus on WHAT changed and WHY.
• Be accurate and specific.
• Do NOT invent behavior not supported by the diff.
• If intent is inferred, state it as inference.

• Highlight important technical implications:
  - correctness
  - edge cases
  - performance
  - memory safety
  - API or behavioral changes

• Call out anything reviewers should pay special attention to.

• Prefer clarity over verbosity.
• Do NOT include code blocks.
• Do NOT include raw diff content.
• Do NOT mention file paths unless absolutely necessary for clarity.

Use the following structure exactly:

## Summary

High-level explanation of the change and its purpose.

## Key Changes

• Bullet list of concrete technical changes

## Motivation

Identify the PRIMARY motivation that drove this PR.

Pay special attention to new combinatorial algorithms and treat them as
likely primary motivations.

This should reflect the original idea or feature that initiated the work,
not secondary refactors or improvements discovered along the way.

Then include a subsection:

### Secondary Improvements

List improvements that were discovered opportunistically while implementing
the primary change.

Be careful to distinguish clearly between:

• the original driving idea
• incidental fixes, refactors, or optimizations

If uncertain, infer carefully from the diff and state that it is inferred.

## Impact

Describe any impact on:

• User-facing behavior
• Performance
• Correctness
• Backwards compatibility

If none, state "No user-visible impact expected."

## Reviewer Notes

Call out:

• Risky areas
• Edge cases
• Logic changes
• Performance-sensitive code
• CRAN / R API considerations if applicable

Be professional and concise.

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
        max_output_tokens=4000,
    )

    print(resp.output_text.strip())


if __name__ == "__main__":
    main()
