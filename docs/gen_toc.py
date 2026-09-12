"""Generate the README table of contents between marker comments.

Scans README.md for ``## `` headings and rewrites the block between
``<!-- TOC -->`` and ``<!-- /TOC -->`` with a bulleted, anchored list. Run with
``--check`` to verify the TOC is current (used by the test suite); run with no
args to rewrite it in place. GitHub-style anchors: lowercase, spaces to hyphens,
punctuation stripped.
"""

import os
import re
import sys

START = "<!-- TOC -->"
END = "<!-- /TOC -->"


def _readme_path():
    return os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                        "README.md")


def _anchor(title):
    a = title.strip().lower()
    a = a.replace(" ", "-")
    a = re.sub(r"[^\w\-]", "", a)   # strip punctuation (keep word chars and hyphen)
    return a


def build_toc(text):
    lines = text.splitlines()
    items = []
    in_code = False
    for line in lines:
        if line.startswith("```"):
            in_code = not in_code
            continue
        if in_code:
            continue
        if line.startswith("## ") and not line.startswith("### "):
            title = line[3:].strip()
            items.append(f"- [{title}](#{_anchor(title)})")
    return "\n".join(items)


def apply_toc(text, toc):
    if START not in text or END not in text:
        raise ValueError("README missing TOC markers")
    pre, rest = text.split(START, 1)
    _, post = rest.split(END, 1)
    return f"{pre}{START}\n{toc}\n{END}{post}"


def main(argv=None):
    argv = argv if argv is not None else sys.argv[1:]
    path = _readme_path()
    with open(path, encoding="utf-8") as f:
        text = f.read()
    toc = build_toc(text)
    new_text = apply_toc(text, toc)
    if "--check" in argv:
        if new_text != text:
            print("README TOC is stale; run: python docs/gen_toc.py", file=sys.stderr)
            return 1
        print("README TOC is up to date.")
        return 0
    with open(path, "w", encoding="utf-8") as f:
        f.write(new_text)
    print(f"wrote TOC ({toc.count(chr(10)) + 1} entries)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
