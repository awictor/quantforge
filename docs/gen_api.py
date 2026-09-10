"""Generate the API reference (docs/api.md) by introspecting the package.

Zero-dependency: walks ``quantforge.__all__``, groups symbols by the module
they live in, and emits each public function/class with its signature and
docstring. Run from the repo root:

    python docs/gen_api.py            # writes docs/api.md
    python docs/gen_api.py --check    # exits non-zero if api.md is stale

Keeping this in the repo means the reference never drifts from the code: the
test suite runs it in --check mode.
"""

import argparse
import inspect
import os
import sys

# Make the package importable when run from the repo root.
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import quantforge  # noqa: E402


def _module_short(obj):
    mod = getattr(obj, "__module__", "") or ""
    return mod.split(".")[-1] if mod.startswith("quantforge") else mod


def _signature(obj):
    try:
        return str(inspect.signature(obj))
    except (TypeError, ValueError):
        return ""


def _first_doc_line(obj):
    doc = inspect.getdoc(obj) or ""
    return doc.strip()


def build_markdown() -> str:
    lines = [
        "# QuantForge API reference",
        "",
        f"Auto-generated from `quantforge` v{quantforge.__version__} by "
        "`docs/gen_api.py` — do not edit by hand.",
        "",
    ]

    # Group public callables/classes by their defining module. Skip module-level
    # constants/flags (e.g. HAS_NUMPY) that are values, not documented API.
    groups = {}
    for name in quantforge.__all__:
        if name == "__version__":
            continue
        obj = getattr(quantforge, name)
        if not (inspect.isfunction(obj) or inspect.isclass(obj)):
            continue
        mod = _module_short(obj)
        groups.setdefault(mod, []).append((name, obj))

    for mod in sorted(groups):
        lines.append(f"## {mod}")
        lines.append("")
        for name, obj in sorted(groups[mod], key=lambda x: x[0]):
            kind = "class" if inspect.isclass(obj) else "function"
            sig = _signature(obj)
            lines.append(f"### `{name}{sig}`  _{kind}_")
            lines.append("")
            doc = inspect.getdoc(obj)
            if doc:
                # Indent the docstring as a block quote for readability.
                for dl in doc.splitlines():
                    lines.append(f"> {dl}" if dl.strip() else ">")
            else:
                lines.append("> (no docstring)")
            lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def _target_path():
    return os.path.join(os.path.dirname(__file__), "api.md")


def main(argv=None):
    ap = argparse.ArgumentParser(description="Generate the QuantForge API reference.")
    ap.add_argument("--check", action="store_true",
                    help="exit non-zero if docs/api.md is out of date")
    args = ap.parse_args(argv)

    md = build_markdown()
    path = _target_path()

    if args.check:
        if not os.path.exists(path):
            print("docs/api.md missing; run: python docs/gen_api.py", file=sys.stderr)
            return 1
        with open(path, encoding="utf-8") as f:
            current = f.read()
        if current != md:
            print("docs/api.md is stale; run: python docs/gen_api.py", file=sys.stderr)
            return 1
        print("docs/api.md is up to date.")
        return 0

    with open(path, "w", encoding="utf-8") as f:
        f.write(md)
    n_documented = md.count("\n### ")
    print(f"wrote {path} ({n_documented} documented symbols)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
