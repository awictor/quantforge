"""Docs stay in sync with the code.

The API reference is generated from the package by docs/gen_api.py; this test
runs it in --check mode so a public API change without regenerating the docs
fails CI. Also asserts every public callable/class is documented.
"""

import sys
import os
import inspect

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "docs"))

import gen_api  # noqa: E402
import quantforge  # noqa: E402


def test_api_reference_is_up_to_date():
    # Exit code 0 means docs/api.md matches the current package.
    assert gen_api.main(["--check"]) == 0


def test_every_public_callable_is_documented():
    md = gen_api.build_markdown()
    for name in quantforge.__all__:
        if name == "__version__":
            continue
        obj = getattr(quantforge, name)
        if inspect.isfunction(obj) or inspect.isclass(obj):
            assert f"### `{name}" in md, f"{name} missing from the API reference"


def test_markdown_has_module_sections():
    md = gen_api.build_markdown()
    # A few representative modules should appear as sections.
    for mod in ("bsm", "exotics", "heston", "rates", "risk"):
        assert f"## {mod}" in md
