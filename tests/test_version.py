"""Version metadata sanity checks."""

import re

import quantforge


def test_version_is_semver():
    assert re.match(r"^\d+\.\d+\.\d+$", quantforge.__version__)


def test_version_is_at_least_1_0():
    major = int(quantforge.__version__.split(".")[0])
    assert major >= 1


def test_public_api_is_importable():
    # Every name advertised in __all__ must actually resolve.
    for name in quantforge.__all__:
        assert hasattr(quantforge, name), f"missing public symbol: {name}"
