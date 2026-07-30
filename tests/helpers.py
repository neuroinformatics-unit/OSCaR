from dataclasses import asdict
from typing import Any, Mapping

import pytest


def assert_dataclass_equal(actual, expected, abs: float = 1e-3) -> None:
    """Assert two dataclasses are equal.

    This can handle dataclasses with nested items e.g. dictionaries within
    dictionaries. It also matches float values within a given tolerance using
    pytest.approx.

    Parameters
    ----------
    actual : dataclass
        The actual dataclass
    expected : dataclass
        The expected dataclass
    abs : float, optional
        The absolute tolerance for comparing float values inside the
        dataclasses. This is passed to pytest.approx()
    """
    actual_dict = asdict(actual)
    expected_dict = asdict(expected)
    _assert_close(actual_dict, expected_dict, abs=abs)


def _assert_close(
    actual: Any, expected: Any, abs: float, message_prefix: str | None = None
) -> None:
    """Assert two values are close, handling nested dicts.

    Parameters
    ----------
    actual : Any
        The actual value
    expected : Any
        The expected value
    abs : float
        The absolute tolerance for comparing float values - this is passed to
        pytest.approx()
    message_prefix : str | None, optional
        Prefix to add to the assertion message. This is useful for nested
        values to indicate where in the structure the error came from.
    """

    # assertion message that will be printed with pytest report
    message = f"actual is {actual} vs in expected {expected}"
    if message_prefix:
        message = f"{message_prefix} in {message}"

    if isinstance(actual, Mapping):
        assert actual.keys() == expected.keys()
        for key in actual:
            if message_prefix:
                new_prefix = f"{message_prefix} - {key}"
            else:
                new_prefix = key

            _assert_close(
                actual[key], expected[key], abs=abs, message_prefix=new_prefix
            )

    elif isinstance(actual, float):
        assert actual == pytest.approx(expected, abs=abs), message

    else:
        assert actual == expected, message
