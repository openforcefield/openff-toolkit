import pytest

from openff.toolkit._utilities.testing import skip_if_missing, skip_if_missing_exec


def test_skips():
    assert not skip_if_missing("pip").args[0]
    assert skip_if_missing("numpynumpy").args[0]

    assert not skip_if_missing_exec("python").args[0]
    assert skip_if_missing_exec("python4").args[0]
    assert not skip_if_missing_exec(["python", "python4"]).args[0]
    assert skip_if_missing_exec(["python4", "python5"]).args[0]

    with pytest.raises(ValueError, match=r"Bad type.*int"):
        skip_if_missing_exec(0)
