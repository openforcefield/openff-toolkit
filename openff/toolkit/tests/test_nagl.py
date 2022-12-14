from openff.utilities.testing import skip_if_missing

from openff.toolkit.utils._nagl_wrapper import _NAGLToolkitWrapper


@skip_if_missing("openff.nagl")
def test_version():
    assert "0.1.0" in _NAGLToolkitWrapper()._toolkit_version
