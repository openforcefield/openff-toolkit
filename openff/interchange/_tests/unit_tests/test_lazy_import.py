import pytest


def test_lazy_load():
    from openff.interchange import Interchange as I1
    from openff.interchange.components.interchange import Interchange as I2

    assert I1 == I2


def test_lazy_load_top_level_module():
    import openff.interchange

    assert "Interchange" in dir(openff.interchange)


def test_lazy_load_error():
    with pytest.raises(
        ImportError,
        match="cannot import name 'Blah' from 'openff.interchange",
    ):
        from openff.interchange import Blah  # noqa
