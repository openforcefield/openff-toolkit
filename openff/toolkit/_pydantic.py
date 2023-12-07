try:
    from pydantic.v1 import PrivateAttr
except ModuleNotFoundError:
    from pydantic import PrivateAttr
