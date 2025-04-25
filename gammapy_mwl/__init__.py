
from .version import version as __version__
from importlib.metadata import version as _version, PackageNotFoundError

try:
    __version__ = _version(__name__)
except PackageNotFoundError:
    pass
del version, PackageNotFoundError


__all__ = []
