__all__ = ["__version__"]

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version(__name__)
except PackageNotFoundError:
    # package is not installed
    pass
del version, PackageNotFoundError
