from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("oscar-colony")
except PackageNotFoundError:
    # package is not installed
    pass
