try:
    from ._version import __version__
except ImportError:
    __version__ = "unknown"

# we cannot make top level imports until we factor out FireConsts os enviroment variables or else they get set here
from . import FireTypes