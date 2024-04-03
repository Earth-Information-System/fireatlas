try:
    from ._version import __version__
except ImportError:
    __version__ = "unknown"

# do not create top-level imports here
# that import from `.FireConsts` until
# we've moved to pydantic so that os env var
# substitution in FireConsts can work for us
from fireatlas import (
    FireTypes,
    FireTime,
    FireEnums,
)

from fireatlas.FireConsts import Settings

settings = Settings()
