try:
    from ._version import __version__
except ImportError:
    __version__ = "unknown"

from fireatlas import (
    FireTypes,
    FireTime,
)

from fireatlas.FireConsts import Settings

# instantiate pydantic settingss
settings = Settings()
