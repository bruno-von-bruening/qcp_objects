
from qcp_objects.util.import_helper import *


try:
    from ._version import __version__
except ImportError:
    warn(f"Could not import version from _version.py (there appears to be a problem with setuptools-scm)")
