# SPDX-FileCopyrightText: 2023 Paul Romano
# SPDX-License-Identifier: MIT

from importlib.metadata import version, PackageNotFoundError

from .material import *
from .incident_neutron import *
from .function import *
from .product import *
from .reaction import *
from . import ace

try:
    __version__ = version("endf")
except PackageNotFoundError:
    # package is not installed
    pass
