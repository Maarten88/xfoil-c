# -*- coding: utf-8 -*-
#   Copyright (c) 2019 D. de Vries
#
#   This file is part of XFoil.
#
#   XFoil is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   XFoil is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with XFoil.  If not, see <https://www.gnu.org/licenses/>.
"""C-backed port of xfoil-python: drop-in replacement for `import xfoil`
(`XFoil` in `xfoil_c.xfoil`, `Airfoil` in `xfoil_c.model`, as in the original package).

`lib_version()` returns the version compiled into libxfoil_c, `build_info()` the build variant
("exact" = bit-identical to the Fortran reference).
"""
from . import model
from ._lib import XFoilError, build_info, version as lib_version
from .xfoil import XFoil

__version__ = "0.1.0"
