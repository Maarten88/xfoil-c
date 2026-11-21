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
import numpy as np
import ctypes
import os
import glob

from ctypes import c_bool, c_int, c_double, byref, POINTER, cdll
from shutil import copy2
from tempfile import NamedTemporaryFile

from .model import Airfoil

# xfoil_c: the C library (libxfoil_c) instead of xfoil_fortran, located by _lib.py, and loaded with
# _lib.cdll, which raises XFoilError after a hard library error (the Fortran library would end the
# process); otherwise this file is the original xfoil/xfoil.py (drop-in for `import xfoil`).
from ._lib import lib_path, cdll  # noqa: F811  (replaces ctypes.cdll imported above)

here = os.path.abspath(os.path.dirname(__file__))
lib_ext = lib_path[lib_path.rfind('.'):]

fptr = POINTER(c_double)
bptr = POINTER(c_bool)


def _range_count(start, end, step):
    """xfoil_c (clean-up C6c): number of points of range(start, end, step) for floats, and the end
    value that makes the library's step (end - start)/n equal to `step`. The original used
    abs(int((end - start)/step)), which loses a point to rounding (aseq(0, 0.7, 0.1) gave 6 points
    0, 0.1167, ...) and then spread the points over the range."""
    x = abs((end - start) / step)
    n = max(0, int(np.ceil(x - 1e-9 * max(1.0, x))))
    return n, start + n * abs(step) * (1.0 if end >= start else -1.0)


class XFoil(object):
    """Interface to the XFoil Fortran routines.

    Attributes
    ----------
    airfoil
    Re
    M
    xtr
    n_crit
    max_iter
    """

    def __init__(self):
        super().__init__()
        # tmp = NamedTemporaryFile(mode='wb', delete=False, suffix=lib_ext)
        # tmp.close()
        # self._lib_path = tmp.name
        # copy2(lib_path, self._lib_path)
        # self._lib = cdll.LoadLibrary(self._lib_path)
        self._lib = cdll.LoadLibrary(lib_path)
        self._lib.init()
        self._airfoil = None

        self._lib.get_print.restype = c_bool
        self._lib.get_reynolds.restype = c_double
        self._lib.get_mach.restype = c_double
        self._lib.get_n_crit.restype = c_double

    def __del__(self):
        name = self._lib._name
        del self._lib
        try:
            ctypes.windll.kernel32.FreeLibrary(name)
        except AttributeError:
            pass
        # finally:
        #     os.remove(self._lib_path)

    @property
    def print(self):
        """bool: True if console output should be shown."""
        return self._lib.get_print()

    @print.setter
    def print(self, value):
        self._lib.set_print(byref(c_bool(value)))

    @property
    def airfoil(self):
        """Airfoil: Instance of the Airfoil class."""
        n = self._lib.get_n_coords()
        x = np.asfortranarray(np.zeros(n), dtype=c_double)
        y = np.asfortranarray(np.zeros(n), dtype=c_double)
        self._lib.get_airfoil(x.ctypes.data_as(fptr), y.ctypes.data_as(fptr), byref(c_int(n)))
        return Airfoil(x.astype(float), y.astype(float))

    @airfoil.setter
    def airfoil(self, airfoil):
        self._airfoil = airfoil
        self._lib.set_airfoil(
            np.asfortranarray(airfoil.x.flatten(), dtype=c_double).ctypes.data_as(fptr),
            np.asfortranarray(airfoil.y.flatten(), dtype=c_double).ctypes.data_as(fptr),
            byref(c_int(airfoil.n_coords))
        )

    @property
    def Re(self):
        """float: Reynolds number."""
        return float(self._lib.get_reynolds())

    @Re.setter
    def Re(self, value):
        self._lib.set_reynolds(byref(c_double(value)))

    @property
    def M(self):
        """float: Mach number."""
        return float(self._lib.get_mach())

    @M.setter
    def M(self, value):
        self._lib.set_mach(byref(c_double(value)))

    @property
    def xtr(self):
        """tuple(float, float): Top and bottom flow trip x/c locations."""
        xtr_top = c_double()
        xtr_bot = c_double()
        self._lib.get_xtr(byref(xtr_top), byref(xtr_bot))
        return xtr_top.value, xtr_bot.value

    @xtr.setter
    def xtr(self, value):
        self._lib.set_xtr(byref(c_double(value[0])), byref(c_double(value[1])))

    @property
    def n_crit(self):
        """float: Critical amplification ratio."""
        return float(self._lib.get_n_crit())

    @n_crit.setter
    def n_crit(self, value):
        self._lib.set_n_crit(byref(c_double(value)))

    @property
    def max_iter(self):
        """int: Maximum number of iterations."""
        return int(self._lib.get_max_iter())

    @max_iter.setter
    def max_iter(self, max_iter):
        self._lib.set_max_iter(byref(c_int(max_iter)))

    def naca(self, specifier):
        """Set a NACA 4 or 5 series airfoil.

        Parameters
        ----------
        specifier : string
            A NACA 4 or 5 series identifier, such as '2412'.
        """
        self._lib.set_naca(byref(c_int(int(specifier))))

    def reset_bls(self):
        """Reset the boundary layers to be reinitialized on the next analysis."""
        self._lib.reset_bls()

    def repanel(self, n_nodes=160, cv_par=1, cte_ratio=0.15, ctr_ratio=0.2, xt_ref=(1, 1), xb_ref=(1, 1)):
        """Re-panel airfoil.

        Parameters
        ----------
        n_nodes : int
            Number of panel nodes
        cv_par : float
            Panel bunching parameter
        cte_ratio : float
            TE/LE panel density ratio
        ctr_ratio : float
            Refined-area/LE panel density ratio
        xt_ref : tuple of two floats
            Top side refined area x/c limits
        xb_ref : tuple of two floats
            Bottom side refined area x/c limits
        """
        self._lib.repanel(byref(c_int(n_nodes)), byref(c_double(cv_par)),
                          byref(c_double(cte_ratio)), byref(c_double(ctr_ratio)),
                          byref(c_double(xt_ref[0])), byref(c_double(xt_ref[1])),
                          byref(c_double(xb_ref[0])), byref(c_double(xb_ref[1])))

    def a(self, a):
        """Analyze airfoil at a fixed angle of attack.

        Parameters
        ----------
        a : float
            Angle of attack in degrees

        Returns
        -------
        cl, cd, cm, cp : float
            Corresponding values of the lift, drag, moment, and minimum pressure coefficients.
        """
        cl = c_double()
        cd = c_double()
        cm = c_double()
        cp = c_double()
        xtrtop = c_double()
        xtrbot = c_double()
        conv = c_bool()

        self._lib.alfa(byref(c_double(a)), byref(cl), byref(cd), byref(cm), byref(cp), byref(xtrtop), byref(xtrbot), byref(conv))

        return (cl.value, cd.value, cm.value, cp.value, xtrtop.value, xtrbot.value) if conv else (np.nan, np.nan, np.nan, np.nan, np.nan, np.nan)

    def cl(self, cl):
        """"Analyze airfoil at a fixed lift coefficient.

        Parameters
        ----------
        cl : float
            Lift coefficient

        Returns
        -------
        a, cd, cm, cp, xtr_top, xtr_bot : float
            Corresponding values of the angle of attack, drag, moment, and minimum pressure coefficients.
        """
        a = c_double()
        cd = c_double()
        cm = c_double()
        cp = c_double()
        xtrtop = c_double()
        xtrbot = c_double()
        conv = c_bool()

        self._lib.cl(byref(c_double(cl)), byref(a), byref(cd), byref(cm), byref(cp), byref(xtrtop), byref(xtrbot), byref(conv))

        return (a.value, cd.value, cm.value, cp.value, xtrtop.value, xtrbot.value) if conv else (np.nan, np.nan, np.nan, np.nan, np.nan, np.nan)

    def aseq(self, a_start, a_end, a_step):
        """Analyze airfoil at a sequence of angles of attack.

        The analysis is done for the angles of attack given by range(a_start, a_end, a_step).

        Parameters
        ----------
        a_start, a_end, a_step : float
            Start, end, and increment angles for the range.

        Returns
        -------
        a, cl, cd, cm, co : np.ndarray
            Lists of angles of attack and their corresponding lift, drag, moment, and minimum pressure coefficients.
        """
        n, a_end = _range_count(a_start, a_end, a_step)

        a = np.zeros(n, dtype=c_double)
        cl = np.zeros(n, dtype=c_double)
        cd = np.zeros(n, dtype=c_double)
        cm = np.zeros(n, dtype=c_double)
        cp = np.zeros(n, dtype=c_double)
        conv = np.zeros(n, dtype=c_bool)

        self._lib.aseq(byref(c_double(a_start)), byref(c_double(a_end)), byref(c_int(n)),
                       a.ctypes.data_as(fptr), cl.ctypes.data_as(fptr),
                       cd.ctypes.data_as(fptr), cm.ctypes.data_as(fptr),
                       cp.ctypes.data_as(fptr), conv.ctypes.data_as(bptr))

        isnan = np.logical_not(conv)
        a[isnan] = np.nan
        cl[isnan] = np.nan
        cd[isnan] = np.nan
        cm[isnan] = np.nan
        cp[isnan] = np.nan
        return a.astype(float), cl.astype(float), cd.astype(float), cm.astype(float), cp.astype(float)

    def cseq(self, cl_start, cl_end, cl_step):
        """Analyze airfoil at a sequence of lift coefficients.

        The analysis is done for the lift coefficients given by range(cl_start, cl_end, cl_step).

        Parameters
        ----------
        cl_start, cl_end, cl_step : float
            Start, end, and increment lift coefficients for the range.

        Returns
        -------
        a, cl, cd, cm, co : np.ndarray
            Lists of angles of attack and their corresponding lift, drag, moment, and minimum pressure coefficients.
        """
        n, cl_end = _range_count(cl_start, cl_end, cl_step)

        a = np.zeros(n, dtype=c_double)
        cl = np.zeros(n, dtype=c_double)
        cd = np.zeros(n, dtype=c_double)
        cm = np.zeros(n, dtype=c_double)
        cp = np.zeros(n, dtype=c_double)
        conv = np.zeros(n, dtype=c_bool)

        self._lib.cseq(byref(c_double(cl_start)), byref(c_double(cl_end)), byref(c_int(n)),
                       a.ctypes.data_as(fptr), cl.ctypes.data_as(fptr),
                       cd.ctypes.data_as(fptr), cm.ctypes.data_as(fptr),
                       cp.ctypes.data_as(fptr), conv.ctypes.data_as(bptr))

        isnan = np.logical_not(conv)
        a[isnan] = np.nan
        cl[isnan] = np.nan
        cd[isnan] = np.nan
        cm[isnan] = np.nan
        cp[isnan] = np.nan
        return a.astype(float), cl.astype(float), cd.astype(float), cm.astype(float), cp.astype(float)

    def get_cp_distribution(self):
        """Get the Cp distribution from the last converged point.

        xfoil_c: cp is NaN if the last analysed point did not converge, or if the airfoil, paneling,
        flap or Mach number changed since (the original returned those stale/unconverged values).

        Returns
        -------
        x : np.array
            X-coordinates
        y : np.array
            Y-coordinates
        cp : np.ndarray
            Pressure coefficients at the corresponding x-coordinates
        """
        n = self._lib.get_n_cp()
        x = np.zeros(n, dtype=c_double)
        y = np.zeros(n, dtype=c_double)
        cp = np.zeros(n, dtype=c_double)

        self._lib.get_cp(x.ctypes.data_as(fptr), y.ctypes.data_as(fptr), cp.ctypes.data_as(fptr), byref(c_int(n)))
        return x.astype(float), y.astype(float), cp.astype(float)

    # ------------------------------------------------------------------------------------------
    # xfoil_c extensions (not in the original xfoil package; see MIGRATION_PLAN.md Phase 10b).
    # Everything above this line is the original xfoil/xfoil.py.
    # ------------------------------------------------------------------------------------------

    def _geometry(self):
        """Geometry of the current (paneled) airfoil: thick, x_thick, camber, x_camber, le_radius,
        te_angle (radians), area. NaN if no airfoil is set."""
        out = [c_double() for _ in range(7)]
        self._lib.get_geometry(*[byref(o) for o in out])
        return tuple(o.value for o in out)

    @property
    def max_thickness(self):
        """float: Maximum thickness of the current airfoil (read-only).

        As XFoil computes it (TCCALC): the largest distance normal to the chord line between a node
        and the opposite surface, taken at the discrete nodes. In airfoil coordinates, i.e. a fraction
        of the chord for a unit-chord airfoil."""
        return self._geometry()[0]

    @property
    def x_max_thickness(self):
        """float: x location of the maximum thickness (read-only; x of the opposite-surface point, as XFoil)."""
        return self._geometry()[1]

    @property
    def max_camber(self):
        """float: Maximum camber of the current airfoil (read-only; signed, same method as max_thickness)."""
        return self._geometry()[2]

    @property
    def x_max_camber(self):
        """float: x location of the maximum camber (read-only; x of the opposite-surface point, as XFoil)."""
        return self._geometry()[3]

    @property
    def le_radius(self):
        """float: Leading-edge radius of the current airfoil (read-only; 0 if the LE curvature is tiny)."""
        return self._geometry()[4]

    @property
    def te_angle(self):
        """float: Trailing-edge angle of the current airfoil in degrees (read-only)."""
        return float(np.degrees(self._geometry()[5]))

    @property
    def area(self):
        """float: Cross-sectional area of the current airfoil (read-only)."""
        return self._geometry()[6]

    # -- flap and hinge moment (XFOIL GDES FLAP / OPER FMOM, without prompts) -------------------

    def set_hinge(self, x_hinge, y_hinge, y_relative=False):
        """Set the flap hinge point; hinge moments are then computed at every analysed point.

        Parameters
        ----------
        x_hinge : float
            Hinge x location (airfoil coordinates).
        y_hinge : float
            Hinge y location; with y_relative=True the fraction of the local thickness at x_hinge
            (0 = bottom surface, 1 = top surface, 0.5 = halfway), as in XFOIL.
        y_relative : bool
            Interpret y_hinge as y/t.
        """
        if self._lib.get_n_coords() < 2:
            raise ValueError('set an airfoil before setting the hinge')
        self._lib.set_hinge(byref(c_double(x_hinge)), byref(c_double(y_hinge)), byref(c_bool(y_relative)))

    def flap(self, angle):
        """Set the flap deflection about the hinge (set_hinge) and repanel with the current paneling settings.

        The angle is absolute: the undeflected airfoil is deflected once by `angle` (flap(5) then
        flap(5) is still 5 degrees; flap(0) restores the undeflected airfoil). As in XFOIL's FLAP, points
        are added/removed near the surface breaks. xf.airfoil returns the deflected airfoil; setting a new
        airfoil (or naca) clears the hinge and the deflection. Moving the hinge keeps the deflection.

        Parameters
        ----------
        angle : float
            Flap deflection in degrees, positive = trailing edge down.
        """
        if self.hinge is None:
            raise RuntimeError('set the hinge with set_hinge() before deflecting the flap')
        self._lib.flap(byref(c_double(angle)))

    @property
    def flap_angle(self):
        """float: Current flap deflection in degrees, + = trailing edge down (read-only; set with flap())."""
        f = self._lib.get_flap_angle
        f.restype = c_double
        return f()

    @property
    def hinge(self):
        """tuple(float, float) or None: Hinge point (x, y) in airfoil coordinates, None if no hinge is set (read-only)."""
        x, y = c_double(), c_double()
        if not self._lib.get_hinge_point(byref(x), byref(y)):
            return None
        return x.value, y.value

    def _hinge(self):
        out = [c_double() for _ in range(3)]
        self._lib.get_hinge(*[byref(o) for o in out])
        return tuple(o.value for o in out)

    @property
    def hinge_moment(self):
        """float: Flap hinge moment coefficient about the hinge of the last analysed point (read-only).

        Same definition and sign convention as XFOIL's HMOM (integrated surface Cp of the flap part,
        per unit span, in units of q*c^2 for a unit chord). NaN if no hinge is set or no point has been
        analysed since set_hinge/flap."""
        return self._hinge()[0]

    @property
    def hinge_fx(self):
        """float: x force coefficient on the flap of the last analysed point (read-only; NaN as hinge_moment)."""
        return self._hinge()[1]

    @property
    def hinge_fy(self):
        """float: y force coefficient on the flap of the last analysed point (read-only; NaN as hinge_moment)."""
        return self._hinge()[2]

    BL_FIELDS = ('s', 'x', 'y', 'ue', 'dstar', 'theta', 'cf', 'hk', 'hstar', 'n', 'ctau')

    def get_bl_nodes(self):
        """Boundary layer of the last analysed point per node, as XFOIL's OPER DUMP writes it.

        Returns
        -------
        dict of np.ndarray
            Keys BL_FIELDS, one value per node: the airfoil nodes in airfoil order (as the current,
            paneled airfoil and get_cp_distribution), then the wake nodes. s: arc length, x, y: node
            position, ue: Ue/Vinf (signed: negative on the bottom surface, as in DUMP), dstar, theta:
            displacement and momentum thickness, cf: skin friction (0 in the wake), hk: kinematic shape
            factor, hstar: energy shape factor H* (NaN in the wake), n: amplification factor on laminar
            stations (else NaN), ctau: max shear coefficient on turbulent and wake stations (else NaN).
            All but s, x, y are NaN unless the last analysed point is a converged viscous one.
        """
        n = self._lib.get_n_bl()
        out = {k: np.zeros(n, dtype=c_double) for k in self.BL_FIELDS}
        self._lib.get_bl(*[out[k].ctypes.data_as(fptr) for k in self.BL_FIELDS], byref(c_int(n)))
        return {k: v.astype(float) for k, v in out.items()}

    def get_bl(self):
        """Boundary layer of the last analysed point per surface (XFOIL's OPER DUMP, split up).

        Returns
        -------
        dict or None
            None unless the last analysed point is a converged viscous one. Otherwise keys 'top',
            'bottom' and 'wake', each a dict of np.ndarray with the keys of BL_FIELDS (see
            get_bl_nodes), ordered downstream: top and bottom from the stagnation point to the trailing
            edge, the wake from the trailing edge. Here ue is the magnitude Ue/Vinf on both surfaces.
            Transition: the first station with a finite ctau (n is finite before it).
        """
        nodes = self.get_bl_nodes()
        n_air = self._lib.get_n_cp()                # airfoil nodes; the rest is wake
        ue = nodes['ue'][:n_air]
        if np.isnan(ue).all():
            return None
        top = np.nonzero(ue >= 0.0)[0][::-1]
        bottom = np.nonzero(ue < 0.0)[0]
        wake = np.arange(n_air, len(nodes['x']))
        res = {}
        for side, idx in (('top', top), ('bottom', bottom), ('wake', wake)):
            res[side] = {k: v[idx] for k, v in nodes.items()}
            res[side]['ue'] = np.abs(res[side]['ue'])
        return res

    # A new airfoil clears the hinge (the Fortran library would keep it): the original setter and
    # naca() above run unchanged, then the hinge is removed.
    _original_airfoil = airfoil
    _original_naca = naca

    @_original_airfoil.setter
    def airfoil(self, airfoil):
        XFoil._original_airfoil.fset(self, airfoil)
        self._lib.clear_hinge()

    def naca(self, specifier):
        """Set a NACA 4 or 5 series airfoil.

        Parameters
        ----------
        specifier : string
            A NACA 4 or 5 series identifier, such as '2412'.
        """
        XFoil._original_naca(self, specifier)
        self._lib.clear_hinge()
