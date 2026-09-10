"""Swaption volatility cube: a SABR smile per (expiry, tenor) node.

The swaption market quotes implied volatility along three axes -- option expiry,
underlying swap tenor, and strike (moneyness). A **vol cube** stores a calibrated
SABR smile at each ``(expiry, tenor)`` grid node and interpolates across all
three axes: SABR handles the strike/smile dimension analytically, while expiry
and tenor are handled by bilinear interpolation of the total variance
``w = sigma^2 * expiry`` (interpolating variance, not vol, is the
no-calendar-arbitrage-friendly choice).

Build the cube from market smiles with :meth:`VolCube.fit`, then query any
``(expiry, tenor, strike)`` with :meth:`vol`. Pure standard library, built on
:mod:`quantforge.sabr`.
"""

import math
from bisect import bisect_left

from .sabr import calibrate_sabr, sabr_vol, SABRParams


class VolCube:
    """A SABR-per-node swaption vol cube with variance interpolation."""

    def __init__(self, expiries, tenors, node_params, forwards):
        """``node_params[(e, T)]`` is a fitted SABRParams; ``forwards[(e, T)]``
        the forward swap rate for that node."""
        self.expiries = sorted(expiries)
        self.tenors = sorted(tenors)
        self.node_params = node_params
        self.forwards = forwards

    @classmethod
    def fit(cls, market, beta=0.5):
        """Fit a cube from market smiles.

        ``market`` maps ``(expiry, tenor)`` to ``(forward, strikes, vols)``. Each
        node is calibrated with SABR (``beta`` fixed). Returns the cube.
        """
        expiries = sorted({e for e, _ in market})
        tenors = sorted({t for _, t in market})
        params = {}
        forwards = {}
        for (e, T), (fwd, strikes, vols) in market.items():
            p, _rmse = calibrate_sabr(fwd, e, strikes, vols, beta=beta)
            params[(e, T)] = p
            forwards[(e, T)] = fwd
        return cls(expiries, tenors, params, forwards)

    def _bracket(self, axis, x):
        """Return (lo, hi, weight) so x ~ (1-w)*axis[lo] + w*axis[hi]."""
        if x <= axis[0]:
            return 0, 0, 0.0
        if x >= axis[-1]:
            n = len(axis) - 1
            return n, n, 0.0
        i = bisect_left(axis, x)
        if axis[i] == x:
            return i, i, 0.0
        lo, hi = i - 1, i
        w = (x - axis[lo]) / (axis[hi] - axis[lo])
        return lo, hi, w

    def _node_vol(self, e_idx, t_idx, expiry, strike):
        """SABR vol at a grid node for the queried strike, using its own forward."""
        e = self.expiries[e_idx]
        T = self.tenors[t_idx]
        p = self.node_params[(e, T)]
        F = self.forwards[(e, T)]
        return sabr_vol(F, strike, expiry, p.alpha, p.beta, p.rho, p.nu)

    def vol(self, expiry, tenor, strike):
        """Implied vol at ``(expiry, tenor, strike)`` by SABR + variance interp.

        Strike is handled by each node's SABR smile; expiry and tenor by
        bilinear interpolation of the total variance ``sigma^2 * expiry``.
        """
        ei_lo, ei_hi, we = self._bracket(self.expiries, expiry)
        ti_lo, ti_hi, wt = self._bracket(self.tenors, tenor)

        def w_at(ei, ti):
            iv = self._node_vol(ei, ti, expiry, strike)
            return iv * iv * expiry   # total variance at the queried expiry

        w00 = w_at(ei_lo, ti_lo)
        w01 = w_at(ei_lo, ti_hi)
        w10 = w_at(ei_hi, ti_lo)
        w11 = w_at(ei_hi, ti_hi)
        w0 = (1 - wt) * w00 + wt * w01
        w1 = (1 - wt) * w10 + wt * w11
        w = (1 - we) * w0 + we * w1
        if expiry <= 0:
            return math.sqrt(max(w, 0.0))
        return math.sqrt(max(w, 0.0) / expiry)

    def node_smile(self, expiry, tenor, strikes):
        """The fitted SABR smile at an exact ``(expiry, tenor)`` node."""
        p = self.node_params[(expiry, tenor)]
        F = self.forwards[(expiry, tenor)]
        return [(K, sabr_vol(F, K, expiry, p.alpha, p.beta, p.rho, p.nu))
                for K in strikes]
