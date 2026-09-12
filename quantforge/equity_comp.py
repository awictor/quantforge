"""Warrants and employee stock options (ESOs) with dilution and forfeiture.

Warrants issued by the company create new shares on exercise, diluting existing
holders, so their value is the vanilla call scaled by ``M / (M + N)`` for ``M``
existing and ``N`` new (warrant) shares. Employee stock options are long-dated
calls whose holders exercise early and forfeit before vesting; the FASB-endorsed
practical model (Hull-White) values them at an expected life shorter than the
contractual term and multiplies by the probability of surviving to exercise given
a post-vest exit rate. Pure standard library on top of :mod:`quantforge.bsm`.
"""

import math

from .bsm import call_price


def dilution_factor(existing_shares, new_shares):
    """Dilution multiplier ``M / (M + N)`` for issuing ``N`` new shares on ``M``.

    Below one whenever new shares are created; one when ``N = 0``. Applied to a
    warrant's per-option payoff because exercise expands the share count.
    """
    if existing_shares <= 0:
        raise ValueError("existing_shares must be positive")
    if new_shares < 0:
        raise ValueError("new_shares must be non-negative")
    return existing_shares / (existing_shares + new_shares)


def warrant_price(S, K, t, r, sigma, existing_shares, new_shares, b=None):
    """Warrant value: the vanilla call scaled by the :func:`dilution_factor`.

    A single-period dilution adjustment -- the standard textbook approximation
    ``value = M/(M+N) * call(S, K, ...)`` -- which reduces to the plain call when
    no new shares are issued. ``S`` is the current (pre-dilution) share price.
    """
    factor = dilution_factor(existing_shares, new_shares)
    return factor * call_price(S, K, t, r, sigma, b)


def eso_expected_life(vesting, contractual_term, exit_rate):
    """Expected life of an ESO given post-vest exit and the contractual term.

    After vesting, holders leave (and exercise or forfeit) at a constant hazard
    ``exit_rate``; the expected time to exercise, capped at the contractual term,
    is

        vesting + (1 - e^{-exit_rate * (T - vesting)}) / exit_rate

    -- between ``vesting`` (immediate exit) and ``contractual_term`` (no exit).
    """
    if not (0.0 <= vesting <= contractual_term):
        raise ValueError("vesting must lie in [0, contractual_term]")
    if exit_rate < 0:
        raise ValueError("exit_rate must be non-negative")
    post = contractual_term - vesting
    if exit_rate == 0.0:
        return contractual_term
    return vesting + (1.0 - math.exp(-exit_rate * post)) / exit_rate


def eso_value(S, K, contractual_term, r, sigma, vesting, exit_rate,
              forfeiture_rate=0.0, b=None):
    """Employee stock option value (Hull-White practical / FASB 123R style).

    Prices the ESO as a call at the :func:`eso_expected_life` rather than the full
    term, then multiplies by the probability of surviving pre-vest forfeiture
    ``e^{-forfeiture_rate * vesting}``. Cheaper than the vanilla call on the
    contractual term because early exercise shortens the option and forfeiture
    can extinguish it. With ``exit_rate = 0`` and ``forfeiture_rate = 0`` it
    reduces to the vanilla call.
    """
    if forfeiture_rate < 0:
        raise ValueError("forfeiture_rate must be non-negative")
    life = eso_expected_life(vesting, contractual_term, exit_rate)
    survive = math.exp(-forfeiture_rate * vesting)
    return survive * call_price(S, K, life, r, sigma, b)
