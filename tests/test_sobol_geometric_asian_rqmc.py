"""RQMC geometric-average Asian vs exact closed form (sobol_geometric_asian_rqmc)."""

import pytest

from quantforge import sobol_geometric_asian_rqmc, OptionType
from quantforge.montecarlo import _discrete_geometric_asian


S, K, T, R, SIG, NS = 100.0, 100.0, 1.0, 0.05, 0.2, 6


@pytest.mark.parametrize("ot", [OptionType.CALL, OptionType.PUT])
def test_matches_exact_closed_form(ot):
    # The discrete geometric average is exactly lognormal, so there is a closed
    # form -- the tightest possible (deterministic) reference.
    cf = _discrete_geometric_asian(S, K, T, R, SIG, ot, R, NS)
    rq = sobol_geometric_asian_rqmc(S, K, T, R, SIG, ot, n_steps=NS,
                                    n_paths=4096, n_rand=24, seed=1)
    assert rq.price == pytest.approx(cf, abs=3.0 * rq.std_error)


@pytest.mark.parametrize("ns", list(range(2, 13)))
def test_extended_sobol_dims_match_closed_form(ns):
    # Locks the direction-number table (dims 7-12 added in 1.209.0): the RQMC geo
    # Asian must match the exact closed form at every supported n_steps. A wrong
    # primitive polynomial or seed would break the low-discrepancy property and
    # bias this.
    cf = _discrete_geometric_asian(S, K, T, R, SIG, OptionType.CALL, R, ns)
    # n_rand=8 keeps the honest SE small enough that 4*SE still catches a broken
    # direction number, at a fraction of the runtime of the full n_rand=24 batch.
    rq = sobol_geometric_asian_rqmc(S, K, T, R, SIG, OptionType.CALL,
                                    n_steps=ns, n_paths=4096, n_rand=8, seed=7)
    assert rq.price == pytest.approx(cf, abs=4.0 * rq.std_error)


def test_dividend_carry_matches_closed_form():
    b = R - 0.03
    cf = _discrete_geometric_asian(S, K, T, R, SIG, OptionType.CALL, b, NS)
    rq = sobol_geometric_asian_rqmc(S, K, T, R, SIG, OptionType.CALL, b=b,
                                    n_steps=NS, n_paths=4096, n_rand=24, seed=2)
    assert rq.price == pytest.approx(cf, abs=3.0 * rq.std_error)


def test_honest_se_small():
    rq = sobol_geometric_asian_rqmc(S, K, T, R, SIG, OptionType.CALL, n_steps=NS,
                                    n_paths=4096, n_rand=24, seed=3)
    assert rq.std_error < 0.02


def test_n_paths_is_total_points():
    rq = sobol_geometric_asian_rqmc(S, K, T, R, SIG, n_steps=4, n_paths=2048,
                                    n_rand=16, seed=4)
    assert rq.n_paths == 2048 * 16


def test_bad_n_rand_raises():
    with pytest.raises(ValueError):
        sobol_geometric_asian_rqmc(S, K, T, R, SIG, n_rand=1)


def test_bad_n_steps_raises():
    with pytest.raises(ValueError):
        sobol_geometric_asian_rqmc(S, K, T, R, SIG, n_steps=99)
