import math
import random

from quantforge import particle_filter, pcg_gaussian
from quantforge.kalman import kalman_local_level


def _local_level_data(seed=0, n=60, q=0.02, r=0.5):
    random.seed(seed)
    true = [0.0]
    for _ in range(n - 1):
        true.append(true[-1] + random.gauss(0, math.sqrt(q)))
    obs = [t + random.gauss(0, math.sqrt(r)) for t in true]
    return true, obs, q, r


def _pf_local_level(obs, q, r, n_particles=5000, seed=42):
    def init(rng):
        return pcg_gaussian(rng, 0.0, 1.0)

    def trans(x, rng):
        return x + pcg_gaussian(rng, 0.0, math.sqrt(q))

    def loglik(z, x):
        return -0.5 * ((z - x) ** 2 / r + math.log(2 * math.pi * r))

    return particle_filter(obs, trans, loglik, init, n_particles=n_particles, seed=seed)


def test_tracks_kalman_on_linear_gaussian():
    true, obs, q, r = _local_level_data()
    lev, _, _ = kalman_local_level(obs, q, r, x0=0.0, p0=1.0)
    pm = _pf_local_level(obs, q, r)["means"]
    maxdiff = max(abs(pm[i] - lev[i]) for i in range(20, 60))
    assert maxdiff < 0.15


def test_reproducible_with_seed():
    _, obs, q, r = _local_level_data()
    a = _pf_local_level(obs, q, r, seed=42)
    b = _pf_local_level(obs, q, r, seed=42)
    assert a["means"] == b["means"]


def test_seed_changes_result_but_stays_close():
    _, obs, q, r = _local_level_data()
    a = _pf_local_level(obs, q, r, seed=42)
    c = _pf_local_level(obs, q, r, seed=7)
    assert a["means"] != c["means"]
    assert abs(a["means"][-1] - c["means"][-1]) < 0.1


def test_ess_bounded_and_resampling_happens():
    _, obs, q, r = _local_level_data()
    res = _pf_local_level(obs, q, r)
    assert all(0 < e <= 5000 + 1e-6 for e in res["ess"])
    assert res["n_resample"] > 0


def test_vector_state_position_velocity():
    random.seed(1)
    tp = [0.0]
    v = 0.3
    for _ in range(39):
        tp.append(tp[-1] + v)
    z = [[p + random.gauss(0, 0.5)] for p in tp]

    def init(rng):
        return [pcg_gaussian(rng, 0, 1), pcg_gaussian(rng, 0, 1)]

    def trans(s, rng):
        return [s[0] + s[1], s[1] + pcg_gaussian(rng, 0, 0.05)]

    def loglik(zz, s):
        return -0.5 * (zz[0] - s[0]) ** 2 / 0.25

    res = particle_filter(z, trans, loglik, init, n_particles=3000, seed=1)
    est = [m[0] for m in res["means"]]
    rmse = math.sqrt(sum((est[i] - tp[i]) ** 2 for i in range(10, 40)) / 30)
    assert rmse < 1.5
