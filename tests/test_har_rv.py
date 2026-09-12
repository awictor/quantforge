"""HAR-RV realized-variance model (Corsi 2009)."""

import random

import pytest

from quantforge import fit_har_rv, har_rv_forecast


def _har_series():
    random.seed(9)
    rv = [0.0004] * 30
    for t in range(29, 400):
        day = rv[t]
        week = sum(rv[t - 4:t + 1]) / 5
        month = sum(rv[t - 21:t + 1]) / 22
        nxt = 0.0001 + 0.3 * day + 0.4 * week + 0.2 * month + random.gauss(0, 1e-5)
        rv.append(max(nxt, 1e-8))
    return rv


def test_persistence_near_true():
    rv = _har_series()
    _, bd, bw, bm = fit_har_rv(rv)
    # The three lag averages are collinear, so persistence (their sum) is the
    # identified quantity; true value 0.3 + 0.4 + 0.2 = 0.9.
    assert abs((bd + bw + bm) - 0.9) < 0.15


def test_in_sample_fit_correlation_high():
    rv = _har_series()
    b0, bd, bw, bm = fit_har_rv(rv)
    preds, acts = [], []
    for t in range(22, len(rv) - 1):
        d = rv[t]
        w = sum(rv[t - 4:t + 1]) / 5
        m = sum(rv[t - 21:t + 1]) / 22
        preds.append(b0 + bd * d + bw * w + bm * m)
        acts.append(rv[t + 1])
    mp, ma = sum(preds) / len(preds), sum(acts) / len(acts)
    cov = sum((preds[i] - mp) * (acts[i] - ma) for i in range(len(preds)))
    sp = sum((p - mp) ** 2 for p in preds) ** 0.5
    sa = sum((a - ma) ** 2 for a in acts) ** 0.5
    assert cov / (sp * sa) > 0.7


def test_forecast_matches_manual():
    rv = _har_series()
    c = fit_har_rv(rv)
    b0, bd, bw, bm = c
    f = har_rv_forecast(c, rv)
    day, week, month = rv[-1], sum(rv[-5:]) / 5, sum(rv[-22:]) / 22
    assert f == pytest.approx(b0 + bd * day + bw * week + bm * month)


def test_validation():
    with pytest.raises(ValueError):
        fit_har_rv([0.0004] * 10)
    with pytest.raises(ValueError):
        har_rv_forecast((0, 0.3, 0.4, 0.2), [0.0004] * 10)
