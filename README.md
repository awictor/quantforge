# QuantForge

![CI](https://github.com/awictor/quantforge/actions/workflows/ci.yml/badge.svg)
![Python](https://img.shields.io/badge/python-3.8%2B-blue)
![License](https://img.shields.io/badge/license-MIT-green)

Fast, **dependency-free** options pricing and risk engine in pure Python.

QuantForge prices European and American options, computes the full set of
analytic Greeks, and solves for implied volatility — with zero third-party
dependencies. It drops into any Python 3.8+ runtime (serverless, embedded,
notebooks, trading bots) without compiling NumPy or SciPy.

## Why

Most option libraries pull in a heavy scientific stack. QuantForge implements
the math from scratch — a full-precision normal CDF/quantile, closed-form
Black-Scholes-Merton, analytic Greeks, a robust Newton-with-bisection implied
vol solver, and a Cox-Ross-Rubinstein binomial lattice for early exercise. The
whole thing is a few hundred lines of readable, tested code.

- **Zero dependencies** — standard library only.
- **Generalized BSM** — one model covers stocks, dividends, futures (Black-76)
  and FX (Garman-Kohlhagen) via the cost-of-carry parameter.
- **Every Greek analytic** — delta, gamma, vega, theta, rho, verified against
  finite differences to 1e-4 or better.
- **Robust implied vol** — Newton-Raphson with a guaranteed bisection bracket,
  arbitrage-band rejection.
- **American exercise** — binomial tree that converges to BSM on European
  payoffs (used as a self-check in the suite).

## Install

```bash
pip install quantforge          # from PyPI (planned)
pip install -e .                # from source
```

## Quick start

```python
from quantforge import call_price, greeks, implied_volatility, american_price

# Price a European call: spot 100, strike 105, 6 months, 4% rate, 25% vol.
print(call_price(S=100, K=105, t=0.5, r=0.04, sigma=0.25))

# Full risk report in one call.
g = greeks(S=100, K=105, t=0.5, r=0.04, sigma=0.25)
print(g.delta, g.gamma, g.vega, g.theta, g.rho)

# Back out implied vol from a market price.
iv = implied_volatility(target_price=6.12, S=100, K=105, t=0.5, r=0.04)
print(iv)

# American put with a 3% dividend yield (carry b = r - q).
print(american_price(S=100, K=100, t=1.0, r=0.05, sigma=0.2,
                     option_type="put", b=0.05 - 0.03, steps=500))

# Or the closed-form Bjerksund-Stensland (2002) approximation: ~500x faster
# than a 500-step tree, accurate to a few cents.
from quantforge import bjerksund_stensland
print(bjerksund_stensland(S=100, K=100, t=1.0, r=0.05, sigma=0.2,
                          option_type="put", b=0.02))

# Or a trinomial lattice (smoother convergence), with Richardson extrapolation
# for extra accuracy at low step counts.
from quantforge import trinomial_price, richardson_american
print(trinomial_price(S=100, K=100, t=1.0, r=0.05, sigma=0.2,
                      option_type="put", steps=400))
print(richardson_american(S=100, K=100, t=1.0, r=0.05, sigma=0.2,
                          option_type="put", steps=50))
```

## Scenario / stress grid

Revalue a book across a grid of spot and vol shocks — the classic trader "risk
matrix" — and pull out the worst case:

```python
from quantforge import Contract, stress_grid

straddle = [
    Contract(S=100, K=100, t=0.25, r=0.02, sigma=0.3, option_type="call", qty=-1),
    Contract(S=100, K=100, t=0.25, r=0.02, sigma=0.3, option_type="put",  qty=-1),
]

grid = stress_grid(straddle,
                   spot_shocks=[-0.2, -0.1, 0.0, 0.1, 0.2],   # +/-20% spot
                   vol_shocks=[-0.5, 0.0, 0.5],               # +/-50% vol
                   relative=True)

print(grid.worst_case())          # (spot_shock, vol_shock, pnl) of the biggest loss
for ss, row in grid.as_rows():
    print(ss, row)                # P&L per vol shock at each spot shock
```

Shocks are fractional when `relative=True`, or absolute price/vol-point moves
when `relative=False`. `spot_ladder(...)` gives the 1-D spot-only P&L profile.

## Second-order Greeks

For hedging through joint moves in spot, vol, and time:

```python
from quantforge import vanna, vomma, charm, veta, speed, zomma, color

kw = dict(S=100, K=105, t=0.5, r=0.04, sigma=0.25)
vanna(**kw)   # d(delta)/d(sigma) = d(vega)/d(spot)
vomma(**kw)   # d(vega)/d(sigma)  (a.k.a. volga)  -- vol convexity
charm(**kw, option_type="call")   # d(delta)/d(time) -- delta bleed
veta(**kw)    # d(vega)/d(time)   -- vega decay
speed(**kw); zomma(**kw); color(**kw)   # third-order in spot / cross terms
```

Every one is analytic and verified against finite differences of the
first-order Greeks in the test suite.

## Batch pricing and portfolio risk

Value a whole book in one call and get net exposures:

```python
from quantforge import Contract, price_book

book = price_book([
    Contract(S=100, K=105, t=0.5, r=0.04, sigma=0.25, option_type="call",
             qty=10, multiplier=100, label="AAPL 105C"),
    Contract(S=100, K=95,  t=0.5, r=0.04, sigma=0.30, option_type="put",
             qty=-5, multiplier=100, label="AAPL 95P"),
])
print(book.net.delta, book.net.gamma, book.net.vega, book.net.market_value)
for pos in book.positions:
    print(pos.contract.label, pos.position_delta)
```

`qty` is signed (short = negative) and `multiplier` scales to notional
(e.g. 100 for US equity options). Net Greeks are position-scaled sums.

## Implied forward and dividend

Recover the forward price and discount factor directly from a call/put chain
using put-call parity — no volatility assumption:

```python
from quantforge import implied_forward

res = implied_forward(strikes, calls, puts, t=1.0, spot=100.0)
print(res.forward, res.discount_factor)
print(res.implied_rate, res.implied_div_yield)
```

The points `(K, C - K)` vs `P` are linear under parity, so a least-squares fit
returns both the forward (slope) and discount factor (intercept) at once, and
the dividend yield follows from `F = S * exp((r - q) t)`.

## Realized volatility

Estimate historical vol from a price series — close-to-close, EWMA, and the
efficient range-based estimators (Parkinson, Garman-Klass, Rogers-Satchell,
Yang-Zhang). All annualized:

```python
from quantforge import close_to_close, yang_zhang, vol_report

closes = [...]                      # daily closes
print(close_to_close(closes))       # classic
print(ewma_vol(closes, lam=0.94))   # RiskMetrics

# Range-based estimators need OHLC bars and are much lower-variance.
rep = vol_report(opens, highs, lows, closes)
print(rep.parkinson, rep.garman_klass, rep.rogers_satchell, rep.yang_zhang)
```

Yang-Zhang is drift-independent and handles overnight gaps; it is the default
choice when you have clean OHLC data. Compare any of these against the implied
vol from `implied_volatility` to trade realized-vs-implied.

## Portfolio risk (VaR / Expected Shortfall)

Three estimators over a priced book — parametric delta-gamma (Cornish-Fisher),
historical, and full-repricing Monte Carlo:

```python
from quantforge import Contract, price_book, parametric_var, montecarlo_var

positions = [
    Contract(S=100, K=100, t=0.25, r=0.02, sigma=0.3, option_type="call", qty=-1),
    Contract(S=100, K=100, t=0.25, r=0.02, sigma=0.3, option_type="put",  qty=-1),
]
book = price_book(positions)

# 99% 1-day VaR from the book's net delta/gamma (captures convexity).
r = parametric_var(book, sigma_annual=0.3, spot=100, confidence=0.99, horizon_days=1)
print(r.var, r.expected_shortfall, r.method)

# Full-repricing MC VaR: no delta-gamma approximation.
r = montecarlo_var(positions, sigma_annual=0.3, confidence=0.99, horizon_days=5,
                   n_paths=20_000, seed=1)
print(r.var, r.expected_shortfall)
```

VaR is a positive loss number; Expected Shortfall (CVaR) is the mean loss
beyond it. The parametric/historical estimators assume a single underlying;
the MC estimator reprices every position exactly.

## Monte Carlo (with variance reduction)

For payoffs without a closed form. The engine uses the standard-library RNG,
antithetic variates, and a geometric-Asian control variate — and returns a
standard error with every price:

```python
from quantforge import european_mc, arithmetic_asian_mc

# European call: converges to the exact BSM price.
res = european_mc(S=100, K=100, t=1, r=0.05, sigma=0.2, seed=42)
print(res.price, "+/-", res.std_error, res.confidence_interval())

# Arithmetic-average Asian, with the geometric Asian as a control variate.
# control_variate=True cuts the standard error by ~10x for free.
res = arithmetic_asian_mc(S=100, K=100, t=1, r=0.05, sigma=0.3,
                          n_steps=50, n_paths=50_000, control_variate=True, seed=1)
print(res.price, "+/-", res.std_error)
```

## Exotic options (closed form)

Analytic prices for binaries, single barriers, and geometric Asians:

```python
from quantforge import (
    cash_or_nothing, asset_or_nothing, barrier_option, geometric_asian, Barrier,
)

# Digital: pays $10 if the call finishes in the money.
cash_or_nothing(S=100, K=105, t=1, r=0.05, sigma=0.2, option_type="call", cash=10)

# Down-and-out call with a $3 knock-out rebate (Reiner-Rubinstein).
barrier_option(S=100, K=90, H=95, t=0.5, r=0.08, sigma=0.25,
               option_type="call", barrier=Barrier.DOWN_OUT, b=0.04, rebate=3)

# Geometric-average Asian call (Kemna-Vorst closed form).
geometric_asian(S=100, K=100, t=1, r=0.05, sigma=0.3, option_type="call")
```

Barrier kinds: `Barrier.DOWN_IN`, `DOWN_OUT`, `UP_IN`, `UP_OUT`. In/out parity
(`in + out = vanilla`) holds exactly and is enforced by the tests.

## Term-structure surface (calendar-arbitrage aware)

Stitch per-expiry SVI smiles into a full surface, interpolate vol at any
`(log-moneyness, expiry)`, and check for calendar arbitrage:

```python
from quantforge import VolSurface

# Per-expiry quotes: (t, log-moneyness points, total variances w = sigma^2 t).
quotes = [
    (0.5, ks_6m, w_6m),
    (1.0, ks_1y, w_1y),
    (2.0, ks_2y, w_2y),
]
surf = VolSurface.fit(quotes)

surf.implied_vol(k=0.05, t=0.75)          # interpolated vol between expiries
surf.is_calendar_arbitrage_free()         # True if variance rises with maturity
for v in surf.calendar_arbitrage():
    print(v.t_short, v.t_long, v.k, v.w_short, v.w_long)   # any crossing curves
```

Total variance is interpolated linearly in maturity (the standard
no-arbitrage-friendly scheme) and the calendar check enforces that
`w(k, t)` is non-decreasing in `t` at every strike.

## SABR stochastic-vol smile

The market-standard SABR model via Hagan's implied-vol expansion, with
calibration of (alpha, rho, nu) at a fixed beta:

```python
from quantforge import sabr_vol, calibrate_sabr

# Evaluate the smile.
sabr_vol(F=100, K=90, t=1.0, alpha=0.2, beta=0.5, rho=-0.3, nu=0.4)

# Calibrate to a market smile.
strikes = [80, 90, 100, 110, 120]
market  = [0.26, 0.235, 0.22, 0.225, 0.24]
params, rmse = calibrate_sabr(F=100, t=0.5, strikes=strikes,
                              market_vols=market, beta=0.5)
print(params, rmse)
```

SABR is the standard for interest-rate and FX smiles; SVI (below) is the
common equity-index parametrization. Both interpolate/extrapolate a smile and
plug into `implied_volatility`/`price` for consistent surface pricing.

## Volatility surface (SVI)

Fit Gatheral's raw SVI smile to market quotes with a built-in Nelder-Mead
calibrator (no SciPy):

```python
from quantforge import calibrate_svi, SVIParams

# Observed smile: log-moneyness k = log(K/F), total variance w = sigma^2 * t.
ks  = [-0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3]
t   = 0.5
vols = [0.28, 0.24, 0.21, 0.20, 0.205, 0.22, 0.245]
tv  = [(v * v) * t for v in vols]

params, rmse = calibrate_svi(ks, tv)
print(params, "rmse:", rmse)

# Interpolate/extrapolate a vol anywhere on the smile.
print(params.implied_vol(k=0.05, t=t))

# Static no-arbitrage wing check (Lee's slope bound).
print(params.is_arbitrage_free_wings())
```

## Command line

Installing the package exposes a `quantforge` CLI:

```bash
quantforge price   -S 100 -K 105 -t 0.5 -r 0.04 --sigma 0.25 --type call
quantforge greeks  -S 100 -K 105 -t 0.5 -r 0.04 --sigma 0.25
quantforge iv      -S 100 -K 105 -t 0.5 -r 0.04 --price 6.12 --type call
quantforge american -S 100 -K 100 -t 1 -r 0.05 --sigma 0.2 --type put -b 0.02
```

## Model coverage

| Instrument            | Set the carry `b` to | Function            |
|-----------------------|----------------------|---------------------|
| Non-dividend stock    | `b = r` (default)    | `price`, `greeks`   |
| Continuous dividend q | `b = r - q`          | `price`, `greeks`   |
| Option on a future    | `b = 0` (Black-76)   | `price`, `greeks`   |
| FX option             | `b = r - r_foreign`  | `price`, `greeks`   |
| American exercise     | any of the above     | `american_price`    |

## Vectorized fast path (optional NumPy)

The core is zero-dependency. If you have NumPy and need to price whole chains,
`quantforge.vectorized` runs the same formulas over arrays:

```python
import numpy as np
from quantforge import price_array, greeks_array

S = np.full(10_000, 100.0)
K = np.linspace(80, 120, 10_000)          # a strike grid
p = price_array(S, K, t=0.5, r=0.04, sigma=0.25, option_type="call")
g = greeks_array(S, K, 0.5, 0.04, 0.25)   # dict of price/delta/gamma/vega arrays
```

Scalars broadcast against arrays. Results match the scalar engine to machine
precision (1e-12) and run ~6x faster on large batches. Install with:

```bash
pip install "quantforge[fast]"
```

## Performance

Pure-Python, but fast enough for interactive risk work. On a typical laptop
(single core, CPython 3.12):

| Operation                     | Throughput        |
|-------------------------------|-------------------|
| Call price (BSM)              | ~1.4M / sec       |
| Full Greeks (6 outputs)       | ~245K / sec       |
| Implied vol (Newton+bisect)   | ~110K / sec       |
| American (200-step tree)      | ~220 / sec        |

Implied-vol round-trips to a max price error of `1e-8`. Reproduce with:

```bash
python benchmarks/bench.py --n 100000
python benchmarks/bench.py --vollib   # accuracy check vs py_vollib, if installed
```

## Testing

```bash
pip install -e ".[test]"
pytest -q
```

The suite pins prices against Hull's textbook reference values, enforces
put-call parity, checks every Greek against central finite differences, and
round-trips implied volatility.

## License

MIT.
