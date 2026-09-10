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
- **Every Greek analytic** — delta, gamma, vega, theta, rho, epsilon (dividend
  rho), plus second-order (vanna/vomma/charm/…), verified against finite
  differences to 1e-4 or better.
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

# Invert a whole chain to a smile in one call: [(log-moneyness, vol), ...].
from quantforge import implied_vol_smile
print(implied_vol_smile(strikes, prices, S=100, t=0.5, r=0.04))

# American put with a 3% dividend yield (carry b = r - q).
print(american_price(S=100, K=100, t=1.0, r=0.05, sigma=0.2,
                     option_type="put", b=0.05 - 0.03, steps=500))

# Or the closed-form Bjerksund-Stensland (2002) approximation: ~500x faster
# than a 500-step tree, accurate to a few cents.
from quantforge import bjerksund_stensland, bjerksund_stensland_greeks
print(bjerksund_stensland(S=100, K=100, t=1.0, r=0.05, sigma=0.2,
                          option_type="put", b=0.02))
# American Greeks (delta/gamma/vega/theta/rho) by finite differences:
print(bjerksund_stensland_greeks(S=100, K=100, t=1.0, r=0.05, sigma=0.2,
                                 option_type="put", b=0.02))
# Split into European value + early-exercise premium:
from quantforge import early_exercise_premium
print(early_exercise_premium(S=90, K=100, t=1.0, r=0.05, sigma=0.3,
                             option_type="put"))

# Or a trinomial lattice (smoother convergence), with Richardson extrapolation
# for extra accuracy at low step counts.
from quantforge import trinomial_price, richardson_american
print(trinomial_price(S=100, K=100, t=1.0, r=0.05, sigma=0.2,
                      option_type="put", steps=400))
print(richardson_american(S=100, K=100, t=1.0, r=0.05, sigma=0.2,
                          option_type="put", steps=50))

# Or Bermudan/American by Longstaff-Schwartz least-squares Monte Carlo.
from quantforge import bermudan_lsm
print(bermudan_lsm(S=100, K=100, t=1.0, r=0.05, sigma=0.2,
                   option_type="put", n_steps=50, n_paths=40_000, seed=1))

# Perpetual (no-expiry) American, exact closed form (Merton).
from quantforge import perpetual_american
print(perpetual_american(S=100, K=100, r=0.08, sigma=0.3, option_type="put"))
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

## Smile-aware delta (sticky-strike / sticky-delta)

The BS delta assumes vol is fixed as spot moves. When the smile rides with the
underlying (sticky-delta / sticky-moneyness), the effective hedge ratio picks
up a skew term:

```python
from quantforge import smile_delta, smile_delta_from_smile, StickyRule

# Explicit skew slope dsigma/dk at this strike:
smile_delta(S=100, K=105, t=0.5, r=0.04, sigma=0.25, dsigma_dk=-0.4,
            option_type="call", sticky=StickyRule.DELTA)

# Or straight from a smile function sigma(K):
smile = lambda K: 0.25 - 0.2 * (K / 100 - 1)
smile_delta_from_smile(S=100, K=105, t=0.5, r=0.04, smile_fn=smile,
                       option_type="call")
```

Sticky-strike returns the plain BS delta; sticky-delta adds
`-vega * (dsigma/dk) / S`, so a downward equity skew raises the call delta.

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

## Option strategies

Named multi-leg builders return a `Book`, so net price and Greeks come from the
same engine, plus a payoff diagram and break-even solver:

```python
from quantforge import straddle, iron_condor, payoff_profile, break_evens

book = straddle(S=100, K=100, t=0.5, r=0.04, sigma=0.25)
print(book.net.market_value, book.net.delta, book.net.vega)   # net premium/Greeks
print(break_evens(book, lo=50, hi=150))                       # [K-prem, K+prem]

condor = iron_condor(S=100, K_put_long=80, K_put_short=90,
                     K_call_short=110, K_call_long=120,
                     t=0.5, r=0.04, sigma=0.25)
print(payoff_profile(condor, spots=[70, 90, 100, 110, 130]))
```

Also `vertical_spread`, `strangle`, `risk_reversal`, `butterfly`,
`ratio_spread` (net short options), and `backspread` (net long options).

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

`book_second_order` aggregates the net second-order Greeks (vanna, vomma,
charm, veta, speed, zomma, color) the same way:

```python
from quantforge import book_second_order

bk = book_second_order(positions)
print(bk.vanna, bk.vomma, bk.charm)
```

See where the book's vol risk sits along the curve with `vega_buckets`:

```python
from quantforge import vega_buckets

vb = vega_buckets(positions, edges=(0.25, 0.5, 1.0, 2.0, 5.0))
print(vb.buckets)   # net vega per maturity bucket; sums to the net book vega
```

Explain a position's realized P&L over a move via its Greeks:

```python
from quantforge import attribute_pnl

a = attribute_pnl(S=100, K=100, t=0.5, r=0.05, sigma=0.25,
                  dS=3, dsigma=0.02, dt=1/252)
print(a.delta_pnl, a.gamma_pnl, a.vega_pnl, a.theta_pnl, a.unexplained)
```

`total = explained + unexplained`; the residual flags what the first/second-
order Greeks miss.

Carry-roll (roll-down) P&L — "if nothing moves, what do I earn?":

```python
from quantforge import carry_roll_pnl

carry_roll_pnl(S=100, K=100, t=1.0, r=0.05, sigma=0.2, horizon=0.25).theta_roll
```

Split the book's time decay into gamma rent vs financing carry:

```python
from quantforge import theta_carry_report

tc = theta_carry_report(positions)
print(tc.theta, tc.gamma_rent, tc.residual)   # theta = gamma_rent + residual
```

Size a hedge to a target Greek off the book's net exposures:

```python
from quantforge import delta_hedge_shares, vega_neutral_quantity

delta_hedge_shares(book)                                        # shares to zero delta
vega_neutral_quantity(book, S=100, K=110, t=0.5, r=0.04, sigma=0.25)  # option units to zero vega
```

Kelly-criterion bet/leverage sizing:

```python
from quantforge import kelly_fraction_binary, kelly_fraction_continuous

kelly_fraction_binary(win_prob=0.6, win_payoff=1.0)     # 0.2 of bankroll
kelly_fraction_continuous(0.08, 0.04, fraction=0.5)     # half-Kelly leverage
```

## Variance / volatility swaps

Model-free fair strike of a variance swap from an option strip (the log-contract
replication behind the VIX):

```python
from quantforge import variance_swap_strike, volatility_swap_strike

kvar = variance_swap_strike(S0=100, t=1.0, r=0.05,
                            put_strikes=puts, put_prices=pp,
                            call_strikes=calls, call_prices=cp)
kvol = volatility_swap_strike(S0=100, t=1.0, r=0.05,
                              put_strikes=puts, put_prices=pp,
                              call_strikes=calls, call_prices=cp)
```

The strike is `(2/t)` times the log-contract value replicated by OTM options
weighted `1/K^2`. For a flat-vol chain `kvar -> sigma^2`. The vol-swap value is
the `sqrt` proxy (an upper bound before the convexity adjustment).

## Risk-neutral density (Breeden-Litzenberger)

Extract the market-implied probability distribution of the underlying at expiry
from a call-price curve, then reprice any payoff against it:

```python
from quantforge import risk_neutral_density, price_from_density, density_total_mass

mids, pdf = risk_neutral_density(strikes, calls, t=1.0, r=0.05)
print(density_total_mass(strikes, calls, 1.0, 0.05))   # ~ 1.0 for a clean curve

# Price an arbitrary European payoff by integrating against the density.
digital = price_from_density(strikes, calls, 1.0, 0.05,
                             payoff=lambda ST: 1.0 if ST > 100 else 0.0)
```

The density is the discounted second derivative of the call curve in strike;
`risk_neutral_cdf` gives the CDF from the first derivative.

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

`dividend_curve` bootstraps the implied dividend-yield term structure from a
multi-expiry chain:

```python
from quantforge import dividend_curve

curve = dividend_curve(chain_by_expiry, spot=100.0)   # [(t, ForwardResult), ...]
for t, res in curve:
    print(t, res.implied_div_yield)
```

## Index implied correlation

The dispersion-trading measure: back out the common correlation an index vol
implies given its members' weights and vols:

```python
from quantforge import implied_correlation, index_vol_from_correlation

weights = [0.4, 0.35, 0.25]
vols    = [0.25, 0.30, 0.20]
implied_correlation(weights, vols, index_vol=0.19)   # rho consistent with the quote
index_vol_from_correlation(weights, vols, rho=0.3)   # forward map
```

`rho=0` gives the fully-diversified `dispersion_basket_vol`; `rho=1` gives the
weighted-average vol.

`correlation_term_structure` implies the correlation at each expiry from index
and member vol curves:

```python
from quantforge import correlation_term_structure

correlation_term_structure(weights, member_vol_curves, index_vol_curve,
                           expiries=[0.5, 1.0, 2.0])   # [(expiry, rho), ...]
```

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
choice when you have clean OHLC data.

`vol_cone` shows the distribution of rolling realized vol per horizon — the
standard "is current vol high or low" tool:

```python
from quantforge import vol_cone

for pt in vol_cone(closes, windows=[5, 21, 63, 126]):
    print(pt.window, pt.minimum, pt.median, pt.maximum, pt.current)
```

`fit_garch` estimates a GARCH(1,1) model and `garch_forecast` projects the vol
forward, mean-reverting to the long-run level:

```python
from quantforge import fit_garch, garch_forecast

p = fit_garch(returns)
garch_forecast(p, last_return=returns[-1], last_variance=h, horizon=21)
``` Compare any of these against the implied
vol from `implied_volatility` to trade realized-vs-implied.

## Delta-hedge P&L simulator

Monte Carlo a discretely delta-hedged short option and see the hedging-error
distribution — including what happens when you hedge at the wrong vol:

```python
from quantforge import simulate_delta_hedge

# Rehedge 50 times over the option's life; error std shrinks like 1/sqrt(n).
res = simulate_delta_hedge(S=100, K=100, t=1.0, r=0.05, sigma=0.2,
                           n_steps=50, n_paths=20_000, seed=1)
print(res.mean_pnl, res.std_pnl)   # mean ~ 0, std = discrete-hedging error

# Hedge at 15% while the market realizes 30%: short gamma bleeds.
res = simulate_delta_hedge(S=100, K=100, t=1.0, r=0.03, sigma=0.15,
                           hedge_vol=0.15, real_vol=0.30, n_steps=50, seed=11)
print(res.mean_pnl)   # negative
```

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

## Two-asset options

Exchange (Margrabe), spread (Kirk), and basket options on two correlated
assets:

```python
from quantforge import exchange_option, spread_option, basket_option

# Right to swap asset 2 for asset 1: max(S1 - S2, 0). Exact, rate-independent.
exchange_option(S1=100, S2=95, t=1.0, sigma1=0.2, sigma2=0.25, rho=0.3)

# Its Greeks (two deltas, own/cross gammas, correlation vega) by finite diff:
from quantforge import exchange_greeks
exchange_greeks(S1=100, S2=95, t=1.0, sigma1=0.2, sigma2=0.25, rho=0.3)

# Spread: max(S1 - S2 - K, 0) via Kirk's approximation.
spread_option(S1=100, S2=95, K=5, t=1.0, r=0.05, sigma1=0.2, sigma2=0.25, rho=0.5)

# Basket call on w1 S1 + w2 S2 (lognormal moment-match).
basket_option(spots=(100, 100), weights=(0.5, 0.5), K=100, t=1.0, r=0.05,
              sigmas=(0.2, 0.3), corr=0.4)

# Rainbow: option on the best (max) or worst (min) of two assets (Monte Carlo).
from quantforge import best_of_call, worst_of_call
best_of_call(S1=100, S2=100, K=100, t=1.0, r=0.05, sigma1=0.2, sigma2=0.25, rho=0.4)
```

Best-of + worst-of equals `call(S1) + call(S2)` (Stulz identity).

## Forward-start and cliquet options

Options whose strike is fixed at a future date (as a multiple of the then-spot),
and strips of them (cliquets/ratchets):

```python
from quantforge import forward_start_price, cliquet_price

# ATM-at-reset call: strike set at t=0.5 to alpha*S, expiring at t=1.0.
forward_start_price(S=100, t_start=0.5, t_expiry=1.0, r=0.05, sigma=0.25,
                    alpha=1.0, option_type="call")

# A quarterly-reset cliquet over one year.
cliquet_price(S=100, reset_times=[0.25, 0.5, 0.75, 1.0], r=0.05, sigma=0.2,
              alpha=1.0, option_type="call")
```

Under BSM the forward-start value scales with the current spot and is
independent of the absolute future strike (Rubinstein), so `alpha` (moneyness)
is the only strike input.

A simple chooser (pick call or put at a future date) has a closed form:

```python
from quantforge import chooser_option

chooser_option(S=100, K=100, t_choose=0.5, T=1.0, r=0.05, sigma=0.25)
```

A compound option (option on an option, Geske) prices all four kinds:

```python
from quantforge import compound_option

compound_option(S=100, K1=5, K2=100, t1=0.5, t2=1.0, r=0.05, sigma=0.25,
                kind="call-on-call")
```

An autocallable note (early redemption + coupon at observation dates, downside
below a protection barrier) is priced by Monte Carlo:

```python
from quantforge import autocallable_mc

autocallable_mc(S=100, t=3.0, r=0.03, sigma=0.25,
                observation_times=[1.0, 2.0, 3.0], autocall_barrier=100,
                coupon=0.08, protection_barrier=70)
```

A capped cliquet (ratchet note) — periodic returns clipped locally and the sum
clipped globally — is priced by Monte Carlo:

```python
from quantforge import capped_cliquet_mc

capped_cliquet_mc(S=100, t=1.0, r=0.05, sigma=0.3,
                  reset_times=[0.25, 0.5, 0.75, 1.0],
                  local_cap=0.05, global_cap=0.15)
```

## Quasi-Monte Carlo

Deterministic low-discrepancy (Halton) integration converges several times
faster than pseudo-random Monte Carlo for a European payoff:

```python
from quantforge import european_qmc, halton

european_qmc(S=100, K=100, t=1.0, r=0.05, sigma=0.2, n_points=8192)
halton(index=0, dim=2)   # a low-discrepancy point in the unit square
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

# Arithmetic-average Asian call (Turnbull-Wakeman moment matching).
from quantforge import arithmetic_asian
arithmetic_asian(S=100, K=100, t=1, r=0.05, sigma=0.3, option_type="call")

# Asian Greeks (delta/gamma/vega/theta) for either average, by finite diff.
from quantforge import asian_greeks
asian_greeks(S=100, K=100, t=1, r=0.05, sigma=0.3, average="geometric")

# Average-strike Asian (strike = realized average) by Monte Carlo.
from quantforge import average_strike_asian_mc
average_strike_asian_mc(S=100, t=1, r=0.05, sigma=0.3, option_type="call")
```

Barrier kinds: `Barrier.DOWN_IN`, `DOWN_OUT`, `UP_IN`, `UP_OUT`. In/out parity
(`in + out = vanilla`) holds exactly and is enforced by the tests.

`barrier_greeks` returns delta/gamma/vega/theta of a barrier option by finite
differences on the closed form:

```python
from quantforge import barrier_greeks, Barrier

barrier_greeks(S=100, K=100, H=90, t=0.5, r=0.05, sigma=0.25,
               option_type="call", barrier=Barrier.DOWN_OUT)
```

Digitals have unbounded pin risk at the strike, so desks super-replicate them
with a tight vanilla spread:

```python
from quantforge import digital_call_overhedge

oh = digital_call_overhedge(S=100, K=100, t=0.25, r=0.05, sigma=0.2,
                            cash=1.0, width=2)
print(oh.cost, oh.digital_value, oh.cushion)   # spread cost >= fair value
```

The spread payoff dominates the digital everywhere and its cost converges to
the fair digital value as `width -> 0`; the cushion is the pin-risk buffer.

A double-knockout (corridor) option pays the vanilla payoff only if the spot
stays inside two barriers for the whole path:

```python
from quantforge import double_knockout_mc

double_knockout_mc(S=100, K=100, t=1.0, r=0.05, sigma=0.25,
                   lower=90, upper=115, option_type="call")
```

A Parisian barrier activates only after the spot stays past the level for a
consecutive window (robust to brief spikes), priced by Monte Carlo:

```python
from quantforge import parisian_barrier_mc

parisian_barrier_mc(S=100, K=100, H=90, t=1.0, r=0.05, sigma=0.25,
                    window=0.1, option_type="call", barrier="down-out")
```

A barrier-contingent digital pays cash only if it finishes ITM *and* the
barrier condition holds (Monte Carlo):

```python
from quantforge import barrier_digital_mc

# Pays $1 if S_T > 100 and 120 was touched over the path.
barrier_digital_mc(S=100, K=100, H=120, t=1.0, r=0.05, sigma=0.2,
                   option_type="call", barrier="up-in", cash=1.0)
```

One-touch / no-touch binaries pay a fixed cash on (or against) a barrier being
hit:

```python
from quantforge import one_touch, no_touch

# Pays $1 the moment spot touches 120 (FX pay-at-hit convention).
one_touch(S=100, H=120, t=1.0, r=0.05, sigma=0.2, cash=1.0, payoff_at_hit=True)

# Pays $1 at expiry only if 80 is never touched.
no_touch(S=100, H=80, t=1.0, r=0.05, sigma=0.2, cash=1.0)

# Standalone barrier rebate: knock-out (on breach) or knock-in (if never hit).
from quantforge import barrier_rebate
barrier_rebate(S=100, H=120, t=1.0, r=0.05, sigma=0.2, knock="out", cash=1.0)
```

Lookbacks (against the realized path extreme) have closed forms too:

```python
from quantforge import floating_strike_lookback, fixed_strike_lookback

# Floating strike (buy at the low / sell at the high): payoff S_T - S_min.
floating_strike_lookback(S=100, t=1.0, r=0.05, sigma=0.3, option_type="call")

# Fixed strike on the realized maximum: payoff max(S_max - K, 0).
fixed_strike_lookback(S=100, K=100, t=1.0, r=0.05, sigma=0.3, option_type="call")
```

Pass `s_extreme` (the running min/max observed so far) to price a seasoned
lookback; it defaults to the current spot at inception.

Gap options (separate trigger and payoff strikes) and power options (payoff on
`S^p`) have closed forms too:

```python
from quantforge import gap_option, power_option

gap_option(S=100, K_trigger=90, K_payoff=110, t=1.0, r=0.05, sigma=0.25,
           option_type="call")
power_option(S=100, K=10000, t=1.0, r=0.05, sigma=0.2, power=2.0,
             option_type="call")
```

## Local volatility (Dupire)

Extract the Dupire local-volatility function from an implied-vol or call-price
surface:

```python
from quantforge import local_vol_from_implied, dupire_local_vol

# From an implied-vol surface sigma(K, T):
smile = lambda K, T: 0.25 - 0.1 * (K / 100 - 1.0)
local_vol_from_implied(smile, S=100, K=105, T=1.0, r=0.03, q=0.0)

# Or directly from a call-price surface C(K, T):
dupire_local_vol(call_fn, K=100, T=1.0, r=0.03, q=0.0)
```

A flat implied surface returns a constant local vol; a pure term structure
returns the analytic `sqrt(dw/dT)`.

`sabr_local_vol` feeds a SABR smile straight into Dupire at one expiry:

```python
from quantforge import sabr_local_vol

sabr_local_vol(S=100, K=90, T=1.0, r=0.05,
               alpha=0.25, beta=0.5, rho=-0.4, nu=0.5)
```

`local_vol_mc` prices any European payoff under a local-vol surface
`sigma_loc(S, t)`:

```python
from quantforge import local_vol_mc

local_vol_mc(S=100, K=100, t=1.0, r=0.05,
             local_vol_fn=lambda S, tau: 0.2 * (100 / S) ** 0.5,
             option_type="call")
```

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
# ...or fit and repair any calendar arbitrage so w is non-decreasing in t:
surf = VolSurface.fit_arbitrage_free(quotes)

surf.implied_vol(k=0.05, t=0.75)          # interpolated vol between expiries
surf.is_calendar_arbitrage_free()         # True if variance rises with maturity
for v in surf.calendar_arbitrage():
    print(v.t_short, v.t_long, v.k, v.w_short, v.w_long)   # any crossing curves
```

Total variance is interpolated linearly in maturity (the standard
no-arbitrage-friendly scheme) and the calendar check enforces that
`w(k, t)` is non-decreasing in `t` at every strike.

The forward vol between two maturities (the vol a forward-starting option
sees) comes straight off the surface:

```python
surf.forward_vol(k=0.0, t1=1.0, t2=2.0)   # sqrt of the forward variance
```

Export a grid of vols for a heatmap/surface plot (by log-moneyness or strike):

```python
surf.vol_grid(ks=[-0.2, 0.0, 0.2], ts=[0.5, 1.0, 2.0])
surf.strike_vol_grid(spot=100, strikes=[90, 100, 110], ts=[0.5, 1.0])
```

## Skew/kurtosis-adjusted pricing (Corrado-Su)

Add the first skewness and kurtosis corrections to Black-Scholes via a
Gram-Charlier expansion, and estimate the moments from history:

```python
from quantforge import (
    corrado_su_price, realized_skewness, realized_excess_kurtosis,
)

corrado_su_price(S=100, K=110, t=1.0, r=0.05, sigma=0.2,
                 skew=-0.5, excess_kurt=2.0, option_type="call")

skew = realized_skewness(returns)
xk = realized_excess_kurtosis(returns)
```

`skew=0, excess_kurt=0` recovers Black-Scholes exactly.

## CEV (constant elasticity of variance)

Local vol that scales with the spot level (the leverage effect), priced in
closed form via the noncentral chi-square distribution — implemented from
scratch, no SciPy:

```python
from quantforge import cev_price

# beta in [0, 1); lower beta => stronger downside skew. beta -> 1 is BSM.
cev_price(S=100, K=90, t=1.0, r=0.05, sigma=0.2, beta=0.5, option_type="put")
```

## Variance-Gamma (pure-jump)

Brownian motion evaluated at a random gamma business time — a pure-jump model
with skew (`theta`) and kurtosis (`nu`) controls, priced via its characteristic
function:

```python
from quantforge import variance_gamma_price

variance_gamma_price(S=100, K=100, t=1.0, r=0.05, sigma=0.2,
                     nu=0.3, theta=-0.1, option_type="call")
```

`nu -> 0` recovers Black-Scholes; a negative `theta` gives the equity left
skew.

## Bachelier (normal) model

For rates and spread options where the forward can go negative and a lognormal
model breaks down. `sigma` is the normal (absolute) volatility:

```python
from quantforge import bachelier_price, bachelier_implied_vol

# Works with negative forwards/strikes (e.g. a rate spread).
bachelier_price(F=-0.5, K=-1.0, t=1.0, r=0.0, sigma=2.0, option_type="call")

# Back out the normal vol from a market price.
bachelier_implied_vol(target_price=6.12, F=100, K=100, t=1.0, r=0.02)
```

Includes analytic `bachelier_delta`, `bachelier_gamma`, and `bachelier_vega`.

Displaced diffusion (shifted lognormal) interpolates between Black-Scholes and
Bachelier and permits negative strikes:

```python
from quantforge import displaced_diffusion_price

displaced_diffusion_price(S=100, K=100, t=1.0, r=0.05, sigma=0.2, shift=50,
                          option_type="call")
```

`shift=0` is Black-Scholes; a larger shift moves toward normal-model behavior.

Interest-rate caps/floors/collars build on it as Bachelier caplet strips (so
they handle negative rates):

```python
from quantforge import CapletPeriod, cap_price, floor_price, collar_price
import math

periods = [CapletPeriod(forward=0.03 + 0.002*i, expiry=float(i), accrual=1.0,
                        discount=math.exp(-0.03*i), sigma_n=0.01)
           for i in range(1, 5)]
cap_price(periods, strike=0.04)
floor_price(periods, strike=0.03)
collar_price(periods, cap_strike=0.045, floor_strike=0.03)
```

Cap(K) - Floor(K) equals the discounted swap PV, which the tests enforce.

Swaptions price as the annuity times a Bachelier option on the forward swap
rate:

```python
from quantforge import swaption_price

# Payer swaption (call on the swap rate); receiver is payer=False.
swaption_price(swap_rate=0.03, strike=0.035, expiry=2.0, sigma_n=0.01,
               periods=periods, payer=True)
```

Key-rate DV01 buckets the interest-rate sensitivity of any book expressed as a
function of the zero curve:

```python
from quantforge import key_rate_dv01

price = lambda curve: sum(cf * 2.718281828 ** (-curve[t] * t)
                          for t, cf in {1.0: 5, 2.0: 105}.items())
kr = key_rate_dv01(price, base_curve={1.0: 0.02, 2.0: 0.03})
print(kr.buckets, kr.parallel)   # per-tenor DV01s sum to the parallel DV01
```

## Merton jump-diffusion

Adds lognormal jumps to the diffusion; priced as a Poisson-weighted series of
Black-Scholes values:

```python
from quantforge import merton_jump_price

merton_jump_price(S=100, K=100, t=1.0, r=0.05, sigma=0.2,
                  lam=1.0, mu_j=-0.1, sigma_j=0.15, option_type="call")
```

`lam` is the jump intensity (jumps per year), `mu_j`/`sigma_j` the log-jump
mean/std. At `lam=0` it is exactly Black-Scholes; the drift is compensated so
put-call parity holds.

## Heston stochastic volatility

Full stochastic-variance pricing via the characteristic function, integrated
with a built-in Gauss-Legendre rule (no SciPy):

```python
from quantforge import heston_price

heston_price(S=100, K=100, t=1.0, r=0.0,
             v0=0.04, kappa=2.0, theta=0.04, xi=0.3, rho=-0.7,
             option_type="call")
```

`v0` is the initial variance, `kappa`/`theta` the mean-reversion speed/level,
`xi` the vol-of-vol, and `rho` the spot/variance correlation. As `xi -> 0` the
price collapses to Black-Scholes; puts follow from put-call parity.

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

For a model-free single smile, `SmileSpline` fits a natural cubic spline
through the quoted vols with flat extrapolation:

```python
from quantforge import SmileSpline

sm = SmileSpline(strikes=[80, 90, 100, 110, 120],
                 vols=[0.28, 0.24, 0.22, 0.225, 0.24])
sm.vol(95)    # smoothly interpolated
sm.vol(200)   # flat beyond the quoted range
```

For an FX smile from the three market quotes, `VannaVolgaSmile` builds the
25-delta pillars and interpolates:

```python
from quantforge import VannaVolgaSmile

sm = VannaVolgaSmile(S=1.2, t=1.0, r_dom=0.02, r_for=0.01,
                     atm=0.10, rr=0.02, bf=0.005)
sm.vol(1.25)      # smile vol at a strike
sm.pillars()      # [(K_25P, v), (K_ATM, v), (K_25C, v)]
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

# Butterfly (density) no-arbitrage check via the Gatheral-Jacquier g-function.
from quantforge import svi_is_butterfly_free, svi_butterfly_arbitrage
print(svi_is_butterfly_free(params))
print(svi_butterfly_arbitrage(params))   # log-moneyness points that violate
```

## Command line

Installing the package exposes a `quantforge` CLI:

```bash
quantforge price   -S 100 -K 105 -t 0.5 -r 0.04 --sigma 0.25 --type call
quantforge greeks  -S 100 -K 105 -t 0.5 -r 0.04 --sigma 0.25
quantforge iv      -S 100 -K 105 -t 0.5 -r 0.04 --price 6.12 --type call
quantforge american -S 100 -K 100 -t 1 -r 0.05 --sigma 0.2 --type put -b 0.02
```

## Quanto options

A foreign-asset option settled in domestic currency at a fixed FX rate — the
quanto adjustment enters through the carry:

```python
from quantforge import quanto_option

quanto_option(S=100, K=100, t=1.0, r_domestic=0.03, r_foreign=0.05,
              sigma_asset=0.25, sigma_fx=0.1, rho=0.5, option_type="call")
```

The carry shifts by `-rho * sigma_asset * sigma_fx`; `rho=0` removes the
adjustment.

A composite (compo) option instead converts at the *floating* FX, so its vol
combines the asset and FX vols:

```python
from quantforge import compo_option

compo_option(S=100, K=100, t=1.0, r_domestic=0.05, r_foreign=0.05,
             sigma_asset=0.2, sigma_fx=0.1, rho=0.5, option_type="call")
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
