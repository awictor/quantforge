# QuantForge

![CI](https://github.com/awictor/quantforge/actions/workflows/ci.yml/badge.svg)
![Python](https://img.shields.io/badge/python-3.8%2B-blue)
![License](https://img.shields.io/badge/license-MIT-green)

Fast, **dependency-free** options pricing and risk engine in pure Python.

QuantForge prices European and American options, computes the full set of
analytic Greeks, and solves for implied volatility — with zero third-party
dependencies. It drops into any Python 3.8+ runtime (serverless, embedded,
notebooks, trading bots) without compiling NumPy or SciPy.

## Contents

<!-- TOC -->
- [Contents](#contents)
- [Why](#why)
- [Install](#install)
- [Quick start](#quick-start)
- [Scenario / stress grid](#scenario--stress-grid)
- [Smile-aware delta (sticky-strike / sticky-delta)](#smile-aware-delta-sticky-strike--sticky-delta)
- [Second-order Greeks](#second-order-greeks)
- [Option strategies](#option-strategies)
- [Batch pricing and portfolio risk](#batch-pricing-and-portfolio-risk)
- [Variance / volatility swaps](#variance--volatility-swaps)
- [Risk-neutral density (Breeden-Litzenberger)](#risk-neutral-density-breeden-litzenberger)
- [Implied forward and dividend](#implied-forward-and-dividend)
- [Index implied correlation](#index-implied-correlation)
- [Realized volatility](#realized-volatility)
- [Delta-hedge P&L simulator](#delta-hedge-pl-simulator)
- [Portfolio risk (VaR / Expected Shortfall)](#portfolio-risk-var--expected-shortfall)
- [Monte Carlo (with variance reduction)](#monte-carlo-with-variance-reduction)
- [Two-asset options](#two-asset-options)
- [Forward-start and cliquet options](#forward-start-and-cliquet-options)
- [Quasi-Monte Carlo](#quasi-monte-carlo)
- [Exotic options (closed form)](#exotic-options-closed-form)
- [Local volatility (Dupire)](#local-volatility-dupire)
- [Term-structure surface (calendar-arbitrage aware)](#term-structure-surface-calendar-arbitrage-aware)
- [Skew/kurtosis-adjusted pricing (Corrado-Su)](#skewkurtosis-adjusted-pricing-corrado-su)
- [CEV (constant elasticity of variance)](#cev-constant-elasticity-of-variance)
- [Variance-Gamma (pure-jump)](#variance-gamma-pure-jump)
- [Bachelier (normal) model](#bachelier-normal-model)
- [Merton jump-diffusion](#merton-jump-diffusion)
- [Heston stochastic volatility](#heston-stochastic-volatility)
- [SABR stochastic-vol smile](#sabr-stochastic-vol-smile)
- [Volatility surface (SVI)](#volatility-surface-svi)
- [Command line](#command-line)
- [Quanto options](#quanto-options)
- [Fixed income](#fixed-income)
- [Callable bonds and OAS](#callable-bonds-and-oas)
- [Carry and roll-down](#carry-and-roll-down)
- [Inflation-linked bonds and derivatives](#inflation-linked-bonds-and-derivatives)
- [Commodities (cost-of-carry and mean reversion)](#commodities-cost-of-carry-and-mean-reversion)
- [Credit (reduced-form)](#credit-reduced-form)
- [Counterparty valuation adjustments (XVA)](#counterparty-valuation-adjustments-xva)
- [Mortgage-backed securities and CMOs](#mortgage-backed-securities-and-cmos)
- [Weather derivatives](#weather-derivatives)
- [Equity compensation and convertibles](#equity-compensation-and-convertibles)
- [Optimal execution](#optimal-execution)
- [Structured notes](#structured-notes)
- [Actuarial (life contingencies and cat bonds)](#actuarial-life-contingencies-and-cat-bonds)
- [Equity swaps and dispersion](#equity-swaps-and-dispersion)
- [Futures/forward convexity](#futuresforward-convexity)
- [Copulas and portfolio credit](#copulas-and-portfolio-credit)
- [Retirement decumulation](#retirement-decumulation)
- [Liability-driven investing](#liability-driven-investing)
- [Bond futures](#bond-futures)
- [Factor models](#factor-models)
- [Performance attribution (Brinson)](#performance-attribution-brinson)
- [Bootstrap and jackknife](#bootstrap-and-jackknife)
- [Markov chains](#markov-chains)
- [Principal component analysis](#principal-component-analysis)
- [Matrix utilities](#matrix-utilities)
- [Numerical utilities](#numerical-utilities)
- [Nelson-Siegel / Svensson curves](#nelson-siegel--svensson-curves)
- [Sample risk measures](#sample-risk-measures)
- [Student-t fat tails](#student-t-fat-tails)
- [Extreme value theory](#extreme-value-theory)
- [Structural credit (Merton)](#structural-credit-merton)
- [Dual-currency deposits](#dual-currency-deposits)
- [FX forwards (covered interest parity)](#fx-forwards-covered-interest-parity)
- [Portfolio optimization](#portfolio-optimization)
- [Portfolio insurance (CPPI)](#portfolio-insurance-cppi)
- [Volatility targeting](#volatility-targeting)
- [Trend and momentum signals](#trend-and-momentum-signals)
- [Pairs trading](#pairs-trading)
- [Performance metrics](#performance-metrics)
- [GARCH volatility](#garch-volatility)
- [Model coverage](#model-coverage)
- [Vectorized fast path (optional NumPy)](#vectorized-fast-path-optional-numpy)
- [Performance](#performance)
- [Testing](#testing)
- [License](#license)
<!-- /TOC -->

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
# The simpler single-boundary Bjerksund-Stensland (1993) is also available:
from quantforge import bjerksund_stensland_1993
print(bjerksund_stensland_1993(S=100, K=100, t=1.0, r=0.05, sigma=0.2,
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

`strategy_report` summarizes any of them — max profit/loss, unbounded-tail
flags, and break-evens:

```python
from quantforge import strategy_report, iron_condor

r = strategy_report(iron_condor(100, 80, 90, 110, 120, 0.5, 0.04, 0.25))
print(r["max_profit"], r["max_loss"], r["break_evens"])
```

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

`book_bump_greeks` gets net delta/gamma/vega/theta by bumping the shared market
and repricing — model-free, so it works for any instrument in the book:

```python
from quantforge import book_bump_greeks

bg = book_bump_greeks(positions)
print(bg.delta, bg.gamma, bg.vega, bg.theta)
```

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

Covered-call / cash-secured-put income yields:

```python
from quantforge import covered_call, cash_secured_put

m = covered_call(S=100, K=105, t=0.25, r=0.04, sigma=0.25)
print(m.annualized_yield, m.breakeven, m.if_assigned_return)
cash_secured_put(S=100, K=95, t=0.25, r=0.04, sigma=0.3)
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

# Back out the correlation implied by a spread-option quote.
from quantforge import implied_spread_correlation
implied_spread_correlation(target_price=10.5, S1=100, S2=95, K=5, t=1.0, r=0.05,
                           sigma1=0.2, sigma2=0.25)

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

Digitals have unbounded pin risk at the strike — `digital_greeks` shows the
delta spiking as expiry nears:

```python
from quantforge import digital_greeks

digital_greeks(S=100, K=100, t=0.02, r=0.05, sigma=0.25)   # large ATM delta
```

so desks super-replicate them with a tight vanilla spread:

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

Merton jump-diffusion also exposes its implied-vol smile:

```python
from quantforge import merton_smile

merton_smile(S=100, strikes=[80, 90, 100, 110, 120], t=1.0, r=0.05,
             sigma=0.2, lam=1.0, mu_j=-0.1, sigma_j=0.15)   # jump smile/skew
```

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

The Vasicek short-rate model prices zero-coupon bonds and options on them:

```python
from quantforge import zero_coupon_bond, bond_option

zero_coupon_bond(r0=0.03, t=5.0, kappa=0.5, theta=0.05, sigma=0.01)
bond_option(r0=0.03, t_option=1.0, t_bond=5.0, strike=0.85,
            kappa=0.5, theta=0.05, sigma=0.01, option_type="call")
```

The Cox-Ingersoll-Ross model keeps rates non-negative (square-root diffusion):

```python
from quantforge import cir_zero_coupon_bond

cir_zero_coupon_bond(r0=0.03, t=5.0, kappa=0.5, theta=0.05, sigma=0.08)
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

`heston_smile` extracts the implied-vol smile the parameters produce:

```python
from quantforge import heston_smile

heston_smile(S=100, strikes=[80, 90, 100, 110, 120], t=1.0, r=0.0,
             v0=0.04, kappa=2.0, theta=0.04, xi=0.5, rho=-0.7)  # downward skew
```

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
from quantforge import svi_is_butterfly_free, svi_butterfly_arbitrage, svi_repair_butterfly
print(svi_is_butterfly_free(params))
print(svi_butterfly_arbitrage(params))   # log-moneyness points that violate
params = svi_repair_butterfly(params)    # shrink the wings until arbitrage-free
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

## Fixed income

Coupon-bond analytics off an explicit cashflow schedule, plus curve-based risk:

```python
from quantforge import (bond_cashflows, yield_to_maturity, macaulay_duration,
                        convexity, key_rate_durations, clean_price)

cf = bond_cashflows(face=100, coupon_rate=0.05, maturity=5, freq=2)
y = yield_to_maturity(cf, price=104.0)
macaulay_duration(cf, y); convexity(cf, y)
key_rate_durations(cf, [0.5, 1, 2, 3, 5], [0.03, 0.032, 0.035, 0.037, 0.04])
```

The three short-rate models expose analytic rate moments too —
`vasicek_expected_rate` / `vasicek_rate_variance` / `vasicek_stationary_distribution`
(Normal), the CIR equivalents (Gamma), and Ho-Lee (drifted Brownian).

## Callable bonds and OAS

Callable/puttable bond pricing on a short-rate binomial tree and the
option-adjusted spread:

```python
from quantforge import (callable_bond_price, straight_bond_tree_price,
                        option_adjusted_spread)

straight = straight_bond_tree_price(100, coupon_rate=0.05, maturity=5, r0=0.04,
                                    sigma=0.15)
callable_ = callable_bond_price(100, 0.05, 5, 0.04, 0.15, call_price=102)  # <= straight
option_adjusted_spread(market_price=99, face=100, coupon_rate=0.05, maturity=5,
                       r0=0.04, sigma=0.15, call_price=110)
```

## Carry and roll-down

Decompose a bond's expected holding return (unchanged curve) into carry and
roll-down:

```python
from quantforge import (carry_return, rolldown_return, total_carry_rolldown,
                        SplineZeroCurve, bond_cashflows)

curve = SplineZeroCurve([0.5, 1, 2, 3, 5, 10], [0.02, 0.025, 0.03, 0.033, 0.037, 0.04])
cf = bond_cashflows(face=100, coupon_rate=0.04, maturity=5, freq=2)
carry_return(coupon_rate=0.04, yield_now=0.037, horizon=1.0, financing_rate=0.02)
rolldown_return(cf, curve, horizon=1.0)     # zero on a flat curve, positive when upward
total_carry_rolldown(cf, curve, 0.04, 1.0, financing_rate=0.02)
```

## Inflation-linked bonds and derivatives

Index-ratio mechanics, TIPS-style linker pricing/risk, the deflation floor,
zero-coupon and year-on-year inflation swaps, an inflation curve from a swap
strip, YoY caplets (lognormal and Bachelier), CPI seasonality and daily
reference-index interpolation, and the Fisher real/nominal curve bridge:

```python
from quantforge import (breakeven_inflation, linker_price, linker_real_duration,
                        zc_inflation_swap_rate, inflation_curve_from_zc_swaps,
                        yoy_caplet_price, real_zero_curve)

breakeven_inflation(nominal_yield=0.05, real_yield=0.0294)      # ~2%
zc_inflation_swap_rate(index_start=100, index_end=133.1, years=3)  # 10% p.a.
real_zero_curve([0.04, 0.045], [0.02, 0.022])                  # strip breakevens
```

The index ratio multiplies the whole linker price, so it cancels in the
fractional real duration/convexity (which equal the standard bond measures on
the real cashflows) while the dollar DV01 scales with it.

## Commodities (cost-of-carry and mean reversion)

Cost-of-carry forwards with storage and convenience yield, implied-yield
inversion, roll yield, the Schwartz (1997) one-factor mean-reverting model
(forward, spot option, Samuelson futures vol, calibration), crack/spark spread
options (Margrabe, Kirk, Bachelier, Monte Carlo), commodity swaps, and average-
price Asians (1/3-variance, exact geometric, Turnbull-Wakeman, MC):

```python
from quantforge import (commodity_forward, implied_convenience_yield, roll_yield,
                        schwartz_forward, kirk_spread_option, turnbull_wakeman_asian)

commodity_forward(spot=100, r=0.05, maturity=1, storage_cost=0.02,
                  convenience_yield=0.03)
schwartz_forward(spot=50, kappa=1.5, alpha_star=4.0, sigma=0.3, maturity=2)
kirk_spread_option(f1=100, f2=90, strike=5, sigma1=0.3, sigma2=0.25, rho=0.4,
                   r=0.05, expiry=1)      # crack spread
```

## Credit (reduced-form)

Piecewise-constant-hazard survival curves, CDS pricing/greeks/bootstrap, and
defaultable bonds:

```python
from quantforge import (SurvivalCurve, cds_par_spread, bootstrap_survival_curve,
                        risky_bond_price)

curve = bootstrap_survival_curve([1, 3, 5], [0.008, 0.010, 0.012], r=0.03)
cds_par_spread(curve, [0.5 * i for i in range(1, 11)], r=0.03, recovery=0.4)
```

## Counterparty valuation adjustments (XVA)

CVA / DVA / bilateral CVA off a survival curve, funding and margin adjustments
(FVA / MVA), analytic swap exposure profiles (EPE / PFE), wrong-way risk, and CSA
collateralization:

```python
from quantforge import (SurvivalCurve, swap_expected_exposure, swap_cva,
                        collateralized_exposure_profile, bcva, mva)

curve = SurvivalCurve([1, 3, 5], [0.02, 0.03, 0.04])
grid = [1, 2, 3, 4, 5]
epe = swap_expected_exposure(notional=1e6, sigma=0.01, maturity=5, grid_times=grid)
swap_cva(curve, 1e6, 0.01, 5, grid, r=0.03, recovery=0.4)

# Collateral (CSA threshold + MTA) caps the exposure and shrinks the CVA.
col = collateralized_exposure_profile(epe, threshold=5000, min_transfer_amount=1000)
```

The CVA default buckets are the marginal survival drops `Q(t_{i-1}) - Q(t_i)`;
`wrong_way_cva` tilts them toward later, higher-exposure dates.

## Mortgage-backed securities and CMOs

Level-payment amortization, PSA / CPR / SMM prepayment, projected pool cashflows,
WAL, price / yield / Z-spread and effective duration, plus sequential and
PAC/support CMO tranching:

```python
from quantforge import (monthly_payment, mbs_cashflows_psa, weighted_average_life,
                        mbs_zspread, sequential_cmo, pac_schedule, pac_support_split)

monthly_payment(balance=300000, annual_rate=0.05, term_months=360)   # 1610.46
flows = mbs_cashflows_psa(300000, 0.05, 360, psa=100)                 # PSA ramp
weighted_average_life(flows, 300000)

# Sequential CMO: earlier tranches retire first (shorter WAL).
tranches = sequential_cmo(flows, [150000, 100000, 50000])

# PAC bond: stable WAL across prepayment speeds inside the 100-300 PSA collar.
sched = pac_schedule(300000, 0.05, 360, psa_low=100, psa_high=300)
pac, support = pac_support_split(flows, sched)
```

## Weather derivatives

Heating/cooling degree days, temperature-index swaps and options (Bachelier), a
degree-day collar, a mean-reverting temperature model, and a Monte Carlo cross-
check:

```python
from quantforge import (heating_degree_days, cooling_degree_days,
                        degree_day_option, expected_temperature)

cooling_degree_days(temps=[70, 72, 68, 75], base=65)
degree_day_option(expected_index=900, strike=880, sigma=45, r=0.03, expiry=0.5,
                  tick_value=20, is_call=True, cap=100)
```

## Equity compensation and convertibles

Warrants with dilution, employee stock options (FASB 123R expected-life model),
discrete cash dividends (escrowed method), and convertible bonds with their
relative-value metrics:

```python
from quantforge import (warrant_price, eso_value, discrete_dividend_price,
                        convertible_bond_value, conversion_premium)

warrant_price(50, 50, 5, 0.05, 0.3, existing_shares=1_000_000, new_shares=100_000)
eso_value(50, 50, contractual_term=10, r=0.05, sigma=0.3, vesting=2,
          exit_rate=0.15, forfeiture_rate=0.03)
discrete_dividend_price(100, 100, 1.0, 0.05, 0.25, dividends=[(0.25, 2), (0.75, 2)])
convertible_bond_value(50, conversion_ratio=20, face=1000, coupon_rate=0.04,
                       maturity=5, r=0.05, sigma=0.3, credit_spread=0.01)

# Full American optionality (early conversion + issuer call) on an equity lattice.
from quantforge import convertible_bond_lattice
convertible_bond_lattice(50, sigma=0.3, face=1000, conversion_ratio=20,
                         coupon_rate=0.04, maturity=5, r=0.05, call_price=1100)
```

## Optimal execution

Almgren-Chriss optimal liquidation (cost/risk trade-off), market-impact models
(Kyle, square-root law), implementation-shortfall decomposition, and TWAP / VWAP /
POV schedules:

```python
from quantforge import (execution_trajectory, implementation_shortfall,
                        square_root_impact, vwap_schedule)

# Risk-averse liquidation front-loads selling; lambda = 0 gives TWAP.
traj = execution_trajectory(total_shares=1e6, n_intervals=10, horizon=1.0,
                            lam=1e-6, sigma=0.3, eta=1e-6)
square_root_impact(order_size=1e5, sigma=0.3, daily_volume=1e7)
vwap_schedule(1e6, volume_profile=[100, 200, 300, 400])
```

## Structured notes

Principal-protected notes (capped and uncapped), reverse convertibles with a fair
coupon solve, and buffered notes, all by component decomposition:

```python
from quantforge import (principal_protected_note, reverse_convertible,
                        reverse_convertible_fair_coupon, buffered_note)

principal_protected_note(100, 100, maturity=3, r=0.04, sigma=0.25, principal=1000,
                         participation=1.0)               # >= discounted principal
reverse_convertible_fair_coupon(100, 100, 3, 0.04, 0.25, principal=1000)  # prices at par
buffered_note(100, 100, buffer=0.10, maturity=3, r=0.04, sigma=0.25, principal=1000)
```

## Actuarial (life contingencies and cat bonds)

Life-table survival, annuities and insurance EPVs, the equivalence-principle net
premium, the Gompertz-Makeham mortality law, and catastrophe-bond analytics:

```python
from quantforge import (gompertz_makeham_survival_curve, life_annuity_due,
                        whole_life_insurance, net_level_premium,
                        curtate_life_expectancy, cat_bond_price)

px = gompertz_makeham_survival_curve(age=40, n_years=60, a=0.0005, b=7.58e-5, c=1.09)
life_annuity_due(px, i=0.05)
net_level_premium(px, i=0.05)          # A_x / a-due_x (equivalence principle)
curtate_life_expectancy(px)            # e_40 ~ 35 years

cat_bond_price(principal=1000, coupon_rate=0.06, expected_loss_rate=0.02,
               r=0.03, maturity=1.0)
```

## Equity swaps and dispersion

Total-return swaps, dividend swaps, variance/volatility swaps, and dispersion-
trade P&L:

```python
from quantforge import (total_return_swap_value, dividend_swap_fair_strike,
                        variance_swap_payoff, dispersion_trade_pnl)

total_return_swap_value(1e6, start_price=100, end_price=110, dividends=3,
                        funding_rate=0.03, spread=0.005, year_fraction=1.0)
variance_swap_payoff(realized_vol=0.25, strike_vol=0.20, variance_notional=1e6)
```

## Futures/forward convexity

Ho-Lee and Hull-White convexity adjustments converting interest-rate futures to
forwards, plus a futures-strip-to-forward-curve bootstrap:

```python
from quantforge import forward_from_futures, forward_curve_from_futures_strip

forward_from_futures(futures_rate=0.05, sigma=0.01, t1=2.0, t2=2.25)   # < futures
forward_curve_from_futures_strip([(0.25, 0.5, 0.05), (0.5, 0.75, 0.052)], sigma=0.012)
```

## Copulas and portfolio credit

Gaussian, Clayton, Gumbel, and Frank copulas with tail dependence and Kendall's-
tau calibration, plus Gaussian-copula credit applications (joint/first-to-default)
and the Vasicek large-pool CDO tranche loss:

```python
from quantforge import (gaussian_copula, clayton_lower_tail_dependence,
                        first_to_default_probability, vasicek_loss_quantile,
                        cdo_tranche_expected_loss)

gaussian_copula(0.4, 0.7, rho=0.5)                       # joint CDF
clayton_lower_tail_dependence(theta=2)                   # joint-crash dependence
first_to_default_probability(pd1=0.05, pd2=0.08, rho=0.5)
vasicek_loss_quantile(q=0.999, pd=0.02, rho=0.15)        # Basel IRB capital
cdo_tranche_expected_loss(attachment=0.03, detachment=0.07, pd=0.05, rho=0.2)
```

## Retirement decumulation

Sustainable withdrawals, portfolio depletion, target-date glide paths, and a
Monte Carlo ruin probability:

```python
from quantforge import (sustainable_withdrawal, portfolio_depletion_years,
                        glide_path_equity_weight, ruin_probability_mc)

sustainable_withdrawal(balance=1e6, real_return=0.03, years=30)   # depletes in 30y
portfolio_depletion_years(1e6, annual_withdrawal=40000, real_return=0.02)  # ~35y
glide_path_equity_weight(years_to_target=10, glide_years=30, start_equity=0.9,
                         end_equity=0.3)
ruin_probability_mc(1e6, annual_withdrawal=50000, mean_return=0.04, vol=0.12,
                    years=30)
```

## Liability-driven investing

Funding ratios, liability duration/convexity, duration hedging, Redington
immunization, and surplus-at-risk:

```python
from quantforge import (liability_pv, funding_ratio, liability_duration,
                        required_hedge_duration, surplus_at_risk)

liabilities = [(5, 100), (10, 150), (20, 250), (30, 300)]
L = liability_pv(liabilities, discount_rate=0.03)
funding_ratio(assets=1.05 * L, liabilities=L)
D = liability_duration(liabilities, 0.03)
required_hedge_duration(asset_value=1.05 * L, liability_duration_=D, liability_value=L)
surplus_at_risk(1.05 * L, L, surplus_volatility=0.08, confidence=0.95)
```

## Bond futures

Conversion factors, delivery invoice/basis, cheapest-to-deliver selection,
implied repo, and futures DV01 hedging:

```python
from quantforge import (conversion_factor, net_basis, cheapest_to_deliver,
                        futures_dv01, futures_hedge_ratio, bond_future_dv01)

conversion_factor(coupon_rate=0.04, years_to_maturity=10, notional_coupon=0.06)
basket = [{"price": 115, "cf": 0.95, "carry": 0.5},
          {"price": 118, "cf": 0.98, "carry": 0.4}]
cheapest_to_deliver(basket, futures_price=120)

# Hedge a cash bond with futures via DV01 matching.
fut_dv01 = futures_dv01(ctd_dv01=bond_future_dv01(115, 9.0), ctd_conversion_factor=0.95)
futures_hedge_ratio(bond_dv01_=bond_future_dv01(120, 8.5), futures_dv01_=fut_dv01)
```

## Factor models

Multi-factor OLS return regression (alpha, betas, R-squared), factor attribution,
and rolling beta:

```python
from quantforge import factor_regression, factor_attribution, rolling_factor_beta

fit = factor_regression(asset_returns, [market_returns, size_returns, value_returns])
fit["alpha"], fit["betas"], fit["r_squared"]
factor_attribution(total_return=0.05, alpha=fit["alpha"], betas=fit["betas"],
                   factor_realized_returns=[0.02, 0.01, -0.005])
```

## Performance attribution (Brinson)

Brinson-Hood-Beebower allocation/selection/interaction effects with Cariño
multi-period geometric linking:

```python
from quantforge import brinson_attribution, carino_linked_effects

res = brinson_attribution(portfolio_weights=[0.5, 0.3, 0.2],
                          benchmark_weights=[0.4, 0.4, 0.2],
                          portfolio_returns=[0.10, 0.05, 0.08],
                          benchmark_returns=[0.08, 0.06, 0.07])
res["allocation_total"], res["selection_total"], res["active_return"]
```

## Bootstrap and jackknife

Nonparametric confidence intervals: IID, stationary (block) bootstrap for
serially-correlated data, BCa, and the jackknife:

```python
from quantforge import (bootstrap_ci, stationary_bootstrap_ci, bca_bootstrap_ci,
                        jackknife_estimate)

returns = [0.01, -0.02, 0.03, 0.00, 0.015, -0.01, 0.025]
bootstrap_ci(returns, confidence=0.95)          # (lower, point, upper) for the mean
bca_bootstrap_ci(returns)                        # bias-corrected accelerated
stationary_bootstrap_ci(returns, mean_block=3)   # block bootstrap for autocorrelation
jackknife_estimate(returns)                      # (estimate, standard_error)
```

## Markov chains

Finite-state chains: n-step transitions, stationary distribution, hitting times,
and absorbing-chain analytics (credit-rating migration):

```python
from quantforge import (stationary_distribution, n_step_transition,
                        absorption_probabilities, expected_steps_to_absorption)

P = [[0.9, 0.1], [0.2, 0.8]]
stationary_distribution(P)          # long-run state distribution
n_step_transition(P, 5)             # P^5

# Absorbing chain (states 0,1 transient, 2 absorbing).
A = [[0.5, 0.3, 0.2], [0.1, 0.6, 0.3], [0.0, 0.0, 1.0]]
expected_steps_to_absorption(A, transient_states=[0, 1])
absorption_probabilities(A, transient_states=[0, 1], absorbing_states=[2])
```

## Principal component analysis

Covariance-matrix PCA via Jacobi eigendecomposition -- the yield-curve
level/slope/curvature factors -- plus low-rank reconstruction and component
scenarios:

```python
from quantforge import pca, pca_scenario, reconstruct_covariance

cov = [[0.04, 0.02, 0.01], [0.02, 0.03, 0.015], [0.01, 0.015, 0.025]]
res = pca(cov)
res["variances"], res["loadings"], res["cumulative_explained"]

# A 2-sigma move along the first component (a level shift for a rate curve).
pca_scenario(component_index=0, n_sigma=2.0, variances=res["variances"],
             loadings=res["loadings"])
```

## Matrix utilities

Cholesky factorization, positive-definiteness check, correlated-normal draws,
Higham nearest-correlation repair, and an n-asset basket Monte Carlo:

```python
from quantforge import (cholesky, nearest_correlation, correlated_normals,
                        basket_option_mc)

corr = [[1, 0.3, 0.2], [0.3, 1, 0.4], [0.2, 0.4, 1]]
cholesky(corr)
correlated_normals([0.5, -1.0, 2.0], corr)     # IID normals -> correlated
nearest_correlation([[1, 0.9, -0.9], [0.9, 1, 0.9], [-0.9, 0.9, 1]])  # repair indefinite

basket_option_mc([100, 100, 100], [1/3, 1/3, 1/3], strike=100, t=1, r=0.05,
                 sigmas=[0.2, 0.25, 0.3], correlation=corr)
```

## Numerical utilities

Cubic interpolation (natural spline and monotone Hermite), a spline-interpolated
zero curve, and general-purpose root finders (bisection, Brent, Newton):

```python
from quantforge import (natural_cubic_spline, monotone_cubic, SplineZeroCurve,
                        brent, newton)

f = monotone_cubic([0, 1, 2, 3], [0, 0, 0, 1])   # no overshoot
curve = SplineZeroCurve([0.5, 1, 2, 5, 10], [0.02, 0.025, 0.03, 0.035, 0.04])
brent(lambda x: x * x - 2, 0, 2)                  # sqrt(2)
```

## Nelson-Siegel / Svensson curves

Parametric yield curves and least-squares calibration to observed zeros:

```python
from quantforge import (nelson_siegel_zero, svensson_zero, fit_nelson_siegel)

nelson_siegel_zero(t=5, beta0=0.04, beta1=-0.02, beta2=0.01, tau=2.0)
mats = [0.5, 1, 2, 3, 5, 7, 10, 20, 30]
zeros = [0.02, 0.025, 0.028, 0.03, 0.033, 0.035, 0.037, 0.039, 0.04]
beta0, beta1, beta2, tau = fit_nelson_siegel(mats, zeros)   # calibrate
```

## Sample risk measures

VaR, expected shortfall, spectral and entropic risk from a P&L sample, plus
coherence checks and Euler component-ES allocation:

```python
from quantforge import (value_at_risk, sample_expected_shortfall,
                        spectral_risk_exponential, component_expected_shortfall)

pnl = [0.01, -0.02, 0.03, -0.05, 0.02, -0.10, 0.015, -0.03, -0.08, 0.01]
value_at_risk(pnl, confidence=0.95)
sample_expected_shortfall(pnl, confidence=0.95)     # coherent CVaR
component_expected_shortfall([book_a_pnl, book_b_pnl])   # contributions sum to total ES
```

## Student-t fat tails

The Student-t distribution (pdf/cdf/quantile), fat-tailed parametric VaR / ES,
and degrees-of-freedom fitting from a return sample:

```python
from quantforge import (t_ppf, student_t_var, student_t_expected_shortfall,
                        fit_student_t)

student_t_var(mean=0, scale=0.02, df=4, confidence=0.99)   # fatter than normal
student_t_expected_shortfall(0, 0.02, df=4, confidence=0.99)
mean, scale, df = fit_student_t(returns)                    # moment-match the tails
```

## Extreme value theory

Tail-index estimation and peaks-over-threshold Generalized Pareto VaR / ES:

```python
from quantforge import hill_estimator, gpd_var, gpd_expected_shortfall

hill_estimator(losses, k=1000)                        # tail index of the top k
threshold = sorted(losses)[int(0.95 * len(losses))]
gpd_var(losses, threshold, confidence=0.99)           # peaks-over-threshold VaR
gpd_expected_shortfall(losses, threshold, 0.99)
```

## Structural credit (Merton)

Firm equity as a call on assets, distance-to-default, default probability, credit
spread, and the KMV solve for unobservable asset value and volatility:

```python
from quantforge import (distance_to_default, risk_neutral_default_probability,
                        credit_spread, solve_asset_value_and_vol)

distance_to_default(asset_value=120, debt_face=100, r=0.05, sigma=0.25, t=1)
risk_neutral_default_probability(120, 100, 0.05, 0.25, 1)
credit_spread(120, 100, 0.05, 0.25, 1)

# Back out the unobservable assets from observed equity value and equity vol.
solve_asset_value_and_vol(equity_value_obs=27.4, equity_vol_obs=0.94,
                          debt_face=100, r=0.05, t=1)
```

## Dual-currency deposits

Yield-enhanced FX-linked deposits by component decomposition:

```python
from quantforge import (dcd_option_premium_rate, dcd_enhanced_yield,
                        dcd_breakeven_spot)

prem = dcd_option_premium_rate(spot=1.10, strike=1.12, tenor=0.25, r_domestic=0.05,
                               r_foreign=0.03, sigma=0.10)
dcd_enhanced_yield(base_deposit_rate=0.05, option_premium_rate=prem, tenor=0.25)
```

## FX forwards (covered interest parity)

```python
from quantforge import fx_forward, forward_points, implied_base_rate

fx_forward(spot=1.10, r_price=0.05, r_base=0.03, t=1.0)   # EURUSD-style
```

## Portfolio optimization

Mean-variance optimizers and risk decomposition from a covariance matrix, plus
Black-Litterman:

```python
from quantforge import (min_variance_weights, max_sharpe_weights,
                        risk_parity_weights, efficient_frontier,
                        black_litterman_weights, portfolio_var, var_budget)

cov = [[0.04, 0.01, 0.0], [0.01, 0.09, 0.02], [0.0, 0.02, 0.16]]
min_variance_weights(cov)
risk_parity_weights(cov)          # equal risk contributions
var_budget([0.4, 0.4, 0.2], cov)  # % risk per position

# Clustering-based and inverse-vol allocation (no matrix inversion).
from quantforge import inverse_volatility_weights, hierarchical_risk_parity
inverse_volatility_weights(cov)
hierarchical_risk_parity(cov)     # Lopez de Prado HRP
```

## Portfolio insurance (CPPI)

Constant Proportion Portfolio Insurance: dynamic risky/safe allocation with a
guaranteed floor.

```python
from quantforge import risky_exposure, cppi_path

risky_exposure(wealth=100, floor_pv=90, multiplier=3)   # = 3 * cushion, capped
cppi_path(initial_wealth=100, floor=90, multiplier=3,
          risky_returns=[-0.30, -0.30, -0.30, -0.30], r=0.02, dt=0.25)  # protected
```

## Volatility targeting

Scale exposure inversely with realized volatility to run at a constant risk level:

```python
from quantforge import target_leverage, vol_targeted_returns, realized_annualized_vol

target_leverage(target_vol=0.10, realized_vol=0.20)          # 0.5x (de-risk)
overlaid = vol_targeted_returns(returns, target_vol=0.15, lookback=60)
realized_annualized_vol(overlaid)                            # ~ 0.15
```

## Trend and momentum signals

Moving averages, MACD, RSI, rolling z-score, and time-series momentum:

```python
from quantforge import sma, ema, macd, rsi, rolling_zscore, time_series_momentum

ema(prices, span=12)
macd(prices, fast=12, slow=26, signal=9)      # (macd_line, signal_line, histogram)
rsi(prices, window=14)                         # 0-100 oscillator
time_series_momentum(prices, lookback=60)      # +1 / -1 / 0 trend sign

from quantforge import bollinger_bands, average_true_range, donchian_channel
bollinger_bands(prices, window=20, num_std=2.0)          # (lower, middle, upper)
donchian_channel(highs, lows, window=20)                 # breakout channel
```

## Pairs trading

Hedge ratio, spread, Ornstein-Uhlenbeck half-life, and spread z-score:

```python
from quantforge import pairs_hedge_ratio, spread_series, ou_half_life, spread_zscore

beta = pairs_hedge_ratio(y, x)          # OLS slope of y on x
spread = spread_series(y, x, beta)
ou_half_life(spread)                     # mean-reversion half-life (years/periods)
spread_zscore(spread, window=60)         # entry/exit signal
```

## Performance metrics

Track-record statistics from a return series:

```python
from quantforge import (sharpe_ratio, sortino_ratio, max_drawdown,
                        calmar_ratio, information_ratio, up_capture,
                        ulcer_index, ulcer_performance_index,
                        probabilistic_sharpe_ratio, deflated_sharpe_ratio)

sharpe_ratio(returns); max_drawdown(returns); calmar_ratio(returns)

# Drawdown pain and statistical significance of the Sharpe.
ulcer_index(returns); ulcer_performance_index(returns)      # Martin ratio
probabilistic_sharpe_ratio(returns, benchmark_sr=0.0)       # P(true SR > 0)
deflated_sharpe_ratio(returns, n_trials=100)                # selection-bias corrected
```

## GARCH volatility

Fit GARCH(1,1), forecast the term volatility, and price consistently with the
vol mean-reversion:

```python
from quantforge import fit_garch, garch_term_variance, garch_option_price

params = fit_garch(returns)
garch_option_price(params, last_return, last_variance,
                   S=100, K=100, r=0.05, option_type="call", horizon=21)
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

Implied-vol round-trips to a max price error of `1e-8`, and the
Newton + Corrado-Miller solver needs ~4-5 iterations vs ~29 for pure bisection.
Reproduce with:

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
