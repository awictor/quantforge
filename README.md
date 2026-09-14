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
- [Shout and ladder options](#shout-and-ladder-options)
- [Installment options](#installment-options)
- [Double-barrier knock-out](#double-barrier-knock-out)
- [Range-accrual note](#range-accrual-note)
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
- [Equity valuation](#equity-valuation)
- [Capital budgeting](#capital-budgeting)
- [Money-market yields](#money-market-yields)
- [Black-Karasinski short rate](#black-karasinski-short-rate)
- [Hull-White (fitted to a curve)](#hull-white-fitted-to-a-curve)
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
- [Hawkes self-exciting process](#hawkes-self-exciting-process)
- [Structured notes](#structured-notes)
- [Survival analysis (Kaplan-Meier / Nelson-Aalen)](#survival-analysis-kaplan-meier--nelson-aalen)
- [Actuarial (life contingencies and cat bonds)](#actuarial-life-contingencies-and-cat-bonds)
- [Equity swaps and dispersion](#equity-swaps-and-dispersion)
- [Futures/forward convexity](#futuresforward-convexity)
- [Copulas and portfolio credit](#copulas-and-portfolio-credit)
- [kth-to-default basket](#kth-to-default-basket)
- [Retirement decumulation](#retirement-decumulation)
- [Liability-driven investing](#liability-driven-investing)
- [Bond futures](#bond-futures)
- [OLS regression](#ols-regression)
- [Logistic regression](#logistic-regression)
- [Classification metrics](#classification-metrics)
- [Forecast calibration (Brier decomposition)](#forecast-calibration-brier-decomposition)
- [Cross-validation](#cross-validation)
- [Feature scaling](#feature-scaling)
- [k-nearest neighbors](#k-nearest-neighbors)
- [Gaussian naive Bayes](#gaussian-naive-bayes)
- [Decision stump](#decision-stump)
- [Factor models](#factor-models)
- [Performance attribution (Brinson)](#performance-attribution-brinson)
- [Bootstrap and jackknife](#bootstrap-and-jackknife)
- [Hodrick-Prescott filter](#hodrick-prescott-filter)
- [Kalman filter (local level)](#kalman-filter-local-level)
- [Newey-West HAC variance](#newey-west-hac-variance)
- [Theil-Sen robust regression](#theil-sen-robust-regression)
- [Isotonic regression (monotone fit)](#isotonic-regression-monotone-fit)
- [Robust scale and location](#robust-scale-and-location)
- [Hurst exponent (long memory)](#hurst-exponent-long-memory)
- [Entropy (time-series regularity)](#entropy-time-series-regularity)
- [Variance-ratio test](#variance-ratio-test)
- [Serial correlation (Ljung-Box / Durbin-Watson)](#serial-correlation-ljung-box--durbin-watson)
- [Benford's law (first-digit anomaly detection)](#benfords-law-first-digit-anomaly-detection)
- [Goodness of fit (Jarque-Bera / KS)](#goodness-of-fit-jarque-bera--ks)
- [Autocorrelation (ACF / PACF)](#autocorrelation-acf--pacf)
- [Exponential smoothing (Holt / Holt-Winters)](#exponential-smoothing-holt--holt-winters)
- [Forecast accuracy](#forecast-accuracy)
- [Rank dependence (Kendall / Spearman)](#rank-dependence-kendall--spearman)
- [Gaussian-copula sampling](#gaussian-copula-sampling)
- [Directional statistics](#directional-statistics)
- [Geodesy](#geodesy)
- [Quaternions](#quaternions)
- [Color spaces](#color-spaces)
- [Spectral analysis](#spectral-analysis)
- [Wavelet transform (Haar multiresolution)](#wavelet-transform-haar-multiresolution)
- [Structural breaks (CUSUM / Chow)](#structural-breaks-cusum--chow)
- [Cointegration (ADF / Engle-Granger)](#cointegration-adf--engle-granger)
- [Ornstein-Uhlenbeck calibration](#ornstein-uhlenbeck-calibration)
- [Ornstein-Uhlenbeck calibration](#ornstein-uhlenbeck-calibration)
- [Markov chains](#markov-chains)
- [K-means clustering](#k-means-clustering)
- [Hierarchical clustering](#hierarchical-clustering)
- [Gaussian mixture model](#gaussian-mixture-model)
- [Hidden Markov model](#hidden-markov-model)
- [Principal component analysis](#principal-component-analysis)
- [Market stress (turbulence / absorption ratio)](#market-stress-turbulence--absorption-ratio)
- [Matrix utilities](#matrix-utilities)
- [Numerical utilities](#numerical-utilities)
- [Range-query structures](#range-query-structures)
- [Interval arithmetic](#interval-arithmetic)
- [Number theory](#number-theory)
- [Computational geometry](#computational-geometry)
- [Graph algorithms](#graph-algorithms)
- [String algorithms](#string-algorithms)
- [Data compression](#data-compression)
- [Probability distributions](#probability-distributions)
- [Nelson-Siegel / Svensson curves](#nelson-siegel--svensson-curves)
- [Sample risk measures](#sample-risk-measures)
- [Student-t fat tails](#student-t-fat-tails)
- [Extreme value theory](#extreme-value-theory)
- [Structural credit (Merton)](#structural-credit-merton)
- [Dual-currency deposits](#dual-currency-deposits)
- [FX forwards (covered interest parity)](#fx-forwards-covered-interest-parity)
- [Entropy pooling (views on scenarios)](#entropy-pooling-views-on-scenarios)
- [Covariance shrinkage (Ledoit-Wolf)](#covariance-shrinkage-ledoit-wolf)
- [Portfolio optimization](#portfolio-optimization)
- [Rebalancing](#rebalancing)
- [Portfolio insurance (CPPI)](#portfolio-insurance-cppi)
- [Leveraged ETFs](#leveraged-etfs)
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

Across correlated assets the growth-optimal leverage vector is
`f* = Sigma^{-1} mu`. Positive correlation pulls the total leverage below the
naive per-asset sum; feed it the shrunk covariance from `ledoit_wolf_shrinkage`.

```python
from quantforge import kelly_fractions_multivariate

mu = [0.08, 0.05]                       # excess returns
cov = [[0.04, 0.012], [0.012, 0.02]]    # covariance
kelly_fractions_multivariate(mu, cov)                 # -> [1.5244, 1.5854]
kelly_fractions_multivariate(mu, cov, fraction=0.5)   # -> [0.7622, 0.7927]  (half-Kelly)
```

Sizing has a flip side: the chance of going broke. `gamblers_ruin_probability`
gives the classic unit-stake ruin probability, and `ruin_probability_gbm` the
chance a drifting account ever falls by a given fraction — the first-passage law
`(1 - loss)^{2 mu / sigma^2}`:

```python
from quantforge import (gamblers_ruin_probability, risk_of_ruin_units,
                        ruin_probability_gbm)

gamblers_ruin_probability(start=5, target=20, win_prob=0.55)   # ruin before target
risk_of_ruin_units(win_prob=0.55, target_units=20)             # bankroll-in-units
ruin_probability_gbm(loss=0.3, mu=0.15, sigma=0.30)            # 0.305, P(ever -30%)
```

A fair game gives the linear `1 - start/target`; a stronger edge (higher
drift-to-variance) makes a given loss less likely, and a driftless account hits any
loss almost surely.

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

# HAR-RV (Corsi) long-memory realized-variance forecast.
from quantforge import fit_har_rv, har_rv_forecast
coeffs = fit_har_rv(realized_variance_series)          # day/week/month regression
har_rv_forecast(coeffs, realized_variance_series)      # one-step-ahead RV
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

For **intraday** data, decompose the day's quadratic variation into a continuous
part and jumps. Realized variance captures everything; bipower variation is
jump-robust; their difference is the jump contribution.

```python
from quantforge import (realized_variance_from_returns, bipower_variation,
                        jump_variation)

r = [0.001, -0.002, 0.0015, -0.001, 0.05, 0.002,   # one big jump mid-series
     -0.0015, 0.001, -0.002, 0.0018, -0.0012, 0.0009]

realized_variance_from_returns(r)   # -> 0.002525  (total quadratic variation)
bipower_variation(r)                # -> 0.000267  (continuous part, jump-robust)
jump_variation(r)                   # -> 0.002258  (max(RV - BV, 0), the jump)
```

On a jump-free path RV and BV coincide and the jump variation is ~0; a jump lifts
RV by roughly its square while BV barely moves.

`min_realized_variance` and `med_realized_variance` (Andersen-Dobrev-Schaumburg)
are nearest-neighbour jump-robust alternatives to bipower — MinRV keeps the smaller
of each adjacent pair, MedRV the median of each triple — and `realized_quarticity`
estimates the integrated quarticity that sets RV's standard error:

```python
from quantforge import (min_realized_variance, med_realized_variance,
                        realized_quarticity)

min_realized_variance(r)    # jump-robust, discards the lone jump
med_realized_variance(r)    # median of three; also robust to two nearby jumps
realized_quarticity(r)      # (n/3) sum r^4, for jump-test standard errors
```

Like bipower these ignore the jump and track the continuous integrated variance;
MedRV additionally tolerates two close jumps and occasional zero returns.

To decide *whether* a jump occurred, `bns_jump_test` gives the Barndorff-Nielsen-
Shephard ratio statistic — `(RV - BV) / RV` standardized by the `tripower_quarticity`
— which is asymptotically standard normal under the no-jump null:

```python
from quantforge import bns_jump_test, tripower_quarticity

z, p = bns_jump_test(r)      # large z / small p rejects "no jump"
tripower_quarticity(r)       # jump-robust integrated-quarticity scale
```

Under the null the statistic holds its ~5% size; a genuine jump drives `z` well
above the normal critical value (power ~1 in simulation).

At the finest sampling frequencies, microstructure noise biases the naive realized
variance *upward* by `2 * n * Var(noise)` — and faster sampling makes it worse.
`two_scale_realized_variance` (Zhang-Mykland-Aït-Sahalia) removes that bias by
combining a subsampled slow scale with the fast scale:

```python
from quantforge import two_scale_realized_variance, realized_variance_naive

# noisy tick log-prices; true integrated variance is 0.002
realized_variance_naive(prices)        # -> 0.00309  (inflated by noise)
two_scale_realized_variance(prices)    # -> 0.00224  (bias-corrected)
```

`realized_kernel` (Barndorff-Nielsen-Hansen-Lunde-Shephard) is an alternative
noise-robust estimator that adds flat-top Parzen-weighted return autocovariances:

```python
from quantforge import realized_kernel

realized_kernel(prices)                # noise-robust integrated variance
realized_kernel(prices, bandwidth=20)  # or fix the number of lags
```

The autocovariance terms cancel the noise inflation in the sum of squared returns,
so the kernel tracks the true integrated variance where the naive estimator runs
~10x high; the bandwidth defaults to the `n^{3/5}` rule of thumb.

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

Once a model is live, backtest it against realized losses. The coverage tests check
the exception rate and clustering; the Acerbi-Szekely statistic checks the ES:

```python
from quantforge import (kupiec_pof, christoffersen_cc, acerbi_szekely_es)

kupiec_pof(losses, var_forecasts, alpha=0.01)        # (LR, p): unconditional coverage
christoffersen_cc(losses, var_forecasts, alpha=0.01) # (LR, p): coverage + independence
acerbi_szekely_es(losses, var_forecasts, es_forecasts, alpha=0.01)   # ~0 if ES calibrated
```

A well-specified 99% VaR passes Kupiec and Christoffersen (large p-values); too-low
a VaR is rejected. The Acerbi-Szekely statistic sits near zero when ES is right,
turns positive when the model understates the tail, negative when it overstates.
`christoffersen_independence` isolates the exception-clustering piece on its own.

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

For a general diffusion `dX = a(X,t) dt + b(X,t) dW`, the SDE integrators simulate arbitrary
drift/diffusion paths, not just GBM:

```python
from quantforge import euler_maruyama, milstein, gbm_paths

# any 1-D SDE from drift and diffusion callables
euler_maruyama(lambda x, t: 0.1 * x, lambda x, t: 0.2 * x, x0=100, t=1.0, n_steps=50)

# GBM convenience wrapper — terminal mean tracks x0 * exp(mu*t)
paths = gbm_paths(mu=0.1, sigma=0.2, x0=100, t=1.0, n_steps=50, n_paths=3000)
```

`euler_maruyama` is the basic scheme (strong order 0.5); `milstein` adds the `0.5 b b'
(dW^2 - dt)` correction for strong order 1.0, converging faster for state-dependent
diffusions — on a GBM path with shared Brownian increments its error against the exact
solution is ~30x smaller than Euler's. `gbm_paths` wraps them for geometric Brownian
motion (its terminal sample mean approaches the analytic `x0 exp(mu t)`), and all use a
deterministic seeded normal stream so paths reproduce exactly.

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

For a basket of **more than two** assets, `levy_basket_option` applies the same
lognormal moment match across an arbitrary correlation matrix:

```python
from quantforge import levy_basket_option

spots, w, sig = [100, 95, 105], [1/3, 1/3, 1/3], [0.20, 0.25, 0.22]
corr = [[1, 0.6, 0.5], [0.6, 1, 0.55], [0.5, 0.55, 1]]
levy_basket_option(spots, w, strike=100, t=1.0, r=0.04, sigmas=sig, corr=corr)
# -> 9.4216   (matches a Monte Carlo basket to ~0.05)
```

It reduces to Black-Scholes for a single asset and satisfies put-call parity
exactly.

Best-of / worst-of options extend to any number of assets by Monte Carlo:

```python
from quantforge import rainbow_option_mc

spots, sig = [100, 95, 105], [0.2, 0.25, 0.3]
corr = [[1, 0.5, 0.3], [0.5, 1, 0.4], [0.3, 0.4, 1]]
rainbow_option_mc(spots, strike=100, t=1.0, r=0.05, sigmas=sig, corr=corr, best=True)   # best-of call
rainbow_option_mc(spots, strike=100, t=1.0, r=0.05, sigmas=sig, corr=corr, best=False)  # worst-of call
```

The two-asset case matches the closed-form rainbow prices, and always
`worst-of <= single-asset <= best-of`.

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

Where Halton and Sobol generate a *sequence*, Latin hypercube sampling designs a *fixed*
set of `n` points that stratifies every axis — each coordinate is split into `n` bins
with exactly one sample per bin, so no dimension clusters or leaves gaps:

```python
from quantforge import latin_hypercube, maximin_lhs, l2_star_discrepancy

design = latin_hypercube(30, dim=2)          # 30 stratified points in the unit square
l2_star_discrepancy(design)                  # uniformity score — smaller is better
maximin_lhs(15, dim=2, tries=30)             # most space-filling of 30 candidate designs
```

Averaged over many seeds a Latin hypercube's L2 star discrepancy (`~0.032` for 30 points
in 2-D) is about half that of plain random sampling (`~0.068`), so it fills space more
evenly with the same point count — the standard design for Monte-Carlo experiments and
sensitivity analysis. `maximin_lhs` picks the design whose closest pair is farthest apart
(raising the minimum inter-point distance from `~0.10` to `~0.15` here), and
`l2_star_discrepancy` (Warnock's formula) scores any point set, so you can compare an LHS,
Sobol, or Halton design head to head.

For repeatedly drawing from a fixed *discrete* distribution — scenario selection, weighted
resampling — `AliasSampler` makes each draw O(1) after an O(n) setup:

```python
from quantforge import AliasSampler

sampler = AliasSampler([0.1, 0.2, 0.3, 0.4], seed=42)
sampler.sample()             # one index in [0, 3], drawn with those probabilities
sampler.sample_many(100000)  # empirical frequencies converge to [0.1, 0.2, 0.3, 0.4]
```

Walker's alias method builds two tables so each `sample()` is a single table lookup rather
than an O(n) cumulative search — the right structure when you draw many times from the
same weights. Weights need not sum to 1, and a fixed `seed` makes the stream reproducible.

For *continuous* distributions, the variate samplers draw from the common families
directly:

```python
from quantforge import sample_normal, sample_exponential, sample_gamma, sample_poisson

sample_normal(1000, mu=5, sigma=2)      # Box-Muller normals
sample_exponential(1000, rate=0.5)      # mean 1/rate = 2
sample_gamma(1000, shape=2.0, scale=2.0)# mean shape*scale = 4
sample_poisson(1000, lam=4.0)           # integer counts, mean = variance = lam
```

`sample_normal` uses Box-Muller, `sample_exponential` the inverse CDF, `sample_gamma`
Marsaglia-Tsang (with the small-shape boost for `shape < 1`), and `sample_poisson`
Knuth's algorithm. Each runs on a deterministic seeded stream (reproducible per `seed`),
and the sample moments match the distribution's — the normal draws recover `(5, 2)`, the
gamma `(mean 4, var 8)`, and the Poisson `mean = var = 4` over enough samples.

The samplers above run on a simple LCG; when the low-bit quality matters, `PCG32` and
`Xorshift128Plus` are statistically strong drop-ins:

```python
from quantforge import PCG32, Xorshift128Plus

rng = PCG32(seed=42, seq=54)
rng.next_uint32()            # 0xa15c02b7 — matches O'Neill's reference stream
rng.random()                 # uniform [0, 1)
rng.randint(0, 5)            # bias-free integer in [0, 5]
```

`PCG32` is O'Neill's permuted congruential generator (its output reproduces the canonical
reference vector bit for bit), with a rejection-sampled `randint` that has no modulo bias.
`Xorshift128Plus` is Vigna's fast 64-bit generator seeded through a splitmix64 warm-up.
Both are uniform (mean ~0.5, variance ~1/12) and reproducible per seed, and pass the
empirical tests a bare LCG fails.

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

# Sharper discrete arithmetic Asian by Curran's geometric conditioning
# (forward measure; n_avg monitoring dates). Matches Monte Carlo to ~0.001.
from quantforge import curran_asian
curran_asian(forward=100, strike=100, sigma=0.25, r=0.05, expiry=1.0, n_avg=12)

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

## Shout and ladder options

Two path-dependent lock-in calls priced on a Cox-Ross-Rubinstein tree. A **shout**
call lets the holder shout once before expiry to lock in the then-current intrinsic
`S* - K` as a floor while keeping the upside, so the terminal payoff is
`max(S_T - K, S* - K)`. A **ladder** call floors its payoff at the highest preset
rung the underlying touches: `max(S_T - K, max_touched L_i - K, 0)`.

```python
from quantforge import shout_call, ladder_call

# Shout: lock in intrinsic once, keep upside.
shout_call(S=100, K=100, t=1.0, r=0.05, sigma=0.25)          # -> 15.4261

# Ladder: floor payoff at highest rung touched (rungs must exceed the strike).
ladder_call(S=100, K=100, rungs=[110, 120, 130], t=1.0, r=0.05, sigma=0.25)  # -> 17.3279
```

Both are worth at least the vanilla call. With no rungs the ladder reduces to the
plain tree-priced European call; rungs at or below the strike are ignored. More or
higher rungs raise the ladder value up toward the shout-like limit.

## Installment options

An installment call is paid for in a stream of premiums rather than a single
upfront cost. At each installment date the holder may **lapse** — stop paying and
forfeit the option for zero — so it is kept alive only while its continuation value
exceeds the next installment. That abandon-option makes it a compound option,
priced by backward induction on a Cox-Ross-Rubinstein tree. `installment_call`
returns the fair *upfront* value given the agreed schedule.

```python
from quantforge import installment_call

# Pay 2.0 at each of three dates to keep a 1y ATM call alive.
installment_call(S=100, K=100, t=1.0, r=0.05, sigma=0.25,
                 installment=2.0, pay_times=[0.25, 0.5, 0.75])   # -> 7.3523

# Zero installment (or an empty schedule) reduces to the plain tree call.
installment_call(S=100, K=100, t=1.0, r=0.05, sigma=0.25,
                 installment=0.0, pay_times=[0.5])                # -> 12.3236
```

A larger installment lowers the upfront value, more payment dates lower it
further, and a prohibitively large installment drives it to zero (the holder
always lapses).

## Double-barrier knock-out

A double knock-out call pays the vanilla call payoff at expiry only if the spot
stays strictly inside a corridor `(L, U)` for the whole life — touching either
the lower barrier `L` or the upper barrier `U` extinguishes it.
`double_knockout_call` prices it in closed form via the Kunitomo-Ikeda (1992)
image series, which enforces both absorbing boundaries and converges geometrically
(a handful of terms is machine-accurate).

```python
from quantforge import double_knockout_call, call_price

double_knockout_call(S=100, K=100, L=80, U=130, t=1.0, r=0.05, sigma=0.25, b=0.05)
# -> 1.9621   (vs vanilla call 12.336 — most of the value is knocked out)
```

The value never exceeds the vanilla call, approaches it as the barriers move far
away, and a tighter corridor lowers it. This is the continuous-monitoring price;
discrete-monitoring Monte Carlo converges down onto it as the step count rises.

The **knock-in** counterpart comes alive only if the spot touches either barrier,
and is priced by the in-out parity `knock-in + knock-out = vanilla`:

```python
from quantforge import double_knockin_call, double_knockout_call, call_price

ki = double_knockin_call(S=100, K=100, L=80, U=130, t=1.0, r=0.05, sigma=0.25, b=0.05)
ko = double_knockout_call(S=100, K=100, L=80, U=130, t=1.0, r=0.05, sigma=0.25, b=0.05)
ki, ko, ki + ko           # -> (10.3739, 1.9621, 12.336 == vanilla call)
```

The knock-in rises toward the vanilla as the corridor tightens (a breach becomes
certain) and falls to zero as the barriers move far away.

## Range-accrual note

A range-accrual note pays a coupon in proportion to the fraction of observation
dates on which a reference index sits inside a band `[L, U]`. By linearity of
expectation the present value of the coupon leg is a discounted sum of GBM range
probabilities `N(d_L) - N(d_U)` across the observation dates — no simulation
needed. `range_accrual_note` returns that present value.

```python
from quantforge import range_accrual_note

# 6% coupon, index in [90, 110], 12 monthly observations over 1y.
range_accrual_note(S=100, L=90, U=110, t=1.0, r=0.05, sigma=0.25,
                   coupon=0.06, observations=12, b=0.05)          # -> 0.026593

# Scale by notional.
range_accrual_note(S=100, L=90, U=110, t=1.0, r=0.05, sigma=0.25,
                   coupon=0.06, observations=12, b=0.05,
                   notional=1_000_000)                            # -> 26592.90
```

The value widens toward the discounted full coupon as the band grows, falls as
volatility rises (the index escapes the band more often), and scales linearly in
the notional.

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

Rates desks quote both a lognormal (Black-76) and a normal (Bachelier) vol for the
same option; convert between them by price-matching:

```python
from quantforge import black_to_normal_vol, normal_to_black_vol

black_to_normal_vol(forward=100, strike=100, t=1.0, sigma_black=0.2)   # -> ~19.97
normal_to_black_vol(forward=100, strike=100, t=1.0, sigma_normal=20.0) # -> ~0.20
```

Near the money `sigma_N ~ sigma_B * F`; the round-trip is exact and both models
reproduce the same price at the converted vol.

Displaced diffusion (shifted lognormal) interpolates between Black-Scholes and
Bachelier and permits negative strikes:

```python
from quantforge import displaced_diffusion_price

displaced_diffusion_price(S=100, K=100, t=1.0, r=0.05, sigma=0.2, shift=50,
                          option_type="call")

# Back out the displaced-diffusion vol from a price (bisection).
from quantforge import displaced_diffusion_implied_vol
displaced_diffusion_implied_vol(target_price=8.5, S=100, K=100, t=1.0, r=0.05,
                                shift=50)
```

The implied-vol solver round-trips with the pricer and, at `shift = 0`, coincides
with the Black-Scholes implied vol.

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

The raw `(a, b, rho, m, s)` are hard to read; the Gatheral-Jacquier jump-wing form
re-expresses a slice in trader quantities — ATM variance, ATM skew, and wing slopes:

```python
from quantforge import raw_to_jumpwing, jumpwing_to_raw

jw = raw_to_jumpwing(params, t)
jw.v, jw.psi, jw.p, jw.c   # ATM variance, ATM skew, left/right wing slopes
jumpwing_to_raw(jw, t)     # exact inverse -> back to raw SVIParams
```

The map is closed-form and round-trips to machine precision, so you can calibrate in
raw space and quote in jump-wing space (or vice versa).

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

From a discount curve, `par_yield` gives the par coupon rate (the par swap rate)
and `par_bond_price` prices any coupon bond off the curve:

```python
from quantforge import par_yield, par_bond_price

par_yield(discount, maturity=10, freq=2)        # par coupon rate for a 10y bond
par_bond_price(discount, coupon_rate=par, maturity=10, freq=2)   # -> 100 at par
```

On a flat curve the par yield equals the flat rate; a bond bearing the par coupon
prices to exactly par, and a richer coupon trades at a premium.

For a fixed-for-floating swap off a single discount curve, `vanilla_swap_value`
values it, `single_curve_par_swap_rate` gives the fair fixed rate, and
`swap_annuity` the PV01:

```python
from quantforge import (vanilla_swap_value, single_curve_par_swap_rate,
                        swap_annuity)

par = single_curve_par_swap_rate(discount, start=0, maturity=5, freq=2)
vanilla_swap_value(discount, par, start=0, maturity=5, freq=2)   # -> 0 at par
swap_annuity(discount, 0, 5, freq=2)                              # PV01 per unit rate
```

The swap is worth zero at the par rate; a receiver is the negative of the payer,
and the value scales linearly in notional.

A single-period forward rate agreement is `fra_value`, with `fra_forward_rate` the
fair rate:

```python
from quantforge import fra_forward_rate, fra_value

fra_forward_rate(discount, t1=1.0, t2=1.5)          # simple forward rate
fra_value(discount, contract_rate=0.03, t1=1.0, t2=1.5, payer=True)
```

The FRA is worth zero at the forward rate; a payer gains when the realized forward
exceeds the contract rate.

## Equity valuation

Cost of capital and discounted-cashflow valuation:

```python
from quantforge import (capm_cost_of_equity, wacc, gordon_growth_value,
                        two_stage_dcf)

ke = capm_cost_of_equity(risk_free=0.03, beta=1.2, market_premium=0.05)
wacc(equity_value=600, debt_value=400, cost_of_equity=ke, cost_of_debt=0.05,
     tax_rate=0.21)
gordon_growth_value(dividend_next=2, discount_rate=0.08, growth=0.03)
two_stage_dcf([100, 110, 121, 133, 146], discount_rate=0.09, terminal_growth=0.03)
```

## Capital budgeting

Project appraisal: NPV, IRR, MIRR, payback period, and profitability index:

```python
from quantforge import npv, irr, profitability_index, payback_period, mirr

cashflows = [-1000, 300, 400, 500, 300]   # outlay then inflows
npv(rate=0.10, cashflows=cashflows)
irr(cashflows)                             # rate where NPV = 0
profitability_index(0.10, cashflows)       # PV inflows / outlay
mirr(cashflows, finance_rate=0.05, reinvest_rate=0.08)
```

## Money-market yields

Discount instruments and the bank-discount / CD / bond-equivalent conventions:

```python
from quantforge import (price_from_discount, bank_discount_yield,
                        money_market_yield, bond_equivalent_yield)

price = price_from_discount(face=100, discount_rate=0.05, days=90)   # T-bill price
bank_discount_yield(100, price, 90)      # actual/360, on face
money_market_yield(100, price, 90)       # actual/360, on price (higher)
bond_equivalent_yield(100, price, 90)    # actual/365, coupon-comparable
```

## Black-Karasinski short rate

A log-normal short-rate model — the log rate mean-reverts, so the rate itself
stays strictly positive (unlike Vasicek/Hull-White). `bk_zero_coupon_bond` prices
a zero on a Hull-White-style trinomial tree in the log rate.

```python
import math
from quantforge import bk_zero_coupon_bond

bk_zero_coupon_bond(r0=0.03, kappa=0.1, theta=math.log(0.03), sigma=0.2, t=5.0)
# -> 0.8564   (5y zero; theta is the log of the target rate level)
```

The price sits near the flat discount factor, falls as the rate or maturity rises,
approaches 1 at short maturity, and falls as volatility rises (Jensen lifts the
expected rate). Rates cannot go negative.

## Hull-White (fitted to a curve)

The Hull-White (extended Vasicek) model adds a time-dependent drift so it reprices
*any* initial discount curve exactly while keeping analytic bond prices.
`hw_zero_from_curve` takes the observed discount curve as a callable and returns
the fitted `P(t, T)`.

```python
import math
from quantforge import hw_zero_from_curve

P0 = lambda T: math.exp(-0.03 * T)          # observed discount curve
hw_zero_from_curve(P0, r0=0.03, a=0.1, sigma=0.01, t=0.0, T=5.0)   # -> 0.860708 (refits P0)
hw_zero_from_curve(P0, r0=0.03, a=0.1, sigma=0.01, t=2.0, T=10.0)  # forward P(2,10)
```

At `t=0` it reproduces the input curve to machine precision (flat, sloped, or the
`a→0` Ho-Lee limit); the `hw_B` helper gives the `B(t,T)` mean-reversion factor.

`hw_bond_option` prices a European option on a zero-coupon bond in closed form
(Jamshidian) off the same curve:

```python
from quantforge import hw_bond_option

# call on the 5y zero, expiring in 1y, struck at the forward
fwd = P0(5.0) / P0(1.0)
hw_bond_option(P0, a=0.1, sigma=0.01, t_option=1.0, t_bond=5.0, strike=fwd, is_call=True)
```

Put-call parity holds exactly (`c - p = P(t_bond) - K·P(t_option)`); the price
rises with volatility and collapses to intrinsic when vol is zero. Bond options are
the building block for caps, floors, and swaptions.

Caps and floors follow directly — each caplet is a scaled bond put:

```python
from quantforge import hw_cap, hw_floor

dates = [1.0, 2.0, 3.0, 4.0, 5.0]     # reset/pay schedule
hw_cap(P0, a=0.1, sigma=0.01, dates=dates, strike=0.03)
hw_floor(P0, a=0.1, sigma=0.01, dates=dates, strike=0.03)
```

The cap is the sum of its caplets and `cap - floor` equals the underlying
fixed-vs-float swap value, so the three are mutually consistent by put-call parity.

`hw_swaption` prices a European payer/receiver swaption by the Jamshidian
decomposition — a portfolio of bond options at the critical short rate:

```python
from quantforge import hw_swaption

hw_swaption(P0, r0=0.03, a=0.1, sigma=0.01, expiry=1.0,
            pay_times=[2, 3, 4, 5], fixed_rate=0.03, is_payer=True)
```

`payer - receiver` equals the forward swap value and an ATM payer equals its
receiver, so the swaption ties back to the cap/floor and bond-option prices by
parity.

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
                   r=0.05, expiry=1)      # single spread

# Refinery crack spread: weighted product basket minus crude (e.g. 3:2:1),
# via the normal-model Bachelier spread (handles a negative margin).
from quantforge import crack_spread_option
crack_spread_option(crude_forward=80, product_forwards=[95, 90],
                    product_weights=[2/3, 1/3], strike=5,
                    sigma_crude=5.0, sigma_products=[7.0, 6.5],
                    corr_products_crude=[0.7, 0.65], r=0.03, expiry=1.0)
```

A single unit-weight product reduces to the plain Bachelier spread; put-call
parity ties the call and put to the forward margin.

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

# Binary payout on a strike breach (discounted breach probability x payout).
from quantforge import degree_day_digital
degree_day_digital(expected_index=1200, strike=1000, sigma=150, r=0.04,
                   expiry=0.5, payout=100000, is_call=True)
```

The call and put digitals sum to the discounted payout (one always pays); an ATM
digital is worth half.

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

Those impact models are theoretical; to measure liquidity from realized price and
volume data, three low-frequency proxies:

```python
from quantforge import roll_spread, amihud_illiquidity, corwin_schultz_spread

roll_spread(prices)                          # 2 sqrt(-cov(dP, dP_lag)): effective spread
amihud_illiquidity(returns, dollar_volumes)  # |return| per dollar traded
corwin_schultz_spread(highs, lows)           # spread from consecutive high-low ranges
```

`roll_spread` recovers a known bid-ask-bounce spread from the negative
autocovariance of price changes; `amihud_illiquidity` rises as volume falls; and
`corwin_schultz_spread` needs only daily highs and lows. They complement the
theoretical Kyle-lambda impact above.

With signed trades, three order-flow measures gauge impact and toxicity directly:

```python
from quantforge import kyle_lambda_regression, order_flow_imbalance, vpin

kyle_lambda_regression(price_changes, signed_volumes)   # empirical Kyle's lambda
order_flow_imbalance(buy_volumes, sell_volumes)         # net signed volume, [-1, 1]
vpin(buy_volumes, sell_volumes)                         # informed-trading proxy, [0, 1]
```

`kyle_lambda_regression` is the OLS slope of price change on signed order flow (it
recovers a known impact coefficient from simulated data); `vpin` (Easley-Lopez de
Prado-O'Hara) is the mean absolute order imbalance across equal-volume buckets —
zero for balanced flow, one for one-sided.

Those measures need each trade signed, which raw tapes rarely provide; three
classifiers infer the aggressor side:

```python
from quantforge import tick_rule, quote_rule, lee_ready

tick_rule(prices)                    # +1 uptick, -1 downtick, carry on a flat tick
quote_rule(prices, bids, asks)       # side of the bid-ask midpoint (0 exactly at it)
lee_ready(prices, bids, asks)        # Lee-Ready hybrid: quote rule, tick tiebreak
```

`lee_ready` is the standard choice: it signs by the quote midpoint away from it and
falls back to the tick rule at the midpoint, so every trade gets a `+1 / -1` you can
feed straight into `vpin` or `order_flow_imbalance`.

Signed trades and quotes also decompose the trading cost into what the liquidity
provider keeps and the permanent price move:

```python
from quantforge import (quoted_spread, effective_spread, realized_spread,
                        price_impact)

quoted_spread(bids, asks)                              # posted (ask - bid) / mid
effective_spread(prices, mids, signs)                  # cost actually paid vs mid
realized_spread(prices, mids, future_mids, signs)      # provider's kept portion
price_impact(prices, mids, future_mids, signs)         # permanent move
```

The decomposition satisfies `effective = realized + price_impact` term by term:
`future_mids` is the midpoint a short horizon after each trade, so the realized
spread nets out the permanent move that the price impact captures.

## Hawkes self-exciting process

Trades and order arrivals cluster: each event lifts the chance of the next. A
Hawkes process with an exponential kernel captures that, with intensity
`lambda(t) = mu + sum alpha exp(-beta (t - t_i))`:

```python
from quantforge import (hawkes_intensity, hawkes_branching_ratio,
                        hawkes_log_likelihood, hawkes_simulate, hawkes_fit)

events = hawkes_simulate(mu=0.5, alpha=0.9, beta=2.0, t_max=5000)
fit = hawkes_fit(events)
fit["mu"], fit["alpha"], fit["beta"], fit["branching_ratio"]   # recovered params
hawkes_intensity(t=10.0, history=events, mu=0.5, alpha=0.9, beta=2.0)
```

The branching ratio `alpha / beta` is the expected number of events triggered by
each event; below 1 the process is stationary with mean rate `mu / (1 - alpha/beta)`,
which the simulator reproduces. `hawkes_fit` maximizes the exact `O(n)` recursive
log-likelihood, and `hawkes_simulate` uses Ogata thinning.

To check a fit, the time-rescaling theorem turns a correct model's events into
unit-rate Poisson arrivals; `hawkes_residuals` computes those rescaled inter-event
times (i.i.d. `Exp(1)` when the model holds) and `hawkes_gof_test` runs a
Kolmogorov-Smirnov test against `Exp(1)`:

```python
from quantforge import hawkes_residuals, hawkes_gof_test

res = hawkes_residuals(events, fit["mu"], fit["alpha"], fit["beta"])   # ~ Exp(1)
D, p = hawkes_gof_test(events, fit["mu"], fit["alpha"], fit["beta"])   # small p rejects
```

A well-specified model gives residuals with mean one and a large p-value; wrong
parameters push the p-value to zero.

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

A Phoenix autocallable — conditional memory coupons, early autocall, and
down-and-in downside — is priced by Monte Carlo:

```python
from quantforge import phoenix_autocall_mc

phoenix_autocall_mc(S=100, t=3, r=0.03, sigma=0.25, observation_times=[1, 2, 3],
                    autocall_barrier=110, coupon_barrier=80, coupon=0.08,
                    protection_barrier=70, memory=True)
```

The memory feature (paying missed coupons on the next barrier touch) raises the
value; a lower coupon barrier pays more often.

## Survival analysis (Kaplan-Meier / Nelson-Aalen)

Nonparametric survival and hazard from right-censored data — defaults, lapses, or
any time-to-event where some observations are cut short:

```python
from quantforge import kaplan_meier, nelson_aalen, survival_at

times  = [2, 3, 3, 5, 7]
events = [1, 1, 0, 1, 1]          # 0 = right-censored
t, s = kaplan_meier(times, events)      # product-limit survival S(t)
t, h = nelson_aalen(times, events)      # cumulative hazard H(t)
survival_at(times, events, query=4)     # S(4)
```

Censored subjects stay in the risk set until their censoring time, then drop out
without an event. With no censoring the Kaplan-Meier curve is just `1 - ECDF`; the
Nelson-Aalen hazard relates to it by `S(t) ~ exp(-H(t))`.

To test whether two groups' survival differs, `log_rank_test` (Mantel-Cox)
accumulates observed-minus-expected events across the pooled timeline:

```python
from quantforge import log_rank_test

chi2, p = log_rank_test(times_a, events_a, times_b, events_b)   # small p: curves differ
```

It returns a chi-square(1) statistic and p-value — the standard test for equal
hazards between two censored samples, and symmetric in the two groups.

Two summaries condense a curve into a number: `median_survival_time` (where survival
crosses 0.5) and `restricted_mean_survival_time` (area under the curve up to a
horizon):

```python
from quantforge import median_survival_time, restricted_mean_survival_time

median_survival_time(times, events)                    # None if it never reaches 0.5
restricted_mean_survival_time(times, events, tau=5.0)  # expected time capped at tau
```

RMST is defined even when the tail is censored (unlike the plain mean), never
exceeds `tau`, and grows with it; the median can be `None` under heavy censoring.

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

For aggregate losses, `panjer_poisson` gives the exact compound-Poisson
aggregate-loss distribution (Poisson claim count, discrete severity), from which
`stop_loss_premium` and `layer_expected_loss` price stop-loss and excess-of-loss
reinsurance:

```python
from quantforge import (panjer_poisson, aggregate_mean, stop_loss_premium,
                        layer_expected_loss)

g = panjer_poisson(lam=3.0, severity_pmf=[0.0, 0.4, 0.6])   # P(S = k), sums to 1
aggregate_mean(g)                       # = lam * E[X]
stop_loss_premium(g, retention=5)       # E[max(S - 5, 0)]
layer_expected_loss(g, attachment=3, limit=2)   # excess-of-loss layer cost
```

The mean and variance match `lam E[X]` and `lam E[X^2]`; a full-width layer equals
the mean. When claims cluster (over-dispersed frequency), use
`panjer_negative_binomial` instead — same recursion, heavier aggregate tail at the
same mean:

```python
from quantforge import panjer_negative_binomial

panjer_negative_binomial(size=5, prob=0.5, severity_pmf=[0.0, 0.4, 0.6])
```

Price the loss with a risk load using a distortion measure — the Wang transform or
proportional hazard:

```python
from quantforge import wang_premium, proportional_hazard_premium

wang_premium(g, lam=0.5)                 # Phi(Phi^-1(S) + lam) distortion
proportional_hazard_premium(g, rho=2.0)  # S^{1/rho} distortion
```

Both equal the expected loss at zero distortion (`lam = 0`, `rho = 1`) and load the
tail above it as the parameter grows.

Capital measures come off the same aggregate grid — `aggregate_var` (loss quantile)
and `aggregate_tvar` (Tail-VaR / CTE):

```python
from quantforge import aggregate_var, aggregate_tvar

aggregate_var(g, confidence=0.99)     # 99% VaR
aggregate_tvar(g, confidence=0.99)    # expected loss beyond the 99% VaR
```

TVaR is always at least the VaR and both rise with the confidence level.

For experience rating, `buhlmann_premium` blends a risk's own mean with the
collective mean by the credibility factor `Z = n / (n + k)`, `k = EPV / VHM`:

```python
from quantforge import buhlmann_premium, credibility_factor, buhlmann_straub_premium

buhlmann_premium(own_mean=500, collective_mean=400, n=10, epv=100, vhm=20)
buhlmann_straub_premium(claims=[50, 50, 50], exposures=[10, 10, 10],
                        collective_mean=4.0, epv=100, vhm=20)   # unequal exposures
```

More data or a larger between-risk spread pulls `Z` toward 1 (trust the
individual); more within-risk noise pulls it toward 0 (trust the collective).

`chain_ladder` reserves outstanding claims from a cumulative run-off triangle via
volume-weighted development factors:

```python
from quantforge import chain_ladder

tri = [[100, 150, 180], [110, 165], [120]]      # cumulative claims by accident year
r = chain_ladder(tri)
r["factors"]         # age-to-age: [1.5, 1.2]
r["ultimate"]        # [180, 198, 216]
r["total_reserve"]   # 129  (IBNR = ultimate - paid)
```

A fully-developed accident year carries zero reserve, and each ultimate is at
least the latest paid amount.

`bornhuetter_ferguson` blends the development pattern with an a-priori ultimate
(reserve = a-priori × undeveloped fraction) — more stable for green years:

```python
from quantforge import bornhuetter_ferguson

bornhuetter_ferguson(tri, apriori_ultimates=[200, 200, 200])["reserve"]
```

Feeding the chain-ladder ultimate as the a-priori reproduces the chain-ladder
reserves exactly, so the two methods are consistent endpoints of the same blend.

`cape_cod` estimates the a-priori loss ratio from the data instead of assuming it —
total losses over premium weighted by development:

```python
from quantforge import cape_cod

r = cape_cod(tri, premiums=[300, 300, 300])
r["elr"]              # fitted expected loss ratio
r["total_reserve"]
```

Cape Cod equals Bornhuetter-Ferguson with a `premium x ELR` a-priori, so the three
reserving methods share one development pattern and are mutually consistent.

`mack_standard_error` attaches an uncertainty to the chain-ladder point estimate:
Mack's (1993) distribution-free standard error of the reserve. It returns the
per-accident-year and total mean-squared error of prediction — process error plus
estimation error, with the between-year correlation in the total — and the
coefficients of variation:

```python
from quantforge import mack_standard_error

wtri = [[100, 180, 230, 245], [120, 200, 260], [90, 175], [130]]
r = mack_standard_error(wtri)
r["total_reserve"]      # 272.02, matches chain_ladder
r["total_std_error"]    # 27.24, standard error of that reserve
r["total_cv"]           # 0.10, std_error / reserve
```

On the Taylor-Ashe triangle from Mack's paper it reproduces his figures to the
dollar: a reserve of 18,680,856 with a standard error of 2,447,095. The total
standard error exceeds the root of the summed per-year variances (positive
correlation) but stays below their plain sum.

Triangles supplied in either form convert with `incremental_to_cumulative` /
`cumulative_to_incremental`, and `paid_to_date` reads the latest diagonal:

```python
from quantforge import incremental_to_cumulative, paid_to_date

cum = incremental_to_cumulative([[100, 50, 30], [110, 55], [120]])
paid_to_date(cum)      # [180, 165, 120]
```

To capture development past the last observed age, extrapolate a tail factor and
apply it:

```python
from quantforge import exponential_tail_factor, chain_ladder_with_tail

tail = exponential_tail_factor([1.5, 1.2, 1.1, 1.05])   # fits decaying excess-over-1
chain_ladder_with_tail(cum, tail)["total_reserve"]
```

A unit tail reproduces plain chain-ladder; a decaying factor pattern gives a
finite tail above 1, and a non-decaying one is rejected.

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

The large-pool loss distribution is complete: `vasicek_loss_cdf` and
`vasicek_loss_quantile` are joined by the density and the expected shortfall (the
average loss beyond the quantile, Tasche's closed form):

```python
from quantforge import vasicek_loss_pdf, vasicek_loss_expected_shortfall

vasicek_loss_pdf(0.03, pd=0.02, rho=0.15)                        # loss density
vasicek_loss_expected_shortfall(q=0.99, pd=0.02, rho=0.15)       # 0.136, >= the VaR
```

The density integrates to one with mean `pd`; the expected shortfall sits at or
above `vasicek_loss_quantile` and rises with both the confidence level and the
asset correlation.

## kth-to-default basket

The full number-of-defaults distribution of a homogeneous basket under the
one-factor Gaussian copula, and the trigger probability of a kth-to-default swap.

```python
from quantforge import basket_default_distribution, kth_to_default_probability

kth_to_default_probability(n=10, k=1, pd=0.05, rho=0.3)   # first-to-default: 0.307
kth_to_default_probability(n=10, k=5, pd=0.05, rho=0.3)   # fifth-to-default: 0.008
basket_default_distribution(10, 0.05, 0.3)                # P(exactly k defaults), sums to 1
```

The trigger falls as `k` rises; correlation clusters defaults, so senior
(high-`k`) triggers rise with `rho` while first-to-default falls. The mean number
of defaults is `n * pd`, and `rho = 0` recovers `1 - (1-pd)^n` for first-to-default.

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

## OLS regression

General multivariate ordinary least squares with the standard diagnostics —
coefficient standard errors, t-statistics, two-sided p-values and confidence
intervals, R-squared, adjusted R-squared, and the overall F-statistic with its
p-value. An intercept column is added by default.

```python
from quantforge import ols_fit

X = [[1, 2], [2, 1], [3, 4], [4, 3], [5, 5]]
y = [5, 4, 10, 9, 13]
m = ols_fit(X, y, confidence=0.95)
m["coefficients"]     # [intercept, b1, b2]
m["t_stats"]          # significance of each coefficient
m["p_values"]         # two-sided p-value per coefficient
m["conf_int"]         # [low, high] per coefficient at `confidence`
m["r_squared"], m["adj_r_squared"]
m["f_stat"], m["f_pvalue"]   # overall model significance
```

The p-values and intervals come from the same t and F distributions in the library:
for a single regressor the overall F equals the slope t-squared and shares its
p-value, and the intervals widen with the confidence level. An exact linear
relationship gives R-squared 1 and zero residuals; adjusted R-squared never exceeds
R-squared. For the finance-specific alpha/beta return regression, see the next
section.

Those standard errors assume homoskedastic, serially-uncorrelated errors. When that
fails — volatility clustering, autocorrelated residuals — use the sandwich estimator:
`white_hc0` for heteroskedasticity and `newey_west` for heteroskedasticity *and*
autocorrelation:

```python
from quantforge import white_hc0, newey_west

w = white_hc0(X, y)                  # HC0 robust SEs
w["coefficients"], w["std_errors"], w["t_stats"]

nw = newey_west(X, y, lags=8)        # Bartlett-weighted HAC SEs
nw["std_errors"]                     # valid under autocorrelated residuals
```

Both return the same coefficients as `ols_fit` (only the standard errors change).
`newey_west` with `lags=0` reduces exactly to `white_hc0`. Under heteroskedastic errors
the White SE is larger than the (invalid) OLS SE; under AR(1) residuals the Newey-West
SE grows larger still, which is the correction that keeps the t-statistics honest —
the covariance is `(X'X)^{-1} S (X'X)^{-1}` with the meat `S` built from the residual
score vectors.

Robust standard errors keep the *inference* honest but leave the OLS point estimate
inefficient. When you know the error structure, weighting restores efficiency:
`weighted_least_squares` down-weights noisy observations, and
`generalized_least_squares` handles a full error covariance:

```python
from quantforge import weighted_least_squares, generalized_least_squares

# weights proportional to 1 / var(eps_t): noisy points count less
w = weighted_least_squares(X, y, weights=[1 / v for v in error_variances])
w["coefficients"], w["std_errors"], w["r_squared"]

# known error covariance Sigma (e.g. AR(1) residuals): whitened by its Cholesky factor
g = generalized_least_squares(X, y, cov=Sigma)
g["coefficients"], g["std_errors"]
```

Equal weights make `weighted_least_squares` reproduce OLS exactly; a diagonal `cov`
makes `generalized_least_squares` coincide with WLS using `1 / diag(cov)` weights, and
an identity `cov` reduces to OLS. GLS whitens the system with the Cholesky factor of
`Sigma` and then runs OLS on the transformed data, so with the correct covariance it is
the minimum-variance linear unbiased estimator.

For streaming data — or a fit that must adapt as it goes — `RecursiveLeastSquares`
updates the coefficients one observation at a time without re-solving, and a forgetting
factor lets it track slowly-varying coefficients:

```python
from quantforge import RecursiveLeastSquares, recursive_least_squares

# batch wrapper: with forgetting = 1 this matches OLS on the same design
recursive_least_squares([[1.0, 1.0], [1.0, 2.0], [1.0, 3.0]], [3, 5, 7])   # [1.0, 2.0]

# streaming, adaptive: forgetting < 1 down-weights old data
rls = RecursiveLeastSquares(n_features=2, forgetting=0.98)
for x_row, y_t in stream:
    rls.update(x_row, y_t)      # x_row includes a leading 1.0 for the intercept
    rls.beta                    # current coefficient estimate
rls.predict([1.0, x_new])
```

Each `update` is `O(k²)` via the Sherman-Morrison identity — no matrix re-inversion and
no need to store history. With `forgetting = 1` it converges to the exact batch OLS fit;
with `forgetting < 1` it down-weights the past geometrically, so it re-adapts within a
few dozen points when the underlying slope shifts (adaptive filtering, time-varying
betas).

When the predictors are collinear (correlated factors, an over-parameterized curve),
`principal_components_regression` rotates them onto their principal components, keeps
the top `k`, and regresses on those — dropping the low-variance directions that make
OLS unstable:

```python
from quantforge import principal_components_regression

pcr = principal_components_regression(X, y, n_components=2)
pcr["coefficients"], pcr["intercept"]     # in the original predictor space
pcr["explained_variance"]                 # fraction of predictor variance kept
```

Retaining all components reproduces OLS exactly; keeping fewer trades a little bias for
much lower variance and stays well-behaved even when two predictors are nearly
identical (where OLS coefficients blow up). The `explained_variance` field reports how
much predictor variance the retained components capture, so you can pick `k` from a
scree-style cutoff. The coefficients are mapped back to the original variables, so they
plug into the same prediction as any other fit.

When regressors are collinear or numerous, `ridge_regression` adds an L2 penalty
that shrinks the slopes and keeps the system solvable:

```python
from quantforge import ridge_regression

ridge_regression(X, y, alpha=0.0)     # reproduces OLS
ridge_regression(X, y, alpha=10.0)    # slopes shrunk toward zero (intercept kept)
```

A larger `alpha` shrinks the slope coefficients more (the intercept is not
penalized); ridge stays solvable even under perfect collinearity, where OLS is
singular.

`bayesian_linear_regression` gives ridge a probabilistic reading: with a Gaussian prior
it returns a full posterior over the coefficients, so you get credible intervals and a
predictive variance the point estimate can't express:

```python
from quantforge import bayesian_linear_regression, bayesian_predict

m = bayesian_linear_regression(X, y, alpha=1.0, beta_noise=4.0)
m["mean"], m["std"]                 # posterior coefficient means and std devs
mean, var = bayesian_predict(m, x_new)   # predictive mean and variance at a new point
```

`alpha` is the prior precision (shrinkage strength) and `beta_noise` the noise precision;
the posterior mean equals ridge with `lambda = alpha / beta_noise`, a weak prior recovers
OLS, and a stronger prior shrinks the coefficients further. The predictive variance is
`1/beta_noise + xᵀΣx` — observation noise plus coefficient uncertainty — so it widens
where data is sparse or when extrapolating, flagging where the model is guessing.

Where ridge shrinks every slope, `lasso_regression` (L1 penalty, coordinate
descent) drives some to exactly zero — so it also selects features:

Where ridge shrinks every slope, `lasso_regression` (L1 penalty, coordinate
descent) drives some to exactly zero — so it also selects features:

```python
from quantforge import lasso_regression

lasso_regression(X, y, alpha=0.0)    # reproduces OLS
lasso_regression(X, y, alpha=0.1)    # irrelevant features set to exactly 0
```

`alpha = 0` is OLS; raising it zeros out the weakest predictors first (feature
selection) and finally collapses the fit to `mean(y)`. The intercept is never
penalized and coefficients come back on the original scale.

`elastic_net` blends the two penalties (`l1_ratio` sets the mix), keeping LASSO's
feature selection while sharing weight across correlated predictors like ridge:

```python
from quantforge import elastic_net

elastic_net(X, y, alpha=0.1, l1_ratio=1.0)    # == lasso_regression
elastic_net(X, y, alpha=0.1, l1_ratio=0.5)    # L1 + L2 mix
```

With `l1_ratio = 1` it is pure LASSO; lowering it adds L2 shrinkage, which on
near-duplicate features spreads the weight between them (the grouping effect) rather
than arbitrarily dropping one.

To fit a conditional *quantile* rather than the mean, `quantile_regression`
minimizes the pinball loss (via iteratively-reweighted least squares):

```python
from quantforge import quantile_regression

quantile_regression(X, y, tau=0.5)    # median (least-absolute-deviations) fit
quantile_regression(X, y, tau=0.9)    # upper conditional-tail line
```

`tau = 0.5` is the robust median regression; a higher `tau` traces the upper tail.
The fit has the defining quantile property — a fraction `tau` of the points fall
below the line — and minimizes the same pinball loss that `pinball_loss` scores.

## Logistic regression

Binary classification / default-probability modeling by iteratively reweighted
least squares. `fit_logistic` maximizes the logistic likelihood via Newton-Raphson
(ridge-stabilized against separation); `predict_proba` returns class probabilities.

```python
from quantforge import fit_logistic, predict_proba

m = fit_logistic(X, y)          # y is binary 0/1
m["coefficients"]                # [intercept, b1, b2, ...]
predict_proba(m, [[1, 0], [-1, 0]])   # -> [0.899, 0.275] class-1 probabilities
```

It recovers the generating logit coefficients on simulated data, returns
probabilities strictly in (0, 1), and classifies separable data near-perfectly. A
positive coefficient makes the probability rise monotonically in that feature.

For *count* responses — event frequencies, claim counts, arrival rates —
`poisson_regression` is the matching GLM: it models `E[y|x] = exp(x'beta)` with a log
link (so the rate is always positive) and fits by IRLS:

```python
from quantforge import poisson_regression, poisson_predict

m = poisson_regression(X, y)      # y is non-negative counts
m["coefficients"]                  # [intercept, b1, ...] on the log-rate scale
poisson_predict(m, X_query)        # predicted rates exp(x'beta), always > 0
```

It recovers a known log-linear rate from sampled counts, the intercept-only fit gives
`exp(intercept) = ` the mean count, and every prediction is positive. Use it instead of
OLS whenever the response is a count whose variance grows with its mean — the Gaussian
assumption breaks there, but the Poisson mean-variance link is built in.

## Classification metrics

Evaluate a probabilistic classifier's scores against binary labels. `roc_auc` is
the rank AUC; `precision_recall_f1` and `confusion_matrix` summarize a threshold;
`log_loss` and `brier_score` score the probabilities directly.

```python
from quantforge import (roc_auc, precision_recall_f1, confusion_matrix,
                        log_loss, brier_score)

roc_auc(y_true, y_score)                 # 1 perfect, 0.5 random
precision_recall_f1(y_true, y_score, threshold=0.5)   # (p, r, f1)
confusion_matrix(y_true, y_score, 0.5)                # (tp, fp, fn, tn)
log_loss(y_true, y_score)                # cross-entropy, 0 = perfect
brier_score(y_true, y_score)             # mean squared prob error, 0 = perfect
```

AUC uses the Mann-Whitney rank statistic (ties count as half); the Brier score is
0 for exact probabilities and 0.25 for all-0.5 guesses. Pair these with
`fit_logistic` / `predict_proba` above.

## Forecast calibration (Brier decomposition)

A single Brier score conflates two very different failings. Murphy's decomposition
splits it into `reliability - resolution + uncertainty`: *reliability* is how far
each forecast group's observed frequency drifts from the forecast (0 = perfectly
calibrated, lower is better), *resolution* rewards forecasts that separate outcomes
away from the base rate (higher is better), and *uncertainty* is the irreducible
`obar (1 - obar)` variance of the outcome.

```python
from quantforge import brier_decomposition

f = [0.2, 0.2, 0.8, 0.8, 0.6]
o = [0,   1,   1,   1,   0]

d = brier_decomposition(f, o)      # grouped by identical value -> exact identity
d["reliability"]   # 0.124
d["resolution"]    # 0.140
d["uncertainty"]   # 0.240   (base rate 0.6 -> 0.6 * 0.4)
d["brier"]         # 0.224 == reliability - resolution + uncertainty == raw Brier
```

With `n_bins=None` (the default) forecasts are grouped by identical value and the
identity reconstructs the raw Brier score exactly; pass an integer `n_bins` to bin
continuous forecasts (the reconstruction is then approximate). The
`reliability_curve` returns the calibration diagram and `expected_calibration_error`
its scalar summary:

```python
from quantforge import reliability_curve, expected_calibration_error

reliability_curve(forecasts, outcomes, n_bins=5)   # [(mean_forecast, observed, count), ...]
expected_calibration_error(forecasts, outcomes, n_bins=5)   # count-weighted |gap|
```

A well-calibrated forecaster's reliability curve sits on the `observed == forecast`
diagonal, so its ECE is near zero (about 0.005 on a large calibrated sample) while a
forecaster that always says 0.9 when the truth is 0.5 lands near 0.4. Feed the raw
scores through `isotonic_fit` (above) to recalibrate, then re-measure here.

Where isotonic calibration assumes only monotonicity, **Platt scaling** assumes a
sigmoidal miscalibration and fits a two-parameter logistic
`P = 1 / (1 + exp(-(A s + B)))` — far more data-efficient, and the right choice when
the raw scores are roughly logit-shaped but mis-scaled (an overconfident classifier):

```python
from quantforge import platt_fit, platt_predict, platt_calibrate

A, B = platt_fit([-2, -1, 0, 1, 2], [0, 0, 0, 1, 1])
platt_predict([-5, 0, 5], A, B)      # [0.0124, 0.3968, 0.9719]

(params, predict) = platt_calibrate(train_scores, train_labels)
predict(new_scores)                  # calibrated probabilities in [0, 1]
```

`platt_fit` minimizes the regularized logistic loss with Newton's method on Platt's
smoothed targets (guarding the tails), returning `(A, B)`; `A > 0` means probability
rises with the score. `platt_calibrate` bundles the fit and a `predict` closure, or
applies directly if you pass `new_scores`. Choose Platt when data is scarce and the
distortion is smooth; choose `isotonic_fit` when you have enough data and the
distortion may be an arbitrary monotone shape.

## Cross-validation

Model-agnostic out-of-sample evaluation. `k_fold_indices` yields disjoint
train/test index pairs; `train_test_split` does a single split;
`cross_val_score` runs a caller-supplied fit/score over the folds.

```python
from quantforge import k_fold_indices, train_test_split, cross_val_score

for train_idx, test_idx in k_fold_indices(len(y), k=5, shuffle=True, seed=42):
    ...                                    # fit on train_idx, evaluate on test_idx

scores = cross_val_score(X, y, fit_fn, score_fn, k=5)   # one score per fold
```

The test folds tile the data exactly (disjoint, covering every index once) with
sizes differing by at most one; a seeded shuffle is reproducible. Pair with
`fit_logistic`/`ols_fit` and the metrics above for a full evaluation loop.

For an in-sample penalty instead of cross-validation, the information criteria trade
fit against complexity — the model with the lowest value wins:

```python
from quantforge import gaussian_log_likelihood, aic, aicc, bic, hqic

ll = gaussian_log_likelihood(rss=12.5, n=50)   # from a least-squares fit
aic(ll, k=4); bic(ll, k=4, n=50)               # k = params incl. noise variance
aicc(ll, 4, 50); hqic(ll, 4, 50)
```

`aic` rewards fit lightly; `bic` penalizes complexity harder for large `n` (favoring
parsimony); `aicc` corrects AIC when `n` is small; `hqic` sits between. Feed a
least-squares RSS through `gaussian_log_likelihood` first.

## Feature scaling

Fit a scaler on training data, apply it to test data — no leakage.
`fit_standardize` (z-score), `fit_min_max` ([0,1]), and `fit_robust` (median/IQR)
return params; `scale_transform` applies them and `scale_inverse_transform` undoes
them.

```python
from quantforge import fit_standardize, scale_transform, scale_inverse_transform

params = fit_standardize(X_train)        # per-column mean and std
Z_train = scale_transform(params, X_train)
Z_test = scale_transform(params, X_test)   # same params — no test-set leakage
X_back = scale_inverse_transform(params, Z_train)   # round-trips exactly
```

Standardized columns have mean 0 and std 1; min-max maps to [0,1]; the robust
scaler uses the median and IQR, so its center is unmoved by outliers. All three
round-trip through the inverse.

To let a linear model fit curvature and interactions, expand the features first
with `polynomial_features`:

```python
from quantforge import polynomial_features

F, terms = polynomial_features(X, degree=2)   # 1, x1, x2, x1^2, x1*x2, x2^2, ...
polynomial_features(X, degree=2, interaction_only=True)   # drops pure powers
```

It returns the expanded matrix and the index-tuple for each column; the feature
count is `C(p + d, d)`. Feed the result straight into `ols_fit`, `ridge_regression`,
or `fit_logistic`.

Categorical columns need encoding first. `fit_label_encoder` maps categories to
integer codes; `one_hot_encode` expands them to 0/1 indicator columns.

```python
from quantforge import fit_label_encoder, label_encode, one_hot_encode

enc = fit_label_encoder(["b", "a", "c", "a"])   # categories sorted -> codes
label_encode(enc, ["a", "c"])                    # -> [0, 2]
one_hot_encode(enc, ["a", "c"])                  # -> [[1,0,0], [0,0,1]]
```

Fit the encoder on the training categories and reuse it on test data; unseen
categories map to `-1` (label) or an all-zero row (one-hot).

## k-nearest neighbors

A lazy, non-parametric baseline for classification and regression: predict a query
from its `k` closest training points. `knn_classify` votes on their labels;
`knn_regress` averages their targets.

```python
from quantforge import knn_classify, knn_regress

knn_classify(X_train, y_train, X_query, k=3)   # majority vote of the k nearest
knn_regress(X_train, y_train, X_query, k=3)    # mean target of the k nearest
```

With `k=1` it reproduces the training labels (or the nearest target) exactly; a
larger `k` smooths the decision boundary. No training step — the model is the data.

## Gaussian naive Bayes

A fast probabilistic classifier that assumes features are conditionally
independent and normally distributed within each class. `fit_gaussian_nb`
estimates the priors and per-feature moments; `predict_gaussian_nb` and
`predict_proba_gaussian_nb` classify and score queries.

```python
from quantforge import (fit_gaussian_nb, predict_gaussian_nb,
                        predict_proba_gaussian_nb)

m = fit_gaussian_nb(X, y)                       # y = class labels (any hashable)
predict_gaussian_nb(m, [[1.1, 1.0]])            # -> ['lo']
predict_proba_gaussian_nb(m, [[1, 1]])          # -> [{'lo': 1.0, 'hi': 0.0}]
```

Posteriors are computed in log space and softmaxed, so they sum to 1 and stay
stable in high dimension. Despite the independence assumption it is a strong,
cheap baseline.

## Decision stump

A depth-1 decision tree: the single feature and threshold that best split the
classes by Gini impurity. Interpretable on its own and the base learner of
boosting.

```python
from quantforge import fit_decision_stump, predict_decision_stump, gini_impurity

s = fit_decision_stump(X, y)
s["feature"], s["threshold"], s["gini"]     # the chosen split and its weighted Gini
predict_decision_stump(s, [[1.2], [8.8]])   # -> [0, 1]
gini_impurity([0, 0, 1, 1])                  # -> 0.5
```

It scans every feature and midpoint for the split minimizing the size-weighted
child Gini; a cleanly separable set gives Gini 0, and identical rows fall back to
the majority class.

Grow it into a full **CART tree** with `fit_decision_tree` — recursive splits up
to a depth or node-size limit:

```python
from quantforge import fit_decision_tree, predict_decision_tree, tree_depth

t = fit_decision_tree(X, y, max_depth=5)
predict_decision_tree(t, X_query)
tree_depth(t)
```

A depth-1 tree is exactly the stump; deeper trees solve non-linearly separable
problems (e.g. XOR) that a single split cannot, at the cost of overfitting — cap
`max_depth` or `min_samples`, and validate with `cross_val_score`.

To cut a single tree's variance, bag many on bootstrap resamples and vote — a
**random forest**:

```python
from quantforge import fit_random_forest, predict_random_forest

f = fit_random_forest(X, y, n_trees=25, max_depth=6)
predict_random_forest(f, X_query)
```

The forest generalizes at least as well as any one tree on held-out data and is
reproducible for a fixed seed.

The same recursive splitting works for a *numeric* target: `fit_regression_tree` splits
to minimize squared error and predicts the mean of each leaf, a piecewise-constant
regressor for nonlinear relationships:

```python
from quantforge import fit_regression_tree, predict_regression_tree

X = [[float(i)] for i in range(10)]
y = [0.0 if i < 5 else 10.0 for i in range(10)]
t = fit_regression_tree(X, y, max_depth=3)
predict_regression_tree(t, [[2.0], [7.0]])   # [0.0, 10.0] — step recovered exactly
```

It recovers a step function exactly, collapses a constant target to one leaf, and its
error falls monotonically as `max_depth` grows (deeper = finer piecewise fit, so cap
depth and validate to avoid overfitting). A leaf predicts the mean of the training
targets that reach it, so the fit is a step function — pair it with bagging for a smooth
ensemble, the regression analogue of the random forest above.

For a stronger regressor, `fit_gradient_boost` boosts *shallow* trees additively: each
new tree fits the residuals of the ensemble so far, corrected by a shrinkage learning
rate:

```python
from quantforge import fit_gradient_boost, predict_gradient_boost

m = fit_gradient_boost(X, y, n_estimators=100, learning_rate=0.1, max_depth=3)
predict_gradient_boost(m, X_query)
```

With squared-error loss the residual is the negative gradient, so this is gradient
descent in function space — a hundred depth-3 trees fit a smooth nonlinearity far more
accurately than one deep tree while resisting overfit through the small `learning_rate`.
Training error falls monotonically as `n_estimators` grows; a `learning_rate` of 0
leaves the prediction at the target mean. The shallow-tree-plus-small-rate combination
is the classic bias-variance sweet spot — validate `n_estimators` and `learning_rate`
together, since more trees demand a smaller rate to avoid overfitting.

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

Where `factor_regression` allows any (signed, unbounded) loadings, Sharpe's
returns-based `style_analysis` constrains them to a long-only, fully-invested mix —
the implied asset-class weights of a fund from its returns alone:

```python
from quantforge import style_analysis

res = style_analysis(fund_returns, [large_cap, small_cap, bonds])
res["weights"]           # non-negative, sum to 1: the implied style mix
res["r_squared"]         # variance explained by the style
res["tracking_error"]    # selection (active) return the mix cannot explain
```

The simplex constraint makes the weights interpretable as portfolio holdings; the
unexplained residual is the manager's selection return.

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

Nonparametric confidence intervals: IID, moving-block and stationary block
bootstraps for serially-correlated data, BCa, and the jackknife:

```python
from quantforge import (bootstrap_ci, stationary_bootstrap_ci,
                        moving_block_bootstrap_ci, bca_bootstrap_ci,
                        jackknife_estimate)

returns = [0.01, -0.02, 0.03, 0.00, 0.015, -0.01, 0.025]
bootstrap_ci(returns, confidence=0.95)          # (lower, point, upper) for the mean
bca_bootstrap_ci(returns)                        # bias-corrected accelerated
stationary_bootstrap_ci(returns, mean_block=3)   # Politis-Romano random-length blocks
moving_block_bootstrap_ci(returns, block=3)      # Kunsch fixed-length blocks
jackknife_estimate(returns)                      # (estimate, standard_error)
```

The IID bootstrap under-covers autocorrelated data; both block bootstraps widen
the interval to restore coverage (with `block = 1` the moving-block reduces to the
IID interval).

Where the bootstrap builds a confidence interval, a permutation test gives an
assumption-free *p-value* for a group difference — shuffling the labels to build the
null directly. `permutation_test` handles two independent samples and any statistic;
`paired_permutation_test` handles matched pairs:

```python
from quantforge import permutation_test, paired_permutation_test

permutation_test(group_a, group_b)                       # difference in means, two-sided
permutation_test(a, b, statistic=lambda a, b: median(a) - median(b))
paired_permutation_test(before, after, alternative="greater")
```

The default statistic is the difference in means, but any `statistic(a, b)` works
(median difference, correlation, a custom loss) — the test is valid for all of them.
`permutation_test` reshuffles the pooled labels; `paired_permutation_test` flips the
sign of each within-pair difference. Both are seeded (deterministic), report a p-value
that is calibrated under the null, and need no distributional assumption — the
significance companion to the bootstrap's intervals.

## Hodrick-Prescott filter

Split a time series into a smooth trend and a cyclical residual by trading fit
against the curvature of the trend. `hp_filter` solves the linear system
`(I + lambda D'D) tau = y` with an O(n) banded factorization — no dense inverse.

```python
from quantforge import hp_filter

y = [100, 101, 103, 102, 105, 107, 106, 109, 111, 110, 113, 115]
trend, cycle = hp_filter(y, lam=1600)     # 1600 is the standard quarterly lambda
trend[:2]                                  # -> [99.57, 100.89]  (smooth)
cycle[:2]                                  # -> [0.425, 0.109]   (y - trend)
```

`trend + cycle` reconstructs `y` exactly. A larger `lambda` gives a smoother
trend (`lambda -> inf` collapses to the least-squares straight line); `lambda = 0`
returns the data untouched. Common choices: 1600 (quarterly), 129600 (monthly),
6.25 (annual).

The Savitzky-Golay filter smooths (or differentiates) by fitting a local polynomial
in a sliding window, preserving peak shape better than a moving average:

```python
from quantforge import savgol_filter, savgol_coeffs

savgol_filter(series, window=11, degree=3)            # smoothed series
savgol_filter(series, window=11, degree=3, deriv=1)   # smoothed first derivative
savgol_coeffs(11, 3)                                  # the convolution weights
```

It reproduces polynomials up to `degree` exactly, so it introduces no bias on
locally-polynomial signals, and the `deriv` argument gives a smoothed numerical
derivative — useful for estimating slopes or curvature from noisy data.

Where Savitzky-Golay uses a fixed window on evenly-spaced data, `lowess` handles
irregular `x` and arbitrary curvature by fitting a local line in a moving
tricube-weighted neighbourhood:

```python
from quantforge import lowess

lowess(x, y, frac=0.3)                 # smoothed value at each x (span = 30% of points)
lowess(x, y, frac=0.5, iterations=3)   # 3 robustifying passes down-weight outliers
```

`frac` sets the span (larger = smoother); it reproduces a straight line exactly, cuts a
noisy sine's error well below the raw scatter, and keeps the input order on unsorted
data. The Cleveland robustifying `iterations` bisquare-down-weight points with large
residuals, so a handful of outliers don't distort the local fits — pin them back to the
trend rather than bending the curve toward them.

## Kalman filter (local level)

Track a slowly drifting level — a time-varying mean, a dynamic hedge ratio, a
smoothed signal — with the exact scalar Kalman recursion for the
random-walk-plus-noise model. `kalman_local_level` returns the filtered level, its
posterior variance, and the Kalman gain at each step;
`kalman_steady_state_gain` gives the closed-form limiting gain.

```python
from quantforge import kalman_local_level, kalman_steady_state_gain

y = [10.1, 10.3, 9.8, 10.5, 10.9, 11.2, 10.7, 11.5, 11.8, 12.0]
levels, variances, gains = kalman_local_level(y, process_var=0.05, obs_var=1.0)
levels[:3]                                  # -> [10.1, 10.172, 10.064]  (smoothed)

kalman_steady_state_gain(0.05, 1.0)         # -> 0.2  (limiting gain)
```

The gain depends only on the signal-to-noise ratio `Q/R`: it rises toward 1 as the
process noise dominates (trust each observation) and toward 0 as the observation
noise dominates (heavy smoothing). With `Q = 0` and a diffuse prior the estimate
is exactly the running mean (recursive least squares).

For a **time-varying regression slope** — a dynamic hedge ratio, factor loading,
or pairs beta — `kalman_regression_beta` filters `beta_t` from
`y_t = beta_t * x_t + v_t` with the slope following a random walk. A positive
process variance lets the beta drift and track a changing relationship:

```python
from quantforge import kalman_regression_beta

# y's slope on x jumps from 1.0 to 3.0 partway through the sample.
betas, variances = kalman_regression_beta(xs, ys, process_var=0.02, obs_var=0.09,
                                          beta0=1.0, p0=1.0)
# filtered slope: ~0.99 over the first regime, ~3.04 over the second
```

With `process_var = 0` and a diffuse prior it collapses to the static OLS slope
`sum(x*y) / sum(x^2)`.

For a general multivariate state-space model, `kalman_filter` runs the full
predict/update recursion with arbitrary transition `F`, observation `H`, and noise
covariances `Q`, `R`; `kalman_smoother` adds the RTS backward pass:

```python
from quantforge import kalman_filter, kalman_smoother

F = [[1, 1], [0, 1]]      # constant-velocity model: position, velocity
H = [[1, 0]]              # observe position only
Q = [[1e-4, 0], [0, 1e-4]]; R = [[0.5]]
res = kalman_filter(observations, F, H, Q, R, x0=[0.0, 0.0], P0=[[1, 0], [0, 1]])
res["states"], res["log_likelihood"]
sm = kalman_smoother(observations, F, H, Q, R, x0=[0.0, 0.0], P0=[[1, 0], [0, 1]])
```

Each observation is a vector; the filter returns filtered means, covariances and the
Gaussian data log-likelihood (for parameter tuning). The smoother conditions on the
whole series, so its covariances never exceed the filter's. With 1x1 matrices it
reduces exactly to `kalman_local_level`.

## Newey-West HAC variance

The sample variance understates the variance of a mean when observations are
autocorrelated. The Newey-West (1987) estimator corrects for it by adding
Bartlett-weighted autocovariances up to a truncation lag — the triangular weights
are exactly what keep the estimate non-negative. This is the long-run variance
behind HAC standard errors.

```python
from quantforge import (newey_west_variance, newey_west_mean_se,
                        autocorrelation)

x = [0.5, 0.7, 0.6, 0.9, 1.1, 0.8, 1.0, 1.3, 1.2, 1.4, 1.1, 1.5]

newey_west_variance(x, 0)      # -> 0.09243  (lag 0 == sample variance)
newey_west_variance(x, 3)      # -> 0.20573  (autocorrelation-corrected, larger)
newey_west_mean_se(x, 3)       # -> 0.13094  (HAC standard error of the mean)
autocorrelation(x, 1)          # -> 0.459
```

Lag 0 reduces to the sample variance; positive autocorrelation inflates both the
long-run variance and the mean's standard error. On an AR(1) the estimate climbs
toward the analytic `sigma^2 / (1 - phi)^2` as the lag grows.

The same HAC variance backs the Diebold-Mariano test of equal predictive accuracy,
which compares two forecasts by their loss differential:

```python
from quantforge import diebold_mariano

dm, p = diebold_mariano(errors1, errors2, h=1)   # negative dm favors forecast 1
```

`errors1`/`errors2` are the two forecast-error series; the statistic standardizes
their mean squared-error (or absolute-error) gap by its HAC standard error, with the
Harvey small-sample correction. It is antisymmetric — swapping the forecasts flips
the sign — so a large `|dm|` with a small p-value picks a winner.

Rather than pick one forecast, combine them — an average usually beats the best
single model. Three weighting schemes:

```python
from quantforge import (simple_average_forecast, inverse_mse_weights,
                        optimal_combination_weights, combine_forecasts)

simple_average_forecast([f1, f2, f3])          # equal weights
w = inverse_mse_weights([errors1, errors2])    # proportional to 1/MSE
w = optimal_combination_weights([errors1, errors2])   # Bates-Granger min-variance
combine_forecasts([f1, f2], w)                 # apply the weights to the forecasts
```

All schemes sum to one. `inverse_mse_weights` leans on the more accurate model;
`optimal_combination_weights` uses the error covariance, so it can put negative
weight on a model that hedges another's errors and drive the combined error variance
below either individual model's.

## Theil-Sen robust regression

Ordinary least squares is dragged by outliers. The Theil-Sen estimator takes the
median of the slopes of all point pairs, giving a ~29% breakdown point — up to
that fraction of the data can be arbitrarily corrupted before the fit blows up —
while staying exact on clean linear data.

```python
from quantforge import theil_sen

x = list(range(12))
y = [2.0 * xi + 5.0 for xi in x]
y[3] += 40      # two gross outliers
y[8] -= 35

theil_sen(x, y)         # -> (2.0, 5.0)   slope and intercept, unmoved by the outliers
```

The median pairwise slope shrugs off the two corrupted points that would tilt an
OLS line. Pairs sharing an `x` value are skipped.

`repeated_median_regression` (Siegel) pushes the robustness further to a **50%
breakdown point** — the theoretical maximum — by taking, for each point, the median of
its slopes to every other point, then the median of those:

```python
from quantforge import repeated_median_regression

x = list(range(10))
y = [2 * xi + 1 for xi in x]
y[5] = 1000                                  # gross outlier
repeated_median_regression(x, y)             # (2.0, 1.0) — completely unmoved
```

The double median tolerates nearly half the data being corrupted where Theil-Sen (a
single median of slopes, ~29% breakdown) starts to fail — the estimator to reach for
when contamination is severe, at the same `O(n²)` cost.

Between OLS and the fully-robust medians sits `huber_regression`, an M-estimator that
is quadratic (OLS-efficient) for small residuals and linear (bounded influence) beyond
a threshold `delta`, fit by iteratively reweighted least squares:

```python
from quantforge import huber_regression

# 50 points on y = 2x + 1 with four gross vertical outliers
h = huber_regression(X, y)
h["coefficients"][1]     # 1.997 — near the true slope 2 (OLS is dragged to ~1.65)
```

The default `delta = 1.345` gives ~95% efficiency at the normal while down-weighting
points whose robust-scaled residual exceeds it; a large `delta` recovers OLS exactly.
Unlike Theil-Sen/Siegel it extends naturally to multiple regressors and gives standard-
error-ready coefficients, trading a little breakdown resistance for efficiency and
generality.

When outliers exceed even 50% — the point where every median estimator fails —
`ransac_line` still recovers the true line by consensus: it fits many minimal (2-point)
candidates and keeps the one the most points agree with, then refits on that inlier set:

```python
from quantforge import ransac_line

# 100 points on y = 2x + 5, with 60 of them replaced by gross noise
r = ransac_line(x, y, threshold=3.0, n_iterations=500)
r["slope"], r["intercept"], r["n_inliers"]   # 2.001, 4.905, 40 (the clean points)
```

`threshold` is the residual within which a point counts as an inlier; a point is kept
if it lands there. RANSAC pays for its robustness with randomness (seed it for
reproducibility) and needs a threshold matched to the noise scale, but it is the tool
for majority-or-minority linear structure buried in heavy contamination — where
Theil-Sen, Siegel, and Huber all break down.

Theil-Sen still assumes `x` is exact. When *both* variables carry measurement error
(comparing two instruments or assays), OLS biases the slope toward zero;
`deming_regression` accounts for error in both and `orthogonal_regression` is its
symmetric (total-least-squares) special case:

```python
from quantforge import deming_regression, orthogonal_regression

deming_regression([1, 2, 3, 4, 5], [3, 5, 7, 9, 11])   # (2.0, 1.0) — exact on a line

x = [1, 2, 3, 4, 5, 6]
y = [1.1, 2.3, 2.9, 4.2, 5.1, 5.8]
orthogonal_regression(x, y)          # (0.9522, 0.234)  slope, intercept
deming_regression(x, y, lam=4.0)     # lam = var(err_x) / var(err_y)
```

The `lam` argument is the ratio of the two error variances: `lam = 1` (the
`orthogonal_regression` default) minimizes perpendicular distances and equals the
first principal-component slope; `lam -> infinity` recovers ordinary OLS of `y` on
`x`. The fitted line always passes through the sample means, and the orthogonal slope
is symmetric — fitting `x` on `y` gives exactly its reciprocal.

Deming assumes the errors are Gaussian; `passing_bablok_regression` is the
distribution-free method-comparison alternative — the standard choice in clinical
chemistry — combining errors-in-variables handling with full outlier resistance:

```python
from quantforge import passing_bablok_regression

passing_bablok_regression([1, 2, 3, 4, 5], [3, 5, 7, 9, 11])   # (2.0, 1.0)

x = list(range(1, 11))
y = [2 * xi + 1 for xi in x]
y[5] = 100                       # a gross outlier
passing_bablok_regression(x, y)  # (2.0, 1.0) — the median slope shrugs it off
```

It takes the median of all pairwise slopes but discards slopes of `-1` (and vertical
pairs) and *shifts* the median rank by the number of slopes below `-1`, the correction
that makes it consistent for method comparison rather than a plain Theil-Sen fit. It
needs no distributional assumption and tolerates outliers in either variable, at the
cost of the `O(n²)` pair enumeration.

Regression tells you how two methods *relate*; `bland_altman` and
`concordance_correlation` tell you how well they *agree*. A high correlation can hide
poor agreement — a constant offset correlates perfectly but never matches:

```python
from quantforge import bland_altman, concordance_correlation

r = bland_altman([10, 20, 30, 40], [12, 19, 33, 38])
r["bias"], r["lower"], r["upper"]    # -0.5, -5.1657, 4.1657 (95% limits of agreement)
r["means"], r["diffs"]               # per-point arrays for the Bland-Altman plot

concordance_correlation([1, 2, 3, 4, 5], [3, 4, 5, 6, 7])   # 0.5 (Pearson is 1.0)
```

`bland_altman` reports the bias (mean difference) and the limits `bias ± 1.96·sd`,
within which about 95% of method-to-method differences fall. Lin's
`concordance_correlation` folds precision (correlation) and accuracy (closeness to the
`y = x` line) into one `[-1, 1]` index that equals 1 only when every point lies on the
line of identity — so the constant `+2` offset above pulls it down to 0.5 even though
the Pearson correlation is a perfect 1.0.

With three or more raters, the intraclass correlation coefficient generalizes this to
"what fraction of the variance is between subjects rather than measurement noise".
`icc` returns the full Shrout-Fleiss family from a subjects × raters table:

```python
from quantforge import icc

data = [[9, 2, 5, 8], [6, 1, 3, 2], [8, 4, 6, 8],
        [7, 1, 2, 6], [10, 5, 6, 9], [6, 2, 4, 7]]

r = icc(data)
r["icc1"], r["icc2_1"], r["icc3_1"]   # 0.1657, 0.2898, 0.7148 (single-rating forms)
r["icc2_k"], r["icc3_k"]              # 0.6201, 0.9093 (mean-of-k-raters forms)
```

`icc1` is the one-way model (different raters per subject); `icc2_*` treats raters as a
random sample and measures absolute agreement; `icc3_*` treats them as fixed and
measures consistency (ignoring rater bias). The `_k` forms apply to the mean of `k`
raters and are related to the single-rating forms by the Spearman-Brown formula. On a
rater with a pure additive offset the consistency ICC(3) stays at 1 while the
agreement ICC(2) drops — the same distinction Bland-Altman's bias exposes.

For *categorical* ratings, kappa corrects the raw agreement for what chance alone would
produce. `cohen_kappa` handles two raters on nominal labels, `weighted_kappa` handles
ordinal categories (penalizing far-apart disagreements more), and `fleiss_kappa`
extends to any number of raters:

```python
from quantforge import cohen_kappa, weighted_kappa, fleiss_kappa

cohen_kappa(['y', 'y', 'n', 'n', 'y'], ['y', 'n', 'n', 'n', 'y'])   # 0.6154

# ordinal: an off-by-one disagreement is scored far above an off-by-two one
weighted_kappa([1, 2, 3, 2, 1], [1, 2, 2, 2, 1], "linear")   #  0.7059
weighted_kappa([1, 2, 3, 2, 1], [3, 2, 1, 2, 3], "linear")   # -0.3636

# m raters per subject: rows are per-category counts (each row sums to m)
fleiss_kappa([[0, 0, 0, 0, 14], [0, 2, 6, 4, 2], [0, 0, 3, 5, 6],
             [0, 3, 9, 2, 0], [2, 2, 8, 1, 1]])              # 0.2519
```

All three are `(p_observed - p_expected) / (1 - p_expected)`: 1 is perfect agreement,
0 is chance level, negative is worse than chance. On two categories a linear
`weighted_kappa` equals `cohen_kappa` exactly; `weighted_kappa` with `"quadratic"`
weights is the usual choice for ordered scales and links back to the ICC.

## Isotonic regression (monotone fit)

When you know the response only moves one way — a dose-response curve, a calibration
map from model scores to probabilities — but don't want to assume a functional form,
fit the best *monotone* step function instead. Isotonic regression finds the
non-decreasing sequence closest to the data in weighted least squares; the
pool-adjacent-violators algorithm solves it exactly in one linear sweep, pooling any
out-of-order adjacent blocks into their weighted mean.

```python
from quantforge import isotonic_regression

isotonic_regression([1, 2, 4, 2, 5])                     # -> [1, 2, 3, 3, 5]
isotonic_regression([5, 2, 4, 2, 1], increasing=False)   # -> [5, 3, 3, 2, 1]
```

The `4, 2` pair violates the increasing order, so it collapses to their mean `3`;
everything already in order is left untouched. Pass `weights` for a weighted fit and
`increasing=False` for a non-increasing one.

`isotonic_fit` sorts by `x`, fits, realigns the result to the original input order,
and returns an interpolating predictor — the shape used for probability calibration,
where raw model scores are mapped to monotone-increasing calibrated probabilities:

```python
from quantforge import isotonic_fit

scores = [0.1, 0.4, 0.35, 0.8, 0.7, 0.9]
labels = [0,   0,   1,    1,   1,   1]

y_hat, predict = isotonic_fit(scores, labels)
y_hat          # [0.0, 0.5, 0.5, 1.0, 1.0, 1.0]  fit aligned to input order
predict(0.5)   # 0.6667  linear interpolation between fitted knots
predict(0.0)   # 0.0     clamped to the low endpoint
```

The fitted probabilities never decrease as the score rises — the monotonicity a
calibration map should have — and `predict` interpolates linearly between knots,
clamping past the data range. It is the exact least-squares monotone fit, verified
against the independent max-min weighted-average formula.

## Robust scale and location

The standard deviation and the mean both break down under a single extreme point.
These estimators degrade gracefully: `median_absolute_deviation` (50% breakdown,
scaled to match the standard deviation under normality), `interquartile_range`,
`winsorize` (clip the tails), and `trimmed_mean`.

```python
from quantforge import (median_absolute_deviation, interquartile_range,
                        winsorize, trimmed_mean)

x = [10, 11, 9, 12, 10, 11, 9, 10, 13, 8, 1000]   # one gross outlier

median_absolute_deviation(x)          # -> 1.4826  (std would be ~299)
interquartile_range(x, scale=True)    # -> 1.4826  (normal-consistent scale)
trimmed_mean(x, 0.1)                  # -> 10.556  (mean would be ~99)
max(winsorize(x, 0.1))                # -> 13      (the 1000 is clipped in)
```

The scaled MAD and IQR both estimate the same underlying sigma the standard
deviation would give on clean data, but ignore the outlier that inflates it.

For robust *location* with far higher efficiency than the median, the
Hodges-Lehmann estimator is the median of all pairwise averages `(x_i + x_j) / 2` — a
29% breakdown point yet ~96% efficiency at the normal, and the point estimate that
pairs with the Wilcoxon signed-rank test:

```python
from quantforge import hodges_lehmann_location, hodges_lehmann_shift

x = [10, 11, 9, 12, 10, 11, 9, 10, 13, 8, 1000]   # one gross outlier
hodges_lehmann_location(x)      # 10.5  (mean would be ~100)

hodges_lehmann_shift([1, 2, 3, 4, 5], [6, 7, 8, 9, 10])   # 5  — the typical y - x gap
```

`hodges_lehmann_shift` is the two-sample analogue — the median of every pairwise
difference `y_j - x_i` — the robust shift estimate that inverts the Wilcoxon rank-sum
test. On a pure `+5` translation it returns exactly `5`, and it stays near zero when
the two samples share a distribution.

The MAD needs a symmetric distribution and a location estimate; the Rousseeuw-Croux `Qn`
and `Sn` estimators need neither and share the MAD's 50% breakdown:

```python
from quantforge import qn_scale, sn_scale, biweight_midvariance

sn_scale(list(range(1, 11)))          # 3.5778 — location-free spread
qn_scale([1, 2, 1.5, 2.5, 1.8, 1000]) # ~1.56 — the gross outlier barely moves it
biweight_midvariance(data)             # smoothly downweighted robust sigma
```

`qn_scale` is a scaled quantile of the pairwise distances and `sn_scale` a median of
per-point medians — both estimate the normal sigma on clean data yet tolerate up to half
the sample being outliers, and neither assumes symmetry (unlike the MAD). Each carries its
normal-consistency constant. `biweight_midvariance` smoothly zero-weights points beyond a
few MADs for ~87% efficiency at the normal. (`qn_scale` omits the finite-sample correction
some libraries apply, so small-sample values differ slightly.)

To see a sample's *shape* without a histogram's arbitrary bin edges, `kde` builds a
smooth Gaussian kernel density estimate:

```python
from quantforge import kde, kde_function, silverman_bandwidth

kde(returns, x=0.0)                     # density at a point (bandwidth auto-selected)
f = kde_function(returns)               # a callable density estimator
kde(returns, [-0.02, 0.0, 0.02])        # evaluate at several points at once
silverman_bandwidth(returns)            # the rule-of-thumb bandwidth it uses
```

The estimate is non-negative and integrates to one, converges to the true density as
the sample grows (recovering the standard-normal peak of ~0.399 on N(0,1) data), and
resolves multiple modes a single summary statistic would hide. The bandwidth defaults
to Silverman's rule (`0.9 min(std, IQR/1.34) n^{-1/5}`, robust to mild non-normality);
pass `rule="scott"` or an explicit `bandwidth` to override — smaller is spikier, larger
is smoother.

For the raw distributional building blocks, `ecdf` gives the empirical CDF, `quantile`
the sample quantile (the inverse), and `qq_points` a quantile-quantile pairing:

```python
from quantforge import ecdf, quantile, qq_points

ecdf(sample, x=0.0)                    # fraction of the sample at or below 0
quantile(sample, 0.95)                 # the 95th percentile (linear interpolation)
qq_points(sample, reference)           # [(ref_quantile, sample_quantile), ...] for a Q-Q plot
```

`ecdf` is a right-continuous step from 0 to 1 (scalar or vectorized over `x`);
`quantile` supports the usual `linear`/`lower`/`higher`/`nearest` conventions and
matches `statistics.quantiles` on the linear rule. `qq_points` pairs the two samples at
matched plotting positions — points fall on the `y = x` line when the distributions
agree, and a slope or bend exposes a scale or shape difference (the visual companion to
the KS and Cramér-von Mises tests).

## Hurst exponent (long memory)

Measure the persistence of a series with the Hurst exponent via rescaled-range
(R/S) analysis: `H ~ 0.5` is memoryless (a random walk's increments), `H > 0.5`
is persistent/trending, and `H < 0.5` is anti-persistent/mean-reverting.

```python
from quantforge import hurst_exponent, rescaled_range

hurst_exponent(white_noise)      # ~0.5  (no memory)
hurst_exponent(price_levels)     # ~1.0  (a random walk in levels)
hurst_exponent(mean_reverting)   # <0.5  (moves tend to reverse)
```

`H` is the slope of `log(R/S)` against `log(window)` over dyadic window sizes.
A reading well above 0.5 flags a trend-following regime; well below flags a
mean-reversion regime.

For a long-memory series, integer differencing (`d = 1`) removes the unit root but
destroys the memory; fractional differencing `(1 - L)^d` for real `d` makes it
stationary while keeping most of the autocorrelation:

```python
from quantforge import (fractional_difference, fixed_width_fracdiff,
                        fracdiff_weights)

fractional_difference(series, d=0.4)         # (1-L)^0.4 by the full expansion
fixed_width_fracdiff(series, d=0.4)          # Lopez de Prado's fixed-width window
fracdiff_weights(0.4, 10)                    # the binomial filter weights
```

Integer orders recover the ordinary differences (`d = 1` is the first difference,
`d = 2` the second); a fractional `d` in `(0, 1)` is the useful middle ground, with
weights that decay slowly rather than truncating at order `d`.

To pick `d` from the data, `gph_estimate` runs the Geweke-Porter-Hudak
log-periodogram regression; `fractional_integrate` is the inverse operator
`(1 - L)^{-d}`, handy for generating test series:

```python
from quantforge import gph_estimate, fractional_integrate

gph_estimate(series)["d"]                    # estimated memory parameter (+ std_error)
fractional_integrate(white_noise, d=0.4)     # build an ARFIMA(0, 0.4, 0) path
```

`gph_estimate` returns `d` and its asymptotic standard error, so you can test
`d = 0` (short memory); estimate it, then difference by that `d` with
`fractional_difference` to obtain a stationary, memory-preserved series.

`dfa_exponent` is a third estimator, detrended fluctuation analysis, which removes a
polynomial trend in each window and so tolerates slow drift that biases R/S:

```python
from quantforge import dfa_exponent, dfa_fluctuations

dfa_exponent(series)                 # alpha ~ 0.5 white, > 0.5 persistent, < 0.5 anti
dfa_fluctuations(series)             # (scales, F(s)) for the log-log fit
```

`alpha` equals the Hurst exponent for a stationary long-memory series and exceeds it
by one for the integrated version, so a random walk reads `~1.5` where white noise
reads `~0.5`.

## Entropy (time-series regularity)

Measure how unpredictable a series is. `approximate_entropy` and `sample_entropy`
score the log-likelihood that short runs which stay close remain close one step
later; `permutation_entropy` is the entropy of the ordinal patterns, normalized to
`[0, 1]`:

```python
from quantforge import approximate_entropy, sample_entropy, permutation_entropy

sample_entropy(returns)              # 0 = perfectly regular, larger = more complex
approximate_entropy(returns, m=2)    # Pincus ApEn (biased low on short series)
permutation_entropy(returns, m=3)    # 0 monotone, ~1 for random orderings
```

A clean sine scores near zero; white noise scores high. `sample_entropy` drops the
self-matches that bias `approximate_entropy`, and `permutation_entropy` is invariant
to any monotone transform of the data, so it reads the same on prices or their logs.

For dependence *between* two series, `mutual_information` scores shared information
(symmetric, zero iff independent) and `transfer_entropy` scores *directed* flow —
how much one series's past predicts another's future:

```python
from quantforge import mutual_information, transfer_entropy

mutual_information(x, y)             # nats of shared information (symmetric)
transfer_entropy(source, target)    # directed: source's past -> target's future
```

Transfer entropy is asymmetric, so `transfer_entropy(x, y)` and
`transfer_entropy(y, x)` reveal lead-lag: in a coupled system where `x` drives `y`,
the `x -> y` value dominates. Both are near zero for independent series.

To compare two *distributions* directly (histograms, empirical PMFs, model vs market
probabilities), the divergence functions quantify how far apart they are:

```python
from quantforge import (kl_divergence, jensen_shannon_divergence,
                        hellinger_distance, total_variation_distance)

kl_divergence([0.7, 0.3], [0.5, 0.5])              # 0.08228 nats (asymmetric)
jensen_shannon_divergence([0.7, 0.3], [0.5, 0.5])  # 0.02101 (symmetric, <= log 2)
total_variation_distance([0.7, 0.3], [0.5, 0.5])   # 0.2 (largest event-probability gap)
```

Inputs are any non-negative weight vectors (normalized internally). `kl_divergence` is
the asymmetric information loss (infinite where the reference has zero mass);
`jensen_shannon_divergence` is its symmetric, always-finite cousin;
`hellinger_distance` and `total_variation_distance` are bounded `[0, 1]` metrics; and
`bhattacharyya_distance` measures overlap. They satisfy the Pinsker inequality
(`TV <= sqrt(KL/2)`) and all vanish exactly when the two distributions coincide.

Those compare aligned *bins*; to compare two raw *samples* (accounting for how far mass
must move along the axis), `wasserstein_distance` gives the earth-mover distance:

```python
from quantforge import wasserstein_distance

wasserstein_distance([1, 2, 3, 4, 5], [4, 5, 6, 7, 8])   # 3.0 — a pure shift of 3
wasserstein_distance([0.0], [1.0])                       # 1.0 — mass moved one unit
wasserstein_distance(sample_a, sample_b, p=2)            # quadratic transport cost
```

Unlike the divergences it works on samples of *different* sizes (merging their empirical
CDFs), is finite even for disjoint supports, and — being a ground-distance metric —
reflects *how far* probability mass moved, not just that it differs. A constant shift of
the whole sample returns exactly that shift; `wasserstein1_sorted` is a faster path when
the two samples are the same length.

## Variance-ratio test

The Lo-MacKinlay variance ratio tests the random-walk null: under it the variance
of a `q`-period return is `q` times the one-period variance, so `VR(q) = 1`.
`VR > 1` signals momentum (positive serial correlation); `VR < 1` signals mean
reversion. `variance_ratio_zstat` returns a heteroskedasticity-robust `z` (valid
under GARCH-type noise) — `|z| > 1.96` rejects at 5%.

```python
from quantforge import variance_ratio, variance_ratio_zstat

variance_ratio(returns, q=2)         # ~1 random walk, >1 trending, <1 mean-reverting
variance_ratio_zstat(returns, q=2)   # standard-normal z; |z| > 1.96 rejects the walk
```

On white noise `VR` sits near 1 with a small `z`; a mean-reverting series drives
`VR` below 1 with a large negative `z`, a trending one above 1 with a large
positive `z`.

## Serial correlation (Ljung-Box / Durbin-Watson)

Test whether a series — usually model residuals — is white noise. `ljung_box`
(and `box_pierce`) check that the first `h` autocorrelations are jointly zero;
`durbin_watson` checks first-order autocorrelation.

```python
from quantforge import ljung_box, box_pierce, durbin_watson

q, p = ljung_box(residuals, lags=10)    # small p rejects "no autocorrelation"
durbin_watson(residuals)                # ~2 = clean, <2 positive, >2 negative
```

White noise gives a large Ljung-Box p-value and a Durbin-Watson near 2; an AR(1)
gives a huge Q with a vanishing p-value and a Durbin-Watson far from 2. The
p-values use a self-contained chi-square survival function (no SciPy).

The Wald-Wolfowitz runs test checks randomness a different way — by counting streaks
rather than autocorrelation:

```python
from quantforge import runs_test, runs_test_binary

runs_test(values)                # dichotomize about the median, test the sign runs
runs_test_binary([1, 0, 1, 1, 0, 0])   # two-symbol sequence directly
```

Both return a two-sided `(z, p_value)`: `z < 0` (too few runs) flags clustering or a
trend, `z > 0` (too many) flags over-alternation / mean reversion. It is
distribution-free, so it catches nonrandomness a linear autocorrelation test can
miss.

## Benford's law (first-digit anomaly detection)

Numbers spanning several orders of magnitude have leading digits distributed as
`P(d) = log10(1 + 1/d)` — digit 1 leads ~30% of the time. Fabricated or constrained
data often breaks this, so a conformance test flags datasets worth auditing:

```python
from quantforge import (benford_chi_square, benford_mad, first_digit_distribution,
                        benford_expected)

benford_chi_square(values)     # (chi2, p): small p rejects Benford conformance
benford_mad(values)            # Nigrini MAD: < 0.006 close, > 0.015 nonconforming
first_digit_distribution(values)   # observed (counts, proportions) for digits 1..9
```

Fibonacci numbers and powers of two conform (they span many magnitudes); a uniform
sample is strongly rejected. The chi-square test is sample-size sensitive, so pair
it with the size-independent `benford_mad` on large datasets.

## Goodness of fit (Jarque-Bera / KS)

Test distributional assumptions. `jarque_bera_test` checks normality from skew and
kurtosis; `ks_two_sample` checks whether two samples share a distribution via the
largest gap between their empirical CDFs.

```python
from quantforge import jarque_bera_test, ks_two_sample

stat, p = jarque_bera_test(returns)      # small p rejects normality
d, p = ks_two_sample(sample_a, sample_b) # small p rejects "same distribution"
```

A normal sample is not rejected; a heavy-tailed or skewed one is. The KS test
catches any distributional difference — location, scale, or shape — not just a
difference in means, and its p-value uses the asymptotic Kolmogorov distribution
(no SciPy).

For a normality check that weights the tails, `anderson_darling_normal` is more
powerful than KS against the heavy-tail and skew departures typical of returns:

```python
from quantforge import anderson_darling_normal

a2, p = anderson_darling_normal(returns)   # small p rejects normality
```

It estimates the mean and standard deviation from the sample, applies the Stephens
small-sample adjustment, and returns a D'Agostino-Stephens p-value. A clean normal
sample passes; an exponential or fat-tailed one is strongly rejected.

`dagostino_k2` gives an omnibus normality test that reports *why* it rejects — a
skewness Z and a kurtosis Z combined into one chi-square statistic:

```python
from quantforge import dagostino_k2

r = dagostino_k2(returns)
r["k2"], r["p_value"]        # omnibus statistic and p-value
r["z_skew"], r["z_kurt"]     # which shape departure drives the rejection
```

It transforms the sample skewness and kurtosis to standard-normal scores (D'Agostino
and Anscombe-Glynn) and sums their squares, so a large `z_skew` flags asymmetry and a
large positive `z_kurt` flags heavy tails (negative = light tails). It passes clean
normals, has a well-calibrated ~5% null rejection rate, and detects the skew and fat
tails typical of returns — the component Z scores tell you which. Needs a moderate
sample (>= 20).

To test whether *several* samples share one (unspecified) distribution — the
nonparametric analogue of one-way ANOVA, but sensitive to any distributional difference
— `anderson_darling_ksample` runs the Scholz-Stephens k-sample AD test:

```python
from quantforge import anderson_darling_ksample

anderson_darling_ksample(group_a, group_b, group_c)
# {'a2k': 8.3926, 'standardized': 4.4798, 'p_value': 0.0055} on the SS worked example
```

It compares each sample's empirical CDF to the pooled CDF with the tail-weighting
Anderson-Darling metric, handles ties, and reports the raw `a2k`, the `standardized`
statistic `(A2k - (k-1))/sqrt(var)`, and an interpolated `p_value`. It reproduces the
Scholz-Stephens (1987) worked example, stays non-significant when the samples share a
distribution, and flags differences in location, scale, or shape that a means-only test
would miss.

For a fully distribution-free two-sample test that needs no CDF assumption at all,
`energy_distance` and its permutation test `energy_test` compare two samples by the
Székely-Rizzo statistic `2A - B - C` (mean cross-sample distance minus the two mean
within-sample distances):

```python
from quantforge import energy_distance, energy_test

energy_distance([1, 2, 3, 4, 5], [1, 2, 3, 4, 5])   # 0.0 — identical samples

res = energy_test(sample_a, sample_b, n_permutations=299)
res["statistic"], res["p_value"]        # small p rejects "same distribution"
```

The statistic is zero if and only if the empirical distributions coincide and grows
with any distributional gap — location, scale, or shape. `energy_test` pools the two
samples, reshuffles the group labels `n_permutations` times, and reports the fraction
of permutations reaching the observed distance, so it catches differences a
mean-focused t-test can miss (for example two samples with the same mean but different
spread, where the energy test's p-value drops well below 0.05).

Where KS uses only the single largest gap between the two empirical CDFs,
`cramer_von_mises_2samp` integrates the *squared* gap across the whole curve — more
powerful against differences spread through the distribution — and returns an
asymptotic p-value from the limiting Cramér-von Mises law (no permutation needed):

```python
from quantforge import cramer_von_mises_2samp

cramer_von_mises_2samp([1, 3, 5, 7, 9], [2, 4, 6, 8, 10])   # interleaved
# {'statistic': 0.05, 'p_value': 0.8763}  — same distribution, not rejected

cramer_von_mises_2samp([1, 2, 3, 4, 5], [6, 7, 8, 9, 10])   # fully separated
# {'statistic': 0.85, 'p_value': 0.0056}  — clearly different
```

It uses Anderson's (1962) rank form of the statistic and evaluates the limiting
distribution's upper tail through its `K_{1/4}` Bessel-function series; the asymptotic
critical values it produces (0.461 → 0.05, 0.743 → 0.01) match the published tables,
and its p-value tracks a permutation test closely.

The classical parametric tests are built on the distribution CDFs, each returning
`(statistic, p_value)`:

```python
from quantforge import (chi_square_gof_test, chi_square_independence_test,
                        one_way_anova, one_sample_t_test, paired_t_test,
                        two_sample_t_test, mann_whitney_u, binomial_test)

chi_square_gof_test([16, 18, 16, 14, 12, 12])            # fit vs uniform
chi_square_independence_test([[10, 20], [30, 40]])        # contingency table
one_way_anova([5.1, 4.9, 5.5], [6.1, 5.9, 6.3])           # equal group means?
one_sample_t_test(sample, mu0=5.0)                        # mean vs a reference
paired_t_test(before, after)                              # dependent samples
two_sample_t_test(sample_a, sample_b, equal_var=False)    # Welch two-sample t
mann_whitney_u(sample_a, sample_b)                        # distribution-free
binomial_test(k=8, n=10, prob=0.5)                        # exact binomial
```

`one_way_anova` on two groups reproduces the pooled `two_sample_t_test` exactly
(`F = t^2`), `paired_t_test` equals a one-sample t on the within-pair differences,
`mann_whitney_u` is the distribution-free alternative when normality is doubtful,
the binomial test is exact (no normal approximation), and the discrete p-values
come from the same gamma/beta identities behind the distribution CDFs.

`chi_square_independence_test` says *whether* two categorical variables are associated;
the association measures say *how strongly*, on a size-independent `[0, 1]` scale:

```python
from quantforge import cramers_v, phi_coefficient, contingency_coefficient

cramers_v([[50, 0], [0, 50]])          # 1.0 — perfect association
phi_coefficient([[10, 20], [20, 40]])  # 0.0 — proportional rows are independent
cramers_v([[20,5,5],[5,20,5],[5,5,20]])# 0.5 — moderate association on a 3x3 table
```

`cramers_v` is the general r x c effect size (`sqrt(chi2 / (n * min(r-1, c-1)))`), `0`
under independence and `1` for a perfect association; `phi_coefficient` is its 2x2 special
case (equal to Cramer's V there and to the closed-form `(ad-bc)` phi). `tschuprow_t` is an
alternative normalization equal to V on square tables, and `contingency_coefficient`
(`sqrt(chi2/(chi2+n))`) is Pearson's, which never quite reaches 1. All are built from the
same Pearson chi-square as the independence test.

ANOVA assumes equal group variances; `levene_test` and `bartlett_test` check that:

```python
from quantforge import levene_test, bartlett_test

levene_test(group_a, group_b, group_c)                 # median-centered (Brown-Forsythe)
levene_test(group_a, group_b, group_c, center="mean")  # original Levene
bartlett_test(group_a, group_b, group_c)               # likelihood ratio (chi-square)
```

`levene_test` runs one-way ANOVA on each point's absolute deviation from its group
center — median-centered by default, the robust Brown-Forsythe form. `bartlett_test` is
a likelihood-ratio test, more powerful under normality but sensitive to fat tails. Both
stay non-significant on equal-variance groups and reject sharply when the spreads differ
(a 1×/3×/6× split gives `p ~ 0`); pair them with `ansari_bradley_test` / `mood_test`
(the two-sample rank versions) for the fully nonparametric route.

For the *paired* nonparametric case, `wilcoxon_signed_rank_test` is the
distribution-free counterpart of the paired t-test (and the test the Hodges-Lehmann
location estimator inverts), while `sign_test` uses only the signs and so assumes
nothing about symmetry:

```python
from quantforge import wilcoxon_signed_rank_test, sign_test

before = [125, 132, 128, 140, 135]
after  = [120, 128, 125, 138, 130]

r = wilcoxon_signed_rank_test(before, y=after)   # paired differences
r["statistic"], r["p_value"]                     # W+ = 15.0, p = 0.0579

sign_test([1, 1, 1, 1, 1, 1, 1, 1, -1, -1])      # n_plus 8 of 10
# {'n_plus': 8, 'n': 10, 'p_value': 0.10937, ...}  exact two-sided binomial
```

The signed-rank test ranks the absolute differences (average ranks for ties, zero
differences dropped) and sums the positive ranks, comparing to a continuity-corrected
normal approximation; on the classic all-positive textbook sample it gives `W+ = 45`
and `p = 0.0092`. The sign test counts observations above the hypothesized median and
tests that count against `Binomial(n, 1/2)`, returning the exact two-sided p-value —
maximally robust, at the cost of power.

All the rank tests above chase a *location* shift and are blind to a pure difference
in *spread*. `ansari_bradley_test` and `mood_test` fill that gap — nonparametric
two-sample tests of equal dispersion:

```python
from quantforge import ansari_bradley_test, mood_test, mann_whitney_u

# two samples, same median, very different spread (sigma 1 vs 4)
ansari_bradley_test(x, y)["p_value"]   # small — flags the scale difference
mood_test(x, y)["p_value"]             # small — same conclusion

mann_whitney_u(x, y)                   # a location test barely reacts (p ~ 0.1)
```

Ansari-Bradley scores the pooled ranks from the outside in (`min(r, N+1-r)`), so a
tightly-concentrated sample earns higher scores; Mood sums each observation's squared
deviation from the center rank, which a dispersed sample inflates by pushing values to
the extremes. Both assume the two samples share a location and return a
normal-approximation p-value whose null z is calibrated to mean 0 and variance 1.

For *more than two* groups, `kruskal_wallis_test` is the rank analogue of one-way
ANOVA, and `friedman_test` handles repeated measures (a nonparametric two-way ANOVA
with one observation per cell):

```python
from quantforge import kruskal_wallis_test, friedman_test

kruskal_wallis_test([2.9, 3.0, 2.5, 2.6, 3.2],
                    [3.8, 2.7, 4.0, 2.4],
                    [2.8, 3.4, 3.7, 2.2, 2.0])
# {'statistic': 0.7714, 'df': 2, 'p_value': 0.68}  — no group difference

# one row per block (subject), one column per treatment
friedman_test([[1, 2, 3], [2, 3, 1], [3, 1, 2], [1, 2, 3], [2, 3, 1]])
# {'statistic': 0.4, 'df': 2, 'p_value': 0.8187}
```

Kruskal-Wallis pools and ranks all observations, comparing each group's mean rank to
the overall mean (tie-corrected, chi-square with `k - 1` df); on two groups its `H`
equals the Mann-Whitney `z²` exactly. Friedman ranks *within* each block across the
treatments, so it removes between-block variation the way a paired test does — use it
when the same subjects are measured under every condition.

A significant Kruskal-Wallis tells you *some* group differs but not which pair.
`dunn_test` is the post-hoc follow-up: it compares every pair using the one pooled
ranking (so it stays consistent with the H statistic) and adjusts the p-values for the
multiple comparisons:

```python
from quantforge import dunn_test

for r in dunn_test([[1, 2, 3, 4, 5], [6, 7, 8, 9, 10], [11, 12, 13, 14, 15]]):
    r["groups"], r["z"], r["p_value"], r["p_adjusted"]
# (0, 1)  z -1.768  p 0.0771  adj 0.1542
# (0, 2)  z -3.536  p 0.0004  adj 0.0012   <- only the far-apart pair survives
# (1, 2)  z -1.768  p 0.0771  adj 0.1542
```

Each pair's z uses the tie-corrected standard error from the pooled ranking, and
`adjust` selects `"holm"` (default), `"bonferroni"` or `None`. The pairwise `z²`
equals the Kruskal-Wallis `H` on two groups, so Dunn's test is the natural drill-down
after the omnibus test rejects.

When the groups have a *natural order* (dose levels, time buckets, rating tiers) and
you expect a monotone response, `jonckheere_terpstra_test` is a directional trend test
— much more powerful than Kruskal-Wallis, which spends power on differences in every
direction:

```python
from quantforge import jonckheere_terpstra_test

jonckheere_terpstra_test([[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12]])
# {'statistic': 48.0, 'mean': 24.0, 'variance': 46.6667, 'z': 3.5132, 'p_value': 0.0004}
```

Pass the groups in the hypothesized order. It sums the Mann-Whitney concordances over
every ordered pair, so a positive `z` signals an increasing trend across the group
order and a negative one a decreasing trend (a perfectly reversed arrangement gives
the mirror-image z). The null mean and tie-corrected variance are the closed-form
Jonckheere values, and the p-value comes from the normal approximation.

For paired *binary* outcomes, `mcnemar_test` is the analogue of the paired t-test (two
measurements per subject) and `cochran_q_test` extends it to `k` binary treatments —
the binary special case of Friedman:

```python
from quantforge import mcnemar_test, cochran_q_test

mcnemar_test(b=10, c=2)      # b, c are the discordant counts
# {'b': 10, 'c': 2, 'p_value': 0.0386, 'chi2_cc': 4.0833, 'p_value_chi2': 0.0433}

# k binary treatments over shared blocks (one row per block)
cochran_q_test([[0, 0, 1], [0, 1, 1], [0, 0, 1], [1, 0, 1], [0, 0, 1], [0, 1, 1]])
# {'statistic': 7.0, 'df': 2, 'p_value': 0.0302}
```

McNemar looks only at the *discordant* pairs (`b` = 0→1, `c` = 1→0) and tests
`b == c`, returning the exact two-sided binomial p-value plus the continuity-corrected
chi-square; pass a 2×2 `table=[[a, b], [c, d]]` instead of `b`/`c` if you have the full
table. Cochran's Q reduces to the uncorrected McNemar statistic on two treatments, and
its statistic is chi-square with `k - 1` degrees of freedom.

For an *unpaired* 2×2 table with small counts — where the chi-square approximation is
unreliable — `fisher_exact_test` gives the exact p-value from the hypergeometric
distribution, plus the sample odds ratio:

```python
from quantforge import fisher_exact_test

fisher_exact_test([[3, 1], [1, 3]])                 # the classic tea-tasting table
# {'odds_ratio': 9.0, 'p_value': 0.4857}

fisher_exact_test([[3, 1], [1, 3]], "greater")["p_value"]   # 0.2429  one-sided
fisher_exact_test([[8, 2], [1, 5]])                 # {'odds_ratio': 20.0, 'p_value': 0.035}
```

Conditioning on the margins, it sums the probabilities of every table no more likely
than the observed one (two-sided), or the upper/lower tail for `"greater"` / `"less"`.
Use it in place of `chi_square_independence_test` whenever any expected cell count is
small (the usual rule of thumb is below 5).

For a proportion estimate, four confidence intervals span the accuracy/simplicity
trade-off:

```python
from quantforge import (wald_interval, wilson_interval, agresti_coull_interval,
                        clopper_pearson_interval)

wald_interval(k=8, n=25)              # normal approximation (simplest)
wilson_interval(8, 25)                # score interval, good small-sample coverage
agresti_coull_interval(8, 25)
clopper_pearson_interval(8, 25)       # exact, never under-covers
```

`clopper_pearson_interval` inverts the binomial CDF and is the widest (conservative
but exact); `wilson_interval` is the practical default; `wald_interval` is accurate
only for large `n` away from 0 or 1. All stay within `[0, 1]`.

Before running a test, size it. The power routines give the rejection probability
for a planned design and invert it for the required sample size:

```python
from quantforge import (two_sample_t_power, two_sample_t_sample_size,
                        proportion_sample_size)

two_sample_t_power(effect_size=0.6, n_per_group=40)      # power at this design
two_sample_t_sample_size(effect_size=0.5, power=0.80)    # -> ~64 per group
proportion_sample_size(p1=0.50, p2=0.65, power=0.80)     # -> ~170 per group
```

`effect_size` is Cohen's `d` (mean gap in pooled-SD units). The computed power
matches an empirical Monte-Carlo rejection rate; a zero effect gives power equal to
the significance level, and the sample sizes reproduce the standard textbook
values. One-sample analogues are `one_sample_z_power` / `one_sample_z_sample_size`.

Running many tests inflates false positives; the correction routines adjust a set
of raw p-values so a fixed threshold controls the family-wise error or the
false-discovery rate:

```python
from quantforge import (bonferroni, holm, benjamini_hochberg,
                        benjamini_yekutieli)

pvals = [0.01, 0.02, 0.03, 0.04, 0.05]
bonferroni(pvals)            # [0.05, 0.10, 0.15, 0.20, 0.25]  (family-wise)
holm(pvals)                  # step-down, dominates Bonferroni
benjamini_hochberg(pvals)    # [0.05, ...]  (false-discovery rate)
benjamini_yekutieli(pvals)   # FDR under arbitrary dependence
```

Each returns monotone adjusted p-values aligned with the input. `holm` is a strict
improvement on `bonferroni` for family-wise control; `benjamini_hochberg` is the
usual choice when many discoveries are expected and a controlled FDR is acceptable.

## Autocorrelation (ACF / PACF)

Identify ARMA structure from the correlograms. `acf` is the autocorrelation at each
lag; `pacf` is the partial autocorrelation (via the Durbin-Levinson recursion),
which removes the intervening lags.

```python
from quantforge import acf, pacf

acf(series, nlags=20)     # [1.0, r1, r2, ...] — decays geometrically for an AR
pacf(series, nlags=20)    # [1.0, p1, p2, ...] — cuts off after lag p for an AR(p)
```

An AR(p) shows a PACF that drops to ~0 beyond lag `p` with a geometrically
decaying ACF; an MA(q) shows the mirror image — an ACF that cuts off after lag `q`.
White noise is flat in both.

Once the order is chosen, `fit_ar_yule_walker` fits the AR(p) coefficients by the
Yule-Walker equations and `ar_forecast` projects the series forward:

```python
from quantforge import fit_ar_yule_walker, ar_forecast

m = fit_ar_yule_walker(series, order=1)
m["coefficients"]      # [phi_1, ...]
m["noise_variance"]    # innovation variance
ar_forecast(m, series[-1:], steps=5)   # mean-reverts toward m["mean"]
```

The fit returns the coefficients, the mean-derived intercept, and the innovation
variance; the forecast iterates the deterministic recursion, decaying toward the
long-run mean for a stationary process.

To pick the order automatically, `select_ar_order` minimizes AIC or BIC over a
range of candidates:

```python
from quantforge import select_ar_order

order, scores = select_ar_order(series, max_order=8, criterion="bic")
# order = the BIC-minimizing AR order; scores = [(p, aic, bic), ...]
```

BIC's heavier `ln(n)` penalty favors more parsimonious models than AIC, so the
BIC-selected order is never larger — it recovers the true order on simulated
AR(1)/AR(2) data where AIC tends to over-fit.

## Exponential smoothing (Holt / Holt-Winters)

Forecast a series as a decaying blend of its recent level, trend, and season.
`holt_linear` is double smoothing (level + trend); `holt_winters_add` adds an
additive seasonal component of a given period.

```python
from quantforge import holt_linear, holt_winters_add

level, trend, fc = holt_linear(series, alpha=0.6, beta=0.4, horizon=3)
# fc = straight-line forecast l + h*b

level, trend, seasonals, fc = holt_winters_add(
    series, alpha=0.3, beta=0.1, gamma=0.3, period=4, horizon=4)
# fc repeats the learned seasonal pattern around the trend
```

Holt recovers a pure linear trend and forecasts a straight line; Holt-Winters
reproduces a repeating seasonal pattern and tracks any underlying trend. Both are
causal O(n) recursions with no external dependencies.

## Forecast accuracy

Score forecasts against realized values with scale-dependent and scale-free
errors. `mase` divides the MAE by a seasonal-naive benchmark, so `< 1` beats naive.

```python
from quantforge import mae, rmse, mape, smape, mase, theil_u1, theil_u2

actual, fc = [100, 102, 101, 105, 108], [101, 100, 103, 104, 107]
mae(actual, fc)      # 1.4    absolute error
rmse(actual, fc)     # 1.483  penalizes big misses
mape(actual, fc)     # 0.0136 scale-free (fraction)
smape(actual, fc)    # 0.0136 symmetric, bounded [0, 2]
mase(actual, fc, train=[95, 96, 98, 99, 100], season=1)   # vs naive: <1 beats it
theil_u2(actual, fc)  # forecast RMSE / no-change RMSE: <1 beats the random walk
theil_u1(actual, fc)  # inequality coefficient in [0, 1]
```

MAPE is scale-invariant but undefined at zero actuals; sMAPE is bounded and robust
to small values; MASE is the scale-free choice for cross-series comparison.
`theil_u2` benchmarks against the persistence (no-change) forecast — `< 1` means the
model adds value over a random walk — while `theil_u1` is a bounded, symmetric
inequality measure.

Those score a point forecast; for forecasts that carry uncertainty, proper scoring
rules grade the whole distribution:

```python
from quantforge import (pinball_loss, interval_score, interval_coverage,
                        crps_ensemble)

pinball_loss(actual, quantile_fc, tau=0.9)         # quantile-regression loss
interval_score(actual, lower, upper, alpha=0.1)    # Winkler score for a 90% interval
interval_coverage(actual, lower, upper)            # fraction inside (should be ~0.9)
crps_ensemble(y, ensemble_members)                 # CRPS of a sample forecast
```

`pinball_loss` is minimized at the true quantile (the objective behind quantile
regression); `crps_ensemble` generalizes absolute error to a full predictive
distribution and is zero only for a perfect forecast; `interval_score` rewards
narrow intervals but penalizes actuals that fall outside, and `interval_coverage`
checks the nominal level is met.

## Rank dependence (Kendall / Spearman)

Rank-based dependence captures monotone (not just linear) co-movement and is
invariant to any monotone transform of the margins — the natural language of
copulas. `kendall_tau` and `spearman_rho` measure it; `pseudo_observations` maps a
margin to its scaled ranks for empirical-copula work.

```python
from quantforge import kendall_tau, spearman_rho, pseudo_observations

kendall_tau(x, y)              # concordant - discordant, in [-1, 1]
spearman_rho(x, y)             # Pearson correlation of ranks, in [-1, 1]
pseudo_observations([12, 5, 9, 20])   # -> [0.6, 0.2, 0.4, 0.8]  (ranks / (n+1))
```

Both equal +1 for a strictly increasing relationship and -1 for decreasing;
Spearman is unchanged by exponentiating a margin, where Pearson would move. Feed
the pseudo-observations of each margin into a copula fit.

Plain `kendall_tau` is the tau-a variant: ties drag it below `+/-1` even under a
perfect monotone relation. `kendall_tau_b` applies the standard tie correction, and
`kendall_tau_test` adds a significance test of `tau = 0`:

```python
from quantforge import kendall_tau, kendall_tau_b, goodman_kruskal_gamma, kendall_tau_test

x = [1, 2, 2, 3, 4]
y = [10, 20, 20, 30, 40]        # perfectly monotone, with an aligned tie

kendall_tau(x, y)               # 0.9   (tau-a undershoots because of the tie)
kendall_tau_b(x, y)             # 1.0   (tie-corrected — reaches +1)
goodman_kruskal_gamma(x, y)     # ignores tied pairs entirely, (C - D) / (C + D)

res = kendall_tau_test(x, y)
res["tau_b"], res["z"], res["p_value"]   # normal-approx test of tau = 0
```

`kendall_tau_b` normalizes by the geometric mean of the untied-pair counts in each
margin, so an aligned tie no longer caps the coefficient below one;
`goodman_kruskal_gamma` drops tied pairs outright. The test uses the large-sample
normal approximation with `Var(S) = n(n-1)(2n+5)/18`, returning a two-sided p-value
equal to `erfc(|z| / sqrt(2))`.

The *unnormalized* count behind Kendall's tau is the number of discordant pairs — the
Kendall-tau *distance* between two rankings — and `count_inversions` computes it (and the
plain out-of-order-pair count of any sequence) in `O(n log n)` via merge sort:

```python
from quantforge import count_inversions, kendall_tau_distance, is_sorted

count_inversions([3, 1, 4, 1, 5, 9, 2, 6])   # 8 — pairs out of order
count_inversions([5, 4, 3, 2, 1])            # 10 — fully reversed, n(n-1)/2
kendall_tau_distance(['a','b','c','d'], ['b','a','d','c'])   # 2 — pairs ordered oppositely
is_sorted([1, 2, 2, 3])                      # True
```

`count_inversions` counts pairs `i < j` with `a[i] > a[j]` — 0 when sorted, `n(n-1)/2` when
reversed — piggybacking on a merge sort so it handles hundreds of thousands of elements a
quadratic scan could not. `kendall_tau_distance` reindexes one ranking by the other's
positions and counts the resulting inversions, giving the minimum adjacent swaps to turn
one ranking into the other. `is_sorted` is the monotonicity check (`strict=True` for
strictly increasing). Verified against brute pair counts over thousands of sequences and
permutation pairs.

A raw correlation can be entirely an artifact of a shared driver; `partial_correlation`
measures the link between two variables *after* regressing out one or more controls:

```python
from quantforge import partial_correlation, semipartial_correlation

# X and Y both driven by a common Z -> raw corr ~0.92, but nothing once Z is removed
partial_correlation(x, y, z)            # ~0.0 — the association was spurious
semipartial_correlation(x, y, z)        # control Z out of y only
```

`partial_correlation` regresses the controls out of *both* variables (OLS) and correlates
the residuals — matching the textbook three-variable formula and dropping a
common-cause correlation of `0.92` to essentially zero, while leaving a genuine link
untouched when the control is irrelevant. `semipartial_correlation` removes the controls
from only one variable (the unique contribution of `x` to `y`). `controls` is a single
column or a list of columns.

Kendall and Spearman still only see *monotone* co-movement. **Distance correlation**
(Székely-Rizzo) sees *any* dependence: it is zero if and only if the two variables
are independent, so it catches nonlinear structure a rank or linear correlation
misses entirely.

```python
from quantforge import distance_correlation, distance_covariance, distance_variance

distance_correlation([1, 2, 3, 4, 5], [2, 4, 6, 8, 10])   # 1.0  — tight linear

x = [-2, -1, 0, 1, 2]
y = [4, 1, 0, 1, 4]        # y = x^2, symmetric -> Pearson is exactly 0
distance_correlation(x, y)                                # 0.5159 — dependence seen
```

It double-centers each sample's pairwise-distance matrix and takes the mean product;
`distance_correlation` lands in `[0, 1]` (0 under independence, 1 for a tight linear
relation), `distance_covariance` is the unnormalized version, and
`distance_variance` is `dCov(x, x)`. On the symmetric parabola above, Pearson
correlation is exactly zero — it is blind to the U-shape — while distance
correlation reports a clear 0.52, flagging the dependence.

`chatterjee_xi` is a newer (2020) rank coefficient built for the same job but with an
`O(n log n)` cost and an asymmetric reading — it measures how far `Y` is a *function*
of `X`, running from 0 under independence toward 1 as `Y` becomes a noiseless function
of `X`, monotone or not:

```python
from quantforge import chatterjee_xi, blomqvist_beta, kendall_tau_b

x = [1, 2, 3, 4, 5, 6, 7]
y = [3, 2, 1, 0, 1, 2, 3]        # a V — a clean function of x, but non-monotone

kendall_tau_b(x, y)              # 0.0    — rank correlation cancels on the V
chatterjee_xi(x, y)             # 0.25   — sees the functional dependence
chatterjee_xi(list(range(1000)), list(range(1000)))   # 0.997 -> 1 for a monotone map

blomqvist_beta([1, 2, 3, 4, 5], [10, 20, 30, 40, 50])  # 1.0  comonotone
```

Chatterjee's `xi` climbs toward 1 as the sample grows whenever `Y = f(X)` exactly
(the V above would too, at larger `n`), which is what makes it a functional-dependence
detector rather than a monotonicity one. `blomqvist_beta` is the median-based
"medial" correlation — the rescaled fraction of points in the concordant quadrants
around the two medians — a fast, robust cousin of Kendall's tau in `[-1, 1]`.

For the linear (Pearson) correlation with a significance test and interval,
`pearson_correlation_test` returns the coefficient, a two-sided t-test of
`rho = 0`, and a Fisher-z confidence interval:

```python
from quantforge import pearson_r, pearson_correlation_test

pearson_r(x, y)                                  # coefficient only
res = pearson_correlation_test(x, y, confidence=0.95)
res["r"], res["t_stat"], res["p_value"]          # t = r sqrt((n-2)/(1-r^2))
res["conf_int"]                                  # [low, high] via atanh(r)
```

The t-test uses the library's Student-t CDF and the interval Fisher's
variance-stabilizing transform; perfect `+/-1` correlation collapses the interval
to the point.

Correlation says nothing about whether two assets crash *together* — tail
dependence does. `upper_tail_dependence` / `lower_tail_dependence` estimate
`P(U > q | V > q)` on the rank scale, and `exceedance_correlation` is the
correlation on joint-tail observations only.

```python
from quantforge import (upper_tail_dependence, lower_tail_dependence,
                        exceedance_correlation)

upper_tail_dependence(x, y, q=0.9)      # ~1 if extremes cluster, ~(1-q) if independent
lower_tail_dependence(x, y, q=0.1)      # joint-downside clustering
exceedance_correlation(x, y, q=0.9, tail="upper")   # correlation in the tail only
```

A comonotone pair has tail dependence ~1; independent margins sit at the `1 - q`
null and fall toward 0 at more extreme thresholds; a common-shock pair shows
clear positive tail dependence (~0.5) even when its bulk correlation is modest.

## Gaussian-copula sampling

Simulate a dependent portfolio: draw uniforms with a target correlation from a
Gaussian copula, then push each margin through its inverse-CDF to get correlated
draws from any distributions. `gaussian_copula_sample` correlates standard normals
via the Cholesky factor and maps them back through the normal CDF.

```python
from quantforge import gaussian_copula_sample, inverse_transform
from quantforge.mathfns import norm_ppf

R = [[1.0, 0.7], [0.7, 1.0]]
u = gaussian_copula_sample(R, n=5000, seed=42)   # list of [u1, u2] uniforms

# Turn margin 0 into standard-normal draws (or any inverse-CDF).
x = inverse_transform([row[0] for row in u], norm_ppf)
```

Each margin is uniform; the cross-margin rank correlation tracks `R` (≈0.69 for a
0.7 target). A non-positive-definite correlation matrix is rejected by the
Cholesky step.

The Gaussian copula has **zero tail dependence** — extremes decouple, which
famously understates joint-crash risk. `student_t_copula_sample` keeps the same
rank correlation but adds symmetric tail dependence that grows as the degrees of
freedom fall:

```python
from quantforge import student_t_copula_sample

student_t_copula_sample(R, df=4, n=5000)     # heavy joint tails
student_t_copula_sample(R, df=100, n=5000)   # ~ Gaussian copula
```

A lower `df` makes joint extremes markedly more likely at the same correlation; as
`df` grows the sampler converges to the Gaussian copula.

## Directional statistics

Angles need their own statistics — the ordinary mean of 350° and 10° is a meaningless
180°, not 0°. These work with the unit-vector sum instead:

```python
import math
from quantforge import circular_mean, resultant_length, circular_std, rayleigh_test

circular_mean([math.radians(350), math.radians(10)])   # 0.0 rad — wraps correctly
resultant_length([1.0, 1.0, 1.0])                       # 1.0 — perfectly concentrated
rayleigh_test(angles)                                   # (R, p_value) for non-uniformity
```

`circular_mean` averages the `(cos, sin)` vectors and takes their `atan2`, so it wraps
across the `2*pi` boundary (and raises when the vectors cancel to no mean direction).
`resultant_length` is the concentration `R` in `[0, 1]` — `1` for identical angles, `~0`
for a uniform spread — and `circular_variance` (`1 - R`) and `circular_std`
(`sqrt(-2 ln R)`) follow from it. `rayleigh_test` tests whether the angles have a preferred
direction at all: on a tight cluster it returns `R ≈ 0.98` with `p ≈ 0` (reject
uniformity), and on uniform data `p > 0.9`. Angles are in radians.

The von Mises distribution is the circular analogue of the normal — the bell curve on a
circle:

```python
from quantforge import von_mises_pdf, von_mises_fit, bessel_i0

von_mises_pdf(0, mu=0, kappa=2)      # 0.51589 — density at the mean
mu_hat, kappa_hat = von_mises_fit(angles)   # MLE mean direction + concentration
```

`von_mises_pdf(theta, mu, kappa)` is `exp(kappa cos(theta - mu)) / (2 pi I0(kappa))` — it
integrates to 1 over the circle, peaks at `mu`, and flattens to the uniform density at
`kappa = 0`. `von_mises_fit` returns the maximum-likelihood mean (the sample circular
mean) and concentration (solving `I1(kappa)/I0(kappa) = R` by Newton's method); on a
simulated `(mu=0.5, kappa=4)` sample it recovers `(0.51, 4.15)`. `bessel_i0` and
`bessel_i1` (the modified Bessel functions the normalizer needs) are exposed too.

## Geodesy

Great-circle distances and bearings between latitude/longitude points on a spherical
Earth:

```python
from quantforge import (haversine_distance, initial_bearing,
                        destination_point, cross_track_distance)

haversine_distance(51.5074, -0.1278, 48.8566, 2.3522)   # 343.6 km (London -> Paris)
initial_bearing(0, 0, 0, 10)                             # 90.0 — due east
destination_point(51.5, 0, bearing=90, distance=111)     # ~ (51.49, 1.60)
```

`haversine_distance` is the numerically stable great-circle distance (kilometres on the
IUGG mean radius by default), matching known city pairs (NYC-LA ~3936 km).
`initial_bearing` gives the forward compass heading in `[0, 360)`, `destination_point`
follows a bearing and distance to the arrival point (the exact inverse — a round-trip
recovers the input distance and bearing), and `cross_track_distance` is the signed
perpendicular offset of a point from a route (positive to the right of the path), for
"how far off course am I". All angles in degrees.

## Quaternions

Unit quaternions represent 3-D rotations without gimbal lock and interpolate smoothly
where matrices cannot:

```python
import math
from quantforge import axis_angle_to_quat, rotate_vector, slerp, quat_multiply

q = axis_angle_to_quat([0, 0, 1], math.pi / 2)   # 90 deg about z
rotate_vector(q, (1, 0, 0))                        # (0, 1, 0) — x rotates to y
mid = slerp(axis_angle_to_quat([0,0,1], 0), q, 0.5)
rotate_vector(mid, (1, 0, 0))                      # (0.707, 0.707, 0) — halfway (45 deg)
```

A quaternion is `(w, x, y, z)`; `axis_angle_to_quat` builds one from an axis and angle,
and `rotate_vector` applies it as `q v q*` (length-preserving). `quat_multiply` composes
rotations — `quat_multiply(a, b)` applies `b` then `a`, matching the sequential rotation —
and `quat_to_axis_angle` inverts the construction. `slerp` interpolates along the shortest
constant-speed arc between two orientations, staying unit-norm, so its midpoint is exactly
the half-angle rotation. `quat_normalize` and `quat_conjugate` round out the set.

A rotation also has matrix and Euler-angle forms, and the `rotation` conversions move
between all three:

```python
from quantforge import quat_to_matrix, euler_to_quat, quat_to_euler
import math

quat_to_matrix((math.cos(math.pi/4), 0, 0, math.sin(math.pi/4)))
                                # the 90-deg-about-z matrix ([0,-1,0],[1,0,0],[0,0,1])
euler_to_quat(math.pi/2, 0, 0)  # (0.7071, 0, 0, 0.7071) — yaw 90 deg
quat_to_euler((0.7071, 0, 0, 0.7071))   # (1.5708, 0, 0) = (yaw, pitch, roll)
```

`quat_to_matrix`/`matrix_to_quat` (Shepperd's method) convert to and from a 3x3 rotation
matrix; `euler_to_quat`/`quat_to_euler` use the aerospace Z-Y-X (yaw, pitch, roll) intrinsic
convention; and `euler_to_matrix`/`matrix_to_euler` compose the two. Quaternions are the
same `(w, x, y, z)` used above, so `quat_to_matrix(q)` applied to a vector matches
`rotate_vector(q, v)`. Verified over thousands of random rotations that every round-trip
recovers the input and the matrices stay orthogonal.

A rotation can also be written as a single *rotation vector* (axis scaled by angle) — the
SO(3) exponential map. `rodrigues` exponentiates it to a matrix and `so3_log` inverts:

```python
from quantforge import rodrigues, so3_log, hat
import math

rodrigues([0, 0, math.pi/2])       # the 90-deg-about-z matrix
so3_log(rodrigues([0.1, 0.2, 0.3]))# [0.1, 0.2, 0.3] — recovers the rotation vector
hat([1, 2, 3])                     # [[0,-3,2],[3,0,-1],[-2,1,0]] — skew matrix
```

`rodrigues(omega)` applies Rodrigues' formula `R = I + sin(t) K + (1 - cos t) K^2` where `t`
is the rotation angle `|omega|` and `K = hat(omega/t)` is the skew-symmetric cross-product
matrix; `so3_log` takes a matrix back to its rotation vector (with dedicated near-zero and
near-pi branches). `hat`/`unhat` convert between a 3-vector and its skew matrix, with
`hat(v) w == cross(v, w)`. This exp/log pair is the natural way to average, interpolate, or
differentiate rotations on the SO(3) manifold — and it agrees with the quaternion path to
machine precision.

The everyday vector operations back them up — dot and cross products, norms, angles,
projection, and reflection:

```python
from quantforge import cross, angle_between, vector_project, reflect
import math

cross((1, 0, 0), (0, 1, 0))                    # (0, 0, 1)
math.degrees(angle_between((1, 0, 0), (0, 1, 0)))   # 90.0
vector_project((3, 3, 0), (1, 0, 0))           # (3.0, 0.0, 0.0) — component along x
reflect((1, 2, 3), (0, 0, 1))                  # (1, 2, -3) — mirror through the xy-plane
```

`dot`, `norm`, `normalize`, `angle_between`, `vector_project`, and `vector_reject` work in
any dimension; `cross` is 3-D (anticommutative and perpendicular to both inputs).
`angle_between` uses the stable `atan2(|a×b|, a·b)` form in 3-D. `vector_project` and
`vector_reject` split a vector into components along and perpendicular to another (they sum
back to the original), and `reflect` mirrors a vector through the plane with a given
normal, preserving its length. The projection helpers are prefixed `vector_` to avoid
clashing with the PCA `project`.

## Color spaces

Conversions between the common color representations — RGB, the cylindrical HSV/HSL
models, and hex:

```python
from quantforge import rgb_to_hsv, rgb_to_hsl, rgb_to_hex, hex_to_rgb

rgb_to_hsv(1, 0, 0)          # (0.0, 1.0, 1.0) — pure red
rgb_to_hsl(0, 0, 1)          # (240.0, 1.0, 0.5) — pure blue
rgb_to_hex(1, 0, 0)          # "#ff0000"
hex_to_rgb("#ff8800")        # (1.0, 0.533, 0.0)
```

RGB channels are floats in `[0, 1]`; HSV/HSL hue is in degrees `[0, 360)` with the other
channels in `[0, 1]`. `rgb_to_hsv`/`rgb_to_hsl` (and their inverses `hsv_to_rgb`/
`hsl_to_rgb`) match Python's `colorsys` and round-trip exactly, and `rgb_to_hex`/
`hex_to_rgb` convert to and from the `"#rrggbb"` web form.

## Spectral analysis

Find cyclical structure — a seasonal pattern, a dominant trading cycle — with the
discrete Fourier transform and periodogram. `dominant_frequency` returns the
strongest non-DC frequency; its reciprocal is the period in samples.

```python
from quantforge import periodogram, dominant_frequency, spectral_energy

freqs, power = periodogram(series)     # one-sided power spectrum, [0, 0.5]
f = dominant_frequency(series)         # cycles per sample of the strongest peak
1 / f                                  # dominant period in samples
spectral_energy(series)                # equals sum(x^2) by Parseval
```

A pure sinusoid peaks exactly at its frequency; a constant series has power only
at zero frequency; white noise spreads its power evenly with no dominant peak.

Compact signal descriptors summarize a series for classification or monitoring — some in
the time domain, some read off the spectrum:

```python
from quantforge import (zero_crossing_rate, rms, crest_factor,
                        spectral_centroid, spectral_bandwidth, spectral_flatness)

rms(sine)                 # 0.707 — 1/sqrt(2) for a unit sine
crest_factor(sine)        # 1.414 — peak / RMS (~1 for a square wave)
spectral_centroid(tone)   # ~ the tone's frequency (spectral center of mass)
spectral_flatness(tone)   # ~0.003 tonal; ~0.55 for white noise
```

`zero_crossing_rate` (sign changes per sample) is a cheap noisiness/pitch proxy; `rms`
and `crest_factor` describe amplitude and peakiness. In the frequency domain,
`spectral_centroid` is the power-weighted mean frequency, `spectral_bandwidth` its spread
(narrow for a tone, wide for noise), and `spectral_flatness` (Wiener entropy) is near `0`
for tonal signals and near `1` for flat, noise-like spectra.

The periodogram uses a direct DFT; for the fast transform itself, `fft` / `ifft` give
the radix-2 Cooley-Tukey pair (power-of-two lengths):

```python
from quantforge import fft, ifft

X = fft(signal)        # O(n log n) complex spectrum
ifft(X)                # recovers the signal (1/N scaled)
```

`fft` matches the direct DFT to machine precision but scales as `O(n log n)`;
`ifft(fft(x))` round-trips. Handy for fast convolution and any spectral transform
where the length is a power of two.

When the length is *not* a power of two, `dft_any`/`idft_any` still run in `O(n log n)` via
Bluestein's chirp-Z algorithm:

```python
from quantforge import dft_any, idft_any

dft_any([1, 2, 3, 4, 5])           # length-5 DFT: [15, -2.5+3.44j, -2.5+0.81j, ...]
idft_any(dft_any([1, 2, 3, 4, 5])) # recovers [1, 2, 3, 4, 5]
```

Bluestein rewrites each DFT output as a chirp multiply, a convolution (padded up to a power
of two and run through the radix-2 `fft`), and a final chirp multiply — so any length works,
including primes, where `fft` cannot apply. It agrees with `fft` exactly on power-of-two
lengths and with a direct DFT everywhere. Verified against a direct DFT and round-trip over
hundreds of arbitrary-length sequences.

When the inputs are integers and you need the *exact* answer with no rounding, `convolve_mod`
runs the number-theoretic transform — the same FFT butterflies but in a finite field —
returning the convolution with every coefficient reduced modulo a prime:

```python
from quantforge import convolve_mod, ntt, intt, NTT_PRIME

convolve_mod([1, 2, 3], [4, 5, 6])         # [4, 13, 28, 27, 18] — exact, no float error
convolve_mod([1, 1, 1, 1], [1, 1, 1, 1])   # [1, 2, 3, 4, 3, 2, 1]
intt(ntt([3, 1, 4, 1]))                     # [3, 1, 4, 1] — round-trips exactly
NTT_PRIME                                   # 998244353 = 119 * 2^23 + 1
```

`ntt`/`intt` transform a power-of-two-length integer sequence using a primitive root of
unity modulo `998244353` (primitive root 3), so there is no floating-point error at all —
every value is an exact residue. `convolve_mod` pads to a power of two, multiplies pointwise
in the transform domain, and inverts, giving the polynomial product (or big-integer
convolution) mod the prime in `O(n log n)`. Verified against a direct integer convolution
over 3000 random sequence pairs.

The FFT diagonalizes *cyclic* convolution; the Walsh-Hadamard transform does the same for
*XOR* convolution, and companion transforms handle AND and OR. `fwht`/`ifwht` and the three
`*_convolve` helpers cover all of it:

```python
from quantforge import fwht, xor_convolve, and_convolve, or_convolve

fwht([1, 0, 0, 0])                          # [1, 1, 1, 1] — impulse transforms to flat
xor_convolve([1, 2, 3, 4], [5, 6, 7, 8])    # [70, 68, 62, 60]  (sum over i^j == k)
and_convolve([1, 2, 3, 4], [5, 6, 7, 8])    # [103, 52, 73, 32] (sum over i&j == k)
or_convolve([1, 2, 3, 4], [5, 6, 7, 8])     # [5, 28, 43, 184]  (sum over i|j == k)
```

`xor_convolve` computes `c[k] = sum_{i ^ j = k} a[i] b[j]` by transforming both inputs with
the Walsh-Hadamard transform, multiplying pointwise, and inverting — `O(n log n)` for a
power-of-two length, and exact for integer inputs. `and_convolve` and `or_convolve` use the
superset and subset zeta (Mobius) transforms for the AND and OR index rules — the tools
behind subset-sum DP and bitmask counting. Verified against brute-force bitwise-index
enumeration over 4000 random sequence pairs for each operation.

Two FFT-backed conveniences build on it:

```python
from quantforge import convolve, fft_autocorrelation

convolve([1, 2, 3], [4, 5, 6, 7])     # full linear conv / polynomial product
fft_autocorrelation(series, max_lag=20)   # autocorrelation, acf[0] = 1
```

`convolve` is the coefficient product of two polynomials in `O(n log n)` (matching a
direct convolution), and `fft_autocorrelation` computes the whole autocorrelation via
Wiener-Khinchin — far faster than the direct lag-by-lag sum on long series.

Where autocorrelation slides a series against *itself*, cross-correlation slides *two*
series past each other and its peak locates the lag at which they line up — the standard
way to estimate a time delay or a lead-lag relationship:

```python
from quantforge import cross_correlation, lag_at_max_correlation

lags, values = cross_correlation(x, y, max_lag=20)   # raw, over -20..+20
lag, coef = lag_at_max_correlation(x, y)             # lag of best alignment, in [-1,1]
```

`cross_correlation` returns the raw overlap at every lag (matching a direct sum to
floating error, computed via one FFT conjugate-multiply); `normalized_cross_correlation`
divides out the amplitudes to give a correlation coefficient in `[-1, 1]`; and
`lag_at_max_correlation` reports the peak. If `y` is `x` delayed by five samples, the
peak lands at lag `5` with coefficient `≈0.90`; a positive lag means `x` leads `y`.
Identical series peak at lag `0` with coefficient `1`.

When the pattern you are hunting for is a *known* waveform rather than another series, the
matched filter is the optimal detector — it correlates the data against the template and
peaks where the two best align, maximizing signal-to-noise ratio in white noise:

```python
from quantforge import matched_filter, find_peaks, detect_template

mf = matched_filter(x, template)             # detection statistic at each offset
find_peaks(mf, height=10, distance=20)       # local maxima, filtered by height/spacing
detect_template(x, template, threshold=0.9)  # start offsets of normalized matches
```

Plant a template in noise and `matched_filter`'s peak lands exactly at its offset;
`normalized_matched_filter` scales the response to a correlation coefficient that hits `1`
at a perfectly scaled match. `find_peaks` returns strict local maxima, dropping any below
a `height` and keeping the taller of two within `distance` samples. `detect_template`
combines them to report every offset where the template recurs — planting a distinctive
16-sample waveform at offsets `[50, 160, 300]` in noise recovers exactly those three.
Because the correlation is amplitude-normalized, short or featureless templates match many
noise windows by shape alone, so use a longer, distinctive template when false positives
matter.

The real-valued cousin of the FFT, the discrete cosine transform, is `dct` (DCT-II) with
inverse `idct` (DCT-III):

```python
from quantforge import dct, idct

X = dct(signal)        # orthonormal cosine spectrum
idct(X)                # recovers the signal exactly
```

The DCT expresses a real signal as a sum of cosines and, being orthonormal, preserves
energy and round-trips to machine precision. It packs most of a smooth signal's energy
into the first few coefficients (a smooth curve puts >95% in its first four) — the
compaction property behind JPEG and MP3 — so it doubles as a compression and
denoising basis. Unlike the FFT it needs no power-of-two length.

The *cepstrum* is the inverse transform of the log-magnitude spectrum. It turns
periodicity *in the spectrum* — the evenly spaced harmonics of a pitched sound, or the
ripple a single echo adds — into a peak at one "quefrency" (a time, in samples):

```python
from quantforge import real_cepstrum, fundamental_quefrency

real_cepstrum(signal)                                   # ifft(log|fft|), read like a lag
q, peak = fundamental_quefrency(signal, min_quefrency=10)  # dominant quefrency
```

A signal with an echo delayed by 50 samples produces a cepstral peak exactly at
quefrency `50`, and a periodic impulse train (the voiced-speech model) peaks at its
period — divide the sample rate by that quefrency to get pitch in Hz. `power_cepstrum`
squares the log-spectrum before inverting to emphasize those peaks. `min_quefrency`
skips the low-quefrency spectral envelope so the search finds the pitch/echo peak, not
the DC hump. This is the classic tool for fundamental-frequency estimation and echo
detection.

When you only need *one* frequency's strength — tone detection, a known harmonic — the
Goertzel algorithm computes a single DFT bin in `O(n)` without a full transform:

```python
from quantforge import goertzel, goertzel_power

goertzel(signal, k)          # complex X[k], identical to the k-th FFT output
goertzel_power(signal, k)    # |X[k]|^2 — large iff a frequency near bin k is present
```

`goertzel` matches the FFT coefficient exactly (both magnitude and phase) with a single
`O(n)` real recurrence — far cheaper than an FFT when you care about a handful of bins.
On a pure tone at bin 8 the power spikes there (`1024` for a 64-sample cosine) and is
essentially zero elsewhere, which is exactly how DTMF and pilot-tone detectors work.

Cutting a finite segment out of a signal leaks energy across the spectrum; tapering it
first with a window function cuts that leakage. The common windows and an `apply_window`
helper are provided:

```python
from quantforge import hann, apply_window, fft

windowed = apply_window(segment, "hann")   # taper before the FFT
X = fft(windowed)                            # far lower sidelobes than the raw segment
hann(256)                                    # the window itself, if you want it
```

`hann`, `hamming`, `blackman`, `bartlett` and `rectangular` return length-`n` window
vectors (all symmetric, tapering to near zero at the edges except the boxcar);
`apply_window` multiplies a signal by a named window or a custom list. On an off-bin
sinusoid a Hann taper drops a distant sidelobe from ~1.0 to ~0.0003 — the standard
pre-processing before a periodogram, and what `welch_psd` applies internally.

The raw periodogram is noisy; `welch_psd` averages windowed segment periodograms for
a much lower-variance spectral-density estimate:

```python
from quantforge import welch_psd

freqs, power = welch_psd(series, segment_length=256)   # smoothed one-sided PSD
```

Welch splits the series into overlapping Hann-windowed segments and averages their
periodograms — trading some frequency resolution for far less variance, so a
spectral peak stands out cleanly against a noisy background where the raw
periodogram would bury it.

Welch averaging also drives the *cross*-spectrum between two signals: `coherence` is the
frequency-domain squared correlation, `1` where they are perfectly linearly related at a
frequency and `0` where unrelated:

```python
from quantforge import coherence, cross_spectral_density

freqs, coh = coherence(x, y, segment_length=256)   # in [0, 1] per frequency
freqs, csd = cross_spectral_density(x, y)          # complex shared spectrum
```

Two signals sharing a 0.1-cycle component reach coherence `~1.0` in that band while
independent noise averages near `0.06` (a nonzero floor set by the segment count — a
single segment gives a meaningless coherence of 1 everywhere, so averaging is essential).
`cross_spectral_density` returns the complex cross-spectrum whose magnitude peaks at the
shared frequency and whose phase gives the frequency-dependent lead/lag.

The periodogram and Welch estimate a spectrum straight from the data; a *parametric*
estimate instead fits an autoregressive model and reads its spectrum off the
coefficients, which resolves sharp peaks far better on a short record:

```python
from quantforge import burg, ar_spectrum

model = burg(series, order=2)                    # forward-backward AR fit
freqs = [i / 500 for i in range(251)]            # 0..0.5 cycles/sample
psd = ar_spectrum(series, order=20, freqs=freqs) # parametric PSD (Burg by default)
```

`burg` estimates the AR coefficients directly from the data (never forming
autocovariances), so on 2000 samples of a known AR(2) process it recovers the
coefficients to within a few percent, and its reflection coefficients stay inside the
unit circle (a stable model). Fed just 128 samples of two sinusoids in noise,
`ar_spectrum` puts clean peaks exactly at both frequencies where a periodogram of the
same length would smear them together. `ar_psd` evaluates the spectrum from coefficients
you already have, and `method="yule_walker"` selects the autocovariance route instead.

Welch *averages away* time to get a cleaner spectrum; the short-time Fourier transform
*keeps* it, taking an FFT of each overlapping windowed frame to show how the spectrum
evolves. Its squared magnitude is the spectrogram — the standard time-frequency view of
speech, music, or any signal whose content shifts:

```python
from quantforge import stft, spectrogram, istft

frames = stft(x, 256, hop=64)            # per-frame complex spectra
S = spectrogram(x, 256, hop=128)         # |STFT|^2, each frame DC..Nyquist
xr = istft(frames, 256, hop=64)          # overlap-add reconstruction
```

`stft` slides a `frame_size` window (a power of two) in steps of `hop` (default 50%
overlap), so a 2048-sample signal gives 8 frames at `hop=256` and 16 at `hop=128`.
`spectrogram` returns `frame_size//2 + 1` power bins per frame (a pure tone lands
entirely in its own bin; a rising chirp's peak bin climbs frame by frame). `istft`
inverts by weighted overlap-add, dividing out the summed squared window, so
`istft(stft(x))` reconstructs the interior to ~1e-14 regardless of the overlap — which
means you can filter in the time-frequency plane and reconstruct.

Where windows shape a *measurement* of the spectrum, a FIR filter reshapes the
*signal* — keeping some frequencies and rejecting others. The windowed-sinc design
takes the ideal brick-wall impulse response and tapers it with one of the same
windows to get a finite, well-behaved filter:

```python
from quantforge import fir_lowpass, fir_highpass, fir_bandpass, fir_apply

taps = fir_lowpass(51, 0.2)      # 51-tap low-pass, cutoff 0.2 x Nyquist
sum(taps)                         # 1.0 — unit DC gain
filtered = fir_apply(taps, x)     # convolve with the signal

fir_highpass(51, 0.2)             # complementary high-pass (zero DC gain)
fir_bandpass(101, 0.15, 0.35)     # passes only the band between the two cutoffs
```

Cutoffs are normalized to the Nyquist frequency, so `0.2` means one fifth of the way
to half the sample rate. The low-pass has unit DC gain (`sum(taps) == 1.0`) and is
symmetric, so it is linear-phase — it delays every frequency equally rather than
smearing the waveform. It passes what is below the cutoff (gain `0.996` at `0.025`)
and crushes what is above it (`0.0008` at `0.2`); the high-pass is the mirror image,
with zero DC gain and near-unit gain in its passband. Fed a two-tone signal of a slow
and a fast cosine, a low-pass keeps the slow one and removes the fast one, leaving an
amplitude near the surviving tone's `1.0`. `fir_bandpass` builds a band-pass as the
difference of two low-pass filters and `fir_apply` runs the direct convolution.

A FIR filter needs many taps for a sharp cutoff; an *IIR* filter reaches the same
selectivity with a handful of coefficients by feeding output back — trading linear phase
for efficiency. The Butterworth response is maximally flat in the passband (no ripple),
designed here via the bilinear transform as a cascade of second-order sections (biquads)
for numerical stability:

```python
from quantforge import butter_lowpass, sosfilt, iir_frequency_response

sos = butter_lowpass(4, 0.1)          # 4th-order low-pass, -3 dB at 0.1 cycles/sample
y = sosfilt(sos, x)                    # filter the signal
iir_frequency_response(sos, [0.1])    # -> ~0.7071, the -3 dB point
```

A 4th-order low-pass comes out as 2 biquads with unit DC gain, exactly `0.7071` (-3 dB)
at the cutoff, and a monotone rolloff into a deep stopband (a `0.3`-cycle tone comes
through at amplitude `0.003`). `butter_highpass` is the mirror — zero DC gain, unit gain
at Nyquist. Higher order rolls off faster, and every section is stable (poles inside the
unit circle) for both even and odd orders. `sosfilt` runs the cascade in Direct Form II
transposed and `iir_frequency_response` gives the magnitude response for checking a
design.

The same low-pass filters drive *sample-rate conversion*, which cannot be done by naively
dropping or repeating samples — that folds high frequencies back as aliases or leaves
spectral images. The correct operations filter at the new Nyquist limit:

```python
from quantforge import sinc_interp, upsample, downsample, resample_rational

sinc_interp(x, [t + 0.5 for t in range(len(x) - 1)])  # band-limited reconstruction
upsample(x, 4)                                         # 4x rate, anti-imaging filter
downsample(x, 2)                                       # 1/2 rate, anti-alias filter
resample_rational(x, 3, 2)                             # rational 3/2 ratio
```

`sinc_interp` is the ideal Whittaker-Shannon reconstruction — exact at integer positions
and the band-limited interpolation between them. `upsample` inserts zeros and filters
(output `factor` times as long); `downsample` filters then keeps every `factor`-th sample
(so a length-256 signal downsampled by 2 has 128 samples). The anti-alias filter is the
point: a tone at `0.32` cycles/sample, above the new Nyquist of `0.25` after halving the
rate, comes out with rms `0.001` — suppressed, not aliased back to a false low frequency
as naive decimation would. Each operation carries the FIR filter's linear-phase group
delay of `(numtaps-1)/2` samples.

Linear filters smear sharp features and get dragged by outliers; *order-statistic*
filters sort a sliding window and pick one sample by rank, rejecting impulsive spikes
while keeping edges sharp:

```python
from quantforge import median_filter, rank_filter, hampel_filter

median_filter(x, 3)                       # sliding median — de-spikes, keeps edges
rank_filter(x, 5, percentile=100)         # sliding max (0 = min, 50 = median)
clean, outliers = hampel_filter(x, window=5)   # MAD-based outlier replacement
```

A single spike in `[1, 1, 1, 50, 1, 1, 1]` vanishes under `median_filter(x, 3)` while a
`[0,0,0,0,0,10,10,10,10,10]` step stays perfectly sharp — a moving average would round
both off. `rank_filter` generalizes to any percentile (min, max, or a quantile envelope).
`hampel_filter` is the surgical option: it computes the local median and a robust
standard deviation from the MAD, and replaces *only* the points more than `n_sigmas`
(default 3) away, returning both the cleaned signal and the outlier indices — inject one
spike into otherwise clean data and it flags exactly that index, leaving everything else
byte-for-byte unchanged.

For a step-like signal buried in noise, *total-variation denoising* recovers the
piecewise-constant shape by minimizing `(1/2)*sum((x-y)^2) + lam*sum|x_{k+1}-x_k|` — the
L1 penalty flattens noise into constant runs while leaving jumps intact:

```python
from quantforge import tv_denoise, tv_total_variation

clean = tv_denoise(noisy, lam=3.0)   # exact minimizer, O(n) (Condat's method)
tv_total_variation(clean)            # sum |x_{k+1} - x_k|, the penalized quantity
```

`lam=0` returns the input; a large `lam` collapses everything to the mean. On a noisy
0-to-10 step it recovers plateaus at `0.01` and `9.75` and cuts the total variation from
`46.7` to `9.7` while keeping the edge sharp. Condat's algorithm returns the *exact*
minimizer in a single O(n) pass — verified against a brute-force optimality check and an
independent subgradient solver — so unlike an iterative smoother there is nothing to tune
but `lam`.

Filtering keeps and rejects frequencies; the *analytic signal* instead turns a real
oscillation into a rotating phasor, so you can read off its instantaneous amplitude and
frequency. The Hilbert transform is the machinery: a 90-degree phase shift of every
component (cosine becomes sine), computed by zeroing the negative-frequency half of the
spectrum and doubling the positive half:

```python
from quantforge import analytic_signal, envelope, instantaneous_frequency

z = analytic_signal(x)            # complex: z.real == x, z.imag == Hilbert(x)
env = envelope(am_signal)         # instantaneous amplitude — recovers an AM envelope
freq = instantaneous_frequency(x) # cycles/sample, from the phase derivative
```

For an amplitude-modulated carrier `(1 + 0.5*cos(2*pi*fm*t)) * cos(2*pi*fc*t)`,
`envelope` recovers the `1 + 0.5*cos(...)` modulating waveform (the fast carrier drops
out entirely). For a pure tone at three cycles over 256 samples the real part of the
analytic signal is the input to floating error, and `instantaneous_frequency` is flat at
`3/256 ≈ 0.0117` cycles/sample. `hilbert_transform` returns just the imaginary part and
`instantaneous_phase` the unwrapped phase. Inputs must have a power-of-two length (the
FFT constraint). This is the standard demodulation and envelope-detection tool for AM/FM
signals, vibration analysis, and empirical-mode decomposition.

## Wavelet transform (Haar multiresolution)

Where the Fourier transform asks *what frequencies are present*, the wavelet
transform asks *what happens at each scale, and where*. The Haar transform is the
simplest orthonormal wavelet: at each level it replaces adjacent pairs by their
scaled sum (a coarse approximation) and difference (the detail), then recurses on
the approximation — a multiresolution view of trend plus detail at every scale.
Length must be a power of two.

```python
from quantforge import haar_dwt, haar_idwt

approx, details = haar_dwt(series, levels=2)   # coarse trend + detail per level
haar_idwt(approx, details)                     # exact reconstruction
```

`haar_dwt` returns the final coarse approximation and a list of detail-coefficient
lists (finest level first); on a length-32 series with `levels=2` the finest detail
has 16 coefficients, the next has 8, and the approximation has 8. Being orthonormal,
`haar_idwt` reconstructs the input to machine precision (round-trip error ~1e-16).

A single level is just the scaled pairwise average and difference:

```python
from quantforge import haar_dwt

approx, details = haar_dwt([4.0, 2.0, 6.0, 8.0], levels=1)
approx      # [4.2426, 9.8995]  == (x0+x1)/sqrt2, (x2+x3)/sqrt2
details[0]  # [1.4142, -1.4142] == (x0-x1)/sqrt2, (x2-x3)/sqrt2
```

`wavelet_energy` decomposes total signal energy into the fraction sitting in each
detail level and in the coarse approximation — because the transform is orthonormal,
the fractions sum to one (Parseval):

```python
from quantforge import wavelet_energy

e = wavelet_energy(series, levels=2)
e["detail"]   # [0.0096, 0.0377] — energy fraction per level, finest first
e["approx"]   # 0.9527 — coarse-approximation fraction
```

A smooth series (a clean sinusoid) concentrates almost all its energy in the coarse
approximation with negligible detail; a noisy or bursty series spreads energy into
the fine detail levels — so the profile is a compact scale-by-scale summary of how
rough a signal is and where.

That energy split is exactly what makes wavelet *denoising* work. Donoho and
Johnstone's shrinkage transforms a noisy signal, shrinks the small detail
coefficients toward zero (noise spreads thinly across many small ones), and inverts —
the few large coefficients carrying real structure survive:

```python
from quantforge import wavelet_denoise

clean_estimate = wavelet_denoise(noisy, mode="soft")   # VisuShrink threshold, auto
```

With no threshold given it uses the VisuShrink universal threshold
`sigma * sqrt(2 log n)`, estimating the noise scale `sigma` robustly from the finest
detail level. On a length-1024 sinusoid plus N(0, 0.4) noise the mean-squared error
against the clean signal drops from about 0.164 (noisy) to 0.102 (denoised). `mode`
is `"soft"` (shrink) or `"hard"` (keep-or-kill); pass an explicit `threshold` to
override the automatic choice.

The pieces are exposed directly for custom shrinkage:

```python
from quantforge import mad_sigma, universal_threshold, soft_threshold, haar_dwt

approx, details = haar_dwt(noisy)
sigma = mad_sigma(details[0])            # median-absolute-deviation noise scale
lam = universal_threshold(len(noisy), sigma)   # ~1.53 here
soft_threshold(coefficient, lam)         # sign(x) * max(|x| - lam, 0)
```

`mad_sigma` uses the median absolute deviation (`median(|d|) / 0.6745`) so a few
large signal coefficients don't inflate the noise estimate the way a plain standard
deviation would — the robustness is what lets a single finest-level pass calibrate
the threshold for the whole signal.

## Structural breaks (CUSUM / Chow)

Detect when a mean or regime shifts. `cusum_mean` returns the standardized
cumulative sum of deviations from the mean — a stable series stays inside a
Kolmogorov band, a level shift drives it out. `chow_test` is an F test for a break
at a known point.

```python
from quantforge import cusum_break_detected, cusum_mean, chow_test

cusum_break_detected(series)          # True if the CUSUM path leaves its 95% band
cusum, band = cusum_mean(series)      # inspect the path and threshold directly

f, dof1, dof2 = chow_test(series, break_index=len(series)//2)
# large f rejects "no break at that point"
```

On a stable series the CUSUM stays well inside the band and the Chow F is small; a
persistent level shift breaches the band and drives the F statistic large.

## Cointegration (ADF / Engle-Granger)

Before trading a spread, test that it is actually mean-reverting. `adf_test` runs
the augmented Dickey-Fuller unit-root test; `engle_granger` regresses one series
on the other and ADF-tests the residual, returning the hedge ratio and a
cointegration verdict.

```python
from quantforge import adf_test, engle_granger

adf_test(series)["reject_5pct"]        # True if the series is stationary

c = engle_granger(y, x)
c["hedge_ratio"]                        # beta from regressing y on x
c["cointegrated_5pct"]                  # True if the residual spread is stationary
```

A stationary series gives a strongly negative ADF statistic and rejects the unit
root; two independent random walks do not. When a pair is cointegrated, feed the
residual spread straight into `fit_ornstein_uhlenbeck` to size the trade.

`kpss_test` is the complementary confirmatory test — its null is *stationarity*, so
a small p-value rejects it (the opposite of ADF):

```python
from quantforge import kpss_test

eta, p = kpss_test(series)                  # level stationarity ("c")
eta, p = kpss_test(series, regression="ct") # trend stationarity
```

Agreement between the two (ADF rejects the unit root *and* KPSS fails to reject
stationarity) is stronger evidence than either alone; disagreement flags a
borderline or trending series.

## Ornstein-Uhlenbeck calibration

## Ornstein-Uhlenbeck calibration

Calibrate the mean-reverting OU process `dX = kappa(theta - X)dt + sigma dW` from a
sampled path — the standard fit for a pairs spread or any mean-reverting series.
The discrete-time OU is an exact AR(1), so `fit_ornstein_uhlenbeck` fits that by
least squares and inverts the relations to recover speed, level, and volatility in
closed form.

```python
from quantforge import fit_ornstein_uhlenbeck

p = fit_ornstein_uhlenbeck(spread, dt=1/252)
p["kappa"]        # mean-reversion speed (per year)
p["theta"]        # long-run level
p["sigma"]        # instantaneous volatility
p["half_life"]    # ln(2) / kappa — time for a deviation to halve
```

A random walk fits a near-zero `kappa` (very long half-life); an oscillating,
anti-persistent series is rejected. Pair the half-life with `spread_zscore` and
`ou_half_life` to size a mean-reversion trade.

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

In continuous time the chain is specified by a rate generator `Q` (rows summing to
zero); the transition matrix over horizon `t` is `exp(Q t)`:

```python
from quantforge import generator_to_transition, generator_default_probability

Q = [[-0.15, 0.12, 0.03], [0.05, -0.20, 0.15], [0.0, 0.0, 0.0]]   # 2 grades + default
generator_to_transition(Q, t=1.0)                     # 1-year transition matrix
generator_default_probability(Q, default_state=2, horizons=[1, 2, 5, 10])
```

`P(t) = exp(Q t)` has the semigroup property `P(s) P(t) = P(s+t)`, so a single
generator prices migration at every horizon — the continuous-time basis for
rating-transition and default-term-structure modeling.

## K-means clustering

Partition points into `k` clusters by Lloyd's algorithm with k-means++ seeding —
regime grouping, basket construction, or any unsupervised segmentation.

```python
from quantforge import kmeans

r = kmeans(points, k=3, seed=42)
r["labels"]        # cluster index per point
r["centroids"]     # k cluster centers
r["inertia"]       # total within-cluster squared distance
```

Well-separated groups are recovered as pure clusters with centroids at their true
centers; runs are reproducible for a fixed seed and inertia falls as `k` rises
(use the elbow in inertia-vs-`k` to choose `k`).

To score cluster quality without labels, use the silhouette — `(b - a)/max(a, b)`
comparing each point's own-cluster cohesion to the nearest other cluster:

```python
from quantforge import silhouette_score

silhouette_score(points, r["labels"])   # ~1 tight/separated, ~0 overlapping, <0 misassigned
```

Picking the `k` that maximizes the silhouette is a common alternative to the
inertia elbow; a wrong `k` scores strictly lower on well-separated data.

When the cluster count is unknown or the groups aren't blob-shaped, `dbscan` clusters
by *density* instead — no preset `k`, arbitrary shapes, and explicit noise labels:

```python
from quantforge import dbscan

labels = dbscan(points, eps=1.0, min_samples=5)
# 0, 1, ... for clusters; -1 for noise (low-density outliers)
```

A point is a *core* point if at least `min_samples` neighbours lie within `eps`;
clusters grow by connecting core points, and anything unreachable is labelled `-1`.
It separates well-spaced blobs with no noise, flags far outliers as `-1`, and — unlike
k-means — recovers two non-convex concentric rings as distinct clusters. Tune `eps`
(the neighbourhood radius) to the data scale: too small labels everything noise, too
large merges everything into one cluster.

## Hierarchical clustering

Build a bottom-up cluster tree with single, complete, or average linkage, then cut
it into flat clusters — no need to fix `k` in advance.

```python
from quantforge import linkage, fcluster

merges = linkage(points, method="average")      # (a, b, distance, size) per merge
fcluster(points, merges, n_clusters=3)           # cut into 3 clusters
fcluster(points, merges, distance_threshold=5)   # or cut below a distance
```

Merge distances are monotone; cutting into the right number of clusters recovers
separated groups purely. Complete linkage yields compact clusters, single linkage
chains, and average sits between.

## Gaussian mixture model

Soft clustering: fit a return series as a mixture of normal regimes (e.g. calm vs
turbulent) by expectation-maximization, recovering each regime's weight, mean and
variance:

```python
from quantforge import fit_gaussian_mixture

r = fit_gaussian_mixture(returns, k=2)
r["weights"]        # regime probabilities (sum to 1)
r["means"]          # regime mean returns
r["variances"]      # regime variances
r["log_likelihood"]
```

Unlike k-means' hard assignment, each point carries a responsibility to every
component; the log-likelihood rises monotonically to a local optimum, and a
single component reduces to the plain sample mean and variance.

## Hidden Markov model

When the regimes are *sequential* — hidden states with a transition matrix — a
discrete HMM scores an observation sequence and decodes the most-likely state path:

```python
from quantforge import hmm_forward, hmm_viterbi

pi = [0.6, 0.4]                       # initial state distribution
A = [[0.7, 0.3], [0.4, 0.6]]          # state transition matrix
B = [[0.5, 0.5], [0.1, 0.9]]          # emission matrix (state -> symbol)
obs = [0, 1, 1, 0]

hmm_forward(pi, A, B, obs)            # log P(obs | model)
path, logp = hmm_viterbi(pi, A, B, obs)   # most-likely hidden-state sequence
```

`hmm_forward` uses the scaled forward algorithm, so the log-likelihood stays finite
on long sequences; `hmm_viterbi` runs in log space and returns the optimal path with
its log-probability — the standard tools for regime inference where the state
carries over in time.

For the per-time state probabilities (not just the single best path), `hmm_posterior`
runs the forward-backward algorithm:

```python
from quantforge import hmm_posterior

gamma = hmm_posterior(pi, A, B, obs)   # gamma[t][i] = P(state_t = i | obs)
```

Each `gamma[t]` is a distribution over states summing to one, using the *entire*
observation sequence — the smoothed marginal, softer than Viterbi's hard assignment.

When the parameters are unknown, `hmm_baum_welch` estimates them from the
observations by EM, and `hmm_simulate` draws sequences from a model:

```python
from quantforge import hmm_simulate, hmm_baum_welch

states, obs = hmm_simulate(pi, A, B, length=2000, seed=1)
fit = hmm_baum_welch(obs, n_states=2, n_symbols=2)
fit["pi"], fit["A"], fit["B"], fit["log_likelihood"]
```

Baum-Welch iterates forward-backward to a local optimum with a monotone
log-likelihood, recovering the transition and emission matrices up to a relabelling
of the hidden states.

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

## Market stress (turbulence / absorption ratio)

Two Kritzman-Li systemic-risk gauges. `turbulence` is the Mahalanobis distance of a
return vector from its history — it spikes on large, *unusual* cross-asset moves
that a volatility reading misses. `absorption_ratio` is the variance share of the
top principal components — high when risk is concentrated in a few factors:

```python
from quantforge import turbulence_series, absorption_ratio

turbulence_series(returns_panel)     # one turbulence value per period (mean ~ n_assets)
absorption_ratio(cov, n_factors=2)   # fraction of variance in the top 2 PCs, in [0, 1]
```

A turbulence spike marks a stressed regime; a rising absorption ratio marks a
fragile, tightly-coupled market that has historically preceded drawdowns.

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

Householder QR factors a (tall) matrix as `A = Q R` and solves least squares without
the ill-conditioned normal equations:

```python
from quantforge import qr_decomposition, qr_solve

Q, R = qr_decomposition(A)        # Q orthogonal, R upper triangular
x = qr_solve(A, b)                # least-squares solution, stabler than (A'A)^-1 A'b
```

`qr_solve` matches the OLS coefficients but stays accurate where the normal-equations
matrix `A'A` would be near-singular — the preferred path for near-collinear designs.

The singular value decomposition factors `A = U S V'` and gives the Moore-Penrose
pseudo-inverse for rank-deficient or minimum-norm least squares:

```python
from quantforge import svd, pseudo_inverse

U, s, V = svd(A)                  # s: descending singular values
Aplus = pseudo_inverse(A)         # V S^+ U', truncating tiny singular values
```

The singular values are the square roots of the eigenvalues of `A'A`; ranking them
exposes the numerical rank, and `pseudo_inverse` drops the near-zero ones so the
solve stays stable even when `A` is rank-deficient.

For a square system, LU with partial pivoting gives the solve and the determinant:

```python
from quantforge import lu_decomposition, lu_solve, determinant

lu_solve([[2, 1], [1, 3]], [5, 10])       # solve A x = b
determinant([[6, 1, 1], [4, -2, 5], [2, 8, 7]])   # -306, signed product of pivots
```

`lu_decomposition` returns `P A = L U`; `lu_solve` forward/back-substitutes, and
`determinant` is the signed product of the U pivots (zero for a singular matrix).

Those solvers use floating point; for an integer or rational matrix the answer can be
computed with *no round-off at all*. Bareiss's fraction-free elimination keeps every
intermediate an exact integer, and exact `Fraction` arithmetic then solves and inverts:

```python
from quantforge import bareiss_determinant, rational_solve, rational_inverse

bareiss_determinant([[1, 2], [3, 4]])     # -2, exact (no float error)
rational_solve([[2, 1], [1, 3]], [3, 5])  # [Fraction(4,5), Fraction(7,5)]
rational_inverse([[2, 1, 1], [1, 3, 2], [1, 0, 0]])  # exact Fraction inverse
```

`bareiss_determinant` returns an exact integer for integer input — each division in the
elimination is provably exact, so nothing rounds. `rational_solve` and `rational_inverse`
return `Fraction` results whose residual `A x - b` (or `A A⁻¹ - I`) is *exactly* zero. On
a Hilbert matrix — the textbook ill-conditioned case where a float solve loses most of its
digits — the exact solver recovers the true answer (all ones) precisely.

Direct solves cost `O(n³)`; for large or structured systems the iterative solvers only
need matrix-vector products. `conjugate_gradient` is the method of choice for a
symmetric positive-definite `A`, and `gauss_seidel` / `jacobi` handle diagonally
dominant systems:

```python
from quantforge import conjugate_gradient, gauss_seidel, jacobi

A = [[4, 1, 0], [1, 4, 1], [0, 1, 4]]
b = [6, 12, 6]
conjugate_gradient(A, b)["x"]     # [0.8571, 2.5714, 0.8571]
gauss_seidel(A, b)["n_iter"]      # 12
jacobi(A, b)["n_iter"]            # 25 — about twice as many as Gauss-Seidel
```

`conjugate_gradient` converges in at most `n` iterations for an SPD matrix (often far
fewer) and returns the residual norm; `gauss_seidel` reuses each freshly-updated
component within the sweep, so it typically converges in about half the iterations of
`jacobi`. All three return `x`, `residual_norm` and `n_iter`, and agree with the direct
`lu_solve` to machine precision on a well-conditioned system.

A *tridiagonal* system — nonzero only on the diagonal and its two neighbours, as in cubic
splines and implicit PDE steps — solves in `O(n)` rather than `O(n^3)` with the Thomas
algorithm, and `solve_cyclic_tridiagonal` handles the periodic-boundary variant:

```python
from quantforge import solve_tridiagonal, solve_cyclic_tridiagonal

# lower (sub-diag), diag, upper (super-diag), rhs — all length n
solve_tridiagonal([0, 1], [2, 2], [1, 0], [3, 3])          # [1.0, 1.0]
solve_cyclic_tridiagonal([1, 1, 1], [2, 2, 2], [1, 1, 1], [4, 4, 4])  # [1, 1, 1]
```

`solve_tridiagonal` runs the Thomas algorithm — forward elimination then back-substitution
in a single ``O(n)`` pass — where `lower[i]`/`upper[i]` are the sub/super-diagonal entries
of row `i`. `solve_cyclic_tridiagonal` adds the corner couplings of a periodic system
(`lower[0]` links row 0 to the last column, `upper[n-1]` the reverse) and solves it with a
Sherman-Morrison correction over two Thomas passes. Both agree with a dense Gaussian solve
to machine precision across thousands of random systems.

When only *one* eigenpair is needed, power iteration is far cheaper than the full
`jacobi_eigen` spectrum: `power_iteration` finds the dominant (largest-magnitude)
eigenvalue and `inverse_iteration` the one nearest a shift `mu`:

```python
from quantforge import power_iteration, inverse_iteration, rayleigh_quotient

A = [[4, 1, 0], [1, 3, 1], [0, 1, 2]]      # eigenvalues 1.2679, 3, 4.7321
power_iteration(A)["eigenvalue"]           # 4.7321 (dominant)
inverse_iteration(A, mu=0.0)["eigenvalue"] # 1.2679 (closest to 0 = smallest magnitude)
rayleigh_quotient(A, [0.79, 0.58, 0.21])   # eigenvalue estimate for a given vector
```

`power_iteration` repeatedly multiplies by `A` and renormalizes; `inverse_iteration`
runs the same loop on `(A - mu I)^{-1}`, so a shift near a known approximate eigenvalue
refines it (and `mu = 0` targets the smallest magnitude). Both return `eigenvalue`,
`eigenvector` (unit norm), `n_iter` and `converged`, and match the corresponding
`jacobi_eigen` value to `1e-5`. The Rayleigh quotient `x'Ax / x'x` is the best
eigenvalue estimate for any vector and is exact on a true eigenvector.

`jacobi_eigen` and power iteration assume a *symmetric* matrix (real spectrum); a general
real matrix can have complex eigenvalues, which `eigenvalues_general` recovers via the
characteristic polynomial:

```python
from quantforge import (characteristic_polynomial, eigenvalues_general,
                        determinant_from_charpoly)

characteristic_polynomial([[2, 0], [0, 3]])   # [1, -5, 6] = x^2 - 5x + 6
eigenvalues_general([[4, 1], [2, 3]])          # [2.0, 5.0]
eigenvalues_general([[0, -1], [1, 0]])         # [-1j, 1j] — a rotation's complex pair
```

`characteristic_polynomial` uses the Faddeev-LeVerrier recurrence to build `det(xI - A)`
from traces of matrix powers, and `eigenvalues_general` roots it with Durand-Kerner, so
complex-conjugate pairs come back correctly (real eigenvalues as `float`, complex ones as
`complex`). The eigenvalues sum to the trace and multiply to the determinant;
`determinant_from_charpoly` reads that determinant straight off the constant term, a cheap
cross-check against the LU value.

SVD-based diagnostics summarize a matrix's conditioning and size:

```python
from quantforge import condition_number, matrix_rank, spectral_norm, frobenius_norm

condition_number(A)   # sigma_max / sigma_min (inf if singular)
matrix_rank(A)        # count of singular values above a relative tolerance
spectral_norm(A)      # largest singular value (operator 2-norm)
frobenius_norm(A)     # sqrt(sum a_ij^2)
```

A large `condition_number` warns that a solve will amplify errors; `matrix_rank`
detects (numerical) rank deficiency before it corrupts a fit.

`matrix_exp` computes the matrix exponential (scaling-and-squaring with a Pade
approximant) — the operation that turns a continuous-time Markov generator into a
transition matrix:

```python
from quantforge import matrix_exp

Q = [[-0.3, 0.2, 0.1], [0.1, -0.2, 0.1], [0.05, 0.05, -0.1]]   # rate generator
P = matrix_exp(Q)     # one-period transition matrix: rows sum to 1, non-negative
```

`exp(0)` is the identity, `exp(A) exp(-A) = I`, and a generator whose rows sum to
zero maps to a valid stochastic matrix — the basis for continuous-time rating
migration.

The multivariate-normal density and its log-determinant are built stably on the
Cholesky factor — the workhorse of Gaussian likelihoods:

```python
from quantforge import log_determinant, mvn_logpdf, mvn_pdf

cov = [[1.0, 0.5], [0.5, 2.0]]
log_determinant(cov)              # 2 * sum(log L_ii), no overflow
mvn_logpdf([0.3, -0.1], mean=[0, 0], cov=cov)   # log N(x; mu, Sigma)
```

`mvn_logpdf` evaluates the Mahalanobis term by a triangular solve rather than an
explicit inverse, so it stays accurate for high-dimensional or ill-conditioned
covariances; it reduces to the univariate normal in one dimension.

## Numerical utilities

Cubic interpolation (natural spline and monotone Hermite), a spline-interpolated
zero curve, and general-purpose root finders (bisection, Brent, Newton):

```python
from quantforge import (natural_cubic_spline, monotone_cubic, SplineZeroCurve,
                        brent, newton)

f = monotone_cubic([0, 1, 2, 3], [0, 0, 0, 1])   # no overshoot
curve = SplineZeroCurve([0.5, 1, 2, 5, 10], [0.02, 0.025, 0.03, 0.035, 0.04])
brent(lambda x: x * x - 2, 0, 2)                  # sqrt(2)

# Higher-order (Halley, cubic) and derivative-free (secant, superlinear) root finders.
from quantforge import halley, secant
halley(lambda x: x*x - 2, lambda x: 2*x, lambda x: 2, 1.0)   # sqrt(2), ~cubic convergence
secant(lambda x: x*x - 2, 1, 2)                              # sqrt(2), no derivative needed

# All roots of a polynomial at once (real and complex) by Durand-Kerner.
from quantforge import polynomial_roots
polynomial_roots([1, -6, 11, -6])   # (x-1)(x-2)(x-3) -> roots near 1, 2, 3
polynomial_roots([1, 0, 1])         # x^2 + 1 -> +/- 1j

# Dense polynomial algebra on coefficient lists (low-degree-first).
from quantforge import (poly_mul, poly_divmod, poly_derivative, poly_gcd, poly_eval)
poly_mul([1, 1], [-1, 1])           # (1+x)(x-1) -> [-1, 0, 1]  = x^2 - 1
poly_divmod([-1, 0, 1], [-1, 1])    # (x^2-1)/(x-1) -> ([1, 1], [0.0])
poly_derivative([0, 0, 0, 1])       # d/dx x^3 -> [0, 0, 3]
poly_gcd([-2, 5, -4, 1], poly_derivative([-2, 5, -4, 1]))  # (x-1)^2(x-2): gcd = [-1, 1] = x-1
poly_eval([1, 2, 3], 2)             # 1 + 2*2 + 3*4 = 17

# Exact polynomial interpolation through n points (Neville / Newton form).
from quantforge import neville, divided_differences, newton_polynomial
xs = [0, 1, 2, 3]
ys = [2 * x**3 - 3 * x**2 + x - 5 for x in xs]
neville(xs, ys, 1.5)                     # (-3.5, 0.75) -> (value, error estimate)
coef = divided_differences(xs, ys)       # [-5, 0.0, 3.0, 2.0]; leading term = leading coeff
newton_polynomial(xs, coef, 1.5)         # -3.5, evaluate the Newton form anywhere
# Neville doubles as Richardson extrapolation: interpolate a step-size sequence to h=0
import math
cd = lambda h: (math.sin(1 + h) - math.sin(1 - h)) / (2 * h)
hs = [0.4, 0.2, 0.1, 0.05]
neville([h * h for h in hs], [cd(h) for h in hs], 0.0)[0]   # 0.5403023059 ~ cos(1)

# Bezier curves (de Casteljau): parametric curves from control points.
from quantforge import bezier_point, bezier_curve, bezier_tangent
bezier_point([0, 1, 0], 0.5)                       # 0.5 — scalar quadratic at t=0.5
bezier_point([(0, 0), (1, 2), (2, 0)], 0.5)        # (1.0, 1.0) — a 2-D curve
bezier_tangent([(0, 0), (1, 2), (2, 0)], 0.5)      # (2.0, 0.0) — derivative vector
bezier_curve([0, 10], 5)                           # [0, 2.5, 5, 7.5, 10] — sample a line

# Chebyshev approximation: spectral accuracy for smooth functions.
from quantforge import chebyshev_fit, chebyshev_eval, chebyshev_derivative
import math
c = chebyshev_fit(math.exp, 0, 2, degree=15)
chebyshev_eval(c, 0, 2, 1.3)                       # ~ exp(1.3) to machine precision
chebyshev_derivative(c, 0, 2)                      # coefficients of the derivative

# Definite-integral quadrature.
from quantforge import (simpson, gauss_legendre, adaptive_simpson, romberg,
                        clenshaw_curtis, tanh_sinh, gauss_kronrod)
simpson(lambda x: x * x, 0, 1)                    # composite Simpson
gauss_legendre(lambda x: x ** 5, 0, 1, n=3)       # exact to degree 2n-1
adaptive_simpson(math.sin, 0, math.pi)            # error-controlled -> 2.0
romberg(math.sin, 0, math.pi)                     # Richardson extrapolation -> 2.0
clenshaw_curtis(math.sin, 0, math.pi, n=64)       # Chebyshev points -> 2.0
gauss_kronrod(lambda x: 1/((x-0.3)**2+1e-3), 0, 1)  # adaptive, resolves sharp peaks
tanh_sinh(lambda x: 1 / math.sqrt(x), 0, 1)       # endpoint singularity -> 2.0

`romberg` (Richardson extrapolation on the trapezoid rule) and `clenshaw_curtis`
(sampling at Chebyshev points) both reach machine precision on smooth integrands;
Clenshaw-Curtis takes a free order `n`, scaling past the fixed 2-5 point
Gauss-Legendre rule and handling awkward integrands like the Runge function.
`gauss_kronrod` is the general-purpose adaptive choice — a G7-K15 rule with an
embedded error estimate that subdivides where the integrand is hardest, so it
resolves sharp peaks a fixed rule would smear over.
`tanh_sinh` (double-exponential) is the one to reach for when the integrand blows up
at an endpoint — `1/sqrt(x)`, `ln x`, `sqrt(1-x^2)` at `x = ±1` — where Simpson and
Gauss-Legendre lose accuracy; it evaluates strictly inside the interval and
converges on the singular cases to machine precision.

`polynomial_roots` finds every root at once; the `poly_*` helpers are the surrounding
algebra — `poly_mul` / `poly_divmod` (long division returning quotient and remainder),
`poly_derivative` / `poly_integral`, `poly_eval` (Horner), and `poly_gcd`, whose
`gcd(p, p')` recovers the repeated-root factor. Coefficients are low-degree-first, and
`poly_mul` is the exact direct product (use the FFT-based `convolve` for long
polynomials).

`bezier_point` evaluates a Bezier curve of any degree from its control points by de
Casteljau's repeated-interpolation scheme — numerically stable and working on scalars or
n-D coordinate tuples alike. `bezier_curve` samples it at evenly spaced parameters,
`bezier_tangent` returns the derivative vector (the curve's degree-`n-1` control differences
via `bezier_derivative_control`), and `bezier_subdivide` splits a curve at `t` into two
exact sub-curves. The curve interpolates its first and last control points and, at
`t = 0.5` of `[0, 1, 0]`, returns the Bernstein value `0.5`; `bernstein(n, i, t)` exposes the
basis. Verified against the Bernstein-basis sum over thousands of random curves.

Where a Bezier curve only touches its endpoints, a Catmull-Rom spline passes *through* every
control point, so `catmull_rom_point`/`catmull_rom_curve` interpolate a set of waypoints
smoothly:

```python
from quantforge import catmull_rom_point, catmull_rom_curve

pts = [(0, 0), (1, 1), (2, 0), (3, 1)]
catmull_rom_point(pts, 1, 0.0)   # (1.0, 1.0) — start of segment 1 == pts[1]
catmull_rom_point(pts, 1, 1.0)   # (2.0, 0.0) — end == pts[2]
catmull_rom_point(pts, 1, 0.5)   # (1.5, 0.5)
catmull_rom_curve(pts, samples_per_segment=5)   # 13 points along the whole spline
```

`catmull_rom_point(points, seg, t)` evaluates segment `seg` (between `points[seg]` and
`points[seg+1]`) at `t` in `[0, 1]`, with each segment's tangents derived from the
neighbouring points — giving a C1 curve that interpolates every input. The `alpha` argument
picks the parametrization: `0` uniform, `0.5` centripetal (the default, free of cusps and
self-intersections), `1` chordal. `catmull_rom_curve` samples the entire spline; points may
be scalars or n-D tuples. Verified over thousands of random point sets that each segment
interpolates its endpoints and adjacent segments join continuously.

Both Bezier and Catmull-Rom are special cases of the *B-spline*, whose knot vector gives
*local* control — moving one control point changes only a few spans. `bspline_point` /
`bspline_curve` evaluate it, `bspline_basis` is the Cox-de Boor basis, and
`open_uniform_knots` builds a clamped knot vector so the curve touches its endpoints:

```python
from quantforge import bspline_point, bspline_curve, open_uniform_knots

c = [(0, 0), (1, 2), (3, 3), (4, 0)]
bspline_point(c, 3, 0.0)              # (0.0, 0.0) — clamped, interpolates the first point
bspline_point(c, 3, 1.0)              # (4.0, 0.0) — and the last
bspline_point(c, 3, 0.5)              # (2.0, 1.875)
open_uniform_knots(4, 3)              # [0,0,0,0, 1,1,1,1] — no interior knots -> a Bezier
```

`bspline_point(control, degree, t)` blends the control points by the degree-`p` Cox-de Boor
basis over a knot vector (a clamped open-uniform one by default, so the endpoints are
interpolated). With `n = degree + 1` control points and no interior knots the B-spline *is*
the Bezier curve; a degree-1 spline is the piecewise-linear interpolant; and interior knots
add local, low-degree control that Bezier's global blend lacks. `bspline_basis` is a
partition of unity, and `bspline_curve` samples the whole curve. Verified over thousands of
random cases against partition-of-unity, endpoint interpolation, and Bezier equivalence.

`halley` and `secant` extend the scalar root finders past Newton: `halley` folds in the
second derivative for cubic convergence (roughly tripling the correct digits per step), and
`secant` needs no derivative at all, converging superlinearly (order ~1.618) from two
guesses. Both reject a step that stalls away from an actual root, so a bad start raises
rather than returning a non-root. Verified against known roots and, over thousands of random
cubics, agreement with a bracketed `brent`.

Where `bisection`/`brent`/`newton` solve one scalar equation, `newton_system` and
`broyden` solve a *system* `F(x) = 0` in several unknowns:

```python
from quantforge import newton_system, broyden

# circle x^2 + y^2 = 1 intersected with the line x = y
f = lambda v: [v[0]**2 + v[1]**2 - 1, v[0] - v[1]]
newton_system(f, [0.5, 0.9])         # ([0.7071, 0.7071], iters) -> 1/sqrt(2)
broyden(f, [0.5, 0.9])               # same root, quasi-Newton
```

`newton_system` forms the finite-difference Jacobian and solves `J dx = -F(x)` each step
(quadratic convergence, `n` extra evaluations per iteration for the Jacobian), while
`broyden` keeps a rank-1-updated inverse-Jacobian estimate that avoids re-differentiating
— cheaper when `F` is expensive. Both return `(solution, iterations)`, agree on the same
root, and raise on a singular Jacobian; on a 3-variable polynomial system they recover
`(1, 2, 3)` to machine precision.

A single bracketing solve needs a sign change to start; `find_all_roots` sweeps an
interval for *every* sign-changing root:

```python
from quantforge import find_all_roots, count_sign_changes
import math

find_all_roots(lambda x: (x-1)*(x-2)*(x-3), 0, 4)   # [1.0, 2.0, 3.0]
find_all_roots(math.sin, 0.5, 10)                    # [pi, 2pi, 3pi]
count_sign_changes(math.cos, 0, 20)                  # 6
```

`find_all_roots` splits `[a, b]` into `n` subintervals, brackets each sign change, and
refines it with Brent — so it recovers all simple roots of a polynomial or transcendental
equation (e.g. both roots of `e^x = 3x`) to machine precision. `count_sign_changes` is the
cheap count-only version. Even-multiplicity roots (no sign change) are missed; increase
`n` for roots closer than the grid spacing.

For a *polynomial* with real coefficients, Sturm sequences count and isolate the real roots
*exactly* — no root-finding, no missed roots — from the sign changes of a remainder chain:

```python
from quantforge import real_root_count, isolate_real_roots

p = [1, -6, 11, -6]                 # (x-1)(x-2)(x-3), highest-degree first
real_root_count(p, 0, 4)            # 3 — distinct real roots in (0, 4]
real_root_count(p, 1.5, 2.5)        # 1
isolate_real_roots(p, 0, 4)         # [(0, 1.0), (1.0, 2.0), (2.0, 4)] — one root each
real_root_count([1, 0, 1], -10, 10) # 0 — x^2 + 1 has no real roots
```

`real_root_count(coeffs, a, b)` returns the exact number of distinct real roots in the
half-open interval `(a, b]` — Sturm's theorem: the drop in sign changes of the Sturm
sequence between the endpoints — counting a repeated root once and needing no floating
root-finder. `isolate_real_roots` bisects on that count to hand back disjoint intervals each
bracketing a single root, ready for a Brent refine. `sturm_sequence` exposes the underlying
chain. Coefficients are highest-degree first, as with `polynomial_roots`. Verified against
polynomials built from known roots over thousands of cases.

Iterating any map `x -> f(x)` eventually repeats, and `floyd_cycle`/`brent_cycle` find the
tail length and cycle length of that "rho" in `O(1)` memory — no set of visited states:

```python
from quantforge import floyd_cycle, brent_cycle, cycle_elements

nxt = [1, 2, 3, 4, 2]        # 0->1->2->3->4->2: tail 0,1 then cycle 2,3,4
f = lambda x: nxt[x]
floyd_cycle(f, 0)            # (2, 3) — (mu tail length, lam cycle length)
cycle_elements(f, 0)        # [2, 3, 4]

lcg = lambda x: (13 * x + 7) % 1000
floyd_cycle(lcg, 0)         # (0, 200) — the LCG's period from 0
```

Both return `(mu, lam)`: `mu` is how many steps before the sequence enters its cycle and
`lam` is the cycle length. `floyd_cycle` is the classic tortoise-and-hare; `brent_cycle`
finds the same answer with typically fewer function calls. `cycle_elements` lists the cycle
itself. Because they use constant memory, they scale to enormous state spaces — the reason
Pollard's rho factoring and PRNG-period analysis rely on them. Cross-checked against a
seen-set walk over thousands of random functional graphs.

When the answer is an *integer* characterized by a monotone or unimodal property — the
"binary search the answer" pattern — `first_true`/`last_true` and the integer ternary
searches locate it in `O(log n)`:

```python
from quantforge import first_true, last_true, ternary_search_int_max, ternary_search_int_min

first_true(0, 100, lambda x: x * x >= 50)          # 8 — smallest x with x^2 >= 50
last_true(0, 100, lambda x: x * x <= 50)           # 7 — largest x with x^2 <= 50
ternary_search_int_max(0, 20, lambda x: -(x-7)**2) # 7 — argmax of a unimodal f
ternary_search_int_min(0, 20, lambda x: (x-13)**2) # 13 — argmin
```

`first_true` returns the smallest integer in `[lo, hi]` where a monotone predicate first
holds (`hi + 1` if never), `last_true` the largest where a decreasing predicate still holds
(`lo - 1` if never) — the workhorses for feasibility questions like "smallest capacity that
fits the schedule". `ternary_search_int_max`/`_min` optimize a strictly unimodal integer
function. All four are checked against brute linear scans over thousands of random monotone
thresholds and unimodal functions.

To evaluate an arithmetic string *without* Python's `eval` — so untrusted input is safe —
`eval_expression` runs the shunting-yard algorithm:

```python
from quantforge import eval_expression

eval_expression("2 + 3 * 4")        # 14.0 — precedence respected
eval_expression("2^3^2")            # 512.0 — power is right-associative
eval_expression("-2^2")             # -4.0 — power binds tighter than unary minus
eval_expression("sqrt(3^2 + 4^2)")  # 5.0
eval_expression("2 * pi")           # 6.283185...
```

`tokenize` splits the string, `shunting_yard` converts it to Reverse Polish honouring
precedence and associativity, and `eval_rpn` folds the postfix stream to a number;
`eval_expression` chains all three. It supports `+ - * / // % ^`, unary minus, parentheses,
the functions `sin/cos/tan/exp/log/sqrt/abs`, and the constants `pi`/`e`. Checked against
Python's `eval` over thousands of random expressions.

`neville` evaluates the unique degree-`(n-1)` polynomial through `n` points at one
`x`, returning `(value, error_estimate)`; `divided_differences` / `newton_polynomial`
build the Newton form once and evaluate it cheaply at many points. Because Neville
extrapolates a tabulated sequence to any target, feeding it a step-size sequence and
`x = 0` performs Richardson extrapolation — the derivative example above lands on
`cos(1)` to 13 digits, the same idea `romberg` uses on the trapezoid rule.

For repeated evaluation, the barycentric form precomputes weights once and then
evaluates in `O(n)` per point, stably — and on Chebyshev nodes it converges spectrally
instead of oscillating (the Runge phenomenon):

```python
from quantforge import (barycentric_weights, barycentric_eval,
                        chebyshev_nodes, chebyshev_barycentric_weights)

runge = lambda x: 1 / (1 + 25 * x * x)
grid = [i / 100 for i in range(-99, 100)]

cn = chebyshev_nodes(-1, 1, 21)
cw = barycentric_weights(cn)
cy = [runge(x) for x in cn]
max(abs(barycentric_eval(cn, cy, cw, x) - runge(x)) for x in grid)   # ~0.0177

en = [-1 + 2 * i / 20 for i in range(21)]           # 21 equispaced nodes
ew = barycentric_weights(en)
ey = [runge(x) for x in en]
max(abs(barycentric_eval(en, ey, ew, x) - runge(x)) for x in grid)   # ~58.6 — blows up
```

`barycentric_weights` costs `O(n²)` once; each `barycentric_eval` is `O(n)`, so
changing the sampled `ys` (same nodes) is free. `chebyshev_nodes` gives the
Chebyshev-Lobatto points and `chebyshev_barycentric_weights` their closed-form
weights; interpolating a smooth function like `exp` on 25 Chebyshev nodes is accurate
to machine precision, where 21 equispaced nodes on the Runge function are off by ~59.

The 1-D interpolators above have a 2-D cousin for lookup tables — vol surfaces, response
grids, heightmaps — sampled on a rectilinear grid:

```python
from quantforge import bilinear_interp, nearest_interp

xs, ys = [0, 1, 2, 3], [0, 1, 2]
z = [[2*x + 3*y + 1 for x in xs] for y in ys]   # a planar field on the grid
bilinear_interp(xs, ys, z, 1.5, 0.5)             # 5.5 — exact for a planar field
nearest_interp(xs, ys, z, 0.6, 0.6)              # snaps to the closest node
```

`bilinear_interp` blends the four surrounding grid values — exact at the nodes, exact for
any planar field, and equal to the four-corner mean at a cell center — while
`nearest_interp` snaps to the closest node (piecewise constant, preserving the sampled
values). Both clamp out-of-range queries to the grid edge; `z` is indexed `z[row][col]`
with strictly increasing axes.

All of the above fit a *polynomial*. When the data has poles or a flattening tail,
Thiele's method fits a *rational* function (a ratio of polynomials) as a continued
fraction and still passes through every node:

```python
import math
from quantforge import thiele_interpolate

f = lambda x: (2 * x + 1) / (x * x + 1)
xs = [0, 1, 2, 3, 4, -1]
ys = [f(x) for x in xs]
thiele_interpolate(xs, ys, 2.7)     # 0.77201448 — recovers the rational exactly

xs = [0.1 * i for i in range(1, 10)]
thiele_interpolate(xs, [math.tan(x) for x in xs], 0.55)   # 0.61310521 ~ tan(0.55)
```

Thiele builds the reciprocal-difference coefficients (`thiele_coefficients`) once and
`thiele_eval` evaluates the continued fraction anywhere. Because it is rational, six
nodes reproduce a degree-2-over-degree-2 rational to machine precision, and it captures
the pole structure of `tan` that a polynomial of the same node count would badly miss.

For expectations under a normal density, Gauss-Hermite quadrature is exact for
polynomials up to degree `2n-1` and needs only a handful of nodes:

```python
from quantforge import gauss_hermite_expectation

gauss_hermite_expectation(lambda x: x * x, mu=1.5, sigma=0.7)   # E[X^2] = mu^2+sigma^2
gauss_hermite_expectation(math.exp, mu=0.05, sigma=0.3)         # lognormal mean
```

The nodes and weights (probabilists' convention, weights summing to one) come from
the Golub-Welsch eigen-decomposition of the Hermite recurrence; `gauss_hermite_nodes_weights(n)`
exposes them directly. Best for smooth integrands — a kinked payoff converges
slowly, so integrate option payoffs with the density routines instead.

For integrals over the half-line `[0, inf)`, Gauss-Laguerre uses the `e^{-x}`
weight; the wrapper factors out an exponential rate so it also handles a plain
`integral_0^inf g(x) dx`:

```python
from quantforge import gauss_laguerre_integral

gauss_laguerre_integral(lambda x: math.exp(-2 * x), rate=2.0)          # -> 0.5
gauss_laguerre_integral(lambda x: x * x * math.exp(-3 * x), rate=3.0)  # -> 2/27
```

Set `rate` near the integrand's true decay constant for best accuracy; it needs an
exponential falloff, so an integrand like `1/(1+x^2)` (only algebraic decay) will
not converge.

For a double integral over a rectangle, `integrate2d_gauss` is the tensor
Gauss-Legendre rule and `integrate2d_simpson` the composite-Simpson fallback:

```python
from quantforge import integrate2d_gauss, integrate2d_simpson

integrate2d_gauss(lambda x, y: x * y, 0, 1, 0, 1)        # 0.25
integrate2d_gauss(lambda x, y: x + y, 0, 1, 0, 1, n=2)   # 1.0
integrate2d_simpson(lambda x, y: math.exp(-(x*x + y*y)), -6, 6, -6, 6, nx=200, ny=200)  # ~pi
```

`integrate2d_gauss` uses `n^2` evaluations and is exact for polynomials up to degree
`2n-1` in each variable — the right tool for smooth integrands. `integrate2d_simpson`
integrates in `y` for each `x` and then in `x` (Fubini), converging on the wider or
less-smooth cases where the fixed Gauss rule is too coarse (it recovers the 2-D Gaussian
integral `pi` on a large box).

To recover a time function from a Laplace transform known only as a formula,
`laplace_inversion` uses the real-arithmetic Gaver-Stehfest algorithm:

```python
from quantforge import laplace_inversion

laplace_inversion(lambda s: 1.0 / (s * s), t=2.0)     # 2.0 (transform of f(t) = t)
laplace_inversion(lambda s: 1.0 / (s + 2.0), t=1.0)   # 0.1354 ~ exp(-2)
```

It evaluates `F(s)` at real points `k ln2 / t` and combines them with tabulated
coefficients, so no complex contour is needed. It recovers constants, polynomials,
`sqrt(t)`, and decaying exponentials to ~1e-3 or better; like all Gaver-Stehfest
implementations it struggles with oscillatory or fast-growing functions (limited by
double precision), so keep the target smooth and non-oscillatory. Useful for inverting
transform-domain formulas in queueing, diffusion, and arithmetic-Asian / occupation-time
pricing.

# One-dimensional minimizers (line search / 1-D calibration).
from quantforge import golden_section_min, brent_min
golden_section_min(lambda x: (x - 3) ** 2 + 1, -10, 10)   # -> (3.0, 1.0)
brent_min(math.cos, 0, 2 * math.pi)                        # -> (pi, -1.0)

# Multivariate minimizers: Nelder-Mead (local) and differential evolution (global).
from quantforge import nelder_mead, differential_evolution
nelder_mead(lambda v: (v[0] - 3) ** 2 + (v[1] + 1) ** 2, [0.0, 0.0])   # -> ([3.0, -1.0], 0.0)
rosen = lambda v: (1 - v[0]) ** 2 + 100 * (v[1] - v[0] ** 2) ** 2
differential_evolution(rosen, [(-5, 5), (-5, 5)], max_iter=2000)["x"]  # -> [1.0, 1.0]

# Central-difference gradient / Hessian / Jacobian (Greeks, calibration).
from quantforge import gradient, hessian, jacobian
gradient(lambda v: v[0] ** 2 + 3 * v[1] ** 2, [1.0, 2.0])   # -> [2.0, 12.0]
hessian(lambda v: v[0] ** 2 + v[0] * v[1], [1.0, 1.0])       # -> [[2, 1], [1, 0]]

# High-accuracy 1-D derivatives (Ridders' adaptive Richardson extrapolation).
from quantforge import ridders_derivative, ridders_second_derivative, price, delta
ridders_derivative(math.sin, 0.7)             # (0.764842187284, ~2e-14 error estimate)
ridders_second_derivative(math.exp, 1.3)      # (3.6692966676, err) == exp(1.3)
# numerical Black-Scholes delta matches the closed form to 8+ digits
ridders_derivative(lambda S: price(S, 100, 1.0, 0.05, 0.2), 100.0)[0]  # 0.63683065
delta(100, 100, 1.0, 0.05, 0.2)                                        # 0.63683065

# Nonlinear least-squares calibration (Levenberg-Marquardt, finite-diff Jacobian).
from quantforge import levenberg_marquardt
xs = [i * 0.2 for i in range(30)]
ys = [2.0 * math.exp(0.5 * x) for x in xs]
fit = levenberg_marquardt(lambda b, x: b[0] * math.exp(b[1] * x), xs, ys, [1.0, 0.1])
fit["parameters"]        # -> [2.0, 0.5]  (recovered exactly)
```

`nelder_mead` is the derivative-free local simplex minimizer (the engine behind several
of the library's calibrations); `differential_evolution` is its global counterpart,
evolving a population to escape local minima on multimodal, non-convex objectives where
Nelder-Mead would stall. On the multimodal Rastrigin function `differential_evolution`
reaches the global optimum while `nelder_mead` from a poor start settles into a nearby
local pit; it is deterministic for a fixed `seed` and confines the search to the given
box `bounds`. Use Nelder-Mead when you have a good starting guess and a smooth basin,
differential evolution when the landscape is rough or the starting region is unknown.

`simulated_annealing` is a third option — a single point that occasionally accepts a
*worse* move (probability `exp(-delta / T)`) so it can climb out of a local minimum, with
the temperature `T` cooling geometrically toward a settled basin:

```python
from quantforge import simulated_annealing

rosen = lambda v: (1 - v[0]) ** 2 + 100 * (v[1] - v[0] ** 2) ** 2
simulated_annealing(rosen, [-3.0, 3.0], bounds=[(-5, 5), (-5, 5)],
                    T0=2.0, max_iter=50000)["x"]     # ~ [1.0, 1.0]
```

It shares the same reproducibility (`seed`) and box-`bounds` behaviour as
`differential_evolution` but keeps no population, so it is cheaper per iteration and a
good fit when each objective evaluation is expensive. Differential evolution is usually
the more reliable global search; simulated annealing shines on very high-dimensional or
combinatorial-flavoured landscapes where maintaining a population is costly.

When the objective is *smooth*, a gradient method converges far faster than any of the
above. `bfgs` is a quasi-Newton minimizer that builds an inverse-Hessian approximation
from successive numerical gradients — near-Newton speed using only the objective:

```python
from quantforge import bfgs

bfgs(lambda v: (v[0] - 3) ** 2 + (v[1] + 1) ** 2, [0.0, 0.0])["x"]   # [3.0, -1.0], 2 iters
bfgs(lambda v: (1 - v[0]) ** 2 + 100 * (v[1] - v[0] ** 2) ** 2, [-1.2, 1.0])["x"]  # [1.0, 1.0]
```

It reaches a quadratic's minimum in a couple of iterations (superlinear convergence)
and the Rosenbrock minimum in a few dozen, where the derivative-free methods take
thousands of evaluations. An Armijo backtracking line search keeps every step a descent,
and it resets to steepest descent if the Hessian approximation loses positive-
definiteness. Use `bfgs` for smooth local refinement — often as a fast polish after a
global method (`differential_evolution` or `simulated_annealing`) has located the right
basin.

`levenberg_marquardt` needs only the model function -- the residual Jacobian is
taken numerically -- and converges from a poor starting guess, making it the
general calibration engine (vol surface, curve, or any parametric fit).

When the coefficients must be non-negative -- weights, proportions, spectral mixing --
`nnls` solves `min ||A x - b||²` subject to `x >= 0` exactly by the Lawson-Hanson
active-set method:

```python
from quantforge import nnls

nnls([[1, 0], [0, 1], [1, 1]], [2, 3, 5])["x"]     # [2.0, 3.0] — exact non-negative fit
nnls([[1.0], [1.0]], [-2.0, -3.0])["x"]            # [0.0] — the negative OLS root is clamped
r = nnls([[1, 1], [1, 2], [1, 3]], [1, 2, 2])
r["x"], r["residual_norm"]                          # [0.6667, 0.5], 0.4082
```

Where the unconstrained least-squares solution is already non-negative, `nnls`
reproduces it; otherwise it returns the best feasible fit, often with several
coefficients exactly zero (a natural sparsity, no penalty needed). The result satisfies
the Karush-Kuhn-Tucker conditions — the objective gradient is zero along the active
(positive) coefficients and non-positive along the zeroed ones — so it is the global
constrained optimum, not a heuristic.

For a *linear* objective under linear constraints, `linprog` solves the LP exactly by
two-phase simplex:

```python
from quantforge import linprog

# maximize 3x + 2y  s.t.  x + y <= 4,  x + 3y <= 6,  x, y >= 0
linprog([3, 2], [([1, 1], "<=", 4), ([1, 3], "<=", 6)])   # x=[4,0], objective 12

# minimize 2x + 3y  s.t.  x + y >= 10
linprog([2, 3], [([1, 1], ">=", 10)], maximize=False)     # objective 20
```

Constraints are `(coeffs, op, rhs)` with `op` in `"<="`, `">="`, `"="`, and variables are
`>= 0`. A phase-1 artificial-variable stage finds an initial feasible basis for `>=`/`=`
rows, and Bland's rule prevents cycling; the result matches a brute-force vertex
enumeration on small problems. Unbounded and infeasible programs raise rather than
returning a wrong vertex.

Two classic *discrete* optimizers round out the section — matching agents to tasks and
filling a budget:

```python
from quantforge import hungarian, knapsack_01

# minimum-cost assignment of 3 rows to 3 columns
hungarian([[4, 1, 3], [2, 0, 5], [3, 2, 2]])   # ([1, 0, 2], 5)

# maximize value with total weight <= 5
knapsack_01([2, 3, 4, 5], [3, 4, 5, 6], 5)     # (7, [0, 1])
```

`hungarian` solves the assignment problem (minimum-cost perfect matching) in `O(n^3)` via
the Kuhn-Munkres potentials, returning the column chosen for each row plus the total cost
— matching a brute-force permutation search. Negate the cost matrix to maximize.
`knapsack_01` is the 0/1 knapsack dynamic program: it maximizes total value under a weight
budget with each item taken at most once, returning the value and the chosen item indices
(both verified against a brute-force subset search).

Three more small dynamic programs cover coin change and subset sum:

```python
from quantforge import min_coins, count_change, subset_sum

min_coins([1, 3, 4], 6)              # (2, [3, 3]) — fewest coins (greedy would misfire)
count_change([1, 2, 5], 5)           # 4 — distinct multisets summing to 5
subset_sum([3, 34, 4, 12, 5, 2], 9)  # (True, [4, 5]) — a subset that hits the target
```

`min_coins` returns the fewest coins summing to a target (unlimited supply of each) with
one witnessing multiset, `-1` if impossible — the `[1, 3, 4]` case shows why a greedy pick
fails where the DP does not. `count_change` counts the distinct *order-independent*
multisets that reach the target, and `subset_sum` decides whether a subset of a list (each
item used once) sums to the target, returning one such subset. Each runs in
`O(target · n)` and is checked against brute-force references over thousands of instances.

`ridders_derivative` / `ridders_second_derivative` return `(value, error_estimate)`:
they evaluate a central difference at a shrinking step sequence and Richardson-
extrapolate across it, stopping when round-off starts to dominate — so they reach
near machine precision without hand-tuning the step, unlike the fixed-step `gradient`.
The numerical Black-Scholes delta above landing on the closed-form value is the kind
of cross-check that catches a sign or scaling slip in a hand-coded Greek.

To *build* a finite-difference stencil rather than apply a fixed one, `fd_weights` gives
Fornberg's coefficients for any derivative order on an arbitrary grid:

```python
from quantforge import fd_weights

fd_weights(0, [-1, 0, 1], 2)   # [[0,1,0], [-0.5,0,0.5], [1,-2,1]]
fd_weights(0, [0, 1, 2], 1)[1] # [-1.5, 2.0, -0.5] — one-sided first derivative
```

`fd_weights(x0, grid, max_deriv)` returns, for each derivative order `0..max_deriv`, the
weights `c` with `f^(k)(x0) ~ sum_i c[i] f(grid[i])` — exact for polynomials up to the grid
degree. The grid need not be uniform or symmetric, so one call produces central, forward,
backward, or staggered stencils to any order; the classic `[-0.5, 0, 0.5]` and `[1, -2, 1]`
fall out for `{-1, 0, 1}`. Verified against the textbook stencils and polynomial exactness
on every monomial up to the grid degree.

When the function can be evaluated on complex inputs, the complex-step method is even
sharper — it has *no* subtractive cancellation, so the step can be made
arbitrarily small:

```python
import cmath
from quantforge import complex_step_derivative, complex_step_gradient

complex_step_derivative(cmath.sin, 0.7)              # 0.7648421872844885
complex_step_derivative(cmath.sin, 0.7, h=1e-100)    # identical — no round-off floor
complex_step_gradient(lambda v: v[0]**2 + 3*v[1]**2 + v[0]*v[1], [1.0, 2.0])  # [4.0, 13.0]
```

`complex_step_derivative` evaluates `Im(f(x + i·h)) / h`: the imaginary part carries
the derivative with an `O(h²)` error and no subtraction of nearby values, so even
`h = 1e-100` gives a machine-precision result where a real finite difference would
collapse to noise. The one requirement is that `f` be written with complex-safe
operations (pass `cmath.sin` rather than `math.sin`, and avoid `abs`/`max`/comparisons
that break on complex numbers).

`dual_derivative` goes one better with forward-mode automatic differentiation: a `Dual`
number `a + b·eps` (with `eps² = 0`) carries the value and its derivative together, so
the chain rule is exact by construction — no truncation *or* round-off, and no
complex-analyticity requirement:

```python
from quantforge import dual_derivative, Dual
from quantforge.dual import exp, sin

dual_derivative(lambda x: x * x * x, 2.0)          # 12.0 — exact, zero error
dual_derivative(lambda x: exp(x) * sin(x), 0.8)    # 3.14705464

x = Dual(3.0, 1.0)                                 # a variable seeded with derivative 1
r = x * x + 2 * x + 1                              # (x+1)^2
r.value, r.deriv                                   # 16.0, 8.0
```

Write the function with `Dual` arithmetic and the module's dual-aware
`exp`/`log`/`sqrt`/`sin`/`cos`/`tan`/`tanh`; `dual_derivative(f, x)` returns `f'(x)` to
the last bit. It even differentiates a dual exponent like `x**x` correctly. Reach for
it when you control the function's source (the exact-derivative gold standard); use
`complex_step_derivative` or `ridders_derivative` when you only have a black-box
callable.

Two helpers build on it: `dual_gradient` for the exact gradient of a scalar function of
a vector, and `dual_newton` for a Newton root solve that needs *only* the function (its
derivative comes from autodiff):

```python
from quantforge import dual_gradient, dual_newton

dual_gradient(lambda v: v[0]**2 + 3*v[1]**2 + v[0]*v[1], [1.0, 2.0])   # [4.0, 13.0]
dual_newton(lambda x: x*x - 2, 1.0)["root"]                           # 1.414213562373
```

`dual_gradient` seeds each coordinate's derivative to 1 in turn (`n` exact evaluations),
matching the analytic gradient to the last bit including transcendentals. `dual_newton`
reads `f(x)` and `f'(x)` from a single dual evaluation per step, so — unlike the library
`newton` — you never supply a hand-coded derivative; it reaches sqrt(2) to machine
precision in six iterations.

For the *second* derivative, hyperdual numbers carry the value, both first-order slots
and the cross term, so one evaluation returns `f`, `f'` and `f''` exactly:

```python
from quantforge import hyperdual_derivatives, second_derivative
from quantforge.hyperdual import exp, sin

hyperdual_derivatives(lambda x: x * x * x, 2.0)          # (8.0, 12.0, 12.0) = (f, f', f'')
second_derivative(lambda x: exp(x) * sin(x), 0.8)        # 3.10109859
```

Unlike differencing a first-difference (which loses precision to a squared step),
hyperdual second derivatives have no truncation error and agree with
`ridders_second_derivative` to machine precision — the exact way to get gamma-style
convexities when you control the pricing function's source.

Slowly-converging sequences and fixed-point iterations can be accelerated with
`aitken` (delta-squared), `shanks`, and `steffensen`:

```python
import math
from quantforge import aitken, steffensen

# partial sums of the Leibniz series for pi/4 converge painfully slowly
seq = []
acc = 0.0
for k in range(20):
    acc += (-1) ** k / (2 * k + 1)
    seq.append(acc)
seq[-1] * 4          # 3.09162  — still wrong in the 2nd decimal after 20 terms
aitken(seq)[-1] * 4  # 3.1415563 — Aitken extrapolation, three digits better

# Steffensen: quadratic fixed-point convergence with no derivative
steffensen(math.cos, 0.5)   # {'root': 0.7390851332, 'iterations': 5, 'converged': True}
```

`aitken` extrapolates each limit from three consecutive terms (`shanks` is the same
transform); `steffensen` fuses that extrapolation into a `x = g(x)` iteration to get
Newton-like quadratic convergence from a merely linear map — the cos fixed point above
takes 5 steps where plain iteration needs about 69.

For series that even Aitken barely dents, `wynn_epsilon` iterates the Shanks transform
to all orders in one table, and `euler_transform` accelerates alternating series
directly:

```python
from quantforge import wynn_epsilon, euler_transform

ps = []
acc = 0.0
for k in range(20):
    acc += (-1) ** k / (2 * k + 1)
    ps.append(acc)
wynn_epsilon(ps) * 4        # 3.141592653589792 — pi to machine precision from 20 terms

euler_transform([1 / (k + 1) for k in range(30)])   # 0.6931471805  ~ ln 2
```

`wynn_epsilon` takes the *partial sums* and returns the best converged even-column
estimate (it deliberately steps back from later columns once round-off degrades them);
`euler_transform` takes the term *magnitudes* of an alternating series and reweights
their forward differences by powers of one half. Wynn is the stronger general tool —
it drives the slowly-alternating Leibniz series to full `double` precision from just
twenty terms.

When you control the *step size*, `richardson_extrapolate` accelerates a sequence of
estimates `A(h), A(h/2), A(h/4), ...` to the `h -> 0` limit by cancelling the leading
error terms:

```python
from quantforge import richardson_extrapolate

cd = lambda h: (math.sin(1 + h) - math.sin(1 - h)) / (2 * h)   # central diff, O(h^2)
richardson_extrapolate([cd(0.2 / 2**i) for i in range(6)], p=2)   # 0.5403023058681 = cos(1)
```

`p` is the leading error exponent (1 for a forward difference, 2 for a central one), and
each tableau column cancels the next error term — this is exactly the engine behind
Romberg integration (feed it the trapezoid sequence with `p=2` and it returns the
Romberg value). `richardson_table` exposes the full triangular tableau so you can watch
the diagonal stabilize.

For initial-value ODEs `y'(t) = f(t, y)`, `rk4` is the fixed-step fourth-order
Runge-Kutta and `rk45` the adaptive Dormand-Prince method (scalar or vector systems):

```python
from quantforge import rk4, rk45

rk4(lambda t, y: y, 0, 1.0, 1.0, n=100)[1][-1][0]        # 2.71828183 = e
ts, ys = rk45(lambda t, y: [y[1], -y[0]], 0, [1.0, 0.0], 2*math.pi, tol=1e-10)
ys[-1]                                                    # [1.0, -0.0] — oscillator back to start
```

Both return `(ts, ys)` — the time points and the state (a list) at each. `rk4` uses `n`
equal steps and is `O(h^4)` accurate for smooth problems; `rk45` carries an embedded
error estimate to resize the step, taking small steps only where the solution moves
fast, and returns the accepted (non-uniform) points. Pass a vector `y0` for a system —
the harmonic oscillator above integrates `[y, y']` back to `[1, 0]` over a full period.

Where RK4 evaluates the derivative four times per step, `adams_bashforth_moulton` reuses
past derivatives — a 4th-order predictor-corrector (PECE) with just one new evaluation per
step after the corrector:

```python
from quantforge import adams_bashforth_moulton
import math

ts, ys = adams_bashforth_moulton(lambda t, y: y, 1.0, 0, 1, 100)
ys[-1]                              # 2.71828183 ~ e
adams_bashforth_moulton(lambda t, y: [y[1], -y[0]], [0.0, 1.0], 0, math.pi, 200)[1][-1]
                                    # [~0, -1.0] — sin/cos system at t = pi
```

It bootstraps the first three steps with RK4 to build the derivative history, then predicts
with the 4-step Adams-Bashforth formula and corrects with Adams-Moulton — cheaper per step
than RK4 when the derivative is expensive, at the same `O(h^4)` accuracy. Fixed step over
`n_steps`; scalar or vector `y0`. Verified against closed-form solutions and confirmed to
converge at 4th order.

When the solution is pinned at *both* ends instead of given an initial slope,
`shooting_bvp` solves the two-point boundary-value problem `y'' = f(t, y, y')`:

```python
from quantforge import shooting_bvp

# y'' = y, y(0) = 0, y(1) = 1  ->  y = sinh(t)/sinh(1)
r = shooting_bvp(lambda t, y, yp: y, 0, 1, alpha=0.0, beta=1.0, s_lo=-5, s_hi=5)
r["slope"]     # 0.85091813 = 1/sinh(1), the initial slope y'(0)
```

It guesses the initial slope, integrates to the far end with adaptive RK45, and
Brent-root-finds the slope whose terminal value hits `beta` — so `s_lo`/`s_hi` must
bracket a sign change in the boundary residual. The result carries the found `slope` and
the `(ts, ys)` solution; on the linear test problems it recovers the analytic slope and
every interior point to machine precision.

For long-horizon *conservative* systems (orbits, molecular dynamics), a general RK
method slowly drifts the total energy; `velocity_verlet` and `leapfrog` are symplectic
integrators that keep it bounded instead:

```python
from quantforge import velocity_verlet, leapfrog

# harmonic oscillator q'' = -q, mass 1, over 2000 steps
qs, ps = velocity_verlet(lambda q: -q, q0=1.0, p0=0.0, mass=1.0, dt=0.05, n_steps=2000)
# energy 0.5(p^2 + q^2) oscillates in a tiny band around 0.5 — it does not drift

# 2-D circular orbit under central gravity keeps its radius over 10k steps
def gravity(q):
    r = (q[0] ** 2 + q[1] ** 2) ** 0.5
    return [-q[0] / r ** 3, -q[1] / r ** 3]
qs, ps = velocity_verlet(gravity, [1.0, 0.0], [0.0, 1.0], mass=1.0, dt=0.001, n_steps=10000)
```

`velocity_verlet` takes the force (position → force) and evolves position and momentum;
`leapfrog` is the equivalent velocity-form for unit mass. Both are second-order and
time-reversible, and because they preserve phase-space structure the energy error
oscillates rather than accumulating — exactly what you want for a simulation run over
many periods, where RK4's small per-step energy leak would compound.

For a diffusion PDE `u_t = alpha u_xx` (heat conduction, or the transformed
Black-Scholes equation), `heat_equation_cn` marches an initial profile forward by
Crank-Nicolson:

```python
from quantforge import heat_equation_cn

# initial profile u0 on a grid of spacing dx, diffusivity alpha, fixed endpoints
u = heat_equation_cn(u0, alpha=0.5, dx=0.01, dt=0.0005, n_steps=200, left=0.0, right=0.0)
```

Crank-Nicolson averages the explicit and implicit Euler steps, so it is second-order in
both space and time and *unconditionally stable* — no CFL restriction linking `dt` to
`dx²`, so you can take large time steps without the solution blowing up. A sine initial
mode decays to its analytic `exp(-alpha (pi/L)² t)` amplitude, a step in the boundary
values relaxes to the exact linear steady state, and a localized bump spreads and
flattens — the standard diffusive behaviour, solved by a tridiagonal system each step.

For the *hyperbolic* wave equation `u_tt = c² u_xx` (vibrating string, acoustics),
`wave_equation` marches an initial displacement and velocity forward by an explicit
leapfrog scheme:

```python
from quantforge import wave_equation

# plucked string u(x,0) = shape, at rest (v0 = 0), speed c
u = wave_equation(u0, v0=[0.0] * len(u0), c=1.0, dx=0.005, dt=0.0025,
                  n_steps=200, left=0.0, right=0.0)
```

Unlike the diffusion solver this one is *conditionally* stable — the Courant number
`c·dt/dx` must stay `<= 1` (the CFL condition, enforced with a clear error). A standing
mode oscillates and returns to its start after one period `2L/c`, and a localized pulse
splits into two half-height waves travelling in opposite directions, exactly as
d'Alembert's solution predicts.

The third PDE type, the *elliptic* Poisson equation `u_xx + u_yy = f` (steady-state heat,
electrostatics), is a boundary-value problem in 2-D; `poisson2d` solves it by successive
over-relaxation:

```python
from quantforge import poisson2d

# f_grid is the source (zeros for Laplace); boundary holds the fixed edge values
u, n_iter = poisson2d(f_grid, boundary, dx=0.05, dy=0.05, omega=1.8)
```

With `f_grid` all zeros it solves Laplace's equation, whose solution is harmonic — the
maximum principle holds (no interior value exceeds the boundary range) and a linear
boundary is reproduced exactly. `omega` is the over-relaxation factor in `(0, 2)`;
tuning it toward ~1.8 converges several times faster than plain Gauss-Seidel (`omega=1`).
A `sin*sin` source recovers its analytic solution to the `O(h^2)` discretization error.

A Padé approximant turns a Taylor series into a *rational* function that often
converges where the series itself diverges. `pade` builds `[m/n]` from Taylor
coefficients and `pade_eval` evaluates it; `lentz_continued_fraction` evaluates a
general continued fraction by the modified Lentz algorithm:

```python
import math
from quantforge import pade, pade_eval, lentz_continued_fraction

coeffs = [1 / math.factorial(k) for k in range(5)]     # exp Taylor series
num, den = pade(coeffs, 2, 2)                          # [2/2] approximant of exp
pade_eval(num, den, 1.0)                               # 2.7142857 (vs e = 2.7182818)

# ln(1+x): its Taylor series diverges at x = 2, but the [4/4] Pade converges
lc = [0.0] + [(-1) ** (k + 1) / k for k in range(1, 9)]
pade_eval(*pade(lc, 4, 4), 2.0)                        # 1.09857 ~ ln 3

# golden ratio 1 + 1/(1 + 1/(1 + ...))
lentz_continued_fraction(lambda k: 1.0, lambda k: 1.0)  # 1.6180339887
```

`pade` matches the input series through order `m + n` but, being rational, models
poles and wide ranges the polynomial can't — the `ln 3` example lands within `1e-4`
where the raw Taylor series has no hope. `lentz_continued_fraction` takes callables
for the partial numerators `a(k)` and denominators `b(k)` and is the same
numerically-stable engine behind the library's incomplete-gamma and beta functions.

Where Padé approximates a *function*, a continued fraction also gives the best rational
approximation of a single *number*. Every real has an expansion whose truncations — the
convergents — are closer than any fraction with a smaller denominator:

```python
from quantforge import cf_expansion, convergents, best_rational

cf_expansion(math.pi, max_terms=5)   # [3, 7, 15, 1, 292]
convergents(cf_expansion(math.pi, 6))# [(3,1), (22,7), (333,106), (355,113)]
best_rational(math.pi, 1000)         # (355, 113) — closest with denominator <= 1000
```

The convergents of pi are the familiar `22/7`, `333/106`, `355/113`; `sqrt(2)` expands to
the repeating `[1; 2, 2, 2, ...]`. `best_rational` returns the closest fraction whose
denominator stays under a bound (`355/113` for pi at 1000, `22/7` at 50), agreeing with
Python's `Fraction.limit_denominator` across thousands of random values — useful for
recovering a simple fraction from a noisy decimal or designing a ratio to a tolerance.

The same rationals form the Stern-Brocot tree, where every positive fraction sits at a
unique left/right path from `1/1`, and the Farey sequences enumerate them in order:

```python
from quantforge import stern_brocot_path, best_rational_bounded, farey_sequence, mediant
import math

stern_brocot_path(Fraction(2, 3))     # 'LR' — the path to 2/3
best_rational_bounded(math.pi, 113)   # Fraction(355, 113)
best_rational_bounded(math.pi, 7)     # Fraction(22, 7)
farey_sequence(4)                     # [0, 1/4, 1/3, 1/2, 2/3, 3/4, 1]
mediant((1, 2), (1, 3))               # Fraction(2, 5)
```

`stern_brocot_path` gives the `L`/`R` steps to any positive rational (and
`stern_brocot_from_path` walks them back), which is a run-length view of the continued
fraction. `best_rational_bounded` descends the tree toward a real and returns the closest
fraction within a denominator bound — the Stern-Brocot route to the same answer as
`best_rational`. `farey_sequence(n)` lists every reduced fraction in `[0, 1]` with
denominator up to `n`, generated by the neighbour recurrence in one pass, and `mediant`
inserts the fraction between two others. Verified against brute nearest-fraction and Farey
enumerations over thousands of cases.

Long or badly-scaled sums lose low-order bits; the compensated routines recover them.
`neumaier_sum` (and `kahan_sum`) carry a running correction, `accurate_dot` does the
same for a dot product, and `welford` computes a one-pass mean and variance that stays
stable where the textbook formula cancels:

```python
from quantforge import neumaier_sum, accurate_dot, welford

neumaier_sum([1.0, 1e100, 1.0, -1e100])        # 2.0 — naive summation loses the 1s
accurate_dot([1e8, 1, -1e8], [1, 1, 1])        # 1.0 — exact despite the cancellation
welford([2, 4, 4, 4, 5, 5, 7, 9])              # (mean 5.0, variance 4.5714, n 8)
```

`neumaier_sum` is the robust default (correct even when a term dwarfs the running
total); summing a million `0.1`s it stays exact where the naive loop drifts by ~1e-6.
`welford` returns `(mean, variance, n)` in a single pass and recovers the variance of
large-mean data (e.g. values near `1e9`) where `mean(x²) − mean(x)²` collapses to zero
through catastrophic cancellation.

`RunningMoments` extends the one-pass idea to the third and fourth moments — mean,
variance, skewness and excess kurtosis as data streams in — and two accumulators
*merge* exactly, so partial results from parallel chunks combine into the whole:

```python
from quantforge import RunningMoments

rm = RunningMoments([2, 4, 4, 4, 5, 5, 7, 9])
rm.mean, rm.variance()        # 5.0, 4.5714
rm.skewness(), rm.kurtosis()  # 0.6563, -0.2188

a = RunningMoments([2, 4, 4, 4])
b = RunningMoments([5, 5, 7, 9])
(a + b).mean                  # 5.0 — same as the combined sample, computed by merge
```

Feed values with `update(x)` or an iterable to the constructor; read any moment at any
time. The `+` operator uses the Chan/Terriberry parallel-combination formulas, so
splitting a sample across workers and merging their accumulators reproduces the
single-pass result to ~1e-10 — the associativity that makes it a true parallel
reduction for large or distributed data.

`RunningCovariance` extends the same one-pass, mergeable idea to *pairs*, tracking
covariance and Pearson correlation of a stream of `(x, y)`:

```python
from quantforge import RunningCovariance

rc = RunningCovariance([1, 2, 3, 4, 5], [2, 4, 6, 8, 10])
rc.correlation()                     # 1.0 — perfectly linear
rc.covariance()                      # 5.0 (sample, ddof=1)
(RunningCovariance([1,2],[2,4]) + RunningCovariance([3,4,5],[6,8,10])).correlation()  # 1.0
```

`update(x, y)` folds in one pair via the Welford co-moment recurrence (stable, no
catastrophic cancellation), and `covariance`/`correlation`/`variance_x`/`variance_y` read
the current estimate — matching batch formulas. Two accumulators merge exactly with `+`,
so partial covariances from parallel chunks combine into the whole.

`RunningRegression` carries the same one-pass, mergeable machinery one step further and
*fits a line*: ordinary least squares of `y` on `x`, reading slope, intercept,
correlation, R-squared, and predictions at any point:

```python
from quantforge import RunningRegression

r = RunningRegression([1, 2, 3, 4, 5], [2.1, 3.9, 6.2, 7.8, 10.1])
r.slope(), r.intercept()             # 1.99, 0.05
r.r_squared()                        # 0.9973
r.predict(6)                         # 11.99
(RunningRegression([1, 2], [2.1, 3.9]) +
 RunningRegression([3, 4, 5], [6.2, 7.8, 10.1])).slope()   # 1.99 — same as the whole
```

`update(x, y)` folds one pair via the Welford co-moment recurrence, and `slope`,
`intercept`, `predict`, `correlation`, and `r_squared` read the current fit — identical to
`statistics.linear_regression` on the same data, with `r_squared` equal to the squared
correlation. As with the other accumulators, `+` merges two fits exactly, so chunks
regressed in parallel combine into the whole-sample fit. A slope needs at least two points
with spread in `x`; a constant `x` raises rather than dividing by zero.

When each observation carries a weight — reliability, sampling, or scenario probability —
the weighted summaries apply:

```python
from quantforge import weighted_mean, weighted_std, weighted_median

weighted_mean([1, 2, 3], [1, 1, 4])       # 2.5 — the heavy third value pulls it up
weighted_median([1, 2, 3], [1, 1, 10])    # 2.82 — mass concentrated near 3
weighted_std([1, 2, 3], [1, 1, 4])        # reliability-weight bias-corrected
```

`weighted_mean`, `weighted_variance`/`weighted_std`, `weighted_quantile`, and
`weighted_median` normalize the weights internally (they need not sum to 1). Equal weights
reproduce the ordinary statistics, and integer weights match the statistics of the
replicated sample exactly; the unbiased variance uses the reliability-weight correction
`V1 / (V1^2 - V2)`.

Quantiles need sorting — unless you estimate them online. `P2Quantile` tracks a single
quantile in O(1) memory (no data stored), and `reservoir_sample` draws a uniform sample
from a stream of unknown length in one pass:

```python
from quantforge import P2Quantile, reservoir_sample

q = P2Quantile(0.95)                    # track the 95th percentile
for x in stream:
    q.update(x)
q.value()                               # running estimate, within ~1% on a large stream

reservoir_sample(range(1_000_000), k=100)   # 100 items, each equally likely, one pass
```

`P2Quantile` uses the Jain-Chlamtac P-square algorithm — five markers that shift by a
piecewise-parabolic rule — and on a large uniform stream its median lands near 0.5 and
its 95th percentile matches the true order statistic to about 0.01, all without holding
the data. `reservoir_sample` (Vitter's algorithm R) is deterministic for a fixed `seed`
and returns every stream element with equal probability; if the stream is shorter than
`k` it returns all of it.

`P2Quantile` is an *estimate* in constant memory; when the *exact* running median is needed,
`RunningMedian` keeps it with two heaps:

```python
from quantforge import RunningMedian

rm = RunningMedian()
for x in (10, 20, 30):
    rm.push(x)
rm.median()                     # 20.0
RunningMedian([1, 2, 3, 4]).median()   # 2.5 — averages the two middle values
```

`RunningMedian` splits the values into a max-heap of the lower half and a min-heap of the
upper half, rebalanced after each `push` so their sizes differ by at most one — so `push` is
`O(log n)`, `median` is `O(1)`, and the answer is exact (the two middle values are averaged
on an even count). Unlike `P2Quantile` it stores every value, so use it when correctness
matters more than the constant memory. Verified against `statistics.median` after each push
over thousands of random streams.

`P2Quantile` tracks *one* preset quantile; when you need to query *any* quantile after the
fact with a guaranteed error, `DDSketch` keeps a compact histogram of log-spaced buckets
and answers each quantile within a fixed *relative* error:

```python
from quantforge import DDSketch

s = DDSketch(alpha=0.01)             # 1% relative-error accuracy
for i in range(1, 10001):
    s.add(float(i))
s.quantile(0.5)                      # ~4965 (true 5000, within 1%)
s.quantile(0.99)                     # ~9801 (true 9900, within 1%)
s.num_buckets()                      # 315 — memory grows with the log-range, not n
```

`DDSketch` (Masson-Rim-Lee) maps each positive value to bucket
`ceil(log(x) / log(gamma))` with `gamma = (1 + alpha)/(1 - alpha)`, so the reported
`q`-quantile is always within `alpha` *relative* error of the true value — tight for small
values, proportionally looser for large ones, which suits latencies and sizes that span
orders of magnitude. `min`/`max` are tracked exactly, and two sketches built with the same
`alpha` merge exactly with `+`, so `(a + b).quantile(0.9)` on split halves matches the
whole-stream estimate. Positive values only.

When elements should be drawn *proportional to a weight* rather than uniformly,
`weighted_reservoir_sample` selects `k` distinct items in the same single pass, and
`weighted_sample_with_replacement` draws `k` independent weight-proportional items:

```python
from quantforge import weighted_reservoir_sample, weighted_sample_with_replacement

weighted_reservoir_sample(['a', 'b', 'c', 'd'], [1, 2, 3, 4], k=2, seed=1)   # ['c', 'd']
weighted_sample_with_replacement(['a', 'b', 'c'], [1, 1, 8], k=10, seed=2)   # mostly 'c'
```

`weighted_reservoir_sample` is the Efraimidis-Spirakis A-Res algorithm: it assigns each
item the key `u^(1/w)` for a uniform `u` and keeps the `k` largest keys in a size-`k`
min-heap, so it runs in `O(n log k)` time and `O(k)` memory over a stream of unknown
length. For `k = 1` this reduces to a weighted choice — across 30,000 seeds the four
weights `1, 2, 3, 4` are picked about `0.10, 0.20, 0.30, 0.40` of the time, matching the
closed-form `w_i / sum(w)`. Zero-weight items are never selected, and `k` past the number
of positive-weight items returns all of them. `weighted_sample_with_replacement` uses
cumulative-weight bisection, so an item can repeat and its frequency converges to its
weight fraction. Both are deterministic for a fixed `seed`.

When recent data should count more than old — a drifting mean or a changing volatility
— `EWMAStats` keeps an exponentially-weighted mean and variance online (the RiskMetrics
recursion), with no fixed window:

```python
from quantforge import EWMAStats, ewma

vol = EWMAStats(lam=0.94)        # ~ a month of daily memory
for r in returns:
    vol.update(r)
vol.mean, vol.std()              # running EW mean and volatility

ewma([1, 2, 3, 4, 5], lam=0.94)  # [1.0, 1.06, 1.1764, 1.3458, 1.5651]
```

Larger `lam` means longer memory (effective window about `1 / (1 - lam)`, so `0.94` ≈ 17
observations). The variance follows the RiskMetrics form — squared deviation from the
previous mean — so `EWMAStats` reacts within a handful of points when volatility shifts
regime, unlike a long fixed-window estimate that would smear the change. The batch
`ewma` returns the running mean at every step for the whole series.

Counting *distinct* items exactly costs memory proportional to the count; `HyperLogLog`
estimates it from a few kilobytes of registers regardless of stream size:

```python
from quantforge import HyperLogLog

hll = HyperLogLog(p=14)               # 2^14 registers, ~0.8% typical error
for item in stream:                   # any stringifiable value
    hll.add(item)
hll.count()                           # estimated distinct count

a.merge(b)                            # union cardinality of two sketches (same p)
```

Each item is hashed; the leading-zero count of the hash tail feeds a register chosen by
its leading bits, and the harmonic mean across registers estimates the cardinality.
Accuracy is about `1.04 / sqrt(2^p)` — under 1% at `p = 14` — and stays within a few
percent from a hundred to a hundred million distinct items. Duplicates never inflate the
estimate, sketches `merge` into the union count, and the SHA-1 hash makes it
deterministic across runs.

Two more sketches cover *frequency* and *membership*. `CountMinSketch` estimates how
often each item appeared (never under-counting), and `BloomFilter` tests set membership
with no false negatives:

```python
from quantforge import CountMinSketch, BloomFilter

cms = CountMinSketch(width=2048, depth=5)
for item in stream:
    cms.add(item)
cms.estimate("hot_key")           # frequency estimate, an upper bound on the true count

bf = BloomFilter(capacity=10000, error_rate=0.01)
bf.add("seen_id")
"seen_id" in bf                   # True; a never-added id is True only ~1% of the time
```

`CountMinSketch` hashes each item into one counter per row and returns the minimum, so
collisions can only *over*-count — a heavy hitter's estimate is essentially exact while
the memory stays fixed. `BloomFilter` sizes its bit array and hash count from the target
`capacity` and `error_rate`; a member always tests present (no false negatives) and a
non-member tests present only at about the configured false-positive rate. Both are
deterministic through the shared SHA-1 hashing.

When the question is *which items dominate* rather than how many there are,
`MisraGries` and `SpaceSaving` track the frequent items in `O(k)` counters — the
heavy-hitters problem:

```python
from quantforge import MisraGries, SpaceSaving

mg = MisraGries(3)
for x in "aaaaaaabbbbbccccdde":     # a=7 b=5 c=4 d=2 e=1, n=19
    mg.add(x)
mg.counts()                         # {'a': 4, 'b': 2, 'c': 1} — undercounts, but keeps the leaders

ss = SpaceSaving(3)
for x in "aaaaaaabbbbbccccdde":
    ss.add(x)
ss.top(2)                           # [('a', 7), ('d', 6)] — 'd' overcounts from an evicted slot
ss.guaranteed(2)                    # [('a', 7), ('d', 2)] — count - error, a true lower bound
```

`MisraGries` keeps `k` counters and, when a new item arrives with all counters full,
decrements every counter — so any item occurring more than `n / (k + 1)` times is
guaranteed to survive, and a stored count never exceeds the true count. `SpaceSaving`
instead evicts the current minimum slot and lets the new item inherit its count as an
`error`, so its count is an *upper* bound and `count - error` a guaranteed lower bound;
`top(m)` ranks the estimates and `guaranteed(m)` returns those lower bounds. On a skewed
stream Space-Saving recovers the true top items even when their raw counts drift high.

`log(sum(exp(x)))` overflows the instant any `x` is large; `logsumexp` shifts by the
maximum first, so it stays exact where the naive form returns `inf`, and `softmax` /
`log_softmax` build the normalized (log-)probability transforms on it:

```python
from quantforge import logsumexp, softmax, log_softmax

logsumexp([1000.0, 1001.0, 1002.0])   # 1002.407606 — naive exp() would overflow to inf
softmax([1.0, 2.0, 3.0])              # [0.09, 0.2447, 0.6652], sums to 1
log_softmax([1.0, 2.0, 3.0])          # x_i - logsumexp(x), avoids log(0)
```

`logsumexp` takes optional non-negative `weights` for a log-weighted-sum-exp (mixture
likelihoods); `softmax` is stable even on 1000-scale inputs, and `log_softmax` computes
`x_i - logsumexp(x)` directly rather than the overflow-prone `log(softmax(x))` — the
log-likelihood form used in classification and Boltzmann-style weighting.

The special functions behind the distribution routines are public:

```python
from quantforge import gammainc, gammaincc, betainc, digamma, erfinv

gammainc(0.5, x=1.0)        # regularized lower incomplete gamma P(a, x) = erf(1)
gammaincc(0.5, 1.0)         # upper Q(a, x) = 1 - P
betainc(2.0, 3.0, x=0.4)    # regularized incomplete beta I_x(a, b)
digamma(1.0)                # psi(1) = -euler_gamma
erfinv(0.95)                # inverse error function
```

`gammainc`/`gammaincc` give the gamma, chi-square and Poisson CDFs; `betainc`
gives the Student-t, F and binomial CDFs; `erfinv` gives normal quantiles
(`erfinv(y) = norm_ppf((1+y)/2) / sqrt(2)`). They satisfy `P + Q = 1`,
`I_x(a,b) = 1 - I_{1-x}(b,a)`, and `erf(erfinv(y)) = y` to machine precision.

The correlated-normal CDFs behind the multi-asset and compound-option models are public
too:

```python
from quantforge import bivariate_normal_cdf, trivariate_normal_cdf

bivariate_normal_cdf(0.0, 0.0, 0.5)                 # 0.333333 = 1/4 + asin(0.5)/(2 pi)
trivariate_normal_cdf(0.0, 0.0, 0.0, 0.3, 0.3, 0.3) # 0.19775
```

`bivariate_normal_cdf(a, b, rho)` is the Drezner-Wesolowsky single-integral form
(accurate to ~1e-7); `trivariate_normal_cdf(a, b, c, r12, r13, r23)` reduces the third
variable to a one-dimensional integral of the bivariate CDF (Genz). Both reduce to the
product of marginals at zero correlation and are monotone in each correlation — the
ingredients for spread, exchange, and two-factor compound-option pricing.

## Range-query structures

Fast prefix/range queries with updates, for cumulative-frequency and running-total work:

```python
from quantforge import FenwickTree, SegmentTree

ft = FenwickTree([1, 2, 3, 4, 5])
ft.range_sum(1, 3)                    # 9 — sum of elements 1..3
ft.update(0, 10)                       # add 10 to element 0, O(log n)

st = SegmentTree([1, 2, 3, 4, 5], combine=max, identity=float("-inf"))
st.query(0, 2)                         # 3 — range maximum over [0, 2]
```

`FenwickTree` (binary indexed tree) does `O(log n)` point updates and prefix/range sums in
tight space. `SegmentTree` generalizes to any associative `combine` — the default is sum,
and passing `min`/`max` (with the matching `identity`) gives range-minimum/maximum with
`O(log n)` point updates. Both are verified against a brute-force recompute across random
arrays with interleaved updates.

`SegmentTree` updates one index at a time; `LazySegmentTree` updates a whole *range* at once
— add a value to every element of `[lo, hi)` and query the sum/min/max over any range, each
in `O(log n)`, via lazy propagation:

```python
from quantforge import LazySegmentTree

st = LazySegmentTree([1, 2, 3, 4, 5], "sum")
st.query(0, 5)              # 15
st.update(1, 4, 10)         # add 10 to indices 1..3 -> [1, 12, 13, 14, 5]
st.query(0, 5), st.query(1, 4)   # (45, 39)

LazySegmentTree([5, 3, 8, 1, 9], "min").query(0, 5)   # 1
```

A range update tags the interior nodes it fully covers and defers the work; the tag is only
pushed down when a later query descends through it, so both operations stay `O(log n)` no
matter how large the range. `mode` selects the aggregate (`"sum"`, `"min"`, `"max"`), ranges
are half-open, and non-power-of-two lengths are handled. Verified against a brute array over
thousands of interleaved update/query sequences for all three modes.

When the "elements" are *lines* rather than array cells, `LiChaoTree` keeps their lower (or
upper) envelope, answering "smallest `y = m x + b` at this `x`" in `O(log n)` with lines
added in any order:

```python
from quantforge import LiChaoTree

t = LiChaoTree([0, 1, 2, 3, 4])
t.add_line(0, 5)      # y = 5
t.add_line(2, -1)     # y = 2x - 1
t.add_line(-2, 8)     # y = -2x + 8
t.query(0), t.query(4)   # (-1, 0) — the minimum over all three lines

hi = LiChaoTree([0, 1, 2, 3, 4], maximize=True)
for m, b in [(0, 5), (2, -1), (-2, 8)]:
    hi.add_line(m, b)
hi.query(0), hi.query(4) # (8, 7) — the upper envelope
```

Construct with the `x` values you will query; then `add_line(m, b)` inserts a line and
`query(x)` returns the minimum (or maximum, with `maximize=True`) over all inserted lines.
Unlike the convex-hull trick, lines may arrive in any slope order. This is the standard
accelerator for the DP recurrence `dp[i] = min_j (m_j x_i + b_j)` — turning an `O(n^2)` scan
into `O(n log n)`. Verified against a brute minimum/maximum over all lines across thousands
of random cases, insertion order included.

When the array never changes, a `SparseTable` answers *idempotent* range queries
(minimum, maximum, gcd) in `O(1)` — faster per query than a segment tree, at the cost of
being read-only:

```python
from quantforge import range_min_query, range_max_query, range_gcd_query

a = [3, 1, 4, 1, 5, 9, 2, 6]
range_min_query(a).query(1, 4)          # 1 — min over indices 1..4 inclusive
range_max_query(a).query(2, 5)          # 9
range_gcd_query([12, 18, 24, 6]).query(0, 2)   # 6 — gcd(12, 18, 24)
```

`SparseTable` precomputes the answer for every power-of-two block in `O(n log n)`, then a
query combines the two (possibly overlapping) blocks that cover `[lo, hi]`. That overlap is
harmless only for *idempotent* combiners — where `f(x, x) = x` — so min, max, and gcd work
but sum does not (use `FenwickTree` for sums). The convenience builders cover the common
cases, and any idempotent `combine` is accepted (e.g. bitwise-or). Verified against
brute-force recomputation over thousands of random arrays and every possible range.

For plain *sums* over a static array — where a Fenwick tree's update support is unnecessary
— a prefix-sum table answers any range in `O(1)` after one linear pass, in 1-D or 2-D, and
a difference array does the mirror image (many range-adds, one read):

```python
from quantforge import PrefixSum1D, PrefixSum2D, DifferenceArray

PrefixSum1D([1, 2, 3, 4, 5]).range_sum(1, 4)          # 9 — sum of [2, 3, 4]

grid = PrefixSum2D([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
grid.range_sum(1, 1, 3, 3)                            # 28 — the bottom-right 2x2 block

da = DifferenceArray(5)
da.add(1, 4, 10); da.add(0, 2, 1)
da.result()                                           # [1, 11, 10, 10, 0]
```

`PrefixSum1D.range_sum(lo, hi)` subtracts two cumulative entries; `PrefixSum2D` is the
summed-area table (integral image) that answers any axis-aligned rectangle by
inclusion-exclusion of four corners — the trick behind box filters and Haar features. Both
use half-open ranges. `DifferenceArray` inverts the pattern: each `add(lo, hi, delta)` marks
only the two endpoints in `O(1)`, and one `result()` pass accumulates them into the final
array — ideal when a batch of range updates precedes a single read. All three checked
against brute-force sums over thousands of random instances.

When the 2-D grid *changes* — a running heatmap, a collision counter — `FenwickTree2D` keeps
`PrefixSum2D`'s rectangle-sum queries while also allowing point updates, each in
`O(log R · log C)`:

```python
from quantforge import FenwickTree2D

ft = FenwickTree2D(3, 3)
for i in range(3):
    for j in range(3):
        ft.add(i, j, i * 3 + j + 1)   # fill with 1..9
ft.range_sum(0, 0, 3, 3), ft.range_sum(1, 1, 3, 3)   # (45, 28)
ft.add(1, 1, 100)
ft.range_sum(1, 1, 2, 2)              # 105 — reflects the update immediately
```

`add(r, c, delta)` adjusts a cell and `range_sum(r0, c0, r1, c1)` sums a half-open rectangle
by four-corner inclusion-exclusion — the 2-D binary indexed tree, nesting the `i & -i`
traversal in both axes. Build from a shape or an initial grid. Where `PrefixSum2D` is fixed
after construction, `FenwickTree2D` trades its `O(1)` query for `O(log R · log C)` in
exchange for updates. Verified against a brute grid over thousands of interleaved
update/query sequences.

Finding a k-th order statistic or the top-k does not need a full sort — quickselect does
it in linear time:

```python
from quantforge import kth_smallest, median, top_k

kth_smallest([5, 2, 8, 1, 9], 1)      # 2 — the 2nd smallest (0-indexed)
median([1, 2, 3, 4])                   # 2.5
top_k([5, 2, 8, 1, 9, 3], 3)          # [9, 8, 5] — the three largest, sorted
```

`kth_smallest` uses quickselect with the median-of-medians pivot for a guaranteed `O(n)`
worst case (matching a full sort for every k, with or without duplicates), `median` is its
even-length-averaging wrapper, and `top_k` selects the threshold then sorts only the
chosen elements (`largest=False` for the bottom k).

The extremum of every length-`k` window of a sequence comes in linear time — no per-window
rescan — from `sliding_window_min`/`sliding_window_max` (and `sliding_window_sum`):

```python
from quantforge import sliding_window_min, sliding_window_max, sliding_window_sum

v = [1, 3, -1, -3, 5, 3, 6, 7]
sliding_window_max(v, 3)          # [3, 3, 5, 5, 6, 7]
sliding_window_min(v, 3)          # [-1, -3, -3, -3, 3, 3]
sliding_window_sum([1, 2, 3, 4], 2)   # [3, 5, 7]
```

`sliding_window_max`/`min` keep a monotonic deque of candidate indices — discarding any that
a newer, better value dominates and dropping those that fall out of the window — so each
element is pushed and popped once, giving `O(n)` overall instead of the `O(n·k)` of scanning
every window. `sliding_window_sum` slides a running sum. All return `n - k + 1` results and
are verified against a brute per-window scan over thousands of random cases.

Collections of `(start, end)` ranges have their own set operations for scheduling and
coverage analysis:

```python
from quantforge import (merge_intervals, total_covered_length,
                        intervals_intersection, max_overlap)

merge_intervals([(1, 3), (2, 6), (8, 10)])           # [(1, 6), (8, 10)]
total_covered_length([(1, 4), (2, 6), (8, 10)])      # 7 (overlaps counted once)
intervals_intersection([(1, 5), (8, 10)], [(3, 9)])  # [(3, 5), (8, 9)]
max_overlap([(1, 5), (2, 6), (4, 8), (10, 12)])      # 3 — peak simultaneous overlap
```

`merge_intervals` collapses overlapping and touching ranges into a minimal cover;
`intervals_union`/`intervals_intersection` combine two collections (verified against an
integer point-set); `total_covered_length` sums the merged cover; and `max_overlap` sweeps
the endpoints for the most intervals active at once (matching a brute-force point scan).

*Choosing* among overlapping intervals is the scheduling problem, and three routines solve
its standard variants:

```python
from quantforge import weighted_interval_schedule, activity_selection, min_rooms

jobs = [(1, 3, 5), (2, 5, 6), (4, 6, 5), (6, 7, 4)]   # (start, end, weight)
weighted_interval_schedule(jobs)   # (14.0, [(1,3,5), (4,6,5), (6,7,4)]) — max total weight

activity_selection([(1, 4), (3, 5), (0, 6), (5, 7), (8, 11), (12, 16)])
                                   # [(1,4), (5,7), (8,11), (12,16)] — most non-overlapping
min_rooms([(0, 30), (5, 10), (15, 20)])   # 2 — fewest resources to run all at once
```

`weighted_interval_schedule` picks the non-overlapping subset of maximum total weight — DP
over intervals sorted by end time, each looking back (by binary search) to the last one it
does not overlap, `O(n log n)` — and returns both the total and the chosen intervals.
`activity_selection` is the unweighted version: the most intervals you can fit, by the
earliest-finish-time greedy. `min_rooms` sweeps the endpoints for the peak number active
simultaneously — the fewest rooms, machines, or tracks needed. Intervals are half-open, so
one ending exactly as another starts do not clash. All three are checked against
brute-force subset enumeration over thousands of random instances.

For points in space rather than on a line, `KDTree` indexes a fixed set of
k-dimensional coordinates and answers nearest-neighbour, radius, and box queries without
scanning every point:

```python
from quantforge import KDTree

pts = [(2, 3), (5, 4), (9, 6), (4, 7), (8, 1), (7, 2)]
tree = KDTree(pts)
tree.nearest((9, 2))              # (4, 1.414...) — index 4 is point (8, 1)
tree.k_nearest((9, 2), 2)         # [(4, 1.414...), (5, 2.0)]
tree.within_radius((5, 5), 3.0)   # [(1, 1.0), (3, 2.236...)] — indices + distances
tree.range_search((4, 1), (8, 4)) # [1, 4, 5] — indices inside the box
```

The tree partitions on alternating axes and, on each query, skips any subtree whose
bounding slab is already farther than the current best — so a nearest-neighbour lookup is
typically `O(log n)` rather than the `O(n)` of a full scan. Queries return the original
point *indices* (in build order), so you can map results back to your own records:
`nearest`/`k_nearest` give `(index, distance)` sorted by distance, `within_radius` lists
every point inside a Euclidean ball, and `range_search` lists those inside an axis-aligned
box. All four are verified against brute-force distance computation over random point sets
in one to four dimensions.

## Interval arithmetic

An `Interval` tracks a range of possible values and propagates it through arithmetic so
the result is *guaranteed* to enclose the true answer — the basis of tolerance analysis
and worst-case bounds:

```python
from quantforge import Interval

Interval(1, 2) + Interval(3, 4)          # Interval(4, 6)
Interval(-2, 3) ** 2                       # Interval(0, 9) — straddles zero, so min is 0
Interval(2, 8).intersect(Interval(5, 10)) # Interval(5, 8)
Interval(2, 8).width()                     # 6.0 — the uncertainty
```

Arithmetic (`+ - * /`), integer powers, and the monotone functions `exp`/`log`/`sqrt`
all return an interval that contains the true value for every input in the operands —
verified to enclose the result across thousands of random samples. Division by an interval
spanning zero raises, and `width`/`midpoint`/`contains`/`intersect` round out the class.
(It uses ordinary floating point without outward rounding, so the enclosure is
mathematically correct but not a certified ULP-level bound.)

## Number theory

Exact integer primitives — primality, factorization, and the derived quantities — built
on Python's arbitrary-precision integers, so there is no overflow:

```python
from quantforge import is_prime, factorize, euler_totient, divisors, gcd, lcm

is_prime(2**61 - 1)                    # True — a Mersenne prime, far past 64-bit
is_prime(561)                          # False — a Carmichael number (fools weak tests)
factorize(60)                          # [(2, 2), (3, 1), (5, 1)]
factorize(1000000007 * 1000000009)     # [(1000000007, 1), (1000000009, 1)]
euler_totient(36)                      # 12
divisors(28)                           # [1, 2, 4, 7, 14, 28] — 28 is perfect
```

`is_prime` is a deterministic Miller-Rabin test with a fixed witness set that is exact
for every integer below `3.3e24` (well past 64-bit), so the answer is a proof, not a
probability — it correctly certifies the Mersenne prime `2^61-1` and rejects the
Carmichael number `561` that fools naive tests. `factorize` combines trial division with
Pollard's rho and certifies each factor prime, so it splits a `~10^18` semiprime that
trial division alone could not reach. `euler_totient` counts the integers coprime to `n`,
`divisors` lists every divisor (the proper divisors of `28` sum to `28`), and `gcd`/`lcm`
satisfy `gcd(a,b) * lcm(a,b) == a*b`.

When you need *many* primes rather than one test, the sieve is far faster:

```python
from quantforge import primes_up_to, prime_count, nth_prime, smallest_prime_factors

primes_up_to(20)         # [2, 3, 5, 7, 11, 13, 17, 19]
prime_count(100)         # 25 — pi(100)
nth_prime(1000)          # 7919
smallest_prime_factors(30)[12]   # 2 — least prime dividing 12
```

`primes_up_to` is the Sieve of Eratosthenes; `prime_count` is the prime-counting function
pi(n); `nth_prime` grows the sieve bound by the prime-number-theorem estimate until it has
enough primes; and `smallest_prime_factors` returns a table that factorizes any number in
range in `O(log n)` — matching `factorize` and `is_prime` exactly.

Modular ("clock") arithmetic adds the operations behind RSA-style key math and residue
computation:

```python
from quantforge import extended_gcd, mod_inverse, chinese_remainder, mod_pow, discrete_log

extended_gcd(240, 46)                  # (2, -9, 47): 240*-9 + 46*47 == 2 == gcd
mod_inverse(3, 11)                     # 4, since 3*4 == 12 == 1 (mod 11)
chinese_remainder([2, 3, 2], [3, 5, 7])# (23, 105): the unique x mod 105
mod_pow(3, -1, 5)                      # 2 — negative exponents invert first
discrete_log(3, 13, 17)                # 4, since 3^4 == 13 (mod 17)
```

`extended_gcd` returns the Bezout coefficients (`a*x + b*y == gcd`), from which
`mod_inverse` builds the modular inverse (raising when the inputs are not coprime).
`chinese_remainder` stitches congruences with pairwise-coprime moduli into a single one —
the classic `x == 2 (mod 3), 3 (mod 5), 2 (mod 7)` gives `23`. `mod_pow` wraps fast
exponentiation and handles negative exponents through the inverse, and `discrete_log`
solves `base^x == target` by baby-step giant-step in `O(sqrt(mod))`.

The quadratic-residue and multiplicative-group tools go deeper into prime arithmetic:

```python
from quantforge import (legendre_symbol, tonelli_shanks,
                        multiplicative_order, primitive_root)

legendre_symbol(2, 7)         # 1 — 2 is a quadratic residue mod 7
tonelli_shanks(2, 7)          # 4 — a square root (4^2 = 16 = 2 mod 7)
multiplicative_order(2, 7)    # 3 — 2^3 = 1 mod 7
primitive_root(7)             # 3 — a generator of the group mod 7
```

`legendre_symbol` (and its composite generalization `jacobi_symbol`) tells whether `a` is
a square mod `p`; `tonelli_shanks` computes the square root when it is (raising on a
non-residue); `multiplicative_order` is the least exponent giving 1; and `primitive_root`
returns a generator whose powers cover every non-zero residue.

Exact counting functions round out the section — all on big integers, so a coefficient
with thirty digits is exact rather than a float approximation:

```python
from quantforge import (binomial, multinomial, stirling_second,
                        bell, catalan, partition_count, derangements)

binomial(10, 3)                    # 120
multinomial([1, 4, 4, 2])          # 34650 — arrangements of "mississippi"
stirling_second(4, 2)              # 7 — partitions of 4 items into 2 blocks
[catalan(n) for n in range(7)]     # [1, 1, 2, 5, 14, 42, 132]
partition_count(8)                 # 22 — integer partitions of 8
derangements(4)                    # 9 — permutations with no fixed point
```

`binomial` and `multinomial` count subsets and multiset arrangements; `stirling_second`
counts set partitions into a fixed number of non-empty blocks, and `bell(n)` sums that
row over all block counts (`bell(5) == 52`). `catalan` gives the ubiquitous
`1, 1, 2, 5, 14, 42, ...`, `partition_count` the number of ways to write `n` as an
unordered sum, and `derangements(n)` the fixed-point-free permutations — equal to
`round(n!/e)` for every `n`.

Symmetric functions relate a multiset of values to the polynomial with those roots.
`elementary_symmetric`, `power_sums`, and `poly_from_roots` compute the three descriptions,
and Newton's identities convert between them exactly:

```python
from quantforge import (elementary_symmetric, power_sums, poly_from_roots,
                        power_to_elementary, elementary_to_power)

elementary_symmetric([1, 2, 3])   # [1, 6, 11, 6] — e_0..e_3
power_sums([1, 2, 3], 3)          # [3, 6, 14, 36] — p_0..p_3
poly_from_roots([1, 2, 3])        # [1, -6, 11, -6] — (t-1)(t-2)(t-3)
power_to_elementary([6, 14, 36])  # [1, 6, 11, 6] as exact Fractions
elementary_to_power([1, 6, 11, 6])# [6, 14, 36]
```

`elementary_symmetric` returns `e_k`, the sum of all products of `k` distinct values (an
`O(n^2)` DP), and `poly_from_roots` reads the monic coefficients off them via Vieta
(`(-1)^k e_k`, highest-degree first, matching `polynomial_roots`). `power_sums` gives
`p_k = sum x_i^k`. `power_to_elementary` and `elementary_to_power` are Newton's identities
in exact `Fraction` arithmetic, so the two representations round-trip without rounding —
the bridge between a polynomial's coefficients and the moments of its roots. Verified
against combination-expansion references over thousands of multisets.

Polynomials over a prime field are exact — the basis of secret sharing and error-correcting
codes. `lagrange_interpolate_mod` recovers the polynomial through a set of points, and
`poly_eval_mod`/`poly_add_mod`/`poly_mul_mod` are the field arithmetic:

```python
from quantforge import lagrange_interpolate_mod, poly_eval_mod, poly_mul_mod

lagrange_interpolate_mod([(0, 1), (1, 3), (2, 7)], 997)   # [1, 1, 1] — x^2 + x + 1
poly_eval_mod([1, 2, 3], 2, 997)                          # 17 — 1 + 2*2 + 3*4 (Horner)
poly_mul_mod([1, 1], [1, 1], 997)                         # [1, 2, 1]

# Shamir secret sharing: the secret is the constant term; any t shares recover it
mod = 998244353
shares = [(x, poly_eval_mod([123456, 111, 222], x, mod)) for x in (1, 2, 3)]
lagrange_interpolate_mod(shares, mod)[0]                  # 123456
```

Coefficients are lowest-degree first (`c[0]` is the constant term).
`lagrange_interpolate_mod` builds the unique degree-`< n` polynomial through `n` points with
distinct `x`, using exact modular inverses — so splitting a secret into shares and
reconstructing it from any threshold subset round-trips exactly, and Reed-Solomon
encode/decode falls out the same way. `poly_mul_mod` is the schoolbook product mod the prime
(pair it with `convolve_mod` for the NTT-accelerated version). Verified over thousands of
random cases across two primes.

Bernoulli numbers and Faulhaber's formula give the exact closed form for sums of powers:

```python
from quantforge import bernoulli_number, faulhaber

bernoulli_number(2)          # Fraction(1, 6)  (B_1 = +1/2; odd B_n > 1 vanish)
faulhaber(10, 2)             # 385   = 1^2 + 2^2 + ... + 10^2
faulhaber(100, 3)            # 25502500 = (100*101/2)^2
```

`bernoulli_number` returns exact `Fraction` values (via the standard recurrence,
`B_1 = +1/2`), and `faulhaber(m, p)` evaluates `1^p + ... + m^p` in closed form — matching
the direct sum for every degree and reproducing the classic `m(m+1)/2`,
`m(m+1)(2m+1)/6`, and `(m(m+1)/2)^2` identities. `bernoulli_sequence(n)` returns
`B_0 .. B_n`.

Numeral-system conversions round out the section — any base and Roman numerals:

```python
from quantforge import to_base, from_base, to_roman, from_roman

to_base(255, 16)          # "ff"
from_base("1010", 2)      # 10
to_roman(1994)            # "MCMXCIV"
from_roman("MMXXIV")      # 2024
```

`to_base`/`from_base` convert non-negative integers to and from any base `2..36` (digits
`0-9a-z`, case-insensitive), round-tripping and matching Python's `int(s, base)`.
`to_roman`/`from_roman` handle `1..3999` in subtractive notation and round-trip for every
value in range.

Combinatorial *ranking* maps permutations and subsets to integers and back, so an object
can be stored or generated by index without enumerating them all:

```python
from quantforge import (gray_code, permutation_unrank, permutation_rank,
                        combination_unrank, combination_rank)

[gray_code(i) for i in range(4)]     # [0, 1, 3, 2] — successive codes differ by one bit
permutation_unrank(3, 3)             # [1, 2, 0] — the 4th permutation of [0,1,2]
permutation_rank([2, 0, 1])          # 4
combination_unrank(4, 5, 2)          # [1, 2] — the 5th 2-subset of [0..4]
```

`gray_code`/`gray_decode` are the reflected binary Gray code (consecutive values differ in
one bit); `permutation_rank`/`permutation_unrank` use the Lehmer code and
`combination_rank`/`combination_unrank` the combinatorial number system — both matching
`itertools`' lexicographic order and round-tripping exactly.

Integer equations close out the section — linear Diophantine, square-root continued
fractions, and Pell:

```python
from quantforge import linear_diophantine, sqrt_continued_fraction, pell_fundamental

linear_diophantine(3, 6, 9)     # (3, 0, 2, -1): x0,y0 + k*(dx,dy) solves 3x+6y=9
sqrt_continued_fraction(7)      # (2, [1,1,1,4]) — sqrt(7) = [2; 1,1,1,4, ...]
pell_fundamental(61)            # (1766319049, 226153980) — the notorious x^2-61y^2=1
```

`linear_diophantine` returns a particular solution plus the step generating all others
(or `None` when `gcd(a,b)` does not divide `c`). `sqrt_continued_fraction` gives the
eventually-periodic CF of an irrational square root, and `pell_fundamental` reads the
smallest positive `x^2 - n y^2 = 1` solution off its convergents — matching the classic
hard cases like `n = 61`.

Integers form a vector space over GF(2) under XOR, and `XorBasis` maintains a basis of the
span of a set of numbers — answering the maximum-subset-XOR and related questions in time
linear in the bit width:

```python
from quantforge import XorBasis

xb = XorBasis([1, 2, 4])
xb.max_xor()             # 7 — the largest XOR of any subset
xb.rank()                # 3 independent vectors
xb.count_distinct()      # 8 reachable values (2^rank)
[xb.kth_smallest(k) for k in range(4)]   # [0, 1, 2, 3]

xb2 = XorBasis([5, 3, 6])
xb2.rank()               # 2 — 6 == 5 ^ 3, so it adds nothing
xb2.can_represent(6)     # True
```

`insert` reduces each number against the current basis by leading bit (Gaussian
elimination over GF(2)) and reports whether it was independent. `max_xor`/`min_xor` greedily
pick basis vectors that raise or lower the running value (with an optional `start` to XOR
against), `can_represent` tests subset-XOR membership, `rank`/`count_distinct` give the
dimension and the `2^rank` reachable values, and `kth_smallest` indexes the sorted reachable
set. Verified against brute subset-XOR enumeration over 4000 random value sets.

Full linear systems over GF(2) — each equation an XOR of boolean variables equal to 0 or 1
— are solved by `solve_gf2`, with `gf2_rank` and `gf2_nullspace_basis` for the structure:

```python
from quantforge import solve_gf2, gf2_rank, gf2_nullspace_basis

# bit j of an equation mask marks variable j; rhs is the 0/1 right-hand side
solve_gf2([0b011, 0b110, 0b001], [1, 0, 1], 3)   # [1, 0, 0]
solve_gf2([0b1, 0b1], [0, 1], 1)                 # None — 0 == 1 is inconsistent
gf2_rank([0b001, 0b010, 0b100])                  # 3
gf2_nullspace_basis([0b0011, 0b1100], 4)         # [[1,1,0,0], [0,0,1,1]]
```

`solve_gf2` runs Gaussian elimination with rows bit-packed into integers, so each
elimination step is a single XOR; it returns one solution (free variables set to 0) or
`None` when the system is inconsistent. `gf2_rank` counts independent rows, and
`gf2_nullspace_basis` returns `n_vars - rank` vectors spanning `{x : A x = 0}` — the tools
behind lights-out puzzles, linear-code decoding, and cycle-space computations. Verified
against brute-force enumeration over 4000 random systems.

Any linear recurrence advances by a fixed matrix, so `linear_recurrence_nth` reaches the
`n`-th term in `O(log n)` matrix multiplies — `matrix_power` does the exponentiation, and
`fibonacci` is the canonical case:

```python
from quantforge import fibonacci, linear_recurrence_nth, matrix_power

fibonacci(100)                              # 354224848179261915075
fibonacci(1_000_000, mod=1_000_000_007)     # 918091266 — huge index, still instant
linear_recurrence_nth([1, 1, 1], [0, 0, 1], 10)   # 81 — tribonacci
matrix_power([[1, 1], [1, 0]], 5)           # [[8, 5], [5, 3]]
```

`linear_recurrence_nth(coeffs, initial, n)` evaluates `x_i = sum_j coeffs[j] * x_{i-1-j}`
by raising the recurrence's companion matrix to the `n`-th power, so the cost grows with
`log n` rather than `n` — the difference between instant and hopeless for indices in the
millions. `matrix_power` is exponentiation-by-squaring on an integer or real square matrix,
and every routine takes an optional `mod` for exact modular arithmetic. Verified against
direct iteration over thousands of random recurrences and known Fibonacci/tribonacci values.

Where `linear_recurrence_nth` advances a *known* recurrence, `berlekamp_massey` runs the
other way: given a sequence's terms modulo a prime, it *finds* the shortest recurrence they
obey, and `berlekamp_massey_next` extrapolates from it:

```python
from quantforge import berlekamp_massey, berlekamp_massey_next

MOD = 998244353
berlekamp_massey([0, 1, 1, 2, 3, 5, 8, 13], MOD)   # [1, 1] — the Fibonacci rule
berlekamp_massey([1, 2, 4, 8, 16], MOD)            # [2]    — doubling
berlekamp_massey_next([1, 1, 2, 3, 5, 8, 13], MOD, 3)   # [21, 34, 55]
```

`berlekamp_massey` returns the coefficients of the shortest linear recurrence
`s_i = sum_j c_j s_{i-1-j} (mod p)` reproducing the input — the minimal LFSR — in `O(n^2)`;
given at least `2L` terms of a length-`L` recurrence it recovers it exactly. It is the tool
for guessing a closed form from a computed prefix or decoding a linear sequence.
`berlekamp_massey_next` rolls the discovered recurrence forward to predict further terms.
Verified over thousands of random modular recurrences that the recovered rule reproduces the
sequence, is minimal, and predicts the true continuation.

Dates reduce to integer arithmetic once mapped to a Julian day number, so `calendar_math`
answers weekday, day-count, leap-year, and Easter questions with no `datetime` dependency:

```python
from quantforge import (julian_day_number, day_of_week_name, days_between,
                        easter_date, add_days, is_leap_year)

julian_day_number(2000, 1, 1)          # 2451545 — the standard reference epoch
day_of_week_name(2026, 9, 9)           # 'Wednesday'
days_between((2024, 1, 1), (2024, 12, 31))   # 365
easter_date(2025)                      # (2025, 4, 20)
add_days(2024, 2, 28, 2)               # (2024, 3, 1) — leap-year aware
is_leap_year(2000), is_leap_year(1900) # (True, False)
```

`julian_day_number`/`jdn_to_date` are the Fliegel-Van Flandern proleptic-Gregorian
conversion — a continuous integer day count, so `days_between` is one subtraction and
`day_of_week` is `jdn % 7`. `add_days` shifts a date by any offset, `day_of_year` and
`days_in_month` round out the arithmetic, and `easter_date` is the Gregorian Computus. All
verified against Python's `datetime`/`calendar` over 8000 random dates; being proleptic, the
routines are not bounded by `datetime`'s year ≤ 9999.

## Computational geometry

Planar geometry primitives on lists of `(x, y)` points — the shape of a point cloud, its
area, and spatial queries:

```python
from quantforge import (convex_hull, polygon_area, polygon_centroid,
                        point_in_polygon, closest_pair)

convex_hull([(0,0), (1,0), (1,1), (0,1), (0.5,0.5)])  # [(0,0),(1,0),(1,1),(0,1)]
polygon_area([(0,0), (4,0), (0,3)])                    # 6.0 (shoelace)
polygon_centroid([(0,0), (6,0), (0,3)])                # (2.0, 1.0)
point_in_polygon((2,2), [(0,0),(4,0),(4,4),(0,4)])     # True
closest_pair([(0,0), (5,5), (1,1), (9,9)])             # ((0,0), (1,1), 1.414...)
```

`convex_hull` uses Andrew's monotone chain to return the counter-clockwise hull (interior
and collinear points dropped) — so the hull of a square plus interior points is just its
four corners. `polygon_area` and `polygon_centroid` apply the shoelace formula and are
orientation-independent. `point_in_polygon` ray-casts with the odd-crossing rule and
handles concave shapes and on-edge points, and `closest_pair` finds the nearest two points
in `O(n log n)`, matching a brute-force search exactly.

Segment intersection and window clipping round out the toolkit — the primitives behind
collision tests and viewport clipping:

```python
from quantforge import (segments_intersect, segment_intersection,
                        polygon_perimeter, clip_polygon)

segment_intersection((0,0),(1,1),(0,1),(1,0))         # (0.5, 0.5)
segments_intersect((0,0),(1,0),(2,2),(3,3))           # False (disjoint)
polygon_perimeter([(0,0),(4,0),(0,3)])                # 12.0 (3-4-5 triangle)
clip_polygon([(0,0),(4,0),(4,4),(0,4)], [(1,1),(3,1),(3,3),(1,3)])  # the 2x2 window
```

`segment_intersection` returns the single crossing point (or `None` for parallel,
collinear, or disjoint segments), while `segments_intersect` is the boolean test and also
catches shared endpoints and T-junctions. `polygon_perimeter` sums the edge lengths, and
`clip_polygon` is Sutherland-Hodgman clipping of a subject polygon against a convex
window — clipping the `[0,4]²` square by the `[1,3]²` window leaves exactly the `2×2`
interior (area `4`).

Extremal measures answer "how big is this cloud" — its span, box, and tightest circle:

```python
from quantforge import bounding_box, polygon_diameter, min_enclosing_circle

bounding_box([(1,2), (3,-1), (0,5), (4,4)])   # (0, -1, 4, 5)
polygon_diameter([(0,0), (1,0), (1,1), (0,1)])# ((0,0), (1,1), 1.414...) farthest pair
min_enclosing_circle([(0,0),(2,0),(2,2),(0,2)])# (1.0, 1.0, 1.414...) center + radius
```

`bounding_box` is the axis-aligned extent; `polygon_diameter` returns the farthest-apart
pair (searching only the convex-hull vertices, where the diameter always lives); and
`min_enclosing_circle` is Welzl's expected-linear algorithm for the smallest circle
covering every point — for the unit square it is centered at `(1,1)` with radius `√2`, and
it cannot be shrunk without leaving a point outside.

Any simple polygon can be cut into triangles for area, rendering, or meshing —
`ear_clipping_triangulate` does it, alongside orientation and convexity tests:

```python
from quantforge import ear_clipping_triangulate, signed_area, is_convex_polygon

L = [(0,0), (4,0), (4,2), (2,2), (2,4), (0,4)]   # concave L-shape
ear_clipping_triangulate(L)          # 4 triangles (n - 2)
signed_area([(0,0),(4,0),(4,4),(0,4)])# 16.0 (positive = counter-clockwise)
is_convex_polygon(L)                 # False
```

`ear_clipping_triangulate` snips off "ears" until a simple polygon (either winding) is
reduced to `n - 2` triangles whose areas sum to the polygon's own. `signed_area` is the
shoelace area — its sign gives the winding, which `is_clockwise` reports — and
`is_convex_polygon` checks that every turn goes the same way.

Point-proximity queries — hit-testing, snapping, route distance — come from the
point-to-line and point-to-segment distances:

```python
from quantforge import (point_to_line_distance, point_segment_distance,
                        closest_point_on_segment, point_polyline_distance)

point_to_line_distance((0, 3), (0, 0), (1, 0))       # 3.0 — perpendicular to the x-axis
point_segment_distance((5, 0), (0, 0), (3, 0))       # 2.0 — clamps past the endpoint
point_polyline_distance((3, 1), [(0,0),(2,0),(2,2),(4,2)])  # 1.0 — nearest segment
```

`point_to_line_distance` is the perpendicular distance to an *infinite* line;
`point_segment_distance` clamps to a finite segment (so a point past an endpoint measures
to that endpoint), with `closest_point_on_segment` returning the foot itself; and
`point_polyline_distance` takes the minimum over a chain of segments.

Thinning a dense polyline — a GPS track, a map contour — while keeping its shape is what
`douglas_peucker` does:

```python
from quantforge import douglas_peucker

douglas_peucker([(i, 0) for i in range(10)], 0.01)     # [(0,0), (9,0)] — a straight run
douglas_peucker([(0,0),(1,0),(2,5),(3,0),(4,0)], 1)    # [(0,0),(2,5),(4,0)] — keeps the peak
```

It keeps the endpoints, recursively retains the farthest vertex from the current chord
while its perpendicular distance exceeds `epsilon`, and drops the rest — so the output is
a subsequence of the input within tolerance of the original. A larger `epsilon` keeps
fewer points (monotonically), and a collinear run collapses to its two endpoints.

Circle constructions and intersections round out the geometry set:

```python
from quantforge import (circle_from_3points, circle_line_intersection,
                        circle_circle_intersection, point_in_circle)

circle_from_3points((1,0), (0,1), (-1,0))           # (0.0, 0.0, 1.0) — the circumcircle
circle_line_intersection((0,0,5), (-10,0), (10,0))  # [(5,0), (-5,0)]
circle_circle_intersection((0,0,5), (8,0,5))        # [(4,3), (4,-3)]
```

`circle_from_3points` returns the circle `(cx, cy, r)` through three non-collinear points
(collinear input raises). `circle_line_intersection` and `circle_circle_intersection`
return 0, 1 (tangent), or 2 points — every one lying on the circle(s) — and
`point_in_circle` is the inside test.

## Graph algorithms

The core graph routines on a dict adjacency list — shortest paths, traversal, components,
and dependency ordering:

```python
from quantforge import (dijkstra, shortest_path, bfs,
                        connected_components, topological_sort)

g = {'A': {'B': 1, 'C': 4}, 'B': {'C': 2, 'D': 5}, 'C': {'D': 1}, 'D': {}}
shortest_path(g, 'A', 'D')                       # (['A','B','C','D'], 4.0)
bfs({'A': ['B','C'], 'B': ['D'], 'C': ['D'], 'D': []}, 'A')  # {'A':0,'B':1,'C':1,'D':2}
connected_components({1: [2], 2: [1], 3: [4], 4: [], 5: []}) # [[1,2],[3,4],[5]]
topological_sort({'shirt': ['tie','belt'], 'tie': ['jacket'],
                 'belt': ['jacket'], 'jacket': []})          # a valid dress order
```

`dijkstra` returns the shortest-path distances and predecessors from a source over
non-negative weights (verified against Bellman-Ford on random graphs), and `shortest_path`
reconstructs the actual route and its length. `bfs` gives unweighted hop distances,
`connected_components` partitions an undirected graph, and `topological_sort` (Kahn's
algorithm) orders a DAG so every edge points forward — raising if a cycle makes that
impossible. Graphs are plain dicts: `{node: {neighbor: weight}}` for weighted routines,
`{node: [neighbors]}` for the unweighted ones.

On a DAG that same forward order makes path optimization a single linear sweep —
`dag_longest_path`/`dag_shortest_path` handle *negative* weights (unlike Dijkstra) with no
relaxation rounds, and the longest path is the CPM critical path:

```python
from quantforge import dag_longest_path, dag_shortest_path, transitive_closure

g = {0: [(1, 3), (2, 2)], 1: [(3, 4)], 2: [(3, 1)], 3: [(4, 2)], 4: []}
dag_longest_path(g, 0, 4)    # (9.0, [0, 1, 3, 4]) — the critical path
dag_shortest_path(g, 0, 4)   # (5.0, [0, 2, 3, 4])
transitive_closure({'a': ['b'], 'b': ['c'], 'c': []})   # {'a': {'b','c'}, 'b': {'c'}, 'c': set()}
```

These take a weighted DAG as `{node: [(neighbor, weight), ...]}`. `dag_shortest_path`/
`dag_longest_path` relax edges once in topological order, so they run in `O(V + E)` and,
because there is no cycle, tolerate negative weights that would break Dijkstra; passing a
`target` returns `(distance, path)` with the route reconstructed. `transitive_closure`
returns each node's reachable set (a DFS from each, cycles allowed). All raise on a cyclic
graph where acyclicity is required. Cross-checked against exhaustive path enumeration over
thousands of random DAGs.

Structure problems on weighted graphs — spanning trees and network capacity — round out
the toolkit:

```python
from quantforge import UnionFind, minimum_spanning_tree, max_flow

edges = [('A','B',1), ('A','C',4), ('B','C',2), ('B','D',5), ('C','D',1)]
minimum_spanning_tree(['A','B','C','D'], edges)   # (tree_edges, 4.0)

g = {'s': {'v1':16,'v2':13}, 'v1': {'v2':10,'v3':12}, 'v2': {'v1':4,'v4':14},
     'v3': {'v2':9,'t':20}, 'v4': {'v3':7,'t':4}, 't': {}}
max_flow(g, 's', 't')                             # 23 (the CLRS classic)
```

`minimum_spanning_tree` runs Kruskal's algorithm on a `UnionFind` disjoint-set forest and
returns the chosen edges plus total weight (verified against a brute-force minimum over
all spanning trees) — a minimum spanning forest for a disconnected graph. `max_flow` is
Edmonds-Karp (BFS-augmenting Ford-Fulkerson) on a `{node: {neighbor: capacity}}` network,
reproducing the textbook value of `23`. `UnionFind` is exposed directly for clustering and
connectivity work, with path compression and union by rank for near-constant operations.

When edges carry a *cost* as well as a capacity, `min_cost_max_flow` finds the maximum flow
that is cheapest among all maximum flows:

```python
from quantforge import min_cost_max_flow, MinCostMaxFlow

edges = [(0, 1, 2, 1), (0, 2, 2, 3), (1, 3, 2, 1), (2, 3, 2, 1), (1, 2, 1, 1)]
min_cost_max_flow(edges, 0, 3)    # (4, 12) — (max flow, minimum cost)

m = MinCostMaxFlow()
m.add_edge("s", "a", 3, 1).add_edge("a", "t", 3, 2)
m.solve("s", "t")                 # (3, 9)
```

Edges are `(u, v, capacity, cost)`. It augments along a *shortest-cost* residual path each
round (found by Bellman-Ford/SPFA, so negative-cost edges are fine) and pushes the path's
bottleneck, until the sink is unreachable — returning the total flow and its least cost.
`MinCostMaxFlow` is the incremental builder; `min_cost_max_flow` the one-shot convenience.
The plain `max_flow` is the all-equal-cost special case, and the two agree on flow value.
Cross-checked against brute integer-flow enumeration over thousands of small networks.

Pairing up two sides of a graph — jobs to workers, slots to applicants — is the bipartite
matching problem, and `maximum_bipartite_matching` solves it optimally:

```python
from quantforge import (maximum_bipartite_matching, maximum_matching_size,
                        minimum_vertex_cover)

adj = {'a': ['x', 'y'], 'b': ['x'], 'c': ['y', 'z']}   # left -> allowed rights
maximum_bipartite_matching(adj)     # {'a': 'y', 'b': 'x', 'c': 'z'} — all three paired
maximum_matching_size(adj)          # 3
minimum_vertex_cover(adj)           # ({'a', 'b', 'c'}, set()) — size 3, by Konig
```

`maximum_bipartite_matching` runs Hopcroft-Karp — repeatedly BFS-layering the graph and
DFS-augmenting along disjoint shortest paths, `O(E·sqrt(V))` — and returns the pairing as a
`{left: right}` dict (`maximum_matching_size` is just its length). `minimum_vertex_cover`
uses Konig's theorem: in a bipartite graph the smallest set of vertices touching every edge
has exactly the size of the maximum matching, recovered here from alternating-path
reachability. Verified against a brute augmenting-path matching over 4000 random graphs,
with the cover checked to match the matching size and cover every edge.

Centrality measures score how "important" each node is — by link structure, degree, or
position on shortest paths:

```python
from quantforge import (pagerank, degree_centrality,
                        closeness_centrality, betweenness_centrality)

pagerank({'A': ['B','C'], 'B': ['C'], 'C': ['A']})   # {'A':0.388,'B':0.215,'C':0.397}
star = {0: [1,2,3,4], 1: [0], 2: [0], 3: [0], 4: [0]}
degree_centrality(star)[0]                            # 1.0 — the hub touches everyone
betweenness_centrality({0:[1],1:[0,2],2:[1,3],3:[2,4],4:[3]})  # midpoint node 2 peaks
```

`pagerank` is the random-surfer stationary distribution (power iteration with
teleportation), summing to 1 and handling dangling nodes by redistributing their mass —
uniform on a symmetric ring. `degree_centrality` normalizes neighbor counts,
`closeness_centrality` rewards short distances to everyone, and `betweenness_centrality`
(Brandes' algorithm) counts how often a node lies on shortest paths — the hub of a star
dominates all three, and the middle of a path has the highest betweenness (`0.667`).

Dijkstra assumes non-negative weights and a single source; three more shortest-path
routines cover the rest:

```python
from quantforge import bellman_ford, floyd_warshall, a_star

g = {'A': {'B': 4, 'C': 5}, 'B': {'C': -3}, 'C': {'D': 2}, 'D': {}}
bellman_ford(g, 'A')[0]                    # {'A':0,'B':4,'C':1,'D':3} — negative edge OK
floyd_warshall(g)['A']['D']                # 3.0, all-pairs distances
a_star(grid, start, goal, heuristic)       # heuristic-guided point-to-point
```

`bellman_ford` handles negative edge weights and raises on a reachable negative cycle
(and agrees with Dijkstra when all weights are non-negative). `floyd_warshall` computes
every pairwise distance in `O(V³)`, matching per-source Dijkstra. `a_star` speeds up a
single-target search with an admissible heuristic — on a 5×5 grid with the Manhattan
heuristic it finds the optimal cost-`8` path, the same answer Dijkstra gives but exploring
fewer nodes, and it reduces exactly to Dijkstra when the heuristic is zero.

On a rooted tree, `LCA` answers lowest-common-ancestor and path-distance queries in
`O(log n)` after an `O(n log n)` build:

```python
from quantforge import LCA

#   a — b — {d, e},  a — c — f
adj = {'a': ['b','c'], 'b': ['a','d','e'], 'c': ['a','f'],
       'd': ['b'], 'e': ['b'], 'f': ['c']}
tree = LCA(adj, root='a')
tree.query('d', 'e')        # 'b' — deepest common ancestor
tree.query('d', 'f')        # 'a'
tree.distance('d', 'f')     # 4 edges on the path
tree.is_ancestor('a', 'f')  # True
```

`LCA` uses binary lifting: the build precomputes each node's `2^j`-th ancestor, so a query
lifts the deeper node to its partner's depth then lifts both in powers of two until they
meet. `depth`, `distance` (edges between any two nodes, via
`depth(u) + depth(v) - 2·depth(lca)`), and `is_ancestor` follow directly. The tree is built
with an iterative DFS, so even a chain of thousands of nodes will not overflow recursion.
Verified against a brute path-to-root reference over 2000 random trees.

On a *directed* graph, `strongly_connected_components` groups nodes that are mutually
reachable, and `condensation` collapses each group to a single vertex to expose the DAG
underneath:

```python
from quantforge import (strongly_connected_components, condensation,
                        is_strongly_connected)

g = {'a': ['b'], 'b': ['c'], 'c': ['a'], 'd': ['a']}   # cycle a->b->c->a, plus d->a
strongly_connected_components(g)   # [['c', 'b', 'a'], ['d']]
is_strongly_connected(g)           # False
component_of, dag = condensation(g)
dag                                # {0: [], 1: [0]} — d's component points at the cycle
```

`strongly_connected_components` is Tarjan's algorithm — one `O(V + E)` DFS tracking each
node's discovery index and the lowest index reachable from its subtree, emitting a
component whenever a node's low-link equals its own index. Components come out in reverse
topological order of the condensation (each before the ones it can reach). `condensation`
returns a `component_of` map and the condensed graph, which is always acyclic, with
duplicate edges and self-loops removed; `is_strongly_connected` and `number_of_sccs` are
shortcuts. The DFS is iterative, so large graphs will not overflow recursion. Verified
against a brute mutual-reachability reference over 3000 random digraphs.

Strongly connected components are exactly what makes 2-satisfiability tractable, so `TwoSat`
builds directly on them: it solves a boolean formula of two-literal clauses `(a OR b)` in
linear time:

```python
from quantforge import TwoSat

ts = TwoSat(3)                     # variables 0, 1, 2
ts.add_or(0, 1)                    # (x0 OR x1)     — literal v is positive, ~v is negated
ts.add_or(~0, 2)                   # (NOT x0 OR x2)
ts.add_or(~1, ~2)                  # (NOT x1 OR NOT x2)
ts.solve()                         # [True, False, True] — a satisfying assignment

bad = TwoSat(1)
bad.force_true(0); bad.force_true(~0)
bad.solve()                        # None — x and NOT x cannot both hold
```

Each clause `(a OR b)` is equivalent to the implications `NOT a -> b` and `NOT b -> a`, so
the formula becomes an implication graph over the `2n` literals. Running
`strongly_connected_components` on it, the formula is satisfiable exactly when no variable
lands in the same component as its own negation, and a satisfying assignment falls out of
the component order. `add_or` (alias `add_clause`) takes literals as `v` or `~v`, with
`add_implication` and `force_true` for readability; `solve` returns the assignment or
`None`, and `is_satisfiable` just checks. Verified against brute-force assignment
enumeration over 5000 random instances.

To trace a route that uses every *edge* exactly once — an Eulerian trail — `eulerian_path`
runs Hierholzer's algorithm, and the `has_eulerian_*` predicates check whether one exists:

```python
from quantforge import eulerian_path, has_eulerian_path, has_eulerian_circuit

g = {0: [1, 3], 1: [2], 2: [3]}     # undirected: edges 0-1, 0-3, 1-2, 2-3
eulerian_path(g)                    # [0, 1, 2, 3, 0] — a circuit (starts and ends at 0)
has_eulerian_circuit(g, False)      # True

gd = {'a': ['b'], 'b': ['c'], 'c': ['a']}
eulerian_path(gd, directed=True)    # ['a', 'b', 'c', 'a']
```

`eulerian_path` returns the `E + 1` vertex sequence of a trail using every edge once (first
and last coincide for a circuit), or `None` if none exists. Existence is a degree
condition, which the predicates test directly: an undirected connected graph needs zero
odd-degree vertices for a circuit or exactly two for an open path; a directed graph needs
in-degree equal to out-degree everywhere, or off by one at a single source/sink pair. Pass
`directed=True` for digraphs and an optional `start` vertex. Verified over 6000 random
graphs that existence matches the criterion and every returned trail uses each edge exactly
once.

## String algorithms

Sequence comparison and pattern matching — on strings or any lists:

```python
from quantforge import (levenshtein, hamming_distance,
                        longest_common_subsequence, longest_common_substring, kmp_search)

levenshtein("kitten", "sitting")                 # 3 edits
longest_common_subsequence("ABCBDAB", "BDCAB")   # "BCAB" (not contiguous)
longest_common_substring("abcdxyz", "xyzabcd")   # "abcd" (contiguous)
kmp_search("abababab", "abab")                   # [0, 2, 4] — overlapping matches
hamming_distance("karolin", "kathrin")           # 3
```

`levenshtein` is the minimum insert/delete/substitute edit distance (symmetric, verified
against a brute-force recursion), and `hamming_distance` its equal-length special case.
`longest_common_subsequence` finds the longest order-preserving but not necessarily
contiguous match, while `longest_common_substring` requires contiguity.  `kmp_search`
returns every start index of a pattern in linear time via the Knuth-Morris-Pratt failure
function, catching overlapping occurrences that a stride-by-length scan would miss.

Plain edit distance charges two edits for a swapped pair; `damerau_levenshtein` counts an
adjacent transposition as one, and `needleman_wunsch` gives the classic global alignment:

```python
from quantforge import damerau_levenshtein, needleman_wunsch

damerau_levenshtein("ca", "ac")            # 1 (Levenshtein would say 2)
damerau_levenshtein("kitten", "sitting")   # 3
needleman_wunsch("GATTACA", "GCATGCU")     # (score, 'G-ATTACA', 'GCA-TGCU')
```

`damerau_levenshtein` is the true unrestricted distance (verified never to exceed
Levenshtein and to match a brute-force BFS over edit operations), which models real typos
better since a transposition is one keystroke slip. `needleman_wunsch` maximizes a
match/mismatch/gap score and reconstructs the aligned, gapped strings — the score always
equals the alignment it returns, and stripping the gaps recovers the inputs.

`smith_waterman` is the *local* counterpart: it floors scores at zero and traces back from
the highest-scoring cell, so it finds the best-matching *substring* pair rather than
aligning end to end — `smith_waterman("xxxHELLOyyy", "zzHELLOww")` recovers the shared
`HELLO` (score 10). No positively-scoring common region gives score 0.

Where edit distance counts *operations*, fuzzy similarity scores return a value in
`[0, 1]` for approximate matching, deduplication, and record linkage:

```python
from quantforge import jaro, jaro_winkler, dice_coefficient, jaccard_similarity

jaro("MARTHA", "MARHTA")                       # 0.944
jaro_winkler("MARTHA", "MARHTA")               # 0.961 — shared-prefix boost
dice_coefficient("night", "nacht")             # 0.25 — bigram overlap
jaccard_similarity("the cat sat", "the dog sat")  # 0.5 — word-set overlap
```

`jaro` matches characters within a sliding window and penalizes transpositions —
reproducing Winkler's published reference values — and `jaro_winkler` rewards a shared
prefix, both tuned for short strings like names. `dice_coefficient` compares character
bigram multisets (robust to small edits and word order), and `jaccard_similarity` is the
set-overlap ratio over whitespace tokens by default, or any granularity via a `tokenize`
callable (`list` for characters). All are `1` for identical inputs and `0` for wholly
dissimilar ones.

A few classic sequence algorithms round out the section — on numbers or any comparable
list:

```python
from quantforge import longest_increasing_subsequence, maximum_subarray, longest_run

longest_increasing_subsequence([10, 9, 2, 5, 3, 7, 101, 18])  # [2, 3, 7, 18]
maximum_subarray([-2, 1, -3, 4, -1, 2, 1, -5, 4])             # (6, 3, 6) sum, start, end
longest_run([1, 1, 2, 3, 3, 3, 1])                             # (3, 3, 3) value, len, start
```

`longest_increasing_subsequence` uses patience sorting (`O(n log n)`) and reconstructs an
actual subsequence (`strict=False` for non-decreasing); `maximum_subarray` is Kadane's
algorithm returning the best contiguous sum and its bounds (correct even for all-negative
input); and `longest_run` finds the longest streak of a repeated value. Verified against
brute-force references over thousands of random arrays.

`kmp_search` finds one pattern; when you need to match *many* patterns at once,
`AhoCorasick` builds an automaton over the whole set and locates every occurrence in a
single pass over the text:

```python
from quantforge import AhoCorasick

ac = AhoCorasick(["he", "she", "his", "hers"])
sorted(ac.find_all("ushers"))          # [(3, 'he'), (3, 'she'), (5, 'hers')]

ac = AhoCorasick(["cat", "dog"])
ac.count_matches("catdogcat")          # 3
ac.contains_any("the dog ran")         # True
```

Each match is `(end_index, pattern)` — the index of the match's last character. The
automaton is a trie of the patterns with *failure links* (the longest proper suffix that
is also a trie prefix) and *output links*, so a longer match also reports every shorter
pattern ending at the same spot (`"she"` at index 3 also yields `"he"`), and the whole
scan runs in `O(len(text) + total_pattern_length + matches)` regardless of how many
patterns there are. `find_all` returns the list, `find` is a lazy generator, and
`contains_any` short-circuits at the first hit. Verified against a brute-force per-pattern
search over 2000 random pattern-set/text cases.

To *index* a single text for repeated substring queries, `suffix_array` sorts all its
suffixes, and `lcp_array` (Kasai) records how much adjacent sorted suffixes share — the two
structures behind fast search, longest-repeat, and distinct-substring counting:

```python
from quantforge import (suffix_array, lcp_array, substring_search,
                        longest_repeated_substring, count_distinct_substrings)

suffix_array("banana")               # [5, 3, 1, 0, 4, 2]
lcp_array("banana")                  # [0, 1, 3, 0, 0, 2]
substring_search("banana", "ana")    # [1, 3] — all occurrences, sorted
longest_repeated_substring("banana") # 'ana'
count_distinct_substrings("banana")  # 15
```

`suffix_array` is built by prefix doubling in `O(n log n)`; `substring_search` then binary-
searches it for a pattern's occurrence block in `O(m log n)`, returning every start index.
`longest_repeated_substring` reads off the deepest entry of the LCP array (the longest
prefix two suffixes share, hence a substring that occurs at least twice, overlaps
included), and `count_distinct_substrings` uses `sum(n - sa[r]) - sum(lcp)` — total suffix
length minus the prefixes already counted. All cross-checked against brute-force references
(sorting suffixes, naive LCP, the substring set) over thousands of random strings.

To index a whole *set* of strings for prefix lookups — autocomplete, dictionary membership
— `Trie` stores them in a prefix tree, so every query costs only its own length:

```python
from quantforge import Trie

t = Trie(["cat", "car", "card", "dog"])
t.keys_with_prefix("ca")      # ['car', 'card', 'cat']
t.count_prefix("car")         # 2 — 'car' and 'card'
t.starts_with("do")           # True
t.longest_prefix_of("cards")  # 'card' — longest stored key that prefixes the query
t.contains("ca")              # False — a prefix, not a stored key
```

`insert`, `contains`, and `delete` manage membership (delete prunes branches that no
longer lead to a key), while `starts_with`/`count_prefix`/`keys_with_prefix` answer prefix
questions in time proportional to the prefix — `count_prefix` is `O(len(prefix))` thanks to
a cached subtree count at each node — regardless of how many keys the trie holds.
`longest_prefix_of` walks the query one character at a time and remembers the deepest stored
key seen, the operation behind routing-table and tokenizer lookups. Verified against a plain
set/list reference over 3000 random word sets with interleaved deletes.

Three linear-time primitives cover pattern search, periodicity, and palindromes. The
`z_function` and `prefix_function` (the KMP border array) are the workhorses; on top of them
sit `count_occurrences`, `smallest_period`/`is_periodic`, `borders`, and the standalone
`manacher_longest_palindrome`:

```python
from quantforge import (z_function, prefix_function, smallest_period,
                        borders, count_occurrences, manacher_longest_palindrome)

prefix_function("abcabcd")             # [0, 0, 0, 1, 2, 3, 0]
smallest_period("abcabcabc")           # 3 — repeats the unit "abc"
borders("abacaba")                     # [1, 3] — "a" and "aba" are prefix == suffix
count_occurrences("ababab", "ab")      # [0, 2, 4] — overlaps included
manacher_longest_palindrome("babad")   # 'bab'
```

`z_function` gives, at each index, the length of the longest substring there that matches a
prefix; `prefix_function` gives the longest proper prefix that is also a suffix. From the
prefix function, `smallest_period` reads `n - pi[-1]` (the shortest block the string
repeats, whole or partial), `is_periodic` checks whether that block tiles it exactly, and
`borders` walks the failure links to list every prefix==suffix length. `count_occurrences`
finds all matches of a pattern (including overlaps) in linear time via the Z-function over
`pattern\0text`, and `manacher_longest_palindrome` returns a longest palindromic substring.
All cross-checked against brute-force references over 4000 random strings each.

A de Bruijn sequence packs *every* length-`n` string over a `k`-symbol alphabet into one
cyclic string of length `k^n`, each appearing exactly once — the shortest such string.
`de_bruijn_sequence` builds it, and `lyndon_words`/`is_lyndon` expose the Lyndon words it is
assembled from:

```python
from quantforge import de_bruijn_sequence, lyndon_words, is_lyndon

de_bruijn_sequence(2, 3)      # [0,0,0,1,0,1,1,1] — every 3-bit string appears once cyclically
de_bruijn_sequence(2, 2)      # [0,0,1,1]
lyndon_words(2, 3)            # [(0,), (0,0,1), (0,1), (0,1,1), (1,)]
is_lyndon([0,0,1]), is_lyndon([0,1,0,1])   # (True, False)
```

`de_bruijn_sequence(k, n)` uses the Fredricksen-Kessler-Maiorana construction — concatenate,
in lexicographic order, the Lyndon words whose length divides `n` — generated by Duval's
constant-amortized-time algorithm. The result is the minimal-length string containing all
`k^n` windows, the classic tool for PIN-pad cracking and rotary-encoder tracks. `is_lyndon`
tests whether a string is strictly smaller than all its rotations, and `lyndon_words` lists
them all up to a length. Verified that every `n`-tuple occurs exactly once cyclically and the
Lyndon words match a brute rotation test.

## Data compression

Two foundational lossless codes, both exact round-trips:

```python
from quantforge import (huffman_encode, huffman_decode,
                        run_length_encode, run_length_decode)

bits, codebook = huffman_encode("mississippi")   # 21 bits (vs 22 fixed-width)
huffman_decode(bits, codebook)                    # back to the original symbols
run_length_encode("aaabbc")                       # [('a',3), ('b',2), ('c',1)]
```

`huffman_encode` builds the optimal prefix code for the symbol distribution — the
shortest expected code length of any prefix code — and returns the bitstring plus the
codebook needed to decode. Its expected length lands within one bit of the Shannon
entropy and never exceeds fixed-width coding; frequent symbols get the shortest codes
(`s` and `i` in "mississippi" get 1- and 2-bit codes). `run_length_encode` collapses runs
of a repeated symbol into `(symbol, count)` pairs — ideal for sparse or blocky data — and
`run_length_decode` expands them back. Both invert losslessly.

The Burrows-Wheeler transform and move-to-front coding are the front of the bzip2
pipeline — a reversible permutation that clusters repeats, then a code that turns those
clusters into mostly-zero runs:

```python
from quantforge import bwt_transform, bwt_inverse, move_to_front_encode

bwt_transform("banana")                        # ('nnbaaa', 3)
bwt_inverse("nnbaaa", 3)                        # 'banana'
move_to_front_encode("aaaabbbbcccc")[0]         # [0,0,0,0,1,0,0,0,2,0,0,0]
```

`bwt_transform` returns the last column of the sorted-rotations matrix plus the primary
index needed to invert; `bwt_inverse` reconstructs the original exactly via the
LF-mapping (verified round-trip on hundreds of random strings). The transform makes runs
of a character adjacent, so `move_to_front_encode` — which replaces each symbol by its
rank in a running alphabet and moves it to the front — emits long runs of zeros that a
following RLE + Huffman pass compresses well. `move_to_front_decode` and `bwt_inverse`
chain back to the original losslessly.

LZW takes the opposite approach — it builds its dictionary of repeated substrings on the
fly, so nothing extra is transmitted, and delta coding shrinks slowly-varying numeric
series:

```python
from quantforge import lzw_compress, lzw_decompress, delta_encode, delta_decode

codes, alphabet = lzw_compress("abcabcabcabc")   # 7 codes for 12 characters
lzw_decompress(codes, alphabet)                   # 'abcabcabcabc'
delta_encode([100, 101, 102, 103, 104])           # [100, 1, 1, 1, 1]
```

`lzw_compress` emits the code of the longest known prefix and adds each new substring to
the table as it goes, so repetitive text needs far fewer codes than characters;
`lzw_decompress` rebuilds the same table (handling the code-equals-next-entry edge case)
and inverts exactly. `delta_encode` replaces a numeric sequence by its first value and
successive differences — a slowly-varying series becomes small numbers a following
entropy coder packs tightly — and `delta_decode` is the cumulative sum back.

The same byte streams need integrity checks; three standard non-cryptographic hashes cover
that:

```python
from quantforge import crc32, adler32, fnv1a_32

crc32("hello")        # 907060870 — matches zlib.crc32
adler32("Wikipedia")  # 300286872
fnv1a_32("foobar")    # 3214735720 — table/bloom-filter mixer
```

`crc32` implements the IEEE 802.3 reflected polynomial (the checksum in Ethernet, gzip,
and PNG) and matches `zlib.crc32` bit for bit; `adler32` is the cheaper Zlib checksum,
also matching `zlib.adler32`; and `fnv1a_32` is a fast, well-dispersed hash for hash
tables and bloom filters, reproducing the published FNV test vectors. All accept a `str`
(UTF-8 encoded) or raw bytes and return a 32-bit unsigned integer — none are
cryptographic.

Checksums only *detect* corruption; error-correcting codes fix it. The Hamming(7,4) code
adds three parity bits so any single-bit flip is corrected, and the Luhn algorithm guards
account/card numbers with a check digit:

```python
from quantforge import hamming74_encode, hamming74_decode, luhn_checksum, luhn_check_digit

code = hamming74_encode([1, 0, 1, 1])     # 7-bit codeword
code[4] ^= 1                               # flip a bit in transit
hamming74_decode(code)                     # ([1, 0, 1, 1], 5) — data recovered, flip at pos 5

luhn_checksum("79927398713")               # 0 -> valid
luhn_check_digit("7992739871")             # 3 -> the digit that validates it
```

`hamming74_decode` computes the 3-bit syndrome, corrects the flagged bit, and returns the
data plus the 1-indexed error position (`0` if clean) — verified to correct a flip at every
one of the seven positions for all sixteen data words. `luhn_checksum` returns `0` for a
valid number (it validates real card numbers) and `luhn_check_digit` produces the digit
that makes any payload valid.

## Probability distributions

Gamma, chi-square, Poisson, F and binomial distributions built on those special
functions — CDFs, densities or masses, and quantiles:

```python
from quantforge import (gamma_cdf, gamma_ppf, chi2_cdf, chi2_ppf, chi2_sf,
                        poisson_cdf, poisson_pmf, f_ppf, binomial_cdf)

chi2_ppf(0.95, df=10)            # -> 18.307  (textbook critical value)
chi2_sf(18.307, df=10)           # -> 0.05    (upper-tail p-value)
f_ppf(0.95, d1=1, d2=10)         # -> 4.965
poisson_cdf(k=7, lam=4.3)        # P(N <= 7), via the gamma identity
binomial_cdf(k=6, n=20, prob=0.3)
gamma_ppf(0.5, shape=2.5, scale=3.0)   # gamma median
```

The Poisson and binomial CDFs use the incomplete-gamma and incomplete-beta
identities (`P(N<=k) = Q(k+1, lam)`, `P(X<=k) = I_{1-p}(n-k, k+1)`), so they stay
exact and stable at large parameters where summing masses would lose precision.
Every continuous quantile inverts its CDF and round-trips to machine precision.

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

The `expectile` is the only risk measure that is both coherent and elicitable —
VaR is elicitable but not coherent, ES coherent but not elicitable:

```python
from quantforge import expectile

expectile(pnl, tau=0.5)      # = mean loss
expectile(pnl, tau=0.95)     # tail-weighted; monotone in tau
```

It solves the asymmetric-least-squares condition on the loss variable, so unlike a
quantile it is backtestable by a single scoring function.

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

The measures above discount at the risk-free rate. The Moody's-KMV distance to
default and expected default frequency use the firm's real-world asset drift `mu`
instead:

```python
from quantforge import physical_distance_to_default, physical_default_probability

physical_distance_to_default(asset_value=120, debt_face=100, mu=0.10, sigma=0.25, t=1)
physical_default_probability(120, 100, mu=0.10, sigma=0.25, t=1)   # 0.158
```

Setting `mu = r` recovers the risk-neutral figures; because a risky firm earns
`mu > r`, its physical default probability sits below the risk-neutral one.

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
from quantforge import (fx_forward, forward_points, implied_base_rate,
                        cross_rate, triangular_arbitrage)

fx_forward(spot=1.10, r_price=0.05, r_base=0.03, t=1.0)   # EURUSD-style
cross_rate(1.10, 1.25)                                     # EURGBP from EURUSD, GBPUSD
triangular_arbitrage(1.10, 150, 1 / (1.10 * 150))         # == 1 if arbitrage-free
```

## Entropy pooling (views on scenarios)

Impose a view on a set of scenarios — "the mean of this quantity should be 3" —
while changing the scenario probabilities as little as possible.
`entropy_pooling_mean` returns the minimal-relative-entropy posterior (an
exponential tilt of the prior, Meucci 2008); `relative_entropy` measures the KL
cost of the view.

```python
from quantforge import entropy_pooling_mean, relative_entropy

x = [-3, -1, 0, 1, 2, 4, 5, 7]          # scenario values, uniform prior
p = entropy_pooling_mean(x, target=3.0)  # posterior probabilities
sum(x[i] * p[i] for i in range(len(x)))  # -> 3.0  (view satisfied exactly)

relative_entropy(p, [1/8] * 8)           # -> 0.0657  (KL cost of the view)
```

A view equal to the prior mean leaves the probabilities unchanged; a stronger view
costs more relative entropy. The reweighted scenarios flow straight into any
downstream risk or allocation calculation.

## Covariance shrinkage (Ledoit-Wolf)

The sample covariance is noisy when the number of assets is not small relative to
the number of observations — a poor input to a mean-variance optimizer.
`ledoit_wolf_shrinkage` blends it with a constant-correlation target by the
data-driven Ledoit-Wolf (2004) intensity `delta`, which minimizes the expected
distance to the true covariance.

```python
from quantforge import ledoit_wolf_shrinkage

# rows = observations, columns = assets.
returns = [
    [0.01, 0.02, -0.01],
    [-0.02, -0.01, 0.03],
    [0.015, 0.01, -0.02],
    [0.00, -0.015, 0.01],
    [-0.005, 0.005, 0.00],
]
sigma, delta = ledoit_wolf_shrinkage(returns)
delta                                   # -> 0.2008  (shrink 20% toward the target)
sigma[0][0]                             # -> 0.00015 (asset-0 variance, preserved)
```

`delta` sits in `[0, 1]` and falls toward zero as the sample grows (the sample
estimate becomes reliable); the result keeps each asset's own variance, pulls the
correlations toward their average, and is symmetric positive definite. Feed
`sigma` straight into the mean-variance optimizer below.

An alternative to shrinkage is random-matrix-theory denoising, which cleans only the
noise eigenvalues instead of pulling the whole matrix toward a target:

```python
from quantforge import marchenko_pastur_edge, clip_correlation_eigenvalues

marchenko_pastur_edge(n_assets=50, n_obs=250)          # noise-eigenvalue cutoff
clip_correlation_eigenvalues(correlation, n_obs=250)   # clip the noise bulk
```

Eigenvalues below the Marchenko-Pastur edge `(1 + sqrt(N/T))^2` are consistent with
pure noise; `clip_correlation_eigenvalues` replaces them with their average while
keeping the signal eigenvalues, so the trace is preserved and the matrix stays a
correlation matrix — better-conditioned for optimization.

For a *time-varying* estimate that weights recent data more heavily, use the
exponentially-weighted (RiskMetrics) covariance:

```python
from quantforge import ewma_covariance_matrix, ewma_correlation_matrix

ewma_covariance_matrix(returns, lam=0.94)     # RiskMetrics decay
ewma_correlation_matrix(returns, lam=0.94)    # unit diagonal, entries in [-1, 1]
```

Both are symmetric positive definite; the EWMA tracks regime shifts a shrinkage
of the full-sample covariance would smooth over.

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

# Arbitrary risk budgets (generalizes risk parity) and the realized contributions.
from quantforge import risk_budget_weights, risk_contributions
w = risk_budget_weights(cov, budgets=[0.5, 0.3, 0.2])   # target risk shares
risk_contributions(w, cov)        # each asset's variance contribution (sums to var)

# Clustering-based and inverse-vol allocation (no matrix inversion).
from quantforge import inverse_volatility_weights, hierarchical_risk_parity
inverse_volatility_weights(cov)
hierarchical_risk_parity(cov)     # Lopez de Prado HRP
```

Measure how diversified a set of weights actually is — nominal count vs independent
risk sources:

```python
from quantforge import (herfindahl_index, effective_number_of_constituents,
                        effective_number_of_bets)

herfindahl_index(weights)                    # sum of squared weights
effective_number_of_constituents(weights)    # 1/HHI: equal-weight-equivalent count
effective_number_of_bets(weights, cov)        # Meucci: independent risk sources
```

`effective_number_of_constituents` counts positions by weight; Meucci's
`effective_number_of_bets` counts *uncorrelated* risk factors, so a portfolio of
many correlated names can hold dozens of constituents but only a handful of real
bets.

## Rebalancing

Weight drift, turnover, transaction-cost drag, and no-trade-band rebalancing:

```python
from quantforge import (drift_weights, turnover, transaction_cost,
                        no_trade_band_rebalance)

drifted = drift_weights([0.5, 0.5], asset_returns=[0.20, -0.10])
turnover(drifted, target_weights=[0.5, 0.5])              # fraction traded
transaction_cost(drifted, [0.5, 0.5], cost_bps=10)        # round-trip drag
no_trade_band_rebalance(drifted, [0.5, 0.5], band=0.05)   # only trade if outside band
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

## Leveraged ETFs

Daily-rebalanced leveraged/inverse ETF path and the volatility (compounding) drag:

```python
from quantforge import volatility_drag, flat_market_decay, leveraged_etf_path

volatility_drag(leverage=3, sigma=0.20)          # per-period drag vs naive 3x
leveraged_etf_path(daily_returns, leverage=3, expense_ratio=0.0095)
flat_market_decay(3, [0.10, 1/1.1 - 1])          # loses value over a net-flat path
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

More drawdown-adjusted return ratios, each dividing return by a different pain
denominator:

```python
from quantforge import gain_to_pain_ratio, sterling_ratio, burke_ratio

gain_to_pain_ratio(returns)   # sum(returns) / sum(|losses|); inf if never loses
sterling_ratio(returns)       # excess return / (average drawdown + 10% margin)
burke_ratio(returns)          # excess return / sqrt(sum of squared drawdowns)
```

`sterling_ratio` uses the L1 average drawdown; `burke_ratio` the L2 root-sum-of-
squares, so Burke penalizes a few deep drawdowns more heavily than many shallow
ones. `gain_to_pain_ratio` is the scale-free profitability-vs-pain summary.

Depth is only half the drawdown story — `drawdown_analytics` adds duration:

```python
from quantforge import drawdown_analytics

a = drawdown_analytics(returns)
a["max_drawdown_depth"]    # deepest peak-to-trough fractional drop
a["time_to_recovery"]      # periods from trough back to the prior peak (None if never)
a["longest_underwater"]    # most periods spent below a prior high
```

It also returns the `peak_index`, `trough_index`, and `recovery_index` of the
deepest episode (indices into the equity curve). A series that ends underwater
reports `recovery_index=None`; the longest underwater stretch need not coincide
with the deepest drawdown.

For drawdown-based *risk* measures — the drawdown analogues of VaR and expected
shortfall — use the drawdown-at-risk family:

```python
from quantforge import (average_drawdown, drawdown_at_risk,
                        conditional_drawdown_at_risk)

average_drawdown(returns)                              # mean drawdown over the path
drawdown_at_risk(returns, confidence=0.95)             # DaR: the drawdown quantile
conditional_drawdown_at_risk(returns, confidence=0.95) # CDaR: mean of the worst tail
```

`CDaR >= DaR >= average_drawdown >= 0`, all capped by the maximum drawdown. CDaR is
coherent (the drawdown analogue of expected shortfall) and is the objective in
drawdown-constrained portfolio optimization.

For downside-risk-adjusted return, `kappa_ratio` generalizes Sortino and
Omega-Sharpe into one family (order 1 = Omega-Sharpe, order 2 = Sortino), and
`upside_potential_ratio` rewards the chance of beating a target per unit of
downside:

```python
from quantforge import kappa_ratio, upside_potential_ratio

kappa_ratio(returns, tau=0.0, order=2)     # per-period Sortino
kappa_ratio(returns, tau=0.0, order=3)     # penalizes deep shortfalls more
upside_potential_ratio(returns, tau=0.0)   # upside expectation / downside deviation
```

A higher `order` weights deep shortfalls more heavily, so Kappa falls as the order
rises for a left-skewed series.

Against a benchmark, the CAPM measures decompose risk-adjusted return by systematic
exposure:

```python
from quantforge import market_beta, treynor_ratio, jensens_alpha, m_squared

market_beta(returns, market)                 # Cov/Var slope on the market
treynor_ratio(returns, market)               # excess return per unit of beta
jensens_alpha(returns, market)               # return beyond what beta explains
m_squared(returns, market)                   # portfolio rescaled to market risk
```

`treynor_ratio` divides by beta (systematic risk) where Sharpe divides by total
volatility; `jensens_alpha` is the CAPM regression intercept (positive = skill);
`m_squared` restates the Sharpe ranking in return units directly comparable to the
market.

For tail risk that accounts for non-normal returns, the Cornish-Fisher pair adjusts
the Gaussian quantile with the sample skewness and excess kurtosis:

```python
from quantforge import cornish_fisher_var, cornish_fisher_expected_shortfall

r = [0.01, -0.02, 0.015, -0.03, 0.008, 0.02, -0.05, 0.012, -0.01, 0.006, 0.018, -0.04]
cornish_fisher_var(r, confidence=0.95)                  # 0.05, skew/kurt-adjusted VaR
cornish_fisher_expected_shortfall(r, confidence=0.95)   # 0.0582, tail mean beyond it
```

Both reduce to the Gaussian VaR/ES for a normal series; negative skew and fat tails
push each above its Gaussian value, and the expected shortfall never falls below the
VaR.

## GARCH volatility

Fit GARCH(1,1), forecast the term volatility, and price consistently with the
vol mean-reversion:

```python
from quantforge import fit_garch, garch_term_variance, garch_option_price

params = fit_garch(returns)
garch_option_price(params, last_return, last_variance,
                   S=100, K=100, r=0.05, option_type="call", horizon=21)

# Asymmetric (leverage) variants: GJR-GARCH and log-variance EGARCH.
from quantforge import (GJRGarchParams, gjr_garch_forecast,
                        EGarchParams, egarch_forecast)
gjr = GJRGarchParams(omega=1e-6, alpha=0.05, beta=0.90, gamma=0.06)
gjr_garch_forecast(gjr, last_return=-0.02, last_variance=4e-4, horizon=5)
eg = EGarchParams(omega=-0.1, alpha=0.15, beta=0.95, gamma=-0.08)
egarch_forecast(eg, last_return=-0.02, last_variance=4e-4, horizon=5)
```

Before fitting a GARCH model, check that the returns actually show ARCH effects —
Engle's LM test regresses the squared returns on their own lags:

```python
from quantforge import arch_lm_test

lm, p = arch_lm_test(returns, lags=5)   # small p: volatility clustering present
```

A small p-value means the variance depends on the recent past, so a GARCH model is
warranted; iid returns are not rejected.

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
