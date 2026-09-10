# Changelog

All notable changes to QuantForge are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/), and the project adheres to
[Semantic Versioning](https://semver.org/).

## [1.17.0] - 2026-09-09

### Added
- `bjerksund_stensland_greeks`: delta, gamma, vega, theta (calendar), and rho of
  the Bjerksund-Stensland American price by central finite differences (the
  2002 closed form has no simple Greek expressions). For a no-dividend American
  call the Greeks equal the European BSM Greeks, as they must.

## [1.16.1] - 2026-09-09

### Changed
- `calibrate_svi` now uses a deterministic multi-start (several fixed seeds,
  keep the best) so it no longer stalls in the degenerate huge-`b` valley of
  the raw-SVI objective. On the example chain the worst per-expiry fit improved
  from rmse ~9e-4 to ~5e-6 and the assembled surface is calendar-arbitrage free.

### Added
- `examples/vol_surface.py`: end-to-end surface workflow (invert quotes -> fit
  SVI per expiry -> assemble `VolSurface` -> calendar check -> interpolate vol
  and extract Dupire local vol), with a smoke test.

## [1.16.0] - 2026-09-09

### Added
- `localvol.py`: Dupire local volatility. `dupire_local_vol` evaluates the
  Dupire formula from a call-price surface `C(K, T)` by finite differences;
  `local_vol_from_implied` wraps an implied-vol surface via Black-Scholes.
  Recovers a flat implied vol as a constant local vol and matches the analytic
  term-structure local variance `dw/dT`.

## [1.15.0] - 2026-09-09

### Added
- `strategy.py`: multi-leg option-strategy builders returning a `Book` (so net
  price/Greeks come from the existing engine): `vertical_spread`, `straddle`,
  `strangle`, `risk_reversal`, `butterfly`, `iron_condor`. Plus
  `payoff_at_expiry` / `payoff_profile` for the P&L diagram and `break_evens`
  (grid scan + bisection) for the zero-P&L spots.

## [1.14.0] - 2026-09-09

### Added
- `arithmetic_asian` (in `exotics.py`): closed-form arithmetic-average Asian via
  Turnbull-Wakeman moment matching (match the average's first two moments to a
  lognormal, then Black-Scholes). Complements the Monte Carlo pricer and agrees
  with it to a few cents; arithmetic value dominates the geometric Asian.

## [1.13.0] - 2026-09-09

### Added
- `multiasset.py`: two-asset options. `exchange_option` (Margrabe, exact),
  `spread_option` (Kirk approximation, puts via parity), and `basket_option`
  (Levy lognormal moment-match on the two-asset weighted sum). Exchange and
  basket match correlated-GBM Monte Carlo; basket single-asset reduces to
  Black-Scholes exactly.

## [1.12.0] - 2026-09-09

### Added
- `varswap.py`: model-free variance- and volatility-swap fair strikes by static
  option replication (Demeterfi-Derman-Kamani-Zou log-strip). `variance_swap_strike`
  integrates an OTM put/call strip weighted by `1/K^2` around the forward;
  `volatility_swap_strike` returns the `sqrt` proxy. Recovers `sigma^2` from a
  flat-vol Black-Scholes chain to strip-truncation error.

## [1.11.0] - 2026-09-09

### Added
- `density.py`: Breeden-Litzenberger risk-neutral density extraction from a
  call-price curve. `risk_neutral_density` (second strike-derivative, non-
  uniform grid), `risk_neutral_cdf` (first derivative), `price_from_density`
  (integrate any payoff against the recovered density), and
  `density_total_mass` (sanity check ~1). Recovers the lognormal pdf from a
  BSM curve and reprices vanilla and digital payoffs.

## [1.10.0] - 2026-09-09

### Added
- `merton.py`: `merton_jump_price` prices European options under Merton (1976)
  jump-diffusion as a Poisson-weighted sum of Black-Scholes prices with
  jump-adjusted volatility and carry. Collapses to Black-Scholes at zero jump
  intensity; the compensated drift keeps the forward a martingale so put-call
  parity holds exactly. Cross-checked against Monte Carlo.

## [1.9.0] - 2026-09-09

### Added
- `hedgesim.py`: `simulate_delta_hedge` Monte Carlos a discretely delta-hedged
  short option, returning the hedging-error distribution (mean, std, min, max).
  Supports hedging at a different vol than the realized path (`hedge_vol` vs
  `real_vol`) to study vol-mismatch P&L.
- Test suite gains a `slow` marker; run `pytest -m "not slow"` for a ~4s fast
  pass (deep Monte Carlo / tree cross-checks are marked slow).

## [1.8.0] - 2026-09-09

### Added
- `bachelier.py`: the Bachelier (normal) model — `bachelier_price`, analytic
  delta/gamma/vega, and `bachelier_implied_vol` (Newton + bisection). Prices
  options on a forward following arithmetic Brownian motion, so it handles
  negative forwards/strikes (rates and spread options) where the lognormal
  model breaks down. `sigma` here is the normal (absolute) volatility.

## [1.7.0] - 2026-09-09

### Added
- `heston.py`: `heston_price` prices European options under the Heston (1993)
  stochastic-volatility model via its characteristic function (Albrecher
  "little trap" form), integrated with a self-contained 64-point
  Gauss-Legendre rule (no SciPy). Puts follow from parity. Collapses to
  Black-Scholes as the vol-of-vol goes to zero.

## [1.6.0] - 2026-09-09

### Added
- `lookback.py`: continuously-monitored lookback options —
  `floating_strike_lookback` (Goldman-Sosin-Gatto; payoff against the realized
  extreme) and `fixed_strike_lookback` (Conze-Viswanathan; ordinary strike on
  the realized extreme). Both take the running extreme and cost of carry `b`,
  with the `b -> 0` singularity handled by a nudge. Cross-checked against Monte
  Carlo path max/min.

## [1.5.0] - 2026-09-09

### Added
- `hedging.py`: smile-aware delta. `smile_delta` returns the effective delta
  under the sticky-strike rule (equals BS delta) or the sticky-delta /
  sticky-moneyness rule (adds a `-vega * (dsigma/dk) / S` skew term).
  `skew_slope` finite-differences a supplied smile, and
  `smile_delta_from_smile` wires the two together.

## [1.4.0] - 2026-09-09

### Added
- `forwardstart.py`: `forward_start_price` prices forward-start options
  (strike fixed at a future date as a multiple of the then-spot) via
  Rubinstein's closed form, and `cliquet_price` values a cliquet/ratchet as a
  strip of consecutive forward-starts.

## [1.3.0] - 2026-09-09

### Added
- `forward.py`: `implied_forward` extracts the implied forward price and
  discount factor from a call/put chain via a put-call-parity least-squares
  fit (no volatility assumption), and backs out the implied rate and dividend
  yield.

## [1.2.0] - 2026-09-09

### Added
- `trinomial.py`: `trinomial_price` prices American/European options on a Boyle
  trinomial lattice (smoother convergence than the binomial tree), and
  `richardson_american` combines `n`/`2n` solves to cancel the leading O(1/n)
  error for a more accurate American price.

## [1.1.0] - 2026-09-09

### Added
- `surface.py`: `VolSurface` stitches per-expiry SVI smiles into a term
  structure, interpolates total variance linearly in maturity, and reports
  calendar arbitrage (total variance must be non-decreasing in `t` at each
  strike). `VolSurface.fit` calibrates one SVI slice per expiry.

## [1.0.0] - 2026-09-09

First stable release. The public API is now considered stable under SemVer.

### Added
- Single-sourced version: `pyproject.toml` reads `quantforge.__version__`.

### Summary of the 1.0 feature set
- **Pricing** — generalized Black-Scholes-Merton (stock, dividend, Black-76,
  FX via cost-of-carry `b`); European and American exercise.
- **Greeks** — analytic delta, gamma, vega, theta, rho, plus the second-order
  vanna, vomma/volga, charm, veta, speed, zomma, color.
- **Implied volatility** — robust Newton-with-bisection solver with
  arbitrage-band rejection.
- **American options** — Cox-Ross-Rubinstein binomial tree and the
  Bjerksund-Stensland (2002) closed form (~500x faster than the tree).
- **Exotics** — cash/asset-or-nothing digitals, single-barrier options
  (Reiner-Rubinstein, all four kinds, with rebate), geometric-Asian.
- **Volatility surfaces** — Gatheral raw SVI and SABR (Hagan expansion), each
  with calibration.
- **Monte Carlo** — GBM engine with antithetic and control-variate variance
  reduction; standard errors on every estimate.
- **Portfolio** — batch pricing, net Greeks, VaR / Expected Shortfall
  (parametric delta-gamma, historical, full-reprice MC), spot×vol stress grid.
- **Realized volatility** — close-to-close, EWMA, Parkinson, Garman-Klass,
  Rogers-Satchell, Yang-Zhang estimators.
- **Performance** — zero-dependency core (~1.4M prices/sec) with an optional
  NumPy vectorized fast path (~6x on large batches).
- **Tooling** — 182 tests, CI on Python 3.8/3.10/3.12, tag-triggered PyPI
  release via Trusted Publishing, reproducible benchmark suite.

## [0.1.0 - 0.13.0]

Pre-1.0 development. Each minor version added one major capability area:
BSM core (0.1), portfolio/CLI/CI (0.2), SVI surface (0.3), exotics (0.4),
Monte Carlo (0.5), VaR/ES (0.6), realized vol (0.7), second-order Greeks and
the release workflow (0.8), benchmarks (0.9), the NumPy fast path (0.10), the
stress grid (0.11), Bjerksund-Stensland American (0.12), and SABR (0.13).
