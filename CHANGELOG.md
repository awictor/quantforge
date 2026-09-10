# Changelog

All notable changes to QuantForge are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/), and the project adheres to
[Semantic Versioning](https://semver.org/).

## [1.64.0] - 2026-09-09

### Added
- `lookback_greeks` (in `lookback.py`): delta, gamma, vega, and theta of a
  floating- or fixed-strike lookback by central finite differences on the
  closed form (`kind="floating"`/`"fixed"`, with the running extreme).

## [1.63.0] - 2026-09-09

### Added
- `ratio_spread` and `backspread` (in `strategy.py`): a ratio spread (long 1,
  short `ratio` at a further strike; net short options) and its mirror the
  backspread (short 1, long `ratio`; net long options). Both return a leg
  `Book`, so net price/Greeks and the payoff diagram come from the engine.

## [1.62.0] - 2026-09-09

### Added
- GARCH(1,1) volatility forecasting (in `volatility.py`): `fit_garch` estimates
  the variance model by Gaussian quasi-MLE (stationary reparametrization,
  Nelder-Mead), and `garch_forecast` gives the annualized vol `horizon` steps
  ahead, mean-reverting to the long-run level. Recovers persistence on synthetic
  GARCH data.

## [1.61.0] - 2026-09-09

### Added
- `compo_option` (in `quanto.py`): composite (compo) FX option on a foreign
  asset converted at the *floating* exchange rate, so the effective vol combines
  the asset and FX vols with their correlation
  (`sqrt(sa^2 + sfx^2 + 2 rho sa sfx)`). Complements the fixed-FX quanto; higher
  correlation raises the price, and it matches a combined-lognormal Monte Carlo.

## [1.60.0] - 2026-09-09

### Added
- `variancegamma.py`: `variance_gamma_price` prices European options under the
  Variance-Gamma (Madan-Carr-Chang) pure-jump model via its characteristic
  function, integrated with the shared Gauss-Legendre quadrature (no SciPy).
  `nu` controls kurtosis and `theta` the skew; `nu -> 0` recovers Black-Scholes.
  Monte-Carlo verified.

## [1.59.0] - 2026-09-09

### Added
- `average_strike_asian_mc` (in `montecarlo.py`): Monte Carlo an average-strike
  Asian option, where the strike is the realized arithmetic average of the path
  (call pays `max(S_T - A, 0)`). Complements the fixed-strike Asians; a single
  monitoring date collapses the payoff to zero.

## [1.58.0] - 2026-09-09

### Added
- `attribution.py`: `attribute_pnl` explains an option position's realized P&L
  over a market move via a second-order Greek expansion (delta / gamma / vega /
  theta / rho P&L), returning a `PnLAttribution` with the true revaluation, the
  explained sum, and the unexplained residual. The residual is tiny for modest
  moves and grows for large ones (higher-order Greeks).

## [1.57.0] - 2026-09-09

### Changed
- `implied_volatility` now seeds the Newton solve with the Corrado-Miller (1996)
  rational approximation instead of the ATM-only Brenner-Subrahmanyam guess.
  It is accurate away from the money too, cutting the Newton iteration count
  several-fold (~5x fewer across a strike/vol grid), and falls back to the old
  seed in the deep wings. The solver's answers are unchanged.

## [1.56.0] - 2026-09-09

### Added
- `local_vol_mc` (in `montecarlo.py`): Monte Carlo a European option under a
  Dupire local-volatility surface `sigma_loc(S, t)`, evolving the spot in
  log-space with a spot/time-dependent vol at each step. A flat local vol
  reproduces the Black-Scholes price; pairs with `dupire_local_vol` /
  `sabr_local_vol` for surface-consistent pricing.

## [1.55.0] - 2026-09-09

### Added
- `parisian_barrier_mc` (in `montecarlo.py`): Monte Carlo for a Parisian barrier
  option, which activates only after the spot stays past the barrier for a
  *consecutive* window (robust to brief spikes) rather than on a single touch.
  Supports all four in/out, up/down kinds; knock-in + knock-out equals the
  vanilla, and a Parisian knock-out is worth more than the instantaneous one.

## [1.54.0] - 2026-09-09

### Added
- `vol_cone` (in `volatility.py`): the realized-volatility cone — for each
  rolling window length it reports the min / 25th / median / 75th / max and the
  current realized vol (annualized), returning `VolConePoint` per window. The
  cone narrows as the window grows and the median tracks the true vol.

## [1.53.0] - 2026-09-09

### Added
- SVI butterfly-arbitrage check (in `svi.py`): `svi_g` evaluates the
  Gatheral-Jacquier g-function of a slice (>= 0 everywhere iff no butterfly /
  density arbitrage), with `svi_butterfly_arbitrage` listing violating strikes
  and `svi_is_butterfly_free` the boolean. Cross-checked against the sign of the
  Breeden-Litzenberger density.

## [1.52.0] - 2026-09-09

### Added
- `theta_carry_report` (in `bookgreeks.py`): decomposes a book's net theta into
  the gamma-rent term (`-0.5 * Gamma * sigma^2 * S^2`) and a residual
  drift/financing carry, returning a `ThetaCarry`. With zero carry the theta is
  pure gamma rent; a non-zero rate produces the financing residual.

## [1.51.0] - 2026-09-09

### Added
- `vannavolga.py`: `VannaVolgaSmile` builds an FX smile from the three market
  quotes (ATM vol, 25-delta risk reversal, 25-delta butterfly), recovering the
  pillar vols/strikes and interpolating the vol at any strike. `pillar_vols`
  exposes the 25P/ATM/25C vols. Exact at the three pillars.

## [1.50.0] - 2026-09-09

### Added
- Kelly-criterion sizing (in `sizing.py`): `kelly_fraction_binary` (optimal
  stake for a binary bet from win probability and odds),
  `kelly_fraction_continuous` (growth-optimal leverage `mu / sigma^2`, with a
  fractional-Kelly multiplier), and `kelly_growth_rate` (expected log-growth at
  a given leverage, maximized at full Kelly).

## [1.49.0] - 2026-09-09

### Added
- `early_exercise_premium` (in `american.py`): decomposes the American price
  (Bjerksund-Stensland) into the European (BSM) value plus the early-exercise
  premium. The premium is zero for a no-dividend American call and positive for
  ITM puts and dividend-paying calls.

## [1.48.0] - 2026-09-09

### Added
- `displaced.py`: `displaced_diffusion_price` prices under Rubinstein's
  displaced-diffusion (shifted-lognormal) model — a Black-Scholes price on
  `S + shift` / `K + shift` with a rescaled vol. `shift = 0` recovers
  Black-Scholes; a positive shift allows negative strikes/spot and moves the
  skew toward normal-model behavior. Monte-Carlo verified.

## [1.47.0] - 2026-09-09

### Added
- `dividend_curve` (in `forward.py`): bootstraps an implied dividend-yield term
  structure from a multi-expiry option chain, applying `implied_forward` per
  expiry and returning sorted `(t, ForwardResult)` points. Recovers a known
  dividend term structure exactly.

## [1.46.0] - 2026-09-09

### Added
- `quanto.py`: `quanto_option` prices a foreign-asset option settled in domestic
  currency at a fixed exchange rate. The quanto adjustment shifts the carry by
  `-rho * sigma_asset * sigma_fx` and discounts at the domestic rate; `rho=0`
  removes it. Monte-Carlo verified.

## [1.45.0] - 2026-09-09

### Added
- `perpetual.py`: `perpetual_american` prices a no-expiry American option in
  exact closed form (Merton 1973), and `perpetual_exercise_boundary` returns the
  flat optimal-exercise spot. Matches the binomial American price at long
  maturity; a call with carry >= rate is never exercised early.

## [1.44.0] - 2026-09-09

### Added
- `gap_option` and `power_option` (in `exotics.py`): closed forms for a gap
  option (separate trigger and payoff strikes; Reiner-Rubinstein) and a power
  option (payoff on `S^power`; adjusted-drift/vol Black-Scholes). Equal gap
  strikes and `power=1` recover the vanilla option; the power form is
  Monte-Carlo verified.

## [1.43.0] - 2026-09-09

### Added
- `correlation_term_structure` (in `correlation.py`): implied correlation at each
  expiry from index and member vol term structures, returning
  `(expiry, rho)` pairs. Recovers a constant or maturity-varying correlation
  exactly.

## [1.42.0] - 2026-09-09

### Added
- `asian_greeks` (in `exotics.py`): delta, gamma, vega, and theta of a Asian
  option by central finite differences on either closed form
  (`average="geometric"` Kemna-Vorst or `"arithmetic"` Turnbull-Wakeman). The
  arithmetic-average delta exceeds the geometric one, matching the price order.

## [1.41.0] - 2026-09-09

### Added
- `compound.py`: `compound_option` prices Geske (1979) compound options — an
  option on an option — for all four kinds (call/put-on-call/put), reusing the
  bivariate-normal CDF and solving the critical-spot exercise boundary by
  bisection. Cross-checked against Monte Carlo; call-on-call collapses to the
  vanilla as the first strike goes to zero.

## [1.40.0] - 2026-09-09

### Added
- `best_of_call` / `worst_of_call` (in `multiasset.py`): rainbow options on the
  max / min of two correlated assets, priced by Monte Carlo on correlated GBM
  (antithetic, seeded). Best-of + worst-of equals the sum of the two single-name
  calls (Stulz identity), which the tests verify.

## [1.39.0] - 2026-09-09

### Added
- `chooser.py`: `chooser_option` prices a simple chooser (Rubinstein 1991) in
  closed form — the holder picks call or put at a future date. Decomposes into a
  call to expiry plus a put on the discounted-forward strike expiring at the
  choice date; equals a straddle when the choice is at expiry.

## [1.38.0] - 2026-09-09

### Added
- `barrier_digital_mc` (in `montecarlo.py`): Monte Carlo for a cash-or-nothing
  digital contingent on a barrier condition (up/down, knock-in/knock-out) — pays
  the cash only if the option finishes in the money AND the barrier condition
  holds over the path. Knock-in + knock-out sums to the plain digital.

## [1.37.0] - 2026-09-09

### Added
- `sabr_local_vol` (in `localvol.py`): Dupire local volatility of a single SABR
  smile — builds the Hagan implied-vol smile on the forward and feeds it through
  the Dupire formula. A flat (nu->0, beta=1) SABR gives a constant local vol =
  alpha; a skewed smile gives the steeper-than-implied local skew.

## [1.36.0] - 2026-09-09

### Added
- `vegabucket.py`: `vega_buckets` groups a multi-expiry book's position-scaled
  vega into maturity buckets defined by upper-edge tenors, returning a
  `VegaBuckets` whose buckets sum to the net book vega — so a desk can see where
  its vol risk sits along the curve.

## [1.35.0] - 2026-09-09

### Added
- `lsm.py`: `bermudan_lsm` prices Bermudan/American options by Longstaff-Schwartz
  least-squares Monte Carlo — backward induction over exercise dates, regressing
  the discounted continuation value on a polynomial basis of spot (normal
  equations solved in pure Python). Converges to the binomial American value as
  the number of exercise dates grows.

## [1.34.0] - 2026-09-09

### Added
- `gramcharlier.py`: Corrado-Su skew/kurtosis-adjusted pricing. `corrado_su_call`
  / `corrado_su_price` add the first skewness and excess-kurtosis corrections to
  Black-Scholes via a Gram-Charlier expansion (skew=kurt=0 recovers BSM;
  kurtosis fattens the tails). `realized_skewness` and
  `realized_excess_kurtosis` estimate those moments from a return series.

## [1.33.0] - 2026-09-09

### Added
- `capped_cliquet_mc` (in `montecarlo.py`): Monte Carlo pricer for a locally-
  and globally-capped cliquet (ratchet) note — sums clipped periodic returns
  and clips the running total, with antithetic variates and a seed. Tightening
  either cap lowers the price; the global cap bounds the payoff.

## [1.32.0] - 2026-09-09

### Added
- `sizing.py`: hedge-quantity helpers. `delta_hedge_shares` zeros a book's net
  delta with the underlying; `neutralize` solves the units of a hedge option to
  move delta/gamma/vega to a target; `vega_neutral_quantity` and
  `gamma_neutral_quantity` are the zero-target shortcuts. All work off the net
  Greeks from `price_book`.

## [1.31.0] - 2026-09-09

### Added
- `VolSurface.forward_variance` / `forward_vol`: the forward (instantaneous-
  average) variance and volatility between two maturities, from the additive
  total-variance surface — the vol of a forward-starting option. Raises on a
  negative forward variance (calendar arbitrage).

## [1.30.0] - 2026-09-09

### Added
- `barrier_greeks` (in `exotics.py`): delta, gamma, vega, and theta (calendar)
  of a single-barrier option by central finite differences on the
  Reiner-Rubinstein price. A far knock-out matches the vanilla Greeks, and
  knock-in + knock-out delta equals the vanilla delta (in-out parity).

## [1.29.0] - 2026-09-09

### Added
- `cev.py`: Constant-Elasticity-of-Variance pricing (Schroder/Hull) for
  `0 <= beta < 1`, with a from-scratch noncentral chi-square CDF
  (`noncentral_chisq_cdf`) and regularized incomplete gamma. The chi-square
  summation starts at the Poisson mode so it is stable at large noncentrality.
  `cev_price` matches Black-Scholes when scaled, satisfies parity, and produces
  the leverage skew (lower beta -> richer downside puts). Cross-checked vs
  Monte Carlo.

## [1.28.0] - 2026-09-09

### Added
- `correlation.py`: index implied correlation. `implied_correlation` inverts
  the index-variance decomposition for the single common correlation consistent
  with a quoted index vol; `index_vol_from_correlation` is the forward map; and
  `dispersion_basket_vol` is the zero-correlation reference. The standard
  dispersion-trading measure.

## [1.27.0] - 2026-09-09

### Added
- `qmc.py`: quasi-Monte Carlo. `halton` generates low-discrepancy points (van
  der Corput radical inverse per prime base) and `european_qmc` prices a
  European option by deterministic Halton integration of the payoff — it
  converges to the Black-Scholes value several times faster than pseudo-random
  Monte Carlo at the same point count.

## [1.26.0] - 2026-09-09

### Added
- `book_second_order` (in `bookgreeks.py`): position-scaled net second-order
  Greeks across a book — vanna, vomma/volga, charm, veta, speed, zomma, color —
  returned as a `BookSecondOrder`. Complements the first-order net Greeks from
  `price_book`.

## [1.25.0] - 2026-09-09

### Added
- Documentation site. `docs/gen_api.py` generates `docs/api.md` by introspecting
  `quantforge.__all__` (signatures + docstrings, grouped by module) with zero
  dependencies; `docs/index.md` and `mkdocs.yml` wire up an mkdocs site.
- CI now runs `docs/gen_api.py --check`, and a test asserts the reference is in
  sync and every public callable/class is documented, so the docs can't drift.

## [1.24.0] - 2026-09-09

### Added
- `epsilon` (dividend rho) in `bsm.py`: analytic sensitivity of the option
  price to the continuous dividend yield, `dPrice/dq`. Negative for calls,
  positive for puts, and verified against finite differences.

## [1.23.0] - 2026-09-09

### Added
- `dv01.py`: `key_rate_dv01` computes bucketed (key-rate) DV01 for any book
  expressed as `price(zero_curve)`. Bumps each tenor independently (central or
  one-sided), reports per-bucket sensitivities normalized to 1bp plus the
  parallel DV01; the buckets sum to the parallel shift.

## [1.22.0] - 2026-09-09

### Added
- `overhedge.py`: super-replicate a cash-or-nothing digital with a tight
  vanilla spread. `digital_call_overhedge` / `digital_put_overhedge` return an
  `Overhedge` (spread cost as a conservative price, the fair digital value, and
  the cushion between them); `overhedge_payoff` gives the spread's terminal
  payoff, which dominates the digital everywhere and converges to it as the
  spread width shrinks.

## [1.21.0] - 2026-09-09

### Added
- `swaption_price` (in `rates.py`): European payer/receiver swaptions on the
  Bachelier model — the swap's PV annuity times a normal-model option on the
  forward swap rate (handles negative rates). Plus `annuity` (PV01) and
  `swaption_parity` (payer - receiver = annuity * (swap_rate - strike)).

## [1.20.0] - 2026-09-09

### Added
- `rates.py`: interest-rate caps, floors, and collars priced as Bachelier
  (normal-vol) caplet/floorlet strips, so they handle negative rates.
  `CapletPeriod` describes each accrual period; `cap_price`/`floor_price` sum
  the strip, `collar_price` is long-cap/short-floor, and `caplet_floorlet_parity`
  gives the check identity. Cap - floor at one strike equals the swap PV.

## [1.19.0] - 2026-09-09

### Added
- `spline.py`: a pure-stdlib natural cubic spline (`CubicSpline`, Thomas-solved
  tridiagonal moments, C2-continuous, clamped outside its range) and
  `SmileSpline`, a strike->implied-vol interpolator with flat extrapolation — a
  model-free alternative to SVI/SABR for a single smile.

## [1.18.0] - 2026-09-09

### Added
- `one_touch` / `no_touch` (in `exotics.py`): continuously-monitored touch
  binaries. `one_touch` pays cash if the barrier is ever reached (immediately
  on hit, the FX convention, or deferred to expiry); `no_touch` pays if it
  never is. Direction (up/down barrier) is inferred from `H` vs `S`.
  Cross-checked against a barrier-crossing Monte Carlo; touch + no-touch (paid
  at expiry) sum to the discounted cash.

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
