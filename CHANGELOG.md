# Changelog

All notable changes to QuantForge are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/), and the project adheres to
[Semantic Versioning](https://semver.org/).

## [1.312.0] - 2026-09-11

### Added
- `sabr_normal_vol` (in `sabr.py`): the Hagan (2002) normal (Bachelier) implied
  volatility for the SABR model — the absolute-vol `sigma_N` quoted in
  interest-rate markets. Uses the Hagan normal expansion (leading factor
  `nu (F - K) / x(z)`, third-order time bracket shared with the lognormal
  `sabr_vol`), with the ATM `F == K` limit handled separately.
- Verified: matches the reference normal vol (SABR-Black vol -> Black-Scholes
  price -> Bachelier implied vol) to within 5e-3 across strikes; the ATM branch
  joins continuously with the just-off-ATM value; `beta = nu = 0` gives exactly
  `alpha`; the normal-vol smile is positive and rises in both wings.

## [1.311.0] - 2026-09-11

### Added
- `displaced_diffusion_smile` (in `displaced.py`): the Black-Scholes implied-vol
  smile a displaced-diffusion (shifted-lognormal) model produces. Prices a call
  at each strike and inverts to BS implied vol, returning `(log_moneyness, vol)`
  on the forward `F = S e^{b t}`. Mirrors `cev_smile` / `merton_smile`.
- Verified: `shift = 0` gives a flat smile at `sigma`; the ATM vol sits near
  `sigma`; a positive shift produces a downward skew, steeper for larger shift;
  output is sorted by strike; each smile vol reprices to the displaced-diffusion
  price it came from.

## [1.310.0] - 2026-09-11

### Added
- `cev_smile` (in `cev.py`): the Black-Scholes implied-vol smile a CEV model
  produces. Prices a call at each strike under CEV and inverts to BS implied
  vol, returning `(log_moneyness, vol)` pairs on the forward
  `F = S e^{(r-q) t}`. Mirrors `merton_smile` / `kou_smile`.
- Verified: the ATM implied vol sits near `sigma` (calibrated to the ATM
  instantaneous vol); the smile has a downward skew for `beta < 1`, steeper for
  smaller `beta`; output is sorted by strike; each smile vol reprices to the CEV
  price it came from.

## [1.309.0] - 2026-09-11

### Added
- `bjerksund_stensland_boundary` (in `american.py`): the Bjerksund-Stensland
  (2002) flat exercise trigger `I` at inception — a call is exercised for
  `S >= I`, a put for `S <= I`. Returns `None` when early exercise is never
  optimal (an American call with `b >= r`). The put trigger follows from the
  same put-call transformation the pricer uses: the transformed call's
  spot-axis trigger `I2t` maps back to `K^2 / I2t`.

### Changed
- Extracted the BS2002 trigger computation into a shared `_bs2002_triggers`
  helper used by both the call pricer and the new boundary function (prices
  unchanged).
- Verified: the call boundary lies above the strike and the put below; the
  BS2002 price equals the exercise intrinsic exactly at the boundary; a spot
  just inside the boundary has continuation value strictly above intrinsic; the
  boundary sits within ~10% of the Barone-Adesi-Whaley one (both flat
  approximations); the refactor leaves the reference prices unchanged.

## [1.308.0] - 2026-09-11

### Added
- `baw_critical_spot` (in `baw.py`): the Barone-Adesi-Whaley early-exercise
  boundary `S*` at inception — a call is exercised for `S >= S*`, a put for
  `S <= S*`. Returns `None` when early exercise is never optimal (an American
  call with `b >= r`).
- `baw_american_greeks` (in `baw.py`): delta, gamma, vega, theta of the BAW
  American price by finite difference.

### Changed
- Refactored the BAW critical-spot Newton solves into reusable `_critical_call`
  and `_critical_put` helpers shared by the pricer and the new boundary
  function (prices unchanged).
- Verified: the call boundary lies above the strike and the put boundary below;
  the American price equals the exercise intrinsic exactly at `S*`; greeks match
  central finite differences for calls and puts across several spots; the put
  delta lies in `(-1, 0)` with positive gamma; the BAW value dominates the
  European price.

## [1.307.0] - 2026-09-11

### Added
- `bachelier_theta` (in `bachelier.py`): analytic calendar theta of the normal
  (Bachelier) model, `theta = r * price - e^{-rt} sigma phi(d) / (2 sqrt(t))`.
  `bachelier_greeks` now uses it instead of a finite difference.
- `bachelier_cash_or_nothing` and `bachelier_asset_or_nothing` (in
  `bachelier.py`): normal-model digitals. `F_T` is Gaussian, so a cash-or-nothing
  call is `cash e^{-rt} N(d)` and an asset-or-nothing call is
  `e^{-rt} (F N(d) + sigma sqrt(t) phi(d))`, with `d = (F - K)/(sigma sqrt t)`.
- Verified: analytic theta matches central finite differences at, above, and
  below the strike for calls and puts; cash-digital call/put parity sums to the
  discount factor; the vanilla Bachelier price decomposes into asset-or-nothing
  minus `K` cash-or-nothing; both digitals match a Gaussian Monte Carlo.

## [1.306.0] - 2026-09-11

### Added
- `partial_time_end_barrier_call` (in `exotics.py`): partial-time (end)
  single-barrier call (Heynen-Kat 1994), where the down barrier is monitored
  only over `[t1, T2]` (inactive before `t1`). Exact bivariate-normal closed
  form coupling the monitoring-start date to expiry (`rho = sqrt(t1/T2)`).
  Down-out is priced directly; down-in follows from in-out parity. (Up-barrier
  partial-time calls have a distinct form and raise `ValueError`.)
- Verified: as `t1 -> 0` it approaches the standard continuously-monitored
  down-out barrier, and as `t1 -> T2` it approaches the vanilla call; a mid
  window sits strictly between the two; in-out parity holds exactly; matches a
  path Monte Carlo across `t1 in {0.25, 0.5, 0.75}` to within ~2%.

## [1.305.0] - 2026-09-11

### Added
- `discrete_barrier_option` (in `exotics.py`): discretely-monitored
  single-barrier option via the Broadie-Glasserman-Kou (1999) continuity
  correction. A barrier checked at `n_fixings` equally-spaced dates is breached
  less often than a continuous one, so the correction shifts the barrier away
  from the spot by `exp(+/- beta sigma sqrt(dt))` (up for up-barriers, down for
  down-barriers, `beta ~ 0.5826`) and prices with the continuous
  `barrier_option`. Applies to all four knock in/out types.
- Verified: matches a path Monte Carlo at `n = 50` across all four barrier
  types (down/up x in/out) to within ~3%; a discrete knock-out sits above and a
  knock-in below the continuous price; in-out parity holds exactly at the
  shifted barrier (KI + KO = vanilla); converges to the continuous barrier on a
  fine grid.

## [1.304.0] - 2026-09-11

### Added
- `discrete_fixed_strike_lookback` (in `lookback.py`): discretely-monitored
  fixed-strike lookback via the Broadie-Glasserman-Kou (1999) continuity
  correction. The realized extreme is sampled at `n_fixings` equally-spaced
  dates; the correction shifts the spot fed to the continuous
  `fixed_strike_lookback` by `exp(-/+ beta sigma sqrt(dt))` (down for a call on
  the max, up for a put on the min, `beta ~ 0.5826`). As `n_fixings` grows the
  shift vanishes and the price converges to the continuous lookback.
- Verified: matches a path Monte Carlo at `n = 50` (call/put) to within ~1.5%;
  the discrete call sits below the continuous lookback; the price converges
  (monotonically, and the correction decays like `1/sqrt(n)`) to the continuous
  value on a fine grid.

## [1.303.0] - 2026-09-11

### Added
- `holder_extendible_put` and `holder_extendible_put_greeks` (in
  `extendible.py`): holder-extendible put (Longstaff 1990), completing the full
  2x2 extendible family (holder/writer x call/put). At the first expiry `t1` the
  holder takes the best of exercising against `K1`, lapsing, or paying a fee `A`
  to extend to `T2` as a put struck at `K2`. The terminal spot splits into
  exercise (`S < I_low`), extend (`I_low <= S <= I_high`), and lapse
  (`S > I_high`); the extended-put strip is built from the call-strip pieces via
  `N(-d) = 1 - N(d)`. Greeks by finite difference.
- Verified: matches a Monte Carlo at the first expiry across three
  strike/maturity/fee cases to within ~0.03%; a very large fee collapses the
  strip and recovers the vanilla put to `t1`; a finite fee adds value; a lower
  fee is worth more; delta is negative and vega positive.

## [1.302.0] - 2026-09-11

### Added
- `writer_extendible_call` and `writer_extendible_call_greeks` (in
  `extendible.py`): writer-extendible call (Longstaff 1990), completing the
  extendible family. At the first expiry `t1` the call is exercised if in the
  money (`S_{t1} > K1`); otherwise the writer's obligation extends automatically
  to `T2` as a call struck at `K2` (no fee). Closed form: a vanilla call to `t1`
  plus the extended-call value collected on `S_{t1} <= K1`, via bivariate
  normals coupling `t1` and `T2`. Greeks by finite difference.
- Verified: matches a Monte Carlo at the first expiry across three
  strike/maturity cases to within ~0.1%; the automatic extension makes it worth
  more than a plain call to `t1`; delta is positive and vega positive.

## [1.301.0] - 2026-09-11

### Added
- `writer_extendible_put` and `writer_extendible_put_greeks` (in
  `extendible.py`): writer-extendible put (Longstaff 1990). At the first expiry
  `t1` the put is exercised if in the money (`S_{t1} < K1`); otherwise the
  writer's obligation extends automatically to `T2` as a put struck at `K2`
  (no fee). Closed form: a vanilla put to `t1` plus the extended-put value
  collected on `S_{t1} >= K1`, via bivariate normals coupling `t1` and `T2`.
  Greeks by finite difference.
- Verified: matches a Monte Carlo at the first expiry across three
  strike/maturity cases to within ~0.1%; the automatic extension makes it worth
  more than a plain put to `t1`; delta is negative and vega positive.

## [1.300.0] - 2026-09-11

### Added
- `holder_extendible_call` and `holder_extendible_call_greeks` (new module
  `extendible.py`): holder-extendible call (Longstaff 1990). At the first expiry
  `t1` the holder takes the best of exercising against `K1`, lapsing, or paying
  a fee `A` to extend the life to `T2` with strike `K2`. Priced in closed form
  by splitting the terminal spot into three regions (lapse / extend / exercise)
  and valuing the extension strip with bivariate normals coupling `t1` and `T2`;
  the strip boundaries are found by bisection. Greeks by finite difference.
- Verified: matches a Monte Carlo at the first expiry across three
  strike/maturity/fee cases to within ~0.15%; a very large fee collapses the
  extension strip and recovers the vanilla call to `t1`; a finite fee adds value
  over that plain call; a lower fee is worth more.

## [1.299.0] - 2026-09-11

### Added
- `complex_chooser_option` and `complex_chooser_option_greeks` (in
  `chooser.py`): complex chooser (Rubinstein 1991), where at the choice date the
  holder keeps whichever is worth more of a call (strike `Kc`, expiry `Tc`) or a
  put (strike `Kp`, expiry `Tp`) — the two legs may differ in both strike and
  maturity. Priced by Rubinstein's bivariate-normal formula, with the critical
  spot (where the two legs are equal at the choice date) found by bisection.
  Greeks by finite difference.
- Verified: equal strikes and maturities reduce exactly to the simple
  `chooser_option`; matches a Monte Carlo at the choice date to within ~0.14%;
  the chooser is worth more than either leg valued outright today; delta matches
  a finite difference.

## [1.298.0] - 2026-09-11

### Changed
- `supershare_greeks` now computes `delta` and `gamma` analytically instead of
  by finite difference. With price `A S (N(d1_lo) - N(d1_hi))`, `A = e^{(b-r)t}/K_low`,
  the spot derivatives are
  `delta = A (dN + dphi/(sigma sqrt t))` and
  `gamma = A/(S sigma sqrt t) (dphi - (d1_lo phi(d1_lo) - d1_hi phi(d1_hi))/(sigma sqrt t))`,
  where `dN`, `dphi` are the differences of the two `N`/`phi` at the corridor
  edges. `vega` and `theta` remain finite differences.
- Verified: analytic delta and gamma match central finite differences across
  four strike/vol/maturity cases.

## [1.297.0] - 2026-09-11

### Changed
- `gap_option_greeks` now computes `delta` and `gamma` analytically instead of
  by finite difference. Using the identity
  `S carry phi(d1) = K_trigger disc phi(d2)`, the spot sensitivities collapse to
  closed forms in `d1, d2` at the trigger, with a cash-driven correction term
  proportional to the strike gap `K_trigger - K_payoff`. `vega` and `theta`
  remain finite differences.
- Verified: analytic delta and gamma match central finite differences across
  three strike/vol/maturity cases for both calls and puts; gamma is identical
  for calls and puts; equal trigger/payoff strikes recover the vanilla
  Black-Scholes delta and gamma.

## [1.296.0] - 2026-09-11

### Added
- `range_binary_greeks` (in `exotics.py`): analytic `delta` and `gamma` for the
  range binary plus finite-difference `vega`/`theta`. Delta is
  `cash e^{-rt} (phi(d2_lo) - phi(d2_hi))/(S sigma sqrt t)`; gamma carries the
  double-sided pin risk near either corridor edge.
- `supershare_greeks` (in `exotics.py`): finite-difference `delta`, `gamma`,
  `vega`, `theta` for the supershare option.
- Verified: analytic range-binary delta and gamma match central finite
  differences across three strike/vol/maturity cases; both `price` fields match
  their pricers; a forward-centered corridor has positive calendar theta (the
  terminal mass concentrates inside as expiry nears).

## [1.295.0] - 2026-09-11

### Added
- `range_binary` (in `exotics.py`): range binary / double digital paying `cash`
  iff `K_low <= S_T <= K_high` at expiry. Exactly the difference of two
  cash-or-nothing calls: `cash e^{-rt} (N(d2(K_low)) - N(d2(K_high)))`.
- `supershare` (in `exotics.py`): Hakansson (1976) supershare paying
  `S_T / K_low` iff `K_low <= S_T <= K_high`. A scaled difference of two
  asset-or-nothing calls.
- Verified: each equals the corresponding digital difference to machine
  precision; both match a terminal Monte Carlo to within ~0.4%; a wide corridor
  range binary approaches the full discounted cash; the zero-vol limit pays iff
  the forward lands inside the corridor.

## [1.294.0] - 2026-09-11

### Added
- `implied_geometric_basket_correlation` (in `multiasset.py`): backs out the
  uniform pairwise correlation implied by a geometric-basket quote, assuming an
  equicorrelation matrix (`corr[i][j] = rho` off-diagonal). The basket
  log-variance rises with `rho`, so both call and put prices increase in `rho`
  (parity `C - P = disc(F - K)` is strike-independent), and a bisection over
  `[-1/(n-1), 1]` — the PSD-preserving range — recovers it.
- Verified: round-trips the correlation for calls and puts across
  `rho in {-0.3, 0, 0.4, 0.8}` and for a two-asset basket; raises when the quote
  lies outside the `rho`-range; requires at least two assets.

## [1.293.0] - 2026-09-11

### Added
- `geometric_basket_option` and `geometric_basket_greeks` (in `multiasset.py`):
  weighted geometric-average basket option on `prod_i S_i^{w_i}` for any number
  of assets. Unlike the arithmetic `basket_option` (a Levy moment-match), the
  weighted geometric average of correlated lognormals is itself lognormal, so
  this is exact: `log B` is Gaussian with variance
  `t sum_ij w_i w_j corr[i][j] sigma_i sigma_j`, priced by a Black-Scholes
  formula on the basket forward. Greeks return per-asset delta/gamma lists plus
  total vega and theta.
- Verified: a single unit-weight asset recovers the vanilla Black-Scholes
  price; put-call parity equals the discounted basket forward minus strike;
  matches a Cholesky Monte Carlo (3 correlated assets) to within ~0.2%;
  per-asset deltas and vega are positive.

## [1.292.0] - 2026-09-11

### Added
- `average_strike_arithmetic_asian` and `average_strike_arithmetic_asian_greeks`
  (in `exotics.py`): average-strike (floating-strike) discrete arithmetic Asian,
  payoff `max(S_T - A, 0)`. Priced as an exchange option between `S_T` (exactly
  lognormal) and the arithmetic average `A` (Levy two-moment lognormal). The
  cross-moment `E[S_T A] = (S^2/n) sum_i exp(b(t+t_i) + sigma^2 t_i)` is exact,
  giving the log-space covariance and the Margrabe spread variance. Greeks by
  finite difference.
- Verified: put-call parity equals the discounted exchange forward
  `E[S_T] - E[A]` exactly; the average-strike arithmetic call is cheaper than
  the geometric one (AM-GM raises the strike); matches an antithetic Monte
  Carlo (n=12) to ~1% (a small stable Levy approximation bias); explicit times
  match the `n_fixings` grid.

### Changed
- `average_strike_geometric_asian_greeks` test now asserts the exact
  homogeneity property (gamma zero, delta = price/S) instead of a fragile
  gamma-sign check, since the average-strike payoff is homogeneous of degree 1
  in spot.

## [1.291.0] - 2026-09-11

### Added
- `average_strike_geometric_asian` and `average_strike_geometric_asian_greeks`
  (in `exotics.py`): average-strike (floating-strike) discrete geometric Asian.
  The strike is the realized geometric average, so a call pays
  `max(S_T - G, 0)`. `S_T` and `G` are jointly lognormal, making this an
  exchange option with an exact Margrabe-style closed form on the spread
  variance `sigma^2 t + v_G - 2 sigma^2 mean(t_i)`. Greeks by finite
  difference.
- Verified: put-call parity equals the discounted exchange forward
  `E[S_T] - E[G]` exactly; matches an antithetic Monte Carlo (n=12, call/put)
  to within ~0.1%; explicit equally-spaced times match the `n_fixings` grid;
  gamma and vega are positive.

## [1.290.0] - 2026-09-11

### Added
- `seasoned_geometric_asian_greeks` and `seasoned_arithmetic_asian_greeks` (in
  `exotics.py`): finite-difference Greeks (`delta`, `gamma`, `vega`, `theta`)
  for the seasoned (in-progress) Asian pricers. The calendar bump scales the
  remaining fixing schedule with `t` so a `theta` difference stays within the
  option's life and preserves the fixing shape.
- Verified: with no observations the Greeks match the fresh
  `discrete_geometric_asian_greeks` / `discrete_arithmetic_asian_greeks`; the
  `price` field matches the pricer; partial-window delta/gamma/vega are
  positive and theta negative for a call; a partially-seasoned option has
  smaller delta and vega than a fresh one (less of the average is still
  stochastic).

## [1.289.0] - 2026-09-11

### Added
- `seasoned_arithmetic_asian` (in `exotics.py`): prices an in-progress discrete
  arithmetic-average Asian where some fixings are already observed. Writing the
  average as `A = (Q + sum remaining)/n` with `Q = sum(observed)` known, a call
  payoff is an arithmetic-average option on the remaining fixings with the
  shifted strike `K' = nK - Q`, scaled by `1/n`. The remaining sum's two exact
  moments are matched to a lognormal (Levy 1992). Deep-in/out branches (`K' <=
  0`) are priced exactly as `disc*(E[A] - K)` / zero.
- Verified: with no observations it reduces to `discrete_arithmetic_asian`;
  all-observed gives the deterministic arithmetic intrinsic; the `K' <= 0`
  branch equals the discounted expected-average payoff; a partial window
  matches a path Monte Carlo (n=12, 6 observed) to within 1%; the seasoned
  arithmetic price sits above the seasoned geometric one (AM-GM).

## [1.288.0] - 2026-09-11

### Added
- `seasoned_geometric_asian` (in `exotics.py`): prices an in-progress discrete
  geometric-average Asian where some fixings have already been observed. The
  observed prices contribute a known constant `A = sum log(S_obs)` and the `k`
  remaining log-prices stay jointly Gaussian, so `log G` is still
  `Normal(m, v)` with `m = (A + k log S + (b - sigma^2/2) sum tau_j)/n` and
  `v = (sigma^2/n^2) sum_ij min(tau_i, tau_j)` — an exact Black-Scholes-style
  closed form on the remaining window.
- Verified: with no observations it reduces exactly to
  `discrete_geometric_asian`; with all fixings observed the payoff is the
  deterministic discounted geometric intrinsic; a partial-window case matches a
  path Monte Carlo (n=12, 6 observed) to within 1%; locking in high observed
  fixings raises the call versus at-the-money observations.

## [1.287.0] - 2026-09-11

### Added
- `discrete_arithmetic_asian` and `discrete_arithmetic_asian_greeks` (in
  `exotics.py`): discretely-monitored arithmetic-average Asian by Levy (1992)
  two-moment lognormal matching. The average's first two moments over the
  fixing dates are exact — `M1 = (S/n) sum_i exp(b t_i)` and
  `M2 = (S^2/n^2) sum_ij exp(b(t_i+t_j) + sigma^2 min(t_i,t_j))` — and are
  matched to a lognormal priced with a Black-Scholes formula on `M1` with
  variance `log(M2/M1^2)`. Accepts `n_fixings` or explicit `fixing_times`.
  Greeks by finite difference.
- Verified: a single fixing recovers the vanilla Black-Scholes price; the
  arithmetic Asian sits above the exact discrete geometric Asian (AM-GM) and
  below the vanilla call; explicit equally-spaced times match the `n_fixings`
  grid; matches a geometric-control-variate Monte Carlo (n=12, call/put) to
  within 1%.

## [1.286.0] - 2026-09-11

### Added
- `discrete_geometric_asian` and `discrete_geometric_asian_greeks` (in
  `exotics.py`): exact closed form for a discretely-monitored geometric-average
  Asian option. `log G` over the fixing dates is Gaussian with mean
  `log S + (b - sigma^2/2) mean(t_i)` and variance
  `(sigma^2/n^2) sum_ij min(t_i, t_j)`, giving a Black-Scholes-style price on
  the lognormal average. Accepts equally-spaced `n_fixings` or an explicit
  `fixing_times` sequence. Greeks by finite difference.
- Verified: a single fixing recovers the vanilla Black-Scholes price; the price
  converges (monotonically, from above) to the continuous Kemna-Vorst
  `geometric_asian` as `n_fixings` grows; explicit equally-spaced fixing times
  match the `n_fixings` grid; a discrete geo Asian call is cheaper than the
  vanilla call; matches a path Monte Carlo (n=12, call/put) to within 1%.

## [1.285.0] - 2026-09-11

### Added
- `powered_option` and `powered_option_greeks` (in `exotics.py`): the powered
  option, payoff `max(S_T - K, 0)**p` (call) or `max(K - S_T, 0)**p` (put) for
  a positive integer power `p`. Distinct from `power_option` (payoff
  `max(S_T**p - K, 0)`): here the option payoff itself is raised to a power.
  Closed form by binomial expansion of the polynomial payoff into truncated
  moments `E[S_T**j] N(d_j)` (Esser 2003; Heynen-Kat 1996). Greeks by finite
  difference.
- Verified: `p = 1` recovers the vanilla Black-Scholes price and Greeks;
  matches an antithetic Monte Carlo for `p = 2` (call/put) and `p = 3` (call)
  to within statistical error; zero-vol limit equals the discounted powered
  intrinsic; delta/gamma/vega positive for the call.

## [1.284.0] - 2026-09-11

### Added
- `ultima` (in `greeks2.py`): the third-order vega Greek,
  `d(vomma)/d(sigma) = d^3(price)/d(sigma)^3`, useful for the convexity of a
  volga hedge. Analytic: `ultima = -(vega/sigma^2) * [d1 d2 (1 - d1 d2) + d1^2 + d2^2]`.
  Same for calls and puts.
- Verified: matches a central finite difference of `vomma` in sigma to
  6-7 digits; matches a central third difference of the BSM price in sigma;
  call equals put; ATM ultima is negative (volga concave in sigma there).

## [1.283.0] - 2026-09-11

### Added
- `double_knock_in_call_greeks` (in `exotics.py`): Greeks of a double-barrier
  knock-in call by in-out parity. Spot/vol/time Greeks are the vanilla
  Black-Scholes Greek minus the double knock-out Greek; barrier sensitivities
  are the negatives of the knock-out's (the vanilla has no barrier dependence).
- Verified: `knock_in + knock_out` Greeks sum to the vanilla Greek; `dV/dL` and
  `dV/dU` negate the knock-out's; a wider corridor lowers the knock-in
  (`dV/dL > 0`, `dV/dU < 0`); the price field matches the pricer.

## [1.282.0] - 2026-09-11

### Added
- `double_knock_out_call_greeks` (in `exotics.py`): finite-difference Greeks of
  the Ikeda-Kunitomo double knock-out call — `delta`, `gamma`, `vega`, `theta`,
  plus the two barrier sensitivities `dV/dL` and `dV/dU`.
- Verified: `delta` matches an independent bump; `vega < 0` (short vol);
  widening the corridor raises value (`dV/dL < 0`, `dV/dU > 0`), both matching
  direct bumps; the price field matches the pricer.

## [1.281.0] - 2026-09-11

### Added
- `surface_arbitrage_report` and `surface_is_arbitrage_free` (in `rnd.py`):
  combine the butterfly (per-slice density-sign) and calendar (total-variance
  monotonicity) static-arbitrage checks into one surface-level verdict. Takes
  smiles in log-moneyness (one per expiry) and returns
  `{"butterfly": {t: [strikes]}, "calendar": [(k, t_lo, t_hi)]}`.
- Verified: an arbitrage-free SSVI surface passes (agreeing with SSVI's own
  check); a butterfly-violating slice is reported under its expiry; a
  calendar-violating term structure is reported; a flat surface is free.

## [1.280.0] - 2026-09-11

### Added
- `calendar_arbitrage_violations` and `surface_is_calendar_arbitrage_free` (in
  `rnd.py`): model-free calendar-arbitrage detection across a smile term
  structure. Flags `(k, t_lo, t_hi)` points where total implied variance
  `w(k, t) = sigma^2 t` decreases with maturity — a horizontal-spread arbitrage.
  Works for any smiles expressed in log-moneyness (SVI, SABR, vanna-volga, raw).
- Verified: a flat-vol term structure is arbitrage-free; a dropping total
  variance is flagged; the verdict agrees with SSVI's own calendar check;
  length-mismatch and non-increasing expiries raise.

## [1.279.0] - 2026-09-11

### Added
- `smile_arbitrage_violations` and `smile_is_arbitrage_free` (in `rnd.py`):
  model-free butterfly-arbitrage detection for any implied-vol smile. Scans a
  log-moneyness grid and flags strikes where the Breeden-Litzenberger density is
  negative (a negative-cost butterfly). Works for SVI, SABR, vanna-volga, or raw
  quotes.
- Performance: the extended-Sobol geometric-Asian direction-number test now runs
  at `n_rand=8` instead of 24 across the 2..12 step sweep (`4*SE` still catches a
  broken direction number), cutting ~6s off the fast gate.
- Verified: flat and convex-SVI smiles are flagged arbitrage-free; a
  butterfly-violating SVI slice is flagged; the density-sign verdict agrees with
  SVI's analytic g-function across slices.

## [1.278.0] - 2026-09-11

### Added
- `double_knock_in_call` (in `exotics.py`): double-barrier knock-in call priced
  by in-out parity, `vanilla - double_knock_out_call`. A knock-in and knock-out
  with the same strike and corridor partition every path, so they sum to the
  vanilla call.
- Verified: `DKI + DKO = vanilla` exactly; wide barriers give a near-zero
  knock-in; a tighter band raises the knock-in value; the price stays in
  `[0, vanilla]`.

## [1.277.0] - 2026-09-11

### Added
- `double_knock_out_call` (in `exotics.py`): the Ikeda-Kunitomo (1992)
  double-barrier knock-out call — payoff `max(S_T - K, 0)` paid only if the
  continuously-monitored spot stays inside a corridor `(L e^{delta1 s},
  U e^{delta2 s})`. Prices via the truncated Ikeda-Kunitomo image series
  (`delta1 = delta2 = 0` gives flat barriers).
- Verified: wide barriers recover the vanilla call; a tighter band and higher
  vol both lower the value; the price stays in `[0, vanilla]`; and it matches a
  fine-step (6000-step) continuously-monitored Monte Carlo (which sits just above
  the closed form by the usual discrete-monitoring bias).

## [1.276.0] - 2026-09-11

### Added
- `double_no_touch_greeks` (in `exotics.py`): finite-difference Greeks of a
  double-no-touch — `delta`, `gamma`, `vega`, `theta`, plus the two barrier
  sensitivities `dV/dL` and `dV/dU`.
- Verified: `delta` matches an independent bump; `vega < 0` and `gamma < 0`
  (short vol, peaked inside the band); `theta > 0` (less time to knock);
  widening the band raises value (`dV/dL < 0`, `dV/dU > 0`), and both barrier
  sensitivities match direct bumps.

## [1.275.0] - 2026-09-11

### Added
- `double_no_touch` and `double_one_touch` (in `exotics.py`): continuously
  monitored double-barrier binaries paying `cash` if the spot stays inside
  `(L, U)` to expiry (DNT), or if it touches either barrier (DOT). DNT uses the
  Fourier eigenfunction expansion of the driftful survival probability in a
  strip.
- Verified: `U -> infinity` collapses DNT to the lower `no_touch`, `L -> 0` to
  the upper; `DNT + DOT = cash e^{-rt}`; a narrower band and higher vol both
  lower survival; and DNT matches a fine-step (4000-step) continuously-monitored
  Monte Carlo.

## [1.274.0] - 2026-09-10

### Added (tests)
- Cross-model surface density parity: an SSVI surface calibrated to SABR smiles
  at two expiries reproduces the SABR Breeden-Litzenberger risk-neutral density
  to ~4% at the money and ~6% in the wings at a fitted expiry; the fitted surface
  is arbitrage-free and both densities integrate to one. Complements the
  single-slice SVI/SABR parity check.

## [1.273.0] - 2026-09-10

### Added
- `volatility_swap_bounds_from_smile` (in `varswap.py`): brackets the fair
  volatility-swap strike between the at-the-money-forward implied vol (lower --
  the Carr-Lee zero-correlation proxy) and `sqrt(K_var)` (upper -- the Jensen
  bound from the variance-swap strike). The bracket width is the convexity /
  vol-of-vol premium the smile implies.
- Verified: a flat smile collapses the bracket; a convex smile keeps
  `lower < upper`; the upper bound equals `sqrt` of the variance-swap strike and
  the lower equals the ATMF vol; the premium widens with SVI curvature.

## [1.272.0] - 2026-09-10

### Changed
- Trimmed the fast gate: the RQMC knock-in/knock-out partition-identity tests
  (Parisian, barrier, barrier-digital) now run at `n_rand=6` instead of 24. The
  identity is exact path-by-path at a shared seed, so it holds for any `n_rand`;
  this drops ~10s off the non-slow suite with no loss of coverage.

### Added (tests)
- Cross-model density parity: an SVI slice calibrated to a SABR smile reproduces
  the SABR Breeden-Litzenberger risk-neutral density to ~2% at the money and ~5%
  in the wings, both densities are non-negative, and both integrate to one -- two
  independent constructions agreeing end to end.

## [1.271.0] - 2026-09-10

### Added
- `VannaVolgaSmile.vol_at_delta`, `.risk_reversal`, and `.butterfly` (in
  `vannavolga.py`): query the FX smile in delta space. `vol_at_delta` solves the
  fixed point `sigma = vol(strike_from_delta(sigma))`; `risk_reversal` and
  `butterfly` recompute the delta-consistent skew and convexity at any delta.
- Verified: at 25-delta the risk reversal and butterfly recover the inputs used
  to build the smile, and the pillar vols are returned exactly; a ~50-delta
  option sits near the ATM vol; a negative RR makes the 25d put richer than the
  25d call; a convex smile has a larger 10-delta butterfly than 25-delta.

## [1.270.0] - 2026-09-10

### Added
- `corrado_su_implied_vol` and `corrado_su_smile` (in `gramcharlier.py`): the
  Black-Scholes implied vol of a Corrado-Su price at a strike, and the smile over
  a strike grid. Turns Gram-Charlier `(sigma, skew, excess_kurt)` into an
  implied-vol curve for plotting or seeding an SVI/SABR fit.
- Verified: zero moments give a flat smile at `sigma`; negative skew lifts the
  low-strike put wing (monotone-decreasing smile); positive excess kurtosis lifts
  both wings above the at-the-money level; the implied vol reprices the
  Corrado-Su value; the smile round-trips through `calibrate_corrado_su`.

## [1.269.0] - 2026-09-10

### Added
- `calibrate_corrado_su` (in `gramcharlier.py`): fits the Corrado-Su
  `(sigma, skew, excess_kurt)` to a set of market call prices by Nelder-Mead
  least squares, with a smooth reparametrization keeping `sigma > 0`. Returns the
  three parameters plus the price RMSE.
- Verified: recovers known `(sigma, skew, kurt)` from synthetic Corrado-Su
  prices to 1e-2; a flat Black-Scholes surface calibrates to zero skew and
  excess kurtosis at the input vol; reprices every strike within RMSE.

## [1.268.0] - 2026-09-10

### Added
- `implied_spread_correlation_bs` (in `multiasset.py`): the correlation implied
  by a spread-option price under the Bjerksund-Stensland (2014) model, recovered
  by bisection on the price (monotone decreasing in `rho`). Companion to the Kirk
  `implied_spread_correlation`.
- Verified: round-trips a known `rho` for both calls and puts; the price is
  monotone decreasing in correlation; a quote above the `rho = -1` maximum
  raises.

## [1.267.0] - 2026-09-10

### Added
- `spread_option_bs_greeks` (in `multiasset.py`): finite-difference Greeks of the
  Bjerksund-Stensland (2014) spread option — two spot deltas, own-gammas,
  cross-gamma `d2V/dS1 dS2`, and `corr_vega` — matching the `spread_greeks`
  layout.
- Verified: delta1 > 0, delta2 < 0, positive gammas, negative cross-gamma and
  correlation sensitivity; delta1 matches an independent bump; call-minus-put
  deltas are `(+1, -1)`; the cross-gamma tracks `-sqrt(gamma1 gamma2)`; gammas
  are identical for call and put.

## [1.266.0] - 2026-09-10

### Added
- `spread_option_bs` (in `multiasset.py`): the Bjerksund-Stensland (2014)
  three-`d` closed-form spread-option approximation for `max(S1 - S2 - K, 0)`,
  generally more accurate than the Kirk approximation at wide strikes, high
  volatility, or dispersed leg vols, and exact (Margrabe) at `K = 0`.
- Verified: matches Margrabe at `K = 0`; satisfies put-call parity on the
  spread; agrees with 200k-400k-path Latin-hypercube Monte Carlo within a few
  standard errors; sits at least as close to MC as Kirk at a wide, dispersed-vol
  strike; the call decreases in correlation.

## [1.265.0] - 2026-09-10

### Added
- `risk_neutral_var_from_smile` and `risk_neutral_cvar_from_smile` (in `rnd.py`):
  risk-neutral Value-at-Risk and Conditional VaR (expected shortfall) of the
  terminal simple return `L = 1 - S_T/S0`, from an implied-vol smile.
  `VaR_alpha = 1 - Q(1-alpha)/S0` via the inverse smile CDF; CVaR integrates the
  truncated first moment `E[S_T ; S_T <= K]` against the Breeden-Litzenberger
  density.
- Verified against closed-form lognormal VaR/expected-shortfall on a flat smile;
  CVaR >= VaR; VaR rises with confidence; an equity-skew smile raises both tail
  metrics.

## [1.264.0] - 2026-09-10

### Added
- `risk_neutral_cdf_from_smile` and `risk_neutral_quantile_from_smile` (in
  `rnd.py`): the Breeden-Litzenberger risk-neutral CDF `F(K) = 1 + e^{rt} dC/dK`
  of the terminal spot from an implied-vol smile, and its inverse (bisection over
  a forward-standard-deviation bracket).
- Verified: a flat smile recovers the Black-Scholes `N(-d2)`; the CDF is monotone
  in `[0, 1]`; the quantile inverts the CDF; the CDF equals the integral of the
  Breeden-Litzenberger density; the median matches the lognormal
  `F e^{-sigma^2 t / 2}`; an equity-skew smile fattens the left tail.

## [1.263.0] - 2026-09-10

### Added
- SSVI surface smile consumers at a fitted expiry (in `ssvi.py`):
  `ssvi_variance_swap_strike`, `ssvi_vix`, `ssvi_svix`, `ssvi_density`
  (Breeden-Litzenberger), and `ssvi_bkm_moments`. Each maps a strike to the
  surface's Black vol `implied_vol(ln(K/F), t)` on the forward and feeds the
  existing smile machinery.
- Verified: on a flat-ATM (`theta_t = sigma^2 t`) surface the var-swap and
  VIX/SVIX recover `sigma`; the density integrates to 1 and has mean equal to
  the forward and is non-negative on an arbitrage-free slice; BKM skewness
  tracks the sign of the surface `rho`. Requesting an unfitted expiry raises.

## [1.262.0] - 2026-09-10

### Added
- `svi_density` (in `svi.py`): Breeden-Litzenberger risk-neutral density `g(K)`
  implied by a raw-SVI slice, `g(K) = e^{rt} d^2C/dK^2` with the call priced at
  the slice's smile vol on the forward.
- Verified: a flat slice matches the closed-form lognormal density; the density
  is non-negative on a butterfly-arbitrage-free slice, integrates to 1, and has
  mean equal to the forward; a wing-ok but butterfly-violating slice produces a
  negative density exactly where `svi_g < 0`.

## [1.261.0] - 2026-09-10

### Added
- `sabr_bkm_moments` (in `sabr.py`): risk-neutral variance, skewness, and excess
  kurtosis implied by a SABR smile via Bakshi-Kapadia-Madan moment replication.
- Verified: a near-flat slice is symmetric; negative `rho` gives negative
  risk-neutral skew (positive `rho` positive); higher vol-of-vol raises the
  excess kurtosis.

## [1.260.0] - 2026-09-10

### Added
- `svi_bkm_moments` (in `svi.py`): risk-neutral variance, skewness, and excess
  kurtosis implied by a raw-SVI slice via the Bakshi-Kapadia-Madan moment
  replication of the slice's smile.
- Verified: a flat slice is near-symmetric (skew ~ 0, excess kurtosis ~ 0);
  negative SVI `rho` gives negative risk-neutral skew and positive `rho`
  positive skew; a convex slice has positive excess kurtosis.

## [1.259.0] - 2026-09-10

### Added
- `sabr_variance_swap_strike` and `sabr_vix` (in `sabr.py`): the fair
  variance-swap strike (annualized variance) and the VIX-style index implied by
  a SABR smile, each pricing the strike chain at `sabr_vol(F, K, ...)` and
  replicating via `variance_swap_from_smile` / `vix_from_smile`.
- Verified: at `beta = 1`, `nu -> 0` the strike is `alpha^2` and the VIX is
  `100 * alpha`; a skewed slice's variance exceeds the ATM variance; higher
  vol-of-vol raises the strike.

## [1.258.0] - 2026-09-10

### Added
- `svi_vix` and `svi_svix` (in `svi.py`): the VIX-style and Martin SVIX indices
  implied by a raw-SVI slice, mapping each strike to `p.implied_vol(ln(K/F), t)`
  and replicating via `vix_from_smile` / `svix_from_smile` (both reported as
  `100 * index`).
- Verified: a flat slice returns `100 * sigma` for both; a skewed slice lifts the
  VIX above `100 *` the ATM vol; higher convexity (`b`) raises the index.

## [1.257.0] - 2026-09-10

### Added
- `svi_variance_swap_strike` (in `svi.py`): fair variance-swap strike
  (annualized variance) replicated model-consistently from a raw-SVI slice --
  maps each strike to `p.implied_vol(ln(K/F), t)` and feeds the smile to
  `variance_swap_from_smile`.
- Verified: a flat slice (`b = 0`) returns the flat variance `sigma^2`; a
  skewed/convex slice returns a variance above the ATM variance (the convexity
  premium); higher `b` (more convexity) raises the strike.

## [1.256.0] - 2026-09-10

### Added
- `DiscountCurve` analytics methods: `instantaneous_forward` (`-d ln DF/dT` by
  central FD), `annuity` (fixed-leg PV01), `swap_value` (unit-notional payer/
  receiver), and `swap_dv01` (value change for a 1bp parallel curve drop via a
  shifted-curve wrapper).
- Verified: on a flat curve the instantaneous forward equals the zero rate; a
  par-struck swap has zero value; the annuity matches a manual sum; a payer's
  DV01 is negative on a rate drop and a receiver's positive, of order
  annuity x 1bp.

## [1.255.0] - 2026-09-10

### Added
- Ho-Lee caps and floors: `holee_caplet`, `holee_floorlet`, `holee_cap`,
  `holee_floor`, each a strip via the bond-put/call identity. Completes cap/floor
  coverage across every short-rate model (Vasicek, CIR, Ho-Lee, Cheyette) plus
  G2++.
- Verified: `cap - floor` equals the underlying swap value (parity) to 1e-8; the
  cap equals the sum of its caplets; positive and monotone in strike.

## [1.254.0] - 2026-09-10

### Added
- Vasicek and CIR caps and floors: `vasicek_caplet`/`vasicek_floorlet`/
  `vasicek_cap`/`vasicek_floor` and `cir_caplet`/`cir_floorlet`/`cir_cap`/
  `cir_floor`, each a strip of caplets/floorlets via the bond-put/call identity
  on the model's exact bond option.
- Verified: for both models `cap - floor` equals the underlying swap value
  (parity) to 1e-8, and the cap equals the sum of its caplets; positive and
  monotone in strike.

## [1.253.0] - 2026-09-10

### Added
- `cheyette_floorlet`, `cheyette_cap`, `cheyette_floor` (in `cheyette.py`):
  Cheyette (Hull-White) floorlet (bond-call identity) and cap/floor as strips of
  caplets/floorlets.
- Verified: `cap - floor` equals the underlying swap value (parity) to 1e-8; cap
  and floor equal the sums of their legs; both positive and monotone in strike.

## [1.252.0] - 2026-09-10

### Added
- `g2pp_floorlet`, `g2pp_cap`, `g2pp_floor` (in `g2pp.py`): G2++ floorlet (via
  the bond-call identity) and cap/floor as strips of caplets/floorlets over a
  schedule of `(t_i, P(0, t_i))` dates.
- Verified: `cap - floor` equals the underlying swap value (put-call parity) to
  1e-8; the cap and floor equal the sums of their caplets/floorlets; both
  positive; a higher strike lowers the cap and raises the floor.

## [1.251.0] - 2026-09-10

### Added
- `sabr_option_greeks` now also returns the total `gamma` (`d2Price/dF2`),
  computed by a central difference of the SABR-repriced surface so it captures
  the smile backbone's curvature -- not just the vol-fixed Black gamma.
- Verified: the total gamma matches a finite difference that recomputes the SABR
  vol at each bumped forward, and is positive.

## [1.250.0] - 2026-09-10

### Added
- `sabr_option_greeks` (in `sabr.py`): Greeks of a Black-76 option priced at the
  Hagan SABR smile vol. The total `delta = black_delta + vega * dsigma/dF`
  includes the smile backbone (the exact AD `dsigma/dF` from
  `sabr_sensitivities`), so it differs from the vol-fixed Black delta. Returns
  `price`, `vol`, `delta` (backbone-adjusted), `black_delta`, and `vega`, with an
  optional `discount` to scale the forward figures to present value.
- Verified: the total delta matches a finite difference that recomputes the SABR
  vol at each bumped forward; it differs measurably from the Black delta; vega
  positive; `vol` equals `sabr_vol`; price and delta scale linearly in `discount`.

## [1.249.0] - 2026-09-10

### Added
- `cir_coupon_bond_option` and `cir_swaption` (in `cir.py`): exact European
  coupon-bond option and swaption under CIR by Jamshidian decomposition. The CIR
  bond is monotone in `r0`, so the critical-rate `r*` split into per-cashflow CIR
  zero-coupon-bond options (`cir_bond_option`) is exact; a payer swaption is a put
  on the fixed-leg coupon bond struck at the notional. Completes the full option
  chain (bond -> coupon-bond option -> swaption) for all three affine short-rate
  models: Vasicek, Ho-Lee, CIR.
- Verified: a single cashflow reduces to the scaled `cir_bond_option`; swaption
  parity `payer - receiver = annuity (swap_rate - strike)` holds to 1e-7; both
  legs positive; a higher strike lowers the payer.

## [1.248.0] - 2026-09-10

### Added
- `cir_bond_option` (in `cir.py`): exact European option on a CIR zero-coupon
  bond (Cox-Ingersoll-Ross 1985), using the noncentral chi-square CDF already in
  the library. The call is `P(0,t_bond) X2(...) - strike P(0,t_option) X2(...)`
  with critical rate `r* = ln(A/strike)/B`; puts follow from put-call parity.
- Verified: put-call parity `C - P = P(0,t_bond) - strike P(0,t_option)` to
  1e-9; the call matches a fine-step CIR Monte Carlo (0.02260, slow test); both
  prices positive; a higher strike lowers the call.

## [1.247.0] - 2026-09-10

### Added
- `holee_bond_option`, `holee_coupon_bond_option`, `holee_swaption` (in
  `holee.py`): exact Ho-Lee zero-coupon-bond option (Black-style with bond vol
  `sigma (t_bond - t_option) sqrt(t_option)`), coupon-bond option by Jamshidian
  decomposition, and European swaption via the coupon-bond-option identity (payer
  = put on the fixed-leg coupon bond struck at notional).
- Verified: a single cashflow reduces to the scaled zero-coupon bond option;
  the bond option satisfies put-call parity; swaption parity
  `payer - receiver = annuity (swap_rate - strike)` holds to 1e-8; both swaption
  legs positive.

## [1.246.0] - 2026-09-10

### Added
- `vasicek_swaption` (in `vasicek.py`): exact European swaption under Vasicek via
  the coupon-bond-option identity. The fixed leg plus notional is a coupon bond;
  a payer swaption is a put on it struck at the notional and a receiver a call,
  both priced by the Jamshidian `coupon_bond_option` -- no normal-model
  approximation.
- Verified: swaption parity `payer - receiver = annuity (swap_rate - strike)`
  holds to 1e-8; at the forward swap rate payer equals receiver; a higher strike
  lowers the payer; both legs positive.

## [1.245.0] - 2026-09-10

### Added
- `vasicek_coupon_bond_option` (in `vasicek.py`): European option on a
  coupon bond under Vasicek by Jamshidian's decomposition. Solves for the
  critical short rate `r*` where the bond value at expiry equals the strike, then
  prices the option as the `c_i`-weighted sum of zero-coupon-bond options struck
  at `K_i = P(t_option, t_i | r*)` -- exact, no simulation.
- Verified: a single cashflow reduces to the scaled zero-coupon `bond_option`;
  put-call parity `C - P = sum c_i P(0,t_i) - K P(0,t_option)` holds to 1e-8;
  prices positive; a higher strike lowers the call.

## [1.244.0] - 2026-09-10

### Added
- `cms_adjustment_greeks` (in `cms.py`): sensitivities of the standard-model CMS
  convexity adjustment -- `d_forward` and `d_sigma` by central finite differences
  on `cms_adjustment_standard`.
- Verified: both match finite differences; `d_sigma > 0` (more vol means more
  convexity, a larger adjustment); a zero-vol adjustment has zero sensitivity.

## [1.243.0] - 2026-09-10

### Added
- `swaption_greeks` (in `rates.py`): analytic Greeks of a normal-model European
  swaption. The value is `annuity * Bachelier(swap_rate, strike, expiry, 0,
  sigma_n)`, so `rate_delta`, `rate_gamma`, and `vega` are the Bachelier Greeks
  scaled by the annuity (also returned).
- Verified: `rate_delta` matches a finite difference of `swaption_price`; a payer
  swaption's rate delta is positive and a receiver's negative; the ATM payer
  delta is exactly `annuity/2` (Bachelier ATM delta 0.5); vega positive.

## [1.242.0] - 2026-09-10

### Added
- `caplet_greeks`, `cap_greeks`, `floor_greeks` (in `rates.py`): analytic Greeks
  of normal-model caplets/floorlets and their strips. Each caplet's value is
  `discount * accrual * Bachelier(F, K, ...)`, so its `rate_delta`, `rate_gamma`,
  and `vega` are the Bachelier Greeks in the forward rate scaled by the same
  factor; `cap_greeks`/`floor_greeks` sum them across periods.
- Verified: the cap `rate_delta` matches a finite difference of `cap_price`; a
  cap's rate delta is positive and a floor's negative; both vegas positive; the
  cap Greeks equal the sum of the caplet Greeks.

## [1.241.0] - 2026-09-10

### Added
- `dual_swap_dv01` (in `dualcurve.py`): risk of a dual-curve (OIS-discounted,
  projection-forward) swap -- the exact fixed-leg `pv01` (annuity), the total
  `dv01` for a 1bp parallel drop of both curves, and the split `ois_dv01` /
  `proj_dv01` from shifting only the discount or only the projection curve.
  Curve shifts use a light `_ShiftedCurve` wrapper (`df(T) e^{-dr T}`).
- Verified: `pv01` equals `-dV/d(fixed_rate)` for a payer to machine precision;
  the total `dv01` is the sum of the per-curve DV01s to first order; a payer's
  `dv01` is negative on a rate drop and a receiver's positive.

## [1.240.0] - 2026-09-10

### Added
- `bermudan_swaption_g2pp_greeks` (in `bermudan_swaption.py`): Greeks of a G2++
  Bermudan swaption by common-random-number bumps on `bermudan_swaption_g2pp` --
  `d_fixed` (dV/d fixed_rate) and `curve_dv01` (value change for a 1bp parallel
  curve drop). Same-seed repricing shares the G2++ state paths for low-variance
  differences.
- Verified: `d_fixed` matches a CRN finite difference; it is negative for a payer
  (paying a higher fixed rate is worth less) and positive for a receiver;
  reproducible under a fixed seed.

## [1.239.0] - 2026-09-10

### Added
- `vasicek_bond_greeks` (in `vasicek.py`) and `holee_bond_greeks` (in
  `holee.py`): exact rate sensitivities of the Vasicek and Ho-Lee zero-coupon
  bonds (`rho_r`, `gamma_r`, `duration`, `convexity`). Vasicek uses the affine
  `B(t)`; Ho-Lee is linear in `r0` inside the exponent, so its duration is
  exactly the maturity `t` and convexity `t^2`.
- Verified: both `rho_r` match finite differences of their bond prices;
  `convexity = duration^2`; the Ho-Lee duration equals the maturity; zero-
  maturity bonds are flat.

## [1.238.0] - 2026-09-10

### Added
- `cir_bond_greeks` (in `cir.py`): exact rate sensitivities of a CIR
  zero-coupon bond. Since `P = A(t) e^{-B(t) r0}`, the short-rate delta is
  `rho_r = -B P`, gamma `gamma_r = B^2 P`, rate duration `B`, and convexity
  `B^2` -- all closed form.
- Verified: `rho_r` and `gamma_r` match finite differences of
  `cir_zero_coupon_bond`; `rho_r < 0`, `gamma_r > 0`, `convexity = duration^2`,
  `duration = -rho_r / price`; a zero-maturity bond is flat.

## [1.237.0] - 2026-09-10

### Added
- `cheyette_bond_option_greeks` (in `cheyette.py`): exact discount-factor deltas
  `delta_T` (`dV/dP0T = N(d1)`, underlying-bond delta) and `delta_S` (`dV/dP0S`),
  plus the vol sensitivity `vega`, of a Cheyette (Hull-White) zero-coupon-bond
  option.
- Verified: both deltas match finite differences of `cheyette_bond_option`; the
  call's `delta_T` is in (0,1) and `delta_S` negative; vega positive; the put's
  `delta_T` negative.

## [1.236.0] - 2026-09-10

### Added
- `vasicek_bond_option_greeks` (in `vasicek.py`): short-rate sensitivities
  `rho_r` (dV/dr0) and `gamma_r` (d2V/dr0^2), plus the vol sensitivity `vega`
  (dV/dsigma) of a Vasicek zero-coupon-bond option by central finite differences
  on `bond_option`.
- Verified: `rho_r` matches a finite difference; a bond call's `rho_r` is
  negative (it loses value as the short rate rises) and a put's is positive;
  vega is positive.

## [1.235.0] - 2026-09-10

### Added
- `g2pp_bond_option_greeks` (in `g2pp.py`): Greeks of a G2++ zero-coupon-bond
  option -- the exact discount-factor deltas `delta_T` (`dV/dP0T = N(d1)`, the
  underlying-bond delta) and `delta_S` (`dV/dP0S`), plus the two factor-vol
  vegas `vega_sigma` and `vega_eta` by finite difference.
- Verified: both deltas match finite differences of `g2pp_bond_option` (delta_T
  equals the analytic `N(d1)`); the call's `delta_T` is in (0,1) and `delta_S`
  negative; both vegas positive; the put's `delta_T` is negative.

## [1.234.1] - 2026-09-10

### Tests
- Trimmed the fast suite ~31s -> ~20s. The two-asset LSM-Greeks and
  rBergomi-Greeks tests each run 5+ re-prices; the `test_reproducible`
  determinism checks (which price twice) dropped to tiny sizes (n_steps 6-10,
  n_paths 800-1500), and the sign/ordering checks (min-put `delta2 < delta1`,
  basket `delta1 > delta2`, both-deltas-negative) to 10 steps / 3000-4000 paths
  -- all robust across seeds at the smaller size. Accuracy cross-checks stay
  under `-m slow`.

## [1.234.0] - 2026-09-10

### Added
- `andreasen_huge_strike_greeks` (in `andreasenhuge.py`): strike-space Greeks of
  the arbitrage-free Andreasen-Huge call surface -- the `dual_delta` (`dC/dK`,
  discounted; equals `-e^{-rT}` times the risk-neutral exceedance probability)
  and the `rnd` (Breeden-Litzenberger risk-neutral density `e^{rT} d2C/dK2`) at
  each interior strike.
- Verified: the density is non-negative for any positive local vols (the
  scheme's convexity), the dual delta is monotone in `[-e^{-rT}, 0]`, and the
  density integrates to ~1 -- rising toward 1 as the strike grid widens
  (0.89 -> 0.994 -> 0.9997), confirming the tail-truncation deficit rather than a
  normalization error.

## [1.233.0] - 2026-09-10

### Added
- `rbergomi_greeks_cv` (in `rbergomi.py`): delta, gamma, and `vega_xi0`
  (forward-variance-level sensitivity) of a rough-Bergomi call by
  common-random-number bumps on the conditional control-variate estimator
  `rbergomi_price_cv`. Same-seed repricing shares the volatility-driving Brownian
  paths, so the finite differences are low-variance.
- Verified: delta matches a common-random-number finite difference of
  `rbergomi_price_cv` to machine precision; call delta in (0,1) with positive
  gamma and `vega_xi0`; reproducible under a fixed seed.

## [1.232.0] - 2026-09-10

### Added
- `rough_heston_greeks` (in `rough_heston.py`): spot Greeks (`delta`, `gamma`)
  and the initial-variance sensitivity `vega_v0` of a rough-Heston option by
  central finite differences on `rough_heston_price`. Each re-price runs the
  O(n_grid^2) fractional-Riccati solve, so it is comparatively slow; the default
  `n_grid` matches the pricer's (small `H` needs a fine grid to stay stable).
- Verified: at `H = 0.5` the delta matches a finite difference of the classical
  Heston price (`xi = kappa * nu`); at `H = 0.3` (rough) the delta matches a
  finite difference of `rough_heston_price`, with call delta in (0,1), positive
  gamma and `vega_v0`.

## [1.231.0] - 2026-09-10

### Added
- `double_heston_greeks` (in `double_heston.py`): spot Greeks (`delta`,
  `gamma`) plus a per-factor initial-variance sensitivity `vega_v01` and
  `vega_v02` of a double-Heston option by central finite differences on
  `double_heston_price` (one-sided at each variance floor).
- Verified: delta matches a finite difference; call delta in (0,1) with positive
  gamma; both variance vegas are positive.

## [1.230.0] - 2026-09-10

### Added
- `meixner_greeks` (in `meixner.py`): spot Greeks (`delta`, `gamma`, `theta`)
  plus the asymmetry/skew sensitivity `d_b` (dV/db) of a Meixner option by
  central finite differences on `meixner_price` (the `d_b` bump is clipped to
  keep `b` in `(-pi, pi)`). Completes the Levy-model Greek family (Merton, Kou,
  VG, NIG, CGMY, Meixner).
- Verified: delta matches a finite difference; call delta in (0,1) with positive
  gamma, put delta negative; the skew sensitivity is finite.

## [1.229.0] - 2026-09-10

### Added
- `cgmy_greeks` (in `cgmy.py`): spot Greeks (`delta`, `gamma`, `theta`) plus the
  tail-activity sensitivity `d_Y` (dV/dY) of a CGMY option by central finite
  differences on `cgmy_price` (`d_Y` falls back to a one-sided bump near `Y = 2`).
- Verified: delta matches a finite difference; call delta in (0,1) with positive
  gamma, put delta negative; the tail sensitivity is finite.

## [1.228.0] - 2026-09-10

### Added
- `bates_greeks` (in `bates.py`): spot Greeks (`delta`, `gamma`), the
  initial-variance sensitivity `vega_v0`, and the jump-intensity sensitivity
  `d_lambda` of a Bates (Heston + Merton jumps) option by central finite
  differences on `bates_price` (one-sided at the `v0`/`lam` floors).
- Verified: at `lam = 0` the delta matches a finite difference of the exact
  Heston price; with jumps the delta matches a finite difference of `bates_price`;
  call delta in (0,1) with positive gamma and `vega_v0`.

## [1.227.0] - 2026-09-10

### Added
- `nig_greeks` (in `nig.py`): spot Greeks (`delta`, `gamma`, `theta`) plus
  process-parameter sensitivities (`d_alpha` tail steepness, `d_beta` skew) of a
  NIG option by central finite differences on `nig_price`. The `d_beta` bump is
  clipped to keep `|beta| < alpha` on both sides.
- Verified: delta matches a finite difference; call delta in (0,1) with positive
  gamma, put delta negative; the parameter sensitivities are finite.

## [1.226.0] - 2026-09-10

### Added
- `variance_gamma_greeks` (in `variancegamma.py`): delta, gamma, vega
  (Brownian-vol sensitivity), and `theta_greek` (calendar decay) of a
  Variance-Gamma option by central finite differences on `variance_gamma_price`.
  The calendar Greek is named `theta_greek` to avoid clashing with the VG skew
  *parameter* `theta`.
- Verified: a small `nu` gives a delta close to the Black-Scholes delta (the
  `nu -> 0` limit); the delta matches a finite difference; call delta in (0,1)
  with positive gamma/vega, put delta negative.

## [1.225.0] - 2026-09-10

### Added
- `kou_greeks` (in `kou.py`): delta, gamma, vega (diffusion-vol sensitivity),
  and theta of a Kou double-exponential jump-diffusion option by central finite
  differences on `kou_price`.
- Verified: at `lam = 0` (no jumps) the Greeks reduce exactly to the vanilla
  Black-Scholes Greeks; with jumps the delta matches a finite difference, the
  call delta is in (0,1) with positive gamma/vega, and the put delta is negative.

## [1.224.0] - 2026-09-10

### Added
- `bachelier_greeks` (in `bachelier.py`): bundles the Bachelier (normal-model)
  Greeks -- `delta`, `gamma`, `vega` from the existing exact closed forms plus
  `theta` (calendar decay) by central finite difference.
- Verified: delta/gamma/vega equal their standalone closed forms; theta matches
  a finite difference; the ATM call delta is exactly 0.5 at zero rate; the put
  delta is negative.

## [1.223.0] - 2026-09-10

### Added
- `merton_jump_greeks` (in `merton.py`): delta, gamma, vega (diffusion-vol
  sensitivity), and theta of a Merton jump-diffusion option by central finite
  differences on `merton_jump_price`.
- Verified: at `lam = 0` (no jumps) the Greeks reduce exactly to the vanilla
  Black-Scholes Greeks; with jumps the delta matches a finite difference, the
  call delta is in (0,1) with positive gamma/vega, and the put delta is negative.

## [1.222.0] - 2026-09-10

### Added
- `cev_greeks` (in `cev.py`): delta, gamma, vega, and theta of a CEV option by
  central finite differences on `cev_price`.
- Verified: delta matches a finite difference across `beta` in {0.3, 0.5, 0.7};
  the call delta is in (0,1) with positive gamma/vega and the put delta is
  negative; a lower `beta` (steeper local-vol skew) lifts the ATM call delta.
  (The `beta -> 1` Black-Scholes limit is *not* asserted -- the Schroder
  noncentral-chi-square parameters are numerically unstable as `beta -> 1`, so
  the checks stay in the well-conditioned interior.)

## [1.221.0] - 2026-09-10

### Added
- `cliquet_greeks` (in `forwardstart.py`): delta, gamma, vega, and theta of a
  cliquet (ratchet) by central finite differences on `cliquet_price` (theta
  shifts every reset date together).
- Verified: every strike scales with the spot (`alpha * S_reset`), so the
  cliquet is homogeneous of degree 1 in `S` -- `delta` equals `price / S`
  exactly and `gamma` is zero. Delta and vega also match finite differences of
  `cliquet_price`.

## [1.220.0] - 2026-09-10

### Added
- `forward_start_greeks` (in `forwardstart.py`): Greeks of a forward-start
  option. Since the Rubinstein price `S e^{(b-r) t_start} u` is exactly linear in
  the spot (the unit price `u` does not depend on `S`), `delta = e^{(b-r) t_start} u`
  is constant in the spot and `gamma = 0` -- a forward-start has no spot gamma
  until its strike is fixed. `vega` and `theta` are central finite differences of
  the closed form.
- Verified: delta matches a finite difference and is constant across spot levels;
  gamma is exactly zero; vega matches a finite difference; theta is zero at
  `b = r` (only the maturity gap matters) and nonzero under a carry.

## [1.219.0] - 2026-09-10

### Added
- `displaced_diffusion_greeks` (in `displaced.py`): delta, gamma, vega, and
  theta of a displaced-diffusion option by central finite differences on
  `displaced_diffusion_price`.
- Verified: at `shift = 0` the Greeks reduce exactly to the vanilla
  Black-Scholes Greeks; the delta matches a finite difference at a positive
  shift; call delta in (0,1) with positive gamma/vega, put delta negative.

## [1.218.0] - 2026-09-10

### Added
- `gap_option_greeks` (in `exotics.py`): delta, gamma, vega, and theta of a gap
  option (Reiner-Rubinstein) by central finite differences on `gap_option`.
- Verified: setting `K_trigger = K_payoff` reduces the Greeks exactly to the
  vanilla Black-Scholes Greeks; the gap delta matches a finite difference for
  both call (trigger above payoff) and put.

## [1.217.0] - 2026-09-10

### Added
- `compo_option_greeks` (in `quanto.py`): Greeks of a composite (compo) FX
  option. The price is a Black-Scholes price on the domestic-currency asset with
  the combined vol `sqrt(sigma_asset^2 + sigma_fx^2 + 2 rho sigma_asset sigma_fx)`
  and carry `r_domestic - q`, so `delta`/`gamma` are exact BSM Greeks; `vega`,
  `fx_vega`, and `corr_vega` are finite differences.
- Verified: delta matches a finite difference; at `sigma_fx = 0` the delta
  equals the plain BSM delta at carry `r_domestic - q`; a compo is long FX
  volatility (`fx_vega > 0`), unlike a quanto; call gamma and vega positive.

### Fixed
- The `fx_vega` finite difference in `quanto_option_greeks` /
  `compo_option_greeks` probed a negative `sigma_fx` when `sigma_fx` sat below
  the bump size (e.g. `sigma_fx = 0`), raising `ValueError`. It now falls back to
  a one-sided difference at the volatility floor.

## [1.216.0] - 2026-09-10

### Added
- `quanto_option_greeks` (in `quanto.py`): Greeks of a quanto option. The price
  is a Black-Scholes price on the foreign asset with the quanto-adjusted carry
  `b_q = r_foreign - q_asset - rho sigma_asset sigma_fx`, so the spot enters only
  through that BSM price -- `delta` and `gamma` are the exact BSM Greeks at `b_q`
  (no finite difference). `vega` (asset vol), `fx_vega` (FX vol), and `corr_vega`
  (correlation) are central finite differences of the closed form.
- Verified: delta matches a finite difference; at `rho = 0` the delta equals the
  plain BSM delta at carry `r_foreign - q`; `corr_vega` is negative for a
  positive-rho call (higher correlation lowers the quanto carry); call gamma and
  asset vega are positive.

## [1.215.0] - 2026-09-10

### Added
- `power_option_greeks` (in `exotics.py`): delta, gamma, vega, and theta of a
  power option (`max(S_T^power - K, 0)`) by central finite differences on the
  closed-form `power_option`.
- Verified: at `power = 1` the Greeks reduce exactly to the vanilla
  Black-Scholes Greeks; at `power = 2` delta matches a finite difference and the
  call gamma/vega are positive; the put delta is negative.

## [1.214.0] - 2026-09-10

### Added
- `chooser_option_greeks` (in `chooser.py`): delta, gamma, and vega of a simple
  chooser option, exact by decomposition. The chooser is exactly a call to `T`
  plus a put struck at the discounted-forward level expiring at the choice date,
  and both legs are Black-Scholes prices in `S` and `sigma`, so the Greeks are
  the exact sums of the two legs' BSM Greeks -- no finite difference.
- Verified: delta, gamma, and vega match finite differences of `chooser_option`
  to step precision; gamma and vega are positive (long both a call and a put).

## [1.213.0] - 2026-09-10

### Added
- `compound_option_greeks` (in `compound.py`): delta, gamma, vega, and theta of
  a Geske compound option by central finite differences on `compound_option`
  (all four kinds: call/put on call/put). Theta shifts both expiries together.
- Verified: delta matches a finite difference of the closed form; a
  call-on-call has positive delta, gamma, and vega; a put-on-call has negative
  delta. Noted (via cross-check) that vega is *not* universally positive -- a put
  on an option is short the compound optionality, so `put-on-call` and
  `put-on-put` can have negative vega; only the call compounds are asserted
  positive.

## [1.212.0] - 2026-09-10

### Added
- `geometric_asian_greeks` (in `exotics.py`): Greeks of the continuously-
  monitored geometric-average Asian. Since the Kemna-Vorst price is exactly a
  Black-Scholes price at the adjusted vol `sigma_A = sigma/sqrt(3)` and carry
  `b_A = (b - sigma^2/6)/2`, the spot enters only through that BSM price, so
  `delta` and `gamma` are the exact BSM Greeks there (no finite difference);
  `vega`, `theta`, and `rho` are central differences of the exact closed form.
- Verified: delta, gamma, and vega match finite differences of `geometric_asian`
  to machine / step precision; call delta in (0,1), gamma and vega positive; put
  delta negative; the `price` field equals `geometric_asian`.

## [1.211.0] - 2026-09-10

### Added
- `bermudan_basket_lsm_greeks` (in `lsm.py`): deltas, own-gammas, and
  cross-gamma of an American basket option by common-random-number bumps on
  `bermudan_basket_lsm`. Same-seed repricing shares the Brownian shocks so the
  finite differences are low-variance; the regression is re-fit at each bump.
  Completes the two-asset LSM Greek family (max-call, spread, min-put, basket).
- Verified: without dividends the two spot deltas match the European
  `basket_greeks` deltas within noise; the call deltas are positive and the
  larger-weight asset carries the larger delta; the put deltas are negative;
  reproducible.

## [1.210.0] - 2026-09-10

### Added
- `sobol_arithmetic_asian_rqmc` (in `sobol.py`): randomized-QMC arithmetic Asian
  with a geometric control variate -- the two strongest variance-reduction
  techniques together. With `control_variate=True` each path's estimator is
  `arith - geo + E[geo]`, where `E[geo]` is the exact discrete-geometric closed
  form over the same `n_steps` dates and the two averages (nearly perfectly
  correlated) share the path. Bridge construction + per-dimension
  Cranley-Patterson rotation.
- Verified: matches the control-variate `arithmetic_asian_mc` for call and put;
  the control cuts the standard error to ~0.11x the RQMC-only estimator (on top
  of the QMC gain), so ~9x below RQMC alone.

## [1.209.0] - 2026-09-10

### Changed
- Extended the Sobol generator from 6 to 12 dimensions, adding the canonical
  Joe-Kuo primitive-polynomial coefficients and direction-number seeds for dims
  7-12. Validated end-to-end: the RQMC geometric-Asian price matches the exact
  discrete-geometric closed form at every `n_steps` from 2 to 12 (a wrong
  polynomial/seed would break the low-discrepancy property and bias it),
  regression-tested by `test_extended_sobol_dims_match_closed_form`. This lets
  every `sobol_*_rqmc` routine use up to 12 monitoring dates.

### Added
- `sobol_parisian_rqmc` (in `sobol.py`): randomized-QMC Parisian barrier option
  with an honest standard error -- knock-out/knock-in triggered only after the
  spot spends a *consecutive* `window` on the barrier's far side. Bridge
  construction + per-dimension Cranley-Patterson rotation; the discrete analogue
  of `parisian_barrier_mc`.
- Verified: matches `parisian_barrier_mc` at equal `n_steps` for down-out and
  down-in; a knock-in plus its knock-out sum to the vanilla at the same seed; a
  longer required window raises the knock-out value.

## [1.208.0] - 2026-09-10

### Added
- `sobol_geometric_asian_rqmc` (in `sobol.py`): randomized-QMC geometric-average
  Asian option (call `max(G - K, 0)`, put `max(K - G, 0)`) with an honest
  standard error. Bridge construction + per-dimension Cranley-Patterson
  rotation. Because the discrete geometric average is exactly lognormal it has a
  closed form (`_discrete_geometric_asian`), so this is cross-checked against a
  *deterministic* reference -- the tightest possible.
- Verified: matches the exact closed form for call and put (also under a
  dividend carry) within the tiny RQMC standard error (< 0.02).

## [1.207.0] - 2026-09-10

### Added
- `bermudan_basket_lsm` (in `lsm.py`): American basket option on
  `w1 S1 + w2 S2` (call or put) by Longstaff-Schwartz. Two correlated GBMs; the
  continuation value is regressed on a quadratic basis in both spots plus the
  basket `{1, S1, S2, S1^2, S2^2, S1 S2, w1 S1 + w2 S2}` over the in-the-money
  paths.
- Verified: without dividends the call matches the European moment-matched
  `basket_option`; a 6% dividend on both assets gives a positive early-exercise
  premium on the call (~0.38); the put carries an early-exercise premium even
  without dividends (~0.26); reproducible.

## [1.206.0] - 2026-09-10

### Added
- `sobol_average_strike_rqmc` (in `sobol.py`): randomized-QMC average-strike
  Asian option with an honest standard error. The strike is the realized
  arithmetic average of the path (call `max(S_T - A, 0)`, put `max(A - S_T, 0)`).
  Bridge construction + per-dimension Cranley-Patterson rotation; the discrete
  analogue of `average_strike_asian_mc`.
- Verified: matches `average_strike_asian_mc` for call and put within MC error,
  at roughly a seventh of the plain-MC standard error at equal points; prices
  positive; SE < 0.02.

## [1.205.1] - 2026-09-10

### Tests
- Trimmed the two `bermudan_min_put_lsm_greeks` sign-check fast tests (each
  re-prices five times) from 20-25 steps / 15000-30000 paths to 12-15 / 5000-8000;
  the delta signs and the `delta2 < delta1` ordering are robust across seeds at
  the smaller size. Fast gate ~60s -> ~30s.

## [1.205.0] - 2026-09-10

### Added
- `bermudan_min_put_lsm_greeks` (in `lsm.py`): deltas, own-gammas, and
  cross-gamma of an American min-put (worst-of protective put) by
  common-random-number bumps on `bermudan_min_put_lsm`. Same-seed repricing
  shares the Brownian shocks so the finite differences are low-variance; the LSM
  regression is re-fit at each bump. Completes the two-asset LSM Greek family
  (max-call, spread, min-put).
- Verified: both spot deltas are negative (a higher spot lifts the min and
  shrinks the put); the lower-starting asset carries the larger-magnitude delta
  (it is more often the min); the price exceeds the European
  `worst_of_put_closed` by the early-exercise premium; reproducible.

## [1.204.0] - 2026-09-10

### Added
- `sobol_cliquet_rqmc` (in `sobol.py`): randomized-QMC capped cliquet (ratchet)
  note with an honest standard error -- the same product as `capped_cliquet_mc`
  (per-period returns clipped to `[local_floor, local_cap]`, running sum clipped
  to `[global_floor, global_cap]`). Brownian motion at the reset times comes from
  one Sobol point via a bridge on the reset grid; the per-period standardized
  shock is the bridge increment over `sqrt(dt_i)`, randomized by a per-dimension
  Cranley-Patterson rotation.
- Verified: matches `capped_cliquet_mc` with full caps and with uncapped-local /
  no-global settings within MC error; a tighter local cap lowers the value; the
  across-randomization SE is tight (< 0.01).

## [1.203.0] - 2026-09-10

### Added
- `sobol_double_knockout_rqmc` (in `sobol.py`): randomized-QMC double-knockout
  (corridor) option with an honest standard error. Pays the vanilla payoff only
  if the spot stays strictly inside `(lower, upper)` at every monitoring date,
  else the cash `rebate`. Bridge construction + per-dimension Cranley-Patterson
  rotation; the discrete analogue of `double_knockout_mc`.
- Verified: matches `double_knockout_mc` within MC error; a very wide corridor
  approaches the vanilla call; a narrower corridor is worth strictly less.

## [1.202.0] - 2026-09-10

### Added
- `sobol_autocallable_rqmc` (in `sobol.py`): randomized-QMC autocallable
  structured note with an honest standard error -- the same product as
  `autocallable_mc` (early redemption with accrued coupons at each observation,
  down-and-in protection at maturity). Each path's Brownian values at the
  observation dates come from one Sobol point via a Brownian bridge built
  directly on the (possibly non-uniform) observation grid (`_bridge_on_times`),
  randomized by a per-dimension Cranley-Patterson rotation.
- Verified: matches `autocallable_mc` with and without a protection barrier
  within MC error; the price rises with the coupon; the across-randomization SE
  is tight (< 0.01).

## [1.201.0] - 2026-09-10

### Added
- `sobol_barrier_digital_rqmc` (in `sobol.py`): randomized-QMC barrier-contingent
  cash-or-nothing digital with an honest standard error. Pays `cash` iff the
  option finishes in the money AND the barrier condition holds over the `n_steps`
  monitoring dates (`up-in`/`down-in`/`up-out`/`down-out`). Bridge construction +
  per-dimension Cranley-Patterson rotation; the discrete analogue of
  `barrier_digital_mc`.
- Verified: matches `barrier_digital_mc` across all four barrier types within MC
  error; a knock-in plus its knock-out sum to the plain `cash_or_nothing` digital
  at the same seed; the value scales linearly in `cash`.

## [1.200.1] - 2026-09-10

### Tests
- Trimmed the two `bermudan_spread_lsm_greeks` sign-check fast tests (each
  re-prices five times) from 20 steps / 15000 paths to 12 / 5000; they assert
  only delta signs, which are robust at the smaller size. Fast gate ~35s -> ~20s.

## [1.200.0] - 2026-09-10

### Added
- `bermudan_spread_lsm_greeks` (in `lsm.py`): deltas, own-gammas, and
  cross-gamma of an American spread option by common-random-number bumps on
  `bermudan_spread_lsm`. Same-seed repricing shares the Brownian shocks so the
  finite differences are low-variance; the LSM regression is re-fit at each bump.
- Verified: without dividends the two spot deltas match the European Kirk
  `spread_greeks` deltas within noise (slow test); the spread-call long-leg
  delta is positive and the short-leg delta negative, with the signs flipping
  for the put; reproducible under a fixed seed.

## [1.199.0] - 2026-09-10

### Added
- `sobol_barrier_rqmc` (in `sobol.py`): randomized-QMC discretely-monitored
  single-barrier vanilla (down/up, in/out, with rebate) and an honest standard
  error. Normals come from an `n_steps`-dim Sobol point through the Brownian
  bridge, randomized by a per-dimension Cranley-Patterson rotation; the discrete
  analogue of `barrier_mc(brownian_bridge=False)`.
- Verified: matches the discretely-monitored `barrier_mc` across all four
  barrier types within MC error; a knock-in plus its knock-out sum to the
  vanilla call at the same seed (they partition every path); the knock-out is
  below the vanilla.

## [1.198.0] - 2026-09-10

### Added
- `two_asset_digital_greeks` (in `multiasset.py`): Greeks of a two-asset
  correlated digital by finite differences on the exact `two_asset_digital`
  closed form (no Monte Carlo noise) -- the two spot deltas, own-gammas,
  cross-gamma, and the correlation sensitivity `corr_vega`.
- Verified: a both-above digital has positive spot deltas and a positive
  `corr_vega` (the two in-the-money events move together), a mixed above/below
  digital has negative `corr_vega`, and summed over the four exhaustive quadrants
  the correlation sensitivity is zero (total probability is `rho`-independent);
  Greeks scale linearly in `cash`; a `below` condition flips that spot delta's
  sign.

## [1.197.0] - 2026-09-10

### Added
- `bermudan_min_put_lsm` (in `lsm.py`): American put on the minimum of two assets
  `max(K - min(S1_T, S2_T), 0)` by Longstaff-Schwartz -- the worst-of protective
  put. Two correlated GBMs; the continuation value is regressed on a quadratic
  basis in both spots plus the running min
  `{1, S1, S2, S1^2, S2^2, S1 S2, min(S1,S2)}` over the in-the-money paths.
- Verified: exceeds the European Stulz `worst_of_put_closed` by a positive
  early-exercise premium (puts carry early-exercise value even without
  dividends, ~0.42 here); is worth at least its intrinsic on the min when deep
  in the money; reproducible under a fixed seed.

### Audit
- Swept every test for Monte Carlo functions called without a `seed` whose
  result feeds an assertion (the class of the `test_pde2d_american` flake fixed
  in 1.196.1). All remaining seedless calls are `pytest.raises(ValueError)`
  argument-validation checks (no sampling) or pass the seed via `**kw`; no
  further nondeterministic references found.

## [1.196.1] - 2026-09-10

### Fixed
- `test_pde2d_american.py` compared the deterministic two-asset ADI price to the
  *Monte Carlo* `best_of_call` / `worst_of_call` (seed unset), so the reference
  itself was random and its noise occasionally breached the tolerance -- the
  source of the rare "1 failed" seen in full-suite runs. Switched both to the
  exact Stulz `best_of_call_closed` / `worst_of_call_closed`; the checks are now
  deterministic (and the worst-of one no longer needs `-m slow`).

## [1.196.0] - 2026-09-10

### Changed
- The test suite now runs across all cores by default via `pytest-xdist`
  (`addopts = "-n auto"` in `pyproject.toml`; `pytest-xdist>=3` added to the
  `test` extra). The ~1700 tests are independent with no shared state, so
  distribution is safe; the fast gate drops from ~55s to ~12s on a 24-core box
  and CI gets a proportional speedup. Run serially with `-n0`.

### Tests
- Reduced the six RQMC lookback structural checks (floating and fixed strike)
  from 4096 paths / 24 randomizations to 2048 / 12; the sign and ordering
  assertions and the SE bound (actual SE ~0.006 vs a 0.05 threshold) hold
  comfortably at the smaller size.

## [1.195.0] - 2026-09-10

### Added
- `sobol_fixed_lookback_rqmc` (in `sobol.py`): randomized-QMC discretely-
  monitored fixed-strike lookback with an honest standard error. Call payoff
  `max(max_i S_{t_i} - K, 0)`, put payoff `max(K - min_i S_{t_i}, 0)` over `S_0`
  and the `n_steps` monitoring dates; normals from an `n_steps`-dim Sobol point
  through the Brownian bridge, randomized by a per-dimension Cranley-Patterson
  rotation.
- Verified: the discrete price is below the continuously-monitored Conze-
  Viswanathan `fixed_strike_lookback` for call and put, and rises toward it as
  monitoring frequency grows (n=2 -> 6: 12.2 -> 14.6 vs 19.2 continuous); the
  call dominates the plain vanilla.

## [1.194.1] - 2026-09-10

### Tests
- Trimmed the fast suite ~66s -> ~50s by shrinking the heavy two-asset LSM/RQMC
  structural checks (each stores full two-asset paths and re-prices several
  times): the max-call cross-gamma check 15 steps / 8000 paths -> 12 / 5000 (its
  sign is robust across seeds there), the delta-interval check 12 / 3000 ->
  10 / 2500, and the American-spread lower-bound check 25 / 15000 -> 15 / 6000.
  The accuracy cross-checks stay under `-m slow`.

## [1.194.0] - 2026-09-10

### Added
- `sobol_lookback_rqmc` (in `sobol.py`): randomized-QMC discretely-monitored
  floating-strike lookback with an honest standard error. Call payoff
  `S_T - min_i S_{t_i}`, put payoff `max_i S_{t_i} - S_T` over the `n_steps`
  monitoring dates; normals come from an `n_steps`-dim Sobol point through the
  Brownian bridge, randomized by a per-dimension Cranley-Patterson rotation so
  `n_rand` shifts give a genuine SE.
- Verified: the discrete price is below the continuously-monitored
  Goldman-Sosin-Gatto `floating_strike_lookback` (fewer sampling dates see less
  extreme highs/lows), for both call and put, and rises toward it as the
  monitoring frequency increases (n=2 -> 6: 11.9 -> 13.8 vs 17.2 continuous);
  the across-randomization SE is tight (< 0.05).

## [1.193.0] - 2026-09-10

### Added
- `bermudan_spread_lsm` (in `lsm.py`): American spread option
  `max(S1_T - S2_T - K, 0)` (call) or `max(K - (S1_T - S2_T), 0)` (put) by
  Longstaff-Schwartz. Two correlated GBMs; the continuation value is regressed
  on a quadratic basis in both spots plus the spread
  `{1, S1, S2, S1^2, S2^2, S1 S2, S1 - S2}` over the in-the-money paths.
- Verified: without dividends it matches the European Kirk `spread_option`
  (11.83); an 8% dividend on the long leg produces a positive early-exercise
  premium (~0.75); it is never cheaper than the European; call and put prices
  are positive and reproducible.

## [1.192.0] - 2026-09-10

### Added
- `bermudan_max_call_lsm_greeks` (in `lsm.py`): deltas, own-gammas, and
  cross-gamma of an American max-call by common-random-number bumps on
  `bermudan_max_call_lsm`. Repricing at bumped spots with the same seed shares
  the Brownian shocks, so the finite differences are low-variance; the LSM
  regression is re-fit at each bump. Deltas are reliable; the gammas (second
  differences over a re-fit regression) are indicative and need many paths.
- Verified: without dividends the two spot deltas match the exact European
  `rainbow_greeks` deltas within noise; both deltas stay in `(0,1)`; the
  cross-gamma is robustly negative across seeds (the two spots are substitutes
  in a max payoff); reproducible under a fixed seed.

## [1.191.1] - 2026-09-10

### Tests
- Trimmed the `bermudan_max_call_lsm` single-asset-lower-bound fast test from 30
  steps / 20000 paths to 20 / 8000; the bound has a ~7-point margin, so the
  smaller sample is comfortably safe. Fast gate back to ~51s.

## [1.191.0] - 2026-09-10

### Added
- `bermudan_max_call_lsm` (in `lsm.py`): American call on the maximum of two
  assets by Longstaff-Schwartz -- the classic two-asset early-exercise
  benchmark. Simulates two correlated GBMs and regresses the continuation value
  on a quadratic basis in both spots plus the running max
  `{1, S1, S2, S1^2, S2^2, S1 S2, max(S1,S2)}` over the in-the-money paths at
  each date.
- Verified: without dividends it matches the European Stulz
  `best_of_call_closed` (little early-exercise value); with a 6% dividend on
  both assets it exceeds the European by a positive early-exercise premium; it
  dominates a single-asset vanilla and is reproducible under a fixed seed.

## [1.190.1] - 2026-09-10

### Tests
- Trimmed the fast suite: the three `heston_mc_greeks` structural checks (each
  runs 9 QE re-prices per Greek set, ~14s combined) drop to smaller sizes -- the
  reproducibility check to 12 steps / 1500 paths and the two positivity/interval
  checks to 15 / 3000. They only assert determinism, a sign, or a `(0,1)` range,
  not accuracy; the Fourier-FD accuracy cross-checks stay under `-m slow`.

## [1.190.0] - 2026-09-10

### Added
- `rainbow_greeks` (in `multiasset.py`): Greeks of a best-of/worst-of rainbow
  option by finite differences on the exact Stulz closed forms (no Monte Carlo
  noise). Returns the two spot deltas, the two own-gammas, the cross-gamma
  `d2V/dS1 dS2`, and the correlation sensitivity `dV/drho`, for calls and puts.
- Verified: the best-of and worst-of call deltas in each asset sum to the
  single-asset Black-Scholes delta (differentiating the Stulz identity
  `C_max + C_min = c1 + c2`), to 1e-4; the max-call `corr_vega` is negative and
  the min-call's positive with the two summing to zero; own-gammas positive;
  worst-of put deltas negative.

## [1.189.0] - 2026-09-10

### Added
- `heston_mc_greeks` (in `heston_mc.py`): Heston Greeks by common-random-number
  finite differences on the QE Monte Carlo. Repricing at bumped inputs with the
  same seed shares the random draws, so the bumped difference is dominated by the
  true sensitivity rather than MC noise. Returns `price`, `delta`, `gamma`,
  `vega_v0` (initial-variance), `vega_theta` (long-variance), `volvol` (dV/dxi),
  and `rho_sens` (dV/drho) -- the sensitivities that matter for a stochastic-vol
  book.
- Verified: at 120000-180000 paths every Greek matches a finite difference of the
  exact Fourier `heston_price` (delta 0.703, gamma 0.0191, vega_v0 47.6,
  vega_theta 48.6, volvol -2.08 vs -2.17, rho 0.20); delta stays in (0,1), gamma
  and both variance vegas are positive.

## [1.188.0] - 2026-09-10

### Added
- `two_asset_gap_option` (in `multiasset.py`): a gap option on asset 1 gated by
  asset 2, separating the trigger strike from the payoff strike. The asset-1 gap
  payoff `(S1_T - K_payoff) 1[S1_T > K_trigger]` (call) fires only if asset 2
  clears its barrier. Decomposes into the two two-asset digitals with the
  trigger strike setting the asset-1 condition and the payoff strike scaling the
  cash leg, so it is an exact bivariate-normal closed form.
- Verified: matches a correlated-GBM Monte Carlo for call and put; reduces to
  `correlation_option` when `K_payoff = K_trigger`; and a larger payoff strike
  lowers the call.

## [1.187.0] - 2026-09-10

### Added
- `correlation_option` (in `multiasset.py`): a two-asset correlation option -- a
  vanilla on asset 1 that pays only if asset 2 satisfies a barrier condition
  (`above`/`below` K2). Decomposes exactly into the two two-asset digitals:
  `call = AoN(S1>K1, cond2) - K1 CoN(S1>K1, cond2)` and
  `put = K1 CoN(S1<K1, cond2) - AoN(S1<K1, cond2)`, where `AoN` is
  `two_asset_asset_or_nothing` and `CoN` is `two_asset_digital`, so it is a
  closed form.
- Verified: matches a correlated-GBM Monte Carlo for call, put, and a below
  barrier; reduces to the plain Black-Scholes vanilla on S1 when the barrier is
  always met; and the `above` + `below` gated calls sum to that vanilla (1e-9).

## [1.186.0] - 2026-09-10

### Added
- `two_asset_asset_or_nothing` (in `multiasset.py`): asset-or-nothing digital
  paying the first asset's terminal value `S1_T` iff both conditions hold (each
  `above`/`below` its strike). Priced under the asset-1 (share) measure, where
  asset 1's drift gains `sigma1^2` and asset 2's shock inherits an extra
  `rho sigma1 sqrt(t)`: `S1 e^{-q1 t} M(s1 a1, s2 a2; s1 s2 rho)`.
- Verified: all four quadrants match a correlated-GBM Monte Carlo; they sum to
  the discounted forward `S1 e^{-q1 t}` (asset 1 is always delivered on some
  quadrant) to 1e-8, including under a dividend yield; and at `rho = 0` the
  price factorizes into the single-asset `asset_or_nothing` on S1 times the
  risk-neutral probability that S2 clears K2.

## [1.185.0] - 2026-09-10

### Added
- `two_asset_digital` (in `multiasset.py`): exact closed form for a cash-or-
  nothing digital on two correlated assets. Pays `cash` iff both single-asset
  conditions hold (each `above` or `below` its strike). Price is
  `cash e^{-rt} M(s1 d1, s2 d2; s1 s2 rho)` with `di` the usual `d2`,
  `si = +/-1` for above/below, and `M` the standardized bivariate-normal CDF;
  flipping a condition flips that `d`'s sign and the correlation.
- Verified: all four quadrant prices match a correlated-GBM Monte Carlo, sum to
  `e^{-rt}` (exhaustive) to 1e-9, factorize into the product of two single-asset
  `cash_or_nothing` digitals at `rho = 0`, and scale linearly in `cash`.

## [1.184.0] - 2026-09-10

### Added
- `best_of_put_closed` / `worst_of_put_closed` (in `multiasset.py`): exact
  closed forms for rainbow puts on the maximum / minimum of two assets. By
  rainbow put-call parity `P = C - disc E[chosen] + K e^{-rt}`, built on the
  call closed forms and `disc E[min] = S2 e^{-q2 t} - exchange_option(S2, S1)`
  (`min(a,b) = b - max(b-a,0)`), with `disc E[max]` the two forwards less the
  discounted expected min.
- Verified: the put identity `P_min + P_max = p(S1) + p(S2)` holds to 1e-9;
  both match their Monte Carlo counterparts within MC error; put-on-min exceeds
  put-on-max, and put-on-max is below each single-asset vanilla put.

## [1.183.0] - 2026-09-10

### Added
- `best_of_call_closed` / `worst_of_call_closed` (in `multiasset.py`): exact
  Stulz (1982) closed forms for rainbow calls on the maximum / minimum of two
  assets, replacing Monte Carlo (`best_of_call` / `worst_of_call`) with an
  analytic price. Built on the spread vol
  `sigma = sqrt(sigma1^2 - 2 rho sigma1 sigma2 + sigma2^2)` and the standardized
  bivariate-normal CDF; the call-on-min is priced directly and the call-on-max
  from the Stulz identity `C_max + C_min = c(S1) + c(S2)`.
- Verified: the identity holds to 1e-9; both match their Monte Carlo
  counterparts within MC error (best 17.150 vs 17.16, worst 4.574 vs 4.58), also
  under negative correlation with dividends; the max-call dominates each vanilla
  and the min-call.

## [1.182.0] - 2026-09-10

### Added
- `basket_option_lhs_mc` (in `montecarlo.py`): two-asset basket option
  `max(w1 S1 + w2 S2 - K, 0)` by Latin hypercube sampling -- the analogue of
  `spread_option_lhs_mc` for a weighted-sum payoff. Each driving normal is
  stratified into `n_paths` equiprobable bins with independently permuted bin
  orders, mapped through the inverse normal CDF; asset 2's shock is correlated
  by `corr z1 + sqrt(1 - corr^2) z2`.
- Verified: agrees with the Levy moment-matched (approximate) `basket_option`
  for call and put (to a few cents), reduces to a Black-Scholes call when all
  weight is on one asset, and -- wrapped in `replicated_mc` for an honest
  across-seed SE -- places the moment-match price within a few SEs. As with any
  LHS estimator the per-run `std_error` overstates the true error; documented.

## [1.181.1] - 2026-09-10

### Tests
- Trimmed the fast suite ~50s -> ~47s. The Heston-CV standard-error-ordering
  check drops from 30000 paths / 100 steps to 8000 / 30 (a robust inequality,
  not an accuracy claim), and the dividend-yield Fourier match moved under
  `-m slow` with the other `heston_cv_mc` accuracy cross-checks. Together these
  removed ~19s of QE simulation from the fast gate.

## [1.181.0] - 2026-09-10

### Added
- `heston_cv_mc` (in `heston_mc.py`): Heston QE Monte Carlo with the underlying
  as a control variate. Andersen's QE step is martingale-corrected, so the
  discounted terminal spot `Y = e^{-rt} S_T` has the known mean `S0 e^{-qt}`;
  the controlled estimator `X - beta (Y - E[Y])` with the regression-optimal
  `beta = Cov(X,Y)/Var(Y)` cuts the standard error at no bias. Same QE variance
  step, `K0..K4` asset constants, and antithetic draws as `heston_qe_mc`.
- Verified: matches the Fourier `heston_price` for call and put (also under a
  dividend yield) within MC error; at equal path count the standard error is
  ~0.56x the plain `heston_qe_mc` (the spot/payoff correlation is weaker under
  stochastic vol than in Black-Scholes, so the reduction is milder than the
  lognormal `european_cv_mc`).

## [1.180.1] - 2026-09-10

### Tests
- Trimmed the fast suite ~50s -> ~46s. The two RQMC-Asian cross-checks use a
  120000-path (was 400000) control-variate reference and n_rand=16 (was 24); the
  vol-surface example builds its SVI calibration once via a module-scoped fixture
  instead of three times; and the local-vol-MC positivity/smoke checks drop from
  100 steps / 20000 paths to 40 / 8000. Accuracy cross-checks stay under `-m slow`.

## [1.180.0] - 2026-09-10

### Fixed
- `arithmetic_asian_mc` control variate was biased. It used the *continuous*
  geometric-Asian closed form (`geometric_asian`, `sigma/sqrt(3)`) as the
  control, but the Monte Carlo averages over `n_steps` *discrete* dates, so the
  control's known expectation did not match the simulated geometric payoff. The
  CV price was biased low -- e.g. 5.76 vs a true ~6.55 at `n_steps=6` -- while
  reporting a tiny SE. Existing tests used `n_steps=50`/`250` (near-continuous)
  and a loose `3*plain_SE` band, so they missed it. Added
  `_discrete_geometric_asian` (exact closed form of the discretely-monitored
  geometric average, `E[ln G] = ln S0 + (b-sig^2/2) dt (n+1)/2`,
  `Var[ln G] = sig^2 dt (n+1)(2n+1)/(6n)`) and use it as the control. CV and
  no-CV prices now agree at every `n_steps`; a tight `n_steps=6` regression test
  guards it.

### Added
- `sobol_asian_rqmc` (in `sobol.py`): randomized-QMC arithmetic Asian with an
  honest standard error -- the multi-dimensional analogue of
  `sobol_european_rqmc`. Each path's normals come from an `n_steps`-dimensional
  Sobol point through the Brownian bridge; a per-dimension Cranley-Patterson
  rotation randomizes the point set, so `n_rand` shifts give i.i.d. QMC estimates
  whose spread is a genuine SE. Cross-checks the (now-fixed) control-variate
  `arithmetic_asian_mc` for call and put; ~0.06x the plain-MC SE at equal points.

## [1.179.0] - 2026-09-10

### Added
- `sobol_european_rqmc` (in `sobol.py`): randomized-QMC European price with an
  honest standard error. Plain Sobol QMC returns a single deterministic number
  with no error bar; this applies a Cranley-Patterson rotation -- shifting the
  whole Sobol point set by a random `U ~ Uniform[0,1)` mod 1 -- which preserves
  the low discrepancy but makes each estimate an unbiased draw. `n_rand`
  independent shifts give i.i.d. QMC estimates whose spread is a genuine SE.
  Returns the mean price, the across-randomization SE, and `n_paths` = total
  points (`n_rand * n_paths`).
- Verified: matches the closed-form Black-Scholes call and put (also under a
  dividend carry); at equal total points the RQMC standard error is ~0.05x the
  plain `european_mc` for this smooth 1-D integral.

## [1.178.0] - 2026-09-10

### Added
- `replicated_mc` (in `montecarlo.py`): batched-replication honest standard
  error for variance-reduced estimators. Stratified, Latin hypercube, and QMC
  estimators draw dependent samples, so the plain i.i.d. `std_error` they report
  is not the truth (see the note on `spread_option_lhs_mc`). Running the whole
  estimator `n_batches` times with distinct seeds gives independent batch means
  whose spread is an unbiased standard error. Takes any callable
  `seed -> MCResult | float`; returns the mean price, the across-batch SE
  `s/sqrt(n_batches)`, and `n_paths = n_batches`.
- Verified: for a plain `european_mc` the across-batch SE matches a single-run
  i.i.d. SE / sqrt(n_batches) (0.032 vs 0.031, calibrated); the mean is unbiased
  vs Black-Scholes for a stratified estimator; and it exposes the LHS spread
  reduction (LHS honest SE < 0.6x a plain two-asset MC's honest SE) that the
  naive per-run SE hides.

## [1.177.0] - 2026-09-10

### Added
- `spread_option_lhs_mc` (in `montecarlo.py`): two-asset spread option
  `max(S1 - S2 - K, 0)` by Latin hypercube sampling. Each of the two driving
  normals is stratified into `n_paths` equiprobable bins with one draw per bin,
  the two dimensions' bin orders are independently permuted, and the stratified
  uniforms are mapped through the inverse normal CDF; the second asset's shock
  is correlated by `rho z1 + sqrt(1 - rho^2) z2`.
- Verified: matches the Kirk `spread_option` and, at `K = 0`, the exact Margrabe
  `exchange_option`. The LHS samples are not independent, so the reported
  `std_error` (plain i.i.d. formula) overstates the true error; the genuine gain
  is in the across-seed spread -- at 4000 paths the RMSE vs Kirk is ~0.10 versus
  ~0.27 for a plain two-asset draw (~2.6x), documented on the function.

## [1.176.0] - 2026-09-10

### Added
- `european_stratified_mc` (in `montecarlo.py`): European price by stratified
  sampling of the terminal normal. The driving normal is split into `n_strata`
  equiprobable strata in probability space; `n_per` uniforms are drawn within
  each and mapped through the inverse normal CDF, spreading draws evenly and
  removing the clustering that inflates plain MC variance. With equiprobable
  strata the estimate is the average of the per-stratum means and its variance
  is `(1/n_strata^2) sum_i s_i^2 / n_per`.
- Verified: matches the closed-form Black-Scholes call and put (also deep OTM
  and under a dividend carry); at 10000 paths (400 strata x 25) the standard
  error is ~0.05x the plain `european_mc` for the smooth call payoff.

## [1.175.0] - 2026-09-10

### Added
- `european_is_adaptive_mc` (in `montecarlo.py`): importance sampling with a
  pilot-tuned optimal shift. The strike-centring shift of `european_is_mc` is
  near-optimal for a digital but not for a vanilla, whose payoff keeps growing
  past the strike and pulls the best shift further OTM. A change of measure
  gives the shifted second moment from plain draws,
  `M(mu) = E_0[payoff(z)^2 exp(-mu z + mu^2/2)]`, so one `N(0,1)` pilot sample
  scores every candidate `mu` on a grid at negligible cost; the main run samples
  at the variance-minimising `mu` and stays unbiased via the likelihood ratio.
- Verified: matches the closed-form Black-Scholes call and put (deep OTM, ATM,
  under carry); for a K=160 deep-OTM call the tuned shift gives ~0.87x the
  standard error of the fixed strike-centring shift at equal main-run paths.

## [1.174.0] - 2026-09-10

### Added
- `digital_is_mc` (in `montecarlo.py`): cash-or-nothing digital price by
  importance sampling. A deep-OTM digital is even harder to simulate plainly
  than a vanilla -- the 0/`cash` indicator has relative SE that blows up like
  `sqrt((1-p)/p)` for a small hit probability `p`. Shifting the terminal normal
  to `N(mu, 1)` and reweighting by `L(z) = exp(-mu z + mu^2/2)` moves paths into
  the money unbiased; the default shift puts the mean draw on the strike
  boundary so about half the paths pay.
- Verified: matches the closed-form `cash_or_nothing` (deep OTM call and put,
  ATM, cash scaling, dividend carry); for a K=160 deep-OTM digital the standard
  error is ~0.17x a plain indicator estimator at equal paths.

## [1.173.0] - 2026-09-10

### Added
- `european_is_mc` (in `montecarlo.py`): European price by importance sampling,
  for deep out-of-the-money options where a plain simulation wastes almost every
  path. Draws the terminal normal from a shifted mean `N(mu, 1)` to push mass
  into the money and corrects with the likelihood ratio
  `L(z) = exp(-mu z + mu^2/2)`. The default shift centres the terminal log-spot
  on the strike, `mu* = (ln(K/S0) - (b - sig^2/2) t) / (sig sqrt(t))`.
- Verified: matches the closed-form Black-Scholes call and put (deep OTM and
  ATM, under a dividend carry, and with an explicit shift); for a K=160 deep-OTM
  call the standard error is ~0.086x the plain `european_mc` at equal paths.

## [1.172.0] - 2026-09-10

### Added
- `lr_digital_greeks` (in `mc_greeks.py`): delta, vega, and gamma of a
  cash-or-nothing digital by the likelihood-ratio method. The digital payoff is
  discontinuous, so the pathwise method is undefined for every Greek; the LR
  method differentiates the log-normal density, so the same one-step
  Black-Scholes score weights as `lr_greeks` (delta `Z/(S0 sig sqrt t)`, vega
  `(Z^2-1)/sig - Z sqrt t`, gamma `(Z^2 - Z sig sqrt t - 1)/(S0^2 sig^2 t)`)
  apply unchanged to the digital.
- Verified: delta and gamma match the analytic `digital_greeks`, vega matches a
  sigma-bump of `cash_or_nothing`, within MC error for both call and put (e.g.
  call delta 0.0274 vs 0.0274, vega -0.478 vs -0.479, gamma -0.00048 vs
  -0.00048); call delta positive, put delta negative.

## [1.171.1] - 2026-09-10

### Tests
- Trimmed the fast suite ~49s -> ~35s. The LSM-Greeks sign/price-field checks
  drop from 25 steps / 8000 paths to 20 / 4000, and the Asian pathwise-vega
  positivity check from 40 steps / 40000 paths to 20 / 15000. These assert
  only a sign or an exact price-field match, not accuracy; the accuracy
  cross-checks stay under `-m slow`, unchanged.

## [1.171.0] - 2026-09-10

### Added
- `barrier_lr_delta` (in `mc_greeks.py`): delta of a discretely-monitored
  single-barrier option by the likelihood-ratio method. The knock-out/knock-in
  payoff is discontinuous in the spot, so the pathwise method fails; the LR
  method uses the score of the path density. The initial spot enters only
  through the mean of the first log-increment, so the score is
  `Z_1 / (S0 sig sqrt(dt))` and `delta = E[disc_payoff * Z_1 / (S0 sig sqrt(dt))]`.
- Verified: matches a common-random-number finite difference of the
  discretely-monitored `barrier_mc` (brownian_bridge=False) within MC error for
  a down-out call (0.784 vs 0.788) and an up-in put (0.089 vs 0.092); down-out
  call delta is positive.

## [1.170.0] - 2026-09-10

### Added
- `european_cv_mc` (in `montecarlo.py`): European Monte Carlo combining both
  variance-reduction techniques -- antithetic sampling and a control variate.
  The discounted terminal spot `Y = e^{-rt} S_T` has the known mean
  `S0 e^{(b-r)t}` and is correlated with the payoff, so the estimator
  `X - beta (Y - E[Y])` with the regression-optimal `beta = Cov(X,Y)/Var(Y)`
  cuts the variance. Antithetic pairs are averaged into one sample before the
  control is applied.
- Verified: matches the closed-form Black-Scholes call and put within MC error
  (also under a dividend carry `b = r - q`); at equal path count the standard
  error is ~0.24x the plain `european_mc` for the ITM call and ~0.42x for the
  put, and under half at the money.

## [1.169.1] - 2026-09-10

### Tests
- Trimmed the fast suite ~49s -> ~35s: the two Heston pathwise-delta
  sign/range checks run at n_steps=40 / n_paths=6000 instead of 50/20000 (1.8s
  each -> under 0.7s). The Fourier-FD accuracy cross-checks stay under `-m slow`.

## [1.169.0] - 2026-09-10

### Added
- `heston_pathwise_delta` (in `heston_mc.py`): Heston delta by the pathwise
  method. In the Andersen-QE simulation the initial spot enters only as the
  additive `ln S0` in the terminal log-price, so `S_T = S0 e^Y` with `Y`
  independent of `S0` and the pathwise delta is exact:
  `delta = e^{-rt} E[1_{S_T>K} S_T/S0]` (put: `-1_{S_T<K}`).
- Verified: the call and put deltas match a finite-difference of the exact
  Fourier `heston_price` within Monte Carlo error; the call delta is in (0, 1)
  and the put delta is negative.

## [1.168.0] - 2026-09-10

### Added
- `smoothed_digital_delta` (in `mc_greeks.py`): pathwise delta of a
  cash-or-nothing digital via a call-spread smoothing of the discontinuous
  indicator (ramp of relative width `eps_rel`), making the payoff Lipschitz so
  the pathwise derivative is well-defined. A single-pass, model-agnostic
  alternative to the likelihood-ratio digital delta.
- Verified: it matches the analytic digital delta; a narrower spread lowers the
  smoothing bias (converging as `eps_rel -> 0`, at the cost of higher variance
  from the `1/eps` ramp -- the classic bias/variance tradeoff, documented); call
  delta positive, put delta negative.

## [1.167.1] - 2026-09-10

### Tests
- Trimmed the fast suite: the two fast LSM-Greeks checks (put-delta sign, price
  field consistency) run at n_steps=25 / n_paths=8000 instead of 40/40000 (9s
  each -> under 2s), cutting ~14s off the fast run. The delta/gamma accuracy
  cross-checks against the binomial tree still run under `-m slow`.

## [1.167.0] - 2026-09-10

### Added
- `bermudan_lsm_greeks` (in `lsm.py`): delta and gamma of a Bermudan/American
  LSM price by common-random-number bumps -- the option is repriced at
  `S(1 +/- h)` on the same seeded path stream, so the finite differences are
  low-variance. Returns price, delta and gamma.
- Verified: delta matches a 2000-step binomial tree (~1e-2); gamma is positive
  and in a loose band around the tree gamma (a second difference over a re-fit
  regression is only indicative and needs many paths -- documented). Delta
  converges to the tree value as paths grow (-0.40 -> -0.412).

## [1.166.0] - 2026-09-10

### Added
- `mixed_gamma` (in `mc_greeks.py`): European gamma by the mixed
  pathwise-likelihood-ratio estimator. It differentiates the pathwise delta
  payoff `D = e^{-rt} 1_{S_T>K} S_T/S0` w.r.t. `S0` -- through both the density
  (LR weight `Z/(S0 sigma sqrt t)`) and the explicit `1/S0` factor -- giving
  `gamma = E[D (Z/(S0 sigma sqrt t) - 1/S0)]`. Well-defined despite the
  indicator (pure pathwise gamma is not) and lower variance than a double
  likelihood-ratio.
- Verified: matches the Black-Scholes gamma; its standard error is ~4x smaller
  than the LR gamma at equal paths; put gamma equals call gamma; positive.

## [1.165.0] - 2026-09-10

### Added
- Multi-level Monte Carlo (new `mlmc.py`): `mlmc_asian` prices a fixed-strike
  arithmetic Asian by the Giles (2008) telescoping estimator
  `E[P_L] = E[P_0] + sum_l E[P_l - P_{l-1}]`, with each correction from coupled
  fine/coarse paths sharing Brownian increments (each coarse step sums `M` fine
  ones). Deeper levels use fewer paths since their corrections have lower
  variance, cutting the cost to a target accuracy.
- Verified: the level corrections decay geometrically; the estimate converges
  to the Turnbull-Wakeman value as levels grow (err ~0.01 at 6 levels) and more
  levels reduce the discretisation bias; call and put prices are positive.

## [1.164.0] - 2026-09-10

### Added
- `asian_pathwise_vega` (in `mc_greeks.py`): pathwise vega of a fixed-strike
  arithmetic-average Asian option -- differentiating the payoff along each path
  w.r.t. `sigma`, with `dS_i/dsigma = S_i (W_i - sigma t_i)` and
  `dA/dsigma = mean_i dS_i/dsigma`. A lower-variance alternative to bumping for
  this path-dependent, Lipschitz payoff.
- Verified: it matches a central-difference bump of the arithmetic-Asian Monte
  Carlo within MC error, and agrees with a finite-difference of the
  Turnbull-Wakeman continuous-average closed form; call and put vega are
  positive.

## [1.163.0] - 2026-09-10

### Added
- `density_var_es` (new `density_var.py`): risk-neutral Value-at-Risk and
  Expected Shortfall of an option position whose horizon P&L is `pnl(S_T)`,
  computed by integrating the smile-implied density into a P&L distribution and
  taking the loss quantile and tail mean. VaR/ES returned as positive losses.
- Verified: a long-forward VaR matches the lognormal quantile; ES >= VaR; a long
  call's VaR is capped at the premium (max loss); a short call has a positive
  VaR with a larger ES; higher confidence gives a larger VaR.
- Bug fixed during validation: the ES accumulation broke on the first
  zero-density deep-tail atom (returning ~0); it now fills the whole tail mass
  regardless of individual atom probabilities.

## [1.162.0] - 2026-09-10

### Added
- `wasserstein_smiles` (in `density_metrics.py`): the 1-Wasserstein
  (earth-mover) distance between two smile-implied densities, computed as the L1
  gap between their CDFs `integral |F_p - F_q| dK` on a shared strike grid.
  Unlike the KL divergence it is a true metric (symmetric, triangle inequality)
  measured in price units.
- Verified: identical smiles give ~0; it is symmetric; positive for differing
  vols and growing with the vol gap; a skewed smile is a positive distance from
  a flat one; the triangle inequality holds.

## [1.161.0] - 2026-09-10

### Added
- `kl_divergence_smiles` (in `density_metrics.py`): the Kullback-Leibler
  divergence `KL(g_p || g_q) = integral g_p ln(g_p/g_q) dK` between two
  smile-implied risk-neutral densities (the second interpolated onto the first's
  grid). Measures how far one implied distribution sits from another -- two
  dates, two models, or implied vs a reference.
- Verified: identical smiles give ~0; differing vols give a positive,
  asymmetric, non-negative divergence that grows with the vol gap; a skewed
  smile diverges positively from a flat one.

## [1.160.0] - 2026-09-10

### Added
- Risk-neutral density tail/shape metrics (new `density_metrics.py`):
  `tail_probability` (`Q(S_T < L)` / `> U`, the undiscounted digital price),
  `density_entropy` (differential entropy of the terminal-spot density) and
  `expected_shortfall` (`E^Q[S_T | S_T` in a tail]), all from the
  Breeden-Litzenberger smile density.
- Verified: the upper-tail probability matches the lognormal `N(d2)`; the two
  tails sum to 1; a downward skew fattens the left tail; entropy increases with
  vol; the tail expected-shortfalls sit beyond their thresholds; an empty tail
  returns NaN.

### Fixed
- `density_grid_from_smile` could evaluate a call at a non-positive strike at the
  low end of a wide/high-vol grid (the finite-difference step exceeded the first
  strike). The FD step is now capped at `0.5 * K` so `K - h > 0` everywhere.

## [1.159.0] - 2026-09-10

### Added
- Breeden-Litzenberger risk-neutral density from a smile (new `rnd.py`):
  `risk_neutral_density_from_smile` gives `g(K) = e^{r t} d^2 C/dK^2`,
  `density_grid_from_smile` returns the density on a strike grid, and
  `price_payoff_from_density` prices any European payoff model-free by
  integrating it against the density.
- Verified: the density integrates to 1; it reprices calls, puts and digitals
  to Black-Scholes; the payoff `K` recovers the spot (discounted forward); the
  ATM density is positive; and a concave (arbitraging) smile produces a negative
  density in the wings.

## [1.158.1] - 2026-09-10

### Tests
- Trimmed the fast-suite runtime from ~45s to ~30s: the variance-term-structure
  tests build the strip at n_strikes=101, the LSM local-vol American>=European
  check runs at n_steps=25/n_paths=8000, the rBergomi auto/forced identity uses
  150 paths, and the Levy VG parameter recovery is marked slow. All accuracy
  claims still run under `-m slow`.

## [1.158.0] - 2026-09-10

### Added
- `variance_term_structure` (in `varswap.py`): the term structure of
  variance-swap strikes across expiries plus the forward (instantaneous)
  variance curve between them. Returns `(spot_var, forward_var)` where
  `forward_var[i] = (K_i t_i - K_{i-1} t_{i-1})/(t_i - t_{i-1})` by
  total-variance additivity.
- Verified: a flat term structure gives a constant curve (to the ~1e-6 strip
  truncation level); a rising vol term structure gives positive, increasing
  forward variances; the forward curve matches
  `forward_variance_swap_from_smile`; and `forward_var[0]` equals the first spot
  strike.

## [1.157.0] - 2026-09-10

### Added
- `moment_risk_premia` (new `moment_premium.py`): skewness- and
  kurtosis-risk premia -- the risk-neutral (BKM) skewness/excess-kurtosis the
  option smile implies versus what realized in the price history. Returns
  realized, implied and premium (implied - realized) for both moments, wiring
  `bkm_moments_from_smile` to the realized-moment estimators.
- Verified: a downward-skewed smile against a near-symmetric GBM history gives a
  negative skew premium (crash-protection demand) and a positive kurtosis
  premium; a flat smile gives a near-zero skew premium; the components are
  internally consistent.

## [1.156.0] - 2026-09-10

### Added
- `equity_premium_lower_bound` (in `vix.py`): Martin's (2013) model-free lower
  bound on the expected equity excess return, `Rf * SVIX^2`, computed from the
  option smile via the simple-variance index. Under the negative-correlation
  condition the annualized expected market excess return is bounded below by
  this quantity.
- Verified: a flat 20% vol gives a bound near `Rf * sigma^2`; a higher vol
  raises it; and it equals `Rf` times the SVIX variance exactly.

## [1.155.0] - 2026-09-10

### Added
- `bermudan_lsm_local_vol` (in `lsm.py`): Longstaff-Schwartz American/Bermudan
  pricing under a local-volatility surface. Same backward-induction regression
  as `bermudan_lsm`, but each Euler step uses `local_vol_fn(S, tau)` -- so it
  prices early-exercise options directly on a calibrated Dupire/SVI local-vol
  surface. Carry `b = r - q`.
- Verified: a flat `local_vol_fn` reproduces the constant-vol LSM price and a
  binomial tree; a genuinely skewed local vol matches an independent
  Crank-Nicolson American PDE (~0.03).

## [1.154.0] - 2026-09-10

### Added
- Monte Carlo Greeks without bumping (new `mc_greeks.py`): `lr_greeks` returns
  European delta, gamma and vega by the likelihood-ratio (Malliavin-flavoured)
  method -- `E[payoff * weight]` with density-derivative weights -- which works
  even for discontinuous payoffs; `pathwise_delta` is the lower-variance
  pathwise estimator for smooth payoffs; and `lr_digital_delta` gives the delta
  of a cash-or-nothing digital, where the pathwise method fails.
- Verified: LR delta/gamma/vega match Black-Scholes within Monte Carlo error;
  the pathwise delta matches and has a smaller standard error than the LR delta;
  and the LR digital delta matches the analytic (finite-difference) digital
  delta.

## [1.153.1] - 2026-09-10

### Tests
- Trimmed the fast-suite runtime from ~44s to ~32s: the Kim-Greeks structural
  checks (price equals the direct price, sign checks) run at n_steps=50 instead
  of 120; the Bermudan-swaption sign/monotonicity checks at n_paths=6000 instead
  of 20000; and the Heston-calibration parameter-validity check is marked slow.
  All accuracy claims still run under `-m slow`.

## [1.153.0] - 2026-09-10

### Added
- `svix_from_smile` (in `vix.py`): Martin's (2013) simple-variance index (SVIX).
  Unlike the VIX log-contract (`1/K^2` weights), the simple variance swap weights
  the OTM strip by the constant `1/F^2`, corresponding to the payoff
  `(S_T - F)^2 / F^2` with no log approximation -- jump-robust and a genuine
  lower bound on the equity premium. Reported as `100 * SVIX`.
- Verified: a flat smile returns ~`100 * sigma`; `SVIX = 100 * sqrt(variance)`;
  it is positive; and it differs from the VIX under a skew (the two coincide only
  to leading order).

## [1.152.0] - 2026-09-10

### Added
- Variance risk premium (new `vrp.py`): `realized_variance` computes annualized
  realized variance from a close series, and `variance_risk_premium` combines it
  with a supplied implied (variance-swap) variance into the additive premium
  (realized - implied), the ratio, and the vol-point premium
  (implied vol - realized vol).
- Verified: realized variance recovers `sigma^2` from a GBM path; implied above
  realized gives a negative VRP with a positive vol premium and ratio < 1;
  realized above implied flips the signs; the components are internally
  consistent.

## [1.151.0] - 2026-09-10

### Added
- Bakshi-Kapadia-Madan risk-neutral moments (new `bkm.py`):
  `bkm_moments_from_smile` extracts the model-free risk-neutral variance,
  skewness and excess kurtosis of the log-return from a smile `vol_fn(K)` via the
  quadratic/cubic/quartic option-spanning contracts (strike-weighted OTM
  strips). `skew_swap_from_smile` returns the risk-neutral skewness (a skew
  swap's fair value).
- Verified: a flat smile gives variance `sigma^2 t`, zero skewness and zero
  excess kurtosis; a downward (equity) skew gives negative skewness and positive
  excess kurtosis; an upward skew gives positive skewness; a steeper skew is more
  negative.

## [1.150.0] - 2026-09-10

### Added
- CBOE VIX-style fair volatility index (new `vix.py`): `vix_from_chain`
  implements the exact discrete CBOE formula -- the `dK/K^2`-weighted OTM strip
  minus the `(F/K0 - 1)^2` forward-correction term, reported as
  `100 * sqrt(variance)` -- and `vix_from_smile` builds the chain from a smile
  `vol_fn(K)`.
- Verified: a flat smile returns `VIX = 100 * sigma` exactly; a downward skew
  lifts the index above the ATM level (the fear premium); a manually-built
  flat-vol chain matches; and `VIX = 100 * sqrt(variance)`.

## [1.149.0] - 2026-09-10

### Added
- `forward_variance_swap_from_smile` (in `varswap.py`): fair forward-start
  variance-swap strike over `[t1, t2]` from the two expiries' smiles. Total
  variance is additive in time, so the accrued variance is
  `(K_var(t2) t2 - K_var(t1) t1)/(t2 - t1)` with each spot-starting leg replicated
  from its smile.
- Verified: a flat term structure returns `sigma^2`; `t1 = 0` recovers the spot
  variance swap exactly; an upward-sloping vol term structure gives a forward
  variance above the front; and the additivity identity holds to 1e-6.

## [1.148.0] - 2026-09-10

### Added
- `gamma_swap_from_smile` (in `varswap.py`): fair gamma-swap (price-weighted
  variance) strike from a smile `vol_fn(K)`. A gamma swap accrues `(S_t/S0)`-
  weighted realized variance, so by Carr-Lewis it is replicated by an option
  strip weighted `1/K` (vs the variance swap's `1/K^2`), scaled by
  `2 e^{rt}/(S0 t)`, with no forward remainder (the price-weighted log contract
  has none).
- Verified: a flat smile returns `sigma^2`; a denser strip reduces the error; a
  downward skew makes the gamma swap worth less than the same smile's variance
  swap (spot-weighting down-weights the high-vol low-strike puts).

### Fixed
- An initial gamma-swap formula carried a spurious forward "drift" term (copied
  from the variance-swap log contract); the price-weighted log contract has no
  such remainder. Removing it makes the flat case equal `sigma^2` and restores
  the correct gamma < variance ordering under a downward skew.

## [1.147.0] - 2026-09-10

### Added
- `corridor_variance_swap_from_smile` (in `varswap.py`): fair corridor
  variance-swap strike from a smile `vol_fn(K)`. A corridor variance swap
  accrues realized variance only while the spot is in `[lower, upper]`; by
  Carr-Lewis static replication this restricts the `1/K^2`-weighted option strip
  to strikes inside the corridor. Each strip option is priced at its smile vol
  with Black-Scholes.
- Verified: a narrower corridor accrues less variance; nested corridors are
  monotone at fixed strike density; a wide corridor recovers the full
  variance-swap strike (~5e-3); the strike is positive.

## [1.146.0] - 2026-09-10

### Added
- `variance_swap_from_smile` (in `varswap.py`): fair variance-swap strike
  replicated directly from a volatility smile `vol_fn(K)`. Builds the OTM
  option strip (puts below the forward split, calls above) by pricing each
  strike at its smile vol with Black-Scholes and delegates to
  `variance_swap_strike`. Convenient for marking a swap off a fitted smile
  (SVI/SABR/vanna-volga).
- Verified: a flat smile returns its variance `sigma^2` (to ~3e-3, converging as
  the strip widens/densifies -- the residual is finite-strip truncation); a
  downward skew adds convexity variance over the same-grid flat level; and it
  agrees with a manually-built flat-vol chain.

## [1.145.1] - 2026-09-10

### Tests
- Trimmed the fast-suite runtime back from ~45s to ~37s: the smooth `K -> 0`
  Asian-PDE sanity check runs on a modest grid (it needs no strike-kink
  resolution), the ATM/OTM Asian-PDE put check is marked slow, and the
  Heston-calibration parameter-validity check uses a short optimiser run (it
  tests the reparametrization, not fit accuracy). Full accuracy coverage is
  unchanged under `-m slow`.

## [1.145.0] - 2026-09-10

### Added
- `asian_pde_price` (new `pde_asian.py`): continuously-averaged arithmetic Asian
  option by an augmented-state 2D PDE in spot and running integral. The integral
  state is pure transport (`dI = S dt`, no diffusion), so each backward step is
  operator-split: a semi-Lagrangian transport in `I` (interpolate along the
  characteristic) then a Crank-Nicolson diffusion in `S` on each `I`-line. The
  `I`-grid is bounded by the average's reachable range so it resolves the strike
  kink.
- Verified: a `K -> 0` average-price call equals the discounted expected average
  (~0.05); prices agree with Turnbull-Wakeman to ~0.25 and converge toward a
  continuous-monitoring Monte Carlo as the `I`-grid refines.

## [1.144.0] - 2026-09-10

### Added
- `adi_two_asset_american` (in `pde2d.py`): two-asset options with optional
  early exercise on the Peaceman-Rachford ADI grid. When `american=True` the
  value grid is floored at the immediate-exercise payoff after each time step
  (the 2D analogue of the vanilla PSOR floor), for American best-of / worst-of /
  spread payoffs.
- Verified: the European limit (`american=False`) matches the best-of and
  worst-of closed forms and converges in the grid; the American value is at
  least the European; a dividend yield produces a positive early-exercise
  premium; and the American value is floored at intrinsic deep in the money.

## [1.143.0] - 2026-09-10

### Added
- `adi_two_asset_cs` (in `pde2d.py`): the Craig-Sneyd ADI scheme for two-asset
  options. Peaceman-Rachford is only first-order in time when a mixed
  (correlation) derivative is present; Craig-Sneyd restores second order with a
  Douglas predictor followed by a corrector that re-applies the explicit cross
  term at the predicted value. The directional operators are solved implicitly
  (Thomas), the cross term explicitly.
- Verified: matches the Margrabe exchange closed form; at few time steps it is
  markedly more accurate than Peaceman-Rachford (0.0004 vs 0.043 at n_time=10),
  its purpose; it converges in space; and the zero-correlation case matches
  Margrabe.

## [1.142.0] - 2026-09-10

### Added
- Two-asset ADI PDE solver (new `pde2d.py`): `adi_two_asset` prices a European
  two-asset option for any terminal payoff by the Peaceman-Rachford
  alternating-direction-implicit scheme in log-prices -- each time step is two
  half-steps (implicit in x1, then x2) with the mixed correlation derivative
  explicit, each half a set of tridiagonal Thomas solves. `adi_spread_option`
  wraps it for `max(S1 - S2 - K, 0)`.
- Verified: the zero-strike spread matches the Margrabe exchange closed form
  (~2e-3, converging in the grid), a non-zero-strike spread matches the Kirk
  approximation (~3e-3), higher asset correlation lowers the spread price
  monotonically, and a generic max-of-two payoff prices sanely.

## [1.141.0] - 2026-09-10

### Added
- `VannaVolgaSmile.vol_cm` (in `vannavolga.py`): the Castagna-Mercurio analytic
  vanna-volga implied-vol approximation from the three delta pillars. `order=1`
  is the vega/vanna/volga-weighted first-order average
  `sigma_atm + sum_i x_i (sigma_i - sigma_atm)`; `order=2` adds the standard
  second-order convexity correction
  `(-sigma_atm + sqrt(sigma_atm^2 + d1 d2 (2 sigma_atm D1 + D2)))/(d1 d2)`.
- Verified: both orders are exact at the three pillars; the second order matches
  the price-corrected vol off-pillar to ~1e-3; the two orders differ off-pillar;
  and a negative risk reversal gives a downward skew.

## [1.140.0] - 2026-09-10

### Added
- Digital and no-touch binaries via the Crank-Nicolson PDE (`pde.py`):
  `crank_nicolson_digital` prices a cash-or-nothing digital (digital-specific
  boundary conditions, Rannacher damping on by default for the payoff jump), and
  `crank_nicolson_no_touch` prices a no-touch binary as a knock-out of a
  constant cash payoff with an absorbing barrier (a node placed exactly on `H`).
  A pay-at-expiry one-touch is `cash * e^{-rt} - no_touch`. Constant `sigma` or
  `local_vol_fn`, carry `b = r - q`.
- Verified: the digital matches the closed-form `cash_or_nothing` (~4e-3),
  digital call + put equals the discount factor exactly, the no-touch matches
  its closed form (~5e-3), and the one-touch complement matches the analytic
  pay-at-expiry one-touch.

## [1.139.0] - 2026-09-10

### Added
- FX delta-space quoting conventions (new `fxdelta.py`): `strike_from_delta` /
  `delta_from_strike` convert between strike and delta in the forward- or
  spot-delta convention, with an optional premium adjustment; `atm_dns_strike`
  gives the delta-neutral-straddle ATM strike `F exp(0.5 sigma^2 t)`; and
  `rr_bf_to_pillars` turns (ATM, 25d risk reversal, 25d butterfly) quotes into
  the put/ATM/call pillar strikes and vols.
- Verified: strike<->delta round-trips to ~1e-16 (forward, spot and
  premium-adjusted); the ATM DNS strike makes a straddle delta-neutral; the
  pillars reproduce the RR and BF quotes exactly and are strike-ordered; a
  negative risk reversal makes the put vol exceed the call vol.

## [1.138.0] - 2026-09-10

### Added
- `calibrate_ssvi` gains a `vega_weighted` option: each quote's total-variance
  error is weighted by an approximate Black vega
  (`sqrt(w) exp(-d1^2/2)`), so near-the-money quotes -- largest vega, deepest
  liquidity -- dominate the fit rather than the deep wings.
- Verified: a clean surface is still recovered exactly; under wing noise the
  vega-weighted fit has a smaller at-the-money vol error than the unweighted
  fit; the fitted surface stays valid and arbitrage-free; and the default
  (unweighted) behaviour is unchanged.

## [1.137.0] - 2026-09-10

### Added
- `calibrate_svi_from_prices` (in `svi.py`): calibrate a raw SVI slice directly
  from market call *prices* rather than pre-inverted vols. Inverts each call to
  its Black-Scholes implied vol, converts to total variance, and fits raw SVI,
  vega-weighting the quotes by default (so a price error in the deep wings,
  where vega is tiny, does not dominate the fit). Returns
  `(params, iv_rmse, price_rmse)`.
- Verified: recovers a synthetic SVI slice from clean prices exactly (iv/price
  RMSE ~0); fits noisy prices; vega-weighted and unweighted both fit clean
  prices; and the fitted slice reprices each strike to within the reported price
  RMSE.

## [1.136.0] - 2026-09-10

### Added
- Rannacher time-stepping in the Crank-Nicolson solver (`pde.py`): the first
  `rannacher` steps out of expiry are taken fully implicit (backward Euler)
  before switching to Crank-Nicolson, damping the spurious oscillation the
  non-smooth payoff kink induces in CN (which otherwise corrupts gamma/theta
  near the strike). Exposed as a `rannacher` argument on `crank_nicolson_price`
  and `crank_nicolson_greeks` (default 2; `0` recovers pure CN).
- The theta-scheme was generalised to a per-step weight so the same code runs
  backward Euler and Crank-Nicolson. Verified: prices still match Black-Scholes;
  `rannacher=0` recovers pure CN; and on a coarse time grid the Rannacher gamma
  near the strike is no worse than (here better than) pure CN.

## [1.135.0] - 2026-09-10

### Added
- `crank_nicolson_barrier` (in `pde.py`): continuously-monitored single-barrier
  options by a Crank-Nicolson PDE with an absorbing boundary -- the value is
  pinned to the rebate on the dead side of the barrier at every time step. The
  grid is anchored so a node lands exactly on `H` (otherwise the barrier is
  applied at an offset node). `down-out`/`up-out` are solved directly;
  `down-in`/`up-in` come from in + out = vanilla. Supports a constant `sigma` or
  a `local_vol_fn` and a carry `b = r - q`.
- Verified: matches the Reiner-Rubinstein closed form on all four barrier kinds
  and a put (within the barrier's O(ds) discretisation bias on a fine grid);
  knock-in + knock-out reproduces the vanilla; a knock-out is below the vanilla;
  spot at the barrier returns the rebate.

## [1.134.0] - 2026-09-10

### Added
- `crank_nicolson_greeks` (in `pde.py`): delta, gamma and theta read straight
  off the Crank-Nicolson grid -- no extra solves. Delta and gamma are node-level
  central differences interpolated to the spot (so the grid offset does not bias
  them), and theta is the difference between the `t=0` grid and the grid one
  time step earlier. The solver was refactored into a shared `_cn_solve` used by
  both the price and Greeks entry points.
- Verified against Black-Scholes: delta to ~1e-4, gamma to ~5e-6, theta to
  ~3e-3; the put delta matches; American-put Greeks have the right signs.

## [1.133.0] - 2026-09-10

### Added
- Crank-Nicolson PDE solver (new `pde.py`): `crank_nicolson_price` solves the
  Black-Scholes-Merton PDE on a spot grid with the (unconditionally stable,
  second-order) Crank-Nicolson scheme. Supports a constant `sigma` or a
  `local_vol_fn(S, t)`, a cost-of-carry `b = r - q`, and American exercise via
  projected SOR (PSOR) enforcing `V >= payoff`. The interior tridiagonal system
  is solved by the Thomas algorithm.
- Verified: European calls and puts match Black-Scholes; the error converges
  `O(1/n^2)` in the grid; the dividend-carry and flat-local-vol cases match BS;
  the American put matches a 3000-step binomial tree to ~1e-2 and exceeds the
  European value.

## [1.132.0] - 2026-09-10

### Added
- Arbitrage-free SSVI calibration (in `ssvi.py`): `calibrate_ssvi` gains an
  `arb_weight` that penalises butterfly and calendar no-arbitrage violations,
  and `calibrate_ssvi_arbitrage_free` ramps that weight (warm-starting each fit
  from the last) until the surface passes `ssvi_is_arbitrage_free`. Trades a
  little fit RMSE for a guaranteed no-arbitrage surface.
- Verified: an intentionally arbitraging SSVI surface is flagged and the plain
  fit reproduces it, while the penalised fit returns an arbitrage-free surface
  (pulling `eta` 8.0 -> 1.8); on a clean surface the arb-free fit matches the
  plain fit's RMSE and removing arbitrage never improves the fit to an
  arbitraging market.

## [1.131.0] - 2026-09-10

### Added
- LSV leverage-function calibration (new `lsv.py`): `calibrate_lsv_leverage`
  calibrates a local-stochastic-volatility leverage surface `L(S, t)` by the
  Guyon-Henry-Labordere particle method -- Gyongy's condition
  `L(K,t)^2 = sigma_Dupire(K,t)^2 / E[V_t | S_t = K]`, with the conditional
  expectation estimated by simulating the Heston variance and binning `V_t` by
  the spot level, marched forward with Euler sub-steps between the calibration
  expiries. Returns the leverage grid and an interpolating `lev_fn(spot, t)`.
- Verified: the front-expiry ATM leverage equals `sigma_local / sqrt(v0)`
  exactly (Gyongy at `t -> 0`, `E[V|S] = v0`); leverage is positive everywhere;
  a higher target vol scales it up; and the flat-target front-expiry leverage is
  uniform at `sigma_local / sqrt(v0)`.

## [1.130.0] - 2026-09-10

### Added
- `calibrate_heston` (new `heston_calib.py`): fit the five Heston parameters
  `(v0, kappa, theta, xi, rho)` to an implied-vol surface by least squares on
  Black vol over `(expiry, strike, vol)` quotes, pricing each candidate with the
  Fourier `heston_price` and optimising with Nelder-Mead under a smooth
  reparametrization. An optional `feller_weight` penalises Feller-condition
  violations (`2 kappa theta < xi^2`) to prefer a strictly-positive variance
  process.
- Verified: recovers a synthetic Heston surface's parameters (rmse ~0); the
  Feller penalty pushes `2 kappa theta` up to the `xi^2` boundary; all fitted
  parameters stay in their valid regions.

## [1.129.0] - 2026-09-10

### Added
- `calibrate_double_heston` (new `double_heston_calib.py`): fit all ten
  double-Heston parameters (two variance factors) to an implied-vol surface by
  least squares on Black vol over `(expiry, strike, vol)` quotes, pricing each
  candidate with the Fourier `double_heston_price` and optimising with
  Nelder-Mead under a smooth reparametrization (variances/vol-of-vols positive,
  correlations in `(-1, 1)`). A two-scale (fast + slow factor) seed is used by
  default.
- Verified: calibrating to a synthetic double-Heston surface (4 expiries x 5
  strikes) fits to ~2e-5 RMSE with all fitted parameters in their valid regions.

## [1.128.1] - 2026-09-10

### Fixed
- `rough_heston_price` now raises a clear error when the fractional Riccati
  diverges (a too-coarse `n_grid` for small `H`) instead of silently returning
  NaN -- found while trimming the test suite, where an H=0.2 parity check on a
  60-step grid produced NaN.

### Tests
- Trimmed the fast-suite runtime from ~66s to ~28s: marked the heaviest
  Fourier/calibration cross-checks `slow` (rough-Heston Hurst-reduction and
  short-skew, Levy NIG/Meixner/cross-model recovery, the Levy-surface calendar
  scan, the 256-step rBergomi FFT identity) and cut redundant path counts on the
  FFT-identity checks (numerical-identity, not statistical). Full coverage is
  unchanged; the heavy checks still run under `-m slow`.

## [1.128.0] - 2026-09-10

### Added
- SABR butterfly-arbitrage detection and repair (in `sabr.py`): `sabr_density`
  is the Breeden-Litzenberger implied density of a SABR smile (second strike
  derivative of the Black call at the SABR vol); `sabr_butterfly_arbitrage`
  returns the strikes where that density goes negative (the Hagan expansion is
  not guaranteed arb-free in the wings); `sabr_is_arbitrage_free` is the
  boolean; and `sabr_repair_butterfly` shrinks the vol-of-vol `nu` until the
  density is non-negative, preserving alpha/beta/rho.
- Verified: a benign smile's density integrates to 1 and is arbitrage-free; an
  extreme vol-of-vol flags negative-density strikes; repair restores arbitrage
  freedom (`nu` 2.5 -> 1.66) while leaving a clean smile untouched.



### Added
- CMS convexity adjustment (new `cms.py`): `cms_adjustment_standard` gives the
  linear-TSR / Hagan standard-model adjustment `G * Var_A(S_T)` (using the exact
  lognormal variance), `cms_rate_convexity_replication` computes it by model-free
  static replication over a swaption strip with a caller-supplied volatility
  smile (Carr-Madan second-moment replication times the level factor `G`), and
  `cms_rate` returns the convexity-adjusted expected CMS rate.
- Verified: zero vol gives no adjustment; the adjustment is positive; the
  flat-smile replication matches the standard model to ~1e-5; both match a
  linear-TSR Monte Carlo to ~2e-4; a convex smile raises the adjustment; and a
  payment lag increases it.



### Added
- Swaption volatility cube (new `volcube.py`): `VolCube` stores a calibrated
  SABR smile at each `(expiry, tenor)` node and interpolates across all three
  market axes -- SABR handles strike/smile analytically, while expiry and tenor
  use bilinear interpolation of the total variance `sigma^2 * expiry` (the
  no-calendar-arbitrage-friendly choice). `VolCube.fit` builds it from market
  smiles; `vol(expiry, tenor, strike)` queries any point.
- Verified: each node reprices its own SABR smile, an interpolated `(3, 7)` node
  lands between its neighbours, expiry interpolation stays in range, the SABR
  smile shape is preserved, and beyond the grid the nearest node is used.

## [1.125.0] - 2026-09-10

### Added
- OIS/LIBOR dual-curve discounting (new `dualcurve.py`): value swaps under the
  post-2008 convention where forward rates are projected off a forward (LIBOR)
  curve while every cashflow is discounted on the collateral (OIS) curve.
  `dual_forward_rate`, `dual_float_leg_value`, `dual_par_swap_rate` and
  `dual_swap_value` take separate OIS and projection `DiscountCurve`s plus an
  additive basis; `dual_calibrate_basis` solves the constant basis spread that
  reprices a set of par swaps.
- Verified: with a single shared curve the dual-curve par rate collapses to the
  single-curve one exactly; the flat-curve forward equals `e^r - 1`; a known
  basis is recovered to ~1e-6; a payer swap is worth zero at par and
  payer + receiver = 0; adding basis raises the float leg; and a projection
  curve above OIS raises the par rate.

## [1.124.0] - 2026-09-10

### Added
- Discount curve (new `discount_curve.py`): `DiscountCurve` builds a log-linear
  (piecewise-constant instantaneous-forward) discount curve from pillar
  `(T, DF)` points or continuously-compounded zero rates, with zero-rate,
  forward-rate and par-swap-rate accessors. `bootstrap_from_swaps` bootstraps
  pillar discount factors from par swap rates (each pillar solves a linear
  equation given the shorter ones). This is exactly the initial-curve object the
  Gaussian short-rate models (`cheyette`, `g2pp`) consume.
- Verified: a flat zero-rate curve reproduces `e^{-zT}`; interpolation is exact
  at pillars; a consecutive-annual-tenor bootstrap reproduces its par rates to
  1e-10; and G2++ reprices the bootstrapped curve's own bonds exactly.

## [1.123.0] - 2026-09-10

### Added
- Rough-Heston model (new `rough_heston.py`): `rough_heston_price` and
  `rough_heston_smile`. The variance is driven by a fractional kernel with Hurst
  `H in (0, 0.5]`; its characteristic function comes from the fractional Riccati
  equation, solved here by the fractional Adams predictor-corrector
  (Diethelm-Ford-Freed) and inverted by the same Gil-Pelaez two-probability
  Gauss-Legendre integral as Heston.
- Verified: `H = 0.5` reduces to classical Heston (with `xi = kappa * nu`, the
  El Euch-Rosenbaum convention) to <5e-3 and converges in the grid size; the
  characteristic function is a martingale; put-call parity holds; and `H = 0.1`
  gives a much steeper short-dated skew (-0.88) than `H = 0.5` (-0.42) -- the
  rough-vol signature.

## [1.122.0] - 2026-09-10

### Added
- Double-Heston two-factor stochastic-volatility model (new `double_heston.py`):
  `double_heston_price` and `double_heston_smile`. Two independent Heston
  variance factors (a fast- and a slow-reverting one) let the short- and
  long-dated skew move more independently than single-factor Heston allows.
  Since the factors are independent, the log-spot characteristic function is the
  product of the two single-factor pieces; priced by the same Gil-Pelaez
  two-probability Gauss-Legendre integral as Heston. The `xi -> 0` factor limit
  is handled analytically (deterministic-variance contribution).
- Verified: zeroing the second factor recovers single-factor Heston exactly,
  put-call parity holds, two factors are worth more than one, negative
  correlations give a downward skew, and the Fourier price matches an
  independent two-factor QE Monte Carlo to ~0.6 SE.

## [1.121.0] - 2026-09-10

### Added
- `bermudan_swaption_g2pp` (new `bermudan_swaption.py`): Bermudan swaption
  pricing by Longstaff-Schwartz on the G2++ state. Simulates the correlated
  ``(x, y)`` factors on the exercise schedule (exact OU transitions), carries a
  discretely-compounded money-market numeraire from the one-period G2++ bonds,
  values the co-terminal swap in closed form from the state, and regresses the
  discounted continuation on a quadratic basis in ``(x, y)`` for the
  exercise decision.
- Verified: price is positive; a multi-date Bermudan is worth at least the
  single-exercise value; payer and receiver are both positive; a payer's value
  falls as the fixed rate rises; payer/receiver swap values are antisymmetric;
  the simulated factors have ~zero mean.

## [1.120.0] - 2026-09-10

### Added
- Two-factor G2++ Gaussian short-rate model (new `g2pp.py`): two correlated
  mean-reverting factors, `r = x + y + phi(t)`. `g2pp_zero_bond` is the
  exponential-affine bond off the initial curve (with the closed-form variance
  `g2pp_V`), `g2pp_bond_option` the exact Gaussian bond-option price, and
  `g2pp_caplet` via the bond-put identity.
- Verified: the bond option matches a forward-measure Monte Carlo to ~1e-5,
  put-call parity holds, turning the second factor off (`eta -> 0`, `rho = 0`)
  reproduces the single-factor Cheyette price exactly, factor correlation raises
  the bond vol and option price monotonically, and `sigma = eta = 0` gives the
  discounted intrinsic.

## [1.119.0] - 2026-09-10

### Added
- Single-factor Cheyette / quasi-Gaussian short-rate model (new `cheyette.py`):
  the Markovian HJM representation with state `(x, y)` and short rate
  `r = f(0,t) + x`. `cheyette_zero_bond` reconstitutes `P(t,T)` off the initial
  curve via the `G(t,T)` function, and `cheyette_bond_option` / `cheyette_caplet`
  give the exact prices for a constant short-rate vol (where the factor
  coincides with Hull-White), the caplet using the bond-put identity.
- Verified: the bond option matches a forward-measure Monte Carlo to ~1e-5,
  put-call parity holds, the caplet equals its `(1 + K tau)` bond-put
  composition exactly, `kappa -> 0` stays finite, and `sigma -> 0` collapses to
  the discounted intrinsic.

## [1.118.0] - 2026-09-10

### Added
- Andreasen-Huge single-step arbitrage-free local-vol smile (new
  `andreasenhuge.py`): `andreasen_huge_prices` solves one implicit Dupire
  finite-difference step (a tridiagonal M-matrix system, Thomas algorithm) to
  produce call prices that are monotone-decreasing and convex in strike -- hence
  arbitrage-free -- for any positive local-vol grid. `andreasen_huge_smile`
  inverts them to implied vols, and `andreasen_huge_calibrate` bootstraps the
  per-strike local vols that reprice a market smile.
- Verified: the prices are monotone and convex for both flat and steeply skewed
  local vols; calibration reproduces a market skew to <5e-3 RMSE (interior
  strikes to ~1e-4); the calibrated smile is downward-sloping.

## [1.117.0] - 2026-09-10

### Added
- `kim_put_greeks` (in `kim.py`): delta, gamma and theta of a Kim American put.
  The early-exercise boundary is spot-independent, so it is solved once (the
  expensive step) and the spot bumps only re-run the cheap European-plus-premium
  evaluation on that fixed boundary; theta uses a maturity bump.
- Verified against binomial bump Greeks: delta matches to ~5e-3, gamma to ~2e-3
  (stable in the grid size -- the small residual is Kim's boundary
  approximation, not quadrature noise), put delta is negative, gamma positive,
  and theta negative.

## [1.116.0] - 2026-09-10

### Added
- `baw_american` (new `baw.py`): the Barone-Adesi-Whaley (1987) quadratic
  approximation for American options. Splits the price into the European value
  plus an early-exercise premium (a power of spot anchored at the critical price
  ``S*``, found by a 1-D Newton solve of the value-matching condition). Fast
  closed-form, dividends via ``b = r - q``; a no-dividend call returns the
  European value.
- Verified: within ~0.1 of a 4000-step binomial tree for puts and dividend
  calls (BAW's known approximation accuracy), within ~0.05 of the accurate Kim
  integral-equation price, and correctly floored at intrinsic.

## [1.115.0] - 2026-09-10

### Added
- American pricing via Kim's (1990) integral equation (new `kim.py`):
  `kim_american_put`, `kim_american_call`, `kim_exercise_boundary`. The
  early-exercise boundary solves a Volterra integral equation by backward
  marching from expiry, bisecting the value-matching condition at each step; the
  price is the European value plus the early-exercise premium integrated at
  spot. The call uses the McDonald-Schroder put-call symmetry
  `C(S,K,r,q) = P(K,S,q,r)`.
- Verified: puts and dividend calls match a 4000-step binomial tree to ~1e-3 at
  200 grid steps; the boundary rises to K at expiry; a no-dividend call equals
  the European value.

### Fixed
- Two bugs caught against the binomial reference: the premium integral indexed
  the boundary backward (`B[n-m]` instead of `B[m]`, the boundary `s` years
  ahead), and the `s -> 0` integrand endpoint was forced to zero when on the
  boundary it is `rK/2 - qS/2` (from `N(0)=1/2`). With both fixed the price
  converges to the tree instead of ~1.8% low.

## [1.114.0] - 2026-09-10

### Added
- Sobol sequence + Brownian-bridge QMC (new `sobol.py`): a `Sobol` generator
  built from primitive-polynomial direction numbers (Joe-Kuo initial values) via
  the Gray-code recurrence, `brownian_bridge_path` which loads a path's dominant
  variance onto the leading (most uniform) Sobol coordinates, and two pricers --
  `sobol_european` (1-D) and `sobol_asian` (Brownian-bridge, multi-step).
- Verified: Sobol points cover the unit cube evenly (each quarter of the 1-D
  coordinate within 2%), `sobol_european` matches Black-Scholes and its error is
  materially smaller than pseudo-random MC at matched N (~0.007 vs ~0.09 at
  N=8192), the Brownian bridge has the right terminal moments (mean ~0, variance
  ~t), and `sobol_asian` matches a discrete-monitoring pseudo-random MC to ~0.01.

## [1.113.1] - 2026-09-10

### Changed
- The rBergomi conditional (turbocharged) estimator's `I1`/`QV` path loop now
  shares the same cached-kernel FFT Volterra convolution as `rbergomi_paths`
  (factored into module-level `_build_far_kernel` / `_far_sums_fft`), with the
  same `fast="auto"` switch at `n_steps >= 200`. Verified the FFT and direct
  `I1`/`QV` agree to ~1e-14 and the CV price/tests are unchanged.

## [1.113.0] - 2026-09-10

### Added
- `LevySurface` (new `levysurface.py`): a unified implied-vol surface generated
  by a single exponential-Levy parameter set (`vg`, `nig`, `meixner`, `cgmy`).
  `implied_vol(k, t)` / `total_variance(k, t)` evaluate any point, `grid(...)`
  builds an expiry x log-moneyness surface with one Carr-Madan FFT strip per
  maturity, and `calendar_violations` / `is_calendar_arbitrage_free` check that
  total variance is non-decreasing in maturity (the same test as the SVI
  `VolSurface`).
- `levy_psi(model, params)` factory (in `levycalib.py`): returns the
  characteristic exponent for a named model + raw parameters, for pricing or
  surface-building outside calibration.
- Verified: the surface reprices the NIG and VG per-model smiles to ~1e-3, a
  genuine Levy law is calendar-arbitrage-free across five expiries, total
  variance grows with maturity, and `beta < 0` gives a downward skew at each
  expiry.

## [1.112.0] - 2026-09-10

### Changed
- rBergomi path generation (`rbergomi_paths`) now evaluates the hybrid-scheme
  Volterra convolution with a cached-kernel FFT instead of the `O(n^2)` inner
  double loop. The far-cell contribution `F[i] = sum_j dW[j] g[i-j]` is a
  discrete convolution of the fixed weight kernel with the Brownian increments;
  the kernel's FFT is precomputed once and reused across every path (one forward
  + one inverse transform of the increments per path), giving `O(n log n)`. A
  new `fast` argument (`"auto"` default) uses the FFT when `n_steps >= 200` and
  the direct loop below that (where pure-Python FFT overhead dominates); the two
  agree to ~1e-12. Measured ~2.3x faster at 512 steps, growing with `n_steps`.

## [1.111.0] - 2026-09-10

### Added
- `calibrate_levy_smile` (new `levycalib.py`): fit an exponential-Levy model
  (`vg`, `nig`, `meixner`, `cgmy`) to a one-expiry market implied-vol smile by
  least squares on vol. Each candidate smile is priced with a single Carr-Madan
  FFT strip (interpolated to the market strikes) rather than an integral per
  strike, and a per-model smooth parameter transform keeps the optimiser
  unconstrained while the raw parameters stay in their valid region.
- Verified: recovers synthetic VG, NIG and Meixner parameters to <1e-3 RMSE;
  CGMY (weakly identified from few strikes) is checked on fit quality (RMSE
  < 5e-3); a cross-model NIG fit to a VG smile is close (RMSE ~1e-3).

## [1.110.0] - 2026-09-10

### Added
- `variance_gamma_smile` (in `variancegamma.py`): the Black-Scholes implied-vol
  smile the VG model produces, `theta < 0` giving a downward skew.

### Changed
- `variance_gamma_price` refactored onto the shared `carrmadan.levy_price`
  engine via the closed-form VG exponent
  `psi(u) = -(1/nu) log(1 - i theta nu u + 0.5 sigma^2 nu u^2)`, replacing its
  own two-probability Gauss-Legendre loop. Same prices, and it now inherits the
  COS-method cross-check like the other Levy models. Verified: Carr-Madan and
  COS agree to ~4.5e-6 on three parameter sets, the existing VG tests still
  pass, parity holds, and small `nu` recovers Black-Scholes.

## [1.109.0] - 2026-09-10

### Added
- `cos_greeks` (in `carrmadan.py`): analytic delta and gamma of the COS-method
  price for any Levy model. Spot enters the cosine sum only through the
  characteristic function's `e^{i u (x + mu)}` factor with `x = ln(S/K)`, so
  differentiating term-by-term gives delta (`i u / S` per term) and gamma
  (`i u (i u - 1)/S^2`) with no re-pricing and no finite differences; the put
  follows by parity.
- Verified: matches Black-Scholes delta/gamma to 1e-6 for the GBM exponent, and
  matches finite-difference COS Greeks on CGMY and NIG to ~1e-4.

## [1.108.0] - 2026-09-10

### Added
- Meixner Levy model (new `meixner.py`): `meixner_price` and `meixner_smile`.
  Schoutens' Meixner process has the analytic exponent
  `psi(u) = 2 d (ln cos(b/2) - ln cosh((a u - i b)/2))` (scale `a`, asymmetry
  `b in (-pi, pi)`, activity `d`) and is priced through the shared Carr-Madan
  engine, so it also gets the COS cross-check for free.
- Verified: Carr-Madan and COS prices agree to ~1.7e-6 on three parameter sets,
  put-call parity holds, `b < 0` gives a downward skew and `b > 0` an upward
  skew, and more activity `d` raises the price.

## [1.107.0] - 2026-09-10

### Added
- COS method (in `carrmadan.py`): `cos_price` implements the Fang-Oosterlee
  (2008) Fourier-cosine expansion of the risk-neutral log-return density on a
  cumulant-based truncation range, with closed-form call payoff coefficients
  (chi/psi). Exponentially convergent and independent of the Carr-Madan
  transform, so it is a genuine cross-check for the whole Levy family.
- Verified: matches Black-Scholes to 1e-8 (GBM exponent), and matches the
  Carr-Madan pricer to ~1e-6/1e-4 on CGMY and three NIG parameter sets; parity
  holds and the error shrinks with the term count.

### Fixed
- The COS characteristic function omitted the `x = ln(S/K)` shift, so only ATM
  strikes priced correctly; and the fourth-difference c4 cumulant estimate blew
  up on non-smooth exponents (CGMY's Gamma(-Y) power law), ballooning the
  truncation range until the density undersampled. c4 now uses a larger step and
  is clamped to a sane multiple of c2^2 (it only fine-tunes the range).

## [1.106.0] - 2026-09-10

### Added
- Carr-Madan FFT strip (in `carrmadan.py`): `carr_madan_strip` prices a whole
  log-strike grid of European calls in a single transform, and
  `carr_madan_smile_strip` returns the implied-vol smile over a log-moneyness
  window from that one pass. Backed by a pure-Python radix-2 Cooley-Tukey `_fft`
  (bit-reversal + butterflies, no NumPy) and Simpson-weighted frequency
  sampling; works for any Levy model via its characteristic exponent (CGMY, NIG,
  ...). The log-strike spacing is `2 pi / (n_fft * eta)`.
- Verified: `_fft` matches a naive DFT and round-trips to ~1e-12; on the native
  FFT grid the CGMY and NIG strip prices match the per-strike Gauss-Legendre
  pricer to ~2e-6, and the smile strip is sorted, in-window, and shows the
  expected skew.

## [1.105.0] - 2026-09-10

### Added
- Shared Carr-Madan engine (new `carrmadan.py`): `levy_price` / `carr_madan_call`
  price any exponential-Levy model from its characteristic exponent `psi(u)`,
  applying the martingale correction `omega = -psi(-i)` and the alpha-damped
  Fourier inversion over the shared Gauss-Legendre nodes.
- Normal Inverse Gaussian model (new `nig.py`): `nig_price` and `nig_smile`.
  Barndorff-Nielsen's NIG has the analytic exponent
  `psi(u) = delta (sqrt(alpha^2 - beta^2) - sqrt(alpha^2 - (beta + i u)^2))`
  (tail `alpha`, asymmetry `beta`, scale `delta`) and is priced through the
  shared engine.
- Verified NIG against an independent Gil-Pelaez inversion to ~1e-6 on three
  parameter sets, with put-call parity, `beta < 0` -> downward skew and
  `beta > 0` -> upward skew.

### Changed
- `cgmy_price` now routes through the shared `carrmadan.levy_price` instead of
  its own inlined Fourier loop (identical prices; the `_cgmy_char_logspot`
  helper is kept for the Gil-Pelaez cross-check).

## [1.104.0] - 2026-09-10

### Added
- CGMY tempered-stable Levy model (new `cgmy.py`): `cgmy_price` and `cgmy_smile`.
  A pure-jump process with a stable Levy density tempered independently on each
  tail (`C` activity, `G`/`M` down/up tempering, fine-structure `Y < 2`), whose
  closed-form characteristic exponent uses `Gamma(-Y)` and complex powers.
  Priced by Carr-Madan Fourier inversion of the damped call transform (damping
  `alpha`, requiring `alpha + 1 < M`) over the shared Gauss-Legendre nodes, with
  a martingale drift correction `omega = -psi(-i)`.
- Verified against an independent Gil-Pelaez inversion of the same
  characteristic function to ~1e-6 across three regimes (`Y = 0.5, 0.8, 1.2`);
  put-call parity holds, `G < M` produces a downward skew and `G = M` a
  symmetric smile, and more activity (`C`) raises the price.

## [1.103.0] - 2026-09-10

### Added
- `leisen_reimer_american_accel` (in `leisen_reimer.py`): the Broadie-Detemple
  (1996) two-point Richardson extrapolation of the American Leisen-Reimer price,
  `V_ext = 2*V(2n) - V(n)`, which cancels the leading `1/n` error that the
  American tree carries (its exercise boundary breaks the smooth-payoff
  assumption behind Peizer-Pratt). Measured against a 6000-step CRR tree it is
  1.6-20x more accurate than a single 101-step LR tree for the same work,
  reaching a few mils; `b = r - q` supports dividends.

## [1.102.0] - 2026-09-10

### Added
- Leisen-Reimer binomial tree (new `leisen_reimer.py`): `leisen_reimer_price`
  and `leisen_reimer_greeks`. Centres the tree on the strike and sets the
  up-move and probability from a Peizer-Pratt inversion of the Black-Scholes
  d1/d2, giving smooth `O(1/n^2)` convergence instead of CRR's slow oscillation.
  Handles European and American exercise and a continuous dividend yield via
  `b = r - q`; the step count is forced odd so the strike sits at the tree
  centre.
- Verified: European prices match Black-Scholes to <2e-3 at ~51 steps and beat
  CRR at the same step count; American prices with dividends reach a 3000-step
  CRR tree to ~1 cent (the American convergence is slower than the European
  `O(1/n^2)`, since the smooth-payoff assumption breaks at the exercise
  boundary, so a larger step count is used).

## [1.101.0] - 2026-09-10

### Added
- Kou (2002) double-exponential jump-diffusion (new `kou.py`): `kou_price` and
  `kou_smile`. The jump size is an asymmetric double exponential (up-jump tail
  rate `eta1 > 1`, down-jump tail rate `eta2 > 0`, up-probability `p`), giving
  fatter, asymmetric tails than Merton's Gaussian jumps with a still-analytic
  Levy characteristic function. Priced with the same two-probability
  Gauss-Legendre Fourier integral as the Heston/Bates pricers, with the jump
  drift compensator keeping the discounted spot a martingale.
- Verified: `lambda = 0` recovers Black-Scholes to 1e-14, put-call parity holds,
  jumps raise the price, asymmetric jumps produce a downward skew, and the
  Fourier price matches an independent double-exponential-jump Monte Carlo to
  <0.9 standard errors on three parameter sets.

## [1.100.0] - 2026-09-10

### Added
- `svi_local_variance` and `svi_surface_local_vol` (in `svi.py`): analytic
  Dupire local volatility for raw SVI, mirroring the SSVI local-vol pair. The
  strike derivatives w_k, w_kk come in closed form from the SVI parametrization
  (reusing `_svi_derivs`); `svi_local_variance` takes a caller-supplied dw/dt
  (one slice carries no maturity), while `svi_surface_local_vol` interpolates
  total variance linearly across a term structure of slices to supply dw/dt.
  Both raise on a non-positive Dupire denominator (butterfly-arbitrage flag).
- Verified: the surface local vol matches a finite-difference Dupire on the same
  slice term structure to ~1e-7 across strikes and maturities.

## [1.99.0] - 2026-09-10

### Added
- Bates (1996) model (new `bates.py`): Heston stochastic volatility plus Merton
  lognormal jumps. `bates_price` reuses the Heston two-probability Fourier
  integral and Gauss-Legendre machinery, multiplying the Heston characteristic
  function by the independent compound-Poisson jump factor with a martingale
  drift compensator. `bates_smile` returns the implied-vol smile. `lambda = 0`
  recovers the Heston price exactly.
- Verified: `lambda = 0` matches Heston to 1e-8, put-call parity holds, jumps
  raise the price vs Heston, a negative mean jump steepens the downward skew,
  and -- the decisive check -- the Fourier price matches an independent
  Andersen-QE variance simulation with a compound-Poisson jump overlay to <0.35
  standard errors on three parameter sets.

### Fixed
- The jump characteristic factor initially used a `+-1/2` measure shift and an
  `i*phi` compensator; the Monte Carlo cross-check exposed a 12-86 sigma error.
  The correct Heston-decomposition shift is `s = 1` (share measure, argument
  `phi - i`) for P1 and `s = 0` for P2, with the compensator `-a*lambda*t*k`;
  after the fix the MC agreement is <0.35 sigma. Parity and the `lambda = 0`
  limit alone did not catch this -- only the independent jump-aware MC did.

## [1.98.1] - 2026-09-10

### Tests
- Cross-check the Turnbull-Wakeman arithmetic-Asian Greeks against an
  independent method: bump Greeks of the arithmetic-Asian Monte Carlo under
  common random numbers. Delta, gamma and vega agree to ~1-2% (TW is a
  two-moment approximation), and the arithmetic and geometric Greeks converge at
  low vol. The prior Asian-Greeks tests only finite-differenced the same closed
  form; these confirm the analytic values against a separate pricer.

## [1.98.0] - 2026-09-10

### Added
- `ssvi_local_vol_fn` and `ssvi_reprice_mc` (in `ssvi.py`): close the calibrate
  -> local-vol -> reprice loop. `ssvi_local_vol_fn` turns a fitted SSVI surface
  into a `(spot, tau) -> local vol` callable (mapping spot to SSVI log-moneyness
  on the forward), and `ssvi_reprice_mc` simulates that surface through
  `local_vol_mc`. A correct local-vol construction reprices the SSVI *implied*
  smile, and it does: the Monte Carlo price matches the closed-form price at the
  SSVI implied vol across the smile to within ~1.2 standard errors, provided the
  fitted surface includes short maturities.

### Fixed
- `ssvi_local_vol_from_params` / `ssvi_local_vol_fn` froze the ATM variance below
  the shortest fitted expiry, which biased the local-vol Monte Carlo more as the
  step count grew (the path integrates `tau` from 0 but saw a frozen short-end
  vol). The short end now linearly extrapolates `theta -> 0` toward the origin,
  matching SSVI's small-time behaviour; the bias-vs-steps drift is gone.

## [1.97.0] - 2026-09-10

### Added
- `ssvi_local_variance` and `ssvi_local_vol_from_params` (in `ssvi.py`): the
  Dupire local variance of an SSVI surface in fully analytic form. The total
  variance's strike derivatives (w_k, w_kk) and its theta-derivative come from
  the SSVI parametrization in closed form -- no finite differences -- and feed
  Gatheral's total-variance Dupire formula. `ssvi_local_vol_from_params` builds
  theta(t) and theta'(t) by piecewise-linear interpolation of the fitted ATM
  variances and returns the local vol at any (k, t) in range; it also raises on
  a non-positive Dupire denominator (a butterfly-arbitrage flag).
- Verified: the analytic local variance matches a finite-difference Dupire on
  the same surface to ~1e-9 across strikes and maturities.

## [1.96.0] - 2026-09-10

### Added
- Surface SVI (new `ssvi.py`): the Gatheral-Jacquier (2014) arbitrage-free
  whole-surface parametrization, tying every expiry's smile together through a
  shared skew function and the ATM total-variance term structure `theta_t`.
  `calibrate_ssvi` fits the global `(rho, eta, gamma)` power-law skew plus one
  `theta` per expiry to a `(t, k, iv)` market surface; `ssvi_butterfly_free`,
  `ssvi_calendar_free` and `ssvi_is_arbitrage_free` implement the sufficient
  no-static-arbitrage conditions (density positivity within a slice, total
  variance non-decreasing in maturity across slices).
- Verified: recovers a synthetic surface's parameters to ~1e-8 RMSE, the
  butterfly condition flags an excessive-skew slice, the calendar condition
  flags a decreasing-`theta` term structure, and `phi` decays with maturity.

## [1.95.0] - 2026-09-10

### Added
- `rbergomi_price_cv` and `rbergomi_smile_cv` (in `rbergomi.py`): the conditional
  ("turbocharged") rough Bergomi estimator of McCrickerd & Pakkanen (2018).
  Conditioning on the volatility-driving Brownian motion makes the terminal
  log-spot Gaussian, so each path contributes a smooth Black-Scholes conditional
  price (effective spot `S exp(rho I1 - rho^2 QV/2)`, effective variance
  `(1 - rho^2) QV`) instead of a noisy indicator payoff. Same price as the plain
  estimator within Monte Carlo error, with the standard error cut severalfold --
  ~2x at rho=-0.7 up to ~6-7x as `|rho|` shrinks (the integrated-out noise
  fraction is `1 - rho^2`).

## [1.94.0] - 2026-09-10

### Added
- `rbergomi_price` and `rbergomi_smile` (new `rbergomi.py`): the rough Bergomi
  stochastic-volatility model (Bayer-Friz-Gatheral 2016), whose variance is
  driven by a rough fractional process with Hurst exponent `H < 1/2`. Simulated
  with the Bennedsen-Lunde-Pakkanen (2017) hybrid scheme (kappa=1): the singular
  cell nearest each step is sampled exactly from the joint law of the Brownian
  increment and its kernel integral, the far cells use the optimally-placed
  discretised kernel. `rbergomi_smile` reprices every strike on common terminal
  spots and inverts to a Black-Scholes vol.
- Verified: the discounted spot is a martingale (`E[S_T] e^{-rt} = S` to <5e-3),
  `eta = 0` collapses to a flat smile at `sqrt(xi0)`, `rho < 0` gives a downward
  skew, and -- the point of rough vol -- `H = 0.1` yields a steeper short-dated
  ATM skew (-0.58) than the `H = 0.5` diffusive case (-0.41).

### Fixed
- Hybrid-scheme weights divided by `alpha = H - 1/2`, blowing up at `H = 1/2`;
  that case now returns the constant-kernel (standard Brownian) weights.

## [1.93.0] - 2026-09-10

### Added
- `VannaVolgaSmile.price` and `.vol_price_corrected` (in `vannavolga.py`): the
  exact second-order Castagna-Mercurio vanna-volga construction. Starts from the
  flat-ATM Black-Scholes price and adds the three pillar options' market-minus-
  ATM price gaps, weighted (via a 3x3 vega/vanna/volga solve) so the hedging
  portfolio matches the target's vega, vanna and volga. Unlike the existing
  quadratic vol interpolation (`.vol`), this reprices the three market
  instruments *exactly in price*, the standard FX smile pricer.
- Verified: reprices all three market pillars to <1e-10, price-corrected vol
  equals the pillar vols, put-call parity holds to 1e-16, a negative risk
  reversal puts the put wing above the call wing, and a flat market gives a flat
  smile.

## [1.92.0] - 2026-09-10

### Added
- `calibrate_sabr_lm` (in `sabr.py`): a Levenberg-Marquardt SABR calibrator that
  reuses the exact `sabr_jacobian` from v1.91.0 instead of the derivative-free
  Nelder-Mead. Solves the damped 3x3 normal equations (new `_solve3` Gaussian
  elimination), adapts the damping to guarantee a downhill step, and box-clamps
  (alpha, rho, nu) to the valid region. Returns `(params, rmse, n_iter)`.
- Recovers synthetic (alpha, rho, nu) to ~1e-12 RMSE in ~4 iterations (vs
  Nelder-Mead's ~2e-9) and lands on the same optimum as `calibrate_sabr` on a
  noisy market; also verified for the beta=1 lognormal smile.

## [1.91.0] - 2026-09-10

### Added
- `sabr_sensitivities` and `sabr_jacobian` (in `sabr.py`): exact partial
  derivatives of the Hagan SABR implied vol via a small forward-mode dual-number
  type, with no finite-difference truncation error. `sabr_sensitivities` returns
  the vol plus its partials w.r.t. F, K, alpha, rho and nu; `sabr_jacobian`
  stacks the (alpha, rho, nu) columns into the calibration Jacobian a
  Gauss-Newton / Levenberg-Marquardt step (and the parameter covariance) needs.
- Verified every partial against central finite differences across strikes to a
  max error of ~1e-10.

### Fixed
- The dual-number vol reused a scalar `log(F/K)`, which zeroed the F/K partials;
  `logFK` is now carried as a dual so the smile backbone/skew slopes propagate.

## [1.90.0] - 2026-09-10

### Added
- `heston_qe_mc` (new `heston_mc.py`): Heston Monte Carlo via Andersen's (2008)
  Quadratic-Exponential scheme. The CIR variance is advanced by moment-matching
  to a shifted squared-Gaussian (low vol-of-vol) or an exponential-with-atom
  (high vol-of-vol), so variances stay non-negative by construction; the
  log-asset step uses Andersen's K0..K4 constants with a martingale correction.
  Includes an Acklam inverse-normal `_norm_ppf` for the QE branch draw.
- Verified against the Fourier `heston_price` on four parameter sets (including
  a long-dated, Feller-violating xi=1 case and deep OTM) to within ~1.3 SE, the
  discounted spot is a martingale (forward recovered to 8e-5), and puts match.

## [1.89.0] - 2026-09-10

### Added
- `barrier_mc` (in `montecarlo.py`): Monte Carlo for single-barrier vanilla
  options with a Brownian-bridge crossing correction, the natural cross-check
  for the `barrier_option` closed form including a continuous dividend yield
  (`b = r - q`). Each step contributes the exact conditional probability the
  bridge between its endpoints touched the barrier, removing the discrete-
  monitoring bias that otherwise over-prices knock-outs. Agrees with the
  Reiner-Rubinstein closed form on all four barrier kinds (with q) to a few
  basis points, and knock-in + knock-out reproduces the vanilla.

## [1.88.0] - 2026-09-10

### Added
- `displaced_implied_shift` (in `displaced.py`): calibrate the displaced-diffusion
  shift that reproduces an observed Black-Scholes vol skew. For each trial shift
  the local `sigma` is re-solved to match the ATM quote exactly, so the shift is
  driven purely by the off-ATM skew and is not biased by the local-vs-implied
  vol convention; golden-section search minimises the RMS vol error. Recovers a
  known shift exactly on synthetic quotes (40 -> 40, 100 -> 100, 0 -> 0).

## [1.87.0] - 2026-09-10

### Added
- `merton_smile` (in `merton.py`): the Black-Scholes implied-vol smile a Merton
  jump-diffusion produces — prices calls across strikes and inverts each to a
  BSM vol. Flat at `sigma` with no jumps; symmetric jumps lift the wings into a
  smile and a negative mean jump tilts it into a downward skew.

## [1.86.1] - 2026-09-10

### Added
- Test hardening: a barrier in-out parity property test (`KI + KO = vanilla`
  across all four barrier kinds, calls/puts, and many strike/barrier
  combinations, with and without dividends) and a Bachelier-to-Black-Scholes
  low-vol ATM convergence test.

## [1.86.0] - 2026-09-10

### Added
- `heston_smile` (in `heston.py`): the Black-Scholes implied-vol smile a Heston
  model produces — prices calls across strikes and inverts each to a BSM vol,
  returning `(log_moneyness, vol)` pairs. Flat at `sqrt(v0)` when the vol-of-vol
  is zero; a negative `rho` gives a downward equity skew, positive an upward
  one.

## [1.85.1] - 2026-09-10

### Added
- Benchmark suite now reports implied-vol solver iterations — Newton with the
  Corrado-Miller seed (~4-5 per quote) vs pure bisection (~29), quantifying the
  seed's benefit. Covered by a smoke test.

## [1.85.0] - 2026-09-10

### Added
- `digital_greeks` (in `exotics.py`): delta and gamma of a cash-or-nothing
  digital by finite differences, demonstrating the pin-risk delta spike as
  expiry nears the strike (the motivation for the call-spread over-hedge).
- Regression guard test: the Kirk spread at strike 0 matches the exact Margrabe
  exchange price.

## [1.84.0] - 2026-09-10

### Added
- `implied_spread_correlation` (in `multiasset.py`): back out the correlation
  implied by a spread-option market price, bisecting the Kirk price (monotone
  decreasing in rho) over `(-1, 1)`. Round-trips exactly and raises for quotes
  outside the rho-spanned price range.

## [1.83.0] - 2026-09-10

### Added
- `svi_repair_butterfly` (in `svi.py`): repair a single SVI slice's butterfly
  arbitrage by geometrically shrinking the wing angle `b` (flattening the smile,
  lifting the g-function) until `svi_is_butterfly_free` passes, preserving the
  level/skew/shift/curvature. A clean slice is returned unchanged.

## [1.82.0] - 2026-09-10

### Added
- `strategy_report` (in `strategy.py`): summarize any strategy `Book`'s expiry
  P&L — net premium, max profit / max loss over a terminal-spot grid, whether
  the profit/loss tail is unbounded, and the break-even spots. Works for every
  builder (verticals, straddles, condors, ratio/backspreads, ...).

## [1.81.0] - 2026-09-10

### Added
- `holee.py`: the Ho-Lee (1986) short-rate model — `holee_zero_coupon_bond` and
  `holee_zero_coupon_yield` via the affine closed form (constant-drift case),
  the simplest no-mean-reversion model. Yield is linear in drift with a
  `sigma^2 t^2 / 6` convexity pull-down; MC-verified.

## [1.80.0] - 2026-09-10

### Added
- `basket_greeks` (in `multiasset.py`): Greeks of a two-asset basket option by
  finite differences — both spot deltas, own-gammas, cross-gamma, and
  correlation sensitivity. Both deltas positive, correlation vega positive
  (higher correlation raises the basket vol), and weight rebalancing shifts the
  delta split.

## [1.79.0] - 2026-09-10

### Added
- `bjerksund_stensland_1993` (in `american.py`): the single-flat-boundary
  Bjerksund-Stensland (1993) American approximation — simpler and slightly less
  accurate than the 2002 two-region version, within a few cents of the binomial
  tree. Calls direct, puts via the exact transformation.

## [1.78.0] - 2026-09-10

### Added
- `lee_wing_slopes` and `lee_bounds_ok` (in `svi.py`): the asymptotic total-
  variance wing slopes of an SVI slice, `b(1-rho)` (left) and `b(1+rho)`
  (right), and a check that both satisfy Lee's moment bound (slope <= 2). Agrees
  with `SVIParams.is_arbitrage_free_wings` and matches the empirical
  large-|k| slope.

## [1.77.0] - 2026-09-10

### Added
- `spread_greeks` (in `multiasset.py`): Greeks of a Kirk spread option by finite
  differences — the two spot deltas, own-gammas, cross-gamma, and correlation
  sensitivity. delta1 > 0 / delta2 < 0, cross-gamma near minus the own-gamma,
  correlation vega negative.

## [1.76.0] - 2026-09-10

### Added
- `cir.py`: the Cox-Ingersoll-Ross (1985) short-rate model — `cir_zero_coupon_bond`
  and `cir_zero_coupon_yield` via the affine closed form with the square-root
  diffusion that keeps the rate non-negative. Short yield equals the short rate,
  the long yield approaches the CIR limit `2 kappa theta / (gamma + kappa)`, and
  bond prices match a floored-Euler Monte Carlo.

## [1.75.0] - 2026-09-10

### Added
- `book_bump_greeks` (in `bookgreeks.py`): net book delta/gamma/vega/theta by
  bumping the shared spot/vol/time across every leg and repricing via
  `price_book` — model-free, so it works for any instrument in the book, not
  only ones with analytic Greeks. Matches the analytic net Greeks on a plain
  option book.

## [1.74.0] - 2026-09-10

### Added
- `income.py`: covered-call and cash-secured-put income analytics —
  `covered_call` and `cash_secured_put` return an `IncomeMetrics` with the
  premium (model or supplied), static and annualized yield, the if-assigned
  return, and the breakeven (`S - premium` / `K - premium`).

## [1.73.0] - 2026-09-10

### Added
- `vasicek.py`: the Vasicek (1977) short-rate model — `zero_coupon_bond` and
  `zero_coupon_yield` (affine closed form) and `bond_option` (Jamshidian
  closed-form European option on a zero-coupon bond). The short yield equals the
  short rate, the long yield approaches the mean-reversion level less the
  convexity term, and bond-option parity holds; MC-verified.

## [1.72.0] - 2026-09-10

### Added
- `VolSurface.vol_grid` and `VolSurface.strike_vol_grid`: export a
  `(expiries x strikes)` grid of implied vols from the surface for charting —
  indexed by log-moneyness or by strike (using the carry-implied forward per
  expiry). Returns the grid plus its axes.

## [1.71.0] - 2026-09-10

### Added
- `double_knockout_mc` (in `montecarlo.py`): Monte Carlo a double-knockout
  (corridor) barrier option that pays the vanilla payoff only if the spot stays
  inside `(lower, upper)` for the whole path, else the cash rebate. Wide
  barriers approach the vanilla; a tighter corridor is cheaper.

## [1.70.0] - 2026-09-10

### Added
- `autocallable_mc` (in `montecarlo.py`): Monte Carlo an autocallable note —
  early redemption with accrued coupon when the spot is at or above the autocall
  barrier at an observation date, otherwise notional at maturity with a
  down-and-in downside below the protection barrier. Higher coupons raise the
  value; a protection barrier lowers it.

## [1.69.0] - 2026-09-10

### Added
- `exchange_greeks` (in `multiasset.py`): Greeks of a Margrabe exchange option
  by finite differences — the two spot deltas, own-gammas, the cross-gamma
  (`d2V/dS1 dS2`), and the correlation sensitivity. Cross-gamma equals minus the
  own-gamma (degree-1 homogeneity) and correlation vega is negative.

## [1.68.0] - 2026-09-10

### Added
- `carry_roll_pnl` (in `attribution.py`): the roll-down / carry-roll P&L of an
  option over a horizon at constant vol — rolls the spot to its forward and
  reprices at the shorter maturity, returning a `CarryRoll` with the value now,
  the rolled value, and the roll P&L (theta bleed net of carry drift). Negative
  for a long option, positive for a short.

## [1.67.0] - 2026-09-10

### Added
- `implied_vol_smile` (in `implied.py`): invert a whole option chain to an
  implied-vol smile in one call, returning `(log_moneyness, vol)` pairs sorted
  by strike and dropping any quote outside the no-arbitrage band.

### Changed
- Marked several heavy Monte Carlo tests `slow`, cutting the default
  (`pytest -m "not slow"`) run back to ~12s.

## [1.66.0] - 2026-09-10

### Added
- `VolSurface.fit_arbitrage_free`: fits each expiry with SVI, then walks from
  the short end up and lifts each slice's level just enough to keep total
  variance non-decreasing in maturity — repairing calendar arbitrage. The
  result passes `is_calendar_arbitrage_free`; a consistent surface is left as
  fitted.

## [1.65.0] - 2026-09-09

### Added
- `barrier_rebate` (in `exotics.py`): the standalone rebate cashflow on a
  barrier — a knock-out rebate (pays cash on breach, at-hit or at-expiry) is a
  one-touch, a knock-in rebate (pays at expiry if never breached) is a
  no-touch. Knock-out (expiry) + knock-in sums to the discounted cash.

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
