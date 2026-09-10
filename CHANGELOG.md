# Changelog

All notable changes to QuantForge are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/), and the project adheres to
[Semantic Versioning](https://semver.org/).

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
