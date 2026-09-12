# QuantForge API reference

Auto-generated from `quantforge` v1.449.0 by `docs/gen_api.py` — do not edit by hand.

## actuarial

### `cat_bond_price(principal, coupon_rate, expected_loss_rate, r, maturity, risk_free_spread=0.0)`  _function_

> Present value of a single-period cat bond.
>
> Pays ``coupon_rate`` on the principal and returns the principal at maturity
> unless a triggering event erodes it by the expected loss. Discounts the
> expected principal repayment ``principal * (1 - expected_loss_rate)`` and the
> coupon at ``r + risk_free_spread``. Falls as the expected loss rises.

### `cat_bond_spread(expected_loss_rate, risk_load=1.0)`  _function_

> Fair coupon spread of a cat bond: expected loss rate times a risk load.
>
> Investors demand a spread above the expected loss to bear the (undiversifiable,
> fat-tailed) catastrophe risk: ``spread = (1 + risk_load) * expected_loss_rate``
> (``risk_load = 0`` is the actuarially fair spread). At or above the expected
> loss rate.

### `cat_expected_loss(loss_scenarios, attachment, exhaustion)`  _function_

> Expected layer loss over equally-likely loss scenarios (as a fraction).
>
> Averages :func:`cat_layer_loss` across ``loss_scenarios`` and divides by the
> layer width, giving the expected loss as a fraction of the layer notional in
> ``[0, 1]`` -- the cat bond's expected loss rate.

### `cat_layer_loss(gross_loss, attachment, exhaustion)`  _function_

> Loss ceded to a reinsurance / cat-bond layer ``[attachment, exhaustion]``.
>
> ``min(max(gross_loss - attachment, 0), exhaustion - attachment)`` -- zero below
> the attachment point, rising one-for-one through the layer, capped at the layer
> width above exhaustion.

### `curtate_life_expectancy(one_year_survival)`  _function_

> Curtate expectation of life ``e_x = sum_{k>=1} kp_x`` (whole years).
>
> The expected number of complete future years lived, the sum of the cumulative
> survival probabilities beyond time zero.

### `endowment_insurance(one_year_survival, i, term)`  _function_

> EPV of an endowment: term insurance plus a pure endowment at ``term``.

### `gompertz_makeham_hazard(age, a, b, c)`  _function_

> Gompertz-Makeham force of mortality ``mu(x) = a + b * c^x``.
>
> ``a`` is the age-independent (accident) component and ``b c^x`` the
> exponentially-rising Gompertz term. Increasing in age for ``c > 1``.

### `gompertz_makeham_survival(age, years, a, b, c)`  _function_

> Survival probability over ``years`` under Gompertz-Makeham mortality.
>
> Integrates the force of mortality from ``age`` to ``age + years``:
>
>     tp_x = exp(-a t - (b / ln c) c^x (c^t - 1)),   t = years
>
> (the closed-form integral of ``a + b c^s``). Falls monotonically with the
> horizon; the ``c -> 1`` limit uses the exponential (Makeham-only) form.

### `gompertz_makeham_survival_curve(age, n_years, a, b, c)`  _function_

> One-year survival probabilities ``[p_x, p_{x+1}, ...]`` for ``n_years``.
>
> Each entry is the one-year Gompertz-Makeham survival at successive ages, ready
> to feed the life-table functions (:func:`life_annuity_due`, etc.).

### `life_annuity_due(one_year_survival, i)`  _function_

> Expected present value of a unit life annuity-due.
>
> ``a-due = sum_k v^k * kp_x`` paying 1 at the start of each year while alive,
> ``v = 1/(1+i)``. Falls as interest or mortality rises.

### `net_level_premium(one_year_survival, i, term=None)`  _function_

> Net annual premium for a (term or whole-life) unit insurance.
>
> By the equivalence principle the level premium equates the EPV of premiums
> (a life annuity-due) to the EPV of benefits (the insurance):
>
>     P = A / a-due
>
> Uses :func:`whole_life_insurance` over ``a-due`` for whole life (``term`` None)
> or :func:`endowment_insurance` over the temporary annuity for an ``term``-year
> endowment. The premium the insurer must charge to break even.

### `pure_endowment(one_year_survival, i, term)`  _function_

> EPV of a unit pure endowment: ``v^n * np_x`` (pays 1 iff alive at ``n``).

### `survival_probabilities(one_year_survival)`  _function_

> Cumulative survival ``[0p_x, 1p_x, 2p_x, ...]`` from one-year ``p_x`` values.
>
> ``kp_x = prod_{j<k} p_{x+j}``, starting at ``0p_x = 1``. The returned list has
> one more entry than the input (the leading 1). Survival is non-increasing.

### `temporary_life_annuity_due(one_year_survival, i, term)`  _function_

> EPV of an ``n``-year temporary life annuity-due.
>
> ``a-due_{x:n} = sum_{k<n} v^k * kp_x`` -- pays 1 at the start of each year while
> alive, for at most ``term`` years. Below the whole-life
> :func:`life_annuity_due` and rising to it as ``term`` grows.

### `term_insurance(one_year_survival, i, term=None)`  _function_

> EPV of a unit term insurance paying 1 at the end of the year of death.
>
> ``A = sum_k v^{k+1} * kp_x * q_{x+k}`` over the first ``term`` years (default:
> the whole table). ``q = 1 - p`` is the one-year death probability.

### `whole_life_insurance(one_year_survival, i)`  _function_

> EPV of whole-life insurance: :func:`term_insurance` over the whole table.

## american

### `bjerksund_stensland(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None)`  _function_

> American option price via Bjerksund-Stensland (2002), closed form.
>
> Args mirror the rest of the engine. ``b`` is the cost of carry (defaults to
> ``r``); dividend yield q enters as b = r - q. American puts are priced via
> the exact put-call transformation P(S,K,r,b) = C(K,S,r-b,-b).

### `bjerksund_stensland_1993(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None)`  _function_

> American option price via Bjerksund-Stensland (1993), single flat boundary.
>
> A simpler and slightly less accurate predecessor to the 2002 two-region
> model (:func:`bjerksund_stensland`): it uses one flat exercise boundary. Calls
> are priced directly; puts via the exact transformation
> ``P(S,K,r,b) = C(K,S,r-b,-b)``.

### `bjerksund_stensland_boundary(K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None)`  _function_

> Bjerksund-Stensland (2002) flat exercise trigger ``I`` at inception.
>
> Immediate exercise is optimal for a call at ``S >= I`` and for a put at
> ``S <= I``. Returns ``None`` when early exercise is never optimal (an
> American call with ``b >= r``). This is the model's flat-boundary
> approximation to the true (curved) early-exercise frontier -- the level
> where the BS2002 price equals the exercise intrinsic.
>
> The put trigger follows from the same put-call transformation used by the
> pricer, ``P(S, K, r, b) = C(K, S, r - b, -b)``: the transformed call's
> spot-axis trigger ``I2t`` maps back to the put boundary ``K^2 / I2t``.

### `bjerksund_stensland_greeks(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None)`  _function_

> Greeks of the Bjerksund-Stensland American price by finite differences.
>
> The 2002 price is a closed form but its Greeks have no simple expression
> (the exercise boundary and the bivariate-normal term move with the inputs),
> so we central-difference the price. Returns a dict with delta, gamma, vega,
> theta (per year, calendar), and rho.
>
> Bumps are chosen small relative to each input; because the BS2002 price is a
> smooth function of its arguments (away from t=0) central differences are
> accurate to a few basis points, plenty for hedging.

### `early_exercise_premium(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None)`  _function_

> Decompose the American price into European value + early-exercise premium.
>
> Returns a dict with ``american`` (Bjerksund-Stensland), ``european`` (BSM),
> and ``premium`` = american - european, the extra value from the right to
> exercise early. The premium is non-negative and is (near) zero for an
> American call with no dividends (``b >= r``), where early exercise is never
> optimal.

## andreasenhuge

### `andreasen_huge_calibrate(F, strikes, T, market_vols, r=0.0, max_iter=60, tol=1e-08)`  _function_

> Calibrate per-strike local vols so the AH smile matches market vols.
>
> Bootstraps the local vol at each interior strike by a 1-D bisection on the
> single-strike implied-vol error (the implicit step couples neighbours only
> weakly, so a few sweeps converge). Returns ``(local_vols, rmse)`` with
> ``rmse`` the root-mean-square implied-vol error over the interior strikes.

### `andreasen_huge_prices(F, strikes, T, local_vols)`  _function_

> Arbitrage-free forward call prices from one implicit Dupire step.
>
> ``strikes`` is an increasing grid, ``local_vols`` the per-strike local vol
> (same length). Returns the list of undiscounted forward call prices
> ``C(K_i)`` at expiry ``T``. Multiply by ``e^{-rT}`` for the discounted price
> if the forward already embeds the carry.
>
> The prices are monotone decreasing and convex in strike by construction, for
> any positive ``local_vols`` -- the whole point of the scheme.

### `andreasen_huge_smile(F, strikes, T, local_vols, r=0.0)`  _function_

> Implied-vol smile from the Andreasen-Huge arbitrage-free call prices.
>
> Prices the grid with :func:`andreasen_huge_prices` (spot ``S = F e^{-rT}``,
> carry ``b = r`` so the pricing forward is ``F``) and inverts each to a
> Black-Scholes implied vol. Returns ``(log_moneyness, vol)`` pairs sorted by
> strike on the forward ``F``.

### `andreasen_huge_strike_greeks(F, strikes, T, local_vols, r=0.0)`  _function_

> Strike-space Greeks of the Andreasen-Huge call surface.
>
> From the arbitrage-free forward call prices (:func:`andreasen_huge_prices`),
> computes at each interior strike:
>
>   * ``dual_delta`` = ``dC/dK`` (discounted), which equals ``-e^{-rT}`` times
>     the risk-neutral probability of finishing above ``K``; monotone in
>     ``[-e^{-rT}, 0]``;
>   * ``rnd`` = ``e^{rT} d2C/dK2``, the Breeden-Litzenberger risk-neutral
>     density, non-negative by the scheme's convexity.
>
> Returns ``(interior_strikes, dual_delta, rnd)`` as three equal-length lists
> (the two Dirichlet edge strikes are dropped). The density is non-negative for
> any positive ``local_vols`` and integrates to approximately 1 over the grid.

## attribution

### `CarryRoll(horizon: float, forward_spot: float, value_now: float, value_rolled_static: float, theta_roll: float) -> None`  _class_

> CarryRoll(horizon: float, forward_spot: float, value_now: float, value_rolled_static: float, theta_roll: float)

### `PnLAttribution(total: float, delta_pnl: float, gamma_pnl: float, vega_pnl: float, theta_pnl: float, rho_pnl: float, explained: float, unexplained: float) -> None`  _class_

> PnLAttribution(total: float, delta_pnl: float, gamma_pnl: float, vega_pnl: float, theta_pnl: float, rho_pnl: float, explained: float, unexplained: float)

### `attribute_pnl(S, K, t, r, sigma, dS, dsigma, dt, dr=0.0, option_type=<OptionType.CALL: 'call'>, b=None, qty=1.0) -> quantforge.attribution.PnLAttribution`  _function_

> Attribute an option position's P&L over a move to its Greeks.
>
> Args:
>     dS: change in spot. dsigma: change in vol. dt: elapsed calendar time
>         (years). dr: change in rate.
>     qty: signed position size (scales every P&L component).
>
> theta here is the calendar theta (per year) from the engine, so the theta
> P&L is ``theta * dt`` (value lost as time passes). vega is per 1.0 vol, rho
> per 1.0 rate; pass dsigma / dr in those units.

### `carry_roll_pnl(S, K, t, r, sigma, horizon, option_type=<OptionType.CALL: 'call'>, b=None, qty=1.0)`  _function_

> Roll-down / carry-roll P&L of an option over ``horizon`` at constant vol.
>
> Rolls the position forward by ``horizon`` years assuming the spot drifts to
> its forward ``S e^{b*horizon}`` and volatility is unchanged, then reprices at
> the shorter remaining maturity. The roll P&L is the change in value -- the
> theta bleed net of the forward drift the carry earns. This is the standard
> "if nothing moves, what do I earn/pay" carry number.
>
> Returns a :class:`CarryRoll`.

## bachelier

### `bachelier_asset_or_nothing(F, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>) -> float`  _function_

> Bachelier asset-or-nothing digital: pays the forward ``F_T`` if in the money.
>
> With ``F_T`` Gaussian, ``E[F_T 1_{F_T > K}] = F N(d) + sigma sqrt(t) phi(d)``
> (call), discounted at ``r``. Note the vanilla Bachelier call equals this
> asset-or-nothing minus ``K`` times the cash-or-nothing, mirroring the
> Black-Scholes decomposition.

### `bachelier_cash_or_nothing(F, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, cash=1.0) -> float`  _function_

> Bachelier cash-or-nothing digital: pays ``cash`` if in the money.
>
> In the normal model ``F_T`` is Gaussian, so with ``d = (F - K)/(sigma sqrt t)``
> a call (pays when ``F_T > K``) is ``cash e^{-rt} N(d)`` and a put is
> ``cash e^{-rt} N(-d)``.

### `bachelier_delta(F, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>) -> float`  _function_

> dPrice/dF (in the forward). Call delta is e^{-rt} N(d).

### `bachelier_gamma(F, K, t, r, sigma) -> float`  _function_

> d2Price/dF2. Same for calls and puts.

### `bachelier_greeks(F, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>)`  _function_

> Bundle the Bachelier Greeks: delta, gamma, vega (analytic) plus theta.
>
> ``delta``, ``gamma``, and ``vega`` reuse the exact closed forms
> :func:`bachelier_delta`, :func:`bachelier_gamma`, :func:`bachelier_vega`;
> ``theta`` (calendar decay, ``-dV/dt``) is a central finite difference of
> :func:`bachelier_price`. Returns a dict with ``price``, ``delta``, ``gamma``,
> ``vega``, ``theta``. All are in normal-model (absolute-vol) terms.

### `bachelier_implied_vol(target_price, F, K, t, r, option_type=<OptionType.CALL: 'call'>, tol=1e-10, max_iter=100) -> float`  _function_

> Solve for the normal volatility that reproduces ``target_price``.
>
> Newton's method on vega with a bisection fallback. Rejects prices outside
> the no-arbitrage band [intrinsic, forward-bound].

### `bachelier_price(F, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>) -> float`  _function_

> Bachelier price of a European option on a forward ``F``.
>
> ``sigma`` is the normal (absolute) volatility. ``r`` discounts the payoff
> from expiry; pass ``r=0`` to price on the forward directly.

### `bachelier_theta(F, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>) -> float`  _function_

> Calendar theta ``-dPrice/dt`` in the Bachelier (normal) model, analytic.
>
> Differentiating the discounted normal price gives, for a call,
>
>     theta = r * price - e^{-rt} [ sigma phi(d) / (2 sqrt(t)) ],
>
> where ``d = (F - K)/(sigma sqrt(t))``. The normal-model time value grows
> with maturity (the ``sigma phi/(2 sqrt t)`` piece, same for calls and puts,
> enters ``-dP/dt`` with a minus sign); the ``r * price`` piece is the
> discount drift. At ``r = 0`` this is the pure decay
> ``-sigma phi(d)/(2 sqrt(t))``.

### `bachelier_vega(F, K, t, r, sigma) -> float`  _function_

> dPrice/dsigma_N (per unit of normal vol). Same for calls and puts.

## bates

### `bates_greeks(S, K, t, r, v0, kappa, theta, xi, rho, lam, mu_j, sigma_j, option_type=<OptionType.CALL: 'call'>, q=0.0)`  _function_

> Greeks of a Bates (Heston + jumps) option by central finite differences.
>
> Central differences of :func:`bates_price` for the spot Greeks ``delta``
> (dV/dS) and ``gamma`` (d2V/dS2), the initial-variance sensitivity ``vega_v0``
> (dV/dv0 -- the stochastic-vol analogue of vega), and the jump-intensity
> sensitivity ``d_lambda`` (dV/dlam). At ``lam = 0`` the Greeks reduce to the
> Heston Greeks. Returns a dict with ``price``, ``delta``, ``gamma``,
> ``vega_v0``, ``d_lambda``.

### `bates_price(S, K, t, r, v0, kappa, theta, xi, rho, lam, mu_j, sigma_j, option_type=<OptionType.CALL: 'call'>, q=0.0, upper=200.0) -> float`  _function_

> Price a European option under the Bates (Heston + Merton jumps) model.
>
> Args:
>     v0, kappa, theta, xi, rho: Heston stochastic-variance parameters.
>     lam: jump intensity (expected jumps per year, >= 0).
>     mu_j, sigma_j: mean and std of the log jump size.
>     q: continuous dividend yield.
>     upper: Fourier-integral truncation (raise for long maturities/large xi).
>
> ``lam = 0`` recovers the Heston price exactly. Puts use put-call parity.

### `bates_smile(S, strikes, t, r, v0, kappa, theta, xi, rho, lam, mu_j, sigma_j, q=0.0)`  _function_

> Black-Scholes implied-vol smile the Bates model produces.
>
> Prices a call at each strike and inverts to a Black-Scholes implied vol,
> returning ``(log_moneyness, vol)`` pairs sorted by strike on the forward
> ``F = S e^{(r-q) t}``. Jumps steepen the short-dated skew beyond what the
> Heston diffusion alone can produce.

## baw

### `baw_american(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None) -> float`  _function_

> American option price by the Barone-Adesi-Whaley quadratic approximation.
>
> ``b`` is the cost of carry (defaults to ``r``); dividend yield ``q`` enters
> as ``b = r - q``. A no-dividend American call (``b = r``) returns the
> European value. Falls back to intrinsic below/above the critical spot.

### `baw_american_greeks(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None)`  _function_

> Greeks of a BAW American option by central finite differences of
> :func:`baw_american`: ``delta`` (dV/dS), ``gamma`` (d2V/dS2), ``vega``
> (dV/dsigma), ``theta`` (calendar decay). Returns a dict with ``price`` and
> those fields.

### `baw_critical_spot(K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None)`  _function_

> Barone-Adesi-Whaley early-exercise boundary S* at inception.
>
> The spot at which immediate exercise becomes optimal: a call is exercised
> for ``S >= S*`` and a put for ``S <= S*``. Returns ``None`` when early
> exercise is never optimal (an American call with ``b >= r`` equals its
> European value, so there is no finite boundary).

## bermudan_swaption

### `bermudan_swaption_g2pp(P0, exercise_times, fixed_rate, a, b, sigma, eta, rho, payer=True, n_paths=20000, seed=None)`  _function_

> Price a Bermudan swaption under G2++ by Longstaff-Schwartz.
>
> Args:
>     P0: today's discount-factor function ``P(0, T)``.
>     exercise_times: increasing list of exercise dates; at date ``t_k`` the
>         holder may enter the co-terminal swap paying/receiving ``fixed_rate``
>         on the remaining schedule ``exercise_times[k:]`` (annual periods).
>     payer: True for a payer swaption (pay fixed), else receiver.
>     a, b, sigma, eta, rho: G2++ parameters.
>
> Returns the Bermudan swaption price (today's value). Uses a
> money-market numeraire built from one-period bonds along each path.

### `bermudan_swaption_g2pp_greeks(P0, exercise_times, fixed_rate, a, b, sigma, eta, rho, payer=True, n_paths=20000, seed=None)`  _function_

> Greeks of a G2++ Bermudan swaption by common-random-number bumps.
>
> Reprices :func:`bermudan_swaption_g2pp` on the *same* seed at bumped inputs,
> so the two simulations share their G2++ state paths and the finite
> differences are low-variance. Returns a dict with ``price`` and:
>
>   * ``d_fixed`` = dV/d(fixed_rate) -- negative for a payer (paying a higher
>     fixed rate is worth less), positive for a receiver;
>   * ``curve_dv01`` = the value change for a 1bp parallel *drop* in the
>     zero curve (applied as ``P0(T) -> P0(T) e^{+1e-4 T}``), the standard
>     DV01 sign convention (positive for a payer).
>
> A fixed ``seed`` is required for the CRN differences to be meaningful.

## binomial

### `american_price(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None, steps=500)`  _function_

> Price an American option via a CRR binomial tree.
>
> Args:
>     steps: number of time steps. Higher = more accurate, O(steps^2) work.
>     b: cost of carry (defaults to r). Dividend yield q enters as b = r - q.

## bkm

### `bkm_moments_from_smile(S0, t, r, vol_fn, q=0.0, n_strikes=401, width=8.0)`  _function_

> Risk-neutral (variance, skewness, excess kurtosis) via Bakshi-Kapadia-Madan.
>
> ``vol_fn(K)`` is the implied-vol smile; OTM options are priced with
> Black-Scholes and the three moment contracts integrated by the trapezoidal
> rule. Returns ``(variance, skewness, excess_kurtosis)`` of the ``t``-horizon
> log-return under the risk-neutral measure.

### `skew_swap_from_smile(S0, t, r, vol_fn, q=0.0, n_strikes=401, width=8.0)`  _function_

> Fair skew-swap value: the risk-neutral skewness from the BKM moments.

## bond_future

### `bond_future_dv01(price, modified_duration)`  _function_

> Dollar value of a 1bp yield rise for a bond: ``duration * price * 1e-4``.
>
> Positive magnitude of the price move per basis point (a price *fall* for a
> yield rise). Reported as a positive number for hedge sizing.

### `cheapest_to_deliver(bonds, futures_price)`  _function_

> Index of the cheapest-to-deliver bond (minimum net basis).
>
> ``bonds`` is a list of dicts with keys ``price``, ``cf`` (conversion factor),
> and ``carry``. Returns the index minimizing :func:`net_basis`. Ties break to
> the first.

### `conversion_factor(coupon_rate, years_to_maturity, notional_coupon=0.06, freq=2)`  _function_

> Conversion factor: price per unit face of a bond yielding the notional coupon.
>
> The bond's price discounted at the exchange's notional coupon
> (``notional_coupon``, e.g. 6% for CBOT Treasury futures), per unit face:
>
>     CF = sum_i (c/freq) / (1 + y/freq)^i + 1 / (1 + y/freq)^n,
>
> with ``c = coupon_rate``, ``y = notional_coupon``, ``n = years * freq``. Above
> one for bonds with coupons over the notional, below one otherwise; exactly one
> when the coupon equals the notional.

### `futures_dv01(ctd_dv01, ctd_conversion_factor)`  _function_

> DV01 of a bond future from the CTD bond's DV01.
>
> ``ctd_dv01 / conversion_factor`` -- the futures price moves ``1/CF`` of the
> CTD price per unit yield (the conversion factor gears the delivery), so the
> futures DV01 is the CTD DV01 divided by its conversion factor.

### `futures_hedge_ratio(bond_dv01_, futures_dv01_)`  _function_

> Number of futures to hedge a cash bond's rate risk.
>
> ``bond_dv01 / futures_dv01`` -- the contract count whose DV01 offsets the
> bond's. Positive; a long bond position is hedged by selling this many futures.

### `gross_basis(bond_price, futures_price, conversion_factor_)`  _function_

> Gross basis ``bond_price - futures_price * CF``.
>
> The raw richness of the cash bond over its futures-implied forward price
> (before carry). Positive means the bond trades above its converted futures
> price.

### `implied_repo_rate(bond_price, accrued_now, futures_price, conversion_factor_, accrued_delivery, days, day_count=360)`  _function_

> Implied repo rate of a deliverable bond.
>
> The financing rate that makes buying the bond and delivering into the future
> break even:
>
>     repo = (invoice - dirty_now) / dirty_now * (day_count / days),
>
> with ``invoice = futures * CF + accrued_delivery`` and ``dirty_now =
> bond_price + accrued_now``. The cheapest-to-deliver bond has the highest
> implied repo.

### `invoice_price(futures_price, conversion_factor_, accrued_interest=0.0)`  _function_

> Invoice price the short receives: ``futures * CF + accrued``.

### `net_basis(bond_price, futures_price, conversion_factor_, carry)`  _function_

> Net basis: gross basis less the carry to delivery.
>
> ``gross_basis - carry`` where ``carry`` is coupon income minus financing over
> the holding period. The cheapest-to-deliver bond minimizes the net basis.

## bondmath

### `accrued_interest(face, coupon_rate, freq, fraction_elapsed) -> float`  _function_

> Accrued interest since the last coupon, straight-line within the period.
>
> ``fraction_elapsed`` in ``[0, 1]`` is the share of the current coupon period
> that has passed at settlement. The accrual is
> ``face * coupon_rate / freq * fraction_elapsed`` -- the linear (actual/
> actual-in-period) convention. The buyer pays this on top of the quoted
> clean price.

### `bond_cashflows(face, coupon_rate, maturity, freq=2) -> List[Tuple[float, float]]`  _function_

> Level-coupon schedule ``[(t, amount), ...]`` for a vanilla bond.
>
> ``freq`` coupons per year of ``face * coupon_rate / freq`` each, with the
> face repaid alongside the final coupon at ``maturity``.

### `bond_dv01(cashflows: Sequence[Tuple[float, float]], y) -> float`  _function_

> Dollar value of a 1bp yield rise (negative: prices fall as yields rise).
>
> ``DV01 = dPrice/dy * 1e-4 = -modified_duration * price * 1e-4``.

### `bond_price_from_yield(cashflows: Sequence[Tuple[float, float]], y) -> float`  _function_

> Present value of the cashflows at continuously-compounded yield ``y``.

### `clean_price(cashflows: Sequence[Tuple[float, float]], y, face, coupon_rate, freq, fraction_elapsed) -> float`  _function_

> Clean (quoted) price: dirty price minus accrued interest.
>
> ``clean = dirty - accrued``. At a coupon date (``fraction_elapsed = 0``) the
> clean and dirty prices coincide; mid-period the clean price strips out the
> accrued coupon so the quote does not saw-tooth across coupon dates.

### `convexity(cashflows: Sequence[Tuple[float, float]], y) -> float`  _function_

> Convexity ``1/price d2Price/dy2``: PV-weighted average of squared time.

### `dated_accrued_interest(settle, prev_coupon_date, next_coupon_date, face, coupon_rate, freq, convention='30/360') -> float`  _function_

> Accrued interest at ``settle`` between two coupon dates, day-count based.
>
> Accrues the current coupon ``face * coupon_rate / freq`` by the fraction of
> the period elapsed, ``year_fraction(prev, settle) / year_fraction(prev,
> next)``, under the day-count ``convention``. Zero at the coupon date, the
> full coupon just before the next.

### `dated_bond_cashflows(start, maturity_years, face, coupon_rate, freq=2, convention='30/360', end_of_month=False, business_day='unadjusted')`  _function_

> Coupon-bond cashflows on a real calendar with day-count accruals.
>
> Generates the coupon schedule from ``start`` (a ``(y, m, d)`` date) with
> :func:`quantforge.generate_schedule`, computes each period's accrual factor
> with :func:`quantforge.year_fraction` under ``convention``, and pays the
> day-count-weighted coupon ``face * coupon_rate * tau_i`` each period plus the
> face at maturity. Returns ``[(pay_date, year_fraction, amount), ...]``.
>
> Unlike :func:`bond_cashflows` (which assumes uniform ``1/freq`` periods),
> this reflects the actual day counts, month-end roll, and business-day
> adjustment -- the difference that matters for act/360 and stub periods.

### `dated_bond_price(settle, dated_cashflows, y, convention='30/360') -> float`  _function_

> Present value of dated cashflows discounted from a settlement date.
>
> ``dated_cashflows`` is the ``[(pay_date, tau, amount), ...]`` output of
> :func:`dated_bond_cashflows`. Each amount is discounted by
> ``exp(-y * T_i)`` where ``T_i`` is the :func:`quantforge.year_fraction` from
> ``settle`` to the pay date under ``convention`` (continuously-compounded
> yield). Cashflows on or before ``settle`` are dropped.

### `dated_bond_yield(settle, dated_cashflows, price, convention='30/360', tol=1e-10, max_iter=100) -> float`  _function_

> Continuously-compounded yield reproducing a dated bond ``price``.
>
> Bisection on the yield (price is monotone decreasing in it), discounting the
> :func:`dated_bond_price` cashflows from ``settle``. Inverse of
> :func:`dated_bond_price`.

### `dated_clean_price(settle, dated_cashflows, y, prev_coupon_date, next_coupon_date, face, coupon_rate, freq, convention='30/360') -> float`  _function_

> Clean (quoted) dated price: dirty price minus dated accrued interest.
>
> Combines :func:`dated_bond_price` (the dirty/invoice price discounted from
> ``settle``) with :func:`dated_accrued_interest`. At a coupon date the clean
> and dirty prices coincide.

### `dirty_price(cashflows: Sequence[Tuple[float, float]], y) -> float`  _function_

> Dirty (invoice) price: the full present value of the remaining cashflows.
>
> Alias of :func:`bond_price_from_yield` -- the cash amount actually paid at
> settlement, before subtracting accrued interest to get the clean quote.

### `effective_duration_from_curve(cashflows, pillar_times, zero_rates, bump=0.0001)`  _function_

> Effective duration under a parallel shift of the whole zero curve.
>
> Shifts every pillar zero rate by +/- ``bump`` and central-differences the
> fractional price change. Equals the sum of the :func:`key_rate_durations` to
> first order.

### `key_rate_durations(cashflows, pillar_times, zero_rates, bump=0.0001)`  _function_

> Key-rate (partial) durations of a bond against a zero-rate curve.
>
> Builds a log-linear :class:`quantforge.DiscountCurve` from the pillar zero
> rates, then bumps each pillar's zero rate up by ``bump`` in turn and measures
> the fractional price change ``-dP/P / bump``. Returns a list of key-rate
> durations aligned with ``pillar_times``; their sum approximates the bond's
> effective duration (a parallel shift is the sum of the pillar bumps).
>
> Continuously-compounded zero rates; ``DF = e^{-z t}`` at each pillar.

### `macaulay_duration(cashflows: Sequence[Tuple[float, float]], y) -> float`  _function_

> Macaulay duration (years): PV-weighted average cashflow time.

### `modified_duration(cashflows: Sequence[Tuple[float, float]], y) -> float`  _function_

> Modified duration ``-1/price dPrice/dy``.
>
> Under continuous compounding this equals the Macaulay duration.

### `price_from_curve(cashflows: Sequence[Tuple[float, float]], curve) -> float`  _function_

> Present value of the cashflows off a discount curve.
>
> ``curve`` is anything callable as ``curve.df(t)`` (e.g.
> :class:`quantforge.DiscountCurve`) or a plain ``curve(t)`` returning the
> discount factor ``P(0, t)``.

### `yield_to_maturity(cashflows: Sequence[Tuple[float, float]], price, tol=1e-10, max_iter=100) -> float`  _function_

> Solve the continuously-compounded yield reproducing ``price``.
>
> Newton on the price/yield relation (derivative is ``-D * price``) with a
> bracketing bisection fallback, since price is monotone decreasing in yield.

## bookgreeks

### `BookSecondOrder(vanna: float = 0.0, vomma: float = 0.0, charm: float = 0.0, veta: float = 0.0, speed: float = 0.0, zomma: float = 0.0, color: float = 0.0) -> None`  _class_

> BookSecondOrder(vanna: float = 0.0, vomma: float = 0.0, charm: float = 0.0, veta: float = 0.0, speed: float = 0.0, zomma: float = 0.0, color: float = 0.0)

### `BumpGreeks(price: float, delta: float, gamma: float, vega: float, theta: float) -> None`  _class_

> BumpGreeks(price: float, delta: float, gamma: float, vega: float, theta: float)

### `ThetaCarry(theta: float, gamma_rent: float, residual: float) -> None`  _class_

> ThetaCarry(theta: float, gamma_rent: float, residual: float)

### `book_bump_greeks(contracts, dS_frac=0.001, dvol=0.0001, dt=0.0001)`  _function_

> Net book Greeks by bumping the shared market and repricing (model-free).
>
> Applies a common shock to every leg's spot, volatility, and time-to-expiry
> and reprices the whole book via :func:`quantforge.price_book`, so the net
> delta/gamma/vega/theta come out numerically without needing analytic Greeks
> for each instrument. Assumes all legs share one underlying and vol (a
> single-name book), the usual case for this kind of check.
>
> ``dS_frac`` is the relative spot bump; ``dvol`` and ``dt`` are absolute.
> Returns a :class:`BumpGreeks`.

### `book_second_order(contracts: Iterable[quantforge.portfolio.Contract]) -> quantforge.bookgreeks.BookSecondOrder`  _function_

> Aggregate position-scaled second-order Greeks across a book.
>
> charm/veta/color are calendar-convention (per year), matching the scalar
> functions in :mod:`quantforge.greeks2`.

### `theta_carry_report(contracts)`  _function_

> Decompose a book's net theta into gamma-rent and a residual carry term.
>
> The theta-gamma relationship says an option's time decay is dominated by the
> "gamma rent" paid on convexity: for a single underlying with volatility
> ``sigma``, ``theta ~= -0.5 * Gamma * sigma^2 * S^2`` plus a smaller
> drift/financing residual (rho- and dividend-carry effects).
>
> Sums the position-scaled net theta and the gamma-rent term across the book;
> the residual is the difference. When every leg shares one spot/vol (a
> single-name book) the gamma-rent uses that common S and sigma.

## bsm

### `Greeks(price: float, delta: float, gamma: float, vega: float, theta: float, rho: float) -> None`  _class_

> Greeks(price: float, delta: float, gamma: float, vega: float, theta: float, rho: float)

### `OptionType(*values)`  _class_

> str(object='') -> str
> str(bytes_or_buffer[, encoding[, errors]]) -> str
>
> Create a new string object from the given object. If encoding or
> errors is specified, then the object must expose a data buffer
> that will be decoded using the given encoding and error handler.
> Otherwise, returns the result of object.__str__() (if defined)
> or repr(object).
> encoding defaults to 'utf-8'.
> errors defaults to 'strict'.

### `box_spread_implied_rate(box_price, K1, K2, t) -> float`  _function_

> Continuously-compounded rate implied by a box-spread price.
>
> A box (bull call spread + bear put spread on strikes ``K1 < K2``) has the
> riskless terminal payoff ``K2 - K1``, so its fair value is
> ``e^{-rt} (K2 - K1)`` and the implied financing rate is
>
>     r = -ln(box_price / (K2 - K1)) / t.
>
> This is the synthetic-lending rate the options market prices, independent of
> the underlying. Requires ``0 < box_price < K2 - K1`` (a positive rate) or
> allows a negative rate if the box trades above its notional.

### `call_price(S, K, t, r, sigma, b=None) -> float`  _function_

> (no docstring)

### `delta(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None) -> float`  _function_

> dPrice/dS.

### `epsilon(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None) -> float`  _function_

> Dividend rho (a.k.a. epsilon / psi): dPrice/dq, per 1.0 change in the
> continuous dividend yield.
>
> The dividend yield enters through the carry ``b = r - q``, so raising q
> lowers the forward. For a call ``dPrice/dq = -S t e^{(b-r)t} N(d1)``; for a
> put ``+S t e^{(b-r)t} N(-d1)``. Assumes ``b`` moves with ``q`` (the standard
> dividend-yield case).

### `forward_price(S, t, r, b=None) -> float`  _function_

> Forward price of the underlying ``F = S e^{b t}`` (carry ``b``, default r).
>
> With ``b = r`` this is the cost-of-carry forward on a non-dividend stock;
> ``b = r - q`` handles a continuous dividend yield and ``b = 0`` a future.

### `gamma(S, K, t, r, sigma, b=None) -> float`  _function_

> d2Price/dS2. Identical for calls and puts.

### `greeks(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None) -> quantforge.bsm.Greeks`  _function_

> Compute price and all first/second-order Greeks in one call.

### `implied_discount_factor(call1, put1, K1, call2, put2, K2) -> float`  _function_

> Discount factor implied by call-put pairs at two strikes.
>
> Put-call parity at each strike is ``C_i - P_i = DF (F - K_i)``; subtracting
> the two eliminates the forward and leaves
>
>     DF = [ (C1 - P1) - (C2 - P2) ] / (K2 - K1).
>
> A model-free read of the discount factor to expiry straight off two
> same-expiry option pairs, needing neither a rate nor a volatility. Combine
> with :func:`implied_forward_from_parity` to also back out the forward.

### `implied_forward_from_parity(call, put, K, t, r) -> float`  _function_

> Forward price implied by a call-put pair at strike ``K`` (invert parity).
>
> Put-call parity ``C - P = e^{-rt}(F - K)`` solves for the forward
> ``F = K + e^{rt}(C - P)`` -- the market's forward read straight off a
> same-strike call and put, with no volatility input. Combined across two
> strikes it also pins the implied discount factor.

### `price(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None) -> float`  _function_

> Price a European option under the generalized BSM model.
>
> Args:
>     S: spot price of the underlying.
>     K: strike price.
>     t: time to expiry in years.
>     r: continuously-compounded risk-free rate.
>     sigma: annualized volatility.
>     option_type: CALL or PUT (also accepts "call"/"put"/"c"/"p").
>     b: cost of carry. Defaults to ``r`` (non-dividend stock).

### `put_call_parity_residual(call, put, S, K, t, r, b=None) -> float`  _function_

> Put-call parity residual ``(C - P) - e^{-rt}(F - K)``, ``F = S e^{b t}``.
>
> Zero (up to rounding) when the call and put are arbitrage-consistent. A
> non-zero value is the parity violation in price terms -- useful for
> validating quotes or a pricer. For a non-dividend stock (``b = r``) the
> forward term reduces to ``S - K e^{-rt}``.

### `put_price(S, K, t, r, sigma, b=None) -> float`  _function_

> (no docstring)

### `rho(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None) -> float`  _function_

> dPrice/dr, per 1.0 change in rate.
>
> Assumes carry moves with the rate (the plain BSM stock case). For models
> where ``b`` is fixed independently of ``r`` (e.g. Black-76), pass ``b`` and
> interpret accordingly.

### `rho_discount(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None) -> float`  _function_

> Discount-only rho ``dPrice/dr`` holding the carry ``b`` fixed.
>
> When the cost of carry is independent of the funding rate (Black-76 on a
> future, an FX or commodity forward), a change in ``r`` moves only the
> discount factor, not the forward. Then ``dPrice/dr = -t * price`` for both
> calls and puts. Contrast :func:`rho`, which assumes ``b`` moves with ``r``
> (the plain stock case).

### `theta(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None) -> float`  _function_

> Calendar-time theta, dPrice/d(calendar time) per year.
>
> This is the market convention: equal to ``-dPrice/dt_expiry``, so long
> options usually show negative theta (value decays as the clock advances).
> Divide by 365 for per-calendar-day decay.

### `vega(S, K, t, r, sigma, b=None) -> float`  _function_

> dPrice/dSigma, per 1.0 change in vol (divide by 100 for per-vol-point).

## carrmadan

### `carr_madan_smile_strip(S, t, r, q, psi, k_lo=-0.5, k_hi=0.5, alpha=1.5, n_fft=4096, eta=0.25)`  _function_

> Implied-vol smile over a log-moneyness window from one Carr-Madan FFT.
>
> Runs :func:`carr_madan_strip` once, keeps the grid strikes whose forward
> log-moneyness ``ln(K / F)`` lies in ``[k_lo, k_hi]``, and inverts each call
> to a Black-Scholes implied vol. Returns ``(log_moneyness, vol)`` pairs sorted
> by strike -- the whole smile from a single transform.

### `carr_madan_strip(S, t, r, q, psi, alpha=1.5, n_fft=4096, eta=0.25)`  _function_

> Price a whole strip of European calls in one FFT (Carr-Madan 1999).
>
> Returns ``(strikes, calls)`` on a log-strike grid centred on the forward.
> The Carr-Madan damped-call transform is sampled at ``n_fft`` frequency points
> spaced ``eta`` apart, Simpson-weighted, and inverted with a single radix-2
> FFT -- so the entire smile costs one transform instead of one Gauss-Legendre
> integral per strike. The log-strike spacing is ``lambda = 2 pi / (n_fft eta)``.
>
> Args:
>     psi: characteristic exponent ``psi(u)`` (as in :func:`levy_price`).
>     alpha: damping factor (> 0); needs ``psi(-(alpha+1) i)`` finite.
>     n_fft: FFT length (power of two).
>     eta: frequency-grid spacing; smaller = finer strikes over a wider range.

### `cos_greeks(S, K, t, r, q, psi, option_type=<OptionType.CALL: 'call'>, n_terms=256, L=12.0, cumulants=None)`  _function_

> Delta and gamma of a COS-method price, analytic in the cosine series.
>
> The COS price is ``disc * sum_k Re(cf_k e^{-i u_k a}) V_k`` where the model
> dependence on spot enters only through ``cf_k = e^{i u_k (x + mu) + t psi}``
> with ``x = ln(S/K)``. Differentiating that sum term-by-term w.r.t. ``S`` (so
> ``dx/dS = 1/S``) gives delta and gamma with no re-pricing and no finite
> differences: each term picks up ``i u_k / S`` for delta and
> ``i u_k (i u_k - 1)/S^2`` for gamma. Returns ``{price, delta, gamma}``.

### `cos_price(S, K, t, r, q, psi, option_type=<OptionType.CALL: 'call'>, n_terms=256, L=12.0, cumulants=None) -> float`  _function_

> Price a European option by the COS method (Fang & Oosterlee, 2008).
>
> Expands the risk-neutral density of the log-return in a Fourier-cosine series
> on a truncation range ``[a, b]``, so the price is a finite sum of the
> characteristic function sampled at ``k pi / (b - a)`` against closed-form
> payoff coefficients. Exponentially convergent in ``n_terms`` for smooth
> densities -- an independent Fourier method to cross-check the Carr-Madan
> pricer.
>
> Args:
>     psi: characteristic exponent ``psi(u)`` (as in :func:`levy_price`).
>     n_terms: number of cosine terms.
>     L: truncation-range width in standard deviations (10-12 is ample).
>     cumulants: optional ``(c1, c2, c4)`` of the log-return to set ``[a, b]``;
>         if omitted they are estimated by differencing ``psi`` numerically.
>
> Puts use put-call parity.

### `levy_price(S, K, t, r, q, psi, option_type=<OptionType.CALL: 'call'>, alpha=1.5, upper=200.0) -> float`  _function_

> Price a European call/put for a Levy model via Carr-Madan + parity.

## cev

### `cev_greeks(S, K, t, r, sigma, beta, option_type=<OptionType.CALL: 'call'>, q=0.0)`  _function_

> Greeks of a CEV option by central finite differences.
>
> Differentiates :func:`cev_price` for ``delta`` (dV/dS), ``gamma`` (d2V/dS2),
> ``vega`` (dV/dsigma), and ``theta`` (calendar decay). As ``beta -> 1`` the
> Greeks approach the Black-Scholes Greeks (``sigma`` is calibrated so the
> ATM instantaneous vol matches). Lower ``beta`` steepens the local-vol skew,
> lifting put deltas and gammas in the left wing. Returns a dict with ``price``
> and those fields.

### `cev_price(S, K, t, r, sigma, beta, option_type=<OptionType.CALL: 'call'>, q=0.0)`  _function_

> European CEV option price (Schroder 1989), for 0 <= beta < 1.
>
> Args:
>     sigma: the CEV volatility level, calibrated so that at ``S`` the
>         instantaneous lognormal vol equals ``sigma`` — i.e. the scale
>         ``delta = sigma * S^{1-beta}``. This makes ``beta`` control only the
>         skew, with ``sigma`` comparable to a Black-Scholes vol.
>     beta: elasticity in [0, 1). beta = 1 recovers Black-Scholes (handled by
>         a limit); beta = 0 is the Bachelier-like absolute-diffusion case.
>     q: continuous dividend yield.
>
> Puts are obtained by put-call parity.

### `cev_smile(S, strikes, t, r, sigma, beta, q=0.0)`  _function_

> The Black-Scholes implied-vol smile a CEV model produces.
>
> Prices a European call at each strike under CEV, then inverts each price to
> its Black-Scholes implied volatility, returning ``(log_moneyness, vol)``
> pairs sorted by strike (log-moneyness on the forward ``F = S e^{(r-q) t}``).
> Because ``sigma`` is calibrated to the ATM instantaneous vol, the smile
> passes near ``sigma`` at the money; ``beta < 1`` makes local volatility fall
> as spot rises, producing a downward skew (steeper for smaller ``beta``).

### `noncentral_chisq_cdf(x: float, k: float, lam: float) -> float`  _function_

> Noncentral chi-square CDF at ``x`` with ``k`` dof and noncentrality ``lam``.
>
> A Poisson(lam/2)-weighted sum of central chi-square CDFs:
> ``F(x; k, lam) = sum_j pois(j; lam/2) * P((k+2j)/2, x/2)``. The summation
> starts at the Poisson mode ``j0 = floor(lam/2)`` and expands outward, so it
> stays numerically stable even when ``lam`` is large (the naive j=0 start
> underflows because ``e^{-lam/2}`` is zero to machine precision).

## cgmy

### `cgmy_greeks(S, K, t, r, C, G, M, Y, option_type=<OptionType.CALL: 'call'>, q=0.0, alpha=1.5)`  _function_

> Greeks of a CGMY option by central finite differences.
>
> Central differences of :func:`cgmy_price` for the spot Greeks ``delta``
> (dV/dS), ``gamma`` (d2V/dS2), and ``theta`` (calendar decay), plus the
> tail-activity sensitivity ``d_Y`` (dV/dY). ``d_Y`` uses a one-sided bump if a
> central one would push ``Y`` to or past 2. Returns a dict with ``price``,
> ``delta``, ``gamma``, ``theta``, ``d_Y``.

### `cgmy_price(S, K, t, r, C, G, M, Y, option_type=<OptionType.CALL: 'call'>, q=0.0, alpha=1.5, upper=200.0) -> float`  _function_

> Price a European option under the CGMY model by Carr-Madan inversion.
>
> Args:
>     C, G, M, Y: CGMY parameters (``C > 0``; ``G, M > 0``; ``Y < 2`` and not a
>         non-negative integer -- ``Y = 0`` is Variance Gamma, handled by the
>         limit only approximately here, so pass a small ``Y`` instead).
>     alpha: Carr-Madan damping factor (> 0); the call transform needs
>         ``E[S_T^{alpha+1}] < infinity``, i.e. ``alpha + 1 < M``.
>     upper: Fourier-integral truncation.
>
> ``C = 0`` gives a degenerate (deterministic-forward) payoff. Puts use
> put-call parity.

### `cgmy_smile(S, strikes, t, r, C, G, M, Y, q=0.0, alpha=1.5)`  _function_

> Black-Scholes implied-vol smile the CGMY model produces.
>
> Prices a call at each strike and inverts to a Black-Scholes implied vol,
> returning ``(log_moneyness, vol)`` pairs sorted by strike on the forward
> ``F = S e^{(r-q) t}``. Tail asymmetry (``G != M``) tilts the smile into a
> skew; smaller ``Y`` fattens the wings.

## cheyette

### `cheyette_G(kappa, tau)`  _function_

> The Cheyette/Hull-White G-function ``(1 - e^{-kappa tau}) / kappa``.

### `cheyette_bond_option(P0S, P0T, kappa, sigma, expiry, maturity, strike, is_call=True)`  _function_

> Price a European option on a zero-coupon bond under Cheyette (const sigma).
>
> Option expires at ``expiry`` (= ``t``) on a bond maturing at ``maturity``
> (= ``T``), struck at ``strike``. ``P0S = P(0, expiry)`` and
> ``P0T = P(0, maturity)`` are today's discount factors. This is the exact
> Hull-White bond-option formula (Cheyette with constant sigma coincides with
> Hull-White), used as the analytic anchor for the model.

### `cheyette_bond_option_greeks(P0S, P0T, kappa, sigma, expiry, maturity, strike, is_call=True)`  _function_

> Greeks of a Cheyette (Hull-White) zero-coupon-bond option.
>
> Sensitivities of :func:`bond_option` to the two discount factors -- exact,
> from the Black-style form -- plus the vol sensitivity by finite difference:
>
>   * ``delta_T`` = dV/dP0T = ``N(d1)`` (call), the underlying-bond delta;
>   * ``delta_S`` = dV/dP0S (discount-leg delta);
>   * ``vega`` = dV/dsigma.
>
> Returns a dict with ``price``, ``delta_T``, ``delta_S``, ``vega``.

### `cheyette_cap(discounts, kappa, sigma, strike, notional=1.0)`  _function_

> Cheyette cap: strip of caplets over successive periods.
>
> ``discounts`` is an increasing list of ``(t_i, P(0, t_i))``; each consecutive
> pair is one caplet ``[t_{i-1}, t_i]``. Returns the summed caplet value.

### `cheyette_caplet(P0_reset, P0_pay, kappa, sigma, reset, pay, strike, notional=1.0)`  _function_

> Price a caplet under Cheyette (constant sigma) on ``[reset, pay]``.
>
> A caplet paying ``tau (L - strike)^+`` at ``pay`` (with ``L`` the simple
> forward rate and ``tau = pay - reset``) equals ``notional (1 + strike tau)``
> put options on the zero-coupon bond ``P(reset, pay)`` struck at
> ``1 / (1 + strike tau)`` (the standard caplet<->bond-put identity).

### `cheyette_floor(discounts, kappa, sigma, strike, notional=1.0)`  _function_

> Cheyette floor: strip of floorlets over successive periods (see :func:`cap`).

### `cheyette_floorlet(P0_reset, P0_pay, kappa, sigma, reset, pay, strike, notional=1.0)`  _function_

> Floorlet on ``[reset, pay]`` under Cheyette via the bond-call identity.

### `cheyette_y(kappa, sigma, t)`  _function_

> The auxiliary state ``y(t)`` for constant sigma (= Var[x(t)]).

### `cheyette_zero_bond(P0T, P0t, x, y, kappa, t, T)`  _function_

> Cheyette zero-coupon bond ``P(t, T)`` given the state ``(x, y)``.
>
> ``P0T = P(0, T)`` and ``P0t = P(0, t)`` are today's discount factors.

## chooser

### `chooser_option(S, K, t_choose, T, r, sigma, b=None) -> float`  _function_

> Price a simple chooser option (Rubinstein 1991).
>
> Args:
>     t_choose: time (years) until the call/put choice is made.
>     T: total time (years) to the underlying option's expiry (>= t_choose).
>     b: cost of carry (defaults to r).

### `chooser_option_greeks(S, K, t_choose, T, r, sigma, b=None)`  _function_

> Greeks of a simple chooser option, exact by decomposition.
>
> The chooser is exactly ``C(S, K, T) + P(S, K e^{-b(T - t_choose)}, t_choose)``
> -- a call to ``T`` plus a put struck at the discounted-forward level expiring
> at the choice date. Both legs are Black-Scholes prices in ``S`` and ``sigma``
> (the put's strike does not depend on either), so ``delta``, ``gamma``, and
> ``vega`` are the exact sums of the two legs' BSM Greeks -- no finite
> difference. Returns a dict with ``price``, ``delta``, ``gamma``, ``vega``.

### `complex_chooser_option(S, Kc, Kp, t_choose, Tc, Tp, r, sigma, b=None) -> float`  _function_

> Complex chooser option (Rubinstein 1991), closed form.
>
> At the choice date ``t_choose`` the holder keeps whichever is worth more: a
> call struck at ``Kc`` expiring at ``Tc``, or a put struck at ``Kp`` expiring
> at ``Tp`` (the two legs may differ in both strike and maturity). Rubinstein's
> formula prices this with bivariate normals coupling the choice date to each
> leg's expiry:
>
>     V = S e^{(b-r)Tc} M(d1, y1; rho_c) - Kc e^{-r Tc} M(d2, y1 - sig sqrt Tc; rho_c)
>         - S e^{(b-r)Tp} M(-d1, -y2; rho_p) + Kp e^{-r Tp} M(-d2, -y2 + sig sqrt Tp; rho_p)
>
> where ``d1,d2`` use the critical spot ``I`` (the level where the two legs are
> equal at ``t_choose``), ``y1,y2`` use ``Kc,Kp``, and
> ``rho_c = sqrt(t_choose/Tc)``, ``rho_p = sqrt(t_choose/Tp)``.

### `complex_chooser_option_greeks(S, Kc, Kp, t_choose, Tc, Tp, r, sigma, b=None)`  _function_

> Greeks of a complex chooser by central finite differences of
> :func:`complex_chooser_option`: ``delta`` (dV/dS), ``gamma`` (d2V/dS2),
> ``vega`` (dV/dsigma), ``theta`` (calendar decay -- all maturities and the
> choice date shrink together). Returns a dict with ``price`` and those
> fields.

## cir

### `cir_bond_greeks(r0, t, kappa, theta, sigma)`  _function_

> Rate sensitivities of a CIR zero-coupon bond, exact.
>
> The bond is ``P = A(t) e^{-B(t) r0}``, so its short-rate sensitivities are
> closed form: ``rho_r = dP/dr0 = -B P`` and ``gamma_r = d2P/dr0^2 = B^2 P``.
> The rate ``duration`` is ``-1/P dP/dr0 = B`` and ``convexity`` is
> ``1/P d2P/dr0^2 = B^2``. Returns a dict with ``price``, ``rho_r``,
> ``gamma_r``, ``duration``, ``convexity``.

### `cir_bond_option(r0, t_option, t_bond, strike, kappa, theta, sigma, option_type=<OptionType.CALL: 'call'>)`  _function_

> European option on a CIR zero-coupon bond (CIR 1985, exact).
>
> Option expires at ``t_option`` on a bond maturing at ``t_bond`` (``> t_option``),
> struck at ``strike`` on the bond price. Uses the noncentral chi-square
> formula: with ``g = sqrt(kappa^2 + 2 sigma^2)`` and the CIR affine
> ``A, B`` over ``t_bond - t_option``, the call is
>
>     P(0,t_bond) X2(...; nc1) - strike P(0,t_option) X2(...; nc2),
>
> with critical rate ``r* = ln(A/strike)/B``. Puts follow from put-call parity.

### `cir_cap(r0, dates, strike, kappa, theta, sigma, notional=1.0)`  _function_

> CIR cap: strip of caplets over successive ``dates`` (increasing times).

### `cir_caplet(r0, reset, pay, strike, kappa, theta, sigma, notional=1.0)`  _function_

> Caplet on ``[reset, pay]`` under CIR via the bond-put identity.

### `cir_coupon_bond_option(r0, t_option, cashflows, strike, kappa, theta, sigma, option_type=<OptionType.CALL: 'call'>)`  _function_

> European option on a coupon bond under CIR (Jamshidian decomposition).
>
> ``cashflows`` is ``[(t_i, c_i), ...]`` with ``t_i > t_option``. The CIR bond
> is monotone decreasing in ``r0``, so Jamshidian applies: solve for the
> critical rate ``r*`` where the coupon bond's value at expiry equals
> ``strike``, then sum the ``c_i``-weighted CIR zero-coupon-bond options
> (:func:`cir_bond_option`) struck at ``K_i = P(t_option, t_i | r*)``. Exact.

### `cir_expected_rate(r0, t, kappa, theta, sigma=0.0)`  _function_

> Expected CIR short rate ``E[r_t] = theta + (r0 - theta) e^{-kappa t}``.
>
> The square-root diffusion does not change the mean, so it matches the
> Vasicek/OU mean (``sigma`` is accepted only for a uniform signature).

### `cir_floor(r0, dates, strike, kappa, theta, sigma, notional=1.0)`  _function_

> CIR floor: strip of floorlets over successive ``dates``.

### `cir_floorlet(r0, reset, pay, strike, kappa, theta, sigma, notional=1.0)`  _function_

> Floorlet on ``[reset, pay]`` under CIR via the bond-call identity.

### `cir_rate_variance(r0, t, kappa, theta, sigma)`  _function_

> Variance of the CIR short rate at horizon ``t``.
>
>     Var[r_t] = r0 (sigma^2/kappa)(e^{-kappa t} - e^{-2 kappa t})
>              + theta (sigma^2/(2 kappa))(1 - e^{-kappa t})^2.
>
> Unlike Vasicek the variance depends on ``r0`` (state-dependent diffusion);
> as ``t -> infinity`` it approaches the stationary variance
> ``theta sigma^2/(2 kappa)``.

### `cir_stationary_distribution(kappa, theta, sigma)`  _function_

> Long-run (stationary) distribution of the CIR rate as a Gamma law.
>
> Returns ``(shape, scale, mean, variance)``. As ``t -> infinity`` the rate is
> Gamma with shape ``2 kappa theta / sigma^2`` and scale ``sigma^2/(2 kappa)``,
> hence mean ``theta`` and variance ``theta sigma^2/(2 kappa)``. The Feller
> condition ``2 kappa theta >= sigma^2`` (shape >= 1) keeps the rate strictly
> positive.

### `cir_swaption(r0, expiry, pay_times, fixed_rate, kappa, theta, sigma, payer=True, notional=1.0)`  _function_

> European swaption under CIR via the coupon-bond-option identity (exact).
>
> A payer swaption is a put on the fixed-leg coupon bond struck at the
> notional; a receiver is a call. Priced by :func:`cir_coupon_bond_option`.
> ``pay_times`` are the fixed-leg payment dates (all ``> expiry``); accruals
> are the gaps, the first measured from ``expiry``.

### `cir_zero_coupon_bond(r0, t, kappa, theta, sigma)`  _function_

> CIR zero-coupon bond price P(0, t) for a unit face, given r(0)=r0.
>
> Requires ``r0 >= 0`` and positive parameters.

### `cir_zero_coupon_yield(r0, t, kappa, theta, sigma)`  _function_

> Continuously-compounded yield of the CIR zero-coupon bond to ``t``.

## cms

### `cms_adjustment_greeks(forward, sigma, expiry, tenor, freq=1.0, pay_lag=0.0)`  _function_

> Sensitivities of the standard-model CMS convexity adjustment.
>
> Central finite differences of :func:`cms_adjustment_standard` for
> ``d_forward`` (d(CA)/d forward) and ``d_sigma`` (d(CA)/d sigma). The
> adjustment is monotone increasing in the vol (more convexity), so
> ``d_sigma > 0``. Returns a dict with ``adjustment``, ``d_forward``,
> ``d_sigma``.

### `cms_adjustment_standard(forward, sigma, expiry, tenor, freq=1.0, pay_lag=0.0)`  _function_

> Standard-model (linear-TSR) CMS convexity adjustment.
>
> Args:
>     forward: forward swap rate S0.
>     sigma: lognormal (Black) swap-rate volatility.
>     expiry: fixing time in years.
>     tenor: swap tenor in years.
>     freq: payment frequency of the underlying swap (per year).
>     pay_lag: payment delay in years (0 = natural payment).
>
> Under the linear terminal-swap-rate model the adjustment is
> ``CA = G * Var_A(S_T) = G * S0^2 (e^{sigma^2 T} - 1)`` -- the exact
> lognormal variance, not just its ``sigma^2 T`` leading term. Returns the
> additive adjustment so ``E_pay[S_T] = forward + CA``.

### `cms_rate(forward, sigma, expiry, tenor, freq=1.0, pay_lag=0.0)`  _function_

> Convexity-adjusted expected CMS rate under the standard model.

### `cms_rate_convexity_replication(forward, expiry, tenor, vol_fn, freq=1.0, pay_lag=0.0, width=8.0, n=800)`  _function_

> CMS convexity adjustment by static replication over a swaption strip.
>
> Under the linear-TSR model the adjustment is ``G * E_A[(S_T - S0)^2]``, and
> the second moment of the swap rate is replicated model-free by a strip of
> swaptions (Carr-Madan variance replication):
>
>     E_A[(S_T - S0)^2] = 2 * integral_0^inf swaption(K) dK,
>
> with payer swaptions for ``K >= S0`` and receiver swaptions for ``K < S0``,
> each priced at the smile vol ``vol_fn(K)``. This captures the whole smile,
> and with a flat ``vol_fn`` it reproduces :func:`cms_adjustment_standard`.

## commodity

### `asian_commodity_option(avg_forward, strike, sigma, r, expiry, reset_var_frac=0.3333333333333333, is_call=True)`  _function_

> Average-price (Asian) commodity option, Black on the average forward.
>
> Prices an option on the arithmetic average of a commodity's price over the
> averaging window. The average of lognormals is not lognormal, so the average
> forward's variance is reduced by ``reset_var_frac`` (the continuous-averaging
> limit is ``1/3`` of the terminal variance): effective total variance
> ``v = reset_var_frac * sigma^2 * T``. Then Black on ``avg_forward``:
>
>     call = e^{-r T} [F Phi(d1) - K Phi(d2)],  d1,2 = (ln(F/K) +/- 0.5 v)/sqrt(v)
>
> The variance reduction makes the Asian cheaper than the vanilla on the same
> forward. Put and call satisfy ``C - P = e^{-r T}(F - K)``.

### `asian_commodity_option_mc(forward, strike, sigma, r, expiry, n_avg, n_paths=100000, seed=2024, is_call=True, geometric=False)`  _function_

> Monte Carlo Asian commodity option over a discrete monitoring path.
>
> Simulates a driftless (forward-measure) GBM ``F(t) = forward
> exp(-0.5 sigma^2 t + sigma W_t)`` on ``n_avg`` equally-spaced dates, averages
> (arithmetic by default, geometric if ``geometric=True``), and discounts the
> payoff. Independent reference validating both the arithmetic 1/3-variance
> :func:`asian_commodity_option` and the exact
> :func:`geometric_asian_option`. Deterministic per seed.

### `bachelier_spread_option(f1, f2, strike, sigma1, sigma2, rho, r, expiry, is_call=True)`  _function_

> Bachelier (normal-model) spread option on two forwards.
>
> Models each forward as arithmetic Brownian motion, so the spread ``F1 - F2``
> is normal with volatility
> ``sigma = sqrt(sigma1^2 - 2 rho sigma1 sigma2 + sigma2^2)`` (absolute, price
> units). Unlike the lognormal :func:`kirk_spread_option` this prices spreads
> that are or can go negative -- the norm for crack and location spreads. With
> ``m = F1 - F2 - K`` and ``s = sigma sqrt(T)``:
>
>     call = e^{-r T} [m Phi(m/s) + s phi(m/s)]
>     put  = e^{-r T} [-m Phi(-m/s) + s phi(m/s)]
>
> Put and call satisfy ``C - P = e^{-r T} (F1 - F2 - K)``.

### `carry_roll_yield(r, storage_cost=0.0, convenience_yield=0.0)`  _function_

> Roll yield implied by the cost-of-carry model ``y - r - u`` (= -net carry).
>
> Under constant carry ``F(T) = S e^{(r+u-y)T}`` the roll yield is exactly the
> negative net carry :func:`net_cost_of_carry`, so it equals the convenience
> yield net of financing and storage. Positive precisely in backwardation.

### `commodity_calendar_spread(spot, r, t_near, t_far, storage_cost=0.0, convenience_yield=0.0)`  _function_

> Far-minus-near forward spread under one carry rate.
>
> ``F(t_far) - F(t_near)`` from :func:`commodity_forward`. Positive in contango
> (net carry ``r + u - y > 0``, the far contract richer) and negative in
> backwardation, so its sign matches :func:`net_cost_of_carry`.

### `commodity_forward(spot, r, maturity, storage_cost=0.0, convenience_yield=0.0)`  _function_

> Cost-of-carry forward ``S * exp((r + u - y) * T)``.
>
> ``storage_cost`` (u) and ``convenience_yield`` (y) are continuous proportional
> rates. Storage lifts the forward (a cost of carrying the physical), the
> convenience yield lowers it (a benefit of holding it).

### `commodity_forward_curve(spot, r, maturities, storage_cost=0.0, convenience_yield=0.0)`  _function_

> Forward prices across a list of maturities under one carry rate.
>
> Returns ``[(T, F(T)), ...]`` from :func:`commodity_forward` at each maturity.

### `commodity_swap_rate(forward_quotes, discount_factors)`  _function_

> Fair fixed price of a commodity swap: DF-weighted average of the forwards.
>
> A commodity swap exchanges a fixed price for the floating settlement (the
> forward) on each reset. The fair fixed price zeroing the swap is the
> discount-factor-weighted average of the reset forwards
> ``sum_i DF_i F_i / sum_i DF_i`` -- so paying this fixed against the floating
> forwards has zero present value.

### `commodity_swap_value(forward_quotes, discount_factors, fixed_price, notional=1.0, pay_fixed=True)`  _function_

> Present value of a commodity swap versus a fixed price.
>
> Fixed-payer value ``notional * sum_i DF_i (F_i - fixed_price)`` (receiver is
> the negative). Zero when ``fixed_price`` equals :func:`commodity_swap_rate`.

### `convenience_yield_curve(spot, r, forward_quotes, storage_cost=0.0)`  _function_

> Per-tenor convenience yields implied by a forward strip.
>
> ``forward_quotes`` is ``[(T, F(T)), ...]``. Inverts each quote with
> :func:`implied_convenience_yield` at a common ``storage_cost``, so recomputing
> the forward at each ``(T, y_T)`` reprices the input strip exactly.

### `geometric_asian_option(forward, strike, sigma, r, expiry, n_avg, is_call=True)`  _function_

> Exact geometric-average Asian commodity option (discrete monitoring).
>
> For ``n_avg`` equally-spaced monitoring dates the geometric average of
> lognormals is itself lognormal, giving an exact Black-style price. With
> monitoring at ``t_i = i T / n`` (``i = 1..n``), the geometric-average forward
> and its variance are
>
>     F_G  = forward * exp(-0.5 sigma^2 (T - t_bar)),  t_bar = mean(t_i)
>     v_G  = sigma^2 / n^2 * sum_i (2 i - 1) (n - i + 1) * (T / n)   [= adj var]
>
> Because the geometric mean is below the arithmetic mean (AM-GM), this is a
> lower bound for the arithmetic :func:`asian_commodity_option`. Put and call
> satisfy ``C - P = e^{-r T}(F_G - K)``.

### `implied_convenience_yield(spot, forward, r, maturity, storage_cost=0.0)`  _function_

> Convenience yield implied by a market forward (inverts the carry formula).
>
> Solves ``F = S exp((r + u - y) T)`` for ``y``:
> ``y = r + u - ln(F/S)/T``. Inverse of :func:`commodity_forward`.

### `implied_storage_cost(spot, forward, r, maturity, convenience_yield=0.0)`  _function_

> Storage cost implied by a market forward (inverts the carry formula).
>
> ``u = ln(F/S)/T - r + y``. Inverse of :func:`commodity_forward` in ``u``.

### `is_backwardation(r, storage_cost=0.0, convenience_yield=0.0)`  _function_

> True when the curve is in backwardation (``y > r + u``, forwards below spot).
>
> Backwardation occurs when the convenience yield exceeds the financing-plus-
> storage carry, so the net carry :func:`net_cost_of_carry` is negative and
> forwards fall with maturity.

### `kirk_spread_option(f1, f2, strike, sigma1, sigma2, rho, r, expiry, is_call=True)`  _function_

> Kirk (1995) approximation for a spread option on two forwards.
>
> Prices ``max(F1 - F2 - K, 0)`` (call) or ``max(K - (F1 - F2), 0)`` (put) by
> treating ``F2 + K`` as a single lognormal asset and applying Black with an
> effective spread volatility
>
>     sigma_K = sqrt(sigma1^2 - 2 rho sigma1 sigma2 w + sigma2^2 w^2),
>     w = F2 / (F2 + K).
>
> At ``K = 0`` (``w = 1``) it collapses exactly to the
> :func:`margrabe_exchange_option`. Call and put satisfy
> ``C - P = e^{-r T} (F1 - F2 - K)``.

### `margrabe_exchange_option(f1, f2, sigma1, sigma2, rho, r, expiry)`  _function_

> Margrabe (1978) option to exchange asset 2 for asset 1, on forwards.
>
> Exact closed form for the payoff ``max(F1 - F2, 0)`` (a zero-strike spread
> option). With the spread volatility
> ``sigma = sqrt(sigma1^2 - 2 rho sigma1 sigma2 + sigma2^2)``:
>
>     price = e^{-r T} [F1 Phi(d1) - F2 Phi(d2)]
>     d1 = (ln(F1/F2) + 0.5 sigma^2 T) / (sigma sqrt(T)),  d2 = d1 - sigma sqrt(T)
>
> The building block the :func:`kirk_spread_option` approximation reduces to at
> zero strike.

### `mean_reversion_half_life(kappa)`  _function_

> Half-life of mean reversion ``ln(2) / kappa`` (years).
>
> Time for a shock to log-spot to decay to half its size under the Schwartz OU
> dynamics; falls as the mean-reversion speed ``kappa`` rises.

### `net_cost_of_carry(r, storage_cost=0.0, convenience_yield=0.0)`  _function_

> Net proportional carry rate ``r + u - y`` (the forward's growth rate).

### `roll_yield(near_forward, far_forward, t_near, t_far)`  _function_

> Annualized roll yield between two forwards ``ln(F_near/F_far)/(t_far-t_near)``.
>
> The return earned rolling a long position from the far to the near contract as
> time passes, assuming spot is unchanged. Positive in backwardation (near above
> far) and negative in contango, so its sign is the opposite of the
> :func:`commodity_calendar_spread` sign.

### `schwartz_forward(spot, kappa, alpha_star, sigma, maturity)`  _function_

> Commodity forward under the Schwartz (1997) one-factor model.
>
> With log-spot lognormal, ``F(T) = exp(E[X_T] + 0.5 Var[X_T])`` from
> :func:`schwartz_log_mean` and :func:`schwartz_log_variance`. Equals the spot
> at ``T = 0`` and converges to the risk-neutral long-run forward
> ``exp(alpha_star + sigma^2/(4 kappa))`` as ``T -> inf`` -- the mean-reverting
> alternative to the constant-carry :func:`commodity_forward`.

### `schwartz_futures_volatility(sigma, kappa, maturity)`  _function_

> Instantaneous return volatility of the ``maturity``-future under Schwartz.
>
> Because log-spot mean-reverts, the futures return volatility decays with time
> to maturity: ``sigma_F(T) = sigma e^{-kappa T}``. It equals the spot vol
> ``sigma`` for the front (``T = 0``) and falls for longer maturities -- the
> Samuelson effect (near contracts more volatile than deferred).

### `schwartz_implied_alpha(spot, forward, kappa, sigma, maturity)`  _function_

> Risk-neutral long-run log level implied by a single forward quote.
>
> Inverts :func:`schwartz_forward` for ``alpha_star``:
>
>     alpha_star = [ln F - 0.5 Var[X_T] - e^{-kappa T} ln S] / (1 - e^{-kappa T}).
>
> Requires ``maturity > 0`` (at ``T = 0`` the forward carries no information
> about the long-run level). Round-trips with :func:`schwartz_forward`.

### `schwartz_log_mean(spot, kappa, alpha_star, maturity)`  _function_

> Risk-neutral mean of log-spot under the Schwartz one-factor model.
>
> Log-spot ``X = ln S`` follows a mean-reverting Ornstein-Uhlenbeck process
> ``dX = kappa (alpha_star - X) dt + sigma dW`` under the pricing measure, where
> ``alpha_star`` is the risk-adjusted long-run log level. The conditional mean is
>
>     E[X_T] = e^{-kappa T} ln S + (1 - e^{-kappa T}) alpha_star,
>
> starting at ``ln S`` for ``T = 0`` and relaxing to ``alpha_star`` as
> ``T -> inf``.

### `schwartz_log_variance(sigma, kappa, maturity)`  _function_

> Variance of log-spot under the Schwartz one-factor model.
>
> ``Var[X_T] = sigma^2 (1 - e^{-2 kappa T}) / (2 kappa)`` -- zero at ``T = 0``,
> rising monotonically to the stationary ``sigma^2 / (2 kappa)`` as
> ``T -> inf``.

### `schwartz_option(spot, kappa, alpha_star, sigma, strike, r, expiry, is_call=True)`  _function_

> European spot option under the Schwartz one-factor model.
>
> At expiry the spot is lognormal with mean :func:`schwartz_log_mean` and
> variance :func:`schwartz_log_variance`, so the option is a Black-style price
> off the model forward ``F* = schwartz_forward`` and total variance
> ``v = Var[X_T]``:
>
>     d1 = (ln(F*/K) + 0.5 v) / sqrt(v),  d2 = d1 - sqrt(v)
>     call = e^{-r T} [F* Phi(d1) - K Phi(d2)]
>     put  = e^{-r T} [K Phi(-d2) - F* Phi(-d1)]
>
> Put and call satisfy ``C - P = e^{-r T} (F* - K)``. At zero variance the price
> is the discounted intrinsic on ``F*``.

### `schwartz_smith_forward(chi0, xi0, kappa, mu_xi, sigma_chi, sigma_xi, rho, maturity)`  _function_

> Commodity forward under the Schwartz-Smith (2000) two-factor model.
>
> With log-spot lognormal, ``F(T) = exp(E[ln S_T] + 0.5 Var[ln S_T])`` from
> :func:`schwartz_smith_log_mean` and :func:`schwartz_smith_log_variance`.
> Equals the spot ``exp(chi0 + xi0)`` at ``T = 0``. The two factors let the
> front of the curve move on short-term shocks while the back tracks the drifting
> equilibrium -- richer than the single-factor :func:`schwartz_forward`.

### `schwartz_smith_futures_volatility(kappa, sigma_chi, sigma_xi, rho, maturity)`  _function_

> Instantaneous return volatility of the ``maturity``-future (Schwartz-Smith).
>
> The two-factor analogue of :func:`schwartz_futures_volatility`. A future's log
> return loads fully on the persistent long-term factor and on the decaying
> short-term factor with weight ``e^{-kappa T}``, so
>
>     sigma_F(T) = sqrt(e^{-2 kappa T} sigma_chi^2 + sigma_xi^2
>                       + 2 e^{-kappa T} rho sigma_chi sigma_xi).
>
> Falls from the front (short-term shocks fully felt) toward the long-term floor
> ``sigma_xi`` as ``T -> inf`` -- the Samuelson effect with a non-zero long-end
> asymptote. Reduces to the one-factor ``sigma_chi e^{-kappa T}`` when
> ``sigma_xi = 0``.

### `schwartz_smith_log_mean(chi0, xi0, kappa, mu_xi, maturity)`  _function_

> Risk-neutral mean of log-spot under the Schwartz-Smith two-factor model.
>
> Log-spot decomposes as ``ln S = chi + xi``: a short-term deviation ``chi``
> that mean-reverts to zero at speed ``kappa`` (Ornstein-Uhlenbeck) and a
> long-term equilibrium ``xi`` that drifts as arithmetic Brownian motion with
> risk-neutral rate ``mu_xi``. The conditional mean is
>
>     E[ln S_T] = e^{-kappa T} chi0 + xi0 + mu_xi T,
>
> starting at ``chi0 + xi0`` for ``T = 0``; the short-term contribution decays
> while the long-term level persists and drifts.

### `schwartz_smith_log_variance(kappa, sigma_chi, sigma_xi, rho, maturity)`  _function_

> Variance of log-spot under the Schwartz-Smith two-factor model.
>
> Combines the mean-reverting short-term variance, the linearly-growing
> long-term variance, and their covariance:
>
>     Var[ln S_T] = (1 - e^{-2 kappa T}) sigma_chi^2 / (2 kappa)
>                   + sigma_xi^2 T
>                   + 2 (1 - e^{-kappa T}) / kappa * rho sigma_chi sigma_xi.
>
> Zero at ``T = 0``; as ``T -> inf`` the short-term term saturates and the total
> variance grows like ``sigma_xi^2 T`` (long-term factor dominates).

### `seasonal_forward(spot, r, maturity, seasonal_factor, storage_cost=0.0, convenience_yield=0.0)`  _function_

> Cost-of-carry forward scaled by a multiplicative seasonal factor.
>
> ``seasonal_factor * commodity_forward(...)`` -- lifts or discounts the carry
> forward for the delivery month's seasonal pattern (e.g. gas in winter). With
> a :func:`quantforge.normalize_seasonal_factors` factor the annual average is
> unchanged.

### `spark_spread_option(power_forward, fuel_forward, heat_rate, strike, sigma_power, sigma_fuel, rho, r, expiry, is_call=True, emissions_rate=0.0, carbon_forward=0.0)`  _function_

> Spark-spread (or dark-spread) option on a power generator's margin.
>
> A gas/coal plant's clean generation margin per MWh is
>
>     spread = power - heat_rate * fuel - emissions_rate * carbon,
>
> where ``heat_rate`` converts fuel price into fuel cost per MWh (MMBtu/MWh) and
> the optional emissions term charges carbon. This prices ``max(spread - K, 0)``
> (a call, the option to run the plant) via the normal
> :func:`bachelier_spread_option`, treating the scaled fuel-plus-carbon leg as a
> single lognormal-free asset with volatility ``heat_rate * sigma_fuel`` (the
> carbon leg is taken as a deterministic adder). Because generation margins are
> routinely negative or near zero, the Bachelier model is used rather than Kirk.
> Put and call satisfy ``C - P = e^{-r T}(power - heat_rate*fuel -
> emissions_rate*carbon - K)``.

### `spread_option_mc(f1, f2, strike, sigma1, sigma2, rho, r, expiry, n_paths=100000, seed=12345, is_call=True)`  _function_

> Monte Carlo price of a spread option under bivariate lognormal forwards.
>
> Simulates ``F1 e^{-0.5 sigma1^2 T + sigma1 sqrt(T) Z1}`` and the analogous
> ``F2`` with correlated normals ``corr(Z1, Z2) = rho`` (Cholesky), averaging the
> discounted payoff ``max(F1 - F2 - K, 0)`` (call) or its put. An independent
> reference for the :func:`kirk_spread_option` approximation. Uses a
> deterministic linear-congruential stream so results are reproducible.

### `tolling_value(power_forwards, fuel_forwards, heat_rate, strike, sigma_power, sigma_fuel, rho, r, expiries, discount_factors=None, emissions_rate=0.0, carbon_forwards=None)`  _function_

> Value of a tolling agreement as a strip of daily spark-spread call options.
>
> A tolling deal grants the right (not obligation) to run a plant each delivery
> period, so its value is the sum of :func:`spark_spread_option` calls over the
> periods -- one per ``(power_forward, fuel_forward, expiry)`` triple at a common
> ``heat_rate``, ``strike`` (variable O&M) and vols. ``strike`` here is the
> per-MWh variable cost; each spark option already discounts by ``e^{-r T}``, so
> ``discount_factors`` (if given) rescales that if a separate curve is wanted --
> by default the internal ``e^{-r expiry}`` is used. Increasing with the number
> of run periods.

### `turnbull_wakeman_asian(forward, strike, sigma, r, expiry, n_avg, is_call=True)`  _function_

> Turnbull-Wakeman arithmetic-average Asian option (two-moment match).
>
> The arithmetic average of lognormals is not lognormal, so Turnbull-Wakeman
> (1991) matches its first two moments to a lognormal and applies Black. For
> ``n_avg`` equally-spaced driftless (forward-measure) monitoring dates
> ``t_i = i T / n``:
>
>     M1 = forward
>     M2 = forward^2 / n^2 * sum_i sum_j exp(sigma^2 min(t_i, t_j))
>
> The effective total variance is ``v = ln(M2 / M1^2)`` and the price is Black on
> ``forward`` with variance ``v``. More accurate than the fixed 1/3-variance
> :func:`asian_commodity_option`; at ``n = 1`` (``v = sigma^2 T``) it reduces to
> the vanilla. Put and call satisfy ``C - P = e^{-r T}(forward - K)``.

## compound

### `compound_option(S, K1, K2, t1, t2, r, sigma, kind='call-on-call', b=None) -> float`  _function_

> Price a compound option (Geske 1979).
>
> Args:
>     K1: strike of the compound (paid at ``t1`` to obtain the underlying).
>     K2: strike of the underlying option (expiring at ``t2``).
>     t1: expiry of the compound (the decision date), ``0 < t1 < t2``.
>     t2: expiry of the underlying option.
>     kind: "call-on-call", "call-on-put", "put-on-call", "put-on-put".
>     b: cost of carry (defaults to r).

### `compound_option_greeks(S, K1, K2, t1, t2, r, sigma, kind='call-on-call', b=None)`  _function_

> Greeks of a compound option (Geske) by central finite differences.
>
> Differentiates :func:`compound_option` for ``delta`` (dV/dS), ``gamma``
> (d2V/dS2), ``vega`` (dV/dsigma), and ``theta`` (calendar decay, ``-dV/dt``
> shifting both expiries together). Returns a dict with ``price`` and those
> fields. ``kind`` is one of ``call-on-call``/``call-on-put``/``put-on-call``/
> ``put-on-put``.

## copula

### `cdo_tranche_expected_loss(attachment, detachment, pd, rho, n_steps=2000)`  _function_

> Expected loss of a CDO tranche in the Vasicek large-pool limit.
>
> Integrates the portfolio loss distribution over the tranche
> ``[attachment, detachment]`` and normalizes by the tranche width, giving the
> expected tranche loss as a fraction of the tranche notional. Equity (low
> attachment) tranches lose more than senior tranches at the same correlation.
> Trapezoidal integration of ``E[min(max(L - a, 0), d - a)] / (d - a)`` using
> the survival ``1 - F(l)``.

### `cdo_tranche_expected_loss_mc(attachment, detachment, pd, rho, n_names=100, n_paths=20000, seed=8675309)`  _function_

> Monte Carlo CDO tranche expected loss under the single-factor model.
>
> Simulates a finite pool of ``n_names``: a common factor ``M`` and idiosyncratic
> shocks give each name's asset value ``sqrt(rho) M + sqrt(1 - rho) Z_i``; a name
> defaults when it falls below ``Phi^{-1}(pd)``. Averages the tranche loss over
> the portfolio-loss fraction across paths. An independent finite-pool reference
> for the large-pool :func:`cdo_tranche_expected_loss` (they agree as
> ``n_names -> inf``). Deterministic per seed.

### `clayton_copula(u, v, theta)`  _function_

> Clayton copula ``(u^{-theta} + v^{-theta} - 1)^{-1/theta}`` (``theta > 0``).
>
> Lower-tail dependent (assets crash together); reduces to independence as
> ``theta -> 0``.

### `clayton_lower_tail_dependence(theta)`  _function_

> Lower-tail dependence of the Clayton copula ``2^{-1/theta}``.
>
> In ``(0, 1)`` for ``theta > 0`` -- rising toward 1 as ``theta`` grows (stronger
> joint-crash dependence).

### `clayton_theta_from_tau(tau)`  _function_

> Clayton ``theta`` from Kendall's tau: ``2 tau / (1 - tau)``.
>
> Inverse of ``tau = theta / (theta + 2)``. Requires ``0 <= tau < 1``.

### `first_to_default_probability(pd1, pd2, rho)`  _function_

> Probability that at least one of two names defaults (Gaussian copula).
>
> ``P(A or B) = pd1 + pd2 - C(pd1, pd2; rho)`` by inclusion-exclusion. Lies
> between ``max(pd1, pd2)`` and ``min(pd1 + pd2, 1)``, and falls as correlation
> rises (correlated defaults overlap more, so fewer *distinct* default events).

### `frank_copula(u, v, theta)`  _function_

> Frank copula (``theta != 0``), symmetric with no tail dependence.
>
> ``C(u, v) = -1/theta * ln(1 + (e^{-theta u} - 1)(e^{-theta v} - 1) /
> (e^{-theta} - 1))``. Positive dependence for ``theta > 0``, negative for
> ``theta < 0``; reduces to independence as ``theta -> 0``.

### `gaussian_copula(u, v, rho)`  _function_

> Gaussian copula ``C(u, v) = Phi_rho(Phi^{-1}(u), Phi^{-1}(v))``.
>
> The dependence structure of a bivariate normal with correlation ``rho``. Zero
> tail dependence for ``|rho| < 1``; reduces to ``u v`` at ``rho = 0``.

### `gaussian_copula_joint_default(pd1, pd2, rho)`  _function_

> Joint default probability of two names under the Gaussian copula.
>
> Both default when their latent normals fall below their default thresholds
> ``Phi^{-1}(pd_i)``; the joint probability is the Gaussian copula
> ``C(pd1, pd2; rho)``. Rises above the independent product ``pd1 * pd2`` for
> ``rho > 0`` and equals it at ``rho = 0``.

### `gumbel_copula(u, v, theta)`  _function_

> Gumbel copula ``exp(-((-ln u)^theta + (-ln v)^theta)^{1/theta})`` (``theta >= 1``).
>
> Upper-tail dependent (assets rally together); reduces to independence at
> ``theta = 1``.

### `gumbel_theta_from_tau(tau)`  _function_

> Gumbel ``theta`` from Kendall's tau: ``1 / (1 - tau)``.
>
> Inverse of ``tau = 1 - 1/theta``. Requires ``0 <= tau < 1``.

### `gumbel_upper_tail_dependence(theta)`  _function_

> Upper-tail dependence of the Gumbel copula ``2 - 2^{1/theta}``.
>
> Zero at ``theta = 1`` (independence) rising toward 1 as ``theta -> inf``.

### `vasicek_loss_cdf(loss, pd, rho)`  _function_

> CDF of the large-homogeneous-portfolio loss fraction (Vasicek limit).
>
> In the single-factor Gaussian-copula limit of an infinitely granular pool with
> default probability ``pd`` and asset correlation ``rho``, the fractional loss
> ``L`` has closed-form CDF
>
>     P(L <= x) = Phi( (sqrt(1 - rho) Phi^{-1}(x) - Phi^{-1}(pd)) / sqrt(rho) ).
>
> ``loss`` is a fraction in ``[0, 1]`` (LGD assumed 1). Increasing in ``loss``.

### `vasicek_loss_quantile(q, pd, rho)`  _function_

> Portfolio loss at confidence ``q`` (the Vasicek/Basel capital formula).
>
> Inverse of :func:`vasicek_loss_cdf`:
>
>     L(q) = Phi( (Phi^{-1}(pd) + sqrt(rho) Phi^{-1}(q)) / sqrt(1 - rho) ).
>
> The worst-case loss not exceeded with probability ``q`` -- the basis of the
> Basel IRB capital charge. Increasing in ``q``, ``pd`` and ``rho``.

## correlation

### `correlation_term_structure(weights, member_vol_curves, index_vol_curve, expiries)`  _function_

> Implied correlation at each expiry across a term structure.
>
> Args:
>     weights: index member weights (constant across expiries).
>     member_vol_curves: list per member of that member's vol at each expiry,
>         i.e. member_vol_curves[i][j] is member i's vol at expiries[j].
>     index_vol_curve: the index's implied vol at each expiry.
>     expiries: the tenors (used only as labels in the returned pairs).
>
> Returns a list of ``(expiry, implied_correlation)`` pairs, applying
> :func:`implied_correlation` slice by slice.

### `dispersion_basket_vol(weights: Sequence[float], vols: Sequence[float]) -> float`  _function_

> The zero-correlation ("fully diversified") index vol, sqrt(sum w^2 sig^2).
>
> A useful lower reference: the index vol if the members were uncorrelated.

### `dispersion_trade_pnl(weights, realized_member_vols, realized_index_vol, strike_member_vols, strike_index_vol, variance_notional=1.0)`  _function_

> P&L of a variance dispersion trade (short index var, long member var).
>
> A dispersion trade sells index variance and buys the weighted member
> variances. Its variance P&L per unit notional is
>
>     (sum_i w_i (rv_i^2 - k_i^2))  -  (rv_index^2 - k_index^2),
>
> the long member-variance legs minus the short index-variance leg (strikes
> ``k``). Because index variance carries the correlation, the trade profits when
> realized correlation comes in *below* what was implied (index realizes calmer
> than the members would imply), and is zero when realized matches strikes.

### `ewma_correlation(returns_x, returns_y, lam=0.94)`  _function_

> Exponentially-weighted correlation of two aligned return series.
>
> The EWMA covariance divided by the product of the EWMA volatilities (all on
> the same decay), so it stays in ``[-1, 1]``.

### `ewma_covariance(returns_x, returns_y, lam=0.94)`  _function_

> Exponentially-weighted covariance of two aligned return series.
>
> RiskMetrics-style recursion ``s_t = lam s_{t-1} + (1-lam) x_t y_t`` seeded
> from the first product, giving more weight to recent observations. ``lam``
> is the decay (0.94 for daily data). The series must be equal length and
> zero-mean is assumed (the RiskMetrics convention for returns).

### `implied_correlation(weights: Sequence[float], vols: Sequence[float], index_vol: float) -> float`  _function_

> Common implied correlation consistent with the quoted ``index_vol``.
>
> Returns rho in principle within [-1, 1]; a value outside that band signals
> an index vol inconsistent with the member vols (arbitrage or stale quotes)
> and is returned unclamped so the caller can see it.

### `index_vol_from_correlation(weights: Sequence[float], vols: Sequence[float], rho: float) -> float`  _function_

> Index volatility implied by member weights/vols and a common correlation.

### `realized_beta(asset_returns, market_returns)`  _function_

> Realized beta of an asset to the market: ``Cov(a, m) / Var(m)``.
>
> Ordinary (equal-weight) sample covariance over variance, the slope of a
> regression of asset returns on market returns. Series must be equal length.

## credit

### `SurvivalCurve(times: Sequence[float], hazards: Sequence[float])`  _class_

> Piecewise-constant hazard-rate survival curve.
>
> Built from pillar times and the *forward* hazard rate on each segment
> ``[t_{i-1}, t_i]``. ``survival(t)`` returns ``Q(t) = exp(-integral h)`` and
> ``default_density(t)`` returns ``h(t) Q(t)``.

### `bootstrap_survival_curve(quote_maturities, quote_spreads, r, recovery=0.4, freq=4, n_steps_per_year=100, tol=1e-10, max_iter=100)`  _function_

> Bootstrap a piecewise-constant hazard curve from par CDS quotes.
>
> Given increasing ``quote_maturities`` and their par ``quote_spreads``, solve
> each tenor's forward hazard in turn (holding earlier segments fixed) so that
> the model par spread of :func:`cds_par_spread` reproduces the quote. Uses a
> bisection on the hazard, which is monotone in the par spread. Premium legs
> pay ``freq`` times a year; the protection-leg grid uses
> ``n_steps_per_year`` points per year. Returns the calibrated
> :class:`SurvivalCurve`.

### `cds_accrual_on_default(curve: quantforge.credit.SurvivalCurve, pay_times, r, n_steps=400)`  _function_

> Accrued-premium annuity paid on default between coupon dates.
>
> A protection buyer who defaults mid-period still owes the premium accrued
> since the last coupon. This returns the accrual factor (to be multiplied by
> the spread): ``integral (t - t_prev) DF(t) (-dQ)`` over each coupon interval,
> the default time approximated on a uniform sub-grid. Adding this to the
> :func:`risky_annuity` gives the full premium-leg annuity.

### `cds_greeks(curve: quantforge.credit.SurvivalCurve, spread, pay_times, r, recovery=0.4, n_steps=400, protection_buyer=True, bump=0.0001)`  _function_

> Risk sensitivities of a CDS mark-to-market by finite difference.
>
> Returns a dict with:
>
>   * ``value``       -- the mark-to-market :func:`cds_value`;
>   * ``credit01``    -- value change for a 1bp parallel bump of the hazard
>     curve (credit spread risk);
>   * ``ir01``        -- value change for a 1bp parallel bump of the discount
>     rate ``r``;
>   * ``recovery01``  -- value change for a 1-point (0.01) rise in recovery;
>   * ``risky_annuity`` -- the survival-weighted premium annuity.
>
> A protection buyer gains when spreads widen (``credit01 > 0``) and loses as
> recovery rises. Bumps are one-sided by ``bump`` (hazard/rate) or 0.01
> (recovery).

### `cds_par_spread(curve: quantforge.credit.SurvivalCurve, pay_times, r, recovery=0.4, n_steps=400, accrual_on_default=False)`  _function_

> Fair (par) CDS spread: protection-leg PV divided by the premium annuity.
>
> With ``accrual_on_default=True`` the annuity includes the accrued premium
> paid on a mid-period default, which lowers the par spread slightly.

### `cds_premium_leg(curve: quantforge.credit.SurvivalCurve, spread, pay_times, r, accrual=None, accrual_on_default=False, n_steps=400)`  _function_

> PV of the CDS premium leg at a given ``spread`` (annualized).
>
> With ``accrual_on_default=True`` the accrued premium paid on a mid-period
> default (:func:`cds_accrual_on_default`) is added to the survival-weighted
> coupon annuity, the market-standard convention.

### `cds_protection_leg(curve: quantforge.credit.SurvivalCurve, maturity, r, recovery=0.4, n_steps=400)`  _function_

> PV of the CDS protection (default) leg, ``(1-R) integral DF(t) (-dQ)``.
>
> Numerically integrates the loss payment over ``[0, maturity]`` on a uniform
> grid, paying ``(1 - recovery)`` at the (grid-approximated) default time.

### `cds_value(curve: quantforge.credit.SurvivalCurve, spread, pay_times, r, recovery=0.4, n_steps=400, protection_buyer=True)`  _function_

> Mark-to-market value of a CDS at a contractual ``spread``.
>
> Protection buyer is long the protection leg and short the premium leg:
> ``V = protection - spread * annuity``. Positive when the par spread has
> widened beyond the contractual spread.

### `risky_annuity(curve: quantforge.credit.SurvivalCurve, pay_times, r, accrual=None)`  _function_

> Risky (survival-weighted) annuity ``sum_i tau_i DF(t_i) Q(t_i)``.
>
> ``pay_times`` are the premium payment dates; ``accrual`` is the per-period
> year fractions (defaults to the gaps between pay times, starting from 0).
> Discount by flat rate ``r`` or a supplied ``r(t)`` function.

### `risky_bond_price(curve: quantforge.credit.SurvivalCurve, cashflows, r, recovery=0.4, face=100.0, n_steps=400)`  _function_

> Price a defaultable coupon bond under a hazard-rate survival curve.
>
> Each scheduled cashflow ``(t, amount)`` is received only if the issuer
> survives to ``t``, so its PV is ``amount DF(t) Q(t)``. On default the holder
> recovers ``recovery * face``, modelled as a payment at the (grid-approximated)
> default time over ``[0, last cashflow]``:
>
>     price = sum_i CF_i DF(t_i) Q(t_i)
>           + recovery * face * integral DF(t) (-dQ).
>
> With ``recovery = 0`` and no defaults this collapses to the survival-weighted
> cashflow PV; with a zero hazard it recovers the risk-free bond price. Uses a
> flat rate ``r`` or a supplied discount function.

### `risky_bond_yield_spread(curve: quantforge.credit.SurvivalCurve, cashflows, r, recovery=0.4, face=100.0, n_steps=400, tol=1e-10, max_iter=100)`  _function_

> Constant credit spread ``s`` over ``r`` that reproduces the risky price.
>
> Prices the bond with :func:`risky_bond_price`, then finds the flat spread
> such that discounting the *promised* cashflows at ``r + s`` (no explicit
> default/recovery) gives the same value -- the bond's z-spread-like quote.
> Solved by bisection (price is monotone decreasing in the spread).

## daycount

### `day_count(start, end) -> int`  _function_

> Actual number of days between two ``(y, m, d)`` dates.

### `year_fraction(start, end, convention='act/365') -> float`  _function_

> Year fraction between ``start`` and ``end`` under ``convention``.
>
> Both dates are ``(year, month, day)`` tuples with ``end >= start``. Negative
> intervals raise.

## density

### `density_total_mass(strikes: Sequence[float], calls: Sequence[float], t: float, r: float) -> float`  _function_

> Integrate the extracted density; should be close to 1 for a good curve.

### `price_from_density(strikes: Sequence[float], calls: Sequence[float], t: float, r: float, payoff: Callable[[float], float]) -> float`  _function_

> Price a European payoff by integrating it against the extracted density.
>
> ``payoff`` maps a terminal underlying value to its cash payoff. The integral
> uses the trapezoidal rule on the recovered density grid; the result is
> discounted at ``r``.

### `risk_neutral_cdf(strikes: Sequence[float], calls: Sequence[float], t: float, r: float) -> Tuple[List[float], List[float]]`  _function_

> Estimate the risk-neutral CDF F(K) = 1 + e^{rt} dC/dK at midpoints.
>
> Uses central first differences; returns ``(mid_strikes, cdf_values)``.

### `risk_neutral_density(strikes: Sequence[float], calls: Sequence[float], t: float, r: float) -> Tuple[List[float], List[float]]`  _function_

> Estimate the risk-neutral pdf at the interior strikes.
>
> Uses a non-uniform central second difference of the call curve, so strikes
> need not be equally spaced. Returns ``(mid_strikes, densities)`` for the
> interior points (the two endpoints have no central second difference).

## density_metrics

### `density_entropy(S0, t, r, vol_fn, q=0.0, n=600, width=8.0)`  _function_

> Differential entropy ``-integral g ln g dK`` of the terminal-spot density.

### `expected_shortfall(S0, t, r, vol_fn, level, lower=True, q=0.0, n=600, width=8.0)`  _function_

> Risk-neutral tail mean ``E^Q[S_T | S_T < level]`` (or ``> level``).
>
> Returns the conditional expectation of the terminal spot in the tail beyond
> ``level``; ``nan`` if that tail has zero probability.

### `kl_divergence_smiles(S0, t, r, vol_fn_p, vol_fn_q, q=0.0, n=600, width=8.0)`  _function_

> Kullback-Leibler divergence ``KL(g_p || g_q)`` of two smile densities.
>
> Both risk-neutral densities are built on the same strike grid (from
> ``vol_fn_p`` and ``vol_fn_q``), renormalised to unit mass, and
>
>     KL = integral g_p(K) ln( g_p(K) / g_q(K) ) dK
>
> is integrated by the trapezoidal rule. Zero iff the two densities coincide,
> always non-negative, and asymmetric in its arguments. Useful for measuring
> how far one implied distribution sits from another (two dates, two models,
> or implied vs a reference).

### `tail_probability(S0, t, r, vol_fn, level, lower=True, q=0.0, n=600, width=8.0)`  _function_

> Risk-neutral tail probability ``Q(S_T < level)`` (or ``> level``).
>
> Equals the undiscounted price of a cash-or-nothing binary struck at
> ``level``. ``lower=True`` returns the downside probability.

### `wasserstein_smiles(S0, t, r, vol_fn_p, vol_fn_q, q=0.0, n=600, width=8.0)`  _function_

> Wasserstein-1 distance between two smile-implied densities.
>
> For one-dimensional distributions the 1-Wasserstein (earth-mover) distance
> equals the L1 gap between their CDFs,
>
>     W1 = integral |F_p(K) - F_q(K)| dK,
>
> computed here on a shared strike grid from the two Breeden-Litzenberger
> densities. Unlike :func:`kl_divergence_smiles` it is a true metric (symmetric,
> satisfies the triangle inequality) and is measured in price units, so it is a
> robust "how far apart are these distributions" number even when their
> supports differ.

## density_var

### `density_var_es(S0, t, r, vol_fn, pnl, confidence=0.99, q=0.0, n=800, width=8.0)`  _function_

> Risk-neutral VaR and Expected Shortfall of a payoff ``pnl(S_T)``.
>
> Args:
>     pnl: horizon profit-and-loss as a function of the terminal spot (losses
>         negative).
>     confidence: e.g. 0.99 for the 99% level.
>
> Returns ``(var, es)`` with both as positive loss numbers: ``var`` is the loss
> the P&L does not exceed with probability ``confidence`` under the
> risk-neutral density, and ``es`` the mean loss beyond it.

## discount_curve

### `DiscountCurve(times, dfs)`  _class_

> Log-linear discount curve from pillar ``(T, DF)`` points.

### `bootstrap_from_swaps(swap_maturities, par_rates, freq=1.0)`  _function_

> Bootstrap a :class:`DiscountCurve` from par swap rates.
>
> Args:
>     swap_maturities: increasing swap tenors in years (each an integer number
>         of ``1/freq``-year periods).
>     par_rates: the fair fixed rate for each tenor.
>     freq: fixed-leg payments per year (1 = annual).
>
> Solves pillar by pillar: with all shorter discount factors known, each new
> par-rate equation is linear in the final ``DF(T_n)``. Returns the curve whose
> par-swap rates reproduce the inputs.

## displaced

### `displaced_diffusion_greeks(S, K, t, r, sigma, shift=0.0, option_type=<OptionType.CALL: 'call'>, b=None)`  _function_

> Greeks of a displaced-diffusion option by central finite differences.
>
> Differentiates :func:`displaced_diffusion_price` for ``delta`` (dV/dS),
> ``gamma`` (d2V/dS2), ``vega`` (dV/dsigma), and ``theta`` (calendar decay).
> At ``shift = 0`` these reduce to the vanilla Black-Scholes Greeks; a positive
> shift flattens the smile toward normal-model behaviour. Returns a dict with
> ``price`` and those fields.

### `displaced_diffusion_price(S, K, t, r, sigma, shift=0.0, option_type=<OptionType.CALL: 'call'>, b=None) -> float`  _function_

> Price a European option under the displaced-diffusion model.
>
> Args:
>     shift: the displacement added to spot and strike. ``shift = 0`` is
>         Black-Scholes; larger positive shifts push toward normal-model
>         behavior and allow the underlying to fall below zero (down to
>         ``-shift``).
>     b: cost of carry (defaults to r).
>
> The payoff is unchanged (``max(S_T - K, 0)`` etc.); only the diffusion is
> displaced, so the price equals a BSM price on ``S + shift`` / ``K + shift``.

### `displaced_diffusion_smile(S, strikes, t, r, sigma, shift=0.0, b=None)`  _function_

> The Black-Scholes implied-vol smile a displaced-diffusion model produces.
>
> Prices a European call at each strike under the displaced diffusion, then
> inverts each price to its Black-Scholes implied volatility, returning
> ``(log_moneyness, vol)`` pairs sorted by strike (log-moneyness on the forward
> ``F = S e^{b t}``). Because ``sigma`` is calibrated to the ATM instantaneous
> vol, the smile passes near ``sigma`` at the money; a positive ``shift`` makes
> the process partly normal, producing a downward skew (steeper for larger
> shift), while ``shift = 0`` returns a flat Black-Scholes smile.

### `displaced_implied_shift(S, t, r, quotes, b=None, shift_lo=None, shift_hi=None)`  _function_

> Calibrate the displacement that reproduces an observed vol skew.
>
> Displaced diffusion has one skew knob, the ``shift``: a positive shift
> lowers the low-strike wing relative to the high-strike wing (a downward
> skew), while ``shift = 0`` is flat Black-Scholes. Given a set of Black-Scholes
> implied-vol quotes ``quotes = [(K, iv), ...]`` this finds the single shift
> whose displaced-diffusion smile best fits them.
>
> The at-the-money volatility is not a free skew knob here, but it is not
> fixed blindly either: for every trial shift the model's local ``sigma`` is
> re-solved so the displaced smile reproduces the ATM quote (the one whose
> strike is closest to the forward ``F = S e^{b t}``) exactly. That decouples
> level from skew, so the shift is driven purely by the off-ATM quotes and the
> fit is not biased by the local-vs-implied vol convention.
>
> Returns ``(shift, sigma_atm, rmse)`` where ``sigma_atm`` is the local vol at
> the fitted shift and ``rmse`` is the root-mean-square implied-vol error
> across the quotes. Minimises the squared vol error over the shift by
> golden-section search on ``[shift_lo, shift_hi]`` (defaults scale with spot:
> ``[-0.9 S, 20 S]``, staying above the ``-shift`` floor).

## double_heston

### `double_heston_greeks(S, K, t, r, v01, kappa1, theta1, xi1, rho1, v02, kappa2, theta2, xi2, rho2, option_type=<OptionType.CALL: 'call'>, q=0.0)`  _function_

> Greeks of a double-Heston option by central finite differences.
>
> Central differences of :func:`double_heston_price` for the spot Greeks
> ``delta`` (dV/dS) and ``gamma`` (d2V/dS2), plus a per-factor
> initial-variance sensitivity ``vega_v01`` and ``vega_v02`` (dV/dv0 for each
> variance factor -- the stochastic-vol analogue of vega). Returns a dict with
> ``price``, ``delta``, ``gamma``, ``vega_v01``, ``vega_v02``.

### `double_heston_price(S, K, t, r, v01, kappa1, theta1, xi1, rho1, v02, kappa2, theta2, xi2, rho2, option_type=<OptionType.CALL: 'call'>, q=0.0, upper=200.0) -> float`  _function_

> Price a European option under the double-Heston model.
>
> Factor 1 is ``(v01, kappa1, theta1, xi1, rho1)`` and factor 2
> ``(v02, kappa2, theta2, xi2, rho2)`` -- typically a fast- and a slow-reverting
> variance. ``q`` is the dividend yield. Zeroing the second factor's ``xi2``
> and ``v02`` recovers single-factor Heston. Puts use put-call parity.

### `double_heston_smile(S, strikes, t, r, v01, kappa1, theta1, xi1, rho1, v02, kappa2, theta2, xi2, rho2, q=0.0)`  _function_

> Black-Scholes implied-vol smile the double-Heston model produces.
>
> Returns ``(log_moneyness, vol)`` pairs sorted by strike on the forward
> ``F = S e^{(r-q) t}``. The two mean-reversion speeds let the short- and
> long-dated skew move more independently than single-factor Heston allows.

## double_heston_calib

### `calibrate_double_heston(S, r, quotes, q=0.0, initial=None, max_iter=6000)`  _function_

> Fit the ten double-Heston parameters to an implied-vol surface.
>
> Args:
>     quotes: iterable of ``(expiry, strike, market_vol)`` Black implied-vol
>         points.
>     initial: optional starting 10-tuple ``(v01, kappa1, theta1, xi1, rho1,
>         v02, kappa2, theta2, xi2, rho2)``; a sensible two-scale seed is used
>         otherwise (a fast- and a slow-reverting factor).
>
> Returns ``(params, rmse)`` -- the fitted 10-tuple and the root-mean-square
> implied-vol error over the quotes.

## dualcurve

### `dual_calibrate_basis(ois_curve, proj_curve, swap_maturities, par_rates, freq=1.0)`  _function_

> Solve the constant basis spread that reprices the given par swap rates.
>
> Minimises the squared par-rate error over a single additive basis on the
> projected forwards, with the built-in Nelder-Mead. Returns ``(basis, rmse)``.

### `dual_float_leg_value(ois_curve, proj_curve, pay_times, basis=0.0)`  _function_

> Value of a unit-notional float leg: projected forwards, OIS-discounted.
>
> ``basis`` is a constant spread (in rate units) added to each projected
> forward -- the tenor/currency basis.

### `dual_forward_rate(proj_curve, T0, T1)`  _function_

> Simply-compounded forward rate over ``[T0, T1]`` off the projection curve.

### `dual_par_swap_rate(ois_curve, proj_curve, pay_times, basis=0.0)`  _function_

> Dual-curve par (fair fixed) rate: float-leg value over the OIS annuity.

### `dual_swap_dv01(ois_curve, proj_curve, pay_times, fixed_rate, payer=True, basis=0.0, bump=0.0001)`  _function_

> Risk of a dual-curve swap by finite differences.
>
> Returns a dict with:
>
>   * ``pv01`` = the fixed-leg annuity ``sum tau_i P_ois(T_i)`` -- the exact
>     ``|dV/d(fixed_rate)|``; a payer's ``dV/d(fixed_rate)`` is ``-pv01``;
>   * ``dv01`` = the value change for a 1bp parallel *drop* in **both** curves
>     (OIS and projection shifted together), the total delta risk;
>   * ``ois_dv01`` / ``proj_dv01`` = the same 1bp-drop risk from shifting only
>     the discount (OIS) or only the projection curve.
>
> ``bump`` is the parallel shift (default 1bp). All figures are per unit
> notional.

### `dual_swap_value(ois_curve, proj_curve, pay_times, fixed_rate, payer=True, basis=0.0)`  _function_

> Value of a unit-notional swap (payer = pay fixed, receive float).

## dv01

### `KeyRateDV01(buckets: Dict[float, float], parallel: float, total_bucketed: float) -> None`  _class_

> KeyRateDV01(buckets: Dict[float, float], parallel: float, total_bucketed: float)

### `key_rate_dv01(price_fn: Callable[[Dict[float, float]], float], base_curve: Dict[float, float], bump: float = 0.0001, one_sided: bool = False) -> quantforge.dv01.KeyRateDV01`  _function_

> Compute key-rate DV01s of a book about ``base_curve``.
>
> Args:
>     price_fn: reprices the book given a ``{tenor: zero_rate}`` curve.
>     base_curve: the current zero curve.
>     bump: the rate shift per bucket (default 1bp). DV01 is reported as the
>         PV change for a +1bp move, scaled from the actual bump.
>     one_sided: use a forward difference instead of the default central one
>         (cheaper, slightly less accurate).
>
> Returns a :class:`KeyRateDV01`. By convention DV01 is negative for a long
> bond-like position (rates up -> PV down).

## equity_comp

### `conversion_premium(convertible_price, S, conversion_ratio)`  _function_

> Conversion (equity) premium: how far the convertible trades above parity.
>
> ``convertible_price / conversion_value - 1`` -- the fractional premium an
> investor pays over the value of the underlying shares for the bond's downside
> protection. Non-negative when the convertible trades at or above parity.

### `conversion_value(S, conversion_ratio)`  _function_

> Parity (conversion) value of a convertible: ``conversion_ratio * S``.
>
> The worth of the shares the bond converts into -- the equity floor of the
> convertible.

### `convertible_bond_value(S, conversion_ratio, face, coupon_rate, maturity, r, sigma, credit_spread=0.0, freq=2)`  _function_

> Convertible bond value via the component (bond floor + call) approximation.
>
> Values the convertible as its :func:`straight_bond_floor` plus a call option on
> the ``conversion_ratio`` shares struck at the floor's per-share equivalent --
> the standard decomposition ``CB = bond floor + conversion_ratio *
> call(S, K=face/ratio adjusted)``. Here the call strike is set so that at
> maturity the holder converts when ``conversion_ratio * S > face``, i.e. strike
> ``= face / conversion_ratio`` on ``conversion_ratio`` shares. The result is at
> least the bond floor and at least the conversion value.

### `convertible_breakeven_years(convertible_price, S, conversion_ratio, bond_coupon_income, dividend_income)`  _function_

> Years for extra income to recoup the conversion premium.
>
> ``(convertible_price - conversion_value) / (bond_coupon_income -
> dividend_income)`` -- the time for the convertible's income advantage over the
> equivalent shares to pay back the dollar conversion premium. Requires the bond
> to yield more than the shares (positive net income); returns ``inf`` if not.

### `dilution_factor(existing_shares, new_shares)`  _function_

> Dilution multiplier ``M / (M + N)`` for issuing ``N`` new shares on ``M``.
>
> Below one whenever new shares are created; one when ``N = 0``. Applied to a
> warrant's per-option payoff because exercise expands the share count.

### `discrete_dividend_price(S, K, t, r, sigma, dividends, is_call=True)`  _function_

> European option on a stock paying known discrete cash dividends.
>
> Uses the escrowed-dividend (spot-minus-PV-of-dividends) approximation: the
> option is priced with Black-Scholes on the dividend-adjusted spot
> ``S - PV(dividends up to expiry)`` and carry ``b = r`` (the adjusted spot grows
> at the risk-free rate). Reduces to the plain BSM call/put when there are no
> dividends before expiry. Put and call satisfy
> ``C - P = (S - PV_div) - K e^{-r t}``.

### `eso_expected_life(vesting, contractual_term, exit_rate)`  _function_

> Expected life of an ESO given post-vest exit and the contractual term.
>
> After vesting, holders leave (and exercise or forfeit) at a constant hazard
> ``exit_rate``; the expected time to exercise, capped at the contractual term,
> is
>
>     vesting + (1 - e^{-exit_rate * (T - vesting)}) / exit_rate
>
> -- between ``vesting`` (immediate exit) and ``contractual_term`` (no exit).

### `eso_value(S, K, contractual_term, r, sigma, vesting, exit_rate, forfeiture_rate=0.0, b=None)`  _function_

> Employee stock option value (Hull-White practical / FASB 123R style).
>
> Prices the ESO as a call at the :func:`eso_expected_life` rather than the full
> term, then multiplies by the probability of surviving pre-vest forfeiture
> ``e^{-forfeiture_rate * vesting}``. Cheaper than the vanilla call on the
> contractual term because early exercise shortens the option and forfeiture
> can extinguish it. With ``exit_rate = 0`` and ``forfeiture_rate = 0`` it
> reduces to the vanilla call.

### `forward_with_dividends(S, t, r, dividends)`  _function_

> Forward price of a stock paying discrete dividends: ``(S - PV_div) e^{r t}``.
>
> The escrowed-dividend forward: the dividend-stripped spot compounded at the
> risk-free rate. Equals ``S e^{r t}`` when no dividends fall before ``t``.

### `investment_premium(convertible_price, bond_floor)`  _function_

> Investment premium: how far the convertible trades above its bond floor.
>
> ``convertible_price / bond_floor - 1`` -- the fractional premium over the
> straight-debt value, paid for the equity upside. Non-negative when the
> convertible trades at or above its floor.

### `pv_dividends(dividends, r)`  _function_

> Present value of a discrete dividend schedule ``[(t, amount), ...]``.
>
> Each cash dividend is discounted at the continuously-compounded rate ``r``:
> ``sum_i D_i e^{-r t_i}``. Dividends at or before time zero (``t <= 0``) are
> excluded (already paid).

### `straight_bond_floor(face, coupon_rate, maturity, r, credit_spread=0.0, freq=2)`  _function_

> Investment (bond) value of a convertible ignoring the conversion option.
>
> Discounts the straight bond's coupons and principal at the risk-free rate plus
> a ``credit_spread`` (continuously compounded). This is the debt floor: the
> convertible cannot be worth less than this if held to maturity without
> converting.

### `warrant_price(S, K, t, r, sigma, existing_shares, new_shares, b=None)`  _function_

> Warrant value: the vanilla call scaled by the :func:`dilution_factor`.
>
> A single-period dilution adjustment -- the standard textbook approximation
> ``value = M/(M+N) * call(S, K, ...)`` -- which reduces to the plain call when
> no new shares are issued. ``S`` is the current (pre-dilution) share price.

## equity_swap

### `dividend_swap_fair_strike(expected_dividends, discount_factors)`  _function_

> Fair strike of a dividend swap: PV of the expected dividend stream.
>
> ``sum_i DF_i * D_i`` -- the present value of the expected dividends the
> floating leg will pay, which the fixed strike must match for the swap to have
> zero value at inception.

### `dividend_swap_value(realized_dividends, strike, discount_factors, notional=1.0)`  _function_

> Value of a dividend swap to the fixed-strike payer (dividend receiver).
>
> ``notional * (sum_i DF_i * D_i - strike)`` -- the PV of realized dividends less
> the fixed strike. Zero when the strike equals the
> :func:`dividend_swap_fair_strike` for the realized stream.

### `financing_leg(notional, funding_rate, spread, year_fraction)`  _function_

> Financing leg of a TRS: ``notional * (funding_rate + spread) * tau``.
>
> The interest the total-return receiver pays on the notional over the accrual
> period ``year_fraction``.

### `total_return_leg(notional, start_price, end_price, dividends)`  _function_

> Equity total-return leg: price return plus dividends on the notional.
>
> ``notional * ((end_price - start_price + dividends) / start_price)`` -- the
> cash the total-return receiver collects (negative if the equity fell more than
> its dividends).

### `total_return_swap_value(notional, start_price, end_price, dividends, funding_rate, spread, year_fraction)`  _function_

> Net value to the total-return receiver: equity leg minus financing leg.
>
> Positive when the equity total return beats the financing cost.

### `trs_fair_spread(start_price, expected_end_price, expected_dividends, funding_rate, year_fraction)`  _function_

> Financing spread that zeroes the expected TRS value.
>
> Solves ``E[equity return] = (funding_rate + spread) * tau`` for the spread:
>
>     spread = E[total return] / tau - funding_rate,
>
> with ``E[total return] = (E[end] - start + E[div]) / start``. The spread the
> financing leg must carry so the swap is fair at inception.

### `variance_swap_mtm(accrued_variance, expected_future_variance, elapsed, total_time, strike_vol, variance_notional, discount_factor)`  _function_

> Mark-to-market of a seasoned variance swap.
>
> Blends the realized (accrued) variance over the elapsed fraction with the
> expected future variance over the remainder to get the expected terminal
> realized variance, then discounts the variance-swap payoff:
>
>     E[realized_var] = (elapsed * accrued + (total - elapsed) * future) / total
>     MTM = DF * variance_notional * (E[realized_var] - strike_vol^2).
>
> At inception (``elapsed = 0``) it is the discounted expected-minus-strike
> variance; at expiry it is the fully-realized payoff.

### `variance_swap_payoff(realized_vol, strike_vol, variance_notional)`  _function_

> Variance-swap payoff ``variance_notional * (realized_vol^2 - strike_vol^2)``.
>
> Settles on the difference between realized and strike *variance* (vols entered
> as decimals, e.g. 0.20 for 20%). Convex in realized vol -- the variance
> convention penalizes large moves more than a volatility swap.

### `vega_notional_to_variance_notional(vega_notional, strike_vol)`  _function_

> Convert a vega notional to the equivalent variance notional.
>
> ``variance_notional = vega_notional / (2 * strike_vol)`` -- the market quotes
> variance swaps in vega terms (P&L per vol point at the strike); this is the
> variance notional that reproduces that sensitivity.

### `volatility_swap_payoff(realized_vol, strike_vol, vega_notional)`  _function_

> Volatility-swap payoff ``vega_notional * (realized_vol - strike_vol)``.
>
> Linear in realized vol (in vol points). Unlike the variance swap it has no
> convexity, so it prices below a variance swap struck at the same vol (the
> convexity value / vol-of-vol adjustment).

## execution

### `cost_variance(trajectory, horizon, sigma)`  _function_

> Timing-risk variance ``sigma^2 tau * sum_k x_k^2`` of the holdings path.
>
> The variance of execution cost from price moves while shares are still held;
> ``x_k`` are the holdings during each interval (using the end-of-interval
> holdings ``x_1..x_N``). Falls as liquidation is front-loaded.

### `efficient_frontier_point(total_shares, n_intervals, horizon, lam, sigma, eta, gamma)`  _function_

> One ``(expected_cost, variance)`` point for a given risk aversion ``lam``.
>
> Builds the :func:`execution_trajectory` at ``lam`` and returns its
> :func:`expected_cost` and :func:`cost_variance` -- sweeping ``lam`` traces the
> Almgren-Chriss efficient frontier (cost rises as variance falls).

### `execution_trades(trajectory)`  _function_

> Per-interval trade sizes ``n_k = x_{k-1} - x_k`` from a holdings trajectory.
>
> Positive sells that sum to the initial holdings.

### `execution_trajectory(total_shares, n_intervals, horizon, lam, sigma, eta)`  _function_

> Optimal holdings trajectory ``[x_0, x_1, ..., x_N]`` (Almgren-Chriss).
>
> ``x_k`` is the shares still held after interval ``k``. Starts at
> ``total_shares`` and ends at zero. For ``lam = 0`` the schedule is linear
> (TWAP); for ``lam > 0`` it is the ``sinh`` profile that liquidates faster
> early. ``eta`` is the temporary-impact coefficient.

### `expected_cost(trajectory, horizon, gamma, eta)`  _function_

> Expected implementation-shortfall cost of a trajectory.
>
> Permanent impact contributes ``0.5 gamma X^2`` (independent of the path); the
> temporary impact contributes ``eta / tau * sum_k n_k^2`` for trades ``n_k``
> over intervals of length ``tau``. Returned in cash (price * shares) units.

### `implementation_shortfall(trajectory, horizon, gamma, eta, sigma)`  _function_

> Decompose expected implementation shortfall into its components.
>
> Returns ``(permanent, temporary, timing_std, total_expected)``: the permanent
> impact ``0.5 gamma X^2``, the temporary impact from :func:`expected_cost` net of
> the permanent part, the timing-risk standard deviation
> ``sqrt(cost_variance)``, and the total *expected* cost (permanent +
> temporary). The expected cost excludes timing risk (mean-zero); the timing std
> is reported separately for the risk budget.

### `kyle_impact(order_size, sigma, daily_volume, liquidity_constant=1.0)`  _function_

> Linear (Kyle) price impact of an order: ``kyle_lambda * order_size``.
>
> Proportional to order size -- doubling the order doubles the impact.

### `kyle_lambda(sigma, daily_volume, liquidity_constant=1.0)`  _function_

> Kyle's lambda: linear price impact per unit of signed order flow.
>
> ``lambda = liquidity_constant * sigma / daily_volume`` -- the slope of price in
> net order flow in Kyle's model, so trading ``q`` shares moves the price by
> ``lambda * q``. Rises with volatility and falls with liquidity (volume).

### `pov_schedule(volume_profile, participation_rate, total_shares=None)`  _function_

> Percentage-of-volume schedule: trade a fixed fraction of each interval.
>
> Each interval's child order is ``participation_rate * volume_profile[i]``. If
> ``total_shares`` is given, trading stops once cumulative fills reach it (the
> last slice is truncated). Returns the per-interval sizes.

### `square_root_impact(order_size, sigma, daily_volume, coefficient=1.0)`  _function_

> Square-root market-impact law ``coefficient * sigma * sqrt(Q / ADV)``.
>
> The empirically-observed concave impact: cost per share grows with the square
> root of participation ``order_size / daily_volume``, so total impact scales
> like ``sqrt(order_size)`` rather than linearly. ``order_size`` is taken as a
> magnitude (absolute value used).

### `twap_schedule(total_shares, n_intervals)`  _function_

> Time-weighted average price schedule: equal-size slices.
>
> Splits ``total_shares`` into ``n_intervals`` equal child orders. The
> benchmark-neutral schedule when volume is uniform.

### `vwap_schedule(total_shares, volume_profile)`  _function_

> Volume-weighted average price schedule: slices proportional to volume.
>
> Allocates ``total_shares`` across intervals in proportion to the expected
> ``volume_profile`` (per-interval volumes), so the fill tracks the day's volume
> curve. With a flat profile it reduces to :func:`twap_schedule`.

## exotics

### `Barrier(*values)`  _class_

> str(object='') -> str
> str(bytes_or_buffer[, encoding[, errors]]) -> str
>
> Create a new string object from the given object. If encoding or
> errors is specified, then the object must expose a data buffer
> that will be decoded using the given encoding and error handler.
> Otherwise, returns the result of object.__str__() (if defined)
> or repr(object).
> encoding defaults to 'utf-8'.
> errors defaults to 'strict'.

### `arithmetic_asian(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None)`  _function_

> Continuously-monitored arithmetic-average-price Asian (Turnbull-Wakeman).
>
> The arithmetic average of a lognormal is not lognormal, so there is no exact
> closed form. Turnbull-Wakeman (1991) matches the first two moments of the
> average to a lognormal and prices with a Black-Scholes-style formula on the
> average's forward. Fast and accurate for typical vols; agrees with the
> arithmetic-Asian Monte Carlo (:func:`quantforge.arithmetic_asian_mc`) to a
> few basis points. Averaging runs over the full life ``[0, t]``.

### `asian_greeks(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None, average='geometric')`  _function_

> Greeks of an Asian option by central finite differences.
>
> ``average`` selects the closed form to differentiate: "geometric"
> (Kemna-Vorst, exact) or "arithmetic" (Turnbull-Wakeman moment match).
> Returns a dict with delta, gamma, vega, and theta (calendar, per year).

### `asset_or_nothing(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None)`  _function_

> Pays the asset value S_T if in the money, else 0.

### `average_strike_arithmetic_asian(S, t, r, sigma, n_fixings=None, fixing_times=None, option_type=<OptionType.CALL: 'call'>, b=None)`  _function_

> Average-strike (floating-strike) discrete arithmetic Asian option.
>
> The strike is the realized arithmetic average: a call pays
> ``max(S_T - A, 0)`` and a put ``max(A - S_T, 0)``, where
> ``A = (1/n) sum_i S_{t_i}``. This is an exchange option between the terminal
> price ``S_T`` (exactly lognormal) and the average ``A`` (matched to a
> lognormal by its first two moments, Levy 1992). The cross-moment
> ``E[S_T A]`` is exact,
>
>     E[S_T A] = (S^2/n) sum_i exp(b (t + t_i) + sigma^2 min(t, t_i))
>              = (S^2/n) sum_i exp(b (t + t_i) + sigma^2 t_i),
>
> so the log-space covariance is ``rho_log = log(E[S_T A]/(E[S_T] E[A]))`` and
> the spread variance is ``Var = sigma^2 t + V_A - 2 rho_log`` with
> ``V_A = log(M2/M1^2)``. A call is then
> ``e^{-rt} (E[S_T] N(d1) - E[A] N(d2))`` in the usual Margrabe form.
>
> Provide either ``n_fixings`` (equally-spaced dates ``t*i/n``, last at
> expiry) or an explicit ``fixing_times`` sequence in ``(0, t]``.

### `average_strike_arithmetic_asian_greeks(S, t, r, sigma, n_fixings=None, fixing_times=None, option_type=<OptionType.CALL: 'call'>, b=None)`  _function_

> Greeks of an average-strike arithmetic Asian by central finite
> differences of :func:`average_strike_arithmetic_asian`: ``delta`` (dV/dS),
> ``gamma`` (d2V/dS2), ``vega`` (dV/dsigma), ``theta`` (calendar decay).
> Returns a dict with ``price`` and those fields.

### `average_strike_geometric_asian(S, t, r, sigma, n_fixings=None, fixing_times=None, option_type=<OptionType.CALL: 'call'>, b=None)`  _function_

> Average-strike (floating-strike) discrete geometric Asian option (exact).
>
> The strike is the realized geometric average rather than a fixed level:
> a call pays ``max(S_T - G, 0)`` and a put ``max(G - S_T, 0)``, where
> ``G = (prod_i S_{t_i})^{1/n}`` is the geometric average over the fixing
> dates. The terminal price ``S_T`` and the average ``G`` are jointly
> lognormal, so this is an exchange option between two lognormal assets and
> has an exact Margrabe-style closed form.
>
> With ``E[S_T] = S e^{bt}``, ``E[G] = exp(m_G + v_G/2)`` (the lognormal
> average forward), and the variance of ``log S_T - log G``
>
>     Var = sigma^2 t + v_G - 2 sigma^2 mean(t_i),
>     v_G = (sigma^2/n^2) sum_ij min(t_i, t_j),
>
> a call is ``e^{-rt} (E[S_T] N(d1) - E[G] N(d2))`` with
> ``d1 = (log(E[S_T]/E[G]) + Var/2)/sqrt(Var)`` and ``d2 = d1 - sqrt(Var)``.
>
> Provide either ``n_fixings`` (equally-spaced dates ``t*i/n``, last at
> expiry) or an explicit ``fixing_times`` sequence in ``(0, t]``.

### `average_strike_geometric_asian_greeks(S, t, r, sigma, n_fixings=None, fixing_times=None, option_type=<OptionType.CALL: 'call'>, b=None)`  _function_

> Greeks of an average-strike geometric Asian by central finite
> differences of :func:`average_strike_geometric_asian`: ``delta`` (dV/dS),
> ``gamma`` (d2V/dS2), ``vega`` (dV/dsigma), ``theta`` (calendar decay).
> Returns a dict with ``price`` and those fields.

### `barrier_greeks(S, K, H, t, r, sigma, option_type=<OptionType.CALL: 'call'>, barrier=<Barrier.DOWN_OUT: 'down-out'>, b=None, rebate=0.0)`  _function_

> Greeks of a single-barrier option by central finite differences.
>
> The Reiner-Rubinstein price is a closed form, but its Greeks are messy and
> change character across the barrier, so we central-difference the price.
> Returns a dict with delta, gamma, vega, and theta (calendar, per year).
>
> Near the barrier the true delta/gamma are large and discontinuous; the
> finite differences there are indicative rather than exact - use a spot bump
> well away from ``H`` when a smooth number is needed.

### `barrier_option(S, K, H, t, r, sigma, option_type=<OptionType.CALL: 'call'>, barrier=<Barrier.DOWN_OUT: 'down-out'>, b=None, rebate=0.0)`  _function_

> Price a single-barrier option with an optional cash rebate.
>
> Args:
>     H: barrier level.
>     barrier: one of the four Barrier kinds.
>     rebate: cash paid if the option is knocked out (out types) or never
>         knocked in (in types), paid at expiry.
>
> Implements the standard Reiner-Rubinstein decomposition. Validated in the
> suite against in-out parity (knock-in + knock-out = vanilla + rebate term).

### `barrier_rebate(S, H, t, r, sigma, knock='out', b=None, cash=1.0, payoff_at_hit=True)`  _function_

> Standalone rebate cashflow attached to a barrier.
>
> A **knock-out rebate** pays ``cash`` if the barrier ``H`` is breached (the
> consolation for the option knocking out); a **knock-in rebate** pays ``cash``
> at expiry if the barrier is *never* breached (the option failed to knock in).
>
> ``payoff_at_hit`` (knock-out only) pays on touch vs at expiry. This reuses
> the touch-option machinery: a knock-out rebate is a one-touch, a knock-in
> rebate is a no-touch.

### `cash_or_nothing(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None, cash=1.0)`  _function_

> Pays ``cash`` if the option finishes in the money, else 0.
>
> Call pays when S_T > K; put pays when S_T < K.

### `contingent_premium_option(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None)`  _function_

> Fair premium of a pay-later (contingent-premium) option.
>
> The holder pays no premium up front; instead a fixed premium is paid at expiry
> *only if* the option finishes in the money. For the deal to be fair at
> inception the premium's expected discounted value must equal the vanilla price:
>
>     vanilla = premium * cash_or_nothing(cash=1),
>
> so ``premium = vanilla / cash_or_nothing_unit``. The premium exceeds the
> vanilla price (it is only collected in the ITM states). Reduces toward the
> vanilla as the option goes deep in the money (ITM probability -> 1).

### `digital_greeks(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None, cash=1.0)`  _function_

> Delta and gamma of a cash-or-nothing digital by finite differences.
>
> Returns a dict with price, delta, and gamma. Near the strike as expiry
> approaches, the digital's delta spikes (and gamma flips sign across the
> strike) -- the "pin risk" that makes digitals hard to hedge and motivates
> the call-spread over-hedge in :mod:`quantforge.overhedge`.

### `discrete_arithmetic_asian(S, K, t, r, sigma, n_fixings=None, fixing_times=None, option_type=<OptionType.CALL: 'call'>, b=None)`  _function_

> Discretely-monitored arithmetic-average-price Asian option (Levy
> moment-matching approximation).
>
> The arithmetic average ``A = (1/n) sum_i S_{t_i}`` is not lognormal, but its
> first two moments over the fixing dates have exact closed forms:
>
>     M1 = (S/n) sum_i exp(b t_i)
>     M2 = (S^2/n^2) sum_i sum_j exp(b (t_i + t_j) + sigma^2 min(t_i, t_j)).
>
> Levy (1992) matches these to a lognormal and prices with a Black-Scholes
> formula on the average's forward ``M1`` and effective variance
> ``V = log(M2/M1^2)``:
>
>     d1 = (log(M1/K) + V/2) / sqrt(V),  d2 = d1 - sqrt(V)
>     call = e^{-rt} (M1 N(d1) - K N(d2)).
>
> Provide either ``n_fixings`` (equally-spaced dates ``t*i/n``, last at expiry)
> or an explicit ``fixing_times`` sequence in ``(0, t]``. A single fixing at
> ``t`` recovers the vanilla Black-Scholes price. The geometric-average
> Asian is an exact lower bound; this arithmetic price sits above it.

### `discrete_arithmetic_asian_greeks(S, K, t, r, sigma, n_fixings=None, fixing_times=None, option_type=<OptionType.CALL: 'call'>, b=None)`  _function_

> Greeks of a discrete arithmetic-average Asian option by central finite
> differences of :func:`discrete_arithmetic_asian`: ``delta`` (dV/dS),
> ``gamma`` (d2V/dS2), ``vega`` (dV/dsigma), ``theta`` (calendar decay).
> Returns a dict with ``price`` and those fields.

### `discrete_barrier_option(S, K, H, t, r, sigma, n_fixings, option_type=<OptionType.CALL: 'call'>, barrier=<Barrier.DOWN_OUT: 'down-out'>, b=None, rebate=0.0)`  _function_

> Discretely-monitored single-barrier option (Broadie-Glasserman-Kou 1999).
>
> A barrier checked at ``n_fixings`` equally-spaced dates is breached less
> often than one monitored continuously, so a knock-out is worth more and a
> knock-in worth less. Broadie-Glasserman-Kou give an asymptotic continuity
> correction: price with the continuous :func:`barrier_option` but shift the
> barrier away from the spot by ``exp(+/- beta sigma sqrt(dt))`` -- up for an
> up-barrier, down for a down-barrier -- with ``beta ~ 0.5826`` and
> ``dt = t / n_fixings``. As ``n_fixings -> infinity`` the shift vanishes and
> the price converges to the continuous barrier.
>
> Accurate to a fraction of a percent for ``n_fixings`` of ~50 or more; the
> correction is asymptotic, so a handful of monitoring dates carries a larger
> error. The shift is applied identically to knock-in and knock-out (they sum
> to the vanilla with the *same* shifted barrier via in-out parity).

### `discrete_geometric_asian(S, K, t, r, sigma, n_fixings=None, fixing_times=None, option_type=<OptionType.CALL: 'call'>, b=None)`  _function_

> Discretely-monitored geometric-average-price Asian option (exact).
>
> The geometric average ``G = (prod_i S_{t_i})^{1/n}`` over the monitoring
> dates ``t_i`` is lognormal, because ``log G`` is a linear combination of the
> jointly-Gaussian log-prices. With ``m = log S + (b - sigma^2/2) * mean(t_i)``
> and ``v = (sigma^2 / n^2) * sum_i sum_j min(t_i, t_j)``, ``log G`` is
> ``Normal(m, v)`` and the price is a Black-Scholes-style closed form on the
> forward ``F = exp(m + v/2)`` discounted at ``r``:
>
>     d1 = (m + v - log K) / sqrt(v),  d2 = d1 - sqrt(v)
>     call = e^{-rt} (F N(d1) - K N(d2)).
>
> Provide either ``n_fixings`` (equally-spaced dates ``t*i/n``, last at expiry)
> or an explicit ``fixing_times`` sequence in ``(0, t]``. A single fixing at
> ``t`` recovers the vanilla Black-Scholes price; as ``n_fixings -> infinity``
> the price converges to the continuous Kemna-Vorst :func:`geometric_asian`.

### `discrete_geometric_asian_greeks(S, K, t, r, sigma, n_fixings=None, fixing_times=None, option_type=<OptionType.CALL: 'call'>, b=None)`  _function_

> Greeks of a discrete geometric-average Asian option by central finite
> differences of :func:`discrete_geometric_asian`: ``delta`` (dV/dS),
> ``gamma`` (d2V/dS2), ``vega`` (dV/dsigma), ``theta`` (calendar decay). When
> ``fixing_times`` is given it is held fixed; with ``n_fixings`` the equally-
> spaced grid rescales with ``t`` (matching the continuous convention).
> Returns a dict with ``price`` and those fields.

### `double_knock_in_call(S, K, L, U, t, r, sigma, b=None, delta1=0.0, delta2=0.0, n_terms=10)`  _function_

> Double-barrier knock-in call: pays the call only if a barrier is touched.
>
> The in-out complement of :func:`double_knock_out_call`: a knock-in and a
> knock-out with the same strike and corridor partition every path, so at expiry
>
>     double_knock_in_call + double_knock_out_call = vanilla call.
>
> Priced as ``vanilla - double_knock_out_call`` with the vanilla evaluated on the
> same carry ``b``. Requires ``L < S < U``.

### `double_knock_in_call_greeks(S, K, L, U, t, r, sigma, b=None, delta1=0.0, delta2=0.0, n_terms=10)`  _function_

> Greeks of a double-barrier knock-in call by in-out parity.
>
> Differentiating ``knock_in = vanilla - knock_out`` term by term, the spot,
> vol and time Greeks are the vanilla Black-Scholes Greek minus the double
> knock-out Greek (:func:`double_knock_out_call_greeks`). The vanilla has no
> barrier dependence, so the knock-in's barrier sensitivities are the negatives
> of the knock-out's (widening the corridor lowers the knock-in). Returns the
> same dict layout: ``price``, ``delta``, ``gamma``, ``vega``, ``theta``,
> ``dV_dL``, ``dV_dU``.

### `double_knock_out_call(S, K, L, U, t, r, sigma, b=None, delta1=0.0, delta2=0.0, n_terms=10)`  _function_

> Ikeda-Kunitomo (1992) double-barrier knock-out call: payoff max(S_T - K, 0).
>
> Pays the vanilla call payoff only if the continuously-monitored spot stays
> strictly inside the (possibly exponentially curved) corridor bounded below by
> ``L e^{delta1 s}`` and above by ``U e^{delta2 s}`` over ``[0, t]``. With
> ``delta1 = delta2 = 0`` the barriers are flat at ``L`` and ``U``. The price is
> the Ikeda-Kunitomo image series
>
>     C = S e^{(b-r)t} sum_n [ (U^n/L^n)^{mu1} (L^n/S)^{mu2} (N(d1)-N(d2))
>                              - (L^{n+1}/(U^n S))^{mu3} (N(d3)-N(d4)) ]
>         - K e^{-rt} sum_n [ ... same with mu-2 and d-sigma sqrt(t) ... ],
>
> truncated at ``|n| <= n_terms`` (the series converges geometrically). Requires
> ``L < S < U`` and ``K < U`` for a non-trivial payoff.

### `double_knock_out_call_greeks(S, K, L, U, t, r, sigma, b=None, delta1=0.0, delta2=0.0, n_terms=10)`  _function_

> Greeks of an Ikeda-Kunitomo double knock-out call by finite differences.
>
> Central differences of :func:`double_knock_out_call` for ``delta`` (dV/dS),
> ``gamma`` (d2V/dS2), ``vega`` (dV/dsigma), ``theta`` (calendar decay
> ``-dV/dt``), and the two barrier sensitivities ``dV/dL`` and ``dV/dU``. The
> value is a knock-out, so more volatility raises the knock probability and
> ``vega < 0`` near the middle of the corridor, and widening either barrier
> raises the value (``dV/dL < 0``, ``dV/dU > 0``). Returns a dict with those
> fields.

### `double_no_touch(S, L, U, t, r, sigma, b=None, cash=1.0, n_terms=200)`  _function_

> Double-no-touch: pays ``cash`` at expiry if spot stays inside ``(L, U)``.
>
> Continuously monitored: the option survives only if the spot never touches
> either the lower barrier ``L`` or the upper barrier ``U`` before expiry. The
> survival probability of driftful Brownian motion in a strip has the classic
> Fourier (eigenfunction) expansion; with ``x = ln(S/L)``, ``Z = ln(U/L)``,
> ``m = b - sigma^2/2`` and ``beta = m/sigma^2``,
>
>     P(survive) = (2/Z) e^{-beta x - m^2 t / (2 sigma^2)}
>         * sum_{n>=1} sin(k_n x) e^{-k_n^2 sigma^2 t / 2}
>                      * k_n (1 - (-1)^n e^{beta Z}) / (beta^2 + k_n^2),
>
> with ``k_n = n pi / Z``. The value is ``cash e^{-rt} P(survive)``. As
> ``U -> infinity`` it collapses to the single lower :func:`no_touch`, and as
> ``L -> 0`` to the upper one. Requires ``L < S < U``.

### `double_no_touch_greeks(S, L, U, t, r, sigma, b=None, cash=1.0, n_terms=200)`  _function_

> Greeks of a double-no-touch by finite differences on the closed form.
>
> Central differences of :func:`double_no_touch` for ``delta`` (dV/dS),
> ``gamma`` (d2V/dS2), ``vega`` (dV/dsigma), ``theta`` (calendar decay
> ``-dV/dt``), and the two barrier sensitivities ``dV/dL`` and ``dV/dU``. A DNT
> is a bet on low realized range, so ``vega < 0`` (more vol -> more likely to
> knock) and widening either barrier raises the value (``dV/dL < 0`` since a
> lower ``L`` widens the band, ``dV/dU > 0``). Returns a dict with those fields.

### `double_one_touch(S, L, U, t, r, sigma, b=None, cash=1.0, n_terms=200)`  _function_

> Double-one-touch: pays ``cash`` at expiry if spot touches ``L`` or ``U``.
>
> The expiry-settled complement of :func:`double_no_touch`:
> ``double_one_touch = cash e^{-rt} - double_no_touch``. Requires ``L < S < U``.

### `gap_option(S, K_trigger, K_payoff, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None)`  _function_

> Gap option: pays off against ``K_payoff`` but is triggered by ``K_trigger``.
>
> A gap call pays ``S_T - K_payoff`` (which may be negative) whenever
> ``S_T > K_trigger``; a gap put pays ``K_payoff - S_T`` whenever
> ``S_T < K_trigger``. Setting the two strikes equal recovers the vanilla
> option. Closed form (Reiner-Rubinstein).

### `gap_option_greeks(S, K_trigger, K_payoff, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None)`  _function_

> Greeks of a gap option (Reiner-Rubinstein).
>
> ``delta`` and ``gamma`` are analytic; ``vega`` and ``theta`` (calendar
> decay) are central finite differences. A gap option is an asset-or-nothing
> minus ``K_payoff`` cash-or-nothings, both triggered at ``K_trigger``, so
> with ``d1, d2`` at the trigger and the identity
> ``S carry phi(d1) = K_trigger disc phi(d2)`` the spot sensitivities collapse
> to (call)
>
>     delta = carry N(d1) + disc (K_trigger - K_payoff) phi(d2)/(S sigma sqrt t)
>     gamma = carry phi(d1)/(S sigma sqrt t)
>             - disc (K_trigger - K_payoff) phi(d2) (d2/(sigma sqrt t) + 1)/(S^2 sigma sqrt t).
>
> Setting ``K_trigger = K_payoff`` recovers the vanilla Black-Scholes delta
> and gamma. Returns a dict with ``price`` and those fields.

### `geometric_asian(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None)`  _function_

> Continuously-monitored geometric-average-price Asian option.
>
> The geometric average of a lognormal is itself lognormal, so the price is a
> BSM price with adjusted volatility and carry:
>
>     sigma_A = sigma / sqrt(3)
>     b_A     = 0.5 * (b - sigma^2 / 6)
>
> (Kemna-Vorst). This gives an exact closed form and is a standard control
> variate for the arithmetic-average Asian priced by simulation.

### `geometric_asian_greeks(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None)`  _function_

> Greeks of a continuously-monitored geometric-average Asian option.
>
> The Kemna-Vorst price is exactly a Black-Scholes price with the adjusted
> volatility ``sigma_A = sigma / sqrt(3)`` and carry ``b_A = (b - sigma^2/6)/2``,
> so the spot ``S`` enters only through the BSM price at ``(sigma_A, b_A)``:
> ``delta`` and ``gamma`` are the exact BSM Greeks evaluated there (no finite
> difference). ``vega``, ``theta``, and ``rho`` do depend on ``sigma``/``t``/
> ``r`` through the adjusted parameters, so they are taken as central finite
> differences of the exact closed form. Returns a dict with ``price``,
> ``delta``, ``gamma``, ``vega``, ``theta``, ``rho``.

### `log_contract(S, t, r, sigma, b=None)`  _function_

> Log contract paying ``ln(S_T / F)`` at expiry, discounted to today.
>
> With ``F = S e^{b t}`` the forward, ``ln(S_T / F)`` is normal with mean
> ``-0.5 sigma^2 t`` under the risk-neutral measure, so the present value is
>
>     value = e^{-r t} * (-0.5 sigma^2 t).
>
> The log contract is the theoretical building block of the variance swap: a
> static log-contract position replicates the payoff of realized variance. The
> value is negative (the holder pays for the guaranteed negative drift of the
> log return).

### `log_contract_fair_variance(S, t, r, sigma, b=None)`  _function_

> Fair variance implied by the log contract: ``-2/t * e^{r t} * value``.
>
> Inverts :func:`log_contract` via the variance-swap replication identity
> ``sigma^2 = -2/t * E[ln(S_T / F)]``. Recovers the input ``sigma^2`` exactly in
> the Black-Scholes world -- the sanity check behind model-free variance-swap
> pricing.

### `no_touch(S, H, t, r, sigma, b=None, cash=1.0)`  _function_

> No-touch binary: pays ``cash`` at expiry if the barrier is never reached.
>
> Complementary to :func:`one_touch` with payment at expiry:
> ``no_touch = cash * e^{-rt} - one_touch(payoff_at_hit=False)``.

### `one_touch(S, H, t, r, sigma, b=None, cash=1.0, payoff_at_hit=True)`  _function_

> One-touch binary: pays ``cash`` if the spot ever reaches barrier ``H``.
>
> A continuously-monitored American digital. ``payoff_at_hit=True`` pays the
> cash immediately when the barrier is touched (the FX-market convention);
> ``False`` defers the payment to expiry. Works for an upper barrier
> (``H > S``) or a lower barrier (``H < S``); the direction is inferred.
>
> Uses the standard Rubinstein-Reiner touch formulas.

### `partial_time_end_barrier_call(S, K, H, t1, T2, r, sigma, barrier=<Barrier.DOWN_OUT: 'down-out'>, b=None) -> float`  _function_

> Partial-time (end) single-barrier call (Heynen-Kat 1994), closed form.
>
> The knock-out barrier is monitored only over ``[t1, T2]`` -- it is inactive
> before ``t1`` and live from ``t1`` to expiry ``T2``. Because the barrier
> watches a shorter window than a full-life barrier, a knock-out is worth more
> than the continuously-monitored one and less than the vanilla; as
> ``t1 -> 0`` it approaches the standard barrier and as ``t1 -> T2`` it
> approaches the vanilla call.
>
> Heynen-Kat's bivariate-normal formula couples the monitoring-start date
> ``t1`` (correlation ``rho = sqrt(t1/T2)``) to expiry. The down-out call is
> supported directly; the down-in value follows from in-out parity
> ``KI = vanilla - KO``. (Up-barrier partial-time calls have a distinct
> Heynen-Kat form and are not handled here.)

### `partial_time_start_barrier_call(S, K, H, t1, T2, r, sigma, barrier=<Barrier.DOWN_OUT: 'down-out'>, b=None) -> float`  _function_

> Partial-time (start) single-barrier call (Heynen-Kat 1994), closed form.
>
> The knock-out barrier is monitored only over ``[0, t1]`` -- it is live from
> inception to ``t1`` and inactive afterwards, with the option paying off at
> ``T2 > t1``. Because the barrier watches a shorter window than a full-life
> barrier, a knock-out is worth more than the continuously-monitored one and
> less than the vanilla; as ``t1 -> 0`` it approaches the vanilla call and as
> ``t1 -> T2`` it approaches the standard barrier.
>
> Heynen-Kat's bivariate-normal formula couples the monitoring-end date ``t1``
> (correlation ``rho = sqrt(t1/T2)``) to expiry. The down-out call is priced
> directly; the down-in value follows from in-out parity ``KI = vanilla - KO``.
> (Up-barrier partial-time calls have a distinct form and are not handled.)

### `pay_later_option_value(S, K, t, r, sigma, premium, option_type=<OptionType.CALL: 'call'>, b=None)`  _function_

> Value to the holder of a pay-later option with a contracted ``premium``.
>
> ``vanilla - premium * cash_or_nothing(cash=1)`` -- the option payoff net of the
> contingent premium collected only in the in-the-money states. Zero at the
> :func:`contingent_premium_option` fair premium, positive below it, negative
> above.

### `power_option(S, K, t, r, sigma, power, option_type=<OptionType.CALL: 'call'>, b=None)`  _function_

> Power option with payoff ``max(S_T^power - K, 0)`` (call) / ``max(K - S_T^power, 0)``.
>
> S_T^power is lognormal, so this has a closed form: an adjusted-drift,
> adjusted-vol Black-Scholes on the transformed underlying. ``power = 1``
> recovers the vanilla option.

### `power_option_greeks(S, K, t, r, sigma, power, option_type=<OptionType.CALL: 'call'>, b=None)`  _function_

> Greeks of a power option (payoff ``max(S_T^power - K, 0)``) by FD.
>
> Central finite differences of the closed-form :func:`power_option` for
> ``delta`` (dV/dS), ``gamma`` (d2V/dS2), ``vega`` (dV/dsigma), and ``theta``
> (calendar decay, ``-dV/dt``). At ``power = 1`` these reduce to the vanilla
> Black-Scholes Greeks. Returns a dict with ``price`` and those fields.

### `powered_option(S, K, t, r, sigma, power, option_type=<OptionType.CALL: 'call'>, b=None)`  _function_

> Powered option: payoff ``max(S_T - K, 0)**power`` (call) or
> ``max(K - S_T, 0)**power`` (put), for a positive **integer** ``power``.
>
> Distinct from :func:`power_option` (whose payoff is ``max(S_T**power - K, 0)``):
> here the *option payoff itself* is raised to a power, so the payoff has a
> higher-order convexity in the terminal spot. Because the payoff is a
> polynomial in ``S_T`` on the exercise region, it decomposes by the binomial
> theorem into a sum of ``S_T**j`` truncated moments, each of which has a
> closed form (Esser 2003; Heynen-Kat 1996). ``power = 1`` recovers the
> vanilla Black-Scholes option.
>
> With ``F_j = E[S_T**j] = S**j exp(j b t + 0.5 j (j-1) sigma^2 t)`` and
> ``d_j = (ln(S/K) + (b + (j - 0.5) sigma^2) t) / (sigma sqrt(t))``, a call is
> ``disc * sum_j C(p,j) (-K)^{p-j} F_j N(d_j)`` and a put is
> ``disc * sum_j C(p,j) K^{p-j} (-1)^j F_j N(-d_j)``.

### `powered_option_greeks(S, K, t, r, sigma, power, option_type=<OptionType.CALL: 'call'>, b=None)`  _function_

> Greeks of a powered option by central finite differences of
> :func:`powered_option`: ``delta`` (dV/dS), ``gamma`` (d2V/dS2), ``vega``
> (dV/dsigma), ``theta`` (calendar decay, ``-dV/dt``). Returns a dict with
> ``price`` and those fields.

### `range_binary(S, K_low, K_high, t, r, sigma, b=None, cash=1.0)`  _function_

> Range binary (double digital): pays ``cash`` iff ``K_low <= S_T <= K_high``.
>
> A bet that the terminal price lands inside a corridor at expiry (no path
> monitoring). It is exactly the difference of two cash-or-nothing calls
> struck at the two levels: paying when ``S_T > K_low`` but not when
> ``S_T > K_high``. With ``d2(K) = (log(S/K) + (b - sigma^2/2) t)/(sigma sqrt t)``,
>
>     price = cash e^{-rt} (N(d2(K_low)) - N(d2(K_high))).

### `range_binary_greeks(S, K_low, K_high, t, r, sigma, b=None, cash=1.0)`  _function_

> Greeks of a range binary (double digital).
>
> ``delta`` and ``gamma`` are analytic. With ``d2(K)`` the digital exponent
> and ``e^{-rt}`` the discount, the corridor value is
> ``cash e^{-rt} (N(d2_lo) - N(d2_hi))``, so differentiating in spot,
>
>     delta = cash e^{-rt} (phi(d2_lo) - phi(d2_hi)) / (S sigma sqrt(t))
>     gamma = -cash e^{-rt} / (S^2 sigma sqrt(t))
>             * ((phi(d2_lo) d2_lo - phi(d2_hi) d2_hi) / (sigma sqrt(t))
>                + phi(d2_lo) - phi(d2_hi)).
>
> ``vega`` and ``theta`` (calendar decay) are central finite differences of
> :func:`range_binary`. Returns a dict with ``price`` and those fields. The
> delta changes sign across the middle of the corridor and gamma is large near
> either edge as expiry approaches (double-sided pin risk).

### `seasoned_arithmetic_asian(S, K, t, r, sigma, observed_prices, n_total, remaining_fixing_times=None, option_type=<OptionType.CALL: 'call'>, b=None)`  _function_

> Seasoned (in-progress) discrete arithmetic-average Asian (Levy match).
>
> Prices an arithmetic Asian partway through its averaging window, when some
> fixings are already observed. The final average is
> ``A = (Q + sum_j S_{tau_j}) / n`` where ``Q = sum(observed_prices)`` is a
> known constant and the ``k`` remaining prices are lognormal. A call payoff
> ``max(A - K, 0) = (1/n) max(sum_j S_{tau_j} - (n K - Q), 0)`` is therefore an
> arithmetic-average option on the *remaining* fixings with the shifted strike
> ``K' = n K - Q``, scaled by ``1/n``. The remaining sum's first two moments
> are exact,
>
>     m1 = S sum_j exp(b tau_j)
>     m2 = S^2 sum_ij exp(b (tau_i + tau_j) + sigma^2 min(tau_i, tau_j)),
>
> and Levy (1992) matches them to a lognormal priced by Black-Scholes.
>
> Special cases handled exactly: if ``K' <= 0`` the call is always in the
> money and worth ``e^{-rt} (E[A] - K)`` (the put is worthless), and vice
> versa. With no observations this reduces to :func:`discrete_arithmetic_asian`;
> with all fixings observed the payoff is the deterministic arithmetic
> intrinsic.

### `seasoned_arithmetic_asian_greeks(S, K, t, r, sigma, observed_prices, n_total, remaining_fixing_times=None, option_type=<OptionType.CALL: 'call'>, b=None)`  _function_

> Greeks of a seasoned arithmetic Asian by central finite differences of
> :func:`seasoned_arithmetic_asian`: ``delta`` (dV/dS), ``gamma`` (d2V/dS2),
> ``vega`` (dV/dsigma), ``theta`` (calendar decay). The observed fixings and
> the remaining schedule are held fixed. Returns a dict with ``price`` and
> those fields.

### `seasoned_geometric_asian(S, K, t, r, sigma, observed_prices, n_total, remaining_fixing_times=None, option_type=<OptionType.CALL: 'call'>, b=None)`  _function_

> Seasoned (in-progress) discrete geometric-average Asian option (exact).
>
> Prices a geometric Asian partway through its averaging window, when some
> fixings have already been observed. With ``n_total`` fixings in all, of
> which ``observed_prices`` are already fixed and ``k`` remain at future
> times ``remaining_fixing_times`` (in ``(0, t]``, measured from now), the
> final average ``G = (prod_{i=1}^{n} S_{t_i})^{1/n}`` is still lognormal: the
> observed factors contribute a known constant ``A = sum log(S_obs)`` and the
> ``k`` future log-prices are jointly Gaussian. Hence ``log G`` is
> ``Normal(m, v)`` with
>
>     m = (A + k log S + (b - sigma^2/2) sum_j tau_j) / n
>     v = (sigma^2 / n^2) sum_i sum_j min(tau_i, tau_j),
>
> and the price is a Black-Scholes-style closed form on ``F = exp(m + v/2)``
> discounted at ``r``. When no fixings are observed this reduces exactly to
> :func:`discrete_geometric_asian`. When all ``n_total`` fixings are observed
> the average is known and the payoff is deterministic.

### `seasoned_geometric_asian_greeks(S, K, t, r, sigma, observed_prices, n_total, remaining_fixing_times=None, option_type=<OptionType.CALL: 'call'>, b=None)`  _function_

> Greeks of a seasoned geometric Asian by central finite differences of
> :func:`seasoned_geometric_asian`: ``delta`` (dV/dS), ``gamma`` (d2V/dS2),
> ``vega`` (dV/dsigma), ``theta`` (calendar decay). The observed fixings and
> the remaining schedule are held fixed. Returns a dict with ``price`` and
> those fields.

### `supershare(S, K_low, K_high, t, r, sigma, b=None)`  _function_

> Supershare option: pays ``S_T / K_low`` iff ``K_low <= S_T <= K_high``.
>
> Introduced by Hakansson (1976) as a building block for mutual-fund payoffs.
> It is a scaled difference of two asset-or-nothing calls struck at the two
> levels, so with ``d1(K) = (log(S/K) + (b + sigma^2/2) t)/(sigma sqrt t)`` and
> carry ``e^{(b-r)t}``,
>
>     price = (S/K_low) e^{(b-r)t} (N(d1(K_low)) - N(d1(K_high))).

### `supershare_greeks(S, K_low, K_high, t, r, sigma, b=None)`  _function_

> Greeks of a supershare option by central finite differences of
> :func:`supershare`: ``delta`` (dV/dS), ``gamma`` (d2V/dS2), ``vega``
> (dV/dsigma), ``theta`` (calendar decay). Returns a dict with ``price`` and
> those fields.

## extendible

### `holder_extendible_call(S, K1, K2, t1, T2, r, sigma, A, b=None) -> float`  _function_

> Holder-extendible call (Longstaff 1990), closed form.
>
> Args:
>     K1: strike applying at the first expiry ``t1``.
>     K2: strike of the extended option (life ``T2``).
>     t1: first expiry (years); T2: extended expiry (years, ``> t1``).
>     A: fee paid at ``t1`` to extend.
>
> Reduces to a vanilla call struck at ``K1`` expiring at ``t1`` when extending
> is never worthwhile (very large ``A``).

### `holder_extendible_call_greeks(S, K1, K2, t1, T2, r, sigma, A, b=None)`  _function_

> Greeks of a holder-extendible call by central finite differences of
> :func:`holder_extendible_call`: ``delta``, ``gamma``, ``vega``, ``theta``
> (calendar decay, both expiries shrinking together). Returns a dict with
> ``price`` and those fields.

### `holder_extendible_put(S, K1, K2, t1, T2, r, sigma, A, b=None) -> float`  _function_

> Holder-extendible put (Longstaff 1990), closed form.
>
> At the first expiry ``t1`` the holder takes the best of exercising against
> ``K1`` for ``K1 - S_{t1}``, lapsing, or paying a fee ``A`` to extend to
> ``T2`` as a put struck at ``K2``:
>
>     payoff(t1) = max(K1 - S_{t1}, P(S_{t1}, K2, T2 - t1) - A, 0).
>
> The terminal spot splits into exercise (``S < I_low``), extend
> (``I_low <= S <= I_high``), and lapse (``S > I_high``). The strip is valued
> with bivariate normals coupling ``t1`` and ``T2``. A very large fee collapses
> the strip and recovers the vanilla put struck at ``K1`` expiring at ``t1``.

### `holder_extendible_put_greeks(S, K1, K2, t1, T2, r, sigma, A, b=None)`  _function_

> Greeks of a holder-extendible put by central finite differences of
> :func:`holder_extendible_put`: ``delta``, ``gamma``, ``vega``, ``theta``
> (calendar decay, both expiries shrinking together). Returns a dict with
> ``price`` and those fields.

### `writer_extendible_call(S, K1, K2, t1, T2, r, sigma, b=None) -> float`  _function_

> Writer-extendible call (Longstaff 1990), closed form.
>
> At the first expiry ``t1`` the call is exercised if it finishes in the money
> (``S_{t1} > K1``, paying ``S_{t1} - K1``); otherwise the writer's obligation
> is automatically extended to ``T2`` as a call struck at ``K2`` (no fee):
>
>     payoff(t1) = (S_{t1} - K1)            if S_{t1} > K1
>                = C(S_{t1}, K2, T2 - t1)   if S_{t1} <= K1.
>
> The value is a vanilla call to ``t1`` plus the extended-call value collected
> on ``S_{t1} <= K1``, via bivariate normals coupling ``t1`` and ``T2``:
>
>     W = c(S, K1, t1)
>         + S e^{(b-r)T2} M(-z1, y1; -rho) - K2 e^{-r T2} M(-z2, y2; -rho),
>
> with ``z1, z2`` the ``d1/d2`` arguments at ``K1`` over ``t1``, ``y1, y2`` the
> same at ``K2`` over ``T2``, and ``rho = sqrt(t1/T2)``.

### `writer_extendible_call_greeks(S, K1, K2, t1, T2, r, sigma, b=None)`  _function_

> Greeks of a writer-extendible call by central finite differences of
> :func:`writer_extendible_call`: ``delta``, ``gamma``, ``vega``, ``theta``
> (calendar decay, both expiries shrinking together). Returns a dict with
> ``price`` and those fields.

### `writer_extendible_put(S, K1, K2, t1, T2, r, sigma, b=None) -> float`  _function_

> Writer-extendible put (Longstaff 1990), closed form.
>
> At the first expiry ``t1`` the put is exercised if it finishes in the money
> (``S_{t1} < K1``, paying ``K1 - S_{t1}``); otherwise the writer's obligation
> is automatically extended to ``T2`` as a put struck at ``K2`` (no fee). The
> terminal-``t1`` payoff is therefore
>
>     payoff(t1) = (K1 - S_{t1})           if S_{t1} < K1
>                = P(S_{t1}, K2, T2 - t1)  if S_{t1} >= K1.
>
> The value is a vanilla put to ``t1`` plus the extended-put value collected on
> ``S_{t1} >= K1``, expressed with bivariate normals coupling ``t1`` and
> ``T2`` (correlation ``rho = sqrt(t1/T2)``):
>
>     W = p(S, K1, t1)
>         + K2 e^{-r T2} M(z2, -y2; -rho) - S e^{(b-r)T2} M(z1, -y1; -rho),
>
> with ``z2, y2`` the ``d2``-type arguments at ``K1`` (over ``t1``) and ``K2``
> (over ``T2``) and ``z1 = z2 + sigma sqrt(t1)``, ``y1 = y2 + sigma sqrt(T2)``.

### `writer_extendible_put_greeks(S, K1, K2, t1, T2, r, sigma, b=None)`  _function_

> Greeks of a writer-extendible put by central finite differences of
> :func:`writer_extendible_put`: ``delta``, ``gamma``, ``vega``, ``theta``
> (calendar decay, both expiries shrinking together). Returns a dict with
> ``price`` and those fields.

## forward

### `ForwardResult(forward: float, discount_factor: float, implied_rate: float, implied_div_yield: float, n_strikes: int, rmse: float) -> None`  _class_

> ForwardResult(forward: float, discount_factor: float, implied_rate: float, implied_div_yield: float, n_strikes: int, rmse: float)

### `dividend_curve(chain_by_expiry, spot)`  _function_

> Bootstrap an implied dividend-yield term structure from a multi-expiry chain.
>
> Args:
>     chain_by_expiry: iterable of ``(t, strikes, calls, puts)`` tuples, one
>         per expiry.
>     spot: current underlying spot.
>
> Returns a list of ``(t, ForwardResult)`` pairs sorted by expiry, each from
> :func:`implied_forward`. The ``implied_div_yield`` field of each result is
> the continuous dividend yield to that expiry (a point on the dividend curve).

### `implied_forward(strikes: Sequence[float], calls: Sequence[float], puts: Sequence[float], t: float, spot: float = None)`  _function_

> Extract the implied forward and discount factor from a parity fit.
>
> Solves ``C - P = D*F - D*K`` as a straight line in ``K`` by ordinary least
> squares: the slope is ``-D`` and the intercept is ``D*F``.
>
> Args:
>     strikes, calls, puts: equal-length chains at a single expiry.
>     t: time to expiry in years (used to annualize the implied rate).
>     spot: if given, also returns the implied continuous dividend yield.
>
> Returns a :class:`ForwardResult`.

## forwardstart

### `cliquet_greeks(S, reset_times: Sequence[float], r, sigma, alpha=1.0, option_type=<OptionType.CALL: 'call'>, b=None)`  _function_

> Greeks of a cliquet (ratchet) by central finite differences.
>
> A cliquet is a strip of consecutive forward-start options. Only the first
> (spot-strike) period carries spot gamma; every later forward-start period is
> linear in the current spot, so the cliquet's ``gamma`` comes entirely from
> the first period and is small relative to a single vanilla. ``delta``,
> ``gamma``, ``vega``, and ``theta`` are central finite differences of
> :func:`cliquet_price`; ``theta`` shifts every reset date together. Returns a
> dict with ``price``, ``delta``, ``gamma``, ``vega``, ``theta``.

### `cliquet_price(S, reset_times: Sequence[float], r, sigma, alpha=1.0, option_type=<OptionType.CALL: 'call'>, b=None) -> float`  _function_

> Price a cliquet (ratchet) as a strip of forward-start options.
>
> Args:
>     reset_times: increasing schedule of reset/expiry dates in years, e.g.
>         [0.25, 0.5, 0.75, 1.0]. Each consecutive pair (t_i, t_{i+1}) is one
>         forward-start period that starts at t_i and expires at t_{i+1}. The
>         first period starts now (t=0) and ends at reset_times[0].
>
> Returns the total present value of the strip.

### `forward_start_greeks(S, t_start, t_expiry, r, sigma, alpha=1.0, option_type=<OptionType.CALL: 'call'>, b=None)`  _function_

> Greeks of a forward-start option (Rubinstein), exact where possible.
>
> The price is ``FS = S e^{(b-r) t_start} * u`` where ``u`` is a unit
> Black-Scholes price on a unit underlying and does **not** depend on ``S``.
> So the value is exactly linear in the spot: ``delta = e^{(b-r) t_start} u``
> (constant in ``S``) and ``gamma = 0`` -- a forward-start has no spot gamma
> until its strike is fixed. ``vega`` and ``theta`` (calendar decay, both
> ``t_start`` and ``t_expiry`` shifting together) are central finite
> differences of the closed form. Returns a dict with ``price``, ``delta``,
> ``gamma``, ``vega``, ``theta``.

### `forward_start_price(S, t_start, t_expiry, r, sigma, alpha=1.0, option_type=<OptionType.CALL: 'call'>, b=None) -> float`  _function_

> Price a forward-start option whose strike is set at ``t_start``.
>
> Args:
>     S: current spot.
>     t_start: time (years) until the strike is fixed. 0 reduces to a vanilla.
>     t_expiry: total time (years) to the option's expiry (> t_start).
>     alpha: moneyness multiple; strike = alpha * S_{t_start}. alpha=1 is ATM.
>     b: cost of carry (defaults to r).
>
> Returns the present value.

## futures_convexity

### `forward_curve_from_futures_strip(futures_quotes, sigma, a=0.0)`  _function_

> Convert a strip of futures quotes to forward rates via convexity.
>
> ``futures_quotes`` is ``[(t1, t2, futures_rate), ...]`` for consecutive
> contracts. Applies :func:`forward_from_futures` to each, returning
> ``[(t1, t2, forward_rate), ...]``. Every forward sits below its futures rate,
> and the adjustment grows down the curve.

### `forward_from_futures(futures_rate, sigma, t1, t2, a=0.0)`  _function_

> Forward rate from a futures rate, subtracting the convexity adjustment.
>
> Uses Hull-White (or Ho-Lee when ``a = 0``). The forward is below the futures
> rate by the (non-negative) adjustment.

### `futures_from_forward(forward_rate, sigma, t1, t2, a=0.0)`  _function_

> Futures rate from a forward rate, adding the convexity adjustment.
>
> Inverse of :func:`forward_from_futures`; the futures rate is above the forward.

### `ho_lee_convexity_adjustment(sigma, t1, t2)`  _function_

> Ho-Lee convexity adjustment ``0.5 sigma^2 t1 t2`` (continuous-comp rates).
>
> ``t1`` is the time to the rate fixing, ``t2`` the time to the end of the
> underlying accrual period (``t2 >= t1``). Non-negative, zero at zero vol, and
> growing with both horizons and the volatility.

### `hull_white_convexity_adjustment(sigma, a, t1, t2)`  _function_

> Hull-White (mean-reverting) convexity adjustment.
>
> With mean-reversion speed ``a`` the adjustment is
>
>     adj = (sigma^2 / (2 a^2)) * (1 - e^{-a(t2 - t1)}) *
>           [ (1 - e^{-a(t2 - t1)}) * (1 - e^{-2 a t1}) / (2 a)
>             + a * B(0, t1)^2 ... ]  (standard HW futures-forward formula)
>
> Implemented in the common compact form
>     adj = (sigma^2 / (2 a)) * B(t1, t2) * [ B(t1, t2) (1 - e^{-2 a t1})
>           + 2 a B(0, t1)^2 ] / 2,
> with ``B(u, v) = (1 - e^{-a (v - u)}) / a``. Reduces to
> :func:`ho_lee_convexity_adjustment` as ``a -> 0``.

### `stub_discount_factors_from_forwards(forward_quotes, df0=1.0)`  _function_

> Bootstrap discount factors from a strip of forward rates.
>
> ``forward_quotes`` is ``[(t1, t2, forward_rate), ...]`` of consecutive simple
> forward rates over ``[t1, t2]``. Chains ``DF(t2) = DF(t1) / (1 + f * (t2 -
> t1))`` starting from ``df0`` at the first ``t1``. Returns ``[(t2, DF(t2)),
> ...]``; discount factors are decreasing for positive rates.

## fxdelta

### `atm_dns_strike(F, t, sigma)`  _function_

> Delta-neutral-straddle ATM strike ``K = F exp(0.5 sigma^2 t)``.
>
> The strike at which a straddle has zero (forward) delta -- the market's
> standard ATM quote for most currency pairs.

### `delta_from_strike(F, t, sigma, K, is_call, premium_adjusted=False, spot_delta=False, r_for=0.0)`  _function_

> Forward (or spot) delta of an option struck at ``K``.
>
> The inverse of :func:`strike_from_delta`; returns the signed delta.

### `rr_bf_to_pillars(F, t, atm, rr, bf, call_delta=0.25)`  _function_

> Convert (ATM, risk-reversal, butterfly) quotes to smile pillars.
>
> Returns ``(K_put, sigma_put, K_atm, sigma_atm, K_call, sigma_call)`` -- the
> three market pillar strikes and vols. Uses the standard smile-implied-from-
> quotes relations ``sigma_25c = atm + bf + rr/2``, ``sigma_25p = atm + bf -
> rr/2`` and the delta-neutral ATM strike.

### `strike_from_delta(F, t, sigma, delta, is_call, premium_adjusted=False, spot_delta=False, r_for=0.0)`  _function_

> Strike with the given delta.
>
> Args:
>     delta: the target delta magnitude convention -- pass the signed delta
>         (call > 0, put < 0), e.g. ``0.25`` for a 25-delta call, ``-0.25``
>         for a 25-delta put.
>     is_call: whether the option is a call.
>     premium_adjusted: use the premium-adjusted delta convention (delta net
>         of the option premium, standard for premium-in-foreign pairs).
>     spot_delta: if True the delta is a spot delta (discounted by ``r_for``);
>         otherwise a forward delta.
>
> Returns the strike ``K``.

## fxforward

### `forward_points(spot, r_price, r_base, t) -> float`  _function_

> Forward points ``F - S`` (positive when the base is at a forward premium).

### `fx_forward(spot, r_price, r_base, t) -> float`  _function_

> Covered-interest-parity forward FX rate ``S exp((r_price - r_base) t)``.
>
> A base currency yielding more than the price currency (``r_base > r_price``)
> trades at a forward discount (``F < S``), and vice versa.

### `fx_forward_from_curves(spot, price_curve, base_curve, t) -> float`  _function_

> FX forward from two discount curves: ``S * DF_base(t) / DF_price(t)``.
>
> Curve-based covered interest parity: a base currency that discounts more
> steeply (higher rates, lower ``DF_base``) trades at a forward discount. Each
> curve is anything callable as ``curve.df(t)`` (e.g.
> :class:`quantforge.DiscountCurve`) or a plain ``curve(t)`` returning the
> discount factor. Consistent with :func:`fx_forward` when the curves are flat
> exponentials.

### `fx_swap_points(spot, r_price, r_base, t_near, t_far) -> float`  _function_

> FX-swap points between two tenors: ``F(t_far) - F(t_near)``.
>
> The pips exchanged in a forward-forward FX swap rolling from the near to the
> far date.

### `implied_base_rate(spot, forward, r_price, t) -> float`  _function_

> Base-currency rate implied by a quoted forward (invert CIP).
>
> ``r_base = r_price - ln(forward / spot) / t``.

### `implied_price_rate(spot, forward, r_base, t) -> float`  _function_

> Price-currency rate implied by a quoted forward (invert CIP).
>
> ``r_price = r_base + ln(forward / spot) / t``.

## g2pp

### `g2pp_V(a, b, sigma, eta, rho, t, T)`  _function_

> The G2++ variance term V(t,T) (Brigo-Mercurio eq. 4.10).

### `g2pp_bond_option(P0S, P0T, a, b, sigma, eta, rho, expiry, maturity, strike, is_call=True)`  _function_

> European option on a zero-coupon bond under G2++ (exact).
>
> Option expires at ``expiry`` on a bond maturing at ``maturity``, struck at
> ``strike``. ``P0S = P(0, expiry)``, ``P0T = P(0, maturity)``. Since
> ``ln P(expiry, maturity)`` is Gaussian the price is a Black-style formula on
> the forward bond ``P0T / P0S`` with the G2++ bond volatility.

### `g2pp_bond_option_greeks(P0S, P0T, a, b, sigma, eta, rho, expiry, maturity, strike, is_call=True)`  _function_

> Greeks of a G2++ zero-coupon-bond option.
>
> Sensitivities of :func:`bond_option` to the two discount factors and the
> model vols, by central finite differences except the two discount-factor
> deltas which are exact (the price is Black-style in ``P0T``/``P0S``):
>
>   * ``delta_T`` = dV/dP0T = ``N(d1)`` (call) -- the underlying-bond delta;
>   * ``delta_S`` = dV/dP0S (the discount-leg delta);
>   * ``vega_sigma`` = dV/dsigma, ``vega_eta`` = dV/deta -- exposure to the two
>     G2++ factor vols.
>
> Returns a dict with ``price``, ``delta_T``, ``delta_S``, ``vega_sigma``,
> ``vega_eta``.

### `g2pp_cap(discounts, a, b, sigma, eta, rho, strike, notional=1.0)`  _function_

> G2++ cap: strip of caplets over successive periods.
>
> ``discounts`` is an increasing list of ``(t_i, P(0, t_i))`` reset/pay dates
> (the first pair is the first reset, then each consecutive pair is a caplet
> ``[t_{i-1}, t_i]``). Returns the summed caplet value.

### `g2pp_caplet(P0_reset, P0_pay, a, b, sigma, eta, rho, reset, pay, strike, notional=1.0)`  _function_

> Caplet on ``[reset, pay]`` under G2++ via the bond-put identity.

### `g2pp_floor(discounts, a, b, sigma, eta, rho, strike, notional=1.0)`  _function_

> G2++ floor: strip of floorlets over successive periods (see :func:`cap`).

### `g2pp_floorlet(P0_reset, P0_pay, a, b, sigma, eta, rho, reset, pay, strike, notional=1.0)`  _function_

> Floorlet on ``[reset, pay]`` under G2++ via the bond-call identity.

### `g2pp_zero_bond(P0T, P0t, x, y, a, b, sigma, eta, rho, t, T)`  _function_

> G2++ zero-coupon bond ``P(t,T)`` given the factor state ``(x, y)``.

## gramcharlier

### `calibrate_corrado_su(S, t, r, strikes, call_prices, b=None, initial=None, max_iter=4000) -> Tuple[float, float, float, float]`  _function_

> Fit Corrado-Su ``(sigma, skew, excess_kurt)`` to market call prices.
>
> Minimizes the sum of squared price errors of :func:`corrado_su_call` over the
> given strikes with Nelder-Mead, using a smooth reparametrization that keeps
> ``sigma > 0`` (the skew and excess-kurtosis coefficients are unconstrained).
> A flat Black-Scholes surface calibrates to ``skew = kurt = 0`` and the input
> vol.
>
> Returns ``(sigma, skew, excess_kurt, price_rmse)``.

### `corrado_su_call(S, K, t, r, sigma, skew=0.0, excess_kurt=0.0, b=None) -> float`  _function_

> Corrado-Su (1996) skew/kurtosis-adjusted European call price.
>
> Args:
>     skew: skewness of the (log) return distribution.
>     excess_kurt: excess kurtosis (kurtosis - 3).
>     b: cost of carry (defaults to r). skew=kurt=0 => Black-Scholes.

### `corrado_su_implied_vol(S, K, t, r, sigma, skew=0.0, excess_kurt=0.0, b=None) -> float`  _function_

> Black-Scholes implied vol of a Corrado-Su price at a single strike.
>
> Prices the option with :func:`corrado_su_price` at the Gram-Charlier
> parameters, then inverts Black-Scholes for the vol that reproduces it. With
> ``skew = excess_kurt = 0`` this returns ``sigma`` at every strike (a flat
> smile); non-zero moments trace the characteristic Gram-Charlier skew/smile:
> negative skew lifts the low-strike (put) wing, positive excess kurtosis lifts
> both wings relative to the at-the-money level.

### `corrado_su_price(S, K, t, r, sigma, skew=0.0, excess_kurt=0.0, option_type=<OptionType.CALL: 'call'>, b=None) -> float`  _function_

> Corrado-Su price for a call or put (put via put-call parity).

### `corrado_su_smile(S, t, r, sigma, strikes, skew=0.0, excess_kurt=0.0, b=None)`  _function_

> Corrado-Su implied-vol smile: ``(strikes, implied_vols)`` over ``strikes``.
>
> Convenience wrapper mapping each strike through
> :func:`corrado_su_implied_vol`. Useful for plotting the skew/kurtosis smile
> or seeding an SVI/SABR fit from Gram-Charlier moments.

### `realized_excess_kurtosis(returns: Sequence[float]) -> float`  _function_

> Sample excess kurtosis (kurtosis - 3) of a return series.

### `realized_skewness(returns: Sequence[float]) -> float`  _function_

> Sample skewness of a return series (bias-corrected denominator n).

## greeks2

### `charm(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None) -> float`  _function_

> Calendar charm = -d(delta)/d(t_expiry): the drift of delta over time.

### `color(S, K, t, r, sigma, b=None) -> float`  _function_

> Calendar color = -d(gamma)/d(t_expiry): the decay of gamma over time.

### `dual_delta(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None) -> float`  _function_

> d(price)/d(strike). For a call ``-e^{-rt} N(d2)``, for a put
> ``e^{-rt} N(-d2)``.
>
> The strike sensitivity. Its negative (for a call) is the discounted
> risk-neutral probability of finishing in the money, so
> ``-dual_delta_call = e^{-rt} N(d2)`` is exactly the cash-or-nothing digital
> value -- the link to :func:`quantforge.risk_neutral_cdf`.

### `dual_gamma(S, K, t, r, sigma, b=None) -> float`  _function_

> d^2(price)/d(strike)^2 = e^{-rt} phi(d2) / (K sigma sqrt(t)).
>
> Same for calls and puts. By Breeden-Litzenberger this is exactly the
> discounted risk-neutral probability density of the terminal price at ``K``,
> so it is always non-negative in an arbitrage-free market.

### `speed(S, K, t, r, sigma, b=None) -> float`  _function_

> d(gamma)/d(spot). Third-order in spot; same for calls and puts.

### `ultima(S, K, t, r, sigma, b=None) -> float`  _function_

> d(vomma)/d(sigma) = d^3(price)/d(sigma)^3. Same for calls and puts.
>
> The third-order vega sensitivity, useful for the convexity of the volga
> hedge. With vega ``v = S e^{(b-r)t} phi(d1) sqrt(t)``,
>
>     ultima = -(v / sigma^2)
>              * [ d1 d2 (1 - d1 d2) + d1^2 + d2^2 ].

### `vanna(S, K, t, r, sigma, b=None) -> float`  _function_

> d(delta)/d(sigma) = d(vega)/d(spot). Same for calls and puts.

### `veta(S, K, t, r, sigma, b=None) -> float`  _function_

> Calendar veta = -d(vega)/d(t_expiry): decay of vega over time.

### `volga(S, K, t, r, sigma, b=None) -> float`  _function_

> d(vega)/d(sigma) (volga). Same for calls and puts.

### `vomma(S, K, t, r, sigma, b=None) -> float`  _function_

> d(vega)/d(sigma) (volga). Same for calls and puts.

### `zomma(S, K, t, r, sigma, b=None) -> float`  _function_

> d(gamma)/d(sigma). Same for calls and puts.

## hedgesim

### `HedgeResult(mean_pnl: float, std_pnl: float, min_pnl: float, max_pnl: float, n_paths: int, n_steps: int) -> None`  _class_

> HedgeResult(mean_pnl: float, std_pnl: float, min_pnl: float, max_pnl: float, n_paths: int, n_steps: int)

### `simulate_delta_hedge(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None, n_steps=50, n_paths=20000, hedge_vol=None, real_vol=None, seed=None, return_samples=False)`  _function_

> Monte Carlo a discretely delta-hedged short option position.
>
> Args:
>     n_steps: rebalancing dates over the option's life.
>     hedge_vol: volatility used to compute the hedge delta (defaults to
>         ``sigma``). Set different from ``real_vol`` to study hedging at the
>         wrong vol.
>     real_vol: volatility of the simulated path (defaults to ``sigma``).
>     return_samples: if True, also return the raw per-path P&L list.
>
> Returns a :class:`HedgeResult` (and the samples if requested).

## hedging

### `StickyRule(*values)`  _class_

> str(object='') -> str
> str(bytes_or_buffer[, encoding[, errors]]) -> str
>
> Create a new string object from the given object. If encoding or
> errors is specified, then the object must expose a data buffer
> that will be decoded using the given encoding and error handler.
> Otherwise, returns the result of object.__str__() (if defined)
> or repr(object).
> encoding defaults to 'utf-8'.
> errors defaults to 'strict'.

### `skew_slope(smile_fn, K, F, h=None)`  _function_

> Estimate d(sigma)/dk at strike ``K`` via central finite difference.
>
> Args:
>     smile_fn: callable ``sigma(K)`` returning implied vol for a strike.
>     K: strike at which to measure the slope.
>     F: forward (used to convert to log-moneyness k = ln(K/F)).
>     h: bump in ``k`` space; defaults to a small fraction.
>
> Returns d(sigma)/dk where k = ln(K/F).

### `smile_delta(S, K, t, r, sigma, dsigma_dk=0.0, option_type=<OptionType.CALL: 'call'>, b=None, sticky=<StickyRule.DELTA: 'delta'>) -> float`  _function_

> Effective (smile-adjusted) delta of an option.
>
> Args:
>     sigma: the option's current implied volatility.
>     dsigma_dk: local skew slope d(sigma)/dk at this strike, where
>         k = ln(K/F). Only used under the sticky-delta rule.
>     sticky: STRIKE (delta == BS delta) or DELTA (add the vega/skew term).
>
> Under sticky-delta the smile is a function of moneyness, so a 1-unit rise in
> spot lowers the log-moneyness of a fixed strike by 1/S, shifting its vol by
> ``-(dsigma/dk)/S``. The effective delta is therefore
>
>     delta_eff = delta_BS + vega * d(sigma)/d(spot)
>               = delta_BS - vega * (dsigma/dk) / S.

### `smile_delta_from_smile(S, K, t, r, smile_fn, F=None, option_type=<OptionType.CALL: 'call'>, b=None) -> float`  _function_

> Sticky-delta effective delta computed directly from a smile function.
>
> Reads the option's vol as ``smile_fn(K)``, estimates the local skew slope by
> finite difference, and returns the adjusted delta. ``F`` defaults to the
> carry-implied forward ``S * exp(b * t)``.

## heston

### `heston_price(S, K, t, r, v0, kappa, theta, xi, rho, option_type=<OptionType.CALL: 'call'>, q=0.0, upper=200.0) -> float`  _function_

> Price a European option under the Heston model.
>
> Args:
>     v0, kappa, theta, xi, rho: Heston parameters (see module docstring).
>     q: continuous dividend yield.
>     upper: truncation of the Fourier integral (200 is ample for typical
>         parameters; raise for very long maturities or large xi).
>
> Returns the option price. Puts are obtained from put-call parity.

### `heston_smile(S, strikes, t, r, v0, kappa, theta, xi, rho, q=0.0)`  _function_

> The Black-Scholes implied-vol smile a Heston model produces.
>
> Prices a European call at each strike under Heston, then inverts each price
> to its Black-Scholes implied volatility, returning ``(log_moneyness, vol)``
> pairs sorted by strike (``log_moneyness = ln(K / F)`` on the forward
> ``F = S e^{(r-q)t}``). This exposes the skew/smile the stochastic-vol
> parameters imply; a negative ``rho`` gives the usual downward equity skew.

## heston_calib

### `calibrate_heston(S, r, quotes, q=0.0, initial=None, feller_weight=0.0, max_iter=6000)`  _function_

> Fit the five Heston parameters to an implied-vol surface.
>
> Args:
>     quotes: iterable of ``(expiry, strike, market_vol)`` Black implied-vol
>         points.
>     initial: optional ``(v0, kappa, theta, xi, rho)`` seed; an ATM-variance
>         based guess is used otherwise.
>     feller_weight: if > 0, adds ``feller_weight * max(0, xi^2 - 2 kappa
>         theta)^2`` to the objective, nudging the fit toward the Feller
>         condition (a strictly positive variance process).
>
> Returns ``(params, rmse)`` -- the fitted 5-tuple and the root-mean-square
> implied-vol error over the quotes.

## heston_mc

### `heston_cv_mc(S, K, t, r, v0, kappa, theta, xi, rho, option_type=<OptionType.CALL: 'call'>, q=0.0, n_steps=100, n_paths=50000, antithetic=True, seed=None, gamma1=0.5) -> quantforge.montecarlo.MCResult`  _function_

> Heston QE Monte Carlo with the underlying as a control variate.
>
> Andersen's QE step is martingale-corrected, so the discounted terminal spot
> ``Y = e^{-r t} S_T`` has the *known* mean ``E[Y] = S0 e^{-q t}`` (the
> discounted forward). ``Y`` is strongly correlated with the option payoff, so
> the controlled estimator ``X - beta (Y - E[Y])`` with the regression-optimal
> ``beta = Cov(X, Y) / Var(Y)`` sharply cuts the standard error at no bias.
> Everything else matches :func:`heston_qe_mc` (same QE variance step,
> ``K0..K4`` asset constants, antithetic draws).
>
> Cross-checks the Fourier :func:`quantforge.heston_price` and reports a
> standard error well below :func:`heston_qe_mc` at equal path count.

### `heston_mc_greeks(S, K, t, r, v0, kappa, theta, xi, rho, option_type=<OptionType.CALL: 'call'>, q=0.0, n_steps=100, n_paths=100000, antithetic=True, seed=None)`  _function_

> Heston Greeks by common-random-number finite differences on the QE MC.
>
> Repricing at bumped inputs with the *same* seed makes the two simulations
> share their random draws, so the bumped price difference is dominated by the
> genuine sensitivity rather than Monte Carlo noise -- far lower variance than
> independent-sample bumps. Returns a dict with ``price`` and
>
>     delta      = dV/dS0
>     gamma      = d2V/dS0^2
>     vega_v0    = dV/dv0        (initial-variance sensitivity)
>     vega_theta = dV/dtheta     (long-variance sensitivity)
>     volvol     = dV/dxi        (vol-of-vol sensitivity)
>     rho_sens   = dV/drho       (spot/vol correlation sensitivity)
>
> The variance-parameter Greeks are the ones that matter for a stochastic-vol
> book; ``vega_v0`` is the closest analogue of Black-Scholes vega. Each is a
> central difference with a relative bump; ``gamma`` reuses the delta re-prices.
> Cross-checks a finite difference of the exact Fourier
> :func:`quantforge.heston_price`.

### `heston_pathwise_delta(S, K, t, r, v0, kappa, theta, xi, rho, option_type=<OptionType.CALL: 'call'>, q=0.0, n_steps=100, n_paths=50000, antithetic=True, seed=None) -> quantforge.montecarlo.MCResult`  _function_

> Heston delta by the pathwise method (spot enters multiplicatively).
>
> In the QE simulation the initial spot appears only as the additive constant
> ``ln S0`` in the terminal log-price, so ``S_T = S0 e^Y`` with ``Y``
> independent of ``S0``. The pathwise delta is therefore exact and simple:
>
>     delta = e^{-r t} E[ 1_{S_T > K} S_T / S0 ]   (put: -1_{S_T < K}),
>
> reusing the same Andersen-QE variance/asset scheme as
> :func:`heston_qe_mc`. Cross-checks a finite-difference bump of the price.

### `heston_qe_mc(S, K, t, r, v0, kappa, theta, xi, rho, option_type=<OptionType.CALL: 'call'>, q=0.0, n_steps=100, n_paths=50000, antithetic=True, seed=None, gamma1=0.5) -> quantforge.montecarlo.MCResult`  _function_

> Price a European option under Heston by Andersen's QE Monte Carlo.
>
> Args:
>     v0, kappa, theta, xi, rho: Heston parameters (initial variance, mean
>         reversion speed, long variance, vol-of-vol, spot/vol correlation).
>     q: continuous dividend yield.
>     n_steps: time steps (the QE variance update is exact in its moments, so
>         the residual bias is only in the asset integral; ~50-200 is plenty).
>     gamma1: weight on the left variance endpoint in the variance integral
>         (``gamma1 = 0.5`` is the central discretisation; ``gamma2`` is set
>         to ``1 - gamma1``).
>
> Returns an :class:`MCResult`. Puts come from simulating the same paths and
> taking the put payoff (parity holds path-by-path at the terminal spot).

## holee

### `holee_bond_greeks(r0, t, theta, sigma)`  _function_

> Exact rate sensitivities of a Ho-Lee zero-coupon bond.
>
> ``P = exp(-r0 t - ...)`` is linear in ``r0`` inside the exponent, so
> ``rho_r = dP/dr0 = -t P``, ``gamma_r = d2P/dr0^2 = t^2 P``, and the rate
> ``duration`` is exactly ``t`` (a Ho-Lee bond has duration equal to its
> maturity), with ``convexity = t^2``. Returns a dict with ``price``,
> ``rho_r``, ``gamma_r``, ``duration``, ``convexity``.

### `holee_bond_option(r0, t_option, t_bond, strike, theta, sigma, is_call=True)`  _function_

> European option on a Ho-Lee zero-coupon bond (exact Black-style).
>
> ``ln P(t_option, t_bond)`` is Gaussian, so the option is a Black formula on
> the forward bond ``P(0, t_bond) / P(0, t_option)`` with bond volatility
> ``sigma_p = sigma * (t_bond - t_option) * sqrt(t_option)`` (the Ho-Lee
> ``B(tau) = tau`` gives the linear maturity factor).

### `holee_cap(r0, dates, strike, theta, sigma, notional=1.0)`  _function_

> Ho-Lee cap: strip of caplets over successive ``dates`` (increasing times).

### `holee_caplet(r0, reset, pay, strike, theta, sigma, notional=1.0)`  _function_

> Caplet on ``[reset, pay]`` under Ho-Lee via the bond-put identity.

### `holee_coupon_bond_option(r0, t_option, cashflows, strike, theta, sigma, is_call=True)`  _function_

> European option on a coupon bond under Ho-Lee (Jamshidian decomposition).
>
> ``cashflows`` is ``[(t_i, c_i), ...]`` with ``t_i > t_option``. The Ho-Lee
> bond is monotone decreasing in ``r0``, so Jamshidian's trick applies: solve
> for the critical rate ``r*`` where the coupon bond's value at expiry equals
> ``strike``, then sum the ``c_i``-weighted zero-coupon-bond options struck at
> ``K_i = P(t_option, t_i | r*)``. Exact.

### `holee_expected_rate(r0, t, theta, sigma=0.0)`  _function_

> Expected Ho-Lee short rate ``E[r_t] = r0 + theta t``.
>
> With constant drift and no mean reversion the rate is
> ``r_t = r0 + theta t + sigma W_t``, so the mean drifts linearly at rate
> ``theta`` (``sigma`` does not enter the mean).

### `holee_floor(r0, dates, strike, theta, sigma, notional=1.0)`  _function_

> Ho-Lee floor: strip of floorlets over successive ``dates``.

### `holee_floorlet(r0, reset, pay, strike, theta, sigma, notional=1.0)`  _function_

> Floorlet on ``[reset, pay]`` under Ho-Lee via the bond-call identity.

### `holee_rate_variance(t, sigma)`  _function_

> Variance of the Ho-Lee short rate ``Var[r_t] = sigma^2 t``.
>
> The rate is a drifted Brownian motion, so its variance grows linearly and
> without bound -- there is no stationary distribution (no mean reversion).

### `holee_swaption(r0, expiry, pay_times, fixed_rate, theta, sigma, payer=True, notional=1.0)`  _function_

> European swaption under Ho-Lee via the coupon-bond-option identity (exact).
>
> A payer swaption is a put on the fixed-leg coupon bond struck at the
> notional; a receiver is a call. Priced by :func:`holee_coupon_bond_option`.
> ``pay_times`` are the fixed-leg payment dates (all ``> expiry``); accruals
> are the gaps, the first measured from ``expiry``.

### `holee_zero_coupon_bond(r0, t, theta, sigma)`  _function_

> Ho-Lee zero-coupon bond price P(0, t) with constant drift ``theta``.
>
> ``P = exp(-r0 t - 0.5 theta t^2 + sigma^2 t^3 / 6)``.

### `holee_zero_coupon_yield(r0, t, theta, sigma)`  _function_

> Continuously-compounded yield of the Ho-Lee zero-coupon bond to ``t``.
>
> ``y(t) = r0 + 0.5 theta t - sigma^2 t^2 / 6`` (linear-in-t drift, quadratic
> convexity pull-down).

## implied

### `implied_vol_smile(strikes, prices, S, t, r, option_type=<OptionType.CALL: 'call'>, b=None, forward=None)`  _function_

> Invert a whole chain of quotes to an implied-vol smile in one call.
>
> Args:
>     strikes, prices: equal-length option-quote arrays at one expiry.
>     forward: optional forward for the log-moneyness output; defaults to the
>         carry-implied forward ``S e^{b t}``.
>
> Returns a list of ``(log_moneyness, implied_vol)`` pairs sorted by strike,
> skipping any quote outside the no-arbitrage band (those cannot be inverted).
> ``log_moneyness = ln(K / forward)``, the standard smile x-axis.

### `implied_volatility(target_price, S, K, t, r, option_type=<OptionType.CALL: 'call'>, b=None, tol=1e-08, max_iter=100, lo=1e-09, hi=10.0)`  _function_

> Solve for the volatility that reproduces ``target_price``.
>
> Returns the implied vol, or raises ValueError if the quote is outside the
> no-arbitrage band (no finite vol can produce it).

## income

### `IncomeMetrics(premium: float, static_yield: float, annualized_yield: float, if_assigned_return: float, breakeven: float) -> None`  _class_

> IncomeMetrics(premium: float, static_yield: float, annualized_yield: float, if_assigned_return: float, breakeven: float)

### `cash_secured_put(S, K, t, r, sigma, b=None, premium=None) -> quantforge.income.IncomeMetrics`  _function_

> Cash-secured put: short a put struck at ``K``, holding ``K`` cash.
>
> Args:
>     premium: option premium; if None, uses the BSM put value.
>
> Yields are relative to the secured cash ``K``. If assigned (S_T <= K) the
> trader buys the stock at ``K`` net of premium, so the effective purchase and
> breakeven price is ``K - premium``.

### `covered_call(S, K, t, r, sigma, b=None, premium=None) -> quantforge.income.IncomeMetrics`  _function_

> Covered call: long stock at ``S``, short a call struck at ``K``.
>
> Args:
>     premium: option premium; if None, uses the BSM call value.
>
> Yields are relative to the stock capital ``S``. If assigned (S_T >= K) the
> return is the capped gain to the strike plus the premium; the breakeven is
> ``S - premium`` (the stock can fall by the premium before a loss).

## inflation

### `apply_seasonality(deseasonalized_index, seasonal_factor)`  _function_

> Add the seasonal component back: ``deseasonalized_index * seasonal_factor``.
>
> Inverse of :func:`deseasonalize`. With a :func:`normalize_seasonal_factors`
> factor this raises or lowers the observed index around its trend without
> shifting the annual average.

### `breakeven_inflation(nominal_yield, real_yield) -> float`  _function_

> Breakeven inflation implied by a nominal and a real yield (Fisher).
>
> ``(1 + nominal)/(1 + real) - 1`` -- the inflation rate at which a nominal and
> an inflation-linked bond of the same maturity have equal return. The market's
> inflation expectation (plus risk premium).

### `deflation_floor_value(face, index_settle, index_base) -> float`  _function_

> Intrinsic value of the deflation floor: floored redemption minus unfloored.
>
> ``face * (max(ratio, 1) - ratio)`` -- zero when the index has risen (the floor
> is out of the money), positive under net deflation. The realized payoff of the
> embedded floor option, ignoring optionality/time value.

### `deflation_floored_redemption(face, index_settle, index_base) -> float`  _function_

> TIPS-style redemption with the deflation floor: principal never below par.
>
> Real (TIPS) principal redeems at ``face * max(index_ratio, 1)`` -- the index
> ratio inflates the principal in inflation, but a cumulative deflation over the
> bond's life cannot pull the redemption below the original face. Equals
> :func:`inflation_adjusted_principal` whenever the index has risen since issue
> (ratio >= 1), and is floored to ``face`` otherwise.

### `deseasonalize(observed_index, seasonal_factor)`  _function_

> Strip the seasonal component: ``observed_index / seasonal_factor``.
>
> Inverse of :func:`apply_seasonality`; recovers the trend index used for
> projecting forward fixings free of the within-year seasonal pattern.

### `fisher_nominal_rate(real, inflation) -> float`  _function_

> Exact Fisher nominal rate ``(1 + real)(1 + inflation) - 1``.

### `fisher_real_rate(nominal, inflation) -> float`  _function_

> Exact Fisher real rate ``(1 + nominal)/(1 + inflation) - 1``.
>
> The rate that, compounded with inflation, reproduces the nominal rate. For
> small rates it is approximately ``nominal - inflation``.

### `forward_inflation_rate(index_start, index_end, t_start, t_end)`  _function_

> Annualized forward inflation between two curve horizons.
>
> ``(I_end / I_start)^(1/(t_end - t_start)) - 1`` -- the constant annual rate
> linking two projected index levels. Chains with the near-leg rate so that
> ``(1 + spot)^t_start (1 + fwd)^(t_end - t_start) = (1 + spot_end)^t_end`` (the
> no-arbitrage forward/spot relation tested against).

### `index_ratio(index_settle, index_base) -> float`  _function_

> Index ratio ``CPI_settle / CPI_base`` used to inflate the principal.
>
> Above 1 when the price index has risen since issue. Both indices must be
> positive.

### `index_ratio_interpolated(cpi_month_start, cpi_next_month, day, days_in_month, cpi_base)`  _function_

> Index ratio using the daily-interpolated :func:`reference_cpi`.
>
> ``reference_cpi(...) / cpi_base`` -- the ratio a linker actually applies to
> its principal on a mid-month settlement, versus the month-boundary
> :func:`index_ratio` which ignores intra-month accrual.

### `inflation_adjusted_principal(face, index_settle, index_base) -> float`  _function_

> Inflation-adjusted principal ``face * index_ratio`` (the linker notional).

### `inflation_curve_from_zc_swaps(index_base, tenors, zc_rates)`  _function_

> Projected index levels implied by a strip of zero-coupon swap rates.
>
> A ZC inflation swap of maturity ``T`` with fair rate ``k_T`` pins the forward
> index to ``I_0 * (1 + k_T)^T`` (the :func:`zc_inflation_swap_rate` identity).
> Given quotes ``(tenors, zc_rates)`` this returns the matching forward index
> levels ``[I_0 (1 + k_T)^T for T in tenors]`` -- the market-implied inflation
> curve, expressed as projected index fixings. By construction reinverting each
> level through :func:`zc_inflation_swap_rate` recovers the input ``zc_rates``.

### `linker_price(real_cashflows, real_yield, index_settle, index_base) -> float`  _function_

> Dirty price of an inflation-linked bond off real cashflows.
>
> ``real_cashflows`` is ``[(t, real_amount), ...]`` in constant (issue-date)
> money. Each flow is discounted at the continuously-compounded ``real_yield``
> and then the whole bond is inflated by the settlement index ratio:
>
>     price = (index_settle / index_base) * sum_i real_amount_i e^{-r t_i}
>
> Because the index ratio multiplies every flow, the price is degree-one
> homogeneous in it -- stripping the ratio recovers a standard real-yield bond
> price (the invariant tested against :mod:`quantforge.bondmath`).

### `linker_real_convexity(real_cashflows, real_yield) -> float`  _function_

> Convexity of a linker w.r.t. its real yield ``1/P d2P/dr2``.
>
> PV-weighted average of squared cashflow time on the real cashflows; the index
> ratio cancels, matching :func:`quantforge.convexity`.

### `linker_real_duration(real_cashflows, real_yield) -> float`  _function_

> Modified duration of a linker w.r.t. its real yield (years).
>
> ``-1/P dP/dr``. Since the settlement index ratio multiplies the whole price it
> cancels in the fractional sensitivity, so the real duration is the PV-weighted
> average cashflow time of the *real* cashflows -- identical to the standard
> :func:`quantforge.modified_duration` on those flows, independent of the index
> level.

### `linker_real_dv01(real_cashflows, real_yield, index_settle, index_base) -> float`  _function_

> Dollar value of a 1bp real-yield rise for a linker (negative).
>
> ``dP/dr * 1e-4 = -duration * price * 1e-4`` on the inflated (dirty) price, so
> unlike the fractional duration this DOES scale with the index ratio.

### `linker_real_yield(real_cashflows, price, index_settle, index_base, tol=1e-10, max_iter=100) -> float`  _function_

> Continuously-compounded real yield reproducing a linker ``price``.
>
> Deflates the quoted price by the index ratio and solves the standard real-
> cashflow bond yield by bisection (price is monotone decreasing in the yield).
> Inverse of :func:`linker_price`.

### `nominal_zero_curve(real_zeros, breakevens)`  _function_

> Per-tenor nominal zero rates from real zeros and breakeven inflation.
>
> Fisher forward direction, ``nominal = (1 + real)(1 + breakeven) - 1``
> (:func:`fisher_nominal_rate`), reconstructing the nominal curve. Inverse of
> :func:`real_zero_curve`, so composing the two is the identity per tenor.

### `normalize_seasonal_factors(raw_factors)`  _function_

> Scale 12 monthly seasonal factors to a geometric mean of one.
>
> Seasonal adjustment must not change the trend level over a full year, so the
> monthly factors are normalized by their geometric mean:
> ``f_i / (prod f_j)^(1/12)``. The result multiplies to one across the year, so
> compounding all twelve leaves the annual index unchanged. Raw factors must be
> positive.

### `real_discount_factor(nominal_df, index_ratio_t)`  _function_

> Real discount factor from a nominal one and the period index growth.
>
> The real (inflation-adjusted) discount factor grows the nominal by the
> realized/projected index ratio over the period: ``nominal_df * index_ratio_t``.
> A real cashflow discounted at the real DF equals its inflated nominal cashflow
> discounted at the nominal DF -- the identity linking the two measures.

### `real_from_breakeven(nominal_yield, breakeven) -> float`  _function_

> Real yield implied by a nominal yield and a breakeven inflation rate.
>
> ``(1 + nominal)/(1 + breakeven) - 1`` -- inverse of
> :func:`breakeven_inflation`.

### `real_zero_curve(nominal_zeros, breakevens)`  _function_

> Per-tenor real zero rates from nominal zeros and breakeven inflation.
>
> Applies the Fisher relation tenor by tenor,
> ``real = (1 + nominal)/(1 + breakeven) - 1`` (:func:`real_from_breakeven`),
> turning a nominal zero curve and a breakeven-inflation curve into the implied
> real zero curve. Inverse of :func:`nominal_zero_curve`.

### `reference_cpi(cpi_month_start, cpi_next_month, day, days_in_month)`  _function_

> Daily reference index by linear interpolation between two monthly fixings.
>
> Inflation-linked bonds accrue off a *reference index* that interpolates
> linearly within the month between the anchor CPI for the first of the month
> and the first of the next month (the standard linker daily-indexation rule,
> applied to the lagged CPIs). For settlement on the ``day``-th of a month with
> ``days_in_month`` days:
>
>     ref = cpi_month_start + (day - 1)/days_in_month
>               * (cpi_next_month - cpi_month_start)
>
> Equals ``cpi_month_start`` on the 1st and approaches ``cpi_next_month`` at
> month end.

### `yoy_cap_implied_vol(price, forward_rates, strike, expiries, discount_factors, notional=1.0, is_cap=True, tol=1e-10, max_iter=100)`  _function_

> Flat Black vol reproducing a year-on-year cap/floor ``price``.
>
> Bisection on the common ``sigma`` (cap value is monotone increasing in vol),
> inverting :func:`yoy_cap_price`. The price must lie between the zero-vol
> intrinsic and the vol -> infinity bound.

### `yoy_cap_price(forward_rates, strike, expiries, sigma, discount_factors, notional=1.0, is_cap=True)`  _function_

> Year-on-year inflation cap/floor: a strip of :func:`yoy_caplet_price`.
>
> Sums the Black-76 caplet (or floorlet) values across each YoY period, one per
> ``(forward_rate, expiry, discount_factor)`` triple, at a common ``strike`` and
> flat ``sigma``. A single-period strip equals the caplet; cap minus floor
> telescopes to ``sum_i DF_i * N * (F_i - K)``.

### `yoy_caplet_implied_normal_vol(price, forward_rate, strike, expiry, discount_factor, notional=1.0, is_cap=True, tol=1e-12, max_iter=100)`  _function_

> Normal (Bachelier) vol reproducing a YoY caplet/floorlet ``price``.
>
> Bisection on ``sigma`` (price is monotone increasing in normal vol), inverting
> :func:`yoy_caplet_price_normal`. Works for any real forward/strike, including
> negative inflation forwards.

### `yoy_caplet_price(forward_rate, strike, expiry, sigma, discount_factor, notional=1.0, is_cap=True)`  _function_

> Black-76 price of a year-on-year inflation cap/floor let.
>
> Prices a single YoY period whose payoff is ``max(YoY - K, 0)`` (caplet) or
> ``max(K - YoY, 0)`` (floorlet), with the year-on-year inflation rate modelled
> as lognormal around its ``forward_rate`` (from :func:`forward_inflation_rate`)
> with volatility ``sigma`` to ``expiry``. Standard Black-76:
>
>     caplet  = DF * N * [F Phi(d1) - K Phi(d2)]
>     floorlet= DF * N * [K Phi(-d2) - F Phi(-d1)]
>     d1,2    = (ln(F/K) +/- 0.5 sigma^2 T) / (sigma sqrt(T))
>
> Requires positive ``forward_rate`` and ``strike`` (lognormal support). At zero
> vol it collapses to the discounted intrinsic ``DF*N*max(F-K,0)`` (cap).

### `yoy_caplet_price_normal(forward_rate, strike, expiry, sigma, discount_factor, notional=1.0, is_cap=True)`  _function_

> Bachelier (normal-model) price of a year-on-year inflation cap/floor let.
>
> Models the YoY rate as *arithmetic* Brownian motion around its forward, so it
> admits zero and negative inflation (where the lognormal
> :func:`yoy_caplet_price` cannot price). Bachelier:
>
>     caplet   = DF * N * [(F - K) Phi(d) + sigma sqrt(T) phi(d)]
>     floorlet = DF * N * [(K - F) Phi(-d) + sigma sqrt(T) phi(d)]
>     d        = (F - K) / (sigma sqrt(T))
>
> ``sigma`` is a normal (absolute-rate) vol. At zero vol it collapses to the
> discounted intrinsic; the ATM caplet equals ``DF*N*sigma*sqrt(T/(2 pi))``.

### `yoy_inflation_rate(index_prev, index_curr) -> float`  _function_

> Year-on-year inflation ``index_curr / index_prev - 1`` between two fixings.

### `yoy_swap_value(notional, fixed_rate, index_levels, discount_factors, index_prev)`  _function_

> Value of a year-on-year inflation swap off a projected index curve.
>
> Each period ``i`` exchanges the realized year-on-year inflation
> ``I_i / I_{i-1} - 1`` (float, received) for ``fixed_rate`` (paid), on
> ``notional``, discounted by ``discount_factors[i]``. ``index_prev`` is the
> fixing one period before the first ``index_levels`` entry (the base for the
> first YoY ratio). Returns the inflation-receiver's value
>
>     notional * sum_i (I_i/I_{i-1} - 1 - fixed_rate) * DF_i.
>
> Unlike the single-payment ZC swap this pays the annual inflation each period.

### `zc_inflation_swap_rate(index_start, index_end, years) -> float`  _function_

> Fair annualized rate of a zero-coupon inflation swap.
>
> A ZC inflation swap exchanges ``(1 + k)^T - 1`` (fixed) for the realized index
> growth ``I_T / I_0 - 1`` (float) at maturity. The par fixed rate that zeroes
> the swap is the annualized index growth ``(I_T / I_0)^(1/T) - 1``, so that
> ``(1 + k)^T * I_0 == I_T`` (the compounding identity tested against).

### `zc_inflation_swap_value(notional, fixed_rate, index_start, index_end, years, discount_factor=1.0) -> float`  _function_

> Value of the inflation leg minus the fixed leg of a ZC inflation swap.
>
> Inflation-leg receiver's value: ``notional * (I_T/I_0 - (1+k)^T)`` at maturity,
> discounted by ``discount_factor``. Zero at the par :func:`zc_inflation_swap_rate`.

## kim

### `kim_american_call(S, K, t, r, sigma, q=0.0, n_steps=80)`  _function_

> American call price via the put-call symmetry for American options.
>
> A dividend-paying American call maps to an American put by the McDonald-
> Schroder symmetry ``C(S, K, r, q) = P(K, S, q, r)`` (spot<->strike,
> rate<->dividend). With ``q = 0`` the call is never exercised early and this
> returns the European call.

### `kim_american_put(S, K, t, r, sigma, q=0.0, n_steps=80)`  _function_

> American put price via Kim's integral equation.
>
> Solves the early-exercise boundary on an ``n_steps`` time grid, then returns
> the European put plus the early-exercise premium integrated at spot ``S``.
> If ``S`` is at or below the current boundary the option is exercised, so the
> intrinsic value is returned.

### `kim_exercise_boundary(K, t, r, sigma, q=0.0, n_steps=80)`  _function_

> Return the American-put early-exercise boundary ``B(t_i)`` on the grid.
>
> ``B[i]`` is the critical spot at time ``i * (t / n_steps)`` below which
> immediate exercise is optimal; ``B[n_steps]`` is the expiry value.

### `kim_put_greeks(S, K, t, r, sigma, q=0.0, n_steps=80)`  _function_

> Delta, gamma, theta of a Kim American put, reusing one boundary solve.
>
> The early-exercise boundary is spot-independent, so it is solved once (the
> expensive step) and the spot/time bumps only re-run the cheap European-plus-
> premium evaluation. Delta and gamma come from central differences in ``S``
> on that fixed boundary; theta from a maturity bump (which does re-solve the
> boundary). Returns ``{price, delta, gamma, theta}``.

## kou

### `kou_greeks(S, K, t, r, sigma, lam, p, eta1, eta2, option_type=<OptionType.CALL: 'call'>, q=0.0)`  _function_

> Greeks of a Kou double-exponential jump-diffusion option by FD.
>
> Central finite differences of :func:`kou_price` for ``delta`` (dV/dS),
> ``gamma`` (d2V/dS2), ``vega`` (dV/dsigma, the *diffusion*-vol sensitivity),
> and ``theta`` (calendar decay). At ``lam = 0`` (no jumps) these reduce to the
> vanilla Black-Scholes Greeks. Returns a dict with ``price`` and those fields.

### `kou_price(S, K, t, r, sigma, lam, p, eta1, eta2, option_type=<OptionType.CALL: 'call'>, q=0.0, upper=200.0) -> float`  _function_

> Price a European option under the Kou double-exponential jump-diffusion.
>
> Args:
>     sigma: diffusion volatility.
>     lam: jump intensity (expected jumps per year, >= 0).
>     p: probability a jump is upward (in [0, 1]).
>     eta1: up-jump tail rate (must be > 1 so E[e^Y] is finite).
>     eta2: down-jump tail rate (> 0).
>     q: continuous dividend yield.
>
> ``lam = 0`` recovers Black-Scholes. Puts use put-call parity.

### `kou_smile(S, strikes, t, r, sigma, lam, p, eta1, eta2, q=0.0)`  _function_

> Black-Scholes implied-vol smile the Kou model produces.
>
> Prices a call at each strike and inverts to a Black-Scholes implied vol,
> returning ``(log_moneyness, vol)`` pairs sorted by strike on the forward
> ``F = S e^{(r-q) t}``. An asymmetric jump distribution (``eta1 != eta2`` or
> ``p != 1/2``) tilts the smile into a skew.

## ldi

### `funded_ratio_return(asset_return, liability_return, funding_ratio_)`  _function_

> Change in funding ratio from asset and liability returns.
>
> ``FR_new / FR_old - 1 = (1 + asset_return) / (1 + liability_return) - 1`` (the
> funding ratio's own return is independent of its level). Positive when assets
> outperform liabilities. ``funding_ratio_`` is accepted for context but the
> fractional change does not depend on it.

### `funding_ratio(assets, liabilities)`  _function_

> Funding ratio ``assets / liabilities`` (above 1 = surplus).

### `hedge_ratio(asset_duration, asset_value, liability_duration_, liability_value)`  _function_

> Fraction of the liability dollar-duration hedged by the assets.
>
> ``(asset_duration * asset_value) / (liability_duration * liability_value)`` --
> the ratio of asset to liability dollar duration (DV01). One means the surplus
> is immune to a parallel rate move; below one leaves residual liability
> interest-rate risk.

### `liability_convexity(cashflows, discount_rate)`  _function_

> Convexity of the liability stream ``sum t^2 PV_i / sum PV_i``.
>
> The second-order interest-rate sensitivity; matching it in addition to
> duration gives the surplus protection against larger, non-parallel rate moves
> (Redington immunization's second condition).

### `liability_duration(cashflows, discount_rate)`  _function_

> Macaulay duration of the liability stream (years).
>
> PV-weighted average cashflow time ``sum t_i PV_i / sum PV_i``. The interest-
> rate sensitivity the asset portfolio must match to immunize the surplus.

### `liability_pv(cashflows, discount_rate)`  _function_

> Present value of a liability stream ``[(t, amount), ...]``.
>
> Continuously-compounded discounting ``sum_i CF_i e^{-r t_i}``.

### `required_hedge_duration(asset_value, liability_duration_, liability_value)`  _function_

> Asset duration that fully immunizes the surplus (hedge ratio = 1).
>
> ``liability_duration * liability_value / asset_value`` -- the duration the
> asset portfolio must carry so its dollar duration matches the liability's.

### `surplus(assets, liabilities)`  _function_

> Plan surplus (deficit if negative): ``assets - liabilities``.

### `surplus_at_risk(assets, liabilities, surplus_volatility, confidence=0.95, horizon=1.0)`  _function_

> Surplus-at-risk: the worst surplus loss at a confidence over a horizon.
>
> ``z * surplus_volatility * sqrt(horizon) * liabilities`` where ``z =
> Phi^{-1}(confidence)`` and ``surplus_volatility`` is the funded-status (surplus/
> liabilities) return volatility. A one-sided downside measure (positive number =
> potential shortfall), analogous to VaR for the plan surplus.

### `surplus_change_under_shock(assets, asset_duration, asset_convexity, liabilities, liability_duration_, liability_convexity_, rate_shock)`  _function_

> Second-order surplus change under a parallel rate shock ``dy``.
>
> Uses the duration-convexity expansion on each side:
>
>     dV = V * (-D dy + 0.5 C dy^2),
>     d(surplus) = dAssets - dLiabilities.
>
> A duration-matched but convexity-mismatched book still moves at second order;
> matching both leaves the surplus (nearly) unchanged.

## leisen_reimer

### `leisen_reimer_american_accel(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None, steps=101) -> float`  _function_

> Richardson-extrapolated Leisen-Reimer American price (Broadie-Detemple).
>
> American LR convergence is only ``O(1/n)`` (the smooth-payoff assumption
> behind the Peizer-Pratt inversion breaks at the early-exercise boundary),
> unlike the ``O(1/n^2)`` European case. Broadie & Detemple (1996) cancel that
> leading ``1/n`` term with a two-point Richardson extrapolation between an
> ``n``-step and a ``2n``-step tree:
>
>     V_ext = 2 * V(2n) - V(n).
>
> For the same work this is several times more accurate than a single tree, so
> a moderate ``steps`` reaches four-figure accuracy. ``b`` is the cost of carry
> (dividend yield ``q`` via ``b = r - q``).

### `leisen_reimer_greeks(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None, steps=101, american=False)`  _function_

> Delta and gamma from Leisen-Reimer tree nodes, plus FD vega/theta.
>
> Delta and gamma are read directly off the first two time steps of the tree
> (no extra pricing passes), while vega and theta use small central
> differences. Returns a dict with price, delta, gamma, vega and theta.

### `leisen_reimer_price(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None, steps=101, american=False) -> float`  _function_

> Price an option with a Leisen-Reimer binomial tree.
>
> Args:
>     steps: number of time steps; forced to the next odd integer so the tree
>         straddles the strike (the source of the fast, monotone convergence).
>     american: if ``True``, apply an early-exercise check at each node;
>         otherwise price the European payoff.
>     b: cost of carry (defaults to ``r``); dividend yield ``q`` enters as
>         ``b = r - q``.

## levycalib

### `calibrate_levy_smile(model, S, t, r, strikes, market_vols, q=0.0, cm_alpha=1.5, max_iter=4000)`  _function_

> Fit a Levy model to a one-expiry market smile by least squares on vol.
>
> Args:
>     model: one of ``"vg"``, ``"nig"``, ``"meixner"``, ``"cgmy"``.
>     strikes, market_vols: matching sequences of strikes and Black-Scholes
>         implied vols at expiry ``t``.
>     cm_alpha: Carr-Madan damping used when pricing each candidate.
>
> Returns ``(params, rmse)`` where ``params`` is the fitted raw-parameter tuple
> for the chosen model and ``rmse`` is the root-mean-square implied-vol error.

### `levy_psi(model, params)`  _function_

> Return the characteristic exponent ``psi(u)`` for a named Levy model.
>
> ``model`` is one of ``"vg"``, ``"nig"``, ``"meixner"``, ``"cgmy"`` and
> ``params`` its raw parameter tuple (the same shape :func:`calibrate_levy_smile`
> returns). Handy for pricing or surface-building a model outside calibration.

## levysurface

### `LevyCalendarViolation(t_short: float, t_long: float, k: float, w_short: float, w_long: float) -> None`  _class_

> LevyCalendarViolation(t_short: float, t_long: float, k: float, w_short: float, w_long: float)

### `LevySurface(model, params, S, r, q=0.0, cm_alpha=1.5)`  _class_

> Implied-vol surface generated by one exponential-Levy parameter set.

## localvol

### `dupire_local_vol(call_fn: Callable[[float, float], float], K: float, T: float, r: float, q: float = 0.0, dK: float = None, dT: float = None) -> float`  _function_

> Dupire local volatility at strike ``K`` and maturity ``T``.
>
> Args:
>     call_fn: ``C(K, T)`` returning the European call price for strike K and
>         maturity T (both positive). Must be evaluable in a neighborhood of
>         (K, T) for the finite differences.
>     r, q: risk-free rate and continuous dividend yield.
>     dK, dT: finite-difference bumps; default to small fractions of K and T.
>
> Returns the local volatility (not variance). Raises if the local variance
> comes out non-positive (a sign of an arbitrageable / too-noisy surface).

### `local_vol_from_implied(implied_vol_fn: Callable[[float, float], float], S: float, K: float, T: float, r: float, q: float = 0.0, dK: float = None, dT: float = None) -> float`  _function_

> Dupire local vol from an implied-vol surface ``sigma_imp(K, T)``.
>
> Wraps :func:`dupire_local_vol` by turning the implied-vol surface into a
> call-price surface with the Black-Scholes-Merton formula (carry ``b = r-q``).

### `sabr_local_vol(S, K, T, r, alpha, beta, rho, nu, q=0.0, dK=None)`  _function_

> Dupire local volatility of a single SABR smile at expiry ``T``.
>
> Builds the SABR implied-vol smile (Hagan) at maturity ``T`` on the forward
> ``F = S e^{(r-q)T}`` and feeds it into the Dupire formula. Only the strike
> derivatives are needed at a fixed expiry, so this reads the SABR smile in
> strike and holds ``T`` fixed for the maturity bump (a flat local term
> structure across the single slice).

## lookback

### `discrete_fixed_strike_lookback(S, K, t, r, sigma, n_fixings, option_type=<OptionType.CALL: 'call'>, b=None) -> float`  _function_

> Discretely-monitored fixed-strike lookback (Broadie-Glasserman-Kou 1999).
>
> The realized extreme is sampled at ``n_fixings`` equally-spaced dates rather
> than continuously, which lowers a call-on-max and raises a put-on-min versus
> continuous monitoring. Broadie-Glasserman-Kou give an asymptotic continuity
> correction: shift the *spot* fed to the continuous
> :func:`fixed_strike_lookback` by ``exp(-/+ beta sigma sqrt(dt))`` (down for a
> call on the max, up for a put on the min), with ``beta ~ 0.5826`` and
> ``dt = t / n_fixings``. As ``n_fixings -> infinity`` the shift vanishes and
> the price converges to the continuous lookback.
>
> Accurate to a few tenths of a percent for ``n_fixings`` of ~50 or more; the
> correction is asymptotic, so coarse monitoring (a handful of dates) carries a
> larger error.

### `discrete_fixed_strike_lookback_greeks(S, K, t, r, sigma, n_fixings, option_type=<OptionType.CALL: 'call'>, b=None)`  _function_

> Greeks of a discretely-monitored fixed-strike lookback by central finite
> differences of :func:`discrete_fixed_strike_lookback`: ``delta`` (dV/dS),
> ``gamma`` (d2V/dS2), ``vega`` (dV/dsigma), ``theta`` (calendar decay). The
> number of monitoring dates ``n_fixings`` is held fixed. Returns a dict with
> ``price`` and those fields.

### `fixed_strike_lookback(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, s_extreme=None, b=None) -> float`  _function_

> Fixed-strike lookback (Conze-Viswanathan).
>
> Call pays ``max(S_max - K, 0)``; put pays ``max(K - S_min, 0)``.
>
> Args:
>     s_extreme: running maximum (call) or minimum (put) so far. Defaults to
>         the current spot.

### `floating_strike_lookback(S, t, r, sigma, option_type=<OptionType.CALL: 'call'>, s_extreme=None, b=None) -> float`  _function_

> Floating-strike lookback (Goldman-Sosin-Gatto).
>
> Args:
>     s_extreme: running minimum (for a call) or maximum (for a put) observed
>         so far. Defaults to the current spot (inception).
>     b: cost of carry (defaults to r).
>
> Call payoff: ``S_T - S_min``. Put payoff: ``S_max - S_T``.

### `lookback_greeks(S, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None, kind='floating', K=None, s_extreme=None)`  _function_

> Greeks of a lookback option by central finite differences.
>
> ``kind`` selects the closed form: ``"floating"``
> (:func:`floating_strike_lookback`) or ``"fixed"``
> (:func:`fixed_strike_lookback`, which needs ``K``). Returns a dict with
> delta, gamma, vega, and theta (calendar, per year). ``s_extreme`` (the
> running min/max) defaults to the current spot.

## lsm

### `bermudan_basket_lsm(S1, S2, w1, w2, K, t, r, sigma1, sigma2, rho, q1=0.0, q2=0.0, option_type=<OptionType.CALL: 'call'>, n_steps=50, n_paths=20000, seed=None) -> float`  _function_

> American basket option on ``w1 S1 + w2 S2`` by Longstaff-Schwartz.
>
> Prices ``max(w1 S1_T + w2 S2_T - K, 0)`` (call) or the put, exercisable at
> ``n_steps`` equally-spaced dates. Two correlated GBMs are simulated and the
> continuation value is regressed on a quadratic basis in both spots plus the
> basket ``B = w1 S1 + w2 S2``: ``{1, S1, S2, S1^2, S2^2, S1 S2, B}`` over the
> in-the-money paths at each date.
>
> Returns the price (in-sample LSM estimate, mildly biased low). It sits at or
> above the European moment-matched :func:`quantforge.basket_option`;
> dividends create an early-exercise premium (and American puts carry one even
> without).

### `bermudan_basket_lsm_greeks(S1, S2, w1, w2, K, t, r, sigma1, sigma2, rho, q1=0.0, q2=0.0, option_type=<OptionType.CALL: 'call'>, n_steps=50, n_paths=40000, seed=None, h_rel=0.01)`  _function_

> Deltas and cross-gamma of an American basket option by common-random bumps.
>
> Reprices :func:`bermudan_basket_lsm` at bumped spots on the *same* seed, so
> the two simulations share their Brownian shocks and the finite differences
> are low-variance. Returns a dict with ``price``, the two spot deltas
> (``delta1`` = dV/dS1, ``delta2`` = dV/dS2), the two own-gammas
> (``gamma1``, ``gamma2``), and the cross-gamma (``cross`` = d2V/dS1 dS2).
>
> For a basket *call* both spot deltas are positive (a higher spot lifts the
> weighted basket); the regression is re-fit at each bump. Deltas are reliable;
> the gammas (second differences over a re-fit regression) are indicative and
> need many paths.

### `bermudan_lsm(S, K, t, r, sigma, option_type=<OptionType.PUT: 'put'>, b=None, n_steps=50, n_paths=20000, degree=3, seed=None) -> float`  _function_

> Price a Bermudan option (exercisable at ``n_steps`` equally-spaced dates).
>
> Args:
>     n_steps: number of exercise opportunities over the life; as it grows the
>         price approaches the continuously-exercisable American value.
>     degree: polynomial degree of the regression basis in spot.
>     b: cost of carry (defaults to r). Dividend yield q enters as b = r - q.
>
> Returns the option price (in-sample LSM estimate, mildly biased low).

### `bermudan_lsm_greeks(S, K, t, r, sigma, option_type=<OptionType.PUT: 'put'>, b=None, n_steps=50, n_paths=40000, degree=3, seed=None, h_rel=0.01)`  _function_

> Delta and gamma of a Bermudan/American LSM price by common-random bumps.
>
> Prices the option at ``S``, ``S(1 +/- h)`` on the *same* random-number
> stream (each call reseeds ``bermudan_lsm`` with the same ``seed``, so the
> Brownian paths coincide up to the spot scaling and the finite differences
> are low-variance). Returns a dict with ``price``, ``delta`` and ``gamma``
> from central differences; ``h_rel`` is the relative spot bump.
>
> Common random numbers make the bump estimator far less noisy than
> independent re-pricing; the LSM regression is re-fit at each bump, which is
> the standard practical scheme. Delta is reliable; gamma (a second difference
> over a re-fit regression) is only indicative and needs many paths.

### `bermudan_lsm_local_vol(S, K, t, r, local_vol_fn, option_type=<OptionType.PUT: 'put'>, q=0.0, n_steps=50, n_paths=20000, degree=3, seed=None) -> float`  _function_

> Bermudan/American option under a local-volatility surface by LSM.
>
> Same Longstaff-Schwartz backward induction as :func:`bermudan_lsm`, but each
> Euler step uses the spot- and time-dependent ``local_vol_fn(S, tau)`` (with
> ``tau`` the elapsed forward time) instead of a constant vol -- so it prices
> early-exercise options directly on a calibrated Dupire / SVI local-vol
> surface. Carry is ``b = r - q``. A flat ``local_vol_fn`` reproduces the
> constant-vol LSM price.

### `bermudan_max_call_lsm(S1, S2, K, t, r, sigma1, sigma2, rho, q1=0.0, q2=0.0, n_steps=50, n_paths=20000, seed=None) -> float`  _function_

> American call on the maximum of two assets by Longstaff-Schwartz.
>
> Prices ``max(max(S1_T, S2_T) - K, 0)`` with early exercise at ``n_steps``
> equally-spaced dates -- the classic two-asset LSM benchmark. Two correlated
> GBMs are simulated (``S2``'s shock is ``rho z1 + sqrt(1-rho^2) z2``) and the
> continuation value is regressed on a quadratic basis in both spots plus the
> running max: ``{1, S1, S2, S1^2, S2^2, S1 S2, max(S1,S2)}``, over the
> in-the-money paths at each date.
>
> Returns the price (in-sample LSM estimate, mildly biased low). It sits at or
> above the European :func:`quantforge.best_of_call_closed`; with dividends the
> gap is the early-exercise premium.

### `bermudan_max_call_lsm_greeks(S1, S2, K, t, r, sigma1, sigma2, rho, q1=0.0, q2=0.0, n_steps=50, n_paths=40000, seed=None, h_rel=0.01)`  _function_

> Deltas and cross-gamma of an American max-call by common-random bumps.
>
> Reprices :func:`bermudan_max_call_lsm` at bumped spots on the *same* seed, so
> the two simulations share their Brownian shocks and the finite differences
> are low-variance. Returns a dict with ``price``, the two spot deltas
> (``delta1`` = dV/dS1, ``delta2`` = dV/dS2), the two own-gammas
> (``gamma1``, ``gamma2``), and the cross-gamma (``cross`` = d2V/dS1 dS2).
>
> The regression is re-fit at each bump (the standard practical scheme). The
> deltas are reliable; the gammas -- second differences over a re-fit
> regression -- are only indicative and need many paths.

### `bermudan_min_put_lsm(S1, S2, K, t, r, sigma1, sigma2, rho, q1=0.0, q2=0.0, n_steps=50, n_paths=20000, seed=None) -> float`  _function_

> American put on the minimum of two assets by Longstaff-Schwartz.
>
> Prices ``max(K - min(S1_T, S2_T), 0)`` with early exercise at ``n_steps``
> equally-spaced dates -- the worst-of protective put, a common structured-note
> hedge. Two correlated GBMs are simulated and the continuation value is
> regressed on a quadratic basis in both spots plus the running min
> ``{1, S1, S2, S1^2, S2^2, S1 S2, min(S1,S2)}`` over the in-the-money paths.
>
> Returns the price (in-sample LSM estimate, mildly biased low). Puts carry
> early-exercise value even without dividends, so it sits above the European
> :func:`quantforge.worst_of_put_closed`.

### `bermudan_min_put_lsm_greeks(S1, S2, K, t, r, sigma1, sigma2, rho, q1=0.0, q2=0.0, n_steps=50, n_paths=40000, seed=None, h_rel=0.01)`  _function_

> Deltas and cross-gamma of an American min-put by common-random bumps.
>
> Reprices :func:`bermudan_min_put_lsm` at bumped spots on the *same* seed, so
> the two simulations share their Brownian shocks and the finite differences
> are low-variance. Returns a dict with ``price``, the two spot deltas
> (``delta1`` = dV/dS1, ``delta2`` = dV/dS2), the two own-gammas
> (``gamma1``, ``gamma2``), and the cross-gamma (``cross`` = d2V/dS1 dS2).
>
> For the worst-of protective put ``max(K - min(S1, S2), 0)`` both spot deltas
> are negative (a higher spot lifts the min, shrinking the put); the LSM
> regression is re-fit at each bump. Deltas are reliable; the gammas (second
> differences over a re-fit regression) are indicative and need many paths.

### `bermudan_spread_lsm(S1, S2, K, t, r, sigma1, sigma2, rho, q1=0.0, q2=0.0, option_type=<OptionType.CALL: 'call'>, n_steps=50, n_paths=20000, seed=None) -> float`  _function_

> American spread option ``max(S1_T - S2_T - K, 0)`` by Longstaff-Schwartz.
>
> Prices an early-exercisable option on the spread ``S1 - S2`` (call) or
> ``K - (S1 - S2)`` (put), exercisable at ``n_steps`` equally-spaced dates. Two
> correlated GBMs are simulated and the continuation value is regressed on a
> quadratic basis in both spots plus the spread ``S1 - S2``:
> ``{1, S1, S2, S1^2, S2^2, S1 S2, S1 - S2}``, over the in-the-money paths at
> each date.
>
> Returns the price (in-sample LSM estimate, mildly biased low). It sits at or
> above the European Kirk :func:`quantforge.spread_option`; dividends on the
> long leg create an early-exercise premium.

### `bermudan_spread_lsm_greeks(S1, S2, K, t, r, sigma1, sigma2, rho, q1=0.0, q2=0.0, option_type=<OptionType.CALL: 'call'>, n_steps=50, n_paths=40000, seed=None, h_rel=0.01)`  _function_

> Deltas and cross-gamma of an American spread option by common-random bumps.
>
> Reprices :func:`bermudan_spread_lsm` at bumped spots on the *same* seed, so
> the two simulations share their Brownian shocks and the finite differences
> are low-variance. Returns a dict with ``price``, the two spot deltas
> (``delta1`` = dV/dS1, ``delta2`` = dV/dS2), the two own-gammas
> (``gamma1``, ``gamma2``), and the cross-gamma (``cross`` = d2V/dS1 dS2).
>
> For a spread *call* ``max(S1 - S2 - K, 0)`` the long-leg delta is positive
> and the short-leg delta negative; the LSM regression is re-fit at each bump.
> Deltas are reliable; the gammas (second differences over a re-fit regression)
> are indicative and need many paths.

## lsv

### `calibrate_lsv_leverage(S0, r, local_vol_fn, kappa, theta, xi, rho, v0, expiries, k_grid, q=0.0, n_paths=20000, seed=None, sub_steps=20)`  _function_

> Calibrate the LSV leverage surface by the particle method.
>
> Args:
>     local_vol_fn: target Dupire local vol ``sigma(K, t)``.
>     kappa, theta, xi, rho, v0: the backbone Heston variance parameters.
>     expiries: increasing calibration times (the leverage is piecewise
>         constant in time between them; ``t=0`` implied at the front).
>     k_grid: log-moneyness bin edges (relative to the forward) for the
>         conditional-expectation binning.
>
> Returns ``leverage`` -- a dict ``{t: {k_center: L}}`` -- and the callable
> ``lev_fn(spot, t)`` that interpolates it, suitable for an LSV Monte Carlo.

## mbs

### `amortization_schedule(balance, annual_rate, term_months)`  _function_

> Level-payment amortization schedule with no prepayment.
>
> Returns ``[(month, interest, principal, ending_balance), ...]``. The balance
> amortizes to (floating-point) zero at the final month.

### `cpr_to_smm(cpr)`  _function_

> Single monthly mortality from an annual CPR: ``1 - (1 - CPR)^{1/12}``.

### `mbs_cashflows(balance, annual_rate, term_months, smm=0.0)`  _function_

> Projected MBS cashflows with a constant SMM prepayment.
>
> Each month pays scheduled interest and principal on the surviving balance,
> plus a prepayment of ``smm`` times the balance remaining after the scheduled
> principal. Returns ``[(month, interest, scheduled_principal, prepayment,
> total_principal, ending_balance), ...]``. With ``smm = 0`` the total principal
> matches :func:`amortization_schedule`.

### `mbs_cashflows_psa(balance, annual_rate, term_months, psa=100.0)`  _function_

> Projected MBS cashflows on the PSA prepayment ramp (age-varying SMM).
>
> Like :func:`mbs_cashflows` but the monthly prepayment uses the age-dependent
> :func:`psa_cpr` converted to SMM at each month, rather than a constant SMM.
> Returns the same ``[(month, interest, scheduled_principal, prepayment,
> total_principal, ending_balance), ...]`` rows. At ``psa = 0`` it reduces to the
> no-prepayment schedule.

### `mbs_effective_convexity(cashflows, annual_yield, bump=0.0001)`  _function_

> Effective convexity of an MBS from a parallel yield bump.
>
> ``(P(y+h) - 2 P(y) + P(y-h)) / (h^2 P(y))`` on :func:`mbs_price` (static
> cashflows).

### `mbs_effective_duration(cashflows, annual_yield, bump=0.0001)`  _function_

> Effective duration of an MBS from a parallel yield bump (central difference).
>
> ``-(P(y+h) - P(y-h)) / (2 h P(y))`` on :func:`mbs_price`. Assumes the cashflows
> are held fixed (a static-duration measure; true option-adjusted duration would
> re-project prepayment at each bumped yield).

### `mbs_price(cashflows, annual_yield)`  _function_

> Present value of projected MBS cashflows at a monthly-compounded yield.
>
> ``cashflows`` are :func:`mbs_cashflows` rows; each month's cash is
> ``interest + total_principal`` discounted by ``(1 + y/12)^{-month}``. Monotone
> decreasing in ``annual_yield``.

### `mbs_price_with_spread(cashflows, zero_rates, spread)`  _function_

> Present value discounting each cashflow at its zero rate plus a spread.
>
> ``zero_rates[i]`` is the monthly-compounded annualized zero rate for the cash
> at ``cashflows[i]``'s month; every flow is discounted at ``zero_rate + spread``
> (a parallel add-on, the static/Z-spread convention). Reduces to
> :func:`mbs_price` at a flat curve.

### `mbs_yield(cashflows, price, tol=1e-10, max_iter=100)`  _function_

> Monthly-compounded annual yield reproducing an MBS ``price``.
>
> Bisection on :func:`mbs_price` (monotone decreasing in yield). Inverse of
> :func:`mbs_price`.

### `mbs_zspread(cashflows, zero_rates, price, tol=1e-12, max_iter=100)`  _function_

> Static (Z-) spread over the zero curve reproducing an MBS ``price``.
>
> Bisection on the constant spread added to every zero rate (price is monotone
> decreasing in the spread). Inverse of :func:`mbs_price_with_spread`.

### `monthly_payment(balance, annual_rate, term_months)`  _function_

> Level fully-amortizing monthly payment for a fixed-rate mortgage.
>
> ``P = B * i / (1 - (1 + i)^{-n})`` with monthly rate ``i = annual_rate / 12``.
> At zero rate this is the straight-line ``balance / term_months``.

### `pac_schedule(balance, annual_rate, term_months, psa_low, psa_high)`  _function_

> Planned-amortization-class principal schedule from a PSA collar.
>
> A PAC bond promises the principal that is available under *both* ends of a PSA
> speed band: at each month the scheduled PAC principal is the minimum of the
> total principal produced at ``psa_low`` and at ``psa_high``
> (:func:`mbs_cashflows_psa`). Returns ``[(month, pac_principal), ...]``. Because
> it is a lower envelope, the PAC schedule is stable for any prepayment speed
> inside the collar -- the support (companion) tranche absorbs the difference.

### `pac_support_split(cashflows, pac_sched)`  _function_

> Allocate pool principal between a PAC band and its support tranche.
>
> At each month the PAC receives its scheduled principal (from
> :func:`pac_schedule`), capped by what the pool actually produces and by the
> PAC's remaining balance; the support tranche receives the remainder. Any PAC
> shortfall in a slow month is made up from later principal before the support is
> paid. Returns ``(pac_rows, support_rows)`` as ``[(month, principal), ...]``.
> The two principal streams sum to the pool principal each month.

### `psa_cpr(month, psa=100.0)`  _function_

> CPR on the PSA ramp at a given loan age, for ``psa`` percent of the model.
>
> Standard 100 PSA: ``CPR = 0.06 * min(month, 30) / 30`` (0.2%/month ramp to 6%
> at month 30, flat after). Scaled by ``psa / 100`` for other speeds.

### `sequential_cmo(cashflows, tranche_sizes)`  _function_

> Split MBS principal across sequential (plain-vanilla) CMO tranches.
>
> Principal from ``cashflows`` (:func:`mbs_cashflows` rows) is paid to tranches
> strictly in order: tranche 0 receives all principal until retired, then
> tranche 1, and so on. ``tranche_sizes`` are the initial tranche balances (must
> sum to the pool's total principal). Returns a list, one per tranche, of
> ``[(month, principal, ending_balance), ...]`` rows. Each tranche's principal
> sums to its size; earlier tranches retire first (shorter WAL).

### `smm_to_cpr(smm)`  _function_

> Annual CPR from a single monthly mortality: ``1 - (1 - SMM)^{12}``.
>
> Inverse of :func:`cpr_to_smm`.

### `tranche_wal(tranche_rows, tranche_size)`  _function_

> Weighted-average life (years) of a single CMO tranche.
>
> ``sum_m (month/12) * principal_m / tranche_size`` over the tranche's principal
> rows from :func:`sequential_cmo`.

### `weighted_average_life(cashflows, balance)`  _function_

> Weighted-average life (years) from projected principal cashflows.
>
> ``WAL = sum_m (month/12) * total_principal_m / balance``. Uses the
> ``total_principal`` column (index 4) of :func:`mbs_cashflows`. Falls as
> prepayment speeds up (principal returns sooner).

## mc_greeks

### `asian_pathwise_vega(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None, n_steps=50, n_paths=100000, antithetic=True, seed=None) -> quantforge.montecarlo.MCResult`  _function_

> Pathwise vega of a fixed-strike arithmetic-average Asian call/put.
>
> Differentiates the payoff along each path with respect to ``sigma``. With
> ``S_i = S0 exp(sum (b - sig^2/2) dt + sig sqrt(dt) Z_j)``, the pathwise
> sensitivity of each monitored spot is ``dS_i/dsig = S_i * (W_i - sig t_i)``
> where ``W_i = sqrt(dt) sum_{j<=i} Z_j`` is the accumulated Brownian motion,
> so the average's derivative is ``dA/dsig = mean_i dS_i/dsig`` and the payoff
> derivative is ``disc * 1_{A>K} * dA/dsig`` (put: ``-1_{A<K}``). Lower
> variance than a bump for this Lipschitz payoff; the kink at ``A = K`` is a
> measure-zero set.

### `barrier_lr_delta(S, K, H, t, r, sigma, option_type=<OptionType.CALL: 'call'>, barrier='down-out', b=None, rebate=0.0, n_steps=100, n_paths=100000, antithetic=True, seed=None) -> quantforge.montecarlo.MCResult`  _function_

> Delta of a discretely-monitored single-barrier option (likelihood ratio).
>
> A knock-out/knock-in payoff is discontinuous in the spot (a path that just
> grazes the barrier pays nothing), so the pathwise method is ill-defined. The
> likelihood-ratio method sidesteps this: in the discrete GBM path the initial
> spot enters only through the mean of the *first* log-increment,
> ``ln S_1 = ln S0 + (b - sig^2/2) dt + sig sqrt(dt) Z_1``, so the score of the
> path density with respect to ``S0`` is ``Z_1 / (S0 sig sqrt(dt))`` and
>
>     delta = E[ discounted_payoff * Z_1 / (S0 sig sqrt(dt)) ].
>
> Monitoring is discrete (hard touch at the ``n_steps`` dates), matching
> :func:`barrier_mc` with ``brownian_bridge=False``; a common-random-number
> finite-difference of that price is the natural cross-check. ``barrier`` is
> ``down-out``/``down-in``/``up-out``/``up-in``; "down" watches ``S <= H``,
> "up" watches ``S >= H``. ``rebate`` is paid at expiry to killed knock-outs
> or never-activated knock-ins.

### `lr_digital_delta(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None, cash=1.0, n_paths=200000, antithetic=True, seed=None) -> quantforge.montecarlo.MCResult`  _function_

> Delta of a cash-or-nothing digital by the likelihood-ratio method.
>
> The digital payoff is discontinuous, so pathwise delta is ill-defined, but
> the LR estimator ``E[payoff * Z/(S0 sigma sqrt t)]`` is fine.

### `lr_digital_greeks(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None, cash=1.0, n_paths=400000, antithetic=True, seed=None)`  _function_

> Delta, vega, gamma of a cash-or-nothing digital by likelihood ratio.
>
> The digital payoff ``cash * 1_{S_T > K}`` (call) is discontinuous, so the
> pathwise method is undefined for *any* of its Greeks. The likelihood-ratio
> method differentiates the log-normal density instead of the payoff, so the
> same one-step Black-Scholes score weights that :func:`lr_greeks` uses for a
> vanilla apply unchanged to the digital:
>
>     delta:  Z / (S0 sigma sqrt(t))
>     vega:   (Z^2 - 1)/sigma - Z sqrt(t)
>     gamma:  (Z^2 - Z sigma sqrt(t) - 1) / (S0^2 sigma^2 t)
>
> Returns a dict with ``price``, ``delta``, ``vega``, ``gamma`` (Monte Carlo
> means) and their ``*_se`` standard errors. Cross-checks the analytic
> :func:`quantforge.digital_greeks` (delta, gamma) and a sigma-bump of
> :func:`quantforge.cash_or_nothing` (vega).

### `lr_greeks(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None, n_paths=100000, antithetic=True, seed=None)`  _function_

> European delta, gamma, vega by the likelihood-ratio method.
>
> Returns a dict with ``price``, ``delta``, ``gamma``, ``vega`` (each a Monte
> Carlo mean) plus their ``*_se`` standard errors. Works for the discontinuous
> digital payoff too (the LR weights do not touch the payoff).

### `mixed_gamma(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None, n_paths=200000, antithetic=True, seed=None) -> quantforge.montecarlo.MCResult`  _function_

> European gamma by the mixed pathwise-likelihood-ratio estimator.
>
> Gamma is ``d(delta)/dS0``. The pathwise delta payoff is
> ``D = e^{-r t} 1_{S_T > K} S_T / S0`` (a call), which depends on ``S0`` both
> through ``S_T`` (density -> LR weight ``Z/(S0 sigma sqrt t)``) and the
> explicit ``1/S0`` factor. Differentiating,
>
>     gamma = E[ D * ( Z/(S0 sigma sqrt t) - 1/S0 ) ].
>
> This "pathwise-then-LR" combination is well-defined even though the pure
> pathwise gamma is not (the delta payoff has an indicator), and it is lower
> variance than a double likelihood-ratio. Cross-checks the Black-Scholes
> gamma.

### `pathwise_delta(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None, n_paths=100000, antithetic=True, seed=None) -> quantforge.montecarlo.MCResult`  _function_

> European delta by the pathwise method (smooth-payoff estimator).
>
> Delta of a call is ``e^{-r t} 1_{S_T > K} S_T / S0`` (put: ``-1_{S_T < K}``).
> Lower variance than the likelihood-ratio delta for these Lipschitz payoffs.

### `smoothed_digital_delta(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None, cash=1.0, eps_rel=0.02, n_paths=200000, antithetic=True, seed=None) -> quantforge.montecarlo.MCResult`  _function_

> Delta of a cash-or-nothing digital by a smoothed (call-spread) payoff.
>
> The digital's indicator is discontinuous, so its pathwise delta is
> undefined. Replacing the indicator with a narrow call-spread ramp of relative
> width ``eps_rel`` -- ``clamp((S_T - (K - eps/2)) / eps, 0, 1)`` for a call --
> makes the payoff Lipschitz, so the pathwise delta
> ``disc * cash * (dramp/dS_T) * (S_T / S0)`` is well-defined. It trades bias
> for variance in ``eps_rel``: a wider spread lowers the variance but adds
> smoothing bias, a narrower one reduces the bias (converging to the true
> digital delta as ``eps_rel -> 0``) but the ``1/eps`` ramp raises the variance.
> A single-pass, model-agnostic alternative to the likelihood-ratio estimator.

## meixner

### `meixner_greeks(S, K, t, r, a, b, d, option_type=<OptionType.CALL: 'call'>, q=0.0, cm_alpha=1.5)`  _function_

> Greeks of a Meixner option by central finite differences.
>
> Central differences of :func:`meixner_price` for the spot Greeks ``delta``
> (dV/dS), ``gamma`` (d2V/dS2), and ``theta`` (calendar decay), plus the
> asymmetry (skew) sensitivity ``d_b`` (dV/db). The ``d_b`` bump is clipped to
> keep ``b`` in ``(-pi, pi)`` on both sides. Returns a dict with ``price``,
> ``delta``, ``gamma``, ``theta``, ``d_b``.

### `meixner_price(S, K, t, r, a, b, d, option_type=<OptionType.CALL: 'call'>, q=0.0, cm_alpha=1.5, upper=200.0) -> float`  _function_

> Price a European option under the Meixner model via Carr-Madan inversion.
>
> Args:
>     a: jump-size scale (> 0).
>     b: asymmetry in ``(-pi, pi)``; ``b < 0`` gives a downward skew.
>     d: activity / tail parameter (> 0).
>     q: continuous dividend yield.
>     cm_alpha: Carr-Madan damping. The transform needs the moment-generating
>         function at ``s = cm_alpha + 1`` to be finite, i.e.
>         ``a (cm_alpha + 1) + b < pi`` (the ``cosh`` argument must stay off its
>         pole).
>
> Puts use put-call parity.

### `meixner_smile(S, strikes, t, r, a, b, d, q=0.0, cm_alpha=1.5)`  _function_

> Black-Scholes implied-vol smile the Meixner model produces.
>
> Prices a call at each strike and inverts to a Black-Scholes implied vol,
> returning ``(log_moneyness, vol)`` pairs sorted by strike on the forward
> ``F = S e^{(r-q) t}``. ``b < 0`` tilts the smile into a downward skew.

## merton

### `merton_jump_greeks(S, K, t, r, sigma, lam, mu_j, sigma_j, option_type=<OptionType.CALL: 'call'>, b=None)`  _function_

> Greeks of a Merton jump-diffusion option by central finite differences.
>
> Differentiates :func:`merton_jump_price` for ``delta`` (dV/dS), ``gamma``
> (d2V/dS2), ``vega`` (dV/dsigma, the *diffusion*-vol sensitivity), and
> ``theta`` (calendar decay). At ``lam = 0`` (no jumps) the Greeks reduce to the
> vanilla Black-Scholes Greeks. Returns a dict with ``price`` and those fields.

### `merton_jump_price(S, K, t, r, sigma, lam, mu_j, sigma_j, option_type=<OptionType.CALL: 'call'>, b=None, max_terms=200, tol=1e-12) -> float`  _function_

> Price a European option under the Merton jump-diffusion model.
>
> Args:
>     sigma: diffusion volatility (the continuous part).
>     lam: jump intensity (expected number of jumps per year, >= 0).
>     mu_j: mean of the log jump size.
>     sigma_j: standard deviation of the log jump size (>= 0).
>     b: cost of carry (defaults to r). The drift is compensated so the
>         discounted asset is a martingale under the given carry.
>
> Returns the option price as the Poisson-weighted BSM series.

### `merton_smile(S, strikes, t, r, sigma, lam, mu_j, sigma_j, b=None)`  _function_

> The Black-Scholes implied-vol smile a Merton jump-diffusion produces.
>
> Prices a European call at each strike under the jump-diffusion, then inverts
> each price to its Black-Scholes implied volatility, returning
> ``(log_moneyness, vol)`` pairs sorted by strike (log-moneyness on the forward
> ``F = S e^{b t}``). Jumps fatten the tails, so the smile curves up in the
> wings; a negative mean jump ``mu_j`` tilts it into a downward skew.

## mlmc

### `mlmc_asian(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None, levels=4, M=2, n_paths=20000, seed=None)`  _function_

> Arithmetic-average Asian price by multi-level Monte Carlo.
>
> Runs levels ``0..levels`` with refinement factor ``M`` (fine grid at the top
> has ``M^levels`` steps). ``n_paths`` is the sample count at level 0; deeper
> levels use fewer (``n_paths // M^l``, floored) since their corrections have
> lower variance. Returns an :class:`~quantforge.MCResult` whose ``n_paths`` is
> the total sample count across levels.
>
> A flat run at ``levels=0`` is plain single-grid Monte Carlo; increasing
> ``levels`` refines the time discretisation while sharing the cost across
> coarser levels.

## moment_premium

### `moment_risk_premia(closes: Sequence[float], S0, t, r, vol_fn, q=0.0, n_strikes=401, width=8.0)`  _function_

> Realized vs risk-neutral (BKM) skewness and excess kurtosis.
>
> Args:
>     closes: realized price history over the measurement window (its log
>         returns give the realized moments).
>     S0, t, r, vol_fn, q: inputs for the BKM risk-neutral moments implied by
>         the option smile at horizon ``t``.
>
> Returns a dict with ``realized_skew``, ``implied_skew``, ``skew_premium``
> (implied - realized), and the analogous ``*_kurt`` excess-kurtosis fields.

## montecarlo

### `MCResult(price: float, std_error: float, n_paths: int) -> None`  _class_

> MCResult(price: float, std_error: float, n_paths: int)

### `arithmetic_asian_mc(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None, n_steps=50, n_paths=50000, antithetic=True, control_variate=True, seed=None) -> quantforge.montecarlo.MCResult`  _function_

> Price a fixed-strike arithmetic-average-price Asian option.
>
> With ``control_variate=True`` the *discretely*-monitored geometric-average
> Asian (known in closed form, over the same ``n_steps`` dates) is used as a
> control, dramatically reducing the standard error since the two averages are
> almost perfectly correlated.

### `autocallable_mc(S, t, r, sigma, observation_times, autocall_barrier, coupon, protection_barrier=None, notional=1.0, b=None, n_paths=50000, antithetic=True, seed=None) -> quantforge.montecarlo.MCResult`  _function_

> Monte Carlo an autocallable structured note.
>
> At each observation date, if the spot is at or above ``autocall_barrier``
> the note redeems early paying ``notional * (1 + coupon * k)`` where ``k`` is
> the observation number (accrued coupons), discounted to today. If it never
> autocalls, at maturity the holder gets the notional back unless the spot
> finished below ``protection_barrier`` (a down-and-in put on the notional),
> in which case they take the downside ``notional * S_T / S``.
>
> Args:
>     observation_times: increasing dates (years); the last is maturity.
>     autocall_barrier / protection_barrier: spot levels (same units as S).
>     coupon: coupon rate paid per elapsed observation on early redemption.

### `average_strike_asian_mc(S, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None, n_steps=50, n_paths=50000, antithetic=True, seed=None) -> quantforge.montecarlo.MCResult`  _function_

> Monte Carlo an average-strike Asian option.
>
> The strike is the realized arithmetic average of the monitored path, so a
> call pays ``max(S_T - A, 0)`` and a put ``max(A - S_T, 0)``, where ``A`` is
> the average over the ``n_steps`` monitoring dates. There is no simple closed
> form; the average and terminal spot come from the same simulated path.

### `barrier_digital_mc(S, K, H, t, r, sigma, option_type=<OptionType.CALL: 'call'>, barrier='up-in', b=None, cash=1.0, n_steps=100, n_paths=50000, antithetic=True, seed=None) -> quantforge.montecarlo.MCResult`  _function_

> Monte Carlo a cash-or-nothing digital contingent on a barrier condition.
>
> Pays ``cash`` at expiry if the option finishes in the money (call: S_T > K;
> put: S_T < K) AND the barrier condition holds over the monitored path:
>
>   * ``"up-in"``   / ``"down-in"``   : the barrier H must be touched;
>   * ``"up-out"``  / ``"down-out"``  : the barrier H must NOT be touched.
>
> "up" barriers watch for S >= H, "down" for S <= H. This is the standard
> barrier-contingent binary; the path dependence has no simple closed form.

### `barrier_mc(S, K, H, t, r, sigma, option_type=<OptionType.CALL: 'call'>, barrier='down-out', b=None, rebate=0.0, n_steps=100, n_paths=50000, antithetic=True, seed=None, brownian_bridge=True) -> quantforge.montecarlo.MCResult`  _function_

> Monte Carlo a single-barrier vanilla option with a Brownian-bridge check.
>
> Prices the continuously-monitored single-barrier option that
> :func:`quantforge.barrier_option` gives in closed form, so it is the natural
> cross-check for that formula (including a continuous dividend yield ``q`` fed
> in as ``b = r - q``). ``barrier`` is ``down-out``/``down-in``/``up-out``/
> ``up-in``; "down" watches for ``S <= H`` and "up" for ``S >= H``.
>
> Naive discrete monitoring misses barrier crossings that happen *between*
> time steps and so systematically over-prices knock-outs. With
> ``brownian_bridge=True`` (the default) each step contributes the exact
> conditional probability that the bridge between its two endpoints touched
> ``H``; a path survives a knock-out only if it dodges the barrier on every
> bridge. This removes the discretisation bias and converges to the
> continuous-monitoring closed form.
>
> ``rebate`` is paid at expiry to knock-outs that are killed, or to knock-ins
> that never activate, matching the closed form's convention.

### `basket_option_lhs_mc(spots, weights, K, t, r, sigmas, corr, q=None, option_type=<OptionType.CALL: 'call'>, n_paths=50000, seed=None) -> quantforge.montecarlo.MCResult`  _function_

> Two-asset basket option ``max(w1 S1 + w2 S2 - K, 0)`` by Latin hypercube MC.
>
> The exact analogue of :func:`spread_option_lhs_mc` for a basket (weighted
> sum) payoff. Each of the two driving normals is stratified into ``n_paths``
> equiprobable bins with one draw per bin, the two bin orders are independently
> permuted, and the stratified uniforms map through the inverse normal CDF;
> asset 2's shock is correlated by ``corr z1 + sqrt(1 - corr^2) z2``.
>
> Args mirror :func:`quantforge.basket_option`: ``spots`` ``(S1, S2)``,
> ``weights`` ``(w1, w2)``, ``sigmas`` ``(s1, s2)``, ``corr`` the correlation,
> and optional ``q`` ``(q1, q2)`` dividend yields.
>
> Note on the reported ``std_error``: as with any LHS estimator the samples are
> dependent, so the returned SE uses the plain i.i.d. formula and overstates
> the true error -- wrap the call in :func:`replicated_mc` for an honest SE.
> This routine is the *unbiased* Monte Carlo reference for the moment-matched
> (approximate) :func:`quantforge.basket_option`.

### `capped_cliquet_mc(S, t, r, sigma, reset_times, local_cap=None, local_floor=0.0, global_cap=None, global_floor=0.0, b=None, n_paths=50000, antithetic=True, seed=None) -> quantforge.montecarlo.MCResult`  _function_

> Monte Carlo a locally- and globally-capped cliquet (ratchet).
>
> The payoff sums the periodic returns of the underlying over consecutive
> reset windows, clipping each period return to ``[local_floor, local_cap]``,
> then clips the running sum to ``[global_floor, global_cap]``. The result is
> discounted at ``r``. This is the standard capped-cliquet structured note;
> the caps make it path-dependent with no closed form.
>
> Args:
>     reset_times: increasing schedule, e.g. [0.25, 0.5, 0.75, 1.0]; the first
>         period runs from now (t=0) to reset_times[0].
>     local_cap / local_floor: per-period return bounds (cap None = uncapped).
>     global_cap / global_floor: bounds on the summed payoff.

### `digital_is_mc(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None, cash=1.0, shift=None, n_paths=100000, seed=None) -> quantforge.montecarlo.MCResult`  _function_

> Cash-or-nothing digital price by importance sampling (deep-OTM friendly).
>
> A deep-OTM digital is even harder to simulate plainly than a vanilla: the
> payoff is a bounded 0/``cash`` indicator, so a far strike gives a tiny hit
> probability ``p`` and a relative standard error that blows up like
> ``sqrt((1-p)/p)``. Sampling the terminal normal from a shifted mean
> ``N(mu, 1)`` and reweighting by the likelihood ratio
> ``L(z) = exp(-mu z + mu^2/2)`` moves paths into the money while staying
> unbiased. The default shift places the mean draw exactly on the strike
> boundary, ``mu* = (ln(K/S0) - (b - sig^2/2) t) / (sig sqrt(t))``, so about
> half the shifted paths pay -- near variance-optimal for the indicator.
>
> Cross-checks the closed-form :func:`quantforge.cash_or_nothing`; for a
> deep-OTM strike the standard error is far below a plain indicator estimator
> at equal paths.

### `double_knockout_mc(S, K, t, r, sigma, lower, upper, option_type=<OptionType.CALL: 'call'>, b=None, rebate=0.0, n_steps=100, n_paths=50000, antithetic=True, seed=None) -> quantforge.montecarlo.MCResult`  _function_

> Monte Carlo a double-knockout barrier option (a corridor).
>
> The option pays the vanilla payoff only if the spot stays strictly inside
> ``(lower, upper)`` for the whole monitored path; if either barrier is
> breached it knocks out and pays the cash ``rebate`` (at expiry, discounted).
> Also known as a double-barrier knock-out or "corridor" option.

### `european_cv_mc(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None, n_paths=100000, antithetic=True, seed=None) -> quantforge.montecarlo.MCResult`  _function_

> European price with the underlying as a control variate (optimal beta).
>
> Combines both variance-reduction techniques: antithetic sampling *and* a
> control variate. The discounted terminal spot ``Y = e^{-r t} S_T`` has the
> known mean ``E[Y] = S0 e^{(b - r) t}`` (the forward, discounted), and it is
> correlated with the discounted call/put payoff ``X``, so the controlled
> estimator ``X - beta (Y - E[Y])`` has lower variance for the regression-
> optimal ``beta = Cov(X, Y) / Var(Y)``. Beta is estimated from the same
> sample; the resulting O(1/N) bias is negligible at these path counts and is
> swamped by the variance reduction. Antithetic pairs are averaged into a
> single sample first so both controls act on the same draws.
>
> Cross-checks the closed-form Black-Scholes value and reports a standard
> error strictly below the plain :func:`european_mc` at equal path count.

### `european_is_adaptive_mc(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None, n_pilot=20000, n_grid=41, n_paths=100000, seed=None) -> quantforge.montecarlo.MCResult`  _function_

> European price by importance sampling with a pilot-tuned optimal shift.
>
> :func:`european_is_mc` centres the sampling shift on the strike, which is
> near-optimal for a digital but not exactly optimal for a vanilla (whose
> payoff keeps growing past the strike, pulling the best shift further OTM).
> This routine tunes the shift ``mu`` from a short pilot instead of guessing.
>
> For an estimator ``payoff(z) * L(z)`` with ``L(z) = exp(-mu z + mu^2/2)``,
> a change of measure gives the second moment under the shifted law in terms
> of plain ``N(0, 1)`` draws:
>
>     M(mu) = E_mu[(payoff L)^2] = E_0[payoff(z)^2 exp(-mu z + mu^2/2)].
>
> So one pilot sample of ``payoff(z)^2`` and ``z`` under ``N(0, 1)`` scores
> *every* candidate ``mu`` on a grid at negligible cost; the variance-minimising
> ``mu`` is the one with the smallest ``M(mu)``. The main run then samples at
> that ``mu`` with the standard likelihood-ratio correction, so the estimate
> stays unbiased regardless of the tuning. The grid spans ``[0, 1.5 |mu0|]``
> (or the mirror for OTM puts) around the strike-centring shift ``mu0``.
>
> Cross-checks the closed-form Black-Scholes value; for a deep-OTM vanilla its
> standard error is at or below :func:`european_is_mc` at equal main-run paths.

### `european_is_mc(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None, shift=None, n_paths=100000, seed=None) -> quantforge.montecarlo.MCResult`  _function_

> European price by importance sampling, for deep out-of-the-money options.
>
> A plain simulation of a far-OTM option wastes almost every path: the payoff
> is zero unless the terminal spot crosses a distant strike, so the estimator
> is dominated by the rare paths that do. Importance sampling draws the
> terminal normal from a *shifted* mean ``N(mu, 1)`` instead of ``N(0, 1)`` to
> push mass into the money, then corrects the bias with the likelihood ratio
>
>     L(z) = exp(-mu z + mu^2 / 2),
>
> so ``E_shifted[payoff * L] = E[payoff]`` is unbiased. The default ``shift``
> centres the terminal log-spot on the strike -- ``mu* = (ln(K/S0) - (b -
> sig^2/2) t) / (sig sqrt(t))`` -- which is near variance-optimal for a digital
> and a large reduction for a deep-OTM vanilla. Pass an explicit ``shift`` to
> override. Antithetic sampling is not used (it would fight the deliberate
> asymmetry of the shift).
>
> Cross-checks the closed-form Black-Scholes value; for a deep-OTM strike its
> standard error is far below the plain :func:`european_mc` at equal paths.

### `european_mc(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None, n_paths=100000, antithetic=True, seed=None) -> quantforge.montecarlo.MCResult`  _function_

> Monte Carlo price of a European option (converges to the BSM value).

### `european_stratified_mc(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None, n_strata=100, n_per=10, seed=None) -> quantforge.montecarlo.MCResult`  _function_

> European price by stratified sampling of the terminal normal.
>
> The single normal that drives the terminal spot is split into ``n_strata``
> equiprobable strata ``[(i)/n, (i+1)/n]`` in probability space. Drawing
> ``n_per`` uniforms *within* each stratum and mapping them through the inverse
> normal CDF spreads the draws evenly across the distribution, removing the
> clustering that inflates plain Monte Carlo variance. Because the strata are
> equiprobable the estimator is the simple average of the per-stratum means,
> and its variance is ``(1/n_strata^2) sum_i s_i^2 / n_per`` from the
> within-stratum sample variances -- always at or below the plain estimator,
> and much lower for a smooth payoff.
>
> Total paths drawn is ``n_strata * n_per``; ``n_per >= 2`` is required so each
> stratum's variance is estimable. Cross-checks the closed-form Black-Scholes
> value and reports a standard error below :func:`european_mc` at equal paths.

### `local_vol_mc(S, K, t, r, local_vol_fn, option_type=<OptionType.CALL: 'call'>, q=0.0, n_steps=100, n_paths=50000, antithetic=True, seed=None) -> quantforge.montecarlo.MCResult`  _function_

> Monte Carlo a European option under a Dupire local-volatility surface.
>
> Args:
>     local_vol_fn: callable ``sigma_loc(S, t_now)`` giving the instantaneous
>         local volatility at spot ``S`` and elapsed time ``t_now``.
>     q: continuous dividend yield (drift is ``r - q``).
>
> Evolves ``dS = (r - q) S dt + sigma_loc(S, t) S dW`` with an Euler step in
> log-space. For a flat local vol this reproduces the Black-Scholes price; for
> a genuine Dupire surface the discretely-simulated price is consistent with
> that surface's vanilla prices.

### `parisian_barrier_mc(S, K, H, t, r, sigma, window, option_type=<OptionType.CALL: 'call'>, barrier='down-out', b=None, n_steps=252, n_paths=40000, antithetic=True, seed=None) -> quantforge.montecarlo.MCResult`  _function_

> Monte Carlo a Parisian barrier option.
>
> Unlike a standard barrier (triggered by a single touch), a Parisian barrier
> triggers only if the spot stays on the barrier's far side for a *consecutive*
> elapsed time of at least ``window`` years. This makes the option robust to
> brief spikes through the level.
>
> ``barrier`` is one of ``down-out``/``down-in``/``up-out``/``up-in``. "down"
> watches for S <= H, "up" for S >= H. Knock-out pays the vanilla payoff
> unless the barrier is activated; knock-in pays only if it is.

### `replicated_mc(estimator, n_batches=30, base_seed=0) -> quantforge.montecarlo.MCResult`  _function_

> Honest standard error for a variance-reduced (dependent-sample) estimator.
>
> Stratified sampling, Latin hypercube, and low-discrepancy (QMC) estimators
> draw *dependent* samples, so the plain i.i.d. ``std_error`` those routines
> report understates -- or, for a negatively-correlated design, overstates --
> the true error and cannot be used to size a run (see the note on
> :func:`spread_option_lhs_mc`). The fix is batched replication: run the whole
> estimator ``n_batches`` times with distinct seeds and treat the batch prices
> as the i.i.d. sample. The batch means *are* independent, so their spread is
> an unbiased estimate of the estimator's true standard error.
>
> ``estimator`` is any callable ``seed -> MCResult | float`` (e.g.
> ``lambda s: european_stratified_mc(..., seed=s)``). Seeds are
> ``base_seed, base_seed + 1, ...``. Returns an :class:`MCResult` whose
> ``price`` is the mean of the batch prices, ``std_error`` is the across-batch
> standard error ``s / sqrt(n_batches)``, and ``n_paths`` is ``n_batches`` (the
> number of independent replications, not the per-batch path count).

### `spread_option_lhs_mc(S1, S2, K, t, r, sigma1, sigma2, rho, q1=0.0, q2=0.0, option_type=<OptionType.CALL: 'call'>, n_paths=50000, seed=None) -> quantforge.montecarlo.MCResult`  _function_

> Two-asset spread option ``max(S1 - S2 - K, 0)`` by Latin hypercube MC.
>
> A plain two-asset simulation draws the two independent normals freely, so
> clumps and gaps in each margin add variance. Latin hypercube sampling
> stratifies *each* dimension into ``n_paths`` equiprobable bins and takes one
> draw per bin, then independently permutes the two dimensions' bin order so
> the pair is decorrelated before the target correlation is imposed. Mapping
> the stratified uniforms through the inverse normal CDF gives two marginally
> well-spread normals ``z1, z2``; the second asset's shock is correlated in the
> usual way ``rho z1 + sqrt(1 - rho^2) z2`` (Cholesky of the 2x2). The margins
> of a spread payoff are close to linear in each normal, so LHS removes most of
> the variance a plain draw leaves in.
>
> Note on the reported ``std_error``: the LHS samples are *not* independent, so
> the returned value is the plain i.i.d. formula and does **not** reflect the
> LHS variance reduction -- it overstates the true error. The genuine gain
> shows up only in the spread of the estimate across independent runs: at 4000
> paths the across-seed RMSE against Kirk is roughly 0.10 versus 0.27 for a
> plain two-asset draw (a ~2.6x reduction), even though both report a similar
> ``std_error``. Use replication, not the reported SE, to size an LHS run.
>
> Cross-checks the Kirk :func:`quantforge.spread_option` (and, at ``K = 0``,
> the exact Margrabe :func:`quantforge.exchange_option`).

## multiasset

### `basket_greeks(spots, weights, K, t, r, sigmas, corr, q=None, option_type=<OptionType.CALL: 'call'>)`  _function_

> Greeks of a two-asset basket option (Levy moment-match) by FD.
>
> Returns a dict with the two spot deltas, own-gammas, the cross-gamma, and
> the correlation sensitivity, all by central finite differences on
> :func:`basket_option`.

### `basket_option(spots, weights, K, t, r, sigmas, corr, q=None, option_type=<OptionType.CALL: 'call'>) -> float`  _function_

> Two-asset basket call/put on ``w1 S1 + w2 S2`` via lognormal moment match.
>
> Args:
>     spots: (S1, S2). weights: (w1, w2). sigmas: (sigma1, sigma2).
>     corr: correlation between the two assets.
>     q: optional (q1, q2) dividend yields; defaults to zeros.
>
> Matches the basket forward's first two moments to a single lognormal (Levy)
> and prices with Black-Scholes. Exact for a single asset; an approximation
> for the sum.

### `best_of_call(S1, S2, K, t, r, sigma1, sigma2, rho, q1=0.0, q2=0.0, option_type=<OptionType.CALL: 'call'>, n_paths=100000, antithetic=True, seed=None)`  _function_

> Option on the maximum of two assets: payoff max(max(S1,S2) - K, 0) (call).
>
> Monte Carlo on correlated GBM. Best-of and worst-of calls satisfy
> ``best + worst = call(S1) + call(S2)`` at the same strike (Stulz), which the
> tests check.

### `best_of_call_closed(S1, S2, K, t, r, sigma1, sigma2, rho, q1=0.0, q2=0.0)`  _function_

> Exact Stulz (1982) price of a call on the maximum of two assets.
>
> ``max(max(S1, S2) - K, 0)``. Uses the Stulz identity
> ``C_max + C_min = c(S1) + c(S2)`` (both vanilla calls at strike ``K``), so
> ``C_max = c(S1) + c(S2) - C_min`` with the exact :func:`_stulz_min_call`.
> This is the closed-form cross-check for the Monte Carlo :func:`best_of_call`.

### `best_of_put_closed(S1, S2, K, t, r, sigma1, sigma2, rho, q1=0.0, q2=0.0)`  _function_

> Exact price of a put on the maximum of two assets: ``max(K - max(S1,S2), 0)``.
>
>     P_max = C_max - disc E[max] + K e^{-r t},
>
> where ``disc E[max] = S1 e^{-q1 t} + S2 e^{-q2 t} - disc E[min]`` (the two
> forwards less the discounted expected min). Uses the exact
> :func:`best_of_call_closed`. Cross-check for the Monte Carlo
> :func:`best_of_call` put.

### `correlation_option(S1, S2, K1, K2, t, r, sigma1, sigma2, rho, option_type=<OptionType.CALL: 'call'>, cond2='above', q1=0.0, q2=0.0)`  _function_

> Two-asset correlation option: a vanilla on asset 1 gated by asset 2.
>
> Pays the asset-1 vanilla payoff at expiry only if asset 2 satisfies its
> barrier condition:
>
>     call: ``max(S1_T - K1, 0) * 1[cond2 on S2]``
>     put:  ``max(K1 - S1_T, 0) * 1[cond2 on S2]``
>
> with ``cond2`` = ``above`` (``S2_T > K2``) or ``below`` (``S2_T < K2``). It
> decomposes exactly into the two two-asset digitals already priced here:
>
>     call = AoN(S1>K1, cond2) - K1 * CoN(S1>K1, cond2)
>     put  = K1 * CoN(S1<K1, cond2) - AoN(S1<K1, cond2)
>
> where ``AoN`` is :func:`two_asset_asset_or_nothing` (pays ``S1_T``) and
> ``CoN`` is :func:`two_asset_digital` (pays 1). Both pieces are exact
> bivariate-normal closed forms, so the correlation option is too.
>
> Cross-checks a correlated-GBM Monte Carlo.

### `exchange_greeks(S1, S2, t, sigma1, sigma2, rho, q1=0.0, q2=0.0)`  _function_

> Greeks of a Margrabe exchange option (payoff max(S1 - S2, 0)) by FD.
>
> Returns a dict with the two spot deltas (``delta1`` = dV/dS1,
> ``delta2`` = dV/dS2), the two own-gammas (``gamma1``, ``gamma2``), the
> cross-gamma (``cross`` = d2V/dS1 dS2), and the correlation sensitivity
> (``corr_vega`` = dV/drho). All by central finite differences on the exact
> Margrabe formula.

### `exchange_option(S1, S2, t, sigma1, sigma2, rho, q1=0.0, q2=0.0) -> float`  _function_

> Margrabe option to exchange asset 2 for asset 1: payoff max(S1 - S2, 0).
>
> Exact closed form; independent of the risk-free rate (the two assets'
> financing cancels), depending only on the dividend yields.

### `geometric_basket_greeks(spots, weights, K, t, r, sigmas, corr, q=None, option_type=<OptionType.CALL: 'call'>)`  _function_

> Greeks of a geometric-basket option by central finite differences of
> :func:`geometric_basket_option`. Returns ``price`` plus per-asset ``delta``
> and ``gamma`` lists (dV/dS_i, d2V/dS_i^2) and the total ``vega`` (bumping all
> sigmas together) and ``theta`` (calendar decay).

### `geometric_basket_option(spots, weights, K, t, r, sigmas, corr, q=None, option_type=<OptionType.CALL: 'call'>) -> float`  _function_

> Weighted geometric-average basket option on ``prod_i S_i^{w_i}`` (exact).
>
> Unlike the arithmetic :func:`basket_option`, the weighted geometric average
> ``B = prod_i S_i^{w_i}`` of correlated lognormal assets is *itself*
> lognormal, so this has an exact closed form for any number of assets ``n``.
> With each asset ``S_i(t) = S_i exp((r - q_i - sigma_i^2/2) t + sigma_i W_i)``
> and correlations ``corr[i][j]``, ``log B`` is Gaussian with
>
>     mean = sum_i w_i (log S_i + (r - q_i - sigma_i^2/2) t)
>     var  = t sum_i sum_j w_i w_j corr[i][j] sigma_i sigma_j,
>
> and the price is a Black-Scholes-style formula on the basket forward
> ``F = exp(mean + var/2)`` discounted at ``r``.
>
> Args:
>     spots, weights, sigmas: length-``n`` sequences.
>     corr: ``n x n`` correlation matrix (list of lists).
>     q: optional length-``n`` dividend yields; defaults to zeros.

### `implied_exchange_correlation(target_price, S1, S2, t, sigma1, sigma2, q1=0.0, q2=0.0, tol=1e-10, max_iter=100)`  _function_

> Back out the correlation implied by a Margrabe exchange-option price.
>
> The exchange price depends on ``rho`` only through the spread vol
> ``sqrt(sigma1^2 - 2 rho sigma1 sigma2 + sigma2^2)``, which falls as ``rho``
> rises, so the price is monotone decreasing in ``rho`` -- a bisection on
> ``rho in (-1, 1)`` recovers it. Raises if the quote lies outside the range
> spanned by ``rho = -1 .. 1``.

### `implied_geometric_basket_correlation(target_price, spots, weights, K, t, r, sigmas, q=None, option_type=<OptionType.CALL: 'call'>, tol=1e-08, max_iter=100)`  _function_

> Back out the uniform pairwise correlation implied by a geometric-basket
> price.
>
> Assumes a single off-diagonal correlation ``rho`` shared by every pair
> (an equicorrelation matrix ``corr[i][j] = rho`` for ``i != j``, ``1`` on the
> diagonal). The basket log-variance
> ``t sum_ij w_i w_j corr[i][j] sigma_i sigma_j`` rises with ``rho``, so a
> geometric-basket call is monotone increasing in ``rho`` (a put decreasing),
> and a bisection recovers the correlation consistent with the quote.
>
> The search is bounded below by ``-1/(n-1)`` (the smallest ``rho`` keeping the
> equicorrelation matrix positive semidefinite) and above by ``1``. Raises if
> the quote lies outside the price range those bounds span.

### `implied_spread_correlation(target_price, S1, S2, K, t, r, sigma1, sigma2, q1=0.0, q2=0.0, option_type=<OptionType.CALL: 'call'>, tol=1e-08, max_iter=100)`  _function_

> Back out the correlation implied by a spread-option market price (Kirk).
>
> The Kirk spread price is monotone decreasing in ``rho`` (higher correlation
> lowers the spread volatility), so a bisection on ``rho in (-1, 1)`` recovers
> the correlation consistent with the quote. Raises if the quote lies outside
> the price range spanned by ``rho = -1 .. 1``.

### `implied_spread_correlation_bs(target_price, S1, S2, K, t, r, sigma1, sigma2, q1=0.0, q2=0.0, option_type=<OptionType.CALL: 'call'>, tol=1e-08, max_iter=100)`  _function_

> Correlation implied by a spread-option price under Bjerksund-Stensland 2014.
>
> Identical bisection to :func:`implied_spread_correlation` but inverts
> :func:`spread_option_bs` instead of the Kirk approximation. The BS spread
> price is likewise monotone decreasing in ``rho`` (higher correlation lowers
> the spread vol), so ``rho in (-1, 1)`` is recovered by bisection. Raises if
> the quote lies outside the price range spanned by ``rho = -1 .. 1``.

### `rainbow_greeks(S1, S2, K, t, r, sigma1, sigma2, rho, kind='best', option_type=<OptionType.CALL: 'call'>, q1=0.0, q2=0.0)`  _function_

> Greeks of a rainbow (best-of/worst-of) option by FD on the Stulz closed form.
>
> ``kind`` is ``"best"`` (option on the maximum) or ``"worst"`` (on the
> minimum). Differentiates the exact :func:`best_of_call_closed` /
> :func:`worst_of_call_closed` / :func:`best_of_put_closed` /
> :func:`worst_of_put_closed` -- no Monte Carlo noise -- for the two spot
> deltas, the two own-gammas, the cross-gamma ``d2V/dS1 dS2``, and the
> correlation sensitivity ``dV/drho``.
>
> A useful check: the best-of and worst-of *call* deltas in each asset sum to
> the corresponding single-asset Black-Scholes delta (differentiate the Stulz
> identity ``C_max + C_min = c(S1) + c(S2)``), and the max-call gains value as
> correlation falls (``corr_vega < 0``) while the min-call gains as it rises.

### `spread_greeks(S1, S2, K, t, r, sigma1, sigma2, rho, q1=0.0, q2=0.0, option_type=<OptionType.CALL: 'call'>)`  _function_

> Greeks of a Kirk spread option (payoff max(S1 - S2 - K, 0)) by FD.
>
> Returns a dict with the two spot deltas, own-gammas, the cross-gamma
> (``d2V/dS1 dS2``), and the correlation sensitivity (``corr_vega``). All by
> central finite differences on the Kirk approximation.

### `spread_option(S1, S2, K, t, r, sigma1, sigma2, rho, q1=0.0, q2=0.0, option_type=<OptionType.CALL: 'call'>) -> float`  _function_

> Kirk (1995) approximation for a spread option: payoff max(S1 - S2 - K, 0).
>
> Reduces to an exact Margrabe formula when K = 0. Puts follow from parity on
> the spread ``S1 - S2``.

### `spread_option_bs(S1, S2, K, t, r, sigma1, sigma2, rho, q1=0.0, q2=0.0, option_type=<OptionType.CALL: 'call'>) -> float`  _function_

> Bjerksund-Stensland (2014) spread-option approximation: max(S1 - S2 - K, 0).
>
> A three-``d`` closed form that is generally more accurate than Kirk
> (:func:`spread_option`) at wide strikes, high volatility, or when the two
> legs' vols differ sharply, while still reducing to the exact Margrabe value
> at ``K = 0``. Treats ``F2 + K`` as the effective second asset with weight
> ``b = F2 / (F2 + K)`` and prices
>
>     C = e^{-rt} [ F1 N(d1) - F2 N(d2) - K N(d3) ],
>
> where each ``d`` uses the blended spread vol
> ``sigma = sqrt(sigma1^2 - 2 b rho sigma1 sigma2 + b^2 sigma2^2)``. Puts follow
> from parity on the spread ``S1 - S2``.

### `spread_option_bs_greeks(S1, S2, K, t, r, sigma1, sigma2, rho, q1=0.0, q2=0.0, option_type=<OptionType.CALL: 'call'>)`  _function_

> Greeks of a Bjerksund-Stensland (2014) spread option by FD.
>
> Same layout as :func:`spread_greeks` (two spot deltas, own-gammas, cross-gamma
> ``d2V/dS1 dS2``, and ``corr_vega``) but differentiates
> :func:`spread_option_bs` instead of the Kirk approximation.

### `two_asset_asset_or_nothing(S1, S2, K1, K2, t, r, sigma1, sigma2, rho, cond1='above', cond2='above', q1=0.0, q2=0.0)`  _function_

> Asset-or-nothing digital paying ``S1_T`` iff both conditions hold.
>
> Pays the *first* asset's terminal value at expiry iff asset 1 is
> ``above``/``below`` ``K1`` and asset 2 is ``above``/``below`` ``K2``. Pricing
> under the asset-1 (share) measure -- where ``S1`` is the numeraire and asset
> 1's drift gains ``sigma1^2`` while asset 2's shock inherits an extra
> ``rho sigma1 sqrt(t)`` -- gives
>
>     S1 e^{-q1 t} * M(s1 a1, s2 a2; s1 s2 rho),
>
>     a1 = (ln(S1/K1) + (r - q1 + sigma1^2/2) t) / (sigma1 sqrt t)
>     a2 = (ln(S2/K2) + (r - q2 - sigma2^2/2) t) / (sigma2 sqrt t)
>          + rho sigma1 sqrt(t)
>
> with ``si = +1`` for an ``above`` condition, ``-1`` for ``below``. The four
> quadrant prices sum to the discounted forward ``S1 e^{-q1 t}`` (asset 1 is
> always delivered, on some quadrant). Cross-checks a correlated-GBM Monte
> Carlo. To pay asset 2 instead, swap the two assets in the call.

### `two_asset_digital(S1, S2, K1, K2, t, r, sigma1, sigma2, rho, cond1='above', cond2='above', q1=0.0, q2=0.0, cash=1.0)`  _function_

> Cash-or-nothing digital on two correlated assets (exact closed form).
>
> Pays ``cash`` at expiry iff both single-asset conditions hold: asset 1 is
> ``above`` (``S1_T > K1``) or ``below`` (``S1_T < K1``) its strike, and
> likewise for asset 2. Under the risk-neutral bivariate lognormal the price is
>
>     cash * e^{-r t} * M(s1 d1, s2 d2; s1 s2 rho)
>
> where ``di = (ln(Si/Ki) + (r - qi - sigma_i^2/2) t) / (sigma_i sqrt(t))`` is
> the usual ``d2``, ``si = +1`` for an ``above`` condition and ``-1`` for a
> ``below`` one, and ``M`` is the standardized bivariate-normal CDF. Flipping a
> condition flips the sign of that ``d`` and of the correlation. The four
> quadrant prices sum to ``cash e^{-r t}`` (the conditions are exhaustive).
>
> Cross-checks a correlated-GBM Monte Carlo.

### `two_asset_digital_greeks(S1, S2, K1, K2, t, r, sigma1, sigma2, rho, cond1='above', cond2='above', q1=0.0, q2=0.0, cash=1.0)`  _function_

> Greeks of a two-asset correlated digital by FD on the exact closed form.
>
> Differentiates :func:`two_asset_digital` -- no Monte Carlo noise -- for the
> two spot deltas (``delta1`` = dV/dS1, ``delta2`` = dV/dS2), the two
> own-gammas, the cross-gamma ``d2V/dS1 dS2``, and the correlation sensitivity
> ``corr_vega`` = dV/drho. Returns a dict with ``price`` and those fields.
>
> The correlation Greek is the interesting one: a both-``above`` (or
> both-``below``) digital *gains* value as correlation rises (the two
> in-the-money events move together), while a mixed above/below digital loses
> it; summed over the four exhaustive quadrants the correlation sensitivity is
> zero (total probability does not depend on ``rho``).

### `two_asset_gap_option(S1, S2, K_trigger, K_payoff, K2, t, r, sigma1, sigma2, rho, option_type=<OptionType.CALL: 'call'>, cond2='above', q1=0.0, q2=0.0)`  _function_

> Two-asset gap option: a gap payoff on asset 1 gated by asset 2.
>
> A gap option separates the *trigger* strike from the *payoff* strike. Here
> the asset-1 gap payoff fires only if asset 2 also clears its barrier:
>
>     call: ``(S1_T - K_payoff) * 1[S1_T > K_trigger] * 1[cond2 on S2]``
>     put:  ``(K_payoff - S1_T) * 1[S1_T < K_trigger] * 1[cond2 on S2]``
>
> (the payoff can be negative when ``K_payoff`` is on the far side of
> ``K_trigger`` -- the defining feature of a gap option). It decomposes into
> the two two-asset digitals with the *trigger* strike setting the asset-1
> condition and the *payoff* strike scaling the cash leg:
>
>     call = AoN(S1>K_trigger, cond2) - K_payoff * CoN(S1>K_trigger, cond2)
>     put  = K_payoff * CoN(S1<K_trigger, cond2) - AoN(S1<K_trigger, cond2)
>
> so it is an exact bivariate-normal closed form. With
> ``K_payoff = K_trigger`` it reduces to :func:`correlation_option`.
> Cross-checks a correlated-GBM Monte Carlo.

### `worst_of_call(S1, S2, K, t, r, sigma1, sigma2, rho, q1=0.0, q2=0.0, option_type=<OptionType.CALL: 'call'>, n_paths=100000, antithetic=True, seed=None)`  _function_

> Option on the minimum of two assets: payoff max(min(S1,S2) - K, 0) (call).

### `worst_of_call_closed(S1, S2, K, t, r, sigma1, sigma2, rho, q1=0.0, q2=0.0)`  _function_

> Exact Stulz (1982) price of a call on the minimum of two assets.
>
> ``max(min(S1, S2) - K, 0)``. Closed-form cross-check for the Monte Carlo
> :func:`worst_of_call`.

### `worst_of_put_closed(S1, S2, K, t, r, sigma1, sigma2, rho, q1=0.0, q2=0.0)`  _function_

> Exact price of a put on the minimum of two assets: ``max(K - min(S1,S2), 0)``.
>
> By put-call parity on the rainbow, a put and call on the same underlying
> (here ``min(S1, S2)``) satisfy ``C - P = disc E[min] - K e^{-r t}``, so
>
>     P_min = C_min - disc E[min] + K e^{-r t}
>
> with the exact :func:`worst_of_call_closed` and :func:`_disc_expected_min`.
> Closed-form cross-check for the Monte Carlo :func:`worst_of_call` put.

## nig

### `nig_greeks(S, K, t, r, alpha, beta, delta, option_type=<OptionType.CALL: 'call'>, q=0.0, cm_alpha=1.5)`  _function_

> Greeks of a NIG option by central finite differences.
>
> Central differences of :func:`nig_price` for the spot Greeks ``delta``
> (dV/dS), ``gamma`` (d2V/dS2), and ``theta`` (calendar decay), plus the
> process-parameter sensitivities ``d_alpha`` (dV/dalpha, tail steepness) and
> ``d_beta`` (dV/dbeta, skew). Returns a dict with ``price``, ``delta``,
> ``gamma``, ``theta``, ``d_alpha``, ``d_beta``.

### `nig_price(S, K, t, r, alpha, beta, delta, option_type=<OptionType.CALL: 'call'>, q=0.0, cm_alpha=1.5, upper=200.0) -> float`  _function_

> Price a European option under the NIG model via Carr-Madan inversion.
>
> Args:
>     alpha: tail-heaviness / steepness (> 0); larger = lighter tails.
>     beta: asymmetry (``|beta| < alpha``); ``beta < 0`` gives a downward skew.
>     delta: scale (> 0).
>     q: continuous dividend yield.
>     cm_alpha: Carr-Madan damping; needs ``alpha - (beta + cm_alpha + 1) > 0``
>         for the martingale transform to stay finite.
>
> Puts use put-call parity.

### `nig_smile(S, strikes, t, r, alpha, beta, delta, q=0.0, cm_alpha=1.5)`  _function_

> Black-Scholes implied-vol smile the NIG model produces.
>
> Prices a call at each strike and inverts to a Black-Scholes implied vol,
> returning ``(log_moneyness, vol)`` pairs sorted by strike on the forward
> ``F = S e^{(r-q) t}``. ``beta < 0`` tilts the smile into a downward skew.

## overhedge

### `Overhedge(cost: float, digital_value: float, cushion: float, long_strike: float, short_strike: float, quantity: float) -> None`  _class_

> Overhedge(cost: float, digital_value: float, cushion: float, long_strike: float, short_strike: float, quantity: float)

### `digital_call_overhedge(S, K, t, r, sigma, cash=1.0, width=None, b=None)`  _function_

> Super-replicate a cash-or-nothing CALL digital with a call spread.
>
> Longs ``cash/width`` calls at ``K - width`` and shorts the same at ``K``, so
> the payoff dominates ``cash * 1{S_T > K}``. Returns an :class:`Overhedge`
> with the spread cost (a conservative price), the fair digital value, and the
> cushion between them.

### `digital_put_overhedge(S, K, t, r, sigma, cash=1.0, width=None, b=None)`  _function_

> Super-replicate a cash-or-nothing PUT digital with a put spread.
>
> Longs ``cash/width`` puts at ``K + width`` and shorts the same at ``K``, so
> the payoff dominates ``cash * 1{S_T < K}``.

### `overhedge_payoff(oh: quantforge.overhedge.Overhedge, spot_at_expiry: float, is_call=True) -> float`  _function_

> Terminal payoff of the replicating spread at ``spot_at_expiry``.

## pde

### `crank_nicolson_barrier(S, K, H, t, r, sigma=None, option_type=<OptionType.CALL: 'call'>, barrier='down-out', b=None, rebate=0.0, local_vol_fn=None, n_space=400, n_time=400)`  _function_

> Price a continuously-monitored single-barrier option by a CN PDE.
>
> Knock-out barriers are imposed as an *absorbing* boundary: at every time
> step the value is set to ``rebate`` on the dead side of ``H`` (V = 0 there
> for a zero rebate), the exact continuous-monitoring condition. ``barrier`` is
> ``down-out``/``up-out`` for the directly-solved knock-outs, or
> ``down-in``/``up-in``, obtained from in + out = vanilla (same rebate handling
> as :func:`quantforge.barrier_option`).
>
> Supports a constant ``sigma`` or a ``local_vol_fn(S, t)`` and a carry
> ``b = r - q``. Returns the value at spot ``S``.

### `crank_nicolson_digital(S, K, t, r, sigma=None, option_type=<OptionType.CALL: 'call'>, b=None, cash=1.0, local_vol_fn=None, n_space=400, n_time=400, rannacher=2)`  _function_

> Price a cash-or-nothing digital by a Crank-Nicolson PDE.
>
> Terminal payoff is ``cash`` if the option finishes in the money (call:
> ``S_T > K``; put: ``S_T < K``), else 0. Boundary conditions are the digital's
> own: a call pays ``cash * e^{-r tau}`` at the top and 0 at 0 (mirror for a
> put). Supports a constant ``sigma`` or a ``local_vol_fn`` and a carry
> ``b = r - q``; the payoff kink at the strike makes Rannacher damping
> especially useful, so it is on by default.

### `crank_nicolson_greeks(S, K, t, r, sigma=None, option_type=<OptionType.CALL: 'call'>, b=None, american=False, local_vol_fn=None, n_space=200, n_time=200, s_max_mult=4.0, psor_tol=1e-08, psor_max_iter=10000, rannacher=2)`  _function_

> Price plus delta, gamma and theta read straight off the CN grid.
>
> Delta and gamma come from central finite differences of the final value
> grid in spot (no extra solves), and theta from the difference between the
> ``t=0`` grid and the grid one time step earlier. ``rannacher`` initial
> fully-implicit steps damp the payoff-kink oscillation that otherwise
> corrupts gamma near the strike. Returns a dict with price, delta, gamma and
> theta (calendar, per year).

### `crank_nicolson_no_touch(S, H, t, r, sigma=None, b=None, cash=1.0, local_vol_fn=None, n_space=400, n_time=400)`  _function_

> Price a no-touch binary (pays ``cash`` at expiry if ``H`` never hit).
>
> Solved as a knock-out of a constant ``cash`` payoff with an absorbing
> barrier at ``H`` (a node is placed exactly on ``H``). A one-touch that pays
> at expiry is ``cash * e^{-r t} - no_touch``.

### `crank_nicolson_price(S, K, t, r, sigma=None, option_type=<OptionType.CALL: 'call'>, b=None, american=False, local_vol_fn=None, n_space=200, n_time=200, s_max_mult=4.0, psor_tol=1e-08, psor_max_iter=10000, rannacher=2)`  _function_

> Price a European or American option by a Crank-Nicolson PDE solve.
>
> Provide either a constant ``sigma`` or a ``local_vol_fn(S, t)`` (time ``t``
> measured forward from today). ``b`` is the cost of carry (defaults to ``r``);
> dividend yield ``q`` enters as ``b = r - q``. American exercise uses PSOR.
> ``rannacher`` sets how many initial fully-implicit steps damp the payoff-kink
> oscillation (0 = pure Crank-Nicolson).
>
> Returns the option value interpolated at spot ``S``.

## pde2d

### `adi_spread_option(S1, S2, K, t, r, sigma1, sigma2, rho, q1=0.0, q2=0.0, n1=60, n2=60, n_time=40, width=5.0)`  _function_

> European spread-option price ``max(S1 - S2 - K, 0)`` by ADI.

### `adi_two_asset(payoff, S1, S2, t, r, sigma1, sigma2, rho, q1=0.0, q2=0.0, n1=60, n2=60, n_time=40, width=5.0)`  _function_

> Price a European two-asset option by Peaceman-Rachford ADI.
>
> Args:
>     payoff: callable ``payoff(s1, s2)`` giving the terminal value.
>     S1, S2: spot prices. sigma1, sigma2, rho: vols and correlation.
>     q1, q2: dividend yields. n1, n2, n_time: grid resolutions.
>     width: half-width of the log-price box in standard deviations.
>
> Returns the discounted option value interpolated at ``(S1, S2)``.

### `adi_two_asset_american(payoff, S1, S2, t, r, sigma1, sigma2, rho, q1=0.0, q2=0.0, n1=60, n2=60, n_time=40, width=5.0, american=True)`  _function_

> Two-asset option with optional early exercise, by Peaceman-Rachford ADI.
>
> Same discretisation as :func:`adi_two_asset` but, when ``american=True``,
> the value grid is floored at the immediate-exercise payoff after each time
> step (the explicit-payoff projection -- the 2D analogue of the vanilla PSOR
> floor). Suitable for American best-of / worst-of / spread payoffs.
>
> Returns the value interpolated at ``(S1, S2)``.

### `adi_two_asset_cs(payoff, S1, S2, t, r, sigma1, sigma2, rho, q1=0.0, q2=0.0, n1=60, n2=60, n_time=40, width=5.0, theta=0.5)`  _function_

> Price a European two-asset option by the Craig-Sneyd ADI scheme.
>
> The Peaceman-Rachford scheme (:func:`adi_two_asset`) is only first-order in
> time when a mixed (correlation) derivative is present. Craig-Sneyd fixes that
> with a Douglas predictor followed by a corrector that re-applies the explicit
> cross term at the predicted value, restoring second-order time accuracy. The
> directional operators ``A1``/``A2`` (each carrying half the ``r`` term) are
> solved implicitly (Thomas); the cross term ``A0`` stays explicit.
>
> Args and return value mirror :func:`adi_two_asset`.

## pde_asian

### `asian_pde_price(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None, n_s=120, n_i=120, n_time=100, s_max_mult=4.0)`  _function_

> Fixed-strike continuously-averaged arithmetic Asian by a 2D PDE.
>
> ``b`` is the cost of carry (defaults to ``r``); dividend yield ``q`` enters
> as ``b = r - q``. Averaging runs over the full life ``[0, t]``; the payoff is
> ``(A_t - K)^+`` for a call and ``(K - A_t)^+`` for a put, with
> ``A_t = (1/t) int_0^t S_u du``.
>
> Returns the value at spot ``S`` (the running integral starts at 0).

## perfmetrics

### `annualized_return(returns: Sequence[float], periods_per_year=252) -> float`  _function_

> Geometric (compound) annualized return.
>
> ``(prod(1 + r))^{periods_per_year / n} - 1`` -- the constant per-year rate
> that compounds to the realized total return over the sample.

### `annualized_volatility(returns: Sequence[float], periods_per_year=252) -> float`  _function_

> Annualized volatility: the sample standard deviation times
> ``sqrt(periods_per_year)``.

### `calmar_ratio(returns: Sequence[float], periods_per_year=252) -> float`  _function_

> Calmar ratio: annualized return divided by the maximum drawdown.
>
> Annualized return is the geometric ``(prod(1+r))^{periods_per_year/n} - 1``.
> Raises if there is no drawdown (undefined ratio).

### `cornish_fisher_var(returns, confidence=0.95, horizon=1.0) -> float`  _function_

> Cornish-Fisher (skew/kurtosis-adjusted) Value-at-Risk, as a positive loss.
>
> Expands the standard-normal quantile ``z`` at ``confidence`` with the sample
> skewness ``S`` and excess kurtosis ``K`` of the returns,
>
>     z_cf = z + (z^2-1) S/6 + (z^3-3z) K/24 - (2z^3-5z) S^2/36,
>
> evaluated at the lower-tail quantile ``z = Phi^{-1}(1-confidence)`` (a
> negative number), then ``VaR = -(mean*horizon + z_cf*sigma*sqrt(horizon))``
> as a positive loss. For a normal series it reduces to the parametric VaR;
> negative skew and fat tails fatten the left tail and push it above the
> Gaussian VaR.

### `cumulative_return(returns: Sequence[float]) -> float`  _function_

> Total compounded return over the series, ``prod(1 + r) - 1``.

### `down_capture(returns, benchmark_returns) -> float`  _function_

> Down-capture ratio: the asset's geometric return in down-benchmark
> periods over the benchmark's. Below 1 means the asset falls less than the
> benchmark in declining markets (good).

### `downside_beta(asset_returns, market_returns) -> float`  _function_

> Beta conditioned on down markets: ``Cov / Var`` over periods where the
> market return is negative.
>
> Measures how much the asset falls with the market on the downside. Computed
> on the subset of periods with ``market < 0`` using the ordinary sample
> covariance and variance. Raises if there are fewer than two down periods.

### `drawdown_curve(returns: Sequence[float]) -> list`  _function_

> Per-period underwater curve: fractional drop from the running peak.
>
> Compounds the returns into an equity curve and returns, for each period, the
> non-negative drawdown ``(peak - equity)/peak`` at that point (0 at a new
> high). The maximum of this curve is :func:`max_drawdown`.

### `historical_cvar(returns, confidence=0.95) -> float`  _function_

> Empirical conditional VaR (expected shortfall), as a positive loss.
>
> The average of the returns at or below the historical-VaR threshold, negated
> -- the mean loss in the worst ``1 - confidence`` of periods. Always at least
> the historical VaR. Falls back to the single worst return when the tail
> holds one observation.

### `historical_var_series(returns, confidence=0.95) -> float`  _function_

> Empirical (historical) Value-at-Risk, as a positive loss.
>
> The ``(1 - confidence)`` percentile of the return distribution, negated to a
> loss. Uses linear-interpolated order statistics -- no distributional
> assumption. A confidence of 0.95 reports the loss the returns exceed 5% of
> the time.

### `hit_rate(returns: Sequence[float]) -> float`  _function_

> Fraction of periods with a strictly positive return.

### `information_ratio(returns, benchmark_returns, periods_per_year=252) -> float`  _function_

> Information ratio: annualized active return over the tracking error.
>
> ``mean(active) * periods_per_year / tracking_error`` where the tracking
> error is itself annualized, so this equals
> ``mean(active) / stdev(active) * sqrt(periods_per_year)`` -- the Sharpe of
> the active-return series. Raises if the active returns have no variance.

### `jarque_bera(returns) -> float`  _function_

> Jarque-Bera test statistic for normality of a return series.
>
> ``JB = n/6 * (skew^2 + excess_kurt^2/4)``, asymptotically chi-squared with 2
> degrees of freedom under normality. Larger values reject normality (the 5%
> critical value is ~5.99). Uses the population skew/kurtosis.

### `longest_drawdown_duration(returns: Sequence[float]) -> int`  _function_

> Longest run of consecutive underwater periods (below a prior peak).
>
> Counts the maximum number of periods between a peak and the point the equity
> curve first recovers to (or exceeds) it. A series that never falls below its
> running peak returns 0.

### `max_drawdown(returns: Sequence[float]) -> float`  _function_

> Maximum peak-to-trough drawdown of the cumulative-return curve.
>
> Compounds the periodic returns into an equity curve and returns the largest
> fractional drop from a running peak, as a non-negative number (0.2 = a 20%
> drawdown). Empty or all-rising series give 0.

### `omega_ratio(returns: Sequence[float], threshold=0.0) -> float`  _function_

> Omega ratio: probability-weighted gains over losses about a threshold.
>
> ``sum(max(r - threshold, 0)) / sum(max(threshold - r, 0))`` -- the ratio of
> upside to downside area relative to ``threshold``. Values above 1 mean more
> gain mass than loss mass. Returns ``inf`` when there is no downside; raises
> if there is neither upside nor downside.

### `profit_factor(returns: Sequence[float]) -> float`  _function_

> Gross profits divided by gross losses (absolute).
>
> Returns ``inf`` when there are no losing periods. Raises if there are no
> gains and no losses.

### `rolling_sharpe(returns: Sequence[float], window: int, risk_free=0.0, periods_per_year=252) -> list`  _function_

> Annualized Sharpe ratio over each trailing window of ``window`` periods.
>
> Returns one Sharpe per window position (``len(returns) - window + 1``
> values), each computed by :func:`sharpe_ratio` on that slice. A
> zero-variance window yields ``float('nan')`` rather than raising, so the
> series stays aligned.

### `sample_kurtosis(returns, excess=True) -> float`  _function_

> Sample kurtosis (fourth standardized moment, population convention).
>
> ``(1/n) sum (x - mean)^4 / sigma^4``; with ``excess=True`` subtracts 3 so a
> normal distribution reads 0 (fat tails positive). Raises on zero variance.

### `sample_skewness(returns) -> float`  _function_

> Sample skewness (third standardized moment, population convention).
>
> ``(1/n) sum (x - mean)^3 / sigma^3`` with the population standard deviation
> (ddof=0). Positive means a longer right tail. Raises on a degenerate
> (zero-variance) series.

### `sharpe_ratio(returns: Sequence[float], risk_free=0.0, periods_per_year=252) -> float`  _function_

> Annualized Sharpe ratio of a periodic return series.
>
> ``(mean_excess / stdev) * sqrt(periods_per_year)`` where ``risk_free`` is the
> per-period risk-free return. Sample standard deviation (ddof=1).

### `sortino_ratio(returns: Sequence[float], risk_free=0.0, target=0.0, periods_per_year=252) -> float`  _function_

> Annualized Sortino ratio: excess mean over downside deviation.
>
> Downside deviation uses only returns below ``target`` (root-mean-square of
> the shortfalls, divided by the full sample count -- the standard
> convention). Raises if there is no downside.

### `tail_ratio(returns: Sequence[float], pct=5.0) -> float`  _function_

> Tail ratio: the right tail's magnitude over the left tail's.
>
> ``|percentile(100 - pct)| / |percentile(pct)|`` -- by default the 95th over
> the 5th percentile (in absolute value). Above 1 means the upside tail is
> fatter than the downside. Raises if the lower tail percentile is zero.

### `tracking_error(returns, benchmark_returns, periods_per_year=252) -> float`  _function_

> Annualized tracking error: stdev of the active (excess) return series.
>
> ``active_t = r_t - b_t``; the sample standard deviation (ddof=1) scaled by
> ``sqrt(periods_per_year)``. Series must be equal length.

### `up_capture(returns, benchmark_returns) -> float`  _function_

> Up-capture ratio: the asset's geometric return in up-benchmark periods
> divided by the benchmark's. Above 1 means the asset outpaces the benchmark
> in rising markets.

## perpetual

### `perpetual_american(S, K, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None)`  _function_

> Price a perpetual American option (Merton 1973), exact closed form.
>
> Args:
>     b: cost of carry (defaults to r). A perpetual call requires ``b < r``
>         (some carry cost / dividend) to be finite and worth exercising;
>         with ``b >= r`` the call is never exercised early and its value
>         tends to the spot.

### `perpetual_exercise_boundary(K, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None)`  _function_

> The optimal-exercise spot ``S*`` for a perpetual American option.

## portfolio

### `Book(positions: list = <factory>, net: quantforge.portfolio.BookRisk = <factory>) -> None`  _class_

> Book(positions: list = <factory>, net: quantforge.portfolio.BookRisk = <factory>)

### `BookRisk(market_value: float = 0.0, delta: float = 0.0, gamma: float = 0.0, vega: float = 0.0, theta: float = 0.0, rho: float = 0.0) -> None`  _class_

> Aggregate book-level exposures (position-scaled sums).

### `Contract(S: float, K: float, t: float, r: float, sigma: float, option_type: quantforge.bsm.OptionType = <OptionType.CALL: 'call'>, b: float = None, qty: float = 1.0, multiplier: float = 1.0, label: str = '') -> None`  _class_

> A single option position.
>
> ``qty`` is signed: positive = long, negative = short. ``multiplier`` scales
> each contract to its notional (e.g. 100 for US equity options).

### `Position(contract: quantforge.portfolio.Contract, greeks: quantforge.bsm.Greeks) -> None`  _class_

> A contract paired with its computed Greeks and position-scaled values.

### `price_book(contracts: Iterable[quantforge.portfolio.Contract]) -> quantforge.portfolio.Book`  _function_

> Value every contract and aggregate net book Greeks.
>
> Returns a ``Book`` with per-position detail and a ``net`` ``BookRisk`` of
> position-scaled (qty * multiplier) sums.

## portopt

### `black_litterman_returns(cov, market_weights, P, Q, tau=0.05, risk_aversion=2.5, omega=None) -> list`  _function_

> Black-Litterman posterior expected returns blending prior and views.
>
> The market-equilibrium prior ``Pi = lambda C w`` is combined with ``k``
> linear views ``P mu = Q`` (each row of ``P`` a portfolio, ``Q`` its expected
> return) of uncertainty ``omega`` (defaults to ``diag(tau P C P^T)``). The
> posterior mean is the standard closed form
>
>     mu = [ (tau C)^{-1} + P^T Omega^{-1} P ]^{-1}
>          [ (tau C)^{-1} Pi + P^T Omega^{-1} Q ].
>
> With no views (empty ``P``) it returns the prior ``Pi``.

### `black_litterman_weights(cov, market_weights, P, Q, tau=0.05, risk_aversion=2.5, omega=None, normalize=True)`  _function_

> Optimal portfolio weights from Black-Litterman posterior returns.
>
> Computes the posterior expected returns (:func:`black_litterman_returns`),
> then the unconstrained mean-variance optimum ``w = (lambda C)^{-1} mu``. With
> ``normalize=True`` the weights are rescaled to sum to 1 (fully invested);
> otherwise the raw utility-optimal holdings are returned. With no views the
> normalized weights reproduce the market weights (the prior is self-
> consistent).

### `component_var(weights, cov) -> list`  _function_

> Component (risk-contribution) VaR: each asset's share of portfolio vol.
>
> The marginal contribution ``(C w)_i / sigma_p`` times ``w_i`` gives the
> component ``w_i (C w)_i / sigma_p``; the components sum to the portfolio
> standard deviation. Scale by the VaR z-quantile to get VaR contributions.

### `diversification_ratio(weights, cov) -> float`  _function_

> Diversification ratio ``(sum_i w_i sigma_i) / sqrt(w^T C w)``.
>
> The weighted average of the assets' standalone volatilities over the
> portfolio volatility. Equals 1 for a single asset (or perfectly correlated
> assets) and rises as diversification lowers the portfolio vol below the
> weighted-average vol.

### `efficient_frontier(mean_returns, cov, targets) -> list`  _function_

> Efficient frontier as ``(target_return, portfolio_std)`` pairs.
>
> For each requested expected return in ``targets`` solves
> :func:`target_return_weights` and reports the achieved return with the
> portfolio standard deviation ``sqrt(w^T C w)``.

### `implied_equilibrium_returns(cov, market_weights, risk_aversion=2.5) -> list`  _function_

> Reverse-optimized (implied) equilibrium excess returns ``Pi = lambda C w``.
>
> Given the market-cap weights and a risk-aversion ``lambda``, the returns
> that make those weights mean-variance optimal are ``lambda C w`` -- the
> Black-Litterman market prior.

### `marginal_var(weights, cov, confidence=0.95, horizon=1.0) -> list`  _function_

> Marginal VaR: sensitivity of the portfolio VaR to each weight.
>
> ``dVaR/dw_i = z sqrt(horizon) (C w)_i / sigma_p`` (the zero-mean parametric
> VaR). Multiplying by ``w_i`` gives the component VaR, and the dot product
> ``sum_i w_i * marginal_i`` recovers the total VaR (VaR is homogeneous of
> degree 1 in the weights).

### `max_diversification_weights(cov) -> list`  _function_

> Most-diversified portfolio: maximizes the diversification ratio.
>
> The maximizer of ``(w^T sigma) / sqrt(w^T C w)`` is proportional to
> ``C^{-1} sigma`` (the tangency portfolio in the assets' own volatilities),
> normalized to sum to 1. Fully invested; may be long/short.

### `max_sharpe_weights(mean_returns, cov, risk_free=0.0) -> list`  _function_

> Tangency (max-Sharpe) weights ``C^{-1} (mu - rf) / sum(...)``.
>
> Maximizes the portfolio Sharpe ratio over fully-invested long/short
> weights. Requires the excess returns not to be orthogonal to ``C^{-1} 1``.

### `min_variance_weights(cov) -> list`  _function_

> Global minimum-variance weights ``C^{-1} 1 / (1^T C^{-1} 1)``.
>
> Fully invested (weights sum to 1); may be long/short. Requires a
> non-singular covariance matrix.

### `portfolio_cvar(weights, cov, mean_returns=None, confidence=0.95, horizon=1.0) -> float`  _function_

> Parametric (Gaussian) Conditional VaR / expected shortfall, as a loss.
>
> ``CVaR = phi(z)/(1-c) * sigma_p * sqrt(horizon) - mu_p * horizon``, the mean
> loss beyond the VaR under normality. Exceeds :func:`portfolio_var`.

### `portfolio_return(weights, mean_returns) -> float`  _function_

> Expected portfolio return ``w^T mu``.

### `portfolio_var(weights, cov, mean_returns=None, confidence=0.95, horizon=1.0) -> float`  _function_

> Parametric (Gaussian) Value-at-Risk of a portfolio, as a positive loss.
>
> ``VaR = z * sigma_p * sqrt(horizon) - mu_p * horizon`` where ``sigma_p`` is
> the portfolio standard deviation, ``mu_p`` the expected return (0 if
> ``mean_returns`` is omitted), and ``z`` the standard-normal quantile at
> ``confidence``. Returned as a non-negative loss figure.

### `portfolio_variance(weights, cov) -> float`  _function_

> Portfolio variance ``w^T C w``.

### `risk_parity_weights(cov, tol=1e-10, max_iter=1000) -> list`  _function_

> Equal-risk-contribution (risk-parity) weights.
>
> Solves for positive weights whose marginal risk contributions
> ``w_i (C w)_i`` are equal, by the standard fixed-point iteration
> ``w_i <- (target / (C w)_i)`` renormalized -- convergent for a positive-
> definite covariance. Weights are long-only and sum to 1.

### `target_return_weights(mean_returns, cov, target) -> list`  _function_

> Minimum-variance weights achieving an exact expected return ``target``.
>
> Solves ``min w^T C w`` subject to ``w^T 1 = 1`` and ``w^T mu = target`` by
> the two-constraint Lagrangian. With the efficient-frontier scalars
> ``A = 1^T C^{-1} 1``, ``B = 1^T C^{-1} mu``, ``C2 = mu^T C^{-1} mu`` and
> ``D = A C2 - B^2``, the weights are
>
>     w = C^{-1} [ (C2 - B target)/D * 1 + (A target - B)/D * mu ].
>
> Fully invested; may be long/short. Sweeping ``target`` traces the efficient
> frontier.

### `var_budget(weights, cov) -> list`  _function_

> Percentage VaR budget: each asset's fractional share of portfolio risk.
>
> ``w_i (C w)_i / (w^T C w)`` -- the component VaRs normalized to sum to 1.
> Independent of the confidence level and horizon (they cancel). Shows how the
> total risk is distributed across positions; equal entries mean risk parity.

## qmc

### `european_qmc(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None, n_points=8192) -> float`  _function_

> Price a European option by 1-D quasi-Monte Carlo (Halton) integration.
>
> Converges to the exact Black-Scholes value much faster than pseudo-random
> Monte Carlo for the same ``n_points``. Deterministic (no seed needed).

### `halton(index: int, dim: int) -> List[float]`  _function_

> The ``index``-th Halton point in ``dim`` dimensions (0-based index).
>
> Skips index 0 (the origin) by convention via a 1-based offset internally,
> so callers can pass 0, 1, 2, ... and get well-spread points.

## quanto

### `compo_option(S, K, t, r_domestic, r_foreign, sigma_asset, sigma_fx, rho, q_asset=0.0, option_type=<OptionType.CALL: 'call'>) -> float`  _function_

> Price a composite (compo) FX option: a foreign asset valued in domestic terms.
>
> Unlike a quanto (fixed FX), a compo option converts the foreign asset to
> domestic currency at the *floating* exchange rate, so the payoff is on the
> domestic-currency asset value ``X = S * FX``. Its volatility combines the
> asset and FX vols with their correlation:
>
>     sigma_compo = sqrt(sigma_asset^2 + sigma_fx^2 + 2 rho sigma_asset sigma_fx)
>
> Both ``S`` and ``K`` are quoted in domestic currency (K is the domestic
> strike on the converted asset). Carry and discounting use the domestic rate;
> the foreign rate enters as the asset's dividend-like yield ``q_asset``.

### `compo_option_greeks(S, K, t, r_domestic, r_foreign, sigma_asset, sigma_fx, rho, q_asset=0.0, option_type=<OptionType.CALL: 'call'>)`  _function_

> Greeks of a composite (compo) FX option.
>
> A compo option is a Black-Scholes price on the domestic-currency asset value
> with the combined volatility
> ``sigma_compo = sqrt(sigma_asset^2 + sigma_fx^2 + 2 rho sigma_asset sigma_fx)``
> and carry ``b = r_domestic - q_asset``. The spot enters only through the BSM
> price, so ``delta`` and ``gamma`` are exact BSM Greeks (no finite difference).
> ``vega`` (dV/dsigma_asset), ``fx_vega`` (dV/dsigma_fx), and ``corr_vega``
> (dV/drho) are central finite differences of the closed form; unlike a quanto,
> a compo is *long* FX volatility (positive ``fx_vega``). Returns a dict with
> ``price``, ``delta``, ``gamma``, ``vega``, ``fx_vega``, ``corr_vega``.

### `quanto_option(S, K, t, r_domestic, r_foreign, sigma_asset, sigma_fx, rho, q_asset=0.0, option_type=<OptionType.CALL: 'call'>) -> float`  _function_

> Price a quanto option (fixed-FX foreign-asset option in domestic terms).
>
> Args:
>     S, K: foreign-asset spot and strike (in foreign-asset units).
>     r_domestic: domestic risk-free rate (used for discounting).
>     r_foreign: foreign risk-free rate.
>     sigma_asset: volatility of the foreign asset.
>     sigma_fx: volatility of the domestic/foreign FX rate.
>     rho: correlation between the asset and the FX rate.
>     q_asset: dividend yield on the foreign asset.
>
> The price is in domestic currency per unit of the fixed exchange rate
> (multiply by the agreed FX level for the cash amount).

### `quanto_option_greeks(S, K, t, r_domestic, r_foreign, sigma_asset, sigma_fx, rho, q_asset=0.0, option_type=<OptionType.CALL: 'call'>)`  _function_

> Greeks of a quanto option.
>
> The quanto price is a Black-Scholes price on the foreign asset with the
> quanto-adjusted carry ``b_q = r_foreign - q_asset - rho sigma_asset
> sigma_fx``, discounted domestically. The spot enters only through that BSM
> price, so ``delta`` and ``gamma`` are the exact BSM Greeks at ``b_q`` (no
> finite difference). ``vega`` (dV/dsigma_asset -- which also moves ``b_q``),
> ``fx_vega`` (dV/dsigma_fx, the quanto's exposure to FX volatility), and
> ``corr_vega`` (dV/drho) are central finite differences of the closed form.
> Returns a dict with ``price``, ``delta``, ``gamma``, ``vega``, ``fx_vega``,
> ``corr_vega``.

## rates

### `CapletPeriod(forward: float, expiry: float, accrual: float, discount: float, sigma_n: float) -> None`  _class_

> CapletPeriod(forward: float, expiry: float, accrual: float, discount: float, sigma_n: float)

### `annuity(periods: Sequence[quantforge.rates.CapletPeriod]) -> float`  _function_

> Present-value annuity (level / PV01) of a swap: sum of accrual*discount.

### `black_swaption_greeks(swap_rate, strike, expiry, sigma_b, periods, payer=True)`  _function_

> Analytic Greeks of a Black (lognormal) European swaption.
>
> The value is ``annuity * Black76(swap_rate, strike, expiry, sigma_b)``, so
> the swap-rate Greeks are the Black-76 Greeks scaled by the annuity:
> ``rate_delta`` (dV/d swap_rate), ``rate_gamma`` (d2V/d swap_rate^2), and
> ``vega`` (dV/dsigma_b). Returns a dict with ``price``, ``rate_delta``,
> ``rate_gamma``, ``vega``, ``annuity``.

### `black_swaption_price(swap_rate, strike, expiry, sigma_b, periods, payer=True) -> float`  _function_

> Black (lognormal) price of a European swaption on the underlying swap.
>
> The market-standard lognormal counterpart to :func:`swaption_price`: the
> forward swap rate is modelled as lognormal with (Black) volatility
> ``sigma_b``, and the swaption is the annuity times a zero-carry Black-76
> option on the rate:
>
>     V = annuity * Black76(swap_rate, strike, expiry, sigma_b).
>
> A payer swaption is a call on the rate, a receiver a put. Requires positive
> ``swap_rate`` and ``strike`` (use :func:`swaption_price` for the normal model
> when rates may be negative).

### `cap_greeks(periods: Sequence[quantforge.rates.CapletPeriod], strike: float)`  _function_

> Aggregate Greeks of a cap: the summed caplet ``price``/``rate_delta``/
> ``rate_gamma``/``vega`` (all per unit notional).

### `cap_price(periods: Sequence[quantforge.rates.CapletPeriod], strike: float) -> float`  _function_

> Price an interest-rate cap as the sum of its caplets.

### `caplet_floorlet_parity(period: quantforge.rates.CapletPeriod, strike: float) -> float`  _function_

> Caplet - floorlet at the same strike = discounted forward-minus-strike.
>
> A put-call-parity identity used to check the pricer:
> ``caplet - floorlet = discount * accrual * (F - K)``.

### `caplet_greeks(period: quantforge.rates.CapletPeriod, strike: float, is_cap: bool = True)`  _function_

> Analytic Greeks of a single caplet/floorlet (normal model).
>
> The value is ``discount * accrual * Bachelier(F, K, expiry, 0, sigma_n)``, so
> its rate Greeks are the Bachelier Greeks in the forward rate scaled by the
> same ``discount * accrual`` factor: ``rate_delta`` (dV/dF), ``rate_gamma``
> (d2V/dF2), and ``vega`` (dV/dsigma_n). A caplet is a call on the forward, so
> its rate delta is positive; a floorlet's is negative. Returns a dict with
> ``price``, ``rate_delta``, ``rate_gamma``, ``vega``.

### `caplet_implied_normal_vol(price, period: quantforge.rates.CapletPeriod, strike, is_cap=True) -> float`  _function_

> Normal (Bachelier) implied vol of a caplet/floorlet from its price.
>
> Divides out the ``discount * accrual`` factor to recover the undiscounted
> Bachelier option value, then inverts it with :func:`bachelier_implied_vol`.
> Inverse of :func:`caplet_price` in the forward rate.

### `caplet_price(period: quantforge.rates.CapletPeriod, strike: float, is_cap: bool = True) -> float`  _function_

> Price a single caplet (cap) or floorlet (floor).
>
> Value = discount * accrual * Bachelier(F, K, expiry, r=0, sigma_n),
> with the option being a call for a caplet and a put for a floorlet. The
> Bachelier price is taken undiscounted (r=0) and discounted explicitly by the
> period's bond factor, which is the market convention.

### `collar_price(periods: Sequence[quantforge.rates.CapletPeriod], cap_strike: float, floor_strike: float) -> float`  _function_

> Price a collar: long a cap at ``cap_strike``, short a floor at ``floor_strike``.
>
> The net value is ``cap - floor``; a zero-cost collar is the pair of strikes
> that makes this zero.

### `compounded_overnight_rate(fixings, accruals) -> float`  _function_

> Annualized rate from daily-compounding overnight fixings (SOFR-style).
>
> The compounded setting-in-arrears rate over a period: multiply the daily
> growth factors ``(1 + r_i tau_i)`` and annualize by the total accrual,
>
>     rate = (prod_i (1 + r_i tau_i) - 1) / sum_i tau_i.
>
> ``fixings`` are the per-day annualized overnight rates and ``accruals`` the
> day-count fractions (typically 1/360). This is how compounded SOFR / SONIA
> coupons are computed.

### `compounded_rate_with_lockout(fixings, accruals, lockout=0) -> float`  _function_

> Compounded overnight rate with a rate lockout of ``k`` days.
>
> The final ``lockout`` business days of the period reuse the last observed
> fixing (the rate is *locked* before period end so the coupon is known early),
> the convention used for compounded fed funds. ``lockout = 0`` reduces to
> :func:`compounded_overnight_rate`.

### `compounded_rate_with_lookback(fixings, accruals, lookback=0) -> float`  _function_

> Compounded overnight rate with a lookback (observation-shift) of ``k`` days.
>
> Each accrual period uses the fixing observed ``lookback`` business days
> earlier (an *observation shift* also shifts the weighting to the earlier
> day). This gives the payment-lag convention used to publish a compounded SOFR
> coupon a few days before period end. ``lookback = 0`` reduces to
> :func:`compounded_overnight_rate`.
>
> ``fixings`` must extend at least ``lookback`` days before the accrual start;
> fixing ``i`` here is the rate applied to accrual ``i``, already shifted by the
> caller when constructing the arrays -- so this compounds ``fixings[i]`` over
> ``accruals[i]`` with an index offset used only to validate coverage.

### `floor_greeks(periods: Sequence[quantforge.rates.CapletPeriod], strike: float)`  _function_

> Aggregate Greeks of a floor: the summed floorlet Greeks.

### `floor_price(periods: Sequence[quantforge.rates.CapletPeriod], strike: float) -> float`  _function_

> Price an interest-rate floor as the sum of its floorlets.

### `sabr_cap_price(periods: Sequence[quantforge.rates.CapletPeriod], strike, alpha, beta, rho, nu) -> float`  _function_

> Price a cap under a single SABR smile (normal model).
>
> Sums caplets, each valued at the SABR normal vol read at its own forward and
> expiry -- so one calibrated ``(alpha, beta, rho, nu)`` prices the whole cap
> consistently across the smile, rather than a flat per-period ``sigma_n``.

### `sabr_floor_price(periods: Sequence[quantforge.rates.CapletPeriod], strike, alpha, beta, rho, nu) -> float`  _function_

> Price a floor under a single SABR smile (normal model): the sum of
> floorlets, each at the SABR normal vol of its own forward and expiry.

### `sabr_swaption_price(swap_rate, strike, expiry, periods, alpha, beta, rho, nu, payer=True, model='black')`  _function_

> Price a European swaption whose smile is described by a SABR model.
>
> Reads the SABR-implied volatility at the (forward swap rate, strike, expiry)
> point and feeds it into the matching swaption pricer:
>
>   * ``model="black"``  -> Hagan lognormal vol :func:`quantforge.sabr_vol`
>     into :func:`black_swaption_price` (requires positive rate and strike);
>   * ``model="normal"`` -> Hagan normal vol :func:`quantforge.sabr_normal_vol`
>     into :func:`swaption_price` (handles negative rates).
>
> This is the standard way SABR is used on swaptions: one calibrated smile
> prices every strike consistently. Returns the swaption present value.

### `simple_average_rate(fixings, accruals) -> float`  _function_

> Accrual-weighted arithmetic average of overnight fixings (Fed-funds style).
>
> ``sum_i r_i tau_i / sum_i tau_i`` -- the simple (non-compounded) averaging
> convention. Lies below the compounded rate when fixings are positive
> (compounding adds interest-on-interest).

### `swaption_greeks(swap_rate, strike, expiry, sigma_n, periods, payer=True)`  _function_

> Analytic Greeks of a European swaption (normal model).
>
> The value is ``annuity * Bachelier(swap_rate, strike, expiry, 0, sigma_n)``,
> so its swap-rate Greeks are the Bachelier Greeks scaled by the annuity:
> ``rate_delta`` (dV/d swap_rate), ``rate_gamma`` (d2V/d swap_rate^2), and
> ``vega`` (dV/dsigma_n). A payer swaption is a call on the swap rate (positive
> rate delta); a receiver is a put (negative). Returns a dict with ``price``,
> ``rate_delta``, ``rate_gamma``, ``vega``, ``annuity``.

### `swaption_implied_black_vol(price, swap_rate, strike, expiry, periods, payer=True) -> float`  _function_

> Black (lognormal) implied vol of a swaption from its price.
>
> Divides out the annuity to recover the undiscounted Black-76 option value,
> then inverts it with :func:`quantforge.implied_volatility` at zero carry.
> Inverse of :func:`black_swaption_price`; requires positive rate and strike.

### `swaption_implied_normal_vol(price, swap_rate, strike, expiry, periods, payer=True) -> float`  _function_

> Normal (Bachelier) implied vol of a swaption from its price.
>
> Divides out the annuity to recover the undiscounted Bachelier option value,
> then inverts it with :func:`bachelier_implied_vol`. Inverse of
> :func:`swaption_price`.

### `swaption_parity(swap_rate, strike, periods) -> float`  _function_

> Payer - receiver at the same strike = annuity * (swap_rate - strike).

### `swaption_price(swap_rate, strike, expiry, sigma_n, periods, payer=True) -> float`  _function_

> Bachelier price of a European swaption on the underlying swap.
>
> A payer swaption is a call on the swap rate; a receiver is a put. The value
> is the swap's PV annuity times a Bachelier option on the forward swap rate:
>
>     V = annuity * Bachelier(swap_rate, strike, expiry, r=0, sigma_n).
>
> Args:
>     swap_rate: current forward swap rate.
>     strike: fixed strike rate.
>     expiry: option expiry (years) — when the swap rate sets.
>     sigma_n: normal (absolute) volatility of the swap rate.
>     periods: the underlying swap's ``CapletPeriod`` legs, used only for the
>         annuity (accrual and discount factors).
>     payer: True for a payer (call), False for a receiver (put).
>
> Rates may be negative; the normal model handles that.

## rbergomi

### `rbergomi_greeks_cv(S, K, t, xi0, eta, H, rho, r=0.0, n_steps=100, n_paths=20000, antithetic=True, seed=None)`  _function_

> Greeks of a rough-Bergomi call by common-random-number bumps.
>
> Reprices the conditional (control-variate) estimator
> :func:`rbergomi_price_cv` at bumped inputs on the *same* seed, so the two
> simulations share their volatility-driving Brownian paths and the finite
> differences are low-variance. Returns a dict with ``price``, ``delta``
> (dV/dS), ``gamma`` (d2V/dS2), and ``vega_xi0`` (dV/dxi0 -- sensitivity to the
> forward-variance level, the rough-Bergomi analogue of vega). Only calls.

### `rbergomi_price(S, K, t, xi0, eta, H, rho, r=0.0, option_type=<OptionType.CALL: 'call'>, n_steps=100, n_paths=20000, antithetic=True, seed=None) -> quantforge.montecarlo.MCResult`  _function_

> Monte Carlo a European option under rough Bergomi.
>
> Args:
>     xi0: the (flat) forward variance curve level; ``sqrt(xi0)`` is the
>         baseline vol.
>     eta: vol-of-vol of the rough driver.
>     H: Hurst exponent in (0, 1); ``H < 0.5`` is rough (steep short skew),
>         ``H = 0.5`` recovers a standard lognormal-vol diffusion.
>     rho: spot/vol correlation (negative gives the equity down-skew).
>
> Returns an :class:`MCResult`.

### `rbergomi_price_cv(S, K, t, xi0, eta, H, rho, r=0.0, n_steps=100, n_paths=20000, antithetic=True, seed=None) -> quantforge.montecarlo.MCResult`  _function_

> Rough Bergomi European call by the conditional (turbocharged) estimator.
>
> Instead of simulating the orthogonal spot noise and averaging noisy payoffs
> (:func:`rbergomi_price`), this conditions on the volatility-driving Brownian
> motion and integrates the orthogonal noise out with a Black-Scholes formula
> (McCrickerd & Pakkanen, 2018). Every path contributes a smooth conditional
> price, so the Monte Carlo standard error drops sharply for the same paths --
> typically several-fold, and more as ``|rho|`` shrinks.
>
> Only calls are provided directly; puts follow from put-call parity on the
> forward ``S e^{r t}``.

### `rbergomi_smile(S, strikes, t, xi0, eta, H, rho, r=0.0, n_steps=100, n_paths=40000, antithetic=True, seed=None)`  _function_

> Black-Scholes implied-vol smile a rough Bergomi model produces.
>
> Simulates one set of terminal spots and reprices every strike on it (common
> random numbers), then inverts each call price to its Black-Scholes implied
> vol. Returns ``(log_moneyness, vol)`` pairs sorted by strike on the forward
> ``F = S e^{r t}``. Rough dynamics (``H < 0.5``) with ``rho < 0`` give the
> steep negative short-maturity skew that motivates the model.

### `rbergomi_smile_cv(S, strikes, t, xi0, eta, H, rho, r=0.0, n_steps=100, n_paths=40000, antithetic=True, seed=None)`  _function_

> Rough Bergomi implied-vol smile via the conditional estimator.
>
> Like :func:`rbergomi_smile` but prices each strike with the low-variance
> conditional call on a shared set of W paths, then inverts to a Black-Scholes
> vol. Returns ``(log_moneyness, vol)`` pairs sorted by strike.

## retirement

### `glide_path_equity_weight(years_to_target, glide_years, start_equity, end_equity)`  _function_

> Linear equity weight along a target-date glide path.
>
> Interpolates the equity allocation from ``start_equity`` (``glide_years`` out)
> down to ``end_equity`` at the target, clamped outside the glide window:
>
>     w = end + (start - end) * clamp(years_to_target / glide_years, 0, 1).
>
> Declines monotonically as the target approaches.

### `portfolio_depletion_years(balance, annual_withdrawal, real_return)`  _function_

> Years a portfolio lasts under a constant real withdrawal.
>
> Solves the annuity exhaustion ``balance = W * (1 - (1+g)^{-n}) / g`` for ``n``
> with real growth ``g = real_return``:
>
>     n = -ln(1 - g * balance / W) / ln(1 + g).
>
> Returns ``inf`` when the withdrawal is at or below the interest earned
> (``W <= g * balance``); at ``g = 0`` it is simply ``balance / W``.

### `ruin_probability_mc(balance, annual_withdrawal, mean_return, vol, years, inflation=0.0, n_paths=10000, seed=20260911)`  _function_

> Probability of portfolio ruin under lognormal returns (Monte Carlo).
>
> Simulates annual real returns as ``N(mean_return, vol^2)``, taking an
> inflation-indexed withdrawal at the start of each year, and reports the
> fraction of paths that hit zero before ``years``. Rises with the withdrawal
> rate and volatility; near zero for withdrawals well below the mean return.
> Deterministic per seed.

### `sustainable_withdrawal(balance, real_return, years)`  _function_

> Largest constant real withdrawal that exactly depletes over ``years``.
>
> The annuity payment ``W = balance * g / (1 - (1+g)^{-n})`` (``= balance / n`` at
> zero real return). The inverse of :func:`portfolio_depletion_years`.

### `withdrawal_balance_path(balance, annual_withdrawal, nominal_return, inflation, years)`  _function_

> Year-end balances under an inflation-indexed withdrawal.
>
> The withdrawal grows with ``inflation`` each year while the portfolio grows at
> the ``nominal_return``; withdrawals happen at year start. Returns the list of
> year-end balances (clipped at zero once depleted).

### `withdrawal_stream_pv(annual_withdrawal, real_discount_rate, years, growth=0.0)`  _function_

> Present value of a (possibly growing) real withdrawal stream.
>
> Discounts ``years`` annual withdrawals -- ``annual_withdrawal`` growing at
> ``growth`` per year -- at ``real_discount_rate``, withdrawals at year start:
>
>     PV = sum_{k=0}^{n-1} W (1+growth)^k / (1+r)^k.
>
> The capital needed to fund the stream. Rises with the withdrawal, the horizon,
> and the growth rate; falls with the discount rate.

## risk

### `VaRResult(var: float, expected_shortfall: float, confidence: float, horizon_days: float, method: str) -> None`  _class_

> VaRResult(var: float, expected_shortfall: float, confidence: float, horizon_days: float, method: str)

### `historical_var(book, return_scenarios, spot, confidence=0.99, horizon_days=1.0)`  _function_

> Historical VaR/ES: apply each realized return to a delta-gamma P&L.
>
> Args:
>     return_scenarios: iterable of simple returns r (e.g. daily), already at
>         the desired horizon.

### `montecarlo_var(contracts, sigma_annual, confidence=0.99, horizon_days=1.0, trading_days=252, n_paths=20000, seed=None)`  _function_

> Full-repricing VaR/ES: shock the spot, reprice every option, measure P&L.
>
> Unlike the parametric/historical estimators this makes no delta-gamma
> approximation — each position is repriced under the shocked spot with time
> advanced by the horizon.

### `parametric_var(book, sigma_annual, spot, confidence=0.99, horizon_days=1.0, trading_days=252)`  _function_

> Delta-gamma VaR/ES via a Cornish-Fisher expansion.
>
> Args:
>     book: a priced :class:`Book` (net.delta, net.gamma in spot terms).
>     sigma_annual: annualized volatility of the underlying's returns.
>     spot: current underlying price (to turn return shocks into price moves).
>     confidence: e.g. 0.99.
>     horizon_days / trading_days: scale vol to the VaR horizon.
>
> P&L ~= delta * dS + 0.5 * gamma * dS^2, with dS = spot * r_h and r_h normal
> with std ``sigma_h``. The gamma term makes P&L non-normal; Cornish-Fisher
> corrects the quantile using the P&L skewness.

## rnd

### `calendar_arbitrage_violations(expiries, vol_fns, ks=None, tol=1e-09)`  _function_

> Log-moneyness/expiry pairs where a smile term structure has calendar arbitrage.
>
> Calendar (horizontal-spread) arbitrage is absent when total implied variance
> ``w(k, t) = sigma(k, t)^2 t`` is non-decreasing in maturity at every fixed
> log-moneyness ``k = ln(K / F_t)``. This scans each adjacent expiry pair and
> reports the ``(k, t_lo, t_hi)`` points where ``w`` *decreases*
> (``w(k, t_hi) < w(k, t_lo) - tol``), which would let one buy the cheaper
> longer-dated variance and sell the richer shorter-dated one for a riskless
> profit.
>
> Args:
>     expiries: increasing list of expiries (years).
>     vol_fns: one smile ``vol_fn(k)`` per expiry, taking *log-moneyness* ``k``
>         and returning the Black implied vol. Order matches ``expiries``.
>     ks: log-moneyness grid to check (default ``[-1.5, 1.5]`` in 0.1 steps).
>
> Returns the list of violating ``(k, t_lo, t_hi)`` tuples; empty means the
> surface is calendar-arbitrage-free on the grid. Model-free: pass any smiles
> (SVI, SABR, vanna-volga, raw quotes) expressed in log-moneyness.

### `density_grid_from_smile(S0, t, r, vol_fn, q=0.0, n=400, width=8.0)`  _function_

> Return ``(strikes, density)`` of the risk-neutral density on a grid.
>
> Strikes span ``width`` standard deviations of log-moneyness around the
> forward. Useful for plotting or integrating custom payoffs.

### `price_payoff_from_density(S0, t, r, vol_fn, payoff, q=0.0, n=400, width=8.0)`  _function_

> Price a European payoff ``payoff(S_T)`` model-free from the smile density.
>
> ``price = e^{-r t} integral payoff(K) g(K) dK`` by the trapezoidal rule over
> the Breeden-Litzenberger density grid. A flat smile reprices vanillas and
> digitals to Black-Scholes.

### `risk_neutral_cdf_from_smile(S0, t, r, vol_fn, K, q=0.0, dK=None)`  _function_

> Risk-neutral CDF ``F(K) = P(S_T <= K)`` implied by an implied-vol smile.
>
> From Breeden-Litzenberger, the digital-put price is ``e^{-rt} P(S_T <= K)`` and
> equals ``-dC/dK`` discounted, so
>
>     F(K) = 1 + e^{r t} dC/dK,
>
> with the call priced at the smile vol ``vol_fn`` and ``dC/dK`` a central
> difference. Clamped to ``[0, 1]`` (a value hitting the clamp flags a smile
> that is not arbitrage-free at ``K``). A flat smile recovers the Black-Scholes
> ``N(-d2)``.

### `risk_neutral_cvar_from_smile(S0, t, r, vol_fn, alpha=0.99, q=0.0, n=8000, width=12.0)`  _function_

> Risk-neutral Conditional VaR (expected shortfall) of the terminal return.
>
> ``CVaR_alpha = E[L | L >= VaR_alpha]``, the average loss in the worst
> ``1 - alpha`` tail. With the loss ``L = 1 - S_T / S0`` and the tail threshold
> ``K = Q(1 - alpha)`` (so ``P(S_T <= K) = 1 - alpha``),
>
>     CVaR_alpha = 1 - E[S_T ; S_T <= K] / (S0 (1 - alpha)),
>
> where the truncated expectation ``E[S_T ; S_T <= K] = integral_0^K S g(S) dS``
> is taken against the Breeden-Litzenberger density ``g``. Returned as a
> positive fraction of ``S0`` and always ``>= VaR_alpha``. A flat smile matches
> the lognormal expected shortfall.

### `risk_neutral_density_from_smile(S0, t, r, vol_fn, K, q=0.0, dK=None)`  _function_

> Breeden-Litzenberger risk-neutral density ``g(K)`` at strike ``K``.
>
> ``g(K) = e^{r t} d^2 C / dK^2`` with the call priced at the smile vol; a
> negative value flags butterfly arbitrage in the smile there.

### `risk_neutral_quantile_from_smile(S0, t, r, vol_fn, p, q=0.0, dK=None, width=12.0, tol=1e-08, max_iter=200)`  _function_

> Inverse risk-neutral CDF: the strike ``K`` with ``P(S_T <= K) = p``.
>
> Bisection on :func:`risk_neutral_cdf_from_smile` over a log-moneyness bracket
> of ``+/- width`` forward standard deviations. ``p`` in ``(0, 1)``. Requires the
> smile CDF to be monotone on the bracket (true for an arbitrage-free smile).

### `risk_neutral_var_from_smile(S0, t, r, vol_fn, alpha=0.99, q=0.0, dK=None, width=12.0)`  _function_

> Risk-neutral Value-at-Risk of the terminal simple return over ``[0, t]``.
>
> Works with the loss ``L = 1 - S_T / S0`` (a positive number is a loss). The
> ``alpha``-VaR is the ``alpha``-quantile of ``L``: with probability ``alpha``
> the loss does not exceed it. Since ``L <= v`` iff ``S_T >= S0 (1 - v)``,
>
>     VaR_alpha = 1 - Q(1 - alpha) / S0,
>
> where ``Q`` is :func:`risk_neutral_quantile_from_smile`. Returned as a
> positive fraction of ``S0`` (e.g. ``0.18`` = an 18% loss). A flat smile
> matches the lognormal VaR.

### `smile_arbitrage_violations(S0, t, r, vol_fn, q=0.0, n=400, width=8.0, tol=1e-09)`  _function_

> Strikes where an implied-vol smile has butterfly (density) arbitrage.
>
> Scans a log-moneyness grid of ``width`` forward standard deviations and returns
> the strikes where the Breeden-Litzenberger density
> (:func:`risk_neutral_density_from_smile`) is negative beyond ``-tol``. A
> negative density means the call price is locally concave in strike -- a
> butterfly spread with a negative cost -- so the smile admits static arbitrage
> there. An empty list means the smile is butterfly-arbitrage-free on the grid.
> Model-free: works for any ``vol_fn`` (SVI, SABR, vanna-volga, raw quotes).

### `smile_is_arbitrage_free(S0, t, r, vol_fn, q=0.0, n=400, width=8.0, tol=1e-09) -> bool`  _function_

> True if the smile has no butterfly arbitrage on the scanned grid.
>
> Convenience wrapper: ``not smile_arbitrage_violations(...)``.

### `surface_arbitrage_report(S0, r, expiries, vol_fns, q=0.0, n=200, width=8.0, ks=None, tol=1e-09)`  _function_

> Full static-arbitrage report for a smile *surface* (butterfly + calendar).
>
> Combines the two model-free static-arbitrage tests:
>
>   * **butterfly** (per slice): the Breeden-Litzenberger density must stay
>     non-negative at every strike (see :func:`smile_arbitrage_violations`);
>   * **calendar** (across slices): total implied variance ``sigma^2 t`` must be
>     non-decreasing in maturity at fixed log-moneyness (see
>     :func:`calendar_arbitrage_violations`).
>
> ``vol_fns`` are smiles in *log-moneyness* ``k = ln(K / F_t)`` (one per expiry,
> matching ``expiries``), the same convention as
> :func:`calendar_arbitrage_violations`. Each is wrapped to a strike-based
> ``vol_fn(K)`` on that expiry's forward for the butterfly scan.
>
> Returns a dict ``{"butterfly": {t: [strikes]}, "calendar": [(k, t_lo, t_hi)]}``
> listing every violation; both empty means the surface is free of static
> arbitrage on the scanned grid.

### `surface_is_arbitrage_free(S0, r, expiries, vol_fns, q=0.0, n=200, width=8.0, ks=None, tol=1e-09) -> bool`  _function_

> True if a smile surface is free of both butterfly and calendar arbitrage.
>
> Convenience wrapper over :func:`surface_arbitrage_report`.

### `surface_is_calendar_arbitrage_free(expiries, vol_fns, ks=None, tol=1e-09) -> bool`  _function_

> True if the smile term structure has no calendar arbitrage on the grid.
>
> Convenience wrapper: ``not calendar_arbitrage_violations(...)``.

## rough_heston

### `rough_heston_greeks(S, K, t, r, v0, kappa, theta, nu, rho, H=0.1, option_type=<OptionType.CALL: 'call'>, q=0.0, n_grid=200)`  _function_

> Greeks of a rough-Heston option by central finite differences.
>
> Central differences of :func:`rough_heston_price` for the spot Greeks
> ``delta`` (dV/dS) and ``gamma`` (d2V/dS2), and the initial-variance
> sensitivity ``vega_v0`` (dV/dv0). Each re-price runs the O(n_grid^2)
> fractional-Riccati solve, so this is comparatively slow. Small ``H`` needs a
> fine grid to stay stable, so the default ``n_grid`` matches the pricer's. At
> ``H = 0.5`` the Greeks approach
> the classical Heston Greeks. Returns a dict with ``price``, ``delta``,
> ``gamma``, ``vega_v0``.

### `rough_heston_price(S, K, t, r, v0, kappa, theta, nu, rho, H=0.1, option_type=<OptionType.CALL: 'call'>, q=0.0, n_grid=200, upper=120.0) -> float`  _function_

> Price a European option under the rough-Heston model.
>
> Args:
>     v0, kappa, theta: initial variance, mean-reversion speed, long variance.
>     nu: volatility of variance *relative to kappa* (El Euch-Rosenbaum
>         convention). At ``H = 0.5`` this reduces to the classical Heston
>         model with vol-of-vol ``xi = kappa * nu``.
>     rho: spot/variance correlation (negative for the equity skew).
>     H: Hurst exponent in (0, 0.5]; ``H = 0.5`` recovers classical Heston.
>     n_grid: fractional-Riccati time-grid resolution (more = more accurate,
>         O(n_grid^2) work per Fourier node).
>
> Puts use put-call parity.

### `rough_heston_smile(S, strikes, t, r, v0, kappa, theta, nu, rho, H=0.1, q=0.0, n_grid=200)`  _function_

> Black-Scholes implied-vol smile the rough-Heston model produces.
>
> Returns ``(log_moneyness, vol)`` pairs sorted by strike on the forward
> ``F = S e^{(r-q) t}``. Small ``H`` steepens the short-dated skew beyond what
> classical Heston can reach.

## sabr

### `SABRParams(alpha: float, beta: float, rho: float, nu: float) -> None`  _class_

> SABRParams(alpha: float, beta: float, rho: float, nu: float)

### `calibrate_sabr(F, t, strikes: Sequence[float], market_vols: Sequence[float], beta: float = 0.5, weights: Sequence[float] = None, initial: quantforge.sabr.SABRParams = None, max_iter: int = 4000) -> Tuple[quantforge.sabr.SABRParams, float]`  _function_

> Fit (alpha, rho, nu) of a SABR smile to market Black vols; ``beta`` fixed.
>
> Returns ``(params, rmse)`` where rmse is the root-mean-square vol error.
> Uses a smooth constrained reparametrization so alpha > 0, nu >= 0 and
> rho in (-1, 1), optimized with the built-in Nelder-Mead.

### `calibrate_sabr_lm(F, t, strikes: Sequence[float], market_vols: Sequence[float], beta: float = 0.5, weights: Sequence[float] = None, initial: quantforge.sabr.SABRParams = None, max_iter: int = 100, tol: float = 1e-14)`  _function_

> Fit (alpha, rho, nu) of a SABR smile by Levenberg-Marquardt.
>
> Uses the exact analytic Jacobian (:func:`sabr_jacobian`) instead of the
> derivative-free Nelder-Mead of :func:`calibrate_sabr`, so it converges in a
> handful of iterations and lands on the same optimum. Parameters are box
> constrained to the valid region (``alpha > 0``, ``nu >= 0``,
> ``-1 < rho < 1``) by clamping each proposed step.
>
> Returns ``(params, rmse, n_iter)``.

### `sabr_bkm_moments(F, t, r, alpha, beta, rho, nu, q=0.0, n_strikes=401, width=8.0)`  _function_

> Risk-neutral (variance, skewness, excess kurtosis) implied by a SABR smile.
>
> Prices the strike chain at ``sabr_vol(F, K, ...)`` and applies the
> Bakshi-Kapadia-Madan moment replication
> (:func:`quantforge.bkm_moments_from_smile`). ``F`` is the forward;
> ``S0 = F e^{-(r-q)t}``. A negative correlation ``rho`` (equity skew) produces
> negative risk-neutral skewness; higher vol-of-vol raises the excess kurtosis.

### `sabr_butterfly_arbitrage(F, t, alpha, beta, rho, nu, strikes=None, r=0.0, tol=1e-08)`  _function_

> Return the strikes where the SABR smile has negative implied density.
>
> Scans ``strikes`` (default a wide grid around the forward) and reports every
> ``K`` whose Breeden-Litzenberger density is below ``-tol`` -- the butterfly-
> arbitrage points of the Hagan expansion (which is not guaranteed arb-free in
> the wings).

### `sabr_density(F, K, t, alpha, beta, rho, nu, r=0.0, dK=None)`  _function_

> Breeden-Litzenberger risk-neutral density of a SABR smile at strike ``K``.
>
> The implied density is ``g(K) = e^{r t} d^2 C / dK^2`` where ``C(K)`` is the
> Black call priced at the SABR vol ``sabr_vol(F, K)`` for each strike. Computed
> by a central second difference in strike. A *negative* density signals
> butterfly (static) arbitrage in the smile.

### `sabr_is_arbitrage_free(F, t, alpha, beta, rho, nu, strikes=None, r=0.0, tol=1e-08)`  _function_

> True if the SABR smile's implied density is non-negative on the grid.

### `sabr_jacobian(F, t, strikes: Sequence[float], alpha, beta, rho, nu)`  _function_

> Calibration Jacobian ``d sabr_vol(K_i) / d (alpha, rho, nu)``.
>
> Returns a list of ``[d_dalpha, d_drho, d_dnu]`` rows, one per strike, using
> the exact dual-number partials. This is the ``J`` a Gauss-Newton or
> Levenberg-Marquardt step needs, and ``(J^T J)^{-1}`` gives the asymptotic
> parameter covariance for standard errors on a fit.

### `sabr_normal_vol(F, K, t, alpha, beta, rho, nu) -> float`  _function_

> Hagan (2002) normal (Bachelier) implied volatility for the SABR model.
>
> Returns the absolute-vol ``sigma_N`` such that a Bachelier option on the
> forward reproduces the SABR price -- the standard quoting convention in
> interest-rate markets (where forwards can be near or below zero). Uses the
> Hagan normal expansion with the ATM limit (``F == K``) handled separately to
> avoid the removable ``z / x(z)`` singularity. The third-order time bracket
> matches the lognormal :func:`sabr_vol`; only the leading factor differs
> (``nu (F - K) / x(z)`` rather than ``alpha z / (denom x(z))``).

### `sabr_option_greeks(F, K, t, alpha, beta, rho, nu, option_type=<OptionType.CALL: 'call'>, discount=1.0)`  _function_

> Greeks of an option priced at the SABR smile volatility.
>
> The option is a Black-76 call/put on the forward ``F`` at the Hagan SABR
> implied vol ``sigma(F, K)``. Its **total** delta includes the smile backbone:
>
>     delta = dPrice/dF = black_delta + black_vega * dsigma/dF,
>
> where ``dsigma/dF`` is the exact (AD) backbone from
> :func:`sabr_sensitivities`. This differs from the pure Black delta because
> moving the forward also moves the SABR vol. Returns a dict with ``price``,
> ``vol``, ``delta`` (total, backbone-adjusted), ``black_delta`` (vol held
> fixed), ``vega`` (dPrice/dsigma), and ``gamma`` (total d2Price/dF2, including
> the backbone curvature, by a central difference of the SABR-repriced
> surface). Prices/greeks are on the *forward* (carry ``b = 0``); pass
> ``discount`` = P(0,T) to scale to present value.

### `sabr_repair_butterfly(F, t, params, r=0.0, strikes=None, max_iter=200, factor=0.95)`  _function_

> Repair a SABR smile's butterfly arbitrage by shrinking the vol-of-vol.
>
> The Hagan expansion loses density positivity when ``nu`` is large relative to
> ``t`` (steep, convex wings). This shrinks ``nu`` geometrically until
> :func:`sabr_is_arbitrage_free` passes, preserving alpha, beta and rho.
> Returns a new :class:`SABRParams` (the input if already arbitrage-free).

### `sabr_sensitivities(F, K, t, alpha, beta, rho, nu)`  _function_

> Exact partial derivatives of the Hagan SABR vol via forward-mode AD.
>
> Returns a dict with the vol itself and its machine-precision partials
>
>     ``vol`` and ``d_dF, d_dK, d_dalpha, d_drho, d_dnu``
>
> computed with dual numbers (no finite-difference truncation error). The
> (alpha, rho, nu) partials are the columns of the calibration Jacobian; the
> ``d_dF`` and ``d_dK`` partials give the smile's backbone and skew slopes.
>
> At exactly ``F == K`` the ATM branch of :func:`sabr_vol` is used, whose
> F/K partials describe that branch (a finite-difference bump moves off ATM);
> the alpha/rho/nu partials are exact everywhere.

### `sabr_variance_swap_strike(F, t, r, alpha, beta, rho, nu, q=0.0, n_strikes=401, width=8.0)`  _function_

> Fair variance-swap strike (annualized *variance*) implied by SABR.
>
> Replicates the variance swap from the SABR smile: each strike is priced at
> ``sabr_vol(F, K, ...)`` and fed to :func:`quantforge.variance_swap_from_smile`.
> ``F`` is the forward; the spot is ``S0 = F e^{-(r-q) t}``. Returns the fair
> variance; take ``sqrt`` for the fair volatility.

### `sabr_vix(F, t, r, alpha, beta, rho, nu, q=0.0, n_strikes=201, width=6.0)`  _function_

> VIX-style index (``~= 100 * sigma``) implied by a SABR smile.
>
> Prices each strike at ``sabr_vol(F, K, ...)`` and applies
> :func:`quantforge.vix_from_smile`. ``F`` is the forward; ``S0 = F e^{-(r-q)t}``.

### `sabr_vol(F, K, t, alpha, beta, rho, nu) -> float`  _function_

> Hagan (2002) lognormal (Black) implied volatility for the SABR model.
>
> Uses the standard expansion with the ATM limit handled separately to avoid
> the removable 0/0 singularity at ``F == K``.

## scenario

### `ScenarioGrid(spot_shocks: tuple, vol_shocks: tuple, pnl: tuple, base_value: float, relative: bool) -> None`  _class_

> ScenarioGrid(spot_shocks: tuple, vol_shocks: tuple, pnl: tuple, base_value: float, relative: bool)

### `spot_ladder(contracts: Sequence[quantforge.portfolio.Contract], spot_shocks, relative=True)`  _function_

> A 1-D price ladder: book P&L vs spot shock only (vol unchanged).
>
> Returns a list of (spot_shock, pnl) pairs.

### `stress_grid(contracts: Sequence[quantforge.portfolio.Contract], spot_shocks, vol_shocks, relative=True) -> quantforge.scenario.ScenarioGrid`  _function_

> Reprice a book across a Cartesian grid of spot and vol shocks.
>
> Args:
>     contracts: the positions (signed qty, multiplier as in ``price_book``).
>     spot_shocks: iterable of shocks to the underlying (e.g. [-0.1, 0, 0.1]).
>     vol_shocks: iterable of shocks to volatility.
>     relative: if True shocks are fractional (0.1 = +10%); if False they are
>         absolute additive moves (spot in price units, vol in vol points).
>
> Returns a :class:`ScenarioGrid` whose ``pnl[i][j]`` is the change in book
> market value under ``(spot_shocks[i], vol_shocks[j])``.

## schedule

### `adjust_business_day(d, convention='following')`  _function_

> Adjust a date off weekends per a business-day convention.
>
> ``following`` rolls forward to the next weekday; ``preceding`` rolls back;
> ``modified_following`` rolls forward unless that crosses into the next month,
> in which case it rolls back. ``unadjusted`` returns the date unchanged.

### `generate_schedule(start, maturity_years, freq_months, end_of_month=False, convention='unadjusted')`  _function_

> Generate period end dates from ``start`` over ``maturity_years``.
>
> Steps ``freq_months`` at a time (1=monthly, 3=quarterly, 6=semiannual,
> 12=annual) until ``maturity_years`` is reached, then applies
> :func:`adjust_business_day` with ``convention``. Returns the list of adjusted
> period end dates (the start date itself is not included).

## sizing

### `delta_hedge_shares(book: quantforge.portfolio.Book) -> float`  _function_

> Shares of the underlying to add to zero the book's net delta.
>
> A share has delta 1, so the hedge is ``-net_delta`` shares (negative = sell).

### `gamma_neutral_quantity(book: quantforge.portfolio.Book, S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None, multiplier=1.0) -> float`  _function_

> Units of a hedge option that zero the book's net gamma.

### `kelly_fraction_binary(win_prob, win_payoff, loss_amount=1.0)`  _function_

> Kelly fraction for a binary bet.
>
> Args:
>     win_prob: probability of winning, p in (0, 1).
>     win_payoff: net amount won per unit staked on a win (the "b" in b-to-1
>         odds).
>     loss_amount: amount lost per unit staked on a loss (default 1).
>
> Returns the fraction of bankroll to wager: ``f = (p*b - q*loss) / (b*loss)``
> where ``q = 1 - p``. A non-positive result means the bet has no edge; the
> optimal stake is then zero (returned as a negative/zero fraction for the
> caller to clamp).

### `kelly_fraction_continuous(expected_excess_return, variance, fraction=1.0)`  _function_

> Continuous Kelly allocation for a normally-distributed return.
>
> For a return with mean excess ``mu`` (over the risk-free rate) and variance
> ``sigma^2``, the growth-optimal leverage is ``f* = mu / sigma^2``. Multiply
> by ``fraction`` for fractional Kelly (e.g. 0.5 for half-Kelly, which trades
> a little growth for much lower drawdown).

### `kelly_growth_rate(expected_excess_return, variance, leverage)`  _function_

> Expected log-growth rate at a given leverage (continuous Kelly).
>
> ``g(f) = f*mu - 0.5 * f^2 * sigma^2``. Maximized at the full-Kelly leverage
> ``f* = mu / sigma^2``; used to compare fractional-Kelly choices.

### `neutralize(book: quantforge.portfolio.Book, greek: str, S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None, multiplier=1.0, target: float = 0.0) -> float`  _function_

> Units of the hedge option to move ``greek`` to ``target`` (default 0).
>
> ``greek`` is one of "delta", "gamma", "vega". Returns the signed quantity
> (in option units, before the multiplier is applied to notionals): solving
> ``net_greek + qty * multiplier * hedge_greek = target``.

### `vega_neutral_quantity(book: quantforge.portfolio.Book, S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None, multiplier=1.0) -> float`  _function_

> Units of a hedge option that zero the book's net vega.

## sobol

### `Sobol(dim: int)`  _class_

> A Sobol sequence generator (Gray-code recurrence).

### `brownian_bridge_path(unifs: List[float], t: float)`  _function_

> Build a Brownian path W_0..W_n at times k*dt from Sobol uniforms.
>
> ``unifs`` has length ``n`` (one per time step). The endpoint is drawn from
> the first coordinate, then successive midpoints from the rest -- so the
> leading (most uniform) Sobol dimensions carry the dominant variance.
> Returns the list ``[W_1, ..., W_n]`` (W_0 = 0).

### `sobol_arithmetic_asian_rqmc(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None, control_variate=True, n_steps=6, n_paths=4096, n_rand=24, seed=None) -> quantforge.montecarlo.MCResult`  _function_

> Randomized-QMC arithmetic Asian with a geometric control variate.
>
> Combines the two strongest variance-reduction techniques available here for a
> path-dependent payoff: low-discrepancy (randomized-QMC) sampling *and* the
> exact discrete-geometric-Asian control variate. With ``control_variate=True``
> each path's estimator is ``arith - geo + E[geo]``, where ``E[geo]`` is the
> closed-form :func:`quantforge.montecarlo._discrete_geometric_asian` over the
> same ``n_steps`` dates and the arithmetic and geometric averages (almost
> perfectly correlated) share the path. Normals come from one Sobol point
> through the Brownian bridge, randomized by a per-dimension Cranley-Patterson
> rotation, so ``n_rand`` shifts give a genuine SE.
>
> Cross-checks the control-variate :func:`quantforge.arithmetic_asian_mc`.
> ``n_steps`` is capped by the Sobol generator's dimension.

### `sobol_asian(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None, n_steps=6, n_paths=8192)`  _function_

> QMC arithmetic-average Asian price with Sobol + a Brownian bridge.
>
> Each path's ``n_steps`` normals come from one Sobol point mapped through the
> Brownian-bridge construction, so the dominant path variance lands on the
> leading (most uniform) Sobol dimensions. ``n_steps`` is capped by the
> generator's dimension.

### `sobol_asian_rqmc(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None, n_steps=6, n_paths=4096, n_rand=24, seed=None) -> quantforge.montecarlo.MCResult`  _function_

> Randomized-QMC arithmetic Asian price with an honest standard error.
>
> The multi-dimensional analogue of :func:`sobol_european_rqmc`. Each path's
> ``n_steps`` normals come from one ``n_steps``-dimensional Sobol point mapped
> through the Brownian bridge (dominant variance on the leading, most uniform
> coordinates). A per-dimension Cranley-Patterson rotation -- shift each Sobol
> coordinate by an independent ``U ~ Uniform[0,1)`` modulo 1 -- randomizes the
> whole point set without disturbing its low discrepancy, so ``n_rand``
> independent shifts give i.i.d. QMC estimates whose spread is a genuine SE.
>
> Returns an :class:`MCResult` with the mean price, the across-randomization
> SE, and ``n_paths`` = total points (``n_rand * n_paths``). Cross-checks the
> geometric-control-variate :func:`quantforge.arithmetic_asian_mc`.

### `sobol_autocallable_rqmc(S, t, r, sigma, observation_times, autocall_barrier, coupon, protection_barrier=None, notional=1.0, b=None, n_paths=4096, n_rand=24, seed=None) -> quantforge.montecarlo.MCResult`  _function_

> Randomized-QMC autocallable structured note with an honest standard error.
>
> The same product as :func:`quantforge.autocallable_mc`: at each observation
> date, if the spot is at or above ``autocall_barrier`` the note redeems early
> paying ``notional (1 + coupon k)`` (``k`` = observation number), discounted;
> if it never autocalls, at maturity the holder gets the notional back unless
> the spot finished below ``protection_barrier`` (a down-and-in put on the
> notional), taking ``notional S_T / S`` instead.
>
> Each path's ``len(observation_times)`` Brownian values come from one Sobol
> point through the Brownian bridge (dominant variance on the leading, most
> uniform coordinates), and a per-dimension Cranley-Patterson rotation
> randomizes the point set, so ``n_rand`` shifts give a genuine SE. The number
> of observations is capped by the Sobol generator's dimension. Cross-checks
> :func:`quantforge.autocallable_mc`.

### `sobol_average_strike_rqmc(S, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None, n_steps=6, n_paths=4096, n_rand=24, seed=None) -> quantforge.montecarlo.MCResult`  _function_

> Randomized-QMC average-strike Asian option with an honest standard error.
>
> The strike is the realized arithmetic average of the monitored path, so a
> call pays ``max(S_T - A, 0)`` and a put ``max(A - S_T, 0)`` with ``A`` the
> average over the ``n_steps`` monitoring dates. Normals come from one
> ``n_steps``-dim Sobol point through the Brownian bridge, randomized by a
> per-dimension Cranley-Patterson rotation, so ``n_rand`` shifts give a genuine
> SE. The discretely-monitored analogue of
> :func:`quantforge.average_strike_asian_mc`, which it cross-checks. ``n_steps``
> is capped by the Sobol generator's dimension.

### `sobol_barrier_digital_rqmc(S, K, H, t, r, sigma, option_type=<OptionType.CALL: 'call'>, barrier='up-in', b=None, cash=1.0, n_steps=6, n_paths=4096, n_rand=24, seed=None) -> quantforge.montecarlo.MCResult`  _function_

> Randomized-QMC barrier-contingent cash-or-nothing digital, honest SE.
>
> Pays ``cash`` at expiry iff the option finishes in the money (call
> ``S_T > K``, put ``S_T < K``) AND the barrier condition holds over the
> ``n_steps`` monitoring dates: ``up-in``/``down-in`` need the barrier touched,
> ``up-out``/``down-out`` need it untouched ("up" watches ``S >= H``, "down"
> ``S <= H``). Normals come from an ``n_steps``-dim Sobol point through the
> Brownian bridge, randomized by a per-dimension Cranley-Patterson rotation, so
> ``n_rand`` shifts give a genuine SE. The discrete analogue of
> :func:`quantforge.barrier_digital_mc`, which it cross-checks.

### `sobol_barrier_rqmc(S, K, H, t, r, sigma, option_type=<OptionType.CALL: 'call'>, barrier='down-out', b=None, rebate=0.0, n_steps=6, n_paths=4096, n_rand=24, seed=None) -> quantforge.montecarlo.MCResult`  _function_

> Randomized-QMC discretely-monitored single-barrier option with honest SE.
>
> Prices the knock-out/knock-in vanilla barrier option monitored at the
> ``n_steps`` dates. ``barrier`` is ``down-out``/``down-in``/``up-out``/
> ``up-in``; "down" watches ``S <= H``, "up" watches ``S >= H``. ``rebate`` is
> paid at expiry to killed knock-outs or never-activated knock-ins.
>
> Each path's normals come from one ``n_steps``-dimensional Sobol point through
> the Brownian bridge, and a per-dimension Cranley-Patterson rotation
> randomizes the point set, so ``n_rand`` shifts give i.i.d. QMC estimates
> whose spread is a genuine SE. This is the discretely-monitored analogue of
> :func:`quantforge.barrier_mc` with ``brownian_bridge=False`` (it does not add
> the continuity correction), and cross-checks it. ``n_steps`` is capped by the
> Sobol generator's dimension.

### `sobol_cliquet_rqmc(S, t, r, sigma, reset_times, local_cap=None, local_floor=0.0, global_cap=None, global_floor=0.0, b=None, n_paths=4096, n_rand=24, seed=None) -> quantforge.montecarlo.MCResult`  _function_

> Randomized-QMC capped cliquet (ratchet) with an honest standard error.
>
> The same product as :func:`quantforge.capped_cliquet_mc`: the payoff sums the
> per-period returns over the consecutive reset windows, each clipped to
> ``[local_floor, local_cap]``, then clips the running sum to
> ``[global_floor, global_cap]``, discounted at ``r``.
>
> Each path's Brownian motion at the reset times comes from one Sobol point via
> a bridge on the reset grid (:func:`_bridge_on_times`), and the per-period
> standardized shock is the bridge increment divided by ``sqrt(dt_i)``. A
> per-dimension Cranley-Patterson rotation randomizes the point set, so
> ``n_rand`` shifts give a genuine SE. The number of reset periods is capped by
> the Sobol generator's dimension. Cross-checks
> :func:`quantforge.capped_cliquet_mc`.

### `sobol_double_knockout_rqmc(S, K, t, r, sigma, lower, upper, option_type=<OptionType.CALL: 'call'>, b=None, rebate=0.0, n_steps=6, n_paths=4096, n_rand=24, seed=None) -> quantforge.montecarlo.MCResult`  _function_

> Randomized-QMC double-knockout (corridor) option with an honest standard error.
>
> Pays the vanilla payoff only if the spot stays strictly inside
> ``(lower, upper)`` at every one of the ``n_steps`` monitoring dates; if
> either barrier is breached it knocks out and pays the cash ``rebate`` at
> expiry. Normals come from an ``n_steps``-dim Sobol point through the Brownian
> bridge, randomized by a per-dimension Cranley-Patterson rotation, so
> ``n_rand`` shifts give a genuine SE. The discretely-monitored analogue of
> :func:`quantforge.double_knockout_mc`, which it cross-checks. ``n_steps`` is
> capped by the Sobol generator's dimension.

### `sobol_european(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None, n_paths=8192)`  _function_

> QMC European price using a 1-D Sobol sequence (single-step payoff).
>
> A direct Sobol analogue of :func:`quantforge.european_qmc`; converges faster
> than pseudo-random Monte Carlo for this smooth one-dimensional integral.

### `sobol_european_rqmc(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None, n_paths=4096, n_rand=24, seed=None) -> quantforge.montecarlo.MCResult`  _function_

> Randomized-QMC European price with an honest standard error.
>
> Plain Sobol QMC (:func:`sobol_european`) returns a single number with no
> error estimate -- the points are deterministic, so there is no variance to
> report. Randomized QMC restores an unbiased error bar by applying a
> Cranley-Patterson rotation: shift the whole Sobol point set by a random
> ``U ~ Uniform[0,1)`` modulo 1. Each shift preserves the sequence's low
> discrepancy but makes the resulting estimate an unbiased draw, so ``n_rand``
> independent shifts give ``n_rand`` i.i.d. QMC estimates whose spread is a
> genuine standard error.
>
> The returned :class:`MCResult` has ``price`` = mean over the randomizations,
> ``std_error`` = their across-randomization SE, and ``n_paths`` = the total
> points evaluated (``n_rand * n_paths``). For this smooth 1-D integral the
> RQMC SE falls off far faster than pseudo-random Monte Carlo's ``1/sqrt(N)``.
>
> Cross-checks the closed-form Black-Scholes value.

### `sobol_fixed_lookback_rqmc(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None, n_steps=6, n_paths=4096, n_rand=24, seed=None) -> quantforge.montecarlo.MCResult`  _function_

> Randomized-QMC discrete fixed-strike lookback with an honest standard error.
>
> Prices the discretely-monitored fixed-strike lookback -- call payoff
> ``max(max_i S_{t_i} - K, 0)``, put payoff ``max(K - min_i S_{t_i}, 0)`` --
> where the running extreme is taken over ``S_0`` and the ``n_steps``
> monitoring dates. Normals come from an ``n_steps``-dimensional Sobol point
> through the Brownian bridge, randomized by a per-dimension Cranley-Patterson
> rotation so ``n_rand`` shifts give a genuine SE.
>
> Discrete monitoring *under*-prices the continuously-monitored Conze-
> Viswanathan :func:`quantforge.fixed_strike_lookback` (fewer sampling dates
> see less extreme highs/lows); the gap shrinks as ``n_steps`` grows.
> ``n_steps`` is capped by the Sobol generator's dimension.

### `sobol_geometric_asian_rqmc(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None, n_steps=6, n_paths=4096, n_rand=24, seed=None) -> quantforge.montecarlo.MCResult`  _function_

> Randomized-QMC geometric-average Asian option with an honest standard error.
>
> Prices the discretely-monitored fixed-strike geometric-average Asian -- call
> ``max(G - K, 0)``, put ``max(K - G, 0)`` with ``G`` the geometric mean of the
> ``n_steps`` monitored spots. Normals come from one ``n_steps``-dim Sobol point
> through the Brownian bridge, randomized by a per-dimension Cranley-Patterson
> rotation, so ``n_rand`` shifts give a genuine SE.
>
> Because the discrete geometric average is exactly lognormal, this has a
> *closed form* -- :func:`quantforge.montecarlo._discrete_geometric_asian` --
> so it is the tightest available cross-check (a deterministic reference, not
> another simulation). ``n_steps`` is capped by the Sobol generator's
> dimension.

### `sobol_lookback_rqmc(S, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None, n_steps=6, n_paths=4096, n_rand=24, seed=None) -> quantforge.montecarlo.MCResult`  _function_

> Randomized-QMC floating-strike lookback with an honest standard error.
>
> Prices the discretely-monitored floating-strike lookback -- call payoff
> ``S_T - min_i S_{t_i}``, put payoff ``max_i S_{t_i} - S_T`` -- where the
> running extreme is taken over the ``n_steps`` monitoring dates (plus the
> known ``S_0``). Each path's normals come from one ``n_steps``-dimensional
> Sobol point through the Brownian bridge, and a per-dimension
> Cranley-Patterson rotation randomizes the point set, so ``n_rand`` shifts
> give i.i.d. QMC estimates whose spread is a genuine SE.
>
> Discrete monitoring always *under*-prices the continuously-monitored
> Goldman-Sosin-Gatto :func:`quantforge.floating_strike_lookback` (fewer
> sampling dates see less extreme highs/lows); the gap shrinks as ``n_steps``
> grows. ``n_steps`` is capped by the Sobol generator's dimension. Returns an
> :class:`MCResult` with the mean price, across-randomization SE, and
> ``n_paths`` = total points.

### `sobol_parisian_rqmc(S, K, H, t, r, sigma, window, option_type=<OptionType.CALL: 'call'>, barrier='down-out', b=None, n_steps=12, n_paths=4096, n_rand=24, seed=None) -> quantforge.montecarlo.MCResult`  _function_

> Randomized-QMC Parisian barrier option with an honest standard error.
>
> A Parisian barrier triggers only if the spot stays on the barrier's far side
> for a *consecutive* elapsed time of at least ``window`` years, so it is
> robust to brief spikes. ``barrier`` is ``down-out``/``down-in``/``up-out``/
> ``up-in`` ("down" watches ``S <= H``, "up" ``S >= H``). Normals come from one
> ``n_steps``-dim Sobol point through the Brownian bridge, randomized by a
> per-dimension Cranley-Patterson rotation, so ``n_rand`` shifts give a genuine
> SE. The discretely-monitored analogue of
> :func:`quantforge.parisian_barrier_mc`, which it cross-checks. ``n_steps`` is
> capped by the Sobol generator's dimension (now 12), so the window is resolved
> to ``round(window / dt)`` consecutive steps.

## spline

### `CubicSpline(xs: Sequence[float], ys: Sequence[float])`  _class_

> Natural cubic spline through ``(xs, ys)`` with ``xs`` strictly increasing.

### `SmileSpline(strikes: Sequence[float], vols: Sequence[float])`  _class_

> Strike -> implied-vol natural cubic spline with flat extrapolation.

## ssvi

### `SSVIParams(rho: float, eta: float, gamma: float, thetas: Dict[float, float] = <factory>) -> None`  _class_

> A fitted SSVI surface: global (rho, eta, gamma) and per-expiry theta.

### `calibrate_ssvi(market: Sequence[Tuple[float, float, float]], initial: quantforge.ssvi.SSVIParams = None, max_iter: int = 8000, arb_weight: float = 0.0, vega_weighted: bool = False) -> Tuple[quantforge.ssvi.SSVIParams, float]`  _function_

> Fit an SSVI surface to market implied vols.
>
> ``market`` is a sequence of ``(t, k, iv)`` points (expiry in years,
> log-moneyness, Black-Scholes implied vol). Fits the global ``(rho, eta,
> gamma)`` and one ``theta`` per distinct expiry by minimising the total-
> variance RMSE, using a smooth reparametrization so ``theta > 0``, ``eta > 0``,
> ``gamma in (0, 1)`` and ``rho in (-1, 1)`` and the built-in Nelder-Mead.
>
> Returns ``(params, rmse)`` where ``rmse`` is the root-mean-square implied-vol
> error across the market points.

### `calibrate_ssvi_arbitrage_free(market: Sequence[Tuple[float, float, float]], initial: quantforge.ssvi.SSVIParams = None, max_iter: int = 8000, weights=(1.0, 10.0, 100.0, 1000.0, 10000.0, 100000.0, 1000000.0)) -> Tuple[quantforge.ssvi.SSVIParams, float]`  _function_

> Fit an arbitrage-free SSVI surface by ramping the no-arb penalty.
>
> Calibrates with :func:`calibrate_ssvi` at increasing ``arb_weight`` values,
> warm-starting each from the previous fit, and returns the first result that
> passes :func:`ssvi_is_arbitrage_free` (or the last, most-penalised fit if
> none does). Trades a little fit RMSE for a guaranteed no-arbitrage surface.

### `ssvi_bkm_moments(params: quantforge.ssvi.SSVIParams, t, S0, r, q=0.0, n_strikes=401, width=8.0)`  _function_

> Risk-neutral (variance, skewness, excess kurtosis) of an SSVI slice at ``t``.
>
> Applies Bakshi-Kapadia-Madan moment replication
> (:func:`quantforge.bkm_moments_from_smile`) to the surface's smile at the
> fitted expiry ``t``. A negative surface ``rho`` (equity skew) yields negative
> risk-neutral skewness; ``t`` must be a fitted expiry.

### `ssvi_butterfly_free(theta: float, rho: float, eta: float, gamma: float, tol: float = 1e-09) -> bool`  _function_

> Gatheral-Jacquier sufficient condition for a butterfly-arbitrage-free slice.
>
> A fixed-``theta`` SSVI slice has no butterfly (density-negative) arbitrage if
>
>     theta * phi * (1 + |rho|) <= 4      and
>     theta * phi^2 * (1 + |rho|) <= 4 .

### `ssvi_calendar_free(params: quantforge.ssvi.SSVIParams, ks: Sequence[float] = None, tol: float = 1e-09) -> bool`  _function_

> Check the surface has no calendar-spread arbitrage on the fitted expiries.
>
> Calendar arbitrage is absent when total variance is non-decreasing in
> maturity at every log-moneyness: ``w(k, t_{i+1}) >= w(k, t_i)``. Checked on a
> grid of ``k`` across each adjacent pair of the fitted expiries.

### `ssvi_density(params: quantforge.ssvi.SSVIParams, t, S0, r, K, q=0.0, dK=None)`  _function_

> Breeden-Litzenberger risk-neutral density ``g(K)`` of an SSVI slice at ``t``.
>
> ``g(K) = e^{rt} d^2C/dK^2`` with the call priced at the surface's smile vol on
> the forward ``F = S0 e^{(r-q)t}``. Non-negative wherever the slice is
> butterfly-arbitrage-free (:func:`ssvi_butterfly_free`); ``t`` must be fitted.

### `ssvi_is_arbitrage_free(params: quantforge.ssvi.SSVIParams, ks: Sequence[float] = None) -> bool`  _function_

> True if every slice is butterfly-free and the surface is calendar-free.

### `ssvi_local_variance(k, t, theta, dtheta_dt, rho, eta, gamma)`  _function_

> Dupire local variance of an SSVI surface, fully analytic (Gatheral).
>
> Given the ATM total variance ``theta = theta(t)`` and its time derivative
> ``dtheta_dt = theta'(t)`` at maturity ``t``, the local variance at
> log-moneyness ``k`` is
>
>     sigma_loc^2 = (dw/dt)
>         / [ 1 - (k/w) w_k + (1/4)(-1/4 - 1/w + k^2/w^2) w_k^2 + (1/2) w_kk ]
>
> with all ``w`` derivatives taken in closed form from the SSVI parametrization
> (no finite differences). ``dw/dt = (dw/dtheta) * theta'(t)``.

### `ssvi_local_vol_fn(params: quantforge.ssvi.SSVIParams, S0, r, q=0.0)`  _function_

> Build a ``(spot, tau) -> local vol`` callable from a fitted SSVI surface.
>
> Converts the running spot and elapsed time into the SSVI log-moneyness
> ``k = log(spot / F_tau)`` on the forward ``F_tau = S0 e^{(r - q) tau}`` and
> returns the analytic Dupire local vol there. Suitable as the ``local_vol_fn``
> argument of :func:`quantforge.local_vol_mc`, which simulates the surface.
>
> Below the shortest fitted expiry the ATM total variance is linearly
> extrapolated toward the origin (``theta -> 0`` as ``tau -> 0``, matching
> SSVI's small-time behaviour) rather than frozen, which is what a Monte Carlo
> path integrating ``tau`` from 0 needs; the long end is clamped to the last
> fitted maturity.

### `ssvi_local_vol_from_params(params: quantforge.ssvi.SSVIParams, k, t)`  _function_

> Local volatility of a fitted SSVI surface at ``(k, t)``.
>
> Builds ``theta(t)`` and ``theta'(t)`` by linear interpolation of the fitted
> per-expiry ATM total variances (piecewise-linear in ``t``), then applies the
> analytic :func:`ssvi_local_variance`. ``t`` must lie within the fitted expiry
> range.

### `ssvi_phi(theta: float, eta: float, gamma: float) -> float`  _function_

> Power-law SSVI skew function phi(theta).

### `ssvi_reprice_mc(params: quantforge.ssvi.SSVIParams, S0, K, t, r, q=0.0, option_type=None, n_steps=100, n_paths=60000, antithetic=True, seed=None)`  _function_

> Monte Carlo a vanilla under the SSVI local-vol surface it calibrates to.
>
> Closes the calibrate -> local-vol -> reprice loop: simulates the analytic
> SSVI Dupire surface (via :func:`quantforge.local_vol_mc`) and returns the
> option's :class:`~quantforge.MCResult`. A correct local-vol construction
> reprices the SSVI *implied* smile, so this Monte Carlo price should match the
> closed-form Black-Scholes price at the SSVI implied vol for that strike.

### `ssvi_svix(params: quantforge.ssvi.SSVIParams, t, S0, r, q=0.0, n_strikes=201, width=6.0)`  _function_

> Martin (2013) SVIX index of an SSVI slice at fitted expiry ``t`` (``100 * SVIX``).

### `ssvi_total_variance(k: float, theta: float, rho: float, eta: float, gamma: float) -> float`  _function_

> SSVI total implied variance w(k, theta).

### `ssvi_variance_swap_strike(params: quantforge.ssvi.SSVIParams, t, S0, r, q=0.0, n_strikes=401, width=8.0)`  _function_

> Fair variance-swap strike (annualized *variance*) of an SSVI slice at ``t``.
>
> Replicates the variance swap from the surface's smile at the fitted expiry
> ``t``: each strike carries the Black vol ``implied_vol(ln(K/F), t)`` on the
> forward ``F = S0 e^{(r-q)t}``, fed to
> :func:`quantforge.variance_swap_from_smile`. Returns the fair *variance*
> (``sqrt`` it back to vol); ``t`` must be a fitted expiry.

### `ssvi_vix(params: quantforge.ssvi.SSVIParams, t, S0, r, q=0.0, n_strikes=201, width=6.0)`  _function_

> VIX-style index (``~= 100 * sigma``) of an SSVI slice at fitted expiry ``t``.

## strategy

### `backspread(S, K_short, K_long, t, r, sigma, kind='call', ratio=2, b=None, mult=1.0)`  _function_

> Backspread: short 1 option at K_short, long ``ratio`` at K_long.
>
> The mirror of a ratio spread -- net long options, so it profits from a large
> move (unlimited upside for a call backspread) and loses a little in the
> middle. Returns the leg :class:`Book`.

### `box_spread(S, K_low, K_high, t, r, sigma, b=None, mult=1.0)`  _function_

> Box spread: a bull call spread plus a bear put spread on the same strikes.
>
> Long call ``K_low`` / short call ``K_high`` (bull call) combined with long
> put ``K_high`` / short put ``K_low`` (bear put). The terminal payoff is the
> constant ``K_high - K_low`` regardless of spot, so the box is a synthetic
> zero-coupon bond: its fair value is the discounted strike width
> ``e^{-rt} (K_high - K_low)``, which the net Greeks confirm are ~zero in spot.
> Returns the leg :class:`Book`.

### `break_evens(book: quantforge.portfolio.Book, lo: float, hi: float, n: int = 2000) -> List[float]`  _function_

> Find terminal spots where total P&L (payoff - net premium) crosses zero.
>
> Scans ``[lo, hi]`` on a grid and refines each sign change by bisection. Net
> premium is the book's market value now (positive = we paid it).

### `butterfly(S, K_low, K_mid, K_high, t, r, sigma, kind='call', b=None, mult=1.0)`  _function_

> Long butterfly: +1 K_low, -2 K_mid, +1 K_high (equally spaced strikes).

### `calendar_spread(S, K, t_near, t_far, r, sigma, kind='call', b=None, mult=1.0)`  _function_

> Calendar (horizontal) spread: short the near expiry, long the far, same
> strike ``K``.
>
> A long calendar is short one near-dated option and long one far-dated option
> at the same strike, financed by the faster time decay of the near leg. The
> two legs carry different maturities (``t_near < t_far``), which the
> per-contract ``t`` supports, so ``price_book`` gives the net debit and the
> net Greeks directly. (The expiry payoff diagram is not well defined by
> intrinsics alone, since the far leg still has time value at the near expiry;
> use the net Greeks and price for analysis.) Returns the leg :class:`Book`.

### `collar(S, K_put, K_call, t, r, sigma, b=None, mult=1.0)`  _function_

> Protective collar on a long share: long a put at ``K_put`` (floor) and
> short a call at ``K_call`` (cap), with ``K_put < K_call``.
>
> Priced here as the two option legs (the underlying share is held
> separately); the net option premium is a small debit or credit depending on
> the skew. The collar caps gains above ``K_call`` and floors losses below
> ``K_put``. Returns the leg :class:`Book`.

### `diagonal_spread(S, K_near, K_far, t_near, t_far, r, sigma, kind='call', b=None, mult=1.0)`  _function_

> Diagonal spread: short the near expiry at ``K_near``, long the far expiry
> at ``K_far`` -- a calendar with different strikes on the two legs.
>
> Combines the horizontal (time) and vertical (strike) spreads. Requires
> ``t_near < t_far``; strikes may differ freely. Returns the leg :class:`Book`;
> net price and Greeks come from ``price_book`` (see :func:`calendar_spread` on
> the expiry-payoff caveat).

### `iron_condor(S, K_put_long, K_put_short, K_call_short, K_call_long, t, r, sigma, b=None, mult=1.0)`  _function_

> Iron condor: sell an OTM put spread and an OTM call spread.
>
> Strikes ordered K_put_long < K_put_short < K_call_short < K_call_long.
> Collects premium; profits if the underlying stays between the short strikes.

### `payoff_at_expiry(book: quantforge.portfolio.Book, spot_at_expiry: float) -> float`  _function_

> Intrinsic payoff of the book's legs at a terminal spot (per multiplier).

### `payoff_profile(book: quantforge.portfolio.Book, spots: Sequence[float]) -> List[float]`  _function_

> Payoff at each terminal spot in ``spots``.

### `ratio_spread(S, K_long, K_short, t, r, sigma, kind='call', ratio=2, b=None, mult=1.0)`  _function_

> Ratio spread: long 1 option at K_long, short ``ratio`` at K_short.
>
> A call ratio spread (K_long < K_short, ratio > 1) is long one lower-strike
> call and short several higher-strike calls -- typically a small credit or
> debit with a capped-profit tent that turns into unlimited downside beyond
> the short strikes. Returns the leg :class:`Book`.

### `risk_reversal(S, K_put, K_call, t, r, sigma, b=None, mult=1.0)`  _function_

> Risk reversal: short an OTM put, long an OTM call (a skew/forward trade).

### `straddle(S, K, t, r, sigma, b=None, mult=1.0, qty=1)`  _function_

> Long straddle: long a call and a put at the same strike.

### `strangle(S, K_put, K_call, t, r, sigma, b=None, mult=1.0, qty=1)`  _function_

> Long strangle: long an OTM put and an OTM call (K_put < K_call).

### `strategy_report(book, lo=None, hi=None, n=4000)`  _function_

> Summarize a strategy's expiry P&L: max profit, max loss, break-evens.
>
> Scans terminal spots on a grid ``[lo, hi]`` (defaults span a wide range
> around the leg strikes) and returns a dict with the net premium, the maximum
> profit and maximum loss seen on the grid (P&L = payoff - premium), whether
> each is bounded (i.e. not still rising/falling at the grid edge), and the
> break-even spots from :func:`break_evens`.
>
> P&L is per multiplier, matching :func:`payoff_at_expiry`.

### `synthetic_forward(S, K, t, r, sigma, b=None, mult=1.0)`  _function_

> Synthetic long forward: long a call and short a put at the same strike.
>
> By put-call parity the position replicates a forward struck at ``K``: its
> present value is ``C - P = e^{-bt} S - e^{-rt} K`` (carry ``b``), its delta is
> ~1, and its gamma/vega net to ~zero. Returns the leg :class:`Book`.

### `vertical_spread(S, K_long, K_short, t, r, sigma, kind='call', b=None, mult=1.0)`  _function_

> Bull/bear vertical: long one option at K_long, short one at K_short.
>
> A call spread with K_long < K_short is a bull spread; a put spread with
> K_long > K_short is a bear spread.

## structured

### `buffered_note(S, K, buffer, maturity, r, sigma, principal, n_shares=None, b=None)`  _function_

> Buffered note: absorbs the first ``buffer`` fraction of downside losses.
>
> The investor is short a put struck at the buffered level ``K * (1 - buffer)``
> rather than at ``K``, so losses only bite once the underlying falls more than
> ``buffer``. Value is ``PV(principal) - n_shares * put(K*(1-buffer))``. A larger
> buffer moves the put further out of the money, raising the note's value.

### `capped_principal_protected_note(S, K, cap_level, maturity, r, sigma, principal, participation=1.0, n_shares=None, b=None)`  _function_

> Principal-protected note with a capped upside (a call spread).
>
> Like :func:`principal_protected_note` but the upside is a call spread -- long
> a call at ``K``, short a call at ``cap_level`` -- so the payoff is capped once
> the underlying passes ``cap_level``. Value is the ZC bond plus
> ``participation * n_shares * (call(K) - call(cap_level))``. At or below the
> uncapped PPN (selling the higher-strike call raises no value).

### `note_embedded_option_value(note_value, principal, r, maturity)`  _function_

> Option component of a note: ``note_value - discounted principal``.
>
> Strips the guaranteed bond leg to isolate the value attributable to the
> embedded option position (positive for a PPN's long call, negative for a
> reverse convertible's short put).

### `note_zero_coupon_bond(principal, r, maturity)`  _function_

> Present value of a zero-coupon bond: ``principal * e^{-r T}``.

### `principal_protected_note(S, K, maturity, r, sigma, principal, participation=1.0, n_shares=None, b=None)`  _function_

> Value of a principal-protected note = ZC bond + participation * call.
>
> Guarantees the ``principal`` at maturity (the ZC-bond leg) and adds
> ``participation`` times a call on the underlying for upside. ``n_shares``
> scales the call to the note's notional (defaults to ``principal / S``, i.e. the
> note buys as many shares as the principal affords at inception). The value is
> always at least the discounted principal (the call leg is non-negative).

### `reverse_convertible(S, K, maturity, r, sigma, principal, coupon_rate, n_shares=None, b=None)`  _function_

> Value of a reverse convertible = ZC bond + coupon PV - short put.
>
> The investor receives an enhanced coupon and the principal, but is short a put
> struck at ``K`` (delivering shares if the underlying falls). Value is
> ``PV(principal + coupon) - n_shares * put``. ``n_shares`` defaults to
> ``principal / K`` (the put covers the principal at the strike). Below the
> plain bond-plus-coupon value because of the short put.

### `reverse_convertible_fair_coupon(S, K, maturity, r, sigma, principal, n_shares=None, b=None, tol=1e-12, max_iter=100)`  _function_

> Coupon rate that prices a reverse convertible at par (its principal).
>
> Solves :func:`reverse_convertible` ``= principal`` for the ``coupon_rate``.
> The short put costs value, so the fair coupon is positive -- the enhanced yield
> that compensates the investor for the downside they sell. Closed form:
> ``coupon = (put_value / disc / principal ... )``; here solved directly since
> the note is linear in the coupon.

## surface

### `CalendarViolation(t_short: float, t_long: float, k: float, w_short: float, w_long: float) -> None`  _class_

> CalendarViolation(t_short: float, t_long: float, k: float, w_short: float, w_long: float)

### `SurfaceSlice(t: float, params: quantforge.svi.SVIParams, rmse: float) -> None`  _class_

> SurfaceSlice(t: float, params: quantforge.svi.SVIParams, rmse: float)

### `VolSurface(slices: List[quantforge.surface.SurfaceSlice])`  _class_

> A term structure of SVI smiles with calendar-arbitrage diagnostics.

## svi

### `SVIParams(a: float, b: float, rho: float, m: float, s: float) -> None`  _class_

> SVIParams(a: float, b: float, rho: float, m: float, s: float)

### `calibrate_svi(ks: Sequence[float], total_variances: Sequence[float], weights: Sequence[float] = None, initial: quantforge.svi.SVIParams = None, max_iter: int = 4000) -> Tuple[quantforge.svi.SVIParams, float]`  _function_

> Fit raw SVI to observed (log-moneyness, total-variance) points.
>
> Returns ``(params, rmse)`` where rmse is the root-mean-square total-variance
> error. Uses an unconstrained Nelder-Mead over a smooth reparametrization
> that enforces ``b >= 0``, ``s > 0``, and ``rho in (-1, 1)``.

### `calibrate_svi_from_prices(F, t, r, strikes, call_prices, q=0.0, vega_weighted=True, initial=None, max_iter=4000)`  _function_

> Calibrate a raw SVI slice directly from market *call prices*.
>
> Inverts each call to its Black-Scholes implied volatility, converts to total
> variance ``w = sigma^2 t``, and fits raw SVI with :func:`calibrate_svi`.
> Quotes are vega-weighted by default (near-the-money prices carry the most
> volatility information, so weighting by Black vega down-weights the deep
> wings where a price error maps to a large vol error).
>
> Args:
>     F: forward. strikes, call_prices: matching market quotes at expiry ``t``.
>     r: discount rate (the calls are priced on the forward, carry ``b = r``
>         relative to spot ``S = F e^{-rt}``... here calls are taken on the
>         forward directly with discounting ``e^{-rt}``).
>     vega_weighted: weight each quote by its Black vega if True.
>
> Returns ``(params, iv_rmse, price_rmse)``.

### `lee_bounds_ok(p: quantforge.svi.SVIParams, tol=1e-09)`  _function_

> True if both SVI wing slopes satisfy Lee's ``slope <= 2`` moment bound.
>
> Equivalent to :meth:`SVIParams.is_arbitrage_free_wings` but exposes the two
> directional slopes explicitly via :func:`lee_wing_slopes`.

### `lee_wing_slopes(p: quantforge.svi.SVIParams)`  _function_

> Asymptotic wing slopes of total variance for a raw-SVI slice.
>
> As ``k -> +/- inf`` the SVI total variance ``w(k)`` is linear with slopes
>
>     right (k -> +inf):  b (1 + rho)
>     left  (k -> -inf):  b (1 - rho)
>
> Lee's moment formula caps the slope of *total variance* at 2 for a valid
> (arbitrage-free-wing) surface, so both slopes must be <= 2. Returns
> ``(left_slope, right_slope)``.

### `svi_bkm_moments(p: quantforge.svi.SVIParams, S0, t, r, q=0.0, n_strikes=401, width=8.0)`  _function_

> Risk-neutral (variance, skewness, excess kurtosis) implied by an SVI slice.
>
> Maps each strike to ``p.implied_vol(ln(K/F), t)`` and applies the
> Bakshi-Kapadia-Madan moment replication
> (:func:`quantforge.bkm_moments_from_smile`). Returns
> ``(variance, skewness, excess_kurtosis)`` of the ``t``-horizon risk-neutral
> log-return. A negative SVI ``rho`` (equity skew) produces negative
> risk-neutral skewness; a flat slice is near-symmetric.

### `svi_butterfly_arbitrage(p: quantforge.svi.SVIParams, ks=None, tol=1e-10)`  _function_

> Return the log-moneyness points where the SVI slice has butterfly arb.
>
> Scans ``ks`` (default a wide grid) and reports those where ``g(k) < -tol``.
> An empty list means the slice is butterfly-arbitrage-free on the grid.

### `svi_density(p: quantforge.svi.SVIParams, S0, t, r, K, q=0.0, dK=None)`  _function_

> Breeden-Litzenberger risk-neutral density ``g(K)`` implied by an SVI slice.
>
> ``g(K) = e^{r t} d^2 C / dK^2`` with the call priced at the slice's smile vol
> ``p.implied_vol(ln(K/F), t)`` on the forward ``F = S0 e^{(r-q)t}``. Non-
> negative wherever the slice is butterfly-arbitrage-free (see
> :func:`svi_is_butterfly_free`); a negative value flags a density violation.

### `svi_g(p: quantforge.svi.SVIParams, k)`  _function_

> Gatheral-Jacquier g-function of a raw-SVI slice at log-moneyness ``k``.
>
> The slice is free of butterfly (static/density) arbitrage iff ``g(k) >= 0``
> for all ``k`` (the implied risk-neutral density is then non-negative). With
> ``w = w(k)``, ``w'`` and ``w''``:
>
>     g(k) = (1 - k w' / (2w))^2 - (w'^2 / 4)(1/w + 1/4) + w''/2.

### `svi_is_butterfly_free(p: quantforge.svi.SVIParams, ks=None) -> bool`  _function_

> (no docstring)

### `svi_local_variance(p: quantforge.svi.SVIParams, k, dw_dt)`  _function_

> Dupire local variance of a single SVI slice, analytic in strike.
>
> Given the slice ``p`` and the total-variance time derivative ``dw_dt =
> dw/dt`` at log-moneyness ``k`` (supplied by the caller, since one slice
> carries no maturity information), the Gatheral total-variance Dupire formula
> gives
>
>     sigma_loc^2 = dw/dt
>         / [ 1 - (k/w) w_k + (1/4)(-1/4 - 1/w + k^2/w^2) w_k^2 + (1/2) w_kk ]
>
> with ``w``, ``w_k = w'(k)`` and ``w_kk = w''(k)`` taken in closed form from
> the SVI parametrization (no finite differences in strike). Raises if the
> Dupire denominator is non-positive (a butterfly-arbitrage flag: the slice's
> ``svi_g`` is negative there).

### `svi_repair_butterfly(p: quantforge.svi.SVIParams, ks=None, max_iter=200, factor=0.98)`  _function_

> Repair a single SVI slice's butterfly arbitrage by shrinking the wings.
>
> If the slice has ``svi_g(k) < 0`` anywhere (a negative density), it reduces
> the wing angle ``b`` geometrically (which flattens the smile and lifts the
> g-function) until :func:`svi_is_butterfly_free` passes or ``max_iter`` is
> reached. Returns a new :class:`SVIParams`; the ATM level, skew, shift and
> curvature are preserved. If already arbitrage-free the input is returned
> unchanged.

### `svi_surface_local_vol(slices, k, t)`  _function_

> Local volatility from a term structure of SVI slices at ``(k, t)``.
>
> ``slices`` maps expiry ``t_i`` (years) to a fitted :class:`SVIParams`. Total
> variance is interpolated *linearly in t* at fixed ``k`` to supply the
> Dupire ``dw/dt`` (the piecewise-constant slope of the bracketing slices),
> while the strike derivatives come analytically from the slice active at
> ``t``. ``t`` must lie within the fitted expiry range.
>
> Returns the local volatility ``sqrt(sigma_loc^2)``.

### `svi_svix(p: quantforge.svi.SVIParams, S0, t, r, q=0.0, n_strikes=201, width=6.0)`  _function_

> Martin (2013) SVIX index implied by a raw-SVI slice.
>
> Maps each strike to ``p.implied_vol(ln(K/F), t)`` and feeds the smile to
> :func:`quantforge.svix_from_smile` (the ``1/F^2``-weighted, put-call-symmetric
> variance index that lower-bounds the equity premium), reported as
> ``100 * SVIX``. A flat slice returns approximately ``100 * sigma``.

### `svi_variance_swap_strike(p: quantforge.svi.SVIParams, S0, t, r, q=0.0, n_strikes=401, width=8.0)`  _function_

> Fair variance-swap strike (annualized *variance*) implied by a raw-SVI slice.
>
> Replicates the variance swap from the SVI smile: at each strike the Black
> implied vol is ``p.implied_vol(k, t)`` with ``k = ln(K / F)`` the
> log-moneyness on the forward ``F = S0 e^{(r-q)t}``. Feeds the smile to
> :func:`quantforge.variance_swap_from_smile`, so the result is model-
> consistent with the fitted slice. Returns the fair *variance* (square it back
> to vol with ``sqrt``); a flat slice (``b = 0``) returns that flat variance
> ``sigma^2``, and a skewed slice returns a variance above the ATM variance (the
> convexity/skew premium).

### `svi_vix(p: quantforge.svi.SVIParams, S0, t, r, q=0.0, n_strikes=201, width=6.0)`  _function_

> VIX-style index (``~= 100 * sigma``) implied by a raw-SVI slice.
>
> Maps each strike to ``p.implied_vol(ln(K/F), t)`` and feeds the smile to
> :func:`quantforge.vix_from_smile`. A flat slice returns ``100 * sigma``.

## trinomial

### `richardson_american(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None, steps=200)`  _function_

> Richardson-extrapolated American price from ``n`` and ``2n`` trinomial solves.
>
> The trinomial price converges to the true value with a leading O(1/n) error,
> so ``2 * P(2n) - P(n)`` cancels that term and converges faster. Returns the
> extrapolated price.

### `trinomial_price(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None, steps=200, american=True)`  _function_

> Price an option on a Boyle trinomial lattice.
>
> Args:
>     american: if True, allow early exercise at every node; if False, price
>         the European payoff (useful as a convergence cross-check).
>     b: cost of carry (defaults to r). Dividend yield q enters as b = r - q.

## vannavolga

### `VannaVolgaSmile(S, t, r_dom, r_for, atm, rr, bf, call_delta=0.25)`  _class_

> A vanna-volga FX smile built from ATM / RR / BF at one expiry.

### `pillar_vols(atm, rr, bf)`  _function_

> Return (sigma_25put, sigma_atm, sigma_25call) from ATM / RR / BF quotes.

## variancegamma

### `variance_gamma_greeks(S, K, t, r, sigma, nu, theta, option_type=<OptionType.CALL: 'call'>, q=0.0, cm_alpha=1.5)`  _function_

> Greeks of a Variance-Gamma option by central finite differences.
>
> Central differences of :func:`variance_gamma_price` for ``delta`` (dV/dS),
> ``gamma`` (d2V/dS2), ``vega`` (dV/dsigma, the Brownian-vol sensitivity), and
> ``theta_greek`` (calendar decay, ``-dV/dt``). As ``nu -> 0`` the Greeks
> approach the Black-Scholes Greeks. ``theta`` is the VG skew *parameter*; the
> calendar Greek is returned as ``theta_greek`` to avoid the name clash.
> Returns a dict with ``price``, ``delta``, ``gamma``, ``vega``, ``theta_greek``.

### `variance_gamma_price(S, K, t, r, sigma, nu, theta, option_type=<OptionType.CALL: 'call'>, q=0.0, cm_alpha=1.5, upper=200.0)`  _function_

> Price a European option under the Variance-Gamma model.
>
> Args:
>     sigma: Brownian volatility. nu: gamma-time variance rate (> 0).
>     theta: Brownian drift (skew; negative for an equity left skew).
>     q: continuous dividend yield.
>     cm_alpha: Carr-Madan damping; the transform needs
>         ``1 - theta nu (cm_alpha+1) - 0.5 sigma^2 nu (cm_alpha+1)^2 > 0``.
>
> Puts follow from put-call parity. As ``nu -> 0`` the price approaches the
> Black-Scholes value.

### `variance_gamma_smile(S, strikes, t, r, sigma, nu, theta, q=0.0, cm_alpha=1.5)`  _function_

> Black-Scholes implied-vol smile the Variance-Gamma model produces.
>
> Prices a call at each strike and inverts to a Black-Scholes implied vol,
> returning ``(log_moneyness, vol)`` pairs sorted by strike on the forward
> ``F = S e^{(r-q) t}``. ``theta < 0`` tilts the smile into a downward skew;
> larger ``nu`` fattens the wings.

## varswap

### `corridor_variance_swap_from_smile(S0, t, r, vol_fn, lower, upper, q=0.0, n_strikes=401, split=None)`  _function_

> Fair corridor variance-swap strike from a smile ``vol_fn(K)``.
>
> A corridor variance swap accrues realized variance only while the spot is in
> the corridor ``[lower, upper]``. By the Carr-Madan static-replication view
> this restricts the ``1/K^2``-weighted option strip to strikes inside the
> corridor (Carr & Lewis): the fair accrued variance is
>
>     K_corr = (2 e^{r t} / t) * ( strip of OTM options with L <= K <= U ).
>
> ``vol_fn(K)`` prices each strip option with Black-Scholes at its smile vol.
> A corridor spanning the whole strip recovers (most of) the plain
> variance-swap strike; a narrower corridor accrues less variance.

### `forward_variance_swap_from_smile(S0, t1, t2, r, vol_fn1, vol_fn2, q=0.0, n_strikes=401, width=8.0)`  _function_

> Fair forward-start variance-swap strike over ``[t1, t2]`` from two smiles.
>
> Total (undiscounted) variance is additive in time, so the fair variance
> accrued between ``t1`` and ``t2`` is
>
>     K_fwd = ( K_var(t2) * t2 - K_var(t1) * t1 ) / (t2 - t1),
>
> where ``K_var(t_i)`` is the spot-starting variance-swap strike replicated
> from the expiry-``t_i`` smile ``vol_fn_i(K)``. Requires ``0 <= t1 < t2``. A
> flat term structure of flat smiles returns that flat variance.

### `gamma_swap_from_smile(S0, t, r, vol_fn, q=0.0, n_strikes=401, width=8.0, split=None)`  _function_

> Fair gamma-swap strike from a smile ``vol_fn(K)``.
>
> A gamma (or "weighted variance") swap accrues ``(S_t / S0) d<ln S>`` -- each
> increment of realized variance weighted by the spot level -- so it is
> replicated by an option strip weighted ``1/K`` (the price-weighted version of
> the variance swap's ``1/K^2``), plus the matching log-contract terms
> (Carr-Lewis). Its fair strike is
>
>     K_gamma = (2 e^{r t} / (S0 t)) * ( 1/K-weighted OTM strip )
>               + (2/t) * (r - q) * (e^{(r-q) t} - 1) / (r - q) ...  [drift term]
>
> Implemented from the price-weighted log contract; ``vol_fn(K)`` prices each
> option with Black-Scholes. A flat smile returns that flat variance.

### `variance_swap_from_smile(S0, t, r, vol_fn, q=0.0, n_strikes=401, width=8.0, split=None)`  _function_

> Fair variance-swap strike from a volatility *smile* ``vol_fn(K)``.
>
> Builds the OTM option strip -- puts below the forward split, calls above --
> by pricing each strike at its smile vol ``vol_fn(K)`` with Black-Scholes,
> then feeds them to :func:`variance_swap_strike`. Convenient for marking a
> variance swap directly off a fitted smile (SVI, SABR, vanna-volga, ...).
>
> Args:
>     vol_fn: callable ``vol_fn(K)`` returning the Black implied vol at strike.
>     q: dividend yield (carry ``b = r - q``).
>     n_strikes: number of strikes on each side; strikes span ``width``
>         standard deviations of log-moneyness around the split.
>
> A flat smile returns exactly that flat variance (the model-free result).

### `variance_swap_strike(S0, t, r, put_strikes: Sequence[float], put_prices: Sequence[float], call_strikes: Sequence[float], call_prices: Sequence[float], split: float = None) -> float`  _function_

> Fair variance-swap strike (annualized variance) by option replication.
>
> Args:
>     S0: current spot.
>     t: swap tenor in years.
>     r: risk-free rate.
>     put_strikes/put_prices: OTM puts, strikes strictly below ``split``.
>     call_strikes/call_prices: OTM calls, strikes strictly above ``split``.
>     split: the forward split level ``K*``. Defaults to the forward
>         ``S0 e^{r t}``.
>
> Returns the fair strike as an annualized variance (multiply tenor and take
> sqrt for a vol number).

### `variance_term_structure(S0, r, expiries, vol_fns, q=0.0, n_strikes=401, width=8.0)`  _function_

> Term structure of variance-swap strikes and the forward-variance curve.
>
> Args:
>     expiries: increasing list of expiries.
>     vol_fns: one smile ``vol_fn(K)`` per expiry (matching order).
>
> Returns ``(spot_var, forward_var)`` where ``spot_var[i]`` is the
> spot-starting variance-swap strike to ``expiries[i]`` and ``forward_var[i]``
> the annualized *forward* variance over ``(expiries[i-1], expiries[i]]``
> (``forward_var[0]`` = ``spot_var[0]``). By total-variance additivity
> ``forward_var[i] = (K_i t_i - K_{i-1} t_{i-1}) / (t_i - t_{i-1})``.
>
> A flat term structure of flat smiles gives a constant curve; a rising
> variance term structure gives positive, increasing forward variances.

### `volatility_swap_bounds_from_smile(S0, t, r, vol_fn, q=0.0, n_strikes=401, width=8.0)`  _function_

> Bracket the fair volatility-swap strike from a smile: ``(lower, upper)``.
>
> The fair vol-swap strike ``E[sqrt(RV)]`` has no model-free replication (unlike
> the variance swap), but it is pinned between two computable levels:
>
>   * upper -- ``sqrt(K_var)`` with ``K_var`` the variance-swap strike, since
>     ``E[sqrt(RV)] <= sqrt(E[RV])`` by Jensen (the convexity/vol-of-vol gap).
>   * lower -- the at-the-money-forward implied vol ``vol_fn(F)``. Carr-Lee
>     (2009) show the fair vol-swap strike equals the ATMF implied vol exactly
>     when spot/vol correlation is zero, and the ATMF vol lies below
>     ``sqrt(K_var)`` whenever the smile is convex; it is the standard
>     first-order proxy.
>
> Returns ``(atmf_vol, sqrt_var_strike)`` with ``lower <= upper``. A flat smile
> collapses the bracket to that flat vol (zero convexity). The true strike sits
> inside; the width ``upper - lower`` is the convexity premium the smile implies.

### `volatility_swap_strike(S0, t, r, put_strikes, put_prices, call_strikes, call_prices, split=None) -> float`  _function_

> Fair volatility-swap strike as ``sqrt(variance strike)``.
>
> This is the standard first-order proxy; it slightly overstates the true
> vol-swap strike because ``E[sqrt(var)] <= sqrt(E[var])`` (Jensen), the
> "convexity" or "vol-of-vol" adjustment, which requires a model to quantify.

## vasicek

### `bond_option(r0, t_option, t_bond, strike, kappa, theta, sigma, option_type=<OptionType.CALL: 'call'>)`  _function_

> European option (Jamshidian) on a zero-coupon bond under Vasicek.
>
> Args:
>     t_option: option expiry. t_bond: the underlying bond's maturity
>         (``t_bond > t_option``). strike: strike on the bond price.
>
> A call pays ``max(P(t_option, t_bond) - strike, 0)`` at the option expiry.
> Uses the closed-form bond-price volatility.

### `vasicek_bond_greeks(r0, t, kappa, theta, sigma)`  _function_

> Exact rate sensitivities of a Vasicek zero-coupon bond.
>
> ``P = A(t) e^{-B(t) r0}`` (with ``B = _B(kappa, t)``), so ``rho_r = -B P``,
> ``gamma_r = B^2 P``, rate ``duration = B``, and ``convexity = B^2``. Returns
> a dict with ``price``, ``rho_r``, ``gamma_r``, ``duration``, ``convexity``.

### `vasicek_bond_option_greeks(r0, t_option, t_bond, strike, kappa, theta, sigma, option_type=<OptionType.CALL: 'call'>)`  _function_

> Greeks of a Vasicek zero-coupon-bond option by central finite differences.
>
> Central differences of :func:`bond_option` for the short-rate sensitivities
> ``rho_r`` (dV/dr0) and ``gamma_r`` (d2V/dr0^2), and the vol sensitivity
> ``vega`` (dV/dsigma). A bond call *falls* as the short rate rises (higher
> rates discount the bond harder), so ``rho_r < 0`` for a call. Returns a dict
> with ``price``, ``rho_r``, ``gamma_r``, ``vega``.

### `vasicek_cap(r0, dates, strike, kappa, theta, sigma, notional=1.0)`  _function_

> Vasicek cap: strip of caplets over successive ``dates`` (increasing times).

### `vasicek_caplet(r0, reset, pay, strike, kappa, theta, sigma, notional=1.0)`  _function_

> Caplet on ``[reset, pay]`` under Vasicek via the bond-put identity.
>
> Pays ``tau (L - strike)^+`` at ``pay``; equals ``notional (1 + strike tau)``
> puts on the zero-coupon bond ``P(reset, pay)`` struck at ``1/(1 + strike tau)``.

### `vasicek_coupon_bond_option(r0, t_option, cashflows, strike, kappa, theta, sigma, option_type=<OptionType.CALL: 'call'>)`  _function_

> European option on a coupon bond under Vasicek (Jamshidian decomposition).
>
> ``cashflows`` is a list of ``(t_i, c_i)`` pairs with ``t_i > t_option``: the
> underlying coupon bond pays ``c_i`` at each ``t_i`` (the last usually
> includes the principal). The option pays ``max(B(t_option) - strike, 0)``
> (call) on the bond's value ``B``.
>
> Since the Vasicek short rate is one-factor and every zero-coupon bond is
> monotone decreasing in ``r``, Jamshidian's trick applies: find the critical
> rate ``r*`` where the bond value at expiry equals ``strike``, split ``strike``
> into per-cashflow strikes ``K_i = P(t_option, t_i | r*)``, and the coupon-bond
> option is the ``c_i``-weighted sum of zero-coupon-bond options struck at each
> ``K_i``. Exact (no simulation).

### `vasicek_expected_rate(r0, t, kappa, theta, sigma=0.0)`  _function_

> Expected short rate ``E[r_t] = theta + (r0 - theta) e^{-kappa t}``.
>
> The mean of the Ornstein-Uhlenbeck process; it decays from ``r0`` toward the
> long-run level ``theta`` at speed ``kappa`` (``sigma`` does not enter the
> mean and is accepted only for a uniform signature).

### `vasicek_floor(r0, dates, strike, kappa, theta, sigma, notional=1.0)`  _function_

> Vasicek floor: strip of floorlets over successive ``dates``.

### `vasicek_floorlet(r0, reset, pay, strike, kappa, theta, sigma, notional=1.0)`  _function_

> Floorlet on ``[reset, pay]`` under Vasicek via the bond-call identity.

### `vasicek_rate_variance(t, kappa, sigma)`  _function_

> Variance of the short rate ``Var[r_t] = sigma^2/(2 kappa) (1 - e^{-2 kappa t})``.
>
> Grows from 0 to the stationary variance ``sigma^2/(2 kappa)`` as ``t`` rises;
> at ``kappa -> 0`` it degenerates to the Brownian ``sigma^2 t``.

### `vasicek_stationary_distribution(kappa, theta, sigma)`  _function_

> Long-run (stationary) distribution of the short rate as ``(mean, variance)``.
>
> As ``t -> infinity`` the OU rate is Normal with mean ``theta`` and variance
> ``sigma^2 / (2 kappa)``. Requires ``kappa > 0`` (otherwise no stationary law).

### `vasicek_swaption(r0, expiry, pay_times, fixed_rate, kappa, theta, sigma, payer=True, notional=1.0)`  _function_

> European swaption under Vasicek via the coupon-bond-option identity (exact).
>
> A physically-settled European swaption is an option on the underlying swap.
> The fixed leg plus notional at maturity is a coupon bond with cashflows
> ``fixed_rate * tau_i`` at each ``pay_times[i]`` and the notional at the last
> date. Entering a *payer* swap (pay fixed, receive float) at ``expiry`` is
> worth ``notional - couponbond``, so a payer swaption is a *put* on that coupon
> bond struck at the notional, and a receiver swaption a *call* -- both priced
> exactly by :func:`coupon_bond_option` (Jamshidian), no approximation.
>
> ``pay_times`` are the fixed-leg payment dates (all ``> expiry``); accruals
> ``tau_i`` are the gaps between them, with the first gap measured from
> ``expiry``.

### `zero_coupon_bond(r0, t, kappa, theta, sigma)`  _function_

> Vasicek zero-coupon bond price P(0, t) for a unit face, given r(0)=r0.
>
> ``P = A(t) * exp(-B(t) * r0)`` with the standard affine coefficients.

### `zero_coupon_yield(r0, t, kappa, theta, sigma)`  _function_

> Continuously-compounded yield of the Vasicek zero-coupon bond to ``t``.

## vectorized

### `delta_array(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None)`  _function_

> (no docstring)

### `gamma_array(S, K, t, r, sigma, b=None)`  _function_

> (no docstring)

### `greeks_array(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None)`  _function_

> Return a dict of vectorized price, delta, gamma, vega arrays.

### `price_array(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None)`  _function_

> Vectorized European price. Inputs may be scalars or NumPy arrays.

### `vega_array(S, K, t, r, sigma, b=None)`  _function_

> (no docstring)

## vegabucket

### `VegaBuckets(buckets: Dict[str, float], total: float) -> None`  _class_

> VegaBuckets(buckets: Dict[str, float], total: float)

### `vega_buckets(contracts: Sequence[quantforge.portfolio.Contract], edges: Sequence[float] = (0.25, 0.5, 1.0, 2.0, 5.0)) -> quantforge.vegabucket.VegaBuckets`  _function_

> Bucket a book's position-scaled vega by expiry.
>
> Args:
>     contracts: the book's legs (signed qty, multiplier as in ``price_book``).
>     edges: sorted upper-edge tenors in years. A contract with expiry ``t``
>         falls in the first bucket whose edge is >= ``t``; longer expiries go
>         to the final ">last" bucket.
>
> Returns a :class:`VegaBuckets`. Vega is per 1.0 change in vol (divide by 100
> for per-vol-point), scaled by ``qty * multiplier``.

## vix

### `equity_premium_lower_bound(S0, t, r, vol_fn, q=0.0, n_strikes=201, width=6.0)`  _function_

> Martin's (2013) lower bound on the expected equity excess return.
>
> Martin shows that, under the (empirically mild) negative-correlation
> condition, the expected simple excess return of the market over ``[0, t]`` is
> bounded below by the risk-neutral *simple variance*:
>
>     (1/t) E_0[ (R_market - R_f) ] >= Rf * SVIX^2,
>
> where ``SVIX^2`` is the annualized simple-variance index
> (:func:`svix_from_smile`) and ``Rf = e^{r t}`` the gross risk-free return.
> This returns the annualized lower bound ``Rf * SVIX^2`` -- a model-free floor
> on the equity premium computable purely from option prices.

### `svix_from_smile(S0, t, r, vol_fn, q=0.0, n_strikes=201, width=6.0)`  _function_

> Martin (2013) "simple variance" index (SVIX) from a smile ``vol_fn(K)``.
>
> Unlike the VIX log-contract (``1/K^2`` weights), the simple variance swap
> weights the OTM strip by ``1/F^2`` -- a constant -- so it corresponds to the
> payoff ``(S_T - F)^2 / F^2`` and needs no log approximation, making it robust
> to large moves/jumps and giving a genuine lower bound on the equity premium
> (Martin). The fair simple variance is
>
>     SVIX^2 = (2 e^{r t} / (t F^2)) * ( OTM option strip ),
>
> reported as ``100 * SVIX``. A flat smile returns approximately
> ``100 * sigma`` (equal to VIX only to leading order; the two differ at higher
> order in vol).

### `vix_from_chain(strikes, q_prices, F, t, r)`  _function_

> CBOE variance/VIX from a discrete OTM option chain.
>
> Args:
>     strikes: increasing list of strikes.
>     q_prices: the OTM option mid-price at each strike (put below the forward,
>         call above; at the money use the average of the two).
>     F: implied forward. t: tenor in years. r: risk-free rate.
>
> Returns ``(variance, vix)`` -- the annualized fair variance and
> ``100 * sqrt(variance)``.

### `vix_from_smile(S0, t, r, vol_fn, q=0.0, n_strikes=201, width=6.0)`  _function_

> VIX-style fair index built from a smile ``vol_fn(K)``.
>
> Samples an OTM chain around the forward, prices each option at its smile vol
> with Black-Scholes, and applies :func:`vix_from_chain`. A flat smile returns
> ``VIX ~= 100 * sigma``.

## volatility

### `GarchParams(omega: float, alpha: float, beta: float) -> None`  _class_

> GarchParams(omega: float, alpha: float, beta: float)

### `VolConePoint(window: int, minimum: float, p25: float, median: float, p75: float, maximum: float, current: float) -> None`  _class_

> VolConePoint(window: int, minimum: float, p25: float, median: float, p75: float, maximum: float, current: float)

### `VolReport(close_to_close: float, parkinson: float, garman_klass: float, rogers_satchell: float, yang_zhang: float, ewma: float) -> None`  _class_

> VolReport(close_to_close: float, parkinson: float, garman_klass: float, rogers_satchell: float, yang_zhang: float, ewma: float)

### `close_to_close(closes: Sequence[float], periods_per_year: int = 252, ddof: int = 1) -> float`  _function_

> Classic close-to-close realized volatility (annualized).
>
> Uses the sample standard deviation of log returns with ``ddof`` degrees of
> freedom removed (1 = unbiased sample variance).

### `ewma_vol(closes: Sequence[float], lam: float = 0.94, periods_per_year: int = 252) -> float`  _function_

> RiskMetrics-style exponentially weighted volatility (annualized).
>
> Variance_t = lam * Variance_{t-1} + (1 - lam) * r_t^2, seeded with the
> first squared return. ``lam=0.94`` is the RiskMetrics daily default.

### `fit_garch(returns: Sequence[float], periods_per_year: int = 252)`  _function_

> Fit a GARCH(1,1) variance model to a return series by quasi-MLE.
>
> Model: ``h_t = omega + alpha * r_{t-1}^2 + beta * h_{t-1}`` with the returns
> assumed conditionally normal (Gaussian quasi-likelihood). Fitted with the
> built-in Nelder-Mead over a smooth reparametrization that keeps
> ``omega > 0``, ``alpha, beta >= 0`` and ``alpha + beta < 1`` (stationary).
>
> Returns ``GarchParams`` (per-period variance parameters).

### `garch_forecast(params: quantforge.volatility.GarchParams, last_return, last_variance, horizon=1, periods_per_year: int = 252)`  _function_

> Forecast annualized volatility ``horizon`` periods ahead under GARCH(1,1).
>
> The one-step-ahead variance is ``h_1 = omega + alpha r^2 + beta h``. Beyond
> that the expected variance mean-reverts toward the long-run level at rate
> ``persistence`` per step: ``E[h_k] = LR + persistence^{k-1} (h_1 - LR)``.
> Returns the annualized volatility for the ``horizon``-step-ahead period.

### `garch_option_price(params: quantforge.volatility.GarchParams, last_return, last_variance, S, K, r, option_type='call', horizon=None, t=None, periods_per_year: int = 252, b=None)`  _function_

> Black-Scholes price using the GARCH term volatility for the maturity.
>
> Bridges the GARCH variance forecast to an option price: the annualized term
> (average) volatility over ``horizon`` steps -- :func:`garch_term_variance` --
> is fed into the Black-Scholes formula. ``horizon`` is the number of GARCH
> steps to expiry; the option's year fraction ``t`` defaults to
> ``horizon / periods_per_year`` but may be passed explicitly (e.g. to use
> calendar rather than trading time). This lets a fitted GARCH model price
> options consistently with its own vol term structure -- capturing the vol
> mean-reversion that a single spot vol misses.

### `garch_term_variance(params: quantforge.volatility.GarchParams, last_return, last_variance, horizon, periods_per_year: int = 252)`  _function_

> Annualized GARCH term (average) volatility over the next ``horizon`` steps.
>
> An option maturing in ``horizon`` periods is priced off the *average* of the
> per-step conditional variances, not a single step. Summing the mean-reverting
> forecasts ``E[h_k] = LR + persistence^{k-1} (h_1 - LR)`` gives, for
> persistence ``p < 1``,
>
>     avg_var = LR + (h_1 - LR)/horizon * (1 - p^horizon)/(1 - p),
>
> the closed form of the geometric-series average. This is the volatility to
> feed a Black-Scholes price for that maturity. Returns the annualized term
> volatility ``sqrt(avg_var * periods_per_year)``.

### `garman_klass(opens, highs, lows, closes, periods_per_year: int = 252) -> float`  _function_

> Garman-Klass OHLC estimator (annualized).
>
> var = mean( 0.5*ln(H/L)^2 - (2ln2 - 1)*ln(C/O)^2 ). Uses the full bar; more
> efficient than Parkinson, still assumes no overnight jump or drift.

### `parkinson(highs: Sequence[float], lows: Sequence[float], periods_per_year: int = 252) -> float`  _function_

> Parkinson high-low range estimator (annualized).
>
> var = (1 / (4 ln2)) * mean( ln(H/L)^2 ). ~5x more efficient than
> close-to-close but ignores drift and overnight moves.

### `rogers_satchell(opens, highs, lows, closes, periods_per_year: int = 252) -> float`  _function_

> Rogers-Satchell OHLC estimator (annualized).
>
> var = mean( ln(H/C)ln(H/O) + ln(L/C)ln(L/O) ). Drift-independent: stays
> unbiased even when the underlying has a non-zero mean return.

### `vol_cone(closes: Sequence[float], windows: Sequence[int], periods_per_year: int = 252)`  _function_

> Realized-volatility cone: the distribution of rolling realized vol per window.
>
> For each window length, computes the annualized close-to-close realized vol
> over every rolling block of returns of that length, then reports the min,
> 25th/50th/75th percentiles, max, and the most-recent (current) value. This
> is the standard "vol cone" used to judge whether current realized vol is
> high or low versus its own history at each horizon.
>
> Args:
>     closes: the price series.
>     windows: rolling window lengths in *returns* (e.g. [5, 21, 63, 126]).
>     periods_per_year: annualization factor.
>
> Returns a list of :class:`VolConePoint`, one per window (skipping windows
> too long for the data).

### `vol_report(opens, highs, lows, closes, periods_per_year: int = 252, ewma_lambda: float = 0.94) -> quantforge.volatility.VolReport`  _function_

> Compute every estimator at once for an OHLC series.

### `yang_zhang(opens, highs, lows, closes, periods_per_year: int = 252) -> float`  _function_

> Yang-Zhang estimator (annualized): drift-independent and jump-robust.
>
> Combines overnight (close-to-open) variance, open-to-close variance, and
> the Rogers-Satchell term:
>
>     var = var_overnight + k * var_open_to_close + (1 - k) * var_RS
>     k   = 0.34 / (1.34 + (N+1)/(N-1))
>
> Requires the previous close, so bars are chained: overnight return uses
> ln(O_t / C_{t-1}).

## volcube

### `VolCube(expiries, tenors, node_params, forwards)`  _class_

> A SABR-per-node swaption vol cube with variance interpolation.

## vrp

### `realized_variance(closes: Sequence[float], periods_per_year: int = 252)`  _function_

> Annualized realized variance of log returns from a close series.

### `variance_risk_premium(closes, implied_variance, periods_per_year=252)`  _function_

> Realized-minus-implied variance risk premium and its components.
>
> Args:
>     closes: realized price history over the measurement window.
>     implied_variance: the annualized implied (variance-swap) variance for
>         the same horizon -- e.g. from a smile replication.
>     periods_per_year: sampling frequency of the closes.
>
> Returns a dict with ``realized_variance``, ``implied_variance``,
> ``vrp`` (realized - implied; usually negative), ``ratio``
> (realized / implied), and ``vol_premium`` (implied vol - realized vol, the
> usual positive number quoted in vol points).

## weather

### `cooling_degree_days(temps, base=65.0)`  _function_

> Accumulated cooling degree days ``sum_d max(T_d - base, 0)`` over the period.

### `degree_day_collar(expected_index, cap_strike, floor_strike, sigma, r, expiry, tick_value)`  _function_

> Zero-cost-style degree-day collar: long a call, short a put.
>
> Buys protection above ``cap_strike`` (a call) and finances it by selling a put
> struck at ``floor_strike``. Value is
> ``degree_day_option(call, K=cap) - degree_day_option(put, K=floor)``. When both
> strikes coincide the collar reduces to the discounted forward payoff
> ``e^{-r T} tick (expected_index - strike)`` by put-call parity.

### `degree_day_index(temps, base=65.0, kind='HDD')`  _function_

> Accumulated degree-day index of the requested ``kind`` ("HDD" or "CDD").

### `degree_day_option(expected_index, strike, sigma, r, expiry, tick_value, is_call=True, cap=None)`  _function_

> Bachelier price of an option on an accumulated degree-day index.
>
> The seasonal degree-day total is modelled as normal with mean
> ``expected_index`` and standard deviation ``sigma`` (in index points), so a
> call (protection against a high index) or put (low index) is priced by the
> Bachelier formula and scaled by ``tick_value``. With ``m = expected_index -
> strike`` for a call (``strike - expected_index`` for a put) and ``s = sigma``:
>
>     undiscounted = m Phi(m/s) + s phi(m/s)
>     price = e^{-r T} * tick_value * undiscounted
>
> An optional ``cap`` limits the maximum payoff (points), pricing the capped leg
> as a call spread: ``value(strike) - value(strike + cap)``. Put and call satisfy
> ``C - P = e^{-r T} tick (expected_index - strike)`` when uncapped.

### `degree_day_option_mc(daily_means, daily_sigma, base, strike, r, expiry, tick_value, kind='HDD', is_call=True, n_paths=20000, seed=4321)`  _function_

> Monte Carlo degree-day option over simulated daily temperatures.
>
> Simulates each day's average temperature as independent normal
> ``N(daily_means[d], daily_sigma^2)``, accumulates the HDD/CDD index over the
> period, and averages the discounted option payoff. An independent reference for
> the Bachelier :func:`degree_day_option` (which approximates the accumulated
> index as normal). Deterministic per seed.

### `degree_day_swap_payoff(index, strike, tick_value, notional_side=1.0)`  _function_

> Linear (swap) payoff on a degree-day index: ``side * tick * (index - strike)``.
>
> ``tick_value`` is the currency amount per index point; ``notional_side`` is
> ``+1`` for the long-index side (gains when the index exceeds the strike) and
> ``-1`` for the short.

### `degree_day_swap_rate(expected_index)`  _function_

> Fair fixed strike of a degree-day swap: the expected accumulated index.
>
> A degree-day swap pays ``tick * (index - strike)``; its expected value is zero
> when the strike equals the expected index, so the fair strike is
> ``expected_index`` itself.

### `expected_temperature(current_temp, seasonal_now, seasonal_future, kappa, horizon)`  _function_

> Expected temperature under a mean-reverting (OU) temperature model.
>
> Temperature reverts to its seasonal mean at speed ``kappa``; the deviation from
> the seasonal curve decays exponentially:
>
>     E[T_h] = seasonal_future + e^{-kappa h} (current_temp - seasonal_now).
>
> Equals ``current_temp`` at ``horizon = 0`` and relaxes to the future seasonal
> mean as ``horizon -> inf``.

### `heating_degree_days(temps, base=65.0)`  _function_

> Accumulated heating degree days ``sum_d max(base - T_d, 0)`` over the period.
>
> ``temps`` is the sequence of daily average temperatures. Each cold day (below
> ``base``) contributes its shortfall; warm days contribute nothing.

### `seasonal_mean_temperature(day, a, b, amplitude, phase, period=365.0)`  _function_

> Deterministic seasonal mean temperature on a given day.
>
> The Alaton-Djehiche-Stillberger seasonal trend
> ``a + b * day + amplitude * sin(2 pi (day - phase) / period)`` -- a linear
> warming/cooling trend ``b`` plus an annual sinusoid. Used as the reversion
> level of the mean-reverting :func:`expected_temperature`.

### `temperature_variance(sigma, kappa, horizon)`  _function_

> Variance of temperature under the OU model ``sigma^2 (1 - e^{-2 kappa h})/(2 kappa)``.
>
> Zero at ``horizon = 0``, rising monotonically to the stationary variance
> ``sigma^2 / (2 kappa)`` as ``horizon -> inf`` (same form as the Schwartz
> commodity model). ``sigma`` is the daily temperature volatility.

## xva

### `bcva(cpty_curve, own_curve, grid_times, epe, ene, r, cpty_recovery=0.4, own_recovery=0.4)`  _function_

> Bilateral CVA ``BCVA = CVA - DVA``.
>
> ``epe`` is the expected positive exposure profile (counterparty default risk)
> and ``ene`` the expected negative exposure profile (our default benefit).
> Returns the net adjustment to the risk-free value; positive when counterparty
> risk dominates.

### `collateralized_exposure(uncollateralized_exposure, threshold, min_transfer_amount=0.0, independent_amount=0.0)`  _function_

> Exposure remaining after a CSA collateral agreement.
>
> Under a credit-support annex the counterparty posts collateral once the
> uncollateralized exposure exceeds a ``threshold`` (plus the minimum transfer
> amount ``min_transfer_amount``), and an ``independent_amount`` of collateral is
> held unconditionally. The residual exposure is
>
>     max(min(E, threshold + MTA) - independent_amount, 0)
>
> -- exposure below the call level is uncollateralized, and above it only the
> threshold + MTA remains at risk (before default-time gap risk). An infinite
> threshold recovers the uncollateralized exposure; a zero threshold with no MTA
> leaves only the independent-amount offset.

### `collateralized_exposure_profile(uncollateralized_profile, threshold, min_transfer_amount=0.0, independent_amount=0.0)`  _function_

> Apply :func:`collateralized_exposure` across an exposure profile.

### `cva(curve, grid_times, expected_exposure, r, recovery=0.4)`  _function_

> Unilateral CVA from an expected-exposure profile and a survival curve.
>
> ``expected_exposure[i]`` is the positive expected exposure at ``grid_times[i]``
> (the representative exposure over bucket ``i``), discounted by ``DF`` (flat
> rate ``r`` or a callable ``r(t)``) and weighted by the counterparty's marginal
> default probability :func:`marginal_default_probs`. Scaled by
> ``LGD = 1 - recovery``. Non-negative, increasing in exposure and in hazard.

### `dva(own_curve, grid_times, negative_expected_exposure, r, recovery=0.4)`  _function_

> Debit valuation adjustment: the mirror of :func:`cva` on our own default.
>
> ``negative_expected_exposure[i]`` is the expected exposure of the counterparty
> to us (our negative exposure, entered as a non-negative magnitude). Weighted by
> *our* marginal default probability from ``own_curve`` and ``LGD``. A benefit to
> us, so it is subtracted from CVA in the bilateral adjustment.

### `fva(grid_times, expected_exposure, funding_spread, r, survival=None)`  _function_

> Funding valuation adjustment on an uncollateralized exposure.
>
> The cost of funding the expected positive exposure at a ``funding_spread`` over
> the risk-free rate:
>
>     FVA = funding_spread * sum_i EE(t_i) DF(t_i) [S(t_{i-1}) - ... ] approx
>         = funding_spread * sum_i EE(t_i) DF(t_i) dt_i * survival(t_i)
>
> Here it is discretized as ``funding_spread * sum_i EE_i DF_i dt_i`` optionally
> weighted by a survival probability ``survival(t)`` (both counterparties alive).
> Proportional to the spread and the exposure; zero at zero spread.

### `marginal_default_probs(curve, grid_times)`  _function_

> Marginal default probability in each grid bucket ``Q(t_{i-1}) - Q(t_i)``.
>
> ``grid_times`` are the bucket end points (strictly increasing, positive); the
> first bucket runs from 0. Returns one probability per bucket, each in
> ``[0, 1]`` and summing to ``1 - Q(t_last)`` (the total default probability by
> the horizon).

### `mva(grid_times, initial_margin, funding_spread, r, survival=None)`  _function_

> Margin valuation adjustment: funding cost of posted initial margin.
>
> Initial margin posted to a CCP or under uncleared-margin rules must be funded
> at a spread over the risk-free rate for the life of the trade. Discretized as
>
>     MVA = funding_spread * sum_i IM(t_i) DF(t_i) dt_i [* survival(t_i)],
>
> where ``initial_margin[i]`` is the IM held over bucket ``i`` (often set to a
> high-quantile :func:`swap_potential_future_exposure`). Proportional to the
> spread and the margin; zero at zero spread. Same shape as :func:`fva` but on
> the margin rather than the net exposure.

### `swap_cva(curve, notional, sigma, maturity, grid_times, r, recovery=0.4)`  _function_

> One-shot unilateral CVA of a par swap from its analytic exposure profile.
>
> Convenience wrapper: builds the :func:`swap_expected_exposure` profile and
> feeds it to :func:`cva` against the counterparty ``curve``. Equivalent to
> composing the two calls by hand.

### `swap_expected_exposure(notional, sigma, maturity, grid_times)`  _function_

> Expected positive exposure profile of a single-rate swap/forward.
>
> A par swap starts at zero value and matures at zero, with its mark-to-market
> diffusing in between. Modelling the value as a driftless Brownian motion with
> per-year volatility ``sigma`` (in value units per unit notional), the value at
> ``t`` is normal with standard deviation ``notional * sigma * sqrt(t)`` scaled by
> the remaining life ``(maturity - t)/maturity`` (linear amortization of the
> remaining risk). The expected positive exposure of a mean-zero normal is
> ``EPE(t) = std(t) / sqrt(2 pi)``. Returns one EPE per grid time; zero at
> ``t = 0`` and ``t = maturity``, humped in between.

### `swap_potential_future_exposure(notional, sigma, maturity, grid_times, quantile=0.95)`  _function_

> Potential future exposure (PFE) profile of a par swap at a high quantile.
>
> Same diffusing-then-amortizing value model as
> :func:`swap_expected_exposure`, but reports the ``quantile`` (e.g. 95th
> percentile) of the positive exposure rather than its mean. For a mean-zero
> normal value with standard deviation ``std(t)`` the upper-tail exposure
> quantile is ``std(t) * Phi^{-1}(quantile)``. Since ``Phi^{-1}(q) > 1/sqrt(2 pi)``
> for ``q`` above ~0.69, the PFE sits above the EPE at usual regulatory
> quantiles.

### `wrong_way_cva(curve, grid_times, expected_exposure, r, recovery=0.4, alpha=0.0)`  _function_

> CVA with a linear wrong-way-risk scaling of the default buckets.
>
> Wrong-way risk is the tendency of exposure to rise as the counterparty's
> credit deteriorates. This applies a simple multiplicative tilt to the marginal
> default probabilities that grows with time,
> ``weight_i = 1 + alpha * (t_i / t_last - 0.5)``, renormalized to preserve the
> total default probability. ``alpha > 0`` shifts default mass toward the later,
> higher-exposure buckets (wrong-way), raising the CVA above the independent
> ``alpha = 0`` case; ``alpha < 0`` is right-way risk. Reduces to :func:`cva`
> at ``alpha = 0``.
