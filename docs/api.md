# QuantForge API reference

Auto-generated from `quantforge` v1.76.0 by `docs/gen_api.py` — do not edit by hand.

## american

### `bjerksund_stensland(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None)`  _function_

> American option price via Bjerksund-Stensland (2002), closed form.
>
> Args mirror the rest of the engine. ``b`` is the cost of carry (defaults to
> ``r``); dividend yield q enters as b = r - q. American puts are priced via
> the exact put-call transformation P(S,K,r,b) = C(K,S,r-b,-b).

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

### `bachelier_delta(F, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>) -> float`  _function_

> dPrice/dF (in the forward). Call delta is e^{-rt} N(d).

### `bachelier_gamma(F, K, t, r, sigma) -> float`  _function_

> d2Price/dF2. Same for calls and puts.

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

### `bachelier_vega(F, K, t, r, sigma) -> float`  _function_

> dPrice/dsigma_N (per unit of normal vol). Same for calls and puts.

## binomial

### `american_price(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None, steps=500)`  _function_

> Price an American option via a CRR binomial tree.
>
> Args:
>     steps: number of time steps. Higher = more accurate, O(steps^2) work.
>     b: cost of carry (defaults to r). Dividend yield q enters as b = r - q.

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

### `gamma(S, K, t, r, sigma, b=None) -> float`  _function_

> d2Price/dS2. Identical for calls and puts.

### `greeks(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None) -> quantforge.bsm.Greeks`  _function_

> Compute price and all first/second-order Greeks in one call.

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

### `put_price(S, K, t, r, sigma, b=None) -> float`  _function_

> (no docstring)

### `rho(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None) -> float`  _function_

> dPrice/dr, per 1.0 change in rate.
>
> Assumes carry moves with the rate (the plain BSM stock case). For models
> where ``b`` is fixed independently of ``r`` (e.g. Black-76), pass ``b`` and
> interpret accordingly.

### `theta(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None) -> float`  _function_

> Calendar-time theta, dPrice/d(calendar time) per year.
>
> This is the market convention: equal to ``-dPrice/dt_expiry``, so long
> options usually show negative theta (value decays as the clock advances).
> Divide by 365 for per-calendar-day decay.

### `vega(S, K, t, r, sigma, b=None) -> float`  _function_

> dPrice/dSigma, per 1.0 change in vol (divide by 100 for per-vol-point).

## cev

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

### `noncentral_chisq_cdf(x: float, k: float, lam: float) -> float`  _function_

> Noncentral chi-square CDF at ``x`` with ``k`` dof and noncentrality ``lam``.
>
> A Poisson(lam/2)-weighted sum of central chi-square CDFs:
> ``F(x; k, lam) = sum_j pois(j; lam/2) * P((k+2j)/2, x/2)``. The summation
> starts at the Poisson mode ``j0 = floor(lam/2)`` and expands outward, so it
> stays numerically stable even when ``lam`` is large (the naive j=0 start
> underflows because ``e^{-lam/2}`` is zero to machine precision).

## chooser

### `chooser_option(S, K, t_choose, T, r, sigma, b=None) -> float`  _function_

> Price a simple chooser option (Rubinstein 1991).
>
> Args:
>     t_choose: time (years) until the call/put choice is made.
>     T: total time (years) to the underlying option's expiry (>= t_choose).
>     b: cost of carry (defaults to r).

## cir

### `cir_zero_coupon_bond(r0, t, kappa, theta, sigma)`  _function_

> CIR zero-coupon bond price P(0, t) for a unit face, given r(0)=r0.
>
> Requires ``r0 >= 0`` and positive parameters.

### `cir_zero_coupon_yield(r0, t, kappa, theta, sigma)`  _function_

> Continuously-compounded yield of the CIR zero-coupon bond to ``t``.

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

### `implied_correlation(weights: Sequence[float], vols: Sequence[float], index_vol: float) -> float`  _function_

> Common implied correlation consistent with the quoted ``index_vol``.
>
> Returns rho in principle within [-1, 1]; a value outside that band signals
> an index vol inconsistent with the member vols (arbitrage or stale quotes)
> and is returned unclamped so the caller can see it.

### `index_vol_from_correlation(weights: Sequence[float], vols: Sequence[float], rho: float) -> float`  _function_

> Index volatility implied by member weights/vols and a common correlation.

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

## displaced

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

### `gap_option(S, K_trigger, K_payoff, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None)`  _function_

> Gap option: pays off against ``K_payoff`` but is triggered by ``K_trigger``.
>
> A gap call pays ``S_T - K_payoff`` (which may be negative) whenever
> ``S_T > K_trigger``; a gap put pays ``K_payoff - S_T`` whenever
> ``S_T < K_trigger``. Setting the two strikes equal recovers the vanilla
> option. Closed form (Reiner-Rubinstein).

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

### `power_option(S, K, t, r, sigma, power, option_type=<OptionType.CALL: 'call'>, b=None)`  _function_

> Power option with payoff ``max(S_T^power - K, 0)`` (call) / ``max(K - S_T^power, 0)``.
>
> S_T^power is lognormal, so this has a closed form: an adjusted-drift,
> adjusted-vol Black-Scholes on the transformed underlying. ``power = 1``
> recovers the vanilla option.

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

## gramcharlier

### `corrado_su_call(S, K, t, r, sigma, skew=0.0, excess_kurt=0.0, b=None) -> float`  _function_

> Corrado-Su (1996) skew/kurtosis-adjusted European call price.
>
> Args:
>     skew: skewness of the (log) return distribution.
>     excess_kurt: excess kurtosis (kurtosis - 3).
>     b: cost of carry (defaults to r). skew=kurt=0 => Black-Scholes.

### `corrado_su_price(S, K, t, r, sigma, skew=0.0, excess_kurt=0.0, option_type=<OptionType.CALL: 'call'>, b=None) -> float`  _function_

> Corrado-Su price for a call or put (put via put-call parity).

### `realized_excess_kurtosis(returns: Sequence[float]) -> float`  _function_

> Sample excess kurtosis (kurtosis - 3) of a return series.

### `realized_skewness(returns: Sequence[float]) -> float`  _function_

> Sample skewness of a return series (bias-corrected denominator n).

## greeks2

### `charm(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None) -> float`  _function_

> Calendar charm = -d(delta)/d(t_expiry): the drift of delta over time.

### `color(S, K, t, r, sigma, b=None) -> float`  _function_

> Calendar color = -d(gamma)/d(t_expiry): the decay of gamma over time.

### `speed(S, K, t, r, sigma, b=None) -> float`  _function_

> d(gamma)/d(spot). Third-order in spot; same for calls and puts.

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

## merton

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

## montecarlo

### `MCResult(price: float, std_error: float, n_paths: int) -> None`  _class_

> MCResult(price: float, std_error: float, n_paths: int)

### `arithmetic_asian_mc(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None, n_steps=50, n_paths=50000, antithetic=True, control_variate=True, seed=None) -> quantforge.montecarlo.MCResult`  _function_

> Price a fixed-strike arithmetic-average-price Asian option.
>
> With ``control_variate=True`` the geometric-average Asian (known in closed
> form) is used as a control, dramatically reducing the standard error since
> the two averages are almost perfectly correlated.

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

### `double_knockout_mc(S, K, t, r, sigma, lower, upper, option_type=<OptionType.CALL: 'call'>, b=None, rebate=0.0, n_steps=100, n_paths=50000, antithetic=True, seed=None) -> quantforge.montecarlo.MCResult`  _function_

> Monte Carlo a double-knockout barrier option (a corridor).
>
> The option pays the vanilla payoff only if the spot stays strictly inside
> ``(lower, upper)`` for the whole monitored path; if either barrier is
> breached it knocks out and pays the cash ``rebate`` (at expiry, discounted).
> Also known as a double-barrier knock-out or "corridor" option.

### `european_mc(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None, n_paths=100000, antithetic=True, seed=None) -> quantforge.montecarlo.MCResult`  _function_

> Monte Carlo price of a European option (converges to the BSM value).

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

## multiasset

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

### `spread_option(S1, S2, K, t, r, sigma1, sigma2, rho, q1=0.0, q2=0.0, option_type=<OptionType.CALL: 'call'>) -> float`  _function_

> Kirk (1995) approximation for a spread option: payoff max(S1 - S2 - K, 0).
>
> Reduces to an exact Margrabe formula when K = 0. Puts follow from parity on
> the spread ``S1 - S2``.

### `worst_of_call(S1, S2, K, t, r, sigma1, sigma2, rho, q1=0.0, q2=0.0, option_type=<OptionType.CALL: 'call'>, n_paths=100000, antithetic=True, seed=None)`  _function_

> Option on the minimum of two assets: payoff max(min(S1,S2) - K, 0) (call).

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

## rates

### `CapletPeriod(forward: float, expiry: float, accrual: float, discount: float, sigma_n: float) -> None`  _class_

> CapletPeriod(forward: float, expiry: float, accrual: float, discount: float, sigma_n: float)

### `annuity(periods: Sequence[quantforge.rates.CapletPeriod]) -> float`  _function_

> Present-value annuity (level / PV01) of a swap: sum of accrual*discount.

### `cap_price(periods: Sequence[quantforge.rates.CapletPeriod], strike: float) -> float`  _function_

> Price an interest-rate cap as the sum of its caplets.

### `caplet_floorlet_parity(period: quantforge.rates.CapletPeriod, strike: float) -> float`  _function_

> Caplet - floorlet at the same strike = discounted forward-minus-strike.
>
> A put-call-parity identity used to check the pricer:
> ``caplet - floorlet = discount * accrual * (F - K)``.

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

### `floor_price(periods: Sequence[quantforge.rates.CapletPeriod], strike: float) -> float`  _function_

> Price an interest-rate floor as the sum of its floorlets.

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

## sabr

### `SABRParams(alpha: float, beta: float, rho: float, nu: float) -> None`  _class_

> SABRParams(alpha: float, beta: float, rho: float, nu: float)

### `calibrate_sabr(F, t, strikes: Sequence[float], market_vols: Sequence[float], beta: float = 0.5, weights: Sequence[float] = None, initial: quantforge.sabr.SABRParams = None, max_iter: int = 4000) -> Tuple[quantforge.sabr.SABRParams, float]`  _function_

> Fit (alpha, rho, nu) of a SABR smile to market Black vols; ``beta`` fixed.
>
> Returns ``(params, rmse)`` where rmse is the root-mean-square vol error.
> Uses a smooth constrained reparametrization so alpha > 0, nu >= 0 and
> rho in (-1, 1), optimized with the built-in Nelder-Mead.

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

## spline

### `CubicSpline(xs: Sequence[float], ys: Sequence[float])`  _class_

> Natural cubic spline through ``(xs, ys)`` with ``xs`` strictly increasing.

### `SmileSpline(strikes: Sequence[float], vols: Sequence[float])`  _class_

> Strike -> implied-vol natural cubic spline with flat extrapolation.

## strategy

### `backspread(S, K_short, K_long, t, r, sigma, kind='call', ratio=2, b=None, mult=1.0)`  _function_

> Backspread: short 1 option at K_short, long ``ratio`` at K_long.
>
> The mirror of a ratio spread -- net long options, so it profits from a large
> move (unlimited upside for a call backspread) and loses a little in the
> middle. Returns the leg :class:`Book`.

### `break_evens(book: quantforge.portfolio.Book, lo: float, hi: float, n: int = 2000) -> List[float]`  _function_

> Find terminal spots where total P&L (payoff - net premium) crosses zero.
>
> Scans ``[lo, hi]`` on a grid and refines each sign change by bisection. Net
> premium is the book's market value now (positive = we paid it).

### `butterfly(S, K_low, K_mid, K_high, t, r, sigma, kind='call', b=None, mult=1.0)`  _function_

> Long butterfly: +1 K_low, -2 K_mid, +1 K_high (equally spaced strikes).

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

### `vertical_spread(S, K_long, K_short, t, r, sigma, kind='call', b=None, mult=1.0)`  _function_

> Bull/bear vertical: long one option at K_long, short one at K_short.
>
> A call spread with K_long < K_short is a bull spread; a put spread with
> K_long > K_short is a bear spread.

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

### `svi_butterfly_arbitrage(p: quantforge.svi.SVIParams, ks=None, tol=1e-10)`  _function_

> Return the log-moneyness points where the SVI slice has butterfly arb.
>
> Scans ``ks`` (default a wide grid) and reports those where ``g(k) < -tol``.
> An empty list means the slice is butterfly-arbitrage-free on the grid.

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

### `variance_gamma_price(S, K, t, r, sigma, nu, theta, option_type=<OptionType.CALL: 'call'>, q=0.0, upper=200.0)`  _function_

> Price a European option under the Variance-Gamma model.
>
> Args:
>     sigma: Brownian volatility. nu: gamma-time variance rate (> 0).
>     theta: Brownian drift (skew; negative for an equity left skew).
>     q: continuous dividend yield.
>
> Puts follow from put-call parity. As ``nu -> 0`` the price approaches the
> Black-Scholes value.

## varswap

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
