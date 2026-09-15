# QuantForge API reference

Auto-generated from `quantforge` v5.34.0 by `docs/gen_api.py` — do not edit by hand.

## acf

### `acf(x, nlags=20)`  _function_

> Autocorrelation function up to ``nlags`` (lag 0 included, always 1).
>
> Uses the biased (divisor ``n``) autocovariance, giving a positive-semidefinite
> sequence. Returns a list of length ``nlags + 1``.

### `pacf(x, nlags=20)`  _function_

> Partial autocorrelation function via the Durbin-Levinson recursion.
>
> Returns a list of length ``nlags + 1`` with ``pacf[0] = 1`` by convention and
> ``pacf[k]`` the order-``k`` partial autocorrelation. For an AR(p) process the
> PACF is ~0 beyond lag ``p``.

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

## adaptive_mcmc

### `adaptive_metropolis(log_prob, x0, n_samples, seed=12345, burn_in=1000, init_scale=0.1, adapt_start=200, epsilon=1e-06)`  _function_

> Haario adaptive-Metropolis sampler for a target ``log_prob`` (list of floats -> logp).
>
> The proposal is Gaussian with covariance ``(2.38^2 / d) * Cov(history) + epsilon I`` once at
> least ``adapt_start`` samples have accumulated; before that it is isotropic with standard
> deviation ``init_scale``. Returns a dict with ``samples`` (after ``burn_in``) and
> ``accept_rate``.

### `slice_sample(log_prob, x0, n_samples, w=1.0, seed=12345, burn_in=0, max_steps=50)`  _function_

> Univariate slice sampler (Neal 2003) for a scalar target ``log_prob(x) -> logp``.
>
> Draws an auxiliary height under the density, steps out an interval of initial width ``w``
> (up to ``max_steps`` expansions each side), then samples uniformly from the interval,
> shrinking on rejection. No proposal scale to tune. Returns a list of ``n_samples`` scalars.

## agglomerative

### `fcluster(X, merges, n_clusters=None, distance_threshold=None)`  _function_

> Flatten a linkage into cluster labels.
>
> Provide either ``n_clusters`` (cut so that many clusters remain) or
> ``distance_threshold`` (merge only below that distance). Returns a label per
> original point, relabeled to ``0..k-1`` in order of first appearance.

### `linkage(X, method='average')`  _function_

> Agglomerative linkage over points ``X``.
>
> Returns a list of ``n-1`` merges, each ``(cluster_a, cluster_b, distance,
> size)``. Cluster ids ``0..n-1`` are the singletons; merge ``m`` creates the new
> id ``n + m``. Merge distances are non-decreasing for single/complete/average
> linkage.

## aho_corasick

### `AhoCorasick(patterns=None)`  _class_

> Aho-Corasick automaton over a set of patterns.
>
> Build once with the patterns, then call :meth:`find` on any number of texts. Each match
> is reported as ``(end_index, pattern)`` where ``end_index`` is the index of the match's
> last character in the text (so the match is ``text[end_index - len(pattern) + 1 :
> end_index + 1]``). Duplicate patterns are collapsed; empty patterns are ignored.

## airy

### `airy_ai(x)`  _function_

> Airy function of the first kind ``Ai(x)``, solving ``y'' = x y`` with ``Ai -> 0`` as ``x -> inf``.

### `airy_bi(x)`  _function_

> Airy function of the second kind ``Bi(x)`` (the growing solution).

## alias_sampler

### `AliasSampler(weights, seed=1234567)`  _class_

> O(1)-per-draw sampler for a fixed categorical distribution (Walker's alias method).
>
> Construct from a list of non-negative ``weights`` (need not sum to 1; normalized
> internally). ``sample()`` returns one index; ``sample_many(k)`` returns ``k`` indices.
> A ``seed`` fixes the deterministic random stream for reproducibility.

## alignment

### `damerau_levenshtein(a, b)`  _function_

> True (unrestricted) Damerau-Levenshtein distance between two sequences.
>
> Counts insertions, deletions, substitutions, and transpositions of adjacent elements,
> each as a single edit -- so ``"ca" -> "ac"`` is distance 1, not 2. Uses the full
> dynamic-programming table with a last-seen index (unlike the restricted OSA variant).

### `needleman_wunsch(a, b, match=1, mismatch=-1, gap=-1)`  _function_

> Needleman-Wunsch global alignment: ``(score, aligned_a, aligned_b)``.
>
> Maximizes the total score with ``match``/``mismatch`` for aligned pairs and ``gap``
> per inserted gap ``'-'``. Returns the optimal score and the two gapped strings (as
> ``str`` when the inputs are strings, else lists). Standard ``O(len(a)*len(b))`` DP with
> traceback.

### `smith_waterman(a, b, match=2, mismatch=-1, gap=-1)`  _function_

> Smith-Waterman local alignment: the best-scoring pair of substrings.
>
> Like Needleman-Wunsch but scores never go below zero (a fresh local alignment can
> start anywhere), and the traceback runs from the highest-scoring cell back to the
> first zero. Returns ``(score, aligned_a, aligned_b)`` for the best local region.
> ``score = 0`` means no positively-scoring common substring.

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

## anderson_darling

### `anderson_darling_normal(values)`  _function_

> Anderson-Darling test that ``values`` are normal (mean/sd estimated).
>
> Returns ``(a2_star, p_value)`` where ``a2_star`` is the sample-size-adjusted
> statistic. A small p-value rejects normality; the test is especially sensitive to
> heavy tails and skew. Requires at least 8 observations for a meaningful p-value.

## anderson_darling_ksample

### `anderson_darling_ksample(*samples)`  _function_

> Scholz-Stephens k-sample Anderson-Darling test.
>
> Pass two or more samples as separate sequence arguments. Returns a dict with the
> raw statistic ``a2k``, the ``standardized`` statistic ``(A2k - (k-1)) / sqrt(var)``,
> and an approximate ``p_value``. A small p-value rejects the null that all samples
> share one distribution.

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

## ar_model

### `ar_forecast(model, history, steps=1)`  _function_

> Forecast ``steps`` ahead from a fitted AR model and recent ``history``.
>
> ``history`` must hold at least ``len(coefficients)`` most-recent observations
> (oldest first). Iterates the deterministic AR recursion (innovations set to
> their zero mean), appending each forecast to drive the next.

### `fit_ar_yule_walker(x, order)`  _function_

> Fit an AR(``order``) model by the Yule-Walker / Durbin-Levinson method.
>
> Parameters
> ----------
> x : sequence of float
>     The series.
> order : int
>     AR order ``p`` (>= 1).
>
> Returns
> -------
> dict
>     ``{"coefficients": [phi_1, ..., phi_p], "intercept": c,
>     "noise_variance": sigma2, "mean": mu}``. The coefficients solve the
>     Yule-Walker equations; ``noise_variance`` is the innovation variance from
>     the final Levinson step.

## ar_select

### `ar_information_criteria(x, order)`  _function_

> AIC and BIC of an AR(``order``) fit to ``x``.
>
> Returns ``(aic, bic, loglik, noise_variance)`` using the concentrated Gaussian
> log-likelihood and ``k = order + 1`` parameters.

### `select_ar_order(x, max_order=10, criterion='aic')`  _function_

> Select the AR order minimizing AIC or BIC over ``1..max_order``.
>
> Parameters
> ----------
> x : sequence of float
>     The series.
> max_order : int
>     Largest candidate order to try.
> criterion : str
>     ``"aic"`` or ``"bic"``.
>
> Returns
> -------
> (best_order, scores) : (int, list[tuple])
>     The selected order and a list of ``(order, aic, bic)`` for every candidate.
>     BIC tends to pick an order no larger than AIC.

## ar_spectrum

### `ar_psd(coeffs, noise_variance, freqs)`  _function_

> AR power spectral density ``sigma^2 / |1 - sum phi_k e^{-i 2 pi f k}|^2``.
>
> ``freqs`` are normalized frequencies in cycles/sample (``0`` to ``0.5`` is DC to
> Nyquist). ``coeffs`` are the AR coefficients ``phi_1..phi_p`` and
> ``noise_variance`` the innovation variance. Returns the PSD at each frequency.

### `ar_spectrum(x, order, freqs, method='burg')`  _function_

> Parametric PSD of ``x`` from an AR(``order``) fit, evaluated at ``freqs``.
>
> ``method`` is ``"burg"`` (default, best for short records) or ``"yule_walker"``.
> Returns the PSD at each normalized frequency in ``freqs``.

### `burg(x, order)`  _function_

> Fit AR coefficients by Burg's method (minimizes forward+backward error).
>
> Returns ``{"coefficients": [phi_1..phi_p], "noise_variance": sigma2,
> "reflection": [k_1..k_p]}``. Burg estimates the reflection coefficients directly
> from the data (never forming autocovariances), which gives sharper, more stable
> spectra than Yule-Walker on short records. ``order`` must be ``>= 1`` and less than
> ``len(x)``.

## arch_test

### `arch_lm_test(residuals, lags=1)`  _function_

> Engle ARCH-LM test on a residual (or return) series.
>
> Regresses the squared series on ``lags`` of its own past and returns
> ``(LM, p_value)`` with ``LM = m * R^2`` (``m`` the number of regression rows)
> referenced to a chi-square with ``lags`` degrees of freedom. A small p-value
> rejects "no ARCH effect", i.e. detects volatility clustering. Subtracts the mean
> first, so it works on returns directly. Requires ``len > 2 * lags + 1``.

## assignment

### `hungarian(cost)`  _function_

> Minimum-cost assignment of ``n`` rows to ``n`` columns (Hungarian algorithm).
>
> ``cost`` is an ``n x n`` matrix. Returns ``(assignment, total_cost)`` where
> ``assignment[i]`` is the column assigned to row ``i`` and each column is used once.
> Uses the ``O(n^3)`` potentials (Kuhn-Munkres) formulation. For ``maximize``, negate
> the costs before calling.

### `knapsack_01(weights, values, capacity)`  _function_

> 0/1 knapsack: maximize total value with total weight ``<= capacity``.
>
> Each item (``weights[i]``, ``values[i]``) is taken at most once. Returns
> ``(best_value, chosen_indices)`` via the standard ``O(n * capacity)`` DP. ``capacity``
> and all weights must be non-negative integers.

## association

### `contingency_coefficient(table)`  _function_

> Pearson's contingency coefficient ``sqrt(chi2 / (chi2 + n))`` in ``[0, 1)``.
>
> Always below 1 (it cannot reach it), which is its known limitation; ``0`` under
> independence. Provided for completeness alongside Cramer's V.

### `cramers_v(table)`  _function_

> Cramer's V association strength in ``[0, 1]`` for an r x c contingency table.
>
> ``sqrt(chi2 / (n * min(r-1, c-1)))``. ``0`` for independent variables, ``1`` for a
> perfect association. The general-purpose categorical effect size.

### `phi_coefficient(table)`  _function_

> Phi coefficient for a 2x2 table: ``sqrt(chi2 / n)`` (equals Cramer's V here).
>
> Ranges ``[0, 1]`` in magnitude; the signed Pearson-correlation form of a 2x2 table.
> Raises unless the table is 2x2.

### `tschuprow_t(table)`  _function_

> Tschuprow's T association measure: ``sqrt(chi2 / (n * sqrt((r-1)(c-1))))``.
>
> Like Cramer's V but reaches ``1`` only for square tables; ``0`` under independence.

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

## bareiss

### `bareiss_determinant(matrix)`  _function_

> Exact determinant of a square integer (or rational) matrix by Bareiss elimination.
>
> Fraction-free Gaussian elimination: each step divides by the previous pivot and the
> division is always exact, so with integer input the result is an exact integer (no
> round-off, no overflow on Python big integers). Returns an ``int`` for integer input.

### `rational_inverse(A)`  _function_

> Exact inverse of a square integer/rational matrix as a matrix of ``Fraction``.
>
> Solves ``A X = I`` column by column with :func:`rational_solve`. Raises if ``A`` is
> singular. Multiplying the result by ``A`` returns the exact identity.

### `rational_solve(A, b)`  _function_

> Exact solution ``x`` of ``A x = b`` over the rationals (Gaussian elimination).
>
> ``A`` is a square integer/rational matrix and ``b`` a right-hand-side vector. Returns
> the exact solution as a list of :class:`fractions.Fraction`. Raises if ``A`` is
> singular. No round-off: an ill-conditioned system that a float solver botches comes
> out exact.

## barycentric

### `barycentric_eval(xs, ys, weights, x)`  _function_

> Evaluate the barycentric interpolant at ``x`` in ``O(n)``.
>
> ``weights`` come from :func:`barycentric_weights`. Handles the case where ``x``
> coincides with a node exactly (returns that node's value).

### `barycentric_weights(xs)`  _function_

> Barycentric weights ``w_j = 1 / prod_{k != j} (x_j - x_k)`` for nodes ``xs``.
>
> ``O(n^2)`` one-time cost; feed the result to :func:`barycentric_eval`. Raises on
> duplicate nodes.

### `chebyshev_barycentric_weights(n)`  _function_

> Closed-form barycentric weights for the ``n`` Chebyshev-Lobatto nodes.
>
> ``w_i = (-1)^i delta_i`` with the endpoints halved (``delta = 1/2`` at ``i = 0,
> n-1`` and ``1`` inside). These pair with :func:`chebyshev_nodes` and are far more
> stable than recomputing the general product. Returned in increasing-``x`` order to
> match :func:`chebyshev_nodes`.

### `chebyshev_nodes(a, b, n)`  _function_

> The ``n`` Chebyshev points of the second kind (Chebyshev-Lobatto) on ``[a, b]``.
>
> ``x_i = (a+b)/2 + (b-a)/2 * cos(i pi / (n-1))`` -- the nodes that make polynomial
> interpolation converge spectrally for smooth functions and avoid Runge
> oscillation. Returns them in increasing order.

## basket_default

### `basket_default_distribution(n, pd, rho, n_quad=200)`  _function_

> Distribution of the number of defaults in a homogeneous basket.
>
> Returns a list ``p`` of length ``n + 1`` with ``p[k] = P(exactly k defaults)``,
> integrating the conditional binomial over the common factor by midpoint
> quadrature on ``[-8, 8]``. The list sums to 1.

### `kth_to_default_probability(n, k, pd, rho, n_quad=200)`  _function_

> Probability of at least ``k`` defaults in a homogeneous basket.
>
> The trigger probability of a kth-to-default swap. Monotone: decreasing in
> ``k`` (harder to reach more defaults) and, for a fixed ``k > 1``, increasing
> in correlation (defaults cluster). ``k = 1`` is first-to-default.

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

## bayesian_regression

### `bayesian_linear_regression(X, y, alpha=1.0, beta_noise=1.0, add_intercept=True)`  _function_

> Conjugate Bayesian linear regression posterior.
>
> ``alpha`` is the prior precision (larger = stronger shrinkage toward zero),
> ``beta_noise`` the noise precision (``1 / variance``). Returns a dict with the
> posterior ``mean`` (coefficient vector, intercept first if added), the posterior
> ``covariance`` matrix, and ``std`` (per-coefficient posterior standard deviations).
> The mean equals ridge regression with ``lambda = alpha / beta_noise`` that also
> penalizes the intercept (the prior shrinks *every* coefficient toward zero).

### `bayesian_predict(model, x_row)`  _function_

> Predictive mean and variance for a new point ``x_row``.
>
> Returns ``(mean, variance)`` where the variance is the *predictive* variance
> ``1/beta_noise + x' S x`` -- observation noise plus the posterior uncertainty in
> the coefficients (so it widens where data is sparse). ``x_row`` excludes the
> intercept if the model was fit with one.

## benford

### `benford_chi_square(values)`  _function_

> Chi-square goodness-of-fit of first digits against Benford's law.
>
> ``sum (O_d - E_d)^2 / E_d`` over digits ``1..9`` with ``E_d = n P_benford(d)``,
> referenced to a chi-square with 8 degrees of freedom. Returns
> ``(statistic, p_value)``; a small p-value rejects Benford conformance.

### `benford_expected()`  _function_

> First-digit probabilities under Benford's law, ``d = 1..9``.
>
> Returns a list of 9 probabilities ``log10(1 + 1/d)`` summing to one.

### `benford_mad(values)`  _function_

> Nigrini's mean absolute deviation from Benford's first-digit law.
>
> ``MAD = (1/9) sum_d |observed_prop(d) - benford_prop(d)|``. Nigrini's rule of
> thumb: below ~0.006 is close conformance, above ~0.015 is nonconformity.
> Independent of sample size, unlike the chi-square statistic.

### `first_digit(x)`  _function_

> Leading (most significant) decimal digit of ``x``, ignoring sign and zeros.
>
> Returns an integer ``1..9``; raises for zero (no leading digit).

### `first_digit_distribution(values)`  _function_

> Observed first-digit counts and proportions for a dataset.
>
> Returns ``(counts, proportions)``, each a list of length 9 for digits ``1..9``.
> Zero values are skipped. Requires at least one non-zero value.

## berlekamp_massey

### `berlekamp_massey(sequence, mod)`  _function_

> Return the coefficients of the shortest recurrence for ``sequence`` modulo ``mod``.
>
> ``mod`` must be prime. The returned list ``c`` has length ``L`` (the recurrence order)
> and satisfies ``sequence[i] == sum_j c[j] * sequence[i-1-j] (mod mod)`` for all
> ``i >= L``. An all-zero sequence returns ``[]`` (order 0).

### `berlekamp_massey_next(sequence, mod, count=1)`  _function_

> Extend ``sequence`` by ``count`` terms using its shortest recurrence (mod ``mod``).
>
> Runs Berlekamp-Massey, then rolls the recurrence forward. Returns the list of the next
> ``count`` terms. Raises if the sequence is too short to determine any recurrence and a
> prediction is requested.

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

## bernoulli

### `bernoulli_number(n)`  _function_

> The ``n``-th Bernoulli number ``B_n`` as an exact :class:`fractions.Fraction`.
>
> Uses the convention ``B_1 = +1/2``. Computed by the recurrence
> ``sum_{k=0}^{n} C(n+1, k) B_k = 0`` and cached. ``B_n = 0`` for odd ``n > 1``.

### `bernoulli_sequence(n)`  _function_

> The Bernoulli numbers ``B_0 .. B_n`` as a list of :class:`Fraction`.

### `faulhaber(m, p)`  _function_

> Sum ``1^p + 2^p + ... + m^p`` in closed form (Faulhaber), exact.
>
> Evaluates ``(1/(p+1)) sum_{j=0}^{p} C(p+1, j) B_j m^{p+1-j}`` with ``B_1 = +1/2``.
> Returns an ``int`` (the sum is always an integer). Handles ``p = 0`` (returns ``m``).

## bessel

### `bessel_j0(x)`  _function_

> Bessel function of the first kind, order 0.

### `bessel_j1(x)`  _function_

> Bessel function of the first kind, order 1.

### `bessel_jn(n, x)`  _function_

> Bessel function of the first kind, integer order ``n`` (``n >= 0``).
>
> ``n = 0, 1`` dispatch to the direct approximations; higher orders use Miller's stable
> downward recurrence, normalized by ``J0`` (or the ``sum`` identity).

### `bessel_y0(x)`  _function_

> Bessel function of the second kind, order 0 (``x > 0``).

### `bessel_y1(x)`  _function_

> Bessel function of the second kind, order 1 (``x > 0``).

## bessel2

### `bessel_in(n, x)`  _function_

> Modified Bessel function of the first kind ``I_n(x)`` for integer ``n >= 0``.
>
> ``n = 0, 1`` use the A&S approximations in :mod:`quantforge.von_mises`; higher orders use
> downward recurrence ``I_{k-1} = I_{k+1} + (2k/x) I_k`` normalized against ``I0``.
> ``I_n(x) = I_n(-x)`` for even ``n`` and ``-`` for odd; ``I_n(0) = [n == 0]``.

### `bessel_kn(n, x)`  _function_

> Modified Bessel function of the second kind ``K_n(x)`` for integer ``n >= 0``, ``x > 0``.
>
> ``K0``/``K1`` use the A&S polynomial approximations; higher orders use the stable upward
> recurrence ``K_{k+1} = K_{k-1} + (2k/x) K_k``. Singular as ``x -> 0``.

### `spherical_bessel_j(n, x)`  _function_

> Spherical Bessel function of the first kind ``j_n(x)``.
>
> ``j_0 = sin(x)/x``, ``j_1 = sin(x)/x^2 - cos(x)/x``, and
> ``j_{k+1} = (2k+1)/x * j_k - j_{k-1}``. Uses that upward recurrence when ``n <= x`` (stable
> there) and downward Miller recurrence otherwise. ``j_n(0) = 1`` if ``n == 0`` else ``0``.

### `spherical_bessel_y(n, x)`  _function_

> Spherical Bessel function of the second kind ``y_n(x)`` (a.k.a. spherical Neumann).
>
> ``y_0 = -cos(x)/x``, ``y_1 = -cos(x)/x^2 - sin(x)/x``, then stable upward recurrence
> ``y_{k+1} = (2k+1)/x * y_k - y_{k-1}``. Singular at ``x = 0``.

## bezier

### `bernstein(n, i, t)`  _function_

> Bernstein basis polynomial ``b_{i,n}(t) = C(n, i) t^i (1-t)^(n-i)``.

### `bezier_curve(control, samples)`  _function_

> Sample the Bezier curve at ``samples`` equally spaced ``t`` in ``[0, 1]`` (inclusive).

### `bezier_derivative_control(control)`  _function_

> Control points of the curve's derivative: a degree ``n-1`` Bezier.
>
> The derivative of a degree-``n`` Bezier is degree ``n-1`` with control points
> ``n (P_{i+1} - P_i)``. Returns those points (empty for a single control point).

### `bezier_point(control, t)`  _function_

> Evaluate the Bezier curve with the given ``control`` points at parameter ``t``.
>
> Uses de Casteljau's repeated-interpolation scheme. Control points may be scalars or
> coordinate tuples (all the same length).

### `bezier_subdivide(control, t)`  _function_

> Split the curve at ``t`` into two control-point lists (left, right).
>
> The de Casteljau intermediate points give both halves exactly: the left curve uses the
> first point of each interpolation level, the right curve the last.

### `bezier_tangent(control, t)`  _function_

> The tangent vector (derivative) of the curve at ``t``.

## bfgs

### `bfgs(func, x0, tol=1e-08, max_iter=500)`  _function_

> Minimize ``func`` from ``x0`` by BFGS with a backtracking line search.
>
> ``func`` takes a length-``n`` list and returns a scalar; the gradient is computed by
> central differences. Returns a dict with ``x`` (minimizer), ``fun`` (its value),
> ``n_iter``, ``converged`` (gradient norm below ``tol``) and ``grad_norm``. Best for
> smooth objectives; use a global method first if the landscape is multimodal.

## binary_search

### `first_true(lo, hi, predicate)`  _function_

> Smallest integer ``x`` in ``[lo, hi]`` with ``predicate(x)`` true, or ``hi + 1`` if none.
>
> ``predicate`` must be monotone: once true it stays true. ``O(log(hi - lo))`` calls.

### `last_true(lo, hi, predicate)`  _function_

> Largest integer ``x`` in ``[lo, hi]`` with ``predicate(x)`` true, or ``lo - 1`` if none.
>
> ``predicate`` must be monotone-decreasing: once false it stays false.

### `ternary_search_int_max(lo, hi, f)`  _function_

> Integer in ``[lo, hi]`` maximizing a strictly unimodal ``f`` (increasing then decreasing).
>
> Narrows the range by thirds until three or fewer candidates remain, then takes the best.
> ``O(log(hi - lo))`` evaluations.

### `ternary_search_int_min(lo, hi, f)`  _function_

> Integer in ``[lo, hi]`` minimizing a strictly unimodal ``f`` (decreasing then increasing).

## binomial

### `american_price(S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None, steps=500)`  _function_

> Price an American option via a CRR binomial tree.
>
> Args:
>     steps: number of time steps. Higher = more accurate, O(steps^2) work.
>     b: cost of carry (defaults to r). Dividend yield q enters as b = r - q.

## bipartite_matching

### `maximum_bipartite_matching(adjacency)`  _function_

> Return a maximum matching as a ``{left: right}`` dict (Hopcroft-Karp).
>
> ``adjacency`` maps each left vertex to an iterable of the right vertices it can pair
> with. Left and right vertex labels live in separate namespaces (they may overlap in
> value without conflict). Only left vertices present as keys are matched.

### `maximum_matching_size(adjacency)`  _function_

> Size of a maximum bipartite matching.

### `minimum_vertex_cover(adjacency)`  _function_

> Return a minimum vertex cover ``(left_set, right_set)`` via Konig's theorem.
>
> Its total size equals the maximum matching size. Computed from the maximum matching by
> marking left vertices reachable by alternating paths from unmatched left vertices: the
> cover is the unmarked left vertices plus the marked right vertices.

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

## black_karasinski

### `bk_zero_coupon_bond(r0, kappa, theta, sigma, t, steps=50)`  _function_

> Zero-coupon bond price under Black-Karasinski via a trinomial tree.
>
> Parameters
> ----------
> r0 : float
>     Current short rate (> 0).
> kappa : float
>     Mean-reversion speed of the log rate (> 0).
> theta : float
>     Long-run mean of the *log* rate (``ln`` of the target rate level).
> sigma : float
>     Volatility of the log rate (> 0).
> t : float
>     Bond maturity (years).
> steps : int
>     Number of tree steps.
>
> Returns
> -------
> float
>     Price of a unit zero-coupon bond maturing at ``t``. Strictly in ``(0, 1]``
>     for positive rates; falls as ``r0`` or ``sigma`` rises.

## bland_altman

### `bland_altman(x, y, k=1.96)`  _function_

> Bland-Altman agreement statistics for paired measurements.
>
> Returns a dict with the ``bias`` (mean of ``x - y``), ``sd`` of the differences,
> the ``lower`` and ``upper`` limits of agreement (``bias +/- k * sd``), and the
> per-point ``means`` and ``diffs`` for plotting. ``k`` defaults to 1.96 (95% limits).

### `concordance_correlation(x, y)`  _function_

> Lin's concordance correlation coefficient (CCC) of paired measurements.
>
> ``CCC = 2 s_xy / (s_x^2 + s_y^2 + (mx - my)^2)`` -- Pearson correlation penalized
> for any departure from the line of identity ``y = x``. Lies in ``[-1, 1]``: 1 only
> when the points fall exactly on ``y = x``, and it drops below the Pearson value
> whenever there is a systematic offset or scale difference.

## bluestein

### `dft_any(x)`  _function_

> Forward DFT of a sequence of *any* length ``n`` (Bluestein), returning complex output.
>
> ``X_k = sum_n x_n exp(-2 pi i k n / N)``. Matches :func:`quantforge.fft.fft` on
> power-of-two lengths and a direct DFT everywhere, in ``O(n log n)``.

### `idft_any(x)`  _function_

> Inverse DFT of any length (Bluestein), with the ``1/N`` scaling.
>
> ``idft(dft(x)) == x`` up to floating error.

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

## brinson

### `allocation_effect(portfolio_weights, benchmark_weights, benchmark_returns)`  _function_

> Per-segment allocation effect ``(w_p - w_b) * r_b``.
>
> Positive when the portfolio overweights segments that outperformed the overall
> benchmark (using the segment benchmark return). Returns one value per segment.

### `brinson_attribution(portfolio_weights, benchmark_weights, portfolio_returns, benchmark_returns)`  _function_

> Full Brinson attribution: allocation, selection, interaction, and totals.
>
> Returns a dict with per-segment ``allocation``, ``selection``, ``interaction``
> lists, their totals, the ``active_return`` (portfolio minus benchmark total
> return), and ``total_effect`` (allocation + selection + interaction totals).
> The total effect equals the active return by construction.

### `carino_factor(portfolio_return, benchmark_return)`  _function_

> Cariño (1999) smoothing factor linking arithmetic effects across periods.
>
> ``k = ln(1 + r_p) - ln(1 + r_b)) / (r_p - r_b)`` when the returns differ, else
> ``1 / (1 + r_p)``. Scaling each period's arithmetic effects by ``k_t`` and
> dividing by the total-period factor makes the smoothed effects compound
> exactly to the geometric active return -- resolving the residual that plain
> arithmetic summation leaves across multiple periods.

### `carino_linked_effects(period_effects, portfolio_returns_by_period, benchmark_returns_by_period)`  _function_

> Cariño-smoothed multi-period effect totals that link geometrically.
>
> ``period_effects`` is a list of per-period arithmetic effect totals (e.g. the
> allocation totals from :func:`brinson_attribution` each period). Each is scaled
> by its Cariño factor and divided by the total-period factor
> ``k = (ln(1 + R_p) - ln(1 + R_b)) / (R_p - R_b)`` on the compounded returns, so
> the smoothed effects across periods sum to the geometrically-linked active
> return. Returns the smoothed per-period effect list.

### `interaction_effect(portfolio_weights, benchmark_weights, portfolio_returns, benchmark_returns)`  _function_

> Per-segment interaction effect ``(w_p - w_b) * (r_p - r_b)``.
>
> The cross term -- the combined effect of active weighting and active selection
> in the same segment.

### `linked_active_return(portfolio_returns_by_period, benchmark_returns_by_period)`  _function_

> Geometrically-linked active return over multiple periods.
>
> Compounds each side's total return across periods and returns the difference
> of the geometric returns:
>
>     (prod(1 + r_p,t) - 1) - (prod(1 + r_b,t) - 1).
>
> The quantity multi-period Brinson effects must sum to under Cariño linking.

### `selection_effect(benchmark_weights, portfolio_returns, benchmark_returns)`  _function_

> Per-segment selection effect ``w_b * (r_p - r_b)``.
>
> Positive when the portfolio's holdings within a segment beat that segment's
> benchmark, weighted at the benchmark weight.

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

## bspline

### `bspline_basis(i, p, knots, t)`  _function_

> Cox-de Boor basis function ``N_{i,p}(t)`` for knot vector ``knots``.
>
> ``i`` is the basis index, ``p`` the degree. Uses the standard recursion with the
> ``0/0 = 0`` convention for repeated knots.

### `bspline_curve(control, degree, samples, knots=None)`  _function_

> Sample the B-spline at ``samples`` equally spaced parameters over its knot range.

### `bspline_point(control, degree, t, knots=None)`  _function_

> Evaluate the B-spline curve at parameter ``t`` (basis-weighted control points).
>
> ``degree`` is the polynomial degree ``p``; ``knots`` defaults to a clamped open-uniform
> vector (so the curve interpolates the endpoints). ``t`` runs over the knot range
> ``[knots[p], knots[-p-1]]`` (``[0, 1]`` for the default knots).

### `open_uniform_knots(n_control, degree)`  _function_

> Clamped (open-uniform) knot vector for ``n_control`` points of the given ``degree``.
>
> Length ``n_control + degree + 1``: the first and last ``degree + 1`` knots are repeated
> (0 and 1), the interior knots equally spaced -- so the curve interpolates its endpoints.

## bvp

### `shooting_bvp(f, a, b, alpha, beta, s_lo, s_hi, tol=1e-08)`  _function_

> Solve ``y'' = f(t, y, y')`` with ``y(a)=alpha``, ``y(b)=beta`` by shooting.
>
> ``f(t, y, yp)`` returns ``y''``. ``s_lo, s_hi`` bracket the unknown initial slope
> ``y'(a)`` (the terminal residual must change sign across them). Returns a dict with
> the found initial ``slope``, and the solution ``ts`` / ``ys`` (state = ``[y, y']``)
> from the accepted RK45 steps.

## bwt

### `bwt_inverse(transformed, primary_index)`  _function_

> Invert the Burrows-Wheeler transform back to the original sequence.
>
> ``transformed`` is the last column and ``primary_index`` the row of the original
> string (both from :func:`bwt_transform`). Reconstructs via the standard LF-mapping.
> Returns the same type as ``transformed``.

### `bwt_transform(data)`  _function_

> Burrows-Wheeler transform of a sequence: ``(transformed, primary_index)``.
>
> Builds the sorted matrix of all rotations and returns the last column plus the index
> of the original string among the sorted rotations (needed to invert). ``data`` is a
> string or list; the transformed output has the same type. Empty input returns
> ``("" or [], 0)``.

### `move_to_front_decode(codes, alphabet)`  _function_

> Invert move-to-front coding given the ``codes`` and the initial ``alphabet``.
>
> Returns the reconstructed sequence: a ``str`` if the alphabet symbols are all
> single characters, else a list.

### `move_to_front_encode(data, alphabet=None)`  _function_

> Move-to-front encode a sequence into a list of integer ranks.
>
> Each symbol is replaced by its current index in a running alphabet list, then moved
> to the front. Clustered inputs (like BWT output) produce many zeros. ``alphabet`` is
> the initial ordered symbol list; by default the sorted set of symbols in ``data``.
> Returns ``(codes, alphabet)`` -- the alphabet is needed to decode.

## calendar_math

### `add_days(year, month, day, n)`  _function_

> Return the date ``n`` days after (or before, if ``n < 0``) the given date.

### `day_of_week(year, month, day)`  _function_

> Weekday index, 0 = Monday .. 6 = Sunday.

### `day_of_week_name(year, month, day)`  _function_

> Weekday name, e.g. ``"Wednesday"``.

### `day_of_year(year, month, day)`  _function_

> Ordinal day within the year (Jan 1 = 1).

### `days_between(date1, date2)`  _function_

> Signed day count ``date2 - date1`` (positive when ``date2`` is later).

### `days_in_month(year, month)`  _function_

> Number of days in ``month`` of ``year`` (1-12).

### `easter_date(year)`  _function_

> Gregorian Easter Sunday of ``year`` as ``(year, month, day)`` (Anonymous Computus).

### `is_leap_year(year)`  _function_

> True if ``year`` is a Gregorian leap year (divisible by 4, not 100 unless 400).

### `jdn_to_date(jdn)`  _function_

> Inverse of :func:`julian_day_number`: the ``(year, month, day)`` for a JDN.

### `julian_day_number(year, month, day)`  _function_

> Julian day number of a proleptic-Gregorian date (Fliegel-Van Flandern).
>
> A monotone integer count of days, so ``jdn(b) - jdn(a)`` is the day difference.

## calibration

### `brier_decomposition(forecasts, outcomes, n_bins=None)`  _function_

> Murphy decomposition of the Brier score into reliability/resolution/uncertainty.
>
> Returns a dict with ``reliability``, ``resolution``, ``uncertainty``, ``brier``
> (the reconstructed ``reliability - resolution + uncertainty``), and ``base_rate``.
> With ``n_bins=None`` forecasts are grouped by identical value and the identity is
> exact; with an integer ``n_bins`` they are put into equal-width bins on ``[0, 1]``
> and the reconstruction is approximate.

### `expected_calibration_error(forecasts, outcomes, n_bins=10)`  _function_

> Expected calibration error: count-weighted mean ``|forecast - observed|`` per bin.
>
> Zero for a perfectly calibrated forecaster; the standard scalar summary of a
> reliability diagram's departure from the diagonal.

### `reliability_curve(forecasts, outcomes, n_bins=10)`  _function_

> Calibration diagram: per-bin mean forecast vs observed frequency.
>
> Returns a list of ``(mean_forecast, observed_frequency, count)`` tuples, one per
> non-empty equal-width bin. A perfectly calibrated forecaster lies on the diagonal
> ``observed == forecast``.

## callable_bond

### `call_option_value(face, coupon_rate, maturity, r0, sigma, call_price, freq=1, p=0.5)`  _function_

> Value of the embedded call to the issuer: straight minus callable price.
>
> ``straight_bond - callable_bond`` -- the value the call option strips from the
> bondholder (non-negative). The callable bond trades cheaper by this amount.

### `callable_bond_price(face, coupon_rate, maturity, r0, sigma, freq=1, call_price=None, put_price=None, p=0.5)`  _function_

> Price a (possibly callable/puttable) bond on a binomial short-rate tree.
>
> ``freq`` coupons per year over ``maturity`` years; the lattice has one step per
> coupon. ``call_price`` caps the value at each pre-maturity node (issuer's call);
> ``put_price`` floors it (holder's put); either may be ``None``. With both
> ``None`` this is the straight bond on the tree. Risk-neutral up-probability
> ``p``. Returns the time-zero price.

### `callable_bond_price_with_spread(face, coupon_rate, maturity, r0, sigma, spread, freq=1, call_price=None, put_price=None, p=0.5)`  _function_

> Callable-bond price with a constant spread added to every tree rate.
>
> Shifts the whole short-rate lattice up by ``spread`` before backward
> induction, so a positive spread discounts harder and lowers the price. The
> building block for the option-adjusted spread solve.

### `option_adjusted_spread(market_price, face, coupon_rate, maturity, r0, sigma, freq=1, call_price=None, put_price=None, p=0.5, tol=1e-08, max_iter=100)`  _function_

> Option-adjusted spread: constant rate spread repricing the bond to market.
>
> Bisection on the :func:`callable_bond_price_with_spread` (monotone decreasing
> in the spread) to hit ``market_price``. Positive when the market price is below
> the zero-spread model price. Strips out the embedded option so the spread
> reflects credit/liquidity risk on a like-for-like basis.

### `straight_bond_tree_price(face, coupon_rate, maturity, r0, sigma, freq=1, p=0.5)`  _function_

> Straight (option-free) bond on the same tree -- the no-optionality baseline.

## capital_budgeting

### `irr(cashflows, tol=1e-10, max_iter=200)`  _function_

> Internal rate of return: the rate at which :func:`npv` is zero.
>
> Bisection on ``[-0.999, 10]`` (NPV is monotone decreasing in the rate for a
> conventional outlay-then-inflows project). Requires a sign change in the
> cashflows. Recovers the discount rate that generated a set of flows.

### `mirr(cashflows, finance_rate, reinvest_rate)`  _function_

> Modified internal rate of return.
>
> Compounds positive cashflows forward at ``reinvest_rate`` to the terminal date,
> discounts negative cashflows back at ``finance_rate`` to time zero, then
> ``MIRR = (FV_pos / -PV_neg)^{1/n} - 1``. Avoids the multiple-IRR problem and
> uses realistic reinvestment.

### `npv(rate, cashflows)`  _function_

> Net present value at a per-period discount ``rate``.
>
> ``sum_t CF_t / (1 + rate)^t`` with ``t = 0`` the first entry. Falls as the
> discount rate rises for a conventional project.

### `payback_period(cashflows)`  _function_

> Payback period: fractional periods to recover the initial outlay.
>
> Accumulates undiscounted cashflows until the running total turns non-negative,
> interpolating within the crossing period. Returns ``inf`` if never recovered.

### `profitability_index(rate, cashflows)`  _function_

> Profitability index: PV of inflows over the initial outlay.
>
> ``PV(CF_1..) / -CF_0`` -- above one for a value-adding project, exactly one at
> the :func:`irr`. Requires a negative time-0 outlay.

## carlson

### `carlson_rc(x, y, tol=1e-12)`  _function_

> Carlson's degenerate ``R_C(x, y) = R_F(x, y, y)``.
>
> ``R_C(x, x) = 1/sqrt(x)``; equals ``arctan``/``arctanh``-type elementary functions.

### `carlson_rd(x, y, z, tol=1e-12)`  _function_

> Carlson's ``R_D(x, y, z) = R_J(x, y, z, z)``, symmetric in ``x, y`` only.
>
> ``= 3/2 integral_0^inf dt / [(t+z) sqrt((t+x)(t+y)(t+z))]``. ``z`` must be positive; ``x, y``
> non-negative with at most one zero. Homogeneous of degree ``-3/2``.

### `carlson_rf(x, y, z, tol=1e-12)`  _function_

> Carlson's ``R_F(x, y, z) = 1/2 integral_0^inf dt / sqrt((t+x)(t+y)(t+z))``.
>
> Symmetric and homogeneous of degree ``-1/2``. At most one argument may be zero; all must be
> non-negative. ``R_F(x, x, x) = 1/sqrt(x)``.

### `carlson_rj(x, y, z, p, tol=1e-12)`  _function_

> Carlson's ``R_J(x, y, z, p) = 3/2 integral_0^inf dt / [(t+p) sqrt((t+x)(t+y)(t+z))]``.
>
> Symmetric in ``x, y, z``; all non-negative with at most one zero, and ``p != 0``. This
> implementation covers ``p > 0``.

### `elliptic_e_incomplete(phi, m, tol=1e-12)`  _function_

> Incomplete elliptic integral of the second kind ``E(phi | m)`` via ``R_F`` and ``R_D``.
>
> ``E(phi | m) = sin(phi) R_F(c, 1-m s^2, 1) - (m/3) sin^3(phi) R_D(c, 1-m s^2, 1)``.
> Reduces to the complete ``E(m)`` at ``phi = pi/2``.

### `elliptic_f(phi, m, tol=1e-12)`  _function_

> Incomplete elliptic integral of the first kind ``F(phi | m)`` via ``R_F``.
>
> ``F(phi | m) = integral_0^phi dtheta / sqrt(1 - m sin^2 theta) = sin(phi) R_F(c, 1 - m s^2, 1)``
> with ``s = sin phi``, ``c = cos^2 phi``. ``m = k^2``; reduces to ``F(pi/2, m) = K(m)``.

### `elliptic_pi(n, phi, m, tol=1e-12)`  _function_

> Incomplete elliptic integral of the third kind ``Pi(n; phi | m)`` via ``R_F`` and ``R_J``.
>
> ``Pi(n; phi | m) = integral_0^phi dtheta / [(1 - n sin^2 theta) sqrt(1 - m sin^2 theta)]``,
> ``= sin(phi) R_F(c, 1-m s^2, 1) + (n/3) sin^3(phi) R_J(c, 1-m s^2, 1, 1 - n s^2)``.
> Implemented for ``n < 1`` (so the ``R_J`` fourth argument stays positive).

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

## carry_rolldown

### `carry_return(coupon_rate, yield_now, horizon, financing_rate=0.0)`  _function_

> Carry over a horizon: coupon income plus financing, per unit face.
>
> ``(coupon_rate - financing_rate) * horizon`` -- the running yield earned net
> of the cost of funding the position, holding prices fixed. Positive when the
> coupon exceeds the financing rate.

### `rolldown_return(cashflows, curve, horizon)`  _function_

> Roll-down return: the price gain purely from the yield rolling down the curve.
>
> Isolates the yield-change effect from the time-value growth. Values the
> surviving cashflows (those maturing after the horizon), each at its *shortened*
> maturity ``t - horizon``, under two curves: the rolled yield ``y(t - horizon)``
> versus the unchanged-maturity yield ``y(t)``. The fractional difference is the
> roll-down -- zero on a flat curve (the yield does not change as the bond rolls)
> and positive on an upward-sloping curve (the bond rolls to a lower yield).

### `total_carry_rolldown(cashflows, curve, coupon_rate, horizon, financing_rate=0.0)`  _function_

> Total expected return = carry + roll-down over the horizon.
>
> Sums :func:`carry_return` and :func:`rolldown_return`. The expected holding-
> period return if the curve is unchanged; the standard relative-value carry-
> and-roll number.

## catmull_rom

### `catmull_rom_curve(points, samples_per_segment=16, alpha=0.5)`  _function_

> Sample the whole spline. Returns a flat list of points along all segments.
>
> Each segment contributes ``samples_per_segment`` points; the shared endpoints appear
> once (the start of each segment except the first is skipped). The curve passes through
> every input point.

### `catmull_rom_point(points, seg, t, alpha=0.5)`  _function_

> Evaluate segment ``seg`` (between ``points[seg]`` and ``points[seg+1]``) at ``t`` in [0,1].
>
> Needs at least two points; endpoints are handled by duplicating the boundary point as the
> phantom neighbour. ``alpha``: 0 uniform, 0.5 centripetal, 1 chordal.

## cdar

### `average_drawdown(returns)`  _function_

> Mean fractional drawdown over the return path (non-negative).

### `conditional_drawdown_at_risk(returns, confidence=0.95)`  _function_

> Conditional drawdown-at-risk (CDaR): mean of the worst ``1 - confidence`` drawdowns.
>
> The coherent drawdown analogue of expected shortfall (Chekhlov-Uryasev-
> Zabarankin). Averages the deepest tail of the drawdown distribution beyond the
> :func:`drawdown_at_risk` threshold, so ``CDaR >= DaR``. Falls back to the single
> worst drawdown when the tail holds one observation.

### `drawdown_at_risk(returns, confidence=0.95)`  _function_

> Drawdown-at-risk: the ``confidence``-quantile of the drawdown distribution.
>
> The drawdown depth that is exceeded only ``1 - confidence`` of the time. A larger
> confidence gives a deeper (more conservative) threshold. Uses the upper-tail
> order statistic of the drawdown series.

## cepstrum

### `fundamental_quefrency(x, min_quefrency=1, max_quefrency=None)`  _function_

> Quefrency (in samples) of the largest real-cepstrum peak in a search band.
>
> For a voiced/pitched or echoed signal this is the fundamental period (or echo
> delay) in samples: divide the sample rate by it to get the pitch in Hz. The search
> ignores the low-quefrency region below ``min_quefrency`` (the spectral envelope) and
> is capped at ``max_quefrency`` (default: half the cepstrum length, the useful range
> of a real cepstrum). Returns ``(quefrency, peak_value)``.

### `power_cepstrum(x)`  _function_

> Power cepstrum of ``x``: ``|ifft(log|fft(x)|^2)|^2``.
>
> Squares the log-magnitude spectrum before inverting, then takes the squared
> magnitude -- emphasizing quefrency peaks. Same length convention as
> :func:`real_cepstrum`.

### `real_cepstrum(x)`  _function_

> Real cepstrum of ``x``: ``ifft(log|fft(x)|)``, real part.
>
> The inverse FFT of the log-magnitude spectrum. Returned length equals the padded
> power-of-two length used internally (``>= len(x)``). The independent variable is
> *quefrency* -- an index in samples, read like a lag. A tiny floor is added inside
> the log to keep spectral nulls finite.

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

## chain_ladder

### `bornhuetter_ferguson(triangle, apriori_ultimates)`  _function_

> Bornhuetter-Ferguson reserves blending development with a-priori ultimates.
>
> For each accident year the BF reserve is
> ``apriori_ultimate * (1 - pct_developed)``, where ``pct_developed`` comes from
> the chain-ladder development pattern. The BF ultimate is the latest paid plus
> that reserve -- an a-priori-anchored estimate that is robust for green
> (little-developed) years where chain-ladder is volatile.
>
> Parameters
> ----------
> triangle : list[list[float]]
>     Cumulative-claims triangle (as in :func:`chain_ladder`).
> apriori_ultimates : sequence of float
>     A-priori ultimate loss per accident year (e.g. premium x expected loss
>     ratio).
>
> Returns
> -------
> dict
>     ``pattern`` (% developed by age), ``reserve`` and ``ultimate`` per year,
>     and ``total_reserve``.

### `cape_cod(triangle, premiums)`  _function_

> Cape Cod (Stanard-Buhlmann) reserving.
>
> Like Bornhuetter-Ferguson but the a-priori loss ratio is estimated from the
> data rather than assumed: the expected loss ratio is
>
>     ELR = sum_i latest_i / sum_i (premium_i * pct_developed_i),
>
> the total observed losses over the total "used-up" premium (premium weighted
> by how developed each year is). Each year's a-priori ultimate is then
> ``premium_i * ELR`` and its reserve ``apriori * (1 - pct_developed_i)``.
>
> Parameters
> ----------
> triangle : list[list[float]]
>     Cumulative-claims triangle (as in :func:`chain_ladder`).
> premiums : sequence of float
>     Earned premium per accident year.
>
> Returns
> -------
> dict
>     ``elr``, ``pattern``, ``reserve`` / ``ultimate`` per year, and
>     ``total_reserve``.

### `chain_ladder(triangle)`  _function_

> Project a claims triangle to ultimate losses and reserves.
>
> Returns a dict with ``factors`` (age-to-age), ``ultimate`` (per accident year),
> ``reserve`` (IBNR per year = ultimate - latest observed), and
> ``total_reserve``. A fully-developed row has zero reserve.

### `chain_ladder_with_tail(triangle, tail_factor)`  _function_

> Chain-ladder projection with an extra tail development factor.
>
> Multiplies each accident year's ultimate by ``tail_factor`` after the usual
> chain-ladder projection, capturing development beyond the triangle. A
> ``tail_factor`` of 1.0 reproduces :func:`chain_ladder`.

### `cumulative_to_incremental(triangle)`  _function_

> Convert a cumulative-claims triangle to incremental (successive differences).

### `development_factors(triangle)`  _function_

> Volume-weighted age-to-age development factors from a cumulative triangle.
>
> ``triangle[i]`` is the observed cumulative claims for accident year ``i`` at
> development ages ``0 .. len(triangle[i]) - 1``. Returns ``n - 1`` factors
> ``f_0 .. f_{n-2}`` linking successive development ages.

### `development_pattern(factors)`  _function_

> Cumulative development pattern (% reported) from age-to-age factors.
>
> Returns ``pct[j]`` = fraction of ultimate developed by age ``j``, computed as
> the reciprocal of the cumulative product of the remaining factors. The final
> age is fully developed (``1.0``).

### `exponential_tail_factor(factors, n_extrapolate=100)`  _function_

> Extrapolate a tail development factor by exponential decay of ``f - 1``.
>
> Fits ``ln(f_j - 1) = a + b j`` to the observed age-to-age factors with
> ``f_j > 1`` (the excess-over-one decays geometrically), projects the excess
> forward ``n_extrapolate`` ages, and returns the product ``prod (1 + excess_k)``
> as a single tail factor applied beyond the last observed age. Returns 1.0 when
> no factor exceeds 1 (fully developed). Requires ``b < 0`` (a decaying tail).

### `incremental_to_cumulative(triangle)`  _function_

> Convert an incremental-claims triangle to cumulative.
>
> Each row's cumulative entry is the running sum of its incremental entries.
> Ragged rows (shorter for recent accident years) are preserved.

### `paid_to_date(cumulative_triangle)`  _function_

> Latest (diagonal) paid amount per accident year of a cumulative triangle.

## chatterjee

### `blomqvist_beta(x, y)`  _function_

> Blomqvist's beta (medial correlation).
>
> ``beta = (n_concordant - n_discordant) / n_used`` where a point is concordant if
> it sits in the same direction from both medians (both above or both below) and
> discordant otherwise; points exactly on a median are dropped (``n_used`` counts
> only the points kept). Lies in ``[-1, 1]``:
> +1 comonotone, -1 countermonotone, 0 under independence.

### `chatterjee_xi(x, y)`  _function_

> Chatterjee's xi rank correlation of ``Y`` on ``X`` (asymmetric).
>
> Sorts by ``x`` (ties in ``x`` broken by stable order), then with
> ``r_i`` = #{j : y_j <= y_(i)} and ``l_i`` = #{j : y_j >= y_(i)},
>
>     xi = 1 - n * sum_i |r_{i+1} - r_i| / (2 * sum_i l_i (n - l_i)).
>
> Zero under independence, approaching 1 when ``Y`` is a noiseless function of
> ``X``. Handles ties in ``y`` via the general (Azadkia-Chatterjee) form.

## chebyshev

### `chebyshev_derivative(coeffs, a, b)`  _function_

> Chebyshev coefficients of the derivative of a series on ``[a, b]``.
>
> Applies the standard Chebyshev differentiation recurrence and the chain-rule
> factor ``2 / (b - a)``. The returned coefficients evaluate (via
> :func:`chebyshev_eval`) to ``f'`` on the same interval.

### `chebyshev_eval(coeffs, a, b, x)`  _function_

> Evaluate a Chebyshev series at ``x`` by the Clenshaw recurrence.
>
> ``coeffs`` are from :func:`chebyshev_fit`; the ``c_0`` term is taken at half
> weight (the standard convention). Stable across ``[a, b]``.

### `chebyshev_fit(f, a, b, degree)`  _function_

> Chebyshev coefficients of ``f`` on ``[a, b]`` up to ``degree``.
>
> Samples ``f`` at the ``degree + 1`` Chebyshev extrema and takes the discrete
> cosine transform. Returns the coefficient list ``[c_0, ..., c_degree]`` (the
> ``c_0`` term is the mean level, used at half weight by :func:`chebyshev_eval`).
> Exact for polynomials of degree ``<= degree``.

## checksums

### `adler32(data)`  _function_

> Adler-32 checksum (as in zlib), matching :func:`zlib.adler32`.
>
> Two running sums modulo 65521 combined into a 32-bit value. Cheaper than CRC-32 but
> weaker on short inputs. Returns a 32-bit unsigned integer.

### `crc32(data)`  _function_

> CRC-32 checksum (IEEE 802.3 reflected polynomial), matching :func:`zlib.crc32`.
>
> Table-driven byte-at-a-time computation over the bytes of ``data`` (a ``str`` is
> UTF-8 encoded). Returns a 32-bit unsigned integer.

### `fnv1a_32(data)`  _function_

> 32-bit FNV-1a hash: XOR then multiply per byte.
>
> A fast, well-dispersed non-cryptographic hash for tables and bloom filters. Returns a
> 32-bit unsigned integer; ``fnv1a_32("")`` is the FNV offset basis ``2166136261``.

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

## circle

### `circle_circle_intersection(c1, c2)`  _function_

> Intersection points of two circles ``(cx, cy, r)``.
>
> Returns a list of 0, 1 (tangent), or 2 ``(x, y)`` points. Coincident circles raise
> (infinitely many intersections).

### `circle_from_3points(a, b, c)`  _function_

> Circle ``(cx, cy, r)`` through three non-collinear points.
>
> Raises ``ValueError`` if the points are collinear (no finite circle).

### `circle_line_intersection(circle, a, b)`  _function_

> Intersection points of a circle with the *infinite* line through ``a`` and ``b``.
>
> Returns a list of 0, 1 (tangent), or 2 ``(x, y)`` points. ``a`` and ``b`` must be
> distinct.

### `point_in_circle(p, circle)`  _function_

> Whether point ``p`` lies inside or on the circle ``(cx, cy, r)``.

## circular_stats

### `circular_mean(angles)`  _function_

> Mean direction of a set of angles (radians), in ``(-pi, pi]``.
>
> Averages the unit vectors ``(cos, sin)`` and takes the atan2 of the result, so it
> wraps correctly across the ``2*pi`` boundary. Undefined (raises) when the vectors
> cancel to the origin (no mean direction).

### `circular_std(angles)`  _function_

> Circular standard deviation ``sqrt(-2 ln R)`` (radians).
>
> Grows without bound as the angles spread (``R -> 0``); ``0`` when all identical.

### `circular_variance(angles)`  _function_

> Circular variance ``1 - R`` in ``[0, 1]`` (0 = concentrated, 1 = dispersed).

### `rayleigh_test(angles)`  _function_

> Rayleigh test for a uniform circular distribution: returns ``(R, p_value)``.
>
> Tests the null hypothesis that the angles are uniformly spread around the circle
> against the alternative of a single preferred direction. A small ``p_value`` rejects
> uniformity. Uses the standard ``Z = n R^2`` statistic with the Zar small-sample
> correction. Needs at least two angles.

### `resultant_length(angles)`  _function_

> Mean resultant length ``R`` in ``[0, 1]``: how concentrated the angles are.
>
> ``R = |sum e^{i theta}| / n``. ``1`` means all angles identical; ``0`` means they are
> spread so the unit vectors cancel. The basis for circular variance and the Rayleigh
> test.

## classification_metrics

### `brier_score(y_true, y_score)`  _function_

> Mean squared error of the predicted probabilities. 0 = perfect.

### `confusion_matrix(y_true, y_score, threshold=0.5)`  _function_

> Confusion counts ``(tp, fp, fn, tn)`` at a probability ``threshold``.

### `log_loss(y_true, y_score, eps=1e-15)`  _function_

> Mean binary cross-entropy. 0 for a perfect confident classifier.

### `precision_recall_f1(y_true, y_score, threshold=0.5)`  _function_

> Precision, recall, and F1 at a threshold, as a ``(p, r, f1)`` tuple.
>
> Precision and recall are 0 when their denominators vanish (no predicted or no
> actual positives).

### `roc_auc(y_true, y_score)`  _function_

> Area under the ROC curve via the Mann-Whitney rank statistic.
>
> Equals the probability that a randomly chosen positive scores higher than a
> randomly chosen negative (ties = 0.5). Raises if the labels are all one class
> (AUC undefined). 1 = perfect ranking, 0.5 = random.

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

## cochran_mcnemar

### `cochran_q_test(blocks)`  _function_

> Cochran's Q test for ``k`` binary treatments over ``b`` blocks.
>
> ``blocks`` is a sequence of rows, each a length-``k`` sequence of 0/1 outcomes for
> one block across the treatments. Returns a dict with the ``statistic`` Q, ``df``
> (``k - 1``) and the chi-square ``p_value``. Blocks whose outcomes are all-0 or
> all-1 contribute nothing (as in the standard formulation).

### `mcnemar_test(table=None, b=None, c=None)`  _function_

> McNemar's test for a paired 2x2 table.
>
> Provide either the 2x2 ``table`` ``[[a, b], [c, d]]`` (a/d concordant, b/c
> discordant) or the two discordant counts ``b`` and ``c`` directly. Returns a dict
> with the discordant counts, the exact two-sided binomial ``p_value``, and the
> continuity-corrected chi-square statistic ``chi2_cc`` with its ``p_value_chi2``.

## cohen_kappa

### `cohen_kappa(rater_a, rater_b)`  _function_

> Cohen's kappa for two raters over nominal categories.
>
> ``rater_a`` and ``rater_b`` are equal-length label sequences. Returns the
> chance-corrected agreement in ``[-1, 1]``.

### `fleiss_kappa(table)`  _function_

> Fleiss' kappa for ``m`` raters per subject over fixed categories.
>
> ``table`` is a list of rows, one per subject, each giving the count of raters
> assigning each category (every row sums to the same ``m``). Returns the
> chance-corrected multi-rater agreement.

### `weighted_kappa(rater_a, rater_b, weights='linear')`  _function_

> Weighted kappa for two raters over ordinal categories.
>
> Disagreements are penalized by category distance: ``weights="linear"`` uses
> ``|i - j| / (k - 1)`` and ``"quadratic"`` uses ``(i - j)^2 / (k - 1)^2``. The
> labels must be sortable into their ordinal order. Quadratic weighting is the common
> choice and, for a square table, coincides with an ICC-style measure.

## coherence

### `coherence(x, y, segment_length=None, overlap=0.5)`  _function_

> Magnitude-squared coherence ``|Pxy|^2 / (Pxx Pyy)`` in ``[0, 1]`` per frequency.
>
> The frequency-domain squared correlation: ``1`` where the signals are perfectly
> linearly related at that frequency, ``0`` where unrelated. Returns ``(freqs, coh)``.
> Must average several segments (default ``segment_length = n // 8``) to be meaningful.

### `cross_spectral_density(x, y, segment_length=None, overlap=0.5)`  _function_

> Welch cross-spectral density of two equal-length signals.
>
> Returns ``(freqs, cross)`` where ``cross[k]`` is the complex averaged cross-spectrum
> at one-sided normalized frequency ``freqs[k]`` in ``[0, 0.5]``. Its magnitude shows
> shared power; its phase, the frequency-dependent lead/lag.

## coin_change

### `count_change(coins, target)`  _function_

> Number of distinct multisets of ``coins`` (unlimited supply) summing to ``target``.
>
> Order does not matter (``1 + 2`` and ``2 + 1`` count once). ``target == 0`` returns 1
> (the empty multiset).

### `min_coins(coins, target)`  _function_

> Fewest coins summing to ``target`` with unlimited supply. Returns ``(count, multiset)``.
>
> ``coins`` are positive integer denominations. Returns ``(-1, [])`` if the target cannot
> be made. ``target == 0`` returns ``(0, [])``. The multiset is sorted ascending.

### `subset_sum(values, target)`  _function_

> Whether a subset of ``values`` (each used once) sums to ``target``. Returns ``(bool, subset)``.
>
> ``values`` are non-negative integers. On success the second element is one witnessing
> subset (as a list of the chosen values); on failure it is ``[]``. ``target == 0`` is
> always reachable by the empty subset.

## cointegration

### `adf_test(y, lags=0)`  _function_

> Augmented Dickey-Fuller test statistic for a unit root (constant, no trend).
>
> Parameters
> ----------
> y : sequence of float
>     The series to test.
> lags : int
>     Number of lagged differences to include (the ``augmented`` part), to soak
>     up serial correlation in the residuals.
>
> Returns
> -------
> dict
>     ``{"statistic", "rho", "reject_1pct", "reject_5pct", "reject_10pct"}``.
>     A more negative statistic is stronger evidence of stationarity; the
>     booleans compare it to the Dickey-Fuller critical values.

### `engle_granger(y, x, lags=0)`  _function_

> Engle-Granger cointegration test between two series.
>
> Regresses ``y`` on ``x`` (with an intercept), then runs :func:`adf_test` on
> the residual spread. Rejecting the unit root in the residual means ``y`` and
> ``x`` are cointegrated.
>
> Returns
> -------
> dict
>     ``{"hedge_ratio", "intercept", "adf", "cointegrated_5pct"}`` where ``adf``
>     is the residual ADF result and ``cointegrated_5pct`` is its 5% rejection.

## color

### `hex_to_rgb(code)`  _function_

> Convert a ``"#rrggbb"`` (or ``"rrggbb"``) hex string to RGB floats in ``[0, 1]``.

### `hsl_to_rgb(h, s, l)`  _function_

> Convert ``(hue_deg, saturation, lightness)`` to RGB (each in ``[0, 1]``).

### `hsv_to_rgb(h, s, v)`  _function_

> Convert ``(hue_deg, saturation, value)`` to RGB (each in ``[0, 1]``).

### `rgb_to_hex(r, g, b)`  _function_

> Convert RGB (each in ``[0, 1]``) to a ``"#rrggbb"`` hex string (rounded, clamped).

### `rgb_to_hsl(r, g, b)`  _function_

> Convert RGB (each in ``[0, 1]``) to ``(hue_deg, saturation, lightness)``.

### `rgb_to_hsv(r, g, b)`  _function_

> Convert RGB (each in ``[0, 1]``) to ``(hue_deg, saturation, value)``.
>
> Hue in ``[0, 360)`` (0 for gray), saturation and value in ``[0, 1]``.

## combinatorics

### `bell(n)`  _function_

> Bell number ``B(n)``: the total number of partitions of an ``n``-element set.
>
> Computed by summing Stirling numbers of the second kind, ``B(n) = sum_k S(n, k)``,
> via the Bell triangle. ``B(0) = 1``.

### `binomial(n, k)`  _function_

> Binomial coefficient ``C(n, k)`` -- the number of ``k``-subsets of ``n`` items.
>
> Zero when ``k < 0`` or ``k > n``. Exact for all non-negative ``n`` (big integers).

### `catalan(n)`  _function_

> Catalan number ``C_n = binomial(2n, n) / (n + 1)``.
>
> Counts balanced-parenthesis strings, binary trees, monotone lattice paths, and many
> other structures. Exact.

### `derangements(n)`  _function_

> Number of derangements ``D(n)`` -- permutations of ``n`` items with no fixed point.
>
> Uses the recurrence ``D(n) = (n-1) * (D(n-1) + D(n-2))`` with ``D(0) = 1``,
> ``D(1) = 0``. Exact.

### `multinomial(counts)`  _function_

> Multinomial coefficient ``(sum counts)! / prod(counts!)``.
>
> The number of distinct arrangements of a multiset with the given group ``counts``
> (e.g. ``multinomial([1, 4, 4]) == 630`` for the letters of "mississippi" is built
> from these). All counts must be non-negative.

### `partition_count(n)`  _function_

> Number of integer partitions of ``n`` (unordered sums of positive integers).
>
> ``p(0) = 1``; ``p(4) = 5`` (4, 3+1, 2+2, 2+1+1, 1+1+1+1). Uses the standard
> coin-change DP over parts ``1..n``, so it is exact and ``O(n^2)``.

### `stirling_second(n, k)`  _function_

> Stirling number of the second kind ``S(n, k)``: partitions of ``n`` items into ``k`` non-empty blocks.
>
> Uses the recurrence ``S(n, k) = k*S(n-1, k) + S(n-1, k-1)``. ``S(0, 0) = 1``.

## combinatorics_rank

### `combination_rank(combo, n)`  _function_

> Lexicographic rank of a ``k``-subset of ``[0, n)`` (inverse of :func:`combination_unrank`).
>
> ``combo`` is a sorted (or sortable) iterable of distinct indices in ``[0, n)``. Returns
> an integer in ``[0, C(n, k))``.

### `combination_unrank(rank, n, k)`  _function_

> The ``rank``-th ``k``-subset of ``[0, n)`` in lexicographic order.
>
> Uses the combinatorial number system. ``rank`` in ``[0, C(n, k))``. Returns a sorted
> list of ``k`` indices.

### `gray_code(n)`  _function_

> Reflected binary Gray code of a non-negative integer ``n`` (``n XOR (n >> 1)``).
>
> Consecutive Gray codes differ in exactly one bit. The inverse is
> :func:`gray_decode`.

### `gray_decode(g)`  _function_

> Inverse of :func:`gray_code`: recover the integer from its Gray code.

### `permutation_rank(perm)`  _function_

> Lexicographic rank of a permutation of ``[0, n)`` (inverse of :func:`permutation_unrank`).
>
> ``perm`` must be a permutation of ``0..n-1``. Returns an integer in ``[0, n!)``.

### `permutation_unrank(rank, n)`  _function_

> The ``rank``-th permutation of ``[0, n)`` in lexicographic order (Lehmer code).
>
> ``rank`` in ``[0, n!)``. Returns a list. ``permutation_unrank(0, n)`` is the identity.

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

### `crack_spread_option(crude_forward, product_forwards, product_weights, strike, sigma_crude, sigma_products, corr_products_crude, r, expiry, is_call=True)`  _function_

> Refinery crack-spread option (normal-model, handles negative spreads).
>
> The crack spread is ``sum_i w_i P_i - crude`` -- the refining margin from a
> barrel of crude yielding weighted refined products (e.g. 3:2:1 = 3 crude ->
> 2 gasoline + 1 heating oil). The weighted product basket is aggregated to a
> single forward and volatility (lognormal moment-free: treat the sum as one
> normal leg with variance ``sum_ij w_i w_j sigma_i sigma_j rho_ij``), then a
> Bachelier spread option is priced against the crude leg.
>
> Parameters
> ----------
> crude_forward, sigma_crude : the crude leg forward and (normal) volatility.
> product_forwards, product_weights, sigma_products : per-product forwards,
>     yield weights, and normal volatilities.
> corr_products_crude : correlation of each product with crude (list).
> strike : the strike on the crack spread.
> r, expiry : discount rate and maturity.
> is_call : call on the refining margin if True.
>
> Assumes products are mutually perfectly correlated within the basket (a common
> simplification for a refinery's co-moving product slate); returns the Bachelier
> spread-option value. Reduces to :func:`bachelier_spread_option` for a single
> unit-weight product.

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

## compensated

### `accurate_dot(a, b)`  _function_

> Compensated dot product of two equal-length vectors.

### `kahan_sum(values)`  _function_

> Kahan compensated summation of ``values``.
>
> Maintains a compensation term for the low-order bits dropped at each addition, so
> the result is far more accurate than the naive running sum for long or
> poorly-scaled sequences.

### `neumaier_sum(values)`  _function_

> Neumaier (improved Kahan) summation.
>
> Like :func:`kahan_sum` but also correct when an individual term exceeds the running
> total in magnitude -- it accumulates the correction from whichever operand is
> larger. The most robust simple compensated sum.

### `welford(values)`  _function_

> Welford's stable one-pass mean and sample variance.
>
> Returns ``(mean, variance, n)`` with the sample (``n-1``) variance. Numerically
> stable for large-mean, small-spread data where ``mean(x^2) - mean(x)^2`` loses all
> precision to cancellation. Raises on an empty input.

## complex_step

### `complex_step_derivative(f, x, h=1e-20)`  _function_

> First derivative ``f'(x)`` by the complex-step method.
>
> ``f`` must accept a complex argument and be analytic near ``x``. Returns the
> derivative ``Im(f(x + i h)) / h``; because there is no subtraction of nearby
> values, the tiny default ``h = 1e-20`` gives essentially machine-precision
> accuracy. Raises ``TypeError`` (propagated) if ``f`` is not complex-safe.

### `complex_step_gradient(f, x, h=1e-20)`  _function_

> Gradient of a scalar ``f`` of a vector ``x`` by the complex-step method.
>
> Perturbs each coordinate by ``i h`` in turn. ``f`` must accept a list/sequence of
> (possibly complex) components and return a complex-safe scalar. Returns the list of
> partial derivatives, each to near machine precision.

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

## compression

### `huffman_codebook(data)`  _function_

> Optimal Huffman prefix-code table ``{symbol: bitstring}`` for a sequence.
>
> Builds the code by repeatedly merging the two least-frequent nodes. A single distinct
> symbol maps to ``"0"`` (a one-bit code). Raises on empty input. The codes are
> prefix-free, so no code is a prefix of another.

### `huffman_decode(bits, codebook)`  _function_

> Decode a Huffman ``bits`` string given its ``codebook`` back to the symbol list.
>
> Inverts :func:`huffman_encode` exactly. Raises if the bitstring is not a valid
> concatenation of codes.

### `huffman_encode(data)`  _function_

> Huffman-encode a sequence: return ``(bitstring, codebook)``.
>
> ``bitstring`` is a ``str`` of ``'0'``/``'1'``; ``codebook`` maps each symbol to its
> code (needed to decode). Expected length is minimal among prefix codes.

### `run_length_decode(pairs)`  _function_

> Expand ``(symbol, count)`` pairs back into the original sequence (a list).
>
> Inverts :func:`run_length_encode`. Raises on a non-positive count.

### `run_length_encode(data)`  _function_

> Run-length encode a sequence into a list of ``(symbol, count)`` pairs.
>
> Consecutive equal symbols collapse into one pair; ``count >= 1``. Empty input yields
> an empty list.

## concentration

### `effective_number_of_bets(weights, cov)`  _function_

> Meucci's effective number of bets from the weights and covariance.
>
> Diagonalizes the covariance into uncorrelated principal-component factors, splits
> the portfolio variance into each factor's contribution ``p_k`` (summing to one),
> and returns ``exp(-sum p_k ln p_k)`` -- the exponential of the entropy of those
> contributions. Counts independent risk sources: ``n`` when the variance is spread
> evenly across uncorrelated factors, and down toward 1 when one factor dominates.

### `effective_number_of_constituents(weights)`  _function_

> Effective number of holdings ``1 / HHI``.
>
> The count of equally-weighted positions with the same concentration. Equals the
> number of holdings at equal weight and approaches 1 as weight concentrates.

### `herfindahl_index(weights)`  _function_

> Herfindahl-Hirschman concentration index ``sum w_i^2``.
>
> Weights need not be normalized; they are normalized to sum to one first
> (absolute values, for long-short books). Ranges from ``1/n`` (equal weight) to
> ``1`` (a single holding).

## conjugate

### `beta_binomial_posterior(prior_alpha, prior_beta, successes, trials)`  _function_

> Beta-Binomial conjugate update for a success probability ``p``.
>
> Prior ``p ~ Beta(prior_alpha, prior_beta)`` and ``successes`` out of ``trials`` Bernoulli
> outcomes give posterior ``Beta(alpha + s, beta + (n - s))``. Returns the posterior ``alpha``,
> ``beta``, ``mean`` and ``var``.

### `gamma_poisson_posterior(prior_shape, prior_rate, total_count, n_obs)`  _function_

> Gamma-Poisson conjugate update for a Poisson rate ``lambda``.
>
> Prior ``lambda ~ Gamma(shape, rate)`` (rate = inverse scale) and ``n_obs`` observations with
> summed count ``total_count`` give posterior ``Gamma(shape + total_count, rate + n_obs)``.
> Returns posterior ``shape``, ``rate``, ``mean`` and ``var``.

### `hpd_interval(density, lo, hi, mass=0.95, n_grid=10000)`  _function_

> Highest-posterior-density interval of a 1-D ``density`` over ``[lo, hi]``.
>
> Evaluates ``density`` on a uniform grid, then finds the shortest interval (by lowering a
> horizontal threshold on the density until the enclosed probability reaches ``mass``). Returns
> ``(low, high)``. The density need not be normalized. Suits any of the conjugate posteriors
> (pass the corresponding pdf) or an arbitrary unimodal 1-D posterior.

### `normal_normal_posterior(prior_mean, prior_var, data, data_var)`  _function_

> Normal-Normal conjugate update for a mean ``mu`` with known observation variance.
>
> Prior ``mu ~ Normal(prior_mean, prior_var)`` and ``data`` (list) drawn with known variance
> ``data_var`` give a Normal posterior. Precisions add: ``1/post_var = 1/prior_var + n/data_var``
> and the posterior mean is the precision-weighted average of prior mean and sample mean.
> Returns posterior ``mean`` and ``var`` (``var`` is the variance of ``mu``, not of the data).

## conjugate_gradient

### `conjugate_gradient(A, b, x0=None, tol=1e-10, max_iter=None)`  _function_

> Solve a symmetric positive-definite system ``A x = b`` by conjugate gradient.
>
> Returns a dict with ``x`` (solution), ``residual_norm`` (``||b - A x||``) and
> ``n_iter``. ``A`` must be symmetric positive-definite for convergence; converges in
> at most ``n`` iterations in exact arithmetic.

### `gauss_seidel(A, b, tol=1e-10, max_iter=None)`  _function_

> Gauss-Seidel iteration for ``A x = b`` (uses freshly-updated components).
>
> Converges for diagonally dominant or SPD ``A``. Returns ``x``, ``residual_norm``
> and ``n_iter``.

### `jacobi(A, b, tol=1e-10, max_iter=None)`  _function_

> Jacobi iteration for ``A x = b`` (all components updated from the old iterate).
>
> Converges for diagonally dominant ``A``; slower than Gauss-Seidel. Returns ``x``,
> ``residual_norm`` and ``n_iter``.

## continued_fraction

### `best_rational(x, max_denominator=1000000)`  _function_

> Closest fraction ``(p, q)`` to ``x`` with ``q <= max_denominator``.
>
> Returns the convergent (or the appropriate *semiconvergent*) with the largest
> admissible denominator -- the best rational approximation under the bound, matching
> the classic Stern-Brocot / ``limit_denominator`` result. ``q >= 1`` always.

### `cf_expansion(x, max_terms=64, tol=1e-15)`  _function_

> Continued-fraction coefficients ``[a0, a1, a2, ...]`` of a real ``x``.
>
> Repeatedly takes the integer part and inverts the remainder. Stops after
> ``max_terms`` terms or when the fractional remainder falls below ``tol`` (a rational
> ``x`` terminates). ``a0`` may be negative or zero; every later term is a positive
> integer.

### `convergents(terms)`  _function_

> Successive convergents ``[(p0, q0), (p1, q1), ...]`` from CF coefficients.
>
> Uses the standard recurrence ``p_k = a_k p_{k-1} + p_{k-2}``,
> ``q_k = a_k q_{k-1} + q_{k-2}``. Each ``(p, q)`` is the best rational approximation
> with denominator ``<= q`` -- the fractions converge to the value, alternating above
> and below it.

## convertible_lattice

### `convertible_bond_lattice(S, sigma, face, conversion_ratio, coupon_rate, maturity, r, steps=200, call_price=None, put_price=None, credit_spread=0.0, q=0.0)`  _function_

> Convertible bond value on a CRR equity tree by backward induction.
>
> ``S`` current stock, ``conversion_ratio`` shares per bond, coupons at rate
> ``coupon_rate`` on ``face`` spread evenly across the steps. At each node the
> holder value is ``max(continuation, conversion_ratio * S_node, put_price)`` and
> the issuer caps it at ``call_price`` (both optional). Continuation discounts at
> ``r + credit_spread`` (risky). ``q`` is the dividend yield. Returns the
> time-zero convertible price.

## convolution

### `convolve(a, b)`  _function_

> Full linear convolution of ``a`` and ``b`` via the FFT.
>
> Returns a real list of length ``len(a) + len(b) - 1`` -- equivalently the
> coefficients of the product of the two polynomials with coefficients ``a`` and
> ``b``. ``O(N log N)`` where ``N`` is the padded power-of-two length.

### `fft_autocorrelation(x, max_lag=None)`  _function_

> Autocorrelation of ``x`` at lags ``0..max_lag`` via the FFT.
>
> Uses the Wiener-Khinchin route (spectrum times its conjugate) on the
> mean-subtracted series, then normalizes by lag-0 so ``acf[0] = 1``. ``max_lag``
> defaults to ``len(x) - 1``. Matches a direct autocovariance sum but in
> ``O(n log n)``.

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

### `vasicek_loss_expected_shortfall(q, pd, rho)`  _function_

> Expected shortfall (average loss beyond the ``q`` quantile) of the pool loss.
>
> The mean fractional loss conditional on exceeding the ``q``-quantile. In the
> Vasicek limit this has the closed form (Tasche 2002)
>
>     ES_q = Phi_2( Phi^{-1}(pd), -Phi^{-1}(q); sqrt(rho) ) / (1 - q),
>
> where ``Phi_2(., .; r)`` is the standard bivariate normal CDF with correlation
> ``r``. Always at least the quantile :func:`vasicek_loss_quantile`, and bounded by
> the mean loss ``pd`` from below and 1 from above. Increasing in ``q`` and ``rho``.

### `vasicek_loss_pdf(loss, pd, rho)`  _function_

> Density of the large-pool loss fraction (Vasicek limit).
>
> Differentiating :func:`vasicek_loss_cdf` gives
>
>     f(x) = sqrt((1 - rho) / rho)
>            * exp( 0.5 y^2 - 0.5 ((sqrt(1 - rho) y - Phi^{-1}(pd)) / sqrt(rho))^2 ),
>
> where ``y = Phi^{-1}(x)``. The distribution is bimodal for ``rho > 0.5`` and
> concentrates at ``x = pd`` as ``rho -> 0``. ``loss`` is a fraction in ``(0, 1)``.

### `vasicek_loss_quantile(q, pd, rho)`  _function_

> Portfolio loss at confidence ``q`` (the Vasicek/Basel capital formula).
>
> Inverse of :func:`vasicek_loss_cdf`:
>
>     L(q) = Phi( (Phi^{-1}(pd) + sqrt(rho) Phi^{-1}(q)) / sqrt(1 - rho) ).
>
> The worst-case loss not exceeded with probability ``q`` -- the basis of the
> Basel IRB capital charge. Increasing in ``q``, ``pd`` and ``rho``.

## copula_sample

### `gaussian_copula_sample(correlation, n, seed=1234567)`  _function_

> Draw ``n`` samples from a Gaussian copula with the given correlation matrix.
>
> Parameters
> ----------
> correlation : list[list[float]]
>     A symmetric positive-definite correlation matrix (unit diagonal).
> n : int
>     Number of sample vectors to draw.
> seed : int
>     Seed for the deterministic normal stream.
>
> Returns
> -------
> list[list[float]]
>     ``n`` vectors of uniforms in (0, 1); each margin is uniform and the
>     cross-margin rank correlation approximates ``correlation``.

### `inverse_transform(u, ppf)`  _function_

> Map copula uniforms to a target margin via its inverse CDF ``ppf``.
>
> ``u`` is a sequence of uniforms in (0, 1); ``ppf`` maps a probability to a
> quantile (e.g. :func:`quantforge.norm_ppf` for a normal margin). Returns the
> transformed sample.

## copula_stats

### `kendall_tau(x, y)`  _function_

> Kendall's rank correlation tau-a between paired samples.
>
> Counts concordant minus discordant pairs over all ``n(n-1)/2`` pairs. Equals
> +1 for a strictly increasing relationship, -1 for strictly decreasing, and ~0
> under independence. This is the tau-a variant (no tie correction).

### `pseudo_observations(x)`  _function_

> Empirical-copula pseudo-observations: scaled ranks in (0, 1).
>
> Maps each value to ``rank / (n + 1)``, the standard normalization that keeps
> the transformed sample strictly inside the open unit interval. Applying this
> to each margin and viewing the joint gives the empirical copula, free of the
> marginal distributions.

### `spearman_rho(x, y)`  _function_

> Spearman's rank correlation: the Pearson correlation of the ranks.
>
> Invariant to any monotone transform of either margin; +1/-1 for a perfectly
> monotone relationship, ~0 under independence.

## cordic

### `cordic_atan2(y, x)`  _function_

> Return ``atan2(y, x)`` in ``(-pi, pi]`` by CORDIC vectoring mode.

### `cordic_hypot(x, y)`  _function_

> Return ``sqrt(x^2 + y^2)`` by CORDIC vectoring mode (no square root).

### `cordic_sincos(theta)`  _function_

> Return ``(cos(theta), sin(theta))`` by CORDIC rotation mode.
>
> ``theta`` is reduced into ``[-pi/2, pi/2]`` (the algorithm's convergence range) using the
> identities for the outer quadrants. Accurate to ~1e-10 with 40 iterations.

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

## correlation_test

### `pearson_correlation_test(x, y, confidence=0.95)`  _function_

> Test ``H0: rho = 0`` and give a Fisher-z confidence interval for ``rho``.
>
> Returns a dict with ``r`` (the sample correlation), ``t_stat`` and ``df`` of the
> two-sided t-test, ``p_value``, and ``conf_int`` ``[low, high]`` at ``confidence``
> from Fisher's z-transform. The p-value is small when the correlation is unlikely
> to be zero; the interval is clamped to ``[-1, 1]`` and, at ``|r| = 1``, collapses
> to the point.

### `pearson_r(x, y)`  _function_

> Pearson product-moment correlation coefficient of paired samples.
>
> ``r = cov(x, y) / (sd(x) sd(y))`` in ``[-1, 1]``. Raises if either sample has
> zero variance.

## count_min

### `BloomFilter(capacity=1000, error_rate=0.01)`  _class_

> Bloom filter for approximate set membership (no false negatives).
>
> ``capacity`` expected items and target ``error_rate`` size the bit array and hash
> count optimally. :meth:`add` inserts, :meth:`contains` tests: a member always
> returns True; a non-member returns True only with probability ~``error_rate``.

### `CountMinSketch(width=1024, depth=5)`  _class_

> Count-Min sketch for streaming frequency estimation.
>
> ``width`` counters per row, ``depth`` rows. The estimate never under-counts; the
> over-count is at most about ``total_added / width`` with probability
> ``1 - (1/2)^depth``. Add occurrences with :meth:`add` (optionally a count), query
> with :meth:`estimate`.

## cppi

### `cppi_path(initial_wealth, floor, multiplier, risky_returns, r, dt)`  _function_

> Simulate a CPPI wealth path over a sequence of risky-asset returns.
>
> Each step: allocate :func:`risky_exposure` to the risky asset (rest at the safe
> rate ``r`` over ``dt``), apply that period's ``risky_return``, and roll forward.
> The floor is discounted to each step's remaining horizon. Returns the list of
> period-end wealths. Wealth stays at or above the floor for a multiplier within
> the gap-risk limit.

### `cushion(wealth, floor_pv)`  _function_

> Cushion: wealth above the (discounted) floor, ``max(wealth - floor_pv, 0)``.

### `discounted_floor(floor, r, horizon)`  _function_

> Present value of the guaranteed floor: ``floor * e^{-r * horizon}``.
>
> The bond floor CPPI must stay above today so the terminal wealth is at least
> ``floor``.

### `risky_exposure(wealth, floor_pv, multiplier)`  _function_

> CPPI risky-asset exposure ``clamp(multiplier * cushion, 0, wealth)``.
>
> The dollar amount in the risky asset: the multiplier times the
> :func:`cushion`, capped at total wealth (no leverage) and floored at zero.
> Zero once wealth hits the floor, protecting the guarantee.

## cramer_von_mises

### `cramer_von_mises_2samp(a, b)`  _function_

> Two-sample Cramer-von Mises test.
>
> Returns a dict with the ``statistic`` (Anderson's ``T``) and the asymptotic
> ``p_value`` from the limiting Cramer-von Mises distribution. A small p-value
> rejects the null that ``a`` and ``b`` are drawn from the same distribution.

## credibility

### `buhlmann_k(expected_process_variance, variance_of_hypothetical_means)`  _function_

> Buhlmann stiffness ``k = EPV / VHM``.
>
> ``EPV`` is the mean within-risk variance; ``VHM`` the between-risk variance of
> the true means. Smaller ``k`` -> more credibility to the individual.

### `buhlmann_premium(own_mean, collective_mean, n, epv, vhm)`  _function_

> Buhlmann credibility premium.
>
> Blends the risk's own mean with the collective mean using ``Z = n/(n+k)``,
> ``k = EPV/VHM``. Lies between the two means; approaches the own mean as ``n``
> grows or the between-risk spread dominates the noise.

### `buhlmann_straub_premium(claims, exposures, collective_mean, epv, vhm)`  _function_

> Buhlmann-Straub premium for a risk with per-period exposures.
>
> Generalizes Buhlmann to unequal exposures ``m_i``: the credibility uses total
> exposure ``m = sum m_i`` with ``Z = m / (m + k)`` and the own estimate is the
> exposure-weighted claim rate ``sum claims_i / m``.

### `credibility_factor(n, k)`  _function_

> Buhlmann credibility ``Z = n / (n + k)`` for ``n`` observations.

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

## cross_correlation

### `cross_correlation(x, y, max_lag=None)`  _function_

> Cross-correlation of ``x`` and ``y`` at lags ``-max_lag..max_lag`` via the FFT.
>
> Returns ``(lags, values)`` where ``values[i]`` is the raw cross-correlation
> ``sum_n x[n + lags[i]] y[n]`` at ``lags[i]``. A positive lag slides ``x`` forward
> relative to ``y``; the peak lag is where the two series line up best. ``max_lag``
> defaults to ``len - 1``. Both inputs must be the same non-empty length.

### `lag_at_max_correlation(x, y, max_lag=None)`  _function_

> Lag (in samples) at which ``x`` and ``y`` are most positively correlated.
>
> The argmax of :func:`normalized_cross_correlation`. A positive result means ``x``
> *leads* ``y`` by that many samples (shift ``x`` forward to align them); a negative
> result means ``x`` lags ``y``. Also returns the peak correlation coefficient as
> ``(lag, coefficient)``.

### `normalized_cross_correlation(x, y, max_lag=None)`  _function_

> Cross-correlation coefficient at each lag: dimensionless, in ``[-1, 1]``.
>
> Subtracts each series' mean and divides by ``sqrt(var_x * var_y) * n`` so the value
> is a correlation coefficient -- ``1`` at the lag of perfect alignment, ``0`` for
> uncorrelated series -- independent of the signals' amplitudes. Returns
> ``(lags, values)``.

## cross_validation

### `cross_val_score(X, y, fit_fn, score_fn, k=5, shuffle=False, seed=1234567)`  _function_

> K-fold cross-validated scores for a caller-supplied model.
>
> ``fit_fn(X_train, y_train)`` returns a fitted model or predictor; ``score_fn``
> is called as ``score_fn(model, X_test, y_test)`` and returns a scalar. Returns
> the list of ``k`` per-fold scores.

### `k_fold_indices(n, k=5, shuffle=False, seed=1234567)`  _function_

> Yield ``(train_indices, test_indices)`` for ``k``-fold cross-validation.
>
> The ``n`` observations are split into ``k`` contiguous folds (of sizes
> differing by at most one); with ``shuffle=True`` the order is permuted first.
> The test folds are disjoint and cover every index exactly once.

### `train_test_split(n, test_fraction=0.2, shuffle=False, seed=1234567)`  _function_

> Single train/test split of ``n`` indices by ``test_fraction``.
>
> Returns ``(train_indices, test_indices)``. The test set gets
> ``round(n * test_fraction)`` indices (at least 1, at most ``n-1``).

## curran_asian

### `curran_asian(forward, strike, sigma, r, expiry, n_avg, is_call=True)`  _function_

> Curran arithmetic-average Asian option price.
>
> Parameters mirror :func:`~quantforge.turnbull_wakeman_asian`. Returns the
> discounted option value; call and put satisfy
> ``C - P = e^{-rT}(forward - strike)``.

## cycle_detection

### `brent_cycle(f, x0, max_iter=10000000)`  _function_

> Return ``(mu, lam)`` for ``x -> f(x)`` from ``x0`` (Brent's algorithm).
>
> Same result as :func:`floyd_cycle` but typically fewer function evaluations: it compares
> against a checkpoint whose distance doubles each phase to find ``lam`` first, then ``mu``.

### `cycle_elements(f, x0)`  _function_

> Return the list of states forming the cycle reached by iterating ``f`` from ``x0``.

### `floyd_cycle(f, x0, max_iter=10000000)`  _function_

> Return ``(mu, lam)`` for the iteration ``x -> f(x)`` from ``x0`` (Floyd's algorithm).
>
> ``mu`` is the index of the first element on the cycle (tail length); ``lam`` is the cycle
> length. Uses two pointers at speed 1 and 2. Raises if no cycle is found within
> ``max_iter`` steps (only possible for an unbounded state space).

## dag_paths

### `dag_longest_path(graph, source, target=None)`  _function_

> Longest-path distances from ``source`` on a DAG (the critical path).
>
> Same arguments as :func:`dag_shortest_path`. Unreachable nodes have distance
> ``-inf``.

### `dag_shortest_path(graph, source, target=None)`  _function_

> Shortest-path distances from ``source`` on a DAG (negative weights allowed).
>
> ``graph`` is ``{node: [(neighbor, weight), ...]}``. Returns the distance dict; if
> ``target`` is given, returns ``(distance, path)`` for that target instead
> (``distance`` is ``inf`` and ``path`` empty if unreachable).

### `transitive_closure(graph)`  _function_

> Return ``{node: set(reachable nodes)}`` (excluding the node itself unless it loops).
>
> ``graph`` may be ``{node: [(neighbor, weight), ...]}`` or ``{node: [neighbor, ...]}``.
> Computed by a DFS from each node; works on any directed graph (cycles allowed).

## dagostino

### `dagostino_k2(values)`  _function_

> D'Agostino-Pearson K^2 omnibus normality test.
>
> Returns a dict with the ``k2`` statistic (chi-square, 2 df), its ``p_value``, and
> the component ``z_skew`` and ``z_kurt`` standard scores. A small p-value rejects
> normality; the components show whether skew, tails, or both drive it. Needs at least
> 20 observations for the transforms to be reliable.

## daycount

### `day_count(start, end) -> int`  _function_

> Actual number of days between two ``(y, m, d)`` dates.

### `year_fraction(start, end, convention='act/365') -> float`  _function_

> Year fraction between ``start`` and ``end`` under ``convention``.
>
> Both dates are ``(year, month, day)`` tuples with ``end >= start``. Negative
> intervals raise.

## dbscan

### `dbscan(X, eps, min_samples)`  _function_

> Cluster points ``X`` by DBSCAN; return a label per point.
>
> ``eps`` is the neighbourhood radius, ``min_samples`` the core-point threshold
> (counting the point itself). Labels are ``0, 1, ...`` for clusters and ``-1`` for
> noise. No cluster count is required; clusters may be non-convex.

## dct

### `dct(x)`  _function_

> Orthonormal DCT-II of a real sequence.
>
> ``X[k] = s(k) sum_n x[n] cos(pi (2n+1) k / (2N))`` with the orthonormal scaling
> ``s(0) = sqrt(1/N)``, ``s(k>0) = sqrt(2/N)``. Energy-preserving; pair with
> :func:`idct`.

### `idct(X)`  _function_

> Orthonormal inverse DCT (DCT-III) -- exact inverse of :func:`dct`.
>
> ``x[n] = sum_k s(k) X[k] cos(pi (2n+1) k / (2N))`` with the same orthonormal scaling.
> Recovers the original sequence to machine precision.

## ddsketch

### `DDSketch(alpha=0.01)`  _class_

> Relative-error quantile sketch over positive values.
>
> ``add(x)`` folds one positive value; ``quantile(q)`` returns the estimated
> ``q``-quantile with relative error at most ``alpha``. ``merge`` / ``+`` combine two
> sketches built with the same ``alpha``. Memory grows only with the log-ratio of the
> largest to smallest value seen, not with the number of points.

## debruijn

### `de_bruijn_sequence(alphabet_size, n)`  _function_

> Return a de Bruijn sequence ``B(k, n)`` as a list of symbols (length ``k^n``).
>
> Every length-``n`` string over ``{0..k-1}`` occurs exactly once as a cyclic substring.
> Built by the FKM algorithm: concatenate the Lyndon words whose length divides ``n``, in
> lexicographic order.

### `is_lyndon(word)`  _function_

> True if ``word`` (a sequence) is a Lyndon word: strictly smaller than all rotations.

### `lyndon_words(alphabet_size, max_len)`  _function_

> Generate all Lyndon words over ``{0..alphabet_size-1}`` of length ``<= max_len``.
>
> Returns a list of tuples in lexicographic order (Duval's algorithm). A Lyndon word is a
> nonempty string strictly smaller than all its proper rotations.

## decision_stump

### `fit_decision_stump(X, y)`  _function_

> Fit a decision stump: the best single-feature threshold split.
>
> Scans every feature and every midpoint between adjacent sorted values, scoring
> each candidate by the size-weighted Gini of the two sides. Returns
> ``{"feature", "threshold", "left_label", "right_label", "gini"}`` where points
> with ``x[feature] <= threshold`` take ``left_label``.

### `gini_impurity(labels)`  _function_

> Gini impurity ``1 - sum_c p_c^2`` of a label list (0 = pure).

### `predict_decision_stump(stump, X_query)`  _function_

> Predict labels for query rows under a fitted stump.

## decision_tree

### `fit_decision_tree(X, y, max_depth=5, min_samples=2)`  _function_

> Fit a CART classification tree.
>
> Parameters
> ----------
> X : list[list[float]]
>     ``n`` rows of ``d`` features.
> y : list
>     Class labels.
> max_depth : int
>     Maximum tree depth (>= 1). Depth 1 is a decision stump.
> min_samples : int
>     Minimum samples required to attempt a split.
>
> Returns
> -------
> dict
>     Nested tree; internal nodes have ``feature``/``threshold``/``left``/
>     ``right``, leaves have only ``prediction``.

### `predict_decision_tree(tree, X_query)`  _function_

> Predict labels for query rows by walking the tree to a leaf.

### `tree_depth(tree)`  _function_

> Depth of a fitted tree (a single leaf has depth 0).

## delaunay

### `delaunay_triangulation(points)`  _function_

> Return the Delaunay triangulation of ``points`` as a list of index triples.
>
> Each triple ``(i, j, k)`` indexes into ``points``. Needs at least three
> non-collinear points. Duplicate points are ignored for the triangulation but indices
> refer to the original list.

## deming

### `deming_regression(x, y, lam=1.0)`  _function_

> Deming (errors-in-variables) regression of ``y`` on ``x``.
>
> ``lam`` is the ratio of the error variances ``var(eps_x) / var(eps_y)`` (default 1,
> i.e. orthogonal regression). Returns ``(slope, intercept)`` for the fitted line
> ``y = slope * x + intercept``. Requires a non-zero ``Sxy``.

### `orthogonal_regression(x, y)`  _function_

> Orthogonal (total-least-squares) regression: Deming with ``lambda = 1``.
>
> Minimizes the perpendicular distances to the line, treating ``x`` and ``y``
> symmetrically. Returns ``(slope, intercept)``.

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

## dfa

### `dfa_exponent(x, scales=None)`  _function_

> Detrended-fluctuation scaling exponent ``alpha``.
>
> The slope of ``log F(s)`` regressed on ``log s`` from :func:`dfa_fluctuations`.
> ``alpha ~ 0.5`` for white noise, ``> 0.5`` for persistent series, ``< 0.5`` for
> anti-persistent; it equals the Hurst exponent for a stationary long-memory
> series and exceeds it by 1 for the integrated (random-walk) version.

### `dfa_fluctuations(x, scales=None)`  _function_

> Fluctuation function ``F(s)`` of DFA over a set of window sizes ``scales``.
>
> Returns ``(scales, fluctuations)``. ``scales`` defaults to a dyadic-ish grid
> between 4 and ``len(x)//4``. Windows that do not divide the series exactly drop
> the remainder. Requires at least 16 points.

## differential_evolution

### `differential_evolution(func, bounds, pop_size=None, F=0.8, cr=0.9, max_iter=1000, tol=1e-12, seed=1234567)`  _function_

> Minimize ``func`` over box ``bounds`` by differential evolution.
>
> ``bounds`` is a list of ``(lo, hi)`` per dimension. ``func`` takes a length-``d``
> list and returns a scalar. ``pop_size`` defaults to ``max(15, 10*d)``; ``F`` is the
> differential weight, ``cr`` the crossover rate. Returns a dict with ``x`` (best
> vector), ``fun`` (its objective), ``n_iter`` and ``converged`` (population spread
> below ``tol``). Deterministic for a fixed ``seed``.

## diophantine

### `linear_diophantine(a, b, c)`  _function_

> Solve ``a x + b y = c`` in integers, or return ``None`` if unsolvable.
>
> Returns ``(x0, y0, dx, dy)``: one particular solution ``(x0, y0)`` and the step
> ``(dx, dy)`` so that ``(x0 + k*dx, y0 + k*dy)`` is a solution for every integer ``k``.
> Solvable iff ``gcd(a, b)`` divides ``c``. ``a`` and ``b`` must not both be zero.

### `pell_fundamental(n)`  _function_

> Smallest positive ``(x, y)`` solving Pell's equation ``x^2 - n y^2 = 1``.
>
> Built from the convergents of ``sqrt(n)``'s continued fraction. ``n`` must be a
> non-square positive integer.

### `sqrt_continued_fraction(n)`  _function_

> Continued-fraction expansion of ``sqrt(n)`` for a non-square ``n``: ``(a0, period)``.
>
> ``sqrt(n) = [a0; a1, a2, ...]`` is eventually periodic; returns the integer part
> ``a0`` and the repeating block ``period`` (a list). Raises if ``n`` is a perfect
> square or negative.

## discount_curve

### `DiscountCurve(times, dfs)`  _class_

> Log-linear discount curve from pillar ``(T, DF)`` points.

### `SplineZeroCurve(times, zero_rates)`  _class_

> Discount curve with cubic-spline-interpolated continuously-compounded zeros.
>
> Interpolates the pillar zero rates with a natural cubic spline (smooth C2 zero
> curve) rather than log-linear discount factors, so the instantaneous forward
> rates are smooth. Exposes the same ``df``/``zero_rate``/``forward_rate``
> interface as :class:`DiscountCurve` and reprices the pillars exactly.

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

### `displaced_diffusion_implied_vol(target_price, S, K, t, r, shift=0.0, option_type=<OptionType.CALL: 'call'>, b=None, tol=1e-10, max_iter=200)`  _function_

> Implied displaced-diffusion volatility from a market price.
>
> Inverts :func:`displaced_diffusion_price` for the ``sigma`` reproducing
> ``target_price`` by bisection (the price is monotone increasing in ``sigma``).
> With ``shift = 0`` this coincides with the Black-Scholes implied vol. Raises if
> the quote lies outside the attainable ``[intrinsic, forward]`` band.

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

## distance_correlation

### `distance_correlation(x, y)`  _function_

> Distance correlation ``dCor(x, y)`` in ``[0, 1]``.
>
> ``dCov(x, y) / sqrt(dVar(x) dVar(y))``. Zero exactly under independence (unlike
> Pearson, this holds for nonlinear dependence too), 1 for a tight linear relation.
> Returns 0 if either sample is constant (distance variance zero).

### `distance_covariance(x, y)`  _function_

> Distance covariance ``dCov(x, y)`` (the square root of the mean product).
>
> Zero if and only if ``x`` and ``y`` are independent (in the population). Returns a
> non-negative value on the same scale as the data.

### `distance_variance(x)`  _function_

> Distance variance ``dVar(x) = dCov(x, x)``; zero only for a constant sample.

## distortion

### `expected_loss(pmf)`  _function_

> Expected loss ``sum_k k * pmf[k]`` (= sum of the survival function).

### `proportional_hazard_premium(pmf, rho)`  _function_

> Proportional-hazard distortion premium ``sum_k S_k^{1/rho}``.
>
> ``rho = 1`` gives the expected loss; ``rho > 1`` loads the tail. Requires
> ``rho >= 1``.

### `wang_premium(pmf, lam)`  _function_

> Wang-transform premium of a discrete loss distribution.
>
> Distorts the survival function by ``g(u) = Phi(Phi^{-1}(u) + lam)`` and sums
> it over the grid (unit spacing). ``lam = 0`` gives the expected loss; positive
> ``lam`` adds a risk load. Monotone increasing in ``lam``.

## distributions

### `binomial_cdf(k, n, prob)`  _function_

> Binomial CDF ``P(X <= k)`` for ``n`` trials with success probability ``prob``.
>
> Uses the incomplete-beta identity ``P(X <= k) = I_{1-p}(n - k, k + 1)``, exact
> and stable for large ``n``. ``k`` is floored to an integer.

### `binomial_pmf(k, n, prob)`  _function_

> Binomial probability mass ``C(n, k) p^k (1-p)^{n-k}``.

### `chi2_cdf(x, df)`  _function_

> CDF of the chi-square distribution with ``df`` degrees of freedom.
>
> A gamma with ``shape = df/2`` and ``scale = 2``: ``F(x) = P(df/2, x/2)``.

### `chi2_ppf(p, df)`  _function_

> Quantile of the chi-square distribution for ``p`` in ``(0, 1)``.

### `chi2_sf(x, df)`  _function_

> Survival function ``1 - chi2_cdf(x, df)`` (upper tail, for p-values).

### `f_cdf(x, d1, d2)`  _function_

> CDF of the F distribution with ``(d1, d2)`` degrees of freedom.
>
> ``F(x) = I_{d1 x / (d1 x + d2)}(d1/2, d2/2)`` via the regularized incomplete
> beta.

### `f_ppf(p, d1, d2)`  _function_

> Quantile of the F distribution for ``p`` in ``(0, 1)``.

### `gamma_cdf(x, shape, scale=1.0)`  _function_

> CDF of the gamma distribution ``Gamma(shape, scale)`` at ``x``.
>
> ``F(x) = P(shape, x / scale)`` via the regularized lower incomplete gamma.
> Reduces to the exponential CDF ``1 - e^{-x/scale}`` when ``shape = 1``.

### `gamma_pdf(x, shape, scale=1.0)`  _function_

> Density of the gamma distribution at ``x >= 0``.

### `gamma_ppf(p, shape, scale=1.0)`  _function_

> Quantile (inverse CDF) of the gamma distribution for ``p`` in ``(0, 1)``.

### `poisson_cdf(k, lam)`  _function_

> Poisson CDF ``P(N <= k)`` via the gamma relation ``= Q(k+1, lam)``.
>
> Uses ``P(N <= k) = gammaincc(k + 1, lam)`` (the regularized upper incomplete
> gamma), exact and stable for large ``lam`` where summing masses would lose
> precision. ``k`` is floored to an integer.

### `poisson_pmf(k, lam)`  _function_

> Poisson probability mass ``P(N = k) = e^{-lam} lam^k / k!``.

## divergences

### `bhattacharyya_distance(p, q)`  _function_

> Bhattacharyya distance ``-log(sum sqrt(p_i q_i))`` (0 for identical, grows apart).

### `hellinger_distance(p, q)`  _function_

> Hellinger distance ``(1/sqrt2) sqrt(sum (sqrt(p_i) - sqrt(q_i))^2)`` in ``[0, 1]``.

### `jensen_shannon_divergence(p, q)`  _function_

> Jensen-Shannon divergence (symmetric, bounded by ``log 2``) in nats.
>
> ``0.5 KL(p || m) + 0.5 KL(q || m)`` with ``m = (p + q)/2``. Always finite; its
> square root is a metric.

### `kl_divergence(p, q)`  _function_

> Kullback-Leibler divergence ``sum p_i log(p_i / q_i)`` in nats.
>
> Asymmetric and non-negative (zero iff ``p == q``). Raises if some ``q_i == 0`` where
> ``p_i > 0`` (the divergence is infinite there).

### `total_variation_distance(p, q)`  _function_

> Total-variation distance ``(1/2) sum |p_i - q_i|`` in ``[0, 1]``.

## double_barrier

### `double_knockin_call(S, K, L, U, t, r, sigma, b=None, q=0.0, terms=8)`  _function_

> Price a double-barrier knock-in call in closed form.
>
> A double knock-in call comes alive only if the spot touches either barrier
> (``L`` or ``U``) at some point before expiry; otherwise it expires worthless.
> It is valued by the in-out parity
>
>     knock-in + knock-out = vanilla,
>
> since exactly one of "the corridor is breached" and "the corridor is never
> breached" occurs on every path. Parameters and defaults match
> :func:`double_knockout_call`.
>
> Returns
> -------
> float
>     Value of the knock-in call. Non-negative, never exceeds the vanilla call,
>     approaches zero as the barriers move far away (a breach becomes rare), and
>     rises toward the vanilla as the corridor tightens (a breach becomes
>     certain).

### `double_knockout_call(S, K, L, U, t, r, sigma, b=None, q=0.0, terms=8)`  _function_

> Price a double-barrier knock-out call in closed form (Kunitomo-Ikeda).
>
> Parameters
> ----------
> S, K : float
>     Spot and strike.
> L, U : float
>     Lower and upper knock-out barriers with ``L < S < U``. Both must be
>     positive and ``L < U``.
> t, r, sigma : float
>     Time to expiry (years), risk-free rate, volatility.
> b : float, optional
>     Cost of carry. Defaults to ``r - q``.
> q : float
>     Continuous dividend yield, used only when ``b`` is not given.
> terms : int
>     Number of image terms on each side of the series (total ``2*terms+1``).
>
> Returns
> -------
> float
>     Value of the down-and-out-and-up-and-out call. Non-negative and never
>     exceeds the vanilla call; it approaches the vanilla as the barriers move
>     far away.

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

## downside_ratios

### `kappa_ratio(returns, tau=0.0, order=2)`  _function_

> Kaplan-Knowles Kappa of the given ``order`` about target ``tau``.
>
> ``kappa_n = (mean - tau) / LPM_n^(1/n)``. ``order=1`` is the Omega-Sharpe
> ratio, ``order=2`` the per-period Sortino ratio. Raises if there is no
> downside (zero lower partial moment).

### `lower_partial_moment(returns, tau=0.0, order=2)`  _function_

> ``n``-th lower partial moment about ``tau``: mean(max(tau - r, 0)^n).

### `upside_potential_ratio(returns, tau=0.0)`  _function_

> Upside-potential ratio: mean upside over downside deviation about ``tau``.
>
> ``UPR = mean(max(r - tau, 0)) / sqrt(LPM_2(tau))`` -- expected outperformance
> of the target per unit of downside risk. Raises if there is no downside.

## drawdown

### `drawdown_analytics(returns)`  _function_

> Depth-and-duration summary of a return series' drawdowns.
>
> Returns a dict with ``max_drawdown_depth`` (fractional, non-negative),
> ``peak_index`` / ``trough_index`` / ``recovery_index`` of the deepest episode
> (indices into the equity curve, which has ``len(returns)+1`` points),
> ``time_to_recovery`` (periods from trough to recovery, or ``None`` if never),
> and ``longest_underwater`` (the most periods spent below a prior peak).
> ``recovery_index`` is ``None`` if the deepest drawdown never recovers.

## dual

### `Dual(value, deriv=0.0)`  _class_

> A dual number ``value + deriv * eps`` for forward-mode autodiff.
>
> Construct a variable as ``Dual(x, 1.0)`` and a constant as ``Dual(c)``. Arithmetic
> and the module's elementary functions propagate the derivative exactly via the
> chain rule; read ``.value`` and ``.deriv`` from the result.

### `dual_derivative(f, x)`  _function_

> Exact derivative ``f'(x)`` by forward-mode autodiff.
>
> ``f`` must accept a :class:`Dual` and return a :class:`Dual`, built from ``Dual``
> arithmetic and this module's elementary functions (``exp``, ``log``, ``sin``, ...).
> Returns ``f'(x)`` with no truncation error.

## dual_calculus

### `dual_gradient(f, x)`  _function_

> Exact gradient of a scalar ``f`` of a vector ``x`` by forward-mode autodiff.
>
> ``f`` takes a list of (possibly :class:`Dual`) components and returns a scalar
> ``Dual``. Seeds each coordinate's derivative to 1 in turn, so the result is the
> vector of partial derivatives with no truncation error. ``n`` function evaluations.

### `dual_newton(f, x0, tol=1e-12, max_iter=100)`  _function_

> Newton's method for ``f(x) = 0`` using autodiff for the derivative.
>
> ``f`` must accept a :class:`Dual` and return a :class:`Dual`. Each step evaluates
> ``f`` once on a seeded dual, reading both ``f(x)`` and ``f'(x)`` at no extra cost,
> so no separate derivative function is needed. Returns a dict with ``root``,
> ``iterations`` and ``converged``. Raises if the derivative vanishes.

## dual_currency

### `dcd_breakeven_spot(coupon_rate, base_deposit_rate, tenor, strike, deposit_ccy_is_base=True)`  _function_

> Spot at which the DCD matches a plain deposit (breakeven).
>
> The converted DCD value equals the plain-deposit value when the extra coupon
> exactly offsets the conversion loss. For a base deposit (converted above
> ``strike``):
>
>     strike * (1 + coupon*tenor) / breakeven = 1 + base*tenor
>     => breakeven = strike * (1 + coupon*tenor) / (1 + base*tenor).
>
> Above the strike (in the converted region) for a positive yield pickup.

### `dcd_enhanced_yield(base_deposit_rate, option_premium_rate, tenor)`  _function_

> Enhanced annual yield of a DCD: base rate plus annualized premium.
>
> ``base_deposit_rate + option_premium_rate / tenor`` -- the option premium
> (as a fraction of notional, earned once) spread over the ``tenor`` in years and
> added to the plain deposit rate. Above the base deposit rate whenever the sold
> option has value.

### `dcd_maturity_payoff(notional, coupon_rate, tenor, spot_at_maturity, strike, deposit_ccy_is_base=True)`  _function_

> Base-currency value of a DCD at maturity, including conversion.
>
> The depositor always earns ``notional * (1 + coupon_rate * tenor)`` of
> principal-plus-coupon; if the option the depositor sold finishes in the money
> the repayment is converted at ``strike`` rather than ``spot_at_maturity``,
> costing the depositor the shortfall. Returns the base-currency value received.
>
> For a base-currency deposit (sold call), conversion bites when
> ``spot_at_maturity > strike``: repayment happens at the worse ``strike``, so
> value = principal_plus_coupon * strike / spot when converted.

### `dcd_option_premium_rate(spot, strike, tenor, r_domestic, r_foreign, sigma, deposit_ccy_is_base=True)`  _function_

> Premium rate (fraction of notional) of the option embedded in a DCD.
>
> The depositor sells a call on the base currency (if depositing the base
> currency) or a put (if depositing the alternate), struck at the conversion
> ``strike``, priced with Garman-Kohlhagen (carry ``b = r_domestic -
> r_foreign``). Returned per unit notional. Higher vol and a nearer strike raise
> the premium, and thus the DCD's yield pickup.

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

## dunn_test

### `dunn_test(groups, adjust='holm')`  _function_

> Dunn's post-hoc pairwise test after Kruskal-Wallis.
>
> ``groups`` is a sequence of samples. Returns a list of pairwise result dicts, one
> per unordered pair ``(i, j)``, each with ``groups`` ``(i, j)``, the ``z`` statistic,
> the raw ``p_value`` and the ``p_adjusted`` value. ``adjust`` is ``"holm"``
> (default), ``"bonferroni"`` or ``None`` for no correction. Uses one pooled ranking
> with the tie correction, so it is consistent with the Kruskal-Wallis H.

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

## ecdf

### `ecdf(data, x)`  _function_

> Empirical CDF ``F_n(x) = #{data_i <= x} / n`` at point(s) ``x``.
>
> ``x`` may be a scalar (returns a float in ``[0, 1]``) or an iterable (returns a
> list). Right-continuous step function.

### `qq_points(sample, reference)`  _function_

> Quantile-quantile pairing of two samples for a Q-Q plot.
>
> Returns a list of ``(reference_quantile, sample_quantile)`` pairs evaluated at the
> plotting positions ``(i - 0.5)/m`` of the smaller sample size ``m``. On the y = x
> line the two samples share a distribution; a slope or curvature departure reveals a
> scale or shape difference. ``reference`` may be another sample.

### `quantile(data, p, method='linear')`  _function_

> Sample quantile at probability ``p`` in ``[0, 1]``.
>
> ``method`` selects the interpolation on the order statistics: ``"linear"`` (the
> NumPy default, ``(n-1)p`` position), ``"lower"``, ``"higher"``, or ``"nearest"``.
> ``p`` may be a scalar or iterable. Matches the standard percentile conventions.

## eigen_general

### `characteristic_polynomial(A)`  _function_

> Coefficients of ``det(xI - A)`` (monic, highest degree first) by Faddeev-LeVerrier.
>
> Returns ``[1, c_{n-1}, ..., c_0]`` -- the characteristic polynomial with leading
> coefficient 1. Also the input to :func:`polynomial_roots` for the eigenvalues. ``A``
> must be square.

### `determinant_from_charpoly(A)`  _function_

> Determinant of ``A`` from its characteristic polynomial's constant term.
>
> ``det(A) = (-1)^n * c_0`` where ``c_0`` is the constant term of ``det(xI - A)``. A
> cheap exact-by-construction cross-check against the LU determinant.

### `eigenvalues_general(A, tol=1e-10)`  _function_

> All eigenvalues of a general real (or complex) square matrix.
>
> Builds the characteristic polynomial via Faddeev-LeVerrier and finds its roots with
> Durand-Kerner, so complex-conjugate eigenvalue pairs are returned correctly. Real
> eigenvalues come back as Python ``float``; genuinely complex ones as ``complex``.
> Sorted by real part then imaginary part.

## ekf

### `extended_kalman_filter(observations, f, h, Q, R, x0, P0)`  _function_

> Extended Kalman filter over ``observations`` for nonlinear ``f`` and ``h``.
>
> ``f`` maps a length-``n`` list of :class:`Var` (the state) to a length-``n`` list -- the
> deterministic state transition. ``h`` maps the state to a length-``m`` list -- the
> measurement model. ``Q`` (n x n), ``R`` (m x m) are the process and observation covariances,
> ``x0`` (n) and ``P0`` (n x n) the initial mean and covariance. Returns a dict with
> ``filtered_means`` and ``filtered_covariances``.
>
> The transition/observation Jacobians are computed exactly at each step by reverse-mode
> autodiff, so no analytic linearization is required.

## elliptic

### `agm(a, b, tol=1e-15, max_iter=100)`  _function_

> Arithmetic-geometric mean of ``a, b >= 0``: the common limit of the AM/GM iteration.
>
> Each step replaces ``(a, b)`` with ``((a+b)/2, sqrt(a b))``; the two converge
> quadratically. ``agm(a, b) == agm(b, a)`` and ``agm(a, a) == a``.

### `elliptic_e(m)`  _function_

> Complete elliptic integral of the second kind ``E(m)``, ``0 <= m <= 1``.
>
> ``E(0) = pi/2`` and ``E(1) = 1``. Uses the AGM descent, accumulating the geometric-step
> corrections that give ``E`` from the same iteration.

### `elliptic_k(m)`  _function_

> Complete elliptic integral of the first kind ``K(m)``, ``0 <= m < 1``.
>
> ``K(0) = pi/2`` and ``K(m) -> inf`` as ``m -> 1``. Computed as
> ``pi / (2 AGM(1, sqrt(1 - m)))``.

## encoding

### `fit_label_encoder(values)`  _function_

> Fit a label encoder: sorted distinct categories -> integer indices.
>
> Returns ``{"categories": [...]}`` where the list position is the code.

### `label_decode(encoder, codes)`  _function_

> Invert label codes back to categories; -1 (unseen) maps to ``None``.

### `label_encode(encoder, values)`  _function_

> Encode ``values`` to integer codes; unseen categories map to -1.

### `one_hot_encode(encoder, values)`  _function_

> One-hot encode ``values`` against a fitted label encoder.
>
> Each row is a length-``len(categories)`` list of 0/1; a known category sets one
> entry to 1 (rows sum to 1), an unseen category yields an all-zero row.

## energy_distance

### `energy_distance(a, b)`  _function_

> Energy distance ``2 A - B - C`` between two samples.
>
> Non-negative; zero only when the empirical distributions coincide. ``A`` is the
> mean cross-sample absolute distance, ``B`` and ``C`` the mean within-sample
> distances (each normalized by the squared sample size, i.e. including the zero
> diagonal, per Szekely-Rizzo).

### `energy_test(a, b, n_permutations=999, seed=1234567)`  _function_

> Permutation test of equal distributions via the energy distance.
>
> Pools the two samples, reshuffles the group labels ``n_permutations`` times, and
> returns a dict with the observed ``statistic``, the ``p_value``
> ``(1 + #{perm >= observed}) / (1 + n_permutations)``, and ``n_permutations``.
> A small p-value is evidence the two samples come from different distributions.

## entropy_pooling

### `entropy_pooling_mean(x, target, prior=None, tol=1e-12, max_iter=200)`  _function_

> Posterior scenario probabilities matching a target mean, min relative entropy.
>
> Parameters
> ----------
> x : sequence of float
>     Scenario values of the quantity being viewed.
> target : float
>     Desired posterior mean ``E_p[x]``. Must lie strictly inside
>     ``(min(x), max(x))`` -- a mean outside the scenario range is infeasible.
> prior : sequence of float, optional
>     Prior probabilities (non-negative, summing to 1). Defaults to uniform.
> tol, max_iter : float, int
>     Bisection tolerance on the achieved mean and its iteration cap.
>
> Returns
> -------
> list[float]
>     Posterior probabilities ``p`` (an exponential tilt of the prior). They sum
>     to 1, reproduce ``target`` as their mean, and reduce to the prior when
>     ``target`` equals the prior mean.

### `relative_entropy(p, q)`  _function_

> Kullback-Leibler divergence ``sum_i p_i log(p_i / q_i)`` (nats).
>
> Non-negative, and zero exactly when ``p == q``. Terms with ``p_i == 0`` are
> dropped (limit ``0 log 0 = 0``).

## entropy_ts

### `approximate_entropy(series, m=2, r=None)`  _function_

> Approximate entropy ``ApEn(m, r)`` (Pincus).
>
> ``r`` is the tolerance for a match (defaults to ``0.2 * std(series)``). Returns
> ``phi_m - phi_{m+1}``; near zero for a perfectly regular series and larger for an
> irregular one. Includes self-matches, so it is biased low on short series -- use
> :func:`sample_entropy` when that matters.

### `permutation_entropy(series, m=3, normalize=True)`  _function_

> Permutation entropy (Bandt-Pompe) from the ordinal patterns of length ``m``.
>
> Slides a window of ``m`` points, records which permutation sorts each window, and
> returns the Shannon entropy of the pattern distribution. With ``normalize=True``
> it is divided by ``log(m!)`` to land in ``[0, 1]`` -- 0 for a monotone series, 1
> for one whose orderings are uniformly random. Invariant to any monotone transform
> of the series.

### `sample_entropy(series, m=2, r=None)`  _function_

> Sample entropy ``SampEn(m, r)`` (Richman-Moorman).
>
> ``-log(A / B)`` where ``B`` counts template matches of length ``m`` and ``A`` of
> length ``m + 1``, both excluding self-matches. Larger means less regular / more
> complex. ``r`` defaults to ``0.2 * std(series)``. Raises if no longer-template
> matches occur (entropy would be infinite).

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

## eulerian

### `eulerian_path(graph, directed=False, start=None)`  _function_

> Return a list of vertices tracing an Eulerian path, or ``None`` if none exists.
>
> The returned list has ``E + 1`` entries (each consecutive pair is a traversed edge).
> For a circuit the first and last vertices coincide. ``graph`` is ``{u: [v, ...]}``; for
> an undirected graph list each edge once (both directions are inferred).

### `has_eulerian_circuit(graph, directed=False)`  _function_

> True if ``graph`` has an Eulerian circuit (closed Eulerian path).

### `has_eulerian_path(graph, directed=False)`  _function_

> True if ``graph`` has an Eulerian path (trail using every edge once).

## evt

### `gev_cdf(x, loc, scale, shape)`  _function_

> Generalized extreme value CDF (block-maxima limit distribution).
>
> ``F(x) = exp(-(1 + shape (x-loc)/scale)^{-1/shape})`` for ``shape != 0`` (Frechet
> ``shape > 0`` / Weibull ``shape < 0``), and the Gumbel limit
> ``exp(-e^{-(x-loc)/scale})`` as ``shape -> 0``. Defined where
> ``1 + shape (x-loc)/scale > 0``.

### `gev_fit_block_maxima(block_maxima)`  _function_

> Method-of-moments GEV fit assuming the Gumbel (shape = 0) limit.
>
> Fits location and scale of a Gumbel to the block maxima by moment matching:
> ``scale = std * sqrt(6)/pi``, ``loc = mean - gamma * scale`` (``gamma`` the
> Euler-Mascheroni constant). Returns ``(loc, scale, shape=0.0)`` -- a simple,
> robust baseline when the shape is not separately estimated.

### `gev_return_level(period, loc, scale, shape)`  _function_

> Return level: the block maximum exceeded once per ``period`` blocks.
>
> Inverts :func:`gev_cdf` at ``p = 1 - 1/period``:
>
>     level = loc + (scale/shape) [ (-ln(1 - 1/period))^{-shape} - 1 ],
>
> with the Gumbel limit ``loc - scale ln(-ln(1 - 1/period))`` as ``shape -> 0``.
> The T-block return level rises with the return period.

### `gpd_expected_shortfall(losses, threshold, confidence=0.99)`  _function_

> Peaks-over-threshold expected shortfall from the GPD tail.
>
> ``ES = VaR/(1 - xi) + (beta - xi threshold)/(1 - xi)`` for ``xi < 1``. Always
> at least the :func:`gpd_var`; finite only for a tail with ``xi < 1``.

### `gpd_fit_pot(losses, threshold)`  _function_

> Fit a Generalized Pareto to peaks over ``threshold`` (method of moments).
>
> Returns ``(xi, beta, n_exceed, n_total)``: the GPD shape ``xi`` and scale
> ``beta`` matching the mean and variance of the exceedances ``X - threshold``,
> the number of exceedances, and the sample size. ``xi > 0`` is a heavy
> (Pareto-type) tail.

### `gpd_var(losses, threshold, confidence=0.99)`  _function_

> Peaks-over-threshold VaR from a fitted Generalized Pareto tail.
>
> ``VaR = threshold + (beta/xi) [ (n/N_u (1 - confidence))^{-xi} - 1 ]`` with
> ``N_u`` exceedances of ``n`` observations. Reduces to the exponential-tail
> limit as ``xi -> 0``. A positive loss quantile deep in the tail.

### `hill_estimator(losses, k)`  _function_

> Hill estimator of the tail index from the top ``k`` order statistics.
>
> ``xi = (1/k) sum_{i=1}^{k} ln(X_(n-i+1) / X_(n-k))`` -- the mean log-excess of
> the ``k`` largest losses over the ``(k+1)``-th. Estimates the shape ``xi`` of a
> heavy power-law tail (tail exponent ``alpha = 1/xi``); larger ``xi`` means a
> heavier tail.

## ewma

### `EWMAStats(lam=0.94, values=None)`  _class_

> Streaming exponentially-weighted mean, variance and volatility.
>
> ``lam`` is the decay in ``(0, 1)`` (larger = longer memory). Feed points with
> :meth:`update`; read :attr:`mean`, :meth:`variance`, :meth:`std`. The first point
> seeds the mean; variance builds from the second. Effective window is about
> ``1 / (1 - lam)`` observations.

### `ewma(values, lam=0.94)`  _function_

> Exponentially-weighted moving mean of a whole series (list output).
>
> Returns the running EWMA at each step, seeded from the first value. A convenience
> wrapper over :class:`EWMAStats` for batch use.

## ewma_cov

### `ewma_correlation_matrix(returns, lam=0.94)`  _function_

> EWMA correlation matrix: the EWMA covariance normalized by its diagonal.
>
> Unit diagonal, off-diagonals in ``[-1, 1]``.

### `ewma_covariance_matrix(returns, lam=0.94)`  _function_

> EWMA covariance matrix of a return panel (rows = periods, cols = assets).
>
> Recurses ``Sigma_t = lam Sigma_{t-1} + (1-lam) r_t r_t'`` from the sample
> covariance seed, treating returns as mean-zero (the RiskMetrics convention).
> Returns the final ``p x p`` symmetric positive-semidefinite matrix.

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

## expint

### `e1(x)`  _function_

> Exponential integral ``E1(x) = integral_x^inf e^-t / t dt`` for ``x > 0``.

### `ei(x)`  _function_

> Exponential integral ``Ei(x)`` (principal value) for real ``x != 0``.

### `en(n, x)`  _function_

> Generalized exponential integral ``E_n(x) = integral_1^inf e^{-x t} / t^n dt``.
>
> ``n >= 0`` integer, ``x >= 0`` (``x > 0`` when ``n <= 1``). Uses ``E_0 = e^-x / x`` and
> the upward recurrence ``E_{n}(x) = (e^-x - x E_{n-1}(x)) / (n - 1)`` seeded from ``E_1``.

## expression

### `eval_expression(expr)`  _function_

> Evaluate an infix arithmetic string safely (no Python ``eval``).

### `eval_rpn(rpn)`  _function_

> Evaluate a postfix (RPN) token list to a number.

### `shunting_yard(tokens)`  _function_

> Convert an infix token list to a postfix (RPN) list (Dijkstra's shunting-yard).

### `tokenize(expr)`  _function_

> Split ``expr`` into a token list of numbers (float), operators, parens, and names.

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

## factor_model

### `factor_attribution(total_return, alpha, betas, factor_realized_returns)`  _function_

> Decompose a realized return into factor contributions plus a residual.
>
> Each factor contributes ``beta_k * factor_realized_return_k``; the residual is
> ``total_return - alpha - sum(contributions)`` (the part not explained by the
> factors, i.e. the model residual for the period). Returns a dict with
> ``factor_contributions`` (list), ``alpha`` (the intercept as its own term), and
> ``residual``. The alpha, factor contributions, and residual sum to
> ``total_return`` by construction.

### `factor_expected_return(alpha, betas, factor_premia)`  _function_

> Expected return from a fitted factor model: ``alpha + sum beta_k premium_k``.

### `factor_regression(asset_returns, factor_returns)`  _function_

> OLS regression of ``asset_returns`` on ``factor_returns`` (list of series).
>
> ``factor_returns`` is a list of equal-length factor return series. Returns a
> dict with ``alpha`` (intercept), ``betas`` (one loading per factor),
> ``r_squared``, and ``residual_vol`` (standard deviation of the residuals).
> Fits ``r_t = alpha + sum_k beta_k f_{k,t} + eps_t`` by minimizing the squared
> residuals.

### `rolling_factor_beta(asset_returns, factor_returns, window)`  _function_

> Rolling single-factor beta over a trailing ``window``.
>
> Runs :func:`factor_regression` on each trailing window of length ``window``,
> returning the list of first-factor betas (one per window end, from index
> ``window - 1`` onward). Tracks how the factor loading drifts through time.

## fd_weights

### `fd_weights(x0, grid, max_deriv)`  _function_

> Fornberg weights for derivatives ``0..max_deriv`` at ``x0`` on ``grid``.
>
> Returns a list ``c`` of length ``max_deriv + 1`` where ``c[k]`` is the weight list (same
> length as ``grid``) for the ``k``-th derivative: ``f^(k)(x0) ~ sum_i c[k][i] f(grid[i])``.
> Requires ``len(grid) > max_deriv`` and distinct grid points.

## fenwick

### `FenwickTree(size_or_values)`  _class_

> Binary indexed tree for point updates and prefix/range sums (``O(log n)`` each).
>
> Construct from a size (all zeros) or an initial list. ``update(i, delta)`` adds to
> element ``i``; ``prefix_sum(i)`` is the sum of elements ``0..i``; ``range_sum(lo, hi)``
> the inclusive sum over ``[lo, hi]``. Zero-indexed.

### `SegmentTree(values, combine=None, identity=0.0)`  _class_

> Segment tree for range queries under an associative ``combine`` with point updates.
>
> ``combine`` defaults to ``+`` (range sum); pass ``min``/``max`` (with the matching
> ``identity``) for range-minimum/maximum. ``update(i, value)`` sets element ``i``;
> ``query(lo, hi)`` folds ``combine`` over ``[lo, hi]`` inclusive. ``O(log n)`` each.

## fenwick2d

### `FenwickTree2D(rows=None, cols=None, grid=None)`  _class_

> Dynamic 2-D prefix sums: point-add and rectangle-sum in ``O(log R * log C)``.
>
> Construct with the grid shape ``(rows, cols)`` (optionally an initial grid). ``add(r, c,
> delta)`` adds to a cell; ``prefix_sum(r, c)`` sums the rectangle ``[0, r) x [0, c)``; and
> ``range_sum(r0, c0, r1, c1)`` sums ``[r0, r1) x [c0, c1)`` by inclusion-exclusion.

## fft

### `fft(x)`  _function_

> Forward FFT of a sequence whose length is a power of two.
>
> Returns the complex DFT ``X_k = sum_n x_n exp(-2 pi i k n / N)``. Accepts real or
> complex input; raises unless the length is a positive power of two.

### `ifft(x)`  _function_

> Inverse FFT: recovers the sequence from its DFT (with the ``1/N`` scaling).
>
> ``ifft(fft(x)) == x`` up to floating error. Length must be a power of two.

## fir_filter

### `fir_apply(taps, x)`  _function_

> Apply an FIR filter (tap list) to signal ``x`` by direct convolution.
>
> Returns the filtered signal, same length as ``x`` (each output is the tap-weighted
> sum of the current and preceding samples; the first ``len(taps)-1`` outputs are the
> transient start-up).

### `fir_bandpass(numtaps, low, high, window='hamming')`  _function_

> Windowed-sinc bandpass FIR taps passing ``[low, high]`` (normalized cutoffs).
>
> ``numtaps`` must be odd; the band-pass is the difference of two low-pass filters.

### `fir_highpass(numtaps, cutoff, window='hamming')`  _function_

> Windowed-sinc highpass FIR taps (spectral inversion of a lowpass).
>
> ``numtaps`` must be odd. Unit gain at Nyquist, zero at DC.

### `fir_lowpass(numtaps, cutoff, window='hamming')`  _function_

> Windowed-sinc lowpass FIR taps.
>
> ``cutoff`` is the normalized cutoff (0 to 1, fraction of Nyquist). ``numtaps`` should
> be odd for a symmetric (linear-phase, zero-delay-at-center) filter. Returns the tap
> list, normalized to unit DC gain.

## fisher_exact

### `fisher_exact_test(table, alternative='two-sided')`  _function_

> Fisher's exact test for a 2x2 ``table`` ``[[a, b], [c, d]]``.
>
> ``alternative`` is ``"two-sided"`` (default), ``"greater"`` or ``"less"`` (one-sided
> on the odds ratio). Returns a dict with the sample ``odds_ratio``
> (``a d / (b c)``, ``inf`` if a denominator is zero) and the ``p_value``. The
> two-sided p-value sums the probabilities of every table (given the margins) no more
> probable than the observed one.

## forecast_combine

### `combine_forecasts(forecasts, weights)`  _function_

> Weighted combination of aligned forecast series by ``weights``.

### `inverse_mse_weights(error_series)`  _function_

> Combination weights proportional to ``1 / MSE`` of each model's errors.
>
> ``error_series`` is a list of forecast-error series. A more accurate model (lower
> MSE) gets more weight; weights sum to one. Ignores cross-model error correlation.

### `optimal_combination_weights(error_series)`  _function_

> Bates-Granger minimum-variance combination weights.
>
> ``w = Sigma^{-1} 1 / (1' Sigma^{-1} 1)`` where ``Sigma`` is the covariance of the
> models' forecast errors. Minimizes the variance of the combined error and can put
> negative weight on a model that hedges another's errors. Weights sum to one.

### `simple_average_forecast(forecasts)`  _function_

> Equal-weight combination of aligned forecast series.
>
> ``forecasts`` is a list of series (one per model); returns the point-wise mean.

## forecast_metrics

### `mae(actual, forecast)`  _function_

> Mean absolute error.

### `mape(actual, forecast)`  _function_

> Mean absolute percentage error (as a fraction; 0.1 = 10%).
>
> Raises if any actual value is zero (the percentage is undefined there).

### `mase(actual, forecast, train=None, season=1)`  _function_

> Mean absolute scaled error.
>
> Scales the forecast MAE by the in-sample MAE of a seasonal-naive forecast
> (differences at lag ``season``) computed on ``train`` (defaults to ``actual``).
> ``< 1`` means the forecast beats naive; ``= 1`` matches it. Raises if the naive
> benchmark has zero error (a perfectly predictable training series).

### `rmse(actual, forecast)`  _function_

> Root mean squared error.

### `smape(actual, forecast)`  _function_

> Symmetric mean absolute percentage error, in [0, 2].
>
> ``mean( |a - f| / ((|a| + |f|) / 2) )``. Terms with both ``a`` and ``f`` zero
> contribute 0. Robust to small actuals and bounded, unlike plain MAPE.

### `theil_u1(actual, forecast)`  _function_

> Theil's U1 inequality coefficient, bounded in ``[0, 1]``.
>
> ``U1 = RMSE / (rms(actual) + rms(forecast))``. Zero for a perfect forecast and
> one in the worst case; scale-free and symmetric. Distinct from
> :func:`theil_u2`, which benchmarks against the naive forecast.

### `theil_u2(actual, forecast, last_actual=None)`  _function_

> Theil's U2 statistic: forecast RMSE over the no-change naive RMSE.
>
> The naive forecast for period ``t`` is the previous actual ``actual[t-1]`` (or
> ``last_actual`` for the first point, if given). ``U2 < 1`` means the forecast
> beats the random walk, ``= 1`` matches it, ``> 1`` is worse. Needs at least two
> points (or one point plus ``last_actual``).

## forecast_test

### `diebold_mariano(errors1, errors2, h=1, power=2)`  _function_

> Diebold-Mariano statistic and two-sided p-value for equal predictive accuracy.
>
> ``errors1`` / ``errors2`` are the forecast-error series of the two models
> (target minus forecast), ``h`` the forecast horizon (sets the number of
> autocovariance lags ``h - 1`` in the HAC variance), and ``power`` the loss
> exponent (2 = squared error, 1 = absolute error). Returns ``(DM, p_value)`` with
> the Harvey-Leybourne-Newbold small-sample correction and a Student-t reference on
> ``n - 1`` degrees of freedom. A negative ``DM`` favors the first forecast.

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

## fra

### `fra_forward_rate(discount, t1, t2)`  _function_

> Simple forward rate over ``[t1, t2]`` implied by the discount curve.

### `fra_value(discount, contract_rate, t1, t2, notional=1.0, payer=True)`  _function_

> Value today of a forward rate agreement.
>
> A payer (long the FRA, paying fixed ``contract_rate``) gains when the forward
> rate exceeds the contract rate: ``notional * tau * (f - K) * P(t2)``. A
> receiver is the negative. Zero at the fair (forward) rate.

## fracdiff

### `fixed_width_fracdiff(series, d, threshold=1e-05)`  _function_

> Fixed-width fractional differencing (Lopez de Prado).
>
> Truncates the weight window at the first ``|w_k| < threshold`` and applies that
> fixed-length filter, so every output point uses the same weights. Returns the
> valid portion (length ``len(series) - width + 1``), where ``width`` is the number
> of retained weights. Preserves stationarity with a constant memory window.

### `fracdiff_weights(d, n)`  _function_

> Binomial weights ``w_0..w_{n-1}`` of the operator ``(1 - L)^d``.
>
> ``w_0 = 1`` and ``w_k = w_{k-1} * -(d - k + 1) / k``. For integer ``d`` the
> weights vanish beyond ``k = d`` (e.g. ``d = 1`` gives ``[1, -1]`` then zeros).

### `fractional_difference(series, d)`  _function_

> Fractionally difference ``series`` by the full expansion of ``(1 - L)^d``.
>
> Returns a list the same length as ``series``; entry ``t`` is
> ``sum_{k=0}^{t} w_k * series[t - k]`` (a growing backward window, so early
> entries use fewer weights). ``d = 0`` returns the series unchanged and ``d = 1``
> returns the first difference (with the first entry equal to ``series[0]``).

## fresnel

### `dawson(x)`  _function_

> Dawson function ``D(x) = e^{-x^2} integral_0^x e^{t^2} dt``.

### `fresnel_c(x)`  _function_

> Fresnel cosine integral ``C(x) = integral_0^x cos(pi t^2 / 2) dt``.

### `fresnel_s(x)`  _function_

> Fresnel sine integral ``S(x) = integral_0^x sin(pi t^2 / 2) dt``.

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

## fuzzy_match

### `dice_coefficient(a, b)`  _function_

> Sørensen-Dice coefficient over character bigrams: ``2|A∩B| / (|A|+|B|)``.
>
> A multiset overlap of adjacent character pairs -- robust to word-order and small
> edits. Identical strings score ``1``; strings sharing no bigram score ``0``. Strings
> shorter than two characters fall back to exact equality.

### `jaccard_similarity(a, b, tokenize=None)`  _function_

> Jaccard similarity ``|A∩B| / |A∪B|`` over token (or character) sets.
>
> By default splits on whitespace into a word set; pass ``tokenize`` (a callable
> returning an iterable of tokens, e.g. ``list`` for characters) to change the
> granularity. Two empty token sets score ``1``.

### `jaro(a, b)`  _function_

> Jaro similarity of two strings in ``[0, 1]`` (1 = identical).
>
> Counts matching characters within a sliding window and penalizes transpositions,
> the classic short-string metric. Empty-vs-empty is ``1``; empty-vs-nonempty is ``0``.

### `jaro_winkler(a, b, prefix_weight=0.1, max_prefix=4)`  _function_

> Jaro-Winkler similarity: Jaro boosted for a shared prefix.
>
> Adds ``prefix_weight * L * (1 - jaro)`` where ``L`` is the common prefix length (up
> to ``max_prefix``), rewarding strings that agree at the start -- the standard tweak
> for names. ``prefix_weight`` must satisfy ``prefix_weight * max_prefix <= 1``.

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

### `cross_rate(ab, cb, via_is_quote=True)`  _function_

> Cross rate from two rates sharing a common currency.
>
> Given ``A/B`` and ``C/B`` (both quoted against the common currency ``B``), the
> cross ``A/C = (A/B) / (C/B)`` when ``via_is_quote`` (the common currency is the
> quote of both), else ``A/C = (A/B) * (B/C)`` for ``A/B`` and ``B/C``. Positive
> rates required.

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

### `is_arbitrage_free(ab, bc, ca, tol=1e-09)`  _function_

> True if the triangular loop is arbitrage-free within ``tol``.
>
> Checks ``|triangular_arbitrage - 1| <= tol``; the cross rates are mutually
> consistent when the round-trip product is one.

### `triangular_arbitrage(ab, bc, ca)`  _function_

> Triangular-arbitrage profit factor around a currency loop ``A->B->C->A``.
>
> Converting one unit of A through ``A/B``... actually multiplying the three
> quoted legs ``(A per B) (B per C) (C per A)`` returns the units of A after a
> round trip; it equals one in an arbitrage-free market. Returns the product;
> values above one (net of costs) signal a profitable loop, below one the reverse
> direction.

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

## gauss_hermite

### `gauss_hermite_expectation(g, mu=0.0, sigma=1.0, n=16)`  _function_

> Approximate ``E[g(X)]`` for ``X ~ N(mu, sigma^2)`` by Gauss-Hermite quadrature.
>
> Shifts and scales the standard-normal nodes to ``x_i = mu + sigma * z_i`` and
> forms ``sum_i w_i g(x_i)``. Exact when ``g`` is a polynomial of degree at most
> ``2n - 1``. With ``sigma = 0`` it collapses to ``g(mu)``.

### `gauss_hermite_nodes_weights(n)`  _function_

> Probabilists' Gauss-Hermite nodes and weights for ``n`` points.
>
> Returns ``(nodes, weights)`` for the weight function ``phi(z)`` (standard
> normal density), so ``sum_i w_i = 1`` and ``sum_i w_i x_i^{2m}`` reproduces the
> standard-normal moments. Exact for polynomials up to degree ``2n - 1``. Nodes
> are symmetric about zero and returned in increasing order.

## gauss_kronrod

### `gauss_kronrod(f, a, b, tol=1e-10, max_intervals=1000)`  _function_

> Adaptive Gauss-Kronrod (G7-K15) integral of ``f`` over ``[a, b]``.
>
> Starts with one panel and repeatedly bisects the panel with the largest local
> error estimate until the total estimated error falls below ``tol`` (or
> ``max_intervals`` panels are used). Returns the integral. Concentrates
> evaluations where the integrand is hardest, so it handles peaks and mild
> endpoint behaviour that a fixed rule would miss.

## gauss_laguerre

### `gauss_laguerre_integral(g, n=32, rate=1.0)`  _function_

> Approximate ``integral_0^inf g(x) dx`` by Gauss-Laguerre quadrature.
>
> Writes the integrand as ``g(x) = [g(x) e^{rate * x}] e^{-rate * x}`` and applies
> the ``e^{-x}`` rule after the substitution ``u = rate * x``:
>
>     integral_0^inf g(x) dx = (1/rate) sum_i w_i g(x_i / rate) e^{x_i}.
>
> ``rate`` should roughly match the integrand's exponential decay for best
> accuracy (choose ``rate`` near the true decay constant). Exact when
> ``g(x) e^{rate x}`` is a polynomial of degree up to ``2n - 1``.

### `gauss_laguerre_nodes_weights(n)`  _function_

> Gauss-Laguerre nodes and weights for the weight ``e^{-x}`` on ``[0, inf)``.
>
> Returns ``(nodes, weights)`` with ``sum_i w_i = 1`` and
> ``sum_i w_i x_i^m = m!`` (the moments of the ``e^{-x}`` density). Exact for
> polynomials up to degree ``2n - 1``. Nodes are positive and returned in
> increasing order.

## gauss_newton

### `gauss_newton(residual, p0, tol=1e-10, max_iter=200)`  _function_

> Minimize ``sum residual(p)^2`` by Levenberg-Marquardt with an exact autodiff Jacobian.
>
> ``residual`` takes a length-``k`` list of :class:`quantforge.reverse_ad.Var` and returns a
> list of ``m`` residual :class:`Var` (``m >= k``). Returns a dict with ``p`` (the fitted
> parameters), ``cost`` (half the sum of squared residuals ``0.5 * ||r||^2``), ``n_iter``,
> ``converged`` and ``grad_norm`` (norm of ``J^T r``).

## geo

### `cross_track_distance(lat, lon, lat1, lon1, lat2, lon2, radius=6371.0088)`  _function_

> Signed distance of point ``(lat, lon)`` from the great circle through 1 and 2.
>
> Positive when the point lies to the *right* of the path direction (1 -> 2), negative
> to the left; magnitude is the perpendicular great-circle distance. Useful for "how far
> off the route am I".

### `destination_point(lat, lon, bearing, distance, radius=6371.0088)`  _function_

> Destination ``(lat, lon)`` reached from a start point on a given bearing/distance.
>
> Follows the great circle from ``(lat, lon)`` heading ``bearing`` degrees for
> ``distance`` (same units as ``radius``). Returns degrees; the inverse of
> :func:`haversine_distance` / :func:`initial_bearing`.

### `haversine_distance(lat1, lon1, lat2, lon2, radius=6371.0088)`  _function_

> Great-circle distance between two lat/lon points (haversine), in ``radius`` units.
>
> Numerically stable for small distances (unlike the spherical law of cosines). Inputs
> in degrees; returns kilometres by default.

### `initial_bearing(lat1, lon1, lat2, lon2)`  _function_

> Initial (forward) bearing from point 1 to point 2 along the great circle, degrees.
>
> The compass bearing to steer at the start of the path, in ``[0, 360)`` clockwise from
> north. It changes along a great circle, so this is the *initial* heading.

## geometry

### `closest_pair(points)`  _function_

> Closest pair of points and their distance: ``(p, q, distance)``.
>
> Divide-and-conquer in ``O(n log n)``. Needs at least two points; ties break to the
> first pair found.

### `convex_hull(points)`  _function_

> Convex hull of a set of 2-D points (Andrew's monotone chain).
>
> Returns the hull vertices in counter-clockwise order, starting from the lowest-then-
> leftmost point, without repeating the first point. Collinear interior points are
> dropped. Needs at least one point; duplicates are ignored.

### `point_in_polygon(point, polygon)`  _function_

> Whether ``point`` lies inside a simple ``polygon`` (ray-casting, odd-crossing rule).
>
> Casts a ray to the right and counts edge crossings; an odd count means inside. Points
> exactly on an edge are reported as inside. ``polygon`` is an ordered ``(x, y)`` vertex
> list.

### `polygon_area(polygon)`  _function_

> Area of a simple polygon by the shoelace formula (unsigned).
>
> ``polygon`` is a list of ``(x, y)`` vertices in order (open ring; the last vertex is
> joined back to the first). Returns the absolute area. Needs at least three vertices.

### `polygon_centroid(polygon)`  _function_

> Centroid ``(cx, cy)`` of a simple polygon (area-weighted, shoelace form).
>
> Falls back to the vertex average for a degenerate (zero-area) polygon. ``polygon`` is
> an ordered ``(x, y)`` vertex list.

## geometry2

### `clip_polygon(subject, clip)`  _function_

> Clip ``subject`` polygon against a convex ``clip`` polygon (Sutherland-Hodgman).
>
> Both are ordered ``(x, y)`` vertex lists; ``clip`` must be convex and given
> counter-clockwise. Returns the clipped polygon (possibly empty) as a vertex list --
> the intersection of the subject with the clip window. The classic viewport/window
> clipping algorithm.

### `polygon_perimeter(polygon)`  _function_

> Perimeter of a polygon: the sum of its edge lengths (closed ring).

### `segment_intersection(p1, p2, p3, p4)`  _function_

> Intersection point of two segments, or ``None`` if they do not cross at a single point.
>
> Returns the ``(x, y)`` crossing for segments that meet at exactly one point; returns
> ``None`` if they are parallel, collinear, or disjoint. Endpoints count as
> intersections.

### `segments_intersect(p1, p2, p3, p4)`  _function_

> Whether segment ``p1-p2`` intersects segment ``p3-p4`` (including endpoints/collinear).
>
> Uses the standard four-orientation test with collinear-overlap handling. Returns a
> bool.

## geometry3

### `bounding_box(points)`  _function_

> Axis-aligned bounding box of a point set as ``(min_x, min_y, max_x, max_y)``.

### `min_enclosing_circle(points, seed=1234567)`  _function_

> Smallest circle containing all ``points`` as ``(center_x, center_y, radius)``.
>
> Welzl's algorithm with a deterministic shuffle (an LCG seeded by ``seed``), running in
> expected linear time. Every input point lies inside or on the returned circle, and the
> circle is the unique smallest such -- no smaller circle contains them all.

### `polygon_diameter(points)`  _function_

> Diameter of a point set: the farthest-apart pair and their distance.
>
> Returns ``(p, q, distance)``. Restricts the search to the convex-hull vertices (the
> diameter is always realized by two hull points), then does an all-pairs scan over the
> hull -- exact, and cheap once the hull has few vertices.

## geometry_dist

### `closest_point_on_segment(p, a, b)`  _function_

> Closest point on the finite segment ``a-b`` to ``p`` (clamped to the endpoints).

### `point_polyline_distance(p, polyline)`  _function_

> Minimum distance from ``p`` to a polyline (a list of ``>= 2`` vertices).
>
> Takes the smallest :func:`point_segment_distance` over consecutive segments.

### `point_segment_distance(p, a, b)`  _function_

> Distance from ``p`` to the nearest point of the finite segment ``a-b``.
>
> Equals the perpendicular distance when the foot of the perpendicular lands on the
> segment, otherwise the distance to the nearer endpoint.

### `point_to_line_distance(p, a, b)`  _function_

> Perpendicular distance from point ``p`` to the *infinite* line through ``a`` and ``b``.
>
> Raises if ``a == b`` (no line defined).

## gf2_linalg

### `gf2_nullspace_basis(equations, n_vars)`  _function_

> Return a basis for the nullspace ``{x : A x = 0}`` over GF(2).
>
> Each basis vector is a list of ``n_vars`` bits. The number of vectors is
> ``n_vars - rank``. An empty list means only the zero vector satisfies the system.

### `gf2_rank(rows)`  _function_

> Rank over GF(2) of a matrix given as a list of bit-packed integer rows.

### `solve_gf2(equations, rhs, n_vars)`  _function_

> Solve an XOR linear system over GF(2). Returns one solution list, or ``None``.
>
> ``equations[i]`` is a bit-packed integer whose bit ``j`` (value ``1 << j``) marks that
> variable ``j`` appears in equation ``i``; ``rhs[i]`` is that equation's right-hand side
> (0 or 1). ``n_vars`` is the number of variables. Returns a list of ``n_vars`` bits
> (free variables set to 0), or ``None`` if the system is inconsistent.

## gmm

### `fit_gaussian_mixture(data, k=2, max_iter=200, tol=1e-08, seed=1234567)`  _function_

> Fit a ``k``-component 1-D Gaussian mixture by EM.
>
> Returns a dict with ``weights`` (mixing proportions summing to one), ``means``,
> ``variances``, ``log_likelihood`` (of the final fit) and ``n_iter``. Components
> are initialized at spread-out data quantiles (deterministic given ``seed``, which
> only jitters the initial means). The log-likelihood is non-decreasing across
> iterations. Requires at least ``k`` distinct points.

## goertzel

### `goertzel(x, k)`  _function_

> DFT coefficient ``X[k]`` of ``x`` at integer bin ``k`` by the Goertzel recurrence.
>
> Returns the complex ``X[k] = sum_n x[n] exp(-2 pi i k n / N)`` -- identical to the
> ``k``-th FFT output, computed with a single ``O(n)`` real recurrence. Integer ``k``
> in ``[0, N)``.

### `goertzel_power(x, k)`  _function_

> Power ``|X[k]|^2`` at bin ``k`` (the Goertzel magnitude-squared).
>
> The efficient tone-detection quantity: large when a frequency near bin ``k`` is
> present, small otherwise. Avoids the final trig of :func:`goertzel`.

## gof_tests

### `jarque_bera_test(returns)`  _function_

> Jarque-Bera normality test: ``(statistic, p_value)``.
>
> The statistic is chi-square(2) under the normal null, so the p-value is its
> upper-tail probability. A normal sample gives a small statistic and a large
> p-value; a heavy-tailed or skewed sample gives a large statistic and a small
> p-value.

### `ks_two_sample(a, b)`  _function_

> Two-sample Kolmogorov-Smirnov test: ``(D, p_value)``.
>
> ``D`` is the maximum absolute difference between the two empirical CDFs,
> evaluated at every observed point. The p-value uses the asymptotic Kolmogorov
> distribution with the effective sample size ``sqrt(n m / (n + m))``. A small p
> rejects "the two samples come from the same distribution".

## gph

### `fractional_integrate(noise, d)`  _function_

> Fractionally integrate white ``noise`` by order ``d``: apply ``(1 - L)^{-d}``.
>
> The inverse of :func:`quantforge.fracdiff.fractional_difference`; useful for
> generating ARFIMA(0, d, 0) test series. Uses the binomial weights of
> ``(1 - L)^{-d}``, ``w_0 = 1``, ``w_k = w_{k-1} (k - 1 + d) / k``.

### `gph_estimate(x, m=None)`  _function_

> Estimate the fractional-integration order ``d`` by the GPH regression.
>
> ``x`` is the series; ``m`` is the number of low frequencies used (defaults to
> ``floor(sqrt(n))``, the common bandwidth). Returns a dict with ``d`` (the memory
> parameter), ``std_error`` (asymptotic, from the ``pi^2/6`` log-periodogram
> variance), ``m`` and ``n``. ``d ~ 0`` indicates short memory, ``d > 0`` long
> memory, ``d < 0`` anti-persistence.

## gradient_boost

### `fit_gradient_boost(X, y, n_estimators=100, learning_rate=0.1, max_depth=3, min_samples=2)`  _function_

> Fit a gradient-boosted regression-tree ensemble (squared-error loss).
>
> Returns a dict with the initial ``base`` prediction (the mean of ``y``), the list
> of ``trees``, and ``learning_rate``. Each tree is fit to the current residuals and
> contributes ``learning_rate * its prediction``. Feed the result to
> :func:`predict_gradient_boost`.

### `predict_gradient_boost(model, X_query)`  _function_

> Predict targets for ``X_query`` with a fitted gradient-boosted ensemble.

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

## graph

### `bfs(graph, source)`  _function_

> Breadth-first hop distances from ``source`` (unweighted).
>
> ``graph`` maps each node to an iterable of neighbors. Returns ``{node: hops}`` for
> every node reachable from ``source`` (``0`` at the source), in increasing hop order.

### `connected_components(graph)`  _function_

> Connected components of an undirected graph as a list of node sets.
>
> ``graph`` maps each node to its neighbors; edges are treated as undirected (a node
> listed as a neighbor is joined even if the reverse edge is absent). Returns the
> components as a list of sorted lists.

### `dijkstra(graph, source)`  _function_

> Shortest-path distances from ``source`` over non-negative edge weights.
>
> ``graph`` is ``{node: {neighbor: weight}}`` with weights ``>= 0``. Returns
> ``(distances, predecessors)``: ``distances[node]`` is the shortest total weight from
> ``source`` (``inf`` if unreachable) and ``predecessors[node]`` the previous node on a
> shortest path (``None`` at the source or if unreachable). Raises on a negative weight.

### `shortest_path(graph, source, target)`  _function_

> Shortest path ``[source, ..., target]`` and its total weight, via Dijkstra.
>
> Returns ``(path, distance)``; ``(None, inf)`` if ``target`` is unreachable.

### `topological_sort(graph)`  _function_

> Topological ordering of a DAG (Kahn's algorithm); raises if the graph has a cycle.
>
> ``graph`` is ``{node: [successors]}`` (or ``{node: {succ: weight}}``). Returns a list
> of nodes such that every edge points forward. Raises ``ValueError`` if a cycle makes
> a valid ordering impossible.

## graph2

### `UnionFind()`  _class_

> Disjoint-set forest with path compression and union by rank.
>
> ``find(x)`` returns the representative of ``x``'s set; ``union(a, b)`` merges two
> sets and returns whether they were previously disjoint; ``connected(a, b)`` tests
> membership. Elements are created on first reference. Near-constant amortized cost.

### `max_flow(graph, source, sink)`  _function_

> Maximum flow from ``source`` to ``sink`` (Edmonds-Karp).
>
> ``graph`` is ``{node: {neighbor: capacity}}`` with capacities ``>= 0``. Returns the
> maximum total flow value. Uses BFS to find shortest augmenting paths in the residual
> network, so it runs in ``O(V E^2)`` and terminates on rational capacities.

### `minimum_spanning_tree(nodes, edges)`  _function_

> Minimum spanning tree (or forest) by Kruskal's algorithm.
>
> ``nodes`` is an iterable of vertices; ``edges`` a list of ``(u, v, weight)``. Returns
> ``(tree_edges, total_weight)`` where ``tree_edges`` is the chosen subset (sorted by
> weight then endpoints). For a disconnected graph this is the minimum spanning forest.

## graph3

### `betweenness_centrality(graph, normalized=True)`  _function_

> Betweenness centrality by Brandes' algorithm (unweighted, undirected/directed).
>
> The fraction of shortest paths (between all pairs) that pass through each node.
> ``normalized`` divides by the number of pairs so scores are comparable across graph
> sizes. Returns ``{node: centrality}``.

### `closeness_centrality(graph)`  _function_

> Closeness centrality: inverse of the mean shortest-path distance (unweighted).
>
> For each node, ``(reachable) / sum(distances)`` scaled by the fraction reachable
> (Wasserman-Faust), so isolated or unreachable nodes score low. Uses BFS from every
> node. Returns ``{node: centrality}``.

### `degree_centrality(graph)`  _function_

> Degree centrality: each node's neighbor count normalized by ``n - 1``.
>
> For an undirected graph pass a symmetric adjacency. Returns ``{node: centrality}`` in
> ``[0, 1]`` (a node linked to all others scores 1).

### `pagerank(graph, damping=0.85, tol=1e-10, max_iter=1000)`  _function_

> PageRank scores by power iteration with teleportation.
>
> ``graph`` maps each node to its out-neighbors (list or ``{neighbor: weight}``).
> ``damping`` is the follow-a-link probability (``1 - damping`` teleports uniformly).
> Dangling nodes (no out-links) redistribute their mass uniformly. Returns
> ``{node: score}`` summing to 1; converges when the L1 change drops below ``tol``.

## graph4

### `a_star(graph, source, target, heuristic)`  _function_

> Point-to-point shortest path by A* search with an admissible ``heuristic``.
>
> ``graph`` is ``{node: {neighbor: weight}}`` with non-negative weights. ``heuristic``
> is a callable ``h(node)`` estimating the remaining cost to ``target``; it must never
> overestimate (admissible) for the result to be optimal. Returns ``(path, cost)``, or
> ``(None, inf)`` if ``target`` is unreachable.

### `bellman_ford(graph, source)`  _function_

> Single-source shortest paths allowing negative edge weights (Bellman-Ford).
>
> ``graph`` is ``{node: {neighbor: weight}}``; weights may be negative. Returns
> ``(distances, predecessors)`` after ``V - 1`` relaxation rounds. Raises
> ``ValueError`` if a negative-weight cycle is reachable (no finite shortest path).

### `floyd_warshall(graph)`  _function_

> All-pairs shortest distances by Floyd-Warshall (``O(V^3)``).
>
> ``graph`` is ``{node: {neighbor: weight}}`` (negative edges allowed, no negative
> cycles). Returns a nested dict ``dist[u][v]`` of shortest path weights (``inf`` if
> unreachable, ``0`` on the diagonal). Raises if a negative cycle is present (a
> diagonal entry goes negative).

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

## grid_interp

### `bilinear_interp(xs, ys, z, x, y)`  _function_

> Bilinear interpolation of grid ``z`` at ``(x, y)``.
>
> ``xs`` are the column coordinates (length ``ncols``), ``ys`` the row coordinates
> (length ``nrows``), and ``z`` an ``nrows x ncols`` value table (``z[row][col]``). Both
> axes must be strictly increasing. Out-of-range queries clamp to the nearest edge.
> Exact at grid nodes; reduces to linear interpolation along a grid line.

### `nearest_interp(xs, ys, z, x, y)`  _function_

> Nearest-neighbour lookup on the grid: the value at the closest node to ``(x, y)``.
>
> Snaps to the nearest grid coordinate on each axis (clamping out-of-range). Piecewise
> constant, so it preserves the exact sampled values -- useful for categorical or
> quantized grids.

## hac

### `autocorrelation(x, lag)`  _function_

> Sample autocorrelation ``rho_k = gamma_k / gamma_0`` at ``lag``.

### `autocovariance(x, lag)`  _function_

> Biased (divisor ``n``) sample autocovariance of ``x`` at ``lag``.
>
> ``gamma_k = (1/n) sum_{t=k}^{n-1} (x_t - xbar)(x_{t-k} - xbar)``. The divisor
> ``n`` (not ``n - k``) is what makes the autocovariance sequence positive
> semidefinite, which the Newey-West weighting relies on.

### `newey_west_mean_se(x, lags)`  _function_

> HAC (Newey-West) standard error of the sample mean of ``x``.
>
> ``se = sqrt(sigma_LR^2 / n)`` with ``sigma_LR^2`` the Newey-West long-run
> variance. For serially uncorrelated data this matches the usual
> ``sqrt(var / n)``; positive autocorrelation inflates it.

### `newey_west_variance(x, lags)`  _function_

> Newey-West long-run variance of a series (Bartlett-weighted HAC).
>
> Parameters
> ----------
> x : sequence of float
>     The series (e.g. demeaned returns or a moment condition).
> lags : int
>     Truncation lag ``L`` (>= 0). ``lags = 0`` reduces to the sample variance.
>
> Returns
> -------
> float
>     The long-run variance estimate. Always non-negative thanks to the Bartlett
>     weights.

## halley

### `halley(f, fprime, fdoubleprime, x0, tol=1e-14, max_iter=100)`  _function_

> Halley's method: cubically-convergent root of ``f`` from ``x0``.
>
> Needs ``f``, its first derivative ``fprime``, and second derivative ``fdoubleprime``.
> The update is ``x - 2 f f' / (2 f'^2 - f f'')``. Raises if the denominator vanishes or
> it fails to converge in ``max_iter`` steps.

### `secant(f, x0, x1, tol=1e-14, max_iter=200)`  _function_

> Secant method: derivative-free root of ``f`` from two initial guesses ``x0``, ``x1``.
>
> Approximates the derivative by the last two points; converges superlinearly
> (order ~1.618). Raises if the points coincide or it fails to converge.

## hamming_code

### `hamming74_decode(code)`  _function_

> Decode a 7-bit Hamming(7,4) codeword, correcting any single-bit error.
>
> Computes the 3-bit syndrome; if non-zero it gives the 1-indexed position of the
> flipped bit, which is corrected before extracting the data. Returns
> ``(data_bits, error_position)`` where ``error_position`` is ``0`` if no error was
> found or the 1-indexed position that was corrected.

### `hamming74_encode(bits)`  _function_

> Encode 4 data bits into a 7-bit Hamming(7,4) codeword.
>
> ``bits`` is a length-4 sequence of 0/1 (``d1 d2 d3 d4``). Returns a length-7 list
> ``[p1, p2, d1, p3, d2, d3, d4]`` where the parity bits ``p1, p2, p3`` cover the
> standard bit-position groups. Any single bit flip in the result is later correctable.

### `luhn_check_digit(digits)`  _function_

> Compute the Luhn check digit to append to a payload ``digits`` (without one).
>
> Returns the single digit ``0-9`` that makes ``digits + [check]`` a valid Luhn number
> (checksum zero).

### `luhn_checksum(digits)`  _function_

> Luhn checksum of a digit sequence: ``0`` iff the number is valid.
>
> ``digits`` is a string or list of decimal digits *including* the trailing check digit.
> Doubles every second digit from the right, sums the digits of the results, and returns
> the total modulo 10. A valid Luhn number gives ``0``.

## hawkes

### `hawkes_branching_ratio(alpha, beta)`  _function_

> Branching ratio ``alpha / beta`` -- expected offspring per event.
>
> Below 1 the process is stationary; at or above 1 it explodes.

### `hawkes_fit(events, t_end=None, x0=None)`  _function_

> Maximum-likelihood fit of ``(mu, alpha, beta)`` to observed event times.
>
> Maximizes :func:`log_likelihood` with Nelder-Mead over a log-parameterization
> (keeping ``mu, alpha, beta`` positive). Returns a dict with ``mu``, ``alpha``,
> ``beta``, ``branching_ratio`` and the attained ``log_likelihood``. Provide
> ``x0 = (mu, alpha, beta)`` to seed the search.

### `hawkes_gof_test(events, mu, alpha, beta)`  _function_

> Kolmogorov-Smirnov goodness-of-fit of a fitted Hawkes model.
>
> Applies the time-rescaling theorem: under a correct model the :func:`residuals`
> are i.i.d. unit-rate exponentials. Returns ``(D, p_value)`` from a one-sample KS
> test of the residuals against the ``Exp(1)`` CDF; a small p-value rejects the
> fitted model. Requires at least three events.

### `hawkes_intensity(t, history, mu, alpha, beta)`  _function_

> Conditional intensity ``lambda(t)`` given past event times ``history``.
>
> Sums the exponential kernel over events strictly before ``t``.

### `hawkes_log_likelihood(events, mu, alpha, beta, t_end=None)`  _function_

> Exact log-likelihood of an exponential-kernel Hawkes process.
>
> Uses the standard recursion for the excitation term ``A_i = sum_{j<i}
> exp(-beta (t_i - t_j))`` (updated as ``A_i = exp(-beta dt)(1 + A_{i-1})``), so it
> is linear in the number of events. ``events`` is a sorted list of event times on
> ``[0, t_end]`` (``t_end`` defaults to the last event). The log-likelihood is
>
>     sum_i log(mu + alpha A_i) - mu T - (alpha/beta) sum_i (1 - exp(-beta (T-t_i))).

### `hawkes_residuals(events, mu, alpha, beta)`  _function_

> Time-rescaling residuals of a fitted Hawkes process.
>
> By the time-rescaling theorem the integrated intensity (compensator) between
> consecutive events, ``tau_i = integral_{t_{i-1}}^{t_i} lambda(u) du``, is a
> sequence of i.i.d. unit-rate exponentials when the model is correctly specified.
> Computes those ``tau_i`` for ``i >= 1`` using the exponential-kernel recursion,
>
>     tau_i = mu (t_i - t_{i-1})
>             + (alpha/beta) sum_{j<i} [exp(-beta (t_{i-1}-t_j)) - exp(-beta (t_i-t_j))],
>
> accumulated via a running kernel sum. Returns the list of ``n - 1`` residuals.

### `hawkes_simulate(mu, alpha, beta, t_max, seed=1234567)`  _function_

> Simulate a Hawkes process on ``[0, t_max]`` by Ogata's thinning algorithm.
>
> Returns the sorted list of event times. Requires ``alpha < beta`` (stationarity)
> for a well-behaved simulation. Deterministic given ``seed`` (an LCG uniform
> stream).

## heat_equation

### `heat_equation_cn(u0, alpha, dx, dt, n_steps, left=None, right=None)`  _function_

> Crank-Nicolson evolution of ``u_t = alpha u_xx``.
>
> ``u0`` is the initial profile (interior + boundary grid values). ``dx``/``dt`` the
> space/time steps, ``alpha`` the diffusivity. ``left``/``right`` fix the Dirichlet
> boundary values (default: hold the initial endpoints). Returns the profile after
> ``n_steps`` steps. Unconditionally stable, second-order in space and time.

## heavy_hitters

### `MisraGries(k)`  _class_

> Misra-Gries frequent-items summary with ``k`` counters.
>
> ``add(item)`` folds one occurrence; ``counts()`` returns the surviving item -> count
> map. Any item with true frequency greater than ``n / (k + 1)`` is guaranteed present,
> and stored counts never exceed the true count (they may undercount by up to the number
> of decrement rounds).

### `SpaceSaving(k)`  _class_

> Space-Saving top-k summary (Metwally-Agarwal-Abbadi) with ``k`` counters.
>
> ``add(item)`` folds one occurrence; ``top(m)`` returns the ``m`` highest
> ``(item, count)`` pairs. Each slot tracks an ``error`` -- the maximum it may overcount
> -- so ``count - error`` is a guaranteed lower bound on the true frequency. When a new
> item arrives and every slot is used, the slot with the smallest count is evicted and
> its count becomes the new item's starting error.

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

## hilbert

### `analytic_signal(x)`  _function_

> Analytic signal of a real sequence ``x`` (length a power of two).
>
> Returns a list of complex values ``z`` with ``z.real == x`` (to floating error)
> and ``z.imag`` the Hilbert transform of ``x``. The one-sided spectrum is formed by
> keeping the DC and Nyquist bins, doubling the positive-frequency bins, and zeroing
> the negative-frequency bins, then inverting the FFT.

### `envelope(x)`  _function_

> Instantaneous-amplitude envelope ``|analytic_signal(x)|``.
>
> For an amplitude-modulated carrier this recovers the modulating envelope; for a
> pure sinusoid it is a constant equal to the amplitude. Length a power of two.

### `hilbert_transform(x)`  _function_

> Discrete Hilbert transform of ``x``: the imaginary part of its analytic signal.
>
> A 90-degree phase shift of every frequency component (cosine -> sine). Length must
> be a power of two.

### `instantaneous_frequency(x, dt=1.0)`  _function_

> Instantaneous frequency (cycles per unit time) from the analytic signal.
>
> The derivative of the unwrapped phase, ``d(phase)/dt / (2*pi)``, by central
> differences on the interior and one-sided differences at the ends. For a pure tone
> this is flat at the tone's frequency. ``dt`` is the sample spacing; length a power
> of two.

### `instantaneous_phase(x)`  _function_

> Unwrapped instantaneous phase (radians) of the analytic signal of ``x``.
>
> ``atan2(imag, real)`` per sample, unwrapped so it is continuous rather than
> jumping by ``2*pi``. Length a power of two.

## hmm

### `hmm_baum_welch(obs, n_states, n_symbols, max_iter=100, tol=1e-06, seed=1234567)`  _function_

> Baum-Welch (EM) estimate of HMM parameters from an observation sequence.
>
> Iterates the forward-backward E-step and the re-estimation M-step from a
> near-uniform (seed-jittered) start, returning a dict with ``pi``, ``A``, ``B``,
> ``log_likelihood`` and ``n_iter``. The log-likelihood is non-decreasing; the
> labelling of states is arbitrary (identifiable only up to a permutation). Pure
> standard library.

### `hmm_forward(pi, A, B, obs)`  _function_

> Log-likelihood ``log P(obs | model)`` by the scaled forward algorithm.
>
> Rescales the forward variables at each step (dividing by their sum) and
> accumulates the log of the scale factors, so the likelihood is exact without
> underflow on long sequences. ``pi`` initial distribution, ``A`` transition, ``B``
> emission matrices; ``obs`` a list of symbol indices.

### `hmm_posterior(pi, A, B, obs)`  _function_

> Smoothed posterior state probabilities ``P(state_t = i | obs)`` (gamma).
>
> Runs the scaled forward-backward algorithm and returns a list of length-``n``
> distributions, one per time step, each summing to one -- the probability of being
> in each hidden state at that time given the *entire* observation sequence (unlike
> Viterbi's single best path, this is the per-time marginal). Pure standard library.

### `hmm_simulate(pi, A, B, length, seed=1234567)`  _function_

> Simulate ``length`` observations from an HMM. Returns ``(states, obs)``.

### `hmm_viterbi(pi, A, B, obs)`  _function_

> Most-likely hidden-state path (Viterbi) and its log-probability.
>
> Dynamic programming in log space; returns ``(path, log_prob)`` where ``path`` is
> the list of state indices maximizing the joint probability of states and
> observations. Ties are broken toward the lower state index.

## hodges_lehmann

### `hodges_lehmann_location(x)`  _function_

> One-sample Hodges-Lehmann estimator: median of the Walsh averages.
>
> Median of ``(x_i + x_j) / 2`` over all ``i <= j`` (including ``i == j``). Robust
> (29% breakdown) and highly efficient at the normal; estimates the center of a
> symmetric distribution.

### `hodges_lehmann_shift(x, y)`  _function_

> Two-sample Hodges-Lehmann shift: median of all pairwise differences.
>
> Median of ``y_j - x_i`` over every pair. The robust estimate of the location
> shift between the two samples, consistent with the Wilcoxon rank-sum test (the
> shift for which the test would not reject). Positive means ``y`` is shifted above
> ``x``.

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

## holt_winters

### `holt_linear(y, alpha, beta, horizon=1)`  _function_

> Holt's linear-trend (double) exponential smoothing.
>
> Parameters
> ----------
> y : sequence of float
>     The series (length >= 2).
> alpha, beta : float
>     Level and trend smoothing parameters in [0, 1].
> horizon : int
>     Forecast horizon (number of steps beyond the last observation).
>
> Returns
> -------
> (level, trend, forecast) : (float, float, list[float])
>     Final level and trend, and the ``horizon``-step-ahead forecasts
>     ``l + h*b`` (a straight line in ``h``).

### `holt_winters_add(y, alpha, beta, gamma, period, horizon=1)`  _function_

> Additive Holt-Winters (triple) exponential smoothing.
>
> Parameters
> ----------
> y : sequence of float
>     The series; length must exceed ``2 * period``.
> alpha, beta, gamma : float
>     Level, trend, and seasonal smoothing parameters in [0, 1].
> period : int
>     Season length ``m`` (>= 2).
> horizon : int
>     Forecast horizon.
>
> Returns
> -------
> (level, trend, seasonals, forecast) : (float, float, list[float], list[float])
>     Final level, trend, the last ``period`` seasonal factors, and the
>     ``horizon``-step forecasts ``l + h*b + s[(h-1) mod m]``.

## hp_filter

### `hp_filter(y, lam=1600.0)`  _function_

> Hodrick-Prescott trend/cycle decomposition.
>
> Parameters
> ----------
> y : sequence of float
>     The time series.
> lam : float
>     Smoothing parameter (>= 0). Larger values yield a smoother trend; 1600 is
>     the standard quarterly value.
>
> Returns
> -------
> (trend, cycle) : (list[float], list[float])
>     The smooth trend and the cyclical residual ``cycle = y - trend``. Their
>     sum reconstructs ``y`` exactly.

## hrp

### `hierarchical_risk_parity(cov)`  _function_

> Hierarchical Risk Parity weights (López de Prado).
>
> Clusters assets by correlation distance, quasi-diagonalizes, and allocates by
> recursive bisection: each split gives the two sub-clusters weights inversely
> proportional to their inverse-variance-weighted variances. Returns positive
> weights summing to one, aligned with the original asset order. Needs no matrix
> inversion, so it is stable for near-singular covariances.

### `inverse_volatility_weights(cov)`  _function_

> Inverse-volatility weights ``(1/sigma_i) / sum_j (1/sigma_j)``.
>
> Each asset weighted by the reciprocal of its standard deviation, normalized to
> sum to one. Higher-volatility assets get less capital; ignores correlations.

## huber_regression

### `huber_regression(X, y, delta=1.345, add_intercept=True, max_iter=50, tol=1e-08)`  _function_

> Huber robust regression by IRLS.
>
> ``delta`` is the residual threshold (in robust-scale units) between the quadratic
> and linear regions; the default 1.345 gives ~95% efficiency at the normal. Returns
> a dict with ``coefficients``, ``n_iter`` and ``scale`` (the final robust residual
> scale). Reduces toward OLS as ``delta`` grows.

## hull_white

### `hw_B(a, tau)`  _function_

> Hull-White ``B(t, T) = (1 - e^{-a*tau}) / a`` for ``tau = T - t``.

### `hw_bond_option(P0, a, sigma, t_option, t_bond, strike, is_call=True)`  _function_

> European option on a zero-coupon bond under Hull-White (analytic).
>
> Prices an option expiring at ``t_option`` on a zero maturing at ``t_bond``,
> struck at ``strike``, using the initial discount curve ``P0``. The forward
> bond price is log-normal with volatility
>
>     sigma_P = sigma * B(a, t_bond - t_option) * sqrt((1 - e^{-2 a t_option}) / (2 a)),
>
> giving a Black-style formula in the discount factors ``P0(t_bond)`` and
> ``P0(t_option)`` (Jamshidian / Hull-White). Returns the option value today.
> Requires ``0 < t_option < t_bond``.

### `hw_cap(P0, a, sigma, dates, strike, notional=1.0)`  _function_

> Hull-White cap: sum of caplets over consecutive ``dates`` (reset, pay pairs).
>
> ``dates`` is the schedule ``[t_0, t_1, ..., t_n]``; caplet ``i`` covers
> ``[t_i, t_{i+1}]``. A floor is the analogous sum of floorlets. By put-call
> parity ``cap - floor`` equals the value of the fixed-vs-float swap.

### `hw_caplet(P0, a, sigma, reset, pay, strike, notional=1.0)`  _function_

> Hull-White caplet: an option on the simple forward rate over ``[reset, pay]``.
>
> A caplet paying ``notional * tau * max(L - strike, 0)`` at ``pay`` (where ``L``
> is the simple rate set at ``reset`` for accrual ``tau = pay - reset``) equals
> ``notional * (1 + strike*tau)`` puts on the ``pay``-zero struck at
> ``1/(1 + strike*tau)``, expiring at ``reset`` -- the standard bond-option
> representation. Priced analytically off the initial curve.

### `hw_floor(P0, a, sigma, dates, strike, notional=1.0)`  _function_

> Hull-White floor: sum of floorlets over consecutive ``dates``.

### `hw_floorlet(P0, a, sigma, reset, pay, strike, notional=1.0)`  _function_

> Hull-White floorlet: ``notional * (1 + strike*tau)`` calls on the pay-zero.

### `hw_swaption(P0, r0, a, sigma, expiry, pay_times, fixed_rate, is_payer=True, notional=1.0)`  _function_

> European swaption under Hull-White via Jamshidian decomposition.
>
> At ``expiry`` the holder enters a swap paying (payer) or receiving (receiver)
> ``fixed_rate`` on ``pay_times``. Jamshidian's trick: the underlying coupon bond
> is monotone in the short rate, so find the critical rate ``r*`` at which its
> value equals par, then the swaption is a portfolio of options on each
> zero-coupon cashflow struck at that cashflow's value at ``r*``. Priced off the
> initial curve; ``r0`` is only used to seed the ``r*`` search.
>
> A payer swaption is a put on the coupon bond (a portfolio of zero puts); a
> receiver is the corresponding call portfolio.

### `hw_zero_from_curve(P0, r0, a, sigma, t, T, f0=None, eps=1e-05)`  _function_

> Hull-White bond price ``P(t, T)`` fitted to an initial curve ``P0``.
>
> Parameters
> ----------
> P0 : callable
>     Initial discount factor ``P0(T)`` observed today (``P0(0) = 1``).
> r0 : float
>     Current short rate; for consistency it should equal the initial instant
>     forward ``f(0,0)``.
> a, sigma : float
>     Mean reversion and volatility (``a`` may be 0 for the Ho-Lee limit).
> t, T : float
>     Valuation and maturity times, ``0 <= t <= T``.
> f0 : callable, optional
>     Initial instantaneous forward ``f(0, t)``; defaults to a finite-difference
>     of ``-ln P0``.
> eps : float
>     Step for the forward finite difference.
>
> Returns
> -------
> float
>     The fitted ``P(t, T)``. At ``t = 0`` it reproduces ``P0(T)`` exactly.

## hurst

### `hurst_exponent(x, min_window=8, max_window=None)`  _function_

> Estimate the Hurst exponent of ``x`` by rescaled-range analysis.
>
> Parameters
> ----------
> x : sequence of float
>     The series (levels for a walk, or increments -- interpret ``H`` relative
>     to what you feed in).
> min_window : int
>     Smallest window length used in the log-log regression (>= 2).
> max_window : int, optional
>     Largest window length. Defaults to ``len(x) // 2``.
>
> Returns
> -------
> float
>     The Hurst exponent, the slope of ``log(R/S)`` versus ``log(n)`` over
>     dyadic window sizes. ~0.5 for a memoryless series, >0.5 persistent,
>     <0.5 mean-reverting.

### `rescaled_range(window)`  _function_

> Rescaled range ``R/S`` of a single window.
>
> ``R`` is the range of the cumulative demeaned series; ``S`` is the window's
> population standard deviation. Returns 0 when the window is constant.

## hurwitz

### `hurwitz_zeta(s, a, terms=20, corrections=10)`  _function_

> Hurwitz zeta ``zeta(s, a) = sum_{n>=0} (n + a)^-s`` for real ``s > 1``, ``a > 0``.
>
> Euler-Maclaurin summation: the first ``terms`` shifted terms are summed exactly, the tail
> is replaced by its integral, and ``corrections`` Bernoulli terms recover the remainder.
> With ``a = 1`` this reduces to the Riemann zeta.

### `polygamma(m, x)`  _function_

> Polygamma ``psi^(m)(x)``, the ``m``-th derivative of the digamma, for ``x > 0``.
>
> ``m = 0`` returns the digamma :func:`quantforge.special.digamma`. For ``m >= 1`` uses
> ``psi^(m)(x) = (-1)^(m+1) m! * hurwitz_zeta(m+1, x)``. ``polygamma(1, .)`` is the trigamma.

## hyperdual

### `HyperDual(f0, f1=0.0, f2=0.0, f12=0.0)`  _class_

> Hyperdual number ``f0 + f1 e1 + f2 e2 + f12 e1e2`` for 1st/2nd derivatives.
>
> Seed a variable as ``HyperDual(x, 1.0, 1.0, 0.0)``; arithmetic and this module's
> elementary functions propagate all four parts exactly. After evaluating ``f``, the
> result's ``f1`` (or ``f2``) is ``f'(x)`` and ``f12`` is ``f''(x)``.

### `hyperdual_derivatives(f, x)`  _function_

> Return ``(f(x), f'(x), f''(x))`` exactly by one hyperdual evaluation.
>
> ``f`` must accept a :class:`HyperDual` built from hyperdual arithmetic and this
> module's elementary functions. No finite-difference error in either derivative.

### `second_derivative(f, x)`  _function_

> Exact ``f''(x)`` by hyperdual autodiff.

## hyperloglog

### `HyperLogLog(p=14)`  _class_

> HyperLogLog distinct-count estimator with ``2^p`` registers.
>
> ``p`` in ``[4, 16]`` trades memory for accuracy: ``m = 2^p`` registers give a
> relative error near ``1.04 / sqrt(m)`` (``p = 14`` -> ~0.8%). Add items with
> :meth:`add` (any hashable stringifiable value); read the estimate from
> :meth:`count`. Two sketches with the same ``p`` merge via :meth:`merge` (union
> cardinality), the property that makes it distributable.

## hypothesis

### `binomial_test(k, n, prob=0.5, alternative='two-sided')`  _function_

> Exact binomial test that the success probability equals ``prob``.
>
> ``alternative`` is ``"greater"`` (``P(X >= k)``), ``"less"`` (``P(X <= k)``), or
> ``"two-sided"`` (sum of all outcome probabilities no larger than the observed
> one). Returns ``(proportion, p_value)`` where ``proportion = k / n``.

### `chi_square_gof_test(observed, expected=None)`  _function_

> Pearson chi-square goodness-of-fit test.
>
> Compares observed counts against ``expected`` (defaults to a uniform
> distribution over the categories). The statistic ``sum (O - E)^2 / E`` is
> referenced to a chi-square with ``k - 1`` degrees of freedom. Returns
> ``(statistic, p_value)``; a small p-value rejects the fit.

### `chi_square_independence_test(table)`  _function_

> Chi-square test of independence for a contingency ``table`` (rows x cols).
>
> Uses the row/column marginals to form the expected counts under independence
> and references ``sum (O - E)^2 / E`` to a chi-square with
> ``(rows - 1)(cols - 1)`` degrees of freedom. Returns ``(statistic, p_value)``.

### `mann_whitney_u(a, b)`  _function_

> Mann-Whitney U rank-sum test (two-sided, normal approximation with ties).
>
> Ranks the pooled samples (average ranks for ties) and forms the smaller of the
> two U statistics; the p-value uses the normal approximation with a tie
> correction to the variance and a continuity correction. Returns
> ``(u, p_value)``, where ``u`` is ``min(U_a, U_b)``. A distribution-free
> alternative to the two-sample t when normality is doubtful.

### `one_sample_t_test(sample, mu0=0.0)`  _function_

> One-sample two-sided Student-t test that the mean equals ``mu0``.
>
> ``t = (xbar - mu0) / (s / sqrt(n))`` on ``n - 1`` degrees of freedom. Returns
> ``(t, p_value)``.

### `one_way_anova(*groups)`  _function_

> One-way ANOVA F-test across two or more samples.
>
> Partitions the total variation into between-group and within-group sums of
> squares and forms ``F = MS_between / MS_within``, referenced to an F with
> ``(k - 1, N - k)`` degrees of freedom. Returns ``(F, p_value)``; a small
> p-value rejects equality of the group means.

### `paired_t_test(a, b)`  _function_

> Paired (dependent) two-sided Student-t test on the within-pair differences.
>
> Equivalent to a one-sample t-test of ``a[i] - b[i]`` against zero, on ``n - 1``
> degrees of freedom. Returns ``(t, p_value)``.

### `two_sample_t_test(a, b, equal_var=True)`  _function_

> Two-sample Student-t test of equal means, two-sided.
>
> With ``equal_var=True`` uses the pooled-variance t-test (``n_a + n_b - 2``
> degrees of freedom); with ``equal_var=False`` uses Welch's t-test with the
> Welch-Satterthwaite degrees of freedom. Returns ``(t, p_value)``.

## icc

### `icc(data)`  _function_

> Shrout-Fleiss intraclass correlation coefficients.
>
> ``data`` is a list of rows, one per subject, each a length-``k`` list of the
> ``k`` raters' scores (a balanced subjects x raters table). Returns a dict with
> ``icc1``, ``icc2_1``, ``icc2_k``, ``icc3_1`` and ``icc3_k``. Requires at least two
> subjects and two raters.

## iir_filter

### `butter_highpass(order, cutoff)`  _function_

> Design a Butterworth high-pass as a cascade of biquad sections.
>
> ``cutoff`` is the -3 dB frequency in cycles/sample, in ``(0, 0.5)``. Returns a list
> of ``(b, a)`` second-order sections with the Nyquist gain normalized to 1.

### `butter_lowpass(order, cutoff)`  _function_

> Design a Butterworth low-pass as a cascade of biquad sections.
>
> ``cutoff`` is the -3 dB frequency in cycles/sample, in ``(0, 0.5)``. Returns a list
> of ``(b, a)`` second-order sections (each a length-3 numerator and denominator with
> ``a[0] == 1``), with the overall DC gain normalized to 1.

### `iir_frequency_response(sections, freqs)`  _function_

> Magnitude response ``|H(f)|`` of a biquad cascade at normalized ``freqs``.
>
> ``freqs`` in cycles/sample (``0`` to ``0.5``). Returns the product of the section
> magnitudes at each frequency -- useful for verifying a design's passband/stopband.

### `sosfilt(sections, x)`  _function_

> Filter ``x`` through a cascade of biquad ``(b, a)`` sections (Direct Form II transposed).
>
> Runs the signal through each second-order section in turn. Returns a list the same
> length as ``x``.

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

## info_criteria

### `aic(log_likelihood, k)`  _function_

> Akaike information criterion ``-2 logL + 2k`` (``k`` = number of parameters).

### `aicc(log_likelihood, k, n)`  _function_

> Small-sample corrected AIC ``AIC + 2k(k+1)/(n-k-1)``.
>
> Approaches :func:`aic` as ``n`` grows; use it when ``n`` is not much larger than
> ``k``. Requires ``n > k + 1``.

### `bic(log_likelihood, k, n)`  _function_

> Bayesian (Schwarz) information criterion ``-2 logL + k ln n``.
>
> Penalizes complexity more than AIC for ``n >= 8`` (``ln n > 2``), so it selects
> more parsimonious models.

### `gaussian_log_likelihood(rss, n)`  _function_

> Concentrated Gaussian log-likelihood from a residual sum of squares.
>
> ``logL = -n/2 (ln(2 pi) + ln(rss/n) + 1)`` at the MLE noise variance
> ``sigma^2 = rss/n``. Increases (toward zero from below) as the fit improves.

### `hqic(log_likelihood, k, n)`  _function_

> Hannan-Quinn information criterion ``-2 logL + 2k ln(ln n)``.
>
> A penalty between AIC and BIC for moderate ``n``. Requires ``n >= 3`` so that
> ``ln(ln n)`` is defined and positive.

## installment

### `installment_call(S, K, t, r, sigma, installment, pay_times, steps=200, q=0.0)`  _function_

> Fair upfront value of a European installment call on a CRR tree.
>
> Parameters
> ----------
> S, K, t, r, sigma : float
>     Spot, strike, time to expiry (years), risk-free rate, volatility.
> installment : float
>     Amount paid at each date in ``pay_times`` to keep the option alive.
>     Must be non-negative.
> pay_times : sequence of float
>     Installment payment dates in years, each in ``(0, t)``. Duplicates and
>     out-of-range values are rejected.
> steps : int
>     Number of tree steps.
> q : float
>     Continuous dividend yield (carry ``b = r - q``).
>
> Returns
> -------
> float
>     Fair upfront premium at inception. At least zero; equals the plain
>     European call value when ``installment`` is zero or ``pay_times`` empty.

## integrate2d

### `integrate2d_gauss(f, a, b, c, d, n=5)`  _function_

> Tensor Gauss-Legendre integral of ``f(x, y)`` over ``[a, b] x [c, d]``.
>
> ``n`` in {2, 3, 4, 5} is the order per axis; exact for polynomials up to degree
> ``2n - 1`` in each variable. Uses ``n^2`` function evaluations.

### `integrate2d_simpson(f, a, b, c, d, nx=50, ny=50)`  _function_

> Composite-Simpson double integral of ``f(x, y)`` over ``[a, b] x [c, d]``.
>
> Integrates in ``y`` for each ``x`` then in ``x`` (Fubini), each by composite
> Simpson with ``ny`` / ``nx`` panels. Slower but robust for integrands the smooth
> Gauss rule would miss.

## interpolation

### `monotone_cubic(xs, ys)`  _function_

> Fritsch-Carlson monotone cubic Hermite interpolant, returning ``f(x)``.
>
> Chooses the Hermite tangents so the interpolant preserves the monotonicity of
> the data: no overshoot between knots. Exact at the knots. The standard choice
> for discount-factor / survival curves that must not wiggle below/above the data.

### `natural_cubic_spline(xs, ys)`  _function_

> Build a natural cubic spline interpolant, returning a callable ``f(x)``.
>
> Solves the tridiagonal system for the second derivatives with zero-curvature
> (natural) end conditions. The returned function evaluates the piecewise cubic
> and is exact at the knots, C2 in between. Clamps to the end segments outside
> ``[xs[0], xs[-1]]``.

## interval

### `Interval(lo, hi=None)`  _class_

> A closed interval ``[lo, hi]`` with guaranteed-enclosure arithmetic.
>
> Construct from two endpoints (order-normalized) or a single number (a degenerate
> point interval). Supports ``+ - * /`` with intervals and scalars, ``width``,
> ``midpoint``, ``contains``, ``intersect``, and monotone functions ``exp``/``log``/
> ``sqrt``/``__pow__`` (integer powers).

## interval_set

### `intervals_intersection(a, b)`  _function_

> Intersection of two interval collections: the ranges covered by *both*.
>
> Merges each side first, then sweeps for overlaps. Returns a sorted
> non-overlapping list.

### `intervals_union(a, b)`  _function_

> Union of two interval collections as a merged non-overlapping list.

### `max_overlap(intervals)`  _function_

> Maximum number of intervals overlapping at any single point.
>
> Sweeps the endpoints (starts before ends at the same coordinate, so touching
> intervals count as overlapping). Returns the peak count; ``0`` for no intervals.

### `merge_intervals(intervals)`  _function_

> Merge overlapping/adjacent intervals into a minimal sorted non-overlapping cover.
>
> ``[(1,3),(2,6),(8,10)] -> [(1,6),(8,10)]``. Touching intervals (``(1,2),(2,3)``) merge
> into ``(1,3)``. Returns a new sorted list; empty input yields ``[]``.

### `total_covered_length(intervals)`  _function_

> Total length covered by the union of the intervals (overlaps counted once).

## inversions

### `count_inversions(values)`  _function_

> Number of inversions in ``values`` (pairs ``i < j`` with ``values[i] > values[j]``).
>
> Uses a merge sort, so ``O(n log n)`` rather than the ``O(n^2)`` brute count. A strictly
> increasing sequence has 0; a strictly decreasing one has ``n(n-1)/2``.

### `is_sorted(values, strict=False)`  _function_

> True if ``values`` is non-decreasing (or strictly increasing when ``strict``).

### `kendall_tau_distance(rank_a, rank_b)`  _function_

> Kendall-tau distance: number of pairs ordered oppositely in two rankings.
>
> ``rank_a`` and ``rank_b`` are sequences of the same items (equal length, same set). The
> distance is the count of pairs ``(x, y)`` whose relative order differs between the two,
> i.e. the inversions of ``rank_b`` reindexed by ``rank_a``'s positions.

## isotonic

### `isotonic_fit(x, y, weights=None, increasing=True)`  _function_

> Fit an isotonic step function of ``x`` and return a predictor callable.
>
> Sorts by ``x``, runs :func:`isotonic_regression` on the reordered response, and
> returns ``(y_hat, predict)`` where ``y_hat`` is the fit aligned to the *original*
> input order and ``predict(x_new)`` interpolates the monotone step fit at a new
> point (linear between fitted knots, clamped to the endpoints). Ties in ``x`` share
> a fitted value.

### `isotonic_regression(y, weights=None, increasing=True)`  _function_

> Weighted isotonic regression by pool-adjacent-violators.
>
> Returns the fitted values ``y_hat`` (same length as ``y``) forming the monotone
> sequence closest to ``y`` in weighted least squares. ``weights`` defaults to all
> ones; ``increasing=False`` fits a non-increasing sequence. The fit is a step
> function: tied blocks share their common weighted mean.

## jacobi_elliptic

### `jacobi_am(u, m, tol=1e-15, max_iter=64)`  _function_

> Jacobi amplitude ``am(u, m)`` -- the angle ``phi`` with ``u = F(phi | m)``.
>
> AGM descent: build the sequence ``a_0 = 1, b_0 = sqrt(1-m), c_0 = sqrt(m)`` down to
> ``c_n ~ 0``, then climb back, halving the accumulated angle. ``m`` is the parameter
> ``k^2`` in ``[0, 1]``. Reduces to ``am(u, 0) = u``.

### `jacobi_cn(u, m)`  _function_

> Jacobi elliptic ``cn(u, m) = cos(am(u, m))``. Reduces to ``cos u`` at ``m = 0``.

### `jacobi_dn(u, m)`  _function_

> Jacobi elliptic ``dn(u, m) = sqrt(1 - m sin^2 am(u, m))``. Reduces to ``1`` at ``m = 0``.

### `jacobi_sn(u, m)`  _function_

> Jacobi elliptic ``sn(u, m) = sin(am(u, m))``. Reduces to ``sin u`` at ``m = 0``.

## jonckheere

### `jonckheere_terpstra_test(groups)`  _function_

> Jonckheere-Terpstra trend test across ordered ``groups``.
>
> ``groups`` is a sequence of samples given in the hypothesized order (e.g. from
> lowest dose to highest). Returns a dict with the ``statistic`` J, its null ``mean``
> and ``variance`` (tie-corrected), the ``z`` normal approximation and the two-sided
> ``p_value``. A positive ``z`` indicates an increasing trend across the group order.

## jump_test

### `bns_jump_test(returns)`  _function_

> Barndorff-Nielsen-Shephard ratio jump-test statistic and p-value.
>
> Returns ``(z, p_value)`` where ``z`` is the ratio statistic (asymptotically
> ``N(0, 1)`` under no jump) and ``p_value`` is the upper-tail ``P(Z > z)``. A
> small p-value rejects the no-jump null in favour of a jump having occurred
> during the sampled period. Requires at least three returns.

### `tripower_quarticity(returns)`  _function_

> Tripower quarticity ``TQ``, a jump-robust estimator of integrated quarticity.
>
> ``TQ = n * mu_{4/3}^{-3} * sum |r_{i-2}|^{4/3} |r_{i-1}|^{4/3} |r_i|^{4/3}`` with
> ``mu_{4/3} = 2^{2/3} Gamma(7/6) / Gamma(1/2)``. Robust to jumps (each enters only
> one triple), it estimates ``integral sigma^4`` and sets the scale of the jump
> test. Requires at least three returns.

## kabsch

### `kabsch(P, Q)`  _function_

> Optimal rotation/translation mapping points ``P`` onto ``Q`` (least-squares).
>
> ``P`` and ``Q`` are equal-length lists of 3-D points (paired). Returns
> ``(R, t, rmsd)``: the 3x3 rotation matrix and translation ``t`` such that
> ``R @ (p - centroid_P) + centroid_Q`` best matches ``q``, and the residual RMSD.

## kalman

### `kalman_local_level(observations, process_var, obs_var, x0=None, p0=None)`  _function_

> Run a scalar local-level Kalman filter over ``observations``.
>
> Parameters
> ----------
> observations : sequence of float
>     The measured series ``y_t``.
> process_var : float
>     State (process) noise variance ``Q`` (>= 0).
> obs_var : float
>     Observation noise variance ``R`` (> 0).
> x0 : float, optional
>     Prior mean for the level. Defaults to the first observation.
> p0 : float, optional
>     Prior variance for the level. Defaults to ``obs_var`` (a diffuse-ish
>     start). Larger values weight the data more at the outset.
>
> Returns
> -------
> (levels, variances, gains) : (list[float], list[float], list[float])
>     The filtered level estimate, its posterior variance, and the Kalman gain
>     at each step.

### `kalman_steady_state_gain(process_var, obs_var)`  _function_

> Steady-state Kalman gain for the local-level model.
>
> Solves the algebraic Riccati fixed point in closed form; depends only on the
> signal-to-noise ratio ``q = process_var / obs_var``. Returns a gain in
> ``(0, 1)`` that rises toward 1 as the process noise dominates and toward 0 as
> the observation noise dominates.

## kalman_beta

### `kalman_regression_beta(x, y, process_var, obs_var, beta0=0.0, p0=1000000.0)`  _function_

> Filter a time-varying regression slope ``beta_t`` of ``y`` on ``x``.
>
> Parameters
> ----------
> x, y : sequence of float
>     Regressor and response of equal length.
> process_var : float
>     Slope random-walk variance ``Q`` (>= 0). ``0`` gives a static slope
>     (recursive least squares); larger values let it drift faster.
> obs_var : float
>     Observation noise variance ``R`` (> 0).
> beta0 : float
>     Prior mean of the slope.
> p0 : float
>     Prior variance of the slope. A large default gives a near-diffuse start so
>     early data dominates.
>
> Returns
> -------
> (betas, variances) : (list[float], list[float])
>     The filtered slope and its posterior variance at each step.

## kalman_filter

### `kalman_filter(observations, F, H, Q, R, x0, P0)`  _function_

> Kalman filter over a sequence of observation vectors.
>
> ``F`` (state x state), ``H`` (obs x state), ``Q`` (state x state process cov),
> ``R`` (obs x obs measurement cov), ``x0``/``P0`` the initial state mean and
> covariance. Each entry of ``observations`` is an observation vector. Returns a
> dict with ``states`` (filtered means), ``covariances``, and ``log_likelihood``
> (the Gaussian log-likelihood of the observations from the prediction errors).

### `kalman_smoother(observations, F, H, Q, R, x0, P0)`  _function_

> Rauch-Tung-Striebel smoother: filter then a backward refinement pass.
>
> Returns a dict with ``states`` (smoothed means) and ``covariances``. Each smoothed
> covariance is no larger (in the positive-definite sense; here checked on the
> trace) than the corresponding filtered covariance, since the smoother conditions
> on the whole series rather than only the past.

## kde

### `kde(data, x, bandwidth=None, rule='silverman')`  _function_

> Evaluate the Gaussian KDE of ``data`` at point(s) ``x``.
>
> ``bandwidth`` overrides the rule (``"silverman"`` or ``"scott"``). ``x`` may be a
> scalar (returns a float) or an iterable (returns a list). The estimate is
> non-negative everywhere and integrates to one.

### `kde_function(data, bandwidth=None, rule='silverman')`  _function_

> Return a callable density estimator ``f(x)`` for ``data`` (bandwidth fixed once).

### `scott_bandwidth(data)`  _function_

> Scott's rule-of-thumb bandwidth ``h = std * n^{-1/5}``.

### `silverman_bandwidth(data)`  _function_

> Silverman's rule-of-thumb bandwidth.
>
> ``h = 0.9 * min(std, IQR/1.34) * n^{-1/5}`` -- robust to mild non-normality via the
> IQR term. The standard default for a unimodal, roughly-Gaussian sample.

## kdtree

### `KDTree(points)`  _class_

> Static k-d tree over a fixed set of points for spatial queries.
>
> Build once from a list of equal-length coordinate tuples; then run `nearest`,
> `k_nearest`, `within_radius`, and `range_search`. Query methods return the original
> point *indices* (in build order), so callers can map back to their own payloads.

## kendall_test

### `goodman_kruskal_gamma(x, y)`  _function_

> Goodman-Kruskal gamma: ``(C - D) / (C + D)``, ignoring all tied pairs.

### `kendall_tau_b(x, y)`  _function_

> Tie-corrected Kendall's tau-b.
>
> ``(C - D) / sqrt((C + D + Tx)(C + D + Ty))`` where ``C``/``D`` are concordant/
> discordant pairs and ``Tx``/``Ty`` are pairs tied only in ``x`` / only in ``y``.
> Reaches +/-1 for a perfect monotone relationship even when ties are present.

### `kendall_tau_test(x, y)`  _function_

> Test ``tau = 0`` by the large-sample normal approximation.
>
> Returns a dict with ``tau_b``, the ``z`` statistic and the two-sided ``p_value``.
> The variance of the concordance statistic ``S = C - D`` under independence is
> ``n(n-1)(2n+5)/18`` (no tie correction to the variance -- adequate away from heavy
> ties); ``z = S / sqrt(Var S)``.

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

## kmeans

### `kmeans(X, k, max_iter=100, seed=1234567)`  _function_

> Cluster ``X`` into ``k`` groups with Lloyd's algorithm (k-means++ seed).
>
> Parameters
> ----------
> X : list[list[float]]
>     ``n`` points in ``d`` dimensions.
> k : int
>     Number of clusters (>= 1, <= n).
> max_iter : int
>     Maximum Lloyd iterations.
> seed : int
>     Seed for the deterministic k-means++ sampler.
>
> Returns
> -------
> dict
>     ``labels`` (cluster index per point), ``centroids``, ``inertia`` (total
>     within-cluster squared distance), and ``iterations``.

## knn

### `knn_classify(X_train, y_train, X_query, k=3)`  _function_

> Classify each query row by majority vote of its ``k`` nearest neighbors.
>
> Ties are broken toward the label that is closest on average (the first
> encountered at minimum total distance). Returns a predicted label per query.

### `knn_regress(X_train, y_train, X_query, k=3)`  _function_

> Predict each query as the mean target of its ``k`` nearest neighbors.

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

## kpss

### `kpss_test(y, lags=None, regression='c')`  _function_

> KPSS stationarity test statistic and approximate p-value.
>
> ``regression`` is ``"c"`` for level stationarity (residuals about the mean) or
> ``"ct"`` for trend stationarity (residuals about a fitted linear trend).
> ``lags`` sets the Newey-West bandwidth for the long-run variance (defaults to
> ``floor(4 (n/100)^{1/4})``). Returns ``(eta, p_value)``; a *small* p-value
> rejects stationarity (unlike ADF, where a small p-value supports it). The p-value
> is interpolated/clamped against the asymptotic critical values, so it is reported
> within ``[0.01, 0.10]`` at the bounds.

## kruskal_wallis

### `friedman_test(blocks)`  _function_

> Friedman test for ``k`` related treatments over ``b`` blocks.
>
> ``blocks`` is a sequence of rows, each a length-``k`` sequence giving one block's
> measurements across the treatments (e.g. one subject rated under every condition).
> Ranks within each block, then compares treatment rank sums. Returns a dict with
> the ``statistic``, ``df`` (``k - 1``) and the chi-square ``p_value``.

### `kruskal_wallis_test(*groups)`  _function_

> Kruskal-Wallis H test that ``k`` groups share a distribution.
>
> Pass each group as a separate sequence argument. Returns a dict with the
> tie-corrected ``statistic`` H, the ``df`` (``k - 1``) and the chi-square upper-tail
> ``p_value``. A small p-value rejects the null of equal distributions (specifically,
> equal medians for similarly-shaped groups).

## lambert

### `lambert_w0(x)`  _function_

> Principal branch ``W0(x)`` solving ``w e^w = x`` with ``w >= -1``. Defined for ``x >= -1/e``.

### `lambert_wm1(x)`  _function_

> Secondary branch ``W_{-1}(x)`` with ``w <= -1``. Defined for ``-1/e <= x < 0``.

## laplace_inversion

### `laplace_inversion(F, t, N=12)`  _function_

> Invert the Laplace transform ``F(s)`` at time ``t`` by Gaver-Stehfest.
>
> ``F`` is a callable of a single real argument ``s > 0``. ``N`` (even, default 12) is
> the number of terms. Returns the estimated ``f(t)``. Best for smooth,
> non-oscillatory functions; oscillatory or discontinuous ``f`` invert poorly.

## lasso

### `elastic_net(X, y, alpha=1.0, l1_ratio=0.5, max_iter=1000, tol=1e-08)`  _function_

> Fit an elastic-net regression: a mix of L1 (LASSO) and L2 (ridge) penalties.
>
> Minimizes ``(1/2n) ||y - X beta||^2 + alpha (l1_ratio ||beta||_1 +
> 0.5 (1 - l1_ratio) ||beta||^2)`` by coordinate descent. ``l1_ratio = 1`` reduces
> to :func:`lasso_regression` (pure L1, sparse); ``l1_ratio = 0`` is a ridge-style
> L2 shrinkage. The L2 part groups correlated predictors while the L1 part still
> selects, which is more stable than pure LASSO when features are collinear.
> Returns ``[intercept, b_1, ..., b_p]`` on the original scale; the intercept is
> unpenalized.

### `lasso_regression(X, y, alpha=1.0, max_iter=1000, tol=1e-08)`  _function_

> Fit a LASSO regression ``y ~ X beta`` by coordinate descent.
>
> Returns ``[intercept, b_1, ..., b_p]`` on the original feature scale. ``alpha``
> is the L1 penalty strength: ``alpha = 0`` recovers ordinary least squares, and a
> large ``alpha`` drives all slopes to zero (the fit collapses to the mean of
> ``y``). Features are standardized internally so the penalty applies evenly; the
> intercept is never penalized. A slope set to exactly zero has been selected out.

## lazy_segtree

### `LazySegmentTree(values, mode='sum')`  _class_

> Segment tree with range-add updates and a range aggregate query.
>
> ``mode`` is ``"sum"``, ``"min"``, or ``"max"``. Build from an initial list; then
> ``update(lo, hi, delta)`` adds ``delta`` to every index in ``[lo, hi)`` and
> ``query(lo, hi)`` returns the aggregate over ``[lo, hi)`` -- both ``O(log n)``.

## lbfgs

### `lbfgs(func, x0, grad=None, m=10, tol=1e-08, max_iter=500)`  _function_

> Minimize ``func`` from ``x0`` by L-BFGS with a strong-Wolfe line search.
>
> ``func`` maps a length-``n`` list to a scalar. ``grad`` is an optional gradient function
> (same signature, returning a length-``n`` list); if omitted, central differences are used.
> ``m`` is the history size. Returns a dict with ``x`` (minimizer), ``fun``, ``n_iter``,
> ``converged`` (gradient norm below ``tol``) and ``grad_norm``.

## lca

### `LCA(adjacency, root)`  _class_

> Lowest-common-ancestor index over a rooted tree.
>
> Construct from an adjacency map ``{node: [neighbors]}`` (an undirected tree) and a root.
> ``query(u, v)`` returns their LCA, ``depth(u)`` the edge count from the root, and
> ``distance(u, v)`` the number of edges on the path between two nodes. Nodes may be any
> hashable label.

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

## legendre_harmonics

### `assoc_legendre(l, m, x)`  _function_

> Associated Legendre function ``P_l^m(x)`` for ``0 <= m <= l`` and ``|x| <= 1``.
>
> Uses the Condon-Shortley phase ``(-1)^m``. Seeds ``P_m^m = (-1)^m (2m-1)!! (1-x^2)^{m/2}``
> then climbs in ``l`` with ``(l-m) P_l^m = x (2l-1) P_{l-1}^m - (l+m-1) P_{l-2}^m``.
> ``assoc_legendre(l, 0, x) == legendre_p(l, x)``.

### `legendre_p(l, x)`  _function_

> Legendre polynomial ``P_l(x)`` via the three-term recurrence.
>
> ``P_0 = 1``, ``P_1 = x``, ``(l+1) P_{l+1} = (2l+1) x P_l - l P_{l-1}``. Valid for any real
> ``x`` (orthogonal on ``[-1, 1]``). ``P_l(1) = 1`` and ``P_l(-1) = (-1)^l``.

### `spherical_harmonic_real(l, m, theta, phi)`  _function_

> Real spherical harmonic ``Y_l^m(theta, phi)`` (orthonormal, ``-l <= m <= l``).
>
> ``theta`` is the polar (colatitude) angle in ``[0, pi]``, ``phi`` the azimuth. The real
> (tesseral) convention:
> ``m > 0``: ``sqrt(2) N_l^m P_l^m(cos theta) cos(m phi)``;
> ``m = 0``: ``N_l^0 P_l^0(cos theta)``;
> ``m < 0``: ``sqrt(2) N_l^|m| P_l^|m|(cos theta) sin(|m| phi)``,
> with ``N_l^m = sqrt((2l+1)/(4 pi) (l-m)!/(l+m)!)``. These are orthonormal over the sphere.

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

## levenberg

### `levenberg_marquardt(model, xs, ys, beta0, max_iter=100, tol=1e-10, lam0=0.001)`  _function_

> Fit ``beta`` minimizing ``sum_i (model(beta, xs[i]) - ys[i])^2``.
>
> Parameters
> ----------
> model : callable
>     ``model(beta, x)`` returning a scalar prediction.
> xs, ys : sequences
>     Inputs and targets of equal length.
> beta0 : sequence of float
>     Initial parameter guess.
> max_iter, tol : int, float
>     Iteration cap and convergence tolerance on the parameter step.
> lam0 : float
>     Initial damping.
>
> Returns
> -------
> dict
>     ``parameters``, ``residual`` (sum of squared errors), ``iterations``,
>     ``converged``.

## leveraged_etf

### `expected_leveraged_return(underlying_return, leverage, sigma, periods)`  _function_

> Approximate multi-period leveraged-ETF log return with drag.
>
> ``periods * (leverage * mu - drag)`` where ``mu`` is the per-period underlying
> log return and ``drag`` is :func:`volatility_drag`. Below naive
> ``periods * leverage * mu`` whenever there is drag.

### `flat_market_decay(leverage, returns)`  _function_

> Cumulative leveraged return over a *round-trip* (net-flat) return path.
>
> A market that ends where it started but moved in between: the leveraged ETF
> still loses value to the drag. Returns the final growth factor minus one; below
> zero for ``|leverage| > 1`` over a volatile flat path.

### `leveraged_etf_path(underlying_returns, leverage, expense_ratio=0.0, periods_per_year=252)`  _function_

> Daily-rebalanced leveraged ETF cumulative return path.
>
> Each period the ETF returns ``leverage * r - expense_ratio/periods_per_year``;
> the path compounds those. ``leverage`` may be negative (inverse ETFs). Returns
> the list of cumulative growth factors (starting after the first period).

### `volatility_drag(leverage, sigma)`  _function_

> Approximate per-period volatility drag of a leveraged ETF.
>
> ``0.5 * leverage * (leverage - 1) * sigma^2`` -- the expected log-return
> shortfall versus naive ``leverage`` times the underlying's log return, from the
> daily-rebalancing compounding. Zero at ``leverage`` 0 or 1; positive (a drag)
> for ``leverage > 1`` or ``leverage < 0``.

## levinson

### `levinson_durbin(autocorr)`  _function_

> Fit an AR model from an autocorrelation sequence ``[r_0, r_1, ..., r_p]``.
>
> Returns ``(ar_coeffs, error, reflection)``: the ``p`` AR coefficients ``a_1..a_p`` such
> that the model predicts ``x_n = sum_k a_k x_{n-k}``, the final prediction error variance,
> and the ``p`` reflection (PARCOR) coefficients. ``O(p^2)``.

### `solve_toeplitz(r, b)`  _function_

> Solve ``T x = b`` where ``T`` is the symmetric Toeplitz matrix with first row ``r``.
>
> ``r[0]`` is the diagonal, ``r[k]`` the ``k``-th off-diagonal (both directions, symmetric).
> ``len(r) == len(b) == n``. ``O(n^2)`` via the Levinson recursion. Raises if a leading
> principal minor is singular (``r[0] == 0`` or a zero prediction error).

## levy_basket

### `levy_basket_option(spots, weights, strike, t, r, sigmas, corr, q=None, is_call=True)`  _function_

> Levy moment-matched basket call/put on ``sum_i w_i S_i``.
>
> Parameters
> ----------
> spots, weights, sigmas : sequences of length ``n``.
> strike, t, r : option strike, maturity, risk-free rate.
> corr : ``n x n`` correlation matrix.
> q : optional per-asset dividend yields (defaults to zeros).
> is_call : call if True, else put.
>
> Returns
> -------
> float
>     Basket option value. Reduces to Black-Scholes for a single asset.

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

## lhs

### `l2_star_discrepancy(points)`  _function_

> L2 star discrepancy of a point set in ``[0, 1]^dim`` (Warnock's formula).
>
> A single number measuring departure from perfect uniformity -- smaller is more
> uniform. Uses Warnock's closed form, so it is exact and ``O(n^2 * dim)``. Useful for
> comparing a Latin hypercube (or Sobol/Halton set) against plain random sampling.

### `latin_hypercube(n, dim, seed=1234567, centered=False)`  _function_

> Latin hypercube design of ``n`` points in ``[0, 1)^dim``.
>
> Each of the ``dim`` axes is split into ``n`` equal bins with exactly one sample per
> bin, then the per-axis bin assignments are independently permuted. With
> ``centered=True`` each sample sits at its bin center (deterministic given the
> permutation); otherwise it is placed uniformly at random within its bin. Returns a
> list of ``n`` points (each a length-``dim`` list).

### `maximin_lhs(n, dim, seed=1234567, tries=20)`  _function_

> Latin hypercube maximizing the minimum inter-point distance over ``tries`` designs.
>
> Generates ``tries`` Latin hypercubes (distinct seeds) and returns the one whose
> closest pair is farthest apart -- a more space-filling design than a single random
> LHS, at ``tries`` times the cost. Returns the chosen point list.

## li_chao

### `LiChaoTree(xs, maximize=False)`  _class_

> Lower/upper envelope of lines over a fixed sorted list of query points.
>
> Construct with the ``xs`` at which queries will be made (any reals, deduplicated and
> sorted internally). ``add_line(m, b)`` inserts ``y = m x + b``; ``query(x)`` returns the
> minimum (or maximum, if ``maximize``) ``y`` over all inserted lines at that ``x``.
> ``x`` must be one of the construction points.

## linalg

### `basket_option_mc(spots, weights, strike, t, r, sigmas, correlation, option_type='call', n_paths=50000, seed=987654321)`  _function_

> Monte Carlo basket option on ``sum_i w_i S_i`` for any number of assets.
>
> Simulates correlated terminal asset prices under geometric Brownian motion
> (correlated normals via :func:`correlated_normals` / Cholesky) and averages the
> discounted basket payoff. A general n-asset reference for the two-asset
> analytic :func:`quantforge.basket_option`. Deterministic per seed.

### `cholesky(matrix)`  _function_

> Lower-triangular Cholesky factor ``L`` with ``L L^T = matrix``.
>
> Requires a symmetric positive-definite input; raises if a non-positive pivot
> is encountered (matrix not PD).

### `correlated_normals(independent, correlation)`  _function_

> Turn independent standard normals into correlated ones via Cholesky.
>
> Returns ``L z`` where ``L`` is the Cholesky factor of ``correlation`` -- the
> resulting vector has the target correlation structure. ``independent`` is a
> vector of IID standard normal draws.

### `is_positive_definite(matrix)`  _function_

> True if the symmetric matrix is positive definite (Cholesky succeeds).

### `nearest_correlation(matrix, max_iter=100, tol=1e-10)`  _function_

> Nearest positive-semidefinite correlation matrix (Higham alternating projection).
>
> Repairs an indefinite estimated correlation matrix to the closest valid one:
> alternately projects onto the PSD cone (clip negative eigenvalues to zero) and
> onto the unit-diagonal set, iterating to convergence. Returns a symmetric PSD
> matrix with unit diagonal; a matrix that is already a valid correlation is
> returned essentially unchanged.

## linear_recurrence

### `fibonacci(n, mod=None)`  _function_

> The ``n``-th Fibonacci number (``F_0 = 0, F_1 = 1``) via matrix power.
>
> Optional modular reduction. ``O(log n)`` matrix multiplies.

### `linear_recurrence_nth(coeffs, initial, n, mod=None)`  _function_

> The ``n``-th term (0-indexed) of ``x_i = sum_j coeffs[j] * x_{i-1-j}``.
>
> ``coeffs`` are the ``k`` recurrence coefficients (``coeffs[0]`` multiplies the most
> recent term); ``initial`` are the first ``k`` terms ``x_0 .. x_{k-1}``. For ``n < k`` the
> initial term is returned directly. With ``mod`` set, the result is reduced modulo it.
> ``O(k^3 log n)``.

### `matrix_power(matrix, power, mod=None)`  _function_

> Raise a square ``matrix`` to a non-negative integer ``power`` by squaring.
>
> ``power == 0`` returns the identity. With ``mod`` set, every entry is reduced modulo it
> (integer matrices only). ``O(k^3 log power)``.

## linprog

### `linprog(c, constraints, maximize=True, tol=1e-09, max_iter=10000)`  _function_

> Solve a linear program by two-phase simplex.
>
> ``c`` is the objective coefficient vector (length ``n``). ``constraints`` is a list of
> ``(coeffs, op, rhs)`` with ``op`` one of ``"<="``, ``">="``, ``"="``. Variables are
> assumed ``>= 0``. Set ``maximize=False`` to minimize. Returns
> ``{"x": [...], "objective": value, "status": "optimal"}``. Raises ``ValueError`` if
> the program is infeasible or unbounded.

## liquidity

### `amihud_illiquidity(returns, dollar_volumes)`  _function_

> Amihud illiquidity: average of ``|return| / dollar_volume`` over the period.
>
> A larger value means a given dollar of trading moves the price more (a less
> liquid asset). ``returns`` and ``dollar_volumes`` are aligned daily series.

### `corwin_schultz_spread(highs, lows)`  _function_

> Corwin-Schultz high-low bid-ask spread estimator.
>
> Uses consecutive daily high/low ranges: the two-day range reflects both
> volatility and the spread, while single-day ranges scale with volatility alone,
> so their combination isolates the spread. Returns the average estimated
> proportional spread over the sample, floored at zero each day (negative daily
> estimates, a known small-sample artefact, are set to zero). ``highs`` and
> ``lows`` are aligned and at least two long.

### `roll_spread(prices)`  _function_

> Roll's implied effective spread from a series of transaction prices.
>
> Bid-ask bounce makes successive price changes negatively autocovaried; Roll's
> estimator is ``spread = 2 sqrt(-cov(dP_t, dP_{t-1}))``. When the sample
> autocovariance is non-negative (no detectable bounce) the estimate is zero.
> ``prices`` are levels; at least three are needed.

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

## logistic

### `fit_logistic(X, y, add_intercept=True, max_iter=100, tol=1e-08, ridge=1e-08)`  _function_

> Fit a logistic regression by IRLS / Newton-Raphson.
>
> Parameters
> ----------
> X : list[list[float]]
>     Design matrix, ``n`` rows of features.
> y : list[float]
>     Binary outcomes (0 or 1).
> add_intercept : bool
>     Prepend an intercept column.
> max_iter, tol : int, float
>     Newton iteration cap and convergence tolerance on the coefficient step.
> ridge : float
>     Small L2 penalty on the Hessian for numerical stability under separation.
>
> Returns
> -------
> dict
>     ``coefficients`` (intercept first if added), ``iterations``,
>     ``log_likelihood``, ``converged``.

### `predict_proba(model, X, add_intercept=True)`  _function_

> Predicted P(y=1) for each row of ``X`` under a fitted logistic model.

## logsumexp

### `log_softmax(x)`  _function_

> Stable log-softmax: ``x_i - logsumexp(x)`` (avoids the overflow of ``log(softmax)``).

### `logsumexp(x, weights=None)`  _function_

> Stable ``log(sum_i w_i exp(x_i))`` (weights default to 1).
>
> Shifts by the maximum so it never overflows; returns ``-inf`` if every weighted
> term is zero. With ``weights`` it is the log of a weighted sum of exponentials
> (weights must be non-negative).

### `softmax(x)`  _function_

> Stable softmax: ``exp(x_i) / sum_j exp(x_j)``, a probability vector summing to 1.

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

## lowess

### `lowess(x, y, frac=0.3, iterations=3)`  _function_

> LOWESS smooth of ``(x, y)``; returns the fitted value at each ``x``.
>
> ``frac`` is the span -- the fraction of points in each local neighbourhood (larger
> = smoother). ``iterations`` robustifying passes (Cleveland) down-weight outliers by
> a bisquare of their residuals; ``iterations=1`` disables robustifying. Points need
> not be sorted. Returns a list aligned to the input order.

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

## lu

### `determinant(A)`  _function_

> Determinant of a square matrix as the signed product of the LU pivots.
>
> ``det = sign * prod(U[i][i])``. Returns 0 for a singular matrix.

### `lu_decomposition(A)`  _function_

> LU decomposition with partial pivoting of a square matrix.
>
> Returns ``(L, U, piv, sign)`` where ``L`` is unit lower-triangular, ``U`` upper-
> triangular, ``piv`` the row-permutation (as a list mapping output row -> source
> row) so that ``A[piv] = L U``, and ``sign`` the permutation parity (+1/-1) used
> by :func:`determinant`. Raises on a singular matrix.

### `lu_solve(A, b)`  _function_

> Solve ``A x = b`` via LU with partial pivoting.
>
> Factorizes ``A``, permutes ``b``, then forward- and back-substitutes. Returns the
> solution vector. Raises on a singular matrix.

## lzw

### `delta_decode(deltas)`  _function_

> Invert :func:`delta_encode`: cumulative sum of the deltas.

### `delta_encode(data)`  _function_

> Delta-encode a numeric sequence: first value, then successive differences.
>
> ``[a, b, c, ...] -> [a, b-a, c-b, ...]``. Slowly-varying data becomes small numbers
> that compress better. Empty input yields an empty list.

### `lzw_compress(data)`  _function_

> LZW-compress a string into a list of integer codes.
>
> Starts with a dictionary of the distinct single characters (assigned codes in sorted
> order) and grows it with each new substring seen, emitting the code of the longest
> known prefix. Returns ``(codes, alphabet)`` -- the sorted initial alphabet is needed
> to seed the decoder. Empty input yields ``([], [])``.

### `lzw_decompress(codes, alphabet)`  _function_

> Decompress LZW ``codes`` given the initial ``alphabet`` back to the string.
>
> Rebuilds the same dictionary the compressor grew, handling the special case where a
> code refers to an entry being defined this step. Inverts :func:`lzw_compress`
> exactly.

## mack

### `mack_standard_error(triangle)`  _function_

> Chain-ladder reserves with Mack's prediction standard errors.
>
> ``triangle`` is a cumulative run-off triangle: ``triangle[i]`` holds the observed
> cumulative claims for accident year ``i`` at development ages ``0 ..``, with the
> most recent accident year the shortest row. Returns a dict with
>
> - ``factors``, ``sigma2``: estimated age-to-age factors and variance parameters,
> - ``ultimate``, ``reserve``: per accident year,
> - ``std_error``: Mack standard error of each accident-year reserve,
> - ``cv``: coefficient of variation ``std_error / reserve`` (0 for a zero reserve),
> - ``total_reserve``, ``total_std_error``, ``total_cv``: for the reserve sum,
>   including the between-year correlation term.
>
> Raises ``ValueError`` for fewer than two accident years or non-positive entries.

## markov

### `absorption_probabilities(P, transient_states, absorbing_states)`  _function_

> Probability of ending in each absorbing state from each transient state.
>
> ``B = N R`` where ``N`` is the :func:`fundamental_matrix` and ``R`` is the
> transient-to-absorbing transition block. Row ``i`` (a transient state) is a
> distribution over ``absorbing_states`` summing to one.

### `cumulative_default_term_structure(P, default_state, horizons, start_state=0)`  _function_

> Cumulative default probability by horizon from a rating-migration matrix.
>
> ``P`` is a one-period rating transition matrix with ``default_state`` an
> absorbing default row. For each ``n`` in ``horizons`` the cumulative default
> probability from ``start_state`` is the default-column entry of ``P^n``, i.e.
> ``(P^n)[start_state][default_state]``. Non-decreasing in the horizon (default
> is absorbing) and rising toward one if default is reachable.

### `expected_hitting_time(P, target)`  _function_

> Expected number of steps to first reach ``target`` from each state.
>
> Solves the linear system ``h_i = 0`` for ``i = target`` and
> ``h_i = 1 + sum_j P_ij h_j`` otherwise, by Gaussian elimination. Returns the
> vector of expected hitting times (``0`` at the target, ``inf`` conceptually if
> unreachable -- the solver raises on a singular system in that case).

### `expected_steps_to_absorption(P, transient_states)`  _function_

> Expected steps to absorption from each transient state.
>
> Row sums of the :func:`fundamental_matrix` ``N`` -- the total expected visits
> across all transient states before hitting an absorbing state. Positive for
> every transient state.

### `fundamental_matrix(P, transient_states)`  _function_

> Fundamental matrix ``N = (I - Q)^{-1}`` of an absorbing chain.
>
> ``transient_states`` lists the indices of the non-absorbing states; ``Q`` is
> their sub-transition block. ``N_ij`` is the expected number of visits to
> transient state ``j`` starting from ``i`` before absorption. Requires the chain
> to be absorbing (every transient state eventually reaches an absorbing one).

### `generator_default_probability(Q, default_state, horizons, start_state=0)`  _function_

> Cumulative default probability at each horizon from a rating generator.
>
> Exponentiates ``Q`` to each horizon and reads the transition probability from
> ``start_state`` into the absorbing ``default_state``. Returns one probability per
> horizon; non-decreasing when the default state is absorbing.

### `generator_to_transition(Q, t=1.0)`  _function_

> Transition matrix ``P(t) = exp(Q t)`` of a continuous-time Markov chain.
>
> ``Q`` is a rate generator (rows summing to zero, non-negative off-diagonals).
> Returns the ``t``-horizon transition matrix, whose rows sum to one with
> non-negative entries. ``P(0)`` is the identity and ``P(s) P(t) = P(s + t)``
> (the semigroup property). The basis for continuous-time rating migration.

### `marginal_default_probabilities(P, default_state, horizons, start_state=0)`  _function_

> Marginal (per-period) default probabilities between successive horizons.
>
> Differences of the :func:`cumulative_default_term_structure`; each is the
> probability of defaulting in ``(horizons[k-1], horizons[k]]`` having survived
> to ``horizons[k-1]``. Non-negative because the cumulative curve is
> non-decreasing.

### `n_step_transition(P, n)`  _function_

> ``n``-step transition matrix ``P^n`` by repeated squaring.
>
> Row ``i`` gives the distribution over states after ``n`` steps starting from
> ``i``. ``P^0`` is the identity; every power is again row-stochastic.

### `stationary_distribution(P, tol=1e-14, max_iter=100000)`  _function_

> Stationary distribution ``pi`` with ``pi P = pi`` (power iteration).
>
> Iterates a uniform start under ``P`` until convergence. For an irreducible
> aperiodic chain this is the unique long-run state distribution; the returned
> vector is non-negative and sums to one.

## matched_filter

### `detect_template(x, template, threshold=0.8, distance=None)`  _function_

> Offsets where ``template`` occurs in ``x`` (normalized correlation above ``threshold``).
>
> Runs :func:`normalized_matched_filter` and returns the peak offsets whose correlation
> exceeds ``threshold`` (a coefficient in ``[-1, 1]``), separated by at least
> ``distance`` samples (default ``len(template)``, so overlapping detections of the same
> hit collapse to one). Each returned offset is the start index of a match.
>
> Because the correlation is amplitude-normalized, short or featureless templates match
> many noise windows by shape alone; use a longer, distinctive template (and a higher
> ``threshold``) when false positives matter.

### `find_peaks(x, height=None, distance=1)`  _function_

> Indices of local maxima in ``x``, optionally filtered by height and spacing.
>
> A point is a peak if it is strictly greater than both neighbors. ``height`` drops
> peaks below that threshold; ``distance`` enforces a minimum gap between reported peaks
> by keeping the taller peak when two fall within ``distance`` samples. Returns a sorted
> list of indices.

### `matched_filter(x, template)`  _function_

> Matched-filter response: correlate ``x`` against ``template`` at every offset.
>
> Returns a list of length ``len(x) - len(template) + 1`` whose entry ``i`` is the dot
> product of ``template`` with the window ``x[i:i+len(template)]``. The maximum locates
> the offset where the template best aligns -- the maximum-SNR detector for a known
> shape in white noise.

### `normalized_matched_filter(x, template)`  _function_

> Normalized matched filter: the correlation coefficient at each offset, in ``[-1, 1]``.
>
> Removes each window's mean and the template mean, then divides by the norms, so the
> response is a dimensionless correlation independent of signal amplitude and offset --
> ``1`` where the window is a scaled, shifted copy of the template. Returns a list of
> length ``len(x) - len(template) + 1``.

## matrix_exp

### `matrix_exp(A)`  _function_

> Matrix exponential ``exp(A)`` by scaling-and-squaring with Pade(6,6).
>
> ``A`` is a square matrix. Returns ``exp(A)``; ``exp(0) = I``, ``exp`` of a
> diagonal matrix is the diagonal of exponentials, and it satisfies the defining
> series. Accurate across a wide norm range thanks to the scaling step.

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

## mcmc

### `hamiltonian_monte_carlo(log_prob, x0, n_samples, step=0.1, n_leapfrog=20, seed=12345, burn_in=0)`  _function_

> Hamiltonian Monte Carlo sampler using autodiff gradients of ``log_prob``.
>
> ``log_prob`` maps a length-``d`` list of :class:`quantforge.reverse_ad.Var` to a single
> ``Var`` (the log target density up to a constant). Each iteration samples a Gaussian momentum,
> runs ``n_leapfrog`` leapfrog steps of step size ``step`` on the Hamiltonian
> ``H = -log_prob(x) + |p|^2 / 2``, and Metropolis-accepts on the total energy. Returns a dict
> with ``samples`` and ``accept_rate``.

### `metropolis_hastings(log_prob, x0, n_samples, step=0.5, seed=12345, burn_in=0, thin=1)`  _function_

> Random-walk Metropolis sampler for a target ``log_prob``.
>
> ``log_prob`` maps a length-``d`` list of floats to the log target density (up to an additive
> constant). Proposals are isotropic Gaussian with standard deviation ``step``. Returns a dict
> with ``samples`` (list of accepted states after burn-in/thinning) and ``accept_rate``.

### `sample_cov(samples)`  _function_

> Component-wise covariance matrix of a list of sample vectors (population form).

### `sample_mean(samples)`  _function_

> Component-wise mean of a list of sample vectors.

## mcmc_diagnostics

### `effective_sample_size(x, max_lag=None)`  _function_

> Effective sample size ``n / tau`` -- independent-draw equivalent of a correlated chain.
>
> Uses :func:`integrated_autocorrelation_time`. Equals ``n`` for white noise and shrinks as
> the chain's autocorrelation grows. Never exceeds ``n``.

### `gelman_rubin(chains)`  _function_

> Gelman-Rubin potential-scale-reduction ``R-hat`` for a list of scalar chains.
>
> ``chains`` is a list of ``m >= 2`` equal-length sequences (one per independent run). Returns
> ``R-hat = sqrt(V_hat / W)`` where ``W`` is the mean within-chain variance and ``V_hat`` the
> variance estimate mixing in the between-chain spread. ``-> 1`` at convergence.

### `integrated_autocorrelation_time(x, max_lag=None)`  _function_

> Integrated autocorrelation time ``tau = 1 + 2 sum_k rho_k`` (Geyer initial positive seq).
>
> Sums the autocorrelations, truncating at the first lag where consecutive-pair sums turn
> negative (Geyer's rule), which keeps the estimate stable. ``tau >= 1``; larger means more
> correlated samples. Returns ``1.0`` for a zero-variance series.

## median_filter

### `hampel_filter(x, window=7, n_sigmas=3.0)`  _function_

> Hampel outlier filter: replace points far from the local median with that median.
>
> In each centered window of length ``2*window+1`` it computes the median and the
> median absolute deviation (MAD), scales the MAD to a robust standard deviation
> (``1.4826 * MAD``), and replaces the center point only if it lies more than
> ``n_sigmas`` robust deviations away. Clean data passes through untouched. Returns
> ``(filtered, outlier_indices)``.

### `median_filter(x, window)`  _function_

> Sliding-window median of ``x`` with an odd ``window`` length.
>
> Replaces each point by the median of the ``window`` samples centered on it (the
> window is clipped at the ends). Removes impulsive spikes while keeping step edges
> sharp, unlike a moving average. Returns a list the same length as ``x``.

### `rank_filter(x, window, percentile)`  _function_

> Sliding-window order-statistic filter at a given ``percentile`` (0-100).
>
> ``percentile=50`` is the median; ``0`` a min filter, ``100`` a max filter. Picks the
> order statistic nearest that percentile in each centered window. Returns a list the
> same length as ``x``.

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

## microstructure

### `effective_spread(trade_prices, mids, signs)`  _function_

> Average effective (proportional) spread ``2 * sign * (price - mid) / mid``.
>
> The cost actually paid relative to the midpoint at the time of the trade:
> ``sign`` is ``+1`` for buys and ``-1`` for sells. Returned as the mean over the
> trades, in the same units as ``price / mid`` (a fraction). Wider than the quoted
> spread when trades walk the book, tighter when they occur inside it.

### `kyle_lambda_regression(price_changes, signed_volumes)`  _function_

> Empirical Kyle's lambda: slope of price change on signed order flow.
>
> Fits ``dP = alpha + lambda * signed_volume`` by ordinary least squares and
> returns ``lambda`` (price move per unit of net signed volume). ``signed_volumes``
> should be signed by trade direction (buys positive, sells negative). A larger
> lambda means a less liquid, higher-impact market. Requires at least three
> aligned observations with non-constant order flow.

### `order_flow_imbalance(buy_volumes, sell_volumes)`  _function_

> Net signed volume over total volume: ``sum(buy - sell) / sum(buy + sell)``.
>
> In ``[-1, 1]``: positive when buys dominate, negative when sells do. Aligned
> non-negative volume series.

### `price_impact(trade_prices, mids, future_mids, signs)`  _function_

> Average (proportional) price impact ``2 * sign * (mid_future - mid) / mid``.
>
> The permanent midpoint move in the trade's direction over the impact horizon --
> the informational half of the spread. By construction
> ``effective = realized + price_impact`` term by term (both defined with the same
> ``2 * sign / mid`` scaling), so this equals the gap between the effective and
> realized spreads.

### `quoted_spread(bids, asks)`  _function_

> Average proportional quoted spread ``(ask - bid) / midpoint``.
>
> The posted cost of a round trip, independent of where trades actually print.
> Aligned bid/ask series with positive midpoints.

### `realized_spread(trade_prices, mids, future_mids, signs)`  _function_

> Average realized (proportional) spread ``2 * sign * (price - mid_future) / mid``.
>
> The portion of the effective spread the liquidity provider *keeps* -- the trade
> price against the midpoint a short horizon later (``future_mids``), so it nets
> out the permanent price move. ``mids`` is the quote midpoint at the trade,
> ``future_mids`` the midpoint after the impact horizon. Mean over the trades.

### `vpin(buy_volumes, sell_volumes)`  _function_

> Volume-synchronized probability of informed trading (VPIN).
>
> Given per-bucket buy and sell volumes (equal-volume buckets), VPIN is the mean
> of ``|buy - sell| / (buy + sell)`` across buckets -- the average absolute order
> imbalance. In ``[0, 1]``: near zero when buys and sells balance, near one when
> trading is one-sided (a proxy for informed order flow / toxicity). Aligned
> non-negative series with at least one bucket.

## min_area_rect

### `min_area_rectangle(points)`  _function_

> Return the minimum-area enclosing rectangle of ``points``.
>
> Result is a dict with ``area``, ``width``, ``height`` (the two side lengths, width the
> longer), ``angle`` (edge direction in radians), and ``corners`` (four ``(x, y)`` points
> in order). Needs at least one point; a degenerate (collinear) set gives a zero-area
> rectangle.

## min_cost_flow

### `MinCostMaxFlow()`  _class_

> Min-cost max-flow on an integer-capacity network built edge by edge.
>
> ``add_edge(u, v, capacity, cost)`` adds a directed edge (and its residual). ``solve(s,
> t)`` returns ``(max_flow, min_cost)``: the maximum flow value and the least total cost
> achieving it. Nodes are any hashable labels.

### `min_cost_max_flow(edges, source, sink)`  _function_

> Convenience: build the network from ``edges`` = ``[(u, v, capacity, cost), ...]``.
>
> Returns ``(max_flow, min_cost)``.

## minimize1d

### `brent_min(f, lo, hi, tol=1e-10, max_iter=200)`  _function_

> Minimize a scalar ``f`` on ``[lo, hi]`` by Brent's method.
>
> Parabolic interpolation with a golden-section fallback. Returns
> ``(x_min, f_min)``.

### `golden_section_min(f, lo, hi, tol=1e-10, max_iter=200)`  _function_

> Minimize a unimodal ``f`` on ``[lo, hi]`` by golden-section search.

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

## modular

### `chinese_remainder(remainders, moduli)`  _function_

> Chinese Remainder Theorem: the unique ``x`` in ``[0, M)`` matching all congruences.
>
> Given ``x == remainders[i] (mod moduli[i])`` with pairwise-coprime ``moduli``,
> returns ``(x, M)`` where ``M`` is the product of the moduli. Raises if the moduli are
> not pairwise coprime.

### `discrete_log(base, target, mod)`  _function_

> Smallest non-negative ``x`` with ``base^x == target (mod mod)``, or ``None``.
>
> Baby-step giant-step: ``O(sqrt(mod))`` time and space. Searches exponents in
> ``[0, mod)``. ``mod`` must be positive; ``target`` is reduced mod ``mod``.

### `extended_gcd(a, b)`  _function_

> Extended Euclidean algorithm: return ``(g, x, y)`` with ``a*x + b*y == g = gcd(a, b)``.
>
> The Bezout coefficients ``x, y`` are the basis of the modular inverse and the CRT.

### `mod_inverse(a, m)`  _function_

> Modular inverse ``a^{-1} mod m``: the ``x`` with ``a*x == 1 (mod m)``.
>
> Exists iff ``gcd(a, m) == 1``; raises ``ValueError`` otherwise. Result is in
> ``[0, m)``.

### `mod_pow(base, exp, mod)`  _function_

> Modular exponentiation ``base^exp mod mod`` (supports negative ``exp`` via inverse).
>
> A thin, explicit wrapper over fast binary exponentiation; a negative exponent inverts
> the base first (requires ``gcd(base, mod) == 1``).

## modular2

### `jacobi_symbol(a, n)`  _function_

> Jacobi symbol ``(a/n)`` for an odd positive ``n``: 0, 1, or -1.
>
> Generalizes the Legendre symbol to composite (odd) ``n`` by multiplicativity in the
> denominator. Equals the Legendre symbol when ``n`` is prime.

### `legendre_symbol(a, p)`  _function_

> Legendre symbol ``(a/p)`` for an odd prime ``p``: 0, 1, or -1.
>
> ``1`` if ``a`` is a non-zero quadratic residue mod ``p``, ``-1`` if a non-residue,
> ``0`` if ``a`` is divisible by ``p``. Uses Euler's criterion ``a^((p-1)/2) mod p``.

### `multiplicative_order(a, n)`  _function_

> Multiplicative order of ``a`` modulo ``n``: the least ``k > 0`` with ``a^k == 1``.
>
> Requires ``gcd(a, n) == 1``. Computes it from the factorization of Euler's totient of
> a prime (here ``n`` is required prime for an exact ``phi = n - 1``); raises otherwise.

### `primitive_root(p)`  _function_

> A primitive root modulo the prime ``p`` (a generator of the multiplicative group).
>
> Returns the smallest ``g`` whose multiplicative order is ``p - 1``. Every power of a
> primitive root covers all non-zero residues exactly once.

### `tonelli_shanks(a, p)`  _function_

> A square root of ``a`` modulo an odd prime ``p`` (Tonelli-Shanks).
>
> Returns ``r`` with ``r*r == a (mod p)`` (the other root is ``p - r``). Raises if ``a``
> is a non-residue. ``a`` is reduced mod ``p`` first.

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

## money_market

### `bank_discount_yield(face, price, days, year_days=360)`  _function_

> Bank discount yield from the price: ``(face - price)/face * year_days/days``.
>
> Inverse of :func:`price_from_discount`. Understates the true return because it
> divides the discount by face rather than price.

### `bond_equivalent_yield(face, price, days, year_days=365)`  _function_

> Bond-equivalent yield: return on price on an actual/365 basis.
>
> ``(face - price)/price * 365/days`` -- puts a discount instrument on the same
> footing as a coupon bond (actual/365), so it exceeds both the discount and the
> actual/360 money-market yield.

### `discount_to_bond_equivalent(discount_rate, days)`  _function_

> Convert a bank discount rate directly to a bond-equivalent yield.
>
> Prices at par 100 off the discount, then takes the actual/365 return on price.
> Always above the input discount rate.

### `holding_period_return(buy_price, sell_price, income=0.0)`  _function_

> Holding-period return ``(sell - buy + income) / buy``.
>
> The total return over the holding period including any interim income; not
> annualized.

### `money_market_yield(face, price, days, year_days=360)`  _function_

> Money-market (CD-equivalent) yield: return on *price*, actual/360.
>
> ``(face - price)/price * year_days/days`` -- the actual return per invested
> dollar, higher than the :func:`bank_discount_yield` (which divides by face).

### `price_from_discount(face, discount_rate, days, year_days=360)`  _function_

> Price of a discount instrument from its bank discount rate.
>
> ``price = face * (1 - discount_rate * days / year_days)`` -- the bank discount
> convention prices off the *face*, not the price, on an actual/360 basis.

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

## multiple_testing

### `benjamini_hochberg(pvals)`  _function_

> Benjamini-Hochberg FDR-adjusted p-values (step-up).
>
> Sorts ascending, scales the ``k``-th smallest (1-based) by ``m / k``, takes a
> running minimum from the largest down so the sequence is monotone, and unshuffles
> to the input order. Controls the false-discovery rate under independence or
> positive dependence; less conservative than family-wise methods.

### `benjamini_yekutieli(pvals)`  _function_

> Benjamini-Yekutieli FDR-adjusted p-values (arbitrary dependence).
>
> As Benjamini-Hochberg but with the extra factor ``c(m) = sum_{i=1}^m 1/i``, which
> makes the procedure valid under any dependence structure at the cost of power.

### `bonferroni(pvals)`  _function_

> Bonferroni-adjusted p-values: ``min(1, m * p)`` for ``m`` tests.

### `holm(pvals)`  _function_

> Holm step-down family-wise adjusted p-values.
>
> Sorts ascending, scales the ``k``-th smallest (0-based) by ``m - k``, then takes
> a running maximum so the sequence is monotone, and unshuffles to the input
> order. Controls the family-wise error rate and dominates :func:`bonferroni`.

## multistep_ode

### `adams_bashforth_moulton(f, y0, t0, t1, n_steps)`  _function_

> Integrate ``y' = f(t, y)`` from ``t0`` to ``t1`` in ``n_steps`` using ABM4 (PECE).
>
> ``y0`` is scalar or a list (vector system). Returns ``(ts, ys)``: the ``n_steps + 1``
> time points and the solution at each. The first three steps use RK4 to build history.

## multivariate_normal_cdf

### `bivariate_normal_cdf(a, b, rho)`  _function_

> Standard bivariate normal CDF ``P(X1 <= a, X2 <= b; corr=rho)``.
>
> ``rho`` in ``[-1, 1]``. Accurate to ~1e-7 across the usable correlation range.

### `trivariate_normal_cdf(a, b, c, r12, r13, r23, n=24)`  _function_

> Standard trivariate normal CDF ``P(X1<=a, X2<=b, X3<=c)``.
>
> ``r12, r13, r23`` are the pairwise correlations of a valid 3x3 correlation matrix.
> Uses the reduction ``P3 = P(X3<=c) * ...`` via a 1-D integral over the third
> variable of a conditional bivariate CDF (Genz). ``n`` sets the quadrature panels.
> Falls back to the product/independent forms when correlations vanish.

## mvn

### `log_determinant(cov)`  _function_

> Log-determinant of a symmetric positive-definite matrix via Cholesky.
>
> ``log det Sigma = 2 sum_i log L_ii``. Numerically stable where a direct product
> of eigenvalues (or the determinant) would under/overflow. Raises if ``cov`` is
> not positive definite.

### `mvn_logpdf(x, mean, cov)`  _function_

> Log-density of the multivariate normal ``N(mean, cov)`` at ``x``.
>
> ``-0.5 [ k ln(2 pi) + ln|Sigma| + (x-mu)' Sigma^{-1} (x-mu) ]``. The quadratic
> form is evaluated as ``||L^{-1}(x-mu)||^2`` from the Cholesky factor, avoiding an
> explicit inverse. Reduces to the univariate normal log-density for ``k = 1``.

### `mvn_pdf(x, mean, cov)`  _function_

> Density of the multivariate normal ``N(mean, cov)`` at ``x`` (``exp`` of the log).

## naive_bayes

### `fit_gaussian_nb(X, y, var_smoothing=1e-09)`  _function_

> Fit a Gaussian naive Bayes model.
>
> Parameters
> ----------
> X : list[list[float]]
>     ``n`` rows of ``d`` features.
> y : list
>     Class labels (any hashable).
> var_smoothing : float
>     Added to every variance for numerical stability (avoids zero variance on
>     constant features).
>
> Returns
> -------
> dict
>     ``classes`` (sorted), ``priors``, ``means`` and ``variances`` (per class,
>     per feature), keyed by class label.

### `predict_gaussian_nb(model, X_query)`  _function_

> Predict the most probable class for each query row.

### `predict_proba_gaussian_nb(model, X_query)`  _function_

> Posterior class probabilities per query (softmax of the log posteriors).

## nelson_siegel

### `fit_nelson_siegel(maturities, zero_rates, tau_grid=None)`  _function_

> Least-squares fit of Nelson-Siegel parameters to observed zero rates.
>
> For a fixed decay ``tau`` the three betas enter linearly (the level/slope/
> curvature loadings), so they are solved by ordinary least squares; ``tau`` is
> chosen by a grid search minimizing the residual sum of squares. Returns
> ``(beta0, beta1, beta2, tau)``. Recovers the true parameters exactly on
> noiseless data whose ``tau`` is in the grid.

### `nelson_siegel_discount(t, beta0, beta1, beta2, tau)`  _function_

> Discount factor ``exp(-z(t) t)`` from the Nelson-Siegel zero rate.

### `nelson_siegel_forward(t, beta0, beta1, beta2, tau)`  _function_

> Instantaneous forward rate under Nelson-Siegel.
>
> ``f(t) = beta0 + beta1 e^{-t/tau} + beta2 (t/tau) e^{-t/tau}``. Equals the
> short rate ``beta0 + beta1`` at ``t = 0`` and the long level ``beta0`` as
> ``t -> inf``.

### `nelson_siegel_zero(t, beta0, beta1, beta2, tau)`  _function_

> Nelson-Siegel zero rate at maturity ``t``.
>
> ``z(t) = beta0 + (beta1 + beta2) (1 - e^{-t/tau}) / (t/tau) - beta2 e^{-t/tau}``.
> ``beta0`` is the long-run level, ``beta0 + beta1`` the short rate (``t -> 0``),
> and ``beta2`` scales the medium-term curvature hump with decay ``tau``.

### `svensson_zero(t, beta0, beta1, beta2, beta3, tau1, tau2)`  _function_

> Svensson zero rate: Nelson-Siegel plus a second curvature term.
>
> Adds ``beta3 ((1 - e^{-t/tau2})/(t/tau2) - e^{-t/tau2})`` with its own decay
> ``tau2`` for a second hump. Reduces to :func:`nelson_siegel_zero` when
> ``beta3 = 0``.

## newton_min

### `newton_min(func, x0, tol=1e-09, max_iter=100, h=1e-05)`  _function_

> Minimize ``func`` from ``x0`` by damped Newton with exact autodiff derivatives.
>
> ``func`` takes a list of :class:`quantforge.reverse_ad.Var` and returns a single ``Var``.
> The gradient is exact (reverse mode); the Hessian is the reverse gradient differenced once
> (:func:`quantforge.reverse_jacobian.reverse_hessian`). Returns a dict with ``x``, ``fun``,
> ``n_iter``, ``converged`` (gradient norm below ``tol``) and ``grad_norm``.

## newton_system

### `broyden(f, x0, tol=1e-10, max_iter=200, rel_step=1e-06)`  _function_

> Solve ``f(x) = 0`` by Broyden's (good) quasi-Newton method.
>
> Seeds the inverse-Jacobian estimate from one finite-difference Jacobian, then updates
> it rank-1 from each step's secant equation -- avoiding a fresh Jacobian per iteration.
> Returns ``(solution, iterations)``. Raises if the seed Jacobian is singular or it
> fails to converge.

### `newton_system(f, x0, tol=1e-10, max_iter=100, rel_step=1e-06)`  _function_

> Solve ``f(x) = 0`` for a vector function by Newton's method.
>
> ``f`` maps a length-``n`` list to a length-``n`` list; ``x0`` is the initial guess.
> Each step solves ``J dx = -f(x)`` with the central-difference Jacobian ``J`` and
> updates ``x += dx``. Returns ``(solution, iterations)``. Raises if the Jacobian is
> singular or convergence is not reached within ``max_iter``.

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

## nnls

### `nnls(A, b, max_iter=None, tol=1e-10)`  _function_

> Non-negative least squares: ``min ||A x - b||^2`` with ``x >= 0``.
>
> ``A`` is an ``m x n`` matrix (list of rows), ``b`` a length-``m`` vector. Returns a
> dict with the non-negative solution ``x`` (length ``n``), the ``residual_norm``
> ``||A x - b||``, and ``n_iter``. Uses the Lawson-Hanson active-set method.

## ntt

### `convolve_mod(a, b, mod=998244353)`  _function_

> Exact convolution (polynomial product) of integer sequences ``a`` and ``b`` mod ``mod``.
>
> Returns the ``len(a) + len(b) - 1`` coefficients of the product, each reduced mod
> ``mod`` -- no floating-point error. Empty inputs give an empty result.

### `intt(a, mod=998244353)`  _function_

> Inverse number-theoretic transform (with the ``1/n`` field scaling).

### `ntt(a, mod=998244353)`  _function_

> Forward number-theoretic transform of ``a`` (length must be a power of two).

## number_theory

### `divisors(n)`  _function_

> All positive divisors of ``n`` in sorted order (from its factorization).

### `euler_totient(n)`  _function_

> Euler's totient ``phi(n)``: the count of integers in ``[1, n]`` coprime to ``n``.
>
> Computed from the factorization as ``n * prod (1 - 1/p)`` over distinct primes ``p``.
> ``phi(1) = 1``.

### `factorize(n)`  _function_

> Prime factorization of ``n`` as a sorted list of ``(prime, exponent)`` pairs.
>
> Combines trial division by small primes with Pollard's rho for large factors and
> Miller-Rabin to certify primality, so it factors numbers far beyond what naive trial
> division reaches. ``n`` must be ``>= 1`` (``1`` has an empty factorization).

### `gcd(a, b)`  _function_

> Greatest common divisor of ``a`` and ``b`` (Euclid's algorithm, non-negative).

### `is_prime(n)`  _function_

> Deterministic Miller-Rabin primality test, exact for all ``n < 3.3 * 10^24``.
>
> Uses a fixed set of witness bases that is proven to give no false positives below
> that bound (well past 64-bit), so the answer is exact, not probabilistic, for any
> integer arising in normal use. ``O(k log^3 n)``.

### `lcm(a, b)`  _function_

> Least common multiple of ``a`` and ``b`` (``0`` if either is zero).

## numdiff

### `gradient(f, x, rel_step=1e-06)`  _function_

> Central-difference gradient of a scalar ``f`` at the vector ``x``.
>
> Returns a list of ``df/dx_i``. Each partial uses ``(f(x+h e_i) - f(x-h e_i))
> / (2 h)`` with a per-coordinate step scaled by the argument magnitude.

### `hessian(f, x, rel_step=0.0001)`  _function_

> Central-difference Hessian of a scalar ``f`` at ``x`` (symmetric).
>
> Diagonal terms use the second-difference stencil; off-diagonal terms the
> four-point cross stencil. A larger default step than the gradient keeps the
> second differences well-conditioned. Returns an ``n x n`` list of lists.

### `jacobian(f, x, rel_step=1e-06)`  _function_

> Central-difference Jacobian of a vector function ``f: R^n -> R^m``.
>
> Returns an ``m x n`` list of lists ``df_k/dx_i``.

## numeral

### `from_base(s, base)`  _function_

> Parse a digit string ``s`` in ``base`` (2..36) back to a non-negative integer.
>
> Case-insensitive; raises on a digit that is out of range for the base.

### `from_roman(s)`  _function_

> Convert a Roman-numeral string to an integer (subtractive notation).
>
> Case-insensitive. Raises on characters that are not Roman digits.

### `to_base(n, base)`  _function_

> Represent a non-negative integer ``n`` in ``base`` (2..36) as a digit string.
>
> Digits above 9 use lowercase letters (``a`` = 10 .. ``z`` = 35). ``to_base(0, b)`` is
> ``"0"``.

### `to_roman(n)`  _function_

> Convert an integer in ``1..3999`` to its Roman-numeral string (subtractive form).

## ode

### `rk4(f, t0, y0, t1, n=100)`  _function_

> Fixed-step RK4 from ``t0`` to ``t1`` in ``n`` steps.
>
> ``f(t, y)`` returns the derivative (scalar or list matching ``y0``). Returns
> ``(ts, ys)``: the ``n+1`` time points and the state at each (each state a list).

### `rk45(f, t0, y0, t1, tol=1e-08, h0=None, max_steps=100000)`  _function_

> Adaptive Dormand-Prince (RK45) integration from ``t0`` to ``t1``.
>
> Controls the step to keep the estimated local error near ``tol``. Returns
> ``(ts, ys)`` at the accepted steps (non-uniform). ``h0`` is the initial step
> (defaults to a fraction of the interval).

## ols

### `ols_fit(X, y, add_intercept=True, confidence=0.95)`  _function_

> Fit an OLS regression and return coefficients with diagnostics.
>
> Parameters
> ----------
> X : list[list[float]]
>     Design matrix, ``n`` rows of ``k`` regressors (no intercept column unless
>     ``add_intercept=False`` and you supply your own).
> y : list[float]
>     Response vector of length ``n``.
> add_intercept : bool
>     Prepend a column of ones (the default).
> confidence : float
>     Confidence level for the coefficient intervals (default 0.95).
>
> Returns
> -------
> dict
>     ``coefficients`` (intercept first if added), ``std_errors``, ``t_stats``,
>     ``p_values`` (two-sided, per coefficient), ``conf_int`` (list of
>     ``[low, high]`` at ``confidence``), ``r_squared``, ``adj_r_squared``,
>     ``f_stat``, ``f_pvalue`` (overall significance), ``residuals``, ``n_obs``,
>     ``df_resid``.

## ols_hac

### `newey_west(X, y, lags, add_intercept=True)`  _function_

> Newey-West HAC OLS standard errors (heteroskedasticity + autocorrelation).
>
> Adds Bartlett-weighted cross-products of the score vectors ``e_t x_t`` out to
> ``lags`` lags to the White meat, so the covariance is consistent under both
> heteroskedasticity and serial correlation. Returns the same dict shape as
> :func:`white_hc0`. ``lags = 0`` reduces exactly to White (HC0).

### `white_hc0(X, y, add_intercept=True)`  _function_

> White (HC0) heteroskedasticity-consistent OLS standard errors.
>
> Returns a dict with ``coefficients``, robust ``std_errors``, ``t_stats`` and the
> full ``cov`` matrix. Valid when errors are heteroskedastic but not autocorrelated.

## online_cov

### `RunningCovariance(xs=None, ys=None)`  _class_

> Streaming covariance/correlation of paired observations (Welford co-moment).
>
> ``update(x, y)`` folds in one pair; ``covariance()`` and ``correlation()`` read the
> current estimate at any time. Sample (``ddof=1``) covariance by default. Two
> accumulators combine with ``+`` using the parallel-merge formula.

## online_regression

### `RunningRegression(xs=None, ys=None)`  _class_

> Streaming OLS fit of ``y`` on ``x`` (Welford co-moment, constant memory).
>
> ``update(x, y)`` folds in one pair; ``slope()``, ``intercept()``, ``correlation()``,
> ``r_squared()``, and ``predict(x)`` read the current fit at any time. Needs at least
> two points with variation in ``x`` before a slope is defined. Two accumulators combine
> with ``+`` via the parallel-merge formula, giving the exact same fit as folding every
> pair into one.

## optimize

### `nelder_mead(f: Callable[[List[float]], float], x0: Sequence[float], step: float = 0.1, max_iter: int = 2000, tol: float = 1e-10, alpha: float = 1.0, gamma: float = 2.0, rho: float = 0.5, sigma: float = 0.5)`  _function_

> Minimize ``f`` over R^n from ``x0``. Returns (best_x, best_f).
>
> A textbook Nelder-Mead: build an initial simplex by perturbing each
> coordinate, then reflect/expand/contract/shrink until the spread of
> function values falls below ``tol`` or ``max_iter`` is hit.

## order_statistic_tree

### `OrderStatisticTree(universe)`  _class_

> Dynamic multiset with rank/select over a fixed sorted value universe.
>
> Construct with the universe of possible values (deduplicated and sorted internally).
> ``add``/``remove`` adjust multiplicities; ``rank(x)`` counts stored values ``< x``;
> ``select(k)`` returns the ``k``-th smallest (0-indexed); ``count_less``/``count_range``
> answer threshold and interval counts. All queries are ``O(log U)``.

## orthogonal_polys

### `chebyshev_t(n, x)`  _function_

> Chebyshev polynomial of the first kind ``T_n(x)`` (weight ``1/sqrt(1-x^2)``).
>
> ``T_0 = 1``, ``T_1 = x``, ``T_{k+1} = 2x T_k - T_{k-1}``. On ``[-1, 1]``,
> ``T_n(cos theta) = cos(n theta)``, so ``|T_n| <= 1`` there.

### `chebyshev_u(n, x)`  _function_

> Chebyshev polynomial of the second kind ``U_n(x)`` (weight ``sqrt(1-x^2)``).
>
> ``U_0 = 1``, ``U_1 = 2x``, ``U_{k+1} = 2x U_k - U_{k-1}``. On ``[-1, 1]``,
> ``U_n(cos theta) = sin((n+1) theta) / sin(theta)``.

### `hermite_h(n, x)`  _function_

> Physicists' Hermite polynomial ``H_n(x)`` (weight ``e^{-x^2}`` on the real line).
>
> ``H_0 = 1``, ``H_1 = 2x``, ``H_{k+1} = 2x H_k - 2k H_{k-1}``. ``H_2 = 4x^2 - 2``.

### `hermite_he(n, x)`  _function_

> Probabilists' Hermite polynomial ``He_n(x)`` (weight ``e^{-x^2/2}``).
>
> ``He_0 = 1``, ``He_1 = x``, ``He_{k+1} = x He_k - k He_{k-1}``. Related to the physicists'
> form by ``He_n(x) = 2^{-n/2} H_n(x / sqrt 2)``.

### `laguerre_l(n, x, alpha=0.0)`  _function_

> Generalized Laguerre polynomial ``L_n^{(alpha)}(x)`` (weight ``x^alpha e^{-x}``).
>
> ``L_0 = 1``, ``L_1 = 1 + alpha - x``, and
> ``(k+1) L_{k+1} = (2k+1+alpha-x) L_k - (k+alpha) L_{k-1}``. Default ``alpha = 0`` gives the
> ordinary Laguerre polynomials.

## ou_fit

### `fit_ornstein_uhlenbeck(x, dt=1.0)`  _function_

> Estimate OU parameters ``(kappa, theta, sigma)`` from a sampled path.
>
> Parameters
> ----------
> x : sequence of float
>     Observations sampled at uniform spacing ``dt``.
> dt : float
>     Time between observations (in the same units as ``kappa`` is desired).
>
> Returns
> -------
> dict
>     ``{"kappa", "theta", "sigma", "half_life"}``. ``kappa`` is the
>     mean-reversion speed, ``theta`` the long-run mean, ``sigma`` the
>     instantaneous volatility, and ``half_life = ln(2)/kappa``. Raises if the
>     series is not mean-reverting (fitted ``b`` outside ``(0, 1)``).

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

## pade

### `lentz_continued_fraction(a, b, tol=1e-15, max_iter=1000, tiny=1e-300)`  _function_

> Evaluate a continued fraction by the modified Lentz algorithm.
>
> Computes ``b0 + a1/(b1 + a2/(b2 + ...))`` where ``a(k)`` and ``b(k)`` are callables
> giving the ``k``-th partial numerator and denominator (``b(0)`` is the leading
> term, ``a(0)`` is unused). Returns the converged value. Robust to zero
> intermediate values via the ``tiny`` guard.

### `pade(coeffs, m, n)`  _function_

> Pade ``[m/n]`` approximant from Taylor coefficients ``coeffs``.
>
> ``coeffs[k]`` is the coefficient of ``x^k``; needs at least ``m + n + 1`` of them.
> Returns ``(num, den)`` -- the numerator (length ``m+1``) and denominator (length
> ``n+1``, normalized to ``den[0] = 1``) coefficient lists. The approximant's Taylor
> series matches ``coeffs`` through order ``m + n``.

### `pade_eval(num, den, x)`  _function_

> Evaluate a Pade approximant ``(num, den)`` at ``x`` by Horner's method.

## pairs

### `ou_half_life(spread)`  _function_

> Mean-reversion half-life from an AR(1) fit to the spread.
>
> Regresses ``delta_t = a + b * spread_{t-1}`` (the discretized OU); the
> mean-reversion speed is ``kappa = -b`` and the half-life is ``ln(2)/kappa``.
> Positive and finite only for a mean-reverting (``-1 < b < 0``) spread; raises
> otherwise.

### `pairs_hedge_ratio(y, x)`  _function_

> OLS hedge ratio (slope) of ``y`` on ``x`` through the mean.
>
> ``beta = cov(x, y) / var(x)`` -- the number of units of ``x`` to short against
> one unit of ``y`` so the spread ``y - beta x`` is mean-reverting. Recovers the
> true beta on a linear relationship.

### `spread_series(y, x, beta=None)`  _function_

> Spread ``y - beta x`` (hedge-ratio residual).
>
> ``beta`` defaults to the :func:`hedge_ratio`. The series a pairs trade bets
> reverts to its mean.

### `spread_zscore(spread, window=None)`  _function_

> Latest spread z-score against its mean and std (full history or a window).
>
> ``(spread[-1] - mean) / std`` over the last ``window`` points (all if ``None``).
> The pairs-trade entry signal: large magnitude means the spread is stretched.

## panjer

### `aggregate_mean(g)`  _function_

> Mean of an aggregate distribution ``g`` (grid units).

### `aggregate_tvar(g, confidence=0.99)`  _function_

> Tail Value-at-Risk (CTE / expected shortfall) of a discrete aggregate loss.
>
> The probability-weighted average loss in the tail beyond the VaR level,
> ``E[S | S >= VaR]`` computed on the grid. Always at least the VaR, equals the
> mean at ``confidence -> 0``, and is non-decreasing in ``confidence``.

### `aggregate_var(g, confidence=0.99)`  _function_

> Value-at-Risk of a discrete aggregate distribution ``g`` (grid units).
>
> The smallest grid point ``k`` with ``P(S <= k) >= confidence`` -- the loss
> quantile. Non-decreasing in ``confidence``.

### `layer_expected_loss(g, attachment, limit)`  _function_

> Expected loss to a reinsurance layer ``[attachment, attachment + limit]``.
>
> ``E[min(max(S - attachment, 0), limit)]`` -- the excess-of-loss layer cost.

### `panjer_negative_binomial(size, prob, severity_pmf, max_k=None)`  _function_

> Aggregate-loss distribution for compound negative-binomial claim counts.
>
> The claim count ``N ~ NegBinom(size=r, prob=p)`` (``P(N=n) = C(n+r-1, n)
> p^r (1-p)^n``, mean ``r(1-p)/p``) is over-dispersed relative to Poisson
> (variance > mean), capturing claim contagion. It is a Panjer ``(a, b)`` class
> with ``a = 1 - p`` and ``b = (r - 1)(1 - p)``:
>
>     g_0 = p^r  (if severity has no mass at 0),
>     g_k = 1/(1 - a f_0) * sum_{j=1}^{k} (a + b j / k) f_j g_{k-j}.
>
> Parameters
> ----------
> size : float
>     The NB ``r`` (number of failures); ``r > 0``.
> prob : float
>     The NB success probability ``p`` in ``(0, 1]``.
> severity_pmf : sequence of float
>     Severity probabilities on an integer grid.
> max_k : int, optional
>     Aggregate grid cutoff.
>
> Returns
> -------
> list[float]
>     ``g[k] = P(S = k)``; sums to ~1.

### `panjer_poisson(lam, severity_pmf, max_k=None)`  _function_

> Aggregate-loss distribution for compound Poisson via Panjer recursion.
>
> Parameters
> ----------
> lam : float
>     Poisson claim frequency (mean number of claims).
> severity_pmf : sequence of float
>     Severity probabilities on an integer grid ``0, 1, 2, ...`` (index =
>     severity in grid units); should sum to 1.
> max_k : int, optional
>     Highest aggregate grid point to compute. Defaults to a cutoff capturing
>     essentially all mass (``ceil(lam * n) * 4 + 20``).
>
> Returns
> -------
> list[float]
>     ``g[k] = P(S = k)`` on the aggregate grid; sums to ~1.

### `stop_loss_premium(g, retention)`  _function_

> Stop-loss premium ``E[max(S - retention, 0)]`` from the aggregate grid.

## par_yield

### `par_bond_price(discount, coupon_rate, maturity, freq=1, face=100.0)`  _function_

> Price of a coupon bond off the discount curve (for the par-yield check).
>
> Sums the discounted coupons plus the discounted principal. Equals ``face``
> exactly when ``coupon_rate`` is the :func:`par_yield`.

### `par_yield(discount, maturity, freq=1)`  _function_

> Par coupon rate for ``maturity`` from a discount curve ``discount(t)``.
>
> Parameters
> ----------
> discount : callable
>     Discount factor ``P(t)`` (``P(0) = 1``, decreasing for positive rates).
> maturity : float
>     Bond maturity in years.
> freq : int
>     Coupon payments per year.
>
> Returns
> -------
> float
>     The annualized par coupon rate. On a flat curve it equals the flat rate;
>     a bond bearing this coupon prices to exactly par.

## partial_corr

### `partial_correlation(x, y, controls)`  _function_

> Partial correlation of ``x`` and ``y`` controlling for ``controls``.
>
> ``controls`` is a single control column (list) or a list of control columns. Removes
> the linear effect of the controls from both ``x`` and ``y`` (OLS) and correlates the
> residuals. Returns a coefficient in ``[-1, 1]``; a correlation that is purely due to
> the controls drops toward 0.

### `semipartial_correlation(x, y, controls)`  _function_

> Semi-partial (part) correlation: control the ``controls`` out of ``y`` only.
>
> Correlates raw ``x`` with the residual of ``y`` after regressing out the controls --
> the unique contribution of ``x`` to ``y`` beyond the controls. Returns a coefficient
> in ``[-1, 1]``.

## particle_filter

### `particle_filter(observations, transition, log_likelihood, init_sampler, n_particles=1000, seed=12345, resample_threshold=0.5)`  _function_

> Bootstrap particle filter over ``observations``.
>
> Callbacks (all receive the PCG32 ``rng`` where randomness is needed):
>
> * ``init_sampler(rng)`` -> an initial state (any object the other callbacks understand).
> * ``transition(state, rng)`` -> the next state, sampled from the process model.
> * ``log_likelihood(observation, state)`` -> log ``p(obs | state)``.
>
> Returns a dict with ``means`` (the weighted-mean state estimate at each step; states must
> support scalar or per-component averaging -- floats or equal-length lists), ``ess`` (the
> effective sample size at each step), and ``n_resample`` (how many steps resampled).

### `pcg_gaussian(rng, mu=0.0, sigma=1.0)`  _function_

> Draw a normal sample from a PCG32 stream via the Box-Muller transform.
>
> Convenience for writing ``transition``/``init_sampler`` callbacks without pulling in another
> RNG. Uses two uniforms; returns a single normal deviate.

## passing_bablok

### `passing_bablok_regression(x, y)`  _function_

> Passing-Bablok regression of ``y`` on ``x``.
>
> Returns ``(slope, intercept)`` for ``y = slope * x + intercept``. Robust to
> outliers in either variable and symmetric in ``x`` and ``y`` (up to reciprocal
> slope), the standard nonparametric alternative to Deming regression.

## pca

### `jacobi_eigen(matrix, tol=1e-12, max_sweeps=100)`  _function_

> Eigenvalues and eigenvectors of a symmetric matrix (Jacobi rotations).
>
> Returns ``(eigenvalues, eigenvectors)`` where ``eigenvectors[i]`` is the
> orthonormal eigenvector for ``eigenvalues[i]``, sorted by descending
> eigenvalue. Requires a symmetric input; iteratively zeroes off-diagonal
> entries with plane rotations.

### `pca(covariance)`  _function_

> PCA of a covariance matrix: sorted variances, loadings, variance explained.
>
> Returns a dict with ``variances`` (eigenvalues, descending), ``loadings``
> (orthonormal eigenvectors), ``explained`` (each variance over the total), and
> ``cumulative_explained``. For a yield-curve covariance the first three
> components are the level, slope, and curvature factors.

### `pca_scenario(component_index, n_sigma, variances, loadings)`  _function_

> A stress scenario shocking one principal component by ``n_sigma`` std devs.
>
> Returns the vector move ``n_sigma * sqrt(variance_i) * loading_i`` -- a
> ``n_sigma``-standard-deviation move along principal component
> ``component_index``. For a yield curve, component 0 is a parallel (level)
> shift, 1 a slope twist, 2 a curvature bend.

### `project(data_row, loadings, k=None)`  _function_

> Project a data vector onto the first ``k`` principal components (scores).
>
> ``sum_j data_row_j loadings_i_j`` for each retained component ``i``. ``k``
> defaults to all components. The scores are the coordinates of the observation
> in the principal-component basis.

### `reconstruct_covariance(variances, loadings, k=None)`  _function_

> Rebuild a covariance matrix from the top ``k`` principal components.
>
> ``sum_i variance_i * (loading_i outer loading_i)`` over the first ``k``
> components. With all components it reproduces the original covariance exactly
> (spectral decomposition); with ``k`` below the rank it is the best rank-``k``
> approximation.

## pcg

### `PCG32(seed=1234567, seq=54)`  _class_

> PCG-XSH-RR 32-bit generator (O'Neill). Strong, small, reproducible.
>
> ``next_uint32()`` yields a 32-bit output; ``random()`` a float in ``[0, 1)``;
> ``randint(lo, hi)`` an integer in ``[lo, hi]`` (inclusive) without modulo bias.
> Seeded by ``seed`` and an optional stream ``seq``.

### `Xorshift128Plus(seed=1234567)`  _class_

> xorshift128+ generator (Vigna). 64-bit output, long period, reproducible.

## pcr

### `principal_components_regression(X, y, n_components=None)`  _function_

> Principal components regression of ``y`` on the columns of ``X``.
>
> Centers ``X`` and ``y``, runs PCA on the predictor covariance, keeps the top
> ``n_components`` (default: all), regresses on the component scores, and maps the
> coefficients back to the original variables. Returns a dict with ``coefficients``
> (per original predictor), ``intercept``, ``n_components`` and the
> ``explained_variance`` (cumulative fraction retained). Retaining all components
> reproduces OLS.

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

### `burke_ratio(returns: Sequence[float], risk_free=0.0, periods_per_year=252) -> float`  _function_

> Burke ratio: annualized excess return over the root-sum-of-squared drawdowns.
>
> ``ann_excess / sqrt(sum(drawdown_t^2))`` -- penalizes a few deep drawdowns more
> than many shallow ones (an L2 denominator), unlike the L1 :func:`sterling_ratio`.
> Higher is better; raises when there is no drawdown.

### `calmar_ratio(returns: Sequence[float], periods_per_year=252) -> float`  _function_

> Calmar ratio: annualized return divided by the maximum drawdown.
>
> Annualized return is the geometric ``(prod(1+r))^{periods_per_year/n} - 1``.
> Raises if there is no drawdown (undefined ratio).

### `cornish_fisher_expected_shortfall(returns, confidence=0.95, horizon=1.0, n_steps=2000) -> float`  _function_

> Cornish-Fisher (skew/kurtosis-adjusted) expected shortfall, a positive loss.
>
> The average loss in the worst ``1 - confidence`` of the distribution when the
> quantile is the Cornish-Fisher expansion :func:`cornish_fisher_var` uses. The
> standardized shortfall is the tail mean of the expanded quantile,
>
>     ES_z = (1/(1-c)) integral_0^{1-c} z_cf(Phi^{-1}(p)) dp,
>
> computed by midpoint quadrature, then scaled to the loss
> ``ES = -(mean*horizon + ES_z*sigma*sqrt(horizon))``. For a normal series it
> reduces to the Gaussian expected shortfall; negative skew and fat tails push it
> above both the Gaussian ES and the Cornish-Fisher VaR. Always at least the
> Cornish-Fisher VaR.

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

### `deflated_sharpe_ratio(returns, n_trials, sr_variance=None)`  _function_

> Deflated Sharpe ratio (Bailey-López de Prado): PSR against a trials-adjusted benchmark.
>
> When many strategy variants are tested, the best in-sample Sharpe is inflated by
> selection. The DSR is the :func:`probabilistic_sharpe_ratio` evaluated against a
> benchmark equal to the *expected maximum* of ``n_trials`` independent Sharpe
> estimates with cross-trial variance ``sr_variance``:
>
>     SR* = sqrt(sr_variance) * ((1 - gamma) Phi^{-1}(1 - 1/N)
>           + gamma Phi^{-1}(1 - 1/(N e)))
>
> (``gamma`` the Euler-Mascheroni constant). Lower than the plain PSR for
> ``n_trials > 1``, and falling as more trials are tested. ``sr_variance`` defaults
> to the sampling variance ``1/(n-1)`` of a single per-period Sharpe estimate.

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

### `gain_to_pain_ratio(returns: Sequence[float]) -> float`  _function_

> Gain-to-pain ratio: sum of returns over the sum of the absolute losses.
>
> ``sum(r) / sum(|r| for r < 0)`` (Schwager). A scale-free profitability-vs-pain
> measure -- above 1 means net gains exceed the total loss magnitude. Returns
> ``inf`` when there are no losing periods; raises on an empty series.

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

### `jensens_alpha(returns, market_returns, risk_free=0.0, periods_per_year=252) -> float`  _function_

> Jensen's alpha: annualized CAPM-risk-adjusted excess return.
>
> ``alpha = ann(r - rf) - beta * ann(market - rf)`` -- the intercept of the CAPM
> regression, the return earned beyond what the market beta explains. Positive
> alpha is outperformance. Uses the ordinary :func:`market_beta`.

### `longest_drawdown_duration(returns: Sequence[float]) -> int`  _function_

> Longest run of consecutive underwater periods (below a prior peak).
>
> Counts the maximum number of periods between a peak and the point the equity
> curve first recovers to (or exceeds) it. A series that never falls below its
> running peak returns 0.

### `m_squared(returns, market_returns, risk_free=0.0, periods_per_year=252) -> float`  _function_

> Modigliani M-squared: the portfolio's return rescaled to the market's risk.
>
> Levers/de-levers the portfolio (with the risk-free asset) to match the market's
> volatility, then reports the resulting annualized return -- a Sharpe-based
> measure in return units directly comparable to the market. Equals
> ``rf_ann + sharpe * market_vol``. Raises on zero portfolio variance.

### `market_beta(asset_returns, market_returns) -> float`  _function_

> Ordinary CAPM beta: ``Cov(asset, market) / Var(market)`` over all periods.
>
> The slope of the asset's returns regressed on the market's. Equal-length series
> of at least two points; raises on zero market variance.

### `max_drawdown(returns: Sequence[float]) -> float`  _function_

> Maximum peak-to-trough drawdown of the cumulative-return curve.
>
> Compounds the periodic returns into an equity curve and returns the largest
> fractional drop from a running peak, as a non-negative number (0.2 = a 20%
> drawdown). Empty or all-rising series give 0.

### `minimum_track_record_length(returns, benchmark_sr=0.0, confidence=0.95)`  _function_

> Minimum track record length for the Sharpe ratio to beat a benchmark.
>
> The number of observations at which the :func:`probabilistic_sharpe_ratio`
> would reach ``confidence`` that the true SR exceeds ``benchmark_sr``:
>
>     MinTRL = 1 + (1 - skew*SR + (kurt-1)/4 SR^2) (z_conf / (SR - SR*))^2.
>
> Requires the observed per-period Sharpe to exceed the benchmark. Longer for a
> smaller edge or a more skewed/fat-tailed series.

### `omega_ratio(returns: Sequence[float], threshold=0.0) -> float`  _function_

> Omega ratio: probability-weighted gains over losses about a threshold.
>
> ``sum(max(r - threshold, 0)) / sum(max(threshold - r, 0))`` -- the ratio of
> upside to downside area relative to ``threshold``. Values above 1 mean more
> gain mass than loss mass. Returns ``inf`` when there is no downside; raises
> if there is neither upside nor downside.

### `pain_index(returns: Sequence[float]) -> float`  _function_

> Pain index: the average depth of the underwater drawdown curve.
>
> ``mean(drawdown_t)`` -- the mean fractional distance below the running peak.
> A gentler (L1) cousin of the :func:`ulcer_index` (L2).

### `pain_ratio(returns: Sequence[float], risk_free=0.0, periods_per_year=252) -> float`  _function_

> Pain ratio: annualized excess return over the :func:`pain_index`.
>
> The L1 analogue of the :func:`ulcer_performance_index`. Higher is better;
> raises when there is no drawdown.

### `probabilistic_sharpe_ratio(returns, benchmark_sr=0.0)`  _function_

> Probabilistic Sharpe ratio (Bailey-López de Prado).
>
> The probability that the true per-period Sharpe ratio exceeds a ``benchmark_sr``
> (also per period), correcting the estimator's standard error for the sample's
> skewness and (excess) kurtosis and the sample length ``n``:
>
>     PSR = Phi( (SR - SR*) sqrt(n - 1)
>                / sqrt(1 - skew*SR + (kurt-1)/4 * SR^2) ),
>
> with ``SR`` the per-period Sharpe. Above 0.5 when the observed SR beats the
> benchmark; rises with a longer, less-skewed, thinner-tailed track record.

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

### `sterling_ratio(returns: Sequence[float], risk_free=0.0, periods_per_year=252, excess=0.1) -> float`  _function_

> Sterling ratio: annualized excess return over the average drawdown plus a margin.
>
> ``ann_excess / (average_drawdown + excess)`` with the classic ``excess = 10%``
> margin that keeps the denominator from collapsing on shallow-drawdown series.
> Higher is better; the drawdown is the mean underwater depth (the pain index).

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

### `treynor_ratio(returns, market_returns, risk_free=0.0, periods_per_year=252) -> float`  _function_

> Treynor ratio: annualized excess return per unit of market beta.
>
> ``ann_excess / beta`` where ``beta`` is the CAPM :func:`market_beta`. Like Sharpe
> but dividing by systematic (non-diversifiable) risk instead of total volatility,
> so it rewards return per unit of market exposure. Raises for a non-positive beta.

### `ulcer_index(returns: Sequence[float]) -> float`  _function_

> Ulcer index: RMS of the underwater drawdown curve.
>
> ``sqrt(mean(drawdown_t^2))`` over the :func:`drawdown_curve` -- a downside risk
> measure that penalizes deep and prolonged drawdowns more than shallow ones,
> unlike volatility which treats up and down moves alike. Zero for a series that
> never draws down.

### `ulcer_performance_index(returns: Sequence[float], risk_free=0.0, periods_per_year=252) -> float`  _function_

> Ulcer performance index (Martin ratio): excess return over the Ulcer index.
>
> ``(annualized_excess_return) / ulcer_index`` -- a return-per-unit-of-drawdown-
> pain ratio, the drawdown analogue of the Sharpe ratio. Higher is better;
> raises if there is no drawdown (infinite ratio).

### `up_capture(returns, benchmark_returns) -> float`  _function_

> Up-capture ratio: the asset's geometric return in up-benchmark periods
> divided by the benchmark's. Above 1 means the asset outpaces the benchmark
> in rising markets.

## permutation_test

### `paired_permutation_test(x, y, n_permutations=9999, alternative='two-sided', seed=1234567)`  _function_

> Paired permutation test on the within-pair differences ``x_i - y_i``.
>
> Under the null the sign of each difference is exchangeable, so each permutation
> flips signs at random. The statistic is the mean difference. Returns the same dict
> shape as :func:`permutation_test`.

### `permutation_test(a, b, statistic=None, n_permutations=9999, alternative='two-sided', seed=1234567)`  _function_

> Two-sample permutation test.
>
> ``statistic(a, b)`` defaults to the difference in means. Pools the two samples,
> reshuffles the labels ``n_permutations`` times, and returns a dict with the
> ``observed`` statistic and the ``p_value``. ``alternative`` is ``"two-sided"``,
> ``"greater"`` or ``"less"``. Deterministic for a fixed ``seed``.

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

## platt_scaling

### `platt_calibrate(train_scores, train_labels, new_scores=None, max_iter=100)`  _function_

> Fit Platt scaling and return ``(params, predict)``.
>
> ``params`` is ``(A, B)`` and ``predict(scores)`` maps raw scores to calibrated
> probabilities. If ``new_scores`` is given, returns ``(params, calibrated_new)``
> instead, applying the fit directly.

### `platt_fit(scores, labels, max_iter=100, tol=1e-10)`  _function_

> Fit Platt sigmoid parameters ``(A, B)`` to scores and binary labels.
>
> ``P(y=1|s) = 1 / (1 + exp(-(A s + B)))``. Uses Platt's smoothed targets and
> Newton's method on the regularized logistic loss. Returns ``(A, B)``; a
> well-separated classifier gives ``A > 0`` (probability rises with the score) in
> this logistic-standard sign convention. Pure standard library.

### `platt_predict(scores, A, B)`  _function_

> Apply a fitted Platt sigmoid to scores, returning calibrated probabilities.

## poisson2d

### `poisson2d(f_grid, boundary, dx, dy, omega=1.5, tol=1e-08, max_iter=10000)`  _function_

> Solve ``u_xx + u_yy = f`` on a grid by SOR with Dirichlet boundaries.
>
> ``f_grid`` is the source term ``f`` as a 2-D list (rows = y, cols = x);
> ``boundary`` a same-shaped grid whose *edge* values set the fixed boundary (interior
> entries are the initial guess). ``dx``/``dy`` are the grid spacings, ``omega`` the
> SOR factor in ``(0, 2)`` (1 = Gauss-Seidel). Returns ``(u, n_iter)`` -- the solution
> grid and the sweeps taken. Laplace is the ``f_grid`` all-zero case.

## poisson_regression

### `poisson_predict(model, X_query, add_intercept=True)`  _function_

> Predicted rates ``exp(x' beta)`` for rows ``X_query``.

### `poisson_regression(X, y, add_intercept=True, max_iter=50, tol=1e-08)`  _function_

> Fit a Poisson GLM (log link) by IRLS.
>
> ``y`` are non-negative counts. Returns a dict with ``coefficients`` (intercept first
> if added), ``n_iter`` and the ``log_likelihood`` at convergence. Predicted rate for
> a row is ``exp(x' beta)``. Uses Fisher-scoring IRLS: working response
> ``z = eta + (y - mu)/mu`` with weights ``mu``.

## poly_features

### `polynomial_features(X, degree=2, interaction_only=False, include_bias=True)`  _function_

> Expand rows of ``X`` into polynomial / interaction features.
>
> Parameters
> ----------
> X : list[list[float]]
>     Input matrix, ``n`` rows of ``p`` features.
> degree : int
>     Maximum total degree of the monomials (>= 1).
> interaction_only : bool
>     If True, exclude pure powers (x_i^2, ...), keeping only products of
>     distinct variables.
> include_bias : bool
>     Prepend a constant 1 column.
>
> Returns
> -------
> (features, powers) : (list[list[float]], list[tuple])
>     The expanded matrix and, for each column, the tuple of input indices whose
>     product forms it (``()`` for the bias). With ``degree=1`` and a bias this
>     is the original matrix with a leading ones column.

## poly_mod

### `lagrange_interpolate_mod(points, mod)`  _function_

> Recover the polynomial through ``points`` = ``[(x_i, y_i)]`` modulo prime ``mod``.
>
> Returns coefficients (lowest-degree first) of the unique polynomial of degree ``< n``
> that passes through all ``n`` points; ``x_i`` must be distinct modulo ``mod``. Builds
> ``sum_i y_i * prod_{j!=i} (x - x_j)/(x_i - x_j)`` with exact modular inverses.

### `poly_add_mod(a, b, mod)`  _function_

> Sum of two polynomials (lowest-degree-first) modulo ``mod``.

### `poly_eval_mod(coeffs, x, mod)`  _function_

> Evaluate a polynomial (lowest-degree-first) at ``x`` modulo ``mod`` (Horner).

### `poly_mul_mod(a, b, mod)`  _function_

> Product of two polynomials (lowest-degree-first) modulo ``mod``.

## polyinterp

### `divided_differences(xs, ys)`  _function_

> Newton divided-difference coefficients for the points ``(xs, ys)``.
>
> Returns the list ``[f[x0], f[x0,x1], ...]`` -- the leading coefficients of the
> Newton form, computed in ``O(n^2)``. Feed these to :func:`newton_polynomial`.

### `neville(xs, ys, x)`  _function_

> Evaluate the interpolating polynomial at ``x`` by Neville's algorithm.
>
> Returns ``(value, error_estimate)`` where the error estimate is the magnitude of
> the last correction -- a practical indicator of interpolation accuracy and the
> basis of Richardson extrapolation (interpolating to ``x = 0`` in the step size).

### `newton_polynomial(xs, coef, x)`  _function_

> Evaluate the Newton form with divided-difference ``coef`` at ``x`` (Horner).

## polylog

### `dilog(x)`  _function_

> Dilogarithm (Spence's function) ``Li_2(x) = sum_{k>=1} x^k / k^2`` for real ``x <= 1``.
>
> Reflection and inversion identities fold ``x`` into ``[-1, 1/2]`` where the series
> converges quickly. ``Li_2(1) = pi^2/6``, ``Li_2(-1) = -pi^2/12``, and
> ``Li_2(1/2) = pi^2/12 - (ln 2)^2/2``. Raises for ``x > 1`` (there ``Li_2`` is complex).

### `polylog(s, z, tol=1e-15, max_terms=200000)`  _function_

> Polylogarithm ``Li_s(z) = sum_{k>=1} z^k / k^s`` for real ``s`` and real ``|z| <= 1``.
>
> Uses direct summation, which converges geometrically for ``|z| < 1``. The endpoint
> ``z = 1`` returns ``zeta(s)`` (requires ``s > 1``); ``z = -1`` returns ``-eta(s)``.
> Convergence slows as ``z -> 1`` with small ``s``; ``max_terms`` bounds the work.

## polynomial

### `poly_add(a, b)`  _function_

> Sum of two polynomials (coefficient lists, low-degree-first).

### `poly_derivative(c)`  _function_

> Derivative of a polynomial: ``[c1, 2 c2, 3 c3, ...]``.

### `poly_divmod(num, den)`  _function_

> Polynomial long division: return ``(quotient, remainder)``.
>
> ``num = quotient * den + remainder`` with ``deg(remainder) < deg(den)``. Raises on
> a zero divisor.

### `poly_eval(c, x)`  _function_

> Evaluate a polynomial at ``x`` by Horner's method.

### `poly_gcd(a, b, tol=1e-09)`  _function_

> Monic greatest common divisor of two polynomials (Euclidean algorithm).
>
> Returns the GCD normalized to a monic polynomial (leading coefficient 1); useful
> for detecting and factoring out repeated roots (``gcd(p, p')``). Coefficients below
> ``tol`` in the remainder are treated as zero to tame round-off.

### `poly_integral(c, constant=0.0)`  _function_

> Antiderivative of a polynomial, with integration constant ``constant``.

### `poly_mul(a, b)`  _function_

> Product of two polynomials by direct convolution.
>
> Exact for the small polynomials typical of algebra; for long polynomials the
> FFT-based :func:`quantforge.convolve` is asymptotically faster.

### `poly_sub(a, b)`  _function_

> Difference ``a - b`` of two polynomials.

## polyroots

### `polynomial_roots(coeffs, tol=1e-12, max_iter=500)`  _function_

> All roots of a polynomial by the Durand-Kerner method.
>
> ``coeffs`` are the coefficients from the highest degree down (e.g. ``[1, -3, 2]``
> for ``x^2 - 3x + 2``); real or complex. Returns a list of the ``n`` roots as
> complex numbers (a root with a negligible imaginary part is still returned as
> ``complex`` -- take ``.real`` if you know it is real). Leading zeros are trimmed.

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

### `risk_budget_weights(cov, budgets, tol=1e-12, max_iter=2000) -> list`  _function_

> Long-only weights whose percentage risk contributions match ``budgets``.
>
> Generalizes :func:`risk_parity_weights` (which targets equal budgets) to an
> arbitrary risk-budget vector. Solves the fixed point
> ``w_i <- sqrt(budget_i * w_i / (C w)_i)`` renormalized, so at convergence the
> percentage risk contribution of asset ``i`` equals ``budget_i / sum(budgets)``.
> ``budgets`` must be positive; they are normalized internally. Weights are
> positive and sum to one.

### `risk_contributions(weights, cov) -> list`  _function_

> Contribution of each asset to total portfolio variance.
>
> ``RC_i = w_i (C w)_i``. The contributions sum to the portfolio variance
> ``w' C w``; dividing by that sum gives the percentage risk contributions. Equal
> percentage contributions is the risk-parity condition.

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

## power

### `one_sample_z_power(effect_size, n, alpha=0.05)`  _function_

> Power of a two-sided one-sample mean test (z-test, normal approx).
>
> ``effect_size`` is the mean shift in SD units; the noncentrality is
> ``effect_size sqrt(n)``.

### `one_sample_z_sample_size(effect_size, power=0.8, alpha=0.05)`  _function_

> Sample size for a target ``power`` in a one-sample mean test.
>
> ``n = (z_{alpha/2} + z_{beta})^2 / d^2``, rounded up.

### `proportion_power(p1, p2, n_per_group, alpha=0.05)`  _function_

> Power of a two-sided two-proportion test (normal approximation).
>
> Uses the unpooled standard error at the alternative and the pooled standard
> error under the null. ``p1``, ``p2`` are the two success probabilities.

### `proportion_sample_size(p1, p2, power=0.8, alpha=0.05)`  _function_

> Per-group sample size for a target ``power`` in a two-proportion test.
>
> ``n = (z_{alpha/2} sqrt(2 pbar (1-pbar)) + z_{beta} sqrt(p1(1-p1)+p2(1-p2)))^2 /
> (p1 - p2)^2``, rounded up.

### `two_sample_t_power(effect_size, n_per_group, alpha=0.05)`  _function_

> Power of a two-sided two-sample test for a mean difference (normal approx).
>
> ``effect_size`` is Cohen's ``d`` (mean difference in pooled-SD units). With
> ``n_per_group`` observations in each arm the noncentrality is
> ``d sqrt(n / 2)``, and the power is ``Phi(ncp - z_{alpha/2})`` plus the far
> tail. Increases with the effect size, the sample size, and ``alpha``.

### `two_sample_t_sample_size(effect_size, power=0.8, alpha=0.05)`  _function_

> Per-group sample size for a target ``power`` in a two-sample mean test.
>
> Inverts the normal-approximation power: ``n = 2 (z_{alpha/2} + z_{beta})^2 /
> d^2``, rounded up. Raises for a zero effect size (infinite sample).

## power_iteration

### `inverse_iteration(A, mu=0.0, x0=None, tol=1e-12, max_iter=1000)`  _function_

> Eigenpair of ``A`` whose eigenvalue is closest to the shift ``mu``.
>
> Runs power iteration on ``(A - mu I)^{-1}``. With ``mu = 0`` this finds the
> smallest-magnitude eigenvalue; with ``mu`` near a known approximate eigenvalue it
> refines that one. Returns the same dict shape as :func:`power_iteration`.

### `power_iteration(A, x0=None, tol=1e-12, max_iter=1000)`  _function_

> Dominant eigenpair of ``A`` by power iteration.
>
> Returns a dict with ``eigenvalue`` (largest magnitude, via the Rayleigh quotient),
> ``eigenvector`` (unit norm), ``n_iter`` and ``converged``. Converges when the
> dominant eigenvalue is unique in magnitude; the sign convention makes the first
> non-negligible component positive.

### `rayleigh_quotient(A, x)`  _function_

> Rayleigh quotient ``x' A x / x' x`` -- the least-squares eigenvalue for ``x``.

## prefix_sum

### `DifferenceArray(n_or_values)`  _class_

> Offline range-add / final-read via a difference array.
>
> ``add(lo, hi, delta)`` adds ``delta`` to every index in ``[lo, hi)`` in ``O(1)``; after
> all updates, ``result()`` materializes the final array in ``O(n)``. Ideal when many
> range updates precede a single read.

### `PrefixSum1D(values)`  _class_

> Constant-time range sums over a fixed 1-D array.
>
> ``range_sum(lo, hi)`` returns the sum of ``values[lo:hi]`` (half-open) in ``O(1)`` after
> an ``O(n)`` build. ``prefix(i)`` gives the sum of the first ``i`` elements.

### `PrefixSum2D(grid)`  _class_

> Summed-area table for constant-time rectangle sums over a fixed 2-D grid.
>
> Build from a list of equal-length rows. ``range_sum(r0, c0, r1, c1)`` returns the sum
> over rows ``[r0, r1)`` and columns ``[c0, c1)`` (half-open) in ``O(1)``.

## prob_forecast

### `crps_ensemble(actual, ensemble)`  _function_

> Continuous ranked probability score of an ensemble forecast (single target).
>
> Uses the empirical-CDF form ``CRPS = mean|X - y| - 0.5 mean|X - X'|`` where ``X``,
> ``X'`` are independent ensemble members and ``y`` the realized value. Reduces to
> the absolute error for a deterministic (single-member) forecast, and is zero when
> every member equals the target. ``ensemble`` is the list of member forecasts.

### `interval_coverage(actual, lower, upper)`  _function_

> Empirical coverage: fraction of actuals within ``[lower, upper]``.
>
> Should match the interval's nominal level (e.g. ~0.9 for a 90% interval).

### `interval_score(actual, lower, upper, alpha=0.1)`  _function_

> Winkler interval score for central ``1 - alpha`` prediction intervals.
>
> ``S = (upper - lower) + (2/alpha)(lower - a) if a < lower
>                       + (2/alpha)(a - upper) if a > upper``.
> Rewards narrow intervals and penalizes actuals falling outside, scaled so the
> penalty grows as the nominal coverage tightens. Lower is better. Aligned series.

### `pinball_loss(actual, quantile_forecast, tau)`  _function_

> Average pinball (quantile) loss at level ``tau``.
>
> ``L = mean( tau (a - q)      if a >= q
>             (1 - tau)(q - a)  otherwise )``.
> Minimized in expectation when ``quantile_forecast`` is the true ``tau``-quantile
> of the target. Aligned series; ``tau`` in ``(0, 1)``.

## proportion_ci

### `agresti_coull_interval(k, n, confidence=0.95)`  _function_

> Agresti-Coull interval: a Wald interval on ``z^2``-adjusted counts.

### `clopper_pearson_interval(k, n, confidence=0.95)`  _function_

> Exact Clopper-Pearson interval by inverting the binomial CDF.
>
> The lower limit is the ``p`` with ``P(X >= k) = alpha/2`` and the upper limit the
> ``p`` with ``P(X <= k) = alpha/2`` (``alpha = 1 - confidence``); the boundary
> cases ``k = 0`` and ``k = n`` give a one-sided interval. Guaranteed to cover at
> least ``confidence`` of the time -- conservative but never under-covering.

### `wald_interval(k, n, confidence=0.95)`  _function_

> Normal-approximation (Wald) interval, clamped to ``[0, 1]``.

### `wilson_interval(k, n, confidence=0.95)`  _function_

> Wilson score interval -- stays in ``[0, 1]`` with good small-sample coverage.

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

## qr

### `qr_decomposition(A)`  _function_

> Householder QR of an ``m x n`` matrix ``A`` (``m >= n``).
>
> Returns ``(Q, R)`` with ``Q`` an ``m x m`` orthogonal matrix and ``R`` an
> ``m x n`` upper-triangular matrix such that ``A = Q R``. Pure Python lists.

### `qr_solve(A, b)`  _function_

> Least-squares solution of ``A x = b`` via QR (``x`` minimizes ||A x - b||).
>
> Computes ``Q' b`` and back-substitutes the upper-triangular ``R``. Numerically
> stabler than the normal equations for ill-conditioned ``A``. ``A`` is ``m x n``
> with ``m >= n``; returns the length-``n`` coefficient vector.

## quadrature

### `adaptive_simpson(f, a, b, tol=1e-10, max_depth=50)`  _function_

> Adaptive Simpson quadrature with error control to ``tol``.
>
> Recursively bisects where the Simpson estimate has not converged, so it
> concentrates work on the hard parts of the integrand. Returns the integral.

### `clenshaw_curtis(f, a, b, n=64)`  _function_

> Clenshaw-Curtis quadrature: sample at Chebyshev points, weight by the DCT.
>
> Evaluates ``f`` at the ``n + 1`` Chebyshev extrema
> ``x_j = cos(pi j / n)`` mapped to ``[a, b]`` and combines them with the classic
> Clenshaw-Curtis weights (a discrete cosine sum of the even Chebyshev moments
> ``2 / (1 - k^2)``). Like Gauss-Legendre it is spectrally accurate for smooth
> integrands, but the order ``n`` is a free parameter and the nodes nest, so it is
> a convenient high-order rule where the fixed 2-5 point Gauss rule is too coarse.
>
> ``n`` must be a positive even integer (rounded up). Pure standard library.

### `gauss_legendre(f, a, b, n=5)`  _function_

> Fixed-order Gauss-Legendre quadrature (``n`` in {2,3,4,5}).
>
> Maps the reference nodes to ``[a, b]``. Exact for polynomials up to degree
> ``2n - 1`` -- very accurate for smooth integrands with few evaluations.

### `romberg(f, a, b, max_order=10, tol=1e-12)`  _function_

> Romberg integration: Richardson extrapolation on the trapezoid rule.
>
> Builds the Romberg tableau, refining the composite trapezoid estimate by
> successive interval halvings and extrapolating away the Euler-Maclaurin error
> terms. Row ``T[k][0]`` is the ``2^k``-panel trapezoid; each further column
> cancels the next even power of the step,
>
>     T[k][j] = (4^j T[k][j-1] - T[k-1][j-1]) / (4^j - 1),
>
> so ``T[k][k]`` converges as ``O(h^{2k+2})`` for a smooth integrand. Stops early
> when two successive diagonal estimates agree to ``tol``. Ideal for smooth
> integrands where it reaches machine precision in a handful of halvings; for
> endpoint singularities use :func:`tanh_sinh` instead.

### `simpson(f, a, b, n=1000)`  _function_

> Composite Simpson's rule (``n`` even) -- exact for cubics.
>
> Rounds ``n`` up to the next even number. Fourth-order accurate.

### `tanh_sinh(f, a, b, levels=6, h0=1.0)`  _function_

> Tanh-sinh (double-exponential) quadrature over ``[a, b]``.
>
> Substitutes ``x = (a+b)/2 + (b-a)/2 * tanh((pi/2) sinh(t))`` and integrates the
> transformed integrand in ``t`` on a uniform grid. The change of variables makes
> the abscissae cluster double-exponentially toward the endpoints and the weights
> decay super-fast, so the rule converges even when ``f`` has integrable
> endpoint singularities (``1/sqrt(x)``, ``ln x``, ...) where Gauss-Legendre and
> Simpson struggle. ``levels`` successive grid halvings refine the step from
> ``h0``; the number of function evaluations is ``O(2^levels / h0)``.
>
> Nodes are computed off the endpoints, so ``f`` is never evaluated exactly at
> ``a`` or ``b`` -- an integrable singularity at either end is fine.

### `trapezoid(f, a, b, n=1000)`  _function_

> Composite trapezoid rule with ``n`` sub-intervals over ``[a, b]``.

## quantile_regression

### `quantile_regression(X, y, tau=0.5, add_intercept=True, max_iter=200, tol=1e-08)`  _function_

> Fit a linear ``tau``-quantile regression ``y ~ X beta``.
>
> Returns the coefficient list (intercept first if added). ``tau`` in ``(0, 1)``
> selects the conditional quantile: 0.5 is the median (least-absolute-deviations)
> fit, higher ``tau`` tracks the upper conditional tail. Solved by IRLS on the
> asymmetric absolute loss; a small floor keeps the reweighting stable at zero
> residuals.

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

## quaternion

### `axis_angle_to_quat(axis, angle)`  _function_

> Unit quaternion for a rotation of ``angle`` radians about ``axis`` (a 3-vector).
>
> The axis is normalized internally. ``angle = 0`` gives the identity ``(1, 0, 0, 0)``.

### `quat_conjugate(q)`  _function_

> Conjugate ``(w, -x, -y, -z)`` -- the inverse rotation for a unit quaternion.

### `quat_multiply(a, b)`  _function_

> Hamilton product ``a * b`` of two quaternions ``(w, x, y, z)``.
>
> Composition of rotations: ``quat_multiply(a, b)`` applies ``b`` then ``a``. Not
> commutative.

### `quat_normalize(q)`  _function_

> Return ``q`` scaled to unit norm; raises on the zero quaternion.

### `quat_to_axis_angle(q)`  _function_

> Recover ``(axis, angle)`` from a unit quaternion; ``angle`` in ``[0, pi]``.
>
> Returns a unit ``axis`` and the rotation angle in radians. For the identity the axis
> is arbitrary; ``(1, 0, 0)`` is returned.

### `rotate_vector(q, v)`  _function_

> Rotate a 3-vector ``v`` by unit quaternion ``q`` (``q v q*``).

### `slerp(a, b, t)`  _function_

> Spherical linear interpolation between unit quaternions ``a`` and ``b`` at ``t``.
>
> ``t = 0`` returns ``a``, ``t = 1`` returns ``b``, and intermediate ``t`` traces the
> shortest constant-speed arc on the unit sphere -- the standard smooth rotation blend.
> Chooses the shorter path (negates ``b`` if the dot product is negative).

## rainbow_n

### `rainbow_option_mc(spots, strike, t, r, sigmas, corr, q=None, best=True, is_call=True, n_paths=100000, seed=1234567)`  _function_

> Monte Carlo price of an ``n``-asset best-of / worst-of option.
>
> Parameters
> ----------
> spots, sigmas : per-asset spot and volatility (length ``n``).
> strike, t, r : option strike, maturity, risk-free rate.
> corr : ``n x n`` correlation matrix.
> q : optional per-asset dividend yields.
> best : True for best-of (max), False for worst-of (min).
> is_call : call if True, else put.
> n_paths, seed : simulation controls.
>
> Returns
> -------
> float
>     Discounted Monte Carlo option value.

## random_forest

### `fit_random_forest(X, y, n_trees=10, max_depth=5, min_samples=2, seed=1234567)`  _function_

> Fit a random forest of bootstrap-resampled CART trees.
>
> Parameters
> ----------
> X, y : data and labels.
> n_trees : int
>     Number of trees in the ensemble (>= 1).
> max_depth, min_samples : passed to each tree.
> seed : int
>     Seed for the reproducible bootstrap resampling.
>
> Returns
> -------
> dict
>     ``{"trees": [...], "classes": sorted labels}``.

### `predict_random_forest(forest, X_query)`  _function_

> Majority-vote prediction across the forest's trees.

## range_accrual

### `range_accrual_note(S, L, U, t, r, sigma, coupon, observations, b=None, q=0.0, notional=1.0)`  _function_

> Present value of a range-accrual note's coupon leg.
>
> Parameters
> ----------
> S, L, U : float
>     Spot and the lower/upper edges of the accrual band, ``0 < L < U``.
> t : float
>     Maturity in years; the coupon is paid at ``t``.
> r, sigma : float
>     Risk-free rate and volatility.
> coupon : float
>     Full coupon rate earned if the index is in range on every observation.
> observations : int
>     Number of equally spaced observation dates in ``(0, t]``. Date ``i`` of
>     ``m`` falls at ``t_i = t * i / m``.
> b : float, optional
>     Cost of carry / index drift. Defaults to ``r - q``.
> q : float
>     Dividend yield, used only when ``b`` is not given.
> notional : float
>     Note notional.
>
> Returns
> -------
> float
>     Discounted expected coupon. Non-negative, rises with a wider band, and
>     approaches ``notional * coupon * e^{-r t}`` as the band widens to cover
>     the whole positive axis.

## ransac

### `ransac_line(x, y, threshold, n_iterations=200, seed=1234567)`  _function_

> RANSAC line fit; returns a dict with the consensus model.
>
> ``threshold`` is the maximum residual for a point to count as an inlier. Runs
> ``n_iterations`` minimal (2-point) fits, keeps the model with the most inliers, and
> refits least squares on that inlier set. Returns ``slope``, ``intercept``,
> ``inliers`` (index list) and ``n_inliers``. Deterministic for a fixed ``seed``.

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

## realized

### `bipower_variation(returns)`  _function_

> Barndorff-Nielsen-Shephard bipower variation.
>
> ``BV = (pi/2) sum_{i=2}^{n} |r_i| |r_{i-1}|``. Estimates the integrated
> variance of the continuous part only and is robust to jumps. Requires at
> least two returns.

### `jump_variation(returns)`  _function_

> Jump component of the quadratic variation, ``max(RV - BV, 0)``.
>
> Zero (up to sampling noise) for a purely continuous path; strictly positive
> when the returns contain a jump. Clamped at zero because the estimator can go
> slightly negative on jump-free data.

### `med_realized_variance(returns)`  _function_

> MedRV jump-robust integrated-variance estimator (Andersen-Dobrev-Schaumburg).
>
> ``MedRV = c * (n / (n - 2)) * sum_i median(|r_{i-1}|, |r_i|, |r_{i+1}|)^2`` with
> ``c = pi / (6 - 4 sqrt(3) + pi)``. Taking the median of three neighbouring
> magnitudes discards a lone jump and, unlike MinRV, is also robust to two nearby
> jumps and less sensitive to zero returns. Requires at least three returns.

### `min_realized_variance(returns)`  _function_

> MinRV jump-robust integrated-variance estimator (Andersen-Dobrev-Schaumburg).
>
> ``MinRV = (pi / (pi - 2)) * (n / (n - 1)) * sum_i min(|r_i|, |r_{i+1}|)^2``. Each
> term pairs adjacent returns and keeps the smaller magnitude, so an isolated jump
> (which lands in one return) is discarded by the minimum. Converges to the
> integrated variance of the continuous part; more robust to jumps than bipower
> variation and to occasional zero returns. Requires at least two returns.

### `realized_quarticity(returns)`  _function_

> Realized quarticity ``(n / 3) * sum r_i^4``.
>
> A consistent estimator of the integrated quarticity ``integral sigma^4``, which
> sets the asymptotic variance of realized variance and appears in the standard
> errors of realized-volatility jump tests. Non-negative.

### `realized_variance_from_returns(returns)`  _function_

> Realized variance ``sum r_i^2`` of a return series.
>
> Consistent for the total quadratic variation (diffusion plus jumps) as the
> sampling frequency rises.

### `realized_volatility_signature(returns, annualization=1.0)`  _function_

> Realized volatility ``sqrt(RV)``, optionally annualized.
>
> ``annualization`` multiplies the variance before the square root (e.g. the
> number of periods per year for intraday returns aggregated to one day times
> 252). Defaults to 1 (the raw realized vol of the supplied returns).

## realized_kernel

### `realized_kernel(prices, bandwidth=None)`  _function_

> Realized-kernel integrated-variance estimate from a (log) price series.
>
> Uses the flat-top Parzen kernel over ``bandwidth`` lags of the intraday-return
> autocovariances. If ``bandwidth`` is ``None`` it defaults to the
> Barndorff-Nielsen rule of thumb ``H ~ n^{3/5}`` (capped below the number of
> returns). Robust to i.i.d. microstructure noise, unlike the naive realized
> variance, and always non-negative for the Parzen kernel. Requires at least three
> prices.

## rebalance

### `drift_weights(weights, asset_returns)`  _function_

> Buy-and-hold weights after one period of ``asset_returns``.
>
> Each position grows by ``(1 + r_i)``; the new weights are the grown values
> renormalized to sum to one. The starting point for the next rebalance decision.

### `no_trade_band_rebalance(current_weights, target_weights, band)`  _function_

> Rebalance only positions whose drift exceeds a ``band`` tolerance.
>
> Positions within ``band`` of their target are left untouched (no trade);
> those outside are moved to target. The remaining weight from the traded legs
> is left as-is (the untouched legs keep their drifted weight), so the result is
> renormalized to sum to one. Reduces turnover versus a full rebalance.

### `transaction_cost(current_weights, target_weights, cost_bps)`  _function_

> Transaction-cost drag of a rebalance: ``2 * turnover * cost_bps / 1e4``.
>
> Costs the round-trip (both sides) at ``cost_bps`` basis points of the traded
> notional. Zero when no trade is needed.

### `turnover(current_weights, target_weights)`  _function_

> One-way turnover ``0.5 * sum |target - current|`` (fraction of the book).
>
> The fraction of the portfolio traded to move from current to target weights;
> zero when already on target, up to one for a full turnover.

## regression_tree

### `fit_regression_tree(X, y, max_depth=5, min_samples=2)`  _function_

> Fit a CART regression tree.
>
> ``X`` is a list of feature rows, ``y`` the numeric targets. Splits greedily to
> maximize squared-error reduction until ``max_depth`` or ``min_samples`` stops it;
> leaves store the mean target. Returns a nested-dict tree for
> :func:`predict_regression_tree`.

### `predict_regression_tree(tree, X_query)`  _function_

> Predict targets for rows ``X_query`` with a fitted regression tree.

## resample

### `bca_bootstrap_ci(data, statistic=None, n_boot=2000, confidence=0.95, seed=1234567)`  _function_

> Bias-corrected accelerated (BCa) bootstrap confidence interval.
>
> Efron's BCa improves on the percentile method by correcting for median bias
> (``z0``, from the fraction of bootstrap replicates below the point estimate)
> and skewness (``a``, the acceleration from the jackknife). The percentiles are
> shifted:
>
>     alpha1 = Phi(z0 + (z0 + z_lo)/(1 - a(z0 + z_lo)))
>     alpha2 = Phi(z0 + (z0 + z_hi)/(1 - a(z0 + z_hi))).
>
> Reduces to the plain :func:`bootstrap_ci` when ``z0`` and ``a`` are zero
> (symmetric, unbiased statistic). Returns ``(lower, point, upper)``.

### `bootstrap_ci(data, statistic=None, n_boot=2000, confidence=0.95, seed=1234567)`  _function_

> IID bootstrap confidence interval for a sample ``statistic``.
>
> Resamples ``data`` with replacement ``n_boot`` times, applies ``statistic`` to
> each resample, and returns ``(lower, point, upper)`` -- the percentile-method
> interval at ``confidence`` plus the statistic on the original sample. The
> interval brackets the point estimate and narrows as the sample grows.
> ``statistic`` defaults to the sample mean.

### `jackknife_estimate(data, statistic=None)`  _function_

> Delete-one jackknife estimate and standard error of a ``statistic``.
>
> Recomputes the statistic on each leave-one-out subsample. Returns
> ``(estimate, standard_error)`` with the bias-aware jackknife SE
> ``sqrt((n-1)/n * sum (theta_i - theta_bar)^2)``. ``statistic`` defaults to the
> sample mean.

### `moving_block_bootstrap_ci(data, statistic=None, block=10, n_boot=2000, confidence=0.95, seed=1234567)`  _function_

> Moving-block (Kunsch) bootstrap CI for serially-correlated data.
>
> Resamples fixed-length overlapping blocks of length ``block`` from the series
> (wrapping at the end) and concatenates ceil(n / block) of them, truncated to
> ``n``, preserving within-block dependence. Like the stationary bootstrap it
> gives valid intervals for autocorrelated data -- wider than the IID
> :func:`bootstrap_ci` for a positively autocorrelated mean, and it agrees with
> the IID interval when ``block = 1``. ``statistic`` defaults to the sample mean.

### `stationary_bootstrap_ci(data, statistic=None, mean_block=10, n_boot=2000, confidence=0.95, seed=1234567)`  _function_

> Stationary (Politis-Romano) bootstrap CI for serially-correlated data.
>
> Resamples geometric-length blocks (expected length ``mean_block``) wrapping
> around the series, preserving short-range dependence, then takes the
> percentile interval. Wider than the IID :func:`bootstrap_ci` for positively
> autocorrelated series (it does not spuriously shrink the variance).
> ``statistic`` defaults to the sample mean.

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

## reverse_ad

### `Var(value, _parents=(), _backward=None)`  _class_

> A scalar node on the autodiff tape.
>
> Wrap each independent input in a ``Var``; build an expression with the arithmetic operators
> and the elementary methods (``exp``, ``log``, ``sin``, ``cos``, ``tanh``, ``sqrt``, ...).
> Calling :meth:`backward` on the output populates ``.grad`` on every node reachable from it.

### `reverse_gradient(f, xs)`  _function_

> Gradient of scalar ``f`` at the point ``xs`` (a sequence of floats).
>
> ``f`` takes a list of :class:`Var` and returns a single :class:`Var`. Returns the list of
> partial derivatives ``[df/dx_0, ...]`` via one reverse pass.

## reverse_jacobian

### `reverse_gradient_vector(f, x)`  _function_

> Exact gradient of scalar ``f`` at ``x`` (thin wrapper returning a plain list).
>
> ``f`` takes a list of :class:`Var` and returns a single :class:`Var`.

### `reverse_hessian(f, x, h=1e-05)`  _function_

> Hessian of scalar ``f: R^n -> R`` at ``x`` by differencing the exact reverse gradient.
>
> ``f`` takes a list of :class:`Var` and returns a single :class:`Var`. Each column is a
> central difference of the (exact) gradient, ``(grad f(x + h e_j) - grad f(x - h e_j)) /
> (2 h)``; the result is symmetrized. Only the outer difference carries truncation error, so
> accuracy is far better than second-differencing ``f`` itself.

### `reverse_jacobian(f, x)`  _function_

> Exact Jacobian of ``f: R^n -> R^m`` at ``x`` via reverse-mode autodiff.
>
> ``f`` takes a list of :class:`Var` and returns a list of :class:`Var` (length ``m``).
> Returns the ``m x n`` Jacobian as a list of rows; row ``i`` is ``grad f_i``. Each output
> component is differentiated by its own backward pass.

## richardson

### `richardson_extrapolate(estimates, p=1.0, t=2.0)`  _function_

> Richardson-extrapolate a sequence of step-halved estimates to the ``h -> 0`` limit.
>
> ``estimates[i] = A(h / t^i)``, ordered from coarsest to finest. ``p`` is the leading
> error exponent (1 for first-order, 2 for a central difference, etc.), ``t`` the step
> ratio between consecutive estimates. Returns the best (last-diagonal) extrapolated
> value. Successive columns cancel the ``h^p, h^{p+1}, ...`` error terms.

### `richardson_table(estimates, p=1.0, t=2.0)`  _function_

> Full Richardson tableau (list of rows) for inspecting convergence.
>
> Row ``i`` holds ``T[i][0..i]``; the diagonal ``T[i][i]`` is the order-``i``
> extrapolation. Handy to watch the estimate stabilize down the diagonal.

## richardson_derivative

### `ridders_derivative(f, x, h=None, con=1.4, safe=2.0, ntab=10)`  _function_

> First derivative ``f'(x)`` by Ridders' polynomial extrapolation.
>
> ``h`` is the initial step (defaults to a scale-aware value); ``con`` is the step
> shrink factor per row, ``ntab`` the tableau size. Returns ``(derivative, error)``
> where ``error`` is the estimated absolute error. Raises if ``f`` cannot be improved
> at all.

### `ridders_second_derivative(f, x, h=None, con=1.4, safe=2.0, ntab=10)`  _function_

> Second derivative ``f''(x)`` by Ridders extrapolation of the central formula.
>
> Uses the three-point second difference ``(f(x+h) - 2 f(x) + f(x-h)) / h^2`` at a
> shrinking step sequence with Richardson extrapolation. Returns
> ``(second_derivative, error)``.

## ridge

### `ridge_regression(X, y, alpha=1.0, add_intercept=True)`  _function_

> Fit an L2-penalized (ridge) regression.
>
> Parameters
> ----------
> X : list[list[float]]
>     Design matrix, ``n`` rows of regressors.
> y : list[float]
>     Response vector.
> alpha : float
>     Ridge penalty ``lambda`` (>= 0). 0 reproduces OLS; larger shrinks the
>     slope coefficients toward zero.
> add_intercept : bool
>     Prepend an (unpenalized) intercept column.
>
> Returns
> -------
> dict
>     ``coefficients`` (intercept first if added), ``fitted``, ``residuals``,
>     ``r_squared``.

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

## riskmeasures

### `component_expected_shortfall(component_pnls, confidence=0.95)`  _function_

> Euler (component) expected-shortfall allocation across sub-portfolios.
>
> ``component_pnls`` is a list of aligned P&L series whose sum is the portfolio.
> Each component's ES contribution is the average of its own losses over the
> scenarios where the *total* portfolio is in its worst ``1 - confidence`` tail;
> the contributions sum to the portfolio :func:`expected_shortfall` (Euler's
> theorem for the positively-homogeneous ES). Returns the list of contributions.

### `entropic_risk(pnl, risk_aversion=1.0)`  _function_

> Entropic (exponential) risk measure ``(1/theta) ln E[e^{-theta X}]``.
>
> The exponential-utility certainty equivalent of the loss; ``theta =
> risk_aversion``. Convex and increasing in ``theta``; approaches the mean loss
> ``-E[X]`` as ``theta -> 0`` and the worst loss as ``theta -> inf``.

### `expectile(pnl, tau=0.95)`  _function_

> Expectile of a P&L sample at level ``tau`` (returned as a positive loss).
>
> The ``tau``-expectile ``e`` solves the asymmetric-least-squares first-order
> condition ``tau * E[(X - e)_+] = (1 - tau) * E[(e - X)_-]`` on the loss variable
> ``X = -pnl``. The expectile is the only risk measure that is both coherent (for
> ``tau >= 0.5``) and elicitable, unlike VaR (elicitable, not coherent) and ES
> (coherent, not elicitable). ``tau = 0.5`` gives the mean loss; larger ``tau``
> weights the right (loss) tail more. Solved by bisection on the monotone
> condition.

### `is_subadditive(pnl_a, pnl_b, confidence=0.95)`  _function_

> Check the subadditivity ``rho(A+B) <= rho(A) + rho(B)`` for expected shortfall.
>
> Adds the two P&L series scenario-by-scenario and compares the combined
> expected shortfall to the sum of the standalone ones. Expected shortfall is
> coherent, so this holds (up to a tiny numerical tolerance) for any two aligned
> series -- a diagnostic that diversification never increases ES.

### `sample_expected_shortfall(pnl, confidence=0.95)`  _function_

> Expected shortfall (CVaR): mean loss in the worst ``1 - confidence`` tail.
>
> Averages the losses at or beyond the VaR quantile. A coherent risk measure,
> always at least the :func:`value_at_risk`.

### `spectral_risk_exponential(pnl, risk_aversion=5.0)`  _function_

> Spectral risk measure with an exponential risk-aversion spectrum.
>
> Weights the sorted losses by the normalized decreasing spectrum
> ``phi(p) ~ e^{-k(1-p)}`` (heavier weight on worse losses). Coherent for any
> decreasing non-negative spectrum; larger ``risk_aversion`` concentrates weight
> on the tail, raising the measure toward the worst loss.

### `value_at_risk(pnl, confidence=0.95)`  _function_

> Historical value-at-risk at ``confidence`` (a positive loss magnitude).
>
> The ``confidence`` quantile of the loss distribution (``-pnl``). Uses the
> lower-index empirical quantile so the VaR is a realized sample loss.

## rls

### `RecursiveLeastSquares(n_features, forgetting=1.0, delta=1000000.0)`  _class_

> Online multivariate least squares with an optional forgetting factor.
>
> Feed observations one at a time with :meth:`update`; read the current fit from
> :attr:`beta`. ``n_features`` is the regressor count (include a constant 1 in each
> ``x`` for an intercept). ``forgetting`` in ``(0, 1]`` down-weights past data
> (``1`` = ordinary growing-window OLS); ``delta`` sets the prior ``P = delta * I``
> (large = diffuse prior). With ``forgetting = 1`` the estimate matches batch OLS
> once enough points have arrived.

### `recursive_least_squares(X, y, forgetting=1.0, delta=1000000.0)`  _function_

> Fit RLS over a whole dataset and return the final coefficient vector.
>
> Convenience wrapper: streams the rows of ``X`` (each already including any
> intercept column) through :class:`RecursiveLeastSquares`. With ``forgetting = 1``
> the result matches batch OLS on the same design.

## rmt

### `clip_correlation_eigenvalues(correlation, n_obs)`  _function_

> Denoise a correlation matrix by clipping sub-Marchenko-Pastur eigenvalues.
>
> Eigen-decomposes ``correlation``, replaces every eigenvalue below the
> Marchenko-Pastur edge with the average of those noise eigenvalues (keeping the
> signal eigenvalues), and rebuilds the matrix, then rescales the diagonal back to
> exactly one. The trace is preserved and the result is a valid, better-conditioned
> correlation matrix. ``n_obs`` is the number of observations used to estimate
> ``correlation``.

### `marchenko_pastur_edge(n_assets, n_obs)`  _function_

> Upper edge ``(1 + sqrt(N/T))^2`` of the Marchenko-Pastur spectrum.
>
> Eigenvalues of a noise correlation matrix (unit variances) lie below this; those
> above it are candidate signal. Requires ``n_obs >= n_assets`` for a full-rank
> sample.

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

## robust_scale

### `biweight_midvariance(x, c=9.0)`  _function_

> Biweight midvariance: a robust variance that smoothly downweights outliers.
>
> Returns the square root (a robust standard deviation). Points more than ``c`` MADs
> from the median get zero weight; the tuning constant ``c=9`` gives ~87% efficiency at
> the normal. Reduces to a near-standard-deviation on clean Gaussian data.

### `qn_scale(x)`  _function_

> Rousseeuw-Croux Qn scale estimator (50% breakdown, location-free).
>
> The (roughly) first quartile of the pairwise distances ``|x_i - x_j|`` (``i < j``),
> scaled by ``2.2219`` for asymptotic consistency with the normal standard deviation.
> Needs at least two points; robust to up to half the data being outliers. No
> finite-sample correction is applied, so small-sample values differ slightly from
> implementations that include one.

### `sn_scale(x)`  _function_

> Rousseeuw-Croux Sn scale estimator (50% breakdown, no location needed).
>
> ``1.1926 * median_i( median_j |x_i - x_j| )`` -- the outer/inner medians give a robust
> spread that, unlike the MAD, does not assume a symmetric distribution. Needs at least
> two points.

## robust_stats

### `interquartile_range(x, scale=False)`  _function_

> Interquartile range ``Q3 - Q1`` of ``x``.
>
> With ``scale=True`` it is divided by 1.349 to give a normal-consistent scale
> estimate (the IQR of a standard normal). Robust to outliers in the outer
> quartiles.

### `median_absolute_deviation(x, scale=True)`  _function_

> Median absolute deviation of ``x``.
>
> ``MAD = median(|x_i - median(x)|)``. With ``scale=True`` (default) it is
> multiplied by 1.4826 so it consistently estimates the standard deviation for
> Gaussian data. Breakdown point 50%: up to half the data can be corrupted
> before it blows up. Zero for constant data.

### `trimmed_mean(x, proportion=0.1)`  _function_

> Mean of ``x`` after discarding a ``proportion`` fraction from each tail.
>
> ``proportion`` must be in ``[0, 0.5)``. With ``proportion = 0`` this is the
> ordinary mean; larger values give a more robust central estimate.

### `winsorize(x, limit=0.05)`  _function_

> Clip the tails of ``x`` to the ``limit`` / ``1 - limit`` quantiles.
>
> Returns a new list with values below the lower quantile raised to it and
> values above the upper quantile lowered to it -- bounding the influence of
> extremes without discarding observations. ``limit`` must be in ``[0, 0.5)``.

## root_scan

### `count_sign_changes(f, a, b, n=1000)`  _function_

> Number of sign changes of ``f`` on ``[a, b]`` over an ``n``-point grid.
>
> A quick lower bound on the number of simple roots (each sign change brackets at least
> one). Cheaper than :func:`find_all_roots` when only the count is needed.

### `find_all_roots(f, a, b, n=1000, tol=1e-12)`  _function_

> All sign-changing roots of ``f`` on ``[a, b]`` via a grid scan plus Brent refinement.
>
> Splits ``[a, b]`` into ``n`` subintervals, and wherever ``f`` changes sign (or hits
> exactly zero at a node) brackets and refines a root with Brent's method. Returns the
> roots in increasing order, de-duplicated. Increase ``n`` to catch roots closer than
> the grid spacing.

## rootfind

### `bisection(f, lo, hi, tol=1e-12, max_iter=200)`  _function_

> Bisection root of ``f`` on ``[lo, hi]`` (requires a sign change).
>
> Halves the bracket until it is narrower than ``tol``. Guaranteed to converge
> for a continuous ``f`` with ``f(lo) f(hi) < 0``.

### `brent(f, lo, hi, tol=1e-12, max_iter=200)`  _function_

> Brent's method root of ``f`` on ``[lo, hi]`` (requires a sign change).
>
> Combines bisection with secant and inverse-quadratic interpolation for
> superlinear convergence while retaining bisection's guaranteed bracketing.

### `newton(f, fprime, x0, tol=1e-12, max_iter=100, lo=None, hi=None)`  _function_

> Newton's method with an optional bisection safeguard.
>
> Steps ``x -= f(x)/f'(x)``; if ``lo``/``hi`` bounds are given, a step leaving
> the bracket (or a zero derivative) falls back to a bisection step. Converges
> quadratically near a simple root when the derivative is well-behaved.

## rotation

### `euler_to_matrix(yaw, pitch, roll)`  _function_

> Z-Y-X Euler angles directly to a 3x3 rotation matrix.

### `euler_to_quat(yaw, pitch, roll)`  _function_

> Z-Y-X intrinsic Euler angles (yaw, pitch, roll, radians) to a quaternion ``(w,x,y,z)``.

### `matrix_to_euler(m)`  _function_

> 3x3 rotation matrix to Z-Y-X Euler angles ``(yaw, pitch, roll)``.

### `matrix_to_quat(m)`  _function_

> Convert a 3x3 rotation matrix to a unit quaternion ``(w, x, y, z)`` (Shepperd's method).

### `quat_to_euler(q)`  _function_

> Quaternion ``(w, x, y, z)`` to Z-Y-X Euler angles ``(yaw, pitch, roll)`` in radians.

### `quat_to_matrix(q)`  _function_

> Convert a unit quaternion ``(w, x, y, z)`` to a 3x3 rotation matrix (row lists).

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

## ruin

### `gamblers_ruin_probability(start, target, win_prob)`  _function_

> Probability of hitting 0 before ``target`` in a unit-stake random walk.
>
> Starting from integer wealth ``start`` with per-round win probability
> ``win_prob`` (win +1, lose -1) and an absorbing target ``target``. The classic
> gambler's-ruin formula; for a fair game (``win_prob = 0.5``) it is the linear
> ``1 - start/target``. Ruin is certain against an unfavorable game as the target
> grows.

### `risk_of_ruin_units(win_prob, target_units)`  _function_

> Risk of ruin starting with ``target_units`` units against a one-unit target run.
>
> A trader with a bankroll of ``target_units`` betting one unit at a time with
> edge ``win_prob``: the probability of losing the whole bankroll before doubling
> it is :func:`gamblers_ruin_probability` from ``target_units`` toward
> ``2 * target_units``. A convenient bankroll-in-units risk measure; falls fast in
> the number of units when the game is favorable.

### `ruin_probability_gbm(loss, mu, sigma)`  _function_

> Probability a drifting log-equity ever falls by at least the fraction ``loss``.
>
> For log-equity following a Brownian motion with drift ``mu > 0`` and volatility
> ``sigma``, the chance the account ever drops to ``(1 - loss)`` of its starting
> value is the first-passage law ``(1 - loss)^{2 mu / sigma^2}``. ``loss`` is a
> fraction in ``(0, 1)``. Decreasing in the drift-to-variance ratio: a stronger
> edge makes a given loss less likely. Returns 1 for a non-positive drift (a
> driftless or losing account eventually hits any loss level almost surely).

## running_median

### `RunningMedian(values=None)`  _class_

> Maintains the exact median of all values inserted so far.
>
> ``push(x)`` adds a value in ``O(log n)``; ``median()`` returns the current median in
> ``O(1)``. ``lower`` is a max-heap (stored as negated values) of the smaller half; ``upper``
> a min-heap of the larger half.

## running_moments

### `RunningMoments(values=None)`  _class_

> One-pass accumulator for mean, variance, skewness, and excess kurtosis.
>
> Feed values with :meth:`update` (or a whole iterable to the constructor). Read
> :attr:`mean`, :meth:`variance`, :meth:`skewness`, :meth:`kurtosis` at any time.
> ``+`` merges two accumulators into one covering both samples, exactly.

## runs_test

### `runs_test(values)`  _function_

> Runs test on a numeric series, dichotomized about its median.
>
> Values above the median are one symbol, below the other; values exactly equal to
> the median are dropped. Returns ``(z, p_value)`` for the null that the sequence
> of above/below signs is random. A small p-value rejects randomness: ``z < 0`` for
> trending/clustered data, ``z > 0`` for over-alternating (mean-reverting) data.

### `runs_test_binary(sequence)`  _function_

> Wald-Wolfowitz runs test on a two-symbol sequence.
>
> ``sequence`` is any list of two distinct values (e.g. 0/1, +/-). Returns
> ``(z, p_value)``; a negative ``z`` (few runs) signals clustering, a positive
> ``z`` (many runs) over-alternation. Requires at least one of each symbol.

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

## sample_rate

### `downsample(x, factor, numtaps=65)`  _function_

> Downsample ``x`` by an integer ``factor`` with an anti-alias low-pass filter.
>
> Low-pass filters at the new Nyquist (cutoff ``1/factor``) to prevent aliasing, then
> keeps every ``factor``-th sample. Returns a signal ``ceil(len(x)/factor)`` long.
> ``numtaps`` sizes the anti-alias filter (odd).

### `resample_rational(x, up, down, numtaps=65)`  _function_

> Resample ``x`` by the rational ratio ``up/down`` (upsample then downsample).
>
> Interpolates by ``up`` and decimates by ``down`` through a single anti-alias/anti-
> imaging low-pass at cutoff ``1/max(up, down)``, so both imaging and aliasing are
> suppressed. Returns a signal about ``len(x) * up / down`` long.

### `sinc_interp(x, positions)`  _function_

> Whittaker-Shannon reconstruction of a band-limited signal at arbitrary ``positions``.
>
> Treats ``x`` as samples at integer indices ``0..len(x)-1`` of a signal band-limited
> to the Nyquist frequency and evaluates ``sum_n x[n] sinc(t - n)`` at each ``t`` in
> ``positions``. Exact at integer positions; the ideal interpolation between them.

### `upsample(x, factor, numtaps=65)`  _function_

> Upsample ``x`` by an integer ``factor`` with an anti-imaging low-pass filter.
>
> Inserts ``factor - 1`` zeros between samples and low-pass filters at the original
> Nyquist (cutoff ``1/factor``), scaling by ``factor`` to preserve amplitude. Returns a
> signal ``factor`` times as long. ``numtaps`` sizes the interpolation filter (odd).

## samplers

### `sample_exponential(n, rate=1.0, seed=1234567)`  _function_

> ``n`` exponential draws with the given ``rate`` (mean ``1/rate``), inverse-CDF method.

### `sample_gamma(n, shape, scale=1.0, seed=1234567)`  _function_

> ``n`` gamma draws with the given ``shape`` (``k``) and ``scale`` (``theta``).
>
> Mean ``shape*scale``, variance ``shape*scale^2``. Uses Marsaglia-Tsang (with the
> small-shape boost). ``shape`` and ``scale`` must be positive.

### `sample_normal(n, mu=0.0, sigma=1.0, seed=1234567)`  _function_

> ``n`` normal draws with mean ``mu`` and standard deviation ``sigma`` (Box-Muller).

### `sample_poisson(n, lam, seed=1234567)`  _function_

> ``n`` Poisson draws with mean ``lam`` (Knuth's algorithm).

## savgol

### `savgol_coeffs(window, degree, deriv=0)`  _function_

> Savitzky-Golay convolution coefficients for a window of odd length ``window``.
>
> Returns the ``window`` weights that, dotted with the windowed samples, give the
> fitted value (``deriv = 0``) or the ``deriv``-th derivative at the window center.
> ``window`` must be odd and larger than ``degree``. For unit spacing; scale a
> derivative by ``1 / h^deriv`` for spacing ``h``.

### `savgol_filter(data, window, degree, deriv=0)`  _function_

> Apply a Savitzky-Golay filter to ``data``.
>
> Smooths (``deriv = 0``) or differentiates the series with a length-``window``,
> degree-``degree`` polynomial fit. Interior points use the centered convolution;
> the ``half`` points at each end are fitted from the nearest full window (so the
> output has the same length as the input). Exactly reproduces polynomials up to
> ``degree``.

## scale_tests

### `ansari_bradley_test(x, y)`  _function_

> Ansari-Bradley two-sample test of equal dispersion.
>
> Scores the pooled ranks from the outside in (``min(r, N+1-r)``) and sums the
> scores over ``x``. Returns a dict with the ``statistic`` (that sum), the ``z``
> normal approximation and the two-sided ``p_value``. Assumes the two samples share
> a location; a smaller ``x`` spread pushes ``x`` toward the center (higher scores).

### `mood_test(x, y)`  _function_

> Mood two-sample test of equal dispersion.
>
> Sums the squared rank deviations ``(r - (N+1)/2)^2`` over ``x``. Returns a dict
> with the ``statistic``, the ``z`` normal approximation and the two-sided
> ``p_value``. A more dispersed ``x`` sends its values to the rank extremes and
> raises the statistic. Assumes a common location.

## scaling

### `fit_min_max(X)`  _function_

> Fit a min-max scaler: per-column ``(min, range)`` (zero range -> 1).

### `fit_robust(X)`  _function_

> Fit a robust scaler: per-column ``(median, IQR)`` (zero IQR -> 1).

### `fit_standardize(X)`  _function_

> Fit a z-score scaler: per-column ``(mean, std)`` (zero std -> 1).

### `scale_inverse_transform(params, X)`  _function_

> Undo a fitted scaler: ``x * scale + center`` column-wise.

### `scale_transform(params, X)`  _function_

> Apply a fitted scaler: ``(x - center) / scale`` column-wise.

## scc

### `condensation(graph)`  _function_

> Collapse each SCC to a node, returning ``(component_of, dag)``.
>
> ``component_of`` maps each original node to its component id (``0..k-1``); ``dag`` is the
> condensed graph ``{comp_id: [successor_comp_ids]}`` with no duplicate edges and no
> self-loops. The condensation is always acyclic.

### `is_strongly_connected(graph)`  _function_

> True if the whole graph is a single strongly connected component.

### `number_of_sccs(graph)`  _function_

> Count the strongly connected components of ``graph``.

### `strongly_connected_components(graph)`  _function_

> Return the SCCs of a directed ``graph`` as a list of node lists.
>
> ``graph`` is ``{node: [successors]}``; every node must appear as a key. Components are
> returned in reverse topological order of the condensation (a component appears before
> the components it can reach), which is the natural output order of Tarjan's algorithm.

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

## scheduling

### `activity_selection(intervals)`  _function_

> Maximum number of mutually non-overlapping intervals (earliest-finish greedy).
>
> Returns the chosen ``(start, end)`` list in order. Ignores weights; half-open overlap.

### `min_rooms(intervals)`  _function_

> Minimum resources to run every interval, i.e. the peak simultaneous count.
>
> An interval ending exactly when another begins does not need a separate room
> (half-open). Returns an integer.

### `weighted_interval_schedule(intervals)`  _function_

> Maximum-weight non-overlapping subset. Returns ``(total_weight, chosen)``.
>
> ``intervals`` is a list of ``(start, end, weight)`` with ``start < end`` and
> ``weight >= 0``. ``chosen`` is the selected intervals in start order. Half-open
> intervals: one ending exactly when another starts do not conflict.

## sde

### `euler_maruyama(drift, diffusion, x0, t, n_steps, seed=1234567)`  _function_

> Euler-Maruyama path of ``dX = drift(X, t) dt + diffusion(X, t) dW``.
>
> ``drift`` and ``diffusion`` are callables ``(x, t) -> float``. Returns the path as a
> list of ``n_steps + 1`` states at times ``0, dt, ..., t``. Strong order 0.5.

### `gbm_paths(mu, sigma, x0, t, n_steps, n_paths, seed=1234567, scheme='milstein')`  _function_

> Simulate geometric Brownian motion paths (``dS = mu S dt + sigma S dW``).
>
> Convenience wrapper over :func:`milstein` (default) or :func:`euler_maruyama` with the
> GBM drift and diffusion. Returns a list of ``n_paths`` state paths. The sample mean of
> the terminal value approaches the analytic ``x0 exp(mu t)`` as ``n_paths`` grows.

### `milstein(drift, diffusion, diffusion_prime, x0, t, n_steps, seed=1234567)`  _function_

> Milstein path with the ``0.5 b b' (dW^2 - dt)`` correction (strong order 1.0).
>
> Adds the derivative of the diffusion ``diffusion_prime(x, t) = d b/d x`` to the
> Euler-Maruyama step, giving strong order 1.0 for state-dependent diffusions. Same
> signature otherwise; returns the ``n_steps + 1`` state path.

## selection

### `kth_smallest(values, k)`  _function_

> The ``k``-th smallest element (0-indexed) via linear-time quickselect.
>
> ``k`` in ``[0, len(values) - 1]``; ``k=0`` is the minimum. Uses the median-of-medians
> pivot for a guaranteed ``O(n)`` worst case. Does not modify the input.

### `median(values)`  _function_

> Median of ``values`` (average of the two middle elements for even length).

### `top_k(values, k, largest=True)`  _function_

> The ``k`` largest (or smallest) elements of ``values``, sorted.
>
> ``largest=True`` returns the top ``k`` in descending order; ``largest=False`` the
> bottom ``k`` ascending. ``k`` is clamped to the list length. Uses quickselect to find
> the threshold, then sorts only the ``k`` selected elements.

## sequence_accel

### `aitken(sequence)`  _function_

> Aitken's delta-squared acceleration of a sequence.
>
> Returns the accelerated sequence, two elements shorter than the input. Each output
> term extrapolates the limit from three consecutive input terms; for a linearly
> convergent sequence it converges markedly faster.

### `shanks(sequence)`  _function_

> Shanks transform of a sequence (one pass; same formula as :func:`aitken`).

### `steffensen(g, x0, tol=1e-12, max_iter=100)`  _function_

> Steffensen's method: quadratically-convergent fixed point of ``x = g(x)``.
>
> Applies Aitken acceleration to the fixed-point iterates, so it converges
> quadratically like Newton's method but needs no derivative -- only the map ``g``.
> Returns a dict with the ``root`` (the fixed point), ``iterations`` and
> ``converged``. Raises if a zero denominator stalls the iteration before
> convergence.

## sequences

### `longest_increasing_subsequence(x, strict=True)`  _function_

> Longest (strictly by default) increasing subsequence of ``x``.
>
> Returns an actual longest subsequence as a list (not necessarily contiguous). With
> ``strict=False`` allows equal consecutive values (non-decreasing). ``O(n log n)`` via
> patience sorting with predecessor tracking. Empty input yields ``[]``.

### `longest_run(x)`  _function_

> Longest run of a single repeated value: ``(value, length, start_index)``.
>
> Scans for the longest maximal streak of equal adjacent elements. Ties break to the
> earliest run. Raises on empty input.

### `maximum_subarray(x)`  _function_

> Maximum-sum contiguous subarray (Kadane): ``(sum, start, end)`` inclusive.
>
> Returns the largest achievable sum of a non-empty contiguous slice and its bounds.
> Handles all-negative inputs (returns the single largest element). Raises on empty
> input.

## serial_correlation

### `box_pierce(x, lags=10)`  _function_

> Box-Pierce portmanteau statistic and p-value.
>
> Returns ``(Q, p_value)`` with ``Q = n * sum_{k=1}^{lags} rho_k^2`` and the
> chi-square(``lags``) p-value. Small p rejects the white-noise null.

### `durbin_watson(x)`  _function_

> Durbin-Watson statistic ``d = sum (x_t - x_{t-1})^2 / sum x_t^2``.
>
> Approximately ``2(1 - rho_1)``: near 2 = no first-order autocorrelation, near
> 0 = strong positive, near 4 = strong negative. Computed on the demeaned series.

### `ljung_box(x, lags=10)`  _function_

> Ljung-Box portmanteau statistic and p-value (small-sample refinement).
>
> Returns ``(Q, p_value)`` with ``Q = n(n+2) sum_{k=1}^{lags} rho_k^2 / (n-k)``
> and the chi-square(``lags``) p-value. Large Q / small p rejects "the first
> ``lags`` autocorrelations are jointly zero".

## series_transform

### `euler_transform(terms)`  _function_

> Euler transform of an alternating series ``sum_{k>=0} (-1)^k terms[k]``.
>
> ``terms`` are the non-negative magnitudes ``a_k`` (the alternating sign is applied
> internally). Returns the accelerated estimate of the sum via
> ``sum_k (-1)^k a_k = sum_n (-1)^0 * Delta^n a_0 / 2^{n+1}`` -- forward differences
> reweighted by powers of one half, which converges geometrically even when the raw
> alternating series crawls.

### `wynn_epsilon(partial_sums)`  _function_

> Accelerate a sequence of partial sums with Wynn's epsilon algorithm.
>
> ``partial_sums`` is the list ``[s_0, s_1, ...]``. Returns the best (last stable)
> even-column estimate of the limit. The epsilon table is built with the recurrence
> ``eps[k+1][j] = eps[k-1][j+1] + 1 / (eps[k][j+1] - eps[k][j])``; the even columns
> hold the accelerated limits.

## shout

### `ladder_call(S, K, rungs, t, r, sigma, steps=200, q=0.0)`  _function_

> Ladder call price on a CRR tree.
>
> A ladder call locks in a guaranteed payoff each time the underlying touches a
> preset rung ``L_i > K``: the terminal payoff is
> ``max(S_T - K, max_touched L_i - K, 0)``. Priced by carrying the highest rung
> reached along each tree path (a state variable on the sorted rungs) through
> backward induction. At least the vanilla call; more/higher rungs raise the
> value up to the shout-like limit.

### `shout_call(S, K, t, r, sigma, steps=200, q=0.0)`  _function_

> Shout call price on a CRR tree.
>
> Once shouted at spot ``S*``, the remaining claim pays
> ``max(S_T - K, S* - K)`` -- a guaranteed ``S* - K`` plus a call struck at ``S*``.
> Its value at the shout node is ``(S* - K) e^{-r tau} + call(S*, S*, tau)`` for
> remaining time ``tau``, taken only when positive. Backward induction compares
> shouting versus continuing. At least the vanilla call value.

## shrinkage

### `constant_correlation_target(cov)`  _function_

> Constant-correlation shrinkage target from a covariance matrix.
>
> Keeps each asset's own variance but replaces every pairwise correlation with
> the average sample correlation ``r_bar``: ``F_ij = r_bar * sqrt(S_ii S_jj)``
> for ``i != j`` and ``F_ii = S_ii``.

### `ledoit_wolf_shrinkage(returns)`  _function_

> Ledoit-Wolf (2004) constant-correlation shrinkage covariance estimate.
>
> Parameters
> ----------
> returns : sequence of sequence of float
>     ``n`` observations of ``p`` asset returns.
>
> Returns
> -------
> (sigma_hat, delta) : (list[list[float]], float)
>     The shrunk ``p x p`` covariance matrix and the shrinkage intensity
>     ``delta`` in ``[0, 1]``. ``delta`` rises toward 1 as the sample estimate
>     gets noisier (small ``n``) and falls toward 0 as it gets reliable.

### `sample_covariance(returns)`  _function_

> Maximum-likelihood sample covariance (divisor ``n``) of a return matrix.
>
> ``returns`` is a sequence of ``n`` observations, each a length-``p`` sequence.
> Returns a ``p x p`` list-of-lists.

## siegel_regression

### `repeated_median_regression(x, y)`  _function_

> Siegel repeated-median regression; returns ``(slope, intercept)``.
>
> For each point, the median of its slopes to all other points; the overall slope is
> the median of those. 50% breakdown point -- robust to nearly half the data being
> corrupted. Points sharing an ``x`` value contribute no slope for that pair.

## sieve

### `nth_prime(n)`  _function_

> The ``n``-th prime (1-indexed): ``nth_prime(1) == 2``.
>
> Grows the sieve bound with the prime-number-theorem estimate ``n(ln n + ln ln n)``
> until enough primes are found. ``n >= 1``.

### `prime_count(limit)`  _function_

> Number of primes ``<= limit`` (the prime-counting function ``pi(limit)``).

### `primes_up_to(limit)`  _function_

> All primes ``<= limit`` by the Sieve of Eratosthenes (a sorted list).
>
> ``limit < 2`` yields an empty list. Marks multiples starting from each prime's square.

### `smallest_prime_factors(limit)`  _function_

> Smallest-prime-factor table for ``0..limit`` (a list; ``spf[k]`` = least prime dividing ``k``).
>
> ``spf[0] = spf[1] = 0``. With it, any ``k <= limit`` factorizes in ``O(log k)`` by
> repeatedly dividing by ``spf[k]``. Built by a linear-ish sieve.

## signal_features

### `crest_factor(x)`  _function_

> Crest factor: peak amplitude divided by RMS.
>
> High for impulsive/peaky signals (a lone spike), low for signals that fill their
> range (~1.41 for a sine, 1.0 for a square wave). Zero-signal raises.

### `rms(x)`  _function_

> Root-mean-square amplitude of a signal.

### `spectral_bandwidth(x)`  _function_

> Spectral bandwidth: the power-weighted standard deviation about the centroid.
>
> Measures how spread out the spectrum is (narrow for a pure tone, wide for noise).

### `spectral_centroid(x)`  _function_

> Spectral centroid: the power-weighted mean frequency (cycles/sample).
>
> The spectrum's "center of mass" -- higher for brighter/higher-pitched signals. Zero
> for a DC-only (constant) signal. Computed from the mean-removed periodogram.

### `spectral_flatness(x)`  _function_

> Spectral flatness (Wiener entropy): geometric mean / arithmetic mean of the spectrum.
>
> Near ``1`` for white-noise-like flat spectra, near ``0`` for tonal signals with power
> concentrated in a few bins. Uses the mean-removed positive-frequency power bins.

### `zero_crossing_rate(x)`  _function_

> Fraction of adjacent sample pairs that straddle zero (sign changes / (n-1)).
>
> A rough pitch/noisiness proxy: high for noisy or high-frequency signals, low for
> smooth low-frequency ones. Returns a value in ``[0, 1]``.

## signals

### `average_true_range(highs, lows, closes, window=14)`  _function_

> Average true range (Wilder): mean of the true range over a trailing window.
>
> True range at ``t`` is ``max(high-low, |high-prev_close|, |low-prev_close|)``;
> ATR smooths it with Wilder's moving average. A non-negative volatility measure
> in price units. Returns one value per position from index ``window`` on.

### `bollinger_bands(series, window=20, num_std=2.0)`  _function_

> Bollinger bands: ``(lower, middle, upper)`` lists over a trailing window.
>
> Middle is the :func:`sma`; the bands are ``middle +/- num_std * rolling std``.
> Price closing above the upper / below the lower band flags stretched moves.
> Returns three aligned lists, one value per window position.

### `donchian_channel(highs, lows, window=20)`  _function_

> Donchian channel: rolling ``(lowest_low, highest_high)`` over a window.
>
> The channel a breakout system trades: a close above the prior highest high is
> a long breakout, below the lowest low a short. Returns ``(lower, upper)`` lists
> with ``upper >= lower`` at every position.

### `ema(series, span)`  _function_

> Exponential moving average with smoothing ``alpha = 2/(span+1)``.
>
> Recursive ``e_t = alpha x_t + (1 - alpha) e_{t-1}`` seeded at the first point.
> Reacts faster than the :func:`sma` of the same length (less lag).

### `macd(series, fast=12, slow=26, signal=9)`  _function_

> MACD line, signal line, and histogram.
>
> MACD line = ``EMA(fast) - EMA(slow)``; signal line = ``EMA(signal)`` of the
> MACD line; histogram = MACD - signal. Returns ``(macd_line, signal_line,
> histogram)`` aligned to the series. Positive MACD indicates the fast average
> above the slow (up-momentum).

### `rolling_zscore(series, window)`  _function_

> Rolling z-score ``(x_t - mean) / std`` over each trailing ``window``.
>
> Standardizes the latest point against its window; a mean-reversion / breakout
> signal. Windows with zero variance yield 0. One value per position from index
> ``window - 1`` on.

### `rsi(series, window=14)`  _function_

> Relative strength index over a trailing ``window`` (Wilder's smoothing).
>
> ``RSI = 100 - 100/(1 + avg_gain/avg_loss)`` in ``[0, 100]``. Above 70 is
> conventionally overbought, below 30 oversold; near 100 in a strong uptrend.
> Returns one value per position from index ``window`` on.

### `sma(series, window)`  _function_

> Simple moving average over each trailing ``window`` (list, one per position).
>
> Returns ``len(series) - window + 1`` values, each the mean of that window.

### `time_series_momentum(prices, lookback)`  _function_

> Sign of the trailing ``lookback``-period return: +1 up, -1 down, 0 flat.
>
> The time-series-momentum signal (Moskowitz-Ooi-Pedersen): go long after a
> positive past return, short after a negative one. Returns one signal per
> position from index ``lookback`` on.

## silhouette

### `silhouette_samples(X, labels)`  _function_

> Per-point silhouette values ``s(i)`` in ``[-1, 1]``.
>
> A point alone in its cluster gets ``s = 0`` by convention. Requires at least
> two distinct cluster labels.

### `silhouette_score(X, labels)`  _function_

> Mean silhouette over all points -- an overall clustering-quality score.
>
> Near 1 = dense, well-separated clusters; near 0 = overlapping; negative =
> mostly misassigned. Use it to compare label sets or pick ``k``.

## simplify

### `douglas_peucker(points, epsilon)`  _function_

> Simplify a polyline to a subset of its points within perpendicular tolerance ``epsilon``.
>
> Returns a new list containing the retained points (always including the first and
> last). Larger ``epsilon`` keeps fewer points; ``epsilon = 0`` keeps every point that
> is not exactly collinear. The result is a subsequence of the input in order.

## simulated_annealing

### `simulated_annealing(func, x0, bounds=None, T0=1.0, cooling=0.995, step=1.0, max_iter=10000, seed=1234567)`  _function_

> Minimize ``func`` from start ``x0`` by simulated annealing.
>
> ``func`` takes a length-``d`` list and returns a scalar. ``bounds`` is an optional
> list of ``(lo, hi)`` per dimension (proposals are clamped into the box). ``T0`` is
> the initial temperature, ``cooling`` the per-iteration geometric decay, ``step`` the
> proposal scale (shrinks with temperature). Returns a dict with ``x`` (best point),
> ``fun`` (its value), ``n_iter`` and ``final_temp``. Deterministic for a fixed seed.

## sinkhorn

### `cost_matrix(xs, ys, p=2)`  _function_

> Ground-cost matrix ``C[i][j] = |xs[i] - ys[j]|^p`` for scalar or vector support points.
>
> ``xs`` and ``ys`` are lists of points (each a float or an equal-length coordinate list);
> ``p`` is the exponent of the Euclidean distance (``p = 2`` gives squared distance, the
> 2-Wasserstein ground cost). Returns the ``len(xs) x len(ys)`` matrix.

### `sinkhorn(a, b, C, eps=0.05, max_iter=2000, tol=1e-09)`  _function_

> Entropic-regularized optimal transport plan between marginals ``a`` and ``b``.
>
> ``a`` (length ``n``) and ``b`` (length ``m``) are non-negative weight vectors with equal
> total mass; ``C`` is the ``n x m`` cost matrix; ``eps`` the entropic regularization strength.
> Returns a dict with ``plan`` (the ``n x m`` transport matrix), ``cost`` (``sum P*C``) and
> ``n_iter``. As ``eps -> 0`` the cost approaches the exact optimal-transport cost.
>
> Runs in the log-domain (stabilized potentials ``f, g``) so it stays accurate for small
> ``eps`` where the raw ``exp(-C/eps)`` kernel would underflow to zero.

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

### `kelly_fractions_multivariate(mean_excess_returns, cov, fraction=1.0)`  _function_

> Growth-optimal Kelly allocation across correlated assets.
>
> For a vector of excess returns with mean ``mu`` and covariance ``Sigma``, the
> continuous multivariate Kelly criterion maximizes the expected log-growth
> ``f . mu - 0.5 f . Sigma f``; the optimum is ``f* = Sigma^{-1} mu``. Returns
> the leverage vector, scaled by ``fraction`` for fractional Kelly.
>
> Reduces to the scalar ``mu / sigma^2`` for a single asset, and to the
> per-asset Kelly fractions when the covariance is diagonal (uncorrelated
> assets). Requires a positive-definite ``Sigma``.

### `kelly_growth_rate(expected_excess_return, variance, leverage)`  _function_

> Expected log-growth rate at a given leverage (continuous Kelly).
>
> ``g(f) = f*mu - 0.5 * f^2 * sigma^2``. Maximized at the full-Kelly leverage
> ``f* = mu / sigma^2``; used to compare fractional-Kelly choices.

### `kelly_growth_rate_multivariate(mean_excess_returns, cov, leverages)`  _function_

> Expected log-growth rate of a multivariate allocation.
>
> ``g(f) = f . mu - 0.5 f . Sigma f``. Maximized at
> ``f* = Sigma^{-1} mu`` (:func:`kelly_fractions_multivariate` with
> ``fraction=1``); used to compare fractional-Kelly leverage vectors.

### `neutralize(book: quantforge.portfolio.Book, greek: str, S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None, multiplier=1.0, target: float = 0.0) -> float`  _function_

> Units of the hedge option to move ``greek`` to ``target`` (default 0).
>
> ``greek`` is one of "delta", "gamma", "vega". Returns the signed quantity
> (in option units, before the multiplier is applied to notionals): solving
> ``net_greek + qty * multiplier * hedge_greek = target``.

### `vega_neutral_quantity(book: quantforge.portfolio.Book, S, K, t, r, sigma, option_type=<OptionType.CALL: 'call'>, b=None, multiplier=1.0) -> float`  _function_

> Units of a hedge option that zero the book's net vega.

## sliced_wasserstein

### `sliced_wasserstein(xs, ys, n_projections=200, p=2, seed=12345)`  _function_

> Sliced ``p``-Wasserstein distance between two multivariate point clouds.
>
> ``xs`` and ``ys`` are lists of points (each a length-``d`` coordinate list; ``d = 1`` scalars
> are also accepted as bare floats). Projects both onto ``n_projections`` random unit
> directions, takes the 1-D ``p``-Wasserstein distance along each, and returns the ``p``-mean
> ``(mean_l W_p(proj_l)^p)^{1/p}``. A seeded PCG32 stream makes the projections reproducible.

## sliding_window

### `sliding_window_max(values, k)`  _function_

> Maximum of each length-``k`` window; returns ``n - k + 1`` values.

### `sliding_window_min(values, k)`  _function_

> Minimum of each length-``k`` window; returns ``n - k + 1`` values.

### `sliding_window_sum(values, k)`  _function_

> Sum of each length-``k`` window in ``O(n)`` (running sum), ``n - k + 1`` values.

## so3

### `hat(v)`  _function_

> Skew-symmetric matrix of a 3-vector: ``hat(v) w == cross(v, w)``.

### `rodrigues(omega)`  _function_

> Exponential map: rotation vector ``omega`` -> 3x3 rotation matrix (Rodrigues' formula).

### `so3_log(R)`  _function_

> Logarithm map: rotation matrix ``R`` -> rotation vector ``omega`` (axis * angle).

### `unhat(m)`  _function_

> Inverse of :func:`hat`: the 3-vector of a skew-symmetric matrix.

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

## sparse_table

### `SparseTable(values, combine=<built-in function min>)`  _class_

> Static array supporting ``O(1)`` idempotent range queries (min/max/gcd/...).
>
> Build from a sequence and an associative *idempotent* ``combine`` (default `min`).
> ``query(lo, hi)`` returns the combined value over the inclusive index range
> ``[lo, hi]`` in constant time. The array is fixed after construction.

### `range_gcd_query(values)`  _function_

> Convenience: a `SparseTable` answering range-*gcd* queries.

### `range_max_query(values)`  _function_

> Convenience: a `SparseTable` answering range-*maximum* queries.

### `range_min_query(values)`  _function_

> Convenience: a `SparseTable` answering range-*minimum* queries.

## special

### `betainc(a, b, x)`  _function_

> Regularized incomplete beta ``I_x(a, b)`` via the Lentz continued fraction.
>
> ``I_0 = 0``, ``I_1 = 1``, and the symmetry ``I_x(a, b) = 1 - I_{1-x}(b, a)`` is
> used where the fraction converges slowly. Requires ``a, b > 0`` and
> ``0 <= x <= 1``.

### `digamma(x)`  _function_

> Digamma ``psi(x) = d/dx ln Gamma(x)`` for ``x > 0``.
>
> Recurses up to ``x >= 6`` with ``psi(x) = psi(x+1) - 1/x``, then applies the
> asymptotic (Stirling) series. Accurate to ~1e-12 for positive arguments.

### `erfinv(y)`  _function_

> Inverse error function on ``(-1, 1)``.
>
> A rational approximation (Giles) seeds a Newton-Halley refinement against
> :func:`math.erf`, giving full double precision. ``erfinv(0) = 0`` and
> ``erfinv(erf(x)) = x``.

### `gammainc(a, x)`  _function_

> Regularized lower incomplete gamma ``P(a, x) = gamma(a, x) / Gamma(a)``.
>
> Series for ``x < a + 1`` and the complement of the Lentz continued fraction
> otherwise (Numerical Recipes). ``P(a, 0) = 0`` and ``P(a, x) -> 1`` as
> ``x -> inf``. Requires ``a > 0`` and ``x >= 0``.

### `gammaincc(a, x)`  _function_

> Regularized upper incomplete gamma ``Q(a, x) = 1 - P(a, x)``.
>
> Uses the continued fraction directly for ``x >= a + 1`` (where it converges
> fast) and the series complement otherwise, so ``gammainc(a, x) + gammaincc(a, x)
> == 1`` to machine precision.

## spectral

### `dft(x)`  _function_

> Discrete Fourier transform of a real (or complex) sequence.
>
> Returns the ``n`` complex coefficients ``X_k = sum_t x_t exp(-2 pi i k t / n)``.

### `dominant_frequency(x)`  _function_

> Frequency (cycles per sample) of the largest non-DC periodogram peak.
>
> Ignores the zero-frequency (mean) component. Returns the frequency; its
> reciprocal is the dominant period in samples.

### `periodogram(x)`  _function_

> One-sided periodogram of a real series.
>
> Returns ``(freqs, power)`` where ``freqs`` are normalized frequencies in
> cycles per sample over ``[0, 0.5]`` and ``power[k] = |X_k|^2 / n``. The DC
> term is included at frequency 0. Peaks mark dominant cycles.

### `spectral_energy(x)`  _function_

> Total periodogram energy, ``sum_k |X_k|^2 / n`` over all n frequencies.
>
> By Parseval's theorem this equals ``sum_t x_t^2`` (the time-domain energy),
> which the tests use as a consistency check.

### `welch_psd(x, segment_length=None, overlap=0.5)`  _function_

> Welch's power-spectral-density estimate: averaged windowed periodograms.
>
> Splits ``x`` into overlapping segments of ``segment_length`` (default ``n // 8``,
> clamped to at least 8), applies a Hann window to each, and averages their
> periodograms. Averaging trades frequency resolution for a much lower-variance
> spectral estimate than the raw periodogram. ``overlap`` is the fractional segment
> overlap in ``[0, 1)``. Returns ``(freqs, power)`` with one-sided normalized
> frequencies in ``[0, 0.5]``.

## sphere

### `angular_distance(u, v)`  _function_

> Angle in radians between two 3-vectors (their directions), in ``[0, pi]``.
>
> Uses ``atan2(|u x v|, u . v)``, which stays accurate for both nearly-parallel and
> nearly-opposite directions (unlike ``acos`` of the dot product).

### `slerp_vectors(u, v, t)`  _function_

> Spherical-linear interpolation between unit directions ``u`` and ``v`` at ``t`` in [0,1].
>
> Returns a unit vector on the great-circle arc from ``u`` (``t=0``) to ``v`` (``t=1``),
> at constant angular speed. Falls back to normalized linear interpolation when the two
> directions are nearly identical.

### `spherical_centroid(vectors)`  _function_

> Mean direction of a set of 3-vectors: the normalized vector sum.
>
> Returns the unit vector minimizing the sum of squared chord distances (the resultant
> direction of directional statistics). Raises if the vectors sum to zero (no mean
> direction).

### `spherical_resultant_length(vectors)`  _function_

> Mean resultant length ``R`` in ``[0, 1]``: concentration of a set of directions.
>
> ``R = |sum unit(v)| / n``. ``R = 1`` means all directions coincide; ``R = 0`` means they
> are perfectly spread. The dispersion measure of spherical statistics.

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

## stern_brocot

### `best_rational_bounded(x, max_denominator)`  _function_

> Closest fraction to ``x`` with denominator ``<= max_denominator`` (Stern-Brocot descent).
>
> Walks the tree toward ``x``, stopping when the next mediant would exceed the denominator
> bound, and returns the nearer of the two current boundaries. ``x`` may be any real.

### `farey_sequence(n)`  _function_

> Return the Farey sequence of order ``n``: reduced fractions in ``[0, 1]``, ascending.
>
> Uses the ``O(1)``-per-term neighbour recurrence, so no gcd or sorting is needed.

### `mediant(a, b)`  _function_

> Mediant of two fractions ``a = p/q`` and ``b = r/s``: ``(p + r)/(q + s)``.
>
> ``a`` and ``b`` may be `Fraction` or ``(num, den)`` pairs. The mediant lies strictly
> between them and is the fraction the Stern-Brocot / Farey construction inserts.

### `stern_brocot_from_path(path)`  _function_

> Return the fraction at a given L/R ``path`` from ``1/1``.

### `stern_brocot_path(target, max_steps=10000)`  _function_

> Return the L/R path from ``1/1`` to a positive rational ``target`` in the tree.
>
> ``target`` is a `Fraction` (or something `Fraction` accepts). Returns a string of ``'L'``
> and ``'R'`` steps; the empty string means ``target == 1``. Raises for non-positive
> targets or if ``max_steps`` is exceeded (non-terminating only for irrationals).

## stft

### `istft(frames, frame_size, hop=None, window='hann', length=None)`  _function_

> Invert an :func:`stft` back to the time domain by weighted overlap-add.
>
> Applies the synthesis window to each inverse-FFT frame and overlap-adds, dividing by
> the overlap-added squared window so the reconstruction is exact wherever the window
> coverage is non-zero (independent of the constant-overlap-add condition). ``length``
> truncates the output; by default it keeps the full overlap-added span.

### `spectrogram(x, frame_size, hop=None, window='hann')`  _function_

> Spectrogram: per-frame power ``|STFT|^2`` over the lower (non-redundant) half.
>
> Returns a list of frames, each a list of ``frame_size // 2 + 1`` power values (DC up
> to Nyquist -- the rest of a real signal's spectrum is a mirror image). Each frame is
> a column in time; each entry a frequency bin.

### `stft(x, frame_size, hop=None, window='hann')`  _function_

> Short-time Fourier transform: a list of per-frame complex FFT spectra.
>
> Slides a ``frame_size`` window along ``x`` in steps of ``hop`` (default
> ``frame_size // 2``, i.e. 50% overlap), applies the ``window`` taper, and FFTs each
> frame. Returns a list of frames, each a length-``frame_size`` complex spectrum.
> ``frame_size`` must be a power of two (FFT constraint). Frames that run past the end
> are zero-padded.

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

## streaming_quantile

### `P2Quantile(p)`  _class_

> P-square streaming estimator of a single quantile ``p`` in ``[0, 1]``.
>
> Feed values with :meth:`update`; read the current estimate from :meth:`value`.
> Uses five markers and O(1) memory; accuracy improves as the stream grows. Best for
> smooth, stationary streams -- it approximates the true quantile without storing data.

### `reservoir_sample(stream, k, seed=1234567)`  _function_

> Vitter reservoir sample: ``k`` uniform items from a stream of unknown length.
>
> ``stream`` is any iterable. Returns a list of up to ``k`` items, each element of the
> stream equally likely to be included. One pass, O(k) memory. Deterministic for a
> fixed ``seed``.

## string_periodicity

### `borders(s)`  _function_

> Return the lengths of all borders of ``s`` (proper prefix == suffix), ascending.

### `count_occurrences(text, pattern)`  _function_

> Return the start indices of every occurrence of ``pattern`` in ``text`` (Z-function).

### `is_periodic(s)`  _function_

> True if ``s`` is a whole number (>= 2) of copies of a shorter block.

### `manacher_longest_palindrome(s)`  _function_

> Return a longest palindromic substring of ``s`` (Manacher's algorithm, ``O(n)``).
>
> Ties resolve to the earliest-starting longest palindrome. The empty string returns
> ``""``.

### `prefix_function(s)`  _function_

> Return the KMP prefix function ``pi``: the longest proper prefix = suffix length.

### `smallest_period(s)`  _function_

> Return the length of the smallest period ``p`` such that ``s`` repeats ``s[:p]``.
>
> A period ``p`` means ``s[i] == s[i - p]`` for all ``i >= p`` (the last block may be
> partial). For a string that is a whole number of copies this is the repeating unit; a
> string with no shorter period returns its own length. Empty string returns ``0``.

### `z_function(s)`  _function_

> Return the Z-array: ``z[i]`` = length of the longest common prefix of ``s`` and ``s[i:]``.
>
> ``z[0]`` is conventionally ``0``. Computed in ``O(n)`` by maintaining the rightmost
> matching segment ``[l, r]`` seen so far.

## strings

### `hamming_distance(a, b)`  _function_

> Hamming distance: positions at which equal-length sequences differ.
>
> Raises ``ValueError`` if the lengths differ (Hamming distance is undefined then).

### `kmp_search(text, pattern)`  _function_

> All start indices where ``pattern`` occurs in ``text`` (Knuth-Morris-Pratt).
>
> Linear-time ``O(len(text)+len(pattern))`` search using the failure function, so it
> never re-examines text characters. An empty pattern matches at every position
> ``0..len(text)``. Returns a list of indices.

### `levenshtein(a, b)`  _function_

> Levenshtein edit distance: minimum single-character insert/delete/substitutions.
>
> The number of one-character edits to turn ``a`` into ``b``. ``O(len(a)*len(b))`` time
> with two rolling rows. Symmetric; zero iff the sequences are equal.

### `longest_common_subsequence(a, b)`  _function_

> Longest common subsequence (not necessarily contiguous) as a list/str.
>
> Returns a longest sequence appearing in both ``a`` and ``b`` in order. Returns the
> same type as ``a`` when ``a`` is a ``str``. ``O(len(a)*len(b))``.

### `longest_common_substring(a, b)`  _function_

> Longest contiguous substring common to ``a`` and ``b``.
>
> Returns the substring (same type as ``a`` when ``a`` is a ``str``); an empty result
> if there is no common character. ``O(len(a)*len(b))``.

## structural_break

### `chow_test(x, break_index)`  _function_

> Chow F test for a break in the mean at ``break_index``.
>
> Compares the pooled residual sum of squares (one mean for the whole sample)
> against the sum from fitting separate means before and after the break. Under
> the no-break null the statistic is F(1, n - 2) distributed; a large value
> rejects. Returns ``(F_statistic, dof1, dof2)``.

### `cusum_break_detected(x, confidence=0.95)`  _function_

> True if the standardized CUSUM path breaches its confidence band.

### `cusum_mean(x, confidence=0.95)`  _function_

> Standardized CUSUM of deviations from the sample mean.
>
> Returns ``(cusum, boundary)`` where ``cusum[k]`` is the cumulative sum of
> demeaned observations through index ``k`` divided by ``sigma * sqrt(n)``, and
> ``boundary`` is the confidence threshold. Under a stable mean the standardized
> path is a Brownian bridge (it starts and ends at 0), so the supremum of its
> absolute value is compared to the Kolmogorov critical values -- a ``cusum``
> magnitude exceeding ``boundary`` flags a structural break in the mean.

## structural_credit

### `credit_spread(asset_value, debt_face, r, sigma, t)`  _function_

> Continuously-compounded credit spread of the risky debt over the risk-free.
>
> ``spread = -ln(D_risky / D e^{-r T}) / T`` -- the yield pickup of the risky
> debt over discounting the face at the risk-free rate. Non-negative, zero in the
> no-default limit, and rising with leverage and volatility.

### `distance_to_default(asset_value, debt_face, r, sigma, t)`  _function_

> Distance to default: ``d2`` (standard deviations of asset drift above debt).
>
> ``d2 = (ln(V/D) + (r - sigma^2/2) T) / (sigma sqrt(T))``. The number of
> standard deviations the log assets must fall to hit the default barrier;
> higher is safer, and ``Phi(-DD)`` is the default probability.

### `equity_value(asset_value, debt_face, r, sigma, t)`  _function_

> Equity as a call on the firm's assets struck at the debt face value.
>
> ``E = call(V, K=D, T)`` -- shareholders own the residual after repaying debt,
> a call on the assets. Increases with asset value and volatility.

### `equity_volatility(asset_value, debt_face, r, asset_vol, t)`  _function_

> Equity volatility implied by the asset volatility (Merton).
>
> From Ito's lemma on the equity call, ``sigma_E = (V/E) N(d1) sigma_V`` -- the
> equity is a levered claim, so its volatility exceeds the asset volatility by
> the delta-elasticity factor ``(V/E) N(d1)``. Rises as leverage rises.

### `physical_default_probability(asset_value, debt_face, mu, sigma, t)`  _function_

> Real-world probability of default ``P(V_T < D) = Phi(-DD)`` under drift ``mu``.
>
> The physical-measure analogue of :func:`risk_neutral_default_probability`,
> using the firm's actual expected asset return rather than the risk-free rate.
> In the KMV framework this maps (through an empirical calibration) to the
> expected default frequency. Setting ``mu = r`` recovers the risk-neutral
> probability; because ``mu > r`` for a risky firm, the physical default
> probability is below the risk-neutral one.

### `physical_distance_to_default(asset_value, debt_face, mu, sigma, t)`  _function_

> Real-world distance to default under the physical asset drift ``mu``.
>
> The risk-neutral :func:`distance_to_default` discounts at the risk-free rate;
> the physical measure uses the firm's actual expected asset return ``mu``:
>
>     DD = (ln(V/D) + (mu - sigma^2/2) T) / (sigma sqrt(T)).
>
> This is the Moody's-KMV distance to default -- the number of asset-return
> standard deviations between the current value and the default point. Setting
> ``mu = r`` recovers the risk-neutral figure; a higher expected return moves the
> firm further from default. Higher is safer.

### `risk_neutral_default_probability(asset_value, debt_face, r, sigma, t)`  _function_

> Risk-neutral probability of default ``P(V_T < D) = Phi(-d2)``.
>
> The chance the assets end below the debt face at maturity under the pricing
> measure. Rises with leverage (``D/V``), volatility, and horizon.

### `risky_debt_value(asset_value, debt_face, r, sigma, t)`  _function_

> Value of the risky debt: assets minus equity (``V - E``).
>
> By the accounting identity ``V = E + D_risky`` the debt is the firm value less
> the equity call. Below the risk-free discounted face; the gap is the credit
> risk.

### `solve_asset_value_and_vol(equity_value_obs, equity_vol_obs, debt_face, r, t, tol=1e-10, max_iter=500)`  _function_

> Recover the unobservable asset value and volatility (KMV two-equation solve).
>
> Given the observed equity value and equity volatility, jointly solves the
> Merton system
>
>     E = call(V, D, T),   sigma_E = (V/E) N(d1) sigma_V
>
> for ``(V, sigma_V)`` by fixed-point iteration: invert the equity-call for ``V``
> at the current ``sigma_V``, then update ``sigma_V`` from the equity-vol
> relation. Returns ``(asset_value, asset_vol)``. Round-trips with
> :func:`equity_value` and :func:`equity_volatility`.

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

### `phoenix_autocall_mc(S, t, r, sigma, observation_times, autocall_barrier, coupon_barrier, coupon, protection_barrier=None, principal=1.0, memory=True, q=0.0, n_paths=40000, seed=1234567)`  _function_

> Monte Carlo a Phoenix autocallable note.
>
> A Phoenix pays a coupon at each observation where the spot is at or above the
> ``coupon_barrier`` (typically below the autocall level). With ``memory=True``
> any coupons missed while below the barrier are paid retroactively the next time
> the barrier is met (snowball/memory feature). If the spot reaches
> ``autocall_barrier`` the note redeems early at par plus the coupon due. At
> maturity, unredeemed, the holder gets par unless the spot is below
> ``protection_barrier`` (down-and-in), taking the downside ``principal * S_T/S``.
>
> Returns the discounted Monte Carlo price. Memory raises the value versus no
> memory; a lower coupon barrier pays more often. Pure standard library.

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

## student_t

### `fit_df_from_kurtosis(excess_kurtosis)`  _function_

> Degrees of freedom implied by a sample's excess kurtosis (moment match).
>
> For a Student-t the excess kurtosis is ``6 / (df - 4)`` (finite only for
> ``df > 4``), so ``df = 4 + 6 / excess_kurtosis``. Requires positive excess
> kurtosis (fatter than normal); larger kurtosis implies fewer degrees of
> freedom (heavier tails).

### `fit_student_t(returns)`  _function_

> Fit a location-scale Student-t to a return sample by moment matching.
>
> Matches the sample mean, variance, and excess kurtosis: ``df`` from
> :func:`fit_df_from_kurtosis`, and the scale from ``variance = scale^2 df/(df-2)``.
> Returns ``(mean, scale, df)``. Requires ``df > 4`` (positive excess kurtosis).

### `student_t_expected_shortfall(mean, scale, df, confidence=0.95)`  _function_

> Parametric expected shortfall of a location-scale Student-t.
>
> Closed form ``ES = -mean + scale * (df + q^2)/(df - 1) * f(q)/(1 - confidence)``
> where ``q = t_ppf(1 - confidence, df)`` and ``f`` is the t density. Requires
> ``df > 1`` (finite mean); always at least the :func:`student_t_var`.

### `student_t_var(mean, scale, df, confidence=0.95)`  _function_

> Parametric VaR of a location-scale Student-t (positive loss magnitude).
>
> ``-(mean + scale * t_ppf(1 - confidence, df))`` -- heavier-tailed than the
> normal VaR for finite ``df``, converging to it as ``df -> inf``.

### `t_cdf(x, df)`  _function_

> Student-t cumulative distribution via the regularized incomplete beta.
>
> Uses the identity ``P(T <= x) = 1 - 0.5 I_{df/(df+x^2)}(df/2, 1/2)`` for
> ``x > 0`` and symmetry for ``x < 0``. Converges to :func:`norm_cdf` as
> ``df -> inf``.

### `t_pdf(x, df)`  _function_

> Student-t probability density with ``df`` degrees of freedom.
>
> ``f(x) = Gamma((df+1)/2) / (sqrt(df pi) Gamma(df/2)) (1 + x^2/df)^{-(df+1)/2}``.
> Symmetric about zero, heavier-tailed than the normal for finite ``df``.

### `t_ppf(p, df)`  _function_

> Student-t quantile (inverse CDF) by bisection on :func:`t_cdf`.
>
> Returns the ``x`` with ``t_cdf(x, df) = p``. Symmetric: ``t_ppf(1-p) =
> -t_ppf(p)``.

## sturm

### `isolate_real_roots(coeffs, a, b, max_depth=100)`  _function_

> Return disjoint sub-intervals of ``(a, b]``, each containing exactly one real root.
>
> Bisects the interval, using :func:`real_root_count` to decide which halves hold roots,
> until every returned interval isolates a single distinct real root. Each interval is a
> ``(lo, hi)`` pair.

### `real_root_count(coeffs, a, b)`  _function_

> Number of *distinct* real roots in the half-open interval ``(a, b]`` (Sturm's theorem).
>
> Exact for a square-free polynomial; repeated roots are counted once. ``a < b`` required.

### `sturm_sequence(coeffs)`  _function_

> Return the Sturm sequence of a real polynomial (list of coefficient lists).
>
> ``coeffs`` are highest-degree first. The sequence starts with the (square-free-agnostic)
> polynomial and its derivative; each subsequent term is the negated division remainder.

## style_analysis

### `style_analysis(fund_returns, index_returns, max_iter=5000, lr=None)`  _function_

> Returns-based style analysis: implied long-only index weights of a fund.
>
> ``fund_returns`` is the return series; ``index_returns`` is a list of index
> return series (one per style factor), each aligned with the fund. Returns a dict
> with ``weights`` (non-negative, summing to one), ``r_squared`` (fraction of fund
> variance explained by the style mix), and ``tracking_error`` (stdev of the
> unexplained residual). Solved by projected-gradient descent on the simplex.

## suffix_array

### `count_distinct_substrings(text)`  _function_

> Number of distinct non-empty substrings of ``text``.
>
> Equals ``sum(n - sa[r]) - sum(lcp)`` -- total suffix lengths minus the prefixes shared
> with the previous sorted suffix (which would be double-counted).

### `lcp_array(text, sa=None)`  _function_

> Kasai LCP array: ``lcp[r]`` = longest common prefix of ``sa[r]`` and ``sa[r-1]``.
>
> ``lcp[0]`` is 0 by convention. Runs in ``O(n)`` given the suffix array (computed if not
> supplied), exploiting that adjacent suffixes in text order lose at most one leading
> character of shared prefix between successive positions.

### `longest_repeated_substring(text)`  _function_

> Return a longest substring occurring at least twice (``""`` if none repeats).
>
> The answer is the text slice at the maximum LCP value: the deepest shared prefix
> between two adjacent sorted suffixes. Ties resolve to the first such position.

### `rank_array(sa)`  _function_

> Inverse of the suffix array: ``rank[i]`` is the sorted position of suffix ``i``.

### `substring_search(text, pattern, sa=None)`  _function_

> Return sorted start indices of every occurrence of ``pattern`` in ``text``.
>
> Binary-searches the suffix array for the block of suffixes that start with ``pattern``
> (``O(m log n)``). An empty pattern matches at every position.

### `suffix_array(text)`  _function_

> Return the suffix array of ``text`` -- suffix start positions in sorted order.
>
> ``result[r]`` is the starting index of the ``r``-th smallest suffix. Built by prefix
> doubling: sort by first character, then repeatedly refine ranks using pairs of ranks a
> power-of-two apart, so ``O(log n)`` rounds of ``O(n)`` counting each.

## surface

### `CalendarViolation(t_short: float, t_long: float, k: float, w_short: float, w_long: float) -> None`  _class_

> CalendarViolation(t_short: float, t_long: float, k: float, w_short: float, w_long: float)

### `SurfaceSlice(t: float, params: quantforge.svi.SVIParams, rmse: float) -> None`  _class_

> SurfaceSlice(t: float, params: quantforge.svi.SVIParams, rmse: float)

### `VolSurface(slices: List[quantforge.surface.SurfaceSlice])`  _class_

> A term structure of SVI smiles with calendar-arbitrage diagnostics.

## survival

### `kaplan_meier(times, events)`  _function_

> Kaplan-Meier product-limit survival estimate.
>
> Returns ``(event_times, survival)`` giving the step-function survival ``S(t)``
> at each distinct time where an event occurs; censoring times only shrink the risk
> set. ``S`` starts at 1, is non-increasing, and drops by the factor
> ``1 - d_i / n_i`` at each event time. With no censoring it equals
> ``1 - ECDF(t)``.

### `log_rank_test(times1, events1, times2, events2)`  _function_

> Log-rank (Mantel-Cox) test comparing two survival curves.
>
> At each distinct event time across the pooled sample, compares the observed
> events in group 1 with the number expected under the null of equal hazards
> (proportional to each group's share of the risk set), accumulating the
> observed-minus-expected and its hypergeometric variance. The statistic
>
>     chi2 = (sum (O1 - E1))^2 / sum V1
>
> is asymptotically chi-square(1). Returns ``(chi2, p_value)``; a small p-value
> rejects equal survival between the groups. Uses the chi-square survival function.

### `median_survival_time(times, events)`  _function_

> Median survival: the earliest time at which Kaplan-Meier ``S(t) <= 0.5``.
>
> Returns ``None`` if the curve never falls to ``0.5`` (survival stays above the
> median over the observed range, e.g. under heavy censoring).

### `nelson_aalen(times, events)`  _function_

> Nelson-Aalen cumulative-hazard estimate.
>
> Returns ``(event_times, cumulative_hazard)`` with ``H(t) = sum d_i / n_i`` over
> event times up to ``t``. Non-decreasing from 0; ``exp(-H(t))`` approximates the
> Kaplan-Meier survival (they agree closely when the per-step hazard is small).

### `restricted_mean_survival_time(times, events, tau)`  _function_

> Restricted mean survival time (RMST): area under KM up to horizon ``tau``.
>
> ``RMST(tau) = integral_0^tau S(t) dt`` where ``S`` is the Kaplan-Meier step
> function (``S = 1`` before the first event). The expected event time capped at
> ``tau`` -- a censoring-robust summary that, unlike the mean, is always defined
> even when the tail of the curve is not estimable. ``tau > 0``.

### `survival_at(times, events, query, estimator='km')`  _function_

> Evaluate the survival function at ``query`` from a fitted step curve.
>
> ``estimator`` is ``"km"`` (Kaplan-Meier) or ``"na"`` (``exp(-Nelson-Aalen)``).
> Returns the survival at the largest event time ``<= query`` (1 before the first
> event).

## svd

### `condition_number(A)`  _function_

> Spectral condition number ``sigma_max / sigma_min`` of ``A``.
>
> The ratio of the largest to smallest singular value; large means ``A`` is
> ill-conditioned (small perturbations blow up the solution). One for an orthogonal
> matrix; ``inf`` when ``A`` is singular (a zero singular value).

### `frobenius_norm(A)`  _function_

> Frobenius norm ``sqrt(sum a_ij^2)`` -- equivalently ``sqrt(sum sigma_k^2)``.

### `matrix_rank(A, rcond=1e-12)`  _function_

> Numerical rank: the number of singular values above ``rcond * sigma_max``.
>
> Counts the singular directions that carry real signal; singular values below the
> relative tolerance are treated as numerical zeros.

### `pseudo_inverse(A, rcond=1e-12)`  _function_

> Moore-Penrose pseudo-inverse ``A^+`` via the SVD.
>
> ``A^+ = V S^+ U'`` with the reciprocals of the singular values above
> ``rcond * s_max`` (smaller ones treated as zero). For a full-rank tall ``A`` this
> is ``(A' A)^{-1} A'``; applied to ``b`` it gives the minimum-norm least-squares
> solution. Returns the ``n x m`` pseudo-inverse.

### `spectral_norm(A)`  _function_

> Spectral (operator 2-) norm: the largest singular value of ``A``.

### `svd(A, tol=1e-14, max_sweeps=60)`  _function_

> One-sided Jacobi SVD of an ``m x n`` matrix (``m >= n``).
>
> Returns ``(U, s, V)`` where ``U`` is ``m x n`` with orthonormal columns, ``s`` is
> the length-``n`` list of singular values (descending, non-negative), and ``V`` is
> ``n x n`` orthogonal, such that ``A = U diag(s) V'``. Iteratively rotates column
> pairs until they are orthogonal.

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

## svi_jumpwing

### `SVIJumpWing(v: float, psi: float, p: float, c: float, vtilde: float) -> None`  _class_

> SVIJumpWing(v: float, psi: float, p: float, c: float, vtilde: float)

### `jumpwing_to_raw(jw: quantforge.svi_jumpwing.SVIJumpWing, t: float) -> quantforge.svi.SVIParams`  _function_

> Convert jump-wing parameters back to raw SVI at expiry ``t``.
>
> Inverts :func:`raw_to_jumpwing` with the Gatheral-Jacquier formulas. ``t`` is the
> expiry in years.

### `raw_to_jumpwing(params: quantforge.svi.SVIParams, t: float) -> quantforge.svi_jumpwing.SVIJumpWing`  _function_

> Convert raw SVI parameters to the jump-wing parameterization at expiry ``t``.
>
> Uses the closed-form Gatheral-Jacquier map. ``t`` is the expiry in years.

## swap

### `single_curve_par_swap_rate(discount, start, maturity, freq=2)`  _function_

> Par (fair) fixed rate that gives the swap zero value at inception.

### `swap_annuity(discount, start, maturity, freq=2)`  _function_

> Fixed-leg annuity ``sum_k tau_k P(t_k)`` (PV01 per unit rate).
>
> The present value of receiving 1 unit of rate on the fixed schedule.

### `vanilla_swap_value(discount, fixed_rate, start, maturity, freq=2, notional=1.0, payer=True)`  _function_

> Value of a fixed-for-floating swap off a single discount curve.
>
> A payer pays ``fixed_rate`` and receives float; its value is the float-leg PV
> (``P(start) - P(end)``) minus the fixed-leg PV (``fixed_rate * annuity``),
> times ``notional``. A receiver is the negative. Zero at the par swap rate.

## symmetric

### `elementary_symmetric(values)`  _function_

> Return ``[e_0, e_1, ..., e_n]`` for ``values`` (``e_0 = 1``).
>
> ``e_k`` is the sum of all products of ``k`` distinct entries. Built by the standard
> ``O(n^2)`` DP that multiplies in one value at a time.

### `elementary_to_power(e, kmax=None)`  _function_

> Convert ``[e_0, ..., e_n]`` to power sums ``[p_1, ..., p_kmax]`` (Newton's identities).
>
> ``p_k = sum_{i=1}^{k-1} (-1)^(i-1) e_i p_{k-i} + (-1)^(k-1) k e_k`` (with ``e_k = 0``
> for ``k > n``). Defaults ``kmax`` to ``n``.

### `poly_from_roots(roots)`  _function_

> Monic polynomial coefficients (highest-degree first) of ``prod (t - r)``.
>
> ``coeffs[k] = (-1)^k e_k``. Matches the coefficient order of
> :func:`quantforge.polyroots.polynomial_roots`.

### `power_sums(values, kmax)`  _function_

> Return ``[p_0, p_1, ..., p_kmax]`` where ``p_k = sum x_i^k`` (``p_0 = n``).

### `power_to_elementary(p)`  _function_

> Convert power sums ``[p_1, ..., p_n]`` to ``[e_0, ..., e_n]`` (Newton's identities).
>
> ``k e_k = sum_{i=1}^{k} (-1)^(i-1) e_{k-i} p_i``. Uses exact `Fraction` arithmetic so the
> division by ``k`` is exact for integer power sums.

## symplectic

### `leapfrog(force, q0, v0, dt, n_steps)`  _function_

> Leapfrog (kick-drift-kick) integration of ``q'' = force(q)`` (unit mass, velocity form).
>
> ``force(q)`` returns the acceleration. Returns ``(qs, vs)`` positions and velocities.
> Algebraically equivalent to :func:`velocity_verlet` with unit mass; kept as the
> velocity-space form common in N-body simulation.

### `velocity_verlet(force, q0, p0, mass, dt, n_steps)`  _function_

> Velocity-Verlet integration of ``q'' = force(q) / mass``.
>
> ``force(q)`` returns the force (accel * mass) as the same shape as ``q0``; ``p`` is
> the momentum ``mass * velocity``. Returns ``(qs, ps)``: the ``n_steps + 1`` position
> and momentum states (each a list). Second-order accurate and symplectic, so total
> energy stays bounded over long runs.

## t_copula_sample

### `student_t_copula_sample(correlation, df, n, seed=1234567)`  _function_

> Draw ``n`` samples from a Student-t copula.
>
> Parameters
> ----------
> correlation : list[list[float]]
>     Symmetric positive-definite correlation matrix (unit diagonal).
> df : int
>     Degrees of freedom (>= 1). Smaller ``df`` gives heavier joint tails; large
>     ``df`` approaches the Gaussian copula.
> n : int
>     Number of sample vectors.
> seed : int
>     Seed for the deterministic normal stream.
>
> Returns
> -------
> list[list[float]]
>     ``n`` vectors of uniforms in (0, 1) with the target rank correlation and
>     symmetric tail dependence set by ``df``.

## tail_dependence

### `exceedance_correlation(x, y, q=0.9, tail='upper')`  _function_

> Pearson correlation computed only on joint-tail observations.
>
> Selects the points where both margins breach the ``q`` quantile (upper tail)
> or fall below the ``1 - q`` quantile (lower tail), then returns the ordinary
> correlation of ``x`` and ``y`` on that subset -- the "correlations rise in the
> tails" diagnostic. Raises if fewer than two joint-tail points exist.

### `lower_tail_dependence(x, y, q=0.05)`  _function_

> Empirical lower tail-dependence coefficient at threshold ``q``.
>
> Estimates ``P(U <= q | V <= q)`` on the rank/uniform scale. Near 1 means joint
> downside extremes cluster. Returns 0 when no point falls below the threshold
> in ``y``.

### `upper_tail_dependence(x, y, q=0.95)`  _function_

> Empirical upper tail-dependence coefficient at threshold ``q``.
>
> Estimates ``P(U > q | V > q)`` on the rank/uniform scale: of the points whose
> ``y`` rank exceeds ``q``, the fraction whose ``x`` rank also exceeds ``q``.
> Near 1 means the two crash/spike together; near ``1 - q``-scaled independence
> means they do not. Returns 0 when no point exceeds the threshold in ``y``.

## theil_sen

### `theil_sen(x, y)`  _function_

> Theil-Sen robust regression slope and intercept.
>
> Parameters
> ----------
> x, y : sequence of float
>     Paired observations of equal length (at least 2). Pairs sharing an ``x``
>     value are skipped (their slope is undefined).
>
> Returns
> -------
> (slope, intercept) : (float, float)
>     The median pairwise slope and the median residual intercept. On exactly
>     collinear data this reproduces the generating line; under heavy-tailed
>     contamination it stays close to the clean fit where OLS is dragged away.

## thiele

### `thiele_coefficients(xs, ys)`  _function_

> Reciprocal-difference coefficients for Thiele's continued fraction.
>
> Returns the list ``rho`` whose entries are the leading inverse differences used by
> :func:`thiele_eval`. Length equals ``len(xs)``. Raises on duplicate nodes or a
> degenerate (zero) difference that stalls the recursion.

### `thiele_eval(xs, coeffs, x)`  _function_

> Evaluate Thiele's continued fraction with reciprocal-difference ``coeffs``.
>
> ``coeffs`` come from :func:`thiele_coefficients`. Evaluated from the bottom up; the
> ``rho[0]`` and ``rho[1]`` entries seed the fraction. Returns the rational
> interpolant's value at ``x``.

### `thiele_interpolate(xs, ys, x)`  _function_

> Convenience: build Thiele coefficients and evaluate at ``x`` in one call.

## trade_sign

### `lee_ready(prices, bids, asks)`  _function_

> Lee-Ready trade classification: quote rule with a tick-rule tiebreak.
>
> Trades away from the midpoint are signed by the quote rule; trades exactly at
> the midpoint are signed by the tick rule on the trade-price series. Returns a
> list of ``+1 / -1`` signs (no zeros). Aligned series of at least one trade.

### `quote_rule(prices, bids, asks)`  _function_

> Classify trades by their side of the prevailing bid-ask midpoint.
>
> ``+1`` above the midpoint (buyer-initiated), ``-1`` below, and ``0`` exactly at
> the midpoint (unclassified; :func:`lee_ready` resolves these). Aligned series.

### `tick_rule(prices)`  _function_

> Classify trades by the tick test on the trade-price series.
>
> ``+1`` on an uptick, ``-1`` on a downtick, and the previous sign carried forward
> on a zero tick (the first trade defaults to ``+1``). Needs at least one price.

## transfer_entropy

### `mutual_information(x, y, bins=8)`  _function_

> Mutual information ``I(X;Y)`` in nats, from a 2-D histogram.
>
> Non-negative, zero iff ``X`` and ``Y`` are independent, and symmetric:
> ``mutual_information(x, y) == mutual_information(y, x)``. ``bins`` sets the
> discretization resolution. Aligned series of at least two points.

### `transfer_entropy(source, target, bins=8)`  _function_

> Transfer entropy ``TE_{source -> target}`` in nats (lag-1, Schreiber).
>
> Measures how much the source's present reduces uncertainty about the target's
> next value beyond the target's own present. Built from the joint histogram of
> ``(target_{t+1}, target_t, source_t)``. Directional: run it both ways to see
> which series leads. Non-negative; near zero when the source carries no extra
> information about the target's future. Aligned series of at least three points.

## triangulate

### `ear_clipping_triangulate(polygon)`  _function_

> Triangulate a simple polygon by ear clipping.
>
> Returns a list of ``n - 2`` triangles, each a tuple of three ``(x, y)`` vertices. The
> polygon must be simple (non-self-intersecting) with vertices in order; both windings
> are accepted. The triangles' areas sum to the polygon's area.

### `is_clockwise(polygon)`  _function_

> Whether the polygon vertices are ordered clockwise (negative signed area).

### `is_convex_polygon(polygon)`  _function_

> Whether a simple polygon is convex (all turns the same direction).

### `signed_area(polygon)`  _function_

> Signed area of a polygon (shoelace): positive if counter-clockwise, negative if CW.

## tridiagonal

### `solve_cyclic_tridiagonal(lower, diag, upper, rhs)`  _function_

> Solve a *cyclic* tridiagonal system (periodic boundary) via Sherman-Morrison.
>
> Same layout as :func:`solve_tridiagonal`, but ``lower[0]`` is the corner entry coupling
> row 0 to column ``n-1``, and ``upper[n-1]`` couples row ``n-1`` to column 0. Uses two
> Thomas solves plus a rank-1 correction. Needs ``n >= 3``.

### `solve_tridiagonal(lower, diag, upper, rhs)`  _function_

> Solve a tridiagonal system by the Thomas algorithm; returns the solution list.
>
> ``diag`` is the main diagonal (length ``n``); ``lower[i]`` is the sub-diagonal entry in
> row ``i`` (``lower[0]`` unused); ``upper[i]`` the super-diagonal entry in row ``i``
> (``upper[n-1]`` unused). All four lists have length ``n``. Raises on a zero pivot (the
> matrix is singular or not diagonally dominant enough for plain elimination).

## trie

### `Trie(words=None)`  _class_

> A prefix tree over string keys (a multiplicity-free set with prefix queries).
>
> `insert`, `contains`, and `delete` manage membership; `starts_with` tests whether any
> key has a given prefix; `count_prefix` counts keys under it; `keys_with_prefix` lists
> them; and `longest_prefix_of` finds the longest stored key that prefixes a query string.

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

## turbulence

### `absorption_ratio(cov, n_factors=None)`  _function_

> Absorption ratio: variance share of the top principal components.
>
> Sums the ``n_factors`` largest eigenvalues of ``cov`` and divides by the total
> (the trace). ``n_factors`` defaults to about a fifth of the assets (Kritzman-Li's
> rule of thumb). In ``[0, 1]``: near ``n_factors / n`` when risk is spread evenly,
> toward 1 when a few factors dominate -- a high or rising ratio signals a fragile,
> tightly-coupled market.

### `turbulence(observation, mean, cov)`  _function_

> Financial turbulence: Mahalanobis distance of ``observation`` from ``mean``.
>
> ``d = (r - mu)' C^{-1} (r - mu)`` for a single return vector ``observation`` given
> the historical ``mean`` and covariance ``cov``. Non-negative; larger means a more
> unusual (stressed) cross-asset move. Under multivariate normality ``E[d]`` equals
> the number of assets.

### `turbulence_series(returns, mean=None, cov=None)`  _function_

> Turbulence for each row of a return panel (rows = periods, cols = assets).
>
> If ``mean``/``cov`` are omitted they are estimated in-sample from ``returns``.
> Returns one turbulence value per period; the average is close to the number of
> assets when the data is multivariate normal.

## tv_denoise

### `tv_denoise(y, lam)`  _function_

> Total-variation denoise ``y`` with regularization weight ``lam`` (Condat's method).
>
> Returns the exact minimizer of ``(1/2) sum (x - y)^2 + lam * sum |x_{k+1} - x_k|``.
> Larger ``lam`` produces flatter output with fewer, larger jumps; ``lam = 0`` returns
> ``y`` unchanged. The result is piecewise constant. ``O(n)`` time, ``O(1)`` extra
> state beyond the output.

### `tv_total_variation(x)`  _function_

> Total variation ``sum |x_{k+1} - x_k|`` of a sequence.
>
> The quantity the L1 penalty in :func:`tv_denoise` shrinks. Useful for confirming a
> denoised signal is flatter (lower total variation) than its noisy input.

## two_sat

### `TwoSat(n)`  _class_

> A 2-SAT instance over ``n`` boolean variables ``0 .. n-1``.
>
> Add clauses with `add_clause` (each literal is a variable index for the positive form or
> its bitwise complement ``~v`` for the negation) or the readable `add_or`. `solve` returns
> a satisfying list of booleans, or ``None`` if the formula is unsatisfiable.

## two_scale_rv

### `noise_variance_estimate(prices)`  _function_

> Estimate the microstructure-noise variance ``Var(eps)``.
>
> Under the noise model the fast-scale RV is dominated by noise, so
> ``Var(eps) ~= RV_fast / (2 * (n - 1))`` where ``n - 1`` is the number of
> returns. Consistent as the noise dwarfs the signal at the finest scale.

### `realized_variance_naive(prices)`  _function_

> Naive all-observations realized variance of a (log) price series.
>
> Upward-biased by ``2 n Var(eps)`` under i.i.d. microstructure noise.

### `two_scale_realized_variance(prices, K=None)`  _function_

> Two-scale realized variance (TSRV), robust to microstructure noise.
>
> Parameters
> ----------
> prices : sequence of float
>     Observed log-prices on a fine grid.
> K : int, optional
>     Number of subsampling grids for the slow scale. Defaults to
>     ``max(2, round(n ** (1/3)))``, the rate-optimal choice.
>
> Returns
> -------
> float
>     Bias-corrected estimate of the integrated variance. On noise-free data it
>     essentially reproduces the realized variance; under noise it is far less
>     biased than :func:`realized_variance_naive`.

## ukf

### `unscented_kalman_filter(observations, f, h, Q, R, x0, P0, alpha=0.001, beta=2.0, kappa=0.0)`  _function_

> Unscented Kalman filter over ``observations`` for nonlinear ``f`` and ``h``.
>
> ``f`` maps a length-``n`` state list to a length-``n`` list (transition); ``h`` maps the
> state to a length-``m`` list (measurement). Both take and return plain floats -- no autodiff
> or Jacobians. ``Q`` (n x n), ``R`` (m x m) covariances; ``x0`` (n), ``P0`` (n x n) initial
> mean/covariance. ``alpha``/``beta``/``kappa`` are the van der Merwe scaling parameters.
> Returns a dict with ``filtered_means`` and ``filtered_covariances``.

## valuation

### `capm_cost_of_equity(risk_free, beta, market_premium)`  _function_

> CAPM cost of equity ``risk_free + beta * market_premium``.
>
> The required return on equity given its market beta and the equity risk
> premium. Rises with beta.

### `gordon_growth_value(dividend_next, discount_rate, growth)`  _function_

> Gordon constant-growth value ``D_1 / (r - g)``.
>
> The present value of a perpetually-growing dividend. Requires ``r > g``
> (otherwise the sum diverges).

### `terminal_value(final_cashflow, discount_rate, growth)`  _function_

> Gordon terminal (continuing) value at the end of an explicit forecast.
>
> ``final_cashflow * (1 + growth) / (discount_rate - growth)`` -- the perpetuity
> value of cashflows beyond the forecast horizon. Requires ``r > g``.

### `two_stage_dcf(cashflows, discount_rate, terminal_growth)`  _function_

> Two-stage DCF: explicit cashflows plus a discounted terminal value.
>
> Discounts the explicit ``cashflows`` (periods 1..n) at ``discount_rate`` and
> adds the :func:`terminal_value` of the last cashflow discounted from period
> ``n``. Returns the enterprise/equity value. Requires ``discount_rate >
> terminal_growth``.

### `wacc(equity_value, debt_value, cost_of_equity, cost_of_debt, tax_rate)`  _function_

> Weighted-average cost of capital.
>
> ``E/V * ke + D/V * kd * (1 - tax)`` with ``V = E + D``. The after-tax blended
> discount rate; lies between the after-tax debt cost and the equity cost.

## vannavolga

### `VannaVolgaSmile(S, t, r_dom, r_for, atm, rr, bf, call_delta=0.25)`  _class_

> A vanna-volga FX smile built from ATM / RR / BF at one expiry.

### `pillar_vols(atm, rr, bf)`  _function_

> Return (sigma_25put, sigma_atm, sigma_25call) from ATM / RR / BF quotes.

## var_backtest

### `acerbi_szekely_es(losses, var_forecasts, es_forecasts, alpha=0.01)`  _function_

> Acerbi-Szekely (2014) Expected-Shortfall test statistic (their Z2).
>
> ``Z = (1 / (n alpha)) sum_t hit_t * loss_t / ES_t - 1``, where ``hit_t`` marks a
> VaR exception. Near zero when ES is well calibrated; positive when realized tail
> losses exceed the forecast ES (the model understates risk), negative when it
> overstates. Returns the scalar statistic.

### `christoffersen_cc(losses, var_forecasts, alpha=0.01)`  _function_

> Christoffersen conditional-coverage test: coverage and independence jointly.
>
> The sum of :func:`kupiec_pof` and :func:`christoffersen_independence` LR
> statistics, tested as a chi-square(2). Returns ``(LR, p_value)``.

### `christoffersen_independence(losses, var_forecasts)`  _function_

> Christoffersen independence test: exceptions should not cluster.
>
> Fits a first-order Markov chain to the exception indicator and tests that the
> probability of an exception does not depend on whether the previous period was an
> exception. Returns ``(LR, p_value)``, a chi-square(1) LR test.

### `kupiec_pof(losses, var_forecasts, alpha=0.01)`  _function_

> Kupiec proportion-of-failures (unconditional-coverage) test.
>
> ``alpha`` is the VaR tail probability (e.g. 0.01 for 99% VaR), so the expected
> exception rate is ``alpha``. Returns ``(LR, p_value)`` with ``LR`` a
> chi-square(1) likelihood ratio; a small p-value rejects the model's coverage.

## variance_ratio

### `variance_ratio(returns, q)`  _function_

> Lo-MacKinlay variance ratio ``VR(q)`` of a return series.
>
> Uses overlapping ``q``-period returns and the unbiased scaling factors from
> Lo-MacKinlay (1988). Returns 1 under a random walk, >1 for trending series,
> <1 for mean-reverting series.

### `variance_ratio_zstat(returns, q)`  _function_

> Heteroskedasticity-robust Lo-MacKinlay ``z`` statistic for ``VR(q) = 1``.
>
> Divides ``VR(q) - 1`` by the robust standard error assembled from the
> autocorrelations of squared demeaned returns. Asymptotically standard normal
> under the random-walk null; ``|z| > 1.96`` rejects at 5%.

## variance_tests

### `bartlett_test(*groups)`  _function_

> Bartlett's test for equal variance (likelihood ratio, chi-square).
>
> Returns a dict with the ``statistic``, ``df`` (``k-1``) and the ``p_value``. More
> powerful than Levene under normality but sensitive to non-normal tails. Requires at
> least two observations per group.

### `levene_test(*groups, center='median')`  _function_

> Levene / Brown-Forsythe test for equal variance across groups.
>
> ``center="median"`` (default) is the robust Brown-Forsythe form; ``"mean"`` is the
> original Levene. Returns a dict with the ``statistic`` (F), ``df`` ``(k-1, N-k)`` and
> the ``p_value``. A small p-value rejects equal variance.

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

## vector3

### `angle_between(a, b)`  _function_

> Angle between two non-zero vectors in radians, in ``[0, pi]``.
>
> Uses the numerically stable ``atan2(|a x b|, a . b)`` form in 3-D and the clamped
> ``acos`` of the normalized dot product otherwise.

### `cross(a, b)`  _function_

> Cross product of two 3-vectors: a vector perpendicular to both.

### `dot(a, b)`  _function_

> Dot product of two equal-length vectors.

### `norm(a)`  _function_

> Euclidean (L2) norm of a vector.

### `normalize(a)`  _function_

> Unit vector in the direction of ``a``; raises on the zero vector.

### `reflect(a, normal)`  _function_

> Reflect vector ``a`` about the plane with unit-normalizable ``normal``.
>
> ``a - 2 (a . n_hat) n_hat`` where ``n_hat`` is the normalized normal -- the standard
> mirror reflection (e.g. a ray bouncing off a surface).

### `vector_project(a, b)`  _function_

> Vector projection of ``a`` onto ``b`` (the component of ``a`` along ``b``).

### `vector_reject(a, b)`  _function_

> Vector rejection of ``a`` from ``b``: the component of ``a`` perpendicular to ``b``.

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

## vol_convert

### `black_to_normal_vol(forward, strike, t, sigma_black, is_call=True)`  _function_

> Convert a Black-76 (lognormal) vol to the equivalent Bachelier (normal) vol.
>
> Prices the option with Black-76 at ``sigma_black`` (zero rates, forward
> measure) and inverts the Bachelier model for the normal vol that reproduces it.

### `normal_to_black_vol(forward, strike, t, sigma_normal, is_call=True)`  _function_

> Convert a Bachelier (normal) vol to the equivalent Black-76 (lognormal) vol.
>
> Prices with Bachelier at ``sigma_normal`` and inverts Black-76. Requires
> positive forward and strike (a lognormal vol is undefined otherwise).

## vol_target

### `realized_annualized_vol(returns, periods_per_year=252)`  _function_

> Annualized realized volatility of a return series (sample std).

### `target_leverage(target_vol, realized_vol, max_leverage=None)`  _function_

> Volatility-target leverage ``target_vol / realized_vol`` (capped).
>
> Above one when realized vol is below target (lever up), below one when it is
> above (de-risk). Capped at ``max_leverage`` if given, and floored at zero.

### `vol_targeted_returns(returns, target_vol, lookback, periods_per_year=252, max_leverage=None)`  _function_

> Apply a rolling volatility-targeting overlay to a return series.
>
> For each period past the first ``lookback``, sizes the position at
> :func:`target_leverage` using the trailing ``lookback``-window annualized
> realized vol, and scales that period's return. Returns the overlaid return
> series (length ``len(returns) - lookback``). The overlay's realized vol sits
> near ``target_vol`` when the estimate tracks the true vol.

## volatility

### `EGarchParams(omega: float, alpha: float, beta: float, gamma: float) -> None`  _class_

> EGARCH(1,1) parameters (Nelson) on the log conditional variance.
>
> ``ln h_t = omega + beta ln h_{t-1} + alpha (|z| - E|z|) + gamma z``, with the
> standardized shock ``z = r/sqrt(h)`` and ``E|z| = sqrt(2/pi)`` for a Gaussian.
> Modelling the *log* variance means ``h`` is positive for any parameters (no
> constraints), and ``gamma < 0`` gives the leverage effect (negative shocks
> raise volatility more). Stationary when ``|beta| < 1``.

### `GJRGarchParams(omega: float, alpha: float, beta: float, gamma: float) -> None`  _class_

> GJR-GARCH(1,1,1) parameters with a leverage (asymmetry) term.
>
> ``h_t = omega + (alpha + gamma * I[r_{t-1} < 0]) r_{t-1}^2 + beta h_{t-1}``.
> ``gamma > 0`` makes negative shocks raise volatility more than positive ones
> (the leverage effect). Stationary when ``alpha + beta + 0.5 gamma < 1``.

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

### `egarch_forecast(params: quantforge.volatility.EGarchParams, last_return, last_variance, horizon=1, periods_per_year: int = 252)`  _function_

> Forecast annualized volatility ``horizon`` steps ahead under EGARCH.
>
> One step uses :func:`egarch_variance`; beyond that the *log* variance
> mean-reverts to its unconditional level ``omega/(1-beta)`` at rate ``beta`` per
> step (the shock terms are mean-zero), and the result is exponentiated and
> annualized. Returns the annualized volatility.

### `egarch_variance(params: quantforge.volatility.EGarchParams, last_return, last_variance)`  _function_

> One-step-ahead EGARCH conditional variance (always positive).
>
> Computes the standardized shock from the last return and variance and applies
> the log-variance recursion, exponentiating back to a variance. Positive for
> any parameters.

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

### `fit_har_rv(realized_variance, weekly=5, monthly=22)`  _function_

> Fit the HAR-RV (Corsi 2009) model to a realized-variance series.
>
> Regresses ``RV_{t+1}`` on the previous day's RV, the trailing ``weekly``-day
> average, and the trailing ``monthly``-day average -- a parsimonious long-memory
> model. Returns ``(beta0, beta_day, beta_week, beta_month)``. Recovers the true
> coefficients on data generated from the model.

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

### `gjr_garch_forecast(params: quantforge.volatility.GJRGarchParams, last_return, last_variance, horizon=1, periods_per_year: int = 252)`  _function_

> Forecast annualized volatility ``horizon`` steps ahead under GJR-GARCH.
>
> One step uses :func:`gjr_garch_variance`; beyond that the expected variance
> mean-reverts to the long-run level at the leverage-adjusted persistence
> ``alpha + beta + 0.5 gamma`` per step. Returns the annualized volatility.

### `gjr_garch_variance(params: quantforge.volatility.GJRGarchParams, last_return, last_variance)`  _function_

> One-step-ahead GJR-GARCH conditional variance.
>
> Adds the leverage term ``gamma`` to the ARCH coefficient when the last return
> was negative, so a down move feeds more into next-period variance than an up
> move of the same size.

### `har_rv_forecast(coeffs, recent_rv, weekly=5, monthly=22)`  _function_

> One-step HAR-RV forecast from the fitted coefficients and recent RV.
>
> ``coeffs`` is ``(beta0, beta_day, beta_week, beta_month)``; ``recent_rv`` is the
> trailing realized-variance history (at least ``monthly`` points). Forms the
> day / week / month averages and applies the linear model.

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

## von_mises

### `bessel_i0(x)`  _function_

> Modified Bessel function of the first kind, order 0, ``I0(x)``.
>
> Uses the Abramowitz & Stegun polynomial approximations (accurate to ~1e-7), the same
> ones used for the von Mises normalizing constant.

### `bessel_i1(x)`  _function_

> Modified Bessel function of the first kind, order 1, ``I1(x)`` (A&S approximation).

### `von_mises_fit(angles, tol=1e-10, max_iter=100)`  _function_

> Maximum-likelihood fit of a von Mises distribution to ``angles`` (radians).
>
> The MLE mean direction is the sample :func:`circular_mean`; the concentration
> ``kappa`` solves ``I1(kappa)/I0(kappa) = R`` (the mean resultant length), found here by
> Newton's method with a standard closed-form seed. Returns ``(mu, kappa)``.

### `von_mises_pdf(theta, mu, kappa)`  _function_

> Von Mises density at angle ``theta`` with mean ``mu`` and concentration ``kappa``.
>
> ``exp(kappa cos(theta - mu)) / (2 pi I0(kappa))``. ``kappa >= 0``; ``kappa = 0`` gives
> the uniform density ``1 / (2 pi)``. Integrates to 1 over any ``2*pi`` interval.

## voronoi

### `delaunay_neighbors(points)`  _function_

> Return the site-adjacency graph ``{i: set(neighbor indices)}`` (Delaunay edges).
>
> Two sites are neighbours iff they share a Delaunay edge -- equivalently, iff their
> Voronoi cells share an edge. The relation is symmetric.

### `nearest_site(points, query)`  _function_

> Return the index of the site nearest to ``query`` (the Voronoi cell it falls in).

### `voronoi_vertices(points)`  _function_

> Return the Voronoi vertices: the circumcenter of each Delaunay triangle.
>
> Each vertex is equidistant from the three sites of its triangle. Returns a list of
> ``(x, y)`` points (one per Delaunay triangle).

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

## walsh_hadamard

### `and_convolve(a, b)`  _function_

> AND convolution: ``c[k] = sum_{i & j = k} a[i] b[j]`` (via superset zeta transform).

### `fwht(a)`  _function_

> In-place-style fast Walsh-Hadamard transform (unnormalized). Returns a new list.
>
> Length must be a power of two. ``ifwht(fwht(a)) == a`` after the ``1/n`` scaling.

### `ifwht(a)`  _function_

> Inverse Walsh-Hadamard transform (applies the ``1/n`` normalization).

### `or_convolve(a, b)`  _function_

> OR convolution: ``c[k] = sum_{i | j = k} a[i] b[j]`` (via subset zeta transform).

### `xor_convolve(a, b)`  _function_

> XOR convolution: ``c[k] = sum_{i ^ j = k} a[i] b[j]``.
>
> ``a`` and ``b`` must share a power-of-two length. Integer inputs give exact integer
> output (the ``1/n`` division is exact because the Walsh transform's inverse sum is
> divisible by ``n`` for integer data).

## wasserstein

### `wasserstein1_sorted(x, y)`  _function_

> 1-Wasserstein distance for two equal-length samples (sorted-difference form).
>
> ``W_1 = mean_i |x_(i) - y_(i)|`` after sorting each sample. Requires equal lengths;
> use :func:`wasserstein_distance` for unequal sizes.

### `wasserstein_distance(x, y, p=1)`  _function_

> p-Wasserstein distance between two 1-D samples of any sizes.
>
> Merges the empirical CDFs and integrates ``|F_x^{-1}(u) - F_y^{-1}(u)|^p`` over
> ``u in [0, 1]`` via the standard interval decomposition (SciPy's approach). ``p=1``
> is earth-mover distance; ``p=2`` the quadratic transport cost. Returns a true metric
> for ``p=1``. Both samples must be non-empty.

## wave_equation

### `wave_equation(u0, v0, c, dx, dt, n_steps, left=0.0, right=0.0)`  _function_

> Explicit finite-difference evolution of ``u_tt = c^2 u_xx``.
>
> ``u0`` is the initial displacement profile (grid values), ``v0`` the initial
> velocity (same length; pass zeros for a plucked start). Fixed ends at ``left`` /
> ``right``. Returns the displacement after ``n_steps`` steps. Requires the CFL
> condition ``c dt / dx <= 1`` for stability (raises otherwise).

## wavelet

### `haar_dwt(x, levels=None)`  _function_

> Multilevel Haar wavelet transform.
>
> Returns ``(approx, details)`` where ``approx`` is the final coarse-approximation
> list and ``details`` is a list (finest level first) of the detail-coefficient
> lists at each level. ``levels`` defaults to the maximum ``log2(len(x))``. Length
> must be a power of two.

### `haar_idwt(approx, details)`  _function_

> Invert :func:`haar_dwt`, reconstructing the original signal.
>
> Takes the coarse approximation and the per-level detail lists (finest first) and
> returns the reconstructed series. Exact up to floating error.

### `wavelet_energy(x, levels=None)`  _function_

> Fraction of signal energy in each Haar detail level and the coarse approximation.
>
> Returns a dict with ``detail`` (a list of energy fractions, finest level first)
> and ``approx`` (the coarse-approximation energy fraction). The fractions sum to
> one because the Haar transform is orthonormal (Parseval). A smooth series
> concentrates energy in the approximation; a noisy one spreads it into the fine
> details.

## wavelet_denoise

### `hard_threshold(x, lam)`  _function_

> Hard-threshold (keep or kill): ``x`` if ``|x| > lam`` else ``0``.

### `mad_sigma(detail)`  _function_

> Robust noise-scale estimate from detail coefficients.
>
> ``sigma = median(|d|) / 0.6745`` -- the median absolute deviation rescaled to
> match the standard deviation of a Gaussian. Robust to the few large
> (signal-bearing) coefficients that would inflate a plain standard deviation.

### `soft_threshold(x, lam)`  _function_

> Soft-threshold (shrink toward zero): ``sign(x) * max(|x| - lam, 0)``.

### `universal_threshold(n, sigma)`  _function_

> VisuShrink universal threshold ``sigma * sqrt(2 log n)``.

### `wavelet_denoise(x, levels=None, mode='soft', threshold=None)`  _function_

> Denoise a signal by Haar wavelet shrinkage.
>
> Transforms ``x`` with :func:`haar_dwt`, shrinks every detail coefficient with the
> ``soft`` (default) or ``hard`` rule at the given ``threshold``, and inverts. If
> ``threshold`` is None the VisuShrink universal threshold is used, with the noise
> scale estimated by :func:`mad_sigma` from the finest detail level. The coarse
> approximation is left untouched (it carries the trend, not noise). Length must be
> a power of two. Returns the reconstructed, denoised signal.

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

### `degree_day_digital(expected_index, strike, sigma, r, expiry, payout, is_call=True)`  _function_

> Digital (binary) degree-day option: fixed ``payout`` on a strike breach.
>
> Under the normal (Bachelier) index model with mean ``expected_index`` and
> standard deviation ``sigma``, a call digital pays ``payout`` if the index
> finishes above ``strike``, a put digital if below. The value is the discounted
> breach probability times the payout:
>
>     call = e^{-rT} payout * Phi((mu - K) / sigma),
>     put  = e^{-rT} payout * Phi((K - mu) / sigma).
>
> Call and put digitals sum to ``e^{-rT} payout`` (one of them always pays).

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

## weighted_reservoir

### `weighted_reservoir_sample(items, weights, k, seed=1234567)`  _function_

> Sample ``k`` distinct items with weight-proportional probability (Efraimidis-Spirakis).
>
> Assigns each item the key ``u^(1/w)`` for a uniform ``u`` and keeps the ``k`` largest
> keys -- one pass, ``O(n log k)`` memory ``O(k)``. Returns a list of the chosen items
> (order not significant). ``k`` is clamped to the number of positive-weight items.

### `weighted_sample_with_replacement(items, weights, k, seed=1234567)`  _function_

> Sample ``k`` items *with* replacement, weight-proportional (cumulative search).
>
> Each of the ``k`` draws is independent, so an item can appear multiple times. Returns
> a list of length ``k``. Weights must be non-negative and not all zero.

## weighted_stats

### `weighted_mean(values, weights)`  _function_

> Weighted arithmetic mean ``sum(w x) / sum(w)``.
>
> ``weights`` must be non-negative and not all zero; they need not sum to 1.

### `weighted_median(values, weights)`  _function_

> Weighted median: the weighted 0.5-quantile.

### `weighted_quantile(values, weights, q)`  _function_

> Weighted ``q``-quantile (``q`` in ``[0, 1]``) by the cumulative-weight method.
>
> Sorts by value, forms the normalized cumulative weight at each point (midpoint
> convention), and linearly interpolates to the target ``q``. ``q=0.5`` is the weighted
> median. Reduces to the ordinary quantile when weights are equal.

### `weighted_std(values, weights, unbiased=True)`  _function_

> Weighted standard deviation: square root of :func:`weighted_variance`.

### `weighted_variance(values, weights, unbiased=True)`  _function_

> Weighted variance about the weighted mean.
>
> With ``unbiased=True`` applies the reliability-weight correction
> ``V1 / (V1^2 - V2)`` where ``V1 = sum(w)`` and ``V2 = sum(w^2)`` (reduces to the
> ``1/(n-1)`` factor for equal weights); with ``unbiased=False`` divides by ``V1``
> (the population form). ``weights`` non-negative, not all zero.

## wilcoxon

### `sign_test(x, mu0=0.0, y=None)`  _function_

> Sign test that the (paired) sample median equals ``mu0``.
>
> Counts observations above ``mu0`` (or ``x_i > y_i`` for paired data), dropping
> exact ties, and tests against Binomial(n, 1/2). Returns a dict with ``n_plus``,
> ``n`` (non-tied count), the exact two-sided binomial ``p_value`` and the
> normal-approximation ``z``.

### `wilcoxon_signed_rank_test(x, mu0=0.0, y=None)`  _function_

> Wilcoxon signed-rank test that the (paired) sample is centered at ``mu0``.
>
> If ``y`` is given the test runs on the paired differences ``x_i - y_i``; otherwise
> on ``x_i - mu0``. Zero differences are dropped (Wilcoxon's convention); tied
> absolute values receive average ranks. Returns a dict with ``statistic`` W (the
> positive-rank sum), the ``z`` normal approximation (continuity-corrected, with the
> tie correction to the variance) and the two-sided ``p_value``.

## windows

### `apply_window(x, window='hann')`  _function_

> Multiply signal ``x`` by a named window (or a precomputed window list).
>
> ``window`` is ``"hann"``, ``"hamming"``, ``"blackman"``, ``"bartlett"``,
> ``"rectangular"``, or a list of the same length as ``x``. Returns the tapered
> signal.

### `bartlett(n)`  _function_

> Bartlett (triangular) window, zero at both ends.

### `blackman(n)`  _function_

> Blackman window: three-term cosine, very low sidelobes at a wider main lobe.

### `hamming(n)`  _function_

> Hamming window: ``0.54 - 0.46 cos(2 pi k/(n-1))`` (lower first sidelobe than Hann).

### `hann(n)`  _function_

> Hann (raised-cosine) window of length ``n``: ``0.5(1 - cos(2 pi k/(n-1)))``.

### `rectangular(n)`  _function_

> Rectangular (boxcar) window: all ones -- no tapering.

## wls

### `generalized_least_squares(X, y, cov, add_intercept=True)`  _function_

> Generalized least squares for a known error covariance ``cov`` (Sigma).
>
> Whitens the system with the Cholesky factor of ``Sigma`` (``Sigma = L L'``), so
> ``L^{-1} y = L^{-1} X beta + white noise``, then applies OLS. Returns a dict with
> ``coefficients`` and ``std_errors``. A diagonal ``cov`` reproduces weighted least
> squares with ``weights = 1 / diag(cov)``.

### `weighted_least_squares(X, y, weights, add_intercept=True)`  _function_

> Weighted least squares: minimize ``sum w_t (y_t - x_t beta)^2``.
>
> ``weights`` is a length-``n`` list of non-negative weights (larger = more trusted).
> Returns a dict with ``coefficients``, ``std_errors`` (using the weighted residual
> variance), ``residuals`` and ``r_squared`` (weighted). Equal weights reproduce OLS.

## xor_basis

### `XorBasis(values=None)`  _class_

> A GF(2) linear basis maintained by leading bit.
>
> ``insert(x)`` adds an integer to the span (reducing it against the current basis);
> ``max_xor``/``min_xor`` return the extremal subset-XOR, ``can_represent`` tests
> membership, ``rank`` is the basis size, ``count_distinct`` is ``2^rank``, and
> ``kth_smallest`` indexes the sorted reachable values. Non-negative integers only.

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

## zeta

### `dirichlet_eta(s)`  _function_

> Dirichlet eta ``eta(s) = sum (-1)^{n-1} n^-s = (1 - 2^{1-s}) zeta(s)`` for ``s > 0``.
>
> Converges (conditionally) for all ``s > 0``; accelerated here by van Wijngaarden /
> alternating-series transformation for robustness near ``s -> 0``.

### `riemann_zeta(s, terms=20, corrections=10)`  _function_

> Riemann zeta ``zeta(s)`` for real ``s > 0``, ``s != 1`` (Euler-Maclaurin).
>
> For ``s > 1`` uses Euler-Maclaurin directly; for ``0 < s < 1`` uses the Dirichlet-eta
> relation ``zeta(s) = eta(s) / (1 - 2^{1-s})``. Raises at the pole ``s = 1``.
