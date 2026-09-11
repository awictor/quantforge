"""QuantForge: fast, dependency-free options pricing and risk engine.

Closed-form Black-Scholes-Merton pricing, full analytic Greeks, robust
implied-volatility solving, and a binomial lattice for American exercise.

Everything is pure Python (standard library only) so it drops into any
runtime with no compiled dependencies.
"""

from .bsm import (
    OptionType,
    price,
    call_price,
    put_price,
    delta,
    gamma,
    vega,
    theta,
    rho,
    epsilon,
    greeks,
    Greeks,
)
from .implied import implied_volatility, implied_vol_smile
from .binomial import american_price
from .leisen_reimer import (
    leisen_reimer_price, leisen_reimer_greeks, leisen_reimer_american_accel,
)
from .trinomial import trinomial_price, richardson_american
from .lsm import (
    bermudan_lsm, bermudan_lsm_local_vol, bermudan_lsm_greeks,
    bermudan_max_call_lsm, bermudan_max_call_lsm_greeks, bermudan_spread_lsm,
    bermudan_spread_lsm_greeks, bermudan_min_put_lsm,
    bermudan_min_put_lsm_greeks, bermudan_basket_lsm,
    bermudan_basket_lsm_greeks,
)
from .mlmc import mlmc_asian
from .perpetual import perpetual_american, perpetual_exercise_boundary
from .forward import implied_forward, ForwardResult, dividend_curve
from .forwardstart import (
    forward_start_price, forward_start_greeks, cliquet_price, cliquet_greeks,
)
from .chooser import chooser_option, chooser_option_greeks
from .compound import compound_option, compound_option_greeks
from .quanto import (
    quanto_option, compo_option, quanto_option_greeks, compo_option_greeks,
)
from .displaced import (
    displaced_diffusion_price, displaced_diffusion_greeks, displaced_implied_shift,
)
from .hedging import StickyRule, smile_delta, smile_delta_from_smile, skew_slope
from .lookback import (
    floating_strike_lookback, fixed_strike_lookback, lookback_greeks,
)
from .heston import heston_price, heston_smile
from .heston_mc import (
    heston_qe_mc, heston_cv_mc, heston_pathwise_delta, heston_mc_greeks,
)
from .bates import bates_price, bates_greeks, bates_smile
from .double_heston import (
    double_heston_price, double_heston_greeks, double_heston_smile,
)
from .rough_heston import (
    rough_heston_price, rough_heston_greeks, rough_heston_smile,
)
from .double_heston_calib import calibrate_double_heston
from .heston_calib import calibrate_heston
from .lsv import calibrate_leverage as calibrate_lsv_leverage
from .pde import (
    crank_nicolson_price, crank_nicolson_greeks, crank_nicolson_barrier,
    crank_nicolson_digital, crank_nicolson_no_touch,
)
from .pde2d import (
    adi_two_asset, adi_spread_option, adi_two_asset_cs, adi_two_asset_american,
)
from .pde_asian import asian_pde_price
from .kou import kou_price, kou_greeks, kou_smile
from .cgmy import cgmy_price, cgmy_greeks, cgmy_smile
from .nig import nig_price, nig_greeks, nig_smile
from .meixner import meixner_price, meixner_greeks, meixner_smile
from .levycalib import calibrate_levy_smile, levy_psi
from .levysurface import LevySurface, LevyCalendarViolation
from .carrmadan import (
    levy_price, carr_madan_strip, carr_madan_smile_strip, cos_price, cos_greeks,
)
from .rbergomi import (
    rbergomi_price, rbergomi_smile, rbergomi_price_cv, rbergomi_smile_cv,
    rbergomi_greeks_cv,
)
from .bachelier import (
    bachelier_price, bachelier_delta, bachelier_gamma, bachelier_vega,
    bachelier_implied_vol, bachelier_greeks,
)
from .variancegamma import (
    variance_gamma_price, variance_gamma_greeks, variance_gamma_smile,
)
from .hedgesim import simulate_delta_hedge, HedgeResult
from .merton import merton_jump_price, merton_jump_greeks, merton_smile
from .density import (
    risk_neutral_density, risk_neutral_cdf, price_from_density, density_total_mass,
)
from .varswap import (
    variance_swap_strike, volatility_swap_strike, variance_swap_from_smile,
    corridor_variance_swap_from_smile, gamma_swap_from_smile,
    forward_variance_swap_from_smile, variance_term_structure,
)
from .vix import (
    vix_from_chain, vix_from_smile, svix_from_smile, equity_premium_lower_bound,
)
from .bkm import bkm_moments_from_smile, skew_swap_from_smile
from .vrp import realized_variance, variance_risk_premium
from .mc_greeks import (
    lr_greeks, pathwise_delta, lr_digital_delta, asian_pathwise_vega,
    mixed_gamma, smoothed_digital_delta, barrier_lr_delta, lr_digital_greeks,
)
from .moment_premium import moment_risk_premia
from .rnd import (
    risk_neutral_density_from_smile, density_grid_from_smile,
    price_payoff_from_density, risk_neutral_cdf_from_smile,
    risk_neutral_quantile_from_smile, risk_neutral_var_from_smile,
    risk_neutral_cvar_from_smile,
)
from .density_metrics import (
    tail_probability, density_entropy, expected_shortfall,
    kl_divergence_smiles, wasserstein_smiles,
)
from .density_var import density_var_es
from .multiasset import (
    exchange_option, spread_option, spread_option_bs, basket_option,
    best_of_call, worst_of_call,
    best_of_call_closed, worst_of_call_closed,
    best_of_put_closed, worst_of_put_closed, two_asset_digital,
    two_asset_asset_or_nothing, correlation_option, two_asset_gap_option,
    two_asset_digital_greeks,
    exchange_greeks, spread_greeks, spread_option_bs_greeks, basket_greeks,
    rainbow_greeks,
    implied_spread_correlation, implied_spread_correlation_bs,
)
from .strategy import (
    payoff_at_expiry, payoff_profile, break_evens,
    vertical_spread, straddle, strangle, risk_reversal, butterfly, iron_condor,
    ratio_spread, backspread, strategy_report,
)
from .localvol import dupire_local_vol, local_vol_from_implied, sabr_local_vol
from .spline import CubicSpline, SmileSpline
from .rates import (
    CapletPeriod, caplet_price, cap_price, floor_price, collar_price,
    caplet_floorlet_parity, annuity, swaption_price, swaption_parity,
    caplet_greeks, cap_greeks, floor_greeks, swaption_greeks,
)
from .vasicek import (
    zero_coupon_bond, zero_coupon_yield, bond_option,
    bond_option_greeks as vasicek_bond_option_greeks,
    bond_greeks as vasicek_bond_greeks,
    coupon_bond_option as vasicek_coupon_bond_option,
    swaption as vasicek_swaption,
    caplet as vasicek_caplet, floorlet as vasicek_floorlet,
    cap as vasicek_cap, floor as vasicek_floor,
)
from .cir import (
    cir_zero_coupon_bond, cir_zero_coupon_yield, cir_bond_greeks,
    cir_bond_option, cir_coupon_bond_option, cir_swaption,
    cir_caplet, cir_floorlet, cir_cap, cir_floor,
)
from .holee import (
    holee_zero_coupon_bond, holee_zero_coupon_yield, holee_bond_greeks,
    holee_bond_option, holee_coupon_bond_option, holee_swaption,
    holee_caplet, holee_floorlet, holee_cap, holee_floor,
)
from .overhedge import (
    Overhedge, digital_call_overhedge, digital_put_overhedge, overhedge_payoff,
)
from .dv01 import KeyRateDV01, key_rate_dv01
from .bookgreeks import (
    BookSecondOrder, book_second_order, ThetaCarry, theta_carry_report,
    BumpGreeks, book_bump_greeks,
)
from .attribution import PnLAttribution, attribute_pnl, CarryRoll, carry_roll_pnl
from .income import IncomeMetrics, covered_call, cash_secured_put
from .vegabucket import VegaBuckets, vega_buckets
from .qmc import halton, european_qmc
from .sobol import (
    Sobol, brownian_bridge_path, sobol_european, sobol_european_rqmc, sobol_asian,
    sobol_asian_rqmc, sobol_lookback_rqmc, sobol_fixed_lookback_rqmc,
    sobol_barrier_rqmc, sobol_barrier_digital_rqmc, sobol_autocallable_rqmc,
    sobol_double_knockout_rqmc, sobol_cliquet_rqmc, sobol_average_strike_rqmc,
    sobol_geometric_asian_rqmc, sobol_parisian_rqmc,
    sobol_arithmetic_asian_rqmc,
)
from .correlation import (
    implied_correlation, index_vol_from_correlation, dispersion_basket_vol,
    correlation_term_structure,
)
from .cev import cev_price, cev_greeks, noncentral_chisq_cdf
from .sizing import (
    delta_hedge_shares, neutralize, vega_neutral_quantity, gamma_neutral_quantity,
    kelly_fraction_binary, kelly_fraction_continuous, kelly_growth_rate,
)
from .gramcharlier import (
    corrado_su_call, corrado_su_price, realized_skewness, realized_excess_kurtosis,
    calibrate_corrado_su, corrado_su_implied_vol, corrado_su_smile,
)
from .portfolio import Contract, Position, BookRisk, Book, price_book
from .svi import (
    SVIParams, calibrate_svi, svi_g, svi_butterfly_arbitrage, svi_is_butterfly_free,
    lee_wing_slopes, lee_bounds_ok, svi_repair_butterfly,
    svi_local_variance, svi_surface_local_vol, calibrate_svi_from_prices,
    svi_variance_swap_strike, svi_vix, svi_svix, svi_bkm_moments, svi_density,
)
from .ssvi import (
    SSVIParams, ssvi_phi, ssvi_total_variance, calibrate_ssvi,
    ssvi_butterfly_free, ssvi_calendar_free, ssvi_is_arbitrage_free,
    ssvi_local_variance, ssvi_local_vol_from_params,
    ssvi_local_vol_fn, ssvi_reprice_mc, calibrate_ssvi_arbitrage_free,
    ssvi_variance_swap_strike, ssvi_vix, ssvi_svix, ssvi_density,
    ssvi_bkm_moments,
)
from .sabr import (
    SABRParams, sabr_vol, calibrate_sabr, calibrate_sabr_lm,
    sabr_sensitivities, sabr_jacobian, sabr_option_greeks,
    sabr_variance_swap_strike, sabr_vix, sabr_bkm_moments,
    sabr_density, sabr_butterfly_arbitrage, sabr_is_arbitrage_free,
    sabr_repair_butterfly,
)
from .vannavolga import VannaVolgaSmile, pillar_vols
from .fxdelta import (
    atm_dns_strike, strike_from_delta, delta_from_strike, rr_bf_to_pillars,
)
from .volcube import VolCube
from .cms import (
    cms_adjustment_standard, cms_rate_convexity_replication, cms_rate,
    cms_adjustment_greeks,
)
from .surface import VolSurface, SurfaceSlice, CalendarViolation
from .exotics import (
    cash_or_nothing, asset_or_nothing, digital_greeks, barrier_option, barrier_greeks,
    geometric_asian, geometric_asian_greeks, arithmetic_asian, asian_greeks,
    one_touch, no_touch,
    gap_option, gap_option_greeks, power_option, power_option_greeks,
    barrier_rebate, Barrier,
)
from .montecarlo import (
    MCResult, european_mc, european_cv_mc, european_is_mc,
    european_is_adaptive_mc, european_stratified_mc, spread_option_lhs_mc,
    basket_option_lhs_mc, replicated_mc, digital_is_mc, arithmetic_asian_mc,
    capped_cliquet_mc,
    barrier_mc, barrier_digital_mc, parisian_barrier_mc, local_vol_mc,
    average_strike_asian_mc, autocallable_mc, double_knockout_mc,
)
from .risk import VaRResult, parametric_var, historical_var, montecarlo_var
from .volatility import (
    close_to_close, ewma_vol, parkinson, garman_klass, rogers_satchell,
    yang_zhang, vol_report, VolReport, vol_cone, VolConePoint,
    fit_garch, garch_forecast, GarchParams,
)
from .greeks2 import vanna, vomma, volga, charm, veta, speed, zomma, color
from .scenario import ScenarioGrid, stress_grid, spot_ladder
from .american import (
    bjerksund_stensland, bjerksund_stensland_greeks, early_exercise_premium,
    bjerksund_stensland_1993,
)
from .kim import (
    kim_american_put, kim_american_call, kim_exercise_boundary, kim_put_greeks,
)
from .baw import baw_american
from .andreasenhuge import (
    andreasen_huge_prices, andreasen_huge_smile, andreasen_huge_calibrate,
    andreasen_huge_strike_greeks,
)
from .cheyette import (
    cheyette_G, cheyette_y, zero_bond as cheyette_zero_bond,
    bond_option as cheyette_bond_option, caplet as cheyette_caplet,
    bond_option_greeks as cheyette_bond_option_greeks,
    floorlet as cheyette_floorlet, cap as cheyette_cap, floor as cheyette_floor,
)
from .g2pp import (
    g2pp_V, zero_bond as g2pp_zero_bond, bond_option as g2pp_bond_option,
    caplet as g2pp_caplet, bond_option_greeks as g2pp_bond_option_greeks,
    floorlet as g2pp_floorlet, cap as g2pp_cap, floor as g2pp_floor,
)
from .bermudan_swaption import (
    bermudan_swaption_g2pp, bermudan_swaption_g2pp_greeks,
)
from .discount_curve import DiscountCurve, bootstrap_from_swaps
from .dualcurve import (
    forward_rate as dual_forward_rate, par_swap_rate as dual_par_swap_rate,
    swap_value as dual_swap_value, float_leg_value as dual_float_leg_value,
    calibrate_basis as dual_calibrate_basis, swap_dv01 as dual_swap_dv01,
)

# Vectorized NumPy fast path is optional; only expose it if NumPy is present.
try:  # pragma: no cover - trivial availability branch
    from .vectorized import (
        price_array, delta_array, gamma_array, vega_array, greeks_array,
        HAS_NUMPY,
    )
except ImportError:  # pragma: no cover
    HAS_NUMPY = False

__version__ = "1.270.0"

__all__ = [
    "OptionType",
    "price",
    "call_price",
    "put_price",
    "delta",
    "gamma",
    "vega",
    "theta",
    "rho",
    "epsilon",
    "greeks",
    "Greeks",
    "implied_volatility",
    "implied_vol_smile",
    "american_price",
    "leisen_reimer_price",
    "leisen_reimer_greeks",
    "leisen_reimer_american_accel",
    "trinomial_price",
    "richardson_american",
    "bermudan_lsm",
    "bermudan_lsm_local_vol",
    "bermudan_lsm_greeks",
    "bermudan_max_call_lsm",
    "bermudan_max_call_lsm_greeks",
    "bermudan_spread_lsm",
    "bermudan_spread_lsm_greeks",
    "bermudan_min_put_lsm",
    "bermudan_min_put_lsm_greeks",
    "bermudan_basket_lsm",
    "bermudan_basket_lsm_greeks",
    "mlmc_asian",
    "perpetual_american",
    "perpetual_exercise_boundary",
    "implied_forward",
    "ForwardResult",
    "dividend_curve",
    "forward_start_price",
    "forward_start_greeks",
    "cliquet_price",
    "cliquet_greeks",
    "chooser_option",
    "chooser_option_greeks",
    "compound_option",
    "compound_option_greeks",
    "quanto_option",
    "compo_option",
    "quanto_option_greeks",
    "compo_option_greeks",
    "displaced_diffusion_price",
    "displaced_diffusion_greeks",
    "displaced_implied_shift",
    "StickyRule",
    "smile_delta",
    "smile_delta_from_smile",
    "skew_slope",
    "floating_strike_lookback",
    "fixed_strike_lookback",
    "lookback_greeks",
    "heston_price",
    "heston_smile",
    "heston_qe_mc",
    "heston_cv_mc",
    "heston_pathwise_delta",
    "heston_mc_greeks",
    "bates_price",
    "bates_greeks",
    "bates_smile",
    "double_heston_price",
    "double_heston_greeks",
    "double_heston_smile",
    "rough_heston_price",
    "rough_heston_greeks",
    "rough_heston_smile",
    "calibrate_double_heston",
    "calibrate_heston",
    "calibrate_lsv_leverage",
    "crank_nicolson_price",
    "crank_nicolson_greeks",
    "crank_nicolson_barrier",
    "crank_nicolson_digital",
    "crank_nicolson_no_touch",
    "adi_two_asset",
    "adi_spread_option",
    "adi_two_asset_cs",
    "adi_two_asset_american",
    "asian_pde_price",
    "kou_price",
    "kou_greeks",
    "kou_smile",
    "cgmy_price",
    "cgmy_greeks",
    "cgmy_smile",
    "nig_price",
    "nig_greeks",
    "nig_smile",
    "meixner_price",
    "meixner_greeks",
    "meixner_smile",
    "calibrate_levy_smile",
    "levy_psi",
    "LevySurface",
    "LevyCalendarViolation",
    "levy_price",
    "carr_madan_strip",
    "carr_madan_smile_strip",
    "cos_price",
    "cos_greeks",
    "rbergomi_price",
    "rbergomi_smile",
    "rbergomi_price_cv",
    "rbergomi_smile_cv",
    "rbergomi_greeks_cv",
    "bachelier_price",
    "bachelier_delta",
    "bachelier_gamma",
    "bachelier_vega",
    "bachelier_implied_vol",
    "bachelier_greeks",
    "variance_gamma_price",
    "variance_gamma_greeks",
    "variance_gamma_smile",
    "simulate_delta_hedge",
    "HedgeResult",
    "merton_jump_price",
    "merton_jump_greeks",
    "merton_smile",
    "risk_neutral_density",
    "risk_neutral_cdf",
    "price_from_density",
    "density_total_mass",
    "variance_swap_strike",
    "volatility_swap_strike",
    "variance_swap_from_smile",
    "corridor_variance_swap_from_smile",
    "gamma_swap_from_smile",
    "forward_variance_swap_from_smile",
    "variance_term_structure",
    "vix_from_chain",
    "vix_from_smile",
    "svix_from_smile",
    "equity_premium_lower_bound",
    "bkm_moments_from_smile",
    "skew_swap_from_smile",
    "realized_variance",
    "variance_risk_premium",
    "lr_greeks",
    "pathwise_delta",
    "lr_digital_delta",
    "asian_pathwise_vega",
    "mixed_gamma",
    "smoothed_digital_delta",
    "barrier_lr_delta",
    "lr_digital_greeks",
    "moment_risk_premia",
    "risk_neutral_density_from_smile",
    "density_grid_from_smile",
    "price_payoff_from_density",
    "risk_neutral_cdf_from_smile",
    "risk_neutral_quantile_from_smile",
    "risk_neutral_var_from_smile",
    "risk_neutral_cvar_from_smile",
    "tail_probability",
    "density_entropy",
    "expected_shortfall",
    "kl_divergence_smiles",
    "wasserstein_smiles",
    "density_var_es",
    "exchange_option",
    "spread_option",
    "spread_option_bs",
    "basket_option",
    "best_of_call",
    "worst_of_call",
    "best_of_call_closed",
    "worst_of_call_closed",
    "best_of_put_closed",
    "worst_of_put_closed",
    "two_asset_digital",
    "two_asset_asset_or_nothing",
    "correlation_option",
    "two_asset_gap_option",
    "two_asset_digital_greeks",
    "exchange_greeks",
    "spread_greeks",
    "spread_option_bs_greeks",
    "basket_greeks",
    "rainbow_greeks",
    "implied_spread_correlation",
    "implied_spread_correlation_bs",
    "payoff_at_expiry",
    "payoff_profile",
    "break_evens",
    "vertical_spread",
    "straddle",
    "strangle",
    "risk_reversal",
    "butterfly",
    "iron_condor",
    "ratio_spread",
    "backspread",
    "strategy_report",
    "dupire_local_vol",
    "local_vol_from_implied",
    "sabr_local_vol",
    "CubicSpline",
    "SmileSpline",
    "CapletPeriod",
    "caplet_price",
    "cap_price",
    "floor_price",
    "caplet_greeks",
    "cap_greeks",
    "floor_greeks",
    "collar_price",
    "caplet_floorlet_parity",
    "annuity",
    "swaption_price",
    "swaption_parity",
    "swaption_greeks",
    "zero_coupon_bond",
    "zero_coupon_yield",
    "bond_option",
    "vasicek_bond_option_greeks",
    "vasicek_bond_greeks",
    "vasicek_coupon_bond_option",
    "vasicek_swaption",
    "vasicek_caplet",
    "vasicek_floorlet",
    "vasicek_cap",
    "vasicek_floor",
    "cir_zero_coupon_bond",
    "cir_zero_coupon_yield",
    "cir_bond_greeks",
    "cir_bond_option",
    "cir_coupon_bond_option",
    "cir_swaption",
    "cir_caplet",
    "cir_floorlet",
    "cir_cap",
    "cir_floor",
    "holee_zero_coupon_bond",
    "holee_zero_coupon_yield",
    "holee_bond_greeks",
    "holee_bond_option",
    "holee_coupon_bond_option",
    "holee_swaption",
    "holee_caplet",
    "holee_floorlet",
    "holee_cap",
    "holee_floor",
    "Overhedge",
    "digital_call_overhedge",
    "digital_put_overhedge",
    "overhedge_payoff",
    "KeyRateDV01",
    "key_rate_dv01",
    "BookSecondOrder",
    "book_second_order",
    "ThetaCarry",
    "theta_carry_report",
    "BumpGreeks",
    "book_bump_greeks",
    "PnLAttribution",
    "attribute_pnl",
    "CarryRoll",
    "carry_roll_pnl",
    "IncomeMetrics",
    "covered_call",
    "cash_secured_put",
    "VegaBuckets",
    "vega_buckets",
    "halton",
    "european_qmc",
    "Sobol",
    "brownian_bridge_path",
    "sobol_european",
    "sobol_european_rqmc",
    "sobol_asian",
    "sobol_asian_rqmc",
    "sobol_lookback_rqmc",
    "sobol_fixed_lookback_rqmc",
    "sobol_barrier_rqmc",
    "sobol_barrier_digital_rqmc",
    "sobol_autocallable_rqmc",
    "sobol_double_knockout_rqmc",
    "sobol_cliquet_rqmc",
    "sobol_average_strike_rqmc",
    "sobol_geometric_asian_rqmc",
    "sobol_parisian_rqmc",
    "sobol_arithmetic_asian_rqmc",
    "implied_correlation",
    "index_vol_from_correlation",
    "dispersion_basket_vol",
    "correlation_term_structure",
    "cev_price",
    "cev_greeks",
    "noncentral_chisq_cdf",
    "delta_hedge_shares",
    "neutralize",
    "vega_neutral_quantity",
    "gamma_neutral_quantity",
    "kelly_fraction_binary",
    "kelly_fraction_continuous",
    "kelly_growth_rate",
    "corrado_su_call",
    "corrado_su_price",
    "calibrate_corrado_su",
    "corrado_su_implied_vol",
    "corrado_su_smile",
    "realized_skewness",
    "realized_excess_kurtosis",
    "Contract",
    "Position",
    "BookRisk",
    "Book",
    "price_book",
    "SVIParams",
    "calibrate_svi",
    "svi_g",
    "svi_butterfly_arbitrage",
    "svi_is_butterfly_free",
    "lee_wing_slopes",
    "lee_bounds_ok",
    "svi_repair_butterfly",
    "svi_local_variance",
    "svi_surface_local_vol",
    "calibrate_svi_from_prices",
    "svi_variance_swap_strike",
    "svi_vix",
    "svi_svix",
    "svi_bkm_moments",
    "svi_density",
    "SSVIParams",
    "ssvi_phi",
    "ssvi_total_variance",
    "calibrate_ssvi",
    "calibrate_ssvi_arbitrage_free",
    "ssvi_butterfly_free",
    "ssvi_calendar_free",
    "ssvi_is_arbitrage_free",
    "ssvi_local_variance",
    "ssvi_local_vol_from_params",
    "ssvi_local_vol_fn",
    "ssvi_reprice_mc",
    "ssvi_variance_swap_strike",
    "ssvi_vix",
    "ssvi_svix",
    "ssvi_density",
    "ssvi_bkm_moments",
    "SABRParams",
    "sabr_vol",
    "calibrate_sabr",
    "calibrate_sabr_lm",
    "sabr_sensitivities",
    "sabr_jacobian",
    "sabr_option_greeks",
    "sabr_variance_swap_strike",
    "sabr_vix",
    "sabr_bkm_moments",
    "sabr_density",
    "sabr_butterfly_arbitrage",
    "sabr_is_arbitrage_free",
    "sabr_repair_butterfly",
    "VannaVolgaSmile",
    "pillar_vols",
    "atm_dns_strike",
    "strike_from_delta",
    "delta_from_strike",
    "rr_bf_to_pillars",
    "VolCube",
    "cms_adjustment_standard",
    "cms_adjustment_greeks",
    "cms_rate_convexity_replication",
    "cms_rate",
    "VolSurface",
    "SurfaceSlice",
    "CalendarViolation",
    "cash_or_nothing",
    "asset_or_nothing",
    "digital_greeks",
    "barrier_option",
    "barrier_greeks",
    "geometric_asian",
    "geometric_asian_greeks",
    "arithmetic_asian",
    "asian_greeks",
    "one_touch",
    "no_touch",
    "gap_option",
    "gap_option_greeks",
    "power_option",
    "power_option_greeks",
    "barrier_rebate",
    "Barrier",
    "MCResult",
    "european_mc",
    "european_cv_mc",
    "european_is_mc",
    "european_is_adaptive_mc",
    "european_stratified_mc",
    "spread_option_lhs_mc",
    "basket_option_lhs_mc",
    "replicated_mc",
    "digital_is_mc",
    "arithmetic_asian_mc",
    "capped_cliquet_mc",
    "barrier_mc",
    "barrier_digital_mc",
    "parisian_barrier_mc",
    "local_vol_mc",
    "average_strike_asian_mc",
    "autocallable_mc",
    "double_knockout_mc",
    "VaRResult",
    "parametric_var",
    "historical_var",
    "montecarlo_var",
    "close_to_close",
    "ewma_vol",
    "parkinson",
    "garman_klass",
    "rogers_satchell",
    "yang_zhang",
    "vol_report",
    "VolReport",
    "vol_cone",
    "VolConePoint",
    "fit_garch",
    "garch_forecast",
    "GarchParams",
    "vanna",
    "vomma",
    "volga",
    "charm",
    "veta",
    "speed",
    "zomma",
    "color",
    "price_array",
    "delta_array",
    "gamma_array",
    "vega_array",
    "greeks_array",
    "HAS_NUMPY",
    "ScenarioGrid",
    "stress_grid",
    "spot_ladder",
    "bjerksund_stensland",
    "bjerksund_stensland_greeks",
    "bjerksund_stensland_1993",
    "early_exercise_premium",
    "kim_american_put",
    "kim_american_call",
    "kim_exercise_boundary",
    "kim_put_greeks",
    "baw_american",
    "andreasen_huge_prices",
    "andreasen_huge_smile",
    "andreasen_huge_calibrate",
    "andreasen_huge_strike_greeks",
    "cheyette_G",
    "cheyette_y",
    "cheyette_zero_bond",
    "cheyette_bond_option",
    "cheyette_bond_option_greeks",
    "cheyette_floorlet",
    "cheyette_cap",
    "cheyette_floor",
    "cheyette_caplet",
    "g2pp_V",
    "g2pp_zero_bond",
    "g2pp_bond_option",
    "g2pp_bond_option_greeks",
    "g2pp_caplet",
    "g2pp_floorlet",
    "g2pp_cap",
    "g2pp_floor",
    "bermudan_swaption_g2pp",
    "bermudan_swaption_g2pp_greeks",
    "DiscountCurve",
    "bootstrap_from_swaps",
    "dual_forward_rate",
    "dual_par_swap_rate",
    "dual_swap_value",
    "dual_float_leg_value",
    "dual_calibrate_basis",
    "dual_swap_dv01",
    "__version__",
]
