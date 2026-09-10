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
from .trinomial import trinomial_price, richardson_american
from .lsm import bermudan_lsm
from .perpetual import perpetual_american, perpetual_exercise_boundary
from .forward import implied_forward, ForwardResult, dividend_curve
from .forwardstart import forward_start_price, cliquet_price
from .chooser import chooser_option
from .compound import compound_option
from .quanto import quanto_option, compo_option
from .displaced import displaced_diffusion_price
from .hedging import StickyRule, smile_delta, smile_delta_from_smile, skew_slope
from .lookback import (
    floating_strike_lookback, fixed_strike_lookback, lookback_greeks,
)
from .heston import heston_price
from .bachelier import (
    bachelier_price, bachelier_delta, bachelier_gamma, bachelier_vega,
    bachelier_implied_vol,
)
from .variancegamma import variance_gamma_price
from .hedgesim import simulate_delta_hedge, HedgeResult
from .merton import merton_jump_price
from .density import (
    risk_neutral_density, risk_neutral_cdf, price_from_density, density_total_mass,
)
from .varswap import variance_swap_strike, volatility_swap_strike
from .multiasset import (
    exchange_option, spread_option, basket_option, best_of_call, worst_of_call,
    exchange_greeks, spread_greeks, basket_greeks, implied_spread_correlation,
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
)
from .vasicek import (
    zero_coupon_bond, zero_coupon_yield, bond_option,
)
from .cir import cir_zero_coupon_bond, cir_zero_coupon_yield
from .holee import holee_zero_coupon_bond, holee_zero_coupon_yield
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
from .correlation import (
    implied_correlation, index_vol_from_correlation, dispersion_basket_vol,
    correlation_term_structure,
)
from .cev import cev_price, noncentral_chisq_cdf
from .sizing import (
    delta_hedge_shares, neutralize, vega_neutral_quantity, gamma_neutral_quantity,
    kelly_fraction_binary, kelly_fraction_continuous, kelly_growth_rate,
)
from .gramcharlier import (
    corrado_su_call, corrado_su_price, realized_skewness, realized_excess_kurtosis,
)
from .portfolio import Contract, Position, BookRisk, Book, price_book
from .svi import (
    SVIParams, calibrate_svi, svi_g, svi_butterfly_arbitrage, svi_is_butterfly_free,
    lee_wing_slopes, lee_bounds_ok, svi_repair_butterfly,
)
from .sabr import SABRParams, sabr_vol, calibrate_sabr
from .vannavolga import VannaVolgaSmile, pillar_vols
from .surface import VolSurface, SurfaceSlice, CalendarViolation
from .exotics import (
    cash_or_nothing, asset_or_nothing, digital_greeks, barrier_option, barrier_greeks,
    geometric_asian, arithmetic_asian, asian_greeks, one_touch, no_touch,
    gap_option, power_option, barrier_rebate, Barrier,
)
from .montecarlo import (
    MCResult, european_mc, arithmetic_asian_mc, capped_cliquet_mc,
    barrier_digital_mc, parisian_barrier_mc, local_vol_mc,
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

# Vectorized NumPy fast path is optional; only expose it if NumPy is present.
try:  # pragma: no cover - trivial availability branch
    from .vectorized import (
        price_array, delta_array, gamma_array, vega_array, greeks_array,
        HAS_NUMPY,
    )
except ImportError:  # pragma: no cover
    HAS_NUMPY = False

__version__ = "1.85.0"

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
    "trinomial_price",
    "richardson_american",
    "bermudan_lsm",
    "perpetual_american",
    "perpetual_exercise_boundary",
    "implied_forward",
    "ForwardResult",
    "dividend_curve",
    "forward_start_price",
    "cliquet_price",
    "chooser_option",
    "compound_option",
    "quanto_option",
    "compo_option",
    "displaced_diffusion_price",
    "StickyRule",
    "smile_delta",
    "smile_delta_from_smile",
    "skew_slope",
    "floating_strike_lookback",
    "fixed_strike_lookback",
    "lookback_greeks",
    "heston_price",
    "bachelier_price",
    "bachelier_delta",
    "bachelier_gamma",
    "bachelier_vega",
    "bachelier_implied_vol",
    "variance_gamma_price",
    "simulate_delta_hedge",
    "HedgeResult",
    "merton_jump_price",
    "risk_neutral_density",
    "risk_neutral_cdf",
    "price_from_density",
    "density_total_mass",
    "variance_swap_strike",
    "volatility_swap_strike",
    "exchange_option",
    "spread_option",
    "basket_option",
    "best_of_call",
    "worst_of_call",
    "exchange_greeks",
    "spread_greeks",
    "basket_greeks",
    "implied_spread_correlation",
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
    "collar_price",
    "caplet_floorlet_parity",
    "annuity",
    "swaption_price",
    "swaption_parity",
    "zero_coupon_bond",
    "zero_coupon_yield",
    "bond_option",
    "cir_zero_coupon_bond",
    "cir_zero_coupon_yield",
    "holee_zero_coupon_bond",
    "holee_zero_coupon_yield",
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
    "implied_correlation",
    "index_vol_from_correlation",
    "dispersion_basket_vol",
    "correlation_term_structure",
    "cev_price",
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
    "SABRParams",
    "sabr_vol",
    "calibrate_sabr",
    "VannaVolgaSmile",
    "pillar_vols",
    "VolSurface",
    "SurfaceSlice",
    "CalendarViolation",
    "cash_or_nothing",
    "asset_or_nothing",
    "digital_greeks",
    "barrier_option",
    "barrier_greeks",
    "geometric_asian",
    "arithmetic_asian",
    "asian_greeks",
    "one_touch",
    "no_touch",
    "gap_option",
    "power_option",
    "barrier_rebate",
    "Barrier",
    "MCResult",
    "european_mc",
    "arithmetic_asian_mc",
    "capped_cliquet_mc",
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
    "__version__",
]
