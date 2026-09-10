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
from .implied import implied_volatility
from .binomial import american_price
from .trinomial import trinomial_price, richardson_american
from .forward import implied_forward, ForwardResult
from .forwardstart import forward_start_price, cliquet_price
from .hedging import StickyRule, smile_delta, smile_delta_from_smile, skew_slope
from .lookback import floating_strike_lookback, fixed_strike_lookback
from .heston import heston_price
from .bachelier import (
    bachelier_price, bachelier_delta, bachelier_gamma, bachelier_vega,
    bachelier_implied_vol,
)
from .hedgesim import simulate_delta_hedge, HedgeResult
from .merton import merton_jump_price
from .density import (
    risk_neutral_density, risk_neutral_cdf, price_from_density, density_total_mass,
)
from .varswap import variance_swap_strike, volatility_swap_strike
from .multiasset import exchange_option, spread_option, basket_option
from .strategy import (
    payoff_at_expiry, payoff_profile, break_evens,
    vertical_spread, straddle, strangle, risk_reversal, butterfly, iron_condor,
)
from .localvol import dupire_local_vol, local_vol_from_implied
from .spline import CubicSpline, SmileSpline
from .rates import (
    CapletPeriod, caplet_price, cap_price, floor_price, collar_price,
    caplet_floorlet_parity, annuity, swaption_price, swaption_parity,
)
from .overhedge import (
    Overhedge, digital_call_overhedge, digital_put_overhedge, overhedge_payoff,
)
from .dv01 import KeyRateDV01, key_rate_dv01
from .bookgreeks import BookSecondOrder, book_second_order
from .qmc import halton, european_qmc
from .correlation import (
    implied_correlation, index_vol_from_correlation, dispersion_basket_vol,
)
from .cev import cev_price, noncentral_chisq_cdf
from .sizing import (
    delta_hedge_shares, neutralize, vega_neutral_quantity, gamma_neutral_quantity,
)
from .portfolio import Contract, Position, BookRisk, Book, price_book
from .svi import SVIParams, calibrate_svi
from .sabr import SABRParams, sabr_vol, calibrate_sabr
from .surface import VolSurface, SurfaceSlice, CalendarViolation
from .exotics import (
    cash_or_nothing, asset_or_nothing, barrier_option, barrier_greeks,
    geometric_asian, arithmetic_asian, one_touch, no_touch, Barrier,
)
from .montecarlo import (
    MCResult, european_mc, arithmetic_asian_mc, capped_cliquet_mc,
)
from .risk import VaRResult, parametric_var, historical_var, montecarlo_var
from .volatility import (
    close_to_close, ewma_vol, parkinson, garman_klass, rogers_satchell,
    yang_zhang, vol_report, VolReport,
)
from .greeks2 import vanna, vomma, volga, charm, veta, speed, zomma, color
from .scenario import ScenarioGrid, stress_grid, spot_ladder
from .american import bjerksund_stensland, bjerksund_stensland_greeks

# Vectorized NumPy fast path is optional; only expose it if NumPy is present.
try:  # pragma: no cover - trivial availability branch
    from .vectorized import (
        price_array, delta_array, gamma_array, vega_array, greeks_array,
        HAS_NUMPY,
    )
except ImportError:  # pragma: no cover
    HAS_NUMPY = False

__version__ = "1.33.0"

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
    "american_price",
    "trinomial_price",
    "richardson_american",
    "implied_forward",
    "ForwardResult",
    "forward_start_price",
    "cliquet_price",
    "StickyRule",
    "smile_delta",
    "smile_delta_from_smile",
    "skew_slope",
    "floating_strike_lookback",
    "fixed_strike_lookback",
    "heston_price",
    "bachelier_price",
    "bachelier_delta",
    "bachelier_gamma",
    "bachelier_vega",
    "bachelier_implied_vol",
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
    "payoff_at_expiry",
    "payoff_profile",
    "break_evens",
    "vertical_spread",
    "straddle",
    "strangle",
    "risk_reversal",
    "butterfly",
    "iron_condor",
    "dupire_local_vol",
    "local_vol_from_implied",
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
    "Overhedge",
    "digital_call_overhedge",
    "digital_put_overhedge",
    "overhedge_payoff",
    "KeyRateDV01",
    "key_rate_dv01",
    "BookSecondOrder",
    "book_second_order",
    "halton",
    "european_qmc",
    "implied_correlation",
    "index_vol_from_correlation",
    "dispersion_basket_vol",
    "cev_price",
    "noncentral_chisq_cdf",
    "delta_hedge_shares",
    "neutralize",
    "vega_neutral_quantity",
    "gamma_neutral_quantity",
    "Contract",
    "Position",
    "BookRisk",
    "Book",
    "price_book",
    "SVIParams",
    "calibrate_svi",
    "SABRParams",
    "sabr_vol",
    "calibrate_sabr",
    "VolSurface",
    "SurfaceSlice",
    "CalendarViolation",
    "cash_or_nothing",
    "asset_or_nothing",
    "barrier_option",
    "barrier_greeks",
    "geometric_asian",
    "arithmetic_asian",
    "one_touch",
    "no_touch",
    "Barrier",
    "MCResult",
    "european_mc",
    "arithmetic_asian_mc",
    "capped_cliquet_mc",
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
    "__version__",
]
