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
    greeks,
    Greeks,
)
from .implied import implied_volatility
from .binomial import american_price
from .trinomial import trinomial_price, richardson_american
from .forward import implied_forward, ForwardResult
from .forwardstart import forward_start_price, cliquet_price
from .portfolio import Contract, Position, BookRisk, Book, price_book
from .svi import SVIParams, calibrate_svi
from .sabr import SABRParams, sabr_vol, calibrate_sabr
from .surface import VolSurface, SurfaceSlice, CalendarViolation
from .exotics import (
    cash_or_nothing, asset_or_nothing, barrier_option, geometric_asian, Barrier,
)
from .montecarlo import MCResult, european_mc, arithmetic_asian_mc
from .risk import VaRResult, parametric_var, historical_var, montecarlo_var
from .volatility import (
    close_to_close, ewma_vol, parkinson, garman_klass, rogers_satchell,
    yang_zhang, vol_report, VolReport,
)
from .greeks2 import vanna, vomma, volga, charm, veta, speed, zomma, color
from .scenario import ScenarioGrid, stress_grid, spot_ladder
from .american import bjerksund_stensland

# Vectorized NumPy fast path is optional; only expose it if NumPy is present.
try:  # pragma: no cover - trivial availability branch
    from .vectorized import (
        price_array, delta_array, gamma_array, vega_array, greeks_array,
        HAS_NUMPY,
    )
except ImportError:  # pragma: no cover
    HAS_NUMPY = False

__version__ = "1.4.0"

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
    "geometric_asian",
    "Barrier",
    "MCResult",
    "european_mc",
    "arithmetic_asian_mc",
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
    "__version__",
]
