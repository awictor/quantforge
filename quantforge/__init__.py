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
from .portfolio import Contract, Position, BookRisk, Book, price_book
from .svi import SVIParams, calibrate_svi
from .exotics import (
    cash_or_nothing, asset_or_nothing, barrier_option, geometric_asian, Barrier,
)
from .montecarlo import MCResult, european_mc, arithmetic_asian_mc
from .risk import VaRResult, parametric_var, historical_var, montecarlo_var
from .volatility import (
    close_to_close, ewma_vol, parkinson, garman_klass, rogers_satchell,
    yang_zhang, vol_report, VolReport,
)

__version__ = "0.7.0"

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
    "Contract",
    "Position",
    "BookRisk",
    "Book",
    "price_book",
    "SVIParams",
    "calibrate_svi",
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
    "__version__",
]
