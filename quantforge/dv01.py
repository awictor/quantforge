"""Key-rate (bucketed) DV01 for a book of rate options.

DV01 is the change in present value for a 1 basis-point parallel shift in the
zero curve. **Key-rate** DV01 decomposes that into the sensitivity to each
curve tenor ("bucket") independently: bump one tenor's zero rate by 1bp, hold
the rest, and reprice. The sum of the key-rate DV01s equals the parallel DV01,
which the test suite checks.

The book is expressed as a callable ``price(zero_rates)`` that reprices given a
mapping ``{tenor: zero_rate}``; this module bumps each tenor and forms the
finite-difference sensitivities. It is model-agnostic - the pricing function
can use the Bachelier caps/swaptions in :mod:`quantforge.rates` or anything
else.
"""

from dataclasses import dataclass, field
from typing import Callable, Dict, List

BP = 1e-4  # one basis point


@dataclass(frozen=True)
class KeyRateDV01:
    buckets: Dict[float, float]     # tenor -> DV01 (PV change per +1bp on that tenor)
    parallel: float                 # DV01 for a simultaneous +1bp on all tenors
    total_bucketed: float           # sum of the per-bucket DV01s

    def as_list(self) -> List[tuple]:
        return sorted(self.buckets.items())


def key_rate_dv01(price_fn: Callable[[Dict[float, float]], float],
                  base_curve: Dict[float, float], bump: float = BP,
                  one_sided: bool = False) -> KeyRateDV01:
    """Compute key-rate DV01s of a book about ``base_curve``.

    Args:
        price_fn: reprices the book given a ``{tenor: zero_rate}`` curve.
        base_curve: the current zero curve.
        bump: the rate shift per bucket (default 1bp). DV01 is reported as the
            PV change for a +1bp move, scaled from the actual bump.
        one_sided: use a forward difference instead of the default central one
            (cheaper, slightly less accurate).

    Returns a :class:`KeyRateDV01`. By convention DV01 is negative for a long
    bond-like position (rates up -> PV down).
    """
    if bump <= 0:
        raise ValueError("bump must be positive")
    scale = BP / bump  # normalize the reported sensitivity to a 1bp move

    buckets = {}
    for tenor in base_curve:
        up = dict(base_curve)
        up[tenor] += bump
        if one_sided:
            base_pv = price_fn(base_curve)
            dv = (price_fn(up) - base_pv)
        else:
            dn = dict(base_curve)
            dn[tenor] -= bump
            dv = 0.5 * (price_fn(up) - price_fn(dn))
        buckets[tenor] = dv * scale

    # Parallel bump: all tenors together.
    up_all = {tt: rr + bump for tt, rr in base_curve.items()}
    if one_sided:
        parallel = (price_fn(up_all) - price_fn(base_curve)) * scale
    else:
        dn_all = {tt: rr - bump for tt, rr in base_curve.items()}
        parallel = 0.5 * (price_fn(up_all) - price_fn(dn_all)) * scale

    return KeyRateDV01(buckets=buckets, parallel=parallel,
                       total_bucketed=sum(buckets.values()))
