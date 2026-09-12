"""Forward rate agreement (FRA) valuation from a discount curve.

An FRA locks in a simple rate ``K`` over a future accrual period ``[T1, T2]``. The
buyer (payer of fixed) receives the difference between the realized simple rate and
``K`` on a notional. Off a discount curve ``P(t)``:

    forward rate  f = (P(T1)/P(T2) - 1) / tau,   tau = T2 - T1,
    value_payer   = notional * tau * (f - K) * P(T2),

which is the settlement discounted to today. The fair FRA rate is ``f`` -- the
value is zero there. Pure standard library.
"""


def fra_forward_rate(discount, t1, t2):
    """Simple forward rate over ``[t1, t2]`` implied by the discount curve."""
    if not (0.0 <= t1 < t2):
        raise ValueError("require 0 <= t1 < t2")
    tau = t2 - t1
    return (discount(t1) / discount(t2) - 1.0) / tau


def fra_value(discount, contract_rate, t1, t2, notional=1.0, payer=True):
    """Value today of a forward rate agreement.

    A payer (long the FRA, paying fixed ``contract_rate``) gains when the forward
    rate exceeds the contract rate: ``notional * tau * (f - K) * P(t2)``. A
    receiver is the negative. Zero at the fair (forward) rate.
    """
    if not (0.0 <= t1 < t2):
        raise ValueError("require 0 <= t1 < t2")
    tau = t2 - t1
    f = fra_forward_rate(discount, t1, t2)
    val = notional * tau * (f - contract_rate) * discount(t2)
    return val if payer else -val
