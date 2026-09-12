"""Par yields and par swap rates from a discount curve.

The par yield of a maturity ``T`` is the coupon rate that makes a bond priced off
the given discount factors worth par (100% of face). For annual (or
``freq``-per-year) coupon dates ``t_1..t_n = T`` with discount factors
``P(t_k)``:

    par_rate = (1 - P(T)) / (sum_k tau_k P(t_k)),

the classic result that the par rate equals the fixed rate of a par swap. Takes a
callable discount curve ``P(t)`` with ``P(0) = 1``. Pure standard library.
"""


def par_yield(discount, maturity, freq=1):
    """Par coupon rate for ``maturity`` from a discount curve ``discount(t)``.

    Parameters
    ----------
    discount : callable
        Discount factor ``P(t)`` (``P(0) = 1``, decreasing for positive rates).
    maturity : float
        Bond maturity in years.
    freq : int
        Coupon payments per year.

    Returns
    -------
    float
        The annualized par coupon rate. On a flat curve it equals the flat rate;
        a bond bearing this coupon prices to exactly par.
    """
    if maturity <= 0:
        raise ValueError("maturity must be positive")
    if freq < 1:
        raise ValueError("freq must be >= 1")
    n = int(round(maturity * freq))
    if n < 1:
        raise ValueError("maturity too short for the coupon frequency")
    tau = 1.0 / freq
    annuity = sum(tau * discount(k * tau) for k in range(1, n + 1))
    if annuity <= 0.0:
        raise ValueError("degenerate annuity")
    return (1.0 - discount(n * tau)) / annuity


def par_bond_price(discount, coupon_rate, maturity, freq=1, face=100.0):
    """Price of a coupon bond off the discount curve (for the par-yield check).

    Sums the discounted coupons plus the discounted principal. Equals ``face``
    exactly when ``coupon_rate`` is the :func:`par_yield`.
    """
    if maturity <= 0 or freq < 1:
        raise ValueError("require maturity > 0 and freq >= 1")
    n = int(round(maturity * freq))
    tau = 1.0 / freq
    c = coupon_rate * tau * face
    pv = sum(c * discount(k * tau) for k in range(1, n + 1))
    pv += face * discount(n * tau)
    return pv
