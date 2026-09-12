"""Single-curve vanilla interest-rate swap valuation.

Values a fixed-for-floating swap off one discount curve ``P(t)`` (the classic
single-curve / OIS-discounting-only setup). The floating leg's PV telescopes to
``P(t_0) - P(t_n)`` (par float), so the payer-swap value is

    V_payer = [P(t_start) - P(t_end)] - K * A,   A = sum_k tau_k P(t_k),

where ``A`` is the fixed-leg annuity (PV01 per unit rate). The par swap rate is
``(P(t_start) - P(t_end)) / A`` -- zero value at that fixed rate. Pure standard
library.
"""


def _schedule(start, maturity, freq):
    n = int(round((maturity - start) * freq))
    if n < 1:
        raise ValueError("schedule has no payment periods")
    tau = 1.0 / freq
    return [start + k * tau for k in range(n + 1)], tau


def swap_annuity(discount, start, maturity, freq=2):
    """Fixed-leg annuity ``sum_k tau_k P(t_k)`` (PV01 per unit rate).

    The present value of receiving 1 unit of rate on the fixed schedule.
    """
    times, tau = _schedule(start, maturity, freq)
    return sum(tau * discount(times[k]) for k in range(1, len(times)))


def par_swap_rate(discount, start, maturity, freq=2):
    """Par (fair) fixed rate that gives the swap zero value at inception."""
    times, _ = _schedule(start, maturity, freq)
    a = swap_annuity(discount, start, maturity, freq)
    if a <= 0.0:
        raise ValueError("degenerate annuity")
    return (discount(times[0]) - discount(times[-1])) / a


def swap_value(discount, fixed_rate, start, maturity, freq=2, notional=1.0,
               payer=True):
    """Value of a fixed-for-floating swap off a single discount curve.

    A payer pays ``fixed_rate`` and receives float; its value is the float-leg PV
    (``P(start) - P(end)``) minus the fixed-leg PV (``fixed_rate * annuity``),
    times ``notional``. A receiver is the negative. Zero at the par swap rate.
    """
    times, _ = _schedule(start, maturity, freq)
    float_pv = discount(times[0]) - discount(times[-1])
    fixed_pv = fixed_rate * swap_annuity(discount, start, maturity, freq)
    payer_val = notional * (float_pv - fixed_pv)
    return payer_val if payer else -payer_val
