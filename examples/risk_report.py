"""Print a risk report for a small option book.

Run: python examples/risk_report.py
"""

from quantforge import greeks, implied_volatility, OptionType

BOOK = [
    # (label, S, K, t, r, sigma, type, qty)
    ("SPX 4500C", 4500, 4550, 0.25, 0.045, 0.18, OptionType.CALL, 10),
    ("SPX 4500P", 4500, 4450, 0.25, 0.045, 0.20, OptionType.PUT, -5),
    ("AAPL 190C", 190, 200, 0.5, 0.045, 0.30, OptionType.CALL, 100),
]


def main():
    header = f"{'Position':<12}{'Price':>10}{'Delta':>10}{'Gamma':>10}{'Vega':>10}{'Theta':>10}"
    print(header)
    print("-" * len(header))
    net = {"delta": 0.0, "gamma": 0.0, "vega": 0.0, "theta": 0.0}
    for label, S, K, t, r, sigma, ot, qty in BOOK:
        g = greeks(S, K, t, r, sigma, ot)
        net["delta"] += qty * g.delta
        net["gamma"] += qty * g.gamma
        net["vega"] += qty * g.vega / 100.0     # per vol point
        net["theta"] += qty * g.theta / 365.0   # per calendar day
        print(f"{label:<12}{g.price:>10.2f}{g.delta:>10.4f}{g.gamma:>10.5f}"
              f"{g.vega/100:>10.4f}{g.theta/365:>10.4f}")
    print("-" * len(header))
    print(f"{'NET':<12}{'':>10}{net['delta']:>10.2f}{net['gamma']:>10.4f}"
          f"{net['vega']:>10.2f}{net['theta']:>10.2f}")
    print("\n(vega per 1 vol point, theta per calendar day)")

    # Implied vol from an observed price.
    iv = implied_volatility(120.0, S=4500, K=4550, t=0.25, r=0.045,
                            option_type=OptionType.CALL)
    print(f"\nImplied vol of SPX 4550C @ 120.00: {iv:.2%}")


if __name__ == "__main__":
    main()
