"""End-to-end volatility-surface workflow on a synthetic option chain.

Steps, all with QuantForge and no third-party dependencies:

  1. Take a multi-expiry chain of option mid-prices.
  2. Invert each quote to a Black-Scholes implied volatility.
  3. Fit a raw-SVI slice per expiry (total variance vs log-moneyness).
  4. Assemble the slices into a term-structure surface and check it for
     calendar arbitrage.
  5. Read an interpolated vol at an arbitrary (strike, expiry), and extract the
     Dupire local vol there.

Run: python examples/vol_surface.py
"""

import math

from quantforge import (
    implied_volatility, calibrate_svi, VolSurface, SurfaceSlice,
    local_vol_from_implied, OptionType, call_price,
)


# A synthetic chain: forward ~ spot (r=0), a downward equity skew per expiry.
SPOT = 100.0
RATE = 0.0
EXPIRIES = [0.25, 0.5, 1.0]
STRIKES = [80, 90, 95, 100, 105, 110, 120]


def _synthetic_vol(K, T):
    """A plausible skewed smile: higher vol for low strikes, rising with T."""
    m = math.log(K / SPOT)
    return 0.20 + 0.10 * T - 0.5 * m + 0.8 * m * m


def build_chain():
    """Return {T: [(K, call_mid), ...]} priced at the synthetic smile."""
    chain = {}
    for T in EXPIRIES:
        rows = []
        for K in STRIKES:
            sigma = _synthetic_vol(K, T)
            rows.append((K, call_price(SPOT, K, T, RATE, sigma)))
        chain[T] = rows
    return chain


def main():
    chain = build_chain()

    # 2-3. Invert to implied vol, then fit an SVI slice per expiry.
    slices = []
    for T, rows in chain.items():
        ks, tv = [], []
        for K, price in rows:
            iv = implied_volatility(price, SPOT, K, T, RATE, OptionType.CALL)
            ks.append(math.log(K / SPOT))
            tv.append(iv * iv * T)          # total variance
        params, rmse = calibrate_svi(ks, tv)
        slices.append(SurfaceSlice(t=T, params=params, rmse=rmse))
        print(f"T={T:>4}: SVI rmse={rmse:.2e}  a={params.a:.4f} b={params.b:.4f} "
              f"rho={params.rho:+.3f} m={params.m:+.3f} s={params.s:.4f}")

    # 4. Assemble the surface and check calendar arbitrage.
    surface = VolSurface(slices)
    violations = surface.calendar_arbitrage()
    print(f"\nCalendar-arbitrage violations: {len(violations)}")

    # 5. Interpolated implied vol and Dupire local vol at a sample point.
    k, T = math.log(105 / SPOT), 0.75
    iv = surface.implied_vol(k, T)
    print(f"\nInterpolated implied vol at K=105, T=0.75: {iv:.4f}")

    def smile(K, t):
        return surface.implied_vol(math.log(K / SPOT), t)

    lv = local_vol_from_implied(smile, SPOT, 105, 0.75, RATE)
    print(f"Dupire local vol at K=105, T=0.75:          {lv:.4f}")


if __name__ == "__main__":
    main()
