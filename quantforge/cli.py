"""Command-line interface for QuantForge.

Examples:
    quantforge price -S 100 -K 105 -t 0.5 -r 0.04 --sigma 0.25 --type call
    quantforge greeks -S 100 -K 105 -t 0.5 -r 0.04 --sigma 0.25
    quantforge iv --price 6.12 -S 100 -K 105 -t 0.5 -r 0.04 --type call
    quantforge american -S 100 -K 100 -t 1 -r 0.05 --sigma 0.2 --type put -b 0.02
"""

import argparse
import sys

from . import (
    price, greeks, implied_volatility, american_price, OptionType, __version__,
)


def _common(sub):
    sub.add_argument("-S", "--spot", type=float, required=True, help="spot price")
    sub.add_argument("-K", "--strike", type=float, required=True, help="strike")
    sub.add_argument("-t", "--time", type=float, required=True, help="years to expiry")
    sub.add_argument("-r", "--rate", type=float, required=True, help="risk-free rate")
    sub.add_argument("-b", "--carry", type=float, default=None,
                     help="cost of carry (default = rate)")
    sub.add_argument("--type", choices=["call", "put"], default="call")


def build_parser():
    p = argparse.ArgumentParser(prog="quantforge", description="Options pricing and risk.")
    p.add_argument("--version", action="version", version=f"quantforge {__version__}")
    sub = p.add_subparsers(dest="cmd", required=True)

    pp = sub.add_parser("price", help="price a European option")
    _common(pp)
    pp.add_argument("--sigma", type=float, required=True, help="volatility")

    gp = sub.add_parser("greeks", help="full risk report")
    _common(gp)
    gp.add_argument("--sigma", type=float, required=True, help="volatility")

    ip = sub.add_parser("iv", help="solve implied volatility from a price")
    _common(ip)
    ip.add_argument("--price", type=float, required=True, help="observed option price")

    ap = sub.add_parser("american", help="price an American option (binomial)")
    _common(ap)
    ap.add_argument("--sigma", type=float, required=True, help="volatility")
    ap.add_argument("--steps", type=int, default=500, help="tree steps")
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    ot = OptionType(args.type)

    if args.cmd == "price":
        v = price(args.spot, args.strike, args.time, args.rate, args.sigma, ot, args.carry)
        print(f"{v:.6f}")
    elif args.cmd == "greeks":
        g = greeks(args.spot, args.strike, args.time, args.rate, args.sigma, ot, args.carry)
        print(f"price {g.price:.6f}")
        print(f"delta {g.delta:.6f}")
        print(f"gamma {g.gamma:.6f}")
        print(f"vega  {g.vega:.6f}")
        print(f"theta {g.theta:.6f}")
        print(f"rho   {g.rho:.6f}")
    elif args.cmd == "iv":
        try:
            iv = implied_volatility(args.price, args.spot, args.strike, args.time,
                                    args.rate, ot, args.carry)
        except ValueError as e:
            print(f"error: {e}", file=sys.stderr)
            return 2
        print(f"{iv:.6f}")
    elif args.cmd == "american":
        v = american_price(args.spot, args.strike, args.time, args.rate, args.sigma,
                           ot, args.carry, steps=args.steps)
        print(f"{v:.6f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
