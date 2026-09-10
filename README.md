# QuantForge

Fast, **dependency-free** options pricing and risk engine in pure Python.

QuantForge prices European and American options, computes the full set of
analytic Greeks, and solves for implied volatility — with zero third-party
dependencies. It drops into any Python 3.8+ runtime (serverless, embedded,
notebooks, trading bots) without compiling NumPy or SciPy.

## Why

Most option libraries pull in a heavy scientific stack. QuantForge implements
the math from scratch — a full-precision normal CDF/quantile, closed-form
Black-Scholes-Merton, analytic Greeks, a robust Newton-with-bisection implied
vol solver, and a Cox-Ross-Rubinstein binomial lattice for early exercise. The
whole thing is a few hundred lines of readable, tested code.

- **Zero dependencies** — standard library only.
- **Generalized BSM** — one model covers stocks, dividends, futures (Black-76)
  and FX (Garman-Kohlhagen) via the cost-of-carry parameter.
- **Every Greek analytic** — delta, gamma, vega, theta, rho, verified against
  finite differences to 1e-4 or better.
- **Robust implied vol** — Newton-Raphson with a guaranteed bisection bracket,
  arbitrage-band rejection.
- **American exercise** — binomial tree that converges to BSM on European
  payoffs (used as a self-check in the suite).

## Install

```bash
pip install quantforge          # from PyPI (planned)
pip install -e .                # from source
```

## Quick start

```python
from quantforge import call_price, greeks, implied_volatility, american_price

# Price a European call: spot 100, strike 105, 6 months, 4% rate, 25% vol.
print(call_price(S=100, K=105, t=0.5, r=0.04, sigma=0.25))

# Full risk report in one call.
g = greeks(S=100, K=105, t=0.5, r=0.04, sigma=0.25)
print(g.delta, g.gamma, g.vega, g.theta, g.rho)

# Back out implied vol from a market price.
iv = implied_volatility(target_price=6.12, S=100, K=105, t=0.5, r=0.04)
print(iv)

# American put with a 3% dividend yield (carry b = r - q).
print(american_price(S=100, K=100, t=1.0, r=0.05, sigma=0.2,
                     option_type="put", b=0.05 - 0.03, steps=500))
```

## Batch pricing and portfolio risk

Value a whole book in one call and get net exposures:

```python
from quantforge import Contract, price_book

book = price_book([
    Contract(S=100, K=105, t=0.5, r=0.04, sigma=0.25, option_type="call",
             qty=10, multiplier=100, label="AAPL 105C"),
    Contract(S=100, K=95,  t=0.5, r=0.04, sigma=0.30, option_type="put",
             qty=-5, multiplier=100, label="AAPL 95P"),
])
print(book.net.delta, book.net.gamma, book.net.vega, book.net.market_value)
for pos in book.positions:
    print(pos.contract.label, pos.position_delta)
```

`qty` is signed (short = negative) and `multiplier` scales to notional
(e.g. 100 for US equity options). Net Greeks are position-scaled sums.

## Command line

Installing the package exposes a `quantforge` CLI:

```bash
quantforge price   -S 100 -K 105 -t 0.5 -r 0.04 --sigma 0.25 --type call
quantforge greeks  -S 100 -K 105 -t 0.5 -r 0.04 --sigma 0.25
quantforge iv      -S 100 -K 105 -t 0.5 -r 0.04 --price 6.12 --type call
quantforge american -S 100 -K 100 -t 1 -r 0.05 --sigma 0.2 --type put -b 0.02
```

## Model coverage

| Instrument            | Set the carry `b` to | Function            |
|-----------------------|----------------------|---------------------|
| Non-dividend stock    | `b = r` (default)    | `price`, `greeks`   |
| Continuous dividend q | `b = r - q`          | `price`, `greeks`   |
| Option on a future    | `b = 0` (Black-76)   | `price`, `greeks`   |
| FX option             | `b = r - r_foreign`  | `price`, `greeks`   |
| American exercise     | any of the above     | `american_price`    |

## Testing

```bash
pip install -e ".[test]"
pytest -q
```

The suite pins prices against Hull's textbook reference values, enforces
put-call parity, checks every Greek against central finite differences, and
round-trips implied volatility.

## License

MIT.
