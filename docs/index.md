# QuantForge

Fast, dependency-free options pricing and risk engine in pure Python.

QuantForge prices European, American, and exotic options; fits and interpolates
volatility surfaces; measures portfolio risk; and computes the full set of
Greeks — with zero third-party dependencies in its core.

- **Pricing models** — Black-Scholes-Merton (generalized carry), Bachelier
  (normal), Heston, SABR, and Merton jump-diffusion.
- **American exercise** — Cox-Ross-Rubinstein binomial, Boyle trinomial with
  Richardson extrapolation, and the Bjerksund-Stensland (2002) closed form.
- **Exotics** — digitals, one-touch/no-touch, single barriers, geometric and
  arithmetic Asians, lookbacks, forward-starts, cliquets, and two-asset
  exchange/spread/basket options.
- **Volatility** — SVI and SABR smiles, a calendar-arbitrage-aware term-structure
  surface, cubic-spline interpolation, Dupire local vol, and realized-vol
  estimators.
- **Risk** — batch pricing, net Greeks, VaR / Expected Shortfall, a spot×vol
  stress grid, a delta-hedge P&L simulator, and key-rate DV01.
- **Rates** — Bachelier caps, floors, collars, and swaptions.

See the [API reference](api.md) for every public function and class. Usage
examples live in the project README.

## Install

```bash
pip install quantforge          # core, zero dependencies
pip install "quantforge[fast]"  # optional NumPy vectorized fast path
```

## Building these docs

```bash
python docs/gen_api.py     # regenerate the API reference from the source
mkdocs serve               # preview locally (requires mkdocs)
```
