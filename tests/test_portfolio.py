"""Tests for batch pricing, portfolio aggregation, and the CLI."""

import pytest

from quantforge import (
    Contract, price_book, greeks, OptionType, call_price,
)
from quantforge.cli import main as cli_main


def test_price_book_matches_individual():
    contracts = [
        Contract(S=100, K=105, t=0.5, r=0.04, sigma=0.25, option_type="call", qty=10),
        Contract(S=100, K=95, t=0.5, r=0.04, sigma=0.30, option_type="put", qty=-5),
    ]
    book = price_book(contracts)
    assert len(book.positions) == 2

    g0 = greeks(100, 105, 0.5, 0.04, 0.25, OptionType.CALL)
    assert book.positions[0].greeks.price == pytest.approx(g0.price)
    assert book.positions[0].position_delta == pytest.approx(10 * g0.delta)


def test_net_greeks_are_scaled_sums():
    contracts = [
        Contract(S=100, K=100, t=1.0, r=0.05, sigma=0.2, option_type="call",
                 qty=2, multiplier=100),
        Contract(S=100, K=100, t=1.0, r=0.05, sigma=0.2, option_type="put",
                 qty=-1, multiplier=100),
    ]
    book = price_book(contracts)
    expected_delta = sum(p.position_delta for p in book.positions)
    expected_mv = sum(p.market_value for p in book.positions)
    assert book.net.delta == pytest.approx(expected_delta)
    assert book.net.market_value == pytest.approx(expected_mv)


def test_multiplier_scales_market_value():
    c = Contract(S=100, K=100, t=1.0, r=0.05, sigma=0.2, option_type="call",
                 qty=3, multiplier=100)
    book = price_book([c])
    unit = call_price(100, 100, 1.0, 0.05, 0.2)
    assert book.net.market_value == pytest.approx(3 * 100 * unit)


def test_short_position_flips_delta_sign():
    long = price_book([Contract(S=100, K=100, t=1, r=0.05, sigma=0.2,
                                option_type="call", qty=1)])
    short = price_book([Contract(S=100, K=100, t=1, r=0.05, sigma=0.2,
                                 option_type="call", qty=-1)])
    assert short.net.delta == pytest.approx(-long.net.delta)


# --- CLI ---
def test_cli_price(capsys):
    rc = cli_main(["price", "-S", "42", "-K", "40", "-t", "0.5",
                   "-r", "0.10", "--sigma", "0.20", "--type", "call"])
    assert rc == 0
    out = capsys.readouterr().out.strip()
    assert float(out) == pytest.approx(4.7594, abs=1e-3)


def test_cli_greeks(capsys):
    rc = cli_main(["greeks", "-S", "100", "-K", "100", "-t", "1",
                   "-r", "0.05", "--sigma", "0.2"])
    assert rc == 0
    out = capsys.readouterr().out
    assert "delta" in out and "gamma" in out


def test_cli_iv_roundtrip(capsys):
    price_val = call_price(100, 105, 0.5, 0.04, 0.25)
    rc = cli_main(["iv", "-S", "100", "-K", "105", "-t", "0.5",
                   "-r", "0.04", "--price", str(price_val), "--type", "call"])
    assert rc == 0
    out = capsys.readouterr().out.strip()
    assert float(out) == pytest.approx(0.25, abs=1e-4)


def test_cli_iv_arbitrage_error(capsys):
    rc = cli_main(["iv", "-S", "100", "-K", "100", "-t", "1",
                   "-r", "0.05", "--price", "500", "--type", "call"])
    assert rc == 2
    err = capsys.readouterr().err
    assert "no-arbitrage" in err


def test_cli_american(capsys):
    rc = cli_main(["american", "-S", "100", "-K", "100", "-t", "1",
                   "-r", "0.05", "--sigma", "0.2", "--type", "put", "--steps", "300"])
    assert rc == 0
    assert float(capsys.readouterr().out.strip()) > 0
