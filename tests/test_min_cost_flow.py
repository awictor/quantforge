"""Tests for min-cost max-flow, cross-checked against brute enumeration and Edmonds-Karp."""

import itertools
import random

import pytest

from quantforge.min_cost_flow import min_cost_max_flow, MinCostMaxFlow
from quantforge import max_flow as ek_max_flow


def _brute(edges, s, t, nodes):
    ranges = [range(c + 1) for _, _, c, _ in edges]
    best_flow = -1
    best_cost = None
    for combo in itertools.product(*ranges):
        net = {n: 0 for n in nodes}
        for (u, v, cap, cost), f in zip(edges, combo):
            net[u] -= f
            net[v] += f
        if any(net[n] != 0 for n in nodes if n != s and n != t):
            continue
        if net[s] > 0:
            continue
        flow = net[t]
        if flow < 0:
            continue
        cost = sum(f * cst for (u, v, cap, cst), f in zip(edges, combo))
        if flow > best_flow or (flow == best_flow and cost < best_cost):
            best_flow = flow
            best_cost = cost
    return best_flow, best_cost


@pytest.mark.slow
def test_fuzz_vs_brute():
    rng = random.Random(411)
    for _ in range(1500):
        n = rng.randint(2, 5)
        nodes = list(range(n))
        edges = []
        for u in range(n):
            for v in range(n):
                if u != v and rng.random() < 0.5:
                    edges.append((u, v, rng.randint(0, 3), rng.randint(0, 5)))
        if not edges:
            continue
        s, t = 0, n - 1
        assert min_cost_max_flow(edges, s, t) == _brute(edges, s, t, nodes)


@pytest.mark.slow
def test_fuzz_maxflow_agrees_with_edmonds_karp():
    rng = random.Random(412)
    for _ in range(1000):
        n = rng.randint(2, 6)
        edges = []
        g = {i: {} for i in range(n)}
        for u in range(n):
            for v in range(n):
                if u != v and rng.random() < 0.4:
                    c = rng.randint(0, 8)
                    if c > 0:
                        edges.append((u, v, c, rng.randint(1, 5)))
                        g[u][v] = g[u].get(v, 0) + c
        if not edges:
            continue
        s, t = 0, n - 1
        mf, _ = min_cost_max_flow(edges, s, t)
        assert mf == ek_max_flow(g, s, t)


def test_prefers_cheaper_path():
    # two unit paths s->t: one cost 1, one cost 10; min-cost uses the cheap one first
    edges = [(0, 1, 1, 1), (1, 3, 1, 0), (0, 2, 1, 10), (2, 3, 1, 0)]
    flow, cost = min_cost_max_flow(edges, 0, 3)
    assert flow == 2
    assert cost == 11  # 1 (cheap path) + 10 (forced expensive path)


def test_clrs_example():
    edges = [(0, 1, 2, 1), (0, 2, 2, 3), (1, 3, 2, 1), (2, 3, 2, 1), (1, 2, 1, 1)]
    assert min_cost_max_flow(edges, 0, 3) == (4, 12)


def test_negative_cost_edge():
    # a negative-cost edge should be saturated to reduce total cost
    edges = [(0, 1, 1, 2), (1, 2, 1, -5), (0, 2, 1, 1)]
    flow, cost = min_cost_max_flow(edges, 0, 2)
    assert flow == 2
    assert cost == (2 - 5) + 1  # path 0->1->2 costs -3, path 0->2 costs 1


def test_zero_capacity_edge():
    mf, c = min_cost_max_flow([(0, 1, 0, 5)], 0, 1)
    assert mf == 0
    assert c == 0


def test_disconnected_sink():
    assert min_cost_max_flow([(0, 1, 1, 1)], 0, 9) == (0, 0)


def test_single_edge():
    assert min_cost_max_flow([(0, 1, 5, 3)], 0, 1) == (5, 15)


def test_builder_interface():
    m = MinCostMaxFlow()
    m.add_edge("s", "a", 3, 1).add_edge("a", "t", 3, 2)
    assert m.solve("s", "t") == (3, 9)


def test_negative_capacity_raises():
    with pytest.raises(ValueError):
        MinCostMaxFlow().add_edge(0, 1, -1, 2)


def test_source_equals_sink_raises():
    with pytest.raises(ValueError):
        min_cost_max_flow([(0, 1, 1, 1)], 0, 0)
