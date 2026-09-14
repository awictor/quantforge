"""Tests for strongly connected components (Tarjan), cross-checked against brute reachability."""

import random

from quantforge.scc import (
    strongly_connected_components,
    condensation,
    is_strongly_connected,
    number_of_sccs,
)


def _reach(graph, s):
    seen = {s}
    st = [s]
    while st:
        u = st.pop()
        for w in graph[u]:
            if w not in seen:
                seen.add(w)
                st.append(w)
    return seen


def _brute_sccs(graph):
    nodes = list(graph)
    rg = {u: [] for u in nodes}
    for u in nodes:
        for w in graph[u]:
            rg[w].append(u)
    seen = set()
    comps = set()
    for s in nodes:
        if s in seen:
            continue
        comp = _reach(graph, s) & _reach(rg, s)
        comps.add(frozenset(comp))
        seen |= comp
    return comps


def _dag_is_acyclic(dag):
    indeg = {u: 0 for u in dag}
    for u in dag:
        for w in dag[u]:
            indeg[w] += 1
    q = [u for u in dag if indeg[u] == 0]
    cnt = 0
    while q:
        u = q.pop()
        cnt += 1
        for w in dag[u]:
            indeg[w] -= 1
            if indeg[w] == 0:
                q.append(w)
    return cnt == len(dag)


def test_fuzz_vs_brute_reachability():
    rng = random.Random(61)
    for _ in range(3000):
        n = rng.randint(1, 15)
        nodes = list(range(n))
        graph = {u: [] for u in nodes}
        for u in nodes:
            for v in nodes:
                if u != v and rng.random() < 0.2:
                    graph[u].append(v)
        comps = strongly_connected_components(graph)
        flat = [x for c in comps for x in c]
        assert sorted(flat) == nodes  # exact partition
        assert {frozenset(c) for c in comps} == _brute_sccs(graph)


def test_fuzz_condensation_is_acyclic():
    rng = random.Random(62)
    for _ in range(2000):
        n = rng.randint(1, 15)
        nodes = list(range(n))
        graph = {u: [] for u in nodes}
        for u in nodes:
            for v in nodes:
                if u != v and rng.random() < 0.25:
                    graph[u].append(v)
        _, dag = condensation(graph)
        assert _dag_is_acyclic(dag)


def test_simple_cycle_and_tail():
    g = {"a": ["b"], "b": ["c"], "c": ["a"], "d": ["a"]}
    scc_sets = {frozenset(c) for c in strongly_connected_components(g)}
    assert frozenset({"a", "b", "c"}) in scc_sets
    assert frozenset({"d"}) in scc_sets
    assert not is_strongly_connected(g)


def test_fully_strongly_connected():
    assert is_strongly_connected({"a": ["b"], "b": ["a"]})
    assert is_strongly_connected({"a": ["b"], "b": ["c"], "c": ["a"]})


def test_dag_has_singleton_components():
    g = {"a": ["b", "c"], "b": ["d"], "c": ["d"], "d": []}
    comps = strongly_connected_components(g)
    assert all(len(c) == 1 for c in comps)
    assert number_of_sccs(g) == 4


def test_reverse_topological_order():
    # a->b->c->a is one SCC; d->a. d's component reaches abc's, so must sort after it.
    g = {"a": ["b"], "b": ["c"], "c": ["a"], "d": ["a"]}
    comp_of, dag = condensation(g)
    # edge goes from d's component to abc's component
    assert comp_of["a"] in dag[comp_of["d"]]
    assert dag[comp_of["a"]] == []


def test_condensation_no_duplicate_edges():
    # two parallel routes between the same components -> single condensation edge
    g = {"a": ["b", "c"], "b": ["d"], "c": ["d"], "d": []}
    _, dag = condensation(g)
    for cid, succ in dag.items():
        assert len(succ) == len(set(succ))


def test_single_node():
    assert strongly_connected_components({"x": []}) == [["x"]]
    assert number_of_sccs({"x": []}) == 1


def test_self_loop():
    assert strongly_connected_components({"x": ["x"]}) == [["x"]]
    assert is_strongly_connected({"x": ["x"]})


def test_empty_graph():
    assert strongly_connected_components({}) == []
    assert is_strongly_connected({}) is True
    assert number_of_sccs({}) == 0


def test_disconnected_singletons():
    g = {"a": [], "b": [], "c": []}
    assert number_of_sccs(g) == 3
    assert not is_strongly_connected(g)


def test_two_separate_cycles():
    g = {"a": ["b"], "b": ["a"], "c": ["d"], "d": ["c"]}
    scc_sets = {frozenset(x) for x in strongly_connected_components(g)}
    assert scc_sets == {frozenset({"a", "b"}), frozenset({"c", "d"})}


def test_deep_chain_no_recursion_error():
    # long chain would overflow recursive DFS; iterative build must handle it
    n = 5000
    g = {i: [i + 1] for i in range(n - 1)}
    g[n - 1] = []
    assert number_of_sccs(g) == n  # every node its own SCC
