"""Strongly connected components of a directed graph (Tarjan's algorithm).

Two nodes are in the same strongly connected component (SCC) when each is reachable from
the other. Tarjan's algorithm finds every SCC in a single ``O(V + E)`` depth-first pass,
tracking each node's DFS discovery index and the lowest index reachable from its subtree;
a node whose low-link equals its own index roots an SCC. Collapsing each SCC to a single
vertex yields the *condensation*, always a DAG. Implemented iteratively so large graphs do
not overflow recursion. Pure standard library.
"""


def strongly_connected_components(graph):
    """Return the SCCs of a directed ``graph`` as a list of node lists.

    ``graph`` is ``{node: [successors]}``; every node must appear as a key. Components are
    returned in reverse topological order of the condensation (a component appears before
    the components it can reach), which is the natural output order of Tarjan's algorithm.
    """
    index_of = {}
    low = {}
    on_stack = {}
    stack = []
    components = []
    counter = [0]

    def strongconnect(root):
        # explicit work stack of (node, iterator-position) frames
        work = [(root, 0)]
        while work:
            node, pi = work[-1]
            if pi == 0:
                index_of[node] = counter[0]
                low[node] = counter[0]
                counter[0] += 1
                stack.append(node)
                on_stack[node] = True
            recursed = False
            succ = graph[node]
            i = pi
            while i < len(succ):
                w = succ[i]
                if w not in index_of:
                    work[-1] = (node, i + 1)
                    work.append((w, 0))
                    recursed = True
                    break
                elif on_stack.get(w):
                    low[node] = min(low[node], index_of[w])
                i += 1
            if recursed:
                continue
            # done with node's successors
            if low[node] == index_of[node]:
                comp = []
                while True:
                    w = stack.pop()
                    on_stack[w] = False
                    comp.append(w)
                    if w == node:
                        break
                components.append(comp)
            work.pop()
            if work:
                parent = work[-1][0]
                low[parent] = min(low[parent], low[node])

    for node in graph:
        if node not in index_of:
            strongconnect(node)
    return components


def condensation(graph):
    """Collapse each SCC to a node, returning ``(component_of, dag)``.

    ``component_of`` maps each original node to its component id (``0..k-1``); ``dag`` is the
    condensed graph ``{comp_id: [successor_comp_ids]}`` with no duplicate edges and no
    self-loops. The condensation is always acyclic.
    """
    comps = strongly_connected_components(graph)
    component_of = {}
    for cid, comp in enumerate(comps):
        for node in comp:
            component_of[node] = cid
    dag = {cid: [] for cid in range(len(comps))}
    seen = {cid: set() for cid in range(len(comps))}
    for node, succ in graph.items():
        cu = component_of[node]
        for w in succ:
            cw = component_of[w]
            if cw != cu and cw not in seen[cu]:
                seen[cu].add(cw)
                dag[cu].append(cw)
    return component_of, dag


def is_strongly_connected(graph):
    """True if the whole graph is a single strongly connected component."""
    if not graph:
        return True
    return len(strongly_connected_components(graph)) == 1


def number_of_sccs(graph):
    """Count the strongly connected components of ``graph``."""
    return len(strongly_connected_components(graph))
