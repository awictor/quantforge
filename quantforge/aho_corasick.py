"""Aho-Corasick multi-pattern string matching.

Finding every occurrence of *many* patterns in a text at once. A single pass over the text
locates all matches in ``O(len(text) + total_pattern_length + matches)`` time, versus
running a single-pattern search once per pattern. The automaton is a trie of the patterns
augmented with *failure links* (the longest proper suffix that is also a trie prefix) and
*output links* (so a match of one pattern also reports every shorter pattern ending at the
same position). Pure standard library.
"""

from collections import deque


class AhoCorasick:
    """Aho-Corasick automaton over a set of patterns.

    Build once with the patterns, then call :meth:`find` on any number of texts. Each match
    is reported as ``(end_index, pattern)`` where ``end_index`` is the index of the match's
    last character in the text (so the match is ``text[end_index - len(pattern) + 1 :
    end_index + 1]``). Duplicate patterns are collapsed; empty patterns are ignored.
    """

    __slots__ = ("_goto", "_fail", "_output", "_patterns", "_built")

    def __init__(self, patterns=None):
        # node 0 is the root. _goto[node] maps a char -> child node.
        self._goto = [{}]
        self._fail = [0]
        self._output = [[]]        # list of pattern strings ending at this node
        self._patterns = []
        self._built = False
        if patterns:
            for p in patterns:
                self.add(p)
            self.build()

    def add(self, pattern):
        """Add a pattern to the trie (before :meth:`build`)."""
        if self._built:
            raise RuntimeError("cannot add patterns after build()")
        if pattern == "":
            return self
        node = 0
        for ch in pattern:
            nxt = self._goto[node].get(ch)
            if nxt is None:
                nxt = len(self._goto)
                self._goto.append({})
                self._fail.append(0)
                self._output.append([])
                self._goto[node][ch] = nxt
            node = nxt
        if pattern not in self._output[node]:
            self._output[node].append(pattern)
            self._patterns.append(pattern)
        return self

    def build(self):
        """Compute failure and output links (breadth-first). Call once, after all adds."""
        queue = deque()
        # depth-1 nodes fail to the root
        for ch, child in self._goto[0].items():
            self._fail[child] = 0
            queue.append(child)
        while queue:
            node = queue.popleft()
            for ch, child in self._goto[node].items():
                queue.append(child)
                # follow failure links from node's fail state to find child's fallback
                f = self._fail[node]
                while f != 0 and ch not in self._goto[f]:
                    f = self._fail[f]
                fallback = self._goto[f].get(ch, 0)
                if fallback == child:
                    fallback = 0
                self._fail[child] = fallback
                # inherit outputs from the failure state (shorter patterns ending here)
                self._output[child].extend(self._output[fallback])
        self._built = True
        return self

    def _next(self, node, ch):
        # transition with failure fallback
        while node != 0 and ch not in self._goto[node]:
            node = self._fail[node]
        return self._goto[node].get(ch, 0)

    def find(self, text):
        """Yield ``(end_index, pattern)`` for every pattern occurrence in ``text``."""
        if not self._built:
            self.build()
        node = 0
        for i, ch in enumerate(text):
            node = self._next(node, ch)
            if self._output[node]:
                for pat in self._output[node]:
                    yield (i, pat)

    def find_all(self, text):
        """Return a list of all ``(end_index, pattern)`` matches in ``text``."""
        return list(self.find(text))

    def contains_any(self, text):
        """True if any pattern occurs in ``text`` (stops at the first match)."""
        if not self._built:
            self.build()
        node = 0
        for ch in text:
            node = self._next(node, ch)
            if self._output[node]:
                return True
        return False

    def count_matches(self, text):
        """Total number of pattern occurrences in ``text`` (overlaps counted)."""
        return sum(1 for _ in self.find(text))

    @property
    def patterns(self):
        """The distinct non-empty patterns in the automaton."""
        return list(self._patterns)
