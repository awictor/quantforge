"""Trie (prefix tree) for string keys: prefix search, counting, and autocomplete.

A trie stores a set of strings in a tree where each edge is a character, so every node
corresponds to a prefix shared by all keys beneath it. That makes prefix queries --
"is anything stored under this prefix?", "how many?", "list them" -- run in time
proportional to the query, independent of how many keys the trie holds. Each node caches a
subtree key count so `count_prefix` is ``O(len(prefix))``. Pure standard library.
"""


class _TrieNode:
    __slots__ = ("children", "is_end", "count")

    def __init__(self):
        self.children = {}
        self.is_end = False
        self.count = 0          # number of keys stored in this subtree


class Trie:
    """A prefix tree over string keys (a multiplicity-free set with prefix queries).

    `insert`, `contains`, and `delete` manage membership; `starts_with` tests whether any
    key has a given prefix; `count_prefix` counts keys under it; `keys_with_prefix` lists
    them; and `longest_prefix_of` finds the longest stored key that prefixes a query string.
    """

    __slots__ = ("_root", "_size")

    def __init__(self, words=None):
        self._root = _TrieNode()
        self._size = 0
        if words:
            for w in words:
                self.insert(w)

    def __len__(self):
        return self._size

    def __contains__(self, word):
        return self.contains(word)

    def insert(self, word):
        """Add ``word``. Returns ``True`` if newly inserted, ``False`` if already present."""
        node = self._root
        path = [node]
        for ch in word:
            nxt = node.children.get(ch)
            if nxt is None:
                nxt = _TrieNode()
                node.children[ch] = nxt
            node = nxt
            path.append(node)
        if node.is_end:
            return False
        node.is_end = True
        self._size += 1
        for n in path:
            n.count += 1
        return True

    def contains(self, word):
        """True if ``word`` is a stored key (exact match, not merely a prefix)."""
        node = self._find(word)
        return node is not None and node.is_end

    def _find(self, prefix):
        node = self._root
        for ch in prefix:
            node = node.children.get(ch)
            if node is None:
                return None
        return node

    def starts_with(self, prefix):
        """True if any stored key begins with ``prefix`` (an empty prefix matches any key)."""
        node = self._find(prefix)
        return node is not None and node.count > 0

    def count_prefix(self, prefix):
        """Number of stored keys that begin with ``prefix`` (``O(len(prefix))``)."""
        node = self._find(prefix)
        return node.count if node is not None else 0

    def keys_with_prefix(self, prefix):
        """Return the sorted list of stored keys beginning with ``prefix``."""
        node = self._find(prefix)
        if node is None:
            return []
        out = []
        self._collect(node, list(prefix), out)
        out.sort()
        return out

    def _collect(self, node, buf, out):
        if node.is_end:
            out.append("".join(buf))
        for ch, child in node.children.items():
            buf.append(ch)
            self._collect(child, buf, out)
            buf.pop()

    def keys(self):
        """Return all stored keys, sorted."""
        return self.keys_with_prefix("")

    def longest_prefix_of(self, word):
        """Return the longest stored key that is a prefix of ``word`` (``""`` if none)."""
        node = self._root
        best = ""
        buf = []
        if node.is_end:
            best = ""
        for ch in word:
            node = node.children.get(ch)
            if node is None:
                break
            buf.append(ch)
            if node.is_end:
                best = "".join(buf)
        return best

    def delete(self, word):
        """Remove ``word``. Returns ``True`` if it was present, ``False`` otherwise."""
        node = self._find(word)
        if node is None or not node.is_end:
            return False
        node.is_end = False
        self._size -= 1
        # decrement subtree counts and prune now-empty branches
        node = self._root
        node.count -= 1
        parents = [(None, None, node)]
        for ch in word:
            child = node.children[ch]
            child.count -= 1
            parents.append((node, ch, child))
            node = child
        # prune from the leaf upward while a node has no children and is not a key
        for parent, ch, cur in reversed(parents):
            if parent is None:
                break
            if not cur.children and not cur.is_end:
                del parent.children[ch]
            else:
                break
        return True
