#pyright: basic
from sortedcontainers import SortedList

class FenwickTreeWithSortedLists:
    """
    Data structure needed to compute win ratio Z-score efficiently.
    """

    def __init__(self, max: int):
        self.nodes: list[SortedList | None] = [None]*(max+1)

    def add(self, x, t):
        assert 0 < t < len(self.nodes)

        while t < len(self.nodes):
            if self.nodes[t] is None:
                self.nodes[t] = SortedList()
            
            node: SortedList = self.nodes[t]
            assert node is not None
            node.add(x)
            t += t&-t
    
    def query(self, x, t_lb, t_ub):
        a = self._query_single(x, t_ub)
        b = self._query_single(x, t_lb)
        return a[0]-b[0], a[1]-b[1]

    def _query_single(self, x, t):
        total_a, total_b = 0, 0
        while t > 0:
            node: SortedList | None = self.nodes[t] # pyright: ignore (need stub file to prevent complaints here)
            if node is not None:
                total_a += node.bisect_left(x)
                total_b += len(self.nodes[t]) - node.bisect_right(x)
            t -= t&-t
        
        return total_a, total_b

    # Just for testing
    def _print(self, all=False):
        print()
        for i, node in filter(lambda x: (x[1] != None) or all, enumerate(self.nodes)):
            print(f'{i}: {node}')  