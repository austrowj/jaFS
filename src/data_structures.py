#pyright: basic
from typing import Any

from sortedcontainers import SortedList

class FenwickTreeWithSortedLists:
    """
    Data structure needed to compute win ratio Z-score efficiently.
    """

    def __init__(self, max: int):
        self.nodes: list[SortedList | None] = [None]*(max + 1)
    
    def __len__(self):
        return len(self.nodes)

    def add(self, x: Any, t: int):
        assert 0 < t < len(self.nodes)

        while t < len(self.nodes):
            if self.nodes[t] is None:
                self.nodes[t] = SortedList()
            
            node: SortedList = self.nodes[t]
            assert node is not None
            node.add(x)
            t += t&-t
    
    def query_range(self, x: Any, t_lb: int, t_ub: int):
        """Returns the number of elements less than x and the number of elements greater than x in the nodes between t_lb and t_ub."""
        a = self.query(x, t_ub)
        b = self.query(x, t_lb)
        return a[0]-b[0], a[1]-b[1]

    def query(self, x: Any, t: int):
        """Returns the number of elements less than x and the number of elements greater than x in the first t nodes."""
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