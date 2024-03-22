# This script demonstrates how FSC can be used to find all paths in the network, including the shortest one. An implementation of Dijkstra's algorithm is provided for comparing the shortest path results and computing time. Obviously, Dijkstra is faster as it looks only for the shortest path.

import numpy as np
import matplotlib.pyplot as plt
from fsc import fsc
from typing import List, Dict
import heapq


class Solution:
    # Implementation for Dijkstra's shortest path algorithm
    def shortestPath(self, n: int, edges: List[List[int]], src: int) -> Dict[int, int]:
        adj = {}
        for i in range(n):
            adj[i] = []
            
        # s = src, d = dst, w = weight
        for s, d, w in edges:
            adj[s].append([d, w])

        # Compute shortest paths
        shortest = {}
        minHeap = [[0, src]]
        while minHeap:
            w1, n1 = heapq.heappop(minHeap)
            if n1 in shortest:
                continue
            shortest[n1] = w1

            for n2, w2 in adj[n1]:
                if n2 not in shortest:
                    heapq.heappush(minHeap, [w1 + w2, n2])
        
        # Fill in missing nodes
        for i in range(n):
            if i not in shortest:
                shortest[i] = -1

        return shortest
    
n = 5
edges = [[0,1,10], [0,2,3], [1,3,2], [2,1,4], [2,3,8], [2,4,2], [3,4,5]]
src = 0

a = Solution()
a.shortestPath(n=5, edges=edges, src=src)

R = np.array([[0,10,3,0,0],
            [10,0,4,2,0],
            [3,4,0,8,2],
            [0,2,8,0,5],
            [0,0,2,5,0]])
V = np.zeros((5,5))
V[3,3] = 1000

I = fsc(V=V, R=R).get_I()
plt.imshow(I)
