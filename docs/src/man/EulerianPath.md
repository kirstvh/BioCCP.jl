# Eulerian Path
## Introduction 
An **Eulerian path** (or Eulerian trail) is a trail in a finite graph that visits every edge only once, allowing for revisiting nodes.

An **Eulerian circuit** (or Eulerian cycle) is an Eulerian path that starts and ends on the same node.

Euler proved that the necessary condition for the existance of Eulerian circuits is that all the nodes have an even degree.

**Euler's Theorem** : A connected graph has an Euler cycle if and only if every node has even degree. 

## Properities
- An undirected graph has an Eulerian trail if and only if exactly zero or two vertices have odd degree, and all of its vertices with nonzero degree belong to a single connected component.

- An undirected graph can be decomposed into edge-disjoint cycles if and only if all of its vertices have even degree. So, a graph has an Eulerian cycle if and only if it can be decomposed into edge-disjoint cycles and its nonzero-degree vertices belong to a single connected component.

- A directed graph has an Eulerian cycle if and only if every vertex has equal in degree and out degree, and all of its vertices with nonzero degree belong to a single strongly connected component. 
  
- Equivalently, a directed graph has an Eulerian cycle if and only if it can be decomposed into edge-disjoint directed cycles and all of its vertices with nonzero degree belong to a single strongly connected component.

## Explanation of the code
```@docs
create_adj_list
has_eulerian_cycle
```


## References:
[https://en.wikipedia.org/wiki/Eulerian_path#Applications](https://en.wikipedia.org/wiki/Eulerian_path#Applications)



