## MaximumFlow

Provides an implementation of the Ford-Fulkerson method with 
breadth-first-search (Edmonds-Karp algorithm) to solve a maximum flow problem.
Some additional necessary functions for the algorithm and visualisation are
provided as well.

More information on the problem and method can be found in the accompanying notebook or on
the [wikipedia page](https://en.wikipedia.org/wiki/Maximum_flow_problem).

```@docs
res_network
bfs
augmenting_path
isfeasible
clean!
```