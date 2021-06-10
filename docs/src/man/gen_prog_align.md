# GenProgAlign
The package GenProgAlign contains utilities to perform multiple sequence alignment 
with genetic programming. The algorithm is inspired by SAGA (Notredame, Higgins, 1996)$^1$ 
but contains only a subset of all the suggested operations and doesn't use phylogenetic information
generated with clustal for the objective function nor the operations.. The function ga\_alignment performs the genetic algorithm, all other functions exported
by GenProgAlign are helper functions used by ga\_alignment. These are exported to enable the user to 
explore the behaviour of the algorithm.


```@docs
GenProgAlign.ga_alignment
GenProgAlign.fasta_to_array
GenProgAlign.Alignment
GenProgAlign.AA
GenProgAlign.pairwise_score
GenProgAlign.objective_function
GenProgAlign.one_point_crossover
GenProgAlign.gap_insertion
GenProgAlign.shifting_one_seq
GenProgAlign.shifting_block
```
1. Cédric Notredame, Desmond G. Higgins, SAGA: Sequence Alignment by Genetic Algorithm, Nucleic Acids Research, Volume 24, Issue 8, 1 April 1996, Pages 1515–1524, https://doi.org/10.1093/nar/24.8.1515
