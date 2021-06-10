## Sudoku local search solver

Provide basic functions to solve a sudoku using Local search

Starts with a random solution, then randomly chooses a positon in the grid and change its value for a different number,
if the sudoku cost decrease, the new number is assigned to the current postion.
Additionally, the values of the resulted sudoku can be swapped randomly to decrease the cost of the solution   

```@docs
sudoku_greedydesc
flip
fliprow
makeflip
makefliprow
search
print_sudoku
```