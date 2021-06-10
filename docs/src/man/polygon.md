## polygon

Provides the triangleevolution algorithm, an evolutionary algorithm that tries to replicate a target image with a bunch of triangles.

As an evolutionary algorithm, it makes use of a population of individuals.
An individual is defined as an array consisting of the following:
- An array of triangles (a triangle is simply 3 xy coordinates of points and a color, not yet an image)
- An array of RGB values, forming an image made from drawing all these triangles onto a white canvas (array of triangles converted to an image)
- The score of the individual, which is too say how much its image differs from the target image.

We try to get an individual with an image as close as possible to the target image.
In order to do this, the algorithm follows these steps:
1) Generate a random starting population.
2) Make a parent population from your population ("the survivors")
    - Pass a given fraction consisting of the population's best indivuals (the elite) on to the parent population
    - Make the next part of the parent population by means of tournament selection (select 2 individuals at random, best one gets into the parent population)
    - Randomly generate new individuals for the last part of the parent population.
3) Generate a child population from the parent population, new individuals with genes (= triangles) inherited from their parents
4) Randomly mutate some genes of the child population (change position and colors of the triangles)
5) This child population becomes the new starting population
6) Repeat step 2 -> 5 for a given amount of times, the amount of generations.

```@docs
triangleevolution
```

The algorithm makes use of the following:

```@docs
Shape
Triangle
samesideofline
in
getboundaries
drawtriangle
drawimage
colordiffsum
checktriangle
generatetriangle
generatepopulation
triangletournament
mutatetriangle
mutatepopulation
makechildpopulation
```