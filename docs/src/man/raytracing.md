## Raytracing

Basic raytracing with Dijkstra's algorithm. Provides functions for creating a 2D scene with a circle, computing the shortest path from a source to a sink within the scene, and plotting the scene along with the shortest path(s).

### Scene creation
```@docs
create_scene
add_objects!
draw_circle
```

### Dijkstra's algorithm
```@docs
dijkstra
reconstruct_path
get_neighbors
```

### Plotting
```@docs
plot_pixels!
plot_pixel_edges!
plot_paths!
plot_circle!
```