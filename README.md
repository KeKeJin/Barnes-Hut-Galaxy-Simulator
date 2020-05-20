# Barnes Hut Galaxy Simulator

This is a galaxy simulator implemented with Barnes Hut Algorithm.

## Getting Started

1. Install julia 1.4.0 ``` https://julialang.org/downloads/  ```
2. Clone this repository from Github: ```      git clone https://github.com/KeKeJin/.git ```

### To run the Barnes Hut Simulator
1. In julia, add the following packages ``` using Pkg; Pkg.add("Revise");Pkg.add("CSV");Pkg.add("LinearAlgebra");Pkg.add("Distributions");Pkg.add("DataFrames") ```
2. In julia, use the following packages ``` using Revise; using CSV; using LinearAlgebra; using Distributions; using DataFrames ```
3. To start the simulation, run ``` includet("barnesHut.jl")```
4. The default simulation simulates 80 bodies in a space of 3300pc, in 150 days. Specs can be changed  ``` initialGalaxy() ```

### To run the Visualization
1. In julia, add the following packages ``` using Pkg; Pkg.add("Revise");Pkg.add("CSV");Pkg.add("DataFrames");Pkg.add("Plots") ```
2. In julia, use the following packages ``` using Revise; using CSV; using DataFrames; using Plots ```
3. To start the simulation, run ``` includet("Visualization.jl")```
4. To generate a gif, run ``` plotting()```
