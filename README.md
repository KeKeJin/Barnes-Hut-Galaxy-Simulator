# Barnes Hut Galaxy Simulator

This is a galaxy simulator implemented with Barnes Hut Algorithm.

## Getting Started

1. Install julia 1.4.0 ``` https://julialang.org/downloads/  ```
2. Clone this repository from Github: ```      git clone https://github.com/KeKeJin/Barnes-Hut-Galaxy-Simulator.git ```
3. In julia, add the following packages ``` using Pkg; Pkg.add("Revise");Pkg.add("CSV");Pkg.add("LinearAlgebra");Pkg.add("Distributions");Pkg.add("DataFrames");Pkg.add("Plots") ```
4. In julia, use the following packages ``` using Revise; using CSV; using LinearAlgebra; using Distributions; using DataFrames; using Plots ```
5. To start the simulation, run ``` includet("barnesHut.jl")```
4. The default simulation simulates 800 bodies in a space of 3300pc, in 100 thousand years. Specs can be changed  ``` initialGalaxy() ```
5. To start the simulation, run ``` includet("Visualization.jl")```
6. To generate a gif, run ``` plotting()```
