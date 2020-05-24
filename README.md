# Barnes Hut Galaxy Simulator

This is a galaxy simulator implemented with Barnes Hut Algorithm.

## Background
Barnes-Hut Algorithm is a novel and clever algorithm in solving the famous n-body problem, which is used to predict the motion of individual objects under the influence of other objects in a cosmological setting. By the usage of tree data structure and grouping nearby bodies, Barnes-Hut Algorithm reduces the time complexity from the brute force O(n^2) to O(nlogn).
Read more: [![wikipedia](https://en.wikipedia.org/wiki/Barnes%E2%80%93Hut_simulation)](https://en.wikipedia.org/wiki/Barnes%E2%80%93Hut_simulation)

## Installation


Clone this repository from Github
```      
git clone https://github.com/KeKeJin/Barnes-Hut-Galaxy-Simulator.git
```

## Usage

To start the simulation, run
```
include("simulation.jl")
```
Theta is set to 0.7, to change it
```
changeTheta(newTheta)
```

### A Random Galaxy
The default simulation simulates 100 bodies in a space of 3300pc, while 1 second in the simulation represents 10 thousand years.  
```
randomGalaxy([, num => number of bodies in the simulator])
```
The simulator writes CSV files in ```\data``` and erases them later.
![random.gif](exampleSimulation/random.gif)
### A Disk-Shaped Galaxy
The default simulation simulates 100 bodies in a disk in the plane x+y=0, while 1 second in the simulation represents 10 thousand years.
```
diskGalaxy([, normal => normal vector of the disk plane,
            num => number of bodies in the simulator])
```
The simulator writes CSV files in ```\data``` and erases them later.
![diskGalaxy.gif](exampleSimulation/diskGalaxy.gif)
### A Line-Shaped Galaxy
The default simulation simulates 100 bodies in a line of in the direction of (1,1,0), while 1 second in the simulation represents 10 thousand years.
```
lineGalaxy([, dir => direction vector of the line,
            num => number of bodies in the simulator])
```
The simulator writes CSV files in ```\data``` and erases them later.

### A System of Two Galaxies
The default simulation simulates a disk-shaped galaxy of 100 bodies in the plane of x+y=0 moving towards a disk-shaped galaxy of 100 bodies in the plane of x+y+z=0
```
twoCollapseGalaxy([,num1 -> number of bodies in the first galaxy,
                  num2 -> number of bodies in the second galaxy,
                  velocity1 -> the velocity of the first galaxy,
                  velocity2 -> the velocity of the second galaxy,
                  type1 -> the type of the first galaxy,
                  type2 -> the type of the second galaxy,
                  vector1 -> the normal vector or direction vector of the first galaxy,
                  vector2 -> the normal vector or direction vector of the second galaxy]
```
![twoGalaxies.gif](exampleSimulation/twoGalaxies.gif)

## Todos
1. Two Galaxies collapsing needs the user to define a better initial condition than the default
2. Add progress bar while the program is running
3. Add interactive interface to view the simulation from different angles

## Credits
Thanks to Professor Vatche Sahakian for inspiring this project and providing feedbacks.
Thanks to my friend Jonathan Hayase for supporting and brainstorming.
Thanks to my Lord for the amazing universe.
