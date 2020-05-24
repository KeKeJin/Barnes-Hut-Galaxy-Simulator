# Barnes Hut Galaxy Simulator

This is a galaxy simulator implemented with Barnes Hut Algorithm.

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
