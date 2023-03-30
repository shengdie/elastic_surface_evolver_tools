Tools for create surface evolver surfaces and read surface evolver data (http://facstaff.susqu.edu/brakke/evolver/evolver.html)


# Any Mesh
For create surface meshes of any shape, currently support annular shape only. It has some general algorithms to help generate surfaces, like distinguishing boundary edges, writing mesh according to Surface Evolver syntax etc, so it can be easily extended to generate other shapes. 

## Examples
```
./create_ring.py -ir 1 -or 20 -es 0.1 -ds "lambda p: abs(ic.dist(p))**0.5/10 + 0.06" -kw "{'max_steps': 100000}" ring1_20
```
which creates a ring with inner radius 1 and outer radius 20. And the edge size is increasing with radius in a square root order.

## Requirements
```
dmsh==0.1.10
```


# Rectangular Mesh with Periodic Boundary Conditions
For create rectangular mesh with periodic condition.

# Surface Reader

Can read rectangular and annular surface evolver data generated pervious tools. It can plot surface, stress, strain, cross sections, contour etc.

## Examples
![fig](figs/surface_reader.png?raw=true "example of plotting stress")

## Requirements
```
h5py==3.7.0
```