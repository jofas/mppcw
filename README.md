percolate v0.1.0
================

This program searches clusters in a L x L matrix that
percolate. A cell of such a matrix can be either full or
empty. Every cell has 4 neighbor cells and builds a cluster
with them (and therefore with their neighbors' neighbors,
making it a recursive relationship). A cluster of free
cells percolates, if such a cluster starts from the left
most column and continues all the way to the right most
column of the matrix.

This project contains a serial and a parallel version of
percolate v0.1.0.

The distributed version of percolate v0.1.0 uses MPI
for performing a distributed version of the clustering 
process. MPI must be install in order to run the parallel
version of percolate v0.1.0.


Build
-----

If you have a copy of this project, open up a terminal and
type the following commands, in order to build percolate:

```
cd path/to/the/project
make
```

This will result in two executables ```percolate_par``` and
```percolate_ser```. If only one of the two versions is
desired, run ```make percolate_[par|ser]```. While the use
of ```make``` or ```make all``` automatically cleans up the
project, ```make percolate_[par|ser]``` does not. Run
```make clean``` for cleaning up the project.

By default ```make``` will use Intel's Fortran compiler
(```ifort``` for the serial version and ```mpif08``` as the
MPI wrapper for the parallel version of percolate).
To change this behavior, change the ```FC``` and 
```FC_PAR``` flags in the ```Makefile``` to your preferred
compilers. For example, if you'd like to change to the
GNU Fortran compiler with OpenMPI's ```mpifort``` wrapper,
set ```FC=gfortran``` and ```FC_PAR=mpifort``` before 
running ```make```.


Run percolate
-------------

After percolate is built, you can run it, simply by typing
``./percolate`` in your terminal (``percolate`` if
you installed it system wide)


Command line interface
----------------------

You can see the command line options of percolate if you
run ``./percolate -h``.
