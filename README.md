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
process. MPI must be installed in order to run the parallel
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

* ```percolate_ser```

  After ```percolate_ser``` is built, you can run it, 
  simply by typing ``./percolate_ser`` in your terminal.

* ```percolate_par```

  Depending on your MPI implmentation, but in general
  executing ```percolate_par``` should look something like:
  ```mpiexec -n 4 ./percolate_par```.


Command line interface
----------------------

Running ```./percolate_[par|ser] -h``` will show you in 
detail the command line interface.

The output of the ```-h``` option:

```
 percolate v0.1.0                                         
                                                          
 Program for computing, whether any cluster in a          
 randomly generated matrix percolates.                    
                                                          
 Usage: percolate_[par|ser] [seed] [options ...]          
                                                          
                                                          
 seed: INT Seed for the random number generator.          
           DEFAULT: 1564                                  
                                                          
                                                          
 options:                                                 
                                                          
     -h, --help                    Prints this help mes-  
                                   sage.                  
                                                          
         --version                 Prints the version of  
                                   this program.          
                                                          
     -l, --length            INT   Sets the dimension of  
                                   the matrix.            
                                   DEFAULT: 20.           
                                                          
     -d, --density           FLOAT Sets the density of    
                                   the full cells.        
                                   DEFAULT: 0.4.          
                                                          
     -p, --print_n_clusters  INT   Sets the amount of     
                                   clusters displayed in  
                                   the .pgm file.         
                                   DEFAULT: max.          
                                   Can display a maximum  
                                   of the 9 biggest       
                                   clusters.              
                                                          
         --pgm_file_path     PATH  Sets the path for the  
                                   .pgm file.             
                                   DEFAULT: map.pgm       
                                                          
         --print_iter_factor FLOAT Sets the factor of how 
                                   often the average cell 
                                   value during cluster-  
                                   ing is printed.        
                                   iter % int(L * factor) 
                                   DEFAULT: 0.5           
                                   This flag is irrele-   
                                   vant for the serial    
                                   version of percolate.  
```


Project structure
-----------------

* src/

  contains the source code.

  + ```cart_comm.f90```

    the communicator object for the 2d cartesian domain
    decomposition in ```percolate_par```. Used for all
    communications between the MPI processes (scattering,
    gathering and halo swapping). 

  + ```cli_info.f90```

    utility file containing the version of percolate and
    the output for the ```-h``` flag.

  + ```io.f90```

    for reading the command line arguments and writing the
    .pgm output file.

  + ```uni.f90```

    pseudo random number generator.

  + ```percolate_par.f90```

    program of ```percolate_par```.

  + ```percolate_ser.f90```

    program of ```percolate_ser```.

* ```test.sh```

  test suite for regression testing ```percolate_par```
  against output from ```percolate_ser```.

* ```test.pbs```

  wrapper around ```test.sh``` for executing the tests in
  a job handled by a job scheduler.


Note on the decomposition used in the parallel version
------------------------------------------------------

You can execute ```percolate_par``` with any amount of 
MPI processes. However, ```percolate_par``` does not 
guarantee that the program will run with the defined amount
of processes. 

This is best explained by an example. Imagine you want to
percolate a 2 x 2 matrix on three MPI processes with
executing ```mpiexec -n 3 percolate -l 2```. 
In this case, ```percolate_par``` would decompose its 
dimensions to 3 x 1, which means 3 processes are stacked
horizontically (in the 2d cartesian topology). Now,
```percolate_par``` would distribute the matrix by giving
the first ```n - 1```---```n``` in this case is 
3---```floor(l/n)``` many elements in the x-direction and
```floor(l/1)``` many elements in the y-direction.
In this case, the first and second process would both have
a chunk of the matrix of the size 0 x 2.

```percolate_par``` prevents empty chunks by killing MPI
processes, if a dimension of the cartesian topology has 
more processes than the matrix elements on the 
corresponding axis. So, in this example, 
```percolate_par``` would actually run with a 2 x 1 
decomposition, so every process has at least one element.
