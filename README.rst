percolate v0.1.0
================

This program searches clusters in a L x L matrix that
percolate. A cell of such a matrix can be either full or
empty. Every cell has 4 neighbor cells and builds a cluster
with them (and therefore with their neighbors' neighbors,
making it a recursive relationship). A cluster of free
cells percolates, if such a cluster starts from the top and
continues all the way down to the bottom of the matrix.


Build/install
-------------

If you have a copy of this project, open up a terminal and
type the following commands, in order to build percolate:

.. code:: bash

    cd path/to/the/project
    make

For more options and some help with building run:

.. code:: bash

    make help


Run percolate
-------------

After percolate is built, you can run it, simply by typing
``./percolate`` in your terminal (``percolate`` if
you installed it system wide)


Command line interface
----------------------

You can see the command line options of percolate if you
run ``./percolate -h``.
