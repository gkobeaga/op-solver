# op-solver
#### Algorithms for the Orienteering Problem

In this repository you will find the implementation of two algorithms to solve the Orienteering Problem (OP):
  - **RB&C** (exact):
 ["A revisited branch-and-cut algorithm for large-scale orienteering problems"]() by G. Kobeaga, M. Merino and J.A. Lozano
  - **EA4OP** (heuristic): ["An evolutionary algorithm for the orienteering problem"](https://www.sciencedirect.com/science/article/abs/pii/S0305054817302241) by G. Kobeaga, M. Merino and J.A. Lozano

Both algorithms can be used to solve either small or large OP problems. Choose between the heuristic or the exact algorithm depending on your needs.

Instructions to build
---------------------

Install the dependencies:
```sh
sudo apt-get install libtool m4 libgmp-dev
```
You also need to have installed [IBM CPLEX Optimizer][2] in your system.

Obtain the source code by:
```sh
git clone https://github.com/gkobeaga/op-solver
cd op-solver
```

Set, if needed, the path to the cplex directory. By default, it searches in `/opt/ibm/ILOG/CPLEX_Studio125/cplex`.
```sh
export CPLEX_PATH=<path/to/cplex>
```

Use GNU Autotools to generate the *configure* script:
```sh
./autogen.sh
```

To build and install the binary type:
```sh
mkdir -p build && cd build
../configure
make
make install
```

How to use it
-------------
First, download the benchmark instances for the OP:
```sh
cd build
git clone https://github.com/bcamath-ds/OPLib.git
```

By default the solver uses the revisited Branch-and-Cut algorithm(RB&C):
```sh
./src/op-solver opt OPLib/instances/gen2/rd100-gen2-50.oplib
```

You can increase the verbosity of the RB&C with:
```sh
./src/op-solver opt--op-exact-bac-verbose 1 OPLib/instances/gen2/rd100-gen2-50.oplib
```

In order to solve the problem using the EA4OP algorithm, add the `--op-exact 0` argument:
```sh
./src/op-solver opt --op-exact 0 OPLib/instances/gen2/rd100-gen2-50.oplib
```

Acknowledgments
---------------
The RB&C algorithm for large OP problems would not be possible without the following implementations:
  - TSP solver [Concorde](http://www.math.uwaterloo.ca/tsp/concorde.html)
  - The B&C code for ["Solving the orienteering problem through branch-and-cut"](https://pubsonline.informs.org/doi/abs/10.1287/ijoc.10.2.133) (provided by Prof. JJ Salazar-González)

[1]:http://www.math.uwaterloo.ca/tsp/concorde.html
[2]:https://www.ibm.com/analytics/cplex-optimizer
