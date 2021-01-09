# op-solver
> Algorithms for the Orienteering Problem

[![Actions Status](https://github.com/gkobeaga/op-solver/workflows/build-test/badge.svg)](https://github.com/gkobeaga/op-solver/actions)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://github.com/gkobeaga/op-solver/blob/master/LICENSE)

----

In this repository you will find the implementation of two algorithms to solve the Orienteering Problem (OP):
  - **RB&C** (exact):
 ["A revisited branch-and-cut algorithm for large-scale orienteering problems"](https://arxiv.org/abs/2011.02743) by G. Kobeaga, M. Merino and J.A. Lozano
  - **EA4OP** (heuristic): ["An evolutionary algorithm for the orienteering problem"](https://www.sciencedirect.com/science/article/abs/pii/S0305054817302241) by G. Kobeaga, M. Merino and J.A. Lozano

Both algorithms can be used to solve either small or large OP problems. Choose between the heuristic or the exact algorithm depending on your needs.

Installation
------------

First, obtain the source code,
```sh
git clone https://github.com/gkobeaga/op-solver
cd op-solver
```

install the dependencies,
```sh
sudo apt install autoconf automake libtool m4 libgmp-dev
```


and generate the *configure* script.
```sh
./autogen.sh
mkdir -p build && cd build
```

Since the external LP solver used in the exact algorithm is proprietary software, there are two options to install our software: to build only the heuristic algorithm or to build both the heuristic and the exact algorithms.

#### 1) Install Heuristic
By default, the solver is built only with the heuristic algorithm:
```sh
make clean
../configure
make
```

#### 2) Install Heuristic and Exact

To build the exact algorithm, you need to have the [IBM ILOG CPLEX][2] installed in your system.
To build the `op-solver` with the exact algorithm:

```sh
make clean
../configure --with-cplex=<CPLEX_PATH>
#../configure --with-cplex=/opt/ibm/ILOG/CPLEX_Studio125/cplex/
make
```

Usage
-------------
Download first the benchmark instances for the OP:
```sh
cd build
git clone https://github.com/bcamath-ds/OPLib.git
```

To solve the problem using the EA4OP algorithm:
```sh
./src/op-solver opt --op-exact 0 OPLib/instances/gen3/kroA150-gen3-50.oplib
```

To solve the OP using the revisited Branch-and-Cut algorithm(RB\&C):
```sh
./src/op-solver opt --op-exact 1 OPLib/instances/gen3/kroA150-gen3-50.oplib
```

You can increase the verbosity of the RB\&C with:
```sh
./src/op-solver opt --op-exact 1 --op-exact-bac-verbose 1 OPLib/instances/gen3/kroA150-gen3-50.oplib
```

Acknowledgments
---------------
The RB&C algorithm for large OP problems would not be possible without the following implementations:
  - TSP solver [Concorde](http://www.math.uwaterloo.ca/tsp/concorde.html)
  - The B&C code for ["Solving the orienteering problem through branch-and-cut"](https://pubsonline.informs.org/doi/abs/10.1287/ijoc.10.2.133) (provided by Prof. JJ Salazar-González)

[1]:http://www.math.uwaterloo.ca/tsp/concorde.html
[2]:https://www.ibm.com/analytics/cplex-optimizer
