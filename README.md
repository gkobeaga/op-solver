# op-solver

> Algorithms for the Orienteering Problem

[![Actions Status](https://github.com/gkobeaga/op-solver/workflows/Build-Heur/badge.svg)](https://github.com/gkobeaga/op-solver/actions)
[![Actions Status](https://github.com/gkobeaga/op-solver/workflows/Build-Exact/badge.svg)](https://github.com/gkobeaga/op-solver/actions)
[![Lines](https://tokei.rs/b1/github/gkobeaga/op-solver)](https://github.com/XAMPPRocky/tokei)
[![Files](https://tokei.rs/b1/github/gkobeaga/op-solver?category=files)](https://github.com/XAMPPRocky/tokei)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://github.com/gkobeaga/op-solver/blob/master/LICENSE)

---

In this repository, you will find the implementation of two algorithms to solve
the Orienteering Problem (OP):

- **RB&C** (exact):
  ["A revisited branch-and-cut algorithm for large-scale orienteering problems"](https://arxiv.org/abs/2011.02743)
  by G. Kobeaga, M. Merino and J.A. Lozano
- **EA4OP** (heuristic):
  ["An evolutionary algorithm for the orienteering problem"](https://www.sciencedirect.com/science/article/abs/pii/S0305054817302241)
  by G. Kobeaga, M. Merino and J.A. Lozano

Both algorithms can be used to solve either small or large OP problems. Choose
between the heuristic or the exact algorithm depending on your needs.

<p align="center">
  <img src="https://user-images.githubusercontent.com/11088290/115864685-849d7400-a437-11eb-88d2-db94b76b3693.gif" alt="animated" width="52.5%" />
</p>

## Installation

First, obtain the source code,

```sh
git clone https://github.com/gkobeaga/op-solver
cd op-solver
```

install the dependencies,

```sh
sudo apt install autoconf automake libtool m4 libgmp-dev
```

and generate the _configure_ script.

```sh
./autogen.sh
mkdir -p build && cd build
```

Since the external LP solver used in the exact algorithm is proprietary software,
there are two options to install our software: to build only the heuristic
algorithm or to build both the heuristic and the exact algorithms.

### 1) Install Heuristic

By default, the solver is built only with the heuristic algorithm:

```sh
make clean
../configure
make
```

### 2) Install Heuristic and Exact

To build the exact algorithm, you need to have the [IBM ILOG CPLEX][2] installed
in your system. To build the `op-solver` with the exact algorithm:

```sh
make clean
../configure --with-cplex=<CPLEX_PATH>
#../configure --with-cplex=/opt/ibm/ILOG/CPLEX_Studio125/cplex/
make
```

## Usage

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

## Running on Docker

> Only available for the heuristic algorithm.

```sh
git clone https://github.com/gkobeaga/op-solver
cd op-solver
docker build -t op-solver .
mkdir tmp
docker run -v $PWD/tmp:/tmp -it --rm op-solver opt /OPLib/instances/gen3/kroA150-gen3-50.oplib
cat tmp/stats.json
```

## Output

By default, the results of the runs are written in a common `stats.json` file.
You can specify an alternative file to write the results:

```sh
./src/op-solver opt --stats my-stats.json OPLib/instances/gen3/kroA150-gen3-50.oplib
```


### Heuristic Output

The output of the evolutionary algorithm is split into two parts: the population initialization
part and the evolution part.

```json
{
  "prob": {
    "name": "kroA150",
    "n": 150,
    "d0": 13262
  },
  "sol": {
    "val": 4110,
    "cap": 13146,
    "sol_ns": 70,
    "lb": 4110,
    "ub": 1e+30,
    "cycle": [ 1, 47, 113, 84, 24, 38, 36, 127, 59, 141, 17, 15, 11, 32, 109, 91, 98, 23, 60,
               62, 20, 12, 86, 27, 149, 55, 83, 120, 115, 123, 43, 3, 46, 29, 132, 112, 107,
               30, 121, 101, 39, 78, 96, 52, 5, 37, 103, 146, 76, 13, 33, 95, 82, 116, 50, 73,
               68, 85, 135, 140, 117, 9, 7, 57, 51, 125, 61, 58, 105, 142 ]
  },
  "param": {
    "time_limit": 18000000,
    "init": 2,
    "select": 0,
    "pinit": 0
  },
  "stats": {
    "time": 205
  },
  "timestamp": 1618594637202,
  "event": "stats_summary",
  "env": "cp_init",
  "seed": 996021,
  "pid": 140336
}
{
  "prob": {
    "name": "kroA150",
    "n": 150,
    "d0": 13262
  },
  "sol": {
    "val": 5019,
    "cap": 13197,
    "sol_ns": 79,
    "lb": 5019,
    "ub": 1e+30,
    "cycle": [ 1, 130, 93, 28, 58, 61, 81, 25, 125, 51, 87, 145, 140, 135, 85, 68, 73, 114,
               144, 44, 50, 116, 82, 126, 95, 13, 76, 33, 146, 103, 37, 5, 52, 78, 96, 39,
               101, 121, 30, 107, 112, 132, 29, 46, 3, 14, 48, 100, 71, 41, 136, 128, 43,
               123, 115, 120, 149, 55, 83, 34, 117, 9, 7, 57, 20, 12, 27, 86, 35, 150, 62,
               60, 77, 110, 23, 98, 91, 109, 47 ]
  },
  "param": {
    "time_limit": 18000000,
    "it_lim": 2147483647,
    "pop_size": 100,
    "pop_stop": 25,
    "d2d": 50,
    "nparsel": 10,
    "pmut": 0.01,
    "len_improve1": 1,
    "len_improve2": 0
  },
  "stats": {
    "time": 748,
    "it": 750,
    "time_infeas_recover": 711
  },
  "timestamp": 1618594637951,
  "event": "stats_summary",
  "env": "cp_heur_ea",
  "seed": 996021,
  "pid": 140336
}
```

Note that, the two outputs share the same `pid` and
`seed`, which can be then used to obtain the total running time.

```sh
jq -s 'group_by(.seed, .pid)[] | .[1].stats.time += .[0].stats.time | .[0] * .[1]' stats.json
```

This [jq](https://github.com/stedolan/jq) command will merge the population
initialization parameters and evolutionary algorithm parameters, and sum the
running times.

```json
{
  "prob": {
    "name": "kroA150",
    "n": 150,
    "d0": 13262
  },
  "sol": {
    "val": 5019,
    "cap": 13197,
    "sol_ns": 79,
    "lb": 5019,
    "ub": 1e+30,
    "cycle": [ 1, 130, 93, 28, 58, 61, 81, 25, 125, 51, 87, 145, 140, 135, 85, 68, 73, 114,
               144, 44, 50, 116, 82, 126, 95, 13, 76, 33, 146, 103, 37, 5, 52, 78, 96, 39,
               101, 121, 30, 107, 112, 132, 29, 46, 3, 14, 48, 100, 71, 41, 136, 128, 43,
               123, 115, 120, 149, 55, 83, 34, 117, 9, 7, 57, 20, 12, 27, 86, 35, 150, 62,
               60, 77, 110, 23, 98, 91, 109, 47 ]
  },
  "param": {
    "time_limit": 18000000,
    "init": 2,
    "select": 0,
    "pinit": 0,
    "it_lim": 2147483647,
    "pop_size": 100,
    "pop_stop": 25,
    "d2d": 50,
    "nparsel": 10,
    "pmut": 0.01,
    "len_improve1": 1,
    "len_improve2": 0
  },
  "stats": {
    "time": 953,
    "it": 750,
    "time_infeas_recover": 711
  },
  "timestamp": 1618594637951,
  "event": "stats_summary",
  "env": "cp_heur_ea",
  "seed": 996021,
  "pid": 140336
}
```

### RB&C Output

The RB&C algorithm reports the following stats for each of the separation algorithms:

Stat | Description
---: | :---
 \*_active | Number of times that the separation algorithm was used
 \*_success | Number of times that the separation algorithm found at least a violated cut
 \*_total | Total number of violated cuts found by the separation algorithm
 \*_time | Total running time of the separation algorithm

```json
{
  "prob": {
    "name": "kroA150",
    "n": 150,
    "d0": 13262
  },
  "sol": {
    "val": 5039,
    "cap": 13246,
    "sol_ns": 79,
    "lb": 5039,
    "ub": 5039,
    "cycle": [ 1, 93, 28, 58, 61, 25, 81, 69, 64, 40, 54, 2, 144, 114, 44, 50, 116, 82, 126,
               95, 13, 76, 33, 146, 103, 37, 5, 52, 78, 96, 39, 101, 121, 30, 107, 112, 132,
               29, 46, 3, 14, 48, 100, 71, 41, 136, 128, 43, 123, 115, 120, 149, 55, 83, 34,
               135, 140, 125, 51, 87, 145, 9, 117, 7, 57, 20, 12, 27, 86, 150, 62, 60, 77,
               110, 23, 98, 91, 109, 47 ]
  },
  "param": {
    "sep_logical": 1,
    "sep_sec_comps": 1,
    "sep_sec_exact": 3,
    "sep_sec_cc_2": 0,
    "sep_sec_cc_extra": 1,
    "sep_blossom_fst": 0,
    "sep_blossom_eph": 1,
    "sep_blossom_egh": 1,
    "sep_cover_edge": 1,
    "sep_cover_vertex": 0,
    "sep_cover_cycle": 1,
    "sep_path": 1,
    "sep_loop": 1,
    "sep_srk_rule": 4,
    "sep_srk_s2": 0,
    "sep_srk_s3": 1,
    "sep_srk_extra": 1,
    "xheur_vph": 1,
    "xheur_vph_meta": 1
  },
  "stats": {
    "time": 37247,
    "sep_logical_active": 2267,
    "sep_logical_success": 290,
    "sep_logical_total": 538,
    "sep_logical_time": 910,
    "sep_sec_comps_active": 2267,
    "sep_sec_comps_success": 34,
    "sep_sec_comps_total": 727,
    "sep_sec_comps_time": 173,
    "sep_sec_exact_active": 937,
    "sep_sec_exact_success": 754,
    "sep_sec_exact_total": 6454,
    "sep_sec_exact_time": 3665,
    "sep_blossom_fast_active": 961,
    "sep_blossom_fast_success": 134,
    "sep_blossom_fast_total": 326,
    "sep_blossom_fast_time": 583,
    "sep_blossom_ghfast_active": 948,
    "sep_blossom_ghfast_success": 90,
    "sep_blossom_ghfast_total": 160,
    "sep_blossom_ghfast_time": 326,
    "sep_blossom_mst_active": 0,
    "sep_blossom_mst_success": 0,
    "sep_blossom_mst_total": 0,
    "sep_blossom_mst_time": 0,
    "sep_cover_edge_active": 507,
    "sep_cover_edge_success": 40,
    "sep_cover_edge_total": 40,
    "sep_cover_edge_time": 993,
    "sep_cover_cycle_active": 884,
    "sep_cover_cycle_success": 2,
    "sep_cover_cycle_total": 2,
    "sep_cover_cycle_time": 58,
    "sep_cover_vertex_active": 0,
    "sep_cover_vertex_success": 0,
    "sep_cover_vertex_total": 40,
    "sep_cover_vertex_time": 0,
    "sep_path_active": 504,
    "sep_path_success": 37,
    "sep_path_total": 298,
    "sep_path_time": 264,
    "sep_loop_time": 13407,
    "sep_loop_it_time": 0,
    "sep_loop_inner_time": 3075,
    "sep_loop_inner_it_time": 1167,
    "sep_loop_middle_time": 11909,
    "sep_loop_middle_it_time": 10262,
    "sep_loop_outer_time": 13407,
    "sep_loop_outer_it_time": 2728,
    "age_cut_time": 843,
    "age_vars_time": 1020,
    "add_vars_time": 7045,
    "add_cuts_time": 5834,
    "xheur_branch_time": 26,
    "xheur_sep_time": 2302
  },
  "timestamp": 1618593157989,
  "event": "stats_summary",
  "env": "cp_exact_bac",
  "seed": 696815,
  "pid": 134228
}
```

## Acknowledgments

The RB&C algorithm for large OP problems would not be possible without the
following implementations:

- TSP solver [Concorde](http://www.math.uwaterloo.ca/tsp/concorde.html)
- The B&C code for ["Solving the orienteering problem through branch-and-cut"](https://pubsonline.informs.org/doi/abs/10.1287/ijoc.10.2.133)
  (provided by Prof. JJ Salazar-González)

[1]: http://www.math.uwaterloo.ca/tsp/concorde.html
[2]: https://www.ibm.com/analytics/cplex-optimizer
