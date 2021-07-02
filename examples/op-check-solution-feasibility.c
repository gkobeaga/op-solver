#include <op-solver/op-solver.h>
#include <op-solver/op/op.h>

int
main(void)
{
    solver_data *data;
    data = data_read("test.oplib", 0);

    if (data->prob == SOLVER_PROB_OP)
    {
        op_env *op_env = op_create_env();
        op_prob *op    = op_create_prob(data);

        int ns        = 15;
        int cycle[15] = {0,  13, 46, 20, 35, 27, 7, 18,
                         43, 45, 15, 12, 11, 23, 14};

        op_sol *sol = cp_get_sol_from_cycle(op, ns, cycle);

        cp_print_sol(op, sol);

        if (sol->cap > op->cap)
            printf(
            "\nInfeasible solution: solution cap %.0f > problem cap %0.f\n",
            sol->cap, op->cap);
        else
            printf(
            "\nFeasible solution: solution cap %.0f <= problem cap %0.f\n",
            sol->cap, op->cap);

        op_free_sol(&sol);
        op_free_prob(&op);
        op_free_env(&op_env);
    }

    return 0;
}

/*
NAME : att48
COMMENT : 48 capitals of the US (Padberg/Rinaldi)
TYPE : OP
DIMENSION : 48
COST_LIMIT : 2657
EDGE_WEIGHT_TYPE : ATT
NODE_COORD_SECTION
1 6823 4674
2 7692 2247
3 9135 6748
4 7721 3451
5 8304 8580
6 7501 5899
7 4687 1373
8 5429 1408
9 7877 1716
10 7260 2083
11 7096 7869
12 6539 3513
13 6272 2992
14 6471 4275
15 7110 4369
16 6462 2634
17 8476 2874
18 3961 1370
19 5555 1519
20 4422 1249
21 5584 3081
22 5776 4498
23 8035 2880
24 6963 3782
25 6336 7348
26 8139 8306
27 4326 1426
28 5164 1440
29 8389 5804
30 4639 1629
31 6344 1436
32 5840 5736
33 5972 2555
34 7947 4373
35 6929 8958
36 5366 1733
37 4550 1219
38 6901 1589
39 6316 5497
40 7010 2710
41 9005 3996
42 7576 7065
43 4246 1701
44 5906 1472
45 6469 8971
46 6152 2174
47 5887 3796
48 7203 5958
NODE_SCORE_SECTION
1 0
2 1
3 1
4 1
5 1
6 1
7 1
8 1
9 1
10 1
11 1
12 1
13 1
14 1
15 1
16 1
17 1
18 1
19 1
20 1
21 1
22 1
23 1
24 1
25 1
26 1
27 1
28 1
29 1
30 1
31 1
32 1
33 1
34 1
35 1
36 1
37 1
38 1
39 1
40 1
41 1
42 1
43 1
44 1
45 1
46 1
47 1
48 1
DEPOT_SECTION
 1
 -1
EOF
 */
