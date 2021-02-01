#include <op-solver/op-solver.h>

int
main(void)
{
    solver_data *data;
    data = data_read("OPLib/instances/gen3/kroA150-gen3-50.oplib", 0);

    if (data->prob == SOLVER_PROB_OP)
    {
        op_env *op_env = op_create_env();
        op_prob *op    = op_create_prob(data);

        op_opt(op, op_env, op->sol);

        cp_print_sol(op, op->sol);

        op_free_prob(&op);
        op_free_env(&op_env);
    }

    return 0;
}
