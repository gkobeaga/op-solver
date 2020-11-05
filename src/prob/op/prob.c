#include "op-solver.h"
#include "op/op.h"

static void
op_init_prob_work(op_prob *op, solver_data *data)
{
    op->cap        = 0.0;
    op->data       = data;
    op->sol_status = SOLVER_UNDEF;
    return;
}

op_prob *
op_create_prob(solver_data *data)
{
    op_prob *op = malloc(sizeof(op_prob));
    op_init_prob_work(op, data);
    return op;
}
