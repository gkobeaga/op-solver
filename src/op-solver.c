#include "op-solver.h"

int
main(int argc, char *argv[])
{
    int rval = 0;

    if (!strcmp(argv[1], "opt"))
    {
        __solver_opt(argc, argv);
    }
    else
    {
        exit(1);
    }

    return rval;
}
