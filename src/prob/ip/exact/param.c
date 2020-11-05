#include "ip/exact/exact.h"
#include "op-solver.h"

ip_exact_param *
ip_create_exact_param(void)
{
    ip_exact_param *param = malloc(sizeof(ip_exact_param));
    param->time_limit     = 5 * 60 * 60 * 1000;
    return param;
}

void
ip_free_exact_param(ip_exact_param **param)
{
    if (*param)
    {
        free(*param);
        *param = NULL;
    }
}
