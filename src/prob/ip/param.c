#include "ip/ip.h"
#include "op-solver.h"

ip_param *
ip_create_param(void)
{
    ip_param *param   = malloc(sizeof(ip_param));
    param->time_limit = 5 * 60 * 60 * 1000;
    return param;
}

void
ip_free_param(ip_param **param)
{
    if (*param)
    {
        free(*param);
        *param = NULL;
    }
}
