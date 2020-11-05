#ifndef IP_EXACT_H
#define IP_EXACT_H

#include "../ip.h"

typedef struct ip_exact_param
{
    long time_limit;
    int appr;
#define SOLVER_IP_EXACT_APPR_BAC 0
} ip_exact_param;

typedef struct ip_exact_stats
{
    struct stats_item *total;
} ip_exact_stats;

typedef struct ip_exact_bac_env ip_exact_bac_env;
typedef struct ip_exact_env
{
    int verbosity;
    ip_exact_param *param;
    ip_exact_stats *stats;
    ip_exact_bac_env *bac;
} ip_exact_env;

ip_exact_stats *
ip_create_exact_stats(void);
void
ip_free_exact_stats(ip_exact_stats **stats);

ip_exact_param *
ip_create_exact_param(void);
void
ip_free_exact_param(ip_exact_param **param);

ip_exact_env *
ip_create_exact_env(void);
void
ip_free_exact_env(ip_exact_env **env);

#include "bac/bac.h"

#endif
