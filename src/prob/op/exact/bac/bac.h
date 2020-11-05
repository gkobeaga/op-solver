#ifndef OP_EXACT_BAC_H
#define OP_EXACT_BAC_H

#include "../exact.h"

typedef struct cp_exact_bac_stats op_exact_bac_stats;
typedef struct cp_exact_bac_env op_exact_bac_env;
typedef struct cp_exact_bac_param op_exact_bac_param;

#define op_create_exact_bac_env() cp_create_exact_bac_env()
#define op_erase_exact_bac_env(env) cp_erase_exact_bac_env((env))
#define op_free_exact_bac_env(env) cp_free_exact_bac_env((env))
#define op_parse_exact_bac_args(argc, argv, env)                               \
    cp_parse_exact_bac_args((argc), (argv), (env))
#define op_create_exact_bac_param() cp_create_exact_bac_param()
#define op_erase_exact_bac_param(param) cp_erase_exact_bac_param((param))
#define op_free_exact_bac_param(param) cp_free_exact_bac_param((param))
#define op_create_exact_bac_stats() cp_create_exact_bac_stats()
#define op_erase_exact_bac_stats(stats) cp_erase_exact_bac_stats((stats))
#define op_free_exact_bac_stats(stats) cp_free_exact_bac_stats((stats))

#define op_opt_exact_bac(prob, env, sol) cp_opt_exact_bac((prob), (env), (sol))

int
op_init_exact_bac(op_prob *op, op_env *env);
#endif
