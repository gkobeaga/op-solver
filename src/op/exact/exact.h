#ifndef OP_EXACT_H
#define OP_EXACT_H

#include "../op.h"

typedef struct cp_exact_param op_exact_param;
typedef struct cp_exact_stats op_exact_stats;
typedef struct cp_exact_env op_exact_env;

#define op_create_exact_prob(data) cp_create_exact_prob((data))
#define op_erase_exact_prob(prob) cp_erase_exact_prob((prob))
#define op_free_exact_prob(prob) cp_free_exact_prob((prob))

#define op_opt_exact(prob, env, sol) cp_opt_exact((prob), (env), (sol))

#define op_create_exact_env() cp_create_exact_env()
#define op_erase_exact_env(env) cp_erase_exact_env((env))
#define op_free_exact_env(env) cp_free_exact_env((env))
#define op_create_exact_param() cp_create_exact_param()
#define op_erase_exact_param(param) cp_erase_exact_param((param))
#define op_free_exact_param(param) cp_free_exact_param((param))
#define op_create_exact_stats() cp_create_exact_stats()
#define op_erase_exact_stats(stats) cp_erase_exact_stats((stats))
#define op_free_exact_stats(stats) cp_free_exact_stats((stats))

#include "bac/bac.h"

#endif
