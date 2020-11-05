#ifndef OP_H
#define OP_H

#include "../../op-solver.h"
#include "../cp/cp.h"

typedef struct cp_sol op_sol;

typedef struct cp_pop op_pop;

typedef struct cp_prob op_prob;

typedef struct cp_param op_param;
typedef struct cp_stats op_stats;
typedef struct cp_env op_env;

#define op_create_prob(data) cp_create_prob((data))
#define op_erase_prob(prob) cp_erase_prob((prob))
#define op_free_prob(prob) cp_free_prob((prob))

#define op_opt(prob, env, sol) cp_opt((prob), (env), (sol))

#define op_create_env() cp_create_env()
#define op_erase_env(env) cp_erase_env((env))
#define op_free_env(env) cp_free_env((env))
#define op_create_param() cp_create_param()
#define op_erase_param(param) cp_erase_param((param))
#define op_free_param(param) cp_free_param((param))
#define op_create_stats() cp_create_stats()
#define op_erase_stats(stats) cp_erase_stats((stats))
#define op_free_stats(stats) cp_free_stats((stats))

#define op_create_sol(prob) cp_create_sol((prob))
#define op_copy_sol(insol, outsol) cp_copy_sol((insol), (outsol))
#define op_erase_sol(sol) cp_erase_sol((sol))
#define op_free_sol(sol) cp_free_sol((sol))

#define op_create_pop(prob, size) cp_create_pop((prob), (size))
#define op_set_pop_sol(pop, sol, pos) cp_set_pop_sol((pop), (sol), (pos))
#define op_update_pop(pop) cp_update_pop((pop))
#define op_erase_pop(pop) cp_erase_pop((pop))
#define op_free_pop(pop) cp_free_pop((pop))

#include "exact/exact.h"
#include "heur/heur.h"
#include "init/init.h"

#endif
