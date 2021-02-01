#ifndef OP_HEUR_H
#define OP_HEUR_H

#include "../op.h"

typedef struct cp_heur_stats op_heur_stats;
typedef struct cp_heur_env op_heur_env;
typedef struct cp_heur_param op_heur_param;

#define op_create_heur_env() cp_create_heur_env()
#define op_erase_heur_env(env) cp_erase_heur_env((env))
#define op_free_heur_env(env) cp_free_heur_env((env))
#define op_create_heur_param() cp_create_heur_param()
#define op_erase_heur_param(param) cp_erase_heur_param((param))
#define op_free_heur_param(param) cp_free_heur_param((param))
#define op_create_heur_stats() cp_create_heur_stats()
#define op_erase_heur_stats(stats) cp_erase_heur_stats((stats))
#define op_free_heur_stats(stats) cp_free_heur_stats((stats))

#define op_opt_heur(prob, env, sol) cp_opt_heur((prob), (env), (sol))

#include "ea/ea.h"

#endif
