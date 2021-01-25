#include "op-solver.h"
#include <sys/stat.h>

unsigned long __seed__;

struct cmd_args
{
    const char *data_file;
    const char *stats_file;
    char sol_file[50];
};

static void
print_version_and_args(int argc, char *argv[]);
int
parse_opt_args_(int argc, char *argv[], struct cmd_args *cmd_args);

//#define RNG_SEED 947107

int
__solver_opt(int argc, char *argv[])
{
    int rval                 = 0;
    solver_data *data        = NULL;
    stats_item *total        = NULL;
    struct cmd_args cmd_args = {NULL};

    total = stats_create("Total");
    stats_start(total);

    /*--------------------------------------------------------------------------*/
    /* Set initial seed */
#ifdef RNG_SEED
    __seed__ = RNG_SEED;
#else
    struct timespec tm_seed;
    clock_gettime(CLOCK_REALTIME, &tm_seed);
    __seed__ = tm_seed.tv_nsec / 1000;
#endif

    /*--------------------------------------------------------------------------*/
    /* print version and arguments */
    print_version_and_args(argc, argv);

    /*--------------------------------------------------------------------------*/
    /* parse command-line parameters */

    rval = parse_opt_args_(argc, argv, &cmd_args);
    check_rval(rval, "", done);

    /*--------------------------------------------------------------------------*/
    /* Initialize random generators */
    srand(__seed__);
    srand48(__seed__);

    /* read data from file*/
    if (cmd_args.data_file == NULL)
    {
        printf("No input problem file specified; try %s --help\n", argv[0]);
        rval = EXIT_FAILURE;
        goto done;
    }
    data = data_read(cmd_args.data_file, 0);
    check_null(data, "data file processing error\n", done);

    /* Solve problem */
    switch (data->prob)
    {
    case SOLVER_PROB_OP: {
        op_prob *op    = NULL;
        op_env *op_env = op_create_env();
        if (cmd_args.stats_file)
            op_env->exact->bac->stats->file = cmd_args.stats_file;
        if (cmd_args.sol_file)
            strcpy(op_env->sol_file, cmd_args.sol_file);
        cp_parse_args(argc, argv, op_env);
        op = op_create_prob(data);
        op_opt(op, op_env, op->sol);
        cp_print_sol(op, op->sol);
        op_free_prob(&op);
        op_free_env(&op_env);
        break;
    }
    default:
        printf("Unknown problem %d.\n", data->prob);
        return 1;
    }

    rval = 0;

    stats_stop(total, 1);
    printf("Total Running Time: %ld (milliseconds)\n",
           stats_get_current_time(total));

done:

    if (data)
        data_free(&data);
    if (total)
        free(total);

    return rval;
}

static void
print_version_and_args(int argc, char *argv[])
{ /* print version information */
    printf("Solver v%d.%d\n", SOLVER_MAJOR_VERSION, SOLVER_MINOR_VERSION);
    if (argc > 1)
    {
        int k;
        printf("  parameter(s):");
        for (k = 1; k < argc; k++)
        {
            if (argc == 3 || k != argc - 1)
                printf(" %s", argv[k]);
            else
                printf("\n\033[%dC %s", 19, argv[k]);
        }
        printf("\n");
    }
    printf("\n");
    return;
}

int
parse_opt_args_(int argc, char *argv[], struct cmd_args *cmd_args)
{
    int k;

    for (k = 2; k < argc; k++)
    {
        if (!strcmp(argv[k], "--seed"))
        {
            k++;
            if (k == argc || argv[k][0] == '\0' || argv[k][0] == '-')
            {
                printf("No stats file specified\n");
                return 1;
            }
            if (cmd_args->stats_file != NULL)
            {
                printf("Only one stats file allowed\n");
                return 1;
            }
            __seed__ = (unsigned long)atoi(argv[k]);
        }
        else if (!strcmp(argv[k], "--stats"))
        {
            k++;
            if (k == argc || argv[k][0] == '\0' || argv[k][0] == '-')
            {
                printf("No stats file specified\n");
                return 1;
            }
            if (cmd_args->stats_file != NULL)
            {
                printf("Only one stats file allowed\n");
                return 1;
            }
            cmd_args->stats_file = argv[k];
        }
        else if (!strcmp(argv[k], "--sol"))
        {
            k++;
            if (k == argc || argv[k][0] == '\0' || argv[k][0] == '-')
            {
                printf("No solution file specified\n");
                return 1;
            }
            strcpy(cmd_args->sol_file, argv[k]);
        }
        else if ((argv[k][0] == '-' && argv[k][1] == '-' && argv[k][2] == 'o' &&
                  argv[k][3] == 'p'))
        {
            k++;
        }
        else if (argv[k][0] == '-' || (argv[k][0] == '-' && argv[k][1] == '-'))
        {
            printf("Invalid option '%s'; try %s --help\n", argv[k], argv[0]);
        }
        else
        {
            if (cmd_args->data_file != NULL)
            {
                printf("Only one input problem file allowed (%s)\n",
                       cmd_args->data_file);
                return 1;
            }
            cmd_args->data_file = argv[k];

            // if (cmd_args->sol_file == NULL)
            {
                char base[50];
                char seed[50];
                strcpy(base, "solutions/");
                mkdir(base, S_IRWXU);
                strcat(base, strrchr(argv[k], '/') + 1);
                const int n = snprintf(NULL, 0, "%lu", __seed__);
                if (!snprintf(seed, n + 1, "%lu", __seed__))
                    exit(1);
                memcpy(cmd_args->sol_file, base, strlen(base) - 6);
                strcat(cmd_args->sol_file, "-");
                strcat(cmd_args->sol_file, seed);
                strcat(cmd_args->sol_file, ".sol");
            }
        }
    }

    return 0;
}
