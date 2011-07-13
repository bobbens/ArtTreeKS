


#ifndef _SOLVER_H
#  define _SOLVER_H


#include "synthesis.h"


/**
 * @brief Options for the minpack solver.
 */
typedef struct minpack_options_s {
   double ftol;      /**< F tolerance. */
   double xtol;      /**< X tolerance. */
   double gtol;      /**< G tolerance. */
   int maxfev;       /**< Maximum function calls. */
   double epsfcn;    /**< Maximum forward step. */
   int    mode;      /**< No idea. */
   double factor;    /**< No idea. */
} minpack_options_t;


/**
 * @brief Minpack information.
 */
typedef struct minpack_info_s {
   unsigned long elapsed;  /**< Elapsed time in seconds. */
   int func_calls;         /**< Number of function calls. */
   int term_cond;          /**< Termination condition. */
   double ftol;            /**< F tolerance. */
   double xtol;            /**< X tolerance. */
   double gtol;            /**< G tolerance. */
} minpack_info_t;


/*
 * Minpack functions.
 */
void minpack_options_default( minpack_options_t *opts );
char* minpack_term_string( const minpack_info_t *info );
int syn_solve_minpack( synthesis_t *syn, const minpack_options_t *opts, minpack_info_t *info );
int syn_solve_fminpack( synthesis_t *syn, const minpack_options_t *opts, minpack_info_t *info );


/*
 * Forward declarations.
 */
struct ga_info_s;
typedef struct ga_info_s ga_info_t;
struct ga_options_s;
typedef struct ga_options_s ga_options_t;


/*
 * Genetic algorithm function types.
 */
typedef int (*ga_epoch_t)( ga_options_t *opts, const ga_info_t *info ); /**< Epoch function. */


/**
 * @brief Genetic algorithm options.
 */
struct ga_options_s {
   int verbose;         /**< Verbosity. */
   int threads;         /**< Threads to use. */
   int converge;        /**< Use local minimizer. */
   /* Parameters of search area. */
   /* Genetic algorithm parameters. */
   int population;      /**< Population size. */
   int generations;     /**< Generations to run. */
   double eliteness;    /**< Eliteness factor [0-1]. % Entities that passes over without change
                             directly at the start of the generation, these are not mutated.*/
   double crossover;    /**< Crossover factor [0-1]. Percent of new entities created by crossover. */
   double mutation;     /**< Mutation factor [0-1]. Percent change of mutating each entity. */
   double seed_mul;     /**< Multiplier that defines how much larger the population should be when seeding. */
   /* Stop conditions. */
   int stop_sigusr1;    /**< Stop on SIGUSR1. */
   int stop_sigint;     /**< Stop on SIGINT. */
   unsigned long stop_elapsed; /**< Stop on elapsed time (0 is infinite). */
   double stop_fitness; /**< Fitness to stop at. */
   int sigfpe;          /**< Enable floating point exceptions. */
   /* MINPACK options. */
   minpack_options_t minpack; /**< MINPACK options to use. */
   /* Callback functions. */
   ga_epoch_t f_epoch;  /**< Epoch function pointer. */
};


/**
 * @brief genetic algorithm information.
 */
struct ga_info_s {
   unsigned long elapsed;  /**< Elapsed time in seconds. */
   int generations;        /**< Generations run. */
   int max_generations;    /**< Maximum generations to run. */
   /* Statistics. */
   double fit_best;        /**< Best fitness. */
   double fit_mean;        /**< Mean fitness. */
   double fit_stddev;      /**< Standard deviation fitness. */
};


/*
 * Genetic algorithm functions.
 */
void ga_options_default( ga_options_t *opts );
int syn_solve_ga( synthesis_t *syn, ga_options_t *opts, ga_info_t *info );


typedef struct nlopt_options_s {
   int algorithm;       /**< Algorithm to use. */
} nlopt_options_t;


typedef struct nlopt_info_s {
   int ret_code;           /**< Result code. */
   const char *ret_str;    /**< Result string. */
   double minf;            /**< Minimum fitness found. */
   unsigned long elapsed;  /**< Elapsed time in seconds. */
} nlopt_info_t;


void nlopt_options_default( nlopt_options_t *opts );
int syn_solve_nlopt( synthesis_t *syn, nlopt_options_t *opts, nlopt_info_t *info );


#endif /* _SOLVER_H */


