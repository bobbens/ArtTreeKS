

#include "solver.h"

#include <ctype.h>
#include <stdio.h>
#include <assert.h>
#include <getopt.h>
#include <limits.h>
#include <malloc.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/time.h>
#include <math.h>

#include <dq/vec3.h>

#include "cmaes_interface.h"
#include "kin_misc.h"


/**
 * @brief Sets the default cmaes options.
 *    @param opts Default options to set.
 */   
void cmaes_options_default( cmaes_options_t *opts )
{
   assert( opts != NULL );
   memset( opts, 0, sizeof(cmaes_options_t) );
   opts->converge       = 1; /* Use convergence by default. */
   opts->lambda         = 0; /* 0 defaults to optimal default value defined by the paper. */
   /* Stop conditions. */
   opts->stop_fitness   = 1e-10;
   opts->stop_evals     = UINT_MAX;
   opts->stop_iter      = UINT_MAX;
   /* Minpack. */
   minpack_options_default( &opts->minpack );
}


/**
 * @brief Normalize synthesis.
 */
static void cmaes_normalize( synthesis_t *syn, cmaes_options_t *opts )
{
#if 0
   int i;
   kin_joint_t *kj;
   for (i=0; i<syn->njoints; i++) {
      kj = syn->joints[i];
      if (kj->claim_S != NULL)
         plucker_normalize( &kj->S );
   }
#endif
   if (opts->converge) {
      syn_map_to_x(      syn, NULL, NULL, syn->x );
      syn_solve_minpack( syn, NULL, NULL );
   }
}


/**
 * @brief Calculcates fitness, does not modify anything other than fvec.
 */
static double cmaes_fitness( synthesis_t *syn )
{
   int i;
   double fit;

   syn_calc( syn );

   fit = 0.;
   for (i=0; i<syn->m; i++)
      fit += fabs( syn->fvec[i] );

   return fit;
}


/**
 * @brief CMA-ES Solver.
 *
 * @note thread safe, uses no global variables.
 *
 *    @param syn Synthesis object.
 *    @param[in] opts Options (NULL for default).
 *    @param[out] Results of the solver.
 *    @return 0 on success.
 */
int syn_solve_cmaes( synthesis_t *syn, cmaes_options_t *opts, cmaes_info_t *info )
{
   cmaes_t evo;
   unsigned int iter;
   int p, i, dim, lambda, npop;
   double fit, *fitvals, *stddev, best;
   const double *xfinal;
   double *const *pop;
   struct timeval tstart, tend;
   const char *done;
   cmaes_options_t *opts_use, opts_def;

   cmaes_options_default( &opts_def );
   opts_use = (opts != NULL) ? opts : &opts_def;

   /* Parameters. */
   dim      = syn->n; /**< Dimension of the system. */
   lambda   = opts_use->lambda;
   if (lambda == 0)
      lambda = 4 + floor( 3*log(syn->n) );

   /* Map to use as initial position. */
   syn_map_to_x( syn, NULL, NULL, syn->x );

   /* Standard deviation. */
   stddev = calloc( syn->n, sizeof(double) );
   for (i=0; i<syn->n; i++) {
      syn->x[i] = 0.5;
      stddev[i] = 0.5; /* Todo smarter initialization. */
   }

   /* Initialize stuff. */
   fitvals = cmaes_init( &evo,
         dim,     /* Dimension. */
         syn->x,  /* X starting position. */
         stddev,  /* X starting standard deviation. */
         0,       /* Starting seed? */
         lambda,  /* Starting population. */
         "non"    /* Read from file. */
         );
   printf( "%s\n", cmaes_SayHello(&evo) );
   //cmaes_ReadSignals( &evo, "signals.par" );

   /* Settings. */

   /* Set stop fitness. */
   evo.sp.stopTolFun        = opts_use->stop_fitness;
   evo.sp.stopMaxFunEvals   = (double) opts_use->stop_evals;
   evo.sp.stopMaxIter       = (double) opts_use->stop_iter;

   /* Start iterating fool! */
   gettimeofday( &tstart, NULL );
   iter = 0;
   while ((done = cmaes_TestForTermination( &evo )) == NULL) {
      pop   = cmaes_SamplePopulation( &evo );
      npop  = cmaes_Get( &evo, "popsize" );

      /* Here we must analyze each member in the population. */
      for (p=0; p<npop; p++) {

         /* First map to synthesis. */
         syn_map_from_x( syn, pop[p], syn->n );

         /* Reenforce plucker coordinates here */
         cmaes_normalize( syn, opts_use );

         /* Calculate and store fitness. */
         fitvals[p] = cmaes_fitness( syn );
      }

      /* Update the distribution. */
      cmaes_UpdateDistribution( &evo, fitvals );

      /* Output some stuff if necessary. */
      iter++;
      best = cmaes_Get( &evo, "fbestever" );
      printf( "[%d] Best: %.3e\n", iter, best );
      if (best < opts_use->stop_fitness)
         break;
   }
   gettimeofday( &tend, NULL );
   if (info != NULL)
      info->elapsed = (unsigned long)((tend.tv_sec - tstart.tv_sec)
                    + (tend.tv_usec - tstart.tv_usec)/1000000);
   if (done != NULL)
      printf("%s\n",done);

   /* Map. */
   xfinal = cmaes_GetPtr( &evo, "xbestever" );
   syn_map_from_x(  syn, xfinal,     syn->n );
   cmaes_normalize( syn, opts_use );
   syn_map_to_x(    syn, NULL, NULL, syn->x );
   fit = cmaes_fitness( syn );
   if (info != NULL) {
      info->minf = fit;
      info->iterations = iter;
   }

   /* Clean up. */
   cmaes_exit( &evo );
   //free( xfinal );

   return 0;
}




