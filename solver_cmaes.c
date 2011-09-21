

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


/**
 * @brief Sets the default cmaes options.
 *    @param opts Default options to set.
 */   
void cmaes_options_default( cmaes_options_t *opts )
{
   assert( opts != NULL );
   memset( opts, 0, sizeof(cmaes_options_t) );
   opts->lambda = 100;
}


/**
 * @brief Normalizes a plucker coordinate.
 */
static void plucker_normalize( plucker_t *P )
{
   double c[3];
   vec3_normalize( P->s );
   vec3_cross( c, P->s, P->s0 );
   vec3_cross( P->s0, c, P->s );
}


/**
 * @brief Normalize synthesis.
 */
static void cmaes_normalize( synthesis_t *syn )
{
   int i;
   kin_joint_t *kj;
   for (i=0; i<syn->njoints; i++) {
      kj = syn->joints[i];
      if (kj->claim_S != NULL)
         plucker_normalize( &kj->S );
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
   double fit, *fitvals, *stddev;
   const double *xfinal;
   double *const *pop;
   struct timeval tstart, tend;
   const char *done;

   /* Parameters. */
   dim      = syn->n; /**< Dimension of the system. */
   lambda   = opts->lambda;

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

   /* Set stop fitness. */
   evo.sp.stopTolFun        = 1e-10;
   //evo.sp.stopMaxFunEvals   = ;
   //evo.sp.stopMaxIter       = ;

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
         cmaes_normalize( syn );

         /* Calculate and store fitness. */
         fitvals[p] = cmaes_fitness( syn );
      }

      /* Update the distribution. */
      cmaes_UpdateDistribution( &evo, fitvals );

      /* Output some stuff if necessary. */
      iter++;
      printf( "[%d] Best: %.3e\n", iter, cmaes_Get( &evo, "fbestever" ) );
   }
   gettimeofday( &tend, NULL );
   if (info != NULL)
      info->elapsed = (unsigned long)((tend.tv_sec - tstart.tv_sec)
                    + (tend.tv_usec - tstart.tv_usec)/1000000);
   printf("%s\n",done);

   /* Map. */
   xfinal = cmaes_GetPtr( &evo, "xbestever" );
   syn_map_from_x(  syn, xfinal,     syn->n );
   cmaes_normalize( syn );
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




