

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
 * @brief Does local convergence on a synthesis object.
 *
 * Uses Levenberg-Marquadt algorithm.
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
   int iter, p, i, dim, lambda, npop;
   double fit, *fitvals, *stddev, *xfinal;
   double *const *pop;
   kin_joint_t *kj;

   /* Parameters. */
   dim      = syn->n; /**< Dimension of the system. */
   lambda   = 250;

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

   /* Start iterating fool! */
   iter = 0;
   while (!cmaes_TestForTermination( &evo )) {
      pop   = cmaes_SamplePopulation( &evo );
      npop  = cmaes_Get( &evo, "popsize" );

      /* Here we must analyze each member in the population. */
      for (p=0; p<npop; p++) {

         /* First map to synthesis. */
         syn_map_from_x( syn, pop[p], syn->n );

         /* Reenforce plucker coordinates here */
         for (i=0; i<syn->njoints; i++) {
            kj = syn->joints[i];
            if (kj->claim_S != NULL)
               plucker_normalize( &kj->S );
         }

         /* Process branches. */
         for (i=0; i<syn->nbranches; i++)
            syn_calc_branch( syn, &syn->branches[i] );

         /* Map back. */
         syn_map_to_fvec( syn, NULL, NULL, syn->fvec );

         /* Calculate error. */
         fit = 0.;
         for (i=0; i<syn->m; i++)
            fit += fabs( syn->fvec[i] );

         /* Store fitness value. */
         fitvals[p] = fit;
      }

      /* Update the distribution. */
      cmaes_UpdateDistribution( &evo, fitvals );

      iter++;
      printf( "[%d] Best: %.3e\n", iter, cmaes_Get( &evo, "fbestever" ) );
   }

   /* Map. */
   xfinal = cmaes_GetNew( &evo, "xbest" );
   syn_map_from_x(   syn, xfinal, syn->n );
   syn_map_to_x(     syn, NULL, NULL, syn->x );

   /* Clean up. */
   cmaes_exit( &evo );
   free( xfinal );

   return 0;
}




