/*
 * ArtTreeKS: Finite dimensional kinematic synthesis for articulated trees
 * Copyright(C) 2010-2012 Edgar Simo-Serra <esimo@iri.upc.edu>
 * License: see synthesis.h
*/


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
#include <nlopt.h>


/**
 * @brief Updates the equations stuff.
 */
static double nlopt_eqns( unsigned n, const double *x, double *grad, void *p )
{
   (void) grad; /* We don't actually calculate the gradient ever. */
   int i;
   synthesis_t *syn;
   double ret;

   /* Passing data makes this very clean. */
   syn = (synthesis_t*) p;

   /* Paranoia check.s */
   assert( (unsigned)syn->n == n );

   /* First map to synthesis. */
   syn_map_from_x( syn, x, n );
   /* Process branches. */
   for (i=0; i<syn->nbranches; i++)
      syn_calc_branch( syn, &syn->branches[i] );
   /* Map back. */
   syn_map_to_fvec( syn, NULL, NULL, syn->fvec );

   /* Calculate error. */
   ret = 0.;
   for (i=0; i<syn->m; i++)
      ret += fabs( syn->fvec[i] );
   return ret;
}


/**
 * @brief Sets the default nlopt options.
 *    @param opts Default options to set.
 */   
void nlopt_options_default( nlopt_options_t *opts )
{
   assert( opts != NULL );
   memset( opts, 0, sizeof(nlopt_options_t) );
   opts->algorithm = NLOPT_GN_DIRECT_L;
}


/**
 * @brief Calculates the return string.
 */
static const char* nlopt_retstr( nlopt_result res )
{
   switch (res) {
      /* Sucessful termination. */
      case NLOPT_SUCCESS:           return "Generic success return value.";
      case NLOPT_STOPVAL_REACHED:   return "Optimization stopped because stopval (above) was reached.";
      case NLOPT_FTOL_REACHED:      return "Optimization stopped because ftol_rel or ftol_abs (above) was reached.";
      case NLOPT_XTOL_REACHED:      return "Optimization stopped because xtol_rel or xtol_abs (above) was reached.";
      case NLOPT_MAXEVAL_REACHED:   return "Optimization stopped because maxeval (above) was reached.";
      case NLOPT_MAXTIME_REACHED:   return "Optimization stopped because maxtime (above) was reached.";
      /* Error codes. */
      case NLOPT_FAILURE:           return "Generic failure code";
      case NLOPT_INVALID_ARGS:      return "Invalid arguments (e.g. lower bounds are bigger than upper bounds, an unknown algorithm was specified, etcetera).";
      case NLOPT_OUT_OF_MEMORY:     return "Ran out of memory.";
      case NLOPT_ROUNDOFF_LIMITED:  return "Halted because roundoff errors limited progress. (In this case, the optimization still typically returns a useful result.)";
      case NLOPT_FORCED_STOP:       return "Halted because of a forced termination: the user called nlopt_force_stop(opt) on the optimization’s nlopt_opt object opt from the user’s objective function or constraints.";
      default:                      return "Unknown return code.";
   }
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
int syn_solve_nlopt( synthesis_t *syn, nlopt_options_t *opts, nlopt_info_t *info )
{
   nlopt_result ret;
   nlopt_options_t def_opts;
   const nlopt_options_t *opts_use;
   struct timeval tstart, tend;
   nlopt_opt opt;
   double minf;

   /* Must be finalized. */
   assert( syn->finalized != 0 );

   /* Choose what options to use. */
   nlopt_options_default( &def_opts );
   opts_use = (opts!=NULL) ? opts : &def_opts;

   /* Set up. */
   opt = nlopt_create( opts_use->algorithm, syn->n ); /* algorithm and dimensionality */
   nlopt_set_min_objective( opt, nlopt_eqns, (void*)syn );
   nlopt_set_xtol_rel( opt, 1e-4 );

   /* Set bounds. */
   nlopt_set_lower_bounds( opt, syn->lb );
   nlopt_set_upper_bounds( opt, syn->ub );

   /* Begin the solvation. */
   gettimeofday( &tstart, NULL );
   ret = nlopt_optimize( opt, syn->x, &minf );
   gettimeofday( &tend, NULL );

   /* Results. */
   info->ret_code = ret;
   info->ret_str  = nlopt_retstr( ret );
   info->minf     = minf;
   if (ret < 0) {
      printf( "Error: %s\n", nlopt_retstr( ret ) );
   }
   else {
      printf( "Termination: %s\n", nlopt_retstr( ret ) );
      printf( "   Minimum found at %.3e\n", minf );
   }

   /* Store results. */
   if (info != NULL) {
      memset( info, 0, sizeof(nlopt_info_t) );
      info->elapsed = (long unsigned)((tend.tv_sec - tstart.tv_sec)
                    + (tend.tv_usec - tstart.tv_usec)/1000000L);
   }

   /* Clean up. */
   nlopt_destroy( opt );

   return 0;
}


