

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

#include "sminpack.h"


/**
 * @brief Updates the equations stuff.
 */
static int minpack_eqns( void *p, int m, int n,
      const double *x, double *fvec, int iflag )
{
   (void) iflag;
   int i;
   synthesis_t *syn;

   /* Passing data makes this very clean. */
   syn = (synthesis_t*) p;

   /* Paranoia check.s */
   assert( syn->m == m );
   assert( syn->n == n );

   /* First map to synthesis. */
   syn_map_from_x( syn, x, n );
   /* Process branches. */
   for (i=0; i<syn->nbranches; i++)
      syn_calc_branch( syn, &syn->branches[i] );
   /* Map back. */
   syn_map_to_fvec( syn, NULL, NULL, fvec );


   /*
   printf("x:");
   for (i=0; i<n; i++)
      printf("   %+.3e", x[i] );
   printf("\nf:");
   for (i=0; i<m; i++)
      printf("   %+.3e", fvec[i] );
   printf("\n");
   */

   return 0;
}


/**
 * @brief Sets the default minpack options.
 *    @param opts Default options to set.
 */   
void minpack_options_default( minpack_options_t *opts )
{
   assert( opts != NULL );
   memset( opts, 0, sizeof(minpack_options_t) );
   opts->ftol     = 2.22044604926e-16;
   opts->xtol     = 2.22044604926e-16;
   opts->gtol     = 0.;
   opts->maxfev   = INT_MAX;
   opts->epsfcn   = 0.;
   opts->mode     = 1;
   opts->factor   = 100.;
}


/**
 * @brief Gets a human readable string on the termination condition of the MINPACK solver.
 *
 * String must be freed.
 *
 *    @param info Information of the solver.
 *    @return Human readable string to free.
 */
char* minpack_term_string( const minpack_info_t *info )
{
   char buf[256];

   /* Results. */
   switch (info->term_cond) {
      case 0:
         snprintf( buf, sizeof(buf), "Improper input parameters.\n" );
         break;
      case 1:
         snprintf( buf, sizeof(buf), "Completed with ftol=%.3e (actual and predicted sum of squares))\n", info->ftol );
         break;
      case 2:
         snprintf( buf, sizeof(buf), "Completed with xtol=%.3e (relative error between consecutive iterates)\n", info->xtol );
         break;
      case 3:
         snprintf( buf, sizeof(buf), "Completed with ftol=%.3e and xtol=%.3e (sum of squares and consecutive relative error)\n", info->ftol, info->xtol );
         break;
      case 4:
         snprintf( buf, sizeof(buf), "Completed with gtol=%.3e (absolute value cosine between fvec and column jacobian)\n", info->gtol );
         break;
      case 5:
         snprintf( buf, sizeof(buf), "Function calls exceeded (limit %d).\n", info->func_calls );
         break;
      case 6:
         snprintf( buf, sizeof(buf), "Tolerance (ftol=%.3e) too small. No further reduction in the sum of squares is possible.\n", info->ftol );
         break;
      case 7:
         snprintf( buf, sizeof(buf), "Tolerance (xtol=%.3e) too small. No further improvement in the approximate solution x is possible.\n", info->xtol );
         break;
      case 8:
         snprintf( buf, sizeof(buf), "Tolerance (gtol=%.3e) too small. fvec is orthogonal to the columns of the Jacobian to machine precision.\n", info->gtol );
         break;
      default:
         snprintf( buf, sizeof(buf), "Unknown info status!\n" );
         break;
   }

   return strdup(buf);
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
int syn_solve_minpack( synthesis_t *syn, const minpack_options_t *opts, minpack_info_t *info )
{
   minpack_options_t def_opts;
   const minpack_options_t *opts_use;
   int m, n, ret;
   double ftol, xtol, gtol;
   double epsfcn, factor;
   int maxfev, mode, nprint, nfev;
   struct timeval tstart, tend;

   /* Must be finalized. */
   assert( syn->finalized != 0 );

   /* Choose what options to use. */
   minpack_options_default( &def_opts );
   opts_use = (opts!=NULL) ? opts : &def_opts;

   /* Load some stuff. */
   n  = syn->n;
   m  = syn->m;

   /* Set parameters. */
   ftol     = opts_use->ftol;
   xtol     = opts_use->xtol;
   gtol     = opts_use->gtol;
   maxfev   = opts_use->maxfev;
   epsfcn   = opts_use->epsfcn;
   mode     = opts_use->mode;
   factor   = opts_use->factor;
   nprint   = 0;           /* Enables controlled printing of iterates if it is positive. */

   /* Begin the solvation. */
   ret = sminpack( minpack_eqns, (void*)syn, m, n, syn->x, syn->fvec, /* System definition. */
         ftol, xtol, gtol,          /* When convergence has been reached. */
         maxfev,                    /* Cycles. */
         epsfcn, mode, factor,      /* Strange stuff. */
         nprint, &nfev, &tstart, &tend ); /* Output information. */

   /* Store results. */
   if (info != NULL) {
      memset( info, 0, sizeof(minpack_info_t) );
      info->elapsed = (long unsigned)((tend.tv_sec - tstart.tv_sec)
                    + (tend.tv_usec - tstart.tv_usec)/1000000L);
      info->func_calls = nfev;
      info->term_cond = ret;
      info->ftol = ftol;
      info->xtol = xtol;
      info->gtol = gtol;
   }

   return 0;
}


