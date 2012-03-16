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

#include <minpack.h>


static synthesis_t *global_syn;


/**
 * @brief Updates the equations stuff.
 */
static void minpack_eqns( int *m, int *n,
      double *x, double *fvec, int *iflag )
{
   (void) iflag;
   int i;
   synthesis_t *syn;

   /* Passing data makes this very clean. */
   syn = global_syn;

   /* Paranoia check.s */
   assert( syn->m == *m );
   assert( syn->n == *n );

   /* First map to synthesis. */
   syn_map_from_x( syn, x, *n );
   /* Process branches. */
   for (i=0; i<syn->nbranches; i++)
      syn_calc_branch( syn, &syn->branches[i] );
   /* Map back. */
   syn_map_to_fvec( syn, NULL, NULL, fvec );
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
int syn_solve_fminpack( synthesis_t *syn, const minpack_options_t *opts, minpack_info_t *info )
{
   minpack_options_t def_opts;
   const minpack_options_t *opts_use;
   int m, n, ret;
   double ftol, xtol, gtol;
   double epsfcn, *diag, factor;
   int maxfev, mode, nprint, nfev, ldfjac, *ipvt;
   double *fjac, *qtf, *wa1, *wa2, *wa3, *wa4;
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
   global_syn = syn;
   ftol     = opts_use->ftol;
   xtol     = opts_use->xtol;
   gtol     = opts_use->gtol;
   maxfev   = opts_use->maxfev;
   epsfcn   = opts_use->epsfcn;
   diag     = calloc( (size_t)n, sizeof(double) ); /* Input array of length n. */
   mode     = opts_use->mode;
   factor   = opts_use->factor;
   nprint   = 0;           /* Enables controlled printing of iterates if it is positive. */
   nfev     = 0;           /* Output number of iterations. */
   fjac     = calloc( (size_t)(m*n), sizeof(double) ); /* Output array. Upper triangular matrix (m x n) such that p *(jac *jac)*p = r *r. */
   ldfjac   = m;           /* Leading dimension of fjac. */
   ipvt     = calloc( (size_t)n, sizeof(int) ); /* Output array of dim n such that jac*p = q*r. */
   qtf      = calloc( (size_t)n, sizeof(double) ); /* Output array of dim n with the first 2 elements of (q transpose)*fvec. */
   /* Work arrays. */
   wa1      = calloc( (size_t)n, sizeof(double) );
   wa2      = calloc( (size_t)n, sizeof(double) );
   wa3      = calloc( (size_t)n, sizeof(double) );
   wa4      = calloc( (size_t)m, sizeof(double) );

   /* Begin the solvation. */
   gettimeofday( &tstart, NULL );
   lmdif_( minpack_eqns, &m, &n, syn->x, syn->fvec, /* System definition. */
         &ftol, &xtol, &gtol,       /* When convergence has been reached. */
         &maxfev,                   /* Cycles. */
         &epsfcn, diag, &mode, &factor, /* Strange stuff. */
         &nprint, &ret, &nfev,      /* Output information. */
         fjac, &ldfjac, ipvt, qtf,  /* Output informatoin. */
         wa1, wa2, wa3, wa4 );      /* Work arrays. */
   gettimeofday( &tend, NULL );

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

   /* Clean up work arrays. */
   free( wa1 );
   free( wa2 );
   free( wa3 );
   free( wa4 );

   /* Clean up results and data. */
   free( diag );
   free( fjac );
   free( ipvt );
   free( qtf  );

   return 0;
}


