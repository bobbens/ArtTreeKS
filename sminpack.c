
#include "sminpack.h"

#include <stdlib.h>


/**
 * @brief Does minpack without having to fuss about the work buffers.
 *
 * What we do here is automatically allocate a single large buffer which we
 * chop up into smaller ones to reduce the number of memory management calls
 * for running minpack.
 */
int sminpack( minpack_func_mn fcn, void *p,
      int m, int n, double *x, double *fvec,
      double ftol, double xtol, double gtol, int maxfev,
      double epsfcn, int mode, double factor, int nprint,
      int *nfev, struct timeval *tstart, struct timeval *tend )
{
   int ret, ldfjac, *ipvt;
   double *diag, *fjac, *qtf, *wa1, *wa2, *wa3, *wa4;
   char *min_buffer;
   int offset;

   /* Load primary buffer. */
   min_buffer  = calloc( (size_t)(sizeof(double)*(n+(m*n)+n+n+n+n+m)+sizeof(int)*n), 1 );
   offset      = 0;

   /* Set parameters. */
   diag     = (double*)&min_buffer[ offset ]; /* n */
   offset  += n*sizeof(double);
   *nfev    = 0;           /* Output number of iterations. */
   fjac     = (double*)&min_buffer[ offset ]; /* Output array. Upper triangular matrix (m x n) such that p *(jac *jac)*p = r *r. */
   offset  += (m*n)*sizeof(double);
   ldfjac   = m;           /* Leading dimension of fjac. */
   ipvt     = (int*)&min_buffer[ offset ]; /* Output array of dim n such that jac*p = q*r. */
   offset  += n*sizeof(int);
   qtf      = (double*)&min_buffer[ offset ]; /* Output array of dim n with the first 2 elements of (q transpose)*fvec. */
   offset  += n*sizeof(double);
   /* Work arrays. */
   wa1      = (double*)&min_buffer[ offset ];
   offset  += n*sizeof(double);
   wa2      = (double*)&min_buffer[ offset ];
   offset  += n*sizeof(double);
   wa3      = (double*)&min_buffer[ offset ];
   offset  += n*sizeof(double);
   wa4      = (double*)&min_buffer[ offset ];
   /*offset  += m*sizeof(double);*/

   /* Begin the solvation. */
   if (tstart != NULL)
      gettimeofday( tstart, NULL );
   ret = lmdif( fcn, p, m, n, x, fvec, /* System definition. */
         ftol, xtol, gtol,          /* When convergence has been reached. */
         maxfev,                    /* Cycles. */
         epsfcn, diag, mode, factor, /* Strange stuff. */
         nprint, nfev,              /* Output information. */
         fjac, ldfjac, ipvt, qtf,   /* Output informatoin. */
         wa1, wa2, wa3, wa4 );      /* Work arrays. */
   if (tend != NULL)
      gettimeofday( tend, NULL );

   /* Clean up the buffer. */
   free( min_buffer );

   return ret;
}



