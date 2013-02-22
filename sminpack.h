/*
 * ArtTreeKS: Finite dimensional kinematic synthesis for articulated trees
 * Copyright(C) 2010-2012 Edgar Simo-Serra <esimo@iri.upc.edu>
 * License: see synthesis.h
*/


#ifndef _SMINPACK_H
#  define _SMINPACK_H


#include <sys/time.h>

#include <cminpack.h>


int sminpack( cminpack_func_mn fcn, void *p,
      int m, int n, double *x, double *fvec,
      double ftol, double xtol, double gtol, int maxfev,
      double epsfcn, int mode, double factor, int nprint,
      int *nfev, struct timeval *tstart, struct timeval *tend );


#endif /* _SMINPACK_H */

