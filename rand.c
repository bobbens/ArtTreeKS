/*
 * ArtTreeKS: Finite dimensional kinematic synthesis for articulated trees
 * Copyright(C) 2010-2012 Edgar Simo-Serra <esimo@iri.upc.edu>
 * License: see synthesis.h
*/


#include "rand.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>


/*
 * Globals.
 */
static int rand_fd   = -1; /**< /dev/urandom */


/**
 * @brief Initializes random subsystem.
 *
 * @TODO Possibly make it thread safe or use a different approach.
 */
void rand_init (void)
{
   if (rand_fd < 0)
      rand_fd = open( "/dev/urandom", O_RDONLY );
   assert( rand_fd >= 0 );

   unsigned int seed;
   ssize_t ret = read( rand_fd, &seed, sizeof(seed) );
   assert( ret == sizeof(seed) );
   srand( seed );
}
/**
 * @brief Exits random subsystem.
 */
void rand_exit (void)
{
   if (rand_fd >= 0)
      close( rand_fd );
   rand_fd = -1;
}
/**
 * @brief Gets a random boolean.
 */
int rand_bool (void)
{
   return (rand_double() > 0.5);
}
/**
 * @brief Gets an int in range [low,high].
 */
int rand_int_range( int low, int high )
{
   return (int)round( rand_double_range( (double)low, (double)high ) );
}
/**
 * @brief Gets a random double in the [0,1] range.
 */
double rand_double (void)
{
#ifdef RAND_ACCURATE
   ssize_t ret;
   uint64_t r;
   ret = read( rand_fd, &r, sizeof(r) );
   assert( ret == sizeof(r) );
   return (double)r / (double)UINT64_MAX;
#else /* RAND_ACCURATE */
   return (double)rand() / (double)RAND_MAX;
#endif /* RAND_ACCURATE */
}
/**
 * @brief Gets a random double in the range [low,high].
 */
double rand_double_range( double low, double high )
{
   return (low + (high-low)*rand_double());
}
/**
 * @brief Gets a random double following the exponential function for a given lambda.
 */
double rand_double_exponential( double lambda )
{
   return (-1./lambda) * log(rand_double());
}
/**
 * @brief Gets a random double that follows the gaussian distribution.
 */
double rand_double_gaussian (void)
{
   double theta, rsq;
   theta = rand_double_range( 0., 2.*M_PI );
   rsq   = rand_double_exponential( 0.5 );
   return (sqrt(rsq) * cos(theta));
}

