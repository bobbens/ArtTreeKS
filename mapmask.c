/*
 * ArtTreeKS: Finite dimensional kinematic synthesis for articulated trees
 * Copyright(C) 2010-2012 Edgar Simo-Serra <esimo@iri.upc.edu>
 * License: see synthesis.h
*/


#include "mapmask.h"

#include <string.h>
#include <malloc.h>
#include <stdlib.h>


static void* memdup( const void *base, size_t size )
{
   void *mem;
   if (size == 0)
      return NULL;
   mem = malloc( size );
   memcpy( mem, base, size );
   return mem;
}


int mm_initMap( mm_vec_t *mm, size_t chunk, int map_len, void *map_vec, const int *map_map, int mask_len )
{
   int i;
   /* Duplicate the new data. */
   mm->chunk   = chunk;
   mm->map_len = map_len;
   mm->map_vec = memdup( map_vec, mm->map_len * mm->chunk );
   mm->map_map = memdup( map_map, mm->map_len * sizeof(int) );
   /* Create the mask and map the vector. */
   mm->mask_len  = mask_len;
   mm->mask_vec  = calloc( mm->mask_len, mm->chunk );
   mm->mask_mask = calloc( mm->mask_len, sizeof(int) );
   for (i=0; i<mm->map_len; i++)
      mm->mask_mask[ mm->map_map[i] ] = 1;
   /* Synchronize data. */
   mm_updateMask( mm );
   return 0;
}

int mm_initMask( mm_vec_t *mm, size_t chunk, int mask_len, void *mask_vec, const int *mask_mask )
{
   int i, j;
   /* Duplicate new data. */
   mm->chunk    = chunk;
   mm->mask_len = mask_len;
   mm->mask_vec = memdup( mask_vec, mm->mask_len * mm->chunk );
   if (mask_mask == NULL) {
      mm->mask_mask = malloc( mm->mask_len * sizeof(int) );
      for (i=0; i<mm->mask_len; i++)
         ((int*)mm->mask_mask)[i] = 1;
   }
   else
      mm->mask_mask = memdup( mask_mask, mm->mask_len * sizeof(int) );
   /* Create the map. */
   mm->map_len = 0;
   for (i=0; i<mm->mask_len; i++)
      if (mm->mask_mask[i])
         mm->map_len++;
   mm->map_vec = calloc( mm->map_len, mm->chunk );
   mm->map_map = calloc( mm->map_len, sizeof(int) );
   j = 0;
   for (i=0; i<mm->mask_len; i++)
      if (mm->mask_mask[i])
         mm->map_map[ j++ ] = i;
   /* Synchronize data. */
   mm_updateMap( mm );
   return 0;
}

int mm_initDup( mm_vec_t *mm, const mm_vec_t *in )
{
	mm->chunk = in->chunk;
	mm->mask_len = in->mask_len;
	mm->mask_vec = memdup( in->mask_vec, in->chunk * in->mask_len );
	mm->mask_mask = memdup( in->mask_mask, in->mask_len * sizeof(int) );
	mm->map_len = in->map_len;
	mm->map_vec = memdup( in->map_vec, in->chunk * in->map_len );
	mm->map_map = memdup( in->map_map, in->map_len * sizeof(int) );
   return 0;
}

void mm_cleanup( mm_vec_t *mm )
{
   free( mm->map_vec );
   free( mm->map_map );
   free( mm->mask_vec );
   free( mm->mask_mask );
}

int mm_updateMap( mm_vec_t *mm )
{
   int i;
   if (mm->chunk == 0)
      return 1;
   for (i=0; i<mm->map_len; i++)
      memcpy( mm->map_vec + mm->chunk*i, mm->mask_vec + mm->chunk*mm->map_map[i], mm->chunk );
   return 0;
}

int mm_updateMask( mm_vec_t *mm )
{
   int i;
   if (mm->chunk == 0)
      return 1;
   for (i=0; i<mm->map_len; i++)
      memcpy( mm->mask_vec + mm->chunk*mm->map_map[i], mm->map_vec + mm->chunk*i, mm->chunk );
   return 0;
}


#if 0

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

static int print_mm( mm_vec_t *mm )
{
   int i;
   double *d, *t;
   d = calloc( mm->chunk, mm->mask_len );
   t = (double*) mm->map_vec;
   for (i=0; i<mm->map_len; i++)
      memcpy( &d[ mm->map_map[i] ], mm->map_vec + mm->chunk*i, sizeof(double) );

   printf( "map:  " );
   for (i=0; i<mm->mask_len; i++)
      printf( "%.0f, ", d[i] );
   printf( "\n" );

   memset( d, 0, mm->chunk * mm->mask_len );
   t = (double*) mm->mask_vec;
   for (i=0; i<mm->mask_len; i++)
      if (mm->mask_mask[i])
         memcpy( &d[i], mm->mask_vec + mm->chunk*i, sizeof(double) );

   printf( "mask: " );
   for (i=0; i<mm->mask_len; i++)
      printf( "%.0f, ", d[i] );
   printf( "\n" );

   free(d);
   return 0;
}


static const char* mm_compare( const mm_vec_t *mma, const mm_vec_t *mmb )
{
   if (mma->chunk    != mmb->chunk)
      return "Chunk mismatch";
   if (mma->map_len  != mmb->map_len)
      return "Map len mismatch";
   if (mma->mask_len != mmb->mask_len)
      return "Mask len mismatch";
   /* Map. */
   if (memcmp( mma->map_vec,   mmb->map_vec,   mma->chunk * mma->map_len ))
      return "Map vec mismatch";
   if (memcmp( mma->map_map,   mmb->map_map,   sizeof(int) * mma->map_len ))
      return "Map map misamtch";
   /* Mask. */
   if (memcmp( mma->mask_vec,  mmb->mask_vec,  mma->chunk * mma->mask_len ))
      return "Mask vec mismatch";
   if (memcmp( mma->mask_mask, mmb->mask_mask, sizeof(int) * mma->mask_len ))
      return "Mask mask mismatch";
   return NULL;
}


int main (void)
{
   double *d, map_vec[6] = { 1., 2., 3., 4., 5., 6. };
   int map_map[6]        = { 1,  3,  5,  7,  9, 11  };
   double mask_vec[12]   = { 0., 1., 0., 2., 0., 3., 0., 4., 0., 5., 0., 6. };
   int mask_mask[12]     = { 0,  1,  0,  1,  0,  1,  0,  1,  0,  1,  0,  1  };
   mm_vec_t mma, mmb;
   const char *ret;

   /* Test load. */
   mm_initMap( &mma, sizeof(double), 6, map_vec, map_map, 12 );
   mm_initMask( &mmb, sizeof(double), 12, mask_vec, mask_mask );
   ret = mm_compare( &mma, &mmb );
   if (ret != NULL) {
      printf( "%s\n", ret );
      printf( "MMA\n" );
      print_mm( &mma );
      printf( "MMB\n" );
      print_mm( &mmb );
   }

   /* Test manipulate map. */
   d = (double*) mma.map_vec;
   d[0] = 6.;
   d[5] = 1.;
   mm_updateMask( &mma );
   d = (double*) mmb.mask_vec;
   d[1]  = 6.;
   d[11] = 1.;
   mm_updateMap( &mmb );
   ret = mm_compare( &mma, &mmb );
   if (ret != NULL) {
      printf( "%s\n", ret );
      printf( "MMA\n" );
      print_mm( &mma );
      printf( "MMB\n" );
      print_mm( &mmb );
   }

   mm_cleanup( &mma );
   mm_cleanup( &mmb );
   return 0;
}

#endif




