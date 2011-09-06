

#include "mapmask.h"

#include <string.h>
#include <malloc.h>


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
   int i;
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

int main (void)
{
   double *d, map_vec[6] = { 1., 2., 3., 5., 4., 6. };
   int map_map[6] = { 1, 3, 5, 9, 7, 11 };
   mm_vec_t mm;

   mm_initMap( &mm, sizeof(double), 6, map_vec, map_map, 12 );
   print_mm( &mm );

   printf( "TEST MAP\n" );
   d = (double*) mm.map_vec;
   d[0] = 6.;
   d[5] = 1.;
   print_mm( &mm );
   mm_updateMask( &mm );
   print_mm( &mm );

   printf( "TEST MASK\n" );
   d = (double*) mm.mask_vec;
   d[1]  = 1.;
   d[11] = 6.;
   print_mm( &mm );
   mm_updateMap( &mm );
   print_mm( &mm );

   mm_cleanup( &mm );
   return 0;
}
#endif






