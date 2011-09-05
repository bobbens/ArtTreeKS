

#ifndef MAP_MASK_H
#  define MAP_MASK_H


#include <stddef.h>


/**
 * @brief Structure for converting to and from masked/mapped arrays respectively.
 *
 * We do a tricky thing here. What we do is create a linear array of derivatives and keep a map to
 * the corresponding frame of each index. Then from that we can create the equivalent vector of the
 * frame length with a mask. This allows us to store the data compactly and then change it to masked
 * format to be able to calculate efficiently.
 */
typedef struct mmvec_s {
   /* Information. */
   size_t chunk;
   /* Map. */
   int map_len;
   void *map_vec;
   int *map_map;
   /* Mask. */
   int mask_len;
   void *mask_vec;
   int *mask_mask;
} mm_vec_t;


int mm_initMap( mm_vec_t *mm, size_t chunk, int map_len, void *map_vec, const int *map_map, int mask_len );
int mm_initMask( mm_vec_t *mm, size_t chunk, int mask_len, void *mask_vec, const int *mask_mask );
int mm_initDup( mm_vec_t *mm, const mm_vec_t *in );

void mm_cleanup( mm_vec_t *mm );

int mm_updateMap( mm_vec_t *mm );
int mm_updateMask( mm_vec_t *mm );


#endif /* MAP_MASK_H */

