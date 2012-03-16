/*
 * ArtTreeKS: Finite dimensional kinematic synthesis for articulated trees
 * Copyright(C) 2010-2012 Edgar Simo-Serra <esimo@iri.upc.edu>
 * License: see synthesis.h
*/

#include "mem.h"

#include <assert.h>
#include <malloc.h>
#include <string.h>


/**
 * @brief Duplicates memory.
 */
void* memdup( const void *base, size_t size )
{
   void *mem;
   if (size == 0)
      return NULL;
   mem = memmalloc( size );
   memcpy( mem, base, size );
   return mem;
}


/**
 * @brief Calloc wrapper.
 */
void* memcalloc( size_t nmemb, size_t size )
{
   void *mem;
   if (size == 0)
      return NULL;
   mem = calloc( nmemb, size );
   assert( mem != NULL );
   return mem;
}


/**
 * @brief Malloc wrapper.
 */
void* memmalloc( size_t size )
{
   void *mem;
   if (size == 0)
      return NULL;
   mem = malloc( size );
   assert( mem != NULL );
   return mem;
}


/**
 * @brief Realloc wrapper.
 */
void* memrealloc( void *ptr, size_t size )
{
   void *mem;
   mem = realloc( ptr, size );
   assert( mem != NULL );
   return mem;
}



