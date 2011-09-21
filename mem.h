

#ifndef _MEM_H
#  define _MEM_H


#include <stddef.h>


void* memdup( const void *base, size_t size );
void* memcalloc( size_t nmemb, size_t size );
void* memmalloc( size_t size );
void* memrealloc( void *ptr, size_t size );


#endif /* _MEM_H */


