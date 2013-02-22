/*
 * ArtTreeKS: Finite dimensional kinematic synthesis for articulated trees
 * Copyright(C) 2010-2012 Edgar Simo-Serra <esimo@iri.upc.edu>
 * License: see synthesis.h
*/

#include "kin_misc.h"

#include <dq/vec3.h>

#include <string.h>


/**
 * @brief Multiple and accumulate for Lie algebra elements.
 *    @param[out] dest Destination buffer.
 *    @param[in] val Value to multiply the source by.
 *    @param[in] src Source element.
 */
void lie_joint_mac( plucker_t *dest,
      double val, const plucker_t *src )
{
   dest->s[0]  += val * src->s[0];
   dest->s[1]  += val * src->s[1];
   dest->s[2]  += val * src->s[2];
   dest->s0[0] += val * src->s0[0];
   dest->s0[1] += val * src->s0[1];
   dest->s0[2] += val * src->s0[2];
}


/**
 * @brief Does the Lie bracket operation.
 *    @param[out] O Output Lie algebra element.
 *    @param[in] A Input element 1.
 *    @param[in] B Input element 2.
 */
void lie_joint_bracket( plucker_t *O,
      const plucker_t *A, const plucker_t *B )
{
   plucker_t T;
   double a[3], b[3];

   vec3_cross( T.s,  A->s,    B->s  );
   vec3_cross( a,    A->s,    B->s0 );
   vec3_cross( b,    A->s0,   B->s  );
   vec3_add(   T.s0, a,       b     );

   memcpy( O, &T, sizeof(plucker_t) );
}


/**
 * @brief Zeros a plucker coordinate (or Lie algebra element).
 *    @param p Element to zero.
 */
void plucker_zero( plucker_t *p )
{
   memset( p, 0, sizeof(plucker_t) );
}


/**
 * @brief Normalizes a plucker coordinate.
 */
void plucker_normalize( plucker_t *P )
{
   double c[3];
   vec3_normalize( P->s );
   vec3_cross( c, P->s, P->s0 );
   vec3_cross( P->s0, c, P->s );
}


/**
 * @brief Obtains the core screw element from the dual quaternion.
 *    @param[out] p Screw element obtained from the dual quaternion.
 *    @param[in] Q Dual quaternion to obtain screw element from.
 */
void plucker_from_dq( plucker_t *p, const dq_t Q )
{
   p->s[0]  = Q[1];
   p->s[1]  = Q[2];
   p->s[2]  = Q[3];
   p->s0[0] = Q[4];
   p->s0[1] = Q[5];
   p->s0[2] = Q[6];
}


/**
 * @brief Obtains the screw orientation from the dual quaternion.
 *    @param[out] p Screw element obtained from the dual quaternion.
 *    @param[in] Q Dual quaternion to obtain screw element from.
 */
void plucker_from_dqT( plucker_t *p, const dq_t Q )
{
   p->s[0]  = 0.;
   p->s[1]  = 0.;
   p->s[2]  = 0.;
   p->s0[0] = Q[1];
   p->s0[1] = Q[2];
   p->s0[2] = Q[3];
}


/**
 * @brief Subtracts a plucker coordinate from another.
 *    @param[out] o Output plucker coordinate.
 *    @param[in] p Input plucker coordinate to subtract from.
 *    @param[in] q Input plucker coordinate to subtract.
 */
void plucker_sub( plucker_t *o, const plucker_t *p, const plucker_t *q )
{
   vec3_sub( o->s,  p->s,  q->s );
   vec3_sub( o->s0, p->s0, q->s0 );
}



