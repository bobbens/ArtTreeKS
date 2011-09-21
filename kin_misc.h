

#ifndef _KIN_MISC_H
#  define _KIN_MISC_H


#include <dq/dq.h>


/**
 * @brief Representation of a line using plucker coordinates.
 */
typedef struct plucker_s {
   double s[3];   /**< Direction of plucker coordinates. */
   double s0[3];  /**< Moment of the plucker coordinate line. */
} __attribute__((__packed__)) plucker_t;


/*
 * Lie algebra of se(3) stuff.
 */
void lie_joint_mac( plucker_t *dest,
      double val, const plucker_t *src );
void lie_joint_bracket( plucker_t *O,
      const plucker_t *A, const plucker_t *B );


/*
 * Plucker stuff.
 */
void plucker_zero( plucker_t *p );
void plucker_normalize( plucker_t *P );
void plucker_from_dq( plucker_t *p, const dq_t Q );
void plucker_sub( plucker_t *o, const plucker_t *p, const plucker_t *q );


#endif /* _KIN_MISC_H */

