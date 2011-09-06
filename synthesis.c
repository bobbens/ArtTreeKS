

#include "synthesis.h"

#include <math.h>
#include <errno.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <sys/stat.h>

#include <dq/vec3.h>


#define PF     "%.18le"


/*
 * Prototypes.
 */
/* Branches. */
static int syn_construct_branches_count( synthesis_t *syn, kin_branch_iter_t *branch, void *data );
static int syn_construct_branches_walk( synthesis_t *syn, kin_branch_iter_t *branch, void *data );
static int syn_branch_iter_walk( synthesis_t *syn,
      int (*func)(synthesis_t*, kin_branch_iter_t*, void*),
      void *data, kin_branch_iter_t *head, kin_branch_iter_t *last );
/* Joints. */
static void kin_joint_data_free( kin_joint_data_t *data );
static void kin_joint_data_dup( kin_joint_data_t *dest, const kin_joint_data_t *src );
static int kin_joint_data_claim( synthesis_t *syn, kin_joint_data_t *data, const char *name );
/* Saving. */
static int kin_obj_save( const kin_object_t *obj, FILE *stream, int depth, const char *sym );
/* Finalization. */
static int kin_obj_fin( synthesis_t *syn, kin_object_t *obj );
static int syn_joint_add( synthesis_t *parent, kin_joint_t *joint );
static int syn_tcp_add( synthesis_t *parent, kin_object_t *tcp );
/* Claims. */
static kin_claim_t* syn_claim_add( kin_claim_t *claim,
      size_t len, size_t indep, double *v, double *lb, double *ub, const char *name );
static void syn_claim_destroy( kin_claim_t *claim );
static int syn_map_claim_to_vec( const kin_claim_t *par, int *n, int *ni,
      double *v, double *lb, double *ub );
static kin_claim_t* syn_claim_x( synthesis_t *syn, size_t len, size_t indep,
      double *v, double *lb, double *ub, const char *name );
static kin_claim_t* syn_claim_fvec( synthesis_t *parent, size_t len, size_t indep,
      double *v, const char *name );
static int syn_set_bounds( synthesis_t *syn );
/* Function voodoo. */
static int kin_obj_chain_free( kin_object_t *obj );
static int kin_obj_chain_dup( const kin_object_t *obj, kin_object_t *newobj );
static int kin_obj_chain_fin( kin_object_t *obj, synthesis_t *syn );
static int kin_obj_chain_save( const kin_object_t *obj, FILE *stream, const char *self );
static int kin_obj_tcp_free( kin_object_t *obj );
static int kin_obj_tcp_dup( const kin_object_t *obj, kin_object_t *newobj );
static int kin_obj_tcp_fin( kin_object_t *obj, synthesis_t *syn );
static int kin_obj_tcp_save( const kin_object_t *obj, FILE *stream, const char *self );
static int kin_obj_split_free( kin_object_t *obj );
static int kin_obj_split_dup( const kin_object_t *obj, kin_object_t *newobj );
static int kin_obj_split_fin( kin_object_t *obj, synthesis_t *syn );
static int kin_obj_split_save( const kin_object_t *obj, FILE *stream, const char *self );
/* Helper functions. */
static void* memdup( const void *base, size_t size );
static void* memcalloc( size_t nmemb, size_t size );
static void* memmalloc( size_t size );
static void* memrealloc( void *ptr, size_t size );


/**
 * @brief Frees a branch.
 *    @param branch Branch to free.
 */
static void syn_branch_free( kin_branch_t *branch )
{
   free( branch->joints );
}
/**
 * @brief Small wrapper to count branches.
 *    @param syn Unused.
 *    @param branch Unused.
 *    @param data Unused.
 */
static int syn_construct_branches_count( synthesis_t *syn, kin_branch_iter_t *branch, void *data )
{
   (void) syn;
   (void) branch;
   (void) data;
   return 0;
}
/**
 * @brief Walks through branches to construct the branch data.
 *    @param syn Synthesis object being walked through.
 *    @param branch Linked list of branch path.
 *    @param data Casted int* representing number of branch iterating through.
 */
static int syn_construct_branches_walk( synthesis_t *syn, kin_branch_iter_t *branch, void *data )
{
   int *i, j, k, n;
   kin_branch_iter_t *l;
   kin_branch_t *b;
   kin_object_t *tcp;
   i = (int*)data;

   /* Count joints. */
   n = 0;
   tcp = NULL;
   for (l=branch; l!=NULL; l=l->next) {
      if (l->obj->type == KIN_TYPE_CHAIN)
         n += l->obj->d.chain.njoints;
      else if (l->obj->type == KIN_TYPE_TCP)
         tcp = l->obj;
   }
   assert( tcp != NULL );

   /* Set pointers. */
   b         = &syn->branches[*i];
   b->joints = memcalloc( n, sizeof(kin_joint_t*) );
   b->njoints = n;
   b->tcp    = tcp;
   j = 0;
   for (l=branch; l!=NULL; l=l->next)
      if (l->obj->type == KIN_TYPE_CHAIN)
         for (k=0; k<l->obj->d.chain.njoints; k++)
            b->joints[j++] = &l->obj->d.chain.joints[k];

   /* Sanity checking. */
   if ((b->tcp->d.tcp.claim_acc != NULL) &&
         (b->tcp->d.tcp.claim_vel == NULL))
      assert( "Velocity data can't be NULL if acceleration data is not NULL." );
   for (k=0; k < b->tcp->d.tcp.nP; k++) {
      if (((b->tcp->d.tcp.A.mask_mask == NULL) ||
               (b->tcp->d.tcp.A.mask_mask[k] != 0)) &&
            ((b->tcp->d.tcp.V.mask_mask != NULL) &&
               (b->tcp->d.tcp.V.mask_mask[k] == 0))) {
         assert( "Velocity mask must be a superset of acceleration mask." );
      }
   }

   /* Increment branch. */
   (*i)++;
   return 0;
}
/**
 * @brief Constructs the branches for a synthesis object.
 *    @param syn Synthesis object to create branches for.
 */
static int syn_construct_branches( synthesis_t *syn )
{
   int i, n;

   assert( syn->branches == NULL );
   assert( !syn->finalized );

   /* Count branches. */
   n = syn_branch_iter( syn, syn_construct_branches_count, NULL );
   assert( n > 0 );

   /* Allocate branches. */
   syn->nbranches = n;
   syn->branches  = memcalloc( syn->nbranches, sizeof(kin_branch_t) );

   i = 0;
   syn_branch_iter( syn, syn_construct_branches_walk, &i );
   return 0;
}

/**
 * @brief Walks over the syntehsis iteration branches.
 *    @param syn Synthesis object to walk over branches.
 *    @param func Function to execute at the end of each branch.
 *    @param head Head of tree-like topology.
 *    @param last Last element being processed.
 *    @return Number of branches found.
 */   
static int syn_branch_iter_walk( synthesis_t *syn,
      int (*func)(synthesis_t*, kin_branch_iter_t*, void*),
      void *data, kin_branch_iter_t *head, kin_branch_iter_t *last )
{
   int i, ret;
   kin_branch_iter_t *tail;
   kin_object_t *o;

   /* Here we look at the last object and continue processing.
    * The goal is to find all the TCP and run the function on them. */
   ret   = 0;
   o     = last->obj;
   switch (o->type) {
      case KIN_TYPE_CHAIN:
         /* If the chain is not a dead end (something hanging on the end)
          * we shal continue iterating past it. It's just a joining element. */
         if (o->next != NULL) {
            tail        = memcalloc( 1, sizeof(kin_branch_iter_t) );
            last->next  = tail;
            tail->obj   = o->next;
            ret         = syn_branch_iter_walk( syn, func, data, head, tail );
            last->next  = NULL;
            free( tail );
         }
         break;

      case KIN_TYPE_SPLITTER:
         /* Splitter has multiple branches leading off so we must iterate over
          * each one exactly like the chain. */
         tail        = memcalloc( 1, sizeof(kin_branch_iter_t) );
         last->next  = tail;
         for (i=0; i<o->d.split.nobjs; i++) {
            tail->obj   = o->d.split.objs[i];
            ret        += syn_branch_iter_walk( syn, func, data, head, tail );
         }
         last->next  = NULL;
         free( tail );
         break;

      case KIN_TYPE_TCP:
         /* We found a head. Now all we have to do is process it. */
         func( syn, head, data );
         ret = 1; /* Mark that we found a head. */
         /* Theoretically we can have a TCP in the middle so we must consider
          * this possibility. */
         if (o->next != NULL) {
            tail        = memcalloc( 1, sizeof(kin_branch_iter_t) );
            last->next  = tail;
            tail->obj   = o->next;
            ret        += syn_branch_iter_walk( syn, func, data, head, tail );
            last->next  = NULL;
            free( tail );
         }
         break;

      default:
         /* This shouldn't happen. */
         assert( "type mismatch" == NULL );
         break;
   }
   return ret;
}


/**
 * @brief Multiple and accumulate for Lie algebra elements.
 *    @param[out] dest Destination buffer.
 *    @param[in] val Value to multiply the source by.
 *    @param[in] src Source element.
 */
static void lie_joint_mac( plucker_t *dest,
      double val, const plucker_t *src )
{
   dest->s[0]  = val * src->s[0];
   dest->s[1]  = val * src->s[1];
   dest->s[2]  = val * src->s[2];
   dest->s0[0] = val * src->s0[0];
   dest->s0[1] = val * src->s0[1];
   dest->s0[2] = val * src->s0[2];
}


/**
 * @brief Does the Lie bracket operation.
 *    @param[out] O Output Lie algebra element.
 *    @param[in] A Input element 1.
 *    @param[in] B Input element 2.
 */
static void lie_joint_bracket( plucker_t *O,
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
static void plucker_zero( plucker_t *p )
{
   memset( p, 0, sizeof(plucker_t) );
}


/**
 * @brief Obtains the core screw element from the dual quaternion.
 */
static void plucker_from_dq( plucker_t *p, const dq_t Q )
{
   p->s[0]  = Q[1];
   p->s[1]  = Q[2];
   p->s[2]  = Q[3];
   p->s0[0] = Q[4];
   p->s0[1] = Q[5];
   p->s0[2] = Q[6];
}


/**
 * @brief Subtracts a plucker coordinate from another.
 */
static void plucker_sub( plucker_t *o, const plucker_t *p, const plucker_t *q )
{
   vec3_sub( o->s,  p->s,  q->s );
   vec3_sub( o->s0, p->s0, q->s0 );
}


/**
 * @brief Calculates a branch.
 *
 * The tree should be already updated from the possible data vectors if being used in the solver.
 *
 * What we have to do is:
 *
 * *0) Convert variable vector to function vector.
 *  1) Set up mask data from mapped data (which should be read directly from claims).
 *  2) Update position function vector and screw axis current location.
 *  3) Update derivative function vector from new screw axis locations.
 *  4) Convert masked function vector to mapped data.
 *  5) Convert tree structure to function vector.
 *
 * * indicates this is not done in this function.
 */
int syn_calc_branch( synthesis_t *syn, kin_branch_t *branch )
{
   int i, k, m;
   dq_t T, R;
   kin_joint_t *j, *jk;
   double p[3] = { 0., 0., 0. };

   /* Branch must make sense. */
   assert( syn->branches != NULL );
   assert( branch->joints != NULL );

   /* Calculate plucker conditions for all the joints, we explicitly use them as equations. */
   for (i=0; i<branch->njoints; i++) {
      j          = branch->joints[i];
      j->cond[0] = (vec3_norm( j->S.s ) - 1.); /* ||s|| = 1 */
      j->cond[1] = vec3_dot( j->S.s, j->S.s0 ); /* s.s0 = 0 */

      /* Need to map data from memory to vector. */
      mm_updateMask( &j->pos.values );
      mm_updateMask( &j->vel.values );
      mm_updateMask( &j->acc.values );

      /* Initialize joint position to initial position. */
      memcpy( &j->S_cur, &j->S, sizeof(plucker_t) );
   }

   /* Calculate all the frame data.
    * This is tricky because we have L-1 positions and L velocities/accelerations.
    * So we handle them all at the same time, however we don't process the position
    * 0 which corresponds to the reference system.
    */
   for (m=0; m<syn->L; m++) {
      if (m != 0) {
         dq_cr_point( T, p ); /* Create identity dual quaternion. */

         /* Concatenate relative displacements around axes. */
         for (i=0; i<branch->njoints; i++) {
            dq_t S;
            double *d;
            j = branch->joints[i];

            /* Store current joint position. */
            dq_cr_line_plucker( S, j->S.s, j->S.s0 );
            dq_op_f2g( S, T, S );
            plucker_from_dq( &j->S_cur, S );

            /* Update propagation. */
            d = (double*) j->pos.values.mask_vec;
            dq_cr_rotation_plucker( R, d[m-1], j->S.s, j->S.s0 );
            dq_op_mul( T, T, R );
         }

         /* Since dual quaternions double map SE(3) we make it so all the
          * dual quaternions are using the same sign convention. */
         if (T[3] < 0.)
            dq_op_sign( T, T );

         /* Store in fvec: T-P. */
         dq_op_sub( branch->tcp->d.tcp.fvec_pos[m-1], T, branch->tcp->d.tcp.P[m] );
      }

      /* Update velocities. */
      if ((branch->tcp->d.tcp.claim_vel != NULL) &&
            ((branch->tcp->d.tcp.V.mask_mask == NULL) ||
             ((branch->tcp->d.tcp.V.mask_len > m) &&
              (branch->tcp->d.tcp.V.mask_mask[m] != 0)))) {
         /* Initialize result. */
         plucker_t v;
		   plucker_t *d, *V;
         plucker_zero( &v );
         /* Straight forward:
          *    v = \sum v_i S_i
          */
         for (i=0; i<branch->njoints; i++) {
            double *dv;
            j = branch->joints[i];
            dv = (double*) j->vel.values.mask_vec;
            lie_joint_mac( &v, dv[m], &j->S_cur );
         }
         /* Copy result. */
			d = (plucker_t*) branch->tcp->d.tcp.fvec_vel.mask_vec;
			V = (plucker_t*) branch->tcp->d.tcp.V.mask_vec;
         plucker_sub( &d[m], &v, &V[m] );
      }

      /* Update accelerations. */
      if ((branch->tcp->d.tcp.claim_acc != NULL) &&
            ((branch->tcp->d.tcp.A.mask_mask == NULL) ||
             ((branch->tcp->d.tcp.A.mask_len > m) &&
              (branch->tcp->d.tcp.A.mask_mask[m] != 0)))) {
         /* Initialize result. */
         plucker_t a;
		   plucker_t *d, *A;
         plucker_zero( &a );
         /* Straight forward sum part:
          *    a_1 = \sum a_i S_i
          */
         for (i=0; i<branch->njoints; i++) {
            double *dv;
            j = branch->joints[i];
            dv = (double*) j->acc.values.mask_vec;
            lie_joint_mac( &a, dv[m], &j->S_cur );
         }
         /* More complex coriolis part.
          *    a_2 = \sum a_i \sum a_k [s_i, s_k]
          */
         for (i=0; i<branch->njoints-1; i++) {
            double *dv;
            plucker_t acc;
            j = branch->joints[i];

            plucker_zero( &acc );
            for (k=i+1; k<branch->njoints; k++) {
               plucker_t br;
               jk = branch->joints[k];
               lie_joint_bracket( &br, &j->S_cur, &jk->S_cur );
               dv = (double*) jk->vel.values.mask_vec;
               lie_joint_mac( &acc, dv[m], &br );
            }
            dv = (double*) j->vel.values.mask_vec;
            lie_joint_mac( &a, dv[m], &acc );
         }
         /* Copy result. */
			d = (plucker_t*) branch->tcp->d.tcp.fvec_acc.mask_vec;
			A = (plucker_t*) branch->tcp->d.tcp.A.mask_vec;
         plucker_sub( &d[m], &a, &A[m] );
      }
   }

   /* Convert back to map layout for the solver. */
   mm_updateMap( &branch->tcp->d.tcp.fvec_vel );
   mm_updateMap( &branch->tcp->d.tcp.fvec_acc );

   return 0;
}


/**
 * @brief Calculates all the fvec, conditions and the likes for the kinematic synthesis.
 *
 * This rewrites syn->fvec, but does not load from syn->x.
 *
 *    @param syn Synthesis to calculate all the data from.
 */
int syn_calc( synthesis_t *syn )
{
   int i;

   /* Must be finalized. */
   assert( syn->finalized );

   /* Process all the branches. */
   for (i=0; i<syn->nbranches; i++)
      syn_calc_branch( syn, &syn->branches[i] );

   /* We have to map back to the fvector. */
   syn_map_to_fvec( syn, NULL, NULL, syn->fvec );
   return 0;
}


/**
 * @brief Walking over branches.
 *    @param syn Synthesis to iterate over branches.
 *    @param func Function to run on each branch.
 *    @param data User data to pass to function.
 *    @return Number of branches found.
 */
int syn_branch_iter( synthesis_t *syn, int (*func)(synthesis_t*, kin_branch_iter_t*, void*), void *data )
{
   int ret;
   kin_branch_iter_t *head;

   assert( syn->obj != NULL );

   head = memcalloc( 1, sizeof(kin_branch_iter_t) );
   head->obj = syn->obj;
   ret = syn_branch_iter_walk( syn, func, data, head, head );
   free(head);

   return ret;
}


/**
 * @brief Initializes a kinematic joint.
 *
 *    @param joint Joint to initialize.
 *    @param type Type of the joint.
 * @sa kin_joint_init
 */
int kin_joint_init( kin_joint_t *joint, kin_joint_type_t type )
{
   const double ub[6] = {  1.,  1.,  1.,  100.,  100.,  100. };
   const double lb[6] = { -1., -1., -1., -100., -100., -100. };

   memset( joint, 0, sizeof(kin_joint_t) );
   joint->type = type;

   /* Default bounds. */
   memcpy( &joint->S_lb, lb, 6*sizeof(double) );
   memcpy( &joint->S_ub, ub, 6*sizeof(double) );

   return 0;
}


/**
 * @brief Frees the kinematic joint data.
 */
static void kin_joint_data_free( kin_joint_data_t *data )
{
   mm_cleanup( &data->values );
   mm_cleanup( &data->values_lb );
   mm_cleanup( &data->values_ub );
#ifndef NDEBUG
   memset( data, 0, sizeof(kin_joint_data_t) );
#endif /* NDEBUG */
}


/**
 * @brief Frees the datastructure in joint.
 *
 *    @param joint Joint to free.
 * @sa kin_joint_init
 */
void kin_joint_free( kin_joint_t *joint )
{
   kin_joint_data_free( &joint->pos );
   kin_joint_data_free( &joint->vel );
   kin_joint_data_free( &joint->acc );
}


/**
 * @brief Duplicates data
 */
static void kin_joint_data_dup( kin_joint_data_t *dest, const kin_joint_data_t *src )
{
   dest->nvalues = src->nvalues;
   dest->constant = src->constant;

   if (src->values.chunk != 0)
      mm_initDup( &dest->values, &src->values );

   if (src->values_lb.chunk != 0)
      mm_initDup( &dest->values_lb, &src->values_lb );

   if (src->values_ub.chunk != 0)
      mm_initDup( &dest->values_ub, &src->values_ub );
}


/**
 * @brief Duplicates a joint.
 */
void kin_joint_dupInit( kin_joint_t *nj, const kin_joint_t *oj )
{
   memset( nj, 0, sizeof(kin_joint_t) );

   /* Copy base stuff. */
   nj->type      = oj->type;
   nj->const_S   = oj->const_S;
   nj->const_pos = oj->const_pos;

   /* Copy axis info. */
   memcpy( &nj->S,    &oj->S,    sizeof(plucker_t) );
   memcpy( &nj->S_lb, &oj->S_lb, sizeof(plucker_t) );
   memcpy( &nj->S_ub, &oj->S_ub, sizeof(plucker_t) );

   /* Duplicate position and derivative memory. */
   kin_joint_data_dup( &nj->pos, &oj->pos );
   kin_joint_data_dup( &nj->vel, &oj->vel );
   kin_joint_data_dup( &nj->acc, &oj->acc );

   /* Copy conditional data. */
   memcpy( nj->cond, oj->cond, 2*sizeof(double) );
}


/**
 * @brief Claims joint data.
 */
static int kin_joint_data_claim( synthesis_t *syn, kin_joint_data_t *data, const char *name )
{
   assert( data->claim == NULL );
   if (data->values.chunk != 0) {
      //assert( data->nvalues == syn->L-1 );
      data->claim = syn_claim_x( syn, data->values.map_len, data->values.map_len,
            (double*) data->values.map_vec,
            (double*) data->values_lb.map_vec,
            (double*) data->values_ub.map_vec, name );
   }
   return 0;
}


/**
 * @brief has a joint create claims.
 *    @param syn Synthesis object to which the joint belongs.
 *    @param joint Joint to do the claiming.
 */
int kin_joint_claim( synthesis_t *syn, kin_joint_t *joint )
{
   /* We must make sure it's not previously initialized and that data is valid. */
   assert( joint->claim_S == NULL );
   assert( joint->claim_cond == NULL );
   assert( joint->claim_const == NULL );
   assert( syn->L > 0 );

   /* Claim the data. */
   kin_joint_data_claim( syn, &joint->pos, "Joint Pos" );
   kin_joint_data_claim( syn, &joint->vel, "Joint Vel" );
   kin_joint_data_claim( syn, &joint->acc, "Joint Acc" );

   /* Claim joint axis if necessary. */
   if (!joint->const_S) {
      joint->claim_S    = syn_claim_x(    syn, 6, 4,
            (double*)&joint->S, (double*)&joint->S_lb, (double*)&joint->S_ub,
            "Joint Plucker" );
      /* The conditions are not independent as they are what lowers the
       * dimension from 6->4. */
      joint->claim_cond = syn_claim_fvec( syn, 2, 0, joint->cond,
            "Joint Plucker Cond" );
   }

   return 0;
}


/**
 * @brief Sets the plucker coordinates for a joint.
 *
 *    @param joint Joint to set plucker coordinates of.
 *    @param s Direction vector of the joint line.
 *    @param s0 Moment of the joint line.
 */
void kin_joint_setPlucker( kin_joint_t *joint, double s[3], double s0[3] )
{
   assert( joint->type != JOINT_TYPE_NULL );
   assert( fabs(vec3_dot( s, s0 )) <= DQ_PRECISION );
   memcpy( joint->S.s,  s,  sizeof(double)*3 );
   memcpy( joint->S.s0, s0, sizeof(double)*3 );
}


/**
 * @brief Sets the data for a joint.
 *    @param data Joint data to set.
 *    @param x Vector to set.
 *    @param len Length of vector.
 *    @param lb Default low bound value.
 *    @param ub Default upper bound value.
 */
static void kin_joint_data_set( kin_joint_data_t *data, double *x, int len, double lb, double ub, const int *mask )
{
   int i;
   double *pub, *plb;

   assert( data->values.chunk == 0 );

   /* The values are straight forward. */
   data->nvalues   = len;
   mm_initMask( &data->values, sizeof(double), len, x, mask );

   /* We'll create temporary vectors to create the bounds. */
   plb = memmalloc( len * sizeof(double) );
   pub = memmalloc( len * sizeof(double) );
   for (i=0; i<len; i++) {
      plb[i] = lb;
      pub[i] = ub;
   }
   mm_initMask( &data->values_lb, sizeof(double), len, plb, mask );
   mm_initMask( &data->values_ub, sizeof(double), len, pub, mask );
   free( plb );
   free( pub );
}


/**
 * @brief Sets the positions.
 *    @param joint Joint to set positions of.
 *    @param x X positions to set.
 *    @param len Length of the positions (should be frames-1).
 */
void kin_joint_setPositions( kin_joint_t *joint, double *x, int len )
{
   /* No mask, we always need position. */
   kin_joint_data_set( &joint->pos, x, len, 0., 2.*M_PI, NULL );
}


/**
 * @brief Sets the velocities.
 *    @param joint Joint to set velocities of.
 *    @param v Velocities to set.
 *    @param len Length of velocity vector.
 */
void kin_joint_setVelocities( kin_joint_t *joint, double *v, int len, const int *mask )
{
   kin_joint_data_set( &joint->vel, v, len, -M_PI, M_PI, mask );
}


/**
 * @brief Sets the accelerations.
 *    @param joint Joint to set accelerations of.
 *    @param a Accelerations to set.
 *    @param len Length of the accelerations vector.
 */
void kin_joint_setAccelerations( kin_joint_t *joint, double *a, int len, const int *mask )
{
   kin_joint_data_set( &joint->acc, a, len, -M_PI, M_PI, mask );
}


/**
 * @brief Sets the axis of the joint as static (not updated in solver).
 *    @param joint Joint to set as static.
 *    @param constant Whether or not to set as static (constant).
 */
void kin_joint_setConstS( kin_joint_t *joint, int constant )
{
   assert( joint->claim_S == NULL );
   joint->const_S = constant;
}


/**
 * @brief Sets the movements of the joint as static (not updated in solver).
 *    @param joint Joint to set as static.
 *    @param constant Whether or not to set as static (constant).
 */
void kin_joint_setConstPos( kin_joint_t *joint, int constant )
{
   assert( joint->pos.claim == NULL );
   joint->pos.constant = constant;
}


/**
 * @brief Sets the bounds of the plucker position coordinate.
 *    @param joint Joint to set plucker bounds of.
 *    @param S_lb Lower bound to set the axis orientation.
 *    @param S_ub Upper bound to set the axis orientation.
 *    @param S0_lb Lower bound to set the moment of the axis.
 *    @param S0_lb Upper bound to set the moment of the axis.
 */
void kin_joint_setPluckerBounds( kin_joint_t *joint,
      double *S_lb, double *S_ub, double *S0_lb, double *S0_ub )
{
#ifndef NDEBUG
   int i;
   for (i=0; i<3; i++) {
      assert( S_lb[i]  <= S_ub[i] );
      assert( S0_lb[i] <= S0_ub[i] );
   }
#endif /* NDEBUG */

   if (S_lb != NULL)
      memcpy( joint->S_lb.s,  S_lb,  sizeof(double)*3 );
   if (S_ub != NULL)
      memcpy( joint->S_ub.s,  S_ub,  sizeof(double)*3 );

   if (S0_lb != NULL)
      memcpy( joint->S_lb.s0, S0_lb, sizeof(double)*3 );
   if (S0_ub != NULL)
      memcpy( joint->S_ub.s0, S0_ub, sizeof(double)*3 );
}


/**
 * @brief Sets the data bounds.
 *    @param data Data to set bounds of.
 *    @param lb Lower bounds.
 *    @param ub Upper bonds.
 *    @param len Length.
 */
static void kin_joint_data_setBounds( kin_joint_data_t *data,
      double *lb, double *ub, int len )
{
#ifndef NDEBUG
   int i;
   for (i=0; i<len; i++)
      assert( lb[i] <= ub[i] );
#endif /* NDEBUG */
   assert( len == data->nvalues );
   assert( data->values.chunk != 0 );

   /* Copy values over. */
   mm_initMask( &data->values_lb, sizeof(double), len, lb, data->values.mask_mask );
   mm_initMask( &data->values_ub, sizeof(double), len, ub, data->values.mask_mask );
}


/**
 * @brief Sets the position bounds for the joint parameters in a joint.
 *    @param joint Joint to set position bounds of.
 *    @param lb Lower bound vector.
 *    @param ub Upper bound vector.
 *    @param len Length of both lower and upper bounds which must match length of joint parameter vector.
 */
void kin_joint_setPositionBounds( kin_joint_t *joint,
      double *lb, double *ub, int len )
{
   kin_joint_data_setBounds( &joint->pos, lb, ub, len );
}


/**
 * @brief Sets the velocity bounds for the joint parameters in a joint.
 *    @param joint Joint to set velocity bounds of.
 *    @param lb Lower bound vector.
 *    @param ub Upper bound vector.
 *    @param len Length of both lower and upper bounds which must match length of joint parameter vector.
 */
void kin_joint_setVelocityBounds( kin_joint_t *joint,
      double *lb, double *ub, int len )
{
   kin_joint_data_setBounds( &joint->vel, lb, ub, len );
}


/**
 * @brief Sets the acceleration bounds for the joint parameters in a joint.
 *    @param joint Joint to set acceleration bounds of.
 *    @param lb Lower bound vector.
 *    @param ub Upper bound vector.
 *    @param len Length of both lower and upper bounds which must match length of joint parameter vector.
 */
void kin_joint_setAccelerationBounds( kin_joint_t *joint,
      double *lb, double *ub, int len )
{
   kin_joint_data_setBounds( &joint->acc, lb, ub, len );
}


/**
 * @brief Frees a kinematic chain object data.
 *    @param obj Object to free data of.
 *    @return 0 on success.
 */
static int kin_obj_chain_free( kin_object_t *obj )
{
   int i;
   for (i=0; i<obj->d.chain.njoints; i++)
      kin_joint_free( &obj->d.chain.joints[i] );
   free( obj->d.chain.joints );
   return 0;
}
/**
 * @brief Duplicates the kinematic chain object data.
 *    @param obj Parent object to duplicate data of.
 *    @param newobj Child object to recieve a copy of the data.
 *    @return 0 on success.
 */
static int kin_obj_chain_dup( const kin_object_t *obj, kin_object_t *newobj )
{
   int i;

   /* Copy block over. */
   memset( &newobj->d.chain, 0, sizeof(kin_chain_data_t) );

   /* Copy joints over. */
   newobj->d.chain.joints  = memmalloc( obj->d.chain.njoints*sizeof(kin_joint_t) );
   newobj->d.chain.njoints = obj->d.chain.njoints;

   /* Duplicate joints. */
   for (i=0; i<newobj->d.chain.njoints; i++)
      kin_joint_dupInit( &newobj->d.chain.joints[i], &obj->d.chain.joints[i] );

   return 0;
}
/**
 * @brief Runs claim on a chain object.
 *    @param obj Kinematic object to claim data.
 *    @param syn Parent synthesis object.
 */
static int kin_obj_chain_fin( kin_object_t *obj, synthesis_t *syn )
{
   int i;
   /* Claim all joints. */
   for (i=0; i<obj->d.chain.njoints; i++) {
      kin_joint_claim( syn, &obj->d.chain.joints[i] );
      /* Add joints to synthesis. */
      syn_joint_add( syn, &obj->d.chain.joints[i] );
   }
   return 0;
}


static void kin_obj_chain_save_pos( char *str, char *str_lb, char *str_ub, int max,
      const kin_joint_data_t *data )
{
   int j, p, p_lb, p_ub;
   const char *tail;

   str[0]    = '\0';
   str_lb[0] = '\0';
   str_ub[0] = '\0';

   p     = 0;
   p_lb  = 0;
   p_ub  = 0;
   for (j=0; j<data->nvalues; j++) {
      tail = (j!=data->nvalues-1) ? "," : "";
      if (data->values.mask_mask[j]) {
         p    += snprintf( &str[p], max - p,
               PF"%s ", ((double*)data->values.mask_vec)[j],    tail );
         p_lb += snprintf( &str_lb[p_lb], max - p_lb,
               PF"%s ", ((double*)data->values_lb.mask_vec)[j], tail );
         p_ub += snprintf( &str_ub[p_ub], max - p_ub,
               PF"%s ", ((double*)data->values_ub.mask_vec)[j], tail );
      }
      else {
         p    += snprintf( &str[p],       max - p,    "nil%s ", tail );
         p_lb += snprintf( &str_lb[p_lb], max - p_lb, "nil%s ", tail );
         p_ub += snprintf( &str_ub[p_ub], max - p_ub, "nil%s ", tail );
      }
   }
}

/**
 * @brief Saves a kinematics object.
 */
static int kin_obj_chain_save( const kin_object_t *obj, FILE *stream, const char *self )
{
   int i;
   kin_joint_t *kj;
   const char *str;
   char pos_str[4096], pos_lb_str[4096], pos_ub_str[4096];
   double *s, *s0, *lb, *lb0, *ub, *ub0;

   for (i=0; i<obj->d.chain.njoints; i++) {
      kj = &obj->d.chain.joints[i];

      /* Joint information and creation. */
      str = "revolute";
      fprintf( stream,
            "   -- Joint %d\n"
            "   local j = kin_joint.new( \"%s\" )\n",
            i, str );

      /* Position information .*/
      kin_obj_chain_save_pos( pos_str, pos_lb_str, pos_ub_str, 4096, &kj->pos );
      fprintf( stream,
            "   j:setPositions( { %s } )\n"
            "   j:setPositionBounds( { %s },\n"
            "                        { %s } )\n",
            pos_str, pos_lb_str, pos_ub_str );
  
      /* Velocity information. */
      if (kj->vel.nvalues > 0) {
         kin_obj_chain_save_pos( pos_str, pos_lb_str, pos_ub_str, 4096, &kj->vel );
         fprintf( stream,
               "   j:setVelocities( { %s }, %d )\n"
               "   j:setVelocityBounds( { %s },\n"
               "                        { %s } )\n",
               pos_str, kj->vel.values.mask_len, pos_lb_str, pos_ub_str );
      }

      /* Acceleration information. */
      if (kj->acc.nvalues > 0) {
         kin_obj_chain_save_pos( pos_str, pos_lb_str, pos_ub_str, 4096, &kj->acc );
         fprintf( stream,
               "   j:setAccelerations( { %s }, %d )\n"
               "   j:setAccelerationBounds( { %s },\n"
               "                        { %s } )\n",
               pos_str, kj->acc.values.mask_len, pos_lb_str, pos_ub_str );
      }

      /* Plucker axis. */
      s     = kj->S.s;
      s0    = kj->S.s0;
      lb    = kj->S_lb.s;
      lb0   = kj->S_lb.s0;
      ub    = kj->S_ub.s;
      ub0   = kj->S_ub.s0;
      fprintf( stream,
            "   j:setPlucker( { "PF", "PF", "PF" }, { "PF", "PF", "PF" } )\n"
            "   j:setPluckerBounds( { "PF", "PF", "PF" }, { "PF", "PF", "PF" },\n"
            "                       { "PF", "PF", "PF" }, { "PF", "PF", "PF" } )\n",
            s[0], s[1], s[2], s0[0], s0[1], s0[2],
            lb[0], lb[1], lb[2],
            ub[0], ub[1], ub[2],
            lb0[0], lb0[1], lb0[2],
            ub0[0], ub0[1], ub0[2] );

      fprintf( stream, "   %s:attach( j )\n", self );
   }
   return 0;
}
/**
 * @brief Frees a kinematic TCP object data.
 *    @param obj Object to free data of.
 *    @return 0 on success.
 */
static int kin_obj_tcp_free( kin_object_t *obj )
{
   free( obj->d.tcp.P );
	mm_cleanup( &obj->d.tcp.V );
	mm_cleanup( &obj->d.tcp.A );
   free( obj->d.tcp.fvec_pos );
   mm_cleanup( &obj->d.tcp.fvec_vel );
   mm_cleanup( &obj->d.tcp.fvec_acc );
   return 0;
}
/**
 * @brief Duplicates the kinematic TCP object data.
 *    @param obj Parent object to duplicate data of.
 *    @param newobj Child object to recieve a copy of the data.
 *    @return 0 on success.
 */
static int kin_obj_tcp_dup( const kin_object_t *obj, kin_object_t *newobj )
{
   memset( &newobj->d.tcp, 0, sizeof(kin_tcp_data_t) );
   if (obj->d.tcp.P != NULL) {
      newobj->d.tcp.nP = obj->d.tcp.nP;
      newobj->d.tcp.P = memdup( obj->d.tcp.P, obj->d.tcp.nP * sizeof(dq_t) );
   }
   if (obj->d.tcp.V.chunk != 0)
		mm_initDup( &newobj->d.tcp.V, &obj->d.tcp.V );
   if (obj->d.tcp.A.chunk != 0)
		mm_initDup( &newobj->d.tcp.A, &obj->d.tcp.A );
   return 0;
}
/**
 * @brief Runs claim on a tcp object.
 *    @param obj Kinematic object to claim data.
 *    @param syn Parent synthesis object.
 */
static int kin_obj_tcp_fin( kin_object_t *obj, synthesis_t *syn )
{
	double *d;

   /* Claims. */
   assert( obj->d.tcp.fvec_pos == NULL );
   assert( obj->d.tcp.claim_pos == NULL );

   /* We note that since they are relative to the reference position there's
    * only m-1 error dual quaternions. */
   obj->d.tcp.fvec_pos     = memcalloc( (syn->L-1), sizeof(dq_t) );
   obj->d.tcp.claim_pos    = syn_claim_fvec( syn, 8*(syn->L-1),
         6*(syn->L-1), (double*)obj->d.tcp.fvec_pos,
         "TCP Pos" );

   /* Derivatives. */
	d = memcalloc( syn->L, sizeof(plucker_t) );
   /* Velocity. */
   if (obj->d.tcp.V.chunk != 0) {
      assert( obj->d.tcp.fvec_vel.chunk == 0 );
      assert( obj->d.tcp.claim_vel == NULL );

		mm_initMask( &obj->d.tcp.fvec_vel, sizeof(plucker_t),
				obj->d.tcp.V.mask_len, d, obj->d.tcp.V.mask_mask );
      obj->d.tcp.claim_vel = syn_claim_fvec( syn, 6*obj->d.tcp.V.map_len,
            6*obj->d.tcp.V.map_len, (double*)obj->d.tcp.fvec_vel.map_vec,
            "TCP Vel" );
   }
   /* Acceleration. */
   if (obj->d.tcp.A.chunk != 0) {
      assert( obj->d.tcp.fvec_acc.chunk == 0 );
      assert( obj->d.tcp.claim_acc == NULL );

		mm_initMask( &obj->d.tcp.fvec_acc, sizeof(plucker_t),
				obj->d.tcp.A.mask_len, d, obj->d.tcp.A.mask_mask );
      obj->d.tcp.claim_acc = syn_claim_fvec( syn, 6*obj->d.tcp.A.map_len,
            6*obj->d.tcp.A.map_len, (double*)obj->d.tcp.fvec_acc.map_vec,
            "TCP Acc" );
   }
	free(d);

   /* Add tcp to list. */
   syn_tcp_add( syn, obj );
   return 0;
}
/**
 * @brief Frees a kinematic splitter object data.
 *    @param obj Object to free data of.
 *    @return 0 on success.
 */
static int kin_obj_split_free( kin_object_t *obj )
{
   int i;
   for (i=0; i<obj->d.split.nobjs; i++)
      kin_obj_destroy( obj->d.split.objs[i] );
   free( obj->d.split.objs );
   return 0;
} 
static int kin_obj_tcp_save_derivative( const mm_vec_t *der, FILE *stream )
{
   int i;
   const char *tail;
   plucker_t *p;

   for (i=0; i<der->mask_len; i++) {
      fprintf( stream, "       -- Frame %d\n", i );
      tail = (i==der->mask_len-1) ? "" : ",";
      if (der->mask_mask[i]) {
         p = (plucker_t*) der->mask_vec;
         fprintf( stream, "       { "PF", "PF", "PF", "PF", "PF", "PF" }%s\n",
               p->s[0], p->s[1], p->s[2], p->s0[0], p->s0[1], p->s0[2], tail );
      }
      else
         fprintf( stream, "       nil%s\n", tail );
   }
   fprintf( stream,       "          }, %d )\n", der->mask_len );
   return 0;
}
static int kin_obj_tcp_save( const kin_object_t *obj, FILE *stream, const char *self )
{
   int i;
   dq_t Q;
   double R[3][3], d[3];

   fprintf( stream, "   %s:setFK( {\n", self );
   for (i=0; i<obj->d.tcp.nP; i++) {

      /* We can extract the information directly from each pose. */
      if (i > 0)
         dq_op_mul( Q, obj->d.tcp.P[i], obj->d.tcp.P[0] );
      else
         dq_cr_copy( Q, obj->d.tcp.P[0] );
      dq_op_extract( R, d, Q );

      /* We can just extract it from the end-effector position information .*/
      fprintf( stream, "       -- Frame %d\n"
                       "       { { "PF", "PF", "PF", "PF" },\n"
                       "         { "PF", "PF", "PF", "PF" },\n" 
                       "         { "PF", "PF", "PF", "PF" },\n" 
                       "         { "PF", "PF", "PF", "PF" } }%s\n" , i,
                       R[0][0], R[0][1], R[0][2], d[0],
                       R[1][0], R[1][1], R[1][2], d[1],
                       R[2][0], R[2][1], R[2][2], d[2],
                       0., 0., 0., 1.,
                       (i < obj->d.tcp.nP-1) ? "," : "" );
   }
   fprintf( stream, "             } )\n" );

   if (obj->d.tcp.V.mask_len > 0) {
      fprintf( stream, "   %s:setVel( {\n", self );
      kin_obj_tcp_save_derivative( &obj->d.tcp.V, stream );
   }

   if (obj->d.tcp.A.mask_len > 0) {
      fprintf( stream, "   %s:setAcc( {\n", self );
      kin_obj_tcp_save_derivative( &obj->d.tcp.A, stream );
   }

   return 0;
}
/**
 * @brief Duplicates the kinematic splitter object data.
 *    @param obj Parent object to duplicate data of.
 *    @param newobj Child object to recieve a copy of the data.
 *    @return 0 on success.
 */
static int kin_obj_split_dup( const kin_object_t *obj, kin_object_t *newobj )
{
   int i;
   memset( &newobj->d.split, 0, sizeof(kin_splitter_data_t) );
   newobj->d.split.objs = memmalloc( obj->d.split.nobjs*sizeof(kin_object_t*) );
   for (i=0; i<obj->d.split.nobjs; i++)
      newobj->d.split.objs[i] = kin_obj_dup( obj->d.split.objs[i] );
   newobj->d.split.nobjs = obj->d.split.nobjs;
   return 0;
}
/**
 * @brief Runs claim on a splitter object.
 *    @param obj Kinematic object to claim data.
 *    @param syn Parent synthesis object.
 */
static int kin_obj_split_fin( kin_object_t *obj, synthesis_t *syn )
{
   int i;
   /* Iterate over ends. */
   for (i=0; i<obj->d.split.nobjs; i++)
      kin_obj_fin( syn, obj->d.split.objs[i] );
   return 0;
}


/**
 * @brief Saves a splitter object.
 *    @param obj Kinematic object to save.
 *    @param stream Stream to save object to.
 *    @param self Name to use for self.
 */
static int kin_obj_split_save( const kin_object_t *obj, FILE *stream, const char *self )
{
   int i;
   for (i=0; i<obj->d.split.nobjs; i++) {
      kin_obj_save( obj->d.split.objs[i], stream, 1, "sp" );
      fprintf( stream, "   %s:attach( sp1 )\n", self );
   }
   return 0;
}


/**
 * @brief Creates a new kinematic object.
 *    @param type Type to create.
 *    @return Newly created kinematic object.
 * @sa kin_obj_destroy
 */
kin_object_t* kin_obj_create( kin_object_type_t type )
{
   kin_object_t *obj;
   obj         = memmalloc( sizeof(kin_object_t) );
   kin_obj_init( obj, type );
   return obj;
}


/**
 * @brief Duplicates a kinematic object.
 *    @param obj Object to duplicate.
 *    @return Duplicated object.
 */
kin_object_t* kin_obj_dup( const kin_object_t *obj )
{
   kin_object_t *newobj;
   newobj         = memcalloc( 1, sizeof(kin_object_t) );
   kin_obj_dupInit( obj, newobj );
   return newobj;
}


/**
 * @brief Initializes a duplicate kinematic object.
 *    @param obj Parent object to duplicate.
 *    @param newobj Child object to recieve duplicated data.
 *    @return 0 on success.
 */
int kin_obj_dupInit( const kin_object_t *obj, kin_object_t *newobj )
{
   memset( newobj, 0, sizeof(kin_object_t) );

   memcpy( &newobj->f, &obj->f, sizeof(kin_object_func_t) );
   newobj->type   = obj->type;
   if (obj->f.dup != NULL)
      obj->f.dup( obj, newobj );

   /* Duplicate next. */
   if (obj->next != NULL)
      newobj->next = kin_obj_dup( obj->next );

   return 0;
}


/**
 * @brief Destroys a kinematic object.
 *    @param obj Object to destroy.
 */
void kin_obj_destroy( kin_object_t *obj )
{
   /* Destroy next in linked list. */
   if (obj->next != NULL)
      kin_obj_destroy( obj->next );

   /* Free stuff. */
   kin_obj_free( obj );
   free(obj);
   return;
}


/**
 * @brief Initializes a kinematic object.
 *    @param obj Object to initialize.
 *    @param type Type of the object to initialize.
 *    @return 0 on success.
 */
int kin_obj_init( kin_object_t *obj, kin_object_type_t type )
{
   memset( obj, 0, sizeof(kin_object_t) );
   obj->type   = type;
   switch (type) {
      case KIN_TYPE_CHAIN:
         obj->f.free = kin_obj_chain_free;
         obj->f.dup  = kin_obj_chain_dup;
         obj->f.fin  = kin_obj_chain_fin;
         obj->f.save = kin_obj_chain_save;
         break;
      case KIN_TYPE_TCP:
         obj->f.free = kin_obj_tcp_free;
         obj->f.dup  = kin_obj_tcp_dup;
         obj->f.fin  = kin_obj_tcp_fin;
         obj->f.save = kin_obj_tcp_save;
         break;
      case KIN_TYPE_SPLITTER:
         obj->f.free = kin_obj_split_free;
         obj->f.dup  = kin_obj_split_dup;
         obj->f.fin  = kin_obj_split_fin;
         obj->f.save = kin_obj_split_save;
         break;
      default:
         assert( "Unknown object type!" == NULL );
         break;
   }
   return 0;
}


/**
 * @brief Frees a kinematic object data (does not discard the object itself, but it is invalid).
 *    @param obj Object to free.
 */
void kin_obj_free( kin_object_t *obj )
{
   obj->f.free( obj );
}


/**
 * @brief Runs claims for a kinematic object.
 *    @param syn Parent synthesis to run claims on.
 *    @param obj Object to run claims on.
 */
static int kin_obj_fin( synthesis_t *syn, kin_object_t *obj )
{
   if (obj->f.fin != NULL)
      obj->f.fin( obj, syn );
   if (obj->next != NULL)
      kin_obj_fin( syn, obj->next );
   return 0;
}


/**
 * @brief Gets the human readable name of an object type.
 */
static const char *kin_obj_type_to_str( kin_object_type_t type )
{
   switch (type) {
      case KIN_TYPE_NULL:     return "null";
      case KIN_TYPE_CHAIN:    return "chain";
      case KIN_TYPE_TCP:      return "tcp";
      case KIN_TYPE_SPLITTER: return "splitter";
      default:
         assert( "Object type" == NULL );
         return NULL;
   }
}


/**
 * @brief Saves a kinematic object.
 */
static int kin_obj_save( const kin_object_t *obj, FILE *stream, int depth, const char *sym )
{
   const char *str;
   char buf[32];

   /* Create the object. */
   str = kin_obj_type_to_str( obj->type );
   snprintf( buf, sizeof(buf), "%s%d", sym, depth );
   fprintf( stream, "   local %s = kin_object.new( \"%s\" )\n", buf, str );

   /* Save particular data. */
   if (obj->f.save != NULL)
      obj->f.save( obj, stream, buf );

   /* We'll attach the next ones as they come. */
   if (obj->next != NULL) {
      kin_obj_save( obj->next, stream, depth+1, sym );
      fprintf( stream, "   %s:attach( %s%d )\n", buf, sym, depth+1 );
   }
   return 0;
}


/**
 * @brief Attaches a kinematic object to another.
 *    @param parent Parent to attach the child to.
 *    @param child Child to attach to the parent.
 *    @return 0 on success.
 */
int kin_obj_attach( kin_object_t *parent, kin_object_t *child )
{
   assert( parent->next == NULL );
   parent->next = child;
   return 0;
}


/**
 * @brief Adds a joint to a kinematic chain.
 *    @param chain Chain to add joint to.
 *    @param joint Joint to add.
 *    @return 0 on success.
 */
int kin_obj_chain_joint_add( kin_object_t *chain, kin_joint_t *joint )
{
   assert( chain->type == KIN_TYPE_CHAIN );

   chain->d.chain.joints = memrealloc( chain->d.chain.joints, sizeof(kin_joint_t)*(chain->d.chain.njoints+1) );
   kin_joint_dupInit( &chain->d.chain.joints[chain->d.chain.njoints], joint );
   chain->d.chain.njoints++;

   return 0;
}


/**
 * @brief Sets the forward kinematics related to a TCP.
 *    @param tcp TCP to set FK data to.
 *    @param pos Positions of the FK.
 *    @param num Number of coordinates.
 *    @return 0 on success.
 */
int kin_obj_tcp_fk( kin_object_t *tcp, const dq_t *pos, int num )
{
   int i;
   dq_t Tinv;

   /* Only works on TCP. */
   assert( tcp->type == KIN_TYPE_TCP );

   if (tcp->d.tcp.nP == 0)
      tcp->d.tcp.nP = num;
   else {
      assert( tcp->d.tcp.nP == num );
   }

   /* We must allocate the proper size and duplicate it when loading. */
   tcp->d.tcp.P         = memdup( (void*)pos, num * sizeof(dq_t) );
   tcp->d.tcp.constant  = 1; /* Generally considered static. */

   /* We must do the calculations with the inversion of the reference.
    * Specifically we must construct it so we obtain Pi0 = Pi P0^{-1} */
   dq_cr_inv( Tinv, tcp->d.tcp.P[0] );

   /* Must update the other P dual quaternions. */
   for (i=1; i<tcp->d.tcp.nP; i++) {

      /* Do the calculation of multiplying itself. */
      dq_op_mul( tcp->d.tcp.P[i], tcp->d.tcp.P[i], Tinv );

      /* It is important to flip them as the dual cover SE(3). */
      if (tcp->d.tcp.P[i][3] < 0.)
         dq_op_sign( tcp->d.tcp.P[i], tcp->d.tcp.P[i] );
   }

   return 0;
}


/**
 * @brief Sets the velocity.
 */
int kin_obj_tcp_velocity( kin_object_t *tcp,
      const plucker_t *vel, const int *mask, int num )
{
   /* Only works on TCP. */
   assert( tcp->type == KIN_TYPE_TCP );

   assert( tcp->d.tcp.V.chunk == 0 );
	mm_initMask( &tcp->d.tcp.V, sizeof(plucker_t), num, (void*)vel, mask );
   return 0;
}


/**
 * @brief Sets the acceleration.
 */
int kin_obj_tcp_acceleration( kin_object_t *tcp,
      const plucker_t *acc, const int *mask, int num )
{
   /* Only works on TCP. */
   assert( tcp->type == KIN_TYPE_TCP );

   assert( tcp->d.tcp.A.chunk == 0 );
	mm_initMask( &tcp->d.tcp.A, sizeof(plucker_t), num, (void*)acc, mask );
   return 0;
}


/**
 * @brief Calculates how it has to skip.
 */
int syn_calc_skip( synthesis_t *syn, int total )
{
   int skip;
   
   skip = (int)floor((double)(total - 1) / ((double)syn->L));
   if (skip == 0)
      skip = 1;

   return skip;
}


/**
 * @brief Loads the TCP forward kinematics from an fk file.
 *    @param syn Synthesis object to which the tcp object belongs.
 *    @param obj TCP object.
 *    @param file File to load from.
 */
int kin_obj_tcp_fk_load( synthesis_t *syn, kin_object_t *obj, const char *file )
{
   int n, i;
   double R[3][3], d[3], unused[4];
   dq_t Tinv, Tp;
   FILE *fd;
   int skip, L_total;

   assert( syn->L > 0 );
   assert( obj->type == KIN_TYPE_TCP );

   /* Open file. */
   fd = fopen( file, "r" );
   assert( fd != NULL );

   /** Allocate information. */
   obj->d.tcp.nP  = syn->L;
   obj->d.tcp.P   = memcalloc( obj->d.tcp.nP, sizeof(dq_t) );

   /* Get length. */
   L_total = 0;
   while ((fscanf( fd,
                  "{{%lf, %lf, %lf, %lf},"
                  " {%lf, %lf, %lf, %lf},"
                  " {%lf, %lf, %lf, %lf},"
                  " {%lf, %lf, %lf, %lf}}\n",
                  &R[0][0], &R[0][1], &R[0][2], &d[0],
                  &R[1][0], &R[1][1], &R[1][2], &d[1],
                  &R[2][0], &R[2][1], &R[2][2], &d[2],
                  &unused[0], &unused[1], &unused[2], &unused[3])) == 16)
      L_total++;
   rewind( fd );

   /* Calculate skip. */
   skip = syn_calc_skip( syn, L_total );

   /* Read and parse the entire file. */
   n = 0;
   i = 0;
   while ((fscanf( fd,
                  "{{%lf, %lf, %lf, %lf},"
                  " {%lf, %lf, %lf, %lf},"
                  " {%lf, %lf, %lf, %lf},"
                  " {%lf, %lf, %lf, %lf}}\n",
                  &R[0][0], &R[0][1], &R[0][2], &d[0],
                  &R[1][0], &R[1][1], &R[1][2], &d[1],
                  &R[2][0], &R[2][1], &R[2][2], &d[2],
                  &unused[0], &unused[1], &unused[2], &unused[3])) == 16) {

      /* Hack to skip. */
      if (n % skip) {
         n++;
         continue;
      }

      /* Skip when done. */
      if (i >= syn->L)
         break;

      /* Set the home. */
      if (i==0) {
         /* Create the home position. We have to invert it because the calculation is:
          * T_{01} = T_1 T_0^{-1} */
         dq_cr_homo( obj->d.tcp.P[i], R, d );
         dq_cr_inv( Tinv, obj->d.tcp.P[i] );
      }
      else {
         /* Create the movements from the previous position.
          *
          * This can be done through the tranformation:
          *
          *    T_{0i} = T_i T_0^{-1}
          */
         dq_cr_homo( Tp, R, d );
         dq_op_mul( obj->d.tcp.P[i], Tp, Tinv );

         /* To ensure they are the same, we always convert to positive. */
         if (obj->d.tcp.P[i][3] < 0.)
            dq_op_sign( obj->d.tcp.P[i], obj->d.tcp.P[i] );
      }

      /* Increment line counter. */
      i++;
      n++;
   }

   /* Close file. */
   fclose( fd );

   /* Sanity check. */
   assert( i == syn->L );

   return 0;
}


/**
 * @brief Sets whether or not an end effector is constant.
 *
 * If it is not constant it will be calculated, this can be used for forward kinematics.
 *
 *    @param tcp Tool center point object (end effector).
 *    @param constant Whether or not it's constant (they default to constant).
 *    @return 0 on success.
 */
int kin_obj_tcp_fk_setConst( kin_object_t *tcp, int constant )
{
   assert( tcp->type == KIN_TYPE_TCP );
   tcp->d.tcp.constant = constant;
   return 0;
}


/**
 * @brief Attaches an object to a splitter object.
 *    @param split Splitter object to attach object to.
 *    @param obj Object to attach.
 *    @return 0 on success.
 */
int kin_obj_split_attach( kin_object_t *split, kin_object_t *obj )
{
   assert( split->type == KIN_TYPE_SPLITTER );

   split->d.split.objs = memrealloc( split->d.split.objs,
         sizeof(kin_object_t*)*(split->d.split.nobjs+1) );
   split->d.split.objs[split->d.split.nobjs] = obj;
   split->d.split.nobjs++;

   return 0;
}


/**
 * @brief Creates a synthesis object.
 *    @return newly created synthesis object.
 * @sa syn_destroy
 */
synthesis_t* syn_create (void)
{
   synthesis_t *syn;

   syn = memcalloc( 1, sizeof(synthesis_t) );

   return syn;
}


/**
 * @brief Destroys a synthesis object.
 *    @param syn Synthesis object to destroy.
 * @sa syn_create
 */
void syn_destroy( synthesis_t *syn )
{
   syn_free( syn );
   free( syn );
}


/**
 * @brief Initializes a synthesis object.
 *    @param syn Synthesis object ot initialize.
 *    @return 0 on success.
 */
int syn_init( synthesis_t *syn )
{
   memset( syn, 0, sizeof(synthesis_t) );
   return 0;
}


/**
 * @brief Frees a synthesis object.
 *    @param syn Synthesis object to free.
 */
void syn_free( synthesis_t *syn )
{
   int i;

   free( syn->x );
   free( syn->fvec );
   free( syn->lb );
   free( syn->ub );

   if (syn->obj != NULL)
      kin_obj_destroy( syn->obj );

   if (syn->claim_x != NULL)
      syn_claim_destroy( syn->claim_x );
   if (syn->claim_fvec != NULL)
      syn_claim_destroy( syn->claim_fvec );

   free( syn->joints );
   free( syn->tcp );

   if (syn->branches != NULL) {
      for (i=0; i<syn->nbranches; i++)
         syn_branch_free( &syn->branches[i] );
      free( syn->branches );
   }

   /* MEMORY MASH. */
   memset( syn, 0, sizeof(synthesis_t) );
}


/**
 * @brief Duplicates a synthesis object.
 *    @param syn Synthesis object to duplicate.
 *    @return Duplicated synthesis object.
 */
synthesis_t* syn_dup( const synthesis_t *syn )
{
   synthesis_t *new_syn;
   new_syn = memmalloc( sizeof(synthesis_t) );
   syn_copy( new_syn, syn );
   return new_syn;
}


/**
 * @brief Duplicates memory.
 */
static void* memdup( const void *base, size_t size )
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
static void* memcalloc( size_t nmemb, size_t size )
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
static void* memmalloc( size_t size )
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
static void* memrealloc( void *ptr, size_t size )
{
   void *mem;
   mem = realloc( ptr, size );
   assert( mem != NULL );
   return mem;
}


/**
 * @brief Copies a synthesis object.
 *    @param syn Destination synthesis object.
 *    @param base Base synthesis object.
 */
int syn_copy( synthesis_t *syn, const synthesis_t *base )
{
   memset( syn, 0, sizeof(synthesis_t) );

   /* Base properties. */
   syn->n      = base->n;
   syn->m      = base->m;
   syn->ni     = base->ni;
   syn->mi     = base->mi;
   syn->L      = base->L;

   /* Allocate objects. */
   if (base->obj != NULL)
      syn->obj = kin_obj_dup( base->obj );

   /* Finalize if finalized. */
   if (base->finalized)
      syn_finalize( syn );

   return 0;
}


/**
 * @brief Does forward kinematics on the synthesis object.
 *    @param[in] syn Synthesis object to perform forward kinematics on.
 *    @param[out] fk Forward kinematics results (length of syn->ntcp).
 *    @param[in] angles Angles of each joint (length of syn->njoints).
 */
int syn_fk( const synthesis_t *syn, dq_t *fk, const double *angles )
{
   int i, j;
   kin_branch_t *kb;
   kin_joint_t *kj;
   dq_t R;
   double kz[3] = { 0., 0., 0. };

   assert( syn->finalized );

   /* Load angle information. */
   for (i=0; i<syn->njoints; i++)
      syn->joints[i]->pos_cur = angles[i];

   /* Now we must perform forward kinematics on all the branches. */
   for (i=0; i<syn->nbranches; i++) {
      kb = &syn->branches[i];

      /* Forward kinematics. */
      dq_cr_point( fk[i], kz );
      for (j=0; j<kb->njoints; j++) {
         kj = kb->joints[j];

         /* Must propagate the joint rotation. */
         dq_cr_rotation_plucker( R, kj->pos_cur, kj->S.s, kj->S.s0 );
         dq_op_mul( fk[i], fk[i], R );
      }

      /* Initial travsformation P. */
      dq_op_mul( fk[i], fk[i], kb->tcp->d.tcp.P[0] );
   }

   return 0;
}


/**
 * @brief Creates a new kinematic claim.
 *    @param claim Parent claim to append to (or NULL if it should be standalone).
 *    @param len Length of the claim.
 *    @param v Vector for the data of the claim.
 *    @return The newly created kinematic claim.
 */
static kin_claim_t* syn_claim_add( kin_claim_t *claim,
      size_t len, size_t indep, double *v, double *lb, double *ub,
      const char *name )
{
   kin_claim_t *l, *n;
   size_t p;

   assert( v != NULL );
   if (len <= 0)
      return NULL;

   /* Set data. */
   n        = memcalloc( 1, sizeof(kin_claim_t) );
   n->dat   = v;
   n->lb    = lb;
   n->ub    = ub;
   n->size  = len;
   n->indep = indep;

   /* Get offset. */
   if (claim != NULL) {
      p = 0;
      for (l=claim; l!=NULL; l=l->next)
         p += l->size;
      n->offset = p;

      /* Go to the end. */
      for (l=claim; l->next!=NULL; l=l->next);
      l->next  = n; /* Attach to the end. */
   }
   else
      n->offset = 0;

   /* We store the name as it can highlight some possible issues. */
#ifdef NDEBUG
   (void) name;
#else /* NDEBUG */
   n->name  = name;
#endif /* NDEBUG */

   return n;
}


/**
 * @brief Destroys a kinematic claim.
 *    @param claim Claim to destroy.
 */
static void syn_claim_destroy( kin_claim_t *claim )
{
   if (claim->next != NULL)
      syn_claim_destroy( claim->next );

   free(claim);
}


/**
 * @brief Claims data from the X vector.
 *    @param syn Synthesis object to claim from.
 *    @param len Length to claim (units is doubles).
 *    @param indep Independent length.
 *    @param v Data vector doing the claiming.
 *    @return The newly created claim.
 */
static kin_claim_t* syn_claim_x( synthesis_t *syn, size_t len, size_t indep,
      double *v, double *lb, double *ub, const char *name )
{
   kin_claim_t *n;
   n = syn_claim_add( syn->claim_x, len, indep, v, lb, ub, name );
   if (syn->claim_x == NULL)
      syn->claim_x = n;
   return n;
}


/**
 * @brief Claims data from the fvec vector.
 *    @param syn Synthesis object to claim from.
 *    @param len Length to claim (units is doubles).
 *    @param indep Independent length.
 *    @param v Data vector doing the claiming.
 *    @return The newly created claim.
 */
static kin_claim_t* syn_claim_fvec( synthesis_t *syn, size_t len, size_t indep,
      double *v, const char *name )
{
   kin_claim_t *n;
   n = syn_claim_add( syn->claim_fvec, len, indep, v, NULL, NULL, name );
   if (syn->claim_fvec == NULL)
      syn->claim_fvec = n;
   return n;
}


/**
 * @brief Adds a joint to the list of joints used by the synthesis object.
 *    @param parent Synthesis object to add a joint to.
 *    @param joint Joint to add.
 *    @return 0 on success.
 */
static int syn_joint_add( synthesis_t *parent, kin_joint_t *joint )
{
   parent->joints = memrealloc( parent->joints, sizeof(kin_joint_t*)*(parent->njoints+1) );
   parent->joints[parent->njoints] = joint;
   parent->njoints++;
   return 0;
}


/**
 * @brief Adds a TCP to the list of TCP used by the synthesis object.
 *    @param parent Synthesis object to add the tcp to.
 *    @param tcp TCP object to add.
 *    @return 0 on success.
 */
static int syn_tcp_add( synthesis_t *parent, kin_object_t *tcp )
{
   parent->tcp = memrealloc( parent->tcp, sizeof(kin_object_t*)*(parent->ntcp+1) );
   parent->tcp[parent->ntcp] = tcp;
   parent->ntcp++;
   return 0;
}


/**
 * @brief Adds a kinematic object to a synthesis object.
 *    @param parent Synthesis object to add a kinematic object to.
 *    @param obj Kinematic object to add to synthesis object.
 *    @return 0 on success.
 */
int syn_object_add( synthesis_t *parent, kin_object_t *obj )
{
   kin_object_t *l;
   assert( !parent->finalized );

   if (parent->obj == NULL) {
      parent->obj = obj;
      return 0;
   }

   for (l=parent->obj; l->next!=NULL; l=l->next);
   l->next = obj;
   return 0;
}


/**
 * @brief Sets the synthesis frames.
 *    @param syn Synthesis to set frame number of.
 *    @param frames Frames to set.
 */
int syn_set_frames( synthesis_t *syn, int frames )
{
   assert( frames > 0 );
   assert( !syn->finalized );
   syn->L = frames;
   return 0;
}


/**
 * @brief Finalizes the synthesis object creating intermediate structures.
 *    @param syn Synthesis object to finalize.
 *    @return 0 on success.
 */
int syn_finalize( synthesis_t *syn )
{
   assert( syn->L > 0 );
   assert( !syn->finalized );
   assert( syn->x == NULL );
   assert( syn->fvec == NULL );

   /* Run claims. */
   if (syn->obj != NULL)
      kin_obj_fin( syn, syn->obj );

   assert( syn->claim_x != NULL );
   assert( syn->claim_fvec != NULL );

   /* Get variable size. */
   syn_map_to_x(    syn, &syn->n, &syn->ni, NULL );
   syn_map_to_fvec( syn, &syn->m, &syn->mi, NULL );

   /* Need enough independent variables. */
   assert( syn->ni <= syn->mi );

   /* Allocate. */
   syn->x    = memcalloc( syn->n, sizeof(double) );
   syn->fvec = memcalloc( syn->m, sizeof(double) );

   /* Handle the bindings. */
   syn->lb   = memcalloc( syn->n, sizeof(double) );
   syn->ub   = memcalloc( syn->n, sizeof(double) );
   syn_set_bounds( syn );

   /* Remap. */
   syn_map_to_x(    syn, NULL, NULL, syn->x );
   syn_map_to_fvec( syn, NULL, NULL, syn->fvec );

   /* Sanity check, make sure initial point is within bounds. */
   /* We don't actually want this check, because the results could be totalyl screwed up. */
   int i;
   kin_joint_t *kj;
   for (i=0; i<syn->njoints; i++) {
      kj = syn->joints[i];
      assert( kj->pos.nvalues == syn->L-1 );
#if 0
      /* Do not actually use because the boundries can be outstepped by the solver,
       * so this would cause it to assert and to fail when the kinematic structure
       * is valid. */
      int j;
      for (j=0; j<3; j++) {
         assert( kj->S_lb.s[j]  <= kj->S.s[j] );
         assert( kj->S_ub.s[j]  >= kj->S.s[j] );
         assert( kj->S_lb.s0[j] <= kj->S.s0[j] );
         assert( kj->S_ub.s0[j] >= kj->S.s0[j] );
      }
      for (j=0; j<kj->npos; j++) {
         assert( kj->pos_lb[j] <= kj->pos[j] );
         assert( kj->pos_ub[j] >= kj->pos[j] );
      }
#endif /* NDEBUG */
   }

   /* Create branches. */
   syn_construct_branches( syn );

   /* Mark as finalized. */
   syn->finalized = 1;

   return 0;
}


/**
 * @brief Resets the finalized state of the synthesis object.
 */
int syn_reset( synthesis_t *syn )
{
   int i;

   /* Must be finalized. */
   if (!syn->finalized)
      return 0;

   /* Free vectors. */
   free( syn->x );
   free( syn->fvec );
   free( syn->lb );
   free( syn->ub );
   syn->x      = NULL;
   syn->fvec   = NULL;
   syn->lb     = NULL;
   syn->ub     = NULL;

   /* Free the claims. */
   if (syn->claim_x != NULL)
      syn_claim_destroy( syn->claim_x );
   if (syn->claim_fvec != NULL)
      syn_claim_destroy( syn->claim_fvec );
   syn->claim_x      = NULL;
   syn->claim_fvec   = NULL;

   /* Remove joints and tcp. */
   free( syn->joints );
   free( syn->tcp );

   /* Remove branches. */
   if (syn->branches != NULL) {
      for (i=0; i<syn->nbranches; i++)
         syn_branch_free( &syn->branches[i] );
      free( syn->branches );
   }
   syn->branches  = NULL;
   syn->nbranches = 0;

   /* Mark as reset. */
   syn->finalized = 0;
   return 0;
}


/**
 * @brief Maps a a claim to a double vector.
 *    @param par Claim to map.
 *    @param n Length of mapped vector (NULL to ignore).
 *    @param ni Independent variables (NULL to ignore).
 *    @param v Mapped vector (NULL to ignore).
 *    @return 0 on success.
 */
static int syn_map_claim_to_vec( const kin_claim_t *par, int *n, int *ni, double *v, double *lb, double *ub )
{
   int p, pi;
   const kin_claim_t *l;
   size_t size;

   assert( par != NULL );

   p  = 0;
   pi = 0;
   for (l=par; l!=NULL; l=l->next) {
      size = l->size*sizeof(double);
      /* Set data. */
      if (v != NULL)
         memcpy( &v[p], l->dat, size );
      /* Set bounds. */
      if ((lb != NULL) && (l->lb != NULL))
         memcpy( &lb[p], l->lb, size );
      if ((ub != NULL) && (l->ub != NULL))
         memcpy( &ub[p], l->ub, size );
      /* Increment. */
      p  += l->size;
      pi += l->indep;
   }

   /* Set value. */
   if (n != NULL)
      *n = p;
   if (ni != NULL)
      *ni = pi;

   return 0;
}


/**
 * @brief Updates the internal bounds.
 */
static int syn_set_bounds( synthesis_t *syn )
{
   return syn_map_claim_to_vec( syn->claim_x, NULL, NULL, syn->x, syn->lb, syn->ub );
}


/**
 * @brief Maps the X claims to a vector.
 *    @param syn Synthesis object to map X claims.
 *    @param n Length of mapped vector (NULL to ignore).
 *    @param ni Number of independent variables.
 *    @param v Mapped vector (NULL to ignore).
 *    @return 0 on success.
 */
int syn_map_to_x( const synthesis_t *syn, int *n, int *ni, double *x )
{
   return syn_map_claim_to_vec( syn->claim_x, n, ni, x, NULL, NULL );
}


/**
 * @brief Maps the fvec claims to a vector.
 *    @param syn Synthesis object to map fvec claims.
 *    @param m Length of mapped vector (NULL to ignore).
 *    @param mi Number of independent fvec.
 *    @param v Mapped vector (NULL to ignore).
 *    @return 0 on success.
 */
int syn_map_to_fvec( const synthesis_t *syn, int *m, int *mi, double *fvec )
{
   return syn_map_claim_to_vec( syn->claim_fvec, m, mi, fvec, NULL, NULL );
}


/**
 * @brief Maps to X claims from a double vector.
 *    @param syn Synthesis object to map X vector to.
 *    @param x Vector to map from.
 *    @param n Length of x vector.
 *    @return 0 on success.
 */
int syn_map_from_x( synthesis_t *syn, const double *x, int n )
{
   int p;
   kin_claim_t *l;

   assert( x != NULL );

   p = 0;
   for (l=syn->claim_x; l!=NULL; l=l->next) {
      if (x != NULL)
         memcpy( l->dat, &x[p], sizeof(double)*l->size );
      p += l->size;
   }

   assert( p == n );

   return 0;
}


/**
 * @brief Compares two synthesis structures.
 *
 *    @param syn_a First synthesis structure to compare.
 *    @param syn_b Second synthesis structure to compare.
 *    @return 0 if equal, 1 otherwise.
 */
int syn_cmp( const synthesis_t *syn_a, const synthesis_t *syn_b )
{
   int i;

   if (syn_a->n != syn_b->n)
      return 1;
   if (syn_a->m != syn_b->m)
      return 1;
   if (syn_a->ni != syn_b->ni)
      return 1;
   if (syn_a->mi != syn_b->mi)
      return 1;
   if (syn_a->L != syn_b->L)
      return 1;
   if (syn_a->finalized != syn_b->finalized)
      return 1;

   syn_map_to_x( syn_a, NULL, NULL, syn_a->x );
   syn_map_to_x( syn_b, NULL, NULL, syn_b->x );

   for (i=0; i<syn_a->n; i++)
      if (fabs(syn_a->x[i]-syn_b->x[i]) > DQ_PRECISION)
         return 1;

   return 0;
}


/**
 * @brief Loads a synthesis object from data in a directory.
 *    @param path Path to directory to load synthesis object from.
 *    @return The synthesis loaded or NULL on failure.
 */
synthesis_t* syn_load( const char *path )
{
   (void) path;
   return NULL;
}


/**
 * @brief Saves a synthesis object to a directory.
 *    @param syn Synthesis object to save.
 *    @param path Path to save synthesis object to.
 *    @return 0 on success.
 */
int syn_save( const synthesis_t *syn, const char *path )
{
   FILE *fp;

   fp = fopen( path, "w" );
   assert( fp != NULL );

   fprintf( fp, "require \"synthesis\"\n"
                "\n"
                "function saved_syn ()\n"
                "   -- Create new synthesis object\n"
                "   local s = syn.new( %d )\n"
                "\n",
                syn->L );

   /* Add object. */
   if (syn->obj != NULL) {
      kin_obj_save( syn->obj, fp, 1, "ko" );
      fprintf( fp, "   s:addObject( ko1 )\n" );
   }

   if (syn->finalized)
      fprintf( fp, "   s:finalize()\n" );
   fprintf( fp, "   return s\n"
                "end\n" );

   /*
   fprintf( fp, "if type(package.loaded[myname]) ~= \"userdata\" then\n"
                "   local s = saved_syn ()\n"
                "   s:print()\n"
                "end\n" );
   */

   fclose( fp );

   return 0;
}


/**
 * @brief Prints a synthesis object to stdout.
 *    @param syn Synthesis object to print.
 */
void syn_print( const synthesis_t *syn )
{
   syn_printf( stdout, syn );
}


/**
 * @brief Prints a synthesis object to a stream.
 *    @param stream Stream to print to.
 *    @param syn Synthesis object to print.
 */
void syn_printf( FILE* stream, const synthesis_t *syn )
{
   int i, j;
   double *t, *t0, err;
   fprintf( stream, "n=%d, m=%d, ni=%d, mi=%d\n",
         syn->n, syn->m, syn->ni, syn->mi );
   for (i=0; i<syn->nbranches; i++) {
      fprintf( stream, "   [%d]: %d\n", i, syn->branches[i].njoints );
      for (j=0; j<syn->branches[i].njoints; j++) {
         t  = syn->branches[i].joints[j]->S.s;
         t0 = syn->branches[i].joints[j]->S.s0;
         fprintf( stream, "      %d:  [ %.3f, %.3f, %.3f ] x [ %.3f, %.3f, %.3f ]\n", j,
               t[0], t[1], t[2], t0[0], t0[1], t0[2] );
      }
   }
   err = 0.;
   for (i=0; i<syn->m; i++)
      err += fabs( syn->fvec[i] );
   fprintf( stream, "sum of error: %.3e\n", err );
}


/**
 * @brief Prints the synthesis in detail to stdout.
 */
void syn_printDetail( const synthesis_t *syn )
{
   syn_printfDetail( stdout, syn );
}


/**
 * @brief Prints the synthesis in detail to a FILE stream.
 */
void syn_printfDetail( FILE* stream, const synthesis_t *syn )
{
   int i, j, k;
   double *t, *t0, err;
   kin_joint_t *kj;

   for (i=0; i<syn->nbranches; i++) {
      fprintf( stream, "   [%d]: %d\n", i, syn->branches[i].njoints );
      for (j=0; j<syn->branches[i].njoints; j++) {
         kj = syn->branches[i].joints[j];
         t  = kj->S.s;
         t0 = kj->S.s0;
         fprintf( stream, "      %d:  [ %.3e, %.3e, %.3e ] x [ %.3e, %.3e, %.3e ] (%.3e)\n", j,
               t[0], t[1], t[2], t0[0], t0[1], t0[2],
               vec3_dot( t, t0 ) );
         mm_updateMask( &kj->pos.values );
         mm_updateMask( &kj->vel.values );
         mm_updateMask( &kj->acc.values );
         for (k=0; k<kj->pos.values.mask_len; k++) {
            fprintf( stream, "         " );
            fprintf( stream, " [%02d] %+.3e", k, ((double*)kj->pos.values.mask_vec)[k] );
            if ((k < kj->vel.values.mask_len) && (kj->vel.values.mask_mask[k]))
               fprintf( stream, " [%02d] %+.3e", k, ((double*)kj->vel.values.mask_vec)[k] );
            if ((k < kj->acc.values.mask_len) && (kj->acc.values.mask_mask[k]))
               fprintf( stream, " [%02d] %+.3e", k, ((double*)kj->acc.values.mask_vec)[k] );
            fprintf( stream, "\n" );
         }
      }
   }
   err = 0.;
   for (i=0; i<syn->m; i++)
      err += fabs( syn->fvec[i] );
   fprintf( stream, "sum of error: %.3e\n", err );
}



