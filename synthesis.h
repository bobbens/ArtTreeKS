

#ifndef _SYNTHESIS_H
#  define _SYNTHESIS_H


#include <stdio.h>
#include <malloc.h>

#include <dq/dq.h>

#include "mapmask.h"


/**
 * @brief Representation of a line using plucker coordinates.
 */
typedef struct plucker_s {
   double s[3];   /**< Direction of plucker coordinates. */
   double s0[3];  /**< Moment of the plucker coordinate line. */
} __attribute__((__packed__)) plucker_t;


/* Forward declaration of sythesis type. */
struct synthesis_s;
typedef struct synthesis_s synthesis_t;


/* Forward declaration of claim type. */
struct kin_claim_s;
typedef struct kin_claim_s kin_claim_t;


/* Forward declaration. */
struct kin_object_s;
typedef struct kin_object_s kin_object_t;


/**
 * @brief Joint types.
 */
typedef enum kin_joint_type_e {
   JOINT_TYPE_NULL,        /**< Invalid or non-existant joint. */
   JOINT_TYPE_REVOLUTE,    /**< Revolutoin joint. */
   JOINT_TYPE_PRISMATIC,   /**< Unsupported. */
} kin_joint_type_t;


/**
 * @brief Joint data.
 */
typedef struct kin_joint_data_s {
   int nvalues;         /**< Number of values. */
   mm_vec_t values;     /**< Values of each value. */
   mm_vec_t values_lb;  /**< Lower bounds for each value. */
   mm_vec_t values_ub;  /**< Upper bounds for each value. */
   kin_claim_t *claim;  /**< Claim data. */
   int constant;        /**< Is the data constant? */
} kin_joint_data_t;


/**
 * @brief Kinematic joint structure.
 *
 * Joint data depends on the joint type. In the case of rotation we always use
 * radians, while translations will be in the reference unit being used.
 */
typedef struct kin_joint_s {
   int const_S;            /**< Whether or not the axis of rotation of the joint is constant. */
   int const_pos;          /**< Whether or not the positions of movement of the joint are static. */
   kin_joint_type_t type;  /**< Type of the joint. */

   /* Plucker coordinates in absolute of the joint. */
   plucker_t S;      /**< Plucker coordinates of the axis of rotation. */
   plucker_t S_lb;   /**< Plucker lower bound. */
   plucker_t S_ub;   /**< Plucker upper bound. */

   /* Current data. */
   kin_joint_data_t pos; /**< Position data. */
   kin_joint_data_t vel; /**< Velocity data. */
   kin_joint_data_t acc; /**< Acceleration data. */
   plucker_t S_cur;  /**< Current joint position, used for calculating. */
   double cond[2];   /**< Stores the plucker conditions. */

   /* Claims. */
   kin_claim_t *claim_S;     /**< Claim of plucker variables (structural). */
   kin_claim_t *claim_cond;  /**< Claim of plucker conditions. */
   kin_claim_t *claim_const; /**< Constraint offset. */

   /* Temporary data. */
   double   pos_cur; /**< Used when iterating. */
} kin_joint_t;


/**
 * @brief Type of kinematic object.
 */
typedef enum kin_object_type_e {
   KIN_TYPE_NULL,       /**< Invalid or nonexistant object. */
   KIN_TYPE_CHAIN,      /**< Kinematic serial chain. */
   KIN_TYPE_TCP,        /**< Tool center point. */
   KIN_TYPE_SPLITTER,   /**< Splits into various objects. */
} kin_object_type_t;


/**
 * @brief Data for a kinamitc chain.
 */
typedef struct kin_chain_data_s {
   kin_joint_t *joints;  /**< Joints used. */
   int         njoints;  /**< Number of joints. */
} kin_chain_data_t;


/**
 * @brief Data for a tool center point (TCP).
 */
typedef struct kin_tcp_data_s {
   /* Data itself. */
   dq_t *P;    /**< Poses, the first is absolute (G), the rest are actually relative transformations.
                    To convert to absolute multiply by the first (G). */
   int  nP;    /**< Number of poses to calculate for. */
   mm_vec_t V; /**< Velocity information. */
   mm_vec_t A; /**< Acceleration information. */
   int constant; /**< Whether or not the TCP changes. */

   dq_t *fvec_pos; /**< Position coordinate variables. */
   mm_vec_t fvec_vel; /**< Velocity coordinate variables. */
   mm_vec_t fvec_acc; /**< Acceleration coordinate values. */

   /* Pointer information. */
   kin_claim_t *claim_pos; /**< Offset in fvec for positions. */
   kin_claim_t *claim_vel; /**< Offset in fvec for velocities. */
   kin_claim_t *claim_acc; /**< Offset in fvec for accelerations. */
} kin_tcp_data_t;


/**
 * @brief Data for the kinematic tree splitter.
 */
typedef struct kin_splitter_data_s {
   kin_object_t **objs;  /**< Objects. */
   int           nobjs;  /**< Number of objects. */
} kin_splitter_data_t;


/**
 * @brief Function pointers to copy easily.
 */
typedef struct kin_object_func_s {
   int (*free)( kin_object_t* );  /**< Free self. */
   int (*dup)( const kin_object_t*, kin_object_t* ); /**< Duplicate data. */
   int (*fin)( kin_object_t*, synthesis_t* ); /**< Run data claims. */
   int (*save)( const kin_object_t*, FILE *stream,
         const char *self ); /**< Saves the object to a stream. */
} kin_object_func_t;


/**
 * @brief Represents a large kinematic object.
 */
struct kin_object_s {
   kin_object_t *next;     /**< Linked list. */

   /* Function pointers. */
   kin_object_func_t f;    /**< Function pointers to use. */

   /* Data and type information. */
   kin_object_type_t type; /**< Type of kinematic object. */
   union {
      kin_chain_data_t  chain;      /**< Chain object. */
      kin_tcp_data_t    tcp;        /**< TCP object. */
      kin_splitter_data_t split;    /**< Splitter object. */
   } d;
};


/**
  @brief A memory claim.
 */
struct kin_claim_s {
   struct kin_claim_s *next;  /**< Linked list. */
   double *dat;               /**< Data to connect to the claim. */
   double *lb;                /**< Lower bound of the data. */
   double *ub;                /**< Upper bound of the data. */
   size_t size;               /**< Size of the claim. */
   size_t indep;              /**< Independent length. */
   size_t offset;             /**< Offset. */
};


/**
 * A kinematic branch.
 */
typedef struct kin_branch_s {
   kin_joint_t **joints; /**< Joints in the branch. */
   int          njoints; /**< Number of joints in the branch. */
   kin_object_t *tcp;    /**< TCP of the branch. */
} kin_branch_t;


/**
 * @brief The mother and father of the sheep.
 */
struct synthesis_s {
   /* Base properties. */
   int n;         /**< Size of x. */
   int m;         /**< Size of fvec. */
   int ni;        /**< Independent x. */
   int mi;        /**< Independent fvec. */
   double *x;     /**< Variables, dimension of n. */
   double *ub;    /**< Upper bound of variables, dimension of n. */
   double *lb;    /**< Lower bound of variables, dimension of n. */
   double *fvec;  /**< Function vector. */
   int L;         /**< Number of frames, determines FK and angles needed. */
   int finalized; /**< Is it finalized? */

   /* Topological layout. */
   kin_object_t *obj;      /**< Kinematic chain. */

   /* Claim.s */
   kin_claim_t *claim_x;   /**< x vector claims. */
   kin_claim_t *claim_fvec;/**< fvec vector claims. */

   /* Joints. */
   kin_joint_t **joints;   /**< Pointers to joints available. */
   int          njoints;   /**< Number of joints. */

   /* Tool center points. */
   kin_object_t **tcp;     /**< Pointers to the tool center points. */
   int           ntcp;     /**< Number of tool center points. */

   /* Branches. */
   kin_branch_t *branches; /**< Branch pointers. */
   int          nbranches; /**< Number of branches. */
};


/**
 * @brief Structure for walking over branches.
 */
typedef struct kin_branch_iter_s {
   struct kin_branch_iter_s *next;  /**< Linked list style. */
   kin_object_t *obj;               /**< Object in current node. */
} kin_branch_iter_t;


/*
 * BRANCH WALKING
 */
int syn_calc_branch( synthesis_t *syn, kin_branch_t *branch );
int syn_branch_iter( synthesis_t *syn, int (*func)(synthesis_t*, kin_branch_iter_t*, void*), void *data );


/*
 * KINEMATIC JOINT
 */
/* Creation/cleanup. */
int kin_joint_init( kin_joint_t *joint, kin_joint_type_t type );
void kin_joint_free( kin_joint_t *joint );
/* Set up. */
void kin_joint_dupInit( kin_joint_t *nj, const kin_joint_t *oj );
int kin_joint_claim( synthesis_t *syn, kin_joint_t *joint );
void kin_joint_setPlucker( kin_joint_t *joint, double s[3], double s0[3] );
void kin_joint_setPositions( kin_joint_t *joint, double *x, int len );
void kin_joint_setVelocities( kin_joint_t *joint, double *v, int len, const int *mask );
void kin_joint_setAccelerations( kin_joint_t *joint, double *a, int len, const int *mask );
void kin_joint_setConstS( kin_joint_t *joint, int constant );
void kin_joint_setConstPos( kin_joint_t *joint, int constant );
void kin_joint_setPluckerBounds( kin_joint_t *joint,
      double *S_lb, double *S_ub, double *S0_lb, double *S0_ub );
void kin_joint_setPositionBounds( kin_joint_t *joint,
      double *lb, double *ub, int len );
void kin_joint_setVelocityBounds( kin_joint_t *joint,
      double *lb, double *ub, int len );
void kin_joint_setAccelerationBounds( kin_joint_t *joint,
      double *lb, double *ub, int len );


/*
 * KINEMATIC OBJECT
 */
/* Creation/cleanup. */
kin_object_t* kin_obj_create( kin_object_type_t type );
kin_object_t* kin_obj_dup( const kin_object_t *obj );
int kin_obj_dupInit( const kin_object_t *obj, kin_object_t *newobj );
void kin_obj_destroy( kin_object_t *obj );
int kin_obj_init( kin_object_t *obj, kin_object_type_t type );
void kin_obj_free( kin_object_t *obj );

/* Attaching. */
int kin_obj_attach( kin_object_t *parent, kin_object_t *child );

/* 
 * Type specific.
 */
/* Chain. */
int kin_obj_chain_joint_add( kin_object_t *chain, kin_joint_t *joint );
/* TCP. */
int kin_obj_tcp_fk( kin_object_t *tcp, const dq_t *pos, int num );
int kin_obj_tcp_velocity( kin_object_t *tcp, const plucker_t *vel, const int *mask, int num );
int kin_obj_tcp_acceleration( kin_object_t *tcp, const plucker_t *acc, const int *mask, int num );
int kin_obj_tcp_fk_load( synthesis_t *syn, kin_object_t *obj, const char *file );
int kin_obj_tcp_fk_setConst( kin_object_t *tcp, int constant );
/* Splitter. */
int kin_obj_split_attach( kin_object_t *split, kin_object_t *obj );


/*
 * SYNTHESIS
 */
/* Creation/cleanup. */
synthesis_t* syn_create (void);
void syn_destroy( synthesis_t *syn );
int syn_init( synthesis_t *syn );
void syn_free( synthesis_t *syn );
synthesis_t* syn_dup( const synthesis_t *syn );
int syn_copy( synthesis_t *syn, const synthesis_t *base );

/* Forward kinematics. */
int syn_fk( const synthesis_t *syn, dq_t *fk, const double *angles );

/* Saving/loading. */
synthesis_t* syn_load( const char *path );
int syn_save( const synthesis_t *syn, const char *path );

/* Addition of stuff. */
int syn_object_add( synthesis_t *parent, kin_object_t *obj );

/* Preperation of synthesis. */
int syn_set_frames( synthesis_t *syn, int frames );
int syn_finalize( synthesis_t *syn );
int syn_reset( synthesis_t *syn );

/* Calculation. */
int syn_calc( synthesis_t *syn );
int syn_calc_skip( synthesis_t *syn, int total );

/* Mappings. */
int syn_map_to_x(    const synthesis_t *syn, int *n, int *ni, double *x );
int syn_map_to_fvec( const synthesis_t *syn, int *m, int *mi, double *fvec );
int syn_map_from_x(  synthesis_t *syn, const double *x, int n );
/*int syn_map_from_fvec(   synthesis_t *syn, double *fvec,  int m );*/

/* Compariosn. */
int syn_cmp( const synthesis_t *syn_a, const synthesis_t *syn_b );

/* Display. */
void syn_print( const synthesis_t *syn );
void syn_printf( FILE* stream, const synthesis_t *syn );
void syn_printDetail( const synthesis_t *syn );
void syn_printfDetail( FILE* stream, const synthesis_t *syn );


#endif /* _SYNTHESIS_H */



