/*
 * ArtTreeKS: Finite dimensional kinematic synthesis for articulated trees
 * Copyright(C) 2010-2012 Edgar Simo-Serra <esimo@iri.upc.edu>
 * License: see synthesis.h
*/

#include <assert.h>
#include <string.h>

#define LUA_LIB
#include "lua.h"
#include "lauxlib.h"

#include "synthesis.h"
#include "solver.h"
#include "visualizer.h"


#define DEBUG printf


/*
 * Generic.
 */
static int lua_reg_metatable( lua_State *L, const luaL_reg *reg, const char *name );
static int lua_is_foo( lua_State *L, int ind, const char *foo );
static int lua_table_to_vec( double v[3], lua_State *L, int ind );

/*
 * Kinematic Joint API.
 */
#define JOINT_METATABLE    "kin_joint"
/* Internal API. */
static int lua_loadjoint( lua_State *L );
static int lua_isjoint( lua_State *L, int ind );
static kin_joint_t* lua_tojoint( lua_State *L, int ind );
static kin_joint_t* luaL_checkjoint( lua_State *L, int ind );
static kin_joint_t* lua_pushjoint( lua_State *L, kin_joint_t );
/* Exposed API. */
static int jointL_gc( lua_State *L );
static int jointL_new( lua_State *L );
static int jointL_setPlucker( lua_State *L );
static int jointL_setPositions( lua_State *L );
static int jointL_setVelocities( lua_State *L );
static int jointL_setAccelerations( lua_State *L );
static int jointL_setConstS( lua_State *L );
static int jointL_setConstPos( lua_State *L );
static int jointL_setPluckerBounds( lua_State *L );
static int jointL_setPositionBounds( lua_State *L );
static int jointL_setVelocityBounds( lua_State *L );
static int jointL_setAccelerationBounds( lua_State *L );
static const luaL_reg jointL_methods[] = {
   { "__gc", jointL_gc },
   { "new", jointL_new },
   { "setPlucker", jointL_setPlucker },
   { "setPositions", jointL_setPositions },
   { "setVelocities", jointL_setVelocities },
   { "setAccelerations", jointL_setAccelerations },
   { "setConstS", jointL_setConstS },
   { "setConstPos", jointL_setConstPos },
   { "setPluckerBounds", jointL_setPluckerBounds },
   { "setPositionBounds", jointL_setPositionBounds },
   { "setVelocityBounds", jointL_setVelocityBounds },
   { "setAccelerationBounds", jointL_setAccelerationBounds },
   { 0, 0 }
};


/*
 * Kinematic Object API.
 */
#define OBJECT_METATABLE   "kin_object"
/* Internal API. */
static int lua_loadobject( lua_State *L );
static int lua_isobject( lua_State *L, int ind );
static kin_object_t* lua_toobject( lua_State *L, int ind );
static kin_object_t* luaL_checkobject( lua_State *L, int ind );
static kin_object_t* lua_pushobject( lua_State *L, kin_object_t obj );
/* Exposed API. */
static int objectL_gc( lua_State *L );
static int objectL_new( lua_State *L );
static int objectL_attach( lua_State *L );
static int objectL_setFK( lua_State *L );
static int objectL_setVel( lua_State *L );
static int objectL_setAcc( lua_State *L );
static const luaL_reg objectL_methods[] = {
   { "__gc", objectL_gc },
   { "new", objectL_new },
   { "attach", objectL_attach },
   { "setFK", objectL_setFK },
   { "setVel", objectL_setVel },
   { "setAcc", objectL_setAcc },
   { 0, 0 }
};


/*
 * Synthesis Object API.
 */
#define SYN_METATABLE      "syn"
/* Internal API. */
static int lua_loadsyn( lua_State *L );
static int lua_issyn( lua_State *L, int ind );
static synthesis_t* lua_tosyn( lua_State *L, int ind );
static synthesis_t* luaL_checksyn( lua_State *L, int ind );
static synthesis_t* lua_pushsyn( lua_State *L, synthesis_t syn );
/* Exposed API. */
static int synL_gc( lua_State *L );
static int synL_cmp( lua_State *L );
static int synL_new( lua_State *L );
static int synL_dup( lua_State *L );
static int synL_addObj( lua_State *L );
static int synL_finalize( lua_State *L );
#ifdef HAVE_MINPACK
static int synL_solve_fminpack( lua_State *L );
#endif /* HAVE_MINPACK */
static int synL_solve_minpack( lua_State *L );
#ifdef HAVE_NLOPT
static int synL_solve_nlopt( lua_State *L );
#endif /* HAVE_NLOPT */
static int synL_solve_cmaes( lua_State *L );
static int synL_solve_ga( lua_State *L );
static int synL_print( lua_State *L );
static int synL_printClaim( lua_State *L );
static int synL_printJacobian( lua_State *L );
static int synL_save( lua_State *L );
static int synL_raw_x( lua_State *L );
static int synL_raw_fvec( lua_State *L );
static int synL_visualize( lua_State *L );
static int synL_stats( lua_State *L );
static const luaL_reg synL_methods[] = {
   { "__gc", synL_gc },
   { "__eq", synL_cmp },
   { "new", synL_new },
   { "cmp", synL_cmp },
   { "dup", synL_dup },
   { "addObject", synL_addObj },
   { "finalize", synL_finalize },
   { "solver_minpack", synL_solve_minpack },
#ifdef HAVE_MINPACK
   { "solver_fminpack", synL_solve_fminpack },
#endif /* HAVE_MINPACK */
   { "solver_ga", synL_solve_ga },
#ifdef HAVE_NLOPT
   { "solver_nlopt", synL_solve_nlopt },
#endif /* HAVE_NLOPT */
   { "solver_cmaes", synL_solve_cmaes },
   { "print", synL_print },
   { "printClaim", synL_printClaim },
   { "printJacobian", synL_printJacobian },
   { "save", synL_save },
   { "raw_x", synL_raw_x },
   { "raw_fvec", synL_raw_fvec },
   { "visualize", synL_visualize },
   { "stats", synL_stats },
   { 0, 0 }
};


/**
 * @brief Loads the Lua synthesis library.
 */
int luaopen_synthesis( lua_State *L )
{
   lua_loadjoint(L);
   lua_loadobject(L);
   lua_loadsyn(L);
   return 0;
}


/*
 * Generic.
 */
/**
 * @brief Wrapper to register a metatable.
 *    @param L Lua state to register metatable in.
 *    @param reg Functions to add.
 *    @param name Name to use for the metatable.
 */
static int lua_reg_metatable( lua_State *L, const luaL_reg *reg, const char *name )
{
   luaL_newmetatable(   L, name );
   lua_pushvalue(       L, -1 );
   lua_setfield(        L, -2,   "__index" );
   luaL_register(       L, NULL, reg );
   lua_setfield(        L,       LUA_GLOBALSINDEX, name );
   return 0;
}
/**
 * @brief Checks to see if an index is a metatable.
 *    @param L Lua state to get index from.
 *    @param ind Index within the Lua state.
 *    @param foo Name of the metatable to check.
 */
static int lua_is_foo( lua_State *L, int ind, const char *foo )
{
   int ret;
   if (lua_getmetatable(L,ind)==0)
      return 0;
   lua_getfield(L, LUA_REGISTRYINDEX, foo);
   ret = 0;
   if (lua_rawequal(L, -1, -2)) /* does it have the correct mt? */
      ret = 1;
   lua_pop(L, 2); /* remove both metatables */
   return ret;
}
/**
 * @brief Maps a Lua table to a 3D vector.
 *    @param[out] v Vector to map to.
 *    @param L Lua state to get index from.
 *    @param ind Index to get table from.
 */
static int lua_table_to_vec( double v[3], lua_State *L, int ind )
{
   int i;
   luaL_checktype(L,ind,LUA_TTABLE);
   for (i=0; i<3; i++) {
      lua_pushnumber(L,i+1);
      lua_gettable(L,ind);
      v[i] = lua_tonumber(L,-1);
      lua_pop(L,1);
   }
   return 0;
}


/*
 * Kinematic Joint API.
 */
static int lua_loadjoint( lua_State *L )
{
   return lua_reg_metatable( L, jointL_methods, "kin_joint" );
}
static int lua_isjoint( lua_State *L, int ind )
{
   return lua_is_foo( L, ind, JOINT_METATABLE );
}
static kin_joint_t* lua_tojoint( lua_State *L, int ind )
{
   return (kin_joint_t*) lua_touserdata(L,ind);
}
static kin_joint_t* luaL_checkjoint( lua_State *L, int ind )
{
   if (lua_isjoint(L,ind))
      return lua_tojoint(L,ind);
   luaL_typerror(L, ind, JOINT_METATABLE);
   return NULL;
}
static kin_joint_t* lua_pushjoint( lua_State *L, kin_joint_t kj )
{
   kin_joint_t *j;
   j  = (kin_joint_t*) lua_newuserdata( L, sizeof(kin_joint_t) );
   assert( j != NULL );
   kin_joint_dupInit( j, &kj );
   luaL_getmetatable(L, JOINT_METATABLE);
   lua_setmetatable(L, -2);
   return j;
}
/* Exposed API. */
static int jointL_gc( lua_State *L )
{
   kin_joint_t *kj = luaL_checkjoint(L,1);
   kin_joint_free( kj );
   return 0;
}
static int jointL_new( lua_State *L )
{
   kin_joint_t kj;
   const char *str;
   kin_joint_type_t type;
   /* String is handled for basic creation. */
   if (lua_isstring(L,1)) {
      str = lua_tostring(L,1);
      if (strcasecmp(str,"revolute")==0)
         type = JOINT_TYPE_REVOLUTE;
      else if (strcasecmp(str,"prismatic")==0)
         type = JOINT_TYPE_PRISMATIC;
      else {
         luaL_error(L, "Invalid joint type '%s'.", str);
         return 0;
      }
      kin_joint_init( &kj, type );
      lua_pushjoint( L, kj );
      return 1;
   }
   luaL_error( L, "Wrong parameters when creating kinematic joint." );
   return 0;
}
static int jointL_setPlucker( lua_State *L )
{
   double s[3], s0[3];
   kin_joint_t *kj = luaL_checkjoint(L,1);
   lua_table_to_vec( s,  L, 2 );
   lua_table_to_vec( s0, L, 3 );
   kin_joint_setPlucker( kj, s, s0 );
   return 0;
}
static double* jointL_setParameter( lua_State *L, int *n )
{
   int i;
   double *v;

   /* Check parameter. */
   luaL_checktype(L, 2, LUA_TTABLE);
   if (lua_isnumber(L,3))
      (*n) = lua_tonumber(L,3);
   else
      (*n) = (int)lua_objlen(L,2);

   /* Create and map vector. */
   v = malloc( (size_t)(*n)*sizeof(double) );
   for (i=0; i<*n; i++) {
      lua_pushnumber(L,i+1);
      lua_gettable(L,2);
      v[i] = lua_tonumber(L,-1);
      lua_pop(L,1);
   }
   return v;
}
static int jointL_setPositions( lua_State *L )
{
   int n;
   double *v;

   /* Check parameter. */
   kin_joint_t *kj = luaL_checkjoint(L,1);
   v = jointL_setParameter( L, &n );

   /* Set data and clean up. */
   kin_joint_setPositions( kj, v, n );
   free(v);
   return 0;
}
static int jointL_setVelocities( lua_State *L )
{
   int n;
   double *v;

   /* Check parameter. */
   kin_joint_t *kj = luaL_checkjoint(L,1);
   v = jointL_setParameter( L, &n );

   /* Set data and clean up. */
   kin_joint_setVelocities( kj, v, n, NULL );
   free(v);
   return 0;
}
static int jointL_setAccelerations( lua_State *L )
{
   int n;
   double *v;

   /* Check parameter. */
   kin_joint_t *kj = luaL_checkjoint(L,1);
   v = jointL_setParameter( L, &n );

   /* Set data and clean up. */
   kin_joint_setAccelerations( kj, v, n, NULL );
   free(v);
   return 0;
}
static int jointL_setConstS( lua_State *L )
{
   kin_joint_t *kj = luaL_checkjoint(L,1);
   int b = lua_toboolean(L,2);
   kin_joint_setConstS(kj,b);
   return 0;
}
static int jointL_setConstPos( lua_State *L )
{
   kin_joint_t *kj = luaL_checkjoint(L,1);
   int b = lua_toboolean(L,2);
   kin_joint_setConstPos(kj,b);
   return 0;
}
static int jointL_setPluckerBounds( lua_State *L )
{
   double S_lb[3], S_ub[3], S0_lb[3], S0_ub[3];
   int top;

   /* Handle parameters. */
   kin_joint_t *kj = luaL_checkjoint(L,1);
   lua_table_to_vec( S_lb, L, 2 );
   lua_table_to_vec( S_ub, L, 3 );
   top = lua_gettop(L);
   if (top > 3) {
      lua_table_to_vec( S0_lb, L, 4 );
      lua_table_to_vec( S0_ub, L, 5 );
   }

   /* Set bounds. */
   kin_joint_setPluckerBounds( kj, S_lb, S_ub,
         (top>3) ? S0_lb : NULL,
         (top>3) ? S0_ub : NULL );

   return 0;
}
static void jointL_setBounds( lua_State *L, double **lb, double **ub, int *n )
{
   int i;

   /* Check parameters. */
   luaL_checktype(L, 2, LUA_TTABLE);
   luaL_checktype(L, 3, LUA_TTABLE);

   /* Allocate temporary memory. */
   if (lua_isnumber(L,4))
      (*n) = lua_tonumber(L,4);
   else
      (*n) = (int)lua_objlen(L,2);
   (*lb) = malloc( (size_t)(*n)*sizeof(double) );
   (*ub) = malloc( (size_t)(*n)*sizeof(double) );

   /* Construct vectors. */
   for (i=0; i<(*n); i++) {
      lua_pushnumber(L,i+1);
      lua_gettable(L,2);
      (*lb)[i] = lua_tonumber(L,-1);
      lua_pop(L,1);
      lua_pushnumber(L,i+1);
      lua_gettable(L,3);
      (*ub)[i] = lua_tonumber(L,-1);
      lua_pop(L,1);
   }
}
static int jointL_setPositionBounds( lua_State *L )
{
   int n;
   double *lb, *ub;

   /* Check parameters. */
   kin_joint_t *kj = luaL_checkjoint(L,1);
   jointL_setBounds( L, &lb, &ub, &n );

   /* Set and clean up. */
   kin_joint_setPositionBounds( kj, lb, ub, n );
   free(lb);
   free(ub);
   return 0;
}
static int jointL_setVelocityBounds( lua_State *L )
{
   int n;
   double *lb, *ub;

   /* Check parameters. */
   kin_joint_t *kj = luaL_checkjoint(L,1);
   jointL_setBounds( L, &lb, &ub, &n );

   /* Set and clean up. */
   kin_joint_setVelocityBounds( kj, lb, ub, n );
   free(lb);
   free(ub);
   return 0;
}
static int jointL_setAccelerationBounds( lua_State *L )
{
   int n;
   double *lb, *ub;

   /* Check parameters. */
   kin_joint_t *kj = luaL_checkjoint(L,1);
   jointL_setBounds( L, &lb, &ub, &n );

   /* Set and clean up. */
   kin_joint_setAccelerationBounds( kj, lb, ub, n );
   free(lb);
   free(ub);
   return 0;
}


/*
 * Kinematic Object API.
 */
static int lua_loadobject( lua_State *L )
{
   return lua_reg_metatable( L, objectL_methods, "kin_object" );
}
static int lua_isobject( lua_State *L, int ind )
{
   return lua_is_foo( L, ind, OBJECT_METATABLE );
}
static kin_object_t* lua_toobject( lua_State *L, int ind )
{
   return (kin_object_t*) lua_touserdata(L,ind);
}
static kin_object_t* luaL_checkobject( lua_State *L, int ind )
{
   if (lua_isobject(L,ind))
      return lua_toobject(L,ind);
   luaL_typerror(L, ind, OBJECT_METATABLE);
   return NULL;
}
static kin_object_t* lua_pushobject( lua_State *L, kin_object_t obj )
{
   kin_object_t *o;
   o  = (kin_object_t*) lua_newuserdata( L, sizeof(kin_object_t) );
   assert( o != NULL );
   kin_obj_dupInit( &obj, o );
   luaL_getmetatable(L, OBJECT_METATABLE);
   lua_setmetatable(L, -2);
   return o;
}
static int objectL_gc( lua_State *L )
{
   kin_object_t *obj = luaL_checkobject(L,1);
   kin_obj_free( obj );
   return 0;
}
static int objectL_new( lua_State *L )
{
   kin_object_t obj;
   const char *str;
   kin_object_type_t type;
   /* String is handled for basic creation. */
   if (lua_isstring(L,1)) {
      str = lua_tostring(L,1);
      if (strcasecmp(str,"chain")==0)
         type = KIN_TYPE_CHAIN;
      else if (strcasecmp(str,"tcp")==0)
         type = KIN_TYPE_TCP;
      else if (strcasecmp(str,"splitter")==0)
         type = KIN_TYPE_SPLITTER;
      else {
         luaL_error(L, "Invalid joint type '%s'.", str);
         return 0;
      }
      kin_obj_init( &obj, type );
      lua_pushobject( L, obj );
      return 1;
   }
   luaL_error( L, "Wrong parameters when creating kinematic object." );
   return 0;
}
static int objectL_attach( lua_State *L )
{
   kin_object_t *obj = luaL_checkobject(L,1);
   if (lua_isobject(L,2)) {
      kin_object_t *obj2 = lua_toobject(L,2);
      if (obj->type == KIN_TYPE_SPLITTER)
         kin_obj_split_attach( obj, kin_obj_dup(obj2) );
      else
         kin_obj_attach( obj, kin_obj_dup(obj2) );
      return 0;
   }
   else if (lua_isjoint(L,2)) {
      if (obj->type != KIN_TYPE_CHAIN) {
         luaL_error(L, "Can only add joints to kinematic chain objects.");
         return 0;
      }
      kin_joint_t *kj = lua_tojoint(L,2);
      kin_obj_chain_joint_add( obj, kj );
      return 0;
   }
   luaL_error( L, "Invalid object to attach." );
   return 0;
}
static int objectL_setFK( lua_State *L )
{
   int i, j, k, n;
   double H[4][4], R[3][3], d[3];
   dq_t *Q;
   kin_object_t *obj = luaL_checkobject(L,1);
   if (obj->type != KIN_TYPE_TCP) {
      luaL_error(L, "Object must be of type 'tcp'.");
      return 0;
   }

   /* Load as matrix. */
   luaL_checktype(L,2,LUA_TTABLE);
   n = (int)lua_objlen(L,2);
   Q = calloc( (size_t)n, sizeof(dq_t) );
   /* Get all the matrix. */
   for (k=0; k<n; k++) {
      lua_pushnumber(L,k+1); /* t, i */
      lua_gettable(L,2);     /* t, k */
      /* Get first row. */
      for (i=0; i<4; i++) {
         lua_pushnumber(L,i+1); /* t, i */
         lua_gettable(L,-2);     /* t, k */
         /* Get second row. */
         for (j=0; j<4; j++) {
            lua_pushnumber(L,j+1);  /* t, k, i */
            lua_gettable(L,-2);     /* t, k, n */
            H[i][j] = lua_tonumber(L,-1);
            lua_pop(L,1);           /* t, k */
         }
         lua_pop(L,1);           /* t */
      }
      lua_pop(L,1);           /* t */
      /* Convert to dual quaternion. */
      for (i=0; i<3; i++) {
         for (j=0; j<3; j++)
            R[i][j] = H[i][j];
         d[i] = H[i][3];
      }
      dq_cr_homo( Q[k], R, d );
   }

   /* Add. */
   kin_obj_tcp_fk( obj, (const dq_t*) Q, n );
   free( Q );

   return 0;
}
static int objectL_setDerivative( double **Q, int **mask, int *n, lua_State *L )
{
   int i, k;

   /* Load as matrix. */
   luaL_checktype(L,2,LUA_TTABLE);
   if (lua_isnumber(L,3))
      (*n) = lua_tonumber(L,3);
   else
      (*n) = (int)lua_objlen(L,2);
   (*Q) = calloc( (size_t)(*n), sizeof(plucker_t) );
   (*mask) = calloc( (size_t)(*n), sizeof(int) );
   /* Get all the matrix. */
   for (k=0; k<(*n); k++) {
      lua_pushnumber(L,k+1); /* t, i */
      lua_gettable(L,2);     /* t, k */
      if (lua_isnil( L, -1 )) {
         /* No data. */
         (*mask)[k] = 0;
      }
      else {
         /* Get the data. */
         for (i=0; i<6; i++) {
            lua_pushnumber(L,i+1); /* t, i */
            lua_gettable(L,-2);     /* t, k */
            (*Q)[6*k+i] = lua_tonumber(L,-1); /* t, k */
            lua_pop(L,1);           /* t */
         }
         (*mask)[k] = 1;
      }
      lua_pop(L,1);           /* t */
   }

   return 0;
}
static int objectL_setVel( lua_State *L )
{
   double *Q;
   int n, *mask;
   kin_object_t *obj = luaL_checkobject(L,1);
   if (obj->type != KIN_TYPE_TCP) {
      luaL_error(L, "Object must be of type 'tcp'.");
      return 0;
   }

   objectL_setDerivative( &Q, &mask, &n, L );
   kin_obj_tcp_velocity( obj, (const plucker_t*) Q, mask, n );
   free( Q );
   free( mask );

   return 0;
}
static int objectL_setAcc( lua_State *L )
{
   double *Q;
   int n, *mask;
   kin_object_t *obj = luaL_checkobject(L,1);
   if (obj->type != KIN_TYPE_TCP) {
      luaL_error(L, "Object must be of type 'tcp'.");
      return 0;
   }

   objectL_setDerivative( &Q, &mask, &n, L );
   kin_obj_tcp_acceleration( obj, (const plucker_t*) Q, mask, n );
   free( Q );
   free( mask );

   return 0;
}


/*
 * Synthesis Object API.
 */
/* Internal API. */
static int lua_loadsyn( lua_State *L )
{
   return lua_reg_metatable( L, synL_methods, "syn" );
}
static int lua_issyn( lua_State *L, int ind )
{
   return lua_is_foo( L, ind, SYN_METATABLE );
}
static synthesis_t* lua_tosyn( lua_State *L, int ind )
{
   return (synthesis_t*) lua_touserdata(L,ind);
}
static synthesis_t* luaL_checksyn( lua_State *L, int ind )
{
   if (lua_issyn(L,ind))
      return lua_tosyn(L,ind);
   luaL_typerror(L, ind, SYN_METATABLE);
   return NULL;
}
static synthesis_t* lua_pushsyn( lua_State *L, synthesis_t syn )
{
   synthesis_t *s;
   s = (synthesis_t*) lua_newuserdata( L, sizeof(synthesis_t) );
   assert( s != NULL );
   syn_copy( s, &syn );
   luaL_getmetatable(L, SYN_METATABLE);
   lua_setmetatable(L, -2);
   return s;
}
/* Exposed API. */
static int synL_gc( lua_State *L )
{
   synthesis_t *syn = luaL_checksyn(L,1);
   syn_free( syn );
   return 0;
}
static int synL_cmp( lua_State *L )
{
   synthesis_t *syn_a, *syn_b;
   syn_a = luaL_checksyn(L,1);
   syn_b = luaL_checksyn(L,2);
   lua_pushboolean(L, !syn_cmp( syn_a, syn_b ) );
   return 1;
}
static int synL_new( lua_State *L )
{
   synthesis_t syn;
   int num = luaL_checkint(L,1);
   syn_init( &syn );
   syn_set_frames( &syn, num );
   lua_pushsyn(L, syn);
   return 1;
}
static int synL_dup( lua_State *L )
{
   synthesis_t syn, *psyn;
   psyn = luaL_checksyn(L,1);
   syn_copy( &syn, psyn );
   lua_pushsyn(L, syn);
   return 1;
}
static int synL_addObj( lua_State *L )
{
   synthesis_t *syn = luaL_checksyn(L,1);
   kin_object_t *obj = luaL_checkobject(L,2);
   syn_object_add( syn, kin_obj_dup( obj ) );
   return 0;
}
static int synL_finalize( lua_State *L )
{
   synthesis_t *syn = luaL_checksyn(L,1);
   syn_finalize( syn );
   syn_calc( syn );
   lua_pushnumber( L, syn->n );
   lua_pushnumber( L, syn->m );
   lua_pushnumber( L, syn->ni );
   lua_pushnumber( L, syn->mi );
   return 4;
}
#define TBL_DBL( x, L, i, n ) \
lua_getfield(L,i,n); \
if (!lua_isnil(L,-1)) \
   x = luaL_checknumber( L, -1 ); \
lua_pop(L,1);
#define TBL_INT( x, L, i, n ) \
lua_getfield(L,i,n); \
if (!lua_isnil(L,-1)) \
   x = luaL_checkinteger( L, -1 ); \
lua_pop(L,1);
#define TBL_ULONG( x, L, i, n ) \
lua_getfield(L,i,n); \
if (!lua_isnil(L,-1)) \
   x = (unsigned long)luaL_checknumber( L, -1 ); \
lua_pop(L,1);
#define TBL_BOOL( x, L, i, n ) \
lua_getfield(L,i,n); \
if (!lua_isnil(L,-1)) \
   x = lua_toboolean( L, -1 ); \
lua_pop(L,1);
static int solver_parse_minpack( lua_State *L, minpack_options_t *opts, int ind )
{
   luaL_checktype(L,ind,LUA_TTABLE);
   TBL_DBL( opts->ftol, L, ind, "ftol" );
   TBL_DBL( opts->xtol, L, ind, "xtol" );
   TBL_DBL( opts->gtol, L, ind, "gtol" );
   TBL_INT( opts->maxfev, L, ind, "maxfev" );
   TBL_DBL( opts->epsfcn, L, ind, "epsfcn" );
   TBL_INT( opts->mode, L, ind, "mode" );
   TBL_DBL( opts->factor, L, ind, "factor" );
   return 0;
}
static int solver_parse_ga( lua_State *L, ga_options_t *opts, int ind )
{
   luaL_checktype(L,ind,LUA_TTABLE);

   TBL_INT( opts->verbose, L, ind, "verbose" );
   TBL_INT( opts->threads, L, ind, "threads" );
   TBL_BOOL( opts->converge, L, ind, "converge" );

   TBL_INT( opts->population, L, ind, "population" );
   TBL_INT( opts->generations, L, ind, "generations" );

   TBL_DBL( opts->eliteness, L, ind, "eliteness" );
   TBL_DBL( opts->crossover, L, ind, "crossover" );
   TBL_DBL( opts->mutation, L, ind, "mutation" );
   TBL_DBL( opts->seed_mul, L, ind, "seed_mul" );

   TBL_BOOL( opts->stop_sigusr1, L, ind, "stop_sigusr1" );
   TBL_BOOL( opts->stop_sigint, L, ind, "stop_sigint" );
   TBL_ULONG( opts->stop_elapsed, L, ind, "stop_elapsed" );
   TBL_DBL( opts->stop_fitness, L, ind, "stop_fitness" );
   
   TBL_BOOL( opts->sigfpe, L, ind, "sigfpe" );

   lua_getfield( L, ind, "minpack" );
   if (!lua_isnil(L,-1))
      solver_parse_minpack( L, &opts->minpack, -1 );
   lua_pop(L,1);

   return 0;
}
static int synL_print_info_minpack( lua_State *L, minpack_info_t *info )
{
   char *str;
   lua_newtable(L);
   lua_pushnumber(L,info->elapsed);
   lua_setfield(L,-2,"elapsed");
   lua_pushnumber(L,info->func_calls);
   lua_setfield(L,-2,"func_calls");
   lua_pushnumber(L,info->term_cond);
   lua_setfield(L,-2,"term_cond");
   lua_pushnumber(L,info->ftol);
   lua_setfield(L,-2,"ftol");
   lua_pushnumber(L,info->xtol);
   lua_setfield(L,-2,"xtol");
   lua_pushnumber(L,info->gtol);
   lua_setfield(L,-2,"gtol");
   str = minpack_term_string( info );
   lua_pushstring(L,str);
   lua_setfield(L,-2,"term_str");
   free(str);
   return 1;
}
#ifdef HAVE_MINPACK
static int synL_solve_fminpack( lua_State *L )
{
   minpack_options_t opts;
   minpack_info_t info;
   synthesis_t *syn = luaL_checksyn(L,1);

   /* Parse options. */
   minpack_options_default( &opts );
   if (lua_gettop(L) > 1)
      solver_parse_minpack( L, &opts, 2 );

   /* Solve. */
   syn_solve_fminpack( syn, &opts, &info );

   /* Push info. */
   return synL_print_info_minpack( L, &info );;
}
#endif /* HAVE_MINPACK */
static int synL_solve_minpack( lua_State *L )
{
   minpack_options_t opts;
   minpack_info_t info;
   synthesis_t *syn = luaL_checksyn(L,1);

   /* Parse options. */
   minpack_options_default( &opts );
   if (lua_gettop(L) > 1)
      solver_parse_minpack( L, &opts, 2 );

   /* Solve. */
   syn_solve_minpack( syn, &opts, &info );

   /* Push info. */
   return synL_print_info_minpack( L, &info );;
}
static int synL_solve_ga( lua_State *L )
{
   ga_options_t opts;
   ga_info_t info;
   synthesis_t *syn = luaL_checksyn(L,1);

   /* Parse options. */
   ga_options_default( &opts );
   if (lua_gettop(L) > 1)
      solver_parse_ga( L, &opts, 2 );

   /* Solve. */
   syn_solve_ga( syn, &opts, &info );

   /* Push info. */
   lua_newtable(L);
   lua_pushnumber(L,info.elapsed);
   lua_setfield(L,-2,"elapsed");
   lua_pushnumber(L,info.generations);
   lua_setfield(L,-2,"generations");
   lua_pushnumber(L,info.max_generations);
   lua_setfield(L,-2,"max_generations");
   lua_pushnumber(L,info.fit_best);
   lua_setfield(L,-2,"fit_best");
   lua_pushnumber(L,info.fit_mean);
   lua_setfield(L,-2,"fit_mean");
   lua_pushnumber(L,info.fit_stddev);
   lua_setfield(L,-2,"fit_stddev");
   return 1;
}

#ifdef HAVE_NLOPT
static int synL_solve_nlopt( lua_State *L )
{
   nlopt_options_t opts;
   nlopt_info_t info;
   synthesis_t *syn = luaL_checksyn(L,1);

   /* Parse options. */
   nlopt_options_default( &opts );

   /* Solve. */
   syn_solve_nlopt( syn, &opts, &info );

   /* Push info. */
   return 0;
}
#endif /* HAVE_NLOPT */

static int solver_parse_cmaes( lua_State *L, cmaes_options_t *opts, int ind )
{
   luaL_checktype(L,ind,LUA_TTABLE);
   TBL_BOOL( opts->converge, L, ind, "converge" );
   TBL_DBL( opts->lambda, L, ind, "lambda" );
   TBL_DBL( opts->stop_fitness, L, ind, "stop_fitness" );
   TBL_ULONG( opts->stop_evals, L, ind, "stop_evals" );
   TBL_ULONG( opts->stop_iter, L, ind, "stop_iter" );
   TBL_ULONG( opts->stop_elapsed, L, ind, "stop_elapsed" );
   return 0;
}
static int synL_solve_cmaes( lua_State *L )
{
   cmaes_options_t opts;
   cmaes_info_t info;
   synthesis_t *syn = luaL_checksyn(L,1);

   /* Parse options. */
   cmaes_options_default( &opts );
   if (lua_gettop(L) > 1)
      solver_parse_cmaes( L, &opts, 2 );

   /* Solve. */
   syn_solve_cmaes( syn, &opts, &info );

   /* Push info. */
   lua_newtable(L);
   lua_pushnumber(L,info.elapsed);
   lua_setfield(L,-2,"elapsed");
   lua_pushnumber(L,info.iterations);
   lua_setfield(L,-2,"iterations");
   lua_pushnumber(L,info.minf);
   lua_setfield(L,-2,"fit_best");
   return 1;
}

static int synL_print( lua_State *L )
{
   synthesis_t *syn = luaL_checksyn(L,1);
   syn_printDetail( syn );
   return 0;
}

static int synL_printClaim( lua_State *L )
{
   synthesis_t *syn = luaL_checksyn(L,1);
   syn_printClaim( syn );
   return 0;
}

static int synL_printJacobian( lua_State *L )
{
   synthesis_t *syn = luaL_checksyn(L,1);
   double step = 1e-3;
   if (lua_gettop(L) > 1)
      step = luaL_checknumber(L,2);
   syn_printJacobian( syn, step );
   return 0;
}


static int synL_save( lua_State *L )
{
   synthesis_t *syn = luaL_checksyn(L,1);
   const char *str = luaL_checkstring(L,2);
   syn_save( syn, str );
   return 0;
}

static int synL_raw_x( lua_State *L )
{
   int i;
   synthesis_t *syn = luaL_checksyn(L,1);
   lua_newtable( L );
   lua_newtable( L );
   lua_newtable( L );
   for (i=0; i<syn->n; i++) {
      lua_pushnumber( L, i+1 );
      lua_pushnumber( L, syn->x[i] );
      lua_settable( L, -5 );
      lua_pushnumber( L, i+1 );
      lua_pushnumber( L, syn->lb[i] );
      lua_settable( L, -4 );
      lua_pushnumber( L, i+1 );
      lua_pushnumber( L, syn->ub[i] );
      lua_settable( L, -3 );
   }
   return 3;
}

static int synL_raw_fvec( lua_State *L )
{
   int i;
   synthesis_t *syn = luaL_checksyn(L,1);
   lua_newtable( L );
   for (i=0; i<syn->m; i++) {
      lua_pushnumber( L, i+1 );
      lua_pushnumber( L, syn->fvec[i] );
      lua_settable( L, -3 );
   }
   return 1;
}

static int synL_visualize( lua_State *L )
{
   synthesis_t *syn_a, *syn_b;
   syn_a = luaL_checksyn(L,1);
   if (lua_gettop(L) > 1)
      syn_b = luaL_checksyn(L,2);
   else
      syn_b = NULL;
   visualize( syn_a, syn_b );
   return 0;
}

#define PUSHN( L, n, s ) \
lua_pushnumber( L, n ); \
lua_setfield( L, -2, s )
static int synL_stats( lua_State *L )
{
   synthesis_t *syn = luaL_checksyn(L,1);
   lua_newtable( L );
   PUSHN( L, syn->finalized,  "finalized" );
   if (!syn->finalized)
      return 1;
   PUSHN( L, syn->n,          "n" );
   PUSHN( L, syn->m,          "m" );
   PUSHN( L, syn->ni,         "ni" );
   PUSHN( L, syn->mi,         "mi" );
   PUSHN( L, syn->L,          "L" );
   PUSHN( L, syn->njoints,    "r" );
   PUSHN( L, syn->nbranches,  "b" );
   return 1;
}
#undef PUSHN




