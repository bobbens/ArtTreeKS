

#include "solver.h"

#include <fenv.h>
#include <math.h>
#include <fcntl.h>
#include <assert.h>
#include <signal.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <pthread.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/types.h>

#include <dq_vec3.h>

#include "synthesis.h"


/** Conditional logging. */
#define LOG(v, str, args...)  \
   if (v) { fprintf( stdout, str, ## args ); fflush( stdout ); }


#define MIN(a,b)           (((a)<(b))?(a):(b))     /**< Your standard minimum macro. */
#define MAX(a,b)           (((a)>(b))?(a):(b))     /**< Your standard maximum macro. */
#define CLAMP(min,max,x)   (MAX(MIN(max,x),min))   /**< Your standard clamping macro. */


#define FPE_FLAGS    (FE_DIVBYZERO | FE_INVALID)


/* Forward declarations. */
struct ga_population_s;
typedef struct ga_population_s ga_population_t;


/**
 * @brief Genetic algorithm entityt.
 */
typedef struct ga_entity_s {
   int id;           /**< ID for debugging. */
   /* ---VERY IMPORTANT---
    * The synthesis branches should always be up to date and we should refresh
    * the vector as needed. This is to keep sanity. */
   synthesis_t *syn; /**< Synthesis. */
   double fitness;   /**< Fitness rating for the entity. */
} ga_entity_t;


/**
 * @brief Threading data structure.
 */
typedef struct synthesis_pthread_s {
   /* Local information. */
   pthread_t         thread;  /**< Thread id. */
   /* Parent information. */
   ga_population_t*  pop;     /**< Population thread is running on. */
   int*              processed; /**< Number of population processed. */
   pthread_mutex_t*  lock;    /**< Lock tied to global conditional. */
   /* Current processing data. */
   const ga_options_t* opts;  /**< Options for the genetic algorithm. */
   void (*func)( ga_entity_t*, const synthesis_t*, const ga_options_t* ); /**< Function to run. */
} synthesis_pthread_t;


/**
 * @brief Genetic algorithm population.
 */
struct ga_population_s {
   synthesis_t *parent;       /**< Parent synthesis. */
   int pop_target;            /**< Target population. */
   int pop;                   /**< Current size of the population. */
   ga_entity_t *entities;     /**< Entities in the population. */
   int alloced;               /**< Allocated memory. */
   double fitness_total;      /**< Sum of total fitness for roulette. */
   synthesis_pthread_t *th;   /**< Threading information. */
   int nth;                   /**< Number of threads. */
   pthread_mutex_t lock;      /**< Thread lock. */
};


/*
 * Prototypes.
 */
/* Randomness. */
static void rand_init (void);
static void rand_exit (void);
static int rand_bool (void);
static int rand_int_range( int low, int high );
static double rand_double (void);
static double rand_double_range( double low, double high );
static double rand_double_exponential( double lambda );
static double rand_double_gaussian (void);
/* Entities. */
static void ga_ent_evaluate( ga_entity_t *ent, const ga_options_t *opts );
static void ga_ent_seed( ga_entity_t *ent, const synthesis_t *syn, const ga_options_t *opts );
static void ga_ent_mutate( ga_entity_t *ent, const synthesis_t *parent, const ga_options_t *opts );
static void ga_ent_crossover( ga_entity_t *son, ga_entity_t *daughter,
      const ga_entity_t *father, const ga_entity_t *mother, const ga_options_t *opts );
static int ga_ent_comp( const void *p1, const void *p2 );
static void ga_ent_cleanup( ga_entity_t *ent );
/* Population memory. */
static int ga_pop_grow( ga_population_t *pop );
static int ga_pop_grow_init( ga_population_t *pop, const ga_entity_t *parent );
static void ga_pop_cleanup( ga_population_t *pop );
/* Population misc. */
static ga_entity_t* ga_pop_select_roulette( const ga_population_t *pop );
static int ga_pop_seed( ga_population_t *pop, const ga_options_t *opts );
static int ga_pop_eliteness( ga_population_t *pop, const ga_population_t *old,
      const ga_options_t *opts );
static int ga_pop_crossover( ga_population_t *pop, const ga_population_t *old,
      const ga_options_t *opts );
static int ga_pop_mutate( ga_population_t *pop, const ga_options_t *opts );
static int ga_pop_converge( ga_population_t *pop, const ga_options_t *opts );
static int ga_pop_evaluate( ga_population_t *pop, const ga_options_t *opts, ga_info_t *info );
static int ga_pop_sort( ga_population_t *pop, const ga_options_t *opts );
/* Threading. */
static void ga_ent_evaluate_pthread( ga_entity_t *ent, const synthesis_t *parent, const ga_options_t *opts );
static void ga_ent_seed_pthread( ga_entity_t *ent, const synthesis_t *parent, const ga_options_t *opts );
static void ga_ent_converge_pthread( ga_entity_t *ent, const synthesis_t *parent, const ga_options_t *opts );
static void* ga_ent_pthread( void *data );
static int ga_pop_pthread_init( ga_population_t *pop, const ga_options_t *opts );
static void ga_pop_pthread_free( ga_population_t *pop, const ga_options_t *opts );
static int ga_pop_pthread( ga_population_t *pop, const ga_options_t *opts,
      void (*func)( ga_entity_t*, const synthesis_t*, const ga_options_t* ) );
/* Misc. */
static void ga_signal_handler( int sig, siginfo_t *info, void *unused );
static void ga_signal_fpe( int sig, siginfo_t *info, void *unused );


/*
 * Globals.
 */
static int rand_fd   = -1; /**< /dev/urandom */
static int ga_finish = 0; /**< Early finish? */


/**
 * @brief Initializes random subsystem.
 *
 * @TODO Possibly make it thread safe or use a different approach.
 */
static void rand_init (void)
{
   if (rand_fd < 0)
      rand_fd = open( "/dev/urandom", O_RDONLY );
   assert( rand_fd >= 0 );

   unsigned int seed;
   ssize_t ret = read( rand_fd, &seed, sizeof(seed) );
   assert( ret == sizeof(seed) );
   srand( seed );
}
/**
 * @brief Exits random subsystem.
 */
static void rand_exit (void)
{
   if (rand_fd >= 0)
      close( rand_fd );
   rand_fd = -1;
}
/**
 * @brief Gets a random boolean.
 */
static int rand_bool (void)
{
   return (rand_double() > 0.5);
}
/**
 * @brief Gets an int in range [low,high].
 */
static int rand_int_range( int low, int high )
{
   return (int)round( rand_double_range( (double)low, (double)high ) );
}
/**
 * @brief Gets a random double in the [0,1] range.
 */
static double rand_double (void)
{
#ifdef RAND_ACCURATE
   ssize_t ret;
   uint64_t r;
   ret = read( rand_fd, &r, sizeof(r) );
   assert( ret == sizeof(r) );
   return (double)r / (double)UINT64_MAX;
#else /* RAND_ACCURATE */
   return (double)rand() / (double)RAND_MAX;
#endif /* RAND_ACCURATE */
}
/**
 * @brief Gets a random double in the range [low,high].
 */
static double rand_double_range( double low, double high )
{
   return (low + (high-low)*rand_double());
}
/**
 * @brief Gets a random double following the exponential function for a given lambda.
 */
static double rand_double_exponential( double lambda )
{
   return (-1./lambda) * log(rand_double());
}
/**
 * @brief Gets a random double that follows the gaussian distribution.
 */
static double rand_double_gaussian (void)
{
   double theta, rsq;
   theta = rand_double_range( 0., 2.*M_PI );
   rsq   = rand_double_exponential( 0.5 );
   return (sqrt(rsq) * cos(theta));
}


/**
 * @brief Does statistics on a set of data.
 *
 *    @param[out] mean Mean of set of data.
 *    @param[out] stddev Standard deviation of set of data.
 *    @param[out] min Minimum data value.
 *    @param[out] max Maximum data value.
 *    @param[in] data Set of data to do statistics on.
 *    @param[in] n Length of set of data.
 */
static void rand_stats( double *mean, double *stddev,
      double *min, double *max, const double *data, int n )
{
   int i;
   double acc, m, M;

   /* First loop, calculate mean, min and max. */
   m   = INFINITY;
   M   = 0.;
   acc = 0.;
   for (i=0; i<n; i++) {
      acc  += data[i];
      if (data[i] < m)
         m = data[i];
      if (data[i] > M)
         M = data[i];
   }
   acc /= (double) n;
   *mean = acc;
   if (min != NULL)
      *min  = m;
   if (max != NULL)
      *max  = M;

   /* Calculate standard deviation. */
   acc = 0.;
   for (i=0; i<n; i++)
      acc += pow( *mean - data[i], 2. );
   acc = sqrt( acc / (double)n );
   *stddev = acc;
}


/**
 * @brief Sets the options to the default.
 *    @param[out] opts Options to set.
 */
void ga_options_default( ga_options_t *opts )
{
   assert( opts != NULL );
   memset( opts, 0, sizeof(ga_options_t) );
   opts->verbose     = 0;
   opts->converge    = 1;
   opts->threads     = 0;
   /* GA. */
   opts->population  = 100;
   opts->generations = 100;
   opts->eliteness   = 0.05;
   opts->crossover   = 0.1;
   opts->mutation    = 0.5;
   opts->seed_mul    = 1.0;
   /* Stop conditions. */
   opts->stop_sigusr1 = 1;
   opts->stop_sigint = 0;
   opts->stop_elapsed = 0.;
   opts->stop_fitness = 1e10;
   opts->sigfpe      = 0;
   /* MINPACK. */
   minpack_options_default( &opts->minpack );
}


/**
 * @brief Evaluates the entity to set it's fitness parameter.
 *    @param ent Entity to evaluate.
 */
static void ga_ent_evaluate( ga_entity_t *ent, const ga_options_t *opts )
{
   int i;
   double fit;
   double v;

   v = opts->verbose;

   /* Calculate the synthesis. */
   syn_calc( ent->syn );

   /* Calculate the sum of fitness. */
   fit = 0.;
   for (i=0; i<ent->syn->m; i++)
      fit += fabs( ent->syn->fvec[i] );
   ent->fitness = 1./fit;

   /* Debugging. */
   if (v>9) {
      LOG( v>9, "-Evaluated entity-%03d----------\n", ent->id );
      syn_printDetail( ent->syn );
   }
}


/**
 * @brief Normalizes a plucker coordinate.
 */
static void plucker_normalize( plucker_t *P )
{
   double c[3];
   vec3_normalize( P->s );
   vec3_cross( c, P->s, P->s0 );
   vec3_cross( P->s0, c, P->s );
}


/**
 * @brief Seeds an entity in the genetic algorithm.
 *    @param ent Entity to seed.
 *    @param syn Base synthesis to seed of off.
 *    @param opts Options to use.
 */
static void ga_ent_seed( ga_entity_t *ent, const synthesis_t *syn, const ga_options_t *opts )
{
   int i, j;
   kin_joint_t *kj;
   double v;

   v = opts->verbose;

   /* Create the synthesis. */
   ent->syn = syn_dup( syn );

   assert( ent->syn->njoints == syn->njoints );

   /* Mutate everything within bounds. */
   for (i=0; i<syn->n; i++)
      ent->syn->x[i] = rand_double_range( ent->syn->lb[i], ent->syn->ub[i] );

   /* Now the new crossed X vectors must be put back on the claims. */
   syn_map_from_x( ent->syn, ent->syn->x, ent->syn->n );

   /* Set plucker conditions. */
   for (i=0; i<syn->njoints; i++) {
      kj    = ent->syn->joints[i];
      LOG( v>2&&v<=9, "   seed joint %d\n", i );
      if (kj->claim_S != NULL) {
         plucker_normalize( &kj->S );
         LOG( v>2&&v<=9, "      plucker = %.3e, %.3e, %.3e x %.3e, %.3e, %.3e (dot=%.3e)\n",
               kj->S.s[0], kj->S.s[1], kj->S.s[2],
               kj->S.s0[0], kj->S.s0[1], kj->S.s0[2],
               vec3_dot( kj->S.s, kj->S.s0 ) );
      }
      if (kj->pos.claim != NULL) {
         mm_updateMask( &kj->pos.values );
         for (j=0; j<kj->pos.nvalues; j++) {
            LOG( v>2&&v<=9, "      pos[%d] = %.3e\n", j, ((double*)kj->pos.values.mask_vec)[j] );
         }
      }
   }

   /* Debug. */
   LOG( v>9, "-Seeded entity-%03d-------------\n", ent->id );
   if (v>9)
      syn_printDetail( ent->syn );
}


/**
 * @brief Mutate a single entity.
 *    @param ent Entity to mutate.
 *    @param parent Parent to base off of.
 *    @param opts Options to use.
 */
static void ga_ent_mutate( ga_entity_t *ent, const synthesis_t *parent, const ga_options_t *opts )
{
   double v;

   v = opts->verbose;

#if 1
   int i, n, nn;
   (void) parent;
   n = rand_int_range( 1, (ent->syn->n / 64) ); /* Mutations to do. */
   for (nn=0; nn<n; nn++) {
      kin_claim_t *c;
      int r = rand_int_range( 0, ent->syn->n-1 );
      int p = 0;
      for (c=ent->syn->claim_x; c!=NULL; c=c->next) {
         int o = r-p;
         if (o-(int)c->size < 0) {
            c->dat[o] = rand_double_range( c->lb[o], c->ub[o] );
            break;
         }
         p += c->size;
      }
   }

   /* We want to renormalize all the plucker coordinates, just in case. */
   for (i=0; i<ent->syn->njoints; i++)
      plucker_normalize( &ent->syn->joints[i]->S );
#else
   int i, j, n, rj, ra;
   plucker_t *pl, *pl_lb, *pl_ub;
   double v, lb, ub;

   /* Mutate angles. */
   if (rand_bool()) {
      n = rand_int_range( 1, (ent->syn->n / 64) );
      for (i=0; i<n; i++) {
         /* Get a random mutatable joint. */
         j = 0;
         do {
            if (j>10000)
               return;
            rj = rand_int_range( 0, ent->syn->njoints-1 );
            j++;
         } while (ent->syn->joints[rj]->pos.claim == NULL);
         /* Random position. */
         if (ent->syn->joints[rj]->pos.claim != NULL) {
            ra = rand_int_range( 0, ent->syn->L-2 ); /* It's from 0 to L-1 */
            /* Calculate bounds. */
            lb    = ((double*)parent->joints[rj]->pos.values_lb.mask_vec)[ra];
            ub    = ((double*)parent->joints[rj]->pos.values_ub.mask_vec)[ra];
            ((double*)ent->syn->joints[rj]->pos.values.mask_vec)[ra] = rand_double_range( lb, ub );
            mm_updateMap( &ent->syn->joints[rj]->pos.values );
         }
      }
   }
   /* Mutate plucker coords. */
   else {
      //if (rand_bool()) {
      if (1) {
         /* Get a random mutatable joint. */
         i = 0;
         do {
            if (i>10000)
               return;
            rj = rand_int_range( 0, ent->syn->njoints-1 );
            i++;
         } while (ent->syn->joints[rj]->claim_S == NULL);
         pl    = &ent->syn->joints[rj]->S;
         pl_lb = &parent->joints[rj]->S_lb;
         pl_ub = &parent->joints[rj]->S_ub;

         /* Choose mutate target. */
         c = rand_int_range( 0, 2 );
         /* Mutate real component. */
         if ((c==0) || (c==1)) {
            /* Mutate a single component. */
            if (rand_bool()) {
               i = rand_int_range( 0, 2 );
               pl->s[i] = rand_double_range( pl_lb->s[i], pl_ub->s[i] );
            }
            /* Mutate all three components. */
            else {
               for (i=0; i<3; i++) {
                  pl->s[i] = rand_double_range( pl_lb->s[i], pl_ub->s[i] );
               }
            }
         }
         /* Mutate dual component. */
         if ((c==0) || (c==2)) {
            for (i=0; i<3; i++) {
               pl->s0[i] = rand_double_range( pl_lb->s0[i], pl_ub->s0[i] );
            }
         }
         /* Normalize the vector. */
         plucker_normalize( pl );
      }
   }
#endif

   /* Debug. */
   LOG( v>9, "-Mutated entity-%03d------------\n", ent->id );
   if (v>9)
      syn_printDetail( ent->syn );
}


/**
 * @brief Crosses over a chunk.
 */
static void ga_ent_crossover_chunk( void *a, void *b,
      const void *c, const void *d, size_t size )
{
   if (rand_bool()) {
      memcpy( a, c, size );
      memcpy( b, d, size );
   }
   else {
      memcpy( a, d, size );
      memcpy( b, c, size );
   }
}


/**
 * @brief Randomly chooses an entity from two.
 */
static void ga_ent_crossover_swap( ga_entity_t **eqa, ga_entity_t **eqb,
      ga_entity_t *eq_in_a, ga_entity_t *eq_in_b )
{
   if (rand_bool()) {
      *eqa = eq_in_a;
      *eqb = eq_in_b;
   }
   else {
      *eqa = eq_in_b;
      *eqb = eq_in_a;
   }
}


/**
 * @brief Produces a son and a daughter from a father and a mother.
 */
static void ga_ent_crossover( ga_entity_t *son, ga_entity_t *daughter,
      const ga_entity_t *father, const ga_entity_t *mother, const ga_options_t *opts )
{
   double v;

   v = opts->verbose;

#if 1
   kin_claim_t *c;
   size_t p;

   /* We must map the X claims to the X vector to be able to manipulate them. */
   syn_map_to_x( father->syn, NULL, NULL, father->syn->x );
   syn_map_to_x( mother->syn, NULL, NULL, mother->syn->x );

   /* Here we can map the claims using the father as a reference. */
   p = 0;
   for (c=father->syn->claim_x; c!=NULL; c=c->next) {
      ga_ent_crossover_chunk( &son->syn->x[p], &daughter->syn->x[p],
            &father->syn->x[p], &mother->syn->x[p], c->size*sizeof(double) );
      p += c->size;
   }

   /* Now the new crossed X vectors must be put back on the claims. */
   syn_map_from_x( son->syn,      son->syn->x,      son->syn->n );
   syn_map_from_x( daughter->syn, daughter->syn->x, daughter->syn->n );

   /* Sanity checks. */
   assert( (int)p == father->syn->n );
   assert( (int)p == mother->syn->n );
   assert( (int)p == son->syn->n );
   assert( (int)p == daughter->syn->n );
#else
   ga_entity_t *eqa, *eqb;
   int i, c, n;
   /* The old genetic algorithm would treat the plucker coordinates and the
    * joints as separate chromosomes and chop them up at a random point in
    * the vector. */
   /* Plucker coordinates. */
   c = rand_int_range( 0, father->syn->njoints-1 );
   ga_ent_crossover_swap( &eqa, &eqb, son, daughter );
   for (i=0; i<c; i++) {
      memcpy( &eqa->syn->joints[i]->S, &father->syn->joints[i]->S, sizeof(plucker_t) );
      memcpy( &eqb->syn->joints[i]->S, &mother->syn->joints[i]->S, sizeof(plucker_t) );
   }
   if (rand_bool()) {
      memcpy( eqa->syn->joints[i]->S.s0, mother->syn->joints[i]->S.s0, sizeof(double)*3 );
      memcpy( eqb->syn->joints[i]->S.s0, father->syn->joints[i]->S.s0, sizeof(double)*3 );
   }
   for (   ; i<father->syn->njoints; i++) {
      memcpy( &eqa->syn->joints[i]->S, &mother->syn->joints[i]->S, sizeof(plucker_t) );
      memcpy( &eqb->syn->joints[i]->S, &father->syn->joints[i]->S, sizeof(plucker_t) );
   }
   /* Angles. */
   n = 0;
   for (i=0; i<father->syn->njoints; i++)
      n += father->syn->joints[i]->npos;
   ga_ent_crossover_swap( &eqa, &eqb, son, daughter );
   c = rand_int_range( 0, n-1 );
   n = 0;
   for (i=0; i<father->syn->njoints; i++) {
      if (c < n) {
         memcpy( eqa->syn->joints[i]->pos, father->syn->joints[i]->pos,
               sizeof(double)*father->syn->joints[i]->npos );
         memcpy( eqb->syn->joints[i]->pos, mother->syn->joints[i]->pos,
               sizeof(double)*father->syn->joints[i]->npos );
      }
      else if (c-n < father->syn->joints[i]->npos) {
         memcpy( eqa->syn->joints[i]->pos, father->syn->joints[i]->pos,
               sizeof(double)*(c-n) );
         memcpy( eqb->syn->joints[i]->pos, mother->syn->joints[i]->pos,
               sizeof(double)*(c-n) );
         memcpy( eqa->syn->joints[i]->pos, mother->syn->joints[i]->pos,
               sizeof(double)*(father->syn->joints[i]->npos - (c-n)) );
         memcpy( eqb->syn->joints[i]->pos, father->syn->joints[i]->pos,
               sizeof(double)*(father->syn->joints[i]->npos - (c-n)) );
      }
      else {
         memcpy( eqa->syn->joints[i]->pos, mother->syn->joints[i]->pos,
               sizeof(double)*father->syn->joints[i]->npos );
         memcpy( eqb->syn->joints[i]->pos, father->syn->joints[i]->pos,
               sizeof(double)*father->syn->joints[i]->npos );
      }
      n += father->syn->joints[i]->npos;
   }
#endif

   /* Debug. */
   if (v>9) {
      LOG( v>9, "-Crossover entity--------------\n" );
      LOG( v<9, "Father: %03d\n", father->id);
      syn_printDetail( father->syn );
      LOG( v<9, "Mother: %03d\n", mother->id);
      syn_printDetail( mother->syn );
      LOG( v<9, "Son: %03d\n", son->id);
      syn_printDetail( son->syn );
      LOG( v<9, "Daughter: %03d\n", daughter->id);
      syn_printDetail( daughter->syn );
   }
}


/**
 * @brief Frees an entity.
 *    @param ent Entity to free.
 */
static void ga_ent_cleanup( ga_entity_t *ent )
{
   if (ent->syn != NULL)
      syn_destroy( ent->syn );
   memset( ent, 0, sizeof(ga_entity_t) );
}


/**
 * @brief Grows a population by one member.
 *    @param pop Population to grow.
 *    @return Offset of the new entity in the population.
 */
static int ga_pop_grow( ga_population_t *pop )
{
   int id;

   pop->pop++;
   if (pop->pop <= pop->alloced) {
      id = pop->pop-1;
      pop->entities[id].id = id;
      return id;
   }

   /* Increase allocation. */
   assert( pop->alloced == 0 );
   assert( pop->pop_target > 0 );
   pop->alloced = pop->pop_target;
   pop->entities = calloc( (size_t)pop->alloced, sizeof(ga_entity_t) );

   id = pop->pop-1;
   pop->entities[id].id = id;
   return id;
}


/**
 * @brief Grows a population by one member copied from another population.
 *    @param pop Population to grow.
 *    @param parent Parent to copy from.
 */
static int ga_pop_grow_init( ga_population_t *pop, const ga_entity_t *parent )
{
   ga_entity_t *ent;
   int off;
   off      = ga_pop_grow( pop );
   ent      = &pop->entities[off];
   memset( ent, 0, sizeof(ga_entity_t) );
   ent->syn = syn_dup( parent->syn );
   ent->id  = off;
   return off;
}


/**
 * @brief Frees a population.
 *    @param pop Population to free.
 */
static void ga_pop_cleanup( ga_population_t *pop )
{
   int i;
   for (i=0; i<pop->pop; i++)
      ga_ent_cleanup( &pop->entities[i] );
   free( pop->entities );
   pop->entities  = 0;
   pop->pop       = 0;
   pop->alloced   = 0;
   pop->fitness_total = 0.;
}


/**
 * @brief Selects an entity based on roulette.
 *    @param pop Population to select entity from.
 *    @return Entity selected from the population.
 */
static ga_entity_t* ga_pop_select_roulette( const ga_population_t *pop )
{
   double f, r;
   int i;

   r = rand_double_range( 0., pop->fitness_total );
   f = 0.;
   for (i=0; i<pop->pop; i++) {
      f += pop->entities[i].fitness;
      if (r < f)
         return &pop->entities[i];
   }

   assert( "Roulette failed." == NULL );
   return NULL;
}

/**
 * @brief Evaluation wrapper for pthreads.
 */
static void ga_ent_evaluate_pthread( ga_entity_t *ent, const synthesis_t *parent, const ga_options_t *opts )
{
   (void) parent;
   ga_ent_evaluate( ent, opts );
}
/**
 * @brief Seed wrapper for pthreads.
 */
static void ga_ent_seed_pthread( ga_entity_t *ent, const synthesis_t *parent, const ga_options_t *opts )
{
   ga_ent_seed( ent, parent, opts );
}
#if 0
static void ga_ent_mutate_pthread( ga_entity_t *ent, const synthesis_t *parent, const ga_options_t *opts )
{
   ga_ent_mutate( ent, parent, opts );
}
#endif
/**
 * @brief Convergence wrapper for pthreads.
 *
 * Slowest part due to the QR decomposition of the Levenberg-Marquadt algorithm.
 */
static void ga_ent_converge_pthread( ga_entity_t *ent, const synthesis_t *parent,
      const ga_options_t *opts )
{
   (void) parent;
   minpack_info_t info;
   double v = opts->verbose;
   /* We must map to X so the minpack solver works and then map back. */
   syn_map_to_x(   ent->syn, NULL, NULL, ent->syn->x );
   syn_solve_minpack(  ent->syn, &opts->minpack, &info );
   syn_map_from_x( ent->syn, ent->syn->x, ent->syn->n );

   /* Check for possible errors. */
   if (v>1) {
      if (info.term_cond == 5) {
         LOG( v>1, "Entity %03d exceeded function calls (limit %d)\n", ent->id, info.func_calls );
      }
   }

   /* Conditional logging when we have high verbosity. */
   if (v>9) {
      LOG( v>9, "-Converged entity-%03d----------\n", ent->id );
      syn_printDetail( ent->syn );
   }
}


/**
 * @brief Pthread runner.
 *
 * This is the heart of threading. These spawn workers that attempt to handle all
 *  the data.
 *
 *    @param data Pthread synthesis data.
 */
static void* ga_ent_pthread( void *data )
{
   synthesis_pthread_t *syn_pt;
   ga_entity_t *ent;

   /* Initialize. */
   syn_pt = (synthesis_pthread_t*) data;

   /* Lock and signal done. */
   while (1) {
      pthread_mutex_lock( syn_pt->lock );

      /* See if we're done. */
      if (*syn_pt->processed >= syn_pt->pop->pop)
         break;

      /* Choose entity to process. */
      ent = &syn_pt->pop->entities[ (*syn_pt->processed) ];
      (*syn_pt->processed)++; /* Processed more. */

      pthread_mutex_unlock( syn_pt->lock );

      /* Work - only unlocked bit. */
      syn_pt->func( ent, syn_pt->pop->parent, syn_pt->opts );
   }

   pthread_mutex_unlock( syn_pt->lock );

   /* Clean up. */
   pthread_exit( NULL );
}


/**
 * @brief Initializes a pthread data.
 */
static int ga_pop_pthread_init( ga_population_t *pop, const ga_options_t *opts )
{
   if (opts->threads <= 1)
      return 0;

   pthread_mutex_init( &pop->lock, NULL );
   pop->th  = calloc( (size_t)opts->threads, sizeof(synthesis_pthread_t) );
   assert( pop->th != NULL );
   pop->nth = opts->threads;

   return 0;
}


/**
 * @brief Frees a pthread data.
 */
static void ga_pop_pthread_free( ga_population_t *pop, const ga_options_t *opts )
{
   if (opts->threads <= 1)
      return;

   /* Clean up. */
   free( pop->th );
   pop->th  = NULL;
   pop->nth = 0;
   pthread_mutex_destroy( &pop->lock );
}


/**
 * @brief Converges the buggers.
 */
static int ga_pop_pthread( ga_population_t *pop, const ga_options_t *opts,
      void (*func)( ga_entity_t*, const synthesis_t*, const ga_options_t* ) )
{
   int i, process;
   synthesis_pthread_t *th;
   double v;

   v = opts->verbose;

   /* No threads. */
   if (opts->threads <= 1) {
      LOG( v>9, "----------------------------------------------------\n" );
      for (i=0; i<pop->pop; i++) {
         func( &pop->entities[i], pop->parent, opts );
      }
      return 0;
   }

   /* Some global stuff. */
   process  = 0;

   /* Init and spawn. */
   pthread_mutex_lock( &pop->lock );
   for (i=0; i<opts->threads; i++) {
      th             = &pop->th[i];
      th->lock       = &pop->lock;
      th->processed  = &process;
      th->pop        = pop;
      th->func       = func;
      th->opts       = opts;

      /* Spawn thread. */
      pthread_create( &th->thread, NULL, ga_ent_pthread, (void*) th );
   }
   /* Wait until they're done. */
   pthread_mutex_unlock( &pop->lock );

   /* Wait until done. */
   for (i=0; i<opts->threads; i++) {
      th             = &pop->th[i];
      pthread_join( th->thread, NULL );
   }

   return 0;
}



/**
 * @brief Seeds a population.
 */
static int ga_pop_seed( ga_population_t *pop, const ga_options_t *opts )
{
   int i;

   assert( pop->entities == NULL );

   /* Allocate and seed entities. */
   pop->pop       = (int)round((double)opts->population * opts->seed_mul);
   pop->entities  = calloc( (size_t)pop->pop, sizeof(ga_entity_t) );
   assert( pop->entities != NULL );
   for (i=0; i<pop->pop; i++)
      pop->entities[i].id = i;

   /* Run threaded. */
   return ga_pop_pthread( pop, opts, ga_ent_seed_pthread );
}


/**
 * @brief Runs eliteness.
 */
static int ga_pop_eliteness( ga_population_t *pop, const ga_population_t *old,
      const ga_options_t *opts )
{
   int i, n, p;

   p = MIN( old->pop, opts->population );
   n = (int)round( ((double)p)*opts->eliteness );
   for (i=0; i<n; i++)
      ga_pop_grow_init( pop, &old->entities[i] );
   return 0;
}


/**
 * @brief Crosses over a population.
 */
static int ga_pop_crossover( ga_population_t *pop, const ga_population_t *old,
      const ga_options_t *opts )
{
   int i, n, p, id_son, id_daughter;
   ga_entity_t *son, *daughter, *father, *mother;
   double v;

   v = opts->verbose;

   /* Crossover as needed. */
   p = MIN( old->pop, opts->population );
   n = (int)round( ((double)p)*opts->crossover/2. );
   for (i=0; i<n; i++) {
      father      = ga_pop_select_roulette( old );
      mother      = ga_pop_select_roulette( old );
      id_son      = ga_pop_grow_init( pop, father );
      id_daughter = ga_pop_grow_init( pop, mother );
      son         = &pop->entities[ id_son ];
      daughter    = &pop->entities[ id_daughter ];
      LOG( v>3, "Crossing %d and %d\n", father->id, mother->id );
      ga_ent_crossover( son, daughter, father, mother, opts );
   }

   /* Fill rest of population. */
   n = opts->population - pop->pop;
   for (i=0; i<n; i++) {
      father      = ga_pop_select_roulette( old );
      ga_pop_grow_init( pop, father );
   }

   return 0;
}


/**
 * @brief Mutates the buggers.
 */
static int ga_pop_mutate( ga_population_t *pop, const ga_options_t *opts )
{
   int i, p;
   double r;

   p = (int)round(((double)pop->pop) * opts->eliteness);
   for (i=p; i<pop->pop; i++) {
      r = rand_double();
      if (r < opts->mutation)
         ga_ent_mutate( &pop->entities[i], pop->parent, opts );
   }
   return 0;
}


static int ga_pop_converge( ga_population_t *pop, const ga_options_t *opts )
{
   int ret;

   /* MINPACK triggers them, so we try to ignore minpack exceptions, although
    * this means errors could sneak by through minpack. */
   if (opts->sigfpe)
      fedisableexcept( FPE_FLAGS );
   ret = ga_pop_pthread( pop, opts, ga_ent_converge_pthread );
   if (opts->sigfpe)
      feenableexcept( FPE_FLAGS );

   return ret;
}


/**
 * @brief Evaluates a population.
 */
static int ga_pop_evaluate( ga_population_t *pop, const ga_options_t *opts, ga_info_t *info )
{
   int i;
   double mean, stddev, d, best;
   ga_entity_t *ent;
   (void) opts;

   /* Calculate all the fitness. */
   ga_pop_pthread( pop, opts, ga_ent_evaluate_pthread );

   /* Evaluation and mean. */
   best = 0.;
   mean = 0.;
   pop->fitness_total = 0.;
   for (i=0; i<pop->pop; i++) {
      ent   = &pop->entities[i];

      /* Calculate mean and main. */
      mean += ent->fitness;
      if (ent->fitness > best)
         best = ent->fitness;

      /* Track total fitness. */
      pop->fitness_total += ent->fitness;
   }
   mean /= (double)pop->pop;

   /* Standard deviation. */
   stddev = 0.;
   for (i=0; i<pop->pop; i++) {
      ent      = &pop->entities[i];
      d        = mean - ent->fitness;
      stddev  += d*d;
   }
   stddev /= (double) pop->pop;

   /* Set information. */
   info->fit_best    = best;
   info->fit_mean    = mean;
   info->fit_stddev  = stddev;

   return 0;
}


/**
 * @brief Compares the fitness of two entities.
 */
static int ga_ent_comp( const void *p1, const void *p2 )
{
   ga_entity_t *e1, *e2;

   e1 = (ga_entity_t*) p1;
   e2 = (ga_entity_t*) p2;

   if (e1->fitness < e2->fitness)
      return +1;
   if (e1->fitness > e2->fitness)
      return -1;
   return 0;
}


/**
 * @brief Sorts the population.
 */
static int ga_pop_sort( ga_population_t *pop, const ga_options_t *opts )
{
   (void) opts;
   qsort( pop->entities, (size_t)pop->pop, sizeof(ga_entity_t), ga_ent_comp );
   return 0;
}


/**
 * @brief Genetic algorithm signal handler.
 */
static void ga_signal_handler( int sig, siginfo_t *info, void *unused )
{
   (void) sig;
   (void) info;
   (void) unused;

   assert( ga_finish < 3 );
   ga_finish++;
}


/**
 * @brief SIGFPE signal handler.
 */
static void ga_signal_fpe( int sig, siginfo_t *info, void *unused )
{
   (void) sig;
   (void) info;
   (void) unused;

   fprintf( stderr, "SIGFPE!\n" );
   exit( EXIT_FAILURE );
}


/**
 * @brief see if must stop.
 */
static int must_stop( ga_options_t *opts_use, ga_info_t *ga_info )
{
   /* Check time and epoch. */
   if ((opts_use->f_epoch != NULL) && opts_use->f_epoch( opts_use, ga_info ))
      return 1;

   /* Check the time stop condition. */
   if ((opts_use->stop_elapsed != 0) && (opts_use->stop_elapsed < ga_info->elapsed))
      return 1;

   /* Check the fitness stop condition. */
   if ((opts_use->stop_fitness != 0.) && (opts_use->stop_fitness < ga_info->fit_best))
      return 1;

   /* Finish condition. */
   if (ga_finish)
      return 1;

   return 0;
}

/**
 * @brief Solves the kinematic synthesis using genetic algorithms.
 *
 * The process is:
 *  - Seed
 *  - Local convergence
 *  - Evaluate
 *  - Loop for generations
 *    - Eliteness passes directly.
 *    - Crossover is run.
 *    - Mutation affects all entities.
 *    - Local convergence.
 *    - Evaluate selection.
 *
 *    @param syn Base synthesis to work off of.
 *    @param opts Options to use.
 *    @param[out] info Results of the solving.
 *    @return 0 on success.
 */
int syn_solve_ga( synthesis_t *syn, ga_options_t *opts, ga_info_t *info )
{
   int i, pc, po, v;
   ga_population_t pop[2];
   ga_options_t def_opts;
   ga_options_t *opts_use;
   ga_info_t ga_info;
   struct timeval tstart, tend;
   struct sigaction sa, so_int, so_usr1, so_fpe;

   /* Must be finalized. */
   assert( syn->finalized != 0 );

   /* Handle options. */
   ga_options_default( &def_opts );
   opts_use = (opts != NULL) ? opts : &def_opts;
   assert( opts_use->population > 0 );
   assert( opts_use->generations > 0 );
   assert( (opts_use->crossover >= 0.) && (opts_use->crossover <= 1.) );
   assert( (opts_use->mutation >= 0.) && (opts_use->mutation <= 1.) );
   assert( (opts_use->eliteness >= 0.) && (opts_use->eliteness <= 1.) );
   assert( (opts_use->mutation + opts_use->eliteness) <= 1. );
   assert( opts_use->seed_mul >= 1. );

   /* Verbosity. */
   v = opts_use->verbose;

   /* Set up signal handler if possible. */
   sa.sa_handler   = NULL;
   sa.sa_sigaction = ga_signal_handler;
   sigemptyset( &sa.sa_mask );
   sa.sa_flags     = SA_SIGINFO;
   if (opts_use->stop_sigint) {
      LOG( v>0, "Enabling SIGINT catcher.\n" );
      sigaction(SIGINT, &sa, &so_int);
   }
   if (opts_use->stop_sigusr1) {
      LOG( v>0, "Enabling SIGUSR1 catcher.\n" );
      sigaction(SIGUSR1, &sa, &so_usr1);
   }

   /* FPU Exceptions. */
   if (opts_use->sigfpe) {
      feenableexcept( FPE_FLAGS );
     // | FE_OVERFLOW );
      sa.sa_sigaction = ga_signal_fpe;
      sigaction(SIGFPE, &sa, &so_fpe );
      LOG( v>0, "FPE enabled.\n" );
   }

   /* Misc logging. */
   if (opts_use->stop_fitness != 0.)
      LOG( v>0, "Stopping at fitness %.3e.\n", opts_use->stop_fitness );
   if (opts_use->stop_elapsed != 0)
      LOG( v>0, "Stopping at %02ld:%02ld:%02ld elapsed time.\n",
            opts_use->stop_elapsed/3600,
            (opts_use->stop_elapsed/60)%60,
            opts_use->stop_elapsed%60 );
   LOG( v>0, "Using population of %d for %d generations.\n",
         opts_use->population, opts_use->generations );
   if (opts_use->converge)
      LOG( v>0, "Using MINPACK for convergance.\n" );

   /* Initialize randomness. */
   rand_init();

   /* Initialize output. */
   memset( &ga_info, 0, sizeof(ga_info_t) );
   ga_info.max_generations = opts_use->generations;

   /* Initialize generation. */
   pc = 1;
   memset( pop, 0, sizeof(pop) );
   for (i=0; i<2; i++) {
      pop[i].pop_target = opts_use->population;
      pop[i].parent     = syn;
      ga_pop_pthread_init( &pop[i], opts_use );
   }

   /* Start timer. */
   gettimeofday( &tstart, NULL );

   /* Allocate the population. */
   ga_pop_seed(      &pop[pc], opts_use );
   if (opts_use->converge)
      ga_pop_converge(  &pop[pc], opts_use );
   ga_pop_evaluate(  &pop[pc], opts_use, &ga_info );
   ga_pop_sort(      &pop[pc], opts_use );

   /* Seed information. */
   if (v > 0) {
      double *vec, mean, stddev, min, max;
      int j, n;
      n = pop[pc].pop;
      vec = malloc( n*sizeof(double) );
      for (j=0; j<n; j++)
         vec[j] = pop[pc].entities[j].fitness;
      rand_stats( &mean, &stddev, &min, &max, vec, n );
      LOG( v>0, "   seed best = %.3e, mean = %.3e, stddev = %.3e, range = [%.3e, %.3e]\n",
            pop[pc].entities[0].fitness, mean, stddev, min, max );
      free( vec );
   }

   /* Update info. */
   gettimeofday( &tend, NULL );
   ga_info.elapsed = (unsigned long)((tend.tv_sec - tstart.tv_sec) + (tend.tv_usec - tstart.tv_usec)/1000000);
   if (must_stop( opts_use, &ga_info ))
      goto ga_done;

   /* Core algorithm. */
   for (i=0; i<opts_use->generations; i++) {
      /* Swap old and current populations. */
      po = pc;
      pc = 1-pc;

      /* Free current population. */
      ga_pop_cleanup(   &pop[pc] );
      LOG( v>2, "Generation %d start\n", i+1 );

      /* Eliteness. */
      ga_pop_eliteness( &pop[pc], &pop[po], opts_use );
      LOG( v>2, "   eliteness: %d\n", pop[pc].pop );

      /* Crossover. */
      ga_pop_crossover( &pop[pc], &pop[po], opts_use );
      LOG( v>2, "   crossover: %d\n", pop[pc].pop );

      /* Mutation. */
      ga_pop_mutate(    &pop[pc], opts_use );
      LOG( v>2, "   mutated\n" );

      /* Convergence. */
      if (opts_use->converge) {
         ga_pop_converge( &pop[pc], opts_use );
         LOG( v>2, "   converged\n" );
      }

      /* Evaluate population. */
      ga_pop_evaluate(  &pop[pc], opts_use, &ga_info );
      LOG( v>2, "   evaluated\n" );

      /* Sort population. */
      ga_pop_sort(      &pop[pc], opts_use );

      /* If verbosity is enabled display generation results. */
      if (v > 0) {
         double *vec, mean, stddev, min, max;
         int j, n;
         n = pop[pc].pop;
         vec = malloc( n*sizeof(double) );
         for (j=0; j<n; j++)
            vec[j] = pop[pc].entities[j].fitness;
         rand_stats( &mean, &stddev, &min, &max, vec, n );
         LOG( v==1, "   [%d] best = %.3e, mean = %.3e, stddev = %.3e, range = [%.3e, %.3e]\n",
               i+1, pop[pc].entities[0].fitness, mean, stddev, min, max );
         LOG( v>2, "   best = %.3e, mean = %.3e, stddev = %.3e, range = [%.3e, %.3e]\n",
               pop[pc].entities[0].fitness, mean, stddev, min, max );
         free( vec );
      }

      /* Update info. */
      gettimeofday( &tend, NULL );
      ga_info.elapsed = (unsigned long)((tend.tv_sec - tstart.tv_sec) + (tend.tv_usec - tstart.tv_usec)/1000000L);
      ga_info.generations = i+1;

      if (must_stop( opts_use, &ga_info ))
         goto ga_done;
   }

ga_done:
   /* Restore signal handlers. */
   if (opts_use->stop_sigint)
      sigaction(SIGINT, &so_int, NULL);
   if (opts_use->stop_sigusr1)
      sigaction(SIGUSR1, &so_usr1, NULL);
   if (opts_use->sigfpe)
      sigaction(SIGFPE, &so_fpe, NULL);

   /* Copy info. */
   if (info != NULL)
      memcpy( info, &ga_info, sizeof(ga_info_t) );

   /* Choose best. */
   syn_free( syn );
   syn_copy( syn, pop[pc].entities[0].syn );
   syn_calc( syn ); /* Update it with calculations. */

   /* Clean up. */
   for (i=0; i<2; i++) {
      ga_pop_cleanup( &pop[i] );
      ga_pop_pthread_free( &pop[i], opts_use );
   }
   rand_exit();

   return 0;
}



