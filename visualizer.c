

#include "visualizer.h"

#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <unistd.h>

#include <png.h>

#include "SDL.h"
#include "SDL_opengl.h"

#include <dq/dq.h>
#include <dq/mat3.h>
#include <dq/vec3.h>
#include <dq/homo.h>

#include "sminpack.h"

#define WARN(str, args...) \
(fprintf(stderr,"Warning: [%s] "str, __func__, ## args))

#define MAX(a,b)  (((a)>(b))?(a):(b))
#define MIN(a,b)  (((a)<(b))?(a):(b))

/* Recommended for compatibility and such */
#if HAS_BIGENDIAN
#  define RMASK   0xff000000 /**< Red bit mask. */
#  define GMASK   0x00ff0000 /**< Green bit mask. */
#  define BMASK   0x0000ff00 /**< Blue bit mask. */
#  define AMASK   0x000000ff /**< Alpha bit mask. */
#else
#  define RMASK   0x000000ff /**< Red bit mask. */
#  define GMASK   0x0000ff00 /**< Green bit mask. */
#  define BMASK   0x00ff0000 /**< Blue bit mask. */
#  define AMASK   0xff000000 /**< Alpha bit mask. */
#endif
#define RGBAMASK  RMASK,GMASK,BMASK,AMASK


#define MAX_ANG_VELOCITY   (2.*M_PI/4.)


/*
 * Global stuff.
 */
static int rotating     = 0;
static int moving       = 0;
static double rho       = 1.;
static double theta     = 0.;
static double phi       = 0.;
static double yawvel, pitchvel, rollvel;
static double vis_zoom   = 1.;
static double vis_center[3] = { 0., 0., 0. };
static double vis_eye[3];
static double vis_at[3];
static double vis_up[3] = { 0., 0., 1. };
static int vis_width     = 800;
static int vis_height    = 600;


static const double motion_k    = 0.5;


/*
 * Prototypes.
 */
static void gl_screenshot( const char *filename );
static int write_png( const char *file_name, png_bytep *rows,
      int w, int h, int colourtype, int bitdepth );
static int vis_updateData( GLfloat *data, GLfloat *col, const synthesis_t *syn,
      GLfloat col_base[3], GLfloat col_axes[3], int frame, double state, dq_t *syn_P,
      const double *angles );


static double ang_clean( double a )
{
   if (fabs(a) > 2.*M_PI)
      a = fmod( a, 2.*M_PI );
   if (a > M_PI)
      a -= 2.*M_PI;
   if (a < -M_PI)
      a += 2.*M_PI;
   return a;
}

static double ang_diff( double a, double b )
{
   double diff = ang_clean( ang_clean(b) - ang_clean(a) );
   diff = (diff < +M_PI) ? diff : diff - 2.*M_PI;
   diff = (diff > -M_PI) ? diff : diff + 2.*M_PI;
   return diff;
}


/**
 * @brief Propagates a point to a new axes.
 */
static void vis_propagatePosition( double next[3], double prop[3],
      const double cur[3], const double n[3],
      const double s[3], const double s0[3] )
{
   int i;
   double x1[3], x2[3], dn[3], A[3][3], b[3];
   const double *d1, *d2, *c2;
   double c1[3];
   double px1a[3], px1b[3], px2a[3], px2b[3];
   double u;

   /*
    * Subindex numbering
    * 1 =>  to ( c,   s )
    * 2 => from( cur, n )
    *
    * n = d1 x d2
    * b = c2 - c1
    * A = d1 -d2 -n
    */
   d1 = s;
   vec3_cross( c1, s, s0 );
   d2 = n;
   c2 = cur;
   vec3_cross( dn, d1, d2 );
   vec3_sub(   b, c2, c1 );
   for (i=0; i<3; i++) {
      A[i][0] =  d1[i];
      A[i][1] = -d2[i];
      A[i][2] = -dn[i];
   }

   /* In this case we must handle it automagically. */
   if (fabs(mat3_det(A)) < DQ_PRECISION) {

      /* This represents the magnitude of error on the direction of the s vector. */
      u = vec3_dot( b, s );

      /* Calculate intersection point on line from current point. */
      for (i=0; i<3; i++) {

         /* We displace the point along the line using the magnitude of error. */
         next[i]  = c1[i] + u * s[i];
         prop[i]  = next[i];
      }
      return;
   }

   /* We now have s, t and k. */
   mat3_solve( x1, A, b );
   for (i=0; i<3; i++)
      A[i][2] = -A[i][2];
   mat3_solve( x2, A, b );

   /* Calculate both possibilities. */
   for (i=0; i<3; i++) {
      px1a[i] = c1[i] + x1[0] * d1[i];
      px1b[i] = c2[i] + x1[1] * d2[i];

      px2a[i] = c1[i] + x2[0] * d1[i];
      px2b[i] = c2[i] + x2[1] * d2[i];
   }

   /* Get minimum distance. */
   if (vec3_distance( px1a, px1b ) < vec3_distance( px2a, px2b )) {
      memcpy( next, px1b, sizeof(double)*3 );
      memcpy( prop, px1a, sizeof(double)*3 );
   }
   else {
      memcpy( next, px2b, sizeof(double)*3 );
      memcpy( prop, px2a, sizeof(double)*3 );
   }
}


/**
 * @brief Copies doubles to floats.
 */
static void fcopy( GLfloat out[3], const double in[3] )
{
   int i;
   for (i=0; i<3; i++)
      out[i] = (GLfloat) in[i];
}


/**
 * @brief Allocates the data.
 */
static GLfloat* vis_allocData( int *n, const synthesis_t *syn, GLfloat **col )
{
   int i, dof;
   GLfloat *data;

   /* Load memory. */
   dof = 0;
   for (i=0; i<syn->nbranches; i++)
      dof += syn->branches[i].njoints + 2;
   *n = dof*3*3*2;
   data = calloc( (*n), sizeof(GLfloat) );

   /* Allocate colour. */
   *col = calloc( (*n), sizeof(GLfloat) );

   return data;
}


/**
 * @brief Frees the allocated data.
 */
static void vis_freeData( GLfloat *data, GLfloat *col )
{
   free( data );
   free( col );
}


typedef struct vis_syn_s {
   const synthesis_t *syn;
   dq_t *buf;
   dq_t *target;
} vis_syn_t;
static int vis_eqns( void *p, int m, int n, const double *x, double *fvec, int iflag )
{
   (void) iflag;
   (void) m;
   (void) n;
   int i;
   dq_t T;
   vis_syn_t *vsyn = (vis_syn_t*) p;

   /* Forward kinematics. */
   syn_fk( vsyn->syn, vsyn->buf, x );

   /* Copy results over. */
   for (i=0; i<vsyn->syn->ntcp; i++) {
      dq_op_sub( T, vsyn->buf[i], vsyn->target[i] );
      memcpy( &fvec[8*i], T, sizeof(double)*8 );
   }
   return 0;
}
static int vis_updateDataFrom( const synthesis_t *syn,
      double *angles, dq_t *P, int frame, double state )
{
   vis_syn_t vsyn;
   double *fvec;
   int i, nfev;

   if (fabs(state) < 1e-10) {
      for (i=0; i<syn->njoints; i++) {
         if (frame == 0)
            angles[i] = 0.;
         else
            angles[i] = syn->joints[i]->pos[ frame-1 ];
      }
   }
   else {
      /* Set up. */
      vsyn.syn    = syn;
      vsyn.buf    = calloc( syn->ntcp, sizeof(dq_t) );
      vsyn.target = P;
      fvec        = malloc( syn->ntcp*8*sizeof(double) );

      /* Minpack. */
      sminpack( vis_eqns, (void*)&vsyn, syn->ntcp*8, syn->njoints, angles, fvec,
            0., 0., 0., 10000, 0., 1, 100., 0, &nfev, NULL, NULL );

      /* Clean up. */
      free( vsyn.buf );
      free( fvec );
   }

   /*
   for (i=0; i<syn->njoints; i++)
      ang_clean( angles[i] );
   */

   return 0;
}
static int vis_updateDataInterpolate( GLfloat *data, const synthesis_t *syn,
      double *angles_cur, double *angles, double dt )
{
   int i;

   /* Update position. */
   for (i=0; i<syn->njoints; i++) {
      double off = ang_diff( angles_cur[i], angles[i] );
      if (fabs(off) < motion_k*MAX_ANG_VELOCITY*dt)
         angles_cur[i]  = angles[i];
      else
         angles_cur[i] += (off/fabs(off)) * motion_k*MAX_ANG_VELOCITY*dt;
   }

   /* Update data. */
   vis_updateData( data, NULL, syn, NULL, NULL, 0, 0., NULL, angles_cur );

   return 0;
}


/**
 * @brief Loads visual data.
 */
static int vis_updateData( GLfloat *data, GLfloat *col, const synthesis_t *syn,
      GLfloat col_base[3], GLfloat col_axes[3], int frame, double state, dq_t *syn_P,
      const double *angles )
{
   int i, j, k, p, c, f;
   dq_t T, R, L, P;
   double s[3], s0[3], o[3], a[3], z[3], ang;
   double M[3][3];
   double next[3], prop[3];
   const double pz[3] = { 0., 0., 0. };
   kin_branch_t *br;
   kin_joint_t *kj;

   /* Copy angles. */
   if (angles != NULL)
      for (i=0; i<syn->njoints; i++)
         syn->joints[i]->pos_cur = angles[i];

   /* Fill memory.
    * We're going to use GL_LINES, so we have to store start and end point of every line. */
   p = 0;
   c = 0;
   f = frame;
   for (i=0; i<syn->nbranches; i++) {
      br = &syn->branches[i];
      /* Start at origin. */
      o[0] = o[1] = o[2] = 0.;
      z[0] = z[1] = 0.;
      z[2] = 1.;

      /* Initialize transformation. */
      dq_cr_point( T, pz );

      /* Deal with all the branch joints. */
      for (j=0; j<br->njoints; j++) {
         kj = br->joints[j];

         /* Transform the line using previous transformation. */
         dq_cr_line_plucker( L, kj->S.s, kj->S.s0 );
         dq_op_f2g( L, T, L );

         /* Extract updated screw axis from dual quaternion. */
         s[0]  = L[1];
         s[1]  = L[2];
         s[2]  = L[3];
         s0[0] = L[4];
         s0[1] = L[5];
         s0[2] = L[6];

         /* Update the transform with current line untransformed. */
         if (angles == NULL) {
            if (f==0)
               ang = state * kj->pos[0];
            else if (f<syn->L-1) {
               double diff = ang_diff( kj->pos[f-1], kj->pos[f] );
               ang = kj->pos[f-1] + state * diff;
            }
            else
               ang = kj->pos[ f-1 ];
         }
         else
            ang = kj->pos_cur;
         dq_cr_rotation_plucker( R, ang, kj->S.s, kj->S.s0 );
         dq_op_mul( T, T, R );

         /* Copy start. */
         fcopy( &data[3*(p++)], o );

         /* Calculate displacement from last. */
         vis_propagatePosition( next, prop, o, z, s, s0 );
         memcpy( z, s, sizeof(double)*3 );

         /* Copy end. */
         fcopy( &data[3*(p++)], next );
         fcopy( &data[3*(p++)], next );
         fcopy( &data[3*(p++)], prop );

         /* Create axis. */
         memcpy( o, prop, sizeof(double)*3 );
         memcpy( a, o,    sizeof(double)*3 );
         a[0] += s[0]/2.;
         a[1] += s[1]/2.;
         a[2] += s[2]/2.;
         fcopy( &data[3*(p++)], a );
         vec3_sub( a, a, s );
         fcopy( &data[3*(p++)], a );

         /* Colour. */
         if (col != NULL) {
            memcpy( &col[3*(c++)], col_base, sizeof(GLfloat)*3 );
            memcpy( &col[3*(c++)], col_base, sizeof(GLfloat)*3 );
            memcpy( &col[3*(c++)], col_base, sizeof(GLfloat)*3 );
            memcpy( &col[3*(c++)], col_base, sizeof(GLfloat)*3 );
            memcpy( &col[3*(c++)], col_axes, sizeof(GLfloat)*3 );
            memcpy( &col[3*(c++)], col_axes, sizeof(GLfloat)*3 );
         }
      }

      /* Now we need the FK. */
      fcopy( &data[3*(p++)], o );
      dq_op_mul( P, T, br->tcp->d.tcp.P[0] );
      dq_op_extract( M, o, P );
      fcopy( &data[3*(p++)], o );

      /* Copy over. */
      if (syn_P != NULL)
         dq_cr_copy( syn_P[i], P );

      /* Null point. */
      for (k=0; k<3; k++)
         a[k] = o[k] + M[k][0]/3.;
      fcopy( &data[3*(p++)], o );
      fcopy( &data[3*(p++)], a );
      for (k=0; k<3; k++)
         a[k] = o[k] + M[k][1]/3.;
      fcopy( &data[3*(p++)], o );
      fcopy( &data[3*(p++)], a );
      for (k=0; k<3; k++)
         a[k] = o[k] + M[k][2]/3.;
      fcopy( &data[3*(p++)], o );
      fcopy( &data[3*(p++)], a );

      /* Colour. */
      if (col != NULL) {
         memcpy( &col[3*(c++)], col_base, sizeof(GLfloat)*3 );
         memcpy( &col[3*(c++)], col_base, sizeof(GLfloat)*3 );
         memcpy( &col[3*(c++)], col_axes, sizeof(GLfloat)*3 );
         memcpy( &col[3*(c++)], col_axes, sizeof(GLfloat)*3 );
         memcpy( &col[3*(c++)], col_axes, sizeof(GLfloat)*3 );
         memcpy( &col[3*(c++)], col_axes, sizeof(GLfloat)*3 );
         memcpy( &col[3*(c++)], col_axes, sizeof(GLfloat)*3 );
         memcpy( &col[3*(c++)], col_axes, sizeof(GLfloat)*3 );
      }
   }

   return 0;
}

static void vis_camRotate( double yaw, double pitch, double roll )
{
   (void) roll;
   theta += yaw;
   phi += pitch;
   if (fabs(phi) > M_PI/2.)
      phi = M_PI/2. * (phi / fabs(phi));
}

static void vis_camAbsolute( double yaw, double pitch, double roll )
{
   (void) roll;
   theta = yaw;
   phi = pitch;
   if (fabs(phi) > M_PI/2.)
      phi = fmod( phi, M_PI/2. );
}


static void vis_camMove( double x, double y, double z )
{
   vis_center[0] += x;
   vis_center[1] += y;
   vis_center[2] += z;
}


static void vis_camUpdate( double dt )
{
   vis_camRotate( dt*yawvel, dt*pitchvel, dt*rollvel );
}


static void vis_setCamera (void)
{
   int i;
   double r, z;
   double up[3] = { 0., 0., 1. };

   r        = rho * cos( phi );
   z        = rho * sin( phi );
   vis_eye[0]   = r * cos( theta );
   vis_eye[1]   = r * sin( theta );
   vis_eye[2]   = z;
   for (i=0; i<3; i++)
      vis_at[i] = -vis_eye[i];

   /* Set up the matrix. */
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   glOrtho( vis_width/2, -vis_width/2, vis_height/2, -vis_height/2, -10000, 10000 );
   gluLookAt( 0., 0., 0.,
         vis_eye[0], vis_eye[1], vis_eye[2],
         up[0], up[1], up[2] );
   glScalef(  vis_zoom, vis_zoom, vis_zoom );
   glTranslatef( vis_center[0], vis_center[1], vis_center[2] );
}


/**
 * @brief Gets the largest dimensional size of the data.
 */
static double vis_dataSize( double center[3], const GLfloat *data, int n )
{
   int i, j;
   double max, f;

   if (data == NULL)
      return 0.;

   max = 0.;
   for (i=0; i<n/3; i++) {
      for (j=0; j<3; j++) {
         f = fabs( center[j] - data[i*3+j] );
         if (f > max)
            max = f;
      }
   }
   return (double) max;
}


/**
 * @brief Gets the center of the data.
 */
#if 0
static void vis_dataCenter( double center[3], const GLfloat *data, int n )
{
   int i, j;
   double acc[3];

   acc[0] = acc[1] = acc[2] = 0.;
   for (i=0; i<n/3; i++)
      for (j=0; j<3; j++)
         acc[j] += data[i*3+j];

   for (j=0; j<3; j++)
      center[j] = acc[j] / (double)(n/3);
}
#endif


/**
 * @brief Renders the visual data.
 */
static void vis_renderData( const synthesis_t *syn,
      const GLfloat *data, const GLfloat *col,
      int n, int finger )
{
   (void)syn;
   int i, fing_offset, fing_size;

   glMatrixMode(GL_MODELVIEW);
   glLoadIdentity();
   glRotatef( 180., 0., 1., 0. );

   /* Render. */
   glEnableClientState( GL_VERTEX_ARRAY );
   glEnableClientState( GL_COLOR_ARRAY );

   if (finger < 0) {
      fing_offset = 0;
      fing_size = n / 3;
   }
   else {
      fing_offset = 0;
      for (i=0; i<finger; i++)
         fing_offset += syn->branches[i].njoints + 2;
      fing_offset *= 3*3*2;
      fing_size = (syn->branches[finger].njoints+2)*3*2;
   }
   glVertexPointer( 3, GL_FLOAT, 0, &data[fing_offset] );
   glColorPointer(  3, GL_FLOAT, 0, &col[fing_offset] );

   glDrawArrays( GL_LINES, 0, fing_size );

   glDisableClientState( GL_COLOR_ARRAY );
   glDisableClientState( GL_VERTEX_ARRAY );
}


static int gl_setup (void)
{
   GLfloat col_bg[3]    = { 1.0, 1.0, 1.0 };

   /* Fancy opengl settings. */
   glViewport( 0, 0, vis_width, vis_height );
   glClearColor( col_bg[0], col_bg[1], col_bg[2], 1.0 );
   glEnable(  GL_DEPTH_TEST );
   glDisable( GL_TEXTURE_2D );
   glLineWidth( 3. );
   glHint(    GL_LINE_SMOOTH_HINT, GL_NICEST );
   glEnable(  GL_LINE_SMOOTH );
   return 0;
}


/**
 * @brief Visualizes a set of equations.
 */
int visualize( const synthesis_t *syn_a_in, const synthesis_t *syn_b_in )
{
   int quit, max_frame, interpolate;
   int ndata, ncomp, cur_frame, cur_finger, target_frame;
   double dx, dy, max;
   double state, dt;
   GLfloat *vdata, *vcomp, *vdatacol, *vcompcol;
   SDL_Event evt;
   synthesis_t *syn, *syn_b;
   dq_t *syn_P;
   double *angles_cmp, *angles_cur;

   /* Duplicate so we can mess with it. */
   syn      = NULL;
   syn_b    = NULL;
   if (syn_a_in != NULL)
      syn      = syn_dup( syn_a_in );
   if (syn_b_in != NULL)
      syn_b    = syn_dup( syn_b_in );

   /* Colours. */
#if 0
   GLfloat col_bg[3]    = { 0.8, 0.8, 0.8 };
   GLfloat col[3]       = { 1.0, 1.0, 1.0 };
   GLfloat cola[3]      = { 0.5, 0.5, 0.5 };
   GLfloat col_comp[3]  = { 0.2, 0.8, 0.2 };
   GLfloat col_compa[3] = { 0.2, 0.2, 0.8 };
#else
   GLfloat col[3]       = { 0.8, 0.8, 0.8 };
   GLfloat cola[3]      = { 0.3, 0.3, 0.3 };
   GLfloat col_comp[3]  = { 0.2, 0.8, 0.2 };
   GLfloat col_compa[3] = { 0.2, 0.2, 0.8 };
#endif

   /* Must be finalized. */
   assert( syn->finalized );

   /* Allocate. */
   syn_P       = calloc( syn->ntcp, sizeof(dq_t) );
   angles_cmp  = calloc( syn->njoints, sizeof(double) );
   angles_cur  = calloc( syn->njoints, sizeof(double) );

   /* Construct fingers. */
   vdata = vis_allocData( &ndata, syn, &vdatacol );
   vis_updateData( vdata, vdatacol, syn, col, cola, 0, 0., NULL, NULL );
   vcomp = NULL;
   if (syn_b != NULL) {
      vcomp = vis_allocData( &ncomp, syn_b, &vcompcol );
      vis_updateData( vcomp, vcompcol, syn_b, col_comp, col_compa, 0, 0., NULL, NULL );
   }

   /* Create window. */
   SDL_Init( SDL_INIT_VIDEO );
   SDL_GL_SetAttribute( SDL_GL_DOUBLEBUFFER, 1 );
   SDL_SetVideoMode( vis_width, vis_height, 0, SDL_OPENGL | SDL_RESIZABLE );
   gl_setup();

   /* Center stuff. */
   max = vis_dataSize( vis_center, vdata, ndata );
   if (syn_b != NULL)
      max   = MAX( max, vis_dataSize( vis_center, vcomp, ncomp ) );
   vis_zoom  = 0.75 * MIN(vis_height, vis_width) / max;

   /* Camera. */
   vis_camAbsolute( M_PI/4., M_PI/4., 0. );

   /* Render. */
   quit        = 0;
   cur_frame   = 0;
   target_frame = 0;
   cur_finger  = -1;
   interpolate = 1;
   max_frame   = syn->L;
   if (syn_b != NULL)
      max_frame   = MIN( max_frame, syn_b->L );
   while (!quit) {

      /* Events. */
      while (SDL_PollEvent( &evt )) {
         if (evt.type == SDL_QUIT)
            quit = 1;
         else if (evt.type == SDL_VIDEORESIZE) {
            vis_width   = evt.resize.w;
            vis_height  = evt.resize.h;
            SDL_SetVideoMode( vis_width, vis_height, 0, SDL_OPENGL | SDL_RESIZABLE );
            gl_setup();
         }
         else if (evt.type == SDL_KEYDOWN) {
            if (evt.key.keysym.sym == SDLK_ESCAPE)
               quit = 1;
            else if (evt.key.keysym.sym == SDLK_f) {
               cur_finger += 1;
               if (cur_finger > syn->nbranches-1)
                  cur_finger = -1;
            }
            else if (evt.key.keysym.sym == SDLK_i)
               interpolate = !interpolate;
            else if (evt.key.keysym.sym == SDLK_DOWN)
               pitchvel = M_PI/2.;
            else if (evt.key.keysym.sym == SDLK_UP)
               pitchvel = -M_PI/2.;
            else if (evt.key.keysym.sym == SDLK_LEFT)
               yawvel   = -M_PI/2.;
            else if (evt.key.keysym.sym == SDLK_RIGHT)
               yawvel   = M_PI/2.;
            else if (evt.key.keysym.sym == SDLK_n) {
               target_frame  += 1;
               cur_frame      = target_frame-1;
               state          = 0.0;
               if (target_frame >= max_frame)
                  target_frame = max_frame-1;
               if (!interpolate || (cur_frame==target_frame)) {
                  cur_frame = target_frame;
                  vis_updateData( vdata, NULL, syn, col, cola, cur_frame, state, syn_P, NULL );
                  if (syn_b != NULL)
                     vis_updateDataFrom( syn_b, angles_cmp, syn_P, cur_frame, state );
               }
            }
            else if (evt.key.keysym.sym == SDLK_p) {
               target_frame  -= 1;
               cur_frame      = target_frame+1;
               if (target_frame < 0) {
                  target_frame   = 0;
                  state          = 0.0;
               }
               else
                  state          = 1.0;
               if (!interpolate || (cur_frame==target_frame)) {
                  cur_frame = target_frame;
                  vis_updateData( vdata, NULL, syn, col, cola, cur_frame, state, syn_P, NULL );
                  if (syn_b != NULL)
                     vis_updateDataFrom( syn_b, angles_cmp, syn_P, cur_frame, state );
               }
            }
            else if (evt.key.keysym.sym == SDLK_s)
               gl_screenshot( "out.png" );
         }
         else if (evt.type == SDL_KEYUP) {
            if      (evt.key.keysym.sym == SDLK_DOWN) {
               pitchvel = 0.;
            }
            else if (evt.key.keysym.sym == SDLK_UP) {
               pitchvel = 0.;
            }
            else if (evt.key.keysym.sym == SDLK_LEFT) {
               yawvel   = 0.;
            }
            else if (evt.key.keysym.sym == SDLK_RIGHT) {
               yawvel   = 0.;
            }
         }
         else if (evt.type == SDL_MOUSEBUTTONDOWN) {
            if           (evt.button.button == SDL_BUTTON_LEFT)
               rotating = 1;
            else if      (evt.button.button == SDL_BUTTON_RIGHT)
               moving   = 1;
            else if      (evt.button.button == SDL_BUTTON_WHEELUP)
               vis_zoom *= 1.1;
            else if (evt.button.button == SDL_BUTTON_WHEELDOWN)
               vis_zoom *= 1./1.1;
         }
         else if (evt.type == SDL_MOUSEBUTTONUP) {
            if           (evt.button.button == SDL_BUTTON_LEFT)
               rotating = 0;
            else if      (evt.button.button == SDL_BUTTON_RIGHT)
               moving   = 0;
         }
         else if (evt.type == SDL_MOUSEMOTION) {
            if (rotating) {
               dx = (double) (evt.motion.xrel);
               dy = (double) (evt.motion.yrel);
               vis_camRotate( 0.005*dx, 0.005*dy, 0. );
            }
            if (moving) {
               dx = (double) (evt.motion.xrel);
               dy = (double) (evt.motion.yrel);
               double bx[3], by[3];
               vec3_cross( bx, vis_up, vis_at );
               vec3_cross( by, bx,     vis_at );
               vec3_normalize( bx );
               vec3_normalize( by );
               dx /= -vis_zoom;
               dy /= -vis_zoom;
               vis_camMove( bx[0]*dx + by[0]*dy, bx[1]*dx + by[1]*dy, by[2]*dy );
            }
         }
      }

      dt = 0.01;

      if (cur_frame != target_frame) {
         int f;

         /* Going "up" */
         if (target_frame > cur_frame) {
            if (state < 1.0)
               state += motion_k*dt;
            if (state > 1.0) {
               cur_frame = target_frame;
               state     = 0.0;
            }
            f = cur_frame;
         }
         /* Going down. */
         else {
            f = cur_frame-1;
            if (state > 0.0)
               state -= motion_k*dt;
            if (state < 0.0) {
               cur_frame = target_frame;
               state     = 0.0;
            }
         }

         vis_updateData( vdata, NULL, syn, col, cola, f, state, syn_P, NULL );
         if (syn_b != NULL)
            vis_updateDataFrom( syn_b, angles_cmp, syn_P, f, state );
      }

      /* Always updated. */
      if (syn_b != NULL)
         vis_updateDataInterpolate( vcomp, syn_b, angles_cur, angles_cmp, dt );

      /* Update camera. */
      vis_camUpdate( dt );

      /* Clear buffers. */
      glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

      /* Set up camera. */
      vis_setCamera();

      /* Render. */
      if (vcomp != NULL)
         vis_renderData( syn_b, vcomp, vcompcol, ncomp, cur_finger );
      vis_renderData( syn, vdata, vdatacol, ndata, cur_finger  );

      /* Swap buffers. */
      SDL_GL_SwapBuffers();

      /* Check for OpenGL errors. */
      if (glGetError() != GL_NO_ERROR)
         WARN("OpenGL Error!\n" );

      /* CPU limit. */
      SDL_Delay( 10 );
   }

   /* Clean up. */
   vis_freeData( vdata, vdatacol );
   if (vcomp != NULL)
      vis_freeData( vcomp, vcompcol );
   SDL_Quit();

   /* Clean up synthesis objects. */
   if (syn != NULL)
      syn_free( syn );
   if (syn_b != NULL)
      syn_free( syn_b );

   /* Done. */
   free( syn_P );
   free( angles_cmp );

   return 0;
}


/**
 * @brief Takes a screenshot.
 *
 *    @param filename Name of the file to save screenshot as.
 */
static void gl_screenshot( const char *filename )
{
   GLubyte *screenbuf;
   png_bytep *rows;
   int i, w, h;

   /* Allocate data. */
   w           = vis_width;
   h           = vis_height;
   screenbuf   = malloc( sizeof(GLubyte) * 3 * w*h );
   rows        = malloc( sizeof(png_bytep) * h );

   /* Read pixels from buffer -- SLOW. */
   glPixelStorei(GL_PACK_ALIGNMENT, 1); /* Force them to pack the bytes. */
   glReadPixels( 0, 0, w, h, GL_RGB, GL_UNSIGNED_BYTE, screenbuf );

   /* Convert data. */
   for (i = 0; i < h; i++)
      rows[i] = &screenbuf[ (h - i - 1) * (3*w) ];

   /* Save PNG. */
   write_png( filename, rows, w, h, PNG_COLOR_TYPE_RGB, 8);

   /* Free memory. */
   free( screenbuf );
   free( rows );
}


/**
 * @brief Saves a png.
 *
 *    @param file_name Name of the file to save the png as.
 *    @param rows Rows containing the data.
 *    @param w Width of the png.
 *    @param h Height of the png.
 *    @param colourtype Colour type of the png.
 *    @param bitdepth Bit depth of the png.
 *    @return 0 on success.
 */
static int write_png( const char *file_name, png_bytep *rows,
      int w, int h, int colourtype, int bitdepth )
{
   png_structp png_ptr;
   png_infop info_ptr;
   FILE *fp;

   /* Open file for writing. */
   if (!(fp = fopen(file_name, "wb"))) {
      WARN("Unable to open '%s' for writing.", file_name);
      return -1;
   }

   /* Create working structs. */
   if (!(png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL))) {
      WARN("Unable to create png write struct.");
      goto ERR_FAIL;
   }
   if (!(info_ptr = png_create_info_struct(png_ptr))) {
      WARN("Unable to create PNG info struct.");
      goto ERR_FAIL;
   }

   /* Set image details. */
   png_init_io(png_ptr, fp);
   png_set_compression_level(png_ptr, Z_DEFAULT_COMPRESSION);
   png_set_IHDR(png_ptr, info_ptr, w, h, bitdepth, colourtype,
         PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT,
         PNG_FILTER_TYPE_DEFAULT);

   /* Write image. */
   png_write_info(png_ptr, info_ptr);
   png_write_image(png_ptr, rows);
   png_write_end(png_ptr, NULL);

   /* Clean up. */
   fclose(fp);

   return 0;

ERR_FAIL:
   fclose(fp);
   return -1;
}


