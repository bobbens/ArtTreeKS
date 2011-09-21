


#ifndef _RAND_H
#  define _RAND_H


/*
 * Init/cleanup.
 */
void rand_init (void);
void rand_exit (void);


/*
 * Get values.
 */
int rand_bool (void);
int rand_int_range( int low, int high );
double rand_double (void);
double rand_double_range( double low, double high );
double rand_double_exponential( double lambda );
double rand_double_gaussian (void);


#endif /* _RAND_H */


