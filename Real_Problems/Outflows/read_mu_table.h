#ifndef read_mu_table_h
#define read_mu_table_h
#if MU_CALC == MU_TABLE
/* Make sure if included elsewhere it is 
 * preceded by #include pluto.h */

/* functions */
void ReadMuTable();

/* global variables*/
/* mu_por is p/rho */
extern double *mu_por, *mu_mu;
extern int mu_ndata;

#endif
#endif

