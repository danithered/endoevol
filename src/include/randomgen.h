#ifndef _RANDOMSZAM_
#define _RANDOMSZAM_

extern "C" {
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>
}
	
extern gsl_rng * r;

/*A main-be a kovetkezo sorokat be kell tenni az elejere es a vegere:*/
/*
 * 
//r deklaralasa main elott global env-ben
gsl_rng * r;
 
//randomszam generator inicializalasa main elejen
	time_t timer;
	r = (gsl_rng *) gsl_rng_alloc (gsl_rng_mt19937);
	gsl_rng_set(r, time(&timer));

//randomszam generator lezarasa main vegen
	gsl_rng_free(r);
*/

#endif
