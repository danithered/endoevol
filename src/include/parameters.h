#ifndef _MYPARAMS_
#define _MYPARAMS_
 
#include <iostream>
#include <cstring>
#include <fstream>
 
extern double par_length_dependence;
extern double par_death;
extern double par_substitution;
extern double par_insertion;
extern double par_deletion;
extern double par_k_noasso;
extern double par_decay_rate;
extern double par_diffusion_rate;

extern int par_maxtime;
extern int par_ncol;
extern int par_nrow;
extern int par_movie_interval;
extern int par_output_interval;
 
extern char par_pic_folder[255];
extern char par_ID[255];


//Functions 
int paramsToFile(char* filename);
int Args(int argc, char **argv);

#endif
