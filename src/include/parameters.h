#ifndef _MYPARAMS_
#define _MYPARAMS_
 
#include <iostream>
#include <cstring>
#include <fstream>
 
#define MAXLEN 300
 
extern double par_length_dependence;
extern double par_death;
extern double par_substitution;
extern double par_insertion;
extern double par_deletion;
extern double par_dissotiation;
//extern double par_k_noasso;
//extern double par_decay_rate;
extern double par_diffusion_rate;
extern double par_backmut;
extern double par_exponent;
extern double par_alpha;

extern int par_maxtime;
extern int par_ncol;
extern int par_nrow;
extern int par_movie_interval;
extern int par_output_interval;
 
extern char par_pic_folder[255];
extern char par_ID[255];


extern char **par_init_seqs;
extern char **par_init_seqs_compl;
extern double *par_init_props;
extern int par_init_no;


//Functions 
int paramsToFile(char* filename);
int Args(int argc, char **argv);

#endif
