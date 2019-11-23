#include <iostream>
#include "randomgen.h"
#include "stringreps.h"
#include "ca.h"
#include "parameters.h"


using namespace std;

gsl_rng * r;

int main(int argc, char *argv[]) {
	//Argoments
	if ( Args(argc, argv) ) {
		return(-1);
	}
	
	//randomszam generator inicializalasa
	time_t timer;
	r = (gsl_rng *) gsl_rng_alloc (gsl_rng_mt19937);
	gsl_rng_set(r, time(&timer));
	
	//variables
	double diff = 0.0, decay = 0.0;
	int what;
	char command[512], output_file_name[255], pic_folder[255];
	
	//init staff
	strrep::Strrep::set_length_activity(par_length_dependence);
	initColors();
	strrep::Sca aut(par_ncol,par_nrow);
	aut.neighAdd();
	
	//initialising output library tree
	sprintf(command, "mkdir -p OUT");
	system(command);
	
	sprintf(command, "mkdir -p OUT/%s\0", par_ID);
	system(command);
	
	sprintf(command, "mkdir -p OUT/%s/%s\0", par_ID, par_pic_folder);
	system(command);
	
	sprintf(output_file_name, "OUT/%s/output.txt\0", par_ID);
	sprintf(pic_folder, "OUT/%s/%s\0", par_ID, par_pic_folder);
	sprintf(command, "OUT/%s/parameters.txt\0", par_ID);
	
	//print parameters
	paramsToFile(command);
	
	//start simulation
	for(int gen=0; gen < par_maxtime && strrep::Strrep::no_repl; gen++) {
/**/		std::cout << "gen " << gen << " with " << strrep::Strrep::no_repl << " replicators" << std::endl;		
		
		if(gen % par_output_interval == 0) aut.Output(output_file_name, gen);
		if(gen % par_movie_interval == 0) aut.Picture(pic_folder, gen);

		//Update
		for(int iter=0; iter < aut.size; iter++){
//			std::cout << "iter " << iter << std::endl;			
			aut.Update(gsl_rng_uniform_int(r, aut.size));
		}
		
		//Decay
		for(decay += par_decay_rate; decay >= 1; decay--) {
			what = gsl_rng_uniform_int(r, aut.size);
//			std::cout << "Decay of molecule " << what << endl;			
			if(aut.get(what)->role != strrep::empty && gsl_rng_uniform(r) < par_death  ) {
//				std::cout << "tested" << std::endl;				
				aut.get(what)->del();
//				std::cout << "deleted" << std::endl;				
			}
		}
		
		//Diffusion
		for(diff += par_diffusion_rate; diff >= 1; diff--) {
			what = gsl_rng_uniform_int(r, aut.size);
			aut.get(what) > aut.rneigh(what);
		}
		
	}
		

	//randomszam generator lezarasa
	gsl_rng_free(r);

	return 0;
} 
