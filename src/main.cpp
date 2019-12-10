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
	int what, where, gen=0;
	char command[512], output_file_path[255], pic_folder[255];
	
	//init staff
	strrep::Strrep::set_length_activity(par_length_dependence);
	initColors();
	strrep::Sca aut(par_ncol,par_nrow);
	aut.neighAdd();
	if(par_init_no) for(int c = 0; c < aut.size; c++) {
		aut.get(c)->setSeq(par_init_seqs[ dvtools::brokenStickVals(par_init_props, par_init_no, 1, gsl_rng_uniform(r)) ]);
	}
	
	//initialising output library tree
	sprintf(command, "mkdir -p OUT");
	system(command);
	
	sprintf(command, "mkdir -p OUT/%s\0", par_ID);
	system(command);
	
	sprintf(command, "mkdir -p OUT/%s/%s\0", par_ID, par_pic_folder);
	system(command);
	
	sprintf(output_file_path, "OUT/%s/\0", par_ID);
	sprintf(pic_folder, "OUT/%s/%s\0", par_ID, par_pic_folder);
	
	//print parameters
	sprintf(command, "OUT/%s/parameters.txt\0", par_ID);
	paramsToFile(command);
	
	//write rng state
	sprintf(command, "OUT/%s/rngsave.bin\0", par_ID);
	//gsl_rng_fwrite (command, r);
	std::cout << "RNG seed: " << (int) timer << std::endl; 
//	int single_count, empty_count, r_count, e_count, rt_count, et_count; 
	
	//start simulation
	for(gen=0; gen < par_maxtime && strrep::Strrep::no_repl>1; gen++) {
/**/		std::cout << "gen " << gen << " with " << strrep::Strrep::no_repl << " replicators" << std::endl;		
//		for(int z =single_count=empty_count=r_count=e_count=rt_count= et_count=0; z<aut.size; z++) {
//			switch(aut.get(z)->role){
//				case strrep::empty : empty_count++; break;
//				case strrep::single : single_count++; break;
//				case strrep::repl : r_count++; break;
//				case strrep::endo : e_count++; break;
//				case strrep::repl_template : rt_count++; break;
//				case strrep::endo_template : et_count++; break;
//			}
//		}
//		cout << "counts: empty: " << empty_count << " single: " << single_count << " repl: " << r_count << " endo: " << e_count << " repl_template: " << rt_count << " endo_template: " << et_count << endl;

		if(gen % par_output_interval == 0) aut.Output(output_file_path, gen);
		if(gen % par_movie_interval == 0) aut.Picture(pic_folder, gen);
		
/*		
		char rrtcee[]="RRRTRRRTRRRRRRRTCTCTCTCTCTCTTTCTCTCTCTTEEEEEEEEEEE\0";
		char eee[]="EEEEEEEEEEE\0";
		if(gen==7000) {
			for(int xx=0; xx<50; xx++ ){
				for(int yy=0; yy <50; yy++){
					if(gsl_rng_uniform(r) < 0.8) aut.get(xx + par_ncol*yy)->setSeq(rrtcee);
					else aut.get(xx + par_ncol*yy)->setSeq(eee);
				}
			}
		}
*/		

		//Update
		for(int iter=0; iter < aut.size; iter++){
//			std::cout << "iter " << iter << std::endl;
//			aut.Output(output_file_name, gen);
			aut.Update(gsl_rng_uniform_int(r, aut.size));
//			aut.Output(output_file_name, gen);		
			//Decay
			if(par_death) {
				what = gsl_rng_uniform_int(r, aut.size);
//				std::cout << "Decay of molecule " << what << endl;			
				if(aut.get(what)->role != strrep::empty && gsl_rng_uniform(r) < par_death  ) {
//					std::cout << "tested" << std::endl;				
					aut.get(what)->del();
//					std::cout << "deleted " << what << std::endl;				
				}
			}
			
			//Diffusion
			for(diff += par_diffusion_rate; diff >= 1; diff--) {
				what = gsl_rng_uniform_int(r, aut.size);
				where = aut.neigh(what, gsl_rng_uniform_int(r,8) + 1, 0);
				aut.get(what)->diff(aut.get(where));
//				cout << "diffusion of replicator " << what << " to " << where << endl; 
			}
			
		}
		
		
		
	}
	
	if(par_output_interval) aut.Output(output_file_path, gen);
	if(par_movie_interval) aut.Picture(pic_folder, gen);


	//randomszam generator lezarasa
	gsl_rng_free(r);
	
	cout << "Simulation ended" << endl;

	return 0;
} 
