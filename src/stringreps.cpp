#include <iostream>
#include <cmath>
#include "stringreps.h"
#include "randomgen.h"

int strrep::Strrep::no_repl = 0;
double strrep::Strrep::length_activity[MAXLEN];

strrep::Strrep::Strrep() {
	//initialise arrays
	seq = new bases[MAXLEN];
	
	//init vars
	numbers[0] = numbers[1] = numbers[2] = numbers[3] = numbers[4] = 0;
	complex = NULL;
	
	//generate sequence and det variables
	no_repl++;
	length = gsl_rng_uniform_int(r, MAXLEN);
	
	if(length) {
		//gen random sequence
		for(int b = 0; b < length; b++) {
			seq[b] = static_cast<strrep::bases>(gsl_rng_uniform_int(r, 4) + 1);
		}
		
		if ( align() ) {
			std::cerr << "Something went wrong initialising a random sequence" << std::endl;
			del();
		} 
	}
	else { //if it is 0 length seq
		del();
	}
}

strrep::Strrep::~Strrep() {
	delete [] (seq);
	
}

int strrep::Strrep::align() {
	numbers[1] = numbers[2] = numbers[3] = numbers[4] = 0;
	for(int b = 0; b < length; b++) {
		numbers[seq[b]]++;
	}
	
	if( R() + E() + C() + T() != length ) {
		std::cerr << "Something went wrong during alignment!" << std::endl;
		return(1);
	}
	
	role = single;
	
	//rates
	krepl = R() * length_activity[length];
	kendo = E() * length_activity[length];
	kasso_repl = T() * length_activity[length];
	kasso_endo = C() * length_activity[length];
	
	return 0;
}

void strrep::Strrep::del() {
	no_repl--;
	numbers[1] = numbers[2] = numbers[3] = numbers[4] = 0;
	//for(int b = 0; b < MAXLEN; b++) {
	//	seq[b] = 0;
	//}
	seq[0] = N;
	length = 0;
	role = empty;
	krepl = kendo = kasso_repl = kasso_endo = 0;
	if(complex != 0){
		complex->complex = NULL;
		complex =NULL;
	}

}

void strrep::Strrep::set_length_activity(double beta) {
	for (int l = 0; l < MAXLEN; l++) {
		length_activity[l] = (1 - exp(beta * l)) / l;
	}
}

