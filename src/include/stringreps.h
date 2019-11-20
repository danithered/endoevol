#ifndef _STRREP_
#define _STRREP_

#include "randomgen.h"
#define MAXLEN 300

namespace strrep {
	
	enum repl_types {empty=0, single=1, repl=2, endo=3, repl_template=4, endo_template=5};
	enum bases {N=0, R=1, E=2, T=3, C=4};
	
	class Strrep {
		public:
			//// VARIABLES
			// for everybody
			static int no_repl;
			static double length_activity[MAXLEN]; // (1 - exp(beta * length)) / length for length 1 to MAXLEN
			
			//for individual cells
			bases* seq;
			
			int numbers[5];		// the number of different bases in the sequence. N:0, R:1, E:2, T:3, C:4
			//int R;
			//int E;
			//int T;
			//int C;
			int length; 			// the length of the sequence
			
			Strrep* complex; 			// who is in complex with this replicator
			
			repl_types role;
			
			double krepl;
			double kendo;
			double kasso_repl;
			double kasso_endo;
			
			//// CONSTRUCTORS
			//Constructor empty
			Strrep();
			
			//Constructor seq
			//Strrep(char* init_seq, int init_seq_length);
			
			//Deconstructor
			~Strrep();
			
			//// FUNCTIONS
			//set length_activity;
			static void set_length_activity(double beta);
			
			//align seq
			int align();
			
			//kill cell, make it empty
			void del();
			
			//get number of bases
			inline int R() {return(numbers[1]);};
			inline int E() {return(numbers[2]);};
			inline int T() {return(numbers[3]);};
			inline int C() {return(numbers[4]);};
			
			//random walk
			void operator >( Strrep* target){
				if( target->role != empty ) std::cerr << "ERROR: Strep: during diffusion non-empty target given, target will be overwritten!" << std::endl; 
				*target = *this;
				if(complex == NULL) {
					this->del();
				}
				else {
					*this=*complex; //pull complex to this position
					target->complex->del(); //delete complex
					complex=target; //this cell's new complex is target
					target->complex = this; //target's new complex is this cell'
				}
			}
	};

	
	
	
}
#endif 
