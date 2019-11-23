#ifndef _STRREP_
#define _STRREP_

#include <fstream>
extern "C" {
#include "mypng.h"
}
#include "randomgen.h"
#include "parameters.h"
#include "dv_tools.h"
#include "ca.h"
#include "mycol.h"


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
			
			//dissociation
			void diss();
			
			//test for dissotioation
			void dissotiation(double rand);
			
			//get sequence as string
			std::string getSeq();
			
			//set sequence from character string
			void setSeq(char *charseq);
			
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
			
			//complex formation
			void operator +( Strrep* target){
				double kassos[5] = {par_k_noasso, kasso_repl, kasso_endo, target->kasso_repl, target->kasso_endo};
				
				//roles
				switch( dvtools::brokenStickVals(kassos, 5, -1, gsl_rng_uniform(r)) ) {
					case 0:
						return ;
					case 1: // this connects to target as replicator
						role=repl;
						target->role=repl_template;
						break;
					case 2: // this connects to target as endonuclease
						role=endo;
						target->role=endo_template;
						break;
					case 3: // target connects to this as replicator
						role=repl_template;
						target->role=repl;
						break;
					case 4: // target connects to this as endonuclease
						role=endo_template;
						target->role=endo;
						break;
					default: 
						std::cerr << "ERROR: brokenStickVals gave back unvalid choice in complex formation" << std::endl;
						return ;
				}
				
				//complex pointers
				complex = target;
				target->complex = this;
				
			}
			
			//endonucleation
			void Clevage(Strrep* child);

			//replication
			int Replication(Strrep* child);
			
	};

	class Sca : public cadv::CellAut<Strrep> {
		public:
			void Update(int cell);
			
			int Output(std::string filename, int time);
			
			Sca(int size1=300, int size2=300, cadv::Ca_layout layout_type = cadv::Ca_layout::square){
				nrow=size1;
				ncol=size2;
				
				strrep::Strrep::no_repl = 0;
				
				size = cadv::grid_init(&matrix, size1, size2, layout_type);
				
				if(!size) layout = cadv::Ca_layout::empty;
				if(size1==1 || size2==1){
					if(size1==size2) {
						layout = cadv::Ca_layout::single;
					}
					else {
						layout = cadv::Ca_layout::line;
					}
				}
			}
			Sca(int size1, int size2, cadv::Ca_layout layout_type, Strrep* pool, double* probs, int no_choices){
				CellAut(size1, size2, layout_type);
				init(pool, probs, no_choices);
			}
			
			int Picture(char* folder, int timestep);
	};
	
	
}
#endif 
