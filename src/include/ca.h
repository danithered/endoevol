#ifndef _CADV_
#define _CADV_

#include <vector>
#include "randomgen.h"
#include "dv_tools.h"

#define MOORE_NEIGH 1
#define VON_NEUMANN_NEIGH 0
#define MARGOLUS_NEIGH -2

namespace cadv {
	enum Ca_layout {empty=0, single=1, line=2, hex=6, square = 4};
	
	int torus(int m, int k);
	int* neighInic(int cellaszam_f, int ncol_f, Ca_layout layout, double neigh_tipus_f);
	int szomsz_meret(Ca_layout layout, double tipus_f);

	
	template <typename matrix_type>
	int grid_init(matrix_type **plane, int size1=300, int size2=300, Ca_layout layout=square) {
		int size = 0;
		
		if(size1 <=0 || size2 <= 0) {
			size =0;
			return(0);
		} 
		else if(size1 == 1 ) {
			size = size2;
		}
		else if(size2 == 1){
			size == size1;
		}
		else switch (layout) { //the matrix is definitely more than 1D based on sizes
				case square:
					size = size1 * size2;
		}
		
		if(layout == square) {
			*plane = new matrix_type[size];
			if(! *plane ) {
				std::cerr << "ERROR: cadv::ca_init: plane could not be initialized!" << std::endl;
				return(0);
			}
		}
		
		return(size);
	}
	
	template <class celltype>
	class CellAut {
		public:
			int nrow;
			int ncol;
			int size;
			cadv::Ca_layout layout;
			
			celltype *matrix;
			std::vector<int*> neighbourhoods;
			std::vector<int> no_neighbours;
			
			//Constructor 1
			CellAut(int size1=300, int size2=300, cadv::Ca_layout layout_type=square){
				nrow=size1;
				ncol=size2;
				
				size = cadv::grid_init(&matrix, size1, size2, layout_type);
				
				if(!size) layout = empty;
				if(size1==1 || size2==1){
					if(size1==size2) {
						layout = single;
					}
					else {
						layout = line;
					}
				}
			
/**/				std::cout << "Basic Constructor Called" << std::endl;
			}
			
			//Constructor 2
			CellAut(int size1, int size2, cadv::Ca_layout layout_type, celltype* pool, double* probs, int no_choices){
				CellAut(size1, size2, layout_type);
				init(pool, probs, no_choices);
			}
			
			//Deconstructor
			~CellAut(){
				if(size) delete [] (matrix);
				for(int i =0; i < neighbourhoods.size(); i++) {
					delete [] (neighbourhoods.at(i));
				}
/**/				std::cout << "Deconstructor Called" << std::endl;
			}
			
			virtual void Update(int cell);
			
			//Functions
			
			
			//Positions
			
			///gives back x pos of cell
			
			///gives back y pos of cell
			
			///finds cell in pos [x,y] and gives back a pointer to it
			
			///gives back pointer to cell
			inline celltype* get(int cell) {
				return(matrix + cell);
			}
			
			///initialises matrix with predefined values, randomly
			void init(celltype* pool, double* probs, int no_choices) {
				int i = 0;
				double sum = 0.0;
				
				for(i=0; i < no_choices; i++) {
					sum += probs[i];
				}
				for(i=0; i < size; i++) {
					matrix[i] = pool[dvtools::brokenStickVals(probs, no_choices, sum, gsl_rng_uniform(r) )];
				}
			}
			
			///adds neighbourhood matrix
			int neighAdd(double type) {
				int* new_neigh;
				new_neigh = neighInic(size, ncol, layout, type);
				if( !new_neigh ) {
					std::cerr << "failed to create neighbourhood" << std::endl;
					return(-1);
				}
				neighbourhoods.push_back( new_neigh );
				no_neighbours.push_back( szomsz_meret(layout, type) );
				return( neighbourhoods.size() - 1 );
			}
			
			///gives back number of neighbours in a neighbourhood
			inline int neighNo(int no_nm = 0) {
				return no_neighbours[no_nm];
			}
			
			///gives back pointer to neighbourhood matrix
			inline int* neigh(int no_nm = 0) {
				return neighbourhoods[no_nm];
			}
			
			///gives back pointer to neighbours of a cell
			inline int* neigh(int cell, int no_nm) {
				return neighbourhoods[no_nm] + no_neighbours[no_nm]*cell;
			}
			
			///gives back pos of neighbour no_neigh of cell
			inline int neigh(int cell, int no_neigh, int no_nm) {
				return neighbourhoods[no_nm] + no_neighbours[no_nm]*cell + no_neigh;
			}
			
			///gives back pointer to a random neighbour except the cell itself
			inline celltype* rneigh(int cell, int no_nm=0) {
				return (matrix + *(neighbourhoods[no_nm] + no_neighbours[no_nm]*cell + gsl_rng_uniform_int(r, no_neighbours[no_nm] - 1) + 1) );
			}
			
			///swaps two cells in the matrix
			
			///copies over a value on an other
			
			//Diffusions
			///Random walk step
			
			///Random walk
			
			///Flipping step
			
			///Flipping
			
			///TM inic
			
			///TM step
			
			///TM async
			
			///TM async updorder
			
			///TM classic
	};
	
}


#endif
