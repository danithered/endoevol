#include <iostream>
#include <cmath>
#include "ca.h"

namespace cadv {

	//feltolti a szomszedsagmatrixot: sajat pozicio, N, W, E, S
	int* neighInic(int cellaszam_f, int ncol_f, Ca_layout layout, double neigh_tipus_f) {
		/*
		* cellaszam_f: az alapmatrix cellaszama
		* neigh_tipus_f: a szomsz. merete: ha 0 -> von Neumann; ha több -> Moore szomsz
		* */
		
		//Valtozok deklaralasa
		int metNeigh_matrixmeret_f=0, alapMCellaszamlalo_f=0, colszamlalo_f=0, rowszamlalo_f=0, vonN_f = 0, szomsz_f=0, szomszedszam_f=0;
		static int *matrix_f;
		/*
		* metNeigh_matrixmeret_f: a szomsz. matrix cellaszama
		* alapMCellaszamlalo_f: szamlalo: hanyadik alapmatrix cellanal jarunk a feltoltesben
		* ncol_f: alapmatrix oszlopainak szama =sqrt(cellaszam_f)
		* colszamlalo_f: hanyadik oszlopban tartunk
		* rowszamlalo_f: hanyadik sorban tartunk
		* szomsz_f: egy cella szomszedsagaba hany cella tartozik pl.: 2es szomsz eseten (2*2+1)^2=25
		* *matrix_f: a szomsz. matrix
		* */

		//initialise neighbourhood matrix
		szomsz_f = szomsz_meret(layout, neigh_tipus_f);
		metNeigh_matrixmeret_f=cellaszam_f* szomsz_f; //szomsz. matrix cellameretenek meghatarozasa
		matrix_f=(int *) calloc(metNeigh_matrixmeret_f, sizeof(int)); //szomsz. matrix lefoglalasa
		if(!matrix_f) {
			std::cerr << "cadv::neighInic: ERROR: could not initialise neighbourhood matrix!" << std::endl;
			return(NULL);
		}		
		
		//create neighbourhood matrix
		switch(layout) {
			case square:
				if(neigh_tipus_f >= 0) { //vonNeumann, Moore and hybrid neighs
					if(neigh_tipus_f==0 || ( neigh_tipus_f - (int)neigh_tipus_f )){
						vonN_f = 1;
						neigh_tipus_f = (double)(int)neigh_tipus_f;
					}
						
					for(alapMCellaszamlalo_f=0; alapMCellaszamlalo_f<cellaszam_f; alapMCellaszamlalo_f++) {
						szomszedszam_f=0;					
						for(rowszamlalo_f = 0 - neigh_tipus_f; rowszamlalo_f <= 0 + neigh_tipus_f; rowszamlalo_f++) {
							for(colszamlalo_f = 0 - neigh_tipus_f; colszamlalo_f <= 0 + neigh_tipus_f; colszamlalo_f++) {
								//szomszedszam_f++;
								if(!colszamlalo_f && !rowszamlalo_f) *(matrix_f+ alapMCellaszamlalo_f*szomsz_f)=torus(cellaszam_f, alapMCellaszamlalo_f+rowszamlalo_f*ncol_f+colszamlalo_f);
								else {
									*(matrix_f+ alapMCellaszamlalo_f*szomsz_f+szomszedszam_f+1)=torus(cellaszam_f, alapMCellaszamlalo_f+rowszamlalo_f*ncol_f+colszamlalo_f);
									szomszedszam_f++;
								}
							}
						}
						if(vonN_f) {
							*(matrix_f + alapMCellaszamlalo_f*szomsz_f + szomszedszam_f + 1) = torus(cellaszam_f, alapMCellaszamlalo_f - (1 + neigh_tipus_f)*ncol_f); //1: N
							*(matrix_f + alapMCellaszamlalo_f*szomsz_f + szomszedszam_f + 2) = torus(cellaszam_f, alapMCellaszamlalo_f - (1 + neigh_tipus_f)); //2: W
							*(matrix_f + alapMCellaszamlalo_f*szomsz_f + szomszedszam_f + 3) = torus(cellaszam_f, alapMCellaszamlalo_f + (1 + neigh_tipus_f)); //3: E
							*(matrix_f + alapMCellaszamlalo_f*szomsz_f + szomszedszam_f + 4) = torus(cellaszam_f, alapMCellaszamlalo_f + (1 + neigh_tipus_f)*ncol_f); //4: S
						}
					}
				}
				else { //Margolus type neighs (only 2x2)
					for(alapMCellaszamlalo_f=0; alapMCellaszamlalo_f<cellaszam_f; alapMCellaszamlalo_f++) {
						*(matrix_f + alapMCellaszamlalo_f*szomsz_f) = alapMCellaszamlalo_f;
						*(matrix_f + alapMCellaszamlalo_f*szomsz_f + 1) = torus(cellaszam_f, alapMCellaszamlalo_f + 1); //jobbra
						*(matrix_f + alapMCellaszamlalo_f*szomsz_f + 2) = torus(cellaszam_f, alapMCellaszamlalo_f + ncol_f); //lent
						*(matrix_f + alapMCellaszamlalo_f*szomsz_f + 3) = torus(cellaszam_f, alapMCellaszamlalo_f + ncol_f + 1); //jobbra le
					}	
				}
				break;
			default:
				std::cerr << "cadv::neighInic: ERROR: did not recognize layout" << std::endl;
		} //switch (layout)
		
		return (matrix_f);
	}

	int szomsz_meret(Ca_layout layout, double tipus_f) {
		switch(layout) {
			case square:
				if(tipus_f == 0) return(5); //vonNeumann
				else if(tipus_f < 0) return(tipus_f * tipus_f);
				else {
					if(tipus_f - (int)tipus_f ) return(szomsz_meret(layout, (double)(int)tipus_f) + 4);
					else return (std::pow(((int)tipus_f*2+1), 2)); //Moore szomszedsag: pl.: 2es szomsz eseten (2*2+1)^2=25 
				}
				break;
			default:
				std::cerr << "cadv::szomsz_meret: ERROR: could not recognize layout" << std::endl;
		}
				return(0);
	}
	
	int torus(int m, int k)	{
		/*torus*/
		/** *******************************************************
		k: it is a grid point that is determined by torus()
		m: size of matrix 
		
		by: Tamas Czaran
		***********************************************************/
		return ((m-((m-k)%m)))%m;
	}/*torus*/


} //namespace
