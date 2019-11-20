#include <iostream>
#include "randomgen.h"
#include "stringreps.h"
#include "ca.h"


using namespace std;

class Vmi {
	public:
		int vmi1;
		double vmi2;
		
		Vmi() {    cout<<"Constructor Called"<<endl;   }
		~Vmi() {    cout<<"Deconstructor Called"<<endl;   }
};

gsl_rng * r;

int main(int argc, char *argv[]) {
	//randomszam generator inicializalasa
	time_t timer;
	r = (gsl_rng *) gsl_rng_alloc (gsl_rng_mt19937);
	gsl_rng_set(r, time(&timer));
/*
	int length =0, i=5;
	double *matrix;
	cadv::CellAut<int> aut(10, 10);
	
	cout << aut.nrow << " " << aut.size << endl;
*/	
	//length = cadv::grid_init(&matrix, 10, 10);
	
	//cout << cadv::Max(length, i) << endl;
	
	//cout << length << endl;
	//cout << *(matrix + 15) << endl;
	
	//delete [] (matrix);
	
	strrep::Strrep::set_length_activity(-0.005);
	//strrep::Strrep cella;
	
	/*
	cout << cella.length << "\t" << cella.R() << "\t" << cella.krepl << endl;
	for(int i=0; i<cella.length; i++) {
		cout << cella.seq[i] << " ";
	}
	cout << endl;
	*/
	
	cadv::CellAut<strrep::Strrep> aut(10,10);
	cout << strrep::Strrep::no_repl << endl;
	
	//aut.get(11)->update();
	
	//randomszam generator lezarasa
	gsl_rng_free(r);

	return 0;
} 
