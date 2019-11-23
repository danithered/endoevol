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
	length = gsl_rng_uniform_int(r, MAXLEN/2);
	
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
	if(complex != NULL){
		complex->complex = NULL;
		complex->role = single;
		complex = NULL;
	}

}

void strrep::Strrep::set_length_activity(double beta) {
	for (int l = 0; l < MAXLEN; l++) {
		length_activity[l] = (1 - exp(beta * l)) / l;
	}
}

void strrep::Strrep::diss() {
	//dissoc partner
	complex->complex = NULL;
	if(complex->role != empty) complex->role = single;
	
	//dissoc cell
	complex = NULL;
	if(role != empty) role = single;
	
}

void strrep::Strrep::dissotiation(double rand){
	switch( role ) {
		case repl_template:
			if(rand < 1 - kasso_repl ) diss();
			break;
		case endo_template:
			if(rand < 1 - kasso_endo ) diss();
			break;
		case repl:
			if(rand < 1 - complex->kasso_repl ) diss();
			break;
		case endo:
			if(rand < 1 - complex->kasso_endo ) diss();
			break;
		default:
			std::cerr << "ERROR: attempted dissotiation on a non-complex replicator of role: " << role << std::endl;
			break;
	}
}

void strrep::Strrep::Clevage(Strrep* child) {
	int i = 0, till=0, cut=0;
	
	if(length < 2) {
		seq[0] = N;
		length = 0;
		role = empty;
		return ;
	}
	
	//determine place of cleavage
	if( (length % 2) == 0 ) cut = length/2; //if it is an even length
	else { //if it is uneven length
		cut = ( gsl_rng_uniform(r) < 0.5 ) ? (length-1)/2 : (length-1)/2 + 1;
	}

	//decide which part is being moved to the child cell
	if( gsl_rng_uniform(r) < 0.5 ) { //if the second part is moved
		//for seq
//					printf("clevage: part1->stays, part2->child\n");
		for(i = cut; i < length; i++) {
			child->seq[i-cut] = seq[i];
			seq[i] = N;
		}
		child->length = length - cut;
		length=cut;
	}
	else { // the first part is moved
//					printf("clevage: part1->child, part2->stays(parent)\n");
		for(i = 0; i < cut; i++) {
			child->seq[i] = seq[i];
			//seq[i] = seq[i + cut];
		}
		for(; i < length; i++) {
			seq[i-cut] = seq[i];
			seq[i] = N;
		}
		child->length = cut;
		length = length - cut;
	}
	
	align();
	child->align();
	
}

int strrep::Strrep::Replication(Strrep* child) {
	int rna_site[MAXLEN]; /* We chose one of the elements of this
						array. The value stored in the chosen
						element is the site of mutation */

	int nmut; /* the number of mutations */
	int end; /* These represent the range of the site numbers of RNA
			string from which we chose sites mutating */
	int msite[MAXLEN]; /* A site number where we will insert a
				point mutation */
	int index; /* Index of 'rna_site'. The element of corresponding
			element stores the position of the
			mutation. Index is determined by calling random
			number generator. */
	int mint; /* It will take 1,2,3 uniformely randomly. Think
			about a room such that A->C->U->T->A(repeat). We
			incliment this room from the original character by
			'mint'. This means that we rondomly chose one of
			the 3 characters */
	int i, return_value=0, original_base;
	strrep::Strrep *original, *copy;
	int pos_original, pos_copy;
	
	bases *copy_ptr;

	//deciding that the copy or the original stays in the parents place
	child->length = length;
	if(gsl_rng_uniform(r) < 0.5) { //the original stays
		original = this;
		copy = child;
//		printf("Replication started. The original sequence stays\n%d\t%s\n", (int)strlen(original), original);
	}
	else { //the copy stays
		copy_ptr = child->seq;
		child->seq = seq;
		seq = copy_ptr;
		
		original = child;
		copy = this;
		return_value++;
//		printf("Replication started. The original sequence goes to child cell\n%d\t%s\n", (int) strlen(original), original);
	}
	
	/* Replicate original to copy. We do mutation later. */
	for(pos_original = length-1, pos_copy = 0; pos_original >= 0; pos_original--){
		//is there a deletion
		if(gsl_rng_uniform(r) < par_deletion) {
			if(gsl_rng_uniform(r) < par_insertion) { //deletion and insertion
				if(pos_copy >= MAXLEN) {std::cerr << "ERROR (nonPerfectReplication): reached max length!" << std::endl; break;} //check if copy is too long
				copy->seq[pos_copy++] = static_cast<strrep::bases>(gsl_rng_uniform_int(r, 4) + 1);
			}
			else {
				copy->length--;
			} //length decreases: only deletion
		}
		else {
			if( gsl_rng_uniform(r) < par_insertion) { //length increases: only insertion...
				if(pos_copy + 1 >= MAXLEN) {std::cerr << "ERROR (nonPerfectReplication): reached max length!!" << std::endl; break;} //check if copy is too long
				if( gsl_rng_uniform(r) < 0.5) { // ...to the right
					copy->seq[pos_copy++] = original->seq[pos_original];
//					if(copy[pos_copy-1] == '\0') printf("ERROR: nonPerfectReplication: RNAc2cc has found a non RNA charaster (%c) during inserting right (%d)!\n%d\t%s\n%d\t%s\n", original[pos_original], length, (int)strlen(original), original, (int)strlen(copy), copy);
					copy->seq[pos_copy++] = static_cast<strrep::bases>(gsl_rng_uniform_int(r, 4) + 1);
					copy->length++;
				}
				else { // ...to the left
					copy->seq[pos_copy++] = static_cast<strrep::bases>(gsl_rng_uniform_int(r, 4) + 1);
					copy->seq[pos_copy++] = original->seq[pos_original];
//					if(copy[pos_copy-1] == '\0') printf("ERROR: nonPerfectReplication: RNAc2cc has found a non RNA charaster (%c) during inserting left (%d)!\n%d\t%s\n%d\t%s\n", original[pos_original], length, (int)strlen(original), original, (int)strlen(copy), copy);
					copy->length++;
				}
			} else { //basic copying
				if(pos_copy >= MAXLEN) {std::cerr << "ERROR (nonPerfectReplication): reached max length!!!" << std::endl; break;} //check if copy is too long
				copy->seq[pos_copy++] = original->seq[pos_original];
//				if(copy[pos_copy-1] == '\0') printf("ERROR: nonPerfectReplication: RNAc2cc has found a non RNA charaster (%c) during basic copying (%d)!\n%d\t%s\n%d\t%s\n", original[pos_original], length, (int)strlen(original), original, (int)strlen(copy), copy);
			}
		}
//		if(original[pos_original] == 'N') printf("ERROR: nonPerfectReplication: pos_original is %c!\n\t%s\n\t%s\n", original[pos_original], original, copy);
	}
	
	//test if it is 0 length
	if(copy->length == 0) {
		copy->seq[0] = N;
		copy->role = empty;
		copy->length = 0;
		return(return_value);
	}
		
		
	if(par_substitution == 0.) {
		copy->align();
		return return_value;
	}

	/* Here, we calculate how many mutations occurs */
	//nmut = (int) bnldev(par_substitution, copy->length);
	nmut= gsl_ran_binomial(r, par_substitution, copy->length);

	if(nmut==0) {
		copy->align();
		return return_value;
	}

	/* Before choosing where mutations occur, initialize
		'rna_site'. We use this to chose where mutations occur. */
	for(i = 0; i < copy->length; i++) rna_site[i]=i;

	/* end will be decremented each time we chose a mutation
		site. */
	end = copy->length;

	/* We chose the sites of mutation */
	for(i=0;i<nmut;i++){
		/* Chose a site */
		/* (the fractional part is truncated) */
		index = gsl_rng_uniform(r) * (double)(end);

		msite[i] = rna_site[index];

		/* We put the site of mutation in the last element of
			"rna_site" into the choice range, i.e., rna_site[end-1]
			into rna_site[index] in order that we do not chose the
			site in rna_site[index] twice and that we can chose a site
			in rna_site[end-1] later. */
		rna_site[index] = rna_site[end-1];

		/* decrement "end". We should not chose rna_site[end-1] again
			where "end" is before the decrement because it has been
			already chosen. After the decrement, rna_site[end-1] can
			be chosen. */
		end--;
	}

	/* We decide what kind of mutation occurs at chosen sites where
		mutation will occur. Then substitute the bases in the copy
		sequence. */
	
	for(i=0;i<nmut;i++){
		original_base = (int) copy->seq[msite[i]];
		mint = gsl_rng_uniform_int(r, 3) + 1 ; 
		original_base = (original_base + mint) % 4;
		copy->seq[msite[i]] = static_cast<strrep::bases>( original_base?original_base:4 );
	}
	
	copy->align();
	return return_value;
}

void strrep::Sca::Update(int cell) {
	strrep::Strrep *x, *y;
	
	x = get(cell);
	y = rneigh(cell);
	
//	std::cout << "update " << cell <<std::endl;

	if(x->role == empty){ // x is empty
		//enzymatic reaction
		switch(y->role) {
			case repl: //replication from y's neighbour to x
				if( gsl_rng_uniform(r) < y->krepl ) y->complex->Replication(x);
				break;
			case endo: //endonucleation from y's neighbour to x
				if( gsl_rng_uniform(r) < y->kendo ) y->complex->Clevage(x);
				break;
			case endo_template: //endonucleation from y to x
				if( gsl_rng_uniform(r) < y->complex->kendo ) y->Clevage(x);
				break;
			case repl_template: //replication from y to x
				if( gsl_rng_uniform(r) < y->complex->krepl ) y->Replication(x);
		}
	}
	else { // x is not empty
		if( x->role == single) { // x is a single molecule
			if( y->role == single) { //x and y are single molecules
				//complex formation
				*x + y;
			}
		}
		else { // x is a complex molecule
			//complex dissotiation
			x->dissotiation(gsl_rng_uniform(r));
		}
	}
}

std::string strrep::Strrep::getSeq(){
	std::string charseq;
	static char basechars[] = "NRETC";
	
	if(role == empty) return( (std::string) "N" );
	
	charseq = basechars[(int)seq[0]];
	for(int b = 1; b < length; b++) {
		charseq = charseq + basechars[(int)seq[b]];
	}
	
	return(charseq);
}

int strrep::Sca::Output(std::string filename, int time){
	// open a file in write mode.
	static int output_open=0;
	static std::fstream output;
	if(!output_open) {
//		std::cout << "Output opened" << std::endl;		
		output_open++;
		output.open(filename, std::fstream::out | std::fstream::app);
	}
	
	for(int i = 0; i < size; i++ ){
		output << time << "\t" << i << "\t" << matrix[i].role << "\t" << matrix[i].getSeq() << "\t" << matrix[i].krepl << "\t" << matrix[i].kendo << "\t" << matrix[i].kasso_repl << "\t" << matrix[i].kasso_endo << std::endl;
	}

	output.flush();
	return(0);
}


int strrep::Sca::Picture(char* folder, int timestep) {
	unsigned char *data;

	data = new unsigned char [3*ncol*nrow];
	
	for (int i=0; i < size; i++) {
//		std::cout << i <<std::endl;
		data[ i*3 ] = userCol[(int) matrix[i].role][0];
		data[ i*3 + 1] = userCol[(int) matrix[i].role][1];
		data[ i*3 + 2] = userCol[(int) matrix[i].role][2];
//		std::cout << "pixel " << i << ": " << data[ i*3 ] << " " << data[ i*3 +1] << " " << data[ i*3 +2] << std::endl;	
	}
//	std::cout << "here" << std::endl;
	PlanePNG(folder, timestep, ncol, nrow, size, data);
//	std::cout << "picture captured" << std::endl;
	delete [] (data);
	return(0);
}
