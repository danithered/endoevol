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
	numbers[0] = numbers[1] = numbers[2] = numbers[3] = numbers[4] = numbers[5] = 0;
	complex = NULL;
	
	//generate sequence and det variables
	no_repl++;
	role=single;
	length = gsl_rng_uniform_int(r, MAXLEN/2);
	
	if(length) {
		//gen random sequence
		for(int b = 0; b < length; b++) {
			seq[b] = static_cast<strrep::bases>(gsl_rng_uniform_int(r, 4) + 1);
			seq_compl[b] = static_cast<strrep::bases>(gsl_rng_uniform_int(r, 4) + 1);
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

void strrep::Strrep::setNumsZero() {
	numbers_compl[1] = numbers_compl[2] = numbers_compl[3] = numbers_compl[4] = numbers_compl[5] = numbers[1] = numbers[2] = numbers[3] = numbers[4] = numbers[5] = 0;
}

int strrep::Strrep::align() {
//	std::cout << "align() called" << std::endl;
	double powsum;
	
	setNumsZero();
	if(length==0) {
		role=empty;
		return 0;
	}
	
	for(int b = 0; b < length; b++) {
		numbers[seq[b]]++;
		numbers_compl[seq_compl[b]]++;
	}
	
	if( R() + E() + C() + T() + B() + Rc() + Ec() + Cc() + Tc() + Bc() != 2 * length ) {
		std::cerr << "Something went wrong during alignment!" << std::endl << R() << "+" << E() << "+" << C() << "+" << T() << "+" << B() << "+" << Rc() << "+" << Ec() << "+" << Cc() << "+" << Cc() << "+" << Tc() << "+" << Bc() << " != 2x " << length << " of seq" << std::endl << getSeq() << std::endl;
		return(1);
	}
	
	role = single;
	
	//rates
	powsum = pow(R(), par_alpha) + pow(T(), par_alpha) + pow(E(), par_alpha) + pow(C(), par_alpha) + pow(B(), par_alpha);
	
	krepl = pow(R(), par_alpha) / powsum * length_activity[length];
	kendo = pow(E(), par_alpha) / powsum * length_activity[length];
	kasso_repl = pow(T(), par_alpha) / powsum * length_activity[length];
	kasso_endo = pow(C(), par_alpha) / powsum * length_activity[length];
//	std::cout << "align ended. role: " << role << ", length: " << length << ", krepl: " << krepl << ", kendo: " << kendo << ", kasso_repl: " << kasso_repl << ", kasso_endo: " << kasso_endo << std::endl << getSeq() << std::endl;	
	return 0;
}

void strrep::Strrep::del() {
	if(role != empty) no_repl--;
	setNumsZero();
	//for(int b = 0; b < MAXLEN; b++) {
	//	seq[b] = 0;
	//}
	seq[0] = N;
	seq_compl[0] = N;
	length = 0;
	role = empty;
	krepl = kendo = kasso_repl = kasso_endo = 0;
//	std::cout << "del" << std::endl;
	if(complex != NULL){
		complex->complex = NULL;
		complex->role = single;
		complex = NULL;
	}
//	std::cout << "etion" << std::endl;

}

void strrep::Strrep::set_length_activity(double beta) {
	for (int l = 0; l < MAXLEN; l++) {
		length_activity[l] = 1 - exp(beta * pow((double) l, par_exponent) )  ;
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
	
	//breaking complex
	diss();
	
	if(length < 2) {
		del();
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
			child->seq_compl[i-cut] = seq_compl[i];
			seq[i] = seq_compl[i] = N;
		}
		child->length = length - cut;
		length=cut;
	}
	else { // the first part is moved
//					printf("clevage: part1->child, part2->stays(parent)\n");
		for(i = 0; i < cut; i++) {
			child->seq[i] = seq[i];
			child->seq_compl[i] = seq_compl[i];
			//seq[i] = seq[i + cut];
		}
		for(; i < length; i++) {
			seq[i-cut] = seq[i];
			seq_compl[i-cut] = seq_compl[i];
			seq_compl[i] = seq[i] = N;
		}
		child->length = cut;
		length = length - cut;
	}
	
	align();
	child->align();
	no_repl++;
	
//	std::cout << "cleavage happened" << std::endl;	
	
}

int strrep::Strrep::Replication(Strrep* child) {
	int rna_site[MAXLEN]; /* We chose one of the elements of this
						array. The value stored in the chosen
						element is the site of mutation */

	int nmut=0; /* the number of mutations */
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
	
	//breaking complex - it is neccessary to do here as the function has multiple exit points
	diss();
//	std::cout << "complex broke" << std::endl;

	//deciding that the copy or the original stays in the parents place
	child->length = length;
	
	if(gsl_rng_uniform(r) < 0.5) { //the original stays
		original = this;
		copy = child;
//		std::cout << "Replication started. The original sequence stays" << std::endl;
	}
	else { //the copy stays
		copy_ptr = child->seq;
		child->seq = seq;
		seq = copy_ptr;
		
		original = child;
		copy = this;
		return_value++;
//		std::cout << "Replication started. The original sequence goes to child cell" << std::endl;
	}
	
	/* Replicate original to copy. We do mutation later. */
	for(pos_original = 0, pos_copy = 0; pos_original < original->length ; pos_original++){
		//is there a deletion
		if(gsl_rng_uniform(r) < par_deletion) {
			if(gsl_rng_uniform(r) < par_insertion) { //deletion and insertion
				if(pos_copy >= MAXLEN) {std::cerr << "ERROR (nonPerfectReplication): reached max length!" << std::endl; break;} //check if copy is too long
				copy->seq[pos_copy++] = static_cast<strrep::bases>(gsl_rng_uniform_int(r, 5) + 1);
//				std::cout << "insertion of char " << copy->seq[pos_copy -1] << std::endl;
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
					copy->seq[pos_copy++] = static_cast<strrep::bases>(gsl_rng_uniform_int(r, 5) + 1);
//					std::cout << "insertion of char " << copy->seq[pos_copy -1] << std::endl;
					copy->length++;
				}
				else { // ...to the left
					copy->seq[pos_copy++] = static_cast<strrep::bases>(gsl_rng_uniform_int(r, 5) + 1);
//					std::cout << "insertion of char " << copy->seq[pos_copy -1] << std::endl;
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
//	std::cout << "insertion/deletion ended" << std::endl;	
	
	if(return_value) original->align();
	
	//test if it is 0 length
	if(copy->length == 0) {
		copy->seq[0] = N;
		copy->role = empty;
		copy->length = 0;
		return(return_value);
	}
	
	no_repl++;
	
//	std::cout << "replication happened" << std::endl;		
	if(par_substitution) {
		nmut= gsl_ran_binomial(r, par_substitution, copy->length);
	}
	//nmut = (int) bnldev(par_substitution, copy->length);
	
	if(nmut==0) {
//		std::cout << "no substitution happened" << std::endl;
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
/*/		if(original_base == strrep::B){
			if(gsl_rng_uniform(r) < par_backmut){
				copy->seq[msite[i]] = static_cast<strrep::bases>(gsl_rng_uniform_int(r, 2)*2 + 1);
			}
		}
		else if(original_base == strrep::R) {
			copy->seq[msite[i]] = static_cast<strrep::bases>(gsl_rng_uniform_int(r, 2)*2 + 3);
		}
		else {
			copy->seq[msite[i]] = static_cast<strrep::bases>(gsl_rng_uniform_int(r, 2)*4 + 1);
		}
/*/		
		if(original_base == strrep::B) {
			if(gsl_rng_uniform(r) < par_backmut) {
				copy->seq[msite[i]] = static_cast<strrep::bases>( gsl_rng_uniform_int(r, 4) + 1  );
			}
		}
		else {
			mint = gsl_rng_uniform_int(r, 4) + 1 ; 
			original_base = (original_base + mint) % 5;
			copy->seq[msite[i]] = static_cast<strrep::bases>( original_base?original_base:5 );
		}
		
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
//		strrep::Strrep *co; co = y->complex;

		switch(y->role) {
			case repl: //replication from y's neighbour to x
//				std::cout << strrep::Strrep::no_repl << " maybe copiing " << x->role << y->role << y->complex->role << std::endl;
				if( gsl_rng_uniform(r) < y->krepl ) y->complex->Replication(x);
//				std::cout << strrep::Strrep::no_repl << " replication happened " << x->role << " " << y->role  << " " << co->role << std::endl;				
				break;
			case endo: //endonucleation from y's neighbour to x
//				std::cout << "maybe cleaving " << x->role << y->role << y->complex->role << std::endl;
				if( gsl_rng_uniform(r) < y->kendo ) y->complex->Clevage(x);
//				std::cout << "cleavage happened " << x->role << " " << y->role << " " << co->role << std::endl;
				break;
			case endo_template: //endonucleation from y to x
//				std::cout << "maybe cleaving " << x->role << " " << y->role << " " << y->complex->role << std::endl;
				if( gsl_rng_uniform(r) < y->complex->kendo ) y->Clevage(x);
//				std::cout << "cleavage happened " << x->role << " " << y->role << " " << co->role << std::endl;
				break;
			case repl_template: //replication from y to x
//				std::cout << strrep::Strrep::no_repl << " maybe copiing " << x->role << y->role << y->complex->role << std::endl;
				if( gsl_rng_uniform(r) < y->complex->krepl ) y->Replication(x);
//				std::cout << strrep::Strrep::no_repl << " replication happened " << x->role << " " << y->role << " " << co->role << std::endl;
		}
	}
	else { // x is not empty
		if( x->role == single) { // x is a single molecule
			if( y->role == single) { //x and y are single molecules
				//complex formation
				*x + y;
//				std::cout << "Comlex formation" << std::endl;
			}
		}
		else { // x is a complex molecule
			//complex dissotiation
			if( gsl_rng_uniform(r) < par_dissotiation ) x->diss();
			//x->dissotiation(gsl_rng_uniform(r));
//			std::cout << "Comlex dissotiation" << std::endl;
		}
	}
}

std::string strrep::Strrep::getRole(){
	switch(role){
		case strrep::empty: return("empty\0");
		case strrep::single: return("single\0");
		case strrep::repl: return("repl\0");
		case strrep::endo: return("endo\0");
		case strrep::repl_template: return("repl_template\0");
		case strrep::endo_template: return("endo_template\0");
		default: return("NA\0");
	}
}

std::string strrep::Strrep::getSeq(){
	std::string charseq;
	static char basechars[] = "NRETCB";
	
	if(role == empty) return( (std::string) "N" );
	
	charseq = basechars[(int)seq[0]];
	for(int b = 1; b < length; b++) {
		charseq = charseq + basechars[(int)seq[b]];
	}
	
	return(charseq);
}

std::string strrep::Strrep::getComplSeq(){
	std::string charseq;
	static char basechars[] = "NRETCB";
	
	if(role == empty) return( (std::string) "N" );
	
	charseq = basechars[(int)seq_compl[0]];
	for(int b = 1; b < length; b++) {
		charseq = charseq + basechars[(int)seq_compl[b]];
	}
	
	return(charseq);
}

void strrep::Strrep::setSeq(char* charseq, char* charseq_compl){
	del();
	for(length = 0; charseq[length] != '\0' && charseq_compl[length] != '\0'; length++){
		if(length > MAXLEN){
			std::cerr << "too long initial sequence! Truncated" << std::endl;
			break;
		}
		switch(charseq[length]){
			case 'R':
				seq[length]=strrep::R;
				break;
			case 'E':
				seq[length]=strrep::E;
				break;
			case 'T':
				seq[length]=strrep::T;
				break;
			case 'C':
				seq[length]=strrep::C;
				break;
			case 'B':
				seq[length]=strrep::B;
				break;
			default:
				std::cerr << "ERROR: setSeq: non regular character found!" <<std::endl;
				return;
		}
		switch(charseq_compl[length]){
			case 'R':
				seq_compl[length]=strrep::R;
				break;
			case 'E':
				seq_compl[length]=strrep::E;
				break;
			case 'T':
				seq_compl[length]=strrep::T;
				break;
			case 'C':
				seq_compl[length]=strrep::C;
				break;
			case 'B':
				seq_compl[length]=strrep::B;
				break;
			default:
				std::cerr << "ERROR: setSeq: non regular complementer character found!" <<std::endl;
				return;
		}
	}
	if(length) no_repl++;
	align();
}

int strrep::Sca::Output(char* filepath, int time){
	// open a file in write mode.
	static int output_open=0;
	static std::fstream output, mean_output;
	static char filename[255] = "\0", mean_filename[255] = "\0", rolechars[]="NSRETC"; 
	
	int numbers[6]={0,0,0,0,0,0}, i=0;
	double length[6]={0.0,0.0,0.0,0.0,0.0,0.0}, kR[6]={0.0,0.0,0.0,0.0,0.0,0.0}, kE[6]={0.0,0.0,0.0,0.0,0.0,0.0}, kT[6]={0.0,0.0,0.0,0.0,0.0,0.0}, kC[6]={0.0,0.0,0.0,0.0,0.0,0.0};
	
//	std::cout << filepath << std::endl;	
	
	if(!output_open) {
		sprintf(filename, "%soutput.txt\0", filepath);		
		sprintf(mean_filename, "%satlagadat.txt\0", filepath);		
		
		output.open(filename, std::fstream::out | std::fstream::app);
		mean_output.open(mean_filename, std::fstream::out | std::fstream::app);
		
		output_open++;
		mean_output << "time\tnumber\trole\tlength\tkR\tkE\tkT\tkC" << std::endl;
/**/		std::cout << "Output opened: " << filename << std::endl;
		//mean_output.open(filename, std::fstream::out | std::fstream::app);
	}
	
	for(i = 0; i < size; i++ ){
		if(matrix[i].role) {
			output << "t_" << time << "\t" << i << "\t"  << matrix[i].length << "\t" << matrix[i].getRole() << "\t" << matrix[i].getSeq() << "\t" << matrix[i].krepl << "\t" << matrix[i].kendo << "\t" << matrix[i].kasso_repl << "\t" << matrix[i].kasso_endo << '\t' << matrix[i].getComplSeq() << std::endl;
			numbers[(int)matrix[i].role]++;
			length[(int)matrix[i].role] += matrix[i].length;
			kR[(int)matrix[i].role] += matrix[i].krepl;
			kE[(int)matrix[i].role] += matrix[i].kendo;
			kT[(int)matrix[i].role] += matrix[i].kasso_repl;
			kC[(int)matrix[i].role] += matrix[i].kasso_endo;
		}
	}
	
	for(i=1; i<6; i++){
		length[i] = numbers[i]?length[i]/numbers[i]:0;
		kR[i] = numbers[i]?kR[i]/numbers[i]:0;
		kE[i] = numbers[i]?kE[i]/numbers[i]:0;
		kT[i] = numbers[i]?kT[i]/numbers[i]:0;
		kC[i] = numbers[i]?kC[i]/numbers[i]:0;
		
		mean_output << time << "\t" << numbers[i] << "\t" << rolechars[i] << "\t" << length[i] << "\t" << kR[i] << "\t" << kE[i] << "\t" << kT[i] << "\t" << kC[i] << std::endl;
	}

	output.flush();
	mean_output.flush();
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
