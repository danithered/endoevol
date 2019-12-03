#include "parameters.h" 

double par_length_dependence = -0.05;
double par_death = 0.1;
double par_substitution = 0;
double par_insertion = 0;
double par_deletion = 0;
double par_k_noasso = 0.1;
//double par_decay_rate = 1;
double par_diffusion_rate = 1;

int par_maxtime = 100000;
int par_ncol = 512;
int par_nrow = 512;
int par_movie_interval = 10;
int par_output_interval = 100;

char par_pic_folder[255] = "movie\0";
char par_ID[255] = "test\0";

char **par_init_seqs;
double *par_init_props;
int par_init_no=0;

//output parameters to file
int paramsToFile(char* filename){
	int i=0;
	
	// open a file in write mode.
	std::fstream paramfile(filename, std::fstream::out);
	
//	std::cout << "Printing parameters to file: " << filename << std::endl;	
	
	paramfile << "maxlen " << MAXLEN << std::endl; 
	paramfile << "par_length_dependence " << par_length_dependence << std::endl; 
	paramfile << "par_death " <<  par_death  << std::endl;
	paramfile << "par_substitution " << par_substitution << std::endl;
	paramfile << "par_insertion " << par_insertion << std::endl;
	paramfile << "par_deletion " << par_deletion << std::endl;
	paramfile << "par_k_noasso " << par_k_noasso << std::endl;
	//paramfile << "par_decay_rate " << par_decay_rate << std::endl;
	paramfile << "par_diffusion_rate " << par_diffusion_rate << std::endl;
	paramfile << "par_maxtime " << par_maxtime << std::endl;
	paramfile << "par_ncol " << par_ncol << std::endl;
	paramfile << "par_nrow " << par_nrow << std::endl;
	paramfile << "par_movie_interval " << par_movie_interval << std::endl;
	paramfile << "par_output_interval " << par_output_interval << std::endl;
	paramfile << "par_pic_folder " << par_pic_folder << std::endl;
	paramfile << "par_ID " << par_ID << std::endl;
	
	paramfile << "par_init_no " << par_init_no << std::endl;
	paramfile << "par_init_seqs ";
	for(i=0; i < par_init_no; i++) paramfile << par_init_seqs[i] << " ";
	paramfile << std::endl;
	paramfile << "par_init_props ";
	for(i=0; i < par_init_no; i++) paramfile << par_init_props[i] << " ";
	paramfile << std::endl;
	
	paramfile.close();
	return(0);
}


//parameter modifying function
int Args(int argc, char **argv)
{
  int i, s;
  char option;
  //  char temp[BUFSIZE];
  /* parameters may be overridden here. */
  for (i = 1; i < argc; i++) {
	if(argv[i][0] == '-'){
		option = argv[i][1];
		switch(option){
			// l - par_length_dependence
			case 'l':
				if (++i == argc) return 1;
				par_length_dependence = atof(argv[i]);
				continue;
			// k - par_death
			case 'k':
				if (++i == argc) return 1;
				par_death = atof(argv[i]);
				if(par_death < 0) {
					std::cerr << "ERROR at reading argoments: option " << option << ": par_death cant be negative!" << std::endl;
					return(-1);
				}
				continue;
			
			case 's':
				if (++i == argc) return 1;
				par_substitution = atof(argv[i]);
				if(par_substitution < 0) {
					std::cerr << "ERROR at reading argoments: option " << option << ": par_substitution cant be negative!" << std::endl;
					return(-1);
				}
				continue;
			
			case 'i':
				if (++i == argc) return 1;
				par_insertion = atof(argv[i]);
				if(par_insertion < 0) {
					std::cerr << "ERROR at reading argoments: option " << option << ": par_insertion cant be negative!" << std::endl;
					return(-1);
				}
				continue;
			
			case 'd':
				if (++i == argc) return 1;
				par_deletion = atof(argv[i]);
				if(par_deletion < 0) {
					std::cerr << "ERROR at reading argoments: option " << option << ": par_deletion cant be negative!" << std::endl;
					return(-1);
				}
				continue;
			
			case 'e':
				if (++i == argc) return 1;
				par_k_noasso = atof(argv[i]);
				if(par_k_noasso < 0) {
					std::cerr << "ERROR at reading argoments: option " << option << ": par_k_noasso cant be negative!" << std::endl;
					return(-1);
				}
				continue;
			
// 			case 'x':
// 				if (++i == argc) return 1;
// 				par_decay_rate = atof(argv[i]);
// 				if(par_decay_rate < 0) {
// 					std::cerr << "ERROR at reading argoments: option " << option << ": par_decay_rate cant be negative!" << std::endl;
// 					return(-1);
// 				}
// 				continue;
			
			case 'D':
				if (++i == argc) return 1;
				par_diffusion_rate = atof(argv[i]);
				if(par_diffusion_rate < 0) {
					std::cerr << "ERROR at reading argoments: option " << option << ": par_diffusion_rate cant be negative!" << std::endl;
					return(-1);
				}
				continue;
			
			case 'T':
				if (++i == argc) return 1;
				par_maxtime = atoi(argv[i]);
				continue;
			
			case 'c':
				if (++i == argc) return 1;
				par_ncol = atoi(argv[i]);
				if(par_ncol <= 0) {
					std::cerr << "ERROR at reading argoments: option " << option << ": par_ncol cant be negative!" << std::endl;
					return(-1);
				}
				continue;
			
			case 'r':
				if (++i == argc) return 1;
				par_nrow = atoi(argv[i]);
				if(par_nrow <= 0) {
					std::cerr << "ERROR at reading argoments: option " << option << ": par_nrow cant be negative!" << std::endl;
					return(-1);
				}
				continue;
			
			case 'm':
				if (++i == argc) return 1;
				par_movie_interval = atoi(argv[i]);
				if(par_movie_interval < 0) {
					std::cerr << "ERROR at reading argoments: option " << option << ": cant be negative!" << std::endl;
					return(-1);
				}
				continue;
			
			case 'o':
				if (++i == argc) return 1;
				par_output_interval = atoi(argv[i]);
				if(par_output_interval < 0) {
					std::cerr << "ERROR at reading argoments: option " << option << ": cant be negative!" << std::endl;
					return(-1);
				}
				continue;
			
			case 'I':
				if (++i == argc) return 1;
				if ( strlen(argv[i]) > 1 ) 
					strcpy(par_ID, argv[i]);
				else {
					std::cerr << "ERROR at reading argoments: option " << option << ": ID should be more than 1 char long!" << std::endl;
					return -1;
				}
				continue;
			
			case 'F':
				if (++i == argc) return 1;
				if ( strlen(argv[i]) > 1 ) 
					strcpy(par_pic_folder, argv[i]);
				else {
					std::cerr << "ERROR at reading argoments: option " << option << ": should be more than 1 char long!" << std::endl;
					return -1;
				}
				continue;
				
			case 'S':
				if (++i == argc) return 1;
				par_init_no = atoi(argv[i]);
				if(par_init_no > 0 && i + 2*par_init_no < argc) {
					par_init_props = new double [par_init_no];
					par_init_seqs = new char* [par_init_no];
					for(s = 0; s < par_init_no; s++){
						par_init_props[s] = atof(argv[++i]);
						if(par_init_props[s] < 0.0) {
							std::cerr << "ERROR at reading argoments: trouble with proportions of initial sequences: it has to be non negative!" << std::endl;
							return -2;
						}
						par_init_seqs[s] = new char [MAXLEN+1];
						if(strlen(argv[++i]) > MAXLEN) {
							std::cerr << "ERROR at reading argoments: trouble initial sequences: has to be shorter than " << MAXLEN << std::endl;
						}
						strcpy(par_init_seqs[s], argv[i]);
					}
				}
				else {
					std::cerr << "ERROR at reading argoments: trouble with initial sequences!" << std::endl;
					return(-2);
				}
				continue;
				
			default:
				std::cerr << "ERROR at reading argoments: not valid argoment!" << std::endl;
				return -1;
		}
	}
  }
  return 0;
}
