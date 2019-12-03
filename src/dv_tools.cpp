#include "dv_tools.h"

#include <iostream>

namespace dvtools {
	
int brokenStickVals(double *values, int noChoices, double sum, double random) {
	/* torott palca egyszeru bemeneti ertekekkel
	 * visszateresi ertekek:
	 * 	0-(noChoices-1): melyiket valassztottuk
	 * 	-1: egyik sem (hiba)
	 * values: a pointer a kumulalt ertekekhez
	 * noChoices: hany tagja van a values-nak
	 * random: a random szam
	 * 
	 * choice: szamlalo, vegigmegy a values-on
	 * cumulate: ebbe kumulaljuk az ertkekeket
	 */
	
	int choice=0;
	double cumulate=0.0;
	
 	if(sum < 0) {
 		for(sum = choice = 0; choice < noChoices; choice++) {
 			sum += *(values + choice);
//			std::cout << "calculating sum: + " << *(values + choice) << " = " << sum << std::endl;			
 		}
 	}
	
	if(sum <= 0) return -1;
	
	for (choice = 0; choice < noChoices; choice++) {		
		if (random < ((cumulate += *(values + choice)) / sum) ) {
//			std::cout << cumulate << "/" << sum << ">" << random << ", so choice is: " << choice << std::endl;
			return(choice);
		}
//		else std::cout << cumulate << "/" << sum << "<=" << random << ", so choice is NOT" << choice << std::endl;
	}
//	std::cerr << "brokenStickVals: some error, maybe sum (" << sum << ") is not valid" << std::endl;
//	for(choice=0; choice < noChoices; choice++) std::cerr << choice << ". choice: " << values[choice] << std::endl;

	return(-2);
}
 

 
}
