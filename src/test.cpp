#include <iostream>
#include <string.h>

extern "C" {
#include <ViennaRNA/fold.h>
#include <ViennaRNA/utils/basic.h>
}

using namespace std;


int main(int argc, char *argv[]) {
    /* The RNA sequence */
    char seq[] = "GAGUAGUGGAACCAGGCUAUGUUUGUGACUCGCAGACUAACA";

    // allocate memory for MFE structure (length + 1)
    char *structure = (char *) vrna_alloc(sizeof(char) * (strlen(seq) + 1));
    
    // predict Minmum Free Energy and corresponding secondary structure
    float mfe = vrna_fold(seq, structure);
    
    // print sequence, structure and MFE
    cout << seq << endl;
    cout << structure << mfe << endl;
    printf("%s\n%s [ %6.2f ]\n", seq, structure, mfe);
    //cleanup memory
    free(structure);

   
   
   return 0;
} 

