#include "depression.h"
#include "mt.h"
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <fstream>

using namespace std;

// Random number generator:

MTRand rnd;

// Pointers on input and output files:

FILE * fichierE;
FILE * fichierS;

int main()
{

  cout << "Program initialization" << "\n";

	// Parameters:

	double Rslt, Estart, pM, Pmax, x, alpha, deltaA, deltaJ, hm, epsi;

	// Opens input and output files:

	bool end;
	ouvrirFichierE();
	ouvrirFichierS();
	end = false;

	int no = 1;
	
	do{
	
        //reads parameter values from input file:

		end = lireFichier(pM, Estart, x, Pmax, alpha, deltaA, deltaJ, hm, epsi);
					   	 // end = true if end of input file
					   	 
        if(!end) cout << "\n___________________________________\n" << "\n";

        if(!end)

 			Rslt = recursion(pM, Estart, x, Pmax, alpha, deltaA, deltaJ, hm, epsi);

            no++;

	}while(!end);

	// Closes files:

    fclose(fichierE);
	fclose(fichierS);

    /*
    system("mkdir results");
    system("mv result*.txt results/");
    system("tar -czvf results.tar.gz results");
    system("rm -r results");
    */

	return 0 ;
}
