#include "depression.h"
#include "mt.h"
#include <iostream>
#include <cmath>
#include <cstdlib>


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

	int i, iteration;
	double N, ei, ai, S, dj, da, dr, x, Um, r, mut, step, it, ng, eE, eA;

	result Rslt;

	// Opens input and output files:

	bool end;
	ouvrirFichierE();
	ouvrirFichierS();
	end = false;

	int no = 1;

	do
	{
        //reads parameter values from input file:

		end = lireFichier(N, ei, ai, S, dj, da, dr, x, Um, r, mut, ng, step, eE, eA, it);
					   	 // end = true if end of input file
        if(!end) cout << "\n___________________________________\n" << "\n";

        if(!end)
            for (i = 0; i < it; i++)
            {

            iteration = i;

            cout << "\n" << "Iteration Number : "<< i << "\n";

			// Simulation:

 			Rslt = recursion(N, ei, ai, S, dj, da, dr, x, Um, r, mut, ng, step, eE, eA, iteration);

            }

            no++;

	} while (!end);

	// Closes files:

    fclose(fichierE);
	fclose(fichierS);

	return 0 ;
}


