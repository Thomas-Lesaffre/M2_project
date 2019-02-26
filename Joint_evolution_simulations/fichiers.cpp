// Functions to open input and output files,
// read parameter values and write them in output file

#include "depression.h"
#include <iostream>
#include <fstream>

using namespace std;

extern FILE * fichierE;
extern FILE * fichierS;


//Opens input file:

void ouvrirFichierE()
{
	fichierE = fopen(fichierLecture,"r");
	if (!fichierE)
		cout << "The file " << fichierLecture
			 << " doesn't exist!\n\n";
}



//Opens output file:

void ouvrirFichierS()
{
	fichierS = fopen(fichierEcriture,"a");
	if (!fichierS)
		cout << "Impossible to open " << fichierEcriture << " !\n\n";
}



// Reads parameter values from input file,
// returns 1 if it reaches the end of input file, else returns 0:

bool lireFichier(double &Nr, double &eir, double &air, double &Sr, double &djr, double &dar, double &drr, double &xr, double &Umr, double &rr, double &mutr, double &ngr, double &stepr, double &EvolEr, double &EvolAr, double &itr)
{
	int z;
	bool term;
	do {z = fgetc(fichierE);} while (!((z == '*') || (z == EOF)));
		// Lines with parameter sets must begin with *
	if (z == EOF)
	{
		cout << "\nEnd of input file\n";
		term = true;
	}
	else
	{
		fscanf(fichierE," %lf",&Nr);
		fscanf(fichierE," %lf",&eir);
		fscanf(fichierE," %lf",&air);
		fscanf(fichierE," %lf",&Sr);
		fscanf(fichierE," %lf",&djr);
		fscanf(fichierE," %lf",&dar);
		fscanf(fichierE," %lf",&drr);
		fscanf(fichierE," %lf",&xr);
		fscanf(fichierE," %lf",&Umr);
		fscanf(fichierE," %lf",&rr);
		fscanf(fichierE," %lf",&mutr);
		fscanf(fichierE," %lf",&ngr);		
		fscanf(fichierE," %lf",&stepr);
		fscanf(fichierE," %lf",&EvolEr);
		fscanf(fichierE," %lf",&EvolAr);
		fscanf(fichierE," %lf",&itr);
							
        term = false;
	}
	
	return term;
}

