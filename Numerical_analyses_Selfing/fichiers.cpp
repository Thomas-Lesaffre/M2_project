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

bool lireFichier(double &pMr, double &Astartr, double &dAr, double &delta_jr, double &delta_ar, double &Pr)
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
		fscanf(fichierE," %lf",&pMr);
		fscanf(fichierE," %lf",&Astartr);
		fscanf(fichierE," %lf",&dAr);
		fscanf(fichierE," %lf",&delta_jr);
		fscanf(fichierE," %lf",&delta_ar);
		fscanf(fichierE," %lf",&Pr);

        term = false;
	}
	return term;
}

bool isFileExist(const std::string& fileName)
{
    return std::ifstream(fileName.c_str());
}

