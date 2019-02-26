#ifndef DEPRESSION_H
#define DEPRESSION_H

#include <vector>
#include <iostream>
// #include "MersenneTwister.h"
// #include "mt.h"

using namespace std;

// global variables

#define fichierLecture "parameters.txt"     // names of input
#define fichierEcriture "results.txt"		// and output files


// Prototypes of functions

void ouvrirFichierE();

void ouvrirFichierS();

bool lireFichier(double &pMr, double &Estartr, double &xr, double &Pmaxr, double &alphar, double &deltaAr, double &deltaJr, double &hmr, double &epsir);

double recursion(double pMv, double Estartv, double xv, double Pmaxv, double alphav, double deltaAv, double deltaJv, double hmv, double epsiv);

bool isFileExist(const std::string& fileName);


double gammln(const double xx);
double poisdev(const double xm);
double binldev(const double pp, const int n);

void cntl_c_handler(int bidon);


#endif
