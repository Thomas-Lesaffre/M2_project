#ifndef DEPRESSION_H
#define DEPRESSION_H

#include <vector>
#include <iostream>
// #include "MersenneTwister.h"
#include "mt.h"

using namespace std;

// global variables

#define fichierLecture "parameters.txt"     // names of input
#define fichierEcriture "results.txt"		// and output files

// definition of structure "chr" representing a chromosome:
// "M" is the Modifier locus
// "sel" is a vector containing the positions of deleterious alleles along the chromosome

struct chr
{
      double S; // Selected loci
      double E;
      double A; // Modifier locus

      chr(){}; // Default constructor
      ~chr(){}; // Default destructor
};


struct result
{
       double pouet;
};


// Prototypes of functions

void ouvrirFichierE();

void ouvrirFichierS();

bool lireFichier(double &Nr, double &eir, double &air, double &Sr, double &djr, double &dar, double &drr, double &xr, double &Umr, double &rr, double &mutr, double &ngr, double &stepr, double &EvolEr, double &EvolAr, double &itr);

result recursion(double Nv, double eiv, double aiv, double Sv, double djv, double dav, double drv, double xv, double Umv, double rv, double mutv, double ngv, double stepv, double EvolEv, double EvolAv, double itv);

double gammln(const double xx);
double poisdev(const double xm);
double binldev(const double pp, const int n);
double fitness(chr &c1, chr &c2, double wHe, double wHo);
void rec(chr &res, chr &c1, chr &c2, double sz);
void cntl_c_handler(int bidon);


#endif
