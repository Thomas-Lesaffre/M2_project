#include "depression.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <vector>
#include <cmath>
#include <csignal>
#include <algorithm>
#include "mt.h"
#include <cstring>

using namespace std;

extern MTRand rnd;
extern FILE * fichierE;
extern FILE * fichierS;

// For stopping the program with Ctrl-C:

bool cntl_c_bool = false;
void cntl_c_handler (int bidon)
{
	cntl_c_bool = true;
}

/* 

The purpose of  this program is to perform numerical analyses of the exact recursions on genotypic frequencies for the case of 
a modifier locus affecting reproductive effort with a fixed selfing rate. Here, in the 'upwards' direction, that is,
starting from a low resident reproductive effort and introducing a mutant increasing reproductive effort,
then gradually increasing the reproductive effort coded by the resident, as the mutant invades, until it no longer does.

*/


// Parameters :
	// pMv : initial resident allele frequency
	// Estartv : Intial reproductive effort
	// xv : parameter controlling the survival vs. reproduction trade-off shape
	// Pmaxv : Extrinsic mortality parameter
	// alphav : Selfing rate
	// deltaAv : Adult inbreeding depression
	// deltaJv : Juvenile inbreeding depression
	// hmv : dominance coefficient at the modifier
	// epsiv : effect of the modifier
	
// Recursion :

double recursion(double pMv, double Estartv, double xv, double Pmaxv, double alphav, double deltaAv, double deltaJv, double hmv, double epsiv)
{
	
	// Defining the 'output file' :
	
	// Naming the output file :

	char nomFichier[256];

	stringstream nomF;

	nomF << "Num_Analysis_output.txt";

	nomF >> nomFichier;
	
	ofstream fout;

	fout.open(nomFichier, std::ios::app);

    // Declaring stuff :

    double A, p, Xi, Yi, Zi, sXo, sXs, sYo, sYs, sZo, sZs, sTOT, gamX, gamY, gamZ, X, Y, Z, gamTOT, gamXr, gamYr, gamZr, gamp, gamq, dTOT, dXs, dYs, dZs, dXo, dYo, dZo, dXsr, dYsr, dZsr, dXor, dYor, dZor, diff, Res, ESS, E;
    int i;
    
    vector<double> Gen(6);
    
    // Initialization :

	E = Estartv; // Initial reproductive effort.
	ESS = Estartv; // ESS initialized at E.
	
	do{
	
			Xi = pow(pMv,2) + pMv*(1-pMv)*(alphav/(2-alphav)); // Initial {M,M} frequency
			Yi = 2*pMv*(1-pMv)*(1-(alphav/(2-alphav))); // Initial {M,m} frequency
			Zi = pow(1-pMv,2) + pMv*(1-pMv)*(alphav/(2-alphav)); // Initial {m,m} frequency
	
			Gen[0] = (1 - alphav)*pow(Xi + 0.5*Yi,2);  // Outcrossed {M,M}
			Gen[1] = alphav*pow(Xi + 0.5*Yi,2);	// Selfed {M,M}
	
			Gen[2] = (1 - alphav)*2*(Xi + 0.5*Yi)*(Zi + 0.5*Yi); // Outcrossed {M,m}
			Gen[3] = alphav*2*(Xi + 0.5*Yi)*(Zi + 0.5*Yi); // Selfed {M,m}
	
			Gen[4] = (1 - alphav)*pow(Zi + 0.5*Yi,2); // Outcrossed {m,m}
			Gen[5] = alphav*pow(Zi + 0.5*Yi,2); // Selfed {m,m}
	
		for(i=0; i< 1000000; i++){ // Run recusrions for 10^6 iterations
		
			// Genotypic frequencies
			
			X = Gen[0] + Gen[1];
			Y = Gen[2] + Gen[3];
			Z = Gen[4] + Gen[5];
			
			// Survival according to reproductive effort, maximal survival probability Pmaxv, status (selfed or outcrossed) and thereby adult inbreeding depression deltaAv, and trade-off shape xv.
			
			sXo = Gen[0]*Pmaxv*(1 - pow(E,xv));
			sYo = Gen[2]*Pmaxv*(1 - pow(E + epsiv*hmv,xv));
			sZo = Gen[4]*Pmaxv*(1 - pow(E + epsiv,xv));
	
			sXs = Gen[1]*Pmaxv*(1 - pow(E,xv))*(1-deltaAv);
			sYs = Gen[3]*Pmaxv*(1 - pow(E + epsiv*hmv,xv))*(1-deltaAv);
			sZs = Gen[5]*Pmaxv*(1 - pow(E + epsiv,xv))*(1-deltaAv);
	
			sTOT = sXo + sYo + sZo + sXs + sYs + sZs;
			
			// Reproduction
			
			// Genotypic frequencies in gametes
			
			gamX = X*E;
			gamY = Y*(E + epsiv*hmv);
			gamZ = Z*(E + epsiv);
	
			gamTOT = gamX+gamY+gamZ;
	
			gamXr = gamX/gamTOT;
			gamYr = gamY/gamTOT;
			gamZr = gamZ/gamTOT;
			
			// Allelic frequencies in gametes
	
			gamp = gamXr + 0.5*gamYr;
			gamq = gamZr + 0.5*gamYr;
			
			// Selfed offspring
	
			dXs = alphav*(gamXr + 0.25*gamYr)*(1-deltaJv);
			dYs = alphav*(0.5*gamYr)*(1-deltaJv);
			dZs = alphav*(gamZr + 0.25*gamYr)*(1-deltaJv);
			
			// Outcrossed offspring
	
			dXo = (1-alphav)*pow(gamp,2);
			dYo = (1-alphav)*2*gamp*gamq;
			dZo = (1-alphav)*pow(gamq,2);
			
			dTOT = dXs + dYs + dZs + dXo + dYo + dZo;
			
			// Renormalized frequencies in the offspring
			
			dXsr = dXs/dTOT;
			dYsr = dYs/dTOT;
			dZsr = dZs/dTOT;
			
			dXor = dXo/dTOT;
			dYor = dYo/dTOT;
			dZor = dZo/dTOT;
			
			// Genotypic frequencies in the next timestep.
			
			Gen[0] = sXo + (1-sTOT)*dXor;
			Gen[1] = sXs + (1-sTOT)*dXsr;
	
			Gen[2] = sYo + (1-sTOT)*dYor;
			Gen[3] = sYs + (1-sTOT)*dYsr;
	
			Gen[4] = sZo + (1-sTOT)*dZor;
			Gen[5] = sZs + (1-sTOT)*dZsr;
	
		}
		
		p = Gen[0] + Gen[1] + 0.5*(Gen[2] + Gen[3]); // Resident allele frequency
		
		A = Gen[1] + Gen[3] + Gen[5]; // Proportion of selfed individuals
		
		diff = p - pMv; // Change in frequency of the resident
	
		if(diff <= 0)
		{
			E += epsiv;
			ESS += epsiv;
		}
		else
		{
			ESS = E;
		}
	
	}while((E <= 1) && (E >= 0) && (diff <= 0));
	
	Res = ESS;
	
	fout << alphav << " " << deltaAv << " " << deltaJv << " " << Pmaxv << " " << ESS  << " " << sTOT  << " " << A << endl;
	
	// cout << "ESS = " << ESS << "|" << " alpha = " << alphav << " deltaA = " << deltaAv << " deltaJ = " << deltaJv << " Pmax = " << Pmaxv << " x = " << xv << endl;
	
	fout.close();
    
    return(Res);
}

