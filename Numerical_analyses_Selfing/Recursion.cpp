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
a modifier locus affecting the selfing rate with a fixed reproductive effort. Starting from a low resident selfing rate and introducing a mutant increasing it,
then gradually increasing the selfing rate coded by the resident, as the mutant invades, until it no longer does.

*/


// Parameters:
	// pMv: initial resident allele frequency
	// Astartv: Intial selfing rate
	// dAv: magnitude of mutations
	// delta_jv: juvenile inbreeding depression
	// delta_av: adult inbreeding depression
	// Pv: maximal survival probability.
	
// Recursion:

double recursion(double pMv, double Astartv, double dAv, double delta_jv, double delta_av, double Pv)
{
	
	// Defining the 'output file':
	
	// Naming the output file:

	char nomFichier[256];

	stringstream nomF;

	nomF << "Num_Analysis_output.txt";

	nomF >> nomFichier;
	
	ofstream fout;

	fout.open(nomFichier, std::ios::app);

    // Declaring stuff :

	int  i;
	double Res, X, Y, Z, A, A_m, A_ess, rXs, rYs, rZs, rXo, rYo, rZo, gpM, gqM, tot_g, pMr, qMr, mean_fitr, rsXs, rsYs, rsZs, rsXo, rsYo, rsZo, S, sXs, sXo, sYs, sYo, sZs, sZo, p, diff;
	
	vector<double> Gen(6);
	
	A = Astartv;
	
	do{
	
	// Initialization
		
		Gen[0] = (1-A)*(pow(pMv,2) + (A/(2 - A))*pMv*(1-pMv)); // Initial outcrossed {M,M} frequency
		Gen[1] = A*(pow(pMv,2) + (A/(2 - A))*pMv*(1-pMv));  // Initial selfed {M,M} frequency
		Gen[2] = (1-A)*(2*pMv*(1-pMv) - 2*(A/(2 - A))*pMv*(1 - pMv));  // Initial outcrossed {M,m} frequency
		Gen[3] = A*(2*pMv*(1-pMv) - 2*(A/(2 - A))*pMv*(1 - pMv));  // Initial selfed {M,m} frequency
		Gen[4] = (1-A)*(pow(1 - pMv,2) + (A/(2 - A))*pMv*(1-pMv));  // Initial outcrossed {m,m} frequency
		Gen[5] = A*(pow(1 - pMv,2) + (A/(2 - A))*pMv*(1-pMv));  // Initial selfed {m,m} frequency
	
	// Recursion
	
		for(i = 0; i < 1000000; i++) // Run for 10^6 timesteps
		{
	
			X = Gen[0] + Gen[1]; // Frequency of {M,M}
			Y = Gen[2] + Gen[3]; // Frequency of {M,m}
			Z = Gen[4] + Gen[5]; // Frequency of {m,m}
		
		// Reproduction
	
			//Selfing
			
			rXs = A*X + (A + dAv/2)*(Y/4); // Production of {M,M} by selfing.
			rYs = (A + dAv/2)*(Y/2); // Production of {M,m} by selfing.
			rZs = (A + dAv/2)*(Y/4) + (A + dAv)*Z; // Production of {m,m} by selfing.
		
			A_m = rXs + rYs + rZs;
			
			// Outcrossing
			
			gpM = X*(1-A) + (Y/2)*(1 - (A + dAv/2));
			gqM = Z*(1 - (A + dAv)) + (Y/2)*(1 - (A + dAv/2));
			
			tot_g = gpM + gqM;
			
			pMr = gpM/tot_g;
			qMr = gqM/tot_g;
		
			rXo = (1 - A_m)*pMr*(X + Y/2); // Production of {M,M} by outcrossing.
			rYo = (1 - A_m)*(pMr*(Z + Y/2) + qMr*(X + Y/2)); // Production of {M,m} by outcrossing.
			rZo = (1 - A_m)*qMr*(Z + Y/2); // Production of {m,m} by outcrossing.
			
			// Selection for recruitment
			
			mean_fitr = A_m*(1-delta_jv) + 1 - A_m;
		
			rsXs = rXs*(1 - delta_jv)/mean_fitr;
			rsYs = rYs*(1 - delta_jv)/mean_fitr;
			rsZs = rZs*(1 - delta_jv)/mean_fitr;
		
			rsXo = rXo*(1/mean_fitr);
			rsYo = rYo*(1/mean_fitr);
			rsZo = rZo*(1/mean_fitr);
		
		// Survival
		
			sXs = Gen[1]*Pv*(1 - delta_av);
			sYs = Gen[3]*Pv*(1 - delta_av);
			sZs = Gen[5]*Pv*(1 - delta_av);		
			
			sXo = Gen[0]*Pv;
			sYo = Gen[2]*Pv;
			sZo = Gen[4]*Pv;	
			
			S = sXs + sXo + sYs + sYo + sZs + sZo;
		
		// Next Generation 
		
			Gen[0] = sXo + (1-S)*rsXo;
			Gen[1] = sXs + (1-S)*rsXs;
			Gen[2] = sYo + (1-S)*rsYo;
			Gen[3] = sYs + (1-S)*rsYs;
			Gen[4] = sZo + (1-S)*rsZo;
			Gen[5] = sZs + (1-S)*rsZs;
			
		}
		
		p = Gen[0] + Gen[1] + 0.5*(Gen[3] + Gen[4]); // New frequency of the resident
		
		diff = p - pMv; // change in frequency of the resident
		
		if(diff < 0) // If it was invaded by the mutant...
		{
			A += dAv; // ...Increase A
			A_ess = A;
			
			// cout << "up : " << A << endl;
		}
		else
		{
			A_ess = A; // ...Otherwise, we have found the ESS.
		}
	
	}while((A <= 1) && (A >= -dAv) && (diff <= 0));
	
	Res = A;
	
	fout << Pv << " " << delta_jv << " " << delta_av << " " << A << endl;
	
	// cout << "P = " << Pv << " | delta_j = " << delta_jv << " | delta_a = " << delta_av << " | alpha = " << A << endl;
	
	fout.close();
    
    return(Res);
}

