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

double recursion(double pMv, double Astartv, double dAv, double delta_jv, double delta_av, double Pv)
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

	int  i;
	double Res, X, Y, Z, A, A_m, A_ess, rXs, rYs, rZs, rXo, rYo, rZo, gpM, gqM, tot_g, pMr, qMr, mean_fitr, rsXs, rsYs, rsZs, rsXo, rsYo, rsZo, S, sXs, sXo, sYs, sYo, sZs, sZo, p, diff;
	
	vector<double> Gen(6);
	
	A = Astartv;
	
	do{
	
	// Initialiation
		
		Gen[0] = (1-A)*(pow(pMv,2) + (A/(2 - A))*pMv*(1-pMv));
		Gen[1] = A*(pow(pMv,2) + (A/(2 - A))*pMv*(1-pMv));
		Gen[2] = (1-A)*(2*pMv*(1-pMv) - 2*(A/(2 - A))*pMv*(1 - pMv));
		Gen[3] = A*(2*pMv*(1-pMv) - 2*(A/(2 - A))*pMv*(1 - pMv));
		Gen[4] = (1-A)*(pow(1 - pMv,2) + (A/(2 - A))*pMv*(1-pMv));
		Gen[5] = A*(pow(1 - pMv,2) + (A/(2 - A))*pMv*(1-pMv));
	
	// Recursion
	
		for(i = 0; i < 1000000; i++)
		{
	
			X = Gen[0] + Gen[1];
			Y = Gen[2] + Gen[3];
			Z = Gen[4] + Gen[5];
		
		// Reproduction
	
			rXs = A*X + (A + dAv/2)*(Y/4);
			rYs = (A + dAv/2)*(Y/2);
			rZs = (A + dAv/2)*(Y/4) + (A + dAv)*Z;
		
			A_m = rXs + rYs + rZs;
			
			gpM = X*(1-A) + (Y/2)*(1 - (A + dAv/2));
			gqM = Z*(1 - (A + dAv)) + (Y/2)*(1 - (A + dAv/2));
			
			tot_g = gpM + gqM;
			
			pMr = gpM/tot_g;
			qMr = gqM/tot_g;
		
			rXo = (1 - A_m)*pMr*(X + Y/2);
			rYo = (1 - A_m)*(pMr*(Z + Y/2) + qMr*(X + Y/2));
			rZo = (1 - A_m)*qMr*(Z + Y/2);
	
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
		
		p = Gen[0] + Gen[1] + 0.5*(Gen[3] + Gen[4]);
		
		diff = p - pMv;
		
		if(diff < 0)
		{
			A += dAv;
			A_ess = A;
			
			cout << "up : " << A << endl;
		}
		else
		{
			A_ess = A;
		}
	
	}while((A <= 1) && (A >= -dAv) && (diff <= 0));
	
	Res = A;
	
	fout << Pv << " " << delta_jv << " " << delta_av << " " << A << endl;
	
	cout << "P = " << Pv << " | delta_j = " << delta_jv << " | delta_a = " << delta_av << " | alpha = " << A << endl;
	
	fout.close();
    
    return(Res);
}

