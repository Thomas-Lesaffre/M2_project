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

// Parameters :
    // Nv : Size of pop.
    // NbGenPrelimv : Nb of generations before the M-locus starts mutating.
    // NbGen : Total nb of generations in the simulation after M starts mutating.
    // Estartv : Initial reproductive investment.
    // Pmaxv : Extrinsic mortality --> The higher Pmaxv is, the lower ext. mortality is.
    // xv : Shape of the survival gain curve.
    // sv : fitness effect of a deleterious mutation.
    // hv : dominance level of a deleterious mutation.
    // alphav : selfing rate.
    // Uv : genomic mutation rate for deleterious mutations.
    // Umv : Mutation rate at the Modifier Locus.
    // Lv : Genome map length in cM.
    // stepv : Nb of gen. between two measurements of ID and E.
    // samplev : Number of individuals sampled to estimate ID in the population.

// Recursion :

result recursion(double Nv, double eiv, double aiv, double Sv, double djv, double dav, double drv, double xv, double Umv, double rv, double mutv, double ngv, double stepv, double EvolEv, double EvolAv, double itv)
{

    // Naming the output file :

		char nomFichier[256];

		stringstream nomF;

		nomF << "coevo" << "_N" << Nv << "_S" << Sv << "_dj" << djv << "_da" << dav << "_dr" << drv << "_ai" << aiv  << "_ei" << eiv << "_r" << rv << "_it" << itv << ".txt";

		nomF >> nomFichier;

		ofstream fout(nomFichier);
		
		fout << "m_a" << " " << "m_e" << " " << "m_Q" << endl;
		
   	// For stopping prog. with Ctrl-C :
   	
		int stp;
		
	   	cntl_c_bool = false;
		signal(SIGINT, cntl_c_handler);
	
	// Declaring stuff
	
		int twoN = 2*Nv; 
	
		chr * pop = new chr[twoN];
		chr * par = new chr[twoN];
		
		chr off1, off2;
		
		result Res;
		
		int g, i, k, nb, ch1, ch2;
		double sumE, erand, mA, mE, mQ, smA, smE, smQ, fitj;
		
		vector<double> vE, vmA, vmE, vmQ;
	
	// Fitness effects 
		
		
	// Initialization
	
	for(i=0; i<Nv; i++)
	{
		nb = 2*i;
		
		pop[nb].E = eiv;
		pop[nb+1].E = eiv;
		
		pop[nb].A = aiv;
		pop[nb+1].A = aiv;
		
		if(rnd.rand() < aiv)
		{
			pop[nb].S = 1;
			pop[nb+1].S = 1;
		}
		else
		{
			pop[nb].S = 0;
			pop[nb+1].S = 0;		
		}
	}

	// Recursion
	
	for(g=0; g<ngv; g++)
	{	
		if (cntl_c_bool) // Control for stopping
		{
			stp = 1;
			break;
		}
		
		// Relative contribution vector

		vE.clear();
		sumE = 0;

		for(i=0; i<Nv; i++)
		{
			nb=2*i;
			
			par[nb] = pop[nb];
			par[nb+1] = pop[nb+1];
			
			sumE += ((par[nb].E + par[nb+1].E)/2)*(1 - drv*((par[nb].S + par[nb+1].S)/2));
			vE.push_back(sumE);

		}
		
		// Recursion
		
		for(i=0; i<Nv; i++)
		{
		
			nb = 2*i;
			
			if(rnd.rand() < Sv*(1 - ((pop[nb].S+pop[nb+1].S)/2)*dav)*(1 - pow((pop[nb].E + pop[nb+1].E)/2, xv))) // If the individual survives
			{
				// Nothing happens !
			}
			else
			{	
				do{
					// Select a mother
				
					erand = rnd.rand(sumE);
					k = -1;
					do{
						k++;
					}while(vE[k] < erand);
				
					// Selfing or Outcrossing ?
					
					if(rnd.rand() < (par[2*k].A + par[2*k + 1].A)/2)
					{
						ch1 = 2*k + rnd.randInt(1);
						
						if(ch1 % 2 == 0)
						{					
							ch2 = ch1 + rnd.randInt(1);
							
							if(rnd.rand() < rv)
							{
								off1.E = par[ch1].E;
								off1.A = par[ch1+1].A;
							}
							else
							{
								off1.E = par[ch1].E;
								off1.A = par[ch1].A;
							}
							
							if(ch2 % 2 == 0)
							{
								if(rnd.rand() < rv)
								{
									off2.E = par[ch2].E;
									off2.A = par[ch2+1].A;
								}
								else
								{
									off2.E = par[ch2].E;
									off2.A = par[ch2].A;
								}								
							}
							else
							{
								if(rnd.rand() < rv)
								{
									off2.E = par[ch2].E;
									off2.A = par[ch2-1].A;
								}
								else
								{
									off2.E = par[ch2].E;
									off2.A = par[ch2].A;
								}							
							}
						}
						else
						{						
							ch2 = ch1 - rnd.randInt(1);
							
							if(rnd.rand() < rv)
							{
								off1.E = par[ch1].E;
								off1.A = par[ch1-1].A;
							}
							else
							{
								off1.E = par[ch1].E;
								off1.A = par[ch1].A;
							}

							if(ch2 % 2 == 0)
							{
								if(rnd.rand() < rv)
								{
									off2.E = par[ch2].E;
									off2.A = par[ch2+1].A;
								}
								else
								{
									off2.E = par[ch2].E;
									off2.A = par[ch2].A;
								}								
							}
							else
							{
								if(rnd.rand() < rv)
								{
									off2.E = par[ch2].E;
									off2.A = par[ch2-1].A;
								}
								else
								{
									off2.E = par[ch2].E;
									off2.A = par[ch2].A;
								}							
							}						
						}
						
						off1.S = 1;
						off2.S = 1;
					}
					else
					{
						ch1 = 2*k + rnd.randInt(1);
						
						if(ch1 % 2 == 0)
						{
							if(rnd.rand() < rv)
							{
								off1.E = par[ch1].E;
								off1.A = par[ch1+1].A;
							}
							else
							{
								off1.E = par[ch1].E;
								off1.A = par[ch1].A;
							}
						}
						else
						{
							if(rnd.rand() < rv)
							{
								off1.E = par[ch1].E;
								off1.A = par[ch1-1].A;
							}
							else
							{
								off1.E = par[ch1].E;
								off1.A = par[ch1].A;
							}
						}
						
						// Choose a father
						
						erand = rnd.rand(sumE);
						k = -1;
						do{
							k++;
						}while(vE[k] < erand);
						
						
						ch2 = 2*k + rnd.randInt(1);
						
						if(ch2 % 2 == 0)
						{
							if(rnd.rand() < rv)
							{
								off2.E = par[ch2].E;
								off2.A = par[ch2+1].A;
							}
							else
							{
								off2.E = par[ch2].E;
								off2.A = par[ch2].A;
							}
						}
						else
						{
							if(rnd.rand() < rv)
							{
								off2.E = par[ch2].E;
								off2.A = par[ch2-1].A;
							}
							else
							{
								off2.E = par[ch2].E;
								off2.A = par[ch2].A;
							}
						}
						
						off1.S = 0;
						off2.S = 0;															
					}
					
					// Mutation

						// Rep. Effort
						
						if(rnd.rand() < Umv*EvolEv)
						{
							if(rnd.rand() < 0.5)
							{
								off1.E -= mutv*rnd.rand();
							}
							else
							{
								off1.E += mutv*rnd.rand();					
							}
							
						}
						
						if(rnd.rand() < Umv*EvolEv)
						{
							if(rnd.rand() < 0.5)
							{
								off2.E -= mutv*rnd.rand();
							}
							else
							{
								off2.E += mutv*rnd.rand();						
							}

						}
						
						// Selfing rate
						
						if(rnd.rand() < Umv*EvolAv)
						{
							if(rnd.rand() < 0.5)
							{
								off1.A -= mutv*rnd.rand();
							}
							else
							{
								off1.A += mutv*rnd.rand();						
							}
						}
						
						if(rnd.rand() < Umv*EvolAv)
						{						
							if(rnd.rand() < 0.5)
							{
								off2.A -= mutv*rnd.rand();
							}
							else
							{
								off2.A += mutv*rnd.rand();						
							}
						}
						
						// Keep modifiers within definition domain
						
						if(off1.A < 0)
						{
							off1.A = 0;
						}
						else if(off1.A > 1)
						{
							off1.A = 1;
						}

						if(off1.E < 0)
						{
							off1.E = 0;
						}
						else if(off1.E > 1)
						{
							off1.E = 1;
						}

						if(off2.A < 0)
						{
							off2.A = 0;
						}
						else if(off2.A > 1)
						{
							off2.A = 1;
						}

						if(off2.E < 0)
						{
							off2.E = 0;
						}
						else if(off2.E > 1)
						{
							off2.E = 1;
						}			
													
					fitj = 1 - ((off1.S + off2.S)/2)*djv;
					
				}while(fitj < rnd.rand());
				
				pop[nb] = off1;
				pop[nb+1] = off2;
			}
		}
		
		// Measurements
		
		mA=0;
		mE=0;
		mQ=0;
		
		for(i=0; i<Nv; i++)
		{
			nb = 2*i;
			
			mA += (pop[nb].A + pop[nb+1].A)/2;
			mE += (pop[nb].E + pop[nb+1].E)/2;
			mQ += (pop[nb].S + pop[nb+1].S)/2;
		}
		
		mA /= Nv;
		mE /= Nv;
		mQ /= Nv;
		
		cout << mA << " " << mE << " " << mQ << endl;
		
		vmA.push_back(mA);
		vmE.push_back(mE);
		vmQ.push_back(mQ);

		if(g % (int) stepv == 0)
		{			
			smE = 0;
			smA = 0;
			smQ = 0;
			
			for(i=0; i<vmE.size(); i++)
			{
				smE += vmE[i];
				smA += vmA[i];
				smQ += vmQ[i];
			}
			
			smE /= vmE.size();
			smA /= vmA.size();
			smQ /= vmQ.size();
			
			vmE.clear();
			vmA.clear();
			vmQ.clear();
			
			fout << smE << " " << smA << " " << smQ << endl;
		}
		
	}

    delete[] pop;
	delete[] par;
		
    return(Res);
}

