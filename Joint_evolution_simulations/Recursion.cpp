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

In this program, the joint evolution of modifiers of reproductive effort and selfing rate is investigated. 

Each individual is depicted by two chromosomes. These chromosomes are structures, that each contain three doubles, 'S', 'E' and 'A':
	- 'S' is set to 1 or 0 depending on whether the considered individual was selfed (1) or not (0).
	- 'E' depicts the reproductive effort modifier.
	- 'A' depicts the selfing rate modifier.

Modifiers are allowed to mutate freely in the [0,1] interval. Since one value corresponds to one allele, there can potentially be infinitely many alleles segregating at these modifiers.


*/

/* 

Parameters :
    - Nv: Population size
    - eiv, aiv: initial reproductive effort and selfing rate, respectively
    - Sv: Maximal survival probability
    - djv, dav, drv: inbreeding depression affecting juvenile survival, adult survival and fecundity, respectively.
    - xv: parameter controlling the shape of the survival vs. reproduction trade-off.
    - Umv: mutation rate at the modifiers.
    - rv: recombination rate between the modifiers.
    - mutv: parameter controlling the magnitude of muations
    - ngv: Number of generations for which the simulation is ran.
    - stepv: frequency of measurements.
    - EvolEv, EvolAv: parameters controlling whether or not the evolution of reproductive effort and selfing should be triggered.
    - itv: parameter giving iteration number (this is only used for naming the output file here).
    
*/

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

	// Initialization: 
	
	for(i=0; i<Nv; i++)
	{
		nb = 2*i;
		
		// Individuals all start with the same reproductive effort and selfing rate, 'eiv' and 'aiv'
		
		pop[nb].E = eiv;
		pop[nb+1].E = eiv;
		
		pop[nb].A = aiv;
		pop[nb+1].A = aiv;
		
		// According to the initial selfing rate, individuals are assigned a status at random: S=1 if they are selfed, and S=0 otherwise. 
		
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
	
	for(g=0; g<ngv; g++) // Run recursion for ngv timesteps.
	{	
		if (cntl_c_bool) // Control for stopping
		{
			stp = 1;
			break;
		}
		
		// Relative contribution vector: a cumulated sum vector of reproductive effort is created so that individuals can be drawn from it for reproduction. 
		// The higher the reproductive effort, the higher the probability that a real number taken between 0 and sumE will fall in the individuals' range

		vE.clear();
		sumE = 0;

		for(i=0; i<Nv; i++)
		{
			nb=2*i;
			
			// Parental chromosomes are copied into the 'par' vector, so that they can be used to produce new offspring when one individual dies.
			
			par[nb] = pop[nb];
			par[nb+1] = pop[nb+1];
			
			sumE += ((par[nb].E + par[nb+1].E)/2)*(1 - drv*((par[nb].S + par[nb+1].S)/2));
			vE.push_back(sumE);

		}
		
		// Recursion, iterating through the 'pop' vector.
		
		for(i=0; i<Nv; i++)
		{
		
			nb = 2*i;
			
			// Survival probability depends on the maximal survival probability Sv, its status S, and thereby of adult inbreeding depression dav if S=1, 
			// its reproductive effort E, and the shape of the trade-off xv
			
			if(rnd.rand() < Sv*(1 - ((pop[nb].S+pop[nb+1].S)/2)*dav)*(1 - pow((pop[nb].E + pop[nb+1].E)/2, xv))) // If the individual survives
			{
				// Nothing happens, the individual stays in place !
			}
			else
			{	
				// If the individual dies, we have to build a replacement.
				
				do{
					// Select a mother: the first encountered value in vE that is higher than erand stops the 'do' loop and the corresponding individual is chosen as mother.
				
					erand = rnd.rand(sumE);
					k = -1;
					do{
						k++;
					}while(vE[k] < erand);
				
					// Selfing or Outcrossing ?
					
					if(rnd.rand() < (par[2*k].A + par[2*k + 1].A)/2) // If the chosen individual reproduces by self-fertilisation
					{
						// Two chromosomes are sampled with replacement from its two chromosomes, and they recombine with probability rv
						// The alleles are then transmitted to the offspring's chromosomes off1 and off2.
						
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
						
						// S is set to 1 because the offspring was selfed.
						
						off1.S = 1;
						off2.S = 1;
					}
					else // If the individual outcrosses,
					{
						ch1 = 2*k + rnd.randInt(1); // One of the mother's chromosome is picked, and recombines with probability rv...
						
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
						
						// ... and we choose a father the same way we chose a mother.
						
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
						
						// S is set to 0 since the offspring was outcrossed.
						
						off1.S = 0;
						off2.S = 0;															
					}
					
					// Mutation at the modifiers. 
					// With probability Umv, modifiers are assigned a new value, sampled in a uniform distribution in a [x - mutv, x + mutv] interval, where 'x' is the former modifier value.

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
													
					fitj = 1 - ((off1.S + off2.S)/2)*djv; // The offspring fitness is calculated according to its status (selfed or outcrossed)...
					
				}while(fitj < rnd.rand()); // ...and it makes it through recruitment with probability fitj.
				
				// Once recruited, it is incorporated in the population in replacement of the dead adult.
				
				pop[nb] = off1;
				pop[nb+1] = off2;
			}
		}
		
		// Measurements every stepv
		
		mA=0;
		mE=0;
		mQ=0;
		
		for(i=0; i<Nv; i++) // At each timestep, the average selfing rate and reproductive effort, and the proportion of selfed individuals are computed and stored in vectors.
		{
			nb = 2*i;
			
			mA += (pop[nb].A + pop[nb+1].A)/2;
			mE += (pop[nb].E + pop[nb+1].E)/2;
			mQ += (pop[nb].S + pop[nb+1].S)/2;
		}
		
		mA /= Nv;
		mE /= Nv;
		mQ /= Nv;
		
		// cout << mA << " " << mE << " " << mQ << endl;
		
		vmA.push_back(mA);
		vmE.push_back(mE);
		vmQ.push_back(mQ);

		if(g % (int) stepv == 0) // Every stepv generation, the average of the vectors is printed into the output file, and the vectors are cleared.
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

