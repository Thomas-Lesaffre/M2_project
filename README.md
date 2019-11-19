# The joint evolution of lifespan and selfing 
# (DOI: https://doi.org/10.1101/420877)

The preprint version of the manuscript can be found here: https://doi.org/10.1101/420877. The preprint is recommended in _PeerCommunityIn Evol Biol_, recommendation by Thomas Bataillon can be found here: https://doi.org/10.24072/pci.evolbiol.100070.

This article is now published in _Journal of Evolutionary Biology_ (DOI: https://doi.org/10.1111/jeb.13543)

--------------------------------------------------------------------------------------------------

This project contains four programs in total:

- The 'Joint_evolution_simulations' repository contains a program for individual-centered simulations of the joint evolution of reproductive effort and selfing rate.

- The 'Numerical_analyses_Reproductive_Effort' repository contains two subfolders, 'Upwards' and 'Downwards', which contain two separate programs. Both programs perform iterations of the exact recursions on genotypic frequencies for the evolution of reproductive effort with a fixed selfing rate (see comments in the code for details), but with two different starting points. The 'Upwards' program starts from a very low resident reproductive effort, and increases it gradually until the mutant allele, which increases reproductive effort no longer invades the population. The 'Downwards' program starts from a high reproductive effort, and decreases it until the mutant allele no longer invades the population.

- The 'Numerical_analyses_Selfing' repository contains a program which performs iterations of the exacts recursions on genotypic frequencies for the evolution of self-fertilisation with a fixed reproductive effort.

--------------------------------------------------------------------------------------------------

All programs are coded in C++, and have a similar structure:

- Recursions are coded in the 'Recursion.cpp' file. 

- The 'ranbin.cpp' file contains algorithms used to sample from various probability distributions, such as Gamma, Poisson or Binomial distributions. 

- The 'mt.h' file is the header associated with the Mersenne Twister random numbers generator, which is used in all our programs for stochastic events. 

- The 'depression.h' file is the header that contains all the prototypes of all particular functions and structures used in our programs.

- The 'fichiers.cpp' file contains the function used to read parameters from the 'parameters.txt' file. In the 'parameters.txt' file, each line starting with a * corresponds to a new parameter set.


