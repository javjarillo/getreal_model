# Predictors of food web resistance to environmental change

Code (.R files and .py files) and data to reproduce the results in ``Predictors of food web resistance to environmental change'', by Jarillo et al. (2025)

* `script_fromdb_newTAs_Tom_corrected_v2_PredPreyDiffLinks.py`

  Python script for simmulating the meta-community dynamics of of six representative macroinvertebrate food-web assemblages found across European streams

  This script require the csv data files from data_GETREAL/newTAs

* `estimating_extinctions_random_matrices_v3_random_estimation_lognormal_histogram_include_min_max.R`

  R script for simmulating the community dynamics of random communities


## Simulation of random communities

R script for simmulating the community dynamics of random communities: `estimating_extinctions_random_matrices_v3_random_estimation_lognormal_histogram_include_min_max.R`

In this code:

*  We specify the number of prey (`nV_max`) and predator (`nP_max`) forming the random communities.
* We sepcify the mean values of the parameters
	* prey growth rate `b`
	* prey death rate `m`
	* predator death rate `d`
	* interespecific competiting strength `alpha`
	* attack rate `lambda`
	* conversion efficiency `eta`
	* prey tolerance `tauV`
 	* predator tolerance `tauP` 
	* the coefficient of variation of all parameters `sigma`
	* whetther the predators are generalists or specialists (encoded in the variable `mode_predators`).
	
Then, the code will assign for the different taxa random parameters for the different rates and interaction strength from a log-normal distributions, with the specified mean values and variances, that ensures that all taxa coexist in the absence of environmental change.

Later on, for different strengths of the environmental change driver (`zi`), and based ond the random tolerances, the program will recalculate the taxa growth rates, and compute the equilibrium abundances. The program will also remove from the community those taxa going form a positive to zero abundance, computing for such environmental change the number of prey and rpedators. Finally, the program will store for the different environmental changes the different diversities (number of prey and predators), together with other statistica about the remaining community (that has not been yet removed).


## Simulations of typical stream communities in realistic spatial networks

Python script for simmulating the meta-community dynamics of realistic stream communities: `script_fromdb_newTAs_Tom_corrected_v2_PredPreyDiffLinks.py`

The program receives 
