# Predictors of food web resistance to environmental change

Code (.R files and .py files) and data to reproduce the results in ``Predictors of food web resistance to environmental change'', by Jarillo et al. (2025)

* `script_fromdb_newTAs_Tom_corrected_v2_PredPreyDiffLinks.py`

  Python script for simmulating the meta-community dynamics of of six representative macroinvertebrate food-web assemblages found across European streams

  This script require the csv data files from data_GETREAL/newTAs

* `estimating_extinctions_random_matrices_v3_random_estimation_lognormal_histogram_include_min_max.R`

  R script for estimating the taxon resistance within random communities


## Resistance within random communities

R script for computing the tazon resistances when considering taxon interactions within random communities: `estimating_extinctions_random_matrices_v3_random_estimation_lognormal_histogram_include_min_max.R`

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
	
Then, the code will assign for the different taxa random parameters for the different rates and interaction strength from a log-normal distributions, with the specified mean values and variances, that ensures that all taxa coexist in the absence of environmental change ($\xi = 0$).

Later on, for different strengths of the environmental change driver intensities $\xi$ (`zi`), and based ond the random tolerances, the program will recalculate the taxa growth rates ($\vec{r} (\xi)$), and compute the equilibrium abundances ($\vec{n}_{eq} = A^{-1} \, \vec{r} (\xi)$, with $A$ a matrix containing the interactions between taxa$).
The program will also remove from the community those taxa going form a positive to zero abundance, computing for such environmental change $\xi$ the number of remaining prey and predators. Finally, the program will store for the different environmental changes the different diversities (number of prey and predators), together with other statistics about the remaining community (that has not been yet removed).

## Simulations of typical stream communities in realistic spatial networks

Python script for simmulating the meta-community dynamics of realistic stream communities: `script_fromdb_newTAs_Tom_corrected_v2_PredPreyDiffLinks.py`

The program first load the traits of the taxa forming a specific typical macroinvertebrate community (or *typical assemblage*). 

Then, it loads a spatial river-network of a given number of nodes, which was previously generated in R employing the OCNet package.

Later on, with the function `SpeciesParameters_Troph_v3_method`, and based on the traits of the taxa, it categorizes them:
* as either predators or prey. Predator taxa are those taxa
	* known to feed on other living macroinvertebrates (`food_livingmacroinvertebrates_t >0`)
	* whose main food reasource are actually other macroinvertebrates: 
	`food_livingmacroinvertebrates_t >= 2./3 * max(food_vertebrates_t, food_livingmicroinvertebrates_t, food_deadanimal_t, food_livingmacrophytes_t, food_livingmicrophytes_t, food_deadplant_t, food_detritus_t, food_micro_t)`
* as either flying or not flying taxa (based on the traits `aquatic_pas`, `aquatic_act`, `flying_pas`, `flying_act`. Depending on these traits, we will allow the taxa to move in a differential time step $dt$ just between connected nodes of the spatial network (*aquatic* dispersal), or between any pair of nodes (*flying* dispersal).
	* The difference between *active* and *pasive* dispersal will rely on how the probability of going to a node $k_2$ from a node $k_1$ scales with the distance between the nodes:
		* for active dispersal, the probability of reaching a node $k_2$ for individiuals dispersing from a node $k_1$ scales with the inverse of the distance between the nodes,
		* for passive dispersal, the probability of reaching a node $k_2$ for individiuals dispersing from a node $k_1$ scales with the **square** of inverse of the distance between the nodes. 

	* A flying taxon with active dispersal will be able to get to any other node $k_2$. Among all the individuals of the taxa dispersing from $k_1$, the fraction of individuals going to node $k_2$ will be given by
	
	$$f_{F,Act} (k_1, k_2) = \frac{1/distance(k_1, k_2)}{\sum_{k'} 1/distance(k_1, k')}$$
	
	* A flying taxon with passive dispersal will be able to get to any other node $k_2$. Among all the individuals of the taxa movingg from $k_1$, the fraction of individuals going to node $k_2$ will be given by
	
	$$f_{F,Pas} (k_1, k_2) = \frac{1/(distance(k_1, k_2))^2}{\sum_{k'} 1/(distance(k_1, k'))^2}$$
	
	* An aquatic taxon with active dispersal will be able to move to any node $k_2$ **connected to the initial node** $k_1$. Among all the individuals of the taxa dispersing from $k_1$, the fraction of individuals going to node $k_2$ will be given by 
	
	$$f_{Aq,Act} (k_1, k_2) = \frac{1/(distance(k_1, k_2))^2}{\sum_{k', \textrm{k' connected to }k_1} 1/distance(k_1, k')}$$
	
	* An aquatic taxon with passive dispersal will be able to move to any node $k_2$ **connected to the initial node** $k_1$. Among all the individuals of the taxa dispersing from $k_1$, the fraction of individuals going to node $k_2$ will be given by 
	
	$$f_{Aq,Pas} (k_1, k_2) = \frac{1/(distance(k_1, k_2))^2}{\sum_{k', \textrm{k' connected to }k_1} 1/(distance(k_1, k'))^2}$$
	
	

	
	


