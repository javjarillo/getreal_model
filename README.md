# Predictors of food web resistance to environmental change

Code (.R files and .py files) and data to reproduce the results in "Predictors of food web resistance to environmental change"', by Jarillo et al. (2025)

* `script_fromdb_newTAs_Tom_corrected_v2_PredPreyDiffLinks.py`
  Python script for simulating the meta-community dynamics of of six representative macroinvertebrate food-web assemblages found across European streams.
  This script requires the csv data files from data_GETREAL/newTAs/ and data_GETREAL/river_nws/

* `estimating_extinctions_random_matrices_v3_random_estimation_lognormal_histogram_include_min_max.R` 
  R script for estimating the taxon resistance within random communities.


## Resistance within random communities

R script for computing the taxon resistances when considering taxon interactions within random communities: `estimating_extinctions_random_matrices_v3_random_estimation_lognormal_histogram_include_min_max.R`

In this code:

*  We specify the number of prey (`nV_max`) and predator (`nP_max`) forming the random communities.
* We specify the mean values of the parameters
	* prey growth rate `b`
	* prey death rate `m`
	* predator death rate `d`
	* interspecific competition strength `alpha`
	* attack rate `lambda`
	* conversion efficiency `eta`
	* prey tolerance `tauV`
 	* predator tolerance `tauP` 
	* the coefficient of variation of all parameters `sigma`
	* whether the predators are generalists or specialists (encoded in the variable `mode_predators`).
	
Then, the code will assign for the different taxa random parameters for the different rates and interaction strength from a log-normal distributions, with the specified mean values and variances, that ensures that all taxa coexist in the absence of environmental change ($\xi = 0$).

Later on, for different strengths of the environmental change driver intensities $\xi$ (`zi`), and based on the random tolerances, the program will recalculate the taxa growth rates ($\vec{r} (\xi)$), and compute the equilibrium abundances ($\vec{n}_{eq} = A^{-1} \, \vec{r} (\xi)$, with $A$ a matrix containing the interactions between taxa$).
The program will also remove from the community those taxa going form a positive to zero abundance, computing for such environmental change $\xi$ the number of remaining prey and predators. Finally, the program will store for the different environmental changes the different diversities (number of prey and predators), together with other statistics about the remaining community (that has not been yet removed).

## Simulations of typical stream communities in realistic spatial networks

Python script for simulating the meta-community dynamics of realistic stream communities: `script_fromdb_newTAs_Tom_corrected_v2_PredPreyDiffLinks.py`

The program first load the traits of the taxa forming a specific typical macroinvertebrate community (or *typical assemblage*). 

Then, it loads a spatial river-network of a given number of nodes, which was previously generated in R employing the OCNet package.

Later on, with the function `SpeciesParameters_Troph_v3_method`, and based on the traits of the taxa, it categorizes them:
* As either *predators* or *prey*. Predator taxa are those taxa
	* known to feed on other living macroinvertebrates (`food_livingmacroinvertebrates_t >0`)
	* whose main food resource are actually other macroinvertebrates: 
	`food_livingmacroinvertebrates_t >= 2./3 * max(food_vertebrates_t, food_livingmicroinvertebrates_t, food_deadanimal_t, food_livingmacrophytes_t, food_livingmicrophytes_t, food_deadplant_t, food_detritus_t, food_micro_t)`
* as either *flying* or not flying taxa (*aquatic*) (based on the traits `aquatic_pas`, `aquatic_act`, `flying_pas`, `flying_act`. Depending on these traits, we will allow the taxa to move in a differential time step $dt$ just between connected nodes of the spatial network (*aquatic* dispersal), or between any pair of nodes (*flying* dispersal). 
  **NOTE:** If a taxa has non zero traits for the different categories  `aquatic_pas`, `aquatic_act`, `flying_pas`, `flying_act`, they will be able to disperse using all multiple ways.
	* The difference between *active* and *pasive* dispersal will rely on how the probability of going to a node $x_2$ from a node $x_1$ scales with the distance between the nodes. 
	  For active dispersal, the probability of reaching a node $x_2$ for individuals dispersing from a node $x_1$ scales with the inverse of the distance between the nodes.
	  For passive dispersal, the probability of reaching a node $x_2$ for individuals dispersing from a node $x_1$ scales with the **square** of inverse of the distance between the nodes. 
		* A flying taxon with active dispersal will be able to get to any other node $x_2$. Among all the individuals of the taxa dispersing from $x_1$, the fraction of individuals going to node $x_2$ will be given by 
		  
		  $$f_{F,Act} (x_1, x_2) = \frac{\frac{1}{distance(x_1, x_2)}}{\sum_{x'} \frac{1}{distance(x_1, x')}}$$
		* A flying taxon with passive dispersal will be able to get to any other node $x_2$. Among all the individuals of the taxa moving from $x_1$, the fraction of individuals going to node $x_2$ will be given by 
		  
		  $$f_{F,Pas} (x_1, x_2) = \frac{\frac{1}{(distance(x_1, x_2))^2}}{\sum_{x'} \frac{1}{(distance(x_1, x'))^2}}$$
		* An aquatic taxon with active dispersal will be able to move to any node $x_2$ **connected to the initial node** $x_1$. Among all the individuals of the taxa dispersing from $x_1$, the fraction of individuals going to node $x_2$ will be given by
		  
		  $$f_{Aq,Act} (x_1, x_2) = \frac{\frac{1}{distance(x_1, x_2)}}{\sum_{k' \textrm{ connected to } x_1} \frac{1}{distance(x_1, x')}}$$
		* An aquatic taxon with passive dispersal will be able to move to any node $x_2$ **connected to the initial node** $x_1$. Among all the individuals of the taxa dispersing from $x_1$, the fraction of individuals going to node $x_2$ will be given by 
		  
		  $$f_{Aq,Pas} (x_1, x_2) = \frac{\frac{1}{(distance(x_1, x_2))^2}}{\sum_{x' \textrm{ connected to } x_1} \frac{1}{(distance(x_1, x'))^2}}$$
		  

Now we explain how we set the parameters of the taxa for the simulations.
* **Dispersal rates** $m_i^{\textrm{mode k}}$. For all taxa $i$, we assigned them random dispersal rates $m_i^{\textrm{mode k}} \sim \langle \textrm{AverDispersal}^{\textrm{mode k}} \rangle \left[ 1+ N(0, \sigma) \right]$. By default, we set $\sigma=0.1$, and
	* $\textrm{AverDispersal}^{\textrm{(Aq,Pas)}} = 0.1$
	* $\textrm{AverDispersal}^{\textrm{(Aq,Act)}} = 0.2$
	* $\textrm{AverDispersal}^{\textrm{(F,Pas)}} = 0.1$
	* $\textrm{AverDispersal}^{\textrm{(F,Act)}} = 0.2$
* **Growth rates** $r_i$.
	* For prey taxa $i$, the program will assign a random growth rate $r_i \sim   \langle \textrm{AverGrowth} \rangle \left[ 1+ N(0, \sigma) \right]$ (default values: $\textrm{AverGrowth} \rangle = 1, \sigma = 0.1$)
	* For predator taxa, growth rate $r_i \sim   - \langle \textrm{AverGrowth} \rangle \left[ 1+ N(0, \sigma) \right]$ (default values: $\textrm{AverGrowth} \rangle = 1, \sigma = 0.1$). This growth rate is negative, reflecting that predator populations decreases in the absence of the prey.
* **Interaction strength** $A_{i,j}$
	* *Interactions between prey taxa*.
		* $A_{i,i} = \sim N(1, \sigma)$ (intraspecific competition is set equal to 1). Default value: $\sigma=0.1$
		* For the the inter-specific competition, the traits of the taxa are employed. 
			* If preys $i$ and $j$ are known to feed on very different food resources (obtained from the Traits in the Tachet database), the interaction is set to zero. $A_{i,j} = 0$.
			* Otherwise, if the prey taxa $i,j$ are reported to feed on similar food resources, the interaction strength is modulated based on the sizes of the taxa: 
			  
			  $$\begin{aligned} 
			  size_{i} &\sim N(\mu_i, \sigma_i) , \,
			  size_j \sim N(\mu_j, \sigma_j) \\
			  \textrm{mean} \left( A_{i,j} \right)&= \frac{\omega}{\sqrt{\omega^2 + 2 \sigma_i^2+2 \sigma_j^2}} e^{-\frac{(\mu_i-\mu_j)^2}{\omega^2 + 2 \sigma_i^2+2 \sigma_j^2}} \\ 
			  A_{ij} &\sim \textrm{mean} \left( A_{i,j} \right) \times N(1,\sigma)
			  \end{aligned}$$
			  
			  Thus, if the taxa present very different size, or there exist a high variability in the taxon sizes, the interaction between them tend to zero. Otherwise, if they present really similar sizes and the variability is low, the interaction strength is increased up to $\textrm{mean} \left( A_{i,j} \right) \approx 1 - \frac{\sigma_1^2 + \sigma_2^2}{\omega^2} - \frac{(\mu_1 - \mu_2)}{\omega^2}$. By default, we set $\omega=\frac{1}{2}$ and $\sigma=0.1$
	* *Interactions between the predator taxa.*
		* $A_{i,i} \sim N(1,\sigma)$ (intraspecific competition is set equal to 1).
		* For the inter-specific competition, $A_{i,j}=0$. (there is no direct competition between predator species; they can compete indirectly through consumption of similar resources)
	* *Interactions between the predators and the prey
		* Employing the GLOBI database, we can set which predator taxa have been reported to food on which other prey taxa. (In the files, such information is already encoded in the csv files `Tachet_NewTA_XXX_210317_gen_predation_connections.csv`).
		* If a predator taxon $i$ has not been reported to feed on a prey taxon $j$, the interactions strengths are set to zero ($A_{i,j} = A_{j,i} = 0$)
		* On the contrary, if such interactions has been reported, we modulate the interactions strength based on the size traits between the prey and the predator.
		  It has been set an optimum ratio of predator-to-prey size of $4$. When the fraction of the sizes differs from this optimum ratio, the attack rate is reduced.
		  The files `Tachet_NewTA_XXX_210317_gen_predation_connections.csv` already contains the reduction of the attack rate between the predators and the prey based on the sizes, $\lambda_{i,j}$. 
		  
		  $$\begin{aligned}
		  A_{i,j} &\sim \langle AverPredPrey \rangle \, \lambda_{i,j} \times N(1,\sigma) \\  
		  A_{j,i}& \sim -\langle AverPredPrey \rangle \, \lambda_{i,j} \, \langle EfficiencyPredation \rangle \times N(1,\sigma)
		  \end{aligned}
		  $$
		  
		  By default we set $\langle AverPredPrey \rangle = 1$ and $\langle EfficiencyPredation \rangle = 0.8$
* **Taxon tolerances and effect of environmental change**.
	* The tolerances of the different taxa against the different considered chemical pollutants (a metal, copper; a herbicide, atrazine; and a pesticide, imidacloprid) were obtained from
		* T. Sinclair, P. Craig, and L. L. Maltby. Climate warming shifts riverine macroinvertebrate communities to be more sensitive to chemical pollutants. Global Change Biology, 30 (4):e17254, 2024. https://doi.org/10.1111/gcb.17254
		* J. F. Jupke,  et al. Europe-wide spatial trends in copper and imidacloprid sensitivity of macroinvertebrate assemblages. Environmental Sciences Europe, 36(1):124, 2024. https://doi.org/10.1186/s12302-024-00944-3
		* J. F. Jupke. The capacity of broad-scale aquatic typology systems to capture differences in composition and chemical sensitivity of biological assemblages. PhD thesis, Rheinland-Pfälzische Technische Universität Kaiserslautern-Landau, 2024. https://doi.org/10.26204/KLUEDO/7665
	* From these references we obtained $EC_{50}$ values for the different taxa against the different pollutants.
	* Then, for a given chemical concentration $\xi$ of an specific pollutant, we alter the growth rates of prey and predator taxa.
		* For prey taxa $i$, its growth rate is decreased to $r_i^{prey} = r_i \, \left[e^{- \frac{\xi}{EC_i}} - 10^{-3} \right]$. The $10^{-3}$ term is included to ensure that for really high chemical concentrations the growth rate of the taxa are not zero, but negative.
		  
		  ![[prey_growth_rate_vs_zi.png]]
		  
		* For predator taxa $j$, its growth rate becomes more negative, $r_j^{predator} = r_j \, \left[2-e^{- \frac{\xi}{EC_j}} + 10^{-3} \right]$.
		  
		  ![[predator_growth_rate_vs_zi.png]]
		  
  
Once we have assigned to all prey and predator taxa these parameters, we can write the differential equations that govern the community dynamics. We will set then small time steps $\Delta t$, and we will then compute the local dynamics and the dispersal dynamics. For each time step, we assume that it first occurs the local dynamics, with the interaction of the taxa, and later on the dispersal dynamics of motion between the nodes of the spatial network. 
$$
\begin{aligned}  
\widetilde{N_i}(x,t)  &= N_i (x,t) + \Delta t \, \frac{dN_i^{(local)}(x,t)}{dt} \\  
N_i(x,t+\Delta t)  &= \widetilde{N_i} (x,t) + \Delta t \, \frac{d\widetilde{N_i}^{(disperal)} (x,t)}{dt}
\end{aligned}
$$

Where
$$
\frac{dN_i^{(local)}(x,t)}{dt} = N_{i} (x,t) \, \left[ r_i + \sum_{j=1}^{n_v+n_p} A_{i,j} \, N_{j} (x,t) \right]
$$
and
$$
\begin{aligned}
\frac{dN_i^{(dispersal)}(x,t)}{dt} =& - \left[ m_i^{F,Pas} + m_i^{F,Act} + m_i^{Aq,Pas} + m_i^{Aq,Act}\right] N_i (x,t) \\
&+ \sum_{y \neq x} \left [ m_i^{F,Pas} f_{F,Pas}(y,x) +  m_i^{F,Act} f_{F,Act}(y,x)\right] N_i (y,t) 
\\
&+ \sum_{y \neq x; \, y \textrm{ connected to } x} \left [ m_i^{Aq,Pas} f_{Aq,Pas}(y,x) +  m_i^{Aq,Act} f_{Aq,Act}(y,x)\right] N_i (y,t) 
\end{aligned}
$$

