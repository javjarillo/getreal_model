# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 09:26:56 2020

@author: jjarillo
"""

import datetime
import networkx as nx
import numpy as np
import os
import pandas
import random
import string
import sys
import time           
import tracemalloc

def signif(x, p):
    """
        signif(x, p).
        This function round a real number 'x' to the specific number of significative figures 'p' 
    """
    x = np.asarray(x)
    x_positive = np.where(np.isfinite(x) & (x != 0), np.abs(x), 10**(p-1))
    mags = 10 ** (p - 1 - np.floor(np.log10(x_positive)))
    return np.round(x * mags) / mags


def Aver_Alphas_given_traits_v2(traits_comp, aver_alphas, pars):
    """
        Aver_Alphas_given_traits(traits_comp, aver_alphas, sigma_size_comp = 0.5, beta_size_comp = 0.5,  reduct_alphas_food, reduct_alphas_feed).
        This function returns a modified aver_alphas matrix, with the competition strength coefficients, given a database with the traits which modulate the strength of the interaction. For this moment, the body size of the species; their prefered food sources; and the feeding behaviour.
        
        1) Effect of body_size on competition: \Delta alpha_{ij}_{rel} = \exp{sigma_size_comp^2 \, \beta_size_comp^2} \, \exp{ - \frac{1}{4 sigma_size_comp^2} \, \left( size_i -size_j + 2 \, sigma_size_comp^2 \, \beta_size_comp  \right)^2  }
        2) Effect of food resource on competition: \Delta alpha_{ij}_{rel} = \exp{   - \frac{1}{2 \, \sigma_food} \sum_{f=1}^{\# food_resources}  \left|  food_i - \food_j \right\}. If the species don't share any food resources, \Delta alpha_{ij}_{rel} = \exp{ - \frac{1}{\sigma_food} } = reduct_alphas_food < 1. If the species have the same resources with the same need ratio, then \Delta alpha_{ij}_{rel} = 1
        3) Effect of feeding behaviour on competition: \Delta alpha_{ij}_{rel} = \exp{   - \frac{1}{2 \, \sigma_feed} \sum_{f=1}^{\# feed_behaviours}  \left|  feed_i - \feed_j \right\}. If the species don't have the same feeding behaviours, \Delta alpha_{ij}_{rel} = \exp{ - \frac{1}{\sigma_feed} } = reduct_alphas_feed < 1. If the species have the same feeding behaviours, then \Delta alpha_{ij}_{rel} = 1
    """
    ## Effect of sizes on aver_alphas
    data_size= traits_comp[ [i for i in traits_comp.columns[traits_comp.columns.str.startswith('Size')]] ]
    aver_size = sum([(i+1) * data_size.iloc[:,i] for i in range(data_size.shape[1])]) / sum([ data_size.iloc[:,i] for i in range(data_size.shape[1])])    
    var_size  = sum([(i+1)**2 * data_size.iloc[:,i] for i in range(data_size.shape[1])]) / sum([ data_size.iloc[:,i] for i in range(data_size.shape[1])]) - aver_size**2
    effects_size = Effect_Traits_Continuous(aver_trait = aver_size, var_trait = var_size, omega = pars['omega_size']) 
    aver_alphas = aver_alphas * effects_size
    
    
    ## Effect of food on aver_alphas
    # reduct_alphas_food = 1/10 #average reduction in the alphas if the species don't share any food resource. 1/10 = 10\%
    traits_food = traits_comp[ [i for i in traits_comp.columns[traits_comp.columns.str.startswith('food')]] ] 
    # I normalize the value of the traits, so its sum is 1. I do not do that if that species has all the trait values equal to zero. zero arrays cannot be normalized.
    overlap_food = Effect_Traits_Categorical(traits = np.array(traits_food)) # cosine of the angles formed betwenn the food-resources vectors of the different species
    aver_alphas[overlap_food < pars['thershold_overlap_foods']] = 0.
                 
    return aver_alphas

def Categorize_disp(traits_disp):
    """
        Aver_Disp_Rates_traits(traits_disp, aver_disp_rate, ratio_disp_rate_fliers).
        Note1: aver_disp_rate should contain the average passive dispersal rate
        Note2: ratio_disp_rate_fliers should contain the ratio of passive aerial and aquative dispersal rates.
        Given the the traits that affects the dispersal rates, it returns:
        1) a boolean array, stating if the species move through the river system or not
        2) an array with the average dispersal rates of the species through the river
        3) a boolean array, stating if the species can fly (and move out of the river network system)
        4) an array with the average dispersal rates of the species out of the river network system
    """
        
    flying_act, flying_pas, aquatic_act, aquatic_pas  = np.zeros(traits_disp.shape[0]), np.zeros(traits_disp.shape[0]), np.zeros(traits_disp.shape[0]), np.zeros(traits_disp.shape[0])
    
    for i in range(traits_disp.shape[0]):
        if (traits_disp.disp_air_passive_t[i] + traits_disp.disp_air_active_t[i] + traits_disp.disp_aq_passive_t[i] + traits_disp.disp_aq_active_t[i] > 0):
            flying_act[i] = traits_disp.disp_air_active_t[i] / (traits_disp.disp_air_passive_t[i] + traits_disp.disp_air_active_t[i] + traits_disp.disp_aq_passive_t[i] + traits_disp.disp_aq_active_t[i])
            flying_pas[i] = traits_disp.disp_air_passive_t[i]/ (traits_disp.disp_air_passive_t[i] + traits_disp.disp_air_active_t[i] + traits_disp.disp_aq_passive_t[i] + traits_disp.disp_aq_active_t[i])
            aquatic_act[i]= traits_disp.disp_aq_active_t[i]  / (traits_disp.disp_air_passive_t[i] + traits_disp.disp_air_active_t[i] + traits_disp.disp_aq_passive_t[i] + traits_disp.disp_aq_active_t[i])
            aquatic_pas[i]= traits_disp.disp_aq_passive_t[i] / (traits_disp.disp_air_passive_t[i] + traits_disp.disp_air_active_t[i] + traits_disp.disp_aq_passive_t[i] + traits_disp.disp_aq_active_t[i])
    
    out = {'aquatic_pas': aquatic_pas, 'aquatic_act': aquatic_act,
           'flying_pas' : flying_pas,  'flying_act' : flying_act}
    return out

def Categorize_Traits(data_traits):
    """
        Categorize_Traits(data_traits).
        Given the database with all the traits, it returns two subdatabases:
        1) with the traits important for competition
        2) with the traits important for dispersal
        3) with the traits important for longitudinal distribution    
    """
    col_comp = []
    col_disp = []
    col_long_dist = []
    
    nrow, ncol = data_traits.shape
    
    TRAITS_COMP = load_database(os.path.join(os.getcwd(), 'getreal-model', 'code', 'traits_comp.csv')) 
    TRAITS_DISP = load_database(os.path.join(os.getcwd(), 'getreal-model', 'code', 'traits_disp.csv')) 
    TRAITS_LONG_DIST = load_database(os.path.join(os.getcwd(), 'getreal-model', 'code','traits_long_dist.csv'))
           
    for trait_comp in TRAITS_COMP.columns:
        if trait_comp in data_traits.columns:
            col_comp.append(trait_comp)
    for trait_disp in TRAITS_DISP.columns:
        if trait_disp in data_traits.columns:
            col_disp.append(trait_disp)
    for trait_long_dist in TRAITS_LONG_DIST.columns:
        if trait_long_dist in data_traits.columns:
            col_long_dist.append(trait_long_dist)
    return data_traits[col_comp], data_traits[col_disp], data_traits[col_long_dist]


def ComputeDynamics(pars, pars_nw, pars_dyn):
    ## Simulating the final population sizes of the different species at the 
    ## different nodes, for the times from 0 to time_sim
    if (pars['number_nodes']>1):
        times = float('nan')
        pops  = Compute_Pop_Sizes_v2(pars = pars, pars_nw = pars_nw, pars_dyn = pars_dyn)
    else:
        pops, times = Compute_Pop_Sizes_v2(pars = pars, pars_nw = pars_nw, pars_dyn = pars_dyn)
    
    final_pops = pops[:,:,-1]
    
    pops_corr = float('nan')
    ## Compute the average body_size for each node, and the fraction of flying species respect to the surviving species
    mean_sizes     = np.zeros((1, pars['number_nodes']))
    frac_flying    = np.zeros((1, pars['number_nodes']))
    frac_predators = np.zeros((1, pars['number_nodes']))
    frac_species_pred = np.zeros((1, pars['number_nodes']))
    
    
    flying = (np.array(pars_dyn['flying_pas'])>0) | (np.array(pars_dyn['flying_act'])>0)
    
    for k in range(pars['number_nodes']):
        if max(final_pops[:,k])>0:
            mean_sizes[0,k] = np.sum(final_pops[:,k] * pars_dyn['mean_size_species'][:])/np.sum(final_pops[:,k])
            frac_flying[0,k]= np.sum( flying[:]  * final_pops[:,k])/ np.sum(final_pops[:,k])
            frac_predators[0,k] = np.sum(pars_dyn['predator'][:] * final_pops[:,k])/ np.sum(final_pops[:,k])
            frac_species_pred[0,k] = np.sum( pars_dyn['predator'][:] & (final_pops[:,k]>0)  ) /  np.sum(final_pops[:,k]>0)
            
    
    ##Corrected values, erasing species with pop<limit_detec    
    final_pops_corr = final_pops.copy()
    final_pops_corr [final_pops_corr <= pars['limit_detec'] ] = 0
    if (pars['number_nodes']==1):
        pops_corr = pops.copy()
        pops_corr [pops_corr <= pars['limit_detec']] =0
    
    ## Compute the average body_size for each node (corrected values, erasing species with pop<limit_detec)
    mean_sizes_corr     = np.zeros((1, pars['number_nodes']))
    frac_flying_corr    = np.zeros((1, pars['number_nodes']))
    frac_predators_corr = np.zeros((1, pars['number_nodes']))
    frac_species_pred_corr = np.zeros((1, pars['number_nodes']))
    for k in range(pars['number_nodes']):
        if max(final_pops_corr[:,k]>0):
            mean_sizes_corr[0,k] = np.sum(final_pops_corr[:,k]*pars_dyn['mean_size_species'][:])/np.sum(final_pops_corr[:,k])
            frac_flying_corr[0,k]= np.sum( flying[:]  * final_pops_corr[:,k])/ np.sum(final_pops_corr[:,k])
            frac_predators_corr[0,k] = np.sum(pars_dyn['predator'][:] * final_pops_corr[:,k])/ np.sum(final_pops_corr[:,k])
            frac_species_pred_corr[0,k] = np.sum( pars_dyn['predator'][:] & (final_pops_corr[:,k]>0)  ) / np.sum(final_pops_corr[:,k]>0)
        
    output = {
        'final_pops': final_pops,
        'mean_sizes': mean_sizes,      
        'frac_flying': frac_flying, 
        'frac_predators': frac_predators, #fraction of predators,weithed by the species population sizes. Eg: 1 predator, size=5; 2 prey, size=10. frac = 5/(10+10+5) = 5/25 = 0.2
        'frac_species_pred': frac_species_pred, #fraction of predator species, NOT weigthed by the species pop. sizes. Eg: 1 predator, size 5, 2 prey, size = 10. frac = 1 predator/(3 species) = 0.33
        'final_pops_corr': final_pops_corr, 
        'mean_sizes_corr': mean_sizes_corr, 
        'frac_flying_corr': frac_flying_corr, 
        'frac_predators_corr': frac_predators_corr, 
        'frac_species_pred_corr': frac_species_pred_corr
        }
    
    if (pars['number_nodes']==1):
        output['pops'] = pops
        output['times'] = times
    
    return output


def Compute_Pop_Sizes_v2(pars, pars_nw, pars_dyn):     
    """
        This function integrate the population dynamics of the species inhabitating
        a graph with different nodes, at different distances, and with different edges
        between them.
    """
    
    inv_dist_air = np.zeros(np.shape(pars_nw['dist_air']))
    for k1 in range(pars['number_nodes']):
        for k2 in range(pars['number_nodes']):
            if (pars_nw['dist_air'][k1,k2]>0):
                inv_dist_air[k1,k2] = 1./pars_nw['dist_air'][k1,k2]
    
    pop = np.reshape( np.zeros(pars['number_species'] * pars['number_nodes']) , 
                     (pars['number_species'], pars['number_nodes']) )
    for k in range(pars['number_nodes']):
        for i in range(pars['number_species']):
            pop[i,k] = 1./(np.max(np.abs(pars_dyn['Alphas'][i,i])) + 1e-3)

    times = np.linspace(0, pars['time_sim'], int(pars['time_sim'])+1)
        
    ## Environmental variance of the environmental stochasticity 
    sigmas = pars['sigma_env'] * np.ones((pars['number_species'], pars['number_nodes']))
        
    args = {
            'rs': pars_dyn['rs'],
            'Alphas': pars_dyn['Alphas'],
            'Alphas_x': pars_dyn['Alphas_x'],
            'bool_polluted_node': pars_nw['bool_polluted_node'],
            'aquatic_pas': pars_dyn['aquatic_pas'], 
            'movs_aq_pas': pars_dyn['movs_aq_pas'],
            'aquatic_act': pars_dyn['aquatic_act'], 
            'movs_aq_act': pars_dyn['movs_aq_act'],
            'flying_pas': pars_dyn['flying_pas'], 
            'movs_fly_pas': pars_dyn['movs_fly_pas'],
            'flying_act': pars_dyn['flying_act'], 
            'movs_fly_act': pars_dyn['movs_fly_act'],
            'edges_norm_dist':  pars_nw['edges_norm_dist'],
            'inv_dist_air': inv_dist_air,
            'sigmas': sigmas,
            'width': pars_nw['width'],
            'vel': pars_nw['vel'],
            'MSS': pars_dyn['MSS'],
            'k_chem': pars_dyn['k_chem'],
            'chem_conc':   pars_nw['chem_conc'],
            'EC': pars_dyn['EC']
            }

        
    ## Stochastic integration"
    pops = np.zeros((pars['number_species'], pars['number_nodes'], len(times)))
    cont_t = 0
    Incr_prec_t = int(25)
    dt = (times[1]-times[0])/Incr_prec_t
    reduct_det_growth = np.zeros((pars['number_species'], pars['number_nodes']))
    env_noise_variance_growth = np.zeros((pars['number_species'], pars['number_nodes']))
    for sp in range(pars['number_species']):
        reduct_det_growth[sp,:] = -1. * args['k_chem'][sp] * args['rs'][sp] *  (1. - np.exp(-1. * args['chem_conc'][:]/args['EC'][sp])) + 1E-3 * np.abs(args['rs'][sp])
    
    cont_t_precision = 0    
    for t in np.linspace(0, pars['time_sim'], Incr_prec_t*int(pars['time_sim'])+1):
        pop = pop \
            +  dt * (np.reshape(dNdt_OCNet_v2_local(t=t, N=pop, args = args), (pars['number_species'], pars['number_nodes'])) 
                     - reduct_det_growth*pop) 
        pop[np.isfinite(pop)==False]=0.
        pop[pop < pars['limit_detec']] = 0.
        pop = pop + dt * np.reshape(dNdt_OCNet_v2_disp(t=t, N=pop, args = args), (pars['number_species'], pars['number_nodes']))
        pop[np.isfinite(pop)==False]=0.
        pop[pop < pars['limit_detec']] = 0.
        
        if (cont_t_precision == 0):
            pops[:,:,cont_t] = np.reshape(pop, (pars['number_species'], pars['number_nodes']))
            cont_t = cont_t + 1
        cont_t_precision = (cont_t_precision + 1) % Incr_prec_t

    if (pars['number_nodes']>1):
        return pops
    else:
        return pops, times

def dNdt_OCNet_v2_local(t, N, args):
    """Population dynamics equations.
     dNdt(t, N, args)
     args: dict with elements rs, Alphas, movs_aq_pas, movs_aq_act, movs_fly_pas,
     movs_fly_act, edges_norm_dist, inv_dist_air"""
    N = np.reshape(N,  (np.shape(args['Alphas'])[0], np.shape(args['edges_norm_dist'])[0])) #Bool_Polluted_nodes
    
    dndt = (
        np.einsum("ik,ik -> ik", args['rs'], N)  # Exponential growth rate
        - np.einsum("ijk,ik,jk -> ik",args['Alphas_x'],N,N) #density regulation 
        )
    return np.reshape(dndt, np.shape(args['Alphas'])[0]* np.shape(args['edges_norm_dist'])[0])  
    
def dNdt_OCNet_v2_disp(t, N, args):
    """Population dynamics equations.
     dNdt(t, N, args)
     args: dict with elements rs, Alphas, movs_aq_pas, movs_aq_act, movs_fly_pas,
     movs_fly_act, edges_norm_dist, inv_dist_air"""
    N = np.reshape(N,  (np.shape(args['Alphas'])[0], np.shape(args['edges_norm_dist'])[0])) #Bool_Polluted_nodes
    
    dndt = (
        # moving from the patch to connected patches for aquatic species (passive)
        - np.einsum("ikl,ik->ik", args['movs_aq_pas'], N)
        # moving to the patch from connected patches for aquatic species (passive)              
        + np.einsum("ilk,il->ik",args['movs_aq_pas'], N)
        # moving from the patch to connected patches for aquatic species (active)
        - np.einsum("ikl,ik->ik", args['movs_aq_act'], N)
        # moving to the patch from connected patches for aquatic species (active)              
        + np.einsum("ilk,il->ik", args['movs_aq_act'], N)
        # for flying species, moving from the patch to neigbour patches (passive)           
        - np.einsum("ikl,ik->ik", args['movs_fly_pas'], N)
        # for flying species, moving to the patch from neighbour patches  (passive)    
        + np.einsum("ilk,il->ik", args['movs_fly_pas'], N) 
        # for flying species, moving from the patch to neigbour patches (active)           
        - np.einsum("ikl,ik->ik", args['movs_fly_act'], N)
        # for flying species, moving to the patch from neighbour patches  (active)    
        + np.einsum("ilk,il->ik", args['movs_fly_act'], N)   
        )
    return np.reshape(dndt, np.shape(args['Alphas'])[0]* np.shape(args['edges_norm_dist'])[0])  



def Effect_Traits_Categorical(traits):
    effect = np.ones((np.shape(traits)[0],np.shape(traits)[0]))
    for i in range(np.shape(traits)[0]):
        for j in range(np.shape(traits)[0]):
            effect[i,j] = np.sum(traits[i,:] * traits[j,:])/np.sqrt(np.sum(traits[i,:]**2) * np.sum(traits[j,:]**2))
    return effect 


def Effect_Traits_Continuous(aver_trait, var_trait , omega):
    effect = np.ones((np.shape(aver_trait)[0],np.shape(aver_trait)[0]))
    for i in range(np.shape(aver_trait)[0]):
        for j in range(np.shape(aver_trait)[0]):
            effect[i,j] = omega/np.sqrt(2*var_trait[i] + 2*var_trait[j] + omega**2) * \
                np.exp(-1. * (aver_trait[i] - aver_trait[j])**2 * \
                       1./(2*var_trait[i] + 2*var_trait[j] + omega**2))
    return effect



def GenerateSpatialNetwork_OCNet_MaxReach_presets(pars, seed=np.int(datetime.datetime.now().strftime("%f"))):
    pars['number_nodes'] = pars['approx_number_nodes']
    repl_river = random.randrange(start=1,stop=11, step=1)
    name_river = '_'.join(['river', str(pars['number_nodes']), 'nodes',
                           str(repl_river), 'repl'])
    
    data_nodes =load_database(os.path.join(os.getcwd(), 'getreal-model', 'data_GETREAL','river_nws',
                                           '_'.join([name_river, 'NODES.csv'])))
    
    data_edges =load_database(os.path.join(os.getcwd(), 'getreal-model', 'data_GETREAL','river_nws',
                                           '_'.join([name_river, 'EDGES.csv'])))
    
    nodes     = np.array(data_nodes['nodes']) - 1 #-1, because in r the first node is the node 1, while in python is the node 0
    posoutlet = np.array(data_nodes['posoutlet']) - 1
    posoutlet = posoutlet[0]
    Xnode     = np.array(data_nodes['Xnode'])
    Ynode     = np.array(data_nodes['Ynode'])
    width     = np.array(data_nodes['width'])
    depth     = np.array(data_nodes['depth'])
    vel       = np.array(data_nodes['vel'])
    lengs     = np.array(data_nodes['lengs'])
    edges     = np.array(data_edges)
    
    G = nx.generators.empty_graph()
    G.add_nodes_from(nodes)
    ranknodes = np.zeros(len(nodes), dtype=int )
    dist_air = np.zeros((len(nodes),len(nodes)))
    for node1 in nodes:
        for node2 in range(node1, np.max(nodes)+1):
            dist_air[node1,node2] = np.sqrt( (Xnode[node1]-Xnode[node2])**2 + (Ynode[node1]-Ynode[node2])**2  )
            if (edges[node1,node2]>0):
                G.add_edge(node1, node2)
    dist_air = symmetrize(dist_air)
    for node1 in nodes:
        ranknodes[node1] = nx.shortest_path_length(G,source=posoutlet)[node1]
    pos = np.zeros((2,len(nodes)))
    pos[0,:] = Xnode
    pos[1,:] = Ynode
    
    dist_river = np.zeros(np.shape(dist_air))
    for node1 in nodes:
        for node2 in range(node1+1, np.max(nodes)+1):
            path_nodes = nx.shortest_path(G, source=node1, target=node2)
            min_pos = np.array(path_nodes)[np.min(ranknodes[path_nodes])==ranknodes[path_nodes]]
            path_nodes.remove(min_pos)
            dist_river[node1,node2] = np.sum(lengs[path_nodes])
    dist_river =   symmetrize(dist_river)                               
    
    if len(nodes)==1:
        connectivities = np.array([0])
    else:
        connectivities = np.array(G.degree())[:,1]
        
    bool_polluted_node = np.zeros(len(nodes))
    bool_polluted_node = np.zeros(len(nodes))
    temp = np.random.uniform(0,1,len(nodes))
    bool_polluted_node[temp< pars['frac_polluted_nodes']] = 1
    concs = bool_polluted_node * pars['conc_polluted_nodes']
    
    pars_nw = {
        'G': G, 
        'node': nodes, 
        'ranknode': ranknodes, 
        'connectivity': connectivities , 
        'pos': pos, 
        'edges': edges,
        'width': width,
        'depth': depth,
        'vel': vel,
        'dist_river': dist_river,
        'dist_air': dist_air,
        'bool_polluted_node': bool_polluted_node, 
        'chem_conc': concs, 
        'repl_river': repl_river
        }
    
    pars_nw['vel'] = pars['vel_river'] * pars_nw['vel']
    
    edges_norm_dist = np.zeros( np.shape( pars_nw['edges'] ))
    if len(nodes)==1:
        pars_nw['edges_norm_dist'] = edges_norm_dist
    else:
        for k1 in range(len(nodes)):
            for k2 in range(len(nodes)):
                if (pars_nw['edges'][k1,k2]>0):
                    edges_norm_dist[k1,k2] = edges[k1,k2] / dist_river[k1,k2]
        pars_nw['edges_norm_dist'] = edges_norm_dist
    
    return pars_nw



def is_non_zero_file(fpath):
    """is_non_zero_file(fpath)
    Return True if the file at fpath exists and is not empty"""
    return os.path.isfile(fpath) and os.path.getsize(fpath)>0

def load_database(filecsv, index_col=False):
    """load_database(filecsv). It returns the database contained in the csv file.
    IMPORTANT: filecsv should be a csv file"""
    if is_non_zero_file(filecsv)==False:
        sys.exit(' '.join(['The csv file', filecsv, 'does not exist']))
    else:
        data = pandas.read_csv(filecsv, sep = ';', encoding='latin-1', index_col=index_col)
    return data

def Random_Alphas(aver_alpha, sigma):
    """
        Random_Alphas(aver_alpha, sigma)
        Given values of average strength competitive interaction, alpha_{ij}, it returns a random number around this value with standard deviation sigma
    """
    number_species = aver_alpha.shape[0]
    alphas = np.zeros((number_species,number_species))
    for i in range(number_species):
        for j in range(number_species):
            #np.random.seed(seed)
            alphas[i,j] = np.random.normal(  aver_alpha[i,j], sigma*np.abs(aver_alpha[i,j]), 1  ) 
    
    return alphas

def RandString(k=8):
    return (''.join(random.choice(string.ascii_uppercase + string.ascii_lowercase + string.digits) for _ in range(int(k))))

def ShannonBiodiv(pops):
    pops[pops<0] = 0
    number_species, number_nodes = np.shape(pops)
    H = np.zeros(number_nodes)
    for k in range(number_nodes):
        if (np.sum(pops[:,k])>0):
            pis = pops[:,k]/np.sum(pops[:,k])
            pis = pis[pis>0]
            if (pis.size > 1):
                H[k] = - np.sum(pis*np.log(pis))
            else:
                H[k] = 0
        else:
            H[k] = 0
    return H

def sigmoid(x):
    return 1 / (1 + np.exp(-x))

def symmetrize(a):
    """symmetrize(a).
    Symetrize any matrix array a"""
    return a + a.T - np.diag(a.diagonal())

def SpeciesParameters_Troph_v3_method(pars, pars_nw, pars_chemres_size1, data_traits, method):
    ## traits for competition and dispersal
    traits_comp, traits_disp, traits_long_dist = Categorize_Traits(data_traits = data_traits)
    
    # BoolPredator will contain a boolean array stating if the species is a predator (True) of a Prey (False). 
    # We assume there exist just 2 trophic levels.
    # For being predator, the maximum food_trait should be food_livingmacroinvertebrates_t
    BoolPredator = (traits_comp.food_livingmacroinvertebrates_t > 0)  & (\
        ((traits_comp.food_livingmacroinvertebrates_t >= traits_comp.food_vertebrates_t) & \
         (traits_comp.food_livingmacroinvertebrates_t >= traits_comp.food_livingmicroinvertebrates_t) & \
         (traits_comp.food_livingmacroinvertebrates_t >= traits_comp.food_deadanimal_t) & \
         (traits_comp.food_livingmacroinvertebrates_t >= traits_comp.food_livingmacrophytes_t) & \
         (traits_comp.food_livingmacroinvertebrates_t >= traits_comp.food_livingmicrophytes_t) & \
         (traits_comp.food_livingmacroinvertebrates_t >= traits_comp.food_deadplant_t) & \
         (traits_comp.food_livingmacroinvertebrates_t >= traits_comp.food_detritus_t) & \
         (traits_comp.food_livingmacroinvertebrates_t >= traits_comp.food_micro_t)) | \
        (traits_comp.food_livingmacroinvertebrates_t >= 2./3 * np.max(
            traits_comp[['food_vertebrates_t', 'food_livingmicroinvertebrates_t', 'food_deadanimal_t', 
                         'food_livingmacrophytes_t', 'food_livingmicrophytes_t', 'food_deadplant_t', 
                         'food_detritus_t', 'food_micro_t']], 
            axis=1)
        ))
        
        
    # We distinguish between up and down specialist, based on the long_dist trait.
    # Thet are: up specialist, down specialist, generalist.
    BoolSpecialist_Up   = np.zeros(pars['number_species'], dtype=bool)
    BoolSpecialist_Down = np.zeros(pars['number_species'], dtype=bool)
    for i in range(pars['number_species']):
        suma = (traits_long_dist['long_dist_crenon_t'][i] + 
            traits_long_dist['long_dist_epirithron_t'][i] +
            traits_long_dist['long_dist_metarithron_t'][i] + 
            traits_long_dist['long_dist_hyporithron_t'][i] + 
            traits_long_dist['long_dist_epipotamon_t'][i] + 
            traits_long_dist['long_dist_metapotamon_t'][i] + 
            traits_long_dist['long_dist_estuary_t'][i])
        if (suma>0):
            up = (1/suma) * (traits_long_dist['long_dist_crenon_t'][i] + 
                traits_long_dist['long_dist_epirithron_t'][i] +
                traits_long_dist['long_dist_metarithron_t'][i] + 
                traits_long_dist['long_dist_hyporithron_t'][i]
                )
            down = (1/suma) * (traits_long_dist['long_dist_epipotamon_t'][i] + 
                traits_long_dist['long_dist_metapotamon_t'][i] + 
                traits_long_dist['long_dist_estuary_t'][i]
                )
            if (up>0.80):
                BoolSpecialist_Up[i] = True
            elif (down>0.5):
                BoolSpecialist_Down[i] = True
    
    
    disp_mode = Categorize_disp(traits_disp = traits_disp)
    aquatic_pas, aquatic_act = disp_mode['aquatic_pas'],  disp_mode['aquatic_act']
    flying_pas,flying_act    = disp_mode['flying_pas'], disp_mode['flying_act']
         
    rs = np.ones((pars['number_species'], pars['number_nodes']))
    
    
    inv_dist_air = np.zeros(np.shape(pars_nw['dist_air']))
    for k1 in range(pars['number_nodes']):
        for k2 in range(pars['number_nodes']):
            if (pars_nw['dist_air'][k1,k2]>0):
                inv_dist_air[k1,k2] = 1./pars_nw['dist_air'][k1,k2]
    
    frac_mov_between_nodes_fly_pas = np.zeros((pars['number_nodes'], pars['number_nodes']))
    frac_mov_between_nodes_fly_act = np.zeros((pars['number_nodes'], pars['number_nodes']))
    frac_mov_between_nodes_aq_pas  = np.zeros((pars['number_nodes'], pars['number_nodes']))
    frac_mov_between_nodes_aq_act  = np.zeros((pars['number_nodes'], pars['number_nodes']))
    for node1 in range(pars['number_nodes']):
        possible_end = np.array(range(pars['number_nodes']))
        possible_end = possible_end[possible_end != node1]
        if np.sum(inv_dist_air[node1, possible_end]) > 0:
            for node2 in range(pars['number_nodes']):
                frac_mov_between_nodes_fly_act[node1,node2] = inv_dist_air[node1,node2]/np.sum(inv_dist_air[node1, possible_end])
                frac_mov_between_nodes_fly_pas[node1,node2] = inv_dist_air[node1,node2]**2/np.sum(inv_dist_air[node1, possible_end]**2)
        else:
            for node2 in range(pars['number_nodes']):
                frac_mov_between_nodes_fly_act[node1,node2] = 1.
                frac_mov_between_nodes_fly_pas[node1,node2] = 1.
        
        possible_end = np.array(range(pars['number_nodes']))
        possible_end = possible_end[possible_end != node1]
        if np.sum(pars_nw['edges_norm_dist'][node1, possible_end]) > 0:
            for node2 in range(pars['number_nodes']):
                frac_mov_between_nodes_aq_act[node1,node2] = pars_nw['edges_norm_dist'][node1,node2]/np.sum(pars_nw['edges_norm_dist'][node1, possible_end])
                frac_mov_between_nodes_aq_pas[node1,node2] = pars_nw['edges_norm_dist'][node1,node2]**2/np.sum(pars_nw['edges_norm_dist'][node1, possible_end]**2)
        else:
            for node2 in range(pars['number_nodes']):
                frac_mov_between_nodes_aq_act[node1,node2] = 1.
                frac_mov_between_nodes_aq_pas[node1,node2] = 1.
       
                
                
    
    movs_aq_pas = np.zeros((pars['number_species'], pars['number_nodes'], pars['number_nodes'])) #dispersal rates of the species through river system, from the different nodes, to another nodes
    movs_aq_act = np.zeros((pars['number_species'], pars['number_nodes'], pars['number_nodes'])) #dispersal rates of the species through river system, from the different nodes, to another nodes
    movs_fly_pas= np.zeros((pars['number_species'], pars['number_nodes'], pars['number_nodes'])) #Out-of-Network (flying) dispersial rates of the species, at the different nodes
    movs_fly_act= np.zeros((pars['number_species'], pars['number_nodes'], pars['number_nodes'])) #Out-of-Network (flying) dispersial rates of the species, at the different nodes
    
    sigma_r = pars['sigma_param']
    sigma_mov = pars['sigma_param']
    sigma_mov_fly = pars['sigma_param']
    for i in range(pars['number_species']):
        rs[i,:] = np.random.normal(pars['aver_growth_rate'], 
                                   sigma_r*pars['aver_growth_rate'], 
                                   pars['number_nodes']  ) 
        
        if (BoolPredator[i]==True): # If the species are predators, their explicit growth rate must be negative. Otherwise, the predator would kill all the prey species
            rs[i,:] = -1.0 * np.abs(rs[i,:])
        
        # Here I add that some species are specialist at some river zones (upstream or downstream)
        # The specialist grows faster than the generalist in their prefered habitat, but cannot survive in the non prefered habitat
        if (BoolSpecialist_Up[i]==True):
            for k1 in range( pars['number_nodes']):
                if (pars_nw['dist_to_0'][k1] > 0.75 * np.max(pars_nw['dist_to_0'])): #If the species is in the prefered habitat: it grows at doble rate (prey) or decrease at half rate (predator)
                    if (BoolPredator[i]==False):
                        rs[i,k1] = rs[i,k1] * (1+pars['increase_growth_prefer_hab'])
                    else:
                        rs[i,k1] = rs[i,k1] + np.abs(rs[i,k1])*pars['increase_growth_prefer_hab']
                else: #elif (pars_nw['dist_to_0'][k1] < 0.25 * np.max(pars_nw['dist_to_0'])): 
                    #If the species is not in the prefered habitat: 
                    if (BoolPredator[i]==False):
                        rs[i,k1] = rs[i,k1] * (1 - pars['decrease_growth_not_prefer_hab'])  #-1. * np.abs(rs[i,k1])
                    else:
                        rs[i,k1] = rs[i,k1] - np.abs(rs[i,k1])*pars['decrease_growth_not_prefer_hab']
       
        elif (BoolSpecialist_Down[i]==True):
            for k1 in range( pars['number_nodes']):
                if (pars_nw['dist_to_0'][k1] < 0.25 * np.max(pars_nw['dist_to_0'])): #If the species is in the prefered habitat: it grows at doble rate (prey) or decrease at half rate (predator)
                    if (BoolPredator[i]==False):
                        rs[i,k1] = rs[i,k1] * (1+pars['increase_growth_prefer_hab'])
                    else:
                        rs[i,k1] = rs[i,k1] + np.abs(rs[i,k1])*pars['increase_growth_prefer_hab']
                else: #elif (pars_nw['dist_to_0'][k1] > 0.75 * np.max(pars_nw['dist_to_0'])): 
                    #If the species is not in the prefered habitat
                    if (BoolPredator[i]==False):
                        rs[i,k1] = rs[i,k1] * (1 - pars['decrease_growth_not_prefer_hab']) #-1. * np.abs(rs[i,k1])
                    else:
                        rs[i,k1] = rs[i,k1] - np.abs(rs[i,k1])*pars['decrease_growth_not_prefer_hab']
        
        
        
        
        for k1 in range( pars['number_nodes']):
            disp_rate_fly_pas = np.random.normal(pars['aver_disp_rate_fly_pas'], sigma_mov_fly * pars['aver_disp_rate_fly_pas'], 1)
            disp_rate_fly_act = np.random.normal(pars['aver_disp_rate_fly_act'], sigma_mov_fly * pars['aver_disp_rate_fly_act'], 1)
            
            disp_rate_aq_pas = np.random.normal(pars['aver_disp_rate_aq_pas'], sigma_mov * pars['aver_disp_rate_aq_pas'], 1)
            disp_rate_aq_act = np.random.normal(pars['aver_disp_rate_aq_act'], sigma_mov * pars['aver_disp_rate_aq_act'], 1)           
            
            for k2 in range( pars['number_nodes']):
                if (pars['number_nodes']==1):
                     movs_fly_pas[i,k1,k2] = 0.
                     movs_aq_pas[i,k1,k2]  = 0.
                     movs_fly_act[i,k1,k2] = 0.
                     movs_aq_act[i,k1,k2] = 0.
                else:
                    movs_fly_pas[i,k1,k2] = disp_mode['flying_pas'][i]  * disp_rate_fly_pas * frac_mov_between_nodes_fly_pas[k1,k2]
                    movs_aq_pas[i,k1,k2]  = disp_mode['aquatic_pas'][i] * disp_rate_aq_pas  * frac_mov_between_nodes_aq_pas[k1,k2]                
                    movs_fly_act[i,k1,k2] = disp_mode['flying_act'][i]  * disp_rate_fly_act * frac_mov_between_nodes_fly_act[k1,k2]
                    movs_aq_act[i,k1,k2]  = disp_mode['aquatic_act'][i] * disp_rate_aq_act  * frac_mov_between_nodes_aq_act[k1,k2]
                
                if (pars_nw['dist_to_0'][k1] < pars_nw['dist_to_0'][k2]): # Moving Against the Current. NOW WE DO NOT DISTINGUISH BETWEEN AGAINST AND FOLLOWING THE CURRENT MOVEMENTS
                    movs_aq_pas[i,k1,k2] = movs_aq_pas[i,k1,k2]  #np.max([movs_aq_pas[i,k1,k2] - pars_nw['vel'][k1], 0.]) # If we want to distinguish
                    movs_aq_act[i,k1,k2] = movs_aq_act[i,k1,k2]  #np.max([movs_aq_act[i,k1,k2] - pars_nw['vel'][k1], 0.]) # If we want to distinguish
                    movs_fly_act[i,k1,k2] = movs_fly_act[i,k1,k2] #* (1.0 - pars['ratio_disp_fly_current_orientation'])  # If we want to distinguish
                    movs_fly_act[i,k1,k2] = np.max([movs_fly_act[i,k1,k2], 0.])
                else: # Moving following the current
                    movs_aq_pas[i,k1,k2] = movs_aq_pas[i,k1,k2]  #movs_aq_pas[i,k1,k2] + pars_nw['vel'][k1] # If we want to distinguish
                    movs_aq_act[i,k1,k2] = movs_aq_act[i,k1,k2] #movs_aq_act[i,k1,k2] + pars_nw['vel'][k1] # If we want to distinguish
                    movs_fly_act[i,k1,k2] = movs_fly_act[i,k1,k2] #movs_fly_act[i,k1,k2]* pars['ratio_disp_fly_current_orientation'] # If we want to distinguish
                
                          
    ## Interaction strength of the species. The following ones are computed for a river stream width = 1.
    ## Since the different nodes will have different widths, later the values at the different nodes will be 
    ## modulated by this width, stored at pars_nw['width']
    
      
    # Computation of the average alpha: the averages interaction strengths at prsitine river nodes with width=1
    
    # Competition between prey species       
    aver_alphas = np.ones((pars['number_species'],pars['number_species']))
    aver_alphas = Aver_Alphas_given_traits_v2(traits_comp, aver_alphas, pars)
    
    # there is no direct competition between predator species 
    ## (they can compete indirectly through consumption of similar resources)
    for sp1 in range(pars['number_species']):
        for sp2 in range(sp1+1,pars['number_species']):
            if ((BoolPredator[sp1]==True) | (BoolPredator[sp2]==True)):
                aver_alphas[sp1,sp2] = 0.
                aver_alphas[sp2,sp1] = 0.
                    

    if (method == "GLOBI_Traits"):
        file_pred_inter = os.path.join(os.getcwd(),'getreal-model','data_GETREAL','newTAs',
                          ''.join(['Tachet_NewTA_', pars['TA'], '_210317_gen_predation_connections.csv'])) 
    elif (method == "EatAll_Traits"):
        file_pred_inter = os.path.join(os.getcwd(),'getreal-model','data_GETREAL','newTAs',
                          ''.join(['Tachet_NewTA_', pars['TA'], '_210317_gen_predation_connections_eatall_traits.csv'])) 
    elif (method == "EatAll_NOTRAITS"):
        file_pred_inter = os.path.join(os.getcwd(),'getreal-model','data_GETREAL','newTAs',
                          ''.join(['Tachet_NewTA_', pars['TA'], '_210317_gen_predation_connections_eatall_NOTRAITS.csv'])) 
    
    data_pred_inter = load_database(file_pred_inter, index_col=0) 
    alphas_preda = pandas.DataFrame(np.array(data_pred_inter), 
                                    columns=list(data_pred_inter.columns), 
                                    index=list(data_pred_inter.columns))
    TAXA_NAMES = np.array(traits_comp['Taxa'])
    aver_alphas_pred = np.zeros((pars['number_species'],pars['number_species']))
    for tx1 in range(len(TAXA_NAMES)):
        for tx2 in range(len(TAXA_NAMES)):
            aver_alphas_pred[tx1, tx2] = alphas_preda[TAXA_NAMES[tx2]][TAXA_NAMES[tx1]]
    aver_alphas_pred = pars['aver_pred_prey'] * aver_alphas_pred
    
    aver_alphas = aver_alphas + aver_alphas_pred
                
    # aver_alpha contains the averages interaction strengths, based on the traits, at pristine nodes with width = 1. 
    # However, maybe the interaction strengths is not exactly equal to the one estimated from the traits. 
    # Then we consider this kind of variation.
    Alphas_pristine = np.zeros((pars['number_species'],pars['number_species']))
    Alphas_pristine[:,:] = Random_Alphas(aver_alpha = aver_alphas, sigma = pars['sigma_param'])
    Alphas_x = np.zeros((pars['number_species'],pars['number_species'],pars['number_nodes']))
    for x in range(pars['number_nodes']):
        Alphas_x[:,:,x] = Random_Alphas(aver_alpha = aver_alphas, sigma = pars['sigma_param'])
        # If the intraspecific competition at a given location were negative, because of the spatial heterogeneity stochasticity, we correct it to be just 0. Otherwise, the population of such species will go to infinity.
        for sp in range(pars['number_species']):
            Alphas_x[sp,sp,x] = np.max([Alphas_x[sp,sp,x], 0.])
        
    
    #flying = flying_pas & flying_act
    
    # In the data traits, regarding the species size, we have body_sizes categories of
    # ≤ .25 cm  Size1; > .25-.5 cm Size2; > .5-1 cm Size3; > 1-2 cm Size4; > 2-4 cm Size5; > 4-8 cm Size6; > 8 cm Size7
    catefories          = [1,      2,     3,    4,   5,   6,   7 ]
    sizes_categories_cm = [0.1875, 0.375, 0.75, 1.5, 3.0, 6.0, 12.0]
    
    
    data_size= traits_comp[ [i for i in traits_comp.columns[traits_comp.columns.str.startswith('Size')]] ]
    aver_size_category = sum([(i+1) * data_size.iloc[:,i] for i in range(data_size.shape[1])]) / sum([ data_size.iloc[:,i] for i in range(data_size.shape[1])])
    aver_size_cm       = np.interp(aver_size_category, catefories, sizes_categories_cm)
    
    if (pars['MOA']=='narcosis'):
        MSS_traits = data_traits['MSS_traits_narcosis']
        EC_traits  = data_traits['EC_traits_narcosis']
    elif (pars['MOA']=='AChE_inhibition'):
        MSS_traits = data_traits['MSS_traits_AChE_inhibition']
        EC_traits  = data_traits['EC_traits_AChE_inhibition']
    elif (pars['MOA']=='atz'):
        MSS_traits = data_traits['MSS_atz']
        EC_traits  = data_traits['EC_atz']
    elif (pars['MOA']=='Cu'):
        MSS_traits = data_traits['MSS_Cu']
        EC_traits  = data_traits['EC_Cu']
    elif (pars['MOA']=='Imida'):
        MSS_traits = data_traits['MSS_Imida']
        EC_traits  = data_traits['EC_Imida']


        
    
    
    k_chem = np.ones(pars['number_species'])
    for sp in range(pars['number_species']):
        if BoolPredator[sp]:
            k_chem[sp] = 1 #np.abs(np.mean(rs[sp,:]))/np.array(EC_traits)[sp]
        else:
            k_chem[sp] = -1 #np.mean(rs[sp,:])/np.array(EC_traits)[sp] #Now, at C=EC, r_prey = r_prey-r_prey = 0; r_pred = -|r_pred| - |r_pred| = -2 * |r_predator|. #0.5 * np.mean(rs[sp,:])/np.array(data_traits.EC_traits)[sp]
    
    
    
    pars_dyn = {
        'rs': rs, 'Alphas': Alphas_pristine, 'Alphas_x': Alphas_x,
        'aquatic_pas': aquatic_pas, 'movs_aq_pas': movs_aq_pas,
        'aquatic_act': aquatic_act, 'movs_aq_act': movs_aq_act,
        'flying_pas': flying_pas, 'movs_fly_pas': movs_fly_pas,
        'flying_act': flying_act, 'movs_fly_act': movs_fly_act,
        'mean_size_species': aver_size_cm, 'predator': BoolPredator,
        'MSS': MSS_traits,
        'EC': EC_traits,
        'k_chem': k_chem}
    
    return pars_dyn


def main(
        TA = '18',
        
        mode_network = 'binary_tree', 
        size_spatial_nw = 32, #64
        #thershold_drainage = 20,
        approx_number_nodes = 25, #50
        esti_edge = 2.0, 
        time_sim=500, 
        
        aver_pred_prey = 0.5, #before 0.5.
        thershold_overlap_foods = 0.5,
        #aver_pred_predator = -1,
        eff_consumption = 0.8, #before 0.5. efficiency of converting dead prey individuals into new predator individuals
        connectance_food_web = 0.05, # connectance of the food web. The number of trophic links scales as L = C * N^2, with N the total number of species of the ecosystem. (ASK CAMILLE)
        preference_1_prey = 0.5, # if one predator have multiple prey species, one of them will be its prefered. (Eklof & Ebenman). If =1, solely predates the prefered species.
        
        aver_disp_rate_aq_pas = 0.1,
        aver_disp_rate_aq_act = 0.2,
        aver_disp_rate_fly_pas = 0.1,
        aver_disp_rate_fly_act = 0.2,
        
        ratio_disp_fly_current_orientation = 0.3,
        
        aver_dist_nodes = 1.0,
        aver_growth_rate = 1.0,
        sigma_env = 1.0,
        sigma_param = 0.1, #Before: 0.3. From the db, we compute average values of the parameters (eg., growth rates). Then, for each node is taken a random value around this mean, with standard deviation sigma_param
        frac_polluted_nodes = 1., #0: No pollution
        conc_polluted_nodes = 3.,
        loc_polluted_nodes = 'random',
        MOA = 'AChE_inhibition',
        limit_detec = 1E-10,
        omega_size = 0.5,
        ratio_size_opt = 4.,
        increase_growth_prefer_hab = 1,
        decrease_growth_not_prefer_hab = 3,
        vel_river = 0.5,
        method = "GLOBI_Traits"):
    
    file_traits = os.path.join(os.getcwd(),'getreal-model','data_GETREAL','newTAs',
      ''.join(['Tachet_NewTA_', TA, '_210317_chems_v4_b.csv']))
    #pristine: chem_conc = 0.
    
    data_traits = load_database(file_traits)
    data_traits.fillna(0) # Change nan values on the database to 0's
    data_traits.set_index('Taxa')
    ## CHOOSING JUST A SUBSECT OF ALL THE DATA SET. WE WILL WORK WITH JUST A SUBSJECT OF ALL THE SPECIES
    # data_traits = data_traits.sample(n=number_species) # I just work with a "number_species" species, select randomly from the dataset
    data_traits = data_traits.reset_index(drop=True) # reset the index values to have 0,1,2,3... (needed when trying to pick values of the traits for the species)
    number_species = data_traits.shape[0]
    
    
    
    pars = {
        'TA': TA, 
        'mode_network': mode_network,
        'number_species': number_species,
        'esti_edge':  esti_edge, 
        'time_sim': time_sim, 
        'aver_pred_prey': aver_pred_prey,
        'thershold_overlap_foods': thershold_overlap_foods,
        #'aver_pred_predator': aver_pred_predator,
        
        'eff_consumption': eff_consumption,
        'connectance_food_web': connectance_food_web, 
        'preference_1_prey': preference_1_prey, 
       
        
        'aver_disp_rate_aq_pas': aver_disp_rate_aq_pas,
        'aver_disp_rate_aq_act': aver_disp_rate_aq_act,
        'aver_disp_rate_fly_pas': aver_disp_rate_fly_pas,
        'aver_disp_rate_fly_act': aver_disp_rate_fly_act,
        
        'increase_growth_prefer_hab': increase_growth_prefer_hab,
        'decrease_growth_not_prefer_hab': decrease_growth_not_prefer_hab,
        
        'ratio_disp_fly_current_orientation': ratio_disp_fly_current_orientation,
        
        'ratio_size_opt': ratio_size_opt,
    
        'aver_dist_nodes': aver_dist_nodes,
        'aver_growth_rate': aver_growth_rate,
        'sigma_env': sigma_env,
        'sigma_param': sigma_param, #Before: 0.3. From the db, we compute average values of the parameters (eg., growth rates). Then, for each node is taken a random value around this mean, with standard deviation sigma_param
        'frac_polluted_nodes': frac_polluted_nodes,
        'conc_polluted_nodes': conc_polluted_nodes,
        'MOA': MOA,
        'limit_detec': limit_detec,
        'omega_size': omega_size, #reduction on the alpha parameters by the size
        'size': size_spatial_nw,
        #'thrA': thershold_drainage,
        'approx_number_nodes': approx_number_nodes,
        'vel_river': vel_river
        }
    
    pars_chemres_size1 = {'EC50': 1., 'k': 1., 's': 0.1}
     
    # print('Simulations for:')
    # print('   ', 'mode_network =', pars['mode_network'])
    # if(pars['mode_network']=='random'):
    #     print('   ', 'esti_edge =', pars['esti_edge'])
    # print('   ', 'average distance between connected patches =', pars['aver_dist_nodes'])
    # print('   ', 'species =', pars['number_species'])
    # print('   ', 'growth_rate =', pars['aver_growth_rate'])
    # print('   ', 'dispersal rate of passive aquatic =', pars['aver_disp_rate_aq_pas'])
    # print('   ', 'dispersal rate of active active =', pars['aver_disp_rate_aq_act'])
    # print('   ', 'dispersal rate of passive flying =', pars['aver_disp_rate_fly_pas']),
    # print('   ', 'dispersal rate of active flying =', pars['aver_disp_rate_fly_act'])
    # print('   ', 'environmental standard deviation =', pars['sigma_env'])
    # print('   ', 'parameters sd for values at different nodes =', pars['sigma_param'])
    
    ## Generate the spatial connected network
    pars_nw = GenerateSpatialNetwork_OCNet_MaxReach_presets(pars = pars)
    number_nodes = len(pars_nw['node'])
    pars['number_nodes'] = number_nodes
    dist_to_0 = np.zeros(number_nodes)
    for k1 in range(number_nodes):
        dist_to_0[k1] = nx.shortest_path_length(pars_nw['G'],0,k1)
    pars_nw['dist_to_0'] = dist_to_0
    pars_nw['loc_polluted_nodes'] = loc_polluted_nodes
    
# =============================================================================
#     ## Actually, I may not want to have a random distribution of chemical patches.
#     conc_polluted_nodes_prev = conc_polluted_nodes
#     conc_polluted_nodes = np.random.uniform(0,conc_polluted_nodes, pars['number_nodes']) 
#     if loc_polluted_nodes=='upstream':
#         temp = np.zeros(pars['number_nodes'])
#         temp[range(round( (1-frac_polluted_nodes)*number_nodes ), number_nodes)] = 1
#         pars_nw['bool_polluted_node'] = temp        
#     elif loc_polluted_nodes=='downstream':
#         temp = np.zeros(pars['number_nodes'])
#         temp[range(round( frac_polluted_nodes*number_nodes ))] = 1
#         pars_nw['bool_polluted_node'] = temp
#     
#     conc_polluted_nodes = pars_nw['bool_polluted_node'] * conc_polluted_nodes
#     pars_nw['chem_conc'] = conc_polluted_nodes
#     pars['frac_polluted_nodes'] = np.sum(pars_nw['bool_polluted_node'] )/number_nodes    
# =============================================================================   
    
    ## Compute from the traits the and the chemical concnetration the model 
    ## parameters for all the species
    pars_dyn = SpeciesParameters_Troph_v3_method(pars = pars, pars_nw = pars_nw, 
                                 pars_chemres_size1 = pars_chemres_size1,
                                 data_traits = data_traits, method = method)
    
    
    ## Connected Network with some of the nodes polluted
    output = ComputeDynamics(pars = pars, pars_nw = pars_nw, pars_dyn = pars_dyn)
    output['ShannonBiodiv'] = ShannonBiodiv(output['final_pops_corr']) 
    
    ## Disconnected Network with some of the nodes polluted
    pars_nw_disc = pars_nw.copy()
    pars_nw_disc['edges_norm_dist'] = np.zeros(np.shape(pars_nw['edges_norm_dist']))
    pars_nw_disc['edges'] = np.zeros(np.shape(pars_nw['edges']))
    pars_nw_disc['connectivity'] = np.zeros(np.shape(pars_nw['connectivity']))
    
    pars_dyn_disc = SpeciesParameters_Troph_v3_method(pars = pars, pars_nw = pars_nw_disc, 
                             pars_chemres_size1 = pars_chemres_size1,
                             data_traits = data_traits, method = method)
    pars_dyn_disc['movs_aq_pas']  = 0. * pars_dyn_disc['movs_aq_pas']
    pars_dyn_disc['movs_aq_act']  = 0. * pars_dyn_disc['movs_aq_act']
    pars_dyn_disc['movs_fly_pas'] = 0. * pars_dyn_disc['movs_fly_pas']
    pars_dyn_disc['movs_fly_act'] = 0. * pars_dyn_disc['movs_fly_act']
    
    output_disc = ComputeDynamics(pars = pars, pars_nw = pars_nw_disc, pars_dyn = pars_dyn_disc)
    output_disc['ShannonBiodiv'] = ShannonBiodiv(output_disc['final_pops_corr'])
    
    
    ## SAVING RESULTS
    sufix = '_'.join([
        'MOA', MOA,
        'Nodes', str(approx_number_nodes)#,
        # 'NSpecies', str(number_species),
        # 'Mode',  mode_network, 
        # 'Loc', loc_polluted_nodes,
        # 'MaxChem', str(conc_polluted_nodes)#,
        # 'DispAqPas', str(aver_disp_rate_aq_pas),
        # 'DispAqAct', str(aver_disp_rate_aq_act),
        # 'DispFlyPas', str(aver_disp_rate_fly_pas),
        # 'DispFlyAct', str(aver_disp_rate_fly_act)
        ])
    sufix = ''.join([sufix, '_v4_NewTomInput_corrected_v3'])
    if (method == "EatAll_Traits"):
        sufix = ''.join([sufix, '_eatall_traits'])
    elif (method == "EatAll_NOTRAITS"): 
        sufix = ''.join([sufix, '_eatall_NOTRAITS'])
    elif (method == "GLOBI_Traits"): 
        sufix = ''.join([sufix, '_globi_traits'])  
    sufix = ''.join([sufix, '.csv'])
    
    idx_sim = RandString(k=6)
    date = datetime.datetime.now()
    idx_sim = ''.join([
        date.strftime("%Y.%m.%d-%H:%M:%S"), '-',idx_sim
        ])
    
    file_save_Net  = ''.join(['NewTA_', pars['TA'], '_v17_PolVelRiv_ChangeDispChem05_new_NET_OCNet_', sufix])
    file_save_summ = ''.join(['NewTA_', pars['TA'], '_v17_PolVelRiv_ChangeDispChem05_new_SUMM_OCNet_', sufix])
    
    file_save_Net = os.path.join( os.getcwd(), 'getreal-model', 'results', 'newTAs', file_save_Net) #os.path.join( os.path.dirname(os.getcwd()), 'results', file_save_Net)
    file_save_summ = os.path.join( os.getcwd(), 'getreal-model', 'results', 'newTAs', file_save_summ)
    
    ## Counting the total regional biodiversity
    output['biodiv_net'],      output_disc['biodiv_net'] = 0, 0
    output['biodiv_net_corr'], output_disc['biodiv_net_corr'] = 0, 0
    for i in range(pars['number_species']):
        # Connected and polluted
        if (output['final_pops'][i,:]>limit_detec).any():
            output['biodiv_net'] = output['biodiv_net'] + 1
        if (output['final_pops_corr'][i,:]>limit_detec).any():
            output['biodiv_net_corr'] = output['biodiv_net_corr'] + 1
        # Disconnected and polluted
        if (output_disc['final_pops'][i,:]>limit_detec).any():
            output_disc['biodiv_net'] = output_disc['biodiv_net'] + 1
        if (output_disc['final_pops_corr'][i,:]>limit_detec).any():
            output_disc['biodiv_net_corr'] = output_disc['biodiv_net_corr'] + 1
    
       
    ## Indicators of local and regional biodiversity
    output['gammadiv_corr']           = np.repeat(output['biodiv_net_corr'],           pars['number_nodes'])
    output_disc['gammadiv_corr']      = np.repeat(output_disc['biodiv_net_corr'],      pars['number_nodes'])
    
    output['alphadiv_corr']           = np.sum(output['final_pops_corr']           >limit_detec, axis=0)
    output_disc['alphadiv_corr']      = np.sum(output_disc['final_pops_corr']      >limit_detec, axis=0)
    
    output['betadiv_corr']           = output['gammadiv_corr']           - output['alphadiv_corr']
    output_disc['betadiv_corr']      = output_disc['gammadiv_corr']      - output_disc['alphadiv_corr']
    
         
                    
    data_net = np.array([
      idx_sim, file_traits, MOA, conc_polluted_nodes, size_spatial_nw, number_nodes, number_species,
      aver_pred_prey, omega_size, ratio_size_opt,thershold_overlap_foods,
      sigma_param,
      increase_growth_prefer_hab, decrease_growth_not_prefer_hab,
      vel_river,
      aver_disp_rate_aq_pas, aver_disp_rate_aq_act,
      # Connected NW with some nodes polluted
      np.mean(output['frac_flying_corr']),
      np.mean(output['frac_species_pred_corr']),     
      np.mean(output['ShannonBiodiv']),
      np.mean(output['alphadiv_corr']),
      np.mean(output['betadiv_corr']),
      np.mean(output['gammadiv_corr']), 
      # Disconnected NW with some nodes polluted
      np.mean(output_disc['frac_flying_corr']),
      np.mean(output_disc['frac_species_pred_corr']),     
      np.mean(output_disc['ShannonBiodiv']),
      np.mean(output_disc['alphadiv_corr']),
      np.mean(output_disc['betadiv_corr']),
      np.mean(output_disc['gammadiv_corr']) 
    ])
    
    data_net = np.reshape(  data_net, (1, len(data_net))  )
    
    if is_non_zero_file(file_save_Net)==False:
        with open(file_save_Net, mode='a+', newline='') as csv_file:
            csv_file.write(
                ';'.join(
                    np.hstack([
                        'idx',  'file_traits', 'MOA', 'MaxChem', 'size_spatial_nw', 'number_nodes','number_species',
                        'aver_pred_prey', 'omega_size', 'ratio_size_opt','thershold_overlap_foods',
                        'sigma_param',
                        'increase_growth_prefer_hab', 'decrease_growth_not_prefer_hab',
                        'vel_river',
                        'aver_disp_rate_pas', 'aver_disp_rate_act',
                        # Connected NW with some nodes polluted
                        'FracFlying_conn','FracPredators_conn',
                        'Shannon_conn','Alphadiv_conn','Betadiv_conn','Gammadiv_conn',
                        # Disconnected NW with some nodes polluted
                        'FracFlying_disc','FracPredators_disc',
                        'Shannon_disc','Alphadiv_disc','Betadiv_disc','Gammadiv_disc',
                    ])
                )
            )
            csv_file.write('\n')
    pandas.DataFrame(data_net).to_csv(file_save_Net, 
                                      header = None, index = None, 
                                      mode = 'a', sep = ';')       
    
    taxa = np.array(data_traits['Taxa'])
                                                            
    taxa_surv = taxa[np.sum(output['final_pops_corr'], axis=1)>0]
    taxa_surv_disc = taxa[np.sum(output_disc['final_pops_corr'], axis=1)>0]
    
    
    data_summ = np.concatenate([
        np.array([
            idx_sim, file_traits, MOA, conc_polluted_nodes, size_spatial_nw, number_nodes,
              number_species, aver_pred_prey, omega_size, ratio_size_opt, 
              thershold_overlap_foods, sigma_param, 
              increase_growth_prefer_hab, decrease_growth_not_prefer_hab,
              vel_river,
              aver_disp_rate_aq_pas, aver_disp_rate_aq_act,
              aver_disp_rate_fly_pas, aver_disp_rate_fly_act
            ]),
        np.sum(output['final_pops_corr'], axis=1),
        np.sum(output_disc['final_pops_corr'], axis=1),
        np.array([
            '-'.join(taxa_surv),
            '-'.join(taxa_surv_disc)])
        ])
    data_summ = np.reshape(  data_summ, (1, len(data_summ))  )
    
    if is_non_zero_file(file_save_summ)==False:
        with open(file_save_summ, mode='a+', newline='') as csv_file:
            csv_file.write(
                ';'.join(
                    np.hstack([
                        'idx',  'file_traits', 'MOA', 'MaxChem', 'size_spatial_nw', 'number_nodes','number_species',
                        'aver_pred_prey', 'omega_size', 'ratio_size_opt','thershold_overlap_foods',
                        'sigma_param',
                        'increase_growth_prefer_hab', 'decrease_growth_not_prefer_hab',
                        'vel_river',
                        'disp_rate_aq_pas', 'disp_rate_aq_act',
                        'disp_rate_fly_pas', 'disp_rate_fly_act',
                        # Connected NW with some nodes polluted
                        ['_'.join([taxa[i],'conn']) for i in range(number_species)],
                        # Disconnected NW with some nodes polluted
                        ['_'.join([taxa[i],'disc']) for i in range(number_species)],
                        'TaxaRemainConn','TaxaRemainDisc'
                    ])
                )
            )
            csv_file.write('\n')
    pandas.DataFrame(data_summ).to_csv(file_save_summ, 
                                      header = None, index = None, 
                                      mode = 'a', sep = ';')    
    
    

    
    
    return 0  


start = time.perf_counter()
tracemalloc.start()


ratio_size_opt = 4
frac_polluted_nodes = 1
increase_growth_prefer_hab = 0.
decrease_growth_not_prefer_hab = 0.
thershold_overlap_foods = 0.5

sigma_param =  float(sys.argv[1])
TA          =  str(sys.argv[2])
method      =  str(sys.argv[3])


AVER_PRED_PREY = [1.]
OMEGA_SIZE =     [.5]
DISP_RATE =      [.1]

VELS_RIVER =     [ 0]
MOAS = ['atz', 'Cu', 'Imida']#['narcosis', 'AChE_inhibition'] 
#CHANGE_HABITAT_GROWTH = [0]
#INCREASE_GROWTH_PREFER = [0, 1, 1]
#DECREASE_GROWTH_PREFER = [0, 1, 2]
for cont in range(5):
    for disp_rate in DISP_RATE:
        for aver_pred_prey in AVER_PRED_PREY:
            for omega_size in OMEGA_SIZE:
                for vel_river in VELS_RIVER:
                    for moa in MOAS:
                        if moa == 'atz':
                            chem_ini = 1
                            chem_end = 6
                        elif moa == 'Cu':
                            chem_ini = 0
                            chem_end = 5
                        elif moa == 'Imida':
                            chem_ini = -4.5
                            chem_end = 2
                        CONCS = np.append([0], signif(10**np.arange(chem_ini,chem_end+0.1,0.1), 4))
                        for conc_polluted_nodes in CONCS:
                            print('-'.join([str(disp_rate), str(aver_pred_prey), str(omega_size), str(conc_polluted_nodes),str(cont)]))
                            main(TA = TA, approx_number_nodes = 100, time_sim = 1000,
                                 aver_pred_prey = aver_pred_prey,
                                 MOA = moa,
                                 conc_polluted_nodes = conc_polluted_nodes,
                                 thershold_overlap_foods = thershold_overlap_foods,
                                 aver_disp_rate_aq_pas = disp_rate, aver_disp_rate_aq_act = disp_rate,
                                 aver_disp_rate_fly_pas = disp_rate, aver_disp_rate_fly_act = disp_rate,
                                 sigma_param = sigma_param, ratio_size_opt = ratio_size_opt, 
                                 increase_growth_prefer_hab = increase_growth_prefer_hab,
                                 decrease_growth_not_prefer_hab = decrease_growth_not_prefer_hab,
                                 frac_polluted_nodes = frac_polluted_nodes,
                                 omega_size = omega_size,
                                 vel_river = vel_river, 
                                 method = method)


current, peak = tracemalloc.get_traced_memory()
print(f'Current memory usage is {current / 10**6}MB; Peak was {peak / 10**6}MB')
tracemalloc.stop()

end = time.perf_counter()
print("Elapsed time during the whole program in seconds:", end-start) 