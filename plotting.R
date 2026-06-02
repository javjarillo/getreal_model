################################################################################
#                                                                              #
# R script to reproduce selected figures from the article                      #
# "Predictors of food web resistance to environmental change"                  #
# (currently under review).                                                    #
#                                                                              #
# To reproduce these plots, you need the simulation output files available at  # 
# https://doi.org/10.5281/zenodo.19110901                                      #
# (which were generated using the code provided in this repository)            #
#                                                                              #
################################################################################





library(tidyverse)
library(scales) # to access break formatting functions
library(glue)
library(matlib)
library(styler)


####################### FIGURE 2 #######################################
######################################################################


dg <- read.csv("results_random/extinction_ln_histograms_v3_max_min.csv",
               sep = ",")

dg <- dg %>% mutate(
  sigma_str = factor(
    sigma,
    levels = c(0.05, 0.10, 0.25),
    labels = c(
      "\u03C3 = 5%",
      "\u03C3 = 10%",
      "\u03C3 = 25%"
    )
  ),
  mode_predator_text = factor(
    mode_predators,
    levels = c("generalists", "specialists"),
    labels = c('Generalist Predators', 'Specialist Predators')
  )
)






dg <- dg %>% 
  mutate(k = 3.0) %>% 
  mutate(
    iOmega_prey_estimate = 1 / (b / m - 1) * (1 + sigma^2) / tauV * (
      1 + k * (1 + (nV_max - 1) * alpha) / (1 - alpha) * sigma),
    iOmega_predator_estimate = (
      (1 + sigma^2) / tauV +
        d / m * (1 + (nV_max - 1) * alpha) / lambda / eta * (1 + sigma^2) / tauP
      ) / (
        b / m - 1 - d / m * (1 + (nV_max - 1) * alpha) / lambda / eta
      ) * (
        1 + k * (
          (d / m * (1 + (nV_max - 1) * alpha) / lambda / eta +
             lambda * nP_max / nV_max) / (
               (1 + sigma^2) / tauV + d / m * (1 + (nV_max - 1) * alpha) / lambda / eta
             ) *
            (1 + sigma^2) / tauP * (1 + sigma^2) / tauP * sigma
        )
      ),
    HP = 1/mean_tauPM1,
    HV = 1/mean_tauVM1,
    fract_taus = HP/HV
  )








dF2_predator <- dg %>% group_by(
  mode_predators, mode_predator_text, sigma, sigma_str, idx, nV_max, nP_max, b, m, d,
  alpha, lambda, eta, tauV, tauP, 
  iOmega_prey_estimate,
  iOmega_predator_estimate
) %>% filter(nP < nP_max) %>% 
  summarise(stressor_th = min(stressor, na.rm = T), species = 'predator', fract_taus = mean(fract_taus, na.rm = TRUE))

dF2_prey <- dg %>% group_by(
  mode_predators, mode_predator_text, sigma, sigma_str, idx, nV_max, nP_max, b, m, d,
  alpha, lambda, eta, tauV, tauP, 
  iOmega_prey_estimate,
  iOmega_predator_estimate
) %>% filter(nV < nV_max) %>% 
  summarise(stressor_th = min(stressor, na.rm = T), species = 'prey', fract_taus = mean(fract_taus, na.rm = TRUE))

dF2_1stext <- rbind(dF2_prey, dF2_predator) %>% mutate(
  iOmega_estimate = ifelse(species == 'prey', iOmega_prey_estimate, iOmega_predator_estimate)
)


dF2_1stext_new <- rbind(
  dF2_1stext %>% mutate(method = 'simulations', valor = stressor_th),
  dF2_1stext %>% mutate(method = 'From tau distr.', valor = 1/iOmega_estimate)
) 


rm(dF2_prey, dF2_predator)


ggplot(data = dg %>% ungroup()  %>% 
         filter((sigma == 0.05 & stressor >= 1 & stressor <=300) | 
                  (sigma == 0.1 & stressor >= 1 & stressor <= 300) |
                  (sigma == 0.25 & stressor >= 1 & stressor <= 300) |
                  (sigma == 0.5 & stressor >= 1 & stressor <= 300)
         ),
       mapping = aes(x=stressor, group = interaction(sigma, mode_predators, idx))) +
  facet_grid(factor(mode_predator_text)~factor(sigma_str), scales = 'free') +
  scale_x_continuous(trans='log10', limits = c(2,200)) +
  scale_y_continuous(breaks=c(0,2,4,6,8,10), limits = c(0,12) ) + 
  geom_vline(
    mapping = aes(xintercept=1/iOmega_prey_estimate, group = interaction(sigma, mode_predators)), 
    color ="#4477AA", size=2, alpha = 0.5) +
  geom_vline(
    mapping = aes(xintercept=1/iOmega_predator_estimate, group = interaction(sigma, mode_predators)),
    color ="#EE6677", size=2, alpha = 0.5) +
  geom_line(aes(y=nV, color = 'prey'), alpha = I(0.2)) +
  geom_line(aes(y=nP, color = 'predator'), alpha = I(0.2)) +
  geom_boxplot(data=dF2_1stext_new %>% ungroup() %>% filter(species=='predator') %>% filter(method == "simulations"), 
               mapping = aes(x=valor, y = 6, group = interaction(sigma_str, mode_predators)), 
               color="darkred", fill = "white", notch=FALSE, outlier.size=1) + 
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  geom_boxplot(data=dF2_1stext_new %>% filter(species=='prey'), 
               mapping = aes(x=stressor_th, y = 11, group = interaction(sigma_str, mode_predators)), 
               color="darkblue", fill = "white", notch=FALSE, outlier.size=1) + 
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  scale_colour_manual(values = c("darkblue", "darkred"), breaks = c("prey", "predator")) + 
  labs(x = 'environmental change', y = 'remaining species',
       colour = 'trophic level') +
  theme_bw()


########## FIGURE 3 A #############################
##################################################

dg3a <- read.csv('results_random/extinction_ln_histograms_v4_simplified_newvals_min_max_fixed.csv')



dg3a_summ_prey <- dg3a %>% group_by(idx, nV_max, nP_max, mode_predators) %>% 
  filter(nV < nV_max) %>% summarise(resist = min(stressor), species = 'prey') %>% 
  select(idx, nV_max, nP_max, mode_predators, resist, species)
dg3a_summ_predator <- dg3a %>% group_by(idx, nV_max, nP_max, mode_predators) %>% 
  filter(nP < nP_max) %>% summarise(resist = min(stressor), species = 'predator') %>% 
  select(idx, nV_max, nP_max, mode_predators, resist, species)

dg3a_summ <- rbind(dg3a_summ_prey, dg3a_summ_predator)

  

library(RColorBrewer)


# levels for nV
levels_nV <- c(2, 4, 8, 16,32) 

# color palette for these leves (YlOrRd)
pal <- brewer.pal(n = length(levels_nV), name = "YlOrRd")

# color for the different levels
fill_colors <- setNames(pal, levels_nV)



ggplot(
  dg3a_summ %>% mutate(nV_max = as.factor(nV_max), nP_max = as.factor(nP_max)) %>% group_by(nV_max, nP_max, mode_predators, species), 
  aes(x = nV_max, y = resist, fill = nP_max, 
      group = interaction(nV_max, nP_max, species, mode_predators))) + 
  facet_grid(species ~ mode_predators, scales = "free") + 
  scale_y_continuous(trans='log10') + #limits = c(1.3,5.2)
  scale_fill_manual(values = fill_colors, drop = TRUE) +  # drop = TRUE evita mostrar niveles vacíos
  scale_color_manual(values = fill_colors, drop=TRUE, guide='none') +
  geom_violin(trim = TRUE, color=NA, alpha=0.7,
              position = position_dodge(width = 0.7), width=0.8, show.legend=F) + 
  geom_boxplot(position = position_dodge(width = 0.7), width=0.5, color='black',
               outlier.size=0.5, outlier.alpha = 0.25) + 
  theme_bw() 




##################### FIGURE 3B ###########################################
##########################################################################

dg3b <- read.csv('results_random/extinction_ln_histograms_v3_max_min_fract_taus_small.csv')

dg3b <- dg3b %>% mutate(
  sigma_str = factor(
    sigma,
    levels = c(0.01, 0.05, 0.10),
    labels = c(
      "\u03C3 = 1%",
      "\u03C3 = 5%",
      "\u03C3 = 10%"
    )
  )
)

dg3b_summ <- dg3b %>% group_by(idx, nV_max, nP_max, fract, sigma_str, mode_predators) %>% 
  filter(nV + nP < nV_max + nP_max) %>% summarise(resist = min(stressor))


ggplot(
  dg3b_summ, aes(group = interaction(fract, mode_predators, sigma_str), x=fract, y = resist, color = sigma_str) ) + 
  facet_grid(~mode_predators) + 
  scale_y_continuous(trans='log10')  + scale_x_continuous(trans='log10') +
  geom_boxplot() + 
  scale_color_brewer(palette ='Dark2', type='qual') + 
  theme_bw()





