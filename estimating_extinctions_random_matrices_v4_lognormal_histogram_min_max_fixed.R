library(tidyverse)
library(scales) # to access break formatting functions
library(glue)
library(matlib)




# CREATING THE DATABASE
path <- file.path(getwd(), "results", "simu_random_com")
if (!dir.exists(path)) {
  dir.create(path, recursive = TRUE)
}
setwd(path)


df <- data.frame()
b = 10; m = 1/10; d = 1/10; alpha = 1/20; 
lambda = 1; eta = 0.5; 
tauM = 1; #tauV=1; tauP = 1;
nV_max = 10; nP_max = 5;

sigma = 0.1

mode_predators = 'specialist'

#############################################################################
# Only needed if you want a trad-off between alpha (competition strength of 
# competitors) and lambda (attack rate), so the stronger competitors are also
# more attacked by the predators
# -------------------------------------------
# increase_alpha_prey_weakcompetitors = 4
# tradeoff_lambda_tau = NaN
###########################################################################

############################################################################
# Specify the minimum and maximum tolerance of prey of predators within the 
# community (so the sample minimum and maximum tolerances do not depend on
# the sample size)
# ------------------------------------------------------------------------
# We put these sample minimum and maximum tolerances as the tolerances at
# 5 standard deviations from the population mean.
# 
# mu_0 =  log(1/sqrt(1+sigma**2))
# sd_0 = sqrt(log(1+sigma**2))
# > exp(mu_0 - 5*sd_0)
# [1] 0.6042714
# > exp(mu_0 + 5*sd_0)
# [1] 1.638501

minimo = 0.5
maximo = 2
############################################################################




results <- list()
k <- 1


for (nV_max in 2^seq(1,5)){
  for (nP_max in 2^seq(1, log2(nV_max))){
    
    for (fract in c(1)){ # 0.1, 1, 10
      tauV = 2*tauM/(1+fract); tauP = tauV*fract;
      for (mode_predators in c('specialists', 'generalists')){
        for (tradeoff_lambda_tau in c(1)){ #c(1,  4,  16)
          for (sigma in c(0.1)){ # c(0.01, 0.05, 0.1) #  #c(0, 0.10, 0.25, 0.50, 1., 2., 10.)
            
            
            log_sd <- sqrt(log(1 + sigma^2))
            log_alpha <- log(alpha / sqrt(1 + sigma^2))
            log_lambda <- log(lambda / sqrt(1 + sigma^2))
            log_eta <- log(eta / sqrt(1 + sigma^2))
            log_b <- log(b / sqrt(1 + sigma^2))
            log_m <- log(m / sqrt(1 + sigma^2))
            log_d <- log(d / sqrt(1 + sigma^2))
            log_tauV <- log(tauV / sqrt(1 + sigma^2))
            log_tauP <- log(tauP / sqrt(1 + sigma^2))
            
            if (sigma == 0) {idx_max = 2} else {idx_max = 2} #500
            
            for(idx in 1:idx_max){
              print(c(nV_max, nP_max, fract, mode_predators, tradeoff_lambda_tau, sigma, idx))
              nV = nV_max; nP = nP_max;
              
              n_eq <- matrix(rep(-1, nV+nP), nV+nP, 1)
              while (sum(n_eq < 0)){
                # Matrix A interactions
                #A_VV <- abs(matrix(rnorm(nV*nV, mean = alpha, sd = sigma*alpha), nV, nV)) # If we want N(mu, sigma)
                A_VV <- exp(matrix(rnorm(nV*nV, mean =log_alpha, sd = log_sd), nV, nV)) # If we want LN(mu, sigma)
                diag(A_VV) <- rep(1, nV)
                
                # # trade-off: first np species: strong competitors, actually eaten by predators
                # # nV-nP other species. weak competitors.
                # 
                # A_VV[1:nP, (nP+1):nV] <- exp(matrix(rnorm(nP*(nV-nP), 
                #                                           mean = log(alpha*increase_alpha_prey_weakcompetitors/sqrt(1+sigma**2)), 
                #                                           sd = sqrt(log(1+sigma**2))), nP, nV-nP))
                # A_VV[(nP+1):nV, 1:nP] <- exp(matrix(rnorm(nP*(nV-nP), 
                #                                           mean = log(alpha*(1/increase_alpha_prey_weakcompetitors)/sqrt(1+sigma**2)), 
                #                                           sd = sqrt(log(1+sigma**2))), nV-nP, nP))
                # 
                
                alphas_vector <- A_VV[row(A_VV) != col(A_VV)]
                
                mean_alpha   = mean(alphas_vector);   var_alpha    = var(alphas_vector)
                mean_alphaM1 = mean(1/alphas_vector); var_alphaM1  = var(1/alphas_vector)
                min_alpha = min(1/alphas_vector);     max_alpha = max(1/alphas_vector)
                mean_log10_alpha = mean(log10(alphas_vector)); var_log10_alpha = var(log10(alphas_vector))
                
                A_PP <- diag(nP)
                if (mode_predators == 'generalists'){
                  #A_VP <- (1./nV)*abs(matrix(rnorm(nV*nP, mean = lambda, sd = sigma*lambda), nV, nP)) # If we want N(mu, sigma)
                  A_VP <- (1./nV)*exp(matrix(rnorm(nV*nP, mean = log_lambda, sd = log_sd), nV, nP)) # If we want LN(mu, sigma)
                  
                  lambdas_vector <- c(nV*A_VP);
                  
                  mean_lambda   = mean(lambdas_vector); var_lambda= var(lambdas_vector)
                  min_lambda    = min(lambdas_vector);  max_lambda    = max(lambdas_vector)
                  mean_lambdaM1 = mean(1/lambdas_vector); var_lambdaM1  = var(1/lambdas_vector);
                  mean_log10_lambda = mean(log10(lambdas_vector)); var_log10_lambda = var(log10(lambdas_vector))
                  

                } else if (mode_predators == 'specialists') {
                  #A_VP <- diag(x=1, nrow = nV, ncol = nP) * abs(matrix(rnorm(nV*nP, mean = lambda, sd = sigma*lambda), nV, nP)) # If we want N(mu, sigma)
                  A_VP <- diag(x=1, nrow = nV, ncol = nP) * exp(matrix(rnorm(nV*nP, mean = log_lambda, sd = log_sd), nV, nP)) # If we want LN(mu, sigma)
                  
                  lambdas_vector <- c(diag(A_VP));
                  
                  mean_lambda   = mean(lambdas_vector); var_lambda= var(lambdas_vector)
                  min_lambda    = min(lambdas_vector);  max_lambda    = max(lambdas_vector)
                  mean_lambdaM1 = mean(1/lambdas_vector); var_lambdaM1  = var(1/lambdas_vector);
                  mean_log10_lambda = mean(log10(lambdas_vector)); var_log10_lambda = var(log10(lambdas_vector))
                  
                }
                
                #A_PV <- -1. * t(A_VP) * abs(matrix(rnorm(nV*nP, mean = eta, sd = sigma*eta), nP, nV)) # If we want N(mu, sigma) 
                eta_PV <- exp(matrix(rnorm(nV*nP, mean = log(eta/sqrt(1+sigma**2)), sd = sqrt(log(1+sigma**2))), nP, nV)) # If we want LN(mu, sigma)
                A_PV <- -1. * t(A_VP) * eta_PV
                
                A_total <- rbind(cbind(A_VV, A_VP), cbind(A_PV, A_PP)) 
                
                if (mode_predators == 'generalists'){
                  etas_vector = c(eta_PV)
                } else if (mode_predators == 'specialists') {
                  etas_vector = c(diag(eta_PV))
                }
                
                mean_eta   = mean(etas_vector);   var_eta    = var(etas_vector)
                min_eta    = min(etas_vector);   max_eta    = max(etas_vector)
                mean_etaM1 = mean(1/etas_vector); var_etaM1  = var(1/etas_vector)
                mean_log10_eta   = mean(log10(etas_vector)); var_log10_eta    = var(log10(etas_vector))
                

                
                
                # Growth and mortality rates
                # If we want N(mu, sigma)
                #bS <- abs(rnorm(nV, mean = b, sd = b*sigma))
                #mS <- abs(rnorm(nV, mean = m, sd = m*sigma))
                #dS <- abs(rnorm(nP, mean = d, sd = d*sigma))
                # If we want LN(mu, sigma)
                bS <- exp(rnorm(nV, mean = log_b, sd = log_sd))
                mS <- exp(rnorm(nV, mean = log_m, sd = log_sd))
                dS <- exp(rnorm(nP, mean = log_d, sd = log_sd))
                
                mean_b = mean(c(bS)); var_b = var(c(bS)); 
                mean_bM1 = mean(1/c(bS)); var_bM1 = var(1/c(bS)); 
                min_b = min(c(bS)); max_b = max(c(bS));
                mean_log10_b = mean(log10(c(bS))); var_log10_b = var(log10(c(bS))); 
                
                mean_m = mean(c(mS)); var_m = var(c(mS)); 
                mean_mM1 = mean(1/c(mS)); var_mM1 = var(1/c(mS)); 
                min_m = min(c(mS)); max_m = max(c(mS));
                mean_log10_m = mean(log10(c(mS))); var_log10_m = var(log10(c(mS)));
                
                mean_d = mean(c(dS)); var_d = var(c(dS)); 
                mean_dM1 = mean(1/c(dS)); var_dM1 = var(1/c(dS)); 
                min_d = min(c(dS)); max_d = max(c(dS));
                mean_log10_d = mean(log10(c(dS))); var_log10_d = var(log10(c(dS)));
                
                
                # Sensitivities
                # If we want N(mu, sigma)
                #tauVS <- abs(rnorm(nV, mean = tauV, sd = tauV*sigma)) #rep(tauV, nV)
                #tauPS <- abs(rnorm(nP, mean = tauP, sd = tauP*sigma)) #rep(tauP, nP)
                # If we want LN(mu, sigma)
                tauVS <- exp(rnorm(nV, mean = log_tauV, sd = log_sd)) #rep(tauV, nV)
                tauPS <- exp(rnorm(nP, mean = log_tauP, sd = log_sd)) #rep(tauP, nP)
                
                # We set the sample mimimum and maximum (so they do not depend
                # on the community size)
                tauVS[tauVS < minimo] <- minimo
                tauPS[tauPS < minimo] <- minimo
                
                tauVS[tauVS > maximo] <- maximo
                tauPS[tauPS > maximo] <- maximo
                
                tauVS[which.min(tauVS)] <- minimo
                tauPS[which.min(tauPS)] <- minimo
                
                tauVS[which.max(tauVS)] <- maximo
                tauPS[which.max(tauPS)] <- maximo
                
                
                
                # trade-off: first np species: strong competitors, actually eaten by predators  ->  low tolerances.
                # nV-nP other species. weak competitors. -> high tolerances
                #tauVS[1:nP] <- exp(rnorm(nP, mean = log((tauV/tradeoff_lambda_tau)/sqrt(1+sigma**2)), sd = sqrt(log(1+sigma**2)))) #rep(tauV, nV)
                #tauVS[(nP+1):nV] <- exp(rnorm(nV-nP, mean = log((tauV*tradeoff_lambda_tau)/sqrt(1+sigma**2)), sd = sqrt(log(1+sigma**2))))
                
                
                mean_tauV = mean(c(tauVS)); var_tauV = var(c(tauVS)); 
                mean_tauVM1 = mean(1/c(tauVS)); var_tauVM1 = var(1/c(tauVS)); 
                min_tauV = min(c(tauVS)); max_tauV = max(c(tauVS));
                mean_log10_tauV = mean(log10(c(tauVS))); var_log10_tauV = var(log10(c(tauVS)));
                
                mean_tauP = mean(c(tauPS)); var_tauP = var(c(tauPS)); 
                mean_tauPM1 = mean(1/c(tauPS)); var_tauPM1 = var(1/c(tauPS)); 
                min_tauP = min(c(tauPS)); max_tauP = max(c(tauPS));
                mean_log10_tauP = mean(log10(c(tauPS))); var_log10_tauP = var(log10(c(tauPS)));
                
                
                # Chemical concentration (first, 0 concentration)
                zS <- rep(0, nV+nP)
                
                # Growth rate vector
                R_total <- rbind( 
                  matrix(bS-mS*(1+zS[1:nV]/tauVS), nV, 1),
                  matrix(-dS*(1+zS[(nV+1):(nV+nP)]/tauPS),nP,1)
                )
                
                # Equilibrium densities
                #n_eq <- inv(A_total) %*% R_total
                n_eq <- solve(A_total, R_total)
                
                # Matrix with predator boolean array
                bool_predator = rbind(
                  matrix(rep(0, nV),nV, 1), 
                  matrix(rep(1, nP), nP, 1)
                )
                
                # bool_eaten = rbind(
                #   matrix(rep(1, nP),nP, 1),
                #   matrix(rep(0, nV-nP), nV-nP, 1),
                #   matrix(rep(0, nP), nP, 1)
                #   
                # )
              }
              
              
              results[[k]] <- data.frame(
                idx = idx, stressor = 0, nV = nV, nP = nP,
                nV_max = nV_max, nP_max = nP_max,
                b = b, m = m, d = d, alpha = alpha,
                lambda = lambda, eta = eta,
                tauM = tauM, fract = fract, sigma = sigma,
                mode_predators = mode_predators,
                min_tauV = min_tauV, max_tauV = max_tauV,
                min_tauP = min_tauP, max_tauP = max_tauP
              )
              
              k <- k + 1
              
              
              
              prev_nV <- nV
              prev_nP <- nP
              
              # For the different concentrations, compute the diversity
              zis <- 10^seq(-5, 5, 0.02) # Different concentrations 0.02
              contador = 0
              
              # print(A_total)
              
              for (zi in zis){
                contador = contador + 1
                if (nV > 0 || nP > 0){   # If there are both prey or predator species
                  zS <- rep(zi, length(R_total))
                  
                  if (nP>0){  # If there are predators
                    R_total <- rbind( 
                      matrix(bS-mS*(1+zS[1:nV]/tauVS), nV, 1),
                      matrix(-dS*(1+zS[(nV+1):(nV+nP)]/tauPS),nP,1)
                    )
                  } else { # If there are no predators
                    R_total <- matrix(bS-mS*(1+zS[1:nV]/tauVS), nV, 1)
                  }
                  
                  # Compute equilibrium densities or abundances
                  if (length(A_total)>1){
                    n_eq <- solve(A_total, R_total) #inv(A_total) %*% R_total
                  } else if (length(A_total)==1){
                    n_eq <- 1/(A_total[1]) %*% R_total
                  } else {
                    n_eq = c()
                  }
                  
                  # n_eaten_prey <- sum(bool_eaten == 1)
                  # n_weak_comp <- nV - n_eaten_prey
                  
                  # Count how many species have positive abundances
                  while (sum(n_eq < 0) && nV > 0){
                    if ((nV>0) | (nP>0)){
                      mantener <- n_eq>0
                      mantener_V <- mantener[1:nV]
                      mantener_P <- mantener[(nV+1):(nV+nP)]
                      
                      if (length(n_eq)>1){
                        A_total <- A_total[mantener, mantener]
                      } else {
                        A_total <- A_total[mantener]
                      }
                      
                      
                      R_total <- R_total[mantener]
                      bool_predator <- bool_predator[mantener]
                      # bool_eaten <- bool_eaten[n_eq > 0]
                      bS <- bS[mantener_V]
                      mS <- mS[mantener_V]
                      dS <- dS[mantener_P]
                      
                      zS <- rep(zi, length(R_total))
                      
                      tauVS <- tauVS[mantener_V]
                      tauPS <- tauPS[mantener_P]
                      
                      # tauS <- rbind(matrix(tauVS,nV,1), matrix(tauPS,nP,1))
                      # tauS <- tauS[mantener]
                      
                      nV <- sum(bool_predator==0)
                      nP <- sum(bool_predator==1)
                      # n_eaten_prey <- sum(bool_eaten == 1)
                      # n_weak_comp <- nV - n_eaten_prey
                      
                      if (nV>0){
                        if (nP>0){
                          # tauVS <- tauS[1:nV]
                          # tauPS <- tauS[(nV+1):(nV+nP)]
                          
                          
                          R_total <- rbind( 
                            matrix(bS-mS*(1+zS[1:nV]/tauVS), nV, 1),
                            matrix(-dS*(1+zS[(nV+1):(nV+nP)]/tauPS),nP,1)
                          )
                          
                          if (length(A_total)>1){
                            n_eq <- solve(A_total, R_total)# inv(A_total) %*% R_total
                          } else if (length(A_total)==1){
                            n_eq <- 1/(A_total[1]) %*% R_total
                          } else {
                            n_eq = c()
                          }
                        } else {
                          # tauVS <- tauS[1:nV]
                          R_total <- matrix(bS-mS*(1+zS[1:nV]/tauVS), nV, 1)
                          if (length(A_total)>1){
                            n_eq <- solve(A_total, R_total)# inv(A_total) %*% R_total
                          } else if (length(A_total)==1){
                            n_eq <- 1/(A_total[1]) %*% R_total
                          } else {
                            n_eq = c()
                          }
                        }
                      } else {
                        n_eq <- c()
                      } 
                      
                    }
                  }
                }
                
                if (nV < prev_nV || nP < prev_nP) {
                  
                  results[[k]] <- data.frame(
                    idx = idx, stressor = zi, nV = nV, nP = nP,
                    nV_max = nV_max, nP_max = nP_max,
                    b = b, m = m, d = d, alpha = alpha,
                    lambda = lambda, eta = eta,
                    tauM = tauM, fract = fract, sigma = sigma,
                    mode_predators = mode_predators,
                    min_tauV = min_tauV, max_tauV = max_tauV,
                    min_tauP = min_tauP, max_tauP = max_tauP
                  )
                  
                  k <- k + 1
                  
                }
                prev_nV <- nV
                prev_nP <- nP
                
                
              }
            }
          }       
        }
      }
    }
    
    
  }
}



df <- bind_rows(results[1:(k-1)])

write.csv(
  df,
  "extinction_ln_histograms_NOTRADEOFF_v4_simplified_newvals_min_max_fixed.csv",
  row.names = FALSE
)