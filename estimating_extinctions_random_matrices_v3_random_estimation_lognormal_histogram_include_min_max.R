library(tidyverse)
library(scales) # to access break formatting functions
library(glue)
library(matlib)




# CREATING THE DATABASE
df <- data.frame()
b = 10; m = 1/10; d = 1/10; alpha = 1/20; 
lambda = 1; eta = 0.5; 
tauV=1; tauP = 1;
nV_max = 10; nP_max = 5;

sigma = 0.1

mode_predators = 'generalists'

row <- data.frame(idx = NaN, stressor = NaN, nV = NaN, nP = NaN,
                  nV_max = nV_max, nP_max = nP_max,
                  b = b, m = m, d = d, alpha = alpha,
                  lambda = lambda, eta = eta, 
                  tauV=tauV, tauP = tauP, sigma = sigma,
                  mode_predators = mode_predators, 
                  
                  mean_b = NaN, var_b = NaN, mean_bM1 = NaN,  var_bM1 = NaN, 
                  mean_log10_b = NaN, var_log10_b = NaN, min_b = NaN, max_b = NaN,
                  mean_m = NaN, var_m = NaN, mean_mM1 = NaN,  var_mM1 = NaN, 
                  mean_log10_m = NaN, var_log10_m = NaN, min_m = NaN, max_m = NaN,
                  mean_d = NaN, var_d = NaN, mean_dM1 = NaN,  var_dM1 = NaN, 
                  mean_log10_d = NaN, var_log10_d = NaN, min_d = NaN, max_d = NaN,
                  
                  mean_alpha = NaN, var_alpha = NaN, 
                  mean_alphaM1 = NaN,  var_alphaM1 = NaN,
                  mean_log10_alpha = NaN, var_log10_alpha = NaN, 
                  min_alpha = NaN, max_alpha = NaN,
                  
                  mean_lambda = NaN, var_lambda = NaN, 
                  mean_lambdaM1 = NaN,  var_lambdaM1 = NaN,
                  mean_log10_lambda = NaN, var_log10_lambda = NaN, 
                  min_lambda = NaN, max_lambda = NaN,
                  
                  mean_eta = NaN, var_eta = NaN, 
                  mean_etaM1 = NaN,  var_etaM1 = NaN,
                  mean_log10_eta = NaN, var_log10_eta = NaN, 
                  min_eta = NaN, max_eta = NaN,
                  
                  mean_tauV = NaN, var_tauV = NaN, 
                  mean_tauVM1 = NaN,  var_tauVM1 = NaN,
                  mean_log10_tauV = NaN, var_log10_tauV = NaN, 
                  min_tauV = NaN, max_tauV = NaN,
                  
                  mean_tauP = NaN, var_tauP = NaN, 
                  mean_tauPM1 = NaN,  var_tauPM1 = NaN,
                  mean_log10_tauP = NaN, var_log10_tauP = NaN, 
                  min_tauP = NaN, max_tauP = NaN
                  )


rm(row)

for (mode_predators in c('generalists', 'specialists')){
  for (sigma in c(0.05, 0.10, 0.25)){  #c(0, 0.10, 0.25, 0.50, 1., 2., 10.)
    
    if (sigma == 0) {idx_max = 2} else {idx_max = 500} #500
    
    for(idx in 1:idx_max){
      print(c(mode_predators, sigma, idx))
      nV = nV_max; nP = nP_max;
      
      n_eq <- matrix(rep(-1, nV+nP), nV+nP, 1)
      while (sum(n_eq < 0)){
        # Matrix A interactions
        #A_VV <- abs(matrix(rnorm(nV*nV, mean = alpha, sd = sigma*alpha), nV, nV))  # If we want our alphas to be normally distributed
        A_VV <- exp(matrix(rnorm(nV*nV, mean = log(alpha/sqrt(1+sigma**2)), sd = sqrt(log(1+sigma**2))), nV, nV)) # If we want our alphas to be log-normally distributed
        diag(A_VV) <- rep(1, nV)
        
        alphas_vector <- A_VV[row(A_VV) != col(A_VV)]
        
        mean_alpha   = mean(alphas_vector);   var_alpha    = var(alphas_vector)
        mean_alphaM1 = mean(1/alphas_vector); var_alphaM1  = var(1/alphas_vector)
        min_alpha = min(1/alphas_vector);     max_alpha = max(1/alphas_vector)
        mean_log10_alpha = mean(log10(alphas_vector)); var_log10_alpha = var(log10(alphas_vector))
        
        rm(alphas_vector)
        
        A_PP <- diag(nP)
        if (mode_predators == 'generalists'){
          #A_VP <- (1./nV)*abs(matrix(rnorm(nV*nP, mean = lambda, sd = sigma*lambda), nV, nP))    # If we want our alphas to be normally distributed
          A_VP <- (1./nV)*exp(matrix(rnorm(nV*nP, mean = log(lambda/sqrt(1+sigma**2)), sd = sqrt(log(1+sigma**2))), nV, nP)) # If we want our alphas to be log-normally distributed
          
          lambdas_vector <- c(nV*A_VP);
          
          mean_lambda   = mean(lambdas_vector); var_lambda= var(lambdas_vector)
          min_lambda    = min(lambdas_vector);  max_lambda    = max(lambdas_vector)
          mean_lambdaM1 = mean(1/lambdas_vector); var_lambdaM1  = var(1/lambdas_vector);
          mean_log10_lambda = mean(log10(lambdas_vector)); var_log10_lambda = var(log10(lambdas_vector))
          
          rm(lambdas_vector)
          
        } else if (mode_predators == 'specialists') {
          #A_VP <- diag(x=1, nrow = nV, ncol = nP) * abs(matrix(rnorm(nV*nP, mean = lambda, sd = sigma*lambda), nV, nP)) # If we want our alphas to be normally distributed
          A_VP <- diag(x=1, nrow = nV, ncol = nP) * exp(matrix(rnorm(nV*nP, mean = log(lambda/sqrt(1+sigma**2)), sd = sqrt(log(1+sigma**2))), nV, nP)) # If we want our alphas to be log-normally distributed
          
          lambdas_vector <- c(diag(A_VP));
          
          mean_lambda   = mean(lambdas_vector); var_lambda= var(lambdas_vector)
          min_lambda    = min(lambdas_vector);  max_lambda    = max(lambdas_vector)
          mean_lambdaM1 = mean(1/lambdas_vector); var_lambdaM1  = var(1/lambdas_vector);
          mean_log10_lambda = mean(log10(lambdas_vector)); var_log10_lambda = var(log10(lambdas_vector))
          
          rm(lambdas_vector)
        }
        #A_PV <- -1. * t(A_VP) * abs(matrix(rnorm(nV*nP, mean = eta, sd = sigma*eta), nP, nV))  # If we want our alphas to be normally distributed
        eta_PV <- exp(matrix(rnorm(nV*nP, mean = log(eta/sqrt(1+sigma**2)), sd = sqrt(log(1+sigma**2))), nP, nV)) # If we want our alphas to be log-normally distributed
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
        
        rm(etas_vector)
        
        rm(A_VV, A_VP, A_PV, A_PP)
        
        # Growth and mortality rates
        #bS <- abs(rnorm(nV, mean = b, sd = b*sigma))   # If we want our parameters to be normally distributed
        #mS <- abs(rnorm(nV, mean = m, sd = m*sigma))   # If we want our parameters to be normally distributed
        #dS <- abs(rnorm(nP, mean = d, sd = d*sigma))   # If we want our parameters to be normally distributed
        bS <- exp(rnorm(nV, mean = log(b/sqrt(1+sigma**2)), sd = sqrt(log(1+sigma**2)))) # If we want our parameters to be log-normally distributed
        mS <- exp(rnorm(nV, mean = log(m/sqrt(1+sigma**2)), sd = sqrt(log(1+sigma**2)))) # If we want our parameters to be log-normally distributed
        dS <- exp(rnorm(nP, mean = log(d/sqrt(1+sigma**2)), sd = sqrt(log(1+sigma**2)))) # If we want our parameters to be log-normally distributed
        
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
        #tauVS <- abs(rnorm(nV, mean = tauV, sd = tauV*sigma)) #rep(tauV, nV)   # If we want our taus to be normally distributed
        #tauPS <- abs(rnorm(nP, mean = tauP, sd = tauP*sigma)) #rep(tauP, nP)    # If we want our taus to be normally distributed
        tauVS <- exp(rnorm(nV, mean = log(tauV/sqrt(1+sigma**2)), sd = sqrt(log(1+sigma**2)))) #If we want our taus to be log-normally distributed
        tauPS <- exp(rnorm(nP, mean = log(tauP/sqrt(1+sigma**2)), sd = sqrt(log(1+sigma**2)))) #If we want our taus to be log-normally distributed
        
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
        n_eq <- inv(A_total) %*% R_total
        
        # Matrix with predator boolean array
        bool_predator = rbind(
          matrix(rep(0, nV),nV, 1), 
          matrix(rep(1, nP), nP, 1)
        )
      }
      
      # Save data
      row <- data.frame(idx = idx, stressor = 0, nV = nV, nP = nP,
                        nV_max = nV_max, nP_max = nP_max,
                        b = b, m = m, d = d, alpha = alpha,
                        lambda = lambda, eta = eta, 
                        tauV=tauV, tauP = tauP, sigma = sigma,
                        mode_predators = mode_predators,
                        
                        mean_b = mean_b, var_b = var_b, 
                        mean_bM1 = mean_bM1,  var_bM1 = var_bM1, 
                        min_b = min_b, max_b = max_b,
                        mean_log10_b = mean_log10_b, var_log10_b = var_log10_b, 
                        
                        mean_m = mean_m, var_m = var_m, 
                        mean_mM1 = mean_mM1,  var_mM1 = var_mM1, 
                        min_m = min_m, max_m = max_m,
                        mean_log10_m = mean_log10_m, var_log10_m = var_log10_m, 
                        
                        
                        mean_d = mean_d, var_d = var_d, 
                        mean_dM1 = mean_dM1,  var_dM1 = var_dM1, 
                        min_d = min_d, max_d = max_d,
                        mean_log10_d = mean_log10_d, var_log10_d = var_log10_d, 

                        mean_alpha = mean_alpha, var_alpha = var_alpha, 
                        mean_alphaM1 = mean_alphaM1,  var_alphaM1 = var_alphaM1,
                        min_alpha = min_alpha, max_alpha = max_alpha,
                        mean_log10_alpha = mean_log10_alpha, var_log10_alpha = var_log10_alpha, 
                        
                        mean_lambda = mean_lambda, var_lambda = var_lambda, 
                        mean_lambdaM1 = mean_lambdaM1,  var_lambdaM1 = var_lambdaM1,
                        min_lambda = min_lambda, max_lambda = max_lambda,
                        mean_log10_lambda = mean_log10_lambda, var_log10_lambda = var_log10_lambda, 
                        
                        mean_eta = mean_eta, var_eta = var_eta, 
                        mean_etaM1 = mean_etaM1,  var_etaM1 = var_etaM1,
                        min_eta = min_eta, max_eta = max_eta,
                        mean_log10_eta = mean_log10_eta, var_log10_eta = var_log10_eta, 
                        
                        mean_tauV = mean_tauV, var_tauV = var_tauV, 
                        mean_tauVM1 = mean_tauVM1,  var_tauVM1 = var_tauVM1,
                        min_tauV = min_tauV, max_tauV = max_tauV,
                        mean_log10_tauV = mean_log10_tauV, var_log10_tauV = var_log10_tauV, 
                        
                        mean_tauP = mean_tauP, var_tauP = var_tauP, 
                        mean_tauPM1 = mean_tauPM1,  var_tauPM1 = var_tauPM1,
                        min_tauP = min_tauP, max_tauP = max_tauP,
                        mean_log10_tauP = mean_log10_tauP, var_log10_tauP = var_log10_tauP
                        
                        )
      
      write.table(row, file = 'extinction_ln_histograms_v3_max_min.csv', sep = ',',
                  append = TRUE, row.names = FALSE, col.names=!file.exists("extinction_ln_histograms_v3_max_min.csv"))
      rm(row)
      
      
      
      # For the different concentrations, compute the diversity
      zis <- 10^seq(-2, 3, 0.02) # Different concentrations 0.02
      contador = 0
      for (zi in zis){
        contador = contador + 1
        if ((nV>0)|(nP>0)){   # If there are both prey or predator species
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
            n_eq <- inv(A_total) %*% R_total
          } else if (length(A_total)==1){
            n_eq <- 1/(A_total[1]) %*% R_total
          } else {
            n_eq = c()
          }
          
          # Count how many species have positive abundances
          while (sum(n_eq < 0) && nV > 0){
            if ((nV>0) | (nP>0)){
              
              if (length(n_eq)>1){
                A_total <- A_total[n_eq>0, n_eq>0]
              } else {
                A_total <- A_total[n_eq>0]
              }
              
              R_total <- R_total[n_eq>0]
              bool_predator <- bool_predator[n_eq>0]
              bS <- bS[(n_eq>0)[1:nV]]
              mS <- mS[(n_eq>0)[1:nV]]
              dS <- dS[(n_eq>0)[(nV+1):(nV+nP)]]
              
              zS <- rep(zi, length(R_total))
              
              
              
              tauS <- rbind(matrix(tauVS,nV,1), matrix(tauPS,nP,1))
              tauS <- tauS[n_eq>0]
              
              nV <- sum(bool_predator==0)
              nP <- sum(bool_predator==1)
              
              if (nV>0){
                if (nP>0){
                  tauVS <- tauS[1:nV]
                  tauPS <- tauS[(nV+1):(nV+nP)]
                  
                  
                  R_total <- rbind( 
                    matrix(bS-mS*(1+zS[1:nV]/tauVS), nV, 1),
                    matrix(-dS*(1+zS[(nV+1):(nV+nP)]/tauPS),nP,1)
                  )
                  
                  if (length(A_total)>1){
                    n_eq <- inv(A_total) %*% R_total
                  } else if (length(A_total)==1){
                    n_eq <- 1/(A_total[1]) %*% R_total
                  } else {
                    n_eq = c()
                  }
                } else {
                  tauVS <- tauS[1:nV]
                  R_total <- matrix(bS-mS*(1+zS[1:nV]/tauVS), nV, 1)
                  if (length(A_total)>1){
                    n_eq <- inv(A_total) %*% R_total
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
        row <- data.frame(idx = idx, stressor = zi, nV = nV, nP = nP,
                          nV_max = nV_max, nP_max = nP_max,
                          b = b, m = m, d = d, alpha = alpha,
                          lambda = lambda, eta = eta, 
                          tauV=tauV, tauP = tauP, sigma = sigma,
                          mode_predators = mode_predators,
                          
                          mean_b = mean_b, var_b = var_b, 
                          mean_bM1 = mean_bM1,  var_bM1 = var_bM1, 
                          min_b = min_b, max_b = max_b,
                          mean_log10_b = mean_log10_b, var_log10_b = var_log10_b, 
                          
                          mean_m = mean_m, var_m = var_m, 
                          mean_mM1 = mean_mM1,  var_mM1 = var_mM1, 
                          min_m = min_m, max_m = max_m,
                          mean_log10_m = mean_log10_m, var_log10_m = var_log10_m, 
                          
                          
                          mean_d = mean_d, var_d = var_d, 
                          mean_dM1 = mean_dM1,  var_dM1 = var_dM1, 
                          min_d = min_d, max_d = max_d,
                          mean_log10_d = mean_log10_d, var_log10_d = var_log10_d, 
                          
                          mean_alpha = mean_alpha, var_alpha = var_alpha, 
                          mean_alphaM1 = mean_alphaM1,  var_alphaM1 = var_alphaM1,
                          min_alpha = min_alpha, max_alpha = max_alpha,
                          mean_log10_alpha = mean_log10_alpha, var_log10_alpha = var_log10_alpha, 
                          
                          mean_lambda = mean_lambda, var_lambda = var_lambda, 
                          mean_lambdaM1 = mean_lambdaM1,  var_lambdaM1 = var_lambdaM1,
                          min_lambda = min_lambda, max_lambda = max_lambda,
                          mean_log10_lambda = mean_log10_lambda, var_log10_lambda = var_log10_lambda, 
                          
                          mean_eta = mean_eta, var_eta = var_eta, 
                          mean_etaM1 = mean_etaM1,  var_etaM1 = var_etaM1,
                          min_eta = min_eta, max_eta = max_eta,
                          mean_log10_eta = mean_log10_eta, var_log10_eta = var_log10_eta, 
                          
                          mean_tauV = mean_tauV, var_tauV = var_tauV, 
                          mean_tauVM1 = mean_tauVM1,  var_tauVM1 = var_tauVM1,
                          min_tauV = min_tauV, max_tauV = max_tauV,
                          mean_log10_tauV = mean_log10_tauV, var_log10_tauV = var_log10_tauV, 
                          
                          mean_tauP = mean_tauP, var_tauP = var_tauP, 
                          mean_tauPM1 = mean_tauPM1,  var_tauPM1 = var_tauPM1,
                          min_tauP = min_tauP, max_tauP = max_tauP,
                          mean_log10_tauP = mean_log10_tauP, var_log10_tauP = var_log10_tauP 
                          
        )
        write.table(row, file = 'extinction_ln_histograms_v3_max_min.csv', sep = ',',
                    append = TRUE, row.names = FALSE, col.names=!file.exists("extinction_ln_histograms_v3_max_min.csv"))
        rm(row)
      }
    }
  } 
}



