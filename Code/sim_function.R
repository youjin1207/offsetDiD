sim_function = function(obs.data, case = 1){
  ## obs.data: a dataframe with obs.A, X, obs.Y0, and obs.Y1
  #### obs.A : a treatment group indicator (1: control, 2: treatment, 3: neighbors)
  #### X : a vector of observed covariates, e.g., (X.1, X.2)
  #### obs.Y0 : the observed outcome at t=0 (at the pre-intervention time periods)
  #### obs.Y1 : the observed outcome at t=1 (at the post-intervention time periods)
  
  ## case = 1 when all the propensity scores and the outcome regression are correctly specified 
  ## case = 2 when the outcome regression model is misspecified
  ## case = 3 the propensity scores are misspecified but the eta is correctly specified 
  ## case = 4 all the propensity scores are misspecified 
  
  N = nrow(obs.data)
  long.data = data.frame(obs.Y = c(obs.data$obs.Y0, obs.data$obs.Y1), obs.A = rep(obs.data$obs.A,2), X1 = rep(obs.data$X.1, 2),
                         X2 = rep(obs.data$X.2, 2), time = c(rep(0, N), rep(1, N)), 
                         id = rep(1:N, 2), cm = rep(obs.data$cm, 2)) # make it as a long data 
  long.data$Z2 = ifelse(long.data$time > 0 & long.data$obs.A == 2, 1, 0) # an indicator of direct interventions
  long.data$Z3 = ifelse(long.data$time > 0 & long.data$obs.A == 3, 1, 0) # an indicator of indirect interventions
  long.data$obs.A2 = ifelse(long.data$obs.A == 2, 1, 0) # an indicator of being in a treatment group
  long.data$obs.A3 = ifelse(long.data$obs.A == 3, 1, 0) # an indicator of being in a neighboring control group
  
  cluster.data = aggregate(cbind(X.1, X.2, D) ~ cm, data = obs.data, mean)
  names(cluster.data) = c("cm", "X1.mean", "X2.mean", "D")
  subind.data = obs.data[obs.data$obs.A!=1, ]
  subind.data$treat = (subind.data$obs.A == 2)
  
  if(case %in% c(1,2,3,4)){
    # fit a multinomial propensity score model (correctly specified)
    fit1 = glm(D ~ X1.mean + X2.mean, data = cluster.data, family = binomial())
    fit2 = glm(treat ~ X.1 + X.2, data = subind.data, family = binomial())
    cluster.prop = rep(NA, nrow(obs.data))
    for(j in 1:nrow(cluster.data)) cluster.prop[obs.data$cm %in% cluster.data$cm[j]] = fitted.values(fit1)[j]
    #cluster.prop = rep(fitted.values(fit1), each = n)
    ind.prop = predict(fit2, newdata = obs.data, type = "response") 
    pihat_control = (1-cluster.prop)
    pihat_treat = cluster.prop*ind.prop
    pihat_neighbor = cluster.prop*(1-ind.prop)
    
    pihat_control = ifelse(pihat_control < 0.01, 0.01, pihat_control)
    pihat_treat = ifelse(pihat_treat < 0.01, 0.01, pihat_treat)
    pihat_neighbor = ifelse(pihat_neighbor < 0.01, 0.01, pihat_neighbor)
    
    etahat = pihat_neighbor / pihat_treat
  }
  
  if(case %in% c(3,4)){
    # fit a multinomial propensity score models (misspecified)
    fit1.wrong = glm(D ~ exp(-X2.mean^2), data = cluster.data, family = binomial())
    cluster.prop.wrong = rep(NA, nrow(obs.data))
    for(j in 1:nrow(cluster.data)) cluster.prop.wrong[obs.data$cm %in% cluster.data$cm[j]] = fitted.values(fit1.wrong)[j]
    pihat_control = (1-cluster.prop.wrong)
    pihat_treat = cluster.prop.wrong*ind.prop
    pihat_neighbor = cluster.prop.wrong*(1-ind.prop)
    
    pihat_control = ifelse(pihat_control < 0.01, 0.01, pihat_control)
    pihat_treat = ifelse(pihat_treat < 0.01, 0.01, pihat_treat)
    pihat_neighbor = ifelse(pihat_neighbor < 0.01, 0.01, pihat_neighbor)
    
    etahat = pihat_neighbor / pihat_treat
    
    if(case == 4){
      fit2.wrong = glm(treat ~ exp(-X.2^2), data = subind.data, family = binomial())
      ind.prop.wrong = predict(fit2.wrong, newdata = obs.data, type = "response") 
      pihat_treat = cluster.prop.wrong*ind.prop.wrong
      pihat_neighbor = cluster.prop.wrong*(1-ind.prop.wrong)
      
      pihat_treat = ifelse(pihat_treat < 0.01, 0.01, pihat_treat)
      pihat_neighbor = ifelse(pihat_neighbor < 0.01, 0.01, pihat_neighbor)
      
      etahat = pihat_neighbor / pihat_treat
    }
  }
  
  if(case %in% c(1,3,4)){
    # fit the outcome model (correctly specified)
    fit.outcome = lmer(obs.Y ~ obs.A2 + obs.A3 + time + 
                         X1 + X2 + 
                         time*X1 + time*X2 + 
                         Z2 + Z3 +  
                         Z2*X1 + Z2*X2  + 
                         Z3*X1 + Z3*X2 + as.factor(cm) + (1|id) , data = long.data)
  }
  
  if(case == 2){
    # fit the outcome model (misspecified)
    fit.outcome = lmer(obs.Y ~ obs.A2 + obs.A3 + time + 
                         X1 + #X2 + 
                         time*X1 + #time*X2 + 
                         Z2 + Z3 +  
                         Z2*X1 + #Z2*X2  + 
                         Z3*X1 + #Z3*X2 + 
                         as.factor(cm) + (1|id) , data = long.data)
  }
  
  
  ## predicted outcomes (conditional on the non-neighboring control group)
  long.control = long.data
  long.control$Z2 = 0; long.control$Z3 = 0; long.control$obs.A2 = 0; long.control$obs.A3 = 0
  muhat_control = predict(fit.outcome, newdata = long.control)
  muhat_control0 = muhat_control[1:N]; muhat_control1 = muhat_control[(N+1):(2*N)]
  
  ## predicted outcomes (conditional on the neighboring control group)
  long.neighbor = long.data
  long.neighbor$Z2 = 0; long.neighbor$Z3 = ifelse(long.neighbor$time == 1, 1, 0); long.neighbor$obs.A2 = 0; long.neighbor$obs.A3 = 1
  muhat_neighbor = predict(fit.outcome, newdata = long.neighbor)
  muhat_neighbor0 = muhat_neighbor[1:N]; muhat_neighbor1 = muhat_neighbor[(N+1):(2*N)]
  
  ## predicted outcomes (conditional on the treatment group)
  long.treat = long.data
  long.treat$Z2 = ifelse(long.treat$time == 1, 1, 0); long.treat$Z3 = 0; long.treat$obs.A2 = 1; long.treat$obs.A3 = 0
  muhat_treat = predict(fit.outcome, newdata = long.neighbor)
  muhat_treat0 = muhat_treat[1:N]; muhat_treat1 = muhat_treat[(N+1):(2*N)]
  
  ## -offsetting effect estimation (21 conditionally doubly robust estimators)
  outcome.part = mean(etahat*pihat_treat*((obs.data$obs.A==3)/pihat_neighbor - (obs.data$obs.A == 1)/pihat_control)*(obs.data$obs.Y1 - obs.data$obs.Y0))/mean(obs.data$obs.A == 2)
  
  a1 = mean(etahat*(obs.data$obs.A == 2)*(muhat_neighbor1-muhat_neighbor0))/mean(obs.data$obs.A == 2)
  b1 = mean(etahat*pihat_treat*(muhat_neighbor1-muhat_neighbor0))/mean(obs.data$obs.A == 2)
  c1 = mean(etahat*pihat_treat*((obs.data$obs.A ==3)/pihat_neighbor)*(muhat_neighbor1-muhat_neighbor0))/mean(obs.data$obs.A == 2)
  d1 = mean(etahat*pihat_treat*((obs.data$obs.A ==1)/pihat_control)*(muhat_neighbor1-muhat_neighbor0))/mean(obs.data$obs.A == 2)
  
  a2 = mean(etahat*(obs.data$obs.A == 2)*(muhat_control1-muhat_control0))/mean(obs.data$obs.A == 2)
  b2 = mean(etahat*pihat_treat*(muhat_control1-muhat_control0))/mean(obs.data$obs.A == 2)
  c2 = mean(etahat*pihat_treat*((obs.data$obs.A ==3)/pihat_neighbor)*(muhat_control1-muhat_control0))/mean(obs.data$obs.A == 2)
  d2 = mean(etahat*pihat_treat*((obs.data$obs.A ==1)/pihat_control)*(muhat_control1-muhat_control0))/mean(obs.data$obs.A == 2)
  
  dr.offset = c()
  for(k in 1:27){
    dr.offset[1] = outcome.part - (a1-d2) + (a1-a2)
    dr.offset[2] = outcome.part - (a1-d2) + (a1-b2)
    dr.offset[3] = outcome.part - (a1-d2) + (a1-c2)
    dr.offset[4] = outcome.part - (a1-d2) + (b1-a2)
    dr.offset[5] = outcome.part - (a1-d2) + (b1-b2)
    dr.offset[6] = outcome.part - (a1-d2) + (b1-c2)
    dr.offset[7] = outcome.part - (a1-d2) + (c1-a2)
    dr.offset[8] = outcome.part - (a1-d2) + (c1-b2)
    dr.offset[9] = outcome.part - (a1-d2) + (c1-c2)
    dr.offset[10] = outcome.part - (b1-d2) + (a1-a2)
    dr.offset[11] = outcome.part - (b1-d2) + (a1-b2)
    dr.offset[12] = outcome.part - (b1-d2) + (a1-c2)
    dr.offset[13] = outcome.part - (b1-d2) + (b1-a2)
    dr.offset[14] = outcome.part - (b1-d2) + (b1-b2)
    dr.offset[15] = outcome.part - (b1-d2) + (b1-c2)
    dr.offset[16] = outcome.part - (b1-d2) + (c1-a2)
    dr.offset[17] = outcome.part - (b1-d2) + (c1-b2)
    dr.offset[18] = outcome.part - (b1-d2) + (c1-c2)
    dr.offset[19] = outcome.part - (c1-d2) + (a1-a2)
    dr.offset[20] = outcome.part - (c1-d2) + (a1-b2)
    dr.offset[21] = outcome.part - (c1-d2) + (a1-c2)
    dr.offset[22] = outcome.part - (c1-d2) + (b1-a2)
    dr.offset[23] = outcome.part - (c1-d2) + (b1-b2)
    dr.offset[24] = outcome.part - (c1-d2) + (b1-c2)
    dr.offset[25] = outcome.part - (c1-d2) + (c1-a2)
    dr.offset[26] = outcome.part - (c1-d2) + (c1-b2)
    dr.offset[27] = outcome.part - (c1-d2) + (c1-c2)
  }
  
  ## regression-based estimator
  reg.offset =  a1 - a2
  ## ipw estimator
  ipw.offset = outcome.part
  
  ### ATT
  control.A = (obs.data$obs.A == 1); treat.A = (obs.data$obs.A == 2); neighbor.A =(obs.data$obs.A == 3)
  a1 = mean(obs.data$obs.A == 2)
  a2 =  mean(((obs.data$obs.A == 2) - pihat_treat*(obs.data$obs.A == 1)/pihat_control)*( (obs.data$obs.Y1 - obs.data$obs.Y0) - (muhat_control1 - muhat_control0)) )
  dr.ATT = a2/a1  # doubly robust estimator
  # regression-based estimator
  reg.ATT = mean((obs.data$obs.A == 2)*(obs.data$obs.Y1 - obs.data$obs.Y0))/mean(obs.data$obs.A == 2) -  
    mean((obs.data$obs.A == 2)*(muhat_control1 - muhat_control0))/mean(obs.data$obs.A == 2)
  # ipw estimator
  ipw.ATT =  mean((obs.data$obs.Y1 - obs.data$obs.Y0)*(treat.A/pihat_treat - control.A/pihat_control)*pihat_treat)/a1
  
  ## AOTT ##
  reg.AOTT = reg.ATT + reg.offset
  dr.AOTT = dr.ATT + dr.offset
  ipw.AOTT = ipw.ATT + ipw.offset
  
  return(list(ipw.offset = ipw.offset, reg.offset = reg.offset, dr.offset = dr.offset,
              ipw.AOTT = ipw.AOTT, reg.AOTT = reg.AOTT, dr.AOTT = dr.AOTT))
}
