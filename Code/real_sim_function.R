real_sim_function = function(long.data, theta = 1){
  ### A dataframe of long.data contains the following variables:
  ### ID, city, PeriodEnd, year, Units_TaxInd, Units_TaxedFam, Units_NotTaxInd, 
  ### Units_NotTaxedFam, ind.SS, ind.AS, family.SS, family.AS, IncomePerHousehold, 
  ### MaleProportion, AsianProportion, HispanicProportion, BlackProportion, obs.Y, obs.A2, obs.A3, time, obs.Z2, obs.Z3
  
  results.list = list()
  # fit the outcome model
  fit.outcome = lmer(obs.Y ~ obs.A2 + obs.A3 + time +  # treatment group indicator, 
                       obs.Z2 + obs.Z3 + # treatment effects (main effect)
                       family.SS + family.AS + ind.SS + ind.AS +  # baseline covariates (main effect)
                       BlackProportion + HispanicProportion + AsianProportion + 
                       MaleProportion + IncomePerHousehold + 
                       obs.Z2*ind.SS + obs.Z3*ind.SS + # (treatment effects) x (baseline covariates)
                       obs.Z2*family.SS + obs.Z3*family.SS + 
                       obs.Z2*ind.AS + obs.Z3*ind.AS +
                       obs.Z2*family.AS + obs.Z3*family.AS + 
                       (1|ID), # unit-level random effects
                     data = long.data)
  
  # predicted outcomes conditional on Baltimore
  long.control = long.data
  long.control$obs.Z2 = 0; long.control$obs.Z3 = 0; long.control$obs.A2 = 0; long.control$obs.A3 = 0
  muhat_control = predict(fit.outcome, newdata = long.control, re.form = ~ (1|ID))
  muhat_control0 = muhat_control[long.data$year == 2016]; 
  muhat_control1 = muhat_control[long.data$year == 2017]
  
  # predicted outcomes conditional on neighboring controls
  long.neighbor = long.data
  long.neighbor$obs.Z2 = 0; long.neighbor$obs.Z3 = ifelse(long.neighbor$time == 1, 1, 0); long.neighbor$obs.A2 = 0; long.neighbor$obs.A3 = 1
  muhat_neighbor = predict(fit.outcome, newdata = long.neighbor)
  muhat_neighbor0 = muhat_neighbor[long.data$year == 2016]; 
  muhat_neighbor1 = muhat_neighbor[long.data$year == 2017]
  
  # predicted outcomes conditional on Philadelphia
  long.treat = long.data;
  long.treat$obs.Z2 = ifelse(long.treat$time == 1, 1, 0); long.treat$obs.Z3 = 0; long.treat$obs.A2 = 1; long.treat$obs.A3 = 0
  muhat_treat = predict(fit.outcome, newdata = long.treat)
  muhat_treat0 = muhat_treat[long.data$year == 2016]; 
  muhat_treat1 = muhat_treat[long.data$year == 2017]
  
  ## short data
  short.data = long.data[long.data$year == 2016,]
  short.data$obs.Y0 = long.data$obs.Y[long.data$year == 2016]
  short.data$obs.Y1 = long.data$obs.Y[long.data$year == 2017]
  short.data$obs.A = ifelse(short.data$obs.A2 == 1, 2, ifelse(short.data$obs.A3 == 1, 3, 1))
  
  ## propensity score model
  prop.fit = multinom(obs.A ~ family.SS + family.AS + ind.SS + ind.AS + 
                        BlackProportion + HispanicProportion + AsianProportion + 
                        MaleProportion + IncomePerHousehold, data = short.data, trace = FALSE)
  pihat_control = fitted.values(prop.fit)[,1]
  pihat_treat = fitted.values(prop.fit)[,2]
  pihat_neighbor = fitted.values(prop.fit)[,3]
  
  pihat_control = ifelse(pihat_control < 0.01, 0.01, pihat_control)
  pihat_treat = ifelse(pihat_treat < 0.01, 0.01, pihat_treat)
  pihat_neighbor = ifelse(pihat_neighbor < 0.01, 0.01, pihat_neighbor)
  etahat = pihat_neighbor / pihat_treat
  
  etahat = ifelse(etahat > 20, 20, ifelse(etahat < 0.05, 0.05, etahat))
  
  
  ## ATT ##
  control.A = (short.data$obs.A == 1); treat.A = (short.data$obs.A == 2); neighbor.A =(short.data$obs.A == 3)
  a1 = mean(short.data$obs.A == 2)
  a2 =  mean(((short.data$obs.A == 2) - pihat_treat*(short.data$obs.A == 1)/pihat_control)*( (short.data$obs.Y1 - short.data$obs.Y0) - (muhat_control1 - muhat_control0)) )
  dr.ATT = a2/a1  # doubly robust estimator
  # regression-based estimator
  reg.ATT = mean((short.data$obs.A == 2)*(short.data$obs.Y1 - short.data$obs.Y0))/mean(short.data$obs.A == 2) -  
    mean((short.data$obs.A == 2)*(muhat_control1 - muhat_control0))/mean(short.data$obs.A == 2)
  # ipw estimator
  ipw.ATT =  mean((short.data$obs.Y1 - short.data$obs.Y0)*(treat.A/pihat_treat - control.A/pihat_control)*pihat_treat)/a1
  # estimated efficient influence function
  phi.ATT =  ((short.data$obs.A == 2)*( (short.data$obs.Y1 - short.data$obs.Y0) - (muhat_treat1 - muhat_treat0)))/mean(short.data$obs.A == 2) - 
    (pihat_treat*(short.data$obs.A == 1)/pihat_control*( (short.data$obs.Y1 - short.data$obs.Y0) - (muhat_control1 - muhat_control0)))/mean(short.data$obs.A == 2) + 
    (short.data$obs.A==2)*(muhat_treat1 - muhat_treat0 - (muhat_control1 - muhat_control0))/mean(short.data$obs.A == 2) - (short.data$obs.A==2)*dr.ATT/mean(short.data$obs.A==2) 
  dr.ATT.var = mean(phi.ATT^2)/nrow(short.data)
  
  ATT.results = list(ipw.ATT, reg.ATT, dr.ATT, dr.ATT.var)
  
  
  ## ATN ##
  a1 = mean(short.data$obs.A == 3)
  a2 = mean(((short.data$obs.A == 3) - pihat_neighbor*(short.data$obs.A == 1)/pihat_control)*( (short.data$obs.Y1 - short.data$obs.Y0) - (muhat_control1 - muhat_control0)) )
  dr.ATN = a2/a1 # doubly-robust estimator
  # regression-based estimator
  reg.ATN = mean((short.data$obs.A == 3)*(short.data$obs.Y1 - short.data$obs.Y0)) /mean(short.data$obs.A == 3) -  
    mean((short.data$obs.A == 3)*(muhat_control1 - muhat_control0))/mean(short.data$obs.A == 3)
  # ipw estimator
  ipw.ATN = mean((short.data$obs.Y1 - short.data$obs.Y0)*(neighbor.A/pihat_neighbor - control.A/pihat_control)*pihat_neighbor)/a1
  # estimated efficient influence function
  phi.ATN = ((short.data$obs.A == 3)*( (short.data$obs.Y1 - short.data$obs.Y0) - (muhat_neighbor1 - muhat_neighbor0)))/mean(short.data$obs.A == 3) - 
    (pihat_neighbor*(short.data$obs.A == 1)/pihat_control*( (short.data$obs.Y1 - short.data$obs.Y0) - (muhat_control1 - muhat_control0)))/mean(short.data$obs.A == 3) + 
    (short.data$obs.A==3)*(muhat_neighbor1 - muhat_neighbor0 - (muhat_control1 - muhat_control0))/mean(short.data$obs.A == 3) -   (short.data$obs.A==3)*dr.ATN/mean(short.data$obs.A==3) 
  dr.ATN.var = mean(phi.ATN^2)/nrow(short.data)
  
  ATN.results = list(ipw.ATN, reg.ATN, dr.ATN, dr.ATN.var)
  
  ## -offsetting effect estimation
  outcome.part = mean(etahat*theta*pihat_treat*(neighbor.A/pihat_neighbor - control.A/pihat_control)*(short.data$obs.Y1 - short.data$obs.Y0))/mean(short.data$obs.A == 2)
  
  a1 = mean(etahat*theta*(short.data$obs.A == 2)*(muhat_neighbor1-muhat_neighbor0))/mean(short.data$obs.A == 2)
  b1 = mean(etahat*theta*pihat_treat*(muhat_neighbor1-muhat_neighbor0))/mean(short.data$obs.A == 2)
  c1 = mean(etahat*theta*pihat_treat*((short.data$obs.A ==3)/pihat_neighbor)*(muhat_neighbor1-muhat_neighbor0))/mean(short.data$obs.A == 2)
  d1 = mean(etahat*theta*pihat_treat*((short.data$obs.A ==1)/pihat_control)*(muhat_neighbor1-muhat_neighbor0))/mean(short.data$obs.A == 2)
  
  a2 = mean(etahat*theta*(short.data$obs.A == 2)*(muhat_control1-muhat_control0))/mean(short.data$obs.A == 2)
  b2 = mean(etahat*theta*pihat_treat*(muhat_control1-muhat_control0))/mean(short.data$obs.A == 2)
  c2 = mean(etahat*theta*pihat_treat*((short.data$obs.A ==3)/pihat_neighbor)*(muhat_control1-muhat_control0))/mean(short.data$obs.A == 2)
  d2 = mean(etahat*theta*pihat_treat*((short.data$obs.A ==1)/pihat_control)*(muhat_control1-muhat_control0))/mean(short.data$obs.A == 2)
  
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
  
  # regression-based estimator
  reg.offset =  a1 - a2
  # ipw estimator
  ipw.offset = outcome.part
  offset.results = list(ipw.offset, reg.offset, dr.offset)
  
  ## AOTT ##
  reg.AOTT = reg.ATT + reg.offset
  dr.AOTT = dr.ATT + dr.offset
  ipw.AOTT = ipw.ATT + ipw.offset
  
  AOTT.results = list(ipw.AOTT, reg.AOTT, dr.AOTT)
  results.list = list(ATT.results, ATN.results, offset.results, AOTT.results)
  
  return(results.list)
}