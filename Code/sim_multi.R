library(MASS)
library(nnet)
library(lme4)
library(xtable)
###

## parameters ##
betas_control = c(0,0); betas_treat = c(1, 0.5); betas_neighbor = c(-0.5, -0.5)
gamma = seq(-1.0, 1.5, 0.1) 
# time-varying effects of the baseline covariate 1
lambda1 = c(c(-0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.2, 0.1, 0.0, -0.1, -0.2, -0.3), 
            c(-0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.2, 0.1, 0.0, -0.1, -0.2, -0.3) + 0.1); 
# time-varying effects of the baseline covariate 2
lambda2 = c(-c(-0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.2, 0.1, 0.0, -0.1, -0.2, -0.3),
            -c(-0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.2, 0.1, 0.0, -0.1, -0.2, -0.3)-0.3)
alpha_a = c(0,  0.5,  -1) # treatment group effect: control, treatment, neighbor
theta_neighbor = 1.0; theta_treat = -1.5 
# theta_neighbor = 0.5; theta_treat = -2.0 # scenario 2
time.effect.neighbor = seq(1, 0, -1/12)[1:13] # time-varying spillover effects
time.effect.treat = c(0.8, 0.9, seq(1.0, 0, -0.1)[1:11]) # time-varying treatment effects

## true values ##
NN = 10000*5
X = mvrnorm(NN, c(0,0), matrix(c(1, 0.3, 0.3, 1), 2, 2))
tmp = c(0, -0.5, 0.3) + X %*%  cbind(betas_control, betas_treat, betas_neighbor)
X[,1] = pmin(pmax(-2.0, round(10*X[,1])/10), 2.0)
X[,2] = pmin(pmax(-2.0, round(10*X[,2])/10), 2.0)
obs.A = matrix(0, NN, 3)
for(i in 1:NN){
  obs.A[i,] =  rmultinom(1, 1, c(exp(tmp)[i,1:3]))
}
obs.A = cbind(obs.A, rep(0, NN))
obs.A = apply(obs.A, 1, function(x) which(x == 1)) # 1: control, 2: treatment, 3: neighbors

true.ATT =  mean(mean(theta_treat*(1 - 0.5*X[obs.A==2,1]))*time.effect.treat)
true.ATN =  mean(mean(theta_neighbor*(1 + 0.5*X[obs.A==3, 1]) + theta_neighbor*(1-(X[obs.A==3,2])^2))*time.effect.neighbor)
true.offset = mean(mean(theta_neighbor*(1 + 0.5*X[obs.A==2, 1]) + theta_neighbor*(1-(X[obs.A==2,2])^2))*time.effect.neighbor)
true.AOTT =  true.ATT + true.offset


naive.neighbor = naive.treat = dr.ATT = dr.ATN = dr.offset = dr.AOTT = 
  dr.ATT.var = dr.ATN.var = dr.offset.var = dr.AOTT.var = dr.AOTT.var2 =  
  cr.ATT = cr.ATN = cr.offset = cr.AOTT = cr.AOTT2 = 
  reg.ATT = reg.ATN = reg.offset = reg.AOTT =  
  ipw.ATT = ipw.ATN = ipw.offset = ipw.AOTT = c()

############### Simulation Study ##################
for(r in 1:1000){
  
  set.seed(r)

  N = 2000 # N = 500, 1000
  epsilons = mvrnorm(N, c(0,0), matrix(c(1, 0.2, 0.2, 1), 2, 2))
  X = mvrnorm(N, c(0,0), matrix(c(1, 0.3, 0.3, 1), 2, 2))
  X[,1] = pmin(pmax(-2.0, round(10*X[,1])/10), 2.0)
  X[,2] = pmin(pmax(-2.0, round(10*X[,2])/10), 2.0)
  
  tmp = c(0, 0.0, 0.3) + X %*%  cbind(betas_control, betas_treat, betas_neighbor)
  obs.A = matrix(0, N, 3)
  for(i in 1:N){
    obs.A[i,] =  rmultinom(1, 1, c(exp(tmp)[i,1:3]))
  }
  obs.A = cbind(obs.A, rep(0, N))
  obs.A = apply(obs.A, 1, function(x) which(x == 1)) # 1: control, 2: treatment, 3: neighbors
  
  # true propensities
  #pi_treat = exp(tmp[,2]) / apply(exp(tmp[,1:3]), 1, sum)
  #pi_neighbor = exp(tmp[,3]) / apply(exp(tmp[,1:3]), 1, sum)
  #pi_control = exp(tmp[,1]) / apply(exp(tmp[,1:3]), 1, sum)
  
  obs.Y0 = obs.Y1 = matrix(NA, 13, N);   t.effect = matrix(NA, 13, N) 
  for(t in 1:13){
    # pre-treatment outcome
    obs.Y0[t,] = gamma[t] + X %*% c(lambda1[t], lambda2[t]) + alpha_a[obs.A] + epsilons[,1]
    # time-varying treatment effects
    t.effect[t,] = ifelse(obs.A == 1, 0, ifelse(obs.A == 2, theta_treat*(1 - 0.5*X[,1])*time.effect.treat[t],
                                                ifelse(obs.A == 3, (theta_neighbor*(1 + 0.5*X[, 1]) + theta_neighbor*(1-(X[,2]^2)))*time.effect.neighbor[t], NA)))
    # post-treatment outcome
    obs.Y1[t,] = gamma[t] + X %*% c(lambda1[t], lambda2[t]) + alpha_a[obs.A] + t.effect[t,] + epsilons[,2]
  }

  long.data = c()
  for(i in 1:N){
    long.data = rbind(long.data, cbind(c(obs.Y0[,i], obs.Y1[,i]), rep(obs.A[i], 26), rep(X[i,1], each = 26),
                                       rep(X[i,2], 26),  c(rep(0, 13), rep(1, 13)), c(1:26),
                                       rep(i, 26)))
  }
  colnames(long.data) = c("obs.Y", "obs.A", "X1", "X2", "time", "month", "id")
  long.data = as.data.frame(long.data)
  long.data$month = as.factor(long.data$month)
  
  
  long.data$Z2 = ifelse(long.data$time > 0 & long.data$obs.A == 2, 1, 0)
  long.data$Z3 = ifelse(long.data$time > 0 & long.data$obs.A == 3, 1, 0)
  long.data$obs.A2 = ifelse(long.data$obs.A == 2, 1, 0)
  long.data$obs.A3 = ifelse(long.data$obs.A == 3, 1, 0)
  
  # fit the outcome model (correctly specified)
  fit.outcome = lmer(obs.Y ~ obs.A2 + obs.A3 + time + 
                       X1 + X2 + 
                       month*X1 + month*X2 + 
                       Z2 + Z3 +  
                       Z2*X1 + Z2*X2  + 
                       Z3*X1 + Z3*X2 + Z3*I(X2^2) + 
                       Z2*X1*month + Z2*X2*month  + 
                       Z3*X1*month + Z3*X2*month + Z3*I(X2^2)*month + 
                       (1|id), data = long.data)
  
  # fit the outcome model (misspecified)
  #fit.outcome = lmer(obs.Y ~ obs.A2 + obs.A3 + time + 
  #                     X1 + X2 + 
  #                     #month*X1 + month*X2 + 
  #                     Z2 + Z3 +  
  #                     #Z2*X1 + Z2*X2  + 
  #                     #Z3*X1 + Z3*X2 + #Z3*I(X2^2) + 
  #                     #Z2*X1*month + Z2*X2*month  + 
  #                     #Z3*X1*month + Z3*X2*month + #Z3*I(X2^2)*month + 
  #                     (1|id), data = long.data)
  
  ## predicted outcomes (conditional on the non-neighboring control group)
  long.control = long.data
  long.control$Z2 = 0; long.control$Z3 = 0; long.control$obs.A2 = 0; long.control$obs.A3 = 0
  muhat_control = predict(fit.outcome, newdata = long.control)
  muhat_control0 = muhat_control[long.data$time == 0]; muhat_control1 = muhat_control[long.data$time == 1]
  
  ## predicted outcomes (conditional on the neighboring control group)
  long.neighbor = long.data
  long.neighbor$Z2 = 0; long.neighbor$Z3 = ifelse(long.neighbor$time == 1, 1, 0); long.neighbor$obs.A2 = 0; long.neighbor$obs.A3 = 1
  muhat_neighbor = predict(fit.outcome, newdata = long.neighbor)
  muhat_neighbor0 = muhat_neighbor[long.data$time == 0]; muhat_neighbor1 = muhat_neighbor[long.data$time == 1]
  
  ## predicted outcomes (conditional on the treatment group)
  long.treat = long.data
  long.treat$Z2 = ifelse(long.treat$time == 1, 1, 0); long.treat$Z3 = 0; long.treat$obs.A2 = 1; long.treat$obs.A3 = 0
  muhat_treat = predict(fit.outcome, newdata = long.treat)
  muhat_treat0 = muhat_treat[long.data$time == 0]; muhat_treat1 = muhat_treat[long.data$time == 1]
  
  ### short data ###
  short.data = long.data[long.data$time == 0,]
  short.data$obs.Y0 = long.data$obs.Y[long.data$time == 0]
  short.data$obs.Y1 = long.data$obs.Y[long.data$time == 1]
  short.data$obs.A = ifelse(short.data$obs.A2 == 1, 2, ifelse(short.data$obs.A3 == 1, 3, 1))
  
  ## propensity score model (correctly specifed)
  prop.fit = multinom(obs.A ~ X1 + X2, data = short.data)
  pihat_control = fitted.values(prop.fit)[,1]
  pihat_treat = fitted.values(prop.fit)[,2]
  pihat_neighbor = fitted.values(prop.fit)[,3]
  
  ## propensity score model (misspecified)
  #prop.fit = multinom(obs.A ~ exp(X2), data = short.data)
  #pihat_control = fitted.values(prop.fit)[,1]
  #pihat_treat = fitted.values(prop.fit)[,2]
  #pihat_neighbor = fitted.values(prop.fit)[,3]
  
  
  ## ATT estimation
  control.A = (short.data$obs.A == 1); treat.A = (short.data$obs.A == 2); neighbor.A =(short.data$obs.A == 3)
  a1 = mean(short.data$obs.A == 2)
  a2 =  mean(((short.data$obs.A == 2) - pihat_treat*(short.data$obs.A == 1)/pihat_control)*((short.data$obs.Y1 - short.data$obs.Y0) - (muhat_control1 - muhat_control0)) )
  dr.ATT[r] = a2/a1 # doubly robust estimator
  # regression-based estimator
  reg.ATT[r] = mean((short.data$obs.A == 2)*(short.data$obs.Y1 - short.data$obs.Y0))/mean(short.data$obs.A == 2) -  
    mean((short.data$obs.A == 2)*(muhat_control1 - muhat_control0))/mean(short.data$obs.A == 2)
  # ipw estimator
  ipw.ATT[r] =  mean((short.data$obs.Y1 - short.data$obs.Y0)*(treat.A/pihat_treat - control.A/pihat_control)*pihat_treat)/a1
  # efficient influence function
  phi.ATT =  ((short.data$obs.A == 2)*( (short.data$obs.Y1 - short.data$obs.Y0) - (muhat_treat1 - muhat_treat0)))/mean(short.data$obs.A == 2) - 
    (pihat_treat*(short.data$obs.A == 1)/pihat_control*( (short.data$obs.Y1 - short.data$obs.Y0) - (muhat_control1 - muhat_control0)))/mean(short.data$obs.A == 2) + 
    (short.data$obs.A==2)*(muhat_treat1 - muhat_treat0 - (muhat_control1 - muhat_control0))/mean(short.data$obs.A == 2) -   (short.data$obs.A==2)*dr.ATT[r]/mean(short.data$obs.A==2)                          
  dr.ATT.var[r] = var(phi.ATT)/N # variance estimator
  # 95% confidence interval coverage indicator
  cr.ATT[r] = (dr.ATT[r] - 1.96*sqrt(dr.ATT.var[r]) <= true.ATT & dr.ATT[r] + 1.96*sqrt(dr.ATT.var[r]) >= true.ATT ) 
  
  ## ATN estimation
  a1 = mean(short.data$obs.A == 3)
  a2 = mean(((short.data$obs.A == 3) - pihat_neighbor*(short.data$obs.A == 1)/pihat_control)*( (short.data$obs.Y1 - short.data$obs.Y0) - (muhat_control1 - muhat_control0)) )
  dr.ATN[r] = a2/a1 # doubly robust estimator
  # regression-based estimator
  reg.ATN[r] = mean((short.data$obs.A == 3)*(short.data$obs.Y1 - short.data$obs.Y0)) /mean(short.data$obs.A == 3) -  
    mean((short.data$obs.A == 3)*(muhat_control1 - muhat_control0))/mean(short.data$obs.A == 3)
  # ipw estimator
  ipw.ATN[r] = mean((short.data$obs.Y1 - short.data$obs.Y0)*(neighbor.A/pihat_neighbor - control.A/pihat_control)*pihat_neighbor)/a1
  # efficient influence function
  phi.ATN = ((short.data$obs.A == 3)*( (short.data$obs.Y1 - short.data$obs.Y0) - (muhat_neighbor1 - muhat_neighbor0)))/mean(short.data$obs.A == 3) - 
    (pihat_neighbor*(short.data$obs.A == 1)/pihat_control*( (short.data$obs.Y1 - short.data$obs.Y0) - (muhat_control1 - muhat_control0)))/mean(short.data$obs.A == 3) + 
    (short.data$obs.A==3)*(muhat_neighbor1 - muhat_neighbor0 - (muhat_control1 - muhat_control0))/mean(short.data$obs.A == 3) -   (short.data$obs.A==3)*dr.ATN[r]/mean(short.data$obs.A==3) 
  dr.ATN.var[r] = mean(phi.ATN^2)/N # variance estimator
  # 95% confidence interval coverage indicator
  cr.ATN[r] = (dr.ATN[r] - 1.96*sqrt(dr.ATN.var[r]) <= true.ATN & dr.ATN[r] + 1.96*sqrt(dr.ATN.var[r]) >= true.ATN ) 
  
  
  ## offsetting effect
  a2 = mean(((short.data$obs.A == 3)*pihat_treat/pihat_neighbor)*((short.data$obs.Y1 - short.data$obs.Y0 -  (muhat_neighbor1 - muhat_neighbor0) ))) + mean((short.data$obs.A==2)*(muhat_neighbor1 - muhat_neighbor0)) -  
    mean((control.A*pihat_treat/pihat_control)*((short.data$obs.Y1 - short.data$obs.Y0 -  (muhat_control1 - muhat_control0) ))) -  mean((short.data$obs.A==2)*(muhat_control1 - muhat_control0))
  dr.offset[r] = a2/ mean(short.data$obs.A == 2) # doubly-robust estimator
  # regression-based estimator
  reg.offset[r] =  mean((short.data$obs.A == 2)*(muhat_neighbor1 - muhat_neighbor0))/mean(short.data$obs.A == 2) -  
    mean((short.data$obs.A == 2)*(muhat_control1 - muhat_control0))/mean(short.data$obs.A == 2)
  # ipw estimator
  ipw.offset[r] = mean((short.data$obs.Y1 - short.data$obs.Y0)*((short.data$obs.A == 3)/pihat_neighbor - control.A/pihat_control)*pihat_treat)/mean(short.data$obs.A == 2)
  # efficient influence function
  phi.offset = ( ((short.data$obs.A == 3)*pihat_treat/pihat_neighbor)*((short.data$obs.Y1 - short.data$obs.Y0 -  (muhat_neighbor1 - muhat_neighbor0) )) + (short.data$obs.A==2)*(muhat_neighbor1 - muhat_neighbor0) - 
                   (control.A*pihat_treat/pihat_control)*((short.data$obs.Y1 - short.data$obs.Y0 -  (muhat_control1 - muhat_control0) )) - (short.data$obs.A==2)*(muhat_control1 - muhat_control0))/mean(short.data$obs.A == 2) - 
    dr.offset[r]
  # variance estimator
  dr.offset.var[r] = var(phi.offset)/N
  cr.offset[r] = (dr.offset[r] - 1.96*sqrt(dr.offset.var[r]) <= true.offset & dr.offset[r] + 1.96*sqrt(dr.offset.var[r]) >= true.offset ) 
  
  reg.AOTT[r] = reg.ATT[r] + reg.offset[r]
  dr.AOTT[r] = dr.ATT[r] + dr.offset[r]
  ipw.AOTT[r] = ipw.ATT[r] + ipw.offset[r]
  dr.AOTT.var[r] = var(phi.ATT + phi.offset)
  cr.AOTT[r] = (dr.AOTT[r] - 1.96*sqrt(dr.AOTT.var[r]/N) <= true.AOTT & dr.AOTT[r] + 1.96*sqrt(dr.AOTT.var[r]/N) >= true.AOTT ) 
 
}

summary(naive.treat); summary(naive.neighbor)

summary(dr.ATT); summary(dr.ATN); summary(dr.offset); summary(dr.AOTT)
summary(reg.ATT); summary(reg.ATN); summary(reg.offset); summary(reg.AOTT)
summary(ipw.ATT); summary(ipw.ATN); summary(ipw.offset); summary(ipw.AOTT)

bias.table = matrix(NA, nrow = 3, ncol= 4)
row.names(bias.table) = c("IPW", "Reg", "DR"); colnames(bias.table) = c("ATT", "ATN", "Offset", "AOTT")
bias.table[,1] = c(mean(ipw.ATT) - true.ATT, mean(reg.ATT) - true.ATT, mean(dr.ATT) - true.ATT)
bias.table[,2] = c(mean(ipw.ATN) - true.ATN, mean(reg.ATN) - true.ATN, mean(dr.ATN) - true.ATN)
bias.table[,3] = c(mean(ipw.offset) - true.offset, mean(reg.offset) - true.offset, mean(dr.offset) - true.offset)
bias.table[,4] = c(mean(ipw.AOTT) - true.AOTT, mean(reg.AOTT) - true.AOTT, mean(dr.AOTT) - true.AOTT)

se.table = matrix(NA, nrow = 3, ncol= 4)
row.names(se.table) = c("IPW", "Reg", "DR"); colnames(se.table) = c("ATT", "ATN", "Offset", "AOTT")
se.table[,1] = c(sd(ipw.ATT), sd(reg.ATT), sd(dr.ATT))
se.table[,2] = c(sd(ipw.ATN), sd(reg.ATN), sd(dr.ATN))
se.table[,3] = c(sd(ipw.offset), sd(reg.offset), sd(dr.offset))
se.table[,4] = c(sd(ipw.AOTT), sd(reg.AOTT), sd(dr.AOTT))


cr.table = matrix(NA, nrow = 3, ncol= 4)
row.names(cr.table) = c("IPW", "Reg", "DR"); colnames(cr.table) = c("ATT", "ATN", "Offset", "AOTT")
cr.table[,1] = c(mean(ipw.ATT >= true.ATT - 1.96*sd(ipw.ATT) &  ipw.ATT <= true.ATT + 1.96*sd(ipw.ATT)), 
                 mean(reg.ATT >= true.ATT - 1.96*sd(reg.ATT) &  reg.ATT <= true.ATT + 1.96*sd(reg.ATT)), 
                 mean(dr.ATT >= true.ATT - 1.96*sd(dr.ATT) & dr.ATT <= true.ATT + 1.96*sd(dr.ATT)))
cr.table[,2] = c(mean(ipw.ATN >= true.ATN - 1.96*sd(ipw.ATN) &  ipw.ATN <= true.ATN + 1.96*sd(ipw.ATN)), 
                 mean(reg.ATN >= true.ATN - 1.96*sd(reg.ATN) &  reg.ATN <= true.ATN + 1.96*sd(reg.ATN)), 
                 mean(dr.ATN >= true.ATN - 1.96*sd(dr.ATN) & dr.ATN <= true.ATN + 1.96*sd(dr.ATN)))
cr.table[,3] =  c(mean(ipw.offset >= true.offset - 1.96*sd(ipw.offset) &  ipw.offset <= true.offset + 1.96*sd(ipw.offset)), 
                  mean(reg.offset >= true.offset - 1.96*sd(reg.offset) &  reg.offset <= true.offset + 1.96*sd(reg.offset)), 
                  mean(dr.offset >= true.offset - 1.96*sd(dr.offset) & dr.offset <= true.offset + 1.96*sd(dr.offset)))
cr.table[,4] = c(mean(ipw.AOTT >= true.AOTT - 1.96*sd(ipw.AOTT) &  ipw.AOTT <= true.AOTT + 1.96*sd(ipw.AOTT)), 
                 mean(reg.AOTT >= true.AOTT - 1.96*sd(reg.AOTT) &  reg.AOTT <= true.AOTT + 1.96*sd(reg.AOTT)), 
                 mean(dr.AOTT >= true.AOTT - 1.96*sd(dr.AOTT) & dr.AOTT <= true.AOTT + 1.96*sd(dr.AOTT)))


print(xtable(bias.table, digits = 3))
print(cr.table)
print(se.table)

print(round(c(mean(dr.AOTT - true.AOTT),  mean(sqrt(dr.AOTT.var/N)), mean(cr.AOTT) ),3))
