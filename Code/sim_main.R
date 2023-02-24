library(MASS)
library(nnet)
library(lme4)
library(xtable)
###

## parameters ##
betas_control = c(0,0); betas_treat = c(1, 0.5); betas_neighbor = c(-0.5, -0.5)
gamma_0 = 1; gamma_1 = 1.5
lambda_0 = c(0.5, 0.5); lambda_1 = c(-0.5, 1)
alpha_a = c(0,  0.5,  -1) # group efffect: control, treatment, neighbor
theta_neighbor = 0.5; theta_treat = -1

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


true.ATT =  mean(theta_treat*(1 - 0.5*X[obs.A==2,1]))
true.ATN =  mean(theta_neighbor*(1 + 0.5*X[obs.A==3, 1]) + theta_neighbor*(1-2*(X[obs.A==3,2])^2))
true.offset = mean(theta_neighbor*(1 + 0.5*X[obs.A==2, 1]) + theta_neighbor*(1-2*(X[obs.A==2,2])^2)) # -(offset effects)
true.AOTT =  mean(theta_treat*(1-0.5*X[obs.A==2,1])) + mean(theta_neighbor*(1 + 0.5*X[obs.A==2, 1]) + theta_neighbor*(1-2*(X[obs.A==2,2])^2))




naive.neighbor = naive.treat = dr.ATT = dr.ATN = dr.offset = dr.AOTT = 
  dr.ATT.var = dr.ATN.var = dr.offset.var = dr.AOTT.var = dr.AOTT.var2 =  
  cr.ATT = cr.ATN = cr.offset = cr.AOTT = cr.AOTT2 = 
  reg.ATT = reg.ATN = reg.offset = reg.AOTT =  
  ipw.ATT = ipw.ATN = ipw.offset = ipw.AOTT = c()

############### Simulation Study ##################
for(r in 1:1000){
  
  set.seed(r)
  
  N = 2000 # N=500, 1000 (for Table 3 in the main manuscript; fix N = 2000 for Table S1 in the Supplementary Materials)
  epsilons = mvrnorm(N, c(0,0), matrix(c(1, 0.2, 0.2, 1), 2, 2)) # correlated within unit outcomes (two time points)
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
  
  ## generate outcome at t=0 (pre-treatment)
  obs.Y0 = gamma_0 + X %*% lambda_0 + alpha_a[obs.A] + epsilons[,1]
  ## generate the treatment effects
  t.effect = ifelse(obs.A == 1, 0, ifelse(obs.A == 2, theta_treat*(1 - 0.5*X[,1]),
                                          ifelse(obs.A == 3, theta_neighbor*(1 + 0.5*X[, 1]) + theta_neighbor*(1-2*(X[,2]^2)), NA)))
  # generate the outcome at t=1 (post-treatment)
  obs.Y1 = gamma_1 + X %*% lambda_1 + alpha_a[obs.A] + t.effect + epsilons[,2]
  
  
  obs.data = data.frame(obs.A = obs.A, X = X, obs.Y0 = obs.Y0, obs.Y1 = obs.Y1)
  long.data = data.frame(obs.Y = c(obs.Y0, obs.Y1), obs.A = rep(obs.A,2), X1 = rep(X[,1], 2),
                         X2 = rep(X[,2], 2), time = c(rep(0, N), rep(1, N)), 
                         id = rep(1:N, 2))
  long.data$Z2 = ifelse(long.data$time > 0 & long.data$obs.A == 2, 1, 0) # an indicator of direct interventions
  long.data$Z3 = ifelse(long.data$time > 0 & long.data$obs.A == 3, 1, 0) # an indicator of indirect interventions
  long.data$obs.A2 = ifelse(long.data$obs.A == 2, 1, 0) # an indicator of being in a treatment group
  long.data$obs.A3 = ifelse(long.data$obs.A == 3, 1, 0) # an indicator of being in a neighboring control group
  
  ## naive difference-in-differences method
  diff_control = mean(obs.Y1[obs.A == 1]) - mean(obs.Y0[obs.A == 1]) ## 1.control
  diff_treat = mean(obs.Y1[obs.A == 2]) - mean(obs.Y0[obs.A == 2]) ## 2.treatment
  diff_neighbor = mean(obs.Y1[obs.A == 3]) - mean(obs.Y0[obs.A == 3]) ## 3.neighbor
  naive.treat[r] = diff_treat - diff_control
  naive.neighbor[r] = diff_neighbor - diff_control
  
  # fit a multinomial propensity score model (correctly specified)
  fit = multinom(obs.A ~ X)
  pihat_control = fitted.values(fit)[,1]
  pihat_treat = fitted.values(fit)[,2]
  pihat_neighbor = fitted.values(fit)[,3]
  
  # fit a multinomial propensity score models (misspecified)
  #fit2 = multinom(obs.A ~ exp(X[,2]^2) )
  #pihat_control = fitted.values(fit2)[,1]
  #pihat_treat = fitted.values(fit2)[,2]
  #pihat_neighbor = fitted.values(fit2)[,3]  
  
  # fit the outcome model (correctly specified)
  fit.outcome = lmer(obs.Y ~ obs.A2 + obs.A3 + time + 
                       X1 + X2 + 
                       time*X1 + time*X2 + 
                       Z2 + Z3 +  
                       Z2*X1 + Z2*X2  + 
                       Z3*X1 + Z3*X2 + Z3*I(X2^2) + (1|id), data = long.data)
  
  # fit the outcome model (misspecified)
 # fit.outcome = lmer(obs.Y ~ obs.A2 + obs.A3 + time + 
 #                       X1 + X2 + 
 #                       #time*X1 + time*X2 + 
 #                       Z2 + Z3 +  
 #                       Z2*X1 + Z2*X2  + #Z3*I(X2^2)
 #                       Z3*X1 + Z3*X2 + (1|id), data = long.data)
  
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
  
  ## ATT estimation
  # doubly robust estimator
  dr.ATT[r] =  mean(((obs.A == 2) - pihat_treat*(obs.A == 1)/pihat_control)*((obs.Y1 - obs.Y0) - (muhat_control1 - muhat_control0)))/mean(obs.A == 2)
  # regression-based estimator
  reg.ATT[r] = mean((obs.A == 2)*(obs.Y1 - obs.Y0))/mean(obs.A == 2) -  
    mean((obs.A == 2)*(muhat_control1 - muhat_control0))/mean(obs.A == 2)
  # ipw estimator
  ipw.ATT[r] =  mean((obs.Y1 - obs.Y0)*((obs.A==2)/pihat_treat - (obs.A==1)/pihat_control)*pihat_treat)/mean(obs.A == 2)
  # efficient influence function of the doubly robust estimator  
  phi.ATT =  ((obs.A == 2)*((obs.Y1 - obs.Y0) - (muhat_treat1 - muhat_treat0)))/mean(obs.A == 2) - 
    (pihat_treat*(obs.A == 1)/pihat_control*( (obs.Y1 - obs.Y0) - (muhat_control1 - muhat_control0)))/mean(obs.A == 2) + 
    (obs.A==2)*(muhat_treat1 - muhat_treat0 - (muhat_control1 - muhat_control0))/mean(obs.A == 2) -   (obs.A==2)*dr.ATT[r]/mean(obs.A==2)                          

  dr.ATT.var[r] = var(phi.ATT)/N # variance of the influence function
  
  cr.ATT[r] = (dr.ATT[r] - 1.96*sqrt(dr.ATT.var[r]) <= true.ATT & dr.ATT[r] + 1.96*sqrt(dr.ATT.var[r]) >= true.ATT) # 95% confidence interval coverage indicator
  
  ## ATN estimation 
  # doubly robust estimator
  dr.ATN[r] = mean(((obs.A == 3) - pihat_neighbor*(obs.A == 1)/pihat_control)*( (obs.Y1 - obs.Y0) - (muhat_control1 - muhat_control0)))/mean(obs.A == 3)
  # regression-based estimator
  reg.ATN[r] = mean((obs.A == 3)*(obs.Y1 - obs.Y0)) /mean(obs.A == 3) -  
    mean((obs.A == 3)*(muhat_control1 - muhat_control0))/mean(obs.A == 3)
  # ipw estimator
  ipw.ATN[r] = mean((obs.Y1 - obs.Y0)*((obs.A==3)/pihat_neighbor - (obs.A==1)/pihat_control)*pihat_neighbor)/mean(obs.A == 3)
  # efficient influence function of the doubly robust estimator
  phi.ATN = ((obs.A == 3)*( (obs.Y1 - obs.Y0) - (muhat_neighbor1 - muhat_neighbor0)))/mean(obs.A == 3) - 
    (pihat_neighbor*(obs.A == 1)/pihat_control*( (obs.Y1 - obs.Y0) - (muhat_control1 - muhat_control0)))/mean(obs.A == 3) + 
    (obs.A==3)*(muhat_neighbor1 - muhat_neighbor0 - (muhat_control1 - muhat_control0))/mean(obs.A == 3) -   (obs.A==3)*dr.ATN[r]/mean(obs.A==3) 

  dr.ATN.var[r] = var(phi.ATN)/N # variance of the influence function
  
  cr.ATN[r] = (dr.ATN[r] - 1.96*sqrt(dr.ATN.var[r]) <= true.ATN & dr.ATN[r] + 1.96*sqrt(dr.ATN.var[r]) >= true.ATN ) # 95% confidence interval coverage indicator
  
  
  ## -offsetting effect estimation
  a2 = mean(( (obs.A == 3)*pihat_treat/pihat_neighbor)*((obs.Y1 - obs.Y0 -  (muhat_neighbor1 - muhat_neighbor0) ))) + mean((obs.A==2)*(muhat_neighbor1 - muhat_neighbor0)) -  
    mean(( (obs.A == 1)*pihat_treat/pihat_control)*((obs.Y1 - obs.Y0 -  (muhat_control1 - muhat_control0) ))) -  mean((obs.A==2)*(muhat_control1 - muhat_control0))
  dr.offset[r] = a2/ mean(obs.A == 2)
  # regression-based estimator
  reg.offset[r] =  mean((obs.A == 2)*(muhat_neighbor1 - muhat_neighbor0))/mean(obs.A == 2) -  
    mean((obs.A == 2)*(muhat_control1 - muhat_control0))/mean(obs.A == 2)
  # ipw estimator
  ipw.offset[r] = mean((obs.Y1 - obs.Y0)*((obs.A == 3)/pihat_neighbor - (obs.A == 1)/pihat_control)*pihat_treat)/mean(obs.A == 2)
  # efficient influence function of the doubly robust estimator
  phi.offset = ( ((obs.A == 3)*pihat_treat/pihat_neighbor)*((obs.Y1 - obs.Y0 -  (muhat_neighbor1 - muhat_neighbor0) )) + (obs.A==2)*(muhat_neighbor1 - muhat_neighbor0) - 
                   ((obs.A == 1)*pihat_treat/pihat_control)*((obs.Y1 - obs.Y0 -  (muhat_control1 - muhat_control0) )) - (obs.A==2)*(muhat_control1 - muhat_control0))/mean(obs.A == 2) - 
    dr.offset[r]
  dr.offset.var[r] = var(phi.offset)/N # variance of the influence function
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


print(xtable(bias.table, digits = 2))

print(xtable(100*cr.table, digits = 1))

print(se.table)

print(round(c(mean(dr.AOTT - true.AOTT),  mean(sqrt(dr.AOTT.var/N)), mean(cr.AOTT) ),3))
