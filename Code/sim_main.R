library(MASS)
library(nnet)
library(lme4)
library(xtable)
source("Code/sim_function.R")
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
obs.A = matrix(0, NN, 3)
for(i in 1:NN){
  obs.A[i,] =  rmultinom(1, 1, c(exp(tmp)[i,1:3]))
  eta[i] = exp(tmp[i,3])/exp(tmp[i,2])
}
obs.A = cbind(obs.A, rep(0, NN))
obs.A = apply(obs.A, 1, function(x) which(x == 1)) # 1: control, 2: treatment, 3: neighbors


true.ATT =  mean(theta_treat*(1 - 0.5*X[obs.A==2,1]))
true.ATN =  mean(theta_neighbor*(1 + 0.5*X[obs.A==3, 1]) + theta_neighbor*(1-2*(X[obs.A==3,2])^2))
true.offset = mean((exp(tmp[obs.A == 2,3])/exp(tmp[obs.A == 2,2]))*theta_neighbor*(1 + 0.5*X[obs.A==2, 1]) + theta_neighbor*(1-0.5*(X[obs.A==2,2])^2)) # -(offset effects)
true.AOTT =  true.ATT + true.offset

dr.AOTT = boot.dr.se = dr.AOTT.cr = matrix(NA, 10, 27)
ipw.AOTT = reg.AOTT = boot.ipw.se = boot.reg.se = ipw.AOTT.cr = reg.AOTT.cr =  c()

############### Simulation Study ##################
for(r in 1:10){
  
  set.seed(r)
  
  N = 2000 # N=500, 1000 
  epsilons = mvrnorm(N, c(0,0), matrix(c(1, 0.2, 0.2, 1), 2, 2)) # correlated within unit outcomes (two time points)
  X = mvrnorm(N, c(0,0), matrix(c(1, 0.3, 0.3, 1), 2, 2))

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
                                          ifelse(obs.A == 3, theta_neighbor*(1 + 0.5*X[, 1]) + theta_neighbor*(1-0.5*(X[,2]^2)), NA)))
  # generate the outcome at t=1 (post-treatment)
  obs.Y1 = gamma_1 + X %*% lambda_1 + alpha_a[obs.A] + t.effect + epsilons[,2]
  
  obs.data = data.frame(obs.A = obs.A, X = X, obs.Y0 = obs.Y0, obs.Y1 = obs.Y1)
  
  results0 = sim_function(obs.data = obs.data, case = 1) ## case = 1, 2, 3, 4
  ipw.AOTT[r] = results0[[1]]
  reg.AOTT[r] = results0[[2]]
  dr.AOTT[r,] = results0[[3]]
  
  ipw.AOTT[r] = results0[[4]]
  reg.AOTT[r] = results0[[5]]
  dr.AOTT[r,] = results0[[6]]
  
  ## bootstrap 
  boot.ipw = boot.reg = boot.dr = c()
  for(b in 1:100){
    set.seed(b)
    sampled = sample(1:N, N, replace = TRUE)
    newdata = obs.data[sampled,]
    results = sim_function(obs.data = newdata, case = 1) ## case = 1, 2, 3, 4
    boot.ipw = c(boot.ipw, results[[4]])
    boot.reg = c(boot.reg, results[[5]])
    boot.dr = rbind(boot.dr, results[[6]])
  }
  
  boot.ipw.se[r] = sd(boot.ipw)
  boot.reg.se[r] = sd(boot.reg)
  boot.dr.se[r,] = apply(boot.dr, 2, sd)
  
  ipw.AOTT.cr[r] = (ipw.AOTT[r] - 1.96*boot.ipw.se[r] <= true.AOTT & ipw.AOTT[r] + 1.96*boot.ipw.se[r] >= true.AOTT ) 
  reg.AOTT.cr[r] = (reg.AOTT[r] - 1.96*boot.reg.se[r] <= true.AOTT & reg.AOTT[r] + 1.96*boot.reg.se[r] >= true.AOTT )
  dr.AOTT.cr[r,] = (dr.AOTT[r,] - 1.96*boot.dr.se[r,] <= true.AOTT & dr.AOTT[r,] + 1.96*boot.dr.se[r,] >= true.AOTT ) 

}

result.table = matrix(NA, nrow = 29, ncol= 5)
row.names(result.table) = c("IPW", "Reg", rep("CDR", 27)); colnames(result.table) = c("Bias", "SE", "Coverage rate", "SE (boot)", "Coverage rate (boot)")
result.table[,1] = c(mean(ipw.AOTT), mean(reg.AOTT), colMeans(dr.AOTT)) - true.AOTT
result.table[,2] = c(sd(ipw.AOTT), sd(reg.AOTT), apply(dr.AOTT, 2, sd))
result.table[1:2,3] = c(mean(ipw.AOTT >= true.AOTT - 1.96*sd(ipw.AOTT) &  ipw.AOTT <= true.AOTT + 1.96*sd(ipw.AOTT)), 
                        mean(reg.AOTT >= true.AOTT - 1.96*sd(reg.AOTT) &  reg.AOTT <= true.AOTT + 1.96*sd(reg.AOTT)))
for(k in 1:27){
  result.table[(2+k), 3] = mean(dr.AOTT[,k] >= true.AOTT - 1.96*sd(dr.AOTT[,k]) & dr.AOTT[,k] <= true.AOTT + 1.96*sd(dr.AOTT[,k]))
}

result.table[,4] = c(mean(boot.ipw.se), mean(boot.reg.se), apply(boot.dr.se, 2, mean))
result.table[,5] = c(mean(ipw.AOTT.cr), mean(reg.AOTT.cr), apply(dr.AOTT.cr, 2, mean))

print(result.table)