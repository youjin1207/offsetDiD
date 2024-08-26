library(MASS)
library(nnet)
library(lme4)
library(xtable)
source("Code/sim_function.R")
###

## parameters ##
betas_cluster = c(-0.5, -1.0); betas_ind = c(-1.0, 1.0)
gamma_0 = 1.0; gamma_1 = 1.5
lambda_0 = c(0.5, 0.5); lambda_1 = c(-0.5, 1.0)
alpha_a = c(0.0,0.5,-1.0) # group effect: control, treatment, neighbor
theta_neighbor = 0.5; theta_treat = -1


true.AOTT = -0.4393935
#######
n.rep = 500
n = 100
nb = 100
#######

dr.AOTT = boot.dr.se = dr.AOTT.bootcr = cboot.dr.se = dr.AOTT.cbootcr = matrix(NA, n.rep, 27)
ipw.AOTT = reg.AOTT = boot.ipw.se = boot.reg.se = ipw.AOTT.bootcr = reg.AOTT.bootcr = 
  cboot.ipw.se = cboot.reg.se = ipw.AOTT.cbootcr = reg.AOTT.cbootcr = c()

############### Simulation Study ##################
for(r in 1:n.rep){
  
  I = 100; N = I*n
  cm = rep(1:I, each = n) # cluster memberships
  epsilons = mvrnorm(N, c(0,0), matrix(c(1, 0.2, 0.2, 1), 2, 2)) # correlated within unit outcomes (two time points)
  X = matrix(NA, N, 2)
  for(i in 1:N) X[i,] = mvrnorm(1, c(0.01*(cm[i]-50),-0.01*(cm[i]-50)), matrix(c(1, 0.3, 0.3, 1), 2, 2))
  
  # level 1: cluster allocations 
  cluster_X = aggregate(X ~ cm, data = cbind(X, cm), mean)[,2:3]
  cluster_assignment = rbinom(I, 1, plogis(betas_cluster %*% t(cluster_X) )) # betas_cluster = c(-0.5 -1.0)
  D = rep(cluster_assignment, each = n)
  
  # level 2: individual allocations 
  Z = rep(NA, N)
  Z[D==0] = 0
  Z[D==1] = rbinom(sum(D==1), 1, plogis(betas_ind %*% t(X[D==1,]))) # betas_ind = c(-1.0, 1.0)
  
  obs.A = ifelse(D == 0, 1, ifelse(Z == 1, 2, 3))
  ## generate outcome at t=0 (pre-treatment)
  obs.Y0 = gamma_0 + X %*% lambda_0 + alpha_a[obs.A] + epsilons[,1]
  ## generate the treatment effects
  t.effect = ifelse(obs.A == 1, 0, ifelse(obs.A == 2, theta_treat*(1 - 1.0*X[,1]),
                                          ifelse(obs.A == 3, theta_neighbor*(1 - 2.0*X[, 2]), NA)))
  # generate the outcome at t=1 (post-treatment)
  obs.Y1 = gamma_1 + X %*% lambda_1 + alpha_a[obs.A] + t.effect + epsilons[,2]
  
  obs.data = data.frame(obs.A = obs.A, X = X, obs.Y0 = obs.Y0, obs.Y1 = obs.Y1, cm = cm, D = D)
  
  results0 = sim_function(obs.data = obs.data, case = 1) ## case = 1, 2, 3, 4
  
  ipw.AOTT[r] = results0[[4]]
  reg.AOTT[r] = results0[[5]]
  dr.AOTT[r,] = results0[[6]]
  
  ## cluster bootstrap
  cboot.ipw = cboot.reg = cboot.dr = c(); 
  for(b in 1:nb){
    set.seed(b)
    sampled = sample(1:I, I, replace = TRUE); newdata = c()
    for(k in 1:I) newdata = rbind(newdata, obs.data[obs.data$cm %in% sampled[k] ,])
    newdata$cm = rep(1:I, each = n)
    results = sim_function(obs.data = newdata, case = 1) ## case = 1, 2, 3, 4
    cboot.ipw = c(cboot.ipw, results[[4]])
    cboot.reg = c(cboot.reg, results[[5]])
    cboot.dr = rbind(cboot.dr, results[[6]])
  }
  
  cboot.ipw.se[r] = sd(cboot.ipw)
  cboot.reg.se[r] = sd(cboot.reg)
  cboot.dr.se[r,] = apply(cboot.dr, 2, sd)
  
  ipw.AOTT.cbootcr[r] = (quantile(cboot.ipw, 0.025) <= mean(true.AOTT) & quantile(cboot.ipw, 0.975) >= mean(true.AOTT) ) 
  reg.AOTT.cbootcr[r] = (quantile(cboot.reg, 0.025) <= mean(true.AOTT) &  quantile(cboot.reg, 0.975) >= mean(true.AOTT) )
  dr.AOTT.cbootcr[r,] = (apply(cboot.dr,2, quantile, 0.025) <= mean(true.AOTT) & apply(cboot.dr,2, quantile, 0.975) >= mean(true.AOTT) ) 
  
}