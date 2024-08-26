library(MASS)
library(nnet)
library(lme4)
library(xtable)
library(plotrix)
###

long.data = read.csv("Data/sample_sales.csv", header = T, sep =",") # read the sample data 
long.data$obs.Y = long.data$Units_TaxInd # (1) taxed individual size beverages
#long.data$obs.Y = long.data$Units_TaxedFam # (2) taxed family size beverages
long.data$obs.A2 = ifelse(long.data$city == "Philadelphia", 1, 0) # treatment group indicator
long.data$obs.A3 = ifelse(long.data$city == "Border Counties", 1, 0) # neighboring control group indicator
long.data$time = ifelse(long.data$year == 2017, 1, 0) # post-treatment time indicator
long.data$obs.Z2 =(long.data$obs.A2)*(long.data$time) # (direct) treatment indicator
long.data$obs.Z3 = (long.data$obs.A3)*(long.data$time) # spillover effect indicator

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

# print out the fitting result (Table S1 in the Supplementary Materials)
z = summary(fit.outcome)$coefficients[,3]
p = (1 - pnorm(abs(z), 0, 1)) * 2
print(xtable(cbind(z,p), digits = 3))

# predicted outcomes conditional on Baltimore
long.control = long.data
long.control$obs.Z2 = 0; long.control$obs.Z3 = 0; long.control$obs.A2 = 0; long.control$obs.A3 = 0
muhat_control = predict(fit.outcome, newdata = long.control, re.form = ~ (1|ID), allow.new.levels = TRUE)
muhat_control0 = muhat_control[long.data$year == 2016]; 
muhat_control1 = muhat_control[long.data$year == 2017]

# predicted outcomes conditional on neighboring controls
long.neighbor = long.data
long.neighbor$obs.Z2 = 0; long.neighbor$obs.Z3 = ifelse(long.neighbor$time == 1, 1, 0); long.neighbor$obs.A2 = 0; long.neighbor$obs.A3 = 1
muhat_neighbor = predict(fit.outcome, newdata = long.neighbor, re.form = ~ (1|ID), allow.new.levels = TRUE)
muhat_neighbor0 = muhat_neighbor[long.data$year == 2016]; 
muhat_neighbor1 = muhat_neighbor[long.data$year == 2017]

# predicted outcomes conditional on Philadelphia
long.treat = long.data;
long.treat$obs.Z2 = ifelse(long.treat$time == 1, 1, 0); long.treat$obs.Z3 = 0; long.treat$obs.A2 = 1; long.treat$obs.A3 = 0
muhat_treat = predict(fit.outcome, newdata = long.treat, re.form = ~ (1|ID), allow.new.levels = TRUE)
muhat_treat0 = muhat_treat[long.data$year == 2016]; 
muhat_treat1 = muhat_treat[long.data$year == 2017]

## short data
short.data = long.data[long.data$year == 2016,]
short.data$obs.Y0 = long.data$obs.Y[long.data$year == 2016]
short.data$obs.Y1 = long.data$obs.Y[long.data$year == 2017]
short.data$obs.A = ifelse(short.data$obs.A2 == 1, 2, ifelse(short.data$obs.A3 == 1, 3, 1))
short.data$obs.D = ifelse(short.data$obs.A == 1, 0, 1)
## 
long.data$obs.D = rep(short.data$obs.D, each = 2)

## propensity score model (two-stage)
cluster.data = aggregate(cbind(family.SS, family.AS, ind.SS, ind.AS, 
                               BlackProportion, HispanicProportion, AsianProportion,
                               MaleProportion, IncomePerHousehold,  obs.D) ~ cm, data = short.data, mean)
names(cluster.data) = c("cm", "family.SS.mean", "family.AS.mean", "ind.SS.mean", "ind.AS.mean", 
                        "BlackProportion.mean", "HispanicProportion.mean", "AsianProportion.mean",
                        "MaleProportion.mean", "IncomePerHousehold.mean", "obs.D")
subind.data = short.data[short.data$obs.A!=1, ]
subind.data$treat = (subind.data$obs.A == 2)


fit1 = glm(obs.D ~ family.SS.mean + family.AS.mean + ind.SS.mean + ind.AS.mean +
             BlackProportion.mean + HispanicProportion.mean + AsianProportion.mean + 
             MaleProportion.mean + IncomePerHousehold.mean, data = cluster.data, family = binomial())
fit2 = glm(treat ~ family.SS + family.AS + ind.SS + ind.AS + 
             BlackProportion + HispanicProportion + AsianProportion + 
             MaleProportion + IncomePerHousehold, data = subind.data, family = binomial())
cluster.prop = rep(NA, nrow(short.data))
for(j in 1:nrow(cluster.data)) cluster.prop[short.data$cm %in% cluster.data$cm[j]] = fitted.values(fit1)[j]
#cluster.prop = rep(fitted.values(fit1), each = n)
ind.prop = predict(fit2, newdata = short.data, type = "response") 
pihat_control = (1-cluster.prop)
pihat_treat = cluster.prop*ind.prop
pihat_neighbor = cluster.prop*(1-ind.prop)


Units_TaxInd_results = real_sim_function(long.data, theta = 1)

boot.ipw.ATT = boot.ipw.ATN = boot.ipw.offset = boot.ipw.AOTT = 
  boot.reg.ATT = boot.reg.ATN = boot.reg.offset = boot.reg.AOTT = c()
boot.dr.ATT = boot.dr.ATN = dr.ATT.var = dr.ATN.var = c();
boot.dr.offset = boot.dr.AOTT = matrix(NA, 500, 27)
I = nrow(cluster.data)

for(b in 1:500){
  set.seed(b)
  print(b)
  sampled = sample(1:I, I, replace = TRUE); newdata = c()
  for(k in 1:I) newdata = rbind(newdata, long.data[long.data$cm %in% sampled[k] ,])
  
  results = real_sim_function(newdata)
  
  boot.ipw.ATT[b] = results[[1]][[1]]
  boot.ipw.ATN[b] = results[[2]][[1]]
  boot.ipw.offset[b] = results[[3]][[1]]
  boot.ipw.AOTT[b] = results[[4]][[1]]
  
  boot.reg.ATT[b] = results[[1]][[2]]
  boot.reg.ATN[b] = results[[2]][[2]]
  boot.reg.offset[b] = results[[3]][[2]]
  boot.reg.AOTT[b] = results[[4]][[2]]
  
  boot.dr.ATT[b] = results[[1]][[3]]
  boot.dr.ATN[b] = results[[2]][[3]]
  boot.dr.offset[b,] = results[[3]][[3]]
  boot.dr.AOTT[b,] = results[[4]][[3]]
  
}

## ATT results
c(Units_TaxInd_results[[1]][[1]], Units_TaxInd_results[[1]][[1]] - 1.96*sd(boot.ipw.ATT), Units_TaxInd_results[[1]][[1]] + 1.96*sd(boot.ipw.ATT))
c(Units_TaxInd_results[[1]][[2]], Units_TaxInd_results[[1]][[2]] - 1.96*sd(boot.reg.ATT), Units_TaxInd_results[[1]][[2]] + 1.96*sd(boot.reg.ATT))
#dr.ATT.results = c(Units_TaxInd_results[[1]][[3]], Units_TaxInd_results[[1]][[3]] - 1.96*sd(boot.dr.ATT), Units_TaxInd_results[[1]][[3]] + 1.96*sd(boot.dr.ATT))
dr.ATT.results = c(Units_TaxInd_results[[1]][[3]], quantile(boot.dr.ATT, 0.025), quantile(boot.dr.ATT, 0.975))
c(Units_TaxInd_results[[1]][[3]], Units_TaxInd_results[[1]][[3]] - 1.96*sqrt(Units_TaxInd_results[[1]][[4]]), Units_TaxInd_results[[1]][[3]] + 1.96*sqrt(Units_TaxInd_results[[1]][[4]]))


## ATN results
c(Units_TaxInd_results[[2]][[1]], Units_TaxInd_results[[2]][[1]] - 1.96*sd(boot.ipw.ATN), Units_TaxInd_results[[2]][[1]] + 1.96*sd(boot.ipw.ATN))
c(Units_TaxInd_results[[2]][[2]], Units_TaxInd_results[[2]][[2]] - 1.96*sd(boot.reg.ATN), Units_TaxInd_results[[2]][[2]] + 1.96*sd(boot.reg.ATN))
dr.ATN.results = c(Units_TaxInd_results[[2]][[3]], quantile(boot.dr.ATN, 0.025), quantile(boot.dr.ATN, 0.975))
c(Units_TaxInd_results[[2]][[3]], Units_TaxInd_results[[2]][[3]] - 1.96*sqrt(Units_TaxInd_results[[2]][[4]]), Units_TaxInd_results[[2]][[3]] + 1.96*sqrt(Units_TaxInd_results[[2]][[4]]))

## offsetting results
c(Units_TaxInd_results[[3]][[1]], Units_TaxInd_results[[3]][[1]] - 1.96*sd(boot.ipw.offset), Units_TaxInd_results[[3]][[1]] + 1.96*sd(boot.ipw.offset))
c(Units_TaxInd_results[[3]][[2]], Units_TaxInd_results[[3]][[2]] - 1.96*sd(boot.reg.offset), Units_TaxInd_results[[3]][[2]] + 1.96*sd(boot.reg.offset))
dr.offset.results = cbind(Units_TaxInd_results[[3]][[3]], apply(boot.dr.offset,2, quantile, 0.025), apply(boot.dr.offset, 2, quantile, 0.975))

## AOTT results
c(Units_TaxInd_results[[4]][[1]], Units_TaxInd_results[[4]][[1]] - 1.96*sd(boot.ipw.AOTT), Units_TaxInd_results[[4]][[1]] + 1.96*sd(boot.ipw.AOTT))
c(Units_TaxInd_results[[4]][[2]], Units_TaxInd_results[[4]][[2]] - 1.96*sd(boot.reg.AOTT), Units_TaxInd_results[[4]][[2]] + 1.96*sd(boot.reg.AOTT))
dr.AOTT.results = cbind(Units_TaxInd_results[[4]][[3]],  apply(boot.dr.AOTT,2,quantile, 0.025), apply(boot.dr.AOTT, 2, quantile, 0.975))

## Out of 27, we only consider three conditionally doubly robust estimator with gamma1 = gamma2
index = c(1, 2, 3)


## Figure 3
par(mfrow = c(1,1), cex.lab = 1.7, cex.axis = 1.5, 
    mar=c(5,5,3,5), tcl = 0.5,  xpd = FALSE, cex.main = 2)
plotCI(x = c(1:5),              
       y = c(dr.ATT.results[1], dr.ATN.results[1], dr.AOTT.results[index,1]),
       li = c(dr.ATT.results[2], dr.ATN.results[2], dr.AOTT.results[index,2]),
       ui = c(dr.ATT.results[3], dr.ATN.results[3], dr.AOTT.results[index,3]), xaxt = "n",
       ylab = "Estimates", xlab = "Target Estimand", main = "(1) Taxed individual sized beverages",
       #main = "(2) Taxed family sized beverages",
       col = c("red", "dodgerblue", rep("black", 9)), lwd = 2, pch = 19, cex = 2)
abline(h = 0)
axis(1, at = c(1,2,3,4,5), tick = TRUE, c("ATT", "ATN", "AOTT", "AOTT", "AOTT"), 
     col = "black")
text(c(3,4,5), -2000, labels = c(1,2,3), 
     col = "black")