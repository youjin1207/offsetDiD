library(MASS)
library(nnet)
library(lme4)
library(xtable)
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

# print out the fitting result (Table S5 in the Supplementary Materials)
z = summary(fit.outcome)$coefficients[,3]
p = (1 - pnorm(abs(z), 0, 1)) * 2
print(xtable(cbind(z,p), digits = 3))

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
                      MaleProportion + IncomePerHousehold, data = short.data)
pihat_control = fitted.values(prop.fit)[,1]
pihat_treat = fitted.values(prop.fit)[,2]
pihat_neighbor = fitted.values(prop.fit)[,3]

# print out the fitting result (Table S4 in the Supplementary Materials)
z = summary(prop.fit)$coefficients/summary(prop.fit)$standard.errors
print(xtable(z, digits = 3))
p = (1 - pnorm(abs(z), 0, 1)) * 2
print(xtable(p, digits = 3))

###############################################################
### Reproduce the results in Table 5 in the main manuscript ###
###############################################################

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
print(c(dr.ATT - 1.96*sqrt(dr.ATT.var),  dr.ATT + 1.96*sqrt(dr.ATT.var)))


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
print(c(dr.ATN - 1.96*sqrt(dr.ATN.var),  dr.ATN + 1.96*sqrt(dr.ATN.var)))

## -offsetting effect (delta) ##
a2 = mean(((short.data$obs.A == 3)*pihat_treat/pihat_neighbor)*((short.data$obs.Y1 - short.data$obs.Y0 -  (muhat_neighbor1 - muhat_neighbor0) ))) + mean((short.data$obs.A==2)*(muhat_neighbor1 - muhat_neighbor0)) -  
  mean((control.A*pihat_treat/pihat_control)*((short.data$obs.Y1 - short.data$obs.Y0 -  (muhat_control1 - muhat_control0) ))) -  mean((short.data$obs.A==2)*(muhat_control1 - muhat_control0))
dr.offset = a2/mean(short.data$obs.A == 2) # doubly robust estimator
# regression-based estimator
reg.offset =  mean((short.data$obs.A == 2)*(muhat_neighbor1 - muhat_neighbor0))/mean(short.data$obs.A == 2) -  
  mean((short.data$obs.A == 2)*(muhat_control1 - muhat_control0))/mean(short.data$obs.A == 2)
# ipw estimator
ipw.offset = mean((short.data$obs.Y1 - short.data$obs.Y0)*((short.data$obs.A == 3)/pihat_neighbor - control.A/pihat_control)*pihat_treat)/mean(short.data$obs.A == 2)
# estimated efficient influence function
phi.offset = (pihat_treat*(short.data$obs.A == 3)/pihat_neighbor*( (short.data$obs.Y1 - short.data$obs.Y0) - (muhat_neighbor1 - muhat_neighbor0)))/mean(short.data$obs.A == 2) - 
  (pihat_treat*(short.data$obs.A == 1)/pihat_control*( (short.data$obs.Y1 - short.data$obs.Y0) - (muhat_control1 - muhat_control0)))/mean(short.data$obs.A == 2) + 
  (short.data$obs.A==2)*(muhat_neighbor1 - muhat_neighbor0 - (muhat_control1 - muhat_control0))/mean(short.data$obs.A == 2) - (short.data$obs.A==2)*dr.offset/mean(short.data$obs.A==2) 
dr.offset.var = mean(phi.offset^2)/nrow(short.data)
print(c(dr.offset - 1.96*sqrt(dr.offset.var),  dr.offset + 1.96*sqrt(dr.offset.var)))

## AOTT ##
reg.AOTT = reg.ATT + reg.offset
dr.AOTT = dr.ATT + dr.offset
ipw.AOTT = ipw.ATT + ipw.offset
dr.AOTT.var = var(phi.ATT + phi.offset)/nrow(short.data)
print(c(dr.AOTT - 1.96*sqrt(dr.AOTT.var),  dr.AOTT + 1.96*sqrt(dr.AOTT.var)))
print(c(dr.ATT - 1.96*sqrt(dr.ATT.var),  dr.ATT + 1.96*sqrt(dr.ATT.var)))