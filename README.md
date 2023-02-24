# Policy effect evaluation under counterfactual neighborhood intervention in the presence of spillover

## Overview 

Policy interventions can spill over to  units of a population that are not directly exposed to the policy but are geographically close to the units receiving the intervention. In recent work, investigations of spillover effects on neighboring regions have focused on estimating the average treatment effect of a particular policy in an observed setting. Our research question broadens this scope by asking what policy consequences would the treated units have experienced under hypothetical exposure settings. When we only observe treated unit(s) surrounded by controls -- as is common when a policy intervention is implemented in a single city or state -- this effect inquires about the policy effects under a counterfactual neighborhood policy status that we do not, in actuality, observe. 

In this work, we extend difference-in-differences (DiD) approaches to spillover settings and develop identification conditions required to evaluate policy effects in counterfactual treatment scenarios. These causal quantities are policy-relevant for designing effective policies for populations subject to various neighborhood statuses. 

We develop doubly robust estimators and use extensive numerical experiments to examine their performance under heterogeneous spillover effects. We apply our proposed method to investigate the effect of the Philadelphia beverage tax on volume sales.

This repository provides the code and sample data that can reproduce the results in the manuscript.

## Data

We use the data obtained from major US retailers, which was purchased from Information Resources Inc (IRI). We focus on beverages prices and sales in Philadelphia, the counties within approximately 3 miles of Philadelphia's border, and Baltimore.
The details can be found [here](https://jamanetwork.com/journals/jama/fullarticle/2733208) and [here](https://ageconsearch.umn.edu/record/234905/).

We provide the sample data `Data/sample_sales.csv`, which has the same data structure as the Philadelphia beverage data we evaluated in the manuscript (e.g., panel data with several time points from three separate treatment groups). This data was created at random and is only intended to be used as an example, not to be replicated.

## Code

- `Code/sim_main.R` : generates the main simulation studies with two time points and derives the doubly robust estimates for diverse causal effects in the presence of spillover. This `R` file can be used to reproduce ``Table 3`` in the main manuscript and ``Table S1`` in the Supplementary Materials.

- `Code/sim_multi.R`: generate supplementary simulation studies with multiple time points. This `R` file can be used to reproduce ``Table S2`` and ``Table S3`` in the Supplementary Materials.  

- `Code/realdata.R`:  reads the real data example and estimates diverse causal effects of the policy intervention. This `R` file can be used to reproduce ``Table 5`` in the main manuscript and ``Table S4`` and ``Table S5`` in the Supplementary Materials, but with the pseudo data. 

 
## Instructions for the use of sample data

In `Data/sample.csv`, we provide a hypothetical data with the same data structure as the real data presented in the paper. The data is provided for illustrative purpose only, not for reproducing the results. To access the complete code for working with the real data, refer to ``Code/realdata.R``.

```{r}
library(MASS)
library(nnet)
library(lme4)
library(xtable)

long.data = read.csv("Data/sample_sales.csv", header = T, sep =",") # read the sample data 
```

We first define the outcome of interests: (1) taxed individual sized beverages (``$Units_TaxInd``) or (2) taxed family sized beverages (``$Units_TaxedFam``). 
```{r}
long.data$obs.Y = long.data$Units_TaxInd #(1) taxed individual size beverages
## long.data$obs.Y = long.data$Units_TaxedFam #(2) taxed family size beverages
```

The ``$city`` of each store is used to defined different treatment groups.
The variable ``$obs.A2`` is an indicator of the treatment group, i.e., Philadelphia, and the variable ``$obs.A3`` is an indicator of the neighboring control group, i.e., neighboring counties of Philadelphia. The variable ``$time`` indicates the post-intervention period. In the Philadelphia beverage tax study, ``year = 2017`` means ``t=1``. Then the treatment status ``$obs.Z2`` denotes the direct exposure to the intervention and ``$obs.Z3`` denotes the indirect exposure to the intervention, each of which is defined by the interaction between the group and time.

```{r}
long.data$obs.A2 = ifelse(long.data$city == "Philadelphia", 1, 0) # treatment group indicator
long.data$obs.A3 = ifelse(long.data$city == "Border Counties", 1, 0) # neighboring control group indicator
long.data$time = ifelse(long.data$year == 2017, 1, 0) # post-treatment time indicator
long.data$obs.Z2 =(long.data$obs.A2)*(long.data$time) # (direct) treatment indicator
long.data$obs.Z3 = (long.data$obs.A3)*(long.data$time) # spillover effect indicator
```

Then the structure of the working data is as follows: 
```{r}
> head(long.data)

  ID         city  PeriodEnd year Units_TaxInd Units_TaxedFam Units_NotTaxInd
1  1 Philadelphia 01/31/2016 2016     9801.914       23700.62           12015
2  1 Philadelphia 02/28/2016 2016    12338.405       24744.12           11723
3  1 Philadelphia 03/27/2016 2016    11897.121       27900.32           11481
4  1 Philadelphia 04/24/2016 2016    12167.598       24456.35           11989
5  1 Philadelphia 05/22/2016 2016    12825.630       27386.46           14817
6  1 Philadelphia 06/19/2016 2016    16804.318       27702.85           12528
  Units_NotTaxedFam   ind.SS   ind.AS family.SS family.AS IncomePerHousehold MaleProportion
1           27569.0 8.131427 5.391999  3.814364  5.339303          0.3006434      0.5005797
2           26854.0 8.131427 5.391999  3.814364  5.339303          0.3006434      0.5005797
3           27928.0 8.131427 5.391999  3.814364  5.339303          0.3006434      0.5005797
4           26439.0 8.131427 5.391999  3.814364  5.339303          0.3006434      0.5005797
5           26963.0 8.131427 5.391999  3.814364  5.339303          0.3006434      0.5005797
6           27082.1 8.131427 5.391999  3.814364  5.339303          0.3006434      0.5005797
  AsianProportion HispanicProportion BlackProportion     obs.Y obs.A2 obs.A3 time obs.Z2
1       0.1609484          0.1008612      0.07528983  9801.914      1      0    0      0
2       0.1609484          0.1008612      0.07528983 12338.405      1      0    0      0
3       0.1609484          0.1008612      0.07528983 11897.121      1      0    0      0
4       0.1609484          0.1008612      0.07528983 12167.598      1      0    0      0
5       0.1609484          0.1008612      0.07528983 12825.630      1      0    0      0
6       0.1609484          0.1008612      0.07528983 16804.318      1      0    0      0
  obs.Z3
1      0
2      0
3      0
4      0
5      0
6      0
```

In the outcome regression model, we consider the fixed effect of treatment group indicator, the time indicator, and their interactions with store-level random intercepts (``(1|ID)``). We also include covariates effects for the zip-code level characteristics (``BlackProportion``, ``HispanicProportion``, ``AsianProportion``, ``MaleProportion``, ``IncomePerHousehold``) and the baseline weighted price per ounce of both family- and individual-sized sweetened beverages (``family.SS``, ``family.AS``, ``ind.SS``, ``ind.AS``). We include the interaction between price and treatment effects as well (``obs.Z2*ind.SS``, ``obs.Z3*ind.SS``, ``obs.Z2*family.SS``, ``obs.Z3*family.SS``, ``obs.Z2*ind.AS``, ``obs.Z3*ind.AS``, ``obs.Z2*family.AS``, ``obs.Z3*family.AS``).

```{r}
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
```

We only use the baseline data to fit a propensity score model. We include weighted price ounce over the first 4-week period (with the end date of 01/01/2016) of family-sized sweetened beverages and individual sized sweetened beverages. We also include publicly available zip-code level characteristics. 

```{r}
## short data
short.data = long.data[long.data$year == 2016,]
short.data$obs.Y0 = long.data$obs.Y[long.data$year == 2016]
short.data$obs.Y1 = long.data$obs.Y[long.data$year == 2017]
short.data$obs.A = ifelse(short.data$obs.A2 == 1, 2, ifelse(short.data$obs.A3 == 1, 3, 1))

## propensity score model
prop.fit = multinom(obs.A ~ family.SS + family.AS + ind.SS + ind.AS + 
                      BlackProportion + HispanicProportion + AsianProportion + 
                      MaleProportion + IncomePerHousehold, data = short.data)                
```
