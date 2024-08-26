# Policy effect evaluation under counterfactual neighborhood intervention in the presence of spillover

## Overview 

Policy interventions can spill over to  units of a population that are not directly exposed to the policy but are geographically close to the units receiving the intervention. In recent work, investigations of spillover effects on neighboring regions have focused on estimating the average treatment effect of a particular policy in an observed setting. Our research question broadens this scope by asking what policy consequences would the treated units have experienced under hypothetical exposure settings. When we only observe treated unit(s) surrounded by controls -- as is common when a policy intervention is implemented in a single city or state -- this effect inquires about the policy effects under a counterfactual neighborhood policy status that we do not, in actuality, observe. 

In this work, we extend difference-in-differences (DiD) approaches to spillover settings and develop identification conditions required to evaluate policy effects in counterfactual treatment scenarios. These causal quantities are policy-relevant for designing effective policies for populations subject to various neighborhood statuses. We develop several estimators that have desirable properties. We apply our proposed method to investigate the effect of the Philadelphia beverage tax on unit sales.

This repository provides the code and sample data that can reproduce the results in the manuscript.

## Data

We use the data obtained from major US retailers, which was purchased from Information Resources Inc (IRI). We focus on beverages prices and sales in Philadelphia, the counties within approximately 3 miles of Philadelphia's border, and Baltimore.
The details can be found [here](https://jamanetwork.com/journals/jama/fullarticle/2733208) and [here](https://ageconsearch.umn.edu/record/234905/).

We provide the sample data `Data/sample_sales.csv`, which has the same data structure as the Philadelphia beverage data we evaluated in the manuscript (e.g., panel data with several time points from three separate treatment groups). This data was created at random and is only intended to be used as an example, not to be replicated.

## Code

- `Code/sim_function.R`: is an auxiliary function to estimate the offseting effect using (i) the inverse probability weighted (IPW) estimator, (ii) outcome regression model, and (iii) conditionally doubly robust estimators

- `Code/sim_main.R` : generates the main simulation studies with two time points and derives the doubly robust estimates for the offsetting effect. This `R` file can be used to reproduce ``Figure 4`` in the main manuscript.

- `Code/real_sim_function.R`: is an auxiliary function to estimate the average treatment effect on the treated (ATT), the average treatment effect on the neighboring controls (ATN), the offsetting effect, and the average offset treatment effect on the treated (AOTT). 

- `Code/realdata.R`:  reads the real data example and estimates diverse causal effects of the policy intervention. This `R` file can be used to reproduce ``Figure 5`` in the main manuscript, but with the pseudo data. 

 
## Instructions for the use of sample data

In `Data/sample.csv`, we provide a hypothetical data with the same data structure as the real data presented in the paper. The data is provided for illustrative purpose only, not for reproducing the results. To access the complete code for working with the real data, refer to ``Code/realdata.R``.

```{r}
library(MASS)
library(nnet)
library(lme4)
library(xtable)

long.data = read.csv("Data/sample_sales.csv", header = T, sep =",") # read the sample data 
```


