# Policy effect evaluation under counterfactual neighborhood intervention in the presence of spillover

## Data

We use the data obtained from major US retailers, which was purchased from Information Resources Inc (IRI). We focus on beverages prices and sales in Philadelphia, the counties within approximately 3 miles of Philadelphia's border, and Baltimore.
The details can be found [here](https://jamanetwork.com/journals/jama/fullarticle/2733208) and [here](https://ageconsearch.umn.edu/record/234905/).

We provide the sample data `Data/sample sales.csv`, which has the same data structure as the Philadelphia beverage data we evaluated in the manuscript (e.g., panel data with several time points from three separate treatment groups). This data was created at random and is only intended to be used as an example, not to be replicated.

## Code

- `Code/sim_main.R` : generates the main simulation studies with two time points and derives the doubly robust estimates for diverse causal effects in the presence of spillover.

- `Code/sim_multi.R`: generate supplementary simulation studies with multiple time points.  

- `Code/realdata.R`:  reads the real data example and estimates diverse causal effects of the policy intervention.

 

## Instructions for the use of sample data

In `Data/sample.csv`, we provide a hypothetical data with the same data structure as the real data presented in the paper. The data is provided for illustrative purpose only, not for reproducing the results. 

```{r}

```
