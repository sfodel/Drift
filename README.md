# Drift
Simulations of random ecological drift in microbial metacommunities under different selection and dispersal regimes


These simulations are performed using parameters for the frequency and the magnitude of drift that are derived experimentally. The respective experiments were performed using fluorescence-independent flow cytometry in synthetic microbial communities of up to three bacterial strains. For more details on the experimental methodology please see here: https://www.nature.com/articles/s41467-019-09306-7 

The code used for these simulations is an attempt to generalize our experimental findings for natural bacterial communities. 

In our simulations, 100 communities of 2000 species whose starting populations follow a log-normal distribution (mean = 4, sd=1.1) grow for 1000 generations under varying degrees of dispersal and selection and under the presence or absence of drift. 

At each generation, each species in each community grows with a rate drawn from a normal distribution with a mean of 1 and with a standard deviation that is analogous to the strength of selection (here, 0.071 < sd <= 0.167). That rate changes further according to the experimentally-derived magnitude of drift in 388 randomly chosen species (19.4% of the total 2000 - that is the experimentally-observed drift frequency).

Additionally, dispersal occurs among communities in a circular, unidirectional way and the percentage of total cells that migrate from community to community at each generation represents the rate of dispersal (here, between 2% and 20%).

Communities grow in parallel both with and without drift, and the effect of drift on fundamental ecological community properties (e.g. species richness, community evennes community size) is derived from comparisons between the two. 

Here, we used this model to examine the effect of drift in 50 different scenarios resulting from five levels of selection strength over ten levels of dispersal rate, but the user can provide any positive value for selection and any positive value between 0 and 1 for dispersal. Moreover, the user can change the number of communities and the number of species in the metacommunity, as well as the number of generations.

Finally, the modified script "Drift_simulation_as_function_with_RC.R" runs the same growth model but instead calculates the proportion of community pairs under dispersal limitation every 100 generations, using the script  of Stegen et al. (https://github.com/stegen/Stegen_etal_ISME_2013/blob/master/Raup_Crick_Abundance.r) that is based on Bray-Curtis-based Raup-Crick distances. Careful, it is extremely slow for metacommunities over 500 species and 50 communities.
