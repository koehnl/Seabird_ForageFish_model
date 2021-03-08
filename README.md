# Seabird_ForageFish_model
R Code to run models from Koehn et al. 2021 - A structured seabird population model reveals how alternative forage fish control rules benefit seabirds and fisheries. Ecological Applications

TO RUN SEABIRD RUNS WITH A FORAGE FISH PREY
- First run forage fish code Run_ForageFishmodel_Siple_etal_2019.R for both Sardine and Anchovy - saves files to your project directory for each forage fish run under different control rules. Seabird model code for Koehn et al. 2021 will load and use these files.

PLEASE READ COMMENTS AT THE TOP OF EACH R SCRIPT FILE

This repository includes:

-Forage fish model function and associated functions from or modified from Siple et al. 2019. And code to run this model for different forage fish - Run_ForageFishmodel_Siple_etal_2019.R
Siple, M.C., Essington, T.E. and E. Plagányi, É., 2019. Forage fish fisheries management requires a tailored approach to balance trade‐offs. Fish and Fisheries, 20(1), pp.110-124.
  - saves forage fish model runs to be used in the seabird model. 

-Seabird model function (in Seabird model folder) - both a non-stochastic model and model with stochasticity. Need non-stochastic for stable-age distribution prior to running with forage fish prey

-Code to run seabird model with base parameters and non-fished forage fish abundance (Run_seabird_2lifehistory.R)- can run this separately, or to produce results/figures, each figure code will source this code. 

-Code to run and produce Figures 4-8 (Results figures) from Koehn et al. 2021. Including -Code to run seabird model with various forage fish harvest control rule scenarios (Fig 6) 

-Code to run code to play with functional response shapes and create Figure 2 from Koehn et al. 2021





Special Thanks to: Megsie Siple (for forage fish model code), Christine Stawitz (for guidance), and Andre Punt (for specific code for functional response)
