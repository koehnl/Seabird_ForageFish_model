# Seabird_ForageFish_model
R Code to run models from Koehn et al. 2021 - A structured seabird population model reveals how alternative forage fish control rules benefit seabirds and fisheries. Ecological Applications

Includes:

-Seabird model function

-Forage fish model function and associated functions from or modified from Siple et al. 2019 
Siple, M.C., Essington, T.E. and E. Plagányi, É., 2019. Forage fish fisheries management requires a tailored approach to balance trade‐offs. Fish and Fisheries, 20(1), pp.110-124.

-Code to run seabird model with base parameters and non-fished forage fish abundance

-Code to run seabird model with various forage fish harvest control rule scenarios

Work flow
1. Run Forage fish model (Run_ForageFishmodel_Siple_etal_2019.R) with first one forage fish (anchovy or sardine) and then the other specified. 
2. This will save runs for both anchovy and sardine
3. Load runs for one forage fish 
4. Run Run_seabird_3types.R with that one forage fish
5. Then run scenario runs at the top of the script for any of the figures and save those specific runs (code for saving)
6. Repeat steps 3-5 with the other forage fish runs - save runs but still also have the runs saved for the previous forage fish and can plot both together
(Don't need to repeat if only one forage fish you want to run for but then figures won't work because rely on running for both different forage fish prey)


NOTE - FIGURE 5 code requires running with different forage fish and jumping between both R code files depending on if you just want functional response sensitivity,
life history sensitivity, or both - VERY MESSY
