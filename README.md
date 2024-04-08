## Directory breakdown 

- `simulation_script.R` - This is the primary script for the simulation.  This contains all of the relevant parameters that can be changed (mortality, fecundity, patch information, etc.).  To run this code, all you need to do is uncomment the specific type of simulation you want and run the code.  It will automatically create a folder in the `results` folder with the naming convention of `{simulation type}-{date}`.  The possible simulation types are:
  - `param_search` - This searches through all of the given paramets (can be updated in the relevant section) and runs 20 simulations with each unique combination of parameters.  `param_search_final` is the same, but it will only run through the top 10 results from `param_search`, it will be required to run the relevant section in `data_analysis` beforehand.
  - `green` and `green_manage` - These two will run simulations to calculate the IUCN green status.  The only difference is that `green_manage` always has management actions set into place.
  - `dispersal` and `habitat` - These two run simulations for the corresponding management actions and will have dispersal be none, targeted, random, and translocation and have habitat expand from 0 to 100%
  -  `management` - This will run simulations with varying levels of captive breeding and predator reduction
- `simulation_functions.R` - This contains the primary function that is run every timestep of the simulation.  Nothing in here should be necessary to change.
- `Fecundity_logistic.R` - This contains the function for fecundity as well as all the code to fit it to a given mortality value.
- `Patch_risk.R` - This contains the code to create logistic regression models for each patch comparing extinction risk to patch abundance.  Should be run with the results from a `patch_persistence` model from the simulation script.
- `Data_analysis.Rmd` - This contains all of the code to create the figures as shown in the manuscript.
