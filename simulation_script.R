library(tidyverse)
library(parallel)
library(furrr)
library(truncnorm)
library(scales)
library(matrixcalc)
library(cowplot)
library(lemon)
library(stats)
library(boot)



time <- format(Sys.time(), "%b%d-%H-%M-%S")
prefix <- "dispersal_"
# prefix <- "habitat_"
# prefix <- "management_"
# prefix <- "green_"
# prefix <- "green_conserve_"
# prefix <- "green_future_"
# prefix <- "param_search_"
# prefix <- "param_search_final_"
# prefix <- "patch_persistence_"
path <- paste0(prefix, time)
    
source("simulation_functions.R")
source("Fecundity_logistic.R")
source("Patch_risk.R")
setwd("./results")

dir.create(path)
setwd(path)

plan(multisession(workers = detectCores() - 1))  



# Patch names
patch_names = c('Hyde', 'Macraes', 'Sutton', 'Pukerangi')

# Surface Area of the 4 Patches.
patch_areas_0 = c(130, 80, 20, 14)

# Number of islands/patches
number_of_patches = length(patch_areas_0)

# Total Current Population. The source said its between 2000 to 5000. Here we only count females.
initial_n = 1000

# Distance Between each Patch.
distances = 
   matrix(
      c(
         0  , 16 , 420 ,  38 , 
         16 , 0  , 423 ,  33 ,  
         420, 423, 0   ,  455,
         38 , 33 , 455 ,   0
      ),
      nrow = number_of_patches,
      ncol = number_of_patches
   )

# Mortality rates
mu = c(.65, .6, .35, .1)

# Predation-induced mortality
rho_0_0 = .069
rho_0_M = .069
# rho_0_0 = 0
# rho_0_M = 0

# Default parameters for each model
seeds = 1:50
n.years = 101

muZ.red = 20
muO.red = 20
muT.red = 10
muM.red = 20
avg.eggs = 2.10

rho.red = 0
captive.breeding = 0
dispersal.mode = 'targeted'
habitat.expansion = 0
slow.expansion = F
expansion.years = 0

target.size = 12000
init.n = 2500
patches = 0

functional.extinction = TRUE
vuln = TRUE

green.status = FALSE
green.conserve = FALSE
green.manage = FALSE
low.conservation = TRUE

check.decreasing = FALSE

start <- proc.time()
if(!seawulf){
   
  #### Management ####
  if (prefix == "management_"){ 
    seeds = 1:250
    
    rho.red = seq(0, 60, by = 5)
    captive.breeding = seq(0, 120, by = 10)
  } 
  #### Habitat expansion ####
  else if (prefix == "habitat_"){
    seeds = 1:250
    
    vuln = FALSE
    
    habitat.expansion = seq(0, 100, by = 10)
    
  } 
  #### Dispersal ####
  else if (prefix == "dispersal_"){
    seeds = 1:250
    
    vuln = FALSE
    
    dispersal.mode = c("targeted", "random", "always_success", "none")
  }
  #### Green status ####
  else if (prefix == "green_"){
    seeds = 1:250
    n.years = c(123, 123 + 36, 223)
    
    init.n = 150000
    
    functional.extinction = FALSE
    vuln = FALSE
    
    green.status = TRUE
    green.conserve = c(TRUE, FALSE)
    green.manage = c(TRUE, FALSE)
    low.conservation = c(TRUE, FALSE)
    
  } 
  #### Green status conserve ####
  else if (prefix == "green_conserve_"){
    seeds = 1:250
    n.years = c(36, 100)
    
    init.n = 2500
    
    functional.extinction = FALSE
    vuln = FALSE
    
    green.status = TRUE
    green.conserve = c(TRUE)
    green.manage = c(TRUE)
    low.conservation = c(TRUE, FALSE)
    
  } 
  #### Parameter search ####
  else if (prefix == "param_search_"){
    seeds = 1:20
    n.years = 133
    
    init.n = 150000
    
    functional_extinction = FALSE
    vuln = FALSE
    
    muZ.red = seq(-30, 30, length = 7)
    muO.red = seq(-30, 30, length = 7)
    muT.red = seq(-30, 30, length = 7)
    muM.red = seq(-30, 30, length = 7)
    avg.eggs = seq(1.03, 2.82, length = 11)
    
    green.status = TRUE
    green.conserve = TRUE
    
    check.decreasing = TRUE
  }  
  else if (prefix == "param_search_final_"){
    for (i in 1:10){
      seeds = 1:500
      n.years = 133
      
      init.n = 150000
      
      functional_extinction = FALSE
      vuln = FALSE
      
      muZ.red = param.dif[[i, 1]]
      muO.red = param.dif[[i, 2]]
      muT.red = param.dif[[i, 3]]
      muM.red = param.dif[[i, 4]]
      avg.eggs = seq(1.03, 2.82, length = 11)
      
      green.status = TRUE
      green.conserve = TRUE
      
      check.decreasing = TRUE
      
      params =
        expand_grid(
          # Simulation options
          model = 'stochastic',
          seed = seeds,
          years = n.years,
          
          # Life history parameter changes
          muZ_reduction = muZ.red, 
          muO_reduction = muO.red, 
          muT_reduction = muT.red, 
          muM_reduction = muM.red, 
          average_eggs = avg.eggs,
          
          # Management options
          rho_reduction = rho.red, 
          captive_breeding = captive.breeding,
          dispersal_mode = dispersal.mode,
          habitat_expansion = habitat.expansion, # Additional options to customize habitat expansion
          slow_expansion = slow.expansion,
          expansion_years = expansion.years,
          
          # Carrying capacity and initial population
          target_size = target.size,
          initial_n = init.n,
          patch = patches,
          
          # IUCN criteria
          functional_extinction = functional.extinction, # T/F stop when pop reaches critical endangerment
          vulnerable = vuln, # T/F Stop when pop reaches vuln
          
          # Green status
          green_status = green.status, # T/F should the simulation start at a historic time period
          conserve = green.conserve, # T/F should previous conservation efforts be added (conservation legacy)
          manage = green.manage, # T/F should future management actions be added
          low_conservation = low.conservation, # T/F should there be low or high conservation in the future
          
          # Parameter fit
          check_decreasing = check.decreasing,
          
          # Constants across all simulations
          gamma = 4, 
          nu_0 = .06,
          nu_expansion = 0, 
          rbar = 25,
          collect_timeseries = TRUE,
          environmental_stochasticity = 0,
          allee = c(TRUE, FALSE),
          
        ) |>
        rowid_to_column(var = 'scenario') |>
        mutate(file_name = paste0(i, '_scenario_', scenario, '.rds'))
      
      result =
        params |>
        future_pmap_dfr(
          model,
          save_file = TRUE,
          .options = furrr_options(seed = NULL)
        )
      
    }
    break
  } 
  #### Risk assessment ####
  else if (prefix == "patch_persistence_"){
    dispersal.mode = "none"
    
    seeds = 1:100
    
    init.n = seq(0, 2500, by = 10)
    patches = 1:4
    
    functional.extinction = FALSE
    vuln = FALSE
  }
  
  params =
    expand_grid(
      # Simulation options
      model = 'stochastic',
      seed = seeds,
      years = n.years,
      
      # Life history parameter changes
      muZ_reduction = muZ.red, 
      muO_reduction = muO.red, 
      muT_reduction = muT.red, 
      muM_reduction = muM.red, 
      average_eggs = avg.eggs,
      
      # Management options
      rho_reduction = rho.red, 
      captive_breeding = captive.breeding,
      dispersal_mode = dispersal.mode,
      habitat_expansion = habitat.expansion, # Additional options to customize habitat expansion
        slow_expansion = slow.expansion,
        expansion_years = expansion.years,
      
      # Carrying capacity and initial population
      target_size = target.size,
      initial_n = init.n,
      patch = patches,
      
      # IUCN criteria
      functional_extinction = functional.extinction, # T/F stop when pop reaches critical endangerment
      vulnerable = vuln, # T/F Stop when pop reaches vuln
      
      # Green status
      green_status = green.status, # T/F should the simulation start at a historic time period
      conserve = green.conserve, # T/F should previous conservation efforts be added (conservation legacy)
      manage = green.manage, # T/F should future management actions be added
      low_conservation = low.conservation, # T/F should there be low or high conservation in the future
      
      # Parameter fit
      check_decreasing = check.decreasing,
      
      # Constants across all simulations
      gamma = 4, 
      nu_0 = .06,
      nu_expansion = 0, 
      rbar = 25,
      collect_timeseries = TRUE,
      environmental_stochasticity = 0,
      allee = c(TRUE, FALSE),
      
    ) |>
    rowid_to_column(var = 'scenario') |>
    mutate(file_name = paste0('scenario_', scenario, '.rds'))

    result =
      params |>
      future_pmap_dfr(
         model,
         save_file = TRUE,
         .options = furrr_options(seed = NULL)
      )
   
}
