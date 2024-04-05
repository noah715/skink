library(tidyverse)

threshold = matrix(
  c(
    2462, 1519, 1068, 737,
    1951, 1249,  621, 345,
    1513,  290,  195,  54,
    21,    0,    0,   0,
    0,    0,    0,   0
  ), 
  nrow = 5,
  ncol = 4,
  byrow = T
)

patch_threshold = c(271, 276, 369, 424)

model = 
  \(
    scenario = NULL,
    model,
    rho_reduction, 
    muZ_reduction, 
    muO_reduction, 
    muT_reduction, 
    muM_reduction, 
    f_reduction = 0,
    captive_breeding = 0,
    check_decreasing = F,
    gamma, 
    nu_expansion, 
    rbar,
    target_size,
    years,
    habitat_expansion,
    fecundity = fecundity.fun,
    peak_multiplier = 0.5,
    initial_n = 1000,
    nu_0 = .06,
    seed = 0,
    average_eggs = 2.17,
    slow_expansion = F,
    expansion_years = 10,
    collect_timeseries = TRUE,
    environmental_stochasticity = 0,
    dispersal_mode = 'random',
    fecundity_type = "deterministic",
    allee = TRUE,
    patch = 0,
    functional_extinction,
    vulnerable = FALSE,
    manage = FALSE,
    low_conservation = TRUE,
    green_status = FALSE,
    conserve = FALSE,
    save_file = FALSE,
    file_name = NULL
  ){
    decreasing_5 = FALSE
    decreasing_10 = FALSE
    
    # mortality
    rho_0 = (1 - rho_reduction / 100) * rho_0_0
    rho_M = (1 - rho_reduction / 100) * rho_0_M
    muZ = (1 - muZ_reduction / 100) * mu[1]
    muO = (1 - muO_reduction / 100) * mu[2]
    muT = (1 - muT_reduction / 100) * mu[3]
    muM = (1 - muM_reduction / 100) * mu[4]
    
    # fecundity parameters
    leslie <- matrix(c(   0,       0,       0, average_eggs,
                          1 - muZ,       0,       0,            0,
                          0, 1 - muO,       0,            0,
                          0,       0, 1 - muT,      1 - muM),
                     nrow = 4, ncol = 4, byrow = T)
    
    leslie.eig <- eigen (leslie)
    mean.mu <- Re(weighted.mean(c(muZ, muO, muT, muM), leslie.eig$vectors[,1]))
    
    eta = 4
    n = 2500
    carrying.cap = 12000
    patch.areas = c(130, 80, 20, 140)
    c1 = carrying.cap / sum(patch.areas)
    
    phi.bar = average_eggs / 4
    
    n = 2500
    h = sum(patch.areas)
    
    c3 <- logit(mean.mu)
    c2 <- exp((logit(phi.bar) - c3) / (1 - n/carrying.cap)) / 2500
    
    c_vals <- c(c1, c2, c3)
    
    # dispersal rate
    nu = (1 + nu_expansion / 100) * nu_0
    
    # habitat size
    if (!slow_expansion){
      patch_areas = (1 + habitat_expansion / 100) * patch_areas_0
    } else{
      patch_areas = patch_areas_0
    }
    
    if (green_status){
      patch_areas = patch_areas * 12.5
      patch_areas_0 = patch_areas_0 * 12.5
    }
    
    if(model == 'stochastic') set.seed(seed)
    
    # Carrying Capacity of each Patch.
    K_0 = target_size * (patch_areas_0 / sum(patch_areas_0))
    
    if (!slow_expansion){
      K = K_0 *  (1 + habitat_expansion / 100) #describe target_size
    } else {
      K = K_0
    }
    
    low_pop = 1
    low_pop_v = 1
    low_pop_vector = numeric(12)
    low_pop_vector_v = numeric(24)
    
    criteria = "alive"
    Status = "Endangered"
    
    
    if (patch != 0){
      nu = 0
      patch_areas[-patch] = 0
      patch_areas_0[-patch] = 0
    }
    
    # Initial abundances
    Z = matrix(NA, nrow = years, ncol = number_of_patches) # eggs
    O = matrix(NA, nrow = years, ncol = number_of_patches) # 1-year olds
    T = matrix(NA, nrow = years, ncol = number_of_patches) # 2-year olds
    M = matrix(NA, nrow = years, ncol = number_of_patches) # mature lizards
    
    Z[1, ] = c(0, 0, 0, 0)
    O[1, ] = c(0, 0, 0, 0)
    T[1, ] = c(0, 0, 0, 0)
    M[1, ] = round(initial_n * patch_areas / sum(patch_areas))
    
    year = 1
    extinction = FALSE
    
    while(year < years & extinction == FALSE){
      
      # Green status projecting into the future
      if (manage){ # Simulation starts in 1900, assumes conservation starts in 2024
        if (low_conservation){
          rho_0 = (1 - 10 / 100) * rho_0_0
          rho_M = (1 - 10 / 100) * rho_0_M
          
          captive_breeding = 40
          
          dispersal_mode = "always_success"
        } else if (!low_conservation){
          rho_0 = (1 - 20 / 100) * rho_0_0
          rho_M = (1 - 20 / 100) * rho_0_M
          
          captive_breeding = 80
          
          dispersal_mode = "always_success"
          
          if (year > 123 + 36){ # Increasing habitat area after doing some management for 3 generations (3*12)
            patch_areas = (1 + 50 / 100) * patch_areas_0
          }
        }
      }
      
      # Slowly expand habitat instead of instantly
      if (year <= expansion_years){
        patch_areas = (1 + year * (habitat_expansion / expansion_years) / 100) * patch_areas_0
        K = K_0 *  (1 + year * (habitat_expansion / expansion_years) / 100)
      }
      
      year = year + 1
      z = Z[year - 1, ]
      o = O[year - 1, ]
      t = T[year - 1, ]
      m = M[year - 1, ]
      
      if(model == 'stochastic'){
        EsZ = 1 - muZ - rho_0 + 1e-16
        EsO = 1 - muO - rho_M + 1e-16
        EsT = 1 - muT - rho_M + 1e-16
        EsM = 1 - muM - rho_M + 1e-16
        
        sZ = rtruncnorm(1, mean = EsZ, sd = environmental_stochasticity * EsZ, a = 0, b = 1)
        sO = rtruncnorm(1, mean = EsO, sd = environmental_stochasticity * EsO, a = 0, b = 1)
        sT = rtruncnorm(1, mean = EsT, sd = environmental_stochasticity * EsT, a = 0, b = 1)
        sM = rtruncnorm(1, mean = EsM, sd = environmental_stochasticity * EsM, a = 0, b = 1)
        nu_draw = rtruncnorm(1, mean = nu, sd = environmental_stochasticity * nu, a = 0, b = 1)
        
        if(dispersal_mode == "none"){
          nu = 0
          
          dist_matrix = exp(-distances / rbar) / rbar
          diag(dist_matrix) = 0
          dispersal_term = 
            rpois(
              n = number_of_patches, 
              lambda = nu / (2 * pi) / rbar * sM * colSums(dist_matrix * outer(m, sqrt(patch_areas)))
            ) 
        }
        if(dispersal_mode == 'random'){
          dist_matrix = exp(-distances / rbar) / rbar
          diag(dist_matrix) = 0
          dispersal_term = 
            rpois(
              n = number_of_patches, 
              lambda = nu / (2 * pi) / rbar * sM * colSums(dist_matrix * outer(m, sqrt(patch_areas)))
            ) 
        } 
        if(dispersal_mode == 'targeted'){
          dispersal_term = 
            rpois(
              n = number_of_patches, 
              lambda = nu / sum(patch_areas) / rbar * sM * 
                colSums(exp(-distances / rbar) * outer(m, patch_areas))
            ) 
        }
        if(dispersal_mode == 'always_success'){
          dispersal_term =
            rpois(
              n = number_of_patches,
              lambda = nu / sum(patch_areas) / rbar * sM *
                colSums(exp(-distances / rbar) * outer(m, patch_areas))
            )
        }
        
        
        
        if (green_status){
          if (!conserve){ # Will continue declining habitat
            patch_areas = patch_areas_0 * ((1/12.5)^(1/90)) ^ year
            K = K_0 * ((1/12.5)^(1/90)) ^ year
            
          } else { # Stops habitat decline after 90 years
            patch_areas = patch_areas_0 * pmax(((1/12.5)^(1/90)) ^ year, 1/12.5)
            K = K_0 * pmax(((1/12.5)^(1/90)) ^ year, 1/12.5)
          }
        }
        
        
        O[year, ] = rbinom(n = number_of_patches, prob = sZ, size = z + 
                             round(captive_breeding * patch_areas / sum(patch_areas)))
        T[year, ] = rbinom(n = number_of_patches, prob = sO, size = o)
        M[year, ] = 
          rbinom(n = number_of_patches, prob = sM * (1 - nu), size = m) + # survival of stationary adults
          rbinom(n = number_of_patches, prob = sT, size = t) + # maturation of two-yo
          dispersal_term # immigration 
        
        if (allee){
          fecundity_prob = pmax(fecundity((Z + O + T + M)[year - 1, ], patch_areas, c_vals[1], c_vals[2], c_vals[3]), 0)
          fecundity_prob = pmax(pmin((1 - f_reduction / 100) * fecundity_prob, 1), 0)
          Z[year, ] = sapply(1:4, function(x) sum(rbinom(n = round(EsM * m)[x], size = gamma, prob = fecundity_prob[x])))
          
        } else if (!allee){
          fecundity_prob = pmax(pmin(((mean.mu - 1)/ K)*(Z + O + T + M)[year - 1, ] + 1, average_eggs/4), 0)
          fecundity_prob = pmax(pmin((1 - f_reduction / 100) * fecundity_prob, 1), 0)
          
          Z[year, ] = sapply(1:4, function(x) sum(rbinom(n = round(EsM * m)[x], size = gamma, prob = fecundity_prob[x])))
        }
        
        if (dispersal_mode == "always_success"){
          mean_crowding = sum(M[year,])/(sum(K))
          K.proportion = (M[year,])/K # Gives crowding of mature individuals
          
          most_crowded = which(K.proportion == max(K.proportion))[1]
          least_crowded = which(K.proportion == min(K.proportion))[1]
          
          most_dif = min(abs(M[year, most_crowded] - mean_crowding * K[most_crowded]), 30) # How many individuals to move
          
          M[year, most_crowded] = max(0, M[year, most_crowded] - round(most_dif))
          M[year, least_crowded] = max(0, M[year, least_crowded] + round(most_dif))
        }
        
      }  
      
      if (sum((Z + O + T + M)[year, ]) == 0) {
        extinction = TRUE
        criteria = "extinct"
      }
      
      if (functional_extinction){ # Checks if critically endangered
        
        
        if (year > 36){ # Criteria A
          if (sum((Z + O + T + M)[year, ]) <= 0.1 * sum((Z + O + T + M)[year - 36, ])){
            extinction = TRUE
            criteria = "A"
            Status = "Critically endangered"
          }
        }
        
        if (sum(M[year,]) < 250){ # Criteria C
          low_pop_vector[((low_pop - 1) %% 12) + 1] = sum((Z + O + T + M)[year, ])
          
          if (low_pop >= 12){
            for (i in 1:12){
              if (low_pop_vector[(i + 10) %% 12 + 1] / low_pop_vector[i] < 0.75){
                extinction = TRUE
                criteria = "C"
                Status = "Critically endangered"
              }
            }
          }
          
          low_pop = low_pop + 1
        } else{
          low_pop = 1
        }
        
        if (sum(M[year,]) < 50){ # Criteria D
          extinction = TRUE
          criteria = "D"
          Status = "Critically endangered"
        }
        
        ext.risk = 1
        
        for (i in 1:4){
          pop.size = (Z + O + T + M)[year, i]
          
          ext.risk = ext.risk * (1 - get.prob(patch.models[[i]], pop.size))
          
        }
        
        if (sum((Z + O + T + M)[year, ]) <= threshold[captive_breeding %/% 30 + 1, rho_reduction %/% 20 + 1]){ # Criteria E
          # if (ext.risk > 0.5){
          extinction = TRUE
          criteria = "E"
          Status = "Critically endangered"
        }
        
        if (extinction){
          break 
        }
      } 
      
      if (vulnerable){ # Checks if vulnerable
        
        # All 3 criteria have to be true to pass
        A = FALSE
        C = FALSE
        D = FALSE
        
        
        if (year > 36){ # Criteria A
          if (!(sum((Z + O + T + M)[year, ]) <= 0.3 * sum((Z + O + T + M)[year - 36, ]))){
            A = TRUE
          } else {
            A = FALSE
          }
        }
        
        
        if (!(sum(M[year,]) < 2500)){ # Criteria C
          C = TRUE
          low_pop_v = 1
        } else {
          low_pop_vector_v[((low_pop_v - 1) %% 24) + 1] = sum((Z + O + T + M)[year, ])
          
          if (low_pop_v >= 24){
            if (!(low_pop_vector_v[(low_pop_v - 1) %% 24 + 1] / low_pop_vector_v[(low_pop_v) %% 24 + 1] < 0.8)){
              C = TRUE
            } else {
              C = FALSE
            }
          }
          low_pop_v = low_pop_v + 1
        }
        
        
        if (!(sum(M[year,]) < 250)){ # Criteria D
          D = TRUE
        } else {
          D = FALSE
        }
        
        if (A & C & D){
          extinction = TRUE
          Status = "Vulnerable"
          criteria = "Vulnerable"
        }
        
      }
      
    }
    
    extinction_year = year
    
    green = 0
    
    for (i in 1:4){
      pop.size = (Z + O + T + M)[year, i]
      
      green = green + get.prob(patch.models[[i]], pop.size)
    }
    
    if (check_decreasing & year >= 128){
      if (sum((Z + O + T + M)[128, ]) < sum((Z + O + T + M)[123, ])){
        decreasing_5 = TRUE
      }
    }
    if (check_decreasing & year >= 133){
      if (sum((Z + O + T + M)[133, ]) < sum((Z + O + T + M)[123, ])){
        decreasing_10 = TRUE
      }
    }
    
    if(collect_timeseries == TRUE){
      dtf_outcome =
        rbind(Z, O, T, M) |>
        data.frame() |>
        drop_na() |>
        mutate(
          year = rep(seq(extinction_year), 4),
          stage = rep(c('egg', 'one', 'two', 'mature'), each = extinction_year)
        ) |>
        as_tibble() |>
        rename_with(
          .fn = \(.) return(patch_names),
          .cols = X1:X4
        )
    }
    
    if(collect_timeseries == FALSE){
      dtf_outcome = 
        tibble(
          year = year,
          egg = Z[year, ],
          one = O[year, ],
          two = T[year, ],
          mature = M[year, ],
          patch = patch_names
        ) |>
        pivot_longer(
          egg:mature, 
          names_to = 'stage', 
          values_to = 'abundance'
        )
    }
    
    result = 
      tibble(
        model,
        habitat_expansion,
        rho_reduction, 
        muZ_reduction, 
        muO_reduction, 
        muT_reduction, 
        muM_reduction, 
        average_eggs,
        decreasing_5,
        decreasing_10,
        captive_breeding,
        criteria,
        gamma, 
        nu_expansion, 
        rbar,
        fecundity_type,
        allee,
        target_size,
        years,
        dispersal_mode,
        initial_n,
        nu_0,
        environmental_stochasticity,
        seed,
        patch,
        green, 
        Status,
        functional_extinction,
        conserve,
        manage,
        low_conservation,
        vulnerable,
        dtf_outcome,
        expansion_years,
        slow_expansion,
        peak_multiplier,
        scenario
      )
    
    if(save_file){
      saveRDS(result, file = file_name)
    }
    if(!save_file){
      return(result)
    }
  }

model_summary = 
  \(stages){
    dtf = 
      stages |> 
      pivot_longer(Hyde:Pukerangi, names_to = 'patch', values_to = 'abundance') |>
      group_by(across(model:patch)) |>
      summarize(abundance = 2 * abundance, .groups = 'drop')
  }

