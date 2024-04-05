library(boot)

dir_path <- "results/patch_persistence_Mar25-12-42-53"

if(!file.exists(paste0(dir_path, 'results.rds'))){
  plan(multisession(workers = detectCores() - 1))
  
  files = 
    list.files(
      pattern = 'scenario_*', 
      path = dir_path, 
      full.names = TRUE
    )
  
  myf = \(filename){
    
    df = 
      readRDS(filename) |> 
      slice_max(year) 
    
    return(df)
  }
  
  results = 
    files |>
    future_map_dfr(myf)
  
  saveRDS(results, file = paste0(dir_path, 'results.rds'))
} else{
  results = readRDS(paste0(dir_path, 'results.rds'))
}

patch.models = list()

results$year[results$year > 100] = 100

for (i in 1:4){
  patch_risk <- results |> 
    filter(patch == i, allee)
  
  patch_risk <- patch_risk |> 
    mutate(remainder = 100 - year)
  
  # Performs a logistic regression
  patch.models[[i]] <- glm(cbind(patch_risk$year, patch_risk$remainder) ~ patch_risk$initial_n, family = binomial)
}

get.prob <- function(model, n){
  coeff = coef(model)
  
  p = (inv.logit(coeff[1] + n * coeff[2]) - inv.logit(coeff[1])) / (1 - inv.logit(coeff[1])) 
  
  return(p)
}


