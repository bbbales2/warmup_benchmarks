library(tidyverse)
library(cmdstanr)

cmdstan_path = "/home/bbales2/cmdstan-warmup/"

set_cmdstan_path(cmdstan_path)

model_names = c("kilpisjarvi", "diamonds", "accel_gp", "accel_splines", "prophet", "radon")

num_chains = 4
metric = "dense_e"

runs = list()
run_number = 1

for(i in 1:length(model_names)) {
  model_name = model_names[i]

  model = cmdstan_model(paste0(cmdstan_path, "examples/", model_name, "/", model_name, ".stan"))
  model$compile()
  
  data_env = new.env(parent = baseenv())
  source(paste0(cmdstan_path, "examples/", model_name, "/", model_name, ".dat"), local = data_env)
  
  if(model_name == "kilpisjarvi") {
    max_adaptation = 7
  } else {
    max_adaptation = 11
  }
  
  for(which_adaptation in 0:max_adaptation) {
    for(chain in 1:num_chains) {
      if(run_number > length(runs)) {
        stdout = capture.output(fit <- model$sample(data = as.list.environment(data_env),
                                                   num_chains = 1,
                                                   experimental = 1,
                                                   which_adaptation = which_adaptation,
                                                   metric = metric))

        process_fit = function(files, draws, stdout) {
          stanfit = rstan::read_stan_csv(files)
          sp = rstan::get_sampler_params(stanfit)
          fitstats = monitor(draws, warmup = 0, print = FALSE)
          
          efficiency_df = tibble(rhat = max(fitstats$Rhat),
                                 min_bulk_ess = min(fitstats$Bulk_ESS),
                                 min_tail_ess = min(fitstats$Tail_ESS),
                                 lp_bulk_ess = fitstats[["lp__", "Bulk_ESS"]],
                                 lp_tail_ess = fitstats[["lp__", "Tail_ESS"]],
                                 nleapfrogs = sapply(sp, function(e) {
                                   e[, "n_leapfrog__"]
                                 }) %>% sum,
                                 ndivergences = sapply(sp, function(e) {
                                   e[, "divergent__"]
                                 }) %>% sum)
          
          picked_data = str_match(stdout[grepl("picked:", stdout)],
                                  "picked: ([0-9]+), name: (.*)")[,-1]
          colnames(picked_data) = c("iter", "which")
          picked_df = picked_data %>%
            as_tibble() %>%
            mutate(iter = as.integer(iter))
          
          adaptation_data = str_match(stdout[grepl("adapt:", stdout)],
                                      "adapt: ([0-9]*), which: ([0-9]+), min: ([0-9.]+), median: ([0-9.]+), max: ([0-9.]+), (.*)")[,-1]
          colnames(adaptation_data) = c("iter", "whichi", "min", "median", "max", "which")
          adaptation_df = adaptation_data %>%
            as_tibble() %>%
            mutate(iter = as.integer(iter),
                   whichi = as.integer(whichi),
                   min = as.numeric(min),
                   median = as.numeric(median),
                   max = as.numeric(max))
          
          return(list(efficiency_df = efficiency_df,
                      adaptation_df = adaptation_df,
                      picked_df = picked_df))
}

        out = process_fit(fit$output_files(), fit$draws(), stdout)
        out[["model_name"]] = model_name
        out[["which_adaptation"]] = which_adaptation
        out[["chain"]] = chain
        out[["metric"]] = metric
        
        print(paste0(model_name, ", which: ", which_adaptation, ", chain: ", chain, ", metric: ", metric))

        runs[[run_number]] = out
        run_number = run_number + 1
      }
    }
  }
}
