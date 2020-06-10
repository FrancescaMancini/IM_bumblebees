#########################################
## Script to create list of parameters 
## and job script for JASMIN
## Adapted from: Mark Logie
## Date created: 2019-12-17
## Date modified: 2020-03-10
#########################################

rm(list = ls())

library(dplyr)


## create the parameters file ----

# starting with bombus species
taxa <- read.csv("./Data/common_species.csv", stringsAsFactors = FALSE)[, "species"]

chain <- 100
max_it <- 100000

write.csv(data.frame(species = lapply(taxa,
                                      FUN = function(x){rep(x,chain)}) %>% unlist(),
                     start = rep(seq(1,((max_it/chain)*(chain-1)+1),(max_it/chain)),length(taxa)),
                     end = rep(seq((max_it/chain),max_it,(max_it/chain)),length(taxa)),
                     chain = rep(1:100,length(taxa)),
                     seed = round(runif(n = chain*length(taxa), min = 0, max = 2147483647),0)),
          './Data/parameters.csv',row.names = FALSE)

# chain <- 200
# start_it <- 200000
# max_it <- 400000
# 
# write.csv(data.frame(species = lapply(taxa,
#                                       FUN = function(x){rep(x,chain)}) %>% unlist(),
#                      start = rep(seq(start_it,((max_it/chain)*(chain-1)+1),(max_it/chain)),length(taxa)),
#                      end = rep(seq((max_it/chain),max_it,(max_it/chain)),length(taxa)),
#                      chain = rep(200:400,length(taxa)),
#                      seed = round(runif(n = chain*length(taxa), min = 0, max = 2147483647),0)),
#           'parameters.csv',row.names = FALSE)


job_num <- paste0(seq(1,(chain*length(taxa)-chain+1),chain),'-',
                  seq(chain,(chain*length(taxa)),chain))

## create job script for occ model ----

bsub_lines <- paste0('bsub -q short-serial -n 1 -R "rusage[mem=20000]" -M 20000',
                     ' -W 23:59 -J BWARS_modelling_Job[', job_num,
                     ']%1 -oo Outputs/Log_files/R-%J-%I.o -eo Outputs/Log_files/R-%J-%I.e < multi_array.job')
write.table(bsub_lines, file = 'bsub_lines.sh', sep = '',
            quote = FALSE, row.names = FALSE, col.names = FALSE)

## create job script for integrated model ----

bsub_lines_int <- paste0('bsub -q short-serial -n 1 -R "rusage[mem=20000]" -M 20000',
                         ' -W 23:59 -J Integrated_modelling_Job[', job_num,
                         ']%1 -oo Outputs/Log_files_int/R-%J-%I.o -eo Outputs/Log_files_int/R-%J-%I.e < multi_array_int.job')
write.table(bsub_lines_int, file = 'bsub_lines_int.sh', sep = '',
            quote = FALSE, row.names = FALSE, col.names = FALSE)


# Run create_visit_data.R
# chmod +x bsub_lines.sh (may need to open in WinSCP, then save, to remove windows linebreaks)
# ./bsub_lines.sh

# Restart failures
# From ./functions
# source('restart_failures.R')
# restart_failures(dir = '../../dragonflies',originalJob = '../../dragonflies/bsub_lines.sh')

# Combine chains
# From ./functions
# source('combine_chains.R')
# combine_chains(data_dir = '../../dragonflies', target = 3340)

# Get posteriors
# From ./functions
# source('extract_posteriors.R')
# extract_posteriors(data_dir = '../../dragonflies/output',output_dir = '../../dragonflies/output/posteriors, target = 3340)

# Plot occupancy
# source('plot_occupancy.R')
# summary_occ <- get_summary_occ(output_dir = '../../dragonflies', source_data = '../../dragonflies/data/dragonflies_outputs.rdata')
# NOT RUN create_summary_plots(summary_occ, plots_dir = '../../dragonflies/plots')

# Calculate and plot trend estimates
# Run full calculate_trend.R script

