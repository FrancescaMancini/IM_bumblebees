####################################
## Diagnostics
## Author: Francesca Mancini
##         Nick Isaac
## Date created: 2020-06-18
## Date modified: 
####################################


library(sparta)
library(rjags)
library(R2jags)
library(MCMCvis)
library(reshape2)
library(ggplot2)
library(GGally)
source("./Functions/combine_chains_no_summary.R")


common_species <- read.csv("/data/input-data/common_species.csv")

# combine the last 5000 iterations in one mcmc list

combine_chains(data_dir = "/data/outputs/Occ_model/", 
               output_dir = "/data/outputs/Occ_model/Last5000combined/",
               target = 1665, end_point = 'max',
               verbose = TRUE, overwrite = FALSE,
               species_list = common_species$species,
               model_type = "not_integrated")


combine_chains(data_dir = "/data/outputs/Integrated/", 
               output_dir = "/data/outputs/Integrated/Last5000combined/",
               target = 1665, end_point = 'max',
               verbose = TRUE, overwrite = FALSE,
               species_list = common_species$species,
               model_type = "integrated")

# modify MCMC_plots function so it reads the combined iterations files
# plots and saves as pdf the trace plots for each species
# I want to get the pairs plots of parameters too
# might need to loop through years to incorporate Nick's code


trace_posterior_plots <- function(species, common_parameters = c("a", "psi.fs"), 
                                  all_parameters = c("beewalks.p"),
                                  jags_occmodel_output_path = getwd(), 
                                  jags_integrated_output_path = getwd(),
                                  output_path = getwd()){
  
  load(file.path(jags_occmodel_output_path, paste0(species, "_last5000it_combined.Rdata")))
  comb.samples.occmodel <- comb.samples
  rm(comb.samples)
  
  
  load(file.path(jags_integrated_output_path, paste0(species, "_last5000it_combined.Rdata")))
  comb.samples.integrated <- comb.samples
  rm(comb.samples)
  
  for(i in 1:length(all_parameters)){
    if(all_parameters[i] != "beewalks.p"){
      
      MCMCtrace(comb.samples.occmodel, 
                params = all_parameters[i], 
                ind = TRUE, pdf = TRUE,
                filename = paste0(output_path, paste(species, "occmodel", 
                                                     all_parameters[i], "MCMCtrace.pdf", 
                                                     sep = "_")))
      
      MCMCtrace(comb.samples.integrated, 
                params = all_parameters[i], 
                ind = TRUE, pdf = TRUE,
                filename = paste0(output_path, paste(species, "integrated", 
                                                     all_parameters[i], "MCMCtrace.pdf", 
                                                     sep = "_")))
    }else{
      MCMCtrace(comb.samples.integrated, 
                params = all_parameters[i], 
                ind = TRUE, pdf = TRUE,
                filename = paste0(output_path, paste(species, "integrated", 
                                                     all_parameters[i], "MCMCtrace.pdf", 
                                                     sep = "_")))
    }
  }
  
  for(j in 1:length(common_parameters)){
    
    
    pdf(paste0(output_path, paste(species, common_parameters[j], "MCMCplot.pdf", sep = "_")))
    
    
    MCMCplot(object = comb.samples.occmodel, 
             object2 = comb.samples.integrated,
             params = common_parameters[j], 
             horiz = FALSE,
             rank = FALSE,
             ref_ovl = FALSE, 
             xlab = 'My x-axis label', 
             main = paste(species, common_parameters[j], sep = " - "))
    dev.off()
  }
}

lapply(common_species$species, trace_posterior_plots, 
       common_parameters = c("a", "beta1", "beta2", "beta3", 
                             "dtype1.p", "dtype2.p", "dtype3.p", "psi.fs"),
       all_parameters = c("a", "beewalks.p", "beta1", "beta2", "beta3", 
                          "dtype1.p", "dtype2.p", "dtype3.p", "psi.fs"),
       jags_occmodel_output_path = "/data/outputs/Occ_model/Last5000combined/", 
       jags_integrated_output_path = "/data/outputs/Occ_model/Last5000combined/",
       output_path = "/data/outputs/diagnostics/")


## parameter correlations ----


plotChains <- function(species, jags_occmodel_output_path = getwd(),
                       jags_integrated_output_path = getwd(),
                       output_path = getwd()){
  # function that takes a sims.array and makes a traceplot and optionally a pairs plot
  # optional argument defines which year's data to look at (defaults to first year)
  # NB this currently works only for multi-year datasets. Ideally we'd also allow a range of years
  
  load(file.path(jags_occmodel_output_path, paste0(species, "_last5000it_combined.Rdata")))
  sims.array.occmodel <- as.array(comb.samples)
  rm(comb.samples)
  
  load(file.path(jags_integrated_output_path, paste0(species, "_last5000it_combined.Rdata")))
  sims.array.integrated <- as.array(comb.samples)
  rm(comb.samples)
  
  # get the data in a format that we want it
  sims.occmodel <- melt(sims.array.occmodel, varnames = c("iteration", "param", "chain"))
  sims.integrated <- melt(sims.array.integrated, varnames = c("iteration", "param", "chain"))
  
  temp.occmodel <- dcast(sims.occmodel, iteration+chain~param)
  temp.integrated <- dcast(sims.integrated, iteration+chain~param)
  
  names(temp.occmodel)[grepl(paste0("\\[","[0-9]","\\]$"), 
                             x=names(temp.occmodel))] <- gsub("\\]", "", 
                                                              gsub("\\[", "_", 
                                                                   names(temp.occmodel[,grepl(paste0("\\[","[0-9]","\\]$"), 
                                                                                              x=names(temp.occmodel))])))
  
  names(temp.integrated)[grepl(paste0("\\[","[0-9]","\\]$"), 
                               x=names(temp.integrated))] <- gsub("\\]", "", 
                                                                  gsub("\\[", "_", 
                                                                       names(temp.integrated[,grepl(paste0("\\[","[0-9]","\\]$"), 
                                                                                                    x=names(temp.integrated))])))
  
  temp.occmodel <- pivot_longer(temp.occmodel, cols = contains("_"),
                                names_to = c(".value", "year"), 
                                names_sep = "_")
  
  temp.integrated <- pivot_longer(temp.integrated, cols = contains("_"),
                                  names_to = c(".value", "year"), 
                                  names_sep = "_")
  
  
  #pp1 <- pairs(temp[, -c(1:2)], col=temp$chain)
  pp.occmodel <- ggpairs(temp.occmodel[, -c(1,2,9)], ggplot2::aes(col=factor(temp.occmodel$chain)))
  pp.integrated <- ggpairs(temp.integrated[, -c(1,2,9)], ggplot2::aes(col=factor(temp.integrated$chain)))
  
  ggsave(pp.occmodel, 
         filename = file.path(output_path, paste0(species, "_occmodel_pplot.png")), 
         width = 18, height = 7, units = "in")
  
  ggsave(pp.integrated, 
         filename = file.path(output_path, paste0(species, "_integrated_pplot.png")), 
         width = 18, height = 7, units = "in")
  
  
} 


####################### end of functions

# apply to all species

lapply(common_species$species, plotChains, 
       jags_occmodel_output_path = "/data/outputs/Occ_model/Last5000combined/",
       jags_integrated_output_path = "/data/outputs/Integrated/Last5000combined/",
       output_path = "/data/outputs/diagnostics/")


## detection phenology ----

plot_detection_phenology <- function(species = "bombus lapidarius",
                                     jags_occmodel_output_path = "./Outputs/Occ_model/200000_it/",
                                     jags_integrated_output_path = "./Outputs/Integrated/200000_it/",
                                     output_path = "./Outputs/Diagnostics/"){
  
  require(patchwork)
  
  source("./Functions/DetectionPhenology.R")
  
  load(paste0(jags_occmodel_output_path, species, "_it2e+05.rdata"))
  
  occ_mod_det_phenology <- plot_DetectionPhenology_modified(out_BWARS_model, density_function = TRUE) +
    ggtitle(paste(species, "Occupancy model", sep = " - "))
  
  
  load(paste0(jags_integrated_output_path, species, "_it2e+05.rdata"))
  
  int_mod_det_phenology <- plot_DetectionPhenology_modified(out_integrated_model, density_function = TRUE)+
    ggtitle(paste(species, "Integrated model", sep = " - "))
  
  det_phen_plot <- occ_mod_det_phenology + int_mod_det_phenology
  
  ggsave(det_phen_plot, filename=paste0(output_path, 
                                        paste(species, "detection_phenology.png", sep = "_")),
         width = 10, height = 5, units = "in")
  
  
}

lapply(common_species$species, plot_detection_phenology,
       jags_occmodel_output_path = "/data/outputs/Occ_model/Last5000combined/",
       jags_integrated_output_path = "/data/outputs/Integrated/Last5000combined/",
       output_path = "/data/outputs/diagnostics/")
