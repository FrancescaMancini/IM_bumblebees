# This is a modified version of the combine_chains function
# for model runs where the three chains are in the same model output file
# also contains an extra parameter to specify the model type
# (integrated or not)

# This function will attempt to identify the final daisy
# for each species and then combine the chains, saving
# a summary table. This assumes all species have the SAME
# number of iterations and chains (could be adapted).

# data_dir - the directory containing the output chain files
# output_dir - where the summary files will be written
# target - The desired number of iterations for estimation,
# NOTE: this is AFTER thinning (1000 iterations daisy == 334 with thinning of 3,
# 5000 = 1667)
# end_point - Numeric (or 'max') What iteration is the end of the chain?
# this allows you to assess convergence at a certain point
# in the chain. Defaults to the maximum available. Should be
# equal to the end value of a daisy (if not it is rounded up to the nearest)
# overwrite - When false species will be skipped where results already exist
# species_list - a character vetor giveing hte species to summarise. If NULL (default), all species are run
# model_type - "integrated" or "not_integrated"
library(readr)
library(abind)
combine_chains <- function(data_dir, output_dir = data_dir,
                           target = NULL, end_point = 'max',
                           verbose = TRUE, overwrite = FALSE,
                           species_list = NULL,
                           model_type = "not_integrated"){
  
  # Combine chains output on jasmin
  library(coda)
  
  on.exit(expr = {options(warn = 0)})
  options(warn = 1)
  
  all_files <- list.files(data_dir, pattern = '.rdata$', full.names = TRUE)
  
  split <- strsplit(basename(all_files), split = 'it')
  file_info <- data.frame(file = all_files,
                          final_iteration = sapply(split , function(x) parse_number(x[length(x)])),
                          stringsAsFactors = FALSE)
  # If we are using an alternative end point then find
  # it and forget all the other files
  if(end_point != 'max'){
    if(!inherits(end_point, 'numeric')) stop('end_point should be a numeric if not "max"')
    its <- sort(unique(file_info$final_iteration))
    if(!end_point %in% its) stop(paste('end_point should be the end number of iterations of a daisy but it isnt',
                                       paste0('[', end_point, ']'),
                                       'end_point values available to you are:', paste(its, collapse = ', ')))
    file_info <- file_info[file_info$final_iteration <= end_point, ]
  }
  
  if(verbose) cat('Total iterations:', max(file_info$final_iteration), '\n')
  
  if(is.null(species_list)){
    
    species <- lapply(all_files, FUN = function(x){
      paste(head(strsplit(basename(x), '_it')[[1]], -1), collapse = '_')
    })
    
    species <- unique(unlist(species))
    
  } else {
    
    species <- species_list
    
  }
  
  if(!overwrite){
    csv_outs <- list.files(output_dir, pattern = '.csv$')
    if(any(species %in% regmatches(csv_outs, regexpr('^.+(?=_it[0-9])', csv_outs, perl = TRUE)))){
      warning('Some species already have results csv, these will be skipped. Use overwrite == TRUE if you want to overwrite these species')
      species <- species[!species %in% regmatches(csv_outs, regexpr('^.+(?=_it[0-9])', csv_outs, perl = TRUE))]
    }
  }
  
  for(sp in species){
    
    # Bracket cause problems in reg exprs so they are removed here
    species_files <- all_files[grepl(x = brackOmit(basename(all_files)),
                                     pattern = paste0('^', brackOmit(sp), '_it[0-9]'))]
    
    species_info <- file_info[file_info$file %in% species_files,]
    
    # We will always need to load in the last chains
    last_it_file <- species_info$file[species_info$final_iteration == max(species_info$final_iteration)]
    
    if(verbose) cat(paste0(sp,'\n'))
    
    load(last_it_file)
    if(model_type == "not_integrated"){
      out <- out_BWARS_model
      rm('out_BWARS_model')
    } else {
      out <- out_integrated_model
      rm('out_integrated_model')
    }
    sims.list <- out$BUGSoutput$sims.array
    rm('out')
    if(verbose) cat(paste0('it-', nrow(sims.list), ' '))
    if(!is.null(target)){
      if(nrow(sims.list) >= target){
        sims.list <- sims.list[1:target,,]
        if(verbose) cat(paste0('trimmed to it-', nrow(sims.list), ' '))
      } else {
        daisies <- sort(decreasing = TRUE, unique(species_info$final_iteration))
        d <- 2
        while(nrow(sims.list) < target){
          if(d > length(daisies)) stop('Your target number of iterations is larger than the number available. NOTE: target refers to the number of iterations AFTER thinning.')
          load(species_info$file[species_info$final_iteration == daisies[d]])
          d <- d + 1
          if(model_type == "not_integrated"){
            out <- out_BWARS_model
            rm('out_BWARS_model')
          } else {
            out <- out_integrated_model
            rm('out_integrated_model')
          }
          sims.list_annex <- out$BUGSoutput$sims.array
          sims.list <- abind(sims.list, sims.list_annex, along = 1)
          rm('sims.list_annex')
          if(verbose) cat(paste0('it-', nrow(sims.list), ' '))
          if(nrow(sims.list) > target){
            sims.list <- sims.list[1:target,,]
            if(verbose) cat(paste0('trimmed to it-', nrow(sims.list), ' '))
          }
        }
      }
    }
    if(verbose) cat('\n')
    sims.list <- list(chain_1 = matrix(sims.list[1:target,1, 1:dim(sims.list)[3]], 
                                       nrow = target, ncol =dim(sims.list)[3], 
                                       dimnames = list(NULL, attr(sims.list, "dimnames")[[3]])),
                      chain_2 = matrix(sims.list[1:target,2, 1:dim(sims.list)[3]], 
                                       nrow = target, ncol =dim(sims.list)[3], 
                                       dimnames = list(NULL, attr(sims.list, "dimnames")[[3]])),
                      chain_3 = matrix(sims.list[1:target,3, 1:dim(sims.list)[3]], 
                                       nrow = target, ncol =dim(sims.list)[3], 
                                       dimnames = list(NULL, attr(sims.list, "dimnames")[[3]])))
    sims.list <- lapply(sims.list, as.mcmc)
    
    ####### it stops working at this point because the array object cannot be turned into an mcmc object ######
    
    # convert these to a mcmc.list
    comb.samples <- mcmc.list(sims.list)
    
    save(comb.samples, file = file.path(output_dir,
                               paste0(sp, '_last5000it_combined.Rdata')))
    
    
    if(verbose) cat('Done\n\n')
    
  }
}


# A function to show how rHats change over a moving window of
# iterations. This is run over a set of csv output produced
# by combine changes. For example:

# myseq <- seq(5000, 13000, by = 1000)
# 
# for(i in myseq){
#   combine_chains(data_dir = '../test_moths_1000_long/output/',
#                  target = 1670,
#                  end_point = i,
#                  verbose = TRUE)
# }
# 
# csvfiles <- list.files('../test_moths_1000_long/output/',
#                        pattern = '.csv$',
#                        full.names = TRUE)

# csvfiles is a vector of paths to all the files to use

# inspect_convergence <- function(csvfiles){
#   
#   require(scales)
#   
#   master.psi <- NULL
#   
#   for(file in csvfiles){
#     csv_data <- read.csv(file = file, stringsAsFactors = FALSE)
#     csv_data <- csv_data[grep('^psi.fs', csv_data$X), ]
#     temp <- data.frame(year = csv_data$X,
#                        species = paste(head(strsplit(basename(file), '_')[[1]], -2),
#                                        collapse = '_'),
#                        end_point = as.numeric(gsub('.csv$', '', gsub('^ep', '', strsplit(basename(file), '_')[[1]][4]))),
#                        rhat = csv_data$PSRF)
#     master.psi <- rbind(master.psi, temp)
#   }
#   
#   master.psi$year <- as.numeric(gsub('\\]$', '', gsub('^psi.fs\\[', '', master.psi$year)))
#   
#   
#   master.pdet <- NULL
#   
#   for(file in csvfiles){
#     csv_data <- read.csv(file = file, stringsAsFactors = FALSE)
#     csv_data <- csv_data[grep('^pdet.alpha', csv_data$X), ]
#     temp <- data.frame(year = csv_data$X,
#                        species = paste(head(strsplit(basename(file), '_')[[1]], -2),
#                                        collapse = '_'),
#                        end_point = as.numeric(gsub('.csv$', '', gsub('^ep', '', strsplit(basename(file), '_')[[1]][4]))),
#                        rhat = csv_data$PSRF)
#     master.pdet <- rbind(master.pdet, temp)
#   }
#   
#   master.pdet$year <- as.numeric(gsub('\\]$', '', gsub('^pdet.alpha\\[', '', master.pdet$year)))
#   
#   
#   par(mfrow = c(1,2))
#   plot(x = master.psi$end_point,
#        y = master.psi$rhat,
#        pch = 4,
#        col = alpha(as.numeric(as.factor(master.psi$species)), 0.5),
#        xlab = 'Iterations',
#        ylab = 'rHat',
#        main = 'psi.fs')
#   abline(h = 1.1, col = 'red')
#   plot(x = master.pdet$end_point,
#        y = master.pdet$rhat,
#        pch = 4,
#        col = alpha(as.numeric(as.factor(master.pdet$species)), 0.5),
#        xlab = 'Iterations',
#        ylab = 'rHat',
#        main = 'pdet.alpha')
#   abline(h = 1.1, col = 'red')
#   
#   return(list(master.psi = master.psi,
#               master.pdet = master.pdet))
# }

brackOmit <- function(x){
  
  stopifnot(inherits(x, 'character'))
  
  gsub('\\(', '', gsub('\\)', '', x))
  
}
