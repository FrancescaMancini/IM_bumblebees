# This function will attempt to identify the final daisy
# for each species and then combine the chains, saving
# a summary table. This assumes all species have the SAME
# number of iterations and chains (could be adapted).

# data_dir - the directory containing the output chain files
# output_dir - where the summary files will be written
# dataset - the name of the dataset used, if we're building multiple models from different datasets
# target - The desired number of iterations.
# NOTE: this is AFTER thinning (1000 iterations daisy == 334 with thinning of 3,
# 5000 = 1667)
# Currently we just
# take the last <target> values but it could be changed to take this number
# of values from the last half of the total chain length
# verbose - write stuff to console? T/F
# sample_size - how many samples to take from the posterior, defaults to 1000
# overwrite - When false species will be skipped where results already exist
# column_regex - the regular expression to select the columns from sims.matrix
# that you want to keep, default to all
# species_list = optional list of species to run and only run these species, as a character vector

extract_posterior <- function(data_dir, output_dir = data_dir, dataset = NULL, min_year = 2010,
                              max_year = 2016, target = NULL, verbose = TRUE, sample_size = 1000,
                              overwrite = FALSE, column_regex = NULL, species_list = NULL, model_type = "not_integrated"){
  
  require(dplyr)
  require(readr)
  require('utils')
  if(is.null(target)){
    warning('You have not set a target.  This means that only the last iteration of the daisy chain will be sampled, which is probably not what you want')
  }
  
  # Check and create output directory
  if(!dir.exists(output_dir)){
    dir.create(output_dir)
  }
  
  # Combine chains output on jasmin
  # library(coda)
  brackOmit <- function(x){
    
    stopifnot(inherits(x, 'character'))
    
    gsub('\\(', '', gsub('\\)', '', x))
    
  }
  
  on.exit(expr = {options(warn = 0)})
  options(warn = 1)
  
  source_dir <- data_dir
  if(!is.null(dataset)){source_dir <- file.path(data_dir, dataset)}
  
  all_files <- list.files(source_dir, pattern = '.rdata$', full.names = TRUE, recursive = TRUE)
  
  split <- strsplit(basename(all_files), split = 'it')
  split_length <- unlist(lapply(split, function(x) length(x)))
  if(length(unique(split_length))>1){
    stop('Check your data files.  Some of the data files are not output files in the data directory passed to function')
  }
  file_info <- data.frame(file = all_files,
                          final_iteration = as.numeric(gsub('.rdata', '', sapply(split , function(x) tail(x,1)))),
                          stringsAsFactors = FALSE)
  
  if(verbose) cat('Total iterations:', max(file_info$final_iteration), '\n')
  
  species <- lapply(all_files, FUN = function(x){
    paste(head(strsplit(basename(x), '_it')[[1]], -1), collapse = '_')
  })
  file_info$species <- species %>% unlist()
  species <- unique(file_info$species)
  
  if(!is.null(species_list)){
    if(!all(species_list %in% species)){
      cat('Some species listed are not found in the species list.\n\n Species list =',species,'\n\n Your list =',species_list,'\n\n')
      stop()
    }
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
    species_info <- file_info[grepl(x = brackOmit(basename(file_info$file)),
                                    pattern = paste0('^', brackOmit(sp), '_it[0-9]')),]
    
    # We will always need to load in the last chains
    last_chain <- species_info$file[species_info$final_iteration == max(species_info$final_iteration)]
    
    if(verbose) cat(paste0(ifelse(!is.null(dataset),paste0(dataset,', '),''),sp,'\n'))
    #
    #    sims.list <- lapply(unique(species_info$chain), FUN = function(chain){
    #
    #      if(verbose) cat(paste0(' Chain-', chain, ': '))
    #
    #      chainfiles <- species_info[species_info$chain == chain, ]
    #
    # Load the last chain
    load(last_chain)
    if(model_type == "not_integrated"){
      out <- out_BWARS_model
      rm('out_BWARS_model')
    } else {
      out <- out_integrated_model
      rm('out_integrated_model')
    }
    sims.matrix <- out$BUGSoutput$sims.matrix
    if(!is.null(column_regex)){
      sims.matrix <- sims.matrix[,grepl(column_regex, colnames(sims.matrix))]
    }
    rm('out')
    if(verbose) cat(paste0('it-', nrow(sims.matrix), ' '))
    
    if(!is.null(target)){
      if(nrow(sims.matrix) >= target){
        sims.matrix <- tail(sims.matrix, target)
        if(verbose) cat(paste0('\n trimmed to it-', nrow(sims.matrix), ' '))
      } else {
        daisies <- sort(decreasing = TRUE, unique(species_info$final_iteration))
        d <- 2
        while(nrow(sims.matrix) < target){
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
          sims.matrix_annex <- out$BUGSoutput$sims.matrix
          if(!is.null(column_regex)){
            sims.matrix_annex <- sims.matrix_annex[,grepl(column_regex, colnames(sims.matrix_annex))]
          }
          sims.matrix <- rbind(sims.matrix_annex, sims.matrix)
          rm('sims.matrix_annex')
          if(verbose) cat(paste0('it-', nrow(sims.matrix), ' '))
          if(nrow(sims.matrix) > target){
            sims.matrix <- tail(sims.matrix, target)
            if(verbose) cat(paste0('\n trimmed to it-', nrow(sims.matrix), ' '))
          }
        }
      }
    }
    if(ncol(sims.matrix) == length(min_year:max_year)){
      colnames(sims.matrix) <- paste('year', min_year:max_year, sep = '_')
    }
    
    # rbind and sample the posteriors
    if(sample_size <= nrow(sims.matrix)){
      samples <- sims.matrix[sample(1:nrow(sims.matrix), size = sample_size, replace = FALSE), ]
    } else {
      warning(paste0('Your sample size (',sample_size,
                     ') is greater than the number of rows of posterior data (',nrow(sims.matrix),
                     ').  You should probably try a larger value for target, ',
                     'or smaller value for sample size'))
      samples <- sims.matrix
    }
    samples <- cbind(samples, 1:nrow(samples), sp)
    colnames(samples)[ncol(samples)-1] <- 'iteration'
    colnames(samples)[ncol(samples)] <- 'species'
    if(!is.null(dataset)){
      samples <- cbind(samples, 1:nrow(samples), dataset)
      colnames(samples)[ncol(samples)] <- 'dataset'
      filename <- paste0(paste(dataset, sp, sample_size, 'posterior', sep = '_'), '.csv')
    } else {
      filename <- paste0(paste(sp, sample_size, 'posterior', sep = '_'), '.csv')
    }
    write.csv(samples,
              file = file.path(output_dir, filename),
              row.names = FALSE)
  }
  return('Another job well done')
}