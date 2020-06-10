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
    max.it <- nrow(comb.samples[[1]])
    
    p <- summary(comb.samples)
    stats <- as.data.frame(p$statistics[,1:4])
    quantiles <- as.data.frame(p$quantiles[,1:5])
    
    gelman.hack <- function (x, confidence = 0.95, transform = FALSE, autoburnin = TRUE, 
                             multivariate = TRUE) 
    {
      x <- as.mcmc.list(x)
      if (nchain(x) < 2) 
        stop("You need at least two chains")
      if (autoburnin && start(x) < end(x)/2) 
        x <- window(x, start = end(x)/2 + 1 ) # Hacked here
      Niter <- niter(x)
      Nchain <- nchain(x)
      Nvar <- nvar(x)
      xnames <- varnames(x)
      if (transform) 
        x <- gelman.transform(x)
      x <- lapply(x, as.matrix)
      S2 <- array(sapply(x, var, simplify = TRUE), dim = c(Nvar, 
                                                           Nvar, Nchain))
      W <- apply(S2, c(1, 2), mean)
      xbar <- matrix(sapply(x, apply, 2, mean, simplify = TRUE), 
                     nrow = Nvar, ncol = Nchain)
      B <- Niter * var(t(xbar))
      if (Nvar > 1 && multivariate) {
        if (is.R()) {
          CW <- chol(W)
          emax <- eigen(backsolve(CW, t(backsolve(CW, B, transpose = TRUE)), 
                                  transpose = TRUE), symmetric = TRUE, only.values = TRUE)$values[1]
        }
        else {
          emax <- eigen(qr.solve(W, B), symmetric = FALSE, 
                        only.values = TRUE)$values
        }
        mpsrf <- sqrt((1 - 1/Niter) + (1 + 1/Nvar) * emax/Niter)
      }
      else mpsrf <- NULL
      w <- diag(W)
      b <- diag(B)
      s2 <- matrix(apply(S2, 3, diag), nrow = Nvar, ncol = Nchain)
      muhat <- apply(xbar, 1, mean)
      var.w <- apply(s2, 1, var)/Nchain
      var.b <- (2 * b^2)/(Nchain - 1)
      cov.wb <- (Niter/Nchain) * diag(var(t(s2), t(xbar^2)) - 2 * 
                                        muhat * var(t(s2), t(xbar)))
      V <- (Niter - 1) * w/Niter + (1 + 1/Nchain) * b/Niter
      var.V <- ((Niter - 1)^2 * var.w + (1 + 1/Nchain)^2 * var.b + 
                  2 * (Niter - 1) * (1 + 1/Nchain) * cov.wb)/Niter^2
      df.V <- (2 * V^2)/var.V
      df.adj <- (df.V + 3)/(df.V + 1)
      B.df <- Nchain - 1
      W.df <- (2 * w^2)/var.w
      R2.fixed <- (Niter - 1)/Niter
      R2.random <- (1 + 1/Nchain) * (1/Niter) * (b/w)
      R2.estimate <- R2.fixed + R2.random
      R2.upper <- R2.fixed + qf((1 + confidence)/2, B.df, W.df) * 
        R2.random
      psrf <- cbind(sqrt(df.adj * R2.estimate), sqrt(df.adj * R2.upper))
      dimnames(psrf) <- list(xnames, c("Point est.", "Upper C.I."))
      out <- list(psrf = psrf, mpsrf = mpsrf)
      class(out) <- "gelman.diag"
      out
    }
    
    
    g <- matrix(NA, nrow=nvar(comb.samples), ncol=2)
    for (v in 1:nvar(comb.samples)) {
      g[v,] <- gelman.hack(comb.samples[,v], multivariate = T, autoburnin = F)$psrf
    }
    
    gelman <- as.data.frame(g)
    names(gelman)[names(gelman) == 'V1'] <- 'PSRF'
    names(gelman)[names(gelman) == 'V2'] <- 'PSRF: 97.5% quantile'
    
    # So staple the summary components and the gelman-rubin diagnostic together:
    posterior <- cbind(stats, quantiles, gelman)
    
    write.csv(posterior,
              file = file.path(output_dir,
                               paste0(sp, '_', 'it', max.it,
                                      '_ep', max(species_info$final_iteration),
                                      '.csv')))
    
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
