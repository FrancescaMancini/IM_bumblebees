
library(R2jags)



daisyChain <- function(rdataFile, nDaisy = 10, outDir = getwd(), 
                       by.it = 100, n.thin = 3, quiet = TRUE, 
                       model = "integrated", update = TRUE){
  
  #  tot <- 0
  
  # This while statement could be changed to assess some other
  # parameter such as Rhat
  
  for(i in 1:nDaisy){
    if(model == "integrated"){
      cat("Daisy", i, "\n", file = paste0(outDir, "progress.txt"), append = T)
      
      chain_out <- daisy_integrated(rdataFile = rdataFile, outDir = outDir, 
                                    n.iter = by.it, n.thin = n.thin,
                                    quiet = TRUE, update = update)
      #    tot <- chain_out$tot.it
      rdataFile <- chain_out$out_integrated_model
      
      
      # return(invisible(chain_out$out_integrated_model))
      
      
    } else {
      chain_out <- daisy_not_integrated(rdataFile = rdataFile, outDir = outDir, 
                                        n.iter = by.it, n.thin = n.thin,
                                        quiet = TRUE, update = update)
      #    tot <- chain_out$tot.it
      rdataFile <- chain_out$out_BWARS_model
      
      
      # return(chain_out$out_BWARS_model)
      
      
    }
  }
  return(invisible(rdataFile))
}


daisy_integrated <- function(rdataFile, outDir = getwd(),
                             n.iter = 100, n.thin = 3,
                             quiet = TRUE, update){
  
  
  
  if(is.character(rdataFile)){
    load(file = rdataFile) #out
    obj <- FALSE
  } else { #its an object
    out_integrated_model <- rdataFile
    if('filename' %in% names(out_integrated_model)){
      rdataFile <- out_integrated_model$filename  
    } else { # make up a name
      rdataFile <- paste0(paste(out_integrated_model$SPP_NAME,
                                gsub('.txt$', '', tolower(basename(out_integrated_model$model.file))),
                                sep = '_'),
                          'integrated.rdata')
      rdataFile <- file.path(outDir, rdataFile)
    }
    
    obj <- TRUE
    
    
  }
  
  cat('\nLoaded', out_integrated_model$n.iter, 'iterations... ', 
      file = paste0(outDir, "progress_int.txt"), append = TRUE)
  
  # if((out_integrated_model$n.iter + n.iter) > total.it) n.iter <- total.it - out_integrated_model$n.iter
  
  # didn't work before this module is loaded - sparta versions include deviance which needs this
  load.module("dic")
  
  # If we have loaded the file we need to recompile it.
  if(update == TRUE){
    if(!obj){
      if(quiet){
        cat('compiling... ')
        invisible(capture.output(R2jags:::recompile.rjags(out_integrated_model, progress.bar = 'none')))
      } else {
        R2jags:::recompile.rjags(out_integrated_model)
        cat('\nRecompiled', file = paste0(outDir, "progress.txt"), append = TRUE)
        
      }
    }
    
    
    new_out_integrated_model <- R2jags:::update.rjags(out_integrated_model,
                                                      n.iter = n.iter,
                                                      n.thin = n.thin)
    
    cat('\nUpdated', file = paste0(outDir, "progress.txt"), append = TRUE)
    
    new_out_integrated_model$SPP_NAME <- out_integrated_model$SPP_NAME
    new_out_integrated_model$min_year <- out_integrated_model$min_year
    new_out_integrated_model$max_year <- out_integrated_model$max_year
    out_integrated_model <- new_out_integrated_model
    
  }
  
  
  cat('output', out_integrated_model$n.iter, '\n', 
      file = paste0(outDir, "progress.txt"), append = TRUE)
  
  return(list(out_integrated_model = out_integrated_model))
}


daisy_not_integrated <- function(rdataFile, outDir = getwd(),
                                 n.iter = 100, n.thin = 3,
                                 quiet = TRUE, update){
  
  
  
  if(is.character(rdataFile)){
    load(file = rdataFile) #out
    obj <- FALSE
  } else { #its an object
    out_BWARS_model <- rdataFile
    if('filename' %in% names(out_BWARS_model)){
      rdataFile <- out_BWARS_model$filename  
    } else { # make up a name
      rdataFile <- paste0(paste(out_BWARS_model$SPP_NAME,
                                gsub('.txt$', '', tolower(basename(out_BWARS_model$model.file))),
                                sep = '_'),
                          'HRS.rdata')
      rdataFile <- file.path(outDir, rdataFile)
    }
    
    obj <- TRUE
    
    
  }
  
  cat('\nLoaded', out_BWARS_model$n.iter, 'iterations... ', 
      file = paste0(outDir, "progress.txt"), append = TRUE)
  
  # if((out_BWARS_model$n.iter + n.iter) > total.it) n.iter <- total.it - out_BWARS_model$n.iter
  
  # didn't work before this module is loaded - sparta versions include deviance which needs this
  load.module("dic")
  
  # If we have loaded the file we need to recompile it.
  if(update == TRUE){
    if(!obj){
      if(quiet){
        cat('compiling... ')
        invisible(capture.output(R2jags:::recompile.rjags(out_BWARS_model, progress.bar = 'none')))
      } else {
        R2jags:::recompile.rjags(out_BWARS_model)
        cat('\nRecompiled', file = paste0(outDir, "progress.txt"), append = TRUE)
        
      }
    }
    
    
    new_out_BWARS_model <- R2jags:::update.rjags(out_BWARS_model,
                                                 n.iter = n.iter,
                                                 n.thin = n.thin)
    
    cat('\nUpdated', file = paste0(outDir, "progress.txt"), append = TRUE)
    
    new_out_BWARS_model$SPP_NAME <- out_BWARS_model$SPP_NAME
    new_out_BWARS_model$min_year <- out_BWARS_model$min_year
    new_out_BWARS_model$max_year <- out_BWARS_model$max_year
    out_BWARS_model <- new_out_BWARS_model
    
  }
  
  
  cat('output', out_BWARS_model$n.iter, '\n', 
      file = paste0(outDir, "progress.txt"), append = TRUE)
  
  return(list(out_BWARS_model = out_BWARS_model))
}


#files_list <- paste0(getwd(), "/Output/Rerun/Integrated/", list.files("./Output/Rerun/Integrated/"))
#
#hosts<-as.character(read.table(Sys.getenv('PBS_NODEFILE'),header=FALSE)[,1]) # read the nodes to use
#
#sfSetMaxCPUs(length(hosts)) # ensure that snowfall can cope with this many hosts
#
#sfInit(parallel=TRUE,type="MPI",cpus=length(hosts),useRscript=TRUE) # initialise the connection
#
#sfExportAll()
#sfLibrary(R2jags)
#
#sfClusterApplyLB(files_list, daisyChain, nDaisy = 10, 
#                 outDir = "./Output/Rerun/Integrated", by.it = 20000, n.thin = 10, 
#                 quiet = TRUE, model = "integrated", update = TRUE)
#
#sfStop()
#
#
#files_list <- paste0(getwd(), "/Output/Rerun/HRS/", list.files("./Output/Rerun/HRS/"))
#
#hosts<-as.character(read.table(Sys.getenv('PBS_NODEFILE'),header=FALSE)[,1]) # read the nodes to use
#
#sfSetMaxCPUs(length(hosts)) # ensure that snowfall can cope with this many hosts
#
#sfInit(parallel=TRUE,type="MPI",cpus=length(hosts),useRscript=TRUE) # initialise the connection
#
#sfExportAll()
#sfLibrary(R2jags)
#
#sfClusterApplyLB(files_list, daisyChain, nDaisy = 10, 
#                 outDir = "./Output/Rerun/HRS", by.it = 20000, n.thin = 10, 
#                 quiet = TRUE, model = "HRS", update = TRUE)
#
#sfStop()
#


