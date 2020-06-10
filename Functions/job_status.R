# Job status
# use files to work out progress through a job, identify jobs that have failed mid-chain
# dir - top directory containing console files with sub-dir 'output' inc. .rdata files
# target.it - The number of iterations you are aiming for default to infinite
# res_summary - is a resource summary output required.  This is optional as this step can be very slow.

job_status <- function(dir = getwd(), target.it = Inf, first.it = 0,
                       roster_file = "roster.csv", datafiles_path = "output/",
                       logfiles_path = "console_output", consolefiles_path = "console_output/console"){
  
  require(reshape2)
  
  # Load information
  roster <- read.csv(file = file.path(dir, roster_file), stringsAsFactors = FALSE)
  data_files <- list.files(file.path(dir, datafiles_path), pattern = '.rdata$')
  console_files <- list.files(file.path(dir, consolefiles_path), pattern = '^console')
  
  # Add column for console output to roster
  roster$console <- FALSE
  roster$console[as.numeric(gsub('^console', '', gsub('.Rout$', '', console_files)))] <- TRUE
  
  # Add column for data output to roster
  expected_data <- paste0(roster$species, '_it', roster$end,'.rdata')
  roster$data <- expected_data %in% data_files
  
  # Find species that have failed
  console_counts <- melt(tapply(roster$console, roster$species, FUN = sum))
  data_counts <- melt(tapply(roster$data, roster$species, FUN = sum))
  diff_counts <- console_counts$value - data_counts$value
  jobs_failed <- console_counts[diff_counts > 1, c(1:2)]
  colnames(jobs_failed) <- c('species', 'chain')
  nfailed <- nrow(jobs_failed)
  
  summary_table <- data.frame(species = console_counts$Var1)
  
  # Add max iteration for each job
  for(sp in unique(summary_table$species)){
    
    sp_ch_roster <- roster[roster$species == sp & roster$data, ]
    if(nrow(sp_ch_roster) == 0){
      max.it = first.it
    } else {
      max.it <- max(sp_ch_roster$end)
    }
    summary_table$max.it[summary_table$species == sp] <- max.it
    
  }
  
  
  summary_table$status <- 'running'
  
  for(i in 1:nrow(jobs_failed)){
    
    summary_table$status[summary_table$species == jobs_failed$species[i]] <- 'failed'
    
  }
  
  # mark as complete those that have reached the target
  summary_table$status[summary_table$max.it >= target.it] <- 'Complete'
  
  # order summary table by max.it
  summary_table <- summary_table[order(summary_table$max.it), ]
  
  # Get runtime and memory summaries
  o_files <- list.files(file.path(dir, logfiles_path), pattern = '\\.[o]$', full.names = TRUE)
  
  mem <- NULL
  time <- NULL
  
  mem_time <- sapply(o_files, FUN = function(x){
    
    temp_f <- trimws(readLines(x))
    CPU_time <- as.numeric(trimws(gsub('^CPU time :', '', gsub('sec\\.$', '', temp_f[grepl('^CPU time', temp_f)]))))
    max_mem <- as.numeric(trimws(gsub('^Max Memory :', '', gsub('MB$', '', temp_f[grepl('^Max Memory', temp_f)]))))
    
    if(is.na(max_mem)){
      max_mem <- 0
    }
    
    if(length(CPU_time) != 1 | length(max_mem) != 1){
      mem <- c(mem, NA)
      time <- c(time, NA)
    } else {
      mem <- c(mem, max_mem)
      time <- c(time, CPU_time)
    }
    
    return(c(mem, time))
    
  })
  
  resource_summary <- rbind(summary(mem_time[1,]), summary((mem_time[2,]/60)/60))
  row.names(resource_summary) <- c('max_memory_MB', 'time_hrs')
  
  return(list(summary_table = summary_table,
              resource_summary = resource_summary))
}
