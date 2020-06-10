# A function to look at the failed runs and set up a job script 
# to restart the needed jobs
#  You should be passing this function a minimum of 'dir', 'originalJob' and
#    probably 'res_summary = FALSE' unless you're feeling patient
# dir - top directory of the job, where the roster file resides
# originalJob - the original job file
# remove_logs - should the old logs be removed? Needed if the summary is going
#   to work in the future
# res_summary - is a resource summary output required.  This is optional as this
#   step can be very slow.
#
#  Example usage:
# restart_failures(dir = '../../data_degradation',
#                  originalJob = '../../data_degradation/bsub_lines.sh')

restart_failures <- function(dir = getwd(), originalJob = NULL, remove_logs = TRUE, datafiles_path = "output/", consolefiles_path = "console_output",
                             logfiles_path = "console_output", roster_file = "roster.csv", q = NULL, walltime = NULL, res_summary = FALSE,
                             job_file = "multi_array.job", first.it = 0, target.it = 50000){
  
  # Get idea of failures
  source('Functions/job_status.R')
  
  j_stat <- job_status(dir = dir, roster_file = roster_file, datafiles_path = datafiles_path, first.it = first.it,
                       logfiles_path = logfiles_path, consolefiles_path = consolefiles_path, 
                       target.it =  target.it)
  failures <- j_stat$summary_table[j_stat$summary_table$status == 'failed', ]
  
  # Get the roster
  roster <- read.csv(file = file.path(dir, roster_file), stringsAsFactors = FALSE)
  
  # Get parameters from the original
  if(is.null(originalJob)) stop('No job script passed to function')
  originalJob_eg <- readLines(originalJob, warn=FALSE)[1]
  if(is.null(q)) q <- gsub('\\-q\\s', '', stringr::str_extract(originalJob_eg, '\\-q\\s[[:alpha:][:punct:]]+'))
  MbMemory <- as.numeric(gsub('mem=', '', stringr::str_extract(originalJob_eg, 'mem=[[:digit:]]+')))
  if(is.null(walltime)) walltime <- gsub('\\-W\\s', '', stringr::str_extract(originalJob_eg, '\\-W\\s[[:digit:][:punct:]]+'))
  jobName <- stringr::str_extract(originalJob_eg, '(?<=-J\\s).+(?=\\[)')
  
  doc <- NULL
  consolefiles <- list.files(file.path(dir, consolefiles_path), pattern = '^console',
                             recursive = TRUE, full.names = TRUE)
  oe_files <- list.files(file.path(dir, logfiles_path), pattern = "\\.[oe]$",
                         recursive = TRUE, full.names = TRUE)
  
  if(nrow(failures) == 0) stop('There are no failures')
  
  for(fail in 1:nrow(failures)){
    # Get the bit of roster to re-run
    roster_sec <- roster[roster$species == failures$species[fail] &
                           roster$end > failures$max.it[fail], ]
    if(!is.null(roster$file)){
      roster_sec <-
        roster_sec[roster_sec$file == paste0('records_',failures$source_file[fail],'.rdata'), ]
      ds <- failures$source_file[fail]
    } else {
      ds <- ''
    }
    
    array_s_e <- range(as.numeric(row.names(roster_sec)))
    
    # write the sh script
    job_line <- paste0('-J ', jobName, '[', array_s_e[1], '-', array_s_e[2], ']%1')
    submit_line <- paste('bsub', # submit
                         paste('-q', q), #queue
                         '-n 1', #N cores,
                         '-R', paste0('"rusage[mem=', MbMemory, ']"'),
                         '-M', as.integer(MbMemory),
                         paste('-W', walltime), #wall time
                         job_line, #job array
                         paste0('-oo ', logfiles_path, '/',ds,'R-%J-%I.o'), #log
                         paste0('-eo ', logfiles_path, '/',ds,'R-%J-%I.e'), #log
                         paste0('< ', job_file)) #call to R
    
    doc <- c(doc, submit_line)
    
    if(remove_logs){
      
      # delete failed console files
      con_del_index <- as.numeric(gsub('^console', '', gsub('.Rout$', '', basename(consolefiles)))) %in% array_s_e[1]:array_s_e[2]
      unlink(x = consolefiles[con_del_index])
      
      # delete failed .o and .e files
      oe_del_index <- as.numeric(gsub('\\.[oe]', '', stringr::str_extract(basename(oe_files), '[:digit:]+\\.[oe]$'))) %in% array_s_e[1]:array_s_e[2]
      unlink(x = oe_files[oe_del_index])
      
    }
  }
  
  time_now <- format(Sys.time(),"%Y%m%d_%H%M%S")
  write.table(doc,
              file = file(paste0(dir, "/bsub_lines_", time_now, '.sh'), "wb"),
              append = TRUE, col.names = FALSE,
              row.names = FALSE, quote = FALSE)
  closeAllConnections()
  
  return(paste0(dir, "bsub_lines_", time_now, '.sh'))
}