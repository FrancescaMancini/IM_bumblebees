###################################################
## Script for running integrated model in JASMIN
## Author: Francesca Mancini
## Date created: 2020-03-10
## Date modified: 
###################################################


rm(list = ls())

# MCMC settings
ni <- 1000
nt <- 3
nb <- 0
nc <- 3



library(R2jags)
library(dplyr)
library(sparta)

#location of data input and output

data_dir <- "Data/"

output_folder <- "./Outputs/Integrated"


# Load the data -----

load(paste0(data_dir,"BWARS_forJAGS.rdata"))

load(paste0(data_dir,"BeeWalks_forJAGS.rdata"))


# Get job index
index <- as.numeric(commandArgs(trailingOnly = TRUE))[1]
str(index)


parameters <- read.csv(paste0(data_dir,"parameters_restart.csv"), stringsAsFactors = FALSE)[index,]
print(parameters)

set.seed(parameters$seed)


if(parameters$start == 1){
  
  # Data reformatting ----
  # y1 is the detection history for one species from the BWARS dataset
  # y2 is the detection history fort he same species from the BeeWalks dataset
  
  # Site1 is a list of sites visited in BWARS dataset (no unique)
  # Site2 is a list of sites visited in BeeWalks (no unique)
  # nsite is the total number of sites visited in both datasets (unique)
  
  # nvisit1 is the numbe of visits in the BWARS dataset
  # nvisit2 is the number of visits in the BeeWalks dataset
  
  # DATATYPE2 is a column that indicates if that visit produced a medium list length (1) or not (0) (for BWARS dataset only)
  # DATATYPE3 same as above for long list length
  
  # create Julian day variable from the date in visit
  # The variable visit at the moment is a combination of the 1Km grid ref and date
  # the code extracts the date from these character strings (characters number 7 to 16)
  # formats them as dates and converts them to numbers
  
  JulDate1 <- as.numeric(format(as.POSIXlt(substr(BWARS_formatted$spp_vis$visit, 7, 16),
                                           format = "%Y-%m-%d"), "%j"))
  
  Site1 <- BWARS_formatted$occDetdata$site
  Year1 <- as.numeric(factor(substr(BWARS_formatted$occDetdata$visit, 7, 10)))
  # summary(BWARS_formatted$occDetdata$L)
  DATATYPE2 <- BWARS_formatted$occDetdata %>%
    mutate(L = case_when(L < 2 ~ 0,
                         L < 4 ~ 1,
                         L >= 4 ~ 0)) %>%
    select(L)
  
  
  DATATYPE3 <- BWARS_formatted$occDetdata %>%
    mutate(L = case_when(L < 4 ~ 0,
                         L >= 4 ~ 1)) %>%
    select(L)
  
  
  # create Julian day variable from the date in visit
  # The variable visit at the moment is a combination of the 1Km grid ref and the date
  # the code extracts the date from these character strings (characters number 7 to 16)
  # formats them as dates and converts them to numbers
  JulDate2 <- as.numeric(format(as.POSIXlt(substr(beewalks_formatted$spp_vis$visit, 7, 16),
                                           format = "%Y-%m-%d"), "%j")) 
  
  
  Site2 <- beewalks_formatted$occDetdata$site
  Year2 <- as.numeric(factor(substr(beewalks_formatted$occDetdata$visit, 7, 10)))
  
  
  nsite <- length(unique(c(as.character(Site1), as.character(Site2))))
  
  
  
  write.table(paste("Integrated model for species",parameters$species,"has begun at ",date(),sep=" "),
              paste("./Outputs/progress_int.txt",sep=""),
              append=TRUE,col.names=FALSE,row.names=FALSE,quote=FALSE)
  
  # format the data - BWARS
  
  y1 <- as.integer(BWARS_formatted$spp_vis[, parameters$species]) 
  
  nvisit1 <- length(y1)
  
  # format the data - POMS
  y2 <- as.integer(beewalks_formatted$spp_vis[, parameters$species])
  
  nvisit2 <- length(y2)
  
  nyear <- length(unique(substr(BWARS_formatted$occDetdata$visit, 7, 10)))
  
  
  ## Integrated model - Julian data as gaussian density function ----
  
  jags_integrated_data <- list(y1 = y1, y2 = y2, Site1 = Site1, Site2 = Site2, nvisit1 = nvisit1,
                               nvisit2 = nvisit2, DATATYPE2 = DATATYPE2, DATATYPE3 = DATATYPE3, 
                               nsite = nsite, JulDate1 = JulDate1, JulDate2 = JulDate2,
                               Year1 = Year1, Year2 = Year2, nyear = nyear)
  
  inits_integrated <- function() list(a = rep(2, nyear), dtype1.p = rep(-0.5, nyear), dtype2.p = 1.1, 
                                      dtype3.p = 2, beewalks.p = rep(0.2, nyear), beta1 = 180, 
                                      beta2 = 16, beta3 = 0.1)
  
  params_integrated <- c("a", "dtype1.p", "dtype2.p", "dtype3.p", "beewalks.p", 
                         "psi.fs", "beta1", "beta2", "beta3")
  
  load.module("dic")
  
  out_integrated_model <- jags(data = jags_integrated_data, parameters.to.save = params_integrated,
                               model.file = "Integrated_model.txt", inits = inits_integrated,
                               n.chains=nc, n.iter=ni, n.thin=nt, n.burnin=nb)
  
  
  
  save(out_integrated_model, file = file.path(output_folder,
                                              paste0(paste(parameters$species,
                                                           paste0("it",
                                                                  parameters$end),
                                                           sep = '_'),
                                                     '.rdata')))
  
} else {
  
  source('Functions/daisy_chain_functions.R')
  
  # get the path to the previous output
  filename <- file.path(output_folder,
                        paste0(paste(parameters$species,
                                     paste0("it", 
                                            parameters$start - 1),
                                     sep = '_'),
                               '.rdata'))
  # run the daisy
  out_integrated_model <- daisyChain(rdataFile = filename,
                                     outDir = output_folder,
                                     nDaisy = 1,
                                     by.it = length(parameters$start:parameters$end),
                                     n.thin = 3,
                                     quiet = TRUE, 
                                     model = "integrated", 
                                     update = TRUE)
  
  save(out_integrated_model, file = file.path(output_folder,
                                              paste0(paste(parameters$species,
                                                           paste0("it",
                                                                  parameters$end),
                                                           sep = '_'),
                                                     '.rdata')))
  
}


quit(save = 'no', runLast = FALSE)
