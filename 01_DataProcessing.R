########################################
## Data cleaning script
## Author: Francesca Mancini
## Date created: 2020-03-02
## Date modified: 2020-05-12
########################################

library(BRCmap)
library(dplyr)
library(ggplot2)
library(sparta)
library(openxlsx)
library(stringr)
library(tidyr)


## Load BWARS data ----
bwars <- read.csv("../Data/BWARS/March2020Update/20200305_CEH_Export.csv", 
                  header = F, stringsAsFactors = F)

# assign column names
names(bwars) <- c("CONCEPT", "TO_GRIDREF", "TO_STARTDATE", "TO_ENDDATE")

# format dates and calculate precision
bwars$TO_STARTDATE <- as.Date(bwars$TO_STARTDATE, format = "%Y-%m-%d")

bwars$TO_ENDDATE <- as.Date(bwars$TO_ENDDATE, format = "%Y-%m-%d")

bwars$DT_ID <- bwars$TO_ENDDATE - bwars$TO_STARTDATE

# delete all records that are not at day precision

bwars <- bwars %>%
  filter(DT_ID < 1)

# create year

bwars$YEAR <- format(bwars$TO_STARTDATE, "%Y")

table(bwars$YEAR)

# only keep data 2010-2016

bwars <- bwars %>%
  filter(YEAR >= 2010) %>%
  filter(YEAR < 2017)




## Load BeeWalks data ----

beewalks <- read.xlsx("../Data/BeeWalks/BeeWalk data 2008-17 01022018 for sharing.xlsx",
                      sheet = 1)

str(beewalks)

table(beewalks$DateType)
table(beewalks$Year)

# only keep records 2010 to 2016

beewalks <- beewalks %>%
  filter(Year >= 2010) %>%
  filter(Year < 2017)

# only keep records to species level
# records recorded at genus level are entered as
# for example: "Bombus sp."

beewalks <- beewalks %>%
  filter(!str_detect(latin, " sp."))

# calculate total number of individuals seen
# and filter out absences

beewalks <- beewalks %>%
  mutate(total = rowSums(.[34:37])) %>%
  filter(total > 0)

# only keep columns we need

beewalks <- beewalks %>%
  select(latin, GridReference, StartDate, Year)

# dates have been read as numbers
# openxlsx::convertToDate can deal with it

beewalks$StartDate <- convertToDate(beewalks$StartDate)


# visualise spatial overlap between the sites from the two data sources ----

# convert the BeeWalks data to 1Km res
beewalks$GRIDREF_1KM_PREC <- reformat_gr(beewalks$GridReference, 
                                         prec_out = 1000)

# filter BWARS data to remove records with precision < 1Km

bwars <- bwars %>%
  filter(nchar(bwars$TO_GRIDREF) > 5)
# convert the BWARS data to 1Km res
bwars$GRIDREF_1KM_PREC <- reformat_gr(bwars$TO_GRIDREF, 
                                      prec_out = 1000)


## DROP REPofIRELAND + NI HERE 
cn_info <- read.csv("../Data/sq1km_country_id_border_dropped.csv", header = TRUE)# add region data
bwars <- bwars[bwars$GRIDREF_1KM_PREC %in% cn_info[cn_info$ENGLAND == 1 | cn_info$WALES == 1 | cn_info$SCOTLAND == 1, "SQ1_SQUARE"],]



#plot first all the BWARS sites and then highlight those overlapping with BeeWalks in red
par(mar = c(0.1,0.1,1,0.1))
plot_GIS(UK$britain, new.window = FALSE, show.axis = FALSE, 
         show.grid = FALSE, xlab = "", ylab = "", 
         main = "BWARS sites overlapping with BeeWalks")

plotUK_gr(bwars$GRIDREF_1KM_PREC, gr_prec = 1000, border = "grey")
plotUK_gr(beewalks$GRIDREF_1KM_PREC, gr_prec = 1000, border = "red")


length(which(unique(beewalks$GRIDREF_1KM_PREC) %in% unique(bwars$GRIDREF_1KM_PREC)))
# there are 205 sites sampled in both datasets



## Clean species names ----

table(bwars$CONCEPT)

table(beewalks$latin)

## BWARS
## keep only bees in BWARS

# load bee species list and turn to lowercase
bee_sp <- read.csv("../Data/BWARS/bee_species_list.csv", 
                   header = T, stringsAsFactors = F)

bee_sp$species <- tolower(bee_sp$species)


# separate species name and taxonomy
bwars <- separate(bwars, "CONCEPT", into = c("CONCEPT", "TAXONOMY"), 
                  sep = ": ", remove = T, extra = "merge")

## turn all species names to lowercase
bwars$CONCEPT <- tolower(bwars$CONCEPT)

# filter bwars on bee species list
bwars <- bwars %>%
  filter(CONCEPT %in% bee_sp$species)

## ceck that taxonomy is the same for every secies
# 
# taxonomy <- bwars %>%
#   group_by(CONCEPT) %>%
#   summarise(taxonomies = n_distinct(TAXONOMY))


## BeeWalks
## turn all species names to lowercase
beewalks$latin <- tolower(beewalks$latin)


## delete brackets and content

beewalks$latin <- gsub("\\s*\\([^\\)]+\\)", "", beewalks$latin)

table(beewalks$latin)

write.csv(bwars, "./Data/BWARS_cleaned.csv", row.names = FALSE)

write.csv(beewalks, "./Data/BeeWalks_cleaned.csv", row.names = FALSE)


# How many species have records in both datasets?

common_species <- intersect(beewalks$latin, bwars$CONCEPT)
write.csv(data.frame(species = common_species), 
          "./Data/common_species.csv", row.names = FALSE)
length(common_species)



## Format data for JAGS ----

# y1 is the detection history for one species from the BWARS dataset
# y2 is the detection history fort he same species from the BeeWalks dataset

# Site1 is a list of sites visited in BWARS dataset (no unique)
# Site2 is a list of sites visited in BeeWalks (no unique)
# nsite is the total number of sites visited in both datasets (unique)

# nvisit1 is the numbe of visits in the BWARS dataset
# nvisit2 is the number of visits in the BeeWalks dataset

# DATATYPE2 is a column that indicates if that visit produced a medium list length (1) or not (0) (for BWARS dataset only)
# DATATYPE3 same as above for long list length

## BWARS

BWARS_formatted <- formatOccData(taxa = bwars$CONCEPT, 
                                 survey = bwars$TO_STARTDATE, 
                                 site = bwars$GRIDREF_1KM_PREC,
                                 closure_period = bwars$YEAR)

# save(BWARS_formatted, file = "./Data/BWARS_forJAGS.rdata")

## BeeWalks


beewalks_formatted <- formatOccData(taxa = beewalks$latin,
                                    survey = beewalks$StartDate,
                                    site = beewalks$GRIDREF_1KM_PREC,
                                    closure_period = beewalks$Year)

# save(beewalks_formatted, file = "./Data/BeeWalks_forJAGS.rdata")

