##----------------------------------------------------------------------------------------------
## Ruben Perez Perez
## Science Officer Data Center at VLIZ 
## EurOBIS - EMODnet Biology
## "Fri Oct  7 13:46:40 2022"

## Transforming this Pangaea dataset from GBIF into OBIS-ENV data format


### Data to be transformed stored and some background in https://jira.vliz.be/browse/EUROBIS-533

##----------------------------------------------------------------------------------------------


library(EMODnetBiocheck)
library(tidyverse)
library(obistools)
library(pangaear)

## Read data in from GBIF dataset
data <- read_delim("~/eurobis-dmt/datasets/pangaea/717138_data.tab", delim = "|")


## Cleaning
data <- data %>% mutate(identificationQualifier = case_when(grepl("sp.", scientificName) ~ "sp.",
                                                            TRUE ~ NA_character_),
                        scientificName = gsub("sp.", "", scientificName))


# Taxon match
taxa <-  data %>% select (scientificName) %>%
                  distinct() %>%
                  mutate(scientificNameID = match_taxa(scientificName)) %>%
                  do.call(what = data.frame) %>%
                  transmute(scientificName = scientificName,
                            scientificNameID = scientificNameID.scientificNameID)

data <- data %>% left_join(taxa, by = "scientificName")

#################################################################################################
# BLOCKER. This dataset doesn't have coordinates. Let's try a different method.

#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------


## Accessing the same dataset from pangaear R package

  
## Get data
res <- pg_data(doi = gsub("doi:", "", data$CollectionCode[1]))

## Get taxa from metadata
metadata <- res[[1]]$metadata$parameters 
  
df <- data.frame(matrix(unlist(metadata), byrow=TRUE),stringsAsFactors=FALSE) %>% 
      rename(x = matrix.unlist.metadata...byrow...TRUE.)

mtaxa <- df %>% filter(grepl("[#]", x)) %>%
                separate(col = x, sep = " \\[#] ", into = c("fullname", "shortname")) %>%
                mutate (shortname = gsub("\\)", "", shortname),
                        shortname = gsub("\\(", "", shortname)) %>%
                distinct() %>%
                mutate(scientificNameID = match_taxa(fullname)) %>%
                do.call(what = data.frame) %>%
                mutate(shortname = trimws(shortname))


## Integrate full taxon names into the data table
data2 <- res[[1]]$data %>% pivot_longer(cols = !1:7,
                                       names_to = "scientificName",
                                       values_to = "individualCount",
                                       names_repair = "unique") %>% 
                          mutate(measurementUnit = case_when(grepl("[#]", scientificName) ~ "[#]",
                                                                           TRUE ~ NA_character_),
                                 scientificName = gsub("\\[#]", "", scientificName),
                                 scientificName = sub(" \\(G.*)", "", scientificName),
                                 scientificName = trimws(scientificName)) %>%
                          left_join (mtaxa, 
                                     by = c("scientificName" = "shortname")) %>%
                          transmute (scientificName = ifelse(!is.na(scientificNameID.scientificName),
                                                          scientificNameID.scientificName,
                                                          as.character(fullname)),
                                  scientificNameID = ifelse(scientificNameID.match_type != "near_2", # Change this depending on the dataset or accept only match_type = exact
                                                            scientificNameID.scientificNameID,
                                                            NA_character_),
                                  collectionCode = res[[1]]$doi,
                                  basisOfRecord = "MaterialSample",
                                  occurrenceStatus = case_when(individualCount == 0 ~ "absent", 
                                                               individualCount != 0 ~ "present",
                                                               TRUE ~ NA_character_),
                                  datasetName = gsub(".*\\): (.+) PANGAEA.*", "\\1", res[[1]]$citation),
                                  institutionCode = "PANGAEA",  # Double check this, is PANGAEA the data creator?
                                  eventID = Event,
                                  occurrenceID = paste0(eventID, 
                                                        "_" ,
                                                        ifelse(!is.na(scientificNameID),
                                                               gsub("urn:lsid:marinespecies.org:taxname:", "", scientificNameID),
                                                               gsub(" ", "", scientificName))),
                                  eventDate = `Date/Time`,
                                  decimalLatitude = Latitude,
                                  decimalLongitude = Longitude,
                                  maximumDepthInMeters = gsub("-", "", `Elevation [m]`),
                                  minimumDepthInMeters = maximumDepthInMeters,
                                  `A. aculeata [#]...7`,
                                  `A. aculeata [#]...8`,
                                  individualCount) %>% 
                           distinct()
                                    




#################################################################################################################################
# Still needs to get Gear and person responsible from res/meta/params/ method/device and PI AND ship from res/meta/events/ Basis
#################################################################################################################################

# Create event table
event <- data2 %>% select (eventID, eventDate, institutionCode, datasetName, decimalLatitude, decimalLongitude, maximumDepthInMeters, minimumDepthInMeters) %>% distinct()

# Create occurrence table                                 
occurrence <- data2 %>% select (eventID, occurrenceID, scientificName, scientificNameID, basisOfRecord, occurrenceStatus, collectionCode, `A. aculeata [#]...7`, `A. aculeata [#]...8`) %>% distinct()

# Create emof table                                 
emof <- data2 %>% select (eventID, occurrenceID, individualCount) %>% distinct() %>%
                  pivot_longer(cols = !1:2,
                               names_to = "measurementType",
                               values_to = "measurementValue") %>%
                  mutate(measurementUnit = "Dimensionless",
                         measurementTypeID = "http://vocab.nerc.ac.uk/collection/P01/current/OCOUNT01/",
                         measurementValueID = NA_character_,
                         measurementUnitID = "http://vocab.nerc.ac.uk/collection/P06/current/UUUU/")



