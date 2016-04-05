library(ggplot2)
library(stringr)
library(dplyr)

setwd("~/workspace/MAARS_analysis_scripts/finalised_analysis/new_metadata")

rm(list=ls())
#****************************************************************
# Combine Sequencing metadata with summarized per sample metadata
#****************************************************************
source("helper_scripts/1_combine_sequencing_with_sample_metadata.R")

#****************************************************************
# Validate combination
#****************************************************************
str(metadata)

for(col in colnames(metadata)){
  print(col)
  metadata[[col]] %>% summary %>% print
}

source("helper_scripts/2_validate_metadata_consistency.R")

#****************************************************************
# Combine generated metadata with Patient level metadata
#****************************************************************
source("helper_scripts/3_combine_metadata_with_patient_data.R")

#****************************************************************
# Validate combination
#****************************************************************
str(megadata)

for(col in colnames(megadata)){
  print(col)
  class(megadata[[col]]) %>% print
  megadata[[col]] %>% summary %>% print
}

source("helper_scripts/4_validate_megadata_consistency.R")


#*******************************************************************
# Final Wrangling / Filtering of megadata
#*******************************************************************

## Remove redundant columns
megadata <- megadata[,-grep("patient_clinical_group", colnames(megadata))]
megadata <- megadata[,-grep("clinical.assessment.global.assessment.score", colnames(megadata))]
megadata <- megadata[,-grep("sampled_lesional_skinsite", colnames(megadata))]

#Discard samples that are marked as excluded in patient metadata or excluded from the sample metadata file
megadata.filt <- megadata %>% filter(is_ntc | (!excluded & !excluded_from_sample_metadata))


for(col in colnames(megadata.filt)){
  print(col)
  class(megadata.filt[[col]]) %>% print
  megadata.filt[[col]] %>% summary %>% print
}

#*******************************************************************
# Find clinical samples with no sequenced samples
# Likely samples that had little DNA for sequencing or missing
# tubes
#*******************************************************************
# not_sequenced <- setdiff(  sample_metadata$sample_id,megadata$maars_sample_id)
# #Exclude samples wtih _03 and _04 that were not included in the study
# not_sequenced.filt <- not_sequenced[! grepl("_0(3|4|12)$",not_sequenced)]
# not_sequenced.filt


#*************************************************************************
# Finally, save everything to an Rdata and csv file
#*************************************************************************
dim(megadata)

dir.create("out/")
save(megadata,megadata.filt, file = "out/160405_megadata.Rdata")
write.csv(megadata, file = "out/160405_megadata.csv", quote=FALSE, row.names = FALSE)
