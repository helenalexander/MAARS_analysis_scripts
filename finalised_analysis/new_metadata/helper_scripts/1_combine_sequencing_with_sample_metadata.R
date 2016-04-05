#**********************************
# 1.1 Load sequencing metadata
#**********************************
sequenced_samples <- read.csv("input_files/sequenced_samples.csv", sep=";", stringsAsFactors=FALSE)
sequenced_samples$is_ntc <- FALSE

sequenced_negativeControls <- read.csv("input_files/sequenced_ntcs.csv", sep=";", stringsAsFactors=FALSE)
colnames(sequenced_negativeControls)[grep("SubmittedID",colnames(sequenced_negativeControls))] <- "MAARS_sample_id"
sequenced_negativeControls$is_ntc <- TRUE

#Combine both into a single dataframe
sequencing_data <- rbind(sequenced_samples,sequenced_negativeControls)
rm(sequenced_samples,sequenced_negativeControls)

#Convert selected columns to factors
factor_cols <- c("RunID")
for(col in factor_cols){
  sequencing_data[,col] <- as.factor(sequencing_data[,col]);
}

colnames (sequencing_data)
length(sequencing_data[,"ScilifeID"])
length(unique(sequencing_data[,"ScilifeID"]))


#***********************************************
# 1.2 Add filtered read counts
#***********************************************

load(file="input_files/humanremovalSummary.RData")
humanremovalSummary_df$filtered_reads = humanremovalSummary_df$readPairs - humanremovalSummary_df$pairsMapped

sequencing_data <- merge(x=sequencing_data, y=, select(humanremovalSummary_df,scilifeID,filtered_reads),
                  by.x="ScilifeID", by.y="scilifeID")

rm(humanremovalSummary_df)
#**************************************************************
# 1.2 Load summarized sample metadata
#**************************************************************
sample_metadata <- read.delim("input_files/MAARS_all_Fri_Apr_04_14h_CEST_2014-3.tsv", stringsAsFactors=FALSE)

#Convert selected columns to factors
factor_cols <- c("clinical_group",
                 "anatomical_location",
                 "anatomical_location_label",
                 "lesional",
                 "Institution",
                 "Gender",
                 "Global_Assessment_Score")

for(col in factor_cols){
  sample_metadata[,col] <- as.factor(sample_metadata[,col]);
}

#Note: clinical_group is set for all samples in the sample_metadata
#     Useful for determining which samples merged successfully
stopifnot( sample_metadata$clinical_group %>% is.na %>% sum() == 0 ) 

#***************************************************************
# 2.1 Join sequencing libraries with sequencing metadata
#***************************************************************
# Left outer join both datasets,  we want to annotate all sequenced samples as much as possible
# We also want to know which sequenced samples have no metadata

metadata <- merge(x=sequencing_data, y=sample_metadata, by.x="MAARS_sample_id", by.y="sample_id",all.x=TRUE )
colnames(metadata)

#***************************************************************
# 2.2 Wrangling of integrated data
#*****************************************************************

#**** Define a variable for unannotated data
metadata["excluded_from_sample_metadata"] <- is.na(metadata$clinical_group) & ( !metadata$is_ntc )

#******** Variable renaming ********************
# it's not gender, it's sex.
colnames(metadata)[grep("Gender", colnames(metadata))] <- "sex"

# just to simplify the name.
colnames(metadata)[grep("CUSTOM_Age", colnames(metadata))] <- "age"

# Clarify that MAARS_identifier refers to a patient specific id
colnames(metadata)[grep("MAARS_identifier", colnames(metadata))] <- "maars_patient_id"

# Standardize ScilifeID to scilife_id. Probably better to change to sequenced_sample_id or something
colnames(metadata)[grep("ScilifeID", colnames(metadata))] <- "scilife_id"

# Change RunID (from HiSeq Run) to batch_id
colnames(metadata)[grep("RunID", colnames(metadata))] <- "batch_id"

# Clarify MSequenced to million raw reads sequenced
colnames(metadata)[grep("MSequenced", colnames(metadata))] <- "million_reads_sequenced"

# Change "Fam._hist._Atopic_dermatitis" to match "Family_History_of_Psoriasis"
colnames(metadata)[grep("Fam._hist._Atopic_dermatitis", colnames(metadata))] <- "Family_History_of_Atopic_dermatitis"

#Infer maars_patient_id for unannotated samples (For later joining with patient metadata)
metadata$maars_patient_id[is.na(metadata$maars_patient_id) & !metadata$is_ntc  ] <- sapply(metadata$MAARS_sample_id[is.na(metadata$maars_patient_id) & !metadata$is_ntc ],
                                                                                            function(sample_id) substring(sample_id,1,nchar(sample_id)-3))
#************** Anatomical location wrangling

# Remove anatomical location label
# The label is less informative than the "anatomical_location"
metadata = metadata[,-grep("anatomical_location_label",colnames(metadata))]
colnames(metadata)

# Join Posterior thigh and Thigh into a single label
metadata[!is.na(metadata$anatomical_location) & (metadata$anatomical_location == "posterior_thigh"), "anatomical_location"] <-  "thigh"

# Join lower back and buttocks into a single label
metadata[!is.na(metadata$anatomical_location) & (metadata$anatomical_location == "buttocks"), "anatomical_location"] <-  "lower_back"

# Reduce the factor levels from 5 to 3
metadata$anatomical_location <- factor(metadata$anatomical_location,levels=c("thigh","lower_back","upper_back"))

#**************

#Remove the 'CUSTOM_' prefix from some of the variable names
grep("CUSTOM",colnames(metadata))
colnames(metadata) <- gsub("^CUSTOM_","",colnames(metadata))

#Remove the '_v2..' infix from some of the variable names
colnames(metadata)[grep("v2..",colnames(metadata))]
colnames(metadata) <- gsub("v2..","",colnames(metadata))
colnames(metadata)

#************** Negative Control variable generation

### Annotate Negative Controls(ntcs) with institution and buffer type infered from the ntc's MAARS_sample_id
# mock_1.xx PBS HHU
# mock_2.xx PBS KINGS
# mock_3.xx PBS UH
# mock_4.xx PBS KI from library prep (Scilife)
# H2O - Water KI from library prep (Scilife) (RunIDs 13_08, 14_05, 14_06, 14_07, 14_08, 15_01, 15_02, 15_04)

#Extract the institution of the pbs ntc from the first number of the id . H2O-named ntcs will remain unchanged
ntc_type <- sapply( tolower(metadata$MAARS_sample_id[metadata$is_ntc] ), function(x) gsub("mock_([1-4])(\\..+)?","\\1",x) )

#Create an additional variable to store which type of buffer were the ntc's kept in
metadata[metadata$is_ntc,"sample_buffer"] <- factor(ntc_type == "h2o",levels=c(FALSE,TRUE),labels=c("PBS","H2O"))
#add PBS as sample buffers to samples annotated from the sample_metadata
metadata[is.na(metadata$sample_buffer),"sample_buffer"] <- "PBS" 

#Assign institution based on the codes
ntc_institution <- factor(ntc_type,levels=c(1,2,3,4,"h2o"))
levels(ntc_institution) <- c("HHU","KINGS","UH","KI","KI")

#Add the KI level to the institution factor (for ntc's extracted at library prep stage) 
metadata$Institution <- factor(metadata$Institution,levels=c(levels(metadata$Institution),"KI"))
#Add Institution information to the ntcs
metadata[metadata$is_ntc,"Institution"] <- ntc_institution

#******************************************************
### Make all column names lowercase for standarization
#******************************************************
colnames(metadata) <- tolower(colnames(metadata))

#Clean up tmp variables
rm(ntc_type,ntc_institution)
rm(col,factor_cols)

#******************************************************
### Deal with duplicated samples
#******************************************************

# # {is it duplicated sequencing or mislabels?}
# duplicated_maars_sample_id <- (megadata.filt %>% filter(!is_ntc) %>% filter(duplicated(maars_sample_id)) %>% select(maars_sample_id))$maars_sample_id
# 
# duplicated.df <- megadata.filt %>% filter(maars_sample_id %in% duplicated_maars_sample_id) %>% select(batch_id,scilife_id,maars_sample_id,million_reads_sequenced) %>% arrange(maars_sample_id,batch_id)
# write.csv(duplicated.df,file = "duplicated_samples_maars_shotgun.csv")
# 
# str(megadata.filt)
# 
# #Verify how many patients have more or less than 2 sampless
# duplicated_patient_ids <- (megadata.filt %>% filter(!is_ntc) %>% count(maars_patient_id) %>% filter( n > 2) %>% arrange(desc(n)))$maars_patient_id
# suspicious_patient_ids <- (megadata.filt %>% filter(!is_ntc) %>% count(maars_patient_id) %>% filter( n < 2) %>% arrange(desc(n)))$maars_patient_id
# 
# #List Patients with more than 2 samples
# megadata.filt %>% filter(maars_patient_id %in% duplicated_patient_ids) %>% 
#   select(batch_id,scilife_id,maars_patient_id,maars_sample_id) %>% 
#   arrange(maars_patient_id)
# 
# ###NOTE!
# #Some of the maars_sample_id end in _12, should we filter them out?


