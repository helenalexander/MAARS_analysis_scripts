# Merge clinical metadata (of samples) with sequencing metadata (= clin_seq_metadata)
# Excluded patients and patients with absent metadata
# Merge clin_seq_metadata with patient metadata ( = clin_seq_pt_metadata)
# Adding negative controls

rm(list=ls())

### load list of samples in the BRC Bioinformatics Cluster
 listOfFastq = read.delim2("/Users/helenalexander/Documents/MAARS/data/metaphlan/analysis/metadata/listOfSamples.txt", header=FALSE)
 tmpFun=function(y){tmp1 = strsplit(x = as.character(y), split = "/");   return(tail(tmp1[[1]], n=1))}
 listOfScilifeID_0 = sapply(X=(as.character(listOfFastq$V1)), FUN=tmpFun)
 tmpFun=function(y){tmp1 = strsplit(x = as.character(y), split = "_index");   return(tmp1[[1]][1])}
 listOfScilifeID = sapply(X=listOfScilifeID_0, FUN=tmpFun)

#**********************************
# Load clinical metadata
#**********************************
clinical_data <- read.delim("/Users/helenalexander/Documents/MAARS/data/metaphlan/analysis/metadata/MAARS_all_Fri_Apr_04_14h_CEST_2014-3.tsv", stringsAsFactors=FALSE)
cnames = colnames(clinical_data)
write.csv(x = cnames, file = "/Users/helenalexander/Documents/MAARS/data/metaphlan/analysis/metadata/clinicaldatacolnames.csv")

#Convert selected columns to factors
factor_cols <- c("clinical_group",
                 "anatomical_location",
                 "anatomical_location_label",
                 "lesional",
                 "MAARS_identifier",
                 "Institution",
                 "Gender",
                 "Global_Assessment_Score")

for(col in factor_cols){
  clinical_data[,col] <- as.factor(clinical_data[,col]);
}

#**********************************
# Load sequencing metadata
#**********************************
sequencing_data <- read.csv("/Users/helenalexander/Documents/MAARS/data/metaphlan/analysis/metadata/sequenced_samples.csv", sep=";", stringsAsFactors=FALSE)
sequenced_negativeControls <- read.csv("/Users/helenalexander/Documents/MAARS/data/metaphlan/analysis/metadata/sequenced_ntcs.csv", sep=";", stringsAsFactors=FALSE)

sequencing_data = data.frame(RunID=c(sequencing_data$RunID, sequenced_negativeControls$RunID),
                             ScilifeID=c(sequencing_data$ScilifeID, sequenced_negativeControls$ScilifeID),
                             MAARS_sample_id=c(sequencing_data$MAARS_sample_id, sequenced_negativeControls$SubmittedID),
                             MSequenced=c(sequencing_data$MSequenced, sequenced_negativeControls$MSequenced),
                             BarcodeSeq=c(sequencing_data$BarcodeSeq, sequenced_negativeControls$BarcodeSeq)
)

#Convert selected columns to factors
factor_cols <- c("RunID")
for(col in factor_cols){
  sequencing_data[,col] <- as.factor(sequencing_data[,col]);
}

colnames (sequencing_data)
length(sequencing_data[,"ScilifeID"])
length(unique(sequencing_data[,"ScilifeID"]))

#*********************************
#Load clinical patient data
#*********************************
patient_data.ctrl <- read.table("/Users/helenalexander/Documents/MAARS/data/metaphlan/analysis/metadata/MAARS_Control_full_20150424_15-38-30-4.csv",
                                sep="\t",header=TRUE,quote="\"",check.names=TRUE,comment.char="",stringsAsFactors=FALSE)
patient_data.ctrl$clinical_group="CTRL"

patient_data.ad <- read.table("/Users/helenalexander/Documents/MAARS/data/metaphlan/analysis/metadata/MAARS_AD_full_20150424_15-37-35-3.csv",
                              sep="\t",header=TRUE,quote="\"",check.names=TRUE,comment.char="",stringsAsFactors=FALSE)
patient_data.ad$clinical_group="AD"

patient_data.pso <- read.table("/Users/helenalexander/Documents/MAARS/data/metaphlan/analysis/metadata/MAARS_PSO_full_20150423_18-43-57-3.csv",
                               sep="\t",header=TRUE,quote="\"",check.names=TRUE,comment.char="",stringsAsFactors=FALSE)
patient_data.pso$clinical_group="PSO"

#Make variable names equal throughout patient metadata files for Smoking
colnames(patient_data.pso)[18] = "Smoking"
colnames(patient_data.ctrl)[11] = "Smoking"
colnames(patient_data.ad)[10] = "Smoking"

#Create variable for "YearsSinceDiagnosis" for Psoriasis = YearsSinceDiagPso 
# (= "patient.Identification.Date.of.first.visit" - "patient.Psoriasis.Patients.Year.of.Diagnosis")
# First convert date first visit string to date
dateFV = as.Date(patient_data.pso[ ,2],format = "%d/%m/%y")
#Convert dateFV to YearFV
YearFV = format(dateFV, "%Y")
#Convert YearFV to Integer
IntYearFV = strtoi(YearFV, base = 0L)
#Create varaible YearD = 'Year of diagnosis'
YearD = patient_data.pso[ ,14]
#Calculate YsinceD
YearsSinceDiagPso = IntYearFV - YearD
patient_data.pso[,"YearsSinceDiagPso"] = YearsSinceDiagPso

#Create 'dummy' columns for variables not present in all patient groups
patient_data.ctrl[,"YearsSinceDiagPso"] = NA
patient_data.ctrl[,"patient.HEAD.HEAD...Total..Sum.x.Area.x0.1."] = NA
patient_data.ctrl[,"patient.TRUNK.TRUNK...Total..Sum.x.Area.x0.3."] = NA
patient_data.ctrl[,"patient.UPPER.LIMBS.UPPER.LIMBS...Total..Sum.x.Area.x0.2."] = NA
patient_data.ctrl[,"patient.LOWER.LIMBS.LOWER.LIMBS...Total..Sum.x.Area.x0.4."] = NA
patient_data.ctrl[,"patient.PASI.Score.Total..TotalHead...TotalTrunk...TotalUpperLimbs...TotalLowerLimbs..PASI.Score."] = NA
patient_data.ctrl[,"patient.Psoriasis.Severity.Index..PSI..PSI...Erythema"]= NA
patient_data.ctrl[,"patient.Psoriasis.Severity.Index..PSI..PSI...Induration"] = NA
patient_data.ctrl[,"patient.Psoriasis.Severity.Index..PSI..PSI...Scaling"] = NA
patient_data.ctrl[,"patient.Psoriasis.Severity.Index..PSI..PSI.Total"] = NA
patient_data.ctrl[,"patient.Evidence.of.superinfection..Non.lesional.skin"] = NA
patient_data.ctrl[,"patient.Please.state.the.areas.used.for.sampling.Lesional.skin"] = NA
patient_data.ctrl [,"patient.Clinical.Assessment.Global.Assessment.Score"] = NA

patient_data.ctrl[,"patient.Hanifin.and.Rajka.diagnostic.criteria.6.Palmar.hyperlinearity..Keratosis.pilaris"] = NA
patient_data.ctrl[,"patient.Hanifin.and.Rajka.diagnostic.criteria.7.1.Total.IgE..IU.ml."] = NA
patient_data.ctrl[,"patient.Hanifin.and.Rajka.diagnostic.criteria.7.1.a..150.IU.ml...x...400.IU.ml"] = NA
patient_data.ctrl[,"patient.Hanifin.and.Rajka.diagnostic.criteria.7.1.b...400.IU.ml"] = NA
patient_data.ctrl[,"patient.Hanifin.and.Rajka.diagnostic.criteria.8..Raised.serum.IgE"] = NA
patient_data.ctrl[,"patient.Hanifin.and.Rajka.diagnostic.criteria.9..Early.age.of.onset"] = NA
patient_data.ctrl[,"patient.Hanifin.and.Rajka.diagnostic.criteria.10..ECP.measurement..ug.ml."] = NA
patient_data.ctrl[,"patient.Hanifin.and.Rajka.diagnostic.criteria.11..Tendency.toward.cutaneous.infections.or.impaired.cell.mediated.immunity"] = NA
patient_data.ctrl[,"patient.Hanifin.and.Rajka.diagnostic.criteria.16..Food.intolerance"] = NA
patient_data.ctrl[,"patient.Hanifin.and.Rajka.diagnostic.criteria.33.a..Personal.history.of."] = NA
patient_data.ctrl[,"patient.33.b..Family.history.of.atopy.Fam..hist..Atopic.dermatitis"] = NA
patient_data.ctrl[,"patient.SCORAD.Objective.SCORAD..A.5...7B.2."] = NA
patient_data.ctrl[,"patient.SCORAD.SCORAD.Score"] = NA
patient_data.ctrl[,"patient.A..Please.state.the.areas.used.for.sampling.Lesional.skin"] = NA
patient_data.ctrl[,"patient.B..SCORAD..Intensity..SCORAD..Erythema"] = NA
patient_data.ctrl[,"patient.B..SCORAD..Intensity..SCORAD..Edema.Papulation"] = NA
patient_data.ctrl[,"patient.B..SCORAD..Intensity..SCORAD..Oozing.crusts"] = NA
patient_data.ctrl[,"patient.B..SCORAD..Intensity..SCORAD..Excoriation"] = NA
patient_data.ctrl[,"patient.B..SCORAD..Intensity..SCORAD..Lichenification"] = NA
patient_data.ctrl[,"patient.B..SCORAD..Intensity..SCORAD..Dryness..evaluated.on.uninvolved.areas."] = NA
patient_data.ctrl[,"patient.B..LOCAL..Intensity..LOCAL..Erythema"] = NA
patient_data.ctrl[,"patient.B..LOCAL..Intensity..LOCAL..Edema.Papulation"] = NA                                                                                   
patient_data.ctrl[,"patient.B..LOCAL..Intensity..LOCAL..Oozing.crusts"] = NA                                                                                      
patient_data.ctrl[,"patient.B..LOCAL..Intensity..LOCAL..Excoriation"] = NA                                                                                        
patient_data.ctrl[,"patient.B..LOCAL..Intensity..LOCAL..Lichenification"] = NA                                                                                    
patient_data.ctrl[,"patient.B..LOCAL..Intensity..LOCAL..Dryness..evaluated.on.uninvolved.areas."] = NA                                                            
patient_data.ctrl[,"patient.B..LOCAL..Intensity.LOCAL.SCORAD.Score"] = NA

patient_data.ad[,"YearsSinceDiagPso"] = NA
patient_data.ad[,"patient.HEAD.HEAD...Total..Sum.x.Area.x0.1."] = NA
patient_data.ad[,"patient.TRUNK.TRUNK...Total..Sum.x.Area.x0.3."] = NA
patient_data.ad[,"patient.UPPER.LIMBS.UPPER.LIMBS...Total..Sum.x.Area.x0.2."] = NA
patient_data.ad[,"patient.LOWER.LIMBS.LOWER.LIMBS...Total..Sum.x.Area.x0.4."] = NA
patient_data.ad [,"patient.PASI.Score.Total..TotalHead...TotalTrunk...TotalUpperLimbs...TotalLowerLimbs..PASI.Score."] = NA
patient_data.ad[,"patient.Psoriasis.Severity.Index..PSI..PSI...Erythema"]= NA
patient_data.ad[,"patient.Psoriasis.Severity.Index..PSI..PSI...Induration"] = NA
patient_data.ad[,"patient.Psoriasis.Severity.Index..PSI..PSI...Scaling"] = NA
patient_data.ad[,"patient.Psoriasis.Severity.Index..PSI..PSI.Total"] = NA
patient_data.ad [,"patient.Please.state.the.areas.used.for.sampling.Lesional.skin"] = NA

patient_data.pso[,"patient.Hanifin.and.Rajka.diagnostic.criteria.6.Palmar.hyperlinearity..Keratosis.pilaris"] = NA
patient_data.pso[,"patient.Hanifin.and.Rajka.diagnostic.criteria.7.1.Total.IgE..IU.ml."] = NA
patient_data.pso[,"patient.Hanifin.and.Rajka.diagnostic.criteria.7.1.a..150.IU.ml...x...400.IU.ml"] = NA
patient_data.pso[,"patient.Hanifin.and.Rajka.diagnostic.criteria.7.1.b...400.IU.ml"] = NA
patient_data.pso[,"patient.Hanifin.and.Rajka.diagnostic.criteria.8..Raised.serum.IgE"] = NA
patient_data.pso[,"patient.Hanifin.and.Rajka.diagnostic.criteria.9..Early.age.of.onset"] = NA
patient_data.pso[,"patient.Hanifin.and.Rajka.diagnostic.criteria.10..ECP.measurement..ug.ml."] = NA
patient_data.pso[,"patient.Hanifin.and.Rajka.diagnostic.criteria.11..Tendency.toward.cutaneous.infections.or.impaired.cell.mediated.immunity"] = NA
patient_data.pso[,"patient.Hanifin.and.Rajka.diagnostic.criteria.16..Food.intolerance"] = NA
patient_data.pso[,"patient.Hanifin.and.Rajka.diagnostic.criteria.33.a..Personal.history.of."] = NA
patient_data.pso[,"patient.33.b..Family.history.of.atopy.Fam..hist..Atopic.dermatitis"] = NA
patient_data.pso[,"patient.SCORAD.Objective.SCORAD..A.5...7B.2."] = NA
patient_data.pso[,"patient.SCORAD.SCORAD.Score"] = NA
patient_data.pso[,"patient.A..Please.state.the.areas.used.for.sampling.Lesional.skin"] = NA
patient_data.pso[,"patient.B..SCORAD..Intensity..SCORAD..Erythema"] = NA
patient_data.pso[,"patient.B..SCORAD..Intensity..SCORAD..Edema.Papulation"] = NA
patient_data.pso[,"patient.B..SCORAD..Intensity..SCORAD..Oozing.crusts"] = NA
patient_data.pso[,"patient.B..SCORAD..Intensity..SCORAD..Excoriation"] = NA
patient_data.pso[,"patient.B..SCORAD..Intensity..SCORAD..Lichenification"] = NA
patient_data.pso[,"patient.B..SCORAD..Intensity..SCORAD..Dryness..evaluated.on.uninvolved.areas."] = NA
patient_data.pso[,"patient.B..LOCAL..Intensity..LOCAL..Erythema"] = NA
patient_data.pso[,"patient.B..LOCAL..Intensity..LOCAL..Edema.Papulation"] = NA                                                                                   
patient_data.pso[,"patient.B..LOCAL..Intensity..LOCAL..Oozing.crusts"] = NA                                                                                      
patient_data.pso[,"patient.B..LOCAL..Intensity..LOCAL..Excoriation"] = NA                                                                                        
patient_data.pso[,"patient.B..LOCAL..Intensity..LOCAL..Lichenification"] = NA                                                                                    
patient_data.pso[,"patient.B..LOCAL..Intensity..LOCAL..Dryness..evaluated.on.uninvolved.areas."] = NA                                                            
patient_data.pso[,"patient.B..LOCAL..Intensity.LOCAL.SCORAD.Score"] = NA

#Join all patient data + excluded information
cols_to_get <- c("patient.Identification.MAARS.identifier", "clinical_group",
                 "patient.Curated.Items.Excluded","patient.Curated.Items.Exclusion.comments",
                 "patient.Ethnicity.Family.History.Ethnicity", "Smoking",
                 "YearsSinceDiagPso",  
                 "patient.Evidence.of.superinfection..Non.lesional.skin",
                 "patient.HEAD.HEAD...Total..Sum.x.Area.x0.1.","patient.TRUNK.TRUNK...Total..Sum.x.Area.x0.3.",
                 "patient.UPPER.LIMBS.UPPER.LIMBS...Total..Sum.x.Area.x0.2.","patient.LOWER.LIMBS.LOWER.LIMBS...Total..Sum.x.Area.x0.4.",
                 "patient.PASI.Score.Total..TotalHead...TotalTrunk...TotalUpperLimbs...TotalLowerLimbs..PASI.Score.",
                 "patient.Psoriasis.Severity.Index..PSI..PSI...Erythema","patient.Psoriasis.Severity.Index..PSI..PSI...Induration",
                 "patient.Psoriasis.Severity.Index..PSI..PSI...Scaling","patient.Psoriasis.Severity.Index..PSI..PSI.Total",
                 "patient.Please.state.the.areas.used.for.sampling.Lesional.skin", 
                 "patient.Clinical.Assessment.Global.Assessment.Score",
                 "patient.Hanifin.and.Rajka.diagnostic.criteria.6.Palmar.hyperlinearity..Keratosis.pilaris",
                 "patient.Hanifin.and.Rajka.diagnostic.criteria.7.1.Total.IgE..IU.ml.",
                 "patient.Hanifin.and.Rajka.diagnostic.criteria.7.1.a..150.IU.ml...x...400.IU.ml",
                 "patient.Hanifin.and.Rajka.diagnostic.criteria.7.1.b...400.IU.ml",
                 "patient.Hanifin.and.Rajka.diagnostic.criteria.8..Raised.serum.IgE",
                 "patient.Hanifin.and.Rajka.diagnostic.criteria.9..Early.age.of.onset",
                 "patient.Hanifin.and.Rajka.diagnostic.criteria.10..ECP.measurement..ug.ml.",
                 "patient.Hanifin.and.Rajka.diagnostic.criteria.11..Tendency.toward.cutaneous.infections.or.impaired.cell.mediated.immunity",
                 "patient.Hanifin.and.Rajka.diagnostic.criteria.16..Food.intolerance",
                 "patient.Hanifin.and.Rajka.diagnostic.criteria.33.a..Personal.history.of.",
                 "patient.33.b..Family.history.of.atopy.Fam..hist..Atopic.dermatitis",
                 "patient.SCORAD.Objective.SCORAD..A.5...7B.2.",
                 "patient.SCORAD.SCORAD.Score",
                 "patient.A..Please.state.the.areas.used.for.sampling.Lesional.skin",
                 "patient.B..SCORAD..Intensity..SCORAD..Erythema",
                 "patient.B..SCORAD..Intensity..SCORAD..Edema.Papulation",
                 "patient.B..SCORAD..Intensity..SCORAD..Oozing.crusts",
                 "patient.B..SCORAD..Intensity..SCORAD..Excoriation",
                 "patient.B..SCORAD..Intensity..SCORAD..Lichenification",
                 "patient.B..SCORAD..Intensity..SCORAD..Dryness..evaluated.on.uninvolved.areas.",
                 "patient.B..LOCAL..Intensity..LOCAL..Erythema",
                 "patient.B..LOCAL..Intensity..LOCAL..Edema.Papulation",                                                                                   
                 "patient.B..LOCAL..Intensity..LOCAL..Oozing.crusts",                                                                                      
                 "patient.B..LOCAL..Intensity..LOCAL..Excoriation",                                                                                        
                 "patient.B..LOCAL..Intensity..LOCAL..Lichenification",                                                                                    
                 "patient.B..LOCAL..Intensity..LOCAL..Dryness..evaluated.on.uninvolved.areas.",                                                            
                 "patient.B..LOCAL..Intensity.LOCAL.SCORAD.Score"
)

patient_data <- rbind(patient_data.ad[,cols_to_get],
                      patient_data.ctrl[,cols_to_get],
                      patient_data.pso[,cols_to_get])

#**********************************
# Integrate sequencing and clinical data
#**********************************

# Inner join both datasets
clin_seq_metadata <- merge(x=clinical_data, y=sequencing_data,  by.x="sample_id", by.y="MAARS_sample_id")

dim(clin_seq_metadata)
dim(clinical_data)
colnames(clinical_data)
colnames(sequencing_data)
colnames(clin_seq_metadata)
dim(sequencing_data)
head(sequencing_data)

# Excluded patients
# Logic vector of exclusions
patient_data$excluded <- patient_data$patient.Curated.Items.Excluded == "Yes"

#List of excluded patietns
excluded_patients <- patient_data$patient.Identification.MAARS.identifier[patient_data$excluded]

## integrate patient and sample level information ###

clin_seq_pt_metadata = merge(x=patient_data, y=clin_seq_metadata,
                 by.x="patient.Identification.MAARS.identifier",
                 by.y="MAARS_identifier")

dim(clin_seq_pt_metadata)
head(clin_seq_pt_metadata)
colnames (patient_data)
colnames(clin_seq_metadata)
colnames(clin_seq_pt_metadata)
# it's not gender, it's sex.
names(clin_seq_pt_metadata)[grep("Gender", names(clin_seq_pt_metadata))] = "Sex"

# just to simplify the name.
names(clin_seq_pt_metadata)[grep("CUSTOM_Age", names(clin_seq_pt_metadata))] = "Age"

names(clin_seq_pt_metadata)[grep("clinical_group.x",names(clin_seq_pt_metadata))] = "clinical_group"

### things that are not merged
setdiff(sequencing_data$MAARS_sample_id, clin_seq_metadata$sample_id)

#### OK
length(clin_seq_metadata$ScilifeID)
length(listOfScilifeID)
length(intersect(listOfScilifeID, clin_seq_metadata$ScilifeID))
samplesWithAbsentMetadata = setdiff(x=listOfScilifeID, y=intersect(clin_seq_metadata$ScilifeID, listOfScilifeID))

# the scilife IDs we have in sequencing data
length(sequencing_data$ScilifeID)

# the number of scilife IDs I have in the cluster
length(listOfScilifeID)

# the length of the intersection
length(intersect(listOfScilifeID, sequencing_data$ScilifeID))

# things for which I have a scilife ID (in the cluster), but are not in seq_data
samplesWithAbsentSequencingData = setdiff(x=listOfScilifeID, y=sequencing_data$ScilifeID)#

absentMetaAndNotNegControls_ScilifeIDs = setdiff(samplesWithAbsentMetadata, sequenced_negativeControls$ScilifeID)
absentMetaAndNotNegControls_MAARS_sampleID = sequencing_data$MAARS_sample_id[sequencing_data$ScilifeID %in% absentMetaAndNotNegControls_ScilifeIDs]

tmpFun = function(i){x = strsplit(x=i, split="_"); o=x[[1]]; o=paste(o[1], "_",o[2],"_", o[3], sep=""); return(o)}
absentMetaAndNotNegControls_MAARS_ID = sapply(X = as.character(absentMetaAndNotNegControls_MAARS_sampleID), FUN=tmpFun )

absentMetaAndNotNegControls_MAARS_ID[!absentMetaAndNotNegControls_MAARS_ID %in% excluded_patients]

#********************************************
# 1) Get putative duplicated sequenced samples
#*********************************************
# Q: why are there duplicated samples?
# {is it duplicated sequencing or mislables?}
duplicated_samples <- sequencing_data$MAARS_sample_id[duplicated(sequencing_data$MAARS_sample_id)]

#*******************************************************************
# Find Sequenced samples LACKING associated clinical data
# Likely samples from patients excluded from study at a later
# timepoint + putative typos
#*******************************************************************

no_clinical_data.maars <- setdiff( sequencing_data$MAARS_sample_id, clinical_data$sample_id ) ## where these sequenced?
no_clinical_data.seq <- lapply(X = no_clinical_data.maars,FUN = function(x) sequencing_data[sequencing_data$MAARS_sample_id == x,c("RunID","ScilifeID") ] )

no_clinical_data <- data.frame( RunID=as.character(sapply(no_clinical_data.seq,function(x) as.character(x$RunID) )) ,
                                ScilifeID=as.character(sapply(no_clinical_data.seq,function(x) x$ScilifeID )) ,
                                sample_id=no_clinical_data.maars
)

#Get patient id
no_clinical_data$patient_id <- sapply(no_clinical_data$sample_id,function(x)substr(x,1,11))

#Check if patient has been excluded
intersect(no_clinical_data$patient_id,excluded_patients)
no_clinical_data[ no_clinical_data$patient_id %in% intersect(no_clinical_data$patient_id,excluded_patients),]

#Samples from non excluded patients
setdiff(no_clinical_data$patient_id,excluded_patients)
no_clinical_data[ no_clinical_data$patient_id %in% setdiff(no_clinical_data$patient_id,excluded_patients),]

#List samples with no associated patient data in patient db
setdiff(setdiff(no_clinical_data$patient_id,excluded_patients), patient_data$patient.Identification.MAARS.identifier)

#*******************************************************************
# Find clinical samples with no sequenced samples
# Likely samples that had little DNA for sequencing or missing
# tubes
#*******************************************************************
not_sequenced <- setdiff(  clinical_data$sample_id,sequencing_data$MAARS_sample_id)
#Exclude samples wtih _03 and _04 that were not included in the study
not_sequenced.filt <- not_sequenced[! grepl("_0[34]$",not_sequenced)]
not_sequenced.filt


#**********************************************************************
# Add negative controls (rows) to clin_seq_pt_metadata 
#**********************************************************************
sequenced_negativeControls <- read.csv("/Users/helenalexander/Documents/MAARS/data/metaphlan/analysis/metadata/sequenced_ntcs.csv", sep=";", stringsAsFactors=FALSE)

## Creat NA columns for all columns in clin_seq_pt_metadata not present in sequenced_negativeControls
colnames(clin_seq_pt_metadata)
sequenced_negativeControls$"clinical_group" = NA
sequenced_negativeControls$"patient.Curated.Items.Excluded" = NA
sequenced_negativeControls$"patient.Curated.Items.Exclusion.comments" = NA
sequenced_negativeControls$"patient.Ethnicity.Family.History.Ethnicity" = NA
sequenced_negativeControls$"Smoking" = NA
sequenced_negativeControls$"YearsSinceDiagPso" = NA
sequenced_negativeControls$"patient.Evidence.of.superinfection..Non.lesional.skin" = NA
sequenced_negativeControls$"patient.HEAD.HEAD...Total..Sum.x.Area.x0.1." = NA
sequenced_negativeControls$"patient.TRUNK.TRUNK...Total..Sum.x.Area.x0.3." = NA
sequenced_negativeControls$"patient.UPPER.LIMBS.UPPER.LIMBS...Total..Sum.x.Area.x0.2." = NA
sequenced_negativeControls$"patient.LOWER.LIMBS.LOWER.LIMBS...Total..Sum.x.Area.x0.4." = NA
sequenced_negativeControls$"patient.PASI.Score.Total..TotalHead...TotalTrunk...TotalUpperLimbs...TotalLowerLimbs..PASI.Score." = NA
sequenced_negativeControls$"patient.Psoriasis.Severity.Index..PSI..PSI...Erythema" = NA
sequenced_negativeControls$"patient.Psoriasis.Severity.Index..PSI..PSI...Induration" = NA
sequenced_negativeControls$"patient.Psoriasis.Severity.Index..PSI..PSI...Scaling" = NA
sequenced_negativeControls$"patient.Psoriasis.Severity.Index..PSI..PSI.Total" = NA
sequenced_negativeControls$"patient.Please.state.the.areas.used.for.sampling.Lesional.skin" = NA
sequenced_negativeControls$"patient.Clinical.Assessment.Global.Assessment.Score" = NA
sequenced_negativeControls$"patient.Hanifin.and.Rajka.diagnostic.criteria.6.Palmar.hyperlinearity..Keratosis.pilaris" = NA                                
sequenced_negativeControls$"patient.Hanifin.and.Rajka.diagnostic.criteria.7.1.Total.IgE..IU.ml."= NA                            
sequenced_negativeControls$"patient.Hanifin.and.Rajka.diagnostic.criteria.7.1.a..150.IU.ml...x...400.IU.ml" = NA                        
sequenced_negativeControls$"patient.Hanifin.and.Rajka.diagnostic.criteria.7.1.b...400.IU.ml" = NA                    
sequenced_negativeControls$"patient.Hanifin.and.Rajka.diagnostic.criteria.8..Raised.serum.IgE"= NA                                                        
sequenced_negativeControls$"patient.Hanifin.and.Rajka.diagnostic.criteria.9..Early.age.of.onset"= NA                                                      
sequenced_negativeControls$"patient.Hanifin.and.Rajka.diagnostic.criteria.10..ECP.measurement..ug.ml." = NA                                               
sequenced_negativeControls$"patient.Hanifin.and.Rajka.diagnostic.criteria.11..Tendency.toward.cutaneous.infections.or.impaired.cell.mediated.immunity"= NA
sequenced_negativeControls$"patient.Hanifin.and.Rajka.diagnostic.criteria.16..Food.intolerance"= NA                                                       
sequenced_negativeControls$"patient.Hanifin.and.Rajka.diagnostic.criteria.33.a..Personal.history.of."= NA                                                 
sequenced_negativeControls$"patient.33.b..Family.history.of.atopy.Fam..hist..Atopic.dermatitis"= NA                                                       
sequenced_negativeControls$"patient.SCORAD.Objective.SCORAD..A.5...7B.2."= NA                                                                             
sequenced_negativeControls$"patient.SCORAD.SCORAD.Score"= NA                                                                                              
sequenced_negativeControls$"patient.A..Please.state.the.areas.used.for.sampling.Lesional.skin"= NA                                                        
sequenced_negativeControls$"patient.B..SCORAD..Intensity..SCORAD..Erythema"= NA                                                                           
sequenced_negativeControls$"patient.B..SCORAD..Intensity..SCORAD..Edema.Papulation" = NA                                                                  
sequenced_negativeControls$"patient.B..SCORAD..Intensity..SCORAD..Oozing.crusts" = NA                                                                     
sequenced_negativeControls$"patient.B..SCORAD..Intensity..SCORAD..Excoriation"= NA                                                                        
sequenced_negativeControls$"patient.B..SCORAD..Intensity..SCORAD..Lichenification" = NA                                                                   
sequenced_negativeControls$"patient.B..SCORAD..Intensity..SCORAD..Dryness..evaluated.on.uninvolved.areas."= NA                                            
sequenced_negativeControls$"patient.B..LOCAL..Intensity..LOCAL..Erythema" = NA                                                                            
sequenced_negativeControls$"patient.B..LOCAL..Intensity..LOCAL..Edema.Papulation" = NA                                                                    
sequenced_negativeControls$"patient.B..LOCAL..Intensity..LOCAL..Oozing.crusts" = NA                                                                       
sequenced_negativeControls$"patient.B..LOCAL..Intensity..LOCAL..Excoriation" = NA                                                                         
sequenced_negativeControls$"patient.B..LOCAL..Intensity..LOCAL..Lichenification"= NA                                                                      
sequenced_negativeControls$"patient.B..LOCAL..Intensity..LOCAL..Dryness..evaluated.on.uninvolved.areas."= NA                                              
sequenced_negativeControls$"patient.B..LOCAL..Intensity.LOCAL.SCORAD.Score" = NA                                                                          
sequenced_negativeControls$"excluded"= NA                                                                                                                 
sequenced_negativeControls$"sample_id" = NA                                                                                                               
sequenced_negativeControls$"clinical_group.y" = NA                                                                                                        
sequenced_negativeControls$"anatomical_location" = NA                                                                                                     
sequenced_negativeControls$"anatomical_location_label" = NA                                                                                               
sequenced_negativeControls$"lesional" = NA                                                                                                                
sequenced_negativeControls$"Institution"= NA                                                                                                              
sequenced_negativeControls$"Age" = NA                                                                                                                     
sequenced_negativeControls$"Sex" = NA                                                                                                                     
sequenced_negativeControls$"Known_Allergies_v2..Pseudo_Drug_Allergy"= NA                                                                                  
sequenced_negativeControls$"Known_Allergies_v2..House_dust_mite"  = NA                                                                                    
sequenced_negativeControls$"Known_Allergies_v2..Food"  = NA                                                                                               
sequenced_negativeControls$"Known_Allergies_v2..Pollen"  = NA                                                                                             
sequenced_negativeControls$"Known_Allergies_v2..Contact_Allergy" = NA                                                                                     
sequenced_negativeControls$"Known_Allergies_v2..Drug_Allergy"  = NA                                                                                       
sequenced_negativeControls$"Known_Allergies_v2..Animal"   = NA                                                                                            
sequenced_negativeControls$"Concomitant_Medication_v2..Anti.Hypertensive"    = NA                                                                         
sequenced_negativeControls$"Concomitant_Medication_v2..Anti.Inflammatory.non_steroid."   = NA                                                             
sequenced_negativeControls$"Concomitant_Medication_v2..Other_hormones" = NA                                                                               
sequenced_negativeControls$"Concomitant_Medication_v2..Thyroid_hormones" = NA                                                                             
sequenced_negativeControls$"Concomitant_Medication_v2..Statins"   = NA                                                                                    
sequenced_negativeControls$"Concomitant_Medication_v2..Insulin"   = NA                                                                                    
sequenced_negativeControls$"Concomitant_Medication_v2..Others"    = NA                                                                                    
sequenced_negativeControls$"Other_concurrent_chronic_diseases_v2..Hyperlipidemia"   = NA                                                                  
sequenced_negativeControls$"Other_concurrent_chronic_diseases_v2..Others"           = NA                                                                  
sequenced_negativeControls$"Other_concurrent_chronic_diseases_v2..Diabetes_.non.insulin."   = NA                                                          
sequenced_negativeControls$"Other_concurrent_chronic_diseases_v2..Thyroid_dysfunction"= NA                                                                
sequenced_negativeControls$"Other_concurrent_chronic_diseases_v2..Asthma"   = NA                                                                          
sequenced_negativeControls$"Other_concurrent_chronic_diseases_v2..Hypertension"  = NA                                                                     
sequenced_negativeControls$"Global_Assessment_Score"    = NA                                                                                              
sequenced_negativeControls$"CUSTOM_Malignancies_._skin" = NA                                                                                              
sequenced_negativeControls$"CUSTOM_Malignancies_._other"  = NA                                                                                            
sequenced_negativeControls$"CUSTOM_Fam._hist._Atopic_dermatitis"    = NA                                                                                  
sequenced_negativeControls$"CUSTOM_Family_History_of_Psoriasis"      = NA  
sequenced_negativeControls$"patient.Identification.MAARS.identifier" = NA

## Creat NA columns for all columns in sequenced_negativeControls not present in clin_seq_pt_metadata
clin_seq_pt_metadata$"MOCK_ID" = NA

# make SubmittedID in sequenced_negativeControls = MOCK_ID
colnames(sequenced_negativeControls)[3] = "MOCK_ID"

#check have same number of columns
dim(clin_seq_pt_metadata)
dim(sequenced_negativeControls)

setdiff (colnames(sequenced_negativeControls), colnames(clin_seq_pt_metadata))

# makes order of the columns the same for both DFs
Os = order(colnames(sequenced_negativeControls))
Om = order(colnames(clin_seq_pt_metadata))
clin_seq_pt_metadata = clin_seq_pt_metadata [,Om]
sequenced_negativeControls = sequenced_negativeControls [,Os]                     

# rbind both DFs
dim(clin_seq_pt_metadata)
clin_seq_pt_metadata = rbind(clin_seq_pt_metadata, sequenced_negativeControls)
dim(clin_seq_pt_metadata)

### Add column on number of reads left in each sample
load(file="/Users/helenalexander/Documents/MAARS/data/metaphlan/analysis/metadata/humanremovalSummary.RData")
humanremovalSummary_df$pairsLeft = humanremovalSummary_df$readPairs - humanremovalSummary_df$pairsMapped

colnames(clin_seq_pt_metadata)[82] = "scilifeID"
colnames(humanremovalSummary_df)
colnames(clin_seq_pt_metadata)

clin_seq_pt_metadata = merge(x=clin_seq_pt_metadata, y=humanremovalSummary_df[,c("scilifeID", "pairsLeft")], 
                 by.x="scilifeID", by.y="scilifeID")
dim(clin_seq_pt_metadata)

## Remove clinical_group.y and change names of clinical_group.x to clinical_group
wc_gy = grep ("clinical_group.y", colnames(clin_seq_pt_metadata))
wc_gy
clin_seq_pt_metadata = clin_seq_pt_metadata[,-wc_gy]

## Remove anatomical location
colnames(clin_seq_pt_metadata)
clin_seq_pt_metadata = clin_seq_pt_metadata[,-3]

clin_seq_pt_metadata[,3]

# Make T and PT both = T 
wPT = grep ("PT", clin_seq_pt_metadata[,"anatomical_location_label"])
class (wPT)
clin_seq_pt_metadata [wPT, "anatomical_location_label"] = "T"

### Unify mock samples 
# mock_1.9, _1.15 = PBS HHU
# _2.7, _2.1 = PBS KINGS
# _3.2, _3.5 = PBS UH
# _4.15, _4.29, 4.33 = PBS KI (Scilife)
# H20 = library prep at KI (RunIDs 13_08, 14_05, 14_06, 14_07, 14_08, 15_01, 15_02, 15_04)

w_1 = grep ("_1", clin_seq_pt_metadata[,"MOCK_ID"])
clin_seq_pt_metadata [w_1, "MOCK_ID"] = "Mock_1"
w_2 = grep ("_2", clin_seq_pt_metadata[,"MOCK_ID"])
clin_seq_pt_metadata [w_2, "MOCK_ID"] = "Mock_2"
w_3 = grep ("_3", clin_seq_pt_metadata[,"MOCK_ID"])
clin_seq_pt_metadata [w_3, "MOCK_ID"] = "Mock_3"
w_4 = grep ("_4", clin_seq_pt_metadata[,"MOCK_ID"])
clin_seq_pt_metadata [w_4, "MOCK_ID"] = "Mock_4"

#*******************************************************************
# Save clin_seq_metadata to object
#*******************************************************************
save(list=c("patient_data","clinical_data", "clin_seq_metadata","no_clinical_data","not_sequenced.filt"), 
     file = "/Users/helenalexander/Documents/MAARS/data/metaphlan/analysis/metadata/20160301clin_seq_metadata.RData")

#****************************************
# Save clin_seq_pt_metadata
#****************************************

save (list = c("clin_seq_pt_metadata"), file = "/Users/helenalexander/Documents/MAARS/data/metaphlan/analysis/metadata/20160301clin_seq_pt_metadata.Rdata")
write.csv(x = clin_seq_pt_metadata, file = "/Users/helenalexander/Documents/MAARS/data/metaphlan/analysis/metadata/20160301clin_seq_pt_metadata.csv", quote=FALSE, row.names
          = FALSE)

