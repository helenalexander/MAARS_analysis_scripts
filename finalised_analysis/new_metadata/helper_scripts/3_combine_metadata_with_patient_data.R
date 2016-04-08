stopifnot("metadata" %in% ls())

#***************************************************************
#3.1 Load patient metadata
#***************************************************************
patient_data.ctrl <- read.table("input_files/MAARS_Control_full_20150424_15-38-30-4.csv",
                                sep="\t",header=TRUE,quote="\"",check.names=TRUE,comment.char="",stringsAsFactors=FALSE,na.strings=c("NA",""))
patient_data.ctrl$patient_clinical_group <- "CTRL"

patient_data.ad <- read.table("input_files/MAARS_AD_full_20150424_15-37-35-3.csv",
                              sep="\t",header=TRUE,quote="\"",check.names=TRUE,comment.char="",stringsAsFactors=FALSE,,na.strings=c("NA",""))
patient_data.ad$patient_clinical_group <- "AD"

patient_data.pso <- read.table("input_files/MAARS_PSO_full_20150423_18-43-57-3.csv",
                               sep="\t",header=TRUE,quote="\"",check.names=TRUE,comment.char="",stringsAsFactors=FALSE,,na.strings=c("NA",""))
patient_data.pso$patient_clinical_group <- "PSO"


colnames(patient_data.ctrl) %>% length
colnames(patient_data.ad) %>% length
colnames(patient_data.pso) %>% length

intersect(colnames(patient_data.ctrl), colnames(patient_data.ad)   ) %>% length
intersect(colnames(patient_data.ctrl), colnames(patient_data.pso)   ) %>% length
intersect(colnames(patient_data.pso), colnames(patient_data.ad)   ) %>% length

#***************************************************************
# 3.2 Data wrangling of patient data at each clinical group
#***************************************************************

#Make variable names equal throughout patient metadata files for Smoking
colnames(patient_data.pso)[18] = "smoking"
#TODO!!! 
#Not sure "patient.Ethnicity.Family.History.Smoking" means the person itself smokes
colnames(patient_data.ctrl)[11] = "smoking"
colnames(patient_data.ad)[10] = "smoking"


#Join sampled_lesional_skinsite variables into a single column
colnames(patient_data.pso)[grep("patient.Please.state.the.areas.used.for.sampling.Lesional.skin", colnames(patient_data.pso))] <- "sampled_lesional_skinsite"
colnames(patient_data.ad)[grep("Please.state.the.areas.used.for.sampling.Lesional.skin", colnames(patient_data.ad))] <- "sampled_lesional_skinsite"

# Create variable for "YearsSinceDiagnosis" for Psoriasis = YearsSinceDiagPso 
# (= "patient.Identification.Date.of.first.visit" - "patient.Psoriasis.Patients.Year.of.Diagnosis")
#Extract year from date of first visit
IntYearFV <- as.Date(patient_data.pso[ ,2],format = "%d/%m/%y") %>% format("%Y") %>% as.numeric()
#Calculate YsinceD
YearsSinceDiagPso <- IntYearFV - patient_data.pso[["patient.Psoriasis.Patients.Year.of.Diagnosis"]]
patient_data.pso[,"psoriasis_years_since_diagnosis"] = YearsSinceDiagPso

#*********************************************************************
# 3.3 Select clinical variables of interest and merge all patient data
#       for the three groups into a single dataframe
#*********************************************************************
cols_to_get <- c("patient.Identification.MAARS.identifier", 
                 "patient_clinical_group",
                 "patient.Curated.Items.Excluded",
                 "patient.Curated.Items.Exclusion.comments",
                 "patient.Ethnicity.Family.History.Ethnicity", 
                 "smoking",
                 "psoriasis_years_since_diagnosis", 
                 "sampled_lesional_skinsite",
                 "patient.Evidence.of.superinfection..Non.lesional.skin",
                 "patient.HEAD.HEAD...Total..Sum.x.Area.x0.1.","patient.TRUNK.TRUNK...Total..Sum.x.Area.x0.3.",
                 "patient.UPPER.LIMBS.UPPER.LIMBS...Total..Sum.x.Area.x0.2.","patient.LOWER.LIMBS.LOWER.LIMBS...Total..Sum.x.Area.x0.4.",
                 "patient.PASI.Score.Total..TotalHead...TotalTrunk...TotalUpperLimbs...TotalLowerLimbs..PASI.Score.",
                 "patient.Psoriasis.Severity.Index..PSI..PSI...Erythema",
                 "patient.Psoriasis.Severity.Index..PSI..PSI...Induration",
                 "patient.Psoriasis.Severity.Index..PSI..PSI...Scaling",
                 "patient.Psoriasis.Severity.Index..PSI..PSI.Total",
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

#For rbind to work, all patient data dataframes should have all columns of interest
# Solution: Create dummy columns with NA's for the missing cols in each df
#Fill CTRL df
for(feature in setdiff(cols_to_get,colnames(patient_data.ctrl))){
  patient_data.ctrl[feature] <- NA
}
#Fill PSO df
for(feature in setdiff(cols_to_get,colnames(patient_data.pso))){
  patient_data.pso[feature] <- NA
}
#Fill AD df
for(feature in setdiff(cols_to_get,colnames(patient_data.ad))){
  patient_data.ad[feature] <- NA
}

#Merge the three dataframes into one
patient_data <- rbind(patient_data.ad[,cols_to_get],
                      patient_data.ctrl[,cols_to_get],
                      patient_data.pso[,cols_to_get])



rm(feature, IntYearFV,YearsSinceDiagPso)
rm(patient_data.ad,patient_data.pso,patient_data.ctrl)

#*******************************************************
# 3.4 Data wrangling of combined patient data
#*******************************************************

# Boolean vector of exclusions
patient_data$excluded <- ! is.na(patient_data$patient.Curated.Items.Excluded )
#Remove the Yes/no column in favor of the boolean one
patient_data <- subset(patient_data,select=-patient.Curated.Items.Excluded)

#********Rename variables directly****************
colnames(patient_data)[grep("patient.Curated.Items.Exclusion.comments", colnames(patient_data))] <- "exclusion_comments"
colnames(patient_data)[grep("patient.Ethnicity.Family.History.Ethnicity", colnames(patient_data))] <- "ethnicity"
colnames(patient_data)[grep("patient.PASI.Score.Total..TotalHead...TotalTrunk...TotalUpperLimbs...TotalLowerLimbs..PASI.Score.", colnames(patient_data))] <- "pasi_score"
colnames(patient_data)[grep("33.b..Family.history.of.atopy.Fam..hist..Atopic.dermatitis", colnames(patient_data))] <- "ad_family_history"
colnames(patient_data)[grep("Evidence.of.superinfection..Non.lesional.skin", colnames(patient_data))] <- "evidence_superinfection_nonles_skin"
colnames(patient_data)[grep("patient.Hanifin.and.Rajka.diagnostic.criteria.8..Raised.serum.IgE", colnames(patient_data))] <- "hanifin_and_rajka_8_raised_serum_ig_e"


#********Simplify variable names****************
#Remove patient.B prefix
colnames(patient_data) <- gsub("^patient\\.B\\.\\.","",colnames(patient_data),fixed=FALSE,perl=TRUE)
#Remove patient. prefix
colnames(patient_data) <- gsub("^patient\\.","",colnames(patient_data),fixed=FALSE,perl=TRUE)
colnames(patient_data) <- gsub("Psoriasis.Severity.Index..PSI..PSI...","psoriasis_severity_index_",colnames(patient_data),fixed=TRUE,perl=FALSE)
colnames(patient_data) <- gsub("SCORAD..Intensity..SCORAD..","scorad_intensity_",colnames(patient_data),fixed=TRUE,perl=FALSE)
colnames(patient_data) <- gsub("LOCAL\\.\\.Intensity\\.\\.?LOCAL\\.\\.?","local_intensity_",colnames(patient_data),fixed=FALSE,perl=TRUE)
colnames(patient_data) <- gsub("Hanifin.and.Rajka.diagnostic.criteria.","hanifin_and_rajka_",colnames(patient_data),fixed=TRUE,perl=FALSE)

#***********Turn all column names lowercase*************
colnames(patient_data) <- tolower(colnames(patient_data))

#*******************************************************
# 4.1 Join sequencing and patient data
#*******************************************************
## integrate patient and sample level information ###
megadata = merge(x=metadata, y=patient_data,
                 by.x="maars_subject_id",
                 by.y="identification.maars.identifier",
                 all.x=TRUE
)

#*******************************************************
# 4.2 Convert columns into factors / boolean columns appropriately
#*******************************************************
megadata$in_patient_metadata <- ! is.na(megadata$patient_clinical_group)

for(col in colnames(megadata)){
  if( class(megadata[[col]]) == "character" ){
    tmpfactor <- factor(megadata[[col]] )
#If the column has only values Yes/No convert to Boolean True / False
    if( length(levels(tmpfactor)) == 2 && sum(levels(tmpfactor) == c("No","Yes"))   ==  2 ){
      megadata[[col]] <- megadata[[col]] == "Yes"
    
#If the column is a scale (with values 1 (description) , 2, 3 etc ), convert into an ordered factor  
    }else if( sum(grepl("^[0-9].*",levels(tmpfactor))) == length(levels(tmpfactor))  ){
      megadata[[col]] <- factor(megadata[[col]],ordered = TRUE)
    }
  }
}

rm(col,tmpfactor,cols_to_get)

megadata["ethnicity"] <- factor(megadata[["ethnicity"]] )
megadata["smoking"] <- factor(megadata[["smoking"]] )
megadata["patient_clinical_group"] <- factor(megadata[["patient_clinical_group"]] )
megadata["sampled_lesional_skinsite"] <- factor(megadata[["sampled_lesional_skinsite"]] )
megadata["evidence_superinfection_nonles_skin"] <- factor(megadata[["evidence_superinfection_nonles_skin"]] )