megadata.nonexcluded <- megadata %>% filter(is_ntc | (in_sample_metadata & !excluded ))

#Non-negative controls & non-excluded samples annotated in the first step should also have patient level annotation
print("Check if non-negative controls have patient_metadata annotation")
stopifnot( (megadata.nonexcluded %>% filter(!is_ntc) %>% select(patient_clinical_group) %>% is.na %>% sum) == 0 )
print("OK")

#****************************************************
# Consistency between sample and patient metadata
#****************************************************
#Validate clinical group in sample metadata corresponds to clinical group in patient metadata
print("Check whether patient metadata clinical group matches sample metadata clinical group")
stopifnot ( (megadata.nonexcluded %>% filter(!is_ntc) %>% select(clinical_group,patient_clinical_group) %>% 
  mutate(diff= clinical_group != patient_clinical_group))$diff %>% sum == 0)
print("OK")

print("Check whether sample and patient metadata global assessment scores match")
stopifnot( (factor(megadata.nonexcluded$clinical.assessment.global.assessment.score,
       labels=levels(megadata.nonexcluded$global_assessment_score)) != megadata.nonexcluded$global_assessment_score ) %>% sum(na.rm=T) == 0)
stopifnot( sum( is.na(megadata.nonexcluded$global_assessment_score) != is.na(megadata.nonexcluded$clinical.assessment.global.assessment.score) ) == 0 )
print("OK")

print("Check whether anatomical location for lesional samples matches sampled_lesional_skinsite")
tmp_factor <- factor((filter(megadata,lesional=="LES"))$sampled_lesional_skinsite)
#Replace [ "buttocks",  "lower back" ,"posterior thigh", "upper back" ] with these for comparison
levels(tmp_factor) <- c("lower_back","lower_back","thigh","upper_back")
stopifnot( sum( filter(megadata,lesional=="LES")$anatomical_location != tmp_factor ) == 0 ) 
rm(tmp_factor)
print("OK")


print("Are there samples not annotated in the sample metadata nor in the patient metadata")
tmp.df <-megadata %>% filter(!is_ntc) %>% filter(!in_sample_metadata &  !in_patient_metadata) %>% select(scilife_id,maars_sample_id)
if(nrow(tmp.df) > 0){
  warning("There are samples with no known metadata")
  print(tmp.df)
}else{
  print("OK. No unannotated samples")
}


print("List of samples not present in the sample metadata which are NOT marked as excluded in patient metadata")
# Most samples with not sample metadata are marked as excluded volunteers in patient metadata, 
# except for 4 samples from three patients (1,1,2) which were omitted from the sample_metadata table, 
# but have no indication of having an excluded status
megadata %>% filter(!is_ntc) %>% 
             filter(!in_sample_metadata & in_patient_metadata ) %>% 
             filter(!excluded) %>% select(scilife_id,maars_sample_id,samples_per_subject,patient_clinical_group) %>% print


#*******************************************************************************
#Check completeness of patient metadata columns
#*******************************************************************************
ntc.count <- megadata.nonexcluded %>% filter(is_ntc) %>% nrow
ad.count <-  megadata.nonexcluded %>% filter(clinical_group == "AD") %>% nrow 
pso.count <-  megadata.nonexcluded %>% filter(clinical_group == "PSO") %>% nrow 
ctrl.count <-  megadata.nonexcluded %>% filter(clinical_group == "CTRL") %>% nrow 

#For all volunteers
for( col in c("ethnicity","smoking") ){
  print(paste("Testing column",col,"for incomplete values"))
  stopifnot( (megadata.nonexcluded[col] %>% is.na %>% sum) == ntc.count )
  print("OK")
}

#Should be in all AD and PSO patients
for( col in c("global_assessment_score") ){
  print(paste("Testing column",col,"for incomplete values"))
  stopifnot( (megadata.nonexcluded[col] %>% is.na %>% sum) == (ntc.count+ctrl.count))
  print("OK")
}

for( col in grep("(pasi)|(psoriasis)|(superinfection)",colnames(megadata),value=TRUE) ){
  print(paste("Testing column",col,"for incomplete values"))
  if((megadata.nonexcluded[col] %>% is.na %>% sum) == (ntc.count+ctrl.count+ad.count) ){
    print("OK")    
  }else{
    print("FAIL")
    warning( paste(col," is not complete") ) 
  }
}


for( col in grep("(hanifin)|(scorad)|(LOCAL)|(SCORAD)",colnames(megadata),value=TRUE)){
  print(paste("Testing column",col,"for incomplete values"))
  if((megadata.nonexcluded[col] %>% is.na %>% sum) == (ntc.count+ctrl.count+pso.count) ){
    print("OK")    
  }else{
    print("FAIL")
    warning( paste(col," is not complete") ) 
  }
}


rm(megadata.nonexcluded)
rm(ntc.count,ad.count,pso.count,ctrl.count)

#*******************************************************************************
#Verify list of samples in the Rosalind cluster corresponds to the samples in megadata
#*******************************************************************************
rosalind.samples = read.delim2("input_files/rosalind_sample_list.txt", header=FALSE)$V1
rosalind.samples.scilife_id <- gsub(".+/","",rosalind.samples,perl=TRUE,fixed=FALSE)
rosalind.samples.scilife_id <- gsub("_index.+","",rosalind.samples.scilife_id,perl=TRUE,fixed=FALSE)

print("Verify that all samples in the Rosalind cluster are the same as in the megadata")
if( length(intersect(rosalind.samples.scilife_id, megadata$scilife_id)) ==  nrow(megadata) ){
  print("OK")
}else{
  #The cluster does not have the raw data of the last HiSeq run, but the filtered data is there
  print( "Samples in the megadata not present in Rosalind")
  print( setdiff(megadata$scilife_id,rosalind.samples.scilife_id))
  print( "Samples in the megadata not present in Rosalind")
  print( setdiff(rosalind.samples.scilife_id,megadata$scilife_id))
  warning("There are inconsistencies between samples in Rosalind and megadata")
}
rm(rosalind.samples,rosalind.samples.scilife_id)
