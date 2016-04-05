#All sequenced samples should be in the metadata dataframe
stopifnot(nrow(metadata) == nrow(sequencing_data) )

#Unannotated samples -> non-mock samples with an NA in clinical_group
unannot.sample.count <- metadata %>% filter( is.na(clinical_group) & !is_ntc ) %>% select(scilife_id,maars_sample_id) %>% nrow
ntc.sample.count <- metadata %>% filter(is_ntc) %>% nrow

#Validate completness of variables (sequencing data + sample metadata)

#This variables should not have NAs
for( col in c("maars_sample_id","batch_id","scilife_id","million_reads_sequenced",
              "barcodeseq","filtered_reads","is_ntc","excluded_from_sample_metadata","resequenced") ){
  print(paste("Testing column",col,"for NAs"))
  stopifnot( (metadata[col] %>% is.na %>% sum) == 0 )
  print("OK")
}

#This columns should have NAs only for the unannotated rows and the mocks
for( col in c("lesional","sex","age","anatomical_location") ){
  print(paste("Testing column",col,"for completeness"))
  stopifnot( (metadata[col] %>% is.na %>% sum) == (unannot.sample.count + ntc.sample.count) )
  print("OK")
}

#maars_subject_id should be set for all samples except mocks
print("Testing column maars_subject_id for completeness")
stopifnot( (metadata$maars_subject_id %>% is.na %>% sum ) == ntc.sample.count )
print("OK")

#Institute should be set for all samples except mocks
print("Testing column Institution for completeness")
stopifnot( (metadata$institution %>% is.na %>% sum ) == unannot.sample.count )
print("OK")

#Verify that all duplicated sample_ids are appropriately annotated as resequenced
print("Check all duplicated samples are correctly annotated as such")
duplicated_maars_sample_id <- (metadata %>% filter(!is_ntc) %>% filter(duplicated(maars_sample_id)) %>% select(maars_sample_id))$maars_sample_id
unannot.dups <- metadata %>% filter(maars_sample_id %in% duplicated_maars_sample_id) %>% filter(!resequenced ) %>%
  select(batch_id,maars_sample_id,million_reads_sequenced,resequenced,resequencing_notes,excluded_from_sample_metadata) %>% arrange(maars_sample_id,batch_id)
if(nrow(unannot.dups) == 0){
  print("OK")
}else{
  warning(paste(unique(unannot.dups$maars_sample_id),"not annotated as resequenced\n"))
  print(unannot.dups)
}
rm(unannot.dups,duplicated_maars_sample_id)



#Report all other sample number inconsistencies
#List volunteers with more than 2 samples []
print("Check if any volunteer has more than 2 distinct samples [duplicates were addressed in the previous test]")
volunteers.too_many_samples <- (metadata %>% filter(!is_ntc & ! resequenced) %>% distinct(maars_sample_id) %>%
                                            count(maars_subject_id) %>% 
                                            filter( n > 2))$maars_subject_id
if(length(volunteers.too_many_samples) > 0){
  warning("Some volunteers have more than 2 distinct samples associated")
  metadata %>% filter(maars_subject_id %in% volunteers.too_many_samples) %>% 
    select(batch_id,scilife_id,maars_subject_id,maars_sample_id) %>% 
    arrange(maars_subject_id) %>% print
}else {
  print("OK")
}
rm(volunteers.too_many_samples)

#Cleanup
rm(ntc.sample.count, unannot.sample.count,col)