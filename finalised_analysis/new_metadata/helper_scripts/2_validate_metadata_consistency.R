#All sequenced samples should be in the metadata dataframe
stopifnot(nrow(metadata) == nrow(sequencing_data) )

#Unannotated samples -> non-mock samples with an NA in clinical_group
unannot.sample.count <- metadata %>% filter( is.na(clinical_group) & !is_ntc ) %>% select(scilife_id,maars_sample_id) %>% nrow
ntc.sample.count <- metadata %>% filter(is_ntc) %>% nrow

#Validate completness of variables (sequencing data + sample metadata)

#This variables should not have NAs
for( col in c("maars_sample_id","batch_id","scilife_id","million_reads_sequenced","barcodeseq","filtered_reads","is_ntc","excluded_from_sample_metadata") ){
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

#maars_patient_id should be set for all samples except mocks
print("Testing column maars_patient_id for completeness")
stopifnot( (metadata$maars_patient_id %>% is.na %>% sum ) == ntc.sample.count )
print("OK")

#Institute should be set for all samples except mocks
print("Testing column Institution for completeness")
stopifnot( (metadata$institution %>% is.na %>% sum ) == unannot.sample.count )
print("OK")


#Cleanup
rm(ntc.sample.count, unannot.sample.count,col)
