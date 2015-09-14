---
title: "TCGA PRAD"
author: "Levi Waldron"
date: "September 14, 2015"
output: html_document
---

# Load TCGA data


```r
if(!require(RTCGAToolbox)){
  biocLite("devtools")
  biocLite(c("limma", "RCircos", "data.table", "RCurl", "RJSONIO"))
  biocLite("Link-NY/RTCGAToolbox")
}
```

# Download TCGA PRAD data


```r
library(RTCGAToolbox)
rundates <- getFirehoseRunningDates()
analysisdates <- getFirehoseAnalyzeDates()
prad <- getFirehoseData("PRAD", runDate=rundates[1],
                      gistic2_Date=analysisdates[1], RNAseq_Gene=TRUE, 
        miRNASeq_Gene=TRUE, RNAseq2_Gene_Norm=TRUE, CNA_SNP=TRUE,
        CNV_SNP=TRUE, CNA_Seq=TRUE, CNA_CGH=TRUE,  Methylation=TRUE,
        Mutation=TRUE, mRNA_Array=TRUE, miRNA_Array=TRUE, RPPA=TRUE)
```

```
## gdac.broadinstitute.org_PRAD.Clinical_Pick_Tier1.Level_4.2015040200.0.0.tar.gz
## Using locally cached version of ./20150402-PRAD-Clinical.txt
## 
## gdac.broadinstitute.org_PRAD.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2015040200.0.0.tar.gz
## Using locally cached version of ./20150402-PRAD-RNAseq2GeneNorm.txt
## RNAseq2 data will be imported! This may take a while!
## Start: 2015-09-14 17:45:19
```

```
## Read 48.7% of 20533 rowsRead 97.4% of 20533 rowsRead 20533 rows and 551 (of 551) columns from 0.090 GB file in 00:00:05
```

```
## Done: 2015-09-14 17:45:23
## gdac.broadinstitute.org_PRAD.Merge_mirnaseq__illuminahiseq_mirnaseq__bcgsc_ca__Level_3__miR_gene_expression__data.Level_3.2015040200.0.0.tar.gz
## Using locally cached version of ./20150402-PRAD-miRNAseqGene.txt
## miRNAseq data will be imported! This may take a while!
## Start: 2015-09-14 17:45:36
## Done: 2015-09-14 17:45:36
## gdac.broadinstitute.org_PRAD.Merge_snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_hg19__seg.Level_3.2015040200.0.0.tar.gz
## Using locally cached version of ./20150402-PRAD-CNASNPHg19.txt
## gdac.broadinstitute.org_PRAD.Merge_snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.Level_3.2015040200.0.0.tar.gz
## Using locally cached version of ./20150402-PRAD-CNVSNPHg19.txt
## gdac.broadinstitute.org_PRAD.Merge_cna__illuminahiseq_dnaseqc__hms_harvard_edu__Level_3__segmentation__seg.Level_3.2015040200.0.0.tar.gz
## Using locally cached version of ./20150402-PRAD-CNAseq.txt
## 
## gdac.broadinstitute.org_PRAD.Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2015040200.0.0.tar.gz
## http://gdac.broadinstitute.org/runs/stddata__2015_04_02/data/PRAD/20150402/gdac.broadinstitute.org_PRAD.Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2015040200.0.0.tar.gz
## File Size: ~2066MB
## File above won't be downloaded due to data size, RTCGAToolbox will skip this data!
## 
## 
## 
## 
## 
## gdac.broadinstitute.org_PRAD.Mutation_Packager_Calls.Level_3.2015040200.0.0.tar.gz
## gdac.broadinstitute.org_PRAD-TP.CopyNumber_Gistic2.Level_4.2014101700.0.0.tar.gz
```

# Extract available data types



```r
choices <- tolower(gsub("_", "", c("RNAseq_Gene", "miRNASeq_Gene",
             "RNAseq2_Gene_Norm", "CNA_SNP", "CNV_SNP", "CNA_Seq",
             "CNA_CGH", "Methylation", "Mutation", "mRNA_Array",
             "miRNA_Array", "RPPA")))
```


```r
dses <- lapply(choices, function(choice) try(extract(prad, choice, 
                                            clinic=TRUE),
                                             silent=TRUE))
names(dses) <- choices
dses
```

```
## $rnaseqgene
## [1] "Error in extract(prad, choice, clinic = TRUE) : \n  There is no data for that data type!\n"
## attr(,"class")
## [1] "try-error"
## attr(,"condition")
## <simpleError in extract(prad, choice, clinic = TRUE): There is no data for that data type!>
## 
## $mirnaseqgene
## ExpressionSet (storageMode: lockedEnvironment)
## assayData: 1046 features, 537 samples 
##   element names: exprs 
## protocolData: none
## phenoData
##   sampleNames: tcga-2a-a8vl-01 tcga-2a-a8vo-01 ... tcga-zg-a9ni-01
##     (537 total)
##   varLabels: vital_status days_to_death ... center (242 total)
##   varMetadata: labelDescription
## featureData: none
## experimentData: use 'experimentData(object)'
## Annotation:  
## 
## $rnaseq2genenorm
## ExpressionSet (storageMode: lockedEnvironment)
## assayData: 20501 features, 540 samples 
##   element names: exprs 
## protocolData: none
## phenoData
##   sampleNames: tcga-2a-a8vl-01 tcga-2a-a8vo-01 ... tcga-zg-a9ni-01
##     (540 total)
##   varLabels: vital_status days_to_death ... center (242 total)
##   varMetadata: labelDescription
## featureData: none
## experimentData: use 'experimentData(object)'
## Annotation:  
## 
## $cnasnp
## [1] "Error in extract(prad, choice, clinic = TRUE) : \n  object 'sourceName' not found\n"
## attr(,"class")
## [1] "try-error"
## attr(,"condition")
## <simpleError in extract(prad, choice, clinic = TRUE): object 'sourceName' not found>
## 
## $cnvsnp
## [1] "Error in extract(prad, choice, clinic = TRUE) : \n  object 'sourceName' not found\n"
## attr(,"class")
## [1] "try-error"
## attr(,"condition")
## <simpleError in extract(prad, choice, clinic = TRUE): object 'sourceName' not found>
## 
## $cnaseq
## [1] "Error in extract(prad, choice, clinic = TRUE) : \n  object 'sourceName' not found\n"
## attr(,"class")
## [1] "try-error"
## attr(,"condition")
## <simpleError in extract(prad, choice, clinic = TRUE): object 'sourceName' not found>
## 
## $cnacgh
## [1] "Error in extract(prad, choice, clinic = TRUE) : \n  There is no data for that data type!\n"
## attr(,"class")
## [1] "try-error"
## attr(,"condition")
## <simpleError in extract(prad, choice, clinic = TRUE): There is no data for that data type!>
## 
## $methylation
## [1] "Error in extract(prad, choice, clinic = TRUE) : \n  There is no data for that data type!\n"
## attr(,"class")
## [1] "try-error"
## attr(,"condition")
## <simpleError in extract(prad, choice, clinic = TRUE): There is no data for that data type!>
## 
## $mutation
## [1] "Error in extract(prad, choice, clinic = TRUE) : \n  object 'sourceName' not found\n"
## attr(,"class")
## [1] "try-error"
## attr(,"condition")
## <simpleError in extract(prad, choice, clinic = TRUE): object 'sourceName' not found>
## 
## $mrnaarray
## [1] "Error in extract(prad, choice, clinic = TRUE) : \n  There is no data for that data type!\n"
## attr(,"class")
## [1] "try-error"
## attr(,"condition")
## <simpleError in extract(prad, choice, clinic = TRUE): There is no data for that data type!>
## 
## $mirnaarray
## [1] "Error in extract(prad, choice, clinic = TRUE) : \n  There is no data for that data type!\n"
## attr(,"class")
## [1] "try-error"
## attr(,"condition")
## <simpleError in extract(prad, choice, clinic = TRUE): There is no data for that data type!>
## 
## $rppa
## [1] "Error in extract(prad, choice, clinic = TRUE) : \n  There is no data for that data type!\n"
## attr(,"class")
## [1] "try-error"
## attr(,"condition")
## <simpleError in extract(prad, choice, clinic = TRUE): There is no data for that data type!>
```


```r
eset <- extract(prad, "rnaseq2genenorm")
write.csv(exprs(eset), file="prad_rnaseq.csv")
write.csv(pData(eset), file="prad_clinical.csv")
saveRDS(eset, file="prad_eset.rds")
```

To load the eset again:


```r
eset <- readRDS("prad_eset.rds")
```

# Make a histogram of PSA (KLK3) expression


```r
hist(exprs(eset["KLK3", ]))
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-1.png) 

```r
hist(log(exprs(eset["KLK3", ])))
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-2.png) 

# See what clinical data are available by default:


```r
summary(pData(eset))
```

```
##  vital_status       days_to_death      days_to_last_followup
##  Length:540         Length:540         Length:540           
##  Class :character   Class :character   Class :character     
##  Mode  :character   Mode  :character   Mode  :character     
##                                                             
##                                                             
##                                                             
##                                                             
##  primary_site_of_disease neoplasm_diseasestage pathology_T_stage 
##  Length:540              Length:540            Length:540        
##  Class :character        Class :character      Class :character  
##  Mode  :character        Mode  :character      Mode  :character  
##                                                                  
##                                                                  
##                                                                  
##                                                                  
##  pathology_N_stage  dcc_upload_date       gender         
##  Length:540         Length:540         Length:540        
##  Class :character   Class :character   Class :character  
##  Mode  :character   Mode  :character   Mode  :character  
##                                                          
##                                                          
##                                                          
##                                                          
##  date_of_initial_pathologic_diagnosis radiation_therapy 
##  Length:540                           Length:540        
##  Class :character                     Class :character  
##  Mode  :character                     Mode  :character  
##                                                         
##                                                         
##                                                         
##                                                         
##  histological_type  radiations_radiation_regimenindication
##  Length:540         Length:540                            
##  Class :character   Class :character                      
##  Mode  :character   Mode  :character                      
##                                                           
##                                                           
##                                                           
##                                                           
##  completeness_of_resection number_of_lymph_nodes gleason_score_combined
##  Length:540                Length:540            Length:540            
##  Class :character          Class :character      Class :character      
##  Mode  :character          Mode  :character      Mode  :character      
##                                                                        
##                                                                        
##                                                                        
##                                                                        
##  gleason_score_primary gleason_score_secondary gleason_score     
##  Length:540            Length:540              Length:540        
##  Class :character      Class :character        Class :character  
##  Mode  :character      Mode  :character        Mode  :character  
##                                                                  
##                                                                  
##                                                                  
##                                                                  
##  psa_result_preop   days_to_preop_psa   psa_value        
##  Length:540         Length:540         Length:540        
##  Class :character   Class :character   Class :character  
##  Mode  :character   Mode  :character   Mode  :character  
##                                                          
##                                                          
##                                                          
##                                                          
##  days_to_psa            race            ethnicity        
##  Length:540         Length:540         Length:540        
##  Class :character   Class :character   Class :character  
##  Mode  :character   Mode  :character   Mode  :character  
##                                                          
##                                                          
##                                                          
##                                                          
##  batch_number       admin.batch_number
##  Length:540         184.61.0: 78      
##  Class :character   429.22.0: 68      
##  Mode  :character   91.57.0 : 68      
##                     389.30.0: 45      
##                     370.32.0: 42      
##                     357.34.0: 30      
##                     (Other) :209      
##                           admin.bcr   admin.day_of_dcc_upload
##  nationwide children's hospital:540   12:  4                 
##                                       26:534                 
##                                       29:  2                 
##                                                              
##                                                              
##                                                              
##                                                              
##                              admin.file_uuid admin.month_of_dcc_upload
##  016de7f0-1048-4add-ba29-85ba6a4dc4f7:  2    1:  2                    
##  060644b7-7ff0-481a-9f83-fc23dc16df91:  2    2:  2                    
##  076d02bb-3504-4139-91b7-5d52883cb39a:  2    3:536                    
##  09cc76ff-9742-4af8-8c0f-2998a5561dfc:  2                             
##  0ad65f96-3ab3-491e-8a17-1844474e3886:  2                             
##  0bb87419-ce94-4a3a-94b9-ca797cca10b7:  2                             
##  (Other)                             :528                             
##  patient.additional_studies patient.age_at_initial_pathologic_diagnosis
##  NA's:540                   66     : 36                                
##                             61     : 32                                
##                             64     : 31                                
##                             57     : 30                                
##                             63     : 29                                
##                             62     : 27                                
##                             (Other):355                                
##                          patient.bcr_patient_uuid
##  00e41e6a-9fe7-44f9-978b-7b05b179506a:  2        
##  0730216b-c201-443c-9092-81e23fd13c6c:  2        
##  07ccb4e9-a2a9-4bfd-97bd-13873b01977a:  2        
##  1426b4bf-f4a5-477c-bbbb-e0e656549764:  2        
##  181c6d37-9d8d-4ffa-827b-e3785450ff88:  2        
##  1d7ffbe7-6139-412b-8739-8082a63a2b33:  2        
##  (Other)                             :528        
##  patient.biochemical_recurrence
##  no  :402                      
##  yes : 61                      
##  NA's: 77                      
##                                
##                                
##                                
##                                
##                                  patient.bone_scan_results
##  abnormal (not related to prostate cancer)    : 23        
##  equivocal                                    :  9        
##  normal (no evidence of prostate cancer) [cm0]:231        
##  NA's                                         :277        
##                                                           
##                                                           
##                                                           
##  patient.clinical_cqcf.country patient.clinical_cqcf.days_to_consent
##  united states:441             31     : 14                          
##  germany      : 59             62     : 13                          
##  canada       : 27             42     : 12                          
##  australia    :  8             0      : 10                          
##  moldova      :  1             63     : 10                          
##  romania      :  1             (Other):451                          
##  (Other)      :  3             NA's   : 30                          
##  patient.clinical_cqcf.days_to_preop_psa
##  0      : 71                            
##  -33    :  8                            
##  -22    :  7                            
##  31     :  7                            
##  -18    :  6                            
##  (Other):403                            
##  NA's   : 38                            
##        patient.clinical_cqcf.ethnicity
##  hispanic or latino    :  5           
##  not hispanic or latino:424           
##  NA's                  :111           
##                                       
##                                       
##                                       
##                                       
##  patient.clinical_cqcf.gleason_score_combined
##  08:  1                                      
##  10: 11                                      
##  6 : 55                                      
##  7 :284                                      
##  8 : 75                                      
##  9 :114                                      
##                                              
##  patient.clinical_cqcf.gleason_score_primary
##  05:  1                                     
##  2 :  2                                     
##  3 :235                                     
##  4 :241                                     
##  5 : 61                                     
##                                             
##                                             
##  patient.clinical_cqcf.gleason_score_secondary
##  03:  1                                       
##  3 :167                                       
##  4 :285                                       
##  5 : 87                                       
##                                               
##                                               
##                                               
##                    patient.clinical_cqcf.histological_type
##  prostate adenocarcinoma  other subtype: 16               
##  prostate adenocarcinoma acinar type   :524               
##                                                           
##                                                           
##                                                           
##                                                           
##                                                           
##                                patient.clinical_cqcf.history_of_prior_malignancy
##  no                                                     :509                    
##  yes, history of prior malignancy                       : 28                    
##  yes, history of synchronous and or bilateral malignancy:  3                    
##                                                                                 
##                                                                                 
##                                                                                 
##                                                                                 
##  patient.clinical_cqcf.normal_tissue_anatomic_site
##  other          :170                              
##  seminal vesicle: 58                              
##  NA's           :312                              
##                                                   
##                                                   
##                                                   
##                                                   
##  patient.clinical_cqcf.normal_tissue_proximity
##  adjacent (< or = 2cm): 12                    
##  distal (>2cm)        :208                    
##  NA's                 :320                    
##                                               
##                                               
##                                               
##                                               
##  patient.clinical_cqcf.other_anatomic_site_normal_tissue
##  blood   :  5                                           
##  prostate:115                                           
##  NA's    :420                                           
##                                                         
##                                                         
##                                                         
##                                                         
##  patient.clinical_cqcf.psa_result_preop
##  4.6    : 14                           
##  5.2    : 12                           
##  5.6    : 11                           
##  5.5    : 10                           
##  7.4    : 10                           
##  (Other):480                           
##  NA's   :  3                           
##                     patient.clinical_cqcf.race
##  american indian or alaska native:  1         
##  asian                           : 12         
##  black or african american       : 64         
##  white                           :444         
##  NA's                            : 19         
##                                               
##                                               
##  patient.day_of_form_completion patient.days_to_birth
##  16     : 45                    -20805 :  4          
##  20     : 31                    -23725 :  4          
##  5      : 29                    -23360 :  3          
##  28     : 28                    -23458 :  3          
##  6      : 28                    -15884 :  2          
##  23     : 27                    (Other):513          
##  (Other):352                    NA's   : 11          
##  patient.days_to_bone_scan_performed patient.days_to_death
##  31     : 11                         1329   :  1          
##  10     : 10                         1855   :  1          
##  12     : 10                         3467   :  1          
##  13     : 10                         628    :  1          
##  14     :  9                         728    :  1          
##  (Other):208                         (Other):  2          
##  NA's   :282                         NA's   :533          
##  patient.days_to_diagnostic_computed_tomography_performed
##  11     : 12                                             
##  31     : 12                                             
##  14     : 11                                             
##  10     :  8                                             
##  30     :  8                                             
##  (Other):167                                             
##  NA's   :322                                             
##  patient.days_to_diagnostic_mri_performed
##  0      :  3                             
##  43     :  3                             
##  45     :  3                             
##  -73    :  2                             
##  -8     :  2                             
##  (Other): 69                             
##  NA's   :458                             
##  patient.days_to_first_biochemical_recurrence
##  189    :  3                                 
##  124    :  2                                 
##  131    :  2                                 
##  81     :  2                                 
##  1009   :  1                                 
##  (Other): 62                                 
##  NA's   :468                                 
##  patient.days_to_initial_pathologic_diagnosis
##  0   :510                                    
##  NA's: 30                                    
##                                              
##                                              
##                                              
##                                              
##                                              
##  patient.days_to_last_followup patient.diagnostic_ct_abd_pelvis_performed
##  62     :  6                   no  :205                                  
##  91     :  6                   yes :220                                  
##  31     :  5                   NA's:115                                  
##  822    :  5                                                             
##  285    :  4                                                             
##  (Other):507                                                             
##  NA's   :  7                                                             
##                           patient.diagnostic_ct_abd_pelvis_results.diagnostic_ct_abd_pelvis_result
##  equivocal                                                    :  5                                
##  extraprostatic extension  localized (e.g. seminal vesicles)  :  5                                
##  extraprostatic extension (regional lymphadenopathy)[e.g. cn1]:  6                                
##  no evidence of extraprostatic extension                      :202                                
##  NA's                                                         :322                                
##                                                                                                   
##                                                                                                   
##  patient.diagnostic_mri_performed
##  no  :332                        
##  yes : 85                        
##  NA's:123                        
##                                  
##                                  
##                                  
##                                  
##                                  patient.diagnostic_mri_results.diagnostic_mri_result
##  equivocal                                                 :  5                      
##  extraprostatic extension localized (e.g. seminal vesicles): 23                      
##  no evidence of extraprostatic extension                   : 53                      
##  NA's                                                      :459                      
##                                                                                      
##                                                                                      
##                                                                                      
##                                     patient.diagnostic_mri_results.diagnostic_mri_result-2
##  equivocal                                                     :  1                       
##  extraprostatic extension (regional lymphadenopathy) [e.g. cn1]:  3                       
##  NA's                                                          :536                       
##                                                                                           
##                                                                                           
##                                                                                           
##                                                                                           
##  patient.drugs.drug-2.bcr_drug_barcode
##  tcga-ch-5751-d37730:  1              
##  tcga-ch-5788-d37734:  1              
##  tcga-g9-7521-d48748:  1              
##  tcga-g9-a9s0-d60766:  1              
##  tcga-hi-7171-d29022:  1              
##  (Other)            : 27              
##  NA's               :508              
##                     patient.drugs.drug-2.bcr_drug_uuid
##  03320f22-90b1-4353-ac00-be36880265e5:  1             
##  06103a93-83eb-418a-9522-2ff6efb7c941:  1             
##  0af1ea65-d032-4370-8744-f0016808415b:  1             
##  13e4a7e5-1c0b-4702-95a8-801835da59e4:  1             
##  187f4262-f427-4804-994f-8b01714a052a:  1             
##  (Other)                             : 27             
##  NA's                                :508             
##  patient.drugs.drug-2.clinical_trail_drug_classification
##  n/a :  2                                               
##  NA's:538                                               
##                                                         
##                                                         
##                                                         
##                                                         
##                                                         
##  patient.drugs.drug-2.day_of_form_completion
##  13     :  7                                
##  25     :  5                                
##  27     :  3                                
##  2      :  2                                
##  20     :  2                                
##  (Other): 13                                
##  NA's   :508                                
##  patient.drugs.drug-2.days_to_drug_therapy_end
##  1036   :  1                                  
##  1228   :  1                                  
##  126    :  1                                  
##  1461   :  1                                  
##  1527   :  1                                  
##  (Other): 12                                  
##  NA's   :523                                  
##  patient.drugs.drug-2.days_to_drug_therapy_start
##  207    :  2                                    
##  1078   :  1                                    
##  115    :  1                                    
##  1163   :  1                                    
##  1219   :  1                                    
##  (Other): 25                                    
##  NA's   :509                                    
##  patient.drugs.drug-2.drug_name
##  lupron    :  9                
##  casodex   :  6                
##  zoladex   :  4                
##  eligard   :  2                
##  leuprolide:  2                
##  (Other)   :  9                
##  NA's      :508                
##          patient.drugs.drug-2.measure_of_response
##  clinical progressive disease:  2                
##  complete response           : 11                
##  partial response            :  3                
##  stable disease              :  2                
##  NA's                        :522                
##                                                  
##                                                  
##  patient.drugs.drug-2.month_of_form_completion
##  4      :  9                                  
##  3      :  5                                  
##  12     :  4                                  
##  8      :  4                                  
##  10     :  2                                  
##  (Other):  8                                  
##  NA's   :508                                  
##  patient.drugs.drug-2.regimen_indication
##  recurrence:  1                         
##  NA's      :539                         
##                                         
##                                         
##                                         
##                                         
##                                         
##  patient.drugs.drug-2.therapy_ongoing
##  no  : 18                            
##  yes : 14                            
##  NA's:508                            
##                                      
##                                      
##                                      
##                                      
##  patient.drugs.drug-2.therapy_types.therapy_type
##  hormone therapy: 32                            
##  NA's           :508                            
##                                                 
##                                                 
##                                                 
##                                                 
##                                                 
##  patient.drugs.drug-2.tx_on_clinical_trial
##  no  : 27                                 
##  yes :  4                                 
##  NA's:509                                 
##                                           
##                                           
##                                           
##                                           
##  patient.drugs.drug-2.year_of_form_completion
##  2012:  3                                    
##  2013:  1                                    
##  2014: 27                                    
##  2015:  1                                    
##  NA's:508                                    
##                                              
##                                              
##  patient.drugs.drug-3.bcr_drug_barcode
##  tcga-ch-5751-d37731:  1              
##  tcga-g9-7521-d48749:  1              
##  tcga-g9-a9s0-d60767:  1              
##  tcga-hi-7171-d29023:  1              
##  tcga-kk-a6e0-d62512:  1              
##  (Other)            :  7              
##  NA's               :528              
##                     patient.drugs.drug-3.bcr_drug_uuid
##  03feb539-c786-41e6-9134-3c83d062c466:  1             
##  057173bd-c7f4-4750-882a-b9b51145c98a:  1             
##  1572b2f1-99f4-4b83-b2e3-67c95c87a27e:  1             
##  47e9f8b2-f44d-49a1-aced-7b3ce938ad99:  1             
##  4a16143e-2d4f-4674-971c-3e045626f7aa:  1             
##  (Other)                             :  7             
##  NA's                                :528             
##  patient.drugs.drug-3.clinical_trail_drug_classification
##  n/a :  1                                               
##  NA's:539                                               
##                                                         
##                                                         
##                                                         
##                                                         
##                                                         
##  patient.drugs.drug-3.day_of_form_completion
##  13     :  1                                
##  16     :  1                                
##  18     :  1                                
##  2      :  1                                
##  20     :  1                                
##  (Other):  7                                
##  NA's   :528                                
##  patient.drugs.drug-3.days_to_drug_therapy_end
##  1339   :  1                                  
##  1630   :  1                                  
##  2170   :  1                                  
##  3467   :  1                                  
##  605    :  1                                  
##  (Other):  2                                  
##  NA's   :533                                  
##  patient.drugs.drug-3.days_to_drug_therapy_start
##  1183   :  1                                    
##  1287   :  1                                    
##  1600   :  1                                    
##  1897   :  1                                    
##  2431   :  1                                    
##  (Other):  6                                    
##  NA's   :529                                    
##  patient.drugs.drug-3.drug_name
##  casodex        :  4           
##  degarelix      :  1           
##  lupron         :  5           
##  taxotere       :  1           
##  zoledronic acid:  1           
##  NA's           :528           
##                                
##          patient.drugs.drug-3.measure_of_response
##  clinical progressive disease:  3                
##  complete response           :  1                
##  partial response            :  2                
##  NA's                        :534                
##                                                  
##                                                  
##                                                  
##  patient.drugs.drug-3.month_of_form_completion
##  4      :  2                                  
##  7      :  2                                  
##  8      :  2                                  
##  11     :  1                                  
##  12     :  1                                  
##  (Other):  4                                  
##  NA's   :528                                  
##  patient.drugs.drug-3.therapy_ongoing
##  no  :  8                            
##  yes :  4                            
##  NA's:528                            
##                                      
##                                      
##                                      
##                                      
##  patient.drugs.drug-3.therapy_types.therapy_type
##  ancillary      :  1                            
##  chemotherapy   :  1                            
##  hormone therapy: 10                            
##  NA's           :528                            
##                                                 
##                                                 
##                                                 
##  patient.drugs.drug-3.tx_on_clinical_trial
##  no  : 11                                 
##  NA's:529                                 
##                                           
##                                           
##                                           
##                                           
##                                           
##  patient.drugs.drug-3.year_of_form_completion
##  2012:  2                                    
##  2013:  1                                    
##  2014:  9                                    
##  NA's:528                                    
##                                              
##                                              
##                                              
##  patient.drugs.drug-4.bcr_drug_barcode
##  tcga-g9-a9s0-d60768:  1              
##  tcga-hi-7171-d29024:  1              
##  tcga-kk-a8ii-d62413:  1              
##  tcga-yl-a8sb-d67267:  1              
##  tcga-yl-a8sp-d59011:  1              
##  tcga-yl-a9wk-d58528:  1              
##  NA's               :534              
##                     patient.drugs.drug-4.bcr_drug_uuid
##  17d96cc4-b3a0-46c2-af8a-5454ad7e54af:  1             
##  3780b9f2-c2ca-4db2-a58e-6f6fa6a66a55:  1             
##  4e9d83b0-0e1a-4de9-8edc-6777986d8f98:  1             
##  4f7e50a0-c217-4635-a1e3-677d4b516094:  1             
##  5c3ec5b6-18d2-46be-9fd9-3703a2fa6fbb:  1             
##  770b5c2b-e349-463b-a7d8-ab2d0ee929fc:  1             
##  NA's                                :534             
##  patient.drugs.drug-4.day_of_form_completion
##  13  :  1                                   
##  16  :  1                                   
##  22  :  1                                   
##  25  :  1                                   
##  3   :  1                                   
##  30  :  1                                   
##  NA's:534                                   
##  patient.drugs.drug-4.days_to_drug_therapy_end
##  2170:  1                                     
##  2738:  1                                     
##  727 :  1                                     
##  NA's:537                                     
##                                               
##                                               
##                                               
##  patient.drugs.drug-4.days_to_drug_therapy_start
##  1609:  1                                       
##  1814:  1                                       
##  1891:  1                                       
##  2516:  1                                       
##  559 :  1                                       
##  727 :  1                                       
##  NA's:534                                       
##  patient.drugs.drug-4.drug_name
##  cabazitaxel:  1               
##  casodex    :  3               
##  docetaxel  :  1               
##  lupron     :  1               
##  NA's       :534               
##                                
##                                
##          patient.drugs.drug-4.measure_of_response
##  clinical progressive disease:  2                
##  partial response            :  1                
##  NA's                        :537                
##                                                  
##                                                  
##                                                  
##                                                  
##  patient.drugs.drug-4.month_of_form_completion
##  11  :  1                                     
##  2   :  1                                     
##  4   :  2                                     
##  6   :  1                                     
##  8   :  1                                     
##  NA's:534                                     
##                                               
##  patient.drugs.drug-4.therapy_ongoing
##  no  :  4                            
##  yes :  2                            
##  NA's:534                            
##                                      
##                                      
##                                      
##                                      
##  patient.drugs.drug-4.therapy_types.therapy_type
##  chemotherapy   :  2                            
##  hormone therapy:  4                            
##  NA's           :534                            
##                                                 
##                                                 
##                                                 
##                                                 
##  patient.drugs.drug-4.tx_on_clinical_trial
##  no  :  5                                 
##  NA's:535                                 
##                                           
##                                           
##                                           
##                                           
##                                           
##  patient.drugs.drug-4.year_of_form_completion
##  2012:  1                                    
##  2014:  5                                    
##  NA's:534                                    
##                                              
##                                              
##                                              
##                                              
##  patient.drugs.drug-5.bcr_drug_barcode
##  tcga-kk-a8ii-d62416:  1              
##  NA's               :539              
##                                       
##                                       
##                                       
##                                       
##                                       
##   patient.drugs.drug.bcr_drug_barcode
##  tcga-v1-a9o5-d70126:  2             
##  tcga-4l-aa1f-d67927:  1             
##  tcga-ch-5741-d37729:  1             
##  tcga-ch-5751-d37727:  1             
##  tcga-ch-5752-d37732:  1             
##  (Other)            : 60             
##  NA's               :474             
##                      patient.drugs.drug.bcr_drug_uuid
##  cf129345-f428-483b-bf72-c431f6bba535:  2            
##  00ca7115-7773-4773-9ac4-ffccd1ccf33b:  1            
##  0f689e8f-3720-4bb6-9b8d-3fd7bf4a4164:  1            
##  1383f0a9-3896-485b-8575-e3412a373b34:  1            
##  1855423a-82d5-40a1-8789-9cb1346acf7d:  1            
##  (Other)                             : 60            
##  NA's                                :474            
##  patient.drugs.drug.clinical_trail_drug_classification
##  hormonal treatment:  1                               
##  n/a               :  4                               
##  NA's              :535                               
##                                                       
##                                                       
##                                                       
##                                                       
##  patient.drugs.drug.day_of_form_completion
##  25     : 10                              
##  13     :  8                              
##  27     :  7                              
##  10     :  4                              
##  20     :  4                              
##  (Other): 33                              
##  NA's   :474                              
##  patient.drugs.drug.days_to_drug_therapy_end
##  1106   :  1                                
##  143    :  1                                
##  1455   :  1                                
##  162    :  1                                
##  1862   :  1                                
##  (Other): 21                                
##  NA's   :514                                
##  patient.drugs.drug.days_to_drug_therapy_start
##  151    :  4                                  
##  118    :  3                                  
##  179    :  2                                  
##  207    :  2                                  
##  209    :  2                                  
##  (Other): 53                                  
##  NA's   :474                                  
##  patient.drugs.drug.drug_name
##  lupron    : 16              
##  casodex   : 14              
##  eligard   :  6              
##  leuprolide:  4              
##  zoladex   :  4              
##  (Other)   : 21              
##  NA's      :475              
##           patient.drugs.drug.measure_of_response
##  clinical progressive disease:  2               
##  complete response           : 19               
##  partial response            :  2               
##  stable disease              :  1               
##  NA's                        :516               
##                                                 
##                                                 
##  patient.drugs.drug.month_of_form_completion
##  8      : 12                                
##  4      : 11                                
##  3      :  9                                
##  1      :  8                                
##  12     :  6                                
##  (Other): 20                                
##  NA's   :474                                
##  patient.drugs.drug.number_cycles patient.drugs.drug.prescribed_dose
##  13  :  1                         22.5:  1                          
##  NA's:539                         400 :  1                          
##                                   NA's:538                          
##                                                                     
##                                                                     
##                                                                     
##                                                                     
##  patient.drugs.drug.prescribed_dose_units
##  mg  :  2                                
##  NA's:538                                
##                                          
##                                          
##                                          
##                                          
##                                          
##  patient.drugs.drug.therapy_ongoing
##  no  : 26                          
##  yes : 39                          
##  NA's:475                          
##                                    
##                                    
##                                    
##                                    
##  patient.drugs.drug.therapy_types.therapy_type
##  chemotherapy   :  3                          
##  hormone therapy: 63                          
##  NA's           :474                          
##                                               
##                                               
##                                               
##                                               
##  patient.drugs.drug.tx_on_clinical_trial
##  no  : 58                               
##  yes :  6                               
##  NA's:476                               
##                                         
##                                         
##                                         
##                                         
##  patient.drugs.drug.year_of_form_completion              patient.ethnicity
##  2012:  7                                   not hispanic or latino:197    
##  2013:  4                                   NA's                  :343    
##  2014: 48                                                                 
##  2015:  7                                                                 
##  NA's:474                                                                 
##                                                                           
##                                                                           
##  patient.follow_ups.follow_up-2.additional_pharmaceutical_therapy
##  no  :  4                                                        
##  yes :  1                                                        
##  NA's:535                                                        
##                                                                  
##                                                                  
##                                                                  
##                                                                  
##  patient.follow_ups.follow_up-2.bcr_followup_barcode
##  tcga-ej-7781-f60365:  2                            
##  tcga-ej-7782-f60366:  2                            
##  tcga-ej-7783-f60370:  2                            
##  tcga-ej-7784-f60372:  2                            
##  tcga-ej-7793-f60362:  2                            
##  (Other)            : 26                            
##  NA's               :504                            
##              patient.follow_ups.follow_up-2.bcr_followup_uuid
##  0a531879-7216-4649-9d57-18ce9d6de7e4:  2                    
##  2fe6895e-1d5b-453a-b8ca-c51311f587d2:  2                    
##  45f62940-8821-43ac-bc5b-268a9ba14895:  2                    
##  601cf690-ca85-4f26-ae09-9ac8daba8d86:  2                    
##  84754d49-164b-4e0c-bb40-690d1f151cb3:  2                    
##  (Other)                             : 26                    
##  NA's                                :504                    
##  patient.follow_ups.follow_up-2.day_of_form_completion
##  6      : 19                                          
##  11     :  4                                          
##  5      :  4                                          
##  30     :  3                                          
##  20     :  2                                          
##  (Other):  4                                          
##  NA's   :504                                          
##  patient.follow_ups.follow_up-2.days_to_first_biochemical_recurrence
##  1443:  2                                                           
##  NA's:538                                                           
##                                                                     
##                                                                     
##                                                                     
##                                                                     
##                                                                     
##  patient.follow_ups.follow_up-2.days_to_last_followup
##  1167   :  2                                         
##  1378   :  2                                         
##  1443   :  2                                         
##  2465   :  2                                         
##  606    :  2                                         
##  (Other): 23                                         
##  NA's   :507                                         
##  patient.follow_ups.follow_up-2.days_to_new_tumor_event_after_initial_treatment
##  1285:  1                                                                      
##  425 :  1                                                                      
##  974 :  1                                                                      
##  NA's:537                                                                      
##                                                                                
##                                                                                
##                                                                                
##  patient.follow_ups.follow_up-2.followup_case_report_form_submission_reason
##  additional new tumor event    :  2                                        
##  scheduled follow-up submission: 34                                        
##  NA's                          :504                                        
##                                                                            
##                                                                            
##                                                                            
##                                                                            
##  patient.follow_ups.follow_up-2.followup_treatment_success
##  complete remission/response: 32                          
##  progressive disease        :  4                          
##  NA's                       :504                          
##                                                           
##                                                           
##                                                           
##                                                           
##  patient.follow_ups.follow_up-2.lost_follow_up
##  no  : 30                                     
##  yes :  4                                     
##  NA's:506                                     
##                                               
##                                               
##                                               
##                                               
##  patient.follow_ups.follow_up-2.month_of_form_completion
##  10  :  2                                               
##  2   :  2                                               
##  3   :  1                                               
##  4   :  3                                               
##  5   :  1                                               
##  6   : 27                                               
##  NA's:504                                               
##      patient.follow_ups.follow_up-2.new_neoplasm_event_type
##  biochemical evidence of disease:  2                       
##  new primary tumor              :  3                       
##  NA's                           :535                       
##                                                            
##                                                            
##                                                            
##                                                            
##  patient.follow_ups.follow_up-2.new_tumor_event_after_initial_treatment
##  no  : 27                                                              
##  yes :  5                                                              
##  NA's:508                                                              
##                                                                        
##                                                                        
##                                                                        
##                                                                        
##  patient.follow_ups.follow_up-2.person_neoplasm_cancer_status
##  tumor free: 29                                              
##  with tumor:  3                                              
##  NA's      :508                                              
##                                                              
##                                                              
##                                                              
##                                                              
##  patient.follow_ups.follow_up-2.primary_therapy_outcome_success
##  complete remission/response: 27                               
##  progressive disease        :  2                               
##  NA's                       :511                               
##                                                                
##                                                                
##                                                                
##                                                                
##  patient.follow_ups.follow_up-2.radiation_therapy
##  no  : 29                                        
##  NA's:511                                        
##                                                  
##                                                  
##                                                  
##                                                  
##                                                  
##  patient.follow_ups.follow_up-2.vital_status
##  alive: 32                                  
##  NA's :508                                  
##                                             
##                                             
##                                             
##                                             
##                                             
##  patient.follow_ups.follow_up-2.year_of_form_completion
##  2012:  2                                              
##  2013:  8                                              
##  2014: 26                                              
##  NA's:504                                              
##                                                        
##                                                        
##                                                        
##  patient.follow_ups.follow_up.additional_pharmaceutical_therapy
##  no  : 34                                                      
##  yes : 18                                                      
##  NA's:488                                                      
##                                                                
##                                                                
##                                                                
##                                                                
##  patient.follow_ups.follow_up.additional_radiation_therapy
##  no  : 35                                                 
##  yes : 17                                                 
##  NA's:488                                                 
##                                                           
##                                                           
##                                                           
##                                                           
##  patient.follow_ups.follow_up.bcr_followup_barcode
##  tcga-ch-5767-f28999:  2                          
##  tcga-ch-5768-f28998:  2                          
##  tcga-ej-7115-f64799:  2                          
##  tcga-ej-7123-f64800:  2                          
##  tcga-ej-7125-f64801:  2                          
##  (Other)            :409                          
##  NA's               :121                          
##               patient.follow_ups.follow_up.bcr_followup_uuid
##  00c8611c-dc2f-4a9c-9157-5fbdd688327e:  2                   
##  014bfbba-5bb0-4b05-9065-765c982720f7:  2                   
##  053bb5a0-578f-42d3-87fe-9ab36d9d8826:  2                   
##  09804433-748d-4bcd-a78c-591903061a38:  2                   
##  0cc76757-83da-4ae3-8ab6-328b0ce091a8:  2                   
##  (Other)                             :409                   
##  NA's                                :121                   
##  patient.follow_ups.follow_up.day_of_form_completion
##  24     : 39                                        
##  3      : 36                                        
##  11     : 35                                        
##  12     : 26                                        
##  10     : 24                                        
##  (Other):259                                        
##  NA's   :121                                        
##  patient.follow_ups.follow_up.days_to_death
##  1329:  1                                  
##  1855:  1                                  
##  3467:  1                                  
##  3502:  1                                  
##  728 :  1                                  
##  791 :  1                                  
##  NA's:534                                  
##  patient.follow_ups.follow_up.days_to_first_biochemical_recurrence
##  124    :  2                                                      
##  420    :  2                                                      
##  1009   :  1                                                      
##  1060   :  1                                                      
##  1068   :  1                                                      
##  (Other): 43                                                      
##  NA's   :490                                                      
##  patient.follow_ups.follow_up.days_to_last_followup
##  396    :  4                                       
##  1434   :  3                                       
##  425    :  3                                       
##  442    :  3                                       
##  686    :  3                                       
##  (Other):395                                       
##  NA's   :129                                       
##  patient.follow_ups.follow_up.days_to_new_tumor_event_after_initial_treatment
##  420    :  2                                                                 
##  1006   :  1                                                                 
##  1009   :  1                                                                 
##  1060   :  1                                                                 
##  1068   :  1                                                                 
##  (Other): 36                                                                 
##  NA's   :498                                                                 
##  patient.follow_ups.follow_up.days_to_second_biochemical_recurrence
##  1523   :  1                                                       
##  1687   :  1                                                       
##  1756   :  1                                                       
##  1934   :  1                                                       
##  2019   :  1                                                       
##  (Other):  4                                                       
##  NA's   :531                                                       
##  patient.follow_ups.follow_up.days_to_third_biochemical_recurrence
##  2011:  1                                                         
##  631 :  1                                                         
##  655 :  1                                                         
##  NA's:537                                                         
##                                                                   
##                                                                   
##                                                                   
##  patient.follow_ups.follow_up.followup_case_report_form_submission_reason
##  additional new tumor event    : 18                                      
##  scheduled follow-up submission:399                                      
##  NA's                          :123                                      
##                                                                          
##                                                                          
##                                                                          
##                                                                          
##  patient.follow_ups.follow_up.followup_treatment_success
##  complete remission/response:312                        
##  partial remission/response : 13                        
##  progressive disease        : 15                        
##  stable disease             : 35                        
##  NA's                       :165                        
##                                                         
##                                                         
##  patient.follow_ups.follow_up.lost_follow_up
##  no  :360                                   
##  yes : 33                                   
##  NA's:147                                   
##                                             
##                                             
##                                             
##                                             
##  patient.follow_ups.follow_up.month_of_form_completion
##  9      :117                                          
##  2      : 74                                          
##  8      : 60                                          
##  11     : 44                                          
##  1      : 33                                          
##  (Other): 91                                          
##  NA's   :121                                          
##  patient.follow_ups.follow_up.new_neoplasm_event_occurrence_anatomic_site
##  bone          :  4                                                      
##  lung          :  1                                                      
##  other, specify:  4                                                      
##  NA's          :531                                                      
##                                                                          
##                                                                          
##                                                                          
##       patient.follow_ups.follow_up.new_neoplasm_event_type
##  biochemical evidence of disease: 45                      
##  distant metastasis             :  5                      
##  locoregional recurrence        :  3                      
##  new primary tumor              :  4                      
##  NA's                           :483                      
##                                                           
##                                                           
##  patient.follow_ups.follow_up.new_neoplasm_occurrence_anatomic_site_text
##  biochemical recurrence only:  2                                        
##  lung and liver             :  1                                        
##  neck                       :  1                                        
##  NA's                       :536                                        
##                                                                         
##                                                                         
##                                                                         
##  patient.follow_ups.follow_up.new_tumor_event_after_initial_treatment
##  no  :333                                                            
##  yes : 60                                                            
##  NA's:147                                                            
##                                                                      
##                                                                      
##                                                                      
##                                                                      
##  patient.follow_ups.follow_up.patient_death_reason
##  prostate cancer:  2                              
##  NA's           :538                              
##                                                   
##                                                   
##                                                   
##                                                   
##                                                   
##  patient.follow_ups.follow_up.person_neoplasm_cancer_status
##  tumor free:300                                            
##  with tumor: 41                                            
##  NA's      :199                                            
##                                                            
##                                                            
##                                                            
##                                                            
##  patient.follow_ups.follow_up.primary_therapy_outcome_success
##  complete remission/response:328                             
##  partial remission/response : 27                             
##  progressive disease        : 24                             
##  stable disease             : 22                             
##  NA's                       :139                             
##                                                              
##                                                              
##  patient.follow_ups.follow_up.radiation_therapy
##  no  :353                                      
##  yes : 49                                      
##  NA's:138                                      
##                                                
##                                                
##                                                
##                                                
##  patient.follow_ups.follow_up.source_of_patient_death_reason
##  medical record:  5                                         
##  NA's          :535                                         
##                                                             
##                                                             
##                                                             
##                                                             
##                                                             
##  patient.follow_ups.follow_up.targeted_molecular_therapy
##  no  :353                                               
##  yes : 46                                               
##  NA's:141                                               
##                                                         
##                                                         
##                                                         
##                                                         
##  patient.follow_ups.follow_up.tumor_progression_post_ht
##  no  :  4                                              
##  yes : 12                                              
##  NA's:524                                              
##                                                        
##                                                        
##                                                        
##                                                        
##  patient.follow_ups.follow_up.type_ofprogression_of_disease_ht
##  biochemical recurrence: 10                                   
##  distant metastasis    :  3                                   
##  NA's                  :527                                   
##                                                               
##                                                               
##                                                               
##                                                               
##  patient.follow_ups.follow_up.vital_status
##  alive:411                                
##  dead :  6                                
##  NA's :123                                
##                                           
##                                           
##                                           
##                                           
##  patient.follow_ups.follow_up.year_of_form_completion
##  2012: 29                                            
##  2013:116                                            
##  2014:218                                            
##  2015: 56                                            
##  NA's:121                                            
##                                                      
##                                                      
##                           patient.histological_type
##  prostate adenocarcinoma  other subtype: 16        
##  prostate adenocarcinoma acinar type   :524        
##                                                    
##                                                    
##                                                    
##                                                    
##                                                    
##  patient.history_of_neoadjuvant_treatment patient.icd_10
##  no  :539                                 c61  :255     
##  NA's:  1                                 c61.9:284     
##                                           NA's :  1     
##                                                         
##                                                         
##                                                         
##                                                         
##  patient.icd_o_3_histology
##  8140/3:284               
##  8255/3:  3               
##  8480/3:  1               
##  8500/3: 10               
##  8550/3:242               
##                           
##                           
##           patient.initial_pathologic_diagnosis_method patient.laterality
##  core needle biopsy             :523                  bilateral:474     
##  transurethral resection (turbt): 12                  left     : 20     
##  NA's                           :  5                  right    : 37     
##                                                       NA's     :  9     
##                                                                         
##                                                                         
##                                                                         
##  patient.lymphnodes_examined patient.month_of_form_completion
##  no  : 79                    3      :168                     
##  yes :458                    4      :134                     
##  NA's:  3                    7      : 59                     
##                              8      : 44                     
##                              1      : 35                     
##                              5      : 30                     
##                              (Other): 70                     
##  patient.new_tumor_events.new_tumor_event.additional_pharmaceutical_therapy
##  no  : 23                                                                  
##  yes : 19                                                                  
##  NA's:498                                                                  
##                                                                            
##                                                                            
##                                                                            
##                                                                            
##  patient.new_tumor_events.new_tumor_event.additional_radiation_therapy
##  no  : 21                                                             
##  yes : 23                                                             
##  NA's:496                                                             
##                                                                       
##                                                                       
##                                                                       
##                                                                       
##  patient.new_tumor_events.new_tumor_event.days_to_new_tumor_event_after_initial_treatment
##  124    :  2                                                                             
##  1006   :  1                                                                             
##  1009   :  1                                                                             
##  1022   :  1                                                                             
##  105    :  1                                                                             
##  (Other): 56                                                                             
##  NA's   :478                                                                             
##  patient.new_tumor_events.new_tumor_event.new_neoplasm_event_occurrence_anatomic_site
##  bone                              :  2                                              
##  non-regional / distant lymph nodes:  3                                              
##  other, specify                    :  1                                              
##  NA's                              :534                                              
##                                                                                      
##                                                                                      
##                                                                                      
##  patient.new_tumor_events.new_tumor_event.new_neoplasm_event_type
##  biochemical evidence of disease: 58                             
##  distant metastasis             :  4                             
##  locoregional recurrence        :  4                             
##  new primary tumor              :  3                             
##  NA's                           :471                             
##                                                                  
##                                                                  
##  patient.new_tumor_events.new_tumor_event.new_neoplasm_occurrence_anatomic_site_text
##  melanoma in situ of the right midback:  1                                          
##  NA's                                 :539                                          
##                                                                                     
##                                                                                     
##                                                                                     
##                                                                                     
##                                                                                     
##  patient.new_tumor_events.new_tumor_event.tumor_progression_post_ht
##  no  : 39                                                          
##  yes :  8                                                          
##  NA's:493                                                          
##                                                                    
##                                                                    
##                                                                    
##                                                                    
##  patient.new_tumor_events.new_tumor_event.type_of_progression_after_ht
##  biochemical recurrence:  6                                           
##  distant metastases    :  2                                           
##  NA's                  :532                                           
##                                                                       
##                                                                       
##                                                                       
##                                                                       
##  patient.new_tumor_events.new_tumor_event_after_initial_treatment
##  no  :249                                                        
##  yes : 70                                                        
##  NA's:221                                                        
##                                                                  
##                                                                  
##                                                                  
##                                                                  
##  patient.number_of_lymphnodes_examined
##  4      : 35                          
##  2      : 31                          
##  6      : 28                          
##  10     : 27                          
##  5      : 27                          
##  (Other):315                          
##  NA's   : 77                          
##  patient.number_of_lymphnodes_positive_by_he
##  0      :368                                
##  1      : 42                                
##  2      : 13                                
##  3      : 11                                
##  6      :  4                                
##  (Other):  9                                
##  NA's   : 93                                
##               patient.patient_death_reason patient.patient_id
##  other non-malignant disease:  1           5761   :  2       
##  prostate cancer            :  2           5767   :  2       
##  NA's                       :537           5768   :  2       
##                                            5769   :  2       
##                                            6333   :  2       
##                                            6342   :  2       
##                                            (Other):528       
##  patient.person_neoplasm_cancer_status
##  tumor free:354                       
##  with tumor: 82                       
##  NA's      :104                       
##                                       
##                                       
##                                       
##                                       
##         patient.primary_therapy_outcome_success
##  complete remission/response:161               
##  partial remission/response : 25               
##  progressive disease        : 14               
##  stable disease             : 20               
##  NA's                       :320               
##                                                
##                                                
##                                                 patient.prior_dx
##  no                                                     :510    
##  yes                                                    : 13    
##  yes, history of prior malignancy                       : 16    
##  yes, history of synchronous and or bilateral malignancy:  1    
##                                                                 
##                                                                 
##                                                                 
##                     patient.race patient.radiation_therapy
##  asian                    :  2   no  :241                 
##  black or african american: 11   yes : 40                 
##  white                    :189   NA's:259                 
##  NA's                     :338                            
##                                                           
##                                                           
##                                                           
##  patient.radiations.radiation-2.anatomic_treatment_site
##  distant recurrence :  1                               
##  primary tumor field:  1                               
##  NA's               :538                               
##                                                        
##                                                        
##                                                        
##                                                        
##  patient.radiations.radiation-2.bcr_radiation_barcode
##  tcga-fc-a5ob-r61440:  1                             
##  tcga-g9-a9s0-r60771:  1                             
##  NA's               :538                             
##                                                      
##                                                      
##                                                      
##                                                      
##         patient.radiations.radiation-2.measure_of_response
##  radiographic progressive disease:  1                     
##  NA's                            :539                     
##                                                           
##                                                           
##                                                           
##                                                           
##                                                           
##  patient.radiations.radiation-2.month_of_form_completion
##  6   :  2                                               
##  NA's:538                                               
##                                                         
##                                                         
##                                                         
##                                                         
##                                                         
##  patient.radiations.radiation.anatomic_treatment_site
##  distant recurrence :  1                             
##  distant site       :  1                             
##  local recurrence   : 12                             
##  primary tumor field: 42                             
##  regional site      : 11                             
##  NA's               :473                             
##                                                      
##  patient.radiations.radiation.bcr_radiation_barcode
##  tcga-hc-7747-r60209:  2                           
##  tcga-v1-a9o5-r70127:  2                           
##  tcga-ej-5524-r46761:  1                           
##  tcga-ej-a65f-r57144:  1                           
##  tcga-fc-a4ji-r48399:  1                           
##  (Other)            : 64                           
##  NA's               :469                           
##              patient.radiations.radiation.bcr_radiation_uuid
##  b13f9180-c4f6-412f-b6d0-a9ba658b3d5f:  2                   
##  c8610195-3eb5-434d-b13d-e6ca6891b2b8:  2                   
##  086dc138-5deb-45f5-9566-81593b3ae273:  1                   
##  0b360fc8-c6ac-47ef-a777-bd95c91dd297:  1                   
##  0c19b3e6-777d-4f3a-9a5f-ad29a5cbc7ba:  1                   
##  (Other)                             : 64                   
##  NA's                                :469                   
##  patient.radiations.radiation.day_of_form_completion
##  25     : 14                                        
##  13     :  7                                        
##  14     :  5                                        
##  28     :  5                                        
##  29     :  5                                        
##  (Other): 35                                        
##  NA's   :469                                        
##  patient.radiations.radiation.days_to_radiation_therapy_end
##  210    :  2                                               
##  358    :  2                                               
##  417    :  2                                               
##  851    :  2                                               
##  1128   :  1                                               
##  (Other): 59                                               
##  NA's   :472                                               
##  patient.radiations.radiation.days_to_radiation_therapy_start
##  162    :  2                                                 
##  206    :  2                                                 
##  208    :  2                                                 
##  322    :  2                                                 
##  103    :  1                                                 
##  (Other): 60                                                 
##  NA's   :471                                                 
##          patient.radiations.radiation.measure_of_response
##  complete response               : 38                    
##  partial response                :  9                    
##  radiographic progressive disease:  6                    
##  stable disease                  :  7                    
##  NA's                            :480                    
##                                                          
##                                                          
##  patient.radiations.radiation.month_of_form_completion
##  8      : 15                                          
##  3      : 13                                          
##  1      :  9                                          
##  4      :  9                                          
##  7      :  7                                          
##  (Other): 18                                          
##  NA's   :469                                          
##  patient.radiations.radiation.numfractions
##  33     : 12                              
##  35     : 11                              
##  37     :  7                              
##  38     :  5                              
##  34     :  3                              
##  (Other): 13                              
##  NA's   :489                              
##  patient.radiations.radiation.radiation_dosage
##  66     : 14                                  
##  70     : 12                                  
##  6660   :  3                                  
##  6840   :  3                                  
##  6480   :  2                                  
##  (Other): 20                                  
##  NA's   :486                                  
##  patient.radiations.radiation.radiation_treatment_ongoing
##  no  : 69                                                
##  yes :  2                                                
##  NA's:469                                                
##                                                          
##                                                          
##                                                          
##                                                          
##  patient.radiations.radiation.radiation_type
##  external     : 66                          
##  external beam:  2                          
##  NA's         :472                          
##                                             
##                                             
##                                             
##                                             
##  patient.radiations.radiation.units
##  cgy : 16                          
##  gy  : 36                          
##  NA's:488                          
##                                    
##                                    
##                                    
##                                    
##  patient.radiations.radiation.year_of_form_completion
##  2012:  2                                            
##  2013:  9                                            
##  2014: 52                                            
##  2015:  8                                            
##  NA's:469                                            
##                                                      
##                                                      
##  patient.residual_tumor patient.source_of_patient_death_reason
##  r0  :354               death certificate:  1                 
##  r1  :149               medical record   :  4                 
##  r2  :  6               NA's             :535                 
##  rx  : 15                                                     
##  NA's: 16                                                     
##                                                               
##                                                               
##  patient.stage_event.gleason_grading.gleason_score
##  10:  3                                           
##  6 : 50                                           
##  7 :286                                           
##  8 : 65                                           
##  9 :136                                           
##                                                   
##                                                   
##  patient.stage_event.gleason_grading.primary_pattern
##  2:  2                                              
##  3:225                                              
##  4:264                                              
##  5: 49                                              
##                                                     
##                                                     
##                                                     
##  patient.stage_event.gleason_grading.secondary_pattern
##  3:171                                                
##  4:262                                                
##  5:107                                                
##                                                       
##                                                       
##                                                       
##                                                       
##  patient.stage_event.gleason_grading.tertiary_pattern
##  0   :  1                                            
##  2   :  3                                            
##  3   : 82                                            
##  4   : 41                                            
##  5   : 51                                            
##  NA's:362                                            
##                                                      
##  patient.stage_event.psa.days_to_psa patient.stage_event.psa.psa_value
##  172    :  4                         0.1    :165                      
##  204    :  4                         0.03   : 42                      
##  285    :  4                         0.01   : 40                      
##  91     :  4                         0      : 18                      
##  128    :  3                         .01    : 11                      
##  (Other):466                         (Other):203                      
##  NA's   : 55                         NA's   : 61                      
##  patient.stage_event.tnm_categories.clinical_categories.clinical_m
##  m0  :487                                                         
##  m1a :  1                                                         
##  m1c :  1                                                         
##  NA's: 51                                                         
##                                                                   
##                                                                   
##                                                                   
##  patient.stage_event.tnm_categories.clinical_categories.clinical_t
##  t1c    :193                                                      
##  t2a    : 59                                                      
##  t2b    : 57                                                      
##  t2c    : 57                                                      
##  t3a    : 37                                                      
##  (Other): 34                                                      
##  NA's   :103                                                      
##  patient.stage_event.tnm_categories.pathologic_categories.pathologic_n
##  n0  :386                                                             
##  n1  : 79                                                             
##  NA's: 75                                                             
##                                                                       
##                                                                       
##                                                                       
##                                                                       
##  patient.stage_event.tnm_categories.pathologic_categories.pathologic_t
##  t2a : 15                                                             
##  t2b : 12                                                             
##  t2c :190                                                             
##  t3a :172                                                             
##  t3b :135                                                             
##  t4  : 13                                                             
##  NA's:  3                                                             
##  patient.targeted_molecular_therapy
##  no  :239                          
##  yes : 42                          
##  NA's:259                          
##                                    
##                                    
##                                    
##                                    
##  patient.tissue_prospective_collection_indicator
##  no :268                                        
##  yes:272                                        
##                                                 
##                                                 
##                                                 
##                                                 
##                                                 
##  patient.tissue_retrospective_collection_indicator
##  no :272                                          
##  yes:268                                          
##                                                   
##                                                   
##                                                   
##                                                   
##                                                   
##  patient.tissue_source_site patient.tumor_levels.tumor_level
##  ej     :120                apex  :311                      
##  hc     : 74                base  : 19                      
##  g9     : 55                middle: 49                      
##  kk     : 45                NA's  :161                      
##  ch     : 36                                                
##  v1     : 30                                                
##  (Other):180                                                
##  patient.tumor_levels.tumor_level-2 patient.tumor_levels.tumor_level-3
##  base  : 83                         base:170                          
##  middle:211                         NA's:370                          
##  NA's  :246                                                           
##                                                                       
##                                                                       
##                                                                       
##                                                                       
##  patient.vital_status patient.year_of_form_completion
##  alive:533            2011:  2                       
##  dead :  7            2012:229                       
##                       2013: 62                       
##                       2014:235                       
##                       2015: 12                       
##                                                      
##                                                      
##  patient.year_of_initial_pathologic_diagnosis
##  2011   :131                                 
##  2012   : 71                                 
##  2010   : 63                                 
##  2013   : 62                                 
##  2009   : 61                                 
##  (Other):122                                 
##  NA's   : 30                                 
##                   patient.zone_of_origin sample_type       
##  central zone                :  5        Length:540        
##  overlapping / multiple zones:125        Class :character  
##  peripheral zone             :136        Mode  :character  
##  transition zone             :  8                          
##  NA's                        :266                          
##                                                            
##                                                            
##  sample_code            vial             portion         
##  Length:540         Length:540         Length:540        
##  Class :character   Class :character   Class :character  
##  Mode  :character   Mode  :character   Mode  :character  
##                                                          
##                                                          
##                                                          
##                                                          
##    analyte             plate              center         
##  Length:540         Length:540         Length:540        
##  Class :character   Class :character   Class :character  
##  Mode  :character   Mode  :character   Mode  :character  
##                                                          
##                                                          
##                                                          
## 
```

# Look for association between KLK3 expression and clinical PSA


```r
psadat <- data.frame(psa=as.numeric(as.character(eset$patient.stage_event.psa.psa_value)),
                     klk3=t(exprs(eset["KLK3", ])))
psadat.complete <- psadat[complete.cases(psadat), ]
plot(KLK3 ~ psa, data=psadat.complete, xlab="clinical PSA", ylab="KLK3 tumor expression", log="xy")
```

```
## Warning in xy.coords(x, y, xlabel, ylabel, log): 29 x values <= 0 omitted
## from logarithmic plot
```

```r
fit <- lowess(x=psadat.complete$psa, y=psadat.complete$KLK3)
lines(fit, col="red", lw=3)
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-1.png) 

## Is there an association between PSA in the urine and KLK3 in the tumor?


```r
cor.test(x=psadat$KLK3,  y=psadat$psa, method="spearman")
```

```
## Warning in cor.test.default(x = psadat$KLK3, y = psadat$psa, method =
## "spearman"): Cannot compute exact p-value with ties
```

```
## 
## 	Spearman's rank correlation rho
## 
## data:  psadat$KLK3 and psadat$psa
## S = 21183000, p-value = 0.0005893
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##        rho 
## -0.1564584
```


# Clinical psa as a function of ethnicity

Need to get the better race variable with complete data from Tiffany:

```r
table(eset$race)
```

```
## 
##                     asian black or african american 
##                         2                        11 
##                     white 
##                       189
```

```r
boxplot(psadat$psa ~ eset$race, ylab="PSA")
```

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-11-1.png) 

