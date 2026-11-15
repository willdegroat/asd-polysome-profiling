root <- 'manuscript/followup'
files<-c('new_model_results.csv','additional_model_results.csv','covariates_and_individual_cultures.csv','fixed_dispersion_NB_architecture_models.csv','all_omissions_NB_architecture_models.csv')
tables<-lapply(files,function(f) {a<-read.csv(file.path(root,f),stringsAsFactors=FALSE);a$sourceFile<-f;a})
cols<-unique(unlist(lapply(tables,names)))
tables<-lapply(tables,function(a) {for(k in setdiff(cols,names(a))) a[[k]]<-NA;a[cols]})
z<-do.call(rbind,tables)
z$pValueExploratoryUnadjusted<-2*pt(-abs(z$tStatistic),df=z$df)
z$confidenceIntervalLow<-z$standardizedCoefficient-qt(.975,z$df)*z$standardError
z$confidenceIntervalHigh<-z$standardizedCoefficient+qt(.975,z$df)*z$standardError
z$feature[is.na(z$feature)]<-'logThreePrimeUtrLength'
z$outcomeUnitsPerDoubledFeature<-z$outcomeUnitsPerSD/z$featureSD*log(2)
z$outcomeUnitsPerDoubledFeature_CILow<-z$confidenceIntervalLow*z$outcomeSD/z$featureSD*log(2)
z$outcomeUnitsPerDoubledFeature_CIHigh<-z$confidenceIntervalHigh*z$outcomeSD/z$featureSD*log(2)
z$BHAcrossAllExploratoryModels<-NULL
z$log2EffectPerDoubledUTR<-NULL
write.csv(z,file.path(root,'all_models_with_inference.csv'),row.names=FALSE,na='')
nb<-z[z$sourceFile=='all_omissions_NB_architecture_models.csv' & z$outcome=='NBWald',]
write.csv(nb,file.path(root,'NB_omissions_Wald_summary.csv'),row.names=FALSE,na='')
print(nb[c('model','transcripts','genes','standardizedCoefficient','confidenceIntervalLow','confidenceIntervalHigh','pValueExploratoryUnadjusted')],row.names=FALSE)
