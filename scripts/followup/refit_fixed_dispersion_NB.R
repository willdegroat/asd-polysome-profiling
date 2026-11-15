suppressPackageStartupMessages(library(DESeq2))
root <- 'manuscript/followup'
counts<-as.matrix(read.csv(file.path(root,'count_model_reconstructed_counts.csv'),row.names=1,check.names=FALSE))
eligible<-read.csv(file.path(root,'architecture_eligible_analysis_population.csv'),stringsAsFactors=FALSE)$transcriptId
counts<-counts[eligible,,drop=FALSE]
storage.mode(counts)<-'integer'
sf<-read.delim('results/transcriptRedistributionSizeFactors.tsv',check.names=FALSE)
sfv<-setNames(sf$sizeFactor,paste(sf$participantId,sf$fraction,sep='_'))
diag<-read.delim('results/transcriptRedistribution.diagnostics.tsv',check.names=FALSE)
disp<-setNames(diag$dispersion,diag$transcriptId)[rownames(counts)]
meta<-data.frame(sampleId=colnames(counts),participantId=sub('_(input|monosome|polysome)$','',colnames(counts)),fraction=sub('^[^_]+_','',colnames(counts)),row.names=colnames(counts))
meta$condition<-ifelse(meta$participantId %in% c('1851','1861','1916'),'ASD','Unaffected')
for(omit in c('none','07c06','07c13','07c16','1851','1861','1916')) {
 if(file.exists(file.path(root,paste0('fixed_dispersion_NB_omit_',omit,'.csv')))) {cat('Using completed refit:',omit,'\n');next}
 m<-meta[meta$participantId!=omit,,drop=FALSE];m$participantId<-factor(m$participantId);m$fraction<-factor(m$fraction,levels=c('input','monosome','polysome'))
 mm<-model.matrix(~participantId+fraction,data=m)
 mm<-cbind(mm,ASD_Monosome=as.numeric(m$condition=='ASD' & m$fraction=='monosome'),ASD_Polysome=as.numeric(m$condition=='ASD' & m$fraction=='polysome'))
 stopifnot(qr(mm)$rank==ncol(mm))
 m$ASD_Monosome<-as.numeric(m$condition=='ASD' & m$fraction=='monosome');m$ASD_Polysome<-as.numeric(m$condition=='ASD' & m$fraction=='polysome')
 dds<-DESeqDataSetFromMatrix(countData=counts[,rownames(m),drop=FALSE],colData=m,design=~participantId+fraction+ASD_Monosome+ASD_Polysome)
 sizeFactors(dds)<-sfv[rownames(m)];dispersions(dds)<-disp
 dds<-nbinomWaldTest(dds,betaTol=1e-8,maxit=250,minmu=.5,quiet=TRUE)
 cn<-rep(0,ncol(mm));cn[ncol(mm)-1]<- -1;cn[ncol(mm)]<-1
 r<-as.data.frame(results(dds,contrast=cn,cooksCutoff=FALSE,independentFiltering=FALSE))
 r$transcriptId<-rownames(r);r$betaConverged<-mcols(dds)$betaConv
 write.csv(r,file.path(root,paste0('fixed_dispersion_NB_omit_',omit,'.csv')),row.names=FALSE)
 cat('Completed original-design fixed-dispersion NB fit; omitted culture:',omit,'; nonconverged:',sum(!r$betaConverged),'\n')
}
writeLines(capture.output(sessionInfo()),file.path(root,'NB_sensitivity_R_session.txt'))
