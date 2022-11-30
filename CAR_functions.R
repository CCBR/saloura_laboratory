############################################################
# QC Analysis
############################################################
main_prep_qc_secondary<-function(contrast_id,exclusionlist=""){
  #set conditions
  condition1=strsplit(contrast_id,"_vs_")[[1]][1]
  condition2=strsplit(contrast_id,"_vs_")[[1]][2]
  
  # filter based off of params
  sampleinfo=subset(groups_df,group==condition1 | group==condition2)
  rownames(sampleinfo)=NULL
  sampleinfo$group = relevel(as.factor(sampleinfo$group),condition2)
  
  # generate rawcounts
  rawcounts=format_counts_matrix(contrast_id,sampleinfo)
  
  # filter
  raw_counts=ceiling(rawcounts)
  cpm_counts=edgeR::cpm(as.matrix(rawcounts))
  log_cpm_counts=log2(cpm_counts)
  keep=rowSums(cpm_counts>0.5)>2
  
  filtered=raw_counts[keep,]
  colnames(filtered)=shorten_sample_id(colnames(filtered))
  
  # run deseq
  dds=run_deseq_analysis(filtered,sampleinfo)
  
  #save results as df
  result_dds <- results(dds)
  results_df <- as.data.frame(result_dds)
  fpath=paste0(output_car_dir,"DESEQ2_res_",contrast_id,".csv")
  write.csv(results_df,fpath)
  
  # annotate res
  peakAnno=peak_annotation(result_dds,contrast_id)
  
  return(peakAnno)
}