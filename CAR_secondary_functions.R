############################################################
# QC Analysis
############################################################
format_counts_matrix<-function(contrast_id,sampleinfo){
  # set variables
  rawcountsmatrix=paste0(contrast_subpath,contrast_id,extensions,"/",
                         contrast_id,extensions,"_",method,"countsmatrix.txt")
  # prep counts
  rawcounts = read.csv(rawcountsmatrix,
                       header = TRUE,sep="\t",
                       comment.char = "#", 
                       strip.white = TRUE,
                       check.names = FALSE,
                       colClasses = "character")
  rawcounts = as.data.frame(rawcounts)
  rawcounts %>% column_to_rownames("peakID") -> rawcounts
  
  # convert character to numeric to integer
  x = matrix(as.numeric(as.matrix(rawcounts)),ncol=ncol(rawcounts))
  x = matrix(mapply(x,FUN=as.integer),ncol=ncol(rawcounts))
  x = as.data.frame(x)
  colnames(x) = colnames(rawcounts)
  rownames(x) = rownames(rawcounts)
  rawcounts = x
  
  return(rawcounts)
}

shorten_sample_id<-function(input_id){
  output_id=gsub("_350K_0_35ng_0_25_HCHO_2m","",input_id)
  output_id=gsub("_0.25HCHO_500K","",output_id)
  output_id=gsub("MOC1_","",output_id)
  output_id=gsub("_10","",output_id)
  output_id=gsub("Smyd3_v","v",output_id)
  output_id=gsub("_2_","_",output_id)
  output_id=gsub("_CRISPR_","_",output_id)
  #output_id=strsplit(output_id,"__")[[1]][1]
  return(output_id)
}

run_deseq_analysis<-function(rawcounts,sampleinfo){
  # set variables
  bbpaths=paste0(contrast_subpath,"/bed_bedgraph_paths.tsv")
  
  # run DESEQ
  dds <- DESeqDataSetFromMatrix(countData = as.matrix(rawcounts),
                                colData = sampleinfo[,c("sampleid","group")],
                                design = ~ group)
  
  # set up df of scaling information
  bbpaths_df = read.csv(bbpaths,
                        header = FALSE,sep="\t",
                        comment.char = "#", 
                        strip.white = TRUE)
  colnames(bbpaths_df)=c("replicate",
                         "sample",
                         "dupstatus",
                         "peaktype",
                         "peakfile",
                         "bedgraph",
                         "scalingfactor")
  sf_df=unique(bbpaths_df[,c("replicate","scalingfactor")])
  dds_cols=colnames(dds)
  sfs=c()
  for (i in dds_cols){
    if (i %in% sf_df$replicate){
      sfs=c(sfs,sf_df[sf_df$replicate==i,"scalingfactor"])
    }
  }
  # scaling factor magnitudes are variable and depend on the constant used while scaling using spiked-in reads
  # DESeq2 size factors are generally hovering around 1
  # we try to rescale the scaling factors by dividing them by mean of all scaling factors ... this way they also 
  # start hovering around 1 ... based on suggestion from Sohyoung.
  if (length(sfs)==length(dds_cols)){
    if (scalesfbymean == "Y") {
      sfs = sfs/mean(sfs)
    }
    
    # AUC-based counts are prescaled, but fragmentbased counts are not prescaled
    if (rawcountsprescaled == "N") {
      rawcounts=round(t(t(rawcounts) * sfs))
      dds <- DESeqDataSetFromMatrix(countData = as.matrix(rawcounts),
                                    colData = sampleinfo[,c("sampleid","group")],
                                    design = ~ group)
    }
    
    DESeq2::sizeFactors(dds)=sfs
  } else {
    print("Samples are spiked, but DESeq2 scaling factors used!!")
  }
  
  dds <- DESeq(dds)
  return(dds)
}

peak_annotation<-function(result_dds,contrast_id){
  #https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2021/ChIPSeq/practicals/ChIP_Practical3_DownstreamAnalysis_2021.html
  #https://support.bioconductor.org/p/103135/
  x = as.data.frame(rownames(result_dds)) 
  colnames(x) = c("peakID")
  x %>% separate(col = c("peakID"),into = c("chrom","coord"),sep = ":") %>% 
    separate(col = c("coord"),into = c("start","end"),sep = "-") -> x
  peaks <- GenomicRanges::makeGRangesFromDataFrame(x)
  options(ChIPseeker.downstreamDistance = 0)
  peakAnno <- ChIPseeker::annotatePeak(peaks,
                                       tssRegion = c(-2000,200),
                                       TxDb = txdb,
                                       level = "gene",
                                       overlap = "all",
                                       annoDb = annodb,
                                       genomicAnnotationPriority = c("Promoter", 
                                                                     "5UTR", "3UTR", "Exon",
                                                                     "Intron","Intergenic"))
  pa <- as.data.frame(peakAnno)
  pa$shortAnno=stringr::word(pa$annotation,1)
  pa$shortAnno[pa$shortAnno=="5'"]="5'UTR"
  pa$shortAnno[pa$shortAnno=="3'"]="3'UTR"
  pa$peakID = paste0(pa$seqnames,":",pa$start,"-",pa$end)
  
  fpath=paste0(output_dir,"peak_annotation_",contrast_id,".csv")
  write.csv(pa,fpath)
  
  return(peakAnno)
}

main_prep_qc<-function(contrast_id,exclusionlist=""){
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
  fpath=paste0(output_dir,"DESEQ2_res_",contrast_id,".csv")
  write.csv(results_df,fpath)
  
  # annotate res
  peakAnno=peak_annotation(result_dds,contrast_id)
  
  return(peakAnno)
}

############################################################
# Create sig counts
############################################################
create_sig_df<-function(contrast_id){
  # read in results
  fpath=paste0(contrast_subpath,contrast_id,extensions,
               "/",contrast_id,extensions,"_",method,"based_diffresults.txt")
  raw_df=read.csv(fpath,sep = "\t")
  
  # filter results for signifcant values
  filt_df=subset(raw_df,padj<=padj_cutoff)
  filt_df=subset(filt_df,(log2FoldChange>=log2fc_cutoff_car) | (log2FoldChange<=-log2fc_cutoff_car))
  if (nrow(filt_df)>0){
    #add metadata
    filt_df$sample=contrast_id
    filt_df$dedup=strsplit(extensions,"__")[[1]][2]
    filt_df$type=strsplit(strsplit(extensions,"__")[[1]][3],"[.]")[[1]][2]
    filt_df$type=filt_df$type %>% replace(is.na(.),"narrowPeak")
    filt_df$method=method
    filt_df$total=nrow(filt_df)
    filt_df$uniqueid=paste0(filt_df$sample,"_",filt_df$dedup,"_",filt_df$type)
    
    print(paste0("----total number of significant peaks for contrast ", contrast_id,": ", nrow(filt_df)))
  } else{
    print(paste0("There are no significant peaks for contrast",contrast_id))
  }
  return(filt_df)
}

############################################################
# Collapse counts
############################################################
create_collapsed_df<-function(merged_sig_df){
  collapsed_df=data.frame()
  
  # for each sample collapse annotation and sort by fold change
  for (sampleid in unique(merged_sig_df$sample)){
    sub_df=subset(merged_sig_df,sample==sampleid)
    
    # collapse to get shortAnno counts
    tmp_collapsed=sub_df %>% dplyr::count(sample,shortAnno,dedup,type,method,total,uniqueid)
    
    # get counts for up/down
    tmp_collapsed$up=0
    tmp_collapsed$down=0
    tmp_direction1=(subset(sub_df,log2FoldChange>0) %>%
                      dplyr::count(sample,shortAnno,dedup,type,method,total,uniqueid))[,c("sample","shortAnno","n")]
    rownames(tmp_direction1)=tmp_direction1$shortAnno
    tmp_direction2=(subset(sub_df,log2FoldChange<0) %>%
                      dplyr::count(sample,shortAnno,dedup,type,method,total,uniqueid))[,c("shortAnno","n")]
    rownames(tmp_direction2)=tmp_direction2$shortAnno
    
    for (rowid2 in rownames(tmp_collapsed)){
      tmp_collapsed[rowid2,"up"]=as.numeric(tmp_direction1[tmp_collapsed[rowid2,"shortAnno"],"n"])
      tmp_collapsed[rowid2,"down"]=as.numeric(tmp_direction2[tmp_collapsed[rowid2,"shortAnno"],"n"])
    }
    tmp_collapsed[is.na(tmp_collapsed)] <- 0
    
    #calculate percentages
    tmp_collapsed$perc=round((tmp_collapsed$n/tmp_collapsed$total)*100,2)
    tmp_collapsed$perc_up=round((tmp_collapsed$up/sum(tmp_collapsed$up))*100,2)
    tmp_collapsed$perc_down=round((tmp_collapsed$down/sum(tmp_collapsed$down))*100,2)
    
    #merge dfs
    collapsed_df=rbind(collapsed_df,tmp_collapsed)
  }
  
  # write out df
  fpath=paste0(output_dir,"collapsed_df.csv")
  write.csv(collapsed_df,fpath)
  
  return(collapsed_df) 
}

############################################################
# Summary graphics
############################################################
chipseeker_plots<-function(peakAnno){
  #http://bioconductor.org/packages/devel/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html
  (plotAnnoPie(peakAnno))
  
  print(upsetplot(peakAnno, vennpie=TRUE))

  (plotDistToTSS(peakAnno,
                      title="Distribution of transcription factor-binding loci\nrelative to TSS"))
}

############################################################
# pie charts for collapsed counts
############################################################
# create pie chart for each comparison
plot_pies_collapsed<-function(sub_in,df_in,y_in,percent_in,plot_in,title_in){
  p = ggplot(sub_in, aes(x = "" , y = get(y_in), fill = fct_inorder(shortAnno))) +
    geom_col(width = 1, color = 1) +
    coord_polar(theta = "y") +
    scale_fill_brewer(palette = "Pastel1") +
    geom_label_repel(data = df_in,
                     aes(y = pos, label = paste0(get(percent_in), "%")),
                     size = 4, nudge_x = 1, show.legend = FALSE) +
    guides(fill = guide_legend(title = "Group")) +
    ggtitle(title_in) +
    theme_void()
  
  p_out=ggarrange(p, 
                  labels = c(plot_in),
                  ncol = 1, nrow = 1)
  return(p_out)
}

# main function
main_piecharts_from_collapsed<-function(sample_id){
  # bring in df
  fpath=paste0(output_car_dir,"collapsed_df.csv")
  collapsed_df=read.csv(fpath,sep=",")
  
  # subset for sample
  sub_df=subset(collapsed_df,sample==sample_id)
  
  # get positions, plot
  ## all totals
  tmp_df <- sub_df %>% 
    mutate(csum = rev(cumsum(rev(n))), 
           pos = n/2 + lead(csum, 1),
           pos = if_else(is.na(pos), n/2, pos))
  plot_title=paste0("Significant Peaks by Annotation (ALL):\n",
                    unique(sub_df$sample)," (N=",sum(sub_df$n),")")
  p1 = plot_pies_collapsed(sub_df,tmp_df,"n","perc","A",plot_title)
  
  ## up
  sampleL=strsplit(unique(sub_df$sample),"_vs_")[[1]][1]
  tmp_df <- sub_df %>% 
    mutate(csum = rev(cumsum(rev(up))), 
           pos = up/2 + lead(csum, 1),
           pos = if_else(is.na(pos), up/2, pos))
  plot_title=paste0(unique(sub_df$sample),"\n",
                    "Significant Peaks by Annotation (N=",sum(sub_df$up),")\n",
                    "Increased in ", sampleL)
  p2 = plot_pies_collapsed(sub_df,tmp_df,"up","perc_up","B",plot_title)
  
  ##down
  tmp_df <- sub_df %>% 
    mutate(csum = rev(cumsum(rev(down))), 
           pos = down/2 + lead(csum, 1),
           pos = if_else(is.na(pos), down/2, pos))
  plot_title=paste0(unique(sub_df$sample),"\n",
                    "Significant Peaks by Annotation (N=",sum(sub_df$down),")\n",
                    "Decreased in ", sampleL)
  p3 = plot_pies_collapsed(sub_df,tmp_df,"down","perc_down","C",plot_title)
  
  print(p1)
  print(p2)
  print(p3)
  
  #create formatted table
  out_df=sub_df[,c("sample","shortAnno","dedup","type","method","n","up","down","total")]
  colnames(out_df)=c("sample_id","annotation","dedup_type","peak_type","norm_type",
                     "sig_peaks","sig_peaks_up","sig_peaks_down","total_peaks")
  DT::datatable(out_df)
}
############################################################
# volcano plots for gene lists
############################################################
generate_volcano_plots<-function(contrast_id,gene_list_in="",gene_title_in=""){
  # read in res from DEG merge
  fpath=paste0(output_car_dir,"DESeq2_res_",contrast_id,".csv")
  res1=read.csv(fpath,sep=",")
  colnames(res1)=c("peakID","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")
  
  # read in PEAK ANNO df from DEG merge
  fpath=paste0(output_car_dir,"peak_annotation_",contrast_id,".csv")
  pa=read.csv(fpath)
  
  # merge annotations and peaks
  results_df = full_join(res1,pa,by=c("peakID"))
  
  # calc log10
  results_df$log_pval=-log10(results_df$pvalue)
  
  # if the gene list is used, subset df to only the master gene list
  if (length(gene_list_in)>1){
    results_df=subset(results_df,SYMBOL %in% gene_list_in)
  }
  
  # add colors and annotate
  colors=brewer.pal(7,"Set1")
  anno_types=levels(as.factor(results_df$shortAnno))
  
  # set all vals to NS and grey
  keyvals=rep("grey",times=nrow(results_df))
  names(keyvals)=rep("NS",times=length(keyvals))
  
  #iterate through each annotation type; change color if
  for ( i in seq(1,length(anno_types))) {
    keyvals[ abs(results_df$log2FoldChange) > log2fc_cutoff_car & 
               results_df$pvalue < padj_cutoff & 
               results_df$shortAnno == anno_types[i]] = colors[i]
    names(keyvals)[keyvals == colors[i]] <- anno_types[i]
  }
  
  ## shapes
  #keyvals.shape <- ifelse(results_df$SYMBOL %in% gene_list_select, 17,3)
  #keyvals.shape[is.na(keyvals.shape)] <- 3
  #names(keyvals.shape)[keyvals.shape == 3] <- 'Not in gene list'
  #names(keyvals.shape)[keyvals.shape == 17] <- 'In gene list'
  
  # print volcano without genelist designation
  p = EnhancedVolcano(results_df,
                      lab = results_df$SYMBOL,
                      x = 'log2FoldChange',
                      y = 'pvalue',
                      ylab = bquote(~-Log[10] ~ FDR),
                      pCutoff = padj_cutoff,
                      FCcutoff = 0,
                      labSize = 4,
                      title = paste0(contrast_id,"\n", gene_title_in),
                      subtitle = "",
                      subtitleLabSize = 1,
                      captionLabSize = 10,
                      colCustom = keyvals,
                      colAlpha = 1,
                      legendLabels = c("NS", expression(Log[2] ~ FC), "FDR", expression(FDR ~ and ~ log[2] ~ FC)),
                      legendLabSize = 10,legendPosition = 'right')
  print(p)
  
  # set labels
  log_pval=-log10(results_df$pvalue)
  y_title="-Log10 FDR"
  stat_type="pvalue"
  log_FC=results_df$log2FoldChange
  Significant=rep("1_NotSignificant",length(log_FC))
  Significant[which(results_df$pvalue<=padj_cutoff & abs(results_df$log2FoldChange)>=log2fc_cutoff_car)]=paste0("3_LogFC_and_",stat_type)
  Significant[which(results_df$pvalue<=padj_cutoff & abs(results_df$log2FoldChange)<=log2fc_cutoff_car)]=paste0("2b_",stat_type,"_Only")
  Significant[which(results_df$pvalue>=padj_cutoff & abs(results_df$log2FoldChange)>=log2fc_cutoff_car)]="2a_LogFC_Only"
  gene=results_df$SYMBOL
  volcano_data=as.data.frame(cbind(gene,log_FC,log_pval,Significant))
  
  p <- plot_ly(data = volcano_data, 
               x = log_FC, y = log_pval, 
               text = gene,
               mode = "markers", 
               color = Significant) %>% 
    layout(title=paste0(contrast_id,"\n", gene_title_in),
           xaxis=list(title="Fold Change",
                      range =c(-5,5),
                      tickvals=c(-5,-4,-3,-2,-1,0,1,2,3,4,5),
                      ticktext=c('-32','-16','-8','-4','-2',
                                 '1','2','4','8','16','32')),
           yaxis=list(title=y_title,range =c(0,15)))
  return(p)
}
