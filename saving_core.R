############################################################
# Create sig counts
############################################################
# creates sig df for peaks / contrast
create_sig_df<-function(contrast_id,norm_type_cutandrun){
  # set extension
  extensions=c(paste0("__",dedup_status,"__",norm_type_cutandrun,".bed"))
  
  # read in results
  fpath=paste0(car_subpath,contrast_id,extensions,
               "/",contrast_id,extensions,"_",method,"based_diffresults.txt")
  raw_df=read.csv(fpath,sep = "\t")
  
  # perform manual change of annotations due to viewing on tracks
  ## only for project CS029758 with sample 53_H4K20me3_IFNb_vs_HN6_H4K20me3_IFNb
  if (cs_id=="CS029758" && contrast_id=="53_H4K20me3_IFNb_vs_HN6_H4K20me3_IFNb"){
    raw_df = fix_annotations_main(raw_df)
  }
  
  fpath=paste0(output_car_dir,"DESEQ2_res_",contrast_id,".csv")
  write.csv(raw_df,fpath)
  
  # create a sig_tracking df
  raw_df$flag_padj <- ifelse(raw_df$padj <=padj_cutoff,"Y", "N") # pval
  raw_df$flag_log2fc <- ifelse(abs(raw_df$log2FoldChange) >=log2fc_cutoff_car,"Y", "N") #log2fc
  raw_df$flag_anno <- ifelse(is.na(raw_df$SYMBOL) | raw_df$shortAnno=="Downstream","N", "Y") #SYMBOL avail & not downstream
  raw_df$flag_anno[raw_df$shortAnno=="Distal"]="N"
  raw_df$flag_padj_log2fc <- ifelse((raw_df$flag_padj == "Y") & (raw_df$flag_log2fc == "Y"),"Y", "N") #pval and log2fc
  raw_df$flag_padj_log2fc_anno <- ifelse((raw_df$flag_padj_log2fc == "Y") & (raw_df$flag_anno == "Y"),"Y", "N") #pval log2fc anno
  fpath=paste0(output_car_dir,"sig_tracking_",contrast_id,".csv")
  write.csv(raw_df,fpath)
  
  # filter results for significant values
  filt_df=subset(raw_df,flag_padj_log2fc_anno=="Y")
  
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

# creates nonsig df for peaks / conrtrast
create_nonsig_df<-function(contrast_id,norm_type_cutandrun){
  # set extension
  extensions=c(paste0("__",dedup_status,"__",norm_type_cutandrun,".bed"))
  
  # read file
  filt_df=read.csv(paste0(output_car_dir,"sig_tracking_",contrast_id,".csv"))
  
  #add metadata
  filt_df$sample=contrast_id
  filt_df$dedup=strsplit(extensions,"__")[[1]][2]
  filt_df$type=strsplit(strsplit(extensions,"__")[[1]][3],"[.]")[[1]][2]
  filt_df$type=filt_df$type %>% replace(is.na(.),"narrowPeak")
  filt_df$method=method
  filt_df$total=nrow(filt_df)
  filt_df$uniqueid=paste0(filt_df$sample,"_",filt_df$dedup,"_",filt_df$type)
  
  return(filt_df)
}

# creates gene df of top n_in significant genes per sample
create_sig_gene_df<-function(cntrl_in,treat_in,n_in){
  #set contrasat
  contras=c(treat_in,cntrl_in)
  
  #read in path
  source_path=paste0(input_dir,"DEG_",contras[1],"-",contras[2],"_0.5_0.5/")
  res1=read.csv(paste0(output_rna_dir,"DESeq2_",contras[1],"-", contras[2],"_DEG_allgenes_res1.txt"),sep="\t")
  
  # subset for sig genes
  sub_res1=subset(res1,(log2FoldChange>=log2fc_cutoff) | (log2FoldChange<=-log2fc_cutoff))
  sub_res1=subset(sub_res1,padj<padj_cutoff)
  
  # subset for top X genes
  sub_top_fc=sub_res1[order(sub_res1$padj),][c(1:n_in),]
  
  #select cols, rename
  sub_top_fc=sub_top_fc[,c("X","log2FoldChange","padj")]
  colnames(sub_top_fc)=c("genes",
                         paste0(contras[1],"_vs_",contras[2],"--log2FoldChange"),
                         paste0(contras[1],"_vs_",contras[2],"--padj"))
  
  return(sub_top_fc)
}

# fills in df for other samples
fillin_sig_gene_df<-function(df_in){
  # add rownames, subset
  rownames(df_in)=df_in$genes
  df_in <- subset(df_in, select = -c(genes))
  
  # for each column, search through and fill in any NA's
  for (colid in colnames(df_in)){
    # pull all NA values for this column
    df_in[df_in=='NA'] <- NA
    missing_rows=rownames(df_in[is.na(df_in[,colid]),])
    
    #read in path
    contras=strsplit(strsplit(colid,"--")[[1]],"_vs_")
    source_path=paste0(input_dir,"DEG_",contras[[1]][1],"-",contras[[1]][2],"_0.5_0.5/")
    res1=read.csv(paste0(output_rna_dir,"DESeq2_",
                         contras[[1]][1],"-", contras[[1]][2],"_DEG_allgenes_res1.txt"),sep="\t")
    rownames(res1)=res1$X
    
    # fill in missing information
    for (rowid in missing_rows){
      df_in[rowid,colid]=res1[rowid,"log2FoldChange"]
    }
  }
  
  return(df_in)
}

############################################################
# Collapse counts
############################################################
create_sig_collapsed_df<-function(merged_df){
  collapsed_df=data.frame()
  
  merged_file_out=paste0(output_car_dir,"merged_sig_",cs_id,".csv")
  collapsed_file_out="collapsed_df.csv"
  
  # write out merged file df
  out_df=merged_df[,c("sample","shortAnno","dedup",
                      "type","method","total","uniqueid")]
  write.csv(out_df,merged_file_out)
  
  # for each sample collapse annotation and sort by fold change
  for (sampleid in unique(merged_df$sample)){
    sub_df=subset(merged_df,sample==sampleid)
    
    # collapse to get shortAnno counts
    tmp_collapsed=sub_df %>% dplyr::count(sample,shortAnno,dedup,type,method,total,uniqueid)
    
    # get counts for up/down
    tmp_collapsed$up=0
    tmp_collapsed$down=0
    tmp_direction1=(subset(sub_df,log2FoldChange>=0) %>%
                      dplyr::count(sample,shortAnno,dedup,type,method,total,uniqueid))[,c("sample","shortAnno","n")]
    rownames(tmp_direction1)=tmp_direction1$shortAnno
    tmp_direction2=(subset(sub_df,log2FoldChange<=0) %>%
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
  fpath=paste0(output_car_dir,collapsed_file_out)
  write.csv(collapsed_df,fpath)
  
  return(collapsed_df) 
}

create_nonsig_collapsed_df<-function(merged_df){
  collapsed_df=data.frame()
  
  merged_file_out=paste0(output_car_dir,"merged_nonsig_",cs_id,".csv")
  collapsed_file_out="collapsed_df_nonsig.csv"
  
  # write out merged file df
  out_df=merged_df[,c("sample","shortAnno","dedup",
                      "type","method","total","uniqueid")]
  write.csv(out_df,merged_file_out)
  
  # for each sample collapse annotation and sort by fold change
  for (sampleid in unique(merged_df$sample)){
    sub_df=subset(merged_df,sample==sampleid&shortAnno != "Downstream")
    
    # update totals without downstream
    sub_df$total=nrow(sub_df)
    
    # collapse to get shortAnno counts
    tmp_collapsed=sub_df %>% dplyr::count(sample,shortAnno,dedup,type,method,total,uniqueid)
    
    #calculate percentages
    tmp_collapsed$perc=round((tmp_collapsed$n/tmp_collapsed$total)*100,2)
    
    #merge dfs
    collapsed_df=rbind(collapsed_df,tmp_collapsed)
  }
  
  # write out df
  fpath=paste0(output_car_dir,collapsed_file_out)
  write.csv(collapsed_df,fpath)
  
  return(collapsed_df) 
}

create_collapsed_pi_df<-function(merged_sig_df){
  pi_df=data.frame()
  for (sampleid in unique(merged_sig_df$sample)){
    
    sub_df=subset(merged_sig_df,sample==sampleid)
    
    # filter for PI genes
    tmp_pi1=subset(sub_df,SYMBOL %in% subset(gene_list,Set=="REPRESSORS")$Human)
    if(nrow(tmp_pi1)!=0){tmp_pi1$gene_list="REPRESSORS"}
    
    tmp_pi2=subset(sub_df,SYMBOL %in% subset(gene_list,Set=="ACCELERATORS")$Human)
    if(nrow(tmp_pi2)!=0){tmp_pi2$gene_list="ACCELERATORS"}
    
    tmp_pi3=subset(sub_df,SYMBOL %in% subset(gene_list,Set=="INVASION")$Human)
    if(nrow(tmp_pi3)!=0){tmp_pi3$gene_list="INVASION"}
    
    if(nrow(tmp_pi1)!=0){pi_df=rbind(pi_df,tmp_pi1)}
    if(nrow(tmp_pi2)!=0){pi_df=rbind(pi_df,tmp_pi2)}
    if(nrow(tmp_pi3)!=0){pi_df=rbind(pi_df,tmp_pi3)}
  }
  
  # write out df
  fpath=paste0(output_car_dir,"collapsed_pi_df.csv")
  write.csv(pi_df,fpath)
  
  return(pi_df)
}

############################################################
# volcano plots for gene lists
############################################################
generate_volcano_plotly_rnaseq<-function(cntrl_in,treat_in,type_in){
  # read in res from DEG merge
  contras=c(treat_in,cntrl_in)
  fpath=paste0(output_rna_dir,"DESeq2_",contras[1],"-", contras[2],"_DEG_allgenes_res1.txt")
  res1=read.csv(fpath,sep="\t")
  
  # Volcano Plots
  if (type_in=="pvalue"){
    log_pval=-log10(res1$pvalue)
    y_title="-Log10 pvalue"
  } else{
    ## logfc and FDR
    log_pval=-log10(res1$padj)
    y_title="-Log10 FDR"
  }
  log_FC=res1$log2FoldChange
  Significant=rep("1_NotSignificant",length(log_FC))
  Significant[which(res1$pvalue<padj_cutoff & abs(res1$log2FoldChange)>=log2fc_cutoff)]=paste0("3_LogFC_and_",type_in)
  Significant[which(res1$pvalue<padj_cutoff & abs(res1$log2FoldChange)<log2fc_cutoff)]=paste0("2b_",type_in,"_Only")
  Significant[which(res1$pvalue>=padj_cutoff & abs(res1$log2FoldChange)>=log2fc_cutoff)]="2a_LogFC_Only"
  gene=res1$X
  volcano_data=as.data.frame(cbind(gene,log_FC,log_pval,Significant))
  p <- plot_ly(data = volcano_data, x = log_FC, y = log_pval, text = gene,
               mode = "markers", 
               color = Significant) %>% layout(title =paste0(contras[1]," vs. ", contras[2]),
                                               xaxis=list(title="Fold Change",
                                                          range =c(-5,5),
                                                          tickvals=c(-5,-4,-3,-2,-1,0,1,2,3,4,5),
                                                          ticktext=c('-32','-16','-8','-4','-2',
                                                                     '1','2','4','8','16','32')),
                                               yaxis=list(title=y_title,range =c(0,15)))
  return(p)
}

############################################################
# find differential genes
############################################################
main_differential_genes<-function(contrast_id,contrast_list){
  fpath=paste0(output_car_dir,"sig_tracking_",contrast_id,".csv")
  car_df=read.csv(fpath)
  
  # subset for sig
  car_df=subset(car_df,flag_padj_log2fc_anno=="Y")
  
  # subset for gene_body_list
  car_df_sub=subset(car_df,shortAnno %in% contrast_list)
  
  # print outputs
  print(paste0("The number of peaks annotated: ",nrow(car_df)))
  print(paste0("The number of genes annotated: ",length(unique(car_df$SYMBOL))))
  print(paste0("The number of peaks annotated in the gene_bodies_list: ",
               nrow(car_df_sub)))
  print(paste0("The number of genes annotated in the gene_bodies_list: ",
               length(unique(car_df_sub$SYMBOL))))
  
  fpath=paste0(output_dir,"diff_genes_",contrast_id,".csv")
  print(fpath)
  write.table(car_df_sub,fpath,sep=",",row.names=FALSE)
  
  fpath=paste0(output_dir,"deeptools_sig_",contrast_id,".csv")
  print(fpath)
  write.csv(unique(car_df_sub$SYMBOL),fpath,row.names=FALSE)
  
}

deg_group_add<-function(cntrl_in,treat_in){
  # run DESE2
  res_comp<-deg_comparison(cntrl_in,treat_in)
  
  # rename cols with sampleid
  tmp_deg=as.data.frame(res_comp)
  colnames(tmp_deg)=paste0(treat_in,"_",colnames(tmp_deg))
}

deg_comparison<-function(cntrl_in,treat_in){
  contras=c(treat_in,cntrl_in)
  
  # subset sample info for contrasts
  full_contrast=paste0(cntrl_in,"-",treat_in)
  sub_df=subset(groups_df,group %in% contras)
  sampleinfo= data.frame(sub_df$group,row.names = sub_df$sampleid)
  colnames(sampleinfo)=c("condition")
  sampleinfo$condition=as.factor(sampleinfo$condition)
  
  # read in Raw count file
  # example: DEG_KO-CRISPR_53_without_IFNb_0.5_0.5/RawCountFile_RSEM_genes_filtered.txt
  fpath=paste0(input_dir,"DEG_",treat_in,"-",
               cntrl_in,"_0.5_0.5/RawCountFile_RSEM_genes_filtered.txt")
  x=read.csv(fpath,sep="\t")
  rownames(x)=x$symbol
  x=x[,c(2:ncol(x))]
  
  # run DESEQ
  ddsHTSeq<-DESeqDataSetFromMatrix(countData=x,colData=sampleinfo, design=~condition)
  dds<-DESeq(ddsHTSeq)
  
  # prep df
  mfc=c()
  mpval=c()
  i=1
  
  # pull results for conditions
  res<-results(dds,contrast=c("condition",as.character(contras[1]),as.character(contras[2])))
  res1=as.data.frame(res)
  write.table(res1,
              file=paste(output_rna_dir,"DESeq2_",contras[1],"-",
                         contras[2],"_DEG_allgenes_res1.txt",
                         sep=""),
              sep="\t",col.names=NA) 
  restmp=res1
  
  # calc fc, create inital df
  restmp$FoldChange <- ifelse(restmp$log2FoldChange<0, -1/(2^restmp$log2FoldChange), 2^restmp$log2FoldChange)
  mfc=cbind(mfc,restmp$FoldChange)
  mpval=cbind(mpval,restmp$pvalue)
  
  x=rownames(restmp)
  ensID=apply(array(as.character(x)),1,function(z) unlist(strsplit(z, "\\|"))[1])
  gene=apply(array(as.character(x)),1,function(z) unlist(strsplit(z, "\\|"))[2])
  restmp=cbind(ensID,gene,restmp)
  
  #remove rownames for df merging downstream
  deseq2out=restmp
  deseq2out$X=rownames(deseq2out)
  
  #subselect df, recalc new values
  deseq2out=deseq2out[,which(names(deseq2out) %in% c("X", "gene","log2FoldChange","pvalue"))]
  deseq2out$fc=2^deseq2out$log2FoldChange
  down_reg=deseq2out$log2FoldChange<0
  deseq2out$fc[down_reg]=-1/deseq2out$fc[down_reg]
  
  # pull final values
  deseq2out=deseq2out[,c("X","gene","fc","log2FoldChange","pvalue")]
  colnames(deseq2out)=c("ensid_gene","gene","fc","log2fc","pvalue")
  deseq2out$fdr=p.adjust(deseq2out$pvalue,method='fdr',n=length(deseq2out$pvalue))
  deseq2out$gsea_ranking_score=-log10(deseq2out$pvalue)*sign(deseq2out$log2fc)
  write.table(deseq2out,file=paste(output_rna_dir,"DESeq2_",contras[1],"-",
                                   contras[2],"_DEG_allgenes.txt",sep=""),
              row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")
  return(deseq2out)
}

################################################################################
# heatmaps for deg
################################################################################
create_two_sample_df<-function(contrast_id,extensions,sig_flag){
  # read in results
  fpath=paste0(car_subpath,contrast_id,extensions,
               "/",contrast_id,extensions,"_",method,"based_diffresults.txt")
  car_df_filt=read.csv(fpath,sep = "\t")
  
  # remove missing values
  car_df_filt=car_df_filt[complete.cases(car_df_filt),]
  
  # split pos, neg
  df_in_pos=subset(car_df_filt,log2FoldChange>0)
  df_in_neg=subset(car_df_filt,log2FoldChange<0)
  
  # remove duplicate IDS
  df_in_pos=df_in_pos[!duplicated(df_in_pos$ENSEMBL),]
  df_in_neg=df_in_neg[!duplicated(df_in_neg$ENSEMBL),]
  
  # merge
  df_out=rbind(df_in_pos,
               df_in_neg)
  
  # if sigflag on, subset for only sigs
  if (sig_flag=="ON"){
    df_out=subset(df_out,(padj<padj_cutoff) & (abs(log2FoldChange)>=log2fc_cutoff_car))
  }
  
  # sort and remove duplicates
  df_out=df_out[order(df_out$log2FoldChange),]
  df_out=df_out[!duplicated(df_out$ENSEMBL),]
  
  # set rownames
  rownames(df_out)=df_out$ENSEMBL
  return(df_out)
}

# creates heatmap
generate_heat_map_two_samples<-function(contrast_id1,contrast_id2,sig_flag="Off"){
  ####################
  # process df
  #####################
  # split df, remove duplicates
  car_df_filt=create_two_sample_df(contrast_id1,extensions1,sig_flag)
  car_df_filt2=create_two_sample_df(contrast_id2,extensions2,sig_flag)
  
  # merge df
  car_df_merged=car_df_filt[,c("ENSEMBL","log2FoldChange")]
  car_df_merged=merge(car_df_merged,
                      car_df_filt2[,c("ENSEMBL","log2FoldChange")],by="ENSEMBL",all=TRUE)
  colnames(car_df_merged)=c("ENSEMBL",
                            strsplit(gsub("53_","",contrast_id1),"_vs_")[[1]][1],
                            strsplit(gsub("53_","",contrast_id2),"_vs_")[[1]][1])
  
  # set rownames
  rownames(car_df_merged)=car_df_merged$ENSEMBL
  
  # set NA to 0
  car_df_merged[is.na(car_df_merged)]=0
  
  ####################
  # formatting
  #####################
  # Overwrite pheatmaps default draw_colnames with new version
  assignInNamespace(x="draw_colnames", value="draw_colnames_45",ns=asNamespace("pheatmap")) 
  
  # Heatmap Color Gradients 
  paletteLength <- 1000
  mycolors <- colorRampPalette(c("blue","white","red"), interpolate = "linear")(paletteLength)
  
  ####################
  # plot
  #####################
  if (sig_flag=="OFF"){
    plot_title="All Genes Included"
  } else {
    plot_title="Significant Genes Included"
  }
  heatmap.2(as.matrix(car_df_merged[,c(2:3)]), 
            labRow = FALSE, main=plot_title, ylab="ENSEMBL ID",
            scale="none",cexRow = 0.5, cexCol = 0.8,trace="none")
}

# intensity hatmap
generate_intensity_maps<-function(contrast_list,fpath1,fpath2,status_flag="single"){
  
  # set promoter
  promoter <- getPromoters(TxDb=txdb_var, upstream=3000, downstream=3000)
  
  
  if (status_flag=="dual"){
    # pull plot lists
    taxMatrix_list=list(getTagMatrix(readPeakFile(fpath1), windows=promoter),
                        getTagMatrix(readPeakFile(fpath2), windows=promoter))
    
    # set names
    names(taxMatrix_list)=c(contrast_list[1],contrast_list[2])
    
    #plot
    tagHeatmap(taxMatrix_list, xlim=c(-3000, 3000), color=NULL)
    plotAvgProf(taxMatrix_list[1], xlim=c(-3000, 3000),
                xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
    plotAvgProf(taxMatrix_list[2], xlim=c(-3000, 3000),
                xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
  } else{
    # pull data
    taxMatrix_list=getTagMatrix(readPeakFile(fpath1), windows=promoter)
    
    #plot
    tagHeatmap(taxMatrix_list, xlim=c(-3000, 3000), color=NULL)
    plotAvgProf(taxMatrix_list, xlim=c(-3000, 3000),
                xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
  }
}

# Overwrites the pheatmap defaults
draw_colnames_45 <- function (coln, gaps, ...) {
  "Overwrites body of pheatmap:::draw_colnames, customizing it my liking"
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 45, gp = gpar(...))
  return(res)
}

# save heatmap
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

# creates heatmap
generate_heat_map_differential<-function(contrast_id,n_up,n_down,fshort,gene_list_flag="OFF"){
  ####################
  # process df
  #####################
  # cut id
  cut_id=strsplit(contrast_id,"_vs_")[[1]][1]
  
  # read in data
  fpath=paste0(output_dir,"diff_genes_",contrast_id,".csv")
  car_df=read.table(fpath,sep=",",header=TRUE)
  
  car_df_filt=subset(car_df,shortAnno %in% gene_bodies_list)
  rownames(car_df_filt)<-NULL
  
  #sort 
  car_df_filt=car_df_filt[order(car_df_filt$log2FoldChange,decreasing=TRUE),]
  ####################
  # formatting
  #####################
  # Overwrite pheatmaps default draw_colnames with new version
  assignInNamespace(x="draw_colnames", value="draw_colnames_45",ns=asNamespace("pheatmap")) 
  
  # Heatmap Color Gradients 
  # split above and below threshold because otherwise it wont center on zero
  paletteLength <- 1000
  nHalf=nrow(car_df_filt)/2
  Min = min(car_df_filt$log2FoldChange)
  Max = max(car_df_filt$log2FoldChange)
  Thresh = 0
  
  rc1 = colorRampPalette(colors = c("blue", "white"), space="Lab")(nHalf)   
  rc2 = colorRampPalette(colors = c("white", "red"), space="Lab")(nHalf)
  rampcols = c(rc1, rc2)
  
  rampcols[c(nHalf, nHalf+1)] = rgb(t(col2rgb("white")), maxColorValue=256) 
  rb1 = seq(Min, Thresh, length.out=nHalf+1)
  rb2 = seq(Thresh, Max, length.out=nHalf+1)[-1]
  rampbreaks = c(rb1, rb2)
  
  mycolors <- colorRampPalette(c("blue","white","red"), interpolate = "linear")(paletteLength)
  
  ####################
  # plot "all" heatmap
  #####################
  p1 = pheatmap(car_df_filt[,c("log2FoldChange")], cluster_cols=F, cluster_rows = F,
                scale = "none", main=paste0("Significantly differentiated ",cut_id,
                                            " peaks\nassociated with gene annotated intragenic regions"),
                cellwidth = 30, fontsize = 12, fontsize_row = 4, fontsize_col = 8, color = rampcols,
                breaks=rampbreaks,border_color = "NA")
  
  fpath=paste0(img_dir,"heatmap_all_",cut_id,".pdf")
  save_pheatmap_pdf(p1, fpath)
  ####################
  # prep "subset" df
  #####################
  # sort
  car_df_sub=car_df_filt[order(car_df_filt$log2FoldChange),]
  
  # take top/bottom N genes
  top_df=rbind(car_df_sub[c(1:n_down),],
               car_df_sub[c((nrow(car_df_sub)-(n_up-1)):nrow(car_df_sub)),])
  
  # create a newID
  top_df$newID=paste0(top_df$SYMBOL," (",top_df$shortAnno,")")
  
  #remove duplicated id's 
  car_df_type=(top_df[!duplicated(top_df$newID),])[,c("newID","log2FoldChange")]
  
  # print values to help with subsetting
  print(paste0("There are up:",nrow(subset(car_df_type,log2FoldChange>=0))))
  print(paste0("There are down:",nrow(subset(car_df_type,log2FoldChange<=0))))
  
  # add rownames, set na to 0
  rownames(car_df_type)=car_df_type$newID
  car_df_type[is.na(car_df_type)]=0
  
  # sort
  car_df_type=car_df_type[order(car_df_type$log2FoldChange,decreasing=TRUE),]
  head(car_df_type)
  
  ####################
  # plot
  #####################
  nHalf=nrow(car_df_type)/2
  Min = min(car_df_type$log2FoldChange)
  Max = max(car_df_type$log2FoldChange)
  Thresh = 0
  
  rc1 = colorRampPalette(colors = c("blue", "white"), space="Lab")(nHalf)   
  rc2 = colorRampPalette(colors = c("white", "red"), space="Lab")(nHalf)
  rampcols = c(rc1, rc2)
  
  rampcols[c(nHalf, nHalf+1)] = rgb(t(col2rgb("white")), maxColorValue=256) 
  rb1 = seq(Min, Thresh, length.out=nHalf+1)
  rb2 = seq(Thresh, Max, length.out=nHalf+1)[-1]
  rampbreaks = c(rb1, rb2)
  
  mycolors <- colorRampPalette(c("blue","white","red"), interpolate = "linear")(paletteLength)
  
  p1=pheatmap(car_df_type[,"log2FoldChange"], cluster_rows=F,labels_row=car_df_type$newID,
              scale = "none", main=paste0("Select ",length(unique(car_df_type$newID)),
                                          " significantly differentiated ",cut_id,
                                          " peaks\nassociated with gene annotated intragenic regions"),
              cellwidth = 30, fontsize = 12, fontsize_row = 6, fontsize_col = 8, color = rampcols,
              breaks=rampbreaks,border_color = "NA",cluster_cols=F)
  
  fpath=paste0(img_dir,"heatmap_select_",cut_id,".pdf")
  save_pheatmap_pdf(p1, fpath)
  
  ####################
  # gene list 
  #####################
  # subset
  if (gene_list_flag=="ON"){
    # subset for only genes in list
    car_df_sub=subset(car_df_filt,SYMBOL%in%gene_list)
    
    # create a newID
    car_df_sub$newID=paste0(car_df_sub$SYMBOL," (",car_df_sub$shortAnno,")")
    
    ####################
    # formatting
    #####################
    # Heatmap Color Gradients 
    # split above and below threshold because otherwise it wont center on zero
    paletteLength <- 1000
    nHalf=nrow(car_df_sub)/2
    Min = min(car_df_sub$log2FoldChange)
    Max = max(car_df_sub$log2FoldChange)
    Thresh = 0
    
    rc1 = colorRampPalette(colors = c("blue", "white"), space="Lab")(nHalf)   
    rc2 = colorRampPalette(colors = c("white", "red"), space="Lab")(nHalf)
    rampcols = c(rc1, rc2)
    
    rampcols[c(nHalf, nHalf+1)] = rgb(t(col2rgb("white")), maxColorValue=256) 
    rb1 = seq(Min, Thresh, length.out=nHalf+1)
    rb2 = seq(Thresh, Max, length.out=nHalf+1)[-1]
    rampbreaks = c(rb1, rb2)
    
    mycolors <- colorRampPalette(c("blue","white","red"), interpolate = "linear")(paletteLength)
    
    ####################
    # plot heatmap
    #####################
    p1 = pheatmap(car_df_sub[,c("log2FoldChange")], cluster_rows=F,cluster_cols=F, labels_row=car_df_sub$newID,
                  scale = "none", main=paste0("Significantly differentiated ",cut_id,
                                              " peaks\nassociated with immune gene annotated intragenic regions"),
                  cellwidth = 30, fontsize = 12, fontsize_row = 6, fontsize_col = 8, color = rampcols,
                  breaks=rampbreaks,border_color = "NA")
    
    fpath=paste0(img_dir,"heatmap_immune_",cut_id,".pdf")
    save_pheatmap_pdf(p1, fpath)
  }
  
}

################################################################################
# boxplots
################################################################################
genebodies_boxplots<-function(contrast_id,extensions,sig_flag="Off"){
  fpath=paste0(car_subpath,contrast_id,extensions,
               "/",contrast_id,extensions,"_",method,"based_diffresults.txt")
  car_df_filt=read.csv(fpath,sep = "\t")
  
  # if sigflag on, subset for only sigs
  if (sig_flag=="ON"){
    car_df_filt=subset(car_df_filt,(padj<padj_cutoff) & (abs(log2FoldChange)>=log2fc_cutoff_car))
  }
  
  car_df_tmp=subset(car_df_filt,shortAnno=="Promoter")
  car_df_tmp$anno_type="Promoter"
  car_df_select=car_df_tmp[,c("log2FoldChange","anno_type")]
  
  car_df_tmp=subset(car_df_filt,shortAnno=="Intron")
  car_df_tmp$anno_type="Intron"
  car_df_select=rbind(car_df_select,car_df_tmp[,c("log2FoldChange","anno_type")])
  
  car_df_tmp=subset(car_df_filt,shortAnno=="Exon")
  car_df_tmp$anno_type="Exon"
  car_df_select=rbind(car_df_select,car_df_tmp[,c("log2FoldChange","anno_type")])
  
  boxplot(log2FoldChange~anno_type,data=car_df_select, main=contrast_id,
          xlab="Annotation Type", ylab="Log2FoldChange")
}


############################################################
# find differential overlap between CAR and RNA
############################################################

create_overlap_chrommap<-function(merged_df,subset_type){
  #http://bioconductor.org/packages/release/bioc/vignettes/karyoploteR/inst/doc/karyoploteR.html
  # subset and remove NA's
  #sub_df=merged_df%>%tidyr::separate(peakID,c("seqnames","peakID"),sep=":")
  sub_df=subset(merged_df,overlap_type=="overlap")
  sub_df=subset(sub_df,seqnames %in% paste0("chr",c(1:22,"X","Y")))
  sub_df=sub_df%>% rename(c("Start"=geneStart,"End"=geneEnd))
  
  # create lists of up and down regulated
  both_up_df=subset(sub_df,rna_status=="up" & car_status=="up")[,c("seqnames","Start","End")]
  both_down_df=subset(sub_df,rna_status=="down" & car_status=="down")[,c("seqnames","Start","End")]
  car_up_df= subset(sub_df,rna_status=="down" & car_status=="up")[,c("seqnames","Start","End")]
  car_down_df=subset(sub_df,rna_status=="up" & car_status=="down")[,c("seqnames","Start","End")]
  
  # total sig genes
  total_peaks=nrow(both_up_df) + nrow(both_down_df) + nrow(car_up_df) + nrow(car_down_df)
  total_genes=length(unique(sub_df$SYMBOL))
  
  # plot
  par(mfrow = c(1,1))
  kp=plotKaryotype(genome=genome)
  if (nrow(both_up_df)>0){ kpPlotRegions(kp, GRanges(both_up_df), col="#FFBE33")}
  if (nrow(both_down_df)>0){ kpPlotRegions(kp, GRanges(both_down_df), col="#5BFF33")}
  if (nrow(car_up_df)>0){ kpPlotRegions(kp, GRanges(car_up_df), col="#337DFF")}
  if (nrow(car_down_df)>0){ kpPlotRegions(kp, GRanges(car_down_df), col="#FFAACC")}
  legend(x="bottomright",
         legend=(c("Up_both","Down_both","Up_CAR_only","Down_CAR_only")),
         fill = c("#FFBE33","#5BFF33","#337DFF","#FFAACC"))
  mtext(paste0("Karyoplot: ", fshort, " sig genes (", total_genes,
               ") and peaks (", total_peaks,")."),
        line=3)
}

create_overlap_DT<-function(df_in){
  df_in=df_in[,c("overlap_type","peakID","annotation","ENSEMBL","SYMBOL",
                 "log2FC_car","padj_car","log2FC_rna","padj_rna")]
  
  # round sigfigs
  col_list=c("log2FC_car","padj_car","log2FC_rna","padj_rna")
  for (colid in col_list){
    df_in[,colid]=signif(df_in[,colid], digits=3)
  }
  
  DT::datatable(df_in)
}

main_differential_overlap<-function(subset_list,rna_regulation,car_regulation){
  #peak db
  fpath=paste0(output_car_dir,"peak_annotation_",contrast_id_car,".csv")
  peak_df=read.csv(fpath)
  
  # annotation db
  fpath=paste0(output_car_dir,"sig_tracking_",contrast_id_car,".csv")
  deseq_df=read.csv(fpath)
  deseq_df=deseq_df[,c(2:ncol(deseq_df))]
  
  # remove any row without a gene symbol
  car_df=deseq_df[complete.cases(deseq_df$SYMBOL),]
  car_df$source_CAR="Y"
  
  # rna db
  fpath=paste0(output_rna_dir,"DESeq2_",contrast_id_rna,"_DEG_allgenes_res1.txt")
  rna_df=read.csv(fpath,sep="\t")
  rna_df=separate(rna_df,"X",c("ENSEMBL","SYMBOL"),sep="[|]")
  
  # separate ensembl for merging
  rna_df=tidyr::separate(rna_df,ENSEMBL,c("ENSEMBL","ID"),sep="[.]")
  rna_df$source_RNA="Y"
  
  # clean ENSEMBL ID's before merging
  missing_eid=rna_df[is.na(rna_df$ENSEMBL),]
  clean_eid=car_df[!is.na(car_df$ENSEMBL),]
  for (rowid in rownames(missing_eid)){
    lookup_symbol=missing_eid[rowid,"SYMBOL"]
    new_eid=subset(clean_eid,SYMBOL==lookup_symbol)$ENSEMBL[[1]]
    if(!is.na(new_eid)){
      rna_df[rowid,"ENSEMBL"]=new_eid
    }
  }
  
  missing_eid=car_df[is.na(car_df$ENSEMBL),]
  nrow(missing_eid)
  clean_eid=rna_df[!is.na(rna_df$ENSEMBL),]
  for (rowid in rownames(missing_eid)){
    lookup_symbol=missing_eid[rowid,"SYMBOL"]
    new_eid=subset(clean_eid,SYMBOL==lookup_symbol)$ENSEMBL
    if(length(new_eid)>0){
      car_df[rowid,"ENSEMBL"]=new_eid
    }
  }
  
  # save merged unfiltered
  unfiltered_merge=full_join(rna_df %>% rename(c("log2FC_rna"=log2FoldChange,
                                                 "padj_rna"=padj,
                                                 "pval_rna"=pvalue,
                                                 "stat_rna"=stat,
                                                 "baseMean_rna"=baseMean,
                                                 "lfcSE_rna"=lfcSE)),
                             car_df %>% rename(c("log2FC_car"=log2FoldChange,
                                                 "padj_car"=padj,
                                                 "pval_car"=pvalue,
                                                 "stat_car"=stat,
                                                 "baseMean_car"=baseMean,
                                                 "lfcSE_car"=lfcSE)),
                             by=c("ENSEMBL","SYMBOL"))
  
  #rna status
  unfiltered_merge$rna_status="not_sig"
  unfiltered_merge$rna_status[unfiltered_merge$log2FC_rna>=log2fc_cutoff_rna & 
                                unfiltered_merge$padj_rna<padj_cutoff]="up"
  unfiltered_merge$rna_status[unfiltered_merge$log2FC_rna<=-log2fc_cutoff_rna & 
                                unfiltered_merge$padj_rna<padj_cutoff]="down"
  unique(unfiltered_merge$rna_status)
  
  # car status
  unfiltered_merge$car_status="not_sig"
  unfiltered_merge$car_status[unfiltered_merge$log2FC_car>=log2fc_cutoff_car & 
                                unfiltered_merge$flag_padj_log2fc_anno=="Y"]="up"
  unfiltered_merge$car_status[unfiltered_merge$log2FC_car<=-log2fc_cutoff_car & 
                                unfiltered_merge$flag_padj_log2fc_anno=="Y"]="down"
  unique(unfiltered_merge$car_status)
  
  unfiltered_merge=unfiltered_merge[,c("peakID","ENSEMBL","ID","SYMBOL",
                                       "rna_status", "car_status",
                                       "source_RNA","source_CAR",
                                       "baseMean_rna","baseMean_car",
                                       "log2FC_rna","log2FC_car",
                                       "padj_rna","padj_car",
                                       "lfcSE_rna","lfcSE_car",
                                       "stat_rna","stat_car",
                                       "pval_rna","pval_car",
                                       "seqnames","start","end","width","strand",
                                       "annotation","geneChr",
                                       "geneStart","geneEnd","geneLength","geneStrand",
                                       "geneId","distanceToTSS",
                                       "GENENAME","shortAnno")]
  
  fpath=paste0(output_dir,"unfilt_overlap_",contrast_id_car,"_",contrast_id_rna,".csv")
  write.table(unfiltered_merge,fpath,sep=",",row.names = FALSE)
  
  # determine overlap based on directional
  if (rna_regulation=="all" & car_regulation=="all"){
    merged_df=subset(unfiltered_merge,rna_status=="up" | rna_status=="down"| car_status=="up" | car_status=="down")
    merged_df$overlap_type="none"
    
    merged_df$overlap_type[merged_df$source_CAR=="Y" & merged_df$source_RNA=="Y"]="overlap"
    merged_df$overlap_type[merged_df$source_CAR=="Y" & merged_df$source_RNA=="N"]="only_rna"
    merged_df$overlap_type[merged_df$source_CAR=="N" & merged_df$source_RNA=="Y"]="only_car"
  } else{
    merged_df=subset(unfiltered_merge,rna_status==rna_regulation | car_status==car_regulation)
    merged_df$overlap_type="none"
    
    merged_df$overlap_type[merged_df$rna_status==rna_regulation & merged_df$car_status==car_regulation]="overlap"
    merged_df$overlap_type[merged_df$rna_status!=rna_regulation & merged_df$overlap_type!="overlap"]="only_car"
    merged_df$overlap_type[merged_df$car_status!=car_regulation & merged_df$overlap_type!="overlap"]="only_rna"
    unique(merged_df$overlap_type)
  }
  
  # subset for list
  if ("all" %in% subset_list){
    print("Running ALL annotations")
  } else {
    print ("Running SET of annotations")
    merged_df=subset(merged_df,shortAnno %in% subset_list)
  }
  
  print(paste0("Running CAR peaks ", car_regulation, " (", nrow(subset(merged_df,car_status==car_regulation)), ") | ",
               "CAR genes: ", car_regulation, " (",length(unique(subset(merged_df,car_status==car_regulation)$SYMBOL)),") | ",
               "RNA genes ", rna_regulation, " (",nrow(subset(merged_df,rna_status==rna_regulation)),")"))
  
  fpath=paste0(output_dir,"filt_overlap_",fshort, "_",contrast_id_car,"_",contrast_id_rna,".csv")
  write.table(merged_df,fpath,sep=",",row.names = FALSE)
  
  #if there are sig genes, create venn diagrams, chrom map, DT
  if (nrow(merged_df)!=0){
    print("There are significant genes")
    create_venn_diagrams(fshort,merged_df)
    create_overlap_chrommap(merged_df,subset_list)
    create_overlap_DT(subset(merged_df,overlap_type=="overlap"))
  } else{
    print(paste0("There are no signifcant genes overlapping between datasets for these filters:\n",
                 "-",contrast_id_car,"\n",
                 "-",contrast_id_rna,"\n",
                 "-",subset_list, "genes"))
  }
}

create_overlap_heatmap<-function(df_in,contrast_id){
  ####################
  # input
  #####################
  #sort 
  car_df_filt=df_in[order(df_in$log2FC_car,decreasing=TRUE),]
  
  ####################
  # formatting
  #####################
  # Overwrite pheatmaps default draw_colnames with new version
  assignInNamespace(x="draw_colnames", value="draw_colnames_45",ns=asNamespace("pheatmap")) 
  
  # Heatmap Color Gradients 
  # split above and below threshold because otherwise it wont center on zero
  paletteLength <- 1000
  nHalf=nrow(car_df_filt)/2
  Min = min(car_df_filt$log2FC_car)
  Max = max(car_df_filt$log2FC_car)
  Thresh = 0
  
  rc1 = colorRampPalette(colors = c("blue", "white"), space="Lab")(nHalf)   
  rc2 = colorRampPalette(colors = c("white", "red"), space="Lab")(nHalf)
  rampcols = c(rc1, rc2)
  
  rampcols[c(nHalf, nHalf+1)] = rgb(t(col2rgb("white")), maxColorValue=256) 
  rb1 = seq(Min, Thresh, length.out=nHalf+1)
  rb2 = seq(Thresh, Max, length.out=nHalf+1)[-1]
  rampbreaks = c(rb1, rb2)
  
  mycolors <- colorRampPalette(c("blue","white","red"), interpolate = "linear")(paletteLength)
  
  ####################
  # plot "all" heatmap
  #####################
  p1 = pheatmap(car_df_filt[,c("log2FC_car")], ,cluster_cols=F, cluster_rows = F,
                scale = "none", main=paste0("Significantly differentiated ",contrast_id,
                                            " peaks \nwhich induce differential mRNA expression associated\nwith gene annotated intragenic regions"),
                cellwidth = 30, fontsize = 12, fontsize_row = 6, fontsize_col = 8, color = rampcols,
                breaks=rampbreaks,border_color = "NA")
  
  fpath=paste0(img_dir,"heatmap_mrna_overlap_all_",contrast_id,".pdf")
  save_pheatmap_pdf(p1, fpath)
  ####################
  # prep "subset" df
  #####################
  # sort
  car_df_sub=car_df_filt[order(car_df_filt$padj_car),]
  
  # take top/bottom N genes
  car_df_type=data.frame()
  n=50
  while(nrow(car_df_type)!=50){
    print(paste0("Attempting ",n))
    top_df=car_df_sub[c(1:n),]
    
    # create a newID
    top_df$newID=paste0(top_df$SYMBOL," (",top_df$shortAnno,")")
    
    #remove duplicated id's 
    car_df_type=(top_df[!duplicated(top_df$newID),])[,c("newID","log2FC_car")]
    n=n+1
  }
  
  # print values to help with subsetting
  print(paste0("There are up:",nrow(subset(car_df_type,log2FC_car>=0))))
  print(paste0("There are down:",nrow(subset(car_df_type,log2FC_car<=0))))
  
  # add rownames, set na to 0
  rownames(car_df_type)=car_df_type$newID
  car_df_type[is.na(car_df_type)]=0
  
  # sort
  car_df_type=car_df_type[order(car_df_type$log2FC_car,decreasing=TRUE),]
  head(car_df_type)
  
  ####################
  # plot
  #####################
  nHalf=nrow(car_df_type)/2
  Min = min(car_df_type$log2FC_car)
  Max = max(car_df_type$log2FC_car)
  Thresh = 0
  
  rc1 = colorRampPalette(colors = c("blue", "white"), space="Lab")(nHalf)   
  rc2 = colorRampPalette(colors = c("white", "red"), space="Lab")(nHalf)
  rampcols = c(rc1, rc2)
  
  rampcols[c(nHalf, nHalf+1)] = rgb(t(col2rgb("white")), maxColorValue=256) 
  rb1 = seq(Min, Thresh, length.out=nHalf+1)
  rb2 = seq(Thresh, Max, length.out=nHalf+1)[-1]
  rampbreaks = c(rb1, rb2)
  
  mycolors <- colorRampPalette(c("blue","white","red"), interpolate = "linear")(paletteLength)
  
  ####################
  # plot
  #####################
  p1=pheatmap(car_df_type[,"log2FC_car"], cluster_rows=F,labels_row=car_df_type$newID,
              scale = "none", main=paste0("Select ",length(unique(car_df_type$newID)),
                                          " Significantly differentiated ",contrast_id,
                                          " peaks \nwhich induce differential mRNA expression associated\nwith gene annotated intragenic regions"),
              cellwidth = 30, fontsize = 12, fontsize_row = 6, fontsize_col = 8, color = rampcols,
              breaks=rampbreaks,border_color = "NA",cluster_cols=F)
  fpath=paste0(img_dir,"heatmap_mrna_overlap_subset_",contrast_id,".pdf")
  save_pheatmap_pdf(p1, fpath)
}

main_differential_overlap_all<-function(contrast_id_car,contrast_id_rna){
  # read genebodies df up and down, merge
  fpath_down=paste0(output_dir,"filt_overlap_genebodies_down_",contrast_id_car,"_",contrast_id_rna,".csv")
  fpath_up=paste0(output_dir,"filt_overlap_genebodies_up_",contrast_id_car,"_",contrast_id_rna,".csv")
  merged_df=rbind(read.csv(fpath_down),read.csv(fpath_up))
  
  # cut id
  cut_id=strsplit(contrast_id_car,"_vs_")[[1]][1]
  
  # subset for overlap
  overlap_df=subset(merged_df,overlap_type=="overlap")
  
  fpath=paste0(output_dir,"filter_overlap_",contrast_id_car,".csv")
  write.csv(overlap_df,fpath)
  
  fpath=paste0(output_dir,"deeptools_overlap_",contrast_id_car,".csv")
  write.csv(unique(overlap_df$SYMBOL),fpath,row.names=FALSE)
  
  # create heatmap of df
  create_overlap_heatmap(overlap_df,cut_id)
  
}
############################################################
# find differential overlap between CAR and CAR
############################################################
create_venn_diagrams_car_to_car<-function(cs_id1_short,cs_id2_short,merged_df){
  # create gene lists
  list_of_car1_genes=unique(subset(merged_df,overlap_type != cs_id2_short)$SYMBOL)
  list_of_car2_genes=unique(subset(merged_df,overlap_type != cs_id1_short)$SYMBOL)
  
  # remove NA"s
  list_of_car1_genes=list_of_car1_genes[!is.na(list_of_car1_genes)]
  list_of_car2_genes=list_of_car2_genes[!is.na(list_of_car2_genes)]
  
  # List of genes
  x <- list(A = list_of_car1_genes, B = list_of_car2_genes)
  
  # Venn diagram with custom category names
  p = ggVennDiagram(x, color = 1, lwd = 0.7,
                    category.names = c(paste0(cs_id1_short," Genes"),
                                       paste0(cs_id2_short," Genes"))) + 
    scale_fill_gradient(low = "red", high = "blue")
  full_title="Significantly Differentiated Genes"
  pf = p + ggtitle(full_title)
  
  # save and print
  fpath=paste0(output_car_dir,"venndiagram_",cs_id1_short,"_",cs_id2_short,"_genes.png")
  ggsave(fpath,pf)
  
  print(pf)
}

main_differential_overlap_car_to_car<-function(contrast_list,contrast_id_car1,contrast_id_car2,cs_id1,cs_id2){
  # car db1
  fpath=paste0(contrast_id_car1,"sig_tracking_",cs_id1,".csv")
  car_df1=read.csv(fpath)
  
  # car db2
  fpath=paste0(contrast_id_car2,"sig_tracking_",cs_id2,".csv")
  car_df2=read.csv(fpath)
  
  # filter dbs for sig and shortAnno
  cols_to_save=c("SYMBOL","peakID","log2FoldChange","padj","seqnames","ENSEMBL","shortAnno")
  car_df1=subset(car_df1,shortAnno %in% contrast_list & flag_padj_log2fc_anno=="Y")[,cols_to_save]
  car_df2=subset(car_df2,shortAnno %in% contrast_list & flag_padj_log2fc_anno=="Y")[,cols_to_save]
  
  # create short id
  cs_id1_short=strsplit(cs_id1,"_vs_")[[1]][1]
  cs_id2_short=strsplit(cs_id2,"_vs_")[[1]][1]
  
  #rename cols
  colnames(car_df1)=c("SYMBOL",paste0(colnames(car_df1)[c(2:7)],"_",cs_id1_short))
  colnames(car_df2)=c("SYMBOL",paste0(colnames(car_df2)[c(2:7)],"_",cs_id2_short))
  
  # create merged df
  filtered_merge=full_join(car_df1,
                           car_df2,
                           by=c("SYMBOL"))
  
  #determine sigs
  filtered_merge$cs_id1=ifelse(!is.na(filtered_merge[,paste0("shortAnno_",cs_id1_short)]),"Y","N")
  filtered_merge$cs_id2=ifelse(!is.na(filtered_merge[,paste0("shortAnno_",cs_id2_short)]),"Y","N")
  
  # create cats
  filtered_merge$overlap_type[filtered_merge$cs_id1=="Y" & filtered_merge$cs_id2=="Y"]="overlap"
  filtered_merge$overlap_type[filtered_merge$cs_id1=="Y" & filtered_merge$cs_id2=="N"]=cs_id1_short
  filtered_merge$overlap_type[filtered_merge$cs_id1=="N" & filtered_merge$cs_id2=="Y"]=cs_id2_short
  
  unique(filtered_merge$overlap_type)
  head(filtered_merge)
  
  # save df
  fpath=paste0(output_dir,"filt_overlap_",cs_id1_short,"_",cs_id2_short,".csv")
  write.table(filtered_merge,fpath,sep=",",row.names = FALSE)
  
  #if there are sig genes, create venn diagrams, chrom map, DT
  if (nrow(subset(filtered_merge,overlap_type=="overlap"))!=0){
    print("There are significant genes")
    create_venn_diagrams_car_to_car(cs_id1_short,cs_id2_short,filtered_merge)
  } else{
    print(paste0("There are no signifcant genes overlapping between datasets"))
  }
}

################################################################################
# functions ORA/GSEA
################################################################################
capture_entrezids<-function(input_df){
  sep_df=input_df %>%
    separate(ensid_gene,sep="[|]",c("ENSEMBL","SYMBOL"))%>%
    separate(ENSEMBL,sep="[.]",c("ENSEMBL","ID"))
  
  # search for ENTREZ by ENSEMBL
  gene_df_e <- bitr(sep_df[,c("ENSEMBL")], fromType = "ENSEMBL",
                    toType = c("ENTREZID"),
                    OrgDb = org.Hs.eg.db)
  
  # if any genes did not map, try by SYMBOL
  missing_df=subset(sep_df,ENSEMBL %ni% gene_df_e$ENSEMBL)
  if (nrow(missing_df)>0){
    gene_df_s <- bitr(sep_df[,c("SYMBOL")], fromType = "SYMBOL",
                      toType = c("ENTREZID"),
                      OrgDb = org.Hs.eg.db)
    
    # merge dfs together
    gene_output=merge(gene_df_e,gene_df_s,all=TRUE,by="ENTREZID")
    
  }
  
  # remove duplicated ENSEMBL ID's that have a gene symbol
  dup_list=gene_output[duplicated(gene_output$ENSEMBL),]
  filter_list=dup_list[is.na(dup_list$SYMBOL),]$ENTREZID
  gene_output2=subset(gene_output, ENTREZID %ni% filter_list)
  
  # fill in missing gene symbols
  filter_list=gene_output2[is.na(gene_output2$SYMBOL),]
  for (i in rownames(filter_list)){
    eid=filter_list[i,"ENSEMBL"]
    #pull symbol from deg
    gene_output2[i,"SYMBOL"]=unique(subset(sep_df,ENSEMBL==eid)$SYMBOL)
  }
  
  # fill in missing ENSEMBL
  filter_list=gene_output2[is.na(gene_output2$ENSEMBL),]
  for (i in rownames(filter_list)){
    sym=filter_list[i,"SYMBOL"]
    #pull symbol from deg
    gene_output2[i,"ENSEMBL"]=subset(sep_df,SYMBOL==sym)$ENSEMBL
  }
  
  #remove dups
  output_df=gene_output2[!duplicated(gene_output2[,c("ENSEMBL")]),]
  
  # rename cols to match deg
  colnames(output_df)=c("ENTREZID","ENSEMBL","gene")
  return(output_df)
}

#set genelist for GSEA
deg2geneList<-function(deg,t2g){
  # create refs of entrez:genes
  gene_ref_db=capture_entrezids(deg)
  
  #add entrezs to deg df
  deg_anno_df=merge(gene_ref_db,deg,by="gene")
  
  # create genelist
  gsea_genelist=deg_anno_df$log2fc
  
  if ((t2g=="C1") | (t2g=="C2:BIOCARTA") | (t2g=="H")){
    names(gsea_genelist)=as.character(deg_anno_df$gene)
  } else{
    names(gsea_genelist)=as.character(deg_anno_df$ENTREZID)
  }
  gsea_genelist=sort(gsea_genelist,decreasing=TRUE)
  
  return(gsea_genelist)
}

# set the annotation dbs
db_lookup<-function(t2g){
  print(t2g)
  # generate gene lists for C1
  # generate gene lists for C2 with subtypes biocarta, kegg, reactome, wiki
  # generate gene lists for C5 with subtypes MF, BP, CC
  # generate gene lists for Hallmark
  if (t2g=="C1"){
    db_out=msigdbr(species = species_in, category = "C1") %>% 
      dplyr::select(gs_name,ensembl_gene)
  } else if (t2g=="C2:BIOCARTA"){
    db_out=msigdbr(species = species, category = "C2", subcategory = "BIOCARTA") %>% 
      dplyr::select(gs_name,ensembl_gene)
  } else if (t2g=="C2:KEGG"){
    db_out=msigdbr(species = species, category = "C2", subcategory = "KEGG") %>% 
      dplyr::select(gs_name,ensembl_gene)
  } else if (t2g=="C2:REACTOME"){
    db_out=msigdbr(species = species, category = "C2", subcategory = "REACTOME") %>%
      dplyr::select(gs_name,ensembl_gene)
  } else if (t2g=="C2:WIKIPATHWAYS"){
    db_out=msigdbr(species = species, category = "C2", subcategory = "WIKIPATHWAYS") %>%
      dplyr::select(gs_name,ensembl_gene)
  } else if (t2g=="C5:MF"){
    db_out=msigdbr(species = species,  category = "C5", subcategory = "GO:MF") %>%
      dplyr::select(gs_name,ensembl_gene)
  } else if (t2g=="C5:BP"){
    db_out=msigdbr(species = species,  category = "C5", subcategory = "GO:BP") %>%
      dplyr::select(gs_name,ensembl_gene)
  } else if (t2g=="C5:CC"){
    db_out=msigdbr(species = species,  category = "C5", subcategory = "GO:CC") %>%
      dplyr::select(gs_name,ensembl_gene)
  } else if (t2g=="H"){
    db_out=msigdbr(species = species, category = "H") %>% 
      dplyr::select(gs_name,ensembl_gene)  
  } else{
    print("DB does not exist. Please review")
  }
  
  return(db_out)
}

# create a small legend
addSmallLegend <- function(myPlot, pointSize = 2, textSize = 5, spaceLegend = 0.1) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_text(size = textSize), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
}

# plot ORA
ora_plus_plot_rna <- function(gl,t2g,contrast_in,n_show=3){
  # pull the DB
  pulled_db=db_lookup(t2g)
  
  # run ORA
  result=enricher(gene=gl, TERM2GENE=pulled_db, pvalueCutoff = padj_cutoff)
  resultdf=as.data.frame(result)
  
  # write out pathways file
  ttl_abbrev=sub(" ","_",sub(":","_",t2g))
  fpath=paste0(output_rna_dir,"ORA_",contrast_in[1],"-",contrast_in[2],"_table_",ttl_abbrev,".txt")
  write.table(resultdf,file=fpath,quote=FALSE,sep="\t",row.names = FALSE,col.names = TRUE)
  
  # create dotplot if pathways are sig
  if(nrow(resultdf)==0){
    p1 = ggparagraph( paste0("\n\n\n No Sig Results for ", "ORA:",t2g,"\n-",
                             contrast_in[1],"-",contrast_in[2]), 
                      color = NULL, size = 20, face = "bold", 
                      family = NULL, lineheight = NULL)
  } else{
    p1 = dotplot(result,
                 title=paste0(contrast_in[1],"-",contrast_in[2],"\nORA:",t2g),
                 font.size = 6, showCategory=n_show)
  }
  # add small legend
  pf = addSmallLegend(p1)
  return(pf)
}

ora_plus_plot_combo <- function(sigGeneList,db_id,contrast_in,n_show=5){
  # pull the DB
  pulled_db=db_lookup(db_id)
  
  # check if gene list is in database
  db_check=pulled_db$ensembl_gene[pulled_db$ensembl_gene %in% sigGeneList]
  if (length(db_check)>1){
    # run ORA
    result=clusterProfiler::enricher(gene=sigGeneList, TERM2GENE=pulled_db, 
                                     pvalueCutoff = padj_cutoff, 
                                     minGSSize = minSize_gene_set)
    resultdf=as.data.frame(result)
  } else{
    resultdf=data.frame()
  }
  
  
  # write out pathways file
  fpath=paste0(output_dir,"top_pathways_ORA_",fshort,"_",contrast_in[1],"_table.txt")
  if(nrow(resultdf)>0){
    if (file.exists(fpath)){
      add_to_df=read.csv(fpath,sep="\t")
      write.table(rbind(add_to_df,resultdf[,colnames(add_to_df)]),
                  file=fpath,quote=FALSE,sep="\t",row.names = FALSE,col.names = TRUE)
    } else{
      write.table(resultdf,
                  file=fpath,quote=FALSE,sep="\t",row.names = FALSE,col.names = TRUE)
    }
  }
  
  # create dotplot if pathways are sig
  unique_pvals=length(unique(resultdf[c(1:n_show),]$p.adjust))
  if((nrow(resultdf)==0) | (n_show-1>unique_pvals)){
    p1 = ggparagraph( paste0("\n\n\n No Sig Results for ", "-ORA:",db_id,"\n-",
                             contrast_in[1],"-",contrast_in[2]), 
                      color = NULL, size = 20, face = "bold", 
                      family = NULL, lineheight = NULL)
  } else{
    p1 = dotplot(result,
                 title=paste0(contrast_in[1],"\nORA:",db_id),
                 font.size = 6, showCategory=n_show)
  }
  # add small legend
  pf = addSmallLegend(p1)
  return(pf)
}

# plot GSEA
gsea_plus_plot <- function(gl,t2g,contrast_in,select_flag="OFF"){
  # shorthand species
  if (species_in=="Homo sapiens"){
    species_short="hsa"
    org_db="org.Hs.eg.db"
  } else if (species_in=="Mus musculus"){
    species_short="mmu"
    org_db="org.Mm.eg.db"
  } else{
    print (paste0(species_in,": species not found"))
  }
  
  # run GSEA
  if (t2g=="C2:KEGG"){
    #https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-kegg.html
    result=gseKEGG(geneList=gl,
                   pvalueCutoff=padj_cutoff,
                   eps=0,
                   pAdjustMethod="BH", 
                   organism=species_short,
                   verbose=FALSE)
    
  } else if (t2g=="C2:REACTOME"){
    #https://yulab-smu.top/biomedical-knowledge-mining-book/reactomepa.html
    result=gsePathway(gene=gl, 
                      pvalueCutoff = padj_cutoff,
                      eps=0,
                      pAdjustMethod = "BH", 
                      verbose = FALSE)
    
  } else if (t2g=="C2:WIKIPATHWAYS"){
    #https://yulab-smu.top/biomedical-knowledge-mining-book/wikipathways-analysis.html
    result=gseWP(gene=gl, 
                 pvalueCutoff = padj_cutoff,
                 eps=0,
                 pAdjustMethod = "BH",
                 organism=species_in)
    
  } else if ((t2g=="C5:MF") | (t2g=="C5:BP") | (t2g=="C5:CC")){
    #https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-go.html
    ont_id=strsplit(t2g,":")[[1]][2]
    result=gseGO(geneList=gl,
                 pvalueCutof = padj_cutoff,
                 eps=0,
                 pAdjustMethod = "BH",
                 OrgDb= get(org_db),
                 ont=ont_id,
                 verbose= FALSE)
    
  } else if ((t2g=="C1") | (t2g=="C2:BIOCARTA") | (t2g=="H")){
    pulled_db=db_lookup(t2g)
    result=GSEA(geneList=gl,
                pvalueCutoff = padj_cutoff,
                eps=0,
                pAdjustMethod = "BH",
                TERM2GENE = pulled_db)
  } else {
    print (paste0(t2g, ": DB selected is not valid"))
  }
  
  # breakpoint for the select_pathway function versus standard analysis
  if (select_flag=="ON"){ 
    return(result)
  } else{
    # save datatable
    resultdf=as.data.frame(result)
    ttl_abbrev=sub(" ","_",sub(":","_",t2g))
    fpath=paste0(output_rna_dir,"GSEA_",contrast_in[1],"-",contrast_in[2],"_table_",ttl_abbrev,".txt")
    write.table(resultdf,file=fpath,quote=FALSE,sep="\t",row.names = FALSE,col.names = TRUE)
    
    # create dot plots for all DB, ridgeplots for specific DB's
    if(nrow(result)==0){
      pf = ggparagraph( paste0("\n\n\n No Sig Results for GSEA:",t2g,
                               "\n-",contrast_in[1]), 
                        size = 20, face = "bold")
    } else{
      p1 = dotplot(result,
                   title=paste0(contrast_in,"\nGSEA:",t2g),
                   font.size = 6, showCategory=2, split=".sign",orderBy="p.adjust") +
        facet_grid(.~.sign)
      p2 = ridgeplot(result, label_format = 30, showCategory = 4, orderBy="p.adjust") +
        labs(x = "Enrichment distribution for top 5 pathways") + 
        theme(text = element_text(size=6),
              axis.text.x = element_text(size=6),
              axis.text.y = element_text(size=5.5))
      pcol <- cowplot::plot_grid(
        p1 + theme(legend.position="none"),
        p2 + theme(legend.position="none"),
        nrow = 2
      )
      legend<-get_legend(p1)
      pf=cowplot::plot_grid(pcol,legend,rel_widths = c(3, .4))
    }
    return(pf)
  }
}

# save plots
print_save_plots<-function(plot_list,contrast_id,type_in){
  if (exists(output_rna_dir)){
    output_dir=output_rna_dir
  }
  # set figure caption letter
  p=1
  
  # ORA plots need to be merged two plots per figure
  if (type_in == "ORA"){
    for (i in seq(from=1,to=length(plot_list),by=2)){
      #if it's the final image, no merging needed
      if ((i+1)>length(plot_list)){
        pf=cowplot::plot_grid(plot_list[[i]],
                              ncol=1, labels=LETTERS[p])
      } else{
        p1=plot_list[[i]]
        p2=plot_list[[i+1]]
        pf=cowplot::plot_grid(p1,p2,ncol=1,
                              labels=LETTERS[p])
      }
      print(pf)
      ggsave(filename = paste0(output_dir, type_in,"_",
                               contrast_id[1],"-",contrast_id[2], 
                               "_dotplot_",LETTERS[p],".png"),
             height = 8.90, width = 12.80, device = "png", plot = pf)
      # increase counter
      p=p+1
    }
  } else{
    # GSEA plots are already merged, one plot into two, just print
    for (i in seq(from=1,to=length(plot_list))){
      pf = cowplot::plot_grid(plot_list[[i]],
                              ncol=1, labels=LETTERS[p])
      
      #print and save
      print(pf)
      ggsave(filename = paste0(output_dir, type_in,"_",
                               contrast_id[1],"-",contrast_id[2], 
                               "_dotplot_",LETTERS[p],".png"),
             height = 8.90, width = 12.80, device = "png", plot = pf)
      # increase counter
      p=p+1
    }
  }
}

# create output dt that summarized pathways
create_dts<-function(type_in,t2g,contras,n_in,db_list){
  
  # read in datatable created during GSEA/ORA plotting
  merged_df=data.frame()
  for (t2g in db_list){
    ttl_abbrev=sub(" ","_",sub(":","_",t2g))
    fpath=paste0(output_rna_dir,type_in,"_",contras[1],"-",contras[2],"_table_",ttl_abbrev,".txt")
    tmp_df=read.csv(fpath,sep="\t")
    
    # if there are no signifcant pathways, skip
    if (nrow(tmp_df)==0){ next }
    tmp_df$anno=t2g
    
    # merge
    merged_df=rbind(merged_df,tmp_df)
  }
  
  # split leading_edge col into three
  if (nrow(merged_df)>0){
    if (type_in=="GSEA"){
      merged_df=separate(
        merged_df,
        leading_edge,
        c("percent_included","percent_list","percent_signal"),
        sep=", "
      )
      # remove value= in cols
      merged_df$percent_included=sub("tags=","",merged_df$percent_included)
      merged_df$percent_list=sub("list=","",merged_df$percent_list)
      merged_df$percent_signal=sub("signal=","",merged_df$percent_signal)
      
      # round cols
      col_list=c("p.adjust","enrichmentScore","NES")
      for (colid in col_list){
        merged_df[,colid]=signif(merged_df[,colid], digits=3)
      }
      
      # select cols
      output_df=merged_df[,c("anno","ID","Description","setSize","percent_included","p.adjust",
                             "enrichmentScore","NES","core_enrichment")]
      colnames(output_df)=c("anno_db","ID","Desc","total_genes","percent_included","p.adj",
                            "enrichmentScore","NES","genes_included")
    } else{
      # pull values for genes included, genes in set
      merged_df=separate(
        merged_df,
        BgRatio,
        c("set_total","n_bgratio"),
        sep="/"
      )
      # calculate percent genes included
      merged_df$set_total=as.numeric(merged_df$set_total)
      merged_df$percent_included=(merged_df$Count/merged_df$set_total)*100
      
      #round cols
      col_list=c("pvalue","p.adjust","qvalue","percent_included")
      for (colid in col_list){
        merged_df[,colid]=signif(merged_df[,colid], digits=3)
      }
      merged_df$percent_included=paste0(merged_df$percent_included,"%")
      
      # select cols
      output_df=merged_df[,c("anno","ID","set_total","percent_included","p.adjust","geneID")]
      colnames(output_df)=c("anno_db","ID","total_genes","percent_included","p.adj","genes_included")
    }
    # sort by p.adjust
    output_df=output_df[order(output_df$p.adj),]
    
    # if number of pathways is less than n_in, adjust
    if (nrow(output_df)<n_in){
      n_in=nrow(output_df)
    }
    
    # create DT
    caption_title=paste0("Top ", n_in, " Pathways for all annotation databases (", type_in, ")\n")
    p <- DT::datatable(output_df[1:n_in,], extensions = 'Responsive', 
                       caption=htmltools::tags$caption(paste0(caption_title, contras[1],"_vs_", contras[2]) ,
                                                       style="color:gray; font-size: 18px" ),rownames=F)
  } else{
    p="No Significant Genes"
  }
  return(p)
}

# main function
main_gsea_ora_function_rna<-function(cntrl_in,treat_in,db_list,top_path_value,ORA_flag="",GSEA_flag=""){
  
  # set contrast
  contras=c(treat_in,cntrl_in)
  
  # read in deg
  deg_file = paste0(input_dir, "DEG_", contras[1],"-",contras[2],"_0.5_0.5/",analysis_type,
                    "_DEG_",contras[1],"-",contras[2],"_all_genes.txt")
  deg=read.csv(deg_file,header=TRUE,sep="\t")
  
  # run ORA
  o=list()
  if (ORA_flag=="ON"){
    # create ORA genelist
    siggenes=subset(deg,fdr <= padj_cutoff)
    siggenes=subset(deg,(fc < -log2fc_cutoff) | (fc > log2fc_cutoff))
    sigGeneList=siggenes$gene
    
    # for each annotation db, run ORA, save plots
    i=1
    for (db_id in db_list){
      o[[i]]=ora_plus_plot_rna(gl=sigGeneList,t2g=db_id,contrast_in=contras)
      i=i+1
    }
    
    # print,save plots
    print_save_plots(o,contras,"ORA")
    
    # print DT
    create_dts("ORA",t2g,contras,top_path_value,db_list)
  }
  
  # run GSEA
  g=list()
  if (GSEA_flag=="ON"){
    
    # for each annotation db, create GSEA gene list - must be done within loop
    # to handle differences in naming list 
    i=1
    for (db_id in db_list){
      # create GSEA genelist
      gsea_genelist=deg2geneList(deg,t2g=db_id)
      
      # run GSEA, save plots
      g[[i]]=gsea_plus_plot(gl=gsea_genelist,t2g=db_id,contrast_in=contras)
      i=i+1
    }
    
    # print,save plots
    print_save_plots(g,contras,"GSEA")
    
    # print DT
    create_dts("GSEA",t2g,contras,top_path_value,db_list)
  }
}

# create ranked df
prep_genelist_df<-function(subset_list,sig_type,analysis_type,use_gene_list="NO"){
  # input dfs
  fpath=paste0(output_dir,"overlap_",fshort,"_",contrast_id_car,"_",contrast_id_rna,".csv")
  overlap_df=read.csv(fpath)
  
  fpath=paste0(output_rna_dir,"DESeq2_",contrast_id_rna,"_DEG_allgenes_res1.txt")
  rna_df=read.csv(fpath,sep="\t")
  rna_df=rna_df %>% 
    separate("X",c("ENSEMBL","SYMBOL"),sep="[|]") %>% 
    separate("ENSEMBL",c("ENSEMBL","ID"),sep="[.]")
  
  # subset CAR for annotation of interest
  if ("all" %in% subset_list){
    sub_df=overlap_df
  } else{
    sub_df=subset(overlap_df,annotation %in% subset_list)
  }
  
  # sort and remove dups
  dedup_df=subset(sub_df,overlap_type=="overlap")
  if (sig_type=="both"){
    sig_RNA=subset(overlap_df,abs(log2FC_rna)>log2fc_cutoff_rna & padj_rna<padj_cutoff)$ENSEMBL
    sig_CAR=subset(overlap_df,abs(log2FC_car)>log2fc_cutoff_car & padj_car<padj_cutoff)$ENSEMBL
    dedup_df=subset(dedup_df,ENSEMBL %in% intersect(sig_RNA,sig_CAR))
  }
  dedup_df=dedup_df[order(dedup_df$padj_rna),]
  dedup_df=dedup_df[!duplicated(dedup_df$ENSEMBL),]
  
  # pull significant genes ENSEMBL
  sorted_overlap_EID=dedup_df$ENSEMBL
  
  # for GSEA analysis, output must include a ranked list of overlapping significant genes
  # AND all RNA genes in analysis. for ORA analysis output must only include the
  # overlapping significant genes
  if(analysis_type=="GSEA"){
    # remove all overlap EIDs and create list of RNA only
    sub_rna_df=rna_df[order(rna_df$padj),]
    sub_rna_df=sub_rna_df[!duplicated(sub_rna_df$ENSEMBL),]
    sub_rna_df=subset(sub_rna_df,ENSEMBL %ni% sorted_overlap_EID)
    sub_rna_df=sub_rna_df%>%rename("log2FC_rna"=log2FoldChange)
    sorted_rna_EID=sub_rna_df[order(sub_rna_df$log2FC_rna),]$ENSEMBL
    
    # create final df
    output_df=full_join(sub_rna_df,dedup_df) %>%
      subset(ENSEMBL %in% c(sorted_overlap_EID,sorted_rna_EID))
    rownames(output_df)=output_df$ENSEMBL
    output_df=output_df[c(sorted_overlap_EID,sorted_rna_EID),]
  } else if (analysis_type=="ORA"){
    output_df=dedup_df
  }
  
  if(use_gene_list=="YES"){
    print("Using subset gene list")
    pander(gene_list,style="rmarkdown")
    output_df=subset(output_df,SYMBOL %in% gene_list)
  }
  
  return(output_df)
}

# create top pathways
create_plot_fgsea<-function(t2g,fgseaRes,ranked_list,msigdbr_list){
  # run pathway analysis
  topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=5), pathway]
  topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=5), pathway]
  
  # generate plots, DT's for up and down
  path_name=str_wrap(gsub("_"," ",topPathwaysUp[1]), width = 30)
  p1 = plotEnrichment(msigdbr_list[[topPathwaysUp[1]]], ranked_list) + 
    labs(title=paste0("Up-regulated\n",path_name),cex=.5)
  path_name=str_wrap(gsub("_"," ",topPathwaysDown[1]), width = 30)
  p2 = plotEnrichment(msigdbr_list[[topPathwaysDown[1]]], ranked_list) + 
    labs(title=paste0("Down-regulated\n",path_name))
  title1=text_grob(paste0("Top pathways for ", t2g), size = 15, face = "bold")
  grid.arrange(
    p1,p2,
    top=title1,
    nrow=1
  )
}

# create datatables
create_dts_fgsea<-function(fgseaRes,db_id){
  # subset pathways
  topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=5), pathway]
  topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=5), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  
  # generate table for all
  sub_df=as.data.frame(fgseaRes)
  sub_df=subset(sub_df,pathway %in% topPathways)
  sub_df=sub_df[order(sub_df$ES),]
  col_select=c("padj","ES","NES","log2err")
  
  for (colid in col_select){
    sub_df[,colid]=signif(sub_df[,colid], digits=3)
  }
  sub_df$db=db_id
  
  col_select=c("db","pathway","size","padj","ES","NES","log2err")
  return(sub_df[,col_select])
}

# main function
main_gsea_ora_function_combo<-function(subset_list,sig_type,db_list,analysis_type,use_gene_list="NO"){
  
  # prep ranked gene list for GSEA or subset sig list for ORA
  subset_df=prep_genelist_df(subset_list,sig_type,analysis_type=analysis_type,use_gene_list)
  
  # create annodb top pathway df
  merged_df=data.frame()
  
  # pull db and generate list
  for (db_id in db_list){
    print(paste0("--working on ", db_id))
    anno_db=db_lookup(db_id)
    msigdbr_list=split(anno_db$ensembl_gene,anno_db$gs_name)
    
    if (analysis_type=="GSEA"){
      # create ranked input
      ranked_list=subset_df$log2FC_rna
      names(ranked_list)=subset_df$ENSEMBL  
      
      # run fgsea
      #https://cran.r-project.org/web/packages/msigdbr/vignettes/msigdbr-intro.html
      #https://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html
      fgseaRes <- fgsea(msigdbr_list, ranked_list,minSize  = minSize_gene_set)
      
      # plot 
      create_plot_fgsea(db_id,fgseaRes,ranked_list,msigdbr_list)
      
      #create merged top pathway df
      merged_df=rbind(merged_df,create_dts_fgsea(fgseaRes,db_id))
    } else{
      sigGeneList=subset_df$ENSEMBL
      
      # for each annotation db, run ORA, save plots
      o=list()
      i=1
      for (db_id in db_list){
        o[[i]]=ora_plus_plot_combo(sigGeneList=sigGeneList,db_id=db_id,
                                   contrast_in=contrast_id_rna)
        i=i+1
      }
    }
  }
  
  # print,save ORA plots
  if (analysis_type=="ORA"){
    print_save_plots(o,contrast_id_rna,"ORA")
    fpath=paste0(output_dir,"top_pathways_ORA_",fshort,"_",contrast_id_rna,"_table.txt")
    print(fpath)
    caption_title=paste0("Top pathways across all annotation databases for ",
                         fshort, " genes using ", analysis_type," analysis.")
    DT::datatable(read.csv(fpath,sep="\t"), extensions = 'Responsive', 
                  caption=htmltools::tags$caption(paste0(caption_title),
                                                  style="color:gray; font-size: 18px" ),
                  rownames=F)
  } else{
    # print top pathway df
    fpath=paste0(output_dir,"top_pathways_",analysis_type,"_",fshort,"_",contrast_id_car,"_",contrast_id_rna,".csv")
    print(fpath)
    write.table(merged_df,fpath,sep=",")
    caption_title=paste0("Top pathways across all annotation databases for ",
                         fshort, " genes using ", analysis_type," analysis.")
    DT::datatable(merged_df, extensions = 'Responsive', 
                  caption=htmltools::tags$caption(paste0(caption_title),
                                                  style="color:gray; font-size: 18px" ),
                  rownames=F)
  }
}

######################################################################
# functions secondary pathway analysis
######################################################################
# creates heatmap formatted df with only log2fc values
create_heatmap_df<-function(df_in,extra_filter=""){
  df_out=dplyr::select(df_in,contains("log"))
  
  # cleanup col names
  output_list=sub('--log2FoldChange', '',colnames(df_out))
  cntrl=strsplit(output_list,"_vs_")
  for (ct in cntrl){
    output_list=sub(ct[2], '',output_list)
    output_list=sub('_vs_', '',output_list)
  }
  
  # use extra filter
  if (extra_filter != ""){
    output_list=sub(extra_filter, '',output_list)
  }
  
  colnames(df_out)=output_list
  return(df_out)
}

# creates df of log2fc and pvalues
create_output_df<-function(df_in,n_in,extra_filter=""){
  output_list=colnames(df_in)
  
  # split treat from control
  cntrl=strsplit(output_list,"_vs_")
  for (ct in cntrl){
    metric=strsplit(ct[2],"--")[[1]][2]
    output_list=sub(ct[2], paste0("_",metric),output_list)
    output_list=sub('_vs_', '',output_list)
  }
  
  # cleanup col names
  output_list=sub('log2FoldChange', 'log2FC',output_list)
  
  # use extra filter
  extra_filter="_without_IFNb"
  if (extra_filter != ""){
    output_list=sub(extra_filter, '',output_list)
  }
  
  colnames(df_in)=output_list
  
  # round df
  for (colid in colnames(df_in)){
    df_in[,colid]=signif(df_in[,colid], digits=3)
  }
  
  # if number of pathways is less than n_in, adjust
  if (nrow(df_in)<n_in){
    n_in=nrow(df_in)
  }
  
  # print out table
  caption_title=paste0("Expression values for ",n_in," genes")
  p = DT::datatable(df_in, extensions = 'Responsive', 
                    caption=htmltools::tags$caption(paste0(caption_title) ,
                                                    style="color:gray; font-size: 18px" ))
  return(p)
}

# creates heatmap
generate_heat_map<-function(df_in,show_names="ON",title_in="",cluster_by_rows="ON",fpath=""){
  
  ####################
  # formatting
  #####################
  # Overwrite pheatmaps default draw_colnames with new version
  assignInNamespace(x="draw_colnames", value="draw_colnames_45",ns=asNamespace("pheatmap")) 
  
  # Heatmap Color Gradients 
  paletteLength <- 1000
  mycolors <- colorRampPalette(c("blue","white","red"), interpolate = "linear")(paletteLength)
  
  ####################
  # metadata
  ####################
  # Creating Dataframe to map samplenames to groups
  meta = groups_df
  groups <- data.frame(as.factor(meta$group))
  colnames(groups) <- "Groups"
  rownames(groups) <- meta$sampleid
  
  # Creating Group Column Annotation Colors
  columnColors <- c("lightpink","lightblue","orange","purple","red","green")
  names(columnColors) <- unique(groups$Groups)
  anno_colors <- list(Groups = columnColors)
  
  # set title
  if(title_in==""){
    title_in=paste0("Significant Genes (N=",nrow(df_in),")")
  }
  ####################
  # function
  ####################
  if (show_names=="OFF" && cluster_by_rows=="ON"){
    p=pheatmap(df_in, 
               scale = "none", main=title_in,
               cellwidth = 30, fontsize = 12, fontsize_row = 7, fontsize_col = 8, color = mycolors, 
               border_color = "NA",cluster_cols=F,annotation_colors = anno_colors, show_rownames = FALSE)
  } else if (show_names=="ON" && cluster_by_rows=="ON") {
    p=pheatmap(df_in, 
               scale = "none", main=title_in,
               cellwidth = 30, fontsize = 12, fontsize_row = 5, fontsize_col = 8, color = mycolors, 
               border_color = "NA",cluster_cols=F,annotation_colors = anno_colors, show_rownames = TRUE)
  } else if (show_names=="OFF" && cluster_by_rows=="OFF") {
    p=pheatmap(df_in, 
               scale = "none", main=title_in,
               cellwidth = 30, fontsize = 12, fontsize_row = 7, fontsize_col = 8, color = mycolors, 
               border_color = "NA",cluster_cols=F,cluster_rows=F,annotation_colors = anno_colors, 
               show_rownames = FALSE)
  }
  if (fpath==""){
    fpath=paste0(img_dir,"heatmap_",contrast_id,".pdf")
  }
  save_pheatmap_pdf(p, fpath)
}

# creates heatmap for multiple replicates
generate_replicate_heatmaps<-function(contrast_id,peak_type,scale,gene_list_subset="",sample_subset=""){
  # split contrast
  contrast_1=strsplit(contrast_id,"_vs_")[[1]][1]
  contrast_2=strsplit(contrast_id,"_vs_")[[1]][2]
  
  # define extension
  contrast_extension=paste0(contrast_id,"__",dedup_status,"__",peak_type,"_peaks.bed")
  
  # read in counts matrix
  fpath=paste0(output_car_dir,"DESEQ_norm_counts_",contrast_id,".csv")
  counts_matrix=read.csv(fpath)
  head(counts_matrix)
  colnames(counts_matrix)=c("peakID",colnames(counts_matrix)[2:ncol(counts_matrix)])
  
  # replace X at the beginning of cols
  colnames(counts_matrix)=gsub("X","",colnames(counts_matrix))
  head(counts_matrix)
  
  # subset for samples, if needed
  if (length(sample_subset)>1){
    sample_list=sample_subset
    counts_matrix=counts_matrix[,c("peakID",sample_list)]
  } else{
    sample_list=subset(groups_df,group%in% c(contrast_1,contrast_2))$sampleid
  }
  
  # add annotation information
  fpath=paste0(output_car_dir,"peak_annotation_",contrast_id,".csv")
  peak_df=read.csv(fpath)[,c("peakID","SYMBOL","shortAnno")]
  head(peak_df)
  counts_matrix_anno=merge(counts_matrix,peak_df,by="peakID")
  head(counts_matrix_anno)
  
  # pull rowname
  rownames(counts_matrix_anno)=counts_matrix_anno$peakID
  counts_matrix=counts_matrix_anno[,2:ncol(counts_matrix_anno)]
  head(counts_matrix_anno)
  
  # read in sig gene list
  fpath=paste0(output_dir,"deeptools_sig_",contrast_id,".csv")
  deeptools_gene_list=read.csv(fpath)$x
  head(deeptools_gene_list)
  
  # subset for gene list
  counts_matrix_subset=subset(counts_matrix_anno,SYMBOL %in% deeptools_gene_list)
  
  # subset for non-distal annotated peaks
  counts_matrix_subset=subset(counts_matrix_subset,shortAnno!="Distal")
  
  # if needed, subset for immune gene list
  if (length(gene_list_subset)>1){
    print("Subsetting based on gene_list provided")
    counts_matrix_subset=subset(counts_matrix_subset,SYMBOL %in% gene_list_subset)
    
    print(paste0("Total number of peaks: ", nrow(counts_matrix_subset)))
    print(paste0("Total number unique genes: ", length(unique(counts_matrix_subset$SYMBOL))))
    
    # keep only one instance per symbol
    counts_matrix_subset=counts_matrix_subset[!duplicated(counts_matrix_subset$SYMBOL), ]
    print(paste0("Total number of peaks after subsetting only unique genes: ", nrow(counts_matrix_subset)))
    
    # rename rows for printing
    rownames(counts_matrix_subset)=make.unique(paste0(counts_matrix_subset$SYMBOL," (",counts_matrix_subset$shortAnno,")"))
  } else{
    print(paste0("Total number of peaks: ", nrow(counts_matrix_subset)))
    print(paste0("Total number unique genes: ", length(unique(counts_matrix_subset$SYMBOL))))
  }
  
  # scale if necessary
  if (scale=="ON"){
    print("performing scaling")
    # ztransform df
    counts_matrix_complete=t(scale(t(counts_matrix_subset[,sample_list])))
    (head(counts_matrix_complete))
  } else{
    print("no scaling will be performed")
    (head(counts_matrix_subset[,sample_list]))
    counts_matrix_complete=counts_matrix_subset[,sample_list]
  }
  
  # define output name, run heatmap
  if (scale=="ON"){
    fpath=paste0(img_dir,"heatmap_withscale_")
  } else{
    fpath=paste0(img_dir,"heatmap_withoutscale_")
  }
  if (length(gene_list_subset)>1){
    fpath=paste0(fpath,"subset_",gene_list_name,"_",contrast_id,".pdf")
    print(fpath)
    
    # generate heatmap
    generate_heat_map(counts_matrix_complete,
                      show_names="ON",
                      title_in="",
                      cluster_by_rows="ON",
                      fpath=fpath)
  } else{
    fpath=paste0(fpath,contrast_id,".pdf")
    
    # generate heatmap
    generate_heat_map(counts_matrix_complete,
                      show_names="OFF",
                      title_in="",
                      cluster_by_rows="ON",
                      fpath=fpath)
  }
}


#create heatmap and DT
generate_heat_map_select<-function(select_deg,contras){
  # prep df for heatmap generation
  rownames(select_deg)=select_deg$ensid_gene
  log_name=paste0(contras[1],"_vs_",contras[2],"--log2FoldChange")
  p_name=paste0(contras[1],"_vs_",contras[2],"--padj")
  select_deg=select_deg[,c("log2fc","fdr")]
  colnames(select_deg)=c(log_name,p_name)
  
  # create heatmap df of sig genes in pathway list
  sig_df=subset(select_deg, get(log_name) > log2fc_cutoff | get(log_name) < -log2fc_cutoff)
  sig_df=subset(sig_df, get(p_name) < padj_cutoff )
  heat_df=create_heatmap_df(sig_df,"")
  
  # run functions
  generate_heat_map(heat_df)
  p = create_output_df(select_deg,nrow(select_deg),"")
  return(p)
}

# create gseaplot
gsea_plus_plots_select<-function(deg,t2g,path_id,contras){
  
  # create GSEA genelist
  gsea_genelist=deg2geneList(deg,t2g=t2g)
  
  # run GSEA, save plots
  result=gsea_plus_plot(gl=gsea_genelist,t2g=t2g,contrast_in=contras,select_flag="ON")
  
  #get rowname of pathway
  result_df=as.data.frame(result)
  rownames(result_df) <- NULL
  row_id=as.numeric(rownames(subset(result_df,ID==path_id))[[1]])
  path_desc=subset(result_df,ID==path_id)$Description[[1]]
  
  # create plot
  p1=gseaplot(result, by = "all", title = paste0(path_id,": ",path_desc),
              geneSetID =row_id)
  pf=p1&theme(text = element_text(size=8),
              axis.text.x = element_text(size=8),
              axis.text.y = element_text(size=8),
              axis.text =element_text(size=8),
              axis.title=element_text(size=8,face="bold")) 
  print(pf)
}

main_selectpath_function<-function(cntrl_in,treat_in,type_in,t2g,path_id){
  t2g="C1"
  type_in="GSEA"
  path_id="chr18q21"
  cntrl_in=cntrl
  treat_in=treatment
  
  # set contrast
  contras=c(treat_in,cntrl_in)
  
  # read in datatable created during GSEA/ORA plotting, get list of gene names in pathway selected
  ttl_abbrev=sub(" ","_",sub(":","_",t2g))
  fpath=paste0(output_rna_dir,type_in,"_",contras[1],"-",contras[2],"_table_",ttl_abbrev,".txt")
  path_df=read.csv(fpath,sep="\t")
  
  #check pathway exists
  if (nrow(subset(path_df,ID==path_id))==0){
    stop(paste0("The selected pathway (",path_id,") does not exist in the annotation database (",t2g,
                "). Please select a valid combination"))
  }
  
  # create gene list from pathway
  genes_in_pathway=strsplit(subset(path_df,ID==path_id)$core_enrichment,"/")[[1]]
  
  # read in created deg
  fpath=paste0(output_rna_dir,"DESeq2_",contras[1],"-",contras[2],"_DEG_allgenes.txt")
  deg=read.csv(fpath,header=TRUE,sep="\t")
  
  # subset deg for genes,
  #convert ENTREZID if necessary
  special_dbs=c("C1","C2:BIOCARTA","H")
  if (t2g %ni% special_dbs){
    gene_ref_db = capture_entrezids(deg)
    genes_in_pathway=subset(gene_ref_db,ENTREZID %in% genes_in_pathway)$gene
  } 
  select_deg=subset(deg, gene %in% genes_in_pathway)
  
  # create heatmaps
  p = generate_heat_map_select(select_deg,contras)
  
  # create gseaPlot
  gsea_plus_plots_select(deg,t2g,path_id,contras)
  
  # return DT
  return(p)
}
