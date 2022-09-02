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
# pie charts for gene lists
############################################################
# create pie chart for significant peaks within PI gene list for promoters
plot_pies_genelist<-function(df_in,y_in,percent_in,fill_in,title_in){
  p = ggplot(df_in, aes(x = "" , y = get(y_in), fill = fct_inorder(get(fill_in)))) +
    geom_col(width = 1, color = 1) +
    coord_polar(theta = "y") +
    scale_fill_brewer(palette = "Pastel1") +
    geom_label_repel(data = df_in,
                     aes(y = pos, label = paste0(get(percent_in), "%")),
                     size = 4, nudge_x = 1, show.legend = FALSE) +
    guides(fill = guide_legend(title = "Group")) +
    ggtitle(title_in) +
    theme_void()
  
  return(p)
}

# main function to generate pie charts
generate_piecharts_from_genelist<-function(sample_id,subset_type="Promoter"){
  #sample_id=contrast_id
  # bring in dfs
  fpath=paste0(output_car_dir,"collapsed_df.csv")
  collapsed_df=read.csv(fpath,sep=",")
  
  fpath=paste0(output_car_dir,"collapsed_pi_df.csv")
  collapsed_pi_df=read.csv(fpath,sep=",")
  
  # subset for sample
  sub_col_df=subset(collapsed_df,sample==sample_id)
  sub_pi_df=subset(collapsed_pi_df,sample==sample_id)
  
  # calculate counts
  REPRESSORS=nrow(subset(sub_pi_df,gene_list=="REPRESSORS" & shortAnno==subset_type))
  REPRESSORS_total=nrow(subset(gene_list,Set=="REPRESSORS"))
  REPRESSORS_perc=round((REPRESSORS/REPRESSORS_total)*100,2)
  ACCELERATORS=nrow(subset(sub_pi_df,gene_list=="ACCELERATORS"& shortAnno==subset_type))
  ACCELERATORS_total=nrow(subset(gene_list,Set=="ACCELERATORS"))
  ACCELERATORS_perc=round((ACCELERATORS/ACCELERATORS_total)*100,2)
  INVASION=nrow(subset(sub_pi_df,gene_list=="INVASION"& shortAnno==subset_type))
  INVASION_total=nrow(subset(gene_list,Set=="INVASION"))
  INVASION_perc=round((INVASION/INVASION_total)*100,2)
  total=unique(sub_col_df$total)
  other=total-sum(REPRESSORS,ACCELERATORS,INVASION)
  
  # create df
  tmp_df=data.frame("IDENTIFIED","REPRESSORS",REPRESSORS,round((REPRESSORS/total)*100,2),REPRESSORS_perc)
  tmp_df=rbind(tmp_df,c("IDENTIFIED","ACCELERATORS",ACCELERATORS,round((ACCELERATORS/total)*100,2),ACCELERATORS_perc))
  tmp_df=rbind(tmp_df,c("IDENTIFIED","INVASION",INVASION,round((INVASION/total)*100,2),INVASION_perc))
  tmp_df=rbind(tmp_df,c("IDENTIFIED","OTHER",other,round((other/total)*100,2),0))
  tmp_df=rbind(tmp_df,c("NOT IDENTIFIED","REPRESSORS",REPRESSORS_total-REPRESSORS,0,(100-REPRESSORS_perc)))
  tmp_df=rbind(tmp_df,c("NOT IDENTIFIED","ACCELERATORS",ACCELERATORS_total-ACCELERATORS,0,100-ACCELERATORS_perc))
  tmp_df=rbind(tmp_df,c("NOT IDENTIFIED","INVASION",INVASION_total-INVASION,0,100-INVASION_perc))
  colnames(tmp_df)=c("Search","Category","n","perc","perc_list")
  tmp_df$n=as.numeric(tmp_df$n)
  
  ## all totals
  tmp_sub_df <- subset(tmp_df,Search=="IDENTIFIED") %>%
    mutate(csum = rev(cumsum(rev(n))),
           pos = n/2 + lead(csum, 1),
           pos = if_else(is.na(pos), n/2, pos))
  plot_title=paste0("All Genes (",sum(tmp_sub_df$n),")")
  p1 = plot_pies_genelist(tmp_sub_df,"n","perc","Category",plot_title)
  
  ## ACCELERATORS
  tmp_sub_df <- subset(tmp_df,Category=="ACCELERATORS") %>%
    mutate(csum = rev(cumsum(rev(n))),
           pos = n/2 + lead(csum, 1),
           pos = if_else(is.na(pos), n/2, pos))
  plot_title=paste0("    ACCELERATORS Genes (N=16)")
  p2 = plot_pies_genelist(tmp_sub_df,"n","perc_list","Search",plot_title)
  
  ## REPRESSORS
  tmp_sub_df <- subset(tmp_df,Category=="REPRESSORS") %>%
    mutate(csum = rev(cumsum(rev(n))), 
           pos = n/2 + lead(csum, 1),
           pos = if_else(is.na(pos), n/2, pos))
  plot_title=paste0("    REPRESSORS Genes (N=10)")
  p3 = plot_pies_genelist(tmp_sub_df,"n","perc_list","Search",plot_title)
  
  ## INVASION
  tmp_sub_df <- subset(tmp_df,Category=="INVASION") %>%
    mutate(csum = rev(cumsum(rev(n))), 
           pos = n/2 + lead(csum, 1),
           pos = if_else(is.na(pos), n/2, pos))
  plot_title=paste0("    INVASION Genes (N=16)")
  p4 = plot_pies_genelist(tmp_sub_df,"n","perc_list","Search",plot_title)
  
  p_final=ggarrange(p1,p2,p3,p4,
                    labels = c("A","B","C","D"),
                    ncol = 2, nrow = 2)
  plot_title=paste0("Significant Peaks by Gene Lists in ",subset_type,":\n",unique(sub_col_df$sample))
  p_final=annotate_figure(p_final, top = text_grob(plot_title, face = "bold", size = 14))
  print(p_final)
}

############################################################
# volcano plots for gene lists
############################################################
generate_volcano_plots<-function(contrast_id){
  
  # read in res from DEG merge
  fpath=paste0(output_car_dir,"DESeq2_res_",contrast_id,".csv")
  res1=read.csv(fpath,sep=",")
  colnames(res1)=c("peakID","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")
  
  # read in PEAK ANNO df from DEG merge
  fpath=paste0(output_car_dir,"peak_annotation_",contrast_id,".csv")
  pa=read.csv(fpath)
  
  # merge annotations and peaks
  results_df = full_join(res1,pa,by=c("peakID"))
  
  # filter df
  results_filtered_df = subset(results_df,(log2FoldChange>=log2fc_cutoff) | (log2FoldChange<=-log2fc_cutoff))
  results_filtered_df = subset(results_filtered_df,padj<=padj_cutoff)
  
  # add colors and annotate
  colors=brewer.pal(7,"Set1")
  anno_types=levels(as.factor(results_df$shortAnno))
  
  # set keyvals
  ## colors
  keyvals=rep("grey",times=nrow(results_filtered_df))
  names(keyvals)=rep("NS",times=length(keyvals))
  for ( i in seq(1,length(anno_types))) {
    keyvals[ abs(results_filtered_df$log2FoldChange) > log2fc_cutoff & 
               results_filtered_df$padj < padj_cutoff & results_filtered_df$shortAnno == anno_types[i] ] = colors[i]
    names(keyvals)[keyvals == colors[i]] <- anno_types[i]
  }
  ## shapes
  keyvals.shape <- ifelse(
    results_filtered_df$SYMBOL %in% gene_list$Human, 17,3)
  keyvals.shape[is.na(keyvals.shape)] <- 3
  names(keyvals.shape)[keyvals.shape == 3] <- 'Not in gene list'
  names(keyvals.shape)[keyvals.shape == 17] <- 'In gene list'
  
  # print volcano without genelist designation
  p = EnhancedVolcano(results_filtered_df,
                      lab = results_filtered_df$SYMBOL,
                      x = 'log2FoldChange',
                      y = 'padj',
                      ylab = bquote(~-Log[10] ~ FDR),
                      pCutoff = padj_cutoff,
                      FCcutoff = log2fc_cutoff,
                      labSize = 4,
                      title = contrast_id,
                      subtitle = "",
                      subtitleLabSize = 1,
                      captionLabSize = 10,
                      colCustom = keyvals,
                      colAlpha = 1,
                      legendLabels = c("NS", expression(Log[2] ~ FC), "FDR", expression(FDR ~ and ~ log[2] ~ FC)),
                      legendLabSize = 10,legendPosition = 'right'
  )
  print(p)
  
  # set labels
  log_pval=-log10(results_df$pvalue)
  y_title="-Log10 FDR"
  stat_type="pvalue"
  log_FC=results_df$log2FoldChange
  Significant=rep("1_NotSignificant",length(log_FC))
  Significant[which(results_df$pvalue<padj_cutoff & abs(results_df$log2FoldChange)>=log2fc_cutoff)]=paste0("3_LogFC_and_",stat_type)
  Significant[which(results_df$pvalue<padj_cutoff & abs(results_df$log2FoldChange)<log2fc_cutoff)]=paste0("2b_",stat_type,"_Only")
  Significant[which(results_df$pvalue>=padj_cutoff & abs(results_df$log2FoldChange)>=log2fc_cutoff)]="2a_LogFC_Only"
  gene=results_df$SYMBOL
  volcano_data=as.data.frame(cbind(gene,log_FC,log_pval,Significant))
  
  p <- plot_ly(data = volcano_data, x = log_FC, y = log_pval, text = gene,
               mode = "markers", 
               color = Significant) %>% layout(title =paste0(contrast_id),
                                               xaxis=list(title="Fold Change",
                                                          range =c(-5,5),
                                                          tickvals=c(-5,-4,-3,-2,-1,0,1,2,3,4,5),
                                                          ticktext=c('-32','-16','-8','-4','-2',
                                                                     '1','2','4','8','16','32')),
                                               yaxis=list(title=y_title,range =c(0,15)))
  return(p)
}

############################################################
# find differential overlap
############################################################
create_venn_diagrams<-function(subtitle,rna_df_in,car_df_in){
  # create gene lists
  list_of_rna_genes=rna_df_in$SYMBOL
  list_of_car_genes=car_df_in$SYMBOL
  
  # remove NA"s
  list_of_rna_genes=list_of_rna_genes[!is.na(list_of_rna_genes)]
  list_of_car_genes=list_of_car_genes[!is.na(list_of_car_genes)]
  
  # List of genes
  x <- list(A = list_of_car_genes, B = list_of_rna_genes)
  
  # Venn diagram with custom category names
  p = ggVennDiagram(x, color = 1, lwd = 0.7,
                    category.names = c("Cut&Run Genes",
                                       "RNASeq Genes")) + 
    scale_fill_gradient(low = "red", high = "blue")
  full_title=paste0("Significantly Differentiated Genes: ",subtitle," genes")
  pf = p + ggtitle(full_title)
  
  # save and print
  fpath=paste0(output_car_dir,"venndiagram_",subtitle,"_genes.png")
  ggsave(fpath,pf)
  
  print(pf)
}

create_overlapping_df<-function(rna_df_in,car_df_in,subset_type){
  
  # create gene lists
  list_of_rna_genes=rna_df_in$SYMBOL
  list_of_car_genes=car_df_in$SYMBOL
  
  # remove NA"s
  list_of_rna_genes=list_of_rna_genes[!is.na(list_of_rna_genes)]
  list_of_car_genes=list_of_car_genes[!is.na(list_of_car_genes)]
  
  # create CAR only, overlapped df
  car_col_list=c("peakID","geneChr","geneStart","geneEnd","shortAnno","ENSEMBL",
                 "SYMBOL","GENENAME","log2FoldChange","padj")
  overlap_genes_df=subset(car_df_in,SYMBOL%in%list_of_rna_genes)[,car_col_list]
  overlap_genes_df$overlap_type="overlap"
  
  only_car_df=subset(car_df_in,SYMBOL%ni%list_of_rna_genes)[,car_col_list]
  only_car_df$overlap_type="only_car"
  
  merged_car_df=full_join(overlap_genes_df,only_car_df) %>%
    rename(c("log2FC_car"=log2FoldChange,"padj_car"=padj,
             "gene_description"=GENENAME,"annotation"=shortAnno))
  
  # create RNA only, overlapped df
  rna_col_list=c("ENSEMBL","SYMBOL","log2FoldChange","padj")
  overlap_genes_df=subset(rna_df_in,SYMBOL%in%list_of_car_genes)[,rna_col_list]
  overlap_genes_df$overlap_type="overlap"
  
  only_rna_df=subset(rna_df_in,SYMBOL%ni%list_of_car_genes)[,rna_col_list]
  only_rna_df$overlap_type="only_rna"
  
  merged_rna_df=full_join(overlap_genes_df,only_rna_df) %>%
    rename(c("log2FC_rna"=log2FoldChange,"padj_rna"=padj))
  
  # fix ensemblID issue before merging to avoid duplicates
  for (rowid in rownames(merged_car_df)){
    if (is.na(merged_car_df[rowid,"ENSEMBL"])){
      lookup_SYMBOL=merged_car_df[rowid,"SYMBOL"]
      found_row=subset(merged_rna_df,SYMBOL==lookup_SYMBOL)
      if (nrow(found_row)>0){merged_car_df[rowid,"ENSEMBL"]=found_row$ENSEMBL[[1]]}
    }
  }
  
  # create final merged df
  merged_df=full_join(merged_car_df,merged_rna_df)
  merged_df=merged_df[,c("overlap_type","peakID","annotation","ENSEMBL","SYMBOL",
                         "log2FC_car","padj_car","log2FC_rna","padj_rna",
                         "geneChr","geneStart","geneEnd","gene_description")]
  fpath=paste0(output_dir,"overlap_",subset_type,"_",contrast_id_car,"_",contrast_id_rna,".csv")
  write.table(merged_df,fpath,sep=",")
  return(merged_df)
}

create_overlap_chrommap<-function(merged_df,subset_type){
  #http://bioconductor.org/packages/release/bioc/vignettes/karyoploteR/inst/doc/karyoploteR.html
  # subset and remove NA's
  sub_df=merged_df%>%tidyr::separate(peakID,c("seqnames","peakID"),sep=":")
  sub_df=subset(sub_df,seqnames %in% paste0("chr",c(1:22)))
  sub_df=sub_df%>% rename(c("Start"=geneStart,"End"=geneEnd))
  
  # create lists of up and down regulated
  both_up_df=subset(sub_df,log2FC_rna>log2fc_cutoff & log2FC_car>log2fc_cutoff)[,c("seqnames","Start","End")]
  both_down_df=subset(sub_df,log2FC_rna<log2fc_cutoff & log2FC_car<log2fc_cutoff)[,c("seqnames","Start","End")]
  car_up_df= subset(sub_df,log2FC_rna<log2fc_cutoff & log2FC_car>log2fc_cutoff)[,c("seqnames","Start","End")]
  car_down_df=subset(sub_df,log2FC_rna>log2fc_cutoff & log2FC_car<log2fc_cutoff)[,c("seqnames","Start","End")]
  
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
  mtext(paste0("Karyoplot of ", subset_type, " genes"),
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

main_differential_overlap<-function(subset_type){
  #peak db
  fpath=paste0(output_car_dir,"peak_annotation_",contrast_id_car,".csv")
  peak_df=read.csv(fpath)
  
  # annotation db
  fpath=paste0(output_car_dir,"DESeq2_res_",contrast_id_car,".csv")
  deseq_df=read.csv(fpath)
  deseq_df=rename(deseq_df,peakID="X")
  
  #merge peak and annotation
  car_df=full_join(peak_df,deseq_df,by="peakID")
  
  # remove any row without a gene symbol
  car_df=car_df[complete.cases(car_df$SYMBOL),]
  
  # rna db
  fpath=paste0(output_rna_dir,"DESeq2_",contrast_id_rna,"_DEG_allgenes_res1.txt")
  rna_df=read.csv(fpath,sep="\t")
  rna_df=separate(rna_df,"X",c("ENSEMBL","SYMBOL"),sep="[|]")
  
  # separate ensembl for merging
  rna_df=tidyr::separate(rna_df,ENSEMBL,c("ENSEMBL","ID"),sep="[.]")
  
  # save merged unfiltered
  unfiltered_merge=full_join(rna_df %>% rename(c("log2FC_rna"=log2FoldChange,"padj_rna"=padj)),
                             car_df %>% rename(c("log2FC_car"=log2FoldChange,"padj_car"=padj)))
  fpath=paste0(output_dir,"unfilt_overlap_",subset_type,"_",contrast_id_car,"_",contrast_id_rna,".csv")
  write.table(unfiltered_merge,fpath,sep=",")
  
  # reduce df to sig only, limit cols
  car_df_filt=subset(car_df,(padj<padj_cutoff) & (abs(log2FoldChange)>=log2fc_cutoff))
  rna_df_filt=subset(rna_df,padj<padj_cutoff & (abs(log2FoldChange)>=log2fc_cutoff))
  
  if (subset_type=="all"){
    # find overlaps and create venn diagrams and DT
    create_venn_diagrams(subset_type,rna_df_filt,car_df_filt)
    merged_df=create_overlapping_df(rna_df_filt,car_df_filt,subset_type)
    create_overlap_chrommap(merged_df,subset_type)
    create_overlap_DT(subset(merged_df,overlap_type=="overlap"))
  } else{
    # find all genes that are defined as the shortAnno type
    # filter the RNA df for these genes
    list_of_car_genes=subset(car_df,shortAnno==subset_type)$SYMBOL
    rna_df_filt2=subset(rna_df_filt,SYMBOL %in% list_of_car_genes)
    
    # of these genes, which are in the CAR df
    car_df_filt2=subset(car_df_filt,shortAnno==subset_type)
    
    # plot and DT
    create_venn_diagrams(subset_type,rna_df_filt2,car_df_filt2)
    merged_df=create_overlapping_df(rna_df_filt2,car_df_filt2,subset_type)
    create_overlap_chrommap(merged_df,subset_type)
    create_overlap_DT(subset(merged_df,overlap_type=="overlap"))
  }
}

############################################################
# run GSEA
############################################################
# create ranked df
prep_ranked_df<-function(subset_type,sig_type){
  # input dfs
  fpath=paste0(output_dir,"overlap_",subset_type,"_",contrast_id_car,"_",contrast_id_rna,".csv")
  overlap_df=read.csv(fpath)
  
  fpath=paste0(output_rna_dir,"DESeq2_",contrast_id_rna,"_DEG_allgenes_res1.txt")
  rna_df=read.csv(fpath,sep="\t")
  rna_df=rna_df %>% 
    separate("X",c("ENSEMBL","SYMBOL"),sep="[|]") %>% 
    separate("ENSEMBL",c("ENSEMBL","ID"),sep="[.]")
  
  # subset CAR for annotation of interest
  if (subset_type=="all"){
    sub_df=overlap_df
  } else{
    sub_df=subset(overlap_df,annotation %in% subset_type)
  }
  
  # sort and remove dups
  dedup_df=subset(sub_df,overlap_type=="overlap")
  if (sig_type=="both"){
    sig_RNA=subset(overlap_df,abs(log2FC_rna)>log2fc_cutoff & padj_rna<padj_cutoff)$ENSEMBL
    sig_CAR=subset(overlap_df,abs(log2FC_car)>log2fc_cutoff & padj_car<padj_cutoff)$ENSEMBL
    dedup_df=subset(dedup_df,ENSEMBL %in% intersect(sig_RNA,sig_CAR))
  }
  dedup_df=dedup_df[order(dedup_df$padj_rna),]
  dedup_df=dedup_df[!duplicated(dedup_df$ENSEMBL),]
  
  # pull significant genes ENSEMBL
  sorted_overlap_EID=dedup_df$ENSEMBL
  
  # remove all overlap EIDs and create list of RNA only
  sub_rna_df=rna_df[order(rna_df$padj),]
  sub_rna_df=sub_rna_df[!duplicated(sub_rna_df$ENSEMBL),]
  sub_rna_df=subset(sub_rna_df,ENSEMBL %ni% sorted_overlap_EID)
  sub_rna_df=sub_rna_df%>%rename("log2FC_rna"=log2FoldChange)
  sorted_rna_EID=sub_rna_df[order(sub_rna_df$log2FC_rna),]$ENSEMBL
  
  # create final df
  ranked_df=full_join(sub_rna_df,dedup_df) %>%
    subset(ENSEMBL %in% c(sorted_overlap_EID,sorted_rna_EID))
  rownames(ranked_df)=ranked_df$ENSEMBL
  ranked_df=ranked_df[c(sorted_overlap_EID,sorted_rna_EID),]
  
  return(ranked_df)
}

# set the annotation dbs
db_lookup<-function(t2g){
  # generate gene lists for C1
  # generate gene lists for C2 with subtypes biocarta, kegg, reactome, wiki
  # generate gene lists for C5 with subtypes MF, BP, CC
  # generate gene lists for Hallmark
  if (t2g=="C1"){
    db_out=msigdbr(species = species, category = "C1") %>% 
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

# run main GSEA function
main_gsea_function<-function(subset_type,sig_type,db_list){
  
  # prep ranked gene list
  ranked_df=prep_ranked_df(subset_type,sig_type)
  
  # create annodb top pathway df
  merged_df=data.frame()
  # pull db and generate list
  for (db_id in db_list){
    anno_db=db_lookup(db_id)
    msigdbr_list=split(anno_db$ensembl_gene,anno_db$gs_name)
    
    # create ranked input
    ranked_list=ranked_df$log2FC_rna
    names(ranked_list)=ranked_df$ENSEMBL  
    
    # run fgsea
    #https://cran.r-project.org/web/packages/msigdbr/vignettes/msigdbr-intro.html
    #https://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html
    fgseaRes <- fgsea(msigdbr_list, ranked_list,minSize  = minSize_gene_set)
    
    # plot 
    create_plot_fgsea(db_id,fgseaRes,ranked_list,msigdbr_list)
    
    #create merged top pathway df
    merged_df=rbind(merged_df,create_dts_fgsea(fgseaRes,db_id))
  }
  
  # print top pathway df
  caption_title=paste0("Top pathways across all annotation databases for ", subset_type, " genes.")
  DT::datatable(merged_df, extensions = 'Responsive', 
                      caption=htmltools::tags$caption(paste0(caption_title),
                                                      style="color:gray; font-size: 18px" ),
                rownames=F)
}

############################################################
# heatmaps for gene lists
############################################################
# create heatmaps for samples
generate_heatmaps_samples<-function(sample_id,project_id,gene_id){
  # read in counts matrix
  tmp_df=read.csv(paste0(parent_dir,
                         project_id,
                         "/results/peaks/contrasts/",
                         sample_id, "__dedup__norm.relaxed.bed/",
                         sample_id, "__dedup__norm.relaxed.bed_countsmatrix.txt"),sep="\t")
  
  # create peak list of sig PI genes
  peak_list=subset(pi_df,gene_list==gene_id)$peakID
  
  # subset for peaks in peak_list
  sub_df=subset(tmp_df,peakID %in% peak_list)
  
  # check for peaks, if none exit with message
  if (nrow(sub_df)<1){
    print(paste0("No peaks found within PI gene list for sample ",sample_id, " for ", gene_id))
  } else{
    # set peakID as rownmae, remove, and set df as numeric
    counts_in=sub_df[,c(2:ncol(sub_df))]
    counts_in <- sapply(counts_in, as.numeric)
    rownames(counts_in)=sub_df$peakID
    
    # transform and scale
    tmean.scale = t(scale(t(counts_in)))
    tmean.scale = tmean.scale[!is.infinite(rowSums(tmean.scale)),]
    tmean.scale = na.omit(tmean.scale)
    
    # Creating Dataframe to map samplenames to groups
    meta = subset(master_sample_df,sample_id %in% colnames(counts_in))
    groups <- data.frame(as.factor(meta$group_id))
    colnames(groups) <- "Groups"
    rownames(groups) <- meta$sample_id
    
    # Creating Group Column Annotation Colors
    columnColors <- c("lightpink","lightblue","orange","purple")
    names(columnColors) <- unique(groups$Groups)
    anno_colors <- list(Groups = columnColors)
    paletteLength <- 1000
    mycolors <- colorRampPalette(c("blue","white","red"), interpolate = "linear")(paletteLength)
    
    if (nrow(counts_in)>20){
      pheatmap(tmean.scale, 
               scale = "none", 
               main=paste0("Peaks associated with genes in ", gene_id,":\n",sample_id),
               cellwidth = 30, fontsize = 10, fontsize_row = 8, fontsize_col = 8, 
               color = mycolors, border_color = "NA",
               legend_breaks = c(-3,-2,-1,0,1,2,3), annotation_colors = anno_colors, 
               show_rownames = FALSE)
    } else{
      # generate annotation - SYMBOL,shortAnno (Ly6e,Promoter)
      anno_list=list()
      for (row_id in rownames(tmean.scale)){
        anno_list=append(anno_list,
                         paste(subset(pi_df,peakID==row_id)$SYMBOL,subset(pi_df,peakID==row_id)$shortAnno,sep="-"))
      }
      pheatmap(tmean.scale, 
               scale = "none", 
               main=paste0("Peaks associated with genes in ", gene_id,":\n",sample_id),
               cellwidth = 30, fontsize = 10, fontsize_row = 8, fontsize_col = 8, 
               color = mycolors, border_color = "NA",
               legend_breaks = c(-3,-2,-1,0,1,2,3), annotation_colors = anno_colors,
               labels_row = anno_list)
    }
  }
}

# create heatmaps for contrasts
generate_heatmaps_contrasts<-function(df_in,title_in){
  
  #subset for complete cases
  sub_df=df_in[complete.cases(df_in),]
  
  #for each sample log2foldchange values for each gene
  heatmap_df=data.frame()
  for(rowid in rownames(sub_df)){
    sample_id=sub_df[rowid,"sample"]
    gene_id=sub_df[rowid,"SYMBOL"]
    heatmap_df[sample_id,gene_id]=sub_df[rowid,"log2FoldChange"]
  }
  
  # shorten rownames
  rownames(heatmap_df)=gsub("_relaxed","",rownames(heatmap_df))
  rownames(heatmap_df)=gsub("_stringent","",rownames(heatmap_df))
  rownames(heatmap_df)=gsub("_narrowPeak","",rownames(heatmap_df))
  rownames(heatmap_df)=gsub("_dedup","",rownames(heatmap_df))
  rownames(heatmap_df)=gsub("_no","",rownames(heatmap_df))
  
  #set colors
  paletteLength <- 1000
  mycolors <- colorRampPalette(c("blue","white","red"), interpolate = "linear")(paletteLength)
  
  # scale
  scale_df= t(scale(t(heatmap_df)))
  scale_df=as.matrix(scale_df %>% replace(is.na(.), 0))
  
  if(ncol(scale_df)<90){
    pheatmap::pheatmap(scale_df, 
                       scale = "none", main=title_in,
                       fontsize = 10, fontsize_row = 6, fontsize_col = 4.5, color = mycolors, 
                       border_color = "NA",
                       legend_breaks = c(-5,-4,-3,-2,-1,0,1,2,3,4,5))
  }else{
    pheatmap::pheatmap(scale_df, 
                       scale = "none", main=title_in,
                       fontsize = 10, fontsize_row = 8, fontsize_col = 6, color = mycolors, 
                       border_color = "NA",
                       legend_breaks = c(-5,-4,-3,-2,-1,0,1,2,3,4,5),
                       show_colnames = FALSE)
  }
}
