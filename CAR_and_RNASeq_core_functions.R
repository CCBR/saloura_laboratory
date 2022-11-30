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
                      FCcutoff = 2^log2fc_cutoff_car,,
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
           xaxis=list(title="log2 Fold Change",
                      range =c(-5,5),
                      tickvals=c(-5,-4,-3,-2,-1,0,1,2,3,4,5),
                      ticktext=c('-32','-16','-8','-4','-2',
                                 '0','2','4','8','16','32')),
           yaxis=list(title=y_title,range =c(0,15)))
  p
  return(p)
}


############################################################
# run GSEA
############################################################
# create ranked df
prep_genelist_df<-function(subset_type,sig_type,analysis_type){
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
  
  return(output_df)
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

# create a small legend
addSmallLegend <- function(myPlot, pointSize = 2, textSize = 5, spaceLegend = 0.1) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_text(size = textSize), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
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

# save plots
print_save_plots<-function(plot_list,contrast_id,type_in){
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

# run ORA analysis
ora_plus_plot <- function(sigGeneList,db_id,contrast_in,n_show=5){
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
  ttl_abbrev=sub(" ","_",sub(":","_",db_id))
  fpath=paste0(output_dir,"ORA_",contrast_in[1],"-",contrast_in[2],"_table_",ttl_abbrev,".txt")
  write.table(resultdf,file=fpath,quote=FALSE,sep="\t",row.names = FALSE,col.names = TRUE)
  
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

# main function
main_gsea_ora_function<-function(subset_type,sig_type,db_list,analysis_type){
  
  # prep ranked gene list for GSEA or subset sig list for ORA
  subset_df=prep_genelist_df(subset_type,sig_type,analysis_type=analysis_type)
  
  # create annodb top pathway df
  merged_df=data.frame()
  
  # pull db and generate list
  for (db_id in db_list){
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
        o[[i]]=ora_plus_plot(sigGeneList=sigGeneList,db_id=db_id,
                             contrast_in=contrast_id_rna)
        i=i+1
      }
    }
  }
  
  # print,save ORA plots
  if (analysis_type=="ORA"){
    print_save_plots(o,contrast_id_rna,"ORA")
  }
  
  # print top pathway df
  caption_title=paste0("Top pathways across all annotation databases for ",
                       subset_type, " genes using ", analysis_type," analysis.")
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
