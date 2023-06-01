################################################################################
# formatting
################################################################################
# shorten colnmaes
shorten_names<-function(list_in){
  shortened_names=sub("", '',list_in)
  return(shortened_names)
}

############################################################
# QC Analysis
############################################################
format_counts_matrix<-function(contrast_id,sampleinfo,peak_type){
  # set extension
  extensions=c(paste0("__",dedup_status,"__",peak_type,".bed"))
  
  # set variables
  rawcountsmatrix=paste0(car_subpath,contrast_id,extensions,"/",
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

CREATE_PEAK_MAP<-function(contrast_id,condition1,condition2,peak_type){
  extensions=c(paste0("__",dedup_status,"__",peak_type,".bed"))
  fpath=paste0(car_subpath,contrast_id,extensions,"/peak_mapping/")
  
  # determine sample list
  sample_list=subset(groups_df, 
                     group%in%c(condition1,condition2))$sampleid
  
  # create df of peakIDs and source peak info from gopeaks
  peakid_mapped_df=data.frame()
  for (sid in sample_list){
    tmp_bed=read.csv(paste0(fpath,sid,"_peakmapped.bed"),sep="\t",header=FALSE)
    tmp_bed$peakid_contrast=paste0(tmp_bed$V4,":",tmp_bed$V5,"-",tmp_bed$V6)
    tmp_bed$peakid_sample=paste0(tmp_bed$V1,":",tmp_bed$V2,"-",tmp_bed$V3)
    
    # subset
    tmp_bed=tmp_bed[,c("peakid_contrast","peakid_sample")]
    
    # merge based on contrast ID's
    collapsed_bed=tmp_bed %>% 
      dplyr::group_by(peakid_contrast) %>% 
      dplyr::summarise(sid = toString(peakid_sample))
    
    # rename cols
    colnames(collapsed_bed)=c("peakid_contrast",sid)
    
    # merge
    if (nrow(peakid_mapped_df)==0){
      peakid_mapped_df=collapsed_bed
    } else{
      peakid_mapped_df=merge.data.frame(peakid_mapped_df, 
                                        collapsed_bed, 
                                        by="peakid_contrast",all=TRUE)
    }
    remove(tmp_bed,collapsed_bed)
  }
  
  # set NA to blank
  peakid_mapped_df[is.na(peakid_mapped_df)]=""
  
  return(peakid_mapped_df)
}

generate_RLE_plot<-function(condition1,condition2,input_data,input_title){
  #set colors
  colors <- brewer.pal(6, "Set2")
  x=shorten_sample_id(subset(groups_df,group==condition1 | group==condition2)$group)
  x=as.factor(x)#set pheno data
  
  #Plot results
  par(mfrow=c(1, 2), oma=c(3, 2, 0, 0)+0.1)
  plotRLE(as.matrix(input_data), outline=FALSE, ylim=c(-.5, .5), col=colors[x],las=2, cex.axis = .8)
  plotPCA(as.matrix(input_data), col=colors[x], cex=.8)
  mtext(input_title, side=2, outer=TRUE, adj=0)  
}

generate_boxplots<-function(sampleinfo,filtered,contrast_id){
  #set lib
  #determine lib reduction factor
  if (mean(colSums(filtered))>10000000){
    lib_factor=1e6
  } else if (mean(colSums(filtered))>1000000){
    lib_factor=1e5
  } else if (mean(colSums(filtered))>100000){
    lib_factor=1e4
  } else if (mean(colSums(filtered))>10000){
    lib_factor=1e3
  } else if (mean(colSums(filtered))>1000){
    lib_factor=1e2
  } else if (mean(colSums(filtered))>100){
    lib_factor=1e1
  } else {
    lib_factor=1e1
  }
  print(paste0("the libfactor is ",lib_factor," ",mean(colSums(filtered))))
  
  # filter
  sampleinfo=sampleinfo[sampleinfo$sampleid %in% colnames(filtered),]
  sampleinfo$library_size=colSums(filtered)/lib_factor
  sampleinfodf = as.data.frame(sampleinfo)
  sampleinfodf$dupstatus = dedup_status
  rownames(sampleinfo) = sampleinfo$sampleid
  pander(sampleinfodf,style="rmarkdown")
  
  file_name=paste0(output_dir,contrast_id,"_library_size.csv")
  write.csv(sampleinfodf,file_name)
  
  # melt data
  rawcounts_logcpm = log2(cpm(filtered))
  cpm_melt=reshape2::melt(rawcounts_logcpm)
  colnames(cpm_melt)=c("peakID","sampleid","log2cpm")
  
  # print boxplots
  p = ggplot(cpm_melt,aes(x=sampleid,y=log2cpm)) + 
    geom_boxplot(fill=as.factor(as.numeric(as.factor(sampleinfo$group))+1)) +
    theme_classic() +
    coord_flip()
  print(p)
}

run_deseq_analysis<-function(rawcounts,sampleinfo){
  # set variables
  bbpaths=paste0(car_subpath,"/bed_bedgraph_paths.tsv")
  
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
  
  # print scaling factors
  out_df=sf_df
  out_df=out_df[complete.cases(out_df),]
  rownames(out_df)=NULL
  out_df=subset(out_df,replicate %in% sampleinfo$sampleid)
  pander(out_df[,c("replicate","scalingfactor")],style="rmarkdown")
  
  # create list
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
      adjusted_sfs = sfs/mean(sfs)
      print("Samples are spiked and scaling factor used")
    }
    
    # AUC-based counts are prescaled, but fragmentbased counts are not prescaled
    if (rawcountsprescaled == "Y") {
      rawcounts=round(t(t(rawcounts) / sfs))
      dds <- DESeqDataSetFromMatrix(countData = as.matrix(rawcounts),
                                    colData = sampleinfo[,c("sampleid","group")],
                                    design = ~ group)
    }
    
    DESeq2::sizeFactors(dds)=1/adjusted_sfs
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
  
  return(pa)
}

generate_pca_plots<-function(dds,sampleinfo,exclusion_list,contrast_id){
  # default for nsub is 1000, if there are less than 1000 rows this will error
  if (nrow(dds) < 1000){
    rld <- vst(dds,nsub=nrow(dds))
  } else{
    rld <- vst(dds)
  }
  
  # analysis of variance
  assayrld = as.data.frame(assay(rld))
  assayrld$row_variance = rowVars(as.matrix(assayrld))
  assayrld = arrange(assayrld,desc(row_variance))
  zero_variance_rows=assayrld$row_variance<1e-5
  assayrld$row_variance = NULL
  assayrld = assayrld[!zero_variance_rows,]
  if (nrow(assayrld) > 500){
    assayrld=assayrld[1:500,]
  }
  
  # create title
  if (length(exclusion_list)==0){
    plottitle=paste0("All Samples Normalized\n",contrast_id)
  } else {
    plottitle=paste0("Selected Samples Normalized\n",contrast_id)
  }
  #plot PCA
  pca=prcomp(t(assayrld),scale. = T)
  m.pc1 = round(pca$sdev[1]^2/sum(pca$sdev^2)*100,2)
  m.pc2 = round(pca$sdev[2]^2/sum(pca$sdev^2)*100,2)
  m.pc3 = round(pca$sdev[3]^2/sum(pca$sdev^2)*100,2)
  xlab=paste0("PC1(",m.pc1,"%)")
  ylab=paste0("PC2(",m.pc2,"%)")
  p = ggplot(pca$x,aes(x=PC1,y=PC2,label=rownames(pca$x)))+geom_point(col=as.factor(as.numeric(as.factor(sampleinfo$group))+1))+
    xlab(xlab)+ylab(ylab)+
    ggtitle(plottitle)+
    geom_text_repel(max.overlaps = 10,size=2)+
    theme_light()
  print(p)
  
  # save final image
  fpath=paste0(img_dir,"pcoa_",gsub("-","",contrast_id),".png")
  ggsave(fpath,p)
}

generate_ecoli_plots<-function(contrast_id){
  bam_subpath=gsub("peaks/0.05/contrasts","bam",car_subpath)
  bam_subpath=gsub("peaks/contrasts","bam",bam_subpath)
  sample1=strsplit(contrast_id,"_vs_")[[1]][1]
  sample2=strsplit(contrast_id,"_vs_")[[1]][2]
  
  ecoli_df=data.frame()
  for (sampleid in subset(groups_df,group==sample1 | group==sample2)$sampleid){
    stats=read.table(paste0(bam_subpath,sampleid,".",dedup_status,".bam.idxstats"))
    stats=stats[,c("V1","V3")]
    colnames(stats)=c("location","read_count")
    stats$sampleid=sampleid
    stats$groupid=groups_df[groups_df$sampleid==sampleid,]$group
    
    if(nrow(ecoli_df)==0){
      ecoli_df=subset(stats,location=="NC_000913.3")
    } else{
      ecoli_df=rbind(subset(stats,location=="NC_000913.3"),
                     ecoli_df)
    }
  }
  
  p=ggplot(data=ecoli_df,aes(x=sampleid,y=read_count,fill=groupid)) + 
    geom_bar(stat="identity")
  p_final=p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    ggtitle(paste0("Spike-in control values\n", 
                   gsub("-","",contrast_id)))
  print(p_final)
  
  # save final image
  fpath=paste0(img_dir,"spikein_control_",gsub("-","",contrast_id),".png")
  ggsave(fpath,p_final)
}

############################################################
# Annotations
############################################################
GET_UPDATED_ANNOTATIONS_CHIPSEEK<-function(df_in){
  # df_in=anno_df
  
  # generate list of ENST ID's
  tmp_list=gsub("3' .*","",df_in$annotation)
  tmp_list=gsub("5' .*","",tmp_list)
  tmp_list=gsub("[D,P].*","",tmp_list)
  tmp_list=gsub("Intron [(]","",tmp_list)
  tmp_list=gsub("Exon [(]","",tmp_list)
  tmp_list=gsub("[/].*","",tmp_list)
  tmp_list=gsub("[.].*","",tmp_list)
  
  # reapply this list as a new col
  df_in$ensembl_transcript_id=tmp_list
  
  # create a unique ENST list to re-annotate
  tmp_list=tmp_list[tmp_list!=""]
  
  # pull annotation information
  mart=useEnsembl(biomart = "ensembl", 
                  dataset = "hsapiens_gene_ensembl", 
                  mirror = "useast")
  ensemble2gene <- getBM(attributes=c("ensembl_transcript_id",
                                      "external_gene_name"),
                         filters = "ensembl_transcript_id",
                         values = tmp_list, 
                         mart = mart)
  # merge with original data
  df_out=merge(df_in,ensemble2gene,by="ensembl_transcript_id",all=TRUE)
  
  # create new_symbol col replacing original symbol with the new symbol
  cleared_annotation_list=c("Distal","Promoter","3'UTR","5'UTR","Downstream")
  df_out$final_sym=df_out$SYMBOL
  df_out$final_sym[df_out$shortAnno %ni% cleared_annotation_list] <- df_out$external_gene_name[df_out$shortAnno %ni% cleared_annotation_list]
  df_out$final_sym[df_out$final_sym == ""] = NA
  
  # replace gene names that do not provide info
  remove_gene_list=c("Y_RNA")
  for (gid in remove_gene_list){
    df_out$final_sym[df_out$final_sym==gid]=NA
  }
  
  # replace column names
  colnames(df_out)=gsub("SYMBOL","old_SYMBOL",colnames(df_out))
  colnames(df_out)=gsub("final_sym","SYMBOL",colnames(df_out))
  
  return(df_out)
}

GET_GENE_BODY_BOUNDARY<-function(gene.symbol.list, promoter_span=2000){
  library(biomaRt);
  print("--prepping by gene boundry")
  ensembl <- useEnsembl(biomart = 'genes', dataset = "hsapiens_gene_ensembl")
  mart <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  gene_info.0 <- biomaRt::getBM(attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position", "strand"),
                                filters = "hgnc_symbol",
                                values = gene.symbol.list,
                                mart = mart)
  gene_info<-gene_info.0[gene_info.0$chromosome_name %in% c(1:22,'X','Y'),]
  
  genes<-c();
  chr<-c();
  st<-c();
  ed<-c();
  strand<-c();
  st.bd<-c();
  ed.bd<-c();
  
  for(gg in unique( sort(gene_info$hgnc_symbol))){
    gg.min.loc<-min(gene_info[gene_info$hgnc_symbol == gg,'start_position'])
    gg.max.loc<-min(gene_info[gene_info$hgnc_symbol == gg,'end_position'])
    gg.strand<-unique(gene_info[gene_info$hgnc_symbol == gg,'strand'])
    gg.chr<-unique(gene_info[gene_info$hgnc_symbol == gg,'chromosome_name'])
    
    if(length(gg.strand)>1){ stop('Cannot continue without updating code.')}
    
    genes<-c(genes, gg)
    chr<-c(chr, gg.chr)
    st<-c(st, gg.min.loc)
    ed<-c(ed, gg.max.loc)
    strand<-c(strand, gg.strand)
    
    if(gg.strand==-1){
      st.bd<-c(st.bd, gg.min.loc);
      ed.bd<-c(ed.bd, gg.max.loc+promoter_span);
    } else {
      if(gg.min.loc-promoter_span>0){
        st.bd<-c(st.bd, gg.min.loc-promoter_span);
      } else {
        st.bd<-c(st.bd, 1);
      }
      ed.bd<-c(ed.bd, gg.max.loc);
    }
  }
  
  gene_bd<-data.frame(genes, chr, st, ed, strand, st.bd, ed.bd)
  rownames(gene_bd)<-genes
  
  re<-{}
  re$version<-listDatasets(ensembl)[grep('hsapiens_gene_ensembl',
                                         listDatasets(ensembl)$dataset),]
  re$gene_info<-gene_info
  re$gene_bd<-gene_bd
  re$promoter_span<-promoter_span
  
  # save RE
  fpath=paste0(output_dir,"ifn.app.annot.db.RData")
  save(re, file=fpath)
  return(re)
}

ANNOTATE_BY_GB<-function(re,cnts.raw.master.annot, goi.ifn, goi.apm){
  library(hash)
  library(GenomicRanges)
  
  print("--reannotating")
  p<-cnts.raw.master.annot
  goi<-re$gene_bd
  
  # set GRANGES obj
  gr.p<-GRanges(seqnames=p$seqnames, 
                ranges=IRanges(start=p$start, end=p$end, names=p$peakID));
  gr.goi<-GRanges(seqnames=paste0('chr',goi$chr), 
                  ranges=IRanges(start=goi$st.bd, end=goi$ed.bd, names=goi$genes));
  
  # ifnd overlaps with new immune annotations
  fo.goi.p<-GenomicRanges::findOverlaps(gr.goi, gr.p)
  
  #sep into dfs
  goi.p.goiind<-data.frame(fo.goi.p)$queryHits
  goi.p.pind<-data.frame(fo.goi.p)$subjectHits
  matched.goi.p<-data.frame(goi.p.goiind, goi.p.pind)
  
  # ID matches
  genes<-names(gr.goi)[matched.goi.p$goi.p.goiind]
  peaks<-names(gr.p)[matched.goi.p$goi.p.pind]
  
  # create new df of matches
  p2g<-data.frame(peaks, genes)
  
  # create hash key and reset gene info
  p2g.hs<-hash(keys=p2g$peaks, values=p2g$genes)
  for(pp in p2g$peaks[duplicated(p2g$peaks)]){
    .set(p2g.hs, keys=pp, 
         values=unique(sort( c(p2g[p2g$peaks %in% pp, 'genes'], 
                               p2g.hs[[pp]])) ) )
  }

  # add flags
  cnts.raw.master.annot$IFNalpha.flag<-'NO'
  cnts.raw.master.annot$APM.flag<-'NO'
  cnts.raw.master.annot$Immune.flag<-'NO'
  cnts.raw.master.annot$Immune.genes<-NA
  for(kk in keys(p2g.hs)){
    cnts.raw.master.annot[kk, 'Immune.genes']<-paste(p2g.hs[[kk]], collapse=':');
    if (sum( p2g.hs[[kk]] %in% goi.apm) >0){
      cnts.raw.master.annot[kk,'APM.flag']<-'YES';
      cnts.raw.master.annot[kk,'Immune.flag']<-'YES';
    }
    if( sum( p2g.hs[[kk]] %in% goi.ifn) >0){
      cnts.raw.master.annot[kk,'IFNalpha.flag']<-'YES';
      cnts.raw.master.annot[kk,'Immune.flag']<-'YES';
    }
    if( sum( p2g.hs[[kk]] %in% goi.apm)==0 & sum( p2g.hs[[kk]] %in% goi.ifn)==0){
      stop('Cannot continue. Update code at annotate_by_gb.')
    }
  }
  
  # set output
  anno.re<-{}
  anno.re$p2g.hs<-p2g.hs
  anno.re$p2g<-p2g
  anno.re$cnts.raw.master.annot<-cnts.raw.master.annot
  
  return(anno.re)
}

GET_UPDATED_ANNOTATIONS_IMMUNE<-function(df_in){
  
  # set gene lists
  goi.ifn<-unique(sort(subset(gene_df,Set=='IFNalpha')$Human))
  goi.apm<-unique(sort(subset(gene_df,Set=='APM')$Human))
  gene.symbol.list<-unique(sort(c(goi.ifn, goi.apm)));
  
  #------------
  # retrieve gene boundary (including 2k promoter) 
  # save as RData file
  #------------
  re=GET_GENE_BODY_BOUNDARY(gene.symbol.list, promoter_span=2000);
  
  #------------
  # Use annotated data as an input
  #------------
  # For example, on the script scripts/_diff_markdown.Rmd
  # run upto the command below and use results_df as an input
  # results_df = merge(results_df,pa,by=c("peakID"))
  #---------------------
  # run reannotation
  anno.re<-ANNOTATE_BY_GB(re,df_in, goi.ifn, goi.apm)
  
  annot2<-anno.re$cnts.raw.master.annot[,c('peakID','Immune.genes',
                                           'Immune.flag','IFNalpha.flag','APM.flag')]
  
  # save files of all genes with immune notation, just immune notation
  fpath=paste0(data_dir,"_",contrast_id,"_SKanno.csv")
  write.csv(annot2, file=fpath, row.names=FALSE)
  
  fpath=paste0(data_dir,"_",contrast_id,"_SKanno_immune.csv")
  write.csv(annot2[!is.na(annot2$Immune.genes),], file=fpath, row.names=FALSE)
  
  return(anno.re)
}

############################################################
# MAIN FUNCTION for PREP
############################################################
main_prep_qc_core<-function(contrast_id,exclusion_list="",peak_type){
  
  #set conditions
  condition1=strsplit(contrast_id,"_vs_")[[1]][1]
  condition2=strsplit(contrast_id,"_vs_")[[1]][2]
  condition1=gsub("-","",condition1)
  condition2=gsub("-","",condition2)
  
  # filter based off of params
  sampleinfo=subset(groups_df,group==condition1 | group==condition2)
  rownames(sampleinfo)=NULL
  sampleinfo$group = relevel(as.factor(sampleinfo$group),condition2)
  
  # generate rawcounts
  rawcounts=format_counts_matrix(contrast_id,sampleinfo,peak_type)
  colnames(rawcounts)=gsub("-","",colnames(rawcounts))
  
  # filter raw counts
  raw_counts=ceiling(rawcounts)
  cpm_counts=edgeR::cpm(as.matrix(rawcounts))
  log_cpm_counts=log2(cpm_counts)
  keep=rowSums(cpm_counts>0.5)>1
  filtered=raw_counts[keep,]
  
  # shorten col names
  colnames(filtered)=shorten_sample_id(colnames(filtered))
  
  # save prenormalized peaks
  fpath=paste0(output_dir,"replicate_prenormalized_",gsub("-","",contrast_id),".csv")
  write.csv(filtered,fpath)
  
  # generate RLE plot
  generate_RLE_plot(condition1,condition2,rawcounts,
                    paste0("Fig. Before Normalization \n",
                           gsub("-","",contrast_id)))
  
  # generate boxplots
  generate_boxplots(sampleinfo,rawcounts,gsub("-","",contrast_id))
  
  # run deseq
  dds=run_deseq_analysis(filtered,sampleinfo)
  
  # read contrast df for annotation info
  extensions=c(paste0("__",dedup_status,"__",peak_type,".bed"))
  fpath=paste0(car_subpath,contrast_id,extensions,
               "/",contrast_id,extensions,"_",method,"based_diffresults.txt")
  anno_df=read.csv(fpath,sep = "\t")
  
  # update annotations from IMMUNE issue
  updated_anno=GET_UPDATED_ANNOTATIONS_IMMUNE(anno_df)
  
  # update annotations from CHIPSEEKER issue
  anno_df=GET_UPDATED_ANNOTATIONS_CHIPSEEK(anno_df)
  anno_df=anno_df[,c("peakID","shortAnno","SYMBOL")]
  head(anno_df)
  
  # annotate normalized peaks
  Normalized_counts_matrix<-as.data.frame(counts(dds, norm=TRUE))
  Normalized_counts_matrix$peakID=rownames(Normalized_counts_matrix)
  normalized_annotated=merge(Normalized_counts_matrix,
                             anno_df[c("peakID","shortAnno","SYMBOL")],
                             by="peakID")
  
  # generate peakcontrast to map gopeak source file
  peak_map_df=CREATE_PEAK_MAP(contrast_id,condition1,condition2,peak_type)
  fpath=paste0(output_dir,"peakmapping_",contrast_id,".txt")
  write.csv(peak_map_df,fpath)
  head(peak_map_df)
  
  # prepare condition list
  condition1=strsplit(contrast_id,"_vs_")[[1]][1]
  condition2=strsplit(contrast_id,"_vs_")[[1]][2]
  condition1=gsub("-","",condition1)
  condition2=gsub("-","",condition2)
  condition_list=c(condition1,condition2)
  
  # filter round1
  ## filter based on read count - ensure that at least the min sample thresholds
  ## have the reads expected
  for (contrastID in condition_list){
    sample_list=subset(groups_df,group==contrastID)$sampleid
    sample_list=gsub("-",".",sample_list)
    
    # determine number of reads for each peak; determine if there are more than 
    # threshold number of samples with that number of reads
    sample_count_threshold=round(length(sample_list)*sample_consensus_threshold+.5)-1
    normalized_annotated$read_threshold=rowSums(normalized_annotated[,sample_list] > (read_minimum_threshold-1))
    normalized_annotated$tmp_threshold=ifelse(normalized_annotated$read_threshold > sample_count_threshold,"Y","N")
    colnames(normalized_annotated)=gsub("tmp_threshold",paste0("sample_read_threshold_",contrastID),colnames(normalized_annotated))
    
    head(normalized_annotated)
  }
  
  # filter round 2
  ## determine if source GoPeaks is present in condition2 cells, at least in the
  ## minimum sample threshold
  peak_map_counts=data.frame()
  sample_list=subset(groups_df,group==condition2)$sampleid
  for (sid in sample_list){
    peak_map_sub=peak_map_df[,c("peakid_contrast",sid)]
    peak_map_sub$condition2=ifelse(peak_map_sub[,2] == "",0,1)
    peak_map_sub=peak_map_sub[,c("peakid_contrast","condition2")]
    colnames(peak_map_sub)=c("peakID",sid)
    
    if (nrow(peak_map_counts)==0){
      peak_map_counts=peak_map_sub
    } else{
      peak_map_counts=merge.data.frame(peak_map_counts,peak_map_sub,
                                       by="peakID",
                                       all=TRUE)
    }
    remove(peak_map_sub)
  }
  peak_map_counts$sampleSum=rowSums(peak_map_counts[,c(2:ncol(peak_map_counts))])
  head(peak_map_counts)
  
  # merge counts with norm data
  normalized_annotated=merge.data.frame(normalized_annotated, 
                                        peak_map_counts[,c(1,ncol(peak_map_counts))],
                                        by="peakID",
                                        all=TRUE)
  
  normalized_annotated[,paste0("sample_enrich_threshold_",condition2)]=ifelse(normalized_annotated$sampleSum > 
                                              (sample_count_threshold),"Y","N")
  normalized_annotated=normalized_annotated[!duplicated(normalized_annotated),]
  fpath=paste0(output_dir,"replicate_normalized_",gsub("-","",contrast_id),".csv")
  write.csv(normalized_annotated,fpath)

  # plot deseq2
  generate_RLE_plot(condition1,condition2,counts(dds, normalize=TRUE),
                    paste0("Fig. DESEq2 Normalization\n",
                           gsub("-","",contrast_id)))
  generate_pca_plots(dds,sampleinfo,exclusion_list,
                     gsub("-","",contrast_id))
  
  # plot ecoli
  generate_ecoli_plots(contrast_id)
  
  # run exclusions
  if (length(exclusion_list)!= 0){
    filtered=filtered[,colnames(filtered) %ni% exclusion_list]
    sampleinfo = subset(sampleinfo,!(sampleid %in% exclusion_list))
    
    # run analysis again
    dds=run_deseq_analysis(filtered,sampleinfo)
    generate_pca_plots(dds,sampleinfo,exclusion_list,contrast_id)
  }
  
  return(peakAnno)
}

main_prep_qc_rna_core<-function(flag_type){
  #load counts
  raw_counts= read.csv(paste0(input_dir,"DEG_ALL/RawCountFile_RSEM_genes.txt"),sep="\t")
  raw_counts=raw_counts[,c("symbol",groups_df$sampleid)]
  
  out_raw_counts=separate(raw_counts,col="symbol",into=c("ENSEMBL","SYMBOL"),sep="[|]")
  fpath=paste0(output_dir,"replicate_prenormalized_",csid,".csv")
  write.csv(out_raw_counts,fpath)
  
  ## Filter by CPM
  #CPM is calcualted as "how many counts would I get for a gene 
  #if the sample had a library size of 1M".
  raw_counts=column_to_rownames(raw_counts)
  raw_counts=ceiling(raw_counts)
  cpm_counts=edgeR::cpm(as.matrix(raw_counts))
  log_cpm_counts=log2(cpm_counts)
  
  sample_count_threshold=round(ncol(raw_counts)*sample_consensus_threshold+.5)-1
  keep=rowSums(cpm_counts>0.5)>sample_count_threshold
  
  filtered=raw_counts[keep,]
  colnames(filtered)=shorten_names(colnames(filtered))
  
  #set colors
  colors <- brewer.pal(6, "Set2")
  x=shorten_names(groups_df$group)
  x=as.factor(x)#set pheno data
  
  #merge into object
  set <- newSeqExpressionSet(as.matrix(filtered),
                             phenoData = data.frame(x, row.names=colnames(filtered)))
  set_u <- betweenLaneNormalization(set, which="upper")
  
  #Plot results
  if (flag_type=="pre_norm"){
    par(mfrow=c(1, 2), oma=c(3, 2, 0, 0)+0.1)
    plotRLE(as.matrix(filtered), outline=FALSE, ylim=c(-.5, .5), col=colors[x],las=2, cex.axis = .8)
    plotPCA(as.matrix(filtered), col=colors[x], cex=.8)
    mtext("Fig. Before Normalization", side=2, outer=TRUE, adj=0)
  }
  if (flag_type=="upper"){
    par(mfrow=c(1, 2), oma=c(3, 2, 0, 0)+0.1)
    plotRLE(set_u, outline=FALSE, ylim=c(-.5, .5), col=colors[x],las=2, cex.axis = .8)
    plotPCA(set_u, col=colors[x], cex=.8)
    mtext("Fig. UpperQuant Normalization", side=2, outer=TRUE, adj=0)
  }
  
  #run DESEQ2
  if (flag_type=="DESEQ2"){
    dds <- DESeqDataSetFromMatrix(countData = counts(set),
                                  colData = pData(set),
                                  design = ~ x)
    dds <- DESeq(dds)
    
    for (contrastID in unique(contrast_df$contrast)){
      # add sample filter
      condition1=strsplit(contrastID,"_vs_")[[1]][1]
      condition2=strsplit(contrastID,"_vs_")[[1]][2]
      condition_list=c(condition1,condition2)
      
      sample_list=c(subset(contrast_df,contrast==contrastID)$CNTRL,
                    subset(contrast_df,contrast==contrastID)$TREATMENT)
      sample_list=subset(groups_df,group %in% sample_list)$sampleid
      Normalized_counts_matrix<-as.data.frame(counts(dds, norm=TRUE))[,sample_list]
      
      for (conditionID in condition_list){
        sample_list=subset(groups_df,group==conditionID)$sampleid
        sample_list=gsub("-",".",sample_list)
        
        sample_count_threshold=round(length(sample_list)*sample_consensus_threshold+.5)-1
        Normalized_counts_matrix$read_threshold=rowSums(Normalized_counts_matrix[,sample_list] > (read_minimum_threshold-1))
        Normalized_counts_matrix$tmp_threshold=ifelse(Normalized_counts_matrix$read_threshold > sample_consensus_threshold,"Y","N")
        
        colnames(Normalized_counts_matrix)=gsub("tmp_threshold",
                                                paste0("sample_threshold_",conditionID),
                                                colnames(Normalized_counts_matrix))
      }
      head(Normalized_counts_matrix)
      
      # write out
      fpath=paste0(output_dir,"replicate_normalized_",contrastID,".csv")
      write.csv(Normalized_counts_matrix,fpath)
    }
    
    # plot
    par(mfrow=c(1, 2), oma=c(3, 2, 0, 0)+0.1)
    plotRLE(counts(dds, normalize=TRUE), outline=FALSE, ylim=c(-.5, .5), 
            col=colors[x],las=2, cex.axis = .8)
    plotPCA(counts(dds,normalize=TRUE), col=colors[x], cex=.8)
    mtext("Fig. DESEq2 Normalization", side=2, outer=TRUE, adj=0)
  }
}

# create summary graphics
chipseeker_plots<-function(peakAnno){
  #http://bioconductor.org/packages/devel/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html
  (plotAnnoPie(peakAnno))
  
  print(upsetplot(peakAnno, vennpie=TRUE))
  
  (plotDistToTSS(peakAnno,
                 title="Distribution of transcription factor-binding loci\nrelative to TSS"))
}

############################################################
# Create sig counts for contrasts
############################################################
# creates sig df for peaks / contrast
create_sig_contrast_df<-function(contrast_id,peak_type,gene_list_name="",rna_seq_path=""){
  # read in results
  ## these results are generated within the pipeline via DESEQ2
  extensions=c(paste0("__",dedup_status,"__",peak_type,".bed"))
  fpath=paste0(car_subpath,contrast_id,extensions,
               "/",contrast_id,extensions,"_",method,"based_diffresults.txt")
  contrast_df=read.csv(fpath,sep = "\t")[,c("peakID","log2FoldChange","padj","pvalue",
                                            "ENSEMBL")]
  contrast_id=gsub("-","",contrast_id)
  head(contrast_df)
  
  # pull in normalized read counts
  condition1=strsplit(contrast_id,"_vs_")[[1]][1]
  condition2=strsplit(contrast_id,"_vs_")[[1]][2]
  
  fpath=paste0(output_dir,"replicate_normalized_",contrast_id,".csv")
  norm_df=read.csv(fpath)
  colnames(norm_df)=gsub("X","",colnames(norm_df))
  colnames(norm_df)=gsub("[.]","_",colnames(norm_df))
  head(norm_df)
  
  norm_contrast_df=merge.data.frame(norm_df,contrast_df,by="peakID")
  head(norm_contrast_df)
  nrow(norm_df)
  nrow(norm_contrast_df)
  
  # perform manual change of annotations due to viewing on tracks
  ## only for project CS029758 with sample 53_H4K20me3_IFNb_vs_HN6_H4K20me3_IFNb
  if (cs_id=="CS029758" && contrast_id=="53_H4K20me3_IFNb_vs_HN6_H4K20me3_IFNb"){
    norm_contrast_df = fix_annotations_main(norm_contrast_df)
  }
  
  # if padj is NA set to 1
  norm_contrast_df$padj[is.na(norm_contrast_df$padj)]=1
  
  # create a sig_tracking df
  norm_contrast_df$flag_padj <- ifelse(norm_contrast_df$padj <=padj_cutoff,"Y", "N") # pval
  norm_contrast_df$flag_log2fc <- ifelse(abs(norm_contrast_df$log2FoldChange) >=log2fc_cutoff_car,"Y", "N") #log2fc
  norm_contrast_df$flag_anno <- ifelse(norm_contrast_df$shortAnno%in%gene_bodies_list,"Y", "N") #SYMBOL avail & not downstream
  norm_contrast_df$flag_padj_log2fc <- ifelse((norm_contrast_df$flag_padj == "Y") & (norm_contrast_df$flag_log2fc == "Y"),"Y", "N") #pval and log2fc
  norm_contrast_df$flag_padj_log2fc_anno <- ifelse((norm_contrast_df$flag_padj_log2fc == "Y") & (norm_contrast_df$flag_anno == "Y"),"Y", "N") #pval log2fc anno
  unique(norm_contrast_df$flag_padj_log2fc_anno)
  
  #add metadata
  norm_contrast_df$sample=contrast_id
  norm_contrast_df$dedup=strsplit(extensions,"__")[[1]][2]
  norm_contrast_df$type=strsplit(strsplit(extensions,"__")[[1]][3],"[.]")[[1]][2]
  norm_contrast_df$type=norm_contrast_df$type %>% replace(is.na(.),"narrowPeak")
  norm_contrast_df$method=method
  norm_contrast_df$uniqueid=paste0(norm_contrast_df$sample,"_",norm_contrast_df$dedup,"_",norm_contrast_df$type)
  
  # add RNA related data
  if (rna_seq_path != ""){
    rna_df=read.csv(rna_seq_path)[,c("SYMBOL","log2fc","fdr","significance_rna")]
    colnames(rna_df)=c("SYMBOL","log2FoldChange_rna","padj_rna","flag_rna")
    rna_df=subset(rna_df,flag_rna %in% c("Y","N"))
    head(rna_df)
    
    final_df=merge.data.frame(norm_contrast_df,rna_df,by="SYMBOL",all=TRUE)
    final_df=final_df[!is.na(final_df$flag_padj_log2fc_anno),]
    final_df=final_df[!duplicated(final_df$peakID),]
    
    final_df$mrna_dir=""
    final_df$mrna_dir[final_df$flag_rna=="N"]="neutral"
    final_df$mrna_dir[final_df$flag_rna=="Y" & final_df$log2FoldChange_rna>=0]="increase"
    final_df$mrna_dir[final_df$flag_rna=="Y" & final_df$log2FoldChange_rna<=0]="decrease"
    nrow(final_df)
  } else{
    final_df=norm_contrast_df[!duplicated(norm_contrast_df$peakID),]
  }
  
  # add IMMUNE related data
  final_df$immune_list="N"
  if (length(gene_list_name) >1){
    for (geneID in gene_list_name){
      final_df$immune_list[final_df$SYMBOL %in% subset(gene_df, Set==geneID)$Human]=geneID
    }
  }
  (unique(final_df$immune_list))
  
  # output sig table
  fpath=paste0(output_car_dir,"contrast_sig_",contrast_id,"_",gene_bodies_filename,".csv")
  write.csv(final_df,fpath)
  
  # filter results for significant values, print stats
  filt_df=subset(norm_contrast_df,flag_padj_log2fc_anno=="Y")
  if (nrow(filt_df)>0){
    print(paste0("----total number of significant peaks for contrast ", contrast_id,": ", nrow(filt_df)))
  } else{
    print(paste0("There are no significant peaks for contrast",contrast_id))
  }}

# create df for peaks / replicate
create_collapsed_sample_df<-function(contrast_id){
  #contrast_id="5-3_H3K4me3_IFNb_vs_HN6_H3K4me3_IFNb"
  
  # read in normalized counts matrix
  fpath=paste0(output_car_dir,"contrast_sig_",contrast_id,"_",gene_bodies_filename,".csv")
  peak_df=read.csv(fpath)
  colnames(peak_df)=gsub("X","",colnames(peak_df))
  head(peak_df)
  
  # separate contrast
  condition1=strsplit(contrast_id,"_vs_")[[1]][1]
  condition2=strsplit(contrast_id,"_vs_")[[1]][2]
  condition_list=c(condition1,condition2)
  
  # for each sample collapse annotation and sort by fold change
  collapsed_df=data.frame()
  for (conditionID in condition_list){
    # filter peaks based on thresholds
    sub_df=peak_df[peak_df[,paste0("sample_threshold_",conditionID)]=="Y",]
    if (nrow(sub_df)>0){
      sub_df=subset(sub_df,shortAnno %in% gene_bodies_list)[-1]
      
      # create total of peaks for calcs
      sub_df$total=nrow(sub_df)
      
      # collapse to get shortAnno counts
      collapsed_df=sub_df %>% dplyr::count(shortAnno,total)
      
      #calculate percentages
      collapsed_df$perc=round((collapsed_df$n/collapsed_df$total)*100,2)
    } else{
      collapsed_df=data.frame()
    }
    
    fpath=paste0(output_car_dir,"replicate_collapsed_annotated_",conditionID,"_",gene_bodies_filename,".csv")
    write.csv(collapsed_df,fpath)
  }
}

# collapses across all contrasts
create_collapsed_contrast_df<-function(contrast_id){
  
  # for each sample collapse annotation and sort by fold change
  fpath=paste0(output_car_dir,"contrast_sig_",contrast_id,"_",gene_bodies_filename,".csv")
  peak_df=read.csv(fpath)
  head(peak_df)
  
  # separate contrast
  condition2=strsplit(contrast_id,"_vs_")[[1]][2]
  
  # remove peaks where HN6 threshold is not met
  sub_df=peak_df[peak_df[,paste0("sample_threshold_",condition2)]=="Y",]
  head(sub_df)
  
  collapsed_df=data.frame()
  sig_value_list=c("Y","N")
  for (sig_value in sig_value_list){
    # if Y then create a sig list, otherwise create an ALL annotated peaks list
    if(sig_value=="N"){
      sub_df=subset(sub_df,flag_anno=="Y")
    } else{
      sub_df=subset(sub_df,flag_padj_log2fc_anno==sig_value)
    }
    
    # create total of peaks for calcs
    sub_df$total=nrow(sub_df)
    
    # collapse to get shortAnno counts
    collapsed_df=sub_df %>% dplyr::count(sample,shortAnno,dedup,type,method,total,uniqueid)
    
    # get counts for up/down
    collapsed_df$up=0
    collapsed_df$down=0
    tmp_direction1=(subset(sub_df,log2FoldChange>=0) %>%
                      dplyr::count(sample,shortAnno,dedup,type,method,total,uniqueid))[,c("sample","shortAnno","n")]
    rownames(tmp_direction1)=tmp_direction1$shortAnno
    tmp_direction2=(subset(sub_df,log2FoldChange<=0) %>%
                      dplyr::count(sample,shortAnno,dedup,type,method,total,uniqueid))[,c("shortAnno","n")]
    rownames(tmp_direction2)=tmp_direction2$shortAnno
    
    for (rowid2 in rownames(collapsed_df)){
      collapsed_df[rowid2,"up"]=as.numeric(tmp_direction1[collapsed_df[rowid2,"shortAnno"],"n"])
      collapsed_df[rowid2,"down"]=as.numeric(tmp_direction2[collapsed_df[rowid2,"shortAnno"],"n"])
    }
    collapsed_df[is.na(collapsed_df)] <- 0
    
    #calculate percentages
    collapsed_df$perc=round((collapsed_df$n/collapsed_df$total)*100,2)
    collapsed_df$perc_up=round((collapsed_df$up/sum(collapsed_df$up))*100,2)
    collapsed_df$perc_down=round((collapsed_df$down/sum(collapsed_df$down))*100,2)
    
    # write out df
    if (sig_value=="Y"){
      fpath=paste0(output_car_dir,"contrast_collapsed_sig_",contrast_id,"_",gene_bodies_filename,".csv")
    } else{
      fpath=paste0(output_car_dir,"contrast_collapsed_annotated_",contrast_id,"_",gene_bodies_filename,".csv")
    }
    write.csv(collapsed_df,fpath)
  }
}

# creates sig df for RNA / contrast
create_sig_RNA_df<-function(contrast_id,car_sig_path){
  ################################################
  # handle DESEQ counts
  ################################################
  # rna db
  fpath=paste0(input_dir,
               "DEG_",contrast_id,"_0.5_0.5",
               "/DESeq2_DEG_",contrast_id,"_all_genes.txt")
  deg_df=read.csv(fpath,sep="\t")
  
  # add sig of RNA
  deg_df$significance_rna="N"
  deg_df$significance_rna[abs(deg_df$log2fc)>log2fc_cutoff & deg_df$fdr<padj_cutoff]="Y"
  colnames(deg_df)=gsub("gene","SYMBOL",colnames(deg_df))
  head(deg_df)
  nrow(deg_df)
  unique(deg_df$significance_rna)
  
  ################################################
  # handle norm data
  ################################################
  fpath=paste0(output_dir,"replicate_normalized_",gsub("-","_vs_",contrast_id),".csv")
  norm_df=read.csv(fpath)
  norm_df=separate(norm_df,col="X",into=c("ENSEMBL","SYMBOL"),sep="[|]")
  
  # filter for threshold
  condition2=strsplit(contrast_id,"-")[[1]][2]
  norm_df=norm_df[norm_df[,paste0("sample_threshold_",condition2)]=="Y",]
  head(norm_df)
  
  # add filtered gene list to df for tracking
  failed_gene_df=subset(norm_df,SYMBOL %ni% deg_df$gene)[,c("ENSEMBL","SYMBOL")]
  col_list=c("fc","log2fc","pvalue","fdr","gsea_ranking_score")
  for (colid in col_list){
    failed_gene_df[,colid]<-0
  }
  failed_gene_df$significance_rna="F"
  colnames(failed_gene_df)=c("ensid_SYMBOL","SYMBOL",col_list,"significance_rna")
  head(failed_gene_df)
  unique(failed_gene_df$significance_rna)
  
  # merge deg_df and failed_gene_df
  rna_df=rbind(deg_df,failed_gene_df)
  head(rna_df)
  
  # merge norm_df and rna_df
  rna_df=merge.data.frame(norm_df,rna_df,by="SYMBOL")
  head(rna_df)
  
  ################################################
  # handle cut and run data
  ################################################
  # car db
  car_df=read.csv(car_sig_path)
  car_df=subset(car_df,flag_padj_log2fc_anno=="Y" & log2FoldChange<0 & shortAnno !="Distal")
  sig_gene_list=unique(car_df$SYMBOL)
  car_df=car_df[,c("SYMBOL","flag_padj_log2fc_anno")]
  car_df=car_df[!duplicated(car_df),]
  head(car_df)
  
  # add car data, clean
  final_df=merge.data.frame(rna_df,car_df,by="SYMBOL",all=TRUE)
  final_df=final_df[!is.na(final_df$read_threshold),]
  nrow(rna_df)
  nrow(final_df)
  head(final_df)
  
  # output files
  fpath=paste0(output_rna_dir,"contrast_sig_",contrast_id,".csv")
  write.csv(final_df,fpath)
}

############################################################
# pie charts for collapsed counts
############################################################
# create pie chart for each comparison
plot_pies_collapsed<-function(sub_in,df_in,y_in,percent_in,plot_in,title_in,fpath_in){
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
  
  # save and return
  ggsave(fpath_in,p_out)
  return(p_out)
}

generate_piecharts<-function(contrast_id){
  #set conditions
  condition1=strsplit(contrast_id,"_vs_")[[1]][1]
  condition2=strsplit(contrast_id,"_vs_")[[1]][2]
  condition1=gsub("-","_",condition1)
  condition2=gsub("-","_",condition2)
  
  plot_id=c("1.)","2.)")
  counter=1
  
  ##################### SAMPLE PEAKS
  # for each contrast, process samples separately
  for (conditionID in c(condition1,condition2)){
    # read in filtered collapsed,df
    fpath=paste0(output_car_dir,"replicate_collapsed_annotated_",conditionID,"_",gene_bodies_filename,".csv")
    if(file.size(fpath) > 3){
      collapsed_df=read.csv(fpath)
      # determine position for pie charts
      tmp_df <- collapsed_df %>% 
        mutate(csum = rev(cumsum(rev(n))), 
               pos = n/2 + lead(csum, 1),
               pos = if_else(is.na(pos), n/2, pos))
      tmp_df$pos=if_else(is.na(tmp_df$pos), tmp_df$n/2, tmp_df$pos)
      
      # plot
      fpath=paste0(img_dir,"sample_piechart_",conditionID,"_",gene_bodies_filename,".pdf")
      plot_title=paste0("Annotation of peaks for condition: ",conditionID,
                        " (N=",unique(tmp_df$total),")")
      p1=plot_pies_collapsed(collapsed_df,tmp_df,"n","perc",plot_id[counter],plot_title,fpath)
      print(p1)
      
      print(DT::datatable(tmp_df))
    }
    counter=counter+1
  }
  
  ##################### NONSIG+SIG PEAKS
  # for each contrast, process  nonsig differential peaks
  fpath=paste0(output_car_dir,"contrast_collapsed_annotated_",contrast_id,"_",gene_bodies_filename,".csv")
  collapsed_df=read.csv(fpath)
  head(collapsed_df)
  
  # determine position for pie charts
  tmp_df <- collapsed_df %>% 
    mutate(csum = rev(cumsum(rev(n))), 
           pos = n/2 + lead(csum, 1),
           pos = if_else(is.na(pos), n/2, pos))
  tmp_df$pos=if_else(is.na(tmp_df$pos), tmp_df$n/2, tmp_df$pos)
  
  # plot
  fpath=paste0(img_dir,"contrast_piechart_annotated_",contrast_id,"_",gene_bodies_filename,".pdf")
  plot_title=paste0("All Peaks by Annotation:\n",
                    contrast_id," (N=",sum(collapsed_df$n),")")
  p1 = plot_pies_collapsed(collapsed_df,tmp_df,"n","perc","3.",plot_title,fpath)
  print(p1)
  print(DT::datatable(tmp_df))
  
  ##################### SIG PEAKS
  # for each contrast, process  sig differential peaks
  fpath=paste0(output_car_dir,"contrast_collapsed_sig_",contrast_id,"_",gene_bodies_filename,".csv")
  collapsed_df=read.csv(fpath)
  head(collapsed_df)
  
  # determine position for pie charts
  tmp_df <- collapsed_df %>% 
    mutate(csum = rev(cumsum(rev(n))), 
           pos = n/2 + lead(csum, 1),
           pos = if_else(is.na(pos), n/2, pos))
  tmp_df$pos=if_else(is.na(tmp_df$pos), tmp_df$n/2, tmp_df$pos)
  
  ## all sig peaks
  fpath=paste0(img_dir,"contrast_piechart_sigpeaks_",contrast_id,"_",gene_bodies_filename,".pdf")
  plot_title=paste0("Significant Peaks by Annotation:\n",
                    contrast_id," (N=",sum(collapsed_df$n),")")
  p1 = plot_pies_collapsed(collapsed_df,tmp_df,"n","perc","4.",plot_title,fpath)
  print(p1)
  
  ## up
  tmp_df <- collapsed_df %>% 
    mutate(csum = rev(cumsum(rev(up))), 
           pos = up/2 + lead(csum, 1),
           pos = if_else(is.na(pos), up/2, pos))
  tmp_df$pos=if_else(is.na(tmp_df$pos), tmp_df$n/2, tmp_df$pos)
  
  fpath=paste0(img_dir,"contrast_piechart_sigpeaks_up_",contrast_id,"_",gene_bodies_filename,".pdf")
  plot_title=paste0(contrast_id,"\n",
                    "Significant Peaks by Annotation (N=",sum(collapsed_df$up),")\n",
                    "Increased in ", condition2)
  p2 = plot_pies_collapsed(collapsed_df,tmp_df,"up","perc_up","5.",plot_title,fpath)
  print(p2)
  
  ##down
  tmp_df <- collapsed_df %>% 
    mutate(csum = rev(cumsum(rev(down))), 
           pos = down/2 + lead(csum, 1),
           pos = if_else(is.na(pos), down/2, pos))
  tmp_df$pos=if_else(is.na(tmp_df$pos), tmp_df$n/2, tmp_df$pos)
  
  fpath==paste0(img_dir,"contrast_piechart_sigpeaks_down_",contrast_id,"_",gene_bodies_filename,".pdf")
  plot_title=paste0(contrast_id,"\n",
                    "Significant Peaks by Annotation (N=",sum(collapsed_df$down),")\n",
                    "Decreased in ", condition2)
  p3 = plot_pies_collapsed(collapsed_df,tmp_df,"down","perc_down","6.",plot_title,fpath)
  print(p3)
}

############################################################
# volcano plots for gene lists
############################################################
generate_volcano_plots<-function(data_type,contrast_id,gene_list_name=""){
  if(data_type=="CAR"){
    # read in sig df
    fpath=paste0(output_car_dir,"contrast_sig_",contrast_id,"_",gene_bodies_filename,".csv")
    results_df=read.csv(fpath,sep=",")
    
    # subset for gene_bodies_list
    results_df=subset(results_df,shortAnno %in% gene_bodies_list)
    
    # convert downstream to distal
    results_df$shortAnno=gsub("Downstream", "Distal", results_df$shortAnno)
    
    set_log=log2fc_cutoff_car
  } else{
    # rna db
    fpath=paste0(output_dir,"contrast_sig_",
                 contrast_id,".csv")
    results_df=read.csv(fpath)
    
    colnames(results_df)=gsub("log2fc","log2FoldChange",colnames(results_df))
    colnames(results_df)=gsub("fdr","padj",colnames(results_df))
    colnames(results_df)=gsub("gene","SYMBOL",colnames(results_df))
    (head(results_df))
    
    set_log=log2fc_cutoff
  }
  
  # calc log10
  results_df$log_pval=-log10(results_df$pvalue)
  
  # add colors and shapes
  #http://sape.inf.usi.ch/quick-reference/ggplot2/colour
  #https://www.datanovia.com/en/blog/ggplot-point-shapes-best-tips/
  colors=c("grey70","red","blue")
  shapes=c(20,2,3)
  
  #apply color based on list
  if (length(gene_list_name)==0){
    gene_list_name="allgenes"    
    anno_types=levels(as.factor(results_df$shortAnno))
    
    # set all vals to NS and grey
    keyvals=rep(colors[1],times=nrow(results_df))
    names(keyvals)=rep("NS",times=length(keyvals))
    
    for ( i in seq(1,length(anno_types))) {
      keyvals[ abs(results_df$log2FoldChange) > set_log & 
                 results_df$pvalue < padj_cutoff & 
                 results_df$shortAnno == anno_types[i]] = colors[i]
      names(keyvals)[keyvals == colors[i]] <- anno_types[i]
    }
    
    # colors are significance + annotation type
    p = EnhancedVolcano(results_df,
                        lab = results_df$SYMBOL,
                        x = 'log2FoldChange',
                        y = 'padj',
                        ylab = bquote(~-Log[10] ~ FDR),
                        pCutoff = padj_cutoff,
                        FCcutoff = set_log,
                        labSize = 4,
                        cutoffLineType = 'twodash',
                        cutoffLineWidth = 0.8,
                        title = paste0(contrast_id,"\n"),
                        subtitle = "",
                        subtitleLabSize = 1,
                        captionLabSize = 10,
                        colCustom = keyvals,
                        colAlpha = 1,
                        legendLabels = c("NS", expression(Log[2] ~ FC), "FDR", expression(FDR ~ and ~ log[2] ~ FC)),
                        legendLabSize = 10,legendPosition = 'right')
    print(p)
    
    
    # rename df for write out
    sub_df=results_df
    
  } else{
    # filter NA's,inf
    results_df=results_df[complete.cases(results_df$log_pval),]
    results_df=subset(results_df,log_pval!="Inf")
    
    # apply gene list types
    results_df$gene_annotation="Other Genes"
    for (gid in gene_list_name){
      gene_list=subset(gene_df,Set==gid)$Human
      results_df$gene_annotation[results_df$SYMBOL %in% gene_list]=gid
    }
    anno_types=levels(as.factor(results_df$gene_annotation))
    
    # create color/shape determination
    # for car data only look at the log2fc and pvalue
    # for RNA also consider if the CAR is sig, then review log2 and pvalue
    results_df$factor="Other Genes"
    if (data_type=="CAR"){
      for ( i in seq(1,length(anno_types))){
        if (anno_types[i]=="Other Genes"){
          next
        } else{results_df$factor[abs(results_df$log2FoldChange) >= set_log & 
                                   results_df$padj < padj_cutoff & 
                                   results_df$gene_annotation == anno_types[i]] <- paste0(anno_types[i]," Genes")
        }
      }
    } else{
      for ( i in seq(1,length(anno_types))){
        if (anno_types[i]=="Other Genes"){
          next
        } else{results_df$factor[abs(results_df$log2FoldChange) >= set_log & 
                                   results_df$padj < padj_cutoff & 
                                   results_df$significance_car == "Y" &
                                   results_df$gene_annotation == anno_types[i]] <- paste0(anno_types[i]," Genes")
        }
      }
    }
    
    #Set all distal to Other Genes
    results_df$factor[results_df$shortAnno=="Distal"]="Other Genes"
    
    # sort df
    results_df=results_df[order(results_df$gene_annotation,decreasing=TRUE),]
    
    # create plot
    p=ggplot(data=results_df, 
             aes(x=log2FoldChange, y=-log10(padj), col=factor, shape = factor)) + 
      geom_point(size = 3.5, stroke=1.5, alpha = 0.4) + 
      scale_color_manual("",
                         values=colors[1:length(unique(results_df$factor))],
                         breaks=c(unique(results_df$factor))) +
      scale_shape_manual("",
                         values=shapes[1:length(unique(results_df$factor))],
                         breaks=c(unique(results_df$factor))) +
      geom_vline(xintercept=c(-0.6, 0.6), col="black",linetype="dotted") +
      geom_hline(yintercept=-log10(0.05), col="black",linetype="dotted") +
      ggtitle(contrast_id) + 
      #scale_x_continuous(limits=c(-3,3)) +
      #scale_y_continuous(limits=c(0,80)) +
      xlab(bquote(Log[2]~fold~change))+
      ylab(bquote(-Log[10]~italic(padj)))
    
    p_final=p+theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank(),
                    legend.text=element_text(size=22),
                    axis.text = element_text(size=22, color="black"),
                    axis.title = element_text(size=22), 
                    axis.line = element_line(color = "black", linewidth=1),
                    legend.key=element_rect(fill="white"))
    print(p_final)
    
    # save plot
    fpath=paste0(img_dir,"contrast_volcano_",contrast_id,"_",gene_bodies_filename,".pdf")
    ggsave(fpath,p_final, width=3, height=3, units="in", scale=3)
    
    # add col for significance
    results_df$significance="Y"
    results_df$significance[results_df$factor=="Other Genes"]="N"
    
    # subset df for write
    if (data_type=="CAR"){
      sub_df=results_df[c("peakID","factor","log2FoldChange",
                          "padj","SYMBOL","shortAnno",
                          "gene_annotation","significance")]
      
    } else{
      sub_df=results_df[c("log2FoldChange","factor",
                          "padj","SYMBOL","gene_annotation")]
    }
  }
  
  # save DT
  sub_df$log2FoldChange=signif(sub_df$log2FoldChange,3)
  sub_df$padj=signif(sub_df$padj,3)
  merge_names=paste(gene_list_name,collapse = "_")
  fpath=paste0(output_dir,"volcano_data_",contrast_id,"_",data_type,"_",merge_names,".csv")
  print(fpath)
  write.csv(sub_df,fpath,row.names = FALSE)
  
}

############################################################
# heatmaps for sample replicates
############################################################
# Overwrites the pheatmap defaults
draw_colnames_45 <- function (coln, gaps, ...) {
  "Overwrites body of pheatmap:::draw_colnames, customizing it my liking"
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 25, gp = gpar(...))
  return(res)
}

# save pheatmap
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

# creates heatmap for multiple replicates
plot_heat_map<-function(df_in,show_names="ON",title_in="",cluster_by_rows="ON",fpath=""){
  #df_in=counts_matrix_complete;show_names="OFF";title_in="";cluster_by_rows="ON";fpath=fpath
  
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
    title_in=paste0("Significant Data Points (N=",nrow(df_in),")")
  }
  ####################
  # function
  ####################
  if (show_names=="OFF" && cluster_by_rows=="ON"){
    p=pheatmap(df_in, 
               scale = "none", main=title_in,
               cellwidth = 30, fontsize = 12, fontsize_row = 5, fontsize_col = 8, color = mycolors, 
               border_color = "NA",cluster_cols=F,annotation_colors = anno_colors, show_rownames = FALSE)
  } else if (show_names=="ON" && cluster_by_rows=="ON") {
    p=pheatmap(df_in, 
               scale = "none", main=title_in,
               cellwidth = 30, fontsize = 12, fontsize_row = 5, fontsize_col = 8, color = mycolors, 
               border_color = "NA",cluster_cols=F,annotation_colors = anno_colors, show_rownames = TRUE)
  } else if (show_names=="OFF" && cluster_by_rows=="OFF") {
    p=pheatmap(df_in, 
               scale = "none", main=title_in,
               cellwidth = 30, fontsize = 12, fontsize_row = 5, fontsize_col = 8, color = mycolors, 
               border_color = "NA",cluster_cols=F,cluster_rows=F,annotation_colors = anno_colors, 
               show_rownames = FALSE)
  }
  
  # set a fpath if it's missing
  if (fpath==""){
    fpath=paste0(img_dir,"replicate_heatmap_",contrast_id,".pdf")
  }
  
  save_pheatmap_pdf(p, fpath)
}

# fix scale_flag
generate_replicate_heatmaps<-function(contrast_id,scale_flag,gene_list_name,
                                      sample_subset="",rna_flag="N"){
  #contrast_id=contrast_id;scale_flag="ON";sample_subset=sample_sub_list;gene_list_name=""
  
  # split contrast
  contrast_1=strsplit(contrast_id,"_vs_")[[1]][1]
  contrast_2=strsplit(contrast_id,"_vs_")[[1]][2]
  
  # read in normalized counts matrix and save peaks df
  fpath=paste0(output_dir,
               "contrast_sig_",contrast_id,"_",gene_bodies_filename,".csv")
  counts_matrix=read.csv(fpath)
  colnames(counts_matrix)=gsub("X","",colnames(counts_matrix))
  print("**Pulling normalized data")
  print(paste0("--Total number of peaks: ", nrow(counts_matrix)))
  print(paste0("--Total number unique genes: ", length(unique(counts_matrix$SYMBOL))))
  
  # subset for only HN6
  hn6_only=counts_matrix[counts_matrix[,paste0("sample_threshold_",contrast_2)]=="Y",]
  print(paste0("**Subsetting for ", contrast_2))
  print(paste0("--Total number of peaks: ", nrow(hn6_only)))
  print(paste0("--Total number unique genes: ", length(unique(hn6_only$SYMBOL))))
  
  # subset without downstream
  no_downstream=subset(hn6_only,shortAnno!="Downstream")
  print(paste0("**Subsetting for ", contrast_2))
  print(paste0("--Total number of peaks: ", nrow(no_downstream)))
  print(paste0("--Total number unique genes: ", length(unique(no_downstream$SYMBOL))))
  
  # subset for sig peaks
  sig_peak_list=subset(no_downstream,flag_padj_log2fc_anno=="Y")$peakID
  sig_only=subset(no_downstream,peakID %in% sig_peak_list)
  print("**Subsetting for Sig genes")
  print(paste0("--Total number of peaks: ", nrow(sig_only)))
  print(paste0("--Total number unique genes: ", length(unique(sig_only$SYMBOL))))
  
  # subset for down
  down_only=subset(sig_only,log2FoldChange<0)
  print("**Subsetting for Down only")
  print(paste0("--Total number of peaks: ", nrow(down_only)))
  print(paste0("--Total number unique genes: ", length(unique(down_only$SYMBOL))))
  
  # subset for down, intragenic
  intragen_only=subset(down_only,shortAnno != "Distal")
  print("**Subsetting for intragenic genes")
  print(paste0("--Total number of peaks: ", nrow(intragen_only)))
  print(paste0("--Total number unique genes: ", length(unique(intragen_only$SYMBOL))))
  
  # subset for RNA genes
  if (rna_flag=="Y"){
    print("**Subsetting for RNA genes")
    rna_genes_list=subset(intragen_only,flag_rna %in% c("Y","N"))$SYMBOL
    
    # subset
    counts_matrix_subset=subset(intragen_only,SYMBOL %in% rna_genes_list)
    
    # print
    print(paste0("--Total number of peaks: ", nrow(counts_matrix_subset)))
    print(paste0("--Total number unique genes: ", length(unique(counts_matrix_subset$SYMBOL))))
  } else{
    counts_matrix_subset=intragen_only
  }
  
  # subset for immune gene list
  if (gene_list_name!=""){
    print("**Subsetting based on gene_list provided")
    
    # define genes
    if(gene_list_name=="APM_INFA"){
      gene_list_subset=subset(gene_df,Set %in% c("APM","IFNalpha"))$Human
    } else if (gene_list_name=="INFA"){
      gene_list_subset=subset(gene_df,Set %in% c("IFNalpha"))$Human
    } else if (gene_list_name=="APM"){
      gene_list_subset=subset(gene_df,Set %in% c("APM"))$Human
    }
    
    # subset
    counts_matrix_subset=subset(counts_matrix_subset,SYMBOL %in% gene_list_subset)
    
    # print
    print(paste0("--Total number of peaks: ", nrow(counts_matrix_subset)))
    print(paste0("--Total number unique genes: ", length(unique(counts_matrix_subset$SYMBOL))))
    
  } else{
    print("**Using all genes")
  }
  
  # subset for samples, if needed
  if (length(sample_subset)>1){
    print("**Subsetting based on sample_list provided")
    # create sample_list
    sample_list=sample_subset
  } else{
    print("**All samples are included")
    sample_list=subset(groups_df,group %in% c(contrast_1,contrast_2))$sampleid
  }
  
  # add rownames for immune print
  rownames(counts_matrix_subset)=make.unique(paste0(counts_matrix_subset$SYMBOL," (",counts_matrix_subset$shortAnno,")"))
  
  # prep df
  counts_matrix_subset=counts_matrix_subset[,sample_list]
  head(counts_matrix_subset)
  
  # scale if necessary
  if (scale_flag=="ON"){
    print("**Performing scaling")
    # ztransform df
    counts_matrix_complete=t(scale(t(counts_matrix_subset)))
    
    # fix any nan or inf
    counts_matrix_complete[is.nan(counts_matrix_complete)] <- 0
    counts_matrix_complete[counts_matrix_complete=="Inf"] <- 0
    range(counts_matrix_complete)
    
    # set fpath
    fpath=paste0(img_dir,"replicate_heatmap_withscale_")
    
  } else{
    print("**No scaling will be performed")
    counts_matrix_complete=counts_matrix_subset
    
    # set fpath
    fpath=paste0(img_dir,"replicate_heatmap_withoutscale_")
  }
  
  # shorten labels
  colnames(counts_matrix_complete)=gsub("_IFNb","",colnames(counts_matrix_complete))
  
  # generate heatmaps, showing names for gene_lists only
  print("**Generating heatmaps")
  if (gene_list_name!=""){
    # write out file
    fpath_f=paste0(fpath,gene_list_name,"_",contrast_id,"_",gene_bodies_filename,".csv")
    write.csv(counts_matrix_complete,fpath_f)
    
    fpath=paste0(fpath,gene_list_name,"_",contrast_id,"_",gene_bodies_filename,".pdf")
    
    # generate heatmap
    plot_heat_map(counts_matrix_complete,
                  show_names="ON",
                  title_in="",
                  cluster_by_rows="ON",
                  fpath=fpath)
  } else{
    # write out file
    fpath_f=paste0(fpath,gene_list_name,"_",contrast_id,".csv")
    write.csv(counts_matrix_complete,fpath_f)
    
    fpath=paste0(fpath,contrast_id,"_",gene_bodies_filename,".pdf")
    
    # generate heatmap
    plot_heat_map(counts_matrix_complete,
                  show_names="OFF",
                  title_in="",
                  cluster_by_rows="ON",
                  fpath=fpath)
  }
}

# creates heatmap for multiple replicates in RNASeq data
generate_replicate_heatmaps_rna<-function(contrast_id,scale_flag,gene_list_name,
                                          sample_subset="",car_sig_path){
  # split contrast
  contrast_1=strsplit(contrast_id,"-")[[1]][1]
  contrast_2=strsplit(contrast_id,"-")[[1]][2]
  
  # read in counts matrix
  fpath=paste0(output_rna_dir,"contrast_sig_",contrast_id,".csv")
  counts_matrix=read.csv(fpath,sep=",")
  head(counts_matrix)
  print(paste0("--Total number unique genes: ", length(unique(counts_matrix$SYMBOL))))
  
  # pull rowname
  rownames(counts_matrix)=make.unique(counts_matrix$SYMBOL)
  head(counts_matrix)
  
  # subset for sig gene list
  print("**Filter for CAR sig")
  counts_matrix_subset=subset(counts_matrix,flag_padj_log2fc_anno=="Y")
  print(paste0("--Total number unique genes: ", length(unique(counts_matrix_subset$SYMBOL))))
  
  # subset for pi gene list, if needed
  if (gene_list_name!=""){
    print("**Subsetting based on gene_list provided")
    
    # define genes
    if(gene_list_name=="APM_INFA"){
      gene_list_subset=subset(gene_df,Set %in% c("APM","IFNalpha"))$Human
    } else if (gene_list_name=="INFA"){
      gene_list_subset=subset(gene_df,Set %in% c("IFNalpha"))$Human
    } else if (gene_list_name=="APM"){
      gene_list_subset=subset(gene_df,Set %in% c("APM"))$Human
    }
    
    # subset for these genes
    counts_matrix_subset=subset(counts_matrix_subset,SYMBOL %in% gene_list_subset)
    print(paste0("--Total number unique genes: ", length(unique(counts_matrix_subset$SYMBOL))))
    
    # rename rows for printing
    rownames(counts_matrix_subset)=make.unique(paste0(counts_matrix_subset$SYMBOL))
  } else{
    print("**Using all genes")
    print(paste0("--Total number unique genes: ", length(unique(counts_matrix_subset$SYMBOL))))
  }
  
  # clean duplicate values
  out_dir=data.frame()
  for (symbolid in unique(counts_matrix_subset$SYMBOL)){
    # subset for SYMBOL
    sub_df=subset(counts_matrix_subset,SYMBOL==symbolid)
    
    # pull all sig flags
    sig_flags=sub_df$significance_rna
    
    # keep row based on order Y > N > F > M
    match_y=match(c("Y"),sig_flags)
    match_n=match(c("N"),sig_flags)
    match_f=match(c("F"),sig_flags)
    match_m=match(c("M"),sig_flags)
    
    if (!is.na(match_y)){
      out_dir=rbind(out_dir,sub_df[match_y,])
    } else if (!is.na(match_n)){
      out_dir=rbind(out_dir,sub_df[match_n,])
    }
  }
  nrow(out_dir)
  rownames(out_dir)=out_dir$SYMBOL
  counts_matrix_subset=out_dir
  
  # subset for samples
  if(length(sample_subset)>1){
    print("**Subsetting based on sample_list provided")
    # create sample_list
    sample_list=sample_subset
  } else{
    print("**All samples are included")
    sample_list=subset(groups_df,group %in% c(contrast_1,contrast_2))$sampleid
  }
  counts_matrix_subset=counts_matrix_subset[,sample_list]
  
  # clean df
  counts_matrix_subset=counts_matrix_subset[!duplicated(counts_matrix_subset),]
  nrow(counts_matrix_subset)
  
  # scale, if needed
  if (scale_flag=="ON"){
    print("performing scaling")
    
    # ztransform df
    counts_matrix_complete=t(scale(t(counts_matrix_subset)))
    
    # fix any nan or inf
    counts_matrix_complete[is.nan(counts_matrix_complete)] <- 0
    counts_matrix_complete[counts_matrix_complete=="Inf"] <- 0
    range(counts_matrix_complete)
    
    # set fpath
    fpath=paste0(img_dir,"replicate_heatmap_withscale_")
    
  } else{
    print("**No scaling will be performed")
    counts_matrix_complete=counts_matrix_subset
    
    # set fpath
    fpath=paste0(img_dir,"replicate_heatmap_withoutscale_")
  }
  
  print("**Generating heatmaps")
  if (gene_list_name!=""){
    # write out file
    fpath_f=paste0(fpath,gene_list_name,"_",contrast_id,"_",gene_bodies_filename,".csv")
    write.csv(counts_matrix_complete,fpath_f)
    
    # generate heatmap
    fpath=paste0(fpath,gene_list_name,"_",contrast_id,"_",gene_bodies_filename,".pdf")
    plot_heat_map(counts_matrix_complete,
                  show_names="ON",
                  title_in="",
                  cluster_by_rows="ON",
                  fpath=fpath)
  } else{
    # write out file
    fpath_f=paste0(fpath,contrast_id,"_",gene_bodies_filename,".csv")
    write.csv(counts_matrix_complete,fpath_f)
    
    # generate heatmap
    fpath=paste0(fpath,contrast_id,"_",gene_bodies_filename,".pdf")
    plot_heat_map(counts_matrix_complete,
                  show_names="OFF",
                  title_in="",
                  cluster_by_rows="ON",
                  fpath=fpath)
  }
}

############################################################
# chipintensity boxplots for sample replicates
############################################################
generate_intensity_boxplot<-function(contrast_id,sig_flag="",gene_list_name="",
                                     scale_factor="", rna_flag=""){
  #sig_flag="N";gene_list_name="";scale_factor="Y"
  
  # split contrast
  contrast_1=strsplit(contrast_id,"_vs_")[[1]][1]
  contrast_2=strsplit(contrast_id,"_vs_")[[1]][2]
  contrast_1=gsub("-","_",contrast_1)
  contrast_2=gsub("-","_",contrast_2)
  
  condition_list=c(contrast_1,contrast_2)
  sample_list=subset(groups_df,group%in% c(contrast_1,contrast_2))$sampleid
  sample_list=gsub("-",".",sample_list)
  
  # read in normalized counts matrix, contrast info and save peaks df
  fpath=paste0(output_dir,"contrast_sig_",contrast_id,"_",gene_bodies_filename,".csv")
  counts_matrix=read.csv(fpath)
  colnames(counts_matrix)=gsub("X","",colnames(counts_matrix))
  print(paste0("--Number of peaks included: ", nrow(counts_matrix)))
  
  # subset
  print("**Filtering for intragenic peaks")
  counts_matrix=subset(counts_matrix,shortAnno!="Distal")
  print(paste0("--Number of peaks included: ", nrow(counts_matrix)))
  
  #define rownames
  rownames(counts_matrix)=counts_matrix$peakID
  
  # filter for significant peaks
  if (sig_flag=="Y"){
    print("**Filtering for significant peaks")
    counts_matrix_subset=subset(counts_matrix,flag_padj_log2fc_anno=="Y")
    head(counts_matrix_subset)
    print(paste0("--Number of peaks included: ", nrow(counts_matrix_subset)))
  } else{
    counts_matrix_subset=counts_matrix
  } 
  
  # filter for immune genes
  if (gene_list_name=="APM_INFA"){
    print("**Filtering based on gene_list provided")
    
    # subset for these genes
    gene_list_subset=subset(gene_df,Set %in% c("APM","IFNalpha"))$Human
    counts_matrix_subset=subset(counts_matrix_subset,SYMBOL %in% gene_list_subset)
    print(paste0("--Number of peaks included: ", nrow(counts_matrix_subset)))
    
    # set fpath
    fpath_i=paste0(img_dir,"sample_boxplot_",gene_list_name,"_",
                   contrast_id,"_",gene_bodies_filename,".pdf")
    fpath_f=paste0(output_dir,"sample_boxplot_",gene_list_name,"_",
                   contrast_id,"_",gene_bodies_filename,".csv")
  } else{
    fpath_i=paste0(img_dir,"sample_boxplot_",
                   contrast_id,"_",gene_bodies_filename,".pdf")
    fpath_f=paste0(output_dir,"sample_boxplot_",
                   contrast_id,"_",gene_bodies_filename,".csv")
  }
  
  # filter for rna genes
  if (rna_flag=="Y"){
    print("**Filtering based on RNA provided")
    counts_matrix_subset=subset(counts_matrix_subset,flag_rna %in% c("Y","N"))
    print(paste0("--Number of peaks included: ", nrow(counts_matrix_subset)))
  }
  
  # create average counts by sample
  print("**Analyzing samples")
  for (conditionID in condition_list){
    # create sample_list
    sample_list=subset(groups_df,group%in% c(conditionID))$sampleid
    sample_list=gsub("-",".",sample_list)
    print(sample_list)
    # subset df
    tmp_df=counts_matrix_subset[,sample_list]
    
    # add one to all
    tmp_df=tmp_df+1
    
    # create mean col
    mean_col_name=paste0("mean_",conditionID)
    counts_matrix_subset[,mean_col_name]=rowMeans(tmp_df)
    
    # create log10 cols
    log_col_name=paste0("log10_",conditionID)
    counts_matrix_subset[,log_col_name]=log10(counts_matrix_subset[,mean_col_name])
  }
  head(counts_matrix_subset)
  
  # subset for mean or scale
  if (scale_factor=="Y"){
    print("**Scaling to log10 transformation")
    # create new df with log10 values
    mean_df=counts_matrix_subset[,grepl("log10",colnames(counts_matrix_subset))]
  } else{
    print("**Mean value is being used")
    # create new df with mean values
    mean_df=counts_matrix_subset[,grepl("mean",colnames(counts_matrix_subset))]
  }
  
  # save dt
  write.csv(mean_df,fpath_f)
  
  # pivot the table, convert to df, numeric
  df.long <- as.data.frame(pivot_longer(mean_df, cols=1:ncol(mean_df),
                                        names_to = "Contrast", values_to = "Value"))
  df.long$Contrast=gsub("log10_","",df.long$Contrast)
  df.long$Contrast=gsub("mean_","",df.long$Contrast)
  df.long$Value=as.numeric(df.long$Value)
  head(df.long)
  
  # update contrast name
  df.long$Contrast=gsub("_H4K20me3_IFNb","",df.long$Contrast)
  df.long$Contrast=gsub("_H3K4me3_IFNb","",df.long$Contrast)
  df.long$Contrast=gsub("_H3K9me3_IFNb","",df.long$Contrast)
  df.long$Contrast=gsub("53","SMYD3 KO",df.long$Contrast)
  df.long$Contrast=gsub("5-3","SMYD3 KO",df.long$Contrast)
  df.long$Contrast=gsub("_INFB_SMYD3","",df.long$Contrast)
  
  # swap order of factors
  short_contrast_id=unique(df.long$Contrast)
  df.long$Contrast=ordered(df.long$Contrast,
                           levels=c(short_contrast_id[2],short_contrast_id[1]))
  factor(df.long$Contrast)
  
  # define colors
  color_pal=c("black","blue")
  
  # plot boxplots
  p <- ggboxplot(df.long, x = "Contrast", y = "Value",
                 color = "Contrast", palette = color_pal,
                 short.panel.labs = FALSE)
  pf = p + 
    ylab("Log10(NormalizedTagCounts +1)") + 
    xlab("") + 
    ggtitle(paste0(contrast_id,"\nPeak Count: ",
                   nrow(counts_matrix_subset))) +
    theme(legend.position="top",
          legend.text=element_text(size=28),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_text(size=28, color="black"),
          axis.text.y = element_text(size=28, color="black")) +
    labs(color=NULL)+
    scale_y_continuous(limits=c(0,5)) +
    guides(color=guide_legend(nrow=2,byrow=TRUE)) +
    stat_compare_means(label = "p.format",method = "wilcox.test", paired=TRUE,
                       label.y.npc = "bottom",label.x.npc = "center",
                       size=10)
  print(pf)
  
  # save plot
  ggsave(fpath_i,pf, width=2, height=4, units="in", scale=3)
  
  # plot histograms
  mu <- ddply(df.long[,c(1,2)], "Contrast", summarise, grp.mean=mean(Value))
  p<-ggplot(df.long, aes(x=Value, color=Contrast)) +
    geom_histogram(fill="white", position="dodge")+
    geom_vline(data=mu, aes(xintercept=grp.mean, color=Contrast),
               linetype="dashed")
  #print(p)
}

CREATE_LONG_DF<-function(df_in){
  # pivot the table, convert to df, numeric
  df.long <- as.data.frame(pivot_longer(df_in, cols=1:ncol(df_in),
                                        names_to = "Contrast", values_to = "Value"))
  df.long$Contrast=gsub("log10_","",df.long$Contrast)
  df.long$Contrast=gsub("mean_","",df.long$Contrast)
  df.long$Value=as.numeric(df.long$Value)
  head(df.long)
  
  return(df.long)
}

generate_intensity_boxplot_two_projs<-function(contrast_id1,csid2,csid2_path,contrast_id2,gene_list_name=""){
  #contrast_id1="53_H4K20me3_IFNb_vs_HN6_H4K20me3_IFNb"; 
  #csid2_path="~/../../Volumes/CUTRUN/analysis/CS031308/r_analysis_230307/"
  #contrast_id2="53_INFB_SMYD3_vs_HN6_INFB_SMYD3"; csid2="CS031308";
  #gene_list_name="APM_INFA"
  
  # set input and output names
  if (gene_list_name==""){
    fpath_in1=paste0(output_dir,"sample_boxplot_",contrast_id1,"_",gene_bodies_filename,".csv")
    fpath_in2=paste0(csid2_path,"sample_boxplot_",contrast_id2,"_",gene_bodies_filename,".csv")
    fpath_out=paste0(img_dir,"sample_boxplot_",
                     cs_id,"_and_",csid2,"_",gene_bodies_filename,".pdf")
  } else{
    fpath_in1=paste0(output_dir,"sample_boxplot_",gene_list_name,"_",
                     contrast_id1,"_",gene_bodies_filename,".csv")
    fpath_in2=paste0(csid2_path,"sample_boxplot_",gene_list_name,"_",
                     contrast_id2,"_",gene_bodies_filename,".csv")
    fpath_out=paste0(img_dir,"sample_boxplot_",gene_list_name,"_",
                     cs_id,"_and_",csid2,"_",gene_bodies_filename,".pdf")
  }
  
  # read in first contrast
  boxdata1=read.csv(fpath_in1)[,c(2,3)]
  rownames(boxdata1)=NULL
  boxdata1=CREATE_LONG_DF(boxdata1)
  head(boxdata1)
  
  # add second contrast
  boxdata2=read.csv(fpath_in2)[,c(2,3)]
  rownames(boxdata2)=NULL
  boxdata2=CREATE_LONG_DF(boxdata2)
  head(boxdata2)
  
  df.long=rbind(boxdata1,boxdata2)
  head(df.long)
  
  # update contrast name
  df.long$Contrast=gsub("_IFNb","",df.long$Contrast)
  df.long$Contrast=gsub("_INFB","",df.long$Contrast)
  unique(df.long$Contrast)
  
  # swap order of factors
  short_contrast_id=unique(df.long$Contrast)
  df.long$Contrast=ordered(df.long$Contrast,
                           levels=c(short_contrast_id[2],short_contrast_id[1],
                                    short_contrast_id[4],short_contrast_id[3]))
  factor(df.long$Contrast)
  
  # define colors
  color_pal=c("black","blue","red","green")
  
  # plot boxplots
  p <- ggboxplot(df.long, x = "Contrast", y = "Value",
                 color = "Contrast", palette = color_pal,
                 short.panel.labs = FALSE)
  pf = p + 
    ylab("Chip Seq Intensity") + 
    xlab("") + 
    ggtitle(paste0(cs_id," and ", csid2)) +
    theme(legend.position="top",
          legend.text=element_text(size=28),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_text(size=28, color="black"),
          axis.text.y = element_text(size=28, color="black")) +
    labs(color=NULL)+
    scale_y_continuous(limits=c(0, round(max(df.long$Value)))) +
    guides(color=guide_legend(nrow=2,byrow=TRUE))
  print(pf)
  
  # save plot
  ggsave(fpath_out,pf, width=4, height=4, units="in", scale=3)
  
  # plot histograms
  mu <- ddply(df.long[,c(1,2)], "Contrast", summarise, grp.mean=mean(Value))
  p<-ggplot(df.long, aes(x=Value, color=Contrast)) +
    geom_histogram(fill="white", position="dodge")+
    geom_vline(data=mu, aes(xintercept=grp.mean, color=Contrast),
               linetype="dashed")
  #print(p)
}

############################################################
# venn diagrams
############################################################
# read in file, filter for sig genes, in gene list, if necessary
PREP_VENN_DIAGRAM<-function(fpath_in,sig_filter_list,gene_list_in,contrast_id){
  # fpath_in=fpath; gene_list_in=gene_list_to_include; sig_filter_list=sig_filter_list[2]
  
  # split contrast
  contrast_1=strsplit(contrast_id,"_vs_")[[1]][1]
  contrast_2=strsplit(contrast_id,"_vs_")[[1]][2]
  
  # read in sig df
  df_in=read.csv(fpath_in)
  print(paste0("--N genes: ", length(unique(df_in$SYMBOL))))
  
  # filter for sig
  if (sig_filter_list=="N"){
    df_sub=df_in
  } else{
    print("**Filtering for sig")
    df_sub=subset(df_in,flag_padj_log2fc_anno=="Y")
    print(paste0("--N genes: ", length(unique(df_sub$SYMBOL))))
  }
  
  # filter for gene list, if needed
  if (gene_list_in!=""){
    print("**Filtering for gene list")
    # gene list from pi df
    gene_list=subset(gene_df,Set %in% strsplit(gene_list_in,"_")[[1]])$Human
    
    # subset df
    df_sub=subset(df_sub,SYMBOL %in% gene_list)
    print(paste0("--N genes: ", length(unique(df_sub$SYMBOL))))
  }
  
  print("**Filtering for HN6 baseline")
  df_sub=df_sub[df_sub[,paste0("sample_threshold_",contrast_2)]=="Y",]
  print(paste0("--N genes: ", length(unique(df_sub$SYMBOL))))
  
  print("**Filtering for intragenic")
  df_sub=subset(df_sub,shortAnno != "Distal")
  print(paste0("--N genes: ", length(unique(df_sub$SYMBOL))))
  
  # create gene lists
  list_of_genes=unique(df_sub$SYMBOL)
  
  # remove NA"s
  list_of_genes=list_of_genes[!is.na(list_of_genes)]
  
  return(list_of_genes)
}

create_venn_diagrams_mult_projs<-function(cs_id_list,cs_id_path_list,contrast_list,sig_filter_list,gene_list_name="",gene_list_to_include=""){
  #cs_num=cs_id_list[1]; counter=1
  #cs_num=cs_id_list[2]; counter=2
  
  # create significant, unique gene lists
  counter=1; sample_id_merged=""
  for (cs_num in cs_id_list){
    # set contrast and path of sig file
    fpath=paste0(cs_id_path_list[counter],"contrast_sig_",
                 contrast_list[counter],"_",gene_bodies_filename,".csv")
    
    # split file to include 53 only
    var_name=strsplit(contrast_list[counter],"_vs_")[[1]][2]
    var_name=gsub("-","",var_name)
    var_name=gsub("_IFNb","",var_name)
    var_name=gsub("_INFB","",var_name)
    
    # create list, assign to sampleid name
    print(paste0("**Sample ", contrast_list[counter]))
    assign(paste0(var_name,"_list_of_genes"),
           PREP_VENN_DIAGRAM(fpath_in=fpath,
                             sig_filter_list=sig_filter_list[counter],
                             gene_list_in=gene_list_to_include,
                             contrast_id=contrast_list[counter]))
    
    # merge id for fpath out
    if (sample_id_merged==""){
      sample_id_merged=var_name
    } else{
      sample_id_merged=paste0(sample_id_merged,"_and_",var_name)
    }
    
    # reset counter
    counter=counter+1
  }
  
  # create merged list of all genes
  gene_list_variables=paste0(strsplit(sample_id_merged,"_and_")[[1]],"_list_of_genes")
  if (length(cs_id_list)==2){
    x <- list(A = get(gene_list_variables[1]), B = get(gene_list_variables[2]))
  } else if (length(cs_id_list)==3){
    x <- list(A = get(gene_list_variables[1]), B = get(gene_list_variables[2]), C=get(gene_list_variables[3]))
  }
  
  # prep sample names
  cat_name=strsplit(sample_id_merged,"_and_")[[1]]
  sig_filter_list=gsub("N","ALL",sig_filter_list)
  sig_filter_list=gsub("Y","SIG ONLY",sig_filter_list)
  
  # prep output df
  ## create list of unique genes
  tmp_list=""
  for (i in (1:length(x))){
    tmp_list=unique(c(tmp_list,x[[i]]))
  }
  tmp_list=tmp_list[tmp_list != ""]
  print(paste0("Union of genes found in all three samples: ",length(tmp_list)))
  ## create df matching genes to samples
  out_df=data.frame(tmp_list=tmp_list)
  for (i in (1:length(x))){
    out_df[,cat_name[[i]]] <- with(out_df, ifelse(tmp_list %in% x[[i]], 'Y', 'N'))
  }
  ## create merged df
  for (rowid in rownames(out_df)){
    rowcount=rowSums(out_df[rowid,2:(length(x)+1)] == "Y")
    if (rowcount==3){
      flag="ALL"
    } else if (rowcount==2){
      if(out_df[rowid,2]=="Y" & out_df[rowid,3]=="Y"){flag=paste0(colnames(out_df)[2],"_AND_",colnames(out_df)[2])}
      if(length(x)>2){
        if(out_df[rowid,2]=="Y" & out_df[rowid,4]=="Y"){flag=paste0(colnames(out_df)[2],"_AND_",colnames(out_df)[4])}
        if(out_df[rowid,3]=="Y" & out_df[rowid,4]=="Y"){flag=paste0(colnames(out_df)[3],"_AND_",colnames(out_df)[4])}
      }
    } else if (rowcount==1){
      if(out_df[rowid,2]=="Y"){flag=colnames(out_df)[2]}
      if(out_df[rowid,3]=="Y"){flag=colnames(out_df)[3]}
      if(length(x)>2){
        if(out_df[rowid,4]=="Y"){flag=colnames(out_df)[4]}
      }
    }
    out_df[rowid,"overlap_type"]=flag
  }
  
  # prep venn diagram labels
  counter=1
  for (catid in cat_name){
    cat_name=gsub(cat_name[counter],paste0(cat_name[counter],"\n(",sig_filter_list[counter],")"),cat_name)
    counter=counter+1
  }
  
  # prep output file name
  if (gene_list_name==""){
    fpath_out=paste0(img_dir,"venndiagram_",
                     sample_id_merged,"_",gene_bodies_filename,".pdf")
    dfpath_out=paste0(img_dir,"venndiagram_",
                      sample_id_merged,"_",gene_bodies_filename,".csv")
    full_title=paste0("Baseline Genes")
  } else{
    fpath_out=paste0(img_dir,"venndiagram_",
                     sample_id_merged,"_",gene_list_name,"_",gene_bodies_filename,".pdf")
    dfpath_out=paste0(img_dir,"venndiagram_",
                      sample_id_merged,"_",gene_list_name,"_",gene_bodies_filename,".csv")
    full_title=paste0("Baseline Genes in ",gene_list_name, " genes")
  }
  
  # write out df
  write.csv(out_df,dfpath_out)
  
  # Venn diagram with custom category names
  p = ggVennDiagram(x, color = 1, lwd = 0.7,
                    category.names = cat_name,
                    label_size = 4) + 
    scale_fill_gradient(low = "blue", high = "red")
  pf = p + ggtitle(full_title) + scale_x_continuous(expand = expansion(mult = .2))
  pf
  
  # save and print
  print(pf)
  ggsave(fpath_out,pf)
}
