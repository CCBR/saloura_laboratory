library(rlang, lib.loc="/usr/local/apps/R/4.2/site-library_4.2.0")


get_gene_body_boundary<-function(gene.symbol.list, promoter_span=2000)
{
   # gene.symbol.list<-unique(sort(c(goi.apm,goi.ifn)))
   library(biomaRt);
   ensembl <- useEnsembl(biomart = 'genes', dataset = "hsapiens_gene_ensembl");
   mart <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
   gene_info.0 <- biomaRt::getBM(attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position", "strand"),
   	filters = "hgnc_symbol",
   	values = gene.symbol.list,
   	mart = mart)
   gene_info<-gene_info.0[gene_info.0$chromosome_name %in% c(1:22,'X','Y'),];
   
      genes<-c();
      chr<-c();
      st<-c();
      ed<-c();
      strand<-c();
      st.bd<-c();
      ed.bd<-c();
   
   for(gg in unique( sort(gene_info$hgnc_symbol) ) )
   {
      gg.min.loc<-min(gene_info[gene_info$hgnc_symbol == gg,'start_position']);
      gg.max.loc<-min(gene_info[gene_info$hgnc_symbol == gg,'end_position']);
      gg.strand<-unique(gene_info[gene_info$hgnc_symbol == gg,'strand']);
      gg.chr<-unique(gene_info[gene_info$hgnc_symbol == gg,'chromosome_name']);
      if(length(gg.strand)>1){
          stop('Cannot continue without updating code.'); 
      }
      genes<-c(genes, gg);
      chr<-c(chr, gg.chr);
      st<-c(st, gg.min.loc);
      ed<-c(ed, gg.max.loc);
      strand<-c(strand, gg.strand);
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
   
   gene_bd<-data.frame(genes, chr, st, ed, strand, st.bd, ed.bd);
   rownames(gene_bd)<-genes;
   
   re<-{};
   re$version<-listDatasets(ensembl)[grep('hsapiens_gene_ensembl',listDatasets(ensembl)$dataset),];
   re$gene_info<-gene_info;
   re$gene_bd<-gene_bd;
   re$promoter_span<-promoter_span;
   re;
}


annotate_by_gb<-function(rdat.file, cnts.raw.master.annot, goi.ifn, goi.apm)
{
library(hash)
# load gene boundary info
load(rdat.file)

suppressPackageStartupMessages(library(GenomicRanges));
p<-cnts.raw.master.annot;
goi<-re$gene_bd;

      gr.p<-GRanges(seqnames=p$seqnames, ranges=IRanges(start=p$start, end=p$end, names=p$peakID));
      gr.goi<-GRanges(seqnames=paste0('chr',goi$chr), ranges=IRanges(start=goi$st.bd, end=goi$ed.bd, names=goi$genes));

      fo.goi.p<-GenomicRanges::findOverlaps(gr.goi, gr.p);
      goi.p.goiind<-data.frame(fo.goi.p)$queryHits;
      goi.p.pind<-data.frame(fo.goi.p)$subjectHits;
      matched.goi.p<-data.frame(goi.p.goiind, goi.p.pind);
      genes<-names(gr.goi)[matched.goi.p$goi.p.goiind];
      peaks<-names(gr.p)[matched.goi.p$goi.p.pind];
      p2g<-data.frame(peaks, genes);

    p2g.hs<-hash(keys=p2g$peaks, values=p2g$genes);
    for(pp in p2g$peaks[duplicated(p2g$peaks)]){
        .set(p2g.hs, keys=pp, values=unique(sort( c(p2g[p2g$peaks %in% pp, 'genes'], p2g.hs[[pp]])) ) );
    }
    msg<-NA;

    cnts.raw.master.annot$IFNalpha.flag<-'NO';
    cnts.raw.master.annot$APM.flag<-'NO';
    cnts.raw.master.annot$Immune.flag<-'NO';
    cnts.raw.master.annot$Immune.genes<-NA;
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

    anno.re<-{};
    anno.re$p2g.hs<-p2g.hs;
    anno.re$p2g<-p2g;
    anno.re$cnts.raw.master.annot<-cnts.raw.master.annot;
    anno.re;
}




# below is an example gene symbols. I passed all the IFN.alpha and apm gene together
imm.fname<-'/data/CUTRUN/docs/gene_list_master.csv';
goi.imm<-read.csv(imm.fname);
goi.ifn<-unique(sort(goi.imm$Human[goi.imm$Set=='IFNalpha']))
goi.apm<-unique(sort(goi.imm$Human[goi.imm$Set=='APM']))
goi.imm<-unique(sort(c(goi.ifn, goi.apm)));
gene.symbol.list<-goi.imm;

#------------
# retrieve gene boundary (including 2k promoter) 
# this needs to be done just once
# Since biomart is called, if the server is not responding on time, it can throw errors. Retry should work.
#------------
re<-get_gene_body_boundary(gene.symbol.list, promoter_span=2000);
save(re, file='ifn.app.annot.db.RData');

#------------
# Use annotated data as an input
#------------
# For example, on the script scripts/_diff_markdown.Rmd
# run upto the command below and use results_df as an input
# results_df = merge(results_df,pa,by=c("peakID"))
#---------------------
# save(results_df, file='/data/CUTRUN/tmp_scratch/tmp_results_df.RData')

load('/data/CUTRUN/tmp_scratch/tmp_results_df.RData')
#-----------

anno.re<-annotate_by_gb(rdat.file='ifn.app.annot.db.RData', results_df, goi.ifn, goi.apm)
annot2<-anno.re$cnts.raw.master.annot[,c('peakID','Immune.genes','Immune.flag','IFNalpha.flag','APM.flag')]

fname<-'53_H4K20me3_IFNb_vs_HN6_H4K20me3_IFNb__dedup__broadGo_peaks.bed_fragments.annot.v2.csv'
write.csv(annot2, file=fname, row.names=FALSE);

annot2.s<-annot2[!is.na(annot2$Immune.genes),];
fname<-'53_H4K20me3_IFNb_vs_HN6_H4K20me3_IFNb__dedup__broadGo_peaks.bed_fragments.annot.v2.immu.only.csv'
write.csv(annot2.s, file=fname, row.names=FALSE);

