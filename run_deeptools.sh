
##################################################################
# docs
##################################################################
# deep tools
# https://hbctraining.github.io/Intro-to-ChIPseq/lessons/10_data_visualization.html
# https://deeptools.readthedocs.io/en/develop/content/tools/computeMatrix.html
##################################################################
cs_id=$1
sampleid=$2
gene_shorthand=$3
min_max=$4

# center, TSS
center_type="center"

# load modules
module load deeptools bedtools

# Example Usage
# sh run_deeptools.sh CS029758 NA genomic_regions 15000
# sh run_deeptools.sh CS029758 NA immune_genes 15000
##################################################################
# Set samples
##################################################################
if [[ $cs_id == "CS029758" ]]; then
	carlisle_id="carlisle_230111"
    deep_id="deeptools_230317"
	deeptools_genomic="/data/CUTRUN/analysis/CS029758/r_analysis_230316/contrast_deeptools_rna_53_H4K20me3_IFNb_vs_HN6_H4K20me3_IFNb_PIE53.csv"
	deeptools_immunogenic="/data/CUTRUN/analysis/CS029758/r_analysis_230316/contrast_deeptools_rna_pi_53_H4K20me3_IFNb_vs_HN6_H4K20me3_IFNb_PIE53.csv"
	
    f_list=("HN6_H4K20me3_IFNb_1" "HN6_H4K20me3_IFNb_2" "HN6_H4K20me3_IFNb_3" "53_H4K20me3_IFNb_1" "53_H4K20me3_IFNb_2" "53_H4K20me3_IFNb_3" )
	l_list=("HN6_Rep1" "HN6_Rep2" "HN6_Rep3" "5-3_Rep1" "5-3_Rep2" "5-3_Rep3")
	type="broadGo"
	f_id="53_H4K20me3_IFNb_vs_HN6_H4K20me3_IFNb"
	f_id1="HN6_H4K20me3_IFNb"
	f_id2="53_H4K20me3_IFNb"
fi

# check gene_list
## gene lists allow for subsetting the heatmaps by specific locations, gene annotation type, etc
if [[ $gene_shorthand == "genomic_regions" ]]; then
	deeptools_input=$deeptools_genomic
	echo "Running genomic regions"

elif [[ $gene_shorthand == "immune_genes" ]]; then
	deeptools_input=$deeptools_immunogenic
	echo "Running immune genes"
else
	echo "Check shorthand"
	exit
fi
##################################################################
# set dirs
##################################################################
code_dir="/data/CUTRUN/analysis/github"
ref_dir="/data/CUTRUN/analysis/$cs_id/$deep_id/${f_id}/refs"
deep_dir="/data/CUTRUN/analysis/$cs_id/$deep_id/${f_id}/$min_max/$gene_shorthand"
peak_dir="/data/CUTRUN/analysis/$cs_id/$carlisle_id/results/peaks/gopeaks"
echo "$deep_dir"

project_dir="/data/CUTRUN/analysis/$cs_id/$carlisle_id/results"
bigwig_dir="$project_dir/bigwig"

# make dirs
dir_list=($ref_dir/tmp $deep_dir)
for newdir in ${dir_list[@]}; do
	if [[ ! -d $newdir ]]; then mkdir -p $newdir; fi
done

dir_list=(sh logs plots gz bed summary counts)
for newdir in ${dir_list[@]}; do
	if [[ ! -d $deep_dir/$newdir ]]; then mkdir -p $deep_dir/$newdir; fi
done
##################################################################
# create gene GTF needed for all anlysis
##################################################################
# create bed files  - genes only
## source gtf - complete annotation
## genes_bed is genes only
## target is in corrected format
source_gtf="/data/CCBR_Pipeliner/db/PipeDB/Indices/GTFs/hg38/gencode.v30.annotation.gtf"
subset_bed="/data/CCBR_Pipeliner/Pipelines/CARLISLE_dev/resources/gencode_genes.v30.bed"
gene_bed="/data/CCBR_Pipeliner/Pipelines/CARLISLE_dev/resources/gencode.v30.bed"
if [[ ! -f $gene_bed ]]; then
	echo "creating DEEPTOOLS bedfile"
    
    # subset for only genes, format
    awk '$3 == "gene" { print $0 }' $source_gtf > $subset_bed
	awk '{ print $1"\t"$4"\t"$5 }' $subset_bed > $gene_bed
	
    echo "Bed here: $gene_bed"
fi

# create bed file specifc to gene list
genelist_specific_bed="$ref_dir/${f_id}_${gene_shorthand}.bed"
if [[ ! -f $genelist_specific_bed ]]; then

	echo "creating DEEPTOOLS filtered bedfile"

    # create loc only bedfile
    awk '{ print $1"\t"$4"\t"$5 }' $subset_bed > tmp.bed

    # convert output of R to bed file
    awk -F"," '{ print $2"\t"$3"\t"$4 }' $deeptools_input | tail -n +2 | sed -s "s/\"//g" > tmp_list.bed

    # interset for gene list
    bedtools intersect -wa -a tmp.bed -b tmp_list.bed > $genelist_specific_bed

    # preview and remove tmps
	echo "Bed here: $genelist_specific_bed"
    head $genelist_specific_bed
	rm tmp.bed tmp_list.bed
fi

# for center algorithm, create genelist_specific_bed relative to the sample
sample_genelist_specific_bed="$ref_dir/${f_id}_${gene_shorthand}_samplespecific.bed"
if [[ $center_type == "center" && ! -f $sample_genelist_specific_bed ]]; then
	
    # for each sample, create bed file of locations of interst only
    # for genomic samples this will include significant differential locations, intragenic
    # for immune samples this will include the above set, subset for only immune-related genes
	for s_id in ${f_list[@]}; do
		peak_bed="$peak_dir/${s_id}*broadGo_peaks.bed"
		sample_bed="$ref_dir/tmp/${s_id}_broadGo_peaks_genes.bed"

		if [[ ! -f $sample_bed ]]; then
			echo "----subsetting bed $s_id"
			bedtools intersect -a $peak_bed -b $genelist_specific_bed > $sample_bed
		fi
	done

	# merge bedfile together to create one contrast specific bed file
	cat $ref_dir/tmp/* > $ref_dir/tmp/allsamples_genesonly.bed

	# remove duplicated values
	sort $ref_dir/tmp/allsamples_genesonly.bed | uniq > $ref_dir/tmp/allsamples_genesonly_unique.bed

    # # sort
    bedtools sort -i $ref_dir/tmp/allsamples_genesonly_unique.bed > $sample_genelist_specific_bed

	# #cleanup
    echo "Bed here: $sample_genelist_specific_bed"
	head $sample_genelist_specific_bed
	rm $ref_dir/tmp/* 
fi

##################################################################
# Compute Matrix
##################################################################
submit_batch(){
	if [[ $gene_shorthand == "immune_genes" ]]; then
		sbatch --cpus-per-task=8 --verbose \
		--output=$deep_dir/sh/%j.out \
		--mem=2g --gres=lscratch:450 --time 00:30:00 \
		--error=$deep_dir/sh/%j.err $sh
	else
		sbatch --cpus-per-task=32 --verbose \
		--output=$deep_dir/sh/%j.out \
		--mem=30g --gres=lscratch:450 --time 2:00:00 \
		--error=$deep_dir/sh/%j.err $sh
	fi
}
sample_genelist_specific_bed="$ref_dir/${f_id}_${gene_shorthand}_samplespecific.bed"
# if [[ -f $bigwig_dir/${f_id1}_1.dedup.bigwig && ! -f $deep_dir/gz/${f_id}__dedup__${type}.gz ]]; then
# 	echo "--matrix $f_id full"
# 	sh="$deep_dir/sh/matrix_${f_id}.sh"

# 	echo "#!/bin/sh
# 	module load deeptools
# 	computeMatrix reference-point \
# 		--referencePoint $center_type -b $min_max -a $min_max -p 6 \
# 		-R $sample_genelist_specific_bed \
# 		-S $bigwig_dir/${f_id1}*.dedup.bigwig $bigwig_dir/${f_id2}*.dedup.bigwig  \
# 		-o $deep_dir/gz/${f_id}__dedup__${type}.gz \
# 		--missingDataAsZero --skipZeros --sortRegions descend \
# 		--outFileSortedRegions $deep_dir/bed/${f_id}_regions_$center_type__dedup__${type}.bed"  > $sh
# 	submit_batch $sh
# fi

if [[ -f $bigwig_dir/${f_id1}_1.dedup.bigwig && ! -f $deep_dir/gz/${f_id}__dedup__${type}_single.gz ]]; then
	echo "--matrix $f_id single"
	sh="$deep_dir/sh/matrix_${f_id}_single.sh"
		
	echo "#!/bin/sh
	module load deeptools
	computeMatrix reference-point \
		--referencePoint $center_type -b $min_max -a $min_max -p 6 \
		-R $sample_genelist_specific_bed \
		-S $bigwig_dir/${f_id1}_1.dedup.bigwig $bigwig_dir/${f_id2}_1.dedup.bigwig  \
		-o $deep_dir/gz/${f_id}__dedup__${type}_single.gz \
		--missingDataAsZero --skipZeros --sortRegions descend \
		--outFileSortedRegions $deep_dir/bed/${f_id}_regions_$center_type__dedup__${type}_single.bed"  > $sh
	# submit_batch $sh
	sh $sh
fi


##################################################################
# Summary
##################################################################
if [[ -f $deep_dir/gz/${f_id}__dedup__${type}.gz && ! -f $deep_dir/summary/${f_id}_scores_per_bin.npz ]]; then
	echo "--summary $f_id full"
	sh="$deep_dir/sh/summary_${f_id}.sh"
			
	echo "multiBigwigSummary bins \
		-b $bigwig_dir/${f_id1}*.dedup.bigwig $bigwig_dir/${f_id2}*.dedup.bigwig \
		--labels ${l_list[@]} \
		-out $deep_dir/summary/${f_id}_scores_per_bin.npz \
		--outRawCounts $deep_dir/counts/${f_id}_scores_per_bin.tab" > $sh
	sh $sh
fi
	
if [[ -f $deep_dir/gz/${f_id}__dedup__${type}_single.gz && ! -f $deep_dir/summary/${f_id}_scores_per_bin_single.npz ]]; then
	echo "--summary $f_id single"
	sh="$deep_dir/sh/summary_${f_id}_single.sh"
			
	echo "multiBigwigSummary bins \
		-b $bigwig_dir/${f_id1}_1.dedup.bigwig $bigwig_dir/${f_id2}_1.dedup.bigwig \
		--labels ${l_list[0]} ${l_list[3]} \
		-out $deep_dir/summary/${f_id}_scores_per_bin_single.npz \
		--outRawCounts $deep_dir/counts/${f_id}_scores_per_bin_single.tab" > $sh
	sh $sh
fi
##################################################################
# plot heatmap, profile
##################################################################
#https://deeptools.readthedocs.io/en/develop/content/tools/plotHeatmap.html
if [[ -f $deep_dir/summary/${f_id}_scores_per_bin.npz && ! -f $deep_dir/plots/${f_id}_profile.png ]]; then
	echo "--plot1 $f_id full"
	sh="$deep_dir/sh/plot1_${f_id}.sh"
	echo "
    plotProfile \
	-m $deep_dir/gz/${f_id}__dedup__${type}.gz \
	-out $deep_dir/plots/${f_id}_profile.png \
	--perGroup \
	--colors green purple red blue orange yellow \
	--samplesLabel "${l_list[@]}" \
	--refPointLabel "$center_type" \
	-T \"Read density\" \
	-z \"\"" > $sh
	sh $sh
	
	echo "--plot2 $f_id full"
	sh="$deep_dir/sh/plot2_${f_id}.sh"
	echo "
    plotHeatmap \
	-m $deep_dir/gz/${f_id}__dedup__${type}.gz \
	-out $deep_dir/plots/${f_id}_heatmap.png \
	--colorList white,red \
	--samplesLabel "${l_list[@]}""  > $sh
	sh $sh
fi 

if [[ -f $deep_dir/summary/${f_id}_scores_per_bin_single.npz ]]; then
	if [[ ! -f $deep_dir/plots/${f_id}_profile_single.png ]]; then
        echo "--plot1 $f_id single"
        sh="$deep_dir/sh/plot1_${f_id}_single.sh"
        echo "
        plotProfile \
        -m $deep_dir/gz/${f_id}__dedup__${type}_single.gz \
        -out $deep_dir/plots/${f_id}_profile_single.png \
        --perGroup \
        --colors green red \
        --samplesLabel "${l_list[0]} ${l_list[3]}" \
        --refPointLabel "$center_type" \
        -T \"Read density\" \
        -z \"\"" > $sh
        sh $sh
    fi
	
	if [[ ! -f $deep_dir/plots/${f_id}_heatmap_single.png ]]; then
        echo "--plot2 $f_id single"
        sh="$deep_dir/sh/plot2_${f_id}_single.sh"
        #  --heatmapHeight 5 --heatmapWidth 3.5 \
        echo "
        plotHeatmap \
        -m $deep_dir/gz/${f_id}__dedup__${type}_single.gz \
        --yMax 40 40 \
        --samplesLabel "${l_list[0]} ${l_list[3]}" \
        --colorList white,red \
        -out $deep_dir/plots/${f_id}_heatmap_single.png"  > $sh
        sh $sh
    fi
fi