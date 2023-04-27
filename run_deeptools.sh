
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
sample_type=$5

# center, TSS
center_type="center"

# load modules
module load deeptools bedtools

# Example Usage
# sh run_deeptools.sh CS029758 NA genomic_regions 15000 single
# sh run_deeptools.sh CS029758 NA immune_genes 15000 single
##################################################################
# Set samples
##################################################################
if [[ $cs_id == "CS029758" ]]; then
	carlisle_id="carlisle_230111"
    deep_id="deeptools_230324"
	deeptools_genomic="/data/CUTRUN/analysis/CS029758/r_analysis_230321_newannotations/contrast_deeptools_53_H4K20me3_IFNb_vs_HN6_H4K20me3_IFNb_rna_PIE53.bed"
	deeptools_immunogenic="/data/CUTRUN/analysis/CS029758/r_analysis_230321_newannotations/contrast_deeptools_53_H4K20me3_IFNb_vs_HN6_H4K20me3_IFNb_rna_immune_PIE53.bed"
	
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
ref_dir="/data/CUTRUN/analysis/$cs_id/$deep_id/${f_id}/refs"
deep_dir="/data/CUTRUN/analysis/$cs_id/$deep_id/${f_id}/$min_max/$gene_shorthand"
project_dir="/data/CUTRUN/analysis/$cs_id/$carlisle_id/results"
bigwig_dir="$project_dir/bigwig"

# make dirs
dir_list=($deep_dir)
for newdir in ${dir_list[@]}; do
	if [[ ! -d $newdir ]]; then mkdir -p $newdir; fi
done

dir_list=(sh logs plots gz bed summary counts)
for newdir in ${dir_list[@]}; do
	if [[ ! -d $deep_dir/$newdir ]]; then mkdir -p $deep_dir/$newdir; fi
done
##################################################################
# prep bed
##################################################################
sample_genelist_specific_bed="$deep_dir/bed/${f_id}_${gene_shorthand}_samplespecific.bed"
awk '{print $2"\t"$3"\t"$4}' $deeptools_input | sed "s/\"//g" | grep -v "start" | sort > $sample_genelist_specific_bed
head $sample_genelist_specific_bed
cat $sample_genelist_specific_bed | wc -l

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
if [[ $sample_type == "multiple" ]]; then
	echo "--Running all samples"
	if [[ -f $bigwig_dir/${f_id1}_1.dedup.bigwig && ! -f $deep_dir/gz/${f_id}__dedup__${type}.gz ]]; then
		echo "--matrix $f_id full"
		sh="$deep_dir/sh/matrix_${f_id}.sh"

		echo "#!/bin/sh
		module load deeptools
		computeMatrix reference-point \
			--referencePoint $center_type -b $min_max -a $min_max -p 6 \
			-R $sample_genelist_specific_bed \
			-S $bigwig_dir/${f_id1}*.dedup.bigwig $bigwig_dir/${f_id2}*.dedup.bigwig  \
			-o $deep_dir/gz/${f_id}__dedup__${type}.gz \
			--missingDataAsZero --skipZeros --sortRegions descend \
			--outFileSortedRegions $deep_dir/bed/${f_id}_regions_$center_type__dedup__${type}.bed"  > $sh
		# submit_batch $sh
		sh $sh
	fi
else
	echo "--Running single samples"
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
		submit_batch $sh
		# sh $sh
	fi
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
        --yMax 20 20 \
        --samplesLabel "${l_list[0]} ${l_list[3]}" \
        --colorList white,red \
        -out $deep_dir/plots/${f_id}_heatmap_single.png"  > $sh
        sh $sh
    fi
fi