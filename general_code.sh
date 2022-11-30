cs_id="CS029758"
carlisle_dir="carlisle_221013"
project_dir="/data/CCBR/projects/ccbr1155/$cs_id/$carlisle_dir/results/"
code_dir="/data/CCBR/projects/ccbr1155/github"

##################################################################
flag=$1
sampleid=$2

##################################################################
# #Scaling factor code
# if [[ $flag_scaling == "Y" ]]; then
# 	cd $project_dir
# 	# create spike in summary file
# 	echo "sampleid,dedup_spikein_reads,dedup_total_reads,spike_in_factor" > $project_dir/spikein_summary_${cs_id}.txt
# 	for f in bam/*.dedup*idx*; do 
# 		raw=`cat $f | grep "NC" | awk '{print $3}'`
# 		total=`awk '{sum+=$3; sum+=$4}END{print sum;}' $f`
# 		name=`echo $f | cut -f2 -d"/" | cut -f1 -d"."`
# 		spikein=`cat bedgraph/$name*.yaml|sort|uniq`
# 		echo "$name,$raw,$total,$spikein"| sed -s "s/spikein_scaling_factor=//g" >> $project_dir/spikein_summary_${cs_id}.txt
# 	done
# 	cd $code_dir

# 	echo "$project_dir/spikein_summary_${cs_id}.txt"
# fi

# # deduplication summary code
# if [[ $flag_deeptools == "Y" ]]; then
# 	cd $project_dir
# 	#create dedup summary file
# 	echo "sampleid,dedup_nreads_genome,dedup_nreads_spikein,nodedup_nreads_genome,nodedup_nreads_spikein,nreads,raw_nreads_genome,raw_nreads_spikein" > $project_dir/dedup_summary_${cs_id}.txt
# 	for f in alignment_stats/*yaml; do 
# 		name=`echo $f | cut -f2 -d"/" | cut -f1 -d"."`
# 		vals=`cat $f | awk -F": " '{print $2}' | awk 'BEGIN { ORS = "," } { print }'`
# 		echo "$name,$vals">> $project_dir/dedup_summary_${cs_id}.txt
# 	done
# 	cd $code_dir

# 	echo "$project_dir/dedup_summary_${cs_id}.txt"
# fi

# if [[ $flag_gopeaks == "Y" ]]; then
# 	echo "Running GoPeaks"

# 	gopeaks_dir="/home/sevillas2/git/gopeaks"
# 	proj_dir="/data/CCBR/projects/ccbr1155/CS031014/carlisle_220920/results/bam"
# 	igg="igG_1.dedup.bam"
# 	sample_list=("53_H3K4me3_1" "53_H3K4me3_2" "53_H3K4me3_3" "HN6_H3K4me3_1" "HN6_H3K4me3_2" "HN6_H3K4me3_3" "HN6_H4K20me3_1" "HN6_H4K20me3_2" "HN6_H4K20me3_3" "53_H4K20m3_1" "53_H4K20m3_2" "53_H4K20m3_3")
# 	ext="dedup.bam"
# 	output_dir="/data/CCBR/projects/ccbr1155/CS031014/gopeaks"

# 	# create dirs
# 	if [[ ! -d $output_dir/logs/sh ]]; then mkdir -p $output_dir/logs/sh; fi

# 	# run each sample
# 	for sample_name in ${sample_list[@]}; do
# 		echo "--$sample_name"
		
# 		#### NARROW
# 		# create sbatch file
# 		sbatch_file="$output_dir/logs/sh/${sample_name}_narrow.sh"
# 		if [[ ! -f $sbatch_file ]]; then	
# 			echo "----creating sbatch narrow"
# 			echo "#!/bin/sh"  > $sbatch_file
# 			echo "cd $gopeaks_dir" >> $sbatch_file
# 			echo "./gopeaks -b $proj_dir/${sample_name}.${ext} -c $proj_dir/$igg -o $output_dir/${sample_name}_narrow" >> $sbatch_file
#     	fi
    
# 		# submit sbatch files to download SRA's
# 		output_file=$output_dir/${sample_name}_narrow_peaks.bed
# 		if [[ ! -f $output_file ]]; then
# 			echo "----submitting SBATCH narrow"
# 			#echo "----missing $output_file"
# 			# cat $sbatch_file
# 			#sbatch --cpus-per-task=12 --verbose --output=$output_dir/logs/%j.out --mem=200g --gres=lscratch:450 --time 10:00:00 --error=$output_dir/logs/%j.err $sbatch_file
# 		fi

# 		#### BROAD
# 		sbatch_file="$output_dir/logs/sh/${sample_name}_broad.sh"
# 		if [[ ! -f $sbatch_file ]]; then	
# 			echo "----creating sbatch broad"
# 			echo "#!/bin/sh"  > $sbatch_file
# 			echo "cd $gopeaks_dir" >> $sbatch_file
# 			echo "./gopeaks -b $proj_dir/${sample_name}.${ext} -c $proj_dir/$igg -o $output_dir/${sample_name}_broad --broad" >> $sbatch_file
#     	fi
    
# 		# submit sbatch files to download SRA's
# 		output_file=$output_dir/${sample_name}_broad_peaks.bed
# 		if [[ ! -f $output_file ]]; then
# 			echo "----submitting SBATCH broad"
# 			#echo "----missing $output_file"
# 			# cat $sbatch_file
# 			sbatch --cpus-per-task=12 --verbose --output=$output_dir/logs/%j.out --mem=200g --gres=lscratch:450 --time 10:00:00 --error=$output_dir/logs/%j.err $sbatch_file
# 		fi
# 	done
# fi

if [[ $flag == "create_bed" ]]; then
	echo "creating DEEPTOOLS bedfile"

	source_gtf="/data/CCBR_Pipeliner/db/PipeDB/Indices/GTFs/hg38/gencode.v30.annotation.gtf"
	sub_bed="/data/CCBR_Pipeliner/Pipelines/CARLISLE_dev/resources/sub_gencode.v30.bed"
	target_bed="/data/CCBR_Pipeliner/Pipelines/CARLISLE_dev/resources/gencode.v30.bed"
	# subset GTF for only genes
	if [[ ! -f $sub_bed ]]; then
		echo "----subsetting bed"
		awk '$3 == "gene" { print $0 }' $source_gtf > $sub_bed
		echo "Bed here: $sub_bed"
	fi

	# create bed with chrom start stop only
	if [[ ! -f $target_bed ]]; then
		echo "----creating target bed"
		awk '{ print $1"\t"$4"\t"$5 }' $sub_bed > $target_bed
		echo "Bed here: $target_bed"
		head $target_bed
		rm $sub_bed
	fi
fi

#https://hbctraining.github.io/Intro-to-ChIPseq/lessons/10_data_visualization.html
#https://deeptools.readthedocs.io/en/develop/content/tools/computeMatrix.html
if [[ $flag == "deeptools_CS031014" ]]; then
	module load deeptools

	target_bed="/data/CCBR_Pipeliner/Pipelines/CARLISLE_dev/resources/gencode.v30.bed"
	deep_dir="/data/CCBR/projects/ccbr1155/CS031014_CS028891_complete/deeptools"
	r_dir="/data/CCBR/projects/ccbr1155/CS031014_CS028891_complete/r_analysis_220825"
	results_dir="/data/CCBR/projects/ccbr1155/CS031014/carlisle_220920/results"

	if [[ $sampleid == "one" ]]; then
		f_list=("53_H3K4me3_1" "53_H3K4me3_2" "53_H3K4me3_3" "HN6_H3K4me3_1" "HN6_H3K4me3_2" "HN6_H3K4me3_3")
		l_list=("5-3_Rep1" "5-3_Rep2" "5-3_Rep3" "HN6_Rep1" "HN6_Rep2" "HN6_Rep3")
		type="narrowGo"
		f_id="53_H3K4me3_vs_HN6_H3K4me3"
		f_id1="53_H3K4me3"
		f_id2="HN6_H3K4me3"
	else
		f_list=("53_H4K20m3_1" "53_H4K20m3_2" "53_H4K20m3_3" "HN6_H4K20me3_1" "HN6_H4K20me3_2" "HN6_H4K20me3_3")
		l_list=("5-3_Rep1" "5-3_Rep2" "5-3_Rep3" "HN6_Rep1" "HN6_Rep2" "HN6_Rep3")
		type="broadGo"
		f_id="53_H4K20m3_vs_HN6_H4K20me3"
		f_id1="53_H4K20m3"
		f_id2="HN6_H4K20me3"
	fi

	# Normalize the bigwigs
	#https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html
	for s_id in ${f_list[@]}; do
		if [[ ! -f $deep_dir/${s_id}.bw  ]]; then
			echo "--normalize $s_id"
			sh="$deep_dir/sh/norm_${s_id}.sh"
			echo "bamCoverage \
			-b $results_dir/bam/${s_id}.dedup.bam \
			-o $deep_dir/${s_id}.bw \
			--normalizeUsing RPGC \
			--effectiveGenomeSize 2913022398 \
			--centerReads \
			-p 6 2> $deep_dir/logs/${s_id}_bamCoverage.log" > $sh
			sh $sh
		fi
	done

	# Compute Matrix
	if [[ -f $deep_dir/${f_id1}_1.bw && ! -f $deep_dir/${f_id}__dedup__${type}.gz ]]; then
		echo "--matrix $f_id"
		sh="$deep_dir/sh/matrix_${f_id}.sh"
		
		echo "computeMatrix reference-point \
			--referencePoint TSS -b 3000 -a 3000 -p 6 \
			-R $target_bed \
			-S $deep_dir/${f_id1}*.bw $deep_dir/${f_id2}*.bw  \
			-o $deep_dir/${f_id}__dedup__${type}.gz \
			--outFileSortedRegions $deep_dir/${f_id}_regions_TSS__dedup__${type}.bed"  > $sh
		sh $sh
	fi
	if [[ -f $deep_dir/${f_id1}_1.bw && ! -f $deep_dir/${f_id}__dedup__${type}_single.gz ]]; then
		echo "--matrix $f_id"
		sh="$deep_dir/sh/matrix_${f_id}.sh"
		
		echo "computeMatrix reference-point \
			--referencePoint TSS -b 3000 -a 3000 -p 6 \
			-R $target_bed \
			-S $deep_dir/${f_id1}_1.bw $deep_dir/${f_id2}_1.bw  \
			-o $deep_dir/${f_id}__dedup__${type}_single.gz \
			--outFileSortedRegions $deep_dir/${f_id}_regions_TSS__dedup__${type}_single.bed"  > $sh
		sh $sh
	fi
	
	# Summary
	if [[ -f $deep_dir/${f_id}__dedup__${type}.gz && ! -f $deep_dir/${f_id}_scores_per_bin.npz ]]; then
		echo "--summary $f_id"
		sh="$deep_dir/sh/summary_${f_id}.sh"
			
		echo "multiBigwigSummary bins \
			-b $deep_dir/${f_id1}*.bw $deep_dir/${f_id2}*.bw \
			--labels ${l_list[@]} \
			-out $deep_dir/${f_id}_scores_per_bin.npz \
			--outRawCounts $r_dir/${f_id}_scores_per_bin.tab" > $sh
		sh $sh
	fi

	# plots
	if [[ -f $deep_dir/${f_id}_scores_per_bin.npz && ! -f $deep_dir/${f_id}_profile.png ]]; then
		echo "--plot1 $f_id"
		sh="$deep_dir/sh/plot1_${f_id}.sh"
		echo "plotProfile \
		-m $deep_dir/${f_id}__dedup__${type}.gz \
		-out $deep_dir/${f_id}_profile.png \
		--perGroup \
		--colors green purple red blue orange brown \
		--samplesLabel "${l_list[@]}" \
		--refPointLabel "TSS" \
		-T \"Read density\" \
		-z \"\"" > $sh
		sh $sh
	
		echo "--plot2 $f_id"
		sh="$deep_dir/sh/plot2_${f_id}.sh"
		echo "plotHeatmap \
		-m $deep_dir/${f_id}__dedup__${type}.gz \
		-out $deep_dir/${f_id}_heatmap.png \
		--colorMap RdBu \
		--zMin -4 --zMax 4"  > $sh
		sh $sh
	fi 
	if [[ -f $deep_dir/${f_id}_scores_per_bin_single.npz && ! -f $deep_dir/${f_id}_profile_single.png ]]; then
		echo "--plot1 $f_id"
		sh="$deep_dir/sh/plot1_${f_id}.sh"
		echo "plotProfile \
		-m $deep_dir/${f_id}__dedup__${type}_single.gz \
		-out $deep_dir/${f_id}_profile_single.png \
		--perGroup \
		--colors green red \
		--samplesLabel "${l_list[0]} ${l_list[3]}" \
		--refPointLabel "TSS" \
		-T \"Read density\" \
		-z \"\"" > $sh
		sh $sh
	
		echo "--plot2 $f_id"
		sh="$deep_dir/sh/plot2_${f_id}.sh"
		echo "plotHeatmap \
		-m $deep_dir/${f_id}__dedup__${type}_single.gz \
		-out $deep_dir/${f_id}_heatmap_single.png \
		--colorMap RdBu \
		--zMin -4 --zMax 4"  > $sh
		sh $sh
	fi 

	# correlation plot
	if [[ -f $deep_dir/${f_id}_scores_per_bin.npz && ! -f $deep_dir/${f_id}__dedup__${type}_correlation.png ]]; then
		echo "--correlation $f_id"
		sh="$deep_dir/sh/correlation_${f_id}.sh"
		
		echo "plotCorrelation \
		--corData $deep_dir/${f_id}_scores_per_bin.npz \
		--corMethod spearman \
		--whatToPlot heatmap \
		--skipZeros \
		--plotNumbers \
		-o $deep_dir/${f_id}__dedup__${type}_correlation.png" > $sh
		sh $sh
	fi
fi

if [[ $flag == "deeptools_CS031308" ]]; then
	module load deeptools

	target_bed="/data/CCBR_Pipeliner/Pipelines/CARLISLE_dev/resources/gencode.v30.bed"
	deep_dir="/data/CCBR/projects/ccbr1155/CS031308_CS028891_complete/deeptools"
	r_dir="/data/CCBR/projects/ccbr1155/CS031308_CS028891_complete/r_analysis_221129"
	results_dir="/data/CCBR/projects/ccbr1155/CS031308/carlisle_221116/results"

	if [[ $sampleid == "one" ]]; then
		f_list=("53_INFB_SMYD3_1" "53_INFB_SMYD3_2" "53_INFB_SMYD3_3" "HN6_INFB_SMYD3_1" "HN6_INFB_SMYD3_2" "HN6_INFB_SMYD3_3")
		l_list=("5-3_Rep1" "5-3_Rep2" "5-3_Rep3" "HN6_Rep1" "HN6_Rep2" "HN6_Rep3")
		type="narrowGo"
		f_id="53_INFB_SMYD3_vs_HN6_INFB_SMYD3"
		f_id1="53_INFB_SMYD3"
		f_id2="HN6_INFB_SMYD3"
	else
		todo="ok"
	fi

	# Normalize the bigwigs
	#https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html
	for s_id in ${f_list[@]}; do
		if [[ ! -f $deep_dir/${s_id}.bw  ]]; then
			echo "--normalize $s_id"
			sh="$deep_dir/sh/norm_${s_id}.sh"
			echo "bamCoverage \
			-b $results_dir/bam/${s_id}.dedup.bam \
			-o $deep_dir/${s_id}.bw \
			--normalizeUsing RPGC \
			--effectiveGenomeSize 2913022398 \
			--centerReads \
			-p 6 2> $deep_dir/logs/${s_id}_bamCoverage.log" > $sh
			sh $sh
		fi
	done

	# Compute Matrix
	if [[ -f $deep_dir/${f_id1}_1.bw && ! -f $deep_dir/${f_id}__dedup__${type}.gz ]]; then
		echo "--matrix $f_id"
		sh="$deep_dir/sh/matrix_${f_id}.sh"
		
		echo "computeMatrix reference-point \
			--referencePoint TSS -b 3000 -a 3000 -p 6 \
			-R $target_bed \
			-S $deep_dir/${f_id1}*.bw $deep_dir/${f_id2}*.bw  \
			-o $deep_dir/${f_id}__dedup__${type}.gz \
			--outFileSortedRegions $deep_dir/${f_id}_regions_TSS__dedup__${type}.bed"  > $sh
		sh $sh
	fi
	if [[ -f $deep_dir/${f_id1}_1.bw && ! -f $deep_dir/${f_id}__dedup__${type}_single.gz ]]; then
		echo "--matrix $f_id"
		sh="$deep_dir/sh/matrix_${f_id}.sh"
		
		echo "computeMatrix reference-point \
			--referencePoint TSS -b 3000 -a 3000 -p 6 \
			-R $target_bed \
			-S $deep_dir/${f_id1}_1.bw $deep_dir/${f_id2}_1.bw  \
			-o $deep_dir/${f_id}__dedup__${type}_single.gz \
			--outFileSortedRegions $deep_dir/${f_id}_regions_TSS__dedup__${type}_single.bed"  > $sh
		sh $sh
	fi
	
	# Summary
	if [[ -f $deep_dir/${f_id}__dedup__${type}.gz && ! -f $deep_dir/${f_id}_scores_per_bin.npz ]]; then
		echo "--summary $f_id"
		sh="$deep_dir/sh/summary_${f_id}.sh"
			
		echo "multiBigwigSummary bins \
			-b $deep_dir/${f_id1}*.bw $deep_dir/${f_id2}*.bw \
			--labels ${l_list[@]} \
			-out $deep_dir/${f_id}_scores_per_bin.npz \
			--outRawCounts $r_dir/${f_id}_scores_per_bin.tab" > $sh
		sh $sh
	fi

	# plots
	if [[ -f $deep_dir/${f_id}_scores_per_bin.npz && ! -f $deep_dir/${f_id}_profile.png ]]; then
		echo "--plot1 $f_id"
		sh="$deep_dir/sh/plot1_${f_id}.sh"
		echo "plotProfile \
		-m $deep_dir/${f_id}__dedup__${type}.gz \
		-out $deep_dir/${f_id}_profile.png \
		--perGroup \
		--colors green purple red blue orange brown \
		--samplesLabel "${l_list[@]}" \
		--refPointLabel "TSS" \
		-T \"Read density\" \
		-z \"\"" > $sh
		sh $sh
	
		echo "--plot2 $f_id"
		sh="$deep_dir/sh/plot2_${f_id}.sh"
		echo "plotHeatmap \
		-m $deep_dir/${f_id}__dedup__${type}.gz \
		-out $deep_dir/${f_id}_heatmap.png \
		--colorMap RdBu \
		--zMin -4 --zMax 4"  > $sh
		sh $sh
	fi 
	if [[ -f $deep_dir/${f_id}_scores_per_bin_single.npz && ! -f $deep_dir/${f_id}_profile_single.png ]]; then
		echo "--plot1 $f_id"
		sh="$deep_dir/sh/plot1_${f_id}.sh"
		echo "plotProfile \
		-m $deep_dir/${f_id}__dedup__${type}_single.gz \
		-out $deep_dir/${f_id}_profile_single.png \
		--perGroup \
		--colors green red \
		--samplesLabel "${l_list[0]} ${l_list[3]}" \
		--refPointLabel "TSS" \
		-T \"Read density\" \
		-z \"\"" > $sh
		sh $sh
	
		echo "--plot2 $f_id"
		sh="$deep_dir/sh/plot2_${f_id}.sh"
		echo "plotHeatmap \
		-m $deep_dir/${f_id}__dedup__${type}_single.gz \
		-out $deep_dir/${f_id}_heatmap_single.png \
		--colorMap RdBu \
		--zMin -4 --zMax 4"  > $sh
		sh $sh
	fi 

	# correlation plot
	if [[ -f $deep_dir/${f_id}_scores_per_bin.npz && ! -f $deep_dir/${f_id}__dedup__${type}_correlation.png ]]; then
		echo "--correlation $f_id"
		sh="$deep_dir/sh/correlation_${f_id}.sh"
		
		echo "plotCorrelation \
		--corData $deep_dir/${f_id}_scores_per_bin.npz \
		--corMethod spearman \
		--whatToPlot heatmap \
		--skipZeros \
		--plotNumbers \
		-o $deep_dir/${f_id}__dedup__${type}_correlation.png" > $sh
		sh $sh
	fi
fi