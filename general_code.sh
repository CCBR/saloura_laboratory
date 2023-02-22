##################################################################
flag=$1

cs_id=$2
sampleid=$3

#overlap, sig
sig_name=$4

# sh general_code.sh create_filt_bed CS031014 two sig; sh general_code.sh deeptools CS031014 two sig
# sh general_code.sh create_filt_bed CS031308 one sig; sh general_code.sh deeptools CS031308 one sig
# sh general_code.sh create_filt_bed CS031308 two sig; sh general_code.sh deeptools CS031308 two sig
##################################################################
if [[ $cs_id == "CS029758" ]]; then
	carlisle_id="carlisle_221013"
	r_id=""
elif [[ $cs_id == "CS031014" ]]; then
	r_id="CS031014_CS028891_complete/r_analysis_221223"
	carlisle_id="carlisle_221206"

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
elif [[ $cs_id == "CS031308" ]]; then
	r_id="CS031308_CS028891_complete/r_analysis_221226"
	carlisle_id="carlisle_221226"

	if [[ $sampleid == "one" ]]; then
		f_list=("53_INFB_SMYD3_1" "53_INFB_SMYD3_2" "53_INFB_SMYD3_3" "HN6_INFB_SMYD3_1" "HN6_INFB_SMYD3_2" "HN6_INFB_SMYD3_3")
		l_list=("5-3_Rep1" "5-3_Rep2" "5-3_Rep3" "HN6_Rep1" "HN6_Rep2" "HN6_Rep3")
		type="narrowGo"
		f_id="53_INFB_SMYD3_vs_HN6_INFB_SMYD3"
		f_id1="53_INFB_SMYD3"
		f_id2="HN6_INFB_SMYD3"
	else
		f_list=("53_INFB_UHRF1_1" "53_INFB_UHRF1_1" "53_INFB_UHRF1_1" "HN6_INFB_UHRF1_1" "HN6_INFB_UHRF1_1" "HN6_INFB_UHRF1_1")
		l_list=("5-3_Rep1" "5-3_Rep2" "5-3_Rep3" "HN6_Rep1" "HN6_Rep2" "HN6_Rep3")
		type="narrowGo"
		f_id="53_INFB_UHRF1_vs_HN6_INFB_UHRF1"
		f_id1="53_INFB_UHRF1"
		f_id2="HN6_INFB_UHRF1"
	fi
elif [[ $cs_id == "CS033351" ]]; then
	r_id="CS033351/r_analysis_230101"
	carlisle_id="carlisle_221231"
fi

echo "Running $cs_id"

##################################################################
code_dir="/data/CUTRUN/analysis/github"

# set dirs
project_dir="/data/CUTRUN/analysis/$cs_id/$carlisle_id/results/"
deep_dir="/data/CUTRUN/analysis/$cs_id/deeptools_221226/${f_id}/$sig_name"
r_dir="/data/CUTRUN/analysis/$r_id"
bigwig_dir="$project_dir/bigwig"

# can create a sub bed with only genes
# source gtf - complete annotation
# sub_bed is genes only
# target is in corrected format
source_gtf="/data/CCBR_Pipeliner/db/PipeDB/Indices/GTFs/hg38/gencode.v30.annotation.gtf"
sub_bed="/data/CCBR_Pipeliner/Pipelines/CARLISLE_dev/resources/sub_gencode.v30.bed"
target_bed="/data/CCBR_Pipeliner/Pipelines/CARLISLE_dev/resources/gencode.v30.bed"

# can create target bed based on filtered list
source_gtf="/data/CCBR_Pipeliner/db/PipeDB/Indices/GTFs/hg38/gencode.v30.annotation.gtf"
target_bed="$deep_dir/refs/${f_id}_${sig_name}.bed"

# center, TSS
center_type="center"

#######################################################################
# new project through rename file
if [[ $flag == "new_project" ]]; then
	raw_dir="/data/CCBR/rawdata/ccbr1155/$cs_id"; \
	proj_dir="/data/CUTRUN/analysis/$cs_id/$carlisle_id/fastq"
	
	#creates soft links
	if [[ ! -d $proj_dir ]]; then mkdir -p $proj_dir; fi; \
	for f in $raw_dir/*/*/*fastq.gz; do \
		ln -s $f $proj_dir/$(basename "$f");\
	done;\
	#ls $proj_dir/*

	#remove _001 and _S[#] from files
	for f in $proj_dir/*; do \
		new_name=`echo $f | sed -s "s/_S[0-9][0-9]//g"`
		new_name=`echo $new_name | sed -s "s/_S[0-9]//g"`
		new_name=`echo $new_name | sed -s "s/_001//g"`
		new_name=`echo $new_name | sed -s "s/_R/.R/g"`
		mv $f $new_name
	done;\
	#ls $proj_dir/*

	#create rename file
	ls $proj_dir/*fastq.gz | cut -f8 -d "/" > $proj_dir/file_rename.csv
	head $proj_dir/file_rename.csv
fi

# after editing rename file, confinue preparing new project
if [[ $flag == "rename" ]]; then
	proj_dir="/data/CUTRUN/analysis/$cs_id/$carlisle_id/fastq"; \

	#rename from file col1 is old and col2 is new
	echo "" >> $proj_dir/file_rename.csv
	sed -e "s/\r//g" $proj_dir/file_rename.csv > $proj_dir/file_rename_clean.csv
	while IFS=',' read -a files 
	do
		mv "$proj_dir/${files[0]}" "$proj_dir/${files[1]}"
	done < $proj_dir/file_rename_clean.csv
	rm $proj_dir/file_rename.csv
	head $proj_dir/file_rename_clean.csv

fi

if [[ $flag == "groups" ]]; then
	proj_dir="/data/CUTRUN/analysis/$cs_id/$carlisle_id/fastq"; \

	#create groups file
	#col1 is base of fastq.gz; col2 is groups for coloring; col3 for sample labels
	#filenames cannot have "-"
	ls $proj_dir/* | cut -f8 -d "/" | sed 's/.R[1-2].fastq.gz//g' | uniq > $proj_dir/col1.txt
	cat $proj_dir/col1.txt | sed 's/_[1-5]//g' > $proj_dir/col2.txt
	cat $proj_dir/col1.txt | sed 's/_IfnB/_/g' | sed 's/CRISPR_//g'| sed 's/parental_//g' | sed 's/HN6_//g' > $proj_dir/col3.txt
	paste $proj_dir/col1.txt $proj_dir/col2.txt $proj_dir/col3.txt > $proj_dir/col_groups.tab
	sed -e "s/\r//g" $proj_dir/col_groups.tab > $proj_dir/groups.tab
	cat $proj_dir/groups.tab
	rm $proj_dir/col*

	# #example
	# #fastqname  group   shorthand name
	# Cntrl_S62	Cntrl	Cntrl_1
	# Cntrl_S63	Cntrl	Cntrl_2
	# Cntrl_S64	Cntrl	Cntrl_3
	# Cntrl_S62	TREAT	TREAT_1
	# Cntrl_S63	TREAT	TREAT_2
	# Cntrl_S64	TREAT	TREAT_3

	#create contasts file
	# both cols are groups from above
	# col 1 is the variation col 2 is the control
	#colcol2 is groups for coloring; col3 for sample labels
	awk '{print $2}' $proj_dir/groups.tab | sort | uniq > $proj_dir/contrasts.tab

	# #example
	# TREAT	CNTRL
fi

if [[ $flag == "sample_list" ]]; then
	proj_dir="/data/CUTRUN/analysis/$cs_id"; \
	awk '{ print $1 }' $proj_dir/$carlisle_id/samples.tsv > $proj_dir/col1.txt
	awk '{ print $2 }' $proj_dir/$carlisle_id/samples.tsv > $proj_dir/col2.txt
	paste $proj_dir/col1.txt $proj_dir/col2.txt -d "_" > /data/CUTRUN/sampleinfo/sample_list_CS033351.txt
fi

##################################################################
#Scaling factor code
if [[ $flag == "scaling" ]]; then
	cd $project_dir
	# create spike in summary file
	echo "sampleid,dedup_spikein_reads,dedup_total_reads,spike_in_factor" > $project_dir/spikein_summary_${cs_id}.txt
	for f in bam/*.dedup*idx*; do 
		raw=`cat $f | grep "NC" | awk '{print $3}'`
		total=`awk '{sum+=$3; sum+=$4}END{print sum;}' $f`
		name=`echo $f | cut -f2 -d"/" | cut -f1 -d"."`
		spikein=`cat bedgraph/$name*.yaml|sort|uniq`
		echo "$name,$raw,$total,$spikein"| sed -s "s/spikein_scaling_factor=//g" >> $project_dir/spikein_summary_${cs_id}.txt
	done
	cd $code_dir

	echo "$project_dir/spikein_summary_${cs_id}.txt"
fi
##################################################################
# deduplication summary code
if [[ $flag == "dedup" ]]; then
	cd $project_dir
	#create dedup summary file
	echo "sampleid,dedup_nreads_genome,dedup_nreads_spikein,nodedup_nreads_genome,nodedup_nreads_spikein,nreads,raw_nreads_genome,raw_nreads_spikein" > $project_dir/dedup_summary_${cs_id}.txt
	for f in alignment_stats/*yaml; do 
		name=`echo $f | cut -f2 -d"/" | cut -f1 -d"."`
		vals=`cat $f | awk -F": " '{print $2}' | awk 'BEGIN { ORS = "," } { print }'`
		echo "$name,$vals">> $project_dir/dedup_summary_${cs_id}.txt
	done
	cd $code_dir

	echo "$project_dir/dedup_summary_${cs_id}.txt"
fi
##################################################################
# go peaks code
if [[ $flag == "gopeaks" ]]; then
	echo "Running GoPeaks"

	gopeaks_dir="/home/sevillas2/git/gopeaks"
	proj_dir="/data/CCBR/projects/ccbr1155/CS031014/carlisle_220920/results/bam"
	igg="igG_1.dedup.bam"
	sample_list=("53_H3K4me3_1" "53_H3K4me3_2" "53_H3K4me3_3" "HN6_H3K4me3_1" "HN6_H3K4me3_2" "HN6_H3K4me3_3" "HN6_H4K20me3_1" "HN6_H4K20me3_2" "HN6_H4K20me3_3" "53_H4K20m3_1" "53_H4K20m3_2" "53_H4K20m3_3")
	ext="dedup.bam"
	output_dir="/data/CCBR/projects/ccbr1155/CS031014/gopeaks"

	# create dirs
	if [[ ! -d $output_dir/logs/sh ]]; then mkdir -p $output_dir/logs/sh; fi

	# run each sample
	for sample_name in ${sample_list[@]}; do
		echo "--$sample_name"
		
		#### NARROW
		# create sbatch file
		sbatch_file="$output_dir/logs/sh/${sample_name}_narrow.sh"
		if [[ ! -f $sbatch_file ]]; then	
			echo "----creating sbatch narrow"
			echo "#!/bin/sh"  > $sbatch_file
			echo "cd $gopeaks_dir" >> $sbatch_file
			echo "./gopeaks -b $proj_dir/${sample_name}.${ext} -c $proj_dir/$igg -o $output_dir/${sample_name}_narrow" >> $sbatch_file
    	fi
    
		# submit sbatch files to download SRA's
		output_file=$output_dir/${sample_name}_narrow_peaks.bed
		if [[ ! -f $output_file ]]; then
			echo "----submitting SBATCH narrow"
			#echo "----missing $output_file"
			# cat $sbatch_file
			#sbatch --cpus-per-task=12 --verbose --output=$output_dir/logs/%j.out --mem=200g --gres=lscratch:450 --time 10:00:00 --error=$output_dir/logs/%j.err $sbatch_file
		fi

		#### BROAD
		sbatch_file="$output_dir/logs/sh/${sample_name}_broad.sh"
		if [[ ! -f $sbatch_file ]]; then	
			echo "----creating sbatch broad"
			echo "#!/bin/sh"  > $sbatch_file
			echo "cd $gopeaks_dir" >> $sbatch_file
			echo "./gopeaks -b $proj_dir/${sample_name}.${ext} -c $proj_dir/$igg -o $output_dir/${sample_name}_broad --broad" >> $sbatch_file
    	fi
    
		# submit sbatch files to download SRA's
		output_file=$output_dir/${sample_name}_broad_peaks.bed
		if [[ ! -f $output_file ]]; then
			echo "----submitting SBATCH broad"
			#echo "----missing $output_file"
			# cat $sbatch_file
			sbatch --cpus-per-task=12 --verbose --output=$output_dir/logs/%j.out --mem=200g --gres=lscratch:450 --time 10:00:00 --error=$output_dir/logs/%j.err $sbatch_file
		fi
	done
fi
