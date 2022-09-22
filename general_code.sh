cs_id="CS031014"
#cs_id="CS029689"
carlisle_dir="carlisle_220920"
project_dir="/data/CCBR/projects/ccbr1155/$cs_id/$carlisle_dir/results/"
code_dir="/data/CCBR/projects/ccbr1155/sam/github/"

##################################################################
flag_scaling="Y"
flag_dedup="Y"

##################################################################
#Scaling factor code
if [[ $flag_scaling == "Y" ]]; then
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
fi

# deduplication summary code
if [[ $flag_dedup == "Y" ]]; then
	cd $project_dir
	#create dedup summary file
	echo "sampleid,dedup_nreads_genome,dedup_nreads_spikein,nodedup_nreads_genome,nodedup_nreads_spikein,nreads,raw_nreads_genome,raw_nreads_spikein" > $project_dir/dedup_summary_${cs_id}.txt
	for f in alignment_stats/*yaml; do 
		name=`echo $f | cut -f2 -d"/" | cut -f1 -d"."`
		vals=`cat $f | awk -F": " '{print $2}' | awk 'BEGIN { ORS = "," } { print }'`
		echo "$name,$vals">> $project_dir/dedup_summary_${cs_id}.txt
	done
	cd $code_dir
fi