#!/usr/bin/env bash
module load R/4.0

# dir 
project_dir="/data/CCBR/projects/ccbr1155"
R_script="${project_dir}/sam/diff_markdown_wrapper.R"

# project files
#project_id="CS030586" #original
#project_id="CS029901"
#project_id="CS030031"
#project_id="CS031188" #CRISPR vishal ran 1 sample
#project_id="CS031188" #CRISPR
#project_id="CS030666" #human
project_id="CS031014" #human

#type of analysis
# deseq or tracks
analysis_type="tracks"
pass_id="pass1"

if [ $project_id == "CS030586" ]; then
    dir_loc="CS030586_CARAP"
    deg_list=("siSmyd3_2m_Smyd3_0.25HCHO_500K_vs_siNC_2m_Smyd3_0.25HCHO_500K" "siSmyd3_5m_Smyd3_0.25HCHO_500K_vs_siNC_5m_Smyd3_0.25HCHO_500K")
elif [ $project_id == "CS029901" ]; then
    dir_loc="CS029901_CS030031_complete"
    deg_list=("siSmyd3_H3K4me3_vs_siNC_H3K4me3"  "siSmyd3_H3K9me3_vs_siNC_H3K9me3" "siSmyd3_H3K27me3_vs_siNC_H3K27me3")
elif [ $project_id == "CS030031" ]; then
    dir_loc="CS029901_CS030031_complete"
    deg_list=("siSmyd3_H4K20me3_vs_siNC_H4K20me3" "siSmyd3_H3K4me1_vs_siNC_H3K4me1" "siSmyd3_H3K27Ac_vs_siNC_H3K27Ac")
elif [ $project_id == "CS031188_v1" ]; then
    dir_loc="CS031188_Symd3_carlisle"
    deg_list=("MOC1_CRISPR_Smyd3_KO_10_Smyd3_350K_0_35ng_0_25_HCHO_2m_vs_MOC1_CRISPR_NC_2_Smyd3_350K_0_35ng_0_25_HCHO_2m")
elif [ $project_id == "CS031188" ]; then
    dir_loc="CS031188_complete"
    deg_list=("CRISPR_NC_2_H3K4me3_vs_CRISPR_KO_10_H3K4me3" "CRISPR_NC_2_H3K27me3_vs_CRISPR_KO_10_H3K27me3" "CRISPR_NC_2_H4K20me3_vs_CRISPR_KO_10_H4K20me3" "CRISPR_NC_2_Smyd3_vs_CRISPR_KO_10_Smyd3")
elif [ $project_id == "CS030666" ]; then
    dir_loc="CS030666/analysis"
    deg_list=("siNC_H3K27me3_vs_siSmyd3_H3K27me3" "siNC_H3K4me3_vs_siSmyd3_H3K4me3" "siNC_H4K20me3_vs_siSmyd3_H4K20me3")
elif [ $project_id == "CS031014" ]; then
    dir_loc="CS031014/analysis_v3"
    deg_list=("53_H3K4me3_vs_HN6_H3K4me3" "53_H4K20m3_vs_HN6_H4K20me3")
fi

#set args
IFS=$'\n' read -d '' -r -a sample_list < "sample_list_${project_id}.txt"
analysis_dir="${project_dir}/${dir_loc}/results/peaks/contrasts"
output_dir="${analysis_dir}/sam_${pass_id}"
track_dir="/data/CCBR/datashare/ccbr1155/${dir_loc}/"
project_exclusion_file="${project_dir}/sam/sample_exclusion_${project_id}.csv"
col_value="peakID" # comma separated list of indexing columns eg. gene_id \gene_name
fdr_value="0.05"
log2fc_value="0.59"
bed_graph_file="${analysis_dir}/bed_bedgraph_paths.tsv"

# peak_list=("narrowPeak" "norm.relaxed.bed" "norm.stringent.bed")
# method_list=("AUCbased" "fragmentsbased")
peak_list=("norm.relaxed.bed")
method_list=("fragmentsbased")
dedup_list=("dedup")

# make output dir
if [[ ! -d $output_dir/reports ]]; then mkdir -p $output_dir/reports; fi
if [[ ! -d $output_dir/results ]]; then mkdir -p $output_dir/results; fi
if [[ ! -d $output_dir/exclusions ]]; then mkdir -p $output_dir/exclusions; fi
if [[ ! -d $track_dir/bigwig ]]; then mkdir -p $track_dir/bigwig; fi
if [[ ! -d $track_dir/bigbed ]]; then mkdir -p $track_dir/bigbed; fi
if [[ -f $output_dir/results/track_info.txt ]]; then rm $output_dir/results/track_info.txt; fi

# for each sample \ run R markdown 
run_deseq (){
    peak_type=$1
    method_type=$2
    dedup_type=$3
    sample_id=$4
    
    # sample name
    # eg siSmyd3_2m_Smyd3_0.25HCHO_500K_vs_siNC_2m_Smyd3_0.25HCHO_500K__no_dedup__norm.relaxed
    complete_sample_id="${sample_id}__${dedup_type}__${peak_type}"

    # add extensions for bed and counts file
    sample_info_file="${analysis_dir}/${complete_sample_id}/${complete_sample_id}_sampleinfo.txt"
    sample_exclusion_file="${output_dir}/exclusions/${complete_sample_id}_exclusion.tsv"
    counts_file="${analysis_dir}/${complete_sample_id}/${complete_sample_id}_${method_type}countsmatrix.txt"
    elbow_file="${analysis_dir}/${complete_sample_id}/${complete_sample_id}_${method_type}_diffanalysis_elbowlimits.tmp.yaml"
    results_file="${output_dir}/results/${complete_sample_id}_${method_type}_diffresults.txt"
    report_file="${output_dir}/reports/${complete_sample_id}_${method_type}_report.html"

    # create sample specific exclusions for pass 2
    if [[ $pass_id == "pass_2" ]]; then 
        cat $project_exclusion_file | grep "${complete_sample_id}_${method_type}" | cut -f2 -d"," > $sample_exclusion_file
    fi

    # fix naming for counts matrix
    rawcountsmatrix=`echo $counts_file | sed "s/fragmentsbased/fragments/g" | sed "s/AUCbased//g"`

    # set inputs for RMarkdown
    results=$results_file
    report=$report_file
    rawcountsmatrix=$rawcountsmatrix
    sampleinfo=$sample_info_file
    exclusionlist=$sample_exclusion_file
    bbpaths=$bed_graph_file
    elbowlimits=$elbow_file
    condition1=`echo $sample_id | awk -F"_vs_" '{ print $1 }'` # contrasts is condition1 vs condition2 ... pay attention to the order of conditions
    condition2=`echo $sample_id | awk -F"_vs_" '{ print $2 }'`
    dupstatus=$dedup_type
    indexcols=$col_value
    fdr_cutoff=$fdr_value
    log2fc_cutoff=$log2fc_value
    spiked="Y"
    rawcountsprescaled="N"
    scalesfbymean="Y"
    htsfilter="Y" # Use HTSFilter (CPM filter does not work well for this type of data)


    if [[ ! -f $report_file ]]; then
        echo "*********************************************"
        echo "Sample PROCESSING: ${complete_sample_id}_${method_id}"
        echo
        Rscript $R_script \
        --countsmatrix $rawcountsmatrix \
        --sampleinfo $sampleinfo \
        --exclusionlist $exclusionlist \
        --dupstatus $dupstatus \
        --condition1 $condition1 \
        --condition2 $condition2 \
        --indexcols $indexcols \
        --bbpaths $bbpaths \
        --elbowlimits $elbowlimits \
        --fdr_cutoff $fdr_cutoff \
        --log2fc_cutoff $log2fc_cutoff \
        --spiked $spiked \
        --rawcountsprescaled $rawcountsprescaled \
        --scalesfbymean $scalesfbymean \
        --htsfilter $htsfilter \
        --results $results \
        --report $report
    else
        echo "*********************************************"
        echo "Sample COMPLETE: ${complete_sample_id}_${method_id}"
    fi
}

run_comparison_tracks (){
    peak_type=$1
    method_type=$2
    dedup_type=$3
    sample_id=$4
    
    # sample name
    # eg siSmyd3_2m_Smyd3_0.25HCHO_500K_vs_siNC_2m_Smyd3_0.25HCHO_500K__no_dedup__norm.relaxed
    complete_sample_id="${sample_id}__${dedup_type}__${peak_type}"
    
    #create hard links
    source_loc="${analysis_dir}/${complete_sample_id}/${complete_sample_id}_fragmentsbased_diffresults.bigbed "
    link_loc="${track_dir}/bigbed/${complete_sample_id}_fragmentsbased_diffresults.bigbed"

    if [[ ! -f $link_loc ]]; then
        ln $source_loc $link_loc
    fi

    # echo track info
    echo "track name=${sample_id} bigDataUrl=https://hpc.nih.gov/~CCBR/ccbr1155/${dir_loc}/bigbed/${complete_sample_id}_fragmentsbased_diffresults.bigbed type=bigBed itemRgb=On" >> $output_dir/results/track_info.txt
}

run_sample_tracks (){
    dedup_type=$1
    sample_id=$2
    
    # sample name
    # eg siNC_H3K27Ac_1.dedup.bigwig
    complete_sample_id="${sample_id}.${dedup_type}"
    
    #create hard links
    source_loc="/data/CCBR/projects/ccbr1155/${dir_loc}/results/bigwig/${complete_sample_id}.bigwig "
    link_loc="${track_dir}/bigwig/${complete_sample_id}.bigwig"

    if [[ ! -f $link_loc ]]; then
        ln $source_loc $link_loc
    fi

    # echo track info
    echo "track type=bigWig bigDataUrl=https://hpc.nih.gov/~CCBR/ccbr1155/${dir_loc}/bigwig/${complete_sample_id}.bigwig name=${sample_id} description=${sample_id} visibility=full autoScale=off maxHeightPixels=128:30:1 viewLimits=1:120 color=65,105,225" >> $output_dir/results/track_info.txt
}

# iterate through samples / peaks / methods / dedup
for sample_id in ${deg_list[@]}; do
    for peak_id in ${peak_list[@]}; do
        for method_id in ${method_list[@]}; do
            for dedup_id in ${dedup_list[@]}; do
                if [[ $analysis_type == "deseq" ]]; then
                    run_deseq $peak_id $method_id $dedup_id $sample_id
                elif [[ $analysis_type == "tracks" ]]; then
                    run_comparison_tracks $peak_id $method_id $dedup_id $sample_id
                else
                    echo "FIX analysis_type"
                fi
            done
        done
    done
done

#iterate through samples / dedup
for sample_id in ${sample_list[@]}; do
    for dedup_id in ${dedup_list[@]}; do
        run_sample_tracks $dedup_id $sample_id
    done
done

cat $output_dir/results/track_info.txt


# store R params
#   params$species="hg38"
#   params$rawcountsmatrix="~/../../Volumes/ccbr1155/CS031014/analysis_v3/results/peaks/contrasts/53_H3K4me3_vs_HN6_H3K4me3__dedup__norm.relaxed.bed/53_H3K4me3_vs_HN6_H3K4me3__dedup__norm.relaxed.bed_fragmentscountsmatrix.txt"
#   params$coldata="~/../../Volumes/ccbr1155/CS031014/analysis_v3/results/peaks/contrasts/53_H3K4me3_vs_HN6_H3K4me3__dedup__norm.relaxed.bed/53_H3K4me3_vs_HN6_H3K4me3__dedup__norm.relaxed.bed_sampleinfo.txt"
#   params$dupstatus="dedup" # dedup or no_dedup
#   params$condition1="53_H3K4me3"
#   params$condition2="HN6_H3K4me3" # contrasts is condition1 vs condition2 ... pay attention to the order of conditions
#   params$results="~/../../Volumes/ccbr1155/CS031014/analysis_v3/results/peaks/contrasts/53_H3K4me3_vs_HN6_H3K4me3__dedup__norm.relaxed.bed/53_H3K4me3_vs_HN6_H3K4me3__dedup__norm.relaxed.bed_fragmentsbased_diffresults.txt.tmp"
#   params$bbpaths="~/../../Volumes/ccbr1155/CS031014/analysis_v3/results/peaks/contrasts/bed_bedgraph_paths.tsv"
#   params$elbowlimits="~/../../Volumes/ccbr1155/CS031014/analysis_v3/results/peaks/contrasts/53_H3K4me3_vs_HN6_H3K4me3__dedup__norm.relaxed.bed/53_H3K4me3_vs_HN6_H3K4me3__dedup__norm.relaxed.bed_fragmentsbased_diff_elbowlimits.tmp.yaml"
