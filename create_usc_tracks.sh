#!/usr/bin/env bash

############################################################################################################
# project info
############################################################################################################
# dir 
project_dir="/data/CCBR/projects/ccbr1155"

# project files
project_id="CS031308"

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
    dir_loc="CS031014/carlisle_220920"
    sample_list="samples_list.tsv"
    deg_list=("53_H3K4me3_vs_HN6_H3K4me3" "53_H4K20m3_vs_HN6_H4K20me3")
elif [ $project_id == "CS029689" ]; then
    dir_loc="CS029689/carlisle_220920"
    sample_list="samples_list.tsv"
    deg_list=("5-3_H3K4me3_IFNb_vs_HN6_H3K4me3_IFNb" "5-3_H3K9me3_IFNb_vs_HN6_H3K9me3_IFNb")
elif [ $project_id == "CS029758" ]; then
    dir_loc="CS029758/carlisle_221013"
    sample_list="samples_list.txt"
    deg_list=("53_H3K27me3_IFNb_vs_HN6_H3K27me3_IFNb" "53_H4K20me3_IFNb_vs_HN6_H4K20me3_IFNb")
elif [ $project_id == "CS031308" ]; then
    dir_loc="CS031308/carlisle_221116"
    sample_list="samples_list.txt"
    deg_list=("53_INFB_SMYD3_vs_HN6_INFB_SMYD3" "53_INFB_UHRF1_vs_HN6_INFB_UHRF1" "53_no_INFB_SMYD3_vs_HN6_no_INFB_SMYD3")
fi

# peak_list=("narrowPeak" "norm.relaxed.bed" "norm.stringent.bed")
# method_list=("AUCbased" "fragmentsbased")
peak_list=("norm.relaxed.bed" "narrowPeak" "broadGo_peaks.bed" "narrowGo_peaks.bed")
method_list=("fragmentsbased")
dedup_list=("dedup")

############################################################################################################
#set args
############################################################################################################
analysis_dir="${project_dir}/${dir_loc}/results/peaks/contrasts"
bed_graph_file="${analysis_dir}/bed_bedgraph_paths.tsv"
track_dir="/data/CCBR/datashare/ccbr1155/${dir_loc}/"
track_info="$project_dir/$dir_loc/results/track_info.txt"

# make output dir
if [[ ! -d $track_dir/bigwig ]]; then mkdir -p $track_dir/bigwig; fi
if [[ ! -d $track_dir/bigbed ]]; then mkdir -p $track_dir/bigbed; fi
if [[ -f $project_dir/results/track_info.txt ]]; then rm $project_dir/results/track_info.txt; fi

# read sample file
IFS=$'\n' read -d '' -r -a sample_list < "$project_dir/$project_id/samples_list.txt"

# remove previous track info
if [[ -f $track_info ]]; then rm $track_info; touch $track_info; fi

############################################################################################################
# run sample 
############################################################################################################
run_sample_tracks (){
    sample_id=$1
    dedup_id=$2
    
    # sample name
    # eg siNC_H3K27Ac_1.dedup.bigwig
    complete_sample_id="${sample_id}.${dedup_id}"

    #create hard links
    source_loc="/data/CCBR/projects/ccbr1155/${dir_loc}/results/bigwig/${complete_sample_id}.bigwig "
    link_loc="${track_dir}/bigwig/${complete_sample_id}.bigwig"

    if [[ ! -f $link_loc ]]; then ln $source_loc $link_loc; fi

    # echo track info
    echo "track type=bigWig bigDataUrl=https://hpc.nih.gov/~CCBR/ccbr1155/${dir_loc}/bigwig/${complete_sample_id}.bigwig name=${sample_id} description=${sample_id} visibility=full autoScale=off maxHeightPixels=128:30:1 viewLimits=1:120 color=65,105,225" >> $track_info
}

# iterate through samples
# at the sample level only DEDUP matters
for sample_id in ${sample_list[@]}; do
    for dedup_id in ${dedup_list[@]}; do
        run_sample_tracks $sample_id $dedup_id
    done
done

############################################################################################################
# run contrasts
############################################################################################################
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

    if [[ ! -f $link_loc ]]; then ln $source_loc $link_loc; fi

    # echo track info
    echo "track name=${sample_id} bigDataUrl=https://hpc.nih.gov/~CCBR/ccbr1155/${dir_loc}/bigbed/${complete_sample_id}_fragmentsbased_diffresults.bigbed type=bigBed itemRgb=On" >> $track_info
}

# iterate through samples / peaks / methods / dedup
for sample_id in ${deg_list[@]}; do
    for peak_id in ${peak_list[@]}; do
        for method_id in ${method_list[@]}; do
            for dedup_id in ${dedup_list[@]}; do
                run_comparison_tracks $peak_id $method_id $dedup_id $sample_id
            done
        done
    done
done

chmod -R 777 $track_dir
echo $track_info