#!/usr/bin/env bash

############################################################################################################
# project info
############################################################################################################
# project files
project_id=$1

#CS028891, CS031038

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
    dir_loc="CS030666/carlisle_230118"
    deg_list=("siNC_H3K27me3_vs_siSmyd3_H3K27me3" "siNC_H3K4me3_vs_siSmyd3_H3K4me3" "siNC_H4K20me3_vs_siSmyd3_H4K20me3")
elif [ $project_id == "CS031014" ]; then
    dir_loc="CS031014/carlisle_221206"
    deg_list=("53_H3K4me3_vs_HN6_H3K4me3" "53_H4K20m3_vs_HN6_H4K20me3")
elif [ $project_id == "CS029689" ]; then
    dir_loc="CS029689/carlisle_230111"
    deg_list=("5-3_H3K4me3_IFNb_vs_HN6_H3K4me3_IFNb" "5-3_H3K9me3_IFNb_vs_HN6_H3K9me3_IFNb")
elif [ $project_id == "CS029758" ]; then
    dir_loc="CS029758/carlisle_230111"
    deg_list=("53_H3K27me3_IFNb_vs_HN6_H3K27me3_IFNb" "53_H4K20me3_IFNb_vs_HN6_H4K20me3_IFNb")
elif [ $project_id == "CS031308" ]; then
    dir_loc="CS031308/carlisle_221226"
    deg_list=("53_INFB_SMYD3_vs_HN6_INFB_SMYD3" "53_INFB_UHRF1_vs_HN6_INFB_UHRF1" "53_no_INFB_SMYD3_vs_HN6_no_INFB_SMYD3")
elif [ $project_id == "CS033351" ]; then
    dir_loc="CS033351/carlisle_221231"
    deg_list=("HA_Mock_IFNb_antiHA_1min_vs_HA_UHRF1_IFNb_antiHA_1min" "HA_Mock_IFNb_antiHA_2min_vs_HA_UHRF1_IFNb_antiHA_2min" "HA_UHRF1_IFNb_antiHA_1min_vs_HA_Mock_IFNb_antiHA_1min" "HA_UHRF1_IFNb_antiHA_2min_vs_HA_Mock_IFNb_antiHA_2min")
fi

# peak_list=("narrowPeak" "norm.relaxed.bed" "norm.stringent.bed")
# method_list=("AUCbased" "fragmentsbased")
#peak_list=("norm.relaxed.bed" "narrowPeak" "broadGo_peaks.bed" "narrowGo_peaks.bed")
peak_list=("broadGo_peaks.bed" "narrowGo_peaks.bed")
method_list=("fragmentsbased")
dedup_list=("dedup")

############################################################################################################
#set args
############################################################################################################
#project_dir="/data/CCBR/projects/ccbr1155"
project_dir="/data/CUTRUN/analysis/$project_id"
track_info="$project_dir/track_info.txt"

analysis_dir="/data/CUTRUN/analysis/$dir_loc/results"
bed_graph_file="${analysis_dir}/peaks/.05/contrasts/bed_bedgraph_paths.tsv"
bigwig_dir="$analysis_dir/bigwig"

track_dir="/data/CUTRUN/datashare/CCBR/${project_id}"

# make output dirs
if [[ ! -d $track_dir/bigwig ]]; then echo "making bigwig dir"; mkdir -p $track_dir/bigwig; fi
if [[ ! -d $track_dir/bigbed ]]; then echo "making bigbed dir"; mkdir -p $track_dir/bigbed; fi

# read sample file
sample_list="samples_list.txt"
IFS=$'\n' read -d '' -r -a sample_list < "$project_dir/$sample_list"
#example
#HA_Mock_IFNb_IgG_neg_control_1min_1
#HA_Mock_IFNb_antiHA_1min_1
#HA_Mock_IFNb_antiHA_1min_2
#HA_Mock_IFNb_antiHA_1min_3
#awk '{print $1"_"$2}' CS029758/samples.tsv | sort | uniq> CS029758/samples_list.txt

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
    source_loc="$bigwig_dir/${complete_sample_id}.bigwig "
    link_loc="${track_dir}/bigwig/${complete_sample_id}.bigwig"
    if [[ ! -f $link_loc ]]; then ln $source_loc $link_loc; fi

    # echo track info
    echo "track type=bigWig bigDataUrl=https://hpc.nih.gov/~CUTRUN/CCBR/${project_id}/bigwig/${complete_sample_id}.bigwig name=${sample_id} description=${sample_id} visibility=full autoScale=off maxHeightPixels=128:30:1 viewLimits=1:120 color=65,105,225" >> $track_info
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
    source_loc="${analysis_dir}/peaks/0.05/contrasts/${complete_sample_id}/${complete_sample_id}_fragmentsbased_diffresults.bigbed "
    link_loc="${track_dir}/bigbed/${complete_sample_id}_fragmentsbased_diffresults.bigbed"
    if [[ ! -f $link_loc ]]; then ln $source_loc $link_loc; fi

    # echo track info
    echo "track name=${sample_id}_${peak_type} bigDataUrl=https://hpc.nih.gov/~CUTRUN/CCBR/${project_id}/bigbed/${complete_sample_id}_fragmentsbased_diffresults.bigbed type=bigBed itemRgb=On" >> $track_info
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