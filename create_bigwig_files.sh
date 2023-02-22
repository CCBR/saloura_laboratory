# refs
#https://research.stowers.org/cws/CompGenomics/Tutorial/bam2bigwig.html

module load samtools bedtools ucsc
project_id=$1

pipeliner_dir="/data/CCBR/projects/ccbr1155/$project_id/pipeliner_220919"

track_dir="/data/CCBR/datashare/ccbr1155/${project_id}/bigwig"
track_info="$output_dir/track_info.txt"

if [ $project_id == "CS028891" ]; then
    dir_loc="pipeliner_220713"
    reference="hg38"
    deg_list=("CRISPR_53_with_IFNb-parental_HN6_with_IFNb" "CRISPR_53_without_IFN-bparental_HN6_without_IFNb" "CRISPR_52_with_IFNb-parental_HN6_with_IFNb" "CRISPR_52_without_IFNb-parental_HN6_without_IFNb")
fi

frag_length_filter="1000"

# set dirs
output_dir="/data/CUTRUN/analysis/$project_id/bigwig"
project_dir="/data/CUTRUN/analysis/$project_id"
pipeliner_dir=$project_dir/$dir_loc

track_dir="/data/CUTRUN/datashare/CCBR/${project_id}/bigwig"
track_info="$project_dir/track_info.txt"

# create dirs
if [[ ! -d $output_dir ]]; then mkdir -p $output_dir; fi
if [[ ! -d $track_dir ]]; then mkdir -p $track_dir; fi
if [[ ! -f $track_info ]]; then touch $track_info; fi

# references:
if [[ $reference == "hg38" ]]; then
    fa="/data/CCBR_Pipeliner/db/PipeDB/Indices/hg38_basic/hg38.fa"
    regions="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"
else
    fa="/data/CCBR_Pipeliner/db/PipeDB/Indices/mm10_basic/mm10.fa"
    regions="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY"
fi

# create the length file if it doesnt exist
ref_length="$output_dir/ref_length"
if [[ ! -f $ref_length ]]; then
    echo "--creating ref length files"
    ln -s $fa $output_dir/$reference.fa
    samtools faidx $output_dir/$reference.fa
    cut -f1,2 $output_dir/$reference.fa.fai > $ref_length
else
    echo "Length file was already created"
fi

############################################################################################################
# run sample 
############################################################################################################
# create bigwig files, track files
sample_list="samples_list.txt"
IFS=$'\n' read -d '' -r -a sample_list < "$project_dir/$sample_list"
for sample_id in ${sample_list[@]}; do
    echo "--$sample_id"
    input_bam="$pipeliner_dir/bams/$sample_id.star_rg_added.sorted.dmark"
    sorted_bam="$output_dir/$sample_id.si.bam"
    
    bam_filtered="$output_dir/$sample_id.filt.bam"
    bed_file="$output_dir/$sample_id.bed"
    clean_bed_file="$output_dir/$sample_id.clean.bed"
    fragment_bed="$output_dir/${sample_id}.fragment.bed"
    bed_graph="$output_dir/${sample_id}.bg"
    bw="$output_dir/${sample_id}_fragment.bw"

    link_loc="${track_dir}/${sample_id}.bigwig"

    if [[ ! -f $sorted_bam ]]; then
        echo "----creating bam_links"
        ln -s $input_bam.bam $sorted_bam
        ln -s $input_bam.bai $sorted_bam.bai
    fi

    if [[ ! -f $bw ]]; then
        # filter bam for chrom regions
        echo "----filtering bam"
        samtools view -b $sorted_bam $regions | samtools sort -n -o $bam_filtered
        
        # convert
        echo "----creating bed"
        bedtools bamtobed -bedpe -i $bam_filtered > $bed_file

        # filter
        echo "----filtering bed"
        awk -v fl=$frag_length_filter '{{ if ($1==$4 && $6-$2 < fl) {{print $0}}}}' $bed_file > $clean_bed_file
        
        # format
        echo "----formatting bed"
        cut -f 1,2,6 $clean_bed_file | LC_ALL=C sort -k1,1 -k2,2n -k3,3n > $fragment_bed

        # scale bed
        echo "----scale bed"
        bedtools genomecov -i $fragment_bed -bg -scale 1 -g $ref_length > $bed_graph

        # create bigwig
        echo "----create bg"
        bedGraphToBigWig $bed_graph $ref_length $bw

    elif [[ -f $bed_graph ]]; then
        echo "----performing cleanup"
        rm $sorted_bam $sorted_bam.bai $bam_filtered $bed_file $clean_bed_file $fragment_bed $bed_graph
    fi
    
    #create hard links
    if [[ ! -f $link_loc ]]; then
        echo "----creating links"
        echo "track type=bigWig bigDataUrl=https://hpc.nih.gov/~CUTRUN/CCBR/${project_id}/bigwig/${sample_id}.bigwig name=${sample_id} description=${sample_id} visibility=full autoScale=off maxHeightPixels=128:30:1 viewLimits=1:120 color=65,105,225" >> $track_info
        ln $bw $link_loc
    fi
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
    complete_sample_id="DEG_${sample_id}_0.5_0.5"
    
    #create hard links
    source_loc="${analysis_dir}/peaks/contrasts/${complete_sample_id}/${complete_sample_id}_fragmentsbased_diffresults.bigbed "
    link_loc="${track_dir}/bigbed/${complete_sample_id}_fragmentsbased_diffresults.bigbed"
    if [[ ! -f $link_loc ]]; then ln $source_loc $link_loc; fi

    # echo track info
    echo "track name=${sample_id}_${peak_type} bigDataUrl=https://hpc.nih.gov/~CUTRUN/CCBR/${project_id}/bigbed/${complete_sample_id}_fragmentsbased_diffresults.bigbed type=bigBed itemRgb=On" >> $track_info
}

# for sample_id in ${deg_list[@]}; do
#     run_comparison_tracks $sample_id
# done
# change permissions
chmod -R 777 $track_dir