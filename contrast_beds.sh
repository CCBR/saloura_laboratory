# input files
sample_file="/data/CUTRUN/analysis/CS029758/r_analysis_230129/HN6_H4K20me3_IFNb_single_sample_peaks.txt"
contrast_file="/data/CUTRUN/analysis/CS029758/r_analysis_230129/DESEQ2_res_53_H4K20me3_IFNb_vs_HN6_H4K20me3_IFNb.csv"
sig_file="/data/CUTRUN/analysis/CS029758/r_analysis_230129/sig_HN6_only.csv"

# output files
sample_output="/data/CUTRUN/analysis/CS029758/r_analysis_230129/beds/HN6_H4K20me3_IFNb.bed"
contrast_output="/data/CUTRUN/analysis/CS029758/r_analysis_230129/beds/53_H4K20me3_IFNb_vs_HN6_H4K20me3_IFNb.bed"
sig_output="/data/CUTRUN/analysis/CS029758/r_analysis_230129/beds/HN6_H4K20me3_IFNb_sig.bed"

# bed files
hn6_bed="/data/CUTRUN/analysis/CS029758/r_analysis_230129/beds/HN6_only.bed"
s53_bed="/data/CUTRUN/analysis/CS029758/r_analysis_230129/beds/53_only.bed"
extra_bed="/data/CUTRUN/analysis/CS029758/r_analysis_230129/beds/extra_sig_HN6.bed"
sig_bed="/data/CUTRUN/analysis/CS029758/r_analysis_230129/beds/sig_HN6.bed"

# prep
rm /data/CUTRUN/analysis/CS029758/r_analysis_230129/beds/*

# sample file
echo "**sample file**"
awk '{print $1,$2,$3}' $sample_file > /data/CUTRUN/analysis/CS029758/r_analysis_230129/beds/tmp
tail -n+2 /data/CUTRUN/analysis/CS029758/r_analysis_230129/beds/tmp > $sample_output
sed -ie "s/\s/\t/g" $sample_output
head $sample_output

# contraast file
echo "**contrast file**"
awk -F"," '{print $9,$10,$11}' $contrast_file > /data/CUTRUN/analysis/CS029758/r_analysis_230129/beds/tmp
tail -n+2 /data/CUTRUN/analysis/CS029758/r_analysis_230129/beds/tmp > $contrast_output
sed -i "s/\"//g" $contrast_output
sed -ie "s/\s/\t/g" $contrast_output
head $contrast_output

# extra sig peaks
echo "**sig file**"
awk -F"," '{print $9,$10,$11}' $sig_file > /data/CUTRUN/analysis/CS029758/r_analysis_230129/beds/tmp
tail -n+2 /data/CUTRUN/analysis/CS029758/r_analysis_230129/beds/tmp > $sig_output
sed -i "s/\"//g" $sig_output
sed -i "s/\s/\t/g" $sig_output
head $sig_output

# run contrasts
module load bedtools
#find what peaks from the contrast are in HN6 - oututs HN6 peakids
bedtools intersect -a $contrast_output -b $sample_output > $hn6_bed
#find what peaks from the contrast are in 53 - outputs 53 peakids
bedtools intersect -a $contrast_output -b $sample_output -v > $53_bed
#find what sig downreg peaks are from HN6 - outputs contrast peakids
bedtools intersect -wa -a $sig_output -b $hn6_bed | sort | uniq > $sig_bed
#find what sig downreg peaks are not HN6 - outputs contrast peakids
bedtools intersect -v -a $sig_output -b $hn6_bed | sort | uniq > $extra_bed

# summarize
total=`cat $sample_file | wc -l`
echo "Total number of peaks in idiv sample: $total"
total=`cat $contrast_file | wc -l`
echo "Total number of peaks in contrast: $total"
total=`cat $sig_file | wc -l`
echo "Total number of peaks in downreg in contrast: $total"

total=`cat $sig_bed | wc -l`
echo "Total number of peaks in contrast, indiv sample, and down sig: $total"
total=`cat $extra_bed | wc -l`
echo "Total number of peaks in contrast, not in indiv sample, and down sig: $total"

if [[ -f /data/CUTRUN/analysis/CS029758/r_analysis_230129/beds/tmp* ]]; then rm /data/CUTRUN/analysis/CS029758/r_analysis_230129/beds/tmp*; fi
if [[ -f $sample_output ]]; then rm $sample_output*; fi
if [[ -f $contrast_output ]]; then rm $contrast_output*; fi
