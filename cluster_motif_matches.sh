#!/bin/bash

proj_dir=/home/hugheslab1/iyellan/codebook_proj/tfbs_conservation
tfbs_set=$1

# just conserved
for bed in "${proj_dir}"/analysis/computematrix_out/all_PWM_hits/hmaps/*_NA/*_"${tfbs_set}"_conserved_sites.bed*; do
    if [[ "${bed}" =~ ".gz" ]]; then
        gunzip -f "${bed}"
        bed="${bed%%.gz}"
    fi
    TF=$(basename "${bed}" | sed -E 's/([A-Z0-9a-z]+)_trip_opt.*/\1/g')
    awk -v t="$TF" 'BEGIN{OFS=FS="\t"} {print $0, t"_"$1"_"$2"_"$3"_"$6}' "$bed"
done | sort -k 1,1 -k2,2n \
> "${proj_dir}"/analysis/computematrix_out/all_PWM_hits/hmaps/all_"${tfbs_set}"_conserved.bed
bedtools merge -i "${proj_dir}"/analysis/computematrix_out/all_PWM_hits/hmaps/all_"${tfbs_set}"_conserved.bed \
-c 6,7 -o collapse,collapse -delim ";" | awk 'BEGIN{OFS=FS="\t"} {print $0, NR}' \
> "${proj_dir}"/analysis/computematrix_out/all_PWM_hits/hmaps/all_"${tfbs_set}"_conserved_clustered.bed
## to find motif matches entirely overlapped by another
bedtools intersect -a "${proj_dir}"/analysis/computematrix_out/all_PWM_hits/hmaps/all_"${tfbs_set}"_conserved.bed \
-b "${proj_dir}"/analysis/computematrix_out/all_PWM_hits/hmaps/all_"${tfbs_set}"_conserved.bed \
-f 1.0 -wao | awk 'BEGIN{OFS=FS="\t"} $7!=$14' \
> "${proj_dir}"/analysis/computematrix_out/all_PWM_hits/hmaps/all_"${tfbs_set}"_conserved_complete_overlaps.bed
## for 100 bp distances
bedtools merge -i "${proj_dir}"/analysis/computematrix_out/all_PWM_hits/hmaps/all_"${tfbs_set}"_conserved.bed \
-c 6,7 -d 100 -o collapse,collapse -delim ";" | awk 'BEGIN{OFS=FS="\t"} {print $0, NR}' \
> "${proj_dir}"/analysis/computematrix_out/all_PWM_hits/hmaps/all_"${tfbs_set}"_conserved_clustered_100.bed
## for 1000 bp distances
bedtools merge -i "${proj_dir}"/analysis/computematrix_out/all_PWM_hits/hmaps/all_"${tfbs_set}"_conserved.bed \
-c 6,7 -d 1000 -o collapse,collapse -delim ";" | awk 'BEGIN{OFS=FS="\t"} {print $0, NR}' \
> "${proj_dir}"/analysis/computematrix_out/all_PWM_hits/hmaps/all_"${tfbs_set}"_conserved_clustered_1000.bed

# conserved and non-conserved
for bed in "${proj_dir}"/analysis/computematrix_out/all_PWM_hits/hmaps/*_NA/*_"${tfbs_set}"_*conserved_sites.bed*; do
    if [[ "${bed}" =~ ".gz" ]]; then
        gunzip -f "${bed}"
        bed="${bed%%.gz}"
    fi
    TF=$(basename "${bed}" | sed -E 's/([A-Z0-9a-z]+)_trip_opt.*/\1/g')
    awk -v t="$TF" 'BEGIN{OFS=FS="\t"} {print $0, t"_"$1"_"$2"_"$3"_"$6}' "$bed"
done | sort -k 1,1 -k2,2n \
> "${proj_dir}"/analysis/computematrix_out/all_PWM_hits/hmaps/all_"${tfbs_set}"_conserved+unconserved.bed
bedtools merge \
-i "${proj_dir}"/analysis/computematrix_out/all_PWM_hits/hmaps/all_"${tfbs_set}"_conserved+unconserved.bed \
-c 6,7 -o collapse,collapse -delim ";" | awk 'BEGIN{OFS=FS="\t"} {print $0, NR}' \
> "${proj_dir}"/analysis/computematrix_out/all_PWM_hits/hmaps/all_"${tfbs_set}"_conserved+unconserved_clustered.bed
bedtools merge -d 200 \
-i "${proj_dir}"/analysis/computematrix_out/all_PWM_hits/hmaps/all_"${tfbs_set}"_conserved+unconserved.bed \
-c 6,7 -o collapse,collapse -delim ";" | awk 'BEGIN{OFS=FS="\t"} {print $0, NR}' \
> "${proj_dir}"/analysis/computematrix_out/all_PWM_hits/hmaps/all_"${tfbs_set}"_conserved+unconserved_clustered_200.bed
bedtools merge -d 1000 \
-i "${proj_dir}"/analysis/computematrix_out/all_PWM_hits/hmaps/all_"${tfbs_set}"_conserved+unconserved.bed \
-c 6,7 -o collapse,collapse -delim ";" | awk 'BEGIN{OFS=FS="\t"} {print $0, NR}' \
> "${proj_dir}"/analysis/computematrix_out/all_PWM_hits/hmaps/all_"${tfbs_set}"_conserved+unconserved_clustered_1000.bed

