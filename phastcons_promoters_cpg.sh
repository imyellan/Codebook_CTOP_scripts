#!/bin/bash

export proj_dir=/home/hugheslab1/iyellan/codebook_proj/tfbs_conservation
phastcons="${proj_dir}"/data/hg38.phastConsElements470way.bed
# all_TSSs=/home/hugheslab1/afathi/UCSC_Tracks/Human_Promoters/Promoters/hg38_promoters.bed
protein_coding_genes=/home/hugheslab1/afathi/UCSC_Tracks/Human_Promoters/protein_coding_genes.txt
# promoters=/home/hugheslab1/afathi/UCSC_Tracks/Human_Promoters/hg38_protein_coding_promoters.bed
promoters=/home/hughespub/Codebook_ALL_data_combed/GenomeTracks/Annotations/gencodeV44_promoters.bed # updated promoter set
cpgs_masked="${proj_dir}"/data/cpgIslandExt.bed
cpgs_unmasked="${proj_dir}"/data/cpgIslandExtUnmasked.bed
repeat_file=/home/hugheslab1/iyellan/data/hg38_ucsc_rmsk.bed
enhancers=/home/hugheslab1/afathi/EnrichmentAnalysis/ChromHMM/HEK293_s6.bed
TF_list="${proj_dir}"/data/TF_list.txt
ctop_clusters="${proj_dir}"/analysis/computematrix_out/all_PWM_hits/hmaps/all_trip_opt_moods_HC_conserved_clustered_100.bed

# prepare TSS list - changing to be +/- 50 bp from TSS to deal with CpG island annotations really close to TSS
awk -F"\t" 'BEGIN{OFS="\t"} {if($6=="+"){
    print $1,$2+949,$3-450,$4,$5,$6} 
    else{print $1,$2+451,$3-950,$4,$5,$6} 
    }' "${promoters}" > "${proj_dir}"/data/gencodeV44_TSSs.bed
all_TSSs="${proj_dir}"/data/gencodeV44_TSSs.bed

mkdir -p "${proj_dir}/analysis/phastcons_promoter_intersect"

# create sets of protein coding gene promoter, non-protein coding gene promoter, and other CpG islands
# bedtools intersect -a "${cpgs_masked}" -b "${all_TSSs}" -wao \
# > "${proj_dir}"/data/cpgIslands_masked_TSS_intersect.txt
bedtools intersect -a "${cpgs_unmasked}" -b "${all_TSSs}" -wao \
> "${proj_dir}"/data/cpgIslands_unmasked_TSS_intersect.txt
bedtools intersect -a "${all_TSSs}" -b "${cpgs_unmasked}" -wao \
> ${all_TSSs%%.bed}"_cpgs_unmasked_intersect.txt"
# classify promoters
awk -F"\t" 'BEGIN{OFS="\t"} {if ($7!="."){
    print $1,$2,$3,$4,$5,$6,"CpG_promoter"}
    else{print $1,$2,$3,$4,$5,$6,"Non_CpG_promoter"}
    }' ${all_TSSs%%.bed}"_cpgs_unmasked_intersect.txt" \
    > "${all_TSSs%%.bed}"_cpg_classified.bed
Rscript "${proj_dir}"/scripts/add_cpgs_promoters.R
promoters_cpgs=$(dirname "${all_TSSs%%.bed}")/$(basename "${promoters%%.bed}")_cpg_classified.bed

# divide CpG islands into protein-coding TSS overlapping and not
# Rscript "$proj_dir"/scripts/divide_cpg_isles.R \
# "${proj_dir}"/data/cpgIslands_masked_TSS_intersect.txt
Rscript "$proj_dir"/scripts/divide_cpg_isles.R \
"${proj_dir}"/data/cpgIslands_unmasked_TSS_intersect.txt
cpgs_masked="${proj_dir}"/data/cpgIslandExt_masked_annotated.bed
cpgs_unmasked="${proj_dir}"/data/cpgIslandExt_unmasked_annotated.bed

intersect_phastcons(){
    TF="$1"
    phastcons="$2"
    conserved_sites="${proj_dir}"/analysis/computematrix_out/all_PWM_hits/hmaps/"${TF}"_NA/"${TF}"_trip_opt_moods_HC_conserved_sites.bed.gz

    if [[ -f "${conserved_sites}" || -f "${conserved_sites%%.gz}" ]]; then
        if [[ -f "${conserved_sites}" ]]; then
            gunzip "$conserved_sites"
            conserved_sites="${conserved_sites%%.gz}"
        else
            conserved_sites="${conserved_sites%%.gz}"
        fi
        bedtools intersect -a "${conserved_sites}" -b "${phastcons}" -wao \
        > "${proj_dir}/analysis/phastcons_promoter_intersect/${TF}_phastcons_intersect.bed"
    fi

    # intersect directly with C-TOP clusters
    if [[ "${TF}" == "${ctop_clusters}" ]]; then 
        bedtools intersect -a "${ctop_clusters}" -b "${phastcons}" -wao \
        > "${proj_dir}/analysis/phastcons_promoter_intersect/phastcons_intersect_ctop_clust.bed"
    fi
}

intersect_promoters(){
    TF="$1"
    promoters="$2"
    conserved_sites="${proj_dir}"/analysis/computematrix_out/all_PWM_hits/hmaps/"${TF}"_NA/"${TF}"_trip_opt_moods_HC_conserved_sites.bed
    unconserved_sites="${proj_dir}"/analysis/computematrix_out/all_PWM_hits/hmaps/"${TF}"_NA/"${TF}"_trip_opt_moods_HC_unconserved_sites.bed

    if [[ -f "${conserved_sites}" || -f "${conserved_sites%%.gz}" ]]; then
        if [[ -f "${conserved_sites}" ]]; then
            gunzip "$conserved_sites"
            conserved_sites="${conserved_sites%%.gz}"
        else
            conserved_sites="${conserved_sites%%.gz}"
        fi
        bedtools intersect -a "${conserved_sites}" -b "${promoters}" -wao \
        > "${proj_dir}/analysis/phastcons_promoter_intersect/${TF}_promoters_intersect.bed"
    fi

    # intersect directly with C-TOP clusters
    if [[ "${TF}" == "${ctop_clusters}" ]]; then 
        bedtools intersect -a "${ctop_clusters}" -b "${promoters}" -wao \
        > "${proj_dir}/analysis/phastcons_promoter_intersect/promoters_intersect_ctop_clust.bed"
    fi

    if [[ -f "${unconserved_sites}" || -f "${unconserved_sites%%.gz}" ]]; then
        if [[ -f "$unconserved_sites" ]]; then
            gunzip "$unconserved_sites"
            unconserved_sites="${unconserved_sites%%.gz}"
        else
            unconserved_sites="${unconserved_sites%%.gz}"
        fi
        bedtools intersect -a "${unconserved_sites}" -b "${promoters}" -wao \
        > "${proj_dir}/analysis/phastcons_promoter_intersect/${TF}_promoters_intersect_uncons.bed"
    fi
}

intersect_cpgs(){
    TF="$1"
    cpgs="$2"
    conserved_sites="${proj_dir}"/analysis/computematrix_out/all_PWM_hits/hmaps/"${TF}"_NA/"${TF}"_trip_opt_moods_HC_conserved_sites.bed.gz
    
    if [[ "${cpgs}" =~ "_masked" ]]; then
        masking="masked"
    else
        masking="unmasked"
    fi

    if [[ -f "${conserved_sites}" || -f "${conserved_sites%%.gz}" ]]; then
        if [[ -f "$conserved_sites" ]]; then
            gunzip "$conserved_sites"
            conserved_sites="${conserved_sites%%.gz}"
        else
            conserved_sites="${conserved_sites%%.gz}"
        fi
        bedtools intersect -a "${conserved_sites}" -b "${cpgs}" -wao \
        > "${proj_dir}"/analysis/phastcons_promoter_intersect/"${TF}"_cpgs_"${masking}"_intersect.bed
    fi
    # intersect directly with C-TOP clusters
    if [[ "${TF}" == "${ctop_clusters}" ]]; then 
        bedtools intersect -a "${ctop_clusters}" -b "${cpgs}" -wao \
        > ${proj_dir}/analysis/phastcons_promoter_intersect/cpgs_"${masking}"_intersect_ctop_clust.bed
    fi
}

intersect_repeats() {
    # intersect tbfs with repeats
    TF="$1"
    repeat_file="$2"
    conserved_sites="${proj_dir}"/analysis/computematrix_out/all_PWM_hits/hmaps/"${TF}"_NA/"${TF}"_trip_opt_moods_HC_conserved_sites.bed.gz
    unconserved_sites="${proj_dir}"/analysis/computematrix_out/all_PWM_hits/hmaps/"${TF}"_NA/"${TF}"_trip_opt_moods_HC_unconserved_sites.bed.gz

    if [[ -f "${conserved_sites}" || -f "${conserved_sites%%.gz}" ]]; then
        if [[ -f "$conserved_sites" ]]; then
            gunzip "$conserved_sites"
            conserved_sites="${conserved_sites%%.gz}"
        else
            conserved_sites="${conserved_sites%%.gz}"
        fi
        bedtools intersect -a "${conserved_sites}" -b "${repeat_file}" -wao \
        > "${proj_dir}"/analysis/phastcons_promoter_intersect/"${TF}"_repeats_intersect.bed
        if [[ ! -s "${proj_dir}/analysis/phastcons_promoter_intersect/${TF}_repeats_intersect.bed" ]]; then
            rm "${proj_dir}"/analysis/phastcons_promoter_intersect/"${TF}"_repeats_intersect.bed
        fi
    fi

    if [[ -f "${unconserved_sites}" || -f "${unconserved_sites%%.gz}" ]]; then
        if [[ -f "$unconserved_sites" ]]; then
            gunzip "$unconserved_sites"
            unconserved_sites="${unconserved_sites%%.gz}"
        else
            unconserved_sites="${unconserved_sites%%.gz}"
        fi
        bedtools intersect -a "${unconserved_sites}" -b "${repeat_file}" -wao \
        > "${proj_dir}"/analysis/phastcons_promoter_intersect/"${TF}"_repeats_intersect_unconserved.bed
        if [[ ! -s "${proj_dir}/analysis/phastcons_promoter_intersect/${TF}_repeats_intersect_unconserved.bed" ]]; then
            rm "${proj_dir}"/analysis/phastcons_promoter_intersect/"${TF}"_repeats_intersect_unconserved.bed
        fi
    fi
    
    # intersect directly with C-TOP clusters
    if [[ "${TF}" == "${ctop_clusters}" ]]; then 
        bedtools intersect -a "${ctop_clusters}" -b "${repeat_file}" -wao \
        > "${proj_dir}"/analysis/phastcons_promoter_intersect/repeats_intersect_ctop_clust.bed
    fi
}

intersect_enhancers(){
    TF="$1"
    enhancers="$2"
    conserved_sites="${proj_dir}"/analysis/computematrix_out/all_PWM_hits/hmaps/"${TF}"_NA/"${TF}"_trip_opt_moods_HC_conserved_sites.bed.gz

    if [[ -f "${conserved_sites}" || -f "${conserved_sites%%.gz}" ]]; then
        if [[ -f "${conserved_sites}" ]]; then
            gunzip "$conserved_sites"
            conserved_sites="${conserved_sites%%.gz}"
        else
            conserved_sites="${conserved_sites%%.gz}"
        fi
        bedtools intersect -a "${conserved_sites}" -b "${enhancers}" -wao \
        > "${proj_dir}/analysis/phastcons_promoter_intersect/${TF}_enhancer_intersect.bed"
    fi

    # intersect directly with C-TOP clusters
    if [[ "${TF}" == "${ctop_clusters}" ]]; then 
        bedtools intersect -a "${ctop_clusters}" -b "${enhancers}" -wao \
        > "${proj_dir}"/analysis/phastcons_promoter_intersect/enhancer_intersect_ctop_clust.bed
    fi
}

export -f intersect_phastcons intersect_promoters intersect_cpgs intersect_repeats intersect_enhancers

if [[ ! -f "${phastcons}" ]]; then
    phastcons="${phastcons%%.gz}"
fi

# if [[ ! -f "${phastcons}" ]]; then
#     if [[ ! -f "${phastcons%%.gz}" ]]; then
#         echo "PhastCons file not found. Exiting."
#         exit 1
#     else
#         phastcons="${phastcons%%.gz}"
#     fils
# elif [[ "${phastcons}" =~ ".gz" ]]; then
#     gunzip "${phastcons}"
#     phastcons="${phastcons%%.gz}"
# fi

# intersect with individual TOPs
# parallel -j 6 -n 1 intersect_phastcons {} "${phastcons}" ::: $(cat "${TF_list}")
parallel -j 6 -n 1 intersect_promoters {} "${promoters_cpgs}" ::: $(cat "${TF_list}")
# parallel -j 6 -n 1 intersect_cpgs {} "${cpgs_masked}" ::: $(cat "${TF_list}")
# parallel -j 6 -n 1 intersect_cpgs {} "${cpgs_unmasked}" ::: $(cat "${TF_list}")
parallel -j 6 -n 1 intersect_repeats {} "${repeat_file}" ::: $(cat "${TF_list}")
# parallel -j 6 -n 1 intersect_enhancers {} "${enhancers}" ::: $(cat "${TF_list}")

# intersect with C-TOP clusters
# intersect_phastcons "${ctop_clusters}" "${phastcons}"
# intersect_promoters "${ctop_clusters}" "${promoters_cpgs}"
# intersect_cpgs "${ctop_clusters}" "${cpgs_masked}"
# intersect_cpgs "${ctop_clusters}" "${cpgs_unmasked}"
# intersect_repeats "${ctop_clusters}" "${repeat_file}"
# intersect_enhancers "${ctop_clusters}" "${enhancers}"


# also calculate CG content means for each site
bigWigAverageOverBed "${proj_dir}"/data/gc5Base.bw \
"${proj_dir}"/analysis/computematrix_out/all_PWM_hits/hmaps/conservation_all_sites_trip_opt_moods_HC.bed \
"${proj_dir}"/analysis/phastcons_promoter_intersect/conservation_all_sites_trip_opt_moods_HC_cg_mean.txt \
-bedOut="${proj_dir}"/analysis/phastcons_promoter_intersect/conservation_all_sites_trip_opt_moods_HC_cg_mean.bed

# and each cluster
bigWigAverageOverBed "${proj_dir}"/data/gc5Base.bw \
"${proj_dir}"/analysis/computematrix_out/all_PWM_hits/hmaps/all_conserved_clustered_fixed.bed \
"${proj_dir}"/analysis/phastcons_promoter_intersect/all_conserved_clustered_cg_mean.txt \
-bedOut="${proj_dir}"/analysis/phastcons_promoter_intersect/all_conserved_clustered_cg_mean.bed

# delete any empy files
find "${proj_dir}"/analysis/phastcons_promoter_intersect -type f -size 0 -delete