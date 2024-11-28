#!/usr/bin/env bash
set -euo pipefail

ght_intersect(){
	chip_moods_fil="${1}"
	ght_fil="${2}"

	out_base=$(basename "${chip_moods_fil%%.bed}")

	## intersect chip_pwm bed with ght bed
	bedtools intersect -a <(cut -f1,2,3,4,5,6 "${chip_moods_fil}") \
	-b <(tail -n +2 "${ght_fil}") -f 1 -wao \
	> "${analysis_dir}"/pwm_chip_ght_intersect/"${out_base}"_ALL.bed

	## split ght results into: 
	### PWMs with ChIP + GHT overlap
	awk -F"\t" 'BEGIN{OFS="\t"} $7!="." {print $0}' \
	"${analysis_dir}"/pwm_chip_ght_intersect/"${out_base}"_ALL.bed \
	> "${analysis_dir}"/pwm_chip_ght_intersect/"${out_base}"_ght.bed
	### PWMs with ChIP but no GHT overlap
	awk -F"\t" 'BEGIN{OFS="\t"} $7=="." {print $0}' \
	"${analysis_dir}"/pwm_chip_ght_intersect/"${out_base}"_ALL.bed \
	> "${analysis_dir}"/pwm_chip_ght_intersect/"${out_base}"_noght.bed
}

pick_top_motif(){
	in_bed="${1}"
	min_summit="${2}" # T or F
	in_bed_dir="$(dirname "${in_bed}")"
	in_bed_base=$(basename "${in_bed%%.bed}")

	# group by column 10, pick row with highest score in column 5
	awk 'BEGIN{FS=OFS="\t"} {if ($5 > max[$10]) {max[$10]=$5; row[$10]=$0}} END {
		for (i in row) print row[i]}' "${in_bed}" \
		> "${in_bed_dir}"/"${in_bed_base}"_maxpwm.bed

	# alternatively, group by column 10, and pick the row with the smallest difference
	# between column 2 and 16, or between 3 and 16
	if [[ "$min_summit" == T ]]; then
		Rscript "${proj_dir}"/scripts/pwm_near_summit.R "${in_bed}"
	fi
}

run_deeptools(){
	in_bed="${1}"
	in_bed_base=$(basename "${in_bed%%.bed}")
	bws=("${@:2}")
	
	# check that the bed file has any rows; if not, skip
	if [[ $(wc -l < "${in_bed}") -eq 0 ]]; then
		echo "${in_bed}" has no rows, skipping
		return
	fi

	# check that the merged outfile doesn't already exist or that redo isn't set; if so, skip
	if [[ -f "${analysis_dir}"/"${computematrix_out_dir}"/"${in_bed_base}"_Matrixout.gz && "${redo}" != T ]]; then
		echo "${analysis_dir}/${computematrix_out_dir}/${in_bed_base}_Matrixout.gz already exists, skipping"
	else 
		# split in_bed into multiple files 5000 lines long, put in a tmp dir 
		mkdir -p "${analysis_dir}"/"${computematrix_out_dir}"/"${in_bed_base}"_tmp
		split -l 5000 "${in_bed}" "${analysis_dir}"/"${computematrix_out_dir}"/"${in_bed_base}"_tmp/

		for bw in "${bws[@]}"; do
			bw_base=$(basename "${bw%%.bw}")
			if [[ $bw_base =~ 241-mammalian-2020v2.bigWig ]]; then
				bw_base=""
			fi
			mkdir -p "${analysis_dir}"/"${computematrix_out_dir}"/"${bw_base}"/hmaps
		done

		process_chunk() {
   			chunk="$1"
			extension="$2"
			motif_length="$3"
			bws=("${@:4}")
   			chunk_base=$(basename "${chunk}")
			outdir=$(dirname "${chunk}")
			
			# first remove output file if it already exists or else 
			# computeMatrix will silently fail
			# extract PhyloP scores for each type
			for bw in "${bws[@]}"; do
				bw_base=_$(basename "${bw%%.bw}")
				if [[ $bw =~ 241-mammalian-2020v2.bigWig ]]; then
					bw_base=""
				fi

				if [[ "${extension}" -ge 1 ]]; then
					# computeMatrix reference-point --referencePoint center \
					# -S "${bw}" -R "${chunk}" \
					# --blackListFileName "${data_dir}"/hg38-blacklist.v2.bed \
					# --beforeRegionStartLength "${extension}" \
					# --afterRegionStartLength "${extension}" \
					# --binSize 1 -p max \
					# -o "${outdir}"/"${chunk_base}"_"${extension}""${bw_base}".gz 
					# rerunning with blacklisted regions included
					computeMatrix reference-point --referencePoint center \
					-S "${bw}" -R "${chunk}" \
					--beforeRegionStartLength "${extension}" \
					--afterRegionStartLength "${extension}" \
					--binSize 1 -p max \
					-o "${outdir}"/"${chunk_base}"_"${extension}""${bw_base}".gz 
				
				elif [[ "${extension}" -eq 0 ]]; then
					# computeMatrix scale-regions -S "${bw}" \
					# -R "${chunk}" -b 0 --beforeRegionStartLength 0 \
					# --afterRegionStartLength 0 \
					# --blackListFileName "${data_dir}"/hg38-blacklist.v2.bed \
					# --binSize 1 -p max \
					# --regionBodyLength "${motif_length}" \
					# -o "${outdir}"/"${chunk_base}"_"${extension}""${bw_base}".gz
					computeMatrix scale-regions -S "${bw}" \
					-R "${chunk}" -b 0 --beforeRegionStartLength 0 \
					--afterRegionStartLength 0 \
					--binSize 1 -p max \
					--regionBodyLength "${motif_length}" \
					-o "${outdir}"/"${chunk_base}"_"${extension}""${bw_base}".gz
				fi
			done
		}

		# Export the function for parallel processing
		export -f process_chunk

		# run process_chunk for 50bp and 0bp extensions
		# Use find and parallel to process each chunk in parallel
		motif_length=$(awk -F"\t" '{print $3-$2}' "${in_bed}" | sort | uniq)
		find "${analysis_dir}"/"${computematrix_out_dir}"/"${in_bed_base}"_tmp/ \
		-regex ".*/[a-z]+" -regextype posix-extended -type f \
		| parallel --jobs $((th_per_job/2)) process_chunk {} 50 \
		"${motif_length}" "${bws[@]}"
		find "${analysis_dir}"/"${computematrix_out_dir}"/"${in_bed_base}"_tmp/ \
		-regex ".*/[a-z]+" -regextype posix-extended -type f \
		| parallel --jobs $((th_per_job/2)) process_chunk {} 0 \
		"${motif_length}" "${bws[@]}"

		# merge the chunks back together
		for bw in "${bws[@]}"; do
			bw_base_=_"$(basename "${bw%%.bw}")"
			bw_base=$(basename "${bw%%.bw}")
			if [[ $bw_base =~ 241-mammalian-2020v2.bigWig ]]; then
					bw_base=""
					bw_base_=""
				fi
			computeMatrixOperations rbind \
			-m "${analysis_dir}"/"${computematrix_out_dir}"/"${in_bed_base}"_tmp/*_50"${bw_base_}".gz \
			-o "${analysis_dir}"/"${computematrix_out_dir}"/"${bw_base}"/"${in_bed_base}"_50_Matrixout.gz
			computeMatrixOperations rbind \
			-m "${analysis_dir}"/"${computematrix_out_dir}"/"${in_bed_base}"_tmp/*_0"${bw_base_}".gz \
			-o "${analysis_dir}"/"${computematrix_out_dir}"/"${bw_base}"/"${in_bed_base}"_0_Matrixout.gz

			# add column with TF name to the matrix
			gunzip -c "${analysis_dir}"/"${computematrix_out_dir}"/"${bw_base}"/"${in_bed_base}"_50_Matrixout.gz \
			| awk 'BEGIN{FS=OFS="\t"} {print $0, "'${TF}'"}' | gzip \
			> "${analysis_dir}"/"${computematrix_out_dir}"/"${bw_base}"/"${in_bed_base}"_50_Matrixout_TF.gz
			gunzip -c "${analysis_dir}"/"${computematrix_out_dir}"/"${bw_base}"/"${in_bed_base}"_0_Matrixout.gz \
			| awk 'BEGIN{FS=OFS="\t"} {print $0, "'${TF}'"}' | gzip \
			> "${analysis_dir}"/"${computematrix_out_dir}"/"${bw_base}"/"${in_bed_base}"_0_Matrixout_TF.gz

			# clean up
			# 
		done
		# rm -r "${analysis_dir}"/"${computematrix_out_dir}"/"${in_bed_base}"_tmp
	fi

	# if [[ ! -f "${analysis_dir}"/computematrix_out/hmaps/"${TF}"_"${construct}"/"${in_bed_base}"_50_hmap.pdf \
	# && ! -f "${analysis_dir}"/computematrix_out/hmaps/"${TF}"_"${construct}"/"${in_bed_base}"_0_hmap.pdf ]]; then
	for bw in "${bws[@]}"; do
		bw_base=$(basename "${bw%%.bw}")
		if [[ $bw_base =~ 241-mammalian-2020v2.bigWig ]]; then
			bw_base=""
		fi
		# plot the heatmap
		plotHeatmap \
		-m "${analysis_dir}"/"${computematrix_out_dir}"/"${bw_base}"/"${in_bed_base}"_50_Matrixout.gz \
		-out "${analysis_dir}"/"${computematrix_out_dir}"/"${bw_base}"/hmaps/"${TF}"_"${construct}"/"${in_bed_base}"_50_hmap.pdf \
		--outFileSortedRegions "${analysis_dir}"/"${computematrix_out_dir}"/"${bw_base}"/hmaps/"${in_bed_base}"_50_SortedRegions.bed \
		--outFileNameMatrix "${analysis_dir}"/"${computematrix_out_dir}"/"${bw_base}"/hmaps/"${in_bed_base}"_50_NameMatrix.gz \
		--dpi 400 --colorMap PiYG --zMin -8 --zMax 8 --xAxisLabel "motif distance (bp)" \
		--refPointLabel "motif center"
		plotHeatmap \
		-m "${analysis_dir}"/"${computematrix_out_dir}"/"${bw_base}"/"${in_bed_base}"_0_Matrixout.gz \
		-out "${analysis_dir}"/"${computematrix_out_dir}"/"${bw_base}"/hmaps/"${TF}"_"${construct}"/"${in_bed_base}"_0_hmap.pdf \
		--outFileSortedRegions "${analysis_dir}"/"${computematrix_out_dir}"/"${bw_base}"/hmaps/"${in_bed_base}"_0_SortedRegions.bed \
		--outFileNameMatrix "${analysis_dir}"/"${computematrix_out_dir}"/"${bw_base}"/hmaps/"${in_bed_base}"_0_NameMatrix.gz \
		--dpi 400 --colorMap PiYG --zMin -8 --zMax 8 --xAxisLabel "motif position" \
		--startLabel "Motif Start" --endLabel "Motif End" --linesAtTickMarks
		gzip -f "${analysis_dir}"/"${computematrix_out_dir}"/"${bw_base}"/hmaps/"${in_bed_base}"_50_SortedRegions.bed
		gzip -f "${analysis_dir}"/"${computematrix_out_dir}"/"${bw_base}"/hmaps/"${in_bed_base}"_0_SortedRegions.bed
	done
}

intersect_fils(){
	moods_fil="${1}"
	TF="${2}"
	construct="${3}"
	type="${4}"
	ID="${5}"

	#moods_dir=/home/hughespub/Codebook_ALL_data_combed/PWMs/ChIP_Optimized_MOODS ## intersecting myself, just ChIP optimized
	# moods_fil="${moods_dir}"/"${ID}"HC_"${TF}"_"${construct}"_MOODS_COpt_hg38.bed
	moods_fil="${moods_dir}"/"${ID}".bed

	chip_fil="${opt_chip_dir}"/"${TF}.bed"
	# if chip file doesn't exist, use file from hard_chip_dir
	mkdir -p "${data_dir}"/tmp_chip_bed
	if [ ! -f "${chip_fil}" ]; then
		tmp_chip_bed="${data_dir}/tmp_chip_bed/${TF}_merged.narrowPeak.bed"
		tail -n +2 "${hard_chip_dir}"/"${TF}"_merged.narrowPeak \
		| cut -f1,2,3,4,5,10 > "${tmp_chip_bed}" # coordinates, name, score, summit
		chip_fil="${tmp_chip_bed}"
	else
		chip_fil_base=$(basename "${chip_fil}")
		# coordinates, name, score, summit
		cut -f1,2,3,4,5,10 "${chip_fil}" > "${data_dir}"/tmp_chip_bed/"${chip_fil_base}"
		chip_fil="${data_dir}"/tmp_chip_bed/"${chip_fil_base}"
	fi
	
	ght_fil="${opt_ght_dir}"/"${TF}".bed
	# if ght file doesn't exist, use file from hard_ght_dir 
	if [ ! -f "${ght_fil}" ]; then
		ght_fil="${hard_ght_dir}"/"${TF}"_"${construct}".bed
	fi

	# check if any of the files are empty or exist; if so, write to a file
	# if moods or chip file doesn't exist, skip this TF
	if [[ "$(wc -l "${moods_fil}")" -eq 0 || ! -f "${moods_fil}" \
	|| $(wc -l "${opt_chip_dir}"/"${TF}.bed") -eq 0 || ! -f "${opt_chip_dir}"/"${TF}.bed" ]]; then
		if [ "$(wc -l "${moods_fil}")" -eq 0 ] || [ ! -f "${moods_fil}" ]; then
		moods_status="missing_moods"
		fi
		if [ "$(wc -l "${opt_chip_dir}"/"${TF}.bed")" -eq 0 ] || [ ! -f "${opt_chip_dir}"/"${TF}.bed" ]; then
		chip_status="missing_chip"
		fi
		
		printf "%s\t%s\t%s\t%s\t%s\n" "${moods_fil}" "${opt_chip_dir}"/"${TF}.bed" "${ght_fil}" \
		"${moods_status}" "${chip_status}" \
		>> "${analysis_dir}"/TFs_with_missing_empty_files.txt
		return
	fi

	# write the files being used to a file
	printf "%s\t%s\t%s\t%s\t%s\n" "${TF}" "${ID}" "${moods_fil}" "${chip_fil}" "${ght_fil}" \
	>> "${analysis_dir}"/files_used.txt

	# intersect MOODS hits with ChIP peaks and GHT peaks
	bedtools intersect -a "${moods_fil}" -b "${chip_fil}" -b "${ght_fil}" -f 1 -wao -filenames \
	> "${analysis_dir}"/pwm_chip_ght_intersect/"${TF}"_"${construct}".tsv

	# parse results, output into files with different intersection results
	Rscript "${proj_dir}"/scripts/process_intersections.R \
	"${analysis_dir}"/pwm_chip_ght_intersect/"${TF}"_"${construct}".tsv \
	"${analysis_dir}"
	# bedtools intersect -a "${moods_fil}" -b "${chip_fil}" -f 1 -wa -wb \
	# > ${analysis_dir}/pwm_chip_intersect/"${TF}"_"${construct}"_moods_chip.bed
	
	# bedtools intersect -a "${moods_fil}" -b <(cut -f1,2,3,4,5,6 "${chip_fil}") \
	# -b "${ght_fil}" -f 1 -wa -wb \
	# -f 1 -wao \
	# > ${analysis_dir}/pwm_chip_intersect/"${TF}"_"${construct}"_moods_chip.bed
	
	# pick top motif (ChIP has summit, GHT does not, so can use minsummit)
	# pick_top_motif ${analysis_dir}/pwm_chip_intersect/"${TF}"_"${construct}"_moods_chip.bed T
	
	# intersect both peak sets with GHT peaks
	# ght_intersect "${analysis_dir}/pwm_chip_intersect/${TF}_${construct}_moods_chip_maxpwm.bed" \
	# "${ght_fil}"
	# ght_intersect "${analysis_dir}/pwm_chip_intersect/${TF}_${construct}_moods_chip_minsummit.bed" \
	# "${ght_fil}"
	
	# also intersect just the motif hits with GHT peaks
	# bedtools intersect -a "${moods_fil}" -b "${ght_fil}" -f 1 -wa -wb \
	# > ${analysis_dir}/pwm_ght_intersect/"${TF}"_"${construct}"_moods_ght.bed
	# pick_top_motif ${analysis_dir}/pwm_ght_intersect/"${TF}"_"${construct}"_moods_ght.bed F
}

main(){
	# set stuff up
	if [[ $PWD =~ "hugheslab1" ]]; then
		export base_dir=/home/hugheslab1
	elif [[ $PWD =~ "home1" ]]; then
		export base_dir=/home/home1/hugheslab	 # for dc compatibility
	fi
	export proj_dir="${base_dir}"/iyellan/codebook_proj/tfbs_conservation
	# proj_dir=/home/hugheslab1/iyellan/codebook_proj/tfbs_conservation
	# moods_dir=/home/hugheslab1/afathi/ChIP_GHT_Jaccard/Cheat_Set/Triple_Optimized_MOODS_BED ## updated to triple optimized results
	opt_chip_dir=/home/hughespub/Codebook_ALL_data_combed/ChIP_seq/ChIP_Peaks_Jaccard_Optimized
	hard_chip_dir=/home/hughespub/Codebook_ALL_data_combed/ChIP_seq/Merged_P100_narrowPeaks
	opt_ght_dir=/home/hughespub/Codebook_ALL_data_combed/GHT_SELEX/GHT_Peaks_Jaccard_Optimized
	hard_ght_dir=/home/hughespub/Codebook_ALL_data_combed/GHT_SELEX/HamedPeaks/Filtered_S30
	export analysis_dir="${proj_dir}"/analysis
	export data_dir="${proj_dir}"/data

	# make variables visible to other functions
	export proj_dir moods_dir opt_chip_dir hard_chip_dir opt_ght_dir hard_ght_dir \
	data_dir analysis_dir

	TF="${1}"
	construct="${2}"
	construct=NA # the new Jaccard optimized sets of everything don't have the construct names in them
	type="${3}"
	ID="${4}"
	export th_per_job="${5}"
	# assign argument 6 and beyond to bw variable
	bws=("${@:6}")

	# deal with the multiple versions of the moods results
	if [[ "$moods_dir" =~ Triple_Optimized_MOODS_BED_HC ]]; then
		export computematrix_out_dir=computematrix_out/trip_opt_moods_HC
		mkdir -p "${analysis_dir}"/"${computematrix_out_dir}"
		moods_tmp_dir="${analysis_dir}"/trip_opt_moods_HC
		mkdir -p "${moods_tmp_dir}"
		moods_suffix=MOODS_tripOpt_HC_hg38
	elif [[ "$moods_dir" =~ Best_MOODS_BED ]]; then
		export computematrix_out_dir=computematrix_out/all_PWM_hits
		mkdir -p "${analysis_dir}"/"${computematrix_out_dir}"
		moods_tmp_dir="${analysis_dir}"/all_PWM_hits
		mkdir -p "${moods_tmp_dir}"
		moods_suffix=MOODS_all_PWM_hits_hg38
	else
		moods_tmp_dir="${analysis_dir}"/trip_opt_moods
		moods_suffix=MOODS_tripOpt_hg38
		export computematrix_out_dir=computematrix_out
	fi

	#echo "${TF}" "${construct}" "${type}" "${ID}"

	# locate all the relevant files
	## triple optimized moods hits
	# moods_fil="${moods_dir}"/"${ID}"_"${TF}".bed
	moods_fil="${moods_dir}"/"${ID}".bed # updated to all motif hits
	if [[ -f "${moods_fil}" ]]; then
		cp "${moods_fil}" "${moods_tmp_dir}"/"${ID}"_"${TF}"_"${construct}"_"${moods_suffix}".bed
		moods_fil="${moods_tmp_dir}"/"${ID}"_"${TF}"_"${construct}"_"${moods_suffix}".bed
	else
		echo "${moods_fil}" does not exist, skipping
		return
	fi
	
	if [[ "${run_intersects}" == T ]]; then
		intersect_fils "${moods_fil}" "${TF}" "${construct}" "${type}" "${ID}"

		# run deeptools computematrix to extract PhyloP scores, using both motif 
		# selection methods and ght and non-ght intersecting peaks
		# NOTE: i think the way to go here is to use maxpwm, since ght peaks don't have summits
		# and broadly there's little difference between the two methods
	fi
	# intersect triple optimized file with repeatmasker annotations
	if [[ $PWD =~ "hugheslab1" ]]; then
		repeat_fil="${base_dir}"/abrechalov/genome_annotations/hg38_repeats.bed
	elif [[ $PWD =~ "home1" ]]; then
		repeat_fil="${base_dir}"/iyellan/codebook_proj/hg38_repeats.bed # for dc compatibility
	fi
	bedtools intersect -a "${moods_fil}" -b "${repeat_fil}" -f 0.1 \
	-wao \
	> "${analysis_dir}/${computematrix_out_dir}/${TF}_${construct}_moods_all_PWMs_repeatmasker.bed"
	# > "${analysis_dir}/${computematrix_out_dir}/${TF}_${construct}_moods_tripOpt_repeatmasker.bed"

	# remove motif hits that overlap with coding regions
	cds_fil="${base_dir}"/iyellan/data/hg38.knownGene.CDS.bed
	
	cp "${moods_fil}" "${moods_fil%%.bed}"_w_coding.bed
	bedtools intersect -a "${moods_fil%%.bed}"_w_coding.bed -b "${cds_fil}" \
	-v > "${moods_fil}"

	# run deeptools with the special triple optimized file
	for bw in "${bws[@]}"; do
		bw_base=$(basename "${bw%%.bw}")
		if [[ $bw_base =~ 241-mammalian-2020v2.bigWig ]]; then
			bw_base=""
		fi
		mkdir -p "${analysis_dir}"/"${computematrix_out_dir}"/"${bw_base}"/hmaps/"${TF}"_"${construct}"
	done
	# run_deeptools "${analysis_dir}/pwm_chip_intersect/${TF}_${construct}_moods_chip_maxpwm.bed"
	# run_deeptools "${analysis_dir}/pwm_chip_ght_intersect/${TF}_${construct}_moods_chip_maxpwm_ght.bed"
	# run_deeptools "${analysis_dir}/pwm_chip_ght_intersect/${TF}_${construct}_moods_chip_maxpwm_noght.bed"
	# run_deeptools "${analysis_dir}/pwm_ght_intersect/${TF}_${construct}_moods_nochip_ght_maxpwm.bed"
	# run_deeptools "${analysis_dir}/pwm_no_intersect/${TF}_${construct}_moods_nochip_noght.bed"
	run_deeptools "${moods_fil}" "${bws[@]}"
}

get_motif_sizes(){
	proj_dir="${base_dir}"/iyellan/codebook_proj/tfbs_conservation
	analysis_dir="${proj_dir}"/analysis
	# moods_dir=/home/hugheslab1/afathi/ChIP_GHT_Jaccard/Cheat_Set/Triple_Optimized_MOODS_BED ## updated triple optimized results
	# moods_dir_alt=/home/hughespub/Codebook_ALL_data_combed/PWMs/ChIP_Optimized_MOODS ## updated to Jaccard optimized results
	
	 ## all motif hits
	if [[ $PWD =~ "hugheslab1" ]]; then
		export moods_dir=/home/hugheslab1/afathi/ChIP_GHT_Jaccard/Cheat_Set/MOODS_job/Best_MOODS_BED
	elif [[ $PWD =~ "home1" ]]; then
		export moods_dir=/home/home1/hugheslab/iyellan/codebook_proj/Best_MOODS_BED	 # for dc compatibility
	fi
	
	TF="${1-default}"
	#construct="${2}"
	construct=NA
	ID="${2-default}"
	## NOTE: THE MOODS COORDINATES ARE 0-BASED BUT THIS ISN'T AN ISSUE FOR WHERE DEEPTOOLS EXTRACTS PHYLOP SCORES
	# moods_fil="${moods_dir}"/"${ID}"_"${TF}".bed
	# moods_fil_alt="${moods_dir_alt}"/"${ID}"HC_"${TF}"_"${construct}"_MOODS_COpt_hg38.bed
	moods_fil="${moods_dir}"/"${ID}".bed
	# moods_fil_alt="${moods_dir_alt}"/"${ID}"_"${TF}".bed
	motif_len=NA
	if [ -f "${moods_fil}" ]; then
		motif_len=$(head "${moods_fil}" | awk 'BEGIN{FS=OFS="\t"} {print $3-$2}' | sort | uniq)
	# elif [ -f "${moods_fil_alt}" ]; then
	# 	motif_len=$(awk 'BEGIN{FS=OFS="\t"} {print $3-$2}' "${moods_fil_alt}" | sort | uniq)
	fi
	printf "%s\t%s\t%s\n" "${TF}" "${ID}" "${motif_len}" \
	>> "${analysis_dir}"/motif_sizes.txt
}

# 
# export run_intersects=${1-default}
# export redo=${2-default}
# use case statement to parse arguments
while getopts ":i:r" opt; do
	case ${opt} in
		i )
			export run_intersects=T
			;;
		r )
			export redo=T
			;;
		\? )
			echo "Usage: cmd [-i] [-r]"
			exit 1
			;;
	esac
done

# parse tsv of best motifs with while loop, assign variables
## set up various out directories
if [[ $PWD =~ "hugheslab1" ]]; then
	export base_dir=/home/hugheslab1
elif [[ $PWD =~ "home1" ]]; then
	export base_dir=/home/home1/hugheslab	 # for dc compatibility
fi
proj_dir="${base_dir}"/iyellan/codebook_proj/tfbs_conservation
data_dir="${proj_dir}"/data
export analysis_dir="${proj_dir}"/analysis
# motifs_list="${proj_dir}"/data/best_motifs_061023.txt
motifs_list="${proj_dir}"/data/best_motifs_061023_tripOpt.txt
# motifs_list="${proj_dir}"/data/best_motifs_061023_tripOpt_failed.txt

mkdir -p ${analysis_dir}/pwm_chip_intersect
mkdir -p "${data_dir}"/tmp_chip_bed
mkdir -p "${analysis_dir}"/pwm_chip_ght_intersect
mkdir -p ${analysis_dir}/pwm_ght_intersect
mkdir -p ${analysis_dir}/pwm_no_intersect
mkdir -p "${analysis_dir}"/trip_opt_moods

if [[ -f "${analysis_dir}"/TFs_with_missing_empty_files.txt ]]; then
	rm "${analysis_dir}"/TFs_with_missing_empty_files.txt
fi
if [[ -f "${analysis_dir}"/files_used.txt ]]; then
	printf "%s\t%s\t%s\t%s\t%s\n" "TF" "ID" "MOODS" "ChIP" "GHT" \
	> "${analysis_dir}"/files_used.txt
fi
if [[ -f "${analysis_dir}"/motif_sizes.txt ]]; then
	printf "%s\t%s\t%s\n" "TF" "ID" "motif_len" \
	> "${analysis_dir}"/motif_sizes.txt
fi
#th=28
export th_per_job=12
export -f ght_intersect pick_top_motif run_deeptools intersect_fils main

# dealing with multiple moods scan types now
# moods_dir=/home/hugheslab1/afathi/ChIP_GHT_Jaccard/Cheat_Set/Triple_Optimized_MOODS_BED ## original triple optimized results
# export moods_dir=/home/hugheslab1/afathi/ChIP_GHT_Jaccard/Cheat_Set/Triple_Optimized_MOODS_BED_HC ## new triple optimized results
if [[ $PWD =~ "hugheslab1" ]]; then
		export moods_dir=/home/hugheslab1/afathi/ChIP_GHT_Jaccard/Cheat_Set/MOODS_job/Best_MOODS_BED
	elif [[ $PWD =~ "home1" ]]; then
		export moods_dir=/home/home1/hugheslab/iyellan/codebook_proj/Best_MOODS_BED	 # for dc compatibility
fi

if [[ "${PWD}" =~ "hugheslab1" ]]; then
	mem=6
elif [[ "${PWD}" =~ "home1" ]]; then
	mem=80
fi

while IFS=$'\t' read -r TF construct type ID; do
	# if [[ -f "${moods_dir}"/"${ID}"_"${TF}".bed ]]; then
	if [[ -f "${moods_dir}"/"${ID}".bed ]]; then
		# moods_fil="${moods_dir}"/"${ID}"_"${TF}".bed
		moods_fil="${moods_dir}"/"${ID}".bed # updated to all motif hits
		echo main "${TF}" "${construct}" \
		"'${type}'" "${ID}" "${th_per_job}" $(find "${data_dir}" -name 241-mammalian-2020v2.bigWig)
		# 
	fi
	# get_motif_sizes "${TF}" "${ID}"
done < <(tail -n +2 "${motifs_list}") > ~/conservation_commands.txt
if [ "$(wc -l ~/conservation_commands.txt | cut -f 1 -d " ")" -gt 10 ]; then
	split -l 10 ~/conservation_commands.txt ~/conservation_commands_
else
	mv ~/conservation_commands.txt ~/conservation_commands_aa
fi 
for f in ~/conservation_commands_*; do
	submitjob -E isaac.yellan95@gmail.com -w 24 -m "${mem}" -c "${th_per_job}" -f "${f}" 
	sleep 10
done
# > ~/conservation_jobid.txt; ~/scripts/failed_job_reporter_v2.sh ~/conservation_jobid.txt
#\
# #\
#

# cat "${motifs_list}" > tmp.txt; \
# submitjob -w 15 -m 5cd  -E isaac.yellan95@gmail.com -c "${th}" parallel \
# -j $((th/th_per_job)) --colsep \'\\t\' main {1} {2} {3} {4} :::: tmp.txt \
# > ~/conservation_jobid.txt; ~/scripts/failed_job_reporter_v2.sh ~/conservation_jobid.txt

#grep -E "CAMTA2|CTCF|CGGBP1" "${motifs_list}" > tmp.txt; \
# submitjob -w 3 -m 10 -c 10 parallel -j $((th/4)) \
# --colsep '\t' -a <(tail -n +2 "${motifs_list}") main {1} {2} {3} {4}
