#!/usr/bin/env bash
set -euo pipefail

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
			
			# first remove output file if it already exists or else computeMatrix
			# will silently fail, then extract PhyloP scores for each type
			for bw in "${bws[@]}"; do
				bw_base=_$(basename "${bw%%.bw}")
				if [[ $bw =~ 241-mammalian-2020v2.bigWig ]]; then
					bw_base=""
				fi

				if [[ "${extension}" -ge 1 ]]; then
					computeMatrix reference-point --referencePoint center \
					-S "${bw}" -R "${chunk}" \
					--blackListFileName "${data_dir}"/hg38-blacklist.v2.bed \
					--beforeRegionStartLength "${extension}" \
					--afterRegionStartLength "${extension}" \
					--binSize 1 -p max \
					-o "${outdir}"/"${chunk_base}"_"${extension}""${bw_base}".gz 
					rerunning with blacklisted regions included
				
				elif [[ "${extension}" -eq 0 ]]; then
					computeMatrix scale-regions -S "${bw}" \
					-R "${chunk}" -b 0 --beforeRegionStartLength 0 \
					--afterRegionStartLength 0 \
					--blackListFileName "${data_dir}"/hg38-blacklist.v2.bed \
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
		rm -r "${analysis_dir}"/"${computematrix_out_dir}"/"${in_bed_base}"_tmp
	fi

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

main(){
	# set stuff up
	if [[ $PWD =~ "hugheslab1" ]]; then
		export base_dir=/home/hugheslab1
	elif [[ $PWD =~ "home1" ]]; then
		export base_dir=/home/home1/hugheslab	 # for dc compatibility
	fi
	export proj_dir="${base_dir}"/iyellan/codebook_proj/tfbs_conservation

	export analysis_dir="${proj_dir}"/analysis
	export data_dir="${proj_dir}"/data

	# make variables visible to other functions
	export proj_dir data_dir analysis_dir

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

	# locate all the relevant files
	## triple optimized moods hits
	moods_fil="${moods_dir}"/"${ID}".bed # updated to all motif hits
	if [[ -f "${moods_fil}" ]]; then
		cp "${moods_fil}" "${moods_tmp_dir}"/"${ID}"_"${TF}"_"${construct}"_"${moods_suffix}".bed
		moods_fil="${moods_tmp_dir}"/"${ID}"_"${TF}"_"${construct}"_"${moods_suffix}".bed
	else
		echo "${moods_fil}" does not exist, skipping
		return
	fi
	
	# intersect triple optimized file with repeatmasker annotations
	if [[ $PWD =~ "hugheslab1" ]]; then
		repeat_fil="${base_dir}"/abrechalov/genome_annotations/hg38_repeats.bed
	elif [[ $PWD =~ "home1" ]]; then
		repeat_fil="${base_dir}"/iyellan/codebook_proj/hg38_repeats.bed # for dc compatibility
	fi
	bedtools intersect -a "${moods_fil}" -b "${repeat_fil}" -f 0.1 -wao \
	> "${analysis_dir}/${computematrix_out_dir}/${TF}_${construct}_moods_all_PWMs_repeatmasker.bed"

	# remove motif hits that overlap with coding regions
	cds_fil="${base_dir}"/iyellan/data/hg38.knownGene.CDS.bed
	
	cp "${moods_fil}" "${moods_fil%%.bed}"_w_coding.bed
	bedtools intersect -a "${moods_fil%%.bed}"_w_coding.bed -b "${cds_fil}" \
	-v > "${moods_fil}"

	# run deeptools with all pwm hits
	for bw in "${bws[@]}"; do
		bw_base=$(basename "${bw%%.bw}")
		if [[ $bw_base =~ 241-mammalian-2020v2.bigWig ]]; then
			bw_base=""
		fi
		mkdir -p "${analysis_dir}"/"${computematrix_out_dir}"/"${bw_base}"/hmaps/"${TF}"_"${construct}"
	done

	run_deeptools "${moods_fil}" "${bws[@]}"
}

# use case statement to parse arguments
while getopts ":i:r" opt; do
	case ${opt} in
		r )
			export redo=T
			;;
		\? )
			echo "Usage: cmd [-r]"
			exit 1
			;;
	esac
done

######### SETUP
## parse tsv of best motifs with while loop, assign variables
## set up various out directories
if [[ $PWD =~ "hugheslab1" ]]; then
	export base_dir=/home/hugheslab1
elif [[ $PWD =~ "home1" ]]; then
	export base_dir=/home/home1/hugheslab	 # for dc compatibility
fi
proj_dir="${base_dir}"/iyellan/codebook_proj/tfbs_conservation
data_dir="${proj_dir}"/data
export analysis_dir="${proj_dir}"/analysis
motifs_list="${proj_dir}"/data/best_motifs_061023_tripOpt.txt

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

export th_per_job=12
export -f run_deeptools main

# dealing with multiple moods scan types now
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

# go through list of TF motifs, submit deeptools extraction jobs in batches
while IFS=$'\t' read -r TF construct type ID; do
	if [[ -f "${moods_dir}"/"${ID}".bed ]]; then
		moods_fil="${moods_dir}"/"${ID}".bed # updated to all motif hits
		echo main "${TF}" "${construct}" \
		"'${type}'" "${ID}" "${th_per_job}" $(find "${data_dir}" -name 241-mammalian-2020v2.bigWig)
	fi
done < <(tail -n +2 "${motifs_list}") > ~/conservation_commands.txt
if [ "$(wc -l ~/conservation_commands.txt | cut -f 1 -d " ")" -gt 10 ]; then
	split -l 10 ~/conservation_commands.txt ~/conservation_commands_
else
	mv ~/conservation_commands.txt ~/conservation_commands_aa
fi 
for f in ~/conservation_commands_*; do
	submitjob -E isaac.yellan95@gmail.com -w 24 -m "${mem}" -c "${th_per_job}" -f "${f}" 
	sleep 10 # I don't remebmer exactly why there's a sleep here
done
