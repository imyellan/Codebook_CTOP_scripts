library(readr)
library(dplyr)
library(purrr)
library(ggplot2)
library(ggrepel)
library(furrr)
library(tidyr)
library(tibble)
library(dtplyr)
# options(future.fork.multithreading.enable = FALSE)
proj_dir <-
  ifelse(grepl("MacBook", Sys.info()["nodename"]),
    "~/iyellan_bc2-clust/codebook_proj/tfbs_conservation",
    "/home/hugheslab1/iyellan/codebook_proj/tfbs_conservation"
  )
all_hmaps <- F # don't make heatmaps for just conserved & unconserved separately

source(file.path(proj_dir, "scripts/conservation_analysis_lrt_fns.R"))
source(file.path(proj_dir, "scripts/conservation_analysis_fns.R"))

# hmap_outs <- ifelse(grepl("MacBook", Sys.info()["nodename"]),
#                     file.path(proj_dir,"analysis/computematrix_out/all_PWM_hits/hmaps"),
#                     commandArgs(trailingOnly = T)[1]
#                     )
hmap_outs <- file.path(proj_dir, "analysis/computematrix_out/all_PWM_hits/hmaps")
tfbs_subset_path <- ifelse(grepl("MacBook", Sys.info()["nodename"]),
  file.path(proj_dir, "data/trip_opt_moods_hc_nms.txt"),
  commandArgs(trailingOnly = T)[1]
)
if (!grepl("all_moods", tfbs_subset_path)) {
  tfbs_subset_nms <<- read_tsv(tfbs_subset_path) %>%
    mutate(TF_tfbs_nm = paste0(TF, tfbs_nm)) %>%
    pull(TF_tfbs_nm)
} else {
  tfbs_subset_nms <<- NA
}

TF_nm <- commandArgs(trailingOnly = T)[2]
pipeline_section <- commandArgs(trailingOnly = T)[3]
filt_set <- gsub(
  "_nms.txt", "",
  basename(tfbs_subset_path)
) # "trip_opt_moods_HC"

compmat_fils <-
  list.files(file.path(hmap_outs, ".."),
    pattern = "PWM.*_NA_MOODS_all_PWM_hits_hg38_50_Matrixout_TF.gz",
    full.names = TRUE
  )
dir.create(file.path(hmap_outs, "figs"), recursive = T)

lrt_sig_mats <- list.files(hmap_outs,
  recursive = T,
  # pattern = "_sig_motifs_wilcoxon_nonnorm.txt.gz",
  pattern = "_lrt_thresh_-10_2nd_d.*icVals_pass.*tsv.gz",
  full.names = T
)
lrt_nonsig_mats <- list.files(hmap_outs,
  recursive = T,
  # pattern = "_nonsig_motifs_wilcoxon_nonnorm.txt.gz",
  pattern = "_lrt_thresh_-10_2nd_d.*icVals.*fail.*tsv.gz",
  full.names = T
)
corr_sig_mats <- list.files(hmap_outs,
  recursive = T,
  # pattern = "_sig_motifs_wilcoxon_nonnorm.txt.gz",
  pattern = "_corr_thresh_.*_pass.tsv.gz",
  full.names = T
)
corr_nonsig_mats <- list.files(hmap_outs,
  recursive = T,
  # pattern = "_nonsig_motifs_wilcoxon_nonnorm.txt.gz",
  pattern = "_corr_thresh_.*_fail.tsv.gz",
  full.names = T
)
wilc_sig_mats <- list.files(hmap_outs,
  recursive = T,
  pattern = "_wilcoxon_thresh_pass.tsv.gz",
  full.names = T
)
wilc_nonsig_mats <- list.files(hmap_outs,
  recursive = T,
  pattern = "_wilcoxon_thresh_fail.tsv.gz",
  full.names = T
)

tf_nms <- read_tsv(file.path(proj_dir, "data/best_motifs_061023_tripOpt.txt"))$TF
# tf_nms <- unique(gsub("_.*", "", basename(c(corr_sig_mats, lrt_sig_mats,
#                                             lrt_nonsig_mats, corr_nonsig_mats,
#                                             wilc_sig_mats, wilc_nonsig_mats))))

tbl_import <- function(tf_nm, tbl_fil) {
  tbl_tmp <- read_tsv(tbl_fil) %>%
    mutate(across(-c(chr, id, strand), as.numeric),
      TF_tfbs_nm = paste0(tf_nm, paste(chr, start, end, sep = "_"))
    )
  if (all(!is.na(tfbs_subset_nms))) {
    tbl_tmp %>%
      filter(TF_tfbs_nm %in% tfbs_subset_nms) %>%
      select(-TF_tfbs_nm)
  } else {
    tbl_tmp %>% select(-TF_tfbs_nm)
  }
}

summarize_conservation <- function(tf_nm) {
  # print(tf_nm)
  non_pos_cols <<- c("chr", "start", "end", "id", "pwm_score", "strand", "TF")
  if (pipeline_section == "full") {
    # import results
    lrt_sig_mat <- 
      tbl_import(tf_nm, grep(paste0(tf_nm, "_"), lrt_sig_mats, value = T))
    lrt_nonsig_mat <- 
      tbl_import(tf_nm, grep(paste0(tf_nm, "_"), lrt_nonsig_mats, value = T))
    corr_sig_mat <- 
      tbl_import(tf_nm, grep(paste0(tf_nm, "_"), corr_sig_mats, value = T))
    corr_nonsig_mat <- 
      tbl_import(tf_nm, grep(paste0(tf_nm, "_"), corr_nonsig_mats, value = T))
    wilc_sig_mat <- 
      tbl_import(tf_nm, grep(paste0(tf_nm, "_"), wilc_sig_mats, value = T))
    wilc_nonsig_mat <- 
      tbl_import(tf_nm, grep(paste0(tf_nm, "_"), wilc_nonsig_mats, value = T))
    motif_positions <- colnames(wilc_nonsig_mat)[-c(1:8)]

    # summarize which sites pass which test
    print(tf_nm)
    non_TF_non_pos_cols <- non_pos_cols[non_pos_cols!="TF"]
    lrt_mat_comb <- rbind(lrt_sig_mat, lrt_nonsig_mat) %>% 
      ## transform all non_pos_cols values with signif
      # mutate(across(-all_of(non_TF_non_pos_cols), ~signif(., 9))) %>%
      lazy_dt()
    corr_mat_comb <- rbind(corr_sig_mat, corr_nonsig_mat) %>% 
      # mutate(across(-all_of(non_TF_non_pos_cols), ~signif(., 9))) %>%
      lazy_dt()
    wilc_mat_comb <- rbind(wilc_sig_mat, wilc_nonsig_mat) %>% 
      # mutate(across(-all_of(non_TF_non_pos_cols), ~signif(., 9))) %>%
      lazy_dt()
    
    fully_comb_mat <-
      unique(full_join(full_join(lrt_mat_comb, corr_mat_comb), wilc_mat_comb)) %>%
      mutate(
        lrt_pass = ifelse(!is.na(lrt_score), lrt_score < -10, F),
        corr_pass = ifelse(!is.na(phylop_IC_corr),
          phylop_IC_corr >= corr_thresh, F
        ),
        wilc_pass = ifelse(!is.na(left_p),
          left_p < 0.01 | right_p < 0.01, F
        )
      ) %>%
      as_tibble()

    motif_pos <- as.character(na.omit(as.numeric(colnames(fully_comb_mat))))
    motif_len <- length(motif_pos)
    total_sites <- nrow(fully_comb_mat)
    total_positions <- total_sites * motif_len
    fully_comb_mat_cons <- fully_comb_mat %>%
      # dplyr::select(-chr, -start, -end, -id, -pwm_score, -strand) %>%
      # rowwise() %>%
      mutate(
        n_positions = motif_len,
        n_conserved = rowSums(across(
          all_of(motif_pos),
          ~ . >= -log10(0.1)
        ), na.rm = T),
        n_positions_in_conserved_sites = ifelse(n_conserved > 0, n_positions, 0),
        n_positive = rowSums(across(
          all_of(motif_pos),
          ~ . <= log10(0.1)
        ), na.rm = T)
      )

    # write bed files of conserved and nonconserved sites
    fully_comb_mat_cons <- fully_comb_mat_cons %>%
      mutate(
        conserved_tfbs = (n_conserved >= 1) & (lrt_pass | corr_pass | wilc_pass)
      )
    fully_comb_mat_cons %>%
      write_csv(
        file.path(
          hmap_outs, paste0(tf_nm, "_NA"),
          paste0(tf_nm, "_", filt_set, "_all_sites.csv.gz")
        ) 
      )
      
    fully_comb_mat_cons %>%
      filter(conserved_tfbs) %>%
      dplyr::select(all_of(1:6)) %>%
      unique() %>%
      write_tsv(
        file.path(
          hmap_outs, paste0(tf_nm, "_NA"),
          paste0(tf_nm, "_", filt_set, "_conserved_sites.bed.gz")
        ),
        col_names = F
      )

    fully_comb_mat_cons %>%
      filter(!conserved_tfbs) %>%
      dplyr::select(all_of(1:6)) %>%
      write_tsv(
        file.path(
          hmap_outs, paste0(tf_nm, "_NA"),
          paste0(tf_nm, "_", filt_set, "_unconserved_sites.bed.gz")
        ),
        col_names = F
      )
  
    # import the full phyloP matrix
    full_mat <- 
      import_matrix(tf_nm,
                    grep(paste0(tf_nm, "_"), compmat_fils, value = T),
                    motif_positions = motif_positions,
                    filter_list = tfbs_subset_nms
    )
    write_csv(full_mat, file.path(
      hmap_outs, paste0(tf_nm, "_NA"),
      paste0(tf_nm, "_", filt_set, "_full_mat.csv.gz")
    ))
  } else{
    full_mat <- read_csv(
      file.path(
        hmap_outs, paste0(tf_nm, "_NA"),
        paste0(tf_nm, "_", filt_set, "_full_mat.csv.gz")
      )
    )
    fully_comb_mat_cons <- read_csv(
        file.path(
          hmap_outs, paste0(tf_nm, "_NA"),
          paste0(tf_nm, "_", filt_set, "_all_sites.csv.gz")
        ) 
      )
    motif_positions <- as.character(na.omit(as.numeric(colnames(fully_comb_mat_cons))))
    motif_pos <- motif_positions
  }
  # make heatmaps
  q_01 <- round(-log10(0.1), digits = 2)
  if (max(full_mat[as.character(1:100)], na.rm = T) >= 6) {
    breaks <- c(seq(-6, -q_01, 1), 0, seq(q_01, 6, 1))
    lg_breaks <- c(-6, -1.3, -q_01, 0, q_01, 1.3, 6)
  } else {
    breaks <- c(seq(-2, -q_01, 1), 0, seq(q_01, 2, 1))
    lg_breaks <- c(-2, -1.3, -q_01, 0, q_01, 1.3, 2)
  }

  full_mat_cons <- left_join(
    full_mat, fully_comb_mat_cons %>%
      select(
        all_of(non_pos_cols[1:6]),
        lrt_pass, corr_pass, wilc_pass, conserved_tfbs
      )
  )

  ## all sites - only if using triple opt/overlapped sets
  if (filt_set != "all_moods") {
    # import promoter data
    prom_intersect_df <-
      list.files(file.path(hmap_outs, "../../../phastcons_promoter_intersect/"),
        pattern = paste0(tf_nm, ".*_promoters_intersect.*.bed"),
        full.names = T
      ) %>%
      map(~ read_tsv(
        .x,
        col_names = c(
          "chr", "start", "end", "id", "score", "strand", "promoter_chr",
          "promoter_start", "promoter_end", "promoter_gene", "promoter_score", 
          "promoter_strand", "promoter_category",
          "promoter_intersect_length"
        ),
        col_types = c("chr" = "c", "promoter_score" = "c")
      ) %>%
        mutate(TF = gsub("_.*", "", basename(.x)))) %>%
      list_rbind() %>%
      ## filter to just binding intersects that are full length
      select(-promoter_score) %>%
      mutate(TF_tfbs_nm = paste(TF, chr, start, end, strand, sep = "_")) %>%
      filter((end - start) / promoter_intersect_length == 1 |
        promoter_intersect_length == 0) %>%
      group_by(TF_tfbs_nm, promoter_gene) %>%
      mutate(promoter_genes = paste(unique(promoter_gene), collapse = ";"))

    ## and create hmap
    row_annots <- full_mat_cons %>%
      unique() %>%
      left_join(prom_intersect_df %>%
        group_by(TF_tfbs_nm) %>%
        slice_max(promoter_intersect_length, with_ties = F) %>%
        mutate(promoter_intersect = promoter_intersect_length > 0)) %>%
      mutate(
        promoter_intersect =
          ifelse(is.na(promoter_intersect), FALSE, promoter_intersect),
        TF_tfbs_nm = paste(TF, chr, start, end, strand, sep = "_")
      ) %>%
      select(TF_tfbs_nm, lrt_pass, corr_pass, wilc_pass, conserved_tfbs, promoter_intersect) %>%
      mutate(across(c(lrt_pass, corr_pass, wilc_pass, promoter_intersect), 
             ~ifelse(is.na(.x), FALSE, .x)))
    ## add column for row_annots showing presence or absence of "main"
    ## (most enriched TE class containing the most enriched family)
    main_types <- read_csv(
      file.path(
        proj_dir,
        "../repeat_enrichment/analysis/optimized_peak_sets/CHS_Fig7_repeat_enrichment.csv.gz"
      )
    ) %>%
      filter(TF == TF_nm, grepl("LINE|LTR|SINE|DNA", Type)) %>%
      slice_max(`subtype -log10(p-val)`, with_ties = F) #
    main_type <<- main_types %>% pull(Subtype)
    main_type2 <<- main_types %>% pull(Type)
    if (length(main_type) > 0) {
      repeat_intersection <-
        read_tsv(
          file.path(
            proj_dir, "analysis/computematrix_out/trip_opt_moods_HC/",
            paste0(TF_nm, "_NA_moods_tripOpt_repeatmasker.bed")
          ),
          col_names = c(
            "tf_chr", "tf_start", "tf_end", "tf_nm", "tf_score",
            "tf_strand", "repeat_chr", "repeat_start",
            "repeat_end", "repeat_strand", "repeat_nm",
            "repeat_class", "overlap_length"
          )
        ) %>%
        mutate(TF_tfbs_nm = paste(TF_nm, tf_chr, tf_start, tf_end, tf_strand, sep = "_")) %>%
        filter(repeat_nm == main_type)
      row_annots <-
        left_join(
          row_annots,
          repeat_intersection %>% select(TF_tfbs_nm, repeat_nm)
        ) %>%
        unique() %>%
        mutate("{main_type2}-{main_type}" :=
          ifelse(is.na(repeat_nm), FALSE, TRUE)) %>%
        select(-repeat_nm, -TF_tfbs_nm) %>%
        dplyr::rename(L = lrt_pass, C = corr_pass, W = wilc_pass, P = promoter_intersect)
    } else {
      # repeat_intersection <- tibble(TF_tfbs_nm = F, repeat_nm = F)
      row_annots <-
        row_annots %>%
        mutate("No repeats bound" = F) %>%
        select(-TF_tfbs_nm) %>%
        dplyr::rename(L = lrt_pass, C = corr_pass, W = wilc_pass, P = promoter_intersect)
    }

    conserved_nm <- paste0(
      tf_nm,
      " n = ",
      nrow(full_mat_cons %>% filter(conserved_tfbs)), "/", nrow(full_mat_cons)
    )
    out_nm <- file.path(
      hmap_outs,
      paste0(tf_nm, "_NA"),
      paste0(tf_nm, "_", filt_set, "_cons+uncons_sites_hmap.pdf")
    )
    
    ## save table for dissemination
    if(length(main_type) > 0){
      TE_col_nm <- paste0("Intersects ", main_type2,"-", main_type)
    } else{TE_col_nm <- "No repeats bound"}
    annot_sub <- row_annots %>% select(all_of(5:6))
    colnames(annot_sub) <- c("Intersects promoter", TE_col_nm)
    full_mat_cons %>% select(-id, -pwm_score) %>% 
      dplyr::rename(`Passes LTR Test` = lrt_pass, `Passes Correlation Test` = corr_pass, 
             `Passes Wilcoxon Test` = wilc_pass, CTOP = conserved_tfbs) %>% 
      cbind(annot_sub) %>%
      write_csv(file.path(hmap_outs, paste0(tf_nm, "_NA"), 
                          paste0(tf_nm, "_", filt_set, "_full_mat_annots.csv.gz")))
    
    ## make heatmap
    produce_hmap(
      full_mat_cons %>%
        select(all_of(as.character(1:100))),
      mat_nm = conserved_nm, breaks = breaks,
      lg_breaks = lg_breaks,
      column_gap_pos = c(
        as.numeric(motif_positions[1]) - 1,
        as.numeric(motif_positions[length(motif_positions)])
      ),
      row_annots = row_annots, outfile = out_nm
    )

    ## conserved sites
    if (all_hmaps) {
      row_annots <- full_mat_cons %>%
        filter(conserved_tfbs) %>%
        select(lrt_pass, corr_pass, wilc_pass)
      conserved_nm <- paste0(
        "n = ",
        nrow(full_mat_cons[full_mat_cons$conserved_tfbs, ]), "/", nrow(full_mat_cons)
      )
      out_nm <- file.path(
        hmap_outs,
        paste0(tf_nm, "_NA"),
        paste0(tf_nm, "_", filt_set, "_conserved_sites_hmap.pdf")
      )
      if (nrow(full_mat_cons %>% filter(conserved_tfbs)) > 0) {
        produce_hmap(
          full_mat_cons %>% filter(conserved_tfbs) %>%
            select(all_of(as.character(1:100))),
          mat_nm = conserved_nm, breaks = breaks,
          lg_breaks = lg_breaks,
          column_gap_pos = c(
            as.numeric(motif_positions[1]) - 1,
            as.numeric(motif_positions[length(motif_positions)])
          ),
          row_annots = row_annots, outfile = out_nm
        )
      }
      # unconserved sites
      row_annots <- full_mat_cons %>%
        filter(!conserved_tfbs) %>%
        select(lrt_pass, corr_pass, wilc_pass)
      unconserved_nm <- paste0(
        "n = ",
        nrow(full_mat_cons[!full_mat_cons$conserved_tfbs, ]), "/", nrow(full_mat_cons)
      )
      out_nm <- file.path(
        hmap_outs,
        paste0(tf_nm, "_NA"),
        paste0(tf_nm, "_", filt_set, "_unconserved_sites_hmap.pdf")
      )
      if (nrow(full_mat_cons %>% filter(!conserved_tfbs)) > 0) {
        produce_hmap(
          full_mat_cons %>% filter(!conserved_tfbs) %>%
            select(all_of(as.character(1:100))),
          mat_nm = unconserved_nm, breaks = breaks,
          lg_breaks = lg_breaks,
          column_gap_pos = c(
            as.numeric(motif_positions[1]) - 1,
            as.numeric(motif_positions[length(motif_positions)])
          ),
          row_annots = row_annots, outfile = out_nm
        )
      }
    }
  }

  ## output summary table
  if(pipeline_section == "full"){
    conditional_mean <- function(nums_vec) {
      mean(nums_vec[nums_vec >= -log10(0.1)], na.rm = T)
    }
    fully_comb_mat_cons_stats <- fully_comb_mat_cons %>%
      # rowwise() %>%
      # mutate(max_phyloP =
      #          ifelse(any(!is.na(c_across(all_of(motif_pos)))),
      #                 max(c_across(all_of(motif_pos)), na.rm = T), NA),
      #        mean_cons_phyloP = conditional_mean(c_across(all_of(motif_pos)))) %>%
      ungroup() %>%
      mutate(
        TF = tf_nm,
        mean_phyloP = rowMeans(across(all_of(motif_pos)), na.rm = T),
        total_sites = total_sites,
        total_positions = total_positions,
        motif_len = motif_len
      ) %>%
      select(-all_of(motif_pos))
    fully_comb_mat_cons_stats <- fully_comb_mat_cons_stats %>%
      relocate(TF, .before = chr)
    fully_comb_mat_cons_stats %>%
      write_csv(file.path(
        hmap_outs, paste0(tf_nm, "_NA"),
        paste0(
          tf_nm, "_conservation_all_sites_",
          filt_set, ".csv.gz"
        )
      ))
  }
}
if (sys.nframe() == 0) {
  summarize_conservation(TF_nm)
}
