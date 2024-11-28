import_matrix <- function(TF_nm, mat_path, motif_positions, filter_list = NA) {
  require(dtplyr)
  print("run read_tsv")
  matrixout <- read_tsv(
    mat_path,
    col_names = c("chr", "start", "end", "id", "pwm_score", 
                  "strand", 1:100, "TF"),
    skip = 1, num_threads = 1
    ) #%>% lazy_dt()
  if (is.character(filter_list)) {
    print("filter down to predefined list of binding sites")
    matrixout <- matrixout %>%
      mutate(TF_tfbs_nm = paste0(TF, paste(chr, start, end, sep = "_"))) %>%
      filter(TF_tfbs_nm %in% filter_list) %>% select(-TF_tfbs_nm)
  }
  matrixout <- matrixout %>% lazy_dt()
  print("filter out motif hits that are NA")
  matrixout_mat <- matrixout %>%
    # select(as.character(1):as.character(100)) %>%
    filter(if_all(as.character(motif_positions), ~ !is.na(.)))
  coords <- matrixout_mat %>% select(1:6) %>% as.data.frame()
  
  print("make separate matrix for NA motif hits")
  gc()
  NA_motifs <- matrixout %>%
    # select(as.character(1):as.character(100)) %>%
    filter(if_any(as.character(motif_positions), ~ is.na(.))) %>% 
    select(-TF) %>% as.data.frame()

  # convert the PhyloP values to BH FDR PhyloPs
  print("convert the PhyloP values to BH FDR PhyloPs")
  qval_matrixout_mat <- 
    matrixout_mat %>% select(-all_of(non_pos_cols)) %>% 
    mutate(row_id = row_number()) %>% 
    pivot_longer(-row_id, names_to = "position", values_to = "phyloP") %>% 
    mutate(sign = sign(phyloP), 
           qval_phyloP = -log10(p.adjust(10^-abs(phyloP), "BH"))*sign) %>% 
    select(-phyloP, -sign) %>%
    pivot_wider(names_from = position, values_from = qval_phyloP) %>%
    arrange(row_id) %>% select(-row_id) %>% 
    as.data.frame() %>% relocate(as.character(1:100))
    
  # signs <- matrixout_mat %>% select(-all_of(non_pos_cols)) %>% as.data.frame() %>% sign() 
  # matrixout_mat_abs <- matrixout_mat %>% as.data.frame() %>% select(-all_of(non_pos_cols)) %>% abs()
  # rm(matrixout_mat)
  gc()
  # qval_matrixout_mat <- 
  #   -log10(p.adjust(10^-as.matrix(matrixout_mat_abs), "BH")) * signs
    # -log10(p.adjust(10^ - as.matrix(matrixout_mat), "BH")) *
    # signs
  # gc()
  as_tibble(rbind(cbind(coords, qval_matrixout_mat), NA_motifs))
}

# phyloP_IC_weighting <- function(
#     matrixout_mat, motif_positions, TF_nm,
#     alt, flank_length, norm = FALSE) {
#   IC_vals <- motif_IC_vals %>%
#     filter(TF == TF_nm) %>%
#     pull(ics) %>%
#     unlist()
#   IC_vals_norm <- IC_vals / 2 # weighted by max possible bit score
#   motif_size <- length(motif_positions)
#   matrixout_mat_pos <- matrixout_mat %>%
#     select(-all_of(non_pos_cols[1:6])) %>%
#     filter(if_all(everything(), ~ !is.na(.)))
#   motif_phylop_all <- matrixout_mat_pos %>%
#     select(as.character(motif_positions))
#   motif_phylop_all_weighted <-
#     sweep(motif_phylop_all, 2, IC_vals_norm, "*")
#   # stitch together the weighted motif and flanking regions
#   matrixout_mat_pos_weighted_motif <-
#     cbind(
#       matrixout_mat_pos[, 1:(motif_positions[1] - 1)],
#       motif_phylop_all_weighted,
#       matrixout_mat_pos[, (motif_positions[motif_size] + 1):ncol
#       (matrixout_mat_pos)]
#     )
#   # plot the stiched together matrix
#   q_01 <- round(-log10(0.1), digits = 2)
#   breaks <- c(seq(-2, -q_01, 1), seq(q_01, 2, 1))
#   pheatmap(
#     matrixout_mat_pos_weighted_motif %>%
#       slice_sample(n = 5000, replace = F),
#     cluster_cols = F,
#     # use_raster = F,
#     color = colorRampPalette(brewer.pal(
#       n = 7,
#       name = "RdBu"
#     ))(16),
#     breaks = breaks
#   )
# 
# 
#   # use this matrix as input to run_wilcox_ks
#   run_wilcox_ks(matrixout_mat_pos_weighted_motif, motif_positions,
#     alt = "two.sided", flank_length = flank_length, norm = FALSE, NA_mat
#   )
# }

run_wilcox_ks <- function(in_mat, motif_positions, alt, flank_length, norm = FALSE,
                          NA_mat, fig_out) {
  motif_length <- length(motif_positions)

  if (norm) {
    adj_motif_positions <- motif_positions - (motif_length / 2) + 1
    motif_mat <- in_mat[, adj_motif_positions]
    left_flank <- in_mat[, (adj_motif_positions[1] - flank_length):(adj_motif_positions[1] - 1)]
    right_flank <- in_mat[, (adj_motif_positions[motif_length] + 1):(adj_motif_positions[motif_length] + flank_length)]
  } else {
    motif_mat <- as.matrix(in_mat[, as.character(motif_positions)])
    left_flank <-
      as.matrix(in_mat[, as.character(
        (motif_positions[1] - flank_length):(motif_positions[1] - 1)
      )])
    right_flank <-
      as.matrix(in_mat[, as.character(
        (motif_positions[motif_length] + 1):(motif_positions[motif_length] +
          flank_length)
      )])
  }

  rowwise_wilcox <- function(i, motif_mat, left_flank, right_flank, alt) {
    if (all(is.na(left_flank[i, ]))) {
      left <- list(p.value = 1)
    } else{
      left <- wilcox.test(motif_mat[i, ], left_flank[i, ],
                          alternative = alt, exact = F)
    }
    if (all(is.na(right_flank[i, ]))) {
      right <- list(p.value = 1)
    } else{
      right <- wilcox.test(motif_mat[i, ], right_flank[i, ],
                           alternative = alt, exact = F)
    }
    tibble(site_n = i, left_p = left$p.value, right_p = right$p.value)
  }
  
  rowwise_wilcox2 <- function(motif_mat, left_flank, right_flank, alt) {
    left <- row_wilcoxon_twosample(motif_mat, left_flank,
                        alternative = alt, exact = F) %>%
      mutate(left_p = ifelse(obs.y < 25 | obs.x < motif_length, 1, pvalue))
    right <- row_wilcoxon_twosample(motif_mat, right_flank,
                         alternative = alt, exact = F) %>%
      mutate(right_p = ifelse(obs.y < 25 | obs.x < motif_length, 1, pvalue))
    tibble(left_p = left$left_p, right_p = right$right_p)
  }

  rowwise_ks <- function(i, motif_mat, left_flank, right_flank, alt) {
    left <- ks.test(motif_mat[i, ], left_flank[i, ],
      alternative = alt # , exact = F
    )
    right <- ks.test(motif_mat[i, ], right_flank[i, ],
      alternative = alt # , exact = F
    )
    tibble(site_n = i, left_p = left$p.value, right_p = right$p.value)
  }

  ## run wilcoxon/mann whitney test
  # plan(multicore, workers = 3)
  # wilcox_df <- map(
  #   1:nrow(motif_mat),
  #   ~ rowwise_wilcox(.x, motif_mat, left_flank, right_flank, alt)) %>% 
  #   list_rbind()
  # plan(sequential)
  wilcox_df <- rowwise_wilcox2(motif_mat, left_flank, right_flank, alt)
  p_thresh <- 0.01
  # passing_pos <- wilcox_df %>% 
  #   filter(left_p < p_thresh & right_p < p_thresh) %>%
  #   pull(site_n)
  # failing_pos <- wilcox_df %>% 
  #   filter(left_p >= p_thresh | right_p >= p_thresh) %>%
  #   pull(site_n)
  
  in_mat_wilcox <- cbind(in_mat, wilcox_df)
  in_mat_wilcox_pass <- in_mat_wilcox %>% 
    filter(left_p < p_thresh & right_p < p_thresh)
  in_mat_wilcox_fail <- in_mat_wilcox %>% 
    filter(left_p >= p_thresh | right_p >= p_thresh)
  
  # make heatmaps
  q_01 <- round(-log10(0.1), digits = 2)
  breaks <- c(seq(-6, -q_01, 1), 0, seq(q_01, 6, 1))
  lg_breaks <- c(-6, -1.3, -q_01, 0, q_01, 6, 2)
  max_length <- min(
    100 - length(motif_positions),
    motif_positions[length(motif_positions)]
  )
  
  # pass hmap
  pass_nm <- paste0(
    "n = ",
    nrow(in_mat_wilcox_pass), "/",
    (nrow(in_mat) + nrow(NA_mat))
  )
  pass_fil_out <- file.path(fig_out, paste0(
    TF_nm,
    "_wilcoxon_thresh_pass.pdf"
  ))
  produce_hmap(in_mat_wilcox_pass[, as.character(1:100)] , 
               pass_nm, breaks = breaks,
               lg_breaks = lg_breaks,
               column_gap_pos = c(motif_positions[1] - 1, max_length),
               row_annots = NA, outfile = pass_fil_out)
  
  # fail hmap
  fail_nm <- paste0(
    "n = ",
    nrow(in_mat_wilcox_fail), "/",
    (nrow(in_mat) + nrow(NA_mat))
  )
  fail_fil_out <- file.path(fig_out, paste0(
    TF_nm,
    "_wilcoxon_thresh_fail.pdf"
  ))
  produce_hmap(in_mat_wilcox_fail[, as.character(1:100)], fail_nm, 
               breaks = breaks, lg_breaks = lg_breaks,
               column_gap_pos = c(motif_positions[1] - 1, max_length),
               row_annots = NA, outfile = fail_fil_out)
  
  ## run KS test
  ### change alt if running KS test, since "less" and "greater" are reversed
  # if (alt == "greater" & norm == FALSE) {
  #   alt <- "less"
  # }
  # pval_vec <- c()
  # i <- 207
  # for(j in 18:40){
  #   ktest <- ks.test(motif_mat[i, 1:j], left_flank[i, 1:25],
  #                    alternative = alt#, exact = F
  #   )
  #   pval_vec <- c(pval_vec, ktest$p.value)
  #   }

  # ks_df <- map(
  #   1:nrow(motif_mat),
  #   ~ rowwise_ks(.x, motif_mat, left_flank, right_flank, alt)
  # ) %>%
  #   list_rbind()
  
  print(file.path(fig_out, paste0(TF_nm, "_wilcoxon_thresh_pass.tsv.gz")))
  write_tsv(in_mat_wilcox_pass %>% select(
                   all_of(c("chr", "start", "end", "id", "pwm_score", "strand", 
                     "left_p", "right_p", as.character(motif_positions)))), 
            file.path(fig_out, paste0(TF_nm, "_wilcoxon_thresh_pass.tsv.gz")),
            num_threads = 1
            )
  write_tsv(rbind(in_mat_wilcox_fail %>% select(
    all_of(c("chr", "start", "end", "id", "pwm_score", "strand", 
             "left_p", "right_p", as.character(motif_positions)))), 
    NA_mat[,c("chr", "start", "end", "id", "pwm_score", "strand",
              as.character(motif_positions))] %>%
      mutate(left_p = NA, right_p = NA)), 
    file.path(fig_out, paste0(TF_nm, "_wilcoxon_thresh_fail.tsv.gz")),
    num_threads = 1
    )
}

plot_filtered_hmaps <- function(in_mat, wilcox_res, type = "nonnorm",
                                out_dir, stat, TF_nm, motif_start,
                                motif_end) {
  wilcox_res_fdr <- wilcox_res %>%
    mutate(
      left_p_FDR = p.adjust(left_p, method = "BH"),
      right_p_FDR = p.adjust(right_p, method = "BH")
    )
  sig_thresh <- 0.1
  wilcox_res_sig <- wilcox_res_fdr %>%
    filter(left_p_FDR < sig_thresh & right_p_FDR < sig_thresh)
  wilcox_res_nonsig <- wilcox_res_fdr %>%
    filter(left_p_FDR > sig_thresh | right_p_FDR > sig_thresh)

  ## heatmap setup
  q_01 <- round(-log10(0.1), digits = 2)
  # breaks <- c(seq(-10, -q_01, 1), seq(q_01, 10, 1))
  breaks <- c(seq(-2, -q_01, 1), seq(q_01, 2, 1))
  run_pheatmap <- function(mat) {
    ## this runs ultra-slow, so speeding it up by sampling to 2000 rows
    total_rows <- nrow(mat)
    max_rows <- 5000
    if (nrow(mat) > max_rows) {
      mat <- mat[sample(1:nrow(mat), max_rows), ]
    }
    pheatmap(mat,
      cluster_rows = TRUE, cluster_cols = FALSE,
      name = paste0("n = ", nrow(mat), "/", total_rows),
      color =
        colorRampPalette(brewer.pal(
          n = 7,
          name = "RdBu"
        ))(20),
      breaks = breaks,
      # legend_breaks = c(-10, -1.3, -q_01, q_01, 1.3, 10),
      legend_breaks = c(-2, -1.3, -q_01, q_01, 1.3, 2),
      gaps_col = c(motif_start - 1, motif_end),
      use_raster = F,
      fontsize_col = 6
    )
  }
  sig_hmap <- run_pheatmap(in_mat[wilcox_res_sig$site_n, 7:ncol(in_mat)])
  non_sighmap <- run_pheatmap(in_mat[wilcox_res_nonsig$site_n, 7:ncol(in_mat)])

  # pdf(file.path(out_dir, paste0(TF_nm, "merged_conservation_wilcoxon_hmap_", type, ".pdf")))
  # wrap_elements(grid.grabExpr(draw(sig_hmap))) +
  #   wrap_elements(grid.grabExpr(draw(non_sighmap)))
  # dev.off()

  pdf(file.path(out_dir, paste0(
    TF_nm, "_sig_conservation_",
    stat, "_hmap_", type, ".pdf"
  )))
  draw(sig_hmap)
  dev.off()

  pdf(file.path(out_dir, paste0(
    TF_nm, "_nonsig_conservation_",
    stat, "_hmap_", type, ".pdf"
  )))
  draw(non_sighmap)
  dev.off()

  gc()
  # write the sig motif sites to a file
  write_tsv(
    in_mat[
      wilcox_res_sig$site_n,
      c(non_pos_cols[1:6], as.character(motif_start:motif_end))
    ],
    file.path(out_dir, paste0(
      TF_nm, "_sig_motifs_",
      stat, "_", type, ".txt.gz"
    )),
    num_threads = 1
  )
  # write the nonsig motif sites to a file
  write_tsv(
    in_mat[
      wilcox_res_nonsig$site_n,
      c(non_pos_cols[1:6], as.character(motif_start:motif_end))
    ],
    file.path(out_dir, paste0(
      TF_nm, "_nonsig_motifs_",
      stat, "_", type, ".txt.gz"
    )),
    num_threads = 1
  )
}

get_conserved_sites <- function(TF_nm, mat_fil) {
  TF_nm <<- TF_nm
  mat_fil <<- unique(mat_fil)
  non_pos_cols <<- c("chr", "start", "end", "id", "pwm_score", "strand", "TF")

  # get motif size and positions
  motif_start_pos <- motif_sizes %>%
    filter(TF == TF_nm) %>%
    pull(motif_start_pos)
  motif_end_pos <- (motif_sizes %>% filter(TF == TF_nm) %>% 
                      pull(motif_end_pos))
  motif_positions <- as.character(motif_start_pos):as.character(motif_end_pos)
  
  # import data - NOTE THAT THIS CONVERTS PHYLOP VALUES TO BH FDR PhyloPs
  print("import matrix")
  matrixout_mat <- import_matrix(TF_nm, mat_fil, motif_positions)
  NA_mat <- matrixout_mat %>%
    filter(if_any(as.character(motif_positions), ~is.na(.)))
  matrixout_mat <- matrixout_mat %>%
    filter(if_all(as.character(motif_positions), ~!is.na(.)))
  
  ## run mann-whitney and Kolmogorov-Smirnov test on each row, comparing motif
  ## phylop to each flank phyloP separately
  ### non-normalized
  print("run wilcox tests")
  fig_out <- unique(file.path(dirname(mat_fil), "hmaps", paste0(TF_nm, "_NA")))
  wilcox_ks_res <- run_wilcox_ks(matrixout_mat, motif_positions,
    "greater",
    flank_length = 25,
    norm = FALSE,
    NA_mat, fig_out
  )
  # 
  # plot_filtered_hmaps(matrixout_mat, wilcox_ks_res[[1]],
  #   out_dir = fig_out, stat = "wilcoxon",
  #   TF_nm = TF_nm,
  #   motif_start = motif_start_pos,
  #   motif_end = motif_end_pos
  # )
  # plot_filtered_hmaps(matrixout_mat, wilcox_ks_res[[2]],
  #   out_dir = fig_out, stat = "ks",
  #   TF_nm = TF_nm,
  #   motif_start = motif_start_pos,
  #   motif_end = motif_end_pos
  # )

  # testing different methods of separating conserved vs non-conserved
  # Ali's lrt method - calculation
  print("calculate LRT values")
  lrt_mat_list <- lrt_mat(matrixout_mat, motif_positions, fig_out)
  lrt_mat_phylop_df <- lrt_mat_list[[1]]
  lrt_mat_IC_df <- lrt_mat_list[[2]]
  matrixout_mat_pos <- lrt_mat_list[[3]]
  ## analyze lrt results, divide into binding site sets
  # phylop
  # analyze_lrt_mats(lrt_mat_phylop_df, "phylopMean", motif_positions, fig_out)
  # # IC
  # if lrt_mat_IC_df is null because of no conservation, create a null df
  if (is.null(lrt_mat_IC_df)) {
    lrt_mat_IC_df <- as_tibble(matrix(nrow = nrow(matrixout_mat), 
                                       ncol = 100 - length(motif_positions)))
    lrt_mat_IC_df[1:nrow(matrixout_mat), 
                  1:(100 - length(motif_positions))] <- 1e100
  }
  
  print("analyze LRT values")
  analyze_lrt_mats(lrt_mat_IC_df, matrixout_mat, 
                   matrixout_mat_pos, "icVals", 
                   TF_nm, motif_positions, fig_out, NA_mat)
  
  # using correlation + permutations
  print("calculate IC correlation values")
  corr_perm_tfbs_conservation(
    matrixout_mat,
    motif_positions, fig_out, lrt_mat_IC_df, NA_mat
  )

  # Tim's IC weighting scheme
  # wilcox_ks_res_weighted_motif <-
  #   phyloP_IC_weighting(matrixout_mat, motif_positions, TF_nm,
  #                       alt = "greater", flank_length = 25,
  #                       norm = FALSE)
  # plot_filtered_hmaps(matrixout_mat, wilcox_ks_res_weighted_motif[[1]],
  #                     out_dir = fig_out, stat = "wilcoxon_IC_weighted",
  #                     TF_nm = TF_nm,
  #                     motif_start = motif_start_pos,
  #                     motif_end = motif_end_pos)
  # plot_filtered_hmaps(matrixout_mat, wilcox_ks_res_weighted_motif[[2]],
  #                     out_dir = fig_out, stat = "ks_IC_weighted",
  #                     TF_nm = TF_nm,
  #                     motif_start = motif_start_pos,
  #                     motif_end = motif_end_pos)

  #   bc_res <- microbenchmark::microbenchmark(
  #     "import" = {
  #       matrixout_mat <- import_matrix(TF_nm, mat_fil)},
  #     "get motif sizes" = {
  #       # get motif size and positions
  #       motif_start_pos <- motif_sizes %>%
  #         filter(TF == TF_nm) %>%
  #         pull(motif_start_pos)
  #       motif_end_pos <- (motif_sizes %>% filter(TF == TF_nm) %>% pull(motif_end_pos))
  #       motif_positions <- as.character(motif_start_pos):as.character(motif_end_pos)
  #     },
  #     "run lrt" = {
  #       matrixout_mat_lrt <- lrt_mat(matrixout_mat, motif_positions, fig_out)
  #     },
  #
  #     "run MW and KS" = {
  #       ## run mann-whitney and Kolmogorov-Smirnov test on each row, comparing motif
  #       ## phylop to each flank phyloP separately
  #       ### non-normalized
  #       fig_out <- file.path(dirname(mat_fil), "hmaps", paste0(TF_nm, "_NA"))
  #       wilcox_ks_res <- run_wilcox_ks(matrixout_mat, motif_positions,
  #                                      "greater", flank_length = 25,
  #                                      norm = FALSE)
  #     },
  #     "make plots" = {
  #       plot_filtered_hmaps(matrixout_mat, wilcox_ks_res[[1]],
  #                           out_dir = fig_out, stat = "wilcoxon",
  #                           TF_nm = TF_nm,
  #                           motif_start = motif_start_pos,
  #                           motif_end = motif_end_pos)
  #       plot_filtered_hmaps(matrixout_mat, wilcox_ks_res[[2]],
  #                           out_dir = fig_out, stat = "ks",
  #                           TF_nm = TF_nm,
  #                           motif_start = motif_start_pos,
  #                           motif_end = motif_end_pos)
  #     }, times = 2)
  #   print(bc_res)
  #   p <- autoplot(bc_res)
  #   ggsave(file.path(fig_out, "benchmarking.pdf"), p)
  ### normalized
  # wilcox_ks_res_norm <- run_wilcox(matrixout_mat_norm,
  #                               motif_positions, "less", 10, TRUE)
  # plot_filtered_hmaps(matrixout_mat, wilcox_res_norm[[1]], "norm", out_dir = fig_out)
}

print(Sys.info()["nodename"])
proj_dir <-
  ifelse(grepl("MacBook", Sys.info()["nodename"]),
    "~/iyellan_bc2-clust/codebook_proj/tfbs_conservation",
    ifelse(grepl("gra", Sys.info()["nodename"]),
           "/scratch/iyellan/codebook_proj/tfbs_conservation",
           "/home/hugheslab1/iyellan/codebook_proj/tfbs_conservation"
  ))

# only runs if run as Rscript, not if sourced
if (sys.nframe() == 0) {
  # require(tidyverse, quietly = T)
  require(dplyr, quietly = T)
  require(readr, quietly = T)
  require(tibble, quietly = T)
  require(tidyr, quietly = T)
  require(purrr, quietly = T)
  require(ggplot2, quietly = T)
  require(stringr, quietly = T)
  require(data.table, quietly = T)
  require(dtplyr)
  require(furrr, quietly = T)
  require(ComplexHeatmap, quietly = T)
  require(patchwork, quietly = T)
  require(grid, quietly = T)
  require(RColorBrewer, quietly = T)
  require(parallel, quietly = T)
  # require(smoother, quietly = T)
  require(pracma, quietly = T)
  require(matrixTests, quietly = T)
  source(file.path(proj_dir, "scripts/conservation_analysis_lrt_fns.R"))
  options(future.globals.maxSize = 9000*1024^2)
  # memory.limit(999999999)
  gc()
  
  # three different sets of computematrix paths: the original from the 200 
  # something mammals, then the 444 something with the new primates, and the 
  # just primates one
  # computematrix_mat_path <- file.path(proj_dir, "analysis/computematrix_out")
  # computematrix_mat_path <- file.path(proj_dir, "analysis/computematrix_out/trip_opt_moods_HC") ## Ali's new affinity-filtered triple optimized ones

  if (is.na(commandArgs(trailingOnly = TRUE)[1])) {
    if (!grepl("MacBook", Sys.info()["nodename"])) {
      # plan(tweak(multicore, workers = 6))
      # plan(multicore)
    }
    # options(future.globals.onReference = "error")
    TF_nms <- gsub(".*PWM[0-9]+_(.*)_NA.*Matrixout_TF.*", "\\1", compmat_fils)
  }

  motif_IC_vals <<- readRDS(file.path(proj_dir, "TE_conservation/IC_correlation/codebook_motif_positionwise_ICs.rds"))
  motif_sizes <<- read_tsv(file.path(proj_dir, "analysis/motif_sizes.txt"), 
                           num_threads = 1
                           ) %>%
    mutate(
      motif_start_pos = 50 - ceiling(motif_len / 2) + 1, # fixed coordinates
      motif_end_pos = 50 + floor(motif_len / 2)
    )
  ic_phylop_corrs <<- read_csv(
    file.path(
      proj_dir,
      "TE_conservation/IC_correlation/motif_IC_correlations_trip_opt_moods_HC.csv"),
    num_threads = 1
    )
  if (is.na(commandArgs(trailingOnly = TRUE)[1])) {
    future_walk2(TF_nms, compmat_fils, \(x, y) get_conserved_sites(x, y))
    # get_conserved_sites("CTCF", grep("CTCF", compmat_fils, value = T))
    # # get_conserved_sites("ZBTB40", grep("ZBTB40", compmat_fils, value = T))
    # # # # get_conserved_sites("ZNF721", grep("ZNF721", compmat_fils, value = T))
    # get_conserved_sites("SP140", grep("SP140_", compmat_fils, value = T))
    # get_conserved_sites("LEUTX", grep("LEUTX_", compmat_fils, value = T))
    get_conserved_sites("CGGBP1", grep("CGGBP1", compmat_fils, value = T))
    # # # get_conserved_sites("ZNF623", grep("ZNF623", compmat_fils, value = T))
    # get_conserved_sites("ZNF696", grep("ZNF696", compmat_fils, value = T))
    # get_conserved_sites("SALL3", grep("SALL3", compmat_fils, value = T))
    # get_conserved_sites("VDR", grep("VDR", compmat_fils, value = T))
    # get_conserved_sites("LEF1", grep("LEF1", compmat_fils, value = T))
    # get_conserved_sites("ZBED2", grep("ZBED2", compmat_fils, value = T))
    # get_conserved_sites("ZNF233", grep("ZNF233", compmat_fils, value = T))
    # get_conserved_sites("CXXC4", grep("CXXC4", compmat_fils, value = T))    #
    # get_conserved_sites("ZNF251", grep("ZNF251", compmat_fils, value = T))
    get_conserved_sites("ZNF384", grep("ZNF384_", compmat_fils, value = T))
    # get_conserved_sites("ZNF57", grep("ZNF57_", compmat_fils, value = T))
    # get_conserved_sites("ZNF66", grep("ZNF66_", compmat_fils, value = T))
    # get_conserved_sites("ZSCAN2", grep("ZSCAN2_", compmat_fils, value = T))
    
    get_conserved_sites("MKX", grep("MKX_", compmat_fils, value = T))
  } else {
    TF_nm <- commandArgs(trailingOnly = TRUE)[1]
    computematrix_mat_path <- commandArgs(trailingOnly = T)[2]
    compmat_fils <-
      list.files(computematrix_mat_path,
                 pattern = "PWM.*_50_Matrixout_TF",
                 full.names = TRUE
      )
    print(TF_nm)
    get_conserved_sites(unique(TF_nm), 
                        unique(grep(paste0(TF_nm, "_"), compmat_fils, value = T))
                        )
  }
}
