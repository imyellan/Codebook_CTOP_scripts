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
hmap_outs <- file.path(proj_dir, "analysis/computematrix_out/all_PWM_hits/hmaps")
# filt_set <- "trip_opt_moods_HC"
filt_set <- commandArgs(trailingOnly = T)[1]
confidence <- ifelse(grepl("_HC", filt_set), "_HC", "")

## table to differentiate controls from Codebook TFs
codebook_metadat <- read_csv(
  file.path(proj_dir, "../data/all_codebook_TFs_ages_DBDs_darks.csv")
)
codebook_metadat <- codebook_metadat %>%
  mutate(
    DBD_simple = ifelse(grepl("C2H2", DBD_new), DBD_new, "Other"),
    TF = `Gene name`, `Role in study`, `TF category`, DBD_simple, max_age
  )

# import conservation stats
stat_pass_summ_cons_df <- read_csv(
  list.files(hmap_outs,
    pattern = paste0("_conservation_all_sites_", filt_set, ".csv.gz"),
    full.names = T, recursive = T
  )
) %>%
  lazy_dt() %>%
  mutate(
    tfbs_nm = paste(chr, start, end, strand, sep = "_"),
    TF_tfbs_nm = paste(TF, tfbs_nm, sep = "_")
  )

# import bedtools clustering results
clustering_results_import <- function(in_path, clust_dist) {
  clust_pref <- paste0("cluster_", clust_dist, "_")
  read_tsv(
    in_path,
    col_names =
      paste0(clust_pref, c("chr", "start", "end", "strands", "tfbs_nm", "id"))
  ) %>%
    rowwise() %>%
    mutate(
      "{clust_pref}TFs" := paste(
        unique(gsub(
          "_.*", "",
          unlist(strsplit(.data[[paste0(clust_pref, "tfbs_nm")]], ";"))
        )),
        collapse = ";"
      )
    ) %>%
    ungroup() %>%
    mutate(
      "{clust_pref}length" := .data[[paste0(clust_pref, "end")]] -
        .data[[paste0(clust_pref, "start")]],
      TF_tfbs_nm = strsplit(.data[[paste0(clust_pref, "tfbs_nm")]], ";")
    ) %>%
    unnest(c(TF_tfbs_nm)) %>%
    lazy_dt() %>%
    mutate(
      "unambiguous_tfbs_{clust_dist}" := !grepl(";", .data[[paste0(clust_pref, "TFs")]])
    ) %>%
    unique()
}

## 100bp clusters
clustering_filt_set <- "trip_opt_moods_HC"
bedtools_merge_df_100 <- clustering_results_import(file.path(
  proj_dir,
  paste0(
    "analysis/computematrix_out/all_PWM_hits/hmaps/all_",
    clustering_filt_set, "_conserved_clustered_100.bed"
  )
), 100)

# 1000 bp clusters
bedtools_merge_df_1000 <- clustering_results_import(file.path(
  proj_dir,
  paste0(
    "analysis/computematrix_out/all_PWM_hits/hmaps/all_",
    clustering_filt_set, "_conserved_clustered_1000.bed"
  )
), 1000)



# clusters with 0bp distance
bedtools_merge_df_0 <- clustering_results_import(file.path(
  proj_dir,
  paste0(
    "analysis/computematrix_out/all_PWM_hits/hmaps/all_",
    clustering_filt_set, "_conserved_clustered.bed"
  )
), 0)

## completely overlapped tfbs's
complete_overlaps_df <- read_tsv(
  file.path(
    proj_dir,
    paste0(
      "analysis/computematrix_out/all_PWM_hits/hmaps/all_", clustering_filt_set,
      "_conserved_complete_overlaps.bed"
    )
  ),
  col_names = c(
    "site1_chr", "site1_start", "site1_end",
    "site1_id", "site1_score", "site1_strand", "site1_tfbs_nm",
    "site2_chr", "site2_start", "site2_end",
    "site2_id", "site2_score", "site2_strand", "site2_tfbs_nm",
    "overlap_length"
  )
) %>%
  rowwise() %>%
  mutate(overlap_id = paste(sort(c(site1_tfbs_nm, site2_tfbs_nm)), collapse = ";"))


# join conservation stats to clustering results
stat_pass_summ_cons_df_clusts <- stat_pass_summ_cons_df %>%
  left_join(bedtools_merge_df_0) %>%
  left_join(bedtools_merge_df_100) %>%
  left_join(bedtools_merge_df_1000)

# incorporate intersection with protein coding gene promoters and PhastCons elements
promoter_intersect_df <-
  list.files(file.path(proj_dir, "analysis/phastcons_promoter_intersect"),
    pattern = "promoters_intersect.bed", full.names = T
  ) %>%
  map(~ read_tsv(.x,
    col_names = c(
      "chr", "start", "end", "id", "score", "strand",
      "promoter_chr", "promoter_start", "promoter_end",
      "promoter_gene", "promoter_score", "promoter_strand",
      "promoter_category",
      "promoter_intersect_length"
    ), 
    col_types = "c"
  ) %>%
    mutate(TF = gsub("_.*", "", basename(.x))) %>% 
    select(-promoter_score)
  ) %>%
  list_rbind() %>%
  lazy_dt() %>%
  ## filter to just binding intersects that are full length
  mutate(TF_tfbs_nm = paste(TF, chr, start, end, strand, sep = "_")) %>%
  filter((end - start) / promoter_intersect_length == 1 |
    promoter_intersect_length == 0) %>%
  group_by(TF_tfbs_nm, promoter_gene) %>%
  mutate(promoter_genes = paste(unique(promoter_gene), collapse = ";"))

enhancer_intersect_df <-
  list.files(file.path(proj_dir, "analysis/phastcons_promoter_intersect"),
    pattern = "enhancer_intersect.bed", full.names = T
  ) %>%
  map(~ read_tsv(.x,
    col_names = c(
      "chr", "start", "end", "id", "score", "strand",
      "enhancer_chr", "ehancer_start", "enhancer_end",
      "enhancer_score1", "enhancer_score2",
      "enhancer_score3", "enhancer_intersect_length"
    ),
    col_types = "c"
  ) %>%
    mutate(TF = gsub("_.*", "", basename(.x)))) %>%
  list_rbind() %>%
  lazy_dt() %>%
  mutate(TF_tfbs_nm = paste(TF, chr, start, end, strand, sep = "_"))
phastcons_intersect_df <-
  list.files(file.path(proj_dir, "analysis/phastcons_promoter_intersect"),
    pattern = "phastcons_intersect.bed", full.names = T
  ) %>%
  map(~ read_tsv(.x,
    col_names = c(
      "chr", "start", "end", "id", "score", "strand",
      "phastcons_chr", "phastcons_start", "phastcons_end",
      "phastcons_lod", "phastcons_score",
      "phastcons_intersect_length"
    ),
    col_types = "c"
  ) %>%
    mutate(TF = gsub("_.*", "", basename(.x)))) %>%
  list_rbind() %>%
  lazy_dt() %>%
  mutate(TF_tfbs_nm = paste(TF, chr, start, end, strand, sep = "_"))

# intersection with masked CpG islands
masked_cpg_intersect_df <-
  list.files(file.path(proj_dir, "analysis/phastcons_promoter_intersect"),
    pattern = "cpgs_masked_intersect.bed", full.names = T
  ) %>%
  map(~ read_tsv(.x,
    col_names = c(
      "chr", "start", "end", "id", "score", "strand",
      "masked_cpg_chr", "masked_cpg_start", "masked_cpg_end",
      "masked_cpg_nm", "masked_cpg_gene", "masked_cpg_category",
      "masked_cpg_intersect_length"
    ),
    col_types = "c"
  ) %>%
    mutate(TF = gsub("_.*", "", basename(.x)))) %>%
  list_rbind() %>%
  lazy_dt() %>%
  mutate(
    masked_cpg_nm = paste0(masked_cpg_chr, masked_cpg_start, masked_cpg_end),
    TF_tfbs_nm = paste(TF, chr, start, end, strand, sep = "_")
  )

# intersection with unmasked CpG islands
unmasked_cpg_intersect_df <-
  list.files(file.path(proj_dir, "analysis/phastcons_promoter_intersect"),
    pattern = "cpgs_unmasked_intersect.bed", full.names = T
  ) %>%
  map(~ read_tsv(.x,
    col_names = c(
      "chr", "start", "end", "id", "score", "strand",
      "unmasked_cpg_chr", "unmasked_cpg_start", "unmasked_cpg_end",
      "unmasked_cpg_nm", "unmasked_cpg_gene", "unmasked_cpg_category",
      "unmasked_cpg_intersect_length"
    ),
    col_types = "c"
  ) %>%
    mutate(TF = gsub("_.*", "", basename(.x)))) %>%
  list_rbind() %>%
  lazy_dt() %>%
  mutate(
    unmasked_cpg_nm = paste0(unmasked_cpg_chr, unmasked_cpg_start, unmasked_cpg_end),
    TF_tfbs_nm = paste(TF, chr, start, end, strand, sep = "_")
  )

cg_content_df <- read_tsv(
  file.path(
    proj_dir,
    paste0(
      "analysis/phastcons_promoter_intersect/conservation_all_sites_trip_opt_moods_HC_cg_mean.bed"
    )
  ),
  col_names = c("chr", "start", "end", "TF_tfbs_nm", "CG_perc_mean")
) %>% lazy_dt()
# cluster_cg_content_df <-
#   read_tsv(file.path(proj_dir,
#                      paste0("analysis/phastcons_promoter_intersect/all_conserved_clustered_cg_mean.bed")),
#     col_names = c(
#       "cluster_chr", "cluster_start", "cluster_end",
#       "cluster_tfbs_nm", "cluster_CG_perc_mean"
#     )
#   ) %>% lazy_dt()

# intersection with repeats
repeat_intersect_df <-
  list.files(file.path(proj_dir, "analysis/phastcons_promoter_intersect"),
    pattern = "repeats_intersect(_unconserved)?.bed", full.names = T
  ) %>%
  map(~ read_tsv(.x,
    col_names = c(
      "chr", "start", "end", "id", "score", "strand",
      "repeat_chr", "repeat_start", "repeat_end",
      "repeat_family", "repeat_score", "repeat_strand", 
      "repeat_class", "repeat_id", "repeat_intersect_length"
    ), col_types = c("repeat_id" = "c", "repeat_score" = "c")
  ) %>%
    mutate(TF = gsub("_.*", "", basename(.x)))) %>%
  list_rbind() %>%
  lazy_dt() %>%
  mutate(TF_tfbs_nm = paste(TF, chr, start, end, strand, sep = "_"))

stat_pass_summ_cons_df_clusts_prom_phast <- stat_pass_summ_cons_df_clusts %>%
  left_join(promoter_intersect_df %>%
    select(TF_tfbs_nm, promoter_intersect_length, promoter_genes, 
           promoter_category) %>%
    unique()) %>%
  left_join(enhancer_intersect_df %>%
    select(TF_tfbs_nm, enhancer_intersect_length) %>%
    group_by(TF_tfbs_nm) %>%
    unique() %>%
    slice_max(enhancer_intersect_length, with_ties = F)) %>%
  left_join(phastcons_intersect_df %>%
    select(TF_tfbs_nm, phastcons_intersect_length) %>%
    group_by(TF_tfbs_nm) %>%
    unique() %>%
    slice_max(phastcons_intersect_length, with_ties = F)) %>%
  left_join(masked_cpg_intersect_df %>%
    select(
      TF_tfbs_nm, masked_cpg_nm, masked_cpg_gene, masked_cpg_category,
      masked_cpg_intersect_length
    ) %>%
    group_by(TF_tfbs_nm, masked_cpg_nm) %>%
    unique() %>%
    slice_max(masked_cpg_intersect_length, with_ties = F)) %>%
  left_join(unmasked_cpg_intersect_df %>%
    select(
      TF_tfbs_nm, unmasked_cpg_nm, unmasked_cpg_gene, unmasked_cpg_category,
      unmasked_cpg_intersect_length
    ) %>%
    group_by(TF_tfbs_nm, unmasked_cpg_nm) %>%
    unique() %>%
    slice_max(unmasked_cpg_intersect_length, with_ties = F)) %>%
  left_join(cg_content_df %>%
    select(TF_tfbs_nm, CG_perc_mean) %>%
    group_by(TF_tfbs_nm) %>% unique() %>%
    slice_max(CG_perc_mean, with_ties = F)) %>%
  # left_join(cluster_cg_content_df %>%
  #   select(cluster_tfbs_nm, cluster_CG_perc_mean) %>%
  #   group_by(cluster_tfbs_nm) %>% unique() %>%
  #   slice_max(cluster_CG_perc_mean, with_ties = F)) %>%
  left_join(repeat_intersect_df %>%
    select(TF_tfbs_nm, repeat_family, repeat_class, repeat_intersect_length) %>%
    group_by(TF_tfbs_nm) %>%
    unique() %>%
    slice_max(repeat_intersect_length, with_ties = F))

# write merged results
stat_pass_summ_cons_df_clusts_prom_phast <-
  stat_pass_summ_cons_df_clusts_prom_phast %>%
  left_join(codebook_metadat %>% select(TF, `Role in study`, `TF category`) %>%
    lazy_dt()) %>%
  as_tibble()
write_csv(
  stat_pass_summ_cons_df_clusts_prom_phast,
  file.path(hmap_outs, paste0("conservation_all_sites_", filt_set, ".csv.gz"))
)

# create different table using the conserved + unconserved clustering results
bedtools_merge_all_df_200 <- clustering_results_import(file.path(
  proj_dir,
  paste0(
    "analysis/computematrix_out/all_PWM_hits/hmaps/all_",
    clustering_filt_set, "_conserved+unconserved_clustered_200.bed"
  )
), 200)
bedtools_merge_all_df_1000 <- clustering_results_import(file.path(
  proj_dir,
  paste0(
    "analysis/computematrix_out/all_PWM_hits/hmaps/all_",
    clustering_filt_set, "_conserved+unconserved_clustered_1000.bed"
  )
), 1000)
stat_pass_summ_cons_df_uncons_cons_clusts_prom_phast <-
  stat_pass_summ_cons_df_clusts_prom_phast %>%
  lazy_dt() %>%
  select(-contains("cluster")) %>%
  left_join(bedtools_merge_all_df_200) %>%
  left_join(bedtools_merge_all_df_1000) %>%
  as_tibble()
write_csv(
  stat_pass_summ_cons_df_uncons_cons_clusts_prom_phast,
  file.path(hmap_outs, paste0(
    "conservation_all_sites_", filt_set,
    "_full_clustering.csv.gz"
  ))
)

# create table summarizing conservation results -
# inclusive of ambiguous sites and exclusive
for (ambig in c("ambig+unambig", "unambiguous")) {
  if (ambig == "unambiguous") {
    df <- stat_pass_summ_cons_df_clusts_prom_phast %>%
      filter(unambiguous_tfbs_0)
  } else {
    df <- stat_pass_summ_cons_df_clusts_prom_phast
  }
  sig_mat_cons_sum <- df %>%
    group_by(TF) %>%
    mutate(
      n_conserved_cons_tfbs = ifelse(conserved_tfbs, n_conserved, 0),
      n_positions_in_conserved_sites_cons_tfbs =
        ifelse(n_positions_in_conserved_sites, n_conserved, 0)
    ) %>%
    summarize(
      total_sites = n_distinct(TF_tfbs_nm),
      total_positions = sum(n_positions),
      motif_len = unique(motif_len),
      n_conserved_pos = sum(n_conserved_cons_tfbs),
      n_conserved_sites =
        as.numeric(length(which(conserved_tfbs))),
      n_positions_in_conserved_sites =
        sum(n_positions_in_conserved_sites_cons_tfbs),
      n_positive_pos = sum(n_positive)
    )

  sig_mat_cons_sum_metadat <- sig_mat_cons_sum %>%
    mutate(
      prop_conserved_pos = n_conserved_pos / total_positions,
      prop_conserved_sites = n_conserved_sites / total_sites
    ) %>%
    left_join(codebook_metadat %>%
      transmute(TF = `Gene name`, `Role in study`, DBD_simple, max_age))

  # write to file
  write_csv(
    sig_mat_cons_sum_metadat,
    file.path(
      hmap_outs,
      paste0("conservation_summary_", filt_set, "_", ambig, ".csv")
    )
  )
}

for (col_aes in c("DBD_simple", "max_age")) {
  if (col_aes == "DBD_simple") {
    legend_fn <- guide_legend(title = "TF type")
    file_suffix <- paste0("TF_type", confidence)
  } else {
    legend_fn <- guide_colourbar(title = "# of conserved positions", barheight = 7)
    file_suffix <- paste0("TF_age", confidence)
  }
  p <- sig_mat_cons_sum_metadat %>%
    ggplot() +
    geom_point(aes(
      x = prop_conserved_sites, y = !!sym(col_aes),
      size = prop_conserved_pos, colour = !!sym(col_aes)
    )) +
    facet_wrap(vars(`Role in study`)) +
    geom_label_repel(
      aes(
        x = prop_conserved_sites, y = !!sym(col_aes),
        label = TF
      ),
      data = sig_mat_cons_sum_metadat %>%
        filter(prop_conserved_sites > 0.05),
      min.segment.length = 0.1,
      box.padding = 0.15,
      label.padding = 0.15,
      max.overlaps = 7,
      label.size = 0.1,
      size = 2.5
    ) +
    guides(
      size = guide_legend(title = "Proportion of conserved bases"),
      colour = legend_fn
    ) +
    labs(
      x = "Proportion of binding sites that are conserved",
      y = "TF Age (Mya)"
    )
  if (col_aes == "max_age") {
    p <- p + scale_colour_viridis_c(
      option = "cividis",
      breaks = c(15.2, 73.8, 180, 351.7, 563.4),
      labels = c(
        "Great Apes (15.20)", "Primates (73.8)",
        "Mammals (180)",
        "Tetrapods (351.7)", "Vertebrates (563.4)"
      )
    )
  }
  ggsave(
    file.path(
      hmap_outs, "figs",
      paste0(
        "sig_conservation_scatterplot_", file_suffix, "_", filt_set,
        ".pdf"
      )
    ),
    p,
    width = 8, height = 3.5
  )
}

# write out counts of various relevant categories - just conserved sites
for (i in c("cluster_100_tfbs_nm", "TF_tfbs_nm")) {
  tab_nm <- ifelse(i == "cluster_100_tfbs_nm", "clusters", "sites")
  stat_pass_summ_cons_df_clusts_prom_phast %>%
    filter(conserved_tfbs) %>%
    group_by(cluster_100_tfbs_nm) %>%
    mutate(`TF category` = ifelse(is.na(`TF category`), "Other", `TF category`)) %>%
    # slice_sample(n = 1) %>%
    transmute(
      .data[[i]],
      `Cluster type` = ifelse(grepl(";", cluster_100_tfbs_nm), "Cluster",
        "Solo binding site"
      ),
      `Ambiguity type` = ifelse(unambiguous_tfbs_100, "Single-TF", "Multi-TF"),
      `TF category` = paste(sort(unique(`TF category`)), collapse = ",")
    ) %>%
    group_by(`Cluster type`, `Ambiguity type`, `TF category`) %>%
    unique() %>%
    summarize("# of {tab_nm}" := n_distinct(.data[[i]])) %>%
    write_tsv(file.path(
      proj_dir, "analysis/computematrix_out/all_PWM_hits/hmaps",
      paste0(filt_set, "_conserved_cluster_100_", tab_nm, ".tsv")
    ))
}
# write out counts of various relevant categories - all sites
for (j in c("cluster_1000_tfbs_nm", "cluster_200_tfbs_nm")) {
  for (i in c(j, "TF_tfbs_nm")) {
    clust_dist <- gsub("cluster_([0-9]+)_tfbs_nm", "\\1", j)
    tab_nm <- ifelse(grepl("cluster", i), "clusters", "sites")
    unamb <- ifelse(grepl("1000", j), "unambiguous_tfbs_1000", "unambiguous_tfbs_200")
    stat_pass_summ_cons_df_uncons_cons_clusts_prom_phast %>%
      # lazy_dt() %>%
      group_by(.data[[j]]) %>%
      mutate(`TF category` = ifelse(is.na(`TF category`), "Other", `TF category`)) %>%
      # slice_sample(n = 1) %>%
      transmute(
        .data[[i]],
        `Cluster type` = ifelse(grepl(";", .data[[j]]), "Cluster", "Solo binding site"),
        `Ambiguity type` = ifelse(.data[[unamb]], "Single-TF", "Multi-TF"),
        `TF category` = paste(sort(unique(`TF category`)), collapse = ",")
      ) %>%
      group_by(`Cluster type`, `Ambiguity type`, `TF category`) %>%
      unique() %>%
      summarize("# of {tab_nm}" := n_distinct(.data[[i]])) %>%
      # as_tibble() %>%
      write_tsv(file.path(
        proj_dir, "analysis/computematrix_out/all_PWM_hits/hmaps",
        paste0(
          filt_set, "_conserved+unconserved_cluster_",
          clust_dist, "_", tab_nm, ".tsv"
        )
      ))
  }
}


THE_BIG_NUMBER <- sig_mat_cons_sum_metadat %>%
  filter(`Role in study` == "Codebook TF") %>%
  pull(n_conserved_pos) %>%
  sum()
THE_BIGGER_NUMBER <- sig_mat_cons_sum_metadat %>%
  filter(`Role in study` == "Codebook TF") %>%
  pull(n_positions_in_conserved_sites) %>%
  sum()

for (tf_type in c("Codebook TF", "Control")) {
  filt_set <- sig_mat_cons_sum_metadat %>%
    filter(`Role in study` == tf_type)
  n_conserved_pos_sum <- filt_set %>%
    pull(n_conserved_pos) %>%
    sum()
  n_conserved_sites_sum <- filt_set %>%
    pull(n_conserved_sites) %>%
    sum()
  n_positions_in_conserved_sites_sum <- filt_set %>%
    pull(n_positions_in_conserved_sites) %>%
    sum()
  print(
    paste0(
      tf_type, " conservation: ", n_conserved_pos_sum,
      " conserved positions in ", n_conserved_sites_sum,
      " binding sites containing ", n_positions_in_conserved_sites_sum,
      " positions total"
    )
  )
}
