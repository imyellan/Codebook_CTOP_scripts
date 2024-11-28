library(dplyr)
library(ggplot2)
library(readr)
proj_dir <-
  ifelse(grepl("MacBook", Sys.info()["nodename"]),
         "~/iyellan_bc2-clust/codebook_proj/tfbs_conservation",
         "/home/hugheslab1/iyellan/codebook_proj/tfbs_conservation"
  )

## annotate just the clusters based on cluster only intersections with features
stat_pass_summ_cons_df_clusts_prom_phast <- 
  read_csv(file.path(proj_dir, "analysis/computematrix_out/all_PWM_hits/hmaps/conservation_all_sites_trip_opt_moods_HC.csv.gz"))
cluster_info <- stat_pass_summ_cons_df_clusts_prom_phast %>%
  group_by(cluster_100_tfbs_nm) %>%
  summarise(cluster_100_chr = unique(cluster_100_chr), 
            cluster_100_start = unique(cluster_100_start), 
            cluster_100_end = unique(cluster_100_end), 
            cluster_100_strands = unique(cluster_100_strands), 
            cluster_100_length = unique(cluster_100_length),
            cluster_n_sites = n_distinct(TF_tfbs_nm),
            mean_phylop = mean(mean_phyloP),
            max_phyloP = max(max_phyloP),
            TFs = paste(unique(TF), collapse = ";"),
            n_TFs = n_distinct(TF_tfbs_nm),
            n_distinct_TFs = n_distinct(TF))

# cluster intersect results 
cluster_cols <- paste0("cluster_100_", 
                       c("chr", "start", "end", "strands", "tfbs_nm", "id"))
clust_prom_df <- 
  read_tsv(file.path(proj_dir, 
                     "analysis/phastcons_promoter_intersect/promoters_intersect_ctop_clust.bed"), 
           col_names = c(cluster_cols, "promoter_chr", 
                         "promoter_start", "promoter_end", 
                         "promoter_gene", "promoter_score", "promoter_strand", "promoter_category", "promoter_intersect_length")) %>%
  select(cluster_100_tfbs_nm, promoter_intersect_length, promoter_gene, 
         promoter_category) %>%
  unique()

clust_cpg_df <- 
  read_tsv(file.path(proj_dir, 
                     "analysis/phastcons_promoter_intersect/cpgs_unmasked_intersect_ctop_clust.bed"), 
           col_names = c(cluster_cols, "cpg_chr", 
                         "cpg_start", "cpg_end", 
                         "cpg_name", "cpg_gene", 
                         "cpg_category", "cpg_intersect_length")) %>%
  group_by(cluster_100_tfbs_nm) %>%
  slice_max(cpg_intersect_length, with_ties = F) %>%
  select(
    cluster_100_tfbs_nm, cpg_name, cpg_gene, cpg_category, cpg_intersect_length
  ) %>% unique()

clust_phastcon_df <- 
  read_tsv(file.path(proj_dir, 
                     "analysis/phastcons_promoter_intersect/phastcons_intersect_ctop_clust.bed"), 
           col_names = c(cluster_cols, "phastcons_chr", 
                         "phastcons_start", "phastcons_end", 
                         "phastcons_lod", "phastcons_score", 
                         "phastcons_intersect_length")) %>%
  group_by(cluster_100_tfbs_nm) %>%
  transmute(
    cluster_100_tfbs_nm,
    phastcons_intersect_length_sum = sum(phastcons_intersect_length)
  ) %>% unique()
  
clust_enhancer_df <- 
  read_tsv(file.path(proj_dir, 
                     "analysis/phastcons_promoter_intersect/enhancer_intersect_ctop_clust.bed"), 
           col_names = c(cluster_cols, "enhancer_chr", 
                         "enhancer_start", "enhancer_end", 
                         "enhancer_id", "enhancer_score", 
                         "enhancer_strand", "enhancer_intersect_length")) %>%
  transmute(cluster_100_tfbs_nm,
    enhancer_nm = paste(enhancer_chr, enhancer_start, enhancer_end, sep = "_"), enhancer_intersect_length) %>%
  group_by(cluster_100_tfbs_nm) %>%
  slice_max(enhancer_intersect_length, with_ties = F)

simple_sat_repeats <- c("Simple_repeat", 
                        "Low_complexity", "Satellite", 
                        "RNA", "rRNA", "scRNA", "snRNA", "srpRNA", "tRNA")
clust_repeat_df <- 
  read_tsv(file.path(proj_dir, 
                     "analysis/phastcons_promoter_intersect/repeats_intersect_ctop_clust.bed"), 
           col_names = c(cluster_cols, "repeat_chr", 
                         "repeat_start", "repeat_end", 
                         "repeat_family", "repeat_score",
                         "repeat_strand",
                         "repeat_class", "repeat_id", 
                         "repeat_intersect_length")) %>%
  group_by(cluster_100_tfbs_nm) %>%
  mutate(repeat_type_priority = 
           case_when(all(repeat_class %in% c(simple_sat_repeats, ".")) &
                       all(repeat_class != ".") ~ 1,
                     any(!repeat_class %in% c(simple_sat_repeats, ".")) ~ 2,
                     TRUE ~ 0)) %>%
  slice_max(repeat_type_priority) %>%
  slice_max(repeat_intersect_length, with_ties = F) %>%
  select(cluster_100_tfbs_nm, repeat_family, repeat_class, 
         repeat_id, repeat_type_priority,
         repeat_intersect_length)

## join together
cluster_info_annot <- cluster_info %>%
  left_join(clust_prom_df) %>%
  left_join(clust_cpg_df) %>%
  left_join(clust_phastcon_df) %>%
  left_join(clust_enhancer_df) %>%
  left_join(clust_repeat_df) %>% 
  group_by(cluster_100_tfbs_nm) %>%
  mutate(promoter_genes = paste(unique(promoter_gene), collapse = ";"),
         promoter_categories = paste(unique(promoter_category), collapse = ";"),
         promoter_intersect_lengths = paste(promoter_intersect_length, collapse = ";"),
         .keep = "unused") %>%
  unique() %>%
  filter(!is.na(cluster_100_tfbs_nm))

cluster_info_annot %>% write_csv(file.path(proj_dir, "analysis/computematrix_out/all_PWM_hits/hmaps/CTOP_clust_100_trip_opt_moods_HC_annots.csv.gz"))
  