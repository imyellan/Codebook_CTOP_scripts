lrt_mat <- function(matrixout_mat, motif_positions, fig_out) {
  IC_vals <- motif_IC_vals %>%
    filter(TF == TF_nm) %>%
    pull(ics) %>%
    unlist()
  TF_IC_phylo_coeff <- ic_phylop_corrs %>%
    filter(TF == TF_nm) %>%
    pull(lm_coeff)
  motif_size <- length(motif_positions)
  
  # annotate binding sites with overlapping repeats
  # matrixout_mat_reps <- matrixout_mat %>% select(chr, start, end, strand) %>%
  #   left_join(repeat_intersection,
  #             by = c("chr" = "tf_chr", "start" = "tf_start", "end" = "tf_end",
  #                    "strand" = "tf_strand")) %>%
  #   group_by(chr, start, end, strand) %>%
  #   slice_max(overlap_length, with_ties = F) %>% ungroup() %>%
  #   mutate(simp_low_comp = grepl("Simple_repeat|Low_complexity", repeat_class)) %>%
  #   select(-c(repeat_chr, repeat_start, repeat_end, repeat_strand))
  # simp_low_comp_sites <- matrixout_mat_reps$simp_low_comp
  
  # prepare matrices
  matrixout_mat_pos <- matrixout_mat %>%
    select(-all_of(non_pos_cols[1:6])) #%>%
    #filter(if_all(everything(), ~ !is.na(.)))
  
  # create a matrix of the highest phyloP scoring positions
  motif_phylop_all <- matrixout_mat_pos %>%
    select(as.character(motif_positions))
  IC_corrs <- apply(motif_phylop_all, 1, function(x) cor(x, IC_vals))
  motif_phylop <- cbind(motif_phylop_all, IC_corrs = IC_corrs) %>%
    filter(if_any(-IC_corrs, ~ . > -log10(0.1)), !is.na(IC_corrs)) %>%
    slice_max(IC_corrs, # prop = 0.1,
              n = 100
    ) %>%
    select(-IC_corrs)
  
  # if nrows(motif_phylop) == 0, set motif_phylop to 0
  if (nrow(motif_phylop) == 0) {
    # return(NULL)
    motif_phylop <- matrix(0, nrow = 2, ncol = motif_size)
  }
  else{
    # visual motif_phylop
    pheatmap(
      as.matrix(motif_phylop),
      cluster_cols = F,
      # cluster_rows = F,
      # use_raster = F,
      color = colorRampPalette(brewer.pal(
        n = 7,
        name = "RdBu"
      ))(16)
    )
  }
    
  motif_phylop_mean <- colMeans(motif_phylop, na.rm = T)
  motif_phylop_sd <- apply(motif_phylop, 2, function(x) sd(x, na.rm = T))
  
  if (any(is.na(motif_phylop_sd))) {
    # print(paste0("ERROR: NA in motif_phylop_sd for ", TF_nm))
    # return(NULL)
    # set motif_phylop_sd to 0
    motif_phylop_sd <- rep(0, motif_size)
  }
  
  # score each position matrix using the profile of motif phyloP means
  lrt_mat_phylop <- calc_lrt(
    matrixout_mat_pos,
    motif_phylop_mean, motif_phylop_sd, motif_size,
    IC_vals,
    TF_IC_phylo_coeff, "phylopMean"
  )
  lrt_mat_IC <- calc_lrt(
    matrixout_mat_pos,
    motif_phylop_mean, motif_phylop_sd, motif_size,
    IC_vals,
    TF_IC_phylo_coeff, "icVals"
  )
  
  # calculate the 2nd derivative of the LRT score at the motif start position
  lrt_mat_phylop_df <- calc_2nd_deriv(lrt_mat_phylop, motif_positions)
  lrt_mat_IC_df <- calc_2nd_deriv(lrt_mat_IC, motif_positions)
  
  return(list(lrt_mat_phylop_df, lrt_mat_IC_df, matrixout_mat_pos))
}

analyze_lrt_mats <- function(df, matrixout_mat, 
                             matrixout_mat_pos, i, TF_nm, 
                             motif_positions, fig_out, NA_mat) {
  # plots of LRT scores and other things
  statistical_lrt_cutoffs <- c()
    ## histogram of LRT scores in motif start position
    df %>%
      ggplot(
        aes(x = .data[[paste0("V", as.character(motif_positions[1]))]])
      ) +
      geom_histogram(fill = "grey", colour = "black", binwidth = 20) +
      xlab("LRT score") +
      theme_bw()
    ggsave(file.path(fig_out, paste0(TF_nm, "_lrt_score_hist_", i, ".pdf")))
    
    ## hexplot of LRT scores vs 2nd derivatives
    if (nrow(df) > 1) {
      df %>%
        ggplot(
          aes(
            x = .data[[paste0("V", as.character(motif_positions[1]))]],
            y = .data[["second_derivative"]]
          )
        ) +
        geom_hex() +
        xlab("LRT score") +
        theme_bw()
      ggsave(file.path(fig_out, paste0(TF_nm, "_lrt_score_v_2nd_d_", i, ".pdf")))
    }
    
    ## plot the mean LRT score across all motif positions
    qplot(1:length(colMeans(df)), colMeans(df))
    ggsave(file.path(fig_out, paste0(TF_nm, "_lrt_score_means_", i, ".pdf")))
    
    ## number of binding sites at different LRT score and 2nd derivative thresholds
    # lrt_opt_summary <- map_dfr(
    #   seq(-50, 500, by = 0.1),
    #   # -0.1,
    #   ~ df %>%
    #     filter(.data[[paste0(
    #       "V",
    #       as.character(motif_positions[1])
    #     )]] < .x) %>%
    #     summarize(
    #       lrt_thresh = .x,
    #       n_passing = nrow(.),
    #       mean_sec_dev = mean(second_derivative, na.rm = T),
    #       # median_sec_dev = median(second_derivative, na.rm = T),
    #       # geom_mean_sec_dev = psych::geometric.mean(second_derivative),
    #       # harm_mean_sec_dev = psych::harmonic.mean(second_derivative),
    #       any_sec_dev_0 = any(second_derivative == 0),
    #       # any_sec_dev_negative = any(second_derivative < 0)
    #     )
    # ) %>%
    #   filter(!is.nan(mean_sec_dev))
    
    # figure comparing tfbs 2nd derivative to lrt score
    # lrt_opt_summary %>%
    #   ggplot(aes(x = lrt_thresh, y = mean_sec_dev, colour = n_passing)) +
    #   geom_point()
    # ggsave(file.path(fig_out, paste0(TF_nm, "_lrt_opt_summary_", i, ".pdf")))
    
    # smoothed_peaks_bullshit()
    
    # make heatmaps at different LRT threshes
    sec_dev_threshes <- -10e50 # because negative sec devs are possible and i'm too lazy to remove this check right now
    lrt_thresh <- -10 ## final hard threshold??
    lrt_hmaps(
      df, matrixout_mat_pos,
      i, motif_positions,
      # smoothed_thresh,
      lrt_thresh,
      # c(smoothed_thresh, inflec_point_thresh),
      # c(-25, smoothed_thresh, inflec_point_thresh, seq(0, 300, by = 50)),
      sec_dev_threshes, fig_out, 
      crit_value = 100, NA_mat
    )
    
    # write table
    # write pass and fails to files
    matrixout_mat_lrt <- 
      cbind(matrixout_mat, lrt_score = df[[paste0("V", as.character(motif_positions[1]))]],
            second_derivative = df$second_derivative)
    
    matrixout_mat_lrt %>%
      filter(lrt_score < lrt_thresh) %>%
      select(chr, start, end, id, pwm_score, strand, lrt_score, second_derivative,
             all_of(as.character(motif_positions))) %>% 
      write_tsv(
        file = file.path(fig_out, paste0(
          TF_nm, "_lrt_thresh_", lrt_thresh,
          "_2nd_d_", sec_dev_threshes, "_", i, "_pass.tsv.gz")),
        num_threads = 1
        )
    
    rbind(matrixout_mat_lrt %>%
      filter(lrt_score >= lrt_thresh),
      NA_mat %>% mutate(lrt_score = NA, second_derivative = NA)) %>%
      select(chr, start, end, id, pwm_score, strand, lrt_score, second_derivative,
             all_of(as.character(motif_positions))) %>% 
      write_tsv(
        file = file.path(fig_out, paste0(
          TF_nm, "_lrt_thresh_", lrt_thresh,
          "_2nd_d_", sec_dev_threshes, "_", i, "_fail.tsv.gz")),
        num_threads = 1
        )
    
    ## first try fitting sigmoid curve to the data
    # scaled_x <- scale(x)
    # scaled_y <- scale(y)
    # fit <- try(nls(scaled_y ~ SSlogis(scaled_x, Asym, xmid, scal),
    #            data = data.frame(scaled_x = scaled_x,
    #                              scaled_y = scaled_y),
    #            control = list(algorithm = "port")))
    # if (inherits(fit, "try-error")) {
    #   # just use 0 which is already a cutoff being examined
    #   statistical_lrt_cutoffs <- c(statistical_lrt_cutoffs, 0)
    #   } else {
    #   # else if the fit succeeds, use the x value where the curve is 0.5,
    #   # determine the bend points of the fitted curve
    #     upper_asym <- unname(coef(fit)["Asym"])
    #     xmid <- unname(coef(fit)["xmid"])
    #     scal <- unname(coef(fit)["scal"])
    #     ## get a threshold value for getting x as close to the asymptote as possible
    #     threshold <- 0.999 * upper_asym
    #     ## get the value of x at threshold of the asymptote
    #     x_near_upper_asymptote <-
    #       uniroot(function(x) SSlogis(x, upper_asym, xmid, scal) - threshold,
    #               interval = c(-1000, 1000))$root
    #     ## get the difference between the xmid and x_near_asymptote
    #     xmid_diff <- xmid - x_near_upper_asymptote
    #     ## add xmid_diff to xmid to get a value of x that is close to the lower
    #     ## asymptote
    #     x_near_lower_asymptote <- xmid + xmid_diff
    #     ## get the value of y at x_lower - this represents an estimate of the
    #     ## lower asymptote
    #     lower_asym <- SSlogis(x_near_lower_asymptote, upper_asym, xmid, scal)
    #     plot(scaled_y ~ scaled_x)
    #     curve(SSlogis(x, upper_asym, xmid, scal), add = TRUE, col = "red")
    #     ## next use the formula for finding the bend points of the sigmoid from
    #     ## Sebaugh and McCray (2003)
    #     k <- 4.6805
    #     y_bend_lower = ((upper_asym - lower_asym) / (1 + (1/k))) + lower_asym
    #     y_bend_higher = ((upper_asym - lower_asym) / (1 + k)) + lower_asym
    #     ## get the x values at the bend points
    #     x_bend_lower <- uniroot(
    #       function(x) SSlogis(x, upper_asym, xmid, scal) - y_bend_lower,
    #       interval = c(-1000, 1000))$root
    #     x_bend_higher <- uniroot(
    #       function(x) SSlogis(x, upper_asym, xmid, scal) - y_bend_higher,
    #       interval = c(-1000, 1000))$root
    #
    #     ## unscale the bend values, plot on the original scale
    #     x_bend_lower_unscaled <- x_bend_lower * sd(x) + mean(x)
    #     x_bend_higher_unscaled <- x_bend_higher * sd(x) + mean(x)
    #     plot(x, y)
    #     abline(v = x_bend_lower_unscaled, col = "blue")
    #     abline(v = x_bend_higher_unscaled, col = "blue")
    #     # save the upper bend point as the LRT cutoff
    #     statistical_lrt_cutoffs <-
    #       c(statistical_lrt_cutoffs, x_bend_higher_unscaled)
    #   }
    # ## also try to partition into linear segments
    # x_cuts <- RcppDynProg::solve_for_partition(x, y, penalty = 0.6)
    # # get slopes of each group
    # x_cuts <- x_cuts %>% group_by(group) %>%
    #   mutate(slope = (pred[2] - pred[1]) / (x[2] - x[1]))
    # # get minimum x value where slope is > -1e-3, use as LRT cutoff
    # statistical_lrt_cutoffs <-
    #   c(statistical_lrt_cutoffs, min(x_cuts$x[x_cuts$slope > -5e-2 &
    #                                             sign(x_cuts$slope) == -1 &
    #                                             x_cuts$x > -47]))
  
  # try different thresholds for the LRT score and 2nd derivative (representing
  # peak sharpness), make heatmaps of the results
  # lrt_threshes <- c(-25, seq(0, 300, by = 50))
  # sec_dev_threshes <- seq(0, 200, by = 50)
  # lrt_bendpoints: attempt to determine good lrt cutoff
  # lrt_hmaps(
  #   lrt_mat_phylop_df, matrixout_mat_pos,
  #   "phylopMean", motif_positions,
  #   c(-25, statistical_lrt_cutoffs[1:2], seq(0, 300, by = 50)),
  #   sec_dev_threshes, fig_out
  # )
  # lrt_hmaps(
  #   lrt_mat_IC_df, matrixout_mat_pos,
  #   "icVals", motif_positions,
  #   c(-25, statistical_lrt_cutoffs[3:4],
  #     seq(0, 500, by = 50)), sec_dev_threshes, fig_out
  # )
}

calc_lrt <- function(matrixout_mat_pos,
                     motif_phylop_mean, motif_phylop_sd,
                     motif_size, IC_vals, TF_IC_phylo_coeff,
                     stat) {
  
  # if a value in motif_phylop_sd == 0, replace it with 1 to prevent dividing by 0
  motif_phylop_sd[motif_phylop_sd == 0] <- 1
  
  # for the null model, take the mean of all the phylop scores in the matrix
  null_mean <- mean(as.matrix(matrixout_mat_pos %>% 
                                select(names(motif_phylop_mean))), na.rm = T)
  null_sd <- sd(as.matrix(matrixout_mat_pos %>% 
                            select(names(motif_phylop_mean))), na.rm = T)
  # create a matrix to store the LRT scores for each position in the motif
  matrixout_mat_norm <-
    matrix(NA, nrow(matrixout_mat_pos), 100 - motif_size)

  # also create a null/neutral model that assumes the mean == 0 and sd is the same as the motif model
  matrixout_mat_null <-
    matrix(NA, nrow(matrixout_mat_pos), 100 - motif_size)

  # score each position in the matrices using the motif model
  for (pos in 1:(100 - motif_size)) {
    window_mat <- as.matrix(matrixout_mat_pos[, pos:(pos + motif_size - 1)])
    # score window using either the PhyloP profile or the IC profile
    if (stat == "phylopMean") {
      matrixout_mat_norm[, pos] <-
        rowSums((sweep(window_mat, 2, motif_phylop_mean, "-"))^2 %>%
          sweep(2, motif_phylop_sd^2, "/"))
    } else if (stat == "icVals") {
      # score window
      matrixout_mat_norm[, pos] <-
        rowSums((sweep(window_mat, 2, IC_vals * TF_IC_phylo_coeff, "-"))^2 %>%
          sweep(2, motif_phylop_sd^2, "/"))
    }
    # score window under null model
    null_means <- rep(0, ncol(window_mat))
    matrixout_mat_null[, pos] <-
      rowSums((sweep(window_mat, 2, null_means, "-"))^2 %>%
        sweep(2, motif_phylop_sd^2, "/"))
    # matrixout_mat_null[, pos] <-
    #   rowSums((window_mat - null_mean)^2/(null_sd^2))
  }

  # calculate the LRT statistic by comparing the null model to the motif model
  lrt_mat <- -2 * (matrixout_mat_null - matrixout_mat_norm)
  # change any Inf values to 0
  lrt_mat[is.infinite(lrt_mat)] <- 0
  lrt_mat
}

smoothed_peaks_bullshit <- function(){
  # examining this relationship, specifically first derivative of the spline
  spline_fn <- splinefun(
    x = lrt_opt_summary$lrt_thresh,
    y = lrt_opt_summary$mean_sec_dev
  )
  smth_window <- 0.3
  smoothed_df <- data.frame(
    lrt_thresh = lrt_opt_summary$lrt_thresh,
    mean_sec_dev = lrt_opt_summary$mean_sec_dev) %>%
    mutate(
      spline = spline_fn(lrt_thresh),
      spline_first_deriv = spline_fn(lrt_thresh, deriv = 1),
      spline_second_deriv = spline_fn(lrt_thresh, deriv = 2),
      # sma_smoothed_spline_first_deriv = 
      #   smth(spline_first_deriv, method = "sma", window = 1000),
      gauss_smoothed_spline_first_deriv = 
        smth.gaussian(spline_first_deriv, smth_window)
    )
  
  # find peaks
  smoothed_peaks <-
    findpeaks(-1*smoothed_df$gauss_smoothed_spline_first_deriv,
              nups = 5, sortstr = TRUE)
  
  # optimize smoothing window size for peak size
  top_score <- 0
  furthest_peak_position <- 0
  for (smth_window in seq(0.9, 1e-4, by = -0.0005)) {
    smoothed_peaks <- findpeaks(-1 * smoothed_df %>%
                                  mutate(gauss_smoothed_spline_first_deriv = 
                                           smth.gaussian(spline_first_deriv, smth_window)) %>%
                                  pull(gauss_smoothed_spline_first_deriv),
                                nups = 10, sortstr = TRUE)
    if (is.null(smoothed_peaks)) {
      next
    }
    if (smoothed_peaks[1,1] > top_score & 
        smoothed_peaks[1,2] > furthest_peak_position) {
      top_score <- smoothed_peaks[1,1]
      furthest_peak_position <- smoothed_peaks[1,2]
      top_smth_window <- smth_window
    }
  }
  
  smoothed_peaks <- 
    findpeaks(-1 * smoothed_df %>%
                mutate(gauss_smoothed_spline_first_deriv = 
                         smth.gaussian(spline_first_deriv, 
                                       top_smth_window)) %>%
                pull(gauss_smoothed_spline_first_deriv),
              nups = 10, sortstr = TRUE)
  
  if(is.null(smoothed_peaks)) {
    for(j in 4:3){
      smoothed_peaks <-
        findpeaks(-1*smoothed_df$gauss_smoothed_spline_first_deriv,
                  nups = j, sortstr = TRUE)
    }
    if(is.null(smoothed_peaks)){
      # if at nups == 2 there are still no peaks, try adjusting the smoothing
      smth_window <- 1e-2
      smoothed_df <- smoothed_df %>%
        mutate(gauss_smoothed_spline_first_deriv = 
                 smth.gaussian(spline_first_deriv, smth_window))
      for (j in 5:1) {
        smoothed_peaks <-
          findpeaks(-1*smoothed_df$gauss_smoothed_spline_first_deriv,
                    nups = j, sortstr = TRUE)
        if(!is.null(smoothed_peaks)){
          break
        }
      }
    }
  }
  
  if(is.null(smoothed_peaks)){
    print(paste0("peak detection failed for ", TF_nm, " ", i))
    next
  }
  # select best peak
  if (nrow(smoothed_peaks) > 1) {
    smoothed_thresh <- 
      ifelse(smoothed_peaks[2,1]/smoothed_peaks[1,1] > 0.3 &
               smoothed_peaks[2,2] > smoothed_peaks[1,2],
             smoothed_df$lrt_thresh[smoothed_peaks[2,2]],
             smoothed_df$lrt_thresh[smoothed_peaks[1,2]])
  } else{
    smoothed_thresh <- smoothed_df$lrt_thresh[smoothed_peaks[1,2]]
  }
  # the resmoothed stuff
  if (smth_window == 1e-2) {
    smoothed_thresh <- smoothed_df$lrt_thresh[smoothed_peaks[1,4]]
  }
  smoothed_df %>% 
    ggplot(aes(x = lrt_thresh, 
               y = smth.gaussian(spline_first_deriv, top_smth_window),
               # y = smth(spline_first_deriv, method = "sma", window = 1000)
               # y = gauss_smoothed_spline_first_deriv
    )) +
    geom_point(size = 0.5) +
    # xlim(c(-5,20)) +
    theme_bw()
  ggsave(file.path(fig_out, paste0(TF_nm, "_", i, "_smoothing.pdf")))
  
  # Attempt to find good LRT threshold by defining an upper boundary LRT value
  # where the relationship between LRT and 2nd dev are still linear and
  # negatively sloped
  x <- lrt_opt_summary$lrt_thresh
  y <- lrt_opt_summary$mean_sec_dev
  
  ## try using first point where a second derivative is 0 as threshold -
  ## this represents inflection point
  # empirical_fn <- splinefun(x = lrt_opt_summary$lrt_thresh,
  #                           y = lrt_opt_summary$mean_sec_dev)
  # empirical_fn(x, deriv = 1)
  # plot(x, y)
  # plot(x, empirical_fn(x))
  # plot(x, empirical_fn(x, deriv = 2))
  inflec_point_thresh <-
    min(lrt_opt_summary$lrt_thresh[lrt_opt_summary$any_sec_dev_0])
}


calc_2nd_deriv <- function(lrt_mat, motif_positions) {
  motif_start <- paste0("V", as.character(motif_positions[1]))
  start_plus_2 <- paste0("V", as.character(motif_positions[1] + 2))
  start_minus_2 <- paste0("V", as.character(motif_positions[1] - 2))

  lrt_mat_df <- as.data.frame(lrt_mat) %>%
    mutate(second_derivative = !!sym(start_minus_2) +
      !!sym(start_plus_2) - 2 * !!sym(motif_start))
  lrt_mat_df
}

lrt_hmaps <- function(lrt_mat_df, matrixout_mat_pos,
                      stat,
                      motif_positions,
                      lrt_threshes, sec_dev_threshes,
                      fig_out, crit_value, NA_mat) {
  q_01 <- round(-log10(0.1), digits = 2)
  breaks <- c(seq(-6, -q_01, 1), 0, seq(q_01, 6, 1))
  lg_breaks <- c(-6, -1.3, -q_01, 0, q_01, 6, 2)
  max_length <- min(
    100 - length(motif_positions),
    motif_positions[length(motif_positions)]
  )

  for (lrt_thresh in lrt_threshes) {
    for (second_derivative_thresh in sec_dev_threshes) {
      lrt_mat_df <- lrt_mat_df %>% 
        mutate(
        meets_lrt_thresh =
          .data[[paste0("V", as.character((motif_positions[1]):(motif_positions[1])))]] < lrt_thresh,
        meets_second_derivative_thresh =
          second_derivative >= second_derivative_thresh,
        meets_both_thresh = meets_lrt_thresh & meets_second_derivative_thresh
      )
      passing_pos <- lrt_mat_df %>% pull(meets_both_thresh)
      
      # make annotation df
      annotation_df <- lrt_mat_df %>% 
        transmute(LRT_score = .data[[paste0("V", as.character((motif_positions[1]):(motif_positions[1])))]],
               second_derivative)
      
      # make pass hmap
      if (nrow(matrixout_mat_pos[passing_pos, ]) >= 1) {
        pass_fil_out <- file.path(fig_out, paste0(
          TF_nm,
          "_lrt_thresh_", lrt_thresh,
          "_2nd_d_", second_derivative_thresh, "_",
          stat, ifelse(lrt_thresh == crit_value, "_smoothed_peak_thresh", ""),
          "_pass.pdf"
        ))
        produce_hmap(
          in_mat = matrixout_mat_pos[passing_pos, ],
          mat_nm = paste0(
            "n = ",
            nrow(matrixout_mat_pos[passing_pos, ]), "/",
            (nrow(matrixout_mat_pos) + nrow(NA_mat))
          ),
          breaks = breaks,
          lg_breaks = lg_breaks,
          column_gap_pos = c(motif_positions[1] - 1, max_length),
          row_annots = annotation_df[passing_pos, ],
          outfile = pass_fil_out
      )
      }
     
      # make fail hmap
      if (nrow(matrixout_mat_pos[!passing_pos, ]) >= 1){
        fail_fil_out <- file.path(fig_out, paste0(
          TF_nm,
          "_lrt_thresh_", lrt_thresh,
          "_2nd_d_", second_derivative_thresh, "_",
          stat, ifelse(lrt_thresh == crit_value, "_smoothed_peak_thresh", ""),
          "_fail.pdf"
        ))
        produce_hmap(
          in_mat = rbind(matrixout_mat_pos[!passing_pos, ], 
                         NA_mat %>% select(-c(1:6))),
          mat_nm = paste0(
            "n = ",
            nrow(matrixout_mat_pos[!passing_pos, ]), "/",
            (nrow(matrixout_mat_pos) + nrow(NA_mat))
          ),
          breaks = breaks,
          lg_breaks = lg_breaks,
          column_gap_pos = c(motif_positions[1] - 1, max_length),
          row_annots = rbind(annotation_df[!passing_pos, ], 
                             data.frame(LRT_score = rep(NA, nrow(NA_mat)),
                                        second_derivative = rep(NA, nrow(NA_mat)))),
          outfile = fail_fil_out
        )
      }
    }
  }
}

corr_perm_tfbs_conservation <- 
  function(matrixout_mat, motif_positions, fig_out, lrt_mat_IC_df, NA_mat) {
  require(shuffle)
  require(permute)
  set.seed(1234)

  IC_vals <- motif_IC_vals %>%
    filter(TF == TF_nm) %>%
    pull(ics) %>%
    unlist()
  max_length <- min(
    100 - length(motif_positions),
    motif_positions[length(motif_positions)]
  )
  lrt_mat_IC_df_pos <- lrt_mat_IC_df %>% 
    select(paste0("V", as.character(motif_positions[1])),
           second_derivative)
  colnames(lrt_mat_IC_df_pos) <- c("LRT_score", "second_derivative")
  
  matrixout_mat_pos <- matrixout_mat %>% select(-all_of(non_pos_cols[1:6])) 
  motif_phylop_all <- matrixout_mat %>% select(as.character(motif_positions))

  null_dist <- list()
  n_perms <- ifelse(numPerms(IC_vals) >= 1000, 1000, numPerms(IC_vals))
  permuted_IC_vals <- lapply(1:1000, function(x) IC_vals[shuffle(length(IC_vals))])
  
  correlate_to_phyloPs <- function(motif_phylop_all, IC_vals_shuff){
    if (nrow(motif_phylop_all) > 1500) {
      motif_phylop_all_samp <- motif_phylop_all[sample(nrow(motif_phylop_all),
                                                       size = 1500,
                                                       replace = F), ]
    } else {
      motif_phylop_all_samp <- motif_phylop_all
    }
    row_cor_pearson(motif_phylop_all_samp, IC_vals_shuff, 
                    alternative = "two.sided", conf.level = 0.95)$cor
    # sapply(
    #     lapply(seq_len(nrow(motif_phylop_all_samp)),
    #            function(i) unlist(motif_phylop_all_samp[i,])),
    #     function(x) cor(x, IC_vals_shuff))
  }
  
  plan(multicore)
  null_dist <- permuted_IC_vals %>% 
    future_map(~correlate_to_phyloPs(motif_phylop_all, .x))
  plan(sequential)
  
  # based on null distribution, get score that corresponds to 0.05 level of significance
  null_dist <- unlist(null_dist)
  null_dist <- null_dist[!is.na(null_dist)]
  alpha <- 0.05
  corr_thresh <<-
    pnorm(alpha, mean = mean(null_dist, na.rm = T), 
          sd = sd(null_dist, na.rm = T), lower.tail = F)

  # take correlation of phylop score with IC score
  cors <- row_cor_pearson(motif_phylop_all, IC_vals, 
                          alternative = "two.sided", conf.level = 0.95)$cor
  mean_phyloPs <- rowMeans(motif_phylop_all, na.rm = T)
  motif_phylop_all_corrs <- motif_phylop_all %>%
    # rowwise() %>%
    mutate(
      mean_phyloP = mean_phyloPs,
      phylop_IC_corr = cors,
      corr_above_thresh = phylop_IC_corr >= corr_thresh
    ) %>%
    ungroup() %>%
    mutate(corr_above_thresh = ifelse(is.na(phylop_IC_corr), F, corr_above_thresh))

  motif_phylop_all_corrs_ltr <- cbind(motif_phylop_all_corrs, lrt_mat_IC_df_pos)
  
  if (length(unlist(null_dist)) > 0) {
    pdf(file.path(fig_out, paste0(TF_nm, "_corr_thresh_hist.pdf")))
    hist(unlist(null_dist))
    abline(v = corr_thresh, col = "red")
    mtext(paste0("Correlation threshold = ", round(corr_thresh, 3)), side = 1, 
          line = 2)
    dev.off()
  }
  
  
  motif_phylop_all_corrs_ltr %>% 
    ggplot(aes(x = LRT_score, y = phylop_IC_corr)) +
    geom_point(aes(color = corr_above_thresh), size = 0.1)
  # write table of binding sites with correlations
  binding_sit_corrs_lrts <- rbind(
    cbind(matrixout_mat, motif_phylop_all_corrs_ltr %>% 
            select(mean_phyloP, phylop_IC_corr, corr_above_thresh, 
                   LRT_score, second_derivative)),
    tibble(NA_mat, mean_phyloP = rep(NA, nrow(NA_mat)),
               phylop_IC_corr = rep(NA, nrow(NA_mat)),
               corr_above_thresh = rep(NA, nrow(NA_mat)),
               LRT_score = rep(NA, nrow(NA_mat)),
               second_derivative = rep(NA, nrow(NA_mat))))
  
  write_csv(binding_sit_corrs_lrts,
            file.path(fig_out, paste0(TF_nm, "_corrvals.csv.gz")),
            num_threads = 1
            )
  
  # save conserved and non-conserved to file
  ## first remove existing files
  file.remove(list.files(file.path(fig_out), 
                         pattern = paste0(TF_nm, "_corr_thresh_0"),
                         full.names = T))
  matrixout_mat_corr <- cbind(matrixout_mat, 
                              motif_phylop_all_corrs_ltr %>% 
                                select(phylop_IC_corr, corr_above_thresh)) %>%
    mutate(corr_thresh = corr_thresh)
  
  write_tsv(matrixout_mat_corr %>% 
              filter(corr_above_thresh) %>%
              select(chr, start, end, id, pwm_score, strand, phylop_IC_corr, 
                     corr_thresh,
                            all_of(as.character(motif_positions))),
            file.path(fig_out, 
                      paste0(TF_nm, "_corr_thresh_", round(corr_thresh, 2),
                             "_", alpha, "_pass.tsv.gz")),
            num_threads = 1
            )
  write_tsv(
    rbind(matrixout_mat_corr %>% filter(!corr_above_thresh) %>%
            select(chr, start, end, id, pwm_score, strand, 
                   phylop_IC_corr, corr_thresh,
                   all_of(as.character(motif_positions))),
                  NA_mat[,c("chr", "start", "end", "id", "pwm_score", "strand",
                            as.character(motif_positions))] %>%
                    mutate(phylop_IC_corr = NA, corr_thresh = NA)),
            file.path(fig_out, paste0(TF_nm, "_corr_thresh_", round(corr_thresh, 2),
                                      "_", alpha, "_fail.tsv.gz")),
    num_threads = 1
    )

  # make hmap of pass and fails
  q_01 <- round(-log10(0.1), digits = 2)
  breaks <- c(seq(-6, -q_01, 1), 0, seq(q_01, 6, 1))
  lg_breaks <- c(-6, -1.3, -q_01, 0, q_01, 1.3, 6)
  
  ## passes
  if(nrow(matrixout_mat_pos[motif_phylop_all_corrs_ltr$corr_above_thresh, ]) > 0){
    pass_fil_out <-
      file.path(fig_out, paste0(TF_nm, "_corr_thresh_", round(corr_thresh, 2),
                                "_", alpha, "_pass.pdf"))
    produce_hmap(matrixout_mat_pos[motif_phylop_all_corrs_ltr$corr_above_thresh, ],
      paste0(
        "n = ",
        nrow(matrixout_mat_pos[motif_phylop_all_corrs_ltr$corr_above_thresh, ]), "/",
        nrow(matrixout_mat_pos)
      ),
      breaks = breaks,
      lg_breaks = lg_breaks,
      column_gap_pos = c(motif_positions[1] - 1, max_length),
      row_annots = motif_phylop_all_corrs_ltr %>% 
        filter(corr_above_thresh) %>% 
        transmute(corr_coeff = phylop_IC_corr, LRT_score, second_derivative),
      outfile = pass_fil_out
    )
  }
  
  # fails
  fail_fil_out <-
    file.path(fig_out, paste0(TF_nm, "_corr_thresh_", round(corr_thresh, 2),
                              "_", alpha, "_fail.pdf"))
  produce_hmap(matrixout_mat_pos[!motif_phylop_all_corrs_ltr$corr_above_thresh, ],
    paste0(
      "n = ",
      nrow(matrixout_mat_pos[!motif_phylop_all_corrs_ltr$corr_above_thresh, ]),
      "/",
      nrow(matrixout_mat_pos)
    ),
    breaks = breaks,
    lg_breaks = lg_breaks,
    column_gap_pos = c(motif_positions[1] - 1, max_length),
    row_annots = motif_phylop_all_corrs_ltr %>% 
      filter(!corr_above_thresh) %>% 
      transmute(corr_coeff = phylop_IC_corr, LRT_score, second_derivative),
    outfile = fail_fil_out
  )
}

produce_hmap <- function(in_mat, mat_nm, breaks = NA, 
                         lg_breaks = NA, column_gap_pos = NULL, 
                         row_annots = NA, outfile) {
  require(rsvg)
  require(ComplexHeatmap)
  require(RColorBrewer)
  require(grImport2)
  #in_mat <- in_mat %>% filter(if_all(everything(), ~!is.na(.)))
  split_cons <- any(grepl("conserved_tfbs", colnames(row_annots)))
  if (nrow(in_mat) > 5000) {
    sample_inds <- sample(1:nrow(in_mat), 5000, replace = F)
    in_mat <- in_mat[sample_inds,]
    if (is.object(row_annots)) {
      row_annots <- row_annots[sample_inds,]
    }
  }
  
  # create row annotation
  library(circlize)
  if (is.object(row_annots)) {
    if (any(grepl("corr_coeff", colnames(row_annots)))) {
      LRT_scores <- row_annots$LRT_score
      sec_devs <- row_annots$second_derivative
      
      col_fun1 <- colorRamp2(c(-1, 0, 1), c("green", "white", "purple"))
      col_fun2 <- colorRamp2(c(min(LRT_scores, na.rm = T), 0,
                               max(LRT_scores, na.rm = T)), 
                             c("cyan", "white", "red"))
      
      if(length(unique(c(min(sec_devs, na.rm = T), 0, 
                         max(sec_devs, na.rm = T)))) >= 2){
        col_fun3 <- colorRamp2(c(min(sec_devs, na.rm = T), 0, 
                                 max(sec_devs, na.rm = T)), 
                               c("blue", "white", "orange"))
      } else{
        col_fun3 <- setNames("white", min(sec_devs, na.rm = T))
      }
      
      rowAnnot <- rowAnnotation(
        corr_coeff = row_annots$corr_coeff,
        LRT_score = row_annots$LRT_score,
        second_derivative = row_annots$second_derivative,
        col = 
          list(corr_coeff = col_fun1, 
               LRT_score = col_fun2, 
               second_derivative = col_fun3),
        annotation_legend_param = list(
          corr_coeff = list(
            at = c(-1, 0, corr_thresh, 1),
            labels = c(-1, 0, round(corr_thresh, 3), 1)
            )
        )
      )
    } else if (any(grepl("LRT_score", colnames(row_annots)))){
      LRT_scores <- row_annots$LRT_score
      sec_devs <- row_annots$second_derivative
      
      col_fun1 <- colorRamp2(c(min(LRT_scores, na.rm = T), 0, 
                               max(LRT_scores, na.rm = T)), 
                             c("green", "white", "red"))
      
      if(length(unique(c(min(sec_devs, na.rm = T), 0, 
                         max(sec_devs, na.rm = T)))) >= 2){
        col_fun2 <- colorRamp2(c(min(sec_devs, na.rm = T), 0, 
                                 max(sec_devs, na.rm = T)), 
                               c("blue", "white", "orange"))
      } else{
        col_fun2 <- setNames("white", min(sec_devs, na.rm = T))
      }

      rowAnnot <- rowAnnotation(
        LRT_score = row_annots$LRT_score,
        second_derivative = row_annots$second_derivative,
        col = list(LRT_score = col_fun1,
                   second_derivative = col_fun2)
        )
    } else if (
      any(
        "L" %in% colnames(row_annots)
        # grepl("lrt_pass", colnames(row_annots))
        )) {
      if(split_cons){
         annots_cols <-
           list(conserved_tfbs = c("FALSE" = "grey","TRUE" = "black"))
         row_annots_cons <- row_annots %>% 
           mutate(conserved_tfbs = 
                    ifelse(conserved_tfbs, "Conserved", "Unconserved")) %>% 
           pull(conserved_tfbs)
         row_annots <- row_annots %>% select(-conserved_tfbs)
      } else{annots_cols <- list()}
    
    annot_col_list <- vector(mode = "list", length = 5) 
    names(annot_col_list) <- c("L", "C", "W", "P", 
                               paste0(main_type2, "-", main_type))
    for(annot in names(annot_col_list)){
      annot_col_list[[annot]] <- 
        c("FALSE" = "white","TRUE" = "black")
    }
    rowAnnot <- rowAnnotation(
      df = row_annots,
      col = 
        # c(annots_cols, 
          # list(
          #    lrt_pass = c("FALSE" = "#ada8e7","TRUE" = "#350cb0"),
          #    corr_pass = c("FALSE" = "#a5b8b7",
          #                  "TRUE" = "#204e4c"),
          #    wilc_pass = c("FALSE" = "#d6c6c4",
          #                  "TRUE" = "#79443B"))
        annot_col_list,
      show_legend = FALSE
      )
    } 
  } else {rowAnnot <- NULL}
  
  
  # for heatmap, replace NA with 0
  in_mat[is.na(in_mat)] <- 0
  # if (all(in_mat != 0)) {
  row_ttl <- character(0)
  if (split_cons) {
    r_split <- row_annots_cons
    if(length(unique(r_split)) == 1){
      r_split <- NULL
      row_ttl <- unique(r_split)
    }
  } else {r_split <- NULL}
  for (i in 1:2) {
    raster <- c(T, F)
    raster_by_magick <- c(TRUE, F)
    raster_suff <- c("_raster_", "_")
    
    for (j in 2) {
      row_order_method <- c("cluster", "means")[j]
      row_clust <- ifelse(row_order_method == "cluster", T, F)
      in_mat_mat <- in_mat %>% as.matrix()
      rownames(in_mat_mat) <- 1:nrow(in_mat_mat)
      if (nrow(in_mat_mat) > 1) {
        row_order <- 
          order(rowMeans(
            in_mat_mat[,(column_gap_pos[1]+1):column_gap_pos[2]]), 
              decreasing = T)
      } else{row_order <- 1}
      if(row_order_method == "cluster"){
        row_order <- NULL
      }
      ## import the motif logo, create an annotation with it
      logo_path_svg <- 
        list.files(file.path(proj_dir, "../data/best_motifs_061023/logos_SVG"),
                   pattern = paste0(TF_nm, "_.*logo.svg"), full.names = T)
      pwm_id <- gsub(".*_(PWM[0-9]+)_.*", "\\1", logo_path_svg)
      logo_svg <- rsvg_svg(logo_path_svg)
      logo_svg_read <- grImport2::readPicture(logo_path_svg)
      panel_fun = function(index, nm) {
        # grid.rect(height = unit(2, "cm"), gp = gpar(fill = NA, col = "black"))
        grid.picture(logo_svg_read, expansion = 0,
                     height = unit(2, "cm"), distort = T)
      }
      logo_block_annot <- columnAnnotation(logo = anno_block(
        align_to = (column_gap_pos[1] + 1):(column_gap_pos[2]),
        panel_fun = panel_fun,
        height = unit(1.2, "cm")
      ))
      
      ht_opt$COLUMN_ANNO_PADDING = unit(0.45, "cm")
      pdf(gsub(".pdf", 
               paste0(raster_suff[i], "_",
                      row_order_method, ".pdf"), outfile), paper = "letter"
          # height = 6
          )
      ht <- pheatmap(in_mat_mat,
               main = mat_nm,
               name = "phyloP",
               # annotation_legend = F,
               cluster_cols = F,
               use_raster = raster[i],
               raster_by_magick = raster_by_magick[i],
               raster_device = "CairoPNG",
               breaks = breaks,
               row_split = r_split,
               row_gap = unit(3, "mm"),
               cluster_row_slices	= F,
               cluster_rows = row_clust,
               row_order = row_order,
               show_rownames = F,
               show_colnames = F,
               # legend_breaks = c(-10, -1.3, -q_01, q_01, 1.3, 10),
               legend_breaks = lg_breaks,
               gaps_col = column_gap_pos,
               color = colorRampPalette(brewer.pal(
                 n = 7,
                 name = "RdBu"
               ))(20),
               right_annotation = rowAnnot,
               top_annotation = logo_block_annot,
               # run_draw = T,
               na_col = "black",
               # column_title = mat_nm,
               column_title_side = "bottom",
               row_title = row_ttl
      )
      draw(ht, padding = unit(c(2, 2, 15, 2), "mm"))
      decorate_annotation("logo", {
        grid.text(pwm_id, y = unit(0.2, "npc"), just = "right", x = unit(1, "npc"))
      })
      dev.off()
    }
  }  # }
}
