# exp3 functions

# gen complete data ====
# B samples & n sizes
gen_data_list = function(def = def, n = n, B = B, seednum = seednum) {
  data_list = vector("list", length(n))
  for (i in 1:length(n)) {
    set.seed(seednum)  # to reset at each n size iteration, so comparable with one step lower n
    # it also allows snsp_data_gen to be run independently per n size
    data_list_each_n = vector("list", B)  # reset to init for each loop
    for (j in 1:B) {
      data_list_each_n[j] = list(snsp_data_gen(n[i], def))
    }
    data_list[i] = list(data_list_each_n)
  }
  return(data_list)
}

# gen data with pvb, mar ====
# B samples & n sizes
gen_data_mar_list = function(data_list, n = n, B = B, seednum = seednum) {
  data_mar_list = vector("list", length(n))
  for (i in 1:length(n)) {
    set.seed(seednum)  # to reset at each n size iteration, so comparable with one step lower n
    # it also allows snsp_data_gen to be run independently per n size
    data_mar_list_each_n = vector("list", B)  # reset to init for each loop
    for (j in 1:B) {
      data_mar_list_each_n[j] = list(pvb_mar_gen(data_list[[i]][[j]], p1 = p_sel_t1, p0 = p_sel_t0))
    }
    data_mar_list[i] = list(data_mar_list_each_n)
  }
  return(data_mar_list)
}

# get missing rate ====
get_missing_rate_sim = function(data_list) {
  missing_list = lapply(data_list, function(x) {t(sapply(x, function(y) 1 - mean(y$verify)))})
  missing_means = as.matrix(sapply(missing_list, mean))
  colnames(missing_means) = "Missing"; rownames(missing_means) = n
  return(missing_means)
}
# 1 - P(V=1)
# MAR: P(V=1) = sum{P(V=1|T=t)P(T=t)}
# MNAR: # P(V=1) = ?, need to obtain by simulation
get_missing_rate_sim_detail = function(data_list) {
  missing_list = lapply(data_list, function(x) {t(sapply(x, function(y) 1 - mean(y$verify)))})
  return(missing_list)
}

# get snsp ====
get_snsp = function(data_list, snsp_fun = "snsp", n = n, new_B = new_B,
                    mi_based = FALSE, m = 5, seednum = NA,
                    sample_based = FALSE, m1 = 1, seednum1 = NULL,
                    em_based = FALSE, t_max = 500, cutoff = 0.0001) {
  # snsp
  if (mi_based == TRUE) {
    snsp_list = lapply(data_list, function(x) {t(sapply(x, snsp_fun, m = m, seed = seednum))})
  }
  if (sample_based == TRUE) {
    snsp_list = lapply(data_list, function(x) {t(sapply(x, snsp_fun, m = m1, seed = seednum1))})
  }
  if (em_based == TRUE) {
    snsp_list = lapply(data_list, function(x) {t(sapply(x, snsp_fun, t_max = t_max, cutoff = cutoff))})
  }
  else {
    snsp_list = lapply(data_list, function(x) {t(sapply(x, snsp_fun))})
  }
  for (i in 1:length(n)) {colnames(snsp_list[[i]]) = c("Sn", "Sp")}

  # mean
  snsp_means = t(sapply(snsp_list, colMeans))  # mean
  colnames(snsp_means) = c("Sn", "Sp"); rownames(snsp_means) = n

  # SE
  snsp_ses = t(sapply(snsp_list, function(x) apply(x, 2, sd)))  # SD/SE
  colnames(snsp_ses) = c("Sn", "Sp"); rownames(snsp_ses) = n

  # bias, MSE
  # snsp_pop_mat = matrix(rep(c(Sn_pop, Sp_pop), B), ncol = 2, byrow = T)  # pop values
  # snsp_rawdiffs = lapply(snsp_list, function(x) (x - snsp_pop_mat))  # raw diff, estimate - true value
  snsp_rawdiffs = lapply(snsp_list, function(x) {
    pop_mat = matrix(rep(c(Sn_pop, Sp_pop), length(unlist(x))/2), ncol = 2, byrow = T)
    x - pop_mat}
  )  # raw diff, estimate - true value
  snsp_biases = t(sapply(snsp_list, function(x) {colMeans(x) - c(Sn_pop, Sp_pop)}))  # Bias
  rownames(snsp_biases) = n
  snsp_mses = t(sapply(snsp_rawdiffs, function(x) colMeans(x^2)))  # MSE 
  # -- not reported in paper 3, focus on bias, SE
  rownames(snsp_mses) = n

  # tibble output
  snsp_df = tibble::tibble(n, B = new_B,
                           "Sn mean"=snsp_means[,"Sn"], "Sn bias"=snsp_biases[,"Sn"],
                           "Sn SE"=snsp_ses[,"Sn"], "Sn MSE"=snsp_mses[,"Sn"],
                           "Sp mean"=snsp_means[,"Sp"], "Sp bias"=snsp_biases[,"Sp"],
                           "Sp SE"=snsp_ses[,"Sp"], "Sp MSE"=snsp_mses[,"Sp"])
  return(snsp_df)  # returns a tibble
}

# create plot from kable ====
library(ggplot2)
comparison_plt = function(comparison_tbl_sub, shp = 5, lty = 3) {
  tbl = comparison_tbl_sub[,-c(2,3)]  # exclude n & B columns
  n = comparison_tbl_sub[1,2]  # get n from column 2

  df_sn = data.frame(grp = reorder(as.character(tbl$Methods), nrow(tbl):1), fit = tbl$`Sn mean`, se = tbl$`Sn SE`)
  k = ggplot(df_sn, aes(grp, fit, ymin = fit-se, ymax = fit+se)) + scale_y_continuous(breaks = seq(-0.2,1.2,0.1))
  k = k + geom_pointrange(shape=shp) + labs(y="Sn ± SE", x=paste0("Methods (n = ", n, ")")) + geom_hline(yintercept=Sn_pop,linetype=lty) + coord_flip() + theme_bw()

  df_sp = data.frame(grp = reorder(as.character(tbl$Methods), nrow(tbl):1), fit = tbl$`Sp mean`, se = tbl$`Sp SE`)
  l = ggplot(df_sp, aes(grp, fit, ymin = fit-se, ymax = fit+se)) + scale_y_continuous(breaks = seq(-0.2,1.2,0.1))
  l = l + geom_pointrange(shape=shp) + labs(y="Sp ± SE", x="") + geom_hline(yintercept=Sp_pop,linetype=lty) + coord_flip() + theme_bw()
  l = l + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())

  return(list(k,l))
}

