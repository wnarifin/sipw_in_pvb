# Generic experiment template for ease of setup
# create a folder named exp#, then set as working dir
# use the file to run the experiment

# Functions ====

## Load functions ====
source("../../functions/simulated_functions.R")
source("../../functions/snsp_functions.R")
source("../../functions/misc_functions.R")
source("../../exp3_functions.R")

## Load settings ====
source("../../functions/simulated_inits.R")

## This folder settings ====
source("par.R")  # do not load this for template_add, only load .Rdata
exp_num = stringr::str_split(getwd(), "/", simplify = T)
exp_num = exp_num[,length(exp_num)]

# p, Sn, Sp
p_pop = p_pop_list[p_pop]          # p
Sn_pop = SnSp_pop_df[SnSp_pop, 1]  # Sn
Sp_pop = SnSp_pop_df[SnSp_pop, 2]  # Sp
# mar, reasonable P(V=1|T=t)
p_sel_t1 = p_sel_t_df[p_sel_t, 1]  # P(V=1|T=1)
p_sel_t0 = p_sel_t_df[p_sel_t, 2]  # P(V=1|T=0)
# mnar
p_sel_d1 = p_sel_d_df[p_sel_d, 1]  # P(V=1|D=1)
p_sel_d0 = p_sel_d_df[p_sel_d, 2]  # P(V=1|D=0)

# Set Ps for simple 2x2 table
# multinomial distribution
# cell Ps calculated based on Bayesian rules
# p_td
# p_+1 = p, p_+0 = 1-p
# p_1+ = t+, p_0+ = t-
p11 = Sn_pop * p_pop  # s1: P(T=1, D=1) = P(T=1|D=1)P(D=1)
p01 = (1-Sn_pop) * p_pop # s0: P(T=0, D=1) = P(T=0|D=1)P(D=1)
p10 = (1-Sp_pop) * (1-p_pop)  # r1: P(T=1, D=0) = P(T=1|D=0)P(D=0)
p00 = Sp_pop * (1-p_pop)  # r0: P(T=0, D=0) = P(T=0|D=0)P(D=0)

# Define data generating distribution
# cells by col, i.e. by D=d
def = defData(varname = "td", formula = "p11;p01;p10;p00", dist = "categorical")
exp_settings = rbind(p_pop = p_pop, Sn_pop = Sn_pop, Sp_pop = Sp_pop,
                    p_sel_t1 = p_sel_t1, p_sel_t0 = p_sel_t0,
                    p_sel_d1 = p_sel_d1, p_sel_d0 = p_sel_d0,
                    p11 = p11, p01 = p01, p10 = p10, p00 = p00)
colnames(exp_settings) = "params"
capture.output(cat("This experiment folder settings:\n\n"),
               exp_settings,
               file = paste0(exp_num, "_settings.txt"))

# Gen m samples for each n size: PVB data ====

## gen complete data ====
# B samples & n sizes
multip = 25  # gen data 25 times larger, invalid rate around 96%!
data_list = gen_data_list(def, n, B * multip, seednum)

## gen data with pvb, mar ====
# B samples & n sizes
data_mar_list = gen_data_mar_list(data_list, n, B * multip, seednum)
# = Exclude invalid samples =
# - Ensure 2x2 table -
# get the dimension sum, should be 4 for 2x2 epid table
# then, return boolean status, TRUE if 4
data_mar_list_dim_check = lapply(data_mar_list,
                                 function(x) lapply(x,
                                                    function(y) {
                                                      z = na.omit(y)
                                                      z = table(z$test, z$disease)
                                                      z = sum(dim(z))
                                                      z > 3}))  # TRUE if > 3
# exclude the ones with sum < 4
for (i in 1:length(data_mar_list)) {
  data_mar_list[[i]] = data_mar_list[[i]][unlist(data_mar_list_dim_check[[i]])]
}
# - Ensure no cell with count = 0 -
data_mar_list_zero_check = lapply(data_mar_list,
                                  function(x) lapply(x,
                                                     function(y) {
                                                       z = na.omit(y)
                                                       z = table(z$test, z$disease)
                                                       z = as.numeric(z)
                                                       !any(z < 1)}))  # TRUE if does not contain 0
# exclude the ones with any cell count < 1
for (i in 1:length(data_mar_list)) {
  data_mar_list[[i]] = data_mar_list[[i]][unlist(data_mar_list_zero_check[[i]])]
}
# = Inform current number of replicates, i.e. minus invalid samples =
# important here to then include new number of replicates info, i.e. new B
# bcs it changes after excluding invalid samples
# e.g. cur_B = length(data_mar_list[[i]])
cur_B_mar = sapply(data_mar_list, length); cur_B_mar
# for pop1, 6, 1, barely make it...
# = reassign data_mar_list to first 500 =
for (i in 1:length(data_mar_list)) {
  data_mar_list[[i]] = data_mar_list[[i]][1:B]
}
# this will throw error if length < B
# = implement same process to data_list for sample-to-sample comparison =
for (i in 1:length(data_list)) {
  data_list[[i]] = data_list[[i]][unlist(data_mar_list_dim_check[[i]])]
}
for (i in 1:length(data_list)) {
  data_list[[i]] = data_list[[i]][unlist(data_mar_list_zero_check[[i]])]
}
cur_B = sapply(data_list, length); cur_B
for (i in 1:length(data_list)) {
  data_list[[i]] = data_list[[i]][1:B]
}
new_B = sapply(data_list, length); new_B
# = inform the right number of replicates =
# important to check
new_B_mar = sapply(data_mar_list, length); new_B_mar
capture.output(cat("B_mar\n"), cur_B_mar, new_B_mar,
               cat("B\n"), cur_B, new_B,
               file = paste0(exp_num, "_cur_new_B_mar.txt"))

# get missing rate & complete n size
missing_mar_means = get_missing_rate_sim(data_mar_list); round(missing_mar_means, 3)
n_complete_means = n * (1 - as.numeric(missing_mar_means)); round(n_complete_means, 1)
missing_n_mar_means = cbind(missing_mar_means, n_complete_means)
colnames(missing_n_mar_means) = c("Missing rate", "n complete")
capture.output(missing_n_mar_means,
               file = paste0(exp_num, "_missing_n_mar_means.txt"))

# ## gen data with pvb, mnar ====
# # B samples & n sizes
# data_mnar_list = gen_data_mnar_list(data_list , n, B, seednum)
#
# # get missing rate
# missing_mnar_means = get_missing_rate_sim(data_mnar_list); missing_mnar_means
#
# # need to add validity check first before perform mar correction
# # may be better in a separate experiment

# Get estimates - MAR - ====

## --- Gold standard ---

## Unbiased, complete data ====
snsp_df = get_snsp(data_list, "snsp", n, new_B_mar)
write.csv(snsp_df, "snsp_df.csv", row.names = F)
snsp_txt = format_snsp(snsp_df, label = "Complete data, CCA estimates")
capture.output(snsp_txt, file = "snsp.txt")

## --- Comparison methods ---

## CCA ====
ptm0 = proc.time()
snsp_cca_df = get_snsp(data_mar_list, "snsp", n, new_B_mar)
write.csv(snsp_cca_df, "snsp_cca_df.csv", row.names = F)
snsp_cca_txt = format_snsp(snsp_cca_df, label = "PVB data, CCA estimates")
ptm1 = proc.time() - ptm0; ptm1 = calc_time_hms(ptm1["elapsed"])
capture.output(snsp_cca_txt, ptm1, file = "snsp_cca.txt")

## BG/EBG ====
ptm0 = proc.time()
snsp_ebg_df = get_snsp(data_mar_list, "snsp_ebg", n, new_B_mar)
write.csv(snsp_ebg_df, "snsp_ebg_df.csv", row.names = F)
snsp_ebg_txt = format_snsp(snsp_ebg_df, label = "PVB data, EBG estimates")
ptm1 = proc.time() - ptm0; ptm1 = calc_time_hms(ptm1["elapsed"])
capture.output(snsp_ebg_txt, ptm1, file = "snsp_ebg.txt")

## IPW ====
ptm0 = proc.time()
snsp_ipw_df = get_snsp(data_mar_list, "snsp_ipw", n, new_B_mar)
write.csv(snsp_ipw_df, "snsp_ipw_df.csv", row.names = F)
snsp_ipw_txt = format_snsp(snsp_ipw_df, label = "PVB data, IPW estimates")
ptm1 = proc.time() - ptm0; ptm1 = calc_time_hms(ptm1["elapsed"])
capture.output(snsp_ipw_txt, ptm1, file = "snsp_ipw.txt")

## MI, LogReg ====
ptm0 = proc.time()
snsp_mi_df = get_snsp(data_mar_list, "snsp_mi", n, new_B_mar,
                      mi_based = TRUE, m = m, seednum = seednum)
write.csv(snsp_mi_df, "snsp_mi_df.csv", row.names = F)
snsp_mi_txt = format_snsp(snsp_mi_df, label = "PVB data, MI LogReg estimates")
ptm1 = proc.time() - ptm0; ptm1 = calc_time_hms(ptm1["elapsed"])
capture.output(snsp_mi_txt, ptm1, file = "snsp_mi.txt")

## --- New methods ---

## Resampling 9 / IPB ====
# bootstrap resampling IPB up to size of complete case, original IPB implementation
ptm0 = proc.time()
snsp_ipb_df = get_snsp(data_mar_list, "snsp_ipb", n, new_B_mar,
                              sample_based = TRUE, m1 = b, seednum1 = seednum)
write.csv(snsp_ipb_df, "snsp_ipb_df.csv", row.names = F)
snsp_ipb_txt = format_snsp(snsp_ipb_df, label = "PVB data, IPB estimates")
ptm1 = proc.time() - ptm0; ptm1 = calc_time_hms(ptm1["elapsed"])
capture.output(snsp_ipb_txt, ptm1, file = "snsp_ipb.txt")

## Resampling 8 / SIPW ====
# bootstrap resampling IPB up to full sample size
ptm0 = proc.time()
snsp_resample8_df = get_snsp(data_mar_list, "snsp_resample8", n, new_B_mar,
                             sample_based = TRUE, m1 = b, seednum1 = seednum)
write.csv(snsp_resample8_df, "snsp_resample8_df.csv", row.names = F)
snsp_resample8_txt = format_snsp(snsp_resample8_df, label = "PVB data, Resampling 8 estimates")
ptm1 = proc.time() - ptm0; ptm1 = calc_time_hms(ptm1["elapsed"])
capture.output(snsp_resample8_txt, ptm1, file = "snsp_resample8.txt")

## Resampling 15 / SIPW-B ====
# bootstrap resampling IPB up to full sample size with adjusted upsampling
ptm0 = proc.time()
snsp_resample15_df = get_snsp(data_mar_list, "snsp_resample15", n, new_B_mar,
                              sample_based = TRUE, m1 = b, seednum1 = seednum)
write.csv(snsp_resample15_df, "snsp_resample15_df.csv", row.names = F)
snsp_resample15_txt = format_snsp(snsp_resample15_df, label = "PVB data, Resampling 15 estimates")
ptm1 = proc.time() - ptm0; ptm1 = calc_time_hms(ptm1["elapsed"])
capture.output(snsp_resample15_txt, ptm1, file = "snsp_resample15.txt")

# Save workspace ====
time_tag = format(Sys.time(), "%Y%m%d_%H%M%S")
save.image(file = paste0(exp_num, "_workspace_", time_tag, ".Rdata"))
save.image(file = "latest.Rdata")
# load("name.Rdata")  # replace with file name.Rdata

# Comparison table ====
pvb_methods_ = c("Complete", "CCA", "EBG", "IPW", 
                "MI",
                "IPB", "SIPW Resample", "SIPW Resample Bal")
snsp_dfs_ = rbind(snsp_df, snsp_cca_df, snsp_ebg_df, snsp_ipw_df,
                 snsp_mi_df,
                 snsp_ipb_df, snsp_resample8_df, snsp_resample15_df)
# above lines will change with _add
pvb_methods = pvb_methods_
snsp_dfs = cbind(Methods = rep(pvb_methods, each = length(n)), snsp_dfs_)
n200_index = seq(1, 3*length(pvb_methods), 3)
n500_index = seq(2, 3*length(pvb_methods), 3)
n1000_index = seq(3, 3*length(pvb_methods), 3)
n200 = snsp_dfs[n200_index, ]
n500 = snsp_dfs[n500_index, ]
n1000 = snsp_dfs[n1000_index, ]

comparison_tbl = list("n = 200" = n200,
                      "n = 500" = n500, # not reported in paper 3, no contrasting results
                      "n = 1000" = n1000)
capture.output(
  cat(paste0("Comparison between PVB methods (MAR). Sn = ", Sn_pop,
             ", Sp = ", Sp_pop, ", p = ", p_pop, ".\n\n",
             "Verification probabilities: P(V=1|T=1) = ", p_sel_t1,
             ", P(V=1|T=0) = ", p_sel_t0, ".")),
  knitr::kable(comparison_tbl, format = "simple", row.names = F),
  file = "snsp_all_tbl.txt")
capture.output(
  cat(paste0("Comparison between PVB methods (MAR). Sn = ", Sn_pop,
             ", Sp = ", Sp_pop, ", p = ", p_pop, ".\n\n",
             "Verification probabilities: P(V=1|T=1) = ", p_sel_t1,
             ", P(V=1|T=0) = ", p_sel_t0, ".")),
  knitr::kable(comparison_tbl, format = "latex", row.names = F),
  file = "snsp_all_tbl_latex.txt")
