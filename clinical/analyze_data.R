# Required Libraries ====
# install PVBcorrect package from https://github.com/wnarifin/PVBcorrect
# ^^^^
# install.packages("boot", "mice")
# install.packages("devtools")
# devtools::install_github("wnarifin/PVBcorrect")
library(PVBcorrect)

# Load data ====

# 1 Hepatic Scintigraphy Test
data_hepatic = read.csv("hepatic.csv")
str(data_hepatic)
sapply(data_hepatic, unique)
view_table(data_hepatic, test = "test", disease = "disease", 
           show_unverified = T, show_total = T)

# 2 Diaphanography Test
data_diapha = read.csv("diaphanography.csv")
str(data_diapha)
sapply(data_diapha, unique)
view_table(data_diapha, test = "test", disease = "disease", 
           show_unverified = T, show_total = T)

# CCA
out_hepatic = acc_cca(data_hepatic, "test", "disease", ci = TRUE)
out_diapha = acc_cca(data_diapha, "test", "disease", ci = TRUE)
out_hepatic1 = round(out_hepatic$acc_results[1, c(1,3:4)],3)
out_hepatic2 = round(out_hepatic$acc_results[2, c(1,3:4)],3)
txt_hepatic1 = paste0(out_hepatic1[1], " (", out_hepatic1[2], ", ", out_hepatic1[3], ")")
txt_hepatic2 = paste0(out_hepatic2[1], " (", out_hepatic2[2], ", ", out_hepatic2[3], ")")
out_diapha1 = round(out_diapha$acc_results[1, c(1,3:4)],3)
out_diapha2 = round(out_diapha$acc_results[2, c(1,3:4)],3)
txt_diapha1 = paste0(out_diapha1[1], " (", out_diapha1[2], ", ", out_diapha1[3], ")")
txt_diapha2 = paste0(out_diapha2[1], " (", out_diapha2[2], ", ", out_diapha2[3], ")")
txt_cca = tibble::tibble(Methods = "CCA", Sn1 = txt_hepatic1, Sp1 = txt_hepatic2, Sn2 = txt_diapha1, Sp2 = txt_diapha2)

# BG
out_hepatic = acc_bg(data_hepatic, "test", "disease", ci = TRUE)
out_diapha = acc_bg(data_diapha, "test", "disease", ci = TRUE)
out_hepatic1 = round(out_hepatic$acc_results[1, c(1,3:4)],3)
out_hepatic2 = round(out_hepatic$acc_results[2, c(1,3:4)],3)
txt_hepatic1 = paste0(out_hepatic1[1], " (", out_hepatic1[2], ", ", out_hepatic1[3], ")")
txt_hepatic2 = paste0(out_hepatic2[1], " (", out_hepatic2[2], ", ", out_hepatic2[3], ")")
out_diapha1 = round(out_diapha$acc_results[1, c(1,3:4)],3)
out_diapha2 = round(out_diapha$acc_results[2, c(1,3:4)],3)
txt_diapha1 = paste0(out_diapha1[1], " (", out_diapha1[2], ", ", out_diapha1[3], ")")
txt_diapha2 = paste0(out_diapha2[1], " (", out_diapha2[2], ", ", out_diapha2[3], ")")
txt_bg = tibble::tibble(Methods = "BG", Sn1 = txt_hepatic1, Sp1 = txt_hepatic2, Sn2 = txt_diapha1, Sp2 = txt_diapha2)

# MI
out_hepatic = acc_mi(data_hepatic, "test", "disease", ci = TRUE, seednum = 3209673)
out_diapha = acc_mi(data_diapha, "test", "disease", ci = TRUE, seednum = 3209673)
out_hepatic1 = round(out_hepatic$acc_results[1, c(1,3:4)],3)
out_hepatic2 = round(out_hepatic$acc_results[2, c(1,3:4)],3)
txt_hepatic1 = paste0(out_hepatic1[1], " (", out_hepatic1[2], ", ", out_hepatic1[3], ")")
txt_hepatic2 = paste0(out_hepatic2[1], " (", out_hepatic2[2], ", ", out_hepatic2[3], ")")
out_diapha1 = round(out_diapha$acc_results[1, c(1,3:4)],3)
out_diapha2 = round(out_diapha$acc_results[2, c(1,3:4)],3)
txt_diapha1 = paste0(out_diapha1[1], " (", out_diapha1[2], ", ", out_diapha1[3], ")")
txt_diapha2 = paste0(out_diapha2[1], " (", out_diapha2[2], ", ", out_diapha2[3], ")")
txt_mi = tibble::tibble(Methods = "MI", Sn1 = txt_hepatic1, Sp1 = txt_hepatic2, Sn2 = txt_diapha1, Sp2 = txt_diapha2)

# IPW
out_hepatic = acc_ipw(data_hepatic, "test", "disease", ci = TRUE, ci_type = "perc", seednum = 3209673, R = 1000)
out_diapha = acc_ipw(data_diapha, "test", "disease", ci = TRUE, ci_type = "perc", seednum = 3209673, R = 1000)
out_hepatic1 = round(out_hepatic$acc_results[1, c(1,3:4)],3)
out_hepatic2 = round(out_hepatic$acc_results[2, c(1,3:4)],3)
txt_hepatic1 = paste0(out_hepatic1[1], " (", out_hepatic1[2], ", ", out_hepatic1[3], ")")
txt_hepatic2 = paste0(out_hepatic2[1], " (", out_hepatic2[2], ", ", out_hepatic2[3], ")")
out_diapha1 = round(out_diapha$acc_results[1, c(1,3:4)],3)
out_diapha2 = round(out_diapha$acc_results[2, c(1,3:4)],3)
txt_diapha1 = paste0(out_diapha1[1], " (", out_diapha1[2], ", ", out_diapha1[3], ")")
txt_diapha2 = paste0(out_diapha2[1], " (", out_diapha2[2], ", ", out_diapha2[3], ")")
txt_ipw = tibble::tibble(Methods = "IPWE", Sn1 = txt_hepatic1, Sp1 = txt_hepatic2, Sn2 = txt_diapha1, Sp2 = txt_diapha2)

# IPB
out_hepatic = acc_ipb(data_hepatic, "test", "disease", ci = TRUE, ci_perc = T, seed = 3209673, b = 1000)
out_diapha = acc_ipb(data_diapha, "test", "disease", ci = TRUE, ci_perc = T, seed = 3209673, b = 1000)
out_hepatic1 = round(out_hepatic$acc_results[1, c(1,3:4)],3)
out_hepatic2 = round(out_hepatic$acc_results[2, c(1,3:4)],3)
txt_hepatic1 = paste0(out_hepatic1[1], " (", out_hepatic1[2], ", ", out_hepatic1[3], ")")
txt_hepatic2 = paste0(out_hepatic2[1], " (", out_hepatic2[2], ", ", out_hepatic2[3], ")")
out_diapha1 = round(out_diapha$acc_results[1, c(1,3:4)],3)
out_diapha2 = round(out_diapha$acc_results[2, c(1,3:4)],3)
txt_diapha1 = paste0(out_diapha1[1], " (", out_diapha1[2], ", ", out_diapha1[3], ")")
txt_diapha2 = paste0(out_diapha2[1], " (", out_diapha2[2], ", ", out_diapha2[3], ")")
txt_ipb = tibble::tibble(Methods = "IPB", Sn1 = txt_hepatic1, Sp1 = txt_hepatic2, Sn2 = txt_diapha1, Sp2 = txt_diapha2)

# SIPW
start_t = proc.time()
out_hepatic = acc_sipw(data_hepatic, "test", "disease", ci = TRUE, ci_type = "perc", seednum = 3209673, b = 1000, R = 1000)
time_sipw_hepatic = calc_time_hms(proc.time() - start_t)
start_t = proc.time()
out_diapha = acc_sipw(data_diapha, "test", "disease", ci = TRUE, ci_type = "perc", seednum = 3209673, b = 1000, R = 1000)
time_sipw_diapha = calc_time_hms(proc.time() - start_t)
out_hepatic1 = round(out_hepatic$acc_results[1, c(1,3:4)],3)
out_hepatic2 = round(out_hepatic$acc_results[2, c(1,3:4)],3)
txt_hepatic1 = paste0(out_hepatic1[1], " (", out_hepatic1[2], ", ", out_hepatic1[3], ")")
txt_hepatic2 = paste0(out_hepatic2[1], " (", out_hepatic2[2], ", ", out_hepatic2[3], ")")
out_diapha1 = round(out_diapha$acc_results[1, c(1,3:4)],3)
out_diapha2 = round(out_diapha$acc_results[2, c(1,3:4)],3)
txt_diapha1 = paste0(out_diapha1[1], " (", out_diapha1[2], ", ", out_diapha1[3], ")")
txt_diapha2 = paste0(out_diapha2[1], " (", out_diapha2[2], ", ", out_diapha2[3], ")")
txt_sipw = tibble::tibble(Methods = "SIPW", Sn1 = txt_hepatic1, Sp1 = txt_hepatic2, Sn2 = txt_diapha1, Sp2 = txt_diapha2)

# SIPW-B
start_t = proc.time()
out_hepatic = acc_sipwb(data_hepatic, "test", "disease", ci = TRUE, ci_type = "perc", seednum = 3209673, b = 1000, R = 1000)
time_sipwb_hepatic = calc_time_hms(proc.time() - start_t)
start_t = proc.time()
out_diapha = acc_sipwb(data_diapha, "test", "disease", ci = TRUE, ci_type = "perc", seednum = 3209673, b = 1000, R = 1000)
time_sipwb_diapha = calc_time_hms(proc.time() - start_t)
out_hepatic1 = round(out_hepatic$acc_results[1, c(1,3:4)],3)
out_hepatic2 = round(out_hepatic$acc_results[2, c(1,3:4)],3)
txt_hepatic1 = paste0(out_hepatic1[1], " (", out_hepatic1[2], ", ", out_hepatic1[3], ")")
txt_hepatic2 = paste0(out_hepatic2[1], " (", out_hepatic2[2], ", ", out_hepatic2[3], ")")
out_diapha1 = round(out_diapha$acc_results[1, c(1,3:4)],3)
out_diapha2 = round(out_diapha$acc_results[2, c(1,3:4)],3)
txt_diapha1 = paste0(out_diapha1[1], " (", out_diapha1[2], ", ", out_diapha1[3], ")")
txt_diapha2 = paste0(out_diapha2[1], " (", out_diapha2[2], ", ", out_diapha2[3], ")")
txt_sipwb = tibble::tibble(Methods = "SIPW-B", Sn1 = txt_hepatic1, Sp1 = txt_hepatic2, Sn2 = txt_diapha1, Sp2 = txt_diapha2)

# Combine
txt_all = rbind(txt_cca,
                txt_bg,
                txt_ipw,
                txt_mi,
                txt_ipb,
                txt_sipw,
                txt_sipwb)
txt_all
saveRDS(txt_all, file = "all_results.Rdata")
cat(knitr::kable(txt_all, format = "simple"), file = "snsp_all_simple_m3.txt")
cat(knitr::kable(txt_all, format = "latex"), file = "snsp_all_latex_m3.txt")
