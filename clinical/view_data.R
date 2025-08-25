# data only

# 1 Hepatic Scintigraphy Test
data_hepatic = read.csv("hepatic.csv")
str(data_hepatic)
sapply(data_hepatic, unique)
summary(data_hepatic)
with(data_hepatic, table(test, disease, useNA = "always")) |> addmargins() -> tbl_hepatic
tbl_hepatic
with(data_hepatic, table(test, disease, useNA = "no")) |> addmargins()
# formatted table
tbl_hepatic |> knitr::kable(format = "simple")
# missing %
1 - mean(data_hepatic$verified)

# 2 Diaphanography Test
data_diapha = read.csv("diaphanography.csv")
str(data_diapha)
sapply(data_diapha, unique)
summary(data_diapha)
with(data_diapha, table(test, disease, useNA = "always")) |> addmargins() -> tbl_diapha
tbl_diapha
with(data_diapha, table(test, disease, useNA = "no")) |> addmargins()
# formatted table
tbl_diapha |> knitr::kable(format = "simple")
# missing %
1 - mean(data_diapha$verified)  
