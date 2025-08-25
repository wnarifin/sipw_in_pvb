# A collection of misc functions used in project

#' View table
view_tbl = function(data, show_unverified = TRUE, show_total = TRUE) {
  # assign variables
  # genData return data.table and data.frame
  # so data[, "name] will not work
  # either data$name, or data[, name] or
  # name = "name"; data[, name] will work
  test = data$test
  disease = data$disease
  
  # option: show verified/not
  if (show_unverified == TRUE) {
    tbl = table(Test = 2 - test, Disease = 2 - disease, useNA = "ifany")  # rearrange for epidemiology view
    if (ncol(tbl) < 3) {
      tbl = cbind(tbl, na = c(0, 0))
    }  # whenever there are no missing disease status
    dimnames(tbl) = list(Test = c("yes","no"),
                         Disease = c("yes","no","unverified"))
  } else if (show_unverified == FALSE) {
    tbl = table(Test = 2 - test, Disease = 2 - disease)  # rearrange for epidemiology view
    rownames(tbl) = colnames(tbl) = c("yes","no")
  }
  
  if (show_total == TRUE) {
    tbl = addmargins(tbl)
    colnames(tbl)[ncol(tbl)] = "Total"
    rownames(tbl)[nrow(tbl)] = "Total"
  }
  
  return(tbl)
}

#' Odds ratio
odds_ratio = function(x1, x2) {
  tbl1 = table(2 - x1, 2 - x2,
               dnn = c(deparse(substitute(x1)), deparse(substitute(x2))))  # rearrange for epid view
  rownames(tbl1) = colnames(tbl1) = c("yes","no")
  # if (sum(tbl1 > 0) < 4) {tbl1 = tbl1 + .5}  # correction for 0 cell
  odds1 = prod(diag(tbl1))
  odds0 = prod(tbl1[upper.tri(tbl1)], tbl1[lower.tri(tbl1)])
  odds_ratio = odds1 / odds0
  return(list(Table = tbl1, OR = odds_ratio))
}

#' Get missing rate on data_list
get_missing_rate = function(data_list) {
  list = lapply(data_list, function(data) {1 - mean(data$verify)})
  matrix = as.matrix(list)
  colnames(matrix) = "Missing"; rownames(matrix) = sapply(data_list, nrow)
  return(matrix)
}

#' Get info on data_list
get_info = function(data_list) {
  lapply(data_list, function(data) {
  list(
    dim = dim(data),
    str = ls.str(data),
    desc = desc_cat(as.data.frame(data)[,-1])
  )})
}

# get Pr(D=1) indirectly
# 2.2 He & McDermott (2012)
# PPV * P(T=1) + (1 - NPV) * P(T = 0)
get_p_indirect = function(data, test, disease) {
  acc = PVBcorrect::acc_cca(data, test, disease, description = F)$acc_results
  ppv = acc[3, ]
  npv = acc[4, ]
  ppv * mean(data[, test] == 1, na.rm = T) + (1 - npv) * mean(data[, test] == 0, na.rm = T)
}

#' Descriptive for categorical variables
#' @description A convenient function to come out with n (%) side by side.
#' @param data Data frame as input.
#' @param percent Can select to display percent (TRUE) or proportion (FALSE).
#' @return Return a list of data frame(s)
desc_cat = function(data, percent = TRUE) {
  list = vector("list", ncol(data))
  for (i in 1:ncol(data)) {
    n = as.data.frame(table(data[i]))
    prop = as.data.frame(prop.table(table(data[i])))
    var_name = c(names(data[i]), rep("-", length(table(data[i]))-1))
    df = data.frame(Variable = var_name,
                    Label = n$Var1,
                    n = n$Freq)
    if(!percent) {df$Proportion = prop$Freq} else {df$Percent = prop$Freq*100}
    list[[i]] = df
  }
  names(list) = names(data)
  return(list)
}

#' Time in seconds to hours, minutes, seconds.
#'
#' @description Convert time in seconds from e.g. proc.time() function to hours, minutes, seconds.
#  This was based on https://stackoverflow.com/questions/6118922/convert-seconds-value-to-hours-minutes-seconds
#' @param time_in_secs A data frame, with at least "Test" and "Disease" variables.
#' @return Time in hours, minutes and seconds.
#' @details See examples in \code{\link{acc_em}} for use case.
#' @export
calc_time_hms = function(time_in_secs) {
  hms = data.frame(hours = floor(time_in_secs / 3600), minutes = floor(time_in_secs %% 3600 / 60), seconds = time_in_secs %% 60)
  return(hms)
}

