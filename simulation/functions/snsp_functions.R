# Contains methods to get Sn Sp point estimates

# Libraries ====
library(mice)

# Baselines ====

# CCA
snsp_cca = function(data) {
  snsp_values = c(0, 0)
  tbl = table(data$test, data$disease)
  snsp_values[1] = tbl[2,2]/sum(tbl[,2])  # Sn
  snsp_values[2] = tbl[1,1]/sum(tbl[,1])  # Sp
  return(snsp_values)
}
snsp = snsp_cca  # alt. known as snsp

# EBG
snsp_ebg = function(data, covariate = NULL) {
  snsp_values = c(0, 0)
  # data = complete case only for modeling, automatic by glm
  if (is.null(covariate)) {
    fit = glm(disease ~ test, data = data, family = "binomial")
  }
  else {
    z = reformulate(c("test", covariate), "disease")
    fit = glm(z, data = data, family = "binomial")
  }
  
  preds = predict(fit, data, type = "response")  # Predicted for complete & incomplete data
  snsp_values[1] = sum(data$test*preds) / sum(preds)  # P(T=1|D=1)
  snsp_values[2] = sum((1-data$test)*(1-preds)) / sum(1-preds)  # P(T=0|D=0)
  
  return(snsp_values)
}

# EBG, saturated model
snsp_ebg_saturated = function(data, covariate = NULL) {
  snsp_values = c(0, 0)
  # data = complete case only for modeling, automatic by glm
  if (is.null(covariate)) {
    fit = glm(disease ~ test, data = data, family = "binomial")
  }
  else {
    z = reformulate(paste(c("test", covariate), collapse = " * "), "disease")
    fit = glm(z, data = data, family = "binomial")
  }
  
  preds = predict(fit, data, type = "response")  # Predicted for complete & incomplete data
  snsp_values[1] = sum(data$test*preds) / sum(preds)  # P(T=1|D=1)
  snsp_values[2] = sum((1-data$test)*(1-preds)) / sum(1-preds)  # P(T=0|D=0)
  
  return(snsp_values)
}

# MI based ====

# MI, LogReg
# covariate not added yet
snsp_mi = function(data, m = 5, method = "logreg", seed = NA) {  # when using mice, set NA, else use NULL
  snsp_values = c(0, 0)
  # data = complete case only for modeling, automatic by glm
  data1 = data[, c("test", "disease")]
  data1$disease = as.factor(data1$disease)  # D as factor
  
  data_mi = mice(data1, m = m, method = method, seed = seed, print = F)  # MIDS class
  data_mi_data = complete(data_mi, "all")  # imputed data
  
  snsp_mi_imp = t(sapply(data_mi_data, snsp))
  # return(snsp_mi_imp)  # for debug purpose
  
  snsp_values = apply(snsp_mi_imp, 2, mean)
  return(snsp_values)
}

# PS based ====

# PS
# advantage:
# no need to obtain D|T i.e. PVs, just need V|T then
# get T|D directly

# IPW
snsp_ipw = function(data, covariate = NULL) {
  snsp_values = c(0, 0)
  
  # gen ps/pi
  # P(V = 1 | T = t)
  if (is.null(covariate)) {
    model_ps = glm(verify ~ test, data = data, family = "binomial")
  }
  else {
    z = reformulate(c("test", covariate), "verify")
    model_ps = glm(z, data = data, family = "binomial")
  }
  data$ps = predict(model_ps, type = "response")
  
  # recode NA to -1
  data[is.na(data$disease), "disease"] = -1  # so as data$verify*data$disease == 0 for NA
  
  snsp_values[1] = sum((data$test * data$verify * data$disease) / data$ps) /
    sum((data$verify * data$disease) / data$ps)
  snsp_values[2] = sum(((1 - data$test) * data$verify * (1 - data$disease)) / data$ps) /
    sum((data$verify * (1 - data$disease)) / data$ps)
  
  return(snsp_values)
}

# IPW Saturated
snsp_ipw_saturated = function(data, covariate = NULL) {
  snsp_values = c(0, 0)
  
  # gen ps/pi
  # P(V = 1 | T = t)
  if (is.null(covariate)) {
    model_ps = glm(verify ~ test, data = data, family = "binomial")
  }
  else {
    z = reformulate(paste(c("test", covariate), collapse = " * "), "verify")
    model_ps = glm(z, data = data, family = "binomial")
  }
  data$ps = predict(model_ps, type = "response")
  
  # recode NA to -1
  data[is.na(data$disease), "disease"] = -1  # so as data$verify*data$disease == 0 for NA
  
  snsp_values[1] = sum((data$test * data$verify * data$disease) / data$ps) /
    sum((data$verify * data$disease) / data$ps)
  snsp_values[2] = sum(((1 - data$test) * data$verify * (1 - data$disease)) / data$ps) /
    sum((data$verify * (1 - data$disease)) / data$ps)
  
  return(snsp_values)
}

# Bootstrap + resampling ====

# Bootstrap resampling 8, Mod-IPB
# based on IPB Nahorniak 2015
# utilizes weight to get resampling probability
# modified to resample to original size
snsp_resample8 = function(data, covariate = NULL, option = 1, interaction = FALSE,
                          m = 1, seed = NULL, return_data = FALSE) {
  # data = original data
  N_data = nrow(data)
  
  # gen ps
  # P(V = 1 | T = t)
  if (is.null(covariate)) {
    model_ps = glm(verify ~ test, data = data, family = "binomial")
  }
  else {
    if (interaction == TRUE) {
      z = reformulate(paste(c("test", covariate), collapse = " * "), "verify")
      model_ps = glm(z, data = data, family = "binomial")
    }
    else {
      z = reformulate(c("test", covariate), "verify")
      model_ps = glm(z, data = data, family = "binomial")
    }
  }
  # save ps to original data
  data$ps = predict(model_ps, type = "response")
  # transform ps to weight
  if (option == 1) {
    # IPW weight
    data$ps_w = ifelse(data$verify == 1, 1/data$ps, 1/(1-data$ps))
  }
  if (option == 2) {
    # W_h weight, Krautenbacher 2017
    data$ps_w = ifelse(data$verify == 1, max(unique(data$ps)) / data$ps,
                       max(unique((1-data$ps))) / (1-data$ps))
  }
  
  # data1 = complete case
  data1 = na.omit(data)
  n_data = nrow(data1)
  data1$p_ipb = data1$ps_w / sum(data1$ps_w)
  
  # resample
  data_resample_list = NULL
  set.seed(seed)
  i = 1
  while (i < m + 1) {
    # resample
    sampled = sample(1:n_data, N_data, TRUE, prob = data1$p_ipb)  # bootstrap w weight
    data_resample = data1[sampled, ]
    # exclude invalid samples: get the dimension sum, should be 4 for 2x2 epid table
    sum_tbl = sum(dim(table(data_resample$test, data_resample$disease)))
    # exclude invalid sample with sum < 4
    if (sum_tbl < 4) {
      # delete invalid sample
      data_resample_list[i] = NULL
      i = i
    } else {
      # save valid sample in list
      data_resample_list[i] = list(data_resample)
      i = i + 1
    }
  }
  
  snsp_resample_list = t(sapply(data_resample_list, snsp))
  
  snsp_values = apply(snsp_resample_list, 2, mean)
  
  if (return_data == TRUE) {
    return(list(data = data_resample_list, snsp = snsp_values))
  }
  else {
    return(snsp_values)
  }
}
snsp_sipw_resample = snsp_resample8  # scaled inverse probability weighting

# Bootstrap resampling 9
# based on IPB Nahorniak 2015
# utilizes weight to get resampling probability
# resample up to size of complete case, original IPB implementation
snsp_resample9 = function(data, covariate = NULL, option = 1, interaction = FALSE,
                          m = 1, seed = NULL, return_data = FALSE) {
  # data = original data

  # gen ps
  # P(V = 1 | T = t)
  if (is.null(covariate)) {
    model_ps = glm(verify ~ test, data = data, family = "binomial")
  }
  else {
    if (interaction == TRUE) {
      z = reformulate(paste(c("test", covariate), collapse = " * "), "verify")
      model_ps = glm(z, data = data, family = "binomial")
    }
    else {
      z = reformulate(c("test", covariate), "verify")
      model_ps = glm(z, data = data, family = "binomial")
    }
  }
  # save ps to original data
  data$ps = predict(model_ps, type = "response")
  # transform ps to weight
  if (option == 1) {
    # IPW weight
    data$ps_w = ifelse(data$verify == 1, 1/data$ps, 1/(1-data$ps))
  }
  if (option == 2) {
    # W_h weight, Krautenbacher 2017
    data$ps_w = ifelse(data$verify == 1, max(unique(data$ps)) / data$ps,
                       max(unique((1-data$ps))) / (1-data$ps))
  }
  
  # data1 = complete case
  data1 = na.omit(data)
  n_data = nrow(data1)
  data1$p_ipb = data1$ps_w / sum(data1$ps_w)
  
  # resample
  data_resample_list = NULL
  set.seed(seed)
  i = 1
  while (i < m + 1) {
    # resample
    sampled = sample(1:n_data, n_data, TRUE, prob = data1$p_ipb)  # bootstrap w weight
    data_resample = data1[sampled, ]
    # exclude invalid samples: get the dimension sum, should be 4 for 2x2 epid table
    sum_tbl = sum(dim(table(data_resample$test, data_resample$disease)))
    # exclude invalid sample with sum < 4
    if (sum_tbl < 4) {
      # delete invalid sample
      data_resample_list[i] = NULL
      i = i
    } else {
      # save valid sample in list
      data_resample_list[i] = list(data_resample)
      i = i + 1
    }
  }
  
  snsp_resample_list = t(sapply(data_resample_list, snsp))
  
  snsp_values = apply(snsp_resample_list, 2, mean)
  
  if (return_data == TRUE) {
    return(list(data = data_resample_list, snsp = snsp_values))
  }
  else {
    return(snsp_values)
  }
}
snsp_ipb = snsp_resample9

# Mod-IPB with higher prob for D=1, i.e. full sized sample
# auto multip, ratio of D=0 (larger group) to D=1 (smaller group)
# with more refined multip, almost equal d1:d0 ratio
# original sample size
snsp_resample15 = function(data, covariate = NULL, option = 1, interaction = FALSE,
                           m = 1, seed = NULL, rel_size = 1, return_data = FALSE) {
  # data = original data
  # rel_size = relative size of d0 to d1 in case-ctrl, default 1 i.e. 1:1
  N_data = nrow(data)
  
  # gen ps
  # P(V = 1 | T = t)
  if (is.null(covariate)) {
    model_ps = glm(verify ~ test, data = data, family = "binomial")
  }
  else {
    if (interaction == TRUE) {
      z = reformulate(paste(c("test", covariate), collapse = " * "), "verify")
      model_ps = glm(z, data = data, family = "binomial")
    }
    else {
      z = reformulate(c("test", covariate), "verify")
      model_ps = glm(z, data = data, family = "binomial")
    }
  }
  # save ps to original data
  data$ps = predict(model_ps, type = "response")
  # transform ps to weight
  if (option == 1) {
    # IPW weight
    data$ps_w = ifelse(data$verify == 1, 1/data$ps, 1/(1-data$ps))
  }
  if (option == 2) {
    # W_h weight, Krautenbacher 2017
    data$ps_w = ifelse(data$verify == 1, max(unique(data$ps)) / data$ps,
                       max(unique((1-data$ps))) / (1-data$ps))
  }
  
  # data1 = complete case
  data1 = na.omit(data)
  n_data = nrow(data1)
  n_d = table(data1$disease)
  multip = n_d[1] / n_d[2]  # d0:d1 ratio to inflate d1
  data1$ps_w_d1 = data1$ps_w
  data1$ps_w_d1[data1$disease == 1] = data1$ps_w_d1[data1$disease == 1] * multip
  k = sum(data1$ps_w_d1[data1$disease == 0])/  # correction factor
    sum(data1$ps_w_d1[data1$disease == 1])     # to allow almost equal d1:d0 ratio
  data1$ps_w_d1[data1$disease == 1] = data1$ps_w_d1[data1$disease == 1] * (k / rel_size)
  data1$p_ipb_d1 = data1$ps_w_d1 / sum(data1$ps_w_d1)
  
  # resample
  data_resample_list = NULL
  set.seed(seed)
  i = 1
  while (i < m + 1) {
    # resample
    sampled = sample(1:n_data, N_data, TRUE, prob = data1$p_ipb_d1)  # bootstrap w weight
    data_resample = data1[sampled, ]
    # exclude invalid samples: get the dimension sum, should be 4 for 2x2 epid table
    sum_tbl = sum(dim(table(data_resample$test, data_resample$disease)))
    # exclude invalid sample with sum < 4
    if (sum_tbl < 4) {
      # delete invalid sample
      data_resample_list[i] = NULL
      i = i
    } else {
      # save valid sample in list
      data_resample_list[i] = list(data_resample)
      i = i + 1
    }
  }
  
  snsp_resample_list = t(sapply(data_resample_list, snsp))
  
  snsp_values = apply(snsp_resample_list, 2, mean)
  
  if (return_data == TRUE) {
    return(list(data = data_resample_list, snsp = snsp_values))
  }
  else {
    return(snsp_values)
  }
}
snsp_sipw_resample_bal = snsp_resample15
