#' Summarize one-way anova and Tukey test results.
#' Gives the mean value, grouping label, probability, MSD value
#'
#' @param df A data.frame
#' @param Vari. A character vacter, column names of the variables
#' @param Factor The column name of the factor
#' @param DeciNum A numeric vector of decimal numbers, default is 2
#' @return A data.frame
#' @examples
#' Omitted
#' @export
Aov1TukeySummary <- function(df, Vari., Factor, DeciNum = rep(2, length(Vari.))){
  HSD_fun <- function(q,s,n_vector,n_group){
    return(q * sqrt(s * sum(1/n_vector)/n_group))
  }

  tmpdf <- mutate(df, Grp1 = df[,Factor])
  FinalTab <- matrix(ncol = length(Vari.)+1, nrow = tmpdf$Grp1 %>% unique() %>% length() %>% +2) %>% data.frame() %>% setNames(., c("Group",Vari.))
  FinalTab$Group <- unique(tmpdf$Grp1) %>% c(., "P", "MSD")
  for (y in Vari.){
    n <- which(Vari. == y)
    # Prob
    Prob <- aov(as.formula(paste(y,"~", "Grp1")), data = tmpdf) %>% summary() %>% .[[1]] %>% .[1,"Pr(>F)"]
    FinalTab[FinalTab$Group == "P",y] <- ifelse(Prob < 0.0001, "<0.0001", sprintf("%.4f",Prob))
    # MSD
    tmp_model <- aov(as.formula(paste(y,"~", "Grp1")), data = tmpdf)
    tmp_Tukey_result <- HSD.test(tmp_model,"Grp1", group = T)
    FinalTab[FinalTab$Group == "MSD", y] <- sprintf(HSD_fun(q = as.numeric(tmp_Tukey_result$parameters["StudentizedRange"]),
                                                            s = as.numeric(tmp_Tukey_result$statistics["MSerror"]),
                                                            n_vector = c(tmp_Tukey_result$means$r),
                                                            n_group =  tmp_Tukey_result$parameters$ntr), fmt = paste0("%.",DeciNum[n]+1,"f"))
    # mean value and group
    for (i in unique(tmpdf$Grp1)){
      FinalTab[FinalTab$Group == i, y] <-
        paste(sprintf(tmp_Tukey_result$groups[i,1], fmt = paste0("%.",DeciNum[n],"f")), tmp_Tukey_result$groups[i,"groups"])
    }
  }
  return(FinalTab)
}

#' Summarize two-way anova and Tukey test results
#' Gives the mean value, grouping label, probability including interactions, MSD value, Eta_square
#'
#' @param df A data.frame
#' @param Vari. A character vacter, column names of the variables
#' @param Factor1 The column name of the factor1
#' @param Factor2 The column name of the factor2
#' @param DeciNum A numeric vector of decimal numbers, default is 2
#' @return A data.frame
#' @examples
#' Omitted
#' @export
Aov2TukeySummary <- function(df, Vari., Factor1, Factor2, DeciNum = rep(2, length(Vari.))){
  HSD_fun <- function(q,s,n_vector,n_group){
    return(q * sqrt(s * sum(1/n_vector)/n_group))
  }

  tmpdf <- mutate(df, Grp1 = df[,Factor1], Grp2 = df[,Factor2])
  RowNm <- c(tmpdf$Grp1 %>% unique(), paste0("P_", Factor1), paste0("MSD_", Factor1),
             "Factor2", tmpdf$Grp2 %>% unique(), paste0("P_", Factor2), paste0("MSD_", Factor2),
             "Interaction", paste0("P_",Factor1,":",Factor2),
             "Eta2", paste0("Eta2_",Factor1),paste0("Eta2_",Factor2), paste0("Eta2_",Factor1,":",Factor2), "Eta2_Residuals")
  FinalTab <- matrix(ncol = length(Vari.)+1,
                     nrow = length(RowNm)) %>% data.frame() %>% setNames(., c("Factor1",Vari.))
  FinalTab$Factor1 <- RowNm
  FinalTab[FinalTab$Factor1 == "Factor2",-1] <-
    FinalTab[FinalTab$Factor1 == "Interaction",-1] <-
    FinalTab[FinalTab$Factor1 == "Eta2",-1] <-
    colnames(FinalTab)[-1]

  for (y in Vari.){
    n <- which(Vari. == y)
    # Prob
    Prob <- aov(as.formula(paste0(y,"~", Factor1, "*", Factor2)), data = tmpdf) %>% summary() %>% .[[1]] %>% as.data.frame()
    rownames(Prob) <- rownames(Prob) %>% ifelse(grepl(" ", .), gsub(" ","",.), .)
    FinalTab[FinalTab$Factor1 == paste0("P_", Factor1), y] <- ifelse(Prob[Factor1,"Pr(>F)"] < 0.0001, "<0.0001", sprintf("%.4f",Prob[Factor1,"Pr(>F)"]))
    FinalTab[FinalTab$Factor1 == paste0("P_", Factor2), y] <- ifelse(Prob[Factor2,"Pr(>F)"] < 0.0001, "<0.0001", sprintf("%.4f",Prob[Factor2,"Pr(>F)"]))
    FinalTab[FinalTab$Factor1 == paste0("P_", Factor1, ":",Factor2), y] <- ifelse(Prob[paste0(Factor1, ":",Factor2),"Pr(>F)"] < 0.0001, "<0.0001", sprintf("%.4f",Prob[paste0(Factor1, ":",Factor2),"Pr(>F)"]))

    # Eta2
    FinalTab[FinalTab$Factor1 == paste0("Eta2_", Factor1), y] <- (Prob[Factor1,"Sum Sq"] / sum(Prob[,"Sum Sq"])) %>% sprintf("%.4f",.)
    FinalTab[FinalTab$Factor1 == paste0("Eta2_", Factor2), y] <- (Prob[Factor2,"Sum Sq"] / sum(Prob[,"Sum Sq"])) %>% sprintf("%.4f",.)
    FinalTab[FinalTab$Factor1 == paste0("Eta2_",Factor1, ":",Factor2), y] <- (Prob[paste0(Factor1, ":",Factor2),"Sum Sq"] / sum(Prob[,"Sum Sq"])) %>% sprintf("%.4f",.)
    FinalTab[FinalTab$Factor1 == paste0("Eta2_Residuals"), y] <- (Prob["Residuals","Sum Sq"] / sum(Prob[,"Sum Sq"])) %>% sprintf("%.4f",.)
    # MSD and mean values
    tmp_model1 <- aov(as.formula(paste0(y,"~", Factor1, "*", Factor2)), data = tmpdf)
    tmp_Tukey_result <- HSD.test(tmp_model1,Factor1, group = T)
    FinalTab[FinalTab$Factor1 == paste0("MSD_",Factor1), y] <- sprintf(HSD_fun(q = as.numeric(tmp_Tukey_result$parameters["StudentizedRange"]),
                                                                               s = as.numeric(tmp_Tukey_result$statistics["MSerror"]),
                                                                               n_vector = c(tmp_Tukey_result$means$r),
                                                                               n_group =  tmp_Tukey_result$parameters$ntr), fmt = paste0("%.",DeciNum[n]+1,"f"))
    for (i in unique(tmpdf[,Factor1])){
      FinalTab[FinalTab$Factor1 == i, y] <-
        paste(sprintf(tmp_Tukey_result$groups[i,1], fmt = paste0("%.",DeciNum[n],"f")), tmp_Tukey_result$groups[i,"groups"])
    }

    tmp_model2 <- aov(as.formula(paste0(y,"~", Factor1, "*", Factor2)), data = tmpdf)
    tmp_Tukey_result <- HSD.test(tmp_model2,Factor2, group = T)
    FinalTab[FinalTab$Factor1 == paste0("MSD_",Factor2), y] <- sprintf(HSD_fun(q = as.numeric(tmp_Tukey_result$parameters["StudentizedRange"]),
                                                                               s = as.numeric(tmp_Tukey_result$statistics["MSerror"]),
                                                                               n_vector = c(tmp_Tukey_result$means$r),
                                                                               n_group =  tmp_Tukey_result$parameters$ntr), fmt = paste0("%.",DeciNum[n]+1,"f"))
    for (j in unique(tmpdf[,Factor2])){
      FinalTab[FinalTab$Factor1 == j, y] <-
        paste(sprintf(tmp_Tukey_result$groups[j,1], fmt = paste0("%.",DeciNum[n],"f")), tmp_Tukey_result$groups[j,"groups"])
    }
  }
  return(FinalTab)
}

#' Summarize two-way anova and Tukey test results
#' Gives the mean value, grouping label, probability including interactions, MSD value, Eta_square
#'
#' @param df A data.frame
#' @param Vari. A character vacter, column names of the variables
#' @param Factor1 The column name of the factor1
#' @param Factor2 The column name of the factor2
#' @param Factor3 The column name of the factor3
#' @param DeciNum A numeric vector of decimal numbers, default is 2
#' @return A data.frame
#' @examples
#' Omitted
#' @export
Aov3TukeySummary <- function(df, Vari., Factor1, Factor2, Factor3, DeciNum = rep(2, length(Vari.))){
  HSD_fun <- function(q,s,n_vector,n_group){
    return(q * sqrt(s * sum(1/n_vector)/n_group))
  }

  tmpdf <- mutate(df, Grp1 = df[,Factor1], Grp2 = df[,Factor2], Grp3 = df[,Factor3])
  RowNm <- c(tmpdf$Grp1 %>% unique(), paste0("P_", Factor1), paste0("MSD_", Factor1),
             "Factor2", tmpdf$Grp2 %>% unique(), paste0("P_", Factor2), paste0("MSD_", Factor2),
             "Factor3", tmpdf$Grp3 %>% unique(), paste0("P_", Factor3), paste0("MSD_", Factor3),
             "Interaction", paste0("P_",Factor1,":",Factor2), paste0("P_",Factor1,":",Factor3), paste0("P_",Factor2,":",Factor3), paste0("P_",Factor1,":",Factor2,":",Factor3),
             "Eta2", paste0("Eta2_",Factor1),paste0("Eta2_",Factor2), paste0("Eta2_",Factor3), paste0("Eta2_",Factor1,":",Factor2), paste0("Eta2_",Factor1,":",Factor3), paste0("Eta2_",Factor2,":",Factor3),paste0("Eta2_",Factor1,":",Factor2,":",Factor3), "Eta2_Residuals")
  FinalTab <- matrix(ncol = length(Vari.)+1,
                     nrow = length(RowNm)) %>% data.frame() %>% setNames(., c("Factor1",Vari.))
  FinalTab$Factor1 <- RowNm
  FinalTab[FinalTab$Factor1 == "Factor2",-1] <-
    FinalTab[FinalTab$Factor1 == "Factor3",-1] <-
    FinalTab[FinalTab$Factor1 == "Interaction",-1] <-
    FinalTab[FinalTab$Factor1 == "Eta2",-1] <-
    colnames(FinalTab)[-1]

  for (y in Vari.){
    n <- which(Vari. == y)
    # Prob
    Prob <- aov(as.formula(paste0(y,"~", Factor1, "*", Factor2, "*", Factor3)), data = tmpdf) %>% summary() %>% .[[1]] %>% as.data.frame()
    rownames(Prob) <- rownames(Prob) %>% ifelse(grepl(" ", .), gsub(" ","",.), .)
    FinalTab[FinalTab$Factor1 == paste0("P_", Factor1), y] <- ifelse(Prob[Factor1,"Pr(>F)"] < 0.0001, "<0.0001", sprintf("%.4f",Prob[Factor1,"Pr(>F)"]))
    FinalTab[FinalTab$Factor1 == paste0("P_", Factor2), y] <- ifelse(Prob[Factor2,"Pr(>F)"] < 0.0001, "<0.0001", sprintf("%.4f",Prob[Factor2,"Pr(>F)"]))
    FinalTab[FinalTab$Factor1 == paste0("P_", Factor3), y] <- ifelse(Prob[Factor3,"Pr(>F)"] < 0.0001, "<0.0001", sprintf("%.4f",Prob[Factor3,"Pr(>F)"]))
    FinalTab[FinalTab$Factor1 == paste0("P_", Factor1, ":",Factor2), y] <- ifelse(Prob[paste0(Factor1, ":",Factor2),"Pr(>F)"] < 0.0001, "<0.0001", sprintf("%.4f",Prob[paste0(Factor1, ":",Factor2),"Pr(>F)"]))
    FinalTab[FinalTab$Factor1 == paste0("P_", Factor1, ":",Factor3), y] <- ifelse(Prob[paste0(Factor1, ":",Factor3),"Pr(>F)"] < 0.0001, "<0.0001", sprintf("%.4f",Prob[paste0(Factor1, ":",Factor3),"Pr(>F)"]))
    FinalTab[FinalTab$Factor1 == paste0("P_", Factor2, ":",Factor3), y] <- ifelse(Prob[paste0(Factor2, ":",Factor3),"Pr(>F)"] < 0.0001, "<0.0001", sprintf("%.4f",Prob[paste0(Factor2, ":",Factor3),"Pr(>F)"]))
    FinalTab[FinalTab$Factor1 == paste0("P_", Factor1, ":", Factor2, ":", Factor3), y] <- ifelse(Prob[paste0(Factor1, ":", Factor2, ":", Factor3),"Pr(>F)"] < 0.0001, "<0.0001", sprintf("%.4f",Prob[paste0(Factor1, ":", Factor2, ":", Factor3),"Pr(>F)"]))
    # Eta2
    FinalTab[FinalTab$Factor1 == paste0("Eta2_", Factor1), y] <- (Prob[Factor1,"Sum Sq"] / sum(Prob[,"Sum Sq"])) %>% sprintf("%.4f",.)
    FinalTab[FinalTab$Factor1 == paste0("Eta2_", Factor2), y] <- (Prob[Factor2,"Sum Sq"] / sum(Prob[,"Sum Sq"])) %>% sprintf("%.4f",.)
    FinalTab[FinalTab$Factor1 == paste0("Eta2_", Factor3), y] <- (Prob[Factor3,"Sum Sq"] / sum(Prob[,"Sum Sq"])) %>% sprintf("%.4f",.)
    FinalTab[FinalTab$Factor1 == paste0("Eta2_",Factor1, ":",Factor2), y] <- (Prob[paste0(Factor1, ":",Factor2),"Sum Sq"] / sum(Prob[,"Sum Sq"])) %>% sprintf("%.4f",.)
    FinalTab[FinalTab$Factor1 == paste0("Eta2_",Factor1, ":",Factor3), y] <- (Prob[paste0(Factor1, ":",Factor3),"Sum Sq"] / sum(Prob[,"Sum Sq"])) %>% sprintf("%.4f",.)
    FinalTab[FinalTab$Factor1 == paste0("Eta2_",Factor2, ":",Factor3), y] <- (Prob[paste0(Factor2, ":",Factor3),"Sum Sq"] / sum(Prob[,"Sum Sq"])) %>% sprintf("%.4f",.)
    FinalTab[FinalTab$Factor1 == paste0("Eta2_",Factor1, ":", Factor2, ":", Factor3), y] <- (Prob[paste0(Factor1, ":",Factor2, ":", Factor3),"Sum Sq"] / sum(Prob[,"Sum Sq"])) %>% sprintf("%.4f",.)
    FinalTab[FinalTab$Factor1 == paste0("Eta2_Residuals"), y] <- (Prob["Residuals","Sum Sq"] / sum(Prob[,"Sum Sq"])) %>% sprintf("%.4f",.)
    # MSD and mean values
    tmp_model1 <- aov(as.formula(paste0(y,"~", Factor1, "*", Factor2,"*", Factor3)), data = tmpdf)
    tmp_Tukey_result <- HSD.test(tmp_model1,Factor1, group = T)
    FinalTab[FinalTab$Factor1 == paste0("MSD_",Factor1), y] <- sprintf(HSD_fun(q = as.numeric(tmp_Tukey_result$parameters["StudentizedRange"]),
                                                                               s = as.numeric(tmp_Tukey_result$statistics["MSerror"]),
                                                                               n_vector = c(tmp_Tukey_result$means$r),
                                                                               n_group =  tmp_Tukey_result$parameters$ntr), fmt = paste0("%.",DeciNum[n]+1,"f"))
    for (i in unique(tmpdf[,Factor1])){
      FinalTab[FinalTab$Factor1 == i, y] <-
        paste(sprintf(tmp_Tukey_result$groups[i,1], fmt = paste0("%.",DeciNum[n],"f")), tmp_Tukey_result$groups[i,"groups"])
    }

    tmp_model2 <- aov(as.formula(paste0(y,"~", Factor1, "*", Factor2,"*", Factor3)), data = tmpdf)
    tmp_Tukey_result <- HSD.test(tmp_model2,Factor2, group = T)
    FinalTab[FinalTab$Factor1 == paste0("MSD_",Factor2), y] <- sprintf(HSD_fun(q = as.numeric(tmp_Tukey_result$parameters["StudentizedRange"]),
                                                                               s = as.numeric(tmp_Tukey_result$statistics["MSerror"]),
                                                                               n_vector = c(tmp_Tukey_result$means$r),
                                                                               n_group =  tmp_Tukey_result$parameters$ntr), fmt = paste0("%.",DeciNum[n]+1,"f"))
    for (j in unique(tmpdf[,Factor2])){
      FinalTab[FinalTab$Factor1 == j, y] <-
        paste(sprintf(tmp_Tukey_result$groups[j,1], fmt = paste0("%.",DeciNum[n],"f")), tmp_Tukey_result$groups[j,"groups"])
    }

    tmp_model3 <- aov(as.formula(paste0(y,"~", Factor1, "*", Factor2,"*", Factor3)), data = tmpdf)
    tmp_Tukey_result <- HSD.test(tmp_model3,Factor3, group = T)
    FinalTab[FinalTab$Factor1 == paste0("MSD_",Factor3), y] <- sprintf(HSD_fun(q = as.numeric(tmp_Tukey_result$parameters["StudentizedRange"]),
                                                                               s = as.numeric(tmp_Tukey_result$statistics["MSerror"]),
                                                                               n_vector = c(tmp_Tukey_result$means$r),
                                                                               n_group =  tmp_Tukey_result$parameters$ntr), fmt = paste0("%.",DeciNum[n]+1,"f"))
    for (k in unique(tmpdf[,Factor3])){
      FinalTab[FinalTab$Factor1 == k, y] <-
        paste(sprintf(tmp_Tukey_result$groups[k,1], fmt = paste0("%.",DeciNum[n],"f")), tmp_Tukey_result$groups[k,"groups"])
    }
  }
  return(FinalTab)
}
