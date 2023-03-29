mcc <- function(TP, FP, TN, FN, convert_negative_to_zero = FALSE) {
  # Convert types to double for better precision
  TP <- as.double(TP)
  FP <- as.double(FP)
  TN <- as.double(TN)
  FN <- as.double(FN)
  # Calculate MCC
  numerator <- (TP * TN - FP * FN)
  denominator <- sqrt((TP + FP)*(TP + FN)*(TN + FP)*(TN + FN))
  if (!convert_negative_to_zero) {
    return(numerator/denominator)
  } else {
    res <- numerator/denominator
    res[res < 0] <- 0
    return(res)
  }
}

meltTruthTableGroup <- function(dt, col = "codons", grp_cols = c("simulation", "classifier"), cutoff = Inf) {
  dt_temp <- dt[, .(TP = sum(TP), FP = sum(FP), FN = sum(FN), TN = sum(TN), size = .N),
                by = .(simulation, classifier, get(col))]; setnames(dt_temp, "get", col)

  setorderv(dt_temp, c(grp_cols, col), order = c(rep(1, length(grp_cols)), 1))
  dt_temp <- cbind(dt_temp[,c(grp_cols, col), with = FALSE],
                   round(dt_temp[, .(precision = TP / (TP + FP), recall = TP / (TP + FN),
                                     accuracy = (TP + TN) / (TP + FP + FN + TN),
                                     mcc = mcc(TP, FP, TN, FN))], 3)) #size = size
  dt_temp[cbind(FALSE, is.nan(as.matrix(dt_temp[,-1])))] <- 0
  dt_temp <- dt_temp[get(col) < cutoff,]
  dt_melt <- data.table::melt.data.table(dt_temp, id.vars = c(grp_cols, col), variable.name = "metric")
  if ("classifier" %in% grp_cols) dt_melt[, classifier := as.factor(classifier)]
  dt_melt[]
  return(dt_melt)
}
