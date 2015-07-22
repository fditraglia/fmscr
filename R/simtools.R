params <- expand.grid(p = 1:6/100, r = 0:5/10, B = 1:3, N = c(10, 20, 30))
results <- data.frame(params, size = rowSums(params))

rtables <- function(table_formula, results_frame){
  # xtabs sums over all unique combinations of factors
  # here there is only one value per combination
  tables <- xtabs(table_formula, results_frame)
  dim_vals <- dimnames(tables)
  dim_names <- names(dim_vals)
  dim_labels <- lapply(seq_along(dim_vals),
                     function(i) paste(dim_names[i], "=", dim_vals[[i]]))
  list_names <- apply(expand.grid(dim_labels[3:length(dim(tables))]),
                    1, paste, collapse = ', ')
  table_list <- plyr::alply(tables, 3:length(dim(tables)))
  names(table_list) <- list_names
  return(table_list)
}

TeXtable <- function(xtab, tab_name = "", row_lab = NULL, col_lab = NULL){
  nRow <- nrow(xtab)
  nCol <- ncol(xtab)
  if(is.null(row_lab)) row_lab <- names(dimnames(xtab))[1]
  if(is.null(col_lab)) col_lab <- names(dimnames(xtab))[2]
  body <- rbind(colnames(xtab), format(xtab))
  body <- cbind(c(tab_name, rownames(xtab)), body)
  body <- apply(body, c(1,2), function(x) paste0('$', x, '$'))
  rownames(body) <- NULL
  colnames(body) <- NULL
  if(tab_name == "") body[1,1] <- ""
  body <- apply(body, 1, function(x) paste(x, collapse = ' & '))
  header <- body[1]
  body <- body[-1]
  row_lab_i <- nRow %/% 2
  body[row_lab_i] <- paste0("$", row_lab, "\\;\\;\\;$ ", body[row_lab_i])
  body <- paste(body, collapse = "\\\\ \n")
  colLegend <- paste0("&\\multicolumn{", nCol, "}{c}{$", col_lab, "$}")
  header <- paste("\\hline\\hline\n", colLegend, "\\\\ \n",
                  header, "\\\\ \n \\hline")
  header <- paste0("\\begin{tabular}{r|",
                   paste(rep("r", nCol), collapse = ""), "}\n", header)
  footer <- "\\\\ \n \\hline \n \\end{tabular}"
  return(paste0(header, body, footer))
}

