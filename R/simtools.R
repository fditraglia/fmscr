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

#
# cover.tables <- lapply(names(cover.list), function(x) unlist(lapply(cover.list[[x]], TeXtable)))
# names(cover.tables) <- names(cover.list)
#
# cover.panels <- lapply(names(cover.tables), function(x) paste(cover.tables[[x]], collapse = '\n \n \\vspace{2em} \n \n'))
# names(cover.panels) <- names(cover.tables)
#
# lapply(names(cover.panels), function(x) cat(cover.panels[[x]], file = paste0("./Results/coverage_", x, ".tex")))
#
#
# width.tables <- lapply(names(width.list), function(x) unlist(lapply(width.list[[x]], TeXtable)))
# names(width.tables) <- names(width.list)
#
# width.panels <- lapply(names(width.tables), function(x) paste(width.tables[[x]], collapse = '\n \n \\vspace{2em} \n \n'))
# names(width.panels) <- names(width.tables)
#
# lapply(names(width.panels), function(x) cat(width.panels[[x]], file = paste0("./Results/width_", x, ".tex")))
#
#
# #Clean up
# rm(TeXtable, cover.list, width.list)
# rm(cover.panels, cover.tables)
# rm(width.panels, width.tables)
