params <- expand.grid(p = 1:6/10, r = 0:5/10, B = 1:3, N = c(50, 100, 500))
results <- data.frame(params, size = rpois(nrow(params), 75))

rtables <- xtabs(size ~ r + p + N + B, results)
#rtables <- xtabs(size ~ r + p + N, results)
dim_vals <- dimnames(rtables)
dim_names <- names(dim_vals)
dim_labels <- lapply(seq_along(dim_vals),
                     function(i) paste(dim_names[i], "=", dim_vals[[i]]))

list_names <- apply(expand.grid(dim_labels[3:length(dim(rtables))]),
                    1, paste, collapse = ', ')
table_list <- plyr::alply(rtables, 3:length(dim(rtables)))
names(table_list) <- list_names

# Then lapply a slightly modified version of TeXtable
# Optionally escape greek letters



TeXtable <- function(Rtable){
  nRow <- nrow(Rtable)
  nCol <- ncol(Rtable)
  body <- rbind(colnames(Rtable), apply(Rtable, c(1,2), format))
  body <- apply(body, c(1,2), function(x) paste0('$', x, '$'))
  rownames(body) <- NULL
  colnames(body) <- NULL
  body <- apply(body, 1, function(x) paste(x, collapse = ' & '))
  header <- body[1]
  body <- body[-1]
  body[nRow %/% 2] <- paste0("$\\gamma\\quad$", body[nRow %/% 2])#rowLegend
  body <- paste(body, collapse = "\\\\ \n")
  colLegend <- paste0("&\\multicolumn{", nCol - 1, "}{c}{$\\rho$}")
  header <- paste("\\hline\\hline\n", colLegend, "\\\\ \n", header, "\\\\ \n \\hline")
  header <- paste0("\\begin{tabular}{r|", paste(rep("r", nCol - 1), collapse = ""), "}\n", header)
  footer <- "\\\\ \n \\hline \n \\end{tabular}"
  return(paste0(header, body, footer))
}


cover.tables <- lapply(names(cover.list), function(x) unlist(lapply(cover.list[[x]], TeXtable)))
names(cover.tables) <- names(cover.list)

cover.panels <- lapply(names(cover.tables), function(x) paste(cover.tables[[x]], collapse = '\n \n \\vspace{2em} \n \n'))
names(cover.panels) <- names(cover.tables)

lapply(names(cover.panels), function(x) cat(cover.panels[[x]], file = paste0("./Results/coverage_", x, ".tex")))


width.tables <- lapply(names(width.list), function(x) unlist(lapply(width.list[[x]], TeXtable)))
names(width.tables) <- names(width.list)

width.panels <- lapply(names(width.tables), function(x) paste(width.tables[[x]], collapse = '\n \n \\vspace{2em} \n \n'))
names(width.panels) <- names(width.tables)

lapply(names(width.panels), function(x) cat(width.panels[[x]], file = paste0("./Results/width_", x, ".tex")))


#Clean up
rm(TeXtable, cover.list, width.list)
rm(cover.panels, cover.tables)
rm(width.panels, width.tables)
