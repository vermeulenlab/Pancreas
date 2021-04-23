#adapted from plot_genes_jitter; Monocle package

function (cds_subset, grouping = "State", min_expr = NULL, cell_size = 0.75, 
          nrow = NULL, ncol = 1, panel_order = NULL, color_by = NULL, 
          plot_trend = FALSE, label_by_short_name = TRUE, relative_expr = TRUE) 
{
  if (cds_subset@expressionFamily@vfamily %in% c("negbinomial", 
                                                 "negbinomial.size")) {
    integer_expression <- TRUE
  }
  else {
    integer_expression <- FALSE
    relative_expr <- TRUE
  }
  if (integer_expression) {
    cds_exprs <- exprs(cds_subset)
    if (relative_expr) {
      if (is.null(sizeFactors(cds_subset))) {
        stop("Error: to call this function with relative_expr=TRUE, you must call estimateSizeFactors() first")
      }
      cds_exprs <- Matrix::t(Matrix::t(cds_exprs)/sizeFactors(cds_subset))
    }
  }
  else {
    cds_exprs <- exprs(cds_subset)
  }
  mat <- as.matrix(cds_exprs)
  col_fun = colorRamp2(c(-2, 0, 2), c("#0d75c3", "white", 
                                      "#cc0a15"))
  z_mat <- apply(mat, 1, scale)
  rownames(z_mat) <- colnames(mat)
  dim(z_mat)
  Heatmap(t(z_mat), col = col_fun, column_split = acinar_DE_unsuper$yfp)
}


function (cds_subset, grouping = "State", min_expr = NULL, cell_size = 0.75, 
          nrow = NULL, ncol = 1, panel_order = NULL, color_by = NULL, 
          plot_trend = FALSE, label_by_short_name = TRUE, relative_expr = TRUE) 
{
  if (cds_subset@expressionFamily@vfamily %in% c("negbinomial", 
                                                 "negbinomial.size")) {
    integer_expression <- TRUE
  }
  else {
    integer_expression <- FALSE
    relative_expr <- TRUE
  }
  if (integer_expression) {
    cds_exprs <- exprs(cds_subset)
    if (relative_expr) {
      if (is.null(sizeFactors(cds_subset))) {
        stop("Error: to call this function with relative_expr=TRUE, you must call estimateSizeFactors() first")
      }
      cds_exprs <- Matrix::t(Matrix::t(cds_exprs)/sizeFactors(cds_subset))
    }
  }
  else {
    cds_exprs <- exprs(cds_subset)
  }
  mat <- as.matrix(cds_exprs)
  col_fun = colorRamp2(c(-3, 0, 3), c("#0d75c3", "white", 
                                      "#cc0a15"))
  z_mat <- apply(mat, 1, scale)
  rownames(z_mat) <- colnames(mat)
  dim(z_mat)
  Heatmap(t(z_mat), col = col_fun, column_split = acinar_DE_unsuper$yfp, 
          show_column_names = FALSE)
}
