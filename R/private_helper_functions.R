# ---------------------
# make helper functions
# ---------------------

#' CV calc
#' @keywords internal
#' @importFrom stats sd
#'
my_CV <- function(x) sd(x, na.rm = TRUE)/mean(x, na.rm = TRUE)*100

#' FC calc row
#' @keywords internal
#'
my_FC_row <- function(arow, Gone, Gtwo) mean(arow[Gone], na.rm = TRUE)/mean(arow[Gtwo], na.rm = TRUE)

#' FC calc apply
#' @keywords internal
#'
my_FC <- function(some_data,Gone,Gtwo) apply(some_data,1,function(x) my_FC_row(x,Gone,Gtwo))

#------- Function added for sparse matrix calculation ---------------

#' slam private function for matrix apply function
#' @keywords internal
#' @import slam
#'
apply_sparseMatrix <- function(X, MARGIN, FUN, ...) {
  #require(slam)
  stopifnot(MARGIN %in% c(1, 2))
  X2 <- slam::as.simple_triplet_matrix(X) # X could be a dense matrix or sparse matrix (dgCMatrix/dgTMatrix).
  if (MARGIN == 1) {
    res <- slam::rowapply_simple_triplet_matrix(X2, FUN, ...)
  } else {
    res <- slam::colapply_simple_triplet_matrix(X2, FUN, ...)
  }
  return(res)
}

#' slam private function for matrix sweep function
#' @keywords internal
#' @import slam
#' @importFrom methods as
#'
sweep_sparseMatrix <- function(x, margin, stats, fun = "*") {
  stopifnot(margin %in% c(1, 2))
  stopifnot(is.vector(stats))
  f <- match.fun(fun)
  x <- as(x, "dgTMatrix")
  if (margin == 1) {
    stopifnot(nrow(x) == length(stats))
    idx <- x@i + 1
  } else {
    stopifnot(ncol(x) == length(stats))
    idx <- x@j + 1
  }
  x@x <- f(x@x, stats[idx])
  return(x)
}
