#' Get the indices of samples nearest to Ri
#'
#' @param dist.mat Distance matrix
#' @param Ri Sample Ri index.
#' @param nbd.method neighborhood method \code{"multisurf"} or \code{"surf"} (no k) 
#' or \code{"relieff"} (specify k). Used by nearestNeighbors().
#' @param k Integer for the number of neighbors (\code{"relieff"} method).
#' @param Ri.radius Radius of the neighborhood (other methods).
#'
#' @return Numeric vector of nearest indices.
#' @export
#' 
get_Ri_nearest2 <- function(dist.mat, Ri, nbd.method, k = 0, Ri.radius = NULL) {
  #browser()
  
  df <- dist.mat %>%
    select(!!Ri) %>% 
    tibble::rownames_to_column(var = "rowname")
  
  if (nbd.method == "relieff") {
    df <- df %>%
      slice_min(order_by = .data[[Ri]], n = k + 1) %>%
      filter(.data[[Ri]] > 0) %>%
      arrange(.data[[Ri]]) # sort by increasing distance
  } else {
    df <- df %>%
      filter((.data[[Ri]] < Ri.radius[Ri]) & (.data[[Ri]] > 0)) %>%
      arrange(.data[[Ri]]) # sort by increasing distance
  }
  
  neighbors <- df %>% pull(rowname) %>% as.integer()
  
  return(neighbors)
}
