#####                                       #####
##### Additional functions required by      #####
##### gene-dropping program                 #####
#####                                       #####


#' Function to find most common value in a vector
#'
#' A function that finds vhich values in a vector which are most common
#'
#' @param values vector. A vector of values
#' @keywords pedigree complete fix
#' @export
#' @examples
#'
find_most_common <- function(values) {
  unique_values <- unique(values)
  value_freqs <- tabulate(match(values, unique_values))
  num_most_common <- max(tabulate(match(values, unique_values)))
  unique_values[which(value_freqs == num_most_common)]
  list(most_common = unique_values[which(value_freqs == num_most_common)], number = num_most_common)
}


### Check that valid loci are provided to function

.check_loci <- function(loci, gene_drop_out) {

  loci <- tryCatch(as.numeric(loci), error=function(e) e, warning =function(w) w)

  if ((is(loci, 'warning')) || length(loci) != 1) {
    stop("Provided loci is not a numeric value",  call. = FALSE)
  }

  loci_num <- length(get_genotype_matrix(gene_drop_out)[[1]])

  if (!(loci > 0 & loci <= loci_num)) {
    stop(paste0("Please Check Loci Number ", loci), call. = FALSE)
  }
  return(loci)}
