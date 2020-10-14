#####                                                   #####
##### Functions required for performing gene-dropping    #####
#####                                                   #####

### Set up class to store gene-dropping output
# Holds the genotype matrix, haplotype information list
# and input pedigree

#' An S4 class to store gene-dropping output
#'
#'
#'
#' @slot genotype_matrix A list of vectors each storing a haplotype
#' @slot haplotype_info A list containing haplotype information for tracing alleles
#' @slot pedgree The pedigree that was passed to the gene-dropping function
#' @slot map_info A matrix containing mapping information and recombination frequencies used
#' @param gene_drop_object A gene_drop_object
#' @export
#' @examples
#'
gene_drop_object <- setClass("gene_drop_object",
                             slots = c(
                               genotype_matrix = "list",
                               haplotype_info = "list",
                               pedigree = "matrix",
                               map_info = "matrix",
                               model_info = "list"
                             )
)


# Default display of object info
#' @export

setMethod(
  "show",
  c("gene_drop_object"),
  function(object) {
    cat("A gene_drop_object containing:\n")
    num_loci <- ifelse (length(get_genotype_matrix(object)) > 0, length(get_genotype_matrix(object)[[1]]),0)

    cat("@genotype_matrix - A genotype matrix with", num_loci, "loci\n")
    cat("@haplotype_info - A list containing haplotype information for tracing alleles\n")
    cat("@pedigree - A pedigree containing", nrow(get_pedigree(object)), "individuals\n")
    cat("@map_info - A matrix containing mapping and recombination frequency information\n")
    cat("@model_info - A list containing the variable passed to the function\n")
  })




#' @describeIn gene_drop_object-class A method to access the pedigree of the gene-drop object
#' @export

setGeneric("get_pedigree", function(gene_drop_object) standardGeneric("get_pedigree"))
setMethod(
  "get_pedigree",
  c("gene_drop_object"),
  function(gene_drop_object) {
    return(slot(gene_drop_object, "pedigree"))
  }
)

#' @describeIn gene_drop_object-class A method to access the haplotype source information of the gene-drop object
#' @export

setGeneric("get_haplotype_info", function(gene_drop_object) standardGeneric("get_haplotype_info"))
setMethod(
  "get_haplotype_info",
  c("gene_drop_object"),
  function(gene_drop_object) {
    return(slot(gene_drop_object, "haplotype_info"))
  }
)

#' @describeIn gene_drop_object-class A method to access the mapping information of the gene-drop object
#' @export

setGeneric("get_map_info", function(gene_drop_object) standardGeneric("get_map_info"))
setMethod(
  "get_map_info",
  c("gene_drop_object"),
  function(gene_drop_object) {
    return(slot(gene_drop_object, "map_info"))
  }
)

#' @describeIn gene_drop_object-class A method to access the variables passed to the model
#' @export

setGeneric("get_model_info", function(gene_drop_object) standardGeneric("get_model_info"))
setMethod(
  "get_model_info",
  c("gene_drop_object"),
  function(gene_drop_object) {
    return(slot(gene_drop_object, "model_info"))
  }
)

#' @describeIn gene_drop_object-class A method to access the genotype matrix of the gene-drop object
#' @export

setGeneric("get_genotype_matrix", function(gene_drop_object) standardGeneric("get_genotype_matrix"))
setMethod(
  "get_genotype_matrix",
  c("gene_drop_object"),
  function(gene_drop_object) {
    return(slot(gene_drop_object, "genotype_matrix"))
  }
)


#' @describeIn gene_drop_object-class A method to extract the indivdual information from the pedigree
#' @export

setGeneric("id_info", function(gene_drop_object, ids) standardGeneric("id_info"))
setMethod(
  "id_info",
  c("gene_drop_object"),
  function(gene_drop_object, ids) {
    return(get_pedigree(gene_drop_object)[match(ids, get_pedigree(gene_drop_object)[, "ID"]), ])
  }
)



#' @describeIn gene_drop_object-class A method to get and individuals row reference from the pedigree
#' @export


setGeneric("id_ref", function(gene_drop_object, ids) standardGeneric("id_ref"))
setMethod(
  "id_ref",
  c("gene_drop_object"),
  function(gene_drop_object, ids) {
    return(match(ids, get_pedigree(gene_drop_object)[, "ID"]))
  }
)


#' @describeIn gene_drop_object-class A method to return a genotype matrix
#' @param ids A scalar or vector of IDs to include. Defaults to ALL
#' @param loci A scalar or vector of loci to include. Defaults to ALL
#' @export

setGeneric("extract_genotype_mat", function(gene_drop_object, ids = "ALL", loci = "ALL") standardGeneric("extract_genotype_mat"))
setMethod(
  "extract_genotype_mat",
  c("gene_drop_object"),
  function(gene_drop_object, ids = "ALL", loci = "ALL") {
    if ((length(ids) == 1 && ids == "ALL") & (length(loci) == 1 && loci == "ALL")) {
      if(slot(gene_drop_object, "model_info")['to_raw']==TRUE){
        return(matrix(as.numeric(do.call(rbind,slot(gene_drop_object, "genotype_matrix"))),
                      nrow =length(slot(gene_drop_object, "genotype_matrix"))))}
      else{
        return(matrix((do.call(rbind,slot(gene_drop_object, "genotype_matrix"))),
                      nrow =length(slot(gene_drop_object, "genotype_matrix"))))}

    } else {
      if (length(ids) == 1 && ids == "ALL"){ids = get_pedigree(gene_drop_object)[,'ID']}
      if (length(loci) == 1 && loci == "ALL"){loci = c(1:length(slot(gene_drop_object, "genotype_matrix")[[1]]))}
      else{loci<-sapply(loci,function(x).check_loci(x, gene_drop_object))} # TODO change this to be more efficient

      id_list <- id_ref(gene_drop_object, ids)
      if (any(is.na(id_list))) {
        stop("ID provided is not in predigree")
      }

      refs <- c(id_list * 2 - 1, id_list * 2)
      id_names <- paste0("ID_", ids)
      loci_names <- paste0("L_", rep(loci, each = 2), "_", rep(c("01", "02"), length(loci)))

      mat_out <-matrix(sapply(get_genotype_matrix(gene_drop_object)[refs],'[',loci), byrow = TRUE,
                       nrow=length(refs))
      if(slot(gene_drop_object, "model_info")['to_raw']==TRUE){
        return(matrix((as.numeric(mat_out)), nrow = length(id_list), dimnames = list(id_names, loci_names)))}
      else{
        return(matrix(((mat_out)), nrow = length(id_list), dimnames = list(id_names, loci_names)))
      }
    }
  }
)



#' Function to calculate recombination frequency
#'
#' A function to calculate recombination frequency based on kosambi map distances
#'
#' @param map_dist A kosambi map distance in cM
#' @keywords recombination frequency kosambi
#' @export
#' @examples
#'
calc_recom_freq <- function(map_dist) {
  (1 - exp(-4 * map_dist)) / 2 * (1 + exp(-4 * map_dist))
}


#' Function to condense haplotype source information
#'
#' A function to condense haplotype source information.
#' A string of 1's and 2's determinining which parental haplotype
#' each allele of the desecendants  haplotype came from are condensed.
#' i.e. (111112222 becomes 1(5)2(4) but in matrix format)
#' A list of two matrices is returned, one relating to the sire and one to the dam.
#'
#' @param sire_hap_info A vector of 1's and 2's determining sire haplotype source information
#' @param dam_hap_info A vector of 1's and 2's determining dam haplotype source information
#' @keywords rle condense haplotype information
#' @export
#' @examples
#'
reduce_hap_code <- function(sire_hap_info, dam_hap_info) {
  z <- rle(sire_hap_info)
  h1 <- rbind(z[[1]], z[[2]])
  z2 <- rle(dam_hap_info)
  h2 <- rbind(z2[[1]], z2[[2]])
  list(h1, h2)
}


### Function to perform gene_dropping
# Requires a pedigree (use fixpedigree), map distances, chromosome lengths(number of loci)
# and the founder haplotypes
# Returns a list of 3
# 1) The genedropped haplotypes
# 2) A hap code list, allowing tracking back to ancestors
# 3) The pedigree that was used in the gene dropping
#  - This is the same pedigree that was provided to the function but
#    returning it keeps all required information in one object

#' Function to perform gene-dropping
#'
#' This function carries out gene-dropping down a provided pedigree
#'
#' @param pedigree A dataframe or matrix. A pedigree containing 'ID', 'Sire' and 'Dam' columns
#' @param map_dist A vector of map distances in cM. The vector should be the same length as the total number of loci in
#' chr_loci_num.  Each value should be the difference in cM between the two markers, such that the first value in the vector
#' is the difference between marker 1 and marker 2.  The final value in the vector is therefore redundant but should be included
#' (e.g. as 0). Alternatively can be set to  NULL if recombination frequencies are provided using recom_freq.
#' @param chr_loci_num A vector containing the number of loci on each chromosome
#' @param found_hap A vector of founder haplotypes, if sample_hap=FALSE the pair of haplotypes for each individual should
#' in the same order as in the provided pedigree
#' @param founders_unk A vector of founders with no associated haplotypes, these individuals have their haplotypes generated
#' during gene-dropping based on allele frequencies
#' @param founders_unk_cohorts A scalar which determines the cohorts used to generated haplotypes for individuals in the
#' founders_unk vector.  The default is 1 which uses only haplotypes from the same cohort, 2 would use the
#' current same cohort along with the previous one and so on.
#' @param to_raw logical TRUE or FALSE.  If TRUE the genotypes is stored are raw vectors.  This allows quicker gene-dropping
#' and reduces the size of the output.  It must be FALSE if alleles do not all have numeric values between 0 and 255
#' @param sample_hap logical TRUE or FALSE.  If TRUE founder haplotypes are assigned to founders at random.  If FALSE
#' pairs of haplotypes are assigned to founders in the order they appear in the pedigree.  When FALSE the number of haplotypes in found_hap
#' must be twice the numeber of founders present in the pedigree
#' @param recom_freq 'kosambi' is default, this setting calculates recombination frequencies based using the map distances
#' provided in map_dist.  Alternatively a vector of the recombination frequencies between each marker can be provided.
#' The vector should be the same length as the total number of loci in chr_loci_num  The first value in the vector should be
#' the recombination frequency between markers 1 and 2. The final value in the vector is therefore redundant but should be included
#' (e.g. as 0.  Values between chromosomes are best set to 0.5 although this happens automatically during gene-dropping.
#' @param progress logical, default TRUE.  Determines in progress should be displayed will gene-dropping
#' @keywords gene-dropping recombination segregation allele tracking
#' @export
#' @examples
#'
genedrop <- function(pedigree, map_dist, chr_loci_num, found_hap, founders_unk = FALSE, founders_unk_cohorts = 1,
                     to_raw = TRUE, sample_hap = TRUE, recom_freq = "kosambi", progress=TRUE) {

  # Check founder_sample
  if (is.logical(founders_unk) && founders_unk == TRUE) {
    stop("founders_unk should be FALSE or a vector")
  }

  # Check founder_unk_cohort

  if (!(length(founders_unk_cohorts)==1 && is.numeric(founders_unk_cohorts) && founders_unk_cohorts > 0)){
    stop("founders_unk_cohorts is not a single numeric value greater than 0")
  }


  # Check pedigree
  ped_col_present(pedigree)

  # Make sure ID's aren't present multiple times
  if (any(duplicated(pedigree[,'ID']))) {
    stop('Some individuals appear more than one in pedigree')
  }

  ### Get model info

  model_info <-list()
  model_info['model_arguments'] <-list(sys.call())
  model_info<-c(model_info,as.list(environment())[6:11])

  gene_drop_out <- new("gene_drop_object", pedigree = pedigree,
                       model_info = model_info)

  if (is.logical(to_raw) & to_raw == FALSE) {
    convert_to_raw <- function(x) {x}
    convert_from_raw <- function(x){x}
  } else if (is.logical(to_raw) & to_raw == TRUE) {
    convert_to_raw <- function(x) {as.raw(x)}
    convert_from_raw <- function(x){as.numeric(x)}
  } else {
    stop("to_raw should be TRUE or FALSE")
  }

  ### Find number of loci
  loci_num <- sum(chr_loci_num)

  # Check founder haplotypes

  # Deal with a vector being passed and convert to matrix

  if (!(is.matrix(found_hap) | is.data.frame(found_hap))){
    found_hap <- matrix(found_hap, nrow = length(found_hap))
  }

  if (dim(found_hap)[2] != loci_num) {
    stop(paste0(
      "number of loci in founder haplotypes (", dim(found_hap)[2],
      ") is not equal to total number of loci in chr_loci_num (", loci_num, ")"
    ))
  }

  if (to_raw == TRUE && any(!is.numeric(found_hap)) && any(found_hap < 0 | found_hap > 255)) {
    stop("to_raw is TRUE but found_hap contains values that are not between 0 and 255", call.=FALSE)
  }

  ### Check map_dist

  if (is.null(map_dist) & length(recom_freq) != loci_num) {
    stop(paste0(
      "map_dist is NULL but recom_freq is not a vector the",
      "same length as the genotype (", loci_num, ")"
    ))
  }

  if (!is.null(map_dist) & length(map_dist) != loci_num) {
    stop(paste0("map_dist is not the same length as the genotype (", loci_num, ")"))
  }


  ### Calculate recombination frequencies

  if (length(recom_freq) == 1 && recom_freq == "kosambi") {
    recom_freq_vec <- apply(matrix(map_dist / 100), 1, calc_recom_freq)
  } else if ((is.vector(recom_freq) && length(recom_freq) == loci_num && all(recom_freq <= 1 & recom_freq >= 0))) {
    recom_freq_vec <- matrix(recom_freq)
  } else {
    stop(paste0(
      "recom_freq format is wrong. It should be 'kosambi' or a vector ",
      "of probabilities (<=1 and >=0) the same length as the genotype (", loci_num, ")"
    ))
  }


  ### Get parent row references
  sire_ref <- match(pedigree[, "Sire"], pedigree[, "ID"], nomatch = 0)
  dam_ref <- match(pedigree[, "Dam"], pedigree[, "ID"], nomatch = 0)


  ### Get references for founders and non_founders
  gd_founders <- which(pedigree[, "Sire"] == 0 & pedigree[, "Dam"] == 0)
  gd_nonfounders <- c(1:nrow(pedigree))[!1:nrow(pedigree) %in% gd_founders]

  ### Get references for founders without associated haplotypes

  if (!is.logical(founders_unk)){
    gd_samp_found <- match(founders_unk, pedigree[, "ID"])
    gd_founders <- gd_founders[!gd_founders %in% gd_samp_found]
    gd_nonfounders_all <- sort(c(gd_nonfounders, gd_samp_found))}
  else{
    gd_nonfounders_all <- gd_nonfounders
  }

  ### Randomly assign founder haplotypes to matrix if sample=TRUE, otherwise add genotypes to
  ### correct founders

  if (is.logical(sample_hap) && sample_hap == TRUE) {
    found_samp <- found_hap[sample(nrow(found_hap), length(gd_founders) * 2, replace = TRUE),, drop = FALSE ]
  } else if (is.logical(sample_hap) && sample_hap == FALSE) {
    if (nrow(found_hap) != length(gd_founders) * 2) {
      stop(paste0(
        "Number of haplotypes (", nrow(found_hap), ") is not twice the number of founders (",
        length(gd_founders), ") but sample_hap=FALSE"
      ), call. = FALSE)
    }
    found_samp <- found_hap
  } else {
    stop(paste0("sample_hap argument should be TRUE or FALSE"),call. = FALSE)
  }

  ### Set up matrix for haplotypes
  gd_hap <- rep(list(NA), nrow(pedigree) * 2)

  hap_split <- split(convert_to_raw(found_samp), 1:nrow(found_samp))

  gd_hap[sort(c(gd_founders * 2 - 1, gd_founders * 2))] <- hap_split

  ### Store mapping information

  map_info <- cbind(Chromosome = rep(1:length(chr_loci_num), chr_loci_num),
                    cMDiff = map_dist,
                    recom_freq = recom_freq_vec)

  ### Set Recombination frequencies between chromosomes as 0.5
  recom_freq_vec <- c(0.5, recom_freq_vec[1:length(recom_freq_vec) - 1])

  recom_freq_vec[c(cumsum(chr_loci_num[1:length(chr_loci_num) - 1]) + 1)] <- 0.5

  ### Set up list for haplotype codes
  hap_code_list <- rep(list(list()), length(gd_nonfounders) * 2)

  # Basic Progress Bar
  if (progress == TRUE){
    cat("Gene-dropping to", length(gd_nonfounders), "individuals:\n")
    split_50 <- floor(length(gd_nonfounders)/50)
    cat('0%|',paste0(rep('-',50)),'|100%','\n', sep='')
    cat('  |')}
  ### Establish genotype for each non-founder individual
  for (ind in 1:length(gd_nonfounders_all)) {

    n <- gd_nonfounders_all[ind]

    if (!is.logical(founders_unk) && n %in% gd_samp_found){
      # get cohort
      foc_cohort <- pedigree[n,'Cohort']
      samp_cohorts <- foc_cohort:(foc_cohort + 1 - founders_unk_cohorts)
      foc_cohort_refs <- which(pedigree[,'Cohort'] %in% samp_cohorts)

      # remove sampled founders
      foc_cohort_refs <- foc_cohort_refs[!foc_cohort_refs %in% gd_samp_found]
      foc_cohort_refs <- c(foc_cohort_refs * 2 -1, foc_cohort_refs * 2)
      sam_hap <- apply(do.call(rbind,gd_hap[foc_cohort_refs]),2,function(x) sample(x,2))

      x_sire_hap<-sam_hap[1,]
      x_dam_hap<-sam_hap[2,]

      gd_hap[[n * 2 - 1]] <- convert_to_raw(x_sire_hap)
      gd_hap[[n * 2]] <- convert_to_raw(x_dam_hap)

      hap_code_list[c(c(ind, ind) + c(ind - 1, ind))] <- NA


      if (progress == TRUE && ind %% split_50 == 0){
        cat('-')}

    } else{


      ### get sire and dam info for focal individual
      x_sire <- c(sire_ref[n] * 2 - 1, sire_ref[n] * 2)
      if (any(x_sire<1)){cat('\n')
        stop("Sire missing, individuals are present with a single missing parent",call. = FALSE)
      }
      x_dam <- c(dam_ref[n] * 2 - 1, dam_ref[n] * 2)
      if (any(x_dam<1)){cat('\n')
        stop("Dam missing, individuals are present with a single missing parent",call. = FALSE)
      }

      ### Establish where recombination events will occur
      ### and which haplotype comes from each chromosome
      crossover <- rbinom(loci_num, 1, recom_freq_vec)

      ### Code new haplotype as 1 or 2 depending on which sire haplotype they come from

      hap1_code <- ifelse(cumsum(crossover) %% 2 == 0, 1, 2)

      ### Get haplotype from sire
      x_sire_hap <- ifelse(hap1_code == 1, convert_from_raw(gd_hap[[x_sire[1]]]), convert_from_raw(gd_hap[[x_sire[2]]]))

      ### Repeat for dams
      crossover <- rbinom(loci_num, 1, recom_freq_vec)

      hap2_code <- ifelse(cumsum(crossover) %% 2 == 0, 1, 2)

      x_dam_hap <- ifelse(hap2_code == 1, convert_from_raw(gd_hap[[x_dam[1]]]), convert_from_raw(gd_hap[[x_dam[2]]]))
      ### Write the genotype to the matrix
      gd_hap[[n * 2 - 1]] <- convert_to_raw(x_sire_hap)
      gd_hap[[n * 2]] <- convert_to_raw(x_dam_hap)

      ### Write haplotype codes to list (for tracing)
      hap_code_list[c(c(ind, ind) + c(ind - 1, ind))] <- reduce_hap_code(hap1_code, hap2_code)

      if (progress == TRUE && ind %% split_50 == 0){
        cat('-')}
    }}
  if (progress == TRUE){
    cat('|\n')}
  slot(gene_drop_out, "haplotype_info") <- hap_code_list
  slot(gene_drop_out, "genotype_matrix") <- gd_hap
  slot(gene_drop_out, "map_info") <- map_info

  return(gene_drop_out)
}
