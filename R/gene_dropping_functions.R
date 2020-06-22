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
#' @export
#' @examples


gene_drop_object=setClass("gene_drop_object",
                          slots=c(genotype_matrix="list",
                                  haplotype_info ="list",
                                  pedigree="matrix")
)


# Default display of object info

setMethod("show",
          c("object"),
          function(object) {
            cat("A gene_drop_object containing:\n")
            cat("@genotype_matrix - A genotype matrix with", ncol(object@genotype_matrix), "loci\n")
            cat("@haplotype_info - A list containing haplotype information for tracing alleles\n")
            cat("@pedigree - A pedigree containing", nrow(object@pedigree), "individuals\n")
          }
)


#' @describeIn gene_drop_object-class A method to extract the indivdual information from the pedigree

setGeneric("id_info", function(gene_drop_object, ids) standardGeneric("id_info") )
setMethod("id_info",
          c("gene_drop_object"),
          function(gene_drop_object, ids) {
            return (gene_drop_object@pedigree[match(ids, gene_drop_object@pedigree[,'ID']),])
          }
)



#' @describeIn gene_drop_object-class A method to get and individuals row reference from the pedigree


setGeneric("id_ref", function(gene_drop_object, ids) standardGeneric("id_ref") )
setMethod("id_ref",
          c("gene_drop_object"),
          function(gene_drop_object, ids) {
            return (match(ids, gene_drop_object@pedigree[,'ID']))
          }
)


#' @describeIn gene_drop_object-class A method to return a genotype matrix
#' @param ids A scalar or vector of IDs to include. Defaults to ALL
#' @param loci A scalar or vector of loci to include. Defaults to ALL

#TODO Complete this function to completely deal with loci requests
setGeneric("genotype_mat", function(gene_drop_object, ids='ALL', loci='ALL') standardGeneric("genotype_mat") )
setMethod("genotype_mat",
          c("gene_drop_object"),
          function(gene_drop_object, ids='ALL', loci='ALL') {
            if (length(ids)==1 && ids == 'ALL'){
              return(t(sapply(gene_drop_object@genotype_matrix, as.numeric)))
            }else{
            id_list <- id_ref(gene_drop_01, ids)
            if(any(is.na(id_list))){stop('ID provided is not in predigree')}
            refs <- c(id_list*2 -1, id_list*2)
            id_names <- paste0('ID_',ids)
            loci_names <- paste0('L_',rep(loci, each=2),'_', rep(c('01','02'),length(loci)))
            #refs <- c(refs[seq(1,length(refs),2)], refs[seq(2,length(refs),2)])
            mat_out<-(sapply(gene_drop_object@genotype_matrix[c(refs)], function(x)as.numeric(x[loci])))
            return(matrix(t(mat_out), nrow=length(id_list), dimnames= list(id_names, loci_names)))
            }}
)


#' Function to calculate recombination frequency
#'
#' A function to calculate recombination frequency based on kosambi map distances
#'
#' @param map_dist A kosambi map distance in cM
#' @keywords recombination frequency kosambi
#' @export
#' @examples



calc_recom_freq <- function(map_dist) {
  (1-exp(-4*map_dist))/2*(1+exp(-4*map_dist))
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

reduce_hap_code <- function (sire_hap_info, dam_hap_info){
  z<-rle(sire_hap_info)
  h1<-rbind(z[[1]],z[[2]])
  z2<-rle(dam_hap_info)
  h2<-rbind(z2[[1]],z2[[2]])
  list(h1,h2)
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
#' @param to_raw logical TRUE or FALSE.  If TRUE the genotypes is stored are raw vectors.  This allows quicker gene-dropping
#' and reduces the size of the output.  It must be FALSE if alleles do not all have numeric values between 0 and 255
#' @param sample_hap logical TRUE or FALSE.  If TRUE founder haplotypes are assigned to founders at random.  If FALSE
#' pairs of haplotypes are assigned to founders in the order they appear in the pedigree.  When FALSE the number of haplotypes in found_hap
#' must be twice the numeber of founders present in the pedigree
#' @param recom_freq 'kosambi' is default, this setting calculates recombination frequencies based using the map distances
#' provided in map_dist.  Alternatively a vector of the recombination frequencies between each marker can be provided.
#' The vector should be the same length as the total number of loci in chr_loci_num  The first value in the vector should be
#' the recombination frequency between markers 1 and 2. The final value in the vector is therefore redundant but should be included
#' (e.g. as 0.  Values between chromosomes are best set to 0.5 although this happens automaticaly during gene-dropping.
#' @keywords gene-dropping recombination segregation allele tracking
#' @export
#' @examples


genedrop<-function(pedigree, map_dist, chr_loci_num, found_hap,
                   to_raw=TRUE, sample_hap=TRUE, recom_freq='kosambi'){


  # Check pedigree
    ped_col_present(pedigree)


    gene_drop_out <- new('gene_drop_object', pedigree = pedigree)

    if (is.logical(to_raw) & to_raw==FALSE){
        convert_to_raw <- function(x){x}
    }else if ( is.logical(to_raw) & to_raw==TRUE){
        convert_to_raw <- function(x){as.raw(x)}
    }else{
        stop('to_raw should be TRUE or FALSE')
    }

  ### Find number of loci
  loci_num <- sum(chr_loci_num)

  # Check founder haplotypes
  #TODO It would be good to find a way to make this check faster, it is very slow!

    if (dim(found_hap)[2] != loci_num){
      stop(paste0("number of loci in founder haplotypes (", dim(found_hap)[2],
                  ") is not equal to total number of loci in chr_loci_num (", loci_num, ")"))
    }

    if (to_raw==TRUE && any(!is.numeric(found_hap)) && any(found_hap < 0 | found_hap > 255)){
       stop('to_raw is TRUE but found_hap contains values that are not between 0 and 255')
      }








  ### Check map_dist

  if (is.null(map_dist) & length(recom_freq) != loci_num){

    stop(paste0("map_dist is NULL but recom_freq is not a vector the",
                "same length as the genotype (",loci_num, ")"))
  }

  if (!is.null(map_dist) & length(map_dist) != loci_num){
    stop(paste0("map_dist is not the same length as the genotype (",loci_num, ")"))
  }


  ### Calculate reconbination frequencies

  if (length(recom_freq) == 1 && recom_freq == 'kosambi'){
    recom_freq_vec <- apply(matrix(map_dist/100), 1, calc_recom_freq)
  }else if ((is.vector(recom_freq) && length(recom_freq) == loci_num && all(recom_freq <= 1 & recom_freq >= 0))){
    recom_freq_vec <- matrix(recom_freq)
  }else{
    stop(paste0("recom_freq format is wrong. It should be 'kosambi' or a vector ",
                "of probabilities (<=1 and >=0) the same length as the genotype (",loci_num, ")"))
  }


  ### Get parent row references
  sire_ref <- match(pedigree[,'Sire'], pedigree[,'ID'], nomatch=0)
  dam_ref <- match(pedigree[,'Dam'], pedigree[,'ID'], nomatch=0)


  ### Get references for founders and non_founders
  gd_founders <- which(pedigree[,'Sire']==0 & pedigree[,'Dam']==0)
  gd_nonfounders <- c(1:nrow(pedigree))[!1:nrow(pedigree) %in% gd_founders]


  ### Randomly assign founder haplotypes to matrix if sample=TRUE, otherwise add genotypes to
  ### correct founders

  if (is.logical(sample_hap) && sample_hap==TRUE){
    found_samp <- found_hap[sample(nrow(found_hap), length(gd_founders)*2, replace=TRUE),]

  }else if(is.logical(sample_hap) && sample_hap==FALSE){
      if (nrow(found_hap) != length(gd_founders)*2){
        stop(paste0('Number of haplotypes (', nrow(found_hap), ') is not twice the number of founders (',
                    length(gd_founders), ') but sample_hap=FALSE'), call. = FALSE)
      }
    found_samp <- found_hap

  }else{
    stop(paste0('sample_hap argument should be TRUE or FALSE'), call. = FALSE)
  }

  ### Set up matrix for haplotypes
  gd_hap <- rep(list(NA), nrow(pedigree)*2)

  hap_split <- split(convert_to_raw(found_samp), 1:nrow(found_samp))

  gd_hap[sort(c(gd_founders*2-1,gd_founders*2))] <- hap_split

  ### Code for testing - inserts a unique allele (3) for one individual so it can be traced
  id_ref <- which(pedigree[,'ID']==8555)

  gd_hap[[id_ref*2-1]][1024] <- convert_to_raw(3)


  ### Set Recombination frequencies between chromosomes as 0.5
  recom_freq_vec <- c(0.5, recom_freq_vec[1:length(recom_freq_vec)-1])

  recom_freq_vec[c(cumsum(chr_loci_num[1:length(chr_loci_num)-1])+1)] <- 0.5

  ### Set up list for haplotype codes
  hap_code_list <- rep(list(list()), length(gd_nonfounders)*2)


  ### Establish genotype for each non-founder individual
  for (ind in 1:length(gd_nonfounders)){
    n <- gd_nonfounders[ind]

  ### get sire and dam info for focal individual
  x_sire <- c(sire_ref[n]*2-1, sire_ref[n]*2)
  x_dam <- c(dam_ref[n]*2-1, dam_ref[n]*2)


  ### Establish where recombination events will occur
  ### and which haplotype comes from each chromosome
  crossover <- rbinom(loci_num, 1, recom_freq_vec)

  ### Code new haplotype as 1 or 2 depending on which sire haplotype they come from

  hap1_code <- ifelse(cumsum(crossover)%%2 == 0, 1, 2)

  ### Get haplotype from sire
  x_sire_hap <- ifelse(hap1_code==1, as.numeric(gd_hap[[x_sire[1]]]), as.numeric(gd_hap[[x_sire[2]]]))

  ### Repeat for dams
  crossover <- rbinom(loci_num, 1, recom_freq_vec)

  hap2_code <- ifelse(cumsum(crossover)%%2 == 0, 1, 2)

  x_dam_hap <- ifelse(hap2_code==1, as.numeric(gd_hap[[x_dam[1]]]), as.numeric(gd_hap[[x_dam[2]]]))

  ### Write the genotype to the matrix
  gd_hap[[n*2-1]] <- convert_to_raw(x_sire_hap)
  gd_hap[[n*2]] <- convert_to_raw(x_dam_hap)

  ### Write haplotype codes to list (for tracing)
  hap_code_list[c(c(ind, ind) + c(ind-1, ind))]<-reduce_hap_code(hap1_code, hap2_code)
}

  gene_drop_out@haplotype_info <- hap_code_list
  gene_drop_out@genotype_matrix <- gd_hap

  return (gene_drop_out)
}