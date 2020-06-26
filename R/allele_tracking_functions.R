#####                                       #####
##### Functions to track haplotypes         #####
##### through descendants and ancestors     #####
#####                                       #####


### set up a class union to allow allele_track_object to store ID
# as either numeric or character

setClassUnion("numORchr", c("numeric", "character"))

# TODO Is this the best way to store these?

#' An S4 class to store to store allele tracking output
#'
#' @slot id A numeric or character value which corresponds to ID in pedigree of gene_drop_object
#' @slot locus A numeric  value which corresponds to locus position in of gene_drop_object
#' @slot allele_dam A matrix to contain allele tracking information for the dam
#' @slot allele_sire A matrix to contain allele tracking information for the sire
#' @export
#' @examples
#'
allele_track_object <- setClass(
  "allele_track_object",
  slots = c(
    id = "numORchr",
    locus = "numeric",
    allele_dam = "matrix",
    allele_sire = "matrix"
  )
)

#' @export
setMethod("show",
          "allele_track_object",
          function(object) {
            cat("A allele_track_object:\n\n")
            cat("ID : ", slot(object, 'id'), " Locus: ", slot(object, 'locus'), "\n")
            cat("\n  Allele Dam: \n")
            print(slot(object, 'allele_dam'))
            cat("\n  Allele Sire: \n")
            print(slot(object, 'allele_sire'))
          })

#' @describeIn allele_track_object-class A method to access the dam allele section of the allele_track_object
#' @export

setGeneric("get_allele_dam", function(allele_track_object) standardGeneric("get_allele_dam"))
setMethod(
  "get_allele_dam",
  c("allele_track_object"),
  function(allele_track_object) {
    return(slot(allele_track_object, "allele_dam"))
  }
)

#' @describeIn allele_track_object-class A method to access the sire allele section of the allele_track_object
#' @export

setGeneric("get_allele_sire", function(allele_track_object) standardGeneric("get_allele_sire"))
setMethod(
  "get_allele_sire",
  c("allele_track_object"),
  function(allele_track_object) {
    return(slot(allele_track_object, "allele_sire"))
  }
)


###
#' Function to track a gene-dropped allele back through ancestors
#'
#'A function that tracks a gene-dropped individuals alleles at a locus back to the ancestor in
#'which they originated.
#'
#' @param id A single character or numeric value identifying an individual in the pedigree
#' @param loci A single numeric value relating to a locus of interest
#' @param gene_drop_out A gene_drop_object from gene-dropping
#' @keywords allele track ancestor haplotype
#' @export
#' @examples
#'

allele_track_back_sing <- function(id, loci, gene_drop_out) {
  # Check that the loci number is present in the output

  loci<-.check_loci(loci, gene_drop_out)

  # Get the pedigree from the gene-drop output
  ped_use <- get_pedigree(gene_drop_out)

  # Get the row number for the individual
  id_ref <- which(ped_use[, "ID"] == id)

  # Check individual is in the pedigree
  if ((length(id_ref) < 1 || is.na(id_ref))) {
    stop(paste0("ID ", id, " Not in Original Pedigree"))
  }

  # Get the relevant references for the genotype matrix
  id_rows <- c(id_ref * 2 - 1, id_ref * 2)

  # Get sire and dam references
  sire_ref <- match(ped_use[, "Sire"], ped_use[, "ID"], nomatch = 0)
  dam_ref <- match(ped_use[, "Dam"], ped_use[, "ID"], nomatch = 0)

  int_par <- c(sire_ref[id_ref], dam_ref[id_ref])

  # Get founder information
  gd_founders <-
    which(ped_use[, "Sire"] == 0 & ped_use[, "Dam"] == 0)

  #  Set up list of matrices to store output
  allele_track <- rep(list(matrix(NA, nrow = 0, ncol = 2)), 2)
  names(allele_track) <- c("Allele_Sire", "Allele_Dam")

  # Look at each allele in turn
  for (hap in 1:2) {
    # Get parent reference
    par_id_ref <- int_par[hap]
    # Get haplotype row reference
    row_sel <- id_rows[hap]
    # Check if there are parents and proceed if true
    parents <- ifelse(par_id_ref == 0, FALSE, TRUE)
    while (parents == TRUE) {
      # Get the correct haplotype code information from gene-drop output
      par_hap_list <-
        get_haplotype_info(gene_drop_out)[(row_sel - (length(gd_founders)) * 2)][[1]]
      # Take the groups from the haplotype code output
      # groups are stretches between crossovers
      groups <- cumsum(par_hap_list[1,])

      # Find the group that the loci is in
      ref <- unname(table(factor(groups < loci, c("TRUE", "FALSE")))["TRUE"] + 1)
      # Get haplotype info for loci in offspring
      par_hap <- par_hap_list[2, ref]

      # Join information together and add it to the list
      info <- c(par_id_ref, par_hap)
      allele_track[[hap]] <- rbind(allele_track[[hap]], info)

      # Get next sire and dam and continue loop if they are present
      sire <- sire_ref[par_id_ref]
      dam <- dam_ref[par_id_ref]
      if (sire == 0 & dam == 0) {
        parents <- FALSE
        break
      }

      # Get haplotype row reference
      row_sel <- c(par_id_ref * 2 - 1, par_id_ref * 2)[par_hap]
      # Get parent refernce
      par_id_ref <- c(sire, dam)[par_hap]
    }
    # Use cohort info as row names
    row_names <- paste0("Cohort_", ped_use[allele_track[[hap]][, 1], c("Cohort")])
    # Extract info from pedigree and add to output
    allele_track[[hap]] <- cbind( rbind(ped_use[allele_track[[hap]][, 1],
                                  c("ID", "Sex")]),
                                  Hap = allele_track[[hap]][, 2])

    if (length(allele_track[[hap]]) > 0) {
      row.names(allele_track[[hap]]) <- row_names
    }
  }
  new(
    "allele_track_object",
    id = id, locus = loci,
    allele_sire = allele_track[[1]],
    allele_dam = allele_track[[2]]
  )
}

###
#' Function to track gene-dropped alleles back through ancestors
#'
#'A function that tracks provided gene-dropped individuals' alleles at requested loci back
#'to the ancestor in which they originated.  A list of allele_track_objects is returned.
#'
#' @param id A vector containing numerical or character values identifying individuals
#' in the pedigree
#' @param loci A vector containing numerical values relating to the loci of interest
#' @param gene_drop_out A gene_drop_object from gene-dropping
#' @keywords allele track ancestor haplotype
#' @export
#' @examples
#'

allele_track_back <- function(id, loci, gene_drop_out) {
  id_loci <- matrix(c(rep(id, each = length(loci)), rep(loci, length(id))), ncol = 2)
  allele_track_all <- lapply(split(id_loci, 1:nrow(id_loci)),
                             function(x) {allele_track_back_sing(x[[1]], x[[2]], gene_drop_out)})
  # Name list based on ID and loci info
  names(allele_track_all) <- paste("ID", id_loci[, 1], "LOCI", id_loci[, 2], sep = "_")
  if (length(allele_track_all) == 1) {
    return(allele_track_all[[1]])
  } else {
    return(allele_track_all)
  }
}


#' Function to find offspring that have a parents allele
#'
#'A function that finds all an individuals' offspring that received a given alleles.
#'Returns a matrix with containing offspring ID, Sex and which haplotype the allele
#'is on.
#'
#' @param id A vector containing a single numerical or character value identifying the individual
#' of interest
#' @param loci A vector containing single numerical value relating to the loci of interest
#' @param hap The haplotype on which the allele of interest is located (Sire : 1, Dam : 2)
#' @param gene_drop_out A gene_drop_object from gene-dropping
#' @keywords allele track offspring haplotype
#' @export
#' @examples


offspring_with_allele <- function(id, loci, hap, gene_drop_out) {

  loci <- .check_loci(loci, gene_drop_out)

  # Get pedigree from gene-drop output
  ped_use <- get_pedigree(gene_drop_out)

  # Check get ID row reference and check it is valid
  id_ref <- which(ped_use[, "ID"] == id)
  if ((length(id_ref) < 1 || is.na(id_ref))) {
    stop(paste0("ID ", id, " Not in Original Pedigree"))
  }

  #  Find which rows contain offspring of individual
  off_ref <- which(ped_use[, "Sire"] %in% id | ped_use[, "Dam"] %in% id)

  # If there are no offspring return nothing
  if (length(unlist(off_ref)) < 1) {
    return()
  }

  # Get offspring ID and count the number of offspring
  offs_id <- ped_use[unlist(off_ref), c("ID")]
  num_off <- length(off_ref)

  # Get list of parent sex, this method is used over
  # e.g. rep(sex, num) as it allows for cases where parents can
  # be either sex

  sex_list <- c(rep(NA, length(off_ref)))
  sex_list[ped_use[which(ped_use[, "ID"] %in% offs_id), "Sire"] == id] <- 1
  sex_list[ped_use[which(ped_use[, "ID"] %in% offs_id), "Dam"] == id] <- 2

  # Get founder information
  gd_founders <- which(ped_use[, "Sire"] == 0 & ped_use[, "Dam"] == 0)

  ### Get row references
  # TODO Try to find a simpler method.
  # Get row references of genotype matrix
  row_sel <- c(off_ref * 2 - 1, off_ref * 2)
  # Split these into a list of each pair belonging to offspring
  row_sel <- split(row_sel, c(seq(length(row_sel) / 2)))
  # Extract the correct reference from each pair based on parents sex
  # i.e. if individual was a sire we are interested in the first haplotype of offspring
  rowref <- mapply(function(x, y) x[y], row_sel, sex_list)

  # Get the correct haplotype code information from gene-drop output
  par_hap_list <- get_haplotype_info(gene_drop_out)[rowref - (length(gd_founders)) * 2]

  # Take the groups from the haplotype code output
  # groups are stretches between crossovers
  groups <- lapply(par_hap_list, function(x) cumsum(x[1, ]))
  # Find the group that the loci is in
  ref <- lapply(groups, function(x) unname(table(factor(x < loci, c("TRUE", "FALSE")))["TRUE"] + 1))

  # Get haplotype info for loci in offspring
  par_hap <- mapply(function(x, y) y[2, x], ref, par_hap_list)

  # Get list of  individuals that have the parent allele we are looking at
  keep_list <- which(par_hap == rep(hap, num_off))

  # Get IDs
  with_allele <- offs_id[keep_list]

  #  If there are no offspring with the allele of interest return nothing
  if (length(with_allele) < 1) {
    return()
  }

  # Get sex of individual
  with_allele_sex <- ped_use[match(c(with_allele), ped_use[, "ID"]), "Sex"]
  return(cbind(ID = with_allele, Sex = with_allele_sex, Hap = sex_list[keep_list]))
}


###
#' Function to track a gene-dropped allele forward through descendants
#'
#'A function that tracks a gene-dropped individuals alleles at a locus forward through all the descendants
#'which received them.
#'
#' @param id A single character or numeric value identifying an individual in the pedigree
#' @param loci A single numeric value relating to a locus of interest
#' @param gene_drop_out A gene_drop_object from gene-dropping
#' @keywords allele track descendant haplotype
#' @export
#' @examples
#'

allele_track_forw_sing <- function(id, loci, gene_drop_out) {
  loci <- .check_loci(loci, gene_drop_out)
  # Set up list to hold values for both alleles
  follow_hap <- rep(list(matrix(NA, nrow = 0, ncol = 3)), 2)
  names(follow_hap) <- c("Allele_Sire", "Allele_Dam")
  # Get the sex of the individual from gene-drop pedigree
  id_sex <- id_info(gene_drop_out, id)['Sex']

  # For each allele at locus, 1 - sire, 2 - dam
  for (hap in 1:2) {
    # Set up initial information
    new_des <- cbind(id, id_sex, hap)
    gen <- 0
    # Loop until break ( there are no more descendants with allele)
    while (TRUE) {
      # Keep track of how many generations of descendants with the allele there are
      gen <- gen + 1
      # Go through each individual in list and apply the offspring_with_allele  function
      new_des <- lapply(
        split(new_des, 1:nrow(new_des)),
        function(x) offspring_with_allele (x[1], loci, x[3], gene_drop_out = gene_drop_out)
      )
      # Join results together neatly
      new_des <- do.call(rbind, new_des)
      # If there are no descendants stop
      if (is.null(new_des)) {
        break
      }
      # Use generation information as row names
      row.names(new_des) <- rep(paste0("Gen_", gen), nrow(new_des))
      follow_hap[[hap]] <- rbind(follow_hap[[hap]], new_des)
    }
  }
  return(new("allele_track_object",
    id = id, locus = loci, allele_sire = follow_hap[[1]], allele_dam = follow_hap[[2]]
  ))
}

###
#' Function to track gene-dropped alleles forward through descendants
#'
#'A function that tracks provided gene-dropped individuals' alleles at requested loci forward
#'through their descendants.  A list of allele_track_objects is returned.
#'
#' @param id A vector containing numerical or character values identifying individuals
#' in the pedigree
#' @param loci A vector containing numerical values relating to the loci of interest
#' @param gene_drop_out A gene_drop_object from gene-dropping
#' @keywords allele track descendant haplotype
#' @export
#' @examples
#'

allele_track_forw <- function(id, loci, gene_drop_out) {
  id_loci <- matrix(c(rep(id, each = length(loci)), rep(loci, length(id))), ncol = 2)
  allele_track_all <- lapply(split(id_loci, 1:nrow(id_loci)), function(x) {
      allele_track_forw_sing(x[[1]], x[[2]], gene_drop_out)
    })
  # Name list based on ID and loci info
  names(allele_track_all) <-
    paste("ID", id_loci[, 1], "LOCI", id_loci[, 2], sep = "_")
  if (length(allele_track_all) == 1) {
    return(allele_track_all[[1]])
  } else {
    return(allele_track_all)
  }
}



