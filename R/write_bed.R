#' @describeIn gene_drop_object-class A method to write gene-drop objects to plink files.
#'It might take a while and hasn't been fully tested, so use with caution.
#' @param filename A string to be used as the output file names.  Default 'GeneDrop_out'
#' @param cM_positions A vector or NULL.  Provides map distance information to be included in the .bim file
#' Either NULL in which case all values will be 0, or a vector the same length as the number of loci with mapping distances in cM.
#' @param gen_positions A vector or NULL.  Provides physical distance information to be included in the .bim file
#' Either NULL in which case all values will be 0, or a vector the same length as the number of loci with distances in base pairs.
#' @param fast Logical. Default is TRUE.  Determines if the whole genotype matrix should be extracted prior to
#' writing the files.  This is faster but requires more memory.  The function automatically switches the value
#' to FALSE if the genotype matrix can't be extracted whole.
#' @export


setGeneric("write_bed", function(gene_drop_object, filename= 'GeneDrop_out',
                                 cM_positions=NULL,
                                 gen_positions=NULL, fast=TRUE) standardGeneric("write_bed"))
setMethod(
  "write_bed",
  c("gene_drop_object"),
  function(gene_drop_object, filename, cM_positions, gen_positions, fast) {

    # full file name
    fileoutbed = paste0(filename,'.bed')
    fileoutbim = paste0(filename,'.bim')
    fileoutfam = paste0(filename,'.fam')

    #

    ### extract information and calculate dimension values
    num_ind = nrow(slot(gene_drop_object, "pedigree"))
    num_loci = length(slot(gene_drop_object, "genotype_matrix")[[1]])
    ex <- num_ind %% 4
    byte_num <- ceiling(num_ind /4)
    bit_num <- ceiling(num_ind /4) * 4
    ref_1 <- seq(1,(num_ind*2),2)
    ref_2 <- seq(2,(num_ind*2),2)

    if (is.null(cM_positions)){
      cM_pos <- 0
    } else if (is.vector(cM_positions) && length(cM_positions)==num_loci){
      cM_pos <- cM_positions
    } else {
      stop("cM_positions are not NULL or a vector the same length as the number of loci")
    }

    if (is.null(gen_positions)){
      g_pos <- 0
    } else if (is.vector(gen_positions) && length(gen_positions)==num_loci){
      g_pos <- gen_positions
    } else {
      stop("gen_positions are not NULL or vector the same length as the number of loci")
    }

    ## Set up matrix to hold allele information for bim file
    bim<-matrix(0, nrow=num_loci, ncol=2)

    ## Set up list to hold warnings for printing at the end
    warninglist<-vector(mode = "list", length = 1)
    warninglist[[1]] <- matrix('\n')
    ## Set up value to store if missing data is found
    missing_value <- FALSE

    ## Ensure there is enough memory to read whole matrix at once
    if (fast == TRUE) {
      geno <- tryCatch(
        do.call(rbind, slot(gene_drop_object, "genotype_matrix")),
        error = function(e) e,
        warning = function(w)w)
      if (is(geno, 'error')) {
        cat(
          "There isn't enough memory to extract the full genotype matrix.",
          "Attempting to write to the file directy. This may take a while!"
        )
        fast <- FALSE
      }
    }

    ## Set up list to hold bytes
    bed<-vector(mode = "list", length = num_loci+1)
    ## Add magic number
    bed[[1]]<-as.raw(c("0x6c","0x1b","0x01"))

    ## Set up progress bar
    cat("Writing", num_ind, "individuals and", num_loci ," loci to", fileoutbed,"\n")
    split_50 <- floor(num_loci/50)
    cat('0%|',paste0(rep(' ',50)),'|100%','\n', sep='')
    cat('  |')

    for (i in 1:num_loci){

      # If geno isn't available get values directly from gene drop object
      if (fast == TRUE){
        snp<-geno[,i]
      } else {
        snp<-sapply(slot(gene_drop_object,"genotype_matrix"),'[',i)
      }

      # Find values in snp and remove NA/0
      values <- unique(snp)
      if(missing_value == FALSE && any(values ==0)){
        miss_warn<-simpleWarning(paste0("\nWarning: alleles coded as 0 at locus ", i,
                                        ", these will be treated as missing in output files.",
                                        "  This may also be true at subsequent loci but will not be checked.\n"))
        missing_value = TRUE}


      values <- values[!values %in% c(0,NA)]

      if (length(values) < 2){
        NMA <- as.raw(99)
        MA<- values[1]
        warninglist[[1]][length(warninglist[[1]])+1] <- paste0("Note: Only one allele at locus ",
                                                               i, ",second allele coded as 99 in .bim file!\n")
      } else {
        if (tabulate(match(snp,values[[1]],0)) >= 0.5 * (num_ind *2)){
          MA <- values[1]
        }else{
          MA <- values[2]
        }
        NMA<-values[!MA==values]}

      if (length(NMA) > 2){
        cat('\n')
        stop(paste0("Loci ", i, " has more than two different alleles, can't write bed file"))
      }

      # Add alleles to bim file
      bim[i, 1] <- ifelse(is.raw(MA), as.numeric(MA), MA)
      bim[i, 2] <- ifelse(is.raw(NMA), as.numeric(NMA), NMA)

      #recode alleles to 1 and 2
      allele_recode_1 <- match(snp[ref_1], c(MA,NMA))
      allele_recode_2 <- match(snp[ref_2], c(MA,NMA))

      #recode genotype to 1, 2, 3 and 4
      geno_recode <- allele_recode_1 + allele_recode_2

      #recode into homo and hetero
      geno_recode[match(geno_recode, 3, 0)> 0] <- '10'
      geno_recode[match(geno_recode, 2, 0)>0] <- '00'
      geno_recode[match(geno_recode, 4, 0)>0] <- '11'
      geno_recode[match(geno_recode, c('10','00','11'),0)<1] <-'01'

      # Add empty bits
      geno_recode_02<-matrix(ncol=bit_num,'00')
      geno_recode_02[1:length(geno_recode)] <- geno_recode

      # split into bytes
      geno_recode_03<-matrix(geno_recode_02, ncol=4, byrow = TRUE)
      # join each byte into single string
      geno_recode_04 <- paste0(geno_recode_03[,4],geno_recode_03[,3],geno_recode_03[,2],geno_recode_03[,1])

      # recode bytes into hex values
      bed[[i+1]] <- as.raw(strtoi(geno_recode_04, base=2))

      # print progress
      if (i %% split_50 == 0){cat('-')}
    }

    ## Open output file for bed, combine hex values and write binary file
    output.file <- file(fileoutbed, "wb")
    writeBin(do.call(c,bed), output.file)
    close(output.file)

    # set up bin
    bim_int <- matrix(NA, ncol = 6, nrow = num_loci)
    bim_int[,1] <-slot(gene_drop_object,'map_info')[,'Chromosome']
    bim_int[,2]<-1:num_loci
    bim_int[,3]<-cM_pos
    bim_int[,4]<-g_pos
    bim_int[,c(5,6)]<-bim

    fam_int <- slot(gene_drop_object,'pedigree')[,c(1,1,2,3,4,5)]
    fam_int[,1]<-1
    fam_int[,6]<--9
    write.table(bim_int, fileoutbim, row.names = FALSE, col.names = FALSE, quote=FALSE)
    write.table(fam_int, fileoutfam, row.names = FALSE, col.names = FALSE, quote=FALSE)
    cat('|\n\n')
    cat("Finished writing plink files :", fileoutbed, fileoutfam, fileoutbim,'\n')
    if(missing_value){message(miss_warn)}
    cat(do.call(c,warninglist), sep='')}
)
