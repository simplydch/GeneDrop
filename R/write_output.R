#' @describeIn gene_drop_object-class A method to write the genotype matrix to a file.
#' For a large matrix this will take a while and create a large file.
#' @param filename A string to be used as the output file names.  Default 'GeneDrop_out'
#' @param format_short Logical, Default is FALSE.  When FALSE there is two lines (the two haplotypes) per individual,
#' when TRUE the file instead contains one line per individual with two columns per loci.
#' @export


setGeneric("write_file", function(gene_drop_object, filename= 'GeneDrop_out',
                                  format_short = FALSE) standardGeneric("write_file"))
setMethod(
  "write_file",
  c("gene_drop_object"),
  function(gene_drop_object, filename, format_short) {

  # Check format_short is logical
  if (!is.logical(format_short)) {stop('format_short is not TRUE or FALSE')}

  # Calculate values and number of lines in output file
  out_lines <- length(get_genotype_matrix(gene_drop_object))
  num_ind <- out_lines / 2
  num_loci <- length(get_genotype_matrix(gene_drop_object)[[1]])
  out_lines <- ifelse(format_short == TRUE, out_lines/2, out_lines)

  # Calculate steps for process bar
  split_50 <- floor(out_lines/50)

  # Print process bar
  cat("Writing", num_ind, "individuals over", out_lines ,"lines\n")
  cat('0%|',paste0(rep('-',50)),'|100%','\n', sep='')
  cat('  |')
  gc()

  # Open output file
  output.file <- file(file, "wb")
  # When short format is requested join the two haplotypes into a single line
  if(format_short == TRUE) {
    for (i in seq(1, (out_lines*2) -1, 2)){
      j = i + 1
      gen_sub <-do.call(rbind,slot(gene_drop_01,'genotype_matrix')[c(i,j)])
      # Make sure the connection to the file hasn't been closed
      fileopen <- tryCatch(isOpen(output.file), error = function(e) e,
                                            warning = function(w)w)
      if (is(fileopen, 'error')){output.file <- file(file, "ab")}

      # Write line to file
      write.table(rbind(as.numeric(c(gen_sub))), output.file,
                  append=TRUE,row.names=FALSE, col.names=FALSE, quote=FALSE)

      # Print progress
      if (i %% split_50 == 0){cat('-')}
         }

  } else {
      # Print each line from the genotype matrix
    for (i in 1:out_lines){
        gen_sub <-slot(gene_drop_01,'genotype_matrix')[[i]]
      # Make sure the connection to the file hasn't been closed
      fileopen <- tryCatch(isOpen(output.file), error = function(e) e,
                           warning = function(w)w)
      if (is(fileopen, 'error')){output.file <- file(file, "ab")}
      # Write line to file
      write.table(rbind(as.numeric(gen_sub)), output.file,
                  append=TRUE,row.names=FALSE, col.names=FALSE, quote=FALSE)
      # Print progress
      if (i %% split_50 == 0){cat('-')}
    }
  }
  close(output.file)
  cat('|\n')
  cat('Genotype matrix written to', file, '\n')
  })
