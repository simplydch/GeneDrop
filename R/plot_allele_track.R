#####                                       #####
##### Functions to plot allele descent      #####
##### down pedigree                         #####
#####                                       #####


# Function plots blank pedigree and returns all links

.plot_blank_ped <- function(gene_drop_object, cohort_labels) {

  # Get pedigree from gene-drop object
  pedigree <- get_pedigree(gene_drop_object)

  # Extract the cohort (for the x-coordinates)
  x <- pedigree[, "Cohort"]

  # Generate evenly spaced y-coordinates
  y <- do.call(c, lapply(tabulate(x), function(x) seq(1, 10, length = x)))

  # Get sire and dam references
  sire_ref <- match(pedigree[, "Sire"], pedigree[, "ID"])
  dam_ref <- match(pedigree[, "Dam"], pedigree[, "ID"])

  # Get IDs and convert them to references
  ids <- pedigree[, "ID"]
  id_ref <- which(pedigree[, "ID"] %in% ids)

  # Separate coordinates into sires and dams
  x_dam <- x[dam_ref[id_ref]]
  y_dam <- y[dam_ref[id_ref]]
  x_sire <- x[sire_ref[id_ref]]
  y_sire <- y[sire_ref[id_ref]]

  # set plot dimensions, add extra space for cohort labels if required
  if (cohort_labels) {
    par(mar = c(0.5, 4, 0.5, 0.5))
  } else {
    par(mar = c(0.5, 0.5, 0.5, 0.5))
  }

  # set-up plot
  plot(x ~ y,
    ylim = c(max(unique(x)), min(unique(x))),
    col = rgb(0, 0, 0, 0.5), pch = "",
    cex = 1, frame.plot = FALSE, axes = F, ylab = "", xlab = ""
  )
  axis(side = 2, labels = cohort_labels, tick = FALSE, las = 2, at = unique(x), cex.axis = 0.6)

  # return coordinates for plotting
  return(list(x, x_dam, x_sire, y, y_dam, y_sire))
}


# Function to plot links from a single sex to an allele

.plot_allele_desc_sing <- function(gene_drop_object, allele = "sire",
                                   plot_coor, tracking,
                                   col_dam_tran, col_dam_solid,
                                   col_sire_tran, col_sire_solid,
                                   point_cex,
                                   id_labels) {
  # Extract plot coordinates

  x <- plot_coor[[1]]
  x_dam <- plot_coor[[2]]
  x_sire <- plot_coor[[3]]
  y <- plot_coor[[4]]
  y_dam <- plot_coor[[5]]
  y_sire <- plot_coor[[6]]


  # Get the pedigree that is to be plotted

  pedigree <- get_pedigree(gene_drop_object)

  # get the focal id from the tracking object
  id <- slot(tracking, "id")

  # get a reference and sex for the focal id
  id_ref(gene_drop_object, id)
  id_sex <- id_info(gene_drop_object, id)["Sex"]

  # get the tracking of either the sire or the dam

  if (allele == "sire") {
    parent_track <- get_allele_sire(tracking)
  } else {
    parent_track <- get_allele_dam(tracking)
  }

  # get index of which individuals have the sire and dam haplotypes
  hap_sire <- match(c(parent_track[, "Hap"]), 1)[1:(length(parent_track[, "ID"]))]
  hap_dam <- match(c(parent_track[, "Hap"]), 2)[1:(length(parent_track[, "ID"]))]

  # convert indexes to IDs
  sire_links <- hap_sire * c(parent_track[, "ID"])
  dam_links <- hap_dam * c(parent_track[, "ID"])

  # get id references for tracked individuals
  id_refs <- c(id_ref(gene_drop_object, c(id, parent_track[, "ID"])))

  # get coordinates for plotting sire links
  x_sire_plot <- x_sire * match(pedigree[, "ID"], sire_links) / match(pedigree[, "ID"], sire_links)
  y_sire_plot <- y_sire * match(pedigree[, "ID"], sire_links) / match(pedigree[, "ID"], sire_links)


  # plot sire links in appropriate colour
  if (allele == "sire") {
    segments(y, x, y_sire_plot, x_sire_plot, col = col_sire_tran, lwd = 0.5, lty = 1)
  } else {
    segments(y, x, y_sire_plot, x_sire_plot, col = col_dam_tran, lwd = 0.5, lty = 1)
  }

  # get coordinates for plotting dam links

  x_dam_plot <- x_dam * match(pedigree[, "ID"], dam_links) / match(pedigree[, "ID"], dam_links)
  y_dam_plot <- y_dam * match(pedigree[, "ID"], dam_links) / match(pedigree[, "ID"], dam_links)

  # plot dam links in appropriate colour

  if (allele == "sire") {
    segments(y, x, y_dam_plot, x_dam_plot, col = col_sire_tran, lwd = 0.5, lty = 1)
  } else {
    segments(y, x, y_dam_plot, x_dam_plot, col = col_dam_tran, lwd = 0.5, lty = 1)
  }


  # get index of which individuals are male and female

  sexes_sire <- match(c(parent_track[, "Sex"]), 1)[1:(length(parent_track[, "ID"]))]
  sexes_dam <- match(c(parent_track[, "Sex"]), 2)[1:(length(parent_track[, "ID"]))]

  # convert indexes to IDs

  sire_points <- sexes_sire * c(parent_track[, "ID"])
  dam_points <- sexes_dam * c(parent_track[, "ID"])

  # get coordinates of sires

  y_points <- c(y[which(pedigree[, "ID"] %in% sire_points)])
  x_points <- c(x[which(pedigree[, "ID"] %in% sire_points)])

  # plot male symbol in appropriate colour

  points(y_points, x_points, col = col_sire_solid, pch = 15, cex = point_cex)
  if (id_labels) {
    text(y_points, x_points + 0.01, labels = pedigree[, "ID"][which(pedigree[, "ID"] %in% sire_points)], cex = 0.5, pos = 3)
  }

  # get coordinates of dams

  y_points <- c(y[which(pedigree[, "ID"] %in% dam_points)])
  x_points <- c(x[which(pedigree[, "ID"] %in% dam_points)])

  # plot female symbol in appropriate colour

  points(y_points, x_points, col = col_dam_solid, pch = 16, cex = point_cex)
  if (id_labels) {
    text(y_points, x_points - 0.01, labels = pedigree[, "ID"][which(pedigree[, "ID"] %in% dam_points)], cex = 0.5, pos = 3)
  }
}


#' Function to plot inheritance of alleles
#'
#' Function to plot inheritance of both alleles present at a locus in an individual from a gene_drop_object
#'
#' @param gene_drop_object A gene_drop_object from gene-dropping
#' @param id A numerical or character scalar of the focal ID to plot.
#' @param loci A numeric scalar of the locus to plot.
#' @param col_sire A character scalar determining the solid colour used for sires.  Should be a colour name, integer, or hexadecimal string of the form "#rrggbb".
#' The default is "#f1a340".
#' @param col_dam A character scalar determining the solid colour used for dams.  Should be a colour name, integer, or hexadecimal string of the form "#rrggbb".
#' The default is "#998ec3".
#' @param line_trans A numeric scalar between 0 and 1 determining the transparency of the plotted lines. The default is 0.6,
#' setting to 0 hides the points.
#' @param point_trans A numeric scalar between 0 and 1 determining the transparency of the plotted points. The default is 1,
#' setting to 0 hides the points.
#' @param focal_cex A numeric scalar determining the size of the point for the focal individual. The default is 1.5.
#' @param point_cex A numeric scalar determining the size of the points used for the descendants. The default is 0.6, setting
#' this value to 0 allows only the point of the focal individual to be plotted.
#' @param background logical TRUE or FALSE. When TRUE the complete pedigree is plotted in light grey.  Default is TRUE but disabling will
#' make plotting much faster.
#' @param cohort_labels logical TRUE or FALSE. When TRUE cohort labels are printed down the left-hand side. Default is TRUE.
#' @param id_labels logical TRUE or FALSE. When TRUE id labels are printed above each point. Default is FASLE.
#' @keywords allele track plot descendant
#' @export
#' @examples
#'


plot_allele_desc <- function(gene_drop_object,
                             id, loci,
                             col_sire = "#f1a340",
                             col_dam = "#998ec3",
                             line_trans = 0.6,
                             point_trans = 1,
                             focal_cex = 1.5,
                             point_cex = 0.6,
                             background = TRUE,
                             cohort_labels = TRUE,
                             id_labels = FALSE) {


  #TODO allow custom cohort labels to be provided

  # Verify that background, cohort_labels and id_labels are logical

  if (!(is.atomic(background) && length(background) == 1L &&
        is.logical(background))){
    stop('background should be TRUE or FALSE')
  }

  if (!(is.atomic(cohort_labels) && length(cohort_labels) == 1L &&
        is.logical(cohort_labels))){
    stop('cohort_labels should be TRUE or FALSE')
  }

  if (!(is.atomic(id_labels) && length(id_labels) == 1L &&
        is.logical(id_labels))){
    stop('id_labels should be TRUE or FALSE')
  }

  # Verify that focal_cex and point_cex are in correct format

  if (!(is.atomic(focal_cex) && length(focal_cex) == 1L &&
        focal_cex >= 0 && !is.logical(focal_cex))){
    stop('focal_cex is not an integer greater than or equal to 0')
  }

  if (!(is.atomic(point_cex) && length(point_cex) == 1L &&
        point_cex >= 0 && !is.logical(point_cex ))){
    stop('point_cex is not an integer greater than or equal to 0')
  }

  # Verify that point_trans and line_trans are in correct format

  if (!(is.atomic(line_trans) && length(line_trans) == 1L &&
         line_trans >= 0 && line_trans <=1 && !is.logical(line_trans ))){
    stop('line_trans is not a value between 0 and 1')
  }

  if (!(is.atomic(point_trans) && length(point_trans) == 1L &&
        point_trans >= 0 && point_trans <=1 && !is.logical(point_trans ))){
    stop('point_trans is not a value between 0 and 1')
  }

  # Verify that colour codes provided can be interpreted as rgb values

  col_sire_rgb <- tryCatch(if (is.atomic(col_sire) && length(col_sire) == 1L){
    col2rgb(col_sire)
  }else{stop()}, error = function(e) e,
      warning = function(w)w)

  if(is(col_sire_rgb , 'warning') | is(col_sire_rgb , 'error')){
    stop('provided value for col_sire is not a valid colour')
  }

  col_dam_rgb <- tryCatch(if (is.atomic(col_dam) && length(col_dam) == 1L){
    col2rgb(col_dam)
  }else{stop()}, error = function(e) e,
  warning = function(w)w)

  if(is(col_sire_rgb , 'warning') | is(col_dam_rgb , 'error')){
    stop('provided value for col_dam is not a valid colour')
  }


  # Set_up a function that creates a hex colour code, with transparency, from rgb values
  create_rgb <- function(col_code) {

    rgb(col_code[1], col_code[2], col_code[3], col_code[4])
  }

  # set-up colours to be used from provided values

  col_sire_tran <- create_rgb(c(col_sire_rgb / 255, line_trans))
  col_sire_solid <- create_rgb(c(col_sire_rgb / 255, point_trans))
  col_dam_tran <- create_rgb(c(col_dam_rgb / 255, line_trans))
  col_dam_solid <- create_rgb(c(col_dam_rgb / 255, point_trans))


  # Track allele down through descendants

  tracking <- allele_track_forw(gene_drop_object, id = id, loci = loci)

  # Get coordinates for all individuals

  plot_coor <- .plot_blank_ped(gene_drop_object, cohort_labels = cohort_labels)

  x <- plot_coor[[1]]
  x_dam <- plot_coor[[2]]
  x_sire <- plot_coor[[3]]
  y <- plot_coor[[4]]
  y_dam <- plot_coor[[5]]
  y_sire <- plot_coor[[6]]

  # Plot whole pedigree if background is TRUE

  if (background == TRUE) {
    segments(y, x, y_dam, x_dam, col = rgb(0.9, 0.9, 0.9, 0.2), lwd = 0.2)
    segments(y, x, y_sire, x_sire, col = rgb(0.9, 0.9, 0.9, 0.2), lwd = 0.2)
  }

  # Get reference and sex of focal individual

  id_sex <- id_info(gene_drop_object, id)["Sex"]
  focal_id_ref <- id_ref(gene_drop_object, id)

  # plot lines for each allele
  if (nrow(get_allele_sire(tracking)) > nrow(get_allele_dam(tracking))) {
    if (nrow(get_allele_sire(tracking)) > 0) {
      .plot_allele_desc_sing(
        allele = "sire",
        plot_coor = plot_coor,
        tracking = tracking,
        gene_drop_object = gene_drop_object,
        col_sire_tran = col_sire_tran,
        col_sire_solid = col_sire_solid,
        col_dam_tran = col_dam_tran,
        col_dam_solid = col_dam_solid,
        point_cex = point_cex,
        id_labels = id_labels
      )
    }
    if (nrow(get_allele_dam(tracking)) > 0) {
      .plot_allele_desc_sing(
        allele = "dam", plot_coor = plot_coor,
        tracking = tracking,
        gene_drop_object = gene_drop_object,
        col_sire_tran = col_sire_tran,
        col_sire_solid = col_sire_solid,
        col_dam_tran = col_dam_tran,
        col_dam_solid = col_dam_solid,
        point_cex = point_cex,
        id_labels = id_labels
      )
    }
  } else {
    if (nrow(get_allele_dam(tracking)) > 0) {
      .plot_allele_desc_sing(
        allele = "dam", plot_coor = plot_coor,
        tracking = tracking,
        gene_drop_object = gene_drop_object,
        col_sire_tran = col_sire_tran,
        col_sire_solid = col_sire_solid,
        col_dam_tran = col_dam_tran,
        col_dam_solid = col_dam_solid,
        point_cex = point_cex,
        id_labels = id_labels
      )
    }
    if (nrow(get_allele_sire(tracking)) > 0) {
      .plot_allele_desc_sing(
        allele = "sire", plot_coor = plot_coor,
        tracking = tracking,
        gene_drop_object = gene_drop_object,
        col_sire_tran = col_sire_tran,
        col_sire_solid = col_sire_solid,
        col_dam_tran = col_dam_tran,
        col_dam_solid = col_dam_solid,
        point_cex = point_cex,
        id_labels = id_labels
      )
    }
  }

  # plot point for focal individual

  if (id_sex == 1) {
    points(y[focal_id_ref], x[focal_id_ref],
      col = col_sire_solid,
      cex = focal_cex, pch = 15
    )
  } else {
    points(y[focal_id_ref], x[focal_id_ref],
      col = col_dam_solid,
      cex = focal_cex, pch = 16
    )
  }

  # plot ID labels if id_labels is TRUE

  if (id_labels) {
    text(y[focal_id_ref], x[focal_id_ref], labels = id, cex = 0.7, pos = 3)
  }
}
