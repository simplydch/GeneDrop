# Documentation for example data

#' Dataset GeneDropEx_ped
#'
#' This dataframe contains a simulated pedigree for use in gene-dropping
#'
#' This is simulated pedigree data for a population that
#' experiences density dependant selection and has modest migration.
#'
#' It is loosely based on the dynamics seen in the data collected by the
#' Soay Sheep Project.
#'
#' There is a small amount of missing data.
#'
#' It spans 35 years and includes 10000 individuals.  Migrants do not have a known cohort, but
#' the time that they arrived in the population can be found in the 'Arrived' column.
#'
#' Please feel free to play around with this data and let me know how realistic, or otherwise,
#' you think it is.
#'
#' @name GeneDropEx_ped
NULL


#' Dataset GeneDropEx_hap
#'
#' This matrix contains a simulated haplotypes for use in gene-dropping
#'
#' This is simulated haplotype data loosely based on the population history
#' of the soay sheep on St Kilda.
#'
#' It contains 650 haplotypes each with 30000 alleles
#'
#' Please feel free to play around with this data and let me know how realistic, or otherwise,
#' you think it is.
#'
#' @name GeneDropEx_hap
NULL

#' Dataset GeneDropEx_map
#'
#' This dataframe contains a simulated mapping information for use in gene-dropping
#'
#' This is simulated map data that describes the physical and genetic location of the loci
#' in the GeneDropEx_hap file. The distances are loosely based that observed in the
#' St Kilda soay sheep population.
#'
#' Please feel free to play around with this data and let me know how realistic, or otherwise,
#' you think it is.
#'
#' @name GeneDropEx_map
NULL
