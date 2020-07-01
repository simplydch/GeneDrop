### Example code for gene-dropping

## Clear workspace
rm(list = ls())
gc()

library(GeneDrop)

## Inspect the included files

# Pedigree file - A pedigree with 10000 individuals
head(GeneDropEx_ped)
dim(GeneDropEx_ped)

# Founder haplotypes file - 650 founder haplotypes for 30000 loci

dim(GeneDropEx_hap)

GeneDropEx_hap[1:10,1:10]


# Map Files - Map distance information for the loci
head(GeneDropEx_map)



# Extract map distances in cM
map_distances <- GeneDropEx_map[['cMdiff']]

# Find number of loci on each chromosome
chrom_loci_num <- unname(table(GeneDropEx_map[, 'Chr']))

# Extract the founder haplotypes
founder_haplotypes<-GeneDropEx_hap

### Add a unique '3' allele that we can track later

founder_haplotypes[255,1024]<-3

## Run fix_pedigree on the file to ensure it is in the correct order
# and that required columns are present
pedigree <- fix_pedigree(GeneDropEx_ped)

#  Run est_cohort to estimated cohorts for individuals that don't have one
pedigree[, 'Cohort'] <- est_cohort(pedigree)

# Remove individuals that still don't have a cohort
pedigree <- pedigree[!is.na(pedigree[, 'Cohort']), ]

# Create vector containing founder information, here all
# individuals born before 1990 are considered founders
founders <- pedigree[pedigree[, 'Cohort'] < 1990, 'ID']

length(founders)

# Run complete_ped_links to randomly fill in missing pedigree links
ped_completed <- complete_ped_links(pedigree, founders)

# Set a seed so the results are repeatable
set.seed(4562)

# Perform gene-dropping
gene_drop_01 <- genedrop(
  pedigree = ped_completed,
  map_dist = map_distances,
  chr_loci_num = chrom_loci_num,
  found_hap = founder_haplotypes,
  to_raw = TRUE
)


### Tracking alleles back to ancestors

# We can find founders that have the '3' allele at locus 1024 with the following code
# We double the number of founders because there are two haplotypes per individual.
# We use as.raw as the genotype matrix is stored as raw vectors

founder_3 <- unlist(sapply(1:(length(founders)*2), function(x)
{
  if (get_genotype_matrix(gene_drop_01)[[x]][1024] == as.raw(3)){x}}))

# These row references for the genotype matrix can be converted into
# pedigree row references

id_ref = ifelse (founder_3 %% 2 != 0, (founder_3 + 1) / 2, founder_3 / 2)

## And finally we can extract the IDs

founder_id_3 = get_pedigree(gene_drop_01)[id_ref, 'ID']

# These founder individuals have a copy of the '3 allele, every '3' allele should be able to
# be traced back to these individuals

print(founder_id_3)

# Similarly we can find the individual furthest down the pedigree with the '3' allele

rows_allele_3 = tail(unlist(sapply(1:length(get_genotype_matrix(gene_drop_01)), function(x) {
  if (get_genotype_matrix(gene_drop_01)[[x]][1024] == as.raw(3)) {x}})),1)


# And find the individuals that allele belongs to
id_ref = ifelse(rows_allele_3 %% 2 != 0, (rows_allele_3 + 1) / 2, (rows_allele_3 / 2))

id_track = get_pedigree(gene_drop_01)[id_ref, 'ID']

print(id_track)

# We can track the allele back to its original ancestor
# allele_track_back takes an ID and loci number,
# along with the gene-drop output

tr_3 <- allele_track_back(id_track, 1024, gene_drop_01)

# An object is returned which shows tracks the sire and the dam allele
# back to their original ancestor

# If this tracking has worked the allele should have been tracked back
# to one of the founders with the '3' allele

any(c(get_allele_dam(tr_3)[,'ID'],get_allele_sire(tr_3)[,'ID']) %in% founder_id_3)

# And every genotype position should have the same allele for both the dam:

id_ref <- id_ref(gene_drop_01, get_allele_dam(tr_3)[, 'ID'])
row_ref = mapply(function(x, y)
  c(x * 2 - 1, x * 2)[y] , id_ref, as.numeric(get_allele_dam(tr_3)[, 'Hap']))

unlist(lapply(get_genotype_matrix(gene_drop_01)[row_ref], function(x)
  as.numeric(x[1024])))

# And the sire

id_ref <- id_ref(gene_drop_01, get_allele_sire(tr_3)[, 'ID'])

row_ref = mapply(function(x, y)
  c(x * 2 - 1, x * 2)[y] , id_ref , as.numeric(get_allele_sire(tr_3)[, 'Hap']))

unlist(lapply(get_genotype_matrix(gene_drop_01)[row_ref], function(x)
  as.numeric(x[1024])))


# allele_track_back can also take vectors of IDs and loci for the first two arguments
# and returns a list of objects for each combination e.g.

t_multi <- allele_track_back(c('1',id_track), c(1024, 256), gene_drop_01)



### Tracking alleles through descendants
# Use the allele_track_forw to find all the  descendants  from the first founder with a '3'
# each allele at locus 1024


tr_3_for <- allele_track_forw(founder_id_3[1], 1024, gene_drop_01)

# The following output should produce all 3's as output
# It takes the ID info from the descendants with the first allele at the locus 1024
# and finds the relevant row references to extract the alleles of the individuals
# from the gene-drop genotype matrix

id_ref <- id_ref(gene_drop_01, get_allele_dam(tr_3_for)[,'ID'])
row_ref = mapply(function(x, y) c(x * 2 - 1, x * 2)[y],
                 id_ref, as.numeric(get_allele_dam(tr_3_for)[, 'Hap']))
unlist(lapply(get_genotype_matrix(gene_drop_01)[row_ref], function(x)
  as.numeric(x[1024])))

# allele_track_forw can also take vectors of IDs and loci for the first two arguments
# and returns an output list of allele_track_objects which contains all combinations
# e.g.

tr_3_for_2 <- allele_track_forw(c(founder_id_3[1], founder_id_3[1]), c(100, 1024), gene_drop_01)


# ID's and genotypes of interest can be extracted using extract_genotype_mat()

## To investigate the change in the '3' allele over time

extracted_loci <- extract_genotype_mat(gene_drop_01, ids='ALL', loci=c(1024))

extracted_loci_info<-cbind(extracted_loci, Cohort = get_pedigree(gene_drop_01)[,'Cohort'])

allele_count<-rbind(t(table(extracted_loci_info[,'L_1024_01'],extracted_loci_info[,'Cohort'])))
cohort_size <- cbind(table(get_pedigree(gene_drop_01)[,'Cohort']))

matplot(apply(allele_count, 2, function(x) x/cohort_size), type='l',
        ylab = "Allele Frequency", xlab="Years Passed")

# The whole genotype matrix can be written to file although this will take quite a while and
# the output file will be large.

write.table(extract_genotype_mat(gene_drop_out),
            output.file,row.names=FALSE, col.names=FALSE, quote=FALSE)


