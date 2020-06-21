#####                                                   #####  
##### Functions required for performing gene-dropping    #####
#####                                                   #####

### Set up class to store gene-dropping output
# Holds the genotype matrix, haplotype information list
# and input pedigree
gene_drop_object=setClass("gene_drop_object",
                          slots=c(genotype_matrix="list",
                                  haplotype_info ="list",
                                  pedigree="matrix")
)

# Default display of object info

setMethod("show",
          "gene_drop_object",
          function(object) {
            cat("A gene_drop_object containing:\n")
            cat("@genotype_matrix - A genotype matrix with", ncol(object@genotype_matrix), "loci\n")
            cat("@haplotype_info - A list containing haplotype information for tracing alleles\n")
            cat("@pedigree - A pedigree containing", nrow(object@pedigree), "individuals\n")
          }
)

# A method to extract the indivdual information from the pedigree 

setGeneric("ID_info", function(x, ID) standardGeneric("ID_info") )
setMethod("ID_info",
          c("gene_drop_object"),
          function(x,ID) {
            return (x@pedigree[match(ID,x@pedigree[,'ID']),])
          }
)

# A method to get and individuals row reference from the pedigree 

setGeneric("ID_ref", function(x, ID) standardGeneric("ID_ref") )
setMethod("ID_ref",
          c("gene_drop_object"),
          function(x,ID) {
            return (match(ID,x@pedigree[,'ID']))
          }
)


# A method to return the genotype matrix

#TODO Complete this function to completely deal with loci requests
setGeneric("Genotype_mat", function(x, ID='All', loci='ALL') standardGeneric("Genotype_mat") )
setMethod("Genotype_mat",
          c("gene_drop_object"),
          function(x,ID='All', loci='ALL') {
            if (length(ID)==1 && ID == 'All'){
              return(t(sapply(x@genotype_matrix, as.numeric)))
            }else{
            ids <- ID_ref(gene_drop_01,ID)
            if(any(is.na(ids))){stop('ID provided is not in predigree')}
            refs <- c(ids*2 -1, ids*2)
            id_names <- paste0('ID_',ID)
            loci_names <- paste0('L_',rep(loci, each=2),'_', rep(c('01','02'),length(loci)))
            #refs <- c(refs[seq(1,length(refs),2)], refs[seq(2,length(refs),2)])
            mat_out<-(sapply(x@genotype_matrix[c(refs)], function(x)as.numeric(x[loci])))
            return(matrix(t(mat_out), nrow=length(ids), dimnames= list(id_names, loci_names)))
            }}
)


### Function to calculated recombination frequency

calcrecomfreq <- function(M) {
  (1-exp(-4*M))/2*(1+exp(-4*M))
}


###  Function to reduce haploytype code into summarised runs of 1's and 0's
# e.g. 000001111 becomes 0,1:5,4

Reduce_Hap_Code<-function (x,y){
  z<-rle(x)
  h1<-rbind(z[[1]],z[[2]])
  z2<-rle(y)
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



genedrop<-function(peduse,mapdist,chr_loci_num,foundhap, to_raw=TRUE){
  
  gene_drop_out <- new('gene_drop_object', pedigree = peduse)
  
  if (to_raw==FALSE){
    convert_to_raw<-function(x){x}
  }else{
    convert_to_raw<-function(x){as.raw(x)}}
  
  ### Get parent row references
  SireRef<-match(peduse[,'Sire'],peduse[,'ID'],nomatch=0)
  DamRef<-match(peduse[,'Dam'],peduse[,'ID'],nomatch=0)  
  
  
  ### Get references for founders and non_founders
  GD_founders<-which(peduse[,'Sire']==0 & peduse[,'Dam']==0)
  GD_nonfounders<-c(1:nrow(peduse))[!1:nrow(peduse) %in% GD_founders]
  
  ### Find number of loci
  locinum<-sum(chr_loci_num)
  
  ### Set up matrix for haplotypes
  GD_Hap<-rep(list(NA), nrow(peduse)*2)
  
  ### Randomly assign founder haplotypes to matrix
  found_samp <-foundhap[sample(nrow(foundhap),length(GD_founders)*2, replace=TRUE),]
  GD_Hap[sort(c(GD_founders*2-1,GD_founders*2))]<-split(convert_to_raw(found_samp), 1:nrow(found_samp))
  
  ### Code for testing - inserts a unique allele (3) for one individual so it can be traced
  IDref<-which(peduse[,'ID']==8555)
  #IDrows<-c(IDref*2-1,IDref*2)
  GD_Hap[[IDref*2-1]][1024]<-convert_to_raw(3)
  
  
  ### Calculate reconbination frequencies
  recomfreq<-apply(matrix(mapdist/100),1,calcrecomfreq)
  
  ### Set Recombination frequencies between chromosomes as 0.5
  recomfreq[c(1,cumsum(chr_loci_num[1:length(chr_loci_num)-1]))]<-0.5
  
  ### Set up list for haplotype codes
  Hap_Code_List<-rep(list(list()),length(GD_nonfounders)*2)
  
  
  ### Establish genotype for each non-founder individual
  for (ind in 1:length(GD_nonfounders)){
    n<-GD_nonfounders[ind]
    
    ### get sire and dam info for focal individual
    x_sire<-c(SireRef[n]*2-1,SireRef[n]*2)
    x_dam<-c(DamRef[n]*2-1,DamRef[n]*2)
    
    ### get sire and dam genotype
    # x_sire_geno<-GD_Hap[x_sire,]
    # x_dam_geno<-GD_Hap[x_dam,]
    
    ### Establish where recombination events will occur
    ### and which haplotype comes from each chromosome
    crossover<-rbinom(length(mapdist),1,recomfreq)
    
    ### Code new haplotype as 1 or 2 depending on which sire haplotype they come from
    
    hap1_code<-ifelse(cumsum(crossover)%%2==0,1,2)
    
    ### Get haplotype from sire
    x_sire_hap<-ifelse(hap1_code==1,as.numeric(GD_Hap[[x_sire[1]]]),as.numeric(GD_Hap[[x_sire[2]]]))
    
    ### Repeat for dams
    crossover<-rbinom(length(mapdist),1,recomfreq)
    
    hap2_code<-ifelse(cumsum(crossover)%%2==0,1,2)
    
    x_dam_hap<-ifelse(hap2_code==1,as.numeric(GD_Hap[[x_dam[1]]]),as.numeric(GD_Hap[[x_dam[2]]]))
    
    ### Write the genotype to the matrix
    GD_Hap[[n*2-1]]<-convert_to_raw(x_sire_hap)
    GD_Hap[[n*2]]<-convert_to_raw(x_dam_hap)
    
    ### Write haplotype codes to list (for tracing)
    Hap_Code_List[c(n*2-1,n*2)-(length(GD_founders)*2)]<-Reduce_Hap_Code(hap1_code,hap2_code)
  }
  
  gene_drop_out@haplotype_info <- Hap_Code_List
  gene_drop_out@genotype_matrix <- GD_Hap
    
  return (gene_drop_out)
}