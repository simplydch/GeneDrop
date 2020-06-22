#####                                       #####  
##### Functions to track haplotypes         #####
##### through decsendants and ancestors     #####
#####                                       #####

# Class to store allele tracking output
# I'm not sure that this is the best way to present these
# but it is a bit tidier than a list of lists

hap_track_object=setClass("hap_track_object",
                          slots=c(ID="numeric",
                                  Locus ="numeric",
                                  Allele_Dam = "matrix",
                                  Allele_Sire = "matrix")
)

setMethod("show",
          "hap_track_object",
          function(object) {
            cat("A hap_anc_object:\n\n")
            cat("ID : ", object@ID, " Locus: ", object@Locus,"\n")
            cat("\n  Allele Dam: \n")
            print(object@Allele_Dam)
            cat("\n  Allele Sire: \n")
            print(object@Allele_Sire)
          })



### Function to track a single ID and Loci
#


trackhapsing<-function(ID,loci,genedropout){
  
  
  # Check that the loci number is present in tht output
  locinum<-length(genedropout@genotype_matrix[[1]])
  if(!(loci>0 & loci <= locinum)){stop(paste0("Please Check Loci Number ",loci))}
  
  # Get the pedigree from the gene-drop output
  peduse<-genedropout@pedigree
  
  # Get the row number for the individual
  IDref<-which(peduse[,'ID']==ID)
  # Check individual is in the pedigree
  if((length(IDref)<1 ||  is.na(IDref))){stop(paste0("ID ",ID," Not in Original Pedigree"))}
  
  # Get the relevant references for the genotype matrix
  IDrows<-c(IDref*2-1,IDref*2)
  
  # Get sire and dam references
  SireRef<-match(peduse[,'Sire'],peduse[,'ID'],nomatch=0)
  DamRef<-match(peduse[,'Dam'],peduse[,'ID'],nomatch=0) 
  
  intpar<-c(SireRef[IDref],DamRef[IDref])
  
  # Get founder information
  GD_founders<-which(peduse[,'Sire']==0 & peduse[,'Dam']==0)
  
  #  Set up list of matrices to store output
  happed<-rep(list(matrix(NA,nrow=0,ncol=2)),2)
  names(happed)<-c("Allele_Sire","Allele_Dam")
  
  # Look at each allele in turn
  for(s in 1:2){
    # Get parent reference
    parID_ref<-intpar[s]
    # Get haplotype row reference
    rowsel<-IDrows[s]
    # Check if there are parents and proceed if true
    ifelse(parID_ref==0,parents<-FALSE, parents<-TRUE)
    while(parents==TRUE){
      # Get the correct haplotype code information from gene-drop output
      par_hap_list<-genedropout@haplotype_info[(rowsel-(length(GD_founders))*2)][[1]]
      # Take the groups from the haplotype code output
      # groups are stretches between crossovers 
      groups<-cumsum(par_hap_list[1,])
      
      # Find the group that the loci is in
      ref<-unname(table(factor(groups<loci,c("TRUE","FALSE")))['TRUE']+1)
      # Get haplotype info for loci in offspring
      par_hap<-par_hap_list[2,ref]
      
      # Join information together and add it to the list
      info<-c(parID_ref,par_hap)
      happed[[s]]<-rbind(happed[[s]],info)
      
      # Get next sire and dam and continue loop if they are present
      sire<-SireRef[parID_ref]
      dam<-DamRef[parID_ref]
      if(sire==0 & dam==0){parents<-FALSE;break}
      
      # Get haplotype row reference 
      rowsel<-c(parID_ref*2-1,parID_ref*2)[par_hap]
      # Get parent refernce
      parID_ref<-c(sire,dam)[par_hap]
    }
    # Use cohort info as row names
    row_names = paste0('Cohort_', peduse[happed[[s]][,1],c('Cohort')])
    # Extract info from pedigree and add to output
    happed[[s]]<-cbind(rbind(peduse[happed[[s]][,1],c('ID','Sex')]),Hap=happed[[s]][,2])
    if (length(happed[[s]])> 0){row.names(happed[[s]]) <- row_names}
        }
  new('hap_track_object', ID=ID, Locus=loci, Allele_Sire=happed[[1]], Allele_Dam=happed[[2]])
  #happed
}

### Extend trackhapsing function to except multiple input

trackhap<-function(ID,loci,genedropout){
  IDLoci<-matrix(c(rep(ID,each=length(loci)),rep(loci,length(ID))), ncol=2)
  happed_all<-lapply(split(IDLoci,1:nrow(IDLoci)),function(x) trackhapsing(x[[1]],x[[2]],genedropout))
  # Name list based on ID and loci info
  names(happed_all)<-paste('ID',IDLoci[,1],'LOCI',IDLoci[,2],sep='_')
  if (length(happed_all) == 1){
    return(happed_all[[1]])
  }else{
  return(happed_all)
  }}


### Function that finds all offspring of an indiviudal (ID) that have a 
# given parental haplotype(1-Sire, 2-Dam) at a loci

findoffhap<-function(ID,loci,hap,genedropout){
  
  # Check that the loci number is present in the output
  locinum<-length(genedropout@genotype_matrix[[1]])
  if(!(loci>0 & loci <= locinum)){stop(paste0("Please Check Loci Number ",loci))}
  
  # Get pedigree from gene-drop output
  peduse<-genedropout@pedigree
  
  # Check get ID row reference and check it is valid
  IDref<-which(peduse[,'ID']==ID)
  if((length(IDref)<1 ||  is.na(IDref))){stop(paste0("ID ",ID," Not in Original Pedigree"))}
  
  #  Find which rows contain offspring of individual
  Off_ref<-which(peduse[,'Sire'] %in% ID | peduse[,'Dam'] %in% ID)

  # If there are no offspring return nothing
  if(length(unlist(Off_ref))<1){return()}
  
  # Get offspring ID and cound the number of offspring
  OffS_ID<-peduse[unlist(Off_ref),c('ID')]
  Num_Off<-length(Off_ref)
  
  # Get list of parent sex, this method is used over
  # e.g. rep(sex, num) as it allows for cases where parents can
  # be either sex
  
  sex_list = c(rep(NA, length(Off_ref)))
  sex_list[peduse[which(peduse[,'ID'] %in% OffS_ID),'Sire']==ID]<-1
  sex_list[peduse[which(peduse[,'ID'] %in% OffS_ID),'Dam']==ID]<-2
  
  # Get founder information
  GD_founders<-which(peduse[,'Sire']==0 & peduse[,'Dam']==0)
  
  ### Get row references
  # TODO Try to find a simpler method.
  # Get row references of genotype matrix
  rowsel<-c(Off_ref*2-1, Off_ref*2)
  # Split these into a list of each pair belonging to offspring
  rowsel<-split(rowsel, c(seq(length(rowsel)/2)))
  # Extract the correct reference from each pair based on parents sex
  # i.e. if individual was a sire we are interested in the first haplotype of offspring
  rowref <- mapply(function(x,y) x[y], rowsel, sex_list)
  
  # Get the correct haplotype code information from gene-drop output
  par_hap_list<-genedropout@haplotype_info[rowref-(length(GD_founders))*2]
  
  # Take the groups from the haplotype code output
  # groups are stretches between crossovers 
  groups<-lapply(par_hap_list,function(x)cumsum(x[1,]))
  # Find the group that the loci is in
  ref<-lapply(groups,function(x)unname(table(factor(x<loci,c("TRUE","FALSE")))['TRUE']+1))
  
  # Get haplotype info for loci in offspring
  par_hap<-mapply(function(x,y) y[2,x],ref,par_hap_list)
 
  # Get list of  indiviudals that have the parent allele we are looking at
  keeplist<-which(par_hap==rep(hap,Num_Off))
  
  # Get IDs
  WithAllele<-OffS_ID[keeplist]
  
  #  If there are no offspring with the allele of interest return nothing
  if(length(WithAllele)<1){return()}
  
  # Get sex of individual
  WithAllele_Sex<-peduse[match(c(WithAllele), peduse[,'ID']),'Sex']
  return(cbind(ID=WithAllele,Sex=WithAllele_Sex,Hap=sex_list[keeplist]))
}


### Function that finds all descendants of an indiviudal (ID) with each
# allele at a locus

pedfollowhapsing<-function(ID,loci,genedropout){
  # Set up list to hold values for both alleles
  followhap<-rep(list(matrix(NA,nrow=0,ncol=3)),2)
  names(followhap)<-c("Allele_Sire","Allele_Dam")
  # Get the sex of the individual from gene-drop pedigree
  ID_Sex = genedropout@pedigree[genedropout@pedigree[,'ID']==ID, 'Sex']
  
  # For each allele at locus, 1 - sire, 2 - dam
  for (s in 1:2){
    # Set up intial information
    newDes<-cbind(ID,ID_Sex,s)
    gen = 0
    # Loop until break ( there are no more descendants with allele)
    while (TRUE){
      # Keep track of how many generations of descendants with the allele there are
      gen = gen +1
      # Go through each individual in list and apply the findoffhap function
      newDes<-lapply(split(newDes,1:nrow(newDes)),
                     function(x)findoffhap(x[1],loci,x[3],genedropout=genedropout))
      # Join results together neatly
      newDes<-do.call(rbind,newDes)
      # If there are no decendants stop
      if(is.null(newDes)){break}
      # Use generation information as row names
      row.names(newDes) = rep(paste0("Gen_", gen), nrow(newDes))
      followhap[[s]]<-rbind(followhap[[s]],newDes)
    }}
  return(new('hap_track_object', 
          ID=ID, Locus=loci, Allele_Sire=followhap[[1]], Allele_Dam=followhap[[2]]))
  }

### Function extends pedfollowhapsing to allow multiple ID and loci to be passed

pedfollowhap<-function(ID,loci,genedropout){
  IDLoci<-matrix(c(rep(ID,each=length(loci)),rep(loci,length(ID))), ncol=2)
  happed_all<-lapply(split(IDLoci,1:nrow(IDLoci)),function(x) pedfollowhapsing(x[[1]],x[[2]],genedropout))
  # Name list based on ID and loci info
  names(happed_all)<-paste('ID',IDLoci[,1],'LOCI',IDLoci[,2],sep='_')
  if (length(happed_all) == 1){
    return(happed_all[[1]])
  }else{
    return(happed_all)}}