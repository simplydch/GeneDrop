#####                                       #####
##### Functions required to fix pedigree    #####
##### prior to performing gene-dropping     #####
#####                                       #####


#' Function to check  pedigree columns
#'
#' This function checks if ID, Sire and Dam Columns,
#' and any other required columns are present in the pedigree
#' @param peduse A dataframe or matrix. A Pedigree
#' @keywords pedigree structure check
#' @export
#' @examples


ped_col_present<-function(peduse, check_column=c()){
  if (!all(c('ID', 'Dam','Sire', check_column) %in% colnames(peduse))){
      add_comma<-ifelse(length(check_column>0), ', ','')
      stop(paste0("Required columns (ID, Sire, Dam",add_comma,paste0(check_column, collapse=', '),
                  ") are not all present in pedigree. Try using fix_pedigree"), call. = FALSE)
}}



#' Function to check pedigree order
#'
#' This function checks if sires and dams are appearing after offspring
#' and returns TRUE or FALSE.
#' @param peduse A dataframe or matrix. A pedigree containing 'ID', 'Sire' and 'Dam' columns
#' @keywords pedigree sorting order
#' @export
#' @examples


test_order<-function(peduse){
  any(match(peduse[, c("Sire")], peduse[, c("ID")]) > match(peduse[, c("ID")],peduse[, c("ID")]), na.rm = T) ||
    any(match(peduse[, c("Dam")],peduse[, c("ID")]) > match(peduse[, c("ID")],peduse[, c("ID")]), na.rm = T)
}


#' Function to add missing parents
#'
#' This function adds indivdiuals to a pedigree that appear as sires or dams but
#' not in the ID column
#' @param peduse A dataframe or matrix. A pedigree containing 'ID', 'Sire' and 'Dam' columns
#' @param Sex Indicates if a Sex Column is present in provided pedigree. Defaults to TRUE
#' @param Cohort Indicates if a Cohort Column is present in provided pedigree. Defaults to TRUE
#' @keywords pedigree complete fix
#' @export
#' @examples

add_missing_parents<-function(peduse, sex=TRUE, cohort=TRUE){

  names_keep<-c('ID','Sire','Dam','Sex','Cohort')
  ### Find any sire IDs that appear as parents but not in ID column
  # Ignore 0's which indicate foudners
  sire_not_id<-unique(peduse[,'Sire'][is.na(match(peduse[,'Sire'],peduse[,'ID']))])
  sire_not_id<-sire_not_id[!sire_not_id==0]

  ### If there are no missing sires create empty matrix
  # else set-up matrix with sire info
  if(length(sire_not_id)<1){
    matrix(nrow=0,ncol=length(names_keep),dimnames=list(NULL,names_keep))
  }else{
    sire_not_id<-cbind(sire_not_id,0,0,1,NA)
    colnames(sire_not_id)<-names_keep
  }
  ### Find any dam IDs that appear as parents but not in ID column
  # Ignore 0's which indicate foudners
  dam_not_id<-unique(peduse[,'Dam'][is.na(match(peduse[,'Dam'],peduse[,'ID']))])
  dam_not_id<-dam_not_id[!dam_not_id==0]

  ### If there are no missing dams create empty matrix
  # else set-up matrix with dam info
  if(length(dam_not_id)<1){
    matrix(nrow=0,ncol=length(names_keep),dimnames=list(NULL,names_keep))
  }else{
    dam_not_id<-cbind(dam_not_id,0,0,2,NA)
    colnames(dam_not_id)<-names_keep
  }

  ### Create pedigree with added individuals and desired columns
  pedout<-rbind(peduse,sire_not_id,dam_not_id)
  if(is.logical(sex) && sex==FALSE){pedout<-pedout[,!colnames(pedout) %in% 'Sex']}
  if(is.logical(cohort) && cohort==FALSE){pedout<-pedout[,!colnames(pedout) %in% 'Cohort']}
  return(pedout)
}




#' Function to calculate pedigree depth
#'
#' This function calculates pedigree depth each individual in a pedigree by finding the
#' number to descendant generations
#' @param peduse A dataframe or matrix. A pedigree with ID, Sire and Dam columns
#' @keywords pedigree complete fix
#' @export
#' @examples

calc_ped_depth<-function(peduse){
  ped_col_present(peduse) # check correct columns are present
  sire_ref<-match(peduse[,'Sire'],peduse[,'ID'],nomatch=0)
  dam_ref<-match(peduse[,'Dam'],peduse[,'ID'],nomatch=0)
  sapply(matrix(1:nrow(peduse)), function(x){
    d<--1
    n<-x
    # Increase depth (d) until both sire and dam are founders
    # i.e have have both parents as 0
    repeat{
      n<-c(sire_ref[n],dam_ref[n])
      peduse[c(n),c('ID')]
      d<-d+1
      if(sum(n)==0){return(d)}
  }})
}


#' Function to ensure pedigree is correctly structured for gene dropping
#'
#' This function checks the pedigree for potentential problems in structure
#' it reorganises the individuals if neccessary and renames columns to ensure the
#' pedigree is in the correct format for gene-dropping.  Sex and Cohort can be
#' included as columns in the pedigree or provided as seperate vectors.
#' It adds in any individuals that are missing and checks sex coding is correct
#' and consistent. It returns the pedigree as a matrix
#'
#' @param peduse A dataframe or matrix. A pedigree with ID, Sire and Dam columns.
#' @param sex logical or vector. TRUE if sex column is included in pedigree, otherwise
#' a vector containing the sex of each individual. Defaults to TRUE
#' @param cohort logical or vector. TRUE if cohort column is included in pedigree, otherwise
#' a vector containing the cohort of each individual. Defaults to TRUE
#' @keywords pedigree complete fix structure
#' @export
#' @examples


fix_pedigree<-function(peduse, sex=TRUE, cohort=TRUE){

  #TODO Allow function to work if sex is absent?
  if(is.logical(sex) && sex==FALSE){stop('Sex should be TRUE or a vector')}
  if(is.logical(cohort) && cohort==FALSE){stop('Cohort should be TRUE or a vector')}

  # If passed dataframe convert to matrix (faster to work with)
  peduse<-as.matrix(peduse)

  # Check that cohort and sex vectors are correct structure if passed
  if(!is.logical(cohort) & !length(cohort)==nrow(peduse)){stop('Provided cohort vector is wrong length or structure')}
  if(!is.logical(sex) & !length(sex)==nrow(peduse)){stop('Provided cohort vector is wrong length or structure')}

  # Check that sex column exists in pedigree if Sex is set to TRUE
  if(is.logical(sex) && sex==TRUE){
    sex_col<-grep("sex",ignore.case=TRUE,colnames(peduse))
    if(!length(sex_col)==1){stop("Sex column not in provided pedigree but Sex==TRUE")}
  }else{sex_col<-NULL}

  # Check that cohort column exists in pedigree if Cohort is set to TRUE
  if(is.logical(cohort) && cohort==TRUE){
    cohort_col<-grep("cohort|birthyear",ignore.case=TRUE,colnames(peduse))
    if(!length(cohort_col)==1){stop("Cohort column not in provided pedigree but Cohort==TRUE")}
  }else{cohort_col<-NULL}


  # Find id, sire and dam columns
  sire_col<-grep("sire|dad|father",ignore.case=TRUE,colnames(peduse))
  dam_col<-grep("mom|mum|dam|mother",ignore.case=TRUE,colnames(peduse))
  id_col<-grep("id",ignore.case=TRUE,colnames(peduse))

  # Check that ID, Sire and Dam columns can be identified
  if(any(!length(sire_col)==1,!length(dam_col)==1,!length(id_col)==1)){
    stop("Columns ID, Sire, Dam don't appear in pedigree")}

  # Take relevant column names for the provided pedigree
  oldnames<-paste(colnames(peduse)[sort(c(sire_col,dam_col,id_col,sex_col,cohort_col))], collapse=',')

  # Make sure ID Sire, Dam, columns are in the correct order
  peduse<-peduse[,c(id_col,sire_col,dam_col,sex_col,cohort_col)]

  # Add Sex and Cohort columns if required
  if(!is.logical(cohort)){peduse<-cbind(peduse,cohort);cohort_col<-4}
  if(!is.logical(sex)){peduse<-cbind(peduse,sex);sex_col<-5}

  # Re-name Columns
  colnames(peduse)<-c("ID","Sire","Dam","Sex","Cohort")[sapply(list(id_col,sire_col,dam_col,sex_col,cohort_col),function(x) !is.null(x))]
  newnames<-paste(colnames(peduse), collapse=',')

  # Warn if columns have been renamed or re-ordered
  if(!oldnames==newnames){
    cat(paste0("Note: Columns ",oldnames, " reordered and/or renamed as ", newnames,"\n"))}

  # Add any parents that are missing
  peduse<-add_missing_parents(peduse, sex=sex, cohort=cohort)

  # Check Sex #


  # Get parent references
  sire_ref<-match(peduse[,'Sire'],peduse[,'ID'],nomatch=0)
  dam_ref<-match(peduse[,'Dam'],peduse[,'ID'],nomatch=0)

  # Find original sex codings based on parent info
  org_sire_sex_code = find_most_common(peduse[unique(sire_ref), 'Sex'])[['most_common']]
  if (length(org_sire_sex_code) != 1){
    stop(paste0("Male sex code can't be determined. These codes are equally common: ",
                paste(org_sire_sex_code, collapse=',')))
  }
  org_dam_sex_code = find_most_common(peduse[unique(dam_ref), 'Sex'])[['most_common']]
  if (length(org_sire_sex_code) != 1){
    stop(paste0("Female sex code can't be determined. These codes are equally common: ",
                paste(org_dam_sex_code, collapse=',')))
  }

  # Recode to sire = 1 and dam = 2
  if (org_sire_sex_code != 1) {peduse[peduse[, 'Sex']== org_sire_sex_code, 'Sex']<- -98}
  if (org_dam_sex_code != 2) {peduse[peduse[, 'Sex']== org_dam_sex_code, 'Sex']<- -99}

  peduse[peduse[, 'Sex'] == -98, 'Sex'] <- 1
  peduse[peduse[, 'Sex'] == -99, 'Sex'] <- 2

  # Print message if there has been recoding
  if (org_sire_sex_code != 1 | org_dam_sex_code != 2){
  message(paste0("Note: Sex of individuals have been recoded. Sires now 1 (were ", org_sire_sex_code,
                 "). Dams now 2 (were ", org_dam_sex_code, ")"))}

  # Get Sex Coding
  sex_codes = peduse[, 'Sex']

  # Ensure that sires and dams have correct sex code
  peduse[unique(dam_ref), 'Sex'] <- 2
  peduse[unique(sire_ref), 'Sex'] <- 1

  # Print message if any changes have been made i.e dams with male sex code
  if (any(peduse[,'Sex'] != sex_codes)){
    message(paste0("Note: Sex of some individuals have been changed due to parent data"))
  }

  # If pedigree is already in correct order return it
  if(!test_order(peduse)){return(peduse)}

  # Otherwise re-order pedigiree

  order2<-calc_ped_depth(peduse)

  # Warn if there are indivdiuals that do not have a cohort
  if(!is.null(cohort_col) && any(is.na(peduse[,'Cohort']))){message("Note: Individuals added do not have a cohort. This may be required for gene-dropping. Using est_cohort might help")}
  return(peduse[order(order2),])
  #data.frame(peduse)
}


#' Function to estimate missing cohorts
#'
#' A function that estimates cohorts for individuals without a known cohort.
#' It uses only pedigree and cohort information and should be secondary to more
#' accurate predictions.  If offspring information is present this is used, based on a provided
#' to reproduction.  When offspring data is absent a mean of sibling cohorts is used.
#' The function loops using the new estimates until no further cohorts can be
#' estimated.
#'
#' @param peduse A dataframe or matrix. A pedigree with ID, Sire and Dam columns
#' @param tt_rep integer. Time to reproduction.  A value of the length of time between
#' birth and reproduction. Defaults to 2
#' @keywords estimate cohort
#' @export
#' @examples


est_cohort<-function(peduse, tt_rep=2){

  #TODO Add seperate time to reproduction for males and females?
  #     Will this make the function less general?

  # Check columns are present
  ped_col_present(peduse, c('Cohort'))

  # Get sire and dam row numbers
  sire_ref<-match(peduse[,'Sire'],peduse[,'ID'],nomatch=0)
  dam_ref<-match(peduse[,'Dam'],peduse[,'ID'],nomatch=0)

  # Get cohorts
  cohorts<-peduse[,'Cohort']

  # Loop until no new cohorts are being estimated
  while(any(is.na(cohorts))){
    # initial number of cohorts missing (x)
    x<-length(cohorts[is.na(cohorts)])

    cohorts<-sapply(matrix(1:length(cohorts)), function(x){
      # if cohort isn't NA don't do anything
      if(!is.na(cohorts[x])){return(cohorts[x])}

      ##  Find average sibling cohort
      # Row references for siblings from sire
      sib_ref_s <- if(peduse[x,'Sire']==0){NA}else{
        cohorts[peduse[x,'Sire'] == peduse[,'Sire']]}

      # Row references for siblings from dam
      sib_ref_d <- if(peduse[x,'Dam']==0){NA}else{
        cohorts[peduse[x,'Dam'] == peduse[,'Dam']]}

      # Remove NA's and calculate mean across siblings from both sire and dam
      sib_by <- c(sib_ref_s, sib_ref_d)
      sib_by <- sib_by[!is.na(sib_by)]
      est_bs_sib <- ifelse(length(sib_by)>0,round(mean(sib_by,0)),NA)


      ##  Estimate based on first offspring
      # Find offspring cohort as sire and then dam
      by_off_s <- cohorts[peduse[x,'ID'] == peduse[,'Sire']]
      by_off_d <- cohorts[peduse[x,'ID'] == peduse[,'Dam']]

      by_off <- c(by_off_s, by_off_d)
      by_off <- by_off[!is.na(by_off)]
      # Estimated based on first offspring cohort minus TTR

      est_b_off <- ifelse(length(by_off>0),min(by_off)-tt_rep,NA)

      # Use prediction based on offspring as preferences, otherwise used sibling mean
      ifelse(!is.na(est_bs_sib) & est_bs_sib==est_b_off | is.na(est_b_off), est_bs_sib, est_b_off)}) # End of sapply
    # If there are missing cohorts that can't be estimated produce warning
    if(length(cohorts[is.na(cohorts)])==x){message(paste0("Warning: Cohort can not be estimated for ", x,
                                                          " individuals. Remove, or manually add cohort prior to gene-dropping."));break}
  }

  cohorts
}



#' Function to complete pedigree links
#'
#' This fills in missing pedigree links in a random fashion using plausable
#' individuals. It for individuals with missing sires or dams it selects a
#' parent that could have been alive within a suitable time period
#' @param peduse A dataframe or matrix. A pedigree with ID, Sire and Dam, Sex columns
#' @param founders A vector of founder IDs matching those in the pedigree
#' @param rep_years_sire vector. A vector of length 2, the first value is the time until a male is reproductively
#' active and the second value the time a male is no longer reproductively active. Default: c(2,10)
#' @param rep_years_dam vector. A vector of length 2, the first value is the time until a female is reproductively
#' active and the second value the time a female is no longer reproductively active. Default: c(2,10)
#' @param cohort An optional vector, the same length of the pedigree containing cohort information which can be
#' provided if there is not a cohort column in the pedigree
#' @keywords complete pedigree
#' @export
#' @examples

###

complete_ped_links<-function(peduse, founders,
                             rep_years_sire=c(2,10),
                             rep_years_dam=c(2,10),
                             cohorts=NULL){

  # Check columns are present
  ped_col_present(peduse, c('Sex'))

  # Check founders are present in pedigree
  if(any(is.na(match(founders, pedigree[,'ID'])))){
    stop(paste0("Not all founder IDs are present in the pedigree"), call. = FALSE)
  }

  # Set-up founder reference

  founder_ref <- as.logical(match(pedigree[,'ID'], c(founders,'aaa'), nomatch = 0))


  # If cohorts aren't provided seperately extract them from pedigree
  if(is.null(cohorts)){
    ped_col_present(peduse, c('Sex', 'Cohort'))
    cohorts<-peduse[ , 'Cohort']}

  # convert any NA's in Dam or Sire column to 0's
  peduse[is.na(peduse[ , 'Dam']), 'Dam'] <- 0
  peduse[is.na(peduse[ , 'Sire']), 'Sire'] <- 0

  # Set up matrices to contain new sires and dam
  new_sire <- matrix(NA,nrow=nrow(peduse), ncol=1)
  new_dam <- matrix(NA,nrow=nrow(peduse), ncol=1)

  # Add sire and dam info, this approach allows the function to be passed either
  # dataframe or matrix
  new_sire[ ,1] <- peduse[ , 'Sire']
  new_dam[ ,1] <- peduse[ , 'Dam']

  # Set founders as 0
  new_sire[founder_ref, ] <- 0
  new_dam[founder_ref, ] <- 0

  # Find max and min possible cohorts of parents based on reproductive ages
  range_dam_high <- peduse[,'Cohort'] - rep_years_dam[1]
  range_dam_low <- peduse[,'Cohort'] - rep_years_dam[2]

  range_sire_high <- peduse[,'Cohort'] - rep_years_sire[1]
  range_sire_low <- peduse[,'Cohort'] - rep_years_sire[2]


  # Get Sire and Dam row references
  dam_ref <- which(peduse[,'Sex'] == 2)
  sire_ref <- which(peduse[,'Sex'] == 1)

  # Find individuals that need estimated dams and sires (no parent but not founder)
  dam_to_est <- which(!founder_ref & new_dam == 0)
  sire_to_est <- which(!founder_ref & new_sire == 0)

  # Assign a random dam born within calculated window to those without known parent
  new_dam[dam_to_est,] <- sapply(matrix(dam_to_est), function(x){
    t <- which(cohorts >= range_dam_low[x] & cohorts <= range_dam_high[x])
    t2 <- t[t %in% dam_ref]
    peduse[t[sample.int(length(t2),1)], c('ID')]
  })

  # Assign a random sire born within calculated window to those withouth known parent
  new_sire[sire_to_est,]<-sapply(matrix(sire_to_est), function(x){
    t <- which(cohorts>=range_sire_low[x] & cohorts<=range_sire_high[x])
    t2 <- t[t %in% sire_ref]
    peduse[t[sample.int(length(t2),1)],c('ID')]
  })

  # return pedigree in expected format
  peduse <- matrix(c(peduse[,'ID'], new_sire, new_dam, peduse[,'Sex'], peduse[,'Cohort']), ncol=5,
                    dimnames=(list(NULL,c('ID','Sire','Dam','Sex','Cohort'))))

  ## Ensure pedigree is in correct order etc
  fix_pedigree(peduse)
}

