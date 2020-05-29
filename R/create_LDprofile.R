
#' Creates an LD profile
#'
#' An LD (linkage disequilibrium) profile is a look-up table containing the expected correlation between SNPs given the genetic distance between them.
#'
#' The input for \code{dist} and \code{x} can be lists. This is to allow multiple datasets to be used in the creation of the LD profile. For example, using all 22 autosomes from the human genome would involve 22 different distance vectors and SNP matrices.
#' Both lists should be the same length and should correspond exactly to eachother (i.e. the distances in each element of dist should go with the SNPs in the same element of x)
#'
#' In the output, bins represent lower bounds. The first bin contains pairs where the genetic distance is greater than or equal to 0 and less than \code{bin_size}. The final bin contains pairs where the genetic distance is greater than or equal to \code{max_dist}-\code{bin_size} and less than \code{max_dist}.
#' If the \code{max_dist} is not an increment of \code{bin_size}, it will be adjusted to the next highest increment.The maximum bin will be the bin that \code{max_dist} falls into. For example, if the \code{max_dist} is given as 4.5 and the \code{bin_size} is 1, the final bin will be 4.\cr
#' By default, Beta parameters are not calculated. To calcualte Beta parameters, needed for the \code{\link{Zalpha_BetaCDF}} and \code{\link{Zbeta_BetaCDF}} statistics, \code{beta_params} should be set to TRUE and the package \code{fitdistrplus} must be installed.
#'
#' @param dist A numeric vector, or a list of numeric vectors, containing genetic distances.
#' @param x A matrix of SNP values, or a list of matrices. Columns represent chromosomes; rows are SNP locations. Hence, the number of rows should equal the length of the \code{dist} vector. SNPs should all be biallelic.
#' @param bin_size The size of each bin, in the same units as \code{dist}.
#' @param max_dist Optional. The maximum genetic distance to be considered. If this is not supplied, it will default to the maximum distance in the \code{dist} vector.
#' @param beta_params Optional. Beta parameters are calculated if this is set to TRUE. Default is FALSE.
#'
#' @return A data frame containing an LD profile that can be used by other statistics in this package.
#' @references Jacobs, G.S., T.J. Sluckin, and T. Kivisild, \emph{Refining the Use of Linkage Disequilibrium as a Robust Signature of Selective Sweeps.} Genetics, 2016. \strong{203}(4): p. 1807
#' @examples
#' ## load the snps example dataset
#' data(snps)
#' ## Create an LD profile using this data
#' create_LDprofile(snps$distances,as.matrix(snps[,3:12]),0.001)
#' ## To get the Beta distribution parameter estimates, the fitdistrplus package is required
#' if (requireNamespace("fitdistrplus", quietly = TRUE)==TRUE) {
#'   create_LDprofile(snps$distances,as.matrix(snps[,3:12]),0.001,beta_params=TRUE)
#' }
#'
#'
#' @export
#' @seealso \code{\link{Zalpha_expected}} \code{\link{Zalpha_rsq_over_expected}} \code{\link{Zalpha_log_rsq_over_expected}} \code{\link{Zalpha_Zscore}} \code{\link{Zalpha_BetaCDF}} \code{\link{Zbeta_expected}} \code{\link{Zbeta_rsq_over_expected}} \code{\link{Zbeta_log_rsq_over_expected}} \code{\link{Zbeta_Zscore}} \code{\link{Zbeta_BetaCDF}} \code{\link{Zalpha_all}}
#'
create_LDprofile<-function(dist,x,bin_size,max_dist=NULL,beta_params=FALSE){

  #Changes dist into a list if it is not one already
  if (is.list(dist)==FALSE){
    dist<-list(dist)
  }

  #Changes x into a list if it is not one already
  if (is.list(x)==FALSE){
    x<-list(x)
  }

  #Checks
  #Check the dist list and the x list have the same length
  if (length(dist)!=length(x)){
    stop("dist and x should contain the same number of elements")
  }

  #Check for each element in dist and x
  for (el in 1:length(dist)){
    #Check dist is vector
    if (is.numeric(dist[[el]]) ==FALSE || is.vector(dist[[el]])==FALSE){
      stop("dist must be a numeric vector or list of numeric vectors")
    }
    #Check x is a matrix
    if (is.matrix(x[[el]])==FALSE){
      stop("x must be a matrix or list of matrices")
    }
    #Check x has rows equal to the length of dist
    if (length(dist[[el]]) != nrow(x[[el]])){
      stop("The number of rows in x must equal the number of SNP genetic distances given in the corresponding dist")
    }
    #Check SNPs are all biallelic
    if (sum(apply(x[[el]],1,function(x){length(na.omit(unique(x)))}) != 2)>0){
      stop("SNPs must all be biallelic")
    }
    #Change matrix x to numeric if it isn't already
    if (is.numeric(x[[el]])==FALSE){
      x[[el]]<-matrix(as.numeric(factor(x[[el]])),nrow=dim(x[[el]])[1])
    }
  }

  #Check bin_size is a number
  if (is.numeric(bin_size) ==FALSE || bin_size <= 0){
    stop("bin_size must be a number greater than 0")
  }
  #Check max_dist is a number or NULL
  if (is.null(max_dist)==FALSE){
    if (is.numeric(max_dist) ==FALSE || max_dist <= 0){
      stop("max_dist must be a number greater than 0")
    }
  } else {
    #Set max_dist to the maximum distance in the data if it was not supplied
    max_dist<-max(sapply(dist,function(x){x[length(x)]-x[1]}),na.rm = TRUE)
  }
  #Adjusts the Max_dist value so it is equal to an increment of bin_size if it isn't already
  if(!isTRUE(all.equal(max_dist,assign_bins(bin_size,max_dist)))){
    max_dist<-assign_bins(bin_size,max_dist)+bin_size
  }
  #Check beta_params is logical
  if (is.logical(beta_params)==FALSE){
    stop("beta_params must be TRUE or FALSE")
  }
  #If beta_params is TRUE, check for fitdistrplus package
  if (beta_params==TRUE){
    if (requireNamespace("fitdistrplus", quietly = TRUE)==FALSE) {
      stop("Package \"fitdistrplus\" needed for Beta parameters to be calculated. Please install it.")
    }
  }

  diffs<-NULL
  rsq<-NULL
  #for each element in dist and x, get the differences and rsquared values
  for (el in 1:length(dist)){
    #Find the differences in genetic distances between pairs of SNPs
    diffs<-c(diffs,lower_triangle(outer(dist[[el]],dist[[el]],"-")))

    #Find the rsquared value between pairs of SNPs
    rsq<-c(rsq,lower_triangle(cor(t(x[[el]]),use="pairwise.complete.obs")^2))
  }
  #Filter for just those less than the max genetic distance and filter out missing distances
  rsq<-rsq[diffs<max_dist & is.na(diffs)==FALSE]
  diffs<-diffs[diffs<max_dist & is.na(diffs)==FALSE]

  #Assign diffs to bins
  bins<-assign_bins(bin_size,diffs)

  #Create LDprofile data frame
  LDprofile<-data.frame(bin=seq(0,max_dist-bin_size,bin_size),rsq=NA,sd=NA,Beta_a=NA,Beta_b=NA,n=NA)

  #Loop for each bin (i)
  for (i in 1:nrow(LDprofile)){
    LDprofile$n[i]<-sum(bins==LDprofile$bin[i])
    #If there is at least one pair whose genetic distance falls within the bin, calculate stats
    if (LDprofile$n[i]>0){
      #Get the rsquared values for all pairs in this bin
      temprsq<-rsq[bins==LDprofile$bin[i]]
      #Calculate the mean
      LDprofile$rsq[i]<-mean(temprsq)
      #Calculate the standard deviation
      LDprofile$sd[i]<-sd(temprsq)

      #Calculate Beta distribution parameters if required
      #Do not calculate for bins containing less than two pairs or the standatd deviation is zero
      if (beta_params==TRUE & LDprofile$n[i]>1 & LDprofile$sd[i]>0){
        if (sum(temprsq==1 | temprsq==0)>0){
          #If there are any 0s or 1s adjust the data
          temprsq<-(temprsq*(length(temprsq)-1)+0.5)/length(temprsq)
        }
        #Try to fit the data to a Beta distribution
        betafit<-try(fitdistrplus::fitdist(temprsq,"beta"))
        if (class(betafit) != "try-error"){
          LDprofile$Beta_a[i]<-betafit$estimate[1]
          LDprofile$Beta_b[i]<-betafit$estimate[2]
        } else {
          #If failed to fit, try again using estimated beta parameters to initialise
          startBetaParams<-est_Beta_Params(LDprofile$rsq[i], LDprofile$sd[i]^2)
          betafit<-try(fitdistrplus::fitdist(temprsq,"beta",start=list(shape1=startBetaParams$alpha, shape2=startBetaParams$beta)))
          if (class(betafit) != "try-error"){
            LDprofile$Beta_a[i]<-betafit$estimate[1]
            LDprofile$Beta_b[i]<-betafit$estimate[2]
          } else {
            #If Beta parameters cannot be fitted, return NA
            LDprofile$Beta_a[i]<-NA
            LDprofile$Beta_b[i]<-NA
          }
        }
      }
    }
  }
  return(LDprofile)
}
