
#' Creates an LD profile
#'
#' An LD (linkage disequilibrium) profile is a look-up table that tell you the expected correlation between SNPs given the genetic distance between them.
#'
#' In the output, bins represent lower bounds. The first bin contains pairs where the genetic distance is greater than or equal to 0 and less than \code{bin_size}. The final bin contains pairs where the genetic distance is greater than or equal to \code{max_dist}-\code{bin_size} and less than \code{max_dist}.
#' If the \code{max_dist} is not an increment of \code{bin_size}, it will be adjusted to the next highest increment.The maximum bin will be the bin that \code{max_dist} falls into. For example, if the \code{max_dist} is given as 4.5 and the \code{bin_size} is 1, the final bin will be 4.\cr
#' By default, Beta parameters are not calculated. To calcualte Beta parameters, needed for the \code{\link{Zalpha_BetaCDF}} and \code{\link{Zbeta_BetaCDF}} statistics, \code{beta_params} should be set to TRUE and the package \code{fitdistrplus} must be installed.
#'
#' @param dist A numeric vector containing genetic distances.
#' @param x A matrix of SNP values. Columns represent chromosomes; rows are SNP locations. Hence, the number of rows should equal the length of the \code{dist} vector. SNPs should all be biallelic.
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
#'
#' @export
#' @seealso \code{\link{Zalpha_expected}} \code{\link{Zalpha_rsq_over_expected}} \code{\link{Zalpha_log_rsq_over_expected}} \code{\link{Zalpha_Zscore}} \code{\link{Zalpha_BetaCDF}} \code{\link{Zbeta_expected}} \code{\link{Zbeta_rsq_over_expected}} \code{\link{Zbeta_log_rsq_over_expected}} \code{\link{Zbeta_Zscore}} \code{\link{Zbeta_BetaCDF}} \code{\link{Zalpha_all}}
#'
create_LDprofile<-function(dist,x,bin_size,max_dist=NULL,beta_params=FALSE){
  #Checks
  #Check dist is vector
  if (is.numeric(dist) ==FALSE || is.vector(dist)==FALSE){
    stop("dist must be a numeric vector")
  }
  #Check x is a matrix
  if (is.matrix(x)==FALSE){
    stop("x must be a matrix")
  }
  #Check x has rows equal to the length of dist
  if (length(dist) != nrow(x)){
    stop("The number of rows in x must equal the number of SNP genetic distances given in dist")
  }
  #Check SNPs are all biallelic
  if (sum(apply(x,1,function(x){length(na.omit(unique(x)))}) != 2)>0){
    stop("SNPs must all be biallelic")
  }
  #Change matrix x to numeric if it isn't already
  if (is.numeric(x)==FALSE){
    x<-matrix(as.numeric(factor(x)),nrow=dim(x)[1])
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
    max_dist<-dist[length(dist)]-dist[1]
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

  #Find the differences in genetic distances between pairs of SNPs
  diffs<-lower_triangle(outer(dist,dist,"-"))

  #Find the rsquared value between pairs of SNPs
  rsq<-lower_triangle(cor(t(x),use="pairwise.complete.obs")^2)

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
