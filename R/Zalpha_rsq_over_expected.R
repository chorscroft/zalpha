

#' Runs the Zalpha function on the r squared values over the expected r squared values for the region
#'
#' Returns a \eqn{Z_{\alpha}^{r^2/E[r^2]}}{Zalpha} value for each SNP location supplied to the function, based on
#' the expected \eqn{r^2} values given an LD profile and cM distances.
#' For more information about the \eqn{Z_{\alpha}^{r^2/E[r^2]}}{Zalpha} statistic please see Jacobs (2016).
#' The \eqn{Z_{\alpha}^{r^2/E[r^2]}}{Zalpha} statistic is defined as:
#' \deqn{{Z_{\alpha}^{r^2/E[r^2]}}=\frac{{|L| \choose 2}^{-1}\sum_{i,j \in L}r^2_{i,j}/E[r^2_{i,j}] + {|R| \choose 2}^{-1}\sum_{i,j \in R}r^2_{i,j}/E[r^2_{i,j}]}{2}}
#' where \code{|L|} and \code{|R|} are the number of SNPs to the left and right of the current locus within the given window \code{ws}, \eqn{r^2}{r^2} is equal to
#' the squared correlation between a pair of SNPs, and \eqn{E[r^2]}{E[r^2]} is equal to the expected squared correlation between a pair of SNPs, given an LD profile.
#'
#' The LD profile describes the expected correlation between SNPs at a given cM distance, generated using simulations or
#' real data. Care should be taken to utilise an LD profile which is representative of the population in question. The LD
#' profile should consist of evenly-sized bins of distances (for example 0.00001 cM per bin), where the value given is the (inclusive) lower
#' bound of the bin.
#'
#' @importFrom stats cor
#'
#' @param pos A numeric vector of SNP locations
#' @param x A matrix of SNP values. Columns represent chromosomes, rows are SNP locations. Hence, the number of rows should equal the length of the \code{pos} vector. SNPs should all be biallelic.
#' @param cM A mnumeric vector of cM distances. This should be the same length as \code{pos}.
#' @param ws The window size which the \eqn{Z_{\alpha}^{r^2/E[r^2]}}{Zalpha} statistic will be calculated over. This should be on the same scale as the \code{pos} vector.
#' @param LDprofile_cM_bins A numeric vector containing the lower bound of the bins used in the LD profile. These should be of equal size.
#' @param LDprofile_rsq A numeric vector containing the expected \eqn{r^2}{r^2} values for the corresponding bin in the LD profile. Must be between 0 and 1.
#' @param minRandL Minimum number of SNPs in each set R and L for the statistic to be calculated. Default is 4.
#' @param minRL Minimum value for the product of the set sizes for R and L. Default is 25.
#' @param X Optional. Specify a region of the chromosome to calculate \eqn{Z_{\alpha}^{r^2/E[r^2]}}{Zalpha} for in the format \code{c(startposition, endposition)}. The start position and the end position should be within the extremes of the positions given in the \code{pos} vector. If not supplied, the function will calculate \eqn{Z_{\alpha}^{r^2/E[r^2]}}{Zalpha} for every SNP in the \code{pos} vector.
#'
#' @return A list containing the SNP positions and the \eqn{Z_{\alpha}^{E[r^2]}}{Zalpha} values for those SNPs
#' @references Jacobs, G.S., T.J. Sluckin, and T. Kivisild, \emph{Refining the Use of Linkage Disequilibrium as a Robust Signature of Selective Sweeps.} Genetics, 2016. \strong{203}(4): p. 1807
#' @export
Zalpha_rsq_over_expected<-function(pos, x, cM, ws, LDprofile_cM_bins, LDprofile_rsq, minRandL = 4, minRL = 25, X = NULL){

  #Check things are in the correct format

  #Check pos is a numeric vector
  if (is.numeric(pos) ==FALSE || is.vector(pos)==FALSE){
    stop("pos must be a numeric vector")
  }
  #Check x is a matrix
  if (is.matrix(x)==FALSE){
    stop("x must be a matrix")
  }
  #Check x has rows equal to the length of pos
  if (length(pos) != nrow(x)){
    stop("The number of rows in x must equal the number of SNP locations given in pos")
  }
  #Check SNPs are all biallelic
  if(sum(apply(x,1,function(x){length(unique(x))}) != 2)>0){
    stop("SNPs must all be biallelic")
  }
  #Check cM is a numeric vector
  if (is.numeric(cM) ==FALSE || is.vector(cM)==FALSE){
    stop("cM must be a numeric vector")
  }
  #Check cM is the same length as pos
  if (length(pos) != length(cM)){
    stop("The number of values in cM must equal the number of SNP locations given in pos")
  }
  #Check windowsize is a number greater than 0
  if(is.numeric(ws) ==FALSE || ws <= 0){
    stop("ws must be a number greater than 0")
  }
  #Check LDprofile_cM_bins is a numeric vector
  if (is.numeric(LDprofile_cM_bins) ==FALSE || is.vector(LDprofile_cM_bins)==FALSE){
    stop("LDprofile_cM_bins must be a numeric vector")
  }
  #Get bin size from LDprofile_cM_bins
  bin_size<-LDprofile_cM_bins[2]-LDprofile_cM_bins[1]
  #Check LDprofile_cM_bins are of equal size
  if (isTRUE(all.equal(diff(LDprofile_cM_bins),rep(bin_size,length(LDprofile_cM_bins)-1)))==FALSE){
    stop("LDprofile_cM_bins must be of equal size")
  }
  #Check LDprofile_rsq is a numeric vector
  if (is.numeric(LDprofile_rsq) ==FALSE || is.vector(LDprofile_rsq)==FALSE){
    stop("LDprofile_rsq must be a numeric vector")
  }
  #Check values of LDprofile_rsq are between 0 and 1
  if (sum(LDprofile_rsq<0 | LDprofile_rsq>1)>0){
    stop("Values stored in LDprofile_rsq must be between 0 and 1")
  }
  #Check that the LDprofile vectors are the same length
  if (length(LDprofile_cM_bins) != length(LDprofile_rsq)){
    stop("LDprofile_rsq must contain the same number of values as there are bins given in LDprofile_cM_bins")
  }
  #Check minRandL is 0 or greater
  if(is.numeric(minRandL) ==FALSE || minRandL < 0){
    stop("minRandL must be a number greater than or equal to 0")
  }
  #Check minRL is 0 or greater
  if(is.numeric(minRL) ==FALSE || minRL < 0){
    stop("minRL must be a number greater than or equal to 0")
  }
  #If X is specified, check it is in the correct format
  if (is.null(X)==FALSE){
    if(is.numeric(X)==FALSE || is.vector(X)==FALSE){
      stop("X should be a numeric vector of length 2 e.g. c(100,200)")
    } else {
      if (length(X) != 2){
        stop("X should be a numeric vector of length 2 e.g. c(100,200)")
      } else {
        # X is in the correct format
        # Check that X will actually return a result (i.e. that the region specied by X overlaps with pos)
        if ((length(pos[pos>=X[1] & pos <= X[2]])>0) == FALSE){
          stop("The region specified by X is outside the region contained in the pos vector")
        }
      }
    }
  } else {
    # Set X equal to the extremes of pos
    X<-c(pos[1],pos[length(pos)])
  }

  #Change matrix x to numeric if it isn't already
  if (is.numeric(x)==FALSE){
    x<-matrix(as.numeric(factor(x)),nrow=dim(x)[1])
  }

  # Set up output list
  outputLength<-length(pos[pos>=X[1] & pos <= X[2]])
  outputList<-list(position=pos[pos>=X[1] & pos <= X[2]],Zalpha_rsq_over_expected=rep(NA,outputLength))


  # Loop over each position in the output list and calculate the expected Zalpha
  for (i in 1:outputLength){

    # Current physical position in chromosome
    currentPos<-outputList$position[i]

    ## check L, R and LR
    noL <- length(pos[pos>=currentPos-ws/2 & pos < currentPos]) ## Number of SNPs to the left of the current SNP
    noR <- length(pos[pos<=currentPos+ws/2 & pos > currentPos]) ## Number of SNPs to the right of the current SNP
    if  (noL < minRandL || noR < minRandL || noL*noR < minRL){
      #NA
      outputList$Zalpha_rsq_over_expected[i]<-NA
    } else {
      ##Left
      # Find cM distances between each SNP in L and round to bin size
      bins<-sapply(lower_triangle(outer(cM[pos>=currentPos-ws/2 & pos < currentPos],cM[pos>=currentPos-ws/2 & pos < currentPos],"-")),assign_bins,bin_size=bin_size)
      Lrsq<- lower_triangle(cor(t(x[pos>=currentPos-ws/2 & pos < currentPos,]))^2)
      LrsqExp<-merge(data.frame(bins,Lrsq),data.frame(LDprofile_cM_bins,LDprofile_rsq),by.x="bins",by.y="LDprofile_cM_bins",all.x=TRUE,sort=FALSE)
      LrsqSum<-sum(LrsqExp$Lrsq/LrsqExp$LDprofile_rsq)
      ##Right
      bins<-sapply(lower_triangle(outer(cM[pos<=currentPos+ws/2 & pos > currentPos],cM[pos<=currentPos+ws/2 & pos > currentPos],"-")),assign_bins,bin_size=bin_size)
      Rrsq<-lower_triangle(cor(t(x[pos<=currentPos+ws/2 & pos > currentPos,]))^2)
      RrsqExp<-merge(data.frame(bins,Rrsq),data.frame(LDprofile_cM_bins,LDprofile_rsq),by.x="bins",by.y="LDprofile_cM_bins",all.x=TRUE,sort=FALSE)
      RrsqSum<-sum(RrsqExp$Rrsq/RrsqExp$LDprofile_rsq)

      outputList$Zalpha_rsq_over_expected[i]<-(LrsqSum/choose(noL,2)+RrsqSum/choose(noR,2))/2
    }
  }
  if (sum(is.na(outputList$Zalpha_rsq_over_expected))==outputLength){
    warning("No Zalpha_rsq_over_expected values were calculated, try reducing minRandL and minRL or increasing the window size")
  }
  return(outputList)
}

