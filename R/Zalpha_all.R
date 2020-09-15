
#' Runs all the statistics in the zalpha package
#'
#' Returns every statistic for each SNP location, given the appropriate parameters. See Details for more information.
#'
#' Not all statistics will be returned, depending on the parameters supplied to the function.\cr
#' If \code{x} is not supplied, only \code{\link{Zalpha_expected}}, \code{\link{Zbeta_expected}}, \code{\link{LR}} and \code{\link{L_plus_R}} will be calculated.\cr
#' For any of the statistics which use an expected \eqn{r^2}{r^2} value, the parameters \code{dist}, \code{LDprofile_bins} and \code{LDprofile_rsq} must be supplied.
#' This includes the statistics: \code{\link{Zalpha_expected}}, \code{\link{Zalpha_rsq_over_expected}}, \code{\link{Zalpha_log_rsq_over_expected}}, \code{\link{Zalpha_Zscore}}, \code{\link{Zalpha_BetaCDF}}, \code{\link{Zbeta_expected}}, \code{\link{Zbeta_rsq_over_expected}}, \code{\link{Zbeta_log_rsq_over_expected}}, \code{\link{Zbeta_Zscore}} and \code{\link{Zbeta_BetaCDF}}.
#' \itemize{
#'   \item For \code{\link{Zalpha_Zscore}} and \code{\link{Zbeta_Zscore}} to be calculated, the parameter \code{LDprofile_sd} must also be supplied.
#'   \item For \code{\link{Zalpha_BetaCDF}} and \code{\link{Zbeta_BetaCDF}} to be calculated, the parameters \code{LDprofile_Beta_a} and \code{LDprofile_Beta_b} must also be supplied.
#' }
#' The LD profile describes the expected correlation between SNPs at a given genetic distance, generated using simulations or
#' real data. Care should be taken to utilise an LD profile that is representative of the population in question. The LD
#' profile should consist of evenly sized bins of distances (for example 0.0001 cM per bin), where the value given is the (inclusive) lower
#' bound of the bin. Ideally, an LD profile would be generated using data from a null population with no selection, however one can be generated
#' using this data. See the \code{\link{create_LDprofile}} function for more information on how to create an LD profile.
#' For more information about the statistics, please see Jacobs (2016).
#'
#' @importFrom stats cor pbeta na.omit
#'
#' @param pos A numeric vector of SNP locations
#' @param ws The window size which the statistics will be calculated over. This should be on the same scale as the \code{pos} vector.
#' @param x Optional. A matrix of SNP values. Columns represent chromosomes; rows are SNP locations. Hence, the number of rows should equal the length of the \code{pos} vector. SNPs should all be biallelic.
#' @param dist Optional. A numeric vector of genetic distances (e.g. cM, LDU). This should be the same length as \code{pos}.
#' @param LDprofile_bins Optional. A numeric vector containing the lower bound of the bins used in the LD profile. These should be of equal size.
#' @param LDprofile_rsq Optional. A numeric vector containing the expected \eqn{r^2}{r^2} values for the corresponding bin in the LD profile. Must be between 0 and 1.
#' @param LDprofile_sd Optional. A numeric vector containing the standard deviation of the \eqn{r^2}{r^2} values for the corresponding bin in the LD profile.
#' @param LDprofile_Beta_a Optional. A numeric vector containing the first estimated Beta parameter for the corresponding bin in the LD profile.
#' @param LDprofile_Beta_b Optional. A numeric vector containing the second estimated Beta parameter for the corresponding bin in the LD profile.
#' @param minRandL Minimum number of SNPs in each set R and L for the statistics to be calculated. L is the set of SNPs to the left of the target SNP and R to the right, within the given window size \code{ws}. Default is 4.
#' @param minRL Minimum value for the product of the set sizes for R and L. Default is 25.
#' @param X Optional. Specify a region of the chromosome to calculate the statistics for in the format \code{c(startposition, endposition)}. The start position and the end position should be within the extremes of the positions given in the \code{pos} vector. If not supplied, the function will calculate the statistics for every SNP in the \code{pos} vector.
#'
#' @return A list containing the SNP positions and the statistics for those SNPs
#' @references Jacobs, G.S., T.J. Sluckin, and T. Kivisild, \emph{Refining the Use of Linkage Disequilibrium as a Robust Signature of Selective Sweeps.} Genetics, 2016. \strong{203}(4): p. 1807
#' @examples
#' ## load the snps and LDprofile example datasets
#' data(snps)
#' data(LDprofile)
#' ## run Zalpha_all over all the SNPs with a window size of 3000 bp
#' ## will return all 15 statistics
#' Zalpha_all(snps$bp_positions,3000,as.matrix(snps[,3:12]),snps$cM_distances,
#'  LDprofile$bin,LDprofile$rsq,LDprofile$sd,LDprofile$Beta_a,LDprofile$Beta_b)
#' ## only return results for SNPs between locations 600 and 1500 bp
#' Zalpha_all(snps$bp_positions,3000,as.matrix(snps[,3:12]),snps$cM_distances,
#'  LDprofile$bin,LDprofile$rsq,LDprofile$sd,LDprofile$Beta_a,LDprofile$Beta_b,X=c(600,1500))
#' ## will only return statistics not requiring an LD profile
#'Zalpha_all(snps$bp_positions,3000,as.matrix(snps[,3:12]))
#'
#' @export
#' @seealso \code{\link{Zalpha}}, \code{\link{Zalpha_expected}}, \code{\link{Zalpha_rsq_over_expected}}, \code{\link{Zalpha_log_rsq_over_expected}}, \code{\link{Zalpha_Zscore}}, \code{\link{Zalpha_BetaCDF}}, \code{\link{Zbeta}}, \code{\link{Zbeta_expected}}, \code{\link{Zbeta_rsq_over_expected}}, \code{\link{Zbeta_log_rsq_over_expected}}, \code{\link{Zbeta_Zscore}}, \code{\link{Zbeta_BetaCDF}}, \code{\link{LR}}, \code{\link{L_plus_R}}, \code{\link{create_LDprofile}}.

Zalpha_all <- function(pos, ws, x=NULL, dist=NULL, LDprofile_bins=NULL, LDprofile_rsq=NULL, LDprofile_sd=NULL, LDprofile_Beta_a=NULL, LDprofile_Beta_b=NULL, minRandL = 4, minRL = 25, X = NULL) {
  #Check things are in the correct format

  #Check pos is a numeric vector
  if (is.numeric(pos) ==FALSE || is.vector(pos)==FALSE){
    stop("pos must be a numeric vector")
  }
  #If x is supplied
  if (is.null(x)==FALSE){
    #Check x is a matrix
    if (is.matrix(x)==FALSE){
      stop("x must be a matrix")
    }
    #Check x has rows equal to the length of pos
    if (length(pos) != nrow(x)){
      stop("The number of rows in x must equal the number of SNP locations given in pos")
    }
    #Check SNPs are all biallelic
    if(sum(apply(x,1,function(x){length(na.omit(unique(x)))}) != 2)>0){
      stop("SNPs must all be biallelic")
    }
    #Change matrix x to numeric if it isn't already
    if (is.numeric(x)==FALSE){
      x<-matrix(as.numeric(factor(x)),nrow=dim(x)[1])
    }
  }
  #Check windowsize is a number greater than 0
  if(is.numeric(ws) ==FALSE || ws <= 0){
    stop("ws must be a number greater than 0")
  }
  if (is.null(dist)==FALSE){
  #Check dist is a numeric vector
    if (is.numeric(dist) ==FALSE || is.vector(dist)==FALSE){
      stop("dist must be a numeric vector")
    }
    #Check dist is the same length as pos
    if (length(pos) != length(dist)){
      stop("The number of values in dist must equal the number of SNP locations given in pos")
    }
  }
  if (is.null(LDprofile_bins)==FALSE){
  #Check LDprofile_bins is a numeric vector
    if (is.numeric(LDprofile_bins) ==FALSE || is.vector(LDprofile_bins)==FALSE){
      stop("LDprofile_bins must be a numeric vector")
    }
    #Get bin size from LDprofile_bins
    bin_size<-LDprofile_bins[2]-LDprofile_bins[1]
    #Check LDprofile_bins are of equal size
    if (isTRUE(all.equal(diff(LDprofile_bins),rep(bin_size,length(LDprofile_bins)-1)))==FALSE){
      stop("LDprofile_bins must be of equal size")
    }
  }
  if (is.null(LDprofile_rsq)==FALSE & is.null(LDprofile_bins)==FALSE){
  #Check LDprofile_rsq is a numeric vector
    if (is.numeric(LDprofile_rsq) ==FALSE || is.vector(LDprofile_rsq)==FALSE){
      stop("LDprofile_rsq must be a numeric vector")
    }
    #Check values of LDprofile_rsq are between 0 and 1
    if (sum(LDprofile_rsq<0 | LDprofile_rsq>1)>0){
      stop("Values stored in LDprofile_rsq must be between 0 and 1")
    }
    #Check that the LDprofile vectors are the same length
    if (length(LDprofile_bins) != length(LDprofile_rsq)){
      stop("LDprofile_rsq must contain the same number of values as there are bins given in LDprofile_bins")
    }
  }
  if(is.null(LDprofile_sd)==FALSE & is.null(LDprofile_bins)==FALSE){
    #Check LDprofile_sd is a numeric vector
    if (is.numeric(LDprofile_sd) ==FALSE || is.vector(LDprofile_sd)==FALSE){
      stop("LDprofile_sd must be a numeric vector")
    }
    #Check that the LDprofile vectors are the same length
    if (length(LDprofile_bins) != length(LDprofile_sd)){
      stop("LDprofile_sd must contain the same number of values as there are bins given in LDprofile_bins")
    }
  }
  if(is.null(LDprofile_Beta_a)==FALSE & is.null(LDprofile_Beta_b)==FALSE & is.null(LDprofile_bins)==FALSE){
    #Check LDprofile_Beta_a is a numeric vector
    if (is.numeric(LDprofile_Beta_a) ==FALSE || is.vector(LDprofile_Beta_a)==FALSE){
      stop("LDprofile_Beta_a must be a numeric vector")
    }
    #Check LDprofile_Beta_b is a numeric vector
    if (is.numeric(LDprofile_Beta_b) ==FALSE || is.vector(LDprofile_Beta_b)==FALSE){
      stop("LDprofile_Beta_b must be a numeric vector")
    }
    #Check that the LDprofile vectors are the same length
    if (length(LDprofile_bins) != length(LDprofile_Beta_a)){
      stop("LDprofile_Beta_a must contain the same number of values as there are bins given in LDprofile_bins")
    }
    #Check that the LDprofile vectors are the same length
    if (length(LDprofile_bins) != length(LDprofile_Beta_b)){
      stop("LDprofile_Beta_b must contain the same number of values as there are bins given in LDprofile_bins")
    }
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

  # Force the R code to print decimals in full rather than in scientific format
  oldOptions<-options(scipen=999)
  on.exit(options(oldOptions))

  # Set up output list
  outputLength<-length(pos[pos>=X[1] & pos <= X[2]])
  outputList<-list(position=pos[pos>=X[1] & pos <= X[2]],LR=rep(NA,outputLength),L_plus_R=rep(NA,outputLength))
  if (is.null(dist)==FALSE & is.null(LDprofile_bins)==FALSE & is.null(LDprofile_rsq)==FALSE){
    outputList$Zalpha_expected<-rep(NA,outputLength)
    outputList$Zbeta_expected<-rep(NA,outputLength)
  }
  if (is.null(x)==FALSE){
    outputList$Zalpha<-rep(NA,outputLength)
    outputList$Zbeta<-rep(NA,outputLength)
    if (is.null(dist)==FALSE & is.null(LDprofile_bins)==FALSE & is.null(LDprofile_rsq)==FALSE){
      outputList$Zalpha_rsq_over_expected<-rep(NA,outputLength)
      outputList$Zalpha_log_rsq_over_expected<-rep(NA,outputLength)
      outputList$Zbeta_rsq_over_expected<-rep(NA,outputLength)
      outputList$Zbeta_log_rsq_over_expected<-rep(NA,outputLength)
      if (is.null(LDprofile_sd)==FALSE){
        outputList$Zalpha_Zscore<-rep(NA,outputLength)
        outputList$Zbeta_Zscore<-rep(NA,outputLength)
      }
      if (is.null(LDprofile_Beta_a)==FALSE & is.null(LDprofile_Beta_b)==FALSE){
        outputList$Zalpha_BetaCDF<-rep(NA,outputLength)
        outputList$Zbeta_BetaCDF<-rep(NA,outputLength)
      }
    }
  }

  # Loop over each position in the output list and calculate Zalpha
  for (i in 1:outputLength){

    # Current physical position in chromosome
    currentPos<-outputList$position[i]

    ## check L, R and LR
    noL <- length(pos[pos>=currentPos-ws/2 & pos < currentPos]) ## Number of SNPs to the left of the current SNP
    noR <- length(pos[pos<=currentPos+ws/2 & pos > currentPos]) ## Number of SNPs to the right of the current SNP
    outputList$LR[i]<-noL*noR
    outputList$L_plus_R[i]<-choose(noL,2)+choose(noR,2)
    if  (noL < minRandL || noR < minRandL || noL*noR < minRL){
      #NA for everything - leave as is
    } else {
      if (is.null(x)==FALSE){
        ##Left
        Lrsq <- lower_triangle(cor(t(x[pos>=currentPos-ws/2 & pos < currentPos,]),use="pairwise.complete.obs")^2)
        ##Right
        Rrsq<-lower_triangle(cor(t(x[pos<=currentPos+ws/2 & pos > currentPos,]),use="pairwise.complete.obs")^2)
        ##Over
        rsq<-as.vector(t((cor(t(x[pos>=currentPos-ws/2 & pos<=currentPos+ws/2,]),use="pairwise.complete.obs")^2)[1:noL,(noL+2):(noL+noR+1)]))
        outputList$Zalpha[i]<-(sum(Lrsq)/choose(noL,2)+sum(Rrsq)/choose(noR,2))/2
        outputList$Zbeta[i]<-sum(rsq)/(noL*noR)
      } else {
        Lrsq<-NA
        Rrsq<-NA
        rsq<-NA
      }
      if (is.null(dist)==FALSE & is.null(LDprofile_bins)==FALSE & is.null(LDprofile_rsq)==FALSE){
        #Left
        bins<-sapply(lower_triangle(outer(dist[pos>=currentPos-ws/2 & pos < currentPos],dist[pos>=currentPos-ws/2 & pos < currentPos],"-")),assign_bins,bin_size=bin_size)
        bins[bins>max(LDprofile_bins)]<-max(LDprofile_bins)
        LrsqExp<-merge(data.frame(bins=as.character(bins),Lrsq),data.frame(LDprofile_bins=as.character(LDprofile_bins),LDprofile_rsq),by.x="bins",by.y="LDprofile_bins",all.x=TRUE,sort=FALSE)
        #Right
        bins<-sapply(lower_triangle(outer(dist[pos<=currentPos+ws/2 & pos > currentPos],dist[pos<=currentPos+ws/2 & pos > currentPos],"-")),assign_bins,bin_size=bin_size)
        bins[bins>max(LDprofile_bins)]<-max(LDprofile_bins)
        RrsqExp<-merge(data.frame(bins=as.character(bins),Rrsq),data.frame(LDprofile_bins=as.character(LDprofile_bins),LDprofile_rsq),by.x="bins",by.y="LDprofile_bins",all.x=TRUE,sort=FALSE)
        #Over
        bins<-sapply(outer(dist[pos<=currentPos+ws/2 & pos > currentPos],dist[pos>=currentPos-ws/2 & pos < currentPos],"-"),assign_bins,bin_size=bin_size)
        bins[bins>max(LDprofile_bins)]<-max(LDprofile_bins)
        rsqExp<-merge(data.frame(bins=as.character(bins),rsq),data.frame(LDprofile_bins=as.character(LDprofile_bins),LDprofile_rsq),by.x="bins",by.y="LDprofile_bins",all.x=TRUE,sort=FALSE)

        outputList$Zalpha_expected[i]<-(sum(LrsqExp$LDprofile_rsq)/choose(noL,2)+sum(RrsqExp$LDprofile_rsq)/choose(noR,2))/2
        outputList$Zbeta_expected[i]<-sum(rsqExp$LDprofile_rsq)/(noL*noR)
        if (is.null(x)==FALSE){
          outputList$Zalpha_rsq_over_expected[i]<-(sum(LrsqExp$Lrsq/LrsqExp$LDprofile_rsq)/choose(noL,2)+sum(RrsqExp$Rrsq/RrsqExp$LDprofile_rsq)/choose(noR,2))/2
          outputList$Zbeta_rsq_over_expected[i]<-sum(rsqExp$rsq/rsqExp$LDprofile_rsq)/(noL*noR)
          #removes zeros by replacing with lowest correlation greater than zero, for logging
          LrsqExplog <- LrsqExp
          RrsqExplog <- RrsqExp
          rsqExplog <- rsqExp
          LrsqExplog$Lrsq[LrsqExplog$Lrsq==0]<-min(LrsqExplog$Lrsq[LrsqExplog$Lrsq>0])
          RrsqExplog$Rrsq[RrsqExplog$Rrsq==0]<-min(RrsqExplog$Rrsq[RrsqExplog$Rrsq>0])
          rsqExplog$rsq[rsqExplog$rsq==0]<-min(rsqExplog$rsq[rsqExplog$rsq>0])
          outputList$Zalpha_log_rsq_over_expected[i]<-(sum(log10(LrsqExplog$Lrsq/LrsqExp$LDprofile_rsq))/choose(noL,2)+sum(log10(RrsqExplog$Rrsq/RrsqExp$LDprofile_rsq))/choose(noR,2))/2
          outputList$Zbeta_log_rsq_over_expected[i]<-sum(log10(rsqExplog$rsq/rsqExp$LDprofile_rsq))/(noL*noR)

          if (is.null(LDprofile_sd)==FALSE){
            LrsqExp<-merge(LrsqExp,data.frame(LDprofile_bins=as.character(LDprofile_bins),LDprofile_sd),by.x="bins",by.y="LDprofile_bins",all.x=TRUE,sort=FALSE)
            RrsqExp<-merge(RrsqExp,data.frame(LDprofile_bins=as.character(LDprofile_bins),LDprofile_sd),by.x="bins",by.y="LDprofile_bins",all.x=TRUE,sort=FALSE)
            rsqExp<-merge(rsqExp,data.frame(LDprofile_bins=as.character(LDprofile_bins),LDprofile_sd),by.x="bins",by.y="LDprofile_bins",all.x=TRUE,sort=FALSE)
            outputList$Zalpha_Zscore[i]<-(sum((LrsqExp$Lrsq-LrsqExp$LDprofile_rsq)/LrsqExp$LDprofile_sd)/choose(noL,2)+sum((RrsqExp$Rrsq-RrsqExp$LDprofile_rsq)/RrsqExp$LDprofile_sd)/choose(noR,2))/2
            outputList$Zbeta_Zscore[i]<-sum((rsqExp$rsq-rsqExp$LDprofile_rsq)/rsqExp$LDprofile_sd)/(noL*noR)
          }
          if (is.null(LDprofile_Beta_a)==FALSE & is.null(LDprofile_Beta_b)==FALSE){
            LrsqExp<-merge(LrsqExp,data.frame(LDprofile_bins=as.character(LDprofile_bins),LDprofile_Beta_a,LDprofile_Beta_b),by.x="bins",by.y="LDprofile_bins",all.x=TRUE,sort=FALSE)
            RrsqExp<-merge(RrsqExp,data.frame(LDprofile_bins=as.character(LDprofile_bins),LDprofile_Beta_a,LDprofile_Beta_b),by.x="bins",by.y="LDprofile_bins",all.x=TRUE,sort=FALSE)
            rsqExp<-merge(rsqExp,data.frame(LDprofile_bins=as.character(LDprofile_bins),LDprofile_Beta_a,LDprofile_Beta_b),by.x="bins",by.y="LDprofile_bins",all.x=TRUE,sort=FALSE)
            outputList$Zalpha_BetaCDF[i]<-(sum(pbeta(LrsqExp$Lrsq,LrsqExp$LDprofile_Beta_a,LrsqExp$LDprofile_Beta_b))/choose(noL,2)+sum(pbeta(RrsqExp$Rrsq,RrsqExp$LDprofile_Beta_a,RrsqExp$LDprofile_Beta_b))/choose(noR,2))/2
            outputList$Zbeta_BetaCDF[i]<-sum(pbeta(rsqExp$rsq,rsqExp$LDprofile_Beta_a,rsqExp$LDprofile_Beta_b))/(noL*noR)
          }
        }
      }
    }
  }
  if (length(outputList)>3){
    if (sum(sapply(outputList[-c(1:3)],function(x) sum(is.na(x))==outputLength))==length(outputList)-3){
      warning("No statistics were calculated, try reducing minRandL and minRL or increasing the window size")
    }
  }
  return(outputList)
}
