
#' Runs the LR function
#'
#' Returns the \code{|L||R|} value for each SNP location supplied to the function.
#' For more information about the \code{|L||R|} diversity statistic please see Jacobs (2016).
#'
#' @param pos A numeric vector of SNP locations
#' @param x A matrix of SNP values. Columns represent chromosomes, rows are SNP locations. Hence, the number of rows should equal the length of the \code{pos} vector. SNPs should all be biallelic.
#' @param ws The window size which the \code{LR} statistic will be calculated over. This should be on the same scale as the \code{pos} vector.
#' @param X Optional. Specify a region of the chromosome to calculate LR for in the format \code{c(startposition, endposition)}. The start position and the end position should be within the extremes of the positions given in the \code{pos} vector. If not supplied, the function will calculate LR for every SNP in the \code{pos} vector.
#'
#' @return A data frame containing the SNP positions and the LR values for those SNPs
#' @references Jacobs, G.S., T.J. Sluckin, and T. Kivisild, \emph{Refining the Use of Linkage Disequilibrium as a Robust Signature of Selective Sweeps.} Genetics, 2016. \strong{203}(4): p. 1807
#' @export

LR <- function(pos, x, ws, X = NULL) {
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
  #Check windowsize is a number greater than 0
  if(is.numeric(ws) ==FALSE || ws <= 0){
    stop("ws must be a number greater than 0")
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

  # Set up output data frame
  outputDF<-data.frame(POS=pos[pos>=X[1] & pos <= X[2]],LR=NA)


  # Loop over each position in the output data frame and calculate LR
  for (i in 1:nrow(outputDF)){

    # Current physical position in chromosome
    currentPos<-outputDF$POS[i]

    ## get L, R and LR
    noL <- length(pos[pos>=currentPos-ws/2 & pos < currentPos]) ## Number of SNPs to the left of the current SNP
    noR <- length(pos[pos<=currentPos+ws/2 & pos > currentPos]) ## Number of SNPs to the right of the current SNP
    outputDF$LR[i]<-noL*noR
  }
  return(outputDF)
}
