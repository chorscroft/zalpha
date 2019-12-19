
#' Runs the Zalpha function
#'
#' Returns a Zalpha value for each SNP location supplied to the function.
#' For more information about the Zalpha statistic please see Jacobs (2016).
#'
#' @importFrom stats cor
#'
#' @param pos A numeric vector of SNP locations
#' @param x A matrix of SNP values. Columns represent chromosomes, rows are SNP locations. Hence, the number of rows should equal the length of the `pos` vector. SNPs should all be biallelic.
#' @param ws The window size which the `Zalpha` statistic will be calculated over. This should be on the same scale as the `pos` vector.
#' @param minRandL Minimum number of SNPs in each set R and L for the statistic to be calculated. Default is 4.
#' @param minRL Minimum value for the product of the set sizes for R and L. Default is 25.
#' @param X Optional. Specify a region of the chromosome to calculate Zalpha for in the format `c(startposition, endposition)`. The start position and the end position should be within the extremes of the positions given in the `pos` vector. If not supplied, the function will calculate Zalpha for every SNP in the `pos` vector.
#'
#' @return A data frame containing the SNP positions and the Zalpha values for those SNPs
#' @export

Zalpha <- function(pos, x, ws, minRandL = 4, minRL = 25, X = NULL) {
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
    x1<-matrix(as.numeric(factor(x)),nrow=dim(x)[1])
  }

  # Set up output data frame
  outputDF<-data.frame(POS=pos[pos>=X[1] & pos <= X[2]],Zalpha=NA)


  # Loop over each position in the output data frame and calculate Zalpha
  for (i in 1:nrow(outputDF)){

    # Current physical position in chromosome
    currentPos<-outputDF$POS[i]

    ## check L, R and LR
    noL <- length(pos[pos>=currentPos-ws/2 & pos < currentPos]) ## Number of SNPs to the left of the current SNP
    noR <- length(pos[pos<=currentPos+ws/2 & pos > currentPos]) ## Number of SNPs to the right of the current SNP
    if  (noL < minRandL || noR < minRandL || noL*noR < minRL){
      #NA
      outputDF$Zalpha[i]<-NA
    } else {
      ##Left
      LrsqSum<-(sum((cor(t(x[pos>=currentPos-ws/2 & pos < currentPos,]))^2))-noL)/2
      ##Right
      RrsqSum<-(sum((cor(t(x[pos<=currentPos+ws/2 & pos > currentPos,]))^2))-noR)/2
      outputDF$Zalpha[i]<-(LrsqSum/choose(noL,2)+RrsqSum/choose(noR,2))/2
    }
  }
  return(outputDF)
}
