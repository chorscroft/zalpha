
#' Runs the L_plus_R function
#'
#' Returns the \eqn{{|L| \choose 2} + {|R| \choose 2}}{(|L| choose 2) + (|R| choose 2)} value for each SNP location supplied to the function.
#' For more information about the \code{L_plus_R} diversity statistic please see Jacobs (2016).
#'
#'
#'
#' @param pos A numeric vector of SNP locations
#' @param ws The window size which the \code{L_plus_R} statistic will be calculated over. This should be on the same scale as the \code{pos} vector.
#' @param X Optional. Specify a region of the chromosome to calculate \code{L_plus_R} for in the format \code{c(startposition, endposition)}. The start position and the end position should be within the extremes of the positions given in the \code{pos} vector. If not supplied, the function will calculate L_plus_R for every SNP in the \code{pos} vector.
#'
#' @return A list containing the SNP positions and the \code{L_plus_R }values for those SNPs
#' @references Jacobs, G.S., T.J. Sluckin, and T. Kivisild, \emph{Refining the Use of Linkage Disequilibrium as a Robust Signature of Selective Sweeps.} Genetics, 2016. \strong{203}(4): p. 1807
#' @examples
#' ## load the snps example dataset
#' data(snps)
#' ## run L_plus_R over all the SNPs with a window size of 3000 bp
#' L_plus_R(snps$positions,3000)
#' ## only return results for SNPs between locations 600 and 1500 bp
#' L_plus_R(snps$positions,3000,X=c(600,1500))
#'
#' @export

L_plus_R <- function(pos, ws, X = NULL) {
  #Check things are in the correct format

  #Check pos is a numeric vector
  if (is.numeric(pos) ==FALSE || is.vector(pos)==FALSE){
    stop("pos must be a numeric vector")
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

  # Set up output list
  outputLength<-length(pos[pos>=X[1] & pos <= X[2]])
  outputList<-list(position=pos[pos>=X[1] & pos <= X[2]],L_plus_R=rep(NA,outputLength))

  # Loop over each position in the output list and calculate L_plus_R
  for (i in 1:outputLength){

    # Current physical position in chromosome
    currentPos<-outputList$position[i]

    ## get L, R and L_plus_R
    noL <- length(pos[pos>=currentPos-ws/2 & pos < currentPos]) ## Number of SNPs to the left of the current SNP
    noR <- length(pos[pos<=currentPos+ws/2 & pos > currentPos]) ## Number of SNPs to the right of the current SNP
    if(noL < 2 || noR < 2){  #Must be at least 2 to calculate n choose 2
      outputList$L_plus_R[i]<-NA
    } else {
      outputList$L_plus_R[i]<-choose(noL,2)+choose(noR,2)
    }
  }
  return(outputList)
}
