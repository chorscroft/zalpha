
#' Runs the Zalpha function
#'
#' Returns a Zalpha value for each SNP location supplied to the function.
#'
#' @param pos A numeric vector of SNP locations
#' @param x A matrix of SNP values. Columns represent chromosomes, rows are SNP locations. Hence, the number of rows should equal the length of the `pos` vector. SNPs should all be biallelic.
#' @param ws The window size which the `Zalpha` statistic will be calculated over. This should be on the same scale as the `pos` vector.
#' @param minRandL Minimum number of SNPs in each set R and L for the statistic to be calculated. Default is 4.
#' @param minRL Minimum value for the product of the set sizes for R and L. Default is 25.
#' @param X Optional. Specify a region of the chromosome to calculate Zalpha for in the format `c(startposition, endposition)`. The start position and the end position should be within the extremes of the `pos` vector. If not supplied, the function will calculate Zalpha for every SNP in the `pos` vector.
#'
#' @return A dataframe containing the positions and the Zalpha values
#' @export

Zalpha <- function(pos, x, ws, minRandL = 4, minRL = 25, X = NULL) {
  #Check things are in the correct format
  #Change matrix x to numeric if it isn't already


  print("Zalpha function")
}
