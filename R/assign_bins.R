
## Calculates the bin a distance falls into
##
## Due to floating point errors, the floor function cannot be relied upon to
## round down a distance between two SNPs to the nearest bin. Thus it is
## necessary to first test the ceiling to see if it is in fact equal, before
## then using the floor function.
##
## For example, the SNP positions of SNPs 1 and 2 are 0.00235 and 0.00345
## respectively. The difference is 0.0011, and with a bin size of 0.0001, one
## would expect the difference to be assigned to the "0.0011" bin. However,
## the floor function in R instead sets the bin to 0.001, as it stores the
## result of 0.00345-0.00235 as a number a tiny bit smaller than 0.0011 due
## to floating point errors.
##
## @param bin_size a number representing the size of the bins in the LD profile
## @param number the number to be assigned to a bin
##
## @return a number representing the bin the number has been assigned to
##
assign_bins<-function(bin_size,number){
  ceilingTemp<-ceiling(number/bin_size)
  if(isTRUE(all.equal(ceilingTemp,number/bin_size))){
    return(ceilingTemp*bin_size)
  } else {
    return(floor(number/bin_size)*bin_size)
  }
}
