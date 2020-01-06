

## Helper function to get the lower triangle of a matrix
##
## @param x a matrix
##
## @return a vector of values from the lower triangle of the matrix
##
##
lower_triangle<-function(x){
  x[lower.tri(x)]
}
