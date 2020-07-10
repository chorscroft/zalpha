
## Compares a vector to a value and returns a vector of logical values
##
## @param vector The vector to be compared to the value
## @param value The value each item in the vector should be compared to
##
## @return A vector of TRUE and FALSE values
##
equal_vector <- function(vector, value){
  tempFunction<-function(x,value){
    isTRUE(all.equal(x,value))
  }
  return(sapply(vector,tempFunction,value))
}



