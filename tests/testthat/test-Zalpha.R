
## Create basic example data frame to check Zalpha with. 15 SNPs over 5 chromosomes.

df<-data.frame(
  SNP=c("SNP1","SNP2","SNP3","SNP4","SNP5","SNP6","SNP7","SNP8","SNP9","SNP10","SNP11","SNP12","SNP13","SNP14","SNP15"),
  POS=c(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500),
  C1=c(1,1,2,1,2,1,1,2,1,2,1,2,1,1,1),
  C2=c(2,2,1,2,1,2,1,2,2,2,1,2,1,1,2),
  C3=c(2,1,2,2,2,1,1,2,2,1,2,2,1,1,2),
  C4=c(1,1,2,1,2,2,1,1,1,1,1,2,2,2,2),
  C5=c(1,1,2,1,2,1,2,1,1,1,1,1,2,1,1)
)

## test that Zalpha is calculated correctly

test_that("Zalpha calculates Zalpha statistic correctly", {

  expect_equal(Zalpha(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), minRandL = 4, minRL = 25, X = NULL),
               list(
                 position=c(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500),
                 Zalpha=c(NA,NA,NA,NA,
                          ((3+1/2)/6+(11+41/144)/45)/2,
                          ((6+1/4)/10+(9+41/48)/36)/2,
                          ((7+31/72)/15+(7+13/48)/28)/2,
                          ((8+17/144)/21+(4+7/16)/21)/2,
                          ((9+131/144)/28+(2+13/16)/15)/2,
                          ((13+97/144)/36+(1+121/144)/10)/2,
                          ((15+25/48)/45+(1+55/144)/6)/2,
                          NA,NA,NA,NA)
               ),tolerance=0.0001)
})

## Test the function with a different window size

test_that("Zalpha calculates Zalpha statistic correctly with a different window size", {

  expect_equal(Zalpha(pos = df$POS, ws  = 1100, x = as.matrix(df[,3:7]), minRandL = 4, minRL = 25, X = NULL),
               list(
                 position=c(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500),
                 Zalpha=c(NA,NA,NA,NA,NA,
                          ((6+1/4)/10+(2+19/48)/10)/2,
                          ((5+5/18)/10+(2+19/48)/10)/2,
                          ((2+71/72)/10+(2+19/48)/10)/2,
                          ((2+3/16)/10+(2+7/144)/10)/2,
                          ((2+3/16)/10+(1+121/144)/10)/2,
                          NA,NA,NA,NA,NA)
               ),tolerance=0.0001)
})

## Test the function with a character matrix as x

test_that("Zalpha calculates Zalpha statistic correctly with character matrix", {

  df1<-df
  df1[df1==1]<-"A"
  df1[df1==2]<-"B"
  expect_equal(Zalpha(pos = df1$POS, ws  = 3000, x = as.matrix(df1[,3:7]), minRandL = 4, minRL = 25, X = NULL),
               list(
                 position=c(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500),
                 Zalpha=c(NA,NA,NA,NA,
                          ((3+1/2)/6+(11+41/144)/45)/2,
                          ((6+1/4)/10+(9+41/48)/36)/2,
                          ((7+31/72)/15+(7+13/48)/28)/2,
                          ((8+17/144)/21+(4+7/16)/21)/2,
                          ((9+131/144)/28+(2+13/16)/15)/2,
                          ((13+97/144)/36+(1+121/144)/10)/2,
                          ((15+25/48)/45+(1+55/144)/6)/2,
                          NA,NA,NA,NA)
               ),tolerance=0.0001)
})

## Test the function with X supplied as a parameter

test_that("Zalpha calculates Zalpha statistic correctly with X supplied", {

  expect_equal(Zalpha(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), minRandL = 4, minRL = 25, X = c(700,900)),
               list(
                 position=c(700,800,900),
                 Zalpha=c(((7+31/72)/15+(7+13/48)/28)/2,
                          ((8+17/144)/21+(4+7/16)/21)/2,
                          ((9+131/144)/28+(2+13/16)/15)/2)
               ),tolerance=0.0001)
})

## Test the function with X supplied as a parameter outside of the region defined in pos

test_that("Zalpha fails with an X supplied outside of the region defined in pos", {

  expect_error(Zalpha(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), minRandL = 4, minRL = 25, X = c(7000,9000)),
               "The region specified by X is outside the region contained in the pos vector")
})

## Test the function with X supplied as character

test_that("Zalpha fails with an X supplied as a character", {

  expect_error(Zalpha(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), minRandL = 4, minRL = 25, X = c("700bp","900bp")),
               "X should be a numeric vector of length 2 e.g. c(100,200)",
               fixed=TRUE)
})

## Test the function with X supplied as only one number

test_that("Zalpha fails with an X supplied as only one number", {

  expect_error(Zalpha(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), minRandL = 4, minRL = 25, X = 700),
               "X should be a numeric vector of length 2 e.g. c(100,200)",
               fixed=TRUE)
})

## Test the function with X supplied with too many numbers

test_that("Zalpha fails with an X supplied with too many numbers", {

  expect_error(Zalpha(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), minRandL = 4, minRL = 25, X = c(700,900,1100)),
               "X should be a numeric vector of length 2 e.g. c(100,200)",
               fixed=TRUE)
})

## Test the function with ws non-numeric

test_that("Zalpha fails when ws is non-numeric", {

  expect_error(Zalpha(pos = df$POS, ws  = "3000bp", x = as.matrix(df[,3:7]), minRandL = 4, minRL = 25, X = NULL),
               "ws must be a number greater than 0")
})

## Test the function with ws = 0

test_that("Zalpha fails when ws is zero", {

  expect_error(Zalpha(pos = df$POS, ws  = 0, x = as.matrix(df[,3:7]), minRandL = 4, minRL = 25, X = NULL),
               "ws must be a number greater than 0")
})

## Test the function with positions non-numeric

test_that("Zalpha fails when pos is non-numeric", {

  expect_error(Zalpha(pos = paste0(df$POS,"bp"), ws  = 3000, x = as.matrix(df[,3:7]), minRandL = 4, minRL = 25, X = NULL),
               "pos must be a numeric vector")
})

## Test the function with minLandR as non-numeric

test_that("Zalpha fails when minLandR is non-numeric", {

  expect_error(Zalpha(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), minRandL = "4snps", minRL = 25, X = NULL),
               "minRandL must be a number greater than or equal to 0")
})

## Test the function with minLandR as negative

test_that("Zalpha fails when minLandR is negative", {

  expect_error(Zalpha(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), minRandL = -1, minRL = 25, X = NULL),
               "minRandL must be a number greater than or equal to 0")
})

## Test the function with minLR as non-numeric

test_that("Zalpha fails when minLR is non-numeric", {

  expect_error(Zalpha(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), minRandL = 4, minRL = "25b", X = NULL),
               "minRL must be a number greater than or equal to 0")
})

## Test the function with minLR as negative

test_that("Zalpha fails when minLR is negative", {

  expect_error(Zalpha(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), minRandL = 4, minRL = -25, X = NULL),
               "minRL must be a number greater than or equal to 0")
})

## Test the function with x not a matrix

test_that("Zalpha fails when x is not a matrix", {

  expect_error(Zalpha(pos = df$POS, ws  = 3000, x = df[,3:7], minRandL = 4, minRL = 25, X = NULL),
               "x must be a matrix")
})

## Test the function with x not having the correct amount of rows

test_that("Zalpha fails when the number of rows in x is not equal to the length of pos", {

  expect_error(Zalpha(pos = df$POS, ws  = 3000, x = t(as.matrix(df[,3:7])), minRandL = 4, minRL = 25, X = NULL),
               "The number of rows in x must equal the number of SNP locations given in pos")
})

## Test the function with a SNP having only one allele

test_that("Zalpha fails when a SNP has only one allele", {

  df1<-df
  df1[1,3:7]<-1
  expect_error(Zalpha(pos = df1$POS, ws  = 3000, x = as.matrix(df1[,3:7]), minRandL = 4, minRL = 25, X = NULL),
               "SNPs must all be biallelic")
})

## Test the function with a SNP having more than two alleles

test_that("Zalpha fails when a SNP has more than two alleles", {

  df1<-df
  df1[1,7]<-3
  expect_error(Zalpha(pos = df1$POS, ws  = 3000, x = as.matrix(df1[,3:7]), minRandL = 4, minRL = 25, X = NULL),
               "SNPs must all be biallelic")
})

## Test the warning that all values returned are NA by increasing minRL

test_that("Zalpha warns about all NAs", {

  expect_warning(Zalpha(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), minRandL = 4, minRL = 50, X = NULL),
               "No Zalpha values were calculated, try reducing minRandL and minRL or increasing the window size")
})

## test that zalpha works with a missing value

df1<-df
df1$C1[15]<-NA
test_that("Zalpha calculates Zalpha statistic correctly with missing value", {

  expect_equal(Zalpha(pos = df1$POS, ws  = 3000, x = as.matrix(df1[,3:7]), minRandL = 4, minRL = 25, X = NULL),
               list(
                 position=c(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500),
                 Zalpha=c(NA,NA,NA,NA,
                          0.434953703703704,
                          0.473283179012346,
                          0.397114748677249,
                          0.317791005291005,
                          0.300801917989418,
                          0.322897376543210,
                          0.360532407407407,
                          NA,NA,NA,NA)
               ),tolerance=0.0001)
})
