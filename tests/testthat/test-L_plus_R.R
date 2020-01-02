## Create basic example data frame to check LplusR with. 15 SNPs over 5 chromosomes.

df<-data.frame(
  SNP=c("SNP1","SNP2","SNP3","SNP4","SNP5","SNP6","SNP7","SNP8","SNP9","SNP10","SNP11","SNP12","SNP13","SNP14","SNP15"),
  POS=c(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500),
  C1=c(1,1,2,1,2,1,1,2,1,2,1,2,1,1,1),
  C2=c(2,2,1,2,1,2,1,2,2,2,1,2,1,1,2),
  C3=c(2,1,2,2,2,1,1,2,2,1,2,2,1,1,2),
  C4=c(1,1,2,1,2,2,1,1,1,1,1,2,2,2,2),
  C5=c(1,1,2,1,2,1,2,1,1,1,1,1,2,1,1)
)

## test that LplusR is calculated correctly

test_that("LplusR calcualtes LplusR statistic correctly", {

  expect_equal(LplusR(pos = df$POS, x = as.matrix(df[,3:7]), ws  = 3000, X = NULL),
               data.frame(
                 POS=c(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500),
                 LplusR=c(NA,NA,67,58,51,46,43,42,43,46,51,58,67,NA,NA)
               ))
})

## Test the function with a different window size

test_that("LplusR calcualtes LplusR statistic correctly with a different window size", {

  expect_equal(LplusR(pos = df$POS, x = as.matrix(df[,3:7]), ws  = 1100, X = NULL),
               data.frame(
                 POS=c(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500),
                 LplusR=c(NA,NA,11,13,16,20,20,20,20,20,16,13,11,NA,NA)
               ))
})

## Test the function with a character matrix as x

test_that("LplusR calcualtes LplusR statistic correctly with character matrix", {

  df1<-df
  df1[df1==1]<-"A"
  df1[df1==2]<-"B"
  expect_equal(LplusR(pos = df1$POS, x = as.matrix(df1[,3:7]), ws  = 3000, X = NULL),
               data.frame(
                 POS=c(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500),
                 LplusR=c(NA,NA,67,58,51,46,43,42,43,46,51,58,67,NA,NA)
               ))
})

## Test the function with X supplied as a parameter

test_that("LplusR calcualtes LplusR statistic correctly with X supplied", {

  expect_equal(LplusR(pos = df$POS, x = as.matrix(df[,3:7]), ws  = 3000, X = c(700,900)),
               data.frame(
                 POS=c(700,800,900),
                 LplusR=c(43,42,43)
               ))
})

## Test the function with X supplied as a parameter outside of the region defined in pos

test_that("LplusR fails with an X supplied outside of the region defined in pos", {

  expect_error(LplusR(pos = df$POS, x = as.matrix(df[,3:7]), ws  = 3000, X = c(7000,9000)),
               "The region specified by X is outside the region contained in the pos vector")
})

## Test the function with X supplied as character

test_that("LplusR fails with an X supplied as a character", {

  expect_error(LplusR(pos = df$POS, x = as.matrix(df[,3:7]), ws  = 3000, X = c("700bp","900bp")),
               "X should be a numeric vector of length 2 e.g. c(100,200)",
               fixed=TRUE)
})

## Test the function with X supplied as only one number

test_that("LplusR fails with an X supplied as only one number", {

  expect_error(LplusR(pos = df$POS, x = as.matrix(df[,3:7]), ws  = 3000, X = 700),
               "X should be a numeric vector of length 2 e.g. c(100,200)",
               fixed=TRUE)
})

## Test the function with X supplied with too many numbers

test_that("LplusR fails with an X supplied with too many numbers", {

  expect_error(LplusR(pos = df$POS, x = as.matrix(df[,3:7]), ws  = 3000, X = c(700,900,1100)),
               "X should be a numeric vector of length 2 e.g. c(100,200)",
               fixed=TRUE)
})

## Test the function with ws non-numeric

test_that("LplusR fails when ws is non-numeric", {

  expect_error(LplusR(pos = df$POS, x = as.matrix(df[,3:7]), ws  = "3000bp", X = NULL),
               "ws must be a number greater than 0")
})

## Test the function with ws = 0

test_that("LplusR fails when ws is zero", {

  expect_error(LplusR(pos = df$POS, x = as.matrix(df[,3:7]), ws  = 0, X = NULL),
               "ws must be a number greater than 0")
})

## Test the function with positions non-numeric

test_that("LplusR fails when pos is non-numeric", {

  expect_error(LplusR(pos = paste0(df$POS,"bp"), x = as.matrix(df[,3:7]), ws  = 3000, X = NULL),
               "pos must be a numeric vector")
})

## Test the function with x not a matrix

test_that("LplusR fails when x is not a matrix", {

  expect_error(LplusR(pos = df$POS, x = df[,3:7], ws  = 3000, X = NULL),
               "x must be a matrix")
})

## Test the function with x not having the correct amount of rows

test_that("LplusR fails when the number of rows in x is not equal to the length of pos", {

  expect_error(LplusR(pos = df$POS, x = t(as.matrix(df[,3:7])), ws  = 3000, X = NULL),
               "The number of rows in x must equal the number of SNP locations given in pos")
})

## Test the function with a SNP having only one allele

test_that("LplusR fails when a SNP has only one allele", {

  df1<-df
  df1[1,3:7]<-1
  expect_error(LplusR(pos = df1$POS, x = as.matrix(df1[,3:7]), ws  = 3000, X = NULL),
               "SNPs must all be biallelic")
})

## Test the function with a SNP having more than two alleles

test_that("LplusR fails when a SNP has more than two alleles", {

  df1<-df
  df1[1,7]<-3
  expect_error(LplusR(pos = df1$POS, x = as.matrix(df1[,3:7]), ws  = 3000, X = NULL),
               "SNPs must all be biallelic")
})

