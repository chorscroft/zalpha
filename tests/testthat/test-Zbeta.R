
## Create basic example data frame to check Zbeta with. 15 SNPs over 5 chromosomes.

df<-data.frame(
  SNP=c("SNP1","SNP2","SNP3","SNP4","SNP5","SNP6","SNP7","SNP8","SNP9","SNP10","SNP11","SNP12","SNP13","SNP14","SNP15"),
  POS=c(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500),
  C1=c(1,1,2,1,2,1,1,2,1,2,1,2,1,1,1),
  C2=c(2,2,1,2,1,2,1,2,2,2,1,2,1,1,2),
  C3=c(2,1,2,2,2,1,1,2,2,1,2,2,1,1,2),
  C4=c(1,1,2,1,2,2,1,1,1,1,1,2,2,2,2),
  C5=c(1,1,2,1,2,1,2,1,1,1,1,1,2,1,1)
)

## test that Zbeta is calculated correctly

test_that("Zbeta calculates Zbeta statistic correctly", {

  expect_equal(Zbeta(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), minRandL = 4, minRL = 25, X = NULL),
               list(
                 position=c(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500),
                 Zbeta=c(NA,NA,NA,NA,
                          (10+5/18)/40,
                          (10+35/36)/45,
                          (11+103/144)/48,
                          (12+73/144)/49,
                          (11+83/144)/48,
                          (11+17/48)/45,
                          (10+65/144)/40,
                          NA,NA,NA,NA)
               ))
})

## Test the function with a different window size

test_that("Zbeta calculates Zbeta statistic correctly with a different window size", {

  expect_equal(Zbeta(pos = df$POS, ws  = 1100, x = as.matrix(df[,3:7]), minRandL = 4, minRL = 25, X = NULL),
               list(
                 position=c(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500),
                 Zbeta=c(NA,NA,NA,NA,NA,
                          (7+11/72)/25,
                          (5+5/9)/25,
                          (6+41/144)/25,
                          (6+101/144)/25,
                          (7+17/144)/25,
                          NA,NA,NA,NA,NA)
               ))
})

## Test the function with a character matrix as x

test_that("Zbeta calculates Zbeta statistic correctly with character matrix", {

  df1<-df
  df1[df1==1]<-"A"
  df1[df1==2]<-"B"
  expect_equal(Zbeta(pos = df1$POS, ws  = 3000, x = as.matrix(df1[,3:7]), minRandL = 4, minRL = 25, X = NULL),
               list(
                 position=c(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500),
                 Zbeta=c(NA,NA,NA,NA,
                         (10+5/18)/40,
                         (10+35/36)/45,
                         (11+103/144)/48,
                         (12+73/144)/49,
                         (11+83/144)/48,
                         (11+17/48)/45,
                         (10+65/144)/40,
                         NA,NA,NA,NA)
               ))
})

## Test the function with X supplied as a parameter

test_that("Zbeta calculates Zbeta statistic correctly with X supplied", {

  expect_equal(Zbeta(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), minRandL = 4, minRL = 25, X = c(700,900)),
               list(
                 position=c(700,800,900),
                 Zbeta=c((11+103/144)/48,
                         (12+73/144)/49,
                         (11+83/144)/48)
               ))
})

## Test the function with X supplied as a parameter outside of the region defined in pos

test_that("Zbeta fails with an X supplied outside of the region defined in pos", {

  expect_error(Zbeta(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), minRandL = 4, minRL = 25, X = c(7000,9000)),
               "The region specified by X is outside the region contained in the pos vector")
})

## Test the function with X supplied as character

test_that("Zbeta fails with an X supplied as a character", {

  expect_error(Zbeta(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), minRandL = 4, minRL = 25, X = c("700bp","900bp")),
               "X should be a numeric vector of length 2 e.g. c(100,200)",
               fixed=TRUE)
})

## Test the function with X supplied as only one number

test_that("Zbeta fails with an X supplied as only one number", {

  expect_error(Zbeta(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), minRandL = 4, minRL = 25, X = 700),
               "X should be a numeric vector of length 2 e.g. c(100,200)",
               fixed=TRUE)
})

## Test the function with X supplied with too many numbers

test_that("Zbeta fails with an X supplied with too many numbers", {

  expect_error(Zbeta(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), minRandL = 4, minRL = 25, X = c(700,900,1100)),
               "X should be a numeric vector of length 2 e.g. c(100,200)",
               fixed=TRUE)
})

## Test the function with ws non-numeric

test_that("Zbeta fails when ws is non-numeric", {

  expect_error(Zbeta(pos = df$POS, ws  = "3000bp", x = as.matrix(df[,3:7]), minRandL = 4, minRL = 25, X = NULL),
               "ws must be a number greater than 0")
})

## Test the function with ws = 0

test_that("Zbeta fails when ws is zero", {

  expect_error(Zbeta(pos = df$POS, ws  = 0, x = as.matrix(df[,3:7]), minRandL = 4, minRL = 25, X = NULL),
               "ws must be a number greater than 0")
})

## Test the function with positions non-numeric

test_that("Zbeta fails when pos is non-numeric", {

  expect_error(Zbeta(pos = paste0(df$POS,"bp"), ws  = 3000, x = as.matrix(df[,3:7]), minRandL = 4, minRL = 25, X = NULL),
               "pos must be a numeric vector")
})

## Test the function with minLandR as non-numeric

test_that("Zbeta fails when minLandR is non-numeric", {

  expect_error(Zbeta(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), minRandL = "4snps", minRL = 25, X = NULL),
               "minRandL must be a number greater than or equal to 0")
})

## Test the function with minLandR as negative

test_that("Zbeta fails when minLandR is negative", {

  expect_error(Zbeta(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), minRandL = -1, minRL = 25, X = NULL),
               "minRandL must be a number greater than or equal to 0")
})

## Test the function with minLR as non-numeric

test_that("Zbeta fails when minLR is non-numeric", {

  expect_error(Zbeta(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), minRandL = 4, minRL = "25b", X = NULL),
               "minRL must be a number greater than or equal to 0")
})

## Test the function with minLR as negative

test_that("Zbeta fails when minLR is negative", {

  expect_error(Zbeta(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), minRandL = 4, minRL = -25, X = NULL),
               "minRL must be a number greater than or equal to 0")
})

## Test the function with x not a matrix

test_that("Zbeta fails when x is not a matrix", {

  expect_error(Zbeta(pos = df$POS, ws  = 3000, x = df[,3:7], minRandL = 4, minRL = 25, X = NULL),
               "x must be a matrix")
})

## Test the function with x not having the correct amount of rows

test_that("Zbeta fails when the number of rows in x is not equal to the length of pos", {

  expect_error(Zbeta(pos = df$POS, ws  = 3000, x = t(as.matrix(df[,3:7])), minRandL = 4, minRL = 25, X = NULL),
               "The number of rows in x must equal the number of SNP locations given in pos")
})

## Test the function with a SNP having only one allele

test_that("Zbeta fails when a SNP has only one allele", {

  df1<-df
  df1[1,3:7]<-1
  expect_error(Zbeta(pos = df1$POS, ws  = 3000, x = as.matrix(df1[,3:7]), minRandL = 4, minRL = 25, X = NULL),
               "SNPs must all be biallelic")
})

## Test the function with a SNP having more than two alleles

test_that("Zbeta fails when a SNP has more than two alleles", {

  df1<-df
  df1[1,7]<-3
  expect_error(Zbeta(pos = df1$POS, ws  = 3000, x = as.matrix(df1[,3:7]), minRandL = 4, minRL = 25, X = NULL),
               "SNPs must all be biallelic")
})

## Test the warning that all values returned are NA by increasing minRL

test_that("Zbeta warns about all NAs", {

  expect_warning(Zbeta(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), minRandL = 4, minRL = 50, X = NULL),
                 "No Zbeta values were calculated, try reducing minRandL and minRL or increasing the window size")
})

## test that Zbeta works with a missing value
df1<-df
df1$C1[15]<-NA
test_that("Zbeta calculates Zbeta statistic correctly with missing value", {

  expect_equal(Zbeta(pos = df1$POS, ws  = 3000, x = as.matrix(df1[,3:7]), minRandL = 4, minRL = 25, X = NULL),
               list(
                 position=c(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500),
                 Zbeta=c(NA,NA,NA,NA,
                         0.248611111111111,
                         0.235185185185185,
                         0.233651620370370,
                         0.257794784580499,
                         0.250144675925926,
                         0.259413580246914,
                         0.271354166666667,
                         NA,NA,NA,NA)
               ))
})
