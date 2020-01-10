## Create basic example data frame to check LR with. 15 SNPs over 5 chromosomes.

df<-data.frame(
  SNP=c("SNP1","SNP2","SNP3","SNP4","SNP5","SNP6","SNP7","SNP8","SNP9","SNP10","SNP11","SNP12","SNP13","SNP14","SNP15"),
  POS=c(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500)
)

## test that LR is calculated correctly

test_that("LR calcualtes LR statistic correctly", {

  expect_equal(LR(pos = df$POS, ws  = 3000, X = NULL),
               list(
                 position=c(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500),
                 LR=c(0,13,24,33,40,45,48,49,48,45,40,33,24,13,0)
               ))
})

## Test the function with a different window size

test_that("LR calcualtes LR statistic correctly with a different window size", {

  expect_equal(LR(pos = df$POS, ws  = 1100, X = NULL),
               list(
                 position=c(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500),
                 LR=c(0,5,10,15,20,25,25,25,25,25,20,15,10,5,0)
              ))
})

## Test the function with X supplied as a parameter

test_that("LR calcualtes LR statistic correctly with X supplied", {

  expect_equal(LR(pos = df$POS, ws  = 3000, X = c(700,900)),
               list(
                 position=c(700,800,900),
                 LR=c(48,49,48)
               ))
})

## Test the function with X supplied as a parameter outside of the region defined in pos

test_that("LR fails with an X supplied outside of the region defined in pos", {

  expect_error(LR(pos = df$POS, ws  = 3000, X = c(7000,9000)),
               "The region specified by X is outside the region contained in the pos vector")
})

## Test the function with X supplied as character

test_that("LR fails with an X supplied as a character", {

  expect_error(LR(pos = df$POS, ws  = 3000, X = c("700bp","900bp")),
               "X should be a numeric vector of length 2 e.g. c(100,200)",
               fixed=TRUE)
})

## Test the function with X supplied as only one number

test_that("LR fails with an X supplied as only one number", {

  expect_error(LR(pos = df$POS, ws  = 3000, X = 700),
               "X should be a numeric vector of length 2 e.g. c(100,200)",
               fixed=TRUE)
})

## Test the function with X supplied with too many numbers

test_that("LR fails with an X supplied with too many numbers", {

  expect_error(LR(pos = df$POS, ws  = 3000, X = c(700,900,1100)),
               "X should be a numeric vector of length 2 e.g. c(100,200)",
               fixed=TRUE)
})

## Test the function with ws non-numeric

test_that("LR fails when ws is non-numeric", {

  expect_error(LR(pos = df$POS, ws  = "3000bp", X = NULL),
               "ws must be a number greater than 0")
})

## Test the function with ws = 0

test_that("LR fails when ws is zero", {

  expect_error(LR(pos = df$POS, ws  = 0, X = NULL),
               "ws must be a number greater than 0")
})

## Test the function with positions non-numeric

test_that("LR fails when pos is non-numeric", {

  expect_error(LR(pos = paste0(df$POS,"bp"), ws  = 3000, X = NULL),
               "pos must be a numeric vector")
})
