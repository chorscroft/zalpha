## Create basic example data frame to check L_plus_R with. 15 SNPs over 5 chromosomes.

df<-data.frame(
  SNP=c("SNP1","SNP2","SNP3","SNP4","SNP5","SNP6","SNP7","SNP8","SNP9","SNP10","SNP11","SNP12","SNP13","SNP14","SNP15"),
  POS=c(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500)
)

## test that L_plus_R is calculated correctly

test_that("L_plus_R calculates L_plus_R statistic correctly", {

  expect_equal(L_plus_R(pos = df$POS, ws  = 3000, X = NULL),
               list(
                 position=c(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500),
                 L_plus_R=c(91,78,67,58,51,46,43,42,43,46,51,58,67,78,91)
               ))
})

## Test the function with a different window size

test_that("L_plus_R calculates L_plus_R statistic correctly with a different window size", {

  expect_equal(L_plus_R(pos = df$POS, ws  = 1100, X = NULL),
               list(
                 position=c(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500),
                 L_plus_R=c(10,10,11,13,16,20,20,20,20,20,16,13,11,10,10)
               ))
})

## Test the function with X supplied as a parameter

test_that("L_plus_R calculates L_plus_R statistic correctly with X supplied", {

  expect_equal(L_plus_R(pos = df$POS, ws  = 3000, X = c(700,900)),
               list(
                 position=c(700,800,900),
                 L_plus_R=c(43,42,43)
               ))
})

## Test the function with X supplied as a parameter outside of the region defined in pos

test_that("L_plus_R fails with an X supplied outside of the region defined in pos", {

  expect_error(L_plus_R(pos = df$POS, ws  = 3000, X = c(7000,9000)),
               "The region specified by X is outside the region contained in the pos vector")
})

## Test the function with X supplied as character

test_that("L_plus_R fails with an X supplied as a character", {

  expect_error(L_plus_R(pos = df$POS, ws  = 3000, X = c("700bp","900bp")),
               "X should be a numeric vector of length 2 e.g. c(100,200)",
               fixed=TRUE)
})

## Test the function with X supplied as only one number

test_that("L_plus_R fails with an X supplied as only one number", {

  expect_error(L_plus_R(pos = df$POS, ws  = 3000, X = 700),
               "X should be a numeric vector of length 2 e.g. c(100,200)",
               fixed=TRUE)
})

## Test the function with X supplied with too many numbers

test_that("L_plus_R fails with an X supplied with too many numbers", {

  expect_error(L_plus_R(pos = df$POS, ws  = 3000, X = c(700,900,1100)),
               "X should be a numeric vector of length 2 e.g. c(100,200)",
               fixed=TRUE)
})

## Test the function with ws non-numeric

test_that("L_plus_R fails when ws is non-numeric", {

  expect_error(L_plus_R(pos = df$POS, ws  = "3000bp", X = NULL),
               "ws must be a number greater than 0")
})

## Test the function with ws = 0

test_that("L_plus_R fails when ws is zero", {

  expect_error(L_plus_R(pos = df$POS, ws  = 0, X = NULL),
               "ws must be a number greater than 0")
})

## Test the function with positions non-numeric

test_that("L_plus_R fails when pos is non-numeric", {

  expect_error(L_plus_R(pos = paste0(df$POS,"bp"), ws  = 3000, X = NULL),
               "pos must be a numeric vector")
})

