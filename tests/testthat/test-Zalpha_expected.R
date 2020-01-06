
## Create basic example data frame to check Zalpha_expected with.

df<-data.frame(
  SNP=c("SNP1","SNP2","SNP3","SNP4","SNP5","SNP6","SNP7","SNP8","SNP9","SNP10","SNP11","SNP12","SNP13","SNP14","SNP15"),
  POS=c(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500),
  CM=c(0,0.00101,0.00123,0.00207,0.00218,0.00223,0.00235,0.00251,0.00272,0.00289,0.00304,0.00316,0.00335,0.00345,0.00374)
)
LDprofile<-data.frame(
  cM_bin=seq(0,0.0049,0.0001),
  rsq=c(0.495714059385946,0.411014233619574,0.395532914859378,0.44354526861954,0.435550945945946,
        0.419534577303153,0.383410708498866,0.402439897670834,0.395945237932081,0.380436909495495,
        0.384229510773621,0.379494011054621,0.368118044626627,0.358753523652643,0.362330915047976,
        0.372142680938693,0.353703415234045,0.341173431307316,0.345578512726934,0.358159779825909,
        0.334687181997538,0.337442516960342,0.343721563062338,0.336509253287721,0.325271170690837,
        0.325488235114597,0.32078065970396,0.317594707821212,0.314613974963111,0.309774332543378,
        0.307619999017408,0.307004105181405,0.300279349768979,0.305505356875903,0.303179706053309,
        0.312783673753436,0.308540965611869,0.292690196360757,0.299428992380521,0.297144304197462,
        0.286701832971995,0.297894654636997,0.283984127549023,0.283253709389766,0.281330372503553,
        0.289052362087009,0.272730959781483,0.277399161038311,0.285764741944136,0.271195118636169)
)

## test that Zalpha_expected is calculated correctly

test_that("Zalpha_expected calcualtes Zalpha_expected statistic correctly", {

  expect_equal(Zalpha_expected(pos = df$POS, cM = df$CM, ws  = 3000, LDprofile_cM_bins = LDprofile$cM_bin, LDprofile_rsq = LDprofile$rsq, minRandL = 4, minRL = 25, X = NULL),
               list(
                 position=c(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500),
                 Zalpha_expected=c(NA,NA,NA,NA,
                                  0.390457304338967,
                                  0.392014054942343,
                                  0.397339546324536,
                                  0.398874728980465,
                                  0.400715520018796,
                                  0.401718327864356,
                                  0.399526703832832,
                                  NA,NA,NA,NA)
               ))
})

## Test the function with a different window size

test_that("Zalpha_expected calcualtes Zalpha_expected statistic correctly with a different window size", {

  expect_equal(Zalpha_expected(pos = df$POS, cM = df$CM, ws  = 1100, LDprofile_cM_bins = LDprofile$cM_bin, LDprofile_rsq = LDprofile$rsq, minRandL = 4, minRL = 25, X = NULL),
               list(
                 position=c(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500),
                 Zalpha_expected=c(NA,NA,NA,NA,NA,
                                   0.397641067838740,
                                   0.407771198285253,
                                   0.412830128247348,
                                   0.419512588688985,
                                   0.421536156959405,
                                   NA,NA,NA,NA,NA)
               ))
})

## Test the function with X supplied as a parameter

test_that("Zalpha_expected calcualtes Zalpha_expected statistic correctly with X supplied", {

  expect_equal(Zalpha_expected(pos = df$POS, cM = df$CM, ws  = 3000, LDprofile_cM_bins = LDprofile$cM_bin, LDprofile_rsq = LDprofile$rsq, minRandL = 4, minRL = 25, X = c(700,900)),
               list(
                 position=c(700,800,900),
                 Zalpha_expected=c(0.397339546324536,
                                   0.398874728980465,
                                   0.400715520018796)
               ))
})

## Test the function with X supplied as a parameter outside of the region defined in pos

test_that("Zalpha_expected fails with an X supplied outside of the region defined in pos", {

  expect_error(Zalpha_expected(pos = df$POS, cM = df$CM, ws  = 3000, LDprofile_cM_bins = LDprofile$cM_bin, LDprofile_rsq = LDprofile$rsq, minRandL = 4, minRL = 25, X = c(7000,9000)),
               "The region specified by X is outside the region contained in the pos vector")
})

## Test the function with X supplied as character

test_that("Zalpha_expected fails with an X supplied as a character", {

  expect_error(Zalpha_expected(pos = df$POS, cM = df$CM, ws  = 3000, LDprofile_cM_bins = LDprofile$cM_bin, LDprofile_rsq = LDprofile$rsq, minRandL = 4, minRL = 25, X = c("700bp","900bp")),
               "X should be a numeric vector of length 2 e.g. c(100,200)",
               fixed=TRUE)
})

## Test the function with X supplied as only one number

test_that("Zalpha_expected fails with an X supplied as only one number", {

  expect_error(Zalpha_expected(pos = df$POS, cM = df$CM, ws  = 3000, LDprofile_cM_bins = LDprofile$cM_bin, LDprofile_rsq = LDprofile$rsq, minRandL = 4, minRL = 25, X = 700),
               "X should be a numeric vector of length 2 e.g. c(100,200)",
               fixed=TRUE)
})

## Test the function with X supplied with too many numbers

test_that("Zalpha_expected fails with an X supplied with too many numbers", {

  expect_error(Zalpha_expected(pos = df$POS, cM = df$CM, ws  = 3000, LDprofile_cM_bins = LDprofile$cM_bin, LDprofile_rsq = LDprofile$rsq, minRandL = 4, minRL = 25, X = c(700,900,1100)),
               "X should be a numeric vector of length 2 e.g. c(100,200)",
               fixed=TRUE)
})

## Test the function with ws non-numeric

test_that("Zalpha_expected fails when ws is non-numeric", {

  expect_error(Zalpha_expected(pos = df$POS, cM = df$CM, ws  = "3000bp", LDprofile_cM_bins = LDprofile$cM_bin, LDprofile_rsq = LDprofile$rsq, minRandL = 4, minRL = 25, X = NULL),
               "ws must be a number greater than 0")
})

## Test the function with ws = 0

test_that("Zalpha_expected fails when ws is zero", {

  expect_error(Zalpha_expected(pos = df$POS, cM = df$CM, ws  = 0, LDprofile_cM_bins = LDprofile$cM_bin, LDprofile_rsq = LDprofile$rsq, minRandL = 4, minRL = 25, X = NULL),
               "ws must be a number greater than 0")
})

## Test the function with positions non-numeric

test_that("Zalpha_expected fails when pos is non-numeric", {

  expect_error(Zalpha_expected(pos = paste0(df$POS,"bp"), cM = df$CM, ws  = 3000, LDprofile_cM_bins = LDprofile$cM_bin, LDprofile_rsq = LDprofile$rsq, minRandL = 4, minRL = 25, X = NULL),
               "pos must be a numeric vector")
})

## Test the function with minLandR as non-numeric

test_that("Zalpha_expected fails when minLandR is non-numeric", {

  expect_error(Zalpha_expected(pos = df$POS, cM = df$CM, ws  = 3000, LDprofile_cM_bins = LDprofile$cM_bin, LDprofile_rsq = LDprofile$rsq, minRandL = "4snps", minRL = 25, X = NULL),
               "minRandL must be a number greater than or equal to 0")
})

## Test the function with minLandR as negative

test_that("Zalpha_expected fails when minLandR is negative", {

  expect_error(Zalpha_expected(pos = df$POS, cM = df$CM, ws  = 3000, LDprofile_cM_bins = LDprofile$cM_bin, LDprofile_rsq = LDprofile$rsq, minRandL = -1, minRL = 25, X = NULL),
               "minRandL must be a number greater than or equal to 0")
})

## Test the function with minLR as non-numeric

test_that("Zalpha_expected fails when minLR is non-numeric", {

  expect_error(Zalpha_expected(pos = df$POS, cM = df$CM, ws  = 3000, LDprofile_cM_bins = LDprofile$cM_bin, LDprofile_rsq = LDprofile$rsq, minRandL = 4, minRL = "25b", X = NULL),
               "minRL must be a number greater than or equal to 0")
})

## Test the function with minLR as negative

test_that("Zalpha_expected fails when minLR is negative", {

  expect_error(Zalpha_expected(pos = df$POS, cM = df$CM, ws  = 3000, LDprofile_cM_bins = LDprofile$cM_bin, LDprofile_rsq = LDprofile$rsq, minRandL = 4, minRL = -25, X = NULL),
               "minRL must be a number greater than or equal to 0")
})

## Test the warning that all values returned are NA by increasing minRL

test_that("Zalpha_expected warns about all NAs", {

  expect_warning(Zalpha_expected(pos = df$POS, cM = df$CM, ws  = 3000, LDprofile_cM_bins = LDprofile$cM_bin, LDprofile_rsq = LDprofile$rsq, minRandL = 4, minRL = 50, X = NULL),
                 "No Zalpha_expected values were calculated, try reducing minRandL and minRL or increasing the window size")
})

## Test the function with cMs non-numeric

test_that("Zalpha_expected fails when cM is non-numeric", {

  expect_error(Zalpha_expected(pos = df$POS, cM = paste0(df$CM,"cM"), ws  = 3000, LDprofile_cM_bins = LDprofile$cM_bin, LDprofile_rsq = LDprofile$rsq, minRandL = 4, minRL = 25, X = NULL),
               "cM must be a numeric vector")
})

## Test the function with cMs a different length to pos

test_that("Zalpha_expected fails when cM is a different length to pos", {

  expect_error(Zalpha_expected(pos = df$POS, cM = c(df$CM,1), ws  = 3000, LDprofile_cM_bins = LDprofile$cM_bin, LDprofile_rsq = LDprofile$rsq, minRandL = 4, minRL = 25, X = NULL),
               "The number of values in cM must equal the number of SNP locations given in pos")
})

## Test the function with LDprofile_cM_bins non-numeric

test_that("Zalpha_expected fails when LDprofile_cM_bins is non-numeric", {

  expect_error(Zalpha_expected(pos = df$POS, cM = df$CM, ws  = 3000, LDprofile_cM_bins = paste0(LDprofile$cM_bin,"cM"), LDprofile_rsq = LDprofile$rsq, minRandL = 4, minRL = 25, X = NULL),
               "LDprofile_cM_bins must be a numeric vector")
})

## Test the function with LDprofile_cM_bins not of equal size
tempLDprofile<-data.frame(
  cM_bin=c(seq(0,0.0048,0.0001),3),
  rsq=c(0.495714059385946,0.411014233619574,0.395532914859378,0.44354526861954,0.435550945945946,
        0.419534577303153,0.383410708498866,0.402439897670834,0.395945237932081,0.380436909495495,
        0.384229510773621,0.379494011054621,0.368118044626627,0.358753523652643,0.362330915047976,
        0.372142680938693,0.353703415234045,0.341173431307316,0.345578512726934,0.358159779825909,
        0.334687181997538,0.337442516960342,0.343721563062338,0.336509253287721,0.325271170690837,
        0.325488235114597,0.32078065970396,0.317594707821212,0.314613974963111,0.309774332543378,
        0.307619999017408,0.307004105181405,0.300279349768979,0.305505356875903,0.303179706053309,
        0.312783673753436,0.308540965611869,0.292690196360757,0.299428992380521,0.297144304197462,
        0.286701832971995,0.297894654636997,0.283984127549023,0.283253709389766,0.281330372503553,
        0.289052362087009,0.272730959781483,0.277399161038311,0.285764741944136,0.271195118636169)
)
test_that("Zalpha_expected fails when LDprofile_cM_bins are not of equal size", {

  expect_error(Zalpha_expected(pos = df$POS, cM = df$CM, ws  = 3000, LDprofile_cM_bins = tempLDprofile$cM_bin, LDprofile_rsq = LDprofile$rsq, minRandL = 4, minRL = 25, X = NULL),
               "LDprofile_cM_bins must be of equal size")
})

## Test the function with LDprofile_rsq non-numeric

test_that("Zalpha_expected fails when LDprofile_rsq is non-numeric", {

  expect_error(Zalpha_expected(pos = df$POS, cM = df$CM, ws  = 3000, LDprofile_cM_bins = LDprofile$cM_bin, LDprofile_rsq = paste0(LDprofile$rsq,"r"), minRandL = 4, minRL = 25, X = NULL),
               "LDprofile_rsq must be a numeric vector")
})

## Test the function with LDprofile_rsq having values outside (0,1)
tempLDprofile<-data.frame(
  cM_bin=seq(0,0.0049,0.0001),
  rsq=c(3,0.411014233619574,0.395532914859378,0.44354526861954,0.435550945945946,
        0.419534577303153,0.383410708498866,0.402439897670834,0.395945237932081,0.380436909495495,
        0.384229510773621,0.379494011054621,0.368118044626627,0.358753523652643,0.362330915047976,
        0.372142680938693,0.353703415234045,0.341173431307316,0.345578512726934,0.358159779825909,
        0.334687181997538,0.337442516960342,0.343721563062338,0.336509253287721,0.325271170690837,
        0.325488235114597,0.32078065970396,0.317594707821212,0.314613974963111,0.309774332543378,
        0.307619999017408,0.307004105181405,0.300279349768979,0.305505356875903,0.303179706053309,
        0.312783673753436,0.308540965611869,0.292690196360757,0.299428992380521,0.297144304197462,
        0.286701832971995,0.297894654636997,0.283984127549023,0.283253709389766,0.281330372503553,
        0.289052362087009,0.272730959781483,0.277399161038311,0.285764741944136,0.271195118636169)
)
test_that("Zalpha_expected fails when LDprofile_rsq contains values not between 0 and 1", {

  expect_error(Zalpha_expected(pos = df$POS, cM = df$CM, ws  = 3000, LDprofile_cM_bins = LDprofile$cM_bin, LDprofile_rsq = tempLDprofile$rsq, minRandL = 4, minRL = 25, X = NULL),
               "Values stored in LDprofile_rsq must be between 0 and 1")
})

## Test the function with LDprofile_cM_bins is a different length to LDprofile_rsq

test_that("Zalpha_expected fails when LDprofile_cM_bins and LDprofile_rsq are different lengths", {

  expect_error(Zalpha_expected(pos = df$POS, cM = df$CM, ws  = 3000, LDprofile_cM_bins = LDprofile$cM_bin, LDprofile_rsq = c(LDprofile$rsq,1), minRandL = 4, minRL = 25, X = NULL),
               "LDprofile_rsq must contain the same number of values as there are bins given in LDprofile_cM_bins")
})
