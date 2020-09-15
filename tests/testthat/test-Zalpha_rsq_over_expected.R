df<-data.frame(
  SNP=c("SNP1","SNP2","SNP3","SNP4","SNP5","SNP6","SNP7","SNP8","SNP9","SNP10","SNP11","SNP12","SNP13","SNP14","SNP15"),
  POS=c(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500),
  C1=c(1,1,2,1,2,1,1,2,1,2,1,2,1,1,1),
  C2=c(2,2,1,2,1,2,1,2,2,2,1,2,1,1,2),
  C3=c(2,1,2,2,2,1,1,2,2,1,2,2,1,1,2),
  C4=c(1,1,2,1,2,2,1,1,1,1,1,2,2,2,2),
  C5=c(1,1,2,1,2,1,2,1,1,1,1,1,2,1,1),
  dist=c(0,0.00101,0.00123,0.00207,0.00218,0.00223,0.00235,0.00251,0.00272,0.00289,0.00304,0.00316,0.00335,0.00345,0.00374)
)
LDprofile<-data.frame(
  bin=seq(0,0.0049,0.0001),
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
## test that Zalpha_rsq_over_expected is calculated correctly

test_that("Zalpha_rsq_over_expected calculates Zalpha_rsq_over_expected statistic correctly", {

  expect_equal(Zalpha_rsq_over_expected(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, minRandL = 4, minRL = 25, X = NULL),
               list(
                 position=c(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500),
                 Zalpha_rsq_over_expected=c(NA,NA,NA,NA,
                                            1.099392777415350,
                                            1.175423932431810,
                                            0.973208670564282,
                                            0.767152145198758,
                                            0.691051112331195,
                                            0.728230393808371,
                                            0.737847444047460,
                                   NA,NA,NA,NA)
               ),tolerance=0.0001)
})

## Test the function with a different window size

test_that("Zalpha_rsq_over_expected calculates Zalpha_rsq_over_expected statistic correctly with a different window size", {

  expect_equal(Zalpha_rsq_over_expected(pos = df$POS, ws  = 1100, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, minRandL = 4, minRL = 25, X = NULL),
               list(
                 position=c(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500),
                 Zalpha_rsq_over_expected=c(NA,NA,NA,NA,NA,
                                            1.120554321638970,
                                            0.961086582792753,
                                            0.658704548108426,
                                            0.499640647493653,
                                            0.478776705516631,
                                   NA,NA,NA,NA,NA)
               ),tolerance=0.0001)
})
## Test the function with a character matrix as x

test_that("Zalpha_rsq_over_expected calculates Zalpha_rsq_over_expected statistic correctly with character matrix", {

  df1<-df
  df1[df1==1]<-"A"
  df1[df1==2]<-"B"
  expect_equal(Zalpha_rsq_over_expected(pos = df$POS, ws  = 3000, x = as.matrix(df1[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, minRandL = 4, minRL = 25, X = NULL),
               list(
                 position=c(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500),
                 Zalpha_rsq_over_expected=c(NA,NA,NA,NA,
                                            1.099392777415350,
                                            1.175423932431810,
                                            0.973208670564282,
                                            0.767152145198758,
                                            0.691051112331195,
                                            0.728230393808371,
                                            0.737847444047460,
                                            NA,NA,NA,NA)
               ),tolerance=0.0001)
})

## Test the function with X supplied as a parameter

test_that("Zalpha_rsq_over_expected calculates Zalpha_rsq_over_expected statistic correctly with X supplied", {

  expect_equal(Zalpha_rsq_over_expected(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, minRandL = 4, minRL = 25, X = c(700,900)),
               list(
                 position=c(700,800,900),
                 Zalpha_rsq_over_expected=c(0.973208670564282,
                                            0.767152145198758,
                                            0.691051112331195)
               ),tolerance=0.0001)
})

## Test the function with X supplied as a parameter outside of the region defined in pos

test_that("Zalpha_rsq_over_expected fails with an X supplied outside of the region defined in pos", {

  expect_error(Zalpha_rsq_over_expected(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, minRandL = 4, minRL = 25, X = c(7000,9000)),
               "The region specified by X is outside the region contained in the pos vector")
})

## Test the function with X supplied as character

test_that("Zalpha_rsq_over_expected fails with an X supplied as a character", {

  expect_error(Zalpha_rsq_over_expected(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, minRandL = 4, minRL = 25, X = c("700bp","900bp")),
               "X should be a numeric vector of length 2 e.g. c(100,200)",
               fixed=TRUE)
})

## Test the function with X supplied as only one number

test_that("Zalpha_rsq_over_expected fails with an X supplied as only one number", {

  expect_error(Zalpha_rsq_over_expected(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, minRandL = 4, minRL = 25, X = 700),
               "X should be a numeric vector of length 2 e.g. c(100,200)",
               fixed=TRUE)
})

## Test the function with X supplied with too many numbers

test_that("Zalpha_rsq_over_expected fails with an X supplied with too many numbers", {

  expect_error(Zalpha_rsq_over_expected(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, minRandL = 4, minRL = 25, X = c(700,900,1100)),
               "X should be a numeric vector of length 2 e.g. c(100,200)",
               fixed=TRUE)
})

## Test the function with ws non-numeric

test_that("Zalpha_rsq_over_expected fails when ws is non-numeric", {

  expect_error(Zalpha_rsq_over_expected(pos = df$POS, ws  = "3000bp", x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, minRandL = 4, minRL = 25, X = NULL),
               "ws must be a number greater than 0")
})

## Test the function with ws = 0

test_that("Zalpha_rsq_over_expected fails when ws is zero", {

  expect_error(Zalpha_rsq_over_expected(pos = df$POS, ws  = 0, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, minRandL = 4, minRL = 25, X = NULL),
               "ws must be a number greater than 0")
})

## Test the function with positions non-numeric

test_that("Zalpha_rsq_over_expected fails when pos is non-numeric", {

  expect_error(Zalpha_rsq_over_expected(pos = paste0(df$POS,"bp"), ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, minRandL = 4, minRL = 25, X = NULL),
               "pos must be a numeric vector")
})

## Test the function with minLandR as non-numeric

test_that("Zalpha_rsq_over_expected fails when minLandR is non-numeric", {

  expect_error(Zalpha_rsq_over_expected(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, minRandL = "4snps", minRL = 25, X = NULL),
               "minRandL must be a number greater than or equal to 0")
})

## Test the function with minLandR as negative

test_that("Zalpha_rsq_over_expected fails when minLandR is negative", {

  expect_error(Zalpha_rsq_over_expected(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, minRandL = -1, minRL = 25, X = NULL),
               "minRandL must be a number greater than or equal to 0")
})

## Test the function with minLR as non-numeric

test_that("Zalpha_rsq_over_expected fails when minLR is non-numeric", {

  expect_error(Zalpha_rsq_over_expected(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, minRandL = 4, minRL = "25b", X = NULL),
               "minRL must be a number greater than or equal to 0")
})

## Test the function with minLR as negative

test_that("Zalpha_rsq_over_expected fails when minLR is negative", {

  expect_error(Zalpha_rsq_over_expected(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, minRandL = 4, minRL = -25, X = NULL),
               "minRL must be a number greater than or equal to 0")
})

## Test the warning that all values returned are NA by increasing minRL

test_that("Zalpha_rsq_over_expected warns about all NAs", {

  expect_warning(Zalpha_rsq_over_expected(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, minRandL = 4, minRL = 50, X = NULL),
                 "No Zalpha_rsq_over_expected values were calculated, try reducing minRandL and minRL or increasing the window size")
})

## Test the function with dists non-numeric

test_that("Zalpha_rsq_over_expected fails when dist is non-numeric", {

  expect_error(Zalpha_rsq_over_expected(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = paste0(df$dist,"dist"), LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, minRandL = 4, minRL = 25, X = NULL),
               "dist must be a numeric vector")
})

## Test the function with dists a different length to pos

test_that("Zalpha_rsq_over_expected fails when dist is a different length to pos", {

  expect_error(Zalpha_rsq_over_expected(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = c(df$dist,1), LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, minRandL = 4, minRL = 25, X = NULL),
               "The number of values in dist must equal the number of SNP locations given in pos")
})

## Test the function with LDprofile_bins non-numeric

test_that("Zalpha_rsq_over_expected fails when LDprofile_bins is non-numeric", {

  expect_error(Zalpha_rsq_over_expected(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = paste0(LDprofile$bin,"dist"), LDprofile_rsq = LDprofile$rsq, minRandL = 4, minRL = 25, X = NULL),
               "LDprofile_bins must be a numeric vector")
})

## Test the function with LDprofile_bins not of equal size
tempLDprofile<-data.frame(
  bin=c(seq(0,0.0048,0.0001),3),
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
test_that("Zalpha_rsq_over_expected fails when LDprofile_bins are not of equal size", {

  expect_error(Zalpha_rsq_over_expected(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = tempLDprofile$bin, LDprofile_rsq = LDprofile$rsq, minRandL = 4, minRL = 25, X = NULL),
               "LDprofile_bins must be of equal size")
})

## Test the function with LDprofile_rsq non-numeric

test_that("Zalpha_rsq_over_expected fails when LDprofile_rsq is non-numeric", {

  expect_error(Zalpha_rsq_over_expected(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = paste0(LDprofile$rsq,"r"), minRandL = 4, minRL = 25, X = NULL),
               "LDprofile_rsq must be a numeric vector")
})

## Test the function with LDprofile_rsq having values outside (0,1)
tempLDprofile<-data.frame(
  bin=seq(0,0.0049,0.0001),
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
test_that("Zalpha_rsq_over_expected fails when LDprofile_rsq contains values not between 0 and 1", {

  expect_error(Zalpha_rsq_over_expected(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = tempLDprofile$rsq, minRandL = 4, minRL = 25, X = NULL),
               "Values stored in LDprofile_rsq must be between 0 and 1")
})

## Test the function with x not a matrix

test_that("Zalpha_rsq_over_expected fails when x is not a matrix", {

  expect_error(Zalpha_rsq_over_expected(pos = df$POS, ws  = 3000, x = df[,3:7], dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, minRandL = 4, minRL = 25, X = NULL),
               "x must be a matrix")
})

## Test the function with x not having the correct amount of rows

test_that("Zalpha_rsq_over_expected fails when the number of rows in x is not equal to the length of pos", {

  expect_error(Zalpha_rsq_over_expected(pos = df$POS, ws  = 3000, x = t(as.matrix(df[,3:7])), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, minRandL = 4, minRL = 25, X = NULL),
               "The number of rows in x must equal the number of SNP locations given in pos")
})

## Test the function with a SNP having only one allele

test_that("Zalpha_rsq_over_expected fails when a SNP has only one allele", {

  df1<-df
  df1[1,3:7]<-1
  expect_error(Zalpha_rsq_over_expected(pos = df1$POS, ws  = 3000, x = as.matrix(df1[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, minRandL = 4, minRL = 25, X = NULL),
               "SNPs must all be biallelic")
})

## Test the function with a SNP having more than two alleles

test_that("Zalpha_rsq_over_expected fails when a SNP has more than two alleles", {

  df1<-df
  df1[1,7]<-3
  expect_error(Zalpha_rsq_over_expected(pos = df1$POS, ws  = 3000, x = as.matrix(df1[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, minRandL = 4, minRL = 25, X = NULL),
               "SNPs must all be biallelic")
})

## Test the function with LDprofile_bins is a different length to LDprofile_rsq

test_that("Zalpha_rsq_over_expected fails when LDprofile_bins and LDprofile_rsq are different lengths", {

  expect_error(Zalpha_rsq_over_expected(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = c(LDprofile$rsq,1), minRandL = 4, minRL = 25, X = NULL),
               "LDprofile_rsq must contain the same number of values as there are bins given in LDprofile_bins")
})

## test that Zalpha_rsq_over_expected works with a missing value
df1<-df
df1$C1[15]<-NA
test_that("Zalpha_rsq_over_expected calculates Zalpha_rsq_over_expected statistic correctly with missing value", {

  expect_equal(Zalpha_rsq_over_expected(pos = df1$POS, ws  = 3000, x = as.matrix(df1[,3:7]), dist = df1$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, minRandL = 4, minRL = 25, X = NULL),
               list(
                 position=c(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500),
                 Zalpha_rsq_over_expected=c(NA,NA,NA,NA,
                                            1.144893254873100,
                                            1.236446352382830,
                                            1.020556368226920,
                                            0.810519367025429,
                                            0.761404523358315,
                                            0.823237169574888,
                                            0.907695973532459,
                                            NA,NA,NA,NA)
               ),tolerance=0.0001)
})

## test what happens when the biggest bin is bigger than the max_dist in the LDprofile

df1<-df
df1$dist[10:15]<-df1$dist[10:15]+0.1
test_that("Zalpha_rsq_over_expected calculates Zalpha_rsq_over_expected statistic correctly when biggest bin is bigger than LDprofile", {

  expect_equal(Zalpha_rsq_over_expected(pos = df1$POS, ws  = 3000, x = as.matrix(df1[,3:7]), dist = df1$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, minRandL = 4, minRL = 25, X = NULL),
               list(
                 position=c(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500),
                 Zalpha_rsq_over_expected=c(NA,NA,NA,NA,
                                                1.19167472080619,
                                                1.27358283231798,
                                                1.05911709750039,
                                                0.813967663774974,
                                                0.691051112331195,
                                                0.728230393808371,
                                                0.76047497606951,
                                                NA,NA,NA,NA)
               ),tolerance=0.0001)
})
