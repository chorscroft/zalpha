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
        0.289052362087009,0.272730959781483,0.277399161038311,0.285764741944136,0.271195118636169),
  sd=c(0.427539638,0.410029691,0.40440292,0.415539937,0.413150465,
       0.403737933,0.394455514,0.400046196,0.39634811,0.388658831,
       0.391039303,0.388333172,0.377804309,0.375516702,0.37593021,
       0.376077309,0.371055975,0.366087782,0.365512538,0.370981607,
       0.35783996,0.360382871,0.362003257,0.357105337,0.353497616,
       0.352204789,0.353083084,0.348153839,0.343274196,0.343452807,
       0.340024865,0.342316376,0.337826872,0.337178573,0.339827834,
       0.344695499,0.33702497,0.331228682,0.335701618,0.337432909,
       0.327185588,0.332486895,0.323867319,0.323176228,0.319542129,
       0.329335786,0.317310776,0.32009561,0.325978889,0.31411119)
)
## test that Zbeta_Zscore is calculated correctly

test_that("Zbeta_Zscore calculates Zbeta_Zscore statistic correctly", {

  expect_equal(Zbeta_Zscore(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, LDprofile_sd = LDprofile$sd, minRandL = 4, minRL = 25, X = NULL),
               list(
                 position=c(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500),
                 Zbeta_Zscore=c(NA,NA,NA,NA,
                                -0.249622295287423,
                                -0.298443278340203,
                                -0.303226918699104,
                                -0.278319517790552,
                                -0.324158658191404,
                                -0.298557987824230,
                                -0.278980061934724,
                                 NA,NA,NA,NA)
               ))
})

## Test the function with a different window size

test_that("Zbeta_Zscore calculates Zbeta_Zscore statistic correctly with a different window size", {

  expect_equal(Zbeta_Zscore(pos = df$POS, ws  = 1100, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, LDprofile_sd = LDprofile$sd, minRandL = 4, minRL = 25, X = NULL),
               list(
                 position=c(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500),
                 Zbeta_Zscore=c(NA,NA,NA,NA,NA,
                                -0.209270924830470,
                                -0.410246297817478,
                                -0.344123580516901,
                                -0.308094601136825,
                                -0.264295466759451,
                                 NA,NA,NA,NA,NA)
               ))
})
## Test the function with a character matrix as x

test_that("Zbeta_Zscore calculates Zbeta_Zscore statistic correctly with character matrix", {

  df1<-df
  df1[df1==1]<-"A"
  df1[df1==2]<-"B"
  expect_equal(Zbeta_Zscore(pos = df$POS, ws  = 3000, x = as.matrix(df1[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, LDprofile_sd = LDprofile$sd, minRandL = 4, minRL = 25, X = NULL),
               list(
                 position=c(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500),
                 Zbeta_Zscore=c(NA,NA,NA,NA,
                                -0.249622295287423,
                                -0.298443278340203,
                                -0.303226918699104,
                                -0.278319517790552,
                                -0.324158658191404,
                                -0.298557987824230,
                                -0.278980061934724,
                                 NA,NA,NA,NA)
               ))
})

## Test the function with X supplied as a parameter

test_that("Zbeta_Zscore calculates Zbeta_Zscore statistic correctly with X supplied", {

  expect_equal(Zbeta_Zscore(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, LDprofile_sd = LDprofile$sd, minRandL = 4, minRL = 25, X = c(700,900)),
               list(
                 position=c(700,800,900),
                 Zbeta_Zscore=c(-0.303226918699104,
                                -0.278319517790552,
                                -0.324158658191404)
               ))
})

## Test the function with X supplied as a parameter outside of the region defined in pos

test_that("Zbeta_Zscore fails with an X supplied outside of the region defined in pos", {

  expect_error(Zbeta_Zscore(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, LDprofile_sd = LDprofile$sd, minRandL = 4, minRL = 25, X = c(7000,9000)),
               "The region specified by X is outside the region contained in the pos vector")
})

## Test the function with X supplied as character

test_that("Zbeta_Zscore fails with an X supplied as a character", {

  expect_error(Zbeta_Zscore(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, LDprofile_sd = LDprofile$sd, minRandL = 4, minRL = 25, X = c("700bp","900bp")),
               "X should be a numeric vector of length 2 e.g. c(100,200)",
               fixed=TRUE)
})

## Test the function with X supplied as only one number

test_that("Zbeta_Zscore fails with an X supplied as only one number", {

  expect_error(Zbeta_Zscore(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, LDprofile_sd = LDprofile$sd, minRandL = 4, minRL = 25, X = 700),
               "X should be a numeric vector of length 2 e.g. c(100,200)",
               fixed=TRUE)
})

## Test the function with X supplied with too many numbers

test_that("Zbeta_Zscore fails with an X supplied with too many numbers", {

  expect_error(Zbeta_Zscore(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, LDprofile_sd = LDprofile$sd, minRandL = 4, minRL = 25, X = c(700,900,1100)),
               "X should be a numeric vector of length 2 e.g. c(100,200)",
               fixed=TRUE)
})

## Test the function with ws non-numeric

test_that("Zbeta_Zscore fails when ws is non-numeric", {

  expect_error(Zbeta_Zscore(pos = df$POS, ws  = "3000bp", x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, LDprofile_sd = LDprofile$sd, minRandL = 4, minRL = 25, X = NULL),
               "ws must be a number greater than 0")
})

## Test the function with ws = 0

test_that("Zbeta_Zscore fails when ws is zero", {

  expect_error(Zbeta_Zscore(pos = df$POS, ws  = 0, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, LDprofile_sd = LDprofile$sd, minRandL = 4, minRL = 25, X = NULL),
               "ws must be a number greater than 0")
})

## Test the function with positions non-numeric

test_that("Zbeta_Zscore fails when pos is non-numeric", {

  expect_error(Zbeta_Zscore(pos = paste0(df$POS,"bp"), ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, LDprofile_sd = LDprofile$sd, minRandL = 4, minRL = 25, X = NULL),
               "pos must be a numeric vector")
})

## Test the function with minLandR as non-numeric

test_that("Zbeta_Zscore fails when minLandR is non-numeric", {

  expect_error(Zbeta_Zscore(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, LDprofile_sd = LDprofile$sd, minRandL = "4snps", minRL = 25, X = NULL),
               "minRandL must be a number greater than or equal to 0")
})

## Test the function with minLandR as negative

test_that("Zbeta_Zscore fails when minLandR is negative", {

  expect_error(Zbeta_Zscore(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, LDprofile_sd = LDprofile$sd, minRandL = -1, minRL = 25, X = NULL),
               "minRandL must be a number greater than or equal to 0")
})

## Test the function with minLR as non-numeric

test_that("Zbeta_Zscore fails when minLR is non-numeric", {

  expect_error(Zbeta_Zscore(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, LDprofile_sd = LDprofile$sd, minRandL = 4, minRL = "25b", X = NULL),
               "minRL must be a number greater than or equal to 0")
})

## Test the function with minLR as negative

test_that("Zbeta_Zscore fails when minLR is negative", {

  expect_error(Zbeta_Zscore(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, LDprofile_sd = LDprofile$sd, minRandL = 4, minRL = -25, X = NULL),
               "minRL must be a number greater than or equal to 0")
})

## Test the warning that all values returned are NA by increasing minRL

test_that("Zbeta_Zscore warns about all NAs", {

  expect_warning(Zbeta_Zscore(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, LDprofile_sd = LDprofile$sd, minRandL = 4, minRL = 50, X = NULL),
                 "No Zbeta_Zscore values were calculated, try reducing minRandL and minRL or increasing the window size")
})

## Test the function with dists non-numeric

test_that("Zbeta_Zscore fails when dist is non-numeric", {

  expect_error(Zbeta_Zscore(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = paste0(df$dist,"dist"), LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, LDprofile_sd = LDprofile$sd, minRandL = 4, minRL = 25, X = NULL),
               "dist must be a numeric vector")
})

## Test the function with dists a different length to pos

test_that("Zbeta_Zscore fails when dist is a different length to pos", {

  expect_error(Zbeta_Zscore(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = c(df$dist,1), LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, LDprofile_sd = LDprofile$sd, minRandL = 4, minRL = 25, X = NULL),
               "The number of values in dist must equal the number of SNP locations given in pos")
})

## Test the function with LDprofile_bins non-numeric

test_that("Zbeta_Zscore fails when LDprofile_bins is non-numeric", {

  expect_error(Zbeta_Zscore(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = paste0(LDprofile$bin,"dist"), LDprofile_rsq = LDprofile$rsq, LDprofile_sd = LDprofile$sd, minRandL = 4, minRL = 25, X = NULL),
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
test_that("Zbeta_Zscore fails when LDprofile_bins are not of equal size", {

  expect_error(Zbeta_Zscore(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = tempLDprofile$bin, LDprofile_rsq = LDprofile$rsq, LDprofile_sd = LDprofile$sd, minRandL = 4, minRL = 25, X = NULL),
               "LDprofile_bins must be of equal size")
})

## Test the function with LDprofile_rsq non-numeric

test_that("Zbeta_Zscore fails when LDprofile_rsq is non-numeric", {

  expect_error(Zbeta_Zscore(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = paste0(LDprofile$rsq,"r"), LDprofile_sd = LDprofile$sd, minRandL = 4, minRL = 25, X = NULL),
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
test_that("Zbeta_Zscore fails when LDprofile_rsq contains values not between 0 and 1", {

  expect_error(Zbeta_Zscore(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = tempLDprofile$rsq, LDprofile_sd = LDprofile$sd, minRandL = 4, minRL = 25, X = NULL),
               "Values stored in LDprofile_rsq must be between 0 and 1")
})

## Test the function with LDprofile_sd non-numeric

test_that("Zbeta_Zscore fails when LDprofile_sd is non-numeric", {

  expect_error(Zbeta_Zscore(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, LDprofile_sd = paste0(LDprofile$sd,"sd"), minRandL = 4, minRL = 25, X = NULL),
               "LDprofile_sd must be a numeric vector")
})

## Test the function with x not a matrix

test_that("Zbeta_Zscore fails when x is not a matrix", {

  expect_error(Zbeta_Zscore(pos = df$POS, ws  = 3000, x = df[,3:7], dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, LDprofile_sd = LDprofile$sd, minRandL = 4, minRL = 25, X = NULL),
               "x must be a matrix")
})

## Test the function with x not having the correct amount of rows

test_that("Zbeta_Zscore fails when the number of rows in x is not equal to the length of pos", {

  expect_error(Zbeta_Zscore(pos = df$POS, ws  = 3000, x = t(as.matrix(df[,3:7])), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, LDprofile_sd = LDprofile$sd, minRandL = 4, minRL = 25, X = NULL),
               "The number of rows in x must equal the number of SNP locations given in pos")
})

## Test the function with a SNP having only one allele

test_that("Zbeta_Zscore fails when a SNP has only one allele", {

  df1<-df
  df1[1,3:7]<-1
  expect_error(Zbeta_Zscore(pos = df1$POS, ws  = 3000, x = as.matrix(df1[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, LDprofile_sd = LDprofile$sd, minRandL = 4, minRL = 25, X = NULL),
               "SNPs must all be biallelic")
})

## Test the function with a SNP having more than two alleles

test_that("Zbeta_Zscore fails when a SNP has more than two alleles", {

  df1<-df
  df1[1,7]<-3
  expect_error(Zbeta_Zscore(pos = df1$POS, ws  = 3000, x = as.matrix(df1[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, LDprofile_sd = LDprofile$sd, minRandL = 4, minRL = 25, X = NULL),
               "SNPs must all be biallelic")
})

## Test the function with LDprofile_bins is a different length to LDprofile_rsq

test_that("Zbeta_Zscore fails when LDprofile_bins and LDprofile_rsq are different lengths", {

  expect_error(Zbeta_Zscore(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = c(LDprofile$rsq,1), LDprofile_sd = LDprofile$sd, minRandL = 4, minRL = 25, X = NULL),
               "LDprofile_rsq must contain the same number of values as there are bins given in LDprofile_bins")
})

## Test the function with LDprofile_bins as a different length to LDprofile_sd

test_that("Zbeta_Zscore fails when LDprofile_bins and LDprofile_sd are different lengths", {

  expect_error(Zbeta_Zscore(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, LDprofile_sd = c(LDprofile$sd,1), minRandL = 4, minRL = 25, X = NULL),
               "LDprofile_sd must contain the same number of values as there are bins given in LDprofile_bins")
})

## test that Zbeta_Zscore works with a missing value
df1<-df
df1$C1[15]<-NA
test_that("Zbeta_Zscore calculates Zbeta_Zscore statistic correctly with missing value", {

  expect_equal(Zbeta_Zscore(pos = df1$POS, ws  = 3000, x = as.matrix(df1[,3:7]), dist = df1$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, LDprofile_sd = LDprofile$sd, minRandL = 4, minRL = 25, X = NULL),
               list(
                 position=c(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500),
                 Zbeta_Zscore=c(NA,NA,NA,NA,
                                -0.273427430734028,
                                -0.322886149060994,
                                -0.332297266849268,
                                -0.272829789964421,
                                -0.301705254123099,
                                -0.280921980725937,
                                -0.253883231792888,
                                NA,NA,NA,NA)
               ))
})
