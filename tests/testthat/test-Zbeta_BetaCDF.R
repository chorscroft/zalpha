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
  Beta_a=c(0.225722015235706,0.227940311887665,0.230472949251006,0.235459534209174,0.234651102759893,
           0.250181212528427,0.242528673437283,0.241857801570361,0.243170537673904,0.249217742759865,
           0.248994063555682,0.250989836734379,0.265936661528539,0.261424850343736,0.262943726904153,
           0.272660652101016,0.269372770456792,0.269652818810360,0.272769051663880,0.271470095932534,
           0.280916111457310,0.277877819948299,0.276548748443875,0.281495922896857,0.279352276141887,
           0.284129685882995,0.277188108793857,0.278566874127134,0.289081444094226,0.285190045281596,
           0.290109830876961,0.284309968687703,0.288859114749642,0.299277331850620,0.285565479732833,
           0.284533597966335,0.298468126579221,0.297147927567662,0.292133986540405,0.282966380219017,
           0.298606708437661,0.294694288786295,0.299856879919194,0.302510385003750,0.306374031153195,
           0.290677944426393,0.302513945092534,0.304440730466093,0.294541534232723,0.306175724145037),
  Beta_b=c(0.138801044892172,0.188067327979005,0.201703485004672,0.175512844704570,0.178863105223407,
           0.219044703395120,0.229296709619857,0.212561622292768,0.218964966452713,0.243483004339105,
           0.239050159189516,0.246175476488997,0.284782373670816,0.285411205922151,0.290262258698151,
           0.290555194256511,0.306610200383973,0.323507472400767,0.324819589090210,0.312032319182944,
           0.354533038535116,0.348168933585225,0.336995503185295,0.361475463089895,0.361612563213451,
           0.378663668896915,0.376302977538757,0.378070928812591,0.407858651574400,0.408113158062142,
           0.423376460191435,0.407449239724648,0.435888337336229,0.451906955123291,0.431221465443690,
           0.397559355997338,0.449081727734466,0.474523994146667,0.446006870316418,0.424949954468304,
           0.492943986863418,0.458146004175747,0.495794738496397,0.518319998000701,0.534229325863624,
           0.458911554103514,0.543660553044380,0.536132041462793,0.482604851701327,0.565265304615674)
)

## test that Zbeta_BetaCDF is calculated correctly

test_that("Zbeta_BetaCDF calculates Zbeta_BetaCDF statistic correctly", {

  expect_equal(Zbeta_BetaCDF(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_Beta_a = LDprofile$Beta_a, LDprofile_Beta_b = LDprofile$Beta_b, minRandL = 4, minRL = 25, X = NULL),
               list(
                 position=c(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500),
                 Zbeta_BetaCDF=c(NA,NA,NA,NA,
                                 0.418016604924725,
                                 0.406297914560424,
                                 0.401067078642937,
                                 0.410829577389033,
                                 0.397509236264354,
                                 0.405347850057668,
                                 0.409781300588220,
                                 NA,NA,NA,NA)
               ),tolerance=0.0001)
})

## Test the function with a different window size

test_that("Zbeta_BetaCDF calculates Zbeta_BetaCDF statistic correctly with a different window size", {

  expect_equal(Zbeta_BetaCDF(pos = df$POS, ws  = 1100, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_Beta_a = LDprofile$Beta_a, LDprofile_Beta_b = LDprofile$Beta_b, minRandL = 4, minRL = 25, X = NULL),
               list(
                 position=c(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500),
                 Zbeta_BetaCDF=c(NA,NA,NA,NA,NA,
                                 0.422727696221137,
                                 0.368571388831177,
                                 0.391643594518566,
                                 0.398259343601550,
                                 0.408526848230388,
                                 NA,NA,NA,NA,NA)
               ),tolerance=0.0001)
})
## Test the function with a character matrix as x

test_that("Zbeta_BetaCDF calculates Zbeta_BetaCDF statistic correctly with character matrix", {

  df1<-df
  df1[df1==1]<-"A"
  df1[df1==2]<-"B"
  expect_equal(Zbeta_BetaCDF(pos = df$POS, ws  = 3000, x = as.matrix(df1[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_Beta_a = LDprofile$Beta_a, LDprofile_Beta_b = LDprofile$Beta_b, minRandL = 4, minRL = 25, X = NULL),
               list(
                 position=c(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500),
                 Zbeta_BetaCDF=c(NA,NA,NA,NA,
                                 0.418016604924725,
                                 0.406297914560424,
                                 0.401067078642937,
                                 0.410829577389033,
                                 0.397509236264354,
                                 0.405347850057668,
                                 0.409781300588220,
                                 NA,NA,NA,NA)
               ),tolerance=0.0001)
})

## Test the function with X supplied as a parameter

test_that("Zbeta_BetaCDF calculates Zbeta_BetaCDF statistic correctly with X supplied", {

  expect_equal(Zbeta_BetaCDF(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_Beta_a = LDprofile$Beta_a, LDprofile_Beta_b = LDprofile$Beta_b, minRandL = 4, minRL = 25, X = c(700,900)),
               list(
                 position=c(700,800,900),
                 Zbeta_BetaCDF=c(0.401067078642937,
                                 0.410829577389033,
                                 0.397509236264354)
               ),tolerance=0.0001)
})

## Test the function with X supplied as a parameter outside of the region defined in pos

test_that("Zbeta_BetaCDF fails with an X supplied outside of the region defined in pos", {

  expect_error(Zbeta_BetaCDF(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_Beta_a = LDprofile$Beta_a, LDprofile_Beta_b = LDprofile$Beta_b, minRandL = 4, minRL = 25, X = c(7000,9000)),
               "The region specified by X is outside the region contained in the pos vector")
})

## Test the function with X supplied as character

test_that("Zbeta_BetaCDF fails with an X supplied as a character", {

  expect_error(Zbeta_BetaCDF(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_Beta_a = LDprofile$Beta_a, LDprofile_Beta_b = LDprofile$Beta_b, minRandL = 4, minRL = 25, X = c("700bp","900bp")),
               "X should be a numeric vector of length 2 e.g. c(100,200)",
               fixed=TRUE)
})

## Test the function with X supplied as only one number

test_that("Zbeta_BetaCDF fails with an X supplied as only one number", {

  expect_error(Zbeta_BetaCDF(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_Beta_a = LDprofile$Beta_a, LDprofile_Beta_b = LDprofile$Beta_b, minRandL = 4, minRL = 25, X = 700),
               "X should be a numeric vector of length 2 e.g. c(100,200)",
               fixed=TRUE)
})

## Test the function with X supplied with too many numbers

test_that("Zbeta_BetaCDF fails with an X supplied with too many numbers", {

  expect_error(Zbeta_BetaCDF(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_Beta_a = LDprofile$Beta_a, LDprofile_Beta_b = LDprofile$Beta_b, minRandL = 4, minRL = 25, X = c(700,900,1100)),
               "X should be a numeric vector of length 2 e.g. c(100,200)",
               fixed=TRUE)
})

## Test the function with ws non-numeric

test_that("Zbeta_BetaCDF fails when ws is non-numeric", {

  expect_error(Zbeta_BetaCDF(pos = df$POS, ws  = "3000bp", x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_Beta_a = LDprofile$Beta_a, LDprofile_Beta_b = LDprofile$Beta_b, minRandL = 4, minRL = 25, X = NULL),
               "ws must be a number greater than 0")
})

## Test the function with ws = 0

test_that("Zbeta_BetaCDF fails when ws is zero", {

  expect_error(Zbeta_BetaCDF(pos = df$POS, ws  = 0, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_Beta_a = LDprofile$Beta_a, LDprofile_Beta_b = LDprofile$Beta_b, minRandL = 4, minRL = 25, X = NULL),
               "ws must be a number greater than 0")
})

## Test the function with positions non-numeric

test_that("Zbeta_BetaCDF fails when pos is non-numeric", {

  expect_error(Zbeta_BetaCDF(pos = paste0(df$POS,"bp"), ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_Beta_a = LDprofile$Beta_a, LDprofile_Beta_b = LDprofile$Beta_b, minRandL = 4, minRL = 25, X = NULL),
               "pos must be a numeric vector")
})

## Test the function with minLandR as non-numeric

test_that("Zbeta_BetaCDF fails when minLandR is non-numeric", {

  expect_error(Zbeta_BetaCDF(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_Beta_a = LDprofile$Beta_a, LDprofile_Beta_b = LDprofile$Beta_b, minRandL = "4snps", minRL = 25, X = NULL),
               "minRandL must be a number greater than or equal to 0")
})

## Test the function with minLandR as negative

test_that("Zbeta_BetaCDF fails when minLandR is negative", {

  expect_error(Zbeta_BetaCDF(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_Beta_a = LDprofile$Beta_a, LDprofile_Beta_b = LDprofile$Beta_b, minRandL = -1, minRL = 25, X = NULL),
               "minRandL must be a number greater than or equal to 0")
})

## Test the function with minLR as non-numeric

test_that("Zbeta_BetaCDF fails when minLR is non-numeric", {

  expect_error(Zbeta_BetaCDF(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_Beta_a = LDprofile$Beta_a, LDprofile_Beta_b = LDprofile$Beta_b, minRandL = 4, minRL = "25b", X = NULL),
               "minRL must be a number greater than or equal to 0")
})

## Test the function with minLR as negative

test_that("Zbeta_BetaCDF fails when minLR is negative", {

  expect_error(Zbeta_BetaCDF(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_Beta_a = LDprofile$Beta_a, LDprofile_Beta_b = LDprofile$Beta_b, minRandL = 4, minRL = -25, X = NULL),
               "minRL must be a number greater than or equal to 0")
})

## Test the warning that all values returned are NA by increasing minRL

test_that("Zbeta_BetaCDF warns about all NAs", {

  expect_warning(Zbeta_BetaCDF(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_Beta_a = LDprofile$Beta_a, LDprofile_Beta_b = LDprofile$Beta_b, minRandL = 4, minRL = 50, X = NULL),
                 "No Zbeta_BetaCDF values were calculated, try reducing minRandL and minRL or increasing the window size")
})

## Test the function with dists non-numeric

test_that("Zbeta_BetaCDF fails when dist is non-numeric", {

  expect_error(Zbeta_BetaCDF(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = paste0(df$dist,"dist"), LDprofile_bins = LDprofile$bin, LDprofile_Beta_a = LDprofile$Beta_a, LDprofile_Beta_b = LDprofile$Beta_b, minRandL = 4, minRL = 25, X = NULL),
               "dist must be a numeric vector")
})

## Test the function with dists a different length to pos

test_that("Zbeta_BetaCDF fails when dist is a different length to pos", {

  expect_error(Zbeta_BetaCDF(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = c(df$dist,1), LDprofile_bins = LDprofile$bin, LDprofile_Beta_a = LDprofile$Beta_a, LDprofile_Beta_b = LDprofile$Beta_b, minRandL = 4, minRL = 25, X = NULL),
               "The number of values in dist must equal the number of SNP locations given in pos")
})

## Test the function with LDprofile_bins non-numeric

test_that("Zbeta_BetaCDF fails when LDprofile_bins is non-numeric", {

  expect_error(Zbeta_BetaCDF(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = paste0(LDprofile$bin,"dist"), LDprofile_Beta_a = LDprofile$Beta_a, LDprofile_Beta_b = LDprofile$Beta_b, minRandL = 4, minRL = 25, X = NULL),
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
test_that("Zbeta_BetaCDF fails when LDprofile_bins are not of equal size", {

  expect_error(Zbeta_BetaCDF(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = tempLDprofile$bin, LDprofile_Beta_a = LDprofile$Beta_a, LDprofile_Beta_b = LDprofile$Beta_b, minRandL = 4, minRL = 25, X = NULL),
               "LDprofile_bins must be of equal size")
})

## Test the function with LDprofile_Beta_a non-numeric

test_that("Zbeta_BetaCDF fails when LDprofile_Beta_a is non-numeric", {

  expect_error(Zbeta_BetaCDF(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_Beta_a = paste0(LDprofile$Beta_a,"r"), LDprofile_Beta_b = LDprofile$Beta_b, minRandL = 4, minRL = 25, X = NULL),
               "LDprofile_Beta_a must be a numeric vector")
})

## Test the function with LDprofile_Beta_b non-numeric

test_that("Zbeta_BetaCDF fails when LDprofile_Beta_b is non-numeric", {

  expect_error(Zbeta_BetaCDF(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_Beta_a = LDprofile$Beta_a, LDprofile_Beta_b = paste0(LDprofile$Beta_b,"sd"), minRandL = 4, minRL = 25, X = NULL),
               "LDprofile_Beta_b must be a numeric vector")
})

## Test the function with x not a matrix

test_that("Zbeta_BetaCDF fails when x is not a matrix", {

  expect_error(Zbeta_BetaCDF(pos = df$POS, ws  = 3000, x = df[,3:7], dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_Beta_a = LDprofile$Beta_a, LDprofile_Beta_b = LDprofile$Beta_b, minRandL = 4, minRL = 25, X = NULL),
               "x must be a matrix")
})

## Test the function with x not having the correct amount of rows

test_that("Zbeta_BetaCDF fails when the number of rows in x is not equal to the length of pos", {

  expect_error(Zbeta_BetaCDF(pos = df$POS, ws  = 3000, x = t(as.matrix(df[,3:7])), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_Beta_a = LDprofile$Beta_a, LDprofile_Beta_b = LDprofile$Beta_b, minRandL = 4, minRL = 25, X = NULL),
               "The number of rows in x must equal the number of SNP locations given in pos")
})

## Test the function with a SNP having only one allele

test_that("Zbeta_BetaCDF fails when a SNP has only one allele", {

  df1<-df
  df1[1,3:7]<-1
  expect_error(Zbeta_BetaCDF(pos = df1$POS, ws  = 3000, x = as.matrix(df1[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_Beta_a = LDprofile$Beta_a, LDprofile_Beta_b = LDprofile$Beta_b, minRandL = 4, minRL = 25, X = NULL),
               "SNPs must all be biallelic")
})

## Test the function with a SNP having more than two alleles

test_that("Zbeta_BetaCDF fails when a SNP has more than two alleles", {

  df1<-df
  df1[1,7]<-3
  expect_error(Zbeta_BetaCDF(pos = df1$POS, ws  = 3000, x = as.matrix(df1[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_Beta_a = LDprofile$Beta_a, LDprofile_Beta_b = LDprofile$Beta_b, minRandL = 4, minRL = 25, X = NULL),
               "SNPs must all be biallelic")
})

## Test the function with LDprofile_bins is a different length to LDprofile_Beta_a

test_that("Zbeta_BetaCDF fails when LDprofile_bins and LDprofile_Beta_a are different lengths", {

  expect_error(Zbeta_BetaCDF(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_Beta_a = c(LDprofile$Beta_a,1), LDprofile_Beta_b = LDprofile$Beta_b, minRandL = 4, minRL = 25, X = NULL),
               "LDprofile_Beta_a must contain the same number of values as there are bins given in LDprofile_bins")
})

## Test the function with LDprofile_bins as a different length to LDprofile_Beta_b

test_that("Zbeta_BetaCDF fails when LDprofile_bins and LDprofile_Beta_b are different lengths", {

  expect_error(Zbeta_BetaCDF(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_Beta_a = LDprofile$Beta_a, LDprofile_Beta_b = c(LDprofile$Beta_b,1), minRandL = 4, minRL = 25, X = NULL),
               "LDprofile_Beta_b must contain the same number of values as there are bins given in LDprofile_bins")
})

## test that Zbeta_BetaCDF works with a missing value
df1<-df
df1$C1[15]<-NA
test_that("Zbeta_BetaCDF calculates Zbeta_BetaCDF statistic correctly with missing value", {

  expect_equal(Zbeta_BetaCDF(pos = df1$POS, ws  = 3000, x = as.matrix(df1[,3:7]), dist = df1$dist, LDprofile_bins = LDprofile$bin, LDprofile_Beta_a = LDprofile$Beta_a, LDprofile_Beta_b = LDprofile$Beta_b, minRandL = 4, minRL = 25, X = NULL),
               list(
                 position=c(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500),
                 Zbeta_BetaCDF=c(NA,NA,NA,NA,
                                 0.412985057460831,
                                 0.400927939537124,
                                 0.395021661102834,
                                 0.415693421938135,
                                 0.407239348282580,
                                 0.414765530721706,
                                 0.422608554537899,
                                 NA,NA,NA,NA)
               ),tolerance=0.0001)
})
