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

## test that everything is calculated correctly given all parameters

test_that("create_LDprofile calculates the LD profile correctly", {

  expect_equal(create_LDprofile(dist = df$dist, x = as.matrix(df[,3:7]), bin_size = 0.001, max_dist = 0.005, beta_params = TRUE),
      data.frame(
     bin=c(0,0.001,0.002,0.003,0.004),
     rsq=c(0.285622427983539,0.280913978494624,0.263888888888889,0.319444444444444,NA),
     sd=c(0.270862044573862,0.201905775929377,0.321786617161322,0.142318760638328,NA),
     Beta_a=c(0.619957744381906,1.125028692019340,0.635410044952769,3.941019442363900,NA),
     Beta_b=c(1.062459890834270,2.446706389704430,1.149319432462400,8.454825333760550,NA),
     n=c(54,31,15,5,0)
      ),tolerance=0.0001)
})

## Test the function with a different max_dist

test_that("create_LDprofile calculates the LD profile correctly with a different max_dist", {

  expect_equal(create_LDprofile(dist = df$dist, x = as.matrix(df[,3:7]), bin_size = 0.001, max_dist = 0.003, beta_params = TRUE),
      data.frame(
     bin=c(0,0.001,0.002),
     rsq=c(0.285622427983539,0.280913978494624,0.263888888888889),
     sd=c(0.270862044573862,0.201905775929377,0.321786617161322),
     Beta_a=c(0.619957744381906,1.125028692019340,0.635410044952769),
     Beta_b=c(1.062459890834270,2.446706389704430,1.149319432462400),
     n=c(54,31,15)
      ),tolerance=0.0001)
})
## Test the function with no max_dist given

test_that("create_LDprofile calculates the LD profile correctly with no max_dist supplied", {

  expect_equal(create_LDprofile(dist = df$dist, x = as.matrix(df[,3:7]), bin_size = 0.001, beta_params = TRUE),
      data.frame(
     bin=c(0,0.001,0.002,0.003),
     rsq=c(0.285622427983539,0.280913978494624,0.263888888888889,0.319444444444444),
     sd=c(0.270862044573862,0.201905775929377,0.321786617161322,0.142318760638328),
     Beta_a=c(0.619957744381906,1.125028692019340,0.635410044952769,3.941019442363900),
     Beta_b=c(1.062459890834270,2.446706389704430,1.149319432462400,8.454825333760550),
     n=c(54,31,15,5)
      ),tolerance=0.0001)
})


## Test the function with a different bin_size

test_that("create_LDprofile calculates the LD profile correctly with a different bin size", {

  expect_equal(create_LDprofile(dist = df$dist, x = as.matrix(df[,3:7]), bin_size = 0.0005, beta_params = TRUE),
      data.frame(
     bin=c(0,0.0005,0.001,0.0015,0.002,0.0025,0.003,0.0035),
     rsq=c(0.238505747126437,0.340277777777778,0.283459595959596,0.274691358024691,0.215277777777778,0.361111111111111,0.288194444444444,0.444444444444445),
     sd=c(0.211600341827602,0.322468326753589,0.220527561550113,0.158590157477369,0.293808275018141,0.387895557,0.14316339,NA),
     Beta_a=c(0.916070145958307,0.637072700079744,1.046576044485340,1.909812912830260,0.775059123115346,1.088198634018290,3.789877096116000,NA),
     Beta_b=c(2.326350552394540,0.872215477086822,2.166981335251990,5.166454170350760,1.748740564135290,1.488374161884570,9.367197007381050,NA),
     n=c(29,25,22,9,10,5,4,1)
      ),tolerance=0.0001)
})

## Test the function with beta_params not specified

test_that("create_LDprofile calculates the LD profile correctly with beta_params not specified", {

  expect_equal(create_LDprofile(dist = df$dist, x = as.matrix(df[,3:7]), bin_size = 0.001, max_dist = 0.005),
      data.frame(
     bin=c(0,0.001,0.002,0.003,0.004),
     rsq=c(0.285622427983539,0.280913978494624,0.263888888888889,0.319444444444444,NA),
     sd=c(0.270862044573862,0.201905775929377,0.321786617161322,0.142318760638328,NA),
     Beta_a=c(NA,NA,NA,NA,NA),
     Beta_b=c(NA,NA,NA,NA,NA),
     n=c(54,31,15,5,0)
      ),tolerance=0.0001)
})

## Test the function with a character matrix as x

df1<-df
df1[df1==1]<-"A"
df1[df1==2]<-"B"
test_that("create_LDprofile calculates the LD profile correctly with character matrix", {

  expect_equal(create_LDprofile(dist = df1$dist, x = as.matrix(df1[,3:7]), bin_size = 0.001, max_dist = 0.005, beta_params = TRUE),
      data.frame(
     bin=c(0,0.001,0.002,0.003,0.004),
     rsq=c(0.285622427983539,0.280913978494624,0.263888888888889,0.319444444444444,NA),
     sd=c(0.270862044573862,0.201905775929377,0.321786617161322,0.142318760638328,NA),
     Beta_a=c(0.619957744381906,1.125028692019340,0.635410044952769,3.941019442363900,NA),
     Beta_b=c(1.062459890834270,2.446706389704430,1.149319432462400,8.454825333760550,NA),
     n=c(54,31,15,5,0)
      ),tolerance=0.0001)
})

## Test all the checks

## Test the function with dists non-numeric

test_that("create_LDprofile fails when dist is non-numeric", {

  expect_error(create_LDprofile(dist = paste0(df$dist,"dist"), x = as.matrix(df[,3:7]), bin_size = 0.001, max_dist = 0.005, beta_params = TRUE),
      "dist must be a numeric vector or list of numeric vectors")
})

## Test the function with dists not a vector

test_that("create_LDprofile fails when dist is not a vector", {

  expect_error(create_LDprofile(dist = as.matrix(df$dist), x = as.matrix(df[,3:7]), bin_size = 0.001, max_dist = 0.005, beta_params = TRUE),
      "dist must be a numeric vector or list of numeric vectors")
})


## Test the function with x not a matrix

test_that("create_LDprofile fails when x is not a matrix", {

  expect_error(create_LDprofile(dist = df$dist, x = df$C1, bin_size = 0.001, max_dist = 0.005, beta_params = TRUE),
      "x must be a matrix or list of matrices")
})

## Test the function with x not having the correct amount of rows

test_that("create_LDprofile fails when the number of rows in x is not equal to the length of pos", {

  expect_error(create_LDprofile(dist = df$dist, x = t(as.matrix(df[,3:7])), bin_size = 0.001, max_dist = 0.005, beta_params = TRUE),
      "The number of rows in x must equal the number of SNP genetic distances given in the corresponding dist")
})

## Test the function with a SNP having only one allele

test_that("create_LDprofile fails when a SNP has only one allele", {

  df1<-df
  df1[1,3:7]<-1
  expect_error(create_LDprofile(dist = df1$dist, x = as.matrix(df1[,3:7]), bin_size = 0.001, max_dist = 0.005, beta_params = TRUE),
      "SNPs must all be biallelic")
})

## Test the function with a SNP having more than two alleles

test_that("create_LDprofile fails when a SNP has more than two alleles", {

  df1<-df
  df1[1,7]<-3
  expect_error(create_LDprofile(dist = df1$dist, x = as.matrix(df1[,3:7]), bin_size = 0.001, max_dist = 0.005, beta_params = TRUE),
      "SNPs must all be biallelic")
})

## Test the function with dists and x have a different number of elements

test_that("create_LDprofile fails when dist and x are different lengths", {

  expect_error(create_LDprofile(dist = list(df$dist,df$dist), x = as.matrix(df[,3:7]), bin_size = 0.001, max_dist = 0.005, beta_params = TRUE),
      "dist and x should contain the same number of elements")
})


## Test the function with bin_size as non-numeric

test_that("create_LDprofile fails when bin_size is non-numeric", {

  expect_error(create_LDprofile(dist = df$dist, x = as.matrix(df[,3:7]), bin_size = "0.001cM", max_dist = 0.005, beta_params = TRUE),
      "bin_size must be a number greater than 0")
})

## Test the function with bin_size as negative

test_that("create_LDprofile fails when bin_size is negative", {

  expect_error(create_LDprofile(dist = df$dist, x = as.matrix(df[,3:7]), bin_size = -1, max_dist = 0.005, beta_params = TRUE),
      "bin_size must be a number greater than 0")
})

## Test the function with max_dist as non-numeric

test_that("create_LDprofile fails when max_dist is non-numeric", {

  expect_error(create_LDprofile(dist = df$dist, x = as.matrix(df[,3:7]), bin_size = 0.001, max_dist = "0.005cM", beta_params = TRUE),
      "max_dist must be a number greater than 0")
})

## Test the function with max_dist as negative

test_that("create_LDprofile fails when max_dist is negative", {

  expect_error(create_LDprofile(dist = df$dist, x = as.matrix(df[,3:7]), bin_size = 0.001, max_dist = -1, beta_params = TRUE),
      "max_dist must be a number greater than 0")
})

## Test the function with beta_params not logical

test_that("create_LDprofile fails when beta_params is not logical", {

  expect_error(create_LDprofile(dist = df$dist, x = as.matrix(df[,3:7]), bin_size = 0.001, max_dist = 0.005, beta_params = 1),
      "beta_params must be TRUE or FALSE")
})

## Test the function with missing value in x
df1<-df
df1$C1[15]<-NA
test_that("create_LDprofile calculates the LD profile correctly with missing value in x", {

  expect_equal(create_LDprofile(dist = df1$dist, x = as.matrix(df1[,3:7]), bin_size = 0.001, max_dist = 0.005, beta_params = TRUE),
      data.frame(
     bin=c(0,0.001,0.002,0.003,0.004),
     rsq=c(0.302340534979424,0.298387096774194,0.256481481481481,0.297222222222222,NA),
     sd=c(0.285102946872176,0.231243979115867,0.324764476229430,0.125615767273278,NA),
     Beta_a=c(0.606846782972070,0.932888189465068,0.616307318328496,4.643089841249520,NA),
     Beta_b=c(0.942600409939103,1.692309859072240,1.139111486418360,11.012925076592600,NA),
     n=c(54,31,15,5,0)
      ),tolerance=0.0001)
})

## Test the function with missing value in dist
df1<-df
df1$dist[5]<-NA
test_that("create_LDprofile calculates the LD profile correctly with missing value in dist", {

  expect_equal(create_LDprofile(dist = df1$dist, x = as.matrix(df1[,3:7]), bin_size = 0.001, max_dist = 0.005, beta_params = TRUE),
      data.frame(
      bin=c(0,0.001,0.002,0.003,0.004),
      rsq=c(0.279320987654321,0.270833333333333,0.255952380952381,0.319444444444444,NA),
      sd=c(0.268871132983145,0.155580832042849,0.332406755525893,0.142318760638328,NA),
      Beta_a=c(0.634602938184746,1.570771368863860,0.611627386753874,3.941019442363900,NA),
      Beta_b=c(1.133751642734910,4.401600723631340,1.116161049339830,8.454825333760550,NA),
      n=c(45,27,14,5,0)
      ),tolerance=0.0001)
})


## Test the function with fitdistrplus package not loaded and beta_params = TRUE
#Tested manually

## Test the function when beta estimation doesn't work the first try
x1<-matrix(c(1,1,1,2,1,2,1,1,1,2,1,2,2,2,2,
             1,1,1,2,1,2,2,2,2,1,1,2,1,1,2,
             1,2,1,2,1,1,1,2,1,1,2,1,2,2,2,
             1,2,1,2,1,1,2,1,2,2,1,2,2,1,2,
             1,2,1,2,1,2,1,2,1,2,1,1,2,1,2,
             1,2,2,1,1,2,1,2,2,1,2,1,1,1,2,
             1,1,2,2,1,1,2,1,2,2,2,1,1,1,2,
             1,1,2,1,1,1,2,2,1,2,2,1,1,1,2,
             1,1,1,1,1,2,2,2,2,1,2,1,2,2,2,
             1,2,2,2,2,1,2,1,2,1,2,2,2,1,2,
             1,2,1,1,1,1,2,1,1,1,2,1,1,2,2,
             1,1,1,2,2,1,2,1,2,1,1,2,2,1,2,
             1,1,1,1,1,1,1,2,2,2,1,2,2,1,2,
             1,2,2,1,2,1,1,1,2,1,2,1,2,1,2,
             1,1,2,1,1,1,2,2,2,1,1,1,2,1,2),byrow=TRUE,nrow=15)
test_that("create_LDprofile calculates the LD profile correctly when beta calculation fails on first attempt", {
  expect_equal(create_LDprofile(dist = rep(0,15), x = x1, bin_size = 0.001, max_dist = 0.001, beta_params = TRUE),
               data.frame(
                 bin=c(0),
                 rsq=c(0.0754751721676325),
                 sd=c(0.0928525038706074),
                 Beta_a=c(0.21155096488681),
                 Beta_b=c(2.8704915799235),
                 n=c(105)
               ),tolerance=0.0001)
})

## Test the function when beta estimation doesn't work on the second try
