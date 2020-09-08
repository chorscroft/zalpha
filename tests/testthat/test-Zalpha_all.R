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
       0.329335786,0.317310776,0.32009561,0.325978889,0.31411119),
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

## test that everything is calculated correctly given all parameters

test_that("Zalpha_all calculates statistics correctly", {

  expect_equal(Zalpha_all(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, LDprofile_sd = LDprofile$sd, LDprofile_Beta_a = LDprofile$Beta_a, LDprofile_Beta_b = LDprofile$Beta_b, minRandL = 4, minRL = 25, X = NULL),
               list(
                 position=c(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500),
                 LR=c(0,13,24,33,40,45,48,49,48,45,40,33,24,13,0),
                 L_plus_R=c(91,78,67,58,51,46,43,42,43,46,51,58,67,78,91),
                 Zalpha_expected=c(NA,NA,NA,NA,0.390457304338967,0.392014054942343,0.397339546324536,0.398874728980465,0.400715520018796,0.401718327864356,0.399526703832832,NA,NA,NA,NA),
                 Zbeta_expected=c(NA,NA,NA,NA,0.350404001770323,0.357372168444253,0.360647397851440,0.361951221318345,0.362988750761055,0.364752550557575,0.366343120440209,NA,NA,NA,NA),
                 Zalpha=c(NA,NA,NA,NA,((3+1/2)/6+(11+41/144)/45)/2,((6+1/4)/10+(9+41/48)/36)/2,((7+31/72)/15+(7+13/48)/28)/2,((8+17/144)/21+(4+7/16)/21)/2,((9+131/144)/28+(2+13/16)/15)/2,((13+97/144)/36+(1+121/144)/10)/2,((15+25/48)/45+(1+55/144)/6)/2,NA,NA,NA,NA),
                 Zbeta=c(NA,NA,NA,NA,(10+5/18)/40,(10+35/36)/45,(11+103/144)/48,(12+73/144)/49,(11+83/144)/48,(11+17/48)/45,(10+65/144)/40,NA,NA,NA,NA),
                 Zalpha_rsq_over_expected=c(NA,NA,NA,NA,1.099392777415350,1.175423932431810,0.973208670564282,0.767152145198758,0.691051112331195,0.728230393808371,0.737847444047460,NA,NA,NA,NA),
                 Zalpha_log_rsq_over_expected=c(NA,NA,NA,NA, -0.129685954728232, -0.086676970494806,-0.207059296326176,  -0.317561061974900, -0.355497049708911, -0.345346264723395, -0.340516242835506,NA,NA,NA,NA),
                 Zbeta_rsq_over_expected=c(NA,NA,NA,NA,0.741001103918817,0.693986983022872,0.691113499222754,0.716092979819490,0.663090587779225,0.690647913157800,0.710787633520743,NA,NA,NA,NA),
                 Zbeta_log_rsq_over_expected=c(NA,NA,NA,NA,-0.310756050486805,-0.323333602724241,-0.349639186187507,-0.339617245857445,-0.355799393664026,-0.313469388842524,-0.304918066115218,NA,NA,NA,NA),
                 Zalpha_Zscore=c(NA,NA,NA,NA,0.083669516080429,0.160231324773080,-0.038223960672508,-0.240917156668399,-0.316679396483571,-0.285401890680237,-0.270109843273656,NA,NA,NA,NA),
                 Zbeta_Zscore=c(NA,NA,NA,NA,-0.249622295287423,-0.298443278340203,-0.303226918699104,-0.278319517790552,-0.324158658191404,-0.298557987824230,-0.278980061934724,NA,NA,NA,NA),
                 Zalpha_BetaCDF=c(NA,NA,NA,NA,0.498582271149287,0.525043954366501,0.461984646273587,0.409097872871755,0.388665184297966,0.396797658473664,0.393583363498442,NA,NA,NA,NA),
                 Zbeta_BetaCDF=c(NA,NA,NA,NA,0.418016604924725,0.406297914560424,0.401067078642937,0.410829577389033,0.397509236264354,0.405347850057668,0.409781300588220,NA,NA,NA,NA)
               ),tolerance=0.0001)
})

## Test the function with a different window size

test_that("Zalpha_all calculates the statistics correctly with a different window size", {

  expect_equal(Zalpha_all(pos = df$POS, ws  = 1100, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, LDprofile_sd = LDprofile$sd, LDprofile_Beta_a = LDprofile$Beta_a, LDprofile_Beta_b = LDprofile$Beta_b, minRandL = 4, minRL = 25, X = NULL),
               list(
                 position=c(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500),
                 LR=c(0,5,10,15,20,25,25,25,25,25,20,15,10,5,0),
                 L_plus_R=c(10,10,11,13,16,20,20,20,20,20,16,13,11,10,10),
                 Zalpha_expected=c(NA,NA,NA,NA,NA,0.397641067838740,0.407771198285253,0.412830128247348,0.419512588688985,0.421536156959405,NA,NA,NA,NA,NA),
                 Zbeta_expected=c(NA,NA,NA,NA,NA,0.369069693255773,0.381250676471648,0.385948372025294,0.388272359451172,0.387629142774876,NA,NA,NA,NA,NA),
                 Zalpha=c(NA,NA,NA,NA,NA,((6+1/4)/10+(2+19/48)/10)/2,((5+5/18)/10+(2+19/48)/10)/2,((2+71/72)/10+(2+19/48)/10)/2,((2+3/16)/10+(2+7/144)/10)/2,((2+3/16)/10+(1+121/144)/10)/2,NA,NA,NA,NA,NA),
                 Zbeta=c(NA,NA,NA,NA,NA,(7+11/72)/25,(5+5/9)/25,(6+41/144)/25,(6+101/144)/25,(7+17/144)/25,NA,NA,NA,NA,NA),
                 Zalpha_rsq_over_expected=c(NA,NA,NA,NA,NA,1.120554321638970,0.961086582792753,0.658704548108426,0.499640647493653,0.478776705516631,NA,NA,NA,NA,NA),
                 Zalpha_log_rsq_over_expected=c(NA,NA,NA,NA,NA,-0.097503940903089,-0.186457320769577,-0.347591637479257,-0.432639322094622,-0.477272952569813,NA,NA,NA,NA,NA),
                 Zbeta_rsq_over_expected=c(NA,NA,NA,NA,NA,0.792819358044316,0.584585947609949,0.653705466293006,0.686521835223347,0.732759533350118,NA,NA,NA,NA,NA),
                 Zbeta_log_rsq_over_expected=c(NA,NA,NA,NA,NA,-0.298569863567948,-0.437982657397525,-0.409664183061625,-0.364550829001460,-0.315576925527679,NA,NA,NA,NA,NA),
                 Zalpha_Zscore=c(NA,NA,NA,NA,NA,0.104386407929017,-0.049338365858586,-0.347772802475024,-0.507186471020645,-0.535949238344209,NA,NA,NA,NA,NA),
                 Zbeta_Zscore=c(NA,NA,NA,NA,NA,-0.209270924830470,-0.410246297817478,-0.344123580516901,-0.308094601136825,-0.264295466759451,NA,NA,NA,NA,NA),
                 Zalpha_BetaCDF=c(NA,NA,NA,NA,NA,0.504981991655003,0.459025835290601,0.371462610998979,0.325273987256311,0.319208662374142,NA,NA,NA,NA,NA),
                 Zbeta_BetaCDF=c(NA,NA,NA,NA,NA,0.422727696221137,0.368571388831177,0.391643594518566,0.398259343601550,0.408526848230388,NA,NA,NA,NA,NA)
               ),tolerance=0.0001)
})

## Test the function with a character matrix as x

test_that("Zalpha_all calculates all statistics correctly with character matrix", {

  df1<-df
  df1[df1==1]<-"A"
  df1[df1==2]<-"B"
  expect_equal(Zalpha_all(pos = df$POS, ws  = 3000, x = as.matrix(df1[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, LDprofile_sd = LDprofile$sd, LDprofile_Beta_a = LDprofile$Beta_a, LDprofile_Beta_b = LDprofile$Beta_b, minRandL = 4, minRL = 25, X = NULL),
               list(
                 position=c(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500),
                 LR=c(0,13,24,33,40,45,48,49,48,45,40,33,24,13,0),
                 L_plus_R=c(91,78,67,58,51,46,43,42,43,46,51,58,67,78,91),
                 Zalpha_expected=c(NA,NA,NA,NA,0.390457304338967,0.392014054942343,0.397339546324536,0.398874728980465,0.400715520018796,0.401718327864356,0.399526703832832,NA,NA,NA,NA),
                 Zbeta_expected=c(NA,NA,NA,NA,0.350404001770323,0.357372168444253,0.360647397851440,0.361951221318345,0.362988750761055,0.364752550557575,0.366343120440209,NA,NA,NA,NA),
                 Zalpha=c(NA,NA,NA,NA,((3+1/2)/6+(11+41/144)/45)/2,((6+1/4)/10+(9+41/48)/36)/2,((7+31/72)/15+(7+13/48)/28)/2,((8+17/144)/21+(4+7/16)/21)/2,((9+131/144)/28+(2+13/16)/15)/2,((13+97/144)/36+(1+121/144)/10)/2,((15+25/48)/45+(1+55/144)/6)/2,NA,NA,NA,NA),
                 Zbeta=c(NA,NA,NA,NA,(10+5/18)/40,(10+35/36)/45,(11+103/144)/48,(12+73/144)/49,(11+83/144)/48,(11+17/48)/45,(10+65/144)/40,NA,NA,NA,NA),
                 Zalpha_rsq_over_expected=c(NA,NA,NA,NA,1.099392777415350,1.175423932431810,0.973208670564282,0.767152145198758,0.691051112331195,0.728230393808371,0.737847444047460,NA,NA,NA,NA),
                 Zalpha_log_rsq_over_expected=c(NA,NA,NA,NA, -0.129685954728232, -0.086676970494806,-0.207059296326176,  -0.317561061974900, -0.355497049708911, -0.345346264723395, -0.340516242835506,NA,NA,NA,NA),
                 Zbeta_rsq_over_expected=c(NA,NA,NA,NA,0.741001103918817,0.693986983022872,0.691113499222754,0.716092979819490,0.663090587779225,0.690647913157800,0.710787633520743,NA,NA,NA,NA),
                 Zbeta_log_rsq_over_expected=c(NA,NA,NA,NA,-0.310756050486805,-0.323333602724241,-0.349639186187507,-0.339617245857445,-0.355799393664026,-0.313469388842524,-0.304918066115218,NA,NA,NA,NA),
                 Zalpha_Zscore=c(NA,NA,NA,NA,0.083669516080429,0.160231324773080,-0.038223960672508,-0.240917156668399,-0.316679396483571,-0.285401890680237,-0.270109843273656,NA,NA,NA,NA),
                 Zbeta_Zscore=c(NA,NA,NA,NA,-0.249622295287423,-0.298443278340203,-0.303226918699104,-0.278319517790552,-0.324158658191404,-0.298557987824230,-0.278980061934724,NA,NA,NA,NA),
                 Zalpha_BetaCDF=c(NA,NA,NA,NA,0.498582271149287,0.525043954366501,0.461984646273587,0.409097872871755,0.388665184297966,0.396797658473664,0.393583363498442,NA,NA,NA,NA),
                 Zbeta_BetaCDF=c(NA,NA,NA,NA,0.418016604924725,0.406297914560424,0.401067078642937,0.410829577389033,0.397509236264354,0.405347850057668,0.409781300588220,NA,NA,NA,NA)
               ),tolerance=0.0001)
})

## Test the function with X supplied as a parameter

test_that("Zalpha_all calculates the statistics correctly with X supplied", {

  expect_equal(Zalpha_all(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, LDprofile_sd = LDprofile$sd, LDprofile_Beta_a = LDprofile$Beta_a, LDprofile_Beta_b = LDprofile$Beta_b, minRandL = 4, minRL = 25, X = c(700,900)),
               list(
                  position=c(700,800,900),
                  LR=c(48,49,48),
                  L_plus_R=c(43,42,43),
                  Zalpha_expected=c(0.397339546324536,0.398874728980465,0.400715520018796),
                  Zbeta_expected=c(0.360647397851440,0.361951221318345,0.362988750761055),
                  Zalpha=c(((7+31/72)/15+(7+13/48)/28)/2,((8+17/144)/21+(4+7/16)/21)/2,((9+131/144)/28+(2+13/16)/15)/2),
                  Zbeta=c((11+103/144)/48,(12+73/144)/49,(11+83/144)/48),
                  Zalpha_rsq_over_expected=c(0.973208670564282,0.767152145198758,0.691051112331195),
                  Zalpha_log_rsq_over_expected=c(-0.207059296326176,-0.317561061974900,-0.355497049708911),
                  Zbeta_rsq_over_expected=c(0.691113499222754,0.716092979819490,0.663090587779225),
                  Zbeta_log_rsq_over_expected=c(-0.349639186187507,-0.339617245857445,-0.355799393664026),
                  Zalpha_Zscore=c(-0.038223960672508, -0.240917156668399,-0.316679396483571),
                  Zbeta_Zscore=c(-0.303226918699104,-0.278319517790552,-0.324158658191404),
                  Zalpha_BetaCDF=c(0.461984646273587,0.409097872871755,0.388665184297966),
                  Zbeta_BetaCDF=c(0.401067078642937,0.410829577389033,0.397509236264354)
               ),tolerance=0.0001)
})

## Test the function with X supplied as a parameter outside of the region defined in pos

test_that("Zalpha_all fails with an X supplied outside of the region defined in pos", {

  expect_error(Zalpha_all(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, LDprofile_sd = LDprofile$sd, LDprofile_Beta_a = LDprofile$Beta_a, LDprofile_Beta_b = LDprofile$Beta_b, minRandL = 4, minRL = 25, X = c(7000,9000)),
               "The region specified by X is outside the region contained in the pos vector")
})

## Test the function with X supplied as character

test_that("Zalpha_all fails with an X supplied as a character", {

  expect_error(Zalpha_all(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, LDprofile_sd = LDprofile$sd, LDprofile_Beta_a = LDprofile$Beta_a, LDprofile_Beta_b = LDprofile$Beta_b, minRandL = 4, minRL = 25, X = c("700bp","900bp")),
               "X should be a numeric vector of length 2 e.g. c(100,200)",
               fixed=TRUE)
})

## Test the function with X supplied as only one number

test_that("Zalpha_all fails with an X supplied as only one number", {

  expect_error(Zalpha_all(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, LDprofile_sd = LDprofile$sd, LDprofile_Beta_a = LDprofile$Beta_a, LDprofile_Beta_b = LDprofile$Beta_b, minRandL = 4, minRL = 25, X = 700),
               "X should be a numeric vector of length 2 e.g. c(100,200)",
               fixed=TRUE)
})

## Test the function with X supplied with too many numbers

test_that("Zalpha_all fails with an X supplied with too many numbers", {

  expect_error(Zalpha_all(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, LDprofile_sd = LDprofile$sd, LDprofile_Beta_a = LDprofile$Beta_a, LDprofile_Beta_b = LDprofile$Beta_b, minRandL = 4, minRL = 25, X = c(700,900,1100)),
               "X should be a numeric vector of length 2 e.g. c(100,200)",
               fixed=TRUE)
})

## Test the function with ws non-numeric

test_that("Zalpha_all fails when ws is non-numeric", {

  expect_error(Zalpha_all(pos = df$POS, ws  = "3000bp", x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, LDprofile_sd = LDprofile$sd, LDprofile_Beta_a = LDprofile$Beta_a, LDprofile_Beta_b = LDprofile$Beta_b, minRandL = 4, minRL = 25, X = NULL),
               "ws must be a number greater than 0")
})

## Test the function with ws = 0

test_that("Zalpha_all fails when ws is zero", {

  expect_error(Zalpha_all(pos = df$POS, ws  = 0, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, LDprofile_sd = LDprofile$sd, LDprofile_Beta_a = LDprofile$Beta_a, LDprofile_Beta_b = LDprofile$Beta_b, minRandL = 4, minRL = 25, X = NULL),
               "ws must be a number greater than 0")
})

## Test the function with positions non-numeric

test_that("Zalpha_all fails when pos is non-numeric", {

  expect_error(Zalpha_all(pos = paste0(df$POS,"bp"), ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, LDprofile_sd = LDprofile$sd, LDprofile_Beta_a = LDprofile$Beta_a, LDprofile_Beta_b = LDprofile$Beta_b, minRandL = 4, minRL = 25, X = NULL),
               "pos must be a numeric vector")
})

## Test the function with minLandR as non-numeric

test_that("Zalpha_all fails when minLandR is non-numeric", {

  expect_error(Zalpha_all(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, LDprofile_sd = LDprofile$sd, LDprofile_Beta_a = LDprofile$Beta_a, LDprofile_Beta_b = LDprofile$Beta_b, minRandL = "4snps", minRL = 25, X = NULL),
               "minRandL must be a number greater than or equal to 0")
})

## Test the function with minLandR as negative

test_that("Zalpha_all fails when minLandR is negative", {

  expect_error(Zalpha_all(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, LDprofile_sd = LDprofile$sd, LDprofile_Beta_a = LDprofile$Beta_a, LDprofile_Beta_b = LDprofile$Beta_b, minRandL = -1, minRL = 25, X = NULL),
               "minRandL must be a number greater than or equal to 0")
})

## Test the function with minLR as non-numeric

test_that("Zalpha_all fails when minLR is non-numeric", {

  expect_error(Zalpha_all(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, LDprofile_sd = LDprofile$sd, LDprofile_Beta_a = LDprofile$Beta_a, LDprofile_Beta_b = LDprofile$Beta_b, minRandL = 4, minRL = "25b", X = NULL),
               "minRL must be a number greater than or equal to 0")
})

## Test the function with minLR as negative

test_that("Zalpha_all fails when minLR is negative", {

  expect_error(Zalpha_all(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, LDprofile_sd = LDprofile$sd, LDprofile_Beta_a = LDprofile$Beta_a, LDprofile_Beta_b = LDprofile$Beta_b, minRandL = 4, minRL = -25, X = NULL),
               "minRL must be a number greater than or equal to 0")
})

## Test the warning that all values returned are NA by increasing minRL

test_that("Zalpha_all warns about all NAs", {

  expect_warning(Zalpha_all(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, LDprofile_sd = LDprofile$sd, LDprofile_Beta_a = LDprofile$Beta_a, LDprofile_Beta_b = LDprofile$Beta_b, minRandL = 4, minRL = 50, X = NULL),
                 "No statistics were calculated, try reducing minRandL and minRL or increasing the window size")
})

## Test the function with dists non-numeric

test_that("Zalpha_all fails when dist is non-numeric", {

  expect_error(Zalpha_all(pos = df$POS, ws= 3000, x = as.matrix(df[,3:7]), dist = paste0(df$dist,"dist"), LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, LDprofile_sd = LDprofile$sd, LDprofile_Beta_a = LDprofile$Beta_a, LDprofile_Beta_b = LDprofile$Beta_b, minRandL = 4, minRL = 25, X = NULL),
               "dist must be a numeric vector")
})

## Test the function with dists a different length to pos

test_that("Zalpha_all fails when dist is a different length to pos", {

  expect_error(Zalpha_all(pos = df$POS, ws = 3000, x = as.matrix(df[,3:7]), dist = c(df$dist,1), LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, LDprofile_sd = LDprofile$sd, LDprofile_Beta_a = LDprofile$Beta_a, LDprofile_Beta_b = LDprofile$Beta_b, minRandL = 4, minRL = 25, X = NULL),
               "The number of values in dist must equal the number of SNP locations given in pos")
})

## Test the function with LDprofile_bins non-numeric

test_that("Zalpha_all fails when LDprofile_bins is non-numeric", {

  expect_error(Zalpha_all(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = paste0(LDprofile$bin,"dist"), LDprofile_rsq = LDprofile$rsq, LDprofile_sd = LDprofile$sd, LDprofile_Beta_a = LDprofile$Beta_a, LDprofile_Beta_b = LDprofile$Beta_b, minRandL = 4, minRL = 25, X = NULL),
               "LDprofile_bins must be a numeric vector")
})

## Test the function with LDprofile_bins not of equal size
tempLDprofile<-data.frame(
  bin=c(seq(0,0.0048,0.0001),3)
)
test_that("Zalpha_all fails when LDprofile_bins are not of equal size", {

  expect_error(Zalpha_all(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = tempLDprofile$bin, LDprofile_rsq = LDprofile$rsq, LDprofile_sd = LDprofile$sd, LDprofile_Beta_a = LDprofile$Beta_a, LDprofile_Beta_b = LDprofile$Beta_b, minRandL = 4, minRL = 25, X = NULL),
               "LDprofile_bins must be of equal size")
})

## Test the function with LDprofile_rsq non-numeric

test_that("Zalpha_all fails when LDprofile_rsq is non-numeric", {

  expect_error(Zalpha_all(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = paste0(LDprofile$rsq,"r"), LDprofile_sd = LDprofile$sd, LDprofile_Beta_a = LDprofile$Beta_a, LDprofile_Beta_b = LDprofile$Beta_b, minRandL = 4, minRL = 25, X = NULL),
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
test_that("Zalpha_all fails when LDprofile_rsq contains values not between 0 and 1", {

  expect_error(Zalpha_all(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = tempLDprofile$rsq, LDprofile_sd = LDprofile$sd, LDprofile_Beta_a = LDprofile$Beta_a, LDprofile_Beta_b = LDprofile$Beta_b, minRandL = 4, minRL = 25, X = NULL),
               "Values stored in LDprofile_rsq must be between 0 and 1")
})

## Test the function with LDprofile_sd non-numeric

test_that("Zalpha_all fails when LDprofile_sd is non-numeric", {

  expect_error(Zalpha_all(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, LDprofile_sd = paste0(LDprofile$sd,"sd"), LDprofile_Beta_a = LDprofile$Beta_a, LDprofile_Beta_b = LDprofile$Beta_b, minRandL = 4, minRL = 25, X = NULL),
               "LDprofile_sd must be a numeric vector")
})

## Test the function with x not a matrix

test_that("Zalpha_all fails when x is not a matrix", {

  expect_error(Zalpha_all(pos = df$POS, ws  = 3000, x = df[,3:7], dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, LDprofile_sd = LDprofile$sd, LDprofile_Beta_a = LDprofile$Beta_a, LDprofile_Beta_b = LDprofile$Beta_b, minRandL = 4, minRL = 25, X = NULL),
               "x must be a matrix")
})

## Test the function with x not having the correct amount of rows

test_that("Zalpha_all fails when the number of rows in x is not equal to the length of pos", {

  expect_error(Zalpha_all(pos = df$POS, ws  = 3000, x = t(as.matrix(df[,3:7])), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, LDprofile_sd = LDprofile$sd, LDprofile_Beta_a = LDprofile$Beta_a, LDprofile_Beta_b = LDprofile$Beta_b, minRandL = 4, minRL = 25, X = NULL),
               "The number of rows in x must equal the number of SNP locations given in pos")
})

## Test the function with a SNP having only one allele

test_that("Zalpha_all fails when a SNP has only one allele", {

  df1<-df
  df1[1,3:7]<-1
  expect_error(Zalpha_all(pos = df1$POS, ws  = 3000, x = as.matrix(df1[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, LDprofile_sd = LDprofile$sd, LDprofile_Beta_a = LDprofile$Beta_a, LDprofile_Beta_b = LDprofile$Beta_b, minRandL = 4, minRL = 25, X = NULL),
               "SNPs must all be biallelic")
})

## Test the function with a SNP having more than two alleles

test_that("Zalpha_all fails when a SNP has more than two alleles", {

  df1<-df
  df1[1,7]<-3
  expect_error(Zalpha_all(pos = df1$POS, ws  = 3000, x = as.matrix(df1[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, LDprofile_sd = LDprofile$sd, LDprofile_Beta_a = LDprofile$Beta_a, LDprofile_Beta_b = LDprofile$Beta_b, minRandL = 4, minRL = 25, X = NULL),
               "SNPs must all be biallelic")
})

## Test the function with LDprofile_bins is a different length to LDprofile_rsq

test_that("Zalpha_all fails when LDprofile_bins and LDprofile_rsq are different lengths", {

  expect_error(Zalpha_all(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = c(LDprofile$rsq,1), LDprofile_sd = LDprofile$sd, LDprofile_Beta_a = LDprofile$Beta_a, LDprofile_Beta_b = LDprofile$Beta_b, minRandL = 4, minRL = 25, X = NULL),
               "LDprofile_rsq must contain the same number of values as there are bins given in LDprofile_bins")
})

## Test the function with LDprofile_bins as a different length to LDprofile_sd

test_that("Zalpha_all fails when LDprofile_bins and LDprofile_sd are different lengths", {

  expect_error(Zalpha_all(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, LDprofile_sd = c(LDprofile$sd,1), LDprofile_Beta_a = LDprofile$Beta_a, LDprofile_Beta_b = LDprofile$Beta_b, minRandL = 4, minRL = 25, X = NULL),
               "LDprofile_sd must contain the same number of values as there are bins given in LDprofile_bins")
})





## Test the function with LDprofile_Beta_a non-numeric

test_that("Zalpha_all fails when LDprofile_Beta_a is non-numeric", {

  expect_error(Zalpha_all(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, LDprofile_sd = LDprofile$sd, LDprofile_Beta_a = paste0(LDprofile$Beta_a,"r"), LDprofile_Beta_b = LDprofile$Beta_b, minRandL = 4, minRL = 25, X = NULL),
               "LDprofile_Beta_a must be a numeric vector")
})

## Test the function with LDprofile_Beta_b non-numeric

test_that("Zalpha_all fails when LDprofile_Beta_b is non-numeric", {

  expect_error(Zalpha_all(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, LDprofile_sd = LDprofile$sd, LDprofile_Beta_a = LDprofile$Beta_a, LDprofile_Beta_b = paste0(LDprofile$Beta_b,"sd"), minRandL = 4, minRL = 25, X = NULL),
               "LDprofile_Beta_b must be a numeric vector")
})

## Test the function with LDprofile_bins is a different length to LDprofile_Beta_a

test_that("Zalpha_all fails when LDprofile_bins and LDprofile_Beta_a are different lengths", {

  expect_error(Zalpha_all(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, LDprofile_sd = LDprofile$sd, LDprofile_Beta_a = c(LDprofile$Beta_a,1), LDprofile_Beta_b = LDprofile$Beta_b, minRandL = 4, minRL = 25, X = NULL),
               "LDprofile_Beta_a must contain the same number of values as there are bins given in LDprofile_bins")
})

## Test the function with LDprofile_bins as a different length to LDprofile_Beta_b

test_that("Zalpha_all fails when LDprofile_bins and LDprofile_Beta_b are different lengths", {

  expect_error(Zalpha_all(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, LDprofile_sd = LDprofile$sd, LDprofile_Beta_a = LDprofile$Beta_a, LDprofile_Beta_b = c(LDprofile$Beta_b,1), minRandL = 4, minRL = 25, X = NULL),
               "LDprofile_Beta_b must contain the same number of values as there are bins given in LDprofile_bins")
})

## test that only the relevant statistics are calculated with no optional parameters

test_that("Zalpha_all calculates statistics correctly with no optional parameters", {

  expect_equal(Zalpha_all(pos = df$POS, ws  = 3000),
               list(
                 position=c(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500),
                 LR=c(0,13,24,33,40,45,48,49,48,45,40,33,24,13,0),
                 L_plus_R=c(91,78,67,58,51,46,43,42,43,46,51,58,67,78,91)
               ))
})

## test that only the relevant statistics are calculated with only x supplied

test_that("Zalpha_all calculates statistics correctly with only x supplied", {

  expect_equal(Zalpha_all(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7])),
               list(
                 position=c(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500),
                 LR=c(0,13,24,33,40,45,48,49,48,45,40,33,24,13,0),
                 L_plus_R=c(91,78,67,58,51,46,43,42,43,46,51,58,67,78,91),
                 Zalpha=c(NA,NA,NA,NA,((3+1/2)/6+(11+41/144)/45)/2,((6+1/4)/10+(9+41/48)/36)/2,((7+31/72)/15+(7+13/48)/28)/2,((8+17/144)/21+(4+7/16)/21)/2,((9+131/144)/28+(2+13/16)/15)/2,((13+97/144)/36+(1+121/144)/10)/2,((15+25/48)/45+(1+55/144)/6)/2,NA,NA,NA,NA),
                 Zbeta=c(NA,NA,NA,NA,(10+5/18)/40,(10+35/36)/45,(11+103/144)/48,(12+73/144)/49,(11+83/144)/48,(11+17/48)/45,(10+65/144)/40,NA,NA,NA,NA)
               ),tolerance=0.0001)
})

## test that only the relevant statistics are calculated with dist, LDprofile_bins and _rsq supplied

test_that("Zalpha_all calculates statistics correctly with dist, LDprofile_bins and _rsq supplied", {

  expect_equal(Zalpha_all(pos = df$POS, ws  = 3000, dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq),
               list(
                 position=c(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500),
                 LR=c(0,13,24,33,40,45,48,49,48,45,40,33,24,13,0),
                 L_plus_R=c(91,78,67,58,51,46,43,42,43,46,51,58,67,78,91),
                 Zalpha_expected=c(NA,NA,NA,NA,0.390457304338967,0.392014054942343,0.397339546324536,0.398874728980465,0.400715520018796,0.401718327864356,0.399526703832832,NA,NA,NA,NA),
                 Zbeta_expected=c(NA,NA,NA,NA,0.350404001770323,0.357372168444253,0.360647397851440,0.361951221318345,0.362988750761055,0.364752550557575,0.366343120440209,NA,NA,NA,NA)
               ),tolerance=0.0001)
})

## test that only the relevant statistics are calculated with x, dist, LDprofile_bins and _rsq supplied

test_that("Zalpha_all calculates statistics correctly with x, dist, LDprofile_bins and _rsq supplied", {

  expect_equal(Zalpha_all(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq),
               list(
                 position=c(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500),
                 LR=c(0,13,24,33,40,45,48,49,48,45,40,33,24,13,0),
                 L_plus_R=c(91,78,67,58,51,46,43,42,43,46,51,58,67,78,91),
                 Zalpha_expected=c(NA,NA,NA,NA,0.390457304338967,0.392014054942343,0.397339546324536,0.398874728980465,0.400715520018796,0.401718327864356,0.399526703832832,NA,NA,NA,NA),
                 Zbeta_expected=c(NA,NA,NA,NA,0.350404001770323,0.357372168444253,0.360647397851440,0.361951221318345,0.362988750761055,0.364752550557575,0.366343120440209,NA,NA,NA,NA),
                 Zalpha=c(NA,NA,NA,NA,((3+1/2)/6+(11+41/144)/45)/2,((6+1/4)/10+(9+41/48)/36)/2,((7+31/72)/15+(7+13/48)/28)/2,((8+17/144)/21+(4+7/16)/21)/2,((9+131/144)/28+(2+13/16)/15)/2,((13+97/144)/36+(1+121/144)/10)/2,((15+25/48)/45+(1+55/144)/6)/2,NA,NA,NA,NA),
                 Zbeta=c(NA,NA,NA,NA,(10+5/18)/40,(10+35/36)/45,(11+103/144)/48,(12+73/144)/49,(11+83/144)/48,(11+17/48)/45,(10+65/144)/40,NA,NA,NA,NA),
                 Zalpha_rsq_over_expected=c(NA,NA,NA,NA,1.099392777415350,1.175423932431810,0.973208670564282,0.767152145198758,0.691051112331195,0.728230393808371,0.737847444047460,NA,NA,NA,NA),
                 Zalpha_log_rsq_over_expected=c(NA,NA,NA,NA, -0.129685954728232, -0.086676970494806,-0.207059296326176,  -0.317561061974900, -0.355497049708911, -0.345346264723395, -0.340516242835506,NA,NA,NA,NA),
                 Zbeta_rsq_over_expected=c(NA,NA,NA,NA,0.741001103918817,0.693986983022872,0.691113499222754,0.716092979819490,0.663090587779225,0.690647913157800,0.710787633520743,NA,NA,NA,NA),
                 Zbeta_log_rsq_over_expected=c(NA,NA,NA,NA,-0.310756050486805,-0.323333602724241,-0.349639186187507,-0.339617245857445,-0.355799393664026,-0.313469388842524,-0.304918066115218,NA,NA,NA,NA)
                 ),tolerance=0.0001)
})

## test that only the relevant statistics are calculated with only LDprofile_Beta_a and _b not supplied

test_that("Zalpha_all calculates statistics correctly with only LDprofile_Beta_a and _b not supplied", {

  expect_equal(Zalpha_all(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, LDprofile_sd = LDprofile$sd),
               list(
                 position=c(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500),
                 LR=c(0,13,24,33,40,45,48,49,48,45,40,33,24,13,0),
                 L_plus_R=c(91,78,67,58,51,46,43,42,43,46,51,58,67,78,91),
                 Zalpha_expected=c(NA,NA,NA,NA,0.390457304338967,0.392014054942343,0.397339546324536,0.398874728980465,0.400715520018796,0.401718327864356,0.399526703832832,NA,NA,NA,NA),
                 Zbeta_expected=c(NA,NA,NA,NA,0.350404001770323,0.357372168444253,0.360647397851440,0.361951221318345,0.362988750761055,0.364752550557575,0.366343120440209,NA,NA,NA,NA),
                 Zalpha=c(NA,NA,NA,NA,((3+1/2)/6+(11+41/144)/45)/2,((6+1/4)/10+(9+41/48)/36)/2,((7+31/72)/15+(7+13/48)/28)/2,((8+17/144)/21+(4+7/16)/21)/2,((9+131/144)/28+(2+13/16)/15)/2,((13+97/144)/36+(1+121/144)/10)/2,((15+25/48)/45+(1+55/144)/6)/2,NA,NA,NA,NA),
                 Zbeta=c(NA,NA,NA,NA,(10+5/18)/40,(10+35/36)/45,(11+103/144)/48,(12+73/144)/49,(11+83/144)/48,(11+17/48)/45,(10+65/144)/40,NA,NA,NA,NA),
                 Zalpha_rsq_over_expected=c(NA,NA,NA,NA,1.099392777415350,1.175423932431810,0.973208670564282,0.767152145198758,0.691051112331195,0.728230393808371,0.737847444047460,NA,NA,NA,NA),
                 Zalpha_log_rsq_over_expected=c(NA,NA,NA,NA, -0.129685954728232, -0.086676970494806,-0.207059296326176,  -0.317561061974900, -0.355497049708911, -0.345346264723395, -0.340516242835506,NA,NA,NA,NA),
                 Zbeta_rsq_over_expected=c(NA,NA,NA,NA,0.741001103918817,0.693986983022872,0.691113499222754,0.716092979819490,0.663090587779225,0.690647913157800,0.710787633520743,NA,NA,NA,NA),
                 Zbeta_log_rsq_over_expected=c(NA,NA,NA,NA,-0.310756050486805,-0.323333602724241,-0.349639186187507,-0.339617245857445,-0.355799393664026,-0.313469388842524,-0.304918066115218,NA,NA,NA,NA),
                 Zalpha_Zscore=c(NA,NA,NA,NA,0.083669516080429,0.160231324773080,-0.038223960672508,-0.240917156668399,-0.316679396483571,-0.285401890680237,-0.270109843273656,NA,NA,NA,NA),
                 Zbeta_Zscore=c(NA,NA,NA,NA,-0.249622295287423,-0.298443278340203,-0.303226918699104,-0.278319517790552,-0.324158658191404,-0.298557987824230,-0.278980061934724,NA,NA,NA,NA)
               ),tolerance=0.0001)
})

## test that only the relevant statistics are calculated with only LDprofile_sd not supplied

test_that("Zalpha_all calculates statistics correctly with only LDprofile_sd not supplied", {

  expect_equal(Zalpha_all(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, LDprofile_Beta_a = LDprofile$Beta_a, LDprofile_Beta_b = LDprofile$Beta_b),
               list(
                 position=c(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500),
                 LR=c(0,13,24,33,40,45,48,49,48,45,40,33,24,13,0),
                 L_plus_R=c(91,78,67,58,51,46,43,42,43,46,51,58,67,78,91),
                 Zalpha_expected=c(NA,NA,NA,NA,0.390457304338967,0.392014054942343,0.397339546324536,0.398874728980465,0.400715520018796,0.401718327864356,0.399526703832832,NA,NA,NA,NA),
                 Zbeta_expected=c(NA,NA,NA,NA,0.350404001770323,0.357372168444253,0.360647397851440,0.361951221318345,0.362988750761055,0.364752550557575,0.366343120440209,NA,NA,NA,NA),
                 Zalpha=c(NA,NA,NA,NA,((3+1/2)/6+(11+41/144)/45)/2,((6+1/4)/10+(9+41/48)/36)/2,((7+31/72)/15+(7+13/48)/28)/2,((8+17/144)/21+(4+7/16)/21)/2,((9+131/144)/28+(2+13/16)/15)/2,((13+97/144)/36+(1+121/144)/10)/2,((15+25/48)/45+(1+55/144)/6)/2,NA,NA,NA,NA),
                 Zbeta=c(NA,NA,NA,NA,(10+5/18)/40,(10+35/36)/45,(11+103/144)/48,(12+73/144)/49,(11+83/144)/48,(11+17/48)/45,(10+65/144)/40,NA,NA,NA,NA),
                 Zalpha_rsq_over_expected=c(NA,NA,NA,NA,1.099392777415350,1.175423932431810,0.973208670564282,0.767152145198758,0.691051112331195,0.728230393808371,0.737847444047460,NA,NA,NA,NA),
                 Zalpha_log_rsq_over_expected=c(NA,NA,NA,NA, -0.129685954728232, -0.086676970494806,-0.207059296326176,  -0.317561061974900, -0.355497049708911, -0.345346264723395, -0.340516242835506,NA,NA,NA,NA),
                 Zbeta_rsq_over_expected=c(NA,NA,NA,NA,0.741001103918817,0.693986983022872,0.691113499222754,0.716092979819490,0.663090587779225,0.690647913157800,0.710787633520743,NA,NA,NA,NA),
                 Zbeta_log_rsq_over_expected=c(NA,NA,NA,NA,-0.310756050486805,-0.323333602724241,-0.349639186187507,-0.339617245857445,-0.355799393664026,-0.313469388842524,-0.304918066115218,NA,NA,NA,NA),
                 Zalpha_BetaCDF=c(NA,NA,NA,NA,0.498582271149287,0.525043954366501,0.461984646273587,0.409097872871755,0.388665184297966,0.396797658473664,0.393583363498442,NA,NA,NA,NA),
                 Zbeta_BetaCDF=c(NA,NA,NA,NA,0.418016604924725,0.406297914560424,0.401067078642937,0.410829577389033,0.397509236264354,0.405347850057668,0.409781300588220,NA,NA,NA,NA)
                ),tolerance=0.0001)
})

## test that only the relevant statistics are calculated with x not supplied

test_that("Zalpha_all calculates statistics correctly with only x not supplied", {

  expect_equal(Zalpha_all(pos = df$POS, ws  = 3000, dist = df$dist, LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, LDprofile_sd = LDprofile$sd, LDprofile_Beta_a = LDprofile$Beta_a, LDprofile_Beta_b = LDprofile$Beta_b),
               list(
                 position=c(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500),
                 LR=c(0,13,24,33,40,45,48,49,48,45,40,33,24,13,0),
                 L_plus_R=c(91,78,67,58,51,46,43,42,43,46,51,58,67,78,91),
                 Zalpha_expected=c(NA,NA,NA,NA,0.390457304338967,0.392014054942343,0.397339546324536,0.398874728980465,0.400715520018796,0.401718327864356,0.399526703832832,NA,NA,NA,NA),
                 Zbeta_expected=c(NA,NA,NA,NA,0.350404001770323,0.357372168444253,0.360647397851440,0.361951221318345,0.362988750761055,0.364752550557575,0.366343120440209,NA,NA,NA,NA)
               ),tolerance=0.0001)
})

## test that only the relevant statistics are calculated with dist not supplied

test_that("Zalpha_all calculates statistics correctly with dist not supplied", {

  expect_equal(Zalpha_all(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), LDprofile_bins = LDprofile$bin, LDprofile_rsq = LDprofile$rsq, LDprofile_sd = LDprofile$sd, LDprofile_Beta_a = LDprofile$Beta_a, LDprofile_Beta_b = LDprofile$Beta_b),
               list(
                 position=c(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500),
                 LR=c(0,13,24,33,40,45,48,49,48,45,40,33,24,13,0),
                 L_plus_R=c(91,78,67,58,51,46,43,42,43,46,51,58,67,78,91),
                 Zalpha=c(NA,NA,NA,NA,((3+1/2)/6+(11+41/144)/45)/2,((6+1/4)/10+(9+41/48)/36)/2,((7+31/72)/15+(7+13/48)/28)/2,((8+17/144)/21+(4+7/16)/21)/2,((9+131/144)/28+(2+13/16)/15)/2,((13+97/144)/36+(1+121/144)/10)/2,((15+25/48)/45+(1+55/144)/6)/2,NA,NA,NA,NA),
                 Zbeta=c(NA,NA,NA,NA,(10+5/18)/40,(10+35/36)/45,(11+103/144)/48,(12+73/144)/49,(11+83/144)/48,(11+17/48)/45,(10+65/144)/40,NA,NA,NA,NA)
               ),tolerance=0.0001)
})

## test that only the relevant statistics are calculated with LDprofile_bins not supplied

test_that("Zalpha_all calculates statistics correctly with LDprofile_bins not supplied", {

  expect_equal(Zalpha_all(pos = df$POS, ws  = 3000, x = as.matrix(df[,3:7]), dist = df$dist, LDprofile_rsq = LDprofile$rsq, LDprofile_sd = LDprofile$sd, LDprofile_Beta_a = LDprofile$Beta_a, LDprofile_Beta_b = LDprofile$Beta_b),
               list(
                 position=c(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500),
                 LR=c(0,13,24,33,40,45,48,49,48,45,40,33,24,13,0),
                 L_plus_R=c(91,78,67,58,51,46,43,42,43,46,51,58,67,78,91),
                 Zalpha=c(NA,NA,NA,NA,((3+1/2)/6+(11+41/144)/45)/2,((6+1/4)/10+(9+41/48)/36)/2,((7+31/72)/15+(7+13/48)/28)/2,((8+17/144)/21+(4+7/16)/21)/2,((9+131/144)/28+(2+13/16)/15)/2,((13+97/144)/36+(1+121/144)/10)/2,((15+25/48)/45+(1+55/144)/6)/2,NA,NA,NA,NA),
                 Zbeta=c(NA,NA,NA,NA,(10+5/18)/40,(10+35/36)/45,(11+103/144)/48,(12+73/144)/49,(11+83/144)/48,(11+17/48)/45,(10+65/144)/40,NA,NA,NA,NA)
               ),tolerance=0.0001)
})

## test that Zalpha_all works with a missing value
df1<-df
df1$C1[15]<-NA
test_that("Zalpha_all calculates statistics correctly with missing value", {

  expect_equal(Zalpha_all(pos = df1$POS, ws  = 3000, x = as.matrix(df1[,3:7])),
               list(
                 position=c(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500),
                 LR=c(0,13,24,33,40,45,48,49,48,45,40,33,24,13,0),
                 L_plus_R=c(91,78,67,58,51,46,43,42,43,46,51,58,67,78,91),
                 Zalpha=c(NA,NA,NA,NA,0.434953703703704,0.473283179012346,0.397114748677249,0.317791005291005,0.300801917989418,0.322897376543210,0.360532407407407,NA,NA,NA,NA),
                 Zbeta=c(NA,NA,NA,NA,0.248611111111111,0.235185185185185,0.233651620370370,0.257794784580499,0.250144675925926,0.259413580246914,0.271354166666667,NA,NA,NA,NA)
               ),tolerance=0.0001)
})
