# EXAMPLE 1.
#
# Analysis of Underdispersed Word Counts
# 
# Source: http://www.statsci.org/data/oz/wdcount.html
#
# References:
#
#   Modelling_Function_Word_Count.pdf
# 
# Description and Data are available in the website.
#
#Description
#
#In studies aimed at characterising an author's style, 
#samples of n words are taken and the number of function 
#words in each sample counted. Often binomial or Poisson 
#distributions are assumed to hold for the proportions of 
#function words. The table shows the combined frequencies (x) 
#of the articles "the", "a" and "an" in samples from Macauley's 
#"Essay on Milton", taken from the Oxford edition of Macualey's (1923) 
#literary essays. Non-overlapping samples were drawn from 
#opening words of two randomly chosen lines from each of 
#50 pages of printed text, 10 word samples being simply 
#extensions of 5 word samples. 
#
#The data show clear evidence of underdispersion
#
#Data
#
#Occurences      5-Word  10-Word
#0       	45      27
#1     		49      44
#2       	6       26
#3       	0       3
#


# EXAMPLE 2
#
#   References:
#   Pala F. S, Alkaya F., Tabakcioglu K., Tokatli F., Uzal C., Parlar Þ., Algunes C. (2008).
#   The Effects of Micronuclei with Whole Chromosome on Biological Dose Estimation.
#   Turkish Journal of Biology, 32, 283-290.
#
# Pala et.al. (2008), Pg. 286.
#
# Table 2. Distributions and frequencies of Co-60 gamma-ray induced MN of pooled donors’ data.
#
#                         				Micronuclei distribution
# 
# Dose  Cells 	Total 	frequency 	   0   1    2     3     4     5     6    7
#			 MN
#	
#0    20.000 	259 		0.013 	19745 251 	4 	-	- 	- 	- 	-	
#0.1  20.000 	432 		0.021 	19579 410 	11 	-	- 	- 	- 	-
#0.25 20.000 	666 		0.033 	19354 628 	16 	2	- 	- 	- 	-
#0.5 	20.000 	1182 		0.059 	18851 1121 	23 	5 	- 	- 	- 	-
#0.75 20.000 	1660 		0.083 	18402 1549 	36 	13	- 	- 	- 	- 
#1 	20.000 	2376 		0.12 		17769 2119 	84	24 	3 	1 	- 	-
#2 	10.000 	2752 		0.28 		7449 	2373 	89 	48 	8 	5 	- 	- 
#3 	10.000 	4866 		0.49 		5812 	3773 	232	122	44 	15 	2 	- 
#4 	10.000 	7781 		0.78 		3663 	5481 	478 	222 	111 	37 	7 	1 
#5 	10.000 	10950 	1.10 		1652 	6818 	847 	420 	167 	70 	22 	4 
#
#
# In this study,
#
# Separate Co-60 dose response curves for MN in the range 0-5.0 Gy were established using pancentromeric
# FISH probes with lymphocytes from the smoker and non-smoker donors. The data from the 2 donors were pooled and the
# distributions of MN examined is as shown in Table 2
#
# sOME cOMMENTS FROM THIS ARTICLE
# 
# As seen in Table 2, the MN distribution showed
# underdispersion at higher doses. This is a phenomenon
# seen with dicentrics, too. Savage (10) explained this with
# a distortion hypothesis where the total number of
# possible dicentrics is restricted by there being only 46
# chromosomes in each cell. For MN, underdispersion might
# be caused by several micronuclei coalescing into a single
# micronucleus (11).
#  
# where the references:
# (10) 
# Savage JRK. (1970) 
# Sites of radiation induced chromosome exchanges.
# Curr. Top. Radiat. Res. 6:129-194.
#
# (11)
# Littlefield LG, Sayer AM. Frome EL. (1989) 
# Comparison of dose-response parameters for radiation-induced 
# acentric fragments and micronuclei observed in 
# cytokinesis-arrested lymphocytes.
# Mutagenesis 4: 265-270.



source("semi.altham_functions.r")
source("semi.altham_data_rd.R")

#
# Goodness-of-fit test for fitted models.
#

rd=function(x=freq5)
{
# Compute sample dispersion index relative to the binomial 
# distribution of the thea same mean and support.
#
  M=length(x)-1
  Support=0:M
  n=sum(x)               # sample size
  Mean=sum(x*Support)/n
  Var=(sum(x*Support^2)-n*Mean^2)/(n-1)
  
  p.hat=Mean/M
  rd=Var/(Mean*(1-p.hat))

  return(rd)
}



semi.altham.goodness.test=function(ObsFreq, PMF)
{
# Goodness-of-fit test for fitted semi.altham model using MLE.
#
  n=sum(ObsFreq)
  ExpFreq=n*PMF
  Test=sum((ObsFreq-ExpFreq)^2/ExpFreq)
  DF=length(PMF)-1-2
  Pvalue=1-pchisq(Test, DF)

  return(list(testing.statistic=Test, degree.freedom=DF, p.value=Pvalue))
}


semi.altham.example1=function(freq,y,M)
{
testm=matrix(NA,byrow=F,ncol=12)
pm=matrix(NA,byrow=F,ncol=12)
aicm=matrix(NA,byrow=F,ncol=12)
mum=matrix(NA,byrow=T,nrow=12,ncol=1)
sigmam=matrix(NA,byrow=T,nrow=12,ncol=1)
dispersionm=matrix(NA,byrow=T,nrow=12,ncol=1)
rdm=matrix(NA,byrow=T,nrow=12,ncol=1)
beta1m=matrix(NA,byrow=T,nrow=12,ncol=1)
beta2m=matrix(NA,byrow=T,nrow=12,ncol=1)

for (hi in 1:12)
{

mle=mle.optim.semi.altham(hi, y, M)

beta1=mle$MLE[1]; beta2=mle$MLE[2]
  beta1m[hi]=beta1
  beta2m[hi]=beta2

beta=c(beta1, beta2)
PMF=dsemi.altham(hi, M, beta)

moment=msemi.altham(hi, M, beta)
  mum[hi]=moment[1]
  sigmam[hi]=moment[4]
  dispersionm[hi]=moment[5]
  rdm[hi]=rd(freq)
newfreq=rep(0, 4)
newfreq[1]=freq[1]
newfreq[2]=freq[2]
newfreq[3]=freq[3]
newfreq[4]=sum(freq[-(1:3)])

newPMF=rep(0, 4)
newPMF[1]=PMF[1]
newPMF[2]=PMF[2]
newPMF[3]=PMF[3]
newPMF[4]=1-sum(newPMF[1:3])


newgoodness=semi.altham.goodness.test(newfreq, newPMF)

Test=newgoodness$testing.statistic
testm[hi]=Test
Pvalue=newgoodness$p.value
pm[hi]=Pvalue
AIC=mle$AIC
aicm[hi]=AIC

}

list(AIC=aicm, Test.statistics=testm, p.value=pm,
     mu=mum, sigma2=sigmam,Dispersion=dispersionm, RD=rdm,
     beta1=beta1m, beta2=beta2m)

}


##################################################
# PART 1: Data from Macaulay’s work


freq5=c(45, 49,  6,  0,  0,  0)  # 5-word samples from Macaulay’s work

freq10=c(27, 44, 26,  3,  0,  0,  0,  0,  0,  0,  0) # 10-word samples from Macaulay’s work


rd5=rd(freq5)
rd5
#[1] 0.6749975
rd10=rd(freq10)
rd10
#[1] 0.6959728


# Fitting of semi.altham 
#

x5=NULL
M=length(freq5)-1
for(i in 1:length(freq5))
{
  if(freq5[i]!=0) x5=c(x5, rep(i-1, freq5[i]))
}

x5
sum(x5)
length(x5)==sum(freq5)

x10=NULL
M=length(freq10)-1
for(i in 1:length(freq10))
{
  if(freq10[i]!=0) x10=c(x10, rep(i-1, freq10[i]))
}

x10
sum(x10)
length(x10)==sum(freq10)



# Test
#####################################
# 5-word samples from Macaulay’s work


freq=freq5; y=x5; M=5
out=semi.altham.example1(freq,y,M)

#> out
#$AIC
#       [,1]     [,2]     [,3]     [,4]     [,5]     [,6]     [,7]     [,8]
#[1,] 179.69 179.6875 179.7261 179.6604 179.7168 179.6669 179.6866 179.6935
#         [,9]   [,10]    [,11]    [,12]
#[1,] 179.9007 179.598 179.7265 179.6601
#
#$Test.statistics
#           [,1]       [,2]      [,3]       [,4]      [,5]       [,6]       [,7]
#[1,] 0.08188775 0.08088954 0.1026752 0.06591256 0.0974091 0.06949528 0.08036456
#           [,8]      [,9]      [,10]     [,11]      [,12]
#[1,] 0.08425458 0.2072750 0.03236179 0.1029450 0.06574598
#
#$p.value
#          [,1]      [,2]      [,3]      [,4]      [,5]     [,6]     [,7]
#[1,] 0.7747553 0.7760955 0.7486426 0.7973838 0.7549612 0.792073 0.776804
#          [,8]      [,9]    [,10]     [,11]     [,12]
#[1,] 0.7716124 0.6489117 0.857236 0.7483237 0.7976344
#
#$mu
#           [,1]
# [1,] 0.6100196
# [2,] 0.6100031
# [3,] 0.6100001
# [4,] 0.6100000
# [5,] 0.6099999
# [6,] 0.6100000
# [7,] 0.6099999
# [8,] 0.6099999
# [9,] 0.6099999
#[10,] 0.6100002
#[11,] 0.6100001
#[12,] 0.6100020
#
#$sigma2
#           [,1]
# [1,] 0.3580905
# [2,] 0.3578980
# [3,] 0.3578997
# [4,] 0.3578999
# [5,] 0.3578996
# [6,] 0.3578997
# [7,] 0.3578996
# [8,] 0.3578997
# [9,] 0.3578997
#[10,] 0.3578998
#[11,] 0.3578997
#[12,] 0.3578985

#$Dispersion
#           [,1]
# [1,] 0.6685847
# [2,] 0.6682408
# [3,] 0.6682469
# [4,] 0.6682472
# [5,] 0.6682468
# [6,] 0.6682469
# [7,] 0.6682469
# [8,] 0.6682471
# [9,] 0.6682470
#[10,] 0.6682469
#[11,] 0.6682468
#[12,] 0.6682429


#$RD
#           [,1]
# [1,] 0.6749975
# [2,] 0.6749975
# [3,] 0.6749975
# [4,] 0.6749975
# [5,] 0.6749975
# [6,] 0.6749975
# [7,] 0.6749975
# [8,] 0.6749975
# [9,] 0.6749975
#[10,] 0.6749975
#[11,] 0.6749975
#[12,] 0.6749975


#$beta1
#            [,1]
# [1,] -1.2060777
# [2,] -1.2639003
# [3,]  0.8532594
# [4,] -3.2698059
# [5,] -0.3756598
# [6,] -2.0404310
# [7,] -1.3689970
# [8,] -1.0464641
# [9,]  1.0011088
#[10,] -3.4362201
#[11,] -0.5318735
#[12,] -1.8847263

#$beta2
#           [,1]
# [1,] 1.1153263
# [2,] 1.1051931
# [3,] 0.6636073
# [4,] 1.5702289
# [5,] 0.9765712
# [6,] 1.2568985
# [7,] 1.0955567
# [8,] 1.1374076
# [9,] 0.5097133
#[10,] 1.7390593
#[11,] 0.9658730
#[12,] 1.2679990



#######################################
# 10-word samples from Macaulay’s work


freq=freq10; y=x10; M=10
out=semi.altham.example1(freq,y,M)

out
#$AIC
#         [,1]     [,2]     [,3]     [,4]     [,5]     [,6]     [,7]     [,8]
#[1,] 238.7758 238.7733 239.4259 238.3622 239.1435 238.4964 238.7706 238.7812
#        [,9]    [,10]    [,11]   [,12]
#[1,] 245.296 238.8775 239.2745 238.428
#
#$Test.statistics
#         [,1]      [,2]     [,3]      [,4]     [,5]      [,6]      [,7]
#[1,] 0.488678 0.4869394 1.029879 0.1549168 0.795584 0.2599626 0.4847583
#         [,8]     [,9]     [,10]     [,11]     [,12]
#[1,] 0.492587 6.442349 0.8428662 0.9111392 0.2019770
#
#$p.value
#          [,1]      [,2]      [,3]      [,4]     [,5]      [,6]      [,7]
#[1,] 0.4845176 0.4852958 0.3101871 0.6938803 0.372417 0.6101458 0.4862748
#          [,8]       [,9]     [,10]     [,11]     [,12]
#[1,] 0.4827756 0.01114312 0.3585783 0.3398123 0.6531298
#
#$mu
#          [,1]
# [1,] 1.050001
# [2,] 1.050034
# [3,] 1.049905
# [4,] 1.050000
# [5,] 1.050002
# [6,] 1.050000
# [7,] 1.050000
# [8,] 1.050005
# [9,] 1.050004
#[10,] 1.050000
#[11,] 1.050000
#[12,] 1.050001

#$sigma2
#           [,1]
# [1,] 0.6474970
# [2,] 0.6474193
# [3,] 0.6476249
# [4,] 0.6474980
# [5,] 0.6474982
# [6,] 0.6474977
# [7,] 0.6475000
# [8,] 0.6475076
# [9,] 0.6474964
#[10,] 0.6474986
#[11,] 0.6474966
#[12,] 0.6474979

#$Dispersion
#           [,1]
# [1,] 0.6890093
# [2,] 0.6889073
# [3,] 0.6892010
# [4,] 0.6890108
# [5,] 0.6890100
# [6,] 0.6890103
# [7,] 0.6890129
# [8,] 0.6890181
# [9,] 0.6890070
#[10,] 0.6890114
#[11,] 0.6890094
#[12,] 0.6890100

#$beta1
#            [,1]
# [1,] -1.2050598
# [2,] -1.2314711
# [3,]  1.4005659
# [4,] -3.8166836
# [5,] -0.4265813
# [6,] -1.9859198
# [7,] -1.2944406
# [8,] -1.1156390
# [9,]  1.3959013
#[10,] -3.9624624
#[11,] -0.5997132
#[12,] -1.8141261

#$beta2
#           [,1]
# [1,] 0.6262076
# [2,] 0.6241141
# [3,] 0.2821010
# [4,] 0.9738040
# [5,] 0.5158499
# [6,] 0.7378968
# [7,] 0.6206741
# [8,] 0.6317225
# [9,] 0.0750704
#[10,] 1.2658395
#[11,] 0.5145035
#[12,] 0.7400238



##########################################
# PART 2 - Data from Chesterton’s work



Freq5=c(32, 35, 3, 0, 0, 0) # 5-word samples from Chesterton’s work

Freq10=c(14, 38, 16, 2, 0, 0, 0, 0, 0, 0, 0) # 10-word samples from Chesterton’s work


Rd5=rd(Freq5)
Rd5
#[1] 0.6442177
Rd10=rd(Freq10)
Rd10
#0.5613253

# Fitting of semi.altham hi=1 
#

X5=NULL
M=length(Freq5)-1
for(i in 1:length(Freq5))
{
  if(Freq5[i]!=0) X5=c(X5, rep(i-1, Freq5[i]))
}

X5
sum(X5)
length(X5)==sum(Freq5)

X10=NULL
M=length(Freq10)-1
for(i in 1:length(Freq10))
{
  if(Freq10[i]!=0) X10=c(X10, rep(i-1, Freq10[i]))
}

X10
sum(X10)
length(X10)==sum(Freq10)


#######################################
# 5-word samples from Chesterton’s work

freq=Freq5; y=x5; M=5
out=semi.altham.example1(freq,y,M)
out

#$AIC
#       [,1]     [,2]     [,3]     [,4]     [,5]     [,6]     [,7]     [,8]
#[1,] 179.69 179.6875 179.7261 179.6604 179.7168 179.6669 179.6866 179.6935
#         [,9]   [,10]    [,11]    [,12]
#[1,] 179.9007 179.598 179.7265 179.6601
#
#$Test.statistics
#          [,1]      [,2]      [,3]      [,4]      [,5]      [,6]      [,7]
#[1,] 0.3432585 0.3404075 0.3356536 0.3440861 0.3367377 0.3431778 0.3405436
#          [,8]      [,9]     [,10]     [,11]     [,12]
#[1,] 0.3396346 0.3199882 0.3536548 0.3355903 0.3441177

#$p.value
#          [,1]      [,2]      [,3]      [,4]      [,5]      [,6]      [,7]
#[1,] 0.5579544 0.5595941 0.5623488 0.5574801 0.5617184 0.5580006 0.5595156
#          [,8]      [,9]     [,10]     [,11]    [,12]
#[1,] 0.5600402 0.5716147 0.5520515 0.5623857 0.557462

#$mu
#           [,1]
# [1,] 0.6100196
# [2,] 0.6100031
# [3,] 0.6100001
# [4,] 0.6100000
# [5,] 0.6099999
# [6,] 0.6100000
# [7,] 0.6099999
# [8,] 0.6099999
# [9,] 0.6099999
#[10,] 0.6100002
#[11,] 0.6100001
#[12,] 0.6100020

#$sigma2
#           [,1]
# [1,] 0.3580905
# [2,] 0.3578980
# [3,] 0.3578997
# [4,] 0.3578999
# [5,] 0.3578996
# [6,] 0.3578997
# [7,] 0.3578996
# [8,] 0.3578997
# [9,] 0.3578997
#[10,] 0.3578998
#[11,] 0.3578997
#[12,] 0.3578985

#$Dispersion
#           [,1]
# [1,] 0.6685847
# [2,] 0.6682408
# [3,] 0.6682469
# [4,] 0.6682472
# [5,] 0.6682468
# [6,] 0.6682469
# [7,] 0.6682469
# [8,] 0.6682471
# [9,] 0.6682470
#[10,] 0.6682469
#[11,] 0.6682468
#[12,] 0.6682429

#$beta1
#            [,1]
# [1,] -1.2060777
# [2,] -1.2639003
# [3,]  0.8532594
# [4,] -3.2698059
# [5,] -0.3756598
# [6,] -2.0404310
# [7,] -1.3689970
# [8,] -1.0464641
# [9,]  1.0011088
#[10,] -3.4362201
#[11,] -0.5318735
#[12,] -1.8847263

#$beta2
#           [,1]
# [1,] 1.1153263
# [2,] 1.1051931
# [3,] 0.6636073
# [4,] 1.5702289
# [5,] 0.9765712
# [6,] 1.2568985
# [7,] 1.0955567
# [8,] 1.1374076
# [9,] 0.5097133
#[10,] 1.7390593
#[11,] 0.9658730
#[12,] 1.2679990



#########################################
# 10 -word samples from Chesterton’s work

freq=Freq10; y=x10; M=10
out=semi.altham.example1(freq,y,M)

#> out
#$AIC
#         [,1]     [,2]     [,3]     [,4]     [,5]     [,6]     [,7]     [,8]
#[1,] 238.7758 238.7733 239.4259 238.3622 239.1435 238.4964 238.7706 238.7812
#        [,9]    [,10]    [,11]   [,12]
#[1,] 245.296 238.8775 239.2745 238.428
#
#$Test.statistics
#         [,1]     [,2]     [,3]     [,4]     [,5]     [,6]     [,7]     [,8]
#[1,] 2.026831 2.026324 1.693619 2.457102 1.809507 2.280119 2.030739 2.023209
#         [,9]    [,10]    [,11]    [,12]
#[1,] 1.760522 5.418062 1.751944 2.367064
#
#$p.value
#          [,1]      [,2]      [,3]      [,4]      [,5]      [,6]      [,7]
#[1,] 0.1545425 0.1545941 0.1931246 0.1169944 0.1785674 0.1310418 0.1541456
#          [,8]      [,9]      [,10]     [,11]     [,12]
#[1,] 0.1549114 0.1845595 0.01992947 0.1856326 0.1239199
#
#$mu
#          [,1]
# [1,] 1.050001
# [2,] 1.050034
# [3,] 1.049905
# [4,] 1.050000
# [5,] 1.050002
# [6,] 1.050000
# [7,] 1.050000
# [8,] 1.050005
# [9,] 1.050004
#[10,] 1.050000
#[11,] 1.050000
#[12,] 1.050001

#$sigma2
#           [,1]
# [1,] 0.6474970
# [2,] 0.6474193
# [3,] 0.6476249
# [4,] 0.6474980
# [5,] 0.6474982
# [6,] 0.6474977
# [7,] 0.6475000
# [8,] 0.6475076
# [9,] 0.6474964
#[10,] 0.6474986
#[11,] 0.6474966
#[12,] 0.6474979

#$Dispersion
#           [,1]
# [1,] 0.6890093
# [2,] 0.6889073
# [3,] 0.6892010
# [4,] 0.6890108
# [5,] 0.6890100
# [6,] 0.6890103
# [7,] 0.6890129
# [8,] 0.6890181
# [9,] 0.6890070
#[10,] 0.6890114
#[11,] 0.6890094
#[12,] 0.6890100

#$beta1
#            [,1]
# [1,] -1.2050598
# [2,] -1.2314711
# [3,]  1.4005659
# [4,] -3.8166836
# [5,] -0.4265813
# [6,] -1.9859198
# [7,] -1.2944406
# [8,] -1.1156390
# [9,]  1.3959013
#[10,] -3.9624624
#[11,] -0.5997132
#[12,] -1.8141261

#$beta2
#           [,1]
# [1,] 0.6262076
# [2,] 0.6241141
# [3,] 0.2821010
# [4,] 0.9738040
# [5,] 0.5158499
# [6,] 0.7378968
# [7,] 0.6206741
# [8,] 0.6317225
# [9,] 0.0750704
#[10,] 1.2658395
#[11,] 0.5145035
#[12,] 0.7400238


###############################################################################
#
# Analysis of examples in
#
# Chakraborty, S. and Das, K.K. (2006).
# On some properties of a class of weighted quasi-binomial distributions.
# Journal of Statistical Planning and Inference,
# 136, 159--182.
# 
# File: Weighted-Quasi-binomial.pdf
#
# Remark:
#   All data examples are not exact binomial data, but count data.
#   They are just fitted by QBD-I and QBD-II which are constructed
#   for binomial data which have a fixed upper bound.
#   It is ok to fit the truncated count data by the model for
#   binomial data.
#
# Data extracted from the reference.
#

# Observed and expected frequencies of European Corn borer in 1296 Corn
# plants.
# Author's chi-squares: 0.0758 (df=1) for QBD-I, 0.0677 (df=1) for QBD-II.
#
Ex1=matrix(0, nrow=5, ncol=3)
rownames(Ex1)=c(as.character(0:3), "4+")
colnames(Ex1)=c("ObsNo", "QBD-I", "QBD-II")
Ex1[, 1]=c(907, 275, 88, 23, 3)
Ex1[, 2]=c(906.41, 277.40, 85.90, 22.59, 3.70)
Ex1[, 3]=c(906.45, 277.26, 86.01, 22.60, 3.68)


# Distribution of yeast cells per square in a haemacytometer.
# Author's chi-squares: 3.6087 (df=1) for QBD-I, 3.6082 (df=1) for QBD-II.
#
Ex2=matrix(0, nrow=6, ncol=3)
rownames(Ex2)=as.character(0:5)
colnames(Ex2)=c("ObsNo", "QBD-I", "QBD-II")
Ex2[, 1]=c(213, 128, 37, 18, 3, 1)
Ex2[, 2]=c(215.73, 118.28, 47.23, 14.89, 3.43, 0.44)
Ex2[, 3]=c(215.72, 118.29, 47.25, 14.89, 3.42, 0.44)


# Distribution of number of seeds by time of day.
# Author's chi-squares: 0.4575 (df=2) for QBD-I, 0.6225 (df=2) for QBD-II.
#
Ex3=matrix(0, nrow=6, ncol=3)
rownames(Ex3)=as.character(0:5)
colnames(Ex3)=c("ObsNo", "QBD-I", "QBD-II")
Ex3[, 1]=c(7, 4, 5, 5, 4, 7)
Ex3[, 2]=c(6.50, 5.38, 4.55, 4.18, 4.40, 6.99)
Ex3[, 3]=c(6.77, 4.95, 4.28, 4.28, 4.95, 6.77)


# Distribution of number of hits per square.
# Author's chi-squares: 0.9409 (df=2) for QBD-I, 0.9344 (df=2) for QBD-II.
#
Ex4=matrix(0, nrow=6, ncol=3)
rownames(Ex4)=as.character(0:5)
colnames(Ex4)=c("ObsNo", "QBD-I", "QBD-II")
Ex4[, 1]=c(229, 211, 93, 35, 7, 1)
Ex4[, 2]=c(231.35, 203.60, 100.24, 33.01, 7.04, 0.76)
Ex4[, 3]=c(231.35, 203.59, 100.26, 33.01, 7.04, 0.76)


# Functions for dispersion investigation
#

rd=function(x=freq5)
{
# Compute sample dispersion index relative to the binomial 
# distribution of the thea same mean and support.
#
  M=length(x)-1
  Support=0:M
  n=sum(x)               # sample size
  Mean=sum(x*Support)/n
  Var=(sum(x*Support^2)-n*Mean^2)/(n-1)
  
  p.hat=Mean/M
  rd=Var/(Mean*(1-p.hat))

  return(rd)
}

d1=rd(Ex1[, 1])
d1
#[1] 1.391126

d2=rd(Ex2[, 1])
d2
#[1] 1.380748

d3=rd(Ex3[, 1])
d3
#[1] 2.787097

d4=rd(Ex4[, 1])
d4
#[1] 1.237367



Fit=function(hi,freq)
{
  x=NULL
  M=length(freq)-1
  for(i in 1:(M+1))
  {
    if(freq[i]!=0) x=c(x, rep(i-1, freq[i]))
  }

  y=x;
  mle.tmp=mle.optim.semi.altham(hi, y, M)
  beta1=mle.tmp$MLE[1]; beta2=mle.tmp$MLE[2]
  beta=c(beta1, beta2)
  PMF=dsemi.altham(hi, M, beta)
  AIC=mle.tmp$AIC

  list(PMF=PMF, AIC=AIC)
}

#Test

fit1=Fit(hi=1,Ex1[, 1])
fit2=Fit(hi=1,Ex2[, 1])
fit3=Fit(hi=1,Ex3[, 1])
fit4=Fit(hi=1,Ex4[, 1])

#
# Goodness-of-fit test for fitted models.
#

chi2=function(Obsd, Exptd)
{
  sum((Obsd-Exptd)^2/Exptd)
}

Pmf=function(hi, Freq)
{
  M=length(Freq)-1; 
  fit=Fit(hi,Freq)
  PMF=fit$PMF
  PMF
}

# Test
hi=1;
Exptd1=sum(hi, Ex1[, 1])*Pmf(hi, Ex1[, 1])
Exptd2=sum(hi, Ex2[, 1])*Pmf(hi, Ex2[, 1])
Exptd3=sum(hi, Ex3[, 1])*Pmf(hi, Ex3[, 1])
Exptd4=sum(hi, Ex4[, 1])*Pmf(hi, Ex4[, 1])

# Example 1: data grouping--0, 1, 2, 3-4+ ==> df=1
# For QBD-I
Obsd=c(Ex1[1:3, 1], sum(Ex1[4:5, 1]))
Exptd=c(Ex1[1:3, 2], sum(Ex1[4:5, 2]))
chi2(Obsd, Exptd)
#[1] 0.07568598
# For QBD-II
Obsd=c(Ex1[1:3, 1], sum(Ex1[4:5, 1]))
Exptd=c(Ex1[1:3, 3], sum(Ex1[4:5, 3]))
chi2(Obsd, Exptd)
#[1] 0.067781
# For semi.altham
Obsd=c(Ex1[1:3, 1], sum(Ex1[4:5, 1]))
Exptd=c(Exptd1[1:3], sum(Exptd1[4:5]))
chi2(Obsd, Exptd)
#[1] 0.8834984


# Example 2: data grouping--0, 1, 2, 3-5 ==>  df=1
# For QBD-I
Obsd=c(Ex2[1:3, 1], sum(Ex2[4:6, 1]))
Exptd=c(Ex2[1:3, 2], sum(Ex2[4:6, 2]))
chi2(Obsd, Exptd)
#[1] 3.608704
# For QBD-II
Obsd=c(Ex2[1:3, 1], sum(Ex2[4:6, 1]))
Exptd=c(Ex2[1:3, 3], sum(Ex2[4:6, 3]))
chi2(Obsd, Exptd)
#[1] 3.618234
# For semi.altham
Obsd=c(Ex2[1:3, 1], sum(Ex2[4:6, 1]))
Exptd=c(Exptd2[1:3], sum(Exptd2[4:6]))
chi2(Obsd, Exptd)
#[1] 4.063314


# Example 3: data grouping--0, 1, 2, 3-4, 5,  ==>  df=2
# For QBD-I: 
Obsd=c(Ex3[1:3, 1], sum(Ex3[4:5, 1]), Ex3[6, 1])
Exptd=c(Ex3[1:3, 2], sum(Ex3[4:5, 2]), Ex3[6, 2])
chi2(Obsd, Exptd)
#[1] 0.4575185
# For QBD-II:
Obsd=c(Ex3[1:3, 1], sum(Ex3[4:5, 1]), Ex3[6, 1])
Exptd=c(Ex3[1:3, 3], sum(Ex3[4:5, 3]), Ex3[6, 3])
chi2(Obsd, Exptd)
#[1] 0.3248038      # Authors made mistake here: 0.6225
# For semi.altham
Obsd=c(Ex3[1:3, 1], sum(Ex3[4:5, 1]), Ex3[6, 1])
Exptd=c(Exptd3[1:3], sum(Exptd3[4:5]), Exptd3[6])
chi2(Obsd, Exptd)
#[1] 0.3685891


# Example 4: data grouping--0, 1, 2, 3, 4-5,  ==>  df=2
# For QBD-I: 
Obsd=c(Ex4[1:4, 1], sum(Ex4[5:6, 1]))
Exptd=c(Ex4[1:4, 2], sum(Ex4[5:6, 2]))
chi2(Obsd, Exptd)
#[1] 0.9408454
# For QBD-II:
Obsd=c(Ex4[1:4, 1], sum(Ex4[5:6, 1]))
Exptd=c(Ex4[1:4, 3], sum(Ex4[5:6, 3]))
chi2(Obsd, Exptd)
#[1] 0.9443742
# For EDU
Obsd=c(Ex4[1:4, 1], sum(Ex4[5:6, 1]))
Exptd=c(Exptd4[1:4], sum(Exptd4[5:6]))
chi2(Obsd, Exptd)
#[1] 2.118788




Chakraborty.examples=function(ex.no,Ex)
{
## matrix for chi-square test-statistics
#
# testm[i]
# 
# the i.th row is for hi=i

testm=matrix(NA,byrow=T,nrow=12,ncol=1) 

## matrix for all p-values
#
# pm[i]
# 
## the i.th row is for hi=i


pm=matrix(NA,byrow=T,nrow=12,ncol=1)

## matrix for all aic values of the model
#
## aicm[i]
# 
## the i.th value is for hi=i

aicm=matrix(NA,byrow=T,nrow=12,ncol=1)


mum=matrix(NA,byrow=T,nrow=12,ncol=1)
sigmam=matrix(NA,byrow=T,nrow=12,ncol=1)
dispersionm=matrix(NA,byrow=T,nrow=12,ncol=1)

 if(ex.no==1 || ex.no==2) df=1;
 if(ex.no==3 || ex.no==4) df=2;

 #df=length(PMF)-1-2

Ex=Ex[,1]

for (hi in 1:12)
 {
  z=NULL; 
  M=length(Ex)-1
   for (j in 1:length(Ex))
   {
    if(Ex[j]!=0) z=c(z,rep(j-1,Ex[j]))
   }

  y=z; # count data
  mle=mle.optim.semi.altham(hi, y, M)
  aicm[hi,1]=mle$AIC
  beta1=mle$MLE[1]; beta2=mle$MLE[2]
  beta=c(beta1, beta2)
  PMF=dsemi.altham(hi, M, beta)
  moment=msemi.altham(hi, M, beta)
  mum[hi,1]=moment[1]
  sigmam[hi,1]=moment[4]
  dispersionm[hi,1]=moment[5]
  Exptd=sum(Ex)*Pmf(hi, Ex)

  if (ex.no==1) 
  {
   Obsd=c(Ex[1:3], sum(Ex[4:5]))
   Exptd=c(Exptd[1:3], sum(Exptd[4:5]))
  }
  if (ex.no==2) 
  {
   Obsd=c(Ex[1:3], sum(Ex[4:6]))
   Exptd=c(Exptd[1:3], sum(Exptd[4:6]))
   }
   if (ex.no==3) 
   {
    Obsd=c(Ex[1:3], sum(Ex[4:5]), Ex[6])
    Exptd=c(Exptd[1:3], sum(Exptd[4:5]), Exptd[6])
    
   }

  if (ex.no==4) 
  {
   Obsd=c(Ex[1:4], sum(Ex[5:6]))
   Exptd=c(Exptd[1:4], sum(Exptd[5:6]))
  }

    testm[hi,1]=chi2(Obsd, Exptd)
    pm[hi,1]=1-pchisq(testm[hi,1],df)  
 }



list(AIC=aicm, Test.chi=testm, P.value=pm, mu=mum, sigma2=sigmam, Dispersion=dispersionm)


}

# Test

out1=Chakraborty.examples(1,Ex1)
out2=Chakraborty.examples(2,Ex2)
out3=Chakraborty.examples(3,Ex3)
out4=Chakraborty.examples(4,Ex4)

Chakraborty.AIC=matrix(c(out1$AIC,out2$AIC,out3$AIC, out4$AIC), byrow=F,ncol=4)
Chakraborty.chi.test=matrix(c(out1$Test.chi,out2$Test.chi,out3$Test.chi, out4$Test.chi), byrow=F,ncol=4)
Chakraborty.p.value=matrix(c(out1$P.value,out2$P.value,out3$P.value, out4$P.value), byrow=F,ncol=4)
Chakraborty.mu=matrix(c(out1$mu,out2$mu,out3$mu, out4$mu), byrow=F,ncol=4)
Chakraborty.sigma=matrix(c(out1$sigma2,out2$sigma2,out3$sigma2, out4$sigma2), byrow=F,ncol=4)
Chakraborty.dispersion=matrix(c(out1$Dispersion,out2$Dispersion,out3$Dispersion, out4$Dispersion), byrow=F,ncol=4)
Chakraborty.AIC
Chakraborty.chi.test
Chakraborty.p.value
Chakraborty.mu
Chakraborty.sigma
Chakraborty.dispersion

#> Chakraborty.AIC
#          [,1]     [,2]     [,3]     [,4]
# [1,] 2201.245 897.9254 117.6120 1462.090
# [2,] 2200.200 898.3593 117.8345 1462.541
# [3,] 2202.126 896.5845 117.9448 1460.697
# [4,] 2200.690 899.8517 117.3484 1464.565
# [5,] 2203.219 896.7713 117.7708 1460.932
# [6,] 2200.013 899.6723 117.5748 1464.258
# [7,] 2200.309 898.3731 117.7708 1462.586
# [8,] 2202.526 897.6044 117.5748 1461.782
# [9,] 2207.112 896.2945 119.4851 1463.733
#[10,] 2201.169 909.8093 116.9552 1479.432
#[11,] 2203.892 896.5044 117.8345 1460.739
#[12,] 2199.815 900.2706 117.5544 1465.031

#> Chakraborty.chi.test
#            [,1]       [,2]        [,3]       [,4]
# [1,] 0.88340849  4.0709726 0.348857532  2.1207300
# [2,] 0.47138727  4.3125777 0.444382839  2.4330207
# [3,] 1.76614471  2.5235082 0.524337239  0.7481686
# [4,] 0.34293260  6.0471189 0.210029191  4.5916202
# [5,] 2.12302650  2.8038856 0.493675842  0.9950141
# [6,] 0.19981505  5.6710453 0.315787286  4.1748225
# [7,] 0.48709864  4.3679267 0.416088567  2.5023598
# [8,] 1.43546373  3.8007586 0.387100383  1.8306223
# [9,] 6.96372289  0.5060037 1.345994715  3.6989823
#[10,] 0.91930804 15.1880084 0.001628746 20.2011257
#[11,] 2.61029934  2.4775928 0.545101466  0.7993686
#[12,] 0.09885083  6.2147511 0.302599636  4.9262636

#> Chakraborty.p.value
#             [,1]         [,2]      [,3]         [,4]
# [1,] 0.347269817 4.362615e-02 0.8399367 3.463294e-01
# [2,] 0.492349650 3.783160e-02 0.8007621 2.962622e-01
# [3,] 0.183859922 1.121608e-01 0.7693813 6.879189e-01
# [4,] 0.558141343 1.392901e-02 0.9003114 1.006798e-01
# [5,] 0.145099567 9.403617e-02 0.7812673 6.080446e-01
# [6,] 0.654870174 1.724718e-02 0.8539406 1.240077e-01
# [7,] 0.485224399 3.662152e-02 0.8121711 2.861669e-01
# [8,] 0.230874820 5.122937e-02 0.8240285 4.003920e-01
# [9,] 0.008317878 4.768740e-01 0.5101771 1.573172e-01
#[10,] 0.337656724 9.731947e-05 0.9991860 4.105644e-05
#[11,] 0.106171720 1.154789e-01 0.7614348 6.705317e-01
#[12,] 0.753213066 1.266902e-02 0.8595899 8.516780e-02

#> Chakraborty.mu
#
#           [,1]      [,2]     [,3]      [,4]
# [1,] 0.4104970 0.6825019 2.500004 0.9288247
# [2,] 0.4104953 0.6824999 2.499895 0.9288193
# [3,] 0.4104958 0.6825028 2.500000 0.9284720
# [4,] 0.4104611 0.6825018 2.500011 0.9289061
# [5,] 0.4105028 0.6825007 2.499960 0.9288228
# [6,] 0.4104946 0.6825018 2.500068 0.9288234
# [7,] 0.4104937 0.6825002 2.500000 0.9288229
# [8,] 0.4106167 0.6825060 2.500002 0.9289031
# [9,] 0.4104952 0.6825017 2.500000 0.9286678
#[10,] 0.4104987 0.6825021 2.499978 0.9288543
#[11,] 0.4104987 0.6825108 2.499949 0.9288212
#[12,] 0.4104948 0.6825007 2.500024 0.9289485

#> Chakraborty.sigma
#           [,1]      [,2]     [,3]      [,4]
# [1,] 0.5120548 0.8116896 3.374963 0.9341654
# [2,] 0.5120494 0.8116837 3.375446 0.9341586
# [3,] 0.5120469 0.8117147 3.374920 0.9344073
# [4,] 0.5120137 0.8116818 3.374921 0.9340742
# [5,] 0.5120510 0.8116842 3.375119 0.9341620
# [6,] 0.5120493 0.8116876 3.374525 0.9341645
# [7,] 0.5120463 0.8116839 3.374923 0.9341637
# [8,] 0.5119518 0.8116932 3.374925 0.9340793
# [9,] 0.5120488 0.8116976 3.374925 0.9339212
#[10,] 0.5120488 0.8116894 3.375069 0.9341205
#[11,] 0.5120484 0.8116628 3.375190 0.9341592
#[12,] 0.5120485 0.8116860 3.374912 0.9342660

#> Chakraborty.dispersion
#          [,1]     [,2]     [,3]     [,4]
# [1,] 1.390056 1.377286 2.699970 1.235208
# [2,] 1.390046 1.377279 2.700357 1.235205
# [3,] 1.390038 1.377327 2.699936 1.235891
# [4,] 1.390051 1.377272 2.699937 1.235004
# [5,] 1.390028 1.377278 2.700095 1.235206
# [6,] 1.390048 1.377282 2.699620 1.235209
# [7,] 1.390042 1.377279 2.699938 1.235208
# [8,] 1.389417 1.377285 2.699940 1.235014
# [9,] 1.390045 1.377299 2.699940 1.235047
#[10,] 1.390034 1.377285 2.700055 1.235119
#[11,] 1.390033 1.377225 2.700152 1.235204
#[12,] 1.390045 1.377281 2.699930 1.235214


########################################################
########################################################
# 1.
#   Zelterman, Daniel. (2004). 
#   Discrete Distributions - Applications in the Health Sciences.
#   John Wiley & Sons, West Sussex.
# Pg. 145
#
# 

example.zelterman=function(tab)
{

dim.tab=dim(tab)

freq.tab=tab[4:dim.tab[1], 4:dim.tab[2]]
dim.freq=dim(freq.tab)


## matrix to record chi-square test-statistics
#
#testm[i,j]
# 
# the i.th row is for hi=i
# and j.th colomn is for the j.th row of the data
#

testm=matrix(NA,byrow=T,nrow=12,ncol=dim.freq[1]) 

## matrix for all p-values pm[i,j]
# 
## the i.th row is for hi=i
## and j.th colomn is for the j.th row of the data
##
#
pm=matrix(NA,byrow=T,nrow=12,ncol=dim.freq[1])
#
## matrix for all aic values of the model
#
## aicm[i,j]
# 
## the i.th row is for results corresponding to hi
# Since there are 12 hi there are 12 rows.
#
## and j.th colomn is for the j.th row of the data
##
#
#
aicm=matrix(NA,byrow=T,nrow=12,ncol=dim.freq[1])

mum=matrix(NA,byrow=T,nrow=12,ncol=dim.freq[1])
sigmam=matrix(NA,byrow=T,nrow=12,ncol=dim.freq[1])
dispersionm=matrix(NA,byrow=T,nrow=12,ncol=dim.freq[1])
rdm=matrix(NA,byrow=T,nrow=12,ncol=dim.freq[1])
beta1m=matrix(NA,byrow=T,nrow=12,ncol=dim.freq[1])
beta2m=matrix(NA,byrow=T,nrow=12,ncol=dim.freq[1])


for (hi in 1:12)
 {

 for (i in 1:dim.freq[1])
 {
  M=(length(freq.tab[i,])-1)
  Freq=freq.tab[i,]
  rdm[hi,i]=rd(freq.tab[i,])
  mle=mle.optim.semi.altham.freq(hi, Freq, M)
  aicm[hi,i]=mle$AIC
  beta1=mle$MLE[1]; beta2=mle$MLE[2]
    beta1m[hi,i]=beta1
    beta2m[hi,i]=beta2
    beta=c(beta1, beta2)
  PMF=dsemi.altham(hi, M, beta)
 
  moment=msemi.altham(hi, M, beta)
  mum[hi,i]=moment[1]
  sigmam[hi,i]=moment[4]
  dispersionm[hi,i]=moment[5]
  
  expted=PMF*sum(freq.tab[i,])
    
  #f=c(freq.tab[i,1:k], sum(freq.tab[i,(k+1):dim.freq[2]]))
  #e=c(expted[1:k],sum(expted[(k+1):dim.freq[2]]))
  
   f=c(freq.tab[i,])
   e=c(expted)

  testm[hi,i]=sum((f-e)^2/e) 
  # chi-square test statistics
  
  tv=testm[hi,i] 
  pm[hi,i]=1-pchisq(tv,length(e)-2)   
#  k=k+1;
  
  
 }

 }

 list(AIC=aicm, Test.Statistics=testm, P.value=pm, mu=mum, sigma2=sigmam,
      Dispersion=dispersionm,Rd=rdm, beta1=beta1m, beta2=beta2m)

}

###################################################
#############################################
#
cum.example.zelterman=function(tab)
{ 
# This function consider the 
# cumulative sum of the row of the table
# Example: tab6.3
# tab6.3
#     n f_n m_n   0  1  2 3 4 5 6 7 8
#[1,] 1 267  12 255 12  0 0 0 0 0 0 0
#[2,] 2 285  48 239 44  2 0 0 0 0 0 0
#[3,] 3 202  80 143 41 15 3 0 0 0 0 0
#[4,] 4 110  54  69 30  9 2 0 0 0 0 0
#[5,] 5 104 103  43 34 15 9 3 0 0 0 0
#[6,] 6  50  67  15 18  8 5 3 0 1 0 0
#[7,] 7  21  38   4  4  7 4 2 0 0 0 0
#[8,] 8  12  28   1  2  4 3 1 1 0 0 0
#
# tab=tab6.3
#
# x=tab
#       0  1  2 3 4 5 6 7 8
#[1,] 255 12  0 0 0 0 0 0 0
#[2,] 239 44  2 0 0 0 0 0 0
#[3,] 143 41 15 3 0 0 0 0 0
#[4,]  69 30  9 2 0 0 0 0 0
#[5,]  43 34 15 9 3 0 0 0 0
#[6,]  15 18  8 5 3 0 1 0 0
#[7,]   4  4  7 4 2 0 0 0 0
#[8,]   1  2  4 3 1 1 0 0 0
#
# cumx
#     [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9]
#[1,]  255   12    0    0    0    0    0    0    0
#[2,]  494   56    2    0    0    0    0    0    0
#[3,]  637   97   17    3    0    0    0    0    0
#[4,]  706  127   26    5    0    0    0    0    0
#[5,]  749  161   41   14    3    0    0    0    0
#[6,]  764  179   49   19    6    0    1    0    0
#[7,]  768  183   56   23    8    0    1    0    0
#[8,]  769  185   60   26    9    1    1    0    0

dim.tab=dim(tab)
puretab=tab[1:dim.tab[1], 4:dim.tab[2]]
  x=puretab
  dim.x=dim(x)
  m=dim.x[1]
  n=dim.x[2]
  cumx=matrix(NA, nrow=m,ncol=n)
  top=vector(length=n)
  for(i in 1:m)
  {
    top=top+x[i,]
    cumx[i,]=top
  }
 
 freq.tab=cumx
 dim.freq=dim(freq.tab)


## matrix to record chi-square test-statistics
#
#testm[i,j]
# 
# the i.th row is for hi=i
# and j.th colomn is for the j.th row of the data
#

testm=matrix(NA,byrow=T,nrow=12,ncol=dim.freq[1]) 

## matrix for all p-values pm[i,j]
# 
## the i.th row is for hi=i
## and j.th colomn is for the j.th row of the data
##
#
pm=matrix(NA,byrow=T,nrow=12,ncol=dim.freq[1])
#
## matrix for all aic values of the model
#
## aicm[i,j]
# 
## the i.th row is for results corresponding to hi
# Since there are 12 hi there are 12 rows.
#
## and j.th colomn is for the j.th row of the data
##
#
#
aicm=matrix(NA,byrow=T,nrow=12,ncol=dim.freq[1])

mum=matrix(NA,byrow=T,nrow=12,ncol=dim.freq[1])
sigmam=matrix(NA,byrow=T,nrow=12,ncol=dim.freq[1])
dispersionm=matrix(NA,byrow=T,nrow=12,ncol=dim.freq[1])
rdm=matrix(NA,byrow=T,nrow=12,ncol=dim.freq[1])
beta1m=matrix(NA,byrow=T,nrow=12,ncol=dim.freq[1])
beta2m=matrix(NA,byrow=T,nrow=12,ncol=dim.freq[1])


for (hi in 1:12)
 {

 for (i in 1:dim.freq[1])
 {
  M=(length(freq.tab[i,])-1)
  Freq=freq.tab[i,]
  rdm[hi,i]=rd(freq.tab[i,])
  mle=mle.optim.semi.altham.freq(hi, Freq, M)
  aicm[hi,i]=mle$AIC
  beta1=mle$MLE[1]; beta2=mle$MLE[2]
    beta1m[hi,i]=beta1
    beta2m[hi,i]=beta2
    beta=c(beta1, beta2)
  PMF=dsemi.altham(hi, M, beta)
 
  moment=msemi.altham(hi, M, beta)
  mum[hi,i]=moment[1]
  sigmam[hi,i]=moment[4]
  dispersionm[hi,i]=moment[5]
  
  expted=PMF*sum(freq.tab[i,])
    
  #f=c(freq.tab[i,1:k], sum(freq.tab[i,(k+1):dim.freq[2]]))
  #e=c(expted[1:k],sum(expted[(k+1):dim.freq[2]]))
  
   f=c(freq.tab[i,])
   e=c(expted)

  testm[hi,i]=sum((f-e)^2/e) 
  # chi-square test statistics
  
  tv=testm[hi,i] 
  pm[hi,i]=1-pchisq(tv,length(e)-2)   

  
  
 }

 }


 list(AIC=aicm, Test.Statistics=testm, P.value=pm, mu=mum, sigma2=sigmam,
      Dispersion=dispersionm,Rd=rdm, beta1=beta1m, beta2=beta2m)
}





# Test
out.zelt1=example.zelterman(tab6.3)
out.zelt1

#$AIC
#          [,1]     [,2]     [,3]     [,4]     [,5]
# [1,] 207.9221 281.7527 158.6891 70.40504 44.17325
# [2,] 207.9181 281.7245 158.7274 70.38445 44.17102
# [3,] 208.2204 282.3850 157.8071 71.10839 44.21134
# [4,] 207.7969 281.8691 160.1701 69.89362 44.24614
# [5,] 208.0875 282.0801 158.0696 70.83715 44.19593
# [6,] 207.8259 281.7397 159.5955 70.05369 44.20101
# [7,] 207.9151 281.7088 158.7592 70.37116 44.17173
# [8,] 207.9292 281.8004 158.6336 70.44018 44.17683
# [9,] 211.2224 289.5937 158.8368 74.44250 44.99733
#[10,] 208.5224 286.8348 167.5802 69.31670 44.83672
#[11,] 208.1411 282.1246 157.9591 70.93201 44.21171
#[12,] 207.8160 281.8399 159.8013 69.99106 44.19974

#$Test.Statistics
#           [,1]     [,2]      [,3]      [,4]      [,5]
# [1,] 0.3147530 1.795445  4.631586 1.4215201 0.9545620
# [2,] 0.3120557 1.773761  4.692900 1.4044644 0.9534887
# [3,] 0.5310232 2.267535  3.177383 2.0239358 0.9450302
# [4,] 0.2573400 2.080765  7.278837 1.0064429 1.0830028
# [5,] 0.4297079 2.018425  3.628448 1.7884999 0.9473473
# [6,] 0.2631638 1.890514  6.207484 1.1332142 1.0152792
# [7,] 0.3100179 1.762389  4.754134 1.3936526 0.9555241
# [8,] 0.3195577 1.832972  4.546970 1.4509221 0.9561244
# [9,] 3.4693850 9.454096  3.874438 5.2578000 1.7549908
#[10,] 1.1462060 7.736491 25.376262 0.6005265 1.8558298
#[11,] 0.4698456 2.046272  3.475463 1.8721103 0.9594475
#[12,] 0.2636274 2.007549  6.518023 1.0801571 1.0180505

#$P.value
#           [,1]      [,2]         [,3]      [,4]      [,5]
# [1,] 0.9998823 0.9702901 0.7048196842 0.9849030 0.9955286
# [2,] 0.9998857 0.9712954 0.6973814255 0.9854341 0.9955443
# [3,] 0.9993248 0.9435625 0.8681222196 0.9585037 0.9956672
# [4,] 0.9999405 0.9552352 0.4004363067 0.9947239 0.9933762
# [5,] 0.9996653 0.9588136 0.8214391221 0.9706142 0.9956338
# [6,] 0.9999358 0.9656507 0.5157423637 0.9923842 0.9945783
# [7,] 0.9998882 0.9718148 0.6899381218 0.9857649 0.9955144
# [8,] 0.9998762 0.9685039 0.7150541909 0.9839603 0.9955056
# [9,] 0.8384575 0.2216760 0.7941208588 0.6285349 0.9721498
#[10,] 0.9921132 0.3564024 0.0006506151 0.9989889 0.9673871
#[11,] 0.9995496 0.9572349 0.8378174916 0.9665783 0.9954564
#[12,] 0.9999355 0.9594214 0.4807232247 0.9934298 0.9945321

#$mu
#           [,1]      [,2]     [,3]     [,4]     [,5]
# [1,] 0.4909105 0.9903955 1.340166 1.809531 2.333336
# [2,] 0.4909090 0.9903850 1.340033 1.809529 2.333336
# [3,] 0.4909087 0.9903892 1.339968 1.809798 2.333342
# [4,] 0.4909436 0.9903872 1.340031 1.809526 2.333342
# [5,] 0.4909097 0.9903893 1.339983 1.809520 2.333331
# [6,] 0.4909101 0.9904415 1.339999 1.809523 2.333343
# [7,] 0.4909026 0.9905158 1.339992 1.809527 2.333334
# [8,] 0.4909295 0.9902838 1.340057 1.809550 2.333346
# [9,] 0.4909059 0.9904031 1.339982 1.809539 2.333309
#[10,] 0.4909100 0.9904068 1.340006 1.809738 2.333358
#[11,] 0.4909110 0.9903997 1.340006 1.809532 2.333329
#[12,] 0.4908369 0.9903848 1.340000 1.809531 2.333342

#$sigma2
#           [,1]     [,2]     [,3]     [,4]     [,5]
# [1,] 0.5226437 1.163361 1.784376 1.487496 1.722143
# [2,] 0.5226421 1.163343 1.784370 1.487485 1.722146
# [3,] 0.5226402 1.163349 1.784006 1.486336 1.722129
# [4,] 0.5226714 1.163354 1.784361 1.487492 1.722160
# [5,] 0.5226427 1.163344 1.784205 1.487468 1.722197
# [6,] 0.5226438 1.163396 1.784327 1.487529 1.722204
# [7,] 0.5226454 1.163360 1.784444 1.487484 1.722142
# [8,] 0.5226688 1.163266 1.783732 1.487375 1.722198
# [9,] 0.5226130 1.163333 1.784227 1.487468 1.722272
#[10,] 0.5226439 1.163436 1.784370 1.486731 1.722148
#[11,] 0.5226438 1.163363 1.784278 1.487498 1.722236
#[12,] 0.5226779 1.163348 1.784319 1.487496 1.722192

#$Dispersion
#          [,1]     [,2]     [,3]     [,4]     [,5]
# [1,] 1.134243 1.340609 1.599390 1.062322 1.041968
# [2,] 1.134243 1.340601 1.599511 1.062315 1.041970
# [3,] 1.134239 1.340603 1.599247 1.061383 1.041958
# [4,] 1.134232 1.340612 1.599505 1.062321 1.041977
# [5,] 1.134243 1.340598 1.599411 1.062307 1.042002
# [6,] 1.134244 1.340597 1.599505 1.062349 1.042004
# [7,] 1.134264 1.340469 1.599616 1.062315 1.041968
# [8,] 1.134256 1.340630 1.598917 1.062228 1.041999
# [9,] 1.134186 1.340569 1.599432 1.062299 1.042053
#[10,] 1.134244 1.340683 1.599538 1.061689 1.041966
#[11,] 1.134242 1.340607 1.599454 1.062323 1.042027
#[12,] 1.134476 1.340607 1.599497 1.062322 1.041996

#$Rd
#          [,1]     [,2]     [,3]     [,4]     [,5]
# [1,] 1.144654 1.353647 1.632214 1.115466 1.136746
# [2,] 1.144654 1.353647 1.632214 1.115466 1.136746
# [3,] 1.144654 1.353647 1.632214 1.115466 1.136746
# [4,] 1.144654 1.353647 1.632214 1.115466 1.136746
# [5,] 1.144654 1.353647 1.632214 1.115466 1.136746
# [6,] 1.144654 1.353647 1.632214 1.115466 1.136746
# [7,] 1.144654 1.353647 1.632214 1.115466 1.136746
# [8,] 1.144654 1.353647 1.632214 1.115466 1.136746
# [9,] 1.144654 1.353647 1.632214 1.115466 1.136746
#[10,] 1.144654 1.353647 1.632214 1.115466 1.136746
#[11,] 1.144654 1.353647 1.632214 1.115466 1.136746
#[12,] 1.144654 1.353647 1.632214 1.115466 1.136746

#$beta1
#             [,1]        [,2]        [,3]       [,4]        [,5]
# [1,]  0.58604925  0.12301722  0.05325056 -0.8709711 -1.16108662
# [2,]  0.55307365  0.09347444  0.03083708 -0.8979184 -1.18208604
# [3,]  2.98746693  2.41281965  2.27033389  1.3289787  0.95855460
# [4,] -1.82610933 -2.18055913 -2.17284086 -3.0866180 -3.29670237
# [5,]  1.36654230  0.83078178  0.70666569 -0.2193170 -0.56607804
# [6,] -0.19872543 -0.59206396 -0.60816935 -1.5296845 -1.76511667
# [7,]  0.47827942  0.02047129 -0.03906356 -0.9688350 -1.25003125
# [8,]  0.69383246  0.22503623  0.14381369 -0.7736453 -1.07310560
# [9,]  2.93217740  2.05653791  1.74274221  0.6897412  0.09468177
#[10,] -1.95820408 -1.96450149 -1.74151897 -2.6141704 -2.63451360
#[11,]  1.19496596  0.64411014  0.51328451 -0.4175201 -0.77088457
#[12,] -0.02855136 -0.40683266 -0.41559972 -1.3330841 -1.56212372

#$beta2
#            [,1]        [,2]        [,3]        [,4]        [,5]
# [1,]  0.2206062  0.15447729  0.10583195  0.26495217  0.25709855
# [2,]  0.2164599  0.14897980  0.09823305  0.25892733  0.24974626
# [3,] -0.1318766 -0.15356304 -0.18039418 -0.02362266 -0.01396988
# [4,]  0.5799235  0.46847862  0.39522842  0.55756079  0.53114440
# [5,]  0.1144536  0.07749921  0.04522539  0.20024202  0.20483639
# [6,]  0.3294993  0.23448793  0.16905975  0.33154080  0.31114492
# [7,]  0.2113619  0.14310491  0.09145706  0.25260508  0.24275634
# [8,]  0.2298357  0.16603631  0.12077244  0.27744250  0.27166285
# [9,] -0.2249124 -0.14783929 -0.13147660  0.03038228  0.08361334
#[10,]  0.7952024  0.52701214  0.37906811  0.54335462  0.46831756
#[11,]  0.1139205  0.08297744  0.05270005  0.20827782  0.21430188
#[12,]  0.3311568  0.22970509  0.16184887  0.32386042  0.30194883


outcum=cum.example.zelterman(tab6.3)
outcum

#$AIC
#          [,1]     [,2]     [,3]     [,4]     [,5]
# [1,] 207.9221 501.5176 668.2563 749.0312 804.0210
# [2,] 207.9181 501.4890 668.2858 749.0143 803.9695
# [3,] 208.2204 502.3957 668.0800 750.1985 805.7778
# [4,] 207.7969 501.7368 670.3925 750.1757 804.9182
# [5,] 208.0875 501.9544 667.9251 749.5430 804.9221
# [6,] 207.8259 501.5371 669.4471 749.5406 804.3152
# [7,] 207.9151 501.4729 668.3109 749.0123 803.9543
# [8,] 207.9292 501.5667 668.2200 749.0742 804.1223
# [9,] 211.2224 514.0188 681.7909 770.1632 829.3754
#[10,] 208.5224 508.9710 687.3546 765.8611 821.3314
#[11,] 208.1411 502.0495 667.8980 749.6895 805.1438
#[12,] 207.8160 501.6668 669.8388 749.8158 804.5642
#
#$Test.Statistics
#           [,1]      [,2]      [,3]      [,4]      [,5]
# [1,] 0.3147530  1.510345  4.490652  3.813685  1.970786
# [2,] 0.3120557  1.490163  4.584695  3.831977  1.944016
# [3,] 0.5310232  2.153394  3.507977  4.424467  3.483451
# [4,] 0.2573400  1.932940  9.031735  6.469659  3.484877
# [5,] 0.4297079  1.797528  3.555785  3.924210  2.676719
# [6,] 0.2631638  1.664706  6.981055  5.142863  2.603349
# [7,] 0.3100179  1.479524  4.656416  3.857360  1.941838
# [8,] 0.3195577  1.546029  4.361952  3.806091  2.033846
# [9,] 3.4693850 13.785312 17.379173 24.632568 27.148657
#[10,] 1.1462060 10.032943 88.247178 56.060068 35.072841
#[11,] 0.4698456  1.864624  3.464963  4.018345  2.874871
#[12,] 0.2636274  1.815410  7.634854  5.583079  2.914691

#$P.value
#           [,1]       [,2]         [,3]         [,4]         [,5]
# [1,] 0.9998823 0.98194936 7.218427e-01 8.009698e-01 9.614396e-01
# [2,] 0.9998857 0.98264826 7.104960e-01 7.989143e-01 9.628738e-01
# [3,] 0.9993248 0.95086450 8.343798e-01 7.297926e-01 8.369752e-01
# [4,] 0.9999405 0.96345851 2.503856e-01 4.860994e-01 8.368246e-01
# [5,] 0.9996653 0.97019251 8.292820e-01 7.884644e-01 9.132129e-01
# [6,] 0.9999358 0.97605411 4.308551e-01 6.425339e-01 9.191158e-01
# [7,] 0.9998882 0.98301009 7.018094e-01 7.960525e-01 9.629892e-01
# [8,] 0.9998762 0.98067321 7.372704e-01 8.018214e-01 9.579433e-01
# [9,] 0.8384575 0.05513422 1.510819e-02 8.813831e-04 3.133573e-04
#[10,] 0.9921132 0.18671428 3.330669e-16 9.188184e-10 1.083733e-05
#[11,] 0.9995496 0.96695159 8.389225e-01 7.776620e-01 8.963313e-01
#[12,] 0.9999355 0.96934714 3.658899e-01 5.891828e-01 8.927802e-01

#$mu
#           [,1]      [,2]      [,3]      [,4]      [,5]
# [1,] 0.4909105 0.7335244 0.8484847 0.9192979 0.9764358
# [2,] 0.4909090 0.7336473 0.8484889 0.9193086 0.9764670
# [3,] 0.4909087 0.7336446 0.8484860 0.9193080 0.9764296
# [4,] 0.4909436 0.7336500 0.8484875 0.9192989 0.9763059
# [5,] 0.4909097 0.7336391 0.8484863 0.9192999 0.9764449
# [6,] 0.4909101 0.7336444 0.8484941 0.9193027 0.9764338
# [7,] 0.4909026 0.7336453 0.8484922 0.9192965 0.9764325
# [8,] 0.4909295 0.7336391 0.8485000 0.9191753 0.9764524
# [9,] 0.4909059 0.7335907 0.8486897 0.9192985 0.9764848
#[10,] 0.4909100 0.7336345 0.8484880 0.9192421 0.9764699
#[11,] 0.4909110 0.7336526 0.8484847 0.9193038 0.9764503
#[12,] 0.4908369 0.7336447 0.8484845 0.9193121 0.9766711

#$sigma2
#           [,1]      [,2]     [,3]     [,4]     [,5]
# [1,] 0.5226437 0.8963273 1.120948 1.210990 1.309170
# [2,] 0.5226421 0.8963351 1.120949 1.211003 1.309194
# [3,] 0.5226402 0.8963202 1.120945 1.210998 1.309133
# [4,] 0.5226714 0.8963411 1.120963 1.211003 1.309079
# [5,] 0.5226427 0.8963471 1.120941 1.211028 1.309133
# [6,] 0.5226438 0.8963400 1.120975 1.211013 1.309173
# [7,] 0.5226454 0.8963339 1.120965 1.211010 1.309167
# [8,] 0.5226688 0.8963145 1.120978 1.210932 1.309156
# [9,] 0.5226130 0.8962556 1.120661 1.210896 1.309240
#[10,] 0.5226439 0.8963303 1.120964 1.211149 1.309408
#[11,] 0.5226438 0.8963390 1.120936 1.210990 1.309187
#[12,] 0.5226779 0.8963313 1.120955 1.210983 1.309388

#$Dispersion
#          [,1]     [,2]     [,3]     [,4]     [,5]
# [1,] 1.134243 1.345297 1.477861 1.488326 1.527160
# [2,] 1.134243 1.345106 1.477855 1.488326 1.527147
# [3,] 1.134239 1.345088 1.477855 1.488321 1.527126
# [4,] 1.134232 1.345111 1.477875 1.488340 1.527229
# [5,] 1.134243 1.345138 1.477849 1.488370 1.527106
# [6,] 1.134244 1.345118 1.477881 1.488347 1.527167
# [7,] 1.134264 1.345108 1.477872 1.488352 1.527162
# [8,] 1.134256 1.345089 1.477876 1.488427 1.527123
# [9,] 1.134186 1.345080 1.477167 1.488209 1.527177
#[10,] 1.134244 1.345120 1.477877 1.488599 1.527392
#[11,] 1.134242 1.345103 1.477844 1.488317 1.527162
#[12,] 1.134476 1.345105 1.477870 1.488297 1.527098

#$Rd
#          [,1]     [,2]     [,3]     [,4]     [,5]
# [1,] 1.144654 1.351440 1.483525 1.493616 1.532372
# [2,] 1.144654 1.351440 1.483525 1.493616 1.532372
# [3,] 1.144654 1.351440 1.483525 1.493616 1.532372
# [4,] 1.144654 1.351440 1.483525 1.493616 1.532372
# [5,] 1.144654 1.351440 1.483525 1.493616 1.532372
# [6,] 1.144654 1.351440 1.483525 1.493616 1.532372
# [7,] 1.144654 1.351440 1.483525 1.493616 1.532372
# [8,] 1.144654 1.351440 1.483525 1.493616 1.532372
# [9,] 1.144654 1.351440 1.483525 1.493616 1.532372
#[10,] 1.144654 1.351440 1.483525 1.493616 1.532372
#[11,] 1.144654 1.351440 1.483525 1.493616 1.532372
#[12,] 1.144654 1.351440 1.483525 1.493616 1.532372

#$beta1
#             [,1]       [,2]        [,3]       [,4]       [,5]
# [1,]  0.58604925  0.4515826  0.43056165  0.3542027  0.3245841
# [2,]  0.55307365  0.4205509  0.40226259  0.3264196  0.2977727
# [3,]  2.98746693  2.7725559  2.71509927  2.6294238  2.5882320
# [4,] -1.82610933 -1.8854076 -1.86821121 -1.9345871 -1.9514840
# [5,]  1.36654230  1.1797651  1.13272335  1.0500050  1.0117126
# [6,] -0.19872543 -0.2849127 -0.28047133 -0.3502923 -0.3714091
# [7,]  0.47827942  0.3471761  0.32996441  0.2544024  0.2261666
# [8,]  0.69383246  0.5548837  0.53034585  0.4533796  0.4218574
# [9,]  2.93217740  2.5040389  2.36551106  2.2553466  2.1873197
#[10,] -1.95820408 -1.7771952 -1.65078474 -1.6855329 -1.6664288
#[11,]  1.19496596  0.9971853  0.94622559  0.8623814  0.8230213
#[12,] -0.02855136 -0.1043378 -0.09553213 -0.1641972 -0.1842453

#$beta2
#            [,1]        [,2]        [,3]        [,4]        [,5]
# [1,]  0.2206062  0.12791281  0.09629834  0.10113501  0.09643816
# [2,]  0.2164599  0.12283156  0.09027105  0.09496008  0.08995259
# [3,] -0.1318766 -0.18905586 -0.20742640 -0.20009280 -0.20127361
# [4,]  0.5799235  0.45260789  0.40645801  0.40827481  0.39945327
# [5,]  0.1144536  0.04496080  0.02287425  0.02944336  0.02740542
# [6,]  0.3294993  0.21466744  0.17347507  0.17636835  0.16892490
# [7,]  0.2113619  0.11710108  0.08414231  0.08875420  0.08362201
# [8,]  0.2298357  0.13904609  0.10875262  0.11378829  0.10964855
# [9,] -0.2249124 -0.20052204 -0.19429286 -0.18220454 -0.17686314
#[10,]  0.7952024  0.55057187  0.45940120  0.45029536  0.42846136
#[11,]  0.1139205  0.04928995  0.02861625  0.03547934  0.03375809
#[12,]  0.3311568  0.21139902  0.16850489  0.17103548  0.16317451

# Remind the dispersion from the observations data
#x=tab6.3; row0=3
#rd.by.row(x, row0)
#[1] 1.393973 1.224720 1.479030 1.749551 1.164075 1.136746
#[1] 1.265661 1.270871 1.464611 1.602024 1.640585 1.688020


# example.zelterman(tab6.1)

out.zelt2=example.zelterman(tab6.1)
out.zelt2

#$AIC
#          [,1]     [,2]
# [1,] 18.20969 23.45910
# [2,] 18.20876 23.61038
# [3,] 18.23977 23.78514
# [4,] 18.18389 23.15973
# [5,] 18.22995 23.57234
# [6,] 18.19149 23.37561
# [7,] 18.20824 23.57234
# [8,] 18.21117 23.37561
# [9,] 18.38019 24.98340
#[10,] 18.12212 22.37547
#[11,] 18.23701 23.61038
#[12,] 18.18595 23.34431

#$Test.Statistics
#            [,1]     [,2]
# [1,] 0.09761512 2.000000
# [2,] 0.09689905 2.116182
# [3,] 0.12098766 2.233129
# [4,] 0.07844120 1.794612
# [5,] 0.11324467 2.088089
# [6,] 0.08399027 1.949820
# [7,] 0.09649987 2.088089
# [8,] 0.09863173 1.949824
# [9,] 0.24308237 3.167092
#[10,] 0.03614894 1.291963
#[11,] 0.11881105 2.116189
#[12,] 0.07993661 1.929008

#$P.value
#           [,1]      [,2]
# [1,] 0.9998471 0.8491450
# [2,] 0.9998498 0.8328468
# [3,] 0.9997406 0.8160366
# [4,] 0.9999109 0.8767712
# [5,] 0.9997795 0.8368269
# [6,] 0.9998945 0.8560439
# [7,] 0.9998513 0.8368270
# [8,] 0.9998431 0.8560433
# [9,] 0.9985787 0.6742429
#[10,] 0.9999870 0.9357565
#[11,] 0.9997519 0.8328458
#[12,] 0.9999066 0.8588784

#$mu
#           [,1]     [,2]
# [1,] 0.7142880 3.000000
# [2,] 0.7142870 3.000060
# [3,] 0.7142934 3.000000
# [4,] 0.7142897 3.000056
# [5,] 0.7142773 2.999987
# [6,] 0.7142904 3.000001
# [7,] 0.7142866 3.000042
# [8,] 0.7142221 3.000027
# [9,] 0.7142844 3.000193
#[10,] 0.7142862 2.999999
#[11,] 0.7142786 3.000012
#[12,] 0.7142897 3.000003

#$sigma2
#           [,1]     [,2]
# [1,] 0.4897937 4.000000
# [2,] 0.4897948 3.999881
# [3,] 0.4897962 3.999789
# [4,] 0.4897860 3.999802
# [5,] 0.4898021 4.000320
# [6,] 0.4897952 4.000130
# [7,] 0.4897957 3.999791
# [8,] 0.4899654 3.999681
# [9,] 0.4897993 3.999058
#[10,] 0.4897953 3.999752
#[11,] 0.4897755 3.999899
#[12,] 0.4897945 3.999828

#$Dispersion
#           [,1]     [,2]
# [1,] 0.7783727 2.666667
# [2,] 0.7783753 2.666587
# [3,] 0.7783716 2.666526
# [4,] 0.7783589 2.666535
# [5,] 0.7783961 2.666880
# [6,] 0.7783728 2.666753
# [7,] 0.7783771 2.666528
# [8,] 0.7787076 2.666454
# [9,] 0.7783849 2.666038
#[10,] 0.7783770 2.666501
#[11,] 0.7783526 2.666599
#[12,] 0.7783725 2.666552

#$Rd
#           [,1]     [,2]
# [1,] 0.9081081 3.333333
# [2,] 0.9081081 3.333333
# [3,] 0.9081081 3.333333
# [4,] 0.9081081 3.333333
# [5,] 0.9081081 3.333333
# [6,] 0.9081081 3.333333
# [7,] 0.9081081 3.333333
# [8,] 0.9081081 3.333333
# [9,] 0.9081081 3.333333
#[10,] 0.9081081 3.333333
#[11,] 0.9081081 3.333333
#[12,] 0.9081081 3.333333


#$beta1
#             [,1]        [,2]
# [1,] -0.77746420  0.00000000
# [2,] -0.82308899  0.06728445
# [3,]  1.41530100  2.01660731
# [4,] -2.97407737 -2.01642517
# [5,]  0.03464193  0.58578094
# [6,] -1.59133846 -0.58745400
# [7,] -0.91562699 -0.02424090
# [8,] -0.63806938  0.02336125
# [9,]  1.51603077  1.60115931
#[10,] -3.14270034 -1.59741917
#[11,] -0.12838516  0.38481299
#[12,] -1.42976249 -0.38757650

#$beta2
#           [,1]          [,2]
# [1,] 0.7040806  7.180547e-15
# [2,] 0.6963496 -3.767577e-02
# [3,] 0.2916423 -3.361012e-01
# [4,] 1.1192170  3.360685e-01
# [5,] 0.5773919 -4.681756e-02
# [6,] 0.8320360  4.698828e-02
# [7,] 0.6886972 -4.677311e-02
# [8,] 0.7188748  4.702556e-02
# [9,] 0.1434515 -2.668679e-01
#[10,] 1.3154188  2.662366e-01
#[11,] 0.5711815 -3.767680e-02
#[12,] 0.8392506  3.797387e-02

outcum.zelt2=cum.example.zelterman(tab6.1)
outcum.zelt2

#$AIC
#          [,1]     [,2]     [,3]     [,4]     [,5]
# [1,] 57.98424 96.38004 159.9677 175.5814 215.4275
# [2,] 57.98497 96.37993 159.9732 175.5868 216.0936
# [3,] 57.98447 96.38391 159.7084 175.3420 214.4540
# [4,] 57.98487 96.37697 160.3102 175.9073 216.6978
# [5,] 57.98417 96.38261 159.7868 175.4129 214.2157
# [6,] 57.98488 96.37785 160.1916 175.7940 216.8981
# [7,] 57.98489 96.37986 159.9776 175.5910 216.0892
# [8,] 57.98417 96.38022 159.9586 175.5725 214.8188
# [9,] 57.98429 96.40611 159.4979 175.2379 213.2904
#[10,] 57.98489 96.37061 161.9786 177.5408 222.5633
#[11,] 57.98418 96.38355 159.7289 175.3601 214.0139
#[12,] 57.98488 96.37720 160.2776 175.8762 217.2103

#$Test.Statistics
#              [,1]        [,2]     [,3]     [,4]      [,5]
# [1,] 3.256171e-05 0.007432135 2.064911 1.648025  5.657344
# [2,] 3.970618e-04 0.007371015 2.069734 1.654363  7.127987
# [3,] 1.477130e-04 0.009511985 1.784059 1.372811  4.259511
# [4,] 3.483536e-04 0.005807076 2.436194 2.014538  7.510986
# [5,] 5.890744e-08 0.008807900 1.871169 1.455851  3.783919
# [6,] 3.544463e-04 0.006271018 2.306205 1.887229  8.406479
# [7,] 3.557403e-04 0.007338584 2.074312 1.659364  7.092161
# [8,] 5.011667e-09 0.007526809 2.056878 1.641562  4.579169
# [9,] 5.869434e-05 0.022023771 1.528614 1.142920  2.324811
#[10,] 3.567851e-04 0.002509311 4.381677 3.957669 17.276452
#[11,] 7.693053e-07 0.009312696 1.808024 1.394375  3.519947
#[12,] 3.511505e-04 0.005931107 2.400487 1.978418  8.954340

#$P.value
#      [,1]      [,2]      [,3]      [,4]        [,5]
# [1,]    1 0.9999997 0.8400924 0.8953816 0.341000805
# [2,]    1 0.9999998 0.8394143 0.8945983 0.211294418
# [3,]    1 0.9999995 0.8781442 0.9272642 0.512689330
# [4,]    1 0.9999999 0.7860723 0.8471298 0.185325282
# [5,]    1 0.9999996 0.8666686 0.9181077 0.580926589
# [6,]    1 0.9999998 0.8053551 0.8645187 0.135211015
# [7,]    1 0.9999998 0.8387699 0.8939789 0.213875908
# [8,]    1 0.9999997 0.8412203 0.8961782 0.469361228
# [9,]    1 0.9999962 0.9097444 0.9502347 0.802616008
#[10,]    1 1.0000000 0.4958689 0.5555266 0.004004270
#[11,]    1 0.9999996 0.8750193 0.9249276 0.620371532
#[12,]    1 0.9999999 0.7914016 0.8521232 0.110899297

#$mu
#           [,1]      [,2]      [,3]      [,4]      [,5]
# [1,] 0.2499996 0.2957746 0.4545437 0.4737305 0.6000074
# [2,] 0.2499729 0.2957753 0.4545497 0.4736896 0.6000032
# [3,] 0.2499967 0.2957746 0.4545469 0.4736847 0.6000145
# [4,] 0.2499785 0.2957745 0.4545378 0.4736583 0.6000110
# [5,] 0.2500028 0.2957747 0.4545613 0.4737009 0.6000062
# [6,] 0.2500066 0.2957733 0.4545468 0.4736846 0.6000045
# [7,] 0.2499780 0.2957746 0.4545452 0.4736840 0.6000045
# [8,] 0.2499999 0.2957747 0.4545452 0.4735473 0.6000170
# [9,] 0.2499986 0.2957760 0.4545463 0.4736851 0.6000161
#[10,] 0.2499780 0.2957746 0.4545479 0.4736844 0.5999995
#[11,] 0.2500000 0.2957746 0.4544942 0.4737003 0.5999397
#[12,] 0.2499783 0.2957737 0.4544842 0.4737287 0.6001217

#$sigma2
#           [,1]      [,2]      [,3]      [,4]      [,5]
# [1,] 0.1875012 0.2364609 0.4752134 0.4809492 0.9599829
# [2,] 0.1875030 0.2364629 0.4752148 0.4809147 0.9599745
# [3,] 0.1875045 0.2364609 0.4752061 0.4808848 0.9600011
# [4,] 0.1875038 0.2364608 0.4752006 0.4808677 0.9599962
# [5,] 0.1875014 0.2364610 0.4752148 0.4808951 0.9599750
# [6,] 0.1875181 0.2364607 0.4752066 0.4808845 0.9599826
# [7,] 0.1875038 0.2364609 0.4752041 0.4808840 0.9599821
# [8,] 0.1874999 0.2364610 0.4752041 0.4805236 0.9600077
# [9,] 0.1875017 0.2364563 0.4752013 0.4808876 0.9599983
#[10,] 0.1875038 0.2364609 0.4752091 0.4808854 0.9600082
#[11,] 0.1875000 0.2364609 0.4751631 0.4809339 0.9601147
#[12,] 0.1875038 0.2364607 0.4751862 0.4809801 0.9601598

#$Dispersion
#           [,1]      [,2]     [,3]     [,4]     [,5]
# [1,] 0.7826147 0.8409167 1.131168 1.102268 1.777727
# [2,] 0.7827023 0.8409221 1.131157 1.102276 1.777722
# [3,] 0.7826374 0.8409167 1.131143 1.102217 1.777742
# [4,] 0.7826887 0.8409167 1.131151 1.102234 1.777742
# [5,] 0.7826061 0.8409167 1.131131 1.102206 1.777715
# [6,] 0.7826643 0.8409195 1.131144 1.102217 1.777734
# [7,] 0.7826904 0.8409167 1.131142 1.102217 1.777733
# [8,] 0.7826088 0.8409169 1.131142 1.101682 1.777747
# [9,] 0.7826202 0.8408968 1.131133 1.102223 1.777732
#[10,] 0.7826907 0.8409168 1.131148 1.102219 1.777794
#[11,] 0.7826089 0.8409168 1.131161 1.102297 1.778149
#[12,] 0.7826894 0.8409186 1.131239 1.102342 1.777753

#$Rd
#         [,1]      [,2]     [,3]     [,4]     [,5]
# [1,] 0.79926 0.8529302 1.144149 1.113948 1.795735
# [2,] 0.79926 0.8529302 1.144149 1.113948 1.795735
# [3,] 0.79926 0.8529302 1.144149 1.113948 1.795735
# [4,] 0.79926 0.8529302 1.144149 1.113948 1.795735
# [5,] 0.79926 0.8529302 1.144149 1.113948 1.795735
# [6,] 0.79926 0.8529302 1.144149 1.113948 1.795735
# [7,] 0.79926 0.8529302 1.144149 1.113948 1.795735
# [8,] 0.79926 0.8529302 1.144149 1.113948 1.795735
# [9,] 0.79926 0.8529302 1.144149 1.113948 1.795735
#[10,] 0.79926 0.8529302 1.144149 1.113948 1.795735
#[11,] 0.79926 0.8529302 1.144149 1.113948 1.795735
#[12,] 0.79926 0.8529302 1.144149 1.113948 1.795735

#$beta1
#           [,1]         [,2]         [,3]       [,4]       [,5]
# [1,] -4.760716 -0.004111614  0.620137873  0.5147770  1.0220229
# [2,] -3.557236 -0.050824939  0.576513991  0.4709763  1.0015897
# [3,] -1.775139  2.221979473  2.779529286  2.6778562  3.1021004
# [4,] -5.805263 -2.230884468 -1.545445444 -1.6540471 -1.0566230
# [5,] -7.099127  0.830615498  1.407713458  1.3052483  1.7315242
# [6,] -4.403951 -0.839141864 -0.171414084 -0.2792349  0.3100807
# [7,] -3.705345 -0.144021521  0.484804878  0.3790670  0.9125054
# [8,] -9.010535  0.135797576  0.755165896  0.6495645  1.1293453
# [9,] -1.981204  2.458337918  2.761982915  2.6716278  2.8951260
#[10,] -6.048720 -2.480586662 -1.621107988 -1.7393126 -0.8262847
#[11,] -5.950884  0.675348051  1.238237520  1.1364340  1.5474496
#[12,] -4.254311 -0.684126584 -0.003303165 -0.1119145  0.4941999

#$beta2
#           [,1]      [,2]       [,3]       [,4]        [,5]
# [1,]  5.859337 0.9902609  0.2400165  0.2735798 -0.01542613
# [2,]  4.602325 0.9832417  0.2309682  0.2647061 -0.03454985
# [3,]  4.665557 0.5556076 -0.1493377 -0.1181431 -0.36630171
# [4,]  5.112298 1.4255325  0.6336413  0.6694386  0.33572184
# [5,]  8.890871 0.8484665  0.1300934  0.1618337 -0.08866799
# [6,]  4.809450 1.1323355  0.3525504  0.3877905  0.05944858
# [7,]  4.649992 0.9760346  0.2227611  0.2566170 -0.04405474
# [8,] 10.263299 1.0044871  0.2574420  0.2912793  0.01419917
# [9,]  4.871594 0.3174511 -0.2150562 -0.1914662 -0.34217223
#[10,]  5.355759 1.6757785  0.7653362  0.8075402  0.31562779
#[11,]  7.576086 0.8370990  0.1285363  0.1597565 -0.08315513
#[12,]  4.826518 1.1439303  0.3553068  0.3909365  0.05420445



# example.zelterman(tab6.2)
# Error in optim(para0, nlk.semi.altham, gr = NULL, method = "BFGS", hessian = TRUE,  : 
#  non-finite finite-difference value [2]


# example.zelterman(tab6.4)
# Error in optim(para0, nlk.semi.altham, gr = NULL, method = "BFGS", hessian = TRUE,  : 
#  non-finite finite-difference value [2]

out.zelt3=example.zelterman(tab6.9)
out.zelt3

#$AIC
#          [,1]     [,2]     [,3]     [,4]     [,5]
# [1,] 207.9221 281.7527 158.6891 70.40504 44.17325
# [2,] 207.9181 281.7245 158.7274 70.38445 44.17102
# [3,] 208.2204 282.3850 157.8071 71.10839 44.21134
# [4,] 207.7969 281.8691 160.1701 69.89362 44.24614
# [5,] 208.0875 282.0801 158.0696 70.83715 44.19593
# [6,] 207.8259 281.7397 159.5955 70.05369 44.20101
# [7,] 207.9151 281.7088 158.7592 70.37116 44.17173
# [8,] 207.9292 281.8004 158.6336 70.44018 44.17683
# [9,] 211.2224 289.5937 158.8368 74.44250 44.99733
#[10,] 208.5224 286.8348 167.5802 69.31670 44.83672
#[11,] 208.1411 282.1246 157.9591 70.93201 44.21171
#[12,] 207.8160 281.8399 159.8013 69.99106 44.19974

#$Test.Statistics
#           [,1]     [,2]      [,3]      [,4]      [,5]
# [1,] 0.3147530 1.795445  4.631586 1.4215201 0.9545620
# [2,] 0.3120557 1.773761  4.692900 1.4044644 0.9534887
# [3,] 0.5310232 2.267535  3.177383 2.0239358 0.9450302
# [4,] 0.2573400 2.080765  7.278837 1.0064429 1.0830028
# [5,] 0.4297079 2.018425  3.628448 1.7884999 0.9473473
# [6,] 0.2631638 1.890514  6.207484 1.1332142 1.0152792
# [7,] 0.3100179 1.762389  4.754134 1.3936526 0.9555241
# [8,] 0.3195577 1.832972  4.546970 1.4509221 0.9561244
# [9,] 3.4693850 9.454096  3.874438 5.2578000 1.7549908
#[10,] 1.1462060 7.736491 25.376262 0.6005265 1.8558298
#[11,] 0.4698456 2.046272  3.475463 1.8721103 0.9594475
#[12,] 0.2636274 2.007549  6.518023 1.0801571 1.0180505

#$P.value
#           [,1]      [,2]         [,3]      [,4]      [,5]
# [1,] 0.9998823 0.9702901 0.7048196842 0.9849030 0.9955286
# [2,] 0.9998857 0.9712954 0.6973814255 0.9854341 0.9955443
# [3,] 0.9993248 0.9435625 0.8681222196 0.9585037 0.9956672
# [4,] 0.9999405 0.9552352 0.4004363067 0.9947239 0.9933762
# [5,] 0.9996653 0.9588136 0.8214391221 0.9706142 0.9956338
# [6,] 0.9999358 0.9656507 0.5157423637 0.9923842 0.9945783
# [7,] 0.9998882 0.9718148 0.6899381218 0.9857649 0.9955144
# [8,] 0.9998762 0.9685039 0.7150541909 0.9839603 0.9955056
# [9,] 0.8384575 0.2216760 0.7941208588 0.6285349 0.9721498
#[10,] 0.9921132 0.3564024 0.0006506151 0.9989889 0.9673871
#[11,] 0.9995496 0.9572349 0.8378174916 0.9665783 0.9954564
#[12,] 0.9999355 0.9594214 0.4807232247 0.9934298 0.9945321

#$mu
#           [,1]      [,2]     [,3]     [,4]     [,5]
# [1,] 0.4909105 0.9903955 1.340166 1.809531 2.333336
# [2,] 0.4909090 0.9903850 1.340033 1.809529 2.333336
# [3,] 0.4909087 0.9903892 1.339968 1.809798 2.333342
# [4,] 0.4909436 0.9903872 1.340031 1.809526 2.333342
# [5,] 0.4909097 0.9903893 1.339983 1.809520 2.333331
# [6,] 0.4909101 0.9904415 1.339999 1.809523 2.333343
# [7,] 0.4909026 0.9905158 1.339992 1.809527 2.333334
# [8,] 0.4909295 0.9902838 1.340057 1.809550 2.333346
# [9,] 0.4909059 0.9904031 1.339982 1.809539 2.333309
#[10,] 0.4909100 0.9904068 1.340006 1.809738 2.333358
#[11,] 0.4909110 0.9903997 1.340006 1.809532 2.333329
#[12,] 0.4908369 0.9903848 1.340000 1.809531 2.333342

#$sigma2
#           [,1]     [,2]     [,3]     [,4]     [,5]
# [1,] 0.5226437 1.163361 1.784376 1.487496 1.722143
# [2,] 0.5226421 1.163343 1.784370 1.487485 1.722146
# [3,] 0.5226402 1.163349 1.784006 1.486336 1.722129
# [4,] 0.5226714 1.163354 1.784361 1.487492 1.722160
# [5,] 0.5226427 1.163344 1.784205 1.487468 1.722197
# [6,] 0.5226438 1.163396 1.784327 1.487529 1.722204
# [7,] 0.5226454 1.163360 1.784444 1.487484 1.722142
# [8,] 0.5226688 1.163266 1.783732 1.487375 1.722198
# [9,] 0.5226130 1.163333 1.784227 1.487468 1.722272
#[10,] 0.5226439 1.163436 1.784370 1.486731 1.722148
#[11,] 0.5226438 1.163363 1.784278 1.487498 1.722236
#[12,] 0.5226779 1.163348 1.784319 1.487496 1.722192

#$Dispersion
#          [,1]     [,2]     [,3]     [,4]     [,5]
# [1,] 1.134243 1.340609 1.599390 1.062322 1.041968
# [2,] 1.134243 1.340601 1.599511 1.062315 1.041970
# [3,] 1.134239 1.340603 1.599247 1.061383 1.041958
# [4,] 1.134232 1.340612 1.599505 1.062321 1.041977
# [5,] 1.134243 1.340598 1.599411 1.062307 1.042002
# [6,] 1.134244 1.340597 1.599505 1.062349 1.042004
# [7,] 1.134264 1.340469 1.599616 1.062315 1.041968
# [8,] 1.134256 1.340630 1.598917 1.062228 1.041999
# [9,] 1.134186 1.340569 1.599432 1.062299 1.042053
#[10,] 1.134244 1.340683 1.599538 1.061689 1.041966
#[11,] 1.134242 1.340607 1.599454 1.062323 1.042027
#[12,] 1.134476 1.340607 1.599497 1.062322 1.041996
#
#$Rd
#          [,1]     [,2]     [,3]     [,4]     [,5]
# [1,] 1.144654 1.353647 1.632214 1.115466 1.136746
# [2,] 1.144654 1.353647 1.632214 1.115466 1.136746
# [3,] 1.144654 1.353647 1.632214 1.115466 1.136746
# [4,] 1.144654 1.353647 1.632214 1.115466 1.136746
# [5,] 1.144654 1.353647 1.632214 1.115466 1.136746
# [6,] 1.144654 1.353647 1.632214 1.115466 1.136746
# [7,] 1.144654 1.353647 1.632214 1.115466 1.136746
# [8,] 1.144654 1.353647 1.632214 1.115466 1.136746
# [9,] 1.144654 1.353647 1.632214 1.115466 1.136746
#[10,] 1.144654 1.353647 1.632214 1.115466 1.136746
#[11,] 1.144654 1.353647 1.632214 1.115466 1.136746
#[12,] 1.144654 1.353647 1.632214 1.115466 1.136746


#$beta1
#             [,1]        [,2]        [,3]       [,4]        [,5]
# [1,]  0.58604925  0.12301722  0.05325056 -0.8709711 -1.16108662
# [2,]  0.55307365  0.09347444  0.03083708 -0.8979184 -1.18208604
# [3,]  2.98746693  2.41281965  2.27033389  1.3289787  0.95855460
# [4,] -1.82610933 -2.18055913 -2.17284086 -3.0866180 -3.29670237
# [5,]  1.36654230  0.83078178  0.70666569 -0.2193170 -0.56607804
# [6,] -0.19872543 -0.59206396 -0.60816935 -1.5296845 -1.76511667
# [7,]  0.47827942  0.02047129 -0.03906356 -0.9688350 -1.25003125
# [8,]  0.69383246  0.22503623  0.14381369 -0.7736453 -1.07310560
# [9,]  2.93217740  2.05653791  1.74274221  0.6897412  0.09468177
#[10,] -1.95820408 -1.96450149 -1.74151897 -2.6141704 -2.63451360
#[11,]  1.19496596  0.64411014  0.51328451 -0.4175201 -0.77088457
#[12,] -0.02855136 -0.40683266 -0.41559972 -1.3330841 -1.56212372

#$beta2
#            [,1]        [,2]        [,3]        [,4]        [,5]
# [1,]  0.2206062  0.15447729  0.10583195  0.26495217  0.25709855
# [2,]  0.2164599  0.14897980  0.09823305  0.25892733  0.24974626
# [3,] -0.1318766 -0.15356304 -0.18039418 -0.02362266 -0.01396988
# [4,]  0.5799235  0.46847862  0.39522842  0.55756079  0.53114440
# [5,]  0.1144536  0.07749921  0.04522539  0.20024202  0.20483639
# [6,]  0.3294993  0.23448793  0.16905975  0.33154080  0.31114492
# [7,]  0.2113619  0.14310491  0.09145706  0.25260508  0.24275634
# [8,]  0.2298357  0.16603631  0.12077244  0.27744250  0.27166285
# [9,] -0.2249124 -0.14783929 -0.13147660  0.03038228  0.08361334
#[10,]  0.7952024  0.52701214  0.37906811  0.54335462  0.46831756
#[11,]  0.1139205  0.08297744  0.05270005  0.20827782  0.21430188
#[12,]  0.3311568  0.22970509  0.16184887  0.32386042  0.30194883


out.zelt4=example.zelterman(tab7.1)
out.zelt4

#$AIC
#          [,1]     [,2]
# [1,] 18.20969 23.45910
# [2,] 18.20876 23.61038
# [3,] 18.23977 23.78514
# [4,] 18.18389 23.15973
# [5,] 18.22995 23.57234
# [6,] 18.19149 23.37561
# [7,] 18.20824 23.57234
# [8,] 18.21117 23.37561
# [9,] 18.38019 24.98340
#[10,] 18.12212 22.37547
#[11,] 18.23701 23.61038
#[12,] 18.18595 23.34431

#$Test.Statistics
#            [,1]     [,2]
# [1,] 0.09761512 2.000000
# [2,] 0.09689905 2.116182
# [3,] 0.12098766 2.233129
# [4,] 0.07844120 1.794612
# [5,] 0.11324467 2.088089
# [6,] 0.08399027 1.949820
# [7,] 0.09649987 2.088089
# [8,] 0.09863173 1.949824
# [9,] 0.24308237 3.167092
#[10,] 0.03614894 1.291963
#[11,] 0.11881105 2.116189
#[12,] 0.07993661 1.929008

#$P.value
#           [,1]      [,2]
# [1,] 0.9998471 0.8491450
# [2,] 0.9998498 0.8328468
# [3,] 0.9997406 0.8160366
# [4,] 0.9999109 0.8767712
# [5,] 0.9997795 0.8368269
# [6,] 0.9998945 0.8560439
# [7,] 0.9998513 0.8368270
# [8,] 0.9998431 0.8560433
# [9,] 0.9985787 0.6742429
#[10,] 0.9999870 0.9357565
#[11,] 0.9997519 0.8328458
#[12,] 0.9999066 0.8588784

#$mu
#           [,1]     [,2]
# [1,] 0.7142880 3.000000
# [2,] 0.7142870 3.000060
# [3,] 0.7142934 3.000000
# [4,] 0.7142897 3.000056
# [5,] 0.7142773 2.999987
# [6,] 0.7142904 3.000001
# [7,] 0.7142866 3.000042
# [8,] 0.7142221 3.000027
# [9,] 0.7142844 3.000193
#[10,] 0.7142862 2.999999
#[11,] 0.7142786 3.000012
#[12,] 0.7142897 3.000003

#$sigma2
#           [,1]     [,2]
# [1,] 0.4897937 4.000000
# [2,] 0.4897948 3.999881
# [3,] 0.4897962 3.999789
# [4,] 0.4897860 3.999802
# [5,] 0.4898021 4.000320
# [6,] 0.4897952 4.000130
# [7,] 0.4897957 3.999791
# [8,] 0.4899654 3.999681
# [9,] 0.4897993 3.999058
#[10,] 0.4897953 3.999752
#[11,] 0.4897755 3.999899
#[12,] 0.4897945 3.999828

#$Dispersion
#           [,1]     [,2]
# [1,] 0.7783727 2.666667
# [2,] 0.7783753 2.666587
# [3,] 0.7783716 2.666526
# [4,] 0.7783589 2.666535
# [5,] 0.7783961 2.666880
# [6,] 0.7783728 2.666753
# [7,] 0.7783771 2.666528
# [8,] 0.7787076 2.666454
# [9,] 0.7783849 2.666038
#[10,] 0.7783770 2.666501
#[11,] 0.7783526 2.666599
#[12,] 0.7783725 2.666552

#$beta1
#             [,1]        [,2]
# [1,] -0.77746420  0.00000000
# [2,] -0.82308899  0.06728445
# [3,]  1.41530100  2.01660731
# [4,] -2.97407737 -2.01642517
# [5,]  0.03464193  0.58578094
# [6,] -1.59133846 -0.58745400
# [7,] -0.91562699 -0.02424090
# [8,] -0.63806938  0.02336125
# [9,]  1.51603077  1.60115931
#[10,] -3.14270034 -1.59741917
#[11,] -0.12838516  0.38481299
#[12,] -1.42976249 -0.38757650

#$beta2
#           [,1]          [,2]
# [1,] 0.7040806  7.180547e-15
# [2,] 0.6963496 -3.767577e-02
# [3,] 0.2916423 -3.361012e-01
# [4,] 1.1192170  3.360685e-01
# [5,] 0.5773919 -4.681756e-02
# [6,] 0.8320360  4.698828e-02
# [7,] 0.6886972 -4.677311e-02
# [8,] 0.7188748  4.702556e-02
# [9,] 0.1434515 -2.668679e-01
#[10,] 1.3154188  2.662366e-01
#[11,] 0.5711815 -3.767680e-02
#[12,] 0.8392506  3.797387e-02

outcum.zelt4=cum.example.zelterman(tab7.1)
outcum.zelt4

#$AIC
#          [,1]     [,2]     [,3]     [,4]     [,5]
# [1,] 57.98424 96.38004 159.9677 175.5814 215.4275
# [2,] 57.98497 96.37993 159.9732 175.5868 216.0936
# [3,] 57.98447 96.38391 159.7084 175.3420 214.4540
# [4,] 57.98487 96.37697 160.3102 175.9073 216.6978
# [5,] 57.98417 96.38261 159.7868 175.4129 214.2157
# [6,] 57.98488 96.37785 160.1916 175.7940 216.8981
# [7,] 57.98489 96.37986 159.9776 175.5910 216.0892
# [8,] 57.98417 96.38022 159.9586 175.5725 214.8188
# [9,] 57.98429 96.40611 159.4979 175.2379 213.2904
#[10,] 57.98489 96.37061 161.9786 177.5408 222.5633
#[11,] 57.98418 96.38355 159.7289 175.3601 214.0139
#[12,] 57.98488 96.37720 160.2776 175.8762 217.2103

#$Test.Statistics
#              [,1]        [,2]     [,3]     [,4]      [,5]
# [1,] 3.256171e-05 0.007432135 2.064911 1.648025  5.657344
# [2,] 3.970618e-04 0.007371015 2.069734 1.654363  7.127987
# [3,] 1.477130e-04 0.009511985 1.784059 1.372811  4.259511
# [4,] 3.483536e-04 0.005807076 2.436194 2.014538  7.510986
# [5,] 5.890744e-08 0.008807900 1.871169 1.455851  3.783919
# [6,] 3.544463e-04 0.006271018 2.306205 1.887229  8.406479
# [7,] 3.557403e-04 0.007338584 2.074312 1.659364  7.092161
# [8,] 5.011667e-09 0.007526809 2.056878 1.641562  4.579169
# [9,] 5.869434e-05 0.022023771 1.528614 1.142920  2.324811
#[10,] 3.567851e-04 0.002509311 4.381677 3.957669 17.276452
#[11,] 7.693053e-07 0.009312696 1.808024 1.394375  3.519947
#[12,] 3.511505e-04 0.005931107 2.400487 1.978418  8.954340

#$P.value
#      [,1]      [,2]      [,3]      [,4]        [,5]
# [1,]    1 0.9999997 0.8400924 0.8953816 0.341000805
# [2,]    1 0.9999998 0.8394143 0.8945983 0.211294418
# [3,]    1 0.9999995 0.8781442 0.9272642 0.512689330
# [4,]    1 0.9999999 0.7860723 0.8471298 0.185325282
# [5,]    1 0.9999996 0.8666686 0.9181077 0.580926589
# [6,]    1 0.9999998 0.8053551 0.8645187 0.135211015
# [7,]    1 0.9999998 0.8387699 0.8939789 0.213875908
# [8,]    1 0.9999997 0.8412203 0.8961782 0.469361228
# [9,]    1 0.9999962 0.9097444 0.9502347 0.802616008
#[10,]    1 1.0000000 0.4958689 0.5555266 0.004004270
#[11,]    1 0.9999996 0.8750193 0.9249276 0.620371532
#[12,]    1 0.9999999 0.7914016 0.8521232 0.110899297

#$mu
#           [,1]      [,2]      [,3]      [,4]      [,5]
# [1,] 0.2499996 0.2957746 0.4545437 0.4737305 0.6000074
# [2,] 0.2499729 0.2957753 0.4545497 0.4736896 0.6000032
# [3,] 0.2499967 0.2957746 0.4545469 0.4736847 0.6000145
# [4,] 0.2499785 0.2957745 0.4545378 0.4736583 0.6000110
# [5,] 0.2500028 0.2957747 0.4545613 0.4737009 0.6000062
# [6,] 0.2500066 0.2957733 0.4545468 0.4736846 0.6000045
# [7,] 0.2499780 0.2957746 0.4545452 0.4736840 0.6000045
# [8,] 0.2499999 0.2957747 0.4545452 0.4735473 0.6000170
# [9,] 0.2499986 0.2957760 0.4545463 0.4736851 0.6000161
#[10,] 0.2499780 0.2957746 0.4545479 0.4736844 0.5999995
#[11,] 0.2500000 0.2957746 0.4544942 0.4737003 0.5999397
#[12,] 0.2499783 0.2957737 0.4544842 0.4737287 0.6001217

#$sigma2
#           [,1]      [,2]      [,3]      [,4]      [,5]
# [1,] 0.1875012 0.2364609 0.4752134 0.4809492 0.9599829
# [2,] 0.1875030 0.2364629 0.4752148 0.4809147 0.9599745
# [3,] 0.1875045 0.2364609 0.4752061 0.4808848 0.9600011
# [4,] 0.1875038 0.2364608 0.4752006 0.4808677 0.9599962
# [5,] 0.1875014 0.2364610 0.4752148 0.4808951 0.9599750
# [6,] 0.1875181 0.2364607 0.4752066 0.4808845 0.9599826
# [7,] 0.1875038 0.2364609 0.4752041 0.4808840 0.9599821
# [8,] 0.1874999 0.2364610 0.4752041 0.4805236 0.9600077
# [9,] 0.1875017 0.2364563 0.4752013 0.4808876 0.9599983
#[10,] 0.1875038 0.2364609 0.4752091 0.4808854 0.9600082
#[11,] 0.1875000 0.2364609 0.4751631 0.4809339 0.9601147
#[12,] 0.1875038 0.2364607 0.4751862 0.4809801 0.9601598

#$Dispersion
#           [,1]      [,2]     [,3]     [,4]     [,5]
# [1,] 0.7826147 0.8409167 1.131168 1.102268 1.777727
# [2,] 0.7827023 0.8409221 1.131157 1.102276 1.777722
# [3,] 0.7826374 0.8409167 1.131143 1.102217 1.777742
# [4,] 0.7826887 0.8409167 1.131151 1.102234 1.777742
# [5,] 0.7826061 0.8409167 1.131131 1.102206 1.777715
# [6,] 0.7826643 0.8409195 1.131144 1.102217 1.777734
# [7,] 0.7826904 0.8409167 1.131142 1.102217 1.777733
# [8,] 0.7826088 0.8409169 1.131142 1.101682 1.777747
# [9,] 0.7826202 0.8408968 1.131133 1.102223 1.777732
#[10,] 0.7826907 0.8409168 1.131148 1.102219 1.777794
#[11,] 0.7826089 0.8409168 1.131161 1.102297 1.778149
#[12,] 0.7826894 0.8409186 1.131239 1.102342 1.777753

#$Rd
#         [,1]      [,2]     [,3]     [,4]     [,5]
# [1,] 0.79926 0.8529302 1.144149 1.113948 1.795735
# [2,] 0.79926 0.8529302 1.144149 1.113948 1.795735
# [3,] 0.79926 0.8529302 1.144149 1.113948 1.795735
# [4,] 0.79926 0.8529302 1.144149 1.113948 1.795735
# [5,] 0.79926 0.8529302 1.144149 1.113948 1.795735
# [6,] 0.79926 0.8529302 1.144149 1.113948 1.795735
# [7,] 0.79926 0.8529302 1.144149 1.113948 1.795735
# [8,] 0.79926 0.8529302 1.144149 1.113948 1.795735
# [9,] 0.79926 0.8529302 1.144149 1.113948 1.795735
#[10,] 0.79926 0.8529302 1.144149 1.113948 1.795735
#[11,] 0.79926 0.8529302 1.144149 1.113948 1.795735
#[12,] 0.79926 0.8529302 1.144149 1.113948 1.795735

#$beta1
#           [,1]         [,2]         [,3]       [,4]       [,5]
# [1,] -4.760716 -0.004111614  0.620137873  0.5147770  1.0220229
# [2,] -3.557236 -0.050824939  0.576513991  0.4709763  1.0015897
# [3,] -1.775139  2.221979473  2.779529286  2.6778562  3.1021004
# [4,] -5.805263 -2.230884468 -1.545445444 -1.6540471 -1.0566230
# [5,] -7.099127  0.830615498  1.407713458  1.3052483  1.7315242
# [6,] -4.403951 -0.839141864 -0.171414084 -0.2792349  0.3100807
# [7,] -3.705345 -0.144021521  0.484804878  0.3790670  0.9125054
# [8,] -9.010535  0.135797576  0.755165896  0.6495645  1.1293453
# [9,] -1.981204  2.458337918  2.761982915  2.6716278  2.8951260
#[10,] -6.048720 -2.480586662 -1.621107988 -1.7393126 -0.8262847
#[11,] -5.950884  0.675348051  1.238237520  1.1364340  1.5474496
#[12,] -4.254311 -0.684126584 -0.003303165 -0.1119145  0.4941999

#$beta2
#           [,1]      [,2]       [,3]       [,4]        [,5]
# [1,]  5.859337 0.9902609  0.2400165  0.2735798 -0.01542613
# [2,]  4.602325 0.9832417  0.2309682  0.2647061 -0.03454985
# [3,]  4.665557 0.5556076 -0.1493377 -0.1181431 -0.36630171
# [4,]  5.112298 1.4255325  0.6336413  0.6694386  0.33572184
# [5,]  8.890871 0.8484665  0.1300934  0.1618337 -0.08866799
# [6,]  4.809450 1.1323355  0.3525504  0.3877905  0.05944858
# [7,]  4.649992 0.9760346  0.2227611  0.2566170 -0.04405474
# [8,] 10.263299 1.0044871  0.2574420  0.2912793  0.01419917
# [9,]  4.871594 0.3174511 -0.2150562 -0.1914662 -0.34217223
#[10,]  5.355759 1.6757785  0.7653362  0.8075402  0.31562779
#[11,]  7.576086 0.8370990  0.1285363  0.1597565 -0.08315513
#[12,]  4.826518 1.1439303  0.3553068  0.3909365  0.05420445


out.zelt5=example.zelterman(tab8.2)
out.zelt5

#$AIC
#          [,1]     [,2]     [,3]     [,4]     [,5]
# [1,] 207.9221 281.7527 158.6891 70.40504 44.17325
# [2,] 207.9181 281.7245 158.7274 70.38445 44.17102
# [3,] 208.2204 282.3850 157.8071 71.10839 44.21134
# [4,] 207.7969 281.8691 160.1701 69.89362 44.24614
# [5,] 208.0875 282.0801 158.0696 70.83715 44.19593
# [6,] 207.8259 281.7397 159.5955 70.05369 44.20101
# [7,] 207.9151 281.7088 158.7592 70.37116 44.17173
# [8,] 207.9292 281.8004 158.6336 70.44018 44.17683
# [9,] 211.2224 289.5937 158.8368 74.44250 44.99733
#[10,] 208.5224 286.8348 167.5802 69.31670 44.83672
#[11,] 208.1411 282.1246 157.9591 70.93201 44.21171
#[12,] 207.8160 281.8399 159.8013 69.99106 44.19974

#$Test.Statistics
#           [,1]     [,2]      [,3]      [,4]      [,5]
# [1,] 0.3147530 1.795445  4.631586 1.4215201 0.9545620
# [2,] 0.3120557 1.773761  4.692900 1.4044644 0.9534887
# [3,] 0.5310232 2.267535  3.177383 2.0239358 0.9450302
# [4,] 0.2573400 2.080765  7.278837 1.0064429 1.0830028
# [5,] 0.4297079 2.018425  3.628448 1.7884999 0.9473473
# [6,] 0.2631638 1.890514  6.207484 1.1332142 1.0152792
# [7,] 0.3100179 1.762389  4.754134 1.3936526 0.9555241
# [8,] 0.3195577 1.832972  4.546970 1.4509221 0.9561244
# [9,] 3.4693850 9.454096  3.874438 5.2578000 1.7549908
#[10,] 1.1462060 7.736491 25.376262 0.6005265 1.8558298
#[11,] 0.4698456 2.046272  3.475463 1.8721103 0.9594475
#[12,] 0.2636274 2.007549  6.518023 1.0801571 1.0180505

#$P.value
#           [,1]      [,2]         [,3]      [,4]      [,5]
# [1,] 0.9998823 0.9702901 0.7048196842 0.9849030 0.9955286
# [2,] 0.9998857 0.9712954 0.6973814255 0.9854341 0.9955443
# [3,] 0.9993248 0.9435625 0.8681222196 0.9585037 0.9956672
# [4,] 0.9999405 0.9552352 0.4004363067 0.9947239 0.9933762
# [5,] 0.9996653 0.9588136 0.8214391221 0.9706142 0.9956338
# [6,] 0.9999358 0.9656507 0.5157423637 0.9923842 0.9945783
# [7,] 0.9998882 0.9718148 0.6899381218 0.9857649 0.9955144
# [8,] 0.9998762 0.9685039 0.7150541909 0.9839603 0.9955056
# [9,] 0.8384575 0.2216760 0.7941208588 0.6285349 0.9721498
#[10,] 0.9921132 0.3564024 0.0006506151 0.9989889 0.9673871
#[11,] 0.9995496 0.9572349 0.8378174916 0.9665783 0.9954564
#[12,] 0.9999355 0.9594214 0.4807232247 0.9934298 0.9945321

#$mu
#           [,1]      [,2]     [,3]     [,4]     [,5]
# [1,] 0.4909105 0.9903955 1.340166 1.809531 2.333336
# [2,] 0.4909090 0.9903850 1.340033 1.809529 2.333336
# [3,] 0.4909087 0.9903892 1.339968 1.809798 2.333342
# [4,] 0.4909436 0.9903872 1.340031 1.809526 2.333342
# [5,] 0.4909097 0.9903893 1.339983 1.809520 2.333331
# [6,] 0.4909101 0.9904415 1.339999 1.809523 2.333343
# [7,] 0.4909026 0.9905158 1.339992 1.809527 2.333334
# [8,] 0.4909295 0.9902838 1.340057 1.809550 2.333346
# [9,] 0.4909059 0.9904031 1.339982 1.809539 2.333309
#[10,] 0.4909100 0.9904068 1.340006 1.809738 2.333358
#[11,] 0.4909110 0.9903997 1.340006 1.809532 2.333329
#[12,] 0.4908369 0.9903848 1.340000 1.809531 2.333342

#$sigma2
#           [,1]     [,2]     [,3]     [,4]     [,5]
# [1,] 0.5226437 1.163361 1.784376 1.487496 1.722143
# [2,] 0.5226421 1.163343 1.784370 1.487485 1.722146
# [3,] 0.5226402 1.163349 1.784006 1.486336 1.722129
# [4,] 0.5226714 1.163354 1.784361 1.487492 1.722160
# [5,] 0.5226427 1.163344 1.784205 1.487468 1.722197
# [6,] 0.5226438 1.163396 1.784327 1.487529 1.722204
# [7,] 0.5226454 1.163360 1.784444 1.487484 1.722142
# [8,] 0.5226688 1.163266 1.783732 1.487375 1.722198
# [9,] 0.5226130 1.163333 1.784227 1.487468 1.722272
#[10,] 0.5226439 1.163436 1.784370 1.486731 1.722148
#[11,] 0.5226438 1.163363 1.784278 1.487498 1.722236
#[12,] 0.5226779 1.163348 1.784319 1.487496 1.722192

#$Dispersion
#          [,1]     [,2]     [,3]     [,4]     [,5]
# [1,] 1.134243 1.340609 1.599390 1.062322 1.041968
# [2,] 1.134243 1.340601 1.599511 1.062315 1.041970
# [3,] 1.134239 1.340603 1.599247 1.061383 1.041958
# [4,] 1.134232 1.340612 1.599505 1.062321 1.041977
# [5,] 1.134243 1.340598 1.599411 1.062307 1.042002
# [6,] 1.134244 1.340597 1.599505 1.062349 1.042004
# [7,] 1.134264 1.340469 1.599616 1.062315 1.041968
# [8,] 1.134256 1.340630 1.598917 1.062228 1.041999
# [9,] 1.134186 1.340569 1.599432 1.062299 1.042053
#[10,] 1.134244 1.340683 1.599538 1.061689 1.041966
#[11,] 1.134242 1.340607 1.599454 1.062323 1.042027
#[12,] 1.134476 1.340607 1.599497 1.062322 1.041996

#$Rd
#          [,1]     [,2]     [,3]     [,4]     [,5]
# [1,] 1.144654 1.353647 1.632214 1.115466 1.136746
# [2,] 1.144654 1.353647 1.632214 1.115466 1.136746
# [3,] 1.144654 1.353647 1.632214 1.115466 1.136746
# [4,] 1.144654 1.353647 1.632214 1.115466 1.136746
# [5,] 1.144654 1.353647 1.632214 1.115466 1.136746
# [6,] 1.144654 1.353647 1.632214 1.115466 1.136746
# [7,] 1.144654 1.353647 1.632214 1.115466 1.136746
# [8,] 1.144654 1.353647 1.632214 1.115466 1.136746
# [9,] 1.144654 1.353647 1.632214 1.115466 1.136746
#[10,] 1.144654 1.353647 1.632214 1.115466 1.136746
#[11,] 1.144654 1.353647 1.632214 1.115466 1.136746
#[12,] 1.144654 1.353647 1.632214 1.115466 1.136746


#$beta1
#             [,1]        [,2]        [,3]       [,4]        [,5]
# [1,]  0.58604925  0.12301722  0.05325056 -0.8709711 -1.16108662
# [2,]  0.55307365  0.09347444  0.03083708 -0.8979184 -1.18208604
# [3,]  2.98746693  2.41281965  2.27033389  1.3289787  0.95855460
# [4,] -1.82610933 -2.18055913 -2.17284086 -3.0866180 -3.29670237
# [5,]  1.36654230  0.83078178  0.70666569 -0.2193170 -0.56607804
# [6,] -0.19872543 -0.59206396 -0.60816935 -1.5296845 -1.76511667
# [7,]  0.47827942  0.02047129 -0.03906356 -0.9688350 -1.25003125
# [8,]  0.69383246  0.22503623  0.14381369 -0.7736453 -1.07310560
# [9,]  2.93217740  2.05653791  1.74274221  0.6897412  0.09468177
#[10,] -1.95820408 -1.96450149 -1.74151897 -2.6141704 -2.63451360
#[11,]  1.19496596  0.64411014  0.51328451 -0.4175201 -0.77088457
#[12,] -0.02855136 -0.40683266 -0.41559972 -1.3330841 -1.56212372

#$beta2
#            [,1]        [,2]        [,3]        [,4]        [,5]
# [1,]  0.2206062  0.15447729  0.10583195  0.26495217  0.25709855
# [2,]  0.2164599  0.14897980  0.09823305  0.25892733  0.24974626
# [3,] -0.1318766 -0.15356304 -0.18039418 -0.02362266 -0.01396988
# [4,]  0.5799235  0.46847862  0.39522842  0.55756079  0.53114440
# [5,]  0.1144536  0.07749921  0.04522539  0.20024202  0.20483639
# [6,]  0.3294993  0.23448793  0.16905975  0.33154080  0.31114492
# [7,]  0.2113619  0.14310491  0.09145706  0.25260508  0.24275634
# [8,]  0.2298357  0.16603631  0.12077244  0.27744250  0.27166285
# [9,] -0.2249124 -0.14783929 -0.13147660  0.03038228  0.08361334
#[10,]  0.7952024  0.52701214  0.37906811  0.54335462  0.46831756
#[11,]  0.1139205  0.08297744  0.05270005  0.20827782  0.21430188
#[12,]  0.3311568  0.22970509  0.16184887  0.32386042  0.30194883



out.zelt6=example.zelterman(tab9.1.2)
out.zelt6

#$AIC
#          [,1]     [,2]     [,3]
# [1,] 314.0457 75.98499 12.47407
# [2,] 314.0377 75.99217 12.47395
# [3,] 314.2285 75.56162 12.50741
# [4,] 314.8903 76.58985 12.44562
# [5,] 313.9815 75.71971 12.49402
# [6,] 314.4645 76.31256 12.45598
# [7,] 314.0313 75.99899 12.47381
# [8,] 314.0609 75.97111 12.47434
# [9,] 320.0611 75.61470 12.72619
#[10,] 325.1774 79.28284 12.36531
#[11,] 313.8845 75.69347 12.50095
#[12,] 314.7064 76.36399 12.45043
#
#$Test.Statistics
#           [,1]     [,2]       [,3]
# [1,]  3.249940 2.391212 0.10344550
# [2,]  3.245348 2.402660 0.10334440
# [3,]  3.258855 1.733506 0.12968406
# [4,]  4.279795 3.400936 0.08199505
# [5,]  3.085968 1.977661 0.11900404
# [6,]  3.772820 2.921741 0.08970216
# [7,]  3.242083 2.413702 0.10324936
# [8,]  3.258412 2.368646 0.10365640
# [9,]  9.228308 1.551125 0.32629513
#[10,] 16.406329 8.530399 0.02699135
#[11,]  2.984053 1.938513 0.12448858
#[12,]  4.030513 2.997001 0.08554950

#$P.value
#            [,1]      [,2]      [,3]
# [1,] 0.95355589 0.9836674 1.0000000
# [2,] 0.95376821 0.9833877 1.0000000
# [3,] 0.95314204 0.9950156 0.9999999
# [4,] 0.89205115 0.9462606 1.0000000
# [5,] 0.96078878 0.9918159 0.9999999
# [6,] 0.92571866 0.9673169 1.0000000
# [7,] 0.95391882 0.9831149 1.0000000
# [8,] 0.95316266 0.9842099 1.0000000
# [9,] 0.41647202 0.9967494 0.9999952
#[10,] 0.05886616 0.4817001 1.0000000
#[11,] 0.96492385 0.9924030 0.9999999
#[12,] 0.90939274 0.9644138 1.0000000

#$mu
#          [,1]     [,2]      [,3]
# [1,] 1.848495 2.125002 1.0000084
# [2,] 1.848495 2.124997 1.0000114
# [3,] 1.848484 2.125000 1.0000150
# [4,] 1.848488 2.125000 0.9999965
# [5,] 1.848485 2.125000 1.0000010
# [6,] 1.848494 2.125009 1.0000001
# [7,] 1.848477 2.125010 1.0000050
# [8,] 1.848489 2.125003 0.9999995
# [9,] 1.848482 2.124997 1.0000002
#[10,] 1.848486 2.125000 0.9999917
#[11,] 1.848547 2.125003 1.0000001
#[12,] 1.848488 2.125004 0.9999757

#$sigma2
#          [,1]     [,2]      [,3]
# [1,] 1.421447 1.192665 0.5000044
# [2,] 1.421460 1.192689 0.5000060
# [3,] 1.421431 1.192673 0.4999928
# [4,] 1.421458 1.192685 0.4999967
# [5,] 1.421439 1.192677 0.4999983
# [6,] 1.421462 1.192645 0.4999982
# [7,] 1.421423 1.192681 0.4999812
# [8,] 1.421481 1.192717 0.5000010
# [9,] 1.421370 1.192675 0.4999965
#[10,] 1.421466 1.192691 0.5000207
#[11,] 1.421206 1.192676 0.4999988
#[12,] 1.421451 1.192683 0.5000539

#$Dispersion
#           [,1]      [,2]      [,3]
# [1,] 0.9433536 0.7127029 0.5555564
# [2,] 0.9433621 0.7127185 0.5555566
# [3,] 0.9433473 0.7127086 0.5555402
# [4,] 0.9433637 0.7127154 0.5555536
# [5,] 0.9433526 0.7127109 0.5555531
# [6,] 0.9433641 0.7126898 0.5555535
# [7,] 0.9433448 0.7127110 0.5555322
# [8,] 0.9433786 0.7127340 0.5555570
# [9,] 0.9433081 0.7127100 0.5555516
#[10,] 0.9433704 0.7127194 0.5555827
#[11,] 0.9431734 0.7127096 0.5555542
#[12,] 0.9433589 0.7127133 0.5556275
#
#$Rd
#           [,1]      [,2]      [,3]
# [1,] 0.9530111 0.7437178 0.7407407
# [2,] 0.9530111 0.7437178 0.7407407
# [3,] 0.9530111 0.7437178 0.7407407
# [4,] 0.9530111 0.7437178 0.7407407
# [5,] 0.9530111 0.7437178 0.7407407
# [6,] 0.9530111 0.7437178 0.7407407
# [7,] 0.9530111 0.7437178 0.7407407
# [8,] 0.9530111 0.7437178 0.7407407
# [9,] 0.9530111 0.7437178 0.7407407
#[10,] 0.9530111 0.7437178 0.7407407
#[11,] 0.9530111 0.7437178 0.7407407
#[12,] 0.9530111 0.7437178 0.7407407

#$beta1
#            [,1]       [,2]       [,3]
# [1,] -0.9975660 -1.6531200 -1.7227655
# [2,] -1.0214330 -1.6767605 -1.7490063
# [3,]  1.3886578  0.6919062  0.9236744
# [4,] -3.4042785 -4.0203681 -4.3747252
# [5,] -0.3463183 -1.0262189 -0.9191871
# [6,] -1.6557701 -2.2876263 -2.5283920
# [7,] -1.0827661 -1.7378677 -1.8126986
# [8,] -0.9123620 -1.5682862 -1.6329392
# [9,]  0.6375204 -0.2266688  1.0998316
#[10,] -2.8969233 -3.3967566 -4.6989522
#[11,] -0.5455988 -1.2316066 -1.0848569
#[12,] -1.4585465 -2.0845838 -2.3637782

#$beta2
#            [,1]      [,2]      [,3]
# [1,] 0.28933543 0.3940222 0.8922487
# [2,] 0.28638084 0.3910629 0.8901262
# [3,] 0.02217972 0.1327419 0.5261678
# [4,] 0.56186918 0.6602463 1.2616865
# [5,] 0.22418153 0.3324808 0.7684284
# [6,] 0.35630240 0.4572376 1.0173323
# [7,] 0.28239599 0.3870623 0.8870011
# [8,] 0.29627947 0.4009618 0.8975568
# [9,] 0.05252771 0.1876973 0.2263254
#[10,] 0.58817766 0.6638695 1.6519409
#[11,] 0.23238536 0.3417894 0.7630604
#[12,] 0.34854596 0.4483653 1.0234339


out.zelt7=example.zelterman(tab9.1.3)
out.zelt7

#$AIC
#          [,1]     [,2]     [,3]
# [1,] 168.6205 27.45598 4.022225
# [2,] 168.6716 27.46525 4.017892
# [3,] 167.7907 27.24119 4.011769
# [4,] 170.0310 27.71161 4.025957
# [5,] 168.0666 27.33195 4.028108
# [6,] 169.3864 27.59425 4.026309
# [7,] 168.7137 27.47333 4.013674
# [8,] 168.5337 27.43899 4.012307
# [9,] 168.2843 27.14469 4.013484
#[10,] 174.6297 28.12289 4.013979
#[11,] 168.0932 27.34538 4.016725
#[12,] 169.3801 27.58135 4.017656

#$Test.Statistics
#           [,1]     [,2]        [,3]
# [1,]  8.240286 3.711729 0.011174415
# [2,]  8.333392 3.718923 0.008986362
# [3,]  6.673795 3.436465 0.005901916
# [4,] 10.849488 4.074461 0.013062940
# [5,]  7.225647 3.556639 0.014153062
# [6,]  9.648861 3.897762 0.013241410
# [7,]  8.414839 3.726613 0.006860510
# [8,]  8.082858 3.697954 0.006172477
# [9,]  6.798786 3.293339 0.006764766
#[10,] 19.352523 4.841991 0.007013996
#[11,]  7.268361 3.566469 0.008397653
#[12,]  9.605200 3.887846 0.008867220

#$P.value
#            [,1]      [,2] [,3]
# [1,] 0.51012442 0.9293414    1
# [2,] 0.50092871 0.9289198    1
# [3,] 0.67103979 0.9444556    1
# [4,] 0.28616708 0.9064436    1
# [5,] 0.61363915 0.9381024    1
# [6,] 0.37965228 0.9180117    1
# [7,] 0.49294027 0.9284677    1
# [8,] 0.52581823 0.9301450    1
# [9,] 0.65805957 0.9515214    1
#[10,] 0.02235684 0.8478589    1
#[11,] 0.60920119 0.9375658    1
#[12,] 0.38338099 0.9186378    1

#$mu
#          [,1]     [,2]     [,3]
# [1,] 3.155573 3.428579 6.998799
# [2,] 3.155568 3.428575 6.999483
# [3,] 3.155565 3.428552 6.999683
# [4,] 3.155573 3.428572 6.999667
# [5,] 3.155566 3.428584 6.998525
# [6,] 3.155548 3.428572 6.998582
# [7,] 3.155557 3.428572 6.999627
# [8,] 3.155545 3.428573 6.999661
# [9,] 3.155563 3.428582 6.999630
#[10,] 3.155581 3.428591 6.999625
#[11,] 3.155552 3.428586 6.999723
#[12,] 3.155556 3.428574 6.999663

#$sigma2
#          [,1]     [,2]        [,3]
# [1,] 2.308894 1.673324 0.011049492
# [2,] 2.308874 1.673329 0.008906061
# [3,] 2.308881 1.673502 0.005867187
# [4,] 2.308968 1.673361 0.012894400
# [5,] 2.308889 1.673312 0.013953388
# [6,] 2.308907 1.673352 0.013066368
# [7,] 2.308881 1.673332 0.006813626
# [8,] 2.309073 1.673333 0.006134497
# [9,] 2.309049 1.673333 0.006719175
#[10,] 2.308998 1.673378 0.006965003
#[11,] 2.308851 1.673312 0.008327645
#[12,] 2.308911 1.673337 0.008789172

#$Dispersion
#          [,1]      [,2]        [,3]
# [1,] 1.069027 0.7426883 0.005260460
# [2,] 1.069019 0.7426908 0.004240564
# [3,] 1.069022 0.7427702 0.002793730
# [4,] 1.069061 0.7427054 0.006139802
# [5,] 1.069026 0.7426825 0.006642605
# [6,] 1.069038 0.7427015 0.006220401
# [7,] 1.069024 0.7426924 0.003244353
# [8,] 1.069115 0.7426930 0.002921000
# [9,] 1.069100 0.7426922 0.003199382
#[10,] 1.069074 0.7427111 0.003316431
#[11,] 1.069011 0.7426822 0.003965336
#[12,] 1.069038 0.7426947 0.004185051

#$Rd
#          [,1]      [,2] [,3]
# [1,] 1.093441 0.8665459  NaN
# [2,] 1.093441 0.8665459  NaN
# [3,] 1.093441 0.8665459  NaN
# [4,] 1.093441 0.8665459  NaN
# [5,] 1.093441 0.8665459  NaN
# [6,] 1.093441 0.8665459  NaN
# [7,] 1.093441 0.8665459  NaN
# [8,] 1.093441 0.8665459  NaN
# [9,] 1.093441 0.8665459  NaN
#[10,] 1.093441 0.8665459  NaN
#[11,] 1.093441 0.8665459  NaN
#[12,] 1.093441 0.8665459  NaN


#$beta1
#            [,1]        [,2]      [,3]
# [1,] -1.2628964 -2.02406120 -72.59591
# [2,] -1.2780746 -2.03910931 -75.47366
# [3,]  0.8862026  0.05944217 -79.34687
# [4,] -3.4415124 -4.13661836 -72.59286
# [5,] -0.7535929 -1.55016146 -69.05026
# [6,] -1.7844810 -2.50887831 -70.45561
# [7,] -1.3329745 -2.09331381 -79.22155
# [8,] -1.1934582 -1.95514116 -81.09369
# [9,] -0.3501909 -1.28324782 -78.89744
#[10,] -2.4868145 -3.04833186 -79.82550
#[11,] -0.9591385 -1.75920788 -76.50725
#[12,] -1.5800417 -2.30085557 -75.92001

#$beta2
#             [,1]       [,2]     [,3]
# [1,]  0.20277231 0.29555972 5.193215
# [2,]  0.19813472 0.29105170 5.387346
# [3,] -0.01669337 0.08338772 5.616480
# [4,]  0.42630122 0.51154595 5.242063
# [5,]  0.16642807 0.26295032 4.948714
# [6,]  0.24087298 0.32962739 5.031348
# [7,]  0.19287736 0.28579164 5.644351
# [8,]  0.21278274 0.30538144 5.814598
# [9,]  0.10265283 0.21539095 5.625696
#[10,]  0.34345947 0.41196624 5.719444
#[11,]  0.17608433 0.27318162 5.470108
#[12,]  0.23130543 0.31949987 5.422676


#example.zelterman(tab2.cb)
#Error in optim(para0, nlk.semi.altham, gr = NULL, method = "BFGS", hessian = TRUE,  : 
#  non-finite finite-difference value [2]



#example.zelterman(tab3.cb)
#Error in optim(para0, nlk.semi.altham, gr = NULL, method = "BFGS", hessian = TRUE,  : 
#  non-finite finite-difference value [2]


# example.zelterman(tab6.cb)
# Error in solve.default(soln$hessian) : 
#  Lapack routine dgesv: system is exactly singular

# example.zelterman(tab8.cb)
# Error in solve.default(soln$hessian) : 
#  Lapack routine dgesv: system is exactly singular

out.zelt8=example.zelterman(tab.pala1)
out.zelt8

#$AIC
#          [,1]     [,2]     [,3]     [,4]     [,5]     [,6]     [,7]
# [1,] 9106.607 11689.42 15116.70 13034.68 18106.41 22279.65 24652.97
# [2,] 9106.622 11689.46 15116.84 13035.12 18108.01 22283.72 24661.10
# [3,] 9104.302 11683.50 15099.92 12985.27 17948.86 21917.64 24013.39
# [4,] 9109.004 11695.59 15134.63 13086.07 18272.87 22660.63 25323.42
# [5,] 9105.154 11685.69 15106.14 13003.87 18008.12 22054.40 24255.55
# [6,] 9108.096 11693.25 15127.70 13066.25 18208.11 22512.16 25062.16
# [7,] 9106.636 11689.49 15116.97 13035.53 18109.45 22287.31 24668.07
# [8,] 9106.579 11689.35 15116.43 13033.82 18103.38 22272.00 24637.89
# [9,] 9096.785 11664.45 15064.08 12829.40 17491.55 20841.78 22094.30
#[10,] 9119.672 11723.15 15213.04 13298.00 18961.72 24199.63 27960.93
#[11,] 9104.686 11684.48 15103.06 12995.13 17981.67 21995.99 24156.39
#[12,] 9108.591 11694.53 15131.09 13075.55 18236.95 22575.55 25169.21

#$Test.Statistics
#           [,1]      [,2]         [,3]         [,4]         [,5]         [,6]
# [1,]  34.11129  87.83451    234.15799    4822.9888 4.977348e+03 8.447529e+04
# [2,]  34.16712  87.94416    236.98217    4901.5784 5.113515e+03 9.243141e+04
# [3,]  26.17698  69.36482     94.39940    1236.1385 1.395057e+03 7.259842e+03
# [4,]  44.53417 111.83283    857.97641   22379.3375 3.350268e+04 2.057600e+06
# [5,]  28.90724  75.71721    127.76219    2040.7075 2.111825e+03 1.599682e+04
# [6,]  40.25788 102.07916    497.98389   12120.3071 1.503442e+04 5.928528e+05
# [7,]  34.21997  88.07098    239.75789    4993.5441 5.247815e+03 9.903348e+04
# [8,]  33.67365  87.56417    228.25338    4659.5698 4.725501e+03 7.252827e+04
# [9,]  10.45356  35.24656     37.20595     118.0648 2.659916e+02 7.089859e+02
#[10,] 136.15265 322.97264 318185.98685 8162423.7035 1.191052e+08 1.374282e+11
#[11,]  27.37540  72.14291    112.10508    1682.5288 1.816976e+03 1.254505e+04
#[12,]  42.51066 107.29798    612.97668   15343.7286 1.940310e+04 8.411521e+05
#              [,7]
# [1,] 6.772868e+05
# [2,] 7.537960e+05
# [3,] 5.324341e+04
# [4,] 1.203990e+07
# [5,] 1.311022e+05
# [6,] 4.023547e+06
# [7,] 8.055090e+05
# [8,] 5.827120e+05
# [9,] 1.760042e+03
#[10,] 1.013727e+11
#[11,] 1.041452e+05
#[12,] 5.290103e+06

#$P.value
#              [,1]         [,2]         [,3] [,4] [,5] [,6] [,7]
# [1,] 3.878292e-05 1.332268e-15 0.000000e+00    0    0    0    0
# [2,] 3.788921e-05 1.221245e-15 0.000000e+00    0    0    0    0
# [3,] 9.794982e-04 6.574630e-12 1.110223e-16    0    0    0    0
# [4,] 4.508909e-07 0.000000e+00 0.000000e+00    0    0    0    0
# [5,] 3.292044e-04 3.542722e-13 0.000000e+00    0    0    0    0
# [6,] 2.868102e-06 0.000000e+00 0.000000e+00    0    0    0    0
# [7,] 3.706208e-05 1.110223e-15 0.000000e+00    0    0    0    0
# [8,] 4.654677e-05 1.443290e-15 0.000000e+00    0    0    0    0
# [9,] 2.346233e-01 2.411004e-05 1.054719e-05    0    0    0    0
#[10,] 0.000000e+00 0.000000e+00 0.000000e+00    0    0    0    0
#[11,] 6.088156e-04 1.837752e-12 0.000000e+00    0    0    0    0
#[12,] 1.085873e-06 0.000000e+00 0.000000e+00    0    0    0    0

#$mu
#            [,1]       [,2]      [,3]      [,4]      [,5]      [,6]     [,7]
# [1,] 0.05910000 0.08300077 0.1187989 0.2759604 0.4866007 0.7781008 1.094870
# [2,] 0.05909999 0.08300004 0.1188045 0.2759851 0.4865945 0.7780999 1.095005
# [3,] 0.05909994 0.08300017 0.1188000 0.2759728 0.4865999 0.7781592 1.094998
# [4,] 0.05910067 0.08300007 0.1188000 0.2760228 0.4865710 0.7781000 1.094993
# [5,] 0.05909999 0.08299999 0.1188000 0.2759736 0.4866009 0.7781013 1.094982
# [6,] 0.05910001 0.08300016 0.1188001 0.2759752 0.4866002 0.7780987 1.095000
# [7,] 0.05909978 0.08300002 0.1187981 0.2759728 0.4866004 0.7781061 1.095002
# [8,] 0.05910002 0.08300015 0.1188000 0.2759727 0.4866002 0.7780941 1.094996
# [9,] 0.05910259 0.08299998 0.1187960 0.2759733 0.4866155 0.7781013 1.095002
#[10,] 0.05910001 0.08299656 0.1187994 0.2759727 0.4866007 0.7781002 1.095000
#[11,] 0.05910000 0.08299998 0.1187997 0.2759723 0.4866007 0.7780999 1.094993
#[12,] 0.05910005 0.08300001 0.1188000 0.2759759 0.4866009 0.7780957 1.094993

#$sigma2
#            [,1]       [,2]      [,3]      [,4]      [,5]      [,6]      [,7]
# [1,] 0.05940717 0.08361100 0.1230834 0.2661881 0.4582198 0.6338582 0.7408518
# [2,] 0.05940716 0.08361102 0.1230925 0.2662236 0.4582529 0.6338576 0.7405677
# [3,] 0.05940706 0.08361102 0.1230864 0.2661973 0.4582186 0.6338248 0.7405680
# [4,] 0.05940587 0.08361104 0.1230865 0.2662008 0.4582198 0.6338579 0.7405145
# [5,] 0.05940716 0.08361095 0.1230866 0.2662041 0.4582187 0.6338586 0.7405597
# [6,] 0.05940719 0.08361099 0.1230870 0.2661971 0.4582195 0.6338637 0.7405708
# [7,] 0.05940704 0.08361099 0.1230866 0.2661974 0.4582198 0.6337644 0.7405704
# [8,] 0.05942150 0.08361099 0.1230865 0.2661972 0.4582193 0.6338783 0.7405771
# [9,] 0.05941023 0.08360887 0.1230876 0.2661979 0.4581383 0.6338540 0.7405673
#[10,] 0.05940719 0.08360633 0.1230866 0.2661975 0.4582197 0.6338592 0.7405723
#[11,] 0.05940718 0.08361094 0.1230865 0.2661963 0.4582196 0.6338568 0.7405714
#[12,] 0.05940725 0.08361100 0.1230866 0.2661647 0.4582204 0.6338784 0.7405822

#$Dispersion
#          [,1]     [,2]     [,3]      [,4]      [,5]      [,6]      [,7]
# [1,] 1.011842 1.016729 1.049924 0.9951000 0.9954986 0.8917161 0.7703755
# [2,] 1.011842 1.016738 1.049953 0.9951466 0.9955824 0.8917162 0.7699980
# [3,] 1.011841 1.016736 1.049940 0.9950912 0.9954976 0.8916086 0.7700025
# [4,] 1.011808 1.016738 1.049941 0.9949299 0.9955557 0.8917165 0.7699501
# [5,] 1.011842 1.016738 1.049941 0.9951135 0.9954956 0.8917161 0.7700037
# [6,] 1.011842 1.016736 1.049944 0.9950819 0.9954988 0.8917260 0.7700042
# [7,] 1.011843 1.016738 1.049958 0.9950915 0.9954993 0.8915787 0.7700027
# [8,] 1.012086 1.016736 1.049941 0.9950913 0.9954983 0.8917513 0.7700133
# [9,] 1.011850 1.016712 1.049985 0.9950919 0.9952930 0.8917098 0.7699996
#[10,] 1.011842 1.016723 1.049947 0.9950919 0.9954983 0.8917181 0.7700058
#[11,] 1.011842 1.016738 1.049944 0.9950894 0.9954983 0.8917152 0.7700092
#[12,] 1.011842 1.016738 1.049941 0.9949585 0.9954994 0.8917498 0.7700206
#
#$Rd
#          [,1]     [,2]     [,3]      [,4]      [,5]      [,6]      [,7]
# [1,] 1.011893 1.016789 1.049994 0.9951926 0.9956008 0.8918092 0.7700856
# [2,] 1.011893 1.016789 1.049994 0.9951926 0.9956008 0.8918092 0.7700856
# [3,] 1.011893 1.016789 1.049994 0.9951926 0.9956008 0.8918092 0.7700856
# [4,] 1.011893 1.016789 1.049994 0.9951926 0.9956008 0.8918092 0.7700856
# [5,] 1.011893 1.016789 1.049994 0.9951926 0.9956008 0.8918092 0.7700856
# [6,] 1.011893 1.016789 1.049994 0.9951926 0.9956008 0.8918092 0.7700856
# [7,] 1.011893 1.016789 1.049994 0.9951926 0.9956008 0.8918092 0.7700856
# [8,] 1.011893 1.016789 1.049994 0.9951926 0.9956008 0.8918092 0.7700856
# [9,] 1.011893 1.016789 1.049994 0.9951926 0.9956008 0.8918092 0.7700856
#[10,] 1.011893 1.016789 1.049994 0.9951926 0.9956008 0.8918092 0.7700856
#[11,] 1.011893 1.016789 1.049994 0.9951926 0.9956008 0.8918092 0.7700856
#[12,] 1.011893 1.016789 1.049994 0.9951926 0.9956008 0.8918092 0.7700856


#$beta1
#             [,1]       [,2]       [,3]        [,4]       [,5]       [,6]
# [1,]  2.53821651  2.2048825  1.9731606  0.86534531  0.3090823 -0.3902616
# [2,]  2.50839839  2.1751558  1.9435128  0.83573194  0.2798412 -0.4195415
# [3,]  5.13230775  4.7953381  4.5502992  3.43798837  2.8516290  2.1275047
# [4,] -0.05806443 -0.3880153 -0.6089736 -1.71306703 -2.2405568 -2.9157729
# [5,]  3.37009772  3.0346240  2.7950083  1.68431665  1.1097155  0.3952386
# [6,]  1.70561808  1.3742787  1.1495214  0.04440328 -0.4942797 -1.1784701
# [7,]  2.43917878  2.1059393  1.8744627  0.76647795  0.2107148 -0.4886343
# [8,]  2.63943505  2.3039511  2.0720662  0.96406429  0.4074494 -0.2921581
# [9,]  5.46609254  5.0922000  4.7146909  3.66122770  2.9403282  2.1445169
#[10,] -0.47984456 -0.8026112 -0.9943757 -2.08167088 -2.5120298 -3.0870814
#[11,]  3.21382002  2.8775793  2.6353379  1.52365359  0.9435123  0.2245234
#[12,]  1.86136516  1.5306158  1.3078754  0.20326833 -0.3296989 -1.0091720
#            [,7]
# [1,] -1.0062080
# [2,] -1.0364110
# [3,]  1.4751078
# [4,] -3.4978017
# [5,] -0.2438379
# [6,] -1.7737535
# [7,] -1.1048544
# [8,] -0.9099126
# [9,]  1.3663096
#[10,] -3.5242940
#[11,] -0.4209478
#[12,] -1.5979150

#$beta2
#            [,1]       [,2]        [,3]        [,4]         [,5]        [,6]
# [1,]  0.2960267  0.2923943  0.19538648  0.39165723  0.368420193  0.45258880
# [2,]  0.2934255  0.2897223  0.19261241  0.38884331  0.365468458  0.44967080
# [3,] -0.1013538 -0.1018265 -0.18786563  0.01046144  0.005548245  0.10043183
# [4,]  0.6954394  0.6888792  0.58295984  0.77722470  0.736360050  0.80944283
# [5,]  0.1569877  0.1552212  0.06471628  0.26236352  0.250350711  0.34117458
# [6,]  0.4357402  0.4303822  0.32760867  0.52246488  0.488420893  0.56569553
# [7,]  0.2897147  0.2860063  0.18880365  0.38517791  0.361661967  0.44583508
# [8,]  0.3004079  0.2986898  0.20180778  0.39817472  0.375179891  0.45946812
# [9,] -0.4426777 -0.4123978 -0.39380437 -0.25568167 -0.187571556 -0.07326974
#[10,]  1.1182069  1.1054147  0.97522443  1.15951929  1.056974911  1.08021451
#[11,]  0.1465982  0.1454969  0.05715892  0.25528846  0.246627683  0.33939944
#[12,]  0.4466284  0.4407597  0.33631460  0.53094168  0.493324003  0.56841071
#            [,7]
# [1,] 0.51993800
# [2,] 0.51730995
# [3,] 0.18137566
# [4,] 0.86331308
# [5,] 0.41723745
# [6,] 0.62487892
# [7,] 0.51320287
# [8,] 0.52745385
# [9,] 0.04444384
#[10,] 1.07127302
#[11,] 0.41783101
#[12,] 0.62498377

#outcum.zelt8=cum.example.zelterman(tab.pala1)
#outcum.zelt8
#Error in optim(para0, nlk.semi.altham, gr = NULL, method = "BFGS", hessian = TRUE,  : 
#  non-finite finite-difference value [2]


out.zelt9=example.zelterman(tab.ec)
out.zelt9

#$AIC
#          [,1]     [,2]
# [1,] 511.6061 280.1801
# [2,] 510.5708 278.8052
# [3,] 511.7841 279.9387
# [4,] 511.6053 280.5052
# [5,] 512.7478 281.2816
# [6,] 510.7530 279.2575
# [7,] 510.7408 279.0586
# [8,] 512.6132 281.4709
# [9,] 513.5595 279.7287
#[10,] 512.8024 282.1181
#[11,] 513.0742 281.5266
#[12,] 510.6146 279.1112
#
#$Test.Statistics
#          [,1]     [,2]
# [1,] 3.221594 4.402125
# [2,] 2.498911 3.016040
# [3,] 3.405872 4.138912
# [4,] 3.220487 4.752348
# [5,] 4.081931 5.509396
# [6,] 2.627023 3.479669
# [7,] 2.618449 3.271690
# [8,] 3.946790 5.722354
# [9,] 5.248221 3.878113
#[10,] 4.465176 6.480474
#[11,] 4.343365 5.752939
#[12,] 2.542629 3.335187

#$P.value
#           [,1]       [,2]
# [1,] 0.3587052 0.22118846
# [2,] 0.4754880 0.38915867
# [3,] 0.3331770 0.24684957
# [4,] 0.3588636 0.19085599
# [5,] 0.2527519 0.13807773
# [6,] 0.4527717 0.32340863
# [7,] 0.4542644 0.35160234
# [8,] 0.2672675 0.12592797
# [9,] 0.1544978 0.27493012
#[10,] 0.2154176 0.09043576
#[11,] 0.2266945 0.12426887
#[12,] 0.4676363 0.34277526

#$mu
#          [,1]     [,2]
# [1,] 1.510101 1.925532
# [2,] 1.509939 1.925655
# [3,] 1.510102 1.925536
# [4,] 1.510103 1.925538
# [5,] 1.510102 1.925534
# [6,] 1.510103 1.925534
# [7,] 1.510101 1.925607
# [8,] 1.510103 1.925535
# [9,] 1.510105 1.925535
#[10,] 1.510103 1.925538
#[11,] 1.510100 1.925540
#[12,] 1.510102 1.925534

#$sigma2
#           [,1]     [,2]
# [1,] 0.7751435 1.154014
# [2,] 0.7758233 1.153986
# [3,] 0.7751474 1.154030
# [4,] 0.7751424 1.154016
# [5,] 0.7751443 1.154014
# [6,] 0.7751456 1.154015
# [7,] 0.7751444 1.153569
# [8,] 0.7751447 1.154015
# [9,] 0.7751479 1.154017
#[10,] 0.7751432 1.154016
#[11,] 0.7751429 1.154131
#[12,] 0.7751449 1.154015

#$Dispersion
#           [,1]     [,2]
# [1,] 0.8246210 1.155616
# [2,] 0.8253790 1.155583
# [3,] 0.8246250 1.155632
# [4,] 0.8246194 1.155618
# [5,] 0.8246215 1.155616
# [6,] 0.8246228 1.155617
# [7,] 0.8246220 1.155168
# [8,] 0.8246218 1.155617
# [9,] 0.8246247 1.155619
#[10,] 0.8246202 1.155617
#[11,] 0.8246206 1.155733
#[12,] 0.8246222 1.155617

#$Rd
#           [,1]     [,2]
# [1,] 0.8288143 1.168057
# [2,] 0.8288143 1.168057
# [3,] 0.8288143 1.168057
# [4,] 0.8288143 1.168057
# [5,] 0.8288143 1.168057
# [6,] 0.8288143 1.168057
# [7,] 0.8288143 1.168057
# [8,] 0.8288143 1.168057
# [9,] 0.8288143 1.168057
#[10,] 0.8288143 1.168057
#[11,] 0.8288143 1.168057
#[12,] 0.8288143 1.168057


#$beta1
#             [,1]       [,2]
# [1,] -1.75551383 -1.3561340
# [2,] -1.76635880 -1.3112040
# [3,]  0.03889354  0.4429674
# [4,] -3.54993006 -3.1549273
# [5,] -1.03152375 -0.6754437
# [6,] -2.48217398 -2.0390309
# [7,] -1.89269023 -1.4471304
# [8,] -1.62121278 -1.2678781
# [9,] -0.13501154  0.2861524
#[10,] -3.37624484 -2.9934768
#[11,] -1.22018876 -0.8730137
#[12,] -2.29479400 -1.8423667

#$beta2
#           [,1]        [,2]
# [1,] 0.5924578  0.35514968
# [2,] 0.5474392  0.29125099
# [3,] 0.1431957 -0.09462739
# [4,] 1.0417349  0.80484723
# [5,] 0.5018334  0.27965269
# [6,] 0.6840067  0.43114455
# [7,] 0.5362429  0.28324883
# [8,] 0.6497782  0.42784777
# [9,] 0.1846456 -0.05542523
#[10,] 1.0005784  0.76448661
#[11,] 0.5060349  0.28705773
#[12,] 0.6802329  0.42390500


# example.zelterman(tab4.6rao)
# Error in matrix(NA, byrow = T, nrow = 12, ncol = dim.freq[1]) : 
#  non-numeric matrix extent


out.zelt10=example.zelterman(tab5.7rao)
out.zelt10

#$AIC
#          [,1]     [,2]
# [1,] 1211.531 1613.159
# [2,] 1202.095 1613.979
# [3,] 1210.448 1613.322
# [4,] 1212.945 1613.325
# [5,] 1219.073 1613.654
# [6,] 1204.818 1613.702
# [7,] 1203.856 1613.701
# [8,] 1219.822 1613.644
# [9,] 1209.199 1616.047
#[10,] 1219.753 1616.109
#[11,] 1220.772 1613.926
#[12,] 1203.602 1613.980

#$Test.Statistics
#          [,1]         [,2]
# [1,] 56.55053 0.0003487227
# [2,] 47.45492 0.8154185434
# [3,] 54.95794 0.1633096886
# [4,] 58.51520 0.1658775354
# [5,] 63.58440 0.4933060879
# [6,] 50.41758 0.5467733787
# [7,] 49.16507 0.5400046194
# [8,] 64.79537 0.4883092913
# [9,] 52.27413 2.9286868707
#[10,] 67.39822 2.9391536782
#[11,] 65.10824 0.7621091769
#[12,] 49.36111 0.8277867752

#$P.value
#              [,1]      [,2]
# [1,] 3.205214e-12 0.9999983
# [2,] 2.781344e-10 0.8457757
# [3,] 7.009504e-12 0.9832830
# [4,] 1.220135e-12 0.9829002
# [5,] 1.006972e-13 0.9203596
# [6,] 6.509726e-11 0.9085017
# [7,] 1.203061e-10 0.9100187
# [8,] 5.551115e-14 0.9214522
# [9,] 2.618084e-11 0.4027514
#[10,] 1.543210e-14 0.4011019
#[11,] 4.751755e-14 0.8585069
#[12,] 1.092818e-10 0.8428099

#$mu
#          [,1]     [,2]
# [1,] 1.478087 1.961881
# [2,] 1.478225 1.961882
# [3,] 1.478059 1.961884
# [4,] 1.478087 1.961886
# [5,] 1.478083 1.961889
# [6,] 1.477990 1.961881
# [7,] 1.478123 1.961884
# [8,] 1.478088 1.961866
# [9,] 1.478082 1.961888
#[10,] 1.478085 1.961885
#[11,] 1.478082 1.961885
#[12,] 1.478087 1.961886

#$sigma2
#          [,1]     [,2]
# [1,] 1.252972 2.007292
# [2,] 1.252908 2.007283
# [3,] 1.252911 2.007284
# [4,] 1.252971 2.007284
# [5,] 1.252968 2.007260
# [6,] 1.252760 2.007293
# [7,] 1.252931 2.007281
# [8,] 1.252972 2.007302
# [9,] 1.252967 2.007286
#[10,] 1.252970 2.007283
#[11,] 1.252966 2.007276
#[12,] 1.252972 2.007291

#$Dispersion
#          [,1]     [,2]
# [1,] 1.344532 2.008021
# [2,] 1.344411 2.008013
# [3,] 1.344477 2.008013
# [4,] 1.344532 2.008013
# [5,] 1.344529 2.007989
# [6,] 1.344341 2.008022
# [7,] 1.344475 2.008010
# [8,] 1.344532 2.008032
# [9,] 1.344529 2.008015
#[10,] 1.344530 2.008013
#[11,] 1.344528 2.008006
#[12,] 1.344533 2.008020

#$Rd
#          [,1]     [,2]
# [1,] 1.347832 2.012060
# [2,] 1.347832 2.012060
# [3,] 1.347832 2.012060
# [4,] 1.347832 2.012060
# [5,] 1.347832 2.012060
# [6,] 1.347832 2.012060
# [7,] 1.347832 2.012060
# [8,] 1.347832 2.012060
# [9,] 1.347832 2.012060
#[10,] 1.347832 2.012060
#[11,] 1.347832 2.012060
#[12,] 1.347832 2.012060


#$beta1
#            [,1]        [,2]
# [1,] -0.6094146  0.03074770
# [2,] -0.5979418  0.07790620
# [3,]  1.1998173  1.84720752
# [4,] -2.4181556 -1.78506379
# [5,]  0.1074174  0.72596332
# [6,] -1.3271225 -0.66395793
# [7,] -0.7276318 -0.05932937
# [8,] -0.4924395  0.12164025
# [9,]  1.0779384  1.74909877
#[10,] -2.2868637 -1.67564234
#[11,] -0.0809504  0.53367026
#[12,] -1.1386924 -0.47134790

#$beta2
#            [,1]         [,2]
# [1,]  0.2507033 -0.002940878
# [2,]  0.1948130 -0.071250442
# [3,] -0.2017038 -0.457056215
# [4,]  0.7029983  0.451011447
# [5,]  0.1665520 -0.078942814
# [6,]  0.3352293  0.072900586
# [7,]  0.1853120 -0.078220017
# [8,]  0.3166539  0.072177652
# [9,] -0.1715846 -0.432529467
#[10,]  0.6705340  0.423656166
#[11,]  0.1717604 -0.072144770
#[12,]  0.3300060  0.065993469

out.zelt11=example.zelterman(tab.hgb)
out.zelt11

#$AIC
#          [,1]     [,2]
# [1,] 128.3794 67.46623
# [2,] 127.7237 67.12026
# [3,] 128.4769 67.32132
# [4,] 128.3193 67.61691
# [5,] 129.0501 67.66373
# [6,] 127.7907 67.29465
# [7,] 127.8310 67.19294
# [8,] 128.9820 67.76507
# [9,] 129.0433 66.88866
#[10,] 128.3793 68.13973
#[11,] 129.2237 67.69593
#[12,] 127.6667 67.27569

#$Test.Statistics
#          [,1]     [,2]
# [1,] 3.048176 3.182321
# [2,] 2.507943 2.744363
# [3,] 3.127497 3.081722
# [4,] 3.008652 3.312830
# [5,] 3.610522 3.489423
# [6,] 2.564069 2.904311
# [7,] 2.594965 2.831926
# [8,] 3.560063 3.583387
# [9,] 3.655480 2.700743
#[10,] 3.152513 3.781596
#[11,] 3.756939 3.554386
#[12,] 2.464108 2.875876

#$P.value
#           [,1]      [,2]
# [1,] 0.3842569 0.3643600
# [2,] 0.4738573 0.4327408
# [3,] 0.3723823 0.3791960
# [4,] 0.3902931 0.3458609
# [5,] 0.3067082 0.3221366
# [6,] 0.4638236 0.4066151
# [7,] 0.4583730 0.4182712
# [8,] 0.3130553 0.3101071
# [9,] 0.3011497 0.4401010
#[10,] 0.3687031 0.2860341
#[11,] 0.2889347 0.3137766
#[12,] 0.4818125 0.4111613

#$mu
#          [,1]     [,2]
# [1,] 1.454545 1.916665
# [2,] 1.454548 1.916922
# [3,] 1.454545 1.914443
# [4,] 1.454552 1.916667
# [5,] 1.454548 1.916945
# [6,] 1.454547 1.917201
# [7,] 1.454545 1.916669
# [8,] 1.454545 1.916666
# [9,] 1.454548 1.916669
#[10,] 1.454547 1.916675
#[11,] 1.454544 1.916666
#[12,] 1.454548 1.916671

#$sigma2
#          [,1]     [,2]
# [1,] 1.066105 3.243034
# [2,] 1.066106 3.242165
# [3,] 1.066112 3.247832
# [4,] 1.066109 3.243036
# [5,] 1.066107 3.242414
# [6,] 1.066105 3.240201
# [7,] 1.066104 3.243037
# [8,] 1.066117 3.243036
# [9,] 1.066093 3.243037
#[10,] 1.066105 3.243038
#[11,] 1.066097 3.243036
#[12,] 1.066106 3.243038

#$Dispersion
#          [,1]     [,2]
# [1,] 1.151774 3.248674
# [2,] 1.151774 3.247769
# [3,] 1.151782 3.253787
# [4,] 1.151777 3.248676
# [5,] 1.151776 3.248015
# [6,] 1.151774 3.245764
# [7,] 1.151773 3.248677
# [8,] 1.151788 3.248676
# [9,] 1.151760 3.248676
#[10,] 1.151774 3.248677
#[11,] 1.151766 3.248676
#[12,] 1.151774 3.248677

#$Rd
#          [,1]     [,2]
# [1,] 1.178571 3.389943
# [2,] 1.178571 3.389943
# [3,] 1.178571 3.389943
# [4,] 1.178571 3.389943
# [5,] 1.178571 3.389943
# [6,] 1.178571 3.389943
# [7,] 1.178571 3.389943
# [8,] 1.178571 3.389943
# [9,] 1.178571 3.389943
#[10,] 1.178571 3.389943
#[11,] 1.178571 3.389943
#[12,] 1.178571 3.389943

#$beta1
#            [,1]      [,2]
# [1,] -0.9275323 1.9782011
# [2,] -0.9255939 2.0205709
# [3,]  0.8784891 3.8174106
# [4,] -2.7330762 0.1511868
# [5,] -0.2059580 2.6848421
# [6,] -1.6504329 1.2662127
# [7,] -1.0533117 1.8846969
# [8,] -0.8035899 2.0740966
# [9,]  0.7452804 3.7438010
#[10,] -2.5921169 0.2246802
#[11,] -0.3937200 2.4983311
#[12,] -1.4631921 1.4616679

#$beta2
#             [,1]        [,2]
# [1,]  0.35277109 -0.48813444
# [2,]  0.30103075 -0.55774892
# [3,] -0.09897757 -0.94777597
# [4,]  0.80440064 -0.03138099
# [5,]  0.26598093 -0.56499146
# [6,]  0.44008060 -0.41000013
# [7,]  0.29070577 -0.56458136
# [8,]  0.41562523 -0.41227636
# [9,] -0.06643185 -0.92953473
#[10,]  0.76997075 -0.04975495
#[11,]  0.27071021 -0.55915930
#[12,]  0.43555038 -0.41802155


out.zelt12=example.zelterman(tab.hauk1)
out.zelt12

#$AIC
#          [,1]     [,2]     [,3]
# [1,] 35.49075 9.675748 11.29724
# [2,] 35.50348 9.675143 11.29210
# [3,] 35.02428 9.665299 11.28825
# [4,] 36.04206 9.687189 11.30833
# [5,] 35.22753 9.671938 11.29699
# [6,] 35.78391 9.679922 11.29751
# [7,] 35.51544 9.674840 11.29239
# [8,] 35.46673 9.676688 11.30272
# [9,] 34.72552 9.664510 11.28058
#[10,] 37.06305 9.690748 11.32857
#[11,] 35.24183 9.672427 11.29708
#[12,] 35.77234 9.679319 11.29740

#$Test.Statistics
#          [,1]     [,2]      [,3]
# [1,] 5.700217 2.132725 0.5705324
# [2,] 5.716917 2.132102 0.5654160
# [3,] 5.055952 2.122140 0.5615852
# [4,] 6.510816 2.144780 0.5819185
# [5,] 5.331603 2.128877 0.5699221
# [6,] 6.126318 2.137129 0.5705153
# [7,] 5.733066 2.131797 0.5655895
# [8,] 5.668676 2.133706 0.5756774
# [9,] 4.608388 2.121328 0.5542165
#[10,] 8.262404 2.148588 0.6067979
#[11,] 5.344393 2.129341 0.5701660
#[12,] 6.116986 2.136463 0.5707556

#$P.value
#           [,1]      [,2]      [,3]
# [1,] 0.7695063 0.9891905 0.9999464
# [2,] 0.7678820 0.9892020 0.9999485
# [3,] 0.8294038 0.9893854 0.9999499
# [4,] 0.6879032 0.9889656 0.9999417
# [5,] 0.8044968 0.9892617 0.9999467
# [6,] 0.7272123 0.9891087 0.9999464
# [7,] 0.7663083 0.9892077 0.9999484
# [8,] 0.7725657 0.9891723 0.9999444
# [9,] 0.8670239 0.9894003 0.9999527
#[10,] 0.5079340 0.9888939 0.9999304
#[11,] 0.8033130 0.9892531 0.9999466
#[12,] 0.7281552 0.9891211 0.9999464

#$mu
#          [,1]     [,2]     [,3]
# [1,] 3.222207 3.999985 6.999976
# [2,] 3.222246 3.999919 6.999958
# [3,] 3.222278 4.000003 6.999942
# [4,] 3.222216 3.999981 7.000005
# [5,] 3.222222 3.999986 6.999918
# [6,] 3.222250 4.000006 6.999895
# [7,] 3.222239 4.000002 6.999912
# [8,] 3.222188 3.999977 6.999874
# [9,] 3.222216 3.999988 6.999956
#[10,] 3.222216 4.000050 6.999999
#[11,] 3.222222 3.999987 6.999919
#[12,] 3.222222 3.999983 6.999979

#$sigma2
#          [,1]      [,2]      [,3]
# [1,] 1.950430 1.0000894 0.6669011
# [2,] 1.950459 1.0006285 0.6669752
# [3,] 1.950444 1.0000249 0.6675833
# [4,] 1.950464 1.0001344 0.6666772
# [5,] 1.950420 1.0000862 0.6675333
# [6,] 1.950471 0.9999807 0.6673857
# [7,] 1.950455 1.0000172 0.6673082
# [8,] 1.950441 1.0001782 0.6675786
# [9,] 1.950500 1.0000579 0.6670277
#[10,] 1.950505 0.9996609 0.6667713
#[11,] 1.950422 1.0000755 0.6672750
#[12,] 1.950457 0.9999996 0.6667829
#
#$Dispersion
#           [,1]      [,2]      [,3]
# [1,] 0.8930764 0.4167044 0.3175705
# [2,] 0.8930840 0.4169314 0.3176047
# [3,] 0.8930725 0.4166770 0.3178933
# [4,] 0.8930907 0.4167233 0.3174656
# [5,] 0.8930696 0.4167030 0.3178680
# [6,] 0.8930889 0.4166584 0.3177964
# [7,] 0.8930834 0.4166738 0.3177605
# [8,] 0.8930842 0.4167417 0.3178870
# [9,] 0.8931071 0.4166912 0.3176296
#[10,] 0.8931093 0.4165236 0.3175101
#[11,] 0.8930706 0.4166986 0.3177451
#[12,] 0.8930866 0.4166671 0.3175144
#
#$Rd
#          [,1]      [,2]      [,3]
# [1,] 1.004805 0.8333333 0.4761905
# [2,] 1.004805 0.8333333 0.4761905
# [3,] 1.004805 0.8333333 0.4761905
# [4,] 1.004805 0.8333333 0.4761905
# [5,] 1.004805 0.8333333 0.4761905
# [6,] 1.004805 0.8333333 0.4761905
# [7,] 1.004805 0.8333333 0.4761905
# [8,] 1.004805 0.8333333 0.4761905
# [9,] 1.004805 0.8333333 0.4761905
#[10,] 1.004805 0.8333333 0.4761905
#[11,] 1.004805 0.8333333 0.4761905
#[12,] 1.004805 0.8333333 0.4761905


#$beta1
#            [,1]      [,2]       [,3]
# [1,] -1.5909768 -3.999476 -10.494504
# [2,] -1.6070716 -4.006694 -10.277360
# [3,]  0.5394045 -2.056464  -8.222680
# [4,] -3.7520673 -5.950644 -12.784289
# [5,] -1.0901890 -3.607046 -10.245164
# [6,] -2.1035928 -4.395302 -10.726144
# [7,] -1.6622756 -4.058612 -10.246052
# [8,] -1.5200992 -3.940486 -10.733690
# [9,] -0.7387136 -3.490611  -9.604236
#[10,] -2.7576377 -4.558216 -11.636399
#[11,] -1.2979985 -3.810873 -10.398005
#[12,] -1.8970811 -4.190859 -10.586916
#
#$beta2
#           [,1]      [,2]      [,3]
# [1,] 0.2481568 0.4999372 0.7496087
# [2,] 0.2437376 0.4945664 0.7259132
# [3,] 0.0301915 0.3046520 0.5314461
# [4,] 0.4703431 0.6961783 0.9690962
# [5,] 0.2122636 0.4770290 0.7408297
# [6,] 0.2857111 0.5232503 0.7571382
# [7,] 0.2385817 0.4890811 0.7131531
# [8,] 0.2578023 0.5108117 0.7854345
# [9,] 0.1540534 0.4479901 0.6710803
#[10,] 0.3836495 0.5577448 0.8466896
#[11,] 0.2222810 0.4867515 0.7456717
#[12,] 0.2758241 0.5134548 0.7532596

outcum.hauk=cum.example.zelterman(tab.hauk1)
outcum.hauk

#$AIC
#          [,1]     [,2]     [,3]     [,4]     [,5]     [,6]
# [1,] 18217.01 21148.73 21590.71 21643.81 21658.06 21698.22
# [2,] 18218.70 21150.83 21594.11 21647.18 21661.41 21701.73
# [3,] 18205.18 21129.27 21567.12 21621.00 21635.69 21672.28
# [4,] 18324.48 21284.86 21738.22 21791.37 21805.47 21851.11
# [5,] 18197.07 21122.02 21559.71 21613.25 21627.74 21665.63
# [6,] 18271.57 21217.89 21667.12 21720.12 21734.25 21777.53
# [7,] 18220.05 21152.57 21596.17 21649.23 21663.45 21703.91
# [8,] 18214.28 21145.34 21585.81 21638.96 21653.24 21693.14
# [9,] 18926.92 21969.65 22426.77 22486.49 22503.16 22537.17
#[10,] 19044.71 22161.13 22666.61 22722.05 22736.83 22800.79
#[11,] 18197.97 21122.93 21560.48 21614.08 21628.59 21666.31
#[12,] 18282.27 21230.65 21680.48 21733.46 21747.59 21791.12

#$Test.Statistics
#              [,1]         [,2]         [,3]         [,4]         [,5]
# [1,]    29.726184    37.681024 4.398394e+01 4.345069e+01 4.289297e+01
# [2,]    32.465515    40.561330 5.258739e+01 5.178423e+01 5.111504e+01
# [3,]    13.903775    14.626691 1.525800e+01 1.580159e+01 1.578577e+01
# [4,]   234.866857   254.381405 3.262772e+02 3.189922e+02 3.156097e+02
# [5,]     5.996180     7.705323 8.252748e+00 8.476908e+00 8.286310e+00
# [6,]   116.111044   133.299388 1.830416e+02 1.793253e+02 1.774297e+02
# [7,]    34.642294    42.936762 5.659517e+01 5.557212e+01 5.492093e+01
# [8,]    25.676649    33.285918 3.622195e+01 3.590457e+01 3.544147e+01
# [9,]   729.860193   836.075170 8.562893e+02 8.619321e+02 8.638331e+02
#[10,] 97794.826284 39557.893441 2.765542e+05 2.506096e+05 2.413336e+05
#[11,]     6.876993     8.615210 9.042005e+00 9.340259e+00 9.159232e+00
#[12,]   131.509927   149.817824 2.020928e+02 1.980347e+02 1.960238e+02
#              [,6]
# [1,] 4.621369e+01
# [2,] 5.405807e+01
# [3,] 1.517415e+01
# [4,] 3.205104e+02
# [5,] 9.039678e+00
# [6,] 1.806783e+02
# [7,] 5.810048e+01
# [8,] 3.869110e+01
# [9,] 8.603425e+02
#[10,] 1.816360e+05
#[11,] 9.735530e+00
#[12,] 1.993047e+02

#$P.value
#              [,1]         [,2]         [,3]         [,4]         [,5]
# [1,] 4.883417e-04 1.988418e-05 1.421085e-06 1.781446e-06 2.255327e-06
# [2,] 1.653976e-04 6.010336e-06 3.499492e-08 4.965449e-08 6.642181e-08
# [3,] 1.257895e-01 1.017139e-01 8.408860e-02 7.114227e-02 7.149235e-02
# [4,] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
# [5,] 7.403005e-01 5.640936e-01 5.088898e-01 4.868889e-01 5.055705e-01
# [6,] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
# [7,] 6.888892e-05 2.213987e-06 6.040427e-09 9.473879e-09 1.260949e-08
# [8,] 2.307133e-03 1.190759e-04 3.620735e-05 4.122244e-05 4.979098e-05
# [9,] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#[10,] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#[11,] 6.499245e-01 4.735241e-01 4.334058e-01 4.064765e-01 4.227075e-01
#[12,] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#              [,6]
# [1,] 5.496259e-07
# [2,] 1.840441e-08
# [3,] 8.626204e-02
# [4,] 0.000000e+00
# [5,] 4.336197e-01
# [6,] 0.000000e+00
# [7,] 3.109102e-09
# [8,] 1.309570e-05
# [9,] 0.000000e+00
#[10,] 0.000000e+00
#[11,] 3.723168e-01
#[12,] 0.000000e+00

#$mu
#           [,1]      [,2]      [,3]      [,4]      [,5]      [,6]
# [1,] 0.9135448 0.9666487 0.9855053 0.9880568 0.9888198 0.9910967
# [2,] 0.9135076 0.9666392 0.9854590 0.9880075 0.9887704 0.9910641
# [3,] 0.9135124 0.9666284 0.9854299 0.9879708 0.9887340 0.9910223
# [4,] 0.9134907 0.9666276 0.9854694 0.9880193 0.9887825 0.9910759
# [5,] 0.9135100 0.9666285 0.9854460 0.9879885 0.9887520 0.9910849
# [6,] 0.9130129 0.9666229 0.9854507 0.9879954 0.9887566 0.9910344
# [7,] 0.9135259 0.9666153 0.9854521 0.9879734 0.9887491 0.9911737
# [8,] 0.9134966 0.9665224 0.9854600 0.9879877 0.9887504 0.9910311
# [9,] 0.9134896 0.9666570 0.9855971 0.9881346 0.9889060 0.9910993
#[10,] 0.9134499 0.9666243 0.9854450 0.9879899 0.9887512 0.9910295
#[11,] 0.9135793 0.9666771 0.9854675 0.9879910 0.9887512 0.9910489
#[12,] 0.9135136 0.9666111 0.9854484 0.9880006 0.9887615 0.9910391

#$sigma2
#          [,1]     [,2]     [,3]     [,4]     [,5]     [,6]
# [1,] 1.359084 1.452353 1.503450 1.509635 1.511796 1.525151
# [2,] 1.359187 1.452302 1.503461 1.509663 1.511829 1.525233
# [3,] 1.359144 1.452319 1.503240 1.509274 1.511467 1.524976
# [4,] 1.359152 1.452308 1.503605 1.509815 1.511984 1.525392
# [5,] 1.359181 1.452266 1.503371 1.509555 1.511732 1.525263
# [6,] 1.358999 1.452318 1.503471 1.509663 1.511825 1.525190
# [7,] 1.359106 1.452240 1.503454 1.510154 1.512040 1.524628
# [8,] 1.359167 1.452091 1.503455 1.509598 1.511766 1.525163
# [9,] 1.358691 1.451986 1.503760 1.509912 1.512129 1.525054
#[10,] 1.359238 1.452311 1.503491 1.509685 1.511849 1.525218
#[11,] 1.359487 1.452336 1.503609 1.509564 1.511723 1.525146
#[12,] 1.359213 1.452247 1.503489 1.509696 1.511857 1.525217

#$Dispersion
#          [,1]     [,2]     [,3]     [,4]     [,5]     [,6]
# [1,] 1.637276 1.663238 1.692343 1.695398 1.696658 1.708146
# [2,] 1.637460 1.663195 1.692427 1.695504 1.696771 1.708287
# [3,] 1.637400 1.663231 1.692223 1.695123 1.696420 1.708064
# [4,] 1.637446 1.663220 1.692573 1.695657 1.696926 1.708447
# [5,] 1.637450 1.663170 1.692346 1.695413 1.696689 1.708289
# [6,] 1.638032 1.663238 1.692451 1.695523 1.696788 1.708285
# [7,] 1.637333 1.663160 1.692430 1.696108 1.697040 1.707442
# [8,] 1.637453 1.663133 1.692419 1.695462 1.696730 1.708260
# [9,] 1.636892 1.662806 1.692553 1.695590 1.696901 1.708033
#[10,] 1.637615 1.663229 1.692482 1.695556 1.696822 1.708324
#[11,] 1.637706 1.663175 1.692580 1.695419 1.696681 1.708214
#[12,] 1.637481 1.663175 1.692475 1.695552 1.696816 1.708308

#$Rd
#          [,1]     [,2]     [,3]     [,4]     [,5]     [,6]
# [1,] 1.637756 1.663472 1.692732 1.695806 1.697072 1.708576
# [2,] 1.637756 1.663472 1.692732 1.695806 1.697072 1.708576
# [3,] 1.637756 1.663472 1.692732 1.695806 1.697072 1.708576
# [4,] 1.637756 1.663472 1.692732 1.695806 1.697072 1.708576
# [5,] 1.637756 1.663472 1.692732 1.695806 1.697072 1.708576
# [6,] 1.637756 1.663472 1.692732 1.695806 1.697072 1.708576
# [7,] 1.637756 1.663472 1.692732 1.695806 1.697072 1.708576
# [8,] 1.637756 1.663472 1.692732 1.695806 1.697072 1.708576
# [9,] 1.637756 1.663472 1.692732 1.695806 1.697072 1.708576
#[10,] 1.637756 1.663472 1.692732 1.695806 1.697072 1.708576
#[11,] 1.637756 1.663472 1.692732 1.695806 1.697072 1.708576
#[12,] 1.637756 1.663472 1.692732 1.695806 1.697072 1.708576

#$beta1
#             [,1]        [,2]        [,3]        [,4]         [,5]         [,6]
# [1,]  0.51152953  0.47213563  0.47202839  0.47142013  0.471463059  0.476539878
# [2,]  0.48898384  0.44975770  0.45013158  0.44958081  0.449644324  0.454860952
# [3,]  2.93185770  2.88226108  2.87533755  2.87390841  2.873687312  2.876760294
# [4,] -1.93846743 -1.96707833 -1.95962845 -1.95938484 -1.959028044 -1.951690359
# [5,]  1.18129162  1.13522846  1.13053532  1.12942707  1.129283174  1.132872828
# [6,] -0.17028936 -0.20427530 -0.19980509 -0.19989198 -0.199656332 -0.193182306
# [7,]  0.42803747  0.38914952  0.38972872  0.38962713  0.389462646  0.393885864
# [8,]  0.59467254  0.55449500  0.55374910  0.55308396  0.553102353  0.557940811
# [9,]  2.40227263  2.32877814  2.31021585  2.30746054  2.306750254  2.307130995
#[10,] -1.61826483 -1.61628726 -1.59058830 -1.58824675 -1.587132623 -1.574636830
#[11,]  0.99191204  0.94503323  0.94015810  0.93891143  0.938750516  0.942270625
#[12,]  0.01704650 -0.01513591 -0.01015207 -0.01019000 -0.009937005 -0.003354151

#$beta2
#               [,1]         [,2]         [,3]         [,4]         [,5]
# [1,]  0.0566846485  0.057120710  0.053765694  0.053481266  0.053337885
# [2,]  0.0531479955  0.053534195  0.050048957  0.049749143  0.049599986
# [3,] -0.2060723505 -0.203127104 -0.204537579 -0.204599498 -0.204667651
# [4,]  0.3307124990  0.328054597  0.322365943  0.321826608  0.321586835
# [5,] -0.0057571395 -0.003621438 -0.005644835 -0.005788831 -0.005878517
# [6,]  0.1237445286  0.122566201  0.117822246  0.117385499  0.117184209
# [7,]  0.0488714680  0.049173490  0.045624882  0.045211222  0.045119566
# [8,]  0.0646367100  0.065297486  0.062113304  0.061849363  0.061713210
# [9,] -0.1773743063 -0.170790115 -0.169747299 -0.169554486 -0.169518573
#[10,]  0.3913072168  0.378580176  0.366585670  0.365340317  0.364841091
#[11,]  0.0009778857  0.003315536  0.001354147  0.001248311  0.001163858
#[12,]  0.1176769046  0.116084907  0.111175521  0.110722048  0.110515207
#               [,6]
# [1,]  5.162342e-02
# [2,]  4.784117e-02
# [3,] -2.057617e-01
# [4,]  3.191336e-01
# [5,] -7.145040e-03
# [6,]  1.150210e-01
# [7,]  4.354279e-02
# [8,]  6.006956e-02
# [9,] -1.698457e-01
#[10,]  3.605353e-01
#[11,] -7.297024e-05
#[12,]  1.083128e-01




# example.zelterman(tab.hauk2)
# Error in optim(para0, nlk.semi.altham, gr = NULL, method = "BFGS", hessian = TRUE,  : 
#  non-finite finite-difference value [2]


out.zelt13=example.zelterman(tab.omariba)
out.zelt13

#$AIC
#          [,1]     [,2]     [,3]     [,4]     [,5]     [,6]     [,7]
# [1,] 779.8416 756.0680 835.4854 677.7268 626.1789 499.8964 759.3826
# [2,] 779.8374 756.0510 835.4811 677.7100 626.1447 499.8693 759.2676
# [3,] 780.7591 758.2953 836.1792 679.1400 628.2486 501.1356 764.8397
# [4,] 779.9168 755.4518 837.0921 677.9685 626.9825 501.1537 758.3656
# [5,] 780.2190 757.1292 835.5594 678.2622 627.0145 500.3334 762.0501
# [6,] 779.8080 755.5781 836.2366 677.7713 626.3878 500.3878 758.3184
# [7,] 779.8339 756.0362 835.4808 677.6952 626.1215 499.8554 759.1824
# [8,] 779.8500 756.1019 835.4940 677.7597 626.2452 499.9488 759.6003
# [9,] 798.1536 783.8527 861.9294 700.9980 656.9825 523.5586 808.0453
#[10,] 784.6991 761.1884 852.9912 689.9701 643.5860 516.2720 777.0507
#[11,] 780.3569 757.3955 835.6316 678.2971 627.1358 500.4038 762.3520
#[12,] 779.8463 755.5732 836.5164 678.0123 626.6397 500.6262 758.4970
#
#$Test.Statistics
#            [,1]      [,2]      [,3]     [,4]      [,5]      [,6]      [,7]
# [1,]  0.4176408  1.958125  1.400489 13.46919  1.780718  4.320277  5.689091
# [2,]  0.4148792  1.941326  1.398297 13.44186  1.753304  4.274223  5.581012
# [3,]  1.1886459  4.306668  1.874741 16.23680  3.545792  5.568812 10.800509
# [4,]  0.6655332  1.445451  3.380358 13.09458  2.974753  6.205553  5.206793
# [5,]  0.6987443  3.029710  1.337173 14.65070  2.417139  4.748301  8.124005
# [6,]  0.4832649  1.532240  2.327315 13.10904  2.214661  5.092843  4.932394
# [7,]  0.4128168  1.926589  1.400426 13.41863  1.736574  4.255723  5.503667
# [8,]  0.4230524  1.992750  1.405197 13.51997  1.835032  4.411992  5.895135
# [9,] 24.3574592 41.644910 34.030559 58.94825 36.726736 34.376496 56.430684
#[10,] 10.4864541  7.316261 42.169272 23.72039 24.246154 32.534462 27.242773
#[11,]  0.8178798  3.299969  1.379583 14.82070  2.512510  4.767851  8.380559
#[12,]  0.5471841  1.541040  2.656816 13.28417  2.496669  5.420991  5.146760

#$P.value
#             [,1]         [,2]         [,3]         [,4]         [,5]
# [1,] 0.999986000 9.921126e-01 9.978199e-01 1.424993e-01 9.944798e-01
# [2,] 0.999986396 9.923618e-01 9.978333e-01 1.436104e-01 9.947955e-01
# [3,] 0.998865308 8.900945e-01 9.932970e-01 6.209820e-02 9.386916e-01
# [4,] 0.999896940 9.975316e-01 9.472906e-01 1.583730e-01 9.652875e-01
# [5,] 0.999873405 9.631054e-01 9.981841e-01 1.009885e-01 9.830295e-01
# [6,] 0.999973711 9.969005e-01 9.851728e-01 1.577340e-01 9.876003e-01
# [7,] 0.999986687 9.925760e-01 9.978202e-01 1.445603e-01 9.949820e-01
# [8,] 0.999985197 9.915817e-01 9.977909e-01 1.404547e-01 9.938157e-01
# [9,] 0.003770551 3.815912e-06 8.822725e-05 2.136852e-09 2.944311e-05
#[10,] 0.312556568 6.042287e-01 3.060470e-06 4.765965e-03 3.928650e-03
#[11,] 0.999755008 9.512062e-01 9.979454e-01 9.598022e-02 9.805452e-01
#[12,] 0.999955201 9.968307e-01 9.763645e-01 1.501637e-01 9.809729e-01
#              [,6]         [,7]
# [1,] 8.890972e-01 7.705868e-01
# [2,] 8.924547e-01 7.810088e-01
# [3,] 7.821765e-01 2.896314e-01
# [4,] 7.191825e-01 8.159218e-01
# [5,] 8.556744e-01 5.216995e-01
# [6,] 8.261410e-01 8.401614e-01
# [7,] 8.937897e-01 7.883804e-01
# [8,] 8.822677e-01 7.503660e-01
# [9,] 7.671662e-05 6.494182e-09
#[10,] 1.609044e-04 1.275082e-03
#[11,] 8.540575e-01 4.962960e-01
#[12,] 7.961748e-01 8.213316e-01

#$mu
#           [,1]      [,2]      [,3]      [,4]      [,5]      [,6]     [,7]
# [1,] 0.2103748 0.3068145 0.3774329 0.4943874 0.6886647 0.6621627 1.279354
# [2,] 0.2103757 0.3068186 0.3774346 0.4943904 0.6886576 0.6621619 1.279365
# [3,] 0.2103690 0.3068225 0.3774324 0.4943499 0.6886940 0.6621628 1.279402
# [4,] 0.2103750 0.3068187 0.3774316 0.4943848 0.6886424 0.6621682 1.279359
# [5,] 0.2103748 0.3068188 0.3774360 0.4943829 0.6886780 0.6621701 1.279428
# [6,] 0.2103746 0.3068159 0.3774404 0.4944003 0.6887061 0.6622141 1.279352
# [7,] 0.2103750 0.3068309 0.3774369 0.4943819 0.6886722 0.6622381 1.279397
# [8,] 0.2103777 0.3068098 0.3774319 0.4943209 0.6886442 0.6621641 1.279373
# [9,] 0.2103822 0.3068732 0.3773686 0.4943849 0.6886537 0.6621572 1.279336
#[10,] 0.2103746 0.3068200 0.3774317 0.4943829 0.6886569 0.6621720 1.279369
#[11,] 0.2102877 0.3068267 0.3774361 0.4943821 0.6886438 0.6621156 1.279426
#[12,] 0.2103752 0.3068163 0.3774297 0.4943876 0.6886608 0.6621790 1.279366
#
#$sigma2
#           [,1]      [,2]      [,3]      [,4]     [,5]     [,6]     [,7]
# [1,] 0.2612167 0.4058731 0.5190057 0.6488499 1.027596 1.025466 1.828725
# [2,] 0.2612178 0.4058604 0.5190162 0.6488606 1.027590 1.025449 1.828748
# [3,] 0.2612186 0.4058527 0.5189957 0.6488799 1.027619 1.025383 1.828628
# [4,] 0.2612173 0.4058625 0.5190163 0.6488439 1.027570 1.025497 1.828765
# [5,] 0.2612122 0.4058522 0.5190132 0.6488276 1.027637 1.025431 1.828664
# [6,] 0.2612163 0.4058496 0.5190219 0.6488504 1.027544 1.025380 1.828745
# [7,] 0.2612164 0.4058902 0.5190248 0.6488332 1.027587 1.025374 1.828821
# [8,] 0.2612364 0.4058566 0.5190116 0.6490473 1.027553 1.025415 1.828788
# [9,] 0.2611864 0.4060118 0.5187161 0.6488223 1.027494 1.025036 1.828163
#[10,] 0.2612172 0.4058640 0.5190208 0.6488424 1.027539 1.025547 1.828810
#[11,] 0.2610895 0.4058716 0.5190194 0.6488254 1.027512 1.025683 1.828730
#[12,] 0.2612176 0.4058604 0.5190102 0.6488491 1.027604 1.025511 1.828782
#
#$Dispersion
#          [,1]     [,2]     [,3]     [,4]     [,5]     [,6]     [,7]
# [1,] 1.268356 1.364734 1.429030 1.380692 1.602517 1.658480 1.639113
# [2,] 1.268356 1.364673 1.429053 1.380707 1.602523 1.658454 1.639122
# [3,] 1.268399 1.364630 1.429004 1.380855 1.602489 1.658345 1.638974
# [4,] 1.268358 1.364680 1.429064 1.380686 1.602524 1.658517 1.639144
# [5,] 1.268334 1.364645 1.429040 1.380656 1.602552 1.658407 1.638978
# [6,] 1.268355 1.364648 1.429048 1.380659 1.602346 1.658222 1.639133
# [7,] 1.268354 1.364720 1.429068 1.380671 1.602487 1.658156 1.639153
# [8,] 1.268435 1.364698 1.429050 1.381288 1.602493 1.658395 1.639150
# [9,] 1.268165 1.364947 1.428467 1.380640 1.602381 1.657797 1.638629
#[10,] 1.268360 1.364679 1.429076 1.380688 1.602445 1.658589 1.639173
#[11,] 1.268252 1.364676 1.429057 1.380654 1.602430 1.658940 1.639038
#[12,] 1.268358 1.364683 1.429054 1.380690 1.602538 1.658515 1.639151
#
#$Rd
#          [,1]     [,2]     [,3]     [,4]     [,5]     [,6]     [,7]
# [1,] 1.270194 1.367272 1.431869 1.384585 1.608458 1.666049 1.645886
# [2,] 1.270194 1.367272 1.431869 1.384585 1.608458 1.666049 1.645886
# [3,] 1.270194 1.367272 1.431869 1.384585 1.608458 1.666049 1.645886
# [4,] 1.270194 1.367272 1.431869 1.384585 1.608458 1.666049 1.645886
# [5,] 1.270194 1.367272 1.431869 1.384585 1.608458 1.666049 1.645886
# [6,] 1.270194 1.367272 1.431869 1.384585 1.608458 1.666049 1.645886
# [7,] 1.270194 1.367272 1.431869 1.384585 1.608458 1.666049 1.645886
# [8,] 1.270194 1.367272 1.431869 1.384585 1.608458 1.666049 1.645886
# [9,] 1.270194 1.367272 1.431869 1.384585 1.608458 1.666049 1.645886
#[10,] 1.270194 1.367272 1.431869 1.384585 1.608458 1.666049 1.645886
#[11,] 1.270194 1.367272 1.431869 1.384585 1.608458 1.666049 1.645886
#[12,] 1.270194 1.367272 1.431869 1.384585 1.608458 1.666049 1.645886
#
#$beta1
#            [,1]       [,2]       [,3]       [,4]        [,5]       [,6]
# [1,]  1.7933703  1.4657398  1.2926576  0.9479794  0.77857492  0.8559212
# [2,]  1.7678034  1.4406421  1.2680768  0.9228747  0.75535739  0.8334167
# [3,]  4.3563763  3.9948895  3.8007554  3.4710808  3.22629469  3.2938503
# [4,] -0.8139223 -1.1079421 -1.2584319 -1.6065501 -1.70340924 -1.6163336
# [5,]  2.5535948  2.2040104  2.0172577  1.6821072  1.46485626  1.5345903
# [6,]  1.0165080  0.7094617  0.5501111  0.2024304  0.07684394  0.1603889
# [7,]  1.7048006  1.3780546  1.2058047  0.8603004  0.69403249  0.7721729
# [8,]  1.8818920  1.5530824  1.3791932  1.0363153  0.86250886  0.9383755
# [9,]  4.1877445  3.7583371  3.5191036  3.1658082  2.78574225  2.8529659
#[10,] -1.0257752 -1.2051859 -1.2793354 -1.6092491 -1.48957705 -1.3847607
#[11,]  2.3785486  2.0249722  1.8361620  1.5014883  1.27806228  1.3478802
#[12,]  1.1864346  0.8840031  0.7274491  0.3802045  0.26209539  0.3461416
#             [,7]
# [1,]  0.17618774
# [2,]  0.15437288
# [3,]  2.55435045
# [4,] -2.22806288
# [5,]  0.82041353
# [6,] -0.47989057
# [7,]  0.09426605
# [8,]  0.25750474
# [9,]  1.89498447
#[10,] -1.77056560
#[11,]  0.62734075
#[12,] -0.28743793

#$beta2
#             [,1]          [,2]          [,3]        [,4]        [,5]
# [1,] -0.02320504 -0.0074060638  0.0007633297  0.05789178  0.03404258
# [2,] -0.02578534 -0.0102015584 -0.0022356996  0.05516291  0.03064071
# [3,] -0.32720274 -0.2982417435 -0.2829614305 -0.23616593 -0.23400174
# [4,]  0.30856521  0.3082366942  0.3068650365  0.36806075  0.31655918
# [5,] -0.11386111 -0.0888446916 -0.0757214546 -0.02485662 -0.03176147
# [6,]  0.07780958  0.0838910014  0.0863687456  0.14644321  0.10609944
# [7,] -0.02937222 -0.0139917107 -0.0061391478  0.05143233  0.02644396
# [8,] -0.01701152 -0.0006673812  0.0078110603  0.06414294  0.04186330
# [9,] -0.33942968 -0.3021898857 -0.2812048799 -0.24511823 -0.21441766
#[10,]  0.59606881  0.5347233206  0.4968064245  0.55507269  0.41201526
#[11,] -0.11113040 -0.0847554584 -0.0709624471 -0.02074340 -0.02553517
#[12,]  0.07863926  0.0824681094  0.0837504891  0.14395160  0.10071051
#             [,6]        [,7]
# [1,]  0.01844058  0.08398246
# [2,]  0.01477637  0.08031332
# [3,] -0.24540312 -0.17082308
# [4,]  0.29678711  0.34706003
# [5,] -0.04426344  0.02654680
# [6,]  0.08787314  0.14504981
# [7,]  0.01051008  0.07585553
# [8,]  0.02678679  0.09227438
# [9,] -0.22240612 -0.12704655
#[10,]  0.38374371  0.36770124
#[11,] -0.03795817  0.03403702
#[12,]  0.08222618  0.13777646


