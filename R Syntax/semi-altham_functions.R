# (1)
# REPARAMETRIZATION of pmf    (i.e. function repar.semi.altham)
#
# Pr[Y_n=i]=h*exp(beta1*i+beta2*i^2)*A
# where A=1/sum(h*exp(beta1*i+beta2*i^2))
#
# in terms of (mu,sigma2)
#
# (2)
# Calculation of PMF, CDF and random number generator
# in terms of (mu, sigma2) (i.e. dsemi.altham; psemi.altham; rsemi.altham) 
# for
# Pr[Y_n=i]=h*exp(beta1*i+beta2*i^2)
# where i=0,1,...,M and parameters
#
# -oo<beta1<oo 
# -oo<beta2<oo 
#
# and hi may vary as
#
# if (hi==1) h=rep(1,M+1); 	flat
# if (hi==2)  h=log(M-X+1)+1;      decreasing
# if (hi==3) h=choose(M, X); 	convex
# if (hi==4) h=1/choose(M, X); 	concave
# if (hi==5) h=X+1; 		increasing
# if (hi==6) h=1/(X+1);         decreasing
# if (hi==7) h=M-X+1;		decreasing
# if (hi==8) h=1/(M-X+1);       increasing
# if (hi==9)  h=X*(M-X)+1;        convex
# if (hi==10)  h=1/(X*(M-X)+1);   concave
# if (hi==11)  h=log(X+1)+1;        increasing              
# if (hi==12) h=1/(log(X+1)+1);  decreasing
#
# (3)
# Compute the 1st, 2nd and 3rd moments to control (i.e. msemi.altham)
#
# (4)
# Graph of PMF for each hi with different dispersion given as
# 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0
# 

fk.exp=function (beta, hi, G, eta) 
{
# 
# functions of mu and sigma2
#   f1=sum((X-mu)*C)
#   f2=sum((X*X-sigma2-mu*mu)*C)
#   f=c(f1, f2)

X=0:M
   
if (hi==1) h=rep(1,M+1);
if (hi==2)  h=log(M-X+1)+1;
if (hi==3) h=choose(M, X);
if (hi==4) h=1/choose(M, X); 
if (hi==5) h=X+1;
if (hi==6) h=1/(X+1); 
if (hi==7) h=M-X+1;
if (hi==8) h=1/(M-X+1);
if (hi==9) h=X*(M-X)+1; 
if (hi==10)  h=1/(X*(M-X)+1);
if (hi==11)  h=log(X+1)+1;
if (hi==12) h=1/(log(X+1)+1);
 

C=h*exp(-beta%*%G[-1,]) 

# C=h*exp(beta1*i+beta2*i^2)
# where beta1*i+beta2*i^2 
# is considered as matrix multiplication
# of beta=c(beta1, beta2) and
# G=matrix(c(rep(1,M+1),X,X^2),byrow=TRUE,nrow=3) 

     
beta0=log(sum(C))

#The normalizing constant
# exp(log(sum(C)))
# 

f=(h * exp(-beta0)*exp(-beta %*% G[-1,]))%*% t(G) - eta
f

}

repar.semi.altham=function (hi, M, mu, sigma2, initial.beta) 

{
# Reparametrize semi.altham(hi, M, beta1, beta2) in terms of (mu1, mu2).
# Here (mu1, mu2) is specified, and (beta1, beta2) is to be determined.
#
# where semi.altham= Pr[Y_n=i]=h*exp(beta1*i+beta2*i^2)*A
# where A=1/sum(h*exp(beta1*i+beta2*i^2))
# i=0,1,...,M 
# and  parameters
#
# -oo<beta1<oo 
# -oo<beta2<oo 
#
# and h may vary as
#
# if (hi==1)  h=rep(1,M+1);
# if (hi==2)  h=log(M-X+1)+1;
# if (hi==3) h=choose(M, X);
# if (hi==4) h=1/choose(M, X); 
# if (hi==5) h=X+1;
# if (hi==6) h=1/(X+1);
# if (hi==7) h=M-X+1;
# if (hi==8) h=1/(M-X+1);
# if (hi==9) h=X(M-X)+1; 
# if (hi==10)  h=1/(X*(M-X)+1);
# if (hi==11)  h=log(X+1)+1;
# if (hi==12) h=1/(log(X+1)+1); 


beta=initial.beta

X=0:M

if (hi==1) h=rep(1,M+1);
if (hi==2)  h=log(M-X+1)+1;
if (hi==3) h=choose(M, X);
if (hi==4) h=1/choose(M, X); 
if (hi==5) h=X+1;
if (hi==6) h=1/(X+1);
if (hi==7) h=M-X+1;
if (hi==8) h=1/(M-X+1);
if (hi==9) h=X*(M-X)+1; 
if (hi==10)  h=1/(X*(M-X)+1);
if (hi==11)  h=log(X+1)+1;
if (hi==12) h=1/(log(X+1)+1); 


G=matrix(c(rep(1,M+1),X,X^2),byrow=TRUE,nrow=3)

## The dimension of G is 3x(M+1)
#   
#  if M=10 and X=0:M
#
# then
# G
#     [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11]
# [1,]    1    1    1    1    1    1    1    1    1     1     1
# [2,]    0    1    2    3    4    5    6    7    8     9    10
# [3,]    0    1    4    9   16   25   36   49   64    81   100


eta=c(1, mu, mu^2+sigma2)

C=h*exp(-beta%*%G[-1,]) 

# C=h*exp(beta1*i+beta2*i^2)
# where beta1*i+beta2*i^2 
# is considered as matrix multiplication
# of beta=c(beta1, beta2) and
# G=matrix(c(rep(1,M+1),X,X^2),byrow=TRUE,nrow=3) 

beta0 = log(sum(C))

# Newton method
#
   repeat {
        beta.old = beta
	C.old=h * exp(-beta.old %*% G[-1, ])
        beta0 = log(sum(C.old))
        f.old = fk.exp(beta.old, hi = hi, G = G, eta = eta)
        dev.ent = c(h * exp(-beta0) * exp(-beta.old %*% G[-1, 
            ]))
        cov.ent = cov.wt(t(G), wt = dev.ent, method = "ML", cor = FALSE)
  
   #    J= [df1/dbeta1  df1/dbeta2]
   #       [df2/dbeta1  df2/dbeta2]
   #  
   #  is weighted variance-covariance matrix

        hess.ent = cov.ent$cov
  # Hessian matrix

        Hess = hess.ent[-1, -1]
        inv.Hess = solve(Hess)
        beta = beta.old + f.old[, -1] %*% inv.Hess
        
        if (max(abs(beta - beta.old)) < 1e-08) 
            break
    }
    out =list(estimates = beta, infmat = inv.Hess)
    beta = out$estimates
    # beta0 = log(sum(h * exp(-beta %*% G[-1, ])))
   
    pmf.semi.altham = (h * exp(-beta0) * exp(-beta %*% G[-1, ]))
   
     #Pr[Y_n=i]=h*exp(beta1*i+beta2*i^2)

    list(Betas=beta , Estimates=pmf.semi.altham)
}


# TEST

hi=1; M=10; mu=5; sigma2=5; initial.beta=c(0,0)
out=repar.semi.altham(hi, M, mu, sigma2, initial.beta)

#$Betas
#           [,1]       [,2]
#[1,] -0.8665593 0.08665593
#
#$Estimates
#           [,1]       [,2]       [,3]      [,4]      [,5]      [,6]      [,7]
#[1,] 0.01944138 0.04240674 0.07778113 0.1199624 0.1555779 0.1696610 0.1555779
#          [,8]       [,9]      [,10]      [,11]
#[1,] 0.1199624 0.07778113 0.04240674 0.01944138

################################################################################
################################################################################
### PMF ###

dsemi.altham=function(hi, M, beta)
{
# PMF for Pr[Y_n=i]=h*exp(beta1*i+beta2*i^2)*A
# where A=1/sum(h*exp(beta1*i+beta2*i^2))
# i=0,1,...,M 
# and  parameters
#
# -oo<beta1<oo 
# -oo<beta2<oo 
#
# and h may vary as
#
# if (hi==1) h=rep(1,M+1); 	flat
# if (hi=2)  h=log(M-i+1)+1;      decreasing
# if (hi==3) h=choose(M, X); 	convex
# if (hi==4) h=1/choose(M, X); 	concave
# if (hi==5) h=X+1; 		increasing
# if (hi=6) h=1/(i+1);         decreasing
# if (hi==7) h=M-X+1;		decreasing
# if (hi==8) h=1/(M-X+1);	increasing
# if (hi=9)  h=X(M-X)+1;        convex
# if (hi=10)  h=1/(X*(M-i)+1);   concave
# if (hi=11)  h=log(i+1)+1;        increasing              
# if (hi=12) h=1/(log(i+1)+1);  decreasing


X=0:M

if (hi==1) h=rep(1,M+1);
if (hi==2)  h=log(M-X+1)+1;
if (hi==3) h=choose(M, X);
if (hi==4) h=1/choose(M, X);
if (hi==5) h=X+1;
if (hi==6) h=1/(X+1); 
if (hi==7) h=M-X+1;
if (hi==8) h=1/(M-X+1); 
if (hi==9) h=X*(M-X)+1; 
if (hi==10)  h=1/(X*(M-X)+1);
if (hi==11)  h=log(X+1)+1;
if (hi==12) h=1/(log(X+1)+1); 


G=matrix(c(rep(1,M+1),X,X^2),byrow=TRUE,nrow=3)

beta0 = log(sum(h * exp(-beta %*% G[-1, ])))
   
pmf.semi.altham = (h * exp(-beta0) * exp(-beta %*% G[-1, ]))

pmf.semi.altham
}

# Test

hi=7; M=10; beta=c(-0.8665593, 0.08665593)
dsemi.altham(hi, M, beta)


#[1,] 0.03888277 0.07633213 0.1244498 0.1679474 0.1866934 0.1696610 0.1244623
#           [,8]       [,9]       [,10] [,11]
#[1,] 0.07197745 0.03111245 0.008481348     0


##########################################################################
##########################################################################
##### CDF ###### 

psemi.altham=function(hi, M, beta)
{
# Cumulative probability distribution for
# p(i)~ h(i)exp(beta1*i+beta2*i^2)
# where h(i) can be
# if (hi==1) h=rep(1,M+1); 	flat
# if (hi=2)  h=log(M-i+1)+1;      decreasing
# if (hi==3) h=choose(M, X); 	convex
# if (hi==4) h=1/choose(M, X); 	concave
# if (hi==5) h=X+1; 		increasing
# if (hi=6) h=1/(i+1);         decreasing
# if (hi==7) h=M-X+1;		decreasing
# if (hi==8) h=1/(M-X+1);      increasing
# if (hi=9)  h=X(M-X)+1;        convex
# if (hi=10)  h=1/(X*(M-i)+1);   concave
# if (hi=11)  h=log(i+1)+1;        increasing              
# if (hi=12) h=1/(log(i+1)+1);  decreasing

X=0:M

P=dsemi.altham(hi, M, beta)
  CDF=cumsum(P)
  return(CDF)
}

# Test

hi=7; M=10; beta=c(-0.8665593, 0.08665593)
psemi.altham(hi, M, beta)

# [1] 0.03888277 0.11521490 0.23966470 0.40761208 0.59430550 0.76396647
# [7] 0.88842875 0.96040620 0.99151865 1.00000000 1.00000000



##########################################################################
##########################################################################
##### Random Numbers ###### 

rsemi.altham=function(hi, k, M, beta)

{
# Generate k random numbers from 
# the  
# which has the support {0, 1, ..., n}, and parameters beta1 and beta2
# in terms of (mu, sigma2)
# 
# p(i)~ h(i)exp(beta1*i+beta2*i^2)
# where h(i) can be
# if (hi==1) h=rep(1,M+1); 	flat
# if (hi=2)  h=log(M-i+1)+1;      decreasing
# if (hi==3) h=choose(M, X); 	convex
# if (hi==4) h=1/choose(M, X); 	concave
# if (hi==5) h=X+1; 		increasing
# if (hi=6) h=1/(i+1);         decreasing
# if (hi==7) h=M-X+1;		decreasing
# if (hi==8) h=1/(M-X+1);      increasing
# if (hi=9)  h=X(M-X)+1;        convex
# if (hi=10)  h=1/(X*(M-i)+1);   concave
# if (hi=11)  h=log(i+1)+1;        increasing              
# if (hi=12) h=1/(log(i+1)+1);  decreasing



  P=psemi.altham(hi, M, beta)
  x.unif=runif(k)

  x=NULL
  for(i in 1:k)
  {
    tmp=P-x.unif[i]
    flag=(tmp >= 0)
    x.tmp=min((1:(k+1))[flag])-1
    x=c(x, x.tmp)
  }

  return(x)
}

#Test

k=500; 
hi=7; M=10; beta=c(-0.8665593, 0.08665593)

table(rsemi.altham(hi, k, M, beta))

# 0  1  2  3  4  5  6  7  8  9 
# 24 38 54 77 97 91 66 36 14  3 



#######################################################################
#######################################################################
### Compute the 1st, 2nd and 3rd moments ###

msemi.altham=function(hi, M, beta)
{
# Compute the 1st, 2nd and 3rd moments for 
# Pr[Y_n=i]=h*exp(beta1*i+beta2*i^2)
# where i=0,1,...,M and parameters
#
# -oo<beta1<oo 
# -oo<beta2<oo 
#
# and hi may vary as
#
# if (hi==1) h=rep(1,M+1); 	flat
# if (hi=2)  h=log(M-i+1)+1;      decreasing
# if (hi==3) h=choose(M, X); 	convex
# if (hi==4) h=1/choose(M, X); 	concave
# if (hi==5) h=X+1; 		increasing
# if (hi=6) h=1/(i+1);         decreasing
# if (hi==7) h=M-X+1;		decreasing
# if (hi==8) h=1/(M-X+1);      increasing
# if (hi=9)  h=X(M-X)+1;        convex
# if (hi=10)  h=1/(X*(M-i)+1);   concave
# if (hi=11)  h=log(i+1)+1;        increasing        
# if (hi=12) h=1/(log(i+1)+1);  decreasing


#
#  in terms of new parameters (mu, sigma2)

  PMF=dsemi.altham(hi, M, beta)
  m1=sum((0:M)*PMF)
  m2=sum((0:M)^2*PMF)
  m3=sum((0:M)^3*PMF)
  sigma2=m2-m1^2
  dispersion=sigma2/(m1*(1-m1/M))
  return(c(m1, m2, m3, sigma2, dispersion))

}

# Test

hi=7; M=10; mu=5; sigma2=5; initial.beta=c(0,0)
out=repar.semi.altham(hi, M, mu, sigma2, initial.beta)
out
# $Betas
#           [,1]       [,2]
# [1,] -0.7220036 0.04512046
#
# $Estimates
#           [,1]       [,2]       [,3]      [,4]      [,5]      [,6]      [,7]
# [1,] 0.02677724 0.04742147 0.07578756 0.1089407 0.1401641 0.1601982 0.1606048
#           [,8]      [,9]      [,10] [,11]
# [1,] 0.1379239 0.0962004 0.04598164     0



hi=7; M=10; beta=c(-0.7220036, 0.04512046)

msemi.altham(hi, M, beta)
#
#[1]   5.0000  30.0000 197.3642   5.0000   2.0000

###########################################################
mle.optim.semi.altham=function(hi, y, M)
{
#
# Pr[Y_n=i]=h*exp(beta1*i+beta2*i^2)*A
# where A=1/exp(log(sum(h*exp(beta1*i+beta2*i^2))))
#
# in terms of (mu,sigma2)
#
# where i=0,1,...,M and parameters
#
# -oo<beta1<oo 
# -oo<beta2<oo 
#
# and hi may vary as
#
# if (hi==1) h=rep(1,M+1); 	flat
# if (hi==2)  h=log(M-X+1)+1;      decreasing
# if (hi==3) h=choose(M, X); 	convex
# if (hi==4) h=1/choose(M, X); 	concave
# if (hi==5) h=X+1; 		increasing
# if (hi==6) h=1/(X+1);         decreasing
# if (hi==7) h=M-X+1;		decreasing
# if (hi==8) h=1/(M-X+1);	 increasing
# if (hi==9)  h=X*(M-X)+1;        convex
# if (hi==10)  h=1/(X*(M-X)+1);   concave
# if (hi==11)  h=log(X+1)+1;        increasing              
# if (hi==12) h=1/(log(X+1)+1);  decreasing

 X=0:M

if (hi==1) h=rep(1,M+1);
if (hi==2)  h=log(M-X+1)+1;
if (hi==3) h=choose(M, X);
if (hi==4) h=1/choose(M, X); 
if (hi==5) h=X+1;
if (hi==6) h=1/(X+1);
if (hi==7) h=M-X+1;
if (hi==8) h=1/(M-X+1);
if (hi==9) h=X*(M-X)+1; 
if (hi==10)  h=1/(X*(M-X)+1);
if (hi==11)  h=log(X+1)+1;
if (hi==12) h=1/(log(X+1)+1); 



# Count frequencies.
  TMP=sort(unique(y))
  Freq=rep(0, length(TMP))
  for(i in 1:length(TMP))
  {
    flag=(y==TMP[i])
    Freq[i]=sum(flag)
  }


# 1. Define the -loglikelyhood function to be minimized.
#
  nlk.semi.altham=function(para)
  {
# Parameter ranges:
#
#  beta1=para[1]: -oo < beta1 < oo.
#  beta2=para[2]: -oo < beta2 < oo.

#
    print(paste("Parameter values:"))
    print(para)
   
# a1=beta1; a2=beta2
a1=para[1]; a2=para[2]
a=c(a1,a2)
   
      s=sum(Freq*log(dsemi.altham(hi, M, a)[TMP+1]))
      s=-s
      names(s)="negative loglikelihood"
    
    print(s)
    return(s)
  }

# 2. Obtain the MLE by minimizing the -loglikelyhood function.
#    The initial values of parameters are set arbitrary

#  para0=c(a1.hat, a2.hat)
a1.hat=0; a2.hat=0;

para0=c(a1.hat, a2.hat)
 soln=optim(para0, nlk.semi.altham, gr=NULL, method="BFGS", hessian=TRUE,
#  soln=optim(para0, nlk.semi.altham, gr=NULL, method="CG", hessian=TRUE,
#  soln=optim(para0, nlk.semi.altham, gr=NULL, method="SANN", hessian=TRUE,
#  soln=optim(para0, nlk.semi.altham, gr=NULL, method="Nelder-Mead", hessian=TRUE,
             control=list(trace=100, maxit=1000))
#             control=list(trace=100, ndeps=1e-5, maxit=1000))
#  soln=nlm(nlk.semi.altham, para0, hessian=TRUE, print.level=2)
  print(soln)
    
# 3. Obtain the value of -loglikelyhood at MLE, MLE and the SE.
#
# -log-likelihood
        nlk.smallest=soln$value               # from optim()
#        nlk.smallest=soln$minimum            # from nlm()
# MLEs
        MLE=soln$par                          # from optim()
#        MLE=soln$estimate                    # from nlm()
# AIC
        AIC=2*nlk.smallest+2*length(para0)
        print(paste("AIC =", AIC))
# Asymptotic covariance matrix
        hess=solve(soln$hessian)
# Estimated standard erros
        SE=sqrt(diag(hess))

    return(soln, nlk.smallest, AIC, MLE, SE)

}


# Test:
# 
# hi=1; M=10; mu=5; sigma2=5; initial.beta=c(0,0)
# out=repar.semi.altham(hi, M, mu, sigma2, initial.beta)
# beta1=out$Betas[1]; beta2=out$Betas[2]; beta=c(beta1, beta2);
# > beta
# [1] -0.86655933  0.08665593

#
# k=100; 
#

#
# y=rsemi.altham(hi, k, M, beta=beta)
# result=mle.optim.semi.altham(hi, y, M)
#
# [1] "Parameter values:"
# [1] -0.87208375  0.08168664
# negative loglikelihood 
#              221.7115 
# $par
# [1] -0.87208375  0.08368664
#
# $value
# [1] 221.5948
#
# $counts
# function gradient 
#      43       10 
#
# $convergence
# [1] 0
#
# $message
# NULL
# 
# $hessian
#          [,1]     [,2]
# [1,]  510.5601  5242.38
# [2,] 5242.3801 57714.72
#
# [1] "AIC = 447.189674619331"



###########################################################
###########################################################
mle.optim.semi.altham.freq=function(hi, Freq, M)
{
#
# Pr[Y_n=i]=h*exp(beta1*i+beta2*i^2)*A
# where A=1/exp(log(sum(h*exp(beta1*i+beta2*i^2))))
#
# in terms of (mu,sigma2)
#
# where i=0,1,...,M and parameters
#
# -oo<beta1<oo 
# -oo<beta2<oo 
#
# and hi may vary as
#
# if (hi==1) h=rep(1,M+1); 	flat
# if (hi==2)  h=log(M-X+1)+1;      decreasing
# if (hi==3) h=choose(M, X); 	convex
# if (hi==4) h=1/choose(M, X); 	concave
# if (hi==5) h=X+1; 		increasing
# if (hi==6) h=1/(X+1);         decreasing
# if (hi==7) h=M-X+1;		decreasing
# if (hi==8) h=1/(M-X+1);	 increasing
# if (hi==9)  h=X*(M-X)+1;        convex
# if (hi==10)  h=1/(X*(M-X)+1);   concave
# if (hi==11)  h=log(X+1)+1;        increasing              
# if (hi==12) h=1/(log(X+1)+1);  decreasing

 
 X=0:M

if (hi==1) h=rep(1,M+1);
if (hi==2)  h=log(M-X+1)+1;
if (hi==3) h=choose(M, X);
if (hi==4) h=1/choose(M, X); 
if (hi==5) h=X+1;
if (hi==6) h=1/(X+1);
if (hi==7) h=M-X+1;
if (hi==8) h=1/(M-X+1);
if (hi==9) h=X*(M-X)+1; 
if (hi==10)  h=1/(X*(M-X)+1);
if (hi==11)  h=log(X+1)+1;
if (hi==12) h=1/(log(X+1)+1); 

# 1. Define the -loglikelyhood function to be minimized.
#
  nlk.semi.altham=function(para)
  {
# Parameter ranges:
#
#  beta1=para[1]: -oo < beta1 < oo.
#  beta2=para[2]: -oo < beta2 < oo.

#
    print(paste("Parameter values:"))
    print(para)
   
# a1=beta1; a2=beta2
a1=para[1]; a2=para[2]
a=c(a1,a2)
   
       s=sum(Freq*log(dsemi.altham(hi, M, a)))
      s=-s
      names(s)="negative loglikelihood"
    
    print(s)
    return(s)
  }

# 2. Obtain the MLE by minimizing the -loglikelyhood function.
#    The initial values of parameters are set arbitrary

#  para0=c(a1.hat, a2.hat)
a1.hat=0; a2.hat=0;

para0=c(a1.hat, a2.hat)
 soln=optim(para0, nlk.semi.altham, gr=NULL, method="BFGS", hessian=TRUE,
#  soln=optim(para0, nlk.semi.altham, gr=NULL, method="CG", hessian=TRUE,
#  soln=optim(para0, nlk.semi.altham, gr=NULL, method="SANN", hessian=TRUE,
#  soln=optim(para0, nlk.semi.altham, gr=NULL, method="Nelder-Mead", hessian=TRUE,
             control=list(trace=100, maxit=1000))
#             control=list(trace=100, ndeps=1e-5, maxit=1000))
#  soln=nlm(nlk.semi.altham, para0, hessian=TRUE, print.level=2)
  print(soln)
    
# 3. Obtain the value of -loglikelyhood at MLE, MLE and the SE.
#
# -log-likelihood
        nlk.smallest=soln$value               # from optim()
#        nlk.smallest=soln$minimum            # from nlm()
# MLEs
        MLE=soln$par                          # from optim()
#        MLE=soln$estimate                    # from nlm()
# AIC
        AIC=2*nlk.smallest+2*length(para0)
        print(paste("AIC =", AIC))
# Asymptotic covariance matrix
        hess=solve(soln$hessian)
# Estimated standard erros
        SE=sqrt(diag(hess))

    return(soln, nlk.smallest, AIC, MLE, SE)

}




