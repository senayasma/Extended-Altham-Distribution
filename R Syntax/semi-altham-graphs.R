########################################
# GRAPHS FOR 
# Pr[Y_n=i]=h*exp(beta1*i+beta2*i^2)*A
# where A=1/sum(h*exp(beta1*i+beta2*i^2))
#
# in terms of (mu,sigma2)
# where i=0,1,...,M and parameters
#
# -oo<beta1<oo 
# -oo<beta2<oo 
#
# and hi may vary as
#
# if (hi==1) h=rep(1,M+1); 	flat
# if (hi==2)  h=log(M-X+1);      decreasing
# if (hi==3) h=choose(M, X); 	convex
# if (hi==4) h=1/choose(M, X); 	concave
# if (hi==5) h=X+1; 		increasing
# if (hi==6) h=1/(X+1);         decreasing
# if (hi==7) h=M-X+1;		decreasing
# if (hi==8) h=1/(M-X+1);   	 increasing
# if (hi==9)  h=X*(M-X)+1;        convex
# if (hi==10)  h=1/(X*(M-X)+1);   concave
# if (hi==11)  h=log(X+1);        increasing              
# if (hi==12) h=1/(log(X+1)+1);  decreasing

########################################
# 1. Create pmf graphs of the semi.altham(mu, sigma2) of various dispersion.
#
# Dispersion for discrete distribution on the support {0, 1, ..., M}
# is defined as
#   D=variance/[mean*(1-mean/M)]
# For binomial(M, p), the dispersion is always 1:
#   D=M*p*(1-p)/[M*p*(1-M*p/M)]=1
# For semi.altham(mu, sigma2), the dispersion is
#   0 < D=sigma2/[mu*(1-mu/M)] < M
#

# A. 
# GRAPH FOR PMF hi=1-6

source("semi.altham_functions.r")

pdf("semi.altham-pmf1-6.pdf")

M=10; mu=5;

par(mfrow=c(2,3))

Support=0:M

D.tmp=seq(0.10, 1, length.out=10)
D.tmp=c(D.tmp, 2:9)

# [1] 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0

for (hi in 1:6)
{

if (hi==1) h="rep(1,M+1)";
if (hi==2)  h="log(M-X+1)+1";
if (hi==3) h="choose(M, X)";
if (hi==4) h="1/choose(M, X)"; 
if (hi==5) h="X+1";
if (hi==6) h="1/(X+1)";
if (hi==7) h="M-X+1";
if (hi==8) h="1/(M-X+1)";
if (hi==9) h="X*(M-X)+1"; 
if (hi==10)  h="1/(X*(M-X)+1)";
if (hi==11)  h="log(X+1)+1"; 
if (hi==12) h="1/(log(X+1)+1)"; 


plot(Support, seq(0, 0.8, length.out=M+1), 
     xlab="X", ylab="Probability", type="n", main=h)
abline(v=Support, lty=3, col=8)

initial.beta=c(0,0)  
i=0
for(D in D.tmp)
{ 
 i=i+1
   
 sigma2=D*mu*(1-mu/M)
 out=repar.semi.altham(hi, M, mu, sigma2, initial.beta=initial.beta) 
 initial.beta=out$Betas

  PMF=dsemi.altham(hi, M, beta=initial.beta)
  LTY=i
  if(D==1) i=1
  COLOR="red"
  if(D<1) COLOR="green"
  if(D>1) COLOR="blue"
  lines(Support, PMF, col=COLOR, lty=LTY)

 }

}

dev.off()


########################################
# B. 
# GRAPH FOR PMF hi=7-12


pdf("semi.altham-pmf7-12.pdf")

M=10; mu=5;

par(mfrow=c(2,3))

Support=0:M

D.tmp=seq(0.10, 1, length.out=10)
D.tmp=c(D.tmp, 2:9)

# [1] 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0

for (hi in 7:12)
{
if (hi==1) h="rep(1,M+1)";
if (hi==2)  h="log(M-X+1)+1";
if (hi==3) h="choose(M, X)";
if (hi==4) h="1/choose(M, X)"; 
if (hi==5) h="X+1";
if (hi==6) h="1/(X+1)"; 
if (hi==7) h="M-X+1";
if (hi==8) h="1/(M-X+1)";
if (hi==9) h="X*(M-X)+1"; 
if (hi==10)  h="1/(X*(M-X)+1)";
if (hi==11)  h="log(X+1)+1";
if (hi==12) h="1/(log(X+1)+1)"; 

plot(Support, seq(0, 0.8, length.out=M+1), 
     xlab="X", ylab="Probability", type="n", main=h)
abline(v=Support, lty=3, col=8)

initial.beta=c(0,0)  
i=0
for(D in D.tmp)
{ 
 i=i+1
 sigma2=D*mu*(1-mu/M)
 out=repar.semi.altham(hi, M, mu, sigma2, initial.beta=initial.beta) 
 initial.beta=out$Betas
  PMF=dsemi.altham(hi, M, beta=initial.beta)
  LTY=i
  if(D==1) i=1
  COLOR="red"
  if(D<1) COLOR="green"
  if(D>1) COLOR="blue"
  lines(Support, PMF, col=COLOR, lty=LTY)

 }

}

dev.off()


# 2. Create the graph of comparison between the equally-dispersed semi.altham 
#    and the binomial of the same mean.  
#

# 2A. comparison-semi.altham-binomial hi=1-6

#postscript("comparison-semi.altham-binomial.ps")
pdf("comparison-semi.altham-binomial-1-6.pdf")

M=40
Support=0:M

#par(mfrow=c(1,1))
par(mfrow=c(2,3))

MEANS=seq(1, M-1, length.out=as.integer(M/3))

#plot(Support, seq(0, 0.4, length.out=M+1), 
#     xlab="X", ylab="Probability", type="n")
# abline(v=Support, lty=3, col=8)

for(hi in 1:6)
{

if (hi==1) h="rep(1,M+1)";
if (hi==2)  h="log(M-X+1)+1";
if (hi==3) h="choose(M, X)";
if (hi==4) h="1/choose(M, X)"; 
if (hi==5) h="X+1";
if (hi==6) h="1/(X+1)"; 
if (hi==7) h="M-X+1";
if (hi==8) h="1/(M-X+1)";
if (hi==9) h="X*(M-X)+1"; 
if (hi==10)  h="1/(X*(M-X)+1)";
if (hi==11)  h="log(X+1)+1";
if (hi==12) h="1/(log(X+1)+1)"; 


plot(Support, seq(0, 0.4, length.out=M+1), 
     xlab="X", ylab="Probability", type="n", main=h)
 abline(v=Support, lty=3, col=8)

 for(mu in MEANS)
 {
  sigma2=mu*(1-mu/M)
#  print(paste("mu =", mu, " sigma2 =", sigma2))

# initial values

#################### hi=1 h=1 ######################################

  if(hi==1 && mu==1) initial.beta=c(0,0);
  if(hi==1 && mu>1 && mu<4.5) initial.beta=c(-1.017185, 0.1290705);  
  if(hi==1 && mu>7 && mu<8) initial.beta=c(-1.199233, 0.08574957);
  if(hi==1 && mu>10 && mu<=11) initial.beta=c(-1.287356, 0.07153209);
  if(hi==1 && mu>13 && mu<=14) initial.beta=c(-1.428352, 0.05951514);
  if(hi==1 && mu>16 && mu<=17) initial.beta=c(-1.599992, 0.05333307);
  if(hi==1 && mu==20) initial.beta=c(-1.904762, 0.05012531);
  if(hi==1 && mu>22 && mu<=25) initial.beta=c(-2.105263, 0.05012531);
  if(hi==1 && mu>25 && mu<=27) initial.beta=c(-2.666654, 0.05333307);
  if(hi==1 && mu>27 && mu<=30) initial.beta=c(-3.332859, 0.05951514);
  if(hi==1 && mu>30 && mu<=33) initial.beta=c(-4.977422, 0.07776383);
  if(hi==1 && mu>33 && mu<=36) initial.beta=c(-7.70935, 0.1100388);
  if(hi==1 && mu==39) initial.beta=c(-15.02605, 0.196416);

################ hi=3 h=choose(M, X); #################################

  if(hi==3 && mu==1) initial.beta=c(3.637586, 2.079528e-16);
  if(hi==3 && mu>1 && mu<4.5) initial.beta=c(2.197225, 2.572754e-16);  
  if(hi==3 && mu>7 && mu<8) initial.beta=c(1.734601, -1.597501e-15);
  if(hi==3 && mu>10 && mu<=11) initial.beta=c(1.098612, -7.651001e-16);
  if(hi==3 && mu>13 && mu<=14) initial.beta=c(0.8472979, -1.672698e-16);
  if(hi==3 && mu>16 && mu<=17) initial.beta=c(0.4054651, 9.438626e-18);
  if(hi==3 && mu==20) initial.beta=c(0.1000835, 1.377368e-16);
  if(hi==3 && mu>22 && mu<=25) initial.beta=c(-0.3022809, -1.860093e-15);
  if(hi==3 && mu>25 && mu<=27) initial.beta=c(-0.6190392, 2.340299e-16);
  if(hi==3 && mu>27 && mu<=30) initial.beta=c(-0.9694006, -9.076763e-15);
  if(hi==3 && mu>30 && mu<=33) initial.beta=c(-1.386294, -3.049898e-14);
  if(hi==3 && mu>33 && mu<=36) initial.beta=c(-1.94591, -9.775039e-15);
  if(hi==3 && mu==39) initial.beta=c(-3.245193, 3.818461e-13);

#################### hi=5  h=X+1; ########################################

  if(hi==5 && mu==1) initial.beta=c(0.4487824, 0.167705);
  if(hi==5 && mu>4 && mu<4.5) initial.beta=c(-0.5765087, 0.1034802);
  if(hi==5 && mu>7.2 && mu<=7.5) initial.beta=c(-0.8717244, 0.07877348);
  if(hi==5 && mu>10 && mu<=11) initial.beta=c(-1.131407, 0.06147013);
  if(hi==5 && mu>13 && mu<=14) initial.beta=c(-1.362723, 0.05315132);
  if(hi==5 && mu>16 && mu<=17) initial.beta=c(-1.583141, 0.04976643);
  if(hi==5 && mu==20) initial.beta=c(-1.851086, 0.04874663);
  if(hi==5 && mu>22 && mu<=25) initial.beta=c(-2.268400, 0.05023506);
  if(hi==5 && mu>25 && mu<=27) initial.beta=c(-2.782498, 0.05423105);
  if(hi==5 && mu>27 && mu<=30) initial.beta=c(-3.568240, 0.06210076);
  if(hi==5 && mu>30 && mu<=33) initial.beta=c(-5.238844, 0.08104664);
  if(hi==5 && mu>33 && mu<=36) initial.beta=c(-8.387756, 0.1183825);
  if(hi==5 && mu==39) initial.beta=c(-19.81765, 0.2514471);

################ hi=7 h=M-X+1; ##########################################

  if(hi==7 && mu==1) initial.beta=c(-0.2981152, 0.2514471);
  if(hi==7 && mu>1 && mu<4.5) initial.beta=c(-1.041274, 0.1286951);  
  if(hi==7 && mu>7 && mu<8) initial.beta=c(-1.222602, 0.08530752);
  if(hi==7 && mu>10 && mu<=11) initial.beta=c(-1.353852, 0.06607121);
  if(hi==7 && mu>13 && mu<=14) initial.beta=c(-1.527676, 0.05521734);
  if(hi==7 && mu>16 && mu<=17) initial.beta=c(-1.714714, 0.05070316);
  if(hi==7 && mu==20) initial.beta=c(-1.953599, 0.04887204);
  if(hi==7 && mu>22 && mu<=25) initial.beta=c(-2.331937, 0.04944631);
  if(hi==7 && mu>25 && mu<=27) initial.beta=c(-2.794208, 0.05239318);
  if(hi==7 && mu>27 && mu<=30) initial.beta=c(-3.479322, 0.05845432);
  if(hi==7 && mu>30 && mu<=33) initial.beta=c(-4.838769, 0.0724201);
  if(hi==7 && mu>33 && mu<=36) initial.beta=c(-7.121272, 0.09716281);
  if(hi==7 && mu==39)initial.beta=c(-13.86518, 0.167705);

################ hi=4 h=1/choose(M, X); ##########################################


  if(hi==4 && mu==1) initial.beta=c(-4.303162, 0.6112999);
  if(hi==4 && mu>1 && mu<4.5) initial.beta=c(-4.314866, 0.2677898);  
  if(hi==4 && mu>7 && mu<8) initial.beta=c(-4.057696, 0.1787719);
  if(hi==4 && mu>10 && mu<=11) initial.beta=c(-3.850276, 0.1373737);
  if(hi==4 && mu>13 && mu<=14) initial.beta=c(-3.740427, 0.1156765);
  if(hi==4 && mu>16 && mu<=17) initial.beta=c(-3.754178, 0.1046272);
  if(hi==4 && mu==20) initial.beta=c(-3.954904, 0.1001240);
  if(hi==4 && mu>22 && mu<=25) initial.beta=c(-4.415353, 0.1025662);
  if(hi==4 && mu>25 && mu<=27) initial.beta=c(-5.154624, 0.1110587);
  if(hi==4 && mu>27 && mu<=30) initial.beta=c(-6.484727, 0.1285782);
  if(hi==4 && mu>30 && mu<=33) initial.beta=c(-8.960896, 0.1617539);
  if(hi==4 && mu>33 && mu<=36) initial.beta=c(-14.10145, 0.2291585);
  if(hi==4 && mu==39) initial.beta=c(-32.33174, 0.4576186);


################ hi=9 h=X*(M-X)+1; ##########################################


  if(hi==9 && mu==1) initial.beta=c(2.269248, -0.05154036);
  if(hi==9 && mu>1 && mu<4.5) initial.beta=c(-0.04780887, 0.06719773);  
  if(hi==9 && mu>7 && mu<8) initial.beta=c(-0.8405965, 0.06993819);
  if(hi==9 && mu>10 && mu<=11) initial.beta=c(-1.113251, 0.0594664);
  if(hi==9 && mu>13 && mu<=14) initial.beta=c(-1.326611, 0.05271925);
  if(hi==9 && mu>16 && mu<=17) initial.beta=c(-1.543899, 0.04895616);
  if(hi==9 && mu==20) initial.beta=c(-1.795508, 0.04739314);
  if(hi==9 && mu>22 && mu<=25) initial.beta=c(-2.2346, 0.04821576);
  if(hi==9 && mu>25 && mu<=27) initial.beta=c(-2.697669, 0.05118175);
  if(hi==9 && mu>27 && mu<=30) initial.beta=c(-3.359099, 0.05682025);
  if(hi==9 && mu>30 && mu<=33) initial.beta=c(-4.340993, 0.0660842);
  if(hi==9 && mu>33 && mu<=36) initial.beta=c(-5.567611, 0.07615741);
  if(hi==9 && mu==39) initial.beta=c(1.85398, -0.05154036);

######################################

################ hi=10 h=1/(X*(M-X)+1); ##########################################


  if(hi==10 && mu==1) next; #initial.beta=c(-4.035039, 0.9901828);
  if(hi==10 && mu>1 && mu<4.5) initial.beta=c(-2.510631, 0.2594033);  
  if(hi==10 && mu>7 && mu<8) initial.beta=c(-1.916851, 0.1245232);
  if(hi==10 && mu>10 && mu<=11) initial.beta=c(-1.680665, 0.07978527);
  if(hi==10 && mu>13 && mu<=14) initial.beta=c(-1.659662, 0.06209267);
  if(hi==10 && mu>16 && mu<=17) initial.beta=c(-1.791217, 0.05526399);
  if(hi==10 && mu==20) initial.beta=c(-2.059873, 0.05274772);
  if(hi==10 && mu>22 && mu<=25) initial.beta=c(-2.472172, 0.05410599);
  if(hi==10 && mu>25 && mu<=27) initial.beta=c(-3.035043, 0.05907201);
  if(hi==10 && mu>27 && mu<=30) initial.beta=c(-4.106849, 0.07198864);
  if(hi==10 && mu>30 && mu<=33) initial.beta=c(-6.565777, 0.1046964);
  if(hi==10 && mu>33 && mu<=36) initial.beta=c(-13.29024, 0.1942729);
  if(hi==10 && mu==39) next;
  
######################################

################ hi=11 h=log(X+1)+1; ##########################################


  if(hi==11 && mu==1) initial.beta=c(0.2635505, 0.1722706);
  if(hi==11 && mu>2 && mu<4.5) initial.beta=c(-0.7727501, 0.1118146);  
  if(hi==11 && mu>7 && mu<8) initial.beta=c(-1.074471, 0.08038401);
  if(hi==11 && mu>10 && mu<=11) initial.beta=c(-1.257837, 0.06440741);
  if(hi==11 && mu>13 && mu<=14) initial.beta=c(-1.461840, 0.0548908);
  if(hi==11 && mu>16 && mu<=17) initial.beta=c(-1.666313, 0.0509663);
  if(hi==11 && mu==20) initial.beta=c(-1.922688, 0.04962248);
  if(hi==11 && mu>22 && mu<=25) initial.beta=c(-2.262122, 0.05050303);
  if(hi==11 && mu>25 && mu<=27) initial.beta=c(-2.837562, 0.05473732);
  if(hi==11 && mu>27 && mu<=30) initial.beta=c(-3.618037, 0.06251163);
  if(hi==11 && mu>30 && mu<=33) initial.beta=c(-5.283736, 0.0813781);
  if(hi==11 && mu>33 && mu<=36) initial.beta=c(-8.42932, 0.1186652);
  if(hi==11 && mu==39) initial.beta=c(-19.85621, 0.2516892);
 

################# hi=2 h=log(M-X+1)##################################

  if(hi==2 && mu==1) initial.beta=c(-0.2803165, 0.2516607);
  if(hi==2 && mu>2 && mu<4.5) initial.beta=c(-1.023629, 0.128937);  
  if(hi==2 && mu>7 && mu<8) initial.beta=c(-1.205370, 0.08558768);
  if(hi==2 && mu>10 && mu<=11) initial.beta=c(-1.33747, 0.06640245);
  if(hi==2 && mu>13 && mu<=14) initial.beta=c(-1.485527, 0.05671613);
  if(hi==2 && mu>16 && mu<=17) initial.beta=c(-1.66826, 0.05173486);
  if(hi==2 && mu==20) initial.beta=c(-1.946819, 0.04951499);
  if(hi==2 && mu>22 && mu<=25) initial.beta=c( -2.334893, 0.05031695);
  if(hi==2 && mu>25 && mu<=27) initial.beta=c(-2.811033, 0.05354517);
  if(hi==2 && mu>27 && mu<=30) initial.beta=c(-3.513287, 0.05991752);
  if(hi==2 && mu>30 && mu<=33) initial.beta=c(-4.595194, 0.07080732);
  if(hi==2 && mu>33 && mu<=36) initial.beta=c( -6.038489, 0.08411448);
  if(hi==2 && mu==39)initial.beta=c(-11.68763, 0.1494534);

####################hi=6 h=1/(X+1) ##############################

  if(hi==6 && mu==1) initial.beta=c(-1.001554, 0.3380987);
  if(hi==6 && mu>2 && mu<4.5) initial.beta=c(-1.472866, 0.1563312);  
  if(hi==6 && mu>7 && mu<8) initial.beta=c(-1.496566, 0.09691907);
  if(hi==6 && mu>10 && mu<=11) initial.beta=c(-1.53867, 0.07202364);
  if(hi==6 && mu>13 && mu<=14) initial.beta=c(-1.635951, 0.06002705);
  if(hi==6 && mu>16 && mu<=17) initial.beta=c(-1.790113, 0.05403395);
  if(hi==6 && mu==20) initial.beta=c(-2.007826, 0.05148602);
  if(hi==6 && mu>22 && mu<=25) initial.beta=c(-2.507002, 0.05245725);
  if(hi==6 && mu>25 && mu<=27) initial.beta=c(-2.931723, 0.05565776);
  if(hi==6 && mu>27 && mu<=30) initial.beta=c(-3.701797, 0.06324358);
  if(hi==6 && mu>30 && mu<=33) initial.beta=c(-5.037949, 0.07823327);
  if(hi==6 && mu>33 && mu<=36) initial.beta=c(-7.764894, 0.1104342);
  if(hi==6 && mu==39) initial.beta=c(-19.91873, 0.2521018);
  
################ hi==12 h=1/(log(X+1)+1) #############################

  if(hi==12 && mu==1) initial.beta=c(-0.8177765, 0.3342095);
  if(hi==12 && mu>2 && mu<4.5) initial.beta=c(-1.275816, 0.1478506);  
  if(hi==12 && mu>7 && mu<8) initial.beta=c(-1.332427, 0.09165817);
  if(hi==12 && mu>10 && mu<=11) initial.beta=c(-1.408890, 0.06892841);
  if(hi==12 && mu>13 && mu<=14) initial.beta=c(-1.532199, 0.0581129);
  if(hi==12 && mu>16 && mu<=17) initial.beta=c(-1.704156, 0.05274838);
  if(hi==12 && mu==20) initial.beta=c(-1.979779, 0.05044074);
  if(hi==12 && mu>22 && mu<=25) initial.beta=c(-2.375869, 0.05142695);
  if(hi==12 && mu>25 && mu<=27) initial.beta=c(-2.876633, 0.055151);
  if(hi==12 && mu>27 && mu<=30) initial.beta=c(-3.651990, 0.06283253);
  if(hi==12 && mu>30 && mu<=33) initial.beta=c(-4.992424, 0.07789207);
  if(hi==12 && mu>33 && mu<=36) initial.beta=c(-7.722824, 0.1101443);
  if(hi==12 && mu==39) initial.beta=c(-19.88017, 0.2518597);

#####################################################################
################ hi==8 h=1/(M-X+1) #############################

  if(hi==8 && mu==1) initial.beta=c(-0.2494144, 0.2521018);
  if(hi==8 && mu>2 && mu<4.5) initial.beta=c(-0.9930963, 0.1294459);  
  if(hi==8 && mu>7 && mu<8) initial.beta=c(-1.175864, 0.08619166);
  if(hi==8 && mu>10 && mu<=11) initial.beta=c(-1.310169, 0.06713857);
  if(hi==8 && mu>13 && mu<=14) initial.beta=c(-1.491110, 0.05658972);
  if(hi==8 && mu>16 && mu<=17) initial.beta=c(-1.689578, 0.05245725);
  if(hi==8 && mu==20) initial.beta=c(-1.948965, 0.05119379);
  if(hi==8 && mu>22 && mu<=25) initial.beta=c(-2.374481, 0.05286741);
  if(hi==8 && mu>25 && mu<=27) initial.beta=c(-2.922270, 0.05754009);
  if(hi==8 && mu>27 && mu<=30) initial.beta=c(-3.801144, 0.06707397);
  if(hi==8 && mu>30 && mu<=33) initial.beta=c(-5.803752, 0.09130015);
  if(hi==8 && mu>33 && mu<=36) initial.beta=c(-9.871037, 0.1419051);
  if(hi==8 && mu==39) initial.beta=c(-26.04635, 0.3380987);

#####################################################################


#PMFs 
 if(hi==1) PMF1=repar.semi.altham(hi=1, M, mu, sigma2, initial.beta=initial.beta)$Estimates;
 if(hi==2) PMF2=repar.semi.altham(hi=2, M, mu, sigma2, initial.beta=initial.beta)$Estimates;
 if(hi==3) PMF3=repar.semi.altham(hi=3, M, mu, sigma2, initial.beta=initial.beta)$Estimates;
 if(hi==4) PMF4=repar.semi.altham(hi=4, M, mu, sigma2, initial.beta=initial.beta)$Estimates;
 if(hi==5) PMF5=repar.semi.altham(hi=5, M, mu, sigma2, initial.beta=initial.beta)$Estimates;
 if(hi==6) PMF6=repar.semi.altham(hi=6, M, mu, sigma2, initial.beta=initial.beta)$Estimates;
 if(hi==7) PMF7=repar.semi.altham(hi=7, M, mu, sigma2, initial.beta=initial.beta)$Estimates;
 if(hi==8) PMF8=repar.semi.altham(hi=8, M, mu, sigma2, initial.beta=initial.beta)$Estimates;
 if(hi==9) PMF9=repar.semi.altham(hi=9, M, mu, sigma2, initial.beta=initial.beta)$Estimates;
 if(hi==10) PMF10=repar.semi.altham(hi=10, M, mu, sigma2, initial.beta=initial.beta)$Estimates;
 if(hi==11) PMF11=repar.semi.altham(hi=11, M, mu, sigma2, initial.beta=initial.beta)$Estimates;
 if(hi==12) PMF12=repar.semi.altham(hi=12, M, mu, sigma2, initial.beta=initial.beta)$Estimates;

# vs binomial

  PMF=dbinom(0:M, M, mu/M)

# lines to plot

 lines(0:M, PMF, type="l", col="red", lty=2)
 if(hi==1) lines(0:M, PMF1, type="l", col="blue",lty=1, ylim=c(0, max(c(PMF1, PMF))))
 if(hi==2) lines(0:M, PMF2, type="l", col="blue", lty=1, ylim=c(0, max(c(PMF2, PMF))))
 if(hi==3) lines(0:M, PMF3, type="l", col="blue", lty=1, ylim=c(0, max(c(PMF3, PMF))))
 if(hi==4) lines(0:M, PMF4, type="l", col="blue", lty=1, ylim=c(0, max(c(PMF4, PMF))))
 if(hi==5) lines(0:M, PMF5, type="l", col="blue", lty=1, ylim=c(0, max(c(PMF5, PMF))))
 if(hi==6) lines(0:M, PMF6, type="l", col="blue", lty=1, ylim=c(0, max(c(PMF6, PMF))))
 if(hi==7) lines(0:M, PMF7, type="l", col="blue", lty=1, ylim=c(0, max(c(PMF7, PMF))))
 if(hi==8) lines(0:M, PMF8, type="l", col="blue", lty=1, ylim=c(0, max(c(PMF8, PMF))))
 if(hi==9) lines(0:M, PMF9, type="l", col="blue", lty=1, ylim=c(0, max(c(PMF9, PMF))))
 if(hi==10) lines(0:M, PMF10, type="l", col="blue", lty=1, ylim=c(0, max(c(PMF10, PMF))))
 if(hi==11) lines(0:M, PMF11, type="l", col="blue", lty=1, ylim=c(0, max(c(PMF11, PMF))))
 if(hi==12) lines(0:M, PMF12, type="l", col="blue", lty=1, ylim=c(0, max(c(PMF12, PMF))))


 }
}

dev.off()

############################################################
## 2B. comparison-semi.altham-binomial hi=7-12

pdf("comparison-semi.altham-binomial-7-12.pdf")

M=40
Support=0:M

#par(mfrow=c(1,1))
par(mfrow=c(2,3))

MEANS=seq(1, M-1, length.out=as.integer(M/3))

#plot(Support, seq(0, 0.4, length.out=M+1), 
#     xlab="X", ylab="Probability", type="n")
# abline(v=Support, lty=3, col=8)

for(hi in 7:12)
{

if (hi==1) h="rep(1,M+1)";
if (hi==2)  h="log(M-X+1)+1";
if (hi==3) h="choose(M, X)";
if (hi==4) h="1/choose(M, X)"; 
if (hi==5) h="X+1";
if (hi==6) h="1/(X+1)"; 
if (hi==7) h="M-X+1";
if (hi==8) h="1/(M-X+1)";
if (hi==9) h="X*(M-X)+1"; 
if (hi==10)  h="1/(X*(M-X)+1)";
if (hi==11)  h="log(X+1)+1";
if (hi==12) h="1/(log(X+1)+1)"; 


plot(Support, seq(0, 0.4, length.out=M+1), 
     xlab="X", ylab="Probability", type="n", main=h)
 abline(v=Support, lty=3, col=8)

 for(mu in MEANS)
 {
  sigma2=mu*(1-mu/M)
#  print(paste("mu =", mu, " sigma2 =", sigma2))

# initial values

#################### hi=1 h=1 ######################################

  if(hi==1 && mu==1) initial.beta=c(0,0);
  if(hi==1 && mu>1 && mu<4.5) initial.beta=c(-1.017185, 0.1290705);  
  if(hi==1 && mu>7 && mu<8) initial.beta=c(-1.199233, 0.08574957);
  if(hi==1 && mu>10 && mu<=11) initial.beta=c(-1.287356, 0.07153209);
  if(hi==1 && mu>13 && mu<=14) initial.beta=c(-1.428352, 0.05951514);
  if(hi==1 && mu>16 && mu<=17) initial.beta=c(-1.599992, 0.05333307);
  if(hi==1 && mu==20) initial.beta=c(-1.904762, 0.05012531);
  if(hi==1 && mu>22 && mu<=25) initial.beta=c(-2.105263, 0.05012531);
  if(hi==1 && mu>25 && mu<=27) initial.beta=c(-2.666654, 0.05333307);
  if(hi==1 && mu>27 && mu<=30) initial.beta=c(-3.332859, 0.05951514);
  if(hi==1 && mu>30 && mu<=33) initial.beta=c(-4.977422, 0.07776383);
  if(hi==1 && mu>33 && mu<=36) initial.beta=c(-7.70935, 0.1100388);
  if(hi==1 && mu==39) initial.beta=c(-15.02605, 0.196416);

################ hi=3 h=choose(M, X); #################################

  if(hi==3 && mu==1) initial.beta=c(3.637586, 2.079528e-16);
  if(hi==3 && mu>1 && mu<4.5) initial.beta=c(2.197225, 2.572754e-16);  
  if(hi==3 && mu>7 && mu<8) initial.beta=c(1.734601, -1.597501e-15);
  if(hi==3 && mu>10 && mu<=11) initial.beta=c(1.098612, -7.651001e-16);
  if(hi==3 && mu>13 && mu<=14) initial.beta=c(0.8472979, -1.672698e-16);
  if(hi==3 && mu>16 && mu<=17) initial.beta=c(0.4054651, 9.438626e-18);
  if(hi==3 && mu==20) initial.beta=c(0.1000835, 1.377368e-16);
  if(hi==3 && mu>22 && mu<=25) initial.beta=c(-0.3022809, -1.860093e-15);
  if(hi==3 && mu>25 && mu<=27) initial.beta=c(-0.6190392, 2.340299e-16);
  if(hi==3 && mu>27 && mu<=30) initial.beta=c(-0.9694006, -9.076763e-15);
  if(hi==3 && mu>30 && mu<=33) initial.beta=c(-1.386294, -3.049898e-14);
  if(hi==3 && mu>33 && mu<=36) initial.beta=c(-1.94591, -9.775039e-15);
  if(hi==3 && mu==39) initial.beta=c(-3.245193, 3.818461e-13);

#################### hi=5  h=X+1; ########################################

  if(hi==5 && mu==1) initial.beta=c(0.4487824, 0.167705);
  if(hi==5 && mu>4 && mu<4.5) initial.beta=c(-0.5765087, 0.1034802);
  if(hi==5 && mu>7.2 && mu<=7.5) initial.beta=c(-0.8717244, 0.07877348);
  if(hi==5 && mu>10 && mu<=11) initial.beta=c(-1.131407, 0.06147013);
  if(hi==5 && mu>13 && mu<=14) initial.beta=c(-1.362723, 0.05315132);
  if(hi==5 && mu>16 && mu<=17) initial.beta=c(-1.583141, 0.04976643);
  if(hi==5 && mu==20) initial.beta=c(-1.851086, 0.04874663);
  if(hi==5 && mu>22 && mu<=25) initial.beta=c(-2.268400, 0.05023506);
  if(hi==5 && mu>25 && mu<=27) initial.beta=c(-2.782498, 0.05423105);
  if(hi==5 && mu>27 && mu<=30) initial.beta=c(-3.568240, 0.06210076);
  if(hi==5 && mu>30 && mu<=33) initial.beta=c(-5.238844, 0.08104664);
  if(hi==5 && mu>33 && mu<=36) initial.beta=c(-8.387756, 0.1183825);
  if(hi==5 && mu==39) initial.beta=c(-19.81765, 0.2514471);

################ hi=7 h=M-X+1; ##########################################

  if(hi==7 && mu==1) initial.beta=c(-0.2981152, 0.2514471);
  if(hi==7 && mu>1 && mu<4.5) initial.beta=c(-1.041274, 0.1286951);  
  if(hi==7 && mu>7 && mu<8) initial.beta=c(-1.222602, 0.08530752);
  if(hi==7 && mu>10 && mu<=11) initial.beta=c(-1.353852, 0.06607121);
  if(hi==7 && mu>13 && mu<=14) initial.beta=c(-1.527676, 0.05521734);
  if(hi==7 && mu>16 && mu<=17) initial.beta=c(-1.714714, 0.05070316);
  if(hi==7 && mu==20) initial.beta=c(-1.953599, 0.04887204);
  if(hi==7 && mu>22 && mu<=25) initial.beta=c(-2.331937, 0.04944631);
  if(hi==7 && mu>25 && mu<=27) initial.beta=c(-2.794208, 0.05239318);
  if(hi==7 && mu>27 && mu<=30) initial.beta=c(-3.479322, 0.05845432);
  if(hi==7 && mu>30 && mu<=33) initial.beta=c(-4.838769, 0.0724201);
  if(hi==7 && mu>33 && mu<=36) initial.beta=c(-7.121272, 0.09716281);
  if(hi==7 && mu==39)initial.beta=c(-13.86518, 0.167705);

################ hi=4 h=1/choose(M, X); ##########################################


  if(hi==4 && mu==1) initial.beta=c(-4.303162, 0.6112999);
  if(hi==4 && mu>1 && mu<4.5) initial.beta=c(-4.314866, 0.2677898);  
  if(hi==4 && mu>7 && mu<8) initial.beta=c(-4.057696, 0.1787719);
  if(hi==4 && mu>10 && mu<=11) initial.beta=c(-3.850276, 0.1373737);
  if(hi==4 && mu>13 && mu<=14) initial.beta=c(-3.740427, 0.1156765);
  if(hi==4 && mu>16 && mu<=17) initial.beta=c(-3.754178, 0.1046272);
  if(hi==4 && mu==20) initial.beta=c(-3.954904, 0.1001240);
  if(hi==4 && mu>22 && mu<=25) initial.beta=c(-4.415353, 0.1025662);
  if(hi==4 && mu>25 && mu<=27) initial.beta=c(-5.154624, 0.1110587);
  if(hi==4 && mu>27 && mu<=30) initial.beta=c(-6.484727, 0.1285782);
  if(hi==4 && mu>30 && mu<=33) initial.beta=c(-8.960896, 0.1617539);
  if(hi==4 && mu>33 && mu<=36) initial.beta=c(-14.10145, 0.2291585);
  if(hi==4 && mu==39) initial.beta=c(-32.33174, 0.4576186);


################ hi=9 h=X*(M-X)+1; ##########################################


  if(hi==9 && mu==1) initial.beta=c(2.269248, -0.05154036);
  if(hi==9 && mu>1 && mu<4.5) initial.beta=c(-0.04780887, 0.06719773);  
  if(hi==9 && mu>7 && mu<8) initial.beta=c(-0.8405965, 0.06993819);
  if(hi==9 && mu>10 && mu<=11) initial.beta=c(-1.113251, 0.0594664);
  if(hi==9 && mu>13 && mu<=14) initial.beta=c(-1.326611, 0.05271925);
  if(hi==9 && mu>16 && mu<=17) initial.beta=c(-1.543899, 0.04895616);
  if(hi==9 && mu==20) initial.beta=c(-1.795508, 0.04739314);
  if(hi==9 && mu>22 && mu<=25) initial.beta=c(-2.2346, 0.04821576);
  if(hi==9 && mu>25 && mu<=27) initial.beta=c(-2.697669, 0.05118175);
  if(hi==9 && mu>27 && mu<=30) initial.beta=c(-3.359099, 0.05682025);
  if(hi==9 && mu>30 && mu<=33) initial.beta=c(-4.340993, 0.0660842);
  if(hi==9 && mu>33 && mu<=36) initial.beta=c(-5.567611, 0.07615741);
  if(hi==9 && mu==39) initial.beta=c(1.85398, -0.05154036);

######################################

################ hi=10 h=1/(X*(M-X)+1); ##########################################


  if(hi==10 && mu==1) next; #initial.beta=c(-4.035039, 0.9901828);
  if(hi==10 && mu>1 && mu<4.5) initial.beta=c(-2.510631, 0.2594033);  
  if(hi==10 && mu>7 && mu<8) initial.beta=c(-1.916851, 0.1245232);
  if(hi==10 && mu>10 && mu<=11) initial.beta=c(-1.680665, 0.07978527);
  if(hi==10 && mu>13 && mu<=14) initial.beta=c(-1.659662, 0.06209267);
  if(hi==10 && mu>16 && mu<=17) initial.beta=c(-1.791217, 0.05526399);
  if(hi==10 && mu==20) initial.beta=c(-2.059873, 0.05274772);
  if(hi==10 && mu>22 && mu<=25) initial.beta=c(-2.472172, 0.05410599);
  if(hi==10 && mu>25 && mu<=27) initial.beta=c(-3.035043, 0.05907201);
  if(hi==10 && mu>27 && mu<=30) initial.beta=c(-4.106849, 0.07198864);
  if(hi==10 && mu>30 && mu<=33) initial.beta=c(-6.565777, 0.1046964);
  if(hi==10 && mu>33 && mu<=36) initial.beta=c(-13.29024, 0.1942729);
  if(hi==10 && mu==39)
  {
  mu=37.7
  sigma2=mu*(1-mu/M)
  initial.beta=c(-35.10269, 0.478259);
  }
######################################

################ hi=11 h=log(X+1)+1; ##########################################


  if(hi==11 && mu==1) initial.beta=c(0.2635505, 0.1722706);
  if(hi==11 && mu>2 && mu<4.5) initial.beta=c(-0.7727501, 0.1118146);  
  if(hi==11 && mu>7 && mu<8) initial.beta=c(-1.074471, 0.08038401);
  if(hi==11 && mu>10 && mu<=11) initial.beta=c(-1.257837, 0.06440741);
  if(hi==11 && mu>13 && mu<=14) initial.beta=c(-1.461840, 0.0548908);
  if(hi==11 && mu>16 && mu<=17) initial.beta=c(-1.666313, 0.0509663);
  if(hi==11 && mu==20) initial.beta=c(-1.922688, 0.04962248);
  if(hi==11 && mu>22 && mu<=25) initial.beta=c(-2.262122, 0.05050303);
  if(hi==11 && mu>25 && mu<=27) initial.beta=c(-2.837562, 0.05473732);
  if(hi==11 && mu>27 && mu<=30) initial.beta=c(-3.618037, 0.06251163);
  if(hi==11 && mu>30 && mu<=33) initial.beta=c(-5.283736, 0.0813781);
  if(hi==11 && mu>33 && mu<=36) initial.beta=c(-8.42932, 0.1186652);
  if(hi==11 && mu==39) initial.beta=c(-19.85621, 0.2516892);
 

################# hi=2 h=log(M-X+1)##################################

  if(hi==2 && mu==1) initial.beta=c(-0.2803165, 0.2516607);
  if(hi==2 && mu>2 && mu<4.5) initial.beta=c(-1.023629, 0.128937);  
  if(hi==2 && mu>7 && mu<8) initial.beta=c(-1.205370, 0.08558768);
  if(hi==2 && mu>10 && mu<=11) initial.beta=c(-1.33747, 0.06640245);
  if(hi==2 && mu>13 && mu<=14) initial.beta=c(-1.485527, 0.05671613);
  if(hi==2 && mu>16 && mu<=17) initial.beta=c(-1.66826, 0.05173486);
  if(hi==2 && mu==20) initial.beta=c(-1.946819, 0.04951499);
  if(hi==2 && mu>22 && mu<=25) initial.beta=c( -2.334893, 0.05031695);
  if(hi==2 && mu>25 && mu<=27) initial.beta=c(-2.811033, 0.05354517);
  if(hi==2 && mu>27 && mu<=30) initial.beta=c(-3.513287, 0.05991752);
  if(hi==2 && mu>30 && mu<=33) initial.beta=c(-4.595194, 0.07080732);
  if(hi==2 && mu>33 && mu<=36) initial.beta=c( -6.038489, 0.08411448);
  if(hi==2 && mu==39)initial.beta=c(-11.68763, 0.1494534);


####################hi=6 h=1/(X+1) ##############################

  if(hi==6 && mu==1) initial.beta=c(-1.001554, 0.3380987);
  if(hi==6 && mu>2 && mu<4.5) initial.beta=c(-1.472866, 0.1563312);  
  if(hi==6 && mu>7 && mu<8) initial.beta=c(-1.496566, 0.09691907);
  if(hi==6 && mu>10 && mu<=11) initial.beta=c(-1.53867, 0.07202364);
  if(hi==6 && mu>13 && mu<=14) initial.beta=c(-1.635951, 0.06002705);
  if(hi==6 && mu>16 && mu<=17) initial.beta=c(-1.790113, 0.05403395);
  if(hi==6 && mu==20) initial.beta=c(-2.007826, 0.05148602);
  if(hi==6 && mu>22 && mu<=25) initial.beta=c(-2.507002, 0.05245725);
  if(hi==6 && mu>25 && mu<=27) initial.beta=c(-2.931723, 0.05565776);
  if(hi==6 && mu>27 && mu<=30) initial.beta=c(-3.701797, 0.06324358);
  if(hi==6 && mu>30 && mu<=33) initial.beta=c(-5.037949, 0.07823327);
  if(hi==6 && mu>33 && mu<=36) initial.beta=c(-7.764894, 0.1104342);
  if(hi==6 && mu==39) initial.beta=c(-19.91873, 0.2521018);
  
################ hi==12 h=1/(log(X+1)+1) #############################

  if(hi==12 && mu==1) initial.beta=c(-0.8177765, 0.3342095);
  if(hi==12 && mu>2 && mu<4.5) initial.beta=c(-1.275816, 0.1478506);  
  if(hi==12 && mu>7 && mu<8) initial.beta=c(-1.332427, 0.09165817);
  if(hi==12 && mu>10 && mu<=11) initial.beta=c(-1.408890, 0.06892841);
  if(hi==12 && mu>13 && mu<=14) initial.beta=c(-1.532199, 0.0581129);
  if(hi==12 && mu>16 && mu<=17) initial.beta=c(-1.704156, 0.05274838);
  if(hi==12 && mu==20) initial.beta=c(-1.979779, 0.05044074);
  if(hi==12 && mu>22 && mu<=25) initial.beta=c(-2.375869, 0.05142695);
  if(hi==12 && mu>25 && mu<=27) initial.beta=c(-2.876633, 0.055151);
  if(hi==12 && mu>27 && mu<=30) initial.beta=c(-3.651990, 0.06283253);
  if(hi==12 && mu>30 && mu<=33) initial.beta=c(-4.992424, 0.07789207);
  if(hi==12 && mu>33 && mu<=36) initial.beta=c(-7.722824, 0.1101443);
  if(hi==12 && mu==39) initial.beta=c(-19.88017, 0.2518597);

#####################################################################
################ hi==8 h=1/(M-X+1) #############################

  if(hi==8 && mu==1) initial.beta=c(-0.2494144, 0.2521018);
  if(hi==8 && mu>2 && mu<4.5) initial.beta=c(-0.9930963, 0.1294459);  
  if(hi==8 && mu>7 && mu<8) initial.beta=c(-1.175864, 0.08619166);
  if(hi==8 && mu>10 && mu<=11) initial.beta=c(-1.310169, 0.06713857);
  if(hi==8 && mu>13 && mu<=14) initial.beta=c(-1.491110, 0.05658972);
  if(hi==8 && mu>16 && mu<=17) initial.beta=c(-1.689578, 0.05245725);
  if(hi==8 && mu==20) initial.beta=c(-1.948965, 0.05119379);
  if(hi==8 && mu>22 && mu<=25) initial.beta=c(-2.374481, 0.05286741);
  if(hi==8 && mu>25 && mu<=27) initial.beta=c(-2.922270, 0.05754009);
  if(hi==8 && mu>27 && mu<=30) initial.beta=c(-3.801144, 0.06707397);
  if(hi==8 && mu>30 && mu<=33) initial.beta=c(-5.803752, 0.09130015);
  if(hi==8 && mu>33 && mu<=36) initial.beta=c(-9.871037, 0.1419051);
  if(hi==8 && mu==39) initial.beta=c(-26.04635, 0.3380987);

#####################################################################

#PMFs 
 if(hi==1) PMF1=repar.semi.altham(hi=1, M, mu, sigma2, initial.beta=initial.beta)$Estimates;
 if(hi==2) PMF2=repar.semi.altham(hi=2, M, mu, sigma2, initial.beta=initial.beta)$Estimates;
 if(hi==3) PMF3=repar.semi.altham(hi=3, M, mu, sigma2, initial.beta=initial.beta)$Estimates;
 if(hi==4) PMF4=repar.semi.altham(hi=4, M, mu, sigma2, initial.beta=initial.beta)$Estimates;
 if(hi==5) PMF5=repar.semi.altham(hi=5, M, mu, sigma2, initial.beta=initial.beta)$Estimates;
 if(hi==6) PMF6=repar.semi.altham(hi=6, M, mu, sigma2, initial.beta=initial.beta)$Estimates;
 if(hi==7) PMF7=repar.semi.altham(hi=7, M, mu, sigma2, initial.beta=initial.beta)$Estimates;
 if(hi==8) PMF8=repar.semi.altham(hi=8, M, mu, sigma2, initial.beta=initial.beta)$Estimates;
 if(hi==9) PMF9=repar.semi.altham(hi=9, M, mu, sigma2, initial.beta=initial.beta)$Estimates;
 if(hi==10) PMF10=repar.semi.altham(hi=10, M, mu, sigma2, initial.beta=initial.beta)$Estimates;
 if(hi==11) PMF11=repar.semi.altham(hi=11, M, mu, sigma2, initial.beta=initial.beta)$Estimates;
 if(hi==12) PMF12=repar.semi.altham(hi=12, M, mu, sigma2, initial.beta=initial.beta)$Estimates;


# vs binomial

  PMF=dbinom(0:M, M, mu/M)

# lines to plot

 lines(0:M, PMF, type="l", col="red", lty=2)
 if(hi==1) 
 {lines(0:M, PMF1, type="l", col="blue",lty=1, ylim=c(0, max(c(PMF1, PMF))))
 max(c(PMF1, PMF))
 }
 if(hi==2) 
 {
 lines(0:M, PMF2, type="l", col="blue", lty=1, ylim=c(0, max(c(PMF2, PMF))))
 max(c(PMF2, PMF))
 }
 if(hi==3)
 { 
 lines(0:M, PMF3, type="l", col="blue", lty=1, ylim=c(0, max(c(PMF3, PMF))))
 max(c(PMF3, PMF))
 }
 if(hi==4) 
 {
 lines(0:M, PMF4, type="l", col="blue", lty=1, ylim=c(0, max(c(PMF4, PMF))))
  max(c(PMF4, PMF))
 }
 if(hi==5)
 {
 lines(0:M, PMF5, type="l", col="blue", lty=1, ylim=c(0, max(c(PMF5, PMF))))
 max(c(PMF5, PMF))
}
 if(hi==6) 
{
 lines(0:M, PMF6, type="l", col="blue", lty=1, ylim=c(0, max(c(PMF6, PMF))))
 max(c(PMF6, PMF))
}
 if(hi==7)
{ 
lines(0:M, PMF7, type="l", col="blue", lty=1, ylim=c(0, max(c(PMF7, PMF))))
 max(c(PMF7, PMF))
}
 if(hi==8) 
{
lines(0:M, PMF8, type="l", col="blue", lty=1, ylim=c(0, max(c(PMF8, PMF))))
 max(c(PMF8, PMF))
}
 if(hi==9) 
{
lines(0:M, PMF9, type="l", col="blue", lty=1, ylim=c(0, max(c(PMF9, PMF))))
 max(c(PMF9, PMF))
}
 if(hi==10) 
{
lines(0:M, PMF10, type="l", col="blue", lty=1, ylim=c(0, max(c(PMF10, PMF))))
 max(c(PMF10, PMF))
}
 if(hi==11) 
{
lines(0:M, PMF11, type="l", col="blue", lty=1, ylim=c(0, max(c(PMF11, PMF))))
 max(c(PMF11, PMF))
}
 if(hi==12)
{
 lines(0:M, PMF12, type="l", col="blue", lty=1, ylim=c(0, max(c(PMF12, PMF))))
 max(c(PMF12, PMF))
}

#create a matrix of PMF and use its row 

# for(i in 1:12)
# {
#   tmp=paste("PMF", i, sep="")
#   lines(0:M, PMF, type="l", col="red", lty=2)
#   lines(0:M, tmp, type="l", col="blue", lty=1, ylim=c(0, max(c(tmp, PMF))))
#   max(c(tmp, PMF)) 
# }

 }
}

dev.off()