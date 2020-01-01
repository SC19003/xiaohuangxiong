## -----------------------------------------------------------------------------
queue1<-function(lambda, mu, T){
  k<-0; wt<-0; wn<-0; ws<-0;
  tp<-0; nA<-0; n<-0; t<-0
  r<-runif(1); tA<--1/lambda*log(r); tD<-Inf
  repeat{
    k<-k+1; wt[k]<-t; wn[k]<-n
    if (tA < T){
      ws[k]<-min(tA, tD)-t
      if (tA < tD){
        t<-tA; n<-n+1; nA<-nA+1
        r<-runif(1); tA<-t-1/lambda*log(r)
        if (n==1){
          r<-runif(1); tD<-t-1/mu*log(r)
        }
      }else{
        t<-tD; n<-n-1
        if (n==0){
          tD<-Inf
        }else{
          r<-runif(1); tD<-t-1/mu*log(r)
        } }
    }else{
      ws[k]<-if(tD==Inf) 0 else tD-t
      if (n>0){
        t<-tD; n<-n-1
        if (n>0){
          r<-runif(1); tD<-t-1/mu*log(r)
        }
      }else
        tp<-1
    }
    
    if (tp==1) break
  }
  data.frame(Ls=sum(ws*wn)/t, Ws=sum(ws*wn)/nA,
             Pwait=sum(ws[wn>=1])/t)
}

## ----eval=TRUE----------------------------------------------------------------
library(SC19003)
result1 <- queue1(4,10,10000) # we assume that queues are allowed forever so the parameter T is 10000 (a large number).
result1

## -----------------------------------------------------------------------------
queue2<-function(lambda, mu, T, S){
  k<-0; wt<-0; wn<-0; ws<-0
  tp<-0; nA<-0; t<-0
  r<-runif(1); tA<--1/lambda*log(r)
  tD<-rep(Inf, S); SS<-rep(0, S+1)
  repeat{
    t1<-if(SS[1]==0) Inf else min(tD)
    i1<-if(SS[1]==0) 1 else which.min(tD)
    k<-k+1; wt[k]<-t; wn[k]<-SS[1]
    if (tA < T){
      ws[k]<-min(tA, t1)-t
      if (tA < t1){
        t<-tA; nA<-nA+1
        r<-runif(1); tA<-t-1/lambda*log(r)
        n<-SS[1]; SS[1]<-n+1
        for (i in 1:S){
          if (SS[1+i]==0){
            SS[1+i]<-1
            r<-runif(1); tD[i]<-t-1/mu*log(r)
            break
          }
        }
        
      }else{
        t<-t1; n<-SS[1]; SS[1]<-n-1
        if (n==1){
          SS[2:(S+1)]<-0; tD[1:S]<-Inf
        }else if (n<=S){
          SS[1+i1]<-0; tD[i1]<-Inf
        }else{
          r<-runif(1); tD[i1]<-t-1/mu*log(r)
        } }
    }else{
      ws[k]<- if( t1==Inf) 0 else t1-t
      n<-SS[1]
      if (n>0){
        t<-t1; SS[1]<-n-1;
        if (n==1){
          SS[2:(S+1)]<-0; tD[1:S]<-Inf
        }else if (n<=S){
          SS[1+i1]<-0; tD[i1]<-Inf
        }else{
          r<-runif(1); tD[i1]<-t-1/mu*log(r)
        }
      }else
        tp<-1
    }
    if (tp==1) break
  }
  data.frame(Ls=sum(ws*wn)/t, Ws=sum(ws*wn)/nA,
             Pwait=sum(ws[wn>=S])/t)
}

## ----eval=TRUE----------------------------------------------------------------
library(SC19003)
result2 <- queue2(15,6,1000,3) # we assume that queues are allowed forever so the parameter T is 1000 (a large number).
result2

## -----------------------------------------------------------------------------
age <-c(21,30,35,44,45,33,25,27,29,48)
height <-c(160,162,165,158,155,164,162,150,167,161)
group <-gl(2,10,20,labels = c("age","height"))
weight <-c(age,height)
lm.D9 <-lm(weight~group)
summary(lm.D9)$coef

## -----------------------------------------------------------------------------
knitr::kable(head(USJudgeRatings))

## -----------------------------------------------------------------------------
par(mar=c(1,1,1,1))
plot(lm.D9)

## -----------------------------------------------------------------------------
n<-1000
u<-runif(n)
sigma<-1
x<-(-2*(sigma^2)*log(u))^(1/2)
hist(x,probability = T,main = expression(f(x)==(x/(sigma^2))*exp(-(x^2)/(2*(sigma^2)))))
y <- seq(0,4,0.01) 
lines(y,y/(sigma^2)*exp(-(y^2)/(2*(sigma^2))))

## -----------------------------------------------------------------------------
sigma<-2
x<-(-2*(sigma^2)*log(u))^(1/2)
hist(x,probability = T,main = expression(f(x)==(x/(sigma^2))*exp(-(x^2)/(2*(sigma^2)))))
y <- seq(0,8,0.01) 
lines(y,y/(sigma^2)*exp(-(y^2)/(2*(sigma^2))))

## -----------------------------------------------------------------------------
sigma<-4
x<-(-2*(sigma^2)*log(u))^(1/2)
hist(x,probability = T,main = expression(f(x)==(x/(sigma^2))*exp(-(x^2)/(2*(sigma^2)))))
y <- seq(0,15,0.01) 
lines(y,y/(sigma^2)*exp(-(y^2)/(2*(sigma^2))))

## -----------------------------------------------------------------------------
library(MASS)
n<-1000
x1<-rnorm(n,0,1)
x2<-rnorm(n,3,1)
p <- 0.75
r<-sample(c(0,1),n,replace = TRUE,prob = c(p,1-p))
z<-(1-r)*x1+r*x2
hist(z)

## -----------------------------------------------------------------------------
p <- 0.5
r<-sample(c(0,1),n,replace = TRUE,prob = c(p,1-p))
z<-(1-r)*x1+r*x2
hist(z)

## -----------------------------------------------------------------------------
n <- 10
sigma_1 <-4 * diag(4)
Wd<-function(sigma_1,n){
  p<-nrow(sigma_1)
  l<-chol(sigma_1)
  A<-matrix(nrow = p,ncol = p)
  A[upper.tri(A)]<-0
  c<-numeric(p)
  for (i in 1:p){
    c[i]<-rchisq(1,n-i+1)
  }
  diag(A)<-c
  A[lower.tri(A)]<-rnorm((p^2-p)/2)
  s<-l%*%A%*%t(A)%*%t(l)
  return(s)
}
Wd(sigma_1,n)

## -----------------------------------------------------------------------------
set.seed(123)
m<-10000
x<-runif(m,min = 0,max = pi/3)
theta.hat<-mean(sin(x))*(pi/3)
print(c(theta.hat,cos(0)-cos(pi/3)))

## -----------------------------------------------------------------------------
set.seed(12345)
n<-10000
u<-runif(n)
y1<-exp(-u)/(1+u^2)
y2<-exp(u-1)/(1+(1-u)^2)
y0<-c(y1,y2)
Use<-mean(y0)
Use
var0<-(1/4)*var(y1)+(1/4)*var(y2)+(1/2)*cov(y1,y2)
var0

## -----------------------------------------------------------------------------
set.seed(12345)
n<-10000
u<-runif(n)
y<-exp(-u)/(1+u^2)
Ise<-mean(y)
Ise
var(y)

## -----------------------------------------------------------------------------
(var(y)-var0)/(var(y))

## -----------------------------------------------------------------------------
set.seed(0)
m<-numeric(6)
a<-c(0,0.2,0.4,0.6,0.8,1)
for (i in 1:6){
  m[i]<--log(1-a[i]*(1-exp(-1)))
}
print(m)

## stratified importance sampling estimate
M <- 10000 
k <- 5 
r <- M / k 
N <- 50 
T2 <- numeric(k)
estimates <- matrix(0, N, 2)
g<-function(x) {
  exp(-x)/(1+x^2)*(x>0)*(x<1)
}

for (i in 1:N) {
  u<-runif(m)
  x<--log(1-u*(1-exp(-1)))    # inverse transformated method
  g1<-g(x)/(exp(-x)/(1-exp(-1)))
  estimates[i,1]<-mean(g1)
  for (j in 1:k) {
    v<-runif(m/5)
    y<--log(exp(-m[j])-v/5*(1-exp(-1)))
    g2<-g(y)/(5*exp(-y)/(1-exp(-1)))
    T2[j]<-mean(g2)
  }
estimates[i,2]<-sum(T2)
}
print(estimates)

## -----------------------------------------------------------------------------
apply(estimates, 2, mean)
apply(estimates, 2, var)

## -----------------------------------------------------------------------------
print((var(estimates[,1])-var(estimates[,2]))/var(estimates[,1]))

## -----------------------------------------------------------------------------
n <- 20
set.seed(012345678)
alpha <- 0.05
CL<-replicate(1000,expr={
  x <- rchisq(n, df = 2)
  UCL<-mean(x) +qt(1-alpha/2, df = n-1)*sd(x)/sqrt(n)
  LCL<-mean(x) - qt(1-alpha/2, df = n-1)*sd(x)/sqrt(n)
  c(LCL,UCL)
})
cum<-0
a<-2
for (i in 1:1000){
  if(CL[1,i]<a && a < CL[2,i] )
  cum <- cum + 1
  }
cp<- cum/1000
print(cp)

## -----------------------------------------------------------------------------
set.seed(1)
n <- 20
alpha <- .05
UCL <- replicate(1000, expr = {
x <- rchisq(n, df = 2)
(n-1) * var(x) / qchisq(alpha, df = n-1)
} )
sum(UCL > 4)
mean(UCL >4)
## The t-interval should be more robust to departures from normality than the interval for variance.In this experiment, only 797 or 79.7% of the intervals contained the population variance.

## -----------------------------------------------------------------------------
set.seed(123)
N <- 50
sk <- function(x) {
                  xbar <- mean(x)
                  m3 <- mean((x - xbar)^3)
                  m2 <- mean((x - xbar)^2)
                  return( m3 / m2^1.5 )
}
estimates <- matrix(0, N, 1)
for(i in 1:N){
  x <- rnorm(1000)
  estimates[i]<-sk(x)
}

## -----------------------------------------------------------------------------
q <- numeric(4)
a <-c(0.025,0.05,0.95,0.975)
for(i in 1:4){
           q[i]<- quantile(estimates,a[i])
}
print(q)

## -----------------------------------------------------------------------------
var <- matrix(0, 4, 1)
g <- function(x){a[i]*(1-a[i])/(N*dnorm(x,0,1)^2)
}
for(i in 1:4){
            var[i] <- g(q[i])
     }
print(sqrt(var))

## -----------------------------------------------------------------------------
set.seed(123)
n <- 1000
x <- rnorm(n,mean = 0,sd=(6/n)^(1/2))
q2 <- numeric(4)
a <-c(0.025,0.05,0.95,0.975)
for(i in 1:4){
           q2[i]<- quantile(x,a[i])
}
print(q2)

## -----------------------------------------------------------------------------
cbind(q,q2)

## -----------------------------------------------------------------------------
library(energy)
set.seed(1)
alpha <- .05
n <- 30
m <- 2000 
exam1 <- exam2 <- numeric(m)
para1 <- para2 <- 1:10
pwr1 <- pwr2 <- numeric(10)

## Estimate the power of the skewness test of normality against symmetric Beta(α, α) distributions
for (i in 1:10) {
  for (j in 1:m) {
    x <- rbeta(n,para1[i],para1[i]) 
    exam1[j] <- as.integer(
    shapiro.test(x)$p.value <= alpha)
  }
  pwr1[i] <- mean(exam1)
}
plot(para1,pwr1,main = "curve 1",xlab = "α",ylab = "power",type = "b",col = "blue")
print(pwr1)

#Estimate the power of the skewness test of normality against symmetric t(ν) distributions
for (j in 1:10) {
  for (i in 1:m) {
    y <- rt(n,para2[j]) 
    exam2[i] <- as.integer(
    shapiro.test(y)$p.value <= alpha)
  }
  pwr2[j] <- mean(exam2)
}
plot(para2,pwr2,main = "curve 2",xlab = "ν",ylab = "power",type = "b",col = "blue")
print(pwr2)

## -----------------------------------------------------------------------------
set.seed(12345678)
n<-50
alpha <- .05 
a<-1
mu0<-1
m <- 10000
#number of replicates
p <- numeric(m) 
#storage for p-values 
for (j in 1:m) { 
  x <- rchisq(n, 1)
  ttest <- t.test(x, alternative = "two.sided", mu = mu0) 
  p[j] <- ttest$p.value
  }
p.hat <- mean(p < alpha)
print(p.hat)

for (j in 1:m) { 
  x <- runif(n, 0, 2)
  ttest <- t.test(x, alternative = "two.sided", mu = mu0) 
  p[j] <- ttest$p.value
  }
p.hat <- mean(p < alpha)
print(p.hat)

for (j in 1:m) { 
  x <- rexp(n, 1)
  ttest <- t.test(x, alternative = "two.sided", mu = mu0) 
  p[j] <- ttest$p.value
  }
p.hat <- mean(p < alpha)
print(p.hat)

## -----------------------------------------------------------------------------
library(bootstrap)
## Use a panel display to display the scatter plots for each pair of test scores.
library(ggplot2)
pairs(~mec+vec+alg+ana+sta,scor)
##  Compare the plot with the sample correlation matrix.
cor(scor)

## -----------------------------------------------------------------------------
##  Obtain bootstrap estimates of the standard errors for each of the following estimates.
#set up the bootstrap
B <- 200 
n <- nrow(scor)
rho12 <- numeric(B)
set.seed(1234)

## -----------------------------------------------------------------------------
for (b in 1:B) {
i <- sample(1:n, size = n, replace = TRUE)
mec <- scor$mec[i]
vec <- scor$vec[i]
rho12[b] <- cor(mec, vec)
}

#output
print(se.rho12 <- sd(rho12))
hist(rho12, prob = TRUE)
lines(density(rho12))

## -----------------------------------------------------------------------------
rho34 <- numeric(B)
for (b in 1:B) {
i <- sample(1:n, size = n, replace = TRUE)
alg <- scor$alg[i]
ana <- scor$ana[i]
rho34[b] <- cor(alg, ana)
}

#output
print(se.rho34 <- sd(rho34))
hist(rho34, prob = TRUE)
lines(density(rho34))

## -----------------------------------------------------------------------------
rho35 <- numeric(B)
for (b in 1:B) {
i <- sample(1:n, size = n, replace = TRUE)
alg <- scor$alg[i]
sta <- scor$sta[i]
rho35[b] <- cor(alg, sta)
}

#output
print(se.rho35 <- sd(rho35))
hist(rho35, prob = TRUE)
lines(density(rho35))

## -----------------------------------------------------------------------------
rho45 <- numeric(B)
for (b in 1:B) {
i <- sample(1:n, size = n, replace = TRUE)
ana <- scor$ana[i]
sta <- scor$sta[i]
rho45[b] <- cor(ana, sta)
}

#output
print(se.rho45 <- sd(rho45))
hist(rho45, prob = TRUE)
lines(density(rho45))

## -----------------------------------------------------------------------------
## Conclusion
library(knitr)
library(kableExtra)
knitr::kable(cbind(se.rho12,se.rho34,se.rho35,se.rho45),formate = "html",col.names = c("sd of rho12","sd of rho34","sd of rho35","sd of rho45")) %>% kable_styling(position = "center")

## -----------------------------------------------------------------------------
## basic settings
library(boot)
n<-1e2
set.seed(012345678)
mu<-0;b<-1;m<-80
sk<- function(x,i) {
  xbar <- mean(x[i])
  m3 <- mean((x[i] - xbar)^3) 
  m2 <- mean((x [i]- xbar)^2) 
  return( m3 / m2^1.5 )
}
## compute the coverage rates for normal populations
norm_cil<-norm_cir<-basic_cil<-basic_cir<-perc_cil<-perc_cir<-numeric(n)
for(i in 1:m){ 
 U<-rnorm(m,0,1) 
  de <- boot(data=U,statistic=sk, R = 2000) 
  ci <- boot.ci(de ,type=c("norm","basic","perc"))
  norm_cil[i]<-ci$normal[2]
  norm_cir[i]<-ci$normal[3]
  basic_cil[i]<-ci$basic[4]
  basic_cir[i]<-ci$basic[5]
  perc_cil[i]<-ci$percent[4]
  perc_cir[i]<-ci$percent[5]
}
cat('norm =',round(mean(norm_cil<=mu & norm_cir>=mu),4), 
    'basic =',mean(basic_cil<=mu & basic_cir>=mu), 'perc =',mean(perc_cil<=mu & perc_cir>=mu))

## -----------------------------------------------------------------------------
#  compute the coverage rates for chi-square(5) distributions
set.seed(012345678)
sk<- function(x,i) {
  xbar <- mean(x[i])
  m3 <- mean((x[i] - xbar)^3) 
  m2 <- mean((x [i]- xbar)^2) 
  return( m3 / m2^1.5 )
}
n<-500      
m<-200      
norm_cil<-norm_cir<-basic_cil<-basic_cir<-perc_cil<-perc_cir<-numeric(n)

for (i in 1:n) {
  x<-rchisq(m,5)
  boot.obj<-boot(x,statistic=sk,R=2000)
  ci<-boot.ci(boot.obj,type=c('norm','basic','perc'))
  norm_cil[i]<-ci$normal[2]
  norm_cir[i]<-ci$normal[3]
  basic_cil[i]<-ci$basic[4]
  basic_cir[i]<-ci$basic[5]
  perc_cil[i]<-ci$percent[4]
  perc_cir[i]<-ci$percent[5]
}
# compute coverage rates 
a<-sqrt(8/5) 
print(c(mean((norm_cil<a)*(norm_cir>a)),mean((basic_cil<a)*(basic_cir>a)),mean((perc_cil<a)*(perc_cir>a))))

## -----------------------------------------------------------------------------
set.seed(0)
library(bootstrap)
library(knitr)
lamda_hat <- eigen(cov(scor))$values
theta_hat <- lamda_hat[1]/sum(lamda_hat)
print(theta_hat)
n<-nrow(scor)
theta.j<-numeric(n)
ev<-matrix(nrow = 88,ncol = 5)

for (i in 1:n){
  ev[i,]<- eigen(cov(scor[-i,]))$values  
  theta.j[i]<-ev[i,1]/sum(ev[i,])}

bias.j<-(n-1)*(mean(theta.j)-theta_hat)
se.j<-sqrt((n-1)*mean((theta.j-mean(theta.j))^2))
print(bias.j)
print(se.j)

library(kableExtra)
knitr::kable(cbind(theta_hat,bias.j,se.j),formate = "html",col.names = c("theta_hat","bias.j","se.j")) %>% kable_styling(position = "center")

## -----------------------------------------------------------------------------
set.seed(123)
library(DAAG); attach(ironslag)
a <- seq(10, 40, .1) #sequence for plotting fits

par(mar=c(1,1,1,1))
L1 <- lm(magnetic ~ chemical)
plot(chemical, magnetic, main="Linear", pch=16)
yhat1 <- L1$coef[1] + L1$coef[2] * a
lines(a, yhat1, lwd=2)

L2 <- lm(magnetic ~ chemical + I(chemical^2))
plot(chemical, magnetic, main="Quadratic", pch=16)
yhat2 <- L2$coef[1] + L2$coef[2] * a + L2$coef[3] * a^2
lines(a, yhat2, lwd=2)

L3 <- lm(log(magnetic) ~ chemical)
plot(chemical, magnetic, main="Exponential", pch=16)
logyhat3 <- L3$coef[1] + L3$coef[2] * a
yhat3 <- exp(logyhat3)
lines(a, yhat3, lwd=2)

L4 <- lm(magnetic ~ chemical + I(chemical^2) + I(chemical^3))
plot(chemical, magnetic, main="Cubic polynomial", pch=16)
yhat4 <- L4$coef[1] + L4$coef[2] * a + L4$coef[3] * a^2 + L4$coef[4] * a^3
lines(a, yhat4, lwd=2)

## Cross validation
n <- length(magnetic) #in DAAG ironslag
e1 <- e2 <- e3 <- e4 <- numeric(n)
# for n-fold cross validation
# fit models on leave-one-out samples
for (k in 1:n) {
  y <- magnetic[-k]
  x <- chemical[-k]
J1 <- lm(y ~ x)
yhat1 <- J1$coef[1] + J1$coef[2] * chemical[k]
e1[k] <- magnetic[k] - yhat1

J2 <- lm(y ~ x + I(x^2))
yhat2 <- J2$coef[1] + J2$coef[2] * chemical[k] +
J2$coef[3] * chemical[k]^2
e2[k] <- magnetic[k] - yhat2

J3 <- lm(log(y) ~ x)
logyhat3 <- J3$coef[1] + J3$coef[2] * chemical[k]
yhat3 <- exp(logyhat3)
e3[k] <- magnetic[k] - yhat3

J4 <- lm(y ~ x + I(x^2) + I(x^3))
yhat4 <- J4$coef[1] + J4$coef[2] * chemical[k] + J4$coef[3] * chemical[k]^2 + J4$coef[4] * chemical[k]^3
e4[k] <- magnetic[k] - yhat4
}
c(mean(e1^2), mean(e2^2), mean(e3^2), mean(e4^2))
library(kableExtra)
knitr::kable(cbind(mean(e1^2), mean(e2^2), mean(e3^2), mean(e4^2)),formate = "html",col.names = c("mean(e1^2)", "mean(e2^2)", "mean(e3^2)", "mean(e4^2)")) %>% kable_styling(position = "center")

## the values of Adjusted R-squared
knitr::kable(cbind(summary(lm(L1))$adj.r.squared,summary(lm(L2))$adj.r.squared,summary(lm(L3))$adj.r.squared,summary(lm(L4))$adj.r.squared),formate = "html",col.names = c("L1 adj-R^2","L2 adj-R^2","L3 adj-R^2","L4 adj-R^2")) %>% kable_styling(position = "center")

## -----------------------------------------------------------------------------
set.seed(123)
n1 <- 30;n2 <- 50
mu1 <- mu2 <- 0;sigma1 <- sigma2 <- 1;m <- 10000
x <- rnorm(n1, mu1, sigma1)
y <- rnorm(n2, mu2, sigma2)

m1 <-log(.025) / log(n1 / (n1 + n2))
m2 <-log(.025) / log(n2 / (n1 + n2))
count5test <- function(x, y) {
X <- x - mean(x)
Y <- y - mean(y)
outx <- sum(X > max(Y)) + sum(X < min(Y))
outy <- sum(Y > max(X)) + sum(Y < min(X))
return(as.integer((outx > 4) | (outy > 7)))
}

R <- 10000
z <- c(x,y)
K <- n1 + n2
reps <- numeric(R)
t0 <- count5test(x,y)
for (i in 1:R) {
  k <- sample(K, size = n1, replace = FALSE)
  x1 <- z[k];y1 <- z[-k]
  X <- x1 - mean(x1); Y <- y1 - mean(y1)
  reps[i] <- count5test(x1, y1)
}

alphahat <- mean(c(t0, reps) > t0)
print(alphahat)

## -----------------------------------------------------------------------------
# functions
dCov <- function(x, y) {
  x <- as.matrix(x); y <- as.matrix(y)
  n <- nrow(x); m <- nrow(y)
  if (n != m || n < 2) stop("Sample sizes must agree")
  if (! (all(is.finite(c(x, y)))))
  stop("Data contains missing or infinite values")
  Akl <- function(x) {
  d <- as.matrix(dist(x))
  m <- rowMeans(d); M <- mean(d)
  a <- sweep(d, 1, m); b <- sweep(a, 2, m)
  b + M
  }
A<- Akl(x); B <- Akl(y)
sqrt(mean(A * B))
}

ndCov2 <- function(z, ix, dims) {
p <- dims[1]
q <- dims[2]
d <- p + q
x <- z[ , 1:p] #leave x as is
y <- z[ix, -(1:p)] #permute rows of y
return(nrow(z) * dCov(x, y)^2)
}

library(Ball)
library(mvtnorm)
library(boot)
n<-seq(from=10,to=100,by=10)
k<-100
alpha<-0.05
pow_dCor_Model1<-pow_ball_Model1<-pow_dCor_Model2<-pow_ball_Model2<-numeric(length(n))
for (j in 1:length(n)) {

  #storage of temp data
  p_ball1<-numeric(k)
  dcor1<-numeric(k)
  p_ball2<-numeric(k)
  dcor2<-numeric(k)
  dcor1<-dcor2<-p_ball1<-p_ball2<-numeric(k)
  for (i in 1:k) {
    set.seed(i)
    # the function "rmvnorm" is used to 
    # generate the multidimensional normal data
    X<-rmvnorm(n[j],rep(0,2),diag(1,2))
    err<-rmvnorm(n[j],rep(0,2),diag(1,2))
    Y1<-(X/4)+err
    Y2<-(X/4)*err
    Z1<-cbind(X,Y1)
    Z2<-cbind(X,Y2)
    t1<-bcov.test(X,Y2,R=99)
    p_ball1[i]<-t1$p.value
    boot.obj1<-boot(data=Z1,statistic=ndCov2,R=99,sim="permutation",dims=c(2, 2))
    temp1<-c(boot.obj1$t0, boot.obj1$t)
    dcor1[i]<-mean(temp1>=temp1[1])
    
    t2<-bcov.test(X,Y2,R=99)
    p_ball2[i]<-t2$p.value
    boot.obj2<-boot(data=Z2,statistic=ndCov2,R=99,sim="permutation",dims=c(2, 2))
    temp2<-c(boot.obj2$t0, boot.obj2$t)
    dcor2[i]<-mean(temp2>=temp2[1])
    }
  pow_dCor_Model1[j]<-mean(dcor1<alpha)
  pow_ball_Model1[j]<-mean(p_ball1<alpha)
  pow_dCor_Model2[j]<-mean(dcor2<alpha)
  pow_ball_Model2[j]<-mean(p_ball2<alpha)  
}
dat<-data.frame(pow_dCor_Model1,pow_ball_Model1,pow_dCor_Model2,pow_ball_Model2)
# the red one is distance correlation test and the blue one is ballcovariance test
library(ggplot2)
ggplot(dat,aes(n))+geom_point(y=pow_dCor_Model1,fill="white")+geom_line(y=pow_dCor_Model1,colour="purple")+geom_point(y=pow_ball_Model1,fill="white")+geom_line(y=pow_ball_Model1,colour="green")

ggplot(dat,aes(n))+geom_point(y=pow_dCor_Model2,fill="white")+geom_line(y=pow_dCor_Model2,colour="purple")+geom_point(y=pow_ball_Model2,fill="white")+geom_line(y=pow_ball_Model2,colour="green")

## -----------------------------------------------------------------------------
library(VGAM)
library(knitr)
library(kableExtra)
set.seed(123)
    rw.Metropolis <- function(n, sigma, x0, N) {
        x <- numeric(N)
        x[1] <- x0
        u <- runif(N)
        k <- 0
        for (i in 2:N) {
            y <- rnorm(1, x[i-1], sigma)
                if (u[i] <= (dlaplace(y, 0, 1) / dlaplace(x[i-1], 0, 1)))
                    x[i] <- y  
                else {
                    x[i] <- x[i-1]
                    k <- k + 1
                }
            }
        return(list(x=x, k=k))
    }

    N <- 2000
    sigma <- c(.05, .5, 2,  16)

    x0 <- 25
    rw1 <- rw.Metropolis(n, sigma[1], x0, N)
    rw2 <- rw.Metropolis(n, sigma[2], x0, N)
    rw3 <- rw.Metropolis(n, sigma[3], x0, N)
    rw4 <- rw.Metropolis(n, sigma[4], x0, N)

    #number of candidate points rejected
    no.reject <- data.frame(sigma=sigma,no.reject=c(rw1$k, rw2$k, rw3$k, rw4$k))
    knitr::kable(no.reject,format='html')%>% kable_styling(position = "center")
    
    #compute the acceptance rates of each chain
    ar <- c(1-rw1$k/N, 1-rw2$k/N, 1-rw3$k/N, 1-rw4$k/N)
    print(ar)
    
    par(mar=c(1,1,1,1))  #display 4 graphs together
    refline <- qlaplace(c(.025, .975))
    rw <- cbind(rw1$x, rw2$x, rw3$x,  rw4$x)
    for (j in 1:4) {
        plot(rw[,j], type="l",
             xlab=bquote(sigma == .(round(sigma[j],3))),
             ylab="X", ylim=range(rw[,j]))
        abline(h=refline)
    }
    par(mar=c(1,1,1,1)) #reset to default
    
    a <- c(.05, seq(.1, .9, .1), .95)
    Q <- qlaplace(a, 0, 1)
    rw <- cbind(rw1$x, rw2$x, rw3$x, rw4$x)
    mc <- rw[501:N, ]
    Qrw <- apply(mc, 2, function(x) quantile(x, a))
    qq <- data.frame(round(cbind(Q, Qrw), 3))
    names(qq) <- c('True','sigma=0.05','sigma=0.5','sigma=2','sigma=16')
    knitr::kable(qq,format='html')%>% kable_styling(position = "center") 

## -----------------------------------------------------------------------------
x <- log(exp(100))
y <- exp(log(100))
isTRUE(x == y)
isTRUE(all.equal(x,y))
print(x-y)

## -----------------------------------------------------------------------------
f <- function(a, k) 1-pt(sqrt(a^2*k/(k+1-a^2)), k)
k <- c(25, 50, 100, 500)
for (i in k) {
  h <- function(x) {f(x, i-1)-f(x, i)}
    a <- seq(0, sqrt(i), .1)
  plot(a, h(a), type = 'b')
  abline(h = 0)
}


## -----------------------------------------------------------------------------
f <- function(a, k) 1-pt(sqrt(a^2*k/(k+1-a^2)), k)
k <- c(25:50, 100, 250, 500)
Ak <- numeric(length(k))
i <- 1
for (j in k) {
  g <- function(x) f(x, j-1)-f(x, j)
  Ak[i] <- uniroot(g, lower = 1, upper = 2)$root
  i <- i+1
}
knitr::kable(cbind(k,Ak))

## -----------------------------------------------------------------------------
f <- function(k) 2/sqrt(pi*k)*exp(lgamma((k+1)/2)-lgamma(k/2))
ck <- function(a, k) sqrt(a^2*k/(k+1-a^2))
g <- function(u, k) (1+u^2/k)^(-(k+1)/2)
k <- c(25:50, 100, 250, 500)
root <- numeric(length(k))
i <- 1
for (j in k) {
  ff <- function(a) f(j)*integrate(function(u) {g(u, j)}, 0, ck(a, j))$value -       f(j-1)*integrate(function(u) {g(u, j-1)}, 0, ck(a, j-1))$value 
  root[i] <- uniroot(ff, lower = 1, upper = 2)$root
  i <- i+1
}
knitr::kable(cbind(k, Ak, root))

## ----echo=FALSE---------------------------------------------------------------
    dat <- rbind(Genotype=c('AA','BB','OO','AO','BO','AB','Sum'),
                 Frequency=c('p^2','q^2','r^2','2pr','2qr','2pq',1),
                 Count=c('nAA','nBB','nOO','nAO','nBO','nAB','n'))
    knitr::kable(dat,format='html')

## -----------------------------------------------------------------------------
nA <- 28; nB <- 24; nOO <- 41; nAB <- 70

theta_0 <- c(.1, .1)
for (j in 1:1000) {
  E <- function(theta) {
    p <- theta[1]
    q <- theta[2]
    r <- 1-p-q
    p0 <- theta_0[1]
    q0 <- theta_0[2]
    r0 <- 1-p0-q0
    return(2*nA*(log(p)+r0/(p0+2*r0)*log(2*r/p))+2*nB*(log(q)+r0/(q0+2*r0)*log(2*r/q))+2*nOO*log(r)+nAB*log(2*p*q))
  }
  Obj <- function(theta) -E(theta)
  optim <- optim(c(.1, .1), Obj)
  theta_0 <- optim$par
}
print(theta_0)

## ----warning=FALSE------------------------------------------------------------
data <- as.matrix(mtcars)
mpg <- data[,1]
disp <- data[,3]
wt <- data[,6]

# for loops
re <- list()
re1 <- mpg ~ disp
re2 <- mpg ~ I(1 / disp)
re3 <- mpg ~ disp + wt
re4 <- mpg ~ I(1 / disp) + wt
h <- c("re1","re2","re3","re4")

for (i in h){
  re[i] <- lm(i)
}
re

# lapply
f <- list(mpg ~ disp,
       mpg ~ I(1 / disp),
       mpg ~ disp + wt,
       mpg ~ I(1 / disp) + wt)
lapply(f, lm)


## ----warning=FALSE------------------------------------------------------------
set.seed(123)
bootstraps <- lapply(1:10, function(i) {
      rows <- sample(1:nrow(mtcars), rep = TRUE)
      mtcars[rows, ]
      })
# for loops
g <- list()
for (i in 1:10){
  g[i] <- lm(mpg~disp, data = bootstraps[[i]])
}
g
# lapply
#lapply(1:10, FUN = function(x) lm(mpg~disp, data = bootstraps[[x]]) )
lapply(bootstraps, lm, formula = mpg~disp)

## -----------------------------------------------------------------------------
rsq <- function(mod) summary(mod)$r.squared

# exercise 3
m <- lapply(f, lm)
lapply(m, rsq)

# exercise 4
n <- lapply(bootstraps, lm, formula = mpg~disp)
lapply(n, rsq)

## -----------------------------------------------------------------------------
set.seed(1234)
trials <- replicate(
        100,
        t.test(rpois(10, 10), rpois(7, 10)),
        simplify = FALSE
       )

# sapply
l <- function(x) x$p.value
sapply(trials, l)

# Extra challenge
sapply(trials, '[[',3)

## ----eval=FALSE---------------------------------------------------------------
#  # mcsapply()
#  library(parallel)
#  boot_df <- function(x) x[sample(nrow(x), rep = T), ]
#  rsquared <- function(mod) summary(mod)$r.squared
#  boot_lm <- function(i) {
#    rsquared(lm(mpg ~ wt + disp, data = boot_df(mtcars)))
#  }
#  system.time(sapply(1:1e5, boot_lm))
#  system.time(unlist(mclapply(1:1e5, boot_lm, mc.cores = 4)))

## ----warning=FALSE------------------------------------------------------------
#the function written before
library(kableExtra)
library(GeneralizedHyperbolic)
rw.Metropolis <- function(sigma,x0,N) {
x <- numeric(N)
x[1] <- x0
u <- runif(N)
k <- 0
for (i in 2:N) {
y <- rnorm(1, x[i-1], sigma)
if (u[i] <= exp(abs(x[i-1]) - abs(y))){
x[i] <- y 
k <- k + 1}
else {
x[i] <- x[i-1]

} }
return(list(x=x, k=k))
}

#the function using Rcpp
library(Rcpp)
cppFunction('List rw_Metropolis(double sigma, double x0, int N) {
NumericVector x(N);
x[0] = x0;
int k = 0;
for (int i = 1;i < N;i++) {
double u = runif(1)[0];
double y = rnorm(1, x[i-1], sigma)[0];
if (u <= exp(abs(x[i-1]) - abs(y))){
x[i] = y;
k = k + 1;}
else 
x[i] = x[i-1];
}
List result = List::create(x,k);
return(result);
}')

#generate random samples
set.seed(123)
N <- 1000
sigma <- 1
x0 <- 0
sample1 <- rw.Metropolis(sigma,x0,N)$x
sample2 <- rw_Metropolis(sigma,x0,N)[[1]]

#qq plot
library(car)
qqplot(sample1, sample2, xlab = "the samples using R",
       ylab = "the samples using Rcpp")
x <- seq(-4,4,.01)
lines(x,x,col = "red")

#Campare the computation time
library(microbenchmark)
ts <- microbenchmark(rw1 = rw.Metropolis(sigma,x0,N),rw2 = rw_Metropolis(sigma,x0,N))
summary(ts)[,c(1,3,5,6)]

