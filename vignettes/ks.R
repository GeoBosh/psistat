#' ---
#' title: "Computations for Lectures W7a-W8b - Practical Statistics I (2016/2017)"
#' subtitle: "Edited examples from the lectures"
#' author:   "Georgi N. Boshnakov"
#' date:     "March 2017"
#' output:
#'   html_document:
#'     toc: true
#'
#'   beamer_presentation:
#'     toc: true
#'     slide_level: 2
#' ---
#' <!--
#' %\VignetteEngine{knitr::rmarkdown}
#' %\VignetteIndexEntry{Supplementary materials}
#' -->
#'
#+ echo = FALSE
## based on the following documents from 2014/2015:
##    - QQandKSnew.html (R source lost)

#+ echo = FALSE
## package 'psistat' is discussed further below, so don't echo it here.
library(psistat)
set.seed(1234)

#' ## Packages and functions for topic "goodness of fit"
#' ### Functions in base R
#+ eval = FALSE
?ks.test
?ecdf
?qqnorm
?qqline

#' ### Package "psistat"

#' Package psistat provides some useful functions for this topic. You can install it on your
#' own computers from the web page of the course (if you have not done this yet).
#+ eval = FALSE
library(psistat)

#' The names of the functions in "psistat" start with "psi.". You can learn more about them
#' using the help system. For example:
#+ eval = FALSE
apropos("psi.*")
help(package="psistat")
package ? psistat
?psi.eqf
?psi.Dn
?psi.lks.exp.test
?psi.jitter

#' ### Other libraries providing relevant functions:
#+ eval = FALSE
library(nortest) # for Lilliefors test
?lillie.test

#+ eval = FALSE
library(e1071)  # for probplot
?probplot

#' Dataset Oxboys is in package "mlmRev"

#+ eval = FALSE
library(mlmRev)
# ?Oxboys

#' Miscellaneous: use this command to plot in a new window:
#+ eval = FALSE
dev.new()

#' Some random samples for the following examples.
x20 <- rnorm(20)            # a short random sample, to show it on screen
x200 <- rnorm(200)          # larger
x500 <- rnorm(500)          # even larger
x500ms <- rnorm(500, mean = 1, sd = 3)

x2 <- rchisq(100, df = 2)

xe <- rexp(500, rate = 1/1500) # life times of electric bulbs




#' ## Empirical cdf (ecdf)
x20.ecdf <- ecdf(x20)              # its ecdf
x20.ecdf
#' x20.ecdf is a function - we can use it as any other R function:
x20.ecdf(0)
x20.ecdf(10)

#' Here we plot the ecdf and overlay the cdf used to simulate the data. The sample is small,
#' so the two do not match very closely:
plot(x20.ecdf)
curve(pnorm, add = TRUE, col = "blue")

#' At the places were the ecdf jumps, the value is the one after the jump (indicated by the
#' filled points).

#' Here is a histogram with overlayed pdf for comparison.
hist(x20, freq = FALSE, ylim = c(0, 0.5))
curve(dnorm, add = TRUE)

#' For larger samples the ecdf is very good estimator of the underlying cdf:

x200.ecdf <- ecdf(x200)
plot(x200.ecdf)
curve(pnorm, add = TRUE, col = "red")

#' ... and an even larger sample:

x500.ecdf <- ecdf(x500)            # its ecdf
plot(x500.ecdf)                    # ecdf and cdf overlayed
curve(pnorm, add = TRUE, col = "red")
## This shows how to add the curve with non-default values of the
## arguments. This simply illustrates the technique:
norm2 <- function(x) pnorm(x, mean = 2)
curve(norm2, add = TRUE, col = "blue")

#' ## Empirical quantiles

# ?quantile
quantile(x500, probs=0.5)
quantile(x500, probs= c(0.25,0.5, 0.75))

quantile(x20, 0.5)
hist(x20, freq=FALSE, ylim = c(0, 0.5))
curve(dnorm, from = -3, to = 3,  add = TRUE, col = "blue")

#' ### Examples with psi.eqf

#' Empirical quantile function (from package "psistat")

#+ eval = FALSE
?psi.eqf


#' The eqf of a small sample:
x20.eqf <- psi.eqf(x20)
plot(x20.eqf, main = "Eqf of x20")
curve(qnorm, add = TRUE, col = "blue")

# plot(x200)
x200ecdf <- ecdf(x200)
x200ecdf(0)
quantile(x200, c(.25, .5, 0.75))

x200eqf <- psi.eqf(x200)
x200eqf(0.5)
x200eqf(0.25)
plot(x200eqf, main = "Eqf of x200")

x500.eqf <- psi.eqf(x500)          # a larger sample
plot(x500.eqf, main = "Eqf of x500")
curve(qnorm, add = TRUE, col = "blue")

#+ echo =FALSE
## #' Illustrate that quantiles can be plotted even without quantile f. (The "Weekend fun"
## #' question that stood for several weeks.)
## y500  <- pnorm(x500)
## plot(sort(x500), sort(y500))  # plot a cdf
##
## plot(sort(y500), sort(x500))  # plot the inverse cdf without using quantile function

#' ## QQ-plots
#' ### DIY qq-plots

#' We follow the procedure given in the notes: prepare points p,
#' evaluate the scores s,  compute the order statistics and plot.

#' Does x20 come from N(0,1)?
n <- length(x20)
n

p <- (1:n - 0.5)   # to check that it gives 0.5, 1.5, ..., 19.5
p

p <- (1:n - 0.5)/n
p

s <- qnorm(p)            # scores
x20.sorted <- sort(x20)  # order statistics
plot(s, x20.sorted)      # qq-plot
abline(0, 1, col = "blue") # plot the reference line


#' Does x2 come from N(0,1)?
n <- length(x2)
p <- ((1:n) - 0.5) / n
s2 <- qnorm(p)
x2s <- sort(x2)
plot(s2, x2s)
abline(0, 1, col = "red")

#' Does x200 come from N(0,1)?
n <- length(x200)
p <- ((1:n) - 0.5) / n
si <- qnorm(p)
x200.sorted <- sort(x200)
plot(si, x200.sorted)
abline(0, 1, col = "blue")

#' Does x200 come from Student-t with 3 d.f.?
si <- qt(p, df = 3)
plot(si, x200.sorted)
abline(0,1, col = "blue")

#' Does x500 come from N(0,1)?
n <- length(x500)                  # x500
p <- (1:n - 0.5)/n
s <- qnorm(p)
x500.sorted <- sort(x500)
plot(s, x500.sorted)
abline(0,1)                        # various ways to overlay a line
lm(x500.sorted ~ s)

abline(lm(x500.sorted ~ s), col = "red")
qqline(x500.sorted)

#' qq-plots for x500 The following examples are of qq-plots for x500
#' with various non-matching distributions.
s <- qt(p, df = 3)      # H0 is t_3 dist
plot(s, x500.sorted)
abline(0, 1, col="blue")

s <- qnorm(p, mean = 5) # H0 is N(5,1)
plot(s,x500.sorted)
abline(0, 1, col = "blue")

s <- qnorm(p, mean = 5, sd = 2) # H0 is N(5,2^2)
plot(s, x500.sorted)
abline(0, 1, col = "blue")

x500ms.sorted <- sort(x500ms)
p <- ((1:500) - 0.5)/500
s <- qnorm(p)
plot(s, x500ms.sorted)
abline(0, 1, col = "blue")

lmfit500 <- lm(x500ms.sorted ~ s)
summary(lmfit500)
abline(lmfit500, col = "red")


#' qq-plot example with an exponential distribution. First generate some data to work with:

#' xe might represent lifetime of bulbs (incandescent have typical expected life equal to
#' 1500 hours).
hist(xe)

n <- length(xe)
p <- ((1:n) - 0.5) / n
s <- qexp(p, rate = 1/1500)
plot(s, sort(xe))

#' ### Non-diy qq-plots

#' qqnorm produces a normal qq-plot, i.e. a qq-plot for hypothesised normal distribution:
qqnorm(x500ms)
abline(0,1)

qqnorm(x500)                   # qqnorm
qqnorm(xe)     # no match here

#' probplot() is another function for qq-plots. Without further
#' arguments it produces a normal plot:
library(e1071, quiet = TRUE)  # for probplot
probplot(x500)

#' For other hypothesised distributions the relevant quantile function should be given.
#' We often need to define the function ourselves.


#
# this doesn't work since probplot insists on naming the argument 'p'
# qexp1500 <- function(x) qexp(x, rate=1/1500)
#
                                        # here we oblige.
#' This function computes quantiles for exponential distribution with mean 1500:
qexp1500 <- function(p) qexp(p, rate = 1/1500)
#' (Note that probplot insists that the first argument is called 'p'.)
#' Compare:
qexp1500(0.5)
qexp(0.5, rate=1/1500)
#' Create the plot:
probplot(xe, qexp1500)

x2a <- rnorm(500, mean=5, sd=2)
qqnorm(x2a)

mean(x2a)
sd(x2a)

#' example: bulbs.txt
bulbs <- scan("bulbs.txt")
bulbs

qexp1500 <- function(p){ qexp(p, rate=1/1500) }
probplot(bulbs, qdist = qexp1500)

e100 <- rexp(100, rate=1/1500)
p <- ((1:100) - 0.5)/100
s <- qexp(p, rate=1/1500)
e100s <- sort(e100)
plot(s, e100s)
abline(0,1)

e1000 <- rexp(1000, rate=1/1500)
p <- ((1:1000) - 0.5)/1000
s <- qexp(p, rate=1/1500)
e1000s <- sort(e1000)
plot(s, e1000s)
abline(0,1)

s <- qnorm(p)
plot(s, e1000s)

library(e1071)
probplot(e1000, qexp)
probplot(e1000, qnorm)

# ?qexp

qchisq5 <- function(p) qchisq(p, df=5)
probplot(e1000,qchisq5)

#' ## Distributions of order statistics

#' Approx mean and variance of order statistics
#' example in Notes, p. 69
#'
#' first using formulae on pp. 68-69
qnorm(4/(19+1))

ex <- qnorm(4/(19+1))
ex

exvar <- 4*(19-4+1)/((19+1)^2*(19+2))/dnorm(ex)^2
exvar

#' ... then via simulation
y <- numeric(1000)
y[1] <- sort(rnorm(19))[4]
y[2] <- sort(rnorm(19))[4]
#' and so on until y[1000] ... :)

# ... but better use a single command
for(i in 1:1000) y[i] <- sort(rnorm(19))[4]

hist(y, freq=FALSE)  # estimate density

mean(y)              # estimate mean

var(y)               # estimate variance

#' this is an alternative to the 'for' loop
#'      (explain "replicate")
y <- replicate(1000, sort(rnorm(19))[4])
mean(y)
var(y)
quantile(y, c(0.25,0.5,0.75))

N <- 1000
x4 <- numeric(N)
x4[1] <- sort(rnorm(19))[4]
x4[1]
x4[2] <- sort(rnorm(19))[4]
# ... and so on; let's do it with a single command:
for(i in 1:N) x4[i] <- sort(rnorm(19))[4]

#' Explore the distribution of the fourth order statistic:
mean(x4)
var(x4)
hist(x4)

x4ecdf <- ecdf(x4)
curve(x4ecdf, from = -2, to = 2)

plot(x4ecdf)
curve(pnorm, add = TRUE, col = "red")

#' ## Inverse PIT

#' DIY generation of a random sample from distribution Expo(1/2)

#' First generate a sample from the uniform distribution.

u <- runif(8)
u

#' Evaluate the quantiles of the required distribution (Expo(1/2) here) for the values in the
#' U(0,1) random sample:
y <- -2*log(1-u)     #  DIY quantile function of Expo(1/2)
y

y1 <- qexp(u, rate = 1/2) # built-in quantile function
y1

#' The results are the same (so, qexp uses the formula -log(1-u)/lambda):
all(y == y1)

#' ## Kolmogorov-Smirnov tests

u <- runif(100)
x <- -1/2*log(1-u)
ks.test(x, "pexp", rate=2)

ks.test(x, "pexp", rate=10)

pexp2 <- function(x) pexp(x, rate = 2)
ks.test(x, pexp2)

psi.Dn(x, pexp2)

x <- runif(10)
x

y <- -1/2*log(1-x)
y

y1 <- qexp(x, rate=2)
y1

all(y==y1)  # TRUE (so, qexp uses the above formula

#' The cdf of the test statistic in KS test

#' See examples for psi.pks.

# ?psi.pks
psi.pks(0.6239385,4)

psi.pks(0.2940753,20)
psi.pks(0.1340279,100)
psi.pks(0.04294685,1000)

#' Compare the distribution of the test statistic for various sample sizes:

xi <- seq(0,1,length=100)            # some x values
plot(xi,psi.pks(xi,4))
lines(xi,psi.pks(xi,4))              # cdf of D_4
lines(xi,psi.pks(xi,50),col="blue")  # overlay the cdf of  D_{50}
lines(xi,psi.pks(xi,100),col="red")  # overlay the cdf of  D_{100}
abline(h=0.95, col="brown")
lines(xi,psi.pks(xi,1000),col="green")  # overlay the cdf of  D_{1000}

#' The abscissa of the intersection of the brown line with each of the curves gives the
#' corresponding critical value of the KS test. Notice that for larger samples the critical
#' value is smaller. In other words, smaller deviations from the hypothesised distribution
#' function are considered significant.

#' Similar to above using curve():

curve(psi.pks(x,4), from=0, to=1)                             # cdf of D_4
curve(psi.pks(x,10), from=0, to=1, col="blue", add=TRUE)      # cdf of D_10
curve(psi.pks(x,50), from=0, to=1, col="brown", add=TRUE)     # cdf of D_50
curve(psi.pks(x,100), from=0, to=1, col="red", add=TRUE)      # cdf of D_100
curve(psi.pks(x,1000), from=0, to=1, col="green", add=TRUE)   # cdf of D_1000
abline(h=0.95, col="brown")

#' ### DIY Dn

x10 <- rnorm(10)

#' diy Dn for H_0 = cdf of N(0,1)

u10 <- pnorm(x10)
u10

u10 <- sort(pnorm(x10))
u10

n <- 10
x10Dn <- max((1:n)/n - u10, u10 - (0:(n-1)/n))
x10Dn

#' now use psi.Dn to compute Dn, should give the same result.
psi.Dn(x10)

# KS test
#
# ?psi.Dn
x

psi.Dn(x)

#' The value of the statistic is the same as that from ks.test:
ks.test(x, pnorm)

psi.Dn(x) == ks.test(x, pnorm)$statistic

pe <- function(x) pexp(x, rate=1/1500)
ks.test(x, pe)


#' KS test 2013/2014 chunk

x20

x20os <- sort(x20)
pnorm(x20os, mean=3, sd=2)

pnorm(x20os, mean=3, sd=2) - (0:(20-1))/n


wrk1 <- pnorm(x20os, mean=3, sd=2) - (0:(20-1))/20
wrk2 <- (1:20)/20 - pnorm(x20os, mean=3, sd=2)
max(wrk1,wrk2)

Dn.x20 <- max(wrk1,wrk2)
Dn.x20

psi.Dn(x20)

# ?psi.Dn
psi.Dn(x20, cdf=pnorm, mean=3, sd=2)

ks.test(x20, pnorm, mean=3, sd=2)

# ?psi.pks
psi.pks(0.5, n=20)

#' ### Example migraine

#
# file.show("migraine.txt")
datamig <- scan("migraine.txt")
datamig

##  [1]  98  90 155  86  80  84  70 128  93  40 108  90 130  48  55 106 145
## [18] 126 100 115  75  95  38  66  63  32 105 118  21 142

ks.test(datamig, pnorm, mean=90, sd=35)

ks.test(unique(datamig), pnorm, mean=90, sd=35)

psi.Dn(unique(datamig), pnorm, mean=90, sd=35)


#' Example from help page of psi.pks

xi <- seq(0,1,length=100)            # some x values
plot(xi, psi.pks(xi,4))               # cdf of D_4

#' diy Dn

x <- unique(datamig)
x

xs <- sort(x)
xs

n <- length(xs)
max( (1:n)/n - pnorm(xs,mean=90,sd=35), pnorm(xs,mean=90,sd=35) - (0:(n-1))/n)

#' Dn using psi.Dn

psi.Dn(unique(datamig), pnorm, mean=90, sd=35)

ks.test(unique(datamig),pnorm, mean=90, sd=35)

migDn <- max( (1:n)/n - pnorm(xs,mean=90,sd=35), pnorm(xs,mean=90,sd=35) - (0:(n-1))/n)
migDn

psi.pks(migDn, n)

1 - psi.pks(migDn, n)

#' ## Lilliefors test for normality

#+ eval = FALSE
apropos("psi.")
?psi.pkls.exp
?psi.plks.exp
?psi.qlks.exp

#' This package provides lillie.test():

library(nortest)
# ?lillie.test
lillie.test(x20)


#' ## Further examples

#' ### Example: moths

# file.show("mothsontrees.txt")
datamoths <- scan("mothsontrees.txt")
datamoths

# ?punif
ks.test(datamoths, punif, min=0, max=25)

#' Alternatively:

mycdf <- function(q) punif(q, min=0, max=25)
ks.test(datamoths, mycdf)

mothssorted <- sort(datamoths)

val <- punif(mothssorted, min=0, max=25)
length(datamoths)

n <- 14
(1:n)/n

(1:n)/n - val

val - (0:(n-1))/n

max( (1:n)/n - val, val - (0:(n-1))/n )

ks.test(datamoths, "punif", min=0, max=25)

f1 <- function(x) psi.pks(x,14)
curve(f1, from=0.01, 0.99)

xi <- seq(0,1, 0.01)
head(cbind(xi, "f1(xi)" = f1(xi)))

ks.test(datamoths, "punif", min=0, max=25)

d14 <- max( (1:n)/n - val, val - (0:(n-1))/n )
d14

1 - psi.pks(d14,14)

ks.test(datamoths, "punif", min=0, max=25)


#' ### Example: datasales

datasales <- c(2,4,8,18,9,11,13)
ks.test(datasales, "pnorm", mean=10, sd=3)

#' ### Example: barbiturate

# file.show("barbiturate.txt")
databarbi <- scan("barbiturate.txt")
databarbi

lillie.test(databarbi)

ks.test(databarbi,pnorm,mean=mean(databarbi), sd=sd(databarbi))


#' ### Example: bulbs1000

databulbs <- scan("bulbs1000.txt")
databulbs

ks.test(databulbs, "pexp", rate=2)

ks.test(databulbs, "pexp", rate=1/1000)

mean(databulbs)
1/mean(databulbs)

# apropos("psi.")
# ?psi.lks.exp.test
psi.lks.exp.test(databulbs)


#' ### Example: bulbs1500

databulbs15 <- scan("bulbs1500.txt")
psi.lks.exp.test(databulbs15)

psi.lks.exp.test(databulbs15)

psi.lks.exp.test(databulbs15, Nsim=10000)


# example: bulbs: is this the same as bulbs1000?
#
bulbs
mean(bulbs)

ks.test(bulbs, pexp, rate = 1/1500)
ks.test(bulbs, pexp, rate = 1/1000)


#' Carry out Lilliefors test
#' (simulation is used, so slight differences in repeated calculation)
psi.lks.exp.test(bulbs)
psi.lks.exp.test(bulbs)

psi.lks.exp.test(bulbs, Nsim = 10000)   # for more precision

#' (semi-)diy (full diy would also calc. the Dn stat. by diy)
z <- bulbs/mean(bulbs)
zDn <- psi.Dn(z, pexp, rate=1)

#' crit. value at alpha=0.05
DnN0p05 <- psi.qlks.exp(1-0.05, length(z))
#' or, for more precision,
DnN0p05 <- psi.qlks.exp(1-0.05, length(z), Nsim=10000)

zDn > DnN0p05

# p-value
1 - psi.plks.exp(zDn, length(z), Nsim=10000)

#' ### Example: Oxboys (needs clean-up)

library(mlmRev)
summary(Oxboys)          # Oxboys are from library(mlmRev)

Oxboys[1:5,]
Oxboys[1:5, "height"]

dataoxheight <- Oxboys[,"height"]
plot(dataoxheight)

summary(dataoxheight)

hist(dataoxheight, freq=FALSE)  # todo: overlay exp pdf?

boxplot(dataoxheight)

qqnorm(dataoxheight)

library(nortest)
lillie.test(dataoxheight)


ks.test(dataoxheight, pnorm, mean=mean(dataoxheight), sd=sd(dataoxheight)) # for comparison

ks.test(unique(dataoxheight), "pnorm", mean=mean(dataoxheight), sd=sd(dataoxheight))

# add some jitter to remove ties;  ?psi.jitter
# ?psi.jitter
ks.test(psi.jitter(dataoxheight,amount=0.5), pnorm, mean=mean(dataoxheight), sd=sd(dataoxheight))

ks.test(psi.jitter(dataoxheight,amount=0.5), pnorm, mean=mean(dataoxheight), sd=sd(dataoxheight))

#?psi.Dn
psi.Dn(dataoxheight, pnorm, mean=mean(dataoxheight), sd=sd(dataoxheight), ties.jitter=TRUE)

# ?psi.jitter
any(duplicated(dataoxheight))

xj <- psi.jitter(dataoxheight, 0.5)
any(duplicated(xj))

ks.test(xj, "pnorm", mean=mean(dataoxheight), sd = sd(dataoxheight))


###########################################

datamig

ls(pattern="data*")

databarbi
databulbs

#' ### Example: checking normality of residuals from lm()

datamilk <-
   data.frame(
     x = c(42.7,40.2,38.2,37.6,32.2,32.2,28,27.2,26.6,23,22.7,21.8,21.3,20.2),
     y = c(1.2,1.16,1.07,1.13,0.96,1.07,0.85,0.87,0.77,0.74,0.76,0.69,0.72,0.64)
   )
fitmilk <- lm(y~x, data=datamilk)
datamilk

fitmilk
summary(fitmilk)

# plot(fitmilk)
resfitmilk <- residuals(fitmilk)
lillie.test(resfitmilk)

#' todo: also qq-plot? ## Lilliefors test for exponentiality
#+ eval = FALSE
?psi.plks
?psi.plks.exp
?psi.lks.exp.test
?psi.Dn
?psi.pks

#'
databulbs <- scan("bulbs1000.txt")
ks.test(databulbs, "pexp", rate=2)

mean(databulbs)

1/mean(databulbs)

ks.test(databulbs, "pexp", rate=1/1000)

psi.plks.exp(0.203, df = length(databulbs))

1 - psi.plks.exp(0.203, df = length(databulbs))

curve(psi.pks(x,length(databulbs)), from=0, to=1)
curve(psi.plks.exp(x,length(databulbs)), from=0, to=1, col="blue", add=TRUE)
abline(h=0.95)

# ?names
# ?c
databulbs


ks.test(databulbs, "pexp", rate=1/1000)

psi.Dn(databulbs, "pexp", rate=1/1000)

#' q : P(Dn < q) = 0.95 for selected values of n

psi.pks(0.6239385,     4)
psi.pks(0.2940753,    20)
psi.pks(0.1340279,   100)
psi.pks(0.04294685, 1000)   # asymptotic approximation

databulbs
psi.lks.exp.test(databulbs)

length(databulbs)

psi.plks.exp(0.203, df=10)

1 - psi.plks.exp(0.203, df=10)

