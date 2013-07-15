# ------------------------------------------------
# define PSD model 
#
# Additional parameters that are not to be fitted but
# passed to the model functions are stored in the MOD
# array.

 model <- function(par, x, mod=c(0, 1.0)) {

   a.low <- mod[2]

   if (mod[1] == 0) { 
     y.mod <- model0(par, x) 

   } 

   if (mod[1] == 1) {
     y.mod <- model1(par, x, a.low)
   }

   if (mod[1] == 2) {
	y.mod <- model2(par, x)
   }

     if (mod[1] == 3) {
	y.mod <- model3(par, x, a.low)
   }

   if (mod[1] == 4) {
        y.mod <- model4(par, x)
   }

     if (mod[1] == 5) {
        y.mod <- model5(par, x)
   }

   if (mod[1] == 6) {
        y.mod <- model6(par, x)
   }

   if (mod[1] == 7) {
        y.mod <- model7(par, x)
   }


 return(y.mod)
}


#
# define PSD model as QPO
#
#
model4 <- function(par, x){

  gamma <- par[1]
  norm <- par[2]
  nu0 <- par[3]

  alpha = exp(norm)*exp(gamma)/(pi*2.0)
  y = alpha/((x - exp(nu0))**2.0 + exp(gamma)**2.0)

  return(y)
}


model5 <- function(par, x){

  ypl <- model0(c(par[1], par[2], par[3]), x)

  yqpo <- model4(c(par[4], par[5], par[6]), x)

  y <- ypl + yqpo

  return(y)
}

model6 <- function(par, x){

  ybpl <- model0(c(par[1], par[2], par[3], par[4], par[5]), x)

  yqpo <- model4(c(par[6], par[7], par[6]), x)

  y <- ybpl + yqpo

  return(y)
}

model7 <- function(par, x){

  yqpo <- model4(c(par[1], par[2], par[3]), x)

  y <- yqpo + par[4]

  return(y)
}





# ------------------------------------------------
# define PSD model as a power law y = b*x^-a + c

model0 <- function(par, x) { 
  logmin <- -100
  a <- par[1]
  b <- par[2]
  c <- par[3]

  ly.mod <- -a*log(x)+b
  ly.mod <- pmax(ly.mod,logmin)
  y.mod <- exp(ly.mod) + exp(c)
  return(y.mod) 
}

# ------------------------------------------------
# define PSD model as a bending power law 
#
# y(x) = b*x^(-a1) / [1 + (x/x_b)^(a2-a1) ] + c
#
# y(x >> x_b) = A * x^(-a2) * x_b^(a2-a1)
# y(x << x_b) = A * x^(-a1)
#
# This is a power law of slope -a1 at x << x_b
# and a power law of slope -a2 at x >> x_b.
# Between these two extremes is smoothly changes
# from one slope to the other.
# (x_b is the break/bend point.)
#
# The parameters are stored in the one dimensional array par
# par = [a2, x_b, b, c]
#      a1  = power law index at x << x_b [given seperately]
#      a2  = power law index at x >> x_b
#      x_b = break point
#      b   = normalisation (input as log[N])
#      c   = constant
#
# x_b, b and c are passed to the function in log_10 units.

model1 <- function(par, x, a.low=1.0) {

  logmin <- -100


  a1 <- par[1]
  a2 <- par[3]
  x_b <- par[4]
  b <- par[2]
  c <- par[5]

# convert x to log units

  logx <- log(x)

# set up arrays

  logq <- x*0
  y <- x*0

# Calculate bending factor q=1+z with z = (x/x_b)^(a2-a1).
# But to avoid overflow/underflow we calculate
#   log[z] = (a2 - a1) * (log[x] - log[x_b])

  logz <- (a2-a1)*(logx - x_b)

# split into regions of low, medium and high z.
# i.e. where x << x_b, x ~ x_b and x >> x_b
# and treat these seperately to avoid underflow/
# overflow errors

  lo <- (logz < -16)
  hi <- (logz > 16)
  me <- (logz > -16) & (logz < 16)

  if (sum(lo, na.rm=TRUE) > 0) { logq[lo] <- log(1.0) }
  if (sum(hi, na.rm=TRUE) > 0) { logq[hi] <- logz[hi] }
  if (sum(me, na.rm=TRUE) > 0) { logq[me] <- log(exp(logz[me])+1) }

# calculate log(y)

  logy <- -a1*logx - logq + b

# watch for very low/high values (catch over/underflow)

  lo <- (logy < logmin)
  hi <- (logy > -logmin)
  me <- (logy > logmin) & (logz < -logmin)

  if (sum(hi, na.rm=TRUE) > 0) { y[hi] <- exp(-logmin) }
  if (sum(me, na.rm=TRUE) > 0) { y[me] <- exp(logy[me]) }

  y <- y + exp(c)

  return(y)

}


# ------------------------------------------------
# define PSD model as a power law y = b*x^-a 

model2 <- function(par, x) { 
  logmin <- -100
  a <- par[1]
  b <- par[2]

  ly.mod <- -a*log(x)+b
  ly.mod <- pmax(ly.mod,logmin)
  y.mod <- exp(ly.mod) 
  return(y.mod) 
}

# ------------------------------------------------
# define PSD model as a bending power law 
#
# y(x) = b*x^(-a1) / [1 + (x/x_b)^(a2-a1) ] 
#
# y(x >> x_b) = A * x^(-a2) * x_b^(a2-a1)
# y(x << x_b) = A * x^(-a1)
#
# This is a power law of slope -a1 at x << x_b
# and a power law of slope -a2 at x >> x_b.
# Between these two extremes is smoothly changes
# from one slope to the other.
# (x_b is the break/bend point.)
#
# The parameters are stored in the one dimensional array par
# par = [a2, x_b, b]
#      a1  = power law index at x << x_b [given seperately]
#      a2  = power law index at x >> x_b
#      x_b = break point
#      b   = normalisation (input as log[N])
#
# x_b, b are passed to the function in log_10 units.

model3 <- function(par, x, a.low=1.0) {

  logmin <- -100

  a1 <- a.low
  a2 <- par[1]
  x_b <- par[2]
  b <- par[3]

# convert x to log units

  logx <- log(x)

# set up arrays

  logq <- x*0
  y <- x*0

# Calculate bending factor q=1+z with z = (x/x_b)^(a2-a1).
# But to avoid overflow/underflow we calculate
#   log[z] = (a2 - a1) * (log[x] - log[x_b])

  logz <- (a2-a1)*(logx - x_b)

# split into regions of low, medium and high z.
# i.e. where x << x_b, x ~ x_b and x >> x_b
# and treat these seperately to avoid underflow/
# overflow errors

  lo <- (logz < -16)
  hi <- (logz > 16)
  me <- (logz > -16) & (logz < 16)

  if (sum(lo, na.rm=TRUE) > 0) { logq[lo] <- log(1.0) }
  if (sum(hi, na.rm=TRUE) > 0) { logq[hi] <- logz[hi] }
  if (sum(me, na.rm=TRUE) > 0) { logq[me] <- log(exp(logz[me])+1) }

# calculate log(y)

  logy <- -a1*logx - logq + b

# watch for very low/high values (catch over/underflow)

  lo <- (logy < logmin)
  hi <- (logy > -logmin)
  me <- (logy > logmin) & (logz < -logmin)

  if (sum(hi, na.rm=TRUE) > 0) { y[hi] <- exp(-logmin) }
  if (sum(me, na.rm=TRUE) > 0) { y[me] <- exp(logy[me]) }

  return(y)

}


# ------------------------------------------------
# define periodogram log likelihood function
# S = -log(likelihood)

mlogl <- function(par, x, y, mod=c(0, 1.0)) {
  mody <- model(par, x, mod=mod)
  l <- sum( log(mody) + y/mody , na.rm=TRUE)
  return(l)
}

# ------------------------------------------------
# define posterior in terms of likelihood * prior
# for the parameters. Actually we work with  
# -log[posterior] = -log[likelihood] + -log[prior]

lpost <- function(par, x, y, mod=c(0, 1.0)) {
  ml.post <- mlogl(par, x, y, mod=mod) + mlprior(par, mod=mod)
  return(ml.post)
}


# ------------------------------------------------
# Define the prior densities for the parameter
# of the model. Actually, calculate the minus
# log prior density, which can then be combined
# with the minus log likelihood (MLOGL).

mlprior <- function(par, mod=c(0, 1.0)) {

# smallest allowed log(number) - to avoid underflow
# and -log(0) infinities.

  logmin <- -100

# set allowed range of parameters

  alim <- c(-1,8)

# define prior density at parameters = par
# as the product of the prior densities 
# for each parameter

  if ( mod[1] == 0 ) {
     a <- par[1]
     b <- par[2]
     c <- par[3]
     pa <- (a >= min(alim) && a <= max(alim))
     pb <- 1.0
     pc <- 1.0
     pr <- pa*pb*pc
  } 

  if (mod[1] == 1) {
     a <- par[1]
     b <- par[2]
     c <- par[3]
     d <- par[4]
     e <- par[5]
     pa <- (a >= min(alim) && a <= max(alim))
     pb <- 1.0
     pc <- 1.0
     pd <- 1.0
     pe <- 1.0

# Alternative priors using Normal densities (for RE J1034)
#     pa <- dnorm(a,mean=2,sd=2)
#     pb <- dnorm(b,mean=-6.9,sd=2.3)
#     pc <- dnorm(c,mean=-4.6,sd=2.3)
#     pd <- dnorm(d,mean=0,sd=2.3)

     pr <- pa*pb*pc*pd*pe
  }

  if (mod[1] == 2) {
     a <- par[1]
     b <- par[2]
     pa <- (a >= min(alim) && a <= max(alim))
     pb <- 1.0
     pr <- pa*pb
  } 

  if (mod[1] == 3) {
     a <- par[1]
     b <- par[2]
     c <- par[3]
     pa <- (a >= min(alim) && a <= max(alim))
     pb <- 1.0
     pc <- 1.0

     pr <- pa*pb*pc
  }

  if (mod[1] == 4) {
     a <- par[1]
     b <- par[2]
     c <- par[3]
     pa <- 1.0
     pb <- 1.0
     pc <- 1.0

     pr <- pa*pb*pc
  }

  if (mod[1] == 5) {
     a <- par[1]
     b <- par[2]
     c <- par[3]
     d <- par[4]
     e <- par[5]
     f <- par[6]

     pa <- (a >= min(alim) && a <= max(alim))
     pb <- 1.0
     pc <- 1.0
     pd <- 1.0
     pe <- 1.0
     pf <- 1.0
     pr <- pa*pb*pc*pd*pe*pf
  }

  if (mod[1] == 6) {
     a <- par[1]
     b <- par[2]
     c <- par[3]
     d <- par[4]
     e <- par[5]
     f <- par[6]
     g <- par[7]
     h <- par[8]
     pa <- (a >= min(alim) && a <= max(alim))
     pb <- 1.0
     pc <- 1.0
     pd <- 1.0
     pe <- 1.0
     pf <- 1.0
     pg <- 1.0
     ph <- 1.0


# Alternative priors using Normal densities (for RE J1034)
#     pa <- dnorm(a,mean=2,sd=2)
#     pb <- dnorm(b,mean=-6.9,sd=2.3)
#     pc <- dnorm(c,mean=-4.6,sd=2.3)
#     pd <- dnorm(d,mean=0,sd=2.3)

     pr <- pa*pb*pc*pd*pe*pf*pg*ph
  }


# take the log, watching out for zeros

  if (pr > 0) {
    mlp <- -log(pr)
  } else {
    mlp <- -logmin
  }

  return(mlp)

}
1

