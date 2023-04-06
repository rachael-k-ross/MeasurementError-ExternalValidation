#####################################################
#
# Code for generation and analysis of a single 
# simulated dataset from 
# Leveraging external validation data: the challenge of transporting measurement error parameters 
# By ### et al.
#
# Code by ### (2023/01/23)
#
#####################################################

library("tidyverse")
library("broom")
library("purrr")
library("rootSolve")
library("numDeriv")
library("geepack")

####################################################
# DEFINE FUNCTIONS FOR DATA GENERATION AND ANALYSIS
####################################################

# Function to generate data
genbifx <- function(pop,N,Zp,Wp,Ap,Y0p,Y1p,Ystarp){
  tibble( 
    id = c(1:N),
    r  = pop,
    u_a = rbinom(N, size=1, prob=.5),
    u_y = rbinom(N, size=1, prob=.5),
    z  = rnorm(N, mean=Zp[1], sd=Zp[2]),
    w  = rbinom(N, size=1, prob=plogis(Wp[1] + Wp[2]*u_a + Wp[3]*u_y)),
    a  = rbinom(N, size=1, prob=plogis(Ap[1] + Ap[2]*z + Ap[3]*u_a)),
    y0 = rbinom(N, size=1, prob=plogis(Y0p[1] + Y0p[2]*z + Y0p[3]*w)),
    y1 = rbinom(N, size=1, prob=plogis(Y1p[1] + Y1p[2]*z + Y1p[3]*w + Y1p[4]*z + Y1p[5])),
    y = a*y1 + (1-a)*y0,
    pystar = plogis(Ystarp[1] + Ystarp[2]*y + Ystarp[3]*z + + Ystarp[4]*a + Ystarp[5]*w),
    ystar = rbinom(N, size=1, prob=plogis(Ystarp[1] + Ystarp[2]*y + Ystarp[3]*z + + Ystarp[4]*a + Ystarp[5]*w)))
}

# Function to set input parameters for each scenario and generate data
togenbi <- function(scenario,n1size,n0size){ 
  
  if(scenario=="Dalt"){
    u_OR = 4.0
  } else {
    u_OR = 1.2
  }
  
  # Parameters shared by both pops
  if(scenario=="D"){
    wcoef <- c(log(1.2),log(1.2)) 
  } else {
    wcoef <- c(log(1),log(1))
  }
  acoef <- c(log(0.6),log(1.2)) 
  y0coef <- c(log(1.15),log(1.2)) 
  y1coef <- c(y0coef,log(0.8),0.275) 

  # Parameters unique to main study
  n1 <- n1size
  zparm1 <- c(0,1)
  if(scenario=="D"){
    wparm1 <- c(0.514,wcoef) 
  } else {
    wparm1 <- c(0.693,wcoef)
  }
  aparm1 <- c(-1.555, acoef) 
  y0parm1 <- c(-1.364, y0coef) 
  y1parm1 <- c(y0parm1[1], y1coef)
  
  # Parameters unique to validation data
  n0 <- n0size
  zparm0 <- c(1,1)
  if(scenario=="D"){
    wparm0 <- c(-0.878,wcoef) 
  } else {
    wparm0 <- c(-0.693,wcoef)
  }
  aparm0 <- c(0.42, acoef) 
  y0parm0 <- c(-2.333, y0coef) 
  y1parm0 <- c(y0parm0[1], y1coef)
  
  # ME parameters (shared by both pops)
  if(scenario=="A"){
    deltas <- c(-1.735, 3.932, 0, 0, 0)
  } else if(scenario=="B") {
    deltas <- c(-1.853, 3.979, log(0.9), log(1.5), 0)
  } else if(scenario=="C") {
    deltas <- c(-1.946, 3.990, log(0.9), log(1.5), log(1.3))  
  } else {
    deltas <- c(-1.946, 3.988, log(0.9), log(1.5), log(1.3))
  }
  
  main <- genbifx(1,n1,zparm1,wparm1,aparm1,y0parm1,y1parm1,deltas)
  validation <- genbifx(0,n0,zparm0,wparm0,aparm0,y0parm0,y1parm0,deltas)
  full <- rbind(main,validation)
  return(full)
}

#### Functions to implement M-estimation

# Function sums each column of estimating fx (b/c estimating fx produces a n by p matrix)
sumstack <- function(parms, estfx){ 
  sums <- colSums(estfx(parms))
  return(sums)
}

# Function estimaties theta 
solvetheta <- function(estfx, init){
  ests <- multiroot(f=sumstack, start=init, estfx=estfx, atol=1e-12)$root
  return(ests)
}

# Function estimates covariance matrix of theta
mysandwich <- function(estfx, thetahat){
  efhat <- estfx(thetahat)
  n <- nrow(efhat)
  meat <- t(efhat) %*% efhat / n
  bread <- -numDeriv::jacobian(sumstack, thetahat, estfx=estfx)/n
  sandwich <- (solve(bread)%*%meat%*%t(solve(bread)))/n
  stderr <- sqrt(diag(sandwich))
  return(stderr)
}

# Implement full M-estimation
mest <- function(estfx, init, conf=0.95){
  thetahat <- solvetheta(estfx, init=init)
  ses <- mysandwich(estfx, thetahat)
  results <- tibble(est = thetahat,
                    se = ses,
                    lcl = est + qnorm((1-conf)/2)*se,
                    ucl = est + qnorm((1+conf)/2)*se)
  return(results)
}


### Functions for estimating point estimates by MLE

# Function for g-comp estimates 
getests <- function(data,b_est){
  nc <- data %>% filter(r==1) %>% mutate(py=plogis(b_est[1] + b_est[2]*a + b_est[3]*z + b_est[4]*a*z)) %>% pull(py)
  exp <- data %>% filter(r==1) %>% mutate(a=1) %>% mutate(py=plogis(b_est[1] + b_est[2]*a + b_est[3]*z + b_est[4]*a*z)) %>% pull(py)
  unexp <- data %>% filter(r==1) %>% mutate(a=0) %>% mutate(py=plogis(b_est[1] + b_est[2]*a + b_est[3]*z + b_est[4]*a*z)) %>% pull(py)
  c(rd=mean(exp)-mean(unexp), all1=mean(exp), all0=mean(unexp), nc=mean(nc))
}

# Function for modified likelihood
loglik <- function(b, data, se, sp){
  mu <- with(data,plogis(b[1] + b[2]*a + b[3]*z + b[4]*a*z))
  ilogL <- with(data,ifelse(ystar==1, log(se*mu + (1-sp)*(1-mu)), log((1-se)*mu + sp*(1-mu))))
  neglogL <- -1*sum(ilogL)
  return(neglogL)
}


#### Estimating equations for each estimator

# Standard gcomp
gc_ef <- function(theta){
    
    # Extract individual parameters
    risk_diff <- theta[1]
    risk_a1 <- theta[2]
    risk_a0 <- theta[3]
    risk_nc <- theta[4]
    beta <- theta[5:length(theta)]
    
    # Nuisance outcome model (logistic regression)
    ee_outmod <- as.vector(out - plogis(xmat %*% beta))*xmat*r
    
    # Generating predicted outcome values
    ee_a1 <- (plogis(xmat1 %*% beta) - risk_a1)*r
    ee_a0 <- (plogis(xmat0 %*% beta) - risk_a0)*r
    ee_nc <- (plogis(xmat %*% beta) - risk_nc)*r
    
    # Contrast
    ee_rd <- as.vector(rep(risk_a1 - risk_a0, length(out)) - risk_diff)*r
    
    return(cbind(ee_rd, ee_a1, ee_a0, ee_nc, ee_outmod))
}

# Accounting for nondifferential outcome mislcassification
gc_ndme_ef <- function(theta){
    
    colx <- ncol(xmat)
    
    # Extract individual parameters
    risk_diff <- theta[1]
    risk_a1 <- theta[2]
    risk_a0 <- theta[3]
    risk_nc <- theta[4]
    beta <- theta[5:(4 + colx)]
    gamma <- theta[(5 + colx):length(theta)]
    
    # ME process in validation data
    ee_g1 <- (1-r)*y*(ystar - gamma[1])
    ee_g2 <- (1-r)*(1-y)*(ystar - gamma[2])
    
    # Nuisance outcome model (using modified log-likelihood accounting for ME)
    mu <- plogis(xmat %*% beta)
    S1 <- ystar*(gamma[1] - gamma[2])/(gamma[2] + mu*(gamma[1] - gamma[2]))
    S2 <-(1-ystar)*(gamma[1] - gamma[2])/((1-gamma[2]) - mu*(gamma[1] - gamma[2]))
    S <- S1-S2
    ee_outmod <- as.vector(mu*(1-mu)*S)*xmat*r
    
    # Generating predicted outcome values
    ee_a1 <- (plogis(xmat1 %*% beta) - risk_a1)*r
    ee_a0 <- (plogis(xmat0 %*% beta) - risk_a0)*r
    ee_nc <- (plogis(xmat %*% beta) - risk_nc)*r
    
    # Contrast
    ee_rd <- as.vector(rep(risk_a1 - risk_a0, nrow(xmat)) - risk_diff)*r
    
    return(cbind(ee_rd, ee_a1, ee_a0, ee_nc, ee_outmod, ee_g1, ee_g2))
}

# Accounting for differential outcome misclassificaton WRT A and Z
gc_dme_ef <- function(theta){
    
    colx <- ncol(xmat)
    
    # Extract individual parameters
    risk_diff <- theta[1]
    risk_a1 <- theta[2]
    risk_a0 <- theta[3]
    risk_nc <- theta[4]
    beta <- theta[5:(4 + colx)]
    delta <- theta[(5 + colx):length(theta)]
    
    # Fit logistic model for ME process in validation data
    ee_mep <- as.vector(ystar - plogis(xmatdag %*% delta))*xmatdag*(1-r)
    
    # Obtain individual level Se/Sp 
    gamma1 <- plogis(xmatdag1 %*% delta)
    gamma2 <- plogis(xmatdag0 %*% delta)
    
    # Nuisance outcome model (using modified log-likelihood accounting for ME)
    mu <- plogis(xmat %*% beta)
    S1 <- ystar*(gamma1 - gamma2)/(gamma2 + mu*(gamma1 - gamma2))
    S2 <-(1-ystar)*(gamma1 - gamma2)/((1-gamma2) - mu*(gamma1 - gamma2))
    S <- S1-S2
    ee_outmod <- as.vector(mu*(1-mu)*S)*xmat*r
    
    # Generating predicted outcome values
    ee_a1 <- (plogis(xmat1 %*% beta) - risk_a1)*r
    ee_a0 <- (plogis(xmat0 %*% beta) - risk_a0)*r
    ee_nc <- (plogis(xmat %*% beta) - risk_nc)*r
    
    # Contrast
    ee_rd <- as.vector(rep(risk_a1 - risk_a0, nrow(xmat)) - risk_diff)*r
    
    return(cbind(ee_rd, ee_a1, ee_a0, ee_nc, ee_outmod, ee_mep))
}

# Accounting for differential outcome misclassificaton WRT A, Z & W by conditioning on W
# Use function above

# Accounting for differential outcome misclassificaton WRT A, Z & W by weighting misclass parameters
gc_dme_wtme_ef <- function(theta){
    
    colx <- ncol(xmat)
    colxdag <- ncol(xmatdag)
    colxddag <- ncol(xmatddag)
    
    # Extract individual parameters
    risk_diff <- theta[1]
    risk_a1 <- theta[2]
    risk_a0 <- theta[3]
    risk_nc <- theta[4]
    beta <- theta[5:(4 + colx)]
    delta <- theta[(5 + colx):(4 + colx + colxdag)]
    nu <- theta[(5 + colx + colxdag):(4 + colx + colxdag + colxddag)]
    phi <- theta[(5 + colx + colxdag + colxddag):length(theta)]
    
    #Fit logistic selection model for stabilization model
    ee_rmod_stab = as.vector(r - plogis(xmatdia %*% phi))*xmatdia
    
    #Fit logistic selection model for weights
    ee_rmod = as.vector(r - plogis(xmatddag %*% nu))*xmatddag
    
    #Construct weights
    pi = (plogis(xmatddag %*% nu)*(1-plogis(xmatdia %*% phi)))/((1-plogis(xmatddag %*% nu))*plogis(xmatdia %*% phi))
    
    # Fit weighted logistic model for ME process in validation data
    ee_mep <- as.vector((ystar - plogis(xmatdag %*% delta))*pi)*xmatdag*(1-r)
    
    # Obtain individual level Se/Sp 
    gamma1 <- plogis(xmatdag1 %*% delta)
    gamma2 <- plogis(xmatdag0 %*% delta)
    
    # Nuisance outcome model (using modified log-likelihood accounting for ME)
    mu <- plogis(xmat %*% beta)
    S1 <- ystar*(gamma1 - gamma2)/(gamma2 + mu*(gamma1 - gamma2))
    S2 <-(1-ystar)*(gamma1 - gamma2)/((1-gamma2) - mu*(gamma1 - gamma2))
    S <- S1-S2
    ee_outmod <- as.vector(mu*(1-mu)*S)*xmat*r
    
    # Generating predicted outcome values
    ee_a1 <- (plogis(xmat1 %*% beta) - risk_a1)*r
    ee_a0 <- (plogis(xmat0 %*% beta) - risk_a0)*r
    ee_nc <- (plogis(xmat %*% beta) - risk_nc)*r
    
    # Contrast
    ee_rd <- as.vector(rep(risk_a1 - risk_a0, nrow(xmat)) - risk_diff)*r
    
    return(cbind(ee_rd, ee_a1, ee_a0, ee_nc, ee_outmod, ee_mep, ee_rmod, ee_rmod_stab))
}

####################################################
# IMPLEMENTATION
####################################################

#### Generate data
dat <- togenbi("C",100000,2000)

#### Data set up
r <- dat$r
y <- ifelse(dat$r==1,0,dat$y) # make y 0 for everyone in main study
ystar <- dat$ystar
ones <- rep(1,nrow(dat)) 

xmat <- matrix(c(ones, dat$a, dat$z, dat$z*dat$a),nrow=nrow(dat)) # design matrix for out model
xmat1 <- matrix(c(ones, ones, dat$z, dat$z),nrow=nrow(dat)) # x set to 1 for all
xmat0 <- matrix(c(ones, 1-ones, dat$z, 1-ones),nrow=nrow(dat)) # x set to 0 for all

xmatdag <- matrix(c(ones, dat$y, dat$a, dat$z),nrow=nrow(dat)) # design matrix for misclass model
xmatdag1 <- matrix(c(ones, ones, dat$a, dat$z),nrow=nrow(dat)) # y set to 1 for all
xmatdag0 <- matrix(c(ones, 1-ones, dat$a, dat$z),nrow=nrow(dat)) # y set to 0 for all


#### Analysis 0: using true outcome
out <- dat$y
mest(estfx=gc_ef, init=rep(0,8))


#### Analysis 1: using misclassified outcome
out <- dat$ystar
mest(estfx=gc_ef, init=rep(0,8))


#### Analysis 2: accounting for nondifferential misclassification
mest(estfx=gc_ndme_ef,init=c(rep(0,8),.9,.1))


#### Analysis 3: accounting for differential misclassification WRT A & Z

# Fit mislcassification model in validation data to get good starting values
modsesp <- glm(ystar ~ y + a + z, data=dat[dat$r==0,], family='binomial')

# Implement M-estimation
mest(estfx=gc_dme_ef, init=c(rep(0,8),modsesp$coefficients))


#### Analysis 4: accounting for differential misclassification WRT A, Z, & W by conditioning on W

# Update data set up - add W to matrices
xmat <- matrix(c(ones, dat$a, dat$z, dat$z*dat$a, dat$w),nrow=nrow(dat)) # design matrix
xmat1 <- matrix(c(ones, ones, dat$z, dat$z, dat$w),nrow=nrow(dat)) # x set to 1 for all
xmat0 <- matrix(c(ones, 1-ones, dat$z, 1-ones, dat$w),nrow=nrow(dat)) # x set to 0 for all

xmatdag <- matrix(c(ones, dat$y, dat$a, dat$z, dat$w),nrow=nrow(dat)) # design matrix for misclass model
xmatdag1 <- matrix(c(ones, ones, dat$a, dat$z, dat$w),nrow=nrow(dat)) # y set to 1 for all
xmatdag0 <- matrix(c(ones, 1-ones, dat$a, dat$z, dat$w),nrow=nrow(dat)) # y set to 0 for all

# Fit mislcassification model in validation data to get good starting values
modsesp <- glm(ystar ~ y + a + z + w, data=dat[dat$r==0,], family='binomial')

# Implement M-estimation
mest(estfx=gc_dme_ef, init=c(rep(0,9),modsesp$coefficients))


#### Analysis 5: accounting for differential misclassification WRT A, Z, & W by wgting mislcass parameters

# Update data set up - no W in matrices
xmat <- matrix(c(ones, dat$a, dat$z, dat$z*dat$a),nrow=nrow(dat)) # design matrix
xmat1 <- matrix(c(ones, ones, dat$z, dat$z),nrow=nrow(dat)) # x set to 1 for all
xmat0 <- matrix(c(ones, 1-ones, dat$z, 1-ones),nrow=nrow(dat)) # x set to 0 for all

xmatdag <- matrix(c(ones, dat$y, dat$a, dat$z),nrow=nrow(dat)) # design matrix for misclass model
xmatdag1 <- matrix(c(ones, ones, dat$a, dat$z),nrow=nrow(dat)) # y set to 1 for all
xmatdag0 <- matrix(c(ones, 1-ones, dat$a, dat$z),nrow=nrow(dat)) # y set to 0 for all

# Matrices for the selection models
xmatddag <- matrix(c(ones, dat$a, dat$z, dat$w),nrow=nrow(dat))
xmatdia <- matrix(c(ones, dat$a, dat$z),nrow=nrow(dat))

# Use MLE to get starting values
modr <- glm(r ~ a + z + w, data=dat, family='binomial')
modr_stab <- glm(r ~ a + z, data=dat, family='binomial')
withwt <- dat %>% mutate(prr = predict(modr,.,type="response"), 
                      prr_stab = predict(modr_stab,.,type="response"),
                      wt = prr/(1-prr)*(1-prr_stab)/prr_stab)
modsesp <- geeglm(ystar ~ y + a + z, data=withwt[withwt$r==0,], family='binomial', 
                  weights=withwt$wt[withwt$r==0], id=id, corstr="independence")

# Implement M-estimation
mest(estfx=gc_dme_wtme_ef, init=c(rep(0,8),modsesp$coeff,modr$coeff,modr_stab$coeff))


