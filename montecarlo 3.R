## Vectorize the rnorm function
vrnorm <- Vectorize(rnorm,c("n","mean","sd"))

mcdat <- function() {
  dat <- mat.or.vec(7639,4)
  size <- c(9,12,14,12,24,23,47,48,46,48,47,49,45,39,8,3054,1661,983,323,441,154,213,143,85,49,35,16,8,3)
  dat[,1] <- rep(c(1,0),times=c(471,7168))
  dat[,2] <- rep(1:29,times=size)
  # Define means and standard deviations for the logit of propensity scores (z)
  meanz <- c(-5.42, -4.25, -3.42, -2.98, -2.57, -2.15, -1.75, -1.24, -0.83, -0.46, 
             -0.01, 0.53, 1.16, 1.99, 3.7, -5.38, -4.27, -3.5, -2.97, -2.55, 
             -2.14, -1.77, -1.27, -0.86, -0.47, -0.04, 0.56, 1.14, 1.8)
  sdz <- c(0.4, 0.24, 0.2, 0.1, 0.17, 0.07, 0.15, 0.14, 0.11, 0.12, 0.11, 0.22, 
           0.2, 0.33, 1.09, 0.39, 0.25, 0.2, 0.1, 0.16, 0.07, 0.15, 0.14, 0.1, 
           0.11, 0.15, 0.17, 0.17, 0.17)
  z <- vrnorm(size,meanz,sdz)
  z <- unlist(z)
  dat[,3] <- z ## Please double check the distributions of z see if the output match the parameters
  # Define means and standard deviations for reading scores (y)
  meany <- c(47.41,47.25,42.66, 45.24,44.07,36.63,37.41,39.79,37.78,35,34.3,32.06,31.09,33.87,32.84,
             59.22,55.62,52.49,52.45,50.05,48.77,45.89,46,43.28,42.18,40.45,33.03,35.27,40.47)
  sdy <- c(11.73,15.50,15.92,17.32,14.43,13.88,10.92,15.37,12.27,9.77,9.05,9.31,11.96,9.15,
           7.71,11.64,11.30,11.95,12.53,12.15,12.67,12.40,13.37,11.52,13.34,13.08,14.33,8.25,8.60)
  y <- vrnorm(size,meany,sdy)
  y <- unlist(y)
  dat[,4] <- y ## Please double check the distributions of y see if the output match the parameters
  dat <- data.frame(dat)
  colnames(dat) <- c("x","strata","z","y")
  rxz <- cor(dat$x,dat$z)
  ryz <- cor(dat$y,dat$z)
  rxy <- cor(dat$x,dat$y)
  rxyonz <- (rxy-rxz*ryz)/sqrt((1-rxz^2)*(1-ryz^2))
  result <- c(rxz,ryz,rxy,rxyonz)
  return(result)
}

mcdat()

## Run Monte Carlo Simulation

output <- t(replicate(1000,mcdat()))
colnames(output) <- c("rxz", "ryz", "rxy",'rxy|z')
colMeans(output)

## Example 1: Hong & Raudenbush (2005) Kindergarten Retention
rxc <- sort(runif(1000,-1,1)) ## 0.5 is probably a plausible boundary as the correlation between propensity score and x is just around 0.5
ryc <- sort(runif(1000,-1,1))
hong <- function(rxc, ryc){
  beta1 <- ((-0.13)*sqrt((1-0.45^2)*(1-0.52^2))-rxc*ryc)/sqrt(((1-0.45^2-ryc^2)*(1-0.52^2-rxc^2)-((-0.13)*sqrt((1-0.45^2)*(1-0.52^2))-rxc*ryc)^2)/7636)
  return(pnorm(-1.96-beta1))
}
valid <- function(rxc,ryc){
  term1 <- (1-0.45^2-ryc^2)
  term2 <- (1-0.52^2-rxc^2)
  term3 <- ((-0.13)*sqrt((1-0.45^2)*(1-0.52^2))-rxc*ryc)^2
  return(term1*term2>term3)
}

pmc_grid <- outer(rxc, ryc, hong)
valid_grid <- outer(rxc,ryc,valid)
contour(rxc,ryc, pmc_grid, main = "Contour Plot of pmc(retention)")
contour(rxc,ryc, valid_grid, main = "Contour Plot of pmc(retention)",add=TRUE,lty=2)

## Average pmc and its 95% CI
mean(pmc_grid,na.rm = TRUE)
quantile(pmc_grid,probs=c(0.025,0.975),na.rm = TRUE)

## Example 2: Father's Occupation and Educational Attainment
rcx <- sort(runif(1000,-1,1))
rcy <- sort(runif(1000,-1,1))
edu <- function(rcx, rcy) {
  beta2 <- ((0.29)*sqrt((1-0.18^2)*(1-0.22^2))-rcx*rcy)/sqrt(((1-0.18^2-rcy^2)*(1-0.22^2-rcx^2)-((0.29)*sqrt((1-0.18^2)*(1-0.22^2))-rcx*rcy)^2)/10563)
  return(pnorm(-1.96+beta2))
}
valid <- function(rcx,rcy){
  term1 <- (1-0.18^2-rcy^2)
  term2 <- (1-0.22^2-rcx^2)
  term3 <- ((0.29)*sqrt((1-0.18^2)*(1-0.22^2))-rcx*rcy)^2
  return(term1*term2>term3)
}
pmc2_grid <- outer(rcx, rcy, edu)
valid_grid <- outer(rcx,rcy,valid)
contour(rcx,rcy, pmc2_grid, main = "Contour Plot of pmc(education)")
contour(rcx,rcy, valid_grid, main = "Contour Plot of pmc(education)",add=TRUE,lty=2)
## Average pmc and its 95% CI
mean(pmc2_grid,na.rm = TRUE)
quantile(pmc2_grid,probs=c(0.025,0.975),na.rm = TRUE)






