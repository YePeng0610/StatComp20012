#' @title Descriptive statistical function using R
#' @description Output the mean and variance of a set of data, or the median absolute deviation and median
#' @param x A set of data
#' @param parametirc Decide what parameters to output
#' @param print Print the results or not
#' @return if parametirc=TRUE then output mean and variance,else output median absolute deviation and medianche
#' @importFrom stats mad median sd
#' @importFrom Rcpp evalCpp
#' @useDynLib StatComp20012
#' @examples
#' \dontrun{
#'x<-rnorm(500)
#'y<- desstat(x,parametirc=TRUE,print=TRUE)
#' }
#' @export
desstat<- function(x,parametirc=TRUE,print=FALSE) {
  if(parametirc){
    center<-mean(x);spread<-sd(x)
  }else{
    center<- median(x);spread<-mad(x)
  }
  if(print & parametirc){
    cat("Mean=",center,"\n","SD=",spread,"\n")
  }else if(print & !parametirc){
    cat("Median=",center,"\n","MAD=",spread,"\n")
  }
  result<-list(center=center,spread=spread)
  return(result)
}




#' @title Estimation of Confidence interval of population mean
#' @description Estimate the confidence interval of mean of the normal population(The population variance is known)
#' @param x A set of data from a normal population
#' @param alpha Significance level
#' @param sigma Standard deviation of population distribution
#' @return Two critical values and an interval
#' @importFrom stats qt qnorm
#' @examples
#' \dontrun{
#' x<-rnorm(10,mean=0,sd=1)
#' m<-mu(alpha=0.05,x,sigma=1)
#' m
#' }
#' @export
mu <- function(alpha,x,sigma=NA){
  n <- length(x)
  meanx <- mean(x)
  if(is.na(sigma)){
    t1 <- qt(1-alpha/2,n-1)
    t2 <- qt(1-alpha,n-1)
    mu11 <- meanx - t1*sqrt(sum((x-meanx)^2)/(n-1))/sqrt(n)
    mu12 <- meanx + t1*sqrt(sum((x-meanx)^2)/(n-1))/sqrt(n)
    mu21 <- meanx + t2*sqrt(sum((x-meanx)^2)/(n-1))/sqrt(n)
    mu22 <- meanx - t2*sqrt(sum((x-meanx)^2)/(n-1))/sqrt(n)
  }
  else{
    u1 <- qnorm(1-alpha/2,0,1)
    u2 <- qnorm(1-alpha,0,1)
    mu11 <- meanx - u1*sigma/sqrt(n)
    mu12 <- meanx + u1*sigma/sqrt(n)
    mu21 <- meanx + u2*sigma/sqrt(n)
    mu22 <- meanx - u2*sigma/sqrt(n)
  }
  string1 <- paste('The bilateral confidence interval is:[',mu11,', ',mu12,'].',sep='')
  string2 <- paste('Upper bound of unilateral confidence interval:',mu21,'.',sep='')
  string3 <- paste('Lower bound of one side confidence interval:',mu22,'.',sep='')
  string <- data.frame(Confidence_Interval=c(string1,string2,string3))
  return(string)
}



#' @title Estimation of Confidence interval of population variance
#' @description Estimate the confidence interval of variance from the normal population(The population mean is known)
#' @param x A set of data from a normal population
#' @param alpha Significance level
#' @param mu mean of population distribution
#' @return Two critical values and an interval
#' @importFrom stats qchisq
#' @examples
#' \dontrun{
#' x<-rnorm(10,mean=0,sd=1)
#' m<-sigmA(alpha=0.05,x,mu=0)
#' m
#' }
#' @export
sigmA <- function(alpha,x,mu=NA){
  n <- length(x)
  if(is.na(mu)){
    meanx <- mean(x)
    chisq11 <- qchisq(1-alpha/2,n-1)
    chisq12 <- qchisq(alpha/2,n-1)
    chisq21 <- qchisq(alpha,n-1)
    chisq22 <- qchisq(1-alpha,n-1)
    sigma11 <- sqrt(sum((x-meanx)^2)/chisq11)
    sigma12 <- sqrt(sum((x-meanx)^2)/chisq12)
    sigma21 <- sqrt(sum((x-meanx)^2)/chisq21)
    sigma22 <- sqrt(sum((x-meanx)^2)/chisq22)
  }
  else{
    chisq11 <- qchisq(1-alpha/2,n)
    chisq12 <- qchisq(alpha/2,n)
    chisq21 <- qchisq(alpha,n)
    chisq22 <- qchisq(1-alpha,n)
    sigma11 <- sqrt(sum((x-mu)^2)/chisq11)
    sigma12 <- sqrt(sum((x-mu)^2)/chisq12)
    sigma21 <- sqrt(sum((x-mu)^2)/chisq21)
    sigma22 <- sqrt(sum((x-mu)^2)/chisq22)
  }
  string1 <- paste('The bilateral confidence interval is:[',sigma11,', ',sigma12,'].',sep='')
  string2 <- paste('Upper bound of unilateral confidence interval:',sigma21,'.',sep='')
  string3 <- paste('Lower bound of one side confidence interval:',sigma22,'.',sep='')
  string <- data.frame(Confidence_Interval=c(string1,string2,string3))
  return(string)
}





#' @title Estimation of confidence interval of the difference between the mean  of two normal populations
#' @description Estimate the confidence interval of the difference between the mean  of two normal populations(The population variance is known)
#' @param alpha Significance level
#' @param x A set of data from a normal population
#' @param y A set of data from a normal population
#' @param sigmax Standard deviation of population distribution x
#' @param sigmay Standard deviation of population distribution y
#' @return a confidence interval
#' @importFrom stats qt
#' @examples
#' \dontrun{
#' x<-rnorm(10,mean=0,sd=1)
#' y<-rnorm(10,mean=0,sd=1)
#' m<-mux_muy(alpha=0.05,x,y,sigmax=1,sigmay=1)
#' m
#' }
#' @export
mux_muy <- function(alpha,x,y,sigmax=NA,sigmay=NA){
  if(is.na(sigmax)|is.na(sigmay)){
    meanx <- mean(x)
    meany <- mean(y)
    m <- length(x)
    n <- length(y)
    sx <- sqrt(sum((x-meanx)^2)/(m-1))
    sy <- sqrt(sum((y-meany)^2)/(n-1))
    sw <- sqrt((m-1)*sx^2/(m+n-2)+(n-1)*sy^2/(m+n-2))
    mu11 <- (meanx-meany)+qt(1-alpha/2,m+n-2)*sw*sqrt(1/m+1/n)
    mu12 <- (meanx-meany)-qt(1-alpha/2,m+n-2)*sw*sqrt(1/m+1/n)
  }
  else{
    meanx <- mean(x)
    meany <- mean(y)
    m <- length(x)
    n <- length(y)
    sx <- sqrt(sum((x-meanx)^2)/m)
    sy <- sqrt(sum((y-meany)^2)/n)
    sw <- sqrt((m-1)*sx^2/(m+n-2)+(n-1)*sy^2/(m+n-2))
    mu11 <- (meanx-meany)+qt(1-alpha/2,m+n)*sw*sqrt(1/m+1/n)
    mu12 <- (meanx-meany)-qt(1-alpha/2,m+n)*sw*sqrt(1/m+1/n)
  }
  string1 <- paste('The confidence interval is:[',mu12,', ',mu11,'].',sep='')
  return(string1)
}















