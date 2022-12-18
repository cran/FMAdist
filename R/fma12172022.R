# R code for input modeling via frequentist model average
# Xi Jiang and Barry L Nelson
# Last update 12/17/2022

# The following distributions are supported:  'normal', 'lognormal', 'beta', 'exponential',
# 'gamma', 'weibull', 'inverse gaussian', 'logistic', 'loglogistic', 'student t', 'uniform',
# 'cauchy', 'pareto', 'rayleigh', 'ED'.

## Example:
## data<-rlnorm(500,meanlog=0,sdlog=0.25)
## Fset<-c('gamma','weibull','normal','ED')
## type<-'P' #by default type<-'Q'
## J<-5  #by default J<-10
## myfit<-fmafit(data,Fset,J,type)
## n<-10000
## sim_data<-rfma(n,myfit)

#' @import stats
#' @import utils

package_install<-function(packages){
  # make sure not to install if the packages is already there
  # packages = the list of packages need to eb installed to have this package work
  # install the required packages only if it's not installed
  for (i in 1:length(packages)){
    package.need<-packages[i]
    if(length(find.package(package.need,quiet=TRUE))==0){
      install.packages(packages[i])
    }
  }
}

# The following packages will be installed if needed
packages<-c("fitdistrplus","actuar","EnvStats","extraDistr","MASS","quadprog")
package_install(packages)

MLE<-function(dist,data){
  # Estimate the MLE parameters of a distribution using data
  # dist = candidate distribution
  # data = vector of data for estimating MLE
  # return the MLE parameters, max value of loglikelihood, AIC and BIC
  if (dist=='normal'){
    fit.norm<-fitdistrplus::fitdist(data,'norm', method='mle')
    theta=fit.norm$estimate #c(mean,sd)
    ll=fit.norm$loglik
    AIC=fit.norm$aic
    BIC=fit.norm$bic
  } else if (dist=='lognormal'){
    fit.lnorm<-fitdistrplus::fitdist(data,'lnorm',method='mle')
    theta=fit.lnorm$estimate #c(meanlog,sdlog)
    ll=fit.lnorm$loglik
    AIC=fit.lnorm$aic
    BIC=fit.lnorm$bic
  } else if (dist=='beta'){
    eps<-.Machine$double.eps #the smallest positive floating-point number x such that 1 + x != 1
    data<-(data-min(data))/(max(data)-min(data)) # scale data to [0,1]
    data<-eps+(1-2*eps)*(data) # scale again to (0,1)
    fit.beta<-fitdistrplus::fitdist(data,'beta',method='mle')
    theta=fit.beta$estimate #c(shape1,shape2)
    ll=fit.beta$loglik
    AIC=2*length(theta)-2*ll
    BIC=log(length(data))*length(theta)-2*ll
  } else if (dist=='exponential'){
    fit.exp <- MASS::fitdistr(data,'exponential')
    theta=fit.exp$estimate #rate
    ll=fit.exp$loglik
    AIC=2*length(theta)-2*ll
    BIC=log(length(data))*length(theta)-2*ll
  } else if (dist=='gamma'){
    fit.gamma<-fitdistrplus::fitdist(data,'gamma',method='mle')
    theta=fit.gamma$estimate #c(shape,rate)
    ll=fit.gamma$loglik
    AIC=fit.gamma$aic
    BIC=fit.gamma$bic
  } else if (dist=='weibull'){
    fit.weibull <- fitdistrplus::fitdist(data,'weibull',method='mle')
    theta=fit.weibull$estimate #c(shape,scale)
    ll=fit.weibull$loglik
    AIC=fit.weibull$aic
    BIC=fit.weibull$bic
  } else if (dist=='inverse gaussian'){
    MLE_invgauss<-fitdistrplus::fitdist(data, "invgauss", start = list(mean = mean(data), shape = length(data)/sum(1/data-1/mean(data))))
    theta=MLE_invgauss$estimate #c(mean,shape)
    ll=MLE_invgauss$logLik
    AIC=2*length(theta)-2*ll
    BIC=log(length(data))*length(theta)-2*ll
  } else if (dist=='logistic'){
    fit.logis <- MASS::fitdistr(data,'logistic')
    theta=fit.logis$estimate #c(location,scale)
    ll=fit.logis$loglik
    AIC=2*length(theta)-2*ll
    BIC=log(length(data))*length(theta)-2*ll
  } else if (dist=='loglogistic'){
    MLE_loglogis<-fitdistrplus::fitdist(data, "llogis")
    theta=MLE_loglogis$estimate #c(shape,scale)
    ll=MLE_loglogis$logLik
    AIC=2*length(theta)-2*ll
    BIC=log(length(data))*length(theta)-2*ll
  } else if (dist=='student t'){
    fit.t<-fitdistrplus::fitdist(data, "t", method = "mle", start = list(df=1))
    theta=fit.t$estimate #df
    ll=fit.t$loglik
    AIC=fit.t$aic
    BIT=fit.t$bic
  } else if (dist=='uniform'){
    fit.uniform<-fitdistrplus::fitdist(data,'unif',method='mle')
    theta=fit.uniform$estimate #c(min,max)
    ll=fit.uniform$loglik
    AIC=fit.uniform$aic
    BIC=fit.uniform$bic
  } else if (dist=='cauchy'){
    fit.cauchy <- MASS::fitdistr(data,'cauchy')
    theta=fit.cauchy$estimate #c(location,scale)
    ll=fit.cauchy$loglik
    AIC=2*length(theta)-2*ll
    BIC=log(length(data))*length(theta)-2*ll
  } else if (dist=='pareto'){
    fit.pareto <- EnvStats::epareto(data,method="mle")
    theta=fit.pareto$parameters #c(location,shape)
    ll=sum(log(EnvStats::dpareto(data, location = theta[1], shape = theta[2])))
    AIC=2*length(theta)-2*ll
    BIC=log(length(data))*length(theta)-2*ll
  } else if (dist=='rayleigh'){
    theta=sqrt(sum(data^2)/2/length(data)) #scale
    ll=sum(log(extraDistr::drayleigh(data, sigma=theta)))
    AIC=2*length(theta)-2*ll
    BIC=log(length(data))*length(theta)-2*ll
  } else if (dist=='ED'){
    theta='NA'
    ll='NA'
    AIC='NA'
    BIC='NA'
  } else{ stop("MLE: Not a supported distribution")}

  list(MLE_theta=theta, ll=ll, AIC=AIC, BIC=BIC)
}

CDF<-function(dist,x,theta){
  # CDF of a distribution
  # dist = some distribution
  # x = a value within the support of the distribution
  # theta = the parameters of the distribution
  # returns F_X(x) = P(X<=x)
  switch(dist,
         'normal'=pnorm(x,mean=theta[1],sd=theta[2]),
         'lognormal'=plnorm(x,meanlog=theta[1],sdlog=theta[2]),
         'beta'=pbeta(x,shape1=theta[1],shape2=theta[2]),
         'exponential'=pexp(x,rate=theta),
         'gamma'=pgamma(x,shape=theta[1],rate=theta[2]),
         'weibull'=pweibull(x,shape=theta[1],scale=theta[2]),
         'inverse gaussian'=actuar::pinvgauss(x,mean=theta[1],shape=theta[2]),
         'student t'=pt(x,df=theta),
         'uniform'=punif(x,min=theta[1],max=theta[2]),
         'cauchy'=pcauchy(x,location=theta[1],scale=theta[2]),
         'pareto'=EnvStats::ppareto(x, location = theta[1], shape = theta[2]),
         'rayleigh'=extraDistr::prayleigh(x, sigma = theta),
         'logistic'=plogis(x, location = theta[1], scale = theta[2]),
         'loglogistic'=actuar::pllogis(x, shape = theta[1], scale = theta[2]),
         stop("CDF: Not a supported distribution")
  )
}

IT<-function(dist,U,theta){
  # Given a vector of unif[0,1]s, generate random variate from distribution using
  #   inverse transform method
  # dist = some distribution
  # U = the list of RV~unif[0,1]
  # theta = the parameters of the distribution
  # return X = F^{-1}(Ugrid), random variates with the specified distribution
  switch(dist,
         'normal'=qnorm(U,mean=theta[1],sd=theta[2]),
         'lognormal'=qlnorm(U,meanlog=theta[1],sdlog=theta[2]),
         'beta'=qbeta(U,shape1=theta[1],shape2=theta[2]),
         'exponential'=qexp(U,rate=theta),
         'gamma'=qgamma(U,shape=theta[1],rate=theta[2]),
         'weibull'=qweibull(U,shape=theta[1],scale=theta[2]),
         'inverse gaussian'=actuar::qinvgauss(U,mean=theta[1],shape=theta[2]),
         'student t'=qt(U,df=theta),
         'uniform'=qunif(U,min=theta[1],max=theta[2]),
         'cauchy'=qcauchy(U,location=theta[1],scale=theta[2]),
         'pareto'=EnvStats::qpareto(U, location = theta[1], shape = theta[2]),
         'rayleigh'=extraDistr::qrayleigh(U, sigma = theta),
         'logistic'=qlogis(U, location = theta[1], scale = theta[2]),
         'loglogistic'=actuar::qllogis(U, shape = theta[1], scale = theta[2]),
         stop("IT: Not a supported distribution")
  )
}

lappend <- function (lst, ...){
  # a function that appends elements to a list
  # lst = orginial list
  # returns the appended list
  lst <- c(lst, list(...))
  return(lst)
}

D_matrix<-function(data,J,n,Fset){
  # Compute the D_matrix for QP that returns the optimal weight vector
  #   for probability average fitting
  # data = vector of data
  # J = number of groups for cross validation
  # n = size of the vector of data
  # Fset = set of candidate distributions
  # returns matrix D for quadratic term in the objective function of QP
  xsim_matrix <- matrix(data, nrow=J, byrow=T)
  D_matrix<-matrix(0,length(Fset),length(Fset))
  for (rep1 in 1:J){
    dat1<-as.vector(t(xsim_matrix[-rep1,]))
    dat2<-xsim_matrix[rep1,]
    ED_CV2<-ecdf(dat2)
    MLE_theta<-list()
    for (rep2 in 1:length(Fset)){
      if (Fset[rep2]!='ED'){
        MLE_theta<-lappend(MLE_theta,MLE(Fset[rep2],dat1)$MLE_theta)
      } else if (Fset[rep2]=='ED'){
        ED_CV1<-ecdf(dat1)
      }
    }
    for (rep3 in 1:(n/J)){
      c_vector<-vector('numeric')
      for (rep4 in 1:length(Fset)){
        if (Fset[rep4]!='ED' & Fset[rep4]!='beta'){
          c_vector<-append(c_vector, CDF(Fset[rep4],dat2[rep3],unlist(MLE_theta[rep4]))-ED_CV2(dat2[rep3]))
        } else if (Fset[rep4]=='beta'){#normalize the data for beta distribution
          eps<-.Machine$double.eps
          data_pt<-(dat2[rep3]-min(data))/(max(data)-min(data))
          data_pt<-eps+(1-2*eps)*(data_pt)
          c_vector<-append(c_vector, CDF(Fset[rep4],data_pt,unlist(MLE_theta[rep4]))-ED_CV2(dat2[rep3]))
        } else if (Fset[rep4]=='ED'){
          c_vector<-append(c_vector, ED_CV1(dat2[rep3])-ED_CV2(dat2[rep3]))
        }
      }
      D_matrix<-D_matrix+c_vector%*%t(c_vector)
    }
  }
  return(D_mat=D_matrix)
}

DG_matrix<-function(data,J,n,Fset){
  # Compute the D_matrix for QP which returns the optimal weight vector
  # for quantile average fitting
  # data = vector of data
  # J = number of groups for cross validation
  # n = size of the vector of data
  # Fset = set of candidate distributions
  # returns matrix D for quadratic term in the objective function of QP
  xsim_matrix <- matrix(data, nrow=J, byrow=T)
  D_matrix<-matrix(0,length(Fset),length(Fset))
  for (rep1 in 1:J){
    dat1<-as.vector(t(xsim_matrix[-rep1,]))
    dat2<-xsim_matrix[rep1,]
    MLE_theta<-list()
    for (rep2 in 1:length(Fset)){
      if (Fset[rep2]!='ED'){
        MLE_theta<-lappend(MLE_theta,MLE(Fset[rep2],dat1)$MLE_theta)
      } else if (Fset[rep2]=='ED'){
        dat1<-sort(dat1)
      }
    }
    for (rep3 in 1:(n/J)){
      c_vector<-vector('numeric')
      for (rep4 in 1:length(Fset)){
        if (Fset[rep4]!='ED' & Fset[rep4]!='beta'){
          c_vector<-append(c_vector, IT(Fset[rep4],rep3/(n/J+1),unlist(MLE_theta[rep4]))-(sort(dat2))[rep3])
        } else if (Fset[rep4]=='beta'){#denormalize the data for beta distribution
          eps<-.Machine$double.eps
          data_pt<-(IT(Fset[rep4],rep3/(n/J+1),unlist(MLE_theta[rep4]))-eps)/(1-2*eps)*(max(data)-min(data))+min(data)
          c_vector<-append(c_vector, data_pt-(sort(dat2))[rep3])
        } else if (Fset[rep4]=='ED'){
          c_vector<-append(c_vector, QT(dat1,rep3/(n/J+1))-sort(dat2)[rep3])
        }
      }
      D_matrix<-D_matrix+c_vector%*%t(c_vector)
    }
  }
  return(D_mat=D_matrix)
}

Qua_opt<-function(D_mat,Fset){
  # QP that finds the optimal weight vector
  # D_mat = the matrix D for QP in the objective function with the quadratic term
  # Fset = set of candidate distributions
  # returns the weight vector that minimizes the J-fold CV criterion
  dvec <- rep(0,nrow(D_mat))
  A.Equality <- matrix(rep(1,nrow(D_mat)), ncol=1)
  Amat <- cbind(A.Equality, diag(nrow(D_mat)))
  bvec <- c(1, rep(0, nrow(D_mat)))
  qp <- quadprog::solve.QP(D_mat, dvec, Amat, bvec, meq=1)$solution
  qp[qp<=.Machine$double.eps]=0
  qp<-qp/sum(qp)
  return(qp)
}

QT<-function(X,probs,type='non LI'){
  # quantile function of empirical CDF with no linear interpolation for sorted data
  # x = sorted samples of data
  # probs = certain probability
  # type = 'non LI' or 'LI'
  # returns the quantile of a given probability of given samples of data
  n <- length(X)
  if(type == 'LI') {
    # If desired, place linear interpolated quantile function here
  } else if (type=='non LI'){
    nppm <- n * probs
    lo <- floor(nppm)
    if (nppm==lo){
      index<-lo
    } else {
      index<-lo+1
    }
    result <- X[index]
  }
  return(result)
}

EmpCDF<-function(X,data){
  # THIS FUNCTION NOT CURRENTLY USED
  # empirical CDF of sorted data
  # X = sorted samples of data
  # data = one of the data point from X
  # returns the empirical cdf value
  n<-length(X)
  y_step<-1:(n+1)
  f.U<-stepfun(X,y_step,f=0)
  result<-(f.U(data)-1)/n
  return(result)
}

# Complete set of supported distributions
SET = c('normal', 'lognormal', 'beta', 'exponential', 'gamma', 'weibull',
        'inverse gaussian', 'logistic', 'loglogistic', 'student t',
        'uniform', 'cauchy', 'pareto', 'rayleigh', 'ED')

Inputcheck<-function(X,Fset,J,type){
  # check the input
  # X = samples of data for fitting
  # Fset = the list of candidate distributions
  # J = the number of folds for cross-validation
  # type = 'P' (probability) or 'Q' (quantile) model averaging
  # stop the function when inputs are not supported
  if (!is.numeric(X)){
    stop("X: data is not numeric")
  }
  if (!all(Fset %in% SET)){
    stop("Fset: Fset includes distributions that are not supported")
  }
  if ((J-floor(J))>=.Machine$double.eps){
    stop("J: not an integer")
  } else {
    if (J<2) {
      stop("J: J >= 2 required ")
    } else {
      if (length(X)<2*J){
      stop("X: length of X >= 2*J required")
      }
    }
  }
  if (type!='P' & type!='Q'){
    stop("type: not a valid fitting type")
  }
}

#' @export
fmafit<-function(X,Fset,J=10,type='P'){
  # Fit a model average distribution to data
  # X = samples of data for fitting
  # Fset = the list of candidate distributions
  # J = the number of folds for cross-validation
  # type = P (probability) or Q (quantile) model averaging
  # returns weight w, MLE_list, Fset, and data
  Inputcheck(X,Fset,J,type)
  n=length(X)
  set.seed(1)
  data = sample(X) # scramble in case sorted
  n <- floor(n/J)*J
  data = data[1:n]
  MLE_theta<-list()
  for (i in 1:length(Fset)){
    MLE_theta<-lappend(MLE_theta,MLE(Fset[i],data)$MLE_theta)
  }
  if (type=='P'){
    D_mat<-D_matrix(data,J,n,Fset)
    weight <- Qua_opt(D_mat,Fset)
  } else if (type=='Q'){
    DG_mat<-DG_matrix(data,J,n,Fset)
    weight <- Qua_opt(DG_mat,Fset)
  } else {
    stop("fmafit: Unsupported type")
  }
  list(w=weight,MLE_list=MLE_theta,Fset=Fset,data=sort(data))
}

#' @export
rfma<-function(n,myfit){
  # Generate n random samples from model average distribution
  # n = the number of samples to generate
  # Fset = the list of distributions
  # MLE_list = the list of MLE of parameters of Fset
  # w = the weight vector associated with distributions in Fset for MAE
  # data = fitting data (needed for ED)
  # returns vector of samples
  w<-myfit$w
  Fset<-myfit$Fset
  MLE_list<-myfit$MLE_list
  data<-myfit$data
  xsim<-vector('numeric')
  len<-length(Fset)
  for (i in 1:n){
    U<-runif(1)
    if (len==1){
      U_1<-runif(1)
      if (Fset!='ED'){
        x = IT(Fset,U_1,unname(MLE_list))
      } else if (Fset=='ED'){
        x<-unname(QT(data,U_1))
      }
    } else {
      x_step<-cumsum(w)
      x_step[len]<-1
      y_step<-1:(len+1)
      f.U<-stepfun(x_step,y_step,right=FALSE)
      k<-f.U(U)
      U_1<-runif(1)
      if (Fset[k]!= 'ED'){
        x = IT(Fset[k],U_1,unlist(MLE_list[k]))
      } else if (Fset[k]=='ED'){
        x<-unname(QT(data,U_1))
      }
    }
    xsim<-append(xsim,x)
  }
  return(xsim)
}
