\name{fmafit}
\alias{fmafit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Building frequentist model averaging input models
}
\description{
Creation of an input model (fitted distribution) via the frequentist model averaging \dQuote{FMA} approach.
}
\usage{
fmafit(X,Fset,J,type)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{X}{
  a numerical vector of nonzero length containing data values for fitting
}
  \item{Fset}{
  a list of character strings that specifies the set of candidate distributions; supported distributions are \code{'normal'}, \code{'lognormal'}, \code{'exponential'}, \code{'gamma'}, \code{'weibull'}, \code{'inverse gaussian'}, \code{'student t'}, \code{'uniform'}, \code{'cauchy'}, \code{'loglogistic'}, \code{'ED'}
}
  \item{J}{
  number of groups to divide the data into for cross-validation; if not specified, \code{J = 10}
}
  \item{type}{
  a character string specifying the type of model
averaging estimator: \code{'P'} for probability, \code{'Q'} for quantile; if not specified, \code{type = 'P'}
}
}
\details{
\code{fmafit} first fits each candidate parametric distribution in
\code{Fset} to the data \code{X} using maximum likelihood estimation, which yields a set of fitted distributions F = \{F_1, F_2,..., F_q\}. The MLEs for each distribution are returned as \code{MLE_list}. Next a weight vector \code{w} = \{w_1, w_2,..., w_q\} is obtained through cross-validation and also returned. The resulting model-average estimator of the true cumulative distribution of the data is

\deqn{F(x,w) = \sum_{m=1}^{q}(w_m)(F_m(x)).}

The model average fitting can be either in the cumulative probability space or quantile space. The difference between the two types of model averaging is only in the weight vector associated with the candidate distributions, which is obtained through cross-validation in either the probability or quantile space.
}
\value{
\code{fmafit} returns an object which is a list with four components:
  \item{w}{weight vector associated with distributions in
\code{Fset}}
  \item{MLE_list}{list of MLE parameter estimates for each candidate distribution with \code{'NA'} for \code{ED} (the empirical distribution)}
  \item{Fset}{same as the input argument}
  \item{data}{sorted input argument \code{X} (needed for ED)}
}
\references{
W. X. Jiang and B. L. Nelson (2018), "Better Input Modeling via Model Averaging," Proceedings of the 2018 Winter Simulation Conference, IEEE Press, 1575-1586.
}

\seealso{
See \code{\link{rfma}} for random-variate generation from the fitted distribution obtained via \dQuote{FMA} approach.
}

\examples{
data<-rlnorm(500,meanlog=0,sdlog=0.25)
Fset<-c('gamma','weibull','normal','ED')
type<-'P' #by default type<-'Q'
J<-5  #by default J<-10
myfit<-fmafit(data,Fset,J,type)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
