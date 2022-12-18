\name{rfma}
\alias{rfma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
FMA random-variate generation
}
\description{
Generate random variates from the distribution specified by \code{myfit} which is the fitted input model via the frequentist model averaging \dQuote{FMA} approach.
}
\usage{
rfma(n, myfit)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
   number of random variates to generate
}
  \item{myfit}{
  a list object returned by \code{fmafit} containing the four components needed for random-variate generation: \code{w}, \code{MLE_list}, \code{Fset}, and \code{data}
}
}
\details{
 \code{rfma} generates random variates that have the distribution of the model-average estimator. Each time a random variate is needed, a distribution is selected with probability equal to the corresponding weight and then a random variate from the fitted distribution is generated.
}
\value{
\code{rfma} generates random variates from the distribution specified by \code{myfit}
}
\references{
W. X. Jiang and B. L. Nelson (2018), "Better Input Modeling via Model Averaging," Proceedings of the 2018 Winter Simulation Conference, IEEE Press, 1575-1586.

B.L. Nelson and A.T. Wan and G. Zou and X. Zhang and W. X. Jiang (2021), "Reducing simulation input-model risk via input model averaging," INFORMS Journal on Computing, 33(2), 672-684.
}

\seealso{
See \code{\link{fmafit}} for creation of an input model (fitted distribution) via the frequentist model averaging \dQuote{FMA} approach.
}

\examples{
data<-rlnorm(500,meanlog=0,sdlog=0.25)
Fset<-c('gamma','weibull','normal','ED')
type<-'P' #by default type<-'Q'
J<-5  #by default J<-10
myfit<-fmafit(data,Fset,J,type)
n<-100
sim_data<-rfma(n,myfit)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{input modeling}
