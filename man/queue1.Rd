% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/queue1.R
\name{queue1}
\alias{queue1}
\title{calculate the main quantitative metric describing the queuing system in the queuing service system.}
\usage{
queue1(lambda, mu, T)
}
\arguments{
\item{lambda}{The successive arrival intervals of customers to the system are independent and subject to an exponential distribution with a parameter of lambda. (numeric)}

\item{mu}{The service time of the service desk is also iid., and obeys the exponential distribution of parameter mu (numeric)}

\item{T}{System service time (numeric)}
}
\value{
the length Ls, the average stay time Ws and the probability of customer waiting Pwait.
}
\description{
The prediction model is described in http://www.babelgraph.org/wp/?p=358.
}
\examples{
\dontrun{
res <- queue1(4,6,1000)
}
}
