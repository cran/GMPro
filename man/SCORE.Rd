% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SCORE.R
\name{SCORE}
\alias{SCORE}
\title{Spectral Clustering On Ratios-of-Eigenvectors.}
\usage{
SCORE(G, K, itermax = NULL, startn = NULL)
}
\arguments{
\item{G}{A 0/1 adjacency matrix.}

\item{K}{A positive integer, indictaing the number of underlying communities in graph \code{G}.}

\item{itermax}{\code{k-means} parameter, indicating the maximum number of
iterations allowed. The default value is 100.}

\item{startn}{\code{k-means} parameter. If centers is a number, how many
random sets should be chosen? The default value is 10.}
}
\value{
A label vector.
}
\description{
Using ratios-of-eigenvectors to detect underlying communities.
}
\details{
\emph{SCORE} is fully established in \emph{Fast community detection by
  SCORE} of Jin (2015). \emph{SCORE} uses the entry-wise ratios between the
  first leading eigenvector and each of the other leading eigenvectors for
  clustering.
}
\examples{
set.seed(2020)
n = 10; K = 2
P  = matrix(c(1/2, 1/4, 1/4, 1/2), byrow = TRUE, nrow = K)
distribution = c(1, 2)
l = sample(distribution, n, replace=TRUE, prob = c(1/2, 1/2))
Pi = matrix(0, n, 2)
for (i in 1:n){
  Pi[i, l[i]] = 1
  }
### define the expectation of the parent graph's adjacency matrix
Omega = Pi \%*\% P \%*\% t(Pi)
### construct the parent graph G
G = matrix(runif(n*n, 0, 1), nrow = n)
G = G - Omega
temp = G
G[which(temp >0)] = 0
G[which(temp <=0)] = 1
diag(G) = 0
G[lower.tri(G)] = t(G)[lower.tri(G)]
SCORE(G, 2)

}
\references{
Jin, J. (2015) \emph{Fast community detection by score},
  \emph{The Annals of Statistics 43 (1),
  57–89}\cr\url{https://projecteuclid.org/euclid.aos/1416322036}\cr
}
