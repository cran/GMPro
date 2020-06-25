#' @title Spectral Clustering On Ratios-of-Eigenvectors.
#' @description Using ratios-of-eigenvectors to detect underlying communities.
#' @details \emph{SCORE} is fully established in \emph{Fast community detection by
#'   SCORE} of Jin (2015). \emph{SCORE} uses the entry-wise ratios between the
#'   first leading eigenvector and each of the other leading eigenvectors for
#'   clustering.
#' @param G A 0/1 adjacency matrix.
#' @param K A positive integer, indictaing the number of underlying communities in graph \code{G}.
#' @param itermax \code{k-means} parameter, indicating the maximum number of
#'   iterations allowed. The default value is 100.
#' @param startn \code{k-means} parameter. If centers is a number, how many
#'   random sets should be chosen? The default value is 10.
#' @return A label vector.
#' @importFrom stats kmeans runif
#' @references Jin, J. (2015) \emph{Fast community detection by score},
#'   \emph{The Annals of Statistics 43 (1),
#'   57–89}\cr\url{https://projecteuclid.org/euclid.aos/1416322036}\cr
#' @examples
#' set.seed(2020)
#' n = 10; K = 2
#' P  = matrix(c(1/2, 1/4, 1/4, 1/2), byrow = TRUE, nrow = K)
#' distribution = c(1, 2)
#' l = sample(distribution, n, replace=TRUE, prob = c(1/2, 1/2))
#' Pi = matrix(0, n, 2)
#' for (i in 1:n){
#'   Pi[i, l[i]] = 1
#'   }
#' ### define the expectation of the parent graph's adjacency matrix
#' Omega = Pi %*% P %*% t(Pi)
#' ### construct the parent graph G
#' G = matrix(runif(n*n, 0, 1), nrow = n)
#' G = G - Omega
#' temp = G
#' G[which(temp >0)] = 0
#' G[which(temp <=0)] = 1
#' diag(G) = 0
#' G[lower.tri(G)] = t(G)[lower.tri(G)]
#' SCORE(G, 2)
#'
#' @export


####################################################
######## Spectral Clustering Method: SCORE #########
####################################################

# Assume there are n nodes and K communities
# Before applying SCORE, need to:
# 1) transform the network graph into an n by n adjacency matrix. It has following properties:
#    i)   symmetrix
#    ii)  diagonals = 0
#    iii) positive entries = 1.


##### SCORE #####
# spectral clustering On ratios-of-eigenvectors
SCORE = function(G, K, itermax = NULL, startn = NULL){
  # Inputs:
  # 1) G: an n by n symmetric adjacency matrix whose diagonals = 0 and positive entries = 1.
  # 2) K: a positive integer which is no larger than n. This is the predefined number of communities.

  # Optional Arguments for Kmeans:
  # 1) itermax: the maximum number of iterations allowed.
  # 2) nstart: R will try startn different random starting assignments and then select the one with the lowest within cluster variation.

  # Outputs:
  # 1) a factor indicating nodes' labels. Items sharing the same label are in the same community.

  # Remark:
  # SCORE only works on connected graphs, i.e., no isolated node is allowed.

  # exclude all wrong possibilities:
  if(!isSymmetric(G)) stop("Error! G is not symmetric!")
  if(K > dim(G)[1]) stop("Error! More communities than nodes!")
  if(K %% 1 != 0) stop("Error! K is not an integer!")
  if(K  <= 0) stop("Error! Nonpositive K!")

  g.eigen = eigen(G)
  if(sum(g.eigen$vectors[, 1]==0) > 0) stop("Error! Zeroes in the first column")
  R = g.eigen$vectors[, -1]
  R = R[, 1: (K-1)]
  R = R / g.eigen$vectors[, 1]

  # apply Kmeans to assign nodes into communities
  if(!is.null(itermax) & !is.null(startn)){
    result = kmeans(R, K, iter.max = itermax, nstart = startn) #apply kmeans on ratio matrix
  }

  if(!is.null(itermax) & is.null(startn)){
    result = kmeans(R, K, iter.max = itermax, nstart = 10) #apply kmeans on ratio matrix
  }

  if(is.null(itermax) & !is.null(startn)){
    result = kmeans(R, K, iter.max = 100, nstart = startn) #apply kmeans on ratio matrix
  }

  else{
    result = kmeans(R, K, iter.max = 100, nstart = 10) #apply kmeans on ratio matrix
  }

  est = as.factor(result$cluster)
  return(est)
}
