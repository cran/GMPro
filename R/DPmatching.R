#' @title calculate degree profile distances between two graphs and match nodes.
#'
#' @description This function constructs empirical distributions of degree
#'   profiles for each vertex and then calculate distances between each pair of
#'   vertices, one from graph \emph{A} and the other from graph \emph{B}. The
#'   default used is the 1-Wasserstein distance. This function also matches
#'   vertices in \emph{A} with vertices in \emph{B} via the distance matrix
#'   between \emph{A} and \emph{B}. The distance matrix can be null and
#'   \emph{DPmatching} will calculate it. \emph{A} and \emph{B} cannot be null
#'   when the distance matrix is null.
#'
#' @param A,B Two 0/1 adjacency matrices.
#' @param W A distance matrix between \code{A} and \code{B}, which can be null.
#'   If null, this function will calculate it. More details in
#'   \emph{DPdistance}.
#'
#' @return \item{Dist}{The distance matrix between two graphs.} \item{match}{A
#'   vector containing matching results.}
#'
#' @importFrom transport wasserstein1d
#' @importFrom igraph graph.incidence max_bipartite_match
#' @importFrom stats rbinom
#' @examples
#' set.seed(2020)
#' n = 10; q = 1/2; s = 1; p = 1
#' Parent = matrix(rbinom(n*n, 1, q), nrow = n, ncol = n)
#' Parent[lower.tri(Parent)] = t(Parent)[lower.tri(Parent)]
#' diag(Parent) <- 0
#' ### Generate graph A
#' dA = matrix(rbinom(n*n, 1, s), nrow = n, ncol=n);
#' dA[lower.tri(dA)] = t(dA)[lower.tri(dA)]
#' A1 = Parent*dA
#' tmp = rbinom(n, 1, p)
#' n.A = length(which(tmp == 1))
#' indA = sample(1:n, n.A, replace = FALSE)
#' A = A1[indA, indA]
#' ### Generate graph B
#' dB = matrix(rbinom(n*n, 1, s), nrow = n, ncol=n);
#' dB[lower.tri(dB)] = t(dB)[lower.tri(dB)]
#' B1 = Parent*dB
#' tmp = rbinom(n, 1, p)
#' n.B = length(which(tmp == 1))
#' indB = sample(1:n, n.B, replace = FALSE)
#' B = B1[indB, indB]
#' DPmatching(A, B)
#' W = DPdistance(A, B)
#' DPmatching(A, B, W)
#'
#' @export


##############################################################################################################
## DP.wasserstein is the function that matches 2 graphs via the W1 distances between nodes' degree profiles ##
##############################################################################################################

DPmatching = function(A, B, W = NULL){
  # Inputs:
  # A and B: 2 symmetric adjacency matrices of graphs that will be matched with each other.

  # Remark: graph A and graph B may not have same size.
  #         A and B have the following properties:
  #             i) symmetrix
  #            ii) diagonals = 0
  #           iii) positive entries = 1

  # Output:
  # Dist: n.A by n.B matrix, where each entry denotes the wassertein distance between the corresponding node in
  #       graph A and the node in graph B.
  # match: the matching of indices. Represented by a vector with length n.A.

  # Require Packages: transport; igraph.

  # exclude all wrong possibilities
  if(!isSymmetric(A)) stop("Error! A is not symmetric!");
  if(!isSymmetric(B)) stop("Error! B is not symmetric!");

  n.A = dim(A)[1];
  n.B = dim(B)[1];

  if(is.null(W)){
    outa = apply(A, 1, sum); # calculate out degrees
    outb = apply(B, 1, sum); # calculate out degrees

    W = matrix(0, nrow = n.A, ncol = n.B); # distance matrix
    pair.vec = matrix(0, nrow = n.A, ncol = n.B); # partnership indicator matrix

    for (i in 1: n.A){
      temp_a = outa[A[i, ] != 0]
      W[i,] = apply(B, 1, function(x) wasserstein1d(temp_a, outb[x!=0]));
      ind= which(W[i,] == min(W[i, ])) # if dis(i, j) is minimum for all possible j
      pair.vec[i, ind] = 1;
    }
  }
  else{
    pair.vec = matrix(0, nrow = n.A, ncol = n.B); # partnership indicator matrix

    for (i in 1: n.A){
      ind= which(W[i,] == min(W[i, ])) # if dis(i, j) is minimum for all possible j
      pair.vec[i, ind] = 1;
    }
  }

  gm <- graph.incidence(pair.vec)
  m <- max_bipartite_match(gm) # find the maximum bipartite matching
  match = m$matching[1: n.A] - n.A #
  return(list(Dist = W, match = match))
}
