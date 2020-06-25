#' @title The edge-exploited version of DPmatching.
#' @description This functions is based on \emph{DPmatching}. Instead of
#'   allowing each vertex in \emph{A} to connect to one and only one vertex in
#'   \emph{B}, here by introducing parameter \code{d}, this function allows for
#'   \code{d} edges for each vertex in \emph{A}. More details are in
#'   \emph{DPmatching}.
#' @param A,B Two symmetric 0/1 addjacency matrices.
#' @param d A positive integer, indicating the number of candidate matching.
#' @param W A distance matrix between \code{A} and \code{B}. This argumnet can
#'   be null. If \code{W} is null, \code{A} and \code{B} cannot be null.
#' @return \item{Dist}{The distance matrix between two graphs.} \item{Z}{An
#'   indicator matrix. Entry \eqn{Z_{i, j} = 1} indicates a matching between
#'   node \code{i} in graph \code{A} and node \code{j} in graph \code{B}, 0
#'   otherwise.}
#' @importFrom transport  wasserstein1d
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
#' DPedge(A, B, d = 5)
#' @export

# With input networks A and B, and an optional distance matrix W, the number of edges for each nodes d
# returns a (0, 1) matrix Z, whose (i, j) entry is 1 if the edge between node i in A and node j in B is kept.

DPedge = function(A = NULL, B = NULL, d, W = NULL){

  # Input: 1) A and B (optional): 2 graphs to be matched with each other.
  #        2) W (optional): the distance matrix between nodes in A and nodes in B.
  #        3) d: number of edges selected for each node.

  # Remark: graph A and graph B may not have same size.
  #         A and B have the following properties:
  #           i) symmetrix
  #          ii) diagonals = 0
  #         iii) positive entries = 1

  # Output: 1) Dist: the distance matrix between nodes in A and nodes in B.
  #         2) Z: rows representing nodes in A and columns representing nodes in B.
  #               Entries indicate edges.
  #               (i, j) is 1 if there is an edge between node i in A and node j in B.

  if(d %% 1 != 0) stop("Error! d is not an integer!")
  if(d <= 0) stop("Error! Nonpositive d!");

  if(is.null(W)){

    if(is.null(A)|is.null(B)) stop("Error! No values for A, B and W!");
    if(!isSymmetric(A)) stop("Error! A is not symmetric!");
    if(!isSymmetric(B)) stop("Error! B is not symmetric!");

    n.A = dim(A)[1];
    n.B = dim(B)[1];

    outa = apply(A, 1, sum);
    outb = apply(B, 1, sum);

    W = matrix(0, nrow = n.A, ncol = n.B);

    for (i in 1: n.A){
      temp_a = outa[A[i, ] != 0];
      W[i,] = apply(B, 1, function(x) wasserstein1d(temp_a, outb[x!=0]));
    }
  }

  n.A = dim(W)[1]; n.B = dim(W)[2];
  Z = matrix(0, nrow = n.A, ncol = n.B);

  for (i in 1: n.A){
    tmp = sort(W[i,], decreasing = F)
    ind = which(W[i,] <= tmp[d])
    Z[i, ind] = 1
  }
  return(list(Dist = W, Z = Z))
}
