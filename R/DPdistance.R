#' calculate degree profile distances between 2 graphs.
#'
#' This function constructs empirical distributions of degree profiles for each
#' vertex and then calculate distances between each pair of vertices, one from
#' graph \emph{A} and the other from graph \emph{B}. The default distance used is the
#' 1-Wasserstein distance.
#'
#' @param A,B Two 0/1 adjacency matrices.
#' @param fun Optional function that computes distance between two distributions.
#' @return A distance matrix. Rows represent nodes in graph \emph{A} and columns
#'   represent nodes in graph \emph{B}. Its \emph{(i, j)} element is the
#'   distance between \eqn{i \in A} and \eqn{i \in B}.
#' @importFrom transport wasserstein1d
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
#' DPdistance(A, B)
#' @export

##############################################################################################################
## DP.wasserstein is the function that matches 2 graphs via the W1 distances between nodes' degree profiles ##
##############################################################################################################

DPdistance = function(A, B, fun = NULL){
  # Inputs:
  # A and B: 2 symmetric adjacency matrices of graphs that will be matched with each other.

  # Remark: graph A and graph B may not have same size.
  #         A and B have the following properties:
  #           i) symmetrix
  #          ii) diagonals = 0
  #         iii) positive entries = 1

  # Output:
  # W: n.A by n.B matrix, where each entry denotes the wassertein distance between the corresponding node in
  #       graph A and the node in graph B.

  # Require Packages: transport.

  # exclude all wrong possibilities
  if(!isSymmetric(A)) stop("Error! A is not symmetric!");
  if(!isSymmetric(B)) stop("Error! B is not symmetric!");

  n.A = dim(A)[1];
  n.B = dim(B)[1];

  outa = apply(A, 1, sum);
  outb = apply(B, 1, sum);

  W = matrix(0, nrow = n.A, ncol = n.B);

  for (i in 1: n.A){
    temp_a = outa[A[i, ] != 0]
    if(is.null(fun)){
      W[i,] = apply(B, 1, function(x) wasserstein1d(temp_a, outb[x!=0]));
    }
    else{
      W[i,] = apply(B, 1, function(x) fun(temp_a, outb[x!=0]));
    }
  }
  return(W)
}
