#' @title Post-processing step for edge-exploited graph matching.
#' @description Funtions \emph{DPmatching} or \emph{DPedge} can produce a
#'   preliminary graph matching result. This function, \emph{EEPost} works on
#'   refining the result iteratively. In addition, \emph{EEpost} is able to
#'   provide a convergence indicator vector \emph{FLAG} for each matching as a
#'   reference for the certainty about the matching since in practice,it has
#'   been observed that the true matches usually reach the convergence and stay
#'   the same after a few iterations, while the false matches may keep changing
#'   in the iterations.
#' @details Similar to function \emph{EEpre}, \emph{EEpost} uses maximum
#'   bipartite matching to maximize the number of common neighbours for the
#'   matched vertices with the knowledge of a preliminary matching result by
#'   defining the similarity between \eqn{i \in A} and \eqn{j \in B} as the
#'   number of common neighbours between \eqn{i} and \eqn{j} according to the
#'   preliminary matching. Then, given a matching result \eqn{\Pi_t}, post
#'   processing step is to seek a refinement \eqn{\Pi_{t+1}} satisfying
#'   \eqn{\Pi_{t+1} \in} argmax \eqn{\langle \Pi, A \Pi_t B \rangle}, where
#'   \eqn{\Pi} is a permutation matrix of dimension \eqn{(n_A, n_B)}.
#' @param W A distance matrix.
#' @param A,B Two 0/1 adjacency matrices.
#' @param rep A positive integer, indicating the number of iterations.
#' @param d A positive integer, indicating the number of candidate matching. The
#'   default value is 1.
#' @param tau A positive threshold. The default value is \eqn{rep/10}.
#' @param matching A preliminary matching result for \emph{EEpost}. If
#'   \code{matching} is null, \emph{EEpost} will apply \emph{DPedge} accordingly
#'   to generate the initial matching.
#' @importFrom igraph graph.incidence max_bipartite_match
#' @importFrom stats rbinom na.omit
#' @importFrom transport wasserstein1d
#' @return \item{Dist}{The distance matrix between two graphs.} \item{match}{A
#'   vector containing matching results.} \item{FLAG}{An indicator vector
#'   indicating whether the matching result is converged. 0 for No and 1 for
#'   Yes.} \item{converged.match}{Converged match result. \code{NA} indicates
#'   the matching result for a certain node is not v=convergent.}
#'   \item{converged.size}{The number of converged nodes.}
#' @examples
#' set.seed(2020)
#' n = 10;p = 1; q = 1/2; s = 1
#' Parent = matrix(rbinom(n*n, 1, q), nrow = n, ncol = n)
#' Parent[lower.tri(Parent)] = t(Parent)[lower.tri(Parent)]
#' diag(Parent) <- 0
#' ### Generate graph A
#' dA = matrix(rbinom(n*n, 1, s), nrow = n, ncol=n)
#' dA[lower.tri(dA)] = t(dA)[lower.tri(dA)]
#' A1 = Parent*dA;
#' tmp = rbinom(n, 1, p)
#' n.A = length(which(tmp == 1))
#' indA = sample(1:n, n.A, replace = FALSE)
#' A = A1[indA, indA]
#' ### Generate graph B
#' dB = matrix(rbinom(n*n, 1, s), nrow = n, ncol=n)
#' dB[lower.tri(dB)] = t(dB)[lower.tri(dB)]
#' B1 = Parent*dB
#' tmp = rbinom(n, 1, p)
#' n.B = length(which(tmp == 1))
#' indB = sample(1:n, n.B, replace = FALSE)
#' B = B1[indB, indB]
#' matching1= DPmatching(A, B)$Dist
#' EEpost(A = A, B = B, rep = 10, d = 5)
#' EEpost(A = A, B = B, rep = 10, d = 5, matching = matching1)
#'
#' @export

# With input networks A and B, and a matching between A and B (not requested to be perfect match or correct match), DP_plus_seed
# returns a one to one matching between A and B, with a convergence parameter to show whether the result converges or not.

EEpost = function(W = NULL, A, B, rep, tau = NULL, d = NULL, matching = NULL){

  # Input: 1) W (optional): the distance matrix between nodes in A and nodes in B.
  #        2) A and B: 2 graphs to be matched with each other
  #        3) rep: number of iterations
  #        4) tau (optional): a postive threshold to detect convergence. Default value is rep/10.
  #        5) d (optional): number of matches to consider for one node in the intial matching.
  #        6) matching  (optional): the matching can be used as the intial input
  #
  # Output: 1) resultmatch: the matching of indices, indicating for the nodes 1:n.A in graph A, the matching node in graph B.
  #         2) convergence: a vector with length the same as the size of A. convergence[i] = 0 indicates the matching pair
  #             corresponding to node i in A does not converge. convergence[i] = 1 indicates the matching pair converges.
  #         3) converged.match: the matching of indices with convergence as 1. The entries for unconverged indices are NA.
  #         4) converged.size: the number of converged pairs
  #         5) Dist: the wassertein distance between nodes.

  n.A = dim(A)[1]
  n.B = dim(B)[1]

  # exclude all wrong possibilities
  if(!isSymmetric(A)) stop("Error! A is not symmetric!");
  if(!isSymmetric(B)) stop("Error! B is not symmetric!");
  if(rep %% 1 != 0) stop("Error! rep is not an integer!");
  if(rep <= 0) stop("Error! rep is nonpositive!");
  if(!is.null(tau)){
    if (tau <= 0) stop("Error! tau is nonpositive!")
  }

  if (is.null(d) &
      is.null(matching))
    stop("Error! Please input d or matching!")

  ##### Fit with edge exploited
  if(is.null(W)){
    outa = apply(A, 1, sum);
    outb = apply(B, 1, sum);

    W = matrix(0, nrow = n.A, ncol = n.B);

    for (i in 1: n.A){
      temp_a = outa[A[i, ] != 0];
      W[i,] = apply(B, 1, function(x) wasserstein1d(temp_a, outb[x!=0]));
    }
  }

  if (!is.null(d)) {
    Pi_new = matrix(0, nrow = n.A, ncol = n.B)

    for (i in 1:n.A) {
      tmp = W[i, ]

      t = sort(tmp, decreasing = FALSE)[d]

      if (t > 0) {
        ind = which(tmp <= t)
        Pi_new[i, ind] = 1
      }
    }
  }
  else{
    Pi_new = matrix(0, nrow = n.A, ncol = n.B)

    match1 = which(!is.na(matching), arr.ind = TRUE)

    match2 = na.omit(matching)

    for (i in 1:length(match1)) {
      Pi_new[match1[i], match2[i]] = 1
    }
  }

  Pi = A %*% Pi_new %*% B;
  m <- graph.incidence(Pi, weighted = TRUE);
  gm <- max_bipartite_match(m); # Find the new matching with maximum bipartite matching

  match_next = gm$matching[1:n.A] - n.A;

  #### For the second part, we want to find the matching of the whole graphs based on the result from the first part.
  ########################################################
  weightage = 2;
  Pi_old = Pi_new;
  conv = rep(0, rep); flag = rep(0, n.A);
  for(t in 1:rep){
    match_pre = match_next;
    Pi_new = matrix(0, nrow = n.A, ncol = n.B);
    match1 = which( !is.na(match_pre), arr.ind=TRUE);
    match2 = na.omit(match_pre);
    for (i in 1: length(match1)){
      Pi_new[match1[i], match2[i]] = 1
    }
    Pi = A %*% Pi_new %*% B;

    m <- graph.incidence(Pi, weighted = TRUE);
    gm <- max_bipartite_match(m) ;
    match_next = gm$matching[1:n.A] - n.A;

    indflag = which(match_next == match_pre); #To check whether the convergence happens or not
    flag[indflag] = flag[indflag] + 1;
    flag[-indflag] = 0;
    Pi_old = Pi_new;
  }
  if(!is.null(tau)){
    convergence = 1*(flag >= tau);
  }
  else{
    convergence = 1*(flag >= ceiling(rep/10));
  }
  convresult = match_next; convresult[convergence == 0] = NA;
  if(sum(convergence) < 1){warning("The algorithm may not converge!");}
  return(list(Dist = W, match = match_next, FLAG = convergence,
              converged.match = convresult, converged.size = sum(convergence)))
}

