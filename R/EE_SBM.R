#' @title Edge exploited degree profile graph matching with community detection.
#' @description Given two community-structured networks, this function first
#'   applies a spectral clustering method \emph{SCORE} to detect perceivable
#'   communities and then applies a certain graph matching method to match
#'   different communities.
#' @details \emph{EE_SBM} can be regarded as a post processing version of
#'   \emph{DP_SBM} using \emph{EEpost}.
#' @param A,B Two 0/1 addjacency matrices.
#' @param K A positive integer, indicating the number of communities in \code{A}
#'   and \code{B}.
#' @param fun A graph matching algorithm. Choices include \emph{DPmatching} and
#'   \code{EEpost}.
#' @param rep Optional parameter if \emph{EEpost} is the initial graph matching algorithm.
#' @param tau Optional parameter if \emph{EEpost} is the initial graph matching
#'   algorithm. The default value is \eqn{rep/10}.
#' @param d Optional parameter if \emph{EEpost} is the initial graph matching
#'   algorithm. The default value is 1.
#' @return \item{match}{A vector containing matching results.} \item{FLAG}{An
#'   indicator vector indicating whether the matching result is converged, 0 for
#'   No and 1 for Yes.}
#' @importFrom combinat permn
#' @importFrom stats na.omit runif rbinom
#' @examples
#' ### Here we use graphs under stochastic block model(SBM).
#' set.seed(2020)
#' K = 2; n = 30; s = 1;
#' P  = matrix(c(1/2, 1/4, 1/4, 1/2), byrow = TRUE, nrow = K)
#' ### define community label matrix Pi
#' distribution = c(1, 2);
#' l = sample(distribution, n, replace=TRUE, prob = c(1/2, 1/2))
#' Pi = matrix(0, n, 2) # label matrix
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
#' G[lower.tri(G)] = t(G)[lower.tri(G)];
#' ### Sample Graphs Generation
#' ### generate graph A from G
#' dA = matrix(rbinom(n*n, 1, s), nrow = n, ncol=n)
#' dA[lower.tri(dA)] = t(dA)[lower.tri(dA)]
#' A1 = G*dA
#' indA = sample(1:n, n, replace = FALSE)
#' labelA = l[indA]
#' A = A1[indA, indA]
#' ### similarly, generate graph B from G
#' dB = matrix(rbinom(n*n, 1, s), nrow = n, ncol=n)
#' dB[lower.tri(dB)] = t(dB)[lower.tri(dB)]
#' B1 = G*dB
#' indB = sample(1:n, n, replace = FALSE)
#' labelB = l[indB]
#' B = B1[indB, indB]
#' EE_SBM(A = A, B = B, K = 2, fun = "EEpost", rep = 10, d = 3)

#' @export EE_SBM

EE_SBM = function(A,B,K,fun = c("DPmatching", "EEpost"),rep = NULL,tau = NULL,d = NULL){

  # Input: 1) A and B: 2 graphs to be matched with each other.
  #        2) K: the number of communities in A and B.
  #        3) rep: number of iterations.
  #        4) d: number of matches to consider for one node.

  # Output: 1) match: the matching of indices.
  #         2) FLAG: a convergence indicator vector.

  # Require Packages: transport; igraph; combinat.

  n = dim(A)[1]

  # exclude all wrong possibilities
  if(!isSymmetric(A)) stop("Error! A is not symmetric!");
  if(!isSymmetric(B)) stop("Error! B is not symmetric!");
  if(dim(A)[1] != dim(B)[1]) stop("Error! A and B have different sizes");
  if(K %% 1 != 0) stop("Error! K is not an integer!");
  if(K > n) stop("Error! More communities than nodes!");

  # match communities in A and communities in B

  # do community detection in graph A and graph B seperately
  est.SC.A = SCORE(A, K);
  est.SC.B = SCORE(B, K);

  temp = 1:K;
  perm = permn(1:K)

  all_matching = vector("list", K);
  matched = rep(0, length(perm));
  conv = rep(0, length(perm));

  if (fun == "DPmatching"){
    for (i in 1: length(perm)){
      matching = rep(NA, n);
      for (k in 1: K){
        A_com = A[which(est.SC.A == k), which(est.SC.A == k)];
        n.A = dim(A_com)[1];
        temp1 = which(est.SC.A == k);
        temp2 = which(est.SC.B == unlist(perm[i])[k]);
        B_com = B[which(est.SC.B == unlist(perm[i])[k]), which(est.SC.B == unlist(perm[i])[k])];
        n.B = dim(B_com)[1];
        result = DPmatching(A_com, B_com);
        matching[temp1[which(!is.na(result$match), arr.ind = TRUE)]] = temp2[na.omit(result$match)];
        matched[i] = matched[i] + length(which(!is.na(matched)));
      }
      all_matching[[i]] = matching;
    }
  }

  if (fun == "EEpost"){
    if (is.null(rep)) stop("Error! Please input rep!");
    if (rep %% 1 != 0) stop("Error! rep is not an integer!");
    if (rep <= 0) stop("Error! rep is nonpositive!");
    if (is.null(d)) stop("Error! Please input d!");
    if (d %% 1 != 0) stop("Error! d is not an integer!");
    if (d <= 0) stop("Error! d is nonpositive!");
    if (is.null(tau)) tau = ceiling(rep/10);
    if(!is.null(tau)){
      if (tau <= 0) stop("Error! tau is nonpositive!")
    }
    for (i in 1: length(perm)){
      matching = rep(NA, n);
      for (k in 1: K){
        A_com = A[which(est.SC.A == k), which(est.SC.A == k)];
        n.A = dim(A_com)[1];
        temp1 = which(est.SC.A == k);
        temp2 = which(est.SC.B == unlist(perm[i])[k]);
        B_com = B[which(est.SC.B == unlist(perm[i])[k]), which(est.SC.B == unlist(perm[i])[k])];
        n.B = dim(B_com)[1];
        result = EEpost(A = A_com, B = B_com, rep = rep, tau = tau, d = d);
        matching[temp1[which(!is.na(result$match), arr.ind = TRUE)]] = temp2[na.omit(result$match)];
        conv[i] = conv[i] + result$converged.size;
      }
      all_matching[[i]] = matching;
    }
  }

  if (fun == "DPmatching") matching_old = all_matching[[which.max(matched)]];
  if (fun == "EEpost") matching_old = all_matching[[which.max(conv)]];
  result = EEpost(A = A, B = B, rep = rep, matching = matching_old);
  return(list(match = result$match, FLAG = result$FLAG))
}
