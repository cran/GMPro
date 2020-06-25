#' @title Degree profile graph matching with community detection.
#' @description Given two community-structured networks, this function first
#'   applies a spectral clustering method \emph{SCORE} to detect perceivable
#'   communities and then applies \emph{DPmatching} or \emph{EEpost} to match
#'   different communities. More details are in \emph{SCORE}, \emph{DPmatching}
#'   and \emph{EEpost}.
#' @details The graphs to be matched are expected to have community structures.
#'   The result is the collection of all possible permutations on
#'   \code{{1,...,K}}.
#' @param A,B Two 0/1 addjacency matrices.
#' @param K A positive integer, the number of communities in \code{A} and \code{B}.
#' @param fun A graph matching algorithm. Choices include
#'   \emph{DPmatching} and \emph{EEpost}.
#' @param rep A parameter if choosing \emph{EEpost} as the initial graph matching algorithm.
#' @param tau Optional parameter if choosing \emph{EEpost} as the initial graph matching
#'   algorithm. The default value is \eqn{rep/10}.
#' @param d Optional parameter if choosing \emph{EEpost} as the initial graph
#'   matching algorithm. The default value is 1.
#' @return A list of matching results for all possible permutations on \code{{1,...,K}}.
#' @importFrom combinat permn
#' @importFrom stats na.omit rbinom runif
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
#' DP_SBM(A = A, B = B, K = 2, fun = "EEpost", rep = 10, d = 3)
#'
#' @export

DP_SBM=function(A,B,K,fun=c("DPmatching","EEpost"),rep=NULL,tau=NULL,d=NULL){

  # Inputs: 1) A and B: 2 symmetric adjacency matrices of graphs under SBM that will be matched with each other.
  #         2) K: the number of communities in A and B.

  # Remark: graph A and graph B should have identical size.
  #         A and B have the following properties:
  #             i) symmetrix
  #            ii) diagonals = 0
  #           iii) positive entries = 1

  # Output: 1) match: a list of matching results for all possible permutations on {1, ..., K}.
  #            each element represents one kind of matching of indices. Represented by a vector with length n.

  # Require Packages: transport; igraph; combinat.

  n = dim(A)[1];

  # exclude all wrong possibilities
  if(!isSymmetric(A)) stop("Error! A is not symmetric!");
  if(!isSymmetric(B)) stop("Error! B is not symmetric!");
  if(K %% 1 != 0) stop("Error! K is not an integer!");
  if(K > n) stop("Error! More communities than nodes!");

  # do community detection in graph A and graph B seperately
  est.SC.A = SCORE(A, K);
  est.SC.B = SCORE(B, K);

  temp = 1:K;
  perm = permn(1:K);

  # match communities in A and communities in B
  all_matching = vector("list", K);

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
      }
      all_matching[[i]] = matching;
    }
  }
  return(match = all_matching)
}

