#' @title Edge exploited degree profile graph matching with preprocessing.
#' @description This function uses seeds to compute edge-exploited matching
#'   results. Seeds are nodes with high degrees. \emph{EEpre} uses seeds to
#'   extend the matching of seeds to the matching of all nodes.
#' @details The high degree vertices have many neighbours and enjoy ample
#'   information for a successful matching. Thereforem, this function employ
#'   these high degree vertices to match other nodes. If the information of
#'   seeds is unavailable, \emph{EEpre} will conduct a grid search grid search
#'   to find the optimal collection of seeds. These vertices are expected to
#'   have high degress and their distances are supposed to be the smallest among
#'   the pairs in consideration.
#' @param A,B Two 0/1 addjacency matrices.
#' @param d A positive integer, indicating the number of candicate matching.
#' @param seed A matrix indicating pair of seeds. \code{seed} can be null.
#' @param AB_dist A nonnegative distance matrix, which can be null. If
#'   \code{AB_dist} is null, \emph{EEpre} will apply \emph{DPdistance} to find
#'   it.
#' @return \item{Dist}{The distance matrix between two graphs} \item{Z}{An
#'   indicator matrix. Entry \eqn{Z_{i, j} = 1} indicates a matching between
#'   node \eqn{i \in A} and node \eqn{j \in B}, 0
#'   otherwise.}
#' @importFrom transport wasserstein1d
#' @importFrom igraph graph.incidence max_bipartite_match as_edgelist
#' @importFrom stats quantile rbinom
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
#' EEpre(A = A, B = B, d = 5)
#' @export

EEpre = function(A, B, d, seed = NULL, AB_dist = NULL){

  # Input: 1) A and B: 2 graphs to be matched with each other.
  #        2) d: number of edges selected for each node.
  #        3) seed (optional): a matrix indicating pair of seeds. Rows represent pairs of seeds.
  #        4) AB_dist (optional): the distance matrix between nodes in A and nodes in B.

  # Remark: graph A and graph B may not have the same size.
  #         A and B should have the following properties:
  #           i) symmetrix
  #          ii) diagonals = 0
  #         iii) positive entries = 1

  # Output: 1) Dist: the distance matrix between nodes in A and nodes in B.
  #         2) Z: rows representing nodes in A and columns representing nodes in B.
  #               Entries indicate edges.
  #               (i, j) is 1 if there is an edge between node i in A and node j in B.

  n.A = dim(A)[1];
  n.B = dim(B)[1];

  outa = apply(A, 1, sum);
  outb = apply(B, 1, sum);

  # exclude all wrong possibilities
  if(!isSymmetric(A)) stop("Error! A is not symmetric!");
  if(!isSymmetric(B)) stop("Error! B is not symmetric!");
  if(d %% 1 != 0) stop("Error! d is not an integer!")
  if(d <= 0) stop("Error! Nonpositive d!")
  if(!is.null(seed)){
    if(dim(seed)[1] == 1) stop("Error! Only one pair of seeds!")
  }
  if(!is.null(AB_dist)){
    if(min(AB_dist) < 0) stop("Error! Negative elements in distance matrix!")
    if(dim(AB_dist)[1] != n.A | dim(AB_dist)[2] != n.B) stop("Error! Matrices are not conformable!")
  }

  if (!is.null(seed)){
    seedmatch1 = apply(seed, 1, sum); #check the row summation
    seedmatch2 = apply(seed, 2, sum); #check the column summation
    if(max(seedmatch1, seedmatch2) > 1) stop("Error: no perfect matching for the seed")
    else{
      W = seeded.match(A, B, seed);
      Z = matrix(0, nrow = n.A, ncol = n.B);
      for (i in 1: n.A){
        tmp = sort(W[i,], decreasing = T)
        ind = which(W[i,] >= tmp[d])
        Z[i, ind] = 1
      }
      AB_dist = W;
    }
  }
  else{
    if (is.null(AB_dist)){
      AB_dist = DPdistance(A, B);
    }
    tmp = apply(AB_dist, 1, min);

    grid1 = c(0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8); grid2 = c(0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5);
    seed_list = array(0, dim=c(n.A, n.B, length(grid1)*length(grid2)));
    seed_match1 = rep(0, length(grid1)*length(grid2)); seed_match2 = rep(0, length(grid1)*length(grid2));

    for (j in 1: length(grid1)){
      tau1 = quantile(outa, grid1[j]);
      for (i in 1: length(grid2)){
        tau2 = quantile(tmp, grid2[i]);
        seed = seed.generate(A, B, tau1, tau2);
        seedmatch1 = apply(seed, 1, sum); #check the row summation
        seedmatch2 = apply(seed, 2, sum); #check the column summation
        seed_match1[(j-1)*length(grid2) + i] = max(seedmatch1);
        seed_match2[(j-1)*length(grid2) + i] = max(seedmatch2);
        seed_list[, , (j-1)*length(grid2) + i] = seed;
      }
    }

    temp1 = seed_match1 + seed_match2;

    if (2 %in% temp1 == TRUE){
      temp2 = which(temp1 == 2);
      perm_num = c();
      for (i in 1: length(temp2)){
        seed = seed_list[, , temp2[i]];
        gseed <- graph.incidence(seed);
        perm <- as_edgelist(gseed);
        perm_num[i] = dim(perm)[1];
      }
      seed = seed_list[, , temp2[which.max(perm_num)]];
      if (max(perm_num) == 1){
        warning("Warning: Only one pair of seeds. Run EE instead.")
        Z = matrix(0, nrow = n.A, ncol = n.B);
        for (i in 1: n.A){
          tmp = sort(AB_dist[i,], decreasing = F)
          ind = which(AB_dist[i,] <= tmp[d])
          Z[i, ind] = 1
        }
      }
      else{
        W = seeded.match(A, B, seed);
        Z = matrix(0, nrow = n.A, ncol = n.B);
        for (i in 1: n.A){
          tmp = sort(W[i,], decreasing = T)
          ind = which(W[i,] >= tmp[d])
          Z[i, ind] = 1
        }
      }
    }

    if (2 %in% temp1 == FALSE){
      warning("Warning: No perfect matching for the seed. Run EE instead.")
      Z = matrix(0, nrow = n.A, ncol = n.B);
      for (i in 1: n.A){
        tmp = sort(AB_dist[i,], decreasing = F)
        ind = which(AB_dist[i,] <= tmp[d])
        Z[i, ind] = 1
      }
    }
  }
  return(list(Dist = AB_dist, Z = Z))
}

seed.generate <- function(A, B, tau1, tau2){
  # tau1: nodes in graph A and Bwith degree larger than tau1 will be selected as seeds
  # tau2: pairs with distance smaller than tau3 will be regarded as true connections in the seed set

  # Output:
  # seed: the selected seed pairs

  n.A = dim(A)[1]
  n.B = dim(B)[1]

  outa = apply(A, 1, sum);
  outb = apply(B, 1, sum);

  seed = matrix(0, nrow = n.A, ncol = n.B);
  seedA = which(outa >= tau1); seedB = which(outb >= tau1);
  for(i in seedA){
    temp_a = outa[A[i, ] != 0];
    temp = rep(tau2+1, n.B);
    temp[seedB] = apply(B[seedB, ], 1, function(x) wasserstein1d(temp_a, outb[x!=0]));
    seed[i,] = 1*(temp <= tau2);
  }
  return(seed)
}

seeded.match <- function(A, B, seed, n, tau3 = NULL){
  #A: adjacency matrix of graph A;
  #B: adjacency matrix of graph B;
  #seed: a matrix with number of rows same as A and number of columns same as B, entries as 0 or 1,
  ##     where 1 indicates a connection between the two nodes. It defines a perfect matching on the subset of A and B.
  #tau3: the threshold for the common neighbors between two nodes.

  #Output:
  # w: the score between every pair. It can be used to check the correct matching.
  colnames(seed) = 1:dim(B)[1]; rownames(seed) = 1:dim(A)[1]
  gseed <- graph.incidence(seed)
  perm <- as_edgelist(gseed)

  rownames(A) = 1:dim(A)[1]; rownames(B) = 1:dim(B)[1];
  Acomp = A[-as.integer(perm[,1]), as.integer(perm[,1])];
  Bcomp = t(B[-as.integer(perm[,2]), as.integer(perm[,2])])
  H1 = Acomp%*%Bcomp;
  H = H1;
  H[H != 0] = 0;

  if (is.null(tau3)){
    temp = as.vector(H1);
    tau3 = quantile(temp, 1-1/max(dim(A)[1], dim(B)[1]));
  }

  H[H1 >= tau3] <- 1;

  gH <- graph.incidence(H, weighted = TRUE);
  match1 <- max_bipartite_match(gH);

  #Deal with the non-matching pairs
  matcha <- match1$matching[1:dim(Acomp)[1]]
  matchb <- match1$matching[-c(1:dim(Acomp)[1])];
  matcha[is.na(matcha)] <- names(matchb)[is.na(matchb)];

  permformat = as.vector(perm[,2]); names(permformat) = perm[,1];
  perm1 = c(matcha, permformat);
  A1 = A[,as.integer(names(perm1))]; B1 = B[as.integer(perm1),];
  w = A1%*%B1;

  return(w)
}

