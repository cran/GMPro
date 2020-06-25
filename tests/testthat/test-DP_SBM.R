require(combinat)
require(stats)
require(testthat)
require(transport)
require(igraph)

context("DP_SBM")

test_that("DP_SBM stops when it should", {
  expect_error( runningmean(0, c(0,0)) )
})

#################################################################
# Here we generate 2 toy SBM graphs to test if DP_SBM works well.
#################################################################

set.seed(2020)

K = 2; n = 30; s = 1;
P  = matrix(c(1/2, 1/4, 1/4, 1/2), byrow = T, nrow = K);
# define community label matrix Pi
distribution = c(1, 2); # This should be adjusted to K
l = sample(distribution, n, replace=TRUE, prob = c(1/2, 1/2)); # node labels
Pi = matrix(0, n, 2) # label matrix
for (i in 1:n){
  Pi[i, l[i]] = 1
}

# define the expectation of the parent graph's adjacency matrix
Omega = Pi %*% P %*% t(Pi)

# construct the parent graph G
G = matrix(runif(n*n, 0, 1), nrow = n);
G = G - Omega;
temp = G;
G[which(temp >0)] = 0;
G[which(temp <=0)] = 1;
diag(G) = 0;
G[lower.tri(G)] = t(G)[lower.tri(G)];

##### Sample Graphs Generation #####
# generate graph A from G
dA = matrix(rbinom(n*n, 1, s), nrow = n, ncol=n);
dA[lower.tri(dA)] = t(dA)[lower.tri(dA)];
A1 = G*dA;
indA = sample(1:n, n, replace = FALSE);
labelA = l[indA];
A = A1[indA, indA];

# generate graph B from G
dB = matrix(rbinom(n*n, 1, s), nrow = n, ncol=n);
dB[lower.tri(dB)] = t(dB)[lower.tri(dB)];
B1 = G*dB;
indB = sample(1:n, n, replace = FALSE);
labelB = l[indB];
B = B1[indB, indB];

perm = permn(1:K)

test_that("This function returns a list of matching results for all possible permutations on {1, ..., K}", {
  expect_length(DP_SBM(A, B, K, fun = "DPmatching"), length(perm))
})

########################################################################################
# Here we generate 2 larger SBM graphs, which are more practical in real-world scenerio.
########################################################################################
#set.seed(2020)

#K = 2; n = 1000; s = 1;
#P  = matrix(c(1/5, 1/10, 1/10, 1/5), byrow = T, nrow = K);
# define community label matrix Pi
#distribution = c(1, 2); # This should be adjusted to K
#l = sample(distribution, n, replace=TRUE, prob = c(1/2, 1/2)); # node labels
#Pi = matrix(0, n, 2) # label matrix
#for (i in 1:n){
#  Pi[i, l[i]] = 1
#}

# define the expectation of the parent graph's adjacency matrix
#Omega = Pi %*% P %*% t(Pi)

# construct the parent graph G
#G = matrix(runif(n*n, 0, 1), nrow = n);
#G = G - Omega;
#temp = G;
#G[which(temp >0)] = 0;
#G[which(temp <=0)] = 1;
#diag(G) = 0;
#G[lower.tri(G)] = t(G)[lower.tri(G)];

##### Sample Graphs Generation #####
# generate graph A from G
#dA = matrix(rbinom(n*n, 1, s), nrow = n, ncol=n);
#dA[lower.tri(dA)] = t(dA)[lower.tri(dA)];
#A1 = G*dA;
#indA = sample(1:n, n, replace = FALSE);
#labelA = l[indA];
#A = A1[indA, indA];

# generate graph B from G
#dB = matrix(rbinom(n*n, 1, s), nrow = n, ncol=n);
#dB[lower.tri(dB)] = t(dB)[lower.tri(dB)];
#B1 = G*dB;
#indB = sample(1:n, n, replace = FALSE);
#labelB = l[indB];
#B = B1[indB, indB];

#result = DP_SBM(A, B, K, fun = "DPmatching")
