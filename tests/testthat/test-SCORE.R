require(stats)
require(testthat)

context("SCORE")

test_that("SCORE stops when it should", {
  expect_error( runningmean(0, c(0,0)) )
})

###############################################################
# Here we generate a toy SBM graph to test if SCORE works well.
###############################################################

set.seed(2020)
n = 10; K = 2;
P  = matrix(c(1/2, 1/4, 1/4, 1/2), byrow = T, nrow = K);

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

test_that("The result should be label for all nodes in the graph", {
  expect_length(SCORE(G, 2), n)
  expect_length(unique(SCORE(G, 2)), 2)
})

######################################################################################
# Here we generate a larger SBM graph, which is more practical in real-world scenerio.
######################################################################################

#set.seed(2020)
#n = 1000; K = 2;
#P  = matrix(c(1/5, 1/10, 1/10, 1/5), byrow = T, nrow = K);

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

#result = SCORE(G, 2)
