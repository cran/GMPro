require(stats)
require(testthat)
require(transport)
require(igraph)

context("EEpost")

test_that("EEpost stops when it should", {
  expect_error( runningmean(0, c(0,0)) )
})

########################################################################
# Here we generate 2 toy ER graphs to test if EEpost works well.
########################################################################

set.seed(2020)
n = 10;p = 1; q = 1/2; s = 1;

Parent = matrix(rbinom(n*n, 1, q), nrow = n, ncol = n);
Parent[lower.tri(Parent)] = t(Parent)[lower.tri(Parent)];
diag(Parent) <- 0;

#Generate graph A
dA = matrix(rbinom(n*n, 1, s), nrow = n, ncol=n);
dA[lower.tri(dA)] = t(dA)[lower.tri(dA)];
A1 = Parent*dA;
tmp = rbinom(n, 1, p)
n.A = length(which(tmp == 1))
indA = sample(1:n, n.A, replace = FALSE)
A = A1[indA, indA];

#Generate graph B
dB = matrix(rbinom(n*n, 1, s), nrow = n, ncol=n);
dB[lower.tri(dB)] = t(dB)[lower.tri(dB)];
B1 = Parent*dB;
tmp = rbinom(n, 1, p)
n.B = length(which(tmp == 1))
indB = sample(1:n, n.B, replace = FALSE)
B = B1[indB, indB];

matching1= DPmatching(A, B)$Dist

test_that("This function returns a distance matrix and a matching indicator matrix", {
  expect_equal(dim(EEpost(A = A, B = B, d = 5, rep = 10)$Dist), c(n.A, n.B))
  expect_length(EEpost(A = A, B = B, d = 5, rep = 10)$match, n.A)
  expect_length(EEpost(A = A, B = B, d = 5, rep = 10)$FLAG, n.A)
})

test_that("This function returns a distance matrix and a matching indicator matrix", {
  expect_equal(dim(EEpost(A = A, B = B, d = 5, rep = 10, matching = matching1)$Dist), c(n.A, n.B))
  expect_length(EEpost(A = A, B = B, d = 5, rep = 10, matching = matching1)$match, n.A)
  expect_length(EEpost(A = A, B = B, d = 5, rep = 10, matching = matching1)$FLAG, n.A)
})

#######################################################################################
# Here we generate 2 larger ER graphs, which are more practical in real-world scenerio.
#######################################################################################

#set.seed(2020)
#n = 300;p = 1; q = 1/8; s = 0.99;

#Parent = matrix(rbinom(n*n, 1, q), nrow = n, ncol = n);
#Parent[lower.tri(Parent)] = t(Parent)[lower.tri(Parent)];
#diag(Parent) <- 0;

#Generate graph A
#dA = matrix(rbinom(n*n, 1, s), nrow = n, ncol=n);
#dA[lower.tri(dA)] = t(dA)[lower.tri(dA)];
#A1 = Parent*dA;
#tmp = rbinom(n, 1, p)
#n.A = length(which(tmp == 1))
#indA = sample(1:n, n.A, replace = FALSE)
#A = A1[indA, indA];

#Generate graph B
#dB = matrix(rbinom(n*n, 1, s), nrow = n, ncol=n);
#dB[lower.tri(dB)] = t(dB)[lower.tri(dB)];
#B1 = Parent*dB;
#tmp = rbinom(n, 1, p)
#n.B = length(which(tmp == 1))
#indB = sample(1:n, n.B, replace = FALSE)
#B = B1[indB, indB];

#matching1= DPmatching(A, B)$Dist

#result1 = EEpost(A = A, B = B, d = 10, rep = 50)
#result2 = EEpost(A = A, B = B, d = 10, rep = 50, matching = matching1)
