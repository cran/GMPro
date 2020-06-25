require(stats)
require(testthat)
require(transport)

context("DPedge")

test_that("DPedge stops when it should", {

  expect_error( runningmean(0, c(0,0)) )

})

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

W = DPdistance(A, B)

test_that("matching result is the number of nodes in A", {
  expect_equal(dim(DPedge(A, B, 5)$Z), c(n.A, n.B))
})

test_that("produce errors when an error is expected", {
  expect_error(DPedge(d = 5), "Error! No values for A, B and W!")
})

test_that("the dimension of distance matrix is the number of nodes in A and B respectively", {
  expect_equal(dim(DPedge(d = 5, W = W)$Z), c(n.A, n.B))
})

#######################################################################################
# Here we generate 2 larger ER graphs, which are more practical in real-world scenerio.
#######################################################################################

#set.seed(2020)
#n = 300; p = 1; q = 1/8; s = 0.99;

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

#result = DPedge(A, B, 10)
