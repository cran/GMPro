require(testthat)
require(transport)

context("DPdistance")

test_that("DPdistance stops when it should", {

  expect_error( runningmean(0, c(0,0)) )

})

###################################################################
# Here we generate 2 toy ER graphs to test if DPmatching works well.
###################################################################

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

test_that("This function returns a distance matrix", {
  expect_equal(dim(DPdistance(A, B))[1], n.A)
  expect_equal(dim(DPdistance(A, B))[2], n.B)
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

#result = DPdistance(A, B)
