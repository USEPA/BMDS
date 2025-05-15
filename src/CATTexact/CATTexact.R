#' Conditional exact Cochran-Armitage trend test
#'
#' \code{catt_exact} calculates the Cochran-Armitage trend test statistic (Cochran (1954), Armitage (1955)) and the one-sided p-value for the corresponding conditional exact test.
#' The conditional exact test has been established by Williams (1988). The computation of its p-value is performed using an algorithm following an idea by Mehta, et al. (1992).
#'
#' @param dose.ratings A vector of dose ratings, the i-th entry corresponds to the dose-rating of the i-th group. This vector must be strictly monotonically increasing
#' @param totals The vector of total individuals per group, the i-th entry corresponds to the total number of individuals in the i-th group.
#' @param cases The vector of incidences per groups, the i-th entry corresponds to the number of incidences in the i-th group.
#' @return A list containing the value of the Cochran-Armitage Trend Test Statistic, its exact and asymptotic p-value.
#' @references Armitage, P. Tests for linear trends in proportions and frequencies. \emph{Biometrics}, 11 (1955): 375-386.
#' @references Cochran, W. G. Some methods for strengthening the common \eqn{\chi^2} tests, \emph{Biometrics}. 10 (1954): 417-451.
#' @references Mehta, C. R., Nitin P., and Pralay S. Exact stratified linear rank tests for ordered categorical and binary data. \emph{Journal of Computational and Graphical Statistics}, 1 (1992): 21-40.
#' @references Portier, C., and Hoel D. Type 1 error of trend tests in proportions and the design of cancer screens. \emph{Communications in Statistics-Theory and Methods}, 13 (1984): 1-14.
#' @references Williams, D. A. Tests for differences between several small proportions. \emph{Applied Statistics}, 37 (1988): 421-434.
#' @examples
#' d <- c(1,2,3,4)
#' n <- rep(20,4)
#' r <- c(1,4,3,8)
#'
#' catt_exact(d, n, r)
#'
#' @export


catt_exact <- function(dose.ratings,totals,cases) {

  le.d <- length(dose.ratings)
  le.n <- length(totals)
  le.r <- length(cases)


  if (le.d != le.n | le.d != le.r) {
    stop("Length of input is differing!")
  }

  le <- le.d

  if (le < 3) {stop("Need at least three groups")}
  ## Extract Input, calculate total number of cases and individuals
  nk <- totals
  nhat <- sum(nk)

  rk <- cases
  rhat <- sum(rk)

  dk <- dose.ratings

  # dk <- dk/dk[2]
  #
  # rest <- dk - floor(dk)
  # rest <- rest[rest>0]
  #
  # mult <- min(1/(prod(rest)), 1)
  #
  # dk <- round(dk * mult)


  ## Input checks
  if (min(as.numeric(round(c(nk,rk)) == c(nk,rk))) == 0) {
    stop("The number of totals and cases must be integer")
  }

  if (min(as.numeric(nk > 0)) == 0) {
    stop("There must be at least one individual in every dose group")
  }

  if (min(as.numeric(rk >= 0)) == 0) {
    stop("The number of cases in each group must be nonnegative")
  }

  if (min(as.numeric(nk >= rk)) == 0) {
    stop("The number of cases can not exceed the size of the group")
  }

  if (max(as.numeric(rk > 0)) == 0) {
    stop("This test can not be applied, when there is no case")
  }

  if (min(as.numeric(nk == rk)) == 1) {
    stop("This test can not be applied, when the number of cases equals the total number of individuals")
  }

  check.dosemonvec <- rep(1, le - 1)
  for (i in 1:(le - 1))
  {check.dosemonvec[i] <- as.numeric(dk[i + 1] > dk[i])}

  check.dosemon <- min(check.dosemonvec)

  if (check.dosemon == 0) {
    stop("Doses must be strictly monotonically increasing")
  }



  factor <- sqrt(nhat / ( (nhat-rhat) * rhat))

  enum <- sum( (rk- (nk / nhat) * rhat) * dk)

  denom <- sqrt(sum( (nk / nhat) * dk ^ 2)-sum( (nk / nhat) * dk) ^ 2)

  test_statistic <- -factor * enum / denom

  .pval_exact <- .pval_exact(dk, nk, rk)
  pval_asy <- .aspvalue(test_statistic)

  return(list("test.statistic" = test_statistic, "exact.pvalue" = .pval_exact, "asymptotic.pvalue" = pval_asy))
}


#' Asymptotic Cochran-Armitage trend test
#'
#' \code{catt_asy} calculates the Cochran-Armitage trend test statistic (Cochran (1954), Armitage (1955)) and the one-sided p-value for the corresponding asymptotic test.
#' The exact form of used test statistic can be found in the paper by Portier and Hoel (1984).
#'
#' @param dose.ratings A vector of dose ratings, the i-th entry corresponds to the dose-rating of the i-th group. This vector must be strictly monotonically increasing
#' @param totals The vector of total individuals per group, the i-th entry corresponds to the total number of individuals in the i-th group
#' @param cases The vector of incidences per groups, the i-th entry corresponds to the number of incidences in the i-th group
#' @return A list containing the value of the Cochran-Armitage Trend Test Statistic and its asymptotic p-value.
#' @references Armitage, P. Tests for linear trends in proportions and frequencies. \emph{Biometrics}, 11 (1955): 375-386.
#' @references Cochran, W. G. Some methods for strengthening the common \eqn{\chi^2} tests, \emph{Biometrics}. 10 (1954): 417-451.
#' @references Portier, C., and Hoel D. Type 1 error of trend tests in proportions and the design of cancer screens. \emph{Communications in Statistics-Theory and Methods}, 13 (1984): 1-14.
#' @examples
#' d <- c(1,2,3,4)
#' n <- rep(20,4)
#' r <- c(1,4,3,8)
#'
#' catt_asy(d, n, r)
#'
#' @export

catt_asy <- function(dose.ratings, totals, cases) {

  le.d <- length(dose.ratings)
  le.n <- length(totals)
  le.r <- length(cases)

  if (le.d != le.n | le.d != le.r) {
    stop("Length of input is differing!")
  }

  le <- le.d

  if (le < 3) {stop("Need at least three groups")}
  ## Extract Input, calculate total number of cases and individuals
  nk <- totals
  nhat <- sum(nk)

  rk <- cases
  rhat <- sum(rk)

  dk <- dose.ratings

  ## Input checks
  if (min(as.numeric(round(c(nk,rk)) == c(nk,rk))) == 0) {
    stop("The number of totals and cases must be integer")
  }

  if (min(as.numeric(nk > 0)) == 0) {
    stop("There must be at least one individual in every dose group")
  }

  if (min(as.numeric(rk >= 0)) == 0) {
    stop("The number of cases in each group must be nonnegative")
  }

  if (min(as.numeric(nk >= rk)) == 0) {
    stop("The number of cases can not exceed the size of the group")
  }

  if (max(as.numeric(rk > 0)) == 0) {
    stop("This test can not be applied, when there is no case")
  }

  if (min(as.numeric(nk == rk)) == 1) {
    stop("This test can not be applied, when the number of cases equals the total number of individuals")
  }

  check.dosemonvec <- rep(1, le - 1)
  for (i in 1:(le - 1))
  {check.dosemonvec[i] <- as.numeric(dk[i + 1] > dk[i])}

  check.dosemon <- min(check.dosemonvec)

  if (check.dosemon == 0) {
    stop("Doses must be strictly monotonically increasing")
  }




  factor <- sqrt(nhat / ( (nhat-rhat) * rhat))

  enum <- sum( (rk- (nk / nhat) * rhat) * dk)

  denom <- sqrt(sum( (nk / nhat) * dk ^ 2)-sum( (nk / nhat) * dk) ^ 2)

  test_statistic <- -factor * enum / denom

  pval_asy <- .aspvalue(test_statistic)

  return(list("test.statistic" = test_statistic, "asymptotic.pvalue" = pval_asy))
}



.pval_exact <- function(dk, nk, rk) {

  dk <- dk/dk[2]

  rest <- dk - floor(dk)
  rest <- rest[rest>0]

  mult <- min(1/(prod(rest)), 10 ^ 12)

  dk <- round(dk * mult)

  le <- length(dk)
  nodes <- vector("list", le)
  nhat <- sum(nk)
  rhat <- sum(rk)
  nodes[1] <- 0
  a0 <- sum(rk * dk)



  # Nodes are created

  for (i in 1:(le - 1)) {
    lowerbound <- max(0, rhat - sum(nk[ (i + 1):le]))
    upperbound <- min(rhat, sum(nk[1:i]))

    nodes[[i + 1]] <- lowerbound:upperbound
  }

  nodes[[le + 1]] <- rhat
  arcs <- vector("list", le)


  # Arcs are created

  for (i in 1:le) {
   for (j in nodes[[i]]) {
    for (k in max(j, min(nodes[[i + 1]])):min(max(nodes[[i + 1]]), j + nk[i])) {
      arcs[[i]] <- c(arcs[[i]], j, k, dk[i] * (k - j), choose(nk[i], k - j))
    }

   }
   arcs[[i]] <- matrix(arcs[[i]], ncol = 4, byrow = TRUE)
  }


  # Zeros are added in the nodes

  nodes[[le + 1]] <- matrix(c(rhat, 0), ncol = 2)

  for (i in 1:le) {
    nodes[[i]] <- matrix(c(nodes[[i]], rep(0, length(nodes[[i]]))), ncol = 2)
  }

  # Backwards processing for calculating longest paths

  for (i in le:1) {
    for (j in nodes[[i]][,1]) {
      # Choose concurring arcs
      arckonkur <- matrix(arcs[[i]][ (which(arcs[[i]][ ,1] == j)),], ncol = 4)

      # Arcs get "consecutive" longest paths
      for (k in 1:length(arckonkur[ ,1])) {
        arckonkur[k,4] <- nodes[[i + 1]][which(nodes[[i+1]][ ,1]==arckonkur[k,2]),2]
      }

      # LP is calculated
      nodes[[i]][which(nodes[[i]][ ,1] == j),2] <- max(arckonkur[ ,4]+arckonkur[ ,3])
    }
  }


  # Two lists to express the tuples of lambdas over the nodes

  nodes.u <- vector("list",le + 1)


  for (i in 1:(le + 1)){
    nodes.u[[i]] <- vector("list",length(nodes[[i]][ ,1]))
  }

  nodes.u[[1]][[1]] <- 0

  nodes.cu <- vector("list",le + 1)

  for (i in 1:(le + 1)) {
    nodes.cu[[i]] <- vector("list",length(nodes[[i]][ ,1]))
  }

  nodes.cu[[1]][[1]] <- 0


  # Prespecify sets for first nodes

  nodes.u[[2]][1:length(nodes[[2]][ ,1])] <- arcs[[1]][ ,3]
  nodes.cu[[2]][1:length(nodes[[2]][ ,1])] <- arcs[[1]][ ,4]


  nodes.with.paths <- nodes[[2]][ ,1]


  for (i in 2:(le)) {
    nodes.with.paths.new <- numeric(0)
    for (j in nodes.with.paths) {
      # All	successor of a node are evaluated
      succ <- matrix(arcs[[i]][(which(arcs[[i]][ ,1]==j)),], ncol = 4)

      # u and c(u) are copied to the following nodes
      u.candidates <- matrix(c(succ,rep(nodes.u[[i]][[which(nodes[[i]][ ,1] == j)]], rep(length(succ) / 4, length(nodes.u[[i]][[which(nodes[[i]][ ,1] == j)]])))), nrow = length(succ) / 4)
      cu.candidates <- matrix(c(succ,rep(nodes.cu[[i]][[which(nodes[[i]][ ,1] == j)]], rep(length(succ) / 4, length(nodes.cu[[i]][[which(nodes[[i]][ ,1] == j)]])))), nrow = length(succ) / 4)

      # u and c(u) are transformed
      u.candidates[ ,5:ncol(u.candidates)] <- u.candidates[ ,5:ncol(u.candidates)] + succ[ ,3]
      cu.candidates[ ,5:ncol(u.candidates)] <- cu.candidates[ ,5:ncol(u.candidates)] * succ[ ,4]

      for (k in 1:(length(succ) / 4)) {
        candidate <- u.candidates[k,2]
        LP <- nodes[[i + 1]][which(nodes[[i + 1]][ ,1] == candidate),2]

        u.liste <- u.candidates[k,5:ncol(u.candidates)]
        cu.liste <- cu.candidates[k,5:ncol(u.candidates)]
        u.liste <- u.liste[(u.liste >= (a0 - LP - 1E-8))]
        cu.liste <- cu.liste[(u.liste >= (a0  - LP - 1E-8))]
        if (length(u.liste) > 0) {nodes.with.paths.new <-union(nodes.with.paths.new, candidate)}

        existing.u <- intersect(u.liste,nodes.u[[i + 1]][[which(nodes[[i + 1]][ ,1] == candidate)]])

        new.u <- setdiff(u.liste,nodes.u[[i + 1]][[which(nodes[[i + 1]][ ,1]==candidate)]])
        new.cu <-cu.liste[which(is.element(u.liste,new.u))]

        for (l in existing.u) {
          index <- which(nodes.u[[i + 1]][[which(nodes[[i + 1]][,1] == candidate)]] == l)  # index of existing l in nodes.u
          index <- which(nodes.u[[i + 1]][[which(nodes[[i + 1]][,1] == candidate)]] == l)  # index of existing l in nodes.u
          index2 <- which(u.liste == l)                                            # index of existing l in u.liste

          nodes.cu[[i + 1]][[which(nodes[[i + 1]][ ,1] == candidate)]][index] <- nodes.cu[[i + 1]][[which(nodes[[i + 1]][ ,1] == candidate)]][index]+cu.liste[index2]
        }

        nodes.u[[i + 1]][[which(nodes[[i + 1]][ ,1] == candidate)]] <- c(nodes.u[[i + 1]][[which(nodes[[i + 1]][ ,1] == candidate)]],new.u)
        nodes.cu[[i + 1]][[which(nodes[[i + 1]][ ,1] == candidate)]] <- c(nodes.cu[[i + 1]][[which(nodes[[i + 1]][ ,1] == candidate)]],new.cu)
      }

    }
    nodes.with.paths <- nodes.with.paths.new
  }
  pval <- sum(nodes.cu[[le + 1]][[1]]) / choose(nhat, rhat)

  return(pval)
}


#' @importFrom stats pnorm


.aspvalue <- function(statistic) {
  pval <- pnorm(statistic)
  return(pval)
}



