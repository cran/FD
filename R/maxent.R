maxent <- function (c.means, c.mat, prior, tol = 1e-08, lambda = F){
     
     # check input
    if (!is.numeric(c.means)) stop("c.means must be a numeric vector\n")
    if (!is.matrix(c.mat)){
       if (is.numeric(c.mat)){
          s.names <- names(c.mat)
          c.mat <- matrix(c.mat, nrow = 1, ncol = length(c.mat) )
          dimnames(c.mat) <- list("constraint", s.names)
         }
       else stop("if c.mat is not a matrix then it can only be a numeric vector\n")
     }
    s.names <- dimnames(c.mat)[[2]]
    c.names <- dimnames(c.mat)[[1]]
    dim.matrix <- dim(c.mat)
    if (dim.matrix[2] == 1 && dim.matrix[1] > 1){
        c.mat <- t(c.mat)
        dim.matrix <- dim(c.mat)
     }
    n.species <- dim.matrix[2]
    if (missing(prior) ) prob <- rep(1 / n.species, n.species) else prob <- prior
    if (length(prob) != n.species) stop("number of states in prior not equal to number in c.mat\n")
    if (any(is.na(c.means)) || any(is.na(c.mat)) || any(is.na(prob)) ) stop("no NA's allowed\n")
    n.constraints <- length(c.means)
    if (n.constraints != dim.matrix[1]) stop("number of constraint means not equal to number of constraints in c.mat\n")

    # run algorithm
    C.values <- rowSums(c.mat)
    test <- 1e+10
    iter <- 0
    while (test > tol){
        iter <- iter + 1
        denom <- c.mat %*% prob
        gamma.values <- log(c.means / denom) / C.values
        if (any(is.na(gamma.values) ) ) stop("NA's in gamma.values\n")
        unstandardized <- exp(t(gamma.values) %*% c.mat) * prob
        new.prob <- as.vector(unstandardized / sum(unstandardized) )
        if (any(is.na(new.prob) ) ) stop("NA's in new.prob\n")
        test <- max(abs(prob - new.prob))
        prob <- new.prob
     }

  # output
  res <- list()
  names(prob) <- s.names
  res$prob <- prob
  moments <- as.numeric(denom)
  names(moments) <- rownames(denom)
  res$moments <- moments
  res$entropy <- -1 * sum(prob * log(prob) )
  res$iter <- iter
  if (lambda){
    lambda <- coef(lm(log(prob) ~ t(c.mat) ) )
    names(lambda) <- c("intercept", c.names)
    res$lambda <- lambda
   }
  return(res)
}

