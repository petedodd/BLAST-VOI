
## function to compute benefits from partial screening
benefit <- function(grads, # slopes
                    N, # pops
                    n, # total screened
                    rnk = NULL, # to set externally derived ranks
                    verbose = FALSE # speak?
) {
  if (is.null(rnk)) {
    rnk <- order(grads, decreasing = TRUE)
  } # ranking based on true grads
  bnfts <- N * grads # available benefits
  cpops <- cumsum(N[rnk]) # cumulative pops in rnk order
  if (n > rev(cpops)[1]) {
    stop("screening > population!")
  }
  lastcomplete <- sum(n > cpops) # which is last community completed
  if (verbose) {
    cat("...lastcomplete = ", lastcomplete, "...\n")
  }
  if (lastcomplete > 0) { # get through any comms?
    ans <- sum(bnfts[rnk[1:lastcomplete]]) # use complete as base
    nleft <- (n - cpops[lastcomplete]) # n left for next zone
    fracdone <- nleft / N[rnk[lastcomplete + 1]] # fraction of last comm done
  } else {
    if (verbose) {
      cat("...(1st comm termination)\n")
    }
    ans <- 0 # use zero as base
    fracdone <- n / N[rnk[1]] # fraction of last=1st comm done
  }
  if (verbose) {
    cat("...complete zone benefits = ", ans, "...\n")
  }
  if (verbose) {
    cat("...fracdone = ", fracdone, "...\n")
  }
  if (lastcomplete < length(N)) { # add on partial bit for last comm
    if (verbose) {
      cat(
        "...stopping in zone = ", rnk[lastcomplete + 1],
        ", with slope = ", grads[rnk[lastcomplete + 1]], "...\n"
      )
    }
    ans <- ans +
      bnfts[rnk[lastcomplete + 1]] * # total benefit available in last comm
        fracdone # fraction of last comm done
    ## same as adding grads[rnk[lastcomplete+1]] * nleft
  }
  if (verbose) {
    cat("...ans = ", ans, "...\n")
  }
  ans
}


DbenefitFromSlps <- function(slps_true, # true slopes
                             slps_est, # estimated slopes: NOTE only rank used
                             N, # population sizes
                             n, # intervention pop covered
                             verbose = FALSE, # print stuff
                             separate = FALSE) {
  rnk_true <- order(slps_true, decreasing = TRUE) # ranking based on true grads
  rnk_est <- order(slps_est, decreasing = TRUE) # ranking based on true grads
  if (verbose) {
    print(rnk_true)
    print(rnk_est)
    cat("BENEFITS (TRUTH):\n")
  }
  bnft_true <- benefit(slps_true, N, n, rnk_true, verbose) # benefits using true parameters
  if (verbose) {
    print(bnft_true)
  }
  if (verbose) {
    cat("BENEFITS (EST):\n")
  }
  bnft_est <- benefit(slps_true, N, n, rnk_est, verbose) # benefits using estimating parameters
  if (verbose) {
    print(bnft_est)
  }
  if (separate) {
    ans <- list(bnft_opt = bnft_true, bnft_sub = bnft_est)
  } else {
    ans <- bnft_true - bnft_est # diff in benefit from these
  }
  ans
}

## for tables etc
gr <- function(m, l, h) {
  m <- formatC(m, format = "f", digits = 2)
  l <- formatC(l, format = "f", digits = 2)
  h <- formatC(h, format = "f", digits = 2)
  glue::glue("{m} ({l} to {h})")
}

## parameters for converting deaths into DALYs
daly <- list(
  r = 3e-2, #discount rate
  LE = 35   #life-expectancy
)
