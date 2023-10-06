if (return_parameters && is.numeric(continuous)) {
  for (i in 1:vZ) {
    trace <- matrix(NA, nrow = G, ncol = S - burn)
    colnames(trace) <- (burn + 1):(S)
    for (j in 1:G) {
      trace[j, ] <- beta[j, (burn + 1):S, i]
    }
    traceplots[paste0("beta_", i)] <- trace
  }
}

if (return_parameters) {
  return(list(S = S, burn = burn, lambda = lambda, media.lambda = media.lambda,
              alfa = alfa, mu = mu, theta = theta, tau.theta = tau.theta,
              phi = phi, tau.phi = tau.phi, epsilon = epsilon, tau.e = tau.e,
              traceplots = traceplots))
} else {
  return(list(S = S, burn = burn, lambda = lambda, media.lambda = media.lambda,
              traceplots = traceplots))
}
}
