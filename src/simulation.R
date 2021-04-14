library(Rcpp)
sourceCpp("Beast.cpp")
source("Beast.R")

# testing bivariate independence
setting1 <- function(k, nsim=10){
  "Bivarite Normal Scenario"
  count = 0
  for (i in 1:nsim){
    epsilon1 <- rnorm(128)
    epsilon2 <- rnorm(128)
    epsilon3 <- rnorm(128)
    X <- sqrt((1.05 + 0.3 * k) / 2) * epsilon1 + sqrt((0.95 - 0.3 * k) / 2) * epsilon2
    Y <- sqrt((1.05 + 0.3 * k) / 2) * epsilon1 - sqrt((0.95 - 0.3 * k) / 2) * epsilon3
    res <- BEAST(cbind(X, Y), 4, test.independence = TRUE,
                 index = list(c(1), c(2)), test.uniformity = FALSE)
    
    if (res$p.value < 0.05){
      count = count + 1
    }
  }
  
  power = count / nsim
  return(power)
}

setting2 <- function(k, nsim=10){
  "Circle Scenario"
  count = 0
  for (i in 1:nsim){
    epsilon4 = rnorm(128)
    epsilon5 = rnorm(128)
    v = runif(128, -pi, pi)
    X = cos(v) + 4^(2/3) * (0.23 * k + 0.05) * epsilon4
    Y = sin(v) - 4^(2/3) * (0.23 * k + 0.05) * epsilon5
    res <- BEAST(cbind(X, Y), 4, test.independence = TRUE,
                 index = list(c(1), c(2)), test.uniformity = FALSE)
    
    if (res$p.value < 0.05){
      count = count + 1
    }
  }
  
  power = count / nsim
  return(power)
}

setting3 <- function(k, nsim=10){
  "Sine Scenario"
  count = 0
  for (i in 1:nsim){
    U = runif(128, -1, 1)
    epsilon6 = rnorm(128)
    X = U
    Y = sin(4 * pi * X) + (2.4 * k + 0.4) * epsilon6
    res <- BEAST(cbind(X, Y), 4, test.independence = TRUE,
                 index = list(c(1), c(2)), test.uniformity = FALSE)
    
    if (res$p.value < 0.05){
      count = count + 1
    }
  }
  
  power = count / nsim
  return(power)
}

setting4 <- function(k, nsim=10){
  "Checkerboard scenario"
  count = 0
  for (i in 1:nsim){
    epsilon7 = rnorm(128)
    epsilon8 = rnorm(128)
    W = sample(c(1, 2, 3), size=128, replace=TRUE)
    V1 = sample(c(2, 4), size=128, replace=TRUE)
    V2 = sample(c(1, 3, 5), size=128, replace=TRUE)
    X = W + (0.3*k + 0.05) * epsilon7
    Y = V1 * (W - 1) * (3 - W) + V2 * abs((W - 2)) + 4 * (0.3 * k + 0.05) * epsilon8
    res <- BEAST(cbind(X, Y), 4, test.independence = TRUE,
                 index = list(c(1), c(2)), test.uniformity = FALSE)
    
    if (res$p.value < 0.05){
      count = count + 1
    }
  }
  
  power = count / nsim
  return(power)
}

# testing univariate uniformity
setting5 <- function(k, nsim=10){
  "Stochastically ordered"
  count = 0
  for (i in 1:nsim){
    X <- rbeta(128, 0.7 + 0.24 * k, 1.3 - 0.24 * k)
    X <- matrix(X, ncol = 1)
    res <- BEAST(X, 5, test.independence = FALSE, test.uniformity = TRUE)
    if (res$p.value > 0.05){
      count = count + 1
    }
  }
  
  power = count / nsim
  return(power)
}

setting6 <- function(k, nsim=10){
  "Symmetric Unimodal"
  count = 0
  for (i in 1:nsim){
    X <- rbeta(128, 1.5 - 0.4 * k, 1.5 - 0.4 * k)
    X <- matrix(X, ncol = 1)
    res <- BEAST(X, 5, test.independence = FALSE, test.uniformity = TRUE)
    if (res$p.value > 0.05){
      count = count + 1
    }
  }
  power = count / nsim
  return(power)
}

setting7 <- function(k, nsim=10){
  "3-Modal"
  count = 0
  for (i in 1:nsim){
    F1 <- rbeta(128, 1, 4)
    F2 <- rbeta(128, 4, 1)
    F3 <- 0.5 + runif(128, -0.2 + 0.1 * k, 0.2 - 0.1 * k)
    X <- 1/3 * F1 + 1/3 * F2 + 1/3 * F3
    X <- matrix(X, ncol = 1)
    res <- BEAST(X, 5, test.independence = FALSE, test.uniformity = TRUE)
    if (res$p.value > 0.05){
      count = count + 1
    }
  }
  power = count / nsim
  return(power)
}

setting8 <- function(k, nsim=10){
  "8-Modal"
  count = 0
  for (i in 1:nsim){
    g1 = matrix(0, nrow=128, ncol=8)
    g2 = matrix(0, nrow=128, ncol=8)
    for (j in 1:8){
      g1[,j] <- runif(128, (2 * j - 1) / 16, j / 8)
      g2[,j] <- runif(128, (j - 1) / 8, (2 * j - 1) / 16)
    }
    G1 = c()
    G2 = c()
    for (n in 1:128){
      G1 = c(G1, sample(g1[n,], 1))
      G2 = c(G2, sample(g2[n,], 1))
    }
    pi1 = 0.525 + k / 4
    pi2 = 0.475 - k / 4
    X <- pi1 * G1 + pi2 * G2
    X <- matrix(X, ncol = 1)
    res <- BEAST(X, 5, test.independence = FALSE, test.uniformity = TRUE)
    if (res$p.value > 0.05){
      count = count + 1
    }
  }
  power = count / nsim
  return(power)
}

# testing independence of a variable and a vector
setting9 <- function(k, nsim=10){
  "Linear"
  count = 0
  for (i in 1:nsim){
    X1 = rnorm(128)
    X2 = rnorm(128)
    hk = sqrt(0.68 + 0.64 * k - 0.32 * k^2)
    epsilon1 = rnorm(128)
    Y = 0.4 * (1 - k) * (X1 + X2) + hk * epsilon1
    res <- BEAST(cbind(X1, X2, Y), 3, index=list(c(1, 2), c(3)), test.independence = TRUE, test.uniformity = FALSE)
    if (res$p.value < 0.05){
      count = count + 1
    }
  }
  power = count / nsim
  return(power)
}

setting10 <- function(k, nsim=10){
  "Sphere"
  count = 0
  for (i in 1:nsim){
    G1 = rnorm(128)
    G2 = rnorm(128)
    G3 = rnorm(128)
    epsilon2 = rnorm(128)
    G_norm = sqrt(G1^2 + G2^2 + G3^2)
    X1 = G1 / G_norm
    X2 = G2 / G_norm
    Y = G3 / G_norm + (27 * k + 9) * epsilon2 / 40
    res <- BEAST(cbind(X1, X2, Y), 3, index=list(c(1, 2), c(3)), test.independence = TRUE, test.uniformity = FALSE)
    if (res$p.value < 0.05){
      count = count + 1
    }
  }
  
  power = count / nsim
  return(power)
}

setting11 <- function(k, nsim=10){
  "Sine"
  count = 0
  for (i in 1:nsim){
    U1 = runif(128, 0, 1)
    U2 = runif(128, 0, 1)
    X1 = U1
    X2 = U2
    epsilon3 = rnorm(128)
    Y = sin(4 * pi * (X1 + X2)) + (9 * k + 1) * epsilon3 / 5
    res <- BEAST(cbind(X1, X2, Y), 3, index=list(c(1, 2), c(3)), test.independence = TRUE, test.uniformity = FALSE)
    if (res$p.value < 0.05){
      count = count + 1
    }
  }
  
  power = count / nsim
  return(power)
}

setting12 <- function(k, nsim=10){
  "Double Helix"
  count = 0
  for (i in 1:nsim){
    I = sample(c(-1, 1), size=128, replace=TRUE)
    v = runif(128, -pi, pi)
    epsilon4 = rnorm(128)
    epsilon5 = rnorm(128)
    epsilon6 = rnorm(128)
    c4 = (0.36 * k + 0.476) * epsilon4
    c5 = (0.36 * k + 0.476) * epsilon5
    X1 = I * cos(v) + c4
    X2 = I * sin(v) + c5
    Y = v + (0.36 * k + 0.476) * epsilon6
    U1 = runif(128, 0, 1)
    U2 = runif(128, 0, 1)
    X1 = U1
    X2 = U2
    epsilon3 = rnorm(128)
    Y = sin(4 * pi * (X1 + X2)) + (9 * k + 1) * epsilon3 / 5
    res <- BEAST(cbind(X1, X2, Y), 3, index=list(c(1, 2), c(3)), test.independence = TRUE, test.uniformity = FALSE)
    if (res$p.value < 0.05){
      count = count + 1
    }
  }
  
  power = count / nsim
  return(power)
}
