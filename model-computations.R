
# for W1 = W2 = w
Nstar <- function (w, b, alpha, sigma) {
  Q <- (w * b/alpha) ^ (1/sigma)
  w/(b + 2*Q)
}



# for N1 > 0, N2 > 0
N1 <- function (h, e, b, alpha, sigma) {
  Q <- (e * b / alpha)^(1/sigma)
  P <- (h - e * (h/e)^(1/sigma))/b
  R <- (h/e)^(1/sigma)
  
  (e - Q*P)/(b + Q*(1+R))
}

# for W1 = e, W2 = h, N1 > 0, N2 > 0
N2_from_N1 <- function (N1, h, e, b, sigma) {
  P <- (h - e * (h/e)^(1/sigma))/b
  R <- (h/e)^(1/sigma)
  
  P + R * N1
}

# for W1 = e, W2 = h, N1 > 0, N2 > 0
N_total <- function (h, e, b, alpha, sigma) {
  N1 <- N1(h = h, e = e, b = b, alpha = alpha, sigma = sigma) 
  N2 <- N2_from_N1(N1, h = h, e = e, b = b, sigma = sigma) 
  N1 + N2
}

plot_N1_N2 <- function (e, b, alpha, sigma, xlim = NULL, ylim = c(0, 6), ...) {
  if (is.null(xlim)) xlim <- c(e, e*5)
  curve(
          N1(x, e = e, b = b, alpha = alpha, sigma = sigma), 
          xlab = "h",
          ylab = "N",
          xlim = xlim,
          ylim = ylim,
          col = "darkblue",
          sub = glue::glue("e = {e}, b = {b}, alpha = {alpha}, sigma = {sigma}"),
          ...
        )
  curve(N2_from_N1(
    N1(x, e = e, b = b, alpha = alpha, sigma = sigma),
    h = x, e = e, b = b, sigma = sigma), add = TRUE, col = "red")
  curve(N_total(x, e = e, b = b, alpha = alpha, sigma = sigma), add = TRUE)
}


# total utility when subutility is log
# 

N_low_educ <- function (w, alpha, b) {
  2 * w / (b * (1+2*w/alpha))   
}

u_low_educ <- function (w, alpha, b) {
  2 * log(w) + 
  2 * log(1 - w/(1+2*w/alpha)) +
  alpha * log(1/b) +
  alpha * log( 2*w / (1 + 2*w/alpha) )
}

u_high_educ <- function (h, e, alpha, b) {
  log(e) + 
  log( 1 - e/(1 + (e+h)/alpha) ) +
  log(h) + 
  log(1 - h/(1 + (e+h)/alpha)) +
  alpha * log( (e + h) / (b * (1 + (e+h)/alpha) ) )
}

plot_us <- function (w, e, alpha, b, xlim = NULL, ...) {
  if (is.null(xlim)) xlim <- c(e, 5*e)
  curve(u_high_educ(h, e = e, alpha = alpha, b = b), xname = "h", 
          xlim = xlim, ...)
  u_lo <- u_low_educ(w = w, alpha = alpha, b = b)
  # dotted line shows utility of low-educ path. To the R of crossover,
  # people choose high-educ
  abline(h = u_lo, lty = 2)
  # red line shows value of h where higher educ ppl have more children
  abline(v = 2 *w - e, col = "red")
}


plot_us(w = 0.4, e = 0.3, alpha = 0.2, b = 0.8, ylim = c(-4,-2), xlim = c(0.4, 0.6))



# == Single period model ==
# 

Q <- function (W, alpha, b, m, sigma) {
  ((m + b * W) / alpha) ^ (1/sigma)
}

Nstar <- function (W, alpha, b, m, sigma) {
  Q <- Q(W, alpha, b, m, sigma)
  W / (Q + W * b + m)
}

dN_dW <- function (W, alpha, b, m, sigma) {
  Q <- Q(W, alpha, b, m, sigma)
  
  1/(Q + b*W + m)^2 * 
                      (Q + m - W * (1/sigma)*Q*b/(b*W + m))
}

# ==== linear utility in children ====

Nstar <- function (W, a, b, sigma) {
  1/b - b^((1 - sigma)/sigma) / a^(1/sigma) * W^((1 - sigma)/sigma)
}

dN_dW <- function (W, a, b, sigma) {
  -b^((1 - sigma)/sigma) / a^(1/sigma) *
  (1 - sigma) / sigma *
  W^((1 - 2 * sigma)/sigma)
}


# ---- version with education ----

e <- function (gamma, eta, h, sigma) {
  R <- (gamma +  eta * h) ^ (sigma/(1 - 2 * sigma))
  (1 - R) * (1 + R * (gamma + eta * h))
}


# ---- version with w = eh ----

N2 <- function (h, a, b, sigma) {
  N2_pos <- 1/b * (1 - (b/a)^(1/(2 * sigma - 1)) * h^((1-sigma)/(2 * sigma - 1)))
  pmax(N2_pos, 0)
}

N1 <- function (h, a, b, sigma) {
  N1_pos <- N2(h, a, b, sigma) - (1/b) * (b/a)^(1/sigma)
  pmax(N1_pos, 0)
}


N <- function (h, a, b, sigma) {
  N1(h, a, b, sigma) + N2(h, a, b, sigma) 
}

curve(N1(h, 0.1, 0.1, 0.8), xname = "h", xlim = c(0, 1.2), col = "darkred", lwd = 2, lty = 2)
curve(N2(h, 0.1, 0.1, 0.8), xname = "h", add = TRUE, col = "darkgreen", lwd = 2, lty =2)
curve(N(h, 0.1, 0.1, 0.8), xname = "h", add = TRUE)


# ---- when N1 = 0 ----


f <- function (e, h, a, b, sigma) {
 lhs <- (1 - e) * e ^ ((1 - 2 * sigma) / sigma^2)
 rhs <- (a / (b*h)) ^ ((1 - sigma) / sigma^2)
 lhs - rhs
}

e_star <- function (h, a, b, sigma) {
  ur <- uniroot(f, c(1e-6, 2), h = h, a = a, b = b, sigma = sigma)
  ur$root
}

N2_from_e <- function (h, a, b, sigma) {
  e <- Vectorize(e_star, "h")(h, a, b, sigma)
  1/b * (1 - (b/a) ^ (1/sigma) * (e * h) ^ ((1-sigma)/sigma))
}

curve(N2_from_e(h, 0.1, 0.1, 0.8), xname = "h", xlim = c(0.7, 1.2), add = TRUE)
# both 0

# ==== raw utility ====
# 

U <- function (params, a, b, sigma, h) {
  N1 <- params[1]
  N2 <- params[2]
  e  <- params[3]
  
  u1 <- (1 - e - b * N1)
  u1 <- u1^(1 - sigma) / (1 - sigma)
  
  u2 <- h * e * (1 - b * N2)
  u2 <- u2^(1 - sigma) / (1 - sigma)
  
  u1 + u2 + a * (N1 + N2)
}

U_alt <- function(N1, N2, e, a, b, sigma, h) {
  U(c(N1, N2, e), a = a, b = b, sigma = sigma, h = h)
}

