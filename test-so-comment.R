
N <- 1e6

s2_x_star <- 2
x_star <- rnorm(N, sd = sqrt(s2_x_star))

s2_eta <- 2
eta <- rnorm(N, sd = sqrt(s2_eta))

x <- x_star + eta


s2_epsilon <- 2.5
epsilon <- rnorm(N, sd = sqrt(s2_epsilon))

beta <- .3

y <- beta * x_star + epsilon

est <- lm(y ~ x)

beta_hat <- coef(est)[2]

R2_hat <- summary(est)$r.squared
R2 <- summary(lm(y ~ x_star))$r.squared

glue::glue("Estimated beta: {round(beta_hat, 4)}")
glue::glue("Predicted beta_hat: {round(s2_x_star/(s2_x_star+s2_eta) * beta, 4)}")

glue::glue("Estimated R2: {round(R2_hat, 4)}")
glue::glue("Predicted beta_hat: {round(s2_x_star/(s2_x_star+s2_eta) * R2, 4)}")
