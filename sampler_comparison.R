library(tidyverse)
library(PKPDmodels)
library(mvtnorm)
library(rstan)
options(mc.cores = parallel::detectCores())

# prepare data for first patient
theoph <- Theoph %>% as_tibble %>% mutate(DoseMass = Dose*Wt)
theoph %>% ggplot(aes(Time, conc, color = Subject)) + geom_point() + geom_line()
theoph.1 <- theoph %>% filter(Subject == 1)
theoph.1 %>% ggplot(aes(Time, conc)) + geom_point() + geom_line()

dat <- list(Nt = nrow(theoph.1),
            ts = theoph.1 %>% pull(Time) %>% `[`(-1),
            y_init = c(theoph.1 %>% pull(DoseMass) %>% head(1),
                       theoph.1 %>% pull(conc) %>% head(1)),
            y = theoph.1 %>% pull(conc) %>% `[`(-1))

#################################
# Stan NUTS
#################################
fit <- stan(file = "theoph_single_individual.stan",
            data = dat,
            iter = 2000, chains = 2, control = list(adapt_delta = 0.9))

post.samples <- extract(fit)
nuts.samples <- tibble(Ka = post.samples$theta[,1], sigma = post.samples$sigma) %>% mutate(Method = "NUTS")

#################################
# Point Estimate
#################################

model <- stan_model(file = "theoph_single_individual.stan")
opt.object <- optimizing(model, data = dat, draws = 2000, hessian = TRUE)


opt.theta1.sigma <- opt.object$par[c("theta[1]", "sigma")]
cov.theta1.sigma <- opt.object$hessian[c("theta.1", "sigma"),c("theta.1", "sigma")]

point.estimate <- tibble(Ka = rep(opt.theta1.sigma[1],2000), sigma = rep(opt.theta1.sigma[2],2000)) %>%
  mutate(Method = "Point Estimate")

#################################
# Asymptotic Gaussian Approximation
#################################
gaussian.samples <- opt.object$theta_tilde[,c("theta[1]","sigma")] %>%
  as_tibble %>%
  select(Ka = `theta[1]`, sigma) %>%
  mutate(Method = "Asymptotic")

#################################
# Metropolis
#################################
full.post.samples <- tibble(b = post.samples$b,
       Ka = post.samples$theta[,1],
       Cl = post.samples$theta[,2],
       V = post.samples$theta[,3],
       c_init = post.samples$c_init,
       sigma = post.samples$sigma)

q0 <- full.post.samples %>% head(1) %>% as.matrix %>% as.vector %>% log

# log likelihood function
ll <- function(upars) log_prob(fit, upars, adjust_transform = TRUE, gradient = FALSE)

# use NUTS samples to get optimal Metropolis jump distribution
Sigma <- full.post.samples %>% mutate_all(log) %>% cov

RandomStep <- function() {
  rmvnorm(n=1, mean=rep(0,6), sigma=Sigma) %>% as.vector
}

GetMetropolisSamples <- function(ll, q0, Sigma, num.samples) {
  
  D <- length(q0)
  samples <- matrix(NA, nrow = num.samples+1, ncol = D)
  samples[1,] <- q0
  
  for(s in 1:num.samples) {
    
    proposal <- samples[s,] + RandomStep()
    if(runif(1) < exp(ll(proposal)-ll(samples[s,]))) {
      samples[s+1,] <- proposal
    } else{
      samples[s+1,] <- samples[s,]
    }
    
  }
  samples %>% as_tibble %>% set_names(c("b", "Ka", "Cl", "V", "c_init", "sigma"))
}

metropolis.samples <- GetMetropolisSamples(ll, q0, Sigma, num.samples = 2000-1) %>%
  select(Ka, sigma) %>%
  mutate_all(exp) %>%
  mutate(Method = "Metropolis")

#################################
# Bind and Plot
#################################
bind_rows(metropolis.samples, point.estimate, nuts.samples, gaussian.samples) %>%
  mutate(Method = factor(Method, levels = c("Point Estimate", "Asymptotic", "Metropolis", "NUTS"))) %>%
  ggplot(aes(Ka,sigma)) +
  geom_point(alpha = 0.2) +
  facet_wrap( ~ Method) +
  scale_x_log10() +
  scale_y_log10() +
  ylab(expression(sigma)) +
  theme(text = element_text(size=40))
