library(tidyverse)
library(PKPDmodels)
library(deSolve)
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

# get posterior samples
fit <- stan(file = "theoph_single_individual.stan",
            data = dat,
            iter = 2000, chains = 2, control = list(adapt_delta = 0.9))

post.samples <- extract(fit)
post.samples <- tibble(b = post.samples$b,
                       Ka = post.samples$theta[,1],
                       Cl = post.samples$theta[,2],
                       V = post.samples$theta[,3],
                       c_init = post.samples$c_init,
                       sigma = post.samples$sigma)

# forward simulate synthetic data with random parameter values from post.samples
ForwardSymTheo <- function(b, Ka, Cl, V, c_init, sigma) {
  
  pkpd.two.compartment <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      dA <-  -Ka*A
      dc <- (1/V)*(Ka*A - (Cl/V)*(c*V))
      list(c(dA, dc))
    })
  }
  
  parameters <- c(b = b, Ka = Ka, Cl = Cl, V = V)
  state      <- c(A = 319.992*b, c = c_init)
  times      <- theoph.1 %>% pull(Time)
  
  ode.soln <- ode(y = state, times = times,
             func = pkpd.two.compartment, parms = parameters,
             method = "bdf") %>%
    as.data.frame %>%
    as.tibble
  
  synth.data <- ode.soln %>%
    select(time,c) %>%
    mutate(c = rnorm(nrow(.), mean = log(c), sd = sigma) %>% exp)
  
  synth.data
}

GetKaConfIntervalsStanOptim <- function(tib) {
  
  dat <- list(Nt = nrow(theoph.1),
              ts = theoph.1 %>% pull(Time) %>% `[`(-1),
              y_init = c(theoph.1 %>% pull(DoseMass) %>% head(1),
                         tib %>% pull(c) %>% head(1)),
              y = tib %>% pull(c) %>% `[`(-1))
  
  m <- stan_model(file = "theoph_single_individual.stan")
  
  opt.object <- optimizing(model, data = dat, draws = 2000, hessian = TRUE)
  
  quantiles <- opt.object$theta_tilde[,c("theta[1]")] %>%
    quantile(probs = c(0.1,0.5,0.9))
  
  quantiles
}

GetKaConfIntervalsStanNuts <- function(tib) {
  
  dat <- list(Nt = nrow(theoph.1),
              ts = theoph.1 %>% pull(Time) %>% `[`(-1),
              y_init = c(theoph.1 %>% pull(DoseMass) %>% head(1),
                         tib %>% pull(c) %>% head(1)),
              y = tib %>% pull(c) %>% `[`(-1))
  
  m <- stan_model(file = "theoph_single_individual.stan")
  
  fit <- stan(file = "theoph_single_individual.stan",
              data = dat,
              iter = 4000, chains = 1, control = list(adapt_delta = 0.9))
  
  post.samples <- extract(fit)
  
  quantiles <- post.samples$theta[,1] %>%
    quantile(probs = c(0.1,0.5,0.9))
  
  quantiles
}

post.sample.to.sim <- post.samples %>% sample_n(9)

synth.data <- post.sample.to.sim %>% pmap(ForwardSymTheo)

asymp.estimates <- synth.data %>%
  map(GetKaConfIntervalsStanOptim) %>%
  map2(1:9, function(quantiles,sample.num) quantiles %>%
         t %>%
         as_tibble %>%
         set_names(c("Q10", "Q50", "Q90")) %>%
         mutate(sample = sample.num) %>%
         mutate(method = "Asymptotic")) %>%
  bind_rows

nuts.estimates <- synth.data %>%
  map(GetKaConfIntervalsStanNuts) %>%
  map2(1:9, function(quantiles,sample.num) quantiles %>%
         t %>%
         as_tibble %>%
         set_names(c("Q10", "Q50", "Q90")) %>%
         mutate(sample = sample.num) %>%
         mutate(method = "NUTS")) %>%
  bind_rows

#save(asymp.estimates, nuts.estimates, file = "point_estimate_comparison.Rdata")

bind_rows(asymp.estimates, nuts.estimates) %>%
  ggplot(aes(method, Q50)) +
  geom_point() +
  geom_errorbar(aes(ymin = Q10, ymax = Q90)) +
  facet_wrap( ~ sample, scales = "free") +
  geom_hline(aes(x = NULL, y = NULL, yintercept = Ka),
             color = "red",
             data = post.sample.to.sim %>% select(Ka) %>% mutate(sample = 1:9)) +
  xlab("Inference Method") +
  ylab("Ka (Absorption Rate of Drug)") +
  theme(text = element_text(size=30))


