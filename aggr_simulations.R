# Lost in Aggregation: The Causal Interpretation of the IV Estimand
# code for simulations + figures 

library(mvtnorm)
library(ggplot2)
library(dplyr)
library(tidyr)
library(knitr)
library(ivreg)
library(latex2exp)
library(patchwork)
library(gridExtra)
library(car)
library(gmm)

set.seed(0)

# figure 2 code -----------------------------------------------------------

# simulate data from SCM (1) with k=2
sim_scm1 <- function(n, alpha, gamma, delta, beta, gamma_y, variances) {
  # pile all constants to construct A1, A2 as a vector
  ab <- cbind(delta, gamma, c(1,0), c(0,1)) 
  
  # unpack variances
  varI <- variances$I
  varU <- variances$U
  var_eps <- variances$eps_A12
  var_eps_y <- variances$eps_Y
  
  # generate I, U, and epsilons from joint multivariate normal distn
  I.U.eps <- rmvnorm(n, rep(0,4), sigma=diag(x=c(varI, varU, 
                                                 rep(var_eps,2)), ncol=4))
  colnames(I.U.eps) <- c("I", "U", "eps1", "eps2")
  U <- I.U.eps[,"U"] 
  I <- I.U.eps[,"I"] 
  A12 <- I.U.eps %*% t(ab) 
  A <- A12%*%(alpha) # compute aggregate A
  Y <- A12%*%(beta) + gamma_y*U + rnorm(n,mean=0,sd=sqrt(var_eps_y))
  
  #--- ivreg --- 
  ivreg.dat <- data.frame(Y=Y, A=A, I=I)
  tsls_model <- ivreg(Y ~ A | I, data = ivreg.dat)
  beta.tsls <- coef(tsls_model)["A"]
  se.tsls <- sqrt(diag(vcov(tsls_model)))["A"] # stored
  
  return(list(beta.tsls = beta.tsls, se.tsls = se.tsls))
}

#-------- SIM A: prop aggr divergence --------
# set aggregation constants
alpha <- matrix(c(1,1), ncol=1)
beta2 <- 2
delta <- matrix(c(1,1), ncol=1)
gamma <- matrix(c(1,1),ncol=1)
gamma_y <- 1

# setting variances for variables
variances <- list(I = 1, U = 1, eps_A12 = 1, eps_Y = 1)

# ACID mean
d1 <- 2
d <- matrix(c(d1,(1-alpha[1]*d1)/alpha[2]), ncol=1)

# calculate IV (2SLS) estimates
beta1.range <- (-10:40)/10
n.range <- c(10, 100, 1e3)
tsls.matrix <- matrix(data=NA, nrow=length(n.range), ncol=length(beta1.range))
se.matrix <- matrix(data=NA, nrow=length(n.range), ncol=length(beta1.range))
colnames(tsls.matrix)  <- colnames(se.matrix)  <- beta1.range
rownames(tsls.matrix) <- rownames(se.matrix) <- n.range
  
for (nn in 1:length(n.range)) {
  for (beta1 in 1:length(beta1.range)) {
    beta <- matrix(c(beta1.range[beta1],beta2),ncol=1)
    cur <- sim_scm1(n.range[nn], alpha, gamma, delta, beta, gamma_y, variances)
    tsls.matrix[nn, beta1] <- cur$beta.tsls
    se.matrix[nn, beta1] <- cur$se.tsls
  }
}

# SIM A results matrix
tsls.df <- tsls.matrix %>%
  as.data.frame() %>%
  mutate(n = as.numeric(rownames(.))) %>%
  pivot_longer(cols = -n, names_to = "beta1", values_to = "val") %>%
  mutate(beta1 = as.numeric(beta1))

#---- SIM B: instr-tuned intervention divergence ----

# set aggregation constants
alpha <- matrix(c(1,1), ncol=1)
beta <- matrix(c(1,2),ncol=1)
delta <- matrix(c(1,1), ncol=1)
gamma <- matrix(c(1,1),ncol=1)
gamma_y <- 1

# setting variances for variables
variances <- list(I = 1, U = 1, eps_A12 = 1, eps_Y = 1)

# ACID conditional mean
d1.range <- (-20:20)/10
n.range <- c(10, 100, 1e3)

tsls.matrix2 <- matrix(data=NA, nrow=length(n.range), ncol=length(d1.range))
colnames(tsls.matrix2) <- d1.range
rownames(tsls.matrix2) <- n.range
  
for (nn in 1:length(n.range)) {
  for (dd in 1:length(d1.range)) {
    d1 <- d1.range[dd]
    d <- matrix(c(d1,(1-alpha[1]*d1)/alpha[2]), ncol=1)
    cur <- sim_scm1(n.range[nn], alpha, gamma, delta, beta, gamma_y, variances)
    tsls.matrix2[nn, dd] <- cur$beta.tsls
  }
}

# SIM B results matrix
tsls.df2 <- tsls.matrix2 %>%
  as.data.frame() %>%
  mutate(n = as.numeric(rownames(.))) %>%
  pivot_longer(cols = -n, names_to = "d1", values_to = "val") %>%
  mutate(d1 = as.numeric(d1))

#-------- plot --------
all_linetypes <- c("IV estimand" = "solid", 
                   "prop aggr" = "dashed",
                   "ITI" = "dotdash",
                   "true ACE(A,Y)" = "solid")
all_shapes <- c("10" = 20, "100" = 17, "1000" = 18) 
all_labels <- c("prop aggr", "ITI", "IV estimand", "true ACE(A,Y)")

p1 <- ggplot(tsls.df, aes(x = beta1, y = val, color = factor(n), alpha = factor(n))) +
  geom_line(aes(y = (beta1*d[1] + beta[2]*(1-alpha[1]*d[1])/alpha[2]), 
                linetype = "true ACE(A,Y)"),
            color="grey", alpha=0.7, linewidth = 1) + 
  geom_line(aes(y = (beta1*delta[1] + beta[2]*delta[2])/
                          (alpha[1]*delta[1] + alpha[2]*delta[2]),
                linetype = "IV estimand"), 
            color = "goldenrod", linewidth = 2, alpha=0.7) + 
  geom_vline(aes(xintercept = 2,
                 linetype = "prop aggr"),
            color = "violet", linewidth = 1, alpha=0.7) +
  geom_point(aes(shape=factor(n)), size=2) +
  labs(x = TeX(r'($\beta_1$)'), y ="aggr causal effect", color = "n", 
       title="A: proportional aggregation") +
  theme_minimal() +
  scale_color_discrete(name ="sample size n") +
  scale_shape_manual(name = "sample size n", values = all_shapes) +
  scale_alpha_manual(values = c("10"=0.8, "100"=0.8, "1000"=0.8), guide="none") +
  scale_linetype_manual(
    name = "",
    values = all_linetypes,
    breaks = all_labels, 
    drop=FALSE
  ) +
  theme(legend.position = "none")

p2 <- ggplot(tsls.df2, aes(x = d1, y = val, color = factor(n), alpha=factor(n))) +
  geom_line(aes(y = (beta[1]*d1 + beta[2]*(1-alpha[1]*d1)/alpha[2]), linetype = "true ACE(A,Y)"),
            color="grey", alpha=0.7, linewidth = 1) + 
  geom_line(aes(y = (beta[1]*delta[1] + beta[2]*delta[2])/
                          (alpha[1]*delta[1] + alpha[2]*delta[2]),
                linetype = "IV estimand"), 
            color = "goldenrod", linewidth = 2, alpha=0.7) + 
  geom_vline(aes(xintercept = delta[1]/(delta[1]*alpha[1] + delta[2]*alpha[2]), 
                 linetype = "ITI"), 
            color = "coral", linewidth = 1, alpha=0.7) + 
  geom_vline(aes(xintercept = 2.1, 
                 linetype = "prop aggr"), 
             color = "violet", linewidth = 1, alpha=0.7) + 
  # ^just a placeholder to make the legend work, the warning can be ignored
  geom_point(aes(shape=factor(n)), size=2) +
  xlim(c(-2,2)) + 
  labs(x = TeX(r'($d_1$)'), y = "", color = "n", title="B: instrument-tuned intervention") +
  theme_minimal() +
  scale_color_discrete(name = "sample size n") +
  scale_shape_manual(name = "sample size n", values = all_shapes) +
  scale_alpha_manual(values = c("10"=0.8, "100"=0.8, "1000"=0.8), guide="none") +
  scale_linetype_manual(
    name = "",
    values = all_linetypes,
    breaks=all_labels, 
    drop=FALSE
  ) +
  theme(legend.position = "none")

combined <- p1 + p2 
combined
ggsave("fig/fig2_divergence.pdf", width=8, height=4)


# figure 4 code -----------------------------------------------------------

set.seed(0)
n <- 1e3

#-------- functions --------
# g function for Sargan
g_pooled <- function(theta, x) {
  i1 <- x[,1]
  i2 <- x[,2]
  y <- x[,3]
  a <- x[,4]
  residual <- y - theta*a
  return(cbind(i1*residual, i2*residual))
}

# simulate data from SCM (25) with k=2
sim_scm25_data <- function(beta1, n, delta1, delta2) {
  # delta_i: coefficients for I_i
  i1 <- rnorm(n)
  i2 <- rnorm(n)
  d11 <- delta1[1]
  d21 <- delta1[2]
  d12 <- delta2[1]
  d22 <- delta2[2]
  
  alpha1 <- 1
  alpha2 <- 1
  beta2 <- 2 
  u <- rnorm(n)
  a1 <- d11*i1 + d12*i2 + 0.5*u + rnorm(n)
  a2 <- d21*i1 + d22*i2 + 0.5*u + rnorm(n)
  a <- alpha1*a1 + alpha2*a2
  y <- beta1*a1 + beta2*a2 + 2*u + rnorm(n)
  
  return(data.frame(i1=i1, i2=i2, y=y, a=a))
}

# run Sargan; return test statistics and p-values
run_test <- function(i1, i2, y, a) {
  result_pooled <- gmm(g_pooled, x=cbind(i1, i2, y, a), t0=1, optfct="optim", method="BFGS")
  sargan <- summary(result_pooled)$stest$test 
  
  # test for weak instruments
  fs_full <- lm(a ~ i1 + i2)
  fs_restricted <- lm(a ~ 1)
  anova_result <- anova(fs_restricted, fs_full)
  F_stat <- anova_result$F[2]
  
  return(data.frame(test_stat_sargan=as.numeric(sargan[,"J-test"]),
                    pval_sargan=as.numeric(sargan[,"P-value"]),
                    F_stat=F_stat))
}

#-------- run sim --------
# beta_1 values
beta_seq <- (-10:40)/10

# instrument strengths
sw <- c(5,3,0.1,0.2)
ss <- c(5,3,4,2)
ww <- c(.15,.1,.08,.05)
instr_coefs <- rbind(sw, ss, ww)
rownames(instr_coefs) <- c("sw","ss","ww")
colnames(instr_coefs) <- c("d1_1", "d1_2", "d2_1", "d2_2")

# table 1
instr_coefs # columns 1-2
denom <- sqrt((instr_coefs[,"d1_1"] + instr_coefs[,"d1_2"])^2 + 
                (instr_coefs[,"d2_1"] + instr_coefs[,"d2_2"])^2 + 3)
table1 <- cbind(instr_coefs, (instr_coefs[,1]+instr_coefs[,2])/denom, 
      (instr_coefs[,3]+instr_coefs[,4])/denom)
colnames(table1) <- c(colnames(instr_coefs), "cor(I1,A)","cor(I2,A)")
table1

all_results <- list(sw=list(), ss=list(), ww=list())

for (strength_combn in rownames(instr_coefs)) {
  delta1 <- instr_coefs[strength_combn,c("d1_1", "d1_2")]
  delta2 <- instr_coefs[strength_combn,c("d2_1", "d2_2")]
  
  for (beta in beta_seq) {
    print(beta)
    for (rep in 1:100) {
      cur_dat <- sim_scm25_data(beta1=beta,n=n, delta1=delta1, delta2=delta2)
      cur_test <- run_test(cur_dat$i1, cur_dat$i2, cur_dat$y, cur_dat$a)
      
      result <- data.frame(beta1 = beta,
                           sim = rep,
                           stat_sargan = cur_test$test_stat_sargan,
                           pval_sargan = cur_test$pval_sargan,
                           cor_IA1 = cor(cur_dat$i1, cur_dat$a),
                           cor_IA2 = cor(cur_dat$i2, cur_dat$a),
                           F_stat = cur_test$F_stat)
      all_results[[strength_combn]] <- c(all_results[[strength_combn]], list(result))
    }
  }
}

combined_results <- bind_rows(
  bind_rows(all_results$sw) %>% mutate(strength = "strong-weak"),
  bind_rows(all_results$ss) %>% mutate(strength = "strong-strong"),
  bind_rows(all_results$ww) %>% mutate(strength = "weak-weak")
)

plot_data <- combined_results %>%
  group_by(beta1, strength) %>%
  summarise(
    pct_0.01 = mean(pval_sargan<0.01)*100,
    pct_0.5 = mean(pval_sargan<0.5)*100,
    median_cor_IA1 = median(cor_IA1),
    median_cor_IA2 = median(cor_IA2),
    pct_F_below_11 = mean(F_stat < 11)*100,
    .groups = 'drop'
  ) 

#--- plot % p-values < 0.01
p1 <- ggplot(plot_data, aes(x = beta1, y = pct_0.01, 
                            color = strength, group = strength)) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = 2, linetype = "dashed", 
             color = "violet", linewidth = 1) + 
  geom_point() +
  labs(title = "A: type-I-error 0.01",
       x = TeX(r'($\beta_1$)'),
       y = "rejection frequency",
       color = "instr strength") +
  theme_minimal() +
  ylim(c(0,100)) + 
  theme(legend.position = "none")

#--- plot % p-values < 0.5
p2 <- ggplot(plot_data, aes(x = beta1, y = pct_0.5, 
                            color = strength, group = strength)) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = 2, linetype = "dashed", 
             color = "violet", linewidth = 1) + 
  geom_point() +
  labs(title = "B: type-I-error 0.5",
       x = TeX(r'($\beta_1$)'),
       y = "",
       color = "instr strength") +
  theme_minimal() +
  ylim(c(0,100))+
  theme(legend.position = "none")

combined <- p1 + p2
combined
ggsave("fig/fig4_sargan_power.pdf", width=8, height=4)


# fig 5 -------------------------------------------------------------------
vline_data <- data.frame(strength="weak-weak", xint = 11) # rule-of-thumb cut-off

ggplot(combined_results, aes(x = F_stat, fill = strength)) +
  geom_histogram(bins=50) +
  geom_vline(data=vline_data, aes(xintercept=xint),
             color = "goldenrod", linewidth = 0.8) +
  facet_wrap(~strength, scales="free") +
  labs(title = "", x = "F-statistic",
       y = "count") +
  theme_minimal() +
  theme(legend.position = "none")

ggsave("fig/fig5_F_statistics.pdf", width=8, height=4)

