
# incorrect y model

source('metrics.R')

#### population ####

library(sampling); library(MASS);library('xgboost'); library(spatstat.geom)

set.seed(14)

N = 1e4; M = 1e3 # population and simulation size
r1 = 0; r2 = 0
sig = matrix(c(  1, r1, r1, r1,
                r1,  1, r2, r2,
                r1, r2,  1, r2,
                r1, r2, r2,  1), 4, 4)

pop = as.data.frame(mvrnorm(n = N, mu = c(1,1,1,1), Sigma = sig))
colnames(pop) = c('x1', 'x2', 'x3', 'x4')
pop$y = 3*pop$x1 + pop$x2 - 5*pop$x3 + 0.1*pop$x4 + rnorm(N)

includP = 0.05

n.NP = 0.1*N
true.mean = mean(pop$y)

#### models ####
senario = c("Z~x1+x2+x3",    #true model
            "Z~x1+x2+x3+x4", 
            "Z~x1+x2",
            "Z~x3+x4",
            "Z~x1+x3",
            "Z~x2+x4",
            "Z~x1+x4",
            "Z~x2+x3",
            "Z~(x1+x2)^2",
            "Z~(x2+x4)^2",
            "Z~(x1+x3)^2",
            "Z~(x1+x4)^2",
            "Z~(x2+x3)^2",
            "Z~(x3+x4)^2",
            "Z~(x1+x2+x3)^2",
            "Z~(x1+x2+x3+x4)^2",
            "Z~x2+x4+x1*x3",
            "Z~x1*x4",
            "Z~x1*x3",
            "Z~x1+x4*x2",
            "Z~I(x1^2)",
            "Z~I(x2^2)",
            "Z~I(x3^2)",
            "Z~I(x4^2)",
            "Z~x2+I(x3^3)",
            "Z~x1+I(x3^2)",
            "Z~x3+I(x4^3)",
            "Z~I(x2^2)+I(x3^2)",
            "Z~x2+I(x3^2)+I(x4^2)",
            "Z~x1+I(x2^2)+I(x3^2)",
            "Z~x1",
            "Z~x2",
            "Z~x3",
            "Z~x4")

#### function for estimation ####
  
tau = 0

for (m in 1:M) {
  
  #### draw two samples ####
  # srs sample (PS)
  S = sample.int(N, size = includP*N, replace = FALSE)
  sample.P = pop[S, -5]
  n_P = nrow(sample.P)
  
  # nonprobability sample (NS)
  x = pop$x1 + pop$x2 - 0.5 * pop$x3 
  f = function(theta) sum(exp(theta + x) / (1 + exp(theta + x))) - n.NP
  theta = uniroot(f, c(-100, 100))$root  
  includNP = exp(theta + x) / (1 + exp(theta + x))
  Sstar = as.logical(UPrandomsystematic(includNP))
  sample.NP = pop[Sstar,]
  n_NP = nrow(sample.NP)
  d = N/nrow(sample.P)
 
  #### propensity estimation ####
  AB = rbind(sample.NP[,-5], sample.P)
  xs = as.matrix(AB)
  AB$Z = c(rep(1,n_NP), rep(0,n_P))
  
  # for metrics
  y_model = lm(y~x2+x4, data = sample.NP)
  y_hat = predict(y_model, newdata = AB)
  Y = mean(predict(y_model, newdata = sample.P))
  
  tb = matrix(data = NA, nrow = length(senario)+1, ncol = 9) 
  colnames(tb) = c('error', 'cal2', 'cal3', 'msb_0', 'msb_025', 'msb_05', 'msb_075', 'msb_1', 'KS')
  
  for (g in 1:length(senario)) {
  model = senario[g]
  glmO = glm(model, data = AB, family = "binomial")
  O = exp(predict(glmO))
  f18 = d/O
  p = predict(glmO, type="response")
  w = f18[1:n_NP]
  
  # metrics
  cal2 = per.cal2(z = AB$Z, w = f18, y_hat = y_hat, d = rep(d, nrow(AB)))
  cal3 = per.cal3(y = sample.NP$y, w = w, y_hat = y_hat[(n_NP+1):length(y_hat)], d = rep(d, n_P))
  msb_000 = per.msb(phi = 0   , y = sample.NP$y, y_hat = y_hat[1:n_NP], Y = Y, f = n_NP/N, w = w)
  msb_025 = per.msb(phi = 0.25, y = sample.NP$y, y_hat = y_hat[1:n_NP], Y = Y, f = n_NP/N, w = w)
  msb_050 = per.msb(phi = 0.50, y = sample.NP$y, y_hat = y_hat[1:n_NP], Y = Y, f = n_NP/N, w = w)
  msb_075 = per.msb(phi = 0.75, y = sample.NP$y, y_hat = y_hat[1:n_NP], Y = Y, f = n_NP/N, w = w)
  msb_100 = per.msb(phi = 1   , y = sample.NP$y, y_hat = y_hat[1:n_NP], Y = Y, f = n_NP/N, w = w)
  
  enp= ewcdf(sample.NP$y, w = w)
  ep = ewcdf(y_hat[(n_NP+1):length(y_hat)], w = rep(d, n_P))
  yt = c(sample.NP$y, y_hat[(n_NP+1):length(y_hat)])
  KS = max(abs(enp(yt)-ep(yt)))
  
  #### performance ####
  est = sum(sample.NP$y*w/sum(w))
  error =  est - true.mean
  
  tb[g,] = c(error, cal2, cal3, msb_000, msb_025, msb_050, msb_075, msb_100, KS)
  }
  
  # boosting
  bst = xgboost(data = xs, label = AB$Z, nrounds = 1e2, 
                objective = "binary:logistic",verbose = 0) 
  p = predict(bst, xs)
  O = p/(1-p)
  f18 = d/O
  w = f18[1:n_NP]

  cal2 = per.cal2(z = AB$Z, w = f18, y_hat = y_hat, d = rep(d, nrow(AB)))
  cal3 = per.cal3(y = sample.NP$y, w = w, y_hat = y_hat[(n_NP+1):length(y_hat)], d = rep(d, n_P))
  msb_000 = per.msb(phi = 0     , y = sample.NP$y, y_hat = y_hat[1:n_NP], Y = Y, f = n_NP/N, w = w)
  msb_025 = per.msb(phi = 0.25  , y = sample.NP$y, y_hat = y_hat[1:n_NP], Y = Y, f = n_NP/N, w = w)
  msb_050 = per.msb(phi = 0.50  , y = sample.NP$y, y_hat = y_hat[1:n_NP], Y = Y, f = n_NP/N, w = w)
  msb_075 = per.msb(phi = 0.75  , y = sample.NP$y, y_hat = y_hat[1:n_NP], Y = Y, f = n_NP/N, w = w)
  msb_100 = per.msb(phi = 1     , y = sample.NP$y, y_hat = y_hat[1:n_NP], Y = Y, f = n_NP/N, w = w)
  
  enp= ewcdf(sample.NP$y, w = w)
  ep = ewcdf(y_hat[(n_NP+1):length(y_hat)], w = rep(d, n_P))
  yt = c(sample.NP$y, y_hat[(n_NP+1):length(y_hat)])
  KS = max(abs(enp(yt)-ep(yt)))
  
  est = sum(sample.NP$y*w/sum(w))
  error =  est - true.mean
  
  tb[g+1,] =  c(error, cal2, cal3, msb_000, msb_025, msb_050, msb_075, msb_100, KS)

  tau = tau + c(cor(tb[,1], tb[,2], method = "kendall"),
                cor(tb[,1], tb[,3], method = "kendall"),
                cor(tb[,1], tb[,4], method = "kendall"),
                cor(tb[,1], tb[,5], method = "kendall"),
                cor(tb[,1], tb[,6], method = "kendall"),
                cor(tb[,1], tb[,7], method = "kendall"),
                cor(tb[,1], tb[,8], method = "kendall"),
                cor(tb[,1], tb[,9], method = "kendall")) /M
  
 }

tau
library("xtable")
xtable(rbind(c('Cal2','Cal3', 'msb0', 'msb0.25', 'msb0.5', 'msb0.75', 'msb1', 'KS'),round(tau,2)), caption = 'kendall Tau')
