
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


#### function for estimation ####
  
est.fun = function(m) {
  
  #### draw two samples ####
  # srs sample (PS)
  S = sample.int(N, size = includP*N, replace = FALSE)
  sample.P = pop[S, -5]
  n_P = nrow(sample.P)
  
  # nonprobability sample (NS)
  x = pop$x1 + 0.5*pop$x2*pop$x4
  f = function(theta) sum(exp(theta + x) / (1 + exp(theta + x))) - n.NP
  theta = uniroot(f, c(-100, 100))$root  
  includNP = exp(theta + x) / (1 + exp(theta + x))
  Sstar = as.logical(UPrandomsystematic(includNP))
  sample.NP = pop[Sstar,]
  n_NP = nrow(sample.NP)
  d = N/nrow(sample.P)
 
  #### propensity estimation ####
  AB = rbind(sample.NP[,-5], sample.P)
  AB$Z = c(rep(1,n_NP), rep(0,nrow(sample.P)))
  glmO = glm(model, data = AB, family = "binomial")
  O = exp(predict(glmO))
  f18 = d/O
  p = predict(glmO, type="response")
  w = f18[1:n_NP]
  
  # for metrics
  y_model = lm(y~x2+x4, data = sample.NP)
  y_hat = predict(y_model, newdata = AB)
  Y = mean(predict(y_model, newdata = sample.P))
  
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
  true.mean = mean(pop$y)
  est = sum(sample.NP$y*w/sum(w))
  bias =  est - true.mean
  mse =  (est - true.mean)^2
  
 
  return(c(bias, mse, cal2, cal3 , msb_000, msb_025, msb_050, msb_075, msb_100, KS))
}

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

tt = matrix(data = NA, nrow = length(senario)+1, ncol = 10) 

result = list()

for (g in 1:length(senario)) {
  model = senario[g]
  result[[g]] = matrix(data = NA, nrow = M, ncol = 10)
  for (m in 1:M) {result[[g]][m,] = est.fun(m)}
  tt[g,] = colMeans(result[[g]])
}



est.fun_xgboost = function(m) {
  
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
  AB$Z = c(rep(1,n_NP), rep(0,nrow(sample.P)))
  bst = xgboost(data = xs, label = AB$Z, nrounds = 1e2, 
                objective = "binary:logistic",verbose = 0) 
  p = predict(bst, xs)
  O = p/(1-p)
  f18 = d/O
  w = f18[1:n_NP]
  
  # for metrics
  y_model = lm(y~x2+x4, data = sample.NP)
  y_hat = predict(y_model, newdata = AB)
  Y = mean(predict(y_model, newdata = sample.P))
  
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
  true.mean = mean(pop$y)
  est = sum(sample.NP$y*w/sum(w))
  bias =  est - true.mean
  mse =  (est - true.mean)^2
  
  
  return(c(bias, mse, cal2, cal3, msb_000, msb_025, msb_050, msb_075, msb_100, KS))
}

g = length(senario)+1
result[[g]] = matrix(data = NA, nrow = M, ncol = 10)
for (m in 1:M) {result[[g]][m,] = est.fun_xgboost(m)}
tt[g,] = colMeans(result[[g]])

tt = as.data.frame(tt)
tt = cbind(c(senario, "boost"), tt)
colnames(tt) = c('model', 'bias', 'mse', 'cal2', 'cal3', 'msb_000', 'msb_025', 'msb_050', 'msb_075', 'msb_100', 'KS')


pdf("sim1-2_nl.pdf", width = 14)
par(mfrow=c(2,4))

plot(tt$mse[-35], tt$cal2[-35], pch = 20, xlab = 'MSE', ylab = 'Cal2', cex = 1.5, cex.lab=1.5, cex.axis=1.5)
#text(paste("tau =", round(cor(tt$mse, tt$cal2, method = "kendall"), 2)), x = 3, y = 0.25, col = 'blue', cex = 1.5)
points(tt$mse[35], tt$cal2[35], pch=4, cex = 1.5)

plot(tt$mse[-35], tt$cal3[-35], pch = 20, xlab = 'MSE', ylab = 'Cal3', cex = 1.5, cex.lab=1.5, cex.axis=1.5)
#text(paste("tau =", round(cor(tt$mse, tt$cal3, method = "kendall"), 2)), x = 3, y = -1, col = 'blue', cex = 1.5)
points(tt$mse[35], tt$cal3[35], pch=4, cex = 1.5)

plot(tt$mse[-35], tt$msb_000[-35], pch = 20, xlab = 'MSE', ylab = 'msb_0', cex = 1.5, cex.lab=1.5, cex.axis=1.5)
#text(paste("tau =", round(cor(tt$mse, tt$msb_000, method = "kendall"), 2)), x = 3, y = -1, col = 'blue', cex = 1.5)
points(tt$mse[35], tt$msb_000[35], pch=4, cex = 1.5)

plot(tt$mse[-35], tt$msb_025[-35], pch = 20, xlab = 'MSE', ylab = 'msb_025', cex = 1.5, cex.lab=1.5, cex.axis=1.5)
#text(paste("tau =", round(cor(tt$mse, tt$msb_025, method = "kendall"), 2)), x = 3, y = 0, col = 'blue', cex = 1.5)
points(tt$mse[35], tt$msb_025[35], pch=4, cex = 1.5)

plot(tt$mse[-35], tt$msb_050[-35], pch = 20, xlab = 'MSE', ylab = 'msb_05', cex = 1.5, cex.lab=1.5, cex.axis=1.5)
#text(paste("tau =", round(cor(tt$mse, tt$msb_050, method = "kendall"), 2)), x = 3, y = 2, col = 'blue', cex = 1.5)
points(tt$mse[35], tt$msb_050[35], pch=4, cex = 1.5)

plot(tt$mse[-35], tt$msb_075[-35], pch = 20, xlab = 'MSE', ylab = 'msb_075', cex = 1.5, cex.lab=1.5, cex.axis=1.5,ylim=c(2,7))
#text(paste("tau =", round(cor(tt$mse, tt$msb_075, method = "kendall"), 2)), x = 3, y = 7, col = 'blue', cex = 1.5)
points(tt$mse[35], tt$msb_075[35], pch=4, cex = 1.5)

plot(tt$mse[-35], tt$msb_100[-35], pch = 20, xlab = 'MSE', ylab = 'msb_1', cex = 1.5, cex.lab=1.5, cex.axis=1.5,ylim=c(14.5,40))
#text(paste("tau =", round(cor(tt$mse, tt$msb_100, method = "kendall"), 2)), x = 3, y = 40, col = 'blue', cex = 1.5)
points(tt$mse[35], tt$msb_100[35], pch=4, cex = 1.5)

plot(tt$mse[-35], tt$KS[-35], pch = 20, xlab = 'MSE', ylab = 'KS', cex = 1.5, cex.lab=1.5, cex.axis=1.5)
#text(paste("tau =", round(cor(tt$mse, tt$KS, method = "kendall"), 2)), x = 18, y = 0.65, col = 'blue', cex = 1.5)
points(tt$mse[35], tt$KS[35], pch=4, cex = 1.5)

dev.off()

