
# applications

source("metrics.R")
library('xgboost')
library("sampling"); library(spatstat.geom)
set.seed(14)

data(MU284)

#  data from Särndal et al (1992) in sampling package
#  CL as the inclusion of NP
#  P85 as the target variable, which is an extremely skew variable

pop = MU284[,-1] # remove labels 

colnames(pop) = c('y', 'x1', 'x2', 'x3', 'x4', 'x5', 'x6', 'x7', 'x8', 'S_star')
true.mean = mean(pop$y)
N = nrow(pop)

# draw a srs sample
d = 10
S = sample.int(N, size = round(N/d), replace = FALSE)
sample.P = pop[S, -ncol(pop)]
sample.P$Z = 0
n_P = nrow(sample.P)

# non-prob by a categorical data
sample.NP = pop[pop$S_star < 15, -ncol(pop)]
sample.NP$Z = 1
n_NP = nrow(sample.NP)

AB = rbind(sample.NP, sample.P)
AB = AB[-1]
xs = as.matrix(AB[,-ncol(AB)])

tuna = expand.grid(eta = c(0, 0.3, 0.5, 0.7),    
                   max_depth = 6,
                   min_child_weight = c(1, 3, 5, 7),
                   gamma = c(0, 1, 2))

tt = matrix(data = NA, nrow = nrow(tuna), ncol = 10) 

for (i in 1:nrow(tuna)) {
  
par = list(eta = tuna[i,1],     
           max_depth = tuna[i,2],
           min_child_weight = tuna[i,3],
           gamma = tuna[i,4])
           

bst = xgboost(data = xs, label = AB$Z, params = par, nrounds = 1e2, 
              objective = "binary:logistic",verbose = 0) 
              
p = predict(bst, xs)
O = p/(1-p)
f18 = d/O
w = f18[1:n_NP]

# for metrics
y_model = lm(y~., data = sample.NP)
y_hat = predict(y_model, newdata = AB)
Y = mean(predict(y_model, newdata = sample.P))

# metrics
cal1 = per.cal1(z = AB$Z, w = f18, X = xs, d = rep(d, nrow(AB)) )
cal2 = per.cal2(z = AB$Z, w = f18, y_hat = y_hat, d = rep(d, nrow(AB)))
cal3 = per.cal3(y = sample.NP$y, w = w, y_hat = y_hat[(n_NP+1):length(y_hat)], d = rep(d, n_P))
msb_00 = per.msb(phi = 0  , y = sample.NP$y, y_hat = y_hat[1:n_NP], Y = Y, f = n_NP/N, w = w)
msb_05 = per.msb(phi = 0.5, y = sample.NP$y, y_hat = y_hat[1:n_NP], Y = Y, f = n_NP/N, w = w)
msb_10 = per.msb(phi = 1  , y = sample.NP$y, y_hat = y_hat[1:n_NP], Y = Y, f = n_NP/N, w = w)
entropy = per.entropy(truth = AB$Z, ppen = p)
brier = per.brier(truth = AB$Z, ppen = p)

enp= ewcdf(sample.NP$y, w = w)
ep = ewcdf(y_hat[(n_NP+1):length(y_hat)], w = rep(d, n_P))
yt = c(sample.NP$y,y_hat[(n_NP+1):length(y_hat)])
KS = max(abs(enp(yt)-ep(yt)))

est = sum(sample.NP$y*w/sum(w))
error =  (est - true.mean)

tt[i,] = c(error, cal1, cal2, cal3, msb_00, msb_05, msb_10, entropy, brier, KS)

}

tt = as.data.frame(tt)
colnames(tt) = c('error', 'cal1', 'cal2', 'cal3', 'msb_00', 'msb_05', 'msb_10', 'entropy', 'brier', 'KS')


pdf("MU284.pdf", width = 12)
par(mfrow=c(2,4))

plot(tt$`error`, tt$cal1, pch = 20, xlab = 'error', ylab = 'Cal1', cex = 1.5, cex.lab=1.5, cex.axis=1.5)
#text(paste("tau =", round(cor(tt$`error`, tt$cal1, method = "kendall"), 2)),x = 4, y = 2500, col = 'blue', cex = 1.5)

plot(tt$`error`, tt$cal2, pch = 20, xlab = 'error', ylab = 'Cal2', cex = 1.5, cex.lab=1.5, cex.axis=1.5)
#text(paste("tau =", round(cor(tt$`error`, tt$cal2, method = "kendall"), 2)), x = 4, y = 15, col = 'blue', cex = 1.5)

plot(tt$`error`, tt$cal3, pch = 20, xlab = 'error', ylab = 'Cal3', cex = 1.5, cex.lab=1.5, cex.axis=1.5)
#text(paste("tau =", round(cor(tt$`error`, tt$cal3, method = "kendall"), 2)), x = 4, y = 15, col = 'blue', cex = 1.5)

plot(tt$`error`, tt$msb_00, pch = 20, xlab = 'error', ylab = 'msb_0', cex = 1.5, cex.lab=1.5, cex.axis=1.5)
#text(paste("tau =", round(cor(tt$`error`, tt$msb_00, method = "kendall"), 2)), x = 4, y = 15, col = 'blue', cex = 1.5)

# plot(tt$`error`, tt$msb_05, pch = 20, xlab = 'error', ylab = 'msb_05', cex = 1.5, cex.lab=1.5, cex.axis=1.5)
# abline(lm(tt$msb_05~tt$`error`), col = 'blue', lwd = 1)
# text(paste("tau =", round(cor(tt$`error`, tt$msb_05, method = "kendall"), 2)), x = 4, y = 15, col = 'blue')

plot(tt$`error`, tt$msb_10, pch = 20, xlab = 'error', ylab = 'msb_1', cex = 1.5, cex.lab=1.5, cex.axis=1.5)
#text(paste("tau =", round(cor(tt$`error`, tt$msb_10, method = "kendall"), 2)), x = 4, y = 15, col = 'blue', cex = 1.5)

plot(tt$`error`, tt$entropy, pch = 20, xlab = 'error', ylab = 'MXE', cex = 1.5, cex.lab=1.5, cex.axis=1.5)
#text(paste("tau =", round(cor(tt$`error`, tt$entropy, method = "kendall"), 2)), x = 4, y = 0.60, col = 'blue', cex = 1.5)

plot(tt$`error`, tt$brier, pch = 20, xlab = 'error', ylab = 'brier', cex = 1.5, cex.lab=1.5, cex.axis=1.5)
#text(paste("tau =", round(cor(tt$`error`, tt$brier, method = "kendall"), 2)), x = 4, y = 0.45, col = 'blue', cex = 1.5)

plot(tt$`error`, tt$KS, pch = 20, xlab = 'error', ylab = 'KS', cex = 1.5, cex.lab=1.5, cex.axis=1.5)
#text(paste("tau =", round(cor(tt$`error`, tt$KS, method = "kendall"), 2)), x = 4, y = 0.18, col = 'blue', cex = 1.5)

dev.off()
