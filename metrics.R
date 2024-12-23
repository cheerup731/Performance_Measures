
#### function for performance metrics ####

## Calibration ##

per.cal1 = function(z, w, X, d){ 
  # z = inclusion of NP in AB set, w = pseudo weight, X = auxiliaries, d = design weight
  X = as.matrix(X)
  sum(abs((z*w/sum(z*w)-(1-z)*d/sum((1-z)*d))%*%X))
}

per.cal2 = function(z, w, y_hat, d){ 
  # y_hat = estimated y in B set
  (z*w/sum(z*w)-(1-z)*d/sum((1-z)*d))%*%y_hat
}

per.cal3 = function(y, w, y_hat, d){ 
  # y = observed y in NP, y_hat = estimated y in P
  sum(w*y)/sum(w)-(sum(d*y_hat)/sum(d))
}

## Measure of Selection Bias

per.msb = function(phi, y, y_hat, Y, f, w){
  # y = observed y in NP, y_hat = estimated y in NP (proxy) 
  # Y = estimated population mean of y_hat by P
  r_my = cor(y, y_hat)
  MUB = ((phi+(1-phi)*r_my)/(phi*r_my+1-phi))*sqrt(var(y)/var(y_hat))*(mean(y_hat)-Y)
  MUB-mean(y)+sum(w*y)/sum(w)
  
}

## comparing to 

per.entropy = function(truth, ppen, eps = 1e-6){
  ppen[ppen < eps] = eps
  ppen[ppen > 1-eps] = 1-eps
  -sum(truth*log(ppen)+(1-truth)*log(1-ppen))/length(truth)
}

per.brier = function(truth, ppen){
  sqrt(sum((ppen - truth)^2)/length(truth))
}


