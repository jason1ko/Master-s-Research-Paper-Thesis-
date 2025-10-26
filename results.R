

##########################    #################

rm(list = ls())

##############       3.1      plot of Phi(Theta0) when rho = 0.

par(mfrow = c(1,1))



theta0_rho0 <- function(g) { # only alpha is added
  th <- numeric(length(g))  # Initialize vector to store results
  th[g > 0 & g < 1] <- alpha^2 / sqrt(alpha^2 / (1 - g[g > 0 & g < 1]) + g[g > 0 & g < 1] / (1 - g[g > 0 & g < 1]))
  th[g > 1] <- (alpha^2 / g[g > 1]) / sqrt(alpha^2 / (g[g > 1] * (g[g > 1] - 1)) + 1 / (g[g > 1] - 1))
  th[g == 1] <- 0  # Handle g == 1
  th[g <= 0] <- NA
  
  return(pnorm(-th))
}


dev.new()

par(mar = c(5,6,1,1))

alpha <- 3
curve(theta0_rho0, 0, 10, col = 'orange', lwd = 2,
      xlab = expression(gamma), ylab = expression(Phi(-Theta(0))), cex.lab = 2)
abline(v = 1, lty = 'dashed')

alphas <- c(1,2)
colorvec <- c('darkgreen', 'purple')
for(j in 1:length(alphas)){
  alpha <- alphas[j]
  color <- colorvec[j]
  curve(theta0_rho0, 0, 10, col = color, lwd = 2, add = TRUE)
}
legend("bottomright", legend = c(1,2,3),
       title = expression(alpha), cex = 1,
       col=c("darkgreen", "purple", "orange"), pch=20, bg = "white")
#legend = expression(alpha == 1, alpha == 2, alpha == 3)
# legend = c(1,2,3)
text(1, par("usr")[3], labels = "x=1", pos = 1, xpd = TRUE, col = 'red')






####################     4.3.1      check the overview of theta
# for different rho, gamma


library(MASS)

library(Rcpp)
#install.packages("RcppEigen")
#library(RcppEigen)
library(RcppArmadillo)
library(RcppDist)

#sourceCpp("rda_rcpp.cpp")
sourceCpp("rda_rcpp2.cpp")




get.stieltjes.func = function(X, Y) {
  
  M = ((sum(Y == 0) - 1) * var(X[Y == 0,])
       + (sum(Y == 1) - 1) * var(X[Y == 1,])) /
    (length(Y) - 2)
  #xx = rev(eigen_cpp(M)$values)
  xx = svd(M, nu = 0, nv = 0)$d
  gamma = ncol(X) / nrow(X)
  
  
  # defs of m(-lambda), v(-lambda),...
  m = Vectorize(function(lambda) {
    mean(1 / (xx + lambda))
  })
  mp = Vectorize(function(lambda) {
    mean(1 / (xx + lambda)^2)
  })
  v = Vectorize(function(lambda) {
    (m(lambda) - 1/lambda) * gamma + 1/lambda
  })
  vp = Vectorize(function(lambda) {
    (mp(lambda) - 1/lambda^2) * gamma + 1/lambda^2
  })
  
  tau = Vectorize(function(lambda) {
    lambda * m(lambda) * v(lambda)
  })
  eta = Vectorize(function(lambda) {
    (v(lambda) - lambda * vp(lambda)) / gamma
  })
  xi = Vectorize(function(lambda) {
    vp(lambda) / v(lambda)^2 - 1
  })
  
  return(list(m=m, v=v, vp=vp, tau=tau, eta=eta, xi=xi))
}







get.theta.ar1 <- function(n,p,rho) {
  alpha <- 2
  Y = c(rep(0, n/2), rep(1, n/2))
  
  mu = replicate(2, alpha * rnorm(p) / sqrt(p)) # p by 2 matrix
  
  Sigma = outer(1:p, 1:p, function(x, y) rho^abs(x - y))
  #Z = mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
  Z = rmvnorm_rcpp(n = n, mean = matrix(0, nrow = 1,ncol = p), sigma = Sigma)
  X = t(mu)[Y + 1,] + Z # refer to page 15 or 20
  
  stj = get.stieltjes.func(X, Y)
  
  theta = function(lambda) {
    alpha^2 * stj$tau(lambda) / sqrt(alpha^2 * stj$eta(lambda) + stj$xi(lambda))
  }
  return(theta)
}






draw_plots_ar1 <- function(n, r){ # new trial
  gammas <- c(0.5, 1, 3)
  alpha <- 2
  rho <- r
  
  if(rho == 0.1) lmax <- 1000 # 500 for alpha = 1
  if(rho == 0.3) lmax <- 20 # 20 for alpha = 1
  if(rho == 0.5) lmax <- 10  # 5 for alpha = 1
  if(rho == 0.7) lmax <- 2 # 1 for alpha = 1
  if(rho == 0.9) lmax <- 0.5 # 0.2 for alpha = 1
  
  colors <- 1:length(gammas)
  

  for(i in colors){
    gamma <- gammas[i]
    theta.s <- get.theta.ar1(n, n*gamma, rho)
    
    max_th <- max(theta.s(seq(0, lmax, length.out = 101)[-1] ) )
    min_th <- min(theta.s(seq(0, lmax, length.out = 101)[-1] ) )
    
    plot(function(ll)theta.s(ll), 0, lmax, lwd = 1, col = 'brown',
         #xlab = expression(lambda), ylab = expression(Theta), cex.lab = 2,
         xlab = "", ylab = "",
         ylim = c(min_th*0.98,max_th*1.02))
    title(main = bquote(gamma ~ "=" ~ .(gammas[i]) ~ "," ~ rho == .(rho) ),
          cex.main = 3)
    
    for(j in 1:19){
      theta.s <- get.theta.ar1(n, n*gamma, rho)
      plot(function(ll)theta.s(ll), 0, lmax, lwd = 1, col = 'brown', add = TRUE)
    }
  }
  
}


n <- 200
set.seed(2023311161)

dev.new()
#par(mfrow = c(2,3), oma = c(2,3,0,3))
par(mfrow = c(2,3), oma = c(5, 5, 3, 0), mar = c(2, 2, 4, 2)) #

# 0.7 or 0.9 are best for alpha = 2
#draw_plots_ar1(n, 0.1)
draw_plots_ar1(n, 0.3)
#draw_plots_ar1(n, 0.5)
#draw_plots_ar1(n, 0.7)
draw_plots_ar1(n, 0.9)


mtext(expression(lambda), 
      side = 1, outer = TRUE, line = 2, cex = 2) # Common x-axis label
mtext(expression(Theta), 
      side = 2, outer = TRUE, line = 2, cex = 2) # Common y-axis label

#mtext(text = bquote(~ theta ~ "vs" ~ lambda ~ "when" ~ rho ~ "=" ~ .(r) ~ ", AR(1)"), outer = TRUE, cex = 3)



####################     4.3.2   optimal lambda




lambda_generate <- function(rho){
  num <- 1000
  lmin <- 0.01
  
  if(rho == 0.0) ls <- seq(lmin, 2000, length.out = num)
  if(rho == 0.1) ls <- seq(lmin, 1000, length.out = num) # 500 for alpha = 1
  if(rho == 0.3) ls <- seq(lmin, 20, length.out = num)  # 20 for alpha = 1
  if(rho == 0.5) ls <- seq(lmin, 10, length.out = num)  # 5 for alpha = 1
  if(rho == 0.7) ls <- seq(lmin, 2, length.out = num)  # 1 for alpha = 1
  if(rho == 0.9) ls <- seq(lmin, 0.5, length.out = num) # 0.2 for alpha = 1
  
  return(ls)
}



#gamma_vec <- seq(0.01, 1, 0.01)
gamma_vec <- seq(0.05, 4, 0.05)




optimal_lambda_ar1 <- function(n, rho){
  
  ls <- lambda_generate(rho)
  #gamma_vec <- seq(0.10, 4.00, by = 0.02)
  
  opt.ls <- c()
  thetas_list <- list()
  thetas_opt <- c()
  thetas_0 <- c()
  thetas_inf <- c()
  
  ps <- as.integer(n*gamma_vec)
  
  for(p in ps){ # 20, 40, ... 600
    idx <- which(ps == p)
    theta <- get.theta.ar1(n, p, rho)
    opt.ls[idx] <- ls[which.max( theta(ls) ) ] 
    
    thetas_list[[idx]] <- theta
  }
  
  # Theta(optimal lambda) vs gamma
  for(i in 1:length(ps)){
    thetas_opt[i] <- thetas_list[[i]](opt.ls[i])
  }
  # Theta(0) vs gamma
  unt <- which(gamma_vec == 1)
  for(i in 1:length(ps)){
    if(i > unt){
      thetas_0[i] <- NA
    }else{
      thetas_0[i] <- thetas_list[[i]](min(ls)) # min(ls)
    }
  }
  # Theta(infty) vs gamma
  for(i in 1:length(ps)){
    thetas_inf[i] <- thetas_list[[i]](10^6) # 10^6
  }
  
  return(cbind(opt.ls, thetas_opt, thetas_0, thetas_inf))
  
}


#ar1_result_0 <- optimal_lambda_ar1(n, 0.0)
ar1_result_3 <- optimal_lambda_ar1(n, 0.3)
ar1_result_5 <- optimal_lambda_ar1(n, 0.5)
ar1_result_7 <- optimal_lambda_ar1(n, 0.7)
ar1_result_9 <- optimal_lambda_ar1(n, 0.9)


plots_optimal_lambda_ar1 <- function(result, rho){
  opt.ls <- result[,'opt.ls']
  #gamma_vec <- seq(0.10, 4.00, by = 0.02)
  
  # plot of optimal lambda vs gamma
  plot(gamma_vec, opt.ls, type = 'p', cex = 0.5,
       #xlab = expression(gamma), ylab = expression(Optimal ~ lambda),
       main = bquote(rho == .(rho)), cex.main = 3
  )
  lines(smooth.spline(gamma_vec, opt.ls), col = "skyblue", lwd = 2)
}




dev.new()
#par(mfrow = c(1,4), oma = c(0, 0, 5, 0))
par(mfrow = c(1, 4), oma = c(5, 5, 3, 0), mar = c(2, 2, 4, 2)) #

plots_optimal_lambda_ar1(ar1_result_3, 0.3)
plots_optimal_lambda_ar1(ar1_result_5, 0.5)
plots_optimal_lambda_ar1(ar1_result_7, 0.7)
plots_optimal_lambda_ar1(ar1_result_9, 0.9)
#mtext(text = bquote(lambda[opt] ~ "vs" ~ gamma ~ ", AR(1)"), outer = TRUE, cex = 3)


mtext(expression(gamma), 
      side = 1, outer = TRUE, line = 2, cex = 2) # Common x-axis label
mtext(expression(lambda[opt]), 
      side = 2, outer = TRUE, line = 2, cex = 2) # Common y-axis label
#par(mfrow = c(1,1), oma = c(0, 0, 0, 0))



###################    4.3.3    Theta(optimal lambda) vs Theta(inf), gamma > 0





#Thetas_0_ar1 <- ar1_result_0[,-1]
Thetas_3_ar1 <- ar1_result_3[,-1]
Thetas_5_ar1 <- ar1_result_5[,-1]
Thetas_7_ar1 <- ar1_result_7[,-1]
Thetas_9_ar1 <- ar1_result_9[,-1]

#sum(!is.na(ar1_result_3[,'thetas_0'])) # 20


plots_opt_inf_ar1 <- function(Thetas, rho){
  
  thetas_opt <- Thetas[,'thetas_opt']
  #thetas_0   <- Thetas[,'thetas_0']
  thetas_inf <- Thetas[,'thetas_inf']
  
  th_max <- max(thetas_opt, thetas_inf)
  th_min <- min(thetas_opt, thetas_inf)
  
  plot(gamma_vec, log(pnorm(-thetas_opt)), type = 'p', cex = 0.5, col = 'red',
       xlab = expression(gamma), ylab = expression(log(Phi(-Theta))), 
       ylim = c(log(pnorm(-th_max)), log(pnorm(-th_min))) , cex.lab = 1
  )
  title(main = bquote(rho == .(rho)), cex.main = 3)
  
  lines(smooth.spline(gamma_vec, log(pnorm(-thetas_opt))), col = 'red', lwd = 2)
  
  points(gamma_vec, log(pnorm(-thetas_inf)), cex = 0.5, col = 'green')
  lines(smooth.spline(gamma_vec, log(pnorm(-thetas_inf))), col = 'green', lwd = 2)
  
  legend("bottomright", legend = c(expression(lambda[opt]), expression(infinity)), 
         title = expression(lambda), col = c('red','green'), lwd = 2, cex = 1.5)
}



dev.new()
par(mfrow = c(1, 4), oma = c(5, 5, 3, 0), mar = c(2, 2, 4, 2))

plots_opt_inf_ar1(Thetas_3_ar1, 0.3)
plots_opt_inf_ar1(Thetas_5_ar1, 0.5)
plots_opt_inf_ar1(Thetas_7_ar1, 0.7)
plots_opt_inf_ar1(Thetas_9_ar1, 0.9)


mtext(expression(gamma), 
      side = 1, outer = TRUE, line = 2, cex = 2) # Common x-axis label
mtext(expression(log(Phi(-Theta))), 
      side = 2, outer = TRUE, line = 2, cex = 2) # Common y-axis label


#mtext(text = bquote("The plot of" ~ Theta(lambda) ~ "for" ~ lambda == lambda[opt] ~ ","~ 0 ~","~ infinity), outer = TRUE, cex = 3)





############    4.3.4   Theta(optimal lambda) vs Theta(inf) vs Theta(0), 0<gamma<1





plots_opt_0_ar1 <- function(Thetas, rho){
  gamma_vec <- seq(0.05, 1, by = 0.05) # length  == 20
  len <- length(gamma_vec)
  
  thetas_opt <- Thetas[,'thetas_opt'][1:len]
  thetas_0   <- Thetas[,'thetas_0'][1:len]
  thetas_inf <- Thetas[,'thetas_inf'][1:len]
  
  th_max <- max(thetas_opt, thetas_0, thetas_inf)
  th_min <- min(thetas_opt, thetas_0, thetas_inf)
  
  # 
  plot(gamma_vec, log(pnorm(-thetas_opt)), type = 'p', cex = 0.5, col = 'red',
       xlab = expression(gamma), ylab = expression(log(Phi(-Theta))), 
       ylim = c(log(pnorm(-th_max)), log(pnorm(-th_min))),
  )
  title(main = bquote(rho == .(rho)) , cex.main = 3)
  lines(smooth.spline(gamma_vec, log(pnorm(-thetas_opt))), col = 'red', lwd = 2)
  
  
  points(gamma_vec, log(pnorm(-thetas_0)), cex = 0.5, col = 'blue')
  lines(smooth.spline(gamma_vec, log(pnorm(-thetas_0))), col = 'blue', lwd = 2)
  
  points(gamma_vec, log(pnorm(-thetas_inf)), cex = 0.5, col = 'green')
  lines(smooth.spline(gamma_vec, log(pnorm(-thetas_inf))), col = 'green', lwd = 2)
  
  legend("bottomright", 
         legend = c(expression(lambda[opt]), 0, expression(infinity)),
         col = c('red','blue', 'green'), pch = 20, cex = 1.2) #lwd = 2
}

dev.new()
par(mfrow = c(1, 4), oma = c(5, 5, 3, 0), mar = c(2, 2, 4, 2)) #

plots_opt_0_ar1(Thetas_3_ar1, 0.3)
plots_opt_0_ar1(Thetas_5_ar1, 0.5)
plots_opt_0_ar1(Thetas_7_ar1, 0.7)
plots_opt_0_ar1(Thetas_9_ar1, 0.9)

mtext(expression(gamma), 
      side = 1, outer = TRUE, line = 2, cex = 2) # Common x-axis label
mtext(expression(log(Phi(-Theta))), 
      side = 2, outer = TRUE, line = 2, cex = 2) # Common y-axis label

######


###################        4.4  when Sigma = Matern covariance function


####    cf) plot of Matern covariance functions



curve(0.5^x,0,5,col = 'red', main = "Outline of AR(1) curves vs Materns")  
# AR(1), decrease more slowly than Matern
curve(exp(-x/0.5), 0, 5, add = TRUE)                  # Matern
curve((1 + sqrt(3)*x/0.5)*exp(-sqrt(3)*x/0.5) ,0,5, col = 'green', add=TRUE)
curve(exp(-x^2/(2*0.5^2)),0,5, col = 'blue', add=TRUE)




##########       4.4.1      Plot of Theta(lambda)




get.theta.matern <- function(n,p,rho) {
  alpha <- 2
  Y = c(rep(0, n/2), rep(1, n/2))
  
  mu = replicate(2, alpha * rnorm(p) / sqrt(p)) # p by 2 matrix
  
  # nu = 1/2
  Sigma = outer(1:p, 1:p, function(x, y) exp(-abs(x-y)/rho) ) 
  
  Z = rmvnorm_rcpp(n = n, mean = matrix(0, nrow = 1,ncol = p), sigma = Sigma)
  X = t(mu)[Y + 1,] + Z # refer to page 15 or 20
  
  stj = get.stieltjes.func(X, Y)
  
  theta = function(lambda) {
    alpha^2 * stj$tau(lambda) / sqrt(alpha^2 * stj$eta(lambda) + stj$xi(lambda))
  }
  return(theta)
}

#str(200*c(0.5, 1, 3))

draw_plots_matern <- function(n, r){
  gammas <- c(1/2, 1, 3)
  alpha <- 2
  rho <- r
  
  if(rho == 0.1) lmax <- 1000
  if(rho == 0.3) lmax <- 200
  if(rho == 0.5) lmax <- 100
  if(rho == 0.7) lmax <- 30
  if(rho == 0.9) lmax <- 20
  if(rho == 0.99) lmax <- 20

  #fs <- list()
  colors <- 1:length(gammas)
 
  

  for(i in colors){
    gamma <- gammas[i]
    theta.s <- get.theta.matern(n, n*gamma, rho)
    
    max_th <- max(theta.s(seq(0,lmax, length.out = 101)[-1] ) )
    min_th <- min(theta.s(seq(0,lmax, length.out = 101)[-1] ) )

    plot(function(ll)theta.s(ll), 0, lmax, lwd = 1, col = 'brown',
         #xlab = expression(lambda), ylab = expression(Theta), 
         xlab = "", ylab = "",
         ylim = c(min_th*0.95,max_th*1.02))
    title(main = bquote(gamma == .(gammas[i]) ~ "," ~ rho == .(rho) ), cex.main = 3)
    
    for(j in 1:19){
      theta.s <- get.theta.matern(n, n*gamma, rho)
      plot(function(ll)theta.s(ll), 0, lmax, lwd = 1, col = 'brown', add = TRUE)
    }
    
  }
  
  #dev.off()

  
}



n <- 200

set.seed(2023311161)

dev.new()
par(mfrow = c(2,3), oma = c(5, 5, 3, 0), mar = c(2, 2, 4, 2)) #

#draw_plots_matern(n, 0.3)
# 0.1, 0.3 의 경우 (즉 rho가 0에 가까울수록) optimal lambda = inf
draw_plots_matern(n, 0.5)
#draw_plots_matern(n, 0.7)
draw_plots_matern(n, 0.9)
#draw_plots_matern(n, 0.99)


mtext(expression(lambda), 
      side = 1, outer = TRUE, line = 2, cex = 2) # Common x-axis label
mtext(expression(Theta), 
      side = 2, outer = TRUE, line = 2, cex = 2) # Common y-axis label



############       4.4.2   optimal lambda




lambda_generate_matern <- function(rho){
  num <- 1000
  lmin <- 0.01
  
  #if(rho == 0.1) ls <- seq(0.01, 1000, length.out = num)
  if(rho == 0.3) ls <- seq(lmin, 200,  length.out = num)
  if(rho == 0.5) ls <- seq(lmin, 100,   length.out = num)
  if(rho == 0.7) ls <- seq(lmin, 30,   length.out = num)
  if(rho == 0.9) ls <- seq(lmin, 20,   length.out = num)
  if(rho == 0.99) ls <- seq(lmin, 20,  length.out = num)
  
  return(ls)
}



optimal_lambda_matern <- function(n, rho){
  
  ls <- lambda_generate_matern(rho)
  #gamma_vec <- seq(0.10, 4.00, by = 0.02)
  
  opt.ls <- c()
  thetas_list <- list()
  thetas_opt <- c()
  thetas_0 <- c()
  thetas_inf <- c()
  
  ls.max.cand <- c()
  
  ps <- as.integer(n*gamma_vec)
  
  for(p in ps){ 
    idx <- which(ps == p)
    
    if(rho <= 0.5){
      for(j in 1:3){
        theta <- get.theta.matern(n,p, rho)
        ls.max.cand[j] <- ls[which.max(theta(ls))]
      }
      opt.ls[idx] <- min(ls.max.cand)
    } else{
      theta <- get.theta.matern(n,p, rho)
      opt.ls[idx] <- ls[which.max(theta(ls))]
    }

    thetas_list[[idx]] <- theta
  }
  
  # Theta(optimal lambda) vs gamma
  for(i in 1:length(ps)){
    thetas_opt[i] <- thetas_list[[i]](opt.ls[i])
  }
  # Theta(0) vs gamma
  unt <- which(gamma_vec == 1)
  for(i in 1:length(ps)){
    if(i > unt){
      thetas_0[i] <- NA
    }else{
      thetas_0[i] <- thetas_list[[i]](min(ls)) # 0.01
    }
  }
  # Theta(infty) vs gamma
  for(i in 1:length(ps)){
    thetas_inf[i] <- thetas_list[[i]](10^6)
  }
  
  return(cbind(opt.ls, thetas_opt, thetas_0, thetas_inf))
}



set.seed(2023311161)
matern_result_5  <- optimal_lambda_matern(n, 0.5)
matern_result_7  <- optimal_lambda_matern(n, 0.7)
matern_result_9  <- optimal_lambda_matern(n, 0.9)
matern_result_99 <- optimal_lambda_matern(n, 0.99)


#gamma_vec <- seq(0.10, 3.00, by = 2/100)
#gamma_vecnew <- as.integer(100*gamma_vec)/100
#which(gamma_vecnew == 1.50)

plots_optimal_lambda_matern <- function(result, rho){
  opt.ls <- result[,'opt.ls']
  #gamma_vec <- seq(0.10, 4.00, by = 0.02)
  
  ymax <- max(opt.ls)

  plot(gamma_vec, opt.ls, type = 'p', cex = 0.5, ylim = c(0, ymax),
       #xlab = expression(gamma), ylab = expression(Optimal ~ lambda), 
       xlab = "", ylab = "",
       main = bquote(rho == .(rho)) , cex.main = 3)
  lines(smooth.spline(gamma_vec, opt.ls), col = "skyblue", lwd = 2)
}

dev.new()
par(mfrow = c(1, 4), oma = c(5, 5, 3, 0), mar = c(2, 2, 4, 2)) #
#par(mfrow = c(1,4))

plots_optimal_lambda_matern(matern_result_5, 0.5)
plots_optimal_lambda_matern(matern_result_7, 0.7)
plots_optimal_lambda_matern(matern_result_9, 0.9)
plots_optimal_lambda_matern(matern_result_99, 0.99)


mtext(expression(gamma), 
      side = 1, outer = TRUE, line = 2, cex = 2) # Common x-axis label
mtext(expression(lambda[opt]), 
      side = 2, outer = TRUE, line = 2, cex = 2) # Common y-axis label

#mtext(text = bquote(lambda[opt] ~ "vs" ~ gamma ~ ", Matern"), outer = TRUE, cex = 3)

#par(mfrow = c(1,1), oma = c(0, 0, 0, 0))




#######          4.4.3   Theta(optimal lambda) vs Theta(inf), gamma > 0





Thetas_5_matern  <- matern_result_5[,-1]
Thetas_7_matern  <- matern_result_7[,-1]
Thetas_9_matern  <- matern_result_9[,-1]
Thetas_99_matern <- matern_result_99[,-1]


plots_opt_inf_matern <- function(Thetas, rho){
  
  thetas_opt <- Thetas[,'thetas_opt']
  thetas_inf <- Thetas[,'thetas_inf']
  
  th_max <- max(thetas_opt, thetas_inf)
  th_min <- min(thetas_opt, thetas_inf)
  
  plot(gamma_vec, log(pnorm(-thetas_opt)), type = 'p', cex = 0.5, col = 'red',
       #xlab = expression(gamma), ylab = expression(log(Phi(-Theta))), 
       xlab = "", ylab = "",
       ylim = c(log(pnorm(-th_max)), log(pnorm(-th_min))) )
  title(main = bquote(rho == .(rho)), cex.main = 3 )
  
  lines(smooth.spline(gamma_vec, log(pnorm(-thetas_opt))), col = 'red', lwd = 2)
  
  points(gamma_vec, log(pnorm(-thetas_inf)), cex = 0.5, col = 'green')
  lines(smooth.spline(gamma_vec, log(pnorm(-thetas_inf))), col = 'green', lwd = 2)
  
  legend("bottomright", legend = c(expression(lambda[opt]), expression(infinity)), 
         col = c('red','green'), pch = 20, title = expression(lambda), cex = 1.5) #lwd = 2, 
}




dev.new()
par(mfrow = c(1, 4), oma = c(5, 5, 1, 0), mar = c(2, 2, 4, 2)) #
#par(mfrow = c(1,4))


plots_opt_inf_matern(Thetas_5_matern, 0.5)
plots_opt_inf_matern(Thetas_7_matern, 0.7)
plots_opt_inf_matern(Thetas_9_matern, 0.9)
plots_opt_inf_matern(Thetas_99_matern, 0.99)


mtext(text = expression(gamma), 
      side = 1, line = 2, outer = TRUE, cex = 2)
mtext(text= bquote(log ~ Phi(-Theta)),
      side = 2, line = 2, outer = TRUE, cex = 2)




#######  4.4.4 Theta(optimal lambda), Theta(inf), theta(0) for 0<gamma<1




plots_opt_0_matern <- function(Thetas, rho){
  gamma_vec <- seq(0.05, 1, by = 0.05) # length  == 20
  len <- length(gamma_vec)
  
  thetas_opt <- Thetas[,'thetas_opt'][1:len]
  thetas_0   <- Thetas[,'thetas_0'][1:len]
  thetas_inf <- Thetas[,'thetas_inf'][1:len]
  
  th_max <- max(thetas_opt, thetas_0, thetas_inf)
  th_min <- min(thetas_opt, thetas_0, thetas_inf)
  
  # 
  plot(gamma_vec, log(pnorm(-thetas_opt)), type = 'p', cex = 0.5, col = 'red',
       #xlab = expression(gamma), ylab = expression(log(Phi(-Theta))), 
       xlab = "", ylab = "",
       ylim = c(log(pnorm(-th_max)), log(pnorm(-th_min)))
  )
  title(main = bquote(rho == .(rho)) ,cex.main = 3)
  lines(smooth.spline(gamma_vec, log(pnorm(-thetas_opt))), col = 'red', lwd = 2)
  
  
  points(gamma_vec, log(pnorm(-thetas_0)), cex = 0.5, col = 'blue')
  lines(smooth.spline(gamma_vec, log(pnorm(-thetas_0))), col = 'blue', lwd = 2)
  
  points(gamma_vec, log(pnorm(-thetas_inf)), cex = 0.5, col = 'green')
  lines(smooth.spline(gamma_vec, log(pnorm(-thetas_inf))), col = 'green', lwd = 2)
  
  legend("bottomright", legend = c(expression(lambda[opt]), 0, expression(infinity)), 
         col = c('red','blue', 'green'), pch=20, cex = 1.2)
}


dev.new()
par(mfrow = c(1, 4), oma = c(5, 6, 1, 0), mar = c(2, 2, 4, 2)) #


plots_opt_0_matern(Thetas_5_matern, 0.5)
plots_opt_0_matern(Thetas_7_matern, 0.7)
plots_opt_0_matern(Thetas_9_matern, 0.9)
plots_opt_0_matern(Thetas_99_matern, 0.99)

mtext(text = expression(gamma), 
      side = 1, line = 2, outer = TRUE, cex = 2)
mtext(text= bquote(log ~ Phi(-Theta)),
      side = 2, line = 2, outer = TRUE, cex = 2)




#####################      

