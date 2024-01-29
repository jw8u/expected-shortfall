if (!require("pacman")) install.packages("pacman")
pacman::p_load(cvar, TLMoments, ROOPSD)

A_ = function(alpha, y){
  N = length(y)
  A_alpha = (floor(N*alpha)+1-(N*alpha))*sort(y)[floor(N*alpha)] +
    (N*alpha - floor(N*alpha))*sort(y)[floor(N*alpha)+1]
  return(A_alpha)
}
s2_ = function(alpha, y, A){
  z = qnorm(alpha)
  x1 = (z^2 + 1 - (z/((1-alpha)*sqrt(2*pi)))*exp(-(z^2)/2))^-1
  x2 = sum((y-A)^2 * ifelse(y>A,1,0))/sum(ifelse(y>A,1,0))
  return(x1*x2)
}
mu_ = function(alpha, A, s2){
  z = qnorm(alpha)
  return(A - sqrt(s2)*z)
}
skew_ = function(y, A){
  x1 = sum((y-A)^3 * ifelse(y>A,1,0))/sum(ifelse(y>A,1,0))
  x2 = (sum((y-A)^2 * ifelse(y>A,1,0))/sum(ifelse(y>A,1,0)))^-1.5
  return(x1*x2)
}
R_ = function(skew, beta){
  if (beta==0.99){
    return(0.8611 + 0.5191*exp(-0.9747*skew) + 0.6099/skew -
             0.9413/(skew^2))
  }
  if (beta==0.995){
    return(0.9919 + 0.6681*exp(-0.9607*skew) + 0.6022/skew -
             1.4623/(skew^2))
  }
}
ES_adj_ = function(mu, s2, beta, A, R){
  z = qnorm(beta)
  ES_X = mu + sqrt(s2)/((1-beta)*sqrt(2*pi)) * exp(-(z^2)/2)
  ES_adj = (ES_X-A)*R + A
  return(ES_adj)
}
ES_W_ = function(beta, distr, df, shape, mu, s, scale, xi){
  if (distr == "t"){
    ES_W = ES(pt, p_loss = 1-beta, dist.type = "cdf", df=df)
  }
  if (distr == "gpd"){
    #xi != 0
    ES_W = mu + s*(((1-beta)^(-xi))/(1-xi) + ((1-beta)^(-xi) - 1)/xi)
  }
  if (distr == "weibull"){
    f = function(x) scale/(1-beta) * x^((1+1/shape)-1) * exp(-x)
    ES_W = (integrate(f,-log(1-beta),+Inf))$value
  }
  return(ES_W)
}
error_ = function(ES_W, ES_adj){
  return((ES_W - ES_adj)/(ES_W))
}
calculate_ES_adj = function(alpha, beta, y){
  A = A_(alpha,y)
  s2 = s2_(alpha, y, A)
  mu = mu_(alpha, A, s2)
  skew = skew_(y,A)
  R = R_(skew, beta)
  ES_adj = ES_adj_(mu, s2, beta, A, R)
  return(ES_adj)
}
calculate_error = function(alpha, beta, y, distr, ...){
  A = A_(alpha,y)
  s2 = s2_(alpha, y, A)
  mu = mu_(alpha, A, s2)
  skew = skew_(y,A)
  R = R_(skew, beta)
  ES_adj = ES_adj_(mu, s2, beta, A, R)
  ES_W = ES_W_(beta, distr, ...)
  error = error_(ES_W, ES_adj)
  return(error)
}

#Arithmetic Average
AA_ = function(beta, y){
  N = length(y)
  AA = sum(sort(y)[ceiling(N*beta):N])/(N+1-ceiling(N*beta))
  return(AA)
}
#Monte Carlo simulations
print_t = function(M){
  for (df in c(3.5,5,8)){
    for (n in c(250,500)){
      for (beta in c(0.99,0.995)){
        ES_W = ES_W_(beta=beta, distr="t", df=df)
        set.seed(001)
        print(mean(replicate(M, (calculate_ES_adj(alpha=0.95, beta=beta,
                                                  y=rt(n=n, df=df) - ES_W))^2)))
        set.seed(001)
        print(mean(replicate(M, (AA_(beta=beta,y=rt(n=n, df=df) -
                                       ES_W))^2)))
        set.seed(001)
        print(mean(var(replicate(M, calculate_ES_adj(alpha=0.95,
                                                     beta=beta, y=rt(n=n, df=df))))))
        set.seed(001)
        print(mean(var(replicate(M, AA_(beta=beta, y=rt(n=n, df=df))))))
        set.seed(001)
        print(mean(replicate(M, calculate_ES_adj(alpha=0.95, beta=beta,
                                                 y=rt(n=n, df=df) - ES_W))))
        set.seed(001)
        print(mean(replicate(M, AA_(beta=beta, y=rt(n=n, df=df) - ES_W))))
      }
    }
  }
}
print_gpd = function(M){
  for (xi in c(0.3,0.2,0.1)){
    for (n in c(250,500)){
      for (beta in c(0.99,0.995)){
        ES_W = ES_W_(beta=beta, distr="gpd", mu=0, s=1, xi=xi)
        set.seed(001)
        print(mean(replicate(M, (calculate_ES_adj(alpha=0.95, beta=beta,
                                                  y=rgpd(n=n, loc=0, scale=1, shape=xi)) - ES_W)^2)))
        set.seed(001)
        print(mean(replicate(M, (AA_(beta=beta, y=rgpd(n=n, loc=0,
                                                       scale=1, shape=xi)) - ES_W)^2)))
        set.seed(001)
        print(mean(var(replicate(M, calculate_ES_adj(alpha=0.95,
                                                     beta=beta, y=rgpd(n=n, loc=0, scale=1, shape=xi))))))
        set.seed(001)
        print(mean(var(replicate(M, AA_(beta=beta, y=rgpd(n=n, loc=0,
                                                          scale=1, shape=xi))))))
        set.seed(001)
        print(mean(replicate(M, calculate_ES_adj(alpha=0.95, beta=beta,
                                                 y=rgpd(n=n, loc=0, scale=1, shape=xi)) - ES_W)))
        set.seed(001)
        print(mean(replicate(M, AA_(beta=beta, y=rgpd(n=n, loc=0, scale=1,
                                                      shape=xi)) - ES_W)))
      }
    }
  }
}
print_weibull = function(M){
  for (shape in c(0.6,0.9,1.4)){
    for (n in c(250,500)){
      for (beta in c(0.99,0.995)){
        ES_W = ES_W = ES_W_(beta=beta, distr="weibull", shape=shape,
                            scale=1)
        set.seed(001)
        print(mean(replicate(M, (calculate_ES_adj(alpha=0.95, beta=beta,
                                                  y=rweibull(n=n, shape=shape, scale=1)) - ES_W)^2)))
        set.seed(001)
        print(mean(replicate(M, (AA_(beta=beta,y=rweibull(n=n,
                                                          shape=shape, scale=1)) - ES_W)^2)))
        set.seed(001)
        print(mean(var(replicate(M, calculate_ES_adj(alpha=0.95,
                                                     beta=beta, y=rweibull(n=n, shape=shape, scale=1))))))
        set.seed(001)
        print(mean(var(replicate(M, AA_(beta=beta, y=rweibull(n=n,
                                                              shape=shape, scale=1))))))
        set.seed(001)
        print(mean(replicate(M, calculate_ES_adj(alpha=0.95, beta=beta,
                                                 y=rweibull(n=n, shape=shape, scale=1)) - ES_W)))
        set.seed(001)
        print(mean(replicate(M, AA_(beta=beta, y=rweibull(n=n,
                                                          shape=shape, scale=1)) - ES_W)))
      }
    }
  }
}
print_t2 = function(M){
  for (df in c(2.5,3)){
    for (n in c(250,500)){
      for (beta in c(0.99,0.995)){
        ES_W = ES_W_(beta=beta, distr="t", df=df)
        set.seed(001)
        print(mean(replicate(M, (calculate_ES_adj(alpha=0.95, beta=beta,
                                                  y=rt(n=n, df=df) - ES_W))^2)))
        set.seed(001)
        print(mean(replicate(M, (AA_(beta=beta,y=rt(n=n, df=df) -
                                       ES_W))^2)))
        set.seed(001)
        print(mean(var(replicate(M, calculate_ES_adj(alpha=0.95,
                                                     beta=beta, y=rt(n=n, df=df))))))
        set.seed(001)
        print(mean(var(replicate(M, AA_(beta=beta, y=rt(n=n, df=df))))))
        set.seed(001)
        print(mean(replicate(M, calculate_ES_adj(alpha=0.95, beta=beta,
                                                 y=rt(n=n, df=df) - ES_W))))
        set.seed(001)
        print(mean(replicate(M, AA_(beta=beta, y=rt(n=n, df=df) - ES_W))))
      }
    }
  }
}
print_gpd2 = function(M){
  for (xi in c(0.5,0.35)){
    for (n in c(250,500)){
      for (beta in c(0.99,0.995)){
        ES_W = ES_W_(beta=beta, distr="gpd", mu=0, s=1, xi=xi)
        set.seed(001)
        print(mean(replicate(M, (calculate_ES_adj(alpha=0.95, beta=beta,
                                                  y=rgpd(n=n, loc=0, scale=1, shape=xi)) - ES_W)^2)))
        set.seed(001)
        print(mean(replicate(M, (AA_(beta=beta, y=rgpd(n=n, loc=0,
                                                       scale=1, shape=xi)) - ES_W)^2)))
        set.seed(001)
        print(mean(var(replicate(M, calculate_ES_adj(alpha=0.95,
                                                     beta=beta, y=rgpd(n=n, loc=0, scale=1, shape=xi))))))
        set.seed(001)
        print(mean(var(replicate(M, AA_(beta=beta, y=rgpd(n=n, loc=0,
                                                          scale=1, shape=xi))))))
        set.seed(001)
        print(mean(replicate(M, calculate_ES_adj(alpha=0.95, beta=beta,
                                                 y=rgpd(n=n, loc=0, scale=1, shape=xi)) - ES_W)))
        set.seed(001)
        print(mean(replicate(M, AA_(beta=beta, y=rgpd(n=n, loc=0, scale=1,
                                                      shape=xi)) - ES_W)))
      }
    }
  }
}
print_t(2500)
print_gpd(2500)
print_weibull(2500)
print_t2(2500)
print_gpd2(2500)
