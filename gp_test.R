# Create data:
set.seed(1)
f <- function(x) x^2 + x
limits <- function(samples) { #Samples are of size N times J, where N is the number of repetitions
  return(apply(samples, 2, function(s) quantile(s, c(.022, .158, .5, 0.842, 0.978))))
}
plotit <- function(x,y, xt, yt, x_legend, xt_legend, title="") {
  library(data.table)
  fnew = limits(yt)
  fnew <- as.data.frame(fnew)
  setattr(fnew, "row.names", c("1","2","y", "4", "5"))
  s <- data.frame('x'=x, 'y'=y)
  fi <- data.frame('x'=xt, 'id'=t(fnew))
  pl = ggplot() 
  pl = pl +
    geom_line(data=fi, aes(y=id.y, x=x, color=x_legend)) + 
    geom_point(data=s, aes(x=x,y=y, color=xt_legend)) +
    scale_color_manual(name="", values = c("red", "blue")) + 
    geom_ribbon(data=fi, aes(ymin=id.1,ymax=id.5,x=x), fill="blue", alpha=0.1) +
    geom_ribbon(data=fi, aes(ymin=id.2,ymax=id.4,x=x), fill="blue", alpha=0.1) +
    xlab('x') + ylab('y') + ggtitle(title) + theme(legend.position=c(0.85, 0.25))
  return(pl)
}

plotit2 <- function(xt, yt,x=NULL,y=NULL, nsamples=10, title="") {
  library(data.table)
  fnew = limits(yt)
  fnew <- as.data.frame(fnew)
  setattr(fnew, "row.names", c("1","2","y", "4", "5"))
  fi <- data.frame('x'=xt, 'id'=t(fnew))
  pl = ggplot() 
  pl = pl +
    geom_line(data=fi, aes(y=id.y, x=x), color="blue") + 
    geom_ribbon(data=fi, aes(ymin=id.1,ymax=id.5,x=x), fill="blue", alpha=0.1) +
    geom_ribbon(data=fi, aes(ymin=id.2,ymax=id.4,x=x), fill="blue", alpha=0.1)
  for(i in 1:nsamples) {
    ij = sample(1:dim(yt)[1],1,replace=T)
    pl = pl + geom_line(data=data.frame('x'=xt, 'y' = yt[ij,]), aes(y=y, x=x), color="black")
  }
  if(!is.null(x)) {
    pl = pl + geom_point(data=data.frame('x'=x, 'y'=y), aes(y=y,x=x), color="red", size=3)
  }
  pl = pl + xlab('x') + ylab('y') + ggtitle(title)
  return(pl)
}

x <- rnorm(50, mean = 0, sd = 0.5)
y <- f(x) + rnorm(length(x), mean=0, sd=0.1)

plot(x,y)

X = cbind(rep(1, length(x)),x)
b = solve((t(X) %*% X), t(X) %*% y)

xt = seq(from=-1.25, to=1.0, length.out = length(x))
ys = b[1] + b[2]*xt

library(ggplot2)
theme_set(theme_bw())
fr <- data.frame('x'=x, 'y'=y)
fr_l <- data.frame('x'=xt, 'y'=ys)
ggplot() +
  geom_point(data=fr, aes(y=y, x=x, color="data")) +
  geom_line(data=fr_l, aes(y=y, x=x, color="fit")) + 
  scale_color_manual(name="", values = c("blue","red"))

library(rstan)
stan_model <- stan_model("stan/linear.stan")
dat = list(N1=length(x),x1=x,y1=y,N2=length(xt),x2=xt)
fit <- sampling(stan_model, data=dat, chains=4, warmup=1000, iter=2000)
samples <- extract(fit)

plotit(x,y, xt, samples$f2,  "fit", "data", title='Linear model fit')
ggsave('pictures/linear.pdf',device="pdf", width=7,height=4)

stan_model2 <- stan_model("stan/gp.stan")
fit2 <- sampling(stan_model2, data=dat, chains=4, warmup=1000, iter=2000)
samples2 <- extract(fit2)

plotit(x,y, xt, samples2$f2,  "fit", "data", title='Gaussian process fit')
ggsave('pictures/gp.pdf',device="pdf", width=7,height=4)

xt2 = seq(from=min(x)-1, to=max(x)+2, length.out = length(x)*2)
dat2 = list(N1=length(x),x1=x,y1=y,N2=length(xt2),x2=xt2)
fit3 <- sampling(stan_model2, data=dat2, chains=4, warmup=1000, iter=2000)
samples3 <- extract(fit3)

plotit(x,y, xt2, samples3$f2,  "fit", "data", title='Gaussian process fit')
ggsave('pictures/gp_extrapolation.pdf',device="pdf", width=7,height=4)

set.seed(5)
Nxt=0
Nxt= 3
xx = as.array(rep(1,0)) 
xx = sort(as.array(runif(Nxt,-3, 3)))
yx = as.array(rep(1,0))
yx = sin(xx)
xx

N1d=0
x1d = as.array(rep(1,0))
N1df=0
x1df = as.array(rep(1,0))
ydf = as.array(rep(1,0))
nu=100
sigma_d_f = 1e-6

stan_model3 <- stan_model("stan/gp_fixed.stan")
xt3 = seq(from=-3, to=3, length.out = 100)
dat3 = list(N1=Nxt,x1=xx,y1=yx,N2=length(xt3),x2=xt3, rho=1, alpha=1, sigma=0.1,N1d=N1d,N1df=N1df,x1d=x1d,x1df=x1df,nu=nu,sigma_d_f=sigma_d_f,ydf=ydf)
fit4 <- sampling(stan_model3, data=dat3, chains=4, warmup=1000, iter=2000)
samples4 <- extract(fit4)

#plotit2(xt3, samples4$f2, nsamples=3, title="Prior distribution")
#ggsave('pictures/smoother.pdf',device="pdf", width=7,height=4)
plotit2(xt3, samples4$f2,x=xx, y=yx, nsamples=3, title="Posterior fit with with three training samples")
ggsave('pictures/smoother_fit.pdf',device="pdf", width=7,height=4)

dat4 = list(N1=Nxt,x1=xx,y1=yx,N2=length(xt3),x2=xt3, rho=0.5, alpha=1, sigma=0.1)
fit5 <- sampling(stan_model3, data=dat4, chains=4, warmup=1000, iter=2000)
samples5 <- extract(fit5)

#plotit2(xt3, samples5$f2,x=xx, y=yx, nsamples=3, title='')
#plotit2(xt3, samples5$f2, nsamples=3, title='Prior distribution')
#ggsave('pictures/less_smoot.pdf',device="pdf", width=7,height=4)

Nxt=1
xx = as.array(rep(3,Nxt)) 
yx = as.array(rep(0,Nxt))
N1d=30
x1d = as.array(seq(from=-3, to=2.9,length.out = N1d))
N1df=0
x1df = as.array(rep(3,N1df))
ydf = as.array(rep(0,N1df))
nu=100
sigma_d_f = 1e-3
xt3 = seq(from=-3, to=3, length.out = 100)

dat5 = list(N1=Nxt,x1=xx,y1=yx,N2=length(xt3),x2=xt3, rho=1, alpha=1, sigma=0.01,N1d=N1d,N1df=N1df,x1d=x1d,x1df=x1df,nu=nu,sigma_d_f=sigma_d_f,ydf=ydf)
fit6 <- sampling(stan_model3, data=dat5, chains=4, warmup=200, iter=1200)
samples4 <- extract(fit6)

plotit2(xt3, samples4$f2, x=xx, y=yx, nsamples=3, title="Posterior distribution for monotonic GP")
ggsave('pictures/monotonic_post_gp',device="pdf", width=7,height=4)

library(pracma)

x = seq(from=-3, to=3, length.out=100)
y = erf(x)/2 +0.5

fi <- data.frame('x'=x, 'y'=y)
pl = ggplot() 
pl = pl + geom_line(data=fi, aes(y=y, x=x), color="black") + ggtitle("Probit likelihood value as a function of latent function derivative value") + xlab("Latent function derivative value") + ylab("Likelihood")
pl
ggsave('pictures/probit.pdf',device="pdf", width=7,height=4)

h=2
k=1
x = seq(from=-0.7, to=5, length.out=100)
y = (x+0.7)^h / (k^h + (x+0.7)^h) 
fi <- data.frame('x'=x, 'y'=y)
pl = ggplot() 
pl = pl + geom_line(data=fi, aes(y=y, x=x), color="black") + ggtitle("Hill function value as a function of age (h=2, K_{age}=1)") + xlab("Age [y]") + ylab("Maturation")+ geom_vline(xintercept=0, linetype="dashed", color = "black") +
  scale_x_continuous(breaks=c(-0.7,0,1,2,3,4,5), labels=c("-0.7", "0", "1", "2", "3", "4", "5"))

pl
ggsave('pictures/hill.pdf',device="pdf", width=7,height=4)
