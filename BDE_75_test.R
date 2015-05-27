## Brown Durban Evans (1975) 'Homogeneity test' ----

## number of intervals ----
n <- 50

## getting data: change to pull from EIA API using EIAdata ----
data.test <- cbind(as.vector(data$ng.ret),as.vector(data$stor.2.yr),as.vector(data$hdd.dev), as.vector(data$cdd.dev))

p <- floor(dim(data.test)[1]/n)

data.test.spl <- array(NA,dim=c(n,dim(data.test)[2],p))
for (i in 1:p){
  data.test.spl[,,i] <- data.test[(1+(i-1)*n):(i*n),]}

## get the residual sum of squares (rss) for each interval
rss <- rep(NA,p)
for (i in 1:p){
  rss[i] <- t(data.test.spl[,1,i])%*%data.test.spl[,1,i]-t(data.test.spl[,1,i])%*%data.test.spl[,2:4,i]%*%solve(t(data.test.spl[,2:4,i])%*%data.test.spl[,2:4,i])%*%t(data.test.spl[,2:4,i])%*%data.test.spl[,1,i]
}

## get the full sample rss
rss.full.sam <- t(data.test[,1])%*%data.test[,1]-t(data.test[,1])%*%data.test[,2:4]%*%solve(t(data.test[,2:4])%*%data.test[,2:4])%*%t(data.test[,2:4])%*%data.test[,1]

## test statistic
k <- dim(data.test)[2]-1
test.stat <- ((dim(data.test)[1]-(k)*p)/((k)*p-(k)))*((rss.full.sam-sum(rss))/sum(rss))

## under Ho distributed as F(kp-k,T-kp)
## p-value
1-pf(test.stat, k*p-k, dim(data.test)[1]-k*p)

