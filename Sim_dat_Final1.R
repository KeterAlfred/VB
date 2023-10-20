
#
#Script to simulate data mimicking TB prevalence surveys
#
#
#=============================================================================

# Clear the memory
rm(list = ls())


N = 10000


M = 100

varnames = c('id','Y1','Y2',"Y25",'Y3','Y4','Y5','Y6','Y50','Y60','V','D','groups')
sim_dat_Final1 = array(dim = c(N,13,M), dimnames = list(NULL,
                                                        varnames,
                                                        NULL))



for (m in 1:M){
  
  set.seed(m)
  
  
  id = seq(1,N,1)
  
  
  
  
  D = rbinom(N,1,0.02)
  prop.table(table(D))
  
  
  a0 = -2.123
  b0 = -2.845
  
  LP0 = b0*D + a0*(1-D)
  Y0 = rbinom(N, 1, 1/( 1 + exp(- LP0 ) ) )
  table(Y0, D)
  prop.table(table(Y0, D),1)
  
  
  
  #Coefficients for Y0 in the following models for Y1-Y4
  #This causes asymptomatic with normal chest X-ray cases
  #to be generated among the TB cases
  c1 = -2
  c2 = -2
  c3 = -1.52
  c4 = -2.52
  
  
  
  
  #Generate Y1
  a1 = c(-1.823, c1)
  b1 = c(-1.323, c1)
  
  M0 = cbind(1,Y0)
  
  LP1 = (M0 %*% b1)*D + (M0 %*% a1)*(1-D)
  Y1 = rbinom(N, 1, 1/( 1 + exp(- LP1 ) ) )
  table(Y1, D)
  prop.table(table(Y1, D),2)
  
  
  
  
  M1 = cbind(1,Y0,Y1)
  
  a2 = c(-1.47, c2, 1.9)
  b2 = c(0.79, c2, 0)
  
  LP2 = (M1 %*% b2)* D + (M1 %*% a2)*(1-D)
  Y2 = rbinom(N, 1, 1/( 1 + exp(- LP2 ) ) )
  table(Y2, D)
  prop.table(table(Y2, D),2)
  
  
  
  
  
  
  
  
  
  # Simulate CAD4TB score
  kappa = 9.98
  gamma0 =   gamma0 = -0.930 #-0.530
  gamma1 = 0.536
  v1 = 1/ ( 1 + exp(-(gamma1 + c3*Y0 + 1.8*Y2)))
  v0 = 1/ ( 1 + exp(-(gamma0 + c3*Y0 + 1.1*Y1 + 0.75*Y2)))
  gA1 = v1*kappa
  gB1 = (1-v1)*kappa
  gA0 = v0*kappa
  gB0 = (1-v0)*kappa
  
  gA = gA1*D + gA0*(1-D)
  gB = gB1*D + gB0*(1-D)
  
  cad6 = rbeta(N,shape1 = gA, shape2 = gB)
  
  
  hist(cad6[D==1])
  hist(cad6[D==0])
  
  plot(density(cad6[D==1]))
  lines(density(cad6[D==0]),col = 'blue')
  
  
  
  Y3 = (cad6>=.53)*1
  table(Y3)
  prop.table(table(Y3,D),2)
  cor(cbind(Y1[D==1], Y2[D==1], Y3[D==1]))
  cor(cbind(Y1[D==0],Y2[D==0], Y3[D==0]))
  
  
  
  
  
  
  # Simulate CAD4TB score
  kappa.2 = 9.98
  gamma0.2 = -3.330
  gamma1.2 = -1.936
  v1.2 = 1/ ( 1 + exp(-(gamma1.2 + c4*Y0 + 1.2*Y2 + 1.3*cad6)))
  v0.2 = 1/ ( 1 + exp(-(gamma0.2 + c4*Y0 + 1.32*Y1 + 2.65*Y2 + 1.25*cad6)))
  gA1.2 = v1.2*kappa.2
  gB1.2 = (1-v1.2)*kappa.2
  gA0.2 = v0.2*kappa.2
  gB0.2 = (1-v0.2)*kappa.2
  
  gA.2 = gA1.2*D + gA0.2*(1-D)
  gB.2 = gB1.2*D + gB0.2*(1-D)
  
  cad7 = rbeta(N,shape1 = gA.2, shape2 = gB.2)
  
  
  hist(cad7[D==1])
  hist(cad7[D==0])
  
  plot(density(cad7[D==1]))
  lines(density(cad7[D==0]),col = 'blue')
  
  
  
  Y4 = (cad7>=.15)*1
  table(Y4)
  prop.table(table(Y4,D),2)
  cor(cbind(Y1[D==1], Y2[D==1], Y3[D==1], Y4[D==1]))
  cor(cbind(Y1[D==0],Y2[D==0], Y3[D==0], Y4[D==0]))
  
  
  
  
  
  #Sample Y5 & Y6
  
  a5 = -5
  b5 = 0.56
  
  
  LP5 = b5*D + a5*(1-D)
  Y5 = rbinom(N, 1, 1/( 1 + exp(- LP5 ) ) )
  table(Y5, D)
  prop.table(table(Y5, D),2)
  
  
  
  a6 = -5
  b6 = c(0.275, 2.67)
  
  M2 = cbind(1, Y5)
  
  LP6 = (M2 %*% b6)*D + a6*(1-D)
  Y6 = rbinom(N, 1, 1/( 1 + exp(- LP6 ) ) )
  table(Y6, D)
  prop.table(table(Y6, D),2)
  
  
  cor(cbind(Y1[D==1], Y2[D==1], Y3[D==1], Y4[D==1], Y5[D==1], Y6[D==1]))
  cor(cbind(Y1[D==0],Y2[D==0], Y3[D==0], Y4[D==0], Y5[D==0], Y6[D==0]))
  
  
  
  
  Y25 = ( cad6>= .25 )*1
  E = ifelse( Y1 == 1 | Y25 == 1 | Y2 == 1, 1, 0 )
  
  
  
  
  #Flag unverified subjects though eligible
  # exclude.id = ifelse(id %in% sample(id[cad6>=.25 & cad6<.60], 2000),1,0)
  exclude.id = ifelse(id %in% sample(id[cad6>=.25 & cad6<.60], 300),1,0)
  failed.id = ifelse( id %in% sample(id[ E == 1 & exclude.id == 0], 1200), 1, 0)
  
  
  
  
  
  # V = ifelse( E %in% 0 | id %in% failed.id | id %in% exclude.id, 0, 1) 
  # | 
  # cad6>= .25 & !id %in% failed.id & !id %in% exclude.id)*1
  
  
  
  
  groups = ifelse(E %in% 1 & failed.id == 0 & exclude.id == 0, 1,  #Eligible & tested
                  ifelse(E %in% 1 & failed.id == 1 & exclude.id == 0, 2,  #Eligible but didn't test
                         ifelse(E %in% 1 & failed.id == 0 & exclude.id == 1, 3, 4))) #Eligible but were excluded using previous criterion
  
  V = ifelse( groups %in% 1, 1, 0) 
  
  prop.table(table(V))
  prop.table(table(V,D),1)
  
  table(groups)
  
  prop.table( table(groups, D),1)
  
  table(groups, Y1)
  prop.table(table(V,D),1)
  prop.table(table(D))
  table(groups, D)
  prop.table(table(groups, D),1)
  
  prop.table(table(groups))
  
  table(groups,V)
  
  table(groups, Y5)
  table(groups, V, E)
  
  
  
  Y5Y6 = ( Y5 | Y6)*1
  
  
  prop.table(table(Y5Y6))
  prop.table(table(D))
  
  
  
  table(Y1,Y2,Y3,Y4,Y5,Y6,V)
  
  
  
  
  Y50 = ifelse( V==0, 9, Y5 )
  Y60 = ifelse( V==0, 9, Y6 )
  
  
  
  sim_dat_Final1[,,m] = array(cbind(id,Y1,Y2,Y25,Y3,Y4,Y5,Y6,Y50,Y60,V,D,groups))
  
}







#Check the distribution of some variables

prop.table(table(sim_dat_Final1[,"D",1]))

prop.table(table(sim_dat_Final1[,c("Y1"),1],sim_dat_Final1[,c("D"),1]),2)
prop.table(table(Y2,D),2)
prop.table(table(Y4,D),2)
prop.table(table(Y5,D),2)
prop.table(table(Y6,D),2)[2,2]







#Overall prevalence estimate in each of the simulated dataset
M=100
prev = numeric(M)
for( m in 1:M){
  prev[m] = with( as.data.frame(sim_dat_Final1[,,m]), prop.table(table(D)))[2]
}
prev
mean(prev)







#Prevalence estimates in each of the three groups (eligible and verified, 
#eligible but unverified, ineligible and unverified)


prev2 = matrix(NA,nrow=M, ncol=3)
for( m in 1:M){
  sim_dat_Final1[,"groups",m] = ifelse( sim_dat_Final1[,"groups",m] == 3,2,
                                        sim_dat_Final1[,"groups",m])
  prev2[m,] = with( as.data.frame(sim_dat_Final1[,,m]), prop.table(table(D, groups),2))[2,]
}
prev2
apply(prev2,2,mean)




#Overall sensitivity
se.Y = array(dim = c(M,5))
for( m in 1:M){
  se.Y[m,1] = with( as.data.frame(sim_dat_Final1[,,m]), prop.table(table(Y1,D),2))[2,2]
  se.Y[m,2] = with( as.data.frame(sim_dat_Final1[,,m]), prop.table(table(Y2,D),2))[2,2]
  se.Y[m,3] = with( as.data.frame(sim_dat_Final1[,,m]), prop.table(table(Y4,D),2))[2,2]
  se.Y[m,4] = with( as.data.frame(sim_dat_Final1[,,m]), prop.table(table(Y5,D),2))[2,2]
  se.Y[m,5] = with( as.data.frame(sim_dat_Final1[,,m]), prop.table(table(Y6,D),2))[2,2]
}
head(se.Y)
apply(se.Y,2,mean)



#Overall specificity
sp.Y = array(dim = c(M,5))
for( m in 1:M){
  sp.Y[m,1] = with( as.data.frame(sim_dat_Final1[,,m]), prop.table(table(Y1,D),2))[1,1]
  sp.Y[m,2] = with( as.data.frame(sim_dat_Final1[,,m]), prop.table(table(Y2,D),2))[1,1]
  sp.Y[m,3] = with( as.data.frame(sim_dat_Final1[,,m]), prop.table(table(Y4,D),2))[1,1]
  sp.Y[m,4] = with( as.data.frame(sim_dat_Final1[,,m]), prop.table(table(Y5,D),2))[1,1]
  sp.Y[m,5] = with( as.data.frame(sim_dat_Final1[,,m]), prop.table(table(Y6,D),2))[1,1]
}
head(sp.Y)


apply(sp.Y,2,mean)







#Covariances and correlations: All groups

table(sim_dat_Final1[
  sim_dat_Final1[,"D",1] %in% c(0,1),c("Y6"),1])



cov(sim_dat_Final1[
  sim_dat_Final1[,"D",1] == 0,c("Y1","Y2","Y4","Y5","Y6"),1]) #[1,2]

M = 100
K=5

covar <- array(NA,dim = c(K,K,100))
corr <- array(NA,dim = c(K,K,100))

covar.D1 = covar
corr.D1 = corr

covar.D0 = covar
corr.D0 = corr


for(m in 1:M){
  for(i in 1:(K-1)){
    for(j in (i+1):K){
      covar.D1[j,i,m] <- cov(sim_dat_Final1[
        sim_dat_Final1[,"D",m] == 1,c("Y1","Y2","Y4","Y5","Y6"),m])[j,i]
      corr.D1[i,j,m] <- cor(sim_dat_Final1[
        sim_dat_Final1[,"D",m] == 1,c("Y1","Y2","Y4","Y5","Y6"),m])[i,j]
      
      covar.D0[j,i,m] <- cov(sim_dat_Final1[
        sim_dat_Final1[,"D",m] == 0,c("Y1","Y2","Y4","Y5","Y6"),m])[j,i]
      corr.D0[i,j,m] <- cor(sim_dat_Final1[
        sim_dat_Final1[,"D",m] == 0,c("Y1","Y2","Y4","Y5","Y6"),m])[i,j]
      
      
    }
  }
}
covar.D1
corr.D1
covar.D0
corr.D0

pool.covar.D1 = apply(covar.D1,c(1,2),mean)
round( pool.covar.D1, 3)

pool.covar.D0 = apply(covar.D0,c(1,2),mean)
round( pool.covar.D0, 3)


pool.corr.D1 = apply(corr.D1,c(1,2),mean)
round(pool.corr.D1, 3)

pool.corr.D0 = apply(corr.D0,c(1,2),mean)
round(pool.corr.D0,3)




m.D1 <- array(NA,dim = c(K,100))
m.D0 <- array(NA,dim = c(K,100))
var.D1 <- array(NA,dim = c(K,100))
var.D0 <- array(NA,dim = c(K,100))

vars = c("Y1","Y2","Y4","Y5","Y6")

for(m in 1:M){
  for(i in 1:K){
    m.D1[i,m] <- mean(sim_dat_Final1[
      sim_dat_Final1[,"D",m] == 1,vars[i],m])
    
    var.D1[i,m] = m.D1[i,m] * (1 - m.D1[i,m])
    
    m.D0[i,m] <- mean(sim_dat_Final1[
      sim_dat_Final1[,"D",m] == 0,vars[i],m])      
    
    var.D0[i,m] = m.D0[i,m] * (1 - m.D0[i,m])
    
  }
}

pool.var.D1 = apply(var.D1,c(1),mean)
round( pool.var.D1, 3)

pool.var.D0 = apply(var.D0,c(1),mean)
round( pool.var.D0, 3)







#Covariances and correlations: Restricted to verified cases only

M = 100
K=5

covar <- array(NA,dim = c(K,K,100))
corr <- array(NA,dim = c(K,K,100))

covar.D1 = covar
corr.D1 = corr

covar.D0 = covar
corr.D0 = corr



for(m in 1:M){
  for(i in 1:(K-1)){
    for(j in (i+1):K){
      covar.D1[j,i,m] <- cov(sim_dat_Final1[
        sim_dat_Final1[,"D",m] == 1 & 
          sim_dat_Final1[,"V",m] == 1,c("Y1","Y2","Y4","Y5","Y6"),m])[j,i]
      corr.D1[i,j,m] <- cor(sim_dat_Final1[
        sim_dat_Final1[,"D",m] == 1 & 
          sim_dat_Final1[,"V",m] == 1,c("Y1","Y2","Y4","Y5","Y6"),m])[i,j]
      
      covar.D0[j,i,m] <- cov(sim_dat_Final1[
        sim_dat_Final1[,"D",m] == 0 & 
          sim_dat_Final1[,"V",m] == 1,c("Y1","Y2","Y4","Y5","Y6"),m])[j,i]
      corr.D0[i,j,m] <- cor(sim_dat_Final1[
        sim_dat_Final1[,"D",m] == 0 & 
          sim_dat_Final1[,"V",m] == 1,c("Y1","Y2","Y4","Y5","Y6"),m])[i,j]
      
      
    }
  }
}
covar.D1
corr.D1
covar.D0
corr.D0

pool.covar.D1 = apply(covar.D1,c(1,2),mean)
round( pool.covar.D1, 3)

pool.covar.D0 = apply(covar.D0,c(1,2),mean)
round( pool.covar.D0, 3)


pool.corr.D1 = apply(corr.D1,c(1,2),mean)
round(pool.corr.D1, 3)

pool.corr.D0 = apply(corr.D0,c(1,2),mean)
round(pool.corr.D0,3)














#Calculate the average frequency distributions across the 100 simulated datasets
#This give the cross classification of the tests by verification status

# Create the design matrix to hold the output

K=6
g = 3
DH = matrix( NA, nrow = 144, ncol = K  )  #

DH[,1] <- rep(c(1, 0), 72)
DH[,2] <- rep(c(rep(1,2), rep(0,2)),36)
DH[,3] <- rep(c(rep(1,4), rep(0,4)),18)
DH[,4] <- rep(c(rep(1,8), rep(0,8)), 9)
DH[,5] <- rep( c( rep(1,16), rep(0,16), rep(9,16) ),3)
DH[,6] <- c(rep(1,48), rep(0,48), rep(9,48))

var.names = c('V', 'Y60', 'Y50', 'Y25', 'Y2', 'Y1')
var.names2 = paste0('freq',1:100)

freqs = matrix(NA,ncol = ncol(DH), nrow = nrow(DH))
freqs[,1:ncol(DH)] <- DH

colnames(DH) = var.names
colnames(freqs) = c(var.names)

M = 100
N=10000
ones = rep(1, N)
for(m in 1:M){
  freq1 = aggregate( ones ~ V + Y60 + Y50 + Y25 + Y2 + Y1, 
                     data = sim_dat_Final1[,,m], FUN = function(x){NROW(x)} )
  freqs = merge(as.data.frame(freqs), as.data.frame(freq1), 
                by = var.names, all = TRUE)
  new.cols = paste0("freqobs",1:m)
  colnames(freqs) <- c(var.names,new.cols)
  
}
head(freqs)



freqs$na.count = apply(freqs,1,function(x)sum(is.na(x)))
freqs.subset = freqs[ freqs$na.count<100,]
freqs.subset$n = 100 - freqs.subset$na.count
for(m in 1:100){
  freqs.subset[,m+6]= ifelse( is.na(freqs.subset[,m+6]),0,
                              freqs.subset[,m+6])
}

ind1 = which(names(freqs) %in% "freqobs1")
ind2 = which(names(freqs) %in% "freqobs100")

freqs.subset$freqobs.t = apply(freqs.subset[,c(ind1:ind2)],1,sum)
min(freqs.subset$n)
freqs.subset$freqobs = round(freqs.subset$freqobs.t/freqs.subset$n)

freq.dist = freqs.subset[,c(1:6,ind1:ind2)]
freq.dist
sum(freq.dist$freqobs)

freq.dist[,c(1:6,10)]

sum(freq.dist[,10])







