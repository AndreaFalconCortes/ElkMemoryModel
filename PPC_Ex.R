load("~/data_Ex.RData")
load("~/fitposM_Ex.RData")

mprow = vector("numeric", length = 250) #probability vectors
mpchain = vector ("numeric", length = 3)

for (i in 1:length(mprow)){
  mprow[i] = 1/length(mprow)      #joint posterior per sample
}

for (j in 1:length(mpchain)){
  mpchain[j] = 1/length(mpchain)  #joint posterior per chain
}

idpr = 1:length(mprow) #index to choose a random sample from posterior
idpc = 1:length(mpchain)

pos = fitpos[[2]] #posterior
  
datos <- data #data from PPC
  
  z = datos$z
  np = length(datos$areas)
  nsteps = length(datos$z)
  areas = datos$areas
  xv = datos$xv
  x = datos$x
  u = datos$u
  
  area = log(areas)
  idx = 1:np
  
  samples = 100 #number of samples for PPC
  lmem = vector("list", samples)
  
  for (h in 1:samples){
    
    prow = sample (idpr, size = 1, replace = TRUE, prob = mprow) #choose a random row
    pchain = sample (idpc, size = 1, replace = TRUE, prob = mpchain) #choose a random chain
    
    a = pos[prow,pchain,1] #joint posterior 
    b = pos[prow,pchain,2]
    lmbd = pos[prow,pchain,3]
    kp = pos[prow,pchain,4]
    q = pos[prow,pchain,6]
    
    zs = vector ("numeric", length(z)) #simulate the data corresponding to posterior parameters
    ds = vector("numeric", np)
    d = vector("numeric", np)
    co = vector("numeric", np)
    k = vector("numeric", np)
    
    xs = vector("numeric", np)
    mem = vector("numeric", np)
    
    zs[1] = z[1]
    
    co <- plogis(lmbd*area+kp)
    
    for (m in 2:length(zs)){
      
      ds = sqrt ( (xv[zs[m-1],1] - xv[,1])^2 + (xv[zs[m-1],2] - xv[,2])^2 )
      
      xs[zs[m-1]] = xs[zs[m-1]] + 1
      mem = xs/(m-1)
      d = exp(-(ds/a)^b)
      
      rnd = runif(1,0,1)
      
      if (rnd < q){
        k = (d*co*mem)/(sum(d*co*mem))
      }
      else{
        k = (d*co)/(sum(d*co))
      }
      
      n = sample (idx, size = 1, replace = TRUE, prob = k)
      zs[m] = n
      
    }
    
    
    us = numeric(length(zs)) 
    for(j in 1:length(zs)){
      us[j] = length(unique(zs[1:j])) 
    }
    
    lmem[[h]] = us
    
  }
  
  library("coda")
  
  t = 1:nsteps #plot PPC
  
  ix = sort(t,index.return=TRUE)$ix
  
  xsup = vector("numeric", length(us))
  xinf = vector("numeric", length(us))
  
  MA = matrix(nrow = length(us), ncol = samples)
  
  for (k1 in 1:length(us)){
    for (j1 in 1:samples){
      MA[k1,j1] = lmem[[j1]][k1]
    }
  }
  
  for (k2 in 1:length(xsup)){
    xaux = HPDinterval(as.mcmc(MA[k2,]), prob = 0.95)
    xinf[k2] = xaux[1]
    xsup[k2] = xaux[2]
  }
  if (max(xsup)>max(u)){
    ylimsup = max(xsup)
  } else {ylimsup = max(u)}
  
  plot (u, col = "blue", ylim = c(0,ylimsup), xlab = "t", ylab = "uvs", main = paste("Model II"), pch = 1, cex = 1.0)
  
  for (h1 in 1:samples){
    points(lmem[[h1]], col = grey(0.8), pch = 5, cex = 0.5)
  }
  polygon(c(rev(t[ix]),t[ix]), c(rev(xsup[ix]),xinf[ix]), col = grey(0.4), border = NA)
  points(u, col = "blue", pch = 1, cex = 1.0)
  legend(1,ylimsup, legend = c("PPC", "Data"), col = c(grey(0.8), "blue"), pch = c(5,1), cex = c(0.95,0.95))
