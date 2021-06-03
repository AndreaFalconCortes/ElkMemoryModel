Sys.setenv(USE_CXX14 = 1) 
library("rstan") # observe startup messages

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

cat(file = "ProbPatchMem.stan", "

functions{
    
    //This functions return the probability vector k (log) when animal use its memory (DAM_lpmf) or not (DA_lpmf)
    
    real DA_lpmf (int id, int np, vector distance, vector areas, real a, real b, real lambda, real ka){
    
        vector [np] d;
        vector [np] co;
        vector [np] k;

        for (i in 1:np){        
            d[i] = exp(-(distance[i]/a)^b); 
        } 
        
        co = inv_logit(lambda*areas + ka);
        k = (d .*co)/sum(d .*co);
        
        return (log(k[id]));
        
    }
    
    real DAM_lpmf (int id, int np, vector distance, vector m, vector areas, real a, real b, real lambda, real ka){
    
        vector [np] d;
        vector [np] co;
        vector [np] k;
        
        for (i in 1:np){        
            d[i] = exp(-(distance[i]/a)^b);
        } 
        
        co = inv_logit(lambda*areas + ka);
        k = (d .*co .*m)/sum(d .*co .*m);
        
        if (k[id] > 0){
            return (log(k[id]));
        } else {
            return negative_infinity();
        }
    }

}

data {

    int <lower = 0> np;
    int <lower = 0> nsteps; 
    vector <lower = 0> [np] areas;
    matrix [np,2] xv;
    int z[nsteps];
        
}

parameters {

    real <lower=0> a; 
    real <lower=0> b; 
    real <lower=0> lambda; 
    real <lower=0> ka; 
    simplex[2] q; 

}  

model {

    vector [np] ds;
    vector [np] xs;
    vector [np] m;
    vector [np] logareas;
    vector [2] log_q = log(q);
    
    logareas = log(areas);

    xs = rep_vector(0, np);
    m = rep_vector(0, np);
    
    //prior distributions
    
    a ~ normal(0,10);
    b ~ normal(0,1); 
    lambda ~ normal(0,1); 
    ka ~ normal(0,1); 
    q ~ beta(1,1); 
    
    //likelihood compute
    
    for (i in 1:(nsteps - 1)){
        vector[2] lps = log_q;
        
        for (k in 1:np){
    
        ds[k] = sqrt ( (xv[z[i],1] - xv[k,1])^2 + (xv[z[i],2] - xv[k,2])^2 ); 
                                                                                 
        }
        
        xs[z[i]] = xs[z[i]] + 1;
        m = xs/i;
        
        lps[1] += DA_lpmf(z[i+1] | np, ds, logareas, a, b, lambda, ka); //(log) probability to choose the new patch without memory use
        lps[2] += DAM_lpmf(z[i+1] | np, ds, m, logareas, a, b, lambda, ka); //(log) probability to choose the new patch with memory use
        target += log_sum_exp(lps); // loglikelihood
    }
    

}
")



np <- 100     #number of patches to consider
nsteps <- 500 #number of steps for the trajectory
areas <- rweibull(np, 2, 1) #patch areas' distribution 

xv = matrix(0,np,2) #matrix to save the patches' euclidean coordinates whitin a 10x10 continuous space
xv [,1] <- runif(length(xv[,1]), 0, 10) #random X coordinate
xv [,2] <- runif (length(xv[,2]), 0, 10)#random Y coordinate 
z = vector('numeric',nsteps) #vector to save the "animal" trajectory
u = vector('numeric', nsteps) #vector to save the number of uvs at time t

rnd = sample (1:np,1)
z[1] = rnd #random initial position

ds = vector("numeric", np) #vector to save the distance to each patch
xs = vector("numeric", np) #vector to save the number of visits to each patch

d = vector("numeric", np) #vector to save the probability decay with the distance
co = vector("numeric", np) #vector to save the probability increase with the patch area
m = vector("numeric", np) #vector to save the probility to revisit a patch due the linear reinforcement 
k = vector("numeric", np) # d*c or d*c*m

idx <- 1:np #patchs' index

a <- 1 # alpha -> scale parameter for the exponential function that defines the probability decay with distance
b <- 1 # beta -> shape parameter for the exponential function that defines the probability decay with distance
lambda <- 2.0 # lambda -> slope parameter for the logit function that defines the probability increase as function of patch area
ka <- 1.0 # kappa -> intercept parameter for the logit function that defines the probability increase as function of patch area
q <- 0.2 # q -> parameter that defines the memory use frequency


co <- plogis(lambda*log(areas)+ka) #logit probability increase with area

for (i in 2:length(z)){ #simulate trajectory
  
  ds <- sqrt ( (xv[z[i-1],1] - xv[,1])^2 + (xv[z[i-1],2] - xv[,2])^2 ); #distance from the actual position z[i-1] to the other patches  
  
  xs[z[i-1]] <- xs[z[i-1]] + 1 #number of visits to the patch z[i-1]
  
  m <- xs/(i-1) #probability of revisit a patch (number of visits/time)
  d <- exp(-(ds/a)^b) #exponential probability decay with distance
  
  rnd <- runif(1,0,1) #choose the probability to use or not memory
  if (rnd < q){
    k <- (d*co*m)/(sum(d*co*m)) #use memory and revisit a patch cosidering distance, area, and reinforcement 
  }
  else{
    k <- (d*co)/(sum(d*co)) #don't use memory and choose a patch considering distance and area
  }
  
  j <- sample (idx, size = 1, replace = TRUE, prob = k) #choose a patch with probability k
  
  z[i] <- j #save the step
  
}

for(j in 1:length(z)){
  u[j] = length(unique(z[1:j])) # number of uvs until time j
}

data.sim <- list(np = np, 
                 nsteps = nsteps,
                 areas = areas,
                 xv = xv,
                 z = z) #necessary data to perform the bayesian fit

data <- list(np = np, 
             nsteps = nsteps,
             areas = areas,
             xv = xv,
             z = z,
             u = u) #data to perform PPC

fit <- stan(file = 'ProbPatchMem.stan', data = data.sim,
            iter = 500, chains = 3, control = list(adapt_delta = 0.9)) #call the model to perform the bayes fit

print(fit) #show the fit
posterior <- as.array(fit) #convert object fit to array

fitpos = list (fit, posterior) #save in a list fit and posterior
save(fitpos, file = "fitposM_Ex.RData") #save in RDATA file
save(data, file = "data_Ex.RData") #save the data for compute PPC
