################### initialization ########################
#original data: yu, ys, # ordered data: y_u, y_s
library(e1071); 
#library(utils);
library(rBeta2009);

##### delete outlier #49
if (FALSE) { # for ENV-1 only
  y_s[[49]]=NULL
  y_u[[49]]=NULL
  n_s = n_s[-49,]
  n_u = n_u[-49,]
  sub_u = sub_u[-49]
  placebo = placebo[-49]
  I=I-1;
  
  sub_s = sub_s[-49]
  N_s = N_s[-49]
  N_u = N_u[-49]
  
  #32
  y_s[[32]]=NULL
  y_u[[32]]=NULL
  n_s = n_s[-32,]
  n_u = n_u[-32,]
  sub_u = sub_u[-32]
  placebo = placebo[-32]
  I=I-1;
  
  sub_s = sub_s[-32]
  N_s = N_s[-32]
  N_u = N_u[-32]
}
######delete outlier #16 for cast-control####
if(FALSE) {
  y_s[[16]]=NULL
  y_u[[16]]=NULL
  n_s = n_s[-16,]
  n_u = n_u[-16,]
  sub_u = sub_u[-16]
  placebo = placebo[-16]
  I=I-1;
  sub_s = sub_s[-16]
  N_s = N_s[-16]
  N_u = N_u[-16]
}
######fix some gamma's #####
indi = array(1,dim=c(I,K)) # 0 indicate that gamma_ik=0
for (k in 1:K1) {
  l2 = which((n_s[,k]/N_s)-(n_u[,k]/N_u)<=0)
  indi[l2,k]=0;
}
indi[,K] = rowSums(indi[,1:K1])
indi = matrix(as.integer(indi), nrow=I)

#####################find empirical ground mean and variance ##########
m_s = array(0, dim = c(M,1)); #hyper prior mean for mu_s
m_u = array(0, dim = c(M,1)); 
grand_m =m_u; grand_sd = m_u;

for ( mm in 1:M) {
  ys_mm = NULL; yu_mm=NULL;
  for ( i in 1:I) {
    ys_mm = c(ys_mm,y_s[[i]][,mm])
    yu_mm = c(yu_mm,y_u[[i]][,mm])
  }
  ys_mm = ys_mm[-which(ys_mm==0)]
  yu_mm = yu_mm[-which(yu_mm==0)]
  m_s[mm] = mean(ys_mm);
  m_u[mm] = mean(yu_mm);
  grand_m[mm] = mean(c(ys_mm,yu_mm))
  grand_sd[mm] =mad(c(ys_mm,yu_mm))
}
m_s = m_u + 500;
mu_s = array(0, dim=c(T,M)); mu_s[1,] = m_s;
mu_u = array(0, dim=c(T,M)); mu_u[1,] = m_u;

Sigma_s = grand_sd;#hyper prior std for mu_s
Sigma_u = grand_sd;


#############################################################
mk = array(0, dim = c(1,K-1));
Istar = 0;
mKstar = 0;

a =1; # prior for gamma
b =1;

#lambda = lambda_true;
gamma = array(0, dim=c(I,K,T));

alpha_u = array(0, dim=c(T,K));
alpha_s = array(0, dim=c(T,K));

varp_s1 = array(sqrt(5),dim=c(K,1)); # sqrt(var) for alpha_s
varp_s2 = array(sqrt(10),dim=c(K,1)); # sqrt(var)
varp_s2[K] = sqrt(20);
varp_s1[K] = sqrt(10);

pvar_s = array(0.8,dim=c(K,1));
pvar_s[K] = 0.6;
varp_u = array(sqrt(8),dim=c(K,1)); #for propose alpha_u

pp = array(0.65, dim=c(I,1));
pb1 = 0.25;
pb2 = 0.75;

lambda_s=array(100000, dim=c(1,K)); lambda_s[K] =150000; #upper bound for alpha_s
lambda_u=lambda_s;

alpha_u[1,1:(K-1)] =10; alpha_u[1,K] = 100; 
alpha_s[1,1:(K-1)] =10; alpha_s[1,K] = 150;

var_alphas = array(0,dim=c(K,K));
diag(var_alphas)=50;
chol_vars = chol(var_alphas);
inv_cvars = t(solve(chol_vars));


pvar_mus = array(0.6, dim=c(M,1));
sqrt_mus1 = array(sqrt(300),dim=c(M,1));
sqrt_mus2 = array(sqrt(500),dim=c(M,1));

pvar_muu = array(0.6, dim=c(M,1));
sqrt_muu1 = array(sqrt(300),dim=c(M,1));
sqrt_muu2 = array(sqrt(500),dim=c(M,1));

### estimate sigma^2_ik and lambda######
sig2ik = array(0, dim = c(I,K,M)); #hyper prior mean for mu_s
meanuik= sig2ik
dk=1;
for ( i in 1:I){
  if (n_s[i,1]>0 && n_u[i,1]>0) {
    yy = c(y_u[[i]][1:n_u[i,1],dk], y_s[[i]][1:n_s[i,1],dk])
    sig2ik[i,1,dk] = mad(yy)
  }else if (n_s[i,1]==0 && n_u[i,1]>1) {
    yy =  y_u[[i]][1:n_u[i,1],dk]
    sig2ik[i,1,dk] = mad(yy)
  } else if (n_s[i,1]>1 && n_u[i,1]==0) {
    yy = y_s[[i]][1:n_s[i,1],dk]
    sig2ik[i,1,dk] = mad(yy)
  } 
  if (n_u[i,1]>1) {meanuik[i,1,dk] = mean(y_u[[i]][1:n_u[i,1],dk])}
  else if (n_u[i,1]==1) {meanuik[i,1,dk] = (y_u[[i]][1:n_u[i,1],dk])}
}

for ( kk in 2:6) {
  dk = d[kk,1:M]; dk = which(dk==1)
  for ( i in 1:I) {
    if (n_s[i,kk]>0 && n_u[i,kk]>0) {
      yy = c(y_u[[i]][(sum(n_u[i,1:(kk-1)])+1):sum(n_u[i,1:kk]),dk], y_s[[i]][(sum(n_s[i,1:(kk-1)])+1):sum(n_s[i,1:kk]),dk])
      sig2ik[i,kk,dk] = mad(yy)        
    }else if (n_s[i,kk]==0 && n_u[i,kk]>1) {
      yy = y_u[[i]][(sum(n_u[i,1:(kk-1)])+1):sum(n_u[i,1:kk]),dk]
      sig2ik[i,kk,dk] = mad(yy)
    }else if (n_s[i,kk]>1 && n_u[i,kk]==0) {
      yy =  y_s[[i]][(sum(n_s[i,1:(kk-1)])+1):sum(n_s[i,1:kk]),dk]
      sig2ik[i,kk,dk] = mad(yy)
    }
    if (n_u[i,kk]>1){ meanuik[i,kk,dk] = mean(y_u[[i]][(sum(n_u[i,1:(kk-1)])+1):sum(n_u[i,1:kk]),dk])}
    else if (n_u[i,kk]==1) {meanuik[i,kk,dk] =y_u[[i]][(sum(n_u[i,1:(kk-1)])+1):sum(n_u[i,1:kk]),dk]} 
  }    
}

for ( kk in 7:K1) {
  dk = d[kk,1:M]; dk = which(dk==1)
  for ( i in 1:I) {
    if (n_s[i,kk]>0 && n_u[i,kk]>0) {
      yy = rbind(y_u[[i]][(sum(n_u[i,1:(kk-1)])+1):sum(n_u[i,1:kk]),dk], y_s[[i]][(sum(n_s[i,1:(kk-1)])+1):sum(n_s[i,1:kk]),dk])
      sig2ik[i,kk,dk] = apply(yy,2,mad)
    }else if (n_s[i,kk]==0 && n_u[i,kk]>1 ) {   
      yy = y_u[[i]][(sum(n_u[i,1:(kk-1)])+1):sum(n_u[i,1:kk]),dk]
      sig2ik[i,kk,dk] = apply(yy,2,mad)
    } else if (n_s[i,kk]>1 && n_u[i,kk]==0) {
      yy = y_s[[i]][(sum(n_s[i,1:(kk-1)])+1):sum(n_s[i,1:kk]),dk]
      sig2ik[i,kk,dk] = apply(yy,2,mad)
    }
    if (n_u[i,kk]>1) {meanuik[i,kk,dk] = apply(y_u[[i]][(sum(n_u[i,1:(kk-1)])+1):sum(n_u[i,1:kk]),dk],2,mean)}
    else if(n_u[i,kk]==1) {meanuik[i,kk,dk] =y_u[[i]][(sum(n_u[i,1:(kk-1)])+1):sum(n_u[i,1:kk]),dk]} 
  }    
}

mean_sig = array(0, dim=c(M,1))
sig_sig = array(0, dim=c(M,1))
beta1 = array(0,dim=c(M,1))
alpha1 = array(0,dim=c(M,1))
temp_var = beta1;

for ( mm in 1:M) {
  pm = which(d[,mm]==1)
  temp = as.vector(sig2ik[,pm,mm])
  temp = temp[-which(temp==0)]
  mean_sig[mm] = mean(temp)
  sig_sig[mm] = mad(temp)
  alpha1[mm] = (mean_sig[mm])^2/(sig_sig[mm])^2+2
  beta1[mm] = (mean_sig[mm])^3/(sig_sig[mm])^2+mean_sig[mm]
  
  temp1 = as.vector(meanuik[,pm,mm])
  temp1 = temp1[-which(temp1==0)]
  temp_var[mm] = sum((temp1-m_u[mm])^2)/(length(temp1)-1)
}
sig_alpha1 = mad(alpha1)*15
sig_beta1 = mad(beta1)*3

lambda = mean(temp_var/(mean_sig)^2)

alpha = array(5, dim = c(T,1))
beta = array(5, dim = c(T,1))
alpha1=mean(alpha1);
beta1 = mean(beta1);
alpha[1] = alpha1;
beta[1]=beta1;

pvar_alpha = 0.6;
sqrt_alpha1 = sig_alpha1/3;
sqrt_alpha2 = sig_alpha1*3;

pvar_beta = 0.6;
sqrt_beta1 =sig_beta1/3;
sqrt_beta2 = sig_beta1*3;
#################### acceptance rate ###########################
A_gm = array(0, dim=c(I,T)); A_gm=matrix(as.integer(A_gm),nrow=I)
A_alphau = array(0, dim=c(K,T)); A_alphau=matrix(as.integer(A_alphau),nrow=K)
A_alphas = array(0, dim=c(K,T));A_alphas=matrix(as.integer(A_alphas),nrow=K)
A_mus = array(0,dim=c(M,T)); A_mus=matrix(as.integer(A_mus),nrow=M)
A_muu = array(0,dim=c(M,T));A_muu=matrix(as.integer(A_muu),nrow=M)
A_alpha = array(0,dim=c(1,T));A_alpha=as.integer(A_alpha)
A_beta = array(0,dim=c(1,T));A_beta=as.integer(A_beta)

#################### Calculate ybar_sik and ybar_uik #############
#if (FALSE) {#
ybar_s = array(0,dim=c(I,K,M))
ybar_u = array(0,dim=c(I,K,M))
ys2_s = array(0,dim=c(I,K,M))
ys2_u = array(0,dim=c(I,K,M))
for ( i in 1:I) { 
  countu=1; counts=1;
  yui = y_u[[i]]; ysi = y_s[[i]];
  for ( k in 1:K1) {
    ldk = d[k,M+1];
    placed = which(d[k,1:M] ==1);
    if(n_u[i,k]>0) {
      
      if(n_u[i,k]>1) { 
        ybar_u[i,k,] = colMeans(yui[countu:(countu+n_u[i,k]-1),]);
      }else {ybar_u[i,k,]=yui[countu,];}
      
      
      ys2_uik = array(0,dim=c(1,ldk));
      for ( c in countu:(countu+n_u[i,k]-1)) {
        ys2_uik = ys2_uik+ (yui[c,placed]-ybar_u[i,k,placed])^2/n_u[i,k];
      }
      ys2_u[i,k,placed] = ys2_uik;
    }
    if (n_s[i,k]>0 ) {
      if ( n_s[i,k]>1) {
        ybar_s[i,k,] = colMeans(ysi[counts:(counts+n_s[i,k]-1),]);
      }else{ybar_s[i,k,] = ysi[counts,];}
      ys2_sik = array(0,dim=c(1,ldk));
      for ( c in counts:(counts+n_s[i,k]-1)) {
        ys2_sik = ys2_sik+ (ysi[c,placed]-ybar_s[i,k,placed])^2/n_s[i,k];
      }
      ys2_s[i,k,placed] = ys2_sik;
    }
    countu = countu+n_u[i,k];counts = counts+n_s[i,k];
  }
}


#} #
#########################################################################################
############### Calculate ybar_sik and ybar_uik under log transform and standardization ##########
if (FALSE){ #
  K1=K-1;
  mean_log_s = array(0,dim=c(1,M)); sd_log_s = array(0,dim=c(1,M));
  mean_log_u = array(0,dim=c(1,M)); sd_log_u = array(0,dim=c(1,M));
  for ( mm in 1:M) {
    ys_mm = NULL; yu_mm=NULL;
    for ( i in 1:I) {
      ys_mm = c(ys_mm,y_s[[i]][,mm])
      yu_mm = c(yu_mm,y_u[[i]][,mm])
    }
    ys_mm = ys_mm[-which(ys_mm==0)]
    yu_mm = yu_mm[-which(yu_mm==0)]
    log_ys = log(ys_mm); 
    mean_log_s[mm] = mean(log_ys);
    sd_log_s[mm] = sd(log_ys);
    
    log_yu = log(yu_mm); 
    mean_log_u[mm] = mean(log_yu);
    sd_log_u[mm] = sd(log_yu);
  }
  
  ybar_s = array(0,dim=c(I,K,M))
  ybar_u = array(0,dim=c(I,K,M))
  ys2_s = array(0,dim=c(I,K,M))
  ys2_u = array(0,dim=c(I,K,M))
  for ( i in 1:I) { 
    countu=1; counts=1;
    yui = y_u[[i]]; ysi = y_s[[i]];
    yui = (log(yui)-repmat(mean_log_u,sum(n_u[i,1:K1]),1))/repmat(sd_log_u,sum(n_u[i,1:K1]),1)
    ysi = (log(ysi)-repmat(mean_log_s,sum(n_s[i,1:K1]),1))/repmat(sd_log_s,sum(n_s[i,1:K1]),1)
    for ( k in 1:K1) {
      ldk = d[k,M+1];
      placed = which(d[k,1:M] ==1);
      if(n_u[i,k]>0) {
        
        if(n_u[i,k]>1) { 
          ybar_u[i,k,] = colMeans(yui[countu:(countu+n_u[i,k]-1),]);
        }else {ybar_u[i,k,]=yui[countu:(countu+n_u[i,k]-1),];}
        
        
        ys2_uik = array(0,dim=c(1,ldk));
        for ( c in countu:(countu+n_u[i,k]-1)) {
          ys2_uik = ys2_uik+ (yui[c,placed]-ybar_u[i,k,placed])^2/n_u[i,k];
        }
        ys2_u[i,k,placed] = ys2_uik;
      }
      if (n_s[i,k]>0 ) {
        if ( n_s[i,k]>1) {
          ybar_s[i,k,] = colMeans(ysi[counts:(counts+n_s[i,k]-1),]);
        }else{ybar_s[i,k,] = ysi[counts:(counts+n_s[i,k]-1),];}
        ys2_sik = array(0,dim=c(1,ldk));
        for ( c in counts:(counts+n_s[i,k]-1)) {
          ys2_sik = ys2_sik+ (ysi[c,placed]-ybar_s[i,k,placed])^2/n_s[i,k];
        }
        ys2_s[i,k,placed] = ys2_sik;
      }
      countu = countu+n_u[i,k];counts = counts+n_s[i,k];
    }
  }
  
  m_s = array(0, dim = c(M,1)); #hyper prior mean for mu_s
  m_u = array(0, dim = c(M,1)); 
  
  mu_s = array(0, dim=c(T,M)); mu_s[1,] = m_s;
  mu_u = array(0, dim=c(T,M)); mu_u[1,] = m_u;
}#


