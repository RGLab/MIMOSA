##########################case-control CD4 ##############
set.seed(536)

source("realdata_Casecontrol.R")
T = as.integer(4000);
source("init_realdata_Casecontrol.R")
ttt=as.integer(2000);
n_ss = as.integer(n_s); n_uu = as.integer(n_u) 

ybars = apply(ybar_s, 3, function(x) as.vector(x))
ybaru = apply(ybar_u, 3, function(x) as.vector(x))
ys2s = apply(ys2_s, 3, function(x) as.vector(x))
ys2u = apply(ys2_u, 3, function(x) as.vector(x))

rm(ybar_s,ybar_u,ys2_u,ys2_s)
gammat = apply(gamma,3, function(x) as.integer(x))

################## parameter tuning ######################
ptm <- proc.time()
result<-.Call("model",T=T,I=I, K=K, M=M, ttt=ttt, SS=as.integer(1),alpha_u=alpha_u, alpha_s=alpha_s, mu_u=mu_u, mu_s=mu_s, alpha=alpha,
                      beta=beta, gamma=gammat,n_s=n_ss, n_u=n_uu,varp_u=varp_u, lambda_u=lambda_u, indi=indi, d=d, ybar_s = ybars, ybar_u=ybaru, 
                      ys2_s = ys2s, ys2_u = ys2u, a=a, b=b,lambda=lambda,mk = mk, Istar = Istar, mKstar = mKstar, pp = pp, pb1 = pb1, 
                      pb2 = pb2,lambda_s = lambda_s, var_1 = varp_s1,var_2 = varp_s2,p_var = pvar_s, p_vars = pvar_mus,var_1s = sqrt_mus1, var_2s=sqrt_mus2,m_s=m_s,Sigma_s =Sigma_s,
                      p_varu=pvar_muu,var_1u=sqrt_muu1,var_2u=sqrt_muu2,m_u=m_u, Sigma_u=Sigma_u,p_vara=pvar_alpha, var_1a=sqrt_alpha1,var_2a=sqrt_alpha2,sig_alpha1=sig_alpha1,alpha1=alpha1,
                      p_varb=pvar_beta,var_1b=sqrt_beta1,var_2b=sqrt_beta2,sig_beta1 = sig_beta1, beta1= beta1, A_alphau=A_alphau,A_alphas=A_alphas, A_gm=A_gm, A_mus=A_mus, A_muu=A_muu,  
                      A_alpha=A_alpha,A_beta=A_beta, Tune = as.integer(1), pgamma = 0.4)
proc.time() - ptm
Istar = result$Istar
mk = result$mk
mKstar = result$mKstar
var_1a = result$var_1a
var_1b = result$var_1b
var_2a = result$var_2a
var_2b = result$var_2b

alpha_u[1,] = alpha_u[T,]; alpha_s[1,] = alpha_s[T,]; gammat[,1] = gammat[,T];
mu_s[1,] = mu_s[T,]; mu_u[1,]=mu_u[T,];  alpha[1] = alpha[T]; beta[1] = beta[T]; 
sTT =10; 
#######################################################
for (stt in 1:sTT) {
    result<-.Call("model",T=T,I=I, K=K, M=M, ttt=ttt, SS=as.integer(1),alpha_u=alpha_u, alpha_s=alpha_s, mu_u=mu_u, mu_s=mu_s, alpha=alpha,
                      beta=beta, gamma=gammat,n_s=n_ss, n_u=n_uu,varp_u=varp_u, lambda_u=lambda_u, indi=indi, d=d, ybar_s = ybars, ybar_u=ybaru, 
                      ys2_s = ys2s, ys2_u = ys2u, a=a, b=b,lambda=lambda,mk = mk, Istar = Istar, mKstar = mKstar, pp = pp, pb1 = pb1, 
                      pb2 = pb2,lambda_s = lambda_s, var_1 = varp_s1,var_2 = varp_s2,p_var = pvar_s, p_vars = pvar_mus,var_1s = sqrt_mus1, var_2s=sqrt_mus2,m_s=m_s,Sigma_s =Sigma_s,
                      p_varu=pvar_muu,var_1u=sqrt_muu1,var_2u=sqrt_muu2,m_u=m_u, Sigma_u=Sigma_u,p_vara=pvar_alpha, var_1a=sqrt_alpha1,var_2a=sqrt_alpha2,sig_alpha1=sig_alpha1,alpha1=alpha1,
                      p_varb=pvar_beta,var_1b=sqrt_beta1,var_2b=sqrt_beta2,sig_beta1 = sig_beta1, beta1= beta1, A_alphau=A_alphau,A_alphas=A_alphas, A_gm=A_gm, A_mus=A_mus, A_muu=A_muu,  
                      A_alpha=A_alpha,A_beta=A_beta, Tune = as.integer(0),pgamma = 0.5)
    if (stt == sTT) {break}
    Istar = result$Istar
    mk = result$mk
    mKstar = result$mKstar
    alpha_u[1,] = alpha_u[T,]; alpha_s[1,] = alpha_s[T,]; gammat[,1] = gammat[,T];
    mu_s[1,] = mu_s[T,]; mu_u[1,]=mu_u[T,];  alpha[1] = alpha[T]; beta[1] = beta[T]; 
    cat("iteration = ", iter <- stt, "\n")
}

##########################################################
for ( tt in 1:T) {
 gamma[,,tt]=matrix(gammat[,tt], nrow=I)
}

Mgamma = apply(gamma,c(1,2),sum)/T
