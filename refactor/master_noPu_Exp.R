############### RV144 CC ###################
set.seed(786)
T = as.integer(40000);
ttt=2000
SS=1;
source("for_CC_Exp.R")
# dyn.load("updategammak_noPu_Exp.so")
# dyn.load("updatealphau_noPu_Exp.so")
# dyn.load("updatealphas_Exp.so")

######################################## parameter tuning #############################################################
for (tt in 2:T) {
       # update alphau
       res2 <- .Call("updatealphau_noPu_Exp", alphaut = alpha_u[tt-1,],n_s = n_s,n_u=n_u, I=I, K=K, lambda_u = lambda_u, var_p = varp_u, ttt = ttt,gammat =gamma[,,tt-1])
       alpha_u[tt,] = res2$alphau_tt;
       A_alphau[,tt] = res2$Aalphau;

       #update gamma
       res1 <- .Call("updategammak_noPu",n_s = n_s,n_u=n_u,gammat = gamma[,,tt-1],I=I,K=K,SS = SS,alphau = alpha_u[tt,],alphas = alpha_s[tt-1,],alpha=1,mk=mk,Istar = Istar,
                    mKstar = mKstar,pp=pp, pb1 = pb1, pb2 = pb2, indi=indi)
       gamma[,,tt] = res1$gamma_tt;
       A_gm[,tt] = res1$Ag;
       Istar = res1$xIstar;
       mk = res1$xmk;
       mKstar = res1$mKstar;

       # update alphas
       res3 <- .Call("updatealphas_Exp", alphast = alpha_s[tt-1,], n_s = n_s,  K=K, I=I, lambda_s = lambda_s, gammat =gamma[,,tt], var_1 = varp_s1,var_2 = varp_s2,p_var = pvar_s,ttt = ttt)
       alpha_s[tt,] = res3$alphas_tt;
       A_alphas[,tt] = res3$Aalphas;
       cat("iteration = ", iter <- tt, "\n")
       ####
       if (tt>4000&&(tt %% 4000)==0) {
          intr = (tt-1000+1):tt
          for (kk in 1:K) {
             Au = mean(A_alphau[kk,intr])
             if (Au<0.001) {varp_u[kk] = varp_u[kk]*sqrt(0.1);}
             else if (Au <0.05) {varp_u[kk] = varp_u[kk]*sqrt(0.5);}
             else if (Au < 0.2) {varp_u[kk] = varp_u[kk]*sqrt(0.9);}
             else if (Au>0.5) {varp_u[kk] = varp_u[kk]*sqrt(1.1);}
             else if (Au>0.75) {varp_u[kk] = varp_u[kk]*sqrt(2);}
             else if (Au>0.95) {varp_u[kk] = varp_u[kk]*sqrt(10);}
             As = mean(A_alphas[kk,intr])
             if (As<0.001) {varp_s1[kk] = varp_s1[kk]*sqrt(0.1);}
             else if (As <0.05) {varp_s1[kk] = varp_s1[kk]*sqrt(0.5);}
             else if (As < 0.2) {varp_s1[kk] = varp_s1[kk]*sqrt(0.9);}
             else if (As>0.5) {varp_s2[kk] = varp_s2[kk]*sqrt(1.1);}
             else if (As>0.75) {varp_s2[kk] = varp_s2[kk]*sqrt(2);}
             else if (As>0.95) {varp_s2[kk] = varp_s2[kk]*sqrt(10);}
          }
          for ( i in 1:I) {      
             Agm = mean(A_gm[i,intr]);
             if (Agm<0.001) {pp[i] = min(0.9,pp[i]*1.4);}
             else if (Agm<0.05) {pp[i] = min(0.9,pp[i]*1.2);}
             else if (Agm<0.2) {pp[i] = min(0.9,pp[i]*1.1);}
             else if (Agm>0.6) {pp[i] = max(0.1,pp[i]*0.8);}
             else if (Agm>0.75) {pp[i] = max(0.1,pp[i]*0.5);}
             else if (Agm>0.95) {pp[i] = max(0.1,pp[i]*0.2);}
          }
        }
}
alpha_u[1,] = alpha_u[T,]
gamma[,,1] = gamma[,,T]
alpha_s[1,] = alpha_s[T,]
################################################## RUN MCMC ###########################################################
sTT=8;
for (stt in 1:sTT) {
   for (tt in 2:T) {
       # update alphau
       res2 <- .Call("updatealphau_noPu_Exp", alphaut = alpha_u[tt-1,],n_s = n_s,n_u=n_u, I=I, K=K, lambda_u = lambda_u, var_p = varp_u, ttt = ttt,gammat =gamma[,,tt-1])
       alpha_u[tt,] = res2$alphau_tt;
       A_alphau[,tt] = res2$Aalphau;
 
       #update gamma
       res1 <- .Call("updategammak_noPu",n_s = n_s,n_u=n_u,gammat = gamma[,,tt-1],I=I,K=K,SS = SS,alphau = alpha_u[tt,],alphas = alpha_s[tt-1,],alpha=1,mk=mk,Istar = Istar,
                    mKstar = mKstar,pp=pp, pb1 = pb1, pb2 = pb2, indi=indi)
       gamma[,,tt] = res1$gamma_tt;
       A_gm[,tt] = res1$Ag;
       Istar = res1$xIstar;
       mk = res1$xmk;
       mKstar = res1$mKstar;

       # update alphas
       res3 <- .Call("updatealphas_Exp", alphast = alpha_s[tt-1,], n_s = n_s,  K=K, I=I, lambda_s = lambda_s, gammat =gamma[,,tt], var_1 = varp_s1,var_2 = varp_s2,p_var = pvar_s,ttt = ttt)
       alpha_s[tt,] = res3$alphas_tt;
       A_alphas[,tt] = res3$Aalphas;
       cat("iteration = ", iter <- tt, "\n")
    }
       if (stt == sTT) {break}
       alpha_u[1,] = alpha_u[T,];
       gamma[,,1] = gamma[,,T];
       alpha_s[1,] = alpha_s[T,];
}

######################################
Tburn=0;
Mgamma = mat.or.vec(I,K);
Tseq = seq(Tburn+1,T,by=1)
for (ttt in Tseq) {
    Mgamma = Mgamma + gamma[,,ttt]; #thining
}
Mgamma = Mgamma/(T-Tburn);

