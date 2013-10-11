# Load the data
res <-readRDS(system.file(package="MIMOSA","data/res.rds")) # only CD4+
meta <-readRDS(system.file(package="MIMOSA","data/meta.rds"))
cd4_counts <- readRDS(system.file(package="MIMOSA","data/cd4_counts.rds"))
data<-SingleCellData(res,cd4_counts,meta)
rm(res,meta,cd4_counts);gc(reset=TRUE)

#### delete bad ptids
#Filter 
del = c(523845,525881,600929,616848,631338,631772,644891,646477,704809,719186,849938,104173,125948,132409,134551,210617,211497,217402,223868,233547,323142,330034,330272,339553)

#### case-control
orig <- readRDS(system.file(package="MIMOSA","data/PTID_map_case_control.rds"))
rPTID =  unique(orig)
colnames(rPTID) = c("PTID", "PTID_orig")
metadata(data)<-merge(metadata(data),rPTID,by.x="PTID",by.y="PTID")

keep<-subset(metadata(data),!PTID_orig%in%del,select="name")[[1]]
data<-data[(keep)]

#get paired data
#########################################################
#########################################################
sub_negctrl_800 = which(meta$Stim =="negctrl 1" & meta$VISITNO == 800) # unstim & VISITNO

sub_stim_800= which(meta$Stim=="92TH023 Env" & meta$VISITNO == 800) # stim & VISITNO

subsub_800 = intersect(meta$PTID[sub_negctrl_800],meta$PTID[sub_stim_800])

sub_s=NULL; sub_u=NULL;
for ( i in 1:length(subsub_800)) {
  sub_s =c(sub_s,which(meta$PTID[sub_stim_800]==(subsub_800[i])))
  sub_u = c(sub_u, which(meta$PTID[sub_negctrl_800]==(subsub_800[i])))
}

sub_s = sub_stim_800[sub_s]; #index for the stim and unstim sample
sub_u = sub_negctrl_800[sub_u];

#cd4
N_s = as.integer(cd4_counts[sub_s])
N_u = as.integer(cd4_counts[sub_u])

ys= vector(mode="list", length=length(sub_s));
yu =vector(mode="list", length=length(sub_u));
for (i in 1:length(sub_s)) {
  ys[[i]] = res[[sub_s[i]]]
  yu[[i]] = res[[sub_u[i]]]  
}
placebo = meta$vaccine[sub_s]
#rm(sub_negctrl_800,subsub_800,sub_stim_800)
###############################################################
library(caTools)

I = length(yu);
M = dim(yu[[1]])[2]
K = 2^M;
K1=K-1;
y_u=yu; y_s = ys;
n_s = array(as.integer(0),dim=c(I,K));
n_u = array(as.integer(0),dim=c(I,K));


d = array(as.integer(0),dim=c(K,M))
count=1;
for (kk in 1:(M-1)) {
  set = combs(1:M,kk);
  ll = dim(combs(1:M,kk)); 
  
  for ( lll in 1:ll[1]) {
    d[(count-1+lll),set[lll,]] = as.integer(1);
  }
  
  count = count+choose(M,kk);
}
d[count,]=as.integer(1)
d = cbind(d,as.integer(rowSums(d))); # last column the number of 1's

for (i in 1:I) {
  count_s = 1; count_u = 1;
  yui = 1*(yu[[i]]!=0)
  ysi = 1*(ys[[i]]!=0)
  for ( kk in 1:K1) {   
    if( length(yui)>0) {
      sel = apply(yui,1,function(xrow)(identical(as.numeric(xrow),as.numeric(d[kk,1:M])))) 
      n_u[i,kk] = as.integer(sum(sel))
      if(n_u[i,kk]>0) {
        y_u[[i]][count_u:(count_u+n_u[i,kk]-1),] = yu[[i]][sel,]
        count_u = count_u+n_u[i,kk]
      }
    }
    if(length(ysi)>0) {
      sel = apply(ysi,1,function(xrow)(identical(as.numeric(xrow),as.numeric(d[kk,1:M])))) 
      
      n_s[i,kk] = as.integer(sum(sel))
      if(n_s[i,kk]>0) {
        y_s[[i]][count_s:(count_s+n_s[i,kk]-1),] = ys[[i]][sel,]
        count_s = count_s+n_s[i,kk]
      }
    }
  }
  n_s[i,K] = as.integer(N_s[i]-sum(n_s[i,1:K1]))
  n_u[i,K] = as.integer(N_u[i]-sum(n_u[i,1:K1]))
}

rm(ys,yu,count_s,count_u,count)
#######take off zero columns #############
nu0 = which(colSums(n_u)==0)
ns0 = which(colSums(n_s)==0)
n0 = intersect(nu0,ns0)

d = d[-n0,]
n_s = n_s[,-n0]
n_u = n_u[,-n0]

K = K-length(n0)

rm(nu0,ns0,n0)
#### additional delete #####
#del = c(7,17,21,22,23,24,26,29,34,35,37,38,39,40,42,43,46,47);
del = which(colSums(n_s)<=12);

for ( i in 1:I) {
  cs = cumsum(n_s[i,]); cu = cumsum(n_u[i,]); posis = NULL; posiu = NULL;
  for (j in del) {
    if(n_s[i,j]>1) {
      posis = c(posis, (cs[j-1]+1):cs[j])
    } else if (n_s[i,j]==1) {
      posis = c(posis,cs[j])
    }
    if(n_u[i,j]>1) {
      posiu = c(posiu, (cu[j-1]+1):cu[j])
    } else if (n_u[i,j] ==1) {
      posiu = c(posiu,cu[j])
    } 
  }
  if (length(posis) >0) {y_s[[i]] = y_s[[i]][-posis,]}
  if (length(posiu) >0) {y_u[[i]] = y_u[[i]][-posiu,]}
}
d = d[-del,]
n_s = n_s[,-del]
n_u = n_u[,-del]

K = dim(n_s)[2]

K1=K-1;
n_s[,K] = N_s-as.integer(rowSums(n_s[,1:K1]))
n_u[,K] = N_u-as.integer(rowSums(n_u[,1:K1]))

###########################################################
###########################################################
library(e1071); 
######fix some gamma's #####
indi = array(1,dim=c(I,K)) # 0 indicate that gamma_ik=0
for (k in 1:K1) {
  l2 = which((n_s[,k]/N_s)-(n_u[,k]/N_u)<=0)
  indi[l2,k]=0;
}
indi[,K] = rowSums(indi[,1:K1])
indi = matrix(as.integer(indi), nrow=I)

#############################################
SS=as.integer(1);
a = 1; b = 1;
mk = array(as.integer(0), dim = c(1,K1));
Istar = 0;
mKstar = 0;

gamma =array(as.integer(0), dim=c(I,K,T));

alpha_u = array(0, dim=c(T,K));
alpha_s = array(0, dim=c(T,K));

varp_s1 = array(sqrt(3),dim=c(K,1)); # sqrt(var)
varp_s2 = array(sqrt(10),dim=c(K,1)); # sqrt(var)
varp_s2[K] = sqrt(20);
varp_s1[K] = sqrt(10);

pvar_s = array(0.8,dim=c(K,1));
pvar_s[K] = 0.6;
varp_u = array(sqrt(8),dim=c(K,1));

pp = array(0.65, dim=c(I,1));
pb1 = 0.1;
pb2 = 0.6;

lambda_s = rep(0,K);  
lambda_s[1:K1] = (10^-2)*max(N_s,N_u)
lambda_s[K] = max(N_s,N_u)-sum(lambda_s[1:K1])
lambda_u = lambda_s

alpha_u[1,1:(K-1)] =10; alpha_u[1,K] = 150; 
alpha_s[1,1:(K-1)] =10; alpha_s[1,K] = 100;

#################### acceptance rate ###########################
A_gm = array(as.integer(0), dim=c(I,T));
A_alphau = array(as.integer(0), dim=c(K,T));
A_alphas = array(as.integer(0), dim=c(K,T));


