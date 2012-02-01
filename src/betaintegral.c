#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>

/* Compute the integral of the joint distribution of ps~Be(aS,bS), pu~Be(a0,b0), where ps>pu.*/
//TODO monte carlo integration. Will be faster than these loops.
//This is now deprecated in favor of Monte-Carlo integration. Will eventually move everything to C++
void betaintegral_c(double* alphaS,double* betaS,double* alpha0, double* beta0,int* Nu,int* nu,int* Ns,int* ns, double* pu, double* dpu, double* res, int* l, int *p)
{
 for(int i=0;i<*l;i++){
  res[i]=lbeta((double) ns[i]+*alphaS,*betaS+(double) Ns[i])+lbeta((double) nu[i]+*alpha0,*beta0+(double) Nu[i]);
  double acc=0;
  for(int j=0;j<*p;j++){
    acc=acc+dbeta(pu[j],*alpha0+(double) nu[i],*beta0+(double) Nu[i],0)*(1-pbeta(pu[j],*alphaS+(double) ns[i],*betaS+(double) Ns[i],0,0))* dpu[0];
  }
  res[i]=res[i]+log(acc);
  //Rprintf("%f\n",log(acc));
//  if(fpclassify(res[i])==FP_INFINITE)
//	res[i]=-DBL_MAX;
// }
 }
}


