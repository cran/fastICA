#include <R.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

                                         
void rowcentre(float*,int,int);
void colstandard(float*,int,int);
void mmult(float*,int,int,float*,int,int,float*);
void onefy(float*,int,float*);
void bff(float*,int,float*,int,int,float,float*,float*);
void bff2(float*,int,float*,int,int,float,float*,float*);
void gramsch(float*,int,int,int);
void rowstd(float*,int,int,int);
void anfu(float*,int,float*,int,int,float,float*);
void anfu2(float*,int,float*,int,int,float,float*);
void transpose_mat(float*, int*, int*, float*);
int min(int*, int*);
int max(int*, int*);

void F77_NAME(sgesvd)(char*, char*, int*, int*, float*, int*, float*, float*, int*, float*, int*, float*, int*, int*);
void F77_NAME(sgemm)(char*,char*,int*,int*,int*,float*,float*,int*,float*,int*,float*,float*,int*);


void rowcentre(float *ans,int n,int p){
/*  mean centres nxp matrix ans */
  double tmp;
  int i,j;
  for(i=0;i<n;i++){
    tmp=0;
    for(j=0;j<p;j++){
      tmp=tmp+((double) ans[p*i+j])/p;}
    for(j=0;j<p;j++){ans[p*i+j]-=(float) tmp;}}
}

void colstandard(float *ans,int n,int p){
  /*  transform columns of nxp matrix ans to have zero mean and unit variance */
  double tmp[2];
  double tmp1;
  int i,j;
  for(i=0;i<p;i++){
    tmp[0]=0; 
    tmp[1]=0;
    
    for(j=0;j<n;j++){
      tmp[0]+=(double) ans[p*j+i]; 
      tmp[1]+=((double) ans[p*j+i])*((double) ans[p*j+i]);
    }
    
    tmp[0]=tmp[0]/n;
    tmp1=(tmp[1]-n*(tmp[0])*(tmp[0]))/(n-1);
    
    tmp[1]=sqrt(tmp1);
    
    for(j=0;j<n;j++){
      ans[p*j+i]=(float) (( ((double)ans[p*j+i]) -tmp[0])/tmp[1]);
    }
  }
}


void svd(float *mat, int *n, int *p, float *u, float *d, float *v){
  
  /*  calculates svd decomposition of nxp matrix mat */
  /*    mat is a pointer to an nxp array of floats */
  /*    n is a pointer to an integer specifying the no. of rows of mat */
  /*    p is a pointer to an integer specifying the no. of cols of mat */
  /*    u is a pointer to a float array of dimension (n,n) */
  /*    d is a pointer to a float array of dimension min(n,p) */
  /*    v is a pointer to a float array of dimension (p,p) */
  
  
  int info,lwork,i,j;
  float *work, *mat1, *u1, *v1;
  char jobu='A',jobvt='A';
  
  i=3*min(n,p)+max(n,p);
  j=5*min(n,p);
  lwork=10*max(&i,&j);
  
  work=Calloc(lwork,float);
  mat1=Calloc((*n)*(*p),float);
  u1=Calloc((*n)*(*n),float);
  v1=Calloc((*p)*(*p),float);
  
  transpose_mat(mat, n, p, mat1);
  
  F77_CALL(sgesvd)(&jobu, &jobvt, n, p, mat1, n, d, u1, n, v1, p, work, &lwork, &info);
  
  transpose_mat(u1, n, n, u);

  transpose_mat(v1, p, p, v);


  Free(mat1);
  Free(u1);
  Free(v1);
  Free(work);

}

void transpose_mat(float *mat, int *n, int *p, float *ans){
/*    transpose nxp matrix mat */
  int i,j;
  
  for(i=0;i<*n;i++){
    for(j=0;j<*p;j++){
      *(ans+j*(*n)+i)=*(mat+i*(*p)+j);}}
}


int min(int *a, int *b){
/*  find minimum of a and b */
  int ans;

  ans=*b;
  if(*a<*b) ans=*a;
  
  return ans;
}

int max(int *a, int *b){
/*  find maximum of a and b */

  int ans;

  ans=*b;
  if(*a>*b) ans=*a;
  
  return ans;
}

void mmult(float *A,int n,int p,float *B,int q,int r,float *C){
/*    matrix multiplication using FORTRAN BLAS routine SGEMM */
/*    A is (n*p) and B is (q*r), A*B returned to C  */

  float alpha=1.0,beta=0.0;
  float *matA,*matB,*matC;
  int M,K,N;
  char transA='N',transB='N';

  if(p != q){printf("Error, matrices not suitable\nfor multiplication.\n\n");}
  else{
    M=n;
    K=p;
    N=r;
    matA=Calloc(M*K,float);
    matB=Calloc(K*N,float);
    matC=Calloc(M*N,float);

    transpose_mat(A,&M,&K,matA);
    transpose_mat(B,&K,&N,matB);
    F77_CALL(sgemm)(&transA,&transB,&M,&N,&K,&alpha,matA,&M,matB,&K,&beta,matC,&M);
    transpose_mat(matC,&N,&M,C);

    Free(matA);
    Free(matB);
    Free(matC);
  }
}

void onefy(float *ww,int e,float *tmpm){
  /* take ww, (e*e), and return to tmpm "ww1" */
  float *u,*v,*d,*temp;int i;
  u=Calloc(e*e,float);
  d=Calloc(e,float);
  v=Calloc(e*e,float);
  temp=Calloc(e*e,float);
 
  svd(ww,&e,&e,u,d,v);
  for(i=0;i<e;i++){temp[i*e+i]=1/(d[i]);}
  
  mmult(u,e,e,temp,e,e,v);
  transpose_mat(u,&e,&e,temp); 
  mmult(v,e,e,temp,e,e,u);
  mmult(u,e,e,ww,e,e,tmpm);
  
  Free(u);Free(v);Free(d);Free(temp);
}

void bff(float *ww,int e,float *ans,int f,int p,float alpha,float *wwpl,float *Tol){
  float *mat1,*mat2,*mat3,*mat4,*mat5,*mat6;
  int i,j;
  float mean;
  /* ww is W, ans is X, wwpl will take the answer matrix*/
  if(e != f){printf("error in bff, dims dont match\n");}
  else{
    mat1=Calloc(e*p,float);
    mat2=Calloc(e*p,float);
    mat3=Calloc(e*e,float);
    mat4=Calloc(e*e,float);
    mat5=Calloc(e*e,float);
    mat6=Calloc(e*e,float);
    
    mmult(ww,e,e,ans,e,p,mat1);/*mat1 is WX*/
    
    
    for(i=0;i<e;i++){for(j=0;j<p;j++){mat1[i*p+j]=tanh(alpha*mat1[i*p+j]);}}/*mat1 is GWX*/
    transpose_mat(ans,&e,&p,mat2);
    for(i=0;i<e;i++){for(j=0;j<p;j++){mat2[i*p+j]=(mat2[i*p+j])/p;}}/*mat2 is t(X)/p */
    mmult(mat1,e,p,mat2,p,e,mat3); /*mat3 is V1 */
    for(i=0;i<e;i++){for(j=0;j<p;j++){mat1[i*p+j]=(alpha*(1-(mat1[i*p+j])*(mat1[i*p+j])));}}
    /*mat1 is GWX1*/
    for(i=0;i<e;i++){mean=0;for(j=0;j<p;j++){mean+=((mat1[i*p+j])/p);}
    mat4[i*e+i]=mean;}  /*mat4 is D */
    mmult(mat4,e,e,ww,e,e,mat5);  /* mat5 is V2 */
    for(i=0;i<e;i++){for(j=0;j<e;j++){mat4[i*e+j]=(mat3[i*e+j]-mat5[i*e+j]);}}
    /* mat4 is W1 */
    transpose_mat(ww,&e,&e,mat6);
    onefy(mat4,e,wwpl); /* wwpl is W2 (sic) */
    
    
    
    
    
    mmult(wwpl,e,e,mat6,e,e,mat5); /*mat5 is C */
    mean=0;
    for(i=0;i<e;i++){
      if(fabs(1-fabs(mat5[i*e+i]))>mean){
	mean=(fabs(1-fabs(mat5[i*e+i])));}}
    *Tol=mean;
    Free(mat1);Free(mat2);Free(mat3);Free(mat4);Free(mat5);Free(mat6);
  }}

void anfu(float *ww,int e,float *ans,int f,int p,float alpha,float *wwpl){
  float *mat1,*mat2,*mat3,*mat4,*mat5;
  int i,j;
  float mean;
  /* ww is W, ans is X, wwpl will take the answer matrix*/
  if(e != f){printf("error in anfu, dims dont match\n");}
  else{
    mat1=Calloc(1*p,float);
    mat2=Calloc(e*p,float);
    mat3=Calloc(1*e,float);
    mat4=Calloc(1*e,float);
    mat5=Calloc(1*e,float);
    
    mmult(ww,1,e,ans,e,p,mat1);/*mat1 is WX*/
    
    
    for(i=0;i<p;i++){mat1[i]=tanh(alpha*mat1[i]);}/*mat1 is GWX*/
    transpose_mat(ans,&e,&p,mat2);
    for(i=0;i<e;i++){for(j=0;j<p;j++){mat2[i*p+j]=(mat2[i*p+j])/p;}}
    /*mat2 is t(X)/p */
    mmult(mat1,1,p,mat2,p,e,mat3); /*mat3 is V1 */
    for(i=0;i<p;i++){mat1[i]=(alpha*(1-(mat1[i])*(mat1[i])));}
    /*mat1 is GWX1*/
    mean=0;
    for(j=0;j<p;j++){mean+=((mat1[j])/p);}
    for(i=0;i<e;i++){mat5[i]=(ww[i])*mean;}  /* mat5 is V2 */
    for(i=0;i<e;i++){wwpl[i]=(mat3[i]-mat5[i]);}  /* wwpl is W1 */
    
    
    Free(mat1);Free(mat2);Free(mat3);Free(mat4);Free(mat5);
    
  }}

void bff2(float *ww,int e,float *ans,int f,int p,float alpha,float *wwpl,float *Tol){
  float *mat1,*mat2,*mat3,*mat4,*mat5,*mat0,*mat6;
  int i,j;
  float mean;
  /* ww is W, ans is X, wwpl will take the answer matrix*/
  if(e != f){printf("error in bff2, dims dont match\n");}
  else{
    mat0=Calloc(e*p,float);
    mat1=Calloc(e*p,float);
    mat2=Calloc(e*p,float);
    mat3=Calloc(e*e,float);
    mat4=Calloc(e*e,float);
    mat5=Calloc(e*e,float);
    mat6=Calloc(e*e,float);
    mmult(ww,e,e,ans,e,p,mat1);/*mat1 is WX*/
    for(i=0;i<e;i++){for(j=0;j<p;j++){mat0[i*p+j]=(mat1[i*p+j])*exp(-0.5*(mat1[i*p+j])*(mat1[i*p+j]));}}/*mat0 is GWX*/
    transpose_mat(ans,&e,&p,mat2);
    for(i=0;i<e;i++){for(j=0;j<p;j++){mat2[i*p+j]=(mat2[i*p+j])/p;}}/*mat2 is t(X)/p */
    mmult(mat0,e,p,mat2,p,e,mat3); /*mat3 is V1 */
    for(i=0;i<e;i++){for(j=0;j<p;j++){mat1[i*p+j]=((1-(mat1[i*p+j])*(mat1[i*p+j]))*exp(-0.5*(mat1[i*p+j])*(mat1[i*p+j])));}}
    /*mat1 is GWX1*/
    for(i=0;i<e;i++){mean=0;for(j=0;j<p;j++){mean+=((mat1[i*p+j])/p);}
    mat4[i*e+i]=mean;}  /*mat4 is D */
    mmult(mat4,e,e,ww,e,e,mat5);  /* mat5 is V2 */
    for(i=0;i<e;i++){for(j=0;j<e;j++){mat4[i*e+j]=(mat3[i*e+j]-mat5[i*e+j]);}}
    /* mat4 is W1 */
    transpose_mat(ww,&e,&e,mat6);
    onefy(mat4,e,wwpl); /* wwpl is W2 (sic) */
    
    mmult(wwpl,e,e,mat6,e,e,mat5); /*mat5 is C */
    mean=0;
    for(i=0;i<e;i++){
      if(fabs(1-fabs(mat5[i*e+i]))>mean){
	mean=(fabs(1-fabs(mat5[i*e+i])));}}
    *Tol=mean;
    Free(mat1);Free(mat2);Free(mat3);Free(mat4);Free(mat5);Free(mat0);Free(mat6);
  }}

void anfu2(float *ww,int e,float *ans,int f,int p,float alpha,float *wwpl){
  float *mat1,*mat2,*mat3,*mat4,*mat5;
  int i,j;
  float mean;
  /*ww is (a row of) W, ans is X, wwpl will take the answer vector*/
  if(e != f){printf("error in anfu2, dims dont match\n");}
  else{
    mat1=Calloc(1*p,float);
    mat2=Calloc(e*p,float);
    mat3=Calloc(1*e,float);
    mat4=Calloc(1*e,float);
    mat5=Calloc(1*e,float);
    
    mmult(ww,1,e,ans,e,p,mat1); /*mat1 is WX*/
    
    for(i=0;i<p;i++){mat1[i]=((mat1[i])*exp(-0.5*(mat1[i])*(mat1[i])));}
    /*mat1 is GWX*/
    transpose_mat(ans,&e,&p,mat2);
    for(i=0;i<e;i++){for(j=0;j<p;j++){mat2[i*p+j]=(mat2[i*p+j])/p;}}
    /*mat2 is t(X)/p */
    mmult(mat1,1,p,mat2,p,e,mat3); /*mat3 is V1 */
    
    mmult(ww,1,e,ans,e,p,mat1);/*mat1 is WX (again) */
    for(i=0;i<p;i++){mat1[i]=((1-(mat1[i])*(mat1[i]))*exp(-.5*(mat1[i])*(mat1[i])));}  /*mat1 is GWX1*/
    mean=0;
    for(j=0;j<p;j++){mean+=((mat1[j])/p);}
    for(i=0;i<e;i++){mat5[i]=(ww[i])*mean;}  /* mat5 is V2 */
    for(i=0;i<e;i++){wwpl[i]=(mat3[i]-mat5[i]);}  /* wwpl is W1 */
    
    
    Free(mat1);Free(mat2);Free(mat3);Free(mat4);Free(mat5);
    
  }}

void gramsch(float *ww,int n,int m,int k){
  int ip,jp;float tmp;
  /* do Gram-Schmidt on row k of (n*m) matrix ww */ 
  k-=1;
  if(k>n){printf("\nError in gramsch\n");}
  else{
    for(ip=0;ip<k;ip++){tmp=0;
    for(jp=0;jp<m;jp++){tmp+=((ww[m*ip+jp])*(ww[m*k+jp]));}
    for(jp=0;jp<m;jp++){ww[m*k+jp]=(ww[m*k+jp]-((ww[m*ip+jp])*tmp));}}}}

void rowstd(float *ww,int n,int m, int k){
  /* for ww (n*m), make ||ww[k, ]|| equal 1 */
  float tmp=0;int i;k-=1;
  if(k>n){printf("\nError in rowstd\n");}
  else{
    for(i=0;i<m;i++){tmp+=((ww[k*m+i])*(ww[k*m+i]));}
    tmp=sqrt(tmp);
    for(i=0;i<m;i++){ww[k*m+i]=((ww[k*m+i])/tmp);}}
}

void icainc(float *data_matrix,float *w_matrix,int *nn,int *pp,int *ee,float *alpha,int *rowflag,int *colflag,int *funflag,int *maxit,float *lim,int *defflag, int *verbose, float *ansx,float *ansk,float *answ,float *ansa,float *ansx2){ 
  
  int i,j,k,n,p,e;
  float tol;
  float *ans,*uu,*dd,*vv,*pwh,*pwhh,*tmpm,*ww,*wwpl;
  float *mat1,*mat2,*mat3,*mat4,*mat5,*mat6;
  
  n=*nn;p=*pp;e=*ee;
  
  ans=Calloc(n*p,float);
  
  for(i=0;i<n;i++){for(j=0;j<p;j++){ans[i*p+j]=data_matrix[i*p+j];}}
  
  if(*rowflag==1){  
    rowcentre(ans,n,p);
    if(*verbose==1) printf("Centering\n");
  }
  
  if(*colflag==1){  
    colstandard(ans,n,p);
  }
  
  pwh=Calloc(n*n,float);
  pwhh=Calloc(n*p,float);
  
  if(*verbose==1) printf("Whitening\n");
  transpose_mat(ans,&n,&p,pwhh);

   mmult(ans,n,p,pwhh,p,n,pwh);

   Free(pwhh);
   for(i=0;i<n;i++){
     for(j=0;j<n;j++){pwh[n*i+j]=pwh[n*i+j]/p;}}

   uu=Calloc(n*n,float);
   dd=Calloc(n,float);
   vv=Calloc(n*n,float);
   
   svd(pwh,&n,&n,uu,dd,vv);

   pwhh=Calloc(n*n,float);
   for(i=0;i<n;i++){
     pwhh[n*i+i]=1/sqrt(dd[i]);
   }
   
   tmpm=Calloc(n*n,float);
   transpose_mat(uu,&n,&n,tmpm);
   mmult(pwhh,n,n,tmpm,n,n,pwh);
   
   Free(tmpm);Free(pwhh);
   
/*     Have ans as X, preprocessed data, (n*p) */
/*     Have pwh as K, prewhitening matrix, (n*n) */
   
   mat1=Calloc(e*n,float);/* kone*/
   pwhh=Calloc(e*p,float);/* xone*/
   
   for(i=0;i<e;i++){
     for(j=0;j<n;j++){mat1[i*n+j]=pwh[i*n+j];}} 
   mmult(mat1,e,n,ans,n,p,pwhh);
/*     have mat1 as K1 */
   
/*     Have pwhh as X1, (e*p) */
   ww=Calloc(e*e,float);/*  W */
   tmpm=Calloc(e*e,float);
   
   for(i=0;i<e;i++){for(j=0;j<e;j++){ww[i*e+j]=w_matrix[i*e+j];}}
   
   onefy(ww,e,tmpm); /*   Have tmpm as W1 */
   wwpl=Calloc(e*e,float);
 
   /*--------------------------------------------------------------*/
   
   if(*defflag==0){
     if(*funflag==1){

       if(*verbose==1) printf("Symmetric FastICA using tanh approx. to neg-entropy function\n"); 

       i=1;
/*         bff is initially taking W1 and X1 as args */
       bff(tmpm,e,pwhh,e,p,*alpha,wwpl,&tol);
       if(*verbose==1) printf("Iteration %d tol=%f\n",i,tol); 
       i=2;
       
       while((tol>(*lim)) && (i<(*maxit))){
	 bff(wwpl,e,pwhh,e,p,*alpha,wwpl,&tol);
	 if(*verbose==1) printf("Iteration %d tol=%f\n",i,tol); 
	 i+=1;
       } 
     }

     if(*funflag==2){
       if(*verbose==1) printf("Symmetric FastICA using exp approx. to neg-entropy function\n");

       i=1;
/*        bff is initially taking W1 and X1 as args */
       bff2(tmpm,e,pwhh,e,p,*alpha,wwpl,&tol);
       if(*verbose==1) printf("Iteration %d tol=%f\n",i,tol); 
       
       i=2;
       while((tol>(*lim)) && (i<(*maxit))){
	 bff2(wwpl,e,pwhh,e,p,*alpha,wwpl,&tol);
	 if(*verbose==1) printf("Iteration %d tol=%f\n",i,tol); 
	 i+=1;
       } 
     }
   }
   
   if(*defflag==1){
     Free(dd);   dd=Calloc(e,float);
     Free(uu);   uu=Calloc(e,float);
     
     if(*funflag==1){
       if(*verbose==1) printf("Deflation FastICA using tanh approx. to neg-entropy function\n"); 
       
       for(i=0;i<e;i++){k=0;
       gramsch(ww,e,e,i+1);
       rowstd(ww,e,e,i+1);
       tol=1;
       
       while((tol>(*lim)) && (k<(*maxit))){
	 for(j=0;j<e;j++){dd[j]=ww[i*e+j];}  
	 anfu(dd,e,pwhh,e,p,*alpha,uu);
	 for(j=0;j<e;j++){ww[i*e+j]=uu[j];}  
	 gramsch(ww,e,e,i+1);
	 rowstd(ww,e,e,i+1); 
	 tol=0;
	 for(j=0;j<e;j++){tol+=((dd[j])*(ww[i*e+j]));}
	 tol=(fabs(fabs(tol)-1));
	 k+=1;
       }            
       
       if(*verbose==1) printf("Component %d needed %d iterations tol=%f\n",i+1,k,tol); 
       
       }}
     if(*funflag==2){
     
       if(*verbose==1) printf("Deflation FastICA using tanh approx. to neg-entropy function\n"); 
       
       for(i=0;i<e;i++){k=0;
       gramsch(ww,e,e,i+1);
       rowstd(ww,e,e,i+1);
       tol=1;
       
       while((tol>(*lim)) && (k<(*maxit))){
	 for(j=0;j<e;j++){dd[j]=ww[i*e+j];}  
	 anfu2(dd,e,pwhh,e,p,*alpha,uu);
	 for(j=0;j<e;j++){ww[i*e+j]=uu[j];}  
	 gramsch(ww,e,e,i+1);
	 rowstd(ww,e,e,i+1); 
	 tol=0;
	 for(j=0;j<e;j++){tol+=((dd[j])*(ww[i*e+j]));}
	 tol=(fabs(fabs(tol)-1));
	 k+=1;
       }            
       
       if(*verbose==1) printf("Component %d needed %d iterations tol=%f\n",i+1,k,tol); 
       
       }}
     for(i=0;i<e;i++){for(j=0;j<e;j++){wwpl[i*e+j]=ww[i*e+j];}}   
   }
   
   
   
   mat2=Calloc(e*n,float);
   mat3=Calloc(e*p,float);
   mat4=Calloc(n*e,float);
   
   mmult(wwpl,e,e,mat1,e,n,mat2); /*  mat2 is UW (checked) */
   mmult(mat2,e,n,ans,n,p,mat3);  /*  mat3 is X2 (checked) */
   transpose_mat(mat2,&e,&n,mat4); /*  mat4 is t(UW) (checked) */
   
   Free(mat1); mat1=Calloc(e*e,float);
   mmult(mat2,e,n,mat4,n,e,mat1); /*  mat1 to be inverted */
   
   Free(uu);Free(dd);Free(vv);
   uu=Calloc(e*e,float);
   dd=Calloc(e,float);
   vv=Calloc(e*e,float);
   Free(mat2); 
   svd(mat1,&e,&e,uu,dd,vv); mat2=Calloc(e*e,float);
   for(i=0;i<e;i++){mat2[e*i+i]=1/(dd[i]);}
   
   mat5=Calloc(e*e,float);
   mat6=Calloc(e*e,float);
   transpose_mat(vv,&e,&e,mat6);
   mmult(mat6,e,e,mat2,e,e,mat5);
   transpose_mat(uu,&e,&e,vv);
   mmult(mat5,e,e,vv,e,e,uu);
   Free(mat2); mat2=Calloc(n*e,float);
   mmult(mat4,n,e,uu,e,e,mat2); /*  mat2 is A (checked)  */
 
   for(i=0;i<n;i++){for(j=0;j<p;j++){ansx[i*p+j]=ans[i*p+j];}}
   for(i=0;i<n;i++){for(j=0;j<n;j++){ansk[i*n+j]=pwh[i*n+j];}}
   for(i=0;i<e;i++){for(j=0;j<e;j++){answ[i*e+j]=wwpl[i*e+j];}}
   for(i=0;i<n;i++){for(j=0;j<e;j++){ansa[i*e+j]=mat2[i*e+j];}}
   for(i=0;i<e;i++){for(j=0;j<p;j++){ansx2[i*p+j]=mat3[i*p+j];}}
 
 
 
   
   Free(mat2);Free(mat3);Free(mat4);Free(mat5);Free(mat6);
   Free(uu);
   Free(dd);
   Free(vv);
   Free(ans);
   Free(pwh);
   Free(tmpm);
   Free(pwhh); 
   Free(ww);
   Free(wwpl);
   Free(mat1);
}

















