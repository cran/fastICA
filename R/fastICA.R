
FastICA<-function(X,n.comp,alg.typ="deflation",fun="logcosh",alpha=1.0,method="R",col.norm=F,maxit=200,tol=0.0001,verbose=F){

  if(is.matrix(X)==F) return("X is not a matrix")
  if(alpha<1 || alpha>2) return("alpha must be in range [1,2]")
  
  n<-nrow(X)
  p<-ncol(X)
  
  if(n.comp>min(n,p)){
    print("n.comp is too large")
    print(paste("n.comp set to",min(n,p)))
    n.comp<-min(n,p)}

  if(method=="R"){

    if(verbose==T){print("Centering")}
    X<-t(apply(X,1,FUN=function(x){x-mean(x)}))
    if(col.norm==T){X<-apply(X,2,FUN=function(x){(x-mean(x))/sqrt(var(x))})
                  }
    
    
  #whiten the data matrix
    if(verbose==T){print("Whitening")}
    V<-X%*%t(X)/p
    s<-svd(V)
    D<-diag(c(1/sqrt(s$d)))
    K<-D%*%t(s$u)
    K<-matrix(K[1:n.comp,],n.comp,n)   
    X1<-K%*%X
    
    
    if(alg.typ=="deflation"){
      a<-ica.R.def(X1,n.comp,tol=tol,fun=fun,alpha=alpha,maxit=maxit,verbose=verbose)
    }
    else if(alg.typ=="parallel"){
      a<-ica.R.par(X1,n.comp,tol=tol,fun=fun,alpha=alpha,maxit=maxit,verbose=verbose)
    }
    else return("Invalid alg.type")
    
    w<-a%*%K #extract unmixing matrix
    S<-w%*%X #unmix data matrix
    A<-t(w)%*%solve(w%*%t(w)) #calculate mixing matrix
    
    return(list(X=X,K=K,W=a,A=A,S=S))
    
    
  }
  
  else if(method=="C"){
    
    w.matrix<-matrix(rnorm(n.comp^2),n.comp,n.comp)
    col.n<-0;if(col.norm==T){col.n<-1}
    funflag<-1;if(fun=="exp"){funflag<-2}
    defflag<-1;if(alg.typ=="parallel"){defflag<-0}
    verbose.flag<-0;if(verbose==T){verbose.flag<-1}
    
    a<-.C("icainc",
          as.single(t(X)),
          as.single(t(w.matrix)),
          as.integer(n),
          as.integer(p),
          as.integer(n.comp),
          as.single(alpha),
          as.integer(1),
          as.integer(col.n),
          as.integer(funflag),
          as.integer(maxit),
          as.single(tol),
          as.integer(defflag),
          as.integer(verbose.flag),
          X=single(n*p),
          K=single(n*n),
          W=single(n.comp*n.comp),
          A=single(n*n.comp),
          S=single(n.comp*p))


    X1<-matrix(a$X,n,p,byrow=T)
    K<-matrix(a$K,n,n,byrow=T)
    W<-matrix(a$W,n.comp,n.comp,byrow=T)
    A<-matrix(a$A,n,n.comp,byrow=T)
    S<-matrix(a$S,n.comp,p,byrow=T)
    
    return(list(X=X1,K=K,W=W,A=A,S=S))
  }
  else {
    return("Invalid method")
  }

  
  
}

  



ica.R.def<-function(X,n.comp,tol,fun,alpha,maxit,verbose){
  
  n<-nrow(X)
  p<-ncol(X)
  
  #set up unmixing matrix W
  W<-matrix(0,n.comp,n.comp)
  
  if(verbose==T){print("Performing ICA")}
  for(i in 1:n.comp){
    if(verbose==T){print(paste("Component",i))}

    #create random unmixing vector w 
    w<-matrix(rnorm(n.comp),n.comp,1)
    
    #orthonormalise using Gram-Schmidt
    if(i>1){
      t<-w
      t[1:length(t)]<-0
      for(u in 1:(i-1)){
        k<-sum(w*W[u,])
        t<-t+k*W[u,]}
      w<-w-t
    }
    w<-w/sqrt(sum(w^2))

    lim<-rep(1000,maxit)
    it<-1
   
    if(fun=="logcosh"){ #ICA using logcosh approx. to neg. entropy
      while(lim[it]>tol && it<maxit){

        #calculate new unmixing vector
        wx<-t(w)%*%X
        gwx<-tanh(alpha*wx)
        gwx<-matrix(gwx,n.comp,p,byrow=T)
        xgwx<-X*gwx
        v1<-apply(xgwx,1,FUN=mean)
        g.wx<-alpha*(1-(tanh(alpha*wx))^2)
        v2<-mean(g.wx)*w
        w1<-v1-v2
        w1<-matrix(w1,n.comp,1)
      
        it<-it+1

        #orthonormalise using Gram-Schmidt
        if(i>1){
          t<-w1
          t[1:length(t)]<-0
          for(u in 1:(i-1)){
            k<-sum(w1*W[u,])
            t<-t+k*W[u,]}
          w1<-w1-t
        }
        w1<-w1/sqrt(sum(w1^2))
        
        lim[it]<-Mod(Mod(sum((w1*w)))-1)
        if(verbose==T){print(paste("Iteration",it-1,"tol =",lim[it]))}
        w<-matrix(w1,n.comp,1)
      }
    }
    
    
    if(fun=="exp"){ #ICA using exponential approx. to neg. entropy
      while(lim[it]>tol && it<maxit){
        
        #calculate new unmixing vector
        wx<-t(w)%*%X
        gwx<-wx*exp(-(wx^2)/2)
        gwx<-matrix(gwx,n.comp,p,byrow=T)
        xgwx<-X*gwx
        v1<-apply(xgwx,1,FUN=mean)
        g.wx<-(1-wx^2)*exp(-(wx^2)/2)
        v2<-mean(g.wx)*w
        w1<-v1-v2
        w1<-matrix(w1,n.comp,1)
        
        it<-it+1
        
        #orthonormalise using Gram-Schmidt
        if(i>1){
          t<-w1
          t[1:length(t)]<-0
          for(u in 1:(i-1)){
          k<-sum(w1*W[u,])
          t<-t+k*W[u,]}
        w1<-w1-t
      }
      w1<-w1/sqrt(sum(w1^2))
        
      lim[it]<-Mod(Mod(sum((w1*w)))-1)
      if(verbose==T){print(paste("Iteration",it-1,"tol =",lim[it]))}
      w<-matrix(w1,n.comp,1)
      }
    }
    W[i,]<-w
  }

  
  return(W)}


ica.R.par<-function(X,n.comp,tol,fun,alpha,maxit,verbose){
  
  n<-nrow(X)
  p<-ncol(X)
  
  #set up unmixing matrix W
  W<-matrix(rnorm(n.comp*n.comp),n.comp,n.comp)
  sW<-svd(W)
  W<-sW$u%*%diag(1/sW$d)%*%t(sW$u)%*%W
  W1<-W
  
  lim<-rep(1000,maxit)
  it<-1
  
  if(verbose==T){print("Performing ICA")}
  if(fun=="logcosh"){ #ICA using logcosh approx. to neg. entropy
    while(lim[it]>tol && it<maxit){
      
      #calculate new estimate of unmixing matrix
      wx<-W%*%X
      gwx<-tanh(alpha*wx)
      v1<-gwx%*%t(X)/p
      g.wx<-alpha*(1-(gwx)^2)
      v2<-diag(apply(g.wx,1,FUN=mean))%*%W
      W1<-v1-v2
      
      #symmetric orthogonalisation of unmixing matrix
      sW1<-svd(W1) 
      W1<-sW1$u%*%diag(1/sW1$d)%*%t(sW1$u)%*%W1
      
      lim[it+1]<-max(Mod(Mod(diag(W1%*%t(W)))-1))
      W<-W1
      if(verbose==T){print(paste("Iteration",it,"tol =",lim[it+1]))}
      
      it<-it+1 
      
    }
  }
  
  if(fun=="exp"){ #ICA using exponential approx. to neg. entropy
    while(lim[it]>tol && it<maxit){
      
      #calculate new estimate of unmixing matrix
      wx<-W%*%X
      gwx<-wx*exp(-(wx^2)/2)
      v1<-gwx%*%t(X)/p
      g.wx<-(1-wx^2)*exp(-(wx^2)/2)
      v2<-diag(apply(g.wx,1,FUN=mean))%*%W
      W1<-v1-v2
      
      #symmetric orthogonalisation of unmixing matrix
      sW1<-svd(W1)
      W1<-sW1$u%*%diag(1/sW1$d)%*%t(sW1$u)%*%W1
      
      lim[it+1]<-max(Mod(Mod(diag(W1%*%t(W)))-1))
      W<-W1
      if(verbose==T){print(paste("Iteration",it,"tol =",lim[it+1]))}
      
      it<-it+1 
    }
  }
  
  return(W)
}
  










