generate_null_data<-function(X)
{
  n<-dim(X)[1]
  p<-dim(X)[2]
  X.pca<-prcomp(X, scale. = TRUE)
  
  X_null<-matrix(rnorm(n*p), ncol = p, nrow = n)
  
  #X_pca_null<-matrix(rnorm(n*p), ncol = n, nrow = p)
  #X_null<-t(scale(X.pca$rotation)%*%X_pca_null)
  
  rownames(X_null)<-rownames(X)
  colnames(X_null)<-colnames(X)
  
  return(X_null)
}