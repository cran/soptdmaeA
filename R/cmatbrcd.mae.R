#Computation of the information matrix (C-matrix) for a given block/row-column design (des) 
cmatbrcd.mae<-function(trt.N,blk.N,theta,des, dtype){
  k=2
  if(!is.element(dtype,c("blkd","rcd"))){stop("The design type is not correctly specified")}
    trtin<-contrasts(as.factor(des),contrasts=FALSE)[as.factor(des),]
  blk.1<-rep(1:blk.N,each=2)
  blkin<-contrasts(as.factor(blk.1),contrasts=FALSE)[as.factor(blk.1),]
  vec.1<-rep(1,blk.N*2)
  R.trt<-t(trtin)%*%trtin
  N.tb<-t(trtin)%*%blkin
  r.trt<-t(trtin)%*%vec.1
    cmat<-R.trt-(1/k)*(N.tb%*%t(N.tb))+theta*((1/k)*(N.tb%*%t(N.tb))-(1/(blk.N*k))*(r.trt%*%t(r.trt)))
    if(dtype=="rcd"){
  forrow.1=contrasts(as.factor(c(des[1,],seq(1:trt.N))),contrasts=FALSE)[as.factor(c(des[1,],seq(1:trt.N))),]
  forrow.2=contrasts(as.factor(c(des[2,],seq(1:trt.N))),contrasts=FALSE)[as.factor(c(des[2,],seq(1:trt.N))),]
  mmmat.r=t(rbind(t(as.matrix(colSums(forrow.1)-1)), t(as.matrix(colSums(forrow.2)-1))))%*% 
    rbind(t(as.matrix(colSums(forrow.1)-1)), t(as.matrix(colSums(forrow.2)-1)))
  cmat<-cmat-(1/blk.N)*(mmmat.r)+(1/(blk.N*k))*(r.trt%*%t(r.trt))
  }
  cmat
}#End of computation of C-matrix, cmatbrcd.mae
