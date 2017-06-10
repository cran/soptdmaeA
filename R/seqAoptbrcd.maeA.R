#Function for construction of sequential A-optimal block design using array exchange algorithm 
seqAoptbrcd.maeA<-function(trt.N,blk.N,theta,nrep,strt,sary,des0,dtype) {
  cmat<-cmatbrcd.mae(trt.N,blk.N,theta,des0,dtype)
  Ascore0=sum(diag(ginv(cmat)))
  #New parametric combinations 
  blk.N=blk.N+sary
  trt.N=trt.N+strt
  #Candidate arrays
  arrays={}
  if(strt==sary){
    for(i in 1:strt){
      arrays0<-cbind(1:(trt.N-strt+i-1),trt.N-strt+i)
      arrays<-as.matrix(rbind(arrays,arrays0))
    }
  } else {arrays<-t(combn(trt.N,2))}
  if(dtype=="rcd") {arrays<-rbind(arrays,t(rbind(arrays[,2],arrays[,1])))}
  na=dim(arrays)[1]
  #House keeping
  del.1<-matrix(1000,na,3)
  desbest.1<-matrix(0,nrep*2,blk.N)
  aoptbest.1<-matrix(0,nrep,2)
  #Start iteration for search of optimal sequential design given initial design
  for(irep in 1:nrep){
    #select connected initial sequential design 
      Con.egv.check=0.00000001
      while(Con.egv.check<0.000001){
        All.trt.check=trt.N-1
        while(All.trt.check <trt.N){
          sqarray=matrix(0,strt,2);
          if(strt>0) {
            for(i in 1:strt){
              sqary=cbind(1:(trt.N-strt+i-1),trt.N-strt+i)
              if(dtype=="rcd") {sqary=rbind(sqary,t(rbind(sqary[,2],sqary[,1])))}
              sqarray0<-sqary[sample(1:(length(sqary)/2),1),]
              sqarray[i,]=sqarray0
            }
          }
          if((sary-strt)>0) {
            adary0<-if((sary-strt)>na) {c(1:na, sample(1:na,(sary-strt-na),replace=TRUE))} else {
              sample(1:na,(sary-strt),replace=FALSE)}
            sqarray<-rbind(sqarray,arrays[adary0,])} 
          des=t(rbind(t(des0),sqarray))
          trtin<-contrasts(as.factor(des),contrasts=FALSE)[as.factor(des),]
          R.trt<-t(trtin)%*%trtin
          All.trt.check<-rankMatrix(R.trt)
        }
        cmato=cmatbrcd.mae(trt.N,blk.N, 0,des, dtype)
        egv<-sort(eigen(cmato)$values)
        Con.egv.check<-egv[2]
      }
    cmat<-cmatbrcd.mae(trt.N, blk.N, theta, des, dtype)
    aopt=sum(diag(ginv(cmat)))
    acold=aopt
    descold=t(des)
    cdel=100
    for(i in (blk.N-sary+1):blk.N){
      for(j in 1:na){
        temp=descold[i,]
        if(all(descold[i,]==arrays[j,]))  {aopt=acold; del.1[j,]<-c(j,(acold-aopt),aopt); next}
        descold[i,]=arrays[j,]
        trtin<-contrasts(as.factor(t(descold)),contrasts=FALSE)[as.factor(t(descold)),]
        R.trt<-t(trtin)%*%trtin
        if (rankMatrix(R.trt)[1]<trt.N)  {aopt=1000; del.1[j,]<-c(j,(acold-aopt),aopt); next}
        cmato=cmatbrcd.mae(trt.N,blk.N, 0,t(descold), dtype)
        egv<-sort(eigen(cmato)$values)
        if(egv[2]<0.000001) {aopt=1000; del.1[j,]<-c(j,(acold-aopt),aopt); next}
        cmat=cmatbrcd.mae(trt.N,blk.N,theta,t(descold), dtype)
        aopt=sum(diag(ginv(cmat)))
        del.n<-del.1[j,]<-c(j,(acold-aopt),aopt)
        descold[i,]=temp
      }
      del.1<-del.1[order(del.1[,3]),]
      delbest=t(del.1[1,])
      descold[i,]=arrays[delbest[1],]
      acold=delbest[3]
      cdel=delbest[2]
    }
    if (irep==1) {desbest.1=t(descold)} else {desbest.1=rbind(desbest.1,t(descold))}
    aoptbest.1[irep,]=c(irep,acold)
  }
  if(nrep==1){nb=aoptbest.1[1]; Ascore=aoptbest.1[2]}else{
    best=aoptbest.1[order(aoptbest.1[,2]),]
    nb=best[1,1]
    Ascore<-best[1,2]
  }
  Aoptde<- desbest.1[c((nb-1)*2+1,nb*2),]
  tkmessageBox(title="Search completed",message=paste("Search completed",sep=""))
  cnames=paste0("Ary",1:blk.N)
  
  dimnames(Aoptde)=list(if(dtype=="rcd"){c("Dye 1:", "Dye 2:")} else {NULL},cnames)
  dimnames(des0)=list(if(dtype=="rcd"){c("Dye 1:", "Dye 2:")} else {NULL},paste0("Ary",1:(blk.N-sary)))
  Aopt_sum2<-list("v"=trt.N,"b"=blk.N,theta=theta,nrep=nrep, strt=strt, sary=sary, dtype=dtype, "optdes0"=des0,"optcrtsv0" =Ascore0,"soptdesF"=Aoptde,"soptcrtsv" =Ascore)
  return(Aopt_sum2)
}#End of SeqAoptbrcd.maeA function: construction of sequential A-optimal block design using array exchange algorithm
