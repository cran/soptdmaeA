#Function for construction of sequential E-optimal block/row-column design using array exchange algorithm) 
seqEoptbrcd.maeA<-function(trt.N,blk.N,theta,nrep,strt,sary,des0,dtype) {
  cmat<-cmatbrcd.mae(trt.N,blk.N,theta,des0,dtype)
  Eegv<-sort(eigen(cmat)$values)
  Escore0<- Eegv[2]
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
  del.1<-matrix(0,na,3)
  desbest.1<-matrix(0,nrep*2,blk.N)
  eoptbest.1<-matrix(0,nrep,2)
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
    Eegv<-sort(eigen(cmat)$values)
    eopt<- Eegv[2]
    ecold=eopt
    descold=t(des)
    cdel=100
    for(i in (blk.N-sary+1):blk.N){
      for(j in 1:na){
        temp=descold[i,]
        if(all(descold[i,]==arrays[j,]))  {eopt=ecold; del.1[j,]<-c(j,(ecold-eopt),eopt); next}
        descold[i,]=arrays[j,]
        trtin<-contrasts(as.factor(t(descold)),contrasts=FALSE)[as.factor(t(descold)),]
        R.trt<-t(trtin)%*%trtin
        if (rankMatrix(R.trt)[1]<trt.N)  {eopt=0; del.1[j,]<-c(j,(ecold-eopt),eopt); next}
        cmato=cmatbrcd.mae(trt.N,blk.N, 0,t(descold), dtype)
        egv<-sort(eigen(cmato)$values)
        if(egv[2]<0.000001) {eopt=0; del.1[j,]<-c(j,(ecold-eopt),eopt); next}
        cmat=cmatbrcd.mae(trt.N,blk.N,theta,t(descold), dtype)
        eopt=sum(diag(ginv(cmat)))
        del.n<-del.1[j,]<-c(j,(ecold-eopt),eopt)
        descold[i,]=temp
      }
      del.1<-del.1[rev(order(del.1[,3])),]
      delbest=t(del.1[1,])
      descold[i,]=arrays[delbest[1],]
      ecold=delbest[3]
      cdel=delbest[2]
    }
    if (irep==1) {desbest.1=t(descold)} else {desbest.1=rbind(desbest.1,t(descold))}
    eoptbest.1[irep,]=c(irep,ecold)
  }
  if(nrep==1){nb=eoptbest.1[1]; Escore=eoptbest.1[2]}else{
    best=eoptbest.1[rev(order(eoptbest.1[,2])),]
    nb=best[1,1]
    Escore<-best[1,2]
  }
  Eoptde<- desbest.1[c((nb-1)*2+1,nb*2),]
  tkmessageBox(title="Search completed",message=paste("Search completed",sep=""))
  cnames=paste0("Ary",1:blk.N)
  dimnames(Eoptde)=list(if(dtype=="rcd"){c( "Dye 1:", "Dye 2:")} else {NULL},cnames)
  dimnames(des0)=list(if(dtype=="rcd"){c("Dye 1:", "Dye 2:")} else {NULL},paste0("Ary",1:(blk.N-sary)))
  Eopt_sum2<-list("v"=trt.N,"b"=blk.N,theta=theta,nrep=nrep, strt=strt, sary=sary, dtype=dtype, "optdes0"=des0,"optcrtsv0" =Escore0,"soptdesF"=Eoptde,"soptcrtsv" =Escore)
  return(Eopt_sum2)
}#End of SeqEoptbrcd.maeA function: construction of sequential E-optimal block design using array exchange algorithm
