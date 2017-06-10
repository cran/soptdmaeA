soptdmaeA.default<-function(trt.N, blk.N, theta, nrep, strt, sary, des0, dtype, Optcrit = " ",...)
{
  trt.N=as.numeric(trt.N)
  blk.N=as.numeric(blk.N)
  theta=as.numeric(theta)
  nrep=as.numeric(nrep)
  des0<-as.matrix(des0)
  if(dim(des0)[1]!=2){des0=t(des0)}
  if(sary<strt) stop(paste("additional number of arrays should be greater or equal to the additional number of treatments"))
  if(!max(des0)==trt.N||!(length(des0)/2)==blk.N){stop("Number of treatments and/or arrays in the initial design is different from the value you set")}
  for(j in 1:blk.N) {if(des0[1,j]==des0[2,j]) stop("Invalid initial design setup. Treatments within each array should be distinct")}
  trtin<-contrasts(as.factor(t(des0)),contrasts=FALSE)[as.factor(des0),]
  R.trt<-t(trtin)%*%trtin
  if (rankMatrix(R.trt)[1]<trt.N)  stop("Initial design do not consists of all initial treatments, please insert valid connected initial design")
  cmato=cmatbrcd.mae(trt.N,blk.N, 0,des0, dtype)
  egv<-sort(eigen(cmato)$values)
  if(egv[2]<0.000001) stop("Initial design is not connected  ",if(dtype=="rcd") {"row-column"} else{"block"} ," design, please insert valid connected initial design")
  #"===================================================================================================="
  if(is.na(theta)|theta<0|theta>1){
    tkmessageBox(title="Error",message=paste("Please insert correct value of theta, it should be between 0 and 1 inclusive of 0 and 1, click OK to reset.",sep=""),icon = "error"); 
    stop("Wrong value of 'theta', it should be between 0 and 1, inclusive of 0 and 1")
  }#end of if
  if(is.na(trt.N)|is.na(blk.N)|trt.N!=round(trt.N)|blk.N!=round(blk.N)) {
    tkmessageBox(title="Error",message=paste("Please insert correct format of the number of treatments and arrays. The number of treatments and arrays should be an integer, click OK to reset the values.",sep=""),icon = "error"); 
    stop("Wrong format of 'trt.N' and/or 'blk.N', both should be an integer")
  }#end of if
  if(trt.N<3|blk.N<3){ 
    tkmessageBox(title="Error",message=paste("The number of arrays and treatments should be greater than or equal to 3, click Ok to reset.",sep=""),icon = "error"); 
    stop("Very small value of number of treatments and/or arrays, minimum value of the two is 3")
  }#end of if
  if(trt.N-blk.N>1){ 
    tkmessageBox(title="Error",message=paste("The number of arrays should be greater than or equal to the number of treatments minus one, click Ok to reset.",sep=""),icon = "error"); 
    stop("The number of treatments are larger than the number of arrays minus one (trt.N>blk.N-1)")
  }#end of if
  if(trt.N>10|blk.N>10){ 
    tkmessageBox(title="Information",message=paste("This might take some minutes, please be patient...",sep=""))
  }#end of if
  if(is.na(nrep)|nrep<1|nrep!=round(nrep)){
    tkmessageBox(title="Error",message=paste("The number of replications should be a positive integer greater than or equal to one, click OK to reset.",sep=""),icon = "error"); 
    stop("Wrong value of 'nrep', it should be greater than or equal to one (only positive integer values)")
  }#end of if
  #"===================================================================================================="
  if(!is.element(Optcrit,c("A","MV","D","E"))){stop("The optimality criterion 'Optcrit' is not correctly specified")}
  cat("\nSerching for sequential ",Optcrit, "-optimal or near-optimal ",if(dtype=="rcd") {"row-column"} else{"block"} ," design ...\n",sep="")
  if(Optcrit=="A") {
    seqoptbrcd_mae<-seqAoptbrcd.maeA(trt.N,blk.N,theta,nrep,strt,sary,des0,dtype)} else if(
      Optcrit=="MV") {
      seqoptbrcd_mae<-seqMVoptbrcd.maeA(trt.N,blk.N,theta,nrep,strt,sary,des0,dtype)} else if(
        Optcrit=="D") {
        seqoptbrcd_mae<-seqDoptbrcd.maeA(trt.N,blk.N,theta,nrep,strt,sary,des0,dtype)} else if(
          Optcrit=="E") {
          seqoptbrcd_mae<-seqEoptbrcd.maeA(trt.N,blk.N,theta,nrep,strt,sary,des0,dtype)} else{
            stop("The optimality criterion is not specified")}
  seqoptbrcd_mae$call<-match.call()
  seqoptbrcd_mae$Optcrit<-Optcrit
  seqoptbrcd_mae$Cmat<-cmatbrcd.mae(seqoptbrcd_mae$v,seqoptbrcd_mae$b,seqoptbrcd_mae$theta,seqoptbrcd_mae$soptdesF,dtype)
  trtin <- contrasts(as.factor(seqoptbrcd_mae$soptdesF), contrasts = FALSE)[as.factor(seqoptbrcd_mae$soptdesF), ]
  vec1 <- rep(1, seqoptbrcd_mae$b * 2)
  vec_trtr <- t(trtin) %*% vec1
  seqoptbrcd_mae$equireplicate<-all(vec_trtr==vec_trtr[1])
  seqoptbrcd_mae$vtrtrep<-t(vec_trtr)
  trtin0 <- contrasts(as.factor(seqoptbrcd_mae$optdes0), contrasts = FALSE)[as.factor(seqoptbrcd_mae$optdes0), ]
  vec10 <- rep(1, (seqoptbrcd_mae$b-seqoptbrcd_mae$sary) * 2)
  vec_trtr0 <- t(trtin0) %*% vec10
  seqoptbrcd_mae$equireplicate0<-all(vec_trtr0==vec_trtr0[1])
  seqoptbrcd_mae$vtrtrep0<-t(vec_trtr0)
  #"======================================================================================"
  titleoptbd<-list(c("      --------------------------------------- ",paste("Title: Sequential ",Optcrit,"-optimal or near-optimal ",if(dtype=="rcd") {"row-column"} else{"block"} ," design          Date:", format(Sys.time(), "%a %b %d %Y %H:%M:%S"),sep=""),
                     "      --------------------------------------- "))
  write.table(titleoptbd, file = paste("seq",Optcrit,"opt",dtype,"_summary.csv",sep = ""),append=T ,sep = ",", row.names=FALSE, col.names=FALSE)
  parcomb<-list(c("     Parametric combination:", "Number of treatments:", "Number of arrays:", 
                  "Theta value:", "Number of replications:", "Number of added treatments: ","Number of added arrays:","Optimality criterion used:"," ",paste("Initial ",Optcrit,"-optimal or near-optimal ",if(dtype=="rcd") {"row-column"} else{"block"} ," design:",sep="")),
                 c(" ",seqoptbrcd_mae$v,seqoptbrcd_mae$b,seqoptbrcd_mae$theta,seqoptbrcd_mae$nrep,seqoptbrcd_mae$strt,seqoptbrcd_mae$sary, paste(Optcrit,"-optimality criterion",sep="")," "," "))
  write.table(parcomb, file =  paste("seq",Optcrit,"opt",dtype,"_summary.csv",sep = ""),append=T ,sep = ",", row.names=FALSE, col.names=FALSE)
  
  optde0<-list("",cbind(if(dtype=="rcd"){c(" ", "Dye 1:", "Dye 2:")},rbind(paste0("Ary",1:(seqoptbrcd_mae$b-seqoptbrcd_mae$sary)),seqoptbrcd_mae$optdes0)))
  write.table(optde0, file =  paste("seq",Optcrit,"opt",dtype,"_summary.csv",sep = ""),append=T ,sep = ",", row.names=FALSE, col.names=FALSE)
  write.table(list(c("",paste(Optcrit,"-Score value of initial design:  " ,sep=""), "Equireplicate  initial design:",""),c("",seqoptbrcd_mae$optcrtsv0,seqoptbrcd_mae$equireplicate0,"")), file =   paste("seq",Optcrit,"opt",dtype,"_summary.csv",sep = ""),append=T ,sep = ",", row.names=FALSE, col.names=FALSE)
  write.table(list(c("Treatment:", "Treatment replication of initial design:"),rbind(1:(seqoptbrcd_mae$v-seqoptbrcd_mae$strt),seqoptbrcd_mae$vtrtrep0)), file =  paste("seq",Optcrit,"opt",dtype,"_summary.csv",sep = ""),append=T ,sep = ",", row.names=FALSE, col.names=FALSE)
  
  reseqd<-list(c("",paste("Resultant sequential ",Optcrit,"-optimal or near-optimal ",if(dtype=="rcd") {"row-column"} else{"block"} ," design:",sep="")))
  write.table(reseqd, file = paste("seq",Optcrit,"opt",dtype,"_summary.csv",sep = ""),append=T ,sep = ",", row.names=FALSE, col.names=FALSE)
  
  optde<-list("",cbind(if(dtype=="rcd"){c(" ", "Dye 1:", "Dye 2:")}, rbind(paste0("Ary",1:seqoptbrcd_mae$b),seqoptbrcd_mae$soptdesF)))
  write.table(optde, file =  paste("seq",Optcrit,"opt",dtype,"_summary.csv",sep = ""),append=T ,sep = ",", row.names=FALSE, col.names=FALSE)
  write.table(list(c("",paste(Optcrit,"-score value of seq. design:",sep=""), "Equireplicate seq. design:",""),c("",seqoptbrcd_mae$soptcrtsv,seqoptbrcd_mae$equireplicate,"")), file =   paste("seq",Optcrit,"opt",dtype,"_summary.csv",sep = ""),append=T ,sep = ",", row.names=FALSE, col.names=FALSE)
  write.table(list(c("Treatment:", "Treatment replication  of seq. design:"),rbind(1:seqoptbrcd_mae$v,seqoptbrcd_mae$vtrtrep)), file =  paste("seq",Optcrit,"opt",dtype,"_summary.csv",sep = ""),append=T ,sep = ",", row.names=FALSE, col.names=FALSE)
  
  seqoptbrcd_mae$file_loc<-paste(as.character(getwd()),  paste("seq",Optcrit,"opt",dtype,"_summary.csv",sep = ""),sep = "/")
  seqoptbrcd_mae$file_loc2<-paste("Summary of obtained sequential ",Optcrit,"-optimal or near-optimal ",if(dtype=="rcd") {"row-column"} else{"block"} ," design is also saved at:",sep="")
  #"======================================================================================"
  class(seqoptbrcd_mae)<-"soptdmaeA"
  seqoptbrcd_mae
}