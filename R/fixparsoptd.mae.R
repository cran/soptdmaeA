# Function to create  a TCL/TK window to set the parametric values and initial design
globalVariables("des0")
fixparsoptd.mae<-function(Optcrit){
  trt.I<-tclVar(paste("3"));#Number of treatments (default)
  blk.I<-tclVar(paste("3"));#Number of arrays (default)
  theta.I<-tclVar(paste("0.0"));#Theta value (default)
  rep.I<-tclVar(paste("10"));#Number of replications (default)
  strt.I<-tclVar(paste("0"));#Initial convergence value (default)
  sary.I<-tclVar(paste("1"));#Initial convergence value (default)
  cbValue.I<-tclVar("0")
  cbValue2.I<-tclVar("0")  
  cbValue3.I<-tclVar("0")  
  optcrtF<-function(Optcrit){
    nrep<-as.numeric(tclvalue(rep.I))
    trt.N<-as.numeric(tclvalue(trt.I))
    blk.N<-as.numeric(tclvalue(blk.I))
    theta<-as.numeric(tclvalue(theta.I))
    cbVal<-as.numeric(tclvalue(cbValue.I))
    strt<-as.numeric(tclvalue(strt.I))
    sary<-as.numeric(tclvalue(sary.I))
    #if(itr.cvrgval>blk.N) itr.cvrgval<-blk.N
    cbVal2<-as.numeric(tclvalue(cbValue2.I))
    cbVal3<-as.numeric(tclvalue(cbValue3.I))
    if(cbVal2==0) {dtype="blkd"} else if(cbVal2==1) {dtype="rcd"}
    if(Optcrit=="A") optrcdFS=soptdmaeA(trt.N, blk.N, theta, nrep, strt, sary, des0, dtype, Optcrit = "A") 
    if(Optcrit=="MV") optrcdFS=soptdmaeA(trt.N, blk.N, theta, nrep, strt, sary, des0, dtype, Optcrit = "MV")
    if(Optcrit=="D") optrcdFS=soptdmaeA(trt.N, blk.N, theta, nrep, strt, sary, des0, dtype, Optcrit = "D")
    if(Optcrit=="E") optrcdFS=soptdmaeA(trt.N, blk.N, theta, nrep, strt, sary, des0, dtype, Optcrit = "E")
    if(cbVal3=="1") optrcdFS=summary(optrcdFS) 
    print(optrcdFS)
    if(cbVal=="1") {graphsoptd.mae(trt.N, blk.N,theta,optrcdFS$soptdesF,Optcrit,strt,sary,dtype)}
  }
  fiPar<-tktoplevel()
  fontHeading<- tkfont.create(family="times",size=40,weight="bold")
  fontHeading3<-tkfont.create(family="times",size=10,weight="bold")
  fontHeading2<-tkfont.create(family="times",size=14,weight="bold")
  fontHeading4<-tkfont.create(family="times",size=12,weight="bold")
  tkwm.title(fiPar,"Set parameter values")
  fiParF<-tkframe(fiPar)
  fiParFup<- tkframe(fiParF,relief="groove",borderwidth=2)
  fiParFmid<-tkframe(fiParF,relief="sunken",borderwidth=2)
  fiParFlow<-tkframe(fiParF,relief="sunken",borderwidth=2)
  fiParFlow2<-tkframe(fiParF,relief="groove",borderwidth=2)
  trt.in<-trt.I
  blk.in<-blk.I
  theta.in<-theta.I
  rep.in<-rep.I
  strt.in<-strt.I
  sary.in<-sary.I
  cbValue.in<-cbValue.I
  cbValue2.in<-cbValue2.I
  cbValue3.in<-cbValue3.I
  trt.i1<-tkentry(fiParFup,width=11,textvariable=trt.in)
  blk.i1<-tkentry(fiParFup,width=11,textvariable=blk.in)
  theta.i1<-tkentry(fiParFup,width=11,textvariable=theta.in)
  rep.i1<-tkentry(fiParFup,width=11,textvariable=rep.in)
  strt.i1<-tkentry(fiParFup,width=11,textvariable=strt.in)
  sary.i1<-tkentry(fiParFup,width=11,textvariable=sary.in)
  cbValue.i1<-tkcheckbutton(fiPar)
  tkconfigure(cbValue.i1,variable=cbValue.in)
  cbValue.i3<-tkcheckbutton(fiPar)
  tkconfigure(cbValue.i3,variable=cbValue3.in)
  desoI.but<-tkbutton(fiParF,text="Insert",command=function() desoI(matrix(c(1:3,2:3,1),3,2)),width=9)
  TrtEx<-tkradiobutton(fiPar)
  ArEx<-tkradiobutton(fiPar)
  tkconfigure(TrtEx,variable=cbValue2.in,value="0")
  tkconfigure(ArEx,variable=cbValue2.in,value="1")
  tkgrid(tklabel(fiParF,text="       Fix Value of Parameters     ",font=fontHeading2))
  tkgrid(tklabel(fiParFup,text="Number of treatments:                  "),trt.i1)
  tkgrid(tklabel(fiParFup,text="Number of arrays:                           "),blk.i1)
  tkgrid(tklabel(fiParFup,text="Theta value:                                     "),theta.i1)
  tkgrid(tklabel(fiParFup,text="Number of added treatments:      "),strt.i1)
  tkgrid(tklabel(fiParFup,text="Number of added arrays:              "),sary.i1)
  tkgrid(tklabel(fiParFup,text="Number of replications:                "),rep.i1)
  tkgrid(tklabel(fiParFup,text="Insert initial design:                      "),desoI.but)
  tkgrid(tklabel(fiParFlow,text="Type of design:      ", font=fontHeading4))
  tkgrid(tklabel(fiParFlow,text="Block design                    "),TrtEx)
  tkgrid(tklabel(fiParFlow,text="Row-column design       "),ArEx)
  tkgrid(tklabel(fiParFmid,text="Show graphical layout        ", font=fontHeading4),cbValue.i1)
  tkgrid(tklabel(fiParFmid,text="Show Summary result        ", font=fontHeading4),cbValue.i3)
  desoI<-function(des0) {
    des0<-as.matrix(des0)
    dimnames(des0)=list(NULL,c("Dye1","Dye2"))
    des0<-fix(des0)
    if(dim(des0)[1]!=2){des0=t(des0)}
    return(des0)}#end of desoI
  exitFP<-function() {
    closeQ=tkmessageBox(title = "Exit set parameter values", message = "You are leaving set parameter values window",
                        icon = "info", type = "okcancel", default = "cancel")
    if(as.character(closeQ)=="ok") tkdestroy(fiPar)
  }#end of exitFP
  exit.but2a<-tkbutton(fiParF,text="Exit",command=exitFP,width=10)
  serch.but<-function(Optcrit){
    if (Optcrit=="A") but.A<-tkgrid(tkbutton(fiParFlow2,text="Search seq A-opt or near-opt design",command=function()optcrtF("A"),width=29), exit.but2a)
    if (Optcrit=="MV") but.A<-tkgrid(tkbutton(fiParFlow2,text="Search seq MV-opt or near-opt design",command=function()optcrtF("MV"),width=29), exit.but2a)
    if (Optcrit=="D") but.A<-tkgrid(tkbutton(fiParFlow2,text="Search seq D-opt or near-opt design",command=function()optcrtF("D"),width=29), exit.but2a)
    if (Optcrit=="E") but.A<-tkgrid(tkbutton(fiParFlow2,text="Search seq E-opt or near-opt design",command=function()optcrtF("E"),width=29), exit.but2a)
    return(but.A)
  }#end of serch.but
  tkgrid(tklabel(fiParFlow2,text="                                                    ",font=0.00001))
  serch.but(Optcrit)
  tkgrid(tklabel(fiParFlow2,text="                                             ",font=0.00001))
  tkgrid(fiParFup)
  tkgrid(tklabel(fiParF,text="---------------------------------", font=fontHeading3))
  tkgrid(fiParFlow)
  tkgrid(tklabel(fiParF,text="---------------------------------", font=fontHeading3))
  tkgrid(fiParFmid)
  tkgrid(tklabel(fiParF,text="---------------------------------", font=fontHeading3))
  tkgrid(fiParFlow2)
  tkgrid(fiParF) 
}#end of fixparsoptd.mae

