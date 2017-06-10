#Function for main menu TCL/TK window
mmenusoptd.mae<-function(){
  fontHeading<- tkfont.create(family="times",size=40,weight="bold")#,slant="italic")
  fontHeading3<-tkfont.create(family="times",size=10,weight="bold")
  AMVDEopt.top<-tktoplevel()
  tkwm.title(AMVDEopt.top,"Sequential Optimal Designs For Microarray Experiments")
  tkgrid(tklabel(AMVDEopt.top,text="     soptdmaeA 1.0.0     ",font=fontHeading))
  tkgrid(tklabel(AMVDEopt.top,text="",font=fontHeading))
  Fixp.butA<-tkbutton(AMVDEopt.top,text="Sequential A-Optimal Design",font=fontHeading3,command=function()fixparsoptd.mae("A"),width=30)
  Fixp.butMV<-tkbutton(AMVDEopt.top,text="   Sequential MV-Optimal Design",font=fontHeading3,command=function() fixparsoptd.mae("MV"),width=30)
  Fixp.butD<-tkbutton(AMVDEopt.top,text="Sequential D-Optimal Design",font=fontHeading3,command=function() fixparsoptd.mae("D"),width=30)
  Fixp.butE<-tkbutton(AMVDEopt.top,text="Sequential E-Optimal Design",font=fontHeading3,command=function() fixparsoptd.mae("E"),width=30)
  exitMM<-function() {
    closeQ=tkmessageBox(title = "Bye...", message = "Bye..., Enjoy your sequential optimal design",
                        icon = "info", type = "okcancel", default = "cancel")
    if(as.character(closeQ)=="ok") tkdestroy(AMVDEopt.top)
  }#end of exitMM
  ExitWin.but<-tkbutton(AMVDEopt.top,text="Exit",font=fontHeading3,command=exitMM,width=15)
  tkgrid(Fixp.butA)
  tkgrid(Fixp.butMV)
  tkgrid(Fixp.butD)
  tkgrid(Fixp.butE)
  tkgrid(ExitWin.but)
  tkgrid(tklabel(AMVDEopt.top,text="",font=fontHeading))   }#end of mmenusoptd.mae function
mmenusoptd.mae()

