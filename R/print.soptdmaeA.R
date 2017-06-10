print.soptdmaeA<-function(x,...)
{
cat("\n        ---------------------------------------    \n")
cat("Title: Sequential ",x$Optcrit,"-optimal or near-optimal ",if(x$dtype=="rcd") {"row-column"} else{"block"} ," design        ","Date: ", format(Sys.time(), "%a %b %d %Y %H:%M:%S"),"\n",sep="")
cat("        ---------------------------------------    \n")
cat("Call:\n")
print(x$call)
cat("\nResultant Sequential",x$Optcrit,"-optimal or near-optimal ",if(x$dtype=="rcd") {"row-column"} else{"block"} ," design:\n",sep="")
cat("\n")
print(data.frame(x$soptdesF))
cat("\n", x$file_loc2,"\n", x$file_loc,"\n")
cat("\n")
}
