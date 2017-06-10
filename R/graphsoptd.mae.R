#Function for the plot of the graphical layout of resultant sequential optimal designs (graphsoptd.mae)
graphsoptd.mae<-function(trt.N, blk.N,theta,soptdesF,Optcrit,strt,sary, dtype) {
  trt.N<-trt.N+strt;
  blk.N<-blk.N+sary;
  trtblkthetano<-paste("(",paste(trt.N, blk.N, theta,sep=", "),")",sep="")
  NOptcrtr<-paste(Optcrit,"-optimal",sep="")#name of optimality criteria
  NOptcrtrG<-paste("Graph_layout_Seq",Optcrit,"opt",dtype,"_array",sep="")#name of folder where the graphical layout will be saved
  NOptcrtrG2<-paste("_GoutSeq",Optcrit,"opt",dtype,"array.pdf",sep="")
  NgoutT=paste("Sequential ", NOptcrtr,if(dtype=="rcd") {"row-column"} else{"block"}, "design for", paste("(v, b, theta) =",trtblkthetano,sep=" "))
  graph.des<-make_graph(as.numeric(as.factor(soptdesF)), directed = if(dtype=="rcd") {TRUE} else {FALSE}) %>%
    set_edge_attr("color", value = "black")# %>%
  weight0<-rep(c(0.5,1), c(blk.N-sary, sary))
  E(graph.des)$weight<-weight0
  E(graph.des)$width <- 2
  E(graph.des)[weight0 >0.5] $width <- 1
  E(graph.des)[weight0 >0.5] $color <- "red"
  V(graph.des)$color <- ifelse(V(graph.des)>(trt.N-strt),"brown", "cyan")
  graph.desid <- tkplot(graph.des, canvas.width=515, canvas.height=500,layout=layout.kamada.kawai)#,edge.color="black")
  canvas <- tk_canvas(graph.desid)
  padding <- 100
  coords <- norm_coords(layout=layout.kamada.kawai(graph.des), 0+padding, 450-padding,
                        50+padding, 500-padding)
  tk_set_coords(graph.desid, coords)
  width <- as.numeric(tkcget(canvas, "-width"))
  height <- as.numeric(tkcget(canvas, "-height"))
  tkcreate(canvas, "text", width/2, 20, text=NgoutT,
           justify="center", font=tcltk::tkfont.create(family="helvetica",size=10,weight="bold"))
  tkcreate(canvas, "text", width/2, 35, text="using array exchange algorithim",
           justify="center", font=tcltk::tkfont.create(family="helvetica",size=10,weight="bold"))
  graph.OutlayoptBlk<-paste(getwd(), NOptcrtrG,sep="/")
  if(!file.exists(graph.OutlayoptBlk)) dir.create(graph.OutlayoptBlk)
  obtdes.goutloptBlk<-paste(graph.OutlayoptBlk,paste(trtblkthetano,NOptcrtrG2,sep=""),sep="/")
  pdf(file=obtdes.goutloptBlk)
  plot(graph.des,edge.arrow.size=0.6, vertex.size=20, margin=0.5,
       layout=layout.kamada.kawai)
  title(paste("Graphical layout of ",  "Sequential ", Optcrit,"-optimal or near-optimal ",if(dtype=="rcd") {"row-column"} else{"block"}," design for:",sep=""), 
        sub = NULL,cex.main = 0.7,   font.main= 1, col.main= "black")
  mtext(paste("(v, b, theta) =", " (",paste(trt.N, blk.N, theta,sep=", "),")",sep=""), line = 0.50, cex=0.7, col = "blue", font = 1)
  
  grphlt0<-make_graph(as.numeric(as.factor((soptdesF[,1:(blk.N-sary)]))), directed = if(dtype=="rcd") {TRUE} else {FALSE})
  plot(grphlt0,edge.arrow.size=0.6, vertex.size=20, margin=0.5,
       layout=layout.kamada.kawai,vertex.color="cyan",edge.color="black")
  title(paste("Graphical layout of ",  "initial ",Optcrit, "-optimal or near-optimal ",if(dtype=="rcd") {"row-column"} else{"block"}," design",sep=""), 
        sub = NULL,cex.main = 0.7,   font.main= 1, col.main= "black")
  mtext(paste("(v, b, theta) =", " (",paste((trt.N-strt), (blk.N-sary), theta, sep=", "),")",sep=""), line = 0.50, col = "blue", cex=0.7, font = 1)
  dev.off()  
  file_loc<-obtdes.goutloptBlk
  file_loc2<-paste("Graphical layout of obtained",  "Sequential ", NOptcrtr, "or near-optimal ",if(dtype=="rcd") {"row-column"} else{"block"}," design is also saved in .pdf at:",sep=" ")
  cat(file_loc2,"\n",file_loc,"\n\n")
}#End Section 5 (plot of the graphical layout of resultant sequential optimal design, graphsoptd.mae)
