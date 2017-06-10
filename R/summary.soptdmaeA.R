summary.soptdmaeA<-function(object,...)
{
seqoptbrcd_maes<-object
seqoptbrcd_maes$grphlt0<-make_graph(as.numeric(as.factor(object$optdes0)), directed = if(object$dtype=="rcd") {TRUE} else {FALSE})
seqoptbrcd_maes$grphlt<-make_graph(as.numeric(as.factor(object$soptdesF)), directed = if(object$dtype=="rcd") {TRUE} else {FALSE}) %>%
  set_edge_attr("color", value = "black")
weight0<-rep(c(0.5,1), c(object$b-object$sary, object$sary))
E(seqoptbrcd_maes$grphlt)$weight<-weight0
E(seqoptbrcd_maes$grphlt)$width <- 2
E(seqoptbrcd_maes$grphlt)[weight0 >0.5] $width <- 1
E(seqoptbrcd_maes$grphlt)[weight0 >0.5] $color <- "red"
V(seqoptbrcd_maes$grphlt)$color <- ifelse(V(seqoptbrcd_maes$grphlt)>(object$v-object$strt),"brown", "cyan")
class(seqoptbrcd_maes)<-"summary.soptdmaeA"
seqoptbrcd_maes
}
