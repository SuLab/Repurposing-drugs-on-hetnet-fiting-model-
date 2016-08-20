rm(list = ls())
library("igraph")

fileontology <- read.csv("/home/nuria/workspace/repurposing-hetio/rephetio-dhimmelstein/Characterization-Disease-Drug-Scope-of-hetio/out/do_ontologyChildParent.tab", sep = "\t", header = F)
edges <- fileontology[c(1,2)]
nodes<- as.data.frame(unique(c(unique(as.character(edges[, 1])), unique(as.character(edges[,2])))))
dim(nodes)
graph.full = graph.data.frame(edges, directed=T, vertices=nodes)
vcount(graph.full) # number of nodes
ecount(graph.full) # number of edges

graph <- igraph::simplify(graph.full)
vcount(graph)  
ecount(graph) 

resultsPathways  <- data.frame(t(rep(NA,3)))
colnames(resultsPathways ) <- c("child", "parent", "length")

fileparent <- read.csv("/home/nuria/workspace/repurposing-hetio/rephetio-dhimmelstein/Characterization-Disease-Drug-Scope-of-hetio/out/do_ontology1rstBranchParents.tab", sep = "\t", header = F)

parent <- fileparent[c(1)]

filechild <- read.csv("/home/nuria/workspace/repurposing-hetio/rephetio-dhimmelstein/Characterization-Disease-Drug-Scope-of-hetio/out/hetionet-doids-list.tab", sep = "\t", header = F)
child <- filechild[c(1)]
child <- child$V1
parent <- parent$V1
length(unique(child))
# 137
colnames(nodes) <- "node"
child <- intersect(child,nodes$node)
length(unique(child))
#  136

c <- 1  
for (nodo1 in child ){     
  for (nodo2 in parent){       
    sp <- get.shortest.paths(graph,from = nodo1, to = nodo2, mode = "out")
    
    if (length(sp$vpath[[1]]) > 0){
      l <- length(sp$vpath[[1]])
      resultsPathways[c, 1] <- nodo1
      resultsPathways[c, 2] <- nodo2      
      resultsPathways[c, 3] <- l # original row
      c<- c+1
    }
  }
}

resultsPathways
write.csv(resultsPathways,"/home/nuria/workspace/repurposing-hetio/rephetio-dhimmelstein/Characterization-Disease-Drug-Scope-of-hetio/out/hetionet-doidsPathways.csv")
str(resultsPathways)
length(unique(parent))
length(unique(child))
doClass.freq = table(sort(resultsPathways$parent))
sum(doClass.freq) # 136
hetio.freq.df <- as.data.frame(doClass.freq)
colnames(hetio.freq.df)[1] = 'doid'
colnames(fileparent) = c('doid','name')
hetio.freq.df.merge = merge(fileparent,hetio.freq.df, by.x = 'doid', all.x = T)
hetio.freq.df.merge$Freq[is.na(hetio.freq.df.merge$Freq)] = 0




doClassPercent.freq = prop.table(doClass.freq)*100
doClassPercentSort.freq = sort(doClassPercent.freq,decreasing = T)
#sortnames <- c("Disease of anatomical entity","Disease of cellular proliferation","Disease by infectious agent","Disease of metabolism","Genetic disease","Disease of mental health","Syndrome","Physical disorder")
sortnames <- c(fileparent[c(2)])
names(doClassPercentSort.freq) <- sortnames$V2

# Save graphic
library(Cairo)
#CairoSVG("/ibi/users/shared/textmining_testing/befree_filtered/Figure1.svg", width=13, height=6)
#par(mai=c(1.,1.,1.,.2))
png("/home/nuria/workspace/repurposing-hetio/rephetio-dhimmelstein/Characterization-Disease-Drug-Scope-of-hetio/out/hetionet-doids-classes-percent.png",width=1700,height=960, units = "px", pointsize = 18)
par(oma=c(.2,8.,.2,.2))
barplot(doClassPercentSort.freq,las=2,cex.names=.8,xlim=c(0,100),main="DO Classes",xlab="% of diseases",cex.axis=.8,horiz=T)
text(doClassPercentSort.freq+2,f1bp,labels=paste0(as.character(round(doClassPercentSort.freq)),"%"),cex=.7,xpd=T)
f1bp <- barplot(doClassPercentSort.freq,las=2,cex.names=.8,xlim=c(0,100),main="DO Classes",xlab="% of diseases",cex.axis=.8,horiz=T)
text(doClassPercentSort.freq+2,f1bp,labels=paste0(as.character(round(doClassPercentSort.freq)),"%"),cex=.7,xpd=T)
dev.off()

DO1.table <- data.frame(cbind(sort(doClass.freq, decreasing = T)))
DO1percent.table <- data.frame(cbind(sort(doClassPercent.freq, decreasing = T)))
colnames(DO1percent.table) <- c("Percentage")
rownames(DO1percent.table) <- sortnames
write.csv(DO1.table,"/home/nqueralt/Desktop/MedBioinformatics-phenotypic-diseasome-2015/data/out/do-cp-classes-freqTable.csv")
write.csv(DO1percent.table,"/home/nqueralt/Desktop/MedBioinformatics-phenotypic-diseasome-2015/data/out/do-cp-classes-freqPercentTable.csv")


# Add curated data and ontology

rm(list = ls())
library("igraph")


resultsPathways  <- data.frame(t(rep(NA,3)))
colnames(resultsPathways ) <- c("child", "parent", "length")


fileparent <- read.csv("/home/nqueralt/Desktop/MedBioinformatics-phenotypic-diseasome-2015/data/out/do_ontology1rstBranchParents.tab", sep = "\t", header = F)
parent <- fileparent[c(1)]
#filechildcurated <- read.csv("/home/nqueralt/Desktop/MedBioinformatics-phenotypic-diseasome-2015/data/out/stats_hpCuratedDiseases.tab", sep = "\t", header = F)
filechildcurated <- read.csv("/home/nqueralt/Desktop/MedBioinformatics-phenotypic-diseasome-2015/data/out/stats_hpCuratedDiseasesEvidence-refined.tab", sep = "\t", header = F)
child <- filechildcurated[c(1)]
child <- child$V1
parent <- parent$V1
length(unique(child))
# 1317
colnames(nodes) <- "node"
child <- intersect(child,nodes$node)
length(unique(child))
#  1317

c <- 1  
for (nodo1 in child ){     
  for (nodo2 in parent){       
    sp <- get.shortest.paths(graph,from = nodo1, to = nodo2, mode = "out")
    
    if (length(sp$vpath[[1]]) > 0){
      l <- length(sp$vpath[[1]])
      #parent <-labels( graph[[sp[[1]][length(sp[[1]])]]])  
      
      resultsPathways[c, 1] <- nodo1
      resultsPathways[c, 2] <- nodo2      
      #resultsPathways[c, 3] <- l-1 # count of edges
      resultsPathways[c, 3] <- l # original row
      #       camino <- sp$vpath[[1]]
      #       caminolabels <-labels( graph[[camino]])  
      #       caminolabels[l]
      #       paste(caminolabels, collapse = ", ")
      c<- c+1
    }
  }
}

resultsPathwaysCurated <- resultsPathways
resultsPathwaysCurated
write.csv(resultsPathwaysCurated,"/home/nqueralt/Desktop/MedBioinformatics-phenotypic-diseasome-2015/data/out/do_resultsPathwaysCurated.csv")
str(resultsPathwaysCurated)
length(unique(parent))
length(unique(child))
doClass.freq = table(sort(resultsPathways$parent))
sum(doClass.freq) # 1324
doClassPercent.freq = prop.table(doClass.freq)*100
doClassPercentSort.freq = sort(doClassPercent.freq,decreasing = T)
sortnames <- c("Disease of anatomical entity","Disease of cellular proliferation","Disease by infectious agent","Disease of metabolism","Genetic disease","Disease of mental health","Syndrome","Physical disorder")
names(doClassPercentSort.freq) <- sortnames

doClassPercentSortCurated.freq <- doClassPercentSort.freq

# Save graphic
library(Cairo)
#CairoSVG("/ibi/users/shared/textmining_testing/befree_filtered/Figure1.svg", width=13, height=6)
#par(mai=c(1.,1.,1.,.2))
png("/home/nqueralt/Desktop/MedBioinformatics-phenotypic-diseasome-2015/data/out/do-classes-percent2-curated.png",width=1700,height=960, units = "px", pointsize = 18)
par(oma=c(.2,8.,.2,.2))
barplot(doClassPercentSortCurated.freq,las=2,cex.names=.8,xlim=c(0,100),main="DO Classes",xlab="% of diseases",cex.axis=.8,horiz=T)
text(doClassPercentSortCurated.freq+2,f1bp,labels=paste0(as.character(round(doClassPercentSortCurated.freq)),"%"),cex=.7,xpd=T)
f1bp <- barplot(doClassPercentSortCurated.freq,las=2,cex.names=.8,xlim=c(0,100),main="DO Classes",xlab="% of diseases",cex.axis=.8,horiz=T)
text(doClassPercentSortCurated.freq+2,f1bp,labels=paste0(as.character(round(doClassPercentSortCurated.freq)),"%"),cex=.7,xpd=T)
dev.off()


# Add second percentage: percentage regarging total DOID

rm(list = ls())
library("igraph")


resultsPathways  <- data.frame(t(rep(NA,3)))
colnames(resultsPathways ) <- c("child", "parent", "length")


fileparent <- read.csv("/home/nqueralt/Desktop/MedBioinformatics-phenotypic-diseasome-2015/data/out/do_ontology1rstBranchParents.tab", sep = "\t", header = F)
parent <- fileparent[c(1)]
filechildDO <- read.csv("/home/nuria/workspace/repurposing-hetio/rephetio-dhimmelstein/Characterization-Disease-Drug-Scope-of-hetio/out/allActIdDOOntology.tab", sep = "\t", header = F)
child <- filechildDO[c(1)]
child <- child$V1
parent <- parent$V1
length(unique(child))
#  6824
colnames(nodes) <- "node"
child <- intersect(child,nodes$node)
length(unique(child))
# 6824

c <- 1  
for (nodo1 in child ){     
  for (nodo2 in parent){       
    sp <- get.shortest.paths(graph,from = nodo1, to = nodo2, mode = "out")
    
    if (length(sp$vpath[[1]]) > 0){
      l <- length(sp$vpath[[1]])
      #parent <-labels( graph[[sp[[1]][length(sp[[1]])]]])  
      
      resultsPathways[c, 1] <- nodo1
      resultsPathways[c, 2] <- nodo2      
      #resultsPathways[c, 3] <- l-1 # count of edges
      resultsPathways[c, 3] <- l # original row
      #       camino <- sp$vpath[[1]]
      #       caminolabels <-labels( graph[[camino]])  
      #       caminolabels[l]
      #       paste(caminolabels, collapse = ", ")
      c<- c+1
    }
  }
}

resultsPathwaysDO <- resultsPathways
resultsPathwaysDO
write.csv(resultsPathwaysDO,"/home/nqueralt/Desktop/MedBioinformatics-phenotypic-diseasome-2015/data/out/do_resultsPathwaysDO.csv")
str(resultsPathwaysDO)
length(unique(parent))
length(unique(child))
doClass.freq = table(sort(resultsPathways$parent))
sum(doClass.freq) # 6849
doClassPercent.freq = prop.table(doClass.freq)*100
doClassPercentSort.freq = sort(doClassPercent.freq,decreasing = T)
sortnames <- c("Disease of anatomical entity","Disease of cellular proliferation","Disease by infectious agent","Disease of metabolism","Genetic disease","Disease of mental health","Syndrome","Physical disorder")
names(doClassPercentSort.freq) <- sortnames

doClassPercentSortDO.freq <- doClassPercentSort.freq

# Save graphic
library(Cairo)
#CairoSVG("/ibi/users/shared/textmining_testing/befree_filtered/Figure1.svg", width=13, height=6)
#par(mai=c(1.,1.,1.,.2))
png("/home/nqueralt/Desktop/MedBioinformatics-phenotypic-diseasome-2015/data/out/do-classes-percent2-doOntology.png",width=1700,height=960, units = "px", pointsize = 18)
par(oma=c(.2,8.,.2,.2))
barplot(doClassPercentSortDO.freq,las=2,cex.names=.8,xlim=c(0,100),main="DO Classes",xlab="% of diseases",cex.axis=.8,horiz=T)
text(doClassPercentSortDO.freq+2,f1bp,labels=paste0(as.character(round(doClassPercentSortDO.freq)),"%"),cex=.7,xpd=T)
f1bp <- barplot(doClassPercentSortDO.freq,las=2,cex.names=.8,xlim=c(0,100),main="DO Classes",xlab="% of diseases",cex.axis=.8,horiz=T)
text(doClassPercentSortDO.freq+2,f1bp,labels=paste0(as.character(round(doClassPercentSortDO.freq)),"%"),cex=.7,xpd=T)
dev.off()

# Save graphic
library(Cairo)
#CairoSVG("/ibi/users/shared/textmining_testing/befree_filtered/Figure1.svg", width=13, height=6)
#par(mai=c(1.,1.,1.,.2))
png("/home/nqueralt/Desktop/MedBioinformatics-phenotypic-diseasome-2015/data/out/do-classes-percent2.png",width=1700,height=960, units = "px", pointsize = 18)
par(oma=c(.2,12.,.2,.2))
barplot(doClassPercentSortDO.freq,las=2,cex.names=.8,xlim=c(0,100),main="DO Classes",xlab="% of diseases",cex.axis=.8,horiz=T, beside=T)
text(doClassPercentSortDO.freq+2,f1bp,labels=paste0(as.character(round(doClassPercentSortDO.freq)),"%"),cex=.7,xpd=T)
f1bp <- barplot(doClassPercentSortDO.freq,las=2,cex.names=.8,xlim=c(0,100),main="DO Classes",xlab="% of diseases",cex.axis=.8,horiz=T)
text(doClassPercentSortDO.freq+2,f1bp,labels=paste0(as.character(round(doClassPercentSortDO.freq)),"%"),cex=.7,xpd=T)
barplot(doClassPercentSort.freq,col="blue",add=T,las=2,cex.names=.8,xlim=c(0,100),main="DO Classes",xlab="% of diseases",cex.axis=.8,horiz=T)
text(doClassPercentSort.freq+2,f1bp,labels=paste0(as.character(round(doClassPercentSort.freq)),"%"),cex=.7,xpd=T)
f1bp <- barplot(doClassPercentSort.freq,col="blue",add=T,las=2,cex.names=.8,xlim=c(0,100),main="DO Classes",xlab="% of diseases",cex.axis=.8,horiz=T)
text(doClassPercentSort.freq+2,f1bp,labels=paste0(as.character(round(doClassPercentSort.freq)),"%"),cex=.7,xpd=T)
barplot(doClassPercentSortCurated.freq,col="red",add=T,las=2,cex.names=.8,xlim=c(0,100),main="DO Classes",xlab="% of diseases",cex.axis=.8,horiz=T)
text(doClassPercentSortCurated.freq+2,f1bp,labels=paste0(as.character(round(doClassPercentSortCurated.freq)),"%"),cex=.7,xpd=T)
f1bp <- barplot(hpClassPercentSortCurated.freq,col="red",add=T,las=2,cex.names=.8,xlim=c(0,100),main="DO Classes",xlab="% of diseases",cex.axis=.8,horiz=T)
text(doClassPercentSortCurated.freq+2,f1bp,labels=paste0(as.character(round(doClassPercentSortCurated.freq)),"%"),cex=.7,xpd=T)
dev.off()

DO1.table <- data.frame(cbind(sort(doClass.freq, decreasing = T)))
DO1percent.table <- data.frame(cbind(sort(doClassPercent.freq, decreasing = T)))
colnames(DO1percent.table) <- c("Percentage")
rownames(DO1percent.table) <- sortnames
write.csv(DO1.table,"/home/nqueralt/Desktop/MedBioinformatics-phenotypic-diseasome-2015/data/out/do-classes-freqTableDO.csv")
write.csv(DO1percent.table,"/home/nqueralt/Desktop/MedBioinformatics-phenotypic-diseasome-2015/data/out/do-classes-freqPercentTableDO.csv")

