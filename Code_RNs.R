source("svmrfeFeatureRanking.R")
library("e1071")
library("igraph")

CSB_RNs <- function(x = x,y = y,network){
  #calculate the RNs score for each miRNA
  networkNode <- unique(c(network[,1],network[,2]))
  x <- x[networkNode,]
  x <- t(x)
  featureRankedList <- svmrfeFeatureRanking(x,y)
  
  g <- graph_from_data_frame(network,directed = F,vertices = networkNode)
  
  networki <- data.frame(row.names = networkNode,Degree = degree(g),AverageShortestPathLength = 1/closeness(g,normalized = T))
  
  
  geneOrdern<-data.frame(cbind(rank=featureRankedList,ID=colnames(x)))
  geneOrder<-geneOrdern[order(featureRankedList),]
  geneOrder<-cbind(geneOrder,score=round(((1+length(geneOrder[,1]))-as.numeric(as.vector(geneOrder$rank)))/length(geneOrder[,1]),5))
  networki <- networki[geneOrder$ID,]
  
  geneOrder <- data.frame(name=geneOrder$ID,
                          svmScore = geneOrder$score,
                          AVGP = networki$AverageShortestPathLength,
                          Degree = networki$Degree)
  rownames(geneOrder) <- geneOrder$name
  geneOrder$RNs <- geneOrder$Degree * geneOrder$svmScore / geneOrder$AVGP
  geneOrder <- geneOrder[order(geneOrder$RNs,decreasing = T),]
  return(geneOrder)
}

#network <- read.delim(file = "string_interactions_short.tsv",header = T)
network <- read.csv(file = "network_tcga3.csv",header = T)
x <- read.csv(file = "tcga_network_exp3.csv",header = T,row.names = 1)
y <- read.csv(file = "tcga_pdata.csv",header = T,row.names = 1)
x <- x[,rownames(y)]
y <- ifelse(y$phase=="NC",0,1)
CSB_RNs(x = x,y = y,network = network)

write.csv(geneOrder,file = "geneOrder3.csv")
