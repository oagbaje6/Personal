#!/usr/local/bin/Rscript


datafile <- "/courses/example-data/csv/MetaHIT_SangerSamples.genus.txt"
data=read.table(datafile, header=T, row.names=1, dec=".", sep="\t")
data=data[-1,]

clsfile <- "/courses/example-data/csv/MetaHIT_SangerSamples.genus.cls"
datacls=read.table(clsfile, header=T, row.names=1, dec=".", sep=",")


KLD <- function(x,y) sum(x * log(x/y))
JSD <- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))

dist.JSD <- function(inMatrix, pseudocount=0.000001, ...) {
    KLD <- function(x,y) sum(x *log(x/y))
    JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
    matrixColSize <- length(colnames(inMatrix))
    matrixRowSize <- length(rownames(inMatrix))
    colnames <- colnames(inMatrix)
    resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
    inMatrix = apply(inMatrix,1:2,function(x) ifelse (x==0,pseudocount,x))
    for(i in 1:matrixColSize) {
        for(j in 1:matrixColSize) { 
            resultsMatrix[i,j]=JSD(as.vector(inMatrix[,i]), as.vector(inMatrix[,j]))
        }
    }
    colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
    as.dist(resultsMatrix)->resultsMatrix
    attr(resultsMatrix, "method") <- "dist"
    return(resultsMatrix) 
}


data.dist = dist.JSD(data)
pam.clustering = function(x,k) { # x is a distance matrix and k the number of clusters
    require(cluster)
    cluster = as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)
    return(cluster)
}

data.cluster=pam.clustering(data.dist, k=3)
sampleList <- colnames(data)
n <- length(sampleList)

cat("sample,group\n")

for (i in 1:n){
    s = paste0(sampleList[i], ",", data.cluster[i], "\n")
    cat (s)
}
#cat(data.cluster)



