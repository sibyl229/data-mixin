library(lattice)

pDat <- read.table('input/Yeast/combined_features.tsv', sep='\t', header=TRUE, row.names=1)
dim(pDat)
pDat <- subset(pDat, subset=pDat$Half_Life < 300) # 300 marks stable proteins
cn <- colnames(pDat)

# feature distribution
for (colnX in cn){
  colv <- pDat[[colnX]]
  #print(histogram(~colv, xlab=colnX))
}

### remove rows which contains feature values marked as outlier
outlierRows <- rep(T,nrow(pDat))
for (i in cn){
  colv <- pDat[[i]]
  if(is.numeric(colv)){
    outvalues <- boxplot.stats(colv, do.out=T)$out
#     print(i)
#     print(outvalues)
    outlierRows <- outlierRows &&
      is.element(colv,outvalues)                 
  }
}
pDat <- pDat[!outlierRows,]
dim(pDat)
densityplot(~Half_Life, dat=pDat)


xyplot(totalDisorderAA ~ Half_Life, data=pDat)
xyplot(CAI ~ Half_Life, data=pDat)
cnX = colnames(pDat)
corX <- cor(pDat$Half_Life, pDat[cnX[cnX != 'Half_Life' & cnX !='nEnd']],
            method='pearson')
sortedColX = corX[,order(abs(corX))]
barchart(as.table(sortedColX), xlab='Pearson\'s correlation')
cor(pDat$CAI, pDat$CODON_BIAS)

# n-end vs half-life
barchart(table(pDat$nEnd))
countByNEnd  <- aggregate(rownames(pDat), by=list(pDat$nEnd), FUN=length)
colnames(countByNEnd) <-  c('nEnd', 'freq')
commNEndAA <- countByNEnd$nEnd[countByNEnd$freq>=15]
bwplot(~Half_Life|nEnd, 
       data=pDat,
       subset=is.element(pDat$nEnd, commNEndAA),
       xlim=c(0,300))
#bwplot(~Half_Life|nEnd, data=pDat, xlim=c(0,300))

