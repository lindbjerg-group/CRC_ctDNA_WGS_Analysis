#Read and modify mutation data (also determines order of the cases and names)
tmbData <- read.csv("comut_TMB_forcomut-4816-6491-excluded.tsv",sep='\t')
tmbData <- tmbData[order(tmbData$TMB,decreasing=TRUE),]
tmbHeight <- max(log2(tmbData$TMB))-0.95*min(log2(tmbData$TMB))
tmbPlotVal <- (log2(tmbData$TMB)-0.95*min(log2(tmbData$TMB)))/tmbHeight
colTMB <- "#34516C"

#Read msi data
msiData <- read.csv("comut_MSIsensor_forcomut-4816-6491.csv")
msiData <- msiData[unlist(sapply(tmbData$sample,function(x) which(msiData$sample==x))),]
msiNorm <- 10
msiVal <- 1000*round(pmin(msiData$value/msiNorm,1),digits=3)+1
colList <- colorRampPalette(c("#8E383C","black"))
colRamp <- colList(1001)
colRamp2 <- colList(10001)
colMSI <- colRamp[msiVal]

#read MSI IHC status
statusData <- read.csv("comut_MSIclin_forcomut-4816-6492.csv",sep=';')
statusData <- statusData[unlist(sapply(tmbData$sample,function(x) which(statusData$PaperID==x))),]
colIHCList <- c("#395841","#8C4A5E","#B8BFA7")
colIHC <- colIHCList[sapply(statusData$MSI_histology,function(x) which(c(1,0,9)==x))]

#WGS data
wgdData <- read.csv("comut_WGD_forcomut-4816-6491-excluded.csv")
wgdData <- wgdData[unlist(sapply(tmbData$sample,function(x) which(wgdData$sample==x))),]
colWGDList <- c("#9E6A67","#8CB0C8")
colWGD <- colWGDList[sapply(wgdData$value,function(x) which(c("WGD","No WGD")==x))]

#read in mutational signature data
sigData <- read.csv("comut_signatures_forcomut-4816-6491-excluded.csv")
sigCols <- c('#3B99D4','#8ED14B','#F0B849','#b4dbed','#ba6e6e','#786846','#1E90FF','#aba7c7','#93bedf','#9de3cb','#C7A7BC','#bda46c')
sigNames <- c("Ageing","APOBEC","HRD","SBS5","MSI","SBS17","SBS18","ROS/MUTYH","SBS40","E.coli","SBS93","BER(NTHL1)")
sigData <- sigData[unlist(sapply(tmbData$sample,function(x) which(sigData$sample==x))),]
endVals <- sapply(1:nrow(sigData),function(x) cumsum(as.numeric(sigData[x,2:ncol(sigData)])))
startVals <- rbind(rep(0,ncol(endVals)),endVals[1:(nrow(endVals)-1),])

pdf("MOMA_Colon_Publication_Figure2.pdf",height=8,width=20)
plot(c(-20,nrow(tmbData)),c(-0.2,1.2),type='n',xaxt='n',yaxt='n',xlab='',ylab='',bty='n')
nullIn <- sapply(1:length(tmbPlotVal),function(x) rect(x-0.92,0.7,x-0.08,0.3*tmbPlotVal[x]+0.7,col=colTMB))
nullIn <- sapply(1:length(tmbPlotVal),function(x) rect(x-0.92,0.58,x-0.08,0.67,col=colMSI[x],border=NA))
nullIn <- sapply(1:length(tmbPlotVal),function(x) rect(x-0.92,0.46,x-0.08,0.55,col=colIHC[x],border=NA))
nullIn <- sapply(1:length(tmbPlotVal),function(x) rect(x-0.92,0.34,x-0.08,0.43,col=colWGD[x],border=NA))
for(i in 1:ncol(endVals)) {
  for(j in 1:nrow(endVals)) {
    if(endVals[j,i]!=startVals[j,i]) {
      rect(i-1,0.3*startVals[j,i],i,0.3*endVals[j,i],col=sigCols[j])
    }
  }
}
for(i in 1:length(tmbPlotVal)) {
  text(i-0.5,-0.005,labels=gsub("[*#]","",wgdData$sample[i]),adj=c(1,0.5),srt=90,cex=0.7)
  if(gsub("[*#]","",wgdData$sample[grep("[#]",wgdData$sample)])==gsub("[*#]","",wgdData$sample[i]))text(i-0.5,-0.15,"#",cex=0.8);
  if(gsub("[*]","",wgdData$sample[grep("[*]",wgdData$sample)])==gsub("[*#]","",wgdData$sample[i]))text(i-0.5,-0.15,"*",cex=1.8);
}
text(-12,0.2,"COSMIC",adj=c(0.5,0.5),cex=1.5)
text(-12,0.1,"signatures",adj=c(0.5,0.5),cex=1.5)
text(-12,0.34+0.09/2,"WGD",adj=c(0.5,0.5),cex=1.5)
text(-12,0.46+0.09/2,"MSI IHC status",adj=c(0.5,0.5),cex=1.5)
text(-12,0.58+0.09/2,"MSIsensor (WGS)",adj=c(0.5,0.5),cex=1.5)
text(-12,0.85,"TMB",adj=c(0.5,0.5),cex=1.5)
#Legend
text(29,1.2,"MSIsensor (WGS)",adj=c(0,0),cex=1)
nullIn <- sapply(seq(1,length(colRamp)-10,10),function(x) rect(29+((0.009)*(x-1)),1.16,29+((0.009)*(x+20)),1.10,col=colRamp[abs(1001-x)+1],border = NA))
rect(31.25-0.125,1.14,31.25+0.125,1.17,col='white',border=NA);rect(31.25-0.125,1.09,31.25+0.125,1.12,col='white',border=NA)
rect(33.5-0.125,1.14,33.5+0.125,1.17,col='white',border=NA);rect(33.5-0.125,1.09,33.5+0.125,1.12,col='white',border=NA)
rect(35.75-0.125,1.14,35.75+0.125,1.17,col='white',border=NA);rect(35.75-0.125,1.09,35.75+0.125,1.12,col='white',border=NA)
text(29,1.06,">10");text(29+4.5,1.06,"5");text(38,1.06,"0");
text(45,1.2,"MSI IHC status",adj=c(0,0),cex=1)
rect(45,1.15,47.5,1.11,col=colIHCList[1],border=NA);text(49,1.13,"MSI",adj=c(0,0.5))
rect(45,1.09,47.5,1.05,col=colIHCList[2],border=NA);text(49,1.07,"MSS",adj=c(0,0.5))
rect(45,1.03,47.5,0.99,col=colIHCList[3],border=NA);text(49,1.01,"Unknown",adj=c(0,0.5))
text(60,1.2,"WGD status",adj=c(0,0),cex=1)
rect(60,1.15,62.5,1.11,col=colWGDList[1],border=NA);text(64,1.13,"WGD",adj=c(0,0.5))
rect(60,1.09,62.5,1.05,col=colWGDList[2],border=NA);text(64,1.07,"No WGD",adj=c(0,0.5))
text(75,1.2,"COSMIC signatures",adj=c(0,0),cex=1)
rect(75,1.15,77.5,1.11,col=sigCols[1]);text(78.5,1.13,sigNames[1],adj=c(0,0.5))
rect(75,1.09,77.5,1.05,col=sigCols[2]);text(78.5,1.07,sigNames[2],adj=c(0,0.5))
rect(75,1.03,77.5,0.99,col=sigCols[3]);text(78.5,1.01,sigNames[3],adj=c(0,0.5))
rect(86,1.15,88.5,1.11,col=sigCols[4]);text(89.5,1.13,sigNames[4],adj=c(0,0.5))
rect(86,1.09,88.5,1.05,col=sigCols[5]);text(89.5,1.07,sigNames[5],adj=c(0,0.5))
rect(86,1.03,88.5,0.99,col=sigCols[6]);text(89.5,1.01,sigNames[6],adj=c(0,0.5))
rect(97,1.15,99.5,1.11,col=sigCols[7]);text(100.5,1.13,sigNames[7],adj=c(0,0.5))
rect(97,1.09,99.5,1.05,col=sigCols[8]);text(100.5,1.07,sigNames[8],adj=c(0,0.5))
rect(97,1.03,99.5,0.99,col=sigCols[9]);text(100.5,1.01,sigNames[9],adj=c(0,0.5))
rect(111,1.15,113.5,1.11,col=sigCols[10]);text(114.5,1.13,sigNames[10],adj=c(0,0.5))
rect(111,1.09,113.5,1.05,col=sigCols[11]);text(114.5,1.07,sigNames[11],adj=c(0,0.5))
rect(111,1.03,113.5,0.99,col=sigCols[12]);text(114.5,1.01,sigNames[12],adj=c(0,0.5))
text(128,1.2,"Synchronous tumors",adj=c(0,0),cex=1)
text(128,1.13,"*",adj=c(0,0.8),cex=2);text(129.5,1.13,"Patient CRC-043",adj=c(0,0.5))
text(128,1.07,"#",adj=c(0,0.5));text(129.5,1.07,"Patient CRC-047",adj=c(0,0.5))
dev.off()

