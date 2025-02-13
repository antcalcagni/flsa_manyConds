
#setwd("/home/antonio/MEGA/Lavoro_sync/Conferences/2025_IES Bressanone/contribution_2/Rcode/")

# Initial settings --------------------------------------------------------
rm(list=ls()); graphics.off(); options(warn = -1)
source("utils.R"); source("FLSA.R"); source("kneePoint.R")
library(tm); library(quanteda); library(stringr)



# Data --------------------------------------------------------------------
txtin <- read.csv(file = "airfrance_tripadvisor_reviews.csv",header = TRUE)
str(txtin)

yrs <- unlist(lapply(strsplit(x = txtin$publishedDate,split = "-"),function(x)x[1]))
table(yrs)

# Note: we are working with a subset of the original data, which includes two years before and after the pandemic
iid <- yrs=="2018" | yrs=="2019" | yrs=="2023" | yrs=="2024" 
txtin <- txtin[iid,]
corpus <- tm::Corpus(VectorSource(txtin$text)); N=length(corpus)
dtm <- tm::DocumentTermMatrix(corpus); dtm_matrix = as.matrix(dtm) #raw corpus

corpus <- tm::tm_map(corpus,content_transformer(tolower))
corpus <- tm::tm_map(corpus,removeNumbers)
corpus <- tm::tm_map(corpus,removePunctuation) 
corpus <- tm::tm_map(corpus,removeWords, stopwords("english"))
corpus <- tm::tm_map(corpus,stripWhitespace)

dtm <- tm::DocumentTermMatrix(corpus); dtm_matrix = as.matrix(dtm)

out = apply(dtm_matrix,2,sum); 
Vocab = data.frame(term=names(out),freq=as.numeric(out),cums=cumsum(as.numeric(out))); Vocab = Vocab[order(Vocab$freq,decreasing = FALSE),];
Vocab$freq_rel = round(Vocab$freq/sum(Vocab$freq),4); Vocab$cums_rel = round(Vocab$cums/sum(Vocab$freq),4)
Vocab = Vocab[order(Vocab$freq,decreasing = TRUE),]; head(Vocab)

summary(Vocab$freq)
jjd <- c(as.numeric(which(apply(dtm_matrix,2,sum)<9))) #remove rare words
dtm_matrix_filtered = dtm_matrix[,setdiff(1:NROW(Vocab),jjd)] #direction: keep

out = apply(dtm_matrix_filtered,2,sum); 
Vocab_filtered = data.frame(term=names(out),freq=as.numeric(out)); Vocab_filtered = Vocab_filtered[order(Vocab_filtered$freq,decreasing = TRUE),]; head(Vocab_filtered)
M = NROW(Vocab_filtered); 
summary(Vocab_filtered$freq)

dtm_filtered <- as.DocumentTermMatrix(dtm_matrix_filtered,weighting = tm::weightTf)
words <- colnames(dtm_filtered)

out_pcp <- Rdimtools::do.rpca(X = dtm_matrix_filtered,abstol=1e-4,maxiter=3e2)
L_matrix <- out_pcp$L; colnames(L_matrix) <- colnames(dtm_matrix_filtered); rownames(L_matrix) <- rownames(dtm_matrix_filtered)
save(dtm_matrix_filtered,dtm_filtered,words,L_matrix,file = "casestudy_data.RData")



# Results -----------------------------------------------------------------
load("results/dataout_aggregated.RData")


### Figure 1
cls = c("#36648B","#CD4F39")
lwdx=2
nocs <- seq(from=2,to=32); nocs2 <- (seq(from=2,to=32,length=32/2))
mains <- c("(A) tSVD","(B) NNMF","(C) LE")

tikzDevice::tikz(file='../contribution/fig1.tex',width=6,height=4,sanitize = TRUE)
par(mfrow=c(2,3),mar = c(5, 5, 3, 2))

## First row of the plot
Z1 <- cbind(STD_results[1,,"svd","J=2"],STD_results[1,,"nnmf","J=2"],STD_results[1,,"le","J=2"])
Z2 <- cbind(PCP_results[1,,"svd","J=2"],PCP_results[1,,"nnmf","J=2"],PCP_results[1,,"le","J=2"])
Z2[30,3] <- mean(Z2[,2],na.rm=TRUE)
ymin <- min(min(Z1,na.rm=TRUE),min(Z2,na.rm=TRUE)); ymax <- max(max(Z1,na.rm=TRUE),max(Z2,na.rm=TRUE))

for(j in 1:3){
  px1=kneePoint(x = nocs,y = Z1[,j],plot = FALSE,df = 1,bty="n",sign = -1,xQuery = nocs)
  px2=kneePoint(x = nocs,y = Z2[,j],plot = FALSE,df = 1,bty="n",sign = -1,xQuery = nocs)
  
  #plot(nocs,smooth.spline(Z1[,j],df = median(Z1[,j]))$y,bty="n",type="l",lwd=lwdx,ylim = c(ymin,ymax),xlab="",ylab="",col=cls[1],axes = FALSE)
  #points(nocs,smooth.spline(Z2[,j],df = median(Z2[,j]))$y,bty="n",type="l",lwd=lwdx,col=cls[2])
  plot(nocs,Z1[,j],bty="n",type="l",lwd=lwdx,ylim = c(ymin,ymax),xlab="",ylab="",col=cls[1],axes = FALSE)
  points(nocs,Z2[,j],bty="n",type="l",lwd=lwdx,col=cls[2])
  title(main = mains[j],adj=0,line=0.5,cex.main=1.5)
  axis(side = 1,at = nocs2,labels = nocs2,cex.axis=1.25); axis(side = 2,at = round(seq(from=ymin,to=ymax,length=9),2),cex.axis=1.25)
  if(j==1){mtext(side = 2,text = "UMASS-like",cex=1.45,padj = -2)}
  
  abline(v = c(px1,px2),col=cls,lty=2,lwd=1.25); 
  text(c(px1)+2,c(Z1[1,j]),labels = c(px1),cex = 1.25,col=cls[1]); text(c(px2)-2,c(Z2[1,j])+0.2,labels = px2,cex = 1.25,col=cls[2])
}

## Second row of the plot
Z1 <- cbind(STD_results[2,,"svd","J=2"],STD_results[2,,"nnmf","J=2"],STD_results[2,,"le","J=2"])
Z2 <- cbind(PCP_results[2,,"svd","J=2"],PCP_results[2,,"nnmf","J=2"],PCP_results[2,,"le","J=2"])
Z2[30,3] <- mean(Z2[,2],na.rm=TRUE)
ymin <- min(min(Z1,na.rm=TRUE),min(Z2,na.rm=TRUE)); ymax <- max(max(Z1,na.rm=TRUE),max(Z2,na.rm=TRUE))

for(j in 1:3){
  px1=kneePoint(x = nocs,y = Z1[,j],plot = FALSE,df = 1,bty="n",sign = 1,xQuery = nocs)
  px2=kneePoint(x = nocs,y = Z2[,j],plot = FALSE,df = 1,bty="n",sign = 1,xQuery = nocs)
  
  plot(nocs,Z1[,j],bty="n",type="l",lwd=lwdx,ylim = c(ymin,ymax),xlab="",ylab="",col=cls[1],axes = FALSE)
  points(nocs,Z2[,j],bty="n",type="l",lwd=lwdx,col=cls[2])
  #title(main = mains[j],adj=0,line=0.5,cex.main=1.5)
  axis(side = 1,at = nocs2,labels = nocs2,cex.axis=1.25); axis(side = 2,at = round(seq(from=ymin,to=ymax,length=9),2),cex.axis=1.25)
  if(j==1){mtext(side = 2,text = "UCI-like",cex=1.45,padj = -2)}
  
  abline(v = c(px1,px2),col=cls,lty=2,lwd=1.25); 
  text(c(px1)+2,c(Z1[1,j]),labels = c(px1),cex = 1.25,col=cls[1]); text(c(px2)-2,c(Z2[1,j])+0.2,labels = px2,cex = 1.25,col=cls[2])
}
add_legend("bottom",legend = c("STD","PCP"),fill = cls,border = FALSE,bty = "n",cex=1.5,ncol=2)
dev.off()

best_nocs <- c(19,14) #STD,PCP (according to SVD)

### Figure 2
tikzDevice::tikz(file='../contribution/fig2.tex',width=9,height=4,sanitize = TRUE)
par(mfrow=c(1,2))
lwdx <- 2; nocs <- seq(from=2,to=32); nocs2 <- (seq(from=2,to=32,length=32/2))

Z1 <- cbind(STD_results[4,,"svd","J=2"],PCP_results[4,,"svd","J=2"])
plot(nocs,Z1[,1],type="l",bty="n",ylim=c(min(Z1),max(Z1)),xlab="",ylab="",axes=FALSE,col=cls[1],lwd=lwdx)
lines(nocs,Z1[,2],col=cls[2],lwd=lwdx)
axis(side = 1,at = nocs2,labels = nocs2,cex.axis=1.25); axis(side = 2,at = round(seq(from=min(Z1),to=max(Z1),length=10),2),cex.axis=1.25)
abline(v = best_nocs,col=cls,lty=2,lwd=1.25); 
text(best_nocs[1]+1,Z1[best_nocs[1],1]+0.05,labels = best_nocs[1],cex = 1.25,col=cls[1]); text(best_nocs[2]-1,Z1[best_nocs[2],2]+0.05,labels = best_nocs[2],cex = 1.25,col=cls[2])
#add_legend("topright",legend = c("STD","PCP"),fill = cls,border = FALSE,bty = "n",cex=1.5,ncol=1)
title(main = "(A) CL",adj=0,line=0.5,cex.main=1.5)

Z1 <- cbind(STD_results[5,,"svd","J=2"],PCP_results[5,,"svd","J=2"])
plot(nocs,Z1[,1],type="l",bty="n",ylim=c(min(Z1),max(Z1)),xlab="",ylab="",axes=FALSE,col=cls[1],lwd=lwdx)
lines(nocs,Z1[,2],col=cls[2],lwd=lwdx)
axis(side = 1,at = nocs2,labels = nocs2,cex.axis=1.25); axis(side = 2,at = round(seq(from=min(Z1),to=max(Z1),length=10),2),cex.axis=1.25)
abline(v = best_nocs,col=cls,lty=2,lwd=1.25); 
text(best_nocs[1]+1,Z1[best_nocs[1],1]+0.05,labels = best_nocs[1],cex = 1.25,col=cls[1]); text(best_nocs[2]-1,Z1[best_nocs[2],2]+0.05,labels = best_nocs[2],cex = 1.25,col=cls[2])
title(main = "(B) RL",adj=0,line=0.5,cex.main=1.5)
add_legend("topright",legend = c("STD","PCP"),fill = cls,border = FALSE,bty = "n",cex=1.5,ncol=1)
dev.off()


### Run _STD and _PCP based fLSA using the selected number of topics
load("casestudy_data.RData")
mod1_std = FLSA_core(X=dtm_matrix_filtered,N = NROW(dtm_matrix_filtered),M = NCOL(dtm_matrix_filtered),K = 19,type = "svd",J_svd = 2)

# Inter-topic distances via h-clustering
distx <- text2vec::dist2(x = t(mod1_std$PW_T),method="cosine") #higher similarity scores indicate greater redundancy between topics
hc <- hclust(d = as.dist(distx),method = "average")
silh <- sapply(2:19,function(x)clusterSim::index.S(d=distx,cl=cutree(hc,x)))  
which.max(silh)+1 #best number of topics: 2

### Run _STD based fLSA on the aggregated topic solution
txtin <- read.csv(file = "airfrance_tripadvisor_reviews.csv",header = TRUE)
yrs <- unlist(lapply(strsplit(x = txtin$publishedDate,split = "-"),function(x)x[1]))
iid <- yrs=="2018" | yrs=="2019" | yrs=="2023" | yrs=="2024" 
yrs <- yrs[iid]
yrs_num <- yrs; yrs_num[yrs_num=="2018"] <- 1; yrs_num[yrs_num=="2019"] <- 2; yrs_num[yrs_num=="2023"] <- 3; yrs_num[yrs_num=="2024"] <- 4
yrs_num <- as.numeric(yrs_num)

mod2_std = FLSA_core(X=dtm_matrix_filtered,N = NROW(dtm_matrix_filtered),M = NCOL(dtm_matrix_filtered),K = 2,type = "svd",J_svd = 2)

PW_T <- mod2_std$PW_T
PW_T_frex = stm::calcfrex(logbeta = log(t(PW_T)),wordcounts = apply(dtm_matrix_filtered,2,sum)); rownames(PW_T_frex) = rownames(PW_T)
Frex_sorted_std = mapply(function(j)names(PW_T_frex[order(PW_T_frex[,j],decreasing = TRUE),j][1:30]),1:2) 
colnames(Frex_sorted_std)=paste0("topic",1:2)

print(Frex_sorted_std)
# Topic 1: Customer Experience & Travel Preferences
# The words in this topic suggest discussions related to customer experiences with flights, services, and travel preferences. 
# Words like "spend", "details", "happy", and "pleased" indicate aspects of satisfaction and decision-making. 
# Additionally, words such as "pilots", "boeing", "technical", and "bulkhead" point to a focus on aircraft, flight experience, and technical aspects. There are also references to "hungry", "shoes", and "stress", which might indicate personal comfort, food service, and travel-related stress.
# 
# Topic 2: Service & Logistics
# This topic appears to revolve around customer service, logistics, and facilities. 
# Words like "courteous", "rate", "smile", and "impeccable" suggest a focus on service quality. 
# Meanwhile, "transfers", "terminals", "works", "website", and "terminal" indicate discussions about airport logistics, online services, and operational aspects. Additionally, words such as "stayed", "returned", and "copenhagen" might reflect travel experiences related to accommodation and destinations.
# 
# In summary:
#   Topic 1 seems to focus on passenger experiences, flight-related aspects, and travel comfort.
#   Topic 2 emphasizes service quality, airport operations, and logistical considerations.


### Figure 3
Pwts = PW_T; Pwts_n = t(apply(Pwts,1,function(x)x/sum(x)))
Wds = Frex_sorted_std
#Iid = matrix(as.numeric(sapply(as.vector(Wds),function(x)which(x==rownames(Pwts)))),ncol = NCOL(Wds))

iid1 <- c(which(rownames(Pwts_n)=="hungry"),which(rownames(Pwts_n)=="happy"),which(rownames(Pwts_n)=="book"),which(rownames(Pwts_n)=="complain"),which(rownames(Pwts_n)=="report"),which(rownames(Pwts_n)=="pilots"),which(rownames(Pwts_n)=="aircraft"),which(rownames(Pwts_n)=="technical"),which(rownames(Pwts_n)=="bulkhead"),which(rownames(Pwts_n)=="wheel"))
iid2 <- c(which(rownames(Pwts_n)=="courteous"),which(rownames(Pwts_n)=="smile"),which(rownames(Pwts_n)=="transfers"),which(rownames(Pwts_n)=="terminal"),which(rownames(Pwts_n)=="feet"),which(rownames(Pwts_n)=="website"),which(rownames(Pwts_n)=="claims"),which(rownames(Pwts_n)=="works"),which(rownames(Pwts_n)=="mention"),which(rownames(Pwts_n)=="inedible"))
Iid <- cbind(iid1,iid2)

cls <- c("#EEB422","#00B2EE")
X <- data.frame(mod2_std$PD_T,yrs,row.names = NULL); X$tPrev <- apply(X[,1:2],1,which.max)
A <- t(mapply(function(j){x<-table(X$tPrev[X$yrs==j]);x/sum(x)},unique(X$yrs)))

tikzDevice::tikz(file='../contribution/fig3.tex',width=10,height=4.5,sanitize = TRUE)
par(mar=c(5,8,5,2)+0.1,mfrow=c(1,3)) # Doubles left margin.
j=1;barplot(t(round(Pwts_n[Iid[1:10,j],],3)),border = FALSE,beside = FALSE,horiz = TRUE,las=2,col = c("#EEB422","#00B2EE","#B4EEB4"),cex.names = 1.75); title(paste0("Topic ",j), line = 1,adj=0,cex.main=1.55)
j=2;barplot(t(round(Pwts_n[Iid[1:10,j],],3)),border = FALSE,beside = FALSE,horiz = TRUE,las=2,col = c("#EEB422","#00B2EE","#B4EEB4"),cex.names = 1.75); title(paste0("Topic ",j), line = 1,adj=0,cex.main=1.55)
barplot(t(A),col = cls,cex.axis = 1.55,cex.names = 1.55,border=FALSE)
add_legend("bottom",legend = c("Topic 1","Topic 2"),fill = c("#EEB422","#00B2EE"),border = FALSE,bty = "n",cex=1.95,ncol=4)
dev.off()














