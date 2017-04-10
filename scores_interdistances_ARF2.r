rm(list=ls())
library(Biostrings)
nRegion <- 600
library(stringr)
source("PNfonctions.r")                 # fonctions auxiliaires



#-------------------------------------read matrices ------------------------------------------


pfm_ARF2<- read.table("m_ARF5.txt",header=TRUE,sep="\t",skip=1)
pfm_ARF2 <- round((t(as.matrix(pfm_ARF2)))*nRegion)+1 ;pfm_ARF2
maxi_ARF2 <- apply(pfm_ARF2,FUN=max, 2)
maxi_ARF2 <- matrix(nrow=4, rep(maxi_ARF2,4),byrow=TRUE)
pwm_ARF2 <- log(pfm_ARF2/maxi_ARF2)
pwm_ARF2_rev <- pwm_ARF2 - minScore(pwm_ARF2)/dim(pwm_ARF2)[2] 

pwm_ARF2 <-  reverseComplement(pwm_ARF2_rev) ; pwm_ARF2

#-------------------------------------read fasta files-----------------------------------------


ARF2_pos <- readDNAStringSet('MP_1000_pos.fas')#[(1:1000)]
ARF2_neg <- readDNAStringSet('MP_1000_1_neg.fas')#[(1:1000)]
width_pos <- width(ARF2_pos)
width_neg <- width(ARF2_neg)

seq_pos <- as.character(ARF2_pos)
seq_neg <- as.character(ARF2_neg)


#-------------------------------------Compute Scores-----------------------------------------

th <- maxScore(pwm_ARF2) - 9

#pos

scores_ARF2_pos<- mapply(seq_pos,FUN=PWMscoreStartingAt,SIMPLIFY=FALSE,  starting.at=mapply(seq,1,width_pos-dim(pwm_ARF2)[2],SIMPLIFY=FALSE),MoreArgs=list(pwm=pwm_ARF2))

scores_ARF2_rev_pos <- mapply(seq_pos,FUN=PWMscoreStartingAt,SIMPLIFY=FALSE,  starting.at=mapply(seq,1,width_pos-dim(pwm_ARF2_rev)[2],SIMPLIFY=FALSE),MoreArgs=list(pwm=pwm_ARF2_rev))

density_pos <- (sapply(FUN=sum,lapply(FUN=">",scores_ARF2_pos,th))+sapply(FUN=sum,lapply(FUN=">",scores_ARF2_rev_pos,th)))/(width_pos*2)

#neg

scores_ARF2_neg <- mapply(seq_neg,FUN=PWMscoreStartingAt,SIMPLIFY=FALSE,  starting.at=mapply(seq,1,width_neg-dim(pwm_ARF2)[2],SIMPLIFY=FALSE),MoreArgs=list(pwm=pwm_ARF2))

scores_ARF2_rev_neg <- mapply(seq_neg,FUN=PWMscoreStartingAt,SIMPLIFY=FALSE,  starting.at=mapply(seq,1,width_neg-dim(pwm_ARF2_rev)[2],SIMPLIFY=FALSE),MoreArgs=list(pwm=pwm_ARF2_rev))

density_neg <- (sapply(FUN=sum,lapply(FUN=">",scores_ARF2_neg,th))+sapply(FUN=sum,lapply(FUN=">",scores_ARF2_rev_neg,th)))/(width_neg*2)

#-------------------------------------Compute ROC-----------------------------------------

rc1 = ROCcurve(pos,neg) # fait la roc
X <- rc1$XY[1,]
Y <- rc1$XY[2,]
rc2 = ROCcurve(pos_pen,neg_pen) # fait la roc
X_pen <- rc2$XY[1,]
Y_pen <- rc2$XY[2,]
rc3 = ROCcurve(pos_pen2,neg_pen2) # fait la roc
X_pen2 <- rc3$XY[1,]
Y_pen2 <- rc3$XY[2,]
rc4 = ROCcurve(pos_pen3,neg_pen3) # fait la roc
X_pen3 <- rc4$XY[1,]
Y_pen3 <- rc4$XY[2,]
#
AU <- rc1$AUC
A <- as.character(round(AU,4))
AUC <- paste("AUC = ", A,sep="")
#
AU_pen<- rc2$AUC
A_pen<- as.character(round(AU_pen,4))
AUC_pen<- paste("AUC with penalties = ", A_pen,sep="")
#
AU_pen2<- rc3$AUC
A_pen2<- as.character(round(AU_pen2,4))
AUC_pen2<- paste("AUC with penalties and density = ", A_pen2,sep="")
#
AU_pen3<- rc4$AUC
A_pen3<- as.character(round(AU_pen3,4))
AUC_pen3<- paste("AUC with penalties, density and corelation = ", A_pen3,sep="")
#
{plot(X,Y,type="l",col="red",lwd=2,
      ylab="ARF2",xlab="ARF5",
      main="ARF2 vs ARF5")}
lines(X_pen,Y_pen,col='cornflowerblue',lwd=2)
lines(X_pen2,Y_pen2,col='green4',lwd=2)
#lines(X_pen3,Y_pen3,col='orange',lwd=2)
{legend(0.3,0.15,legend=c(AUC,AUC_pen,AUC_pen2),#AUC_pen3),
        col=c("red","cornflowerblue","green4","orange"),lty=rep(1,4),lwd=rep(2,4))}

dev.copy(device = png, filename = 'ROC_ARF2_vs_ARF5_test.png', width = 800, height = 600) 
dev.off()

