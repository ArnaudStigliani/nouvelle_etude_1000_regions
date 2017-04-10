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

############################################### pos ############################################

#-------------------------------------read fasta files-----------------------------------------


ARF2_pos <- readDNAStringSet('MP_1000_pos.fas')#[(1:1000)]
width_pos <- width(ARF2_pos)
seq_pos <- as.character(ARF2_pos)



#-------------------------------------Compute Scores-----------------------------------------


scores_ARF2_pos<- sapply(seq_pos,FUN=PWMscoreStartingAt,starting.at=(1:(width_pos[1]-dim(pwm_ARF2)[2]+1)),pwm=pwm_ARF2)

scores_ARF2_rev_pos <-  sapply(seq_pos,FUN=PWMscoreStartingAt,starting.at=(1:(width_pos[1]-dim(pwm_ARF2)[2]+1)),pwm=pwm_ARF2_rev)

pos <- apply(FUN=max,ifelse(scores_ARF2_pos > scores_ARF2_rev_pos, scores_ARF2_pos,scores_ARF2_rev_pos),2)

############################################### neg ############################################

list_neg <- list('MP_1000_1_neg.fas','MP_1000_2_neg.fas','MP_1000_3_neg.fas','MP_1000_4_neg.fas')

i <- 0
rc <- list()
X <- list()
Y <- list()
AU <- list()
A <- list()
AUC <- list()
for (elt in list_neg)
{
    i <- i+1
    ARF2_neg <- readDNAStringSet(elt)
    width_neg <- width(ARF2_neg)
    seq_neg <- as.character(ARF2_neg)
    #
    scores_ARF2_neg<- sapply(seq_neg,FUN=PWMscoreStartingAt,starting.at=(1:(width_neg[1]-dim(pwm_ARF2)[2]+1)),pwm=pwm_ARF2)
    #
    scores_ARF2_rev_neg <-  sapply(seq_neg,FUN=PWMscoreStartingAt,starting.at=(1:(width_neg[1]-dim(pwm_ARF2)[2]+1)),pwm=pwm_ARF2_rev)
    #
    neg <- apply(FUN=max,ifelse(scores_ARF2_neg > scores_ARF2_rev_neg, scores_ARF2_neg,scores_ARF2_rev_neg),2)
    #
    rc[[i]] = ROCcurve(pos,neg) # fait la roc
    X[[i]] <- rc[[i]]$XY[1,]
    Y[[i]] <- rc[[i]]$XY[2,]
    AU[[i]] <- rc[[i]]$AUC
    A[[i]] <- as.character(round(AU[[i]],4))
    AUC[[i]] <- paste("AUC = ", A[[i]],sep="")

}
#-------------------------------------Compute ROC-----------------------------------------


{plot(X[[i]],Y[[i]],type="l",col="red",lwd=2,
      ylab="ARF2",xlab="ARF5",
      main="ARF2 vs ARF5")}
lines(X[[2]],Y[[2]],col='green4',lwd=2)
lines(X[[3]],Y[[3]],col='orange',lwd=2)
lines(X[[4]],Y[[4]],col='cornflowerblue',lwd=2)
{legend(0.3,0.15,legend=c(AUC[[1]],AUC[[2]],AUC[[3]],AUC[[4]]),
        col=c("red","green4","orange","cornflowerblue"),lty=rep(1,4),lwd=rep(2,4))}

## dev.copy(device = png, filename = 'ROC_ARF2_vs_ARF5_test.png', width = 800, height = 600) 
## dev.off()

