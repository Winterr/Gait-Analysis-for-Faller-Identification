  variables = c('avg[TUG sec]','avg[Feats]','avg[Dists]')
  ROCSF <- ROC(fusionSumTUGSec, TUGs[,4], 200, invert=0)
  plot(ROCSF[,4],ROCSF[,2],t="l",lty=1,lwd=3, xlab="Probability threshold", ylab="Sensitivity X Specificity", cex.lab=1.25, cex.axis=1.2)
  lines(ROCSF[,4],1-ROCSF[,1],t="l",lty=2,lwd=3)
  axis(1,at=seq(0,1,by=0.1),cex.axis=1.2)
  axis(2,at=seq(0,1,by=0.1),cex.axis=1.2)
  segments(.598, 0, .598, 1,lty=3,lwd=3, col="#000077")
  legend(0.01,0.88,inset=.05, lwd=c(2,2),lty=c(1,2),c("Sensitivity","Specificity"),cex=1.25,bg="white",col="#333333")
  dev.print(pdf,"fig_SS_avgTUG.pdf")

  mAUC <- computeAUC(ROCSF)
  title <- paste0(variables[1],', AUC=', round(mAUC,3))
  variables[1] <- title

  ROCFF1 <- ROC(fusionTUGCsAvg[,1], fTUG[,8], 200,invert=1)
  plot(ROCFF1[,4],ROCFF1[,2],t="l",lty=1,lwd=3, xlab="Probability threshold", ylab="Sensitivity X Specificity", cex.lab=1.25, cex.axis=1.2)
  lines(ROCFF1[,4],1-ROCFF1[,1],t="l",lty=2,lwd=3)
  axis(1,at=seq(0,1,by=0.1),cex.axis=1.2)
  axis(2,at=seq(0,1,by=0.1),cex.axis=1.2)
  segments(.513, 0, .513, 1,lty=3,lwd=3, col="#000077")
  legend(0.01,0.88,inset=.05, lwd=c(2,2),lty=c(1,2),c("Sensitivity","Specificity"),cex=1.25,bg="white",col="#333333")
  dev.print(pdf,"fig_SS_avgFeats.pdf")

  mAUC <- computeAUC(ROCFF1)
  title <- paste0(variables[2],', AUC=', round(mAUC,3))
  variables[2] <- title

  ROCFF2 <- ROC(fusionMeanFeats2TUG[,1], fTUG[,8], 200,invert=1)
  plot(ROCFF2[,4],ROCFF2[,2],t="l",lty=1,lwd=3, xlab="Probability threshold", ylab="Sensitivity X Specificity", cex.lab=1.25, cex.axis=1.2)
  lines(ROCFF2[,4],1-ROCFF2[,1],t="l",lty=2,lwd=3)
  axis(1,at=seq(0,1,by=0.1),cex.axis=1.2)
  axis(2,at=seq(0,1,by=0.1),cex.axis=1.2)
  segments(.5, 0, .5, 1,lty=3,lwd=3, col="#000077")
  legend(0.01,0.88,inset=.05, lwd=c(2,2),lty=c(1,2),c("Sensitivity","Specificity"),cex=1.25,bg="white",col="#333333")
  dev.print(pdf,"fig_SS_avgDists.pdf")

  mAUC <- computeAUC(ROCFF2)
  title <- paste0(variables[3],', AUC=', round(mAUC,3))
  variables[3] <- title
   
  plot(ROCSF[,1:2], lwd=3, pch=1, col=2, xlim=c(0,1), ylim=c(0,1), ylab="Sensitivity (TPR)", xlab="[1 - Specificity] (FPR)", cex.lab=1.25, cex.axis=1.25)
  segments(0,0, 1,1, col=1, lwd=1)
  segments(0,1, 1,0, col='darkgray', lwd=1, lty=2)
  lines(ROCSF[,1:2], lwd=3, col=2)
  axis(1,at=seq(0,1,by=0.1),cex.axis=1.25)
  axis(2,at=seq(0,1,by=0.1),cex.axis=1.25)

  points(ROCFF1[,1:2], lwd=3, pch=2, col=3, xlim=c(0,1), ylim=c(0,1), ann=FALSE)
  lines(ROCFF1[,1:2], lwd=3, col=3)

  points(ROCFF2[,1:2], lwd=3, pch=3, col=4, xlim=c(0,1), ylim=c(0,1), ann=FALSE)
  lines(ROCFF2[,1:2], lwd=3, col=4)

  legend("bottomright", inset=.05, cex = 1.2, variables, col=c(2,3,4), pch = c(1,2,3),
               horiz=FALSE, lty=c(1,1,1), lwd=c(2,2), bg="grey96")
  dev.print(pdf,"fig_ROC_fusion.pdf") 
