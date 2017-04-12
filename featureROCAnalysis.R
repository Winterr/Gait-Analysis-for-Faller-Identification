## -------------------------------------------------------
# R functions to run ROC analysis on the  
#	feature extracted from gait data
#
# Moacir A. Ponti / 2016
# in collaboration with: Patricia Bet; Paula Castro
# -------------------------------------------------------


# Read CSV File
dataRead <- function(filename) {
	data <- read.csv(filename, header=FALSE, sep=",")
}

# Normalize variable to 0-1
normalizeZeroOne <- function(datav) {

	a <- (datav-min(datav))/(max(datav)-min(datav))	

	return(a/max(a))
}

#' Generates values for plotting a ROC curve which computes 
#'	(1) sensitivity: the test's ability to correctly 
#'	  detect patients who DO have the condition. 
#'	  Is also known as TPR, true positive rate.
#'	(2) specificity: the test's ability to correctly
#'	  detect patients WITHOUT a condition
#'	(3) false positive rate: which is 1-specificity
#'
ROC <- function(dataraw, label, nsteps=10, invert=1) {
	
	n <- length(dataraw)
	posMask <- (label == 1)

	if (invert == 1) {
		data <- invert-normalizeZeroOne(dataraw)
	} else {
		data <- normalizeZeroOne(dataraw)
	}
		
	step = 1/nsteps

	# computes sensitivity (or TPR), specificity (or TNR) 
	# for several thresholds
	TPR <- matrix(vector(),1,nsteps)
	TNR <- matrix(vector(),1,nsteps)
	ACC <- matrix(vector(),1,nsteps)
	CUT <- matrix(vector(),1,nsteps)

	print(nsteps)
	dshow <- cbind(round(data,2),label)
	dshow <- cbind(dshow[order(dshow[,1]),1], dshow[order(dshow[,1]),2])
	print(dshow) 

	i <- 1
	inits <- step
	ends <- 1-step
	TPR[i] <- 1.0
	TNR[i] <- 0.0
	ACC[i] <- 0.5
	CUT[i] <- 0
	i = i + 1
	for (t in seq(inits,ends,by=step)) {
		predLab <- as.integer(data >= t)
		cat(paste(i,t),"\n")

		print(cbind(label, predLab, data))
		
		TP <- sum(predLab[label==1] == label[label==1])
		TN <- sum(predLab[label==0] == label[label==0])

		FP <- sum(predLab[label==0])
		FN <- sum(!predLab[label==1])

		cat(paste(t)," Conf Matrix:\n")
		cat("pos:", paste(TP,FN),"\n")
		cat("neg:", paste(FP,TN),"\n\n")

		TPR[i] <- TP/(TP+FN) #sensitivity
		TNR[i] <- TN/(TN+FP) #specificity
		ACC[i] <- 2*TP/(2*TP+FP+FN)
		CUT[i] <- t
		
		cat("TPR, TNR, fS:", paste(TPR[i],1-TNR[i],ACC[i]),"\n\n")
		i = i + 1
	}
	TPR[i] <- 0
	TNR[i] <- 1
	ACC[i] <- 0.5
	CUT[i] <- 1
	#print(TPR)
	#print(1-TNR)

	# for ROC we usually output False Positive Rate
	# which is 1 - specificity:
	FPR = 1-TNR

	FPR <- rev(FPR)
	TPR <- rev(TPR)
	ACC <- rev(ACC)
	CUT <- CUT
	# AUC (for debug, remove)
	#he <- (TPR[-1] + TPR[-length(TPR)])/2
	#wi <- diff(FPR)
	#print(sum(he*wi))

	ROC <- cbind(FPR,TPR, ACC, CUT)
}

#' Compute Area Under the Curve ROC (AUC)
#'
computeAUC <- function(ROCvalues) {
	TPR <- ROCvalues[,1]
	FPR <- ROCvalues[,2]
	he <- (TPR[-1] + TPR[-length(TPR)])/2
	wi <- diff(FPR)
	return(1-sum(he*wi))
}

computeROCs <- function(data, nsteps=4, type=1, invert=1) {

	if (type==1) {
		variables <- c('d(PSEs,PSEc)', 'd(PSPs,PSPc)','d(PSPFt,PSPFm)','d(WPSPm,WPSPc)','WPSP_c,3','PSE_c','WPSP_c,2')
		colours <- seq(2,8)
	} else if (type == 2) {
		variables <- c('d(PSEs,PSEc)', 'd(PSPs,PSPc)','d(PSPFt,PSPFm)','d(WPSPm,WPSPc)')
		colours <- seq(2,5)
	} else if (type == 3) {
		variables <- c('WPSP_c,3','PSE_c','WPSP_c,2')
		colours <- seq(2,4)
	} else if (type == 4) {
		variables <- c('TUG','TUG-M','TUG-C')
		colours <- seq(2,4)
	}
	n <- nrow(data)
	m <- ncol(data)

	labels <- data[,m]

	print(paste("Compute ROC analysis for: ", nsteps))
	for (i in seq(1,m-1)) {
		print(i)
		rocv <- ROC(data[,i],labels,nsteps=nsteps, invert=invert)

		mAUC <- computeAUC(rocv)
		title <- paste0(variables[i],', AUC=', round(mAUC,4))
		variables[i] <- title

		print(title)
		if (i == 1) {
			plot(rocv[,1:2], lwd=3, pch=i, col=i+1, xlim=c(0,1), ylim=c(0,1), ylab="Sensitivity (TPR)", xlab="[1 - Specificity] (FPR)", cex.lab=1.25, cex.axis=1.25)
			axis(1,at=seq(0,1,by=0.1),cex.axis=1.25)
			axis(2,at=seq(0,1,by=0.1),cex.axis=1.25)

			lines(rocv[,1:2], lwd=3, col=i+1)
			segments(0,0, 1,1, col=1, lwd=1)
			segments(0,1, 1,0, col='darkgray', lwd=1, lty=2)
		} 
		else {
			points(rocv[,1:2], lwd=3, pch=i, col=i+1, xlim=c(0,1), ylim=c(0,1), ann=FALSE)
			lines(rocv[,1:2], lwd=3, col=i+1)
		
		}
	}

	legend("bottomright",
	       inset=.05, cex = 1.2,
       		variables, col=colours,	pch = colours-1,
	       horiz=FALSE, lty=c(1,1), lwd=c(2,2), bg="grey96")
}
