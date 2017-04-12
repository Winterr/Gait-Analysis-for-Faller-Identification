# -------------------------------------------------------
# R functions to manipulate gait data:
# 	- feature extraction
#	- signal segmentation
#
# Moacir A. Ponti / 2016
# in collaboration with: Patricia Bet; Paula Castro
#
# Code prepared for purpose of paper reviewing
# -------------------------------------------------------

library(signal)

# Function to remove null entries on the beginning of the list
rmNullEntries <- function(x) {
	k <- 0
	while ( is.null(x[[1]]) ) { 
		k <- k + 1
		x[[1]] <- NULL 
	}
	print(paste("total null entries removed", k))
	return (x)
}


# Function that gets a table with 3 columns (axis x,y,z) and returns
# a single column
axisFusion <- function(inputXYZ) {

	# performs the fusion of the triaxial data into one time series
	# using the euclidian sum: square root of the squared sum
	newVector <- sqrt((inputXYZ[,1]^2) + (inputXYZ[,2]^2) + (inputXYZ[,3]^2))

	# remove small values, keeping only observations higher then 0.01
	newVector <- newVector[newVector > 0.01]

} 

# Function to plot the acceleration curve
plotGait <- function(input) {

	# plot an empty window (type="n") and then display gait as a line
	plot(input, type="n")
	lines(input)
}

# Function to read all files and store the data into a list
# it works by getting the list of files to be read as input 
# and returns the list with the data ready to be analysed
# Example:
#    dataList <- readFiles("filelist.txt")
#
#
# 'dataList' is a list in which each data can be accessed by
# the name of the file it contained in the format:
# dataList[['filename.dat']]
#
# Examples:
#    dataList[['C_009.dat']]
#    plotGait(datalist[['C_000_1.dat']]
#
readFilesOld <- function(filelist) {

	print('Deprecated function')

	# read the file containg the name of each datafile to open
	# OBS: uses t() to transpose it into a vector
	filenames <- t( read.table(filelist,header=F) )

	# create a list to store the data, it is a vector with
	# length equal to the number of filenames
	listdata <- vector(mode="list", length=length(filenames))

	# for each file 'f' in the filenames read
	for (f in 1:length(filenames)) {
		# read the current file into 'x'
		x <- read.table(filenames[f])
		# perform fusion
		xfusion <- axisFusion(x)
	
		# include the data into the list
		listdata[[filenames[f]]] <- xfusion
	}
	return(listdata)
}

## Function to read all files and store the data into a list
# It works by getting the list of files to be read as input 
# and returns the list with the data ready to be analysed
# Example:
#    dataList <- readFiles("filelist.txt")
#
#
# 'dataList' is a list in which each data can be accessed by
# the name of the file it contained in the format:
# dataList[['filename.dat']]
#
# Examples:
#    dataList[['C_009.dat']]
#    plotGait(datalist[['C_036.dat']]
#
# OBS: important, it MERGES consecutive files
#      for example, 'C_009.dat' 'C_009_1.dat' and 'C_009_2.dat'
#      become a single entry named 'C_009.dat'
readFiles <- function(filelist) {

	# read the file containg the name of each datafile to open
	# OBS: uses t() to transpose it into a vector
	filenames <- t( read.table(filelist,header=F) )

	# create a list to store the data, it is a vector with
	# length equal to the number of filenames
	listdata <- vector(mode="list", length=length(filenames))

	num <- 0

	# for each file 'f' in the filenames read
	for (f in 1:length(filenames)) {
		# read the current file into 'x'
		x <- read.table(filenames[f])

		print(filenames[f])
		# performs fusion
		xfusion <- axisFusion(x)
		
		## extract the number of the data
		# check if it is a continuation
		snum <- sub(".*_(\\d{1}).*","\\1",filenames[f])

		# extract number with 3 digits
		c_num <- sub("\\D*(\\d{3}).*","\\1",filenames[f])

		if (snum > 0 && c_num == num) {
			print(paste0("...continuation: ",num, " and ", as.numeric(snum)))
	
			f <- f - as.numeric(snum)
			xfusion <- c(listdata[[filenames[f]]], xfusion)

			strfind <- paste0("_",snum)
			filenames[f] <- gsub(strfind, "", filenames[f])
			print(paste0("concatenating to ",filenames[f]))
		} else if (snum > 0) {
			strfind <- paste0("_",snum)
			filenames[f] <- gsub(strfind, "", filenames[f])
		}
		# include the data into the list
		listdata[[filenames[f]]] <- xfusion
		num <- c_num
	}

	listdata <- rmNullEntries(listdata)
	return(listdata)
}

# Prints the length of all data in the list
printLengths <- function(gaitList) {

	for (j in names(gaitList)) {
		print(paste(j, length(gaitList[[j]])))
	}

}

# Compute the Normalized Power Spectrum of the gait data
# it allows to study the frequencies within the signal
# 
# Parameters:
#       norm (default FALSE)
#		TRUE: uses the parameter 'norm.value' in order to normalize
#		FALSE: performs regular normalization (ignoring 'norm.value')
#	norm.value: used only if norm is TRUE
#
# Examples: 
#	fftGait1 <- fourierGait(gait1, norm=TRUE, 1024)
#       fftGait2 <- fourierGait(gait2)
fourierGait <- function(g, power2n=TRUE, norm=FALSE, norm.value) { 
	# compute power of two
	if (power2n == TRUE) {
		g.n <- 2^floor(log2( length(g)) ) 
		#print(paste("n = ", length(g), " n_fft = ", g.n))
	} else {
		g.n <- length(g)
	}

	# calculate fft using g.n samples
	g.F <- fft(g[1:g.n])

	# getting half transform (due to symmetry) and normalizing
	if (norm == TRUE) {
		g.F <- g.F[1:(g.n/2)] / (norm.value)
	} else {
		g.F <- g.F[1:(g.n/2)] / (g.n/2)
	}

	return(abs(g.F))
}

# Compute the frequency analysis for all gait data in a list
# with the option of filtering using butterworth low pass filter
#
# Parameters
#	filter (default TRUE): enable filtering
#	cutoff : cutoff frequencies (default 1/100, i.e. 100 Hz)
#
# Examples:
#	gd2.F <- compute.Fourier.List(gd, filter=TRUE, cutoff=1/50)
#
compute.Fourier.List <- function(gaitList, filter=TRUE, cutoff=(1/100), pow2n=TRUE) {
	gaitList.F <- vector(mode="list")

	require(signal)
	bf <- butter(1, cutoff)

	for (g in names(gaitList)) {
		if (filter) {
			gfilt <- filtfilt(bf, gaitList[[g]])
		} else {
			gfilt <- gaitList[[g]]
		}
		#gaitList.F[[g]] <- fourierGait(gfilt, norm=TRUE, norm.value=32768)
		gaitList.F[[g]] <- fourierGait(gfilt, power2n=pow2n)
	}
	return(gaitList.F)
}

# compute PSE (power spectral entropy) feature
feat.pse <- function(S, cutoff=0, first=1) {
	n <- length(S)
	if (cutoff == 0) cutoff = n
	
	return( -sum(S[first:cutoff] * log2( S[first:cutoff] + 0.001 )) )
}

# compute PSP (power spectrum peak)
feat.psp <- function(S, cutoff=0){
	n <- length(S)
	if (cutoff == 0) cutoff = n

	max1 <- max(S[4:cutoff])
	pmax1 <- which(S == max1)
	
	max2 <- max(S[(pmax1+10):cutoff])
	pmax2 <- which(S == max2)
	
	max3 <- max(S[(pmax2+10):cutoff])
	pmax3 <- which(S == max3)
	
	return (c(max1, max2, max3))
}

# compute PSPF (power spectrum peak frequency)
feat.pspf <- function(S, cutoff=0) {
	n <- length(S)
	if (cutoff == 0) cutoff = n

	max1 <- max(S[4:cutoff])
	pmax1 <- which(S == max1)
	
	max2 <- max(S[(pmax1+10):cutoff])
	pmax2 <- which(S == max2)
	
	max3 <- max(S[(pmax2+10):cutoff])
	pmax3 <- which(S == max3)
	
	return (c(pmax1, pmax2, pmax3))
}


# compute WPSP (weighed power spectrum peak)
feat.wpsp <- function(S, cutoff=0) {
	n <- length(S)
	if (cutoff == 0) cutoff = n

	max1 <- max(S[4:cutoff])
	pmax1 <- which(S == max1)
	
	max2 <- max(S[(pmax1+10):cutoff])
	pmax2 <- which(S == max2)
	
	max3 <- max(S[(pmax2+10):cutoff])
	pmax3 <- which(S == max3)
	
	return (c(pmax1*max1, pmax2*max2, pmax3*max3))
}

#' computes Number of Steps using zerocrossing information
#'
#' Number of Steps is a feature commonly used to characterize gait
#' @param s The signal from which to extract the feture
#' @param labs labels to segment TUGs, will consider only values for which labs>0
#' @return number of steps
#' @export
feat.nsteps <- function(s, labs=NULL) {

	s <- (s - mean(s))^2
	bf <- butter(2, 1/5, type="low")
	s <- filtfilt(bf, s)

	factor <- 2

	#plot(s,t="l")
	if (!(is.null(labs))) {
		s[which(labs==0)] <- 0
		factor <- 6
	}
	#lines(s,col="green")
	#line <- readline()
	
	#maxS <- which(S == max(S))
	#thresh <- mean(S[-maxS])*7
	thresh <- mean(s)*factor

	s[s < thresh] <- 0
	#print(thresh)
	return (round(sum(diff(sign(s)) != 0)/4))
}


# Tries to get features from a list of gaitFourier data
#
compute.frequency.features <- function(fourierList) {

	gaitFeat <- vector(mode="list")

	for (g in names(fourierList)) {

		smooth <- rollapply(fourierList[[g]][1:130], width=5, FUN=mean, align="center", fill=NA)
		smooth <- (smooth / max(smooth[5:125])) * max(fourierList[[g]][1:130])

		peaks <- sort(smooth, decreasing=TRUE)[1:20]

		# get first harmonic
		f.g <- peaks[1]; featA.g <- f.g;
		# get first harmonic index
		featF.g <- which(fourierList[[g]] == f.g)
		print(f.g)
		# go through frequencies
		for (i in (featF.g+10):length(fourierList[[g]])) {
			print(paste(i, " = ", fourierList[[g]][i]))
			diff <- abs(fourierList[[g]][i]) - f.g
			print(paste("diff = ", diff))
		
#			while (  <= 10 ) {
#				i <- i + 1
#				diff <- abs(abs(fourierList[[g]][i]) - f.g)
#				print(paste("diff = ", diff))
#			}
			f.g <- fourierList[[g]][i]
			featA.g <- cbind(featA.g, f.g)
			featF.g <- cbind(featF.g, i)
		}
		print(featA.g)
		print(featF.g)	
	}
}

groupsGait <- function() {
	grps <- c("CR", "NC", "CR", "NC", "NC", "CE", "CE", "CE", "NC", "CE", 
		    "NC", "CE", "NC", "CE", "CR", "CR", "CR", "CR", "NC", "NC",
		    "CR", "CE", "CE", "CE", "CR", "CE", "NC", "CR", "NC", "CR")
	return(grps)
}

	# it considers only the data without errors (note missing numbers)
groupsTUG <- function() {
		#  1     2    3      4    5     6     8      11    13    14
		#  1     2    3      4    5     6     7      8      9   10
	grps <- c("NC", "C" , "NC", "NC", "C" , "NC", "NC", "C" , "NC", "NC", 
		# 15    16    17    18    19    20     21     22    23    24 
		# 11    12    13    14    15    16     17     18    19    20 
	          "NC", "C" , "C" , "C" , "NC", "NC", "NC", "NC", "NC", "C" ,
		# 25    26    27    28    29    30     31    33    34    35  
		# 21    22    23    24    25    26     27    28    29    30  
	          "C" , "C" , "C" , "NC", "C" , "NC", "C" , "C" , "C" , "C",
		# 36,   37    38   40    41    42
		# 31,   32    33   34    35    36
		  "C", "NC", "NC", "NC", "C", "C")
	
	return(grps)
}



# Compute Frequency Features based on peaks
feature.frequency.peaks <- function(g.F, cutoff=0, groups, first=2) {

	features <- matrix(vector(), length(g.F), 10)

	j <- 1
	for (g in g.F) {
		features[j, 1] <- feat.pse(g, cutoff=150, first)
		features[j, 2:4] <- feat.psp(g, cutoff=cutoff)
		features[j, 5:7] <- feat.pspf(g, cutoff=cutoff)
		features[j, 8:10] <- feat.wpsp(g, cutoff=cutoff)
		j <- j + 1
	}

	feats <- data.frame(features)
	feats["id"] <- names(g.F)
	feats["group"] <- groups
	return(feats)
}


feature.frequency.energies <- function(gd2.F, lfr = c(2, 5, 25, 45, 65) , ufr = c(100,35,55,75,95), groups) {

	features <- matrix(vector(), length(gd2.F), length(lfr))

	j <- 1
	for (g in gd2.F) {
		for (i in 1:length(lfr)) {
			if (!exists("freq.f")) {
				freq.f <- sum(g[lfr[i] : ufr[i]])
			} else {
				freq.f <- c(freq.f, sum(g[lfr[i] : ufr[i]]))
			}
		}
		features[j,] <- freq.f
		j <- j + 1
		rm(freq.f)
	}

	feats <- data.frame(features)
	feats["id"] <- names(gd2.F)
	feats["group"] <- groups

	return(feats)
}

plotEnergiesTUG <- function(gd2.F, startF = 2, endF = 100, step=25, stride=0) {

	pltNC_C <- c()
	pltLab <- c()
	pltTic <- c()
	pltPVa <- c()
	pltMax <- c()
	for (lfr in seq(startF, endF, by=step/2)) {
		print(lfr)
		ufr <- lfr+step
		eNC <- c( 
			sum(gd2.F[['TUG_001.dat']][lfr:ufr]),
			sum(gd2.F[['TUG_003.dat']][lfr:ufr]),
			sum(gd2.F[['TUG_004.dat']][lfr:ufr]),
			sum(gd2.F[['TUG_006.dat']][lfr:ufr]),
			sum(gd2.F[['TUG_008.dat']][lfr:ufr]),
			sum(gd2.F[['TUG_013.dat']][lfr:ufr]),
			sum(gd2.F[['TUG_014.dat']][lfr:ufr]),
			sum(gd2.F[['TUG_015.dat']][lfr:ufr]),
			sum(gd2.F[['TUG_019.dat']][lfr:ufr]),
			sum(gd2.F[['TUG_020.dat']][lfr:ufr]),
			sum(gd2.F[['TUG_021.dat']][lfr:ufr]),
			sum(gd2.F[['TUG_022.dat']][lfr:ufr]),
			sum(gd2.F[['TUG_023.dat']][lfr:ufr]),
			sum(gd2.F[['TUG_028.dat']][lfr:ufr]),
			sum(gd2.F[['TUG_030.dat']][lfr:ufr]),
			sum(gd2.F[['TUG_037.dat']][lfr:ufr]),
			sum(gd2.F[['TUG_038.dat']][lfr:ufr]),
			sum(gd2.F[['TUG_040.dat']][lfr:ufr]) )

		eC <- c( 
			sum(gd2.F[['TUG_002.dat']][lfr:ufr]),
			sum(gd2.F[['TUG_005.dat']][lfr:ufr]),
			sum(gd2.F[['TUG_011.dat']][lfr:ufr]),
			sum(gd2.F[['TUG_016.dat']][lfr:ufr]),
			sum(gd2.F[['TUG_017.dat']][lfr:ufr]),
			sum(gd2.F[['TUG_018.dat']][lfr:ufr]),
			sum(gd2.F[['TUG_024.dat']][lfr:ufr]),
			sum(gd2.F[['TUG_025.dat']][lfr:ufr]),
			sum(gd2.F[['TUG_026.dat']][lfr:ufr]),
			sum(gd2.F[['TUG_027.dat']][lfr:ufr]),
			sum(gd2.F[['TUG_029.dat']][lfr:ufr]),
			sum(gd2.F[['TUG_031.dat']][lfr:ufr]),
			sum(gd2.F[['TUG_033.dat']][lfr:ufr]),
			sum(gd2.F[['TUG_034.dat']][lfr:ufr]),
			sum(gd2.F[['TUG_035.dat']][lfr:ufr]),
			sum(gd2.F[['TUG_036.dat']][lfr:ufr]),
			sum(gd2.F[['TUG_041.dat']][lfr:ufr]),
			sum(gd2.F[['TUG_042.dat']][lfr:ufr]) )

		a <- t.test(eNC,eC)
		print(a$p.value)
		pltPVa <- rbind(pltPVa, a$p.value)
		pltMax <- rbind(pltMax, max(cbind(eNC,eC)))
		pltNC_C <- rbind(pltNC_C,eNC)
		pltNC_C <- rbind(pltNC_C,eC)
		pltLab <- cbind(pltLab,paste0(round(lfr-startF+1),'-',round(ufr-startF)))
		pltTic <- cbind(pltTic,paste0('NF'))
		pltTic <- cbind(pltTic,paste0('F'))
	}

	ct <- nrow(pltNC_C)
	print(ct)
	maxv = max(pltNC_C)+0.008
	minv = min(pltNC_C)-0.008
	boxplot(t(pltNC_C), names = pltTic, cex.lab=0.8, xlim=c(0.75,ct+0.25), ylim=c(minv, maxv+0.1))
	
	j <- 1
	print(pltPVa)
	for (i in seq(2.5,ct-1,by=2)) {
		segments(i,-0.15,i,5)
		print(pltLab[j])
		print(i)
		text(i-1 , maxv+0.09, pltLab[j], cex=1.3, col=4)

		if (pltPVa[j] <= 0.05) {
			segments(i-1.5, pltMax[j]+0.03, i-0.5, pltMax[j]+0.03)
			segments(i-1.5, pltMax[j]+0.03, i-1.5, pltMax[j]+0.02)
			segments(i-0.5, pltMax[j]+0.03, i-0.5, pltMax[j]+0.02)
			text(i-1, pltMax[j]+0.07, "*", cex=1.5)
		}
		j <- j + 1
	}
	text(i+1.0 , maxv+0.09, pltLab[j], cex=1.3, col=4)
}

segmentByLabel <- function(data, labels, expression) {


	newdata <- data
	for (i in 1:length(data)) {
		expr <- paste0("labels[[i]]",expression)
		newdata[[i]][which(!eval(parse(text=expr)))] <- 0
		newdata[[i]] <- newdata[[i]][which(newdata[[i]]>0)]
	}

	return(newdata)
}
	


plotEnergies <- function(gd2.F, lfr = 5, ufr = 100) {

	eN2 = c( 
		sum(gd2.F[[2]][lfr:ufr]),
		sum(gd2.F[[4]][lfr:ufr]),
		sum(gd2.F[[5]][lfr:ufr]),
		sum(gd2.F[[9]][lfr:ufr]),
		sum(gd2.F[[11]][lfr:ufr]),
		sum(gd2.F[[13]][lfr:ufr]),
		sum(gd2.F[[19]][lfr:ufr]),
		sum(gd2.F[[20]][lfr:ufr]),
		sum(gd2.F[[27]][lfr:ufr]),
		sum(gd2.F[[29]][lfr:ufr]) )

	eE2 = c(	sum(gd2.F[[6]][lfr:ufr]),
		sum(gd2.F[[7]][lfr:ufr]),
		sum(gd2.F[[8]][lfr:ufr]),
		sum(gd2.F[[10]][lfr:ufr]),
		sum(gd2.F[[12]][lfr:ufr]),
		sum(gd2.F[[14]][lfr:ufr]),
		sum(gd2.F[[22]][lfr:ufr]),
		sum(gd2.F[[23]][lfr:ufr]),
		sum(gd2.F[[24]][lfr:ufr]),
		sum(gd2.F[[26]][lfr:ufr]) )

	eR2 = c(	sum(gd2.F[[1]][lfr:ufr]),
		sum(gd2.F[[3]][lfr:ufr]),
		#sum(gd2.F[[15]][lfr:ufr]),
		sum(gd2.F[[16]][lfr:ufr]),
		sum(gd2.F[[17]][lfr:ufr]),
		sum(gd2.F[[18]][lfr:ufr]),
		sum(gd2.F[[21]][lfr:ufr]),
		sum(gd2.F[[25]][lfr:ufr]),
		sum(gd2.F[[28]][lfr:ufr]),
		sum(gd2.F[[30]][lfr:ufr]) )

	boxplot(eN2, eE2, eR2)
}

boxplotFeats <- function(data, featL=c(), labId=-1, featN=c()) {

	if (length(featL) == 0) {
		featL = seq(1, ncol(data)-1)
	}
	if (labId == -1) {
		labId = ncol(data)-1
	}
	if (length(featN) == 0) {
		featN = seq(1,length(featL))
	}

	pltNC_C <- c()
	pltTic <- c()
	pltPVa <- c()
	pltMax <- c()
	for (f in featL) {
	
		x <- normalizeZeroOne(data[,f])
		x1 <- x[data[,labId]=='C']
		x2 <- x[data[,labId]=='NC']
		a <- wilcox.test( x1, x2 )
		print(a$p.value)

		pltPVa <- rbind(pltPVa, a$p.value)
		pltMax <- rbind(pltMax, max(cbind(x1,x2)))

		pltNC_C <- rbind(pltNC_C,x1)
		pltNC_C <- rbind(pltNC_C,x2)

		pltTic <- cbind(pltTic,paste0('F'))
		pltTic <- cbind(pltTic,paste0('NF'))
	}

	ct <- nrow(pltNC_C)
	print(ct)
	maxv = max(pltNC_C)+0.008
	minv = min(pltNC_C)-0.008
	boxplot(t(pltNC_C), names = pltTic, cex.lab=0.8, xlim=c(0.95,ct+0.05), ylim=c(minv, maxv+0.13), yaxt="n")
	axis(2, labels=F)
	
	j <- 1
	print(pltPVa)
	for (i in seq(2.5,ct-1,by=2)) {
		segments(i,-0.15,i,5)
		print(featN[j])
		print(i)
		text(i-1 , maxv+0.11, featN[j], cex=1.22, col=4)

		if (pltPVa[j] <= 0.05) {
			segments(i-1.5, pltMax[j]+0.03, i-0.5, pltMax[j]+0.03)
			segments(i-1.5, pltMax[j]+0.03, i-1.5, pltMax[j]+0.02)
			segments(i-0.5, pltMax[j]+0.03, i-0.5, pltMax[j]+0.02)
			text(i-1, pltMax[j]+0.065, "*", cex=1.5)
		}
		j <- j + 1
	}
	i <- i+2
	text(i-1.0 , maxv+0.11, featN[j], cex=1.22, col=4)
	if (pltPVa[j] <= 0.05) {
		segments(i-1.5, pltMax[j]+0.03, i-0.5, pltMax[j]+0.03)
		segments(i-1.5, pltMax[j]+0.03, i-1.5, pltMax[j]+0.02)
		segments(i-0.5, pltMax[j]+0.03, i-0.5, pltMax[j]+0.02)
		text(i-1, pltMax[j]+0.065, "*", cex=1.5)
	}
}	

#' Compute NSteps Features for TUG data
#'   repeat it for each TUG, 
#'   with size 4
#' @param dTUG list of lists of TUGS
#' @param lTUG list with TUG labels
#' @param nTUGs number of Tugs in the signal (default = 3)
features.TUG.NSteps <- function(dTUG, lTUG, nTUGs=3) {

	# create matrix of features
	features <- matrix(vector(), length(dTUG), 4)

	j <- 1
	for (s in dTUG) {
		#suppress signal outside labels
		m <- 1	
		features[j, m] <- feat.nsteps(s,lTUG[[j]])

		#print(paste(j,features[j,m]))

		lbls <- unique(lTUG[[j]])
		lbls <- lbls[2:length(lbls)]
		
		for (l in lbls) {
			m <- m + 1
			tug <- s[which(lTUG[[j]]==l)]
			features[j, m] <- feat.nsteps(tug, labs=NULL)
		}
		j <- j + 1
	}

	feats <- data.frame(features)
	feats["id"] <- names(dTUG)
	feats["group"] <- groupsTUG()
	return(feats)
}



#' Compute all Features for TUG data
#' -> pse, psp, pspf, wpsp, nsteps
#'   repeat it for each TUG, 
#'   with size 11 x (nTUGS + 1)
#'
#'   
features.TUG <- function(dTUG, lTUG, nTUGs=3, cutFreq=150, first=2) {

	# create matrix of features
	features <- matrix(vector(), length(dTUG), 11*(nTUGs+1) )

	j <- 1
	for (s in dTUG) {
		print(j)
		scopy <- s
		#suppress signal outside labels
		s[which(lTUG[[j]]==0)] <- 0
		# get FT
		s.F <- fourierGait(s, power2n=TRUE)


		# extract global features
		features[j, 1] <- feat.pse(s.F, cutoff=150, first)
		features[j, 2:4] <- feat.psp(s.F, cutoff=cutFreq)
		features[j, 5:7] <- feat.pspf(s.F, cutoff=cutFreq)
		features[j, 8:10] <- feat.wpsp(s.F, cutoff=cutFreq)
		features[j, 11] <- feat.nsteps(scopy, labs=lTUG[[j]])
		
		lbls <- unique(lTUG[[j]])
		
		m <- 11
		for (l in lbls[2:length(lbls)]) {
			tug <- s[which(lTUG[[j]]==l)]
			tug.F <- fourierGait(tug, power2n=TRUE)
			features[j, m+1] <- feat.pse(tug.F, cutoff=cutFreq, first)
			features[j, (m+2):(m+4)] <- feat.psp(tug.F, cutoff=cutFreq)
			features[j, (m+5):(m+7)] <- feat.pspf(tug.F, cutoff=cutFreq)
			features[j, (m+8):(m+10)] <- feat.wpsp(tug.F, cutoff=cutFreq)
			features[j, (m+11)] <- feat.nsteps(tug, labs=NULL)
			m <- m + 11
		}
		j <- j + 1
	}

	feats <- data.frame(features)
	feats["id"] <- names(dTUG)
	feats["group"] <- groupsTUG()
	return(feats)
}


#' Segments a series of TUG trials
#'
#' When performing several TUGs in a row, we need to segment the
#' signal in order to extract characteristics for each trial
#' @param S The original (1-d) signal with all TUGs
#' @param nTUGs The number of TUG trials in the signal
#' @return A vector of labels with 0 values for no data, and 1+ for each segmented TUG
#' @export
segmentTUG <- function(S, nTUGs=3, freq=100) {

	S <- (S - mean(S))^2
	bf <- butter(2, 1/20, type="low")
	S <- filtfilt(bf, S)
	kmed <- (freq/4)-1
	if (kmed%%2 == 0) kmed <- kmed + 1

	S <- rollmedian(S, kmed)
	meanS <- mean(S)
	#print(meanS)

	M <- floor(length(S)/freq)
	seg <- rep(0, length(S))
	lab <- rep(0, M)
	nlab<- rep(0, M)
	l <- 0
	nl <- 0

	#plot(S, t="l")
	#par(new=T)
	#plot(Sm, col="green")	

	# for each second, sum and check
	for (i in 0:(M-1)) {
		idx <- (i*freq)+1
		sumsec <- sum(S[idx:(idx)+(freq-1)])

		#print(paste(i, idx, (idx)+(freq-1)))
		#print(sumsec)

		if (sumsec >= (meanS/3)) {
			if (i > 1 && lab[i-1] != 0) {
				if (lab[i] == 0) {
					lab[i] = l
				}
			} 

			if (i > 0 && lab[i] == 0) {
				l <- l + 1
			}
			nlab[l] <- nlab[l] + 1
			lab[(i+1)] <- l
			#print(paste("seg", lab[i], lab[i+1]))
		} 
	}
	
	nlab2 <- sort(nlab, decreasing=T)
	nlabt <- nlab2[nTUGs]
	for (i in 2:(M-1)) {

		if (lab[i] != 0) {

			if (nlab[lab[i]] < nlabt) {
				lab[i] <- 0
				next
			}
			if (lab[i-1] == 0 && lab[i+1] == 0) {
				lab[i] <- 0
				next
			}
		}
	}

	dlabs <- unique(lab[lab>0])
	ll <- 1
	for (i in dlabs) {
		lab[lab == i] = ll
		ll <- ll + 1		
	}

	for (i in 0:(M-1)) {
		idxf <- (i*freq)+1
		idxl <- (idx)+(freq-1)
		seg[idxf:idxl] <- lab[i+1]
		
	}

	return(seg)
}

#' Get TUG labels from a segmentation algorithm
getTUGlabels <- function(Ss, n, f, N) {

	labels <- vector(mode="list", length=N)
	for (i in 1:N) {
		print(i)
		labels[[i]] <- segmentTUG(Ss[[i]],nTUGs=n,freq=f)
	}
	return(labels)
}

#' Function for local maxima extraction
localmax <- function(x) {
  # Use -Inf instead if x is numeric (non-integer)
  y <- diff(c(-.Machine$integer.max, x)) > 0L
  rle(y)$lengths
  y <- cumsum(rle(y)$lengths)
  y <- y[seq.int(1L, length(y), 2L)]
  if (x[[1]] == x[[2]]) {
    y <- y[-1]
  }
  y
}

#' Test function for the TUG segmentation
testseg <- function(Ss, n, f, N) {

	for (i in 1:N) {
		print(i)
		seg <- segmentTUG(Ss[[i]],nTUGs=n,freq=f)
		lines(seg, col="red")
		line <- readline()
	}
}


#' Test function to find unique values and plot labeled functions
testLabel <- function(d,l,i) {

	plot(d[[i]],t="l",ylim=c(-1,3.4))
	lines(l[[i]], col="red")
	unique(l[[i]])
}


#' Segments a gait walking test signal into two components: walking and turning
#'
#' @param S The original (1-d) signal
#' @return A vector of labels with 0 values for no data, and 1+ for each segmented TUG
#' @export
segmentWalk <- function(S, freq=100, pthr=0.1) {

	#standarize and center the signal
	S <- (S - mean(S))
	
	#filter signal
	bf <- butter(2, 1/20, type="low")
	S <- filtfilt(bf, S)
	
	kmed <- (freq/4)-1
	if (kmed%%2 == 0) kmed <- kmed + 1

	S <- rollmedian(S, kmed)
	meanS <- mean(S)
	dS <- diff(S)
	#plot(S, t="l")
	#plot(dS, t="l", col=2)

	# the threshold is based on 10% of the max value with respect to the mean
	maxv <- max(abs(dS))
	thr <- ((maxv - (mean(dS)))/maxv)*pthr*maxv
	

	segTurnId <- (abs(dS) <= thr)
	#segTurn <- S[segTurnId]
	#segWalk <- S[!segTurnId]
	#lines(segTurn, t="l", col=3)
	#lines(segWalk, t="l", col=2)
	#lines(segTurnId, col=5)

	return(segTurnId)
}

#' Returns all segmented Walks
segmentWalk.All <- function(SS, N=30, f=100,pthr) {

	walkIds <- vector(mode="list", length=N)
	for (i in 1:N) {
		print(i)
		walkIds[[i]] <- segmentWalk(SS[[i]], freq=f, pthr=pthr)
	}
	return(walkIds)

}

#' Compute all Features for Walking Test data
#' -> pse, psp, pspf, wpsp, nsteps
#'   repeat it for each Gait Data, 
#'   with size 11 x 3
#'   
features.Walk <- function(dWalk, lWalk, cutFreq=150, first=2) {

	# create matrix of features
	features <- matrix(vector(), length(lWalk), 33 )

	for (j in 1:length(lWalk)) {
		s <- dWalk[[j]]
		cat(paste(j, ": "))
		scopy <- s
		#suppress signal outside labels
		sTurn <- s[lWalk[[j]]]
		sWalk <- s[!lWalk[[j]]]
		# get FT
		s.F <- fourierGait(s, power2n=TRUE)
		sT.F <- fourierGait(sTurn, power2n=TRUE)
		sW.F <- fourierGait(sWalk, power2n=TRUE)

		# extract global features
		cat("G ")	
		features[j, 1] <- feat.pse(s.F, cutoff=150, first)
		features[j, 2:4] <- feat.psp(s.F, cutoff=cutFreq)
		features[j, 5:7] <- feat.pspf(s.F, cutoff=cutFreq)
		features[j, 8:10] <- feat.wpsp(s.F, cutoff=cutFreq)
		features[j, 11] <- feat.nsteps(scopy)

		# extract features from walk
		cat("W ")	
		features[j, 12] <- feat.pse(sW.F, cutoff=150, first)
		features[j, 13:15] <- feat.psp(sW.F, cutoff=cutFreq)
		features[j, 16:18] <- feat.pspf(sW.F, cutoff=cutFreq)
		features[j, 19:21] <- feat.wpsp(sW.F, cutoff=cutFreq)
		features[j, 22] <- feat.nsteps(sWalk)

		# extract features from turns
		cat("T\n")	
		features[j, 23] <- feat.pse(sT.F, cutoff=150, first)
		features[j, 24:26] <- feat.psp(sT.F, cutoff=cutFreq)
		features[j, 27:29] <- feat.pspf(sT.F, cutoff=cutFreq)
		features[j, 30:32] <- feat.wpsp(sT.F, cutoff=cutFreq)
		features[j, 33] <- feat.nsteps(sTurn)
		#j <- j + 1
	}

	feats <- data.frame(features)
	#feats["id"] <- names(dTUG)
	feats["group"] <- groupsGait()
	return(feats)
}


#' Compute distances between Features for Walking Test data
#'   repeat it for each Gait Data, 
#'   with size 11 x 3
#'   
features.Dist.Walk <- function(feat, first=1, M=11, seqs=3) {

	# create matrix of features
	features <- matrix(vector(), nrow(feat), 3*11 )

	# 1 vs 2
	for (j in 1:M) {
		features[, j] <- (feat[j]-feat[j+M])^2
	}
	# 1 vs 3
	for (j in 1:M) {
		features[, j+M] <- (feat[j]-feat[j+(M*2)])^2
	}
	# 2 vs 3
	for (j in (M+1):(M*2)) {
		print(paste(j,j+M))
		features[, j+M] <- (feat[j]-feat[j+M])^2
	}

	dfeats <- data.frame(features)
	dfeats["group"] <- groupsGait()
	return(dfeats)
}


test.Features <- function(feat) {

	M <- length(feat)
	l <- feat[M]
	for (i in 1:(M-1)) {
		temp <- cbind(feat[i],l)
		names(temp)[1] <- "x1"
		ai <- aov(temp$x1~temp$group)
		
		pval<- summary(ai)[[1]][["Pr(>F)"]][1] 
		cat(paste(i,pval))
		if (pval <= 0.05) cat("***") 
		cat("\n")
	}
}
