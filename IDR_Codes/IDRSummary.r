#!/usr/bin/env Rscript

#===========================================================
# R script for summarizing the results of IDR analysis between different sample replicates

#Author: Sourya Bhattacharyya
#Vijay-Ay lab, LJI

# usage: Rscript result_summary.r $inpfile
#===========================================================

args <- commandArgs(TRUE)

# file containing the overlapped peak information
CommonPeakFile <- args[1]
npeak1 <- as.integer(args[2])
npeak2 <- as.integer(args[3])

inpdir <- dirname(CommonPeakFile)

## read the common peaks
## by default, it'll remove the row names
CommonPeakInfo <- read.table(CommonPeakFile, header=TRUE)

## error check
## sometimes, the output file may not have proper chromosome or interval information
## match chromosome names, 
CommonPeakInfo <- CommonPeakInfo[which(CommonPeakInfo[,1] == CommonPeakInfo[,5]), ]

# number of overlapped peaks (considering all the IDR values)
ncommonpeak <- nrow(CommonPeakInfo)
fracpeak1 <- (ncommonpeak * 1.0) / npeak1
fracpeak2 <- (ncommonpeak * 1.0) / npeak2

# peaks with IDR <= specified thresholds (0.05 and 0.1) 
NumIDRPass1 <- length(which(CommonPeakInfo[,10] <= 0.05))
FracIDRPass1 <- (NumIDRPass1 * 1.0) / ncommonpeak
NumIDRPass2 <- length(which(CommonPeakInfo[,10] <= 0.1))
FracIDRPass2 <- (NumIDRPass2 * 1.0) / ncommonpeak

# divide the input overlapped peak files into two different structures
# corresponding to the peak information of two different inputs
# the seq() function also includes the row number for every interaction
# this row number serves as the id of peaks
PeakInfoInput1 <- cbind(seq(1:ncommonpeak), CommonPeakInfo[,1:4])
PeakInfoInput2 <- cbind(seq(1:ncommonpeak), CommonPeakInfo[,5:8])

# sort the data according to the significance value (last column of both the data)
# decreasing order is employed
PeakInfoInput1_Sort <- PeakInfoInput1[ order(-PeakInfoInput1[,5]),]
PeakInfoInput2_Sort <- PeakInfoInput2[ order(-PeakInfoInput2[,5]),]

# we check the cumulative percent of samples in both peak sets
# and find out the overlap of peaks
fraction_overlap <- c()

for (x in seq(0, 1, 0.1)) {
	if ((x != 0) && (x != 1)) {
		# number of elements of both peak lists
		nsample <- as.integer(ncommonpeak * x)
		# common elements in both peak lists
		# the common factor is the first column: peak id
		OverlapSet <- PeakInfoInput1_Sort[1:nsample, 1] %in% PeakInfoInput2_Sort[1:nsample, 1]
		ncommon <- length(OverlapSet[OverlapSet==TRUE])
		frac_common <- (ncommon * 1.0 / nsample)
		fraction_overlap <- c(fraction_overlap, frac_common)

		# we also note down two different fraction overlap statistics
		# corresponding to 10\%, 20% and 50% strongest peaks
		if (x == 0.1) {
			frac_overlap_10Pct = frac_common
		}
		if (x == 0.2) {
			frac_overlap_20Pct = frac_common
		}
		if (x == 0.5) {
			frac_overlap_50Pct = frac_common
		}
	}
}

# write the results in a text file
OutFilename <- paste0(inpdir, '/Stat.tab')

fp <- file(OutFilename, open="w")
write(paste0('NPeak1', '\t', 'NPeak2', '\t', 'CommonPeak', '\t', 'FracPeak1', '\t', 'FracPeak2', '\t', 'IDR_0.05_Peak', '\t', 'Frac_IDR_0.05_Peak', '\t', 'IDR_0.1_Peak', '\t', 'Frac_IDR_0.1_Peak', '\t', 'MeanOverlap', '\t', 'Overlap10', '\t', 'Overlap20', '\t', 'Overlap50'), file=fp, append=T)
write(paste(npeak1, '\t', npeak2, '\t', ncommonpeak, '\t', fracpeak1, '\t', fracpeak2, '\t', NumIDRPass1, '\t', FracIDRPass1, '\t', NumIDRPass2, '\t', FracIDRPass2, '\t', mean(fraction_overlap), '\t', frac_overlap_10Pct, '\t', frac_overlap_20Pct, '\t', frac_overlap_50Pct), file=fp, append=T)
close(fp)


##========================
## dump the significant peaks with respect to specific IDR thresholds
##========================

IDROutDir <- paste0(dirname(CommonPeakFile), '/FINAL_IDR_Peaks')
system(paste("mkdir -p", IDROutDir))

for (IDRThr in c(0.05, 0.1)) {
	outfile <- paste0(IDROutDir, '/IDR_Peaks_FDR', IDRThr, '.txt')
	idx <- which(CommonPeakInfo[, ncol(CommonPeakInfo)] <= IDRThr)
	write.table(CommonPeakInfo[idx, ], outfile, row.names=F, col.names=T, sep="\t", quote=F, append=F)
}






