#!/usr/bin/env Rscript

#===========================================================
# R script for finding the ChIP-seq quality metrics
# using the package "phantompeakqualtools"
# from Anshul Kundaje et al.
# and the corresponding source code spp
# this code is invoked from the main ChIP pipeline

# Author: Sourya Bhattacharyya
# Vijay-Ay lab, LJI
#===========================================================

args <- commandArgs(TRUE)
if(length(args)<1)	{
	q("no")
}

# first parameter is the R executable
RPackageExec <- args[1]

# next parameter is the SPP package executable 
# using the package "phantompeakqualtools"
# from Anshul Kundaje et al.
sppexec <- args[2]

# third parameter is the output directory
OutDir <- args[3]

# check if the last character of the output directory is '/'
# unless append that character
if (substr(OutDir,nchar(OutDir),nchar(OutDir)) != '/') {
	OutDir <- paste0(OutDir, '/')
}

# all the subsequent parameters store one or more alignment files
# which are combined and processed as a single file
allInpFiles <- ''
for (i in (4:length(args))) {
	inpbamfile <- args[i]
	# the following file converts the input bam file
	filename <- paste0(OutDir, "chipSampleRep", (i-3), ".tagAlign.gz")
	system(paste("samtools view -F 0x0204 -o -", inpbamfile, "| awk \'BEGIN{OFS=\"\t\"}{if (and($2,16) > 0) {print $3,($4-1),($4-1+length($10)),\"N\",\"1000\",\"-\"} else {print $3,($4-1),($4-1+length($10)),\"N\",\"1000\",\"+\"}}\' | gzip -c >", filename))
	allInpFiles <- paste0(allInpFiles, ' ', filename, ' ')
}

cat(sprintf("\n OutDir : %s ", OutDir))
cat(sprintf("\n allInpFiles : %s ", allInpFiles))

# this file will contain the zipped archieve of all the replicates
combined_filename <- paste0(OutDir, 'chipSampleMaster.tagAlign.gz')
tab_filename <- paste0(OutDir, 'chipSampleMaster.tagAlign.tab')

# combine the input bam files to produce the final alignment file
# for quality
system(paste("zcat", allInpFiles, "|gzip -c >", combined_filename))

cat(sprintf("\n\n Analysing cross-correlation and fragment length\n"))

# sourya - we overwrite the plots for the moment, y=using the -rf option
# it can be removed later
system(paste0(RPackageExec, " ", sppexec, " -rf -s=-100:5:600 -c=", combined_filename, " -savp -out=", tab_filename))

# remove the temporary files
for (i in (4:length(args))) {
	filename <- paste0(OutDir, "chipSampleRep", (i-3), ".tagAlign.gz")
	system(paste("rm", filename))
}




