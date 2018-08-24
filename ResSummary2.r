# sample R script file
# used to determine the summary statistics of various ChIP seq samples
# within a particular experiment (run)

# author: Sourya Bhattacharyya
# Vijay-AY lab

#==============
# Execution command:
# Rscript ResSummary2.r $OutBaseDir $BAMRead $Tagmentation $OldMethod
# input arguments:
# 1) OutBaseDir: Directory under which results of all the different samples are stored
# 2) BAMRead: Boolean variable, indicating if the BAM files were used for analysis (1) or fastq (0). Default 0
# 3) Tagmentation: Binary variable. If 1, denotes that the data samples have been under ChIPMentation method. Default 0.
# 4) OldMethod: Binary variable. If 1, old summary statistics is used for performance evaluation. Default 0 (preferred).
# 5) ControlPeak: 	Variable with value of either 0, 1 or 2. If 1, control samples has been used for peak calling. 
# 					If 0, no control sample is used. If 2, peaks using control and not using control both exist. Default 0.

# sample execution command:
# Rscript ResSummary2.r /home/sourya/ChIPResults/ 0 1 0 0
#==============

args <- commandArgs(TRUE)
if(length(args)<1)	{
	q("no")
}

# directory under which all the results have been stored
# all the subdirectories indicate different samples
# corresponding to the experimental condition mentioned 
# in the main directory
baseresdir <- args[1]

# check if the last character of the output base directory is '/'
# unless append that character
if (substr(baseresdir,nchar(baseresdir),nchar(baseresdir)) != '/') {
	baseresdir <- paste0(baseresdir, '/')
}

# denotes whether the analysis started from BAM files
# in such a case, FASTQ reads are not available
if (length(args) > 1) {
	BAMRead <- as.integer(args[2])	
} else {
	cat(sprintf("\n User did not provide information whether the BAM file was the starting point -- assuming 0, means that fastq files were used"))
	BAMRead <- 0
}

# denotes whether the data samples have been under ChIPMentation method
if (length(args) > 2) {
	Tagmentation <- as.integer(args[3])	
} else {
	cat(sprintf("\n User did not provide information whether the data is ChIPMentation (1) or ChIP seq (0) -- assuming 0"))
	Tagmentation <- 0
}

# denotes, if 1, whether old method of summary statistics would be used
if (length(args) > 3) {
	OldMethod <- as.integer(args[4])	
} else {
	OldMethod <- 0
}

# This variable is applicable for new settings only
# 0: peak is derived without any control (default)
# 1: peak is derived with control
# 2: peak is derived with both control and no control
if (length(args) > 4) {
	ControlPeak <- as.integer(args[4])	
} else {
	ControlPeak <- 0
}

#==========================
# following are list of file and directory name conventions which are used in the ChIP seq pipeline
NRFfilenamefmt <- 'out_NRF_MAPQ30.txt'
QCStrandCorrNameFmt <- 'chipSampleMaster.tagAlign.tab'
FRiPFileNameFmt <- 'out_FRiP.txt'
PeakCountFileFmt <- 'Peak_Statistics.txt'
ReadStatFileNameFmt <- 'Read_Count_Stat.txt'

if (OldMethod == 1) {
	if (Tagmentation == 1) {
		MACS2DirControlNameFmt <- 'MACS2_Tag_with_Control'
		MACS2DirNoControlNameFmt <- 'MACS2_Tag_No_Control'
	} else {
		MACS2DirControlNameFmt <- 'MACS2_with_Control'
		MACS2DirNoControlNameFmt <- 'MACS2_No_Control'	
	}
} else {
	if (Tagmentation == 1) {
		MACS2DirDefControlNameFmt <- 'MACS2_Default_Tag_with_Control'
		MACS2DirExtControlNameFmt <- 'MACS2_Ext_Tag_with_Control'
		MACS2DirDefNoControlNameFmt <- 'MACS2_Default_Tag_No_Control'
		MACS2DirExtNoControlNameFmt <- 'MACS2_Ext_Tag_No_Control'
	} else {
		MACS2DirDefControlNameFmt <- 'MACS2_Default_with_Control'
		MACS2DirExtControlNameFmt <- 'MACS2_Ext_with_Control'
		MACS2DirDefNoControlNameFmt <- 'MACS2_Default_No_Control'
		MACS2DirExtNoControlNameFmt <- 'MACS2_Ext_No_Control'
	}
}

# output file containing the summary results
if (OldMethod == 1) {
	outtextfile <- paste0(baseresdir, 'Results_All_Samples_Summary.txt')
	con <- file(outtextfile, "w")
	cat(paste(c("Dir", "NumRead", "UniqMappRead", "UniqMapPos", "NRF", "NSC", "RSC", "Qtag", "MappedReadPeak(NoCtrl)", "FRiP(NoCtrl)", "nPeakNoCtrl", "nPeak(Q<0.05)NoCtrl", "nPeak(Q<0.01)NoCtrl", "MappedReadPeak(Ctrl)", "FRiP(Ctrl)", "nPeakCtrl", "nPeak(Q<0.05)Ctrl", "nPeak(Q<0.01)Ctrl"), collapse='\t'), file=con)	
} else {
	outtextfile <- paste0(baseresdir, 'Results_Summary_NEW.xls')
	con <- file(outtextfile, "w")

	# initialize the header of the summary excel file
	HeaderVec <- c("Dir", "FastqRead", "MappRead")

	if (BAMRead == 0) {
		HeaderVec <- c(HeaderVec, "%MappRead", "RandomDelRead", "%RandomDelRead", "UniqMappRead", "%UniqMappRead", "MapqThrRead", "%MapqThrRead", "rmDupRead", "%rmDupRead")
	} else {
		HeaderVec <- c(HeaderVec, "%MappRead(w.r.t FastqRead)", "RandomDelRead", "%RandomDelRead (w.r.t MappRead)", "UniqMappRead", "%UniqMappRead (w.r.t MappRead)", "MapqThrRead", "%MapqThrRead (w.r.t MappRead)", "rmDupRead", "%rmDupRead (w.r.t MappRead)")
	}

	HeaderVec <- c(HeaderVec, "UniqMapPos", "NRF", "M1", "M2", "PBC1", "PBC2", "NSC", "RSC", "Qtag")

	if ((ControlPeak == 0) | (ControlPeak == 2)) {
		HeaderVec <- c(HeaderVec, "MappedReadPeak(Def_NoCtrl)", "FRiP(Def_NoCtrl)", "nPeak(Def_NoCtrl)", "nPeak(Q<0.05)(Def_NoCtrl)", "nPeak(Q<0.01)(Def_NoCtrl)", "MappedReadPeak(Ext_NoCtrl)", "FRiP(Ext_NoCtrl)", "nPeak(Ext_NoCtrl)", "nPeak(Q<0.05)(Ext_NoCtrl)", "nPeak(Q<0.01)(Ext_NoCtrl)")
	}
	if ((ControlPeak == 1) | (ControlPeak == 2)) {
		HeaderVec <- c(HeaderVec, "MappedReadPeak(Def_Ctrl)", "FRiP(Def_Ctrl)", "nPeak(Def_Ctrl)", "nPeak(Q<0.05)(Def_Ctrl)", "nPeak(Q<0.01)(Def_Ctrl)", "MappedReadPeak(Ext_Ctrl)", "FRiP(Ext_Ctrl)", "nPeak(Ext_Ctrl)", "nPeak(Q<0.05)(Ext_Ctrl)", "nPeak(Q<0.01)(Ext_Ctrl)")
	}
	cat(paste(HeaderVec, collapse='\t'), file=con)
}

#=================================
# process individual directories under the main results directory
dir.list <- list.dirs(path = baseresdir, full.names = FALSE, recursive = FALSE)
for (dr in dir.list) {
	
	# condition whether old summary method would be used
	if (OldMethod == 0) {

		# this vector stores the output results
		# same length as the header vector
		OutResVec <- rep("NA", length(HeaderVec))

		# directory
		OutResVec[1] <- basename(dr)

		# read the text file containing read count statistic
		ReadStatFile <- paste0(baseresdir, dr, "/", ReadStatFileNameFmt)
		if (file.exists(ReadStatFile) && (file.access(ReadStatFile, 4) == 0)){

			cat(sprintf("\n Found the file: %s \n", ReadStatFile))

			x <- readLines(ReadStatFile)
			lastline <- strsplit(x[length(x)], "\t")[[1]]

			# the line can have 8 or 7 fields
			# 8 fields if fastq file is used in the pipeline
			# 7 fields if already aligned file is used in the pipeline
			nfields <- length(lastline)
			cat(sprintf("\n No of fields in the read statistics file: %s ", nfields))

			if (BAMRead == 0) {
			 	OutResVec[2] <- lastline[nfields-6]	#lastline[2]
			 	OutResVec[3] <- lastline[nfields-5]	#lastline[3]
			 	OutResVec[4] <- as.double(OutResVec[3]) / as.double(OutResVec[2])				
			} else {
				OutResVec[3] <- lastline[nfields-5]
			}
		 	# RandomDelRead <- lastline[nfields-4]
		 	# Percent_RandomDelRead <- as.double(RandomDelRead) / as.double(OutResVec[2])
		 	OutResVec[5] <- lastline[nfields-3]	#lastline[4]
		 	if (BAMRead == 0) {
		 		OutResVec[6] <- as.double(OutResVec[5]) / as.double(OutResVec[2])
	 		} else {
	 			OutResVec[6] <- as.double(OutResVec[5]) / as.double(OutResVec[3])
	 		}
		 	OutResVec[7] <- lastline[nfields-2]	#lastline[5]
		 	if (BAMRead == 0) {
		 		OutResVec[8] <- as.double(OutResVec[7]) / as.double(OutResVec[2])
	 		} else {
	 			OutResVec[8] <- as.double(OutResVec[7]) / as.double(OutResVec[3])
	 		}
		 	OutResVec[9] <- lastline[nfields-1]	#lastline[6]
		 	if (BAMRead == 0) {
		 		OutResVec[10] <- as.double(OutResVec[9]) / as.double(OutResVec[2])	
		 	} else {
		 		OutResVec[10] <- as.double(OutResVec[9]) / as.double(OutResVec[3])
		 	}
		 	OutResVec[11] <- lastline[nfields]	#lastline[7]
		 	if (BAMRead == 0) {
		 		OutResVec[12] <- as.double(OutResVec[11]) / as.double(OutResVec[2])
			} else {
				OutResVec[12] <- as.double(OutResVec[11]) / as.double(OutResVec[3])
			}
		}

		# the following file stores the NRF / library complexity value
		NRF_textfile <- paste0(baseresdir, dr, "/", NRFfilenamefmt)

		if (file.exists(NRF_textfile) && (file.access(NRF_textfile, 4) == 0)){
			# open the file and read line by line
			# extract performance values
			finp <- file(NRF_textfile, "r")
			lineset <- readLines(finp)
			for (i in 1:length(lineset)) {
				curr_line <- lineset[i]
				if (regexpr("Unique_Genome_Pos", curr_line) > 0) {
					OutResVec[13] <- as.integer(strsplit(curr_line,":")[[1]][2])
				} else if (regexpr("NRF", curr_line) > 0) {
					OutResVec[14] <- as.numeric(strsplit(curr_line,":")[[1]][2])
				} else if (regexpr("M1", curr_line) > 0) {
					OutResVec[15] <- as.integer(strsplit(curr_line,":")[[1]][2])
				} else if (regexpr("M2", curr_line) > 0) {
					OutResVec[16] <- as.integer(strsplit(curr_line,":")[[1]][2])
				} else if (regexpr("PBC1", curr_line) > 0) {
					OutResVec[17] <- as.numeric(strsplit(curr_line,":")[[1]][2])
				} else if (regexpr("PBC2", curr_line) > 0) {
					OutResVec[18] <- as.numeric(strsplit(curr_line,":")[[1]][2])
				}
			}
			
			# close the input file
			close(finp)
		}

		# following file stores the NSC, RSC, metrics
	  	tabfile <- paste0(baseresdir, dr, "/", QCStrandCorrNameFmt)
	  	if (file.exists(tabfile) && (file.access(tabfile, 4) == 0)){
	  		#print(sprintf("\n Tab text file exists"))
	  		x <- readLines(tabfile)
	  		lastline <- strsplit(x[length(x)], "\t")[[1]]
	  		OutResVec[19] <- as.numeric(lastline[9])
	  		OutResVec[20] <- as.numeric(lastline[10])
	  		OutResVec[21] <- as.integer(lastline[11])
		}

		#=======================================
		# insert the peak information 
		counter <- 22

		# peaks when no control is used
		if ((ControlPeak == 0) | (ControlPeak == 2)) {
			
			# insert statistics for default MACS2 peaks without control
			# FRIP
			fripfile <- paste0(baseresdir, dr, "/", MACS2DirDefNoControlNameFmt, "/", FRiPFileNameFmt)
			if (file.exists(fripfile) && (file.access(fripfile, 4) == 0)) {
				x <- readLines(fripfile)
				lastline <- strsplit(x[length(x)], "\t")[[1]]
				OutResVec[counter] <- lastline[2]
				OutResVec[counter+1] <- lastline[3]
			}
			counter <- counter + 2

			# peak count
			peakstatfile <- paste0(baseresdir, dr, "/", MACS2DirDefNoControlNameFmt, "/", PeakCountFileFmt)
			if (file.exists(peakstatfile) && (file.access(peakstatfile, 4) == 0)) {
				x <- readLines(peakstatfile)
				lastline <- strsplit(x[length(x)], "\t")[[1]]
				OutResVec[counter] <- lastline[1]
				OutResVec[counter+1] <- lastline[2]
				OutResVec[counter+2] <- lastline[3]
			}
			counter <- counter + 3

			# insert statistics for extsize MACS2 peaks without control
			# FRIP
			fripfile <- paste0(baseresdir, dr, "/", MACS2DirExtNoControlNameFmt, "/", FRiPFileNameFmt)
			if (file.exists(fripfile) && (file.access(fripfile, 4) == 0)) {
				x <- readLines(fripfile)
				lastline <- strsplit(x[length(x)], "\t")[[1]]
				OutResVec[counter] <- lastline[2]
				OutResVec[counter+1] <- lastline[3]
			}
			counter <- counter + 2

			# peak count
			peakstatfile <- paste0(baseresdir, dr, "/", MACS2DirExtNoControlNameFmt, "/", PeakCountFileFmt)
			if (file.exists(peakstatfile) && (file.access(peakstatfile, 4) == 0)) {
				x <- readLines(peakstatfile)
				lastline <- strsplit(x[length(x)], "\t")[[1]]
				OutResVec[counter] <- lastline[1]
				OutResVec[counter+1] <- lastline[2]
				OutResVec[counter+2] <- lastline[3]
			}
			counter <- counter + 3

		}	# end peak check for MACS2 peaks without control

		# peaks when control is used
		if ((ControlPeak == 1) | (ControlPeak == 2)) {
			
			# insert statistics for default MACS2 peaks with control
			# FRIP
			fripfile <- paste0(baseresdir, dr, "/", MACS2DirDefControlNameFmt, "/", FRiPFileNameFmt)
			if (file.exists(fripfile) && (file.access(fripfile, 4) == 0)) {
				x <- readLines(fripfile)
				lastline <- strsplit(x[length(x)], "\t")[[1]]
				OutResVec[counter] <- lastline[2]
				OutResVec[counter+1] <- lastline[3]
			}
			counter <- counter + 2

			# peak count
			peakstatfile <- paste0(baseresdir, dr, "/", MACS2DirDefControlNameFmt, "/", PeakCountFileFmt)
			if (file.exists(peakstatfile) && (file.access(peakstatfile, 4) == 0)) {
				x <- readLines(peakstatfile)
				lastline <- strsplit(x[length(x)], "\t")[[1]]
				OutResVec[counter] <- lastline[1]
				OutResVec[counter+1] <- lastline[2]
				OutResVec[counter+2] <- lastline[3]
			}
			counter <- counter + 3

			# insert statistics for extsize MACS2 peaks with control
			# FRIP
			fripfile <- paste0(baseresdir, dr, "/", MACS2DirExtControlNameFmt, "/", FRiPFileNameFmt)
			if (file.exists(fripfile) && (file.access(fripfile, 4) == 0)) {
				x <- readLines(fripfile)
				lastline <- strsplit(x[length(x)], "\t")[[1]]
				OutResVec[counter] <- lastline[2]
				OutResVec[counter+1] <- lastline[3]
			}
			counter <- counter + 2

			# peak count
			peakstatfile <- paste0(baseresdir, dr, "/", MACS2DirExtControlNameFmt, "/", PeakCountFileFmt)
			if (file.exists(peakstatfile) && (file.access(peakstatfile, 4) == 0)) {
				x <- readLines(peakstatfile)
				lastline <- strsplit(x[length(x)], "\t")[[1]]
				OutResVec[counter] <- lastline[1]
				OutResVec[counter+1] <- lastline[2]
				OutResVec[counter+2] <- lastline[3]
			}
			counter <- counter + 3

		}	# end peak check for MACS2 peaks with control

		# write the statistic in the results log file
		cat('\n', paste(OutResVec, collapse='\t'), file=con)		

	} else {
		#==================
		# the following file stores the NRF / library complexity value
		# at the third column of the last line
		#==================
		NRF_textfile <- paste0(baseresdir, dr, "/", NRFfilenamefmt)

		if (file.exists(NRF_textfile) && (file.access(NRF_textfile, 4) == 0)){
		 	x <- readLines(NRF_textfile)
		 	# the 2nd line in string splitted structure
		 	lastline <- strsplit(x[length(x)], "\t")[[1]]
		 	NumRead <- lastline[1]
		 	UniqMappedRead <- lastline[2]
		 	UniqMappedPos <- lastline[3]
		  	NRF_val <- lastline[4]
	 	} else {
	 		NumRead <- 'NA'
	 	 	UniqMappedRead <- 'NA'
		 	UniqMappedPos <- 'NA'		
		 	NRF_val <- 'NA'
	 	}
	 	#==================
	  	# following file stores the NSC, RSC, metrics
	  	#==================
	  	tabfile <- paste0(baseresdir, dr, "/", QCStrandCorrNameFmt)
	  	if (file.exists(tabfile) && (file.access(tabfile, 4) == 0)){
	  		#print(sprintf("\n Tab text file exists"))
	  		x <- readLines(tabfile)
	  		lastline <- strsplit(x[length(x)], "\t")[[1]]
	  		NSC_val <- lastline[9]
	  		RSC_val <- lastline[10]
	  		Qtag_val <- lastline[11]
		} else {
	  		NSC_val <- 'NA'
	  		RSC_val <- 'NA'
	  		Qtag_val <- 'NA'
		}

		#==================
		# FRiP and Peak count measures
		# when MACS2 is executed with input control BAM files
		#==================	
		FRiP_textfile_ctrl <- paste0(baseresdir, dr, "/", MACS2DirControlNameFmt, "/", FRiPFileNameFmt)
		if (file.exists(FRiP_textfile_ctrl) && (file.access(FRiP_textfile_ctrl, 4) == 0)){
		 	x <- readLines(FRiP_textfile_ctrl)
		 	lastline <- strsplit(x[length(x)], "\t")[[1]]
		 	MappedReadPeakCtrl <- lastline[2]
		  	FRiP_val_ctrl <- lastline[3]
		} else {
			MappedReadPeakCtrl <- 'NA'
			FRiP_val_ctrl <- 'NA'
		}

		PeakCount_TextFile_ctrl <- paste0(baseresdir, dr, "/", MACS2DirControlNameFmt, "/", PeakCountFileFmt)
		if (file.exists(PeakCount_TextFile_ctrl) && (file.access(PeakCount_TextFile_ctrl, 4) == 0)){
		 	x <- readLines(PeakCount_TextFile_ctrl)
		 	lastline <- strsplit(x[length(x)], "\t")[[1]]
		 	TotPeakCtrl <- lastline[1]
		  	TotPeak_Q_Five_Pct_ctrl <- lastline[2]
		  	TotPeak_Q_One_Pct_ctrl <- lastline[3]
		} else {
			TotPeakCtrl <- 'NA'
			TotPeak_Q_Five_Pct_ctrl <- 'NA'
			TotPeak_Q_One_Pct_ctrl <- 'NA'
		}

		#==================
		# FRiP and Peak count measures
		# when MACS2 is executed with default (without any control bam file)
		#==================	
		FRiP_textfile_Noctrl <- paste0(baseresdir, dr, "/", MACS2DirNoControlNameFmt, "/", FRiPFileNameFmt)
		if (file.exists(FRiP_textfile_Noctrl) && (file.access(FRiP_textfile_Noctrl, 4) == 0)){
		 	x <- readLines(FRiP_textfile_Noctrl)
		 	lastline <- strsplit(x[length(x)], "\t")[[1]]
		 	MappedReadPeak <- lastline[2]
		  	FRiP_val <- lastline[3]
		} else {
			MappedReadPeak <- 'NA'
			FRiP_val <- 'NA'
		}

		PeakCount_TextFile_Noctrl <- paste0(baseresdir, dr, "/", MACS2DirNoControlNameFmt, "/", PeakCountFileFmt)
		if (file.exists(PeakCount_TextFile_Noctrl) && (file.access(PeakCount_TextFile_Noctrl, 4) == 0)){
		 	x <- readLines(PeakCount_TextFile_Noctrl)
		 	lastline <- strsplit(x[length(x)], "\t")[[1]]
		 	TotPeak <- lastline[1]
		  	TotPeak_Q_Five_Pct <- lastline[2]
		  	TotPeak_Q_One_Pct <- lastline[3]
		} else {
			TotPeak <- 'NA'
			TotPeak_Q_Five_Pct <- 'NA'
			TotPeak_Q_One_Pct <- 'NA'
		}

		#==============================================
		# write the statistic in the results log file
		cat('\n', paste(c(basename(dr), NumRead, UniqMappedRead, UniqMappedPos, NRF_val, NSC_val, RSC_val, Qtag_val, MappedReadPeak, FRiP_val, TotPeak, TotPeak_Q_Five_Pct, TotPeak_Q_One_Pct, MappedReadPeakCtrl, FRiP_val_ctrl, TotPeakCtrl, TotPeak_Q_Five_Pct_ctrl, TotPeak_Q_One_Pct_ctrl), collapse='\t'), file=con) 

	}	# end process individual directories
}	# end old method condition


if (OldMethod == 0) {
	#============================
	# a few summary statements for the output text file
	#============================
	outtext <- paste0("\n\n\n\n\n *** Important parameters ***** \n\n\n")
	writeLines(outtext, con=con, sep="\n")

	outtext <- paste0("\n\n  FastqRead: number of reads in individual fastq file(s)")
	writeLines(outtext, con=con, sep="\n")

	outtext <- paste0("\n\n  MappRead: number of mappable reads.")
	writeLines(outtext, con=con, sep="\n")

	outtext <- paste0("\n\n  %MappRead: percentage of mappable reads.")
	writeLines(outtext, con=con, sep="\n")

	outtext <- paste0("\n\n  RandomDelRead: number of reads REMAINING after deleting random chromosome reads.")
	writeLines(outtext, con=con, sep="\n")

	outtext <- paste0("\n\n  %RandomDelRead: percentage of reads REMAINING after deleting random chromosome reads.")
	writeLines(outtext, con=con, sep="\n")

	outtext <- paste0("\n\n  UniqMappRead and %UniqMappRead: number (and percentage) of reads uniquely mapped to the reference genome.")
	writeLines(outtext, con=con, sep="\n")

	outtext <- paste0("\n\n  MapqThrRead and %MapqThrRead: number (and percentage) of reads REMAINING after quality thresholding (default MAPQ=30)")
	writeLines(outtext, con=con, sep="\n")

	outtext <- paste0("\n\n rmDupRead and %rmDupRead: number (and percentage) of reads REMAINING after deleting the deplicated reads")
	writeLines(outtext, con=con, sep="\n")

	outtext <- paste0("\n\n  UniqMapPos: number of distinct genome position where at least one read maps uniquely.")
	writeLines(outtext, con=con, sep="\n")

	outtext <- paste0("\n\n  NRF (Non redundant fraction): number of distinct genome positions for uniquely mapped reads / number of uniquely mapped reads =============== < 0.5: concerning,  between 0.5 and 0.8 : acceptable, between 0.8 and 0.9 : compliant, > 0.9: ideal")
	writeLines(outtext, con=con, sep="\n")

	outtext <- paste0("\n\n  M1: number of genomic locations where exactly one read maps uniquely.")
	writeLines(outtext, con=con, sep="\n")

	outtext <- paste0("\n\n  M2: number of genomic locations where exactly two reads map uniquely.")
	writeLines(outtext, con=con, sep="\n")

	outtext <- paste0("\n\n  PBC1 (PCR bottlenecking coefficient 1): M1 / UniqMapPos ===============  < 0.5: severe bottlenecking level,  between 0.5 and 0.8: moderate bottlenecking level, between 0.8 and 0.9 : mild bottlenecking level, >= 0.9 : no bottlenecking ")
	writeLines(outtext, con=con, sep="\n")

	outtext <- paste0("\n\n  PBC2: M1 / M2 ===============  < 1: severe bottlenecking level,  between 1 and 3: moderate bottlenecking level, between 3 and 10 : mild bottlenecking level, >= 10 : no bottlenecking ")
	writeLines(outtext, con=con, sep="\n")


	outtext <- paste0("\n\n  NSC: Normalized strand correlation --- acceptable if NSC >= 1.05 ")
	writeLines(outtext, con=con, sep="\n")

	outtext <- paste0("\n\n  RSC: relative strand correlation --- acceptable if NSC >= 0.8 ")
	writeLines(outtext, con=con, sep="\n")

	outtext <- paste0("\n\n  Quality tag: negative: reject sample, 0: marginally acceptable, 1: very good sample, 2: excellent quality ")
	writeLines(outtext, con=con, sep="\n")

	

	outtext <- paste0("\n\n\n\n\n  MACS2 outputs without control input ----- missing values are replaced by NA \n\n")
	writeLines(outtext, con=con, sep="\n")

	outtext <- paste0("\n\n MappedReadPeak(NoCtrl): mapped reads in peaks ")
	writeLines(outtext, con=con, sep="\n")

	outtext <- paste0("\n\n FRiP(NoCtrl): MappedReadPeak(NoCtrl) / UniqMapRead ")
	writeLines(outtext, con=con, sep="\n")

	outtext <- paste0("\n\n nPeakNoCtrl: number of peaks (determined by p value threshold of 0.01) ")
	writeLines(outtext, con=con, sep="\n")

	outtext <- paste0("\n\n nPeak(Q<0.05)NoCtrl: number of peaks (determined by q value threshold of 0.05) ")
	writeLines(outtext, con=con, sep="\n")

	outtext <- paste0("\n\n nPeak(Q<0.01)NoCtrl: number of peaks (determined by q value threshold of 0.01) ")
	writeLines(outtext, con=con, sep="\n")


	outtext <- paste0("\n\n\n\n\n  MACS2 outputs with control input ----- missing values are replaced by NA \n\n")
	writeLines(outtext, con=con, sep="\n")

	outtext <- paste0("\n\n MappedReadPeak(Ctrl): mapped reads in peaks ")
	writeLines(outtext, con=con, sep="\n")

	outtext <- paste0("\n\n FRiP(Ctrl): MappedReadPeak(Ctrl) / UniqMapRead ")
	writeLines(outtext, con=con, sep="\n")

	outtext <- paste0("\n\n nPeakCtrl: number of peaks (determined by p value threshold of 0.01) ")
	writeLines(outtext, con=con, sep="\n")

	outtext <- paste0("\n\n nPeak(Q<0.05)Ctrl: number of peaks (determined by q value threshold of 0.05) ")
	writeLines(outtext, con=con, sep="\n")

	outtext <- paste0("\n\n nPeak(Q<0.01)Ctrl: number of peaks (determined by q value threshold of 0.01) ")
	writeLines(outtext, con=con, sep="\n")

}	# old method condition

# close the output summary text file
close(con)



