# ChIPLine - a pipeline for ChIP-seq analysis

Developers
--------------
Devloped by : Sourya Bhattacharyya

Supervisors: Dr. Ferhat Ay and Dr. Pandurangan Vijayanand

La Jolla Institute for Allergy and Immunology

La Jolla, San Diego, CA 92037, USA



ChIPLine is a pipeline to analysis ChIP-seq data, starting from input Fastq/BAM files and generating alignment summary, various quality statistics, peak calling, and BigWig formatted tracks ready for visualization in UCSC genome browser.


Required packages
-----------------

ChIPLine requires the following packages / libraries to be installed in the system.

1) Bowtie2 (we have used version 2.3.3.1) http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
2) samtools (we have used version 1.6) http://samtools.sourceforge.net/
3) PICARD tools https://broadinstitute.github.io/picard/
4) Package phantompeakqualtools (Developed by Kundaje et al., for analyzing ChIP-seq quality) https://code.google.com/archive/p/phantompeakqualtools/ 
5) Utilities "bedGraphToBigWig", "bedSort", "bigBedToBed", "hubCheck" and "fetchChromSizes" downloaded from UCSC repository. Executables corresponding to the linux system, for example, is provided in this link: http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/
6) deepTools (we have used version 2.0) https://deeptools.readthedocs.io/en/develop/
7) MACS2 (we have used version 2.1.1) https://github.com/taoliu/MACS
8) HOMER (we recommend using the latest version) http://homer.ucsd.edu/homer/
9) R environment (we have used 3.4.3)


User should include the PATH of above mentioned libraries / packages inside their SYSTEM PATH variable. Some of these PATHS are also to be mentioned in a separate configuration file (mentioned below).


Execution
---------

Current package includes a sample script file "pipeline_exec.sh". It conains sample commands required to invoke the main executable named "pipeline.sh", which is provided within the folder "bin".

In general, ChIP-seq pipeline (the executable "pipeline.sh") involves following command line options:

Options:    

Mandatory parameters:

	-C  ConfigFile		    
                 A configuration file to be separately provided. Mandatory parameter. 
                 Current package includes a sample configuration file named "configfile". 
                 Details of the entries in this file are mentioned later.
                  
	-f  FASTQ1          
                Read 1 (or forward strand) of paired-end sequencing data  [.fq|.gz|.bz2]. 
		Or, even an aligned genome (.bam file) can be provided.
            
	-r  FASTQ2          
                R2 of pair-end sequencing data [.fq|.gz|.bz2]. If not provided, and the -f parameter 
                is not a BAM file, the input is assumed to be single ended.
              
	-n  PREFIX           
                Prefix string of output files. For example, -n "TEST" means that the 
                output filenames start with the string "TEST".

	-g  BOWTIE2_GENOME   
                Bowtie2 indexed reference genome. Basically, the folder containing 
		the bwt2 indices are to be provided. 
		Mandatory parameter if user provides fastq files as input (-f and -r options).
		If user provides .bam files as input (-f option) then no need to provide this value.

	-d  OutDir 			  
                Output directory which will contain all the results.

  	-c  CONTROLBAM		 
             	Control file(s) used for peak calling using MACS2. One or more 
		alignment files can be provided to be used 
		as a control. It may not be specified at all, in which 
		case MACS2 operates without any control. 
		Control file can be either in BAM or in  (tagalign.gz) format. 
		If multiple control files are provided, user needs to ensure that all of the 
		control files follow the same format (i.e. either all BAM or all TAGAlign).
		Example: -c control1.bam -c control2.bam puts two control files for using in MACS2.
		
		Conversion from BAM to TagAlign.gz format can be done using the script "TagAlign.sh" 
		provided within the folder "bin".
		
	-w 	BigWigGenome	 
		Reference genome which is used to convert BAM file to a BigWig file. 
		Used for visualization track creation purpose. 
		If -g option is enabled (i.e. the Bowtie2 index genome is provided) 
		then this option is not required. 
		Otherwise, this is a mandatory parameter. Allowed values are 'hg19' 
		(default), 'mm9', 'hg38', and 'mm10'.

    	-T  Tagmentation	 
		If 1, means that Tagmentation was used for ChIP file creation. 
		Then, forward and reverse strands 
		of the current ChIP signal are shifted by the transposon 
		length, and a tagAlign file is generated. 
		Peaks are called from this tagAlign file. Similar to the ATAC seq principle. 
		Applicable for the ChIPMentation technique (Christian Schmidl et. al., 
		ChIPmentation: fast, robust, low-input ChIP-seq for histones and transcription factors, 
		Nature Methods volume 12, pages 963–965, 2015). Default value of this parameter is 0.			
		
 	-D  DEBUG_TXT		 
		Binary variable. If 1 (recommended), different statistics corresponding to 
		quality metrics and reads are printed. Useful when a summary of a large set 
		of ChIP-seq samples are to be generated.
		
	-q  MAPQ_THR		 
		Quality value threshold, below which the mapped reads are removed (Default 30).
		
	-p  PEAKCALLGENOMESIZE 
		genome size parameter for MACS2 peak calling ("hs", "mm", "ce", "dm": default "hs")
		

Optional parameters:

	-O 	Overwrite		 
		Binary variable. If 1, overwrites the existing files (if any). Default = 0.
						 
	-t  NUMTHREADS              
		Number of sorting, Bowtie2 mapping THREADS [Default = 1]. For parallel processing of Bowtie2, 
		user should specify > 1 value such as 4 ot 8.
		
	-m  MAX_MEM          
		Set max memory used for PICARD duplication removal [Default = 8G].
		
	-a  ALIGNVALIDMAX	 
		Set the number of (max) valid alignments which will be searched [Default = 4] 
		for Bowtie2.
		
	-l  MAXFRAGLEN 		 
		Set the maximum fragment length to be used for Bowtie2 alignment [Default = 2000]















