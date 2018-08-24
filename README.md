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
                Bowtie2 indexed reference genome. Basically, the folder containing the bwt2 indices are to be provided.

	-d  OutDir 			  
                Output directory which will contain all the results.

  -c  CONTROLBAM		 
  
  
             Control file used to call MACS2. 
						 The control file can be a single file, or can be a collection of files. 
						 Even it may not be specified at all, in which case 
						 MACS2 operates without any control.
					 	 Control file can be either in BAM or in tagalign (.gz) format.
					 	 If a set of control files are provided, user needs to 
					 	 ensure that all of the files follow the same format.
             
	-w 	BigWigGenome	 The reference genome which is used to convert BAM file to a BigWig file.
						 If -g option is enabled (i.e. the Bowtie2 index genome is provided) 
						 then this option is not required
 						 as the reference genome will be derived from the genome 
 						 name provided as the Bowtie2 index.
 						 Otherwise, this option needs to be filled with the 
 						 reference genome string (such as 'hg19', 'mm9', etc.)
    -T  Tagmentation	 If 1, corresponds to the adjustment due to Tagmentation. 
						 It shifts the forward and reverse strands by the transposon length, 
						 forms a tagAlign file, and then calls peaks using this tagAlign 
						 file. Similar to the ATAC seq. 
						 The procedure is useful for ChIPMentation type of data. Default 0.	
 	-D  DEBUG_TXT		 this flag signifies whether the read count and other 
 						 statistical paramters are computed or not. 
 						 Can be 1 or 0. If 1, the statistics is generated and 
 						 stored in respective files.
	-O 	Overwrite		 this boolean option signifies whether existing 
						 output files would be overwritten (1) or not (0).
						 Default = 0						 

	-q  MAPQ_THR		 Quality value threshold, below which the mapped reads 
						 are to be removed (Default 30)

	-p  PEAKCALLGENOMESIZE genome size parameter for MACS2 peak calling 
							("hs", "mm", "ce", "dm": default is "hs")

	-t  INT              Set number of sorting, Bowtie2 mapping THREADS [Default = 1].
	-m  MAX_MEM          Set max memory of duplication removal [Default = 8G].
	-a  ALIGNVALIDMAX	 Set the number of (max) valid alignments which will be searched [Default = 4]
	-l  MAXFRAGLEN 		 Set the maximum fragment length to be used for Bowtie2 alignment [Default = 2000]















