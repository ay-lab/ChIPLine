# ChIPLine - ChIP-seq analysis pipeline
---------------------------------------


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

