#!/bin/bash -ex
#PBS -l nodes=1:ppn=1
#PBS -l mem=10GB
#PBS -l walltime=24:00:00
#PBS -m ae
#PBS -j eo
#PBS -V
source ~/.bashrc
#source ~/.bash_profile
hostname
TMPDIR=/scratch
cd $PBS_O_WORKDIR

#=================================
# sample script for IDR execution
# where peaks generated from multiple ChIP-seq replicates are provided as input
#=================================

# main executable of IDR script
# when peak files are used as input
IDRScript='./IDR_Codes/IDRMain.sh'

# main executable of IDR script
# when BAM files are used as input
IDRScriptBAM='./IDR_Codes/IDR_SubSampleBAM_Main.sh'

#******************************
# path containing the IDRCode package by Anshul Kundaje et. al.
# user should replace this path with their custom installation directory
IDRCodePackage='/home/sourya/packages/idrCode/'
#******************************


#====================
# IDR testing script 1
# examining IDR between two peak files
# top 25K common peaks between two samples are experimented
#====================

SampleBaseDir='/mnt/BioAdHoc/Groups/vd-vijay/sourya/ChIP_seq/VIJAY_LAB_DATASETS/004995_Log7_DY_R24_CD4N_ChIP_edited_Merged/'

$IDRScript -I $SampleBaseDir'096_F_CD4N_16_GM12878_CTRL_Final/MACS2_Default_Tag_No_Control/096_F_CD4N_16_GM12878_CTRL_Final.macs2_peaks.narrowPeak_Q0.01filt' -I $SampleBaseDir'080_E_CD4N_16_GM12879_CTRL_Final/MACS2_Default_Tag_No_Control/080_E_CD4N_16_GM12879_CTRL_Final.macs2_peaks.narrowPeak_Q0.01filt' -d '/mnt/BioAdHoc/Groups/vd-vijay/sourya/ChIP_seq/VIJAY_LAB_DATASETS/004995_Log7_DY_R24_CD4N_ChIP_edited_Merged_IDR/Sample_IDR_Peaks' -P $IDRCodePackage


#====================
# IDR testing script 2
# examining IDR between two BAM files
# first these BAM files are subsampled
# and their peaks are estimated using MACS2
# top 25K common peaks between two samples are experimented
# no control BAM file is provided
# user may specify one or more control BAM files using -C option
# like -C control1.bam -C control2.bam etc.
# tagmentation parameter (-T) is provided as 1
# user may reset it to 0, depending on the ChIP sample
#====================

SampleBaseDir='/mnt/BioAdHoc/Groups/vd-vijay/sourya/ChIP_seq/VIJAY_LAB_DATASETS/004995_Log7_DY_R24_CD4N_ChIP_edited_Merged/'

$IDRScriptBAM -I $SampleBaseDir'096_F_CD4N_16_GM12878_CTRL_Final/Alignment_MAPQ30/096_F_CD4N_16_GM12878_CTRL_Final.align.sort.MAPQ30.rmdup.bam' -I $SampleBaseDir'080_E_CD4N_16_GM12879_CTRL_Final/Alignment_MAPQ30/080_E_CD4N_16_GM12879_CTRL_Final.align.sort.MAPQ30.rmdup.bam' -d '/mnt/BioAdHoc/Groups/vd-vijay/sourya/ChIP_seq/VIJAY_LAB_DATASETS/004995_Log7_DY_R24_CD4N_ChIP_edited_Merged_IDR/Sample_IDR_BAMFiles' -P $IDRCodePackage -c 25000 -T 1








