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

#=================
# main executable script of the ChIP seq pipeline
#=================
# developed by - Sourya Bhattacharyya
# Vijay-AY lab
# La Jolla Institute for Allergy and Immunology
#=================

CodeExec=`pwd`'/pipeline.sh'

inpfile='/mnt/NGSAnalyses/ChIP-Seq/Mapping/004997_Log7_DY_R24_NK_ChIP_Merged/All_Uniquely_Mapped_BAM_Files/016_A_NK_16_GM12878_CTRL_Final.bam'
outdir='/mnt/BioAdHoc/Groups/vd-vijay/sourya/ChIP_seq/VIJAY_LAB_DATASETS/004997_Log7_DY_R24_NK_ChIP_Merged/016_A_NK_16_GM12878_CTRL_Final_Copy/'

# refgenome='/mnt/BioHome/ferhatay/data/bowtie2_indexes/hg19-full/hg19'
# $CodeExec -C 'configfile' -f $inpfile -n '001_RunA_SLE10_S1057a_1_Final' -g $refgenome -d $outdir -t 8 -m '16G' -q 30 -T 1 -D 1 -O 1
$CodeExec -C 'configfile' -f $inpfile -n '016_A_NK_16_GM12878_CTRL' -d $outdir -t 8 -m '16G' -q 30 -T 1 -D 0 -O 0 -w "hg38"




