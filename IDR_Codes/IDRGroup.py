#!/usr/bin/env python

"""
Created: August 02,2017

This is a special script file
it groups individual samples according to different criterion (such as donor ID)
and then performs IDR on individual groups

Author: Sourya Bhattacharyya
Vijay-Ay lab, LJI
"""

import os
import sys
import numpy

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

#===============================================
def main():

	basedir_list = ['/mnt/BioAdHoc/Groups/vd-vijay/sourya/ChIP_seq/003696_GrSe95_R24_NK_ChIPAB/', '/mnt/BioAdHoc/Groups/vd-vijay/sourya/ChIP_seq/003705_R24_NK_ChIPCH/']

	ngroup = 10

	# list of sample types which need to be excluded from consideration
	Exclude_Sample_List = ['D8', 'GM12878']

	# max no of samples a group can have
	NsampleLimitPerGrp = 12

	# we create a set of empty lists
	# each of these empty lists will contain the samples for IDR operation
	SampleListbyGrp = [[]] * ngroup

	# navigate through individual base directory
	for basedir_idx in range(len(basedir_list)):
		currbasedir = basedir_list[basedir_idx]
		if 0:
			print '\n Examining the currbasedir: ', currbasedir
		# list of folders within a particular directory
		dirs = [d for d in os.listdir(currbasedir) if os.path.isdir(os.path.join(currbasedir, d))]
		no_of_dirs = len(dirs)
		# randomly permute the contents of this list
		permute_dir_idx_list = numpy.random.permutation(no_of_dirs)

		for d_idx in permute_dir_idx_list:
			d = dirs[d_idx]
			name_contents = str(d).split('_')
			if 0:
				print '\n Examine placement of directory: ', d, '  Its list version: ', name_contents
			# if the 6th field denotes that it is within "Exclude_Sample_List", discard
			if (str(name_contents[5]) in Exclude_Sample_List):
				if 0:
					print '\n The directory is not considered for analysis'
				continue
			donor_ID = str(name_contents[6])
			if 0:
				print '\n donor ID for this sample: ', donor_ID

			# permute thje groups for examine
			Grp_Permute_List_idx = numpy.random.permutation(len(SampleListbyGrp))
			for grp_idx in Grp_Permute_List_idx: 	#range(len(SampleListbyGrp)):
				if (len(SampleListbyGrp[grp_idx]) == 0):
					# the current group does not have any sample  - insert the current sample
					temp = []
					temp.append(os.path.join(currbasedir, d))
					SampleListbyGrp[grp_idx] = temp
					if 0:
						print '\n --- test: grp_idx: ', grp_idx, '   length 0 - auto insert'
						print 'SampleListbyGrp: ', SampleListbyGrp
					break
				elif (len(SampleListbyGrp[grp_idx]) < NsampleLimitPerGrp):
					# check the samples in the current group
					# flag 1 means that there is a matching donor in the samples of current group 
					# in such a case, we cannot insert the current sample in this group
					# for flag 0 we can insert the sample
					if 0:
						print '\n --- test: grp_idx: ', grp_idx, '   length: ', len(SampleListbyGrp[grp_idx])
					flag = 0
					for i in range(len(SampleListbyGrp[grp_idx])):
						currsample_fullname = str(SampleListbyGrp[grp_idx][i])
						if 0:
							print '\n i: ', i, '  currsample_fullname: ', currsample_fullname
						OnlySampleName = str(currsample_fullname.rsplit('/',1)[1])
						if 0:
							print '\n OnlySampleName: ', OnlySampleName
						onlysample_contents = OnlySampleName.split('_')
						if 0:
							print '\n onlysample_contents: ', onlysample_contents
						if (len(onlysample_contents) > 6):
							donor_ID_curr_sample = str(onlysample_contents[6])
							if (donor_ID == donor_ID_curr_sample):
								if 0:
									print '\n There exists donor with the same ID - this group cannot contain the current sample'
								flag = 1
								break
					if (flag == 0):
						if 0:
							print '\n --- test: grp_idx: ', grp_idx, '   insert sample in this group'
						temp = []
						temp.append(os.path.join(currbasedir, d))
						SampleListbyGrp[grp_idx].extend(temp)
						break

	# now we print the groups
	if 1:
		for grp_idx in range(len(SampleListbyGrp)):
			print '\n grp_idx: ', grp_idx, '   no of samples: ', len(SampleListbyGrp[grp_idx])
			print '\n samples list: ', SampleListbyGrp[grp_idx][:]

	#======================================
	# IDR analysis for each group
	#======================================
	# base directory containing all the output results
	BASEOUT = '/mnt/BioAdHoc/Groups/vd-vijay/sourya/ChIP_seq/BEN_July25_2017_ChIP_NK_IDR/Group_Specific_Random_IDR/'
	# main executable of IDR script
	IDRScript='/home/sourya/proj/2017_ChIP_Seq/bin/IDRMain.sh'
	MACS2_Tag='MACS2_Tag_No_Control'

	for grp_idx in range(len(SampleListbyGrp)):	
		print '\n\n **** IDR analysis of the group : ', grp_idx, ' **** \n\n'	
		OutDir = BASEOUT + 'Group_' + str(grp_idx)
		print '\n Output directory: ', OutDir, '\n'
		nsample = len(SampleListbyGrp[grp_idx])
		# command for executing IDR
		exec_cmd = IDRScript
		# include the input samples
		for i in range(nsample):
			currsample_fullname = str(SampleListbyGrp[grp_idx][i])
			OnlySampleName = str(currsample_fullname.rsplit('/',1)[1])
			Macs2_file_name1 = str(currsample_fullname) + '/' + str(MACS2_Tag) + '/' + str(OnlySampleName) + '.macs2_peaks.narrowPeak'
			Macs2_file_name2 = str(currsample_fullname) + '/' + str(MACS2_Tag) + '/' + str(OnlySampleName) + '_peaks.narrowPeak'
			if (os.path.isfile(Macs2_file_name1) == True):
				exec_cmd = exec_cmd + ' -I ' + str(Macs2_file_name1)
			else:
				exec_cmd = exec_cmd + ' -I ' + str(Macs2_file_name2)
		# include the directories
		exec_cmd = exec_cmd + ' -d ' + str(OutDir) + '/OrigBAM/MACS2_Tag_Peak_25k/ -c 25000'
		print '\n Execution command: ', exec_cmd
		os.system(exec_cmd)

#===============================================
if __name__ == "__main__":
    main()


