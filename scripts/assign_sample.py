#!/usr/bin/env python

import glob 
import os
import subprocess
import sys

prefix = 'pilot_run'

os.mkdir('sample_fastq')
bam_list = glob.glob('*.bam')
barcode_map = open(sys.argv[1], 'rU')
barcode_dict = {}
for line in barcode_map:
	line = line.strip('\n')
	key = line.split('\t')[0] + '--' + line.split('\t')[1]
	barcode_dict[key] = line.split('\t')[2]

for key in barcode_dict:
	bam_in_name = prefix + '_ccs_p5_rq0.999_demux.' + str(key) + '.bam'
	fastq_out_name = 'sample_fastq' + '/' + barcode_dict[key] + '_' + prefix
	bam2fastq_cmd = 'bam2fastq -o %s %s' % (fastq_out_name, bam_in_name) 
	process = subprocess.Popen(bam2fastq_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)	
	(out, err) = process.communicate() #the stdout and stderr








			