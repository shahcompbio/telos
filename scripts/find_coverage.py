import os
import pysam
import subprocess
import pandas as pd
import numpy as np


def find_coverage(df):
	print "finding coverage of BAM files"

	df['cov'] = None

	for index, row in df.iterrows():
		bam = row['BAM_path']
		temp_cov_output = str(bam) + ".cov" 
		
		stats = pysam.idxstats(bam)
		nreads_mapped = 0
		nreads_unmapped = 0
		total_bp = 0
		for row in stats.split("\n"):
			row.rstrip("\r")
			fields = row.split("\t")
			if len(fields) > 3 and fields[0] != '*':
				total_bp += int(fields[1])
				nreads_mapped += int(fields[2])
				nreads_unmapped += int(fields[3])

		# print "mapped:", nreads_mapped
		# print "unmapped:", nreads_unmapped
		# print "genome length", total_bp

		cov = (nreads_unmapped + nreads_mapped) * 150. / total_bp

		# print "coverage", cov

		df.loc[index, 'cov'] = cov

	return df
