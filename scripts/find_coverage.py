import os
import pysam
import subprocess
import pandas as pd
import numpy as np


def find_coverage(df):
	print "finding coverage of BAM files"

	df['cov'] = None

	# TODO: use multiprocessing to parallelize this for-loop
	for index, row in df.iterrows():
		bam = row['BAM_path']
		temp_cov_output = str(bam) + ".cov" 
		
		# depths = pysam.depth(bam, "-a")
		# depths = [i for i in depths]
		# print len(depths)
		# print depths[:20]
		# cov = np.mean(depths)

		# num_reads = reduce(lambda x, y: x + y, [eval('+'.join(l.rstrip('\n').split('\t')[2:])) for l in pysam.idxstats(bam)])
		# print num_reads
		# depth = num_reads / 3E9
		# print depth

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

		print "mapped:", nreads_mapped
		print "unmapped:", nreads_unmapped
		print "genome length", total_bp

		cov = (nreads_unmapped + nreads_mapped) * 150. / total_bp

		print "coverage", cov

		# subprocess.check_call(
		# 	"module load samtools ; \
		# 	samtools idxstats {bam} > {out}".
		# 	format(bam=bam, out=temp_cov_output),
		# 	shell=True)

		# subprocess.check_call(
		# 	"module load samtools ; \
		# 	samtools depth -a {bam}  | awk \'{ sum+=$3 } END { print sum/NR } \' > {out}".
		# 	format(bam=bam, out=temp_cov_output),
		# 	shell=True)

		# with open(str(temp_cov_output)) as f:
		# 	num_reads = f.readline().strip()

		# print(num_reads)

		# cov = float(num_reads) / 3E9

		# os.remove(temp_cov_output)

		df.loc[index, 'cov'] = cov

	return df
