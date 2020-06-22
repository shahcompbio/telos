import pysam
import os
import subprocess
import pandas as pd


def find_coverage(df):
	df['cov'] = None

	# TODO: use multiprocessing to parallelize this for-loop
	for index, row in df.iterrows():
		bam = row['BAM_path']
		temp_cov_output = str(bam) + ".cov" 
		
		subprocess.check_call(
			"module load samtools ; \
			samtools depth -a {bam}  | awk '{sum+=$3} END { print sum/NR}' > {out}".
			format(bam=bam, out=temp_cov_output),
			shell=True)

		with open(str(temp_cov_output)) as f:
			cov = f.readline().strip()

		os.remove(temp_cov_output)

		df.loc[index, 'cov'] = cov

