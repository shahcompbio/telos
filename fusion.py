import os
import subprocess
from argparse import ArgumentParser


def get_args():
	p = ArgumentParser()

	p.add_argument('-i', '--input_bams', type=str, help='txt file containing a list of BAM files (one per line)')
	p.add_argument('-o', '--outbam_dir', type=str, help='directory to store output telbams')
	return p.parse_args()


def txt_to_list(input_bams):
	with open(input_bams, 'r') as f:
		out = []
		for line in f:
			entry = line.strip()
			if ".bam" in entry:
				out.append(entry)
	return out


if __name__ == '__main__':
	argv = get_args()

	input_paths = txt_to_list(argv.input_bams)

	if os.path.exists("out.txt"):
		os.remove("out.txt")

	if os.path.exists("err.txt"):
		os.remove("err.txt")

	# remove slash if included as last character
	if argv.outbam_dir[-1] == '/':
		outbam_dir = argv.outbam_dir[:-1]

	# make sure output directory is present with proper permissions
	if not os.path.exists(outbam_dir):
		try:
			os.makedirs(outbam_dir)
		except:
			raise ValueError('Error: can not find outbam_dir path and could not make it either')
	if not os.access(outbam_dir, os.W_OK | os.X_OK):
		raise ValueError('Error: do not have right permission to write into outbam_dir path')

	# submit each bam's job to the cluster, effectively allowing them to run in parallel on different nodes
	for bam in input_paths:
		subprocess.check_call("bsub -W {time} -M{mem} -R\"span[hosts=1] select[mem>{mem}] rusage[mem={mem}]\" \
			-o out.txt -e err.txt \
			python scripts/fusion_finder.py -i {bam} -o {out}".
			format(bam=bam, out=outbam_dir, time="1:00", mem="8"),
			shell=True)

