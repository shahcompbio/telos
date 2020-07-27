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

	if not os.path.exists(argv.outbam_dir):
		try:
			os.makedirs(argv.outbam_dir)
		except:
			# raise ValueError(f'Error: can not find outbam_dir path and could not make it either: \'{outbam_dir}\'.\n')
			raise ValueError('Error: can not find outbam_dir path and could not make it either')
	if not os.access(argv.outbam_dir, os.W_OK | os.X_OK):
		# raise ValueError(f'Error: do not have right permission to write into outbam_dir path: \'{outbam_dir}\'.\n')
		raise ValueError('Error: do not have right permission to write into outbam_dir path')

	# TODO: figure out why process communication errors occur sporadically
	for bam in input_paths:
		subprocess.call("bsub -W {time} -M{mem} -R\"span[hosts=1] select[mem>{mem}] rusage[mem={mem}]\" \
			-o out.txt -e err.txt \
			python scripts/fusion_finder.py -i {bam} -o {out}".
			format(bam=bam, out=argv.outbam_dir, time="6:00", mem="8"),
			shell=True)

