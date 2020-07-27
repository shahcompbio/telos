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

	# TODO: this loop should be submitting an LSF job for each bam instead of waiting for run() to finish before moving onto next bam
	for bam in input_paths:
		print bam
		subprocess.call("bsub -W {time} -M{mem} -R\"span[hosts=1] select[mem>{mem}] rusage[mem={mem}]\" -o out.txt -e err.txt \
			python scripts/fusion_finder.py -i {bam} -o {out}".
			format(bam=bam, out=argv.outbam_dir, time="6:00", mem="24"),
			shell=True)

