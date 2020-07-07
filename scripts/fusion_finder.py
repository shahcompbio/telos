"""
Subset BAM files to only include reads with evidence of direct telomere fusions.
"""

import sys
import os
import gc
import parabam
from argparse import ArgumentParser

def get_args():
	p = ArgumentParser()

	p.add_argument('-i', '--input_bams', type=str, help='txt file containing a list of BAM files (one per line)')
	p.add_argument('-o', '--outbam_dir', type=str, help='directory to store output telbams')
	return p.parse_args()


def txt_to_list(input_bams):
	with open(input_bams, 'r') as f:
		out = [line.strip() for line in f]
	return out


def rule(reads, constants, master):
	tel_pats = constants["tel_pats"]
	telomere_status = False

	for read in iter(reads):
		for pattern in tel_pats:
			if read.seq is not None and pattern in read.seq:
				telomere_status = True
				break

	results = []
	if telomere_status:
		results = [("telbam", reads[0]), ("telbam", reads[1])]
	return results


class Bam2Telbam(parabam.core.Interface):
	"""Interface to interact with telomerecat programatically. This interface
	is also called by the `telomerecat` script"""

	def __init__(
		self,
		temp_dir=None,
		task_size=250000,
		total_procs=8,
		reader_n=2,
		verbose=False,
		announce=True,
		cmd_run=False,
	):
		self.temp_dir = temp_dir
		self.task_size = task_size
		self.total_procs = total_procs
		self.reader_n = reader_n
		self.verbose = verbose


	def run(self, input_paths, outbam_dir=None):
		"""The main function for invoking the part of the 
			program which creates a telbam from a bam
		Arguments:
			bams (list): The BAM files we wish to run telomerecat telbam on
			total_procs (int): The maximum numbers of task that will be run at one time
			task_size (int): The amount of reads that any one task will process concurrently
			keep_in_temp (bool): Files will be kept in temp file after processing. 
				Useful for incorporation into pipelines"""


		subset_types = ["telbam"]
		tel_pats = ["TTAGGGTTAGGG", "CCCTAACCCTAA"]

		# need to define my constants and engine here:
		telbam_constants = {"thresh": 1, "tel_pats": tel_pats}

		final_output_paths = {}

		keep_in_temp = outbam_dir is None
		if not keep_in_temp:
			# check if the folder is writable
			if not os.access(outbam_dir, os.W_OK | os.X_OK):
				raise ValueError(f'Error: do not have right permission to write into outbam_dir path: \'{outbam_dir}\'.\n')
			if not os.path.exists(outbam_dir):
				try:
					os.mkdirs(outbam_dir)
				except:
					raise ValueError(f'Error: can not find outbam_dir path and could not make it either: \'{outbam_dir}\'.\n')

		for input_path in input_paths:

			subset_interface = parabam.Subset(
				temp_dir=self.temp_dir,
				total_procs=self.total_procs,
				task_size=self.task_size,
				reader_n=self.reader_n,
				verbose=self.verbose,
				pair_process=True,
				include_duplicates=True,
				keep_in_temp=keep_in_temp,
			)

			# call to parabam subset
			telbam_paths = subset_interface.run(
				input_paths=[input_path],
				subsets=subset_types,
				constants=telbam_constants,
				rule=rule,
				outbam_dir=outbam_dir
			)

			gc.collect()
			final_output_paths.update(telbam_paths)

		return final_output_paths


if __name__ == '__main__':
	argv = get_args()

	input_paths = txt_to_list(argv.input_bams)

	telbam = Bam2Telbam()
	telbam.run(input_paths, argv.outbam_dir)

