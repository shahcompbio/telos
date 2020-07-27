"""
Subset BAM files to only include reads with evidence of direct telomere fusions.
"""
import random
import sys
import os
import gc
import parabam
import time
import shutil
from argparse import ArgumentParser

def get_args():
	p = ArgumentParser()

	p.add_argument('-i', '--input', type=str, help='BAM file input')
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




def run( 
		input_paths, 
		outbam_dir=None, 
		temp_dir=None, 
		task_size=250000, 
		total_procs=8,
		reader_n=2,
		verbose=False,
	):
	"""The main function for invoking the part of the 
		program which creates a telbam from a bam
	Arguments:
		input_paths (list): The BAM files we wish look for patterns in
		total_procs (int): The maximum numbers of task that will be run at one time
		task_size (int): The amount of reads that any one task will process concurrently
		keep_in_temp (bool): Files will be kept in temp file after processing. 
			Useful for incorporation into pipelines"""

	if temp_dir is None:
		temp_dir = "telos_temp_" + str(time.time()) + str(random.randint(0, 99999))
		temp_dir = temp_dir.replace('.', '_')

	print 'temp_dir', temp_dir
	if not os.path.exists(temp_dir):
		try:
			os.makedirs(temp_dir)
		except:
			# raise ValueError(f'Error: can not find temp_dir path and could not make it either: \'{temp_dir}\'.\n')
			raise ValueError('Error: can not find temp_dir path and could not make it either')

	subset_types = ["telbam"]
	tel_pats = ["TTAGGGTTAGGG", "CCCTAACCCTAA"]

	# need to define my constants and engine here:
	telbam_constants = {"thresh": 1, "tel_pats": tel_pats}

	final_output_paths = {}

	keep_in_temp = outbam_dir is None
	if not keep_in_temp:
		# check if the folder is writable
		if not os.path.exists(outbam_dir):
			try:
				os.makedirs(outbam_dir)
			except:
				# raise ValueError(f'Error: can not find outbam_dir path and could not make it either: \'{outbam_dir}\'.\n')
				raise ValueError('Error: can not find outbam_dir path and could not make it either')
		if not os.access(outbam_dir, os.W_OK | os.X_OK):
			# raise ValueError(f'Error: do not have right permission to write into outbam_dir path: \'{outbam_dir}\'.\n')
			raise ValueError('Error: do not have right permission to write into outbam_dir path')

	for input_path in input_paths:

		subset_interface = parabam.Subset(
			temp_dir=temp_dir,
			total_procs=total_procs,
			task_size=task_size,
			reader_n=reader_n,
			verbose=verbose,
			pair_process=True,
			include_duplicates=True,
			keep_in_temp=keep_in_temp,
		)

		print 'input_path', input_path
		print 'subset_types', subset_types
		print 'telbam_constants', telbam_constants
		print 'outbam_dir', outbam_dir
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

		for k, v in telbam_paths.iteritems():
			for key, item in v.iteritems():
				temp_output_path = item

		print 'telbam_paths', telbam_paths
		print 'temp_output_path', temp_output_path

	basename = os.path.basename(temp_output_path)
	new_basename = basename.replace('_telbam', '_fusions')
	final_output = str(outbam_dir) + '/' + str(new_basename)

	shutil.move(temp_output_path, final_output)

	shutil.rmtree(temp_dir, ignore_errors=True)

	return final_output_paths


if __name__ == '__main__':
	argv = get_args()

	run([argv.input], argv.outbam_dir)

