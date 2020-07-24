import pandas as pd
import os

def get_number_of_telomeres(df, cn):
	""" Estimate the number of telomeres for each cell in df using the hmmcopy information in cn. """
	# add column for number of telomeres
	df['num_tel'] = None
	for key, row in df.iterrows():
		cell_id = os.path.basename(row['BAM_path']).split('.')[0]
		cell_df = cn.loc[cell_id]
		num_tel = 0
		for chrom, group in cell_df.groupby('chr'):
			temp = 2 * group['state'].mode().values[0]
			num_tel += temp  # add chromosome arm's ploidy to running total
		df.loc[key, 'num_tel'] = num_tel

	return df


def find_num_tel(df, hmm_file_paths):
	""" Load in hmmcopy results and estimate the number of telomeres. """
	# make sure we have proper hmmcopy files if number of telomeres isn't already provided
	try:
		hmm_files = []
		with open(hmm_file_paths, 'rb') as f:
			for line in f:
				hmm_files.append(line.strip())
	except:
		raise ValueError("Error when loading in hmmcopy files")

	# read hmmfiles into dataframe
	cn_pieces = []
	for f in hmm_files:
		piece = pd.read_csv(
			f, index_col=['chr', 'start', 'end', 'width', 'cell_id']
		)
		piece = piece[['reads', 'gc', 'map', 'copy', 'state']]
		cn_pieces.append(piece)
	cn = pd.concat(cn_pieces)
	cn = cn.reset_index()
	cn.set_index('cell_id', inplace=True)

	# use hmmcopy info to estimate number of telomeres for each cell in df
	df = get_number_of_telomeres(df, cn)

	return df
