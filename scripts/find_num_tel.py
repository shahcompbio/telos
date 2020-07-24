import pandas as pd
import os

def find_num_tel(df, cn):
	""" Estimate the number of telomeres for each cell in df using the hmmcopy information in cn. """
	# add column for number of telomeres
	df['num_tel'] = None

	for key, row in df.iterrows():
		cell_id = os.path.basename(row['BAM_path']).split('.')[0]
		cell_df = df.loc[cell_id]
		num_tel = 0
		for chrom, group in cell_df.groupby('chr'):
			num_tel += 2 * mode(group['state'])  # add chromosome arm's ploidy to running total
		df.loc[key, 'num_tel'] = num_tel

	return df
