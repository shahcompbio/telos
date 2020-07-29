import pandas as pd
import numpy as np
import os


def find_num_tel(df, hmm_file_paths, gaps_path=None):
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
	cn = cn.query('gc > 0')
	cn = cn.query('map > 0.75')

	# find ploidy for each cell using hmmcopy info
	summary_func = 'median'
	if gaps_path is not None:
		gaps = pd.read_csv(argv.gaps, sep='\t', dtype={'chrom': str})

		centromeres = gaps.loc[
			gaps['type'] == 'centromere', ['chrom', 'chromStart', 'chromEnd']
		]

		cn = find_chr_arm_ploidy(
			cn, summary_func, centromeres
		)
	else:
		cn = find_chr_arm_ploidy(
			cn, summary_func
		)

	# use ploidy info to estimate number of telomeres for each cell in df
	cn.set_index('cell_id', inplace=True)
	df = get_number_of_telomeres(df, cn, gaps_path)

	return df


def get_number_of_telomeres(df, cn, gaps_path=None):
	""" Estimate the number of telomeres for each cell in df using the ploidy information in cn. """
	# add column for number of telomeres
	df['num_tel'] = None
	for key, row in df.iterrows():
		cell_id = os.path.basename(row['BAM_path']).split('.')[0]
		cell_df = cn.loc[cell_id]
		num_tel = 0
		# group by chromosome if no centromere gaps provided
		if gaps_path is None:
			for chrom, group in cell_df.groupby('chr'):
				temp = 2 * group['ploidy'].mode().values[0]  # should all be same value so taking the mode
				num_tel += temp  # add chromosome arm's ploidy to running total
		else:
			for chrom, group in cell_df.groupby(['chr', 'arm']):
				temp = group['ploidy'].mode().values[0]
				num_tel += temp  # add chromosome arm's ploidy to running total

		df.loc[key, 'num_tel'] = num_tel

	return df


def find_chr_arm_ploidy(cn, summary_func, centromeres=None):
	print 'Removing Y chromosome'
	cn = cn[cn['chr'] != 'Y'].copy()

	first_mode = lambda x: mode(x)[0][0]
	if centromeres is not None:
		cn = assign_chr_arm(cn, centromeres)
		cn_grouped = cn.groupby(['cell_id', 'chr', 'arm'])
	else:
		cn_grouped = cn.groupby(['cell_id', 'chr'])

	if summary_func == 'mode':
		cn['ploidy'] = cn_grouped['state'].transform(first_mode)
	elif summary_func == 'median':
		cn['ploidy'] = cn_grouped['state'].transform(np.median)
	else:
		template = 'Warning: bad summary function choice "{}", ' + \
			'resorting to median'
		print(template.format(summary_func))
		cn['ploidy'] = cn_grouped['state'].transform(np.median)

	return cn


def assign_chr_arm(cn, centromeres):
	cn['arm'] = ''

	for c in cn['chr'].unique():
		centromere = centromeres[(centromeres['chrom'] == 'chr' + c)]

		chr_mask = cn['chr'] == c
		p_mask = cn['start'] < centromere['chromStart'].iloc[0]
		q_mask = cn['end'] > centromere['chromEnd'].iloc[0]

		# check we don't still have centromere bins
		n_centromere_bins = sum(chr_mask & ~(p_mask | q_mask))
		if (n_centromere_bins > 0):
			template = 'WARNING: {} centromere bins found in chr {}, removing'
			print(template.format(n_centromere_bins, c))

		cn.loc[chr_mask & p_mask, 'arm'] = 'p'
		cn.loc[chr_mask & q_mask, 'arm'] = 'q'

	cn = cn[cn['arm'].isin(['p', 'q'])]

	return cn

