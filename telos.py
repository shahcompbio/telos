import argparse
import pandas as pd
from scripts.find_coverage import find_coverage
from scripts.find_num_tel import find_num_tel
from scripts.estimate_lengths import estimate_lengths


parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_cells', type=str,
	help='CSV file that stores `BAM_path`, `cov` (coverage, optional), and num_tel (number of telomeres, optional)')
parser.add_argument('-h', '--hmmcopy', type=str, nargs='+',
	help='paths to hmmcopy output for the libraries containing the input_cells')
parser.add_argument('-o', '--output', type=str,
	help='path for writing the output csv file')


def load_input(input_cells_path):
	return pd.read_csv(input_cells_path, index_col=False)


def main():
	args = parser.parse_args()

	df = load_input(args.input_cells)
	print df.head()

	assert 'BAM_path' in df.columns

	if 'cov' not in df.columns:
		print 'coverage not provided'
		df = find_coverage(df)
		print 'done estimating coverage'

	if 'num_tel' not in df.columns:
		# make sure we have proper hmmcopy files if number of telomeres isn't already provided
		try:
			cn_pieces = []
			for f in args.hmmcopy:
				piece = pd.read_csv(
					f, index_col=['chr', 'start', 'end', 'width', 'cell_id']
				)
				piece = piece[['reads', 'gc', 'map', 'copy', 'state']]
				cn_pieces.append(piece)
			cn = pd.concat(cn_pieces)
			cn = cn.reset_index()
		except:
			raise ValueError("Error when loading in hmmcopy files")

		# use hmmcopy info to estimate number of telomeres for each cell in df
		df = find_num_tel(df, cn)
		print df

	estimate_lengths(df, args.output)


if __name__ == "__main__":
	main()