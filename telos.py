import argparse
import pandas as pd
from scripts.find_coverage import find_coverage
from scripts.find_num_tel import find_num_tel
from scripts.estimate_lengths import estimate_lengths


parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_cells', type=str,
	help='CSV file that stores `BAM_path`, `cov` (coverage, optional), and num_tel (number of telomeres, optional)')
parser.add_argument('-hmm', '--hmmcopy', type=str, nargs='?', default=None,
	help='paths to hmmcopy output for the libraries containing the input_cells')
parser.add_argument('-o', '--output', type=str,
	help='path for writing the output csv file')


def load_input(input_cells_path):
	return pd.read_csv(input_cells_path, index_col=False)


def main():
	args = parser.parse_args()

	df = load_input(args.input_cells)
	print df.head()

	# TODO: should I swtich this to a snakemake workflow?
	# not sure if I could execute functions/rules dependent on columns of input file (like I'm doing here) if I switch to snakemake
	assert 'BAM_path' in df.columns

	if 'cov' not in df.columns:
		print 'coverage not provided'
		df = find_coverage(df)
		print 'done estimating coverage'

	if 'num_tel' not in df.columns:
		if args.hmmcopy is None:
			# TODO: run the single cell pipeline to get hmmcopy output with all files stored in a list much like args.hmmcopy would be
			# hmmcopy_paths = run_single_cell_pipeline()
			# df = find_num_tel(df, hmmcopy_paths)
			pass
		else:
			df = find_num_tel(df, args.hmmcopy, '/juno/work/shah/funnellt/projects/rdc-mutsig/data/ref/gap_hg19.txt')
		print df

	estimate_lengths(df, args.output)


if __name__ == "__main__":
	main()