import argparse
import pandas as pd
from scripts.find_coverage import find_coverage
# from scripts.find_num_tel import find_num_tel
from scripts.estimate_lengths import estimate_lengths


parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_cells', type=str,
	help='CSV file that stores `BAM_path`, `cov` (coverage, optional), and num_tel (number of telomeres, optional)')
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
		df = find_coverage(df)

	# if 'num_tel' not in df.columns:
	# 	df = find_num_tel(df)

	estimate_lengths(df, args.output)


if __name__ == "__main__":
	main()