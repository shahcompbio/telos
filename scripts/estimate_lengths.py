import os
import pysam
import subprocess
import csv
import pandas as pd
from collections import OrderedDict


def get_telbam(path_to_bam):
	# telbam file created by telomerecat
	temp_telbam = os.path.splitext(os.path.basename(path_to_bam))[0] + "_telbam.bam"
	output_telbam = path_to_bam.replace(".bam", "_telbam.bam")

	# create appropriate telbam file using telomerecat
	subprocess.check_call(
		"telomerecat bam2telbam {input}".
		format(input=path_to_bam),
		shell=True)
	subprocess.check_call(
		"mv {temp_telbam} {output}".
		format(temp_telbam=temp_telbam, output=output_telbam),
		shell=True)

	return output_telbam


def subset_telbams(temp_telbam_list, temp_good_telbams, temp_bad_telbams, k=5):
	good_telbams = []
	bad_telbams = []
	with open(temp_telbam_list, 'rb') as f:
		for line in f:
			telbam_file = line.strip()
			try:
				file = pysam.AlignmentFile(telbam_file, "rb")
				count = file.count(until_eof=True)  # find number of reads in telbam
				if count < k:
					bad_telbams.append(telbam_file)
				else:
					good_telbams.append(telbam_file)
			except:
				bad_telbams.append(telbam_file)
		f.close()

	# write output for good telbams
	with open(temp_good_telbams, 'w') as f:
		for item in good_telbams:
			f.write(item + '\n')
		f.close()

	# write output for good telbams
	with open(temp_bad_telbams, 'w') as f:
		for item in bad_telbams:
			f.write(item + '\n')
		f.close()


def make_pseudobulk(temp_good_telbams, temp_pseudobulk):
	subprocess.check_call(
		"module load samtools ; samtools merge {output} -b {input}".
		format(input=temp_good_telbams, output=temp_pseudobulk),
		shell=True)
	# throw a process error if the output file doens't exist yet
	subprocess.check_call(
		"ls {output}".
		format(output=temp_pseudobulk),
		shell=True)


def add_cov_ntel_to_good_telbams(temp_good_telbams, new_good_telbams, df):
	df = df.set_index('telbam_name')
	with open(new_good_telbams, 'w') as n:
		with open(temp_good_telbams, 'rb') as f:
			for item in f:
				item = item.strip()
				cov = df.loc[item, 'cov']
				ntel = df.loc[item, 'num_tel']
				n.write(str(item) + ',' + str(cov) + ',' + str(ntel) + '\n')
			f.close()
		n.close()


def get_lengths(temp_good_telbams, temp_bad_telbams, temp_pseudobulk, csv_output):
	subprocess.check_call(
		"telomerecat telbam2length -cnt -e -i {input} -b {bulk} --output {output}".
		format(input=temp_good_telbams, bulk=temp_pseudobulk, output=csv_output),
		shell=True)

	# throw a process error if the output file doesn't exist yet
	connection = subprocess.check_call(
			"ls {output}".
			format(output=csv_output),
			shell=True)

	if connection == 0:
		with open(temp_bad_telbams, 'rb') as bad_files:
			for telbam in bad_files:
				telbam = telbam.strip()
				no_read_dict = OrderedDict([
								("Sample", os.path.splitext(os.path.basename(telbam))[0].split("_")[0] + ".bam"),
								("F1", "NA"),
								("F2", "NA"), 
								("F4", "NA"),
								("Psi", "NA"),
								("Insert_mean", "NA"),
								("Insert_sd", "NA"),
								("Read_length", "NA"),
								("Initial_read_length", "NA"),
								("F2a", "NA"),
								("F2a_c", "NA"),
								("Length", "NA"),
								("Length_std", "NA")])
				with open(csv_output, 'a') as f:  # use 'wb' mode to write a new file
					w = csv.DictWriter(f, no_read_dict.keys())
					w.writerow(no_read_dict)


def estimate_lengths(df, path_to_output):
	df['telbam_name'] = None

	temp_telbam_list = "temp_telbam_list.txt"
	temp_good_telbams = "temp_good_telbams.txt"
	temp_good_telbams_cov_ntel = "temp_good_telbams_cov_ntel.txt"
	temp_bad_telbams = "temp_bad_telbams.txt"
	temp_pseudobulk = "temp_pseudobulk.bam"
	temp_output_csv = "temp_telos_output.csv"

	for index, row in df.iterrows():
		telbam_path = get_telbam(row['BAM_path'])  # make the telbam
		df['telbam_name'][index] = telbam_path

	df["Sample"] = df["telbam_name"].apply(lambda x: os.path.splitext(os.path.basename(x))[0].replace("_telbam", "") + ".bam")

	# make temporary txt file containing one telbam path per line
	df.to_csv(temp_telbam_list, sep='\n', columns=['telbam_name'], index=False, header=False)

	# subset to find the "good" telbams
	subset_telbams(temp_telbam_list, temp_good_telbams, temp_bad_telbams)

	# create a pseudobulk using the good telbams
	make_pseudobulk(temp_good_telbams, temp_pseudobulk)

	# use legend to add coverage and number of telomeres for all the good telbams
	add_cov_ntel_to_good_telbams(temp_good_telbams, temp_good_telbams_cov_ntel, df)

	# estimate telomere lengths, creating temporary output.csv file
	get_lengths(temp_good_telbams_cov_ntel, temp_bad_telbams, temp_pseudobulk, temp_output_csv)

	# merge temporary output.csv with legend dataframe
	telos_df = pd.read_csv(temp_output_csv)
	telos_df.drop(columns=["coverage", "num_tel"], inplace=True)
	df = pd.merge(df, telos_df, on="Sample")
	df.drop(columns=["Sample"], inplace=True)

	# remove temporary files
	os.remove(temp_telbam_list)
	os.remove(temp_bad_telbams)
	os.remove(temp_good_telbams)
	os.remove(temp_good_telbams_cov_ntel)
	os.remove(temp_output_csv)
	os.remove(temp_pseudobulk)

	# save output file containing legend and length estimates
	df.to_csv(path_to_output, index=False)
