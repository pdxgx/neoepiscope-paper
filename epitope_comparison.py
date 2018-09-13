#!/usr/bin/env python

import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-d', '--directory', type=str, required=True, 
					help='path to benchmarking directory')

benchmark_dir = args.directory
patients = {'Mel5': 'Mel5_tumor_v_Mel5_normal', 
			'Mel8': 'Mel8_tumor_v_Mel8_normal', 
			'Mel12': 'Mel12_tumor_v_Mel12_normal'}

peptides = {}

for pat in patients:
	peptides[pat] = {}
	# Parse pVAC-Seq data
	pvac_file = ''.join([benchmark_dir, patients[pat], '.final.tsv'])
	with open(pvac_file, 'r') as f:
		f.readline()
		for line in f:
			tokens = line.strip('\n').split('\t')
			neoepitope = tokens[15]
			if neoepitope not in peptides[pat]:
				peptides[pat][neoepitope] = [1, 0, 0, 0, 0]
	# Parse MuPeXI data
	mupexi_file = ''.join([benchmark_dir, patients[pat], '.mupexi'])
	with open(mupexi_file, 'r') as f:
		for i in range(0,6):
			f.readline()
		for line in f:
			tokens = line.strip('\n').split('\t')
			neoepitope = tokens[7]
			if neoepitope not in peptides[pat]:
				peptides[pat][neoepitope] = [0, 1, 0, 0, 0]
			else:
				peptides[pat][neoepitope][1] = 1
	# Parse neoepiscope tumor+germline data
	neoepiscope_file = ''.join([benchmark_dir, patients[pat], 
								'.neoepiscope.out'])
	with open(neoepiscope_file, 'r') as f:
		f.readline()
		for line in f:
			tokens = line.strip('\n').split('\t')
			neoepitope = tokens[0]
			if neoepitope not in peptides[pat]:
				peptides[pat][neoepitope] = [0, 0, 1, 0, 0]
			else:
				peptides[pat][neoepitope][2] = 1
	# Parse neoepiscope tumor-only data
	neoepiscope_file = ''.join([benchmark_dir, patients[pat], 
								'.neoepiscope.tumor.out'])
	with open(neoepiscope_file, 'r') as f:
		f.readline()
		for line in f:
			tokens = line.strip('\n').split('\t')
			neoepitope = tokens[0]
			if neoepitope not in peptides[pat]:
				peptides[pat][neoepitope] = [0, 0, 0, 1, 0]
			else:
				peptides[pat][neoepitope][3] = 1
	# Parse neoepiscope NMD data
	neoepiscope_file = ''.join([benchmark_dir, patients[pat], 
								'.neoepiscope.comprehensive.out'])
	with open(neoepiscope_file, 'r') as f:
		f.readline()
		for line in f:
			tokens = line.strip('\n').split('\t')
			neoepitope = tokens[0]
			if neoepitope not in peptides[pat]:
				peptides[pat][neoepitope] = [0, 0, 0, 0, 1]
			else:
				peptides[pat][neoepitope][4] = 1
	# Write output
	output_file = ''.join([benchmark_dir, pat, '.peptide_overlap.out'])
	with open(output_file, 'w') as f:
		f.write('\t'.join(['Peptide', 'pVAC_Seq', 'MuPeXI', 'neoepiscope_all', 
						   'neoepiscope_tumor', 
						   'neoepiscope_comprehensive']) + '\n')
		for pep in peptides[pat]:
			outline = [pep]
			for i in range(0,5):
				outline.append(str(peptides[pat][pep][i]).upper())
			f.write('\t'.join(outline) + '\n')

output_file = ''.join([benchmark_dir, 'combined.peptide_overlap.out'])
with open(other_out, 'w') as f:
	f.write('\t'.join(['Peptide', 'pVAC_Seq', 'MuPeXI', 'neoepiscope_all', 
					   'neoepiscope_tumor', 
					   'neoepiscope_comprehensive']) + '\n')
	combined_peps = {}
	for pat in peptides:
		for pep in peptides[pat]:
			if pep not in combined_peps:
				combined_peps[pep] = [0, 0, 0, 0, 0]
			for i in range(0, 5):
				combined_peps[pep][i] += peptides[pat][pep][i]
	for peptide in combined_peps:
		outline = [peptide]
		for i in range(0,5):
			outline.append(str(int(bool(combined_peps[peptide][i]))))
		f.write('\t'.join(outline) + '\n')
