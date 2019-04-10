#!/usr/bin/env python

import os
import glob
import argparse

# Set up command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('-e', '--epitope_dir', type=str, required=True,
					help='path to directory containing neoepiscope results'
					)
parser.add_argument('-o', '--output_dir', type=str, required=True,
					help='path to output directory'
					)
args = parser.parse_args()

# Grab phased epitope files
epitope_file_list = glob.glob(os.path.join(args.epitope_dir, '.neoepiscope.out'))

# Open output file and write header
with open(os.path.join(args.output_dir, 'phasing_epitope_data.tsv'), 'w') as o:
	o.write('\t'.join(['Patient', 'Tumor', 'Total_eps', 'Shared', 'Phased_only', 'Unphased_only']) + '\n')
	for phased_eps in epitope_file_list:
		# Get corresponding unphased epitope file
		unphased_eps = phased_eps.replace('.out', '.somatic_unphased.out')
		# Process phased epitopes
		with open(phased_eps) as f:
			f.readline()
			f.readline()
			for ep_line in f:
				ep_tokens = ep_line.strip().split('\t')
				phased.add(ep_tokens[0])
		# Process unphased epitopes
		with open(unphased_eps) as f:
			f.readline()
			f.readline()
			for ep_line in f:
				ep_tokens = ep_line.strip().split('\t')
				if ep_tokens[0] not in phased:
					# Epitope is unique to unphased calling
					unphased_only += 1
				else:
					# Epitope is common to both phased and unphased calling
					shared += 1 
		# Write output
		o.write('\t'.join([tokens[0], tokens[2], str(len(phased) + unphased_only), str(shared), str(len(phased) - shared), str(unphased_only)]) + '\n')
