#!/usr/bin/env python

from __future__ import print_function
import glob
import pickle
import random
from collections import defaultdict
from numpy import median

random.seed(8211954)

rna_path = '/home/exacloud/lustre1/CompBio/projs/RThompson/immunotherapy/jobs/rna_support/'

dicts = glob.glob(''.join([rna_path, 'rna_dicts/*.rna_support_updated.pickle']))

file_dict = defaultdict(list)
distance_dict = defaultdict(list)

bad_tumors = ['SRR3653958', 'SRR3653960', 'SRR3653962', 'SRR2672864', 'SRR5149302', 
			  'SRR5149303', 'SRR5149304', 'SRR5149305', 'SRR5149306', 'SRR5149307',
			  'Pem10T', 'Pem17T', 'Pem12T', 'Pem13T', 'Pem22T', 'Pem26T', 'Pem27T', 
			  'Pem28T', 'Pem35T', 'Pem1T', 'SRR6450189', 'SRR6450222', 'SRR6450224',
			  'SRR6450226', 'SRR6450227', 'SRR6450228', 'SRR6450229', 'SRR6450230',
			  'SRR6450231', 'SRR6450232', 'SRR6461147', 'SRR6461150', 'SRR6461155',
			  'SRR6461159', 'SRR6461162', 'SRR6461163', 'SRR6461164', 'SRR6461165',
			  'SRR6461166', 'SRR6461167', 'SRR6461168', 'SRR6461169', 'SRR6461170',
			  'SRR6461171', 'SRR6461172', 'SRR6461173', 'SRR6461175', 'SRR6461177',
			  'SRR6461180', 'SRR6461181', 'SRR6468302', 'SRR6468307', 'SRR6468311',
			  'SRR6468316', 'SRR6468321', 'SRR6468326', 'SRR6468328', 'SRR6468330',
			  'SRR6468332', 'SRR6468336', 'SRR6468339', 'SRR6468341', 'SRR6468347',
			  'SRR6468348', 'SRR6468373', 'SRR6468374', 'SRR6468375', 'SRR6491482',
			  'SRR6504580', 'SRR6504581', 'SRR6504582', 'SRR6504583', 'SRR6504586',
			  'SRR6504588', 'SRR6504589', 'SRR6504838', 'SRR6504840', 'SRR6889242']

for d in dicts:
	patient = '.'.join(d.split('/')[-1].replace('.rna_support_updated.pickle', '').split('.')[0:-1])
	tumor = d.split('/')[-1].replace('.rna_support_updated.pickle', '').split('.')[-1]
	if tumor not in bad_tumors:
		with open(d, 'rb') as f:
			rna_support = pickle.load(f)
		data = [len(rna_support['all_relevant_germline']), len(rna_support['all_relevant_somatic']), 
				len(rna_support['covered_germline']), len(rna_support['covered_somatic']),
				len(rna_support['validated_germline']), len(rna_support['validated_somatic']),
				len(rna_support['phased_germline']), len(rna_support['phased_somatic']), len(rna_support['phased_combined']), 
				len(rna_support['total_germline_pairs']), len(rna_support['total_somatic_pairs']),
				len(rna_support['covered_germline_pair']), len(rna_support['covered_somatic_pair']),
				len(rna_support['supported_germline']), len(rna_support['supported_somatic']),
				len(rna_support['not_supported_germline']), len(rna_support['not_supported_somatic']),
				len(rna_support['novel_germline']), len(rna_support['novel_somatic']),
				len(rna_support['covered_across_exon_pairs']), len(rna_support['validated_across_exon_pairs']),
				len(rna_support['novel_across_exon_pairs']), len(rna_support['not_supported_across_exon_pairs'])]
				len(rna_support['predicted_across_exon_pairs'])
		file_dict[(patient, tumor)] = [float(x) for x in data]
		if (len(rna_support['novel_germline']) > 0) or (len(rna_support['novel_somatic']) > 0):
			dis_path = ''.join([rna_path, 'rna_dicts/', patient, '.', tumor, '.pair_distances.pickle'])
			with open(dis_path, 'rb') as f:
				distances = pickle.load(f)
			for pair in rna_support['novel_germline']:
				distance_dict[(patient, tumor)].append(distances[pair])
			for pair in rna_support['novel_somatic']:
				distance_dict[(patient, tumor)].append(distances[pair])

multisample = set([x[0] for x in file_dict if len([y for y in file_dict if x[0] in y]) > 1])
for patient in multisample:
	relevant_keys = [x for x in file_dict if x[0] == patient]
	relevant_entries = [file_dict[x] for x in relevant_keys]
	#print(relevant_entries)
	new_key = (patient, ';'.join([x[1] for x in relevant_keys]))
	new_entry = []
	for i in range(0, len(relevant_entries[0])):
		new_entry.append(str(median([float(x[i]) for x in relevant_entries])))
	file_dict[new_key] = new_entry
	assert len(new_entry) == len(relevant_entries[0])
	for key in relevant_keys:
		del file_dict[key]

with open(''.join([rna_path, 'paired_rna_support_summary.tsv']), 'w') as o:
	header = ['Patient', 'Tumor', 'All_relevant_germline_variants', 'All_relevant_somatic_variants', 
			  'Germline_variants_covered', 'Somatic_variants_covered', 'Germline_variants_validated', 
			  'Somatic_variants_validated', 'Somatic_variants_phased_with_germline', 'Somatic_variants_phased_with_somatic', 
			  'Total_somatic_variants_phased', 'Germline_phased_pairs', 'Somatic_phased_pairs', 
			  'Covered_germline_pairs', 'Covered_somatic_pairs',
			  'Supported_germline_phasing', 'Supported_somatic_phasing', 'Unsupported_germline_phasing', 
			  'Unsupported_somatic_phasing', 'Novel_germline_phasing', 'Novel_somatic_phasing',
			  'Covered_across_exon_pairs', 'Supported_across_exon_pairs', 'Novel_across_exon_pairs', 
			  'Not_supported_across_exons', 'Predicted_across_exon_pairs']
	print('\t'.join(header), file=o)
	with open(''.join([rna_path, 'novel_rna_distance.tsv']), 'w') as o2:
		header2 = ['Patient', 'Tumor', 'Distance']
		print('\t'.join(header2), file=o2)
		for (patient, tumor) in file_dict:
			out_line = [patient, tumor] + [str(x) for x in file_dict[(patient, tumor)]]
			print('\t'.join([str(x) for x in out_line]), file=o)
			if ';' in tumor:
				okay_tumor = random.sample(tumor.split(';'), 1)[0]
				distances = distance_dict[(patient, okay_tumor)]
			else:
				distances = distance_dict[(patient, tumor)]
			for dis in distances:
				print('\t'.join([patient, tumor, str(dis)]), file=o2)



