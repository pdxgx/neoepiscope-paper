#!/usr/bin/env python

import argparse
import glob
import pickle
import os
from operator import itemgetter
from collections import defaultdict
from neoepiscope import bowtie_index
from neoepiscope.transcript import get_transcripts_from_tree
from intervaltree import Interval, IntervalTree

parser = argparse.ArgumentParser()
parser.add_argument('-o', '--outputdir', type=str, required=True,
        help='path to write output_files'
    )
parser.add_argument('-n', '--neoepiscope-data', type=str, required=True,
        help='path to neoepiscope data dir'
    )
parser.add_argument('-c', '--hapcutdir', type=str, required=True,
        help='path to directory containing prepped haplotype files'
    )
parser.add_argument('-d', '--distance', type=str, required=True,
        help='distance to check for phasing (bp)'
    )
args = parser.parse_args()

with open(os.path.join(args.neoepiscope_data,'gencode_v19', 'intervals_to_transcript.pickle'), 'rb') as interval_stream:
    interval_dict = pickle.load(interval_stream)
with open(os.path.join(args.neoepiscope_data, 'gencode_v19', 'transcript_to_CDS.pickle'), 'rb') as interval_stream:
    cds_dict = pickle.load(interval_stream)
reference_index = bowtie_index.BowtieIndexReference(os.path.join('hg19'))

hapcut_results = glob.glob(''.join([args.hapcutdir, '*.hapcut.out.prepped']))

with open(os.path.join(args.outputdir, 'germline_distances.tsv'), 'w') as f:
    f.write('\t'.join(['Patient', 'Tumor', 'Distance']) + '\n')
with open(os.path.join(args.outputdir, 'somatic_distances.tsv'), 'w') as f:
    f.write('\t'.join(['Patient', 'Tumor', 'Distance']) + '\n')
with open(os.path.join(args.outputdir, 'combined_distances.tsv'), 'w') as f:
    f.write('\t'.join(['Patient', 'Tumor', 'Distance']) + '\n')
with open(os.path.join(args.outputdir, ''.join(['phasing_summary_', args.distance, '.tsv'])), 'w') as f:
    f.write('\t'.join(['Patient', 'Tumor', 'Total_somatic_mutations', 'Phased_germline', 
                        'Phased_somatic', 'Combined_phasing', 'No_haplotype', 
                        'Frameshift_affects', 'Total_eps', 'Phased_eps']) + '\n')

for res in hapcut_results:
    patient = '.'.join(res.split('/')[-1].replace('.hapcut.out.prepped', '').split('.')[0:-1])
    tumor = res.split('/')[-1].replace('.hapcut.out.prepped', '').split('.')[-1]
    total_somatic_mutations = 0
    somatic_phasing = []
    germline_phasing = []
    combined_phasing = []
    somatic_distances = []
    germline_distances = []
    combined_distances = []
    frameshift_counts = 0
    frameshift_distances = []
    stop_codon_muts = []
    with open(res, 'r') as r:
        block = []
        for line in r:
            if line.startswith('BLOCK'):
                continue
            elif line[0] == '*':
                if len(block) == 1:
                    if not block[0][6].endswith('*'):
                        total_somatic_mutations += 1
                elif [x[6] for x in block if not x[6].endswith('*')] != []:
                    transcripts = set()
                    strand = set()
                    i = 0
                    while i < (len(block) - 1):
                        if (block[i][7] == 'V' and block[i+1][7] == 'V' and not block[i][6].endswith('*') and not block[i+1][6].endswith('*')):
                            if (block[i][1] == (block[i+1][1] - 1) and block[i][4] == block[i+1][4] and block[i][5] == block[i+1][5]):
                                # They are adjacent to each other and have the same phasing info
                                block[i][2] = block[i][2] + block[i+1][2]
                                block[i][3] = block[i][3] + block[i+1][3]
                                del block[i+1]
                        i += 1                            
                    for i in range(0, len(block)):
                        overlapping_transcripts = get_transcripts_from_tree(''.join(['chr', block[i][0]]), block[i][1], block[i][8], interval_dict)
                        if len(overlapping_transcripts) > 0:
                            for tx in overlapping_transcripts:
                                tx_data = cds_dict[tx]
                                strand_info = set([x[4] for x in tx_data])
                                for st in strand_info:
                                    strand.add(st)
                            if not block[i][6].endswith('*'):
                                total_somatic_mutations += 1
                                germline_phase_count = 0
                                somatic_phase_count = 0
                                combined_phase_count = 0
                                g_distance = []
                                s_distance = []
                                c_distance = []
                                frameshift = False
                                for j in range(0, len(block)):
                                    if j < i:
                                        if block[j][6].endswith('*'):
                                            if block[i][4] == block[j][4] or block[i][5] == block[j][5]:
                                                if block[i][1] - block[j][8] < int(args.Distance):
                                                    germline_phase_count += 1
                                                    combined_phase_count += 1
                                                g_distance.append(block[i][1] - block[j][8])
                                                c_distance.append(block[i][1] - block[j][8])
                                        if (block[j][7] == 'I' and (len(block[j][3]) % 3)) or (block[j][7] == 'D' and (block[j][3] % 3)):
                                            if block[i][4] == block[j][4] or block[i][5] == block[j][5]:
                                                if '+' in strand:
                                                    frameshift = True
                                    elif j > i:
                                        if block[j][6].endswith('*'):
                                            if block[i][4] == block[j][4] or block[i][5] == block[j][5]:
                                                if block[j][1] - block[i][8] < int(args.Distance):
                                                    germline_phase_count += 1
                                                    combined_phase_count += 1
                                                g_distance.append(block[j][1] - block[i][8])
                                                c_distance.append(block[j][1] - block[i][8])
                                        else:
                                            if block[i][4] == block[j][4] or block[i][5] == block[j][5]:
                                                if block[j][1] - block[i][8] < int(args.Distance):
                                                    somatic_phase_count += 1
                                                    combined_phase_count += 1
                                                s_distance.append(block[j][1] - block[i][8])
                                                c_distance.append(block[j][1] - block[i][8])
                                        if (block[j][7] == 'I' and (len(block[j][3]) % 3)) or (block[j][7] == 'D' and (block[j][3] % 3)):
                                            if block[i][4] == block[j][4] or block[i][5] == block[j][5]:
                                                if '-' in strand:
                                                    frameshift = True
                                germline_phasing.append(germline_phase_count)
                                somatic_phasing.append(somatic_phase_count)
                                combined_phasing.append(combined_phase_count)
                                if g_distance != []:
                                    germline_distances.append(min(g_distance))
                                if s_distance != []:
                                    somatic_distances.append(min(s_distance))
                                if c_distance != []:
                                    combined_distances.append(min(c_distance))
                                if frameshift:
                                    frameshift_counts += 1
                        else:
                            total_somatic_mutations += 1
                # Reset
                block = []
            else:
                # Add mutation to block
                tokens = line.strip("\n").split()
                contig = tokens[3]
                if ',' in tokens[6]:
                    alternatives = tokens[6].split(',')
                else:
                    alternatives = [tokens[6]]
                for i in range(0, min(len(alternatives), 2)):
                    if len(tokens[5]) == len(alternatives[i]):
                        mutation_type = 'V'
                        pos = int(tokens[4])
                        ref = tokens[5]
                        alt = alternatives[i]
                        mut_size = len(tokens[5])
                        end = pos + mut_size
                    elif len(tokens[5]) > len(alternatives[i]):
                        mutation_type = 'D'
                        deletion_size = len(tokens[5]) - len(alternatives[i])
                        pos = int(tokens[4]) + (len(tokens[5]) - deletion_size)
                        ref = tokens[5]
                        alt = deletion_size
                        end = pos + deletion_size
                    elif len(tokens[5]) < len(alternatives[i]):
                        mutation_type = 'I'
                        insertion_size = len(alternatives[i]) - len(tokens[5])
                        pos = int(tokens[4]) + len(tokens[5]) - 1
                        ref = ''
                        alt = alternatives[i][len(tokens[5]):]
                        end = pos + 1
                    if len(alternatives) > 1:
                        if i == 0:
                            if tokens[1] == '1':
                                gen1 = '1'
                                gen2 = '0'
                            else:
                                gen1 = '0'
                                gen2 = '1'
                        elif i == 1:
                            if tokens[1] == '2':
                                gen1 = '1'
                                gen2 = '0'
                            else:
                                gen1 = '0'
                                gen2 = '1'
                    else:
                        gen1 = tokens[1]
                        gen2 = tokens[2]
                    block.append([contig, pos, ref, alt, gen1, gen2, tokens[7], mutation_type, end])
    with open(os.path.join(args.outputdir, 'germline_distances.tsv'), 'a') as f:
        refined_germline = sorted([x for x in germline_distances if x > 0])
        for dis in refined_germline:
            f.write('\t'.join([patient, tumor, disease, str(dis)]) + '\n')
    with open(os.path.join(args.outputdir, 'somatic_distances.tsv'), 'a') as f:
        refined_somatic = sorted([x for x in somatic_distances if x > 0])
        for dis in refined_somatic:
            f.write('\t'.join([patient, tumor, disease, str(dis)]) + '\n')
    with open(os.path.join(args.outputdir, 'combined_distances.tsv'), 'a') as f:
        refined_combo = sorted([x for x in combined_distances if x > 0])
        for dis in refined_combo:
            f.write('\t'.join([patient, tumor, disease, str(dis)]) + '\n')
    with open(os.path.join(args.outputdir, 'fs_distances.tsv'), 'a') as f:
        refined_fs = sorted([x for x in frameshift_distances if x > 0])
        for dis in refined_fs:
            f.write('\t'.join([patient, tumor, disease, str(dis)]) + '\n')
    with open(os.path.join(args.outputdir, ''.join(['phasing_summary_', args.distance, '.tsv'])), 'a') as f:
        out_line = [patient, tumor, disease, str(total_somatic_mutations), str(len([x for x in germline_phasing if x > 0])),
                        str(len([x for x in somatic_phasing if x > 0])), str(len([x for x in combined_phasing if x > 0])), 
                        str(frameshift_counts)]
        f.write('\t'.join(out_line) + '\n')
