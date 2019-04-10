#!/usr/bin/env python

from __future__ import print_function
import argparse
import copy
import glob
import os
import pickle
import sys
from operator import itemgetter
from bisect import bisect_left
from collections import defaultdict
from neoepiscope import bowtie_index
from neoepiscope.transcript import get_transcripts_from_tree, Transcript
from intervaltree import Interval, IntervalTree
from string import maketrans

# Reverse complement tool
revcomp_translation_table = maketrans("ATCG", "TAGC")

# Set up command line arguments
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
parser.add_argument('-d', '--distance', type=int, required=True,
        help='distance to check for phasing (bp)'
    )
args = parser.parse_args()

# Load reference data
with open(os.path.join(args.neoepiscope_data,'gencode_v19', 'intervals_to_transcript.pickle'), 'rb') as interval_stream:
    interval_dict = pickle.load(interval_stream)
with open(os.path.join(args.neoepiscope_data, 'gencode_v19', 'transcript_to_CDS.pickle'), 'rb') as interval_stream:
    cds_dict = pickle.load(interval_stream)
reference_index = bowtie_index.BowtieIndexReference(os.path.join(args.neoepiscope_data, 'hg19'))

# Generate interval tree for identifying stop codon overlap
stop_intervals = {}
for tx in cds_dict:
    stop_codons = [x for x in cds_dict[tx] if x[1] == 'stop_codon']
    if len(stop_codons) > 0:
        stop_codon = stop_codons[0]
        if stop_codon[0] not in stop_intervals:
            stop_intervals[stop_codon[0]] = IntervalTree()
        stop_intervals[stop_codon[0]][stop_codon[2]:stop_codon[3]+1] = tx

# Grab HapCUT2 results
hapcut_results = glob.glob(''.join([args.hapcutdir, '*.hapcut.out.prepped']))

# Write headers out to relevant files
with open(os.path.join(args.outputdir, 'germline_tx_distances.tsv'), 'w') as f:
    print('\t'.join(['Patient', 'Tumor', 'Distance']), file=f)
with open(os.path.join(args.outputdir, 'somatic_tx_distances.tsv'), 'w') as f:
    print('\t'.join(['Patient', 'Tumor', 'Distance']), file=f)
with open(os.path.join(args.outputdir, 'combined_tx_distances.tsv'), 'w') as f:
    print('\t'.join(['Patient', 'Tumor', 'Distance']), file=f)
with open(os.path.join(args.outputdir, ''.join(['phasing_summary_', str(args.distance), '.tsv'])), 'w') as f:
    print('\t'.join(['Patient', 'Tumor', 'Total_somatic_mutations', 'Total_germline_mutations', 'Phased_germline_no_tx', 
                        'Phased_somatic_no_tx', 'Combined_phasing_no_tx', 'Frameshift_phasing',
                        'Phased_germline_tx', 'Phased_somatic_tx', 'Combined_phasing_tx', 
                        'Nonstop_phasing']), file=f)

# Store patient variants
patient_variants = {}

# Process each patient's HapCUT2 results
for res in hapcut_results:
    # Get patient and tumor IDs
    patient = '.'.join(res.split('/')[-1].replace('.hapcut.out.prepped', '').split('.')[0:-1])
    tumor = res.split('/')[-1].replace('.hapcut.out.prepped', '').split('.')[-1]
    # Set up lists, counts, sets, and dictionaries
    pat_blocks = []
    homozygous_mutations = defaultdict(list)
    total_somatic_mutations = 0
    total_germline_mutations = 0
    no_tx_germline_phased = 0
    no_tx_somatic_phased = 0
    no_tx_combined_phased = 0
    somatic_tx_phasing = []
    germline_tx_phasing = []
    combined_tx_phasing = []
    somatic_tx_distances = []
    germline_tx_distances = []
    combined_tx_distances = []
    frameshift_counts = 0
    nonstop_phasing = set()
    variants = {'somatic-germline': defaultdict(int), 'somatic-somatic': defaultdict(int), 
                'frameshift': set(), 'nonstop': set()}
    # Open HapCUT2 results and process
    with open(res, 'r') as r:
        block = []
        for line in r:
            if line.startswith('BLOCK'):
                # Skip block header
                continue
            elif line[0] == '*':
                # Reached end of block
                if len(block) == 1:
                    # Only one variant in block - store block and tally variant counts
                    if not block[0][6].endswith('*'):
                        pat_blocks.append([block[0]])
                        total_somatic_mutations += 1
                    else:
                        total_germline_mutations += 1
                elif [x for x in block if not x[6].endswith('*')] != []:
                    # More than one variant and at least one somatic variant in block
                    block.sort(key=itemgetter(1))
                    i = 0
                    # Recompose decomposed substitutions
                    while i < (len(block) - 1):
                        if block[i][7] == 'V' and block[i+1][7] == 'V':
                            if (not block[i][6].endswith('*') and not block[i+1][6].endswith('*')) or (block[i][6].endswith('*') and block[i+1][6].endswith('*')):
                                if (block[i+1][1] == block[i][8] and block[i][4] == block[i+1][4] and block[i][5] == block[i+1][5]):
                                    # They are adjacent to each other and have the same phasing info
                                    block[i][2] = block[i][2] + block[i+1][2]
                                    block[i][3] = block[i][3] + block[i+1][3]
                                    block[i][8] = block[i+1][8]
                                    del block[i+1]
                        i += 1
                    # Store block and tally mutations
                    pat_blocks.append(block)
                    for i in range(len(block)):
                        if not block[i][6].endswith('*'):
                            total_somatic_mutations += 1
                        else:
                            total_germline_mutations += 1
                else:
                    # No somatic variants in block - tally germline mutations
                    total_germline_mutations += len(block)
                # Reset
                block = []
            else:
                # Mutation line - add mutation to block
                tokens = line.strip("\n").split()
                if ''.join(['chr', tokens[3]]) in reference_index.recs.keys():
                    contig = ''.join(['chr', tokens[3]])
                else:
                    contig = tokens[3]
                if ',' in tokens[6]:
                    alternatives = tokens[6].split(',')
                else:
                    alternatives = [tokens[6]]
                # Determine variant type, start/end, ref/alt for each alternative allele
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
                        ref = tokens[5][len(alternatives[i]):]
                        alt = deletion_size
                        end = pos + deletion_size
                    elif len(tokens[5]) < len(alternatives[i]):
                        mutation_type = 'I'
                        insertion_size = len(alternatives[i]) - len(tokens[5])
                        pos = int(tokens[4]) + len(tokens[5]) - 1
                        ref = ''
                        alt = alternatives[i][len(tokens[5]):]
                        end = pos + 1
                    # Process genotyping info
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
                    if gen1 != gen2:
                        # Heterozygous variant - add to block
                        block.append([contig, pos, ref, alt, gen1, gen2, tokens[7], mutation_type, end])
                    else:
                        # Homozygous variant - add to store
                        overlapping_transcripts = get_transcripts_from_tree(contig, pos, end, interval_dict)
                        if not tokens[7].endswith('*'):
                            total_somatic_mutations += 1
                        else:
                            total_germline_mutations += 1
                        for tx in overlapping_transcripts:
                            homozygous_mutations[tx].append([contig, pos, ref, alt, gen1, gen2, tokens[7], mutation_type, end])
    # Process all haplotype blocks
    for x in pat_blocks:
        block = copy.copy(x)
        tx_set = set()
        # Add relevant homozygous variants
        for i in range(0, len(block)):
            overlapping_transcripts = get_transcripts_from_tree(block[i][0], block[i][1], block[i][8], interval_dict)
            for tx in overlapping_transcripts:
                tx_set.add(tx)
                if tx in homozygous_mutations:
                    for homozygous in homozygous_mutations[tx]:
                        block.append(homozygous)
        block.sort(key=itemgetter(1))
        # Process through each somatic variant to determine phasing
        for i in range(len(block)):
            if not block[i][6].endswith('*'):
                germline_tx_phase_count = 0
                somatic_tx_phase_count = 0
                combined_tx_phase_count = 0
                g_tx_distance = []
                s_tx_distance = []
                c_tx_distance = []
                frameshift = False
                # Determine transcripts overlapping variant to check for stop codon overlap
                overlapping_transcripts = get_transcripts_from_tree(block[i][0], block[i][1], block[i][8], interval_dict)
                tx_objects = []
                for tx in overlapping_transcripts:
                    tx_data = cds_dict[tx]
                    # Create transcript object
                    transcript_ob = Transcript(reference_index, [[str(chrom), "blah", seq_type, str(start), str(end), ".", tx_strand]
                                                             for (chrom, seq_type, start, end, tx_strand, tx_type) in tx_data],
                                                    tx)
                    tx_objects.append(transcript_ob)
                    if transcript_ob.stop_codon and tuple([tuple(x) for x in block]) not in nonstop_phasing:
                        stop_overlap = False
                        for k in [4,5]:
                            if block[i][k] == '1':
                                phased_muts = [mut for mut in block if mut[k] == '1']
                                muts_to_apply = []
                                for mut in phased_muts:
                                    if len(get_transcripts_from_tree(mut[0], mut[1], mut[8], interval_dict)) > 0:
                                        muts_to_apply.append(mut)
                                overlapping_muts = []
                                for mut in muts_to_apply:
                                    # Check for stop codon overlaps
                                    stop_query = list(stop_intervals[mut[0]].overlap(mut[1], mut[8]))
                                    stop_options = [x.data for x in stop_query]
                                    if tx in stop_options:
                                        # Get mutation sequence
                                        mutation_seq_index = {}
                                        if mut[7] == 'V':
                                            for l in range(len(mut[3])):
                                                mutation_seq_index[mut[1]+l] = mut[3][l]
                                        elif mut[7] == 'D':
                                            mutation_seq = reference_index.get_stretch(mut[0], mut[1]-1, mut[8]-mut[1])
                                            for l in range(len(mutation_seq)):
                                                mutation_seq_index[mut[1]+l] = ''
                                        elif mut[7] == 'I':
                                            mutation_seq = mut[3]
                                            mutation_seq_index[mut[1]] = mut[3]
                                        # Get stop codon sequence
                                        stop_codon_seq = reference_index.get_stretch(mut[0], transcript_ob._stop_codon, 3)
                                        stop_codon_seq_index = {}
                                        for l in range(len(stop_codon_seq)):
                                            if transcript_ob.stop_codon+l in mutation_seq_index:
                                                if mut[7] == 'V' or mut[7] == 'D':
                                                    stop_codon_seq_index[transcript_ob.stop_codon+l] = mutation_seq_index[transcript_ob.stop_codon+l]
                                                else:
                                                    stop_codon_seq_index[transcript_ob.stop_codon+l] = ''.join([stop_codon_seq[l], mutation_seq_index[transcript_ob.stop_codon+l]])
                                            else:
                                                stop_codon_seq_index[transcript_ob.stop_codon+l] = stop_codon_seq[l]
                                        # Determine new sequence at stop codon position after applying mutation
                                        new_stop_seq = ''
                                        for l in range(len(stop_codon_seq)):
                                            new_stop_seq += stop_codon_seq_index[transcript_ob.stop_codon+l]
                                        new_stop_seq = new_stop_seq[0:3]
                                        if transcript_ob.rev_strand and new_stop_seq[::-1].translate(revcomp_translation_table) not in ['TAA', 'TGA', 'TAG']:
                                            # Mutation changes stop codon sequence on reverse strand
                                            stop_overlap = True
                                            overlapping_muts.append(mut)
                                        elif not transcript_ob.rev_strand and new_stop_seq not in ['TAA', 'TGA', 'TAG']:
                                            # Mutation changes stop codon sequence on forward strand
                                            stop_overlap = True
                                            overlapping_muts.append(mut)
                                if len(muts_to_apply) > 1 and stop_overlap:
                                    # There's more than one relevant mutation beside the nonstop
                                    post_stop_muts = []
                                    for mut in muts_to_apply:
                                        if mut not in overlapping_muts:
                                            if transcript_ob.rev_strand:
                                                if mut[8] < transcript_ob.stop_codon:
                                                    if (transcript_ob.stop_codon - mut[8]) <= args.distance:
                                                        # Mutation after but close to nonstop mutation
                                                        post_stop_muts.append(mut)
                                            else:
                                                if mut[1] > (transcript_ob.stop_codon+2):
                                                    if (mut[1] - (transcript_ob.stop_codon+2)) <= args.distance:
                                                        # Mutation after but close to nonstop mutation
                                                        post_stop_muts.append(mut)
                                    if len(post_stop_muts) > 0:
                                        # There is at least one mutation after the nonstop
                                        if (len([x for x in post_stop_muts if not x[6].endswith('*')]) > 0) or (len([x for x in overlapping_muts if not x[6].endswith('*')]) > 0):
                                            # Either nonstop or post-nonstop mutation is somatic
                                            nonstop_phasing.add(tuple([tuple(x) for x in block]))
                                            variants['nonstop'].add((tuple([tuple(x) for x in overlapping_muts]), tuple([tuple(x) for x in post_stop_muts])))
                # ID other phased variants
                for j in range(len(block)):
                    proceed = False
                    if j < i and block[j][6].endswith('*'):
                        # Second variant is germline variant before somatic
                        if (block[i][4] == block[j][4]) or (block[i][5] == block[j][5]):
                            # They are phased together
                            relevant_obs = []
                            # Check whether second variant overlaps any of the same transcripts as first
                            for tx_ob in tx_objects:
                                if (bisect_left(tx_ob.intervals, (block[j][1]-1)) % 2) or (bisect_left(tx_ob.intervals, (block[j][8]-1)) % 2):
                                    # Second variant overlaps transcript
                                    proceed = True
                                    relevant_obs.append(tx_ob)
                            if proceed:
                                # At least one tx shared
                                transcriptomic_distances = set()
                                for tx_ob in relevant_obs:
                                    intervals = tx_ob.intervals
                                    i_start_index = bisect_left(intervals, block[i][1]-1)
                                    i_end_index = bisect_left(intervals, block[i][8]-1)
                                    j_start_index = bisect_left(intervals, block[j][1]-1)
                                    j_end_index = bisect_left(intervals, block[j][8]-1)
                                    interval_length = len(intervals)
                                    if i_start_index != j_end_index:
                                        # Start and end of variants don't share exon - need to remove introns
                                        intron_to_remove = 0
                                        if (j_end_index % 2) and (i_start_index % 2):
                                            # End and start are inside exons
                                            for k in range(j_end_index, i_start_index, 2):
                                                intron_to_remove += (intervals[k+1] - intervals[k])
                                            tx_distance = block[i][1] - block[j][8] - intron_to_remove
                                            transcriptomic_distances.add(tx_distance)
                                        elif (j_end_index % 2) and (not i_start_index % 2):
                                            # End is in exon, start is in intron
                                            for k in range(j_end_index, i_start_index+1, 2):
                                                intron_to_remove += (intervals[k+1] - intervals[k])
                                            if block[i][7] == 'D':
                                                # Need to remove deleted sequence
                                                intron_to_remove += block[i][8] - (intervals[i_end_index - 1]+2)
                                            tx_distance = block[i][8] - block[j][8] - intron_to_remove
                                            transcriptomic_distances.add(tx_distance)
                                        elif (not j_end_index % 2) and (i_start_index % 2):
                                            # End is in intron, start is in exon
                                            for k in range(j_end_index-1, i_start_index, 2):
                                                intron_to_remove += (intervals[k+1] - intervals[k])
                                            if block[j][7] == 'D':
                                                # Need to remove deleted sequence
                                                intron_to_remove += (intervals[j_start_index] - (block[j][1]-2))
                                            tx_distance = block[i][1] - block[j][1] - intron_to_remove
                                            transcriptomic_distances.add(tx_distance)
                                        elif (not j_end_index % 2) and (not i_start_index % 2):
                                            # Neither is in exon
                                            for k in range(j_start_index-1, i_start_index+1, 2):
                                                intron_to_remove += (intervals[k+1] - intervals[k])
                                            if block[j][7] == 'D':
                                                intron_to_remove += (intervals[j_start_index] - (block[j][1]-2))
                                            if block[i][7] == 'D':
                                                intron_to_remove += (block[i][8] - (intervals[i_end_index - 1]+2))
                                            tx_distance = block[i][8] - block[j][1] - intron_to_remove
                                            transcriptomic_distances.add(tx_distance)
                                    else:
                                        # Don't need to consider introns in distance calculations
                                        transcriptomic_distances.add(block[i][1] - block[j][8])
                                    if block[j][7] == 'I' and (len(block[j][3]) % 3) and not tx_ob.rev_strand:
                                        # Second variant is upstream frameshifting indel
                                        frameshift = True
                                        variants['frameshift'].add((tuple(block[j]), tuple(block[i])))
                                    elif block[j][7] == 'D' and not tx_ob.rev_strand:
                                        # Second variant is upstream deletion
                                        if j_start_index == j_end_index:
                                            # Deletion is contained within one exon
                                            if block[j][3] % 3:
                                                # Frameshifting deletion
                                                frameshift = True
                                                variants['frameshift'].add((tuple(block[j]), tuple(block[i])))
                                        else:
                                            # Deletion spans outside one exon - need to take introns into account
                                            total_indel_size = 0
                                            if j_start_index % 2:
                                                total_indel_size += ((intervals[j_start_index]+1) - block[j][1] + 1)
                                            if j_end_index % 2:
                                                total_indel_size += (block[j][8] - (intervals[j_end_index-1]+2))
                                            if total_indel_size % 3:
                                                # Frameshift deletion
                                                frameshift = True
                                                variants['frameshift'].add((tuple(block[j]), tuple(block[i])))
                                # Find closest distance from all relevant transcripts
                                min_distance = min(transcriptomic_distances)
                                g_tx_distance.append(min_distance)
                                c_tx_distance.append(min_distance)
                                if min_distance < int(args.distance):
                                    # Variants are phased within distance
                                    germline_tx_phase_count += 1
                                    combined_tx_phase_count += 1
                                    pair = (tuple(block[j]), tuple(block[i]))
                                    variants['somatic-germline'][pair] = min_distance
                            else:
                                # Variant is phased, but not in the same transcript
                                no_tx_germline_phased += 1
                                no_tx_combined_phased += 1
                    elif j > i:
                        # Second variant is after the first variant
                        if (block[i][4] == block[j][4]) or (block[i][5] == block[j][5]):
                            # Variants are phased together
                            relevant_obs = []
                            # Check whether second variant overlaps any of the same transcripts as first
                            for tx_ob in tx_objects:
                                if (bisect_left(tx_ob.intervals, (block[j][1] - 1)) % 2) or (bisect_left(tx_ob.intervals, (block[j][8] - 1)) % 2):
                                    # Variant overlaps transcript
                                    proceed = True
                                    relevant_obs.append(tx_ob)
                            if proceed:
                                # At least one tx shared
                                transcriptomic_distances = set()
                                for tx_ob in relevant_obs:
                                    intervals = tx_ob.intervals
                                    i_start_index = bisect_left(intervals, block[i][1]-1)
                                    i_end_index = bisect_left(intervals, block[i][8]-1)
                                    j_start_index = bisect_left(intervals, block[j][1]-1)
                                    j_end_index = bisect_left(intervals, block[j][8]-1)
                                    interval_length = len(intervals)
                                    if i_end_index != j_start_index:
                                        # Start and end of variants don't share exon - need to remove introns
                                        intron_to_remove = 0
                                        if (j_start_index % 2) and (i_end_index % 2):
                                            # Start and end are inside exons
                                            for k in range(i_end_index, j_start_index, 2):
                                                intron_to_remove += (intervals[k+1] - intervals[k])
                                            tx_distance = block[j][1] - block[i][8] - intron_to_remove
                                            transcriptomic_distances.add(tx_distance)
                                        elif (j_start_index % 2) and (not i_end_index % 2):
                                            # Start is in exon but end is not
                                            for k in range(i_end_index-1, j_start_index, 2):
                                                intron_to_remove += (intervals[k+1] - intervals[k])
                                            if block[i][7] == 'D':
                                                # Need to remove deleted sequence
                                                intron_to_remove += ((intervals[i_start_index]+1) - (block[i][1]-1))
                                            tx_distance = block[j][1] - block[i][1] - intron_to_remove
                                            transcriptomic_distances.add(tx_distance)
                                        elif (not j_start_index % 2) and (i_end_index % 2):
                                            # Start is not in exon but end is
                                            for k in range(i_end_index, j_start_index, 2):
                                                intron_to_remove += (intervals[k+1] - intervals[k])
                                            if block[j][7] == 'D':
                                                # Need to remove deleted sequence
                                                intron_to_remove += (block[j][8] - (intervals[j_start_index]+2))
                                            tx_distance = block[j][8] - block[i][8] - intron_to_remove
                                            transcriptomic_distances.add(tx_distance)
                                        elif not j_start_index % 2 and not i_end_index % 2:
                                            # Neither start or end are in exons
                                            for k in range(i_end_index-1, j_start_index, 2):
                                                intron_to_remove += (intervals[k+1] - intervals[k])
                                            if block[j][7] == 'D':
                                                intron_to_remove += (block[j][8] - (intervals[j_start_index]+2))
                                            if block[i][7] == 'D':
                                                intron_to_remove += ((intervals[i_start_index]+1) - (block[i][1]-1))
                                            tx_distance = block[j][8] - block[i][1] - intron_to_remove
                                            transcriptomic_distances.add(tx_distance)
                                    else:
                                        # Start and end are in same exon - don't need to remove introns
                                        transcriptomic_distances.add(block[j][1] - block[i][8])
                                    if block[j][7] == 'I' and len(block[j][3]) % 3 and tx_ob.rev_strand:
                                        # Second variant is upstream frameshifting indel
                                        frameshift = True
                                        variants['frameshift'].add((tuple(block[j]), tuple(block[i])))
                                    elif block[j][7] == 'D' and tx_ob.rev_strand:
                                        # Second variant is upstream deletion
                                        if j_start_index == j_end_index:
                                            # Deletion is contained within one exon
                                            if block[j][3] % 3:
                                                # Frameshifting deletion
                                                frameshift = True
                                                variants['frameshift'].add((tuple(block[j]), tuple(block[i])))
                                        else:
                                            # Deletion spans outside one exon - need to take introns into account
                                            total_indel_size = 0
                                            if j_start_index % 2:
                                                total_indel_size += ((intervals[j_start_index]+1) - block[j][1] + 1)
                                            if j_end_index % 2:
                                                total_indel_size += (block[j][8] - (intervals[j_end_index-1]+2))
                                            if total_indel_size % 3:
                                                # Frameshifting deletion
                                                frameshift = True
                                                variants['frameshift'].add((tuple(block[j]), tuple(block[i])))
                                # Find closest distance from all relevant transcripts
                                min_distance = min(transcriptomic_distances)
                                if block[j][6].endswith('*'):
                                    g_tx_distance.append(min_distance)
                                else:
                                    s_tx_distance.append(min_distance)
                                c_tx_distance.append(min_distance)
                                if min_distance < int(args.distance):
                                    # Variants are phased within distance
                                    if block[j][6].endswith('*'):
                                        germline_tx_phase_count += 1
                                        pair = (tuple(block[i]), tuple(block[j]))
                                        variants['somatic-germline'][pair] = min_distance
                                    else:
                                        somatic_tx_phase_count += 1
                                        pair = (tuple(block[i]), tuple(block[j]))
                                        variants['somatic-somatic'][pair] = min_distance
                                    combined_tx_phase_count += 1
                            else:
                                # Variant is phased, but not in the same transcript
                                if block[j][6].endswith('*'):
                                    no_tx_germline_phased += 1
                                else:
                                    no_tx_somatic_phased += 1
                                no_tx_combined_phased += 1
                # Store data
                germline_tx_phasing.append(germline_tx_phase_count)
                somatic_tx_phasing.append(somatic_tx_phase_count)
                combined_tx_phasing.append(combined_tx_phase_count)
                if g_tx_distance != []:
                    germline_tx_distances.append(min(g_tx_distance))
                if s_tx_distance != []:
                    somatic_tx_distances.append(min(s_tx_distance))
                if c_tx_distance != []:
                    combined_tx_distances.append(min(c_tx_distance))
                if frameshift:
                    frameshift_counts += 1
    # Write to output files
    with open(os.path.join(args.outputdir, 'germline_tx_distances.tsv'), 'a') as f:
        refined_germline = sorted([x for x in germline_tx_distances if x > 0])
        for dis in refined_germline:
            print('\t'.join([patient, tumor, str(dis)]), file=f)
    with open(os.path.join(args.outputdir, 'somatic_tx_distances.tsv'), 'a') as f:
        refined_somatic = sorted([x for x in somatic_tx_distances if x > 0])
        for dis in refined_somatic:
            print('\t'.join([patient, tumor, str(dis)]), file=f)
    with open(os.path.join(args.outputdir, 'combined_tx_distances.tsv'), 'a') as f:
        refined_combo = sorted([x for x in combined_tx_distances if x > 0])
        for dis in refined_combo:
            print('\t'.join([patient, tumor, str(dis)]), file=f)
    with open(os.path.join(args.outputdir, ''.join(['phasing_summary_', str(args.distance), '.tsv'])), 'a') as f:
        out_line = [patient, tumor, str(total_somatic_mutations), str(total_germline_mutations), 
                        str(no_tx_germline_phased), str(no_tx_somatic_phased), str(no_tx_combined_phased),
                        str(frameshift_counts), str(len([x for x in germline_tx_phasing if x > 0])),
                        str(len([x for x in somatic_tx_phasing if x > 0])), str(len([x for x in combined_tx_phasing if x > 0])), 
                        str(len(nonstop_phasing))]
        print('\t'.join(out_line), file=f)

    patient_variants[(patient, tumor)] = variants

# Store pickled dictionary of patient variant pairs and their distances
with open(os.path.join(args.outputdir, ''.join(['patient_variants_', str(args.distance), '.pickle'])), 'wb') as f:
    pickle.dump(patient_variants, f)
