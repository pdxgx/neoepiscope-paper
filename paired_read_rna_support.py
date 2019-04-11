#!/usr/bin/env python

from __future__ import print_function
import argparse
import ast
import bisect
import os
import pickle
import pysam
import re
import subprocess
import sys
from collections import defaultdict
from intervaltree import Interval, IntervalTree
from neoepiscope import bowtie_index
from neoepiscope import paths
from neoepiscope.transcript import get_transcripts_from_tree, Transcript
from operator import itemgetter

# Taken from Rail-RNA: https://github.com/nellore/rail
def parsed_md(md):
    """ Divides an MD string up by boundaries between ^, letters, and numbers
        md: an MD string (example: 33A^CC).
        Return value: MD string split by boundaries described above.
    """
    md_to_parse = []
    md_group = [md[0]]
    for i, char in enumerate(md):
    	if i == 0: continue
    	if (re.match('[A-Za-z]', char) is not None) \
            != (re.match('[A-Za-z]', md[i-1]) is not None) or \
            (re.match('[0-9]', char) is not None) \
            != (re.match('[0-9]', md[i-1]) is not None):
            if md_group:
            	md_to_parse.append(''.join(md_group))
            md_group = [char]
        else:
        	md_group.append(char)
    if md_group:
    	md_to_parse.append(''.join(md_group))
    return [char for char in md_to_parse if char != '0']

# Taken from Rail-RNA: https://github.com/nellore/rail
def indels_junctions_exons_mismatches(
            cigar, md, pos, seq, drop_deletions=False, junctions_only=False
        ):
    """ Finds indels, junctions, exons, mismatches from CIGAR, MD string, POS
    
        cigar: CIGAR string
        md: MD:Z string
        pos: position of first aligned base
        seq: read sequence
        drop_deletions: drops deletions from coverage vectors iff True
        junctions_only: does not populate mismatch list
        Return value: tuple
            (insertions, deletions, junctions, exons, mismatches).
        Insertions is a list of tuples (last genomic position before insertion, 
                                 string of inserted bases). Deletions
            is a list of tuples (first genomic position of deletion,
                                 string of deleted bases). Junctions is a list
            of tuples (intron start position (inclusive),
                       intron end position (exclusive),
                       left_diplacement, right_displacement). Exons is a list
            of tuples (exon start position (inclusive),
                       exon end position (exclusive)). Mismatches is a list
            of tuples (genomic position of mismatch, read base)
    """
    insertions, deletions, junctions, exons, mismatches = [], [], [], [], []
    cigar = re.split(r'([MINDS])', cigar)[:-1]
    md = parsed_md(md)
    seq_size = len(seq)
    cigar_chars, cigar_sizes = [], []
    cigar_index, md_index, seq_index = 0, 0, 0
    max_cigar_index = len(cigar)
    while cigar_index != max_cigar_index:
    	if cigar[cigar_index] == 0:
    		cigar_index += 2
    		continue
    	if cigar[cigar_index+1] == 'M':
    		aligned_base_cap = int(cigar[cigar_index])
    		aligned_bases = 0
    		while True:
    			try:
    				aligned_bases += int(md[md_index])
    				if aligned_bases <= aligned_base_cap:
    					md_index += 1
    			except ValueError:
                    # Not an int, but should not have reached a deletion
    				assert md[md_index] != '^', '\n'.join(
                                                ['cigar and md:',
                                                 ''.join(cigar), ''.join(md)]
                                            )
    				if not junctions_only:
    					mismatches.append(
                                (pos + aligned_bases,
                                    seq[seq_index + aligned_bases])
                            )
    				correction_length = len(md[md_index])
    				m_length = aligned_base_cap - aligned_bases
    				if correction_length > m_length:
    					md[md_index] = md[md_index][:m_length]
    					aligned_bases = aligned_base_cap
    				else:
    					aligned_bases += correction_length
    					md_index += 1
    			if aligned_bases > aligned_base_cap:
    				md[md_index] = aligned_bases - aligned_base_cap
    				break
    			elif aligned_bases == aligned_base_cap:
    				break
            # Add exon
    		exons.append((pos, pos + aligned_base_cap))
    		pos += aligned_base_cap
    		seq_index += aligned_base_cap
    	elif cigar[cigar_index+1] == 'N':
    		skip_increment = int(cigar[cigar_index])
            # Add junction
    		junctions.append((pos, pos + skip_increment,
                            seq_index, seq_size - seq_index))
            # Skip region of reference
    		pos += skip_increment
    	elif cigar[cigar_index+1] == 'I':
            # Insertion
    		insert_size = int(cigar[cigar_index])
    		insertions.append(
                    (pos - 1, seq[seq_index:seq_index+insert_size])
                )
    		seq_index += insert_size
    	elif cigar[cigar_index+1] == 'D':
    		assert md[md_index] == '^', '\n'.join(
                                                ['cigar and md:',
                                                 ''.join(cigar), ''.join(md)]
                                            )
            # Deletion
    		delete_size = int(cigar[cigar_index])
    		md_delete_size = len(md[md_index+1])
    		assert md_delete_size >= delete_size
    		deletions.append((pos, md[md_index+1][:delete_size]))
    		if not drop_deletions: exons.append((pos, pos + delete_size))
    		if md_delete_size > delete_size:
    			md[md_index+1] = md[md_index+1][delete_size:]
    		else:
    			md_index += 2
            # Skip deleted part of reference
    		pos += delete_size
    	else:
    		assert cigar[cigar_index+1] == 'S'
    		seq_index += int(cigar[cigar_index])
    	cigar_index += 2
    # Merge exonic chunks/deletions; insertions/junctions could have chopped them up
    new_exons = []
    last_exon = exons[0]
    for exon in exons[1:]:
    	if exon[0] == last_exon[1]:
            # Merge ECs
    		last_exon = (last_exon[0], exon[1])
    	else:
            # Push last exon to new exon list
    		new_exons.append(last_exon)
    		last_exon = exon
    new_exons.append(last_exon)
    return insertions, deletions, junctions, new_exons, mismatches

def process_haplotypes(haplotypes, reference_index, interval_dict, cds_dict, bed_path, preprocessed_dis, patient):
    # Set up counters, sets, dictionaries
    som_count = 0
    no_tx_count = 0
    recomposed = set()
    abbreviated_muts_germ = {}
    abbreviated_muts_som = {}
    haplotype_dict = defaultdict(list)
    homozygous_mutations = defaultdict(set)
    all_variants = defaultdict(set)
    mutation_intervals = defaultdict(IntervalTree)
    distance_dict = {}
    homozygous_set = set()
    with open(bed_path, 'w') as o:
        with open(haplotypes, 'r') as f:
            block_transcripts = defaultdict(list)
            for line in f:
                if line.startswith('BLOCK'):
                    # Skip block headers
                    continue
                elif line[0] == '*':
                    # Process complete block
                    for transcript_id in block_transcripts:
                        block = block_transcripts[transcript_id]
                        block.sort(key=itemgetter(1))
                        if len(block) == 1:
                            if not block[0][6].endswith('*'):
                                # Add unphased somatic mutation to dictionary
                                all_variants[transcript_id].add(tuple(block[0]))
                                adj = (block[0][0], block[0][1], block[0][2], block[0][3], block[0][8])
                                abbreviated_muts_som[adj] = tuple(block[0])
                                mutation_intervals[block[0][0]][block[0][1]:block[0][8]] = adj
                                if contig in reference_index.recs.keys():
                                    o.write('\t'.join([block[0][0], str(block[0][1]), str(block[0][8])]) + '\n')
                        elif [x[6] for x in block if not x[6].endswith('*')] != []:
                            # If there are somatic mutations in the block, recombine decomposed mutations
                            i = 0
                            while i < (len(block) - 1):
                                if block[i][7] == 'V' and block[i+1][7] == 'V':
                                    if (not block[i][6].endswith('*') and not block[i+1][6].endswith('*')) or (block[i][6].endswith('*') and block[i+1][6].endswith('*')):
                                        if (block[i+1][1] == block[i][8] and block[i][4] == block[i+1][4] and block[i][5] == block[i+1][5]):
                                            block[i][2] = block[i][2] + block[i+1][2]
                                            block[i][3] = block[i][3] + block[i+1][3]
                                            block[i][8] = block[i+1][8]
                                            del block[i+1]
                                            recomposed.add(tuple(block[i]))
                                i += 1
                            for i in range(len(block)):
                                # Store info for all variants in block
                                all_variants[transcript_id].add(tuple(block[i]))
                                adj = (block[i][0], block[i][1], block[i][2], block[i][3], block[i][8])
                                if block[i][6].endswith('*'):
                                    abbreviated_muts_germ[adj] = tuple(block[i])
                                else:
                                    abbreviated_muts_som[adj] = tuple(block[i])
                                mutation_intervals[block[i][0]][block[i][1]:block[i][8]] = adj
                                if contig in reference_index.recs.keys():
                                    o.write('\t'.join([block[i][0], str(block[i][1]), str(block[i][8])]) + '\n')
                    # Reset the block
                    block_transcripts = defaultdict(list)
                else:
                    # Process mutation line
                    tokens = line.strip("\n").split()
                    if not tokens[7].endswith('*'):
                        som_count += 1
                    if ''.join(['chr', tokens[3]]) in reference_index.recs.keys():
                        contig = ''.join(['chr', tokens[3]])
                    else:
                        contig = tokens[3]
                    alternatives = tokens[6].split(',')
                    for i in range(0, min(len(alternatives), 2)):
                        if len(tokens[5]) == len(alternatives[i]):
                            # Subsitution
                            mutation_type = 'V'
                            pos = int(tokens[4])
                            mut_size = len(tokens[5])
                            end = pos + mut_size
                            ref = tokens[5]
                            alt = alternatives[i]
                        elif len(tokens[5]) > len(alternatives[i]):
                            # Deletion
                            mutation_type = 'D'
                            deletion_size = len(tokens[5]) - len(alternatives[i])
                            pos = int(tokens[4]) + (len(tokens[5]) - deletion_size)
                            end = pos + deletion_size
                            ref = tokens[5]
                            alt = deletion_size
                        elif len(tokens[5]) < len(alternatives[i]):
                            # Insertion
                            mutation_type = 'I'
                            pos = int(tokens[4]) + len(tokens[5]) - 1
                            end = pos + 1
                            ref = ''
                            alt = alternatives[i][len(tokens[5]):]
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
                        overlapping_transcripts = get_transcripts_from_tree(contig, pos, end, interval_dict)
                        if len(overlapping_transcripts) == 0 and not tokens[7].endswith('*'):
                            no_tx_count += 1
                        if gen1 != gen2:
                            # Heterozygous variant
                            for transcript_id in overlapping_transcripts:
                                block_transcripts[transcript_id].append([contig, pos, ref, alt, gen1, gen2, tokens[7], mutation_type, end])
                        else:
                            # Homozygous variant
                            for transcript_id in overlapping_transcripts:
                                homozygous_mutations[transcript_id].add((contig, pos, ref, alt, gen1, gen2, tokens[7], mutation_type, end))
                                homozygous_set.add((contig, pos, ref, alt, end))
                                all_variants[transcript_id].add((contig, pos, ref, alt, gen1, gen2, tokens[7], mutation_type, end))
                                if tokens[7].endswith('*'):
                                    abbreviated_muts_germ[(contig, pos, ref, alt, end)] = (contig, pos, ref, alt, gen1, gen2, tokens[7], mutation_type, end)
                                else:
                                    abbreviated_muts_som[(contig, pos, ref, alt, end)] = (contig, pos, ref, alt, gen1, gen2, tokens[7], mutation_type, end)
                                mutation_intervals[contig][pos:end] = (contig, pos, ref, alt, end)
                                if contig in reference_index.recs.keys():
                                    o.write('\t'.join([contig, str(pos), str(end)]) + '\n')
    # Look for variants that cross exons in the same transcript
    potential_across_exon_pairs = set()
    germline_phased_dict = defaultdict(set)
    somatic_phased_dict = defaultdict(set)
    for transcript_id in all_variants:
        # Initialize transcript object
        transcript_object = Transcript(reference_index,
                                       [[str(chrom), "blah", seq_type, str(start), str(end), ".", strand]
                                        for (chrom, seq_type, start, end, strand, tx_type) in cds_dict[transcript_id]],
                                       transcript_id)
        intervals = transcript_object.intervals
        interval_length = len(intervals)
        # Sort associated variants
        variants = list(all_variants[transcript_id])
        variants.sort(key=itemgetter(1))
        # Iterate through variants
        for i in range(len(variants)):
            if not variants[i][6].endswith('*'):
                # Get start and end index of variants
                i_start = bisect.bisect_left(intervals, variants[i][1]-1)
                i_end = bisect.bisect_left(intervals, variants[i][8]-1)
                adj_i = (variants[i][0], variants[i][1], variants[i][2], variants[i][3], variants[i][8])
                for j in range(len(variants)):
                    adj_j = (variants[j][0], variants[j][1], variants[j][2], variants[j][3], variants[j][8])
                    pair = tuple(sorted([adj_i, adj_j]))
                    j_start = bisect.bisect_left(intervals, variants[j][1]-1)
                    j_end = bisect.bisect_left(intervals, variants[j][8]-1)
                    if j < i and variants[j][6].endswith('*'):
                        # Second variant is germline variant upstream of somatic variant
                        if (variants[j], variants[i]) in preprocessed_dis[patient]['somatic-germline']:
                            distance_dict[pair] = preprocessed_dis[patient]['somatic-germline'][(variants[j], variants[i])]
                            germline_phased_dict[adj_i].add(adj_j)
                        else:
                            # Relevant variant
                            if i_start != j_end:
                                intron_to_remove = 0
                                if (j_end % 2) and (i_start % 2):
                                    for k in range(j_end, i_start, 2):
                                        intron_to_remove += (intervals[k+1] - intervals[k])
                                    tx_distance = variants[i][1] - variants[j][8] - intron_to_remove
                                elif (j_end % 2) and (not i_start % 2):
                                    for k in range(j_end, i_start+1, 2):
                                        intron_to_remove += (intervals[k+1] - intervals[k])
                                    if variants[i][7] == 'D':
                                        intron_to_remove += variants[i][8] - (intervals[i_end - 1]+2)
                                    tx_distance = variants[i][8] - variants[j][8] - intron_to_remove
                                elif (not j_end % 2) and (i_start % 2):
                                    for k in range(j_end-1, i_start, 2):
                                        intron_to_remove += (intervals[k+1] - intervals[k])
                                    if variants[j][7] == 'D':
                                        intron_to_remove += (intervals[j_start] - (variants[j][1]-2))
                                    tx_distance = variants[i][1] - variants[j][1] - intron_to_remove
                                elif (not j_end % 2) and (not i_start % 2):
                                    for k in range(j_start-1, i_start+1, 2):
                                        intron_to_remove += (intervals[k+1] - intervals[k])
                                    if variants[j][7] == 'D':
                                        intron_to_remove += (intervals[j_start] - (variants[j][1]-2))
                                    if variants[i][7] == 'D':
                                        intron_to_remove += (variants[i][8] - (intervals[i_end - 1]+2))
                                        tx_distance = variants[i][8] - variants[j][1] - intron_to_remove
                            else:
                                tx_distance = variants[i][1] - variants[j][8]
                            if pair in distance_dict:
                                distance_dict[pair] = min(distance_dict[pair], tx_distance)
                            else:
                                distance_dict[pair] = tx_distance
                    elif j > i:
                        # Second variant is downstream of somatic variant
                        if variants[j][6].endswith('*') and (variants[i], variants[j]) in preprocessed_dis[patient]['somatic-germline']:
                            distance_dict[pair] = preprocessed_dis[patient]['somatic-germline'][(variants[i], variants[j])]
                            germline_phased_dict[adj_i].add(adj_j)
                        elif not variants[j][6].endswith('*') and (variants[i], variants[j]) in preprocessed_dis[patient]['somatic-somatic']:
                            distance_dict[pair] = preprocessed_dis[patient]['somatic-somatic'][(variants[i], variants[j])]
                            somatic_phased_dict[adj_i].add(adj_j)
                            somatic_phased_dict[adj_j].add(adj_i)
                        else:
                            if i_end != j_start:
                                intron_to_remove = 0
                                if (j_start % 2) and (i_end % 2):
                                    for k in range(i_end, j_start, 2):
                                        intron_to_remove += (intervals[k+1] - intervals[k])
                                    tx_distance = variants[j][1] - variants[i][8] - intron_to_remove
                                elif (j_start % 2) and (not i_end % 2):
                                    for k in range(i_end-1, j_start, 2):
                                        intron_to_remove += (intervals[k+1] - intervals[k])
                                    if variants[i][7] == 'D':
                                        intron_to_remove += ((intervals[i_start]+1) - (variants[i][1]-1))
                                    tx_distance = variants[j][1] - variants[i][1] - intron_to_remove
                                elif (not j_start % 2) and (i_end % 2):
                                    for k in range(i_end, j_start, 2):
                                        intron_to_remove += (intervals[k+1] - intervals[k])
                                    if variants[j][7] == 'D':
                                        intron_to_remove += (variants[j][8] - (intervals[j_start]+2))
                                        tx_distance = variants[j][8] - variants[i][8] - intron_to_remove
                                elif (not j_start % 2) and (not i_end % 2):
                                    for k in range(i_end-1, j_start, 2):
                                        intron_to_remove += (intervals[k+1] - intervals[k])
                                    if variants[j][7] == 'D':
                                        intron_to_remove += (variants[j][8] - (intervals[j_start]+2))
                                    if variants[i][7] == 'D':
                                        intron_to_remove += ((intervals[i_start]+1) - (variants[i][1]-1))
                                    tx_distance = variants[j][8] - variants[i][1] - intron_to_remove
                            else:
                                tx_distance = variants[j][1] - variants[i][8]
                            if pair in distance_dict:
                                distance_dict[pair] = min(distance_dict[pair], tx_distance)
                            else:
                                distance_dict[pair] = tx_distance
                    if i_start != j_start and i_start != j_end and i_end != j_start and i_end != j_end:
                        # Variants aren't in the same exon
                        potential_across_exon_pairs.add(pair)
    return germline_phased_dict, somatic_phased_dict, potential_across_exon_pairs, abbreviated_muts_germ, abbreviated_muts_som, mutation_intervals, distance_dict, homozygous_set

### Taken with modifications from https://www.biostars.org/p/306041/
def read_pair_generator(bam, contig):
    """
    Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found.
    """
    read_dict = defaultdict(lambda: [None, None])
    for read in bam.fetch(contig=contig):
        if not read.is_proper_pair or read.is_secondary or read.is_supplementary:
            # Ignore reads that aren't mapped as pairs or from primary alignments
            continue
        if read.query[-2:] in ['.1', '.2']:
            # Strip '.1' or '.2' from query name
            qname = read.query_name[:-2]
        else:
            qname = read.query_name
        if qname not in read_dict:
            # Establishing new read pair
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            # Finish established read pair
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]


def rna_support_dict(pat, STAR_DIR, HAPCUT_DIR, DISTANCE_PICKLE, OUTPUT):
    # Set up input files
    hapcut = os.path.join(HAPCUT_DIR, ''.join([pat[0], pat[1], '.hapcut.out.prepped']))
    bed_file = os.path.join(STAR_DIR, pat[2], ''.join([pat[2], '.variant_regions.bed']))
    RNA_BAM = os.path.join(STAR_DIR, pat[2], ''.join([pat[2], '.sorted.bam']))
    new_RNA_BAM = os.path.join(STAR_DIR, pat[2], ''.join([pat[2], '.reduced.bam']))
    with open(os.path.join(paths.gencode_v19, "intervals_to_transcript.pickle"), "rb") as interval_stream:
        interval_dict = pickle.load(interval_stream)
    with open(os.path.join(paths.gencode_v19, "transcript_to_CDS.pickle"), "rb") as cds_stream:
        cds_dict = pickle.load(cds_stream)
    reference_index = bowtie_index.BowtieIndexReference(paths.bowtie_hg19)
    with open(DISTANCE_PICKLE, 'rb') as f:
        distance_dict = pickle.load(f)
    # Initialize sets to store covered/validated variants
    covered_germline = set()
    validated_germline = set()
    covered_somatic = set()
    validated_somatic = set()
    phased_somatic = set()
    phased_germline = set()
    phased_combined = set()
    # Initialize sets to store info on variant pairs
    supported_germline = set()
    supported_somatic = set()
    novel_germline = set()
    novel_somatic = set()
    covered_germline_pair = set()
    covered_somatic_pair = set()
    not_supported_germline = set()
    not_supported_somatic = set()
    covered_across_exon_pairs = set()
    validated_across_exon_pairs = set()
    novel_across_exon_pairs = set()
    not_supported_across_exon_pairs = set()
    transcript_dict = defaultdict(set)
    # Process haplotypes to create BED file
    germline_dict, somatic_dict, exon_pairs, hapcut_muts_germ, hapcut_muts_som, mut_intervals, distances, homozygous = process_haplotypes(
                                                                                                        hapcut, 
                                                                                                        reference_index, 
                                                                                                        interval_dict, 
                                                                                                        cds_dict, bed_file,
                                                                                                        distance_dict,
                                                                                                        (pat[0], pat[1]))
    print('Variants processed', file=sys.stderr)
    subprocess.call(['samtools', 'view', '-b', '-L', bed_file, '-o', new_RNA_BAM, RNA_BAM])
    subprocess.call(['samtools', 'index', new_RNA_BAM, ''.join([new_RNA_BAM, '.bai'])])
    print('BAM prepared', file=sys.stderr)
    BAM = pysam.AlignmentFile(new_RNA_BAM, 'rb')
    for contig in reference_index.recs.keys():
        print(contig, file=sys.stderr)
        for read1, read2 in read_pair_generator(BAM, contig):
            r1_tokens = str(read1).split('\t')
            r2_tokens = str(read2).split('\t')
            r1_tags = ast.literal_eval(r1_tokens[11])
            r2_tags = ast.literal_eval(r2_tokens[11])
            r1_md = [x[1] for x in r1_tags if x[0] == 'MD'][0]
            r2_md = [x[1] for x in r2_tags if x[0] == 'MD'][0]
            r1_cigar = r1_tokens[5]
            r2_cigar = r2_tokens[5]
            r1_pos = int(r1_tokens[3])+1
            r2_pos = int(r2_tokens[3])+1
            r1_seq = r1_tokens[9]
            r2_seq = r2_tokens[9]
            # Extract data from each read
            r1_insertions, r1_deletions, r1_junctions, r1_exons, r1_mismatches = indels_junctions_exons_mismatches(r1_cigar, r1_md, r1_pos, r1_seq)
            r2_insertions, r2_deletions, r2_junctions, r2_exons, r2_mismatches = indels_junctions_exons_mismatches(r2_cigar, r2_md, r2_pos, r2_seq)
            covered_mutations = set()
            # Add insertions to covered mutations if they are predicted
            for ins in r1_insertions + r2_insertions:
                full_ins = (contig, ins[0], '', ins[1], ins[0]+1)
                if full_ins in hapcut_muts_germ or full_ins in hapcut_muts_som:
                    covered_mutations.add(full_ins)
            # Add deletions to covered mutations is they are predicted
            for dele in r1_deletions + r2_deletions:
                full_dele = (contig, dele[0], dele[1], len(dele[1]), dele[0]+len(dele[1]))
                if full_dele in hapcut_muts_germ or full_dele in hapcut_muts_som:
                    covered_mutations.add(full_dele)
            # Add mismatches to covered mutations if they are predicted
            for mm in r1_mismatches + r2_mismatches:
                normal_base = reference_index.get_stretch(contig, mm[0]-1, len(mm[1]))
                full_mm = (contig, mm[0], normal_base, mm[1], mm[0]+len(mm[1]))
                if full_mm in hapcut_muts_germ or full_mm in hapcut_muts_som:
                    covered_mutations.add(full_mm)
            tx_set = set()
            # Find genomic intervals covered, accounting for junctions
            intervals = []
            for e in r1_exons + r2_exons:
                intervals.extend([e[0], e[1]+1])
            # Find with variants the read overlaps
            all_overlapping = set()
            # Check each interval for overlap
            for i in range(0, len(intervals), 2):
                overlapping_mutations = mut_intervals[contig].overlap(intervals[i], intervals[i+1])
                for mutation in [x.data for x in overlapping_mutations]:
                    all_overlapping.add(mutation)
                    if mutation in hapcut_muts_som:
                        # Somatic variant
                        covered_somatic.add(mutation)
                        if mutation in covered_mutations:
                            validated_somatic.add(mutation)
                    elif mutation in hapcut_muts_germ:
                        # Germline variant
                        covered_germline.add(mutation)
                        if mutation in covered_mutations:
                            validated_germline.add(mutation)
            # Process all variants overlapped by read
            all_overlapping = list(all_overlapping)
            all_overlapping.sort(key=itemgetter(1))
            for i in range(len(all_overlapping)):
                for j in range(len(all_overlapping)):
                    if i == j:
                        continue
                    elif i < j:
                        pair = tuple([all_overlapping[i], all_overlapping[j]])
                    elif j < i:
                        pair = tuple([all_overlapping[j], all_overlapping[i]])
                    if pair in distances:
                        if distances[pair] < 72:
                            if all_overlapping[i] in hapcut_muts_som:
                                # First variant is somatic
                                if all_overlapping[j] in hapcut_muts_som:
                                    # Both variants are somatic
                                    if all_overlapping[i] in somatic_dict[all_overlapping[j]]:
                                        # Variants are predicted to be phased
                                        covered_somatic_pair.add(pair)
                                        if pair in exon_pairs:
                                            covered_across_exon_pairs.add(pair)
                                        if all_overlapping[i] in covered_mutations and all_overlapping[j] in covered_mutations:
                                            # Pair validated by RNA-seq
                                            supported_somatic.add(pair)
                                            if pair in exon_pairs:
                                                validated_across_exon_pairs.add(pair)
                                        elif all_overlapping[i] in covered_mutations and all_overlapping[j] not in homozygous:
                                            # Phasing is not validated
                                            not_supported_somatic.add(pair)
                                            if pair in exon_pairs:
                                                not_supported_across_exon_pairs.add(pair)
                                        elif all_overlapping[j] in covered_mutations and all_overlapping[i] not in homozygous:
                                            # Phasing is not validated
                                            not_supported_somatic.add(pair)
                                            if pair in exon_pairs:
                                                not_supported_across_exon_pairs.add(pair)
                                else:
                                    # Second variant is germline
                                    if all_overlapping[j] in germline_dict[all_overlapping[i]]:
                                        # Variants are predicted to be phased
                                        covered_germline_pair.add(pair)
                                        if pair in exon_pairs:
                                            covered_across_exon_pairs.add(pair)
                                        if all_overlapping[i] in covered_mutations and all_overlapping[j] in covered_mutations:
                                            # Pair validated by RNA-seq
                                            supported_germline.add(pair)
                                            if pair in exon_pairs:
                                                validated_across_exon_pairs.add(pair)
                                        elif all_overlapping[i] in covered_mutations and all_overlapping[j] not in homozygous:
                                            # Phasing is not validated
                                            not_supported_germline.add(pair)
                                            if pair in exon_pairs:
                                                not_supported_across_exon_pairs.add(pair)
                            else:
                                # First variant is germline
                                if all_overlapping[j] in hapcut_muts_som:
                                    # Second variant is somatic
                                    if all_overlapping[i] in germline_dict[all_overlapping[j]]:
                                        # Variants are predicted to be phased
                                        covered_germline_pair.add(pair)
                                        if pair in exon_pairs:
                                            covered_across_exon_pairs.add(pair)
                                        if all_overlapping[i] in covered_mutations and all_overlapping[j] in covered_mutations:
                                            # Pair validated by RNA-seq
                                            supported_germline.add(pair)
                                            if pair in exon_pairs:
                                                validated_across_exon_pairs.add(pair)
                                        elif all_overlapping[j] in covered_mutations and all_overlapping[i] not in homozygous:
                                            # Phasing is not validated
                                            not_supported_germline.add(pair)
                                            if pair in exon_pairs:
                                                not_supported_across_exon_pairs.add(pair)
            # Process the variants present in the read
            covered_mutations = list(covered_mutations)
            covered_mutations.sort(key=itemgetter(1))
            for i in range(len(covered_mutations)):
                if covered_mutations[i] in hapcut_muts_som:
                    overlapping_txs = get_transcripts_from_tree(contig, covered_mutations[i][1], 
                                                                covered_mutations[i][4], interval_dict)
                    for tx in overlapping_txs:
                        transcript_dict[tx].add(covered_mutations[i])
                for j in range(len(covered_mutations)):
                    if i == j:
                        continue
                    elif i < j:
                        pair = tuple([covered_mutations[i], covered_mutations[j]])
                    elif j < i:
                        pair = tuple([covered_mutations[j], covered_mutations[i]])
                    if pair in distances:
                        if distances[pair] < 72:
                            if covered_mutations[i] in hapcut_muts_som:
                                # First variant is somatic
                                if covered_mutations[j] in hapcut_muts_som:
                                    # Second variant is somatic
                                    if covered_mutations[i] in somatic_dict[covered_mutations[j]]:
                                        # Variants are predicted to be phased
                                        supported_somatic.add(pair)
                                        if pair in exon_pairs:
                                            validated_across_exon_pairs.add(pair)
                                    else:
                                        # Novel variant phasing
                                        novel_somatic.add(pair)
                                        if pair in exon_pairs:
                                            # Pair is novelly phased across exons
                                            novel_across_exon_pairs.add(pair)
                                    # Store info about individual variants
                                    if j < i:
                                        phased_somatic.add(covered_mutations[j])
                                        phased_combined.add(covered_mutations[j])
                                    else:
                                        phased_somatic.add(covered_mutations[i])
                                        phased_combined.add(covered_mutations[i])
                                else:
                                    # Second variant is germline
                                    if covered_mutations[j] in germline_dict[covered_mutations[i]]:
                                        # Variants are predicted to be phased
                                        supported_germline.add(pair)
                                        if pair in exon_pairs:
                                            validated_across_exon_pairs.add(pair)
                                    else:
                                        # Novel variant phasing
                                        novel_germline.add(pair)
                                        if pair in exon_pairs:
                                            # Pair is novelly phased across exons
                                            novel_across_exon_pairs.add(pair)
                                    # Store info about individual variants
                                    phased_germline.add(covered_mutations[i])
                                    phased_combined.add(covered_mutations[i])
                            else:
                                # First variant is germline
                                if covered_mutations[j] in hapcut_muts_som:
                                    # Second variant is somatic
                                    if covered_mutations[i] in germline_dict[covered_mutations[j]]:
                                        # Variants are predicted to be phased
                                        supported_germline.add(pair)
                                        if pair in exon_pairs:
                                            validated_across_exon_pairs.add(pair)
                                    else:
                                        # Novel variant phasing
                                        novel_germline.add(pair)
                                        if pair in exon_pairs:
                                            # Pair is novelly phased across exons
                                            novel_across_exon_pairs.add(pair)
                                    phased_germline.add(covered_mutations[j])
                                    phased_combined.add(covered_mutations[j])
    # Remove germline mutations from covered/validated sets if they aren't on tumor-expressed transcripts
    for mut in list(covered_germline):
        somatic_expressed = False
        overlapping_txs = get_transcripts_from_tree(mut[0], mut[1], mut[4], interval_dict)
        for tx in overlapping_txs:
            if len(transcript_dict[tx]) > 0:
                somatic_expressed = True
        if not somatic_expressed:
            covered_germline.remove(mut)
            if mut in validated_germline:
                validated_germline.remove(mut)
    # Remove variant pairs from not supported sets if they aren't on tumor-expressed transcripts
    for pair in list(not_supported_somatic):
        if pair in supported_somatic:
            not_supported_somatic.remove(pair)
        else:
            somatic_expressed = False
            overlapping_txs = get_transcripts_from_tree(pair[0][0], pair[0][1], pair[1][4], interval_dict)
            for tx in overlapping_txs:
                if len(transcript_dict[tx]) > 0:
                    somatic_expressed = True
            if not somatic_expressed:
                not_supported_somatic.remove(pair)
    for pair in list(not_supported_germline):
        if pair in supported_germline:
            not_supported_germline.remove(pair)
        else:
            somatic_expressed = False
            overlapping_txs = get_transcripts_from_tree(pair[0][0], pair[0][1], pair[1][4], interval_dict)
            for tx in overlapping_txs:
                if len(transcript_dict[tx]) > 0:
                    somatic_expressed = True
            if not somatic_expressed:
                not_supported_germline.remove(pair)
    for pair in list(not_supported_across_exon_pairs):
        if pair in validated_across_exon_pairs:
            not_supported_across_exon_pairs.remove(pair)
        else:
            somatic_expressed = False
            overlapping_txs = get_transcripts_from_tree(pair[0][0], pair[0][1], pair[1][4], interval_dict)
            for tx in overlapping_txs:
                if len(transcript_dict[tx]) > 0:
                    somatic_expressed = True
            if not somatic_expressed:
                not_supported_across_exon_pairs.remove(pair)
    # Remove variant pairs from covered sets if they aren't on tumor-expressed transcripts
    for pair in list(covered_germline_pair):
        somatic_expressed = False
        overlapping_txs = get_transcripts_from_tree(pair[0][0], pair[0][1], pair[1][4], interval_dict)
        for tx in overlapping_txs:
            if len(transcript_dict[tx]) > 0:
                somatic_expressed = True
        if not somatic_expressed:
            covered_germline_pair.remove(pair)
    for pair in list(covered_across_exon_pairs):
        somatic_expressed = False
        overlapping_txs = get_transcripts_from_tree(pair[0][0], pair[0][1], pair[1][4], interval_dict)
        for tx in overlapping_txs:
            if len(transcript_dict[tx]) > 0:
                somatic_expressed = True
        if not somatic_expressed:
            covered_across_exon_pairs.remove(pair)
    # Remove pairs from general novel germline/somatic sets if they cross exon-exon junctions
    for pair in list(novel_germline):
        if pair in novel_across_exon_pairs:
            novel_germline.remove(pair)
    for pair in list(novel_somatic):
        if pair in novel_across_exon_pairs:
            novel_somatic.remove(pair)
    # Get total count of predicted phasing across exon-exon junctions
    predicted_across_exon_pairs = set()
    for pair in exon_pairs:
        if pair in distance_dict[(pat[0], pat[1])]['somatic-germline'] or pair in distance_dict[(pat[0], pat[1])]['somatic-somatic']:
            predicted_across_exon_pairs.add(pair)
    # Create dictionary
    rna_support = {}
    rna_support['all_relevant_germline'] = set(hapcut_muts_germ.keys())
    rna_support['all_relevant_somatic'] = set(hapcut_muts_som.keys())
    rna_support['total_germline_pairs'] = distance_dict[(pat[0], pat[1])]['somatic-germline'].keys()
    rna_support['total_somatic_pairs'] = distance_dict[(pat[0], pat[1])]['somatic-somatic'].keys()
    rna_support['covered_germline'] = covered_germline
    rna_support['covered_somatic'] = covered_somatic
    rna_support['validated_germline'] = validated_germline
    rna_support['validated_somatic'] = validated_somatic
    rna_support['phased_somatic'] = phased_somatic
    rna_support['phased_germline'] = phased_germline
    rna_support['phased_combined'] = phased_combined
    rna_support['covered_germline_pair'] = covered_germline_pair
    rna_support['covered_somatic_pair'] = covered_somatic_pair
    rna_support['supported_germline'] = supported_germline
    rna_support['supported_somatic'] = supported_somatic
    rna_support['novel_germline'] = novel_germline
    rna_support['novel_somatic'] = novel_somatic
    rna_support['not_supported_germline'] = not_supported_germline
    rna_support['not_supported_somatic'] = not_supported_somatic
    rna_support['predicted_across_exon_pairs'] = predicted_across_exon_pairs
    rna_support['validated_across_exon_pairs'] = validated_across_exon_pairs
    rna_support['novel_across_exon_pairs'] = novel_across_exon_pairs
    rna_support['covered_across_exon_pairs'] = covered_across_exon_pairs
    rna_support['not_supported_across_exon_pairs'] = not_supported_across_exon_pairs
    # Write pickled dictionary
    with open(os.path.join(OUTPUT, ''.join([pat[0], '.', pat[1],'.rna_support_updated.pickle'])), 'wb') as f:
        pickle.dump(rna_support, f)
    with open(os.path.join(OUTPUT, ''.join([pat[0], '.', pat[1],'.pair_distances.pickle'])), 'wb') as f:
        pickle.dump(distances, f)
    print('FINISHED!', file=sys.stderr)



parser = argparse.ArgumentParser()
parser.add_argument('-p', '--patient_info', type=str, required=True,
                    help='patient ID info'
    )
parser.add_argument('-s', '--star_output_dir', type=str, required=True,
                    help='path to directory containing STAR output subdirectories'
    )
parser.add_argument('-c', '--hapcut2_output_dir', type=str, required=True,
                    help='path to directory containing HapCUT2 output'
    )
parser.add_argument('-o', '--output_dir', type=str, required=True,
                    help='path to output directory'
    )
parser.add_argument('-d', '--pickled_dictionary', type=str, required=True,
                    help='path to 72 bp variant_to_distance pickled dictionary'
    )
args = parser.parse_args()

rna_support_dict(args.patient_info.split(','), args.star_output_dir, args.hapcut2_output_dir, 
                 args.pickled_dictionary, args.output_dir)

