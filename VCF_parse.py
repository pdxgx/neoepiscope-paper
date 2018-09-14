#!/usr/bin/env python

import datetime
import argparse
from intervaltree import Interval, IntervalTree

def combineVariants(*args):
	''' Parses a VCF to obtain variant entries and store as dictionary
		Keys in the dictionary are quadruples:
			(contig, position on contig, reference nucleotide, alternate nucleotide)
		Values are lists of the other data associated with that variant entry in the VCF,
			as well as the caller which produced the vcf
		If more than one caller called the variant, the associated data is added to the entry in the dictionary
		Unique header lines are also retained as a list for subsequent processing
		
		args: tuple of (caller, vcf)
			vcf: path vcf which has been normalized and decomposed using vt
			caller: name of the caller which produced the vcf
		
		Return values: variant dictionary, header info list
	'''
	# Establish dictionary to hold processed variants
	variants = {}
	headers = []
	# Parse each set
	for set in args:
		caller = set[0]
		vcf = set[1]
		with open(vcf, "r") as f:
			# Establish tumor sample as last column by default
			tumor_last = True
			# Filter calls to throw out - REJECT for MuTect, Tier 5 for MuSE
			reject_filters = ["REJECT", "Tier5", "DETP20", "IRC", "MMQS100", "MMQSD50", 
							  "MQD30", "MVC4", "MVF5", "NRC", "PB10", "RLD25", "SB1", "str10"]
			# Parse VCF
			for line in f:
			
				# Extract data from header lines
				if line[0:2] == "##":
					line = line.strip("\n")
					if "ID=" in line and line not in headers:
						headers.append(line)
			
				# Determine sample order from column headers
				elif line[0] == "#":
					line = line.strip("\n").split("\t")
					if "NORM" in line[-1] or "_N" in line[-1]:
						tumor_last = False
			
				# Parse variant entries
				else:
					line = line.strip("\n").split("\t")
					filter = line[6]
					# Check that variant passes filters before proceeding
					if len([item for item in reject_filters if (item in filter)]) == 0:
						contig = line[0]
						if contig.isdigit() == True:
							contig = int(contig)
						position = int(line[1])
						id = line[2]
						ref = line[3]
						alt = line[4]
						if ref > alt and len(ref) > 5000:
							too_large = True
						else:
							too_large = False
						qual = line[5]
						info = line[7]
						format = line[8]
						sample1 = line[9]
						sample2 = line[10]
						if tumor_last == True:
							normal = sample1
							tumor = sample2
						else:
							normal = sample2
							tumor = sample1
						entry = (contig, position, ref, alt)
						data = [id, qual, filter, info, format, normal, tumor, caller]
						if entry not in variants and not too_large:
							variants[entry] = [data]
						elif entry in variants and data not in variants[entry] and not too_large:
							variants[entry].append(data)
	
	return variants, headers


def parseHeaders(header_list):
	''' Parses a list of header lines to produce lists of format, filter, and info definitions
			
		header_list: list of header lines to process
		
		Return values: format, filter, and info lists
	'''
	# Establish output dictionaries
	format_list = []
	filter_list = []
	info_list = []
	
	# Iterate through header lines
	for item in header_list:
		
		# Store format header line	
		if 'FORMAT=<' in item:
			format_list.append(item)
		
		# Store filter header line	
		elif 'FILTER=<' in item:
			filter_list.append(item)
		
		# Store info header line		
		elif 'INFO=<' in item:
			info_list.append(item)
	
	return format_list, filter_list, info_list


def writeVCF(number, variants, outdir, sample, callers, format, filter, info):
	''' Write a consensus VCF from a refined set of variants
		Also produces output file with descriptive data about the call sets
			Number of callers used and which callers were used
		
		variants = dictionary of filtered variants
		outdiir = path to output directory
		sample = sample name to attach to output files
		callers = list of callers used to generate call sets
		format = header format definition dictionary
		filter = header filter definition dictionary
		info = header info definition dictionary
		
		Return value: none
	'''
	# Establish output files
	outputVCF = outdir + "/" + sample + ".consensus.vcf"
	outputCallerCombos = outdir + "/" + sample + ".caller.combos.txt"
	outputCallerCounts = outdir + "/" + sample + ".caller.counts.txt"
	
	# Establish data dictionaries
	caller_combos = {}
	number_caller_counts = {}
	
	# Establish header info
	today = datetime.date.today()
	caller_headers = []
	for caller in callers:
		line = '##FILTER=<ID=' + caller.upper() + ',Description="Variant called by ' + caller.lower() + '">\n'
		caller_headers.append(line)
	
	with open(outputVCF, "w") as f:
		# Write header lines
		f.write("##fileformat=VCF4.2\n")
		f.write("##fileDate=" + str(today) + "\n")
		# Write format header lines
		for ID in format:
			f.write(ID+'\n')
		f.write('##FORMAT=<ID=CALLER,Number=.,Type=String,Description="Caller used to produce data"\n')
		# Write info header lines
		for ID in info:
			f.write(ID+'\n')
		# Write filter header lines
		for ID in filter:
			f.write(ID+'\n')
		for header in caller_headers:
			f.write(header)
		# Write column header lines
		f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR\n")
		
		# Loop through each variant to collect data and write to VCF file if relevant
		pindel_variants = {}
		for entry in variants.keys():
			for callset in variants[entry]:
				caller = callset[7].upper()
				if caller == "PINDEL" and len(entry[2]) > 1:
					chrom = entry[0]
					if chrom not in pindel_variants:
						pindel_variants[chrom] = IntervalTree()
					if len(variants[entry]) == 1:
						pindel_variants[chrom][int(entry[1]):int(entry[1])+len(entry[2])] = [entry[0], entry[1], callset[0], 
																							 entry[2], entry[3], callset[1], 
																							 callset[2], callset[3], callset[4], 
																							 callset[5], callset[6], callset[7], 0]
					else:
						pindel_variants[chrom][int(entry[1]):int(entry[1])+len(entry[2])] = [entry[0], entry[1], callset[0], 
																							 entry[2], entry[3], callset[1], 
																							 callset[2], callset[3], callset[4], 
																							 callset[5], callset[6], callset[7], 1]
						outline = '\t'.join([str(x) for x in pindel_variants[chrom][int(entry[1]):int(entry[1])+len(entry[2])][:-1]]) + '\n'
						f.write(outline)
                
		for entry in sorted(variants.keys()):
			
			# Collect data on which caller(s) called the variant and store data on pindel structural events
			caller_list = []
			for callset in variants[entry]:
				caller = callset[7].upper()
				caller_list.append(caller)

			caller_list.sort()
			callers_used = ",".join(caller_list)
			if callers_used not in caller_combos:
				caller_combos[callers_used] = 1
			else:
				caller_combos[callers_used] += 1
			
			# Collect data on how many callers called the variant
			if len(caller_list) in number_caller_counts:
				number_caller_counts[len(caller_list)] += 1
			else:
				number_caller_counts[len(caller_list)] = 1
			
			if "PINDEL" not in caller_list and entry[0] in pindel_variants:
				overlap = pindel_variants[entry[0]][int(entry[1]):int(entry[1])+max(len(entry[2]), len(entry[3]))]
				if len(overlap) > 0:
					pindel_overlap = True
				else:
					pindel_overlap = False
			else:
				pindel_overlap = False

			# If at least the required number of callers called the variant, write to output VCF
			if len(caller_list) >= number and not pindel_overlap:
				
				# Initialize lists to store each data type
				set_callers = []
				ID_list = []
				qual_list = []
				filter_list = []
				info_list = []
				format_list = []
				tumor_data = []
				normal_data = []
				
				# Loop through data from each caller that called the variant
				for callset in variants[entry]:
					# Store data on callers
					caller = callset[7].upper()
					set_callers.append(caller)
					
					# Store data on IDs
					if callset[0] != "." and callset[0] not in ID_list:
						ID_list.append(callset[0])
				
					# Store data on quality scores
					if callset[1] != "." and callset[1] != 0:
						qual_list.append(float(callset[1]))
					
					# Store data on filters
					if ";" in callset[2]:
						filters = callset[2].split(";")
						for item in filters:
							if item not in filter_list:
								filter_list.append(item)
					elif ";" not in callset[2] and callset[2] not in filter_list and callset[2] != "PASS":
						filter_list.append(callset[2])
				
					# Store data on info
					if callset[3] != "" and callset[3] != ".":
						info_items = callset[3].split(";")
						for item in info_items:
							if item not in info_list:
								info_list.append(item)
							
					# Adjust and store genotype fields
					if "FREQ" not in callset[4]:
						callset[4] = callset[4] + ":FREQ:CALLER"
						format_list.append(callset[4])
						callset[5] = callset[5] + ":.:" + caller
						normal_data.append(callset[5])
						callset[6] = callset[6] + ":.:" + caller
						tumor_data.append(callset[6])
					else:
						callset[4] = callset[4] + ":CALLER"
						format_list.append(callset[4])
						callset[5] = callset[5] + ":" + caller
						normal_data.append(callset[5])
						callset[6] = callset[6] + ":" + caller
						tumor_data.append(callset[6])
			
				# Finalize ID data
				if len(ID_list) == 0:
					ID = "."
				else:
					ID = ";".join(ID_list)
				
				# Finalize quality score data
				if len(qual_list) != 0:
					QUAL = sum(qual_list)/len(qual_list)
				else:
					QUAL = "."
				
				# Finalize filter data
				set_callers = list(set(set_callers))	
				if len(filter_list) != 0:
					FILTER = ";".join(filter_list) + ";" + ";".join(set_callers)
				else:
					FILTER = ";".join(set_callers)
				
				# Finalize info data
				if len(info_list) != 0:
					INFO = ";".join(info_list)
				else:
					INFO = "."
				
				# Finalize genotype fields
				FORMAT = ";".join(format_list)
				NORMAL = ";".join(normal_data)
				TUMOR = ";".join(tumor_data)
				
				# Write data to VCF
				outline = "\t".join([str(entry[0]), str(entry[1]), ID, entry[2], entry[3], str(QUAL), FILTER, INFO, FORMAT, NORMAL, TUMOR]) + "\n"
				f.write(outline)

			elif pindel_overlap:
				for event in overlap:
					if event.data[-1] == 0:
						outline = "\t".join([str(x) for x in event.data[:-1]]) + "\n"
						f.write(outline)
						new_data = event.data[:-1] + [1]
						start, end = event.begin, event.end
						pindel_variants[entry[0]].removei(event.begin, event.end, event.data)
						pindel_variants[entry[0]][start:end] = new_data
	
	# Write count data for caller combos
	with open(outputCallerCombos, "w") as f:
		for combo in caller_combos:
			line = combo + "\t" + str(caller_combos[combo]) + "\n"
			f.write(line)
	
	# Write count data for number of callers per variant
	with open(outputCallerCounts, "w") as f:
		for count in number_caller_counts:
			line = str(count) + "\t" + str(number_caller_counts[count]) + "\n"
			f.write(line)

			
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--vcfs', type=str, required=True,
            help='comma-separated list of paths to input vcf'
        )
    parser.add_argument('-c', '--callers', type=str, required=True,
            help='comma-separated list of callers which generated the input vcfs (in same order)'
        )
    parser.add_argument('-o', '--outputdir', type=str, required=True,
            help='path to output consensus vcf'
        )
    parser.add_argument('-n', '--number', type=int, required=False,
            help='number of callers required to maintain a variant'
        )
    parser.add_argument('-s', '--sample', type=str, required=False,
            help='unique sample name for output'
        )	
    args = parser.parse_args()
    
    # Turn comma separated VCF and caller lists into list objects
    vcf_list = args.vcfs.split(",")
    caller_list = args.callers.split(",")
    
    # Turn VCF and callers lists into a list of tuples for combineVariants()
    set_list = []
    for i in range(0, len(caller_list)):
    	this_set = (caller_list[i], vcf_list[i])
    	set_list.append(this_set)
    
    # Produce a variant dictionary and list of header data	
    variant_dict, header_data = combineVariants(*set_list)
    
    # Produce a dictionary of unique header data
    format_dict, filter_dict, info_dict = parseHeaders(header_data)
    
    # Write consensus VCF and associated metadata
    writeVCF(args.number, variant_dict, args.outputdir, args.sample, caller_list, format_dict, filter_dict, info_dict)
