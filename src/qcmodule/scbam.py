#!/usr/bin/env python3

# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2010–2026 Liguo Wang

# -*- coding: utf-8 -*-
"""
Analzye 10X genomics single cell BAM files.
"""

import pysam
import sys
import collections
import logging
import re
import pandas as pd
import subprocess
import numpy as np

def diff_str(s1, s2):
	'''
	Comparing orignal barcode to the corrected barcode
	find the index and the nucleotide that has been corrected.


	Parameters
	----------
	s1 : str
		the original barcode
	s2 : str
		the corrrected barcode

	'''
	results = []
	if len(s1) != len(s2):
		return results
	diff_positions = [i for i in range(len(s1)) if s1[i] != s2[i]]
	for pos in diff_positions:
		results.append([pos, s1[pos], s2[pos]])
	return results

def read_match_type(cigar_str):
	'''return the matching type between read and ref'''
	match_type = ''
	if bool(re.search(r'\A\d+M\Z', cigar_str)):
		match_type ='Map_consecutively'
	elif bool(re.search(r'\A\d+M\d+N\d+M\Z', cigar_str)):
		match_type = 'Map_with_splicing'
	elif bool(re.search(r'\A\d+S\d+M\Z', cigar_str)):
		match_type = 'Map_with_clipping'
	elif bool(re.search(r'\A\d+M\d+S\Z', cigar_str)):
		match_type = 'Map_with_clipping'
	elif bool(re.search(r'\A\d+M\d+N\d+M\d+S\Z', cigar_str)):
		match_type = 'Map_with_splicing_and_clipping'
	elif bool(re.search(r'\A\d+S\d+M\d+N\d+M\Z', cigar_str)):
		match_type = 'Map_with_splicing_and_clipping'
	else:
		match_type = 'Others'
	return match_type

def list2str (lst):
	'''
	translate samtools returned cigar_list into cigar_string
	'''
	code2Char={'0':'M','1':'I','2':'D','3':'N','4':'S','5':'H','6':'P','7':'=','8':'X'}
	cigar_str=''
	for i in lst:
		cigar_str += str(i[1]) + code2Char[str(i[0])]
	return cigar_str

def barcode_edits(infile, outfile, step_size=10000, limit=2000000, CR_tag = 'CR', CB_tag = 'CB', UR_tag = 'UR', UB_tag = 'UB'):
	'''
	Analzye barcode in BAM file.

	Parameters
	----------
	infile : str
		Input BAM file. Must be sorted and indexed.
	outfile : str
		Prefix of output files.
	step_size: int
		Output progress report when step_size alignments have been processed.
	limit : int
		Only process this number of alignments and stop.
	'''
	logging.info("Reading BAM file \"%s\" ..." % infile)
	samfile = pysam.Samfile(infile,'rb')

	CB_miss = 0 #number of reads without cell barcode
	CB_same = 0 #number of reads whose original cell barcode same as edited barcode
	CB_diff = 0 #number of reads whose cell barcode has been edited
	CB_freq = collections.defaultdict(int) # cell barcode: raw reads
	CB_corrected_bases = collections.defaultdict(dict)


	UMI_miss = 0
	UMI_same = 0
	UMI_diff = 0
	UMI_freq = collections.defaultdict(int) # UMI : raw reads
	UMI_corrected_bases = collections.defaultdict(dict)

	total_alignments = 0
	try:
		while(1):
			total_alignments += 1
			aligned_read = next(samfile)
			tag_dict = dict(aligned_read.tags) #{'NM': 1, 'RG': 'L1'}

			original_CB = ''
			corrected_CB = ''
			if CR_tag in tag_dict and CB_tag in tag_dict:
				original_CB = tag_dict[CR_tag].replace('-1','')
				corrected_CB = tag_dict[CB_tag].replace('-1','')
				CB_freq[corrected_CB] +=1
				if original_CB != corrected_CB:
					CB_diff += 1
					for diff in diff_str(original_CB, corrected_CB):
						try:
							CB_corrected_bases[diff[0]][diff[1]+ ':' + diff[2]] += 1
						except:
							CB_corrected_bases[diff[0]][diff[1]+ ':' + diff[2]] = 1
				else:
					CB_same +=1
			else:
				CB_miss += 1

			original_UMI = ''
			corrected_UMI = ''
			if UR_tag in tag_dict and UB_tag in tag_dict:
				original_UMI = tag_dict[UR_tag].replace('-1','')
				corrected_UMI = tag_dict[UB_tag].replace('-1','')
				UMI_freq[corrected_UMI] += 1
				if original_UMI != corrected_UMI:
					UMI_diff += 1
					for diff in diff_str(original_UMI, corrected_UMI):
						try:
							UMI_corrected_bases[diff[0]][diff[1]+ ':' + diff[2]] += 1
						except:
							UMI_corrected_bases[diff[0]][diff[1]+ ':' + diff[2]] = 1
				else:
					UMI_same += 1
			else:
				UMI_miss += 1

			if total_alignments % step_size == 0:
				print("%d alignments processed.\r" % total_alignments, end=' ', file=sys.stderr)
			if limit is not None:
				if total_alignments >= limit:
					break

	except StopIteration:
		pass
	logging.info ("Total alignments processed: %d" % total_alignments)

	logging.info("Number of alignmenets with <cell barcode> kept AS IS: %d" % CB_same)
	logging.info("Number of alignmenets with <cell barcode> edited: %d" % CB_diff)
	logging.info("Number of alignmenets with <cell barcode> missing: %d" % CB_miss)
	logging.info("Number of alignmenets with UMI kept AS IS: %d" % UMI_same)
	logging.info("Number of alignmenets with UMI edited: %d" % UMI_diff)
	logging.info("Number of alignmenets with UMI missing: %d" % UMI_miss)

	# writing cell barcode
	logging.info ("Writing cell barcode frequencies to \"%s\"" % (outfile + '.CB_freq.tsv'))
	with open(outfile + '.CB_freq.tsv','w') as CB_OUT:
		for bc,count in sorted(CB_freq.items(), key=lambda item: item[1], reverse=True):
			CB_OUT.write(bc + '\t' + str(count) + '\n')

	# writing UMI
	logging.info ("Writing UMI frequencies to \"%s\"" % (outfile + '.UMI_freq.tsv'))
	with open(outfile + '.UMI_freq.tsv','w') as UMI_OUT:
		for bc,count in sorted(UMI_freq.items(), key=lambda item: item[1], reverse=True):
			UMI_OUT.write(bc + '\t' + str(count) + '\n')

	CB_mat_file = outfile + '.CB_edits_count.csv'
	logging.info ("Writing the nucleotide editing matrix (count) of cell barcode to \"%s\"" % CB_mat_file)
	CB_diff_mat = pd.DataFrame.from_dict(CB_corrected_bases)
	CB_diff_mat = CB_diff_mat.fillna(0)
	#CB_diff_mat = CB_diff_mat.T
	CB_diff_mat.sort_index(inplace=True)
	CB_diff_mat = CB_diff_mat.sort_index(axis=1)
	CB_diff_mat.index.name='Edits'
	CB_diff_mat.to_csv(CB_mat_file, index=True, index_label="Index")

	UMI_mat_file = outfile + '.UMI_edits_count.csv'
	logging.info ("Writing the nucleotide editing matrix of molecular barcode (UMI) to \"%s\"" % UMI_mat_file)
	UMI_diff_mat = pd.DataFrame.from_dict(UMI_corrected_bases)
	UMI_diff_mat = UMI_diff_mat.fillna(0)
	#UMI_diff_mat = UMI_diff_mat.T
	UMI_diff_mat.sort_index(inplace=True)
	UMI_diff_mat = UMI_diff_mat.sort_index(axis=1)
	UMI_diff_mat.index.name='Edits'
	UMI_diff_mat.to_csv(UMI_mat_file, index=True, index_label="Index")



def mapping_stat(infile,  step_size=50000, CB_tag = 'CB', UMI_tag = 'UB',RE_tag = 'RE', TX_tag = 'TX', AN_tag = 'AN', xf_tag = 'xf', chrM_id ='chrM', n_thread = 1):
	'''
	Reads mapping statistics

	Parameters
	----------
	infile : str
		Input BAM file. Must be sorted and indexed.
	outfile : str
		Prefix of output files. If outfile is None, only do counting and do not generate BAM files.
	step_size: int
		Output progress report when step_size alignments have been processed.
	limit : int
		Only process this number of alignments and stop.
	'''
	logging.info("Reading BAM file \"%s\" ..." % infile)
	try:
		#older versions of pysam
		samfile = pysam.AlignmentFile(infile, mode = 'rb',require_index = True, thread = n_thread)
	except:
		#latest verion of pysam (v0.19.1)
		samfile = pysam.AlignmentFile(infile, mode = 'rb',require_index = True, threads = n_thread)
	if samfile.check_index():
		pass
	else:
		logging.error("Cannot find the index file")
		sys.exit(0)
	chrom_info = zip(samfile.references, samfile.lengths) #[('chr1', 195471971), ('chr10', 130694993),...]

	total_alignments = 0
	confi_alignments = 0

	total_reads_n = 0
	confi_reads_n = 0

	## for confidently mapped reads
	#PCR duplicate or not
	confi_reads_dup = 0
	confi_reads_nondup = 0

	#Reverse or forward
	confi_reads_rev = 0
	confi_reads_fwd = 0

	# with error-corrected CB or UMI
	confi_CB = 0
	confi_UB = 0

	exon_reads = 0
	intron_reads = 0
	intergenic_reads = 0
	other_reads1 = 0

	#sense or antisense
	sense_reads = 0
	anti_reads = 0
	other_reads2 =0

	chrM_reads = 0
	#read match type
	read_type = collections.defaultdict(int)

	for chr_id, chr_len in chrom_info:
		logging.info("Processing \"%s\" ..." % chr_id)
		chrom_count = 0
		chrom_total_reads = set() #total reads in BAM file
		chrom_confi_reads = set() #reads marked as confidently mapped to transcriptome by xf:i:1 tag

		ALL = open(chr_id + '.all_reads_id.txt','w')
		CONF = open(chr_id + '.confident_reads_id.txt','w')
		for aligned_read in samfile.fetch(chr_id):
			total_alignments += 1
			chrom_count += 1
			read_id = aligned_read.query_name
			tag_dict = dict(aligned_read.tags) #{'NM': 1, 'RG': 'L1'}
			cigar_str = list2str( aligned_read.cigar)
			chrom_total_reads.add(read_id)

			#confident alignments
			if xf_tag in tag_dict and tag_dict[xf_tag]& 0x1 != 0:

				if chr_id == chrM_id:
					chrM_reads += 1
				#with or without CB/UMI barcode
				if CB_tag in tag_dict:
					confi_CB += 1
				if UMI_tag in tag_dict:
					confi_UB += 1

				#duplicate or not
				if aligned_read.is_duplicate:
					confi_reads_dup += 1
				else:
					confi_reads_nondup += 1

				#forward or reverse
				if aligned_read.is_reverse:
					confi_reads_rev += 1
				else:
					confi_reads_fwd += 1

				#Single character indicating the region type of this alignment (E = exonic, N = intronic, I = intergenic).
				if RE_tag in tag_dict:
					if tag_dict[RE_tag] == "E":
						exon_reads += 1
					elif tag_dict[RE_tag] == "I":
						intron_reads += 1
				elif tag_dict[RE_tag] == "N":
					intergenic_reads += 1
				else:
					other_reads1 += 1

				#sense or antisense
				if TX_tag in tag_dict:
					sense_reads += 1
				elif AN_tag in tag_dict:
					anti_reads += 1
				else:
					other_reads2 += 1

				#map type
				cigar_str = list2str( aligned_read.cigar)
				tmp = read_match_type(cigar_str)
				read_type[tmp] += 1

				confi_alignments += 1
				chrom_confi_reads.add(read_id)

			if chrom_count % step_size == 0:
				print("%d alignments processed.\r" % chrom_count, end=' ', file=sys.stderr)
		logging.info("Processed %d alignments from \"%s\"" % (chrom_count, chr_id))
		for i in chrom_total_reads:
			print (i, file=ALL)
		for i in chrom_confi_reads:
			print (i, file=CONF)
		ALL.close()
		CONF.close()

	logging.info("Processing total %d alignments mapped to all chromosomes." % total_alignments)

	logging.info("Count total mapped reads ...")
	output1 = subprocess.check_output("awk '!a[$0]++' *.all_reads_id.txt |tee All_reads_uniqID.txt | wc -l", shell=True)

	logging.info("Count confidently mapped reads ...")
	output2 = subprocess.check_output("awk '!a[$0]++' *.confident_reads_id.txt |tee confident_reads_uniqID.txt | wc -l", shell=True)

	logging.info("Removing intermediate files ...")
	subprocess.run("rm -rf *.all_reads_id.txt", shell=True)
	subprocess.run("rm -rf *.confident_reads_id.txt", shell=True)

	total_reads_n = int(output1.decode('utf-8').strip())
	confi_reads_n = int(output2.decode('utf-8').strip())
	non_confi_reads = total_reads_n - confi_reads_n



	print ('')
	print("\nTotal_alignments: %d" % total_alignments)
	print ("└--Confident_alignments: %d" % confi_alignments)
	print ('')
	print ("Total_mapped_reads:\t%d" % total_reads_n)
	print ("|--Non_confidently_mapped_reads:\t%d\t(%.2f%%)" % (non_confi_reads, non_confi_reads*100.0/total_reads_n))
	print ("└--Confidently_mapped_reads:\t%d\t(%.2f%%)" % (confi_reads_n, confi_reads_n*100.0/total_reads_n))

	print ("   |--Reads_with_PCR_duplicates:\t%d\t(%.2f%%)" % (confi_reads_dup, confi_reads_dup*100.0/confi_reads_n))
	print ("   └--Reads_no_PCR_duplicates:\t%d\t(%.2f%%)" % (confi_reads_nondup, confi_reads_nondup*100.0/confi_reads_n))
	print ('')

	print ("   |--Reads_map_to_forward(Waston)_strand:\t%d\t(%.2f%%)" % (confi_reads_fwd, confi_reads_fwd*100.0/confi_reads_n))
	print ("   └--Reads_map_to_Reverse(Crick)_strand:\t%d\t(%.2f%%)" % (confi_reads_rev, confi_reads_rev*100.0/confi_reads_n))
	print ('')

	print ("   |--Reads_map_to_sense_strand:\t%d\t(%.2f%%)" % (sense_reads, sense_reads*100.0/confi_reads_n))
	print ("   └--Reads_map_to_antisense_strand:\t%d\t(%.2f%%)" % (anti_reads, anti_reads*100.0/confi_reads_n))
	print ("   └--Other:\t%d\t(%.2f%%)" % (other_reads2, other_reads2*100.0/confi_reads_n))
	print ('')

	print ("   |--Reads_map_to_exons:\t%d\t(%.2f%%)" % (exon_reads, exon_reads*100.0/confi_reads_n))
	print ("   └--Reads_map_to_introns:\t%d\t(%.2f%%)" % (intron_reads, intron_reads*100.0/confi_reads_n))
	print ("   └--Reads_map_to_intergenic:\t%d\t(%.2f%%)" % (intergenic_reads, intergenic_reads*100.0/confi_reads_n))
	print ("   └--Other:\t%d\t(%.2f%%)" % (other_reads1, other_reads1*100.0/confi_reads_n))
	print ('')

	print ("   |--Reads_with_error-corrected_barcode:\t%d\t(%.2f%%)" % (confi_CB, confi_CB*100.0/confi_reads_n))
	print ("   └--Reads_no_error-corrected_barcode:\t%d\t(%.2f%%)" % ((confi_reads_n - confi_CB), (confi_reads_n - confi_CB)*100.0/confi_reads_n))
	print ('')
	print ("   |--Reads_with_error-corrected_UMI:\t%d\t(%.2f%%)" % (confi_UB, confi_UB*100.0/confi_reads_n))
	print ("   └--Reads_no_error-corrected_UMI:\t%d\t(%.2f%%)" % ((confi_reads_n - confi_UB), (confi_reads_n - confi_UB)*100.0/confi_reads_n))
	print ('')
	print ("   |--Reads_map_to_mitochonrial_genome:\t%d\t(%.2f%%)" % (chrM_reads, chrM_reads*100.0/confi_reads_n))
	print ("   └--Reads_map_to_nuclear_genome:\t%d\t(%.2f%%)" % ((confi_reads_n - chrM_reads), (confi_reads_n - chrM_reads)*100.0/confi_reads_n))

	print ('')

	for i in sorted(read_type):
		if i == 'Others':continue
		print ("   |--%s:\t%d\t(%.2f%%)" % (i, read_type[i], read_type[i]*100.0/confi_reads_n))
	print ("   └--%s:\t%d\t(%.2f%%)" % ('Others', read_type['Others'], read_type['Others']*100.0/confi_reads_n))
	print ('')



def readCount(infile, outfile, step_size=10000, limit=1000000, csv_out=False):
	'''
	Save reads that confidenlty mapped to transcriptome to a BAM file.

	Parameters
	----------
	infile : str
		Input BAM file. Must be sorted and indexed.
	outfile : str
		Prefix of output files.
	step_size: int
		Output progress report when step_size alignments have been processed.
	limit : int
		Only process this number of alignments and stop.
	'''
	logging.info("Reading BAM file \"%s\" ..." % infile)
	samfile = pysam.AlignmentFile(infile,'rb')

	#OUT = open(outfile,'w')
	total_alignments = 0
	CB_GN_READ = collections.defaultdict(dict)
	#CB_GN_UMI = collections.defaultdict(list)
	try:
		while(1):
			aligned_read = next(samfile)
			read_id = aligned_read.query_name
			tag_dict = dict(aligned_read.tags) #{'NM': 1, 'RG': 'L1'}
			if 'xf' in tag_dict and tag_dict['xf']& 0x1 == 0:
				continue
			if 'CB' in tag_dict:
				CB = tag_dict['CB'].replace('-1','')
			else:
				logging.debug('%s has no cell barcode!' % read_id)
			#if 'UB' in tag_dict:
			#	UMI = tag_dict['UB'].replace('-1','')
			if 'GX' in tag_dict:
				geneID = tag_dict['GX']
			else:
				geneID = 'NA'
			if 'GN' in tag_dict:
				geneSymbol = tag_dict['GN']
			else:
				geneSymbol="NA"

			if geneID == 'NA':
				continue

			try:
				CB_GN_READ[CB][geneID + '|' + geneSymbol] += 1
			except KeyError:
				CB_GN_READ[CB][geneID + '|' + geneSymbol] = 1

			#CB_GN_UMI[CB + ':' + geneID].append(UMI)

			total_alignments += 1
			if total_alignments % step_size == 0:
				print("%d alignments processed.\r" % total_alignments, end=' ', file=sys.stderr)
			if limit is not None:
				if total_alignments >= limit:
					break

	except StopIteration:
		pass
	"""
	logging.info('Total %d alignments processed' % total_alignments)
	logging.info('Convert dict of dict to pandas data frame')
	count_mat = pd.DataFrame(CB_GN_READ)
	(gene_num, CB_num) = count_mat.shape
	logging.info('Data frame contains %d genes (rows) and %d cell barcodes (columns)' % (gene_num, CB_num))

	if csv_out is True:
		logging.info('Save data frame to \"%s\"' % (outfile + '.csv'))
		count_mat.to_csv(outfile + '.csv', index=True, index_label="Gene")


	logging.info('Save data frame to \"%s\"' % (outfile + '.h5'))
	count_mat.to_hdf(outfile + '.h5',key='df',mode='w',complevel=6)
	"""


def CBC_UMIcount(infile, outfile, step_size=50000, CB_tag = 'CB', UMI_tag = 'UB', gene_tag = 'GN', CB_num = 100000, read_num=200, UMI_num = 200, gene_num = 200):
	'''
	Calculate UMI count for each cell barcode.

	Parameters
	----------
	infile : str
		Input BAM file. Must be sorted and indexed.
	outfile : str
		Prefix of output files.
	step_size: int
		Print progress report when step_size alignments have been processed.
	CB_tag : str
		Tag for cell barcode in the input BAM file.
	UMI_tag : str
		Tag for UMI in the input BAM file.
	CB_num : int
		Number of cell barcodes to be considered.
	read_num : int
		Read number threshold to filter invalid cell barcode.
		Cell barcode with reads less than this value will be skipped.
	'''


	CB_read_freq = collections.defaultdict(int) # cell_barcode:read_frequency
	CB_UMI_freq = {} # cell_barcode:UMI
	CB_gene_freq = {} # cell_barcode:gene
	CB_cutoff = CB_num
	logging.info("Top %d cell barcodes (ranked by associated UMI frequency) will be analyzed." % CB_cutoff)
	read_cutoff = read_num
	logging.info("Only count UMIs for cell barcodes with more than %d reads." % read_cutoff)

	#count read number for each cell barcode
	logging.info("Reading BAM file \"%s\". Count reads for each cell barcode ..." % infile)
	samfile = pysam.AlignmentFile(infile,'rb')
	total_alignments = 0
	try:
		while(1):
			aligned_read = next(samfile)
			tag_dict = dict(aligned_read.tags) #{'NM': 1, 'RG': 'L1'}
			if 'xf' in tag_dict and tag_dict['xf']& 0x1 == 0:
				continue

			#count reads
			if CB_tag in tag_dict:
				CB = tag_dict[CB_tag].replace('-1','')
			else:
				continue
			CB_read_freq[CB] += 1
			total_alignments += 1
			if total_alignments % step_size == 0:
				print("%d alignments processed.\r" % total_alignments, end=' ', file=sys.stderr)
	except StopIteration:
		pass
	logging.info('Total %d alignments processed' % total_alignments)

	logging.info("Filtering cell barcodes ...")
	CB_usable = set()
	for k,v in CB_read_freq.items():
		if v >= read_cutoff:
			CB_usable.add(k)
	logging.info("Total cell barcode: %d" % len(CB_read_freq))
	logging.info("Cell barcode with more than %d reads: %d" % (read_cutoff, len(CB_usable)))


	#count UMI for each cell barcode
	logging.info("Reading BAM file \"%s\". Count UMIs for each cell barcode ..." % infile)
	samfile = pysam.AlignmentFile(infile,'rb')
	CB_freq_list = {} #UMI count for each barcode
	gene_freq = {} # gene count for each barcode
	total_alignments = 0
	try:
		while(1):
			aligned_read = next(samfile)
			read_id = aligned_read.query_name
			#chrom = samfile.get_reference_name(aligned_read.reference_id)
			#if aligned_read.is_duplicate:continue
			tag_dict = dict(aligned_read.tags) #{'NM': 1, 'RG': 'L1'}
			if 'xf' in tag_dict and tag_dict['xf']& 0x1 == 0:
				continue

			if CB_tag in tag_dict:
				CB = tag_dict[CB_tag].replace('-1','')
				if CB not in CB_usable: #filter out CB with low number of reads
					continue
			else:
				logging.debug('%s has no cell barcode!' % read_id)
				continue

			# count UMI frequency
			if UMI_tag in tag_dict:
				UMI = tag_dict[UMI_tag].replace('-1','')
				if 	CB not in CB_freq_list:
					CB_freq_list[CB] = set(UMI)
				else:
					CB_freq_list[CB].add(UMI)
			else:
				logging.debug('%s has no UMI!' % read_id)
				continue

			# count gene frequency
			if gene_tag in tag_dict:
				gene_names = tag_dict[gene_tag].split(';')
				if CB not in gene_freq:
					gene_freq[CB] = set()
				else:
					for gene_name in gene_names:
						gene_freq[CB].add(gene_name)
			else:
				continue

			total_alignments += 1
			if total_alignments % step_size == 0:
				print("%d alignments processed.\r" % total_alignments, end=' ', file=sys.stderr)
	except StopIteration:
		pass
	logging.info('Total %d alignments processed' % total_alignments)

	for k,v in CB_freq_list.items():
		temp = len(v)
		CB_UMI_freq[k] = temp

	for k,v in gene_freq.items():
		temp = len(v)
		CB_gene_freq[k] = temp

	OUT = open(outfile + '.Read_UMI_freq.tsv','w')
	logging.info ("Writing cell barcodes' reads and UMI frequencies to \"%s\"" % (outfile + '.Read_UMI_freq.tsv'))
	print ('\t'.join(['Serial','Cell_barcode', 'Read_count', 'UMI_count', 'Gene_count']), file=OUT) #do NOT change the header
	count = 0
	for k in sorted(CB_UMI_freq, key=CB_UMI_freq.get, reverse=True):
		count += 1
		if CB_UMI_freq[k] < UMI_num:
			continue
		if k in CB_gene_freq:
			if CB_gene_freq[k] < gene_num:
				continue
			else:
				print ('\t'.join([str(count), k, str(CB_read_freq[k]), str(CB_UMI_freq[k]), str(CB_gene_freq[k])]), file=OUT)
		else:
			continue
		if count > CB_cutoff:
			break
	OUT.close()
	logging.info ("Done.")
