#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug  2 15:43:45 2020

@author: m102324
"""
import sys
import collections
import pandas as pd
import logomaker
import matplotlib.pyplot as plt

import logging
from qcmodule import ireader


def fasta_iter(infile):
	"""
	Generate seq or qual string.

	Parameters
	----------
	infile : str
		Input fasta file.

	Yields
	------
	str
		String of nucleotides or quality scores.
	"""
	logging.info("Reading FASTA file \"%s\" ..." % infile)
	for l in ireader.reader(infile):
		l = l.strip()
		if len(l) == 0:continue
		if l.startswith('>'):continue
		yield l

def fastq_iter(infile, mode = 'seq'):
	"""
	Generate seq or qual string.

	Parameters
	----------
	infile : str
		Input fastq file.
	mode : str
		Must be 'seq' or 'qual'.

	Yields
	------
	str
		String of nucleotides or quality scores.
	"""
	logging.info("Reading FASTQ file \"%s\" ..." % infile)
	count = 0
	s = ''
	q = ''
	for l in ireader.reader(infile):
		l = l.strip()
		if len(l) == 0:continue
		count += 1
		if count == 2:
			s = l
		if count == 4:
			q = l
			if mode == 'seq':
				yield s
			elif mode == 'qual':
				yield q
			count = 0

def qual2countMat(q_obj,limit, step_size=100000):
	"""
	Generate count data frame.

	Parameters
	----------
	q_obj : iterable
		Generator returned by fastq_iter
	limit : int
		Only read this many sequences.
	step_size : int
		Output progress report every step_size.

	Return
	------
	Data frame
	"""
	dat = collections.defaultdict(dict)
	count = 0
	for qstr in q_obj:
		count += 1
		for i,q in enumerate(qstr):
			q_score = ord(q) - 33
			if q_score not in dat[i]:
				dat[i][q_score] = 1
			else:
				dat[i][q_score] += 1
		if count % step_size == 0:
			print("%d quality sequences finished\r" % count, end=' ', file=sys.stderr)
		if limit is not None:
			if count >= limit:
				break

	logging.info("%d quality sequences finished" % count)

	logging.info("Make data frame from dict of dict ...")
	qual_mat = pd.DataFrame.from_dict(dat)

	logging.info("Filling NA as zero ...")
	qual_mat = qual_mat.fillna(0)
	qual_mat.index.name='pos'
	return qual_mat.T

def seq2countMat(s_obj, limit, step_size=100000,exclude_N = False):
	"""
	Generate count data frame.

	Parameters
	----------
	s_obj : iterable
		Generator returned by fastq_iter or fasta_iter.
	limit : int
		Only read this many sequences.
	step_size : int
		Output progress report when step_size sequences have been processed.
	exclude_N : bool
		If True, sequences containing "N" will be skipped.
	Return
	------
	Data frame
	"""
	mat = collections.defaultdict(dict)
	count = 0
	for s in s_obj:
		count += 1
		if (exclude_N is True) and "N" in s:
			continue
		for indx,base in enumerate(s):
			if base not in mat[indx]:
				mat[indx][base] = 1
			else:
				mat[indx][base] += 1
		if count % step_size == 0:
			print("%d sequences finished\r" % count, end=' ', file=sys.stderr)
		if limit is not None:
			if count >= limit:
				break
	logging.info("%d sequences finished" % count)

	logging.info("Make data frame from dict of dict ...")
	count_mat = pd.DataFrame.from_dict(mat)

	logging.info("Filling NA as zero ...")
	count_mat = count_mat.fillna(0)
	count_mat.index.name='pos'
	return count_mat.T

def make_logo(mat, outfile, exclude_N=False, font_name='sans', stack_order='big_on_top', flip_below=True, shade_below=0.0, fade_below=0.0, highlight_start = 0, highlight_end = 15, oformat='pdf'):
	"""
	Call logomaker (https://github.com/jbkinney/logomaker) to make nucleotide logo.

	Parameters
	----------
	mat : data frame
		Nucleotide frequencies table.
	exclude_N : bool
		If True, sequences containing "N" will be skipped.
	font_name : str
		The character font to use when rendering the logo.
		run logomaker.list_font_names() to get a list of
		valid font names.
	stack_order : str
		Must be 'big_on_top', 'small_on_top', or 'fixed'.
		'big_on_top' : nucleotide with the highest frequency will be on top;
		'small_on_top' : nucleotide with the lowest frequency will be on top;
		'fixed' : nucleotides from top to bottom are in the same order as characters appear in the data frame.
	flip_below : bool
		If True, characters below the x-axis (which correspond to negative
		values in the matrix) will be flipped upside down
	shade_below : float
		The amount of shading to use for characters drawn below the x-axis.
		Must be in [0.0, 1.0]
	fade_below : float
		The amount of fading to use for characters drawn below the x-axis.
		Must be in [0.0, 1.0]
	highlight_start : int
		Highlight logo from this position. Must be within [0, len(logo)-1].
	highlight_end : int
		Highlight logo to this position. Must be within [0, len(logo)-1].
	"""
	logging.info("Making logo ...")
	logging.debug("Font name is: %s" % font_name)
	logging.debug("Stack order is: %s" % stack_order)
	logging.debug("Flip below flag: %s" % flip_below)
	logging.debug("Shade below score: %f" % shade_below)
	logging.debug("Fade below score: %f" % fade_below)
	if exclude_N:
		logging.info("'N' will be excluded.")
		color = {'A':'green','C':'blue','G':'orange','T':'red'}
	else:
		logging.info("'N' will be kept.")
		color = {'A':'green','C':'blue','G':'orange','T':'red','N':'grey'}
	if (highlight_start is not None) and (highlight_end is not None):
		if highlight_start < 0 or highlight_end < 0 or highlight_start > highlight_end:
			logging.error("Incorrect highlight positions")
			sys.exit(0)
	logging.info("Mean-centered logo saved to \"%s\"." % (outfile + '.logo_mean_centered.' + oformat))
	logo = logomaker.Logo(mat,center_values=True,color_scheme=color, font_name = font_name, stack_order = stack_order, flip_below = flip_below, shade_below = shade_below, fade_below = fade_below)
	if isinstance(highlight_start, int) and isinstance(highlight_end, int):
		logging.info("Highlight logo from %d to %d" % (highlight_start, highlight_end))
		logo.highlight_position_range(pmin = highlight_start, pmax = highlight_end)
	plt.savefig(outfile + '.logo.mean_centered.%s' % oformat.lower())

	logging.info("Logo saved to \"%s\"." % (outfile + '.logo.' + oformat))
	logo = logomaker.Logo(mat,center_values=False,color_scheme=color, font_name = font_name, stack_order = stack_order, flip_below = flip_below, shade_below = shade_below, fade_below = fade_below)
	if isinstance(highlight_start, int) and isinstance(highlight_end, int):
		logging.info("Highlight logo from %d to %d" % (highlight_start, highlight_end))
		logo.highlight_position_range(pmin = highlight_start, pmax = highlight_end)
	plt.savefig(outfile + '.logo.%s' % oformat.lower())

