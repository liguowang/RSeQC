#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  5 15:24:22 2020

@author: m102324
"""
import logging
import subprocess

def make_heatmap(infile, outfile, filetype, cell_width, cell_height, col_angle, font_size, text_color, no_numbers = False,log2_scale=False):
	"""
	Use pheatmap to generate heatmap.

	Parameters
	----------
	infile : str
		Input file
	outfile : str
		Output file prefix.
	filetype : str
		Filetype is decided by the extension. Currently following formats are supported: png, pdf, tiff, bmp, jpeg.
	cell_width : int
		Points of cell width in the heatmap.
	cell_height : int
		Points of cell height in the heatmap.
	cell_angle : int
		The angle (must be 0, 45, 90, 270, 315) of column text lables under the heatmap.
	font_size : int
		Font size. If --display-num was set, the font size of number = 0.8 * font_size.
	text_color : str
		The color of numbers in each cell.
	dispaly_num : bool
		If set, numbers will be displayed on each cell
	"""
	logging.info("Writing R code to \"%s\"" % (outfile + '.r'))
	ROUT = open(outfile + '.r','w')
	print ("if(!require(pheatmap)){install.packages(\"pheatmap\")}", file=ROUT)
	print ('library(pheatmap)', file=ROUT)
	print ("dat = read.table(file = '%s',sep=',',header=T,row.names=1,check.names=FALSE)" % infile, file=ROUT)

	if no_numbers:
		logging.info("Does not displayed numerical values on heatmap")
		if log2_scale:
			print("pheatmap(log2(as.matrix(dat)+1), filename='%s', cellwidth = %d, cellheight = %d, display_numbers = FALSE, angle_col=%d, fontsize=%d,cluster_rows=F, cluster_cols=F, scale='none',color = colorRampPalette(c('#f1eef6','#d7b5d8','#df65b0', '#ce1256'))(50))" % ((outfile + '.' + filetype), cell_width, cell_height, col_angle, font_size), file=ROUT)
		else:
			print("pheatmap(as.matrix(dat), filename='%s', cellwidth = %d, cellheight = %d, display_numbers = FALSE, angle_col=%d, fontsize=%d,cluster_rows=F, cluster_cols=F, scale='none',color = colorRampPalette(c('#f1eef6','#d7b5d8','#df65b0', '#ce1256'))(50))" % ((outfile + '.' + filetype), cell_width, cell_height, col_angle, font_size), file=ROUT)
	else:
		logging.info("Displayed numerical values on heatmap")
		if log2_scale:
			logging.info("Numbers will be displayed on log2 scale")
			print("pheatmap(log2(as.matrix(dat)+1), filename='%s', cellwidth = %d, cellheight = %d, display_numbers = TRUE, angle_col=%d, fontsize=%d, number_color=\'%s\',cluster_rows=F, cluster_cols=F, scale='none',color = colorRampPalette(c('#f1eef6','#d7b5d8','#df65b0', '#ce1256'))(50))" % ((outfile + '.' + filetype), cell_width, cell_height, col_angle, font_size, text_color), file=ROUT)
		else:
			print("pheatmap(as.matrix(dat), filename='%s', cellwidth = %d, cellheight = %d, display_numbers = TRUE, angle_col=%d, fontsize=%d, number_color=\'%s\',cluster_rows=F, cluster_cols=F, scale='none',color = colorRampPalette(c('#f1eef6','#d7b5d8','#df65b0', '#ce1256'))(50))" % ((outfile + '.' + filetype), cell_width, cell_height, col_angle, font_size, text_color), file=ROUT)


	ROUT.close()

	logging.info("Running R script file \"%s\"" % (outfile + '.r'))
	try:
		subprocess.call("Rscript " + outfile + '.r', shell=True)
	except:
		logging.error("Failed to run Rscript file \"%s\"" + (outfile + '.r'))
		pass
