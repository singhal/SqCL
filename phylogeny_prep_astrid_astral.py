from rpy2.robjects.packages import importr
import rpy2.robjects as ro

import argparse
import glob
import os
import multiprocessing as mp
import pandas as pd
import re
import subprocess
import random

"""
Sonal Singhal
created on 28 June 2016
Written assuming:
	* R
	* 'ape' in R
	* creates files that are to be used by ASTRID & ASTRAL
"""

def get_args():
	parser = argparse.ArgumentParser(
		description="This creates the files that then get " 
					"aligned in the next script.",
		formatter_class=argparse.ArgumentDefaultsHelpFormatter
		)

	# file
	parser.add_argument(
		'--file',
		type=str,
		default=None,
		help='File with information for phylogeny making.'
		)

	# miss
	parser.add_argument(
		'--miss',
		type=float,
		default=None,
		help='How much missing data will you tolerate?'
		)

	# tolerance for small branch lengths
	# most gene tree programs resolve all branches
	# artificially
	parser.add_argument(
		'--tol',
		type=float,
		default=1e-5,
		help='Collapse branch lengths shorter than this.'
		)

	# collapse low support nodes
	# this is a weird subjective exercise
	parser.add_argument(
		'--collapse',
		type=int,
		default=1,
		help='Collapse nodes with less support than this.'
		)

	# dir
	parser.add_argument(
		'--dir',
		type=str,
		default=None,
		help='Base directory when used in context of '
			 'pipeline.'
		)

	# output dir
	parser.add_argument(
		'--outdir',
		type=str,
		default=None,
		help='Output directory for phylogeny if not '
			 'running in context of pipeline.'
		)

	return parser.parse_args()


def get_sp_loci(args):
	d = pd.read_csv(args.file)
	sp = d['lineage'].unique().tolist()

	if args.dir:
		outdir = os.path.join(args.dir, 'phylogeny')
	else:
		outdir = args.outdir

	loc_file = os.path.join(outdir, 'locus_data.csv')
	d = pd.read_csv(loc_file)

	loci = d.ix[d.missingness >= args.miss, 'locus'].tolist()

	return outdir, sp, loci


def manipulate_gene_tree(ape, phangorn, tree, tol, collapse):
	# polytomize any weak nodes
	tree = ape.di2multi(tree, tol=tol)

	# collapse any low support branches
	# first convert node labels (bootstraps) to numerics
	asnumeric = ro.r('as.numeric')
	tree[4] = asnumeric(tree[4])
	tree = phangorn.pruneTree(tree, collapse)

	return tree


def create_files(args, dir, loci):
	# where the gene trees are
	subdir = os.path.join(dir, 'gene_trees')
	if not os.path.isdir(subdir):
		os.mkdir(subdir)	

	# where to put the files
	outdir = os.path.join(dir, 'astrid_astral')
	if not os.path.isdir(outdir):
		os.mkdir(outdir)

	# start up r
	ape = importr('ape')
	phangorn = importr('phangorn')

	out = os.path.join(outdir, 'best_trees_miss%s_tol%s_collapse%s.trees' % 
						  (args.miss, args.tol, args.collapse))
	bs = os.path.join(outdir, 'bootstrap_files_miss%s_tol%s_collapse%s.txt' % 
						 (args.miss, args.tol, args.collapse))
	bs_out = open(bs, 'w')

	for locus in loci:
		tree = os.path.join(subdir, '%s.bestTree.tre' % locus)
		if os.path.isfile(tree):
			# deal with best tree
			a = ape.read_tree(tree)
			a = manipulate_gene_tree(ape, phangorn, a, args.tol, args.collapse)
			ape.write_tree(a, file=out, append=True)
	
			bs = os.path.join(subdir, '%s.bootstrap.trees' % locus)
			bs_out.write('%s\n' % bs)

	bs_out.close()


def main():
	args = get_args()
	dir, sp, loci = get_sp_loci(args)
	create_files(args, dir, loci)	

if __name__ == "__main__":
	main()
