import argparse
import os
import pandas as pd
import re
import subprocess

"""
Sonal Singhal
created on 23 June 2016
Written assuming nothing!
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
		help='Output directory for alignments if not '
	             'running in context of pipeline.'
		)

	return parser.parse_args()


def get_files(args):
	# make all the directory structure
	if args.outdir:
		outdir = args.outdir
	else:
		outdir = os.path.join(args.dir, 'phylogeny')

	subdir = os.path.join(outdir, 'alignments')
	
	# main result folder
	if not os.path.isdir(outdir):
		os.mkdir(outdir)

	# alignment result folder
	if not os.path.isdir(subdir):
		os.mkdir(subdir)

	# get the genomes
	d = pd.read_csv(args.file)
	gen = d['genome'].unique().tolist()	
	gen = pd.isnull(gen)

	genomes = {}
	# genomes aren't defined
	# define them ourselves
	if True in gen:
		for l in d['lineage'].unique().tolist():
			g = os.path.join(args.dir, 'PRG', '%s.fasta' % l)
			genomes[l] = g 
	else:
		for l in d['lineage'].unique().tolist():
			genomes[l] = d.ix[d['lineage'] == l, 'genome'].unique().tolist()[0] 

	return outdir, subdir, genomes
	

def get_seq(genomes):
	seqs = {}
	ids = {}

	for lin, file in genomes.items():
		seqs[lin] = {}

		seq = os.path.join(file)
                s = open(seq, 'r')
                id = ''

                for l in s:
                        if re.search('>', l):
                                id = re.search('>(\S+)', l.rstrip()).group(1)
                                seqs[lin][id] = ''
				if id not in ids:
					ids[id] = 0
				ids[id] += 1
                        else:
                                seqs[lin][id] += l.rstrip()
                s.close()

	return seqs, ids


def print_loci(dir, subdir, seq, loci):
	sps = sorted(seq.keys())
	n_sp = len(sps)

	d = os.path.join(dir, 'locus_data.csv')
	d = open(d, 'w')
	d.write('locus,n_lineages,missingness,length,PICs\n')

	for locus in loci:
		# count how many sps have the locus
		count = sum([1 for sp in sps if locus in seq[sp]])

		# only print out the locus if in 4 sp
		if count >= 4:
			out = os.path.join(subdir, '%s.fasta' % locus)
			o = open(out, 'w')

			for sp in sps:
				if locus in seq[sp]:
					o.write('>%s\n%s\n' % (sp, seq[sp][locus]))
			o.close()

		d.write('%s,%s,%.3f,NA,NA\n' % (locus, count, count / float(n_sp)))

	d.close()


def main():
	# get arguments
	args = get_args()
	# get genome files, make dirs
	dir, subdir, gen = get_files(args)
	# get sequences
	seq, loci = get_seq(gen)
	# print the loci
	print_loci(dir, subdir, seq, loci)	


if __name__ == "__main__":
	main()
