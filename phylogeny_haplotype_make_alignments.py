import argparse
import os
import pandas as pd
import re
import subprocess
import glob

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

        # miss
        parser.add_argument(
                '--miss',
                type=float,
                default=0.5,
                help='Do not include sequence if more than this value missing (0<1).'
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

	# folders with haplotypes
	d = pd.read_csv(args.file)

	haps = {}
	nonhaps = {}
	lineages = d.lineage.unique().tolist()

	for lineage in lineages:
		hapdir = os.path.join(args.dir, 'variants', lineage)
		if not os.path.isdir(hapdir):
			nonhaps[lineage] = os.path.join(args.dir, 'PRG', '%s.fasta' % lineage)			
		else:
			haps[lineage] = hapdir

	return outdir, subdir, haps, nonhaps
	

def get_seq(haps, nonhaps, args):
	seqs = {}
	ids = {}

	for lin in haps:
		seqfiles = glob.glob(haps[lin] + '/*')
		for seqfile in seqfiles:

                	s = open(seqfile, 'r')
                	id = ''
			seqname = re.sub('^.*/', '', seqfile)
			seqname = re.sub('\..*', '', seqname)
	
			if seqname not in seqs:
				seqs[seqname] = {}

                	for l in s:
                        	if re.search('>', l):
                                	id = re.search('>(\S+)', l.rstrip()).group(1)
                                	seqs[seqname][id] = ''
					if id not in ids:
						ids[id] = 0
					ids[id] += 1
                        	else:
                                	seqs[seqname][id] += l.rstrip()
			
			s.close()

	# get the non haplo sequences in there
	for lin in nonhaps:
		seqfile = nonhaps[lin]
		s = open(seqfile, 'r')
		
		for l in s:
			if re.search('>', l):
				seqname = re.search('>(\S+)', l.rstrip()).group(1)
				if seqname not in seqs:
					seqs[seqname] = {}
				seqs[seqname][lin] = ''
				if lin not in ids:
					ids[id] = 0
				ids[id] += 1
			else:
				seqs[seqname][lin] += l.rstrip()
		s.close()

	# filter
	for locus in seqs:
		for ind in ids:
			if ind in seqs[locus]:
				seq = seqs[locus][ind]
				per_miss = (seq.count('N') + seq.count('n')) / float(len(seq))

				if per_miss > args.miss:
					del seqs[locus][ind]

			
	return seqs, ids


def print_loci(dir, subdir, seq, inds):
	loci = sorted(seq.keys())
	inds = sorted(inds)
	n_sp = len(inds)

	d = os.path.join(dir, 'locus_data.csv')
	d = open(d, 'w')
	d.write('locus,n_lineages,missingness,length,PICs\n')

	for locus in loci:
		# count how many sps have the locus
		count = len(seq[locus])

		# only print out the locus if in 2 sp
		if count >= 2:
			out = os.path.join(subdir, '%s.fasta' % locus)
			o = open(out, 'w')

			for ind in seq[locus]:
				o.write('>%s\n%s\n' % (ind, seq[locus][ind]))
			o.close()

		d.write('%s,%s,%.3f,NA,NA\n' % (locus, count, count / float(n_sp)))

	d.close()


def main():
	# get arguments
	args = get_args()
	# get genome files, make dirs
	dir, subdir, haps, nonhaps = get_files(args)
	# get sequences
	seq, inds = get_seq(haps, nonhaps, args)
	# print the loci
	print_loci(dir, subdir, seq, inds)	


if __name__ == "__main__":
	main()
