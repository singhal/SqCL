import argparse
import glob
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

        # miss
        parser.add_argument(
                '--miss',
                type=float,
                default=None,
                help='How much missing data will you tolerate?'
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


def make_concatenated(args, outdir, sps, loci):
	subdir = os.path.join(outdir, 'concatenated')
	if not os.path.isdir(subdir):
		os.mkdir(subdir)

	# where the alignments are
	seqdir = os.path.join(outdir, 'alignments')

	file = os.path.join(subdir, 'concatenated%s.phy' % args.miss)

	seq = {}
	for sp in sps:
		seq[sp] = ''

	for locus in loci:
		f = os.path.join(seqdir, '%s.fasta.aln' % locus)
		f = open(f, 'r')
		id = ''
		s = {}
		for l in f:
			if re.search('>', l):
				id = re.search('>(\S+)', l).group(1)
				s[id] = ''
			else:
				s[id] += l.rstrip()
		f.close()

		null = '-' * len(s[s.keys()[0]])
		for sp in sps:
			if sp not in s:
				seq[sp] += null
			else:
				seq[sp] += s[sp]

	f = open(file, 'w')
	f.write(' %s %s\n' % (len(sps), len(seq[seq.keys()[0]])))
	for sp, s in seq.items():
		f.write('%s   %s\n' % (sp, s))
	f.close()


def main():
	args = get_args()
	outdir, sp, loci = get_sp_loci(args)
	make_concatenated(args, outdir, sp, loci)

if __name__ == "__main__":
	main()

