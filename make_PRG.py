import argparse
import os
import pandas as pd
import re
import subprocess

"""
Sonal Singhal
created on 21 June 2016
"""

def get_args():
	parser = argparse.ArgumentParser(
		description="Make pseudoreference genomes per lineage.",
        	formatter_class=argparse.ArgumentDefaultsHelpFormatter
		)

	# lineage
	parser.add_argument(
                '--lineage',
                type=str,
                default=None,
                help='Lineage for which to make PRG.'
                )

	# sample file
	parser.add_argument(
                '--file',
                type=str,
                default=None,
                help='File with sample info.'
                )
	
	# base dir
	parser.add_argument(
		'--dir',
		type=str,
		default=None,
		help="Base directory as necessary"
		     " when used with pipeline"
                )

	# loc of matches files
	parser.add_argument(
                '--mdir',
                type=str,
                default=None,
                help='Dir with match data, only '
                     ' necessary if not used with pipeline'
                )

	# loc of assembly files
	parser.add_argument(
                '--adir',
                type=str,
                default=None,
                help='Dir with assembly files, '
                     ' only necessary if not used '
                     ' with pipeline'
                )

	# output dir
	parser.add_argument(
                '--outdir',
                type=str,
                default=None,
                help='Directory to output PRGs, '
                     'only necessary if not used '
                     'with pipeline'
                )

        # matches to keep
        parser.add_argument(
                '--keep',
                type=str,
                default=None,
                help='Tags to keep; comma-delimited. Options are '
                     ' easy_recip_match, complicated_recip_match '
                )

	return parser.parse_args()


def get_samples(args):
	d = pd.read_csv(args.file)
	return d.ix[d['lineage'] == args.lineage, 'sample'].tolist()


def get_sequences(args, samples):
	if args.adir:
		adir = args.adir
	else:
		adir = os.path.join(args.dir, 'trinity_assembly')

	seqs = {}
	for sample in samples:
		seqs[sample] = {}

		seq = os.path.join(adir, '%s.fasta' % sample)
		s = open(seq, 'r')
		id = ''

		for l in s:
			if re.search('>', l):
				id = re.search('>(\S+)', l.rstrip()).group(1)
				seqs[sample][id] = ''
			else:
				seqs[sample][id] += l.rstrip()
		s.close()

	return seqs


def rev_comp(seq):
	seq = seq.upper()
	seq_dict = {'A':'T','T':'A','G':'C','C':'G'}
	return "".join([seq_dict[base] for base in reversed(seq)])


def output(args, samples, seq):
	if args.mdir:
		mdir = args.mdir
	else:
		mdir = os.path.join(args.dir, 'matches')

	keep = re.split(',', args.keep)

	match = {}
	for sample in samples:
		m = os.path.join(mdir, '%s_matches.csv' % sample)
		f = open(m, 'r')
		head = f.next()
		for l in f:
			d = re.split(',', l.rstrip())
			if d[5] in keep:
				c = d[1]
				if c in match:
					# only keep match if evalue diff not too big
					# and it is the longer contig	
					if float(d[6]) / match[c]['eval'] < 1e3:
						if len(seq[sample][d[0]]) > match[c]['len']:
							match[c] = {'sample': sample, 'con': d[0],
                                        		            'len': len(seq[sample][d[0]]),
                                        		            'eval': float(d[6]), 'orr': d[4]}
				else:
					match[c] = {'sample': sample, 'con': d[0],
						    'len': len(seq[sample][d[0]]),
                                	             'eval': float(d[6]), 'orr': d[4]}
		f.close()

	
	if not args.outdir:
		outdir = os.path.join(args.dir, 'PRG')
	else:
		outdir = args.outdir

	if not os.path.isdir(outdir):
		os.mkdir(outdir)

	out = os.path.join(outdir, '%s.fasta' % args.lineage)
	o = open(out, 'w')
	for c in match:
		s = seq[match[c]['sample']][match[c]['con']]
		if match[c]['orr'] == '-':
			s = rev_comp(s)
		o.write('>%s %s\n%s\n' % (c, match[c]['con'], s))
	o.close()


def main():
	# get arguments
	args = get_args()
	# get samples
	samples = get_samples(args)
	# get sequences
	seq = get_sequences(args, samples)
	# output PRG
	output(args, samples, seq)


if __name__ == "__main__":
	main()
