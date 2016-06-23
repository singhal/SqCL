import argparse
import os
import pandas as pd
import re
import subprocess
import zipfile

"""
Sonal Singhal
created on 21 June 2016
Written assuming blat
"""

def get_args():
	parser = argparse.ArgumentParser(
		description="Match contigs to their probes. "
                            "Written assuming blat v36.",
        	formatter_class=argparse.ArgumentDefaultsHelpFormatter
		)

	# blat
	parser.add_argument(
		'--blat',
		type=str,
		default=None,
		help='Full path to blat executable'
		)
	
	# sample
	parser.add_argument(
                '--sample',
                type=str,
                default=None,
                help='Sample for which to run script.'
                )

	# dir
	parser.add_argument(
                '--dir',
                type=str,
                default=None,
                help='Baseline directory if running as pipeline .'
                )

	# evalue
	parser.add_argument(
		'--evalue',
		type=float,
		default=1e-10,
		help="Minimum evalue req'd for match, given in "
	             "Xe-X format."
		)

	# percent matching
	parser.add_argument(
		'--match',
		type=float,
		default=80,
		help="Percent matching req'd for match, given in "
		     "XX format."
		)

	# outdir
	parser.add_argument(
		'--outdir',
		type=str,
		default=None,
		help="Output directory for match info, "
                     " only define if using out of pipeline."
		)	

	# database
	parser.add_argument(
		'--db',
		type=str,
		default=None,
		help="Database to use to search."
		)

	# query
	parser.add_argument(
                '--query',
                type=str,
                default=None,
                help="Query to use to search, "
                     " only define if using out of pipeline."
                )

	return parser.parse_args()


def get_query(args):
	if args.query:
		query = args.query
	else:
		query = os.path.join(args.dir, 'trinity_assembly', 
	                             '%s.fasta' % args.sample)

	return query


def run_blat(args, query):
	if not args.outdir:
		outdir = os.path.join(os.dir, 'matches')
	else:
		outdir = args.outdir

	# make the outdir if it doesn't exist
	if not os.path.isdir(outdir):
		os.mkdir(outdir)
	# make subdir
	subdir = os.path.join(outdir, 'blat_results')
	if not os.path.isdir(subdir):
		os.mkdir(subdir)
	

	outfile1 = os.path.join(subdir, '%s_to_probes' % args.sample)
	subprocess.call("%s %s %s %s -out=blast8" % (args.blat, args.db,
		        query, outfile1), shell=True)

	outfile2 = os.path.join(subdir, 'probes_to_%s' % args.sample)
	subprocess.call("%s %s %s %s -out=blast8" % (args.blat, query,
			args.db, outfile2), shell=True)
	
	return outdir, outfile1, outfile2


def sub_parse_blat(args, out, regex):
	matches = {}

        f = open(out, 'r')
        for l in f:
                d = re.split('\s+', l.rstrip())
                d[regex] = re.sub('_p\d+', '', d[regex])
        	
		# check orientation
		orr = '+'
		if int(d[9]) < int(d[8]):
			orr = '-'
	
	        res = {'match': d[1], 'per': float(d[2]), 'length': int(d[3]),
                       'eval': float(d[10]), 'status': None, 'orr': orr}
		c = d[0]

                # only keep these matches
                if res['eval'] <= args.evalue and res['per'] >= args.match:
                        if c in matches:
				exist = None
				for ix, hash in enumerate(matches[c]):
					if hash['match'] == res['match']:
						exist = ix
				if exist is not None:
					# update hash if better match
					if res['eval'] < matches[c][exist]['eval']:
						matches[c][exist] = res
				else:
					# only keep if within 10 orders a match
                                        mineval = min([x['eval'] for x in matches[c]])

					if res['eval'] / mineval < 1e10:
						matches[c].append(res)
                        else:
                                matches[c] = []
                                matches[c].append(res)
        f.close()

	return matches


def parse_blat(args, dir, query, out1, out2):
	# get contig lengths
	s = open(query, 'r')
	c_len = {}
	id = ''
	for l in s:
		if re.search('>', l):
			id = re.search('>(\S+)', l.rstrip()).group(1)
			c_len[id] = ''
		else:
			c_len[id] += l.rstrip()
	s.close()
	
	for id, s in c_len.items():
		c_len[id] = len(s)

	m1 = sub_parse_blat(args, out1, 1)
	m2 = sub_parse_blat(args, out2, 0)

	for c in m1:
		top1 = m1[c][0]['match']

		if top1 in m2:
			if len(m2[top1]) == 1 and len(m1[c]) == 1:
				if c == m2[top1][0]['match']:
					# yay, recip
					m1[c][0]['status'] = 'easy_recip_match'
				else:
					m1[c][0]['status'] = 'ditched_no_recip_match'
			else:
				if len(m1[c]) == 1:
					mineval = min([x['eval'] for x in m2[top1]])
					contigs = [x['match'] for x in m2[top1] if x['eval'] / mineval < 1e2]
					lengths = dict([(x, c_len[x]) for x in contigs])
					winner = max(lengths, key=lengths.get)
				
					if winner == c:
						m1[c][0]['status'] = 'complicated_recip_match'
					else:
						m1[c][0]['status'] = 'ditched_no_recip_match'
				else:
					m1[c][0]['status'] = 'ditched_too_many_matches'
		else:
			m1[c][0]['status'] = 'ditched_no_match'

	outfile = os.path.join(dir, '%s_matches.csv' % args.sample)
	o = open(outfile, 'w')

	keys = ['match', 'per', 'length', 'orr', 'status', 'eval']
	o.write('contig,%s\n' % (','.join(keys)))
	for c in m1:
		o.write('%s,%s\n' % (c, ','.join([str(m1[c][0][x]) for x in keys])))
	o.close()


def main():
	# get arguments
	args = get_args()
	# get query seq
	query = get_query(args)
	# run blat
	dir, out1, out2 = run_blat(args, query)
	# parse blat
	parse_blat(args, dir, query, out1, out2)


if __name__ == "__main__":
	main()
