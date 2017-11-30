import argparse
import gzip
import os
import pandas as pd
import re
import subprocess
import random
import numpy as np
import math

"""
Sonal Singhal
created on 4 Feb 2017
"""

def get_args():
	parser = argparse.ArgumentParser(
		description="Get D estimates for LD.",
        	formatter_class=argparse.ArgumentDefaultsHelpFormatter
		)
 
	# outdir
	parser.add_argument(
                '--outdir',
                type=str,
                default=None,
                help='Output directory for alignments, only needed '
                     'if not running in context of pipeline.'
                )
                		   
        # vcf
        parser.add_argument(
                '--vcf',
                type=str,
                default=None,
                help="Full path to VCF if "
                     "you aren't running in context of pipeline."
                )
	
	return parser.parse_args()


def get_haplo(outdir, vcf):
	f = gzip.open(vcf, 'r')

	# for haplotype blocks
	blocks = {}

	for l in f:
		if re.search('#CHROM', l):
			samps = re.split('\t', l.rstrip())[9:]
			haplo = {}
			for samp in samps:
				haplo[samp] = {}

		elif not re.search('#', l):
			d = re.split('\t', l.rstrip())
			
			c = d[0]
			pos = int(d[1])
			ref = d[3]
			alt = d[4]
			keep = False

			# only consider biallelic SNPs
			# no indels or multiallelics
			if ref in ['A', 'T', 'C', 'G'] and alt in ['A', 'T', 'C', 'G']:
				genos = d[9:]
				for ix, fullgeno in enumerate(genos):
					geno = re.search('(\S\S\S)', fullgeno).group(1)
					geno = re.split('/', geno)
					# possible phase!
					if re.search('0/1', fullgeno):
						keep = True
						x = re.split(':', fullgeno)
					
						# this is where the phase data live
						phase = x[4]
					
						# tells the orient phase
						phase_pos = re.search('^(\d+)-', phase).group(1)
						phase_pos = int(phase_pos)		
				
						# define the new block
						if c not in blocks:
							blocks[c] = {}
						if phase_pos not in blocks[c]:
							blocks[c][phase_pos] = {}
						if pos not in blocks[c][phase_pos]:
							blocks[c][phase_pos][pos] = 1
		
						# define the SNP orientation within the block
						invert = False
						if re.search('-2.*-1', phase):
							invert = not invert
						if invert:
							geno = [geno[1], geno[0]]

					genos[ix] = geno
		
				if keep:
					if c not in haplo:
						haplo[c] = {}
					haplo[c][pos] = genos

	return blocks, haplo

def calc_D_sub(geno1, geno2):
	new1 = []
	new2 = []
	for g1, g2 in zip(geno1, geno2):
		if g1[0] != '.'  and g2[0] != '.':
			new1.append(g1)
			new2.append(g2)

	if len(new1) >= 4:
		pA = len([y for x in new1 for y in x if y == '0']) / float(len(new1) * 2)
		pB = len([y for x in new2 for y in x if y == '0']) / float(len(new2) * 2)

		d = 0
		for g1, g2 in zip(new1, new2):
			if g1[0] == '0' and g2[0] == '0':
				d += 1
			if g1[1] == '0' and g2[1] == '0':
				d += 1
		d = (d / float(len(new1) * 2)) - pA * pB

		denom = pA * (1 - pA) * pB * (1 - pB)
		if denom == 0:
			return np.nan
		else:
			return d ** 2 / float(denom)
	else:
		return np.nan
		

def calc_D(blocks, haplo):
	for locus in blocks:
		for pos1 in blocks[locus]:
			if len(blocks[locus][pos1]) > 1:
				for pos2 in blocks[locus][pos1]:
					if pos1 != pos2:
						D = calc_D_sub(haplo[locus][pos1], haplo[locus][pos2])
						if not math.isnan(D):
							print('%s,%s,%s,%s' % (locus, pos1, pos2, D))

def main():
	args = get_args()
	out_vcf = args.vcf
	outdir = args.outdir
	blocks, haplo = get_haplo(outdir, out_vcf)
	calc_D(blocks, haplo)

if __name__ == "__main__":
	main()
