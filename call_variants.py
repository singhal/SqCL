import argparse
import os
import pandas as pd
import re
import subprocess

"""
Sonal Singhal
created on 22 June 2016
Written assuming:
	* GATK 3.6
"""

def get_args():
	parser = argparse.ArgumentParser(
		description="Call SNPs. Assumes GATK 3.6",
        	formatter_class=argparse.ArgumentDefaultsHelpFormatter
		)

	# lineage
	parser.add_argument(
                '--lineage',
                type=str,
                default=None,
                help='Lineage for which to run script.'
                )

	# file
	parser.add_argument(
                '--file',
                type=str,
                default=None,
                help='File with sample info.'
                )

	# outdir
	parser.add_argument(
                '--outdir',
                type=str,
                default=None,
                help='Output directory for variants, '
                     ' only needed if running out of pipeline.'
                )

	# basedir
	parser.add_argument(
                '--dir',
                type=str,
                default=None,
                help="Full path to base dir with reads & assemblies & "
                     "everything else."
                )

	# bamfiles
	parser.add_argument(
		'--bamfile',
		type=str,
		default=None,
		help="Full path to file with BAM files, listed one "
		     "per line if running not in context of pipeline. "
		)

	# PRG
	parser.add_argument(
                '--prg',
                type=str,
                default=None,
                help="Full path to pseudoref genome if "
                     "you aren't running in context of pipeline."
                )

	# GATK
        parser.add_argument(
                '--gatk',
                type=str,
                default=None,
                help='GATK executable, full path.'
               )

	# memory
        parser.add_argument(
                '--mem',
                type=int,
                default=1,
                help='Memory available, as an int, in terms of Gb.'
               )

	# CPUs
	parser.add_argument(
                '--CPU',
                type=int,
                default=1,
                help='# of CPUs to use in alignment.'
               )

	# qual
	parser.add_argument(
		'--qual',
		type=int,
		default=20,
		help='Minimum quality to retain variant for '
                     'creating final call set.'
		)

	# depth
	parser.add_argument(
		'--dp',
		type=int,
		default=10,
		help='Minimum depth to retain variant for '
                     'creating final call set.'
		)

	return parser.parse_args()


def get_files(args):
	if args.bamfile:
		f = open(args.bamfile, 'r')
		files = []
		for l in f:
			files.append(l.rstrip())
		f.close()
	else:
		d = pd.read_csv(args.file)
		samps = d.ix[d['lineage'] == args.lineage, 'sample'].tolist()
		files = []
		for samp in samps:
			file = os.path.join(args.dir, 'alignments', 
                                            '%s.realigned.dup.rg.mateFixed.sorted.recal.bam' % samp)
			files.append(file)
	files = sorted(files)

	if args.prg:
		genome = args.prg
	else:
		genome = os.path.join(args.dir, 'PRG', '%s.fasta' % args.lineage)


	if args.outdir:
		outdir = args.outdir
	else:
		outdir = os.path.join(args.dir, 'variants')

	if not os.path.isdir(outdir):
		os.mkdir(outdir)
	
	return files, genome, outdir


def get_vcf(args, files, seq, dir):
	vcf_out1 = os.path.join(dir, '%s.raw.vcf' % args.lineage)
	vcf_out2 = os.path.join(dir, '%s.qual_filtered.vcf' % args.lineage)

	bam = '-I ' + ' -I '.join(files)

	subprocess.call("java -Xmx%sg -jar %s -T UnifiedGenotyper -R %s %s -o %s "
                        "--output_mode EMIT_ALL_SITES -nt %s"
                        % (args.mem, args.gatk, seq, bam, vcf_out1, args.CPU), shell=True) 

	subprocess.call("java -Xmx%sg -jar %s -T VariantFiltration -R %s -V %s "
                        "--filterExpression \"QUAL < %s\" --filterName \"LowQual\""
                        " -o %s -nt %s" % (args.mem, args.gatk, seq, vcf_out1,
                        args.qual, vcf_out2, args.CPU), shell=True)

	os.remove(vcf_out1)
	return vcf_out2


def depth_filter(args, infile, dir):
	out = os.path.join(dir, '%s.qual_filtered.cov_filtered.vcf' % args.lineage)

	f = open(infile, 'r')
	o = open(out, 'w')

	for l in f:
		if re.search('^#', l):
			o.write(l)
		else:
			d = re.split('\t', l.rstrip())
			if d[6] == 'PASS':
				tags = re.split(':', d[8])
				depth = tags.index('DP')
				genos = d[9:]
				miss = 0
				for ix, gen in enumerate(genos):
					info = re.split(':', gen)
					if int(info[depth]) < args.dp:
						d[ix + 9] = re.sub('^\S/\S', './.', d[ix + 9])
						miss += 1
				if miss < len(genos):
					o.write('\t'.join(d) + '\n')
	f.close()
	o.close()

	os.remove(infile)


def main():
	args = get_args()
	files, seq, dir = get_files(args)
	# get vcf and filter for qual
	out = get_vcf(args, files, seq, dir)
	# filter vcf for depth
	depth_filter(args, out, dir)

if __name__ == '__main__':
	main()
