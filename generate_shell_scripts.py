import pandas as pd
import subprocess
import os

outdir = '/scratch/drabosky_flux/sosi/brazil_sub/'
file = '/scratch/drabosky_flux/sosi/brazil_sub/samples.csv'
d = pd.read_csv(file)

name = "align"
nodes = 1
cpu = 8
mem = 24
hours = 3

lineages = d['lineage'].unique().tolist()
# lineages = ['Bothrops_moojeni', 'Corallus_hortulanus']
samps = d.ix[d.lineage.isin(lineages), 'sample']

# for ix, lineage in enumerate(lineages):
for ix, sample in enumerate(samps):
	sh_out = '%s%s.sh' % (name, ix)
	o = open(sh_out, 'w')

	o.write("#PBS -N %s%s\n" % (name, ix))
	o.write("#PBS -M sosi@umich.edu\n")
	o.write("#PBS -A drabosky_flux\n")
	o.write("#PBS -l qos=flux\n")
	o.write("#PBS -q flux\n")
	o.write("#PBS -l nodes=%s:ppn=%s,mem=%sgb\n" % (nodes, cpu, mem))
	o.write("#PBS -l walltime=%s:00:00\n" % hours)
	o.write("#PBS -j oe\n")
	o.write("#PBS -V\n")

	o.write("\n")
	

	# o.write("python ~/squamateUCE/call_variants.py --lineage %s --file %s --dir %s --gatk ~/bin/GenomeAnalysisTK.jar --mem %s --CPU %s" % (lineage, file, outdir, mem, cpu))
	# o.write("python ~/squamateUCE/calculate_divergence_and_Fst.py --lineage %s --file %s --dir %s" % (lineage, file, outdir))
	# o.write("python ~/squamateUCE/calculate_pi_per_species.py --lineage %s --file %s --dir %s" % (lineage, file, outdir))
	# o.write("python ~/squamateUCE/align_reads2.py --lineage %s --file %s --dir %s --samtools ~/bin/samtools-1.3.1/samtools --gatk ~/bin/GenomeAnalysisTK.jar --dp 5 --qual 20 --CPU %s --mem %s" % (lineage, file, outdir, cpu, mem))
	o.write("python ~/squamateUCE/align_reads1.py --sample %s --file %s --dir %s --bwa ~/bin/bwa-0.7.12/bwa --samtools ~/bin/samtools-1.3.1/samtools --gatk ~/bin/GenomeAnalysisTK.jar --picard ~/bin/picard-tools-2.4.1/picard.jar --CPU %s --mem %s" % (sample, file, outdir, cpu, mem))
	# o.write("python ~/squamateUCE/make_PRG.py --lineage %s --file %s --dir %s --keep easy_recip_match,complicated_recip_match" % (lineage, file, outdir))
	# o.write("python ~/squamateUCE/match_contigs_to_probes.py --blat ~/bin/blat --sample %s --dir %s --evalue 1e-20 --db /scratch/drabosky_flux/sosi/brazil/squamate_AHE_UCE_genes_loci.fasta" % (sample, outdir))
	# o.write("python /home/sosi/squamateUCE/trinity_assembly.py --trinity ~/bin/trinityrnaseq-2.2.0/Trinity --sample %s --dir %s --mem %s --CPU %s --normal" % (sample, outdir, mem, cpu))
	# o.write("python ~/squamateUCE/clean_reads.py --trimjar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar --PEAR ~/bin/pear-0.9.10/pear-0.9.10-bin-64 --dir %s --sample %s --file %s\n" % (outdir, sample, file))
	o.close()

	out = '/scratch/drabosky_flux/sosi/brazil_sub/alignments/%s.realigned.dup.rg.mateFixed.sorted.recal.bam' % sample
	if os.path.isfile(out):
	  	os.remove(sh_out)
	# subprocess.call("qsub %s" % sh_out, shell=True)
