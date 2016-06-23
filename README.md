# squamateUCE - target capture pipeline
Scripts to work with UCE data from squamates.

## Areas of future improvements
- use read error correction prior to assembly
	- `http://www.genome.umd.edu/quorum.html`
- use cap3 to make PRG rather than random picking
- better trusted SNP set for BQSR
	- maybe do Cortex intersection with raw GATK call
- check SNP calling filters?
	- `https://www.broadinstitute.org/gatk/guide/article?id=3225`
	- `https://www.broadinstitute.org/gatk/guide/article?id=6925`
- phase SNPs

## Notes before you start
- When running any script, you can see the script arguments using `python script.py --help`.
- This pipeline is best run with these scripts all in a row, but each script can be run on its own. However, this piecemeal functionality hasn't been tested. 

## Steps
1. **Start with demultiplexed data.**
2. **Create a CSV file summarizing data**
	- note, column order does not matter.
	- sample: sample name
	- read1: complete path to read 1 file
	- read2: complete path to read 2 file
	- adaptor1: adaptor 1 to scrub; put asterisk where barcode goes
	- adaptor2: adaptor 2 to scrub; put asterisk where barcode goes
	- barcode1: barcode 1 to scrub
	- barcode2: barcode 2 to scrub (if it exists)
	- lineage: the group to which the sample belongs; used for assembly and SNP calling
3. **Clean reads**
	- Remove adaptors using `Trimmomatic`
	- Merge reads using `PEAR`
	- Lightly remove low quality sequence using `Trimmomatic`
	- Ends up creating three files: 2 paired files with ordering retained, and 1 unpaired read file
	- Assumes: `Trimmomatic v0.36` and `PEAR 0.9.10`
	```
	python ~/squamateUCE/clean_reads.py --trimjar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar \
		--PEAR ~/bin/pear-0.9.10/pear-0.9.10-bin-64 --dir /scratch/drabosky_flux/sosi/uce_test/ \
		--sample Anolis_carolinensis --file /scratch/drabosky_flux/sosi/uce_test/samples.csv
	```
4. **Assemble reads**
	- `Trinity` is meant for transcriptomes, but the authors of `Phyluce` have found it works well with non-transcriptome data sets
	- This script combines read 1 and unpaired reads into one file to use with `Trinity`
	- Then assembles
	- This is very memory intensive - recommended to give 1 Gb of RAM per 1M reads
	- Can use `--normal` flag to normalize big data sets
		- found this to be necessary because otherwise had heap error during `Butterfly` step
	- Requires `Trinity v2.2.0`
	```
	python /home/sosi/squamateUCE/trinity_assembly.py --trinity ~/bin/trinityrnaseq-2.2.0/Trinity \
		--sample Anolis_carolinensis --dir /scratch/drabosky_flux/sosi/uce_test/ \
		--mem 55 --CPU 16
	```
5. **Match assemblies to original targets**
	- Uses `blat` and a soft reciprocal match to identify which probes match to which contigs in the assembly
	- Multiple probes belong to the same locus, and while the probes don't necessarily overlap, they are likely to match to the same contig in the assembly
	- Script uses a regex that is dataset specific to synonymize probes to targeted loci
	- This script only considers matches that are above the e-value provided by the user and that are no more than 10 orders of magnitude worse than the best match
	- Have a few kinds of outcomes
		- *easy_recip_match*: 1-to-1 unique match between contig and targeted locus
		- *complicated_recip_match*: 1-to-1 non-unique match, in which one targeted locus matches to multiple contigs
		- *ditched_no_recip_match*: a case in which the contig matches to the targeted locus, but it isn't the best match
		- *ditched_no_match*: a case in which the contig matches to the targeted locus, but the locus doesn't match to the contig
		- *ditched_too_many_matches*: a case in which one contig has multiple good matches to multiple targeted loci
			- removed because possible paralogs?
	- This requires `blat v36`
	```
	python ~/squamateUCE/match_contigs_to_probes.py --blat ~/bin/blat --sample Anolis_carolinensis \
		--dir /scratch/drabosky_flux/sosi/uce_test/ --evalue 1e-30 \
		--db /scratch/drabosky_flux/sosi/uce_test/uce-5k-probes.fasta
	```
6. **Generate pseudo-reference genomes**
	- Uses lineage designation in sample file created in step 2 to create a pseudo-reference genome (PRG) for each lineage
	- The script identifies all contig matches to all targeted loci across all individuals in that lineage
		- if there are multiple contigs across individuals that match to a given targeted locus, then:
			- either picks the top contig if it is much better than the other matches (>1e3)
			- or if not, picks the longest one
	- This script can be run with only the 'easy' and / or 'complicated' reciprocal matches
	- Also will reverse complement contigs so that they are ordered the same as the contig (important later for phylogeny making!)
	- Requires no additional programs.
	```
	python ~/squamateUCE/make_PRG.py --lineage l1 --file /scratch/drabosky_flux/sosi/uce_test/samples.csv \
		--dir /scratch/drabosky_flux/sosi/uce_test/ --keep easy_recip_match
	```
7. **Align reads**
	- `align_reads1.py`: run by individual
		- Align paired reads & unpaired reads in two separate runs using `bwa`
			- paired read run performs mate rescue
		- Then run fixmate to fix any broken paired-end information
		- Then sort the BAM files
		- Then add the read groups as required by `GATK`
		- Then mark duplicates
		- Then identify indels
		- Then realign around indels
	- `align_reads2.py`: run by lineage first, then by individual
		- Generates raw set of variants using `GATK` for the lineage
		- Filters raw set lightly (`--dp` and `--qual`)
		- Uses filtered set to perform base quality score recalibration (BQSR)
		- Recalibrates individual BAM files
	- These scripts require `BWA 0.7.12`, `samtools 1.3.1`, `GATK 3.6`, `Picard 2.4.1`
	```
	python ~/squamateUCE/align_reads1.py --sample Mus_musculus \ 
		--file /scratch/drabosky_flux/sosi/uce_test/samples.csv \
		--dir /scratch/drabosky_flux/sosi/uce_test/ --bwa ~/bin/bwa-0.7.12/bwa \
		--samtools ~/bin/samtools-1.3.1/samtools --gatk ~/bin/GenomeAnalysisTK.jar \
		--picard ~/bin/picard-tools-2.4.1/picard.jar --CPU 1 --mem 1
		
	python ~/squamateUCE/align_reads2.py --lineage l1 --file /scratch/drabosky_flux/sosi/uce_test/samples.csv \
		--dir /scratch/drabosky_flux/sosi/uce_test/ --samtools ~/bin/samtools-1.3.1/samtools \
		--gatk ~/bin/GenomeAnalysisTK.jar --mem 3 --dp 10 --qual 20 --CPU 4
	```
8. **Call and filter variants**
	- Call variant (and invariant sites!) based on "final" BAM files per lineage
		- Note the importance of calling invariant sites - it is the denominator in all pop gen analyses
	- Variants are filtered based on quality (`--qual`)
	- And then filtered on depth, on a per individual basis (`--dp`)
		- If an individual has too low of coverage, 'ALLELE/ALLELE' becomes missing ('./.')
		- If every individual in the lineage is missing, then the site is dropped
	- This requires `GATK 3.6`
		- Note that we use `UnifiedGenotyper` instead of `HaplotypeCaller` because `HaplotypeCaller` output really odd results
	```
	python ~/squamateUCE/call_variants.py --lineage l1 --file /scratch/drabosky_flux/sosi/uce_test/samples.csv \
		--dir /scratch/drabosky_flux/sosi/uce_test/ --gatk ~/bin/GenomeAnalysisTK.jar --mem 4 \
		--CPU 4 --dp 10 --qual 20
	```
