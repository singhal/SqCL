import re

def get_seq(f):
	id = ''
	seq = {}
	f = open(f, 'r')
	for l in f:
		if re.search('>', l.rstrip()):
			id = re.search('>([^_]+)', l.rstrip()).group(1)
			seq[id] = ''
		else:
			seq[id] += l.strip()
	f.close()
	return seq

f = "/home/sosi/SqCL/squamate_AHE_UCE_genes_loci2.fasta"
seq = get_seq(f)
for id, s in seq.items():
	print('>%s\n%s' % (id, s))

