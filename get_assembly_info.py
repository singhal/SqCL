import argparse
import re
import pandas as pd
import os
import numpy as np

"""
Sonal Singhal
created on 2 August 2016
added to pipeline Oct 2, 2016
Written assuming: NOTHING!
"""

def get_args():
    parser = argparse.ArgumentParser(
        description="Summarizes assemblies & PRGs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )

    # file
    parser.add_argument(
        '--file',
        type=str,
        default=None,
        help='File with sample info.'
        )

    # basedir
    parser.add_argument(
        '--dir',
        type=str,
        default=None,
        help="Full path to base dir with reads & assemblies & "
             "everything else."
        )

    # assembly dir
    parser.add_argument(
            '--gdir',
            type=str,
            default=None,
            help="Full path to trinity assembly dir if "
                 "you aren't running in context of pipeline."
            )

    # PRG dir
    parser.add_argument(
        '--pdir',
            type=str,
            default=None,
            help="Full path to dir with PRG if "
                 "you aren't running in context of pipeline."
            )

    # out dir
    parser.add_argument(
        '--outdir',
            type=str,
            default=None,
            help="Full path to dir for data output "
                 "if you aren't running in context of pipeline."
            )

    return parser.parse_args()


def get_data(args):
    d = pd.read_csv(args.file)

    if args.gdir:
        gdir = args.gdir
    else:
        gdir = os.path.join(args.dir, 'trinity_assembly')

    if args.pdir:
        pdir = args.pdir
    else:
        pdir = os.path.join(args.dir, 'PRG')

    if args.outdir:
        outdir = args.outdir
    else:
        outdir = os.path.join(args.dir, 'metadata')
    
    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    return d, gdir, pdir, outdir

def summarize_contigs(d, gdir, pdir, outdir):
    out = os.path.join(outdir, 'assembly_info.csv')
    o = open(out, 'w')
    o.write('sample,assembled_contigs,assembled_length,assembled_n50,on_target\n')

    for ind, sp in zip(d['sample'], d['lineage']):
        a = os.path.join(gdir, '%s.fasta' % ind)
        
        seq = {}
        f = open(a, 'r')
        for l in f:
            if re.search('>', l):
                id = re.search('>(\S+)', l).group(1)
                seq[id] = ''
            else:
                seq[id] += l.rstrip()
        f.close()

        seq = [len(y) for x,y in seq.items()]
        seq = sorted(seq)
        contig_count = len(seq)
        tot_length = np.sum(seq)
        
        run_count = 0
        for ix, val in enumerate(seq):
            run_count += val
            if run_count >= (tot_length / 2.0):
                n50 = val
                break

        prg = os.path.join(pdir, '%s.fasta' % sp)
        target = 0
        f = open(prg, 'r')
        for l in f:
            if re.search('>', l):
                target += 1
        f.close()
        target = round(target / float(contig_count), 3)

        o.write('%s,%s,%s,%s,%s\n' % (ind, contig_count, tot_length, n50, target))
    o.close()


def main():
    args = get_args()
    d, gdir, pdir, outdir = get_data(args)
    summarize_contigs(d, gdir, pdir, outdir)

if __name__ == "__main__":
    main()


            
