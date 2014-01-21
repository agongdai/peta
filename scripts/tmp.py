#
#line_no = 0
#with open('/home/carl/Projects/peta/rnaseq/hg19/SRX011545/SRR027876_corrected.fa') as ori:
#    with open('/home/carl/Projects/peta/rnaseq/hg19/SRX011545/SRR027876.trim1.fa', 'w') as out:
#        for line in ori:
#            line_no += 1
#            if line_no % 2 == 0:
#                line = line.strip()
#                line = line[1:-1]
#                out.write(line + '\n')
#            else:
#                out.write(line)

import sys, os
from eva import *
from argparse import ArgumentParser

def char2bit(args):
    bin = args.bin
    seq = ''
    if not len(bin) % 2 == 0:
        bin = '0' + bin
    for i in range(len(bin) / 2):
        c = bin[i:i+2]
        if c == '00':
            seq += 'A'
        if c == '01':
            seq += 'C'
        if c == '10':
            seq += 'G'
        if c == '11':
            seq += 'T'
    print seq
    
def to_files(args):
    left = args.fq[:-3] + '.1.fq'
    right = args.fq[:-3] + '.2.fq'
    with open(args.fq) as fq:
        line_no = 0
        with open(left, 'w') as left_mates:
            with open(right, 'w') as right_mates:
                for line in fq:
                    if line_no % 8 < 4:
                        left_mates.write(line)
                    else:
                        right_mates.write(line)
                    line_no += 1
    print 'Check %s and %s' % (left, right)

def count_small_seqs(args):
    fa = FastaFile(args.fa)
    out = FastaFile()
    count = 0
    for name, seq in fa.seqs.iteritems():
        if len(seq) >= 100:
            count += 1
            out.seqs[name] = seq
    out.save_to_disk(args.fa[:-3] + '.100.fa')
    print 'Check %s.100.fa' % (args.fa[:-3])
    print 'Longer than 100 in %s: %d / %d' % (args.fa, count, len(fa.seqs))
    
def get_mates(args):
    hits = read_blat_hits(args.psl, 'query')
#    fa = FastaFile(args.fa)
    out = FastaFile()
    for qname, hits in hits.iteritems():
        for h in hits:
            print h
            if h.n_match >= h.qlen - 1:
                print h
                
def get_bad_template_ids(args):
    hits = read_blat_hits(args.psl, 'ref')
    ids = []
    with open(args.ids) as f:
        for line in f:
            line = line.strip()
            if len(line) > 2:
                ids.append(line)
    for id in ids:
        if id in hits:
            tx_hits = hits[id]
            if len(tx_hits) <= 2:
                for h in tx_hits:
                    print '%s\t%s\t%.2f' % (h.qname.split('_')[0], id, float(h.n_match) / float(h.rlen))

def main():
    parser = ArgumentParser()
    subparsers = parser.add_subparsers(help='sub command help')
    
    parser_char2bit = subparsers.add_parser('b2n', help='Covert binary form to nucleotides')
    parser_char2bit.set_defaults(func=char2bit)
    parser_char2bit.add_argument('bin', help='binary')
    
    parser_2files = subparsers.add_parser('2files', help='Split shuffled FASTQ file into two')
    parser_2files.set_defaults(func=to_files)
    parser_2files.add_argument('fq', help='fastq')

    parser_count = subparsers.add_parser('count', help='counting')
    parser_count.add_argument('fa', help='fasta')
    parser_count.set_defaults(func=count_small_seqs)
    
    parser_tpl = subparsers.add_parser('tpl', help='get bad templates of some ids')
    parser_tpl.add_argument('ids', help='ids')
    parser_tpl.add_argument('psl', help='psl file')
    parser_tpl.set_defaults(func=get_bad_template_ids)
    
    parser_mate = subparsers.add_parser('mate', help='get mates in the psl file')
    parser_mate.add_argument('fa', help='fasta')
    parser_mate.add_argument('psl', help='psl')
    parser_mate.add_argument('out', help='out')
    parser_mate.set_defaults(func=get_mates)
    
    args = parser.parse_args()
    args.func(args)

if __name__ == '__main__':
    sys.exit(main())
    
