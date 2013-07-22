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
    count = 0
    for name, seq in fa.seqs.iteritems():
        if len(seq) >= 100:
            count += 1
    print 'Longer than 100 in %s: %d' % (args.fa, count)

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
    
    args = parser.parse_args()
    args.func(args)

if __name__ == '__main__':
    sys.exit(main())
    
