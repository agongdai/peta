import sys, os
from argparse import ArgumentParser

def summarize_fa(args):
    fa = open(args.fa_file, 'r')
    n_base = 0
    ave_len = 0
    n_seqs = 0
    n50 = 0
    lengths = []
    seq_len = 0
    for line in fa:
        line = line.strip()
        if '>' in line:
            n_seqs += 1
            if seq_len > 0:
                lengths.append(seq_len)
                seq_len = 0
        else:
            n_base += len(line)
            seq_len += len(line)
    fa.close()
    lengths.sort()
    lengths.reverse()
    n50_len = 0
    for l in lengths:
        if n50_len * 2 >= n_base:
            n50 = l
            break
        else:
            n50_len += l
    if n_seqs > 0:
        ave_len = n_base / n_seqs
    print 'Summary of ' + args.fa_file
    print 'Total base: ' + str(n_base)
    print '# of sequences: ' + str(n_seqs)
    print 'Average length: ' + str(ave_len)
    print 'N50 value: ' + str(n50)

def main():
    parser = ArgumentParser()
    subparsers = parser.add_subparsers(help='sub command help')
    parser_differ = subparsers.add_parser('summary', help='get summary information of a fasta file')
    parser_differ.set_defaults(func=summarize_fa)
    parser_differ.add_argument('-f', required=True, help='the files to be summarized', dest='fa_file')

    args = parser.parse_args()
    args.func(args)

if __name__ == '__main__':
    sys.exit(main())
