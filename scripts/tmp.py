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

def main():
    parser = ArgumentParser()
    subparsers = parser.add_subparsers(help='sub command help')
    
    parser_char2bit = subparsers.add_parser('b2n', help='Covert binary form to nucleotides')
    parser_char2bit.set_defaults(func=char2bit)
    parser_char2bit.add_argument('bin', help='binary')
    
    args = parser.parse_args()
    args.func(args)

if __name__ == '__main__':
    sys.exit(main())
    
