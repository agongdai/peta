import sys, os
from argparse import ArgumentParser

def rev_comp(sequence):
    complement = {'A':'T','C':'G','G':'C','T':'A','N':'N'}
    return "".join([complement.get(nt.upper(), '') for nt in sequence[::-1]])

def combine(args):
    files = args.inputs.split(',')
    name = args.output
    combined = open(name, 'w')
    counter = 0
    for f in files:
        openf = open(f, 'r')
        for line in openf:
            if '@' in line:
                combined.write('@' + str(counter) + '\n')
                counter += 1
            else:
                combined.write(line)
        openf.close()
    combined.close()

def merge_fq(args):
    left_file = args.left
    right_file = args.right
    try:
        left = open(left_file, 'r')
        right = open(right_file, 'r')
        filebase, ext = os.path.splitext(left_file)
        combined_file = filebase.split('_')[0] + '.fastq'
        combined = open(combined_file, 'w')
        counter = 0
        seqs = []
        line = left.readline()
        while line:
            if '@' in line:
#                print line
                combined.write('@' + str(counter) + '\n')
                l = left.readline()
                combined.write(l)
                l = left.readline()
                combined.write(l)
                l = left.readline()
                combined.write(l)
                counter += 1
                
                r = right.readline()
                combined.write('@' + str(counter) + '\n')
                r = right.readline()
                r = rev_comp(r) + '\n'
                combined.write(r)
                r = right.readline()
                combined.write(r)
                r = right.readline()
                combined.write(r)
                counter += 1
            line = left.readline()
                
        left.close()
        right.close()
        combined.close()
    except:
        print("IOError")

def main():
    parser = ArgumentParser()
    subparsers = parser.add_subparsers(help='sub command help')
    parser_combine = subparsers.add_parser('combine', help='combine multiple files')
    parser_combine.set_defaults(func=combine)
    parser_combine.add_argument('-i', required=True, help='the files to combine, seperated by a ","', dest='inputs')
    parser_combine.add_argument('-o', required=True, help='combined file name', dest='output')

    parser_merge = subparsers.add_parser('merge', help='merge the shuffled fasta/fastq files, reverse complement the right mates')
    parser_merge.set_defaults(func=merge_fq)
    parser_merge.add_argument('-l', required=True, help='left mate file', dest='left')
    parser_merge.add_argument('-r', required=True, help='right mate file', dest='right')
    args = parser.parse_args()
    args.func(args)

if __name__ == '__main__':
    sys.exit(main())
