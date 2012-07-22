import sys, os
from argparse import ArgumentParser

def rev_comp(sequence):
    complement = {'A':'T','C':'G','G':'C','T':'A','N':'N'}
    return "".join([complement.get(nt.upper(), '') for nt in sequence[::-1]])

def add_file(fq_file, combined, counter):
    openf = open(fq_file, 'r')
    line_no = 0
    for line in openf:
        if line_no % 4 == 0:
            combined.write('@' + str(counter) + '\n')
            counter += 1
        else:
            combined.write(line)
        line_no += 1
    openf.close()
    return counter

def combine(args):
    files = args.inputs.split(',')
    name = args.output
    combined = open(name, 'w')
    counter = 0
    for f in files:
        counter = add_file(f, combined, counter)
        print 'File ' + f + ' combined. # of seqs: ' + str(counter)
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
            if line.startswith('@'):
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
        print 'Unexpected Error:', sys.exc_info()[1]

def rename_ids(args):
    filebase, ext = os.path.splitext(args.fq)
    renamed_file = filebase + '_did.fastq'
    fq_file = open(args.fq, 'r')
    renamed = open(renamed_file, 'w')
    counter = 0
    line = fq_file.readline()
    while line:
        if line.startwith('@'):
            renamed.write('@' + str(counter) + '\n')
            line = fq_file.readline()
            renamed.write(line)
            line = fq_file.readline()
            renamed.write(line)
            line = fq_file.readline()
            renamed.write(line)
            counter += 1
        line = fq_file.readline()
    fq_file.close()
    renamed.close()

def main():
    parser = ArgumentParser()
    subparsers = parser.add_subparsers(help='sub command help')
    parser_combine = subparsers.add_parser('combine', help='combine multiple files')
    parser_combine.set_defaults(func=combine)
    parser_combine.add_argument('-i', required=True, help='the files to combine, seperated by a ","', dest='inputs')
    parser_combine.add_argument('-o', required=True, help='combined file name', dest='output')

    parser_merge = subparsers.add_parser('merge', help='merge the shuffled fasta/fastq files, reverse complement the right mates')
    parser_merge.set_defaults(func=merge_fq)
    parser_merge.add_argument('-l', required=True, help='left mate file', dest='left', metavar='FILE')
    parser_merge.add_argument('-r', required=True, help='right mate file', dest='right', metavar='FILE')

    parser_merge = subparsers.add_parser('rename', help='rename read ids to be 0, 1, 2...')
    parser_merge.set_defaults(func=rename_ids)
    parser_merge.add_argument('-f', required=True, help='fastq file', dest='fq', metavar='FILE')

    args = parser.parse_args()
    args.func(args)

if __name__ == '__main__':
    sys.exit(main())
