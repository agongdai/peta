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
    
def app_pair_suffix(args):
    fq_file = args.fq
    ori = args.ori
    out = args.fq_with_tag
    try:
        fq = open(fq_file, 'r')
        fq_with_tag = open(out, 'w')
        if ori == 'right':
            tag = '/2'
        else:
            tag = '/1'
        counter = 0
        for line in fq:
            if counter % 4 == 0:
                line = line.strip()
                fq_with_tag.write(line + tag + '\n')
            else:
                fq_with_tag.write(line)
            counter += 1
        fq.close()
        fq_with_tag.close()
    except:
        print 'Unexpected Error:', sys.exc_info()[1]

def merge_fq(args):
    left_file = args.left
    right_file = args.right
    try:
        left = open(left_file, 'r')
        right = open(right_file, 'r')
        filebase, ext = os.path.splitext(left_file)
        combined_file = filebase.split('_')[0] + '_merged.fastq'
        combined = open(combined_file, 'w')
        counter = 0
        seqs = []
        line = left.readline()
        while line:
            if line.startswith('@'):
#               print line
                if args.pair:
                    combined.write('@' + str(counter) + '/1\n')
                else:
                    combined.write('@' + str(counter) + '\n')
                l = left.readline()
                combined.write(l)
                l = left.readline()
                combined.write(l)
                l = left.readline()
                combined.write(l)
                if not args.pair:
                    counter += 1
                
                r = right.readline()
                if args.pair:
                    combined.write('@' + str(counter) + '/2\n')
                else:
                    combined.write('@' + str(counter) + '\n')
                r = right.readline()
                if args.shuffle:
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
        if line.startswith('@'):
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

def fq2fa(args):
    fq = open(args.fq_fn)
    fa = open(args.fa_fn, 'w')
    line_no = 0
    for line in fq:
        if line_no % 4 == 0:
            line = line.replace('@', '>')
            fa.write(line)
        elif line_no % 4 == 1:
            fa.write(line)
        line_no += 1
    fq.close()
    fa.close()
    
def decode(args):
    fq = open(args.fq)
    filebase, ext = os.path.splitext(args.fq)
    left = open(filebase + '_1.fastq', 'w')
    right = open(filebase + '_2.fastq', 'w')
    counter = 0
    line_no = 0
    for line in fq:
        line = line.strip()
        if line_no % 4 == 0:
            left.write('@' + str(counter) + '\n')
            right.write('@' + str(counter) + '\n')
            counter += 1
        elif line_no % 4 == 1 or line_no % 4 == 3:
            left.write(line[0:68] + '\n')
            right.write(line[76::] + '\n')
        else:
            left.write(line + '\n')
            right.write(line + '\n')
        line_no += 1
    fq.close()
    left.close()
    right.close()

def main():
    parser = ArgumentParser()
    subparsers = parser.add_subparsers(help='sub command help')
    parser_combine = subparsers.add_parser('combine', help='combine multiple files')
    parser_combine.set_defaults(func=combine)
    parser_combine.add_argument('-i', required=True, help='the files to combine, seperated by a ","', dest='inputs')
    parser_combine.add_argument('-o', required=True, help='combined file name', dest='output')

    parser_merge = subparsers.add_parser('merge', help='merge the shuffled fasta/fastq files')
    parser_merge.set_defaults(func=merge_fq)
    parser_merge.add_argument('-l', required=True, help='left mate file', dest='left', metavar='FILE')
    parser_merge.add_argument('-r', required=True, help='right mate file', dest='right', metavar='FILE')
    parser_merge.add_argument('-s', required=False, help='shuffle or not', dest='shuffle', default=False, metavar='FILE')
    parser_merge.add_argument('-p', required=False, help='pair or not', dest='pair', default=True, metavar='FILE')
    
    parser_merge = subparsers.add_parser('shuffle', help='merge the shuffled fasta/fastq files, reverse complement the right mates')
    parser_merge.set_defaults(func=merge_fq)
    parser_merge.add_argument('-l', required=True, help='left mate file', dest='left', metavar='FILE')
    parser_merge.add_argument('-r', required=True, help='right mate file', dest='right', metavar='FILE')
    parser_merge.add_argument('-s', required=False, help='shuffle or not', dest='shuffle', default=True, metavar='FILE')

    parser_rename = subparsers.add_parser('rename', help='rename read ids to be 0, 1, 2...')
    parser_rename.set_defaults(func=rename_ids)
    parser_rename.add_argument('-f', required=True, help='fastq file', dest='fq', metavar='FILE')
    
    parser_rename = subparsers.add_parser('decode', help='remove the barcode')
    parser_rename.set_defaults(func=decode)
    parser_rename.add_argument('-f', required=True, help='fastq file', dest='fq', metavar='FILE')

    parser_fq2fa = subparsers.add_parser('fq2fa', help='convert fastq file to fasta file')
    parser_fq2fa.set_defaults(func=fq2fa)
    parser_fq2fa.add_argument('-q', required=True, help='fastq file', dest='fq_fn', metavar='FILE')
    parser_fq2fa.add_argument('-a', required=True, help='fasta file', dest='fa_fn', metavar='FILE')

    parser_app_pair_suffix = subparsers.add_parser('pair', help='add /1 or /2 to the fastq record ids')
    parser_app_pair_suffix.set_defaults(func=app_pair_suffix)
    parser_app_pair_suffix.add_argument('-q', required=True, help='fastq file', dest='fq', metavar='FILE')
    parser_app_pair_suffix.add_argument('-o', required=True, help='output fastq file', dest='fq_with_tag', metavar='FILE')
    parser_app_pair_suffix.add_argument('-d', required=True, default='left', help='left or right', dest='ori')

    args = parser.parse_args()
    args.func(args)

if __name__ == '__main__':
    sys.exit(main())
