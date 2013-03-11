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
    seq_file = args.seq_file
    ori = args.ori
    out = args.fq_with_tag
    delimiter = 4
    if args.type == 'fa':
        delimiter = 2
    try:
        fq = open(seq_file, 'r')
        fq_with_tag = open(out, 'w')
        if ori == 'right':
            tag = '/2'
        else:
            tag = '/1'
        counter = 0
        for line in fq:
            if counter % delimiter == 0:
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
    
def trim(args):
    filebase, ext = os.path.splitext(args.fa)
    trimmed_fn = '%s.trim%d.fa' % (filebase, args.n_trim)
    trimmed = open(trimmed_fn, 'w')
    fa = open(args.fa)
    line_no = 0
    for line in fa:
        if line_no % 2 == 0:
            trimmed.write(line)
        else:
            line = line.strip()
            line = line[args.n_trim:0 - args.n_trim]
            trimmed.write(line + '\n')
        line_no += 1
    trimmed.close()
    fa.close()
    print 'Check trimmed file %s ' % trimmed_fn
    
def is_n_rep(line, n):
    if len(line) - n <= 0:
        return 0
    is_n_rep = 1
    for i in range(len(line) - n):
        if not line[i] == line[i + n]:
            is_n_rep = 0
            break
    if is_n_rep:
        return 1
    return 0
    
def view(args):
    f = open(args.fa)
    filebase, ext = os.path.splitext(args.fa)
    filtered_fn = '%s.filtered.fa' % (filebase)
    filtered = open(filtered_fn, 'w')
    n_total = 0
    n_has_n = 0
    n_all_same = 0
    n_biased = 0
    n_2_rep = 0
    n_3_rep = 0
    n_part_rep = 0
    read_id = 0
    for line in f:
        if '>' in line:
            continue
        n_total += 1
        counter = [0 for x in range(5)]
        if n_total % 1000000 == 0:
            print 'Processed %d...' % n_total
        line = line.upper().strip()
        if 'N' in line:
            n_has_n += 1
            continue
        for i in range(len(line)):
            if line[i] == 'A':
                counter[0] += 1
            elif line[i] == 'C':
                counter[1] += 1
            elif line[i] == 'G':
                counter[2] += 1
            elif line[i] == 'T':
                counter[3] += 1
            elif line[i] == 'N':
                counter[4] += 1
        if len(line) in counter:
            n_all_same += 1
            #print 'All same: ' + line
            continue
        if len(line) - 1 in counter or len(line) - 2 in counter or len(line) - 3 in counter or len(line) - 4 in counter or len(line) - 5 in counter:
            n_biased += 1
            #print 'Biased: ' + line
            continue
        if is_n_rep(line, 2):
            n_2_rep += 1
            #print '2 bases: ' + line
            continue
        if is_n_rep(line, 3):
            n_3_rep += 1
            continue
        filtered.write('>%d\n' % read_id)
        filtered.write(line + '\n')
        read_id += 1
        for i in range(len(line) - 8):
            if is_n_rep(line[i:i+10], 2):
                n_part_rep += 1
                #print 'Partial 2 bases: ' + line
                break
    print 'Total reads:              %d' % n_total
    print 'Has N\'s:                  %d' % n_has_n
    print 'All base same:            %d' % n_all_same
    print 'Bases biased:             %d' % n_biased
    print '2 bases repeat:           %d' % n_2_rep
    print '3 bases repeat:           %d' % n_3_rep
    print 'Partial repeat:           %d' % n_part_rep
    print '%d reads remained. Check file %s ' % (read_id, filtered_fn)

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
    parser_merge.add_argument('-p', required=False, help='pair or not', dest='pair', default=False, metavar='FILE')

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

    parser_app_pair_suffix = subparsers.add_parser('pair', help='add /1 or /2 to the fasta/fastq record ids')
    parser_app_pair_suffix.set_defaults(func=app_pair_suffix)
    parser_app_pair_suffix.add_argument('-t', required=True, help='fq or fa', dest='type', default='fq')
    parser_app_pair_suffix.add_argument('-f', required=True, help='file', dest='seq_file', metavar='FILE')
    parser_app_pair_suffix.add_argument('-o', required=True, help='output file', dest='fq_with_tag', metavar='FILE')
    parser_app_pair_suffix.add_argument('-d', required=True, default='left', help='left or right', dest='ori')
    
    parser_trim = subparsers.add_parser('trim', help='trim the reads, head and tail')
    parser_trim.set_defaults(func=trim)
    parser_trim.add_argument('-f', '--fasta', required=True, help='fasta file', dest='fa', metavar='FILE')
    parser_trim.add_argument('-n', required=False, type=int, default=2, help='fasta file', dest='n_trim')
    
    parser_view = subparsers.add_parser('view', help='view the quality of a Fasta file')
    parser_view.set_defaults(func=view)
    parser_view.add_argument('fa', help='fasta file', metavar='FILE')

    args = parser.parse_args()
    args.func(args)

if __name__ == '__main__':
    sys.exit(main())
