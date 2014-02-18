import sys, os, eva
from argparse import ArgumentParser
import re, glob
from eva import *

def rev_comp(sequence):
    complement = {'A':'T','C':'G','G':'C','T':'A','N':'N'}
    return "".join([complement.get(nt.upper(), '') for nt in sequence[::-1]])

def simple_rev(args):
    if os.path.isfile(args.seq):
        ori = eva.FastaFile(args.seq)
        for tx_name, seq in ori.seqs.iteritems():
            ori.seqs[tx_name] = rev_comp(seq)
        ori.save_to_disk(args.seq + '.rev')
    else:
        print rev_comp(args.seq)

def add_file(fq_file, combined, counter):
    openf = open(fq_file, 'r')
    line_no = 0
    for line in openf:
        if line_no % 4 == 0:
            combined.write('@' + str(counter) + '\1\n')
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
    filebase, ext = os.path.splitext(args.input)
    renamed_file = '%s_did%s' % (filebase, ext)
    input_file = open(args.input, 'r')
    renamed = open(renamed_file, 'w')
    counter = 0
    line = input_file.readline()
    while line:
        if line.startswith('@'):
            renamed.write('@' + str(counter) + '\n')
            line = input_file.readline()
            renamed.write(line)
            line = input_file.readline()
            renamed.write(line)
            line = input_file.readline()
            renamed.write(line)
            counter += 1
        if line.startswith('>'):
            renamed.write('>%d\n' % counter)
            line = input_file.readline()
            renamed.write(line)
            counter += 1
        line = input_file.readline()
    input_file.close()
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

def rev(args):
    line_no = 0
    every_n_line = 2
    if args.type == 'fq':
        every_n_line = 4
    with open(args.out, 'w') as output:
        with open(args.input) as input:
            for line in input:
                if not line_no % every_n_line == 1:
                    output.write(line)
                else:
                    line = line.strip()
                    line = rev_comp(line)
                    output.write(line + '\n')
                line_no += 1
    print 'Check file %s' % args.out
    
def merge_peta_fa(args):
    id = 0
    every_n_line = 2
    first_c = '>'
    with open(args.out, 'w') as out:
        with open(args.left) as left:
            with open(args.right) as right:
                line = left.readline().strip()
                while line:
                    out.write('>%d\n' % id)
                    line = left.readline()
                    out.write(line)
                    id += 1
                    
                    line = right.readline()
                    line = right.readline()
                    out.write('>%d\n' % id)
                    out.write(line)
                    id += 1
                    
                    line = left.readline()
    print 'Check file %s' % args.out

def split2(args):
    left = args.fa[:-3] + '.left.fa'
    right = args.fa[:-3] + '.right.fa'
    with open(args.fa) as fa:
        with open(left, 'w') as l:
            with open(right, 'w') as r:
                line_no = 0
                left_id = 0
                right_id = 0
                for line in fa:
                    line = line.strip()
                    line_no += 1
                    if line_no % 4 == 1:
                        l.write('>%d/1\n' % left_id)
                    if line_no % 4 == 2:
                        l.write(rev_comp(line) + '\n')
                        left_id += 1
                    if line_no % 4 == 3:
                        r.write('>%d/2\n' % right_id)
                    if line_no % 4 == 0:
                        r.write(rev_comp(line) + '\n')
                        right_id += 1
    print 'Check left file %s' % left
    print 'Check right file %s' % right

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

def app(args):
    out_fa = FastaFile()
    count = 0
    for infile in glob.glob(os.path.join(args.dir, '*.fa')):
        comp_id = os.path.basename(infile)
        comp_id = comp_id[5:-3]
        fa = FastaFile(infile)
        for id, seq in fa.seqs.iteritems():
            out_fa.seqs[comp_id + '_comp_' + id] = seq
            count += 1
    out_fa.save_to_disk(args.out)
    print 'Check %d seqs in file %s' % (count, args.out)
    
def inspect(args):
    hits = read_blat_hits(args.psl, 'ref')
    with open(args.ids) as ids:
        for line in ids:
            id = line.strip()
            if not id in hits:
                print '======== %s not covered ========'
                continue
            tx_hits = hits[id]
            tx_hits.sort(key=lambda x: x.rstart)
            print '======== %s =======' % id
            for h in tx_hits:
                print h
                
def peta_shuffle(args):
    line_no = 0
    id = 0
    with open(args.fa) as fa:
        with open(args.fa[:-3] + '.shuffled.fa', 'w') as out:
            for line in fa:
                if '>' in line:
                    if line_no % 4 == 0:
                        out.write('>%d/1\n' % id)
                    else:
                        out.write('>%d/2\n' % id)
                        id += 1
                else:
                    line = line.strip()
                    if line_no % 4 == 1:
                        out.write(line + '\n')
                    else:
                        line = rev_comp(line)
                        out.write(line + '\n')
                line_no += 1

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
    parser_rename.add_argument('-t', required=True, help='fq or fa', dest='type', default='fq')
    parser_rename.add_argument('input', help='fastq file')
    
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

    parser_rev = subparsers.add_parser('reverse', help='reverse complement fasta/fastq')
    parser_rev.set_defaults(func=rev)
    parser_rev.add_argument('-t', required=True, help='fq or fa', dest='type', default='fq')
    parser_rev.add_argument('-f', required=True, help='file', dest='input', metavar='FILE')
    parser_rev.add_argument('-o', required=True, help='output file', dest='out', metavar='FILE')
    
    parser_merge_peta = subparsers.add_parser('merge_peta', help='simply merge left and right reads for PETA')
    parser_merge_peta.set_defaults(func=merge_peta_fa)
    parser_merge_peta.add_argument('-l', required=True, help='left mate file', dest='left', metavar='FILE')
    parser_merge_peta.add_argument('-r', required=True, help='right mate file', dest='right', metavar='FILE')
    parser_merge_peta.add_argument('-o', required=True, help='output file', dest='out', metavar='FILE')    
    
    parser_rev_comp = subparsers.add_parser('rev', help='simply get the reverse complement')
    parser_rev_comp.set_defaults(func=simple_rev)
    parser_rev_comp.add_argument('seq', help='ACGTN seq')
    
    parser_split2 = subparsers.add_parser('split', help='split the left and right mates')
    parser_split2.set_defaults(func=split2)
    parser_split2.add_argument('fa', help='FASTA file')
    
    parser_ps = subparsers.add_parser('peta_shuffle', help='shuffle peta FASTA for IDBA-tran')
    parser_ps.set_defaults(func=peta_shuffle)
    parser_ps.add_argument('fa', help='FASTA file')    
    
    parser_2files = subparsers.add_parser('2files', help='Split shuffled FASTQ file into two')
    parser_2files.set_defaults(func=to_files)
    parser_2files.add_argument('fq', help='fastq')
    
    parser_app = subparsers.add_parser('app', help='Append component id and merge FASTA files')
    parser_app.set_defaults(func=app)
    parser_app.add_argument('dir', help='directory')
    parser_app.add_argument('out', help='output file')
    
    parser_inspect = subparsers.add_parser('inspect', help='Append component id and merge FASTA files')
    parser_inspect.set_defaults(func=inspect)
    parser_inspect.add_argument('psl', help='psl file')
    parser_inspect.add_argument('ids', help='ids file')

    args = parser.parse_args()
    args.func(args)

if __name__ == '__main__':
    sys.exit(main())
