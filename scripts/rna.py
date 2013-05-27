from eva import *
from zoom import *
from argparse import ArgumentParser
import operator

class Junction(object):
    def __init__(self, main_id, branch_id, locus = 0, weight = 0, ori = 0):
        self.main = main_id
        self.branch = branch_id
        self.locus = locus
        self.weight = weight
        self.ori = ori

# Split string like '[1, 233]'        
def read_tpl(tpl_str):
    tpl_str = tpl_str[1:-1]
    f = tpl_str.split(',')
    id = int(f[0].strip())
    length = int(f[1].strip())
    return id, length

# Read junctions in format of: [main_id, main_length]    [branch_id, branch_length]    locus    weight    orientation
def read_junctions(junction_file):
    junctions = []
    with open(junction_file) as jf:
        line_no = 0
        for line in jf:
            line_no += 1
            if line_no <= 1:
                continue
            line = line.strip()
            f = line.split('\t')
            new_j = Junction(read_tpl(f[0]), read_tpl(f[1]), int(f[2]), int(f[3]), int(f[4]))
            junctions.append(new_j)
    return junctions

def exam_junctions(args):
    junctions = read_junctions(args.junction)
    hits = read_blat_hits(args.psl, 'ref')
    n_one_contig_full_length = 0
    for tx_name, tx_hits in hits.iteritems():
        tx_hits.sort(key=lambda x: x.rstart)
        stop_here = False
        tx_len = 0
        for h in tx_hits:
            tx_len = h.rlen
            if h.n_match >= h.rlen - 10:
#                 print 'Full length: %s' % tx_name
                n_one_contig_full_length += 1
                stop_here = True
                break
        if not stop_here:
            bases_covered = [0 for _ in range(tx_len)]
            for h in tx_hits:
                for i in range(h.n_blocks):
                    for j in range(h.r_block_starts[i], h.r_block_starts[i] + h.block_sizes[i]):
                        bases_covered[j] = 1
            n_covered = 0
            for c in bases_covered:
                n_covered += c
            if n_covered < tx_len - 10:
                stop_here = True
                print '%s is not fully covered.' % tx_name
        if not stop_here:
            
        
    print 'n_one_contig_full_length: %d' % n_one_contig_full_length

def split_exons(args):
    base_coverage = {}
    base_coverage['chr1'] = [0 for _ in range(5579133)]
    base_coverage['chr2'] = [0 for _ in range(4539804)]
    base_coverage['chr3'] = [0 for _ in range(2452883)]
    line_no = 0
    with open(args.psl) as psl:
        for line in psl:
            line_no += 1
            if line_no <= 5:
                continue
            line = line.strip()
            h = read_psl_line(line)
            for i in range(h.n_blocks):
                for j in range(h.block_sizes[i]):
                    base_coverage[h.rname][h.r_block_starts[i] + j] += 1
            if line_no % 1000000 == 0:
                print >> sys.stderr, '# of reads counted: %d ' % line_no

    exon_start = 0
    exon_end = 0
    weight = 0.0
    for chr, bases in base_coverage.iteritems():
        for i in range(len(bases)):
            cov = bases[i]
            if cov > 0:
                exon_end = i + 1
                weight += cov
            else:
                if exon_end > exon_start:
                    weight /= exon_end - exon_start
                    print '%s\t%d\t%d\t%.2f' % (chr, exon_start, exon_end, weight)
                exon_start = i + 1
                weight = 0.0
    if exon_end > exon_start:
        weight /= exon_end - exon_start
        print '%s\t%d\t%d\t%.2f' % (chr, exon_start, exon_end, weight)

def main():
    parser = ArgumentParser()
    subparsers = parser.add_subparsers(help='sub command help')
    
    parser_exon = subparsers.add_parser('splice', help='Split exons from PSL file')
    parser_exon.set_defaults(func=split_exons)
    parser_exon.add_argument('psl', help='transcript/contigs-to-genome PSL file')
    
    parser_exon = subparsers.add_parser('junction', help='Exam junctions from PSL file')
    parser_exon.set_defaults(func=exam_junctions)
    parser_exon.add_argument('junction', help='PETA junction file')
    parser_exon.add_argument('psl', help='transcript/contigs-to-genome PSL file')

    args = parser.parse_args()
    args.func(args)
                
if __name__ == '__main__':
    sys.exit(main())