from eva import *
from zoom import *
from merge import *
from argparse import ArgumentParser
import operator, random

class Junction(object):
    def __init__(self, main_id, branch_id, main_len = 0, branch_len = 0, locus = 0, weight = 0, ori = 0):
        self.main = main_id
        self.main_len = main_len
        self.branch = branch_id
        self.branch_len = branch_len
        self.locus = locus
        self.weight = weight
        self.ori = ori
        
    def __repr__(self):
        s = '[%s,%d]\t[%s,%d]\t%d\t%d\t%d\n' % (self.main, self.main_len, self.branch, self.branch_len, self.locus, self.weight, self.ori)
        return s

# Split string like '[1, 233]'        
def read_tpl_info(tpl_str):
    tpl_str = tpl_str[1:-1]
    f = tpl_str.split(',')
    tpl_id = f[0].strip()
    length = int(f[1].strip())
    return tpl_id, length

# Read junctions in format of: [main_id, main_length]    [branch_id, branch_length]    locus    weight    orientation
def read_junctions(junction_file):
    junctions = []
    with open(junction_file) as jf:
        line_no = 0
        for line in jf:
            line_no += 1
            line = line.strip()
            if line_no <= 1 or line == '':
                continue
            f = line.split('\t')
            if f[-1] == '-1':
                continue
            main_id, main_len = read_tpl_info(f[0])
            branch_id, branch_len = read_tpl_info(f[1])
            new_j = Junction(main_id, branch_id, main_len, branch_len, int(f[2]), int(f[3]), int(f[4]))
            junctions.append(new_j)
    return junctions

def get_junction_dict(junctions):
    my_junctions = {}
    for j in junctions:
        if not j.main in my_junctions:
            my_junctions[j.main] = []
        if not j.branch in my_junctions:
            my_junctions[j.branch] = []
        my_junctions[j.main].append(j.branch)
        my_junctions[j.branch].append(j.main)
    return my_junctions

visited = []

def exons_are_connected(exon_1, exon_2, j_dict, level = 0):
    if level > 5:
        return False
    if not exon_1 in j_dict or not exon_2 in j_dict:
        return False
#     if len(visited) == 0:
#         visited.append(exon_1)
    connected_1 = j_dict[exon_1]
    if exon_2 in connected_1:
        return True
#     for e in connected_1:
#         if e in visited:
#             continue
#         visited.append(e)
#         # print exon_1, e, exon_2
#         if exons_are_connected(e, exon_2, j_dict, level + 1):
#             return True
    return True
    
def exam_junctions(args):
    junctions = read_junctions(args.junction)
    j_dict = get_junction_dict(junctions)
    existing = []
    with open(args.full) as full:
        for line in full:
            line = line.strip()
            if not line == '':
                existing.append(line)
    print 'Existing full length: %d' % (len(existing))
    
    hits = read_blat_hits(args.psl, 'ref')
    n_one_contig_full_length = 0
    n_covered_by_connected_graph = 0
    full = open(args.junction + '.full', 'w')
    for tx_name, tx_hits in hits.iteritems():
        tx_hits.sort(key=lambda x: x.rstart)
        stop_here = False
        tx_len = 0
        
        if tx_name in existing:
            continue
        
        for h in tx_hits:
            tx_len = h.rlen
        
        # If the transcript is not expressed, ignore
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
            # print '%s is not fully covered.' % tx_name
        if stop_here:
            continue
        
        max_cover_hits = []
        back_bone_h = None
        for h in tx_hits:
            if not back_bone_h:
                back_bone_h = h
                max_cover_hits.append(h)
            else:
                if h.rend < back_bone_h.rend:
                    continue
                if h.rstart < back_bone_h.rend + 10:
                    visited = []
                    if exons_are_connected(back_bone_h.qname, h.qname, j_dict):
                        back_bone_h = h
                        max_cover_hits.append(h)
                        
        bases_covered = [0 for _ in range(tx_len)]
        for h in max_cover_hits:
            for i in range(h.n_blocks):
                for j in range(h.r_block_starts[i], h.r_block_starts[i] + h.block_sizes[i]):
                    bases_covered[j] = 1
        n_covered = 0
        for c in bases_covered:
            n_covered += c
        if n_covered >= tx_len - 10 and len(max_cover_hits) > 1:
            print 'Transcript %s is recovered by connected graph' % tx_name 
            print max_cover_hits
            print '--------------------------------------------------------'
            n_covered_by_connected_graph += 1
            full.write(tx_name + '\n')
    
    full.close()    
    print 'Check file %s.full' % args.junction
    print 'n_one_contig_full_length: %d' % len(existing)
    print 'n_covered_by_connected_graph: %d' % n_covered_by_connected_graph

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

def reads_from_contigs(args):
    ref = FastaFile(args.contigs)
    out = '%s.%d.fa' % (args.contigs[:-3], args.read_len)
    with open(out, 'w') as reads:
        id = 0
        for chr, seq in ref.seqs.iteritems():
            seq = seq.upper()
            #id = 0
            print chr
            for i in range(len(seq) - (args.read_len - 1)):
                #head = '>%s:%d\n' % (chr, id)
                head = '>%d %s\n' % (id, chr)
                reads.write(head)
                reads.write(seq[i:i+args.read_len] + '\n')
                id += 1
    print 'Check file %s' % out

def read_grep_lines(lines, sub):
    ls = lines.split('\n')
    if len(ls) < 3:
        return ''
    line_no = 0
    id = ''
    seq = ''
    l = ls[1]
    readable = ' ' * (len(l) - len(sub)) + sub + '\n'
    for l in ls:
        line_no += 1
        l = l.strip()
        if line_no % 3 == 1:
            id = l[1:]
        if line_no % 3 == 2:
            seq = l
            index = seq.find(sub)
            seq = ' ' * (len(seq) - len(sub) - index) + seq + ' ' * index
            readable += '%s\t%s\n' % (seq, id)
    return readable

def grep_kmer(args):
    fa25 = args.ref
    print '--------- %s ---------' % args.kmer
    cmd = 'grep -B 1 %s %s' % (args.kmer, fa25)
    print read_grep_lines(runInShell(cmd), args.kmer)
    if args.ori == 'right':
        print '>>>>>>>>> To the right >>>>>>>>>'
        first_24 = args.kmer[1:25]
        
        next = first_24 + 'A'
        print '--------- KMER: %s ---------' % args.kmer
        print '--------- NEXT:  %s ---------' % next
        cmd = 'grep -B 1 %s %s' % (next, fa25)
        print read_grep_lines(runInShell(cmd), next)
        
        next = first_24 + 'C'
        print '--------- KMER: %s ---------' % args.kmer
        print '--------- NEXT:  %s ---------' % next
        cmd = 'grep -B 1 %s %s' % (next, fa25)
        print read_grep_lines(runInShell(cmd), next)
        
        next = first_24 + 'G'
        print '--------- KMER: %s ---------' % args.kmer
        print '--------- NEXT:  %s ---------' % next
        cmd = 'grep -B 1 %s %s' % (next, fa25)
        print read_grep_lines(runInShell(cmd), next)
        
        next = first_24 + 'T'
        print '--------- KMER: %s ---------' % args.kmer
        print '--------- NEXT:  %s ---------' % next
        cmd = 'grep -B 1 %s %s' % (next, fa25)
        print read_grep_lines(runInShell(cmd), next)
        
    else:
        print '<<<<<<<<< To the left <<<<<<<<<'
        last_24 = args.kmer[0:24]        
        
        next = 'A' + last_24
        print '--------- KMER:  %s ---------' % args.kmer
        print '--------- NEXT: %s ---------' % next
        cmd = 'grep -B 1 %s %s' % (next, fa25)
        print read_grep_lines(runInShell(cmd), next)
        
        next = 'C' + last_24
        print '--------- KMER:  %s ---------' % args.kmer
        print '--------- NEXT: %s ---------' % next
        cmd = 'grep -B 1 %s %s' % (next, fa25)
        print read_grep_lines(runInShell(cmd), next)
        
        next = 'G' + last_24
        print '--------- KMER:  %s ---------' % args.kmer
        print '--------- NEXT: %s ---------' % next
        cmd = 'grep -B 1 %s %s' % (next, fa25)
        print read_grep_lines(runInShell(cmd), next)
        
        next = 'T' + last_24
        print '--------- KMER:  %s ---------' % args.kmer
        print '--------- NEXT: %s ---------' % next
        cmd = 'grep -B 1 %s %s' % (next, fa25)
        print read_grep_lines(runInShell(cmd), next)    
    
def find_kmer(args):
    grep_kmer(args)
    print '==================== Reverse complement ====================='
    args.kmer = rev_comp(args.kmer)
    if args.ori == 'right':
        args.ori = 'left'
    else:
        args.ori = 'right'
    grep_kmer(args)

def fa2single(args):
    out = args.fa[:-3] + '.single.fa'
    fa = FastaFile(args.fa)
    with open(out, 'w') as single:
        for name, seq in fa.seqs.iteritems():
            single.write('>%s\n' % name)
            single.write('%s\n' % seq)
    print 'Check file %s' % out
    
def simu(args):
    annotated = FastaFile(args.annotated)
    reads = FastaFile()
    read_id = 0
    out = '/home/carl/Projects/peta/rnaseq/Spombe/SRR097897/simu.fa'
    ref = '/home/carl/Projects/peta/rnaseq/Spombe/genome/spombe.broad.tx.fasta.rev'
    for tx in args.transcripts:
        seq = annotated.seqs[tx]
        for i in range(len(seq) - args.read_len + 1):
            reads.seqs[read_id] = seq[i:i+args.read_len]
            read_id += 1
            reads.seqs[read_id] = seq[i:i+args.read_len]
            read_id += 1
    reads.save_to_disk(out)
    cmd = 'blat %s %s %s.psl' % (ref, out, out)
    print cmd
    print runInShell(cmd)
    print 'Check file %s' % out

def main():
    parser = ArgumentParser()
    subparsers = parser.add_subparsers(help='sub command help')
    
    parser_splice = subparsers.add_parser('splice', help='Split exons from PSL file')
    parser_splice.set_defaults(func=split_exons)
    parser_splice.add_argument('psl', help='transcript/contigs-to-genome PSL file')
    
    parser_g2r = subparsers.add_parser('c2r', help='Split the contigs to kmers')
    parser_g2r.set_defaults(func=reads_from_contigs)
    parser_g2r.add_argument('contigs', help='contigs Fasta file')
    parser_g2r.add_argument('read_len', type=int, help='read length')
    
    parser_kmer = subparsers.add_parser('kmer', help='Find kmers in RNA-seq library')
    parser_kmer.set_defaults(func=find_kmer)
    parser_kmer.add_argument('-r', '--ref', required=False, default='/home/carl/Projects/peta/rnaseq/Spombe/SRR097897/SRR097897.fa', help='single line Fasta file to grep kmer')
    parser_kmer.add_argument('kmer', help='kmer sequence')
    parser_kmer.add_argument('ori', help='left or right')
    
    parser_single = subparsers.add_parser('f2s', help='Convert fasta file to single line sequences')
    parser_single.set_defaults(func=fa2single)
    parser_single.add_argument('fa', help='Fasta file')
    
    parser_exon = subparsers.add_parser('junction', help='Exam junctions from PSL file')
    parser_exon.set_defaults(func=exam_junctions)
    parser_exon.add_argument('full', help='Existing full length list (would be skipped)')
    parser_exon.add_argument('junction', help='PETA junction file')
    parser_exon.add_argument('psl', help='transcript/contigs-to-annotation PSL file')
    
    parser_exon = subparsers.add_parser('simu', help='simulate reads from specific contigs for small debugging')
    parser_exon.set_defaults(func=simu)
    parser_exon.add_argument('transcripts', nargs='+', help='transcript names')
    parser_exon.add_argument('-t', '--annotated', default='/home/carl/Projects/peta/rnaseq/Spombe/genome/spombe.broad.tx.fasta.rev', help='annotated transcript file')
    parser_exon.add_argument('-l', '--read_len', default=68, type=int, help='read length')

    args = parser.parse_args()
    args.func(args)
                
if __name__ == '__main__':
    sys.exit(main())
