from eva import *
from zoom import *
from merge import *
from argparse import ArgumentParser
import operator, random
import math

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

def junc(args):
    junctions = read_junctions(args.jun)
    fa = FastaFile(args.fa)
    joint = FastaFile()
    seq_names = []
    for j in junctions:
        if j.ori == 0:
            name = j.main + '_' + j.branch + '_' + str(j.locus) + '_' + str(j.ori)
            left = fa.seqs[j.main][0:j.locus]
            right = fa.seqs[j.branch]
        else:
            name = j.main + '_' + j.branch + '_' + str(j.branch_len) + '_' + str(j.ori)
            left = fa.seqs[j.branch]
            right = fa.seqs[j.main][j.locus:j.main_len]
        joint.seqs[name] = left + right
        seq_names.append(name);
    
    with open(args.fa[:-3] + '.joint.fa', 'w') as out:
        for tx_name in seq_names:
            seq = joint.seqs[tx_name]
            out.write('>%s length: %d\n' % (tx_name, len(seq)))
            l = 0
            for c in seq:
                out.write(c)
                l += 1
                if l % 50 == 0:
                    out.write('\n')
            if not l % 50 == 0:
                out.write('\n')
    
    print 'Check file %s.joint.fa' % (args.fa[:-3])
    
def collect_last_reads(args):
    paired = FastaFile(args.paired)
    reads = FastaFile(args.reads)
    last_reads = FastaFile()
    junctions = read_junctions(args.junc)
    single_juncs = []
    for this_j in junctions:
        multi = False
        for that_j in junctions:
            if not this_j == that_j and this_j.branch == that_j.branch:
                multi = True
                break
        if not multi:
            single_juncs.append(this_j)
    print single_juncs
        
    i = 0
    ctg_read_dict = {}
    for name, seq in paired.seqs.iteritems():
        ids = name.split('_')
        tpl_id = ids[0]
        last_read_id = ids[1]
        ctg_read_dict[last_read_id] = tpl_id
        if not last_read_id == '0':
            for j in junctions:
                if j.branch == tpl_id:
                    i += 1
                    print i, name
                    read = reads.seqs[last_read_id]
                    last_reads.seqs[last_read_id] = read
    last_reads.save_to_disk(args.paired[:-3] + '.last.fa')
    print 'Check file %s.last.fa' % args.paired[:-3]
    
    cmd = 'blat spombe.fa %s.last.fa %s.last.ref.psl' % (args.paired[:-3], args.paired[:-3])
    print runInShell(cmd)
    hits = read_blat_hits('%s.last.ref.psl' % args.paired[:-3], 'query')
    n_half = 0
    n_single_half = 0
    n_total = 0
    for r in last_reads.seqs.iterkeys():
        n_total += 1
        if not r in hits:
            print 'Read not aligned: %s' % r
        else:
            hs = hits[r]
        for h in hs:
            if h.n_match < h.qlen - 2:
                print h.hit_line
                n_half += 1
                for j in single_juncs:
                    if ctg_read_dict[r] == j.branch:
                        n_single_half += 1
                        break
    print '%d/%d' % (n_half, n_total)
    print '%d/%d' % (n_single_half, len(single_juncs))

def simple_format_junctions(args):
    junctions = read_junctions(args.junc)
    with open(args.junc + '.nolen', 'w') as simple:
        for j in junctions:
            simple.write('%d\t%d\t%d\t%d\t%d\n' % (int(j.main), int(j.branch), int(j.locus), int(j.weight), int(j.ori)))

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
    ids = []
    with open(args.id_list) as lines:
        for line in lines:
            ids.append(line.strip())
    print 'Transcripts to check: %d' % len(ids)
    hits = read_blat_hits(args.psl, 'ref')
    
    n_not_touched = 0
    n_should_be_merged = 0
    n_bad_leaf = 0
    n_has_gap = 0
    n_ol_st11 = 0
    for tx in ids:
        if not tx in hits:
            print '---- %s is not touched at all ----' % tx
            n_not_touched += 1
            continue
        tx_hits = hits[tx]
        tx_hits.sort(key=lambda x: x.rstart)
        pre_h = None
        
        is_bad_leaf = False
        should_be_merged = False
        has_gap = False
        st11 = False
        for h in tx_hits:
            if pre_h is None:
                if h.rstart > 10:
                    break
                pre_h = h
            else:
                dis = pre_h.rend - h.rstart
                if dis > 0:
                    if dis >= 11:
                        if pre_h.qend < pre_h.qlen - 2 or h.qstart >= 2:
                            is_bad_leaf = True
                        else:
                            if h.n_blocks == 1 and pre_h.n_blocks == 1:
                                should_be_merged = True
                    else:
                        st11 = True
                else:
                    dis = abs(dis)
                    has_gap = True
        if st11:
            n_ol_st11 += 1
            print '---- %s some overlap smaller than 11bp ----' % tx
            print tx_hits
        elif has_gap:
            n_has_gap += 1
            print '---- %s has gap ----' % tx
            print tx_hits    
        elif is_bad_leaf:
            n_bad_leaf += 1
            print '---- %s has bad leaf ----' % tx
            print tx_hits
        elif should_be_merged:
            n_should_be_merged += 1
            print '---- %s should be merged ----' % tx
            print tx_hits
    
    print 'Should be merged: %d' % n_should_be_merged
    print 'Has bad leaf: %d' % n_bad_leaf 
    print 'Not touched: %d' % n_not_touched 
    print 'Has gap: %d' % n_has_gap 
    print 'Overlap <11bp: %d' % n_ol_st11
    print 'Others: %d' % (len(ids) - (n_should_be_merged + n_bad_leaf + n_not_touched + n_has_gap + n_ol_st11))
    print 'Totally %d' % (n_should_be_merged + n_bad_leaf + n_not_touched + n_has_gap + n_ol_st11)

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
    if args.ins_size > 0:
        simu_pair(args)
        return
    annotated = FastaFile(args.annotated)
    reads = FastaFile()
    read_id = 0
    out = os.path.join(os.path.dirname(args.annotated), 'simu.fa')
    part = os.path.join(os.path.dirname(args.annotated), 'part.fa')
    part_fa = FastaFile()
    # out = '/home/carl/Projects/peta/rnaseq/Spombe/SRR097897/simu.fa'
    for tx in args.transcripts:
        seq = annotated.seqs[tx]
        part_fa.seqs[tx] = seq
        for i in range(len(seq) - args.read_len + 1):
            reads.seqs[read_id] = seq[i:i+args.read_len]
            read_id += 1
            reads.seqs[read_id] = rev_comp(seq[i:i+args.read_len])
            read_id += 1
    reads.save_to_disk(out)
    part_fa.save_to_disk(part)
    cmd = 'blat %s %s %s.psl' % (args.annotated, out, out)
#    print cmd
#    print runInShell(cmd)
    print 'Check %d reads in file %s' % (read_id, out)
    
def simu_pair(args):
    annotated = FastaFile(args.annotated)
    reads = FastaFile()
    read_id = 0
    out = os.path.join(os.path.dirname(args.annotated), 'simu.pe.fa')
    part = os.path.join(os.path.dirname(args.annotated), 'part.fa')
    part_fa = FastaFile()
    # out = '/home/carl/Projects/peta/rnaseq/Spombe/SRR097897/simu.fa'
    for tx in args.transcripts:
        seq = annotated.seqs[tx]
        part_fa.seqs[tx] = seq
        for i in range(len(seq) - args.read_len + 1):
#            if i > 250:
#                break
            if (i + args.read_len + args.ins_size <= len(seq)):
#                if i > 50 and i < 80:
#                    continue
                reads.seqs[read_id] = seq[i:i+args.read_len]
                read_id += 1
                reads.seqs[read_id] = seq[i+args.ins_size:i+args.ins_size+args.read_len]
                read_id += 1
            else:
                break
            
            reads.seqs[read_id + 1] = rev_comp(seq[i:i+args.read_len])
            read_id += 1
            reads.seqs[read_id - 1] = rev_comp(seq[i+args.ins_size:i+args.ins_size+args.read_len])
            read_id += 1
            
    reads.save_to_disk(out)
    part_fa.save_to_disk(part)
    cmd = 'blat %s %s %s.psl' % (args.annotated, out, out)
#    print cmd
#    print runInShell(cmd)
    print 'Check %d reads in file %s' % (read_id, out)
    
def match(args):
    fa = FastaFile(args.fa)
    count = 0
    print args.seq
    for name, seq in fa.seqs.iteritems():
        try:
            index = args.seq.index(seq)
            print ' ' * index + '%s\t%s' % (seq, name)
            count += 1
        except:
            pass
    print '==== HITS: %d ====' % count
    
def find_splicing(args):
    hits = read_blat_hits(args.psl, 'ref')
    tx_w_splicing = {}
    print '# of transcripts with splicing: %d' % len(hits)
    with open(args.psl + '.ids', 'w') as ids:
        for tx in hits.iterkeys():
            ids.write(tx + '\n')
    print 'Check file %s.ids' % (args.psl)
    for ref, tx_hits in hits.iteritems():
        for h in tx_hits:
            if h.n_ref_gap_bases >= 10 or h.n_query_gap_bases >= 10 or abs(h.rlen - h.qlen) >= 100:
                tx_w_splicing[h.rname] = h
                tx_w_splicing[h.qname] = h
    print '# of transcripts not closely similar: %d' % len(tx_w_splicing)
    with open(args.psl + '.nots.ids', 'w') as ids:
        for tx in tx_w_splicing.iterkeys():
            ids.write(tx + '\n')
    print 'Check file %s.nots.ids' % args.psl

def extract(args):
    ref = FastaFile(args.ref)
    reads = FastaFile(args.reads)
    part_reads = FastaFile()
    part_ref = FastaFile()
    read_id = 0
    read_len = len(reads.seqs['0'])
    saved_ids = []
    for t in args.tx:
        part_ref.seqs[t] = ref.seqs[t]
        cmd = 'grep %s %s.psl' % (t, args.reads)
        raw_lines = runInShell(cmd)
        lines = raw_lines.split('\n')
        tx_hits = read_psl_hits(lines, 'ref')
        hits = tx_hits[t]
        print '%d hits on %s' % (len(hits), t)
        for h in hits:
#            if read_id > 30:
#                break
            if h.qname in reads.seqs and not h.qname in saved_ids:
                #print '>%s' % h.qname
                #print reads.seqs[h.qname]
                # if the read is left mate
                saved_ids.append(h.qname)
                if int(h.qname) % 2 == 0:
                    part_reads.seqs[read_id] = reads.seqs[h.qname]
                    mate_id = str(int(h.qname) + 1)
                    for h2 in hits:
                        if h2.qname == mate_id:
                            #print '>%s' % h2.qname
                            #print reads.seqs[h2.qname]
                            read_id += 1
                            part_reads.seqs[read_id] = reads.seqs[mate_id]
                            break
                    if read_id % 2 == 0:
                        read_id += 1
                        part_reads.seqs[read_id] = 'N' * read_len
                    read_id += 1
                    saved_ids.append(mate_id)
                else:
                    read_id += 1
                    part_reads.seqs[read_id] = reads.seqs[h.qname]
                    mate_id = str(int(h.qname) - 1)
                    has_mate = False
                    for h2 in hits:
                        if h2.qname == mate_id:
                            #print '>%s' % h2.qname
                            #print reads.seqs[h2.qname]
                            part_reads.seqs[read_id - 1] = reads.seqs[mate_id]
                            has_mate = True
                            break
                    if not has_mate:
                        part_reads.seqs[read_id - 1] = 'N' * read_len
                    read_id += 1
                    saved_ids.append(mate_id)
    part_ref.save_to_disk(args.ref + '.part.fa')
    print 'Check %d transcripts in file %s.part.fa' % (len(args.tx), args.ref)
    part_reads.save_to_disk(args.reads[:-3] + '.tx.fa')
    print 'Check %d reads in file %s.tx.fa' % (read_id, args.reads[:-3])

def check_overlap(args):
    hits = read_blat_hits(args.tx2tx, 'ref')
    clustered = {}
    n_cluster = 0
    n_tx = 0
    for tx, tx_hits in hits.iteritems():
        prefix = tx.split('.')[0]
        if tx in clustered.iterkeys(): continue
        has_ol = False
        for h in tx_hits:
            if h.qname == tx: continue
            p = h.qname.split('.')[0]
            if p == prefix: continue
            for i in range(h.n_blocks):
                if h.block_sizes[i] > 68:
                    if not has_ol:
                        print '%d: ----------- %s ------------' % (n_cluster, tx)
                    has_ol = True
                    print tx, h.qname
                    clustered[tx] = True
                    clustered[h.qname] = True
        if has_ol: n_cluster += 1
    print '# of cluster: %d' % n_cluster
    print '# of transcripts clustered: %d' % len(clustered.keys())

def get_singleton(args):
    txtx_hits = read_blat_hits(args.tx2tx, 'ref')
    with open(args.out, 'w') as out:
        for tx, tx_hits in txtx_hits.iteritems():
            if len(tx_hits) == 1:
                out.write(tx + '\n')
    print 'Check file %s' % args.out
    
def continuous_paired(args):
    ref = FastaFile(args.ref)
    ids = []
    with open(args.ids) as input:
        for line in input:
            ids.append(line.strip())
    n_100 = 0
    n_200 = 0
    print 'Totally %d ids' % len(ids)
    for tx in ids:
#         tx = 'SPAC4F8.08_T0'
        seq = ref.seqs[tx]
        cmd = 'grep %s %s' % (tx, args.fa2tx)
        hit_lines = runInShell(cmd)
        lines = hit_lines.split('\n')
        if len(lines) == 0: continue
        hits = read_psl_hits(lines, 'ref')
        if not tx in hits: continue
        tx_hits = hits[tx]
        print '------------ %20s length %d: %d hits -----------' % (tx, len(seq), len(lines))
        tx_hits.sort(key=lambda x:x.rstart)
        reads = []
        for h in tx_hits:
            reads.append(h.qname)
        bases = [0 for _ in range(len(seq))]
        
        for h in tx_hits:
            read_id = int(h.qname)
            mname = read_id
            if read_id % 2 == 0:
                mname = str(read_id + 1)
            else:
                mname = str(read_id - 1)
            if mname in reads:
                end = h.rstart + h.qlen
                if end > h.rlen: end = h.rlen
                for i in range(h.rstart, end):
                    bases[i] = 1
        
        sum = 0
        for b in bases:
            sum += b
        print 'Sum covered: %d' % sum
        
        max = 0
        not_paired_len = 0
        start_0 = 0
        end_0 = 0
        cursor = 0
        pre = bases[0]
        for i in range(1, len(bases)):
            b = bases[i]
#             print '%d: %d' % (i, b)
            if b == 0 and bases[pre] == 1:
                cursor = i    
            if (b == 1 and bases[pre] == 0) or (i == len(bases) - 1 and b == 0):
                if  i - cursor > max:
                    start_0 = cursor
                    end_0 = i
                    max = end_0 - start_0
            pre = i
#         if max > 200:
        print '%s: [%d, %d]' % (tx, start_0, end_0)
        if max > 100:
            n_100 += 1
        if max > 200:
            n_200 += 1
#         break
    print '>100bp: %d' % n_100
    print '>200bp: %d' % n_200
    
def avg(l):
    return float(sum(l)) / float(len(l))

def ins_size(args):
    ids = []
    with open(args.ids) as f:
        for line in f:
            ids.append(line.strip())
    sizes = []
    categories = {}
    for i in range(1000):
        categories[i] = 0
    categories[1000] = 0
    c = 0
    for id in ids:
        c += 1
        print >>sys.stderr, '%d -------------- %s -------------' % (c, id)
        cmd = 'grep %s %s' % (id, args.read2ref)
        print >>sys.stderr, cmd
        hit_lines = runInShell(cmd).split("\n")
        hits = read_psl_hits(hit_lines, 'ref')
        if not id in hits:
            continue
        tx_hits = hits[id]
        tx_hits.sort(key=lambda x:int(x.qname))
        pre = None
        n = 0
        for h in tx_hits:
            if pre:
                if int(h.qname) % 2 == 1 and int(h.qname) - int(pre.qname) == 1:
                    s = abs(h.rstart - pre.rstart)
                    #print id, h.qname, pre.qname, 
                    #print abs(h.rstart - pre.rstart)
                    sizes.append(s)
                    n += 2
                    if s >= 1000:
                        categories[1000] += 1
                    else:
                        categories[s] += 1
            pre = h
        print '%s\t%.2f' % (id, float(n) / float(len(tx_hits)))
        #break
    print categories
    a = avg(sizes)
    print >>sys.stderr, 'Average: %.2f' % a
    variance = map(lambda x: (x - a)**2, sizes)
    standard_deviation = math.sqrt(avg(variance))
    print >>sys.stderr, 'Standard deviation: %.2f' % standard_deviation

def oracle(args):
    hits = read_blat_hits(args.read2ref, 'ref')
    for tx, tx_hits in hits.iteritems():
        tx_hits.sort(key=lambda x:x.rstart)

def get_mate_id(name):
    i = int(name)
    if i % 2 == 0:
        return str(i + 1)
    else:
        return str(i - 1)

def mismatch(args):
    fasta = FastaFile(args.ref)
    all_tx_hits = read_psl_hits(args.read2ref, 'ref')
    all_read_hits = read_psl_hits(args.read2ref, 'query')
    n_bad_but_from_another = []
    n_bad_with_mismatches = []
    n_bad_pairs = 0
    n_bad_singles = 0
    n_bad_total = 0
    n_bad_cross_ol = []
    counter = 0
    for tx, seq in fasta.seqs.iteritems():
        counter += 1
        if not tx in all_tx_hits: continue
        tx_hits = all_tx_hits[tx]
        print '%d hits' % len(tx_hits)
        tx_hits.sort(key=lambda x:x.rstart, reverse=False)
        pre = None
        nex = None
        for i in range(len(tx_hits)):
            h = tx_hits[i]
            if i > 0: pre = tx_hits[i - 1]
            if h.n_match + 2 >= h.qlen: continue
            if not h.n_match + h.n_mismatch == h.qlen: continue
            r_hits = all_read_hits[h.qname]
            r_hits.sort(key=lambda x:x.n_match)
            first_h = r_hits[0]
            if first_h.n_match + 2 >= h.qlen:
                n_bad_but_from_another.append(h)
            else:
                left_ol = -1 
                right_ol = -1
                cross_ol = -1
                mate_name = get_mate_id(h.qname)
                
                mate_pos = -1
                if mate_name in all_read_hits:
                    mate_hits = all_read_hits[mate_name]
                    mate_hits.sort(key=lambda x:x.n_match, reverse=True)
                    for mate_h in mate_hits:
                        if mate_h.rname == h.rname:
                            mate_pos = mate_h.rstart
                            break
                if pre:
                    left_ol = pre.rend - h.rstart
                if i < len(tx_hits) - 1:
                    nex = tx_hits[i + 1]
                    right_ol = h.rend - nex.rstart
                if pre and nex:
                    #print pre
                    #print nex
                    cross_ol = pre.rend - nex.rstart
                if mate_pos == -1:
                    n_bad_singles += 1
                else:
                    n_bad_pairs += 1
                n_bad_with_mismatches.append((h, left_ol, right_ol, cross_ol))
                n_bad_total += 1
                if not left_ol == -1 and not right_ol == -1 and abs(cross_ol) <= 25:
                    n_bad_cross_ol.append((tx, h.qname, h.n_mismatch, cross_ol))
                print 'Read %s on %s \t %d mismatches: \t Left %d \t Right %d \t Cross %d \t Pair [%d, %d] %d' % (h.qname, h.rname, h.n_mismatch, left_ol, right_ol, cross_ol, h.rstart, mate_pos, h.rstart - mate_pos)
        #break
    print 'Reads that >2 mismatches, but on another transcript with <=2 mismatches: %d' % (len(n_bad_with_mismatches))
    print 'Reads whose best alignment have >2 mismatches: %d' % (len(n_bad_with_mismatches))
    print 'Bad paired reads: %d' % n_bad_pairs
    print 'Bad single reads: %d' % n_bad_singles
    print 'Bad cross overlap: %d' % len(n_bad_cross_ol)
    for (tx, h.qname, h.n_mismatch, cross_ol) in n_bad_cross_ol:
        print tx, h.qname, h.n_mismatch, cross_ol
    #for (h, left_ol, right_ol, cross_ol) in n_bad_with_mismatches:
    #    print 'Read %s on %s \t %d mismatches \t Left %d \t Right %d \t Cross %d' % (h.qname, h.rname, left_ol, right_ol, cross_ol)

def qc(args):
    tx = FastaFile(args.tx)
    print 'Building read-to-hit hash...'
    readHitHash = [[] for _ in range(7481042)]
    line_no = 0
    with open(args.read2ref) as read2ref:
        for line in read2ref:
            line_no += 1
            if line_no <= 5: continue
            fs = line.split('\t')
            if len(fs) < 9: continue
            readHitHash[int(fs[9])].append(line.strip())
    idList = []
    with open(args.ids) as ids:
        for line in ids:
            idList.append(line.strip())
    print 'Reading transcript to genome hits...'
    tx2refHits = read_blat_hits(args.tx2ref, 'query')
    print 'Transcript\tLength\tHits\tReads\tPaired reads\tPaired reads >2 mismatches\tPaired reads >2 soft-clipped to introns\tMates not aligned\tAvg matches of not aligned\tMates at UTR\tMates far away'
    for tx_name in idList:
        out = tx_name
        cmd = 'grep %s %s' % (tx_name, args.read2tx)
        rawLines = runInShell(cmd)
        lines = rawLines.split('\n')
        if len(lines) < 2: continue
        l = lines[0]
        fs = l.split('\t')
        out = '%s\t%s' % (out, fs[14])
        out = '%s\t%d' % (out, len(lines) - 1)
        tx_hits = read_psl_hits(lines, 'query')
        out = '%s\t%d' % (out, len(tx_hits.keys()))
        pairs = {}
        pairs_but_bad_mis = {}
        pairs_but_bad_mis_intron = {}
        # If pairs are on the same transcripts, they are good
        for qname, read_hits in tx_hits.iteritems():
            mname = get_mate_id(qname)
            if mname in tx_hits:
                pairs[qname] = mname
                read_hits.sort(key=lambda x:x.n_match, reverse=True)
                h = read_hits[0]
                if h.n_match < h.qlen - 2:
                    pairs_but_bad_mis[qname] = h
                    # If no soft clipping, skip
                    if h.n_match + h.n_mismatch == h.qlen: continue
                    
                    ref_hits = tx2refHits[tx_name]
                    ref_hits.sort(key=lambda x:x.n_match, reverse=True)
                    h = ref_hits[0]
                    lines = readHitHash[int(qname)]
                    readRefHits = read_psl_hits(lines, 'ref')
                    if h.rname in readRefHits:
                        read_ref_hits = readRefHits[h.rname]
                        for rh in read_ref_hits:
                            if rh.rstart > h.rstart and rh.r_block_starts[rh.n_blocks - 1] + rh.block_sizes[rh.n_blocks - 1] < h.r_block_starts[h.n_blocks - 1] + h.block_sizes[h.n_blocks - 1]:
                                if rh.n_match >= rh.qlen - 2:
                                    pairs_but_bad_mis_intron[qname] = h
        out = '%s\t%d' % (out, len(pairs))
        out = '%s\t%d' % (out, len(pairs_but_bad_mis.keys()))
        out = '%s\t%d' % (out, len(pairs_but_bad_mis_intron.keys()))
        
        # Those reads whose mates are not aligned to the same transcript
        n_not_aligned = []
        n_mate_at_utr = []
        n_mate_far = []
        not_aligned_matches = []
        for qname, readHits in tx_hits.iteritems():
            if qname in pairs: continue
            mname = get_mate_id(qname)
            lines = readHitHash[int(mname)]
            if len(lines) == 0:
                n_not_aligned.append(mname)
                readHits.sort(key=lambda x: x.n_match, reverse=True)
                h = readHits[0]
                not_aligned_matches.append(h.n_match)
            else:
                ref_hits = tx2refHits[tx_name]
                ref_hits.sort(key=lambda x:x.n_match, reverse=True)
                h = ref_hits[0]
                startTxOnRef = h.rstart - 1000
                endTxOnRef = h.r_block_starts[h.n_blocks - 1] + h.block_sizes[h.n_blocks - 1] + 1000
                read_hits = read_psl_hits(lines, 'ref')
                if h.rname in read_hits:
                    hits = read_hits[h.rname]
                    atUtr = False
                    for mh in hits:
                        if mh.rstart >= startTxOnRef and mh.rstart <= endTxOnRef:
                            n_mate_at_utr.append(mname)
                            atUtr = True
                            #print '%s at %s: %d~%d' % (tx_name, h.rname, h.rstart, h.r_block_starts[h.n_blocks - 1] + h.block_sizes[h.n_blocks - 1])
                            #print '%s at %s: %d M%d -%d' % (mname, h.rname, mh.rstart, mh.n_match, mh.n_mismatch)
                            break
                    if not atUtr:
                        n_mate_far.append(mname)
                else:
                    n_mate_far.append(mname)
                
        out = '%s\t%d' % (out, len(n_not_aligned))
        avg_matches = 68
        if len(not_aligned_matches) > 0: avg_matches = avg(not_aligned_matches)
        out = '%s\t%.2f' % (out, avg_matches)
        out = '%s\t%d' % (out, len(n_mate_at_utr))
        out = '%s\t%d' % (out, len(n_mate_far))
        print out

def resave(args):
    fasta = FastaFile(args.fa)
    fasta.save_to_disk(args.fa + '.50')
    print 'Check file %s.50' % (args.fa)
    
def gene(args):
    fasta = FastaFile(args.tx)
    genes = {}
    for tx, seq in fasta.seqs.iteritems():
        ids = tx.split('_')
        if ids[0] in genes: genes[ids[0]].append(tx)
        else: genes[ids[0]] = [tx]
        
    for gene, transcripts in genes.iteritems():
        print '%s\t%d' % (gene, len(transcripts))
        
def extend(args):
    fasta = FastaFile(args.tx)
    genome = FastaFile(args.ref)
    tx_hits = read_blat_hits(args.psl, 'query')
    extended = FastaFile()
    for tx, hits in tx_hits.iteritems():
        hits.sort(key=lambda x: x.alen, reverse=True)
        if len(hits) <= 0: 
            print '%s not aligned to reference genome' % tx
            continue
        h = hits[0]
        seq = fasta.seqs[tx]
        if h.strand == '-': seq = rev_comp(seq)
        chr_seq = genome.seqs[h.rname]
        s = h.rstart - 1000
        if s < 0: s = 0
        ext_seq = chr_seq[s:h.rstart]
        for i in range(h.n_blocks):
            ext_seq += chr_seq[h.r_block_starts[i]:(h.r_block_starts[i] + h.block_sizes[i])]
        e = h.rend + 1000
        if e > len(chr_seq): e = len(chr_seq)
        ext_seq += chr_seq[h.rend:e]
        extended.seqs[tx] = ext_seq.upper()
    extended.save_to_disk(args.tx + '.extended')

def main():
    parser = ArgumentParser()
    subparsers = parser.add_subparsers(help='sub command help')
    
    parser_mis = subparsers.add_parser('mismatch', help='Pick reads with many mismatches')
    parser_mis.set_defaults(func=mismatch)
    parser_mis.add_argument('ref', help='reference transcript')
    parser_mis.add_argument('read2ref', help='read2ref PSL')
    
    parser_ins = subparsers.add_parser('ins', help='Calculate insert size')
    parser_ins.set_defaults(func=ins_size)
    parser_ins.add_argument('ids', help='singleton id list')
    parser_ins.add_argument('read2ref', help='read2ref PSL')
    
    parser_oracle = subparsers.add_parser('oracle', help='Get oracle set')
    parser_oracle.set_defaults(func=oracle)
    parser_oracle.add_argument('ref', help='Reference transcripts')
    parser_oracle.add_argument('read2ref', help='read2ref PSL')
    
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
    
    parser_match = subparsers.add_parser('match', help='simply reads to a short sequence')
    parser_match.set_defaults(func=match)
    parser_match.add_argument('fa', help='Fasta read file')
    parser_match.add_argument('seq', help='Sequence')
    
    parser_junction = subparsers.add_parser('junction', help='Exam junctions from PSL file')
    parser_junction.set_defaults(func=exam_junctions)
    parser_junction.add_argument('id_list', help='List of ids to inspect')
    parser_junction.add_argument('psl', help='transcript/contigs-to-annotation PSL file')
    
    parser_simple = subparsers.add_parser('simple', help='Convert junctions to simpler format')
    parser_simple.set_defaults(func=simple_format_junctions)
    parser_simple.add_argument('junc', help='somefile.junctions')
    
    parser_exon = subparsers.add_parser('simu', help='simulate reads from specific contigs for small debugging')
    parser_exon.set_defaults(func=simu)
    parser_exon.add_argument('transcripts', nargs='+', help='transcript names')
    parser_exon.add_argument('-t', '--annotated', default='/home/carl/Projects/peta/rnaseq/Spombe/genome/spombe.broad.tx.fasta.rev', help='annotated transcript file')
    parser_exon.add_argument('-l', '--read_len', default=68, type=int, help='read length')
    parser_exon.add_argument('-i', '--ins_size', required=False, default=0, type=int, help='insert size')

    parser_junc = subparsers.add_parser('junc', help='Check whether the junction templates are aligned to the same transcript')
    parser_junc.set_defaults(func=junc)
    parser_junc.add_argument('fa', help='paired.fa')
    parser_junc.add_argument('jun', help='paired.junctions file')
    
    parser_splice = subparsers.add_parser('splice', help='Get transcripts with splicing')
    parser_splice.set_defaults(func=find_splicing)
    parser_splice.add_argument('psl', help='PSL after removing self-to-self alignments')
    
    parser_extract = subparsers.add_parser('extract', help='Extract hits on transcripts')
    parser_extract.set_defaults(func=extract)
    parser_extract.add_argument('ref', help='transcript fasta')
    parser_extract.add_argument('reads', help='reads fasta')
    parser_extract.add_argument("tx", nargs='+', help="[transcript ids]")

    parser_last = subparsers.add_parser('last', help='Get last read of branch templates')
    parser_last.set_defaults(func=collect_last_reads)
    parser_last.add_argument('paired', help='paired.fa')
    parser_last.add_argument('junc', help='paired.junctions')
    parser_last.add_argument('reads', help='reads.fa')
    
    parser_ol = subparsers.add_parser('ol', help='Check overlap between S.pombe transcripts')
    parser_ol.set_defaults(func=check_overlap)
    parser_ol.add_argument('tx', help='transcript fasta')
    parser_ol.add_argument('tx2tx', help='tx-to-tx PSL')
    
    parser_singleton = subparsers.add_parser('single', help='get singleton transcripts')
    parser_singleton.set_defaults(func=get_singleton)
    parser_singleton.add_argument('tx2tx', help='tx-to-tx PSL')
    parser_singleton.add_argument('out', help='output file')
    
    parser_continuous = subparsers.add_parser('blank', help='get max region where no pairs')
    parser_continuous.set_defaults(func=continuous_paired)
    parser_continuous.add_argument('ids', help='ids file')
    parser_continuous.add_argument('ref', help='transcript file')
    parser_continuous.add_argument('fa2tx', help='read-to-tx PSL')
    
    parser_qc = subparsers.add_parser('qc', help='Quality control')
    parser_qc.set_defaults(func=qc)
    parser_qc.add_argument('ids', help='ID list')
    parser_qc.add_argument('tx', help='transcript FASTA')
    parser_qc.add_argument('read2tx', help='reads-to-tx PSL')
    parser_qc.add_argument('read2ref', help='reads-to-genome PSL')
    parser_qc.add_argument('tx2ref', help='transcripts-to-genome PSL')
    
    parser_resave = subparsers.add_parser('resave', help='With different base pair length')
    parser_resave.set_defaults(func=resave)
    parser_resave.add_argument('fa', help='fa file')
    
    parser_gene = subparsers.add_parser('gene', help='Inspect gene-transcript multiplicity')
    parser_gene.set_defaults(func=gene)
    parser_gene.add_argument('tx', help='transcript fasta file')
    
    parser_exd = subparsers.add_parser('extend', help='Add 1000bp to left/right side of a transcript from genome')
    parser_exd.set_defaults(func=extend)
    parser_exd.add_argument('tx', help='transcript FASTA file')
    parser_exd.add_argument('ref', help='genome FASTA file')
    parser_exd.add_argument('psl', help='transcript-to-genome PSL file')
    
    args = parser.parse_args()
    args.func(args)
                
if __name__ == '__main__':
    sys.exit(main())
