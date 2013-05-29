import sys, os, pysam
from subprocess import Popen, PIPE
import eva, merge
from eva import *
from argparse import ArgumentParser

SEP = '    '

BLAT = '/home/carl/Projects/blat/blat'
BLAT_OCC_11 = '/home/carl/Projects/blat/11.ooc'

SPOMBE_GENOME = '/home/carl/Projects/peta/rnaseq/Spombe/genome/spombe.fa'

READ_REF_PSL = '/home/carl/Projects/peta/rnaseq/hg19/SRX011545/SRR027876.psl'
REF_TO_REF = '/home/carl/Projects/peta/rnaseq/hg19/genome/ref.ref.psl'
REF = '/home/carl/Projects/peta/rnaseq/hg19/genome/human.ensembl.cdna.fa'
READ = '/home/carl/Projects/peta/rnaseq/hg19/SRX011545/SRR027876.fa'
CONTIG = '/home/carl/Projects/peta/SRR027876_out/pair_contigs.fa'

READ_REF_PSL_SRR097897 = '/home/carl/Projects/peta/rnaseq/Spombe/SRR097897/SRR097897.rev.psl'
REF_TO_REF_SRR097897 = '/home/carl/Projects/peta/rnaseq/Spombe/genome/ref.ref.rev.psl'
REF_SRR097897 = '/home/carl/Projects/peta/rnaseq/Spombe/genome/spombe.broad.tx.fasta.rev'
READ_SRR097897 = '/home/carl/Projects/peta/rnaseq/Spombe/SRR097897/SRR097897.fa'

READ_SRR027876 = '/home/carl/Projects/peta/rnaseq/hg19/SRX011545/SRR027876.fa'
REF_SRR027876 = '/home/carl/Projects/peta/rnaseq/hg19/genome/human.ensembl.cdna.fa.oracle'
REF_TO_REF_SRR027876 = '/home/carl/Projects/peta/rnaseq/hg19/genome/ref.ref.psl'
READ_REF_PSL_SRR027876 = '/home/carl/Projects/peta/rnaseq/hg19/SRX011545/SRR027876.fa.psl'
HIT_STR_DESCRIPTION = '[read_id]:[read_length]    [qstart]@[rstart]M[n_matches],Ins[n_ref_gap_bases][strand][n_mismatches]'

def runInShell(cmd):
    p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    return p.communicate()[0]
	
def draw_all_dots(args):
	txs = FastaFile(REF_SRR097897)
	for tx_name, seq in txs.seqs.iteritems():
		print tx_name
		cmd = 'py zoom.py draw %s' % tx_name
		runInShell(cmd)
		cmd = 'dot -Tpng %s.dot > ../dots/%s.png' % (tx_name, tx_name)
		runInShell(cmd)
		cmd = 'mv %s.dot ../dots/' % tx_name
		runInShell(cmd)
        
def visu_all_reads_to_ref(args):
    genome = FastaFile(SPOMBE_GENOME)
    base_in_exon = {}
    base_in_exon['chr1'] = [0 for _ in range(genome.get_seq_len('chr1'))]
    base_in_exon['chr2'] = [0 for _ in range(genome.get_seq_len('chr2'))]
    base_in_exon['chr3'] = [0 for _ in range(genome.get_seq_len('chr3'))]
    
    cmd = 'tail -n +6 %s' % args.tx_psl
    print cmd
    lines = runInShell(cmd)
    hit_lines = lines.split('\n')
    hits = eva.read_psl_hits(hit_lines, 'ref')
    for chr_index in range(3):
        chr = 'chr%d' % (chr_index + 1)
        for h in hits[chr]:
            for i in range(h.n_blocks):
                b_start = h.r_block_starts[i]
                for j in range(h.block_sizes[i]):
                    base_in_exon[chr][b_start + j] = 1
    
    cmd = 'tail -n +6 %s' % args.read_psl
    print cmd
    lines = runInShell(cmd)
    hit_lines = lines.split('\n')
    read_hits = eva.read_psl_hits(hit_lines, 'query')
    reads = FastaFile(READ_SRR097897)
    for read_id in reads.seqs.iterkeys():
        hits = read_hits[read_id]
        print 'Read: %s' % read_id
        print hits
        pos_in_exon = [0 for _ in range(reads.get_seq_len(read_id))]
        for h in hits:
            for i in range(h.n_blocks):
                b_start = h.q_block_starts[i]
                for j in range(h.block_sizes[i]):
                    pos_in_exon[j + b_start] = 1
        print pos_in_exon
        break

def sub_read(fasta, query):
    query_read_psl = query + '.psl'
    cmd = '%s %s %s -ooc=%s %s' % (BLAT, fasta, query, BLAT_OCC_11, query_read_psl)
    
def ctg_to_ref(args):
    ref = args.tx
    if args.transcript:
        hits_file_name = args.transcript + '.ctg.hits'
    else:
        hits_file_name = contig + '.ctg.hits'
    hits_file = open(hits_file_name, 'w')
    hits_file.write(HIT_STR_DESCRIPTION + '\n')
    for contig in args.contigs:
        if contig.endswith('fa') or contig.endswith('fasta'):
            print 'Aligning %s to %s ' % (contig, ref)
            ctg_ref_psl = contig + '.psl'
            cmd = '%s %s %s -ooc=%s %s' % (BLAT, ref, contig, BLAT_OCC_11, ctg_ref_psl)
            print cmd
            runInShell(cmd)
        else:
            ctg_ref_psl = contig
        summary, hits = zoom_tx(args.transcript, ref, ctg_ref_psl, 'ctg')
        if len(hits) <= 0:
            continue
        ctg_fa = contig[0:-4]
        hits_file.write(get_align_str(ref, ctg_fa, hits))
        hits_file.write('+' * (hits[0].rlen) + '\n')
    hits_file.close()
    print 'Check file %s' % (hits_file_name)
        
def read_to_ref(args):
    tx_name = args.transcript
    ref = args.tx
    summary, tx_hits = zoom_tx(tx_name, ref, args.psl, 'read')
    align_str = get_align_str(ref, args.reads, tx_hits)
    hits_f = open(tx_name + '.reads.hits', 'w')
    hits_f.write(align_str)
    hits_f.close()
    print summary
    print 'Check text alignment at file %s.reads.hits' % tx_name
    
def pair_to_ref(args):
    tx_name = args.transcript
    ref = args.tx
    summary, tx_hits = zoom_tx(tx_name, ref, args.psl, 'read')
    align_str = get_pair_str(ref, args.reads, tx_hits)
    hits_f = open(tx_name + '.pairs.hits', 'w')
    hits_f.write(align_str)
    hits_f.close()
    print summary
    print 'Check text alignment at file %s.pairs.hits' % tx_name

def zoom_tx(tx_name, ref, blat_psl, ctg_or_read='read'):
    lines = runInShell('grep ' + tx_name + ' ' + blat_psl)
    hit_lines = lines.split('\n')
    tx_hits = []
    if len(hit_lines) <= 0:
        return tx_hits
    hits = eva.read_psl_hits(hit_lines, 'ref')
    tx = FastaFile(ref)
    tx_seq = tx.seqs[tx_name]
    for rname, hs in hits.iteritems():
        if rname == tx_name:
            tx_hits = hs
    summary = ''    
    if ctg_or_read == 'read':
        n_reads = len(tx_hits)
        n_pairs = 0
        n_match = 0
        n_match_minus_1 = 0
        n_match_minus_2 = 0
        n_match_minus_3 = 0
        n_match_minus_4 = 0
        n_match_minus_5 = 0
        n_match_lt_5 = 0
        n_lt_one_block = 0
        all_base_covered = False
        cvr = [0 for x in range(len(tx_seq))]
        ave_cvr = 0.0
        sd_cvr = 0.0
        
        tx_hits.sort(key=lambda x: x.qname, reverse=False)
        h_pre = None
        for h in tx_hits:
            if h.qlen - h.n_match == 0:
                n_match += 1
            if h.qlen - h.n_match == 1:
                n_match_minus_1 += 1
            if h.qlen - h.n_match == 2:
                n_match_minus_2 += 1
            if h.qlen - h.n_match == 3:
                n_match_minus_3 += 1
            if h.qlen - h.n_match == 4:
                n_match_minus_4 += 1
            if h.qlen - h.n_match == 5:
                n_match_minus_5 += 1
            if h.qlen - h.n_match > 5:
                n_match_lt_5 += 1
            if h.n_blocks > 1:
                n_lt_one_block += 1
            try:
                if not h_pre is None:
                    if int(h.qname) - int(h_pre.qname) == 1:
                        n_pairs += 1
            except:
                pass
            h_pre = h
        if len(tx_hits) > 0:
            summary += 'transcript:       %s\n' % tx_name
            summary += 'tx length:        %s\n' % len(tx_seq)
            summary += '# of reads:       %d\n' % n_reads
            summary += '# of pairs:       2 * %d\n' % n_pairs
            summary += 'Full match:       %d\n' % n_match
            summary += '1 mismatch:       %d\n' % n_match_minus_1
            summary += '2 mismatch:       %d\n' % n_match_minus_2
            summary += '3 mismatch:       %d\n' % n_match_minus_3
            summary += '4 mismatch:       %d\n' % n_match_minus_4
            summary += '5 mismatch:       %d\n' % n_match_minus_5
            summary += '>5 mismatch:      %d\n' % n_match_lt_5
            summary += '>1 blocks:        %d\n' % n_lt_one_block
    
    tx_hits.sort(key=lambda x: x.n_blocks, reverse=False)
    tx_hits.sort(key=lambda x: x.rstart, reverse=False)
    return summary, tx_hits

'''
To get a string representing the hit
    hit: a BlastHit instance
    reads: a FastaFile instance with the read sequence
'''
def get_hit_str(hit, reads):
    hit_str = ''
    pre_end = 0
    read_seq = reads.seqs[hit.qname].upper()
    rev_read_seq = merge.rev_comp(read_seq).upper()
    ins_strs = []
    for i in range(0, hit.n_blocks):
        if pre_end > 0:
            hit_str += '-' * (hit.r_block_starts[i] - pre_end)
        if hit.strand == '+':
            hit_str += read_seq[hit.q_block_starts[i]:hit.q_block_starts[i] + hit.block_sizes[i]] 
        else:
            hit_str += rev_read_seq[hit.q_block_starts[i]:hit.q_block_starts[i] + hit.block_sizes[i]] 
        pre_end = hit.r_block_starts[i] + hit.block_sizes[i]
    if (hit.qend < hit.qlen):
        hit_str += '=>'
        if hit.strand == '+':
            hit_str += read_seq[hit.qend:hit.qlen]
        else:
            hit_str += rev_read_seq[hit.qend:hit.qlen]
    hit_str += SEP + hit.qname + ':' + str(len(read_seq)) + ' [' + str(hit.qstart) + '@' + str(hit.rstart) + 'M' + str(hit.n_match) + ',Ins' + str(hit.n_ref_gap_bases) + str(hit.strand) + str(hit.n_mismatch) + ']'
    return hit_str

def get_pair_hit_str(mate_h_1, mate_h_2, reads):
    if mate_h_1.rstart < mate_h_2.rstart:
        left_mate = mate_h_1
        right_mate = mate_h_2
    else:
        left_mate = mate_h_2
        right_mate = mate_h_1
    left_str = get_hit_str(left_mate, reads)
    right_str = get_hit_str(right_mate, reads)
    dis = right_mate.rstart - left_mate.rstart
    left_split = left_str.split('    ')
    right_split = right_str.split('    ')
    
    if left_mate.strand == '+':
        left_ori = '>'
    else:
        left_ori = '<'
    if right_mate.strand == '+':
        right_ori = '>'
    else:
        right_ori = '<'
    
    if dis < left_mate.qlen:
        align_str = ' ' * left_mate.rstart + left_str + SEP + 'OL Pair\n'
        align_str += ' ' * right_mate.rstart + right_str + SEP + 'OL Pair\n'
    elif dis > len(left_str) + 10:
        align_str = ' ' * left_mate.rstart + left_str + SEP + left_ori
        align_str += '=' * (right_mate.rstart - len(align_str) - 1 - len(SEP)) + right_ori + SEP
        align_str += right_str + '\n'
    else:
        align_str = ' ' * left_mate.rstart + left_split[0] + SEP + left_ori
        align_str += '=' * (right_mate.rstart - len(align_str) - 1 - len(SEP)) + right_ori + SEP
        align_str += right_split[0]
        align_str += SEP + left_split[1] + SEP + right_split[1] + '\n'
    return align_str

def get_mate_id(id):
    if int(id) % 2 == 0:
        return str(int(id) + 1)
    else:
        return str(int(id) - 1)

'''
Align the reads to the transcript/contig
Pairs are connected and displayed together
    ref:    transcript/contig file name
    fasta:  reads file name
    hits:   BlastHit instances, hits to display
'''
def get_pair_str(ref, fasta, hits):
    if len(hits) <= 0:
        return ''
    tx_name = hits[0].rname
    algin_str = tx_name + '\n'
    tx = FastaFile(ref)
    reads = FastaFile(fasta)
    seq = tx.seqs[tx_name].upper()
    read_no = 0
    for h in hits:
        if read_no % 15 == 0:
            algin_str += seq + '\n'
        mate_id = get_mate_id(h.qname)
        mate_hit = None
        for m in hits:
            if m.qname == mate_id:
                mate_hit = m
                break
        if mate_hit is None:
            hit_str = ' ' * h.rstart + get_hit_str(h, reads) + '\n'
        else:
            hit_str = get_pair_hit_str(h, mate_hit, reads)
        algin_str += hit_str
        read_no += 1
    return algin_str

# For a list of hits, get the text alignments
def get_align_str(ref, fasta, hits):
    if len(hits) <= 0:
        return ''
    tx_name = hits[0].rname
    algin_str = tx_name + '\n'
    tx = FastaFile(ref)
    reads = FastaFile(fasta)
    seq = tx.seqs[tx_name].upper()
    read_no = 0
    for h in hits:
        if read_no % 15 == 0:
            algin_str += seq + '\n'
        hit_str = get_hit_str(h, reads)
        algin_str += ' ' * h.rstart + hit_str + '\n'
        read_no += 1
    return algin_str

class Block(object):
    def __init__(self, tx, id, start, end, is_covered, matched_block=None):
        self.tx = tx
        self.id = id
        self.start = start
        self.end = end
        self.is_covered = is_covered
        self.matched_block = matched_block
    def __repr__(self):
        s = '%s, %d: %d~%d, %r' % (self.tx, self.id, self.start, self.end, self.is_covered)
        return s

def get_color(transcript, tx, b):
    if not b.is_covered:
        return 'gray'
    if tx == transcript:
        return 'yellow'
    else:
        return 'white'

def nice_name(name):
    return name.replace('.', '_')

def check_covered(tx, hits, start, end):
    for h in hits:
        if h.qname == tx and not h.qname == h.rname:
            for i in range(h.n_blocks):
                b_start = h.q_block_starts[i]
                b_end = h.q_block_starts[i] + h.block_sizes[i]
                if b_start <= start and b_end >= end:
                    return True
    return False

def draw_dot(args):
    tx = FastaFile(args.tx)
    dot = open(args.transcript + '.dot', 'w')
    milestones = []
    
    lines = runInShell('grep ' + args.transcript + ' ' + args.psl)
    hit_lines = lines.split('\n')
    print 'grep ' + args.transcript + ' ' + args.psl
    raw_hits = eva.read_psl_hits(hit_lines, 'query')
    hits = raw_hits[args.transcript]
    
    # Every starting/ending point on the query transcript is a 'milestone'
    m = []
    print hits
    print '-----------------------------' 
    for h in hits:
        if h.qname == args.transcript and not h.qname == h.rname:
            for i in range(h.n_blocks):
                if not h.q_block_starts[i] in m:
                    milestones.append((h.q_block_starts[i], 'start', h.rname, i))
                    m.append(h.q_block_starts[i])
                if not (h.q_block_starts[i] + h.block_sizes[i]) in m:
                    milestones.append((h.q_block_starts[i] + h.block_sizes[i], 'end', h.rname, i))
                    m.append(h.q_block_starts[i] + h.block_sizes[i])
                     
    blocks = {}
    tx_cursors = {}
    tx_node_ids = {}
    milestones.sort()
    blocks[args.transcript] = []
    tx_node_ids[args.transcript] = 0
    tx_cursors[args.transcript] = 0
    dot.write('digraph g { \n\trankdir = LR \n')
    n_solid_blocks = 0
    # Determine all blocks of the query transcript
    print milestones
    for i in range(len(milestones) - 1):
        (next_m, next_start_or_end, next_rname, next_index) = milestones[i + 1]
        (m, start_or_end, rname, index) = milestones[i]
        is_covered = True
        if i == 0 and m > 0:
            b = Block(args.transcript, -1, 0, m, False)
            blocks[args.transcript].append(b)
        is_covered = check_covered(args.transcript, hits, m, next_m)
        if is_covered:
            b = Block(args.transcript, n_solid_blocks, m, next_m, is_covered)
            b.matched_block = b
            n_solid_blocks += 1
        else:
            b = Block(args.transcript, -1, m, next_m, is_covered)
        blocks[args.transcript].append(b)
        #print i, len(milestones)
        #print milestones[i + 1], tx.get_seq_len(args.transcript)
        if i == len(milestones) - 2 and next_m < tx.get_seq_len(args.transcript):
            b = Block(args.transcript, -1, next_m, h.qlen, False)
            blocks[args.transcript].append(b)
            
    print '-----------------------------'
    print blocks[args.transcript]
    print '-----------------------------'

    # Determine the blocks of 'reference' transcripts
    for h in hits:
        rname = h.rname
        if h.qname == args.transcript and not h.qname == h.rname:
            blocks[rname] = []
            tx_cursors[rname] = 0
            tx_node_ids[rname] = 0
            for i in range(h.n_blocks):
                if i == 0:
                    if h.r_block_starts[0] > 0:
                        b = Block(rname, len(blocks[rname]), 0, h.r_block_starts[0], False)
                        blocks[rname].append(b)
                q_start = h.q_block_starts[i]
                q_end = q_start + h.block_sizes[i]
                # For those blocks covered by this transcript, add them
                for j in range(len(blocks[args.transcript])):
                    b = blocks[args.transcript][j]
                    diff = h.r_block_starts[i] - h.q_block_starts[i]
                    if (b.start >= q_start and b.end <= q_end):
                        new_b = Block(rname, b.id, b.start + diff, b.end + diff, True, b)
                        blocks[rname].append(new_b)
                    if b.end >= q_end:
                        break
                # Add the 'not hit' block in between the hit blocks
                if i < h.n_blocks - 1 and h.r_block_starts[i] + h.block_sizes[i] < h.r_block_starts[i + 1]:
                    b = Block(rname, len(blocks[rname]), h.r_block_starts[i] + h.block_sizes[i], h.r_block_starts[i + 1], False)
                    blocks[rname].append(b)
                # If all blocks do not cover till the end, add the remaining area as a 'not hit' block
                if i == h.n_blocks - 1 and h.block_sizes[i] + h.r_block_starts[i] < h.rlen:
                    b = Block(rname, len(blocks[rname]), h.block_sizes[i] + h.r_block_starts[i], h.rlen, False)
                    blocks[rname].append(b)
    
    for t, tx_blocks in blocks.iteritems():
        print t, tx_blocks
        if len(blocks[t]) > 0:
            dot.write('\t%s_0 [fixedsize=true, width=3, style=filled, fillcolor=%s, shape=box, label="%s:%d"] \n' % (nice_name(t), 'orange', t, tx.get_seq_len(t)))
    
    print 'n_solid_blocks ', n_solid_blocks
    
    for i in range(n_solid_blocks + 1):
        # Check whether there is a 'not hit' block before hit block i
        has_dummy = False
        all_covered = True
        
        for t, tx_blocks in blocks.iteritems():
            if tx_cursors[t] < len(tx_blocks):
                b = tx_blocks[tx_cursors[t]]
                matched_b = b.matched_block
                if not b.is_covered:
                    has_dummy = True
                    all_covered = False
                if b.id != i:
                    has_dummy = True
        
        print '----------------------'    
        print 'has_dummy ', has_dummy
        print 'all_covered ', all_covered
        print tx_cursors
        # For each 'referece transcript', draw until the current 'hit block' i
        for t, tx_blocks in blocks.iteritems():
            cursor = tx_cursors[t]
            if cursor >= len(tx_blocks):
                continue
            b = tx_blocks[cursor]
            if has_dummy:
                if b.is_covered:
                    if not all_covered:
                        dot.write('\t%s_%d [fixedsize=true, width=2, style="filled,dashed", fillcolor=%s, shape=box, label=""] \n' % (nice_name(t), tx_node_ids[t] + 1, get_color(args.transcript, t, b)))
                        dot.write('\t%s_%d -> %s_%d \n' % (nice_name(t), tx_node_ids[t], nice_name(t), tx_node_ids[t] + 1))
                        tx_node_ids[t] += 1
                    if b.id == i:
                        dot.write('\t%s_%d [fixedsize=true, width=2, style="filled", fillcolor=%s, shape=box, label="%d: %d~%d"] \n' % (nice_name(t), tx_node_ids[t] + 1, get_color(args.transcript, t, b), b.end - b.start, b.start, b.end))
                        dot.write('\t%s_%d -> %s_%d \n' % (nice_name(t), tx_node_ids[t], nice_name(t), tx_node_ids[t] + 1))
                        tx_cursors[t] += 1
                        tx_node_ids[t] += 1
                    else:
                        dot.write('\t%s_%d [fixedsize=true, width=2, style="filled,dashed", fillcolor=%s, shape=box, label=""] \n' % (nice_name(t), tx_node_ids[t] + 1, get_color(args.transcript, t, b)))
                        dot.write('\t%s_%d -> %s_%d \n' % (nice_name(t), tx_node_ids[t], nice_name(t), tx_node_ids[t] + 1))
                        tx_node_ids[t] += 1
                else:
                    dot.write('\t%s_%d [fixedsize=true, width=2, style="filled", fillcolor=%s, shape=box, label="%d: %d~%d"] \n' % (nice_name(t), tx_node_ids[t] + 1, get_color(args.transcript, t, b), b.end - b.start, b.start, b.end))
                    dot.write('\t%s_%d -> %s_%d \n' % (nice_name(t), tx_node_ids[t], nice_name(t), tx_node_ids[t] + 1))
                    tx_cursors[t] += 1
                    tx_node_ids[t] += 1
                    if tx_cursors[t] < len(tx_blocks):
                        b = tx_blocks[tx_cursors[t]]
                        if b.is_covered and b.id == i:
                            dot.write('\t%s_%d [fixedsize=true, width=2, style="filled", fillcolor=%s, shape=box, label="%d: %d~%d"] \n' % (nice_name(t), tx_node_ids[t] + 1, get_color(args.transcript, t, b), b.end - b.start, b.start, b.end))
                            dot.write('\t%s_%d -> %s_%d \n' % (nice_name(t), tx_node_ids[t], nice_name(t), tx_node_ids[t] + 1))
                            tx_cursors[t] += 1
                            tx_node_ids[t] += 1
                        else:
                            dot.write('\t%s_%d [fixedsize=true, width=2, style="filled,dashed", fillcolor=%s, shape=box, label=""] \n' % (nice_name(t), tx_node_ids[t] + 1, get_color(args.transcript, t, b)))
                            dot.write('\t%s_%d -> %s_%d \n' % (nice_name(t), tx_node_ids[t], nice_name(t), tx_node_ids[t] + 1))
                            tx_node_ids[t] += 1
            else:
                dot.write('\t%s_%d [fixedsize=true, width=2, style="filled", fillcolor=%s, shape=box, label="%d: %d~%d"] \n' % (nice_name(t), tx_node_ids[t] + 1, get_color(args.transcript, t, b), b.end - b.start, b.start, b.end))
                dot.write('\t%s_%d -> %s_%d \n' % (nice_name(t), tx_node_ids[t], nice_name(t), tx_node_ids[t] + 1))
                tx_cursors[t] += 1
                tx_node_ids[t] += 1
        
    dot.write('} \n')
    dot.close()
    
def gen_cov(args):
    tx = eva.FastaFile(args.transcript)
    hits = eva.read_blat_hits(args.psl, 'ref')
    cov = open(args.transcript + '.cov', 'w')
    for tx_name, seq in tx.seqs.iteritems():
        n_reads = 0
        coverage = 0.0
        ref_len = len(seq)
        if tx_name in hits:
            n_reads = len(hits[tx_name])
            coverage = n_reads / ref_len
        cov.write('%s\t%d%.2f' % (tx_name, n_reads, coverage))
    print 'Done. Check file %s.cov' % args.transcript 
    cov.close()
    
def zoom_region(args):
    region_s = args.start
    region_e = args.end
    tx = args.transcript
    
    cmd = 'grep %s %s' % (tx, args.psl)
    print cmd
    lines = runInShell(cmd)
    hit_lines = lines.split('\n')
    hits = eva.read_psl_hits(hit_lines, 'ref')
    
    n_reads = 0
    n_single = 0
    n_in = 0
    n_pairs = 0
    n_pairs_in_out = 0
    n_junction = 0
    n_junction_pairs = 0
    for h in hits[tx]:
        check_mate = False
        n_reads += 1
        if h.rstart >= region_s and h.rend <= region_e:
            n_in += 1
            check_mate = True
        if h.rstart <= region_s and h.rend >= region_s:
            n_junction += 1
        if h.rstart <= region_e and h.rend >= region_e:
            n_junction += 1 
        if check_mate:
            mate_id = get_mate_id(h.qname)
            mate_hit = eva.find_hit(hits, mate_id, tx)
            if mate_hit is None:
                n_single += 1
            else:
                if mate_hit.rend < region_s or mate_hit.rstart > region_e:
                    n_pairs_in_out += 1
                if mate_hit.rstart >= region_s and mate_hit.rend <= region_e:
                    n_pairs += 1
                if mate_hit.rstart <= region_s and mate_hit.rend >= region_s:
                    n_junction_pairs += 1
                if mate_hit.rstart <= region_e and mate_hit.rend >= region_e:
                    n_junction_pairs += 1
    print 'All hit reads:            %d ' % n_reads
    print 'All ingle reads:          %d ' % n_single
    print 'Reads in the region:      %d ' % n_in
    print 'Reads in the junction:    %d ' % n_junction
    print 'Pairs in the region:      %d ' % n_pairs
    print 'Pairs in the junction:    %d ' % n_junction_pairs
    print 'One in region, one out:   %d ' % n_pairs_in_out
    
def pair_regions(args):
    tx_name = args.transcript
    summary, tx_hits = zoom_tx(tx_name, args.tx, args.psl, 'read')
    region_hits = []
    n_pairs = 0
    n_reads_1 = 0
    n_reads_2 = 0
    for h in tx_hits:
        if h.rstart > args.start_1 and h.rstart < args.end_1:
            n_reads_1 += 1
        if h.rstart > args.start_2 and h.rstart < args.end_2:
            n_reads_2 += 1
        
        mate_h = None
        mate_id = get_mate_id(h.qname)
        for m in tx_hits:
            if m.qname == mate_id:
                mate_h = m
                break
        if not mate_h is None:
            if h.rstart > args.start_1 and h.rstart < args.end_1:
                if mate_h.rstart > args.start_2 and mate_h.rstart < args.end_2:
                    region_hits.append(h)
                    region_hits.append(mate_h)
                    n_pairs += 1
    align_str = get_pair_str(args.tx, args.reads, region_hits)
    fn = '%s.%d_%d.%d_%d.pairs.hits' % (tx_name, args.start_1, args.end_1, args.start_2, args.end_2)
    hits_f = open(fn, 'w')
    hits_f.write(' ' * args.start_1)
    hits_f.write('|<' + '-' * (args.end_1 - args.start_1 - 4) + '>|')
    hits_f.write(' ' * (args.start_2 - args.end_1))
    hits_f.write('|<' + '-' * (args.end_2 - args.start_2 - 4) + '>|')
    hits_f.write(align_str)
    hits_f.close()
    
    print 'Hits in region [%d, %d]: %d' % (args.start_1, args.end_1, n_reads_1)
    print 'Hits in region [%d, %d]: %d' % (args.start_2, args.end_2, n_reads_2)
    print 'Pairs:   %d' % n_pairs
    print 'Check file %s' % fn
    
def check_solid(args):
    solid = open(args.solid, 'r')
    for line in solid:
        line = line.strip()
        cmd = 'grep %s %s' % (line, args.transcript)
        result = runInShell(cmd)
        if not result is None and not result == '':
            print '========================== %s ==========================' % line
            r = result.split('\n')
            for read in r:
                print read.strip()
                
def oracle(args):
    anno_tx = FastaFile(args.transcripts)
    oracle_set = FastaFile()
    for tx_name, seq in anno_tx.seqs.iteritems():
        print 'Checking reads on %s: length %d ...' % (tx_name, len(seq))
        lines = runInShell('grep ' + tx_name + ' ' + args.psl)
        lines = lines.strip()
        if lines == '':
            print '===== %s is not expressed =====' % tx_name
            continue
        
        tx_len = len(seq)
        summary, tx_hits = zoom_tx(tx_name, args.transcripts, args.psl, 'ref')
        bin_hit = [0 for _ in range(tx_len)]
        # tx_hits.sort(key=lambda x: int(x.qname), reverse=False)
        tx_hits_valid = []
        # Filter low quality hits
        for i in range(len(tx_hits)):
            h = tx_hits[i]
            if h.rstart <= 1:
                if h.rend - h.rstart - h.n_match <= 1:
                    tx_hits_valid.append(h)
                continue
            if h.rend >= h.qlen - 1:
                if h.rend - h.rstart - h.n_match <= 1:
                    tx_hits_valid.append(h)
                continue
            # If more than 1 mismatches, ignore
            if h.n_match >= h.qlen - 1:
                tx_hits_valid.append(h)
        # Set the bases hit
        for i in range(len(tx_hits_valid) - 1):
            h = tx_hits_valid[i]
            for j in range(h.rstart, h.rend):
                bin_hit[j] = 1
        # Count the bases hit
        hit_bases = 0
        for i in range(len(bin_hit)):
            b = bin_hit[i]
            if b == 0:
                print 'Not covered base: %d' % i
            hit_bases += b
        print '%s Hit bases: %d/%d' % (tx_name, hit_bases, tx_len)
        if hit_bases >= tx_len - 10:
            oracle_set.seqs[tx_name] = seq
        else:
            print '===== %s is not expressed =====' % tx_name
    oracle_set.save_to_disk('oracle.fa')
    print 'Check file oracle.fa'
            
def bedgraph(args):
    cmd = 'tail -n +6 %s' % args.psl
    lines = runInShell(cmd)
    hit_lines = lines.split('\n')
    for i in range(len(hit_lines) / 10000 + 1):
        start_index = i*10000
        end_index = (start_index + 10000)
        if (start_index + 10000) > len(hit_lines):
            end_index = len(hit_lines) 
        hits = eva.read_psl_hits(hit_lines[start_index:end_index], 'ref')
        for chr, chr_hits in hits.iteritems():
            for h in chr_hits:
                for b in range(h.n_blocks):
                    print '%s\t%d\t%d\t1' % (chr, h.r_block_starts[b], h.r_block_starts[b] + h.block_sizes[b])
#        break

def read_id_list(fn):
    ids = []
    with open(fn) as f:
        for line in f:
            line = line.strip()
            fs = line.split('\t')
            ids.append(fs[0])
    return ids

def intersect(args):
    ids_1 = read_id_list(args.list_1)
    ids_2 = read_id_list(args.list_2)
    common = list(set(ids_1) & set(ids_2))
    n_1 = n_2 = n_both = 0
    with open(args.list_1 + '.only', 'w') as only_1:
        for id in ids_1:
            if not id in common:
                only_1.write(id + '\n')
                n_1 += 1
    with open(args.list_2 + '.only', 'w') as only_2:
        for id in ids_2:
            if not id in common:
                only_2.write(id + '\n')
                n_2 += 1
    with open(args.list_1 + '.both', 'w') as both:
        for id in common:
            both.write(id + '\n')
    n_both = len(common)
    print 'Check file %s.only: %d' % (args.list_1, n_1)
    print 'Check file %s.only: %d' % (args.list_2, n_2)
    print 'Check file %s.both: %d' % (args.list_1, n_both)
                
def main():
    parser = ArgumentParser()
    subparsers = parser.add_subparsers(help='sub command help')
    parser_read_to_ref = subparsers.add_parser('rtx', help='align all reads to a transcript and visualize text alignments')
    parser_read_to_ref.set_defaults(func=read_to_ref)
    parser_read_to_ref.add_argument('transcript', help='annotated transcripts ensembl ID')
    parser_read_to_ref.add_argument('-r', '--reads', required=False, default=READ_SRR097897, metavar='FILE', help='reads file', dest='reads')
    parser_read_to_ref.add_argument('-f', '--tx', required=False, default=REF_SRR097897, metavar='FILE', help='annotated transcripts file', dest='tx')
    parser_read_to_ref.add_argument('-p', '--psl', required=False, default=READ_REF_PSL_SRR097897, metavar='FILE', help='reads to annotated transcripts', dest='psl')

    parser_pair_regions = subparsers.add_parser('pr', help='display paired reads of two regions')
    parser_pair_regions.set_defaults(func=pair_regions)
    parser_pair_regions.add_argument('transcript', help='annotated transcripts ensembl ID')
    parser_pair_regions.add_argument('-r', '--reads', required=False, default=READ_SRR097897, metavar='FILE', help='reads file', dest='reads')
    parser_pair_regions.add_argument('-f', '--tx', required=False, default=REF_SRR097897, metavar='FILE', help='annotated transcripts file', dest='tx')
    parser_pair_regions.add_argument('-p', '--psl', required=False, default=READ_REF_PSL_SRR097897, metavar='FILE', help='reads to annotated transcripts', dest='psl')
    parser_pair_regions.add_argument('-s1', '--start_1', required=True, type=int, help='starting point of the region 1', dest='start_1')
    parser_pair_regions.add_argument('-e1', '--end_1', required=True, type=int, help='ending point of the region 1', dest='end_1')
    parser_pair_regions.add_argument('-s2', '--start_2', required=True, type=int, help='starting point of the region 2', dest='start_2')
    parser_pair_regions.add_argument('-e2', '--end_2', required=True, type=int, help='ending point of the region 2', dest='end_2')
    
    parser_pair_to_ref = subparsers.add_parser('ptx', help='align all reads to a transcript and visualize text alignments in pairs')
    parser_pair_to_ref.set_defaults(func=pair_to_ref)
    parser_pair_to_ref.add_argument('transcript', help='annotated transcripts ensembl ID')
    parser_pair_to_ref.add_argument('-r', '--reads', required=False, default=READ_SRR097897, metavar='FILE', help='reads file', dest='reads')
    parser_pair_to_ref.add_argument('-f', '--tx', required=False, default=REF_SRR097897, metavar='FILE', help='annotated transcripts file', dest='tx')
    parser_pair_to_ref.add_argument('-p', '--psl', required=False, default=READ_REF_PSL_SRR097897, metavar='FILE', help='reads to annotated transcripts', dest='psl')
    
    parser_one_read_to_ref = subparsers.add_parser('draw', help='draw transcripts splicing patterns')
    parser_one_read_to_ref.set_defaults(func=draw_dot)
    parser_one_read_to_ref.add_argument('transcript', help='annotated transcripts ensembl ID')
    parser_one_read_to_ref.add_argument('-f', '--tx', required=False, default=REF_SRR097897, help='annotated transcripts file', metavar='FILE', dest='tx')
    parser_one_read_to_ref.add_argument('-p', '--psl', required=False, default=REF_TO_REF_SRR097897, help='ref-to-ref psl file', metavar='FILE', dest='psl')
    
    parser_ctg_to_ref = subparsers.add_parser('ctx', help='align all contigs to a transcript and visualize text alignments')
    parser_ctg_to_ref.set_defaults(func=ctg_to_ref)
    parser_ctg_to_ref.add_argument('contigs', nargs='+', help='contigs file')
    parser_ctg_to_ref.add_argument('-t', '--transcript', required=False, default=None, help='check hits only for one transcript', dest='transcript')
    parser_ctg_to_ref.add_argument('-f', '--tx', required=False, default=REF_SRR097897, metavar='FILE', help='annotated transcripts file', dest='tx')

    parset_cov = subparsers.add_parser('cov', help='calculate the coverage information of some rna-seq library against annotated transcripts')
    parset_cov.set_defaults(func=gen_cov)
    parset_cov.add_argument('transcript', help='annotated transcripts')
    parset_cov.add_argument('-p', '--psl', required=False, default=READ_REF_PSL_SRR027876, metavar='FILE', help='psl file', dest='psl')
    
    parset_oracle = subparsers.add_parser('oracle', help='check the transcripts expressed')
    parset_oracle.set_defaults(func=oracle)
    parset_oracle.add_argument('transcripts', help='annotated transcripts')
    parset_oracle.add_argument('psl', help='read-to-transcript PSL file')    
    
    parset_zoom = subparsers.add_parser('zoom', help='Check the hits within some region of a transcript')
    parset_zoom.set_defaults(func=zoom_region)
    parset_zoom.add_argument('transcript', help='the transcript to zoom in')
    parset_zoom.add_argument('-p', '--psl', required=True, metavar='FILE', help='psl file', dest='psl')
    parset_zoom.add_argument('-s', '--start', required=True, type=int, help='starting point of the region', dest='start')
    parset_zoom.add_argument('-e', '--end', required=True, type=int, help='ending point of the region', dest='end')
    
    parset_check = subparsers.add_parser('check', help='Check whether some reads on a contig are solid reads')
    parset_check.set_defaults(func=check_solid)
    parset_check.add_argument('transcript', help='the reads-to-transcript file to check')
    parset_check.add_argument('solid', help='solid reads file, ids only')
	
    parset_dot = subparsers.add_parser('dot', help='Draw all dot files of the annotated transcripts')
    parset_dot.set_defaults(func=draw_all_dots)
    
    parset_rtg = subparsers.add_parser('rtg', help='Categorize all read-to-genome alignments')
    parset_rtg.add_argument('tx_psl', help='PSL of transcript-to-genome alignment')
    parset_rtg.add_argument('read_psl', help='PSL of read-to-genome alignment')
    parset_rtg.set_defaults(func=visu_all_reads_to_ref)
    
    parset_coverage = subparsers.add_parser('coverage', help='Prepare bedgraph file for PSL hits')
    parset_coverage.add_argument('psl', help='PSL file')
    parset_coverage.set_defaults(func=bedgraph)
    
    parset_set = subparsers.add_parser('set', help='Intersect two files with ids')
    parset_set.add_argument('list_1', help='File with list of ids')
    parset_set.add_argument('list_2', help='File with list of ids')
    parset_set.set_defaults(func=intersect)    
    
    args = parser.parse_args()
    args.func(args)

if __name__ == '__main__':
    sys.exit(main())

