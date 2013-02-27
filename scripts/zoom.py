import sys, os, pysam
from subprocess import Popen, PIPE
import eva, merge
from eva import FastaFile
from argparse import ArgumentParser

BLAT = '/home/carl/Projects/blat/blat'
BLAT_OCC_11 = '/home/carl/Projects/blat/11.ooc'

READ_REF_PSL = '/home/carl/Projects/peta/rnaseq/hg19/SRX011545/SRR027876.psl'
REF_TO_REF = '/home/carl/Projects/peta/rnaseq/hg19/genome/ref.ref.psl'
REF = '/home/carl/Projects/peta/rnaseq/hg19/genome/human.ensembl.cdna.fa'
READ = '/home/carl/Projects/peta/rnaseq/hg19/SRX011545/SRR027876.fa'
CONTIG = '/home/carl/Projects/peta/SRR027876_out/pair_contigs.fa'

READ_REF_PSL_SRR097897 = '/home/carl/Projects/peta/rnaseq/Spombe/SRR097897/SRR097897.psl'
REF_TO_REF_SRR097897 = '/home/carl/Projects/peta/rnaseq/Spombe/genome/ref.ref.psl'
REF_SRR097897 = '/home/carl/Projects/peta/rnaseq/Spombe/genome/spombe.broad.tx.fasta'
READ_SRR097897 = '/home/carl/Projects/peta/rnaseq/Spombe/SRR097897/SRR097897.fa'

READ_SRR027876 = '/home/carl/Projects/peta/rnaseq/hg19/SRX011545/SRR027876.fa'
READ_REF_PSL_SRR027876 = '/home/carl/Projects/peta/rnaseq/hg19/SRX011545/SRR027876.psl'

def runInShell(cmd):
    p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    return p.communicate()[0]

def sub_read(fasta, query):
    query_read_psl = query + '.psl'
    cmd = '%s %s %s -ooc=%s %s' % (BLAT, fasta, query, BLAT_OCC_11, query_read_psl)
    
def ctg_to_ref(args):
    contig = args.contigs
    ref = args.tx
    print 'Aligning %s to %s ' % (contig, ref)
    ctg_ref_psl = contig + '.psl'
    if args.transcript:
        hits_file_name = args.transcript + '.hits'
    else:
        hits_file_name = contig + '.hits'
    hits_file = open(hits_file_name, 'w')
    cmd = '%s %s %s -ooc=%s %s' % (BLAT, ref, contig, BLAT_OCC_11, ctg_ref_psl)
    runInShell(cmd)
    ref_fa = FastaFile(ref)
    no = 0
    open_psl = open(ctg_ref_psl, 'r')
    hit_tx = {}
    for line in open_psl:
        no += 1
        if no <= 5:
            continue
        f = line.split('\t')
        hit_tx[f[13]] = line
    open_psl.close()
    print 'Printing the text alignments to %s...' % hits_file_name
    for tx_name, tx_seq in hit_tx.iteritems():
        if args.transcript and not args.transcript == tx_name:
            continue
        no += 1
        summary, hits = zoom_tx(tx_name, ref, ctg_ref_psl, 'ctg')
        if len(hits) <= 0:
            continue
        hits_file.write(get_align_str(ref, contig, hits))
        hits_file.write('=' * 500 + '\n')
    hits_file.close()
    print 'Check file %s' % (hits_file_name)
        
def read_to_ref(args):
    tx_name = args.transcript
    ref = args.tx
    ctg_or_read = 'read'
    summary, tx_hits = zoom_tx(tx_name, ref, args.psl, 'read')
    align_str = get_align_str(ref, args.reads, tx_hits)
    hits_f = open(tx_name + '.reads.hits', 'w')
    hits_f.write(align_str)
    hits_f.close()
    print summary
    print 'Check text alignment at file %s.reads.hits' % tx_name

def zoom_tx(tx_name, ref, blat_psl, ctg_or_read):
    lines = runInShell('grep ' + tx_name + ' ' + blat_psl)
    hit_lines = lines.split('\n')
    tx_hits = []
    if len(hit_lines) <= 0:
        return tx_hits
    hits = eva.read_psl_hits(hit_lines, 'query')
    tx = FastaFile(ref)
    tx_seq = tx.seqs[tx_name]
    for qname, h in hits.iteritems():
        tx_hits.append(h[0])
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
            summary += '# of pairs:       %d\n' % n_pairs
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


# For a list of hits, get the text alignments
def get_align_str(ref, fasta, hits):
    if len(hits) <= 0:
        return ''
    tx_name = hits[0].rname
    algin_str = tx_name + '\n'
    tx = FastaFile(ref)
    reads = FastaFile(fasta)
    seq = tx.seqs[tx_name]
    read_no = 0
    for h in hits:
        if read_no % 30 == 0:
            algin_str += seq + '\n'
        read_seq = reads.seqs[h.qname].upper()
        rev_read_seq = merge.rev_comp(read_seq).upper()
        pre_end = 0
        for i in range(0, h.n_blocks):
            if pre_end > 0:
                hit_str += '-' * (h.r_block_starts[i] - pre_end)
            else:
                hit_str = ' ' * h.r_block_starts[i]
            if h.strand == '+':
                hit_str += read_seq[h.q_block_starts[i]:h.q_block_starts[i] + h.block_sizes[i]] 
            else:
                hit_str += rev_read_seq[h.q_block_starts[i]:h.q_block_starts[i] + h.block_sizes[i]] 
            pre_end = h.r_block_starts[i] + h.block_sizes[i]
        hit_str += '    ' + h.qname + ':' + str(len(read_seq)) + ' [' + str(h.rstart) + 'M' + str(h.n_match) + ',III' + str(h.n_ref_gap_bases) + str(h.strand) + str(h.n_mismatch) + ']\n'
        algin_str += hit_str
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
    cmd = 'grep %s %s' % (args.transcript, args.tx)
    hit_lines = runInShell(cmd)
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
            if int(h.qname) % 2 == 0:
                mate_id = str(int(h.qname) + 1)
            else:
                mate_id = str(int(h.qname) - 1)
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
            
def main():
    parser = ArgumentParser()
    subparsers = parser.add_subparsers(help='sub command help')
    parser_read_to_ref = subparsers.add_parser('rtx', help='align all reads to a transcript and visualize text alignments')
    parser_read_to_ref.set_defaults(func=read_to_ref)
    parser_read_to_ref.add_argument('transcript', help='annotated transcripts ensembl ID')
    parser_read_to_ref.add_argument('-r', '--reads', required=False, default=READ_SRR097897, metavar='FILE', help='reads file', dest='reads')
    parser_read_to_ref.add_argument('-f', '--tx', required=False, default=REF_SRR097897, metavar='FILE', help='annotated transcripts file', dest='tx')
    parser_read_to_ref.add_argument('-p', '--psl', required=False, default=READ_REF_PSL_SRR097897, metavar='FILE', help='reads to annotated transcripts', dest='psl')
    
    parser_one_read_to_ref = subparsers.add_parser('draw', help='draw transcripts splicing patterns')
    parser_one_read_to_ref.set_defaults(func=draw_dot)
    parser_one_read_to_ref.add_argument('transcript', help='annotated transcripts ensembl ID')
    parser_one_read_to_ref.add_argument('-f', '--tx', required=False, default=REF_SRR097897, help='annotated transcripts file', metavar='FILE', dest='tx')
    parser_one_read_to_ref.add_argument('-p', '--psl', required=False, default=REF_TO_REF_SRR097897, help='ref-to-ref psl file', metavar='FILE', dest='psl')
    
    parser_ctg_to_ref = subparsers.add_parser('ctx', help='align all contigs to a transcript and visualize text alignments')
    parser_ctg_to_ref.set_defaults(func=ctg_to_ref)
    parser_ctg_to_ref.add_argument('contigs', help='contigs file')
    parser_ctg_to_ref.add_argument('-t', '--transcript', required=False, default=None, help='check hits only for one transcript', dest='transcript')
    parser_ctg_to_ref.add_argument('-f', '--tx', required=False, default=REF_SRR097897, metavar='FILE', help='annotated transcripts file', dest='tx')

    parset_cov = subparsers.add_parser('cov', help='calculate the coverage information of some rna-seq library against annotated transcripts')
    parset_cov.set_defaults(func=gen_cov)
    parset_cov.add_argument('transcript', help='annotated transcripts')
    parset_cov.add_argument('-p', '--psl', required=False, default=READ_REF_PSL_SRR027876, metavar='FILE', help='psl file', dest='psl')
    
    parset_zoom = subparsers.add_parser('zoom', help='Check the hits within some region of a transcript')
    parset_zoom.set_defaults(func=zoom_region)
    parset_zoom.add_argument('transcript', help='the transcript to zoom in')
    parset_zoom.add_argument('-p', '--psl', required=True, metavar='FILE', help='psl file', dest='psl')
    parset_zoom.add_argument('-s', '--start', required=True, type=int, help='starting point of the region', dest='start')
    parset_zoom.add_argument('-e', '--end', required=True, type=int, help='ending point of the region', dest='end')
    
    args = parser.parse_args()
    args.func(args)

if __name__ == '__main__':
    sys.exit(main())

