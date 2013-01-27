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
    summary, tx_hits = zoom_tx(tx_name, ref, READ_REF_PSL, 'read')
    print summary

def zoom_tx(tx_name, ref, blat_psl, ctg_or_read):
    lines = runInShell('grep ' + tx_name + ' ' + blat_psl)
    hit_lines = lines.split('\n')
    tx_hits = []
    if len(hit_lines) <= 0:
        return tx_hits
    hits = eva.read_psl_hits(hit_lines)
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
        hit_str += '    ' + h.qname + ':' + str(len(read_seq)) + ' [' + str(h.rstart) + 'M' + str(h.n_match) + ',' + str(h.strand) + str(h.n_mismatch) + ']\n'
        algin_str += hit_str
        read_no += 1
    return algin_str

def draw_dot(args):
    tx = FastaFile(args.tx)
    dot = open(args.transcript + '.dot', 'w')
    cmd = 'grep %s %s' % (args.transcript, args.tx)
    hit_lines = runInShell(cmd)
    milestones = []
    
    lines = runInShell('grep ' + args.transcript + ' ' + args.psl)
    hit_lines = lines.split('\n')
    raw_hits = eva.read_psl_hits(hit_lines)
    hits = raw_hits[args.transcript]
    
    # Every starting/ending point on the query transcript is a 'milestone'
    m = []
    for h in hits:
        if h.qname == args.transcript and not h.qname == h.rname:
            for i in range(h.n_blocks):
                if not h.q_block_starts[i] in m:
                    milestones.append((h.q_block_starts[i], 'start', i))
                    m.append(h.q_block_starts[i])
                if not (h.q_block_starts[i] + h.block_sizes[i]) in m:
                    milestones.append((h.q_block_starts[i] + h.block_sizes[i], 'end', i))
                    m.append(h.q_block_starts[i] + h.block_sizes[i])
                     
    blocks = {}
    milestones.sort()
    blocks[args.transcript] = []
    dot.write('digraph g { \n\trankdir = LR \n')
    # Determine all blocks of the query transcript
    for i in range(len(milestones) - 1):
        (next_m, next_start_or_end, next_index) = milestones[i + 1]
        (m, start_or_end, index) = milestones[i]
        is_covered = True
        if i == 0 and m > 0:
            blocks[args.transcript].append((0, m, False))
        if start_or_end == 'end' and next_start_or_end == 'start' and index < next_index:
            is_covered = False
        blocks[args.transcript].append((m, next_m, is_covered))
        print i, len(milestones)
        print milestones[i + 1], tx.get_seq_len(args.transcript)
        if i == len(milestones) - 2 and next_m < tx.get_seq_len(args.transcript):
            blocks[args.transcript].append((next_m, h.qlen, False))

    # Determine the blocks of 'reference' transcripts
    n_solid_blocks = 0
    for h in hits:
        rname = h.rname
        if h.qname == args.transcript and not h.qname == h.rname:
            blocks[rname] = []
            for i in range(h.n_blocks):
                if i == 0:
                    if h.r_block_starts[0] > 0:
                        blocks[rname].append((0, h.r_block_starts[0], False))
                q_start = h.q_block_starts[i]
                q_end = q_start + h.block_sizes[i]
                # For those blocks covered by this transcript, add them
                for (m, next_m, is_covered) in blocks[args.transcript]:
                    diff = h.r_block_starts[i] - h.q_block_starts[i]
                    if (m >= q_start and next_m <= q_end):
                        blocks[rname].append((m + diff, next_m + diff, True))
                        n_solid_blocks += 1
                    if next_m > q_end:
                        break
                # Add the 'not hit' block in between the hit blocks
                if i < h.n_blocks - 1:
                    blocks[rname].append((h.r_block_starts[i] + h.block_sizes[i], h.r_block_starts[i + 1], False))
                # If all blocks do not cover till the end, add the remaining area as a 'not hit' block
                if i == h.n_blocks - 1 and h.block_sizes[i] + h.r_block_starts[i] < h.rlen:
                    blocks[rname].append((h.block_sizes[i] + h.r_block_starts[i], h.rlen, False))
    
    for t, tx_blocks in blocks.iteritems():
        if len(blocks[t]) > 0:
            (m, next_m, is_covered) = blocks[t][j]
            dot.write('\t%s_0 [fixedsize=true, width=3, style=filled, fillcolor=%s, shape=box, label="%s:%d"] \n' % (t, 'orange', t, tx.get_seq_len(t)))
    
    has_dummy = False
    for t, tx_blocks in blocks.iteritems():
        if len(blocks[t]) > 0:
            (m, next_m, is_covered) = blocks[t][j]
            if is_covered:
                has_dummy = True
                break

    n_has_drawn = 1
    for t, tx_blocks in blocks.iteritems():
        color = 'white'
        if t == args.transcript:
            color = 'yellow'
        if len(blocks[t]) > 0:
            (m, next_m, is_covered) = blocks[t][0]
            if has_dummy:
                if not is_covered:
                    dot.write('\t%s_%d [fixedsize=true, width=2, style=filled, fillcolor=%s, shape=box, label="%s:%d"] \n' % (t, n_has_drawn, color, t, tx.get_seq_len(t)))
                    dot.write('\t%s_%d -> %s_%d \n' % (t, n_has_drawn - 1, t, n_has_drawn))
    
    for t, tx_blocks in blocks.iteritems():
        if t == args.transcript:
            continue        
        for i in range(len(blocks[t])):
            (m, next_m, is_covered) = blocks[t][i]
            color = 'white'
            if not is_covered:
                color = 'gray'
            if i == 0:
                dot.write('\t%s [fixedsize=true, width=3, style=filled, fillcolor=%s, shape=box, label="%s:%d"] \n' % (t, 'orange', t, tx.get_seq_len(t)))
                dot.write('\t%s_%d [fixedsize=true, width=2, style=filled, fillcolor=%s, shape=box, label="%d: %d~%d"] \n' % (t, i, color, next_m - m, m, next_m))
                dot.write('\t%s -> %s_%d \n' % (t, t, 0))
            else:
                dot.write('\t%s_%d [fixedsize=true, width=2, style=filled, fillcolor=%s, shape=box, label="%d: %d~%d"] \n' % (t, i, color, next_m - m, m, next_m))
                dot.write('\t%s_%d -> %s_%d \n' % (t, i - 1, t, i))
                
    for i in range(len(blocks[args.transcript])):
        (m, next_m, is_covered) = blocks[args.transcript][i]
        color = 'yellow'
        if not is_covered:
            color = 'gray'
        if i == 0:
            dot.write('\t%s [fixedsize=true, width=3, style=filled, fillcolor=%s, shape=box, label="%s:%d"] \n' % (args.transcript, 'orange', args.transcript, tx.get_seq_len(args.transcript)))
            dot.write('\t%s_%d [fixedsize=true, width=2, style=filled, fillcolor=%s, shape=box, label="%d: %d~%d"] \n' % (args.transcript, i, color, next_m - m, m, next_m))
            dot.write('\t%s -> %s_%d \n' % (args.transcript, args.transcript, 0))
        else:
            dot.write('\t%s_%d [fixedsize=true, width=2, style=filled, fillcolor=%s, shape=box, label="%d: %d~%d"] \n' % (args.transcript, i, color, next_m - m, m, next_m))
            dot.write('\t%s_%d -> %s_%d \n' % (args.transcript, i - 1, args.transcript, i))

    dot.write('} \n')
    dot.close()

def main():
    parser = ArgumentParser()
    subparsers = parser.add_subparsers(help='sub command help')
    parser_read_to_ref = subparsers.add_parser('rtx', help='align all reads to a transcript and visualize text alignments')
    parser_read_to_ref.set_defaults(func=read_to_ref)
    parser_read_to_ref.add_argument('transcript', help='annotated transcripts ensembl ID')
    parser_read_to_ref.add_argument('-r', '--reads', required=False, default=READ, metavar='FILE', help='reads file', dest='reads')
    parser_read_to_ref.add_argument('-f', '--tx', required=False, default=REF, metavar='FILE', help='annotated transcripts file', dest='tx')
    
    parser_one_read_to_ref = subparsers.add_parser('draw', help='draw transcripts splicing patterns')
    parser_one_read_to_ref.set_defaults(func=draw_dot)
    parser_one_read_to_ref.add_argument('transcript', help='annotated transcripts ensembl ID')
    parser_one_read_to_ref.add_argument('-f', '--tx', required=False, default=REF, help='annotated transcripts file', metavar='FILE', dest='tx')
    parser_one_read_to_ref.add_argument('-p', '--psl', required=False, default=REF_TO_REF, help='ref-to-ref psl file', metavar='FILE', dest='psl')
    
    parser_ctg_to_ref = subparsers.add_parser('ctr', help='align all contigs to a transcript and visualize text alignments')
    parser_ctg_to_ref.set_defaults(func=ctg_to_ref)
    parser_ctg_to_ref.add_argument('contigs', help='contigs file')
    parser_ctg_to_ref.add_argument('-f', '--tx', required=False, default=REF, metavar='FILE', help='annotated transcripts file', dest='tx')

    args = parser.parse_args()
    args.func(args)

if __name__ == '__main__':
    sys.exit(main())

