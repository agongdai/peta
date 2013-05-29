from zoom import *

def psl2gene(psl):
    cmd = r"""tail -n +6 %s | awk -F "\t" 'BEGIN {{ OFS="\t"; }} {{ if ($5 == 0) print $10, $14, $9, $16, $17, $16, $17, $18, $21, $19;}}'""" % (psl)
    print cmd
    lines = runInShell(cmd)
    gene_lines = lines.split('\n')
    real_lines = []
    with open(psl + '.gene', 'w') as genes:
        genes.write('spombe.name\tspombe.chrom\tspombe.strand\tspombe.txStart\tspombe.txEnd\tspombe.cdsStart\tspombe.cdsEnd\tspombe.exonCount\tspombe.exonStarts\tspombe.exonEnds\tspombe.kgXref.geneSymbol\tspombe.kgXref.mRNA\n')
        for l in gene_lines:
            if l == '':
                continue
            fields = l.split('\t')
            sizes = fields[-1].split(',')
            starts = fields[-2].split(',')
            ends = ''
            for f in fields[:-1]:
                ends += f + '\t'
            for i in range(len(starts)):
                if starts[i] == '':
                    continue
                s = int(starts[i])
                sz = int(sizes[i])
                ends += str(s + sz) + ','
            ends = ends[:-1]
            ends += ('\t' + fields[0]) * 2
            genes.write(ends + '\n')
    print 'Check file %s.gene' % psl
    
def fasta2chrs(fasta):
    with open(fasta) as f:
        out = None
        for line in f:
            if '>' in line:
                if out:
                    out.close()
                chr = line[1:-1]
                out = open(chr, 'w')
            out.write(line)

def count_not_align(psl):
    reads_have_hits = {}
    transcripts = {}
    with open(psl) as f:
        line_no = 0
        for line in f:
            line_no += 1
            if line_no <= 5:
                continue
            fields = line.split('\t')
            transcripts[fields[13]] = 1
            reads_have_hits[fields[9]] = 1
    print 'Reads having hits on %d transcripts: %d' % (len(transcripts), len(reads_have_hits))            
    
#psl2gene('/home/carl/Projects/peta/rnaseq/Spombe/genome/tx.rev.genome.psl')
#psl2gene('/home/carl/Projects/peta_dev/scripts/spombe_630.genome.psl')
#psl2gene('/home/carl/Projects/peta/rnaseq/Spombe/genome/tx.630.oracle.genome.psl')
#psl2gene('/home/carl/Projects/peta_dev/spombe_630/paired.630.oracle.psl')
psl2gene('/home/carl/Projects/peta_dev/SRR097897_out/paired.genome.psl')

#count_not_align('/home/carl/Projects/peta/rnaseq/Spombe/630.tx.rev.psl')
