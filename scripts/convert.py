from zoom import *

def psl2gene(psl):
    cmd = r"""tail -n +6 %s | awk -F "\t" 'BEGIN {{ OFS="\t"; }} {{ print $10, $14, $9, $16, $17, $16, $17, $18, $21, $19;}}'""" % (psl)
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
            ends += ('\t' + fields[0]) * 9
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
    
#psl2gene('/home/carl/Projects/peta/rnaseq/Spombe/genome/tx.rev.genome.psl')
#psl2gene('/home/carl/Projects/peta_dev/scripts/spombe_630.genome.psl')
#psl2gene('/home/carl/Projects/peta_dev/SRR097897_out/Trinity.630.genome.psl')