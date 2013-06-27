from zoom import *

def psl2gene(psl, asm):
    cmd = r"""tail -n +6 %s | awk -F "\t" 'BEGIN {{ OFS="\t"; }} {{ if ($5 == 0) print $10, $14, $9, $16, $17, $16, $17, $18, $21, $19;}}'""" % (psl)
    print cmd
    lines = runInShell(cmd)
    gene_lines = lines.split('\n')
    real_lines = []
    with open(psl + '.gene', 'w') as genes:
        genes.write('%s.name\t' % asm)
        genes.write('%s.chrom\t' % asm)
        genes.write('%s.strand\t' % asm)
        genes.write('%s.txStart\t' % asm)
        genes.write('%s.txEnd\t' % asm)
        genes.write('%s.cdsStart\t' % asm)
        genes.write('%s.cdsEnd\t' % asm)
        genes.write('%s.exonCount\t' % asm)
        genes.write('%s.exonStarts\t' % asm)
        genes.write('%s.exonEnds\t' % asm)
        genes.write('%s.kgXref.geneSymbol\t' % asm)
        genes.write('%s.kgXref.mRNA\n')
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
    
def bed2gene(bed, asm):
    with open('%s.gene' % bed, 'w') as genes:
        genes.write('%s.name\t' % asm)
        genes.write('%s.chrom\t' % asm)
        genes.write('%s.strand\t' % asm)
        genes.write('%s.txStart\t' % asm)
        genes.write('%s.txEnd\t' % asm)
        genes.write('%s.cdsStart\t' % asm)
        genes.write('%s.cdsEnd\t' % asm)
        genes.write('%s.exonCount\t' % asm)
        genes.write('%s.exonStarts\t' % asm)
        genes.write('%s.exonEnds\t' % asm)
        genes.write('%s.kgXref.geneSymbol\t' % asm)
        genes.write('%s.kgXref.mRNA\n')
        with open(bed) as b:
            for line in b:
                line = line.strip()
                if line == '':
                    continue
                fields = line.split('\t')
                genes.write(fields[3] + '\t')
                genes.write(fields[0] + '\t')
                genes.write(fields[5] + '\t')
                genes.write(fields[1] + '\t')
                genes.write(fields[2] + '\t')
                genes.write(fields[1] + '\t')
                genes.write(fields[2] + '\t')
                genes.write(fields[9] + '\t')
                #genes.write(fields[11] + '\t')
                s = fields[11].split(',')
                for i in range(int(fields[9])):
                    genes.write('%d,' % (int(s[i]) + int(fields[1])))
                genes.write('\t')
                sz = fields[10].split(',')
                ends = ''
                for i in range(int(fields[9])):
                    size = int(sz[i])
                    start = int(s[i])
                    ends += '%d,' % (start + size + int(fields[1]))
                genes.write(ends + '\t')
                genes.write(fields[3] + '\t')
                genes.write(fields[3] + '\n')
    print 'Check %s.gene' % bed
    
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
psl2gene('/home/ariyaratnep/shaojiang/peta_copies/peta_dev/spombe_630/idba-60.fa.psl', 'hg19')
#bed2gene('/home/ariyaratnep/shaojiang/peta/rnaseq/hg19/ensembl.genes.bed')
#count_not_align('/home/carl/Projects/peta/rnaseq/Spombe/630.tx.rev.psl')
