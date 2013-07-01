from eva import *
from zoom import *
from merge import *
from rna import *
from argparse import ArgumentParser
import operator, random

class Vertex(object):
    def __init__(self, tpl_id='-1', length=0, n_reads=0, seq=''):
        self.tpl_id = tpl_id
        self.n_reads = n_reads
        self.length = length
        self.seq = seq
        
    def __repr__(self):
        s = 'Template [%s, %d]\n' % (self.tpl_id, self.length)
        s += '%d reads: %s\n' % (self.n_reads, self.seq)
        return s

class Graph(object):
    def __init__(self, junction_file='', contigs_file=''):
        self.vertexes = []
        self.junctions = []
        self.contigs = FastaFile(contigs_file)
        if not junction_file == '':
            self.read_from_file(junction_file)
    
    def read_from_file(self, junction_file):
        if self.contigs.seqs:
            seqs = self.contigs.seqs
        with open(junction_file) as jf:
            line_no = 0
            for line in jf:
                line_no += 1
                line = line.strip()
                if line_no <= 1 or line == '':
                    continue
                f = line.split('\t')
                if f[-1] == '-1':
                    tpl_id, length = read_tpl_info(f[0])
                    seq = ''
                    if seqs:
                        seq = seqs[tpl_id]
                    new_tpl = Vertex(tpl_id, length, int(f[3]), seq)
                    self.vertexes.append(new_tpl)
                else:
                    main_id, main_len = read_tpl_info(f[0])
                    branch_id, branch_len = read_tpl_info(f[1])
                    new_j = Junction(main_id, branch_id, main_len, branch_len, int(f[2]), int(f[3]), int(f[4]))
                    self.junctions.append(new_j)
                    
    def break_by_exons(self):
        pass
                    
    def __repr__(self):
        s = '================ Templates ================\n'
        for c in self.vertexes:
            s += '> %s \n' % str(c)
        s += '================ Junctions =================\n'
        for j in self.junctions:
            s += '> %s' % str(j)
        return s

def main():
    parser = ArgumentParser()
    subparsers = parser.add_subparsers(help='sub command help')
    
    parser_splice = subparsers.add_parser('filter', help='Filter not qualified junctions')
    parser_splice.set_defaults(func=split_exons)
    parser_splice.add_argument('psl', help='transcript/contigs-to-genome PSL file')
    
    args = parser.parse_args()
    args.func(args)
                
if __name__ == '__main__':
    g = Graph('/home/carl/Projects/peta_dev/SRR097897_part/paired.junctions', '/home/carl/Projects/peta_dev/SRR097897_part/paired.fa')
    print g
#    sys.exit(main())