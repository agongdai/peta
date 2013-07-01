from eva import *
from zoom import *
from merge import *
from rna import *
from argparse import ArgumentParser
import operator, random

class Vertex(object):
    def __init__(self, id=0, len=0, n_reads=0):
        self.id = id
        self.n_reads = n_reads
        self.len = len
        self.seq = ''

class Graph(object):
    def __init__(self):
        self.vertexes = []
        self.junctions = []
    
    def read_from_file(self, junction_file):
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
                    id, len = read_tpl_info(f[0])
                    new_tpl = Vertex(id, len, int(f[3]))
                    self.vertexes.append(new_tpl)
                else:
                    new_j = Junction(read_tpl(f[0]), read_tpl(f[1]), int(f[2]), int(f[3]), int(f[4]))
                    self.junctions.append(new_j)

# Split string like '[1, 233]'        
def read_tpl_info(tpl_str):
    tpl_str = tpl_str[1:-1]
    f = tpl_str.split(',')
    id = int(f[0].strip())
    length = int(f[1].strip())
    return id, length

def main():
    parser = ArgumentParser()
    subparsers = parser.add_subparsers(help='sub command help')
    
    parser_splice = subparsers.add_parser('filter', help='Filter not qualified junctions')
    parser_splice.set_defaults(func=split_exons)
    parser_splice.add_argument('psl', help='transcript/contigs-to-genome PSL file')
    
    args = parser.parse_args()
    args.func(args)
                
if __name__ == '__main__':
    sys.exit(main())