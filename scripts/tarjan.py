
class Vertex(object):
    def __init__(self, id=0):
        self.id = id
        self.index = 0
        self.lowlink = 0
        self.outs = []
        self.ins = []
    
    def __repr__(self):
        s = '==== Vertex %d ====\n' % self.id
        s += 'Ins: \n'
        for e in self.ins:
            s += '\t'
            s += str(e)
            s += '\n'
        s += 'Outs: \n'
        for e in self.outs:
            s += '\t'
            s += str(e)
            s += '\n'
        return s

class Edge(object):
    def __init__(self, id=0, left=None, right=None):
        self.id = id
        self.left = left
        self.right = right
        self.left.outs.append(self)
        self.right.ins.append(self)
        
    def __repr__(self):
        return 'Edge: %d -> %d' % (self.left.id, self.right.id)
        
class Graph(object):
    def __init__(self):
        self.vertexes = []
        self.edges = []
        self.index = 0
        self.stack = []
    
    def add_vertex(self, v):
        self.vertexes.append(v)
        
    def add_edge(self, e):
        self.edges.append(e)

    def draw(self):
        with open('graph.dot', 'w') as dot:
            dot.write('digraph G {\n')
            dot.write('graph [rankdir=LR];\n')
            for v in self.vertexes:
                dot.write("%d [label=\"%d\" shape=box]; \n" % (v.id, v.id))
            for e in self.edges:
                dot.write('%d -> %d [label=\"\"]; \n' % (e.left.id, e.right.id))
            dot.write('}\n')

def tarjan(g):
    for v in g.vertexes:
        if v.index == 0:
            connect(v, g)

def connect(v, g):
    v.index = g.index
    v.lowlink = g.index
    g.index += 1
    g.stack.append(v)
    
    print '================ Vertex %d ================' % v.id
    for o in v.outs:
        w = o.right
        if w.index == 0:
            connect(w, g)
            v.lowlink = min(v.lowlink, w.lowlink)
        elif w in g.stack:
            print g.stack
            v.lowlink = min(v.lowlink, w.index)
    
    this_scc = []        
    if (v.lowlink == v.index):
        for i in range(len(g.stack)):
            w = g.stack.pop()
            if w == v:
                this_scc.append(w)
                break
            else:
                this_scc.append(w)
#    if len(this_scc):
#        print '=== SCC root: Vertex %d ===' % v.id
#        print this_scc

g = Graph()
for i in range(0, 8):
    g.add_vertex(Vertex(i))
    
e1 = Edge(1, g.vertexes[0], g.vertexes[1])
e2 = Edge(2, g.vertexes[1], g.vertexes[7])
e3 = Edge(3, g.vertexes[7], g.vertexes[1])
e4 = Edge(4, g.vertexes[1], g.vertexes[2])
e5 = Edge(5, g.vertexes[2], g.vertexes[3])
e6 = Edge(6, g.vertexes[2], g.vertexes[5])
e7 = Edge(7, g.vertexes[3], g.vertexes[4])
e8 = Edge(8, g.vertexes[4], g.vertexes[5])
e9 = Edge(9, g.vertexes[5], g.vertexes[6])
e10 = Edge(10, g.vertexes[3], g.vertexes[6])

#for i in range(0, 10):
#    g.add_vertex(Vertex(i))
#
#e1 = Edge(1, g.vertexes[0], g.vertexes[2])
#e2 = Edge(2, g.vertexes[1], g.vertexes[3])
#e3 = Edge(3, g.vertexes[0], g.vertexes[1])
#e4 = Edge(4, g.vertexes[2], g.vertexes[3])
#e5 = Edge(5, g.vertexes[3], g.vertexes[0])
#
#e6 = Edge(6, g.vertexes[2], g.vertexes[4])
#e7 = Edge(7, g.vertexes[3], g.vertexes[5])
#e8 = Edge(8, g.vertexes[4], g.vertexes[5])
#e9 = Edge(9, g.vertexes[4], g.vertexes[6])
#e10 = Edge(10, g.vertexes[6], g.vertexes[5])
#
#e11 = Edge(11, g.vertexes[6], g.vertexes[7])
#e12 = Edge(12, g.vertexes[6], g.vertexes[8])
#e13 = Edge(13, g.vertexes[6], g.vertexes[9])
#e14 = Edge(14, g.vertexes[7], g.vertexes[9])
#e15 = Edge(15, g.vertexes[8], g.vertexes[9])
#e16 = Edge(16, g.vertexes[9], g.vertexes[8])
#e17 = Edge(17, g.vertexes[5], g.vertexes[4])
#
g.add_edge(e1)
g.add_edge(e2)
g.add_edge(e3)
g.add_edge(e4)
g.add_edge(e5)
g.add_edge(e6)
g.add_edge(e7)
g.add_edge(e8)
g.add_edge(e9)
g.add_edge(e10)
#g.add_edge(e11)
#g.add_edge(e12)
#g.add_edge(e13)
#g.add_edge(e14)
#g.add_edge(e15)
#g.add_edge(e16)
#g.add_edge(e17)

g.draw()

tarjan(g)