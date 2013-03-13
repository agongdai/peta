import operator

seq = open('../rnaseq/Spombe/SRR097897/SRR097897.fa')
kmer_fn = '../SRR097897_out/kmer.freq'
kmer = open(kmer_fn, 'w')
k = 25
counter = {}
read_no = 0
for l in seq:
    if '>' in l:
        continue
    read_no += 1
    line = l.upper().strip()
    for i in range(len(line) - k):
        s = line[i:i+k]
        if not s in counter:
            counter[s] = 0
        counter[s] += 1
    if read_no % 1000000 == 0:
        print '%d reads kmer counted...' % read_no

seq.close()
print '%d reads kmer counted.' % read_no
print 'Sorting the kmer frequencies...'
sorted_counter = sorted(counter.iteritems(), key=operator.itemgetter(1), reverse=True)
print 'Saving to %s' % kmer_fn
for s, number in sorted_counter.iteritems():
    kmer.write('>%d\n' % number)
    kmer.write('%s\n' % s)