import operator

seq = open('/home/carl/Projects/peta/rnaseq/Spombe/SRR097897/SRR097897.fa')
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
    for i in range(len(line) - k + 1):
        s = line[i:i+k]
        #print s
        try:
            counter[s] += 1
        except:
            counter[s] = 1
#    if read_no == 1:
#        break
    if read_no % 1000000 == 0:
        print '%d reads kmer counted...' % read_no

#i = 0
#for s, count in counter.iteritems():
#    i += 1
#    print '%d-%s: %d' % (i, s, count)

seq.close()
print '%d reads kmer counted.' % read_no
print 'Saving to %s' % kmer_fn
sorted_counter = sorted(counter.iterkeys())
for s in sorted_counter:
    number = counter[s]
    kmer.write('>%d\n' % number)
    kmer.write('%s\n' % s)