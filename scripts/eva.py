from __future__ import division
import sys, os, pysam
from argparse import ArgumentParser

class ResultSummary(object):
	def __init__(self, contig_fn):
		self.contig_fn = contig_fn
		self.n_bases = 0
		self.n_aligned_bases = 0
		self.n_contigs = 0
		self.n_tx_full_length = 0
		self.n_tx_one_on_one = 0
		self.n_tx_one_covered = 0
		self.n_tx_covered_70 = 0
		self.n_not_aligned = 0
		self.n_fragmented = 0
		self.n_bases_not_aligned = 0
		self.n_not_reached = 0
		self.n_bases_not_reached = 0
		self.n50_aligned = 0
		self.n50_raw = 0
		self.n50_optimal = 0
		self.base_coverage = 0.0

	def set_n_bases(self, n):
		self.n_bases = n

	def report(self):
		print '==================================================================='
		print 'Evaluation result of ' + self.contig_fn
		print '\tAssembled base: \t\t' + str(self.n_bases)
		print '\tBases aligned: \t\t\t' + str(self.n_aligned_bases)
		print '\t# of contigs: \t\t\t' + str(self.n_contigs)
		print '\t# of Full length: \t\t' + str(self.n_tx_full_length)
		print '\t# of one-on-one: \t\t' + str(self.n_tx_one_on_one)
		print '\t# of 70% covered: \t\t' + str(self.n_tx_covered_70)
		print '\tCovered by one contig: \t\t' + str(self.n_tx_one_covered)
		print '\t# of fragmented: \t\t' + str(self.n_fragmented)
		print '\t# of Ctgs not aligned: \t\t' + str(self.n_not_aligned)
		print '\tBases not aligned: \t\t' + str(self.n_bases_not_aligned)
		print '\t# of Ctgs not reached: \t\t' + str(self.n_not_reached)
		print '\tBases not reached: \t\t' + str(self.n_bases_not_reached)
		print '\tOptimal N50: \t\t\t' + str(self.n50_optimal)
		print '\tRaw N50: \t\t\t' + str(self.n50_raw)
		print '\tAligned N50: \t\t\t' + str(self.n50_aligned)
		print '\tBase coverage: \t\t\t%.2f%%' % (self.base_coverage * 100)
		print '==================================================================='

class FastaFile(object):
	def __init__(self, fa=''):
		self.filename = fa
		self.seqs = {}
		self.lengths = []
		self.n_seqs = 0
		self.ave_len = 0
		self.min_len = -1
		self.max_len = 0
		self.n50 = 0
		self.n_bases = 0
		if not fa == '':
			self.read_seqs()
			self.set_lengths()

	def read_seqs(self):
		fa = open(self.filename, 'r')
		name = ''
		s = ''
		for line in fa:
			line = line.strip()
			if '>' in line:
				if not name == '' and not s == '':
					self.seqs[name] = s
					self.lengths.append(len(s))
				line = line[1:-1]
				fields = line.split()
				name = fields[0]
				s = ''
			else:
				s += line
		if not name == '' and not s == '':
			self.seqs[name] = s
			self.lengths.append(len(s))
		fa.close()

	def set_lengths(self):
		self.n_bases = 0
		self.n_seqs = len(self.lengths)
		for l in self.lengths:
			self.n_bases += l
			if l > self.max_len:
				self.max_len = l
			if l < self.min_len or self.min_len == -1:
				self.min_len = l
		if self.n_seqs > 0:
			self.ave_len = self.n_bases / self.n_seqs

	def set_n50(self, n50):
		self.n50 = n50

	def get_seqs(self):
		return self.seqs

	def get_seq_len(self, seq_name):
		return len(self.seqs[seq_name])

	def print_summary(self):
		print '==================================================================='
		print 'Summary of ' + self.filename
		print '\tTotal base: \t\t' + str(self.n_bases)
		print '\t# of transcripts: \t' + str(self.n_seqs)
		print '\tLongest length: \t' + str(self.max_len)
		print '\tShortest length: \t' + str(self.min_len)
		print '\tAverage length: \t%.2f' % self.ave_len
		print '\tN50 value: \t\t' + str(self.n50)
		print '==================================================================='

class BlastHit(object):
	def __init__(self, qname = None, rname = None, alen = 0, pos = 0, astart = 0, aend = 0, tid = 0):
		self.qname = qname
		self.rname = rname
		self.alen = alen
		self.pos = pos
		self.astart = astart
		self.aend = aend
		self.tid = 0

def eva_hits(args, ref, contigs, aligns, summary, hits, aligned_lengths):
	file_full_length = open(os.path.join(args.out_dir, 'full_length.txt'), 'w')
	file_one_on_one = open(os.path.join(args.out_dir, 'one_on_one.txt'), 'w')
	file_covered_70 = open(os.path.join(args.out_dir, 'covered_70.txt'), 'w')
	file_one_covered = open(os.path.join(args.out_dir, 'one_covered.txt'), 'w')

	for tx_name, tx_seq in ref.seqs.iteritems():
		if not tx_name in hits:
			summary.n_not_reached += 1
			summary.n_bases_not_reached += len(tx_seq)
		else:
			hits[tx_name].sort(key=lambda x: x.pos, reverse=True)

	n_obtained_bases = 0
	for tx_name, tx_seq in ref.seqs.iteritems():
		is_set = False
		if tx_name in hits:
			#print tx_name
			for a in hits[tx_name]:
				if a.alen >= len(tx_seq) * 0.9 and a.alen >= contigs.get_seq_len(a.qname) * 0.9:
					summary.n_tx_one_on_one += 1
					file_one_on_one.write(tx_name + '\n')
					is_set = True
					break
			if not is_set:
				for a in hits[tx_name]:
					if a.alen >= len(tx_seq) * 0.9:
						summary.n_tx_full_length += 1
						is_set = True
						file_full_length.write(tx_name + '\n')
						break
			if not is_set:
				for a in hits[tx_name]:
					if a.alen >= len(tx_seq) * 0.7:
						summary.n_tx_covered_70 += 1
						is_set = True
						file_covered_70.write(tx_name + '\n')
						break
			seq_len = len(tx_seq)
			binary_covered = [0 for x in range(seq_len)]
			group_hits = {}
			for a in hits[tx_name]:
				if not a.qname in group_hits:
					group_hits[a.qname] = []
				group_hits[a.qname].append(a)
				end = a.pos + a.alen
				if end > len(tx_seq):
					end = len(tx_seq)
				for i in range(a.pos, end):
					binary_covered[i] = 1
			for i in binary_covered:
				n_obtained_bases += i

			if not is_set:
				for qname, algs in group_hits.iteritems():
					binary_covered = [0 for x in range(seq_len)]
					n_base_one_contig = 0
					for a in algs:
						end = a.aend
						if a.aend > len(tx_seq):
							end = len(tx_seq)
						for i in range(a.pos, end):
							binary_covered[i] = 1
					for i in binary_covered:
						n_base_one_contig += i
					if n_base_one_contig >= seq_len * 0.9:
						summary.n_tx_one_covered += 1
						file_one_covered.write(tx_name + '\n')
						break
	
	summary.n_bases = contigs.n_bases
	summary.n_contigs = contigs.n_seqs
	summary.n_fragmented = summary.n_contigs - summary.n_tx_full_length - summary.n_tx_one_on_one - summary.n_tx_covered_70 - summary.n_not_aligned
	summary.n50_aligned = get_n50(aligned_lengths)
	summary.base_coverage = n_obtained_bases / ref.n_bases
	summary.n50_raw = contigs.n50
	summary.n50_optimal = ref.n50
	summary.report()

	print 'Check %s.'%os.path.join(args.out_dir, 'full_length.txt')
	print 'Check %s.'%os.path.join(args.out_dir, 'one_on_one.txt')
	print 'Check %s.'%os.path.join(args.out_dir, 'covered_70.txt')
	print 'Check %s.'%os.path.join(args.out_dir, 'one_covered.txt')

	file_one_covered.close()
	file_full_length.close()
	file_covered_70.close()
	file_one_on_one.close()

def read_blastn_hits(blastn_fn):
	blastn = open(blastn_fn, 'r')
	hits = []
	qname = ''
	for line in blastn:
		line = line.strip()
		if line.startswith('#'):
			if 'Query' in line:
				f = line.split()
				qname = f[2]
			if '0 hits found' in line:
				hit = BlastHit(qname)
				hits.append(hit)
			continue
		f = line.split('\t')
		pos = f[7]
		if int(f[6]) < int(f[7]):
			pos = f[6]
		hit = BlastHit(f[0], f[1], int(f[3]), int(pos), int(f[9]), int(f[10]))
		hits.append(hit)
	blastn.close()
	return hits

def eva_blastn(args):
	ref = FastaFile(args.ref)
	contigs = FastaFile(args.contigs)
	aligns = read_blastn_hits(args.blastn)

	summarize_fa(args.ref)
	summarize_fa(args.contigs)
	ref.set_n50(get_n50(ref.lengths, ref.n_bases))
	contigs.set_n50(get_n50(contigs.lengths, contigs.n_bases))

	summary = ResultSummary(args.contigs)

	hits = {}
	aligned_lengths = []
	# Set all hits
	for a in aligns: 
		if a.alen == 0:
			summary.n_not_aligned += 1
			summary.n_bases_not_aligned += contigs.get_seq_len(a.qname)
		else:
			rid = a.rname
			if not rid in hits:
				hits[rid] = []
			hits[rid].append(a)
			summary.n_aligned_bases += a.alen
			aligned_lengths.append(a.alen)
			# print a.qname, rid, a.pos, a.alen, a.astart, a.aend
	eva_hits(args, ref, contigs, aligns, summary, hits, aligned_lengths)

def get_n50(arr, total_length=0):
	arr.sort()
	arr.reverse()
	n50_len = 0
	n50 = 0
	if total_length == 0:
		for l in arr:
			total_length += l
	for l in arr:
		if n50_len * 2 >= total_length:
			n50 = l
			break
		else:
			n50_len += l
	return n50

def differ(args):
	files = args.inputs.split(',')
	out_file = open(args.output, 'w')
	genes = {}
	for f in files:
		fp = open(f, 'r')
		for line in fp:
			line = line.strip()
			if line in genes:
				genes[line] += ',' + f
			else:
				genes[line] = f
		fp.close()
	categories = {}
	for key, value in sorted(genes.iteritems(), key=lambda (k,v): (v,k)):
		out_file.write(key + ': ' + value  + '\n')
		if value in categories:
			categories[value] += 1
		else:
			categories[value] = 0
	for key, value in categories.iteritems():
		print key, value
	out_file.close()

def summarize_fa(fa_file):
	fa = FastaFile(fa_file)
	fa.set_n50(get_n50(fa.lengths, fa.n_bases))
	fa.print_summary()	

def eva_bwa(args):
	ref = FastaFile(args.ref)
	contigs = FastaFile(args.contigs)
	sam = pysam.Samfile(args.sam, "r")
	summarize_fa(args.ref)
	summarize_fa(args.contigs)
	ref.set_n50(get_n50(ref.lengths, ref.n_bases))
	contigs.set_n50(get_n50(contigs.lengths, contigs.n_bases))

	summary = ResultSummary(args.contigs)
	aligns = sam.fetch()
	hits = {}
	aligned_lengths = []
	# Set all hits
	for a in aligns:
		if a.tid == -1:
			summary.n_not_aligned += 1
			summary.n_bases_not_aligned += contigs.get_seq_len(a.qname)
		else:
			rid = sam.getrname(a.tid)
			if not rid in hits:
				hits[rid] = []
			hits[rid].append(a)
			summary.n_aligned_bases += a.alen
			aligned_lengths.append(a.alen)
			# print a.qname, rid, a.pos, a.aend
	eva_hits(args, ref, contigs, aligns, summary, hits, aligned_lengths)

def check_dup(args):
	ids = open(args.input, 'r')
	unique = []
	counter = 0
	for line in ids:
		line = line.strip()	
		counter += 1
		if line in unique:
			print 'Line ' + str(counter) + ': Duplicated: ' + line
		else:
			unique.append(line)
	print '\nDone.'
	ids.close()

def main():
    parser = ArgumentParser()
    subparsers = parser.add_subparsers(help='sub command help')
    parser_differ = subparsers.add_parser('diff', help='differ multiple files')
    parser_differ.set_defaults(func=differ)
    parser_differ.add_argument('-i', required=True, help='the files to compare, seperated by a ","', dest='inputs')
    parser_differ.add_argument('-o', required=True, help='result file', dest='output')

    parser_bwa = subparsers.add_parser('bwa', help='evaluate the performance by aligning in BWA')
    parser_bwa.set_defaults(func=eva_bwa)
    parser_bwa.add_argument('-t', required=True, help='reference transcript file', dest='ref')
    parser_bwa.add_argument('-c', required=True, help='transcripts reported to be evaluated', dest='contigs')
    parser_bwa.add_argument('-s', required=True, help='sam file', dest='sam')
    parser_bwa.add_argument('-o', required=True, help='output folder', dest='out_dir')

    parser_bla = subparsers.add_parser('blast', help='evaluate the performance by alinging by Blastn')
    parser_bla.set_defaults(func=eva_blastn)
    parser_bla.add_argument('-t', required=True, help='reference transcript file', dest='ref')
    parser_bla.add_argument('-c', required=True, help='transcripts reported to be evaluated', dest='contigs')
    parser_bla.add_argument('-b', required=True, help='blastn file', dest='blastn')
    parser_bla.add_argument('-o', required=True, help='output folder', dest='out_dir')

    parser_cmp = subparsers.add_parser('hasdup', help='evaluate the performance by alinging by Blastn')
    parser_cmp.set_defaults(func=check_dup)
    parser_cmp.add_argument('-f', required=True, help='file', dest='input')

    args = parser.parse_args()
    args.func(args)

if __name__ == '__main__':
    sys.exit(main())
