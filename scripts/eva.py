from __future__ import division
import sys, os, pysam
from argparse import ArgumentParser
import collections

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
		self.n50_aligned_on_all = 0
		self.n50_raw = 0
		self.n50_optimal = 0
		self.base_coverage = 0.0

	def set_n_bases(self, n):
		self.n_bases = n

	def report(self):
		print '==================================================================='
		print 'Evaluation result of ' + self.contig_fn
		print '\tAssembled base:               ' + str(self.n_bases)
		print '\tBases aligned:                ' + str(self.n_aligned_bases)
		print '\t# of contigs:                 ' + str(self.n_contigs)
		print '\t# of Full length:             ' + str(self.n_tx_full_length)
		print '\t# of one-on-one:              ' + str(self.n_tx_one_on_one)
		print '\t# of 70% covered:             ' + str(self.n_tx_covered_70)
		print '\tCovered by one contig:        ' + str(self.n_tx_one_covered)
#		print '\t# of fragmented:              ' + str(self.n_fragmented)
		print '\t# of Ctgs not aligned:        ' + str(self.n_not_aligned)
		print '\tBases not aligned:            ' + str(self.n_bases_not_aligned)
		print '\t# of Ctgs not reached:        ' + str(self.n_not_reached)
		print '\tBases not reached:            ' + str(self.n_bases_not_reached)
		print '\tOptimal N50:                  ' + str(self.n50_optimal)
		print '\tRaw N50:                      ' + str(self.n50_raw)
		print '\tAligned N50:                  ' + str(self.n50_aligned)
		print '\tAligned N50 based on all:     ' + str(self.n50_aligned_on_all)
		print '\tBase coverage:                %.2f%%' % (self.base_coverage * 100)
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
				line = line[1:]
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
		
	def cal_coverage(self, sam):
		coverage_file = open(sam + '.cov', 'w')
		sam_file = pysam.Samfile(sam, "rb")
		sorted_rnames = sorted(self.seqs.iterkeys())
		
		coverage_file.write('Transcript\tLength\t# of reads\tCoverage\n')
		for rname in sorted_rnames:
			n_reads = 0
			for a in sam_file.fetch(rname):
				n_reads += 1
			coverage_file.write('%s\t%d\t%d\t%.2f\n' % (rname, len(self.seqs[rname]), n_reads, n_reads / len(self.seqs[rname])))
		coverage_file.close()

class BlastHit(object):
	def __init__(self, qname = None, rname = None, alen = 0, pos = 0, astart = 0, aend = 0, tid = 0):
		self.qname = qname
		self.n_mismatch = 0
		self.n_match = 0
		self.n_rep_match = 0
		self.n_count = 0
		self.n_query_gap = 0
		self.n_query_gap_bases = 0
		self.n_ref_gap = 0
		self.n_ref_gap_bases = 0
		self.strand = '+'
		self.qlen = 0
		self.rname = rname
		self.alen = alen
		self.qstart = 0
		self.qend = 0
		self.rlen = 0
		self.rstart = 0
		self.rend = 0
		self.n_blocks = 0
		self.block_sizes = []
		self.q_block_starts = []
		self.r_block_starts = []
		self.similarity = 0.0
		self.identity = 0.0
		self.n_bad_bases = 0
	
	def set_similarity(self):
		self.alen = 0
		for l in self.block_sizes:
			self.alen += int(l)
		if not self.alen == 0:
			self.identity = self.n_match / self.alen
		if not (self.rend - self.rstart == 0):
			self.similarity = self.n_match / (self.rend - self.rstart)
		self.n_bad_bases = self.n_mismatch + self.n_ref_gap_bases

def eva_hits(args, ref, contigs, summary, hits, r_hits, aligned_lengths):
	file_full_length = open(os.path.join(args.out_dir, 'full_length.txt'), 'w')
	file_one_on_one = open(os.path.join(args.out_dir, 'one_on_one.txt'), 'w')
	file_covered_70 = open(os.path.join(args.out_dir, 'covered_70.txt'), 'w')
	file_one_covered = open(os.path.join(args.out_dir, 'one_covered.txt'), 'w')
	similarity = float(args.similarity)

	for tx_name, tx_seq in ref.seqs.iteritems():
		if not tx_name in r_hits:
			summary.n_not_reached += 1
			summary.n_bases_not_reached += len(tx_seq)
		else:
			r_hits[tx_name].sort(key=lambda x: x.alen, reverse=True)

	n_obtained_bases = 0
	for tx_name, tx_seq in ref.seqs.iteritems():
		is_set = False
		if tx_name in r_hits:
			#print tx_name
			for a in r_hits[tx_name]:
				if a.similarity >= similarity and (a.rend - a.rstart) >= len(tx_seq) * 0.9 and a.alen >= contigs.get_seq_len(a.qname) * 0.9 and a.n_bad_bases <= 10:
					summary.n_tx_one_on_one += 1
					file_one_on_one.write(tx_name + '\t' + str(a.similarity) + '\t' + str(a.alen) + '\n')
					is_set = True
					break
			if not is_set:
				for a in r_hits[tx_name]:
					if a.similarity >= similarity and (a.rend - a.rstart) >= len(tx_seq) * 0.9 and a.n_bad_bases <= 10:
						summary.n_tx_full_length += 1
						is_set = True
						file_full_length.write(tx_name + '\t' + str(a.similarity) + '\t' + str(a.alen) + '\n')
						break
			if not is_set:
				for a in r_hits[tx_name]:
					if a.similarity >= similarity and (a.rend - a.rstart) >= len(tx_seq) * 0.7 and a.n_bad_bases <= 10:
						summary.n_tx_covered_70 += 1
						is_set = True
						file_covered_70.write(tx_name + '\t' + str(a.similarity) + '\t' + str(a.alen) + '\n')
						break

			if not is_set: 
				for a in r_hits[tx_name]:
					if (a.rend - a.rstart) >= len(tx_seq) * 0.9:
						summary.n_tx_one_covered += 1
						file_one_covered.write(tx_name + '\t' + str(a.similarity) + '\t' + str(a.alen) + '\n')
						
			seq_len = len(tx_seq)
			binary_covered = [0 for x in range(seq_len)]
			for a in r_hits[tx_name]:
				for i in range(a.n_blocks):
					start = a.r_block_starts[i]
					end = a.r_block_starts[i] + a.block_sizes[i]
					if end > len(tx_seq):
						end = len(tx_seq)
					for i in range(start, end):
						binary_covered[i] = 1
			for i in binary_covered:
				n_obtained_bases += i
				
	summary.n_bases = contigs.n_bases
	summary.n_contigs = contigs.n_seqs
	summary.n_fragmented = summary.n_contigs - summary.n_tx_full_length - summary.n_tx_one_on_one - summary.n_tx_covered_70 - summary.n_not_aligned
	summary.n50_aligned = get_n50(aligned_lengths)
	summary.n50_aligned_on_all = get_n50(aligned_lengths, summary.n_bases)
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

def analyze(args, ref, contigs, aligns, hits):
	report = open(os.path.join(args.out_dir, 'report.txt'), 'w')
	sorted_rnames = sorted(ref.seqs.iterkeys())
	report.write('Transcript\tLength\tContigs\tHits\tCovered\tLargest Covered\t# of reads\tCoverage\n')
	for rname in sorted_rnames:
		ref_len = len(ref.seqs[rname])
		report.write(rname + '\t' + str(ref_len) + '\t')
		if rname in hits:
			hit = hits[rname]
			report.write(str(len(hit)) + '\t') # Number of contigs
			hit.sort(key=lambda x: x.pos)
			longest_hit = None
			covered = 0
			most_end = 0
			for h in hit:
				report.write('[' + h.qname + ':' + str(h.pos) + ',' + str(h.pos + h.alen) + ']; ')
				if longest_hit is None:
					longest_hit = h
					covered = h.alen
					most_end = h.pos + h.alen - 1
				else:
					if longest_hit.alen < h.alen:
						longest_hit = h
					if h.pos + h.alen > most_end:
						if h.pos > most_end:
							covered += h.alen
						else:
							covered += h.pos + h.alen - most_end - 1
						most_end = h.pos + h.alen
			report.write('\t')
			covered_perc = covered / ref_len
			report.write('%.2f\t' % (covered_perc))
			largest_cover_perc = longest_hit.alen / ref_len
			report.write('%.2f\t' % (largest_cover_perc))
		else:
			report.write('0\tNONE\t0\t0\t')
		report.write('\n')
	report.close()

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

def read_blat_hits(blat_fn):
	blat = open(blat_fn, 'r')
	hits = {}
	qname = ''
	reading_hits = False
	line_no = 0
	for line in blat:
		line = line.strip()
		line_no += 1
		if True:
			if '-----' in line:
				reading_hits = True 
				continue
			if not reading_hits:
				continue
			f = line.split('\t')
			hit = BlastHit(f[9])
			hit.n_match = int(f[0])
			hit.n_mismatch = int(f[1])
			hit.n_rep_match = int(f[2])
			hit.n_count = int(f[3])
			hit.n_query_gap = int(f[4])
			hit.n_query_gap_bases = int(f[5])
			hit.n_ref_gap = int(f[6])
			hit.n_ref_gap_bases = int(f[7])
			hit.strand = f[8]
			hit.qlen = int(f[10])
			hit.qstart = int(f[11])
			hit.qend = int(f[12])
			hit.rname = f[13]
			hit.rlen = int(f[14])
			hit.rstart = int(f[15])
			hit.rend = int(f[16])
			hit.n_blocks = int(f[17])
			f[18] = f[18][0:-1]
			hit.block_sizes = f[18].split(',')
			hit.block_sizes = [int(s) for s in hit.block_sizes]
			f[19] = f[19][0:-1]
			hit.q_block_starts = f[19].split(',')
			hit.q_block_starts = [int(s) for s in hit.q_block_starts]
			f[20] = f[20][0:-1]
			hit.r_block_starts = f[20].split(',')
			hit.r_block_starts = [int(s) for s in hit.r_block_starts]
			hit.set_similarity()
			if not hit.qname in hits:
				hits[hit.qname] = []
			hits[hit.qname].append(hit)
		
	blat.close()
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
			# print a.qname, rid, a.pos, a.alen, a.astart, a.aend
	# For the hits on the same transcript, get the longest alignment length
	for rid, a_list in hits.iteritems():
		longest = 0
		a_longest = None
		for a in a_list:
			if a.alen > longest:
				longest = a.alen
				a_longest = a
		if longest > 0:
			aligned_lengths.append(longest)
			# print a_longest.qname, a_longest.rname, a_longest.pos, a_longest.alen, a_longest.astart, a_longest.aend
	eva_hits(args, ref, contigs, aligns, summary, hits, aligned_lengths)
	analyze(args, ref, contigs, aligns, hits)
	if not args.sam is None:
		ref.cal_coverage(args.sam)

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

def eva_blat(args):
	ref = FastaFile(args.ref)
	contigs = FastaFile(args.contigs)
	summarize_fa(args.ref)
	summarize_fa(args.contigs)
	ref.set_n50(get_n50(ref.lengths, ref.n_bases))
	contigs.set_n50(get_n50(contigs.lengths, contigs.n_bases))

	summary = ResultSummary(args.contigs)
	hits = read_blat_hits(args.psl)
	r_hits = {}
	aligned_lengths = {}	
	# r_hits: ref_name->all hits on this reference
	# aligned_lengths: ref_name->longest hit length
	for qname in contigs.get_seqs():
		if not qname in hits:
			summary.n_not_aligned += 1
			summary.n_bases_not_aligned += contigs.get_seq_len(qname)
	for qname, qhits in hits.iteritems():
		longest = 0
		for hit in qhits:
			rname = hit.rname
			for l in hit.block_sizes:
				summary.n_aligned_bases += l
				if int(l) > longest:
					longest = int(l)
			if rname in aligned_lengths:
				if aligned_lengths[rname] < longest:
					aligned_lengths[rname] = longest
			else:
				aligned_lengths[rname] = longest
				
			if not rname in r_hits:
				r_hits[rname] = []
			r_hits[rname].append(hit)

	lengths = []
	for rname, length in aligned_lengths.iteritems():
		lengths.append(length)
	
	eva_hits(args, ref, contigs, summary, hits, r_hits, lengths)

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

    parser_bwa = subparsers.add_parser('blat', help='evaluate the performance by aligning in blat')
    parser_bwa.set_defaults(func=eva_blat)
    parser_bwa.add_argument('-t', required=True, help='reference transcript file', dest='ref')
    parser_bwa.add_argument('-c', required=True, help='transcripts reported to be evaluated', dest='contigs')
    parser_bwa.add_argument('-p', required=True, help='psl file', dest='psl')
    parser_bwa.add_argument('-o', required=True, help='output folder', dest='out_dir')
    parser_bwa.add_argument('-s', required=True, help='similarity', dest='similarity')

    parser_bla = subparsers.add_parser('blast', help='evaluate the performance by alinging by Blastn')
    parser_bla.set_defaults(func=eva_blastn)
    parser_bla.add_argument('-t', required=True, help='reference transcript file', dest='ref')
    parser_bla.add_argument('-c', required=True, help='transcripts reported to be evaluated', dest='contigs')
    parser_bla.add_argument('-b', required=True, help='blastn file', dest='blastn')
    parser_bla.add_argument('-o', required=True, help='output folder', dest='out_dir')
    parser_bla.add_argument('-s', required=False, help='reads sam file', dest='sam')

    parser_ana = subparsers.add_parser('report', help='evaluate the performance by alinging by Blastn')
    parser_ana.set_defaults(func=analyze)
    parser_ana.add_argument('-t', required=True, help='reference transcript file', dest='ref')
    parser_ana.add_argument('-c', required=True, help='transcripts reported to be evaluated', dest='contigs')
    parser_ana.add_argument('-b', required=True, help='blastn file', dest='blastn')
    parser_ana.add_argument('-o', required=True, help='output folder', dest='out_dir')

    parser_cmp = subparsers.add_parser('hasdup', help='evaluate the performance by alinging by Blastn')
    parser_cmp.set_defaults(func=check_dup)
    parser_cmp.add_argument('-f', required=True, help='file', dest='input')

    args = parser.parse_args()
    args.func(args)

if __name__ == '__main__':
    sys.exit(main())
