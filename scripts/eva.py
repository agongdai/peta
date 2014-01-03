from __future__ import division
import sys, os, pysam
from argparse import ArgumentParser
import collections
from subprocess import Popen, PIPE

bad_bases_thre = 10
full_length_perc = 0.99
near_full_length = 0.9

def runInShell(cmd):
    p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    return p.communicate()[0]

class ResultSummary(object):
	def __init__(self, contig_fn):
		self.contig_fn = contig_fn
		self.n_bases = 0 
		self.n_aligned_bases = 0
		self.n_contigs = 0
		self.n_tx_full_length = 0
		self.n_tx_near_full_length = 0
		self.n_tx_one_on_one = 0
		self.n_tx_one_covered = 0
		self.n_tx_covered_70 = 0
		self.n_not_aligned = 0
		self.n_fragmented = 0
		self.n_bases_not_aligned = 0
		self.n_not_reached = 0
		self.n_bases_not_reached = 0
		self.n_chimaeras = 0
		self.n50_aligned = 0
		self.n50_aligned_on_all = 0
		self.n50_raw = 0
		self.n50_optimal = 0
		self.base_coverage = 0.0
		self.similarity_thre = 0
		
		self.accuracy = 0.0
		self.completeness = 0.0
		self.contiguity = 0.0
		self.chimerism = 0.0
		self.variant_resolution = 0.0
		self.threshold = 0.8

	def set_n_bases(self, n):
		self.n_bases = n

	def report(self):
		print '==================================================================='
		print 'Evaluation result of ' + self.contig_fn
		print '[similarity threshold: ' + str(self.similarity_thre) + ']'
		print '[full length threshold: ' + str(full_length_perc) + ']'
		print '[n_mismatches + n_indels < ' + str(bad_bases_thre) + ']'
		print '\tAssembled base:               ' + str(self.n_bases)
		print '\tBases aligned:                ' + str(self.n_aligned_bases)
		print '\t# of contigs:                 ' + str(self.n_contigs)
		print '\t# of Full length:             ' + str(self.n_tx_full_length + self.n_tx_one_on_one)
# 		print '\t# of near full length %.0f%%:   ' % (near_full_length * 100) + ' ' + str(self.n_tx_near_full_length)
# 		print '\t# of 70% covered:             ' + str(self.n_tx_covered_70)
# 		print '\tCovered by one contig:        ' + str(self.n_tx_one_covered)
# 		print '\t# of fragmented:              ' + str(self.n_fragmented)
		print '\t# of Ctgs not aligned:        ' + str(self.n_not_aligned)
		print '\tBases not aligned:            ' + str(self.n_bases_not_aligned)
		print '\t# of Ctgs not reached:        ' + str(self.n_not_reached)
		print '\tBases not reached:            ' + str(self.n_bases_not_reached)
		print '\tOptimal N50:                  ' + str(self.n50_optimal)
		print '\tRaw N50:                      ' + str(self.n50_raw)
		print '\tAligned N50:                  ' + str(self.n50_aligned)
# 		print '\tAligned N50 based on all:     ' + str(self.n50_aligned_on_all)
# 		print '\tBase coverage:                %.2f%%' % (self.base_coverage * 100)
		print '\tAccuracy:                     %.2f%%' % (self.accuracy * 100)
		print '\tCompleteness %.0f%%:             %.2f%%' % (self.threshold * 100, self.completeness * 100)
		print '\tContiguity %.0f%%:               %.2f%%' % (self.threshold * 100, self.contiguity * 100)
		print '\t# of chimaeras:               ' + str(self.n_chimaeras)
		print '\tChimerism:                    %.2f%%' % (self.chimerism * 100)
		print '\tVariant resolution:           %.2f%%' % (self.variant_resolution * 100)
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
			
	def save_to_disk(self, output):
		with open(output, 'w') as out:
			for tx_name, seq in self.seqs.iteritems():
				out.write('>%s length: %d\n' % (tx_name, len(seq)))
				l = 0
				for c in seq:
					out.write(c)
					l += 1
					if l % 100 == 0:
						out.write('\n')
				if not l % 100 == 0:
					out.write('\n')

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

''' Trop a sequence line to multiple lines with equal length of line_len
'''
def trop_seq(seq, line_len):
	line = ''
	c_counter = 0
	for c in seq:
		line += c
		c_counter += 1
		if c_counter % line_len == 0:
			line += '\n'
	return line

def get_tx_str(tx_name, seq, line_len):
	header = '>%s\n' % tx_name
	seq_line = trop_seq(seq, line_len)
	return header + seq_line

class BlastHit(object):
	def __init__(self, qname=None, hit_line='', rname=None, alen=0, pos=0, astart=0, aend=0, tid=0):
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
		self.mate = None
		self.hit_line = hit_line
		
	def print_ori(self):
		line = ''
		line += str(self.qname)
		return line
	
	def set_similarity(self):
		self.alen = 0
		for l in self.block_sizes:
			self.alen += int(l)
		if not self.alen == 0:
			self.similarity = self.n_match / self.alen
		if not (self.rend - self.rstart == 0):
			self.identity = self.n_match / abs(self.rend - self.rstart)
		self.n_bad_bases = self.n_mismatch + self.n_ref_gap_bases
		
	def set_mate(self, mate):
		self.mate = mate
		
	def __repr__(self):
		repr = '%sQuery [%s, %d]: %d->%d\t [' % (self.strand, self.qname, self.qlen, self.qstart, self.qend)
		for s in self.block_sizes:
			repr += str(s) + ','
		repr += '] ['
		for s in self.q_block_starts:
			repr += str(s) + ','
		repr += ']\t Ref [%s, %d]: %d->%d\t [' % (self.rname, self.rlen, self.rstart, self.rend)
		for s in self.r_block_starts:
			repr += str(s) + ','
		repr += ']\t'
		repr += 'Match: %d; Mismatch: %d; Indels: %d/%d\n' % (self.n_match, self.n_mismatch, self.n_query_gap_bases, self.n_ref_gap_bases)
		return repr
	
def find_hit(hits, query, tx):
	for h in hits[tx]:
		if h.qname == query:
			return h
	return None

def cal_base_cover(annotated_tx, hits):
	base_cover = {}
	for tx_name, tx_seq in annotated_tx.iteritems():
		binary_covered = [0 for _ in range(len(tx_seq))]
		if tx_name in hits:
			tx_hits = hits[tx_name]
			for a in tx_hits:
				for i in range(a.n_blocks):
					start = a.r_block_starts[i]
					end = a.r_block_starts[i] + a.block_sizes[i]
					if end > len(tx_seq):
						end = len(tx_seq)
					for i in range(start, end):
						binary_covered[i] = 1
		base_cover[tx_name] = binary_covered
	return base_cover

def pileup_coverage(hits):
	if len(hits) <= 0: return []
	h = hits[0]
	tx_len = h.rlen
	base_cover = [0 for _ in range(tx_len)]
	hits.sort(key=lambda x:x.rstart)
	for h in hits:
		for i in range(h.n_blocks):
			s = h.r_block_starts[i]
			e = s + h.block_sizes[i]
			for j in range(s, e):
				base_cover[j] = base_cover[j] + 1
	return base_cover

def cal_accuracy(annotated_tx, hits):
	n_mis = 0.0
	n_correct = 0.0
	for tx_name, tx_hits in hits.iteritems():
		tx_hits.sort(key=lambda x: x.n_match, reverse=True)
		if len(tx_hits) > 0:
			h = tx_hits[0]
			n_correct += h.n_match
			n_mis += h.n_mismatch
			#print tx_name
			#print h.hit_line
	covered = n_correct / (n_correct + n_mis)
	return covered

def cal_completeness(annotated_tx, hits, threshold=0.99):
	base_cover = cal_base_cover(annotated_tx, hits)
	covered = 0.0
	for tx_name, binary_covered in base_cover.iteritems():
		n_covered = 0
		for i in binary_covered:
			n_covered += i
		tx_cover = n_covered / len(annotated_tx[tx_name])
		if tx_cover >= threshold:
			covered += 1
	covered /= len(annotated_tx)
	return covered

def cal_contiguity(annotated_tx, hits, threshold=0.99):
	covered = 0.0
	for tx_name, tx_seq in annotated_tx.iteritems():
		if tx_name in hits:
			max_block = 0
			tx_hits = hits[tx_name]
			for a in tx_hits:
				for i in range(a.n_blocks):
					if a.block_sizes[i] > max_block:
						max_block = a.block_sizes[i]
			if max_block >= len(annotated_tx[tx_name]) * threshold:
				covered += 1
	covered /= len(annotated_tx)
	return covered

def cal_continuous_n50(annotated_tx, hits):
	n50 = 0
	max_blocks = []
	for tx_name, tx_seq in annotated_tx.iteritems():
		if tx_name in hits:
			max = 0
			tx_hits = hits[tx_name]
			for a in tx_hits:
				for i in range(a.n_blocks):
					if a.block_sizes[i] > max:
						max = a.block_sizes[i]
			max_blocks.append(max)
	n50 = get_n50(max_blocks)
	return n50

'''
	Two regions: [s1, e1] and [s2, e2]
	Get their shared region size
'''
def share_region(s1, e1, s2, e2):
	smaller_s = s1
	first_end = e1
	larger_s = s2
	second_end = e2
	if s1 > s2:
		smaller_s = s2
		first_end = e2
		larger_s = s1
		second_end = e1
	if larger_s >= first_end:
		return 0
	if second_end <= first_end:
		return second_end - larger_s
	else:
		return first_end - smaller_s

def cal_chimerism(hits):
	n_chimaeras = 0
	for contig, ctg_hits in hits.iteritems():
		tx = {}
		for h in ctg_hits:
			if h.n_match > 50:
				bydot = h.rname.split('.')
				gene = bydot[0]
				if not gene in tx:
					tx[gene] = []
				tx[gene].append((h.qstart, h.qend))
		stop = False
		if len(tx) > 1:
			for gene, regions in tx.iteritems():
				if stop: break
				for gene2, regions2 in tx.iteritems():
					if stop: break
					if gene == gene2:
						continue
					for (start, end) in regions:
						if stop: break
						for (start2, end2) in regions2:
							if stop: break
							if share_region(start, end, start2, end2) < 50:
								#print contig
								n_chimaeras += 1
								stop = True
	return n_chimaeras
		
def eva_hits(args, ref, contigs, summary, hits, r_hits, aligned_lengths):
	prefix = os.path.join(args.out_dir, os.path.basename(contigs.filename))
	file_full = open(prefix + '.full', 'w')
	file_full_length = open(prefix + '_full_length.txt', 'w')
	file_one_on_one = open(prefix + '_one_on_one.txt', 'w')
	file_near_full_length = open(prefix + '_near_full_length.txt', 'w')
	file_covered_70 = open(prefix + '_covered_70.txt', 'w')
	file_one_covered = open(prefix + '_one_covered.txt', 'w')
	file_not_aligned = open(prefix + '_not_aligned.txt', 'w')
	file_partial = open(prefix + '_partial.txt', 'w')
	similarity = float(args.similarity)

	for tx_name, tx_seq in ref.seqs.iteritems():
		if not tx_name in r_hits:
			summary.n_not_reached += 1
			summary.n_bases_not_reached += len(tx_seq)
			file_not_aligned.write('%s\t%d\n' % (tx_name, len(tx_seq)))
		else:
			r_hits[tx_name].sort(key=lambda x: x.alen, reverse=True)

	n_obtained_bases = 0
	for tx_name, tx_seq in ref.seqs.iteritems():
		is_set = False
		if tx_name in r_hits:
			for a in r_hits[tx_name]:
				if a.similarity >= similarity and abs(a.rend - a.rstart) >= len(tx_seq) * full_length_perc and a.alen >= contigs.get_seq_len(a.qname) * 0.9 and (a.n_query_gap_bases + a.n_bad_bases <= bad_bases_thre or a.n_blocks == 1):
					summary.n_tx_one_on_one += 1
					file_one_on_one.write(tx_name + '\t' + str(a.similarity) + '\t' + str(a.alen) + '\n')
					file_full.write(tx_name + '\t' + str(a.similarity) + '\t' + str(a.alen) + '\n')
					is_set = True
					break
			if not is_set:
				for a in r_hits[tx_name]:
					# print a.qname, a.rname, a.similarity, abs(a.rend - a.rstart), len(tx_seq), (a.n_query_gap_bases + a.n_bad_bases)
					if a.similarity >= similarity and abs(a.rend - a.rstart) >= len(tx_seq) * full_length_perc and (a.n_query_gap_bases + a.n_bad_bases <= bad_bases_thre or a.n_blocks == 1):
						summary.n_tx_full_length += 1
						is_set = True
						file_full_length.write(tx_name + '\t' + str(a.similarity) + '\t' + str(a.alen) + '\n')
						file_full.write(tx_name + '\t' + str(a.similarity) + '\t' + str(a.alen) + '\n')
						break
			if not is_set:
				for a in r_hits[tx_name]:
					if a.similarity >= similarity and abs(a.rend - a.rstart) >= len(tx_seq) * near_full_length and (a.n_bad_bases <= bad_bases_thre or a.n_blocks == 1):
						summary.n_tx_near_full_length += 1
						is_set = True
						file_near_full_length.write(tx_name + '\t' + str(a.similarity) + '\t' + str(a.alen) + '\n')
						break
			if not is_set:
				for a in r_hits[tx_name]:
					if a.similarity >= similarity and abs(a.rend - a.rstart) >= len(tx_seq) * 0.7 and (a.n_bad_bases <= bad_bases_thre or a.n_blocks == 1):
						summary.n_tx_covered_70 += 1
						is_set = True
						file_covered_70.write(tx_name + '\t' + str(a.similarity) + '\t' + str(a.alen) + '\n')
						break

			if not is_set: 
				for a in r_hits[tx_name]:
					if abs(a.rend - a.rstart) >= len(tx_seq) * full_length_perc:
						summary.n_tx_one_covered += 1
						file_one_covered.write(tx_name + '\t' + str(a.similarity) + '\t' + str(a.alen) + '\n')
						is_set = True
						break
						
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
			tmp = n_obtained_bases
			for i in binary_covered:
				n_obtained_bases += i
			if (len(tx_seq) - (n_obtained_bases - tmp)) >= 50:
				file_partial.write('======= %s: %d \t %d =======\n' % (tx_name, len(tx_seq), (len(tx_seq) - (n_obtained_bases - tmp))))
				for i in range(len(binary_covered)):
					if binary_covered[i] == 0:
						file_partial.write('%d\n' % i)

	file_full.close()
	file_one_covered.close()
	file_full_length.close()
	file_near_full_length.close()
	file_covered_70.close()
	file_one_on_one.close()
	file_not_aligned.close()
				
	summary.n_bases = contigs.n_bases
	summary.n_contigs = contigs.n_seqs
	summary.n_fragmented = summary.n_contigs - summary.n_tx_full_length - summary.n_tx_one_on_one - summary.n_tx_covered_70 - summary.n_not_aligned
	summary.n50_aligned = cal_continuous_n50(ref.seqs, r_hits)
	summary.base_coverage = n_obtained_bases / ref.n_bases
	summary.n50_raw = contigs.n50
	summary.n50_optimal = ref.n50
	summary.accuracy = cal_accuracy(ref.seqs, r_hits)
 	summary.completeness = cal_completeness(ref.seqs, r_hits, summary.threshold)
 	summary.contiguity = cal_contiguity(ref.seqs, r_hits, summary.threshold)
 	summary.n_chimaeras = cal_chimerism(hits)
 	summary.chimerism = summary.n_chimaeras / len(contigs.seqs)
	summary.report()

	print 'Check %s' % os.path.join(args.out_dir, os.path.basename(contigs.filename) + '.full')
	print 'Check %s' % os.path.join(args.out_dir, os.path.basename(contigs.filename) + '_full_length.txt')
	print 'Check %s' % os.path.join(args.out_dir, os.path.basename(contigs.filename) + '_one_on_one.txt')
	print 'Check %s' % os.path.join(args.out_dir, os.path.basename(contigs.filename) + '_near_full_length.txt')
	print 'Check %s' % os.path.join(args.out_dir, os.path.basename(contigs.filename) + '_covered_70.txt')
	print 'Check %s' % os.path.join(args.out_dir, os.path.basename(contigs.filename) + '_one_covered.txt')
	print 'Check %s' % os.path.join(args.out_dir, os.path.basename(contigs.filename) + '_not_aligned.txt')
	print 'Check %s' % os.path.join(args.out_dir, os.path.basename(contigs.filename) + '_partial.txt')

def rm_self(hits):
	for key, h in hits.iteritems():
		clean_hits = []
		for one_hit in h:
			if not one_hit.qname == one_hit.rname:
				clean_hits.append(one_hit)
			else:
				if one_hit.n_blocks > 1:
					clean_hits.append(one_hit)
				else:
					if not one_hit.block_sizes[0] == one_hit.qlen or not one_hit.q_block_starts[0] == 0:
						clean_hits.append(one_hit)
					if not one_hit.block_sizes[0] == one_hit.rlen or not one_hit.r_block_starts[0] == 0:
						clean_hits.append(one_hit)	
		hits[key] = clean_hits
	return hits

def stat(args):
	report = open(os.path.join(args.out_dir, 'report.txt'), 'w')
	hits = read_blat_hits(args.ctg2ref, 'ref')
	ref = FastaFile(args.ref)
	sorted_rnames = sorted(ref.seqs.iterkeys())
	ref2ref_hits = read_blat_hits(args.ref2ref)
	ref2ref_hits = rm_self(ref2ref_hits)
	
	singleton_bin = {}
	n_singleton = 0
	n_cross = 0
	for tx, h in ref2ref_hits.iteritems():
		if len(h) == 0:
			singleton_bin[tx] = True
			n_singleton += 1
		else:
			singleton_bin[tx] = False
			n_cross += 1
			
# 	read2ref_hits = read_blat_hits(args.read2ref, 'ref')
	report.write('Transcript\tCross\tLength\tContigs\tHits\tCovered\tLargest Covered\t# of reads\tCoverage\n')
	for rname in sorted_rnames:
		ref_len = len(ref.seqs[rname])
		report.write(rname + '\t')
		if len(ref2ref_hits[rname]) == 0:
			report.write('None\t')
		else:
			for h in ref2ref_hits[rname]:
				report.write(h.rname + ',')
			report.write('\t')
		report.write(str(ref_len) + '\t')
		if rname in hits:
			hit = hits[rname]
			report.write(str(len(hit)) + '\t')  # Number of contigs
			hit.sort(key=lambda x: x.rstart)
			longest_hit = 0
			covered = [0 for x in range(ref_len)]
			most_end = 0
			for h in hit:
				for i in range(h.n_blocks):
					report.write('[' + h.qname + ':' + str(h.r_block_starts[i]) + ',' + str(h.r_block_starts[i] + h.block_sizes[i]) + ']; ')
					if longest_hit == 0:
						longest_hit = h.block_sizes[i]
						most_end = h.r_block_starts[i] + h.block_sizes[i] - 1
					else:
						if longest_hit < h.block_sizes[i]:
							longest_hit = h.block_sizes[i]
					for j in range(h.r_block_starts[i], h.r_block_starts[i] + h.block_sizes[i]):
						covered[j] = 1
			covered_len = 0
			for i in range(ref_len):
				covered_len += covered[i]
			report.write('\t')
			covered_perc = covered_len / ref_len
			report.write('%.2f\t' % (covered_perc))
			largest_cover_perc = longest_hit / ref_len
			report.write('%.2f\t' % (largest_cover_perc))
# 			if rname in read2ref_hits:
# 				report.write('%d\t' % len(read2ref_hits[rname]))
# 				report.write('%.2f\t' % (len(read2ref_hits[rname]) / ref_len))
# 			else:
# 				report.write('0\t')
# 				report.write('0\t')
		else:
			report.write('0\tNone\t0\t0\t')
		report.write('\n')
	report.close()
	print 'Done. Check report file %s/report.txt' % args.out_dir

def read_psl_line(line):
	f = line.split('\t')
	hit = BlastHit(f[9], line)
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
	return hit

def read_psl_hits(hit_lines, key):
	hits = {}
	qname = ''
	reading_hits = False
	line_no = 0
	for line in hit_lines:
		line = line.strip()
		line_no += 1
		if len(line) < 4:
			continue
		if True:
			if not reading_hits:
				if '-----' in line:
					reading_hits = True
					continue 
				f = line.split('\t')
				try:
					n_match = int(f[0])
					reading_hits = True
				except:
					continue
			if not reading_hits:
				continue
			hit = read_psl_line(line)
			if key == 'query':
				if not hit.qname in hits:
					hits[hit.qname] = []
				hits[hit.qname].append(hit)
			else:
				if not hit.rname in hits:
					hits[hit.rname] = []
				hits[hit.rname].append(hit)
	return hits

def read_blat_hits(blat_fn, key='query'):
	blat = open(blat_fn, 'r')
	hit_lines = []
	for line in blat:
		hit_lines.append(line)
	hits = read_psl_hits(hit_lines, key)
	blat.close()
	return hits

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
	names = args.names.split(',')
	out_file = open(args.output, 'w')
	tx = FastaFile(args.ref)
	genes = {}
	if not len(files) == len(names):
		print 'Number of files and names are not the same. Abort.'
		return
	for i in range(len(files)):
		f = files[i]
		fp = open(f, 'r')
		for line in fp:
			line = line.strip().split('\t')[0]
			if line in genes:
				genes[line] += ',' + names[i]
			else:
				genes[line] = names[i]
		fp.close()
	categories = {}
	for key, value in sorted(genes.iteritems(), key=lambda (k, v): (v, k)):
		out_file.write(key + '\t' + value + '\n')
		if value in categories:
			categories[value] += 1
		else:
			categories[value] = 1
	for tx_name, seq in tx.seqs.iteritems():
		if not tx_name in genes:
			out_file.write(tx_name + '\tNone\n')
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
	summary.similarity_thre = args.similarity
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
	print 'Done.\n'
	ids.close()
	
def get_unaligned(args):
	contigs = FastaFile(args.contigs)
	hits = read_blat_hits(args.psl)
	unaligned = FastaFile()
	for ctg_name, seq in contigs.seqs.iteritems():
		if not ctg_name in hits:
			unaligned.seqs[ctg_name] = seq
	unaligned.save_to_disk('%s.unaligned' % args.psl)
	print 'Check %s.unaligned' % args.psl
	
def paths(args):
	all_hits = read_blat_hits(args.psl, 'ref')
	no = 1
	full_length = 0
	print >> sys.stderr, '%d reference transcripts hit' % len(all_hits)
#	print all_hits['ENST00000475972']
	for tx_name, hits in all_hits.iteritems():
		paths = []
		vertexes = []
		tx_len = 0
		hits.sort(key=lambda x: x.rstart)
		pre_h = None
		for h in hits:
			tx_len = h.rlen
			if h.n_mismatch <= 10 and h.n_match >= h.rlen - 10 and h.n_ref_gap == 0 and h.n_query_gap == 0:
				paths.append((h.qname, h.n_mismatch))
			else:
				if pre_h is None:
					pre_h = h
					vertexes.append(h)
				else:
					if abs(pre_h.rend - h.rstart) < 10:
						vertexes.append(h)
						pre_h = h
		if len(paths) > 0 or len(vertexes) > 0:
			print >> sys.stderr, no, tx_name, 'length:', tx_len
			no+=1
		if len(paths) > 0:
			full_length += 1
			for (p, n_mis) in paths:
				print >> sys.stderr, '\tPath %s: \t\t mismatch %d' % (p, n_mis)
		if len(vertexes) > 0:
			for h in vertexes:
				print >> sys.stderr, '\tVertex %s: [%d, %d] \t mismatch %d indels %d/%d' % (h.qname, h.rstart, h.rend, h.n_mismatch, h.n_ref_gap, h.n_query_gap)
	print >> sys.stderr, "\nFull-length: %d/%d" % (full_length, len(all_hits))

def diff_novo(args):
	tx = FastaFile(args.transcript)
	hits_1 = read_blat_hits(args.psl_1, 'ref')
	hits_2 = read_blat_hits(args.psl_2, 'ref')
	n_novo_base_1 = 0
	n_novo_base_2 = 0
	with open(args.psl_2 + '.novo', 'w') as novo:
		for tx_name, seq in tx.seqs.iteritems():
			if tx_name in hits_1 and not tx_name in hits_2:
				hs = hits_1[tx_name]
				l = 0
				for h in hs:
					l += h.alen
				novo.write('%s: %s %d\n' % (args.psl_1, tx_name, l))
				n_novo_base_1 += l
			if tx_name in hits_2 and not tx_name in hits_1:
				hs = hits_2[tx_name]
				l = 0
				for h in hs:
					l += h.alen
				novo.write('%s: %s\n' % (args.psl_2, tx_name, l))
				n_novo_base_2 += l
	print 'Novo base only 1: %d' % n_novo_base_1
	print 'Novo base only 2: %d' % n_novo_base_2
	print 'Check %s.novo' % args.psl_2
	
def cmp_psl(args):
	tx = FastaFile(args.transcript)
	hits_1 = read_blat_hits(args.psl_1, 'ref')
	hits_2 = read_blat_hits(args.psl_2, 'ref')
	n_base_2_file = open(args.psl_2 + '.base', 'w')
	n_base_1 = 0
	n_base_2 = 0
	n_base_both = 0
	n_valid = 0
	n_none = 0
	for tx_name, seq in tx.seqs.iteritems():
		tx_len = len(seq)
		bit_counter_1 = [0 for x in range(tx_len)]
		if tx_name in hits_1:
			for h in hits_1[tx_name]:
				for i in range(h.n_blocks):
					for j in range(h.r_block_starts[i], h.r_block_starts[i] + h.block_sizes[i]):
						bit_counter_1[j] = 1
					
		bit_counter_2 = [0 for x in range(tx_len)]
		if tx_name in hits_2:
			for h in hits_2[tx_name]:
				for i in range(h.n_blocks):
					for j in range(h.r_block_starts[i], h.r_block_starts[i] + h.block_sizes[i]):
						bit_counter_2[j] = 1
		
		n_base_2_on_tx = 0
		for i in range(tx_len):
			n_valid += 1
			if bit_counter_1[i] == 1 and bit_counter_2[i] == 1:
				n_base_both += 1
			if bit_counter_1[i] == 1 and bit_counter_2[i] == 0:
				n_base_1 += 1
			if bit_counter_1[i] == 0 and bit_counter_2[i] == 1:
				n_base_2 += 1
				n_base_2_on_tx += 1
			if bit_counter_1[i] == 0 and bit_counter_2[i] == 0:
				n_none += 1
				n_valid -= 1
		if n_base_2_on_tx >= 50:
			n_base_2_file.write('============== [%s, %d]: length %d ===============\n' % (tx_name, tx_len, n_base_2_on_tx))
			for i in range(tx_len):
				if bit_counter_1[i] == 0 and bit_counter_2[i] == 1:
					n_base_2_file.write('%d \n' % i)
	
	n_base_2_file.close()			
	print '%s: [psl_1]; %s: [psl_2]' % (args.psl_1, args.psl_2)
	print 'Bases by both:            %d' % n_base_both
	print 'Bases by psl_1 only:      %d' % n_base_1
	print 'Bases by psl_2 only:      %d' % n_base_2
	print 'Bases by none:            %d' % n_none
	if n_base_2 > 0:
		print 'psl_2 bases covered:      %.2f' % (n_base_both / (n_base_2 + n_base_both))
	print 'Check file: %s.base' % args.psl_2

def check_junctions(args):
	tx_hits = read_blat_hits(args.psl, 'query')
	

def comps(args):
    ref_hits = read_blat_hits(args.psl, 'ref')
    not_touched = []
    components = {}
    with open(args.csv) as csv:
        for line in csv:
            line = line.strip()
            fs = line.split("\t")
            components[fs[0]] = fs[1]
            
    n_covered = 0
    n_not_covered = 0
    n_abit_shorter = 0
    with open(args.id_file) as ids:
        for line in ids:
            ref = line.strip()
            hit_str = ref
            if not ref in ref_hits:
                hit_str += '\t0\t0\t0\n'
            else:
                comps = {}
                comp_hits = {}
                hits = ref_hits[ref]
                hit_str += '\t%d' % len(hits)
                n_vertexes = 0
                hits.sort(key=lambda x: x.rstart)
                for h in hits:
                    id_str = h.qname
                    fs = id_str.split("_")
                    comp_id = fs[0]
                    if comp_id in comps:
                        comps[comp_id] += 1
                        comp_hits[comp_id] = []
                    else:
                        n_vertexes += int(components[comp_id])
                        comps[comp_id] = 1
                    comp_hits[comp_id].append(h)
                hit_str += '\t%d' % len(comps)
                hit_str += '\t%d' % n_vertexes
                
                covered = False
                for comp_id, this_hits in comp_hits.iteritems():
                	if covered:
                		break
                	if len(this_hits) == 1:
                		h = this_hits[0]
                		if abs(h.rstart - h.rend) >= h.rlen - 50:
                			n_abit_shorter += 1
                			print h
                			break
                	else:
                		this_hits.sort(key=lambda x: x.rstart)
                		pre_h = None
                		for h in this_hits:
                			if pre_h is None:
                				pre_h = h
                				if h.rstart >= 20:
                					break
                			else:
                				if (h.rstart - pre_h.rend) > 20:
                					break
                			if h == this_hits[len(this_hits) - 1]:
                				if h.rlen - h.rend > 20:
                					break
                				n_covered += 1
                		  		covered = True
                		  		print this_hits
                if not covered:
                	n_not_covered += 1
                continue

                smaller_than_k = False
                gap_in_between = False
                if len(hits) > 0:
                    h = hits[0]
                    if h.rstart > 0:
                        gap_in_between = True

                for i in range(len(hits) - 1):
                    pre = hits[i]
                    post = hits[i + 1]
                    if post.rstart > pre.rend:
                        gap_in_between = True
                    if abs(pre.rend - post.rstart) < 25:
                        smaller_than_k = True 
                if smaller_than_k:
                    hit_str += '\t25'
                else:
                    hit_str += '\tMixed'
                if gap_in_between:
                    hit_str += '\tGap'
                else:
                    hit_str += '\tCovered'

                splicing_error = False
                if len(comps) < 5 and n_vertexes < 5:
                    for h in hits:
                        if h.n_ref_gap > 0 or h.n_query_gap > 0:
                            splicing_error = True
                if splicing_error:
                    hit_str += '\tSplicing error'
                else:
                    hit_str += '\tMixed'
                print hit_str
		print 'n_covered: %d' % (n_covered)
		print 'n_abit_shorter: %d' % (n_abit_shorter)
		print 'n_not_covered: %d' % (n_not_covered)

def analyze(args):
	in_path_ids = []
	with open(args.in_paths, 'r') as p:
		for line in p:
			line = line.strip()
			if len(line) > 0:
				in_path_ids.append(line) 
	ids = []
	with open(args.id_file, 'r') as i:
		for line in i:
			line = line.strip()
			if len(line) > 0:
				ids.append(line)
				
	print 'Totally %d missing transcripts' % len(ids)
	
	singletons = []
	with open(args.singleton) as single:
		for line in single:
			id = line.strip()
			singletons.append(id)
	missing_clusters = list(set(ids) - set(singletons))
	missing_singletons = list(set(ids) - set(missing_clusters))
	print 'Missing singletons: %d' % len(missing_singletons)
	for id in missing_singletons:
		print id
	print ''
    
	all_hits = {}
	counter = 0
	for id in ids:
		#if not id == 'SPBC8D2.01_T0':
		#	continue
		cmd = 'grep %s %s' % (id, args.read2tx)
		counter += 1
		print '%d:\t%s ...' % (counter, cmd)
		raw_lines = runInShell(cmd)
		hit_lines = raw_lines.split('\n')
		all_hits[id] = read_psl_hits(hit_lines, 'ref')[id]
	
	no_enough_coverage = {}
	for id in ids:
		#id = 'SPBC8D2.01_T0'
		all_tx_hits = all_hits[id]
		tx_hits = []
		for h in all_tx_hits:
			if h.n_match >= h.qlen - 2 or h.rstart == 0 or h.rend == h.rlen:
				tx_hits.append(h)
		tx_hits.sort(key=lambda x:x.rstart)
		base_coverage = pileup_coverage(tx_hits)
		h = tx_hits[0]
		for i in range(len(base_coverage)):
			b = base_coverage[i]
			if b == 0:
				no_enough_coverage[id] = 'Not covered: 0 at %d' % i
				break
				
		cursor = 0
		tx_hits.sort(key=lambda x:x.rstart)
		if not id in no_enough_coverage:
			for i in range(len(tx_hits) - 1):
				h = tx_hits[i]
				h_next = tx_hits[i + 1]
				if cursor > 0 and cursor - h_next.rstart < 11:
					no_enough_coverage[id] = 'Small overlapped region / Indels: [%d, %d)' % (h_next.rstart, h.rend)
					break
				else:
					tmp = h_next.rstart + h_next.block_sizes[0]
					if tmp > cursor: cursor = tmp
				#print h
				#print h_next
				#print 'Cursor: %d' % cursor
				#print '----'
		#break
	
	print '------------ No coverage: %d ---------------' % (len(no_enough_coverage))
	for id, des in no_enough_coverage.iteritems():
		print '%20s: \t%s' % (id, des)
	ids = list(set(ids) - set(no_enough_coverage.keys()))

	captured_in_paths = []
	for id in ids:
		if id in in_path_ids:
			captured_in_paths.append(id)
	print '------------ In path: %d ---------------' % len(captured_in_paths)
	for id in captured_in_paths:
		print id
	ids = list(set(ids) - set(in_path_ids))
	
	paired_hits = read_blat_hits(args.paired2tx, 'ref')
	not_covered = []
	for id in ids:
		if not id in paired_hits:
			not_covered.append(id)
	print '----------- No touched: %d ---------------' % len(not_covered)
	for id in not_covered:
		print id
	ids = list(set(ids) - set(not_covered))
	
	captured_in_paired = []
	for id in ids:
		hits = paired_hits[id]
		hits.sort(key=lambda x: x.alen, reverse=True)
		h = hits[0]
		if h.alen >= h.rlen * full_length_perc and h.n_bad_bases + h.n_query_gap_bases < bad_bases_thre:
			captured_in_paired.append(id)
	print '------------ In paired: %d ---------------' % len(captured_in_paired)
	for id in captured_in_paired:
		print id
	ids = list(set(ids) - set(captured_in_paired))
	
	should_jumped = []
	should_jumped_hits = {}
	for id in ids:
		tx_hits = paired_hits[id]
		tx_hits.sort(key=lambda x:x.rstart)
		h = tx_hits[0]
		should_jumped_hits[id] = tx_hits
		pre = h
		n_bad = h.n_bad_bases + h.n_query_gap_bases
		jumped = False
		if h.rstart <= 2 and n_bad <= 10 and h.qend >= h.qlen - 2:
			for i in range(1, len(tx_hits)):
				if id == 'SPBC18H10.16_T0':
					print h
				h = tx_hits[i]
				if h.rend >= h.rlen - 2:
					n_bad += h.n_bad_bases + h.n_query_gap_bases
					jumped = True
					break
				if h.rstart - pre.rend < 0 and h.qend >= h.qlen - 2:
					n_bad += h.n_bad_bases + h.n_query_gap_bases
				else:
					break
				pre = h
			if jumped:
				should_jumped.append(id)
	print '----------- Should jumped: %d ---------------' % len(should_jumped)
	for id in should_jumped:
		s = '%20s: \t%d hits\t' % (id, len(should_jumped_hits[id]))
		for h in should_jumped_hits[id]:
			s += '[%d, %d] ' % (h.rstart, h.rend)
		print s
	ids = list(set(ids) - set(should_jumped))
	
	has_repeat = []
	for id in ids:
		tx_hits = paired_hits[id]
		queries = []
		for h in tx_hits:
			if h.qname in queries:
				has_repeat.append(id)
				break
			else:
				queries.append(h.qname)
	print '----------- Has repeat: %d ---------------' % len(has_repeat)	
	for id in has_repeat:
		print id
	ids = list(set(ids) - set(has_repeat))
	
	full_but_gap_id = []
	for id in ids:
		tx_hits = paired_hits[id]
		tx_hits.sort(key=lambda x:x.alen, reverse=False)
		h = tx_hits[0]
		if h.alen >= h.rlen - 2:
			full_but_gap_id.append(id)
	print '----------- Full but gap in middle: %d ---------------' % len(full_but_gap_id)	
	for id in full_but_gap_id:
		s = '%20s: \t%d hits\t' % (id, len(paired_hits[id]))
		h = paired_hits[id][0]
		s += 'Mismatch: %d\t' % h.n_mismatch
		s += 'RefGap: %d,%d\t' % (h.n_ref_gap, h.n_ref_gap_bases)
		s += 'QueryGap: %d,%d\t' % (h.n_query_gap, h.n_query_gap_bases)
		print s
	ids = list(set(ids) - set(full_but_gap_id))
	
	annotated_tx = FastaFile(args.tx)
	base_cover = cal_base_cover(annotated_tx.seqs, paired_hits)
	partial = []
	ht_partial = {}
	head_tail_missing = []
	ht_missing_hits = {}
	for tx in ids:
		bases = base_cover[tx]
		n = 0.0
		for i in bases:
			n += i
		percentage = n / len(bases)
		if percentage < 0.99:
			partial.append(tx)
			ht_partial[tx] = paired_hits[tx][0]
	gap_in_middle = []
	ht_gap_m = {}
	for tx in partial:
		tx_hits = paired_hits[tx]
		tx_hits.sort(key=lambda x:x.n_match, reverse=True)
		h = tx_hits[0]
		if (h.rstart < 100 and h.rend > (h.rlen - 100)) and h.n_mismatch < 10 and (h.n_query_gap_bases + h.n_bad_bases) < 10:
			head_tail_missing.append(tx)
			ht_missing_hits[tx] = h
		else:
			if h.rend - h.rstart > h.rlen * 0.99:
				gap_in_middle.append(tx)
				ht_gap_m[tx] = h
	ids = list(set(ids) - set(partial))
	print '----------- Covered <0.99: %d ---------------' % len(partial)
	partial = list(set(partial) - set(head_tail_missing))
	partial = list(set(partial) - set(gap_in_middle)) 
	for id in partial:
		h = ht_partial[id]
		r = 'Range [%d, %d] / %d' % (h.rstart, h.rend, h.rlen)
		print '%20s: \tBadBases: %d+%d \t %30s; mismatches: %d; gaps: %d ' % (id, h.n_bad_bases, h.n_query_gap_bases, r, h.n_mismatch, h.n_query_gap + h.n_ref_gap)
	print '\t--- Head/Tail Missing: %d ---' % len(head_tail_missing)
	for id in head_tail_missing:
		h = ht_missing_hits[id]
		print '\t%20s: \tBadBases: %d \t Range [%d, %d] / %d; mismatches: %d; gaps: %d ' % (id, h.n_bad_bases, h.rstart, h.rend, h.rlen, h.n_mismatch, h.n_query_gap + h.n_ref_gap)
	print '\t--- Gap in the middle: %d ---' % len(gap_in_middle)
	for id in gap_in_middle:
		h = ht_gap_m[id]
		print '\t%20s: \tRefGap: %d,%dbp; QueryGap: %d,%d; Mismatches: %d' % (id, h.n_ref_gap, h.n_ref_gap_bases, h.n_query_gap, h.n_query_gap_bases, h.n_mismatch)
	
	missing_junctions = {}
	n = 0
	for tx in ids:
		tx_hits = paired_hits[tx]
		tx_hits.sort(key=lambda x:x.rstart)
		h = tx_hits[0]
		tx_len = h.rlen
		start = -1
		end = 99999999
		missing_junctions[tx] = []
		should_jumped_merged = []
		for h in tx_hits:
			if h.rstart > end + 2:
				continue
			if not end == 99999999 and h.rend < end:
				continue
			if h.n_mismatch < 10 and (h.n_query_gap_bases + h.n_bad_bases) < 10:
				if start == -1:
					start = h.rstart
				end = h.rend
				missing_junctions[tx].append(h)
		n_covered = abs(end - start) 
		if n_covered / tx_len < 0.99:
			missing_junctions[tx] = []
		else:
			n += 1
	
	print '-------------- Missing junctions: %d ---------------' % n
	tmp = []
	for tx, hits in missing_junctions.iteritems():
		if len(hits) > 0:
			n += 1
			tmp.append(tx)
			n_bad = 0
			for h in hits:
				n_bad += h.n_bad_bases
			line = '%20s: \tBadBases: %d \t Ranges: %d\t' % (tx, n_bad, hits[0].rlen)
			for h in hits:
				line += '[%d, %d] ' % (h.rstart, h.rend)
			print line
	
	ids = list(set(ids) - set(tmp))		
	print '-------------- Others: %d -----------------' % len(ids)
	for id in ids:
		print id

def quality(args):
    ids = {}
    with open(args.id_file) as id_f:
        for line in id_f:
            line = line.strip()
            if len(line) > 1:
                ids[line] = 0
    counter = 0
    n_bad_hit = 0
    n_good_hit = 0
    n_not_full = 0
    n_no_hit = 0
    with open(args.start) as starts:
        for line in starts:
            line = line.strip()
            fs = line.split(', ')
            if len(fs) < 4:
                continue
            counter += 1
            read_id = fs[3]
            print '---------- Template %s, Read %s ------------' % (fs[2], read_id)
            cmd = 'grep %s %s' % (read_id, args.read2tx)
            hit_lines = runInShell(cmd).split('\n')
            raw_hits = read_psl_hits(hit_lines, 'query')
            if len(raw_hits) <= 0 or not read_id in raw_hits:
                print 'No hits'
                n_no_hit += 1
                continue
            #print raw_hits
            hits = raw_hits[read_id]
            for h in hits:
                if h.n_match >= h.qlen - 2:
                    print 'Good hit on %s' % h.rname
                    n_good_hit += 1
                else:
                    if h.rstart == 0 and h.n_match >= h.qlen - h.qstart - 2:
                        print 'Good hit on %s' % h.rname
                        n_good_hit += 1
                    elif h.rend == h.rlen and h.n_match >= h.qlen - h.qend - 2:
                        print 'Good hit on %s' % h.rname
                        n_good_hit += 1
                    else:
                        print 'Bad hit on %s' % h.rname
                        n_bad_hit += 1
                        print h
                        if h.rname in ids.keys():
                            print 'Not full-length'
                            cmd = 'grep %s %s' % (h.rname, args.paired2tx)
                            tmp = runInShell(cmd).split('\n')
                            paired_raw_hits = read_psl_hits(tmp, 'ref')
                            paired_hits = paired_raw_hits[h.rname]
                            print paired_hits
                            print h
                            for p in paired_hits:
                                if p.rend == h.rend:
                                    print 'Misleading bad hit' 
                                    cmd = 'zoom ctx -t %s -f spombe.broad.tx.fasta.rev ../SRR097897_half/paired.fa.psl' % h.rname
                                    print cmd
                                    cmd = 'zoom rtx -r SRR097897.half.fa -f spombe.broad.tx.fasta.rev -p SRR097897.half.fa.psl %s' % h.rname
                                    print cmd
                                    break
                            if ids[h.rname] == 0:
                                n_not_full += 1
                                ids[h.rname] = 1
    print 'Total: %d' % counter
    print 'Good hits: %d' % n_good_hit
    print 'No hits: %d' % n_no_hit
    print 'Bad hit: %d' % n_bad_hit
    print 'Not full-length: %d' % n_not_full

def main():
    parser = ArgumentParser()
    subparsers = parser.add_subparsers(help='sub command help')
    parser_differ = subparsers.add_parser('diff', help='differ multiple files')
    parser_differ.set_defaults(func=differ)
    parser_differ.add_argument('ref', help='annotated transcripts')
    parser_differ.add_argument('-i', required=True, help='the files to compare, seperated by a ","', dest='inputs')
    parser_differ.add_argument('-n', required=True, help='name of the compared assemblers, seperated by a ","', dest='names')
    parser_differ.add_argument('-o', required=True, help='result file', dest='output')

    parser_blat = subparsers.add_parser('blat', help='evaluate the performance by aligning in blat')
    parser_blat.set_defaults(func=eva_blat)
    parser_blat.add_argument('-t', required=True, help='reference transcript file', dest='ref')
    parser_blat.add_argument('-c', required=True, help='transcripts reported to be evaluated', dest='contigs')
    parser_blat.add_argument('-p', required=True, help='psl file', dest='psl')
    parser_blat.add_argument('-o', required=True, help='output folder', dest='out_dir')
    parser_blat.add_argument('-s', required=True, help='similarity', dest='similarity')

    parser_cmp = subparsers.add_parser('hasdup', help='evaluate the performance by alinging by Blastn')
    parser_cmp.set_defaults(func=check_dup)
    parser_cmp.add_argument('-f', required=True, help='file', dest='input')
    
    parser_sta = subparsers.add_parser('sta', help='Report statistics about the contigs')
    parser_sta.set_defaults(func=stat)
    parser_sta.add_argument('ref', help='annotated transcripts')
    parser_sta.add_argument('out_dir', help='output directory')
    parser_sta.add_argument('-c', '--ctg2ref', required=True, metavar='FILE', help='PSL file: contigs aligned to annotated transcripts', dest='ctg2ref')
    parser_sta.add_argument('-t', '--ref2ref', required=True, metavar='FILE', help='PSL file: annotated transcripts aligned to itself', dest='ref2ref')
    parser_sta.add_argument('-r', '--read2ref', required=False, metavar='FILE', help='PSL file: raw reads aligned to annotated transcripts', dest='read2ref')
    parser_sta.add_argument('-s', required=True, help='similarity', dest='similarity')
    
    parser_psl = subparsers.add_parser('psl', help='Compare two PSL hit files')
    parser_psl.set_defaults(func=cmp_psl)
    parser_psl.add_argument('transcript', help='annotated transcripts')
    parser_psl.add_argument('-1', '--psl_1', required=True, metavar='FILE', help='PSL file 1', dest='psl_1')
    parser_psl.add_argument('-2', '--psl_2', required=True, metavar='FILE', help='PSL file 2', dest='psl_2')
    
    parser_unaligned = subparsers.add_parser('unaligned', help='Get only unaligned contigs')
    parser_unaligned.set_defaults(func=get_unaligned)
    parser_unaligned.add_argument('contigs', help='assembled contigs')
    parser_unaligned.add_argument('-p', '--psl', required=True, metavar='FILE', help='PSL file', dest='psl')
    
    parser_paths = subparsers.add_parser('paths', help='Check how many full length paths')
    parser_paths.set_defaults(func=paths)
    parser_paths.add_argument('psl', help='path-to-ref PSL file')
    
    parser_comps = subparsers.add_parser('comp', help='Check which transcript from which component')
    parser_comps.set_defaults(func=comps)
    parser_comps.add_argument('psl', help='component-to-ref PSL file')
    parser_comps.add_argument('id_file', help='file of id list to check')
    parser_comps.add_argument('csv', help='components.csv')
    
    parser_quality = subparsers.add_parser('qua', help='Check quality of starting reads when jumping')
    parser_quality.set_defaults(func=quality)
    parser_quality.add_argument('id_file', help='file of id list to check')
    parser_quality.add_argument('read2tx', help='read-to-transcript psl')
    parser_quality.add_argument('paired2tx', help='paired-to-transcript psl')
    parser_quality.add_argument('start', help='file with starting reads id')
    
    parser_analyze = subparsers.add_parser('any', help='Analyse missing transcripts')
    parser_analyze.set_defaults(func=analyze)
    parser_analyze.add_argument('id_file', help='file of id list to check')
    parser_analyze.add_argument('singleton', help='file of id list of singleton transcripts')
    parser_analyze.add_argument('tx', help='reference transcript file')
    parser_analyze.add_argument('in_paths', help='ids captured in paths')
    parser_analyze.add_argument('read2tx', help='read-to-transcript psl')
    parser_analyze.add_argument('paired2tx', help='paired-to-transcript psl')
    parser_analyze.add_argument('junctions', help='paired.junctions')
    
    args = parser.parse_args()
    args.func(args)

if __name__ == '__main__':
    sys.exit(main())
