#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <glib.h>
#include <math.h>
#include "rand.h"
#include "ass.h"
#include "bwase.h"
#include "utils.h"
#include "peseq.h"
#include "pool.h"
#include "pechar.h"
#include "pet.h"
#include "pealn.h"
#include "roadmap.h"
#include "edge.h"
#include "edgelist.h"
#include "readrm.h"
#include "pechar.h"

/**
 * PET is not used.
 *
 * Original contig: abcdefghijklmnopqrst
 * PET: bcd, opqr
 * Process: start from 'bcd', extend to the right;
 * 			start from 'opqr', extend to the left;
 * 			try to combine them.
 * Following methods are different when dealing with the two PET mates:
 * 1. initial query;
 * 2. extension;
 * 3. forwarding;
 * 4. query expansion.
 * The two component are flagged by the global variable left_max_ctg_id.
 */
unsigned int contig_id = 0;
// unsigned int left_max_ctg_id = 0; // Record the largest contig id after assembling the left mate.
roadmap *left_rm = NULL;
roadmap *right_rm = NULL;
edgearray *all_edges = NULL;
ass_opt *opt = 0;
int n_reads_consumed = 0;

ass_opt *init_opt() {
	ass_opt *o = (ass_opt*) calloc(1, sizeof(ass_opt));
	o->nm = 2;
	o->ol = 50;
	o->rl = 100;
	o->n_alg = 30;
	o->mean = 0;
	o->sd = 0;
	o->pair = 0;
	o->nsp = 1;
	o->mode = (BWA_MODE_GAPE | BWA_MODE_COMPREAD);
	o->solid_reads = 0;
	return o;
}

ext_msg *new_msg() {
	ext_msg *m = (ext_msg*) malloc(sizeof(ext_msg));
	m->data = 0;
	m->type = 0;
	m->message = 0;
	m->query = 0;
	return m;
}

void free_msg(ext_msg *m) {
	if (m) {
		free(m->message);
		bwa_free_read_seq(1, m->query);
		free(m);
		m = 0;
	}
}

int check_ign(const int ori, bwa_seq_t *s, bwa_seq_t *contig) {
	int ignore = 0;
	if (ori) {
		if (s->rev_com) {
			if ((s->rseq[s->cursor + 1] != contig->seq[contig->len - 1]))
				ignore = 1;
		} else {
			if (s->seq[s->cursor + 1] != contig->seq[contig->len - 1])
				ignore = 1;
		}
	} else {
		if (s->rev_com) {
			if ((s->rseq[s->cursor - 1] != contig->seq[contig->len - 1]))
				ignore = 1;
		} else {
			if (s->seq[s->cursor - 1] != contig->seq[contig->len - 1])
				ignore = 1;
		}
	}
	return ignore;
}

void filter_pool(pool *p, edge *ass_eg, const hash_table *ht) {
	readarray *reads = p->reads;
	bwa_seq_t *r = 0, *used_mate;
	int i = 0, range_l = 0, range_h = 0;
	for (i = 0; i < reads->len; i++) {
		r = g_ptr_array_index(reads, i);
		// When extending to the left, the sequence is reversed first!
		// So the actual shift value is: length of contig - shift + 1
		range_l = (ass_eg->len - opt->mean - opt->sd * SD_TIMES);
		range_h = (ass_eg->len - 1);
		used_mate = get_mate(r, ht->seqs);
		if (r->used || used_mate->shift < range_l || used_mate->shift > range_h) {
			pool_rm_index(p, i);
			i--;
		}
	}
}

/**
 * To increase the reliability, only if the last bases of current contig and the aligned seq
 * are identical, we extend to next base.
 */
void upd_cur_pool(const alignarray *alns, int *next, pool *cur_pool,
		pool *mate_pool, bwa_seq_t *query, const hash_table *ht, edge *ass_eg,
		const int ori) {
	int i = 0;
	int check_c_1 = 0, check_c_2 = 0, confirm_c = 0, confirm_c_2 = 0;
	index64 index;
	alg *a;
	gboolean removed;
	bwa_seq_t *s, *seqs, *contig = ass_eg->contig, *mate, *q = query;
	seqs = ht->seqs;
	for (i = 0; i < alns->len; i++) {
		a = g_ptr_array_index(alns, i);
		index = a->r_id;
		if (index >= ht->n_seqs)
			continue;
		s = &seqs[index];
		s->rev_com = a->rev_comp;
		mate = get_mate(s, seqs);
		mate->rev_com = s->rev_com;
		if (ass_eg->len > opt->mean + opt->sd * SD_TIMES + opt->rl) {
			if (ori && is_left_mate(s->name)) {
				if ((mate && mate->contig_id != ass_eg->id))
					continue;
			}
			if (!ori && is_right_mate(s->name)) {
				if ((mate && mate->contig_id != ass_eg->id))
					continue;
			}
		}
		if (s->is_in_c_pool || (a->pos + opt->ol - 3) > s->len)
			continue;
		// cursor points to the next char
		if (s->rev_com)
			s->cursor = ori ? (s->len - opt->ol - 1 - a->pos) : (s->len
					- a->pos);
		else
			s->cursor = ori ? (a->pos - 1) : (a->pos + opt->ol);
		// If cursor reaches the end, ignore it
		if (s->cursor >= s->len || s->cursor < 0) {
			s->cursor = 0;
			continue;
		}

		if (check_ign(ori, s, contig)) {
			s->cursor = 0;
			continue;
		}

		// printf("[extend] Inserting alignment %s:%d into current pool. \n", p->r_id, p->pos);
		// p_pool("Before inserting: ", r_pool);
		pool_uni_add(cur_pool, s);
		//pool_add(mate_pool, mate);
		// mate->is_in_c_pool = 0; // mate read is not in current pool
	}

	// Remove those mate reads that would not be useful.
	//	if (ass_eg->len % opt->rl == 0) {
	//		filter_pool(mate_pool, ass_eg, ht);
	//	}
	// In case the alignment only finds partial results,
	//   from the mate pool add reads
	//	for (i = 0; i < mate_pool->n; i++) {
	//		mate = g_ptr_array_index(mate_pool->reads, i);
	//		if (mate->used)
	//			continue;
	//		sub_index = is_sub_seq(q, 0, mate, opt->nm, 0);
	//		if (sub_index != NOT_FOUND) {
	//			mate->cursor = ori ? (sub_index - 1) : (sub_index + opt->ol);
	//			if (mate->cursor >= mate->len || mate->cursor < 0) {
	//				mate->cursor = 0;
	//				continue;
	//			}
	//			// Add the mate to the current pool and remove from the mate pool
	//			pool_uni_add(cur_pool, mate);
	//			pool_rm_index(mate_pool, i);
	//			i--;
	//		}
	//	}
	check_c_1 = ori ? query->seq[0] : query->seq[query->len - 1];
	check_c_2 = ori ? query->seq[1] : query->seq[query->len - 2];
	for (i = 0; i < cur_pool->n; i++) {
		s = g_ptr_array_index(cur_pool->reads, i);
		confirm_c = ori ? s->seq[s->cursor + 1] : s->seq[s->cursor - 1];
		confirm_c_2 = ori ? s->seq[s->cursor + 2] : s->seq[s->cursor - 2];
		if (s->rev_com) {
			confirm_c = ori ? s->rseq[s->cursor + 1] : s->rseq[s->cursor - 1];
			confirm_c_2 = ori ? s->rseq[s->cursor + 2] : s->rseq[s->cursor - 2];
		}
		// Remove those reads probably at the splicing junction
		if ((is_sub_seq(q, 0, s, opt->nm, 0) == NOT_FOUND && is_sub_seq_byte(
				q->rseq, q->len, 0, s, opt->nm, 0) == NOT_FOUND) || (check_c_1
				!= confirm_c && check_c_2 != confirm_c_2)) {
			mate = get_mate(s, seqs);
			removed = pool_rm_index(cur_pool, i);
			//pool_rm_fast(mate_pool, mate);
			if (removed)
				i--;
		} else {
			if (s->rev_com)
				check_c(next, s->rseq[s->cursor]);
			else
				check_c(next, s->seq[s->cursor]);
		}
	}
}

/**
 * Get a read in the pool, which is at least overlapping with ONE other read
 */
readarray *get_parents_reads(const edge *parent, const int ori) {
	readarray *reads = g_ptr_array_sized_new(INIT_EDGE_SPACE);
	int acc_len = parent->len;
	edge *pre = 0;
	if (ori) {
		if (parent->out_egs->len > 0)
			pre = g_ptr_array_index(parent->out_egs, 0);
	} else {
		if (parent->in_egs->len > 0)
			pre = g_ptr_array_index(parent->in_egs, 0);
	}
	g_ptr_array_concat(reads, parent->reads);
	while (pre && acc_len < opt->mean + opt->sd * SD_TIMES) {
		g_ptr_array_concat(reads, pre->reads);
		acc_len += pre->len;
		if (ori) {
			if (pre->out_egs->len > 0)
				pre = g_ptr_array_index(pre->out_egs, 0);
			else
				pre = 0;
		} else {
			if (pre->in_egs->len > 0)
				pre = g_ptr_array_index(pre->in_egs, 0);
			else
				pre = 0;
		}
	}
	if (acc_len <= opt->rl) {
		return 0;
	}
	return reads;
}

bwa_seq_t *get_query_ol(edge *ass_eg, bwa_seq_t *seqs, pool *m_pool,
		const int ori) {
	bwa_seq_t *read, *mate;
	//readarray *p_ra;
	int start, start_copy, i = 0;
	if (!ass_eg)
		return 0;
	start = get_mid_pos(ass_eg->reads, ori, opt->mean);
	start_copy = start;
	//p_readarray(ass_eg->reads, 1);
	show_debug_msg(__func__, "Starting point: %d\n", start);
	// From start to the end
	start = (start < 0) ? 0 : start;
	while (start < ass_eg->reads->len) {
		read = g_ptr_array_index(ass_eg->reads, start);
		if ((ori && is_right_mate(read->name)) || (!ori && is_left_mate(
				read->name))) {
			mate = get_mate(read, seqs);
			if (!mate->used && !has_n(mate))
				return mate;
		}
		start++;
	}
	// From start to the beginning
	start = start_copy;
	while (start >= 0) {
		read = g_ptr_array_index(ass_eg->reads, start);
		if ((ori && is_right_mate(read->name)) || (!ori && is_left_mate(
				read->name))) {
			mate = get_mate(read, seqs);
			if (!mate->used)
				return mate;
		}
		start--;
	}

	for (i = 0; i < m_pool->reads->len; i++) {
		read = g_ptr_array_index(m_pool->reads, i);
		if ((!ori && is_right_mate(read->name)) || (ori && is_left_mate(
				read->name))) {
			if (!read->used)
				return read;
		}
	}
	return 0;
}

/**
 * Get the reads used by its parents/children.
 * NOTE: Only go to the first parent/child
 */
int get_parents_len(edge *parent, const int ori) {
	int acc_len = parent->len;
	edge *pre = 0;
	edgearray *egs = g_ptr_array_sized_new(16);
	if (ori) {
		if (parent->out_egs->len > 0)
			pre = g_ptr_array_index(parent->out_egs, 0);
	} else {
		if (parent->in_egs->len > 0)
			pre = g_ptr_array_index(parent->in_egs, 0);
	}
	g_ptr_array_add(egs, parent);
	while (pre) {
		if (edgearray_find(egs, pre) != NOT_FOUND)
			break;
		g_ptr_array_add(egs, pre);
		acc_len += pre->len;
		if (ori) {
			if (pre->out_egs->len > 0)
				pre = g_ptr_array_index(pre->out_egs, 0);
			else
				pre = 0;
		} else {
			if (pre->in_egs->len > 0)
				pre = g_ptr_array_index(pre->in_egs, 0);
			else
				pre = 0;
		}
	}
	g_ptr_array_free(egs, TRUE);
	return acc_len;
}

/**
 * If within_range is 1: return every read;
 * 		otherwise: return reads within the range;
 *
 * If used is 0: return those not used yet;
 * 		otherwise: return whatever got.
 */
pool *get_mate_pool(const edge *eg, const hash_table *ht, const int ori,
		const int within_range, const int used) {
	pool *mate_pool = new_pool();
	int start = 0, end = 0, i = 0, max_rg = 0;
	bwa_seq_t *seqs, *s, *mate;
	readarray *reads = 0;
	seqs = ht->seqs;
	reads = get_parents_reads(eg, ori);

	if (!reads)
		return mate_pool;
	// p_readarray(reads, 1);
	max_rg = (opt->mean + SD_TIMES * opt->sd);

	if (ori) {
		start = 0;
		end = max_rg;
	} else {
		start = (eg->len > max_rg) ? (eg->len - max_rg) : 0;
		end = eg->len;
	}
	//show_debug_msg(__func__, "start = %d, end = %d \n", start, end);
	// Add the mates of the reads in the middle first
	for (i = 0; i < reads->len; i++) {
		s = g_ptr_array_index(reads, i);
		// If not range requirement, just go into the 'if' block;
		// If range is required (within_range = 1), check the shift value as well
		if (!within_range || (s->shift >= start && s->shift < end)) {
			if (ori)
				mate = get_left_mate(s, seqs);
			else
				mate = get_right_mate(s, seqs);
			if (mate && (used || !mate->used))
				pool_uni_add(mate_pool, mate);
		}
	}
	//p_pool("Mate pool: ", mate_pool);
	return mate_pool;
}

pool *get_init_pool(const hash_table *ht, bwa_seq_t *init_read, const int ori) {
	alignarray *alns = g_ptr_array_sized_new(N_DEFAULT_ALIGNS);
	int i = 0;
	pool *init_pool = new_pool();
	bwa_seq_t *seqs = ht->seqs, *s;
	alg *a;
	p_query(__func__, init_read);
	pe_aln_query(init_read, init_read->seq, ht, opt->nm + 2, opt->ol, 0, alns);
	if (opt->nsp) {
		pe_aln_query(init_read, init_read->rseq, ht, opt->nm + 2, opt->ol, 1,
				alns);
	}
	for (i = 0; i < alns->len; i++) {
		a = (alg*) g_ptr_array_index(alns, i);
		s = &seqs[a->r_id];
		s->rev_com = a->rev_comp;
		if (s->rev_com)
			s->cursor = ori ? (s->len - opt->ol - 1 - a->pos) : (s->len
					- a->pos);
		else
			s->cursor = ori ? (a->pos - 1) : (a->pos + s->len);
		if (s->used || s->cursor < 0 || s->cursor >= s->len) {
			s->cursor = 0;
			continue;
		}
		pool_add(init_pool, s);
	}
	//p_pool("Initial Pool: ", init_pool);
	free_alg(alns);
	return init_pool;
}

int vld_ext(edge *parent, bwa_seq_t *query, const hash_table *ht, const int ori) {
	bwa_seq_t *seqs = ht->seqs, *s, *mate;
	int i = 0, n_on = 0, acc_len = 0, q_pos_on_mate;
	ubyte_t cursor_char, next_char;
	readarray *reads = 0;
	reads = get_parents_reads(parent, ori);
	acc_len = get_parents_len(parent, ori);
	if (!reads || acc_len < opt->mean) {
		show_debug_msg(__func__, "Accumulated length too short: %d \n", acc_len);
		return 1;
	}

	//p_readarray(reads, 1);
	cursor_char = ori ? query->seq[0] : query->seq[query->len - 1];
	p_query(__func__, query);
	for (i = 0; i < reads->len; i++) {
		s = g_ptr_array_index(reads, i);
		//p_query(__func__, s);
		mate = get_mate(s, seqs);
		if (mate && !mate->used) {
			if (s->rev_com) {
				q_pos_on_mate = is_sub_seq_byte(query->rseq, query->len, 0,
						mate, opt->nm, 0);
				if (q_pos_on_mate != NOT_FOUND) {
					next_char = ori ? mate->rseq[q_pos_on_mate]
							: mate->rseq[q_pos_on_mate + mate->len - 1];
					//					show_debug_msg(__func__,
					//							"next_char: %d; cursor_char: %d\n", next_char,
					//							cursor_char);
					if (next_char == cursor_char)
						n_on++;
				}
			} else {
				q_pos_on_mate = is_sub_seq(query, 0, mate, opt->nm, 0);
				if (q_pos_on_mate != NOT_FOUND) {
					next_char = ori ? mate->seq[q_pos_on_mate]
							: mate->seq[q_pos_on_mate + query->len - 1];
					//					show_debug_msg(__func__,
					//							"next_char: %d; cursor_char: %d\n", next_char,
					//							cursor_char);
					if (next_char == cursor_char)
						n_on++;
				}
			}
		}
	}

	show_debug_msg(__func__, "n_on: %d \n", n_on);
	g_ptr_array_free(reads, TRUE);
	if (n_on >= N_MIN_VAL_EXT)
		return 1;
	else
		return 0;
}

void rev_reads_pos(edge *eg) {
	int i = 0;
	readarray *reads = eg->reads;
	bwa_seq_t *r;
	for (i = 0; i < reads->len; i++) {
		r = g_ptr_array_index(reads, i);
		r->shift = eg->len - r->shift + 1;
	}
}

void add_pool_by_ol(pool *p, bwa_seq_t *query, bwa_seq_t *read, const int ori) {
	int confirm_c = 0, confirm_c_2 = 0, index = 0, cursor = 0, check_c = 0,
			check_c_2 = 0;

	if (read->used)
		return;

	confirm_c = ori ? query->seq[0] : query->seq[query->len - 1];
	confirm_c_2 = ori ? query->seq[1] : query->seq[query->len - 2];
	//	show_debug_msg(__func__, "confirm_c: %d; confirm_c_2: %d \n", confirm_c, confirm_c_2);
	//	p_query(__func__, query);
	//	p_query(__func__, read);
	if (read->rev_com) {
		index = is_sub_seq_byte(query->rseq, query->len, 0, read, opt->nm, 0);
		if (index > 0) {
			cursor = ori ? (read->len - opt->ol - index - 1) : (query->len
					- index);
			check_c = ori ? read->rseq[cursor + 1] : read->rseq[cursor - 1];
			check_c_2 = ori ? read->rseq[cursor + 2] : read->rseq[cursor - 2];
			if (cursor >= 0 && cursor < read->len && confirm_c == check_c
					&& confirm_c_2 == check_c_2) {
				read->cursor = cursor;
				pool_uni_add(p, read);
			}
		}
	} else {
		index = is_sub_seq(query, 0, read, opt->nm, 0);
		if (index > 0) {
			cursor = ori ? (index - 1) : (index + query->len);
			check_c = ori ? read->seq[cursor + 1] : read->seq[cursor - 1];
			check_c_2 = ori ? read->seq[cursor + 2] : read->seq[cursor - 2];
			if (cursor >= 0 && cursor < read->len && confirm_c == check_c
					&& confirm_c_2 == check_c_2) {
				read->cursor = cursor;
				pool_uni_add(p, read);
			}
		}
	}
}

void fill_in_pool(pool *p, edgearray *reads, bwa_seq_t *query, const int ori) {
	int i = 0;
	bwa_seq_t *read;
	//p_query(query);
	for (i = 0; i < reads->len; i++) {
		read = g_ptr_array_index(reads, i);
		if (!read->used)
			add_pool_by_ol(p, query, read, ori);
	}
}

void adj_shift(edge *eg, const int trun_len) {
	int i = 0;
	bwa_seq_t *r;
	if (!eg || trun_len < 0)
		return;
	for (i = 0; i < eg->reads->len; i++) {
		r = g_ptr_array_index(eg->reads, i);
		r->shift -= trun_len;
	}
}

void forward_by_ra(edge *ass_eg, pool *cur_pool, const readarray *ra,
		bwa_seq_t *query, const int ori) {
	int next_char = 0, i = 0;
	bwa_seq_t *read;
	next_char = get_next_char(cur_pool, ori, ass_eg);
	ext_que(query, next_char, ori);
	ext_con(ass_eg->contig, next_char, ori);
	ass_eg->len = ass_eg->contig->len;
	//p_query(__func__, query);
	for (i = 0; i < ra->len; i++) {
		read = g_ptr_array_index(ra, i);
		add_pool_by_ol(cur_pool, query, read, ori);
	}
	for (i = 0; i < cur_pool->n; i++) {
		read = g_ptr_array_index(cur_pool->reads, i);
		// If the read has no overlaps with
		if ((!read->rev_com && is_sub_seq(query, 0, read, opt->nm, 0)
				== NOT_FOUND) || (read->rev_com && is_sub_seq_byte(query->rseq,
				query->len, 0, read, opt->nm, 0) == NOT_FOUND)) {
			pool_rm_index(cur_pool, i);
			i--;
		}
	}
	ass_eg->len = ass_eg->contig->len;
}

int est_gap(const edge *left_eg, const edge *right_eg, const hash_table *ht) {
	readarray *l_par_reads, *r_par_reads, *ea;
	bwa_seq_t *read, *mate;
	int i = 0, l_index = 0, r_index = 0, len = 0;
	int p_count = 0, gap_sum = 0, gap_ave = 0;
	show_debug_msg(__func__, "left_eg = [%d, %d], right_eg = [%d, %d] \n",
			left_eg->id, left_eg->len, right_eg->id, right_eg->len);
	l_par_reads = get_parents_reads(left_eg, 0);
	r_par_reads = get_parents_reads(right_eg, 1);
	ea = get_paired_reads(l_par_reads, r_par_reads, ht->seqs, 0);
	if (ea->len == 0) {
		show_debug_msg(__func__, "No reads spanning the two edges\n");
		return INVALID;
	}
	//p_readarray(ea, 1);
	for (i = 0; i < ea->len; i++) {
		read = g_ptr_array_index(ea, i);
		mate = get_mate(read, ht->seqs);
		l_index = read->shift;
		r_index = mate->shift;
		//		l_index = is_sub_seq(read, 0, left_eg->contig, opt->nm, 0);
		//		r_index = is_sub_seq(mate, 0, right_eg->contig, opt->nm, 0);
		//		if (l_index == NOT_FOUND || r_index == NOT_FOUND) {
		//			continue;
		//		}
		len = (left_eg->len - l_index) + r_index;
		//show_debug_msg(__func__, "left = %d, right = %d \n", l_index, r_index);
		gap_sum += opt->mean - len;
		p_count++;
	}
	//show_debug_msg(__func__, "gan_sum = %d, p_count = %d \n", gap_sum, p_count);
	gap_sum = (gap_sum > 0) ? gap_sum : 0;
	if (p_count > 0) {
		gap_ave = gap_sum / p_count;
		if (gap_ave > (opt->mean + opt->sd * SD_TIMES))
			return INVALID;
		return gap_ave;
	} else
		return INVALID; // If no pairs are counted, is because some reads are from parents.
}

/**
 * Fill in the gap between left_eg and right_eg;
 * Mind that when calling this function, left_eg should come first;
 * Parameter ori indicates where to merge the seq to: 0 means to left, 1 means to right.
 *
 * At the very beginning, check whether there are paired-end reads spanning the left_eg and right_eg.
 *
 * First get right mate pool of the left_eg, trying to extend the right_eg to the left;
 * Then get left mate pool of right_eg, trying to extend the left_eg to the right;
 * When there are enough overlapping, merge seqs;
 * When there are no enough overlapping, estimate the gap, fill in with N's.
 * When the gap is too large, abandon.
 */
void fill_in_gap(edge *left_eg, edge *right_eg, const int reason_gap,
		const hash_table *ht, const int ori) {
	readarray *paired_reads, *mates = 0, *left_p_reads, *right_p_reads;
	pool *cur_pool, *ass_mate_pool, *left_pool = 0;
	int olpped = 0, gap = 0, ol_len = 0, to_merge = 0, ori_len = 0;
	int ol = 0;
	bwa_seq_t *query = 0;
	eg_gap *added_gap;
	//	p_query(left_eg->contig);
	//	p_query(right_eg->contig);
	//p_ctg_seq(__func__, left_eg->contig);
	//p_ctg_seq(__func__, right_eg->contig);
	ol_len = find_ol(left_eg->contig, right_eg->contig, opt->nm + 2);
	show_debug_msg(__func__, "Filling the gap...\n");
	// If the overlapped length is above some threshold
	ol = ol_len >= MIN_OL && ol_len < (opt->rl + MIN_OL);
	if (ol || reason_gap < MIN_OL) {
		ol_len = (ol) ? ol_len : (reason_gap);
		trun_seq(right_eg->contig, ol_len);
		// Adjust the shift value of reads after truncating.
		adj_shift(right_eg, ol_len);
		right_eg->len -= ol_len;
		if (ori)
			merge_eg_to_right(left_eg, right_eg, 0);
		else
			merge_eg_to_left(left_eg, right_eg, 0);
		show_debug_msg(__func__, "Overlapped left_eg and right_eg: %d\n",
				ol_len);
		return;
	}

	left_p_reads = get_parents_reads(left_eg, 0);
	right_p_reads = get_parents_reads(right_eg, 1);
	// p_readarray(left_p_reads, 1);
	paired_reads = get_paired_reads(left_p_reads, right_p_reads, ht->seqs, 0);
	// If no paired reads spanning the two edges, just return.
	if (paired_reads->len == 0) {
		g_ptr_array_free(paired_reads, TRUE);
		show_debug_msg(__func__, "No reads spanning left_eg and right_eg \n");
		return;
	}
	g_ptr_array_free(paired_reads, TRUE);
	g_ptr_array_free(left_p_reads, TRUE);
	g_ptr_array_free(right_p_reads, TRUE);

	cur_pool = new_pool();
	ass_mate_pool = get_mate_pool(left_eg, ht, 0, 0, 1);
	mates = ass_mate_pool->reads;
	left_pool = get_mate_pool(right_eg, ht, 1, 0, 1);
	g_ptr_array_concat(mates, left_pool->reads);
	// p_readarray(mates, 0);

	// Trying to extend from the right_eg to the left.
	if (right_eg->len > opt->ol) {
		ori_len = right_eg->len;
		query = new_seq(right_eg->contig, opt->ol, 0);
		rev_reads_pos(right_eg); // To ensure the shift value is correct
		fill_in_pool(cur_pool, mates, query, 1);
		// Extend by mates in the pool
		while (1) {
			//p_query(query);
			if (is_repetitive_q(query))
				break;
			olpped = seq_ol(left_eg->contig, right_eg->contig, MIN_OL, opt->nm);
			if (olpped || cur_pool->n == 0)
				break;
			forward_by_ra(right_eg, cur_pool, mates, query, 1);
			//p_pool("Pool in gap: ", cur_pool);
		}
		show_debug_msg(__func__, "Right edge extended from %d to %d\n",
				ori_len, right_eg->len);
		clear_pool(cur_pool);
		rev_reads_pos(right_eg);
	}

	// If overlapped, stop here.
	if (olpped) {
		trun_seq(right_eg->contig, MIN_OL);
		to_merge = 1;
	} else {
		// Trying to extend from the left_eg to the right.
		if (left_eg->len > opt->ol) {
			ori_len = left_eg->len;
			query = new_seq(left_eg->contig, opt->rl - TLR_LEN, (left_eg->len
					- (opt->rl - TLR_LEN)));
			fill_in_pool(cur_pool, mates, query, 0);
			while (1) {
				if (is_repetitive_q(query)) {
					break;
				}
				olpped = seq_ol(left_eg->contig, right_eg->contig, MIN_OL,
						opt->nm);
				if (olpped || cur_pool->n == 0)
					break;
				forward_by_ra(left_eg, cur_pool, mates, query, 0);
			} // end of extension
			show_debug_msg(__func__, "Left edge extended from %d to %d\n",
					ori_len, left_eg->len);
			if (olpped) {
				trun_seq(right_eg->contig, MIN_OL);
				to_merge = 1;
			} else {
				gap = est_gap(left_eg, right_eg, ht);
				if (gap != INVALID) {
					to_merge = 1;
				}
				if ((gap - reason_gap) >= TLR_LEN)
					gap = reason_gap;
			}
		}
	}
	if (to_merge) {
		show_debug_msg(__func__, "Gap size: %d \n", gap);
		if (ori) {
			merge_eg_to_right(left_eg, right_eg, gap);
			added_gap = init_gap(right_eg->len, gap, ori);
			g_ptr_array_add(right_eg->gaps, added_gap);
		} else {
			merge_eg_to_left(left_eg, right_eg, gap);
			added_gap = init_gap(left_eg->len, gap, ori);
			g_ptr_array_add(left_eg->gaps, added_gap);
		}
	}
	free_pool(cur_pool);
	free_pool(ass_mate_pool);
	free_pool(left_pool);
	bwa_free_read_seq(1, query);
}

/**
 */
bwa_seq_t *get_init_q(edge *ass_eg, bwa_seq_t *init_q, const int ori) {
	bwa_seq_t *query = 0;
	if (ass_eg->len >= opt->ol) {
		if (ori)
			query = new_seq(ass_eg->contig, opt->ol, 0);
		else
			query = new_seq(ass_eg->contig, opt->ol, ass_eg->len - opt->ol);
	} else {
		// If it is a branch, the edge's length is 1, and extend from the initial query
		if (ori)
			query = new_seq(init_q, opt->ol, 0);
		else
			query = new_seq(init_q, opt->ol, init_q->len - opt->ol);
	}
	return query;
}

bwa_seq_t *ext_by_mates(edge *ass_eg, const hash_table *ht, pool *cur_pool,
		bwa_seq_t *init_q, const int ori) {
	pool *mate_pool;
	bwa_seq_t *query;
	int ori_len = 0;

	// @TODO: if ori is 1, the shift value of reads is messy.
	ori_len = ass_eg->len;
	query = get_init_q(ass_eg, init_q, ori);
	mate_pool = get_mate_pool(ass_eg, ht, ori, 0, 0);
	p_query(__func__, query);
	//p_query(__func__, ass_eg->contig);
	//p_pool("Mate pool: ", mate_pool, 0);
	//p_readarray(mate_pool->reads, 1);
	fill_in_pool(cur_pool, mate_pool->reads, query, ori);
	// Only extend by mates for a small length
	//p_pool("In sub path: ", cur_pool, 0);
	while (ass_eg->len < (opt->rl + 2 + ori_len)) {
		//p_query(__func__, query);

		//show_debug_msg(__func__, "Extended by %d \n", (ass_eg->len - ori_len));
		// p_pool("In sub path: ", cur_pool, 0, 0);
		if (is_repetitive_q(query) || cur_pool->n == 0) {
			break;
		}
		forward_by_ra(ass_eg, cur_pool, mate_pool->reads, query, ori);
		//p_ctg_seq(__func__, ass_eg->contig);
		//show_debug_msg(__func__, "[%d, %d] \n", ass_eg->id, ass_eg->len);
		//p_query(__func__, ass_eg->contig);
	}
	show_debug_msg(__func__, "Edge %d extended by mates from %d to %d \n",
			ass_eg->id, ori_len, ass_eg->len);
	if (ass_eg->len > ori_len) {
		upd_reads(ass_eg, opt->nm);
	}
	//free_pool(cur_pool);
	//p_query(__func__, ass_eg->contig);
	return query;
}

ext_msg *single_ext(edge *ass_eg, pool *c_pool, bwa_seq_t *init_q,
		const hash_table *ht, const int ori) {
	bwa_seq_t *query, *contig, *used = 0;
	alignarray *aligns;
	pool *cur_pool = 0, *mate_pool = new_pool();
	int *next = (int*) calloc(5, sizeof(int));
	int *c = (int*) calloc(5, sizeof(int));
	ext_msg *m = new_msg();

	printf("\n");
	show_debug_msg(__func__,
			"******************* Starting of Single Extension ******************\n");
	cur_pool = c_pool ? c_pool : new_pool();
	contig = ass_eg->contig;
	ass_eg->len = contig->len;
	aligns = g_ptr_array_sized_new(N_DEFAULT_ALIGNS);
	// If we are going to extend from an initial read or from some existing edge
	query = get_init_q(ass_eg, init_q, ori);
	if (ori) {
		seq_reverse(ass_eg->len, ass_eg->contig->seq, 0);
		rev_reads_pos(ass_eg);
	}
	p_query(__func__, query);
	while (1) {
		if (ass_eg->len % 50 == 0) {
			show_debug_msg(__func__, "Assembling... [%d, %d] \n", ass_eg->id,
					ass_eg->len);
		}
		reset_c(next, c); // Reset the counter
		// show_debug_msg(__func__, "Current edge: [%d, %d] \n", ass_eg->id,
		// 		ass_eg->len);
		//p_query(__func__, query);
		if (is_repetitive_q(query)) {
			show_debug_msg(__func__, "[%d, %d] Repetitive pattern, stop!\n",
					ass_eg->id, ass_eg->len);
			m->type = REP_QUE;
			p_query(__func__, query);
			break;
		}
		pe_aln_query(query, query->seq, ht, opt->nm, opt->ol, 0, aligns);
		if (opt->nsp) {
			pe_aln_query(query, query->rseq, ht, opt->nm, opt->ol, 1, aligns);
		}
		//p_align(aligns);
		// Extend the contig, update the counter and sequence pool
		upd_cur_pool(aligns, next, cur_pool, mate_pool, query, ht, ass_eg, ori);
		reset_alg(aligns);
		// p_pool("Current pool: ", cur_pool, next);
		c = get_most(next);
		// If cannot extend, or multiple path, just stop here
		if (c[0] == INVALID_CHAR) {
			m->type = NOT_EXTEND;
			show_debug_msg(__func__, "Not extended any more: [%d, %d] \n",
					ass_eg->id, ass_eg->len);
			break;
		}
		if (c[1] != INVALID_CHAR) {
			show_debug_msg(
					__func__,
					"Multiple branching: a:c:g:t:n => %d:%d:%d:%d:%d \t %d:%d:%d:%d:%d [%d, %d]\n",
					next[0], next[1], next[2], next[3], next[4], c[0], c[1],
					c[2], c[3], c[4], ass_eg->id, ass_eg->len);
			m->type = MUL_PATH;
			m->data = c;
			//p_query("Query now: ", query);
			//p_pool("Pool when multi-branching: ", cur_pool, ori, next);
			break;
		}
		// Only consider one branch. If some reads are used, stop
		used = forward(cur_pool, c[0], ass_eg, ori); // The read used.
		if (used && ass_eg->len >= MIN_LEN_BF_CHECK) {
			show_debug_msg(
					__func__,
					"The read has been used before: %s used in %d at %d [%d, %d]\n",
					used->name, used->contig_id, used->shift, ass_eg->id,
					ass_eg->len);
			//tmp_eg = g_ptr_array_index(all_edges, used->contig_id);
			//p_query(__func__, query);
			//p_query("Used reads", used);
			//p_ctg_seq("Contig now", ass_eg->contig);
			if (used->contig_id == ass_eg->id) {
				m->type = REP_EXTEND;
			} else
				m->type = QUE_PFD;
			m->data = used;
			break;
		}
		ext_con(ass_eg->contig, c[0], 0);
		ass_eg->len = ass_eg->contig->len;
		ext_que(query, c[0], ori);
	}
	m->query = new_seq(query, query->len, 0);
	if (ori) {
		seq_reverse(ass_eg->len, ass_eg->contig->seq, 0);
		rev_reads_pos(ass_eg);
	}
	free(next);
	if (m->type != MUL_PATH)
		free(c);
	free_alg(aligns);
	free_pool(cur_pool);
	free_pool(mate_pool);
	bwa_free_read_seq(1, query);
	show_debug_msg(__func__,
			"******************* Ending of Single Extension ******************\n\n");
	return m;
}

int order_edges(edge *left_eg, edge *right_eg) {

}

int linear_ext(edge *ass_eg, const hash_table *ht, bwa_seq_t *cur_query,
		const int type, const int ori) {
	bwa_seq_t *mate = 0, *query = 0;
	int ori_len = 0, extended = 0, len_init = 0, len_le = 0, len_re = 0,
			reason_gap = 0, max_try_times = 4;
	edge *m_eg;
	ext_msg *m = 0, *m2 = 0;
	pool *m_pool = 0, *c_pool = 0;
	ori_len = ass_eg->len;
	eg_gap *exi_gap = 0;
	if (ass_eg->len >= opt->ol && (type != REP_EXTEND && type != QUE_PFD)) {
		query = get_init_q(ass_eg, 0, ori);
		c_pool = new_pool();
		ext_by_mates(ass_eg, ht, c_pool, query, ori);
		free_pool(c_pool);
		if (ass_eg->len > ori_len)
			return 1;
	}

	show_debug_msg(__func__, "Trying to extend [%d, %d] from the mates.\n",
			ass_eg->id, ass_eg->len);
	while (max_try_times-- > 0 && len_le == len_init) {
		m_pool = get_mate_pool(ass_eg, ht, ori, 1, 0);
		mate = get_query_ol(ass_eg, ht->seqs, m_pool, ori);
		free_pool(m_pool);
		if (!mate) {
			show_debug_msg(__func__, "Mates are empty, return \n");
			return 0;
		}

		// Have to create a new copy of the query, not to affect to following extension
		m_eg = new_eg();
		m_eg->contig = mate->rev_com ? new_mem_rev_seq(mate, opt->rl, 0)
				: new_seq(mate, opt->rl, 0);
		m_eg->id = ass_eg->id;
		m_eg->len = m_eg->contig->len;
		len_re = len_le = len_init = m_eg->len;

		show_debug_msg(__func__, "Mate query, %d times [%p]: \n",
				max_try_times, mate);
		p_query(__func__, mate);
		c_pool = get_init_pool(ht, mate, 0);
		//p_pool("Initial pool: ", c_pool, 0);
		m = single_ext(m_eg, c_pool, 0, ht, 0);

		len_re = m_eg->len;
		c_pool = get_init_pool(ht, mate, 1);
		//p_pool("Initial pool: ", c_pool, 1);
		m2 = single_ext(m_eg, c_pool, 0, ht, 1);
		len_le = m_eg->len;

		mate->used = 1;
		mate->shift = len_le - len_re;
		mate->contig_id = m_eg->id;
		g_ptr_array_add(m_eg->reads, mate);
	}

	show_debug_msg(__func__, "Lengths: %d->%d->%d\n", len_init, len_re, len_le);
	reason_gap = ori ? (opt->mean - opt->rl - (len_re - len_init)) : (opt->mean
			- opt->rl - (len_le - len_re));
	show_debug_msg(__func__, "Reasonable gap: %d\n", reason_gap);

	reason_gap = reason_gap > 0 ? reason_gap : 0;

	if (m_eg->len > opt->rl) {
		// If the edge length is not long enough to close the gap, just merge them
		if (ass_eg->len * 2 < opt->mean && m_eg->len * 2 < opt->mean) {
			if (ori)
				merge_eg_to_right(m_eg, ass_eg, 0);
			else
				merge_eg_to_left(ass_eg, m_eg, 0);
			extended = 1;
		}
		// If not extended, try to merge the edges
		if (!extended) {
			show_debug_msg(__func__, "Reasonable Gap: %d \n", reason_gap);
			exi_gap = find_hole(ass_eg, m_eg, ori);
			if (exi_gap) {
				// If should be concatenated directly
				if (exi_gap->s_index == -1) {
					if (ori)
						fill_in_gap(m_eg, ass_eg, reason_gap, ht, ori);
					else
						fill_in_gap(ass_eg, m_eg, reason_gap, ht, ori);
				} else {
					fill_in_hole(ass_eg, m_eg, ori, exi_gap);
				}
			}
		}

		if (ass_eg->len > ori_len) {
			extended = 1;
			show_debug_msg(__func__, "Edge %d extended from %d to %d \n",
					ass_eg->id, ori_len, ass_eg->len);
		}
	}
	upd_reads(ass_eg, opt->nm);
	destroy_eg(m_eg);
	free_msg(m);
	free_msg(m2);
	return extended;
}

/**
 * Construct the roadmap for a RNA-seq read recursively
 */
edge *pe_ass_edge(edge *parent, edge *cur_eg, pool *c_pool,
		bwa_seq_t *init_query, const hash_table *ht, int level, int ori) {
	bwa_seq_t *query, *contig, *used, *sub_query;
	int c_index = 0, no_sub_path = 0;
	int *c = 0, extended = 0;
	edge *ass_eg, *tmp_eg;
	ext_msg *msg = 0;

	if (cur_eg) {
		ass_eg = cur_eg;
	} else {
		ass_eg = new_eg();
		ass_eg->id = contig_id++;
		ass_eg->ori = ori;
		g_ptr_array_add(all_edges, ass_eg);
		if (parent) { // If is not a root edge
			if (ori) {
				contig = new_seq(init_query, 1, 0);
				g_ptr_array_add(parent->in_egs, ass_eg);
				g_ptr_array_add(ass_eg->out_egs, parent);
			} else {
				contig = new_seq(init_query, 1, init_query->len - 1);
				g_ptr_array_add(parent->out_egs, ass_eg);
				g_ptr_array_add(ass_eg->in_egs, parent);
			}
		} else {
			contig = new_seq(init_query, init_query->len, 0);
		}
		init_query->used = 1;
		init_query->contig_id = ass_eg->id;
		init_query->shift = 0;
		ass_eg->contig = contig;
		ass_eg->len = contig->len;
	}
	query = new_seq(init_query, init_query->len, 0);

	// If there is one option to extend, just forward; otherwise, stop
	while (1) {
		free_msg(msg);
		if (parent) {
			show_debug_msg(__func__,
					"[%d, %d] In sub path: extending by mates... \n",
					ass_eg->id, ass_eg->len);
			// p_readarray(p_er, 0);
			c_pool = new_pool(); // be reused by single_ext
			query = ext_by_mates(ass_eg, ht, c_pool, query, ori);
			p_query(__func__, query);
			parent = 0;
		}
		msg = single_ext(ass_eg, c_pool, query, ht, ori); // c_pool freed here.
		bwa_free_read_seq(1, query);
		query = 0;
		c_pool = 0; // Content has been free in function single_extend.
		if (msg->type == REP_QUE || msg->type == NOT_EXTEND || msg->type
				== REP_EXTEND) {
			if (!opt->pair)
				break;
			extended = linear_ext(ass_eg, ht, msg->query, msg->type, ori);
			if (extended) {
				continue;
			} else
				break;
		} else if (msg->type == QUE_PFD) {
			used = (bwa_seq_t*) msg->data;
			// If the ori is reverse, just not set the right_ctg value.
			if (used->contig_id != ass_eg->id && !ori) {
				// Set the shifted pos, such that it is able to recover the original transcript
				tmp_eg = (edge*) g_ptr_array_index(all_edges, used->contig_id);
				ass_eg->r_shift = used->shift + used->cursor - 1;
				ass_eg->right_ctg = tmp_eg;
				ass_eg->out_egs = tmp_eg->out_egs;
				g_ptr_array_add(tmp_eg->in_egs, ass_eg);
			}
			break;
		} else if (msg->type == MUL_PATH) {
			if (level > MAX_BRANCH_LEVEL) {
				show_debug_msg(__func__, "%d Branching too deep, STOP!!! \n",
						ass_eg->id);
				break;
			}
			no_sub_path = 1;
			c = (int*) msg->data;
			c_index = 0;
			while (c[c_index] != INVALID_CHAR) {
				sub_query = new_seq(msg->query, opt->ol, msg->query->len
						- opt->ol);
				ext_que(sub_query, c[c_index], ori);
				show_debug_msg(__func__,
						"[%d, %d] Multi-braching %d, level %d\n", ass_eg->id,
						ass_eg->len, c_index, level);
				if (vld_ext(ass_eg, sub_query, ht, ori)) {
					show_debug_msg(__func__, "Subpath is feasible to go \n");
					tmp_eg = pe_ass_edge(ass_eg, 0, 0, sub_query, ht,
							level + 1, ori);
					no_sub_path = 0;
				}
				c_index++;
			}
			free(c);
			if (no_sub_path && opt->pair) {
				extended = linear_ext(ass_eg, ht, msg->query, msg->type, ori);
				if (extended)
					continue;
				else
					break;
			} else
				break;
		} else {
			break;
		}
	}
	free_msg(msg);
	log_edge(ass_eg);
	return ass_eg;
}

edge *pe_ass_ctg(roadmap *rm, bwa_seq_t *read, hash_table *ht) {
	edge *eg_i, *cur_eg;
	pool *c_pool;
	int i = 0, r_ext_len = 0;
	readarray *reads;
	bwa_seq_t *r;

	p_query(__func__, read);
	c_pool = get_init_pool(ht, read, 0);
	// p_pool("Initial Pool: ", c_pool);
	cur_eg = pe_ass_edge(0, 0, c_pool, read, ht, 0, 0);
	reads = cur_eg->reads;
	r_ext_len = cur_eg->len;
	// reset_pool_info(ht->seqs, ht->n_seqs);
	// The c_pool is freed after extension
	c_pool = get_init_pool(ht, read, 1);
	//p_pool("Initial Pool: ", c_pool);
	show_debug_msg(__func__,
			"Extending to the left ***************************** [%d, %d]\n",
			cur_eg->id, cur_eg->len);
	pe_ass_edge(0, cur_eg, c_pool, read, ht, 0, 1);
	// If the read is not extendable;
	// or it connects to some existing contigs right away, just remove it
	if (cur_eg->len <= opt->rl + 4 && (cur_eg->right_ctg
			|| (cur_eg->out_egs->len == 0 && cur_eg->in_egs->len == 0))) {
		if (cur_eg->right_ctg) {
			for (i = 0; i < cur_eg->in_egs->len; i++) {
				eg_i = g_ptr_array_index(cur_eg->in_egs, i);
				if (eg_i == cur_eg) {
					g_ptr_array_remove_index_fast(cur_eg->in_egs, i);
					break;
				}
			}
		}
		cur_eg->right_ctg = 0;
		cur_eg->alive = 0;
		return cur_eg;
	}
	for (i = 0; i < reads->len; i++) {
		r = g_ptr_array_index(reads, i);
		r->shift += (cur_eg->len - r_ext_len);
	}

	for (i = 0; i < all_edges->len; i++) {
		eg_i = g_ptr_array_index(all_edges, i);
		if (eg_i->in_egs->len == 0 && eg_i->alive) {
			g_ptr_array_add(rm->start_egs, eg_i);
			rm->start_eg_n++;
			eg_i->is_root = 1;
		}
	}
	return cur_eg;
}

void pe_ass_core(const char *starting_reads, const char *fa_fn,
		const char *pet_fn) {
	int counter = 0, index = 0, s_index = 0, e_index = 0, pre_ctg_id = 0, i = 0;
	bwa_seq_t *p; // sequence of RNA-SEQs, RNA-PETs and temp
	char *h = malloc(BUFSIZE), *msg = calloc(BUFSIZE, sizeof(char));
	char line[80];
	FILE *ass_fa = xopen("read/ass_tx.fa", "w");
	FILE *ass_contigs = xopen("read/ass_contigs.fa", "w");
	FILE *all_contigs = xopen("read/contigs.fa", "w");
	FILE *start_reads = xopen("read/start_reads.txt", "w");
	FILE *solid_reads = xopen(starting_reads, "r");

	hash_table *ht;
	edge *eg, *eg_i;
	clock_t t = clock();
	float t_eclipsed = 0;
	all_edges = g_ptr_array_sized_new(INIT_EDGE_SPACE);

	ht = pe_load_hash(fa_fn);
	left_rm = new_rm();

	s_index = 0;
	e_index = 1;
	while (ht->n_seqs * STOP_THRE > n_reads_consumed) {
		if (fgets(line, 80, solid_reads) != NULL && counter <= 10000)
			index = atoi(line);
		else
			index = (int) (rand_f() * ht->n_seqs);
		if (counter < s_index) {
			counter++;
			continue;
		}
		if (counter >= e_index)
			break;
		t_eclipsed = (float) (clock() - t) / CLOCKS_PER_SEC;
		p = &ht->seqs[5226253];
		//p = &ht->seqs[index];
		if (p->used) {
			show_msg(__func__, "Read used: %s\n", p->name);
			continue;
		}
		if (has_rep_pattern(p)) {
			show_debug_msg(__func__, "Read has repeat pattern, skip.\n");
			p_query(__func__, p);
			continue;
		}
		//counter++;
		show_msg(__func__, "Reads Consumed: %d\n", n_reads_consumed);
		sprintf(msg, "Processing read %d: %s... \n", counter, p->name);
		show_msg(__func__, msg);
		show_debug_msg(__func__, msg);
		fputs(p->name, start_reads);
		fputs("\n", start_reads);
		eg = pe_ass_ctg(left_rm, p, ht);
		p->used = 1;
		if (eg->len == p->len) {
			destroy_eg(eg);
			contig_id--;
			n_reads_consumed++;
		} else {
			counter++;
			for (i = pre_ctg_id; i < all_edges->len; i++) {
				eg_i = g_ptr_array_index(all_edges, i);
				if (eg_i->alive)
					n_reads_consumed += eg_i->reads->len;
				//log_edge(eg_i);
			}
		}
		pre_ctg_id = all_edges->len;
		p->contig_id = eg->id;

		sprintf(msg, "-------------------------------------- %.2f sec \n",
				t_eclipsed);
		show_msg(__func__, msg);
		show_debug_msg(__func__, msg);
		//post_pro(left_rm, all_edges, contig_id, opt);
	} // All RNA-PETs are read.

	fprintf(stderr,
			"[pe_ass_core] ------------------------------------------------------ \n");

	// Post process the roadmaps.
	log_reads(all_edges);
	show_msg(__func__, "Post processing the roadmap... \n");
	graph_by_edges(all_edges, "graph/rm_bf_update.dot");
	save_edges(all_edges, all_contigs, 0, 1, opt->rl * 1.5);
	post_pro(left_rm, all_edges, contig_id, opt);
	graph_by_edges(all_edges, "graph/rm_after_update.dot");
	save_edges(all_edges, ass_contigs, 0, 0, opt->rl * 1.5);

	free(h);
	fclose(solid_reads);
	fclose(ass_fa);
	fclose(ass_contigs);
	fclose(all_contigs);
	fclose(start_reads);
}

int ass_usage() {
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:   peta ass [options] <lib.fa> <rnapet.fa>\n\n");
	return 1;
}

int pe_ass(int argc, char *argv[]) {
	int c;
	clock_t t = clock();
	opt = init_opt();

	while ((c = getopt(argc, argv, "n:o:r:a:m:s:p:d:b:")) >= 0) {
		switch (c) {
		case 'n':
			opt->nm = atoi(optarg);
			break;
		case 'o':
			opt->ol = atoi(optarg);
			break;
		case 'r':
			opt->rl = atoi(optarg);
			break;
		case 'a':
			opt->n_alg = atoi(optarg);
			break;
		case 'm':
			opt->mean = atoi(optarg);
			break;
		case 'd':
			opt->sd = atoi(optarg);
			break;
		case 's':
			opt->nsp = !atoi(optarg);
			break;
		case 'p':
			opt->pair = atoi(optarg);
			break;
		case 'b':
			opt->solid_reads = strdup(optarg);
			break;
		default:
			return 1;
		}
	}
	if (optind + 2 > argc) {
		return ass_usage();
	}

	fprintf(stderr, "%s \n", argv[optind]);
	fprintf(stderr, "%s \n", argv[optind + 1]);
	fprintf(stderr, "nm = %d \n", opt->nm);
	fprintf(stderr, "ol = %d \n", opt->ol);
	fprintf(stderr, "rl = %d \n", opt->rl);
	fprintf(stderr, "n_alg = %d \n", opt->n_alg);
	fprintf(stderr, "m = %d \n", opt->mean);
	fprintf(stderr, "sd = %d \n", opt->sd);
	fprintf(stderr, "nsp = %d \n", opt->nsp);

	pe_ass_core(opt->solid_reads, argv[optind], argv[optind + 1]);

	free(opt);

	fprintf(stderr, "[pe_ass] Done: %.2f sec\n", (float) (clock() - t)
			/ CLOCKS_PER_SEC);
	return 0;
}
