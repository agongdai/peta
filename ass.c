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
#include "pealn.h"
#include "roadmap.h"
#include "edge.h"
#include "edgelist.h"
#include "readrm.h"

/**
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
	m->used = NULL;
	m->type = 0;
	m->message = NULL;
	m->query = NULL;
	m->counter = NULL;
	return m;
}

void free_msg(ext_msg *m) {
	if (m) {
		free(m->message);
		free(m->counter);
		bwa_free_read_seq(1, m->query);
		bwa_free_read_seq(1, m->used);
		free(m);
		m = 0;
	}
}

int apply_new_ctg_id() {
	edge *tmp_eg = NULL;
	int next_id = 0, i = 0, largest_id = 0;
	if (!all_edges)
		all_edges = g_ptr_array_sized_new(INIT_EDGE_SPACE);
	if (all_edges->len == 0)
		return 0;
	for (i = 0; i < all_edges->len; i++) {
		tmp_eg = g_ptr_array_index(all_edges, i);
		p_flat_eg(tmp_eg);
		if (tmp_eg->id > largest_id)
			largest_id = tmp_eg->id;
	}
	next_id = largest_id + 1;
	show_debug_msg("NEW_CONTIG", "Next id = %d \n", next_id);
	return next_id;
}

/**
 * Checks whether the position right before the cursor are the same.
 * contig: aaaaaaaac
 * s:      aaaaaaaatggggg
 * The position 'c' is not the same as 't', return False
 *
 * contig: aaaaaaaacccccccccccccg
 * s:      ttttttttccccccccccccctggggggg
 * Starting positions of 's' not the same as contig, return False
 */
int check_ign(const int ori, bwa_seq_t *s, bwa_seq_t *contig) {
	int ignore = 0, to_check_overlap = 0;
	bwa_seq_t *tmp = NULL;
	to_check_overlap = contig->len >= s->len / 4;
	if (ori) {
		seq_reverse(contig->len, contig->seq, 0); // Because single_extension reverse the sequence first
		if (s->rev_com) {
			if ((s->rseq[s->cursor + 1] != contig->seq[contig->len - 1]))
				ignore = 1;
			tmp = new_mem_rev_seq(s, s->len, 0);
			if (to_check_overlap && find_ol(tmp, contig, opt->nm) < s->len / 4) {
				ignore = 1;
			}
			bwa_free_read_seq(1, tmp);
		} else {
			if (s->seq[s->cursor + 1] != contig->seq[contig->len - 1])
				ignore = 1;
			if (to_check_overlap && find_ol(s, contig, opt->nm) < s->len / 4) {
				ignore = 1;
			}
		}
		seq_reverse(contig->len, contig->seq, 0);
	} else {
		if (s->rev_com) {
			if ((s->rseq[s->cursor - 1] != contig->seq[contig->len - 1]))
				ignore = 1;
			tmp = new_mem_rev_seq(s, s->len, 0);
			if (to_check_overlap && find_ol(contig, tmp, opt->nm) < s->len / 4)
				ignore = 1;
			bwa_free_read_seq(1, tmp);
		} else {
			if (s->seq[s->cursor - 1] != contig->seq[contig->len - 1])
				ignore = 1;
			if (to_check_overlap && find_ol(contig, s, opt->nm) < s->len / 4)
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
	index64 index;
	alg *a;
	bwa_seq_t *s = NULL, *seqs, *contig = ass_eg->contig, *mate;
	seqs = ht->seqs;
	for (i = 0; i < alns->len; i++) {
		a = g_ptr_array_index(alns, i);
		index = a->r_id;
		if (index >= ht->n_seqs)
			continue;
		s = &seqs[index];
		mate = get_mate(s, seqs);
		// If the read is used already, but the orientation is different, ignore it.
		if (s->is_in_c_pool || (s->used && (s->rev_com != a->rev_comp))
				|| (mate->used && (mate->rev_com != a->rev_comp)))
			continue;
		s->rev_com = a->rev_comp;
		mate->rev_com = s->rev_com;
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
		pool_add(cur_pool, s);
		if (mate_pool && !mate->is_in_m_pool && !mate->is_in_c_pool) {
			mate_pool_add(mate_pool, mate);
		}
	}
	// To ensure that current reads are similar to the contig
	//p_pool("Before", cur_pool, NULL);
	rm_partial(cur_pool, ori, query, opt->nm);
	//p_pool("After", cur_pool, NULL);
	// Add overlapped reads from the mate_pool, the overlapping length is read length / 4.
	if (mate_pool)
		overlap_mate_pool(cur_pool, mate_pool, contig, ori);
	for (i = 0; i < cur_pool->n; i++) {
		s = g_ptr_array_index(cur_pool->reads, i);
		if (s->rev_com)
			check_c(next, s->rseq[s->cursor]);
		else
			check_c(next, s->seq[s->cursor]);
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
	concate_readarray(reads, parent->reads);
	while (pre && acc_len < opt->mean + opt->sd * SD_TIMES) {
		concate_readarray(reads, pre->reads);
		acc_len += pre->len;
		if (ori) {
			if (pre->out_egs->len > 0 && !pre->right_ctg)
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
		g_ptr_array_free(reads, TRUE);
		return NULL;
	}
	return reads;
}

bwa_seq_t *get_query_ol(edge *ass_eg, bwa_seq_t *seqs, pool *m_pool,
		const int ori) {
	bwa_seq_t *read, *mate = NULL;
	//readarray *p_ra;
	int start, start_copy, i = 0;
	int is_correct_ori = 0;
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
		if (ori)
			is_correct_ori = read->rev_com ? is_left_mate(read->name)
					: is_right_mate(read->name);
		else
			is_correct_ori = read->rev_com ? is_right_mate(read->name)
					: is_left_mate(read->name);
		if (is_correct_ori) {
			mate = get_mate(read, seqs);
			if (!mate->used && !has_n(mate) && mate->contig_id != -1)
				return mate;
		}
		start++;
	}
	// From start to the beginning
	start = start_copy;
	while (start >= 0) {
		read = g_ptr_array_index(ass_eg->reads, start);
		if (ori)
			is_correct_ori = read->rev_com ? is_left_mate(read->name)
					: is_right_mate(read->name);
		else
			is_correct_ori = read->rev_com ? is_right_mate(read->name)
					: is_left_mate(read->name);
		if (is_correct_ori) {
			mate = get_mate(read, seqs);
			if (!mate->used && !has_n(mate) && mate->contig_id != -1)
				return mate;
		}
		start--;
	}

	for (i = 0; i < m_pool->reads->len; i++) {
		read = g_ptr_array_index(m_pool->reads, i);
		if (ori)
			is_correct_ori = read->rev_com ? is_left_mate(read->name)
					: is_right_mate(read->name);
		else
			is_correct_ori = read->rev_com ? is_right_mate(read->name)
					: is_left_mate(read->name);
		if (is_correct_ori) {
			if (!read->used && !has_n(read) && mate->contig_id != -1)
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
	//show_msg(__func__, "start = %d, end = %d \n", start, end);
	// Add the mates of the reads in the middle first
	for (i = 0; i < reads->len; i++) {
		s = g_ptr_array_index(reads, i);
		// If not range requirement, just go into the 'if' block;
		// If range is required (within_range = 1), check the shift value as well
		if (!within_range || (s->shift >= start && s->shift < end)) {
			if (ori)
				mate = s->rev_com ? get_right_mate(s, seqs) : get_left_mate(s,
						seqs);
			else
				mate = s->rev_com ? get_left_mate(s, seqs) : get_right_mate(s,
						seqs);
			if (mate && (used || !mate->used))
				pool_uni_add(mate_pool, mate);
		}
	}
	g_ptr_array_free(reads, TRUE);
	//p_pool("Mate pool: ", mate_pool);
	return mate_pool;
}

pool *get_init_pool(const hash_table *ht, bwa_seq_t *init_read, const int ori) {
	alignarray *alns = NULL;
	int i = 0;
	pool *init_pool = NULL;
	bwa_seq_t *seqs = ht->seqs, *s = NULL, *query = init_read, *mate = NULL;
	alg *a;
	if (is_repetitive_q(init_read) || has_rep_pattern(init_read)) {
		return NULL;
	}
	init_pool = new_pool();
	alns = g_ptr_array_sized_new(N_DEFAULT_ALIGNS);
	if (init_read->rev_com)
		query = new_mem_rev_seq(init_read, init_read->len, 0);
	p_query(__func__, query);

	pe_aln_query(query, query->seq, ht, opt->nm + 2, opt->ol, 0, alns);
	if (opt->nsp) {
		pe_aln_query(query, query->rseq, ht, opt->nm + 2, opt->ol, 1, alns);
	}
	for (i = 0; i < alns->len; i++) {
		a = (alg*) g_ptr_array_index(alns, i);
		s = &seqs[a->r_id];
		mate = get_mate(s, seqs);
		s->rev_com = a->rev_comp;
		mate->rev_com = a->rev_comp;
		if (s->rev_com)
			s->cursor = ori ? (0 - a->pos - 1) : (s->len - a->pos);
		else
			s->cursor = ori ? (a->pos - 1) : (a->pos + s->len);
		if (s->contig_id == -1 || s->used || s->cursor < 0 || s->cursor
				>= s->len) {
			s->cursor = 0;
			continue;
		}
		pool_add(init_pool, s);
	}
	free_alg(alns);
	if (init_read->rev_com)
		bwa_free_read_seq(1, query);
	return init_pool;
}

int vld_ext(edge *parent, bwa_seq_t *query, const hash_table *ht, const int ori) {
	bwa_seq_t *seqs = ht->seqs, *s, *mate;
	int i = 0, n_on = 0, acc_len = 0, q_pos_on_mate;
	ubyte_t cursor_char, next_char;
	readarray *reads = NULL;

	reads = get_parents_reads(parent, ori);
	acc_len = get_parents_len(parent, ori);
	if (!reads || acc_len < opt->mean) {
		if (reads)
			g_ptr_array_free(reads, TRUE);
		show_debug_msg(__func__, "Accumulated length too short: %d \n", acc_len);
		return 1;
	}

	//p_readarray(reads, 1);
	cursor_char = ori ? query->seq[0] : query->seq[query->len - 1];
	p_query(__func__, query);
	for (i = 0; i < reads->len; i++) {
		s = g_ptr_array_index(reads, i);
		mate = get_mate(s, seqs);
		if (mate && !mate->used) {
			if (s->rev_com) {
				q_pos_on_mate = is_sub_seq_byte(query->rseq, query->len, 0,
						mate, opt->nm, 0);
				if (q_pos_on_mate != NOT_FOUND) {
					next_char = ori ? mate->rseq[q_pos_on_mate]
							: mate->rseq[q_pos_on_mate + mate->len - 1];
					//					show_msg(__func__,
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
					//					show_msg(__func__,
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
	//	show_msg(__func__, "confirm_c: %d; confirm_c_2: %d \n", confirm_c, confirm_c_2);
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
	for (i = 0; i < reads->len; i++) {
		read = g_ptr_array_index(reads, i);
		if (!read->used)
			add_pool_by_ol(p, query, read, ori);
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
	ea = get_paired_reads(l_par_reads, r_par_reads, ht->seqs);
	g_ptr_array_free(l_par_reads, TRUE);
	g_ptr_array_free(r_par_reads, TRUE);
	if (ea->len == 0) {
		g_ptr_array_free(ea, TRUE);
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
		//show_msg(__func__, "left = %d, right = %d \n", l_index, r_index);
		gap_sum += opt->mean - len;
		p_count++;
	}
	g_ptr_array_free(ea, TRUE);
	//show_msg(__func__, "gan_sum = %d, p_count = %d \n", gap_sum, p_count);
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
	pool *cur_pool, *ass_mate_pool, *left_pool = NULL;
	int olpped = 0, gap = 0, ol_len = 0, to_merge = 0, ori_len = 0;
	int ol = 0;
	bwa_seq_t *query = NULL;
	eg_gap *added_gap = NULL;
	//	p_query(left_eg->contig);
	//	p_query(right_eg->contig);
	//p_ctg_seq(__func__, left_eg->contig);
	//p_ctg_seq(__func__, right_eg->contig);
	ol_len = find_ol(left_eg->contig, right_eg->contig, opt->nm);
	show_debug_msg(__func__, "Filling the gap...\n");

	show_debug_msg(__func__,
			"Checking mates spanning left and right edges...\n");
	left_p_reads = get_parents_reads(left_eg, 0);
	//p_readarray(left_p_reads, 1);
	right_p_reads = get_parents_reads(right_eg, 1);
	//p_readarray(right_p_reads, 1);
	paired_reads = get_paired_reads(left_p_reads, right_p_reads, ht->seqs);
	free_readarray(left_p_reads);
	free_readarray(right_p_reads);
	// p_readarray(paired_reads, 1);
	// If no paired reads spanning the two edges, just return.
	if (paired_reads->len < MIN_VALID_PAIRS) {
		show_debug_msg(__func__,
				"Not enough pairs spanning left_eg and right_eg: %d \n",
				paired_reads->len);
		free_readarray(paired_reads);
		return;
	}
	free_readarray(paired_reads);

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

	cur_pool = new_pool();
	ass_mate_pool = get_mate_pool(left_eg, ht, 0, 0, 1);
	mates = ass_mate_pool->reads;
	left_pool = get_mate_pool(right_eg, ht, 1, 0, 1);
	g_ptr_array_concat(mates, left_pool->reads);
	// p_readarray(mates, 0);

	// Trying to extend from the right_eg to the left.
	show_debug_msg(__func__,
			"Trying to extend from the right_eg to the left...\n");
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
			olpped = seq_ol(left_eg->contig, right_eg->contig, CLOSE_MIN_OL,
					opt->nm);
			if (olpped || cur_pool->n == 0)
				break;
			forward_by_ra(right_eg, cur_pool, mates, query, 1);
			//p_pool("Pool in gap: ", cur_pool);
		}
		show_debug_msg(__func__, "Right edge extended from %d to %d\n",
				ori_len, right_eg->len);
		clear_pool(cur_pool);
		rev_reads_pos(right_eg);
		bwa_free_read_seq(1, query);
	}

	// If overlapped, stop here.
	if (olpped) {
		trun_seq(right_eg->contig, MIN_OL);
		to_merge = 1;
	} else {
		// Trying to extend from the left_eg to the right.
		show_debug_msg(__func__,
				"Trying to extend from the left_eg to the right...\n");
		if (left_eg->len > opt->ol) {
			ori_len = left_eg->len;
			query = new_seq(left_eg->contig, opt->rl - TLR_LEN, (left_eg->len
					- (opt->rl - TLR_LEN)));
			fill_in_pool(cur_pool, mates, query, 0);
			while (1) {
				if (is_repetitive_q(query)) {
					break;
				}
				olpped = seq_ol(left_eg->contig, right_eg->contig,
						CLOSE_MIN_OL, opt->nm);
				if (olpped || cur_pool->n == 0)
					break;
				forward_by_ra(left_eg, cur_pool, mates, query, 0);
			} // end of extension
			show_debug_msg(__func__, "Left edge extended from %d to %d\n",
					ori_len, left_eg->len);
			if (olpped) {
				trun_seq(right_eg->contig, CLOSE_MIN_OL);
				to_merge = 1;
			} else {
				gap = est_gap(left_eg, right_eg, ht);
				if (gap != INVALID) {
					to_merge = 1;
				}
				if ((gap - reason_gap) >= TLR_LEN)
					gap = reason_gap;
			}
			bwa_free_read_seq(1, query);
		}
	}
	if (to_merge) {
		//		p_ctg_seq("Left: ", left_eg->contig);
		//		p_ctg_seq("Right:", right_eg->contig);
		show_debug_msg(__func__, "Gap size: %d \n", gap);
		if (ori) {
			added_gap = init_gap(right_eg->len, gap, ori);
			g_ptr_array_add(right_eg->gaps, added_gap);
			merge_eg_to_right(left_eg, right_eg, gap);
		} else {
			added_gap = init_gap(left_eg->len, gap, ori);
			g_ptr_array_add(left_eg->gaps, added_gap);
			merge_eg_to_left(left_eg, right_eg, gap);
		}
	}
	free_pool(cur_pool);
	free_mate_pool(ass_mate_pool);
	free_pool(left_pool);
}

/**
 */
bwa_seq_t *get_init_q(edge *ass_eg, bwa_seq_t *init_q, const int ori) {
	bwa_seq_t *query = 0;
	if (ass_eg && ass_eg->len >= opt->ol) {
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

bwa_seq_t *ext_by_pool(edge *ass_eg, const hash_table *ht, pool *cur_pool,
		bwa_seq_t *init_q, const int ori) {
	bwa_seq_t *query = NULL;
	alignarray *aligns = NULL;
	int ori_len = 0;
	int *next = (int*) calloc(5, sizeof(int));
	int most_c = 4;

	aligns = g_ptr_array_sized_new(N_DEFAULT_ALIGNS);
	ori_len = ass_eg->len;
	query = get_init_q(ass_eg, init_q, ori);
	p_query(__func__, query);
	pe_aln_query(query, query->seq, ht, opt->nm, opt->ol, 0, aligns);
	if (opt->nsp) {
		pe_aln_query(query, query->rseq, ht, opt->nm, opt->ol, 1, aligns);
	}
	// p_align(aligns);
	// Only extend for maximum 4 bases
	while (ass_eg->len - ori_len < EXT_BY_POOL_SIZE) {
		reset_c(next, NULL); // Reset the counter
		upd_cur_pool(aligns, next, cur_pool, NULL, query, ht, ass_eg, ori);
		// p_pool(__func__, cur_pool, NULL);
		most_c = get_pure_most(next); // No branching at all
		forward(cur_pool, most_c, ass_eg, ori);
		ext_con(ass_eg->contig, most_c, 0);
		ass_eg->len = ass_eg->contig->len;
		ext_que(query, most_c, ori);
		clean_cur_pool(cur_pool);
		reset_alg(aligns);
		p_ctg_seq(__func__, ass_eg->contig);
	}
	show_debug_msg(__func__, "Edge %d extended by pool reads from %d to %d \n",
			ass_eg->id, ori_len, ass_eg->len);
	if (ass_eg->len > ori_len) {
		upd_reads(ass_eg, opt->nm);
	}
	bwa_free_read_seq(1, init_q);
	free(next);
	free_alg(aligns);
	return query;
}

int *vld_branchs(edge *ass_eg, const int *c, const hash_table *ht,
		const int ori, bwa_seq_t *query) {
	int c_index = 0, valid_count = 0;
	bwa_seq_t *sub_query = NULL;
	int *valid_c = (int*) calloc(5, sizeof(int));
	for (c_index = 0; c_index < 5; c_index++) {
		valid_c[c_index] = INVALID_CHAR;
	}
	c_index = 0;
	while (c[c_index] != INVALID_CHAR) {
		sub_query = new_seq(query, opt->ol, query->len - opt->ol);
		ext_que(sub_query, c[c_index], ori);
		if (vld_ext(ass_eg, sub_query, ht, ori)) {
			valid_c[valid_count++] = c[c_index];
		}
		bwa_free_read_seq(1, sub_query);
		c_index++;
	}
	return valid_c;
}

ext_msg *single_ext(edge *ass_eg, pool *c_pool, bwa_seq_t *init_q,
		const hash_table *ht, const int ori) {
	bwa_seq_t *query, *contig, *used = NULL;
	alignarray *aligns = NULL;
	pool *cur_pool = NULL, *mate_pool = NULL;
	int *next = (int*) calloc(5, sizeof(int));
	int *c = (int*) calloc(5, sizeof(int)), *valid_c = NULL;
	ext_msg *m = new_msg();

	printf("\n");
	show_debug_msg(__func__,
			"******************* Starting of Single Extension ******************\n");
	cur_pool = c_pool ? c_pool : new_pool();
	mate_pool = new_pool();
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
	//p_pool("Init Current pool: ", cur_pool, next);
	while (1) {
		if (ass_eg->len % 50 == 0) {
			show_debug_msg(__func__, "Assembling... [%d, %d] \n", ass_eg->id,
					ass_eg->len);
			clean_mate_pool(mate_pool);
		}
		reset_c(next, c); // Reset the counter
		// show_msg(__func__, "Current edge: [%d, %d] \n", ass_eg->id,
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
		// p_align(aligns);
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
			if (check_next_cursor(cur_pool, c, ori)) { // Check one base ahead the cursor
				valid_c = vld_branchs(ass_eg, c, ht, ori, query); // Validate the branch by pairs
				if (valid_c[0] == INVALID_CHAR) { // If there is no valid branch
					m->type = NOT_EXTEND;
					free(valid_c);
					valid_c = NULL;
					break;
				}
				if (valid_c[1] != INVALID_CHAR) {
					show_debug_msg(
							__func__,
							"Multiple branching: a:c:g:t:n => %d:%d:%d:%d:%d \t %d:%d:%d:%d:%d [%d, %d]\n",
							next[0], next[1], next[2], next[3], next[4], c[0],
							c[1], c[2], c[3], c[4], ass_eg->id, ass_eg->len);
					m->type = MUL_PATH;
					m->counter = c;
					// p_query("Query now: ", query);
					// p_pool("Pool when multi-branching: ", cur_pool, next);
					free(valid_c);
					valid_c = NULL;
					break;
				}
				c[0] = valid_c[0];
			}
		}
		// Only consider one branch. If some reads are used, stop
		used = forward(cur_pool, c[0], ass_eg, ori); // The read used.
		if (used && ass_eg->len >= MIN_LEN_BF_CHECK) {
			show_debug_msg(
					__func__,
					"The read has been used before: %s used in %d at %d [%d, %d]\n",
					used->name, used->contig_id, used->shift, ass_eg->id,
					ass_eg->len);
			p_query("USED", used);
			//			p_pool("Current pool", cur_pool, next);
			if (used->contig_id == ass_eg->id) {
				m->type = REP_EXTEND;
			} else
				m->type = QUE_PFD;
			m->used = new_seq(used, used->len, 0);
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
	free(valid_c);
	free_alg(aligns);
	free_pool(cur_pool);
	free_mate_pool(mate_pool);
	bwa_free_read_seq(1, query);
	show_debug_msg(__func__,
			"******************* Ending of Single Extension ******************\n\n");
	return m;
}

int linear_ext(edge *ass_eg, const hash_table *ht, bwa_seq_t *cur_query,
		const int type, const int ori) {
	bwa_seq_t *mate = NULL, *tmp = NULL;
	int ori_len = 0, extended = 0, len_init = 0, len_le = 0, len_re = 0,
			reason_gap = 0, max_try_times = 4;
	edge *m_eg = NULL;
	ext_msg *m = NULL, *m2 = NULL;
	pool *m_pool = NULL, *c_pool = NULL;
	eg_gap *exi_gap = NULL;

	ori_len = ass_eg->len;
	show_debug_msg(__func__, "Trying to extend [%d, %d] from the mates.\n",
			ass_eg->id, ass_eg->len);
	while (max_try_times-- > 0 && len_le == len_init) {
		destroy_eg(m_eg); // If tried multiple times, the m_eg must be destroyed.
		m_pool = get_mate_pool(ass_eg, ht, ori, 1, 0);
		mate = get_query_ol(ass_eg, ht->seqs, m_pool, ori);
		free_mate_pool(m_pool);
		if (!mate) {
			show_debug_msg(__func__, "Mates are empty, return \n");
			return 0;
		}

		// Have to create a new copy of the query, not to affect to following extension
		m_eg = new_eg();
		m_eg->contig = mate->rev_com ? new_mem_rev_seq(mate, opt->rl, 0)
				: new_seq(mate, opt->rl, 0);
		tmp = get_mate(mate, ht->seqs);
		p_query("USED", tmp);
		p_query("MATE", mate);
		p_ctg_seq("INIT_CTG", m_eg->contig);
		mate->contig_id = -1; // Means used by some hangling template
		m_eg->id = ass_eg->id;
		m_eg->len = m_eg->contig->len;
		len_re = len_le = len_init = m_eg->len;

		show_debug_msg(__func__, "Mate query, %d times [%p]: \n",
				max_try_times, mate);
		p_query(__func__, mate);

		c_pool = get_init_pool(ht, mate, 0);
		if (c_pool) // c_pool is freed in single_ext()
			m = single_ext(m_eg, c_pool, 0, ht, 0);

		len_re = m_eg->len;
		c_pool = get_init_pool(ht, mate, 1);
		if (c_pool) {
			show_debug_msg(__func__, "Linear extension to left");
			m2 = single_ext(m_eg, c_pool, 0, ht, 1);
		}
		len_le = m_eg->len;

		mate->shift = len_le - len_re;
		readarray_add(m_eg, mate);
		free_msg(m);
		free_msg(m2);
		m = NULL;
		m2 = NULL;
	}

	show_debug_msg(__func__, "Lengths: %d->%d->%d\n", len_init, len_re, len_le);
	reason_gap = ori ? (opt->mean - opt->rl - (len_re - len_init)) : (opt->mean
			- opt->rl - (len_le - len_re));

	if (reason_gap < (0 - opt->mean)) {
		show_debug_msg(__func__,
				"Gap size is not reasonable %d, wrong jumping! \n", reason_gap);
		extended = 0;
	} else {
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
				show_debug_msg(__func__, "Reasonable gap: %d \n", reason_gap);
				//			p_readarray(ass_eg->reads, 1);
				//			p_readarray(m_eg->reads, 1);
				exi_gap = find_hole(ass_eg, m_eg, ori);
				if (exi_gap) {
					// If should be concatenated directly
					show_debug_msg(__func__,
							"Existing Gap: %d + %d, ori %d \n",
							exi_gap->s_index, exi_gap->size, exi_gap->ori);
					if (exi_gap->s_index == -1) {
						free_eg_gap(exi_gap);
						if (ori)
							fill_in_gap(m_eg, ass_eg, reason_gap, ht, ori);
						else
							fill_in_gap(ass_eg, m_eg, reason_gap, ht, ori);
					} else {
						fill_in_hole(ass_eg, m_eg, ori, exi_gap, opt->nm,
								opt->rl);
					}
				} else {
					show_debug_msg(__func__,
							"No place reasonable to put [%d, %d]\n", m_eg->id,
							m_eg->len);
				}
			}

			if (ass_eg->len > ori_len) {
				extended = 1;
				show_debug_msg(__func__, "Edge %d extended from %d to %d \n",
						ass_eg->id, ori_len, ass_eg->len);
			}
		}
	}
	if (extended)
		upd_reads(ass_eg, opt->nm);
	destroy_eg(m_eg);
	return extended;
}

int post_vld_branch(edge *ass_eg, edge *child, const hash_table *ht,
		const int ori) {
	readarray *paired_reads = NULL, *main_eg_reads = NULL;
	readarray *branch_eg_reads = NULL;
	int opp_ori = ori ? 0 : 1;
	int valid = 1, main_acc_len = 0, child_acc_len = 0;

	main_eg_reads = get_parents_reads(ass_eg, ori);
	branch_eg_reads = get_parents_reads(child, opp_ori);
	main_acc_len = get_parents_len(ass_eg, ori);
	child_acc_len = get_parents_len(child, opp_ori);
	if (main_acc_len < opt->mean || branch_eg_reads < opt->mean)
			return 1;
	if (ori)
		paired_reads = get_paired_reads(branch_eg_reads, main_eg_reads,
				ht->seqs);
	else
		paired_reads = get_paired_reads(main_eg_reads, branch_eg_reads,
				ht->seqs);
	if (paired_reads->len < MIN_VALID_PAIRS) {
		show_debug_msg(__func__, "Branch [%d, %d] abandoned. \n", child->id,
				child->len);
		free_branch(child, ori, all_edges);
		cut_connection(ass_eg, child, ori);
	} else
		valid = 1;
	free_readarray(main_eg_reads);
	free_readarray(branch_eg_reads);
	free_readarray(paired_reads);
	return valid;
}

edge *init_ctg(edge *parent, edge *cur_eg, bwa_seq_t *init_query, const int ori) {
	edge *ass_eg = NULL;
	bwa_seq_t *contig = NULL;
	if (cur_eg) {
		ass_eg = cur_eg;
	} else {
		ass_eg = new_eg();
		ass_eg->id = apply_new_ctg_id();
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
			readarray_uni_add(ass_eg, init_query); // Only add it for ORIGINAL reads!
		}
		init_query->shift = 0;
		ass_eg->contig = contig;
		ass_eg->len = contig->len;
	}
	return ass_eg;
}

/**
 * Construct the roadmap for a RNA-seq read recursively
 */
edge *pe_ass_edge(edge *parent, edge *cur_eg, pool *c_pool,
		bwa_seq_t *init_query, const hash_table *ht, int level, int ori) {
	bwa_seq_t *query, *used, *sub_query;
	int c_index = 0, no_sub_path = 0;
	int *c = 0, extended = 0;
	edge *ass_eg, *tmp_eg = NULL;
	ext_msg *msg = 0;

	ass_eg = init_ctg(parent, cur_eg, init_query, ori);
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
			query = ext_by_pool(ass_eg, ht, c_pool, query, ori);
			p_query(__func__, query);
			parent = NULL;
		}
		msg = single_ext(ass_eg, c_pool, query, ht, ori); // c_pool freed here.
		bwa_free_read_seq(1, query);
		query = NULL;
		c_pool = NULL; // Content has been free in function single_extend.
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
			used = (bwa_seq_t*) msg->used;
			// If the ori is reverse, just not set the right_ctg value.
			if (used->contig_id >= 0 && used->contig_id != ass_eg->id && !ori) {
				// Set the shifted pos, such that it is able to recover the original transcript
				tmp_eg = (edge*) edgearray_find_id(all_edges, used->contig_id);
				if (!tmp_eg) {
					p_query("ERROR", used);
					err_fatal(
							__func__,
							"There must be something wrong when extending [%d, %d]! Contig %d not found! \n",
							ass_eg->id, ass_eg->len, used->contig_id);
				}
				p_query("USED", used);
				ass_eg->r_shift = used->shift + used->cursor - 1;
				show_debug_msg(__func__,
						"Right connect [%d, %d] to [%d, %d], shift %d \n",
						ass_eg->id, ass_eg->len, tmp_eg->id, tmp_eg->len,
						ass_eg->r_shift);
				ass_eg->right_ctg = tmp_eg;
				g_ptr_array_free(ass_eg->out_egs, TRUE);
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
			c = (int*) msg->counter;
			c_index = 0;
			while (c[c_index] != INVALID_CHAR) {
				sub_query = new_seq(msg->query, opt->ol, msg->query->len
						- opt->ol);
				ext_que(sub_query, c[c_index], ori);
				show_debug_msg(__func__,
						"[%d, %d] Multi-braching %d, level %d\n", ass_eg->id,
						ass_eg->len, c_index, level);
				show_debug_msg(__func__, "Subpath is feasible to go \n");
				tmp_eg = pe_ass_edge(ass_eg, 0, 0, sub_query, ht, level + 1,
						ori);
				p_ctg_seq("BRANCHING", tmp_eg->contig);
				if (post_vld_branch(ass_eg, tmp_eg, ht, ori))
					no_sub_path = 0;
				bwa_free_read_seq(1, sub_query);
				c_index++;
			}
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
	//	FILE *debug = xopen("debug.txt", "a+");

	p_query(__func__, read);
	c_pool = get_init_pool(ht, read, 0);
	// p_pool("Initial Pool: ", c_pool);
	cur_eg = pe_ass_edge(0, 0, c_pool, read, ht, 0, 0);
	reads = cur_eg->reads;
	r_ext_len = cur_eg->len;
	// The c_pool is freed after extension
	c_pool = get_init_pool(ht, read, 1);
	// p_pool("Initial Pool: ", c_pool, NULL);
	show_debug_msg(__func__,
			"Extending to the left ***************************** [%d, %d]\n",
			cur_eg->id, cur_eg->len);
	pe_ass_edge(0, cur_eg, c_pool, read, ht, 0, 1);

	if (cur_eg->len <= opt->rl + 4 && (cur_eg->in_egs->len == 0)
			&& (cur_eg->right_ctg || cur_eg->out_egs->len == 0)) {
		show_msg(__func__, "Not extended. Edge destroyed. \n");
		g_ptr_array_remove_index(all_edges, all_edges->len - 1);
		if (cur_eg->right_ctg) {
			eg_i = cur_eg->right_ctg;
			g_ptr_array_remove_fast(eg_i->in_egs, cur_eg);
		}
		destroy_eg(cur_eg);
		n_reads_consumed++;
		return 0;
	}
	for (i = 0; i < reads->len; i++) {
		r = g_ptr_array_index(reads, i);
		r->shift += (cur_eg->len - r_ext_len);
	}

	for (i = 0; i < all_edges->len; i++) {
		eg_i = g_ptr_array_index(all_edges, i);
		if (eg_i && eg_i->in_egs->len == 0 && eg_i->alive) {
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
	e_index = 10;
	while (fgets(line, 80, solid_reads) != NULL && ht->n_seqs * STOP_THRE
			> n_reads_consumed) {
//		if (counter <= 12000)
			index = atoi(line);
//		else
//			index = (int) (rand_f() * ht->n_seqs);
//		if (counter < s_index) {
//			counter++;
//			continue;
//		}
//		if (counter >= e_index)
//			break;
		t_eclipsed = (float) (clock() - t) / CLOCKS_PER_SEC;
		p = &ht->seqs[index];
		if (p->used || p->contig_id < 0) {
			show_msg(__func__, "Read used: %s\n", p->name);
			continue;
		}
		if (has_rep_pattern(p) || is_repetitive_q(p)) {
			show_msg(__func__, "Read has repeat pattern, skip.\n");
			p_query(__func__, p);
			continue;
		}
		sprintf(msg, "Processing read %d: %s... \n", counter, p->name);
		show_msg(__func__, msg);
		fputs(p->name, start_reads);
		fputs("\n", start_reads);
		eg = pe_ass_ctg(left_rm, p, ht);
		if (eg) {
			readarray_add(eg, p);
			counter++;
			for (i = pre_ctg_id; i < all_edges->len; i++) {
				eg_i = g_ptr_array_index(all_edges, i);
				if (eg_i->alive) {
					n_reads_consumed += eg_i->reads->len;
					show_debug_msg("CONSUMED", "[%d, %d]: %d\n", eg_i->id,
							eg_i->len, eg_i->reads->len);
				}
				//log_edge(eg_i);
			}
			show_msg(__func__, "Reads Consumed: %d; ctg id: %d, %d\n",
					n_reads_consumed, eg->id, eg->len);
		}
		pre_ctg_id = all_edges->len;

		sprintf(msg, "-------------------------------------- %.2f sec \n",
				t_eclipsed);
		show_msg(__func__, msg);
		//break;
	} // All solid reads assembled.

	fprintf(stderr,
			"[pe_ass_core] ------------------------------------------------------ \n");

	// Post process the roadmaps.
	// log_reads(all_edges);
	show_msg(__func__, "Post processing the roadmap... \n");
	graph_by_edges(all_edges, "graph/rm_bf_update.dot");
	save_edges(all_edges, all_contigs, 0, 1, opt->rl * 1.5);
	post_pro(left_rm, all_edges, opt);
	graph_by_edges(all_edges, "graph/rm_after_update.dot");
	save_edges(all_edges, ass_contigs, 0, 0, opt->rl * 1.5);

	free(h);
	free(msg);
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
