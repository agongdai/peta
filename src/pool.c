/*
 * pool.c
 *
 *  Created on: 20-Jun-2011
 *      Author: carl
 */
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>
#include "pool.h"
#include "peseq.h"
#include "edgelist.h"
#include "utils.h"
#include "pelib.h"
#include "pechar.h"
#include "pehash.h"

/**
 * Count the occurrences of next probable chars.
 * Result: set the array 'next' as the counters.
 */
void check_next_char(pool *cur_pool, edge *eg, int *next, const int ori) {
	int i = 0, pre_pos = 0, check_pre = 0;
	bwa_seq_t *s = NULL;
	if (cur_pool->n > 10)
		check_pre = 1;
	for (i = 0; i < cur_pool->n; i++) {
		s = g_ptr_array_index(cur_pool->reads, i);
		pre_pos = ori ? s->cursor + 1 : s->cursor - 1;
		// Only if read's last character is some as the contig's character, count it.
		if (s->rev_com) {
			if (check_pre) {
				if (pre_pos >= 0 && pre_pos < s->len && s->rseq[pre_pos]
						== eg->contig->seq[eg->contig->len - 1])
					check_c(next, s->rseq[s->cursor]);
			} else {
				check_c(next, s->rseq[s->cursor]);
			}
		} else {
			if (check_pre) {
				if (pre_pos >= 0 && pre_pos < s->len && s->seq[pre_pos]
						== eg->contig->seq[eg->contig->len - 1])
					check_c(next, s->seq[s->cursor]);
			} else {
				check_c(next, s->seq[s->cursor]);
			}
		}
	}
}

int should_start(bwa_seq_t *query) {
	int tid = atoi(query->name);
	if (query->status != FRESH)
		return 0;
	// If this read is currently used by another thread
	if (query->is_in_c_pool > 0 && query->is_in_c_pool != tid)
		return 0;
	if (query->is_in_m_pool > 0 && query->is_in_m_pool != tid)
		return 0;
	if (has_n(query, 4) || is_biased_q(query) || has_rep_pattern(query)
			|| is_repetitive_q(query)) {
		query->status = TRIED;
		return 0;
	}
	return 1;
}

void pool_add(pool *p, bwa_seq_t *new_seq, const int tid) {
	if (!new_seq)
		return;
	g_ptr_array_add(p->reads, new_seq);
	new_seq->is_in_c_pool = tid;
	new_seq->tid = tid;
	p->n++;
}

void mate_pool_add(pool *p, bwa_seq_t *new_seq, const int tid) {
	if (!new_seq)
		return;
	g_ptr_array_add(p->reads, new_seq);
	new_seq->is_in_m_pool = tid;
	new_seq->tid = tid;
	p->n = p->reads->len;
}

void pool_uni_add(pool *p, bwa_seq_t *new_seq) {
	if (pool_exists(p, new_seq))
		return;
	pool_add(p, new_seq, 1);
}

void mate_pool_uni_add(pool *p, bwa_seq_t *new_seq) {
	if (pool_exists(p, new_seq))
		return;
	mate_pool_add(p, new_seq, 1);
}

gboolean pool_rm(pool *r_pool, bwa_seq_t *rm_seq) {
	gboolean r;
	r = g_ptr_array_remove(r_pool->reads, rm_seq);
	rm_seq->is_in_c_pool = 0;
	r_pool->n = r_pool->reads->len;
	return r;
}

gboolean pool_rm_fast(pool *p, bwa_seq_t *read) {
	gboolean r;
	r = g_ptr_array_remove(p->reads, read);
	read->is_in_c_pool = 0;
	p->n = p->reads->len;
	return r;
}

gboolean pool_rm_index(pool *p, const int i) {
	gboolean r;
	bwa_seq_t *read = g_ptr_array_index(p->reads, i);
	read->is_in_c_pool = 0;
	read->tid = -1;
	read->pos = -1;
	r = g_ptr_array_remove_index_fast(p->reads, i);
	p->n = p->reads->len;
	return r;
}

gboolean mate_pool_rm(pool *r_pool, bwa_seq_t *rm_seq) {
	gboolean r;
	r = g_ptr_array_remove(r_pool->reads, rm_seq);
	rm_seq->is_in_m_pool = 0;
	rm_seq->tid = -1;
	r_pool->n = r_pool->reads->len;
	return r;
}

gboolean mate_pool_rm_fast(pool *p, bwa_seq_t *read) {
	gboolean r;
	r = g_ptr_array_remove(p->reads, read);
	read->is_in_m_pool = 0;
	read->tid = -1;
	p->n = p->reads->len;
	return r;
}

pool *get_init_mate_pool(pool *cur_pool, bwa_seq_t *seqs, const int ori,
		const int tid) {
	bwa_seq_t *read = NULL, *mate = NULL;
	int i = 0;
	pool *mate_pool = new_pool();
	if (cur_pool) {
		for (i = 0; i < cur_pool->reads->len; i++) {
			read = g_ptr_array_index(cur_pool->reads, i);
			mate = get_mate(read, seqs);
			if (is_paired(read, ori)) {
				if (mate->status != USED && !mate->is_in_m_pool
						&& !mate->is_in_c_pool) {
					mate->rev_com = read->rev_com;
					mate_pool_add(mate_pool, mate, tid);
				}
			}
		}
	}
	return mate_pool;
}

gpointer mate_pool_rm_index(pool *p, const int i) {
	gpointer r;
	bwa_seq_t *read = g_ptr_array_index(p->reads, i);
	read->is_in_m_pool = 0;
	r = g_ptr_array_remove_index_fast(p->reads, i);
	p->n = p->reads->len;
	return r;
}

void pool_get_majority(pool *cur_pool, const char c, edge *ass_eg) {
	bwa_seq_t *s = NULL;
	readarray *reads = cur_pool->reads;
	int i = 0;
	char cursor_char = 0;

	for (i = 0; i < reads->len; i++) {
		s = g_ptr_array_index(reads, i);
		cursor_char = s->rev_com ? s->rseq[s->cursor] : s->seq[s->cursor];
		if (cursor_char != c) {
			pool_rm_index(cur_pool, i);
			s->status = 2; // Will be reused later.
			s->contig_id = ass_eg->id; // Indicates this read will not be used by this edge anymore.
			i--;
		}
	}
}

/**
 * Moving forward by updating the pools
 */
bwa_seq_t *forward(pool *cur_pool, const char c, edge *ass_eg, const int ori) {
	bwa_seq_t *p, *used = 0;
	// Temp reads, holding all valid reads.
	readarray *reads = cur_pool->reads;
	int i;
	for (i = 0; i < cur_pool->n; i++) {
		p = (bwa_seq_t*) g_ptr_array_index(reads, i);
		// Move forward the cursor, update the next character,
		// then assign the read to a temp read.
		// If right mate, the cursor is going backward;
		// otherwise, going forward.
		if (ori)
			p->cursor--;
		else
			p->cursor++;
		// If the cursor char is 'N', just ignore it.
		if (p->seq[p->cursor] == 4) {
			p->status = TRIED;
			p->contig_id = ass_eg->id;
			p->pos = -1;
			pool_rm_index(cur_pool, i);
			// This read is still used by current thread,
			//   only after the pool is freed, the tid is set to be -1
			p->tid = ass_eg->tid;
			i--;
			continue;
		}
		if ((p->cursor < p->len) && p->cursor >= 0) {
			if (p->status == USED) {
				if (!used)
					used = p;
				pool_rm_index(cur_pool, i);
				i--;
			}
		} else { // If the cursor is out of length, mark it as used!
			p->shift = (ass_eg->len - p->cursor + 1);
			pool_rm_index(cur_pool, i);
			readarray_add(ass_eg, p);
			i--;
		}
		p->tid = ass_eg->tid;
	}
	return used;
}

void clean_cur_pool(pool *cur_pool) {
	int i = 0;
	bwa_seq_t *r = NULL;
	readarray *ra = cur_pool->reads;
	for (i = 0; i < ra->len; i++) {
		r = g_ptr_array_index(ra, i);
		if (r->status) {
			if (pool_rm_index(cur_pool, i))
				i--;
		}
	}
}

int get_next_char(pool *cur_pool, const int ori, edge *eg) {
	readarray *reads = cur_pool->reads;
	int *c = (int*) calloc(5, sizeof(int));
	int i = 0, next_char = 0, max_c = 0;
	bwa_seq_t *r;
	c[0] = c[1] = c[2] = c[3] = c[4] = 0;
	for (i = 0; i < reads->len; i++) {
		r = g_ptr_array_index(reads, i);
		if (r->rev_com)
			c[r->rseq[r->cursor]]++;
		else
			c[r->seq[r->cursor]]++;
		if (ori)
			r->cursor--;
		else
			r->cursor++;
		if (r->cursor < 0 || r->cursor >= r->len) {
			pool_rm_fast(cur_pool, r);
			r->shift = (eg->len - r->cursor + 1);
			readarray_add(eg, r);
			i--;
		}
	}
	for (i = 0; i < 5; i++) {
		if (c[i] > max_c) {
			max_c = c[i];
			next_char = i;
		}
	}
	free(c);
	return next_char;
}

int pool_exists(const pool *p, const bwa_seq_t *read) {
	int i = 0;
	readarray *reads = p->reads;
	for (i = 0; i < p->n; i++) {
		if (g_ptr_array_index(reads, i) == read)
			return 1;
	}
	return 0;
}

/**
 * Remove reads from current pool, if the read has two bases inconsistencies at the extension point.
 * The read maybe just share a partial portion
 * contig: aaaaaaaatg
 * read:   aaaaaaaactgggggg
 * 'tg' and not the same as 'ct', remove read from the current pool.
 */

void rm_partial(edge *eg, pool *cur_pool, int ori, bwa_seq_t *seqs,
		bwa_seq_t *query, int nm) {
	int check_c_1 = 0, check_c_2 = 0, confirm_c = 0, confirm_c_2 = 0;
	int removed = 0, is_at_end = 0, i = 0;
	bwa_seq_t *s = NULL;
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
		// If the cursor has reached the end, do not remove it.
		// In case that the read will be removed and not marked as used, which confuses the extending from mates.
		is_at_end = ori ? (s->cursor <= nm) : (s->cursor >= s->len - nm - 1);
		// Remove those reads probably at the splicing junction
		if (!is_at_end) {
			if (s->tid != eg->tid || (check_c_1 != confirm_c && check_c_2
					!= confirm_c_2)) {
				removed = pool_rm_index(cur_pool, i);
				//p_query(__func__, s);
				//p_ctg_seq(__func__, eg->contig);
				s->contig_id = eg->id;
				s->status = FRESH;
				if (removed) {
					s->shift = eg->len - s->cursor + 1;
					i--;
				}
			}
		}
	}
}

void overlap_mate_pool(pool *cur_pool, pool *mate_pool, bwa_seq_t *contig,
		const int ori) {
	int i = 0, overlapped = 0;
	bwa_seq_t *mate = NULL, *tmp = NULL;
	// Add the mate reads which overlap with the tail into the current pool
	for (i = 0; i < mate_pool->n; i++) {
		mate = g_ptr_array_index(mate_pool->reads, i);
		if (mate->is_in_c_pool || (is_right_mate(mate->name) && ori)
				|| (is_left_mate(mate->name) && !ori))
			continue;
		tmp = mate;
		if (mate->rev_com)
			tmp = new_mem_rev_seq(mate, mate->len, 0);
		else
			tmp = new_seq(mate, mate->len, 0);
		if (ori) { // Single extension reverse the template first
			seq_reverse(tmp->len, tmp->seq, 0);
		}
		overlapped = find_ol(contig, tmp, MISMATCHES);
		if (overlapped >= MATE_OVERLAP_THRE) {
			mate->cursor = ori ? (mate->len - overlapped - 1) : overlapped;
			pool_add(cur_pool, mate, 1);
			if (mate_pool_rm_index(mate_pool, i))
				i--;
		}
		bwa_free_read_seq(1, tmp);
	}
}

/**
 * Check whether the character after the current cursor are consistent enough.
 *
 * Two possible branches: c, t
 * aaaaacttaag
 * aaaaatccgtct
 * aaaaatccgtctc
 *
 * The position 6: good count: 2 c's; bad count: 1 't'
 * The position 7: good count: 2 c's; bad count: 1 't'
 * The position 8: good count: 2 g's: bad count: 1 'a'
 * ....
 * Totally: good count: 13; bad: 5. 5/18 > 0.2, return 1, meaning there are actual branches
 */
int bases_sup_branches(pool *cur_pool, const int ori, double threshold) {
	bwa_seq_t *read = NULL;
	int i = 0, index = 0, next_cursor = 0, most_index = 0;
	double good_count = 0, bad_count = 0, more = 1;
	int *sta = (int*) calloc(5, sizeof(int));
	while (more) {
		more = 0;
		index++;
		for (i = 0; i < cur_pool->reads->len; i++) {
			read = g_ptr_array_index(cur_pool->reads, i);
			next_cursor = ori ? read->cursor - index : read->cursor + index;
			if (next_cursor < 0 || next_cursor >= read->len)
				continue;
			more = 1;
			if (read->rev_com) {
				sta[read->rseq[next_cursor]]++;
			} else
				sta[read->seq[next_cursor]]++;
		}
		most_index = get_pure_most(sta);
		for (i = 0; i < 5; i++) {
			if (i == most_index)
				good_count += sta[i];
			else
				bad_count += sta[i];
		}
		reset_c(sta, NULL);
	}
	free(sta);
	//show_debug_msg(__func__, "good_count/bad_count: %.0f/%.0f \n", good_count,
	//		bad_count);
	if (good_count > 0 && bad_count > 0) {
		if (good_count < (good_count + bad_count) * threshold)
			return 1;
	}
	return 0;
}

void clean_mate_pool(reads_ht *rht, pool *mate_pool, edge *eg) {
	bwa_seq_t *mate = NULL;
	int i = 0;
	if (!mate_pool)
		return;
	for (i = 0; i < mate_pool->reads->len; i++) {
		mate = g_ptr_array_index(mate_pool->reads, i);
		if (mate->status == USED || mate->is_in_c_pool == eg->tid
				|| mate->is_in_m_pool != eg->tid || (mate->status == TRIED
				&& mate->contig_id == eg->id)) {
			if (mate_pool_rm_index(mate_pool, i)) {
				i--;
				rm_read_from_ht(rht, mate);
			}
		}
	}
}

pool *new_pool() {
	pool *r_pool = (pool*) malloc(sizeof(pool));
	readarray *reads = g_ptr_array_sized_new(POOLSIZE);
	r_pool->n = 0;
	r_pool->space = POOLSIZE;
	r_pool->reads = reads;
	return r_pool;
}

void clear_pool(pool *r_pool) {
	readarray *reads = r_pool->reads;
	bwa_seq_t *r = 0;
	int i = 0;
	for (i = 0; i < r_pool->n; i++) {
		r = g_ptr_array_index(reads, i);
		r->is_in_c_pool = 0;
		r->is_in_m_pool = 0;
		r->tid = -1;
	}
	while (reads->len > 0)
		g_ptr_array_remove_index(r_pool->reads, 0);
	r_pool->reads->len = 0;
	r_pool->n = 0;
}

void free_pool(pool *r_pool) {
	readarray *reads = NULL;
	bwa_seq_t *r = NULL;
	// p_pool("Pool to free: ", r_pool);
	int i = 0;
	// show_debug_msg(__func__, "Pool: %p \n", r_pool);
	if (r_pool) {
		reads = r_pool->reads;
		if (r_pool->n > 0) {
			for (i = 0; i < r_pool->n; i++) {
				r = g_ptr_array_index(reads, i);
				r->is_in_c_pool = 0;
				r->tid = -1;
			}
			r_pool->n = 0;
			if (reads) {
				g_ptr_array_free(reads, TRUE);
			}
		}
		free(r_pool);
		r_pool = 0;
	}
}

void free_mate_pool(pool *mate_pool) {
	readarray *reads = mate_pool->reads;
	bwa_seq_t *r = 0;
	// p_pool("Pool to free: ", mate_pool);
	int i = 0;
	// show_debug_msg(__func__, "Pool: %p \n", mate_pool);
	if (mate_pool) {
		if (mate_pool->n > 0) {
			for (i = 0; i < mate_pool->n; i++) {
				r = g_ptr_array_index(reads, i);
				r->is_in_m_pool = 0;
			}
			mate_pool->n = 0;
			if (reads) {
				g_ptr_array_free(reads, TRUE);
			}
		}
		free(mate_pool);
		mate_pool = 0;
	}
}

void p_pool_read(gpointer *data, gpointer *user_data) {
	bwa_seq_t *p = (bwa_seq_t*) data, *p2;
	char c;
	int i = 0;
	if ((unsigned int) p->seq[0] <= 4) {
		p2 = new_seq(p, p->len, 0);
		map(p2);
		c = p->rev_com ? p2->rseq[p2->cursor] : p2->seq[p2->cursor];
		if (p2->cursor < 0 || p2->cursor >= p2->len)
			c = 'n';
		for (i = 0; i < p2->len - p2->cursor; i++)
			printf(" ");

		if (p->rev_com) {
			printf("%s", p2->rseq);
			for (i = 0; i < p2->cursor + 2; i++)
				printf(" ");
			printf("%d->%c@%d\t%s\t[pool: %d]\t[tid: %d]\t[rev_com]",
					p2->cursor, c, p->status, p2->name, p->is_in_c_pool, p->tid);
		} else {
			printf("%s", p2->seq);
			for (i = 0; i < p2->cursor + 2; i++)
				printf(" ");
			printf("%d->%c@%d\t%s\t[pool: %d]\t[tid: %d]",
					p2->cursor, c, p->status, p2->name, p->is_in_c_pool, p->tid);
		}
		printf("\n");
		bwa_free_read_seq(1, p2);
	} else {
		c = p->rev_com ? p->rseq[p->cursor] : p->seq[p->cursor];
		if (p->cursor < 0 || p->cursor >= p->len)
			c = 'n';
		printf("[p_pool] %d_%d: %s %s %d->%c\n", p->contig_id, p->shift,
				p->name, p->seq, p->cursor, c);
		//printf("[p_pool] %d_%d: %s %s %d->%c %d\n", p->contig_id, p->shift,
		//		p->name, p->rseq, p->cursor, c, p->rev_com);
	}
}

void p_edge_reads(const edge *eg) {
	int i = 0;
	bwa_seq_t *p;
	show_debug_msg(__func__, "--------------------------------- \n");
	show_debug_msg(__func__, "Reads assemblying edge [%d, %d] \n", eg->id,
			eg->len);
	for (i = 0; i < eg->reads->len; i++) {
		p = g_ptr_array_index(eg->reads, i);
		show_debug_msg(__func__, "%d: %d_%d Read %s: %d \n", i, p->contig_id,
				p->shift, p->name, p->shift);
	}
	show_debug_msg(__func__, "--------------------------------- \n");
}

void p_readarray(const readarray *ra, const int all) {
	int i = 0;
	bwa_seq_t *p;
	show_debug_msg(__func__, "--------------------------------- \n");
	show_debug_msg(__func__, " # of Reads: %d \n", ra->len);
	for (i = 0; i < ra->len; i++) {
		p = g_ptr_array_index(ra, i);
		if ((!all && i % 200 == 0) || all)
			p_query(__func__, p);
		//			show_debug_msg(__func__, "%d: %d_%d Read %s: %d \n", i,
		//					p->contig_id, p->shift, p->name, p->shift);
	}
	show_debug_msg(__func__, "--------------------------------- \n");
}

void p_pool(const char *header, const pool *r_pool, const int *next) {
	readarray *reads = r_pool->reads;
	printf("[p_pool]****************************** \n");
	if (next)
		printf("[p_pool] %s %zd: %d:%d:%d:%d\n", header, r_pool->n, next[0],
				next[1], next[2], next[3]);
	else
		printf("[p_pool] %s %zd:\n", header, r_pool->n);
	g_ptr_array_foreach(reads, (GFunc) p_pool_read, NULL);
	printf("[p_pool]****************************** \n");
}
