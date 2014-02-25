/*
 * pool.c
 *
 *  Created on: Jul 23, 2013
 *      Author: carl
 */
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "pool.hpp"
#include "peseq.h"
#include "utils.h"
#include "pechar.h"
#include "k_hash.h"
#include "tpl.hpp"
#include "junction.hpp"

using namespace std;

pool *new_pool() {
	pool *p = (pool*) malloc(sizeof(pool));
	p->reads = g_ptr_array_sized_new(0);
	return p;
}

/**
 * Destroy a pool
 */
void destroy_pool(pool *p) {
	int i = 0;
	bwa_seq_t *r = NULL;
	if (p) {
		if (p->reads) {
			//show_debug_msg(__func__, "Reads: %d \n", p->reads->len);
			for (i = 0; i < p->reads->len; i++) {
				r = (bwa_seq_t*) g_ptr_array_index(p->reads, i);
				//p_query(__func__, r);
				reset_to_fresh(r);
			}
			g_ptr_array_free(p->reads, TRUE);
		}
		free(p);
	}
}

void p_pool_read(gpointer *data, gpointer *user_data) {
	bwa_seq_t *p = (bwa_seq_t*) data, *p2 = NULL;
	char c;
	int i = 0;
	if ((unsigned int) p->seq[0] <= 4) {
		p2 = new_seq(p, p->len, 0);
		map(p2);
		c = p->rev_com ? p2->rseq[p2->cursor] : p2->seq[p2->cursor];
		if (p2->cursor < 0 || p2->cursor >= p2->len)
			c = 'N';
		for (i = 0; i < p2->len - p2->cursor; i++)
			printf(" ");

		if (p->rev_com) {
			for (i = 0; i < p2->cursor; i++) {
				printf("%c", p2->rseq[i]);
			}
			printf(" %c ", p2->rseq[p2->cursor]);
			for (i = p2->cursor + 1; i < p2->len; i++) {
				printf("%c", p2->rseq[i]);
			}

			for (i = 0; i < p2->cursor + 2; i++)
				printf(" ");
			printf("%d->%c@%s\t[status: %d]\t[    <<<< %d]", p2->cursor, c,
					p2->name, p->status, p->pos);
		} else {
			for (i = 0; i < p2->cursor; i++) {
				printf("%c", p2->seq[i]);
			}
			printf(" %c ", p2->seq[p2->cursor]);
			for (i = p2->cursor + 1; i < p2->len; i++) {
				printf("%c", p2->seq[i]);
			}

			for (i = 0; i < p2->cursor + 2; i++)
				printf(" ");
			printf("%d->%c@%s\t[status: %d]\t[>>>>     %d]", p2->cursor, c,
					p2->name, p->status, p->pos);
		}
		printf("\t\n");
		bwa_free_read_seq(1, p2);
	}
}

void p_pool(const char *header, const pool *p, const int *next) {
	readarray *reads = p->reads;
	g_ptr_array_sort(reads, (GCompareFunc) cmp_reads_by_cursor);
	printf("[p_pool]****************************** \n");
	if (next)
		printf("[p_pool] %s %zd: %d:%d:%d:%d\n", header, reads->len, next[0],
				next[1], next[2], next[3]);
	else
		printf("[p_pool] %s %zd:\n", header, reads->len);
	g_ptr_array_foreach(reads, (GFunc) p_pool_read, NULL);
	printf("[p_pool]****************************** \n");
}

void p_pool_read_mates(char *header, bwa_seq_t *seqs, pool *p) {
	printf("[%s]*********************************\n", header);
	int i = 0;
	bwa_seq_t *r = NULL, *m = NULL;
	for (i = 0; i < p->reads->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(p->reads, i);
		m = get_mate(r, seqs);
		p_query(header, r);
		p_query(header, m);
		printf("[%s]\n", header);
	}
	printf("[%s]*********************************\n", header);
}

/**
 * Add a read to the pool
 */
void add2pool(pool *p, bwa_seq_t *r) {
	g_ptr_array_add(p->reads, r);
	r->status = IN_POOL;
	// Indicates how many mismatches between the read and template
	if (r->pos == IMPOSSIBLE_NEGATIVE)
		r->pos = 0;
}

/**
 * When the extension is terminated in the middle, mark the reads in pool as TRIED.
 */
void mark_pool_reads_tried(pool *p, tpl *t) {
	int i = 0;
	bwa_seq_t *r = NULL;
	for (i = 0; i < p->reads->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(p->reads, i);
		add2tried(t, r);
	}
}

void empty_pool(pool *p) {
	if (!p || !p->reads || p->reads->len <= 0)
		return;
	bwa_seq_t *r = NULL;
	while (p->reads->len > 0) {
		r = (bwa_seq_t*) g_ptr_array_index(p->reads, 0);
		reset_to_fresh(r);
		g_ptr_array_remove_index_fast(p->reads, 0);
	}
}

/**
 * Remove a read from the pool and reset the the read status
 */
void rm_from_pool(pool *p, int index) {
	if (index < 0 || !p || !p->reads || index >= p->reads->len)
		return;
	bwa_seq_t * r = (bwa_seq_t*) g_ptr_array_index(p->reads, index);
	reset_to_fresh(r);
	g_ptr_array_remove_index_fast(p->reads, index);
}

/**
 * Count frequencies of next characters, get the most frequent one.
 * If -1: all of them are 0
 */
int get_next_char(hash_table *ht, pool *p, tpl *t,
		const int ori) {
	readarray *reads = p->reads;
	junction *jun = NULL;
	float *c = NULL, weight = 0.0, max_c = 0.0, multi = INS_SIZE;
	int i = 0, j = 0, next_char = -1;
	bwa_seq_t *r = NULL, *m = NULL;
	int pre_cursor = 0, this_c = 0, pre_c = 0, counted = 0, n_pairs = 0;
	tpl *main_tpl = NULL;
	int pre_t_c = ori ? t->ctg->seq[0] : t->ctg->seq[t->len - 1];

	if (!p->reads || p->reads->len == 0)
		return -1;

	c = (float*) calloc(5, sizeof(float));
	c[0] = c[1] = c[2] = c[3] = c[4] = 0;
	for (i = 0; i < reads->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(reads, i);
		pre_cursor = ori ? (r->cursor + 1) : (r->cursor - 1);
		// Check only if there are more than 2 reads
		// If in the previous round, this read has a mismatch, then does not count it this time
		if (pre_t_c > -1 && pre_cursor >= 0 && pre_cursor < r->len) {
			pre_c = r->rev_com ? r->rseq[pre_cursor] : r->seq[pre_cursor];
			if (pre_c != 4 && pre_c != pre_t_c)
				continue;
		}
		this_c = r->rev_com ? r->rseq[r->cursor] : r->seq[r->cursor];
		// Simply ignore 'N's
		if (this_c == 4)
			continue;
		m = get_mate(r, ht->seqs);
		// The overlap length with template
		weight = ori ? (r->len - r->cursor - 1) : r->cursor;
		// Minus the mismatches with the template
		weight -= r->pos * MISMATCH_WEIGHT;
		// If its mate is USED or IN_POOL, it is guaranteed to be nearby by next_pool().
		if (m->status == IN_POOL || m->status == USED) {
			weight += multi;
		}
		c[this_c] += weight;
		counted++;
		// At most count MAX_POOL_N_READS reads, in case the pool is large
		if (counted >= MAX_POOL_N_READS)
			break;
	}
	// In case only few reads in pool, and all of them not support the previous base
	if (counted == 0) {
		for (i = 0; i < reads->len; i++) {
			r = (bwa_seq_t*) g_ptr_array_index(reads, i);
			this_c = r->rev_com ? r->rseq[r->cursor] : r->seq[r->cursor];
			weight = ori ? (r->len - r->cursor - 1) : r->cursor;
			weight -= r->pos * MISMATCH_WEIGHT;
			c[this_c] += weight;
		}
	}
	// Do not count the 'N's
	for (i = 0; i < 4; i++) {
		if (c[i] > max_c) {
			max_c = c[i];
			next_char = i;
		}
	}
	if (max_c == 0.0)
		next_char = -1;
	free(c);
	return next_char;
}

/**
 * Moving forward by updating the pools
 * If some read is marked as USED, return 1
 */
int forward(pool *p, tpl *t, const int ori) {
	bwa_seq_t *r = NULL;
	readarray *reads = p->reads;
	int i = 0, locus = 0;
	int some_read_used = 0;
	for (i = 0; i < reads->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(reads, i);
		// Move forward the cursor, update the next character,
		// If right mate, the cursor is going backward;
		// otherwise, going forward.
		if (ori)
			r->cursor--;
		else
			r->cursor++;
		// If the cursor is out of length, mark it as used, remove from pool
		if ((r->cursor >= r->len) || r->cursor < 0) {
			// If extending to the left, the locus value is temp, updated by upd_locus_on_tpl
			locus = ori ? t->len : (t->len - r->len);
			add2tpl(t, r, locus);
			g_ptr_array_remove_index_fast(reads, i--);
			some_read_used = 1;
		}
	}
	return some_read_used;
}

int has_n_in_pool(pool *p) {
	int i = 0;
	bwa_seq_t *r = NULL;
	for (i = 0; i < p->reads->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(p->reads, i);
		if (has_n(r, 1))
			return 1;
	}
	return 0;
}

/**
 * If a read is at the junction, maybe partially used and later removed
 * The read attribute 'pos' indicates how many mismatches between read and the template
 */
void rm_half_clip_reads(pool *p, tpl *t, int tpl_c, int mismatches, int ori) {
	index64 i = 0;
	bwa_seq_t *r = NULL;
	int read_c = 0;

	for (i = 0; i < p->reads->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(p->reads, i);
		read_c = r->rev_com ? r->rseq[r->cursor] : r->seq[r->cursor];
		// If the character in read is 'N', do not count it as mismatch
		if (read_c == 4 || read_c != tpl_c)
			r->pos++;
		if (r->pos > mismatches) {
			//p_query("REMOVED", r);
			//p_pool("BEFORE", p, NULL);
			rm_from_pool(p, i--);
			// Set status to TRIED first. Will be reset to FRESH after this template is done.
			add2tried(t, r);
		}
	}
}

/**
 * If 2 bases before current cursor is not the same as template, remove from pool
 */
void rm_bad_ol_reads(pool *p, tpl *t, const int ori) {
	int i = 0, pre_cursor = 0, pre_pre_cursor = 0;
	bwa_seq_t *r = NULL;
	ubyte_t pre_tpl_c = 0, pre_pre_tpl_c = 0, pre_read_c = 0, pre_pre_read_c =
			0;
	if (t->len < 2)
		return;
	pre_tpl_c = ori ? t->ctg->seq[0] : t->ctg->seq[t->len - 1];
	pre_pre_tpl_c = ori ? t->ctg->seq[1] : t->ctg->seq[t->len - 2];
	for (i = 0; i < p->reads->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(p->reads, i);
		pre_cursor = ori ? r->cursor + 1 : r->cursor - 1;
		pre_pre_cursor = ori ? r->cursor + 2 : r->cursor - 2;
		if (pre_cursor < 0 || pre_cursor >= r->len || pre_pre_cursor < 0
				|| pre_pre_cursor >= r->len)
			continue;
		pre_read_c = r->rev_com ? r->rseq[pre_cursor] : r->seq[pre_cursor];
		pre_pre_read_c = r->rev_com ? r->rseq[pre_pre_cursor]
				: r->seq[pre_pre_cursor];
		if (pre_read_c != pre_tpl_c && pre_pre_read_c != pre_pre_tpl_c) {
			rm_from_pool(p, i--);
			add2tried(t, r);
		}
	}
}

/**
 * Align read length template tails to get initial pool
 */
void init_pool(hash_table *ht, pool *p, tpl *t, int tail_len, int mismatches,
		const int ori) {
	bwa_seq_t *r = NULL, *tail = NULL, *read = NULL;
	index64 i = 0, j = 0, start = 0, end = ht->o->read_len;
	GPtrArray *hits = NULL;
	int cursor = 0;

	// Must be some error happens!
	if (!t->start_read || t->start_read->len <= 0)
		return;
	read = get_tail(t, ht->o->read_len, ori);
	//p_tpl(t);
	//p_query("Init read", read);
	start = ori ? 0 : 1;
	end = ori ? read->len - tail_len - 1 : read->len - tail_len;
	for (i = start; i <= end; i++) {
		tail = new_seq(read, tail_len, i);
		//p_query("TAIL", tail);
		hits = align_tpl_tail(ht, t, tail, 0, i, mismatches, FRESH, ori);
		for (j = 0; j < hits->len; j++) {
			r = (bwa_seq_t*) g_ptr_array_index(hits, j);
			//p_query("CANDIDATE", r);
			if (r->cursor >= 0 && r->cursor < read->len) {
				add2pool(p, r);
			} else {
				reset_to_fresh(r);
			}
		}
		bwa_free_read_seq(1, tail);
		g_ptr_array_free(hits, TRUE);
		//show_debug_msg(__func__, "-----------------------------------\n\n");
	}
	bwa_free_read_seq(1, read);
}

/**
 * To check whether a mate read is in range
 */
int mate_is_in_range(int used_locus, int tpl_len, int read_len, int ori) {
	if (ori) {
		return in_range(used_locus - read_len, tpl_len);
	} else {
		return in_range(used_locus, tpl_len);
	}
}

/**
 * Traverse the splicing graph, check whether the distance of a pair is within range.
 */
int is_good_pair(GPtrArray *near_tpls, tpl *t, bwa_seq_t *r, bwa_seq_t *m, int ori) {
	// @TODO: if the mates overlap, treat as good now
	if (m->status == IN_POOL)
		return 1;
	// If they are on the same template, check distance only
	if (m->status == USED && mate_is_in_range(m->contig_locus, t->len, m->len, ori))
		return 1;
	int i = 0, is_near = 0;
	tpl *near = NULL;
	if (m->status == USED) {
		for (i = 0; i < near_tpls->len; i++) {
			near = (tpl*) g_ptr_array_index(near_tpls, i);
			if (near->id == m->contig_id) {
				is_near = 1;
				break;
			}
		}
	}
	// If the mate is not near, not good pair
	if (!is_near) return 0;
	// @TODO: Check their distance on all possible paths
	return 1;
}

/**
 * Align the tail to the RNA-seq reads, add fresh reads to current pool
 */
void next_pool(hash_table *ht, pool *p, GPtrArray *near_tpls, tpl *t, bwa_seq_t *tail,
		int mismatches, const int ori) {
	GPtrArray *fresh_reads = NULL;
	index64 i = 0, shift = 0;
	bwa_seq_t *r = NULL, *m = NULL;
	int limit = 0;

	// Find all fresh reads and add to the pool
	if (t->len == -1 && t->id == -1)
		p_query(__func__, tail);
	shift = ori ? 0 : ht->o->read_len - tail->len;
	fresh_reads = align_tpl_tail(ht, t, tail, limit, shift, mismatches, FRESH,
			ori);
	for (i = 0; i < fresh_reads->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(fresh_reads, i);
		m = get_mate(r, ht->seqs);

		if (t->id == 108220) {
			p_query(__func__, r);
			p_query(__func__, m);
		}
		//add2pool(p, r);
		// Add single read
		if (m->status == FRESH)
			add2pool(p, r);
		else {
			// Add paired-end read
			if (is_good_pair(near_tpls, t, r, m, ori)) {
				add2pool(p, r);
			}
		}
		if (r->status != IN_POOL)
			reset_to_fresh(r);
		if (t->id == 108220) {
			p_query(__func__, r);
			p_query(__func__, m);
		}
	}
	g_ptr_array_free(fresh_reads, TRUE);
}

/**
 * Correct bases on the template to the consensus base
 */
void correct_init_tpl_base(pool *p, tpl *t, int ori) {
	int i = 0, j = 0, pos = 0;
	ubyte_t c = 0, max_c = 0, rev_c = 0;
	bwa_seq_t *r = NULL;
	int counter[5], max = 0;
	int n_counted = 0;
	if (!p || !p->reads || p->reads->len < 3)
		return;
	//p_ctg_seq("ORIGINAL ", t->ctg);
	for (i = 1; i < t->len; i++) {
		max = 0;
		max_c = 0;
		n_counted = 0;
		for (j = 0; j < 5; j++) {
			counter[j] = 0;
		}
		for (j = 0; j < p->reads->len; j++) {
			r = (bwa_seq_t*) g_ptr_array_index(p->reads, j);
			pos = ori ? r->cursor + 1 + i : r->cursor - r->len + i;
			if (pos >= 0 && pos < r->len) {
				c = r->rev_com ? r->rseq[pos] : r->seq[pos];
				// If 'N', ignore
				if (c == 4)
					continue;
				counter[c]++;
				n_counted++;
				if (n_counted >= MAX_POOL_N_READS)
					break;
			}
		}
		for (j = 0; j < 4; j++) {
			if (counter[j] > max) {
				max = counter[j];
				max_c = j;
			}
		}
		if (max <= 0)
			continue;
		t->ctg->seq[i] = max_c;
		rev_c = 3 - max_c;
	} // End of template base correction
	set_rev_com(t->ctg);
	//p_ctg_seq("CORRECTED ", t->ctg);
}

/**
 * Find the fresh mates which overlap with the template tail, at least 11bp with 0 mismatches
 */
void find_match_mates(hash_table *ht, pool *p, GPtrArray *near_tpls, tpl *t,
		bwa_seq_t *to_match, int mismatches, int ori) {
	bwa_seq_t *m = NULL, *r = NULL, *part = NULL, *tail = NULL;
	int i = 0, j = 0;
	int ol = 0, rev_com = 0, n_mis = 0, cursor = 0;
	junction *jun = NULL;
	tpl *near = NULL;
	GPtrArray *existing_reads = NULL;

	if (!to_match)
		tail = get_tail(t, ht->o->read_len, ori);
	else
		tail = new_seq(to_match, to_match->len, 0);

	//p_query(__func__, tail);
	// In case the tail is an biased seq like: TTTTCTTTTTT
	if (!tail || is_biased_q(tail) || has_n(tail, 1)) {
		bwa_free_read_seq(1, tail);
		return;
	}

	// Find all FRESH mate reads
	existing_reads = g_ptr_array_sized_new(t->reads->len);
	for (i = 0; i < near_tpls->len; i++) {
		near = (tpl*) g_ptr_array_index(near_tpls, i);
		for (j = 0; j < near->reads->len; j++) {
			r = (bwa_seq_t*) g_ptr_array_index(near->reads, j);
			m = get_mate(r, ht->seqs);
			if (m->status == FRESH && !is_bad_query(m)) {
				g_ptr_array_add(existing_reads, r);
			}
		}
	}

	//p_tpl(t);
	//show_debug_msg(__func__, "ORI: %d \n", ori);
	//p_query(__func__, tail);

	for (i = 0; i < existing_reads->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(existing_reads, i);
		m = get_mate(r, ht->seqs);
		if (t->id == 108220) {
			p_query(__func__, r);
			p_query(__func__, m);
			p_query("TAIL", tail);
		}
		// If the used read is on this template but not in range, not consider
		if (r->contig_id == t->id && !mate_is_in_range(r->contig_locus, t->len, r->len, ori)) {
			continue;
		}
		// like 'ATATATATAT' or 'AAAAAAAA'
		if (is_bad_query(m))
			continue;
		// Find the overlapping between mate and tail
		ol = find_fr_ol_within_k(m, tail, mismatches, ht->o->k - 1,
				ht->o->read_len, ori, &rev_com, &n_mis);
		if (t->id == 108220)
		show_debug_msg(__func__, "OVERLAP: %d \n", ol);
		if (r->rev_com == rev_com && ol >= ht->o->k - 1 && ol >= n_mis
				* ht->o->k) {
			part = ori ? new_seq(tail, ol, 0) : new_seq(tail, ol,
					ht->o->read_len - ol);
			cursor = ori ? (m->len - ol - 1) : ol;
			// Make sure the overlapped region has no 'N's
			if (has_n(part, 1) || cursor < 0 || cursor >= m->len) {
				bwa_free_read_seq(1, part);
				continue;
			}
			bwa_free_read_seq(1, part);

			m->rev_com = rev_com;
			m->cursor = cursor;
			m->pos = n_mis;
			add2pool(p, m);
		}
	}
	bwa_free_read_seq(1, tail);
	g_ptr_array_free(existing_reads, TRUE);
}

/**
 * If there are reads with zero mismatches, remove others
 * Otherwise, keep those with one mismatches
 * Otherwise, keep all
 */
void keep_fewer_mis_reads(pool *p) {
	int i = 0;
	bwa_seq_t *r = NULL;
	int has_zero = 0, has_one = 0;
	for (i = 0; i < p->reads->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(p->reads, i);
		if (r->pos == 0)
			has_zero = 1;
		else if (r->pos == 1)
			has_one = 1;
		if (has_zero)
			break;
	}
	if (has_zero) {
		for (i = 0; i < p->reads->len; i++) {
			r = (bwa_seq_t*) g_ptr_array_index(p->reads, i);
			if (r->pos != 0) {
				rm_from_pool(p, i--);
			}
		}
	} else if (has_one) {
		for (i = 0; i < p->reads->len; i++) {
			r = (bwa_seq_t*) g_ptr_array_index(p->reads, i);
			if (r->pos != 1) {
				rm_from_pool(p, i--);
			}
		}
	}
}

/**
 * Find those reads with at least 11 * 2 overlap with the template,
 * 	by searching the hash table with LESS_MISMATCH
 */
void find_ol_mates(hash_table *ht, pool *p, GPtrArray *near_tpls, tpl *t,
		bwa_seq_t *tail, int mismatches, int ori) {
	// Base by base overlapping.
	find_match_mates(ht, p, near_tpls, t, NULL, mismatches, ori);
	rm_bad_ol_reads(p, t, ori);
	keep_fewer_mis_reads(p);
}
