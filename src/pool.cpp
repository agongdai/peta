/*
 * pool.c
 *
 *  Created on: Jul 23, 2013
 *      Author: carl
 */
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "pool.hpp"
#include "peseq.h"
#include "utils.h"
#include "pechar.h"
#include "k_hash.h"
#include "tpl.hpp"

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
	uint32_t i = 0;
	bwa_seq_t *r = NULL;
	if (p) {
		if (p->reads) {
			for (i = 0; i < p->reads->len; i++) {
				r = (bwa_seq_t*) g_ptr_array_index(p->reads, i);
				r->pos = IMPOSSIBLE_NEGATIVE;
				r->cursor = -1;
				if (r->status == IN_POOL)
					r->status = FRESH;
			}
		}
		g_ptr_array_free(p->reads, TRUE);
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
		printf("\n");
		bwa_free_read_seq(1, p2);
	}
}

void p_pool(const char *header, const pool *p, const int *next) {
	readarray *reads = p->reads;
	printf("[p_pool]****************************** \n");
	if (next)
		printf("[p_pool] %s %zd: %d:%d:%d:%d\n", header, reads->len, next[0],
				next[1], next[2], next[3]);
	else
		printf("[p_pool] %s %zd:\n", header, reads->len);
	g_ptr_array_foreach(reads, (GFunc) p_pool_read, NULL);
	printf("[p_pool]****************************** \n");
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
int get_next_char(hash_table *ht, pool *p, tpl *t, const int ori) {
	readarray *reads = p->reads;
	float *c = NULL, weight = 0.0, max_c = 0.0;
	int i = 0, next_char = -1;
	bwa_seq_t *r = NULL, *m = NULL;
	int pre_cursor = 0, this_c = 0, pre_c = 0, counted = 0;
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
		// If its mate is on the same template, triple the weight
		if (m->status == USED && m->contig_id == t->id)
			weight *= MATE_MULTI;
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
 */
void forward(pool *p, tpl *t, const int ori) {
	bwa_seq_t *r = NULL;
	readarray *reads = p->reads;
	int i = 0, locus = 0;
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
		}
	}
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
		if (read_c != 4 && read_c != tpl_c)
			r->pos++;
	}

	for (i = 0; i < p->reads->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(p->reads, i);
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
	start = ori ? 0 : 1;
	end = ori ? read->len - tail_len - 1 : read->len - tail_len;
	for (i = start; i <= end; i++) {
		tail = new_seq(read, tail_len, i);
		//p_query("TAIL", tail);
		hits = align_tpl_tail(ht, t, tail, i, mismatches, FRESH, ori);
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
 * Align the tail to the RNA-seq reads, add fresh reads to current pool
 */
void next_pool(hash_table *ht, pool *p, tpl *t, bwa_seq_t *tail,
		int mismatches, const int ori) {
	GPtrArray *fresh_reads = NULL;
	index64 i = 0, shift = 0;
	bwa_seq_t *r = NULL;

	// Find all fresh reads and add to the pool
	shift = ori ? 0 : ht->o->read_len - tail->len;
	fresh_reads = align_tpl_tail(ht, t, tail, shift, mismatches, FRESH, ori);
	for (i = 0; i < fresh_reads->len; i++) {
		//p_query(__func__, r);
		r = (bwa_seq_t*) g_ptr_array_index(fresh_reads, i);
		//show_debug_msg(__func__, "Read %s cursor: %d\n", r->name, r->cursor);
		add2pool(p, r);
	}
	g_ptr_array_free(fresh_reads, TRUE);
}

/**
 * Correct bases on the template to the concensus base
 */
void correct_init_tpl_base(pool *p, tpl *t, int t_len) {
	int i = 0, j = 0, pos = 0;
	ubyte_t c = 0, max_c = 0, rev_c = 0;
	bwa_seq_t *r = NULL;
	int counter[5], max = 0;
	int n_counted = 0;
	if (!p || t_len <= 0 || t->len < t_len || !p->reads || p->reads->len < 3)
		return;
	p_ctg_seq("ORIGINAL ", t->ctg);
	for (i = 1; i < t->len; i++) {
		max = 0;
		max_c = 0;
		n_counted = 0;
		for (j = 0; j < 5; j++) {
			counter[j] = 0;
		}
		for (j = 0; j < p->reads->len; j++) {
			r = (bwa_seq_t*) g_ptr_array_index(p->reads, j);
			pos = r->cursor - t->len + i;
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
		//show_debug_msg(__func__, "BASES at %d \n", i);
		for (j = 0; j < 4; j++) {
			//show_debug_msg(__func__, "\t%d BASE %c: %d\n", i, "ACGTN"[j],
			//		counter[j]);
			if (counter[j] > max) {
				max = counter[j];
				max_c = j;
			}
		}
		//show_debug_msg(__func__, "Max at %d: [%c, %d] \n", i, "ACGTN"[max_c],
		//		max);
		if (max <= 0)
			continue;
		t->ctg->seq[i] = max_c;
		rev_c = 3 - max_c;
		// If more than HIGH_N_READS (50) reads in pool, correct the reads as well
		// This could be slow if too many reads in pool
		if (n_counted >= HIGH_N_READS && n_counted < MAX_POOL_N_READS) {
			for (j = 0; j < p->reads->len; j++) {
				r = (bwa_seq_t*) g_ptr_array_index(p->reads, j);
				// Number of mismatches against the template is 0 now.
				r->pos = 0;
				pos = r->cursor - t->len + i;
				if (pos >= 0 && pos < r->len) {
					c = r->rev_com ? r->rseq[pos] : r->seq[pos];
					if (c != max_c) {
						//show_debug_msg(__func__,
						//		"Read %s at pos %d: %c => %c \n", r->name,
						//		pos, "ACGTN"[c], "ACGTN"[max_c]);
						if (r->rev_com) {
							r->rseq[pos] = max_c;
							r->seq[r->len - pos - 1] = rev_c;
						} else {
							r->seq[pos] = max_c;
							r->rseq[r->len - pos - 1] = rev_c;
						}
					}
				} // End of correction of a read
			} // End of reads loop
		} // End of correcting bases on reads
	} // End of template base correction
	set_rev_com(t->ctg);
	p_ctg_seq("CORRECTED ", t->ctg);
}

/**
 * Find the fresh mates which overlap with the template tail, at least 11bp with 0 mismatches
 */
void find_match_mates(hash_table *ht, pool *p, tpl *t, int tail_len,
		int mismatches, int ori) {
	bwa_seq_t *m = NULL, *r = NULL, *tail = NULL, *part = NULL;
	index64 i = 0;
	int ol = 0, rev_com = 0, n_mis = 0;
	tail = get_tail(t, tail_len, ori);

	// In case the tail is an biased seq like: TTTTCTTTTTT
	if (is_biased_q(tail) || has_n(tail, 1)) {
		bwa_free_read_seq(1, tail);
		return;
	}

	//p_tpl(t);
	//show_debug_msg(__func__, "ORI: %d \n", ori);
	//p_query(__func__, tail);

	for (i = 0; i < t->reads->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(t->reads, i);
		m = get_mate(r, ht->seqs);

		// If the mate is used already
		// If the orientation is not correct
		if (m->status != FRESH || is_paired(r, ori) || is_biased_q(m)) {
			continue;
		}
		// Find the overlapping between mate and tail
		ol = find_fr_ol_within_k(m, tail, mismatches, ht->o->k - 1, tail_len - 1,
				ori, &rev_com, &n_mis);
		//p_query("USED ", r);
		//p_query("FRESH", m);
		//show_debug_msg(__func__, "OVERLAP: %d\n", ol);

		if (ol >= ht->o->k - 1 && ol >= n_mis * ht->o->k) {
			part = ori ? new_seq(tail, ol, 0) : new_seq(tail, ol, tail->len - ol);
			//p_query(__func__, part);
			if (is_biased_q(part) || has_n(part, 1)) {
				bwa_free_read_seq(1, part);
				continue;
			}
			bwa_free_read_seq(1, part);

			m->rev_com = rev_com;
			m->cursor = ori ? (m->len - ol - 1) : ol;
			m->pos = n_mis;
			add2pool(p, m);
		}
	}
	bwa_free_read_seq(1, tail);
}

/**
 * Find those reads with at least 11 * 2 overlap with the template,
 * 	by searching the hash table with LESS_MISMATCH
 */
void find_hashed_mates(hash_table *ht, pool *p, tpl *t, int full_tail_len,
		int mismatches, int ori) {
	int tail_len = ht->o->k * 2;
	int i = 0, ol = 0, rev_com = 0, n_mis = 0, shift = 0;
	int added = 0;
	bwa_seq_t *seqs = ht->seqs;
	// Query tail and overlap tail are used in different length
	bwa_seq_t *ol_tail = 0, *q_tail = NULL, *r = NULL, *m = NULL;
	GPtrArray *mates = NULL;

	if (t->len < tail_len || t->len < 0)
		return;
	//show_debug_msg(__func__, "Looking for shorter; ori %d...\n", ori);
	// Query tail: shorter than a normal tail, just 22bp, query 2 kmers.
	q_tail = get_tail(t, tail_len, ori);
	// In case the tail is an biased seq like: TTTTCTTTTTT
	if (is_biased_q(q_tail) || has_n(q_tail, 1)) {
		bwa_free_read_seq(1, q_tail);
		return;
	}
	// Overlap tail. Used to verify the overlapping between mate and template.
	ol_tail = get_tail(t, full_tail_len, ori);

	shift = ori ? 0 : ht->o->read_len - q_tail->len;
	mates = align_tpl_tail(ht, t, q_tail, shift, mismatches, FRESH, ori);
	for (i = 0; i < mates->len; i++) {
		m = (bwa_seq_t*) g_ptr_array_index(mates, i);
		r = get_mate(m, seqs);
		// Keep those reads whose mate is used by current template
		if (!r || r->contig_id != t->id || r->status != USED || is_biased_q(m)) {
			reset_to_fresh(m);
			continue;
		}
		// Find the overlapping between mate and tail
		//p_query(__func__, m);
		ol = find_fr_ol_within_k(m, ol_tail, mismatches, tail_len - 1, full_tail_len
				- 1, ori, &rev_com, &n_mis);
		//show_debug_msg("TAG", "OL: %d\n", ol);
		if (ol >= ht->o->k && ol >= n_mis * ht->o->k) {
			m->rev_com = rev_com;
			m->cursor = ori ? (m->len - ol - 1) : ol;
			m->pos = n_mis;
			//p_query(__func__, m);
			add2pool(p, m);
			added = 1;
		} else {
			reset_to_fresh(m);
		}
	}
	// With even shorter overlap and less mismatches allow.
	// Base by base overlapping.
	if (!added) {
		find_match_mates(ht, p, t, ht->o->read_len - 1, mismatches, ori);
	}
	g_ptr_array_free(mates, TRUE);
	bwa_free_read_seq(1, q_tail);
	bwa_free_read_seq(1, ol_tail);
}
