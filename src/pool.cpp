/*
 * pool.c
 *
 *  Created on: Jul 23, 2013
 *      Author: carl
 */
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "pool.h"
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
		for (i = 0; i < p->reads->len; i++) {
			r = (bwa_seq_t*) g_ptr_array_index(p->reads, i);
			r->pos = -1;
			r->cursor = -1;
			if (r->status == IN_POOL)
				r->status = FRESH;
		}
		g_ptr_array_free(r->reads, TRUE);
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
			printf("%d->%c@%s\t[status: %d]\t[rev_com]", p2->cursor, c,
					p2->name, p->status);
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
			printf("%d->%c@%s\t[status: %d]\t[forward]", p2->cursor, c,
					p2->name, p->status);
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
}

/**
 * Count frequencies of next characters, get the most frequent one.
 * If -1: all of them are 0
 */
int get_next_char(pool *p, const int ori, tpl *t) {
	readarray *reads = p->reads;
	int *c = (int*) calloc(5, sizeof(int));
	int i = 0, next_char = -1, max_c = 0;
	bwa_seq_t *r;
	c[0] = c[1] = c[2] = c[3] = c[4] = 0;
	for (i = 0; i < reads->len; i++) {
		r = g_ptr_array_index(reads, i);
		if (r->rev_com)
			c[r->rseq[r->cursor]]++;
		else
			c[r->seq[r->cursor]]++;
	}
	for (i = 0; i < 5; i++) {
		if (c[i] > max_c) {
			max_c = c[i];
			next_char = i;
		}
	}
	if (max_c == 0)
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
			locus = ori ? 0 : (t->len - r->len);
			add2tpl(t, r, locus);
			g_ptr_array_remove_index_fast(reads, i--);
		}
	}
}

/**
 * Align the tail to the RNA-seq reads, add fresh reads to current pool
 */
void next_pool(pool *p, tpl *t, hash_table *ht, bwa_seq_t *tail,
		int mismatches, const int ori) {
	GPtrArray *fresh_reads = NULL;
	index64 i = 0;
	bwa_seq_t *r = NULL;

	// Find all fresh reads and add to the pool
	fresh_reads = align_tpl_tail(ht, t, tail, mismatches, FRESH, ori);
	for (i = 0; i < fresh_reads->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(fresh_reads, i);
		add2pool(p, r);
	}
	g_ptr_array_free(fresh_reads, TRUE);
}
