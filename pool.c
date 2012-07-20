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

/**
 * Assume that the new read has been verified not existed using exists().
 */
int get_insert_pos(const pool *r_pool, const bwa_seq_t *read) {
	unsigned int start = 0, end = r_pool->n - 1, middle = 0;
	int read_id, id;
	readarray *reads = r_pool->reads;
	bwa_seq_t *r;

	if (!r_pool->n || !read)
		return 0;
	read_id = atoi(read->name);
	r = g_ptr_array_index(reads, 0);
	id = atoi(r->name);
	if (read_id < id)
		return 0;
	r = g_ptr_array_index(reads, reads->len - 1);
	id = atoi(r->name);
	if (read_id > id)
		return reads->len;

	while (start <= end) {
		middle = (end + start) / 2;
		r = g_ptr_array_index(reads, middle);
		id = atoi(r->name);
		if (id == read_id)
			return -1;
		if (id < read_id)
			start = middle + 1;
		else
			end = middle - 1;
	}
	return end + 1;
}

void insert_fast_index(pool *r_pool, const int index, bwa_seq_t *read) {
	readarray *reads = r_pool->reads;
	g_ptr_array_add_index(reads, read, index);
	r_pool->n = reads->len;
}

/**
 * Add a new sequence to the pool
 */
void pool_sort_ins(pool *r_pool, bwa_seq_t *new_seq) {
	int index = binary_exists(r_pool, new_seq);
	if (index)
		return;
	index = get_insert_pos(r_pool, new_seq);
	assert(index != -1);
	insert_fast_index(r_pool, index, new_seq);
}

void pool_add(pool *p, bwa_seq_t *new_seq) {
	if (!new_seq)
		return;
	g_ptr_array_add(p->reads, new_seq);
	new_seq->is_in_c_pool = 1;
	p->n++;
}

void pool_uni_add(pool *p, bwa_seq_t *new_seq) {
	if (pool_exists(p, new_seq))
		return;
	pool_add(p, new_seq);
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
	read->rev_com = 0;
	p->n = p->reads->len;
	return r;
}

gboolean pool_rm_index(pool *p, const int i) {
	gboolean r;
	bwa_seq_t *read = g_ptr_array_index(p->reads, i);
	read->is_in_c_pool = 0;
	read->rev_com = 0;
	r = g_ptr_array_remove_index_fast(p->reads, i);
	p->n = p->reads->len;
	return r;
}

void syn_pools(pool *cur_pool, pool *mate_pool, const bwa_seq_t *seqs,
		const int ori) {
	int i = 0;
	bwa_seq_t *read, *mate;
	for (i = 0; i < cur_pool->n; i++) {
		read = g_ptr_array_index(cur_pool->reads, i);
		mate = get_mate(read, seqs);
		if (ori && is_right_mate(read->name)) {
			pool_uni_add(mate_pool, mate);
		}
		if (!ori && is_left_mate(read->name)) {
			pool_uni_add(mate_pool, mate);
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
			g_ptr_array_remove_index_fast(reads, i);
			i--;
			cur_pool->n--;
			continue;
		}
		if ((p->cursor < p->len) && p->cursor >= 0) {
			if (p->used) {
				p->is_in_c_pool = 0;
				if (!used)
					used = p;
				g_ptr_array_remove_index_fast(reads, i);
				i--;
				cur_pool->n--;
			}
		} else { // If the cursor is out of length, mark it as used!
			p->used = 1;
			p->shift = (ass_eg->len - p->cursor + 1);
			p->contig_id = ass_eg->id;
			p->is_in_c_pool = 0;
			pool_rm_index(cur_pool, i);
			g_ptr_array_add(ass_eg->reads, p);
			i--;
		}
	}
	return used;
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
			g_ptr_array_add(eg->reads, r);
			r->shift = (eg->len - r->cursor + 1);
			r->used = 1;
			r->contig_id = eg->id;
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
 * Assume the reads are sorted increasingly by read id
 */
int binary_exists(const pool *r_pool, const bwa_seq_t *read) {
	unsigned int start = 0, end = r_pool->n - 1, middle = 0;
	int read_id, id;
	readarray *reads = r_pool->reads;
	bwa_seq_t *r;

	if (!r_pool->n || !read)
		return 0;

	read_id = atoi(read->name);
	r = g_ptr_array_index(reads, 0);
	id = atoi(r->name);
	if (read_id < id)
		return 0;
	r = g_ptr_array_index(reads, reads->len - 1);
	id = atoi(r->name);
	if (read_id > id)
		return 0;

	// Binary search
	//	printf("[exists] Looking for %d \n", read_id);
	while (start <= end) {
		middle = (end + start) / 2;
		r = g_ptr_array_index(reads, middle);
		id = atoi(r->name);
		//		printf("[exists] id = %d, read_id = %d, [%d, %d, %d] \n", id, read_id, start,
		//				middle, end);
		if (id == read_id) {
			return middle + 1;
		} else {
			if (id < read_id)
				start = middle + 1;
			else
				end = middle - 1;
		}
	}
	return 0;
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
	}
	while (reads->len > 0)
		g_ptr_array_remove_index(r_pool->reads, 0);
	r_pool->reads->len = 0;
	r_pool->n = 0;
}

void free_pool(pool *r_pool) {
	readarray *reads = r_pool->reads;
	bwa_seq_t *r = 0;
	// p_pool("Pool to free: ", r_pool);
	int i = 0;
	// show_debug_msg(__func__, "Pool: %p \n", r_pool);
	if (r_pool) {
		if (r_pool->n > 0) {
			for (i = 0; i < r_pool->n; i++) {
				r = g_ptr_array_index(reads, i);
				r->is_in_c_pool = 0;
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
			printf("%d->%c\t%s\t%d\t[rev_com]", p2->cursor, c, p2->name,
					p2->used);
		} else {
			printf("%s", p2->seq);
			for (i = 0; i < p2->cursor + 2; i++)
				printf(" ");
			printf("%d->%c\t%s\t%d", p2->cursor, c, p2->name, p2->used);
		}
		printf("\n");
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
			show_debug_msg(__func__, "%d: %d_%d Read %s: %d \n", i,
					p->contig_id, p->shift, p->name, p->shift);
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

