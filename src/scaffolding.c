/*
 * scaffolding.c
 *
 *  Created on: 04-Feb-2013
 *      Author: carl
 */
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <glib.h>
#include "scaffolding.h"
#include "pepath.h"
#include "hits.h"
#include "edgelist.h"
#include "utils.h"
#include "peseq.h"
#include "merge.h"

comp *new_comp() {
	comp *c = (comp*) malloc(sizeof(comp));
	c->id = 0;
	c->id = 0;
	c->contigs = g_ptr_array_sized_new(4);
	c->edges = g_ptr_array_sized_new(4);
	c->hits = g_ptr_array_sized_new(4);
	c->alive = 1;
	return c;
}

void free_comp(comp *c) {
	g_ptr_array_free(c->contigs, TRUE);
	g_ptr_array_free(c->edges, TRUE);
	g_ptr_array_free(c->hits, TRUE);
	free(c);
}

comp *merge_comps(comp *c, comp *to_merge) {
	int i = 0;
	edge *eg = NULL;
	if (c == to_merge)
		return c;
	show_debug_msg(__func__, "Merging component %d to %d \n",
			to_merge->comp_start, c->comp_start);
	for (i = 0; i < to_merge->edges->len; i++) {
		eg = g_ptr_array_index(to_merge->edges, i);
		eg->comp_id = c->id;
	}
	g_ptr_array_concat(c->edges, to_merge->edges);
	g_ptr_array_concat(c->hits, to_merge->hits);
	g_ptr_array_concat(c->contigs, to_merge->contigs);
	to_merge->alive = 0;
	return c;
}

/**
 * Assumption: the components are with ids: 0, 1, 2, 3, ...
 */
void concat_connected_comps(edgearray *all_edges, const int insert_size, GPtrArray *comps, edge *eg,
		bwa_seq_t *seqs) {
	int i = 0;
	edge *in_out = NULL;
	edgearray *probable_in_out = NULL;
	comp *this_c = NULL, *c = NULL;

	this_c = g_ptr_array_index(comps, eg->comp_id);
	if (!this_c->alive)
		return;
	probable_in_out = get_probable_in_out(all_edges, insert_size, eg, seqs);
	//show_debug_msg(__func__, "Probable in out: %d \n", probable_in_out->len);
	for (i = 0; i < probable_in_out->len; i++) {
		in_out = g_ptr_array_index(probable_in_out, i);
		c = g_ptr_array_index(comps, in_out->comp_id);
		if (c->alive && this_c != c) {
			merge_comps(this_c, c);
		}
	}
	g_ptr_array_free(probable_in_out, TRUE);
}

void merge_roadmap_comps(edgearray *all_edges, const int insert_size, GPtrArray *comps,
		bwa_seq_t *seqs) {
	comp *c = NULL;
	int i = 0, j = 0;
	edge *eg = NULL;

	for (i = 0; i < comps->len; i++) {
		c = g_ptr_array_index(comps, i);
		for (j = 0; j < c->edges->len; j++) {
			eg = g_ptr_array_index(c->edges, j);
			concat_connected_comps(all_edges, insert_size, comps, eg, seqs);
		}
	}
}

GPtrArray *get_components(edgearray *all_edges, char *psl_name) {
	GPtrArray *hits = NULL;
	GPtrArray *all_comps = NULL;
	comp *new_c = NULL;
	int i = 0, eg_id = 0, j = 0, has_more = 1, comp_id = 0, hit_start = 0;
	blat_hit *h = NULL;
	edge *eg = NULL, *eg_q = NULL, *eg_t = NULL;

	all_comps = g_ptr_array_sized_new(16);
	show_debug_msg(__func__, "Reading Blat hits...\n");
	hits = read_blat_hits(psl_name);
	g_ptr_array_sort(hits, (GCompareFunc) cmp_hit_by_qname);

	show_debug_msg(__func__, "Determining components...\n");
	for (i = 0; i < all_edges->len; i++) {
		eg = g_ptr_array_index(all_edges, i);
		if (eg->visited || !eg->alive)
			continue;
		new_c = new_comp();
		new_c->id = comp_id++;
		new_c->comp_start = eg->id;
		eg->comp_id = new_c->id;
		g_ptr_array_add(new_c->edges, eg);
		has_more = 1;
		eg->visited = 1;
		while (has_more) {
			has_more = 0;
			//show_debug_msg(__func__, "Start edge: [%d, %d]; hit_start: %d \n", eg->id, eg->len, hit_start);
			for (j = hit_start; j < hits->len; j++) {
				h = g_ptr_array_index(hits, j);
				if (h->visited)
					continue;
				eg_id = atoi(h->qname);
				if (eg_id == eg->id) {
					hit_start = j + 1;
				}
				eg_q = g_ptr_array_index(all_edges, eg_id);
				eg_id = atoi(h->tname);
				eg_t = g_ptr_array_index(all_edges, eg_id);
				if (edgearray_find(new_c->edges, eg_q) != NOT_FOUND
						|| edgearray_find(new_c->edges, eg_t) != NOT_FOUND) {
					// A->B and B->A, only add either one of the hits
					if (!eg_q->visited || !eg_t->visited) {
						g_ptr_array_add(new_c->hits, h);
						h->visited = 1;
					}
					if (!eg_q->visited) {
						g_ptr_array_add(new_c->edges, eg_q);
						eg_q->visited = 1;
						eg_q->comp_id = new_c->id;
						has_more = 1;
					}
					if (!eg_t->visited) {
						g_ptr_array_add(new_c->edges, eg_t);
						eg_t->visited = 1;
						eg_t->comp_id = new_c->id;
						has_more = 1;
					}
				}
			}
		}
		g_ptr_array_add(all_comps, new_c);
	}
	return all_comps;
}

int post_val_edge_thread(gpointer e, gpointer data) {
	int i = 0, j = 0, start = 0, end = 0;
	int max = 0, min = 1000000;
	int *base_cov = NULL;
	bwa_seq_t *r = NULL, *mate = NULL;
	edge *eg = (edge*) e;
	hash_table *ht = (hash_table*) data;
	bwa_seq_t *seqs = ht->seqs;

	base_cov = (int*) calloc(eg->len, sizeof(int));
	realign_reads_by_ht(ht, eg, MISMATCHES);
	show_debug_msg(__func__, "Edge [%d, %d], reads %d=>%d \n", eg->id, eg->len,
			eg->reads->len, eg->pairs->len);
	for (i = 0; i < eg->reads->len; i++) {
		r = g_ptr_array_index(eg->reads, i);
		mate = get_mate(r, seqs);
		// p_query(__func__, r);
		// p_query(__func__, mate);
		if (mate->contig_id == eg->id && mate->rev_com == r->rev_com) {
			start = (r->shift > mate->shift) ? mate->shift : r->shift;
			end = start + abs(r->shift - mate->shift) + r->len;
			end = (end > eg->len) ? eg->len : end;
			for (j = start; j < end; j++) {
				base_cov[j]++;
			}
		}
	}
	if (eg->len > 2 * seqs->len) {
		for (i = seqs->len; i < eg->len - seqs->len; i++) {
			max = (max < base_cov[i]) ? base_cov[i] : max;
			min = (min > base_cov[i]) ? base_cov[i] : min;
			show_debug_msg(__func__, "%d\t%d\n", i, base_cov[i]);
		}
	}
	show_debug_msg(__func__, "[Min, Max]: [%d, %d] \n", min, max);
	if (min <= 2) {
		eg->alive = 0;
		free(base_cov);
		return 0;
	}
	eg->alive = 1;
	free(base_cov);
	return 1;
}

void post_validation(edgearray *all_edges, hash_table *ht, const int n_threads) {
	GThreadPool *thread_pool = NULL;
	int i = 0, ori_len = all_edges->len;
	edge *eg = NULL;
	thread_pool = g_thread_pool_new((GFunc) post_val_edge_thread, ht,
			n_threads, TRUE, NULL);
	for (i = 0; i < all_edges->len; i++) {
		eg = g_ptr_array_index(all_edges, i);
		g_thread_pool_push(thread_pool, (gpointer) eg, NULL);
	}
	g_thread_pool_free(thread_pool, 0, 1);
	for (i = 0; i < all_edges->len; i++) {
		eg = g_ptr_array_index(all_edges, i);
		if (!eg->alive) {
			show_debug_msg(__func__, "Edge %d,%d is not alive \n", eg->id,
					eg->len);
			g_ptr_array_remove_index_fast(all_edges, i);
			i--;
		}
	}
	show_msg(__func__, "Edge validation shrinks edges from %d to %d \n",
			ori_len, all_edges->len);
}

void scaffolding(edgearray *all_edges, const int insert_size,
		const int sd_insert_size, hash_table *ht, const int n_threads,
		char *psl_name) {
	GPtrArray *all_comps = NULL;
	int i = 0, j = 0;
	blat_hit *h = NULL;
	comp *c = NULL;
	edge *eg = NULL;

	all_comps = get_components(all_edges, psl_name);
	for (i = 0; i < all_edges->len; i++) {
		eg = g_ptr_array_index(all_edges, i);
		g_ptr_array_sort(eg->reads, (GCompareFunc) cmp_read_by_name);
	}
	merge_roadmap_comps(all_edges, insert_size, all_comps, ht->seqs);
	for (i = 0; i < all_comps->len; i++) {
		c = g_ptr_array_index(all_comps, i);
//		if (!c->alive)
//			continue;
		show_debug_msg(
				__func__,
				"======== %d: Hits for component %d ========= Alive: %d =========\n",
				c->id, c->comp_start, c->alive);
		printf("Edges: ");
		for (j = 0; j < c->edges->len; j++) {
			eg = g_ptr_array_index(c->edges, j);
			printf("[%d, %d], ", eg->id, eg->len);
		}
		printf("\n");
		for (j = 0; j < c->hits->len; j++) {
			h = g_ptr_array_index(c->hits, j);
			p_hit(h);
		}
	}
}
