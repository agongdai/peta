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

comp *new_comp() {
	comp *c = (comp*) malloc(sizeof(comp));
	c->contigs = g_ptr_array_sized_new(4);
	c->edges = g_ptr_array_sized_new(4);
	c->hits = g_ptr_array_sized_new(4);
	return c;
}

void free_comp(comp *c) {
	g_ptr_array_free(c->contigs, TRUE);
	g_ptr_array_free(c->edges, TRUE);
	g_ptr_array_free(c->hits, TRUE);
}

GPtrArray *get_components(edgearray *all_edges, char *psl_name) {
	GPtrArray *hits = NULL;
	GPtrArray *all_comps = NULL;
	comp *new_c = NULL;
	int i = 0, eg_id = 0, j = 0, has_more = 1;
	blat_hit *h = NULL;
	edge *eg = NULL, *eg_q = NULL, *eg_t = NULL;

	all_comps = g_ptr_array_sized_new(16);
	show_debug_msg(__func__, "Reading Blat hits...\n");
	hits = read_blat_hits(psl_name);

	show_debug_msg(__func__, "Determining components...\n");
	for (i = 0; i < all_edges->len; i++) {
		eg = g_ptr_array_index(all_edges, i);
		if (eg->visited)
			continue;
		new_c = new_comp();
		g_ptr_array_add(new_c->edges, eg);
		has_more = 1;
		eg->visited = 1;
		while (has_more) {
			has_more = 0;
			for (j = 0; j < hits->len; j++) {
				h = g_ptr_array_index(hits, j);
				if (h->visited)
					continue;
				eg_id = atoi(h->qname);
				show_debug_msg(__func__, "Query Edge %d \n", eg_id);
				eg_q = g_ptr_array_index(all_edges, eg_id);
				eg_id = atoi(h->tname);
				show_debug_msg(__func__, "Ref Edge %d \n", eg_id);
				eg_t = g_ptr_array_index(all_edges, eg_id);
				if (edgearray_find(new_c->edges, eg_q) != NOT_FOUND
						|| edgearray_find(new_c->edges, eg_t) != NOT_FOUND) {
					if (!eg_q->visited) {
						g_ptr_array_add(new_c->edges, eg_q);
						eg_q->visited = 1;
						has_more = 1;
					}
					if (!eg_t->visited) {
						g_ptr_array_add(new_c->edges, eg_t);
						eg_t->visited = 1;
						has_more = 1;
					}
					g_ptr_array_add(new_c->hits, h);
					h->visited = 1;
				}
			}
		}
		g_ptr_array_add(all_comps, new_c);
	}
	return all_comps;
}

void scaffolding(edgearray *all_edges, const int insert_size,
		const int sd_insert_size, hash_table *ht, const int n_threads,
		char *psl_name) {
	GPtrArray *all_comps = NULL;
	int i = 0, j = 0;
	blat_hit *h = NULL;
	comp *c = NULL;


	all_comps = get_components(all_edges, psl_name);
	for (i = 0; i < all_comps->len; i++) {
		c = g_ptr_array_index(all_comps, i);
		for (j = 0; j < c->hits->len; j++) {
			h = g_ptr_array_index(c->hits, j);
			p_hit(h);
		}
	}
}
