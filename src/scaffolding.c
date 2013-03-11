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
#include "readrm.h"

GMutex *comp_mutex = NULL;

gint cmp_break_point(gpointer a, gpointer b) {
	int* x = (int*) a;
	int* y = (int*) b;
	return *y - *x;
}

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

void p_comps(GPtrArray *all_comps) {
	int i = 0, j = 0;
	comp *c = NULL;
	edge *eg = NULL;
	blat_hit *h = NULL;
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

comp *combine_two_comps(comp *c, comp *to_merge) {
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
 * Combine two components if there are paired reads spanning them
 * Assumption: the components are with ids: 0, 1, 2, 3, ...
 */
void combine_connected_comps_thread(gpointer e, gpointer data) {
	int i = 0;
	edge *in_out = NULL, *eg = (edge*) e;
	edgearray *probable_in_out = NULL;
	comp *this_c = NULL, *c = NULL;
	comps_aux_t *d = (comps_aux_t*) data;

	this_c = g_ptr_array_index(d->comps, eg->comp_id);
	if (!this_c->alive || !eg->alive)
		return;
	probable_in_out = get_probable_in_out(d->all_edges, d->insert_size,
			d->sd_insert_size, eg, d->ht->seqs, 1);
	if (probable_in_out->len >= 1) {
		show_debug_msg(__func__, ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");
		show_debug_msg(__func__, "Probable in out of edge %d: %d \n", eg->id,
				probable_in_out->len);
	}
	// If any probable in-out edges are used by another thread, abondan it.
	g_mutex_lock(comp_mutex);
	for (i = 0; i < probable_in_out->len; i++) {
		in_out = g_ptr_array_index(probable_in_out, i);
		if (in_out->level != eg->level) {
			g_mutex_unlock(comp_mutex);
			g_ptr_array_free(probable_in_out, TRUE);
			return;
		}
	}
	g_mutex_unlock(comp_mutex);
	for (i = 0; i < probable_in_out->len; i++) {
		in_out = g_ptr_array_index(probable_in_out, i);
		show_debug_msg(__func__, "\tEdge [%d, %d]\n", in_out->id, in_out->len);
		c = g_ptr_array_index(d->comps, in_out->comp_id);
		if (this_c->alive && c->alive && this_c != c) {
			combine_two_comps(this_c, c);
		}
	}
	if (probable_in_out->len >= 1) {
		show_debug_msg(__func__, "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n");
	}
	g_ptr_array_free(probable_in_out, TRUE);
}

void combine_roadmap_comps(edgearray *all_edges, const int insert_size,
		const int sd_insert_size, GPtrArray *comps, hash_table *ht,
		const int n_threads) {
	comp *c = NULL;
	int i = 0, j = 0;
	edge *eg = NULL;
	GThreadPool *thread_pool = NULL;
	comps_aux_t *data = (comps_aux_t*) malloc(sizeof(comps_aux_t));

	if (!comp_mutex)
		comp_mutex = g_mutex_new();
	data->ht = ht;
	data->insert_size = insert_size;
	data->sd_insert_size = sd_insert_size;
	data->all_edges = all_edges;
	data->comps = comps;

	// A flag is set to avoid recursive calling of get_probable in-out edges
	for (i = 0; i < all_edges->len; i++) {
		eg = g_ptr_array_index(all_edges, i);
		eg->level = eg->comp_id;
	}
	thread_pool = g_thread_pool_new((GFunc) combine_connected_comps_thread,
			data, n_threads, TRUE, NULL);
	for (i = 0; i < comps->len; i++) {
		c = g_ptr_array_index(comps, i);
		if (!c->alive)
			continue;
		for (j = 0; j < c->edges->len; j++) {
			eg = g_ptr_array_index(c->edges, j);
			g_thread_pool_push(thread_pool, (gpointer) eg, NULL);
		}
	}
	g_thread_pool_free(thread_pool, 0, 1);
	for (i = 0; i < all_edges->len; i++) {
		eg = g_ptr_array_index(all_edges, i);
		eg->level = 0;
	}
	free(data);
}

int merge_comps_egs_thread(gpointer c, gpointer data) {
	comp *com = (comp*) c;
	edgearray *comp_edges = com->edges;
	comps_aux_t *d = (comps_aux_t*) data;
	show_debug_msg(__func__,
			"Merging overlapped edges in separated component %d \n");
	merge_ol_comp_edges(comp_edges, d->ht, d->insert_size, d->sd_insert_size);
	return 0;
}

void merge_comps_egs(edgearray *all_edges, GPtrArray *comps, hash_table *ht,
		const int insert_size, const int sd_insert_size, const int n_threads) {
	comp *c = NULL;
	int i = 0;
	GThreadPool *thread_pool = NULL;
	comps_aux_t *data = (comps_aux_t*) malloc(sizeof(comps_aux_t));
	data->ht = ht;
	data->insert_size = insert_size;
	data->sd_insert_size = sd_insert_size;
	data->all_edges = NULL;
	data->comps = NULL;

	for (i = 0; i < comps->len; i++) {
		c = g_ptr_array_index(comps, i);
		if (!c->alive) {
			g_ptr_array_remove_index_fast(comps, i);
			i--;
		}
	}

	thread_pool = g_thread_pool_new((GFunc) merge_comps_egs_thread, data,
				n_threads, TRUE, NULL);

	for (i = 0; i < comps->len; i++) {
		c = g_ptr_array_index(comps, i);
		if (c->edges->len > 1)
			g_thread_pool_push(thread_pool, (gpointer) c, NULL);
	}
	g_thread_pool_free(thread_pool, 0, 1);
	free(data);
}

/**
 * Converge overlapped edges into separated components
 */
GPtrArray *get_components(edgearray *all_edges, char *psl_name) {
	GPtrArray *hits = NULL;
	GPtrArray *all_comps = NULL;
	comp *new_c = NULL;
	int i = 0, eg_id = 0, j = 0, has_more = 1, comp_id = 0, hit_start = 0,
			no_more_hits = 0;
	blat_hit *h = NULL;
	edge *eg = NULL, *eg_q = NULL, *eg_t = NULL;

	all_comps = g_ptr_array_sized_new(16);
	show_debug_msg(__func__, "Reading Blat hits...\n");
	hits = read_blat_hits(psl_name);
	g_ptr_array_sort(hits, (GCompareFunc) cmp_hit_by_qname);

	for (i = 0; i < all_edges->len; i++) {
		eg = g_ptr_array_index(all_edges, i);
		eg->visited = 0;
	}
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
			no_more_hits = 1;
			for (j = hit_start; j < hits->len; j++) {
				h = g_ptr_array_index(hits, j);
				if (!h->visited) {
					no_more_hits = 0;
					break;
				}
			}
			if (no_more_hits)
				break;
			//show_debug_msg(__func__, "Start edge: [%d, %d]; hit_start: %d \n", eg->id, eg->len, hit_start);
			for (j = hit_start; j < hits->len; j++) {
				h = g_ptr_array_index(hits, j);
				if (h->visited)
					continue;
				// If the overlapped length is too short, ignore this hit
				if (h->matches < MIN_OL_TO_SCA)
					continue;
				eg_id = atoi(h->qname);
				if (eg_id == eg->id) {
					hit_start = j + 1;
				}
				eg_q = g_ptr_array_index(all_edges, eg_id);
				eg_id = atoi(h->tname);
				eg_t = g_ptr_array_index(all_edges, eg_id);
				if (!eg_q->alive || !eg_t->alive) {
					h->visited = 1;
					continue;
				}
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
	for (i = 0; i < all_edges->len; i++) {
		eg = g_ptr_array_index(all_edges, i);
		eg->visited = 0;
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

void order_comp_egs(edgearray *con_egs, bwa_seq_t *seqs) {
	int i = 0, j = 0, order = 0;
	edge *eg_i = NULL, *eg_j = 0;

	for (i = 0; i < con_egs->len; i++) {
		eg_i = g_ptr_array_index(con_egs, i);
		for (j = i + 1; j < con_egs->len; j++) {
			eg_j = g_ptr_array_index(con_egs, j);
			show_debug_msg(__func__, "Getting Edge %d and Edge %d order\n",
					eg_i->id, eg_j->id);
			order = order_two_edges(eg_i, eg_j, seqs);
			if (order == 0)
				continue;
			show_debug_msg(__func__, "Edge %d and Edge %d order: %d \n",
					eg_i->id, eg_j->id, order);
			if (order == 1) {
				g_ptr_array_add(eg_i->out_egs, eg_j);
				g_ptr_array_add(eg_j->in_egs, eg_i);
			} else {
				g_ptr_array_add(eg_i->in_egs, eg_j);
				g_ptr_array_add(eg_j->out_egs, eg_i);
			}
		}
	}
}

void get_break_points(edgearray *con_egs, GPtrArray *hits) {
	int i = 0, j = 0, k = 0, n_break_point = 0, point = 0;
	edge *eg = NULL;
	blat_hit *h = NULL;
	int *starts = NULL, *sizes = NULL;
	// For each edge, get its break points.
	// If an edge has 3 break points, it would be broken into 4 shorter edges.
	for (i = 0; i < con_egs->len; i++) {
		eg = g_ptr_array_index(con_egs, i);
		n_break_point = 0;
		starts = NULL;
		for (j = 0; j < hits->len; j++) {
			h = g_ptr_array_index(hits, j);
			if (h->alen < 100)
				continue;
			sizes = h->block_sizes;
			if (eg->id == atoi(h->qname))
				starts = h->q_starts;
			if (eg->id == atoi(h->tname))
				starts = h->t_starts;
			if (starts) {
				for (k = 0; k < h->block_count; k++) {
					g_array_append_val(eg->break_points, starts[k]);
				}
				if (h->q_size - h->q_end >= 10 && h->t_size - h->t_end >= 10) {
					point = starts[h->block_count - 1] + sizes[h->block_count
							- 1];
					g_array_append_val(eg->break_points, point);
				}
			}
		}
		g_array_sort(eg, cmp_break_point);
	}
}

void break_ol_egs(edgearray *con_egs, GPtrArray *hits, edgearray *all_edges) {
	int j = 0, i = 0;
	edge *eg = NULL;
	get_break_points(con_egs, hits);
	for (i = 0; i < con_egs->len; i++) {
		eg = g_ptr_array_index(con_egs, i);
		for (j = 0; j < eg->break_points->len; j++) {
			show_debug_msg(__func__, "Edge %d: %d \n", eg->id,
					g_array_index(eg->break_points, int, j));
		}
	}
}

int scaffold_comp_egs_thread(gpointer component, gpointer data) {
	comps_aux_t *d = (comps_aux_t*) data;
	edgearray *con_egs = NULL, *ol_egs = NULL;
	comp *c = (comp*) component;
	int i = 0;
	blat_hit *h = NULL;
	edge *eg = NULL;
	hash_table *ht = d->ht;

	ol_egs = g_ptr_array_sized_new(0);
	for (i = 0; i < c->hits->len; i++) {
		h = g_ptr_array_index(c->hits, i);
		eg = edgearray_find_id(c->edges, atoi(h->qname));
		if (!eg->alive)
			continue;
		if (eg)
			g_ptr_array_uni_add(ol_egs, eg);
		eg = edgearray_find_id(c->edges, atoi(h->tname));
		if (eg)
			g_ptr_array_uni_add(ol_egs, eg);
	}

	if (c->edges->len - ol_egs->len > 10) {
		g_ptr_array_free(ol_egs, TRUE);
		return 0;
	}

	con_egs = g_ptr_array_sized_new(0);
	for (i = 0; i < c->edges->len; i++) {
		eg = g_ptr_array_index(c->edges, i);
		if (!eg->alive)
			continue;
		if (edgearray_find(ol_egs, eg) == NOT_FOUND)
			g_ptr_array_uni_add(con_egs, eg);
	}

	order_comp_egs(con_egs, ht->seqs);
	//break_ol_egs(con_egs, c->hits, d->all_edges);
	g_ptr_array_free(ol_egs, TRUE);
	g_ptr_array_free(con_egs, TRUE);
	return 0;
}

void scaffold_comp_egs(GPtrArray *all_comps, edgearray *all_edges,
		hash_table *ht, const int n_threads) {
	GThreadPool *thread_pool = NULL;
	int i = 0;
	comp *c = NULL;
	edge *eg = NULL;
	comps_aux_t *data = (comps_aux_t*) malloc(sizeof(comps_aux_t));

	data->all_edges = all_edges;
	data->ht = ht;
	thread_pool = g_thread_pool_new((GFunc) scaffold_comp_egs_thread, data,
			n_threads, TRUE, NULL);
	for (i = 0; i < all_comps->len; i++) {
		c = g_ptr_array_index(all_comps, i);
		g_thread_pool_push(thread_pool, (gpointer) c, NULL);
	}
	g_thread_pool_free(thread_pool, 0, 1);
	for (i = 0; i < all_edges->len; i++) {
		eg = g_ptr_array_index(all_edges, i);
		eg->level = 0;
		if (!eg->in_egs || eg->in_egs->len == 0)
			eg->is_root = 1;
		else
			eg->is_root = 0;
	}
	free(data);
}

/**
 * Assumption: the components are with ids: 0, 1, 2, 3, ...
 */
edgearray *scaffolding(edgearray *all_edges, const int insert_size,
		const int sd_insert_size, hash_table *ht, const int n_threads,
		char *psl_name) {
	GPtrArray *all_comps = NULL;
	int i = 0;
	edge *eg = NULL;
	comp *c = NULL;

	FILE *merged_pair_contigs = NULL;

	all_comps = get_components(all_edges, psl_name);
	// Prepare the edges with refreshed status
	for (i = 0; i < all_edges->len; i++) {
		eg = g_ptr_array_index(all_edges, i);
		eg->tid = -1;
		g_ptr_array_sort(eg->reads, (GCompareFunc) cmp_read_by_name);
	}
	combine_roadmap_comps(all_edges, insert_size, sd_insert_size, all_comps,
			ht, n_threads);
	p_comps(all_comps);
	merge_comps_egs(all_edges, all_comps, ht, insert_size, sd_insert_size,
			n_threads);
	show_msg(__func__, "Scaffolding the edges in components... \n");
	scaffold_comp_egs(all_comps, all_edges, ht, n_threads);

	//	merged_pair_contigs = xopen("../SRR097897_out/validated.3.fa", "w");
	//	reset_edge_ids(all_edges);
	//	save_edges(all_edges, merged_pair_contigs, 0, 0, 0);

	for (i = 0; i < all_comps->len; i++) {
		c = g_ptr_array_index(all_comps, i);
		free_comp(c);
	}
	return all_edges;
}
