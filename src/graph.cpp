/*
 * graph.cpp
 *
 *  Created on: Jun 30, 2013
 *      Author: carl
 */

#include <unordered_map>
#include <stdint.h>
#include <glib.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>
#include "bwtaln.h"
#include "tpl.h"
#include "kmers.hpp"
#include "ass.hpp"
#include "rnaseq.h"
#include "pechar.h"
#include "utils.h"
#include "peseq.h"
#include "junction.hpp"
#include "hash.hpp"
#include "graph.hpp"
#include "path.hpp"

int edge_id = 0;
int vertex_id = 0;

/**
 * Get the reads on a sequence allowing N_MISMATCHES
 */
GPtrArray *reads_on_seq(bwa_seq_t *seq, hash_map *hm, const int n_mismatch) {
	int j = 0, i = 0, read_len = hm->o->read_len;
	bwa_seq_t *part = NULL, *r = NULL;
	GPtrArray *hits = NULL;
	GPtrArray *reads = g_ptr_array_sized_new(32);

	for (i = 0; i <= seq->len - read_len; i++) {
		part = new_seq(seq, read_len, i);
		hits = align_full_seq(part, hm, n_mismatch);
		for (j = 0; j < hits->len; j++) {
			r = (bwa_seq_t*) g_ptr_array_index(hits, j);
			g_ptr_array_add(reads, r);
		}
		bwa_free_read_seq(1, part);
		g_ptr_array_free(hits, TRUE);
	}
	return reads;
}

splice_graph *new_graph() {
	splice_graph *g = (splice_graph*) malloc(sizeof(splice_graph));
	g->edges = g_ptr_array_sized_new(0);
	g->vertexes = g_ptr_array_sized_new(0);
	g->len = 0;
	return g;
}

gint cmp_junc_by_locus(gpointer a, gpointer b) {
	junction *c_a = *((junction**) a);
	junction *c_b = *((junction**) b);
	return ((c_a->locus) - c_b->locus);
}

gint cmp_junc_by_branch_id(gpointer a, gpointer b) {
	junction *c_a = *((junction**) a);
	junction *c_b = *((junction**) b);
	return ((c_a->branch_tpl->id) - c_b->branch_tpl->id);
}

vertex *new_vertex(tpl *t, int start, int len, hash_map *hm) {
	vertex *v = (vertex*) malloc(sizeof(vertex));
	v->ctg = new_seq(t->ctg, len, start);
	v->len = len;
	v->ins = g_ptr_array_sized_new(0);
	v->outs = g_ptr_array_sized_new(0);
	v->reads = reads_on_seq(v->ctg, hm, N_MISMATCHES);
	v->weight = (float) v->reads->len;
	v->status = 0;
	v->id = vertex_id++;
	return v;
}

void p_vertex(vertex *v) {
	show_debug_msg(__func__, "==== Vertex %d: %.2f ====\n", v->id, v->weight);
	p_ctg_seq(__func__, v->ctg);
	//p_readarray(v->reads, 1);
}

void p_edge(edge *e) {
	show_debug_msg(__func__, "---- Edge %d: %.2f ----\n", e->id, e->weight);
	show_debug_msg(__func__, "Left len: %d; right len: %d\n", e->left_len,
			e->right_len);
	p_ctg_seq(__func__, e->junc_seq);
	//p_readarray(e->reads, 1);
}

void p_graph(splice_graph *g) {
	edge *e = NULL;
	vertex *v = NULL, *left = NULL, *right = NULL;
	int i = 0, j = 0;
	GPtrArray *roots = NULL;
	char entry[BUFSIZ];
	FILE *dot = xopen("graph.dot", "w");
	fputs("digraph G {\n", dot);
	fputs("graph [rankdir=LR];", dot);
	for (i = 0; i < g->vertexes->len; i++) {
		v = (vertex*) g_ptr_array_index(g->vertexes, i);
		//		if (v->id == 1 || v->id == 2 || v->id == 3 || v->id == 5 || v->id == 7
		//				|| v->id == 8 || v->id == 10)
		//			sprintf(entry, "%d [label=\"%d: %d\" shape=box color=blue]; \n",
		//					v->id, v->id, v->len);
		//		else
		sprintf(entry, "%d [label=\"%d: %d\" shape=box]; \n", v->id, v->id,
				v->len);
		fputs(entry, dot);
	}
	for (j = 0; j < g->edges->len; j++) {
		e = (edge*) g_ptr_array_index(g->edges, j);
		left = e->left;
		right = e->right;
		sprintf(entry, "%d -> %d [label=\"%.0f\"]; \n", left->id, right->id,
				e->weight);
		fputs(entry, dot);
	}
	fputs("}\n", dot);
	fclose(dot);
}

void save_vertexes(GPtrArray *vertexes) {
	int i = 0;
	vertex *v = NULL;
	char entry[BUFSIZ];
	FILE *v_fp = xopen("vertexes.fa", "w");
	for (i = 0; i < vertexes->len; i++) {
		v = (vertex*) g_ptr_array_index(vertexes, i);
		sprintf(entry, ">%d length: %d\n", v->id, v->len);
		save_con(entry, v->ctg, v_fp);
	}
	fclose(v_fp);
}

edge *new_edge(vertex *left, vertex *right) {
	edge *e = (edge*) malloc(sizeof(edge));
	e->junc_seq = NULL;
	e->left = left;
	e->right = right;
	e->weight = -1;
	e->reads = NULL;
	e->left_len = e->right_len = 0;
	e->len = 0;
	e->status = 0;
	e->id = edge_id++;
	return e;
}

void destroy_edge(edge *eg) {
	if (eg) {
		bwa_free_read_seq(1, eg->junc_seq);
		g_ptr_array_free(eg->reads, TRUE);
		free(eg);
	}
}
void destroy_vertex(vertex *v) {
	if (v) {
		bwa_free_read_seq(1, v->ctg);
		g_ptr_array_free(v->ins, TRUE);
		g_ptr_array_free(v->outs, TRUE);
		g_ptr_array_free(v->reads, TRUE);
		free(v);
	}
}

void destroy_graph(splice_graph *g) {
	int i = 0;
	vertex *v = NULL;
	edge *e = NULL;
	if (g) {
		if (g->vertexes) {
			for (i = 0; i < g->vertexes->len; i++) {
				v = (vertex*) g_ptr_array_index(g->vertexes, i);
				destroy_vertex(v);
			}
			g_ptr_array_free(g->vertexes, TRUE);
		}
		if (g->edges) {
			for (i = 0; i < g->edges->len; i++) {
				e = (edge*) g_ptr_array_index(g->edges, i);
				destroy_edge(e);
			}
			g_ptr_array_free(g->edges, TRUE);
		}
		free(g);
	}
}

/**
 * Get junctions with the template t as main/branch
 */
GPtrArray *tpl_junctions(tpl *t, GPtrArray *all_juncs) {
	junction *junc = NULL;
	uint64_t i = 0, j = 0;
	GPtrArray *tpl_junc = NULL;
	tpl_junc = g_ptr_array_sized_new(4);
	for (j = 0; j < all_juncs->len; j++) {
		junc = (junction*) g_ptr_array_index(all_juncs, j);
		if (junc->main_tpl == t && !(junc->branch_tpl == t))
			g_ptr_array_add(tpl_junc, junc);
	}
	g_ptr_array_sort(tpl_junc, (GCompareFunc) cmp_junc_by_locus);
	return tpl_junc;
}

/**
 * Break a template into a set of edges and vertexes
 */
void break_tpl(tpl *t, GPtrArray *main_juncs, splice_graph *g, hash_map *hm) {
	junction *j = NULL;
	tpl *branch = NULL;
	vertex *v = NULL, *left = NULL, *right = NULL;
	edge *e = NULL;
	int i = 0, pre_start = 0;
	GPtrArray *branch_vs = NULL, *branch_juncs = NULL;
	int max_len = (hm->o->read_len - SHORT_BRANCH_SHIFT) * 2;
	GPtrArray *this_vs = g_ptr_array_sized_new(main_juncs->len + 1);

	// Break the main template into vertexes
	for (i = 0; i < main_juncs->len; i++) {
		j = (junction*) g_ptr_array_index(main_juncs, i);
		if (j->main_tpl != t || (j->main_tpl == t && j->branch_tpl == t)
				|| j->locus == pre_start)
			continue;
		v = new_vertex(t, pre_start, j->locus - pre_start, hm);
		g_ptr_array_add(g->vertexes, v);
		g_ptr_array_add(this_vs, v);
		pre_start = j->locus;
	}
	// Create the last vertex
	v = new_vertex(t, pre_start, t->len - pre_start, hm);
	g_ptr_array_add(g->vertexes, v);
	g_ptr_array_add(this_vs, v);

	show_debug_msg(__func__, "Template %d, Vertex %d \n", t->id, v->id);

	// Create the edges between vertexes
	pre_start = 0;
	for (i = 0; i < this_vs->len - 1; i++) {
		left = (vertex*) g_ptr_array_index(this_vs, i);
		right = (vertex*) g_ptr_array_index(this_vs, i + 1);
		e = new_edge(left, right);
		e->junc_seq = get_junc_seq(t, pre_start, &e->left_len, t, pre_start,
				&e->right_len, max_len);
		e->len = e->left_len + e->right_len;
		e->reads = reads_on_seq(e->junc_seq, hm, N_MISMATCHES);
		e->weight = (float) e->reads->len;
		g_ptr_array_add(left->outs, e);
		g_ptr_array_add(right->ins, e);
		g_ptr_array_add(g->edges, e);
		pre_start += left->len;
	}

	t->vertexes = this_vs;
}

/**
 * Find the vertex with some locus;
 * If ori 0: return the right vertex of locus
 * If ori 1: return the left vertex of locus
 */
vertex *find_vertex_by_locus(tpl *t, int locus, int ori) {
	int sum_shift = 0;
	uint32_t i = 0;
	vertex *this_v = NULL, *next_v = NULL;
	if (!t || !t->vertexes)
		return NULL;
	if (t->vertexes->len == 1) {
		this_v = (vertex*) g_ptr_array_index(t->vertexes, 0);
		if (locus == t->len && ori == 1)
			return this_v;
		else if (locus == 0 && ori == 0)
			return this_v;
		else
			return NULL;
	}
	for (i = 0; i < t->vertexes->len - 1; i++) {
		this_v = (vertex*) g_ptr_array_index(t->vertexes, i);
		next_v = (vertex*) g_ptr_array_index(t->vertexes, i + 1);
		sum_shift += this_v->len;
		if (locus == sum_shift) {
			if (ori == 1)
				return this_v;
			else
				return next_v;
		}
	}
	return NULL;
}

// Special cases:
// 	Main		Branch	Locus	Weight	Direction
//	[1, 1704]	[2, 0]	802		128		1
//	[1, 1704]	[2, 0]	491		98		0
void remove_zero_vertexes(tpl *t, GPtrArray *main_juncs, int read_len,
		splice_graph *g) {
	vertex *left = NULL, *right = NULL;
	edge *e = NULL;
	uint32_t i = 0;
	int left_locus = 0, right_locus = 0;
	junction *j = NULL, *next_j = NULL;
	int max_len = (read_len - SHORT_BRANCH_SHIFT) * 2;
	GPtrArray *zero_juncs = g_ptr_array_sized_new(main_juncs->len);
	// Pick those junctions with branch template length 0
	for (i = 0; i < main_juncs->len; i++) {
		j = (junction*) g_ptr_array_index(main_juncs, i);
		if (j->branch_tpl->len == 0) {
			g_ptr_array_add(zero_juncs, j);
			g_ptr_array_remove_index_fast(main_juncs, i--);
		}
	}
	if (zero_juncs->len >= 2) {
		g_ptr_array_sort(zero_juncs, (GCompareFunc) cmp_junc_by_branch_id);
		for (i = 0; i < zero_juncs->len - 1; i++) {
			j = (junction*) g_ptr_array_index(zero_juncs, i);
			next_j = (junction*) g_ptr_array_index(zero_juncs, i + 1);
			if (j->branch_tpl == next_j->branch_tpl) {
				left_locus = min(j->locus, next_j->locus);
				right_locus = max(j->locus, next_j->locus);
				if (right_locus - left_locus > 0) {
					left = find_vertex_by_locus(t, left_locus, 1);
					right = find_vertex_by_locus(t, right_locus, 0);
					e = new_edge(left, right);
					e->junc_seq = get_junc_seq(t, left_locus, &e->left_len, t,
							right_locus, &e->right_len, max_len);
					g_ptr_array_add(left->outs, e);
					g_ptr_array_add(right->ins, e);
					g_ptr_array_add(g->edges, e);
					show_debug_msg(__func__,
							"New edge between vertexes %d and %d \n", left->id,
							right->id);
				}
			}
		}
	}
	g_ptr_array_free(zero_juncs, TRUE);
}

/**
 * Add edges to connect different vertexes on templates
 */
void connect_tpls(tpl *t, GPtrArray *main_juncs, splice_graph *g, hash_map *hm) {
	junction *j = NULL;
	int i = 0, n_zeros = 0, nth_zero = 0;
	int max_len = (hm->o->read_len - SHORT_BRANCH_SHIFT) * 2;
	edge *e = NULL;
	vertex *branch_v = NULL, *left = NULL, *right = NULL;
	GPtrArray *branch_vs = NULL, *main_vs = NULL;
	// Special cases:
	// 	Main		Branch	Locus	Weight	Direction
	//	[1, 3842]	[2, 66]	736		128		1
	//	[1, 3842]	[2, 66]	736		124		0
	//
	//	[1, 1704]	[2, 0]	802		128		1
	//	[1, 1704]	[2, 0]	491		98		0
	for (i = 0; i < main_juncs->len; i++) {
		j = (junction*) g_ptr_array_index(main_juncs, i);
		if (j->main_tpl != t || (j->main_tpl == t && j->branch_tpl == t))
			continue;

		show_debug_msg(__func__, "Junction %d: nth_zero: %d \n", i, nth_zero);
		p_junction(j);

		branch_vs = j->branch_tpl->vertexes;
		main_vs = j->main_tpl->vertexes;
		if (j->ori) { // Branch is left, main is right
			branch_v = (vertex*) g_ptr_array_index(branch_vs, branch_vs->len
					- 1);
			// Shifting on the main_tpl to find connecting vertex
			right = find_vertex_by_locus(t, j->locus, 0);
			e = new_edge(branch_v, right);
			e->junc_seq = get_junc_seq(j->branch_tpl, j->branch_tpl->len,
					&e->left_len, t, j->locus, &e->right_len, max_len);
			g_ptr_array_add(branch_v->outs, e);
			g_ptr_array_add(right->ins, e);
		} else { // Branch is right, main is left
			branch_v = (vertex*) g_ptr_array_index(branch_vs, 0);
			// Shifting on the main_tpl to find connecting vertex
			left = find_vertex_by_locus(t, j->locus, 1);
			e = new_edge(left, branch_v);
			e->junc_seq = get_junc_seq(t, j->locus, &e->left_len,
					j->branch_tpl, 0, &e->right_len, max_len);
			g_ptr_array_add(left->outs, e);
			g_ptr_array_add(branch_v->ins, e);
		}
		e->len = e->left_len + e->right_len;
		// The weight is not set because of probable alternative paths.
		//e->reads = reads_on_seq(e->junc_seq, hm);
		//e->weight = (float) e->reads->len;
		g_ptr_array_add(g->edges, e);
	}
}

void p_tpl_juncs(tpl *t, GPtrArray *t_juncs) {
	show_debug_msg(__func__, "==== Junctions of Template %d ====\n", t->id);
	uint32_t i = 0;
	junction *j = NULL;
	vertex *v = NULL;
	for (i = 0; i < t_juncs->len; i++) {
		j = (junction*) g_ptr_array_index(t_juncs, i);
		p_junction(j);
	}
	show_debug_msg(__func__, "==== Vertexes of Template %d ====\n", t->id);
	for (i = 0; i < t->vertexes->len; i++) {
		v = (vertex*) g_ptr_array_index(t->vertexes, i);
		printf("\tVertex %d: %d \n", v->id, v->len);
	}
}

splice_graph *build_graph(GPtrArray *all_tpls, GPtrArray *all_juncs,
		hash_map *hm) {
	splice_graph *g = new_graph();
	GPtrArray *t_juncs = NULL;
	vertex *v = NULL;
	int i = 0;
	tpl *t = NULL;
	for (i = 0; i < all_tpls->len; i++) {
		t = (tpl*) g_ptr_array_index(all_tpls, i);
		t_juncs = tpl_junctions(t, all_juncs);
		break_tpl(t, t_juncs, g, hm);
		p_tpl_juncs(t, t_juncs);
		g_ptr_array_free(t_juncs, TRUE);
	}

	for (i = 0; i < all_tpls->len; i++) {
		t = (tpl*) g_ptr_array_index(all_tpls, i);
		t_juncs = tpl_junctions(t, all_juncs);
		printf("\n\n\n==== Building subgraph connected to template %d ===\n ",
				t->id);
		remove_zero_vertexes(t, t_juncs, hm->o->read_len, g);
		connect_tpls(t, t_juncs, g, hm);
		p_tpl_juncs(t, t_juncs);
		g_ptr_array_free(t_juncs, TRUE);
	}

	g->len = 0;
	for (i = 0; i < g->vertexes->len; i++) {
		v = (vertex*) g_ptr_array_index(g->vertexes, i);
		g->len += v->len;
	}
	return g;
}

/**
 * Remove isolated vertexes with length smaller than 10
 */
void clean_graph(splice_graph *g) {
	vertex *v = NULL;
	uint32_t i = 0;
	for (i = 0; i < g->vertexes->len; i++) {
		v = (vertex*) g_ptr_array_index(g->vertexes, i);
		if (v->len < 10 && v->ins->len == 0 && v->outs->len == 0)
			g_ptr_array_remove_index_fast(g->vertexes, i--);
	}
}

void process_graph(GPtrArray *all_tpls, GPtrArray *all_juncs, hash_map *hm) {
	splice_graph *g = NULL;
	g = build_graph(all_tpls, all_juncs, hm);
	clean_graph(g);
	p_graph(g);
	save_vertexes(g->vertexes);
	determine_paths(g, hm);
}