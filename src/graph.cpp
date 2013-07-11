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

int comp_id = 0;
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
	g->stack = g_ptr_array_sized_new(0);
	g->scc = g_ptr_array_sized_new(0);
	g->components = g_ptr_array_sized_new(0);
	g->len = 0;
	g->index = 0;
	return g;
}

comp *new_comp() {
	comp *c = (comp*) malloc(sizeof(comp));
	c->vertexes = g_ptr_array_sized_new(4);
	c->edges = g_ptr_array_sized_new(4);
	c->stack = g_ptr_array_sized_new(0);
	c->scc = g_ptr_array_sized_new(0);
	c->id = ++comp_id;
	return c;
}

void destroy_comp(comp *c) {
	if (c) {
		if (c->vertexes)
			g_ptr_array_free(c->vertexes, TRUE);
		if (c->edges)
			g_ptr_array_free(c->edges, TRUE);
	}
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
	v->id = ++vertex_id;
	v->index = 0;
	v->lowlink = 0;
	return v;
}

void p_vertex(vertex *v) {
	show_debug_msg(__func__, "==== Vertex %d: %.2f ====\n", v->id, v->weight);
	//p_ctg_seq(__func__, v->ctg);
	//p_readarray(v->reads, 1);
}

void p_edge(edge *e) {
	show_debug_msg(__func__, "---- Edge %d: %.2f ----\n", e->id, e->weight);
	show_debug_msg(__func__, "Left len: %d; right len: %d\n", e->left_len,
			e->right_len);
	p_vertex(e->left);
	p_vertex(e->right);
	//p_ctg_seq(__func__, e->junc_seq);
	//p_readarray(e->reads, 1);
}

void p_comp_dot(GPtrArray *vertexes, GPtrArray *edges, FILE *dot) {
	edge *e = NULL;
	vertex *v = NULL, *left = NULL, *right = NULL;
	int i = 0, j = 0;
	char entry[BUFSIZ];

	fputs("digraph G {\n", dot);
	fputs("graph [rankdir=LR];", dot);
	for (i = 0; i < vertexes->len; i++) {
		v = (vertex*) g_ptr_array_index(vertexes, i);
		//		if (v->id == 1 || v->id == 2 || v->id == 3 || v->id == 5 || v->id == 7
		//				|| v->id == 8 || v->id == 10)
		//			sprintf(entry, "%d [label=\"%d: %d\" shape=box color=blue]; \n",
		//					v->id, v->id, v->len);
		//		else
		sprintf(entry, "%d [label=\"%d: %d\" shape=box]; \n", v->id, v->id,
				v->len);
		fputs(entry, dot);
	}
	for (j = 0; j < edges->len; j++) {
		e = (edge*) g_ptr_array_index(edges, j);
		left = e->left;
		right = e->right;
		sprintf(entry, "%d -> %d [label=\"%.0f\"]; \n", left->id, right->id,
				e->weight);
		fputs(entry, dot);
	}
	fputs("}\n", dot);
}

void p_comp(comp *c) {
	char fn[128];
	sprintf(fn, "comp.%d.dot", c->id);
	FILE *dot = xopen(fn, "w");
	p_comp_dot(c->vertexes, c->edges, dot);
	fclose(dot);
}

void p_comps(splice_graph *g) {
	uint32_t i = 0;
	comp *c = NULL;
	for (i = 0; i < g->components->len; i++) {
		c = (comp*) g_ptr_array_index(g->components, i);
		p_comp(c);
	}
}

void p_graph(splice_graph *g) {
	FILE *dot = xopen("graph.dot", "w");
	p_comp_dot(g->vertexes, g->edges, dot);
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
	e->id = ++edge_id;
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
		if (g->stack) {
			g_ptr_array_free(g->stack, TRUE);
		}
		if (g->scc) {
			g_ptr_array_free(g->scc, TRUE);
		}
		if (g->components) {
			g_ptr_array_free(g->components, TRUE);
		}
		free(g);
	}
}

/**
 * Remove an edge from the graph
 */
void remove_edge_by_index(splice_graph *g, uint32_t index) {
	uint32_t i = 0;
	edge *e = NULL, *e2 = NULL;
	vertex *v = NULL;
	e = (edge*) g_ptr_array_index(g->edges, index);
	v = e->left;
	for (i = 0; i < v->outs->len; i++) {
		e2 = (edge*) g_ptr_array_index(v->outs, i);
		if (e == e2) {
			g_ptr_array_remove_index_fast(v->outs, i);
			break;
		}
	}
	v = e->right;
	for (i = 0; i < v->ins->len; i++) {
		e2 = (edge*) g_ptr_array_index(v->ins, i);
		if (e == e2) {
			g_ptr_array_remove_index_fast(v->ins, i);
			break;
		}
	}
	g_ptr_array_remove_index_fast(g->edges, index);
}

/**
 * Remove an edge from the graph
 */
void remove_edge_from_g(splice_graph *g, edge *rmv_e) {
	uint32_t i = 0;
	edge *e = NULL;
	for (i = 0; i < g->edges->len; i++) {
		e = (edge*) g_ptr_array_index(g->edges, i);
		if (e == rmv_e) {
			remove_edge_by_index(g, i);
			break;
		}
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

/**
 * Build the graph from the templates and junctions
 */
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

int tarjan_connect(vertex *v, splice_graph *g) {
	edge *e = NULL;
	vertex *w = NULL, *tmp = NULL;
	uint32_t i = 0, j = 0;
	int removed = 0, to_break = 0;

	// Strongly connected component for this vertex
	GPtrArray *scc = g_ptr_array_sized_new(4);
	v->index = g->index;
	v->lowlink = g->index;
	g->index++;
	g_ptr_array_add(g->stack, v);
	for (i = 0; i < v->outs->len; i++) {
		e = (edge*) g_ptr_array_index(v->outs, i);
		w = e->right;
		if (w->index == 0) {
			tarjan_connect(w, g);
			v->lowlink = min(v->lowlink, w->lowlink);
		} else {
			/**
			 for (j = 0; j < g->stack->len; j++) {
			 tmp = (vertex*) g_ptr_array_index(g->stack, j);
			 p_vertex(tmp);
			 }
			 **/
			if (find_in_array(g->stack, (gpointer) w) >= 0) {
				//p_vertex(v);
				//p_vertex(w);
				v->lowlink = min(v->lowlink, w->index);
			}
		}
	}

	if (v->lowlink == v->index) {
		while (g->stack->len) {
			w = (vertex*) g_ptr_array_index(g->stack, g->stack->len - 1);
			g_ptr_array_remove_index_fast(g->stack, g->stack->len - 1);
			g_ptr_array_add(scc, w);
			if (w == v)
				break;
		}
	}

	// Remove edges to break the strongly connected component
	// The edges to remove:
	//		1. identify a vertex A with incoming edges from outside (not in scc)
	//		2. remove all incoming edges to this vertex A.
	// Since a property of SCC is that all vertexes are reachable by others
	//		In this way, other vertexes cannot reach A.
	if (scc->len >= 2) {
		g_ptr_array_add(g->scc, scc);
		for (j = 0; j < scc->len; j++) {
			to_break = 0;
			v = (vertex*) g_ptr_array_index(scc, j);
			// Identify vertex A
			for (i = 0; i < v->ins->len; i++) {
				e = (edge*) g_ptr_array_index(v->ins, i);
				if (find_in_array(scc, e->left) == -1) {
					to_break = 0;
					break;
				} else
					to_break = 1;
			}
			// Remove all incoming SCC edges to A
			if (to_break) {
				for (i = 0; i < v->ins->len; i++) {
					e = (edge*) g_ptr_array_index(v->ins, i);
					if (find_in_array(scc, (gpointer) e->left) >= 0) {
						p_edge(e);
						remove_edge_from_g(g, e);
						removed = 1;
					}
				}
			}
			if (removed)
				break;
		}
	}
	return removed;
}

int tarjan(splice_graph *g) {
	uint32_t i = 0;
	int removed = 0;
	vertex *v = NULL;
	for (i = 0; i < g->vertexes->len; i++) {
		v = (vertex*) g_ptr_array_index(g->vertexes, i);
		if (v->index == 0)
			removed = tarjan_connect(v, g);
	}
	return removed;
}

/**
 * Check whether a vertex is in a SCC.
 */
int vertex_in_scc(splice_graph *g, vertex *v) {
	GPtrArray *scc = NULL;
	uint32_t i = 0;
	for (i = 0; i < g->scc->len; i++) {
		scc = (GPtrArray*) g_ptr_array_index(g->scc, i);
		if (find_in_array(scc, v) >= 0) {
			return 1;
		}
	}
	return 0;
}

void tmp_filter(splice_graph *g) {
	uint32_t i = 0, j = 0;
	vertex *v = NULL;
	edge *e = NULL, *e2 = NULL;
	for (i = 0; i < g->vertexes->len; i++) {
		v = (vertex*) g_ptr_array_index(g->vertexes, i);
		if (v->id == 10) {
			while (v->outs->len)
				g_ptr_array_remove_index_fast(v->outs, 0);
		}
		if (v->id >= 11 && v->id <= 14)
			g_ptr_array_remove_index_fast(g->vertexes, i--);
	}
	for (i = 0; i < g->edges->len; i++) {
		e = (edge*) g_ptr_array_index(g->edges, i);
		if (e->left->id >= 11 && e->left->id <= 14) {
			g_ptr_array_remove_index_fast(g->edges, i--);
			continue;
		}
		if (e->right->id >= 11 && e->right->id <= 14)
			g_ptr_array_remove_index_fast(g->edges, i--);
	}
}

/**
 * Fill in the vertexes and edges into components recursively
 */
void iterate_comp(splice_graph *g, comp *c, vertex *root) {
	edge *e = NULL;
	uint32_t i = 0, j = 0;
	if (root->status > 0)
		return;
	g_ptr_array_add(c->vertexes, root);
	root->status = c->id;
	// Check all incoming edges and vertexes
	for (i = 0; i < root->ins->len; i++) {
		e = (edge*) g_ptr_array_index(root->ins, i);
		if (e->status == 0) {
			e->status = c->id;
			g_ptr_array_add(c->edges, e);
			if (root == e->left)
				iterate_comp(g, c, e->right);
			if (root == e->right)
				iterate_comp(g, c, e->left);
		}
	}
	// Check all outgoing edges and vertexes
	for (i = 0; i < root->outs->len; i++) {
		e = (edge*) g_ptr_array_index(root->outs, i);
		if (e->status == 0) {
			e->status = c->id;
			g_ptr_array_add(c->edges, e);
			if (root == e->left)
				iterate_comp(g, c, e->right);
			if (root == e->right)
				iterate_comp(g, c, e->left);
		}
	}
}

/**
 * Break the graph into components
 */
void break_to_comps(splice_graph *g) {
	uint32_t i = 0, j = 0;
	vertex *v = NULL;
	edge *e = NULL;
	comp *c = NULL;
	for (i = 0; i < g->vertexes->len; i++) {
		v = (vertex*) g_ptr_array_index(g->vertexes, i);
		if (v->status == 0) {
			c = new_comp();
			iterate_comp(g, c, v);
			g_ptr_array_add(g->components, c);
		}
	}
}

/**
 * Reset the vertex and edge status to be zero
 */
void reset_status(splice_graph *g) {
	uint32_t i = 0;
	vertex *v = NULL;
	edge *e = NULL;
	for (i = 0; i < g->vertexes->len; i++) {
		v = (vertex*) g_ptr_array_index(g->vertexes, i);
		v->status = 0;
	}
	for (i = 0; i < g->edges->len; i++) {
		e = (edge*) g_ptr_array_index(g->edges, i);
		e->status = 0;
	}
}

void calc_comp_stat(splice_graph *g) {
	uint32_t i = 0, j = 0;
	comp *c = NULL;
	FILE *stat = xopen("components.csv", "w");
	char entry[BUFSIZE];
	for (i = 0; i < g->components->len; i++) {
		c = (comp*) g_ptr_array_index(g->components, i);
		sprintf(entry, "%d\t%d\t%d\n", c->id, c->vertexes->len, c->edges->len);
		fputs(entry, stat);
	}
	show_debug_msg(__func__, "Statistics of components are output to components.csv.\n");
}

void process_graph(GPtrArray *all_tpls, GPtrArray *all_juncs, hash_map *hm) {
	splice_graph *g = NULL;
	g = build_graph(all_tpls, all_juncs, hm);
	//tmp_filter(g);
	clean_graph(g);
	save_vertexes(g->vertexes);
	// Temporarily just break the cycles.
	while(tarjan(g)){
	}
	break_to_comps(g);
	p_graph(g);
	p_comps(g);
	calc_comp_stat(g);
	reset_status(g);
	determine_paths(g, hm);
}
