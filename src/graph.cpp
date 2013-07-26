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
#include "tpl.hpp"
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
#include "k_hash.h"

int comp_id = 0;
int edge_id = 0;
int vertex_id = 0;

/**
 * Get the reads on a sequence allowing N_MISMATCHES
 */
GPtrArray *reads_on_seq(bwa_seq_t *seq, hash_table *ht, const int n_mismatch) {
	int j = 0, i = 0, read_len = ht->o->read_len;
	bwa_seq_t *part = NULL, *r = NULL;
	GPtrArray *hits = NULL;
	GPtrArray *reads = g_ptr_array_sized_new(32);

	for (i = 0; i <= seq->len - read_len; i++) {
		part = new_seq(seq, read_len, i);
		hits = g_ptr_array_sized_new(4);
		hits = find_both_fr_full_reads(ht, part, hits, n_mismatch);
		for (j = 0; j < hits->len; j++) {
			r = (bwa_seq_t*) g_ptr_array_index(hits, j);
			g_ptr_array_add(reads, r);
		}
		bwa_free_read_seq(1, part);
		g_ptr_array_free(hits, TRUE);
		hits = NULL;
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
	c->vertexes = g_ptr_array_sized_new(2);
	c->edges = g_ptr_array_sized_new(2);
	c->stack = g_ptr_array_sized_new(0);
	c->scc = g_ptr_array_sized_new(0);
	c->id = ++comp_id;
	c->index = 0;
	return c;
}

void destroy_comp(comp *c) {
	if (c) {
		if (c->vertexes)
			g_ptr_array_free(c->vertexes, TRUE);
		if (c->edges)
			g_ptr_array_free(c->edges, TRUE);
		if (c->stack)
			g_ptr_array_free(c->stack, TRUE);
		if (c->scc)
			g_ptr_array_free(c->scc, TRUE);
		free(c);
	}
}

vertex *new_vertex(tpl *t, int start, int len, hash_table *ht) {
	vertex *v = (vertex*) malloc(sizeof(vertex));
	v->ctg = new_seq(t->ctg, len, start);
	v->len = len;
	v->ins = g_ptr_array_sized_new(0);
	v->outs = g_ptr_array_sized_new(0);
	v->reads = reads_on_seq(v->ctg, ht, N_MISMATCHES);
	v->weight = (float) v->reads->len;
	v->status = 0;
	v->id = ++vertex_id;
	v->index = 0;
	v->lowlink = 0;
    v->from = t;
	return v;
}

void p_vertex(vertex *v) {
	show_debug_msg(__func__, "==== Vertex %d: %.2f ====\n", v->id, v->weight);
	//p_ctg_seq(__func__, v->ctg);
	//p_readarray(v->reads, 1);
}

void p_edge(edge *e) {
	show_debug_msg(__func__, "---- Edge %d: %.2f status: %d ----\n", e->id,
			e->weight, e->status);
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
		sprintf(entry, "%d -> %d [label=\"\"]; \n", left->id, right->id);
		fputs(entry, dot);
	}
	fputs("}\n", dot);
}

void p_comp(comp *c) {
	char fn[128];
	sprintf(fn, "../simu_out/components/comp.%d.dot", c->id);
	FILE *dot = xopen(fn, "w");
	p_comp_dot(c->vertexes, c->edges, dot);
	fclose(dot);
    sprintf(fn, "../simu_out/components/comp.%d.fa", c->id);
    save_vertexes(c->vertexes, fn);
}

void p_comps(splice_graph *g) {
	uint32_t i = 0;
	comp *c = NULL;
	for (i = 0; i < g->components->len; i++) {
		c = (comp*) g_ptr_array_index(g->components, i);
		p_comp(c);
	}
}

void p_graph(splice_graph *g, char *fn) {
	FILE *dot = xopen(fn, "w");
	p_comp_dot(g->vertexes, g->edges, dot);
	fclose(dot);
}

void save_vertexes(GPtrArray *vertexes, char *fn) {
	int i = 0;
	vertex *v = NULL;
	char entry[BUFSIZ];
	FILE *v_fp = xopen(fn, "w");
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

edge *bridge_edge(vertex *left, vertex *right) {
	edge *e = NULL;
	uint32_t i = 0;
	for (i = 0; i < left->outs->len; i++) {
		e = (edge*) g_ptr_array_index(left->outs, i);
		if (e->right == right)
			return NULL;
	}
	return new_edge(left, right);
}

void destroy_edge(edge *eg) {
	if (eg) {
		bwa_free_read_seq(1, eg->junc_seq);
		if (eg->reads)
			g_ptr_array_free(eg->reads, TRUE);
		free(eg);
	}
}

void destroy_vertex(vertex *v) {
	if (v) {
		bwa_free_read_seq(1, v->ctg);
		if (v->ins)
			g_ptr_array_free(v->ins, TRUE);
		if (v->outs)
			g_ptr_array_free(v->outs, TRUE);
		if (v->reads)
			g_ptr_array_free(v->reads, TRUE);
		free(v);
	}
}

void destroy_graph(splice_graph *g) {
	int i = 0;
	vertex *v = NULL;
	edge *e = NULL;
	comp *c = NULL;
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
			for (i = 0; i < g->components->len; i++) {
				c = (comp*) g_ptr_array_index(g->components, i);
				destroy_comp(c);
			}
			g_ptr_array_free(g->components, TRUE);
		}
		free(g);
	}
}

/**
 * Remove an edge from the graph
 */
void remove_edge_by_index(comp *g, uint32_t index) {
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
void remove_edge_from_c(comp *c, edge *rmv_e) {
	uint32_t i = 0;
	edge *e = NULL;
	for (i = 0; i < c->edges->len; i++) {
		e = (edge*) g_ptr_array_index(c->edges, i);
		if (e == rmv_e) {
			remove_edge_by_index(c, i);
			break;
		}
	}
}

/**
 * Break a template into a set of edges and vertexes
 */
void break_tpl(tpl *t, GPtrArray *main_juncs, splice_graph *g, hash_table *ht) {
	junction *j = NULL;
	vertex *v = NULL, *left = NULL, *right = NULL;
	edge *e = NULL;
	int i = 0, pre_start = 0;
	int max_len = (ht->o->read_len - SHORT_BRANCH_SHIFT) * 2;
	GPtrArray *this_vs = NULL;

	if (t->len <= 0)
		return;

	this_vs = g_ptr_array_sized_new(main_juncs->len + 1);

	// Break the main template into vertexes
	// The junctions are sorted by locus already
	for (i = 0; i < main_juncs->len; i++) {
		j = (junction*) g_ptr_array_index(main_juncs, i);
		if (j->locus > t->len || j->locus < 0) {
			continue;
		}
		if (j->main_tpl != t || (j->main_tpl == t && j->branch_tpl == t)
				|| j->locus == pre_start)
			continue;
		v = new_vertex(t, pre_start, j->locus - pre_start, ht);
		g_ptr_array_add(g->vertexes, v);
		g_ptr_array_add(this_vs, v);
		pre_start = j->locus;
		show_debug_msg(__func__, "Template %d, Vertex %d \n", t->id, v->id);
	}

	// Create the last vertex
	v = new_vertex(t, pre_start, t->len - pre_start, ht);
	g_ptr_array_add(g->vertexes, v);
	g_ptr_array_add(this_vs, v);

	show_debug_msg(__func__, "Template %d, Vertex %d \n", t->id, v->id);

	// Create the edges between vertexes
	pre_start = 0;
	for (i = 0; i < this_vs->len - 1; i++) {
		left = (vertex*) g_ptr_array_index(this_vs, i);
		right = (vertex*) g_ptr_array_index(this_vs, i + 1);
		e = bridge_edge(left, right);
		if (!e)
			continue;
		e->junc_seq = get_junc_seq(t, pre_start, &e->left_len, t, pre_start,
				&e->right_len, max_len);
		e->len = e->left_len + e->right_len;
        show_debug_msg(__func__, "Aligning reads to edge %d...\n", e->id);
		e->reads = reads_on_seq(e->junc_seq, ht, N_MISMATCHES);
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
	if (!t || !t->vertexes || t->vertexes->len == 0)
		return NULL;
	if (t->vertexes->len == 1) {
		this_v = (vertex*) g_ptr_array_index(t->vertexes, 0);
		return this_v;
	}
	if (locus == 0 && ori == 0)
		return (vertex*) g_ptr_array_index(t->vertexes, 0);
	if (locus == t->len && ori == 1)
		return (vertex*) g_ptr_array_index(t->vertexes, t->vertexes->len - 1);
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
					if (!left || !right)
						continue;
					e = bridge_edge(left, right);
					if (!e)
						continue;
					e->junc_seq = get_junc_seq(t, left_locus, &e->left_len, t,
							right_locus, &e->right_len, max_len);
					g_ptr_array_add(left->outs, e);
					g_ptr_array_add(right->ins, e);
					g_ptr_array_add(g->edges, e);
					/**
					 show_debug_msg(__func__,
					 "New edge between vertexes %d and %d \n", left->id,
					 right->id);
					 **/
				}
			}
		}
	}
	g_ptr_array_free(zero_juncs, TRUE);
}

/**
 * Add edges to connect different vertexes on templates
 */
void connect_tpls(tpl *t, GPtrArray *main_juncs, splice_graph *g, hash_table *ht) {
	junction *j = NULL;
	int i = 0, n_zeros = 0, nth_zero = 0;
	int max_len = (ht->o->read_len - SHORT_BRANCH_SHIFT) * 2;
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
		if (j->locus < 0 || j->locus >= t->len)
			continue;

		/**
		 show_debug_msg(__func__, "Junction %d: nth_zero: %d \n", i, nth_zero);
		 p_junction(j);
		 **/

		branch_vs = j->branch_tpl->vertexes;
		main_vs = j->main_tpl->vertexes;
		if (j->ori) { // Branch is left, main is right
			branch_v = (vertex*) g_ptr_array_index(branch_vs,
					branch_vs->len - 1);
			// Shifting on the main_tpl to find connecting vertex
			right = find_vertex_by_locus(t, j->locus, 0);
			if (!right)
				continue;
			e = bridge_edge(branch_v, right);
			if (!e)
				continue;
			e->junc_seq = get_junc_seq(j->branch_tpl, j->branch_tpl->len,
					&e->left_len, t, j->locus, &e->right_len, max_len);
			g_ptr_array_add(branch_v->outs, e);
			g_ptr_array_add(right->ins, e);
		} else { // Branch is right, main is left
			branch_v = (vertex*) g_ptr_array_index(branch_vs, 0);
			// Shifting on the main_tpl to find connecting vertex
			left = find_vertex_by_locus(t, j->locus, 1);
			if (!left)
				continue;
			e = bridge_edge(left, branch_v);
			if (!e)
				continue;
			e->junc_seq = get_junc_seq(t, j->locus, &e->left_len,
					j->branch_tpl, 0, &e->right_len, max_len);
			g_ptr_array_add(left->outs, e);
			g_ptr_array_add(branch_v->ins, e);
		}
		e->len = e->left_len + e->right_len;
		// The weight is not set because of probable alternative paths.
		//e->reads = reads_on_seq(e->junc_seq, ht);
		//e->weight = (float) e->reads->len;
		g_ptr_array_add(g->edges, e);
	}
}

/**
 * Build the graph from the templates and junctions
 */
splice_graph *build_graph(GPtrArray *all_tpls, GPtrArray *all_juncs,
		hash_table *ht) {
	splice_graph *g = new_graph();
	GPtrArray *t_juncs = NULL;
	vertex *v = NULL;
	int i = 0;
	tpl *t = NULL;
	int start_index = 0;
	show_debug_msg(__func__, "Breaking templates into vertexes...\n");
	for (i = 0; i < all_tpls->len; i++) {
		t = (tpl*) g_ptr_array_index(all_tpls, i);
		t_juncs = tpl_junctions(t, all_juncs, start_index, 1);
		g_ptr_array_sort(t_juncs, (GCompareFunc) cmp_junc_by_locus);
		start_index += t_juncs->len;
		break_tpl(t, t_juncs, g, ht);
		//p_tpl_juncs(t, t_juncs);
		g_ptr_array_free(t_juncs, TRUE);
	}

	start_index = 0;
	show_debug_msg(__func__, "Building subgraph connected to templates...\n");
	for (i = 0; i < all_tpls->len; i++) {
		t = (tpl*) g_ptr_array_index(all_tpls, i);
		t_juncs = tpl_junctions(t, all_juncs, start_index, 1);
		g_ptr_array_sort(t_juncs, (GCompareFunc) cmp_junc_by_locus);
		start_index += t_juncs->len;
		remove_zero_vertexes(t, t_juncs, ht->o->read_len, g);
		connect_tpls(t, t_juncs, g, ht);
		//p_tpl_juncs(t, t_juncs);
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
 * Remove isolated vertexes with length smaller than 10;
 * Remove edges with status -1
 */
void clean_graph(splice_graph *g) {
	vertex *v = NULL;
	uint32_t i = 0, j = 0;
	edge *e = NULL;
	for (i = 0; i < g->vertexes->len; i++) {
		v = (vertex*) g_ptr_array_index(g->vertexes, i);
		// Consider removing vertexes with length <10
		if (v->len >= 10)
			continue;
		// If it's isolated, remove
		if (v->ins->len == 0 && v->outs->len == 0)
			g_ptr_array_remove_index_fast(g->vertexes, i--);
	}
	for (i = 0; i < g->edges->len; i++) {
		e = (edge*) g_ptr_array_index(g->edges, i);
		if (e->status == -1) {
			g_ptr_array_remove_index_fast(g->edges, i--);
			destroy_edge(e);
		}
	}
}

int tarjan_connect(vertex *v, comp *g) {
	edge *e = NULL;
	vertex *w = NULL, *tmp = NULL;
	uint32_t i = 0, j = 0;
	// A flag to indicate whether some edges are removed at this round
	int removed = 0;
	int to_break = 0;

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
			removed |= tarjan_connect(w, g);
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
	//		2. remove all incoming edges from scc vertexes to this vertex A.
	// Since a property of SCC is that all vertexes are reachable by others
	//		In this way, other vertexes cannot reach A.
	if (scc->len >= 2) {
		g_ptr_array_add(g->scc, scc);
		for (j = 0; j < scc->len; j++) {
			to_break = 0;
			v = (vertex*) g_ptr_array_index(scc, j);
			//p_vertex(v);
			// Identify vertex A
			for (i = 0; i < v->ins->len; i++) {
				e = (edge*) g_ptr_array_index(v->ins, i);
				if (find_in_array(scc, (gpointer) e->left) == -1) {
					to_break = 1;
					break;
				} else
					to_break = 0;
			}
			// Remove all incoming SCC edges from SCC vertexes to A
			if (to_break) {
				//show_debug_msg(__func__, "Removing edges to break the SCC...\n");
				for (i = 0; i < v->ins->len; i++) {
					e = (edge*) g_ptr_array_index(v->ins, i);
					//p_edge(e);
					if (e->right == v && find_in_array(scc, (gpointer) e->left) >= 0) {
						// Will be removed from the global graph later by clean_graph
						e->status = -1;
						//show_debug_msg(__func__, "To remove edge %d...\n",
						//		e->id);
						//p_edge(e);
						remove_edge_from_c(g, e);
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

/**
 * Tarjan's algorithm to detect strongly connected components (SCC), and break them
 */
int tarjan(comp *c) {
	uint32_t i = 0;
	int removed = 0;
	vertex *v = NULL;
	// Reset the index and stack in case they are not empty
	while (c->stack->len > 0)
		g_ptr_array_remove_index_fast(c->stack, 0);
	c->index = 0;
	for (i = 0; i < c->vertexes->len; i++) {
		v = (vertex*) g_ptr_array_index(c->vertexes, i);
		v->index = 0;
		v->lowlink = 0;
	}
	for (i = 0; i < c->vertexes->len; i++) {
		v = (vertex*) g_ptr_array_index(c->vertexes, i);
		if (v->index == 0)
			removed |= tarjan_connect(v, c);
	}
	return removed;
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
	comp *c = NULL;
	for (i = 0; i < g->vertexes->len; i++) {
		v = (vertex*) g_ptr_array_index(g->vertexes, i);
		if (v->status == 0) {
			c = new_comp();
			iterate_comp(g, c, v);
			// Try to break cycles if there are some edges removed from last round
			while (tarjan(c) > 0) {
			}
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

/**
 * Output the counts of edges and vertexes in components
 */
void calc_comp_stat(splice_graph *g) {
	uint32_t len = 0, i = 0, j = 0;
	comp *c = NULL;
	FILE *stat = xopen("components.csv", "w");
	vertex *v = NULL;
	char entry[BUFSIZE];
	for (i = 0; i < g->components->len; i++) {
		c = (comp*) g_ptr_array_index(g->components, i);
		len = 0;
		for (j = 0; j < c->vertexes->len; j++) {
			v = (vertex*) g_ptr_array_index(c->vertexes, j);
			len += v->len;
		}
		sprintf(entry, "%d\t%d\t%d\t%d\n", c->id, c->vertexes->len,
				c->edges->len, len);
		fputs(entry, stat);
	}
	fclose(stat);
	show_debug_msg(__func__,
			"Statistics of components are output to components.csv.\n");
}

void process_graph(GPtrArray *all_tpls, GPtrArray *all_juncs, hash_table *ht) {
	splice_graph *g = NULL;
	show_msg(__func__, "Building the splice graph...\n");
	g = build_graph(all_tpls, all_juncs, ht);
	show_msg(__func__, "Simplifying the splice graph...\n");
	p_graph(g, "graph.ori.dot");
	clean_graph(g);
	save_vertexes(g->vertexes, "../simu_out/vertexes.fa");

	show_msg(__func__, "Breaking into components...\n");
	break_to_comps(g);
	// Some edges may be marked as 'dead' when breaking SCCs
	clean_graph(g);
	p_graph(g, "graph.dot");
	//p_comps(g);
	calc_comp_stat(g);
	reset_status(g);
	show_msg(__func__, "Running EM to get paths...\n");
	determine_paths(g, ht);
	destroy_graph(g);
}
