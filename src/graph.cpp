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
	return g;
}

gint cmp_junc_by_locus(gpointer a, gpointer b) {
	junction *c_a = *((junction**) a);
	junction *c_b = *((junction**) b);
	return ((c_a->locus) - c_b->locus);
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
		if (v->id == 1 || v->id == 2 || v->id == 3 || v->id == 5 || v->id == 7
				|| v->id == 8 || v->id == 10)
			sprintf(entry, "%d [label=\"%d: %d\" shape=box color=blue]; \n",
					v->id, v->id, v->len);
		else
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
		if (j->main_tpl != t || (j->main_tpl == t && j->branch_tpl == t))
			continue;
		v = new_vertex(t, pre_start, j->locus - pre_start, hm);
		g_ptr_array_add(g->vertexes, v);
		g_ptr_array_add(this_vs, v);
		pre_start = j->locus;
	}
	v = new_vertex(t, pre_start, t->len - pre_start, hm);
	g_ptr_array_add(g->vertexes, v);
	g_ptr_array_add(this_vs, v);

	// Create the edges between vertexes
	for (i = 0; i < main_juncs->len; i++) {
		j = (junction*) g_ptr_array_index(main_juncs, i);
		if (j->main_tpl != t || (j->main_tpl == t && j->branch_tpl == t))
			continue;
		left = (vertex*) g_ptr_array_index(this_vs, i);
		right = (vertex*) g_ptr_array_index(this_vs, i + 1);
		e = new_edge(left, right);
		e->junc_seq = get_junc_seq(t, j->locus, &e->left_len, t, j->locus,
				&e->right_len, max_len);
		e->reads = reads_on_seq(e->junc_seq, hm, N_MISMATCHES);
		e->weight = (float) e->reads->len;
		g_ptr_array_add(left->outs, e);
		g_ptr_array_add(right->ins, e);
		g_ptr_array_add(g->edges, e);
	}

	t->vertexes = this_vs;
}

/**
 * Add edges to connect different vertexes on templates
 */
void connect_tpls(tpl *t, GPtrArray *main_juncs, splice_graph *g, hash_map *hm) {
	junction *j = NULL;
	int i = 0;
	int max_len = (hm->o->read_len - SHORT_BRANCH_SHIFT) * 2;
	edge *e = NULL;
	vertex *branch_v = NULL, *left = NULL, *right = NULL;
	GPtrArray *branch_vs = NULL, *main_vs = NULL;
	for (i = 0; i < main_juncs->len; i++) {
		j = (junction*) g_ptr_array_index(main_juncs, i);
		if (j->main_tpl != t || (j->main_tpl == t && j->branch_tpl == t))
			continue;
		branch_vs = j->branch_tpl->vertexes;
		main_vs = j->main_tpl->vertexes;
		if (j->ori) { // Branch is left, main is right
			branch_v = (vertex*) g_ptr_array_index(branch_vs, branch_vs->len
					- 1);
			right = (vertex*) g_ptr_array_index(main_vs, i + 1);
			e = new_edge(branch_v, right);
			e->junc_seq = get_junc_seq(j->branch_tpl, j->branch_tpl->len,
					&e->left_len, t, j->locus, &e->right_len, max_len);
			g_ptr_array_add(branch_v->outs, e);
			g_ptr_array_add(right->ins, e);
		} else { // Branch is right, main is left
			branch_v = (vertex*) g_ptr_array_index(branch_vs, 0);
			left = (vertex*) g_ptr_array_index(main_vs, i);
			e = new_edge(left, branch_v);
			e->junc_seq = get_junc_seq(t, j->locus, &e->left_len,
					j->branch_tpl, 0, &e->right_len, max_len);
			g_ptr_array_add(left->outs, e);
			g_ptr_array_add(branch_v->ins, e);
		}
		// The weight is not set because of probable alternative paths.
		//e->reads = reads_on_seq(e->junc_seq, hm);
		//e->weight = (float) e->reads->len;
		g_ptr_array_add(g->edges, e);
	}
}

splice_graph *build_graph(GPtrArray *all_tpls, GPtrArray *all_juncs,
		hash_map *hm) {
	splice_graph *g = new_graph();
	GPtrArray *t_juncs = NULL;
	int i = 0;
	tpl *t = NULL;
	for (i = 0; i < all_tpls->len; i++) {
		t = (tpl*) g_ptr_array_index(all_tpls, i);
		t_juncs = tpl_junctions(t, all_juncs);
		p_ctg_seq(__func__, t->ctg);
		break_tpl(t, t_juncs, g, hm);
		g_ptr_array_free(t_juncs, TRUE);
	}

	for (i = 0; i < all_tpls->len; i++) {
		t = (tpl*) g_ptr_array_index(all_tpls, i);
		t_juncs = tpl_junctions(t, all_juncs);
		connect_tpls(t, t_juncs, g, hm);
		g_ptr_array_free(t_juncs, TRUE);
	}
	return g;
}

void process_graph(GPtrArray *all_tpls, GPtrArray *all_juncs, hash_map *hm) {
	splice_graph *g = NULL;
	g = build_graph(all_tpls, all_juncs, hm);
	p_graph(g);
	save_vertexes(g->vertexes);
	determine_paths(g, hm);
}
