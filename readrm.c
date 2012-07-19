/*
 * readrm.c
 *
 *  Created on: Nov 30, 2011
 *      Author: carl
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <glib.h>
#include "roadmap.h"
#include "bwase.h"
#include "peseq.h"
#include "pool.h"
#include "pealn.h"
#include "edge.h"
#include "edgelist.h"
#include "readrm.h"
#include "utils.h"

void p_flat_eg(const edge *eg) {
	edge *eg_i;
	int i = 0;
	show_debug_msg(__func__, "[%p] [%d: %d]\n", eg, eg->id, eg->len);
	show_debug_msg(__func__, "\t Out: %p \n", eg->out_egs);
	for (i = 0; i < eg->out_egs->len; i++) {
		eg_i = g_ptr_array_index(eg->out_egs, i);
		show_debug_msg(__func__, "\t Out: [%p] [%d: %d]\n", eg_i, eg_i->id,
				eg_i->len);
	}
	show_debug_msg(__func__, "\t In: %p \n", eg->in_egs);
	for (i = 0; i < eg->in_egs->len; i++) {
		eg_i = g_ptr_array_index(eg->in_egs, i);
		show_debug_msg(__func__, "\t In: [%p] [%d: %d]\n", eg_i, eg_i->id,
				eg_i->len);
	}
	if (eg->right_ctg)
		show_debug_msg(__func__, "\t Right: [%p] [%d: %d]\n", eg->right_ctg,
				eg->right_ctg->id, eg->right_ctg->len);
}

void w_flat_eg(const edge *eg, FILE *debug) {
	edge *eg_i;
	int i = 0;
	char content[BUFSIZE];
	sprintf(content, "[%p] [%d: %d]\n", eg, eg->id, eg->len);
	fputs(content, debug);
	sprintf(content, "\t Out: %p \n", eg->out_egs);
	fputs(content, debug);
	for (i = 0; i < eg->out_egs->len; i++) {
		eg_i = g_ptr_array_index(eg->out_egs, i);
		sprintf(content, "\t Out: [%p] [%d: %d]\n", eg_i, eg_i->id, eg_i->len);
		fputs(content, debug);
	}
	sprintf(content, "\t In: %p \n", eg->in_egs);
	fputs(content, debug);
	for (i = 0; i < eg->in_egs->len; i++) {
		eg_i = g_ptr_array_index(eg->in_egs, i);
		sprintf(content, "\t In: [%p] [%d: %d]\n", eg_i, eg_i->id, eg_i->len);
		fputs(content, debug);
	}
	if (eg->right_ctg) {
		sprintf(content, "\t Right: [%p] [%d: %d]\n", eg->right_ctg,
				eg->right_ctg->id, eg->right_ctg->len);
		fputs(content, debug);
	}
}

void p_eg(const edge *root, int level) {
	int i = 0, next_level = level + 1;
	edgearray *in_egs, *out_egs;
	edge *e;

	if (!root)
		return;
	in_egs = root->in_egs;
	out_egs = root->out_egs;
	printf("[p_edge] %d [%p]: ", root->id, root);
	while (level--) {
		printf("\t");
	}
	for (i = 0; i < in_egs->len; i++) {
		e = g_ptr_array_index(in_egs, i);
		if (e)
			printf("%d,", e->id);
	}
	if (root->right_ctg)
		printf("->%d: (%d, %d) \n", root->len, root->right_ctg->id,
				root->r_shift);
	else {
		printf("->%d \n", root->len);
		for (i = 0; i < out_egs->len; i++) {
			p_eg(g_ptr_array_index(out_egs, i), next_level);
		}
	}
}

void p_rm(const roadmap *rm) {
	int i = 0;
	edge *root;
	if (!rm || !rm->start_eg_n) {
		printf("[p_rm] Empty Roadmap! \n");
	}
	for (; i < rm->start_eg_n; i++) {
		root = g_ptr_array_index(rm->start_egs, i);
		if (root) {
			printf("[p_rm] -------------------------------------------- \n");
			p_eg(root, 0);
		}
	}
}

void graph_by_edges(const edgearray *all_edges, char *dotfile) {
	int i = 0, j = 0;
	char *edge_str = malloc(BUFSIZE);
	FILE *dot;
	edge *eg, *eg_i;
	edgearray *out_egs, *in_egs;
	if (!all_edges || all_edges->len == 0) {
		printf("[graph_by_edges] Empty Roadmap! \n");
		return;
	}
	dot = xopen(dotfile, "w");
	setbuf(dot, NULL);
	fputs("digraph g { \n\trankdir = LR \n", dot);
	for (i = 0; i < all_edges->len; i++) {
		eg = g_ptr_array_index(all_edges, i);

		if (eg && eg->alive) {
			in_egs = eg->in_egs;
			out_egs = eg->out_egs;
			sprintf(edge_str, "\t%d [shape=box, label=\"%d:%d\"] \n", eg->id,
					eg->id, eg->len);
			fputs(edge_str, dot);
			if (eg->right_ctg) {
				sprintf(edge_str, "\t%d -> %d [label=\"s=%d\"] \n", eg->id,
						eg->right_ctg->id, eg->r_shift);
				fputs(edge_str, dot);
			} else {
				for (j = 0; j < out_egs->len; j++) {
					eg_i = g_ptr_array_index(out_egs, j);
					sprintf(edge_str, "\t%d -> %d\n", eg->id, eg_i->id);
					fputs(edge_str, dot);
				}

			}
		}
	}
	fputs("} \n", dot);
	fclose(dot);
	free(edge_str);
}

void save_node(const edge *eg, FILE *ass_fa) {
	char *h;
	int i = 0;
	bwa_seq_t *contig;
	edgearray *out_egs;
	if (!eg)
		return;
	h = malloc(BUFSIZE);
	contig = eg->contig;
	out_egs = eg->out_egs;
	if (contig) {
		//		fprintf(stderr, "[save_node] Writing contig %d to disk... %d\n",
		//				eg->id, contig->len);
		sprintf(h, ">%d len=%d \n", eg->id, contig->len);
		save_con(h, contig, ass_fa);
	}
	free(h);
	for (i = 0; i < out_egs->len; i++) {
		save_node(g_ptr_array_index(out_egs, i), ass_fa);
	}
}

void save_edges(edgearray *pfd_ctg_ids, FILE *ass_fa, const int ori,
		const int p_all, const int min_len) {
	int i = 0;
	edge *eg;
	char *h;
	bwa_seq_t *contig;
	if (!pfd_ctg_ids || pfd_ctg_ids->len <= 0)
		return;
	h = malloc(BUFSIZE);
	for (i = 0; i < pfd_ctg_ids->len; i++) {
		eg = g_ptr_array_index(pfd_ctg_ids, i);
		if (p_all || (eg && eg->alive && eg->contig && eg->len > min_len)) {
			contig = eg->contig;
			if (ori)
				seq_reverse(contig->len, contig->seq, 0);
			sprintf(h, ">%d len=%d \n", eg->id, contig->len);
			save_con(h, contig, ass_fa);
		}
	}
	free(h);
}

void save_tx(const roadmap *rm, FILE *ass_fa) {
	char *h = malloc(BUFSIZE);
	unsigned int i = 0, j = 0;
	edge *root = 0;
	edgearray *out_egs;
	bwa_seq_t *contig;
	if (!rm || !rm->start_eg_n)
		return;
	for (j = 0; j < rm->start_eg_n; j++) {
		root = g_ptr_array_index(rm->start_egs, j);
		contig = root->contig;
		sprintf(h, ">%d len=%d \n", root->id, contig->len);
		// fprintf(stderr, "[save_tx] Writing root contig %d to disk... %d\n",
		//		root->id, contig->len);
		save_con(h, contig, ass_fa);
		out_egs = root->out_egs;
		for (i = 0; i < out_egs->len; i++) {
			save_node(g_ptr_array_index(out_egs, i), ass_fa);
		}
	}
	free(h);
}

int is_start_eg(const roadmap *rm, const edge *eg) {
	edge *tmp;
	unsigned int i;
	if (!rm || !eg || eg->id < 0)
		return 0;
	for (i = rm->start_eg_n; i > 0; i--) {
		tmp = g_ptr_array_index(rm->start_egs, i);
		if (tmp && tmp->id == eg->id)
			return 1;
	}
	return 0;
}

int has_edge(const edge *eg, const int contig_id, const int in_out) {
	edge *tmp;
	unsigned int i, n_edges;
	edgearray *edges;
	if (!eg || contig_id < 0)
		return 0;
	edges = in_out ? eg->in_egs : eg->out_egs;
	n_edges = edges->len;
	if (edges) {
		for (i = 0; i < n_edges; i++) {
			tmp = g_ptr_array_index(edges, i);
			if (tmp->id == contig_id)
				return 1;
		}
	}
	return 0;
}

/**
 * Check whether two edges are siblings.
 */
int is_sbl(const edge *eg_1, const edge *eg_2) {
	edgearray *in_egs, *out_egs;
	edge *eg_i;
	int i = 0, j = 0;
	if (!eg_1 || !eg_2)
		return 0;
	in_egs = eg_1->in_egs;
	for (i = 0; i < in_egs->len; i++) {
		eg_i = g_ptr_array_index(in_egs, i);
		out_egs = eg_i->out_egs;
		for (j = 0; j < out_egs->len; j++) {
			eg_i = g_ptr_array_index(out_egs, j);
			if (eg_i == eg_2)
				return 1;
		}
	}
	in_egs = eg_2->in_egs;
	for (i = 0; i < in_egs->len; i++) {
		eg_i = g_ptr_array_index(in_egs, i);
		out_egs = eg_i->out_egs;
		for (j = 0; j < out_egs->len; j++) {
			eg_i = g_ptr_array_index(out_egs, j);
			if (eg_i == eg_1)
				return 1;
		}
	}
	return 0;
}

edgearray *get_sbls(const edge *eg) {
	edgearray *sbls, *in_egs, *out_egs;
	edge *eg_i, *eg_j;
	int i = 0, j = 0;
	if (!eg || !eg->alive)
		return 0;
	sbls = g_ptr_array_new();
	if (eg->ori) {
		out_egs = eg->out_egs;
		for (i = 0; i < out_egs->len; i++) {
			eg_i = g_ptr_array_index(out_egs, i);
			in_egs = eg_i->in_egs;
			for (j = 0; j < in_egs->len; j++) {
				eg_j = g_ptr_array_index(in_egs, j);
				if (eg_j && eg_j->alive && eg_j->id != eg->id) {
					g_ptr_array_add(sbls, eg_j);
				}
			}
		}
	} else {
		in_egs = eg->in_egs;
		for (i = 0; i < in_egs->len; i++) {
			eg_i = g_ptr_array_index(in_egs, i);
			out_egs = eg_i->out_egs;
			for (j = 0; j < out_egs->len; j++) {
				eg_j = g_ptr_array_index(out_egs, j);
				if (eg_j && eg_j->alive && eg_j->id != eg->id) {
					g_ptr_array_add(sbls, eg_j);
				}
			}
		}
	}
	return sbls;
}
