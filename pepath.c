/*
 * pepath.c
 *
 *  Created on: Sep 6, 2011
 *      Author: Cai Shaojiang
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "pepath.h"
#include "roadmap.h"
#include "bwase.h"
#include "pechar.h"
#include "utils.h"
#include "pehash.h"

int path_id = 0;

extern unsigned char nst_nt4_table[256];

rm_path *new_path() {
	rm_path *p = (rm_path*) malloc(sizeof(rm_path));
	p->edges = g_ptr_array_sized_new(16);
	p->n_ctgs = 0;
	p->len = 0;
	return p;
}

void p_path(const rm_path *p) {
	int i = 0;
	edge *contig = NULL;
	bwa_seq_t *query;
	if (!p) {
		printf("Empty path! \n");
		return;
	}
	printf("[p_path] ----------------------------------------\n");
	printf("[p_path] Path %d (%p): %d \n", p->id, p, p->len);
	for (i = 0; i < p->edges->len; i++) {
		contig = g_ptr_array_index(p->edges, i);
		printf("[p_path] \t Contig %d [%s]: %d (%d, %d)\n", contig->id,
				query->name, contig->len,
				contig->right_ctg ? contig->right_ctg->id : 0, contig->r_shift);
	}
	printf("[p_path] ----------------------------------------\n");
}

/*
 * Recursive method to get all edges reachable by the input edge
 */
void get_block_edges(edge *eg, GPtrArray *block) {
	edgelist *in_out = NULL;
	edge *in_out_eg = NULL;
	int i = 0;
	eg->visited = 1;
	for (i = 0; i < eg->in_egs->len; i++) {
		in_out_eg = g_ptr_array_index(eg->in_egs, i);
		if (in_out_eg->alive && !in_out_eg->visited) {
			g_ptr_array_add(block, in_out_eg);
			in_out_eg->visited = 1;
			get_block_edges(in_out_eg, block);
		}
	}
	for (i = 0; i < eg->out_egs->len; i++) {
		in_out_eg = g_ptr_array_index(eg->out_egs, i);
		if (in_out_eg->alive && !in_out_eg->visited) {
			g_ptr_array_add(block, in_out_eg);
			in_out_eg->visited = 1;
			get_block_edges(in_out_eg, block);
		}
	}
	if (eg->right_ctg && eg->right_ctg->alive && !eg->right_ctg->visited) {
		g_ptr_array_add(block, eg->right_ctg);
		eg->right_ctg->visited = 1;
	}
}

void iterate_block(GPtrArray *block, GPtrArray *paths) {
	edge *eg = NULL;
	rm_path *path = new_path();
	int i = 0, has_fresh = 1;

	if (block->len == 1) {
		g_ptr_array_add(paths, g_ptr_array_index(block, 0));
		return;
	}
	if (block->len >= 32) {
		show_debug_msg(__func__, "Too many edges in this block (>=32) \n");
		for (i = 0; i < block->len; i++) {
			eg = g_ptr_array_index(block, i);
			p_flat_eg(eg);
		}
		return;
	}
	while (has_fresh) {
		has_fresh = 0;
		for (i = 0; i < block->len; i++) {
			eg = g_ptr_array_index(block, i);
			if (eg->level == -1) {

			}
		}
	}
}

GPtrArray *report_path(edgearray *all_edges) {
	int block_size = 0, i = 0, j = 0;
	GPtrArray *block = NULL, *paths = NULL;
	edge *eg = NULL;
	edgelist *in_out = NULL;

	paths = g_ptr_array_sized_new(1024);
	for (i = 0; i < all_edges->len; i++) {
		eg = g_ptr_array_index(all_edges, i);
		if (eg->alive && !eg->visited) {
			block = g_ptr_array_sized_new(128);
			get_block_edges(eg, block);
			show_debug_msg(__func__, "NEW BLOCK ---------------------------\n");
			p_flat_eg(eg);
			for (j = 0; j < block->len; j++) {
				eg = g_ptr_array_index(block, j);
				p_flat_eg(eg);
			}
			iterate_block(block, paths);
			g_ptr_array_free(block, TRUE);
		}
	}
	return paths;
}

edgearray *load_rm(const char *rm_dump_file, const hash_table *ht) {
	edgearray *edges = NULL;
	FILE *dump_fp = NULL;
	int i = 0, j = 0, len = 0, tmp = 0, eg_id = 0;
	edge *eg = NULL, *in_out_eg = NULL;
	bwa_seq_t *r = NULL;

	dump_fp = xopen(rm_dump_file, "r");
	fread(&len, sizeof(int), 1, dump_fp);
	show_debug_msg(__func__, "# of edges: %d \n", len);
	edges = g_ptr_array_sized_new(len);
	for (i = 0; i < len; i++) {
		fread(eg, sizeof(edge), 1, dump_fp);
		fread(&tmp, sizeof(int), 1, dump_fp);
		eg->out_egs = g_ptr_array_sized_new(tmp);
		for (j = 0; j < tmp; j++) {
			fread(&eg_id, sizeof(int), 1, dump_fp);
			g_ptr_array_add(eg->out_egs, eg_id);
		}
		fread(&tmp, sizeof(int), 1, dump_fp);
		eg->in_egs = g_ptr_array_sized_new(tmp);
		for (j = 0; j < tmp; j++) {
			fread(&eg_id, sizeof(int), 1, dump_fp);
			g_ptr_array_add(eg->in_egs, eg_id);
		}
		fread(&tmp, sizeof(int), 1, dump_fp);
		eg->reads = g_ptr_array_sized_new(tmp);
		for (j = 0; j < tmp; j++) {
			fread(&eg_id, sizeof(int), 1, dump_fp);
			g_ptr_array_add(eg->in_egs, eg_id);
		}
		g_ptr_array_add(edges, eg);
		show_debug_msg(__func__, "Edge [%d, %d] \n", eg->id, eg->len);
//		p_flat_eg(eg);
//		fread(eg->id, sizeof(int), 1, dump_fp);
//		fread(tmp, sizeof(int), 1, dump_fp);
//		for (j = 0; j < tmp; j++) {
//			eg->alive = 0;
//		}
	}

	fclose(dump_fp);
}

int pe_path(int argc, char *argv[]) {
	int c;
	clock_t t = clock();
	hash_table *ht = NULL;

	fprintf(stderr, "%s \n", argv[1]);
	fprintf(stderr, "%s \n", argv[2]);
	ht = pe_load_hash(argv[1]);

	load_rm(argv[2], ht);

	fprintf(stderr, "[pe_path] Done: %.2f sec\n", (float) (clock() - t)
			/ CLOCKS_PER_SEC);
	return 0;
}
