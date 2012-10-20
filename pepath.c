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
#include "rnaseq.h"
#include "peseq.h"

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

int count_comma(const char *str, const int len) {
	int commas = 0, i = 0;
	for (i = 0; i < len; i++) {
		if (str[i] == ',')
			commas++;
	}
	return commas;
}

edgearray *load_rm(const hash_table *ht, const char *rm_dump_file,
		const char *rm_reads_file, const char *contig_file) {
	edgearray *edges = NULL;
	FILE *dump_fp = NULL, *reads_fp = NULL, *contigs_fp = NULL;
	int i = 0, j = 0, len = 0, n_ctgs = 0, char_space = BUFSIZ - 1, char_len = 0;
	edge *eg = NULL, *in_out_eg = NULL;
	bwa_seq_t *r = NULL, *contigs = NULL;
	readarray *reads = NULL;
	char buf[BUFSIZ];
	char *read_str = (char*) calloc(char_space, sizeof(char)); // allocate buffer.
	char *attr[16], ch = 0;
	char *read_ids = NULL;

	dump_fp = xopen(rm_dump_file, "r");
	reads_fp = xopen(rm_reads_file, "r");
	contigs = load_reads(contig_file, &n_ctgs);

	while ((ch = fgetc(reads_fp)) != EOF) {
		if (char_len >= char_space) { // time to expand ?
			char_space += 1;
			kroundup32(char_space); // expand to double the current size of anything similar.
			read_str = realloc(read_str, char_space); // re allocate memory.
		}
		if (ch == '\n') {
			i = 0;
			attr[0] = strtok(read_str, "\t");
			while (attr[i] != NULL) { //ensure a pointer was found
				attr[++i] = strtok(NULL, "\t"); //continue to tokenize the string
			}
			eg = new_eg();
			eg->id = atoi(attr[0]);

			len = count_comma(attr[1], char_len);
			char *read_ids[len];
			read_ids[0] = strtok(attr[1], ",");
			show_debug_msg(__func__, "read_ids[0]: %s, len: %d\n", read_ids[0], len);
			i = 0;
			while (read_ids[i] != NULL) { //ensure a pointer was found
				read_ids[++i] = strtok(NULL, ","); //continue to tokenize the string
				if (strcmp(read_ids[i], "") == 0) {
					break;
				}
				r = &ht->seqs[atoi(read_ids[i])];
				g_ptr_array_add(eg->reads, r);
				show_debug_msg(__func__, read_ids[i]);
				p_query(__func__, r);
			}
			p_flat_eg(eg);

			free(read_str);
			read_str = (char*) calloc(BUFSIZ, sizeof(char)); // allocate buffer.
			char_space = BUFSIZ;
			char_len = 0;
		} else {
			read_str[char_len] = ch; // stuff in buffer.
			char_len++;
		}
	}

	fclose(reads_fp);
	fclose(dump_fp);
}

int pe_path(int argc, char *argv[]) {
	int c;
	clock_t t = clock();
	hash_table *ht = NULL;

	fprintf(stderr, "%s \n", argv[1]);
	fprintf(stderr, "%s \n", argv[2]);
	fprintf(stderr, "%s \n", argv[3]);
	fprintf(stderr, "%s \n", argv[4]);
	ht = pe_load_hash(argv[1]);

	load_rm(ht, argv[2], argv[3], argv[4]);

	fprintf(stderr, "[pe_path] Done: %.2f sec\n", (float) (clock() - t)
			/ CLOCKS_PER_SEC);
	return 0;
}
