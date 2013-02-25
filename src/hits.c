/*
 * hits.c
 *
 *  Created on: 18-Feb-2013
 *      Author: carl
 */
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <glib.h>
#include "hits.h"
#include "utils.h"
#include "edgelist.h"

blat_hit *new_hit() {
	blat_hit *h = (blat_hit*) malloc(sizeof(blat_hit));
	h->matches = 0;
	h->mismatches = 0;
	h->n_rep_matches = 0;
	h->n_ns = 0;
	h->q_gap_count = 0;
	h->q_gap_bases = 0;
	h->t_gap_count = 0;
	h->t_gap_bases = 0;
	h->strand = '+';
	h->qname = NULL;
	h->q_size = 0;
	h->q_start = 0;
	h->q_end = 0;
	h->tname = NULL;
	h->t_size = 0;
	h->t_start = 0;
	h->t_end = 0;
	h->block_count = 0;
	h->block_sizes = NULL;
	h->q_starts = NULL;
	h->t_starts = NULL;
	h->visited = 0;
	h->alen = 0;
	return h;
}

void p_hit(blat_hit *h) {
	int i = 0;
	printf("%d\t", h->matches);
	printf("%d\t", h->mismatches);
	printf("%d\t", h->n_rep_matches);
	printf("%d\t", h->n_ns);
	printf("%d\t", h->q_gap_count);
	printf("%d\t", h->q_gap_bases);
	printf("%d\t", h->t_gap_count);
	printf("%d\t", h->t_gap_bases);
	printf("%c\t", h->strand);
	printf("%s\t", h->qname);
	printf("%d\t", h->q_size);
	printf("%d\t", h->q_start);
	printf("%d\t", h->q_end);
	printf("%s\t", h->tname);
	printf("%d\t", h->t_size);
	printf("%d\t", h->t_start);
	printf("%d\t", h->t_end);
	printf("%d\t", h->block_count);
	for (i = 0; i < h->block_count; i++) {
		printf("%d,", h->block_sizes[i]);
	}
	printf("\t");
	for (i = 0; i < h->block_count; i++) {
		printf("%d,", h->q_starts[i]);
	}
	printf("\t");
	for (i = 0; i < h->block_count; i++) {
		printf("%d,", h->t_starts[i]);
	}
	printf("\t\n");
}

gint cmp_hit_by_qname(gpointer a, gpointer b) {
	blat_hit *hit_a = *((blat_hit**) a);
	blat_hit *hit_b = *((blat_hit**) b);
	return (atoi(hit_a->qname) - atoi(hit_b->qname));
}

GPtrArray *read_blat_hits(const char *psl_file) {
	char buf[BUFSIZ];
	char *attr[64], *index[16];
	FILE *psl = NULL;
	int i = 0, line_no = 0, j = 0;
	blat_hit *h = NULL;
	GPtrArray *hits = g_ptr_array_sized_new(32);

	psl = xopen(psl_file, "r");
	while (fgets(buf, sizeof(buf), psl)) {
		line_no += 1;
		if (line_no <= 5)
			continue;
		//show_debug_msg(__func__, "%s", buf);
		i = 0;
		attr[0] = strtok(buf, "\t");
		while (attr[i] != NULL) { //ensure a pointer was found
			attr[++i] = strtok(NULL, "\t"); //continue to tokenize the string
		}
		h = new_hit();
		h->matches = atoi(attr[0]);
		h->mismatches = atoi(attr[1]);
		h->n_rep_matches = atoi(attr[2]);
		h->n_ns = atoi(attr[3]);
		h->q_gap_count = atoi(attr[4]);
		h->q_gap_bases = atoi(attr[5]);
		h->t_gap_count = atoi(attr[6]);
		h->t_gap_bases = atoi(attr[7]);
		h->strand = strdup(attr[8])[0];
		h->qname = strdup(attr[9]);
		h->q_size = atoi(attr[10]);
		h->q_start = atoi(attr[11]);
		h->q_end = atoi(attr[12]);
		h->tname = strdup(attr[13]);
		h->t_size = atoi(attr[14]);
		h->t_start = atoi(attr[15]);
		h->t_end = atoi(attr[16]);
		h->block_count = atoi(attr[17]);
		h->block_sizes = (int*) calloc(h->block_count, sizeof(int));
		j = 0;
		index[0] = strtok(attr[18], ",");
		while (index[j] != NULL) { //ensure a pointer was found
			index[++j] = strtok(NULL, ","); //continue to tokenize the string
		}
		for (j = 0; j < h->block_count; j++) {
			h->block_sizes[j] = atoi(index[j]);
			h->alen += h->block_sizes[j];
		}
		j = 0;
		index[0] = strtok(attr[19], ",");
		while (index[j] != NULL) { //ensure a pointer was found
			index[++j] = strtok(NULL, ","); //continue to tokenize the string
		}
		h->q_starts = (int*) calloc(h->block_count, sizeof(int));
		for (j = 0; j < h->block_count; j++) {
			h->q_starts[j] = atoi(index[j]);
		}
		j = 0;
		index[0] = strtok(attr[20], ",");
		while (index[j] != NULL) { //ensure a pointer was found
			index[++j] = strtok(NULL, ","); //continue to tokenize the string
		}
		h->t_starts = (int*) calloc(h->block_count, sizeof(int));
		for (j = 0; j < h->block_count; j++) {
			h->t_starts[j] = atoi(index[j]);
		}
		if (strcmp(h->qname, h->tname) == 0 && h->block_count == 1
				&& h->block_sizes[0] == h->t_size) {
			free(h);
			continue;
		}
		g_ptr_array_add(hits, h);
		//p_hit(h);
		//break;
	}
	fclose(psl);
	return hits;
}

int realign_thread(gpointer e, gpointer data) {
	edge *eg = (edge*) e;
	hash_table *ht = (hash_table*) data;
	realign_reads_by_ht(ht, eg, MISMATCHES);
	return 0;
}

void realign_by_blat(edgearray *all_edges, hash_table *ht, const int n_threads) {
	GThreadPool *thread_pool = NULL;
	int i = 0;
	edge *eg = NULL;
	thread_pool = g_thread_pool_new((GFunc) realign_thread, ht,
			n_threads, TRUE, NULL);
	show_msg(__func__, "Realigning all reads to %d edges... \n", all_edges->len);
	for (i = 0; i < all_edges->len; i++) {
		eg = g_ptr_array_index(all_edges, i);
		g_thread_pool_push(thread_pool, (gpointer) eg, NULL);
	}
	g_thread_pool_free(thread_pool, 0, 1);
}
