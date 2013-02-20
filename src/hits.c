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
	return h;
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
		i = 0;
		attr[0] = strtok(buf, " ");
		while (attr[i] != NULL) { //ensure a pointer was found
			attr[++i] = strtok(NULL, " "); //continue to tokenize the string
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
		h->strand = attr[8][0];
		h->qname = attr[9];
		h->q_size = atoi(attr[10]);
		h->q_start = atoi(attr[11]);
		h->q_end = atoi(attr[12]);
		h->tname = attr[13];
		h->t_size = atoi(attr[14]);
		h->t_start = atoi(attr[15]);
		h->t_end = atoi(attr[16]);
		h->block_count = atoi(attr[17]);
		h->block_sizes = (int*) calloc(h->block_count, sizeof(int));
		j = 0;
		index[0] = strtok(attr[18], " ");
		while (index[j] != NULL) { //ensure a pointer was found
			index[++j] = strtok(NULL, " "); //continue to tokenize the string
		}
		for (j = 0; j < h->block_count; j++) {
			h->block_sizes[j] = atoi(index[j]);
		}
		j = 0;
		index[0] = strtok(attr[19], " ");
		while (index[j] != NULL) { //ensure a pointer was found
			index[++j] = strtok(NULL, " "); //continue to tokenize the string
		}
		j = 0;
		index[0] = strtok(attr[20], " ");
		while (index[j] != NULL) { //ensure a pointer was found
			index[++j] = strtok(NULL, " "); //continue to tokenize the string
		}
		h->q_starts = (int*) calloc(h->block_count, sizeof(int));
		for (j = 0; j < h->block_count; j++) {
			h->q_starts[j] = atoi(index[j]);
		}
		h->t_starts = (int*) calloc(h->block_count, sizeof(int));
		for (j = 0; j < h->block_count; j++) {
			h->t_starts[j] = atoi(index[j]);
		}
		g_ptr_array_add(hits, h);
	}
	fclose(psl);
	return hits;
}
