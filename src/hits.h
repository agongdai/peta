/*
 * hits.h
 *
 *  Created on: 18-Feb-2013
 *      Author: carl
 */

#ifndef HITS_H_
#define HITS_H_

#include <glib.h>
#include "pehash.h"
#include "edgelist.h"

typedef struct {
	int matches;
	int mismatches;
	int n_rep_matches;
	int n_ns;
	int q_gap_count;
	int q_gap_bases;
	int t_gap_count;
	int t_gap_bases;
	char strand;
	char *qname;
	int q_size;
	int q_start;
	int q_end;
	char *tname;
	int t_size;
	int t_start;
	int t_end;
	int block_count;
	int *block_sizes;
	int *q_starts;
	int *t_starts;
	int visited;
	int alen;
} blat_hit;

GPtrArray *read_blat_hits(const char *psl_file);
void p_hit(blat_hit *h);
gint cmp_hit_by_qname(gpointer a, gpointer b);
void realign_by_blat(edgearray *all_edges, hash_table *ht, const int n_threads);
void free_blat_hits(GPtrArray *hits);

#endif /* HITS_H_ */
