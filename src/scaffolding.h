/*
 * scaffolding.h
 *
 *  Created on: 18-Nov-2012
 *      Author: carl
 */

#ifndef SCAFFOLDING_H_
#define SCAFFOLDING_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include "roadmap.h"
#include "edgelist.h"
#include "bwase.h"
#include <glib.h>

#define MAX_EDGE_NM				11
#define EDGE_OL_THRE			11
#define EDGE_OL_PERC			0.2

int get_edges_ori(edge *eg_left, edge *eg_right, bwa_seq_t *seqs);
int has_reads_in_common(edge *eg_1, edge *eg_2);
edge *merge_edges(edge *eg_1, edge *eg_2, hash_table *ht);
void scaffolding(edgearray *single_edges, const int insert_size,
		const hash_table *ht, const int n_threads);
void merge_ol_edges(edgearray *single_edges, const int insert_size,
		const hash_table *ht, const int n_threads);

#ifdef __cplusplus
}
#endif

#endif /* SCAFFOLDING_H_ */
