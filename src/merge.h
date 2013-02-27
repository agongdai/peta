/*
 * scaffolding.h
 *
 *  Created on: 18-Nov-2012
 *      Author: carl
 */

#ifndef MERGE_H_
#define MERGE_H_

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
#define VAGUE_TAIL_LEN			100
#define EDGE_OL_PERC			0.2

enum ORDER_TYPE {
	HEAD_HEAD,
	HEAD_TAIL,
	TAIL_HEAD,
	TAIL_TAIL,
};

int get_edges_ori(edge *eg_left, edge *eg_right, bwa_seq_t *seqs);
int has_reads_in_common(edge *eg_1, edge *eg_2);
edge *merge_edges(edge *eg_1, edge *eg_2, hash_table *ht);
void merge_ol_edges(edgearray *single_edges, const int insert_size,
		const int sd_insert_size, hash_table *ht, const int n_threads);
void mark_sub_edge(edgearray *all_edges, GPtrArray *hits);
GPtrArray *get_probable_in_out(GPtrArray *all_edges, const int insert_size, edge *eg, bwa_seq_t *seqs);

#ifdef __cplusplus
}
#endif

#endif /* MERGE_H_ */
