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

int get_edges_ori(edge *eg_left, edge *eg_right);
int has_reads_in_common(edge *eg_1, edge *eg_2);
edge *merge_edges(edge *eg_1, edge *eg_2);
GPtrArray *scaffolding(GPtrArray *single_edges, const int insert_size, bwa_seq_t *seqs);
void merge_ol_edges(edgearray *single_edges, const int insert_size,
		bwa_seq_t *seqs);

#ifdef __cplusplus
}
#endif

#endif /* SCAFFOLDING_H_ */
