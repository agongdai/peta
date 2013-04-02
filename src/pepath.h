/*
 * pepath.h
 *
 *  Created on: Sep 6, 2011
 *      Author: xuyiling
 */

#ifndef PEPATH_H_
#define PEPATH_H_

#include <stdio.h>
#include "roadmap.h"
#include "edgelist.h"
#include "bwtaln.h"
#include <glib.h>

#define INIT_PATH_CTGS 		64
#define INIT_PATH_N 		1024
#define MAX_ROADMAP_LEVEL	8
#define MIN_TX_LEN			100
#define PATH_MISMATCHES		4
#define PATH_MAX_GAPS		4
#define SCORE_MATCH			2
#define SCORE_MISMATCH		-1
#define SCORE_GAP			-2

typedef struct {
	int id;
	GPtrArray *edges;
	int n_ctgs;
	int len;
	int alive;
	bwa_seq_t *seq;
	readarray *reads;
} rm_path;

#ifdef __cplusplus
extern "C" {
#endif

rm_path *get_single_edge_path(edge *eg);
GPtrArray *report_paths(edgearray *all_edges, bwa_seq_t *seqs);
int pe_path(int argc, char *argv[]);
edgearray *load_rm(const hash_table *ht, const char *rm_dump_file,
		const char *rm_reads_file, const char *contig_file);
void save_paths(GPtrArray *paths, const char *tx_fn, const int min_len);

#ifdef __cplusplus
}
#endif

#endif /* PEPATH_H_ */
