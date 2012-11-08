/*
 * pepath.h
 *
 *  Created on: Sep 6, 2011
 *      Author: xuyiling
 */

#ifndef PEPATH_H_
#define PEPATH_H_
#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include "roadmap.h"
#include "edgelist.h"
#include "bwase.h"
#include <glib.h>

#define INIT_PATH_CTGS 		64
#define INIT_PATH_N 		1024
#define MAX_ROADMAP_LEVEL	8
#define MIN_TX_LEN			100
#define PATH_MISMATCHES		4
#define SCORE_MATCH			2
#define SCORE_MISMATCH		-1
#define SCORE_GAP			-1

typedef struct {
	int id;
	GPtrArray *edges;
	int n_ctgs;
	int len;
	int alive;
	bwa_seq_t *seq;
	readarray *reads;
} rm_path;

GPtrArray *report_paths(edgearray *all_edges);
int pe_path(int argc, char *argv[]);
edgearray *load_rm(const hash_table *ht, const char *rm_dump_file,
		const char *rm_reads_file, const char *contig_file);

#ifdef __cplusplus
}
#endif
#endif /* PEPATH_H_ */
