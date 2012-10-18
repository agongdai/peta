/*
 * roadmap.h
 *
 *  Created on: Sep 2, 2011
 *      Author: carl
 */

#ifndef ROADMAP_H_
#define ROADMAP_H_
#ifdef __cplusplus
extern "C" {
#endif
#include <glib.h>
#include "bwase.h"
#include "pealn.h"
#include "pool.h"
#include "edgelist.h"

#ifndef EDGE_NODE
#define EDGE_NODE

#define INIT_PET_N 					1023
#define INIT_REMOVE_SIZE 			1023
#define SD_TIMES 					4
#define MEAN_TIMES 					2
#define PATH_RATIO 					2
#define RIGHT_CONTIG_CONNETER_COUNT	4

#endif

#ifndef RM
#define RM
typedef struct {
	edgearray *start_egs;
	int start_eg_n;
	int n_node;
	bwa_seq_t *pet_s;
	bwa_seq_t *pet_e;
} roadmap;
typedef struct {
	edgearray *egs;
} path;
#endif

roadmap *new_rm();
edge *new_eg();
void destroy_eg(edge *eg);
void free_readarray(readarray *ra);
void free_branch(edge *eg, const int ori, edgearray *all_edges);
void cut_connection(edge *ass_eg, edge *tmp_eg, const int ori);
void post_pro(edgearray *all_edges, const ass_opt *opt);
void dump_rm(edgearray *all_edges, char *rm_dump_file);

#ifdef __cplusplus
}
#endif

#endif /* ROADMAP_H_ */
