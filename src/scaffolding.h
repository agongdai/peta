/*
 * scaffolding.h
 *
 *  Created on: 04-Feb-2013
 *      Author: carl
 */

#ifndef SCAFFOLDING_H_
#define SCAFFOLDING_H_

#include "edgelist.h"
#include <glib.h>

#ifdef __cplusplus
extern "C" {
#endif

#define MIN_OL_TO_SCA		100

typedef struct {
	int id;
	int comp_start;
	GPtrArray *edges; 		// Original edges after merging
	GPtrArray *hits;		// Blat hits of all edges
	GPtrArray *contigs;		// Contigs after breaking, ordering
	int alive;
} comp;

typedef struct {
	hash_table *ht;
	edgearray *all_edges;
	GPtrArray *comps;
	int insert_size;
	int sd_insert_size;
	int tid;
} comps_aux_t;

edgearray *scaffolding(edgearray *all_edges, const int insert_size,
		const int sd_insert_size, hash_table *ht, const int n_threads,
		char *psl_name);

#ifdef __cplusplus
}
#endif
#endif /* SCAFFOLDING_H_ */
