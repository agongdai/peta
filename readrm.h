/*
 * readrm.h
 *
 *  Created on: Nov 30, 2011
 *      Author: carl
 */

#ifndef READRM_H_
#define READRM_H_

#include "edge.h"
#include "roadmap.h"
#include "pealn.h"
#include "stdio.h"

void p_flat_eg(const edge *eg);
void p_eg(const edge *root, int level);
void p_rm(const roadmap *rm);
void graph_rm(const roadmap *rm, char *dotfile);
void graph_by_edges(const edgearray *all_edges, char *dotfile);
void save_node(const edge *eg, FILE *ass_fa);
void save_tx(const roadmap *rm, FILE *ass_fa);
void save_edges(edgearray *pfd_ctg_ids, FILE *ass_fa, const int max_ctg_id,
		const int p_all, const int min_len);
int is_start_eg(const roadmap *rm, const edge *eg);
int has_edge(const edge *eg, const int contig_id, const int in_out);
int is_sbl(const edge *eg_1, const edge *eg_2);
edgearray *get_sbls(const edge *eg);

#endif /* READRM_H_ */
