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

#define INIT_PATH_CTGS 64
#define INIT_PATH_N 1024

typedef struct {
	int id;
	edgelist *edges;
	int n_ctgs;
	int len;
} rm_path;

void save_paths(const rm_path *paths, const int n_paths, FILE *ass_fa);
void rtv_back_paths(edge *eg, rm_path *path, rm_path *all_paths, int *n_paths,
		int *space, int level);
rm_path *rtv_rm_paths(const roadmap *rm, int *n_paths,
		const int left_max_ctg_id);
rm_path *rtv_paths(const roadmap *left_rm, const roadmap *right_rm,
		const int left_max_ctg_id, int *n_paths);
int app_ctg_to_path(rm_path *p, edge *contig);
rm_path *new_path();
void add_path(rm_path *all, rm_path *p, int *n_paths, int *space);
void p_path(const rm_path *p);

#ifdef __cplusplus
}
#endif
#endif /* PEPATH_H_ */
