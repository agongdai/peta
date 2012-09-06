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
#include <glib.h>

#define INIT_PATH_CTGS 64
#define INIT_PATH_N 1024

typedef struct {
	int id;
	GPtrArray *edges;
	int n_ctgs;
	int len;
} rm_path;

#ifdef __cplusplus
}
#endif
#endif /* PEPATH_H_ */
