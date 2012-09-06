/*
 * pepath.c
 *
 *  Created on: Sep 6, 2011
 *      Author: Cai Shaojiang
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "pepath.h"
#include "roadmap.h"
#include "bwase.h"
#include "pechar.h"

int path_id = 0;

extern unsigned char nst_nt4_table[256];

void p_path(const rm_path *p) {
	int i = 0;
	edge *contig = NULL;
	bwa_seq_t *query;
	if (!p) {
		printf("Empty path! \n");
		return;
	}
	printf("[p_path] ----------------------------------------\n");
	printf("[p_path] Path %d (%p): %d \n", p->id, p, p->len);
	for (i = 0; i < p->edges->len; i++) {
		contig = g_ptr_array_index(p->edges, i);
		printf("[p_path] \t Contig %d [%s]: %d (%d, %d)\n", contig->id,
							query->name, contig->len,
							contig->right_ctg ? contig->right_ctg->id : 0,
							contig->r_shift);
	}
	printf("[p_path] ----------------------------------------\n");
}

