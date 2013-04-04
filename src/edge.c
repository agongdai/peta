/*
 * edge.c
 *
 *  Created on: 06-Mar-2012
 *      Author: carl
 */

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <glib.h>
#include "utils.h"
#include "edge.h"
#include "peseq.h"

eg_gap *init_gap(int s_index, int size, int ori) {
	eg_gap *gap = (eg_gap*) malloc (sizeof(eg_gap));
	gap->s_index = s_index;
	gap->size = size;
	gap->ori = ori;
	return gap;
}

void free_eg_gap(eg_gap *gap) {
	if (gap)
		free(gap);
}

edge *new_eg() {
	edge *eg = (edge*) malloc(sizeof(edge));
	eg->contig = 0;
	eg->in_egs = g_ptr_array_sized_new(0);
	eg->out_egs = g_ptr_array_sized_new(0);
	eg->reads = g_ptr_array_sized_new(0);
	eg->pairs = g_ptr_array_sized_new(0);
	eg->name = NULL;
	eg->right_ctg = NULL;
	eg->left_ctg = NULL;
	eg->len = 0;
	eg->r_shift = 0;
	eg->l_shift = 0;
	eg->id = 0;
	eg->visited = 0;
	eg->alive = 1;
	eg->is_root = 0;
	eg->ori = 0;
	eg->gaps = g_ptr_array_sized_new(0);
	eg->level = -1;
	eg->comp_id = -1;
	return eg;
}
