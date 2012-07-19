/*
 * edge.h
 *
 *  Created on: 13-Oct-2011
 *      Author: carl
 */

#ifndef EDGE_H_
#define EDGE_H_
#include "bwase.h"
#include "glib.h"

#define MAX_N_EDGE_OUT 7
#define MAX_N_EDGE_IN 7
#define INIT_N_GAPS 3
#define INIT_N_READ_USED 7
#define INIT_N_READ_PAIRED 63

struct edge;

typedef GList edgelist;
typedef GPtrArray edgearray;
typedef GPtrArray readarray;
typedef struct edge edge;

struct edge {
	bwa_seq_t *contig;
	edgearray *in_egs;
	edgearray *out_egs;
	readarray *reads;
	char *name;
	edge *right_ctg; 	// If current contig is done before, record it
	edge *left_ctg;
	int r_shift; 		// shifted position of the right node
	int l_shift;
	int id; 			// contig id
	int len;
	int visited;
	int alive;
	int is_root;
	int ori;			// Orientation
	GPtrArray *gaps;
};

typedef struct {
	int s_index;
	int size;
	int ori;
} eg_gap;

eg_gap *init_gap(int s_index, int size, int ori);

#endif /* EDGE_H_ */
