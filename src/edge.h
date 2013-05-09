/*
 * edge.h
 *
 *  Created on: 13-Oct-2011
 *      Author: carl
 */

#ifndef EDGE_H_
#define EDGE_H_
#include <stdio.h>
#include "bwtaln.h"
#include "glib.h"

#define MAX_N_EDGE_OUT 				7
#define MAX_N_EDGE_IN 				7
#define INIT_N_GAPS 				3
#define INIT_N_READ_USED 			7
#define INIT_N_READ_PAIRED 			63
#define UNUSED_CONTIG_ID			-1
#define INVALID_CONTIG_ID			-2	// Means a reads has been tried to extend, but fail, never use later

struct edge;

typedef GPtrArray edgearray;
typedef GPtrArray readarray;
typedef struct edge edge;

struct edge {
	bwa_seq_t *contig;
	bwa_seq_t *r_tail;
	bwa_seq_t *l_tail;
	edgearray *in_egs;
	edgearray *out_egs;
	int id; 			// contig id
	int32_t len;
	int8_t alive;
	int8_t is_root;
	int8_t ori;				// Orientation
	uint64_t tid;			// Thread id
	int32_t comp_id;		// Component id
	uint64_t start_kmer_int;
};

typedef struct {
	int s_index;
	int size;
	int ori;
} eg_gap;

#ifdef __cplusplus
extern "C" {
#endif

	eg_gap *init_gap(int s_index, int size, int ori);
	void free_eg_gap(eg_gap *gap);
	edge *new_eg();
	void reset_tid(edge *eg);
	void destroy_eg(edge *eg);
	bwa_seq_t *cut_edge_tail(edge *eg, const int tail_len, const int shift, const int ori);
	void set_tail(edge *branch, edge *parent_eg, const int shift, const int tail_len, const int ori);
	void save_edges(edgearray *pfd_ctg_ids, FILE *ass_fa, const int ori,
			const int p_all, const int min_len);

#ifdef __cplusplus
}
#endif

#endif /* EDGE_H_ */
