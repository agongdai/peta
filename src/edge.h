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
	bwa_seq_t *ctg;			// Sequence
	bwa_seq_t *r_tail;		// Right tail, maybe virtual
	bwa_seq_t *l_tail;		// Left tail, maybe virtual
	edgearray *in_egs;		// Incoming edges
	edgearray *out_egs;		// Outgoing edges
	int id; 				// template id
	int32_t len;			// Length
	int8_t alive;			// Whether alive
	int8_t is_root;			// Whether it's a root node in the graph
	int8_t ori;				// Orientation
	int8_t in_connect;		// Indicates whether some template connects to it already
	uint64_t tid;			// Thread id
	int32_t comp_id;		// Component id
	uint64_t start_kmer;	// The starting kmer
	float coverage;			// Its kmer coverage
	uint32_t kmer_freq;		// Sum of all kmer frequencies
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
	void free_eg_seq(edge *eg);
	bwa_seq_t *cut_edge_tail(edge *eg, const int tail_len, const int pos, const int ori);
	void set_tail(edge *branch, edge *parent_eg, const int shift, const int tail_len, const int ori);
	void save_edges(edgearray *pfd_ctg_ids, FILE *ass_fa, const int ori,
			const int p_all, const int min_len);

#ifdef __cplusplus
}
#endif

#endif /* EDGE_H_ */
