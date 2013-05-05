/*
 * edge.h
 *
 *  Created on: 13-Oct-2011
 *      Author: carl
 */

#ifndef EDGE_H_
#define EDGE_H_
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

typedef GList edgelist;
typedef GPtrArray edgearray;
typedef GPtrArray readarray;
typedef struct edge edge;

struct edge {
	bwa_seq_t *contig;
	bwa_seq_t *r_tail;
	bwa_seq_t *l_tail;
	edgearray *in_egs;
	edgearray *out_egs;
	readarray *reads;
	readarray *pairs;
	char *name;
	edge *right_ctg; 	// If current contig is done before, record it
	edge *left_ctg;
	short r_shift; 		// shifted position of the right node
	short l_shift;
	uint64_t id; 			// contig id
	int len;
	int8_t visited;
	int8_t alive;
	int8_t is_root;
	int8_t ori;			// Orientation
	int16_t level;			// For post-processing
	uint64_t tid;			// Thread id
	int16_t comp_id;		// Component id
	uint64_t start_kmer_int;
	GPtrArray *gaps;
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
	void set_tail(edge *eg, const int shift, const int tail_len, const int ori);

#ifdef __cplusplus
}
#endif

#endif /* EDGE_H_ */
