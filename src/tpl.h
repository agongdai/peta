/*
 * tpl.h
 *
 *  Created on: 13-Oct-2011
 *      Author: carl
 */

#ifndef TPL_H_
#define TPL_H_
#include <stdio.h>
#include "bwtaln.h"
#include <glib.h>

#define MAX_N_tpl_OUT 				7
#define MAX_N_tpl_IN 				7
#define INIT_N_GAPS 				3
#define INIT_N_READ_USED 			7
#define INIT_N_READ_PAIRED 			63
#define UNUSED_CONTIG_ID			-1
#define INVALID_CONTIG_ID			-2	// Means a reads has been tried to extend, but fail, never use later

typedef GPtrArray tplarray;
typedef GPtrArray readarray;

typedef struct {
	bwa_seq_t *ctg;			// Sequence
	bwa_seq_t *r_tail;		// Right tail, maybe virtual
	bwa_seq_t *l_tail;		// Left tail, maybe virtual
	int id; 				// template id
	int32_t len;			// Length
	int8_t alive;			// Whether alive
	int8_t is_root;			// Whether it's a root node in the graph
	int8_t ori;				// Orientation
	int8_t in_connect;		// Indicates whether some template connects to it already
	uint64_t tid;			// Thread id
	uint64_t start_kmer;	// The starting kmer
	float coverage;			// Its kmer coverage
	uint32_t kmer_freq;		// Sum of all kmer frequencies
	GPtrArray *vertexes;	// The template be broken into vertexes linearly
	GPtrArray *reads;		// Reads on it. Used for paired validation
} tpl;

typedef struct {
	int s_index;
	int size;
	int ori;
} eg_gap;

#ifdef __cplusplus
extern "C" {
#endif

	gint cmp_tpl_by_id(gpointer a, gpointer b);
	eg_gap *init_gap(int s_index, int size, int ori);
	void free_eg_gap(eg_gap *gap);
	tpl *new_eg();
	void reset_tid(tpl *t);
	void destroy_eg(tpl *t);
	void free_eg_seq(tpl *t);
	bwa_seq_t *cut_tpl_tail(tpl *t, const int tail_len, const int pos, const int ori);
	void set_tail(tpl *branch, tpl *parent_eg, const int shift, const int tail_len, const int ori);
	void save_tpls(tplarray *pfd_ctg_ids, FILE *ass_fa, const int ori,
			const int p_all, const int min_len);
	int vld_tpl_mates(tpl *t1, tpl *t2, int start_2, int end_2,
			const int min_n_pairs);
	void add_read_to_tpl(tpl *t, bwa_seq_t *r, const int locus);
	int find_pairs(GPtrArray *reads_1, GPtrArray *reads_2, int t1_id, int t2_id,
			int start_2, int end_2, const int min_n_pairs);

#ifdef __cplusplus
}
#endif

#endif /* TPL_H_ */
