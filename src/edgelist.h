/*
 * arraylist.h
 *
 *  Created on: 13-Oct-2011
 *      Author: carl
 */

#ifndef EDGELIST_H_
#define EDGELIST_H_
#ifdef __cplusplus
extern "C" {
#endif

#include "edge.h"
#include "bwase.h"
#include "pealn.h"

#define LIST_SIZE				64
#define TRIVIAL_DIFF			4
#define GAP_OL					8
#define NOT_FOUND				-1
#define INVALID					-1
#define STOP_PAIRING			400
#define MIN_VALID_PAIRS			8
#define MISMATCHES				2
#define SHORT_MISMATCH			1
#define UPD_READS_THRE			10000

gint cmp_read_by_name(gpointer a, gpointer b);
void g_ptr_array_add_index(GPtrArray *array, gpointer data, const int index);
void g_ptr_array_replace_index(GPtrArray *array, gpointer data, const int index);
void g_ptr_array_replace_ptr(GPtrArray *array, gpointer data, gpointer olddata);
void g_ptr_array_iterator(gpointer value, gpointer user_data);
void g_ptr_array_uni_add(GPtrArray *array, gpointer data);
void g_ptr_array_concat(GPtrArray *array, GPtrArray *array_2);
void p_edgearray(const edgearray *array);
int edgearray_find(edgearray *array, edge *eg);
edge *edgearray_find_id(edgearray *array, const int ctg_id);
int readarray_find(readarray *array, bwa_seq_t *r);
int edgearray_find_similar(edgearray *array, edge *eg);
void adj_shift(edge *eg, const int trun_len);
readarray *get_paired_reads(readarray *ra_1, readarray *ra_2, bwa_seq_t *seqs);
readarray *find_unconditional_paired_reads(edge *eg_1, edge *eg_2, bwa_seq_t *seqs);
void merge_eg_to_left(edge *left_eg, edge *right_eg, const int gap);
void merge_eg_to_right(edge *left_eg, edge *right_eg, const int gap);
int get_mid_pos(readarray *ra, const int ori, const int lib_mean);
void upd_ctg_id(edge *eg, const int ctg_id, const int status);
void upd_reads_by_ol(bwa_seq_t *seqs, edge *eg, const int mismatches);
void upd_reads_by_ht(const hash_table *ht, edge *eg, const int mismatches, const int stage);
void upd_reads(const hash_table *ht, edge *eg, const int mismatches, const int stage);
int has_pairs_on_edge(edge *eg, bwa_seq_t *seqs, const int n_stop_pairs);
void log_reads(edgearray *ea);
void log_edge(const edge *eg);
void combine_reads(edge *left_eg, edge *right_eg, const int upd_shift, const int gap,
		const int ori);
void concat_readarray(readarray *left_reads, readarray *right_reads);
void clear_used_reads(edge *eg, const int reset_ctg_id);
void fill_in_hole(edge *ass_eg, edge *m_eg, const int ori, eg_gap *gap, const int nm, const int rl);
eg_gap *find_hole(edge *ass_eg, edge *m_eg, const int ori);
double *get_pair_dis_on_edge(edge *eg, int *n_pairs);
void readarray_add(edge *eg, bwa_seq_t *read);
void readarray_uni_add(edge *eg, bwa_seq_t *read);
void readarray_remove(edge *eg, bwa_seq_t *read);
void readarray_frozen(readarray *ra, bwa_seq_t *read);
void readarray_unfrozen(readarray *ra);
void rev_edge(edge *eg);
int binary_exists(const readarray *reads, const bwa_seq_t *read);
int has_most_fresh_reads(readarray *ra, const int max);
void rev_reads_pos(edge *eg);

#ifdef __cplusplus
}
#endif
#endif /* EDGELIST_H_ */
