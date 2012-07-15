/*
 * arraylist.h
 *
 *  Created on: 13-Oct-2011
 *      Author: carl
 */

#ifndef ARRAYLIST_H_
#define ARRAYLIST_H_
#ifdef __cplusplus
extern "C" {
#endif

#include "edge.h"
#include "bwase.h"

#define LIST_SIZE		64
#define TRIVIAL_DIFF	4
#define GAP_OL			8
#define NOT_FOUND		-1
#define INVALID			-1

void g_ptr_array_add_index(GPtrArray *array, gpointer data, const int index);
void g_ptr_array_replace_index(GPtrArray *array, gpointer data, const int index);
void g_ptr_array_replace_ptr(GPtrArray *array, gpointer data, gpointer olddata);
void g_ptr_array_iterator(gpointer value, gpointer user_data);
void g_ptr_array_uni_add(GPtrArray *array, gpointer data);
void g_ptr_array_concat(GPtrArray *array, GPtrArray *array_2);
void p_edgearray(const edgearray *array);
int edgearray_find(edgearray *array, edge *eg);
int readarray_find(readarray *array, bwa_seq_t *r);
int edgearray_find_similar(edgearray *array, edge *eg);
readarray *get_paired_reads(readarray *ra_1, readarray *ra_2, bwa_seq_t *seqs, const int ori);
void merge_eg_to_left(edge *left_eg, edge *right_eg, const int gap);
void merge_eg_to_right(edge *left_eg, edge *right_eg, const int gap);
int get_mid_pos(readarray *ra, const int ori, const int lib_mean);
void upd_reads(edge *eg, const int mismatches);
void log_reads(edgearray *ea);
void log_edge(const edge *eg);
void combine_reads(edge *left_eg, edge *right_eg, const int upd_shift, const int gap,
		const int ori);
void clear_used_reads(edge *eg, const int reset_ctg_id);
void fill_in_hole(edge *ass_eg, edge *m_eg, const int ori, eg_gap *gap, const int nm, const int rl);

#ifdef __cplusplus
}
#endif
#endif /* ARRAYLIST_H_ */
