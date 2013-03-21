#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <glib.h>
#include "list.h"
#include "bwtaln.h"
#include "peseq.h"
#include "edgelist.h"

gint cmp_kmer_by_seq(gpointer a, gpointer b) {
	int i = 0;
	bwa_seq_t *read_a = *((bwa_seq_t**) a);
	bwa_seq_t *read_b = *((bwa_seq_t**) b);
//	return (strcmp(read_a->seq, read_b->seq));
	for (i = 0; i < read_a->len; i++) {
		if (read_a->seq[i] < read_b->seq[i])
			return -1;
		if (read_a->seq[i] > read_b->seq[i])
			return 1;
	}
	return 0;
}

/*
 * Compare the kmer frequency by the frequency desc, which is stored as seq name
 */
gint cmp_kmer_by_freq(gpointer a, gpointer b) {
	int i = 0;
	bwa_seq_t *read_a = *((bwa_seq_t**) a);
	bwa_seq_t *read_b = *((bwa_seq_t**) b);
	return (atoi(read_b->name) - atoi(read_a->name));
}

int cmp_seq_kmer(bwa_seq_t *a, bwa_seq_t *b) {
	int i = 0;
	for (i = 0; i < a->len; i++) {
		if (a->seq[i] < b->seq[i])
			return -1;
		if (a->seq[i] > b->seq[i])
			return 1;
	}
	return 0;
}

slist *new_slist() {
	slist *sl = (slist*) malloc(sizeof(slist));
	sl->order = -1;
	sl->values = g_ptr_array_sized_new(0);
	sl->kmers = g_ptr_array_sized_new(0);
	return sl;
}

void free_slist(slist *sl) {
	if (!sl)
		return;
	if (sl->values)
		g_ptr_array_free(sl->values, TRUE);
	if (sl->kmers)
		g_ptr_array_free(sl->kmers, TRUE);
	free(sl);
}

/**
 * Binary search some value in the sorted list
 * Assumption: the order is increasing
 * Return: position of the hit query on the list;
 * 		   -1 if not found.
 */
int slist_binary(slist *sl, bwa_seq_t *query) {
	int start = 0, end = sl->values->len - 1, middle = 0;
	int cmp_rs = 0;
	bwa_seq_t *r = NULL;
	while (start <= end) {
		middle = (start + end) / 2;
		r = g_ptr_array_index(sl->values, middle);
		cmp_rs = cmp_seq_kmer(query, r);
		//show_debug_msg(__func__, "[Start, Middle, End]: [%d, %d, %d] \n", start, middle, end);
		//p_query(__func__, r);
		if (cmp_rs == 0) {
			//show_debug_msg(__func__, "Hit %d !\n", middle);
			return middle;
		}
		if (cmp_rs > 0)
			start = middle + 1;
		else
			end = middle - 1;
	}
	//show_debug_msg(__func__, "Not Hit!\n");
	return -1;
}

/**
 * Binary search some value in the sorted list
 * Assumption: the order is increasing
 * Return: if query found, return the index;
 *         if query not found, return the index where the query should be inserted into
 */
int slist_ins_pos(slist *sl, bwa_seq_t *query) {
	int start = 0, end = sl->values->len - 1, middle = 0;
	int cmp_rs = 0;
	bwa_seq_t *v = NULL;
	while (start <= end) {
		middle = (start + end) / 2;
		v = g_ptr_array_index(sl->values, middle);
		cmp_rs = cmp_seq_kmer(query, v);
		if (cmp_rs == 0)
			return middle;
		if (cmp_rs > 0)
			start = middle + 1;
		else
			end = middle - 1;
	}
	if (end < 0)
		return 0;
	if (start > sl->values->len - 1)
		return sl->values->len - 1;
	return start;
}

/**
 * Try to insert a value into the sorted list and keep the order
 * Assumption: the order is increasing
 * Return: if already exists, return 0;
 *         if inserted correctly, return 1.
 */
int slist_ins_pt(slist *sl, bwa_seq_t *new_pt) {
	int pos = 0;
	if (slist_binary(sl, new_pt) >= 0)
		return 0;
	pos = slist_ins_pos(sl, new_pt);
	g_ptr_array_add_index(sl->values, new_pt, pos);
	return 1;
}
