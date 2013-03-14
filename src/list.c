#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <glib.h>
#include "list.h"
#include "edgelist.h"
#include "bwase.h"

gint cmp_kmer_by_seq(gpointer a, gpointer b) {
	int i = 0;
	bwa_seq_t *read_a = *((bwa_seq_t**) a);
	bwa_seq_t *read_b = *((bwa_seq_t**) b);
	show_debug_msg(__func__, "read_a: %p \n", a);
	show_debug_msg(__func__, "read_b: %p \n", b);
	show_debug_msg(__func__, "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n");
	show_debug_msg(__func__, "read_a: %p \n", read_a);
	show_debug_msg(__func__, "read_b: %p \n", read_b);
	for (i = 0; i < 30; i++) {
		show_debug_msg(__func__, "i = %d; a: %c \n", i, read_a->seq[i]);
		if (read_a->seq[i] < read_b->seq[i])
			return -1;
		if (read_a->seq[i] > read_b->seq[i])
			return 1;
	}
	return 0;
}

slist *new_slist(GCompareFunc cmp_f) {
	slist *sl = (slist*) malloc(sizeof(slist));
	sl->order = -1;
	sl->values = g_ptr_array_sized_new(0);
	sl->cmp_f = cmp_f;
	return sl;
}

void free_slist(slist *sl) {
	if (!sl)
		return;
	if (sl->values)
		g_ptr_array_free(sl->values, TRUE);
	free(sl);
}

/**
 * Binary search some value in the sorted list
 * Assumption: the order is increasing
 * Return: position of the hit query on the list;
 * 		   -1 if not found.
 */
int slist_binary(slist *sl, gpointer query) {
	int start = 0, end = sl->values->len - 1, middle = 0;
	int cmp_rs = 0;
	GCompareFunc f = &(sl->cmp_f);
	gpointer v = NULL;
	bwa_seq_t *r = NULL;
	if (sl->order == -1)
		g_ptr_array_sort(sl->values, sl->cmp_f);
	while (start <= end) {
		middle = (start + end) / 2;
		show_debug_msg(__func__, "%d, %d, %d \n", start, middle, end);
		r = g_ptr_array_index(sl->values, middle);
		show_debug_msg(__func__, "query: %p\n", query);
		show_debug_msg(__func__, "v: %p\n", r);
		p_query(__func__, query);
		p_query(__func__, v);
		cmp_rs = cmp_kmer_by_seq(query, (gpointer) r);
		if (cmp_rs == 0)
			return middle;
		if (cmp_rs > 0)
			start = middle + 1;
		else
			end = middle - 1;
	}
	return -1;
}

/**
 * Binary search some value in the sorted list
 * Assumption: the order is increasing
 * Return: if query found, return the index;
 *         if query not found, return the index where the query should be inserted into
 */
int slist_ins_pos(slist *sl, gpointer query) {
	int start = 0, end = sl->values->len - 1, middle = 0;
	int cmp_rs = 0;
	GCompareFunc f = sl->cmp_f;
	gpointer v = NULL;
	if (sl->order == -1)
		g_ptr_array_sort(sl->values, sl->cmp_f);
	while (start <= end) {
		middle = (start + end) / 2;
		v = g_ptr_array_index(sl->values, middle);
		cmp_rs = cmp_kmer_by_seq(query, v);
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
int slist_ins_pt(slist *sl, gpointer new_pt) {
	int pos = 0;
	if (slist_binary(sl, new_pt) >= 0)
		return 0;
	pos = slist_ins_pos(sl, new_pt);
	g_ptr_array_add_index(sl->values, new_pt, pos);
	return 1;
}
