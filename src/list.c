#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <ctype.h>
#include <string.h>
#include <glib.h>
#include "list.h"
#include "kmer.h"
#include "peseq.h"
#include "edgelist.h"

gint cmp_kmer_by_seq(gpointer a, gpointer b) {
	int64_t rs = 0;
	mer *read_a = *((mer**) a);
	mer *read_b = *((mer**) b);
	rs = (int64_t) read_a->s - (int64_t) read_b->s;
	if (rs == 0)
		return 0;
	if (rs > 0)
		return 1;
	return -1;
}

/*
 * Compare the kmer frequency by the frequency desc
 */
gint cmp_kmer_by_count(gpointer a, gpointer b) {
	mer *read_a = *((mer**) a);
	mer *read_b = *((mer**) b);
	return ((read_b->count) - (read_a->count));
}

slist *new_slist() {
	slist *sl = (slist*) malloc(sizeof(slist));
	sl->order = -1;
	sl->kmers = g_ptr_array_sized_new(0);
	sl->starts = g_ptr_array_sized_new(0);
	return sl;
}

void free_slist(slist *sl) {
	if (!sl)
		return;
	if (sl->starts)
		g_ptr_array_free(sl->starts, TRUE);
	if (sl->kmers)
		g_ptr_array_free(sl->kmers, TRUE);
	free(sl);
}

void slist_sorted_ins(mer *array, const uint32_t size, mer *data, const int index) {
	mer *s, *r = NULL;
	if (index < 0)
		return;
	s = &array[index];
	r = &array[index + 1];
	memmove(r, s, sizeof(mer) * (size - index));
	array[index] = *data;
}

/**
 * Binary search some value in the sorted list
 * Assumption: the order is increasing
 * Return: position of the hit query on the list;
 * 		   -1 if not found.
 */
int slist_binary(mer *kmers, const uint32_t size, uint64_t s) {
	int start = 0, end = size - 1, middle = 0;
	mer *m = NULL;
	while (start <= end) {
		middle = (start + end) / 2;
		m = &kmers[middle];
		if (s == m->s) {
			return middle;
		}
		if (s > m->s)
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
int slist_ins_pos(mer *kmers, const uint32_t size, mer *query) {
	int start = 0, end = size - 1, middle = 0;
	mer *v = NULL;
	if (size == 0)
		return 0;
	while (start <= end) {
		middle = (start + end) / 2;
		//show_debug_msg(__func__, "[%d, %d, %d] \n", start, middle, end);
		v = &kmers[middle];
		//show_debug_msg(__func__, "q [%" ID64 ", %d] \n", query->s, query->count);
		//show_debug_msg(__func__, "v [%" ID64 ", %d] \n", v->s, v->count);
		if (query->s == v->s)
			return middle;
		if (query->s > v->s)
			start = middle + 1;
		else
			end = middle - 1;
	}
	if (end < 0)
		return 0;
	if (start > size - 1)
		return size - 1;
	return start;
}

/**
 * Try to insert a value into the sorted list and keep the order
 * Assumption: the order is increasing
 * Return: if already exists, return 0;
 *         if inserted correctly, return 1.
 */
mer *slist_ins_pt(mer *kmers, uint32_t *size, uint32_t *full_space, mer *new_mer) {
	int pos = 0;
	if ((*size + 2) > *full_space) {
		*full_space += 2;
		kroundup32(*full_space);
		kmers = (mer*) realloc(kmers, sizeof(mer) * (*full_space));
	}
	pos = slist_ins_pos(kmers, *size, new_mer);
	slist_sorted_ins(kmers, *size, new_mer, pos);
	*size += 1;
	return kmers;
}
