/*
 * arraylist.c
 *
 *  Created on: 13-Oct-2011
 *      Author: carl
 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <glib.h>
#include "edgelist.h"
#include "peseq.h"
#include "utils.h"

gint cmp_read_by_shift(gpointer a, gpointer b) {
	bwa_seq_t *read_a = *((bwa_seq_t**) a);
	bwa_seq_t *read_b = *((bwa_seq_t**) b);
	return (read_a->shift - read_b->shift);
}

gint cmp_read_by_name(gpointer a, gpointer b) {
	bwa_seq_t *read_a = *((bwa_seq_t**) a);
	bwa_seq_t *read_b = *((bwa_seq_t**) b);
	return (atoi(read_a->name) - atoi(read_b->name));
}

void g_ptr_array_add_index(GPtrArray *array, gpointer data, const int index) {
	gpointer s, r;
	if (index < 0 || index > array->len)
		return;
	g_ptr_array_add(array, NULL); // Spare one more space
	if (array->len > 0 && array->len > index) { // If append to the end or the the beginning, not done here.
		s = &array->pdata[index];
		r = &array->pdata[index + 1];
		memmove(r, s, sizeof(gpointer) * (array->len - 1 - index));
	}
	array->pdata[index] = data;
}

void g_ptr_array_replace_index(GPtrArray *array, gpointer data, const int index) {
	if (!array || index < 0 || index >= array->len)
		return;
	array->pdata[index] = data;
}

void g_ptr_array_iterator(gpointer value, gpointer user_data) {
	edge *eg = (edge*) value;
	printf("[iterator] ");
	printf(" => %d \n", eg->id);
}

void p_edgearray(const edgearray *array) {
	edge *eg;
	int i = 0;
	printf("[p_edgearray] ---------------------------- \n");
	printf("[p_edgearray] Array Is: \n");
	for (i = 0; i < array->len; i++) {
		eg = g_ptr_array_index(array, i);
		if (eg)
			printf("[p_edgearray] %d => %d, len=%d \n", i, eg->id, eg->len);
		else
			printf("[p_edgearray] %d => Nil \n", i);
	}
	printf("[p_edgearray] ---------------------------- \n");
}

void g_ptr_array_replace_ptr(GPtrArray *array, gpointer data, gpointer olddata) {
	int i = 0;
	if (!array || !olddata)
		return;
	for (i = 0; i < array->len; i++) {
		if (g_ptr_array_index(array, i) == olddata) {
			g_ptr_array_replace_index(array, data, i);
			return;
		}
	}
}

void g_ptr_array_uni_add(GPtrArray *array, gpointer data) {
	int i = 0;
	if (!array || !data)
		return;
	for (i = 0; i < array->len; i++) {
		if (g_ptr_array_index(array, i) == data)
			return;
	}
	g_ptr_array_add(array, data);
}

void g_ptr_array_concat(GPtrArray *array, GPtrArray *array_2) {
	gpointer data = NULL;
	int i = 0;
	if (!array_2)
		return;
	for (i = 0; i < array_2->len; i++) {
		data = g_ptr_array_index(array_2, i);
		g_ptr_array_uni_add(array, data);
	}
}

int edgearray_find(edgearray *array, edge *eg) {
	int i = 0;
	if (!array || !eg)
		return NOT_FOUND;
	for (i = 0; i < array->len; i++) {
		if (g_ptr_array_index(array, i) == eg) {
			return i;
		}
	}
	return NOT_FOUND;
}

edge *edgearray_find_id(edgearray *array, const int ctg_id) {
	edge *to_return = NULL;
	int i = 0;
	if (!array || array->len <= 0 || ctg_id < 0)
		return NULL;
	for (i = 0; i < array->len; i++) {
		to_return = g_ptr_array_index(array, i);
		if (to_return->id == ctg_id)
			return to_return;
	}
	return NULL;
}

int readarray_find(readarray *array, bwa_seq_t *r) {
	int i = 0;
	if (!array || !r)
		return NOT_FOUND;
	for (i = 0; i < array->len; i++) {
		if (g_ptr_array_index(array, i) == r) {
			return i;
		}
	}
	return NOT_FOUND;
}

int edgearray_find_similar(edgearray *array, edge *eg) {
	edge *eg_i, *eg_right, *eg_i_right;
	int i = 0, diff = 0;
	if (!array || !eg)
		return NOT_FOUND;
	for (i = 0; i < array->len; i++) {
		eg_i = g_ptr_array_index(array, i);
		eg_right = eg->right_ctg;
		eg_i_right = eg_i->right_ctg;
		if ((eg_right && !eg_i_right) || (!eg_right && eg_i_right))
			continue;
		if (eg_right && eg_i_right && eg_right->id != eg_i_right->id)
			continue;
		diff = abs((eg->r_shift - eg->len) - (eg_i->r_shift - eg_i->len));
		if (diff < TRIVIAL_DIFF)
			return i;
		if (eg_right && eg_i_right) {
			diff = abs((eg->r_shift - eg_right->len) - (eg_i->r_shift
					- eg_i_right->len));
			if (diff < TRIVIAL_DIFF)
				return i;
		}
	}
	return NOT_FOUND;
}

void adj_shift(edge *eg, const int trun_len) {
	int i = 0;
	bwa_seq_t *r;
	if (!eg || trun_len < 0)
		return;
	for (i = 0; i < eg->reads->len; i++) {
		r = g_ptr_array_index(eg->reads, i);
		r->shift -= trun_len;
	}
}

double *get_pairs_on_edge(edge *eg, int *n_pairs) {
	int i = 0, index = 0;
	bwa_seq_t *s = NULL, *next_s = NULL;
	readarray *reads = eg->reads;
	double *pairs = (double*) calloc(reads->len / 2 + 1, sizeof(double));
	if (!reads || reads->len <= 0)
		return pairs;
	for (i = 0; i < reads->len - 1; i++) {
		s = g_ptr_array_index(reads, i);
		index = atoi(s->name);
		if (index % 2 == 0) {
			next_s = g_ptr_array_index(reads, i + 1);
			if (atoi(next_s->name) - index == 1) {
				i++;
				//				p_query("Left", s);
				//				p_query("Right", next_s);
				pairs[*n_pairs] = abs(next_s->shift - s->shift);
				*n_pairs += 1;
			}
		}
	}
	return pairs;
}

/**
 * ra_1 is at the left, and ra_2 is at the right
 */
readarray *get_paired_reads(readarray *ra_1, readarray *ra_2, bwa_seq_t *seqs) {
	readarray *paired = g_ptr_array_sized_new(INIT_N_READ_PAIRED);
	int i = 0, correct_dir = 0;
	bwa_seq_t *read_1, *read_2;
	if (!ra_1 || !ra_2 || ra_1->len == 0 || ra_2->len == 0)
		return paired;
	for (i = 0; i < ra_1->len; i++) {
		read_1 = g_ptr_array_index(ra_1, i);
		correct_dir = read_1->rev_com ? is_right_mate(read_1->name)
				: is_left_mate(read_1->name);
		if (correct_dir) {
			read_2 = get_mate(read_1, seqs);
			if (readarray_find(ra_2, read_1) != NOT_FOUND || readarray_find(
					ra_1, read_2) != NOT_FOUND)
				continue;
			//			p_query("PAIR", read_1);
			//			p_query("PAIR", read_2);
			g_ptr_array_add(paired, read_1);
		}
	}
	return paired;
}

/**
 * Assume large_ra and small_ra are sorted by name increasingly.
 */
int has_counter_pair(readarray *large_ra, readarray *small_ra,
		const int lower_bound, const int upper_bound) {
	int i = 0, j = 0, ori = 0, pair_ori = 0, matched_pairs = 0;
	int large_cursor = 0;
	bwa_seq_t *s_read, *l_read;
	for (i = 0; i < small_ra->len; i++) {
		s_read = g_ptr_array_index(small_ra, i);
		for (j = large_cursor; j < large_ra->len; j++) {
			l_read = g_ptr_array_index(large_ra, j);
			if (atoi(l_read->name) - atoi(s_read->name) > 1)
				break;
			large_cursor = j;
			if ((l_read->shift < 0 || l_read->shift >= lower_bound) // Shift could be negative
					&& l_read->shift <= upper_bound && strcmp(s_read->name,
					get_mate_name(l_read->name)) == 0) {
				if (!s_read->rev_com && is_left_mate(s_read->name)) {
					pair_ori = 1; // The overall orientation should be: large_ra.....small_ra
				} else
					pair_ori = 2;
				break;
			}
		}
		if (!ori) { // If the overall ori is not set, set it.
			ori = pair_ori; // If there is no mate, the pair_ori remains 0
		} else {
			if (pair_ori != ori) {
				return 1;
			}
			matched_pairs++;
		}
	}
	return 0;
}

/**
 * Return value of the eg_gap:
 * 1. NULL: no place to put the flying edge;
 * 2. (-1, -1, ori): put the flying edge to the end or the beginning
 * 3. others: put into some existing hole
 */
eg_gap *find_hole(edge *ass_eg, edge *m_eg, const int ori) {
	eg_gap *gap = NULL, *hole = NULL;
	int i = 0, lower = 0, upper = ass_eg->len;
	int left_ol_len = 0, right_ol_len = 0;
	readarray *ra = ass_eg->reads;
	readarray *m_ra = m_eg->reads;
	bwa_seq_t *left_seq = 0;
	bwa_seq_t *right_seq = 0;
	// The edge m_eg could be put at the end (ori == 0) or at the beginning (ori == 1)
	show_debug_msg(__func__, "Looking for holes... \n");
	g_ptr_array_sort(ra, (GCompareFunc) cmp_read_by_name);
	g_ptr_array_sort(m_ra, (GCompareFunc) cmp_read_by_name);
	// If there is no counter read pairs
	if (has_counter_pair(ra, m_ra, lower, upper)) {
		show_debug_msg(__func__, "Has counter pairs:");
		// p_readarray(ra, 1);
		// p_readarray(m_ra, 1);
		// The edge m_eg should be put into some hole
		for (i = ass_eg->gaps->len - 1; i >= 0; i--) {
			gap = g_ptr_array_index(ass_eg->gaps, i);
			show_debug_msg(__func__, "Found a hole: %d + %d, ori %d \n",
					gap->s_index, gap->size, gap->ori);
			lower = ori ? gap->s_index : 0;
			upper = ori ? ass_eg->len : gap->s_index;
			if (!has_counter_pair(ra, m_ra, lower, upper)) {
				hole = gap;
				break;
			}
		}
	} else {
		hole = init_gap(-1, -1, ori);
	}

	// Check overlapping when there is no counter pair and no hole could be put
	// For example, ass_eg is: aaattccgg......tctcagtg
	// 				m_eg is:        ccggagac
	// Here left part of m_eg 'ccgg' overlaps with 'aaattccgg',
	// Then truncate m_eg, put m_eg into the hole
	if (hole && hole->size == -1) {
		for (i = ass_eg->gaps->len - 1; i >= 0; i--) {
			gap = g_ptr_array_index(ass_eg->gaps, i);
			if (ori) {
				left_seq = new_seq(ass_eg->contig, ass_eg->len - gap->s_index
						- gap->size, 0);
				right_seq = new_seq(ass_eg->contig, gap->s_index, ass_eg->len
						- gap->s_index);
			} else {
				left_seq = new_seq(ass_eg->contig, gap->s_index, 0);
				right_seq = new_seq(ass_eg->contig, ass_eg->len - gap->size
						- gap->s_index, gap->size + gap->s_index);
			}
			left_ol_len = find_ol(left_seq, m_eg->contig, 1);
			right_ol_len = find_ol(m_eg->contig, right_seq, 1);
			//			p_query("Left: ", left_seq);
			//			p_query("Middle: ", m_eg->contig);
			//			p_query("Right: ", right_seq);
			bwa_free_read_seq(1, left_seq);
			bwa_free_read_seq(1, right_seq);
			// If overlapped, return current hole
			if (left_ol_len > MIN_OL || right_ol_len > MIN_OL) {
				hole = gap;
				break;
			}
		}
	}
	g_ptr_array_sort(ra, (GCompareFunc) cmp_read_by_shift);
	g_ptr_array_sort(m_ra, (GCompareFunc) cmp_read_by_shift);
	return hole;
}

void fill_in_hole(edge *ass_eg, edge *m_eg, const int ori, eg_gap *gap,
		const int nm, const int rl) {
	int start = 0, new_len = 0, ol_len = 0, left_ol = 0, right_ol = 0;
	bwa_seq_t *ass_seq = ass_eg->contig, *m_seq = m_eg->contig;
	bwa_seq_t *left_seq, *right_seq;
	if (ori) { // Extending to the left
		left_seq
				= new_seq(ass_seq, (ass_eg->len - gap->s_index - gap->size), 0);
		right_seq = new_seq(ass_seq, gap->s_index, ass_eg->len - gap->s_index);
		ol_len = find_ol(left_seq, m_seq, nm);
		if (ol_len > GAP_OL && ol_len < rl + GAP_OL) {
			trun_seq(m_seq, ol_len);
			left_ol = 1;
		}
		ol_len = find_ol(m_seq, right_seq, nm);
		if (ol_len > GAP_OL && ol_len < rl + GAP_OL) {
			m_seq->len -= ol_len;
			m_seq->seq[m_seq->len] = '\0';
			right_ol = 1;
		}
		m_eg->len = m_seq->len;
		// Just replace the gap with the new edge
		if (left_ol && right_ol) {
			if (m_seq->len > gap->size) {
				new_len = (ass_eg->len - gap->size + m_eg->len);
				if (new_len > ass_seq->full_len) {
					kroundup32(new_len);
					ass_seq->seq = (ubyte_t*) realloc(ass_seq->seq,
							sizeof(ubyte_t) * new_len);
					ass_seq->full_len = new_len;
				}
			}
			// Move the right part to the correct position
			memmove(&ass_seq->seq[(ass_seq->len - gap->s_index - gap->size)
					+ m_eg->len], &ass_seq->seq[ass_seq->len - gap->s_index],
					right_seq->len);
			memcpy(&ass_seq->seq[ass_seq->len - gap->s_index - gap->size],
					m_seq->seq, m_seq->len);
			ass_seq->len -= gap->size - m_seq->len;
			ass_seq->seq[ass_seq->len] = '\0';
			ass_eg->len = ass_seq->len;
		} else { // Put the edge somewhere in the middle of the gap
			if (m_eg->len > gap->size) { // If the new edge length is larger than the gap, remove the gap
				// spare some more space for it.
				new_len = ass_seq->len + (m_eg->len - gap->size);
				if (new_len > ass_seq->full_len) {
					ass_seq->len = new_len;
					kroundup32(new_len);
					ass_seq->seq = (ubyte_t*) realloc(ass_seq->seq,
							sizeof(ubyte_t) * new_len);
					ass_seq->full_len = new_len;
					ass_seq->seq[ass_seq->len] = '\0';
				}
				// Here the length of contig is changed, but length of edge is not!
				memmove(&ass_seq->seq[(ass_eg->len - gap->s_index - gap->size)
						+ m_eg->len],
						&ass_seq->seq[ass_eg->len - gap->s_index],
						right_seq->len);
				start = ass_eg->len - gap->s_index - gap->size;
				ass_eg->len = ass_seq->len; // Set them to be correct
			} else {
				if (left_ol)
					start = ass_seq->len - gap->s_index - gap->size;
				if (right_ol)
					start = (ass_seq->len - gap->s_index - gap->size)
							+ (gap->size - m_eg->len);
				if (!left_ol && !right_ol)
					start = (ass_seq->len - gap->s_index - gap->size)
							+ (gap->size - m_eg->len) / 2;
			}
			memcpy(&ass_seq->seq[start], m_seq->seq, m_eg->len);
		}
	} else { // Extending to the right
		left_seq = new_seq(ass_seq, gap->s_index, 0);
		right_seq = new_seq(ass_seq, ass_eg->len - gap->size - gap->s_index,
				gap->size + gap->s_index);
		ol_len = find_ol(left_seq, m_seq, nm);
		if (ol_len > GAP_OL && ol_len < rl + GAP_OL) {
			trun_seq(m_seq, ol_len);
			left_ol = 1;
		}
		ol_len = find_ol(m_seq, right_seq, nm);
		if (ol_len > GAP_OL && ol_len < rl + GAP_OL) {
			m_seq->len -= ol_len;
			m_seq->seq[m_seq->len] = '\0';
			right_ol = 1;
		}
		m_eg->len = m_seq->len;
		// Just replace the gap with the new edge
		if (left_ol && right_ol) {
			if (m_seq->len > gap->size) {
				new_len = (ass_eg->len - gap->size + m_eg->len);
				if (m_seq->full_len < new_len) {
					kroundup32(new_len);
					ass_seq->seq = (ubyte_t*) realloc(ass_seq->seq,
							sizeof(ubyte_t) * new_len);
					m_seq->full_len = new_len;
				}
			}
			// Move the right part to the correct position
			memmove(&ass_seq->seq[gap->s_index + m_eg->len],
					&ass_seq->seq[gap->s_index + gap->size], right_seq->len);
			memcpy(&ass_seq->seq[gap->s_index], m_seq->seq, m_seq->len);
			ass_seq->len -= gap->size - m_seq->len;
			ass_seq->seq[ass_seq->len] = '\0';
			ass_eg->len = ass_seq->len;
		} else {
			if (m_eg->len > gap->size) { // If the new edge length is larger than the gap, remove the gap
				// spare some more space for it.
				new_len = ass_seq->len + (m_eg->len - gap->size);
				if (new_len > ass_seq->full_len) {
					ass_seq->len = new_len;
					kroundup32(new_len);
					ass_seq->seq = (ubyte_t*) realloc(ass_seq->seq,
							sizeof(ubyte_t) * new_len);
					ass_seq->full_len = new_len;
					ass_seq->seq[ass_seq->len] = '\0';
				}
				memmove(&ass_seq->seq[gap->s_index + m_eg->len],
						&ass_seq->seq[gap->s_index + gap->size], right_seq->len);
				start = ass_eg->len - gap->s_index - gap->size;
				ass_eg->len = ass_seq->len; // Set them to be correct
			} else {
				if (left_ol)
					start = gap->s_index;
				if (right_ol)
					start = gap->s_index + (gap->size - m_eg->len);
				if (!left_ol && !right_ol)
					start = gap->s_index + (gap->size - m_eg->len) / 2;
			}
			memcpy(&ass_seq->seq[start], m_seq->seq, m_eg->len);
		}
	}
	bwa_free_read_seq(1, left_seq);
	bwa_free_read_seq(1, right_seq);
}

void clear_used_reads(edge *eg, const int reset_ctg_id) {
	int i = 0;
	bwa_seq_t *r;
	if (!eg || eg->alive)
		return;
	if (reset_ctg_id) {
		for (i = 0; i < eg->reads->len; i++) {
			r = g_ptr_array_index(eg->reads, i);
			r->used = 0;
			r->contig_id = -1;
		}
	}
	while (eg->reads->len)
		g_ptr_array_remove_index(eg->reads, 0);
}

/**
 * ori: if 0 means add all reads to left_eg; vice versa
 */
void combine_reads(edge *left_eg, edge *right_eg, const int upd_shift,
		const int gap, const int ori) {
	int i = 0;
	bwa_seq_t *r;
	if (upd_shift) {
		for (i = 0; i < right_eg->reads->len; i++) {
			r = g_ptr_array_index(right_eg->reads, i);
		}
	}
	if (ori) {
		for (i = 0; i < left_eg->reads->len; i++) {
			r = g_ptr_array_index(left_eg->reads, i);
			r->used = 0;
			r->contig_id = -1; // Reset this read to be not used.
			if (readarray_find(right_eg->reads, r) == NOT_FOUND) {
				readarray_add(right_eg, r); // This read will be marked as used by right_eg
			}
		}
		clear_used_reads(left_eg, 0); // Clear the reads from left_eg without reseting
	} else {
		for (i = 0; i < right_eg->reads->len; i++) {
			r = g_ptr_array_index(right_eg->reads, i);
			r->used = 0;
			r->contig_id = -1; // Reset this read to be not used.
			if (readarray_find(left_eg->reads, r) == NOT_FOUND) {
				readarray_add(left_eg, r);
			}
		}
		clear_used_reads(right_eg, 0);
	}
}

void concate_readarray(readarray *left_reads, readarray *right_reads) {
	int i = 0;
	bwa_seq_t *read = NULL;
	for (i = 0; i < right_reads->len; i++) {
		read = g_ptr_array_index(right_reads, i);
		g_ptr_array_add(left_reads, read);
	}
}

void merge_eg_to_left(edge *left_eg, edge *right_eg, const int gap) {
	combine_reads(left_eg, right_eg, 1, gap, 0);
	left_eg->contig = merge_seq_to_left(left_eg->contig, right_eg->contig, gap);
	left_eg->len = left_eg->contig->len;
	right_eg->alive = 0;
}

void merge_eg_to_right(edge *left_eg, edge *right_eg, const int gap) {
	combine_reads(left_eg, right_eg, 1, gap, 1);
	right_eg->contig = merge_seq_to_right(left_eg->contig, right_eg->contig,
			gap);
	right_eg->len = right_eg->contig->len;
	left_eg->alive = 0;
}

int get_mid_pos(readarray *ra, const int ori, const int lib_mean) {
	int i = 0, max_shift = 0, mid_shift = 0, len = 0;
	bwa_seq_t *read;

	if (!ra || ra->len == 0)
		return INVALID;
	g_ptr_array_sort(ra, (GCompareFunc) cmp_read_by_shift);
	//p_readarray(ra, 1);
	read = g_ptr_array_index(ra, ra->len - 1);
	len = read->shift + read->len;
	max_shift = read->shift;
	if (len < lib_mean) {
		mid_shift = max_shift / 2;
	} else {
		mid_shift = ori ? (lib_mean / 2) : (max_shift - (lib_mean / 2));
	}
	show_debug_msg(__func__, "max_shift = %d, mid_shift = %d \n", max_shift,
			mid_shift);
	for (i = 0; i < ra->len; i++) {
		read = g_ptr_array_index(ra, i);
		if (read->shift >= mid_shift)
			return i;
	}
	return INVALID;
}

void upd_reads(edge *eg, const int mismatches) {
	int i = 0, index = 0, overlap_len = 0;
	bwa_seq_t *read, *rev_read = NULL;
	char *pre_name = NULL;
	if (!eg || !eg->reads || eg->reads->len == 0)
		return;
	read = g_ptr_array_index(eg->reads, 0);
	if (!read || eg->len < read->len)
		return;
	for (i = 0; i < eg->reads->len; i++) {
		read = g_ptr_array_index(eg->reads, i);
		pre_name = read->name;
		if (read->shift < 0) {
			overlap_len = find_ol(read, eg->contig, mismatches);
			if (overlap_len > read->len / 4) {
				read->shift = read->len - overlap_len;
				continue;
			} else {
				overlap_len = find_ol(eg->contig, read, mismatches);
				if (overlap_len > read->len / 4) {
					read->shift = eg->len - overlap_len;
					continue;
				}
			}
		}
		rev_read = new_mem_rev_seq(read, read->len, 0);
		index = read->rev_com ? is_sub_seq(rev_read, 0, eg->contig, mismatches
				+ 2, 0) : is_sub_seq(read, 0, eg->contig, mismatches + 2, 0);
		bwa_free_read_seq(1, rev_read);
		if (index == NOT_FOUND) {
			if (g_ptr_array_remove_index_fast(eg->reads, i)) {
				read->used = 0; // This read is used somewhere and cannot be used again; to avoid infinite loop.
				read->contig_id = -1;
				read->shift = 0;
				i--;
			}
			//p_query("REMOVED", read);
			continue;
		}
		read->shift = index;
		read->contig_id = eg->id;
		read->used = 1;
	}
}

void p_edge_read(bwa_seq_t *p, FILE *log) {
	bwa_seq_t *p2;
	int i = 0, base = 30;
	char buf[BUFSIZ];
	if ((unsigned int) p->seq[0] <= 4) {
		p2 = new_seq(p, p->len, 0);
		map(p2);
		for (i = 0; i < base + p2->shift; i++)
			fputs(" ", log);

		if (p->rev_com) {
			sprintf(buf, "%s\t", p2->rseq);
			fputs(buf, log);
			sprintf(buf, "%d->%d\t%s\t%d\t[rev_com]", p2->contig_id, p2->shift,
					p2->name, p2->used);
			fputs(buf, log);
		} else {
			sprintf(buf, "%s\t", p2->seq);
			fputs(buf, log);
			sprintf(buf, "%d->%d\t%s\t%d", p2->contig_id, p2->shift, p2->name,
					p2->used);
			fputs(buf, log);
		}
		fputs("\n", log);
	}
}

gint cmp_seqs_by_shift(gpointer a, gpointer b) {
	bwa_seq_t* eg_a = *(bwa_seq_t**) a;
	bwa_seq_t* eg_b = *(bwa_seq_t**) b;
	return eg_a->shift - eg_b->shift;
}

void log_reads(edgearray *ea) {
	int i = 0, j = 0, k = 0;
	edge *eg_i;
	bwa_seq_t *r = 0;
	char line[BUFSIZ];
	FILE *reads_fp = xopen("read/reads_used.fa", "w");
	for (i = 0; i < ea->len; i++) {
		eg_i = g_ptr_array_index(ea, i);
		for (j = 0; j < eg_i->reads->len; j++) {
			r = g_ptr_array_index(eg_i->reads, j);
			sprintf(line, ">%s\n", r->name);
			fputs(line, reads_fp);
			for (k = 0; k < r->len; k++) {
				if (r->seq[k] > 4)
					fputc(r->seq[k], reads_fp);
				else
					fputc("acgtn"[(int) r->seq[k]], reads_fp);
			}
			fputc('\n', reads_fp);
		}
	}
	fclose(reads_fp);
}

void log_edge(const edge *eg) {
	FILE *log;
	char buf[BUFSIZ];
	int i = 0, j = 0;
	readarray *reads;
	bwa_seq_t *r;

	if (!eg || eg->id < 0)
		return;
	sprintf(buf, "log/%d.log", eg->id);
	log = xopen(buf, "w");
	reads = eg->reads;

	sprintf(buf, "Edge %d, length %d \n", eg->id, eg->len);
	fputs(buf, log);
	for (i = 0; i < 30; i++)
		fputs(" ", log);
	for (i = 0; i < eg->contig->len; i++) {
		if (eg->contig->seq[i] > 4)
			sprintf(buf, "%c", eg->contig->seq[i]);
		else
			sprintf(buf, "%c", "acgtn"[(int) eg->contig->seq[i]]);
		fputs(buf, log);
	}
	fputs(
			"\n================================================================\n",
			log);
	g_ptr_array_sort(reads, (GCompareFunc) cmp_seqs_by_shift);

	for (j = 0; j < reads->len; j++) {
		r = g_ptr_array_index(reads, j);
		if ((j + 1) % 30 == 0) {
			for (i = 0; i < 30; i++)
				fputs(" ", log);
			for (i = 0; i < eg->contig->len; i++) {
				if (eg->contig->seq[i] > 4)
					sprintf(buf, "%c", eg->contig->seq[i]);
				else
					sprintf(buf, "%c", "acgtn"[(int) eg->contig->seq[i]]);
				fputs(buf, log);
			}
			fputs("\n", log);
		}
		p_edge_read(r, log);
	}

	setbuf(log, NULL);
	fclose(log);
}

void readarray_add(edge *eg, bwa_seq_t *read) {
	if (!eg || !read)
		return;
	g_ptr_array_add(eg->reads, read);
	read->used = 1;
	read->contig_id = eg->id;
}

void readarray_uni_add(edge *eg, bwa_seq_t *read) {
	if (!eg || !read)
		return;
	g_ptr_array_uni_add(eg->reads, read);
	read->used = 1;
	read->contig_id = eg->id;
}

void readarray_remove(edge *eg, bwa_seq_t *read) {
	if (!eg || !read)
		return;
	g_ptr_array_remove(eg->reads, read);
	read->used = 0;
	read->contig_id = -1;
}

void readarray_remove_fast(edge *eg, const int index) {
	bwa_seq_t *read = NULL;
	if (!eg || index < 0 || index >= eg->reads->len)
		return;
	read = g_ptr_array_index(eg->reads, index);
	read->used = 0;
	read->contig_id = -1;
	g_ptr_array_remove_index_fast(eg->reads, index);
}
