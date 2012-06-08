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
	gpointer *data;
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

// param ori: 0 means ra_1 is at the left, and ra_2 is at the right
readarray *get_paired_reads(readarray *ra_1, readarray *ra_2, bwa_seq_t *seqs,
		const int ori) {
	readarray *paired = g_ptr_array_sized_new(INIT_N_READ_PAIRED);
	int i = 0, right_dir = 0;
	bwa_seq_t *read_1, *read_2;
	if (!ra_1 || !ra_2 || ra_1->len == 0 || ra_2->len == 0)
		return paired;
	for (i = 0; i < ra_1->len; i++) {
		read_1 = g_ptr_array_index(ra_1, i);
		if (ori) {
			if (is_left_mate(read_1->name))
				right_dir = 0;
			else
				right_dir = 1;
		} else {
			if (is_left_mate(read_1->name))
				right_dir = 1;
			else
				right_dir = 0;
		}
		if (right_dir) {
			read_2 = get_mate(read_1, seqs);
			if (readarray_find(ra_2, read_2) != NOT_FOUND) {
				g_ptr_array_add(paired, read_1);
			}
		}
	}
	return paired;
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
			r->shift += left_eg->len + gap;
		}
	}
	if (ori) {
		for (i = 0; i < left_eg->reads->len; i++) {
			r = g_ptr_array_index(left_eg->reads, i);
			if (readarray_find(right_eg->reads, r) == NOT_FOUND) {
				r->contig_id = right_eg->id;
				g_ptr_array_add(right_eg->reads, r);
			}
		}
		clear_used_reads(left_eg, 0);
	} else {
		for (i = 0; i < right_eg->reads->len; i++) {
			r = g_ptr_array_index(right_eg->reads, i);
			if (readarray_find(left_eg->reads, r) == NOT_FOUND) {
				r->contig_id = left_eg->id;
				g_ptr_array_add(left_eg->reads, r);
			}
		}
		clear_used_reads(right_eg, 0);
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

gint cmp_read_by_shift(gpointer a, gpointer b) {
	bwa_seq_t *read_a = *((bwa_seq_t**) a);
	bwa_seq_t *read_b = *((bwa_seq_t**) b);
	return (read_a->shift - read_b->shift);
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
	int i = 0, index = 0;
	bwa_seq_t *read;
	for (i = 0; i < eg->reads->len; i++) {
		read = g_ptr_array_index(eg->reads, i);
		if (read->shift < 0)
			continue;
		index = is_sub_seq(read, 0, eg->contig, mismatches, 0);
		if (index == NOT_FOUND) {
			index = is_sub_seq(read, 0, eg->contig, mismatches + 2, 0);
			if (index == NOT_FOUND) {
				g_ptr_array_remove_index_fast(eg->reads, i);
				read->used = 0; // This read is used somewhere and cannot be used again; to avoid infinite loop.
				read->contig_id = -1;
				read->shift = 0;
				i--;
				continue;
			}
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
