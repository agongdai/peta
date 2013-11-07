/*
 * seq_util.c
 *
 *  Created on: May 16, 2011
 *      Author: Carl
 */

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <glib.h>
#include "utils.h"
#include "peseq.h"
#include "bwtaln.h"
#include "pechar.h"

extern unsigned char nst_nt4_table[256];

gint cmp_reads_by_name(gpointer a, gpointer b) {
	bwa_seq_t *seq_a = *((bwa_seq_t**) a);
	bwa_seq_t *seq_b = *((bwa_seq_t**) b);
	return (atoi(seq_a->name) - atoi(seq_b->name));
}

gint cmp_reads_by_cursor(gpointer a, gpointer b) {
	bwa_seq_t *seq_a = *((bwa_seq_t**) a);
	bwa_seq_t *seq_b = *((bwa_seq_t**) b);
	return seq_a->cursor - seq_b->cursor;
}

gint cmp_reads_by_rev_cursor(gpointer a, gpointer b) {
	bwa_seq_t *seq_a = *((bwa_seq_t**) a);
	bwa_seq_t *seq_b = *((bwa_seq_t**) b);
	return seq_b->cursor - seq_a->cursor;
}

gint cmp_reads_by_contig_id(gpointer a, gpointer b) {
	bwa_seq_t *seq_a = *((bwa_seq_t**) a);
	bwa_seq_t *seq_b = *((bwa_seq_t**) b);
	return (seq_a->contig_id - seq_b->contig_id);
}

gint cmp_reads_by_contig_locus(gpointer a, gpointer b) {
	bwa_seq_t *seq_a = *((bwa_seq_t**) a);
	bwa_seq_t *seq_b = *((bwa_seq_t**) b);
	return (seq_a->contig_locus - seq_b->contig_locus);
}

/**
 * Remove duplicates in an array
 */
GPtrArray *rm_duplicates(GPtrArray *reads) {
	int i = 0;
	bwa_seq_t *pre = NULL, *post = NULL;
	GPtrArray *uni_reads = NULL;
	if (!reads || reads->len < 2) {
		return reads;
	}
	uni_reads = g_ptr_array_sized_new(reads->len);
	g_ptr_array_sort(reads, (GCompareFunc) cmp_reads_by_name);
	pre = (bwa_seq_t*) g_ptr_array_index(reads, 0);
	g_ptr_array_add(uni_reads, pre);
	for (i = 1; i < reads->len; i++) {
		post = (bwa_seq_t*) g_ptr_array_index(reads, i);
		if (post != pre) {
			g_ptr_array_add(uni_reads, post);
			pre = post;
		}
	}
	g_ptr_array_free(reads, TRUE);
	return uni_reads;
}

/**
 * Remove reads with the same contig id and contig locus;
 * To avoid duplcate trial in connect_by_full_length
 */
GPtrArray *rm_dup_connectors(GPtrArray *reads) {
	int i = 0, *locus = NULL;
	bwa_seq_t *pre = NULL, *post = NULL;
	GPtrArray *uni_reads = NULL;
	if (!reads || reads->len < 2) {
		return reads;
	}
	uni_reads = g_ptr_array_sized_new(reads->len);
	g_ptr_array_sort(reads, (GCompareFunc) cmp_reads_by_contig_locus);
	pre = (bwa_seq_t*) g_ptr_array_index(reads, 0);
	g_ptr_array_add(uni_reads, pre);
	for (i = 1; i < reads->len; i++) {
		post = (bwa_seq_t*) g_ptr_array_index(reads, i);
		if (post->contig_id == pre->contig_id && post->contig_locus
				== pre->contig_locus) {
		} else {
			g_ptr_array_add(uni_reads, post);
			pre = post;
		}
	}
	g_ptr_array_free(reads, TRUE);
	return uni_reads;
}

/**
 * Find common reads of two lists of reads
 */
GPtrArray *interset_reads(GPtrArray *list_1, GPtrArray *list_2, GPtrArray *set) {
	uint32_t index_1 = 0, index_2 = 0;
	bwa_seq_t *read_1 = NULL, *read_2 = NULL;
	g_ptr_array_sort(list_1, (GCompareFunc) cmp_reads_by_name);
	g_ptr_array_sort(list_2, (GCompareFunc) cmp_reads_by_name);
	if (!set)
		set = g_ptr_array_sized_new(list_1->len);
	while (index_1 < list_1->len && index_2 < list_2->len) {
		read_1 = (bwa_seq_t*) g_ptr_array_index(list_1, index_1);
		//p_query("READ 1", read_1);
		while (index_2 < list_2->len) {
			read_2 = (bwa_seq_t*) g_ptr_array_index(list_2, index_2);
			//p_query("READ 2", read_2);
			if (read_1 == read_2) {
				g_ptr_array_add(set, read_1);
				//show_debug_msg(__func__, "INTERSET\n");
				index_2++;
				index_1++;
				break;
			}
			if (atoi(read_1->name) < atoi(read_2->name)) {
				read_1->pos = IMPOSSIBLE_NEGATIVE;
				index_1++;
				break;
			} else {
				read_1->pos = IMPOSSIBLE_NEGATIVE;
				index_2++;
			}
		}
	}
	return set;
}

/**
 * Save the query into disk.
 */
void save_fq(const bwa_seq_t *seqs, const char *fp_fn, const uint16_t ol) {
	char *header = malloc(FNLEN);
	FILE *fq_fp = xopen(fp_fn, "w+");
	int i;

	sprintf(header, "@%s\n", "0.1");
	fputs(header, fq_fp);
	for (i = seqs->len - ol; i < seqs->len; i++) {
		fputc("ACGTN"[(int) seqs->seq[i]], fq_fp);
	}
	fputs("\n+\n", fq_fp);
	for (i = ol - 1; i >= 0; i--) {
		fputc(':', fq_fp);
	}
	fflush(fq_fp);
	fclose(fq_fp);
	free(header);
}

void p_seq(const char *header, const ubyte_t *seq, const int len) {
	int i = 0;
	printf("[%s] ", header);
	for (i = 0; i < len; i++) {
		if (seq[i] > 4)
			printf("%c", seq[i]);
		else
			printf("%c", "ACGTN"[(int) seq[i]]);
	}
	printf("\n");
}

int trun_seq(bwa_seq_t *s, const int shift) {
	ubyte_t *seq = 0;
	if (!s || shift <= 0 || shift > s->len)
		return 0;
	if (shift > s->len) {
		s->len = 0;
		return 0;
	}
	if (shift == s->len) {
		s->len = 0;
		s->seq = 0;
		return 1;
	}
	seq = s->seq;
	memmove(s->seq, &seq[shift], s->len - shift);
	s->len -= shift;
	s->seq[s->len] = '\0';
	return 1;
}

bwa_seq_t *merge_seq_to_right(bwa_seq_t *s1, bwa_seq_t *s2, const int gap) {
	int i = 0;
	ubyte_t *rev = NULL;
	if (!s1 || !s2)
		return 0;
	if (!s1->seq || !s2->seq)
		return s2;
	s2->full_len = (s1->len + s2->len + 1 + gap);
	kroundup32(s2->full_len);
	s2->seq = (ubyte_t*) realloc(s2->seq, sizeof(ubyte_t) * s2->full_len);
	memmove(&s2->seq[s1->len + gap], s2->seq, sizeof(ubyte_t) * s2->len);
	for (i = 0; i < gap; i++) {
		s2->seq[s1->len + i] = N_NUCLITODE;
	}
	memcpy(s2->seq, s1->seq, sizeof(ubyte_t) * s1->len);
	s2->len += s1->len + gap;
	s2->seq[s2->len] = '\0';

	rev = s2->rseq;
	free(rev);
	rev = (ubyte_t*) malloc(s2->len + 1);
	memcpy(rev, s2->seq, s2->len);
	seq_reverse(s2->len, rev, 1);
	rev[s2->len] = '\0';
	s2->rseq = rev;
	return s2;
}

bwa_seq_t *merge_seq_to_left(bwa_seq_t *s2, bwa_seq_t *s1, const int gap) {
	int i = 0;
	ubyte_t *rev = NULL;
	if (!s1 || !s2)
		return 0;
	s2->full_len = (s1->len + s2->len + 1 + gap);
	kroundup32(s2->full_len);
	s2->seq = (ubyte_t*) realloc(s2->seq, sizeof(ubyte_t) * s2->full_len);
	for (i = 0; i < gap; i++) {
		s2->seq[s2->len + i] = N_NUCLITODE;
	}
	memcpy(&s2->seq[s2->len + gap], s1->seq, sizeof(ubyte_t) * s1->len);
	s2->len += s1->len + gap;
	s2->seq[s2->len] = '\0';

	rev = s2->rseq;
	free(rev);
	rev = (ubyte_t*) malloc(s2->len + 1);
	memcpy(rev, s2->seq, s2->len);
	seq_reverse(s2->len, rev, 1);
	rev[s2->len] = '\0';
	s2->rseq = rev;
	return s2;
}

/**
 * Merge partial s2 to s1.
 * For example, s2->len = 100, shift = 20, concat s2->seq[20:100] to s1->seq
 */
bwa_seq_t *merge_seq(bwa_seq_t *s1, bwa_seq_t *s2, const int shift) {
	if (!s1 || !s2)
		return 0;
	s1->full_len = (s1->len + s2->len + 1 - shift);
	kroundup32(s1->full_len);
	s1->seq = (ubyte_t*) realloc(s1->seq, sizeof(ubyte_t) * s1->full_len);
	memcpy(&s1->seq[s1->len], &s2->seq[shift], sizeof(ubyte_t) * (s2->len
			- shift));
	s1->len += s2->len - shift;
	s1->seq[s1->len] = '\0';
	return s1;
}

void map(bwa_seq_t *bwa_seq) {
	unsigned int i = 0;
	for (i = 0; i < bwa_seq->len; i++) {
		bwa_seq->seq[i] = "ACGTN"[(int) bwa_seq->seq[i]];
	}
	for (i = 0; i < bwa_seq->len; i++) {
		bwa_seq->rseq[i] = "ACGTN"[(int) bwa_seq->rseq[i]];
	}
}

int is_left_mate(const char *seq_id) {
	int id = atoi(seq_id);
	return id % 2 - 1;
}

int is_right_mate(const char *seq_id) {
	int id = atoi(seq_id);
	return id % 2;
}

int is_mates(const char *left, const char *right) {
	char *mate;
	int good = 0;
	if (!left || !right)
		return 0;
	mate = get_mate_name(left);
	if (strcmp(right, mate) == 0)
		good = 1;
	else
		good = 0;
	free(mate);
	return good;
}

char *get_mate_name(const char *seq_id) {
	char *mate = strdup(seq_id);
	if (is_left_mate(seq_id))
		sprintf(mate, "%d", atoi(seq_id) + 1);
	else
		sprintf(mate, "%d", atoi(seq_id) - 1);
	return mate;
}

index64 get_index(const char *seq_id) {
	return atoll(seq_id);
}

index64 get_mate_index(const index64 seq_id) {
	assert(seq_id >= 0);
	if (seq_id % 2 == 0)
		return seq_id + 1;
	else
		return seq_id - 1;
}

bwa_seq_t *get_mate(const bwa_seq_t *s, bwa_seq_t *seqs) {
	index64 mate_id = get_mate_index(get_index(s->name));
	return &seqs[mate_id];
}

bwa_seq_t *get_right_mate(const bwa_seq_t *left, bwa_seq_t *seqs) {
	index64 mate;
	if (is_left_mate(left->name)) {
		mate = get_mate_index(atoll(left->name));
		return &seqs[mate];
	} else
		return 0;
}

bwa_seq_t *get_left_mate(const bwa_seq_t *right, bwa_seq_t *seqs) {
	index64 mate;
	if (!is_left_mate(right->name)) {
		mate = get_mate_index(atoll(right->name));
		return &seqs[mate];
	} else
		return 0;
}

void p_shift_query(const bwa_seq_t *q, const int locus) {
	int i = 0;
	if (!q) {
		printf("Query is NULL \n");
		return;
	}
	for (i = 0; i < q->len - locus; i++) {
		printf(" ");
	}
	for (i = 0; i < q->len; i++) {
		if (q->rev_com) {
			if (q->rseq[i] > 4)
				printf("%c", q->rseq[i]);
			else
				printf("%c", "ACGTN"[(int) q->rseq[i]]);
		} else {
			if (q->seq[i] > 4)
				printf("%c", q->seq[i]);
			else
				printf("%c", "ACGTN"[(int) q->seq[i]]);
		}
	}
	for (i = 0; i < locus + 4; i++) {
		printf(" ");
	}
	if (q->rev_com)
		printf(" %s[reverse@%d]", q->name, q->pos);
	else
		printf(" %s[forward@%d]", q->name, q->pos);
	printf(" [%d: %d]", q->contig_id, q->contig_locus);
	printf("\n");
}

void p_query(const char *header, const bwa_seq_t *q) {
	unsigned int i = 0;
	if (!q) {
		printf("[%s] Query is NULL \n", header);
		return;
	}
	printf("[%s] [%s, %d]: \t", header, q->name, q->len);
	for (i = 0; i < q->len; i++) {
		if (q->rev_com) {
			if (q->rseq[i] > 4)
				printf("%c", q->rseq[i]);
			else
				printf("%c", "ACGTN"[(int) q->rseq[i]]);
		} else {
			if (q->seq[i] > 4)
				printf("%c", q->seq[i]);
			else
				printf("%c", "ACGTN"[(int) q->seq[i]]);
		}
	}
	printf(" [status: %d] ", q->status);
	if (q->rev_com)
		printf(" [    <<<<@%d]", q->pos);
	else
		printf(" [>>>>    @%d]", q->pos);
	printf(" [%d: %d]", q->contig_id, q->contig_locus);
	printf(" [cursor: %d]", q->cursor);
	//	printf("\n[rev_com] ");
	//	for (i = 0; i < q->len; i++) {
	//		if (q->rseq[i] > 4)
	//			printf("%c", q->rseq[i]);
	//		else
	//			printf("%c", "acgtn"[(int) q->rseq[i]]);
	//	}
	printf("\n");
}

void p_readarray(const GPtrArray *ra, const int all) {
	int i = 0;
	bwa_seq_t *p;
	show_debug_msg(__func__, "--------------------------------- \n");
	show_debug_msg(__func__, " # of Reads: %d \n", ra->len);
	for (i = 0; i < ra->len; i++) {
		p = g_ptr_array_index(ra, i);
		if ((!all && i % 200 == 0) || all)
			p_query(__func__, p);
		// show_debug_msg(__func__, "%d: %d_%d Read %s: %d \n", i,
		// p->contig_id, p->shift, p->name, p->shift);
	}
	show_debug_msg(__func__, "--------------------------------- \n");
}

void p_ctg_seq(const char *header, const bwa_seq_t *q) {
	unsigned int i = 0;
	if (!q) {
		printf("[%s] Contig seq is NULL \n", header);
		return;
	}
	printf("[%s] [%s, %d]: \t", header, q->name, q->len);
	for (i = 0; i < q->len; i++) {
		if (q->seq[i] > 4)
			printf("%c", q->seq[i]);
		else
			printf("%c", "ACGTN"[(int) q->seq[i]]);
	}
	printf("\n");
}

void ext_con(bwa_seq_t *contig, const ubyte_t c, const int ori) {
	if (contig->full_len <= contig->len + 2) {
		contig->full_len = contig->len + 2;
		kroundup32(contig->full_len);
		contig->seq = (ubyte_t*) realloc(contig->seq, sizeof(ubyte_t)
				* contig->full_len);
	}
	if (ori) {
		memmove(&contig->seq[1], contig->seq, contig->len);
		contig->seq[0] = c;
	} else {
		contig->seq[contig->len] = c;
	}
	contig->len++;
}

void ext_que(bwa_seq_t *q, const ubyte_t c, const int ori) {
	ubyte_t comp_c = (c >= 4) ? c : 3 - c;
	if (!ori) { // Left mate of the PET
		memmove(q->seq, &q->seq[1], q->len - 1);
		q->seq[q->len - 1] = c;
		memmove(&q->rseq[1], q->rseq, q->len - 1);
		q->rseq[0] = comp_c;
	} else {
		memmove(&q->seq[1], q->seq, q->len - 1);
		q->seq[0] = c;
		memmove(q->rseq, &q->rseq[1], q->len - 1);
		q->rseq[q->len - 1] = comp_c;
	}
}

int has_n(const bwa_seq_t *read, int max) {
	int i = 0;
	for (i = 0; i < read->len; i++) {
		if (read->seq[i] == 4 || read->seq[i] == 'n' || read->seq[i] == 'N')
			max--;
		if (max <= 0)
			return 1;
	}
	return 0;
}

void clear_array(GPtrArray *arr) {
	while (arr->len > 0)
		g_ptr_array_remove_index_fast(arr, 0);
}

/**
 * Create a 'sudo-read' for overlapping alignment
 */
bwa_seq_t *new_seq(const bwa_seq_t *query, const int ol, const int shift) {
	if (ol < 0 || ol + shift > query->len)
		return blank_seq(0);
	bwa_seq_t *p = (bwa_seq_t*) malloc(sizeof(bwa_seq_t));
	p->status = query->status;
	p->contig_id = query->contig_id;
	p->contig_locus = query->contig_locus;
	p->full_len = p->len = ol;
	p->rev_com = query->rev_com;
	p->pos = query->pos;

	if (query->name)
		p->name = strdup(query->name);
	else
		p->name = (char*) calloc(1, sizeof(char));
	p->seq = (ubyte_t*) malloc(ol);
	p->rseq = NULL;
	memcpy(p->seq, query->seq + shift, ol);
	p->rseq = (ubyte_t*) malloc(ol);
	memcpy(p->rseq, p->seq, p->len);
	seq_reverse(p->len, p->rseq, 1);
	p->cursor = query->cursor;

	return p;
}

void copy_partial(bwa_seq_t *s, bwa_seq_t *copied, int start, int len) {
	if (start + len > s->len)
		return;
	memcpy(copied->seq, s->seq + start, sizeof(char) * len);
	memcpy(copied->rseq, s->rseq + (s->len - start - len), sizeof(char) * len);
	copied->len = len;
}

/**
 * Noted: to avoid segfault when calling realloc, 
 * 			here use memcpy, do not simply swap the address
 **/
void switch_fr(bwa_seq_t *s) {
	ubyte_t *b = (ubyte_t*) malloc(sizeof(ubyte_t) * s->len);
	memcpy(b, s->seq, sizeof(ubyte_t) * s->len);
	memcpy(s->seq, s->rseq, sizeof(ubyte_t) * s->len);
	memcpy(s->rseq, b, sizeof(ubyte_t) * s->len);
	free(b);
}

void set_rev_com(bwa_seq_t *s) {
	if (!s)
		return;
	if (s->rseq)
		free(s->rseq);
	s->rseq = (ubyte_t*) malloc(s->len + 1);
	memcpy(s->rseq, s->seq, s->len);
	seq_reverse(s->len, s->rseq, 1);
	s->rseq[s->len] = '\0';
}

void reset_read(bwa_seq_t *r) {
	r->pos = IMPOSSIBLE_NEGATIVE;
	r->cursor = -1;
	r->contig_id = -1;
	r->contig_locus = -1;
	r->rev_com = 0;
}

void reset_to_fresh(bwa_seq_t *r) {
	r->status = FRESH;
	reset_read(r);
}

void reset_to_dead(bwa_seq_t *r) {
	r->status = DEAD;
	reset_read(r);
}

void reset_to_hang(bwa_seq_t *r) {
	r->status = HANG;
	reset_read(r);
}

bwa_seq_t *blank_seq(const int len) {
	bwa_seq_t *p = (bwa_seq_t*) malloc(sizeof(bwa_seq_t));
	p->status = 0;
	p->contig_id = 0;
	p->full_len = len;
	p->len = 0;
	p->rev_com = 0;

	p->name = NULL;
	p->seq = (ubyte_t*) calloc(len, sizeof(ubyte_t));
	p->rseq = (ubyte_t*) calloc(len, sizeof(ubyte_t));

	return p;
}

bwa_seq_t *new_mem_rev_seq(const bwa_seq_t *query, const int ol,
		const int shift) {
	bwa_seq_t *p = NULL;
	ubyte_t *tmp = NULL;
	p = new_seq(query, ol, shift);
	tmp = p->seq;
	p->seq = p->rseq;
	p->rseq = tmp;
	return p;
}

bwa_seq_t *new_rev_seq(const bwa_seq_t *query) {
	bwa_seq_t *p = (bwa_seq_t*) malloc(sizeof(bwa_seq_t));
	ubyte_t *tmp = NULL;

	p->status = FRESH;
	p->contig_id = query->contig_id;
	p->full_len = p->len = query->len;
	p->name = query->name;
	p->rev_com = query->rev_com;

	tmp = p->seq;
	p->seq = query->rseq;
	p->rseq = tmp;
	return p;
}

void save_con(const char *header, const bwa_seq_t *contig, FILE *tx_fp) {
	int i = 0;
	char c = to_upper_lower('N');
	if (!contig)
		return;
	fputs(header, tx_fp);

	for (i = 0; i < contig->len; i++) {
		c = contig->seq[i];
		c = "ACGTN"[(int) c];
		fputc(c, tx_fp);
		if ((i + 1) % LINELEN == 0) {
			fputc('\n', tx_fp);
		}
	}
	if (i % LINELEN != 0)
		fputc('\n', tx_fp);
}

int same_q(const bwa_seq_t *query, const bwa_seq_t *seq) {
	if (query->len != seq->len)
		return 0;
	return seq_ol(query, seq, query->len, 0) == -1 ? 0 : 1;
}

/**
 * Check whether too many N's in the sequence
 */
int too_many_ns(const ubyte_t *s, const int k) {
	int n_n = 0, i = 0;
	for (i = 0; i < k; i++) {
		if (s[i] == 4)
			n_n++;
	}
	if (n_n * 2 >= k)
		return 1;
	else
		return 0;
}

/**
 * Never use the low quality reads: with >= 3 N's in reads; repetitive reads
 */
void mark_low_qua_reads(bwa_seq_t *seqs, index64 n_seqs) {
	index64 i = 0, n_dead = 0, j = 0;
	int max_ns = 0;
	bwa_seq_t *r = NULL, *part = NULL;
	if (!seqs || n_seqs <= 0)
		return;
	max_ns = seqs->len / 4;
	for (i = 0; i < n_seqs; i++) {
		r = &seqs[i];
		if (has_n(r, max_ns) || is_biased_q(r) || is_repetitive_q(r)) {
			r->status = DEAD;
			//p_query(__func__, r);
			n_dead++;
		}
	}
	show_msg(__func__, "N_DEADS: %d\n", n_dead);
}

/**
 * Check whether all bases are the same, allow one base different
 */
int same_bytes(const ubyte_t *s, const int k) {
	int *c = (int*) calloc(5, sizeof(int));
	int i = 0, is_same = 0;
	int thre = k - 1;
	for (i = 0; i < k; i++) {
		c[s[i]]++;
	}
	if (c[4] >= thre)
		is_same = 1;
	for (i = 0; i < 4; i++) {
		if (c[i] + c[4] >= thre) {
			is_same = 1;
			break;
		}
	}
	free(c);
	return is_same;
}

/**
 * Smith-Waterman algorithm to compare two sequences
 * Score for match: 1; mismatch: -1; indel: -2
 * @min_score: the minimal score we accept
 */
int smith_waterman_simple(const bwa_seq_t *seq_1, const bwa_seq_t *seq_2,
		int *seq_1_start, int *seq_1_stop, int *seq_2_start, int *seq_2_stop,
		int min_score) {
	int score_m = 1, score_mis = 2, score_gap = 3;
	int *scores = NULL;
	int columns = seq_1->len + 1, max_score = 0, rows = seq_2->len + 1;
	int i = 0, j = 0, s = 0;
	int up_left = 0, up = 0, left = 0, remain = 0;
	// A simulated 2d-array to store scores.
	scores = (int*) calloc(rows * columns, sizeof(int));
	*seq_1_stop = 0;
	*seq_2_stop = 0;
	for (i = 1; i < rows; i++) {
		for (j = 1; j < columns; j++) {
			up = scores[(i - 1) * columns + j] - score_gap;
			left = scores[i * columns + j - 1] - score_gap;
			if (seq_1->seq[j - 1] == seq_2->seq[i - 1])
				up_left = scores[(i - 1) * columns + j - 1] + score_m;
			else
				up_left = scores[(i - 1) * columns + j - 1] - score_mis;
			scores[i * columns + j] = max3(left, up, up_left);
		}
		for (j = 0; j < columns; j++) {
			if (scores[i * columns + j] > max_score) {
				max_score = scores[i * columns + j];
				*seq_1_stop = j;
				*seq_2_stop = i;
				remain = min(rows - i, columns - j);
				if (max_score + remain < min_score) {
					free(scores);
					return 0;
				}
			}
		}
	}
//	printf("\t\t");
//	for (i = 0; i < seq_1->len; i++) {
//		printf("%c\t", "ACGTN"[seq_1->seq[i]]);
//	}
//	printf("\n");
//	for (i = 0; i < rows; i++) {
//		if (i == 0)
//			printf("\t");
//		else
//			printf("%c\t", "ACGTN"[seq_2->seq[i - 1]]);
//		for (j = 0; j < columns; j++) {
//			printf("%d\t", scores[i * columns + j]);
//		}
//		printf("\n");
//	}
	// Back trace to find the starting points
	*seq_1_start = *seq_1_stop - 1;
	*seq_2_start = 0;
	for (i = *seq_2_stop; i >= 1; i--) {
		s = scores[(i) * columns + *seq_1_start + 1];
		up_left = scores[(i - 1) * columns + *seq_1_start];
		left = scores[(i) * columns + *seq_1_start];
//		printf("%d\n", up_left);
//		printf("%d\t%d\n", left, s);
		if (s == 1) {
			*seq_2_start = i - 1;
			break;
		}
//		show_debug_msg(__func__, "seq_1: %c; seq_2: %c \n",
//				"ACGTN"[seq_1->seq[*seq_1_start]], "ACGTN"[seq_2->seq[i - 1]]);
		if (seq_1->seq[*seq_1_start] == seq_2->seq[i - 1]) {
			if (s == up_left + score_m)
				*seq_1_start = *seq_1_start - 1;
		} else {
			if (s == up_left - score_mis || s == left - score_gap)
				*seq_1_start = *seq_1_start - 1;
		}
//		printf("i: %d; seq_1_start: %d; seq_2_start: %d\n", i, *seq_1_start, i
//				- 1);
	}
//	printf("seq_1: [%d, %d]; seq_2: [%d, %d] \n", *seq_1_start, *seq_1_stop,
//			*seq_2_start, *seq_2_stop);
	free(scores);
	return max_score;
}

/**
 * http://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm
 */
int smith_waterman(const bwa_seq_t *seq_1, const bwa_seq_t *seq_2,
		const int score_mat, const int score_mis, const int score_gap,
		const int min_acceptable_score) {
	int *previous_row = NULL, *current_row = NULL;
	int columns = seq_1->len + 1, max_score = 0, rows = seq_2->len + 1;
	int i = 0, j = 0, max = 0;
	int up_left = 0, up = 0, left = 0;
	previous_row = (int*) calloc(columns + 1, sizeof(int));
	current_row = (int*) calloc(columns + 1, sizeof(int));
	for (i = 1; i < rows; i++) {
		for (j = 1; j < columns; j++) {
			up = previous_row[j] + score_gap; // not updated, so is the value from previous row
			left = current_row[j - 1] + score_gap; // updated already, so is the value from current row
			if (seq_1->seq[j - 1] == seq_2->seq[i - 1])
				up_left = previous_row[j - 1] + score_mat;
			else
				up_left = previous_row[j - 1] + score_mis;
			max = up > left ? up : left;
			max = up_left > max ? up_left : max;
			current_row[j] = max;
		}
		//		printf("Previous row: \n");
		for (j = 0; j < columns; j++) {
			//			printf("%d,", previous_row[j]);
			previous_row[j] = current_row[j];
			max_score = current_row[j] > max_score ? current_row[j] : max_score;
		}
		//		printf("\n");
		//		printf("Row %d: \t", i);
		for (j = 0; j < columns; j++) {
			if (current_row[j] == max_score) {
				//				printf("[%d]\t", current_row[j]);
			} else {
				//				printf("%d\t", current_row[j]);
			}
		}
		//		printf("\t\t");
		//		printf("Max score: %d \n", max_score);
		// If the minimal acceptable score is not reachable, stop and return.
		if ((max_score + (rows - i) * score_mat) <= min_acceptable_score) {
			free(previous_row);
			free(current_row);
			return -1;
		}
	}
	free(previous_row);
	free(current_row);
	return max_score;
}

/**
 * Smith-waterman local alignment algorithm.
 */
int similar_seqs(const bwa_seq_t *query, const bwa_seq_t *seq,
		const int mismatches, const int max_n_gaps, const int score_mat,
		const int score_mis, const int score_gap) {
	int min_acceptable_score = 0, min_len = 0, similarity_score = 0;
	if (!query || !seq || !seq->seq || !query->seq || mismatches < 0)
		return 0;
	if (get_abs(query->len - seq->len) > mismatches)
		return 0;
	min_len = query->len;
	min_len = min_len > seq->len ? seq->len : min_len;
	min_acceptable_score = min_len * score_mat + mismatches * score_mis
			+ get_abs(query->len - seq->len) * score_gap;
	similarity_score = smith_waterman(query, seq, score_mat, score_mis,
			score_gap, min_acceptable_score);
	//show_debug_msg(__func__, "Score: %d\n", similarity_score);
	//show_debug_msg(__func__, "min_acceptable_score: %d\n", min_acceptable_score);
	if (similarity_score >= min_acceptable_score)
		return similarity_score;
	return 0;
}

/**
 * Shift is that position "starting of query on the read"
 * Offset is the position on the read, from which the k-mer is matched.
 *
 * Query: ---------()-------------------------- 35
 * Read:      -----()------------------------------ 35
 * Shift: -4, offset: 5
 * Check: whether from the position of () to the end of query,
 * 		they are different within specified mismatches;
 * 		starting position for query is offset-shift, for read is offset
 *
 * Query:     -----()---------------------- 27
 * Read:  ---------()-------------------------- 35
 * Shift: 4, offset: 9
 * Check: whether from beginning of query to the end of query,
 * 		they are different within specified mismatches;
 * 		starting position for query is offset-shift, for read is offset
 *
 * If query->len + shift >= read->len, meaning that no further base can be
 * 		obtained from the read, just return NOT_FOUND.
 *
 * RETURN: if not found, return NOT_FOUND; otherwise, return how many mismatches remained.
 */

int is_sub_seq_aln(const ubyte_t *query, const int q_len, const int shift,
		const int offset, const bwa_seq_t *seq, int mismatches, const int ol) {
	unsigned int i = 0, start = 0, end = q_len - 1;
	if (!query || !seq || !seq->seq)
		return NOT_FOUND;
	if (q_len + shift > seq->len)
		end = seq->len - shift - 1;
	start = (shift >= 0) ? 0 : get_abs(shift);
	if (ol) {
		end = (end > (start + ol - 1)) ? end : (start + ol - 1);
	}
	// start = offset - shift;
	for (i = start; i <= end; i++) {
		if (query[i] != seq->seq[i + shift]) {
			mismatches--;
		}
		if (mismatches < 0)
			return NOT_FOUND;
	}
	return mismatches;
}

/**
 * To check whether a seq can be found on another seq with certain mismatches.
 * query: ----()------------------------------- 35
 * seq:   ------------------------------------------------------- 50
 *
 * shift = 0; ol = 0; query->len == seq->len: check whether two seqs are identical;
 * shift = 0; ol = 0: check whether query is a subsequence of seq.
 *
 * shift = 4, ol_len = 23: means that check whether the subsequence of query [4, 27) is on seq.
 * If ol_len = 0: the subsequence to compare is from shift to the end.
 */
int is_sub_seq(const bwa_seq_t *query, const int shift, const bwa_seq_t *seq,
		int mismatches, const int ol) {
	unsigned int i = 0, start = 0, offset = 0, end = query->len;
	int nm;
	if (!query || !seq || !seq->seq || !query->seq)
		return NOT_FOUND;
	start = (shift >= 0) ? shift : get_abs(shift);
	if (ol) {
		end = (end > (start + ol)) ? (start + ol) : end;
	}
	if ((end - start) > seq->len)
		return NOT_FOUND;
	for (offset = 0; offset <= seq->len - (end - start); offset++) {
		nm = mismatches;
		for (i = start; i < end; i++) {
			if (query->seq[i] != seq->seq[i + offset - shift]) {
				nm--;
			}
			if (nm < 0)
				break;
		}
		// if 2 mismatches allowed:
		// nm = 0: 2 mis; nm = 1: 1 mis; nm = 2: 0 mis.
		if (nm >= 0)
			return offset;
	}
	return NOT_FOUND;
}

/**
 *
 */
int is_sub_seq_byte(const ubyte_t *query, const int q_len, const int shift,
		const bwa_seq_t *seq, int mismatches, const int ol) {
	unsigned int i = 0, start = 0, offset = 0, end = q_len;
	int nm;
	if (!query || !seq || !seq->seq)
		return NOT_FOUND;
	start = (shift >= 0) ? shift : get_abs(shift);
	if (ol) {
		end = (end > (start + ol)) ? (start + ol) : end;
	}
	if ((end - start) > seq->len)
		return NOT_FOUND;
	for (offset = 0; offset <= seq->len - (end - start); offset++) {
		nm = mismatches;
		for (i = start; i < end; i++) {
			// If either base is 'N', treat them as the same
			//if (query[i] == 4 || seq->seq[i + offset - shift])
			//	continue;
			if (query[i] != seq->seq[i + offset - shift]) {
				nm--;
			}
			if (nm < 0)
				break;
		}
		// if 2 mismatches allowed:
		// nm = 0: 2 mis; nm = 1: 1 mis; nm = 2: 0 mis.
		if (nm >= 0)
			return offset;
	}
	return NOT_FOUND;
}

// Check whether two seqs share some portion whose length is ol
int share_subseq_byte(const ubyte_t *seq_1, const int len,
		const bwa_seq_t *seq_2, const int mismatches, const int ol) {
	int i = 0;
	if (!seq_1 || !seq_2)
		return 0;
	for (i = 0; i < len - ol; i++) {
		if ((ol < seq_2->len) && is_sub_seq_byte(seq_1, len, i, seq_2,
				mismatches, ol) != NOT_FOUND)
			return 1;
	}
	return 0;
}

/**
 * Check whether two sequences share similar regions (length 'ol') with max 'mismatches'.
 * If the a seq is rev_com, use its reverse complement to match
 */
int seq_ol(const bwa_seq_t *left_seq, const bwa_seq_t *right_seq, const int ol,
		int mismatches) {
	int i = 0;
	int n_mis = 0;
	ubyte_t *left = NULL, *right = NULL;
	if (!left_seq || !right_seq || ol <= 0 || ol > left_seq->len || ol
			> right_seq->len)
		return -1;
	left = left_seq->rev_com ? left_seq->rseq : left_seq->seq;
	right = right_seq->rev_com ? right_seq->rseq : right_seq->seq;
	for (i = 0; i < ol; i++) {
		// If there is an 'N' on either seq, just treat them as the same
		//if (left[i + (left_seq->len - ol)] == 4 || right[i] == 4)
		//	continue;
		if (left[i + (left_seq->len - ol)] != right[i]) {
			n_mis++;
		}
		if (n_mis > mismatches)
			return -1;
	}
	return n_mis;
}

/**
 * All 'N's would be ignored.
 */
int similar_bytes(ubyte_t *b_1, ubyte_t *b_2, int len, int mismatches) {
	int i = 0;
	for (i = 0; i < len; i++) {
		//if (b_1[i] == 4 || b_2[i] == 4)
		//	continue;
		if (b_1[i] != b_2[i])
			mismatches--;
		if (mismatches < 0)
			return 0;
	}
	return 1;
}

/**
 * If both the ref and query have similar heads and tails, return true.
 */
int head_tail_similar(bwa_seq_t *ref, bwa_seq_t *query, const int len,
		int mismatches, int *rev_com) {
	int similar = 0;
	if (!ref || !query || ref->len < len || query->len < len)
		return 0;
	// Check whether forward sequences are similar
	//if (query->full_len == 1000) {
	//	p_query("REF", ref);
	//	p_query("QUE", query);
	//}
	// Head: similar_bytes(ref->seq, query->seq, len, mismatches)
	if (similar_bytes(ref->seq, query->seq, len, mismatches) && similar_bytes(
			ref->seq + (ref->len - len), query->seq + (query->len - len), len,
			mismatches)) {
		similar = 1;
		*rev_com = 0;
	} else { // Check reverse sequence
		// Head: similar_bytes(ref->seq, query->rseq, len, mismatches)
		if (similar_bytes(ref->seq, query->rseq, len, mismatches)
				&& similar_bytes(ref->seq + (ref->len - len), query->rseq
						+ (query->len - len), len, mismatches)) {
			similar = 1;
			*rev_com = 1;
		}
	}
	//if (query->full_len == 1000)
	//	show_debug_msg(__func__, "SIMILAR: %d; REV_COM: %d\n", similar,
	//			*rev_com);
	return similar;
}

/**
 * Find the maximum overlap between the mate and template;
 * The overlap length is limited within [min_len, max_len] to save comparison time
 */
int find_ol_within_k(const bwa_seq_t *mate, const bwa_seq_t *tpl,
		const int mismatches, const int min_len, const int max_len,
		const int ori, int *n_mis) {
	int i = 0, olpped = -1;
	int m = 0;
	if (!mate || !tpl || min_len <= 0 || max_len <= 0 || max_len < min_len)
		return 0;
	m = min3(max_len, mate->len, tpl->len);
	for (i = m; i >= min_len; i--) {
		olpped = ori ? seq_ol(mate, tpl, i, mismatches) : seq_ol(tpl, mate, i,
				mismatches);
		if (olpped >= 0) {
			// How many mismatches in the overlap
			*n_mis = olpped;
			olpped = i;
			break;
		}
	}
	return olpped;
}

/**
 * Find the maximum overlapping length between two sequences
 */
int find_fr_ol_within_k(const bwa_seq_t *mate, const bwa_seq_t *tail,
		const int mismatches, const int min_len, const int max_len,
		const int ori, int *rev_com, int *n_mis) {
	int olpped = 0;
	//p_ctg_seq(__func__, mate);
	//p_ctg_seq(__func__, tail);
	olpped = find_ol_within_k(mate, tail, mismatches, min_len, max_len, ori,
			n_mis);
	if (olpped >= min_len && olpped <= max_len) {
		*rev_com = 0;
		return olpped;
	}
	// Switch the forward and reverse sequences temply
	// Make sure to switch back before returning.
	switch_fr(mate);
	//p_ctg_seq(__func__, mate);
	//p_ctg_seq(__func__, tail);
	olpped = find_ol_within_k(mate, tail, mismatches, min_len, max_len, ori,
			n_mis);
	if (olpped >= min_len && olpped <= max_len) {
		*rev_com = 1;
		switch_fr(mate);
		return olpped;
	}
	switch_fr(mate);
	return -1;
}

int find_ol(const bwa_seq_t *left_seq, const bwa_seq_t *right_seq,
		const int mismatches) {
	int i = 0, min_len = 0, olpped = 0;
	if (!left_seq || !right_seq)
		return 0;
	min_len = left_seq->len;
	min_len = (left_seq->len < right_seq->len) ? min_len : right_seq->len;
	for (i = mismatches; i < mismatches * 2; i++) {
		olpped = seq_ol(left_seq, right_seq, i, 1);
		if (olpped >= 0) {
			return i;
		}
	}
	for (i = min_len; i >= mismatches * 2; i--) {
		olpped = seq_ol(left_seq, right_seq, i, mismatches);
		if (olpped >= 0) {
			olpped = i;
			break;
		}
	}
	return olpped;
}

int is_repetitive_q(const bwa_seq_t *query) {
	int i = 0, is_rep = 0, rep_b = 0;
	ubyte_t *seq = 0, c = 0, c2 = 0;
	if (!query || !query->seq)
		return 0;
	seq = query->seq;

	for (rep_b = 1; rep_b <= NO_REPEAT_BASES; rep_b++) {
		is_rep = 1;
		for (i = 0; i < query->len - rep_b; i++) {
			c = seq[i];
			c2 = seq[i + rep_b];
			if (c != c2) {
				is_rep = 0;
				break;
			}
		}
		if (is_rep)
			return 1;
	}
	return 0;
}

int is_bad_query(bwa_seq_t *query) {
	return same_bytes(query->seq, query->len) || is_repetitive_q(query);
}

int is_biased_q(const bwa_seq_t *query) {
	int *c = (int*) calloc(5, sizeof(int));
	int i = 0, is_biased = 0;
	int thre = query->len * 4 / 5;
	for (i = 0; i < query->len; i++) {
		c[query->seq[i]]++;
	}
	for (i = 0; i < 5; i++) {
		if (c[i] >= thre) {
			is_biased = 1;
			break;
		}
	}
	free(c);
	return is_biased;
}

int has_rep_pattern(const bwa_seq_t *read) {
	int i = 0, j = 0, is_rep = 1, n_mis = 0;
	ubyte_t c = 0, c2 = 0;
	for (i = 0; i < read->len - NO_REPEAT_LEN - 1; i++) {
		is_rep = 1;
		n_mis = 0;
		for (j = i; j < i + NO_REPEAT_LEN; j++) {
			c = read->seq[j];
			c2 = read->seq[j + 1];
			if (c != c2) {
				n_mis++;
			}
			if (n_mis > 1) {
				is_rep = 0;
				break;
			}
		}
		if (is_rep)
			return 1;
	}
	return 0;
}

void pe_reverse_seqs(bwa_seq_t *seqs, const int n_seqs) {
	int i = 0;
	bwa_seq_t *s;
	for (i = 0; i < n_seqs; i++) {
		s = &seqs[i];
		seq_reverse(s->len, s->seq, 0);
		//seq_reverse(s->len, s->rseq, 0);
	}
}

/**
 * Return 1:
 * The mate of 'read' should be in the used area;
 * Return 0:
 * The mate of 'read' should be in the future
 */
int is_paired(const bwa_seq_t *read, const int ori) {
	int paired = 0;
	if (ori) {
		if ((is_left_mate(read->name) && !read->rev_com) || (is_right_mate(
				read->name) && read->rev_com)) {
			paired = 1;
		}
	} else {
		if ((is_right_mate(read->name) && !read->rev_com) || (is_left_mate(
				read->name) && read->rev_com)) {
			paired = 1;
		}
	}
	return paired;
}

/**
 * Count how many mismatches on the overlapping region
 * Query: ----------------------------------
 * Seq:                           ----*--*----------------------
 * Return 2.
 */
int get_mismatches_on_ol(const bwa_seq_t *query, const bwa_seq_t *seq,
		const int ol, const int max) {
	int n_mismatches = 0, i = 0;
	if (ol <= 0)
		return 0;
	for (i = 0; i < ol; i++) {
		if (n_mismatches >= max)
			return n_mismatches;
		if (query->seq[query->len - ol + i] != seq->seq[i])
			n_mismatches++;
	}
	return n_mismatches;
}

void free_read_seq(bwa_seq_t *p) {
	if (p) {
		free(p->name);
		free(p->seq);
		free(p->rseq);
		free(p);
	}
}
