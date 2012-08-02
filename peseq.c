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
#include "utils.h"
#include "peseq.h"
#include "bwtaln.h"
#include "bntseq.h"
#include "pechar.h"
#include "pehash.h"

extern unsigned char nst_nt4_table[256];

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
		fputc("acgtn"[(int) seqs->seq[i]], fq_fp);
	}
	fputs("\n+\n", fq_fp);
	for (i = ol - 1; i >= 0; i--) {
		fputc(':', fq_fp);
	}
	fflush(fq_fp);
	fclose(fq_fp);
	free(header);
}

seq *read_seq(const char *fn) {
	unsigned int counter = 0, in_des = 0, h_c = 0;
	char c;
	seq *fa = (seq*) malloc(sizeof(seq));
	fa->h = (char*) malloc(BUFSIZE);
	fa->m = BUFSIZE * 2;
	fa->s = (char*) malloc(fa->m);
	FILE *fp = xopen(fn, "r");

	while ((c = (char) fgetc(fp))) {
		h_c = 0;
		if (c == -1)
			break;
		if (c == '>') {
			fa->h[h_c++] = c;
			in_des = 1;
			continue;
		}
		if (c == '\n') {
			in_des = 0;
			fa->h[h_c++] = c;
			fa->h[h_c++] = '\0';
			continue;
		}
		if (in_des) {
			fa->h[h_c++] = c;
			continue;
		}

		if (1 + fa->l >= fa->m) {
			fa->m = fa->l + 2;
			kroundup32(fa->m); // Adjust the size to nearest 2^k
			//			fprintf(stderr, "[fa->l, fa->m]: [%zd, %zd]\n", fa->l, fa->m);
			fa->s = (char*) realloc(fa->s, fa->m);
		}
		fa->s[fa->l++] = c;
		counter++;
	}
	fclose(fp);
	return fa;
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
	return s2;
}

bwa_seq_t *merge_seq_to_left(bwa_seq_t *s2, bwa_seq_t *s1, const int gap) {
	int i = 0;
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
	return s2;
}

void map(bwa_seq_t *bwa_seq) {
	unsigned int i = 0;
	for (i = 0; i < bwa_seq->len; i++) {
		bwa_seq->seq[i] = "acgtn"[(int) bwa_seq->seq[i]];
	}
	for (i = 0; i < bwa_seq->len; i++) {
		bwa_seq->rseq[i] = "acgtn"[(int) bwa_seq->rseq[i]];
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
	if (!left || !right)
		return 0;
	mate = get_mate_name(left);
	if (!strcmp(right, mate))
		return 1;
	else
		return 0;
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

void p_query(const char *header, const bwa_seq_t *q) {
	unsigned int i = 0;
	if (!q) {
		printf("[%s] Query is NULL \n", header);
		return;
	}
	printf("[%s] [%s, %d]: ", header, q->name, q->len);
	for (i = 0; i < q->len; i++) {
		if (q->seq[i] > 4)
			printf("%c", q->seq[i]);
		else
			printf("%c", "acgtn"[(int) q->seq[i]]);
	}
	if (q->is_in_c_pool)
		printf(" [pool]");
	else
		printf(" [no_pool]");
	if (q->rev_com)
		printf(" [rev_com]");
	else
		printf(" [not_rev_com]");
	if (q->used)
		printf(" [%d, %d]", q->contig_id, q->shift);
	else
		printf(" [not_used]");
	//	printf("\n[rev_com] ");
	//	for (i = 0; i < q->len; i++) {
	//		if (q->rseq[i] > 4)
	//			printf("%c", q->rseq[i]);
	//		else
	//			printf("%c", "acgtn"[(int) q->rseq[i]]);
	//	}
	printf("\n");
}

void p_ctg_seq(const char *header, const bwa_seq_t *q) {
	unsigned int i = 0;
	if (!q) {
		printf("[%s] Contig seq is NULL \n", header);
		return;
	}
	printf("[%s] [%s, %d]: ", header, q->name, q->len);
	for (i = 0; i < q->len; i++) {
		if (q->seq[i] > 4)
			printf("%c", q->seq[i]);
		else
			printf("%c", "acgtn"[(int) q->seq[i]]);
	}
	printf("\n");
}

void ext_con(bwa_seq_t *contig, const ubyte_t c, const int ori) {
	if (contig->full_len <= contig->len + 1) {
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

int has_n(const bwa_seq_t *read) {
	int i = 0;
	for (i = 0; i < read->len; i++) {
		if (read->seq[i] == 4 || read->seq[i] == 'n' || read->seq[i] == 'N')
			return 1;
	}
	return 0;
}

/**
 * Create a 'sudo-read' for overlapping alignment
 */
bwa_seq_t *new_seq(const bwa_seq_t *query, const int ol, const int shift) {
	bwa_seq_t *p = (bwa_seq_t*) malloc(sizeof(bwa_seq_t));
	int i = 0;
	if (ol + shift > query->len)
		return 0;
	for (i = 0; i < 16; i++) {
		p->bc[i] = 0;
	}
	p->tid = -1; // no assigned to a thread
	p->qual = p->strand = p->type = p->dummy = p->extra_flag = p->n_mm
			= p->n_gapo = p->n_gape = p->mapQ = p->score = p->n_aln = p->aln
					= p->n_multi = p->multi = p->sa = p->pos = p->c1 = p->c2
							= p->n_cigar = p->cigar = p->seQ = p->nm = p->md
									= NULL;

	p->shift = p->is_in_c_pool = 0;
	p->used = query->used;
	p->contig_id = query->contig_id;
	p->full_len = p->clip_len = p->len = ol;
	p->cursor = query->cursor;
	p->shift = query->shift;

	p->name = strdup((const char*) query->name);
	p->seq = (ubyte_t*) malloc(ol + 1);
	p->rseq = 0;
	memcpy(p->seq, query->seq + shift, ol);
	p->seq[ol] = '\0';
	p->rseq = (ubyte_t*) malloc(ol + 1);
	memcpy(p->rseq, p->seq, p->len);
	seq_reverse(p->len, p->rseq, 1);
	p->rseq[ol] = '\0';

	return p;
}

bwa_seq_t *new_mem_rev_seq(const bwa_seq_t *query, const int ol,
		const int shift) {
	bwa_seq_t *p = new_seq(query, ol, shift);
	ubyte_t *tmp = p->seq;
	p->seq = p->rseq;
	p->rseq = tmp;
	return p;
}

bwa_seq_t *new_rev_seq(const bwa_seq_t *query) {
	bwa_seq_t *p = (bwa_seq_t*) malloc(sizeof(bwa_seq_t));
	int i = 0;
	for (i = 0; i < 16; i++) {
		p->bc[i] = 0;
	}
	p->tid = -1; // no assigned to a thread
	p->qual = p->strand = p->type = p->dummy = p->extra_flag = p->n_mm
			= p->n_gapo = p->n_gape = p->mapQ = p->score = p->n_aln = p->aln
					= p->n_multi = p->multi = p->sa = p->pos = p->c1 = p->c2
							= p->n_cigar = p->cigar = p->seQ = p->nm = p->md
									= 0;

	p->used = p->shift = p->is_in_c_pool = 0;
	p->contig_id = -1;
	p->full_len = p->clip_len = p->len = query->len;
	p->cursor = query->cursor;
	p->shift = query->shift;
	p->name = query->name;

	p->seq = query->rseq;
	p->rseq = query->seq;
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
		c = "acgtn"[(int) c];
		fputc(c, tx_fp);
		if ((i + 1) % LINELEN == 0) {
			fputc('\n', tx_fp);
		}
	}
	if (i % LINELEN != 0)
		fputc('\n', tx_fp);
}

indexes *load_index(const char *fn) {
	indexes *in = (indexes*) malloc(sizeof(indexes));
	char *str = (char*) calloc(strlen(fn) + 10, 1);
	strcpy(str, fn);
	strcat(str, ".bwt");
	in->bwt[0] = bwt_restore_bwt(str);
	in->bwt[2] = bwt_restore_bwt(str);
	strcpy(str, fn);
	strcat(str, ".sa");
	bwt_restore_sa(str, in->bwt[2]);

	strcpy(str, fn);
	strcat(str, ".rbwt");
	in->bwt[1] = bwt_restore_bwt(str);
	in->bwt[3] = bwt_restore_bwt(str);
	strcpy(str, fn);
	strcat(str, ".rsa");
	bwt_restore_sa(str, in->bwt[3]);
	free(str);

	in->bns = bns_restore(fn);
	return in;
}

int same_q(const bwa_seq_t *query, const bwa_seq_t *seq) {
	unsigned int i = 0;
	if (!query || !seq || !seq->seq || !query->seq)
		return 0;
	if (query->len != seq->len)
		return 0;
	i = query->len;
	while (i--) {
		if (query->seq[i] != seq->seq[i])
			return 0;
	}
	return 1;
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
	start = (shift >= 0) ? 0 : abs(shift);
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
	start = (shift >= 0) ? shift : abs(shift);
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

int is_sub_seq_byte(const ubyte_t *query, const int q_len, const int shift,
		const bwa_seq_t *seq, int mismatches, const int ol) {
	unsigned int i = 0, start = 0, offset = 0, end = q_len;
	int nm;
	if (!query || !seq || !seq->seq)
		return NOT_FOUND;
	start = (shift >= 0) ? shift : abs(shift);
	if (ol) {
		end = (end > (start + ol)) ? (start + ol) : end;
	}
	if ((end - start) > seq->len)
		return NOT_FOUND;
	for (offset = 0; offset <= seq->len - (end - start); offset++) {
		nm = mismatches;
		for (i = start; i < end; i++) {
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

// Check whether two seqs share some small portion whose length is ol
int share_subseq(const bwa_seq_t *seq_1, const bwa_seq_t *seq_2,
		const int mismatches, const int ol) {
	int i = 0;
	if (!seq_1 || !seq_2)
		return 0;
	for (i = 0; i < seq_1->len - ol; i++) {
		if (is_sub_seq(seq_1, i, seq_2, mismatches, ol) != NOT_FOUND)
			return 1;
	}
	return 0;
}

int seq_ol(const bwa_seq_t *left_seq, const bwa_seq_t *right_seq, const int ol,
		int mismatches) {
	int i = 0;
	if (!left_seq || !right_seq || ol <= 0 || ol > left_seq->len || ol
			> right_seq->len)
		return 0;
	for (i = 0; i < ol; i++) {
		if (left_seq->seq[i + (left_seq->len - ol)] != right_seq->seq[i]) {
			mismatches--;
		}
		if (mismatches < 0)
			return 0;
	}
	return 1;
}

int find_ol(const bwa_seq_t *left_seq, const bwa_seq_t *right_seq,
		const int mismatches) {
	int i = 0, min_len = 0, olpped = 0;
	if (!left_seq || !right_seq)
		return 0;
	min_len = left_seq->len;
	min_len = (left_seq->len < right_seq->len) ? min_len : right_seq->len;
	for (i = min_len; i > 0; i--) {
		olpped = seq_ol(left_seq, right_seq, i, mismatches);
		if (olpped) {
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

int has_rep_pattern(bwa_seq_t *read) {
	int i = 0, j = 0, is_rep = 1;
	ubyte_t c = 0, c2 = 0;
	for (i = 0; i < read->len - NO_REPEAT_LEN - 1; i++) {
		is_rep = 1;
		for (j = i; j < i + NO_REPEAT_LEN; j++) {
			c = read->seq[j];
			c2 = read->seq[j + 1];
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

void pe_reverse_seqs(bwa_seq_t *seqs, const int n_seqs) {
	int i = 0;
	bwa_seq_t *s;
	for (i = 0; i < n_seqs; i++) {
		s = &seqs[i];
		seq_reverse(s->len, s->seq, 0);
		//seq_reverse(s->len, s->rseq, 0);
	}
}

void destroy_index(indexes *in) {
	if (in) {
		bwt_destroy(in->bwt[0]);
		bwt_destroy(in->bwt[1]);
		bwt_destroy(in->bwt[2]);
		bwt_destroy(in->bwt[3]);
		bns_destroy(in->bns);
		free(in);
	}
}

void free_read_seq(bwa_seq_t *p) {
	if (p) {
		free(p->name);
		free(p->seq);
		free(p->rseq);
		free(p->qual);
		free(p->aln);
		free(p->md);
		free(p->multi);
		free(p->cigar);
		free(p);
	}
}
