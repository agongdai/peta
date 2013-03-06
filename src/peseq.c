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
			kroundup32(fa->m);
			// Adjust the size to nearest 2^k
			//			fprintf(stderr, "[fa->l, fa->m]: [%zd, %zd]\n", fa->l, fa->m);
			fa->s = (char*) realloc(fa->s, fa->m);
		}
		fa->s[fa->l++] = c;
		counter++;
	}
	fclose(fp);
	return fa;
}

void p_seq(const char *header, const ubyte_t *seq, const int len) {
	int i = 0;
	printf("[%s] ", header);
	for (i = 0; i < len; i++) {
		if (seq[i] > 4)
			printf("%c", seq[i]);
		else
			printf("%c", "acgtn"[(int) seq[i]]);
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
	int i = 0;
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
		if (q->rev_com) {
			if (q->rseq[i] > 4)
				printf("%c", q->rseq[i]);
			else
				printf("%c", "acgtn"[(int) q->rseq[i]]);
		} else {
			if (q->seq[i] > 4)
				printf("%c", q->seq[i]);
			else
				printf("%c", "acgtn"[(int) q->seq[i]]);
		}
	}
	if (q->is_in_c_pool)
		printf(" [pool: %d]", q->is_in_c_pool);
	else
		printf(" [no_pool]");
	if (q->is_in_c_pool)
		printf(" [m_pool]");
	else
		printf(" [no_m_pool]");
	if (q->rev_com)
		printf(" [rev_com]");
	else
		printf(" [not_rev_com]");
	printf(" [%d: %d, %d] [tid: %d]", q->status, q->contig_id, q->shift, q->tid);
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
									= 0;

	p->is_in_c_pool = p->is_in_m_pool = 0;
	p->status = query->status;
	p->contig_id = query->contig_id;
	p->full_len = p->clip_len = p->len = ol;
	p->cursor = query->cursor;
	p->shift = query->shift;
	p->rev_com = query->rev_com;

	if (query->name)
		p->name = strdup((const char*) query->name);
	else
		p->name = (char*) calloc(1, sizeof(char));
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

void set_rev_com(bwa_seq_t *s) {
	free(s->rseq);
	s->rseq = (ubyte_t*) malloc(s->len + 1);
	memcpy(s->rseq, s->seq, s->len);
	seq_reverse(s->len, s->rseq, 1);
	s->rseq[s->len] = '\0';
}

bwa_seq_t *blank_seq() {
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

	p->is_in_c_pool = p->is_in_m_pool = 0;
	p->status = 0;
	p->contig_id = 0;
	p->full_len = p->clip_len = p->len = 0;
	p->cursor = 0;
	p->shift = 0;
	p->rev_com = 0;

	p->name = NULL;
	p->seq = (ubyte_t*) calloc(1, sizeof(ubyte_t));
	p->rseq = (ubyte_t*) calloc(1, sizeof(ubyte_t));

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

	p->status = p->shift = p->is_in_c_pool = p->is_in_m_pool = 0;
	p->contig_id = query->contig_id;
	p->full_len = p->clip_len = p->len = query->len;
	p->cursor = query->cursor;
	p->shift = query->shift;
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
		c = "acgtn"[(int) c];
		fputc(c, tx_fp);
		if ((i + 1) % LINELEN == 0) {
			fputc('\n', tx_fp);
		}
	}
	if (i % LINELEN != 0)
		fputc('\n', tx_fp);
}

void save_read(const char *header, const bwa_seq_t *read, FILE *read_fp) {
	int i = 0;
	char c = to_upper_lower('N');
	if (!read)
		return;
	fputs(header, read_fp);
	for (i = 0; i < read->len; i++) {
		c = read->seq[i];
		c = "ACGTN"[(int) c];
		fputc(c, read_fp);
	}
	fputc('\n', read_fp);
}

int save_unpaired_seqs(const char *part_solid_fn, bwa_seq_t *seqs,
		const int n_seqs) {
	int n_unpaired = 0, i = 0;
	bwa_seq_t *s = NULL;
	FILE *solid = NULL;
	char *h = malloc(BUFSIZ);

	solid = xopen(part_solid_fn, "w");
	for (i = 0; i < n_seqs; i++) {
		s = &seqs[i];
		//		if (s->status != USED && s->status != DEAD) {
		sprintf(h, ">%d\n", n_unpaired);
		save_read(h, s, solid);
		n_unpaired++;
		//		}
	}
	free(h);
	fclose(solid);
	return n_unpaired;
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
		//printf("Previous row: \n");
		for (j = 0; j < columns; j++) {
			//	printf("%d,", previous_row[j]);
			previous_row[j] = current_row[j];
			max_score = current_row[j] > max_score ? current_row[j] : max_score;
		}
		//printf("\n");
		//printf("Current row: \n");
		//for (j = 0; j < columns; j++) {
		//	printf("%d,", current_row[j]);
		//}
		//printf("\n");
		//printf("Max score: %d \n", max_score);
		// If the minimal acceptable score is not reachable, stop and return.
		if ((max_score + (rows - i) * score_mat) < min_acceptable_score) {
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
	if (abs(query->len - seq->len) > mismatches)
		return 0;
	min_len = query->len;
	min_len = min_len > seq->len ? seq->len : min_len;
	min_acceptable_score = min_len * score_mat + mismatches * score_mis + abs(
			query->len - seq->len) * score_gap;
	similarity_score = smith_waterman(query, seq, score_mat, score_mis,
			score_gap, min_acceptable_score);
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

int find_ol_within_k(const bwa_seq_t *mate, const bwa_seq_t *template,
		const int mismatches, const int min_len, const int max_len,
		const int ori) {
	int i = 0, olpped = 0;
	if (!mate || !template)
		return 0;
	for (i = max_len; i > min_len; i--) {
		olpped = ori ? seq_ol(mate, template, i, mismatches) : seq_ol(template,
				mate, i, mismatches);
		if (olpped > 0) {
			olpped = i;
			break;
		}
	}
	return olpped;
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
		if (olpped > 0) {
			return i;
		}
	}
	for (i = min_len; i >= mismatches * 2; i--) {
		olpped = seq_ol(left_seq, right_seq, i, mismatches);
		if (olpped > 0) {
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

int is_biased_q(const bwa_seq_t *query) {
	int *c = (int*) calloc(5, sizeof(int));
	int i = 0, is_biased = 0;
	int thre = query->len / 4 / 3;
	for (i = 0; i < query->len; i++) {
		c[query->seq[i]]++;
	}
	for (i = 0; i < 4; i++) {
		if (c[i] < thre) {
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
