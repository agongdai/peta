/*
 * seq_util.h
 *
 *  Created on: May 16, 2011
 *      Author: Carl
 */
#ifndef PESEQ_H_
#define PESEQ_H_

#include <stdio.h>
#include <glib.h>
#include "bwtaln.h"
#include "utils.h"

enum READ_STATUS {
	FRESH, USED, TRIED, DEAD, HANG, HAS_N, REPETITIVE, MULTI,
};

#ifdef __cplusplus
extern "C" {
#endif

void rescue_reads(bwa_seq_t *seqs, const int n_seqs);
gint cmp_reads_by_name(gpointer a, gpointer b);
gint cmp_reads_by_contig_id(gpointer a, gpointer b);
void save_fq(const bwa_seq_t *seqs, const char *fp_fn, const uint16_t ol);
bwa_seq_t *merge_seq_to_right(bwa_seq_t *s1, bwa_seq_t *s2, const int gap);
int trun_seq(bwa_seq_t *s, const int shift);
bwa_seq_t *merge_seq_to_left(bwa_seq_t *s2, bwa_seq_t *s1, const int gap);
int is_left_mate(const char *seq_id);
int is_right_mate(const char *seq_id);
int is_mates(const char *left, const char *right);
char *get_mate_name(const char *seq_id);
index64 get_mate_index(const index64 seq_id);
bwa_seq_t *get_mate(const bwa_seq_t *s, bwa_seq_t *seqs);
bwa_seq_t *get_right_mate(const bwa_seq_t *left, bwa_seq_t *seqs);
bwa_seq_t *get_left_mate(const bwa_seq_t *right, bwa_seq_t *seqs);
void p_query(const char *header, const bwa_seq_t *q);
void p_readarray(const GPtrArray *ra, const int all);
void p_ctg_seq(const char *header, const bwa_seq_t *q);
void p_seq(const char *header, const ubyte_t *seq, const int len);
void ext_con(bwa_seq_t *contig, const ubyte_t c, const int ori);
void ext_que(bwa_seq_t *q, const ubyte_t c, const int left_max_ctg_id);
int has_n(const bwa_seq_t *read, int max);
bwa_seq_t *new_seq(const bwa_seq_t *query, const int ol, const int shift);
void set_rev_com(bwa_seq_t *s);
bwa_seq_t *blank_seq(const int len);
bwa_seq_t *new_mem_rev_seq(const bwa_seq_t *query, const int ol,
		const int shift);
bwa_seq_t *new_rev_seq(const bwa_seq_t *query);
bwa_seq_t *merge_seq(bwa_seq_t *s1, bwa_seq_t *s2, const int shift);
void map(bwa_seq_t *bwa_seq);
index64 get_index(const char *seq_id);
void save_con(const char *header, const bwa_seq_t *contig, FILE *tx_fp);
int save_unpaired_seqs(const char *part_solid_fn, bwa_seq_t *seqs,
		const int n_seqs);
int same_q(const bwa_seq_t *query, const bwa_seq_t *seq);
int same_bytes(const ubyte_t *s, const int k);
int similar_seqs(const bwa_seq_t *query, const bwa_seq_t *seq,
		const int mismatches, const int max_n_gaps, const int score_mat,
		const int score_mis, const int score_gap);
int is_biased_q(const bwa_seq_t *query);
int is_sub_seq_aln(const ubyte_t *query, const int q_len, const int shift,
		const int offset, const bwa_seq_t *seq, int mismatches, const int ol);
int is_sub_seq(const bwa_seq_t *query, const int shift, const bwa_seq_t *seq,
		int mismatches, const int ol_len);
int is_sub_seq_byte(const ubyte_t *query, const int q_len, const int shift,
		const bwa_seq_t *seq, int mismatches, const int ol);
int share_subseq_byte(const ubyte_t *seq_1, const int len,
		const bwa_seq_t *seq_2, const int mismatches, const int ol);
int seq_ol(const bwa_seq_t *left_seq, const bwa_seq_t *right_seq, const int ol,
		int mismatches);
int find_ol(const bwa_seq_t *left_seq, const bwa_seq_t *right_seq,
		const int mismatches);
void switch_fr(bwa_seq_t *s);
int find_ol_within_k(const bwa_seq_t *mate, const bwa_seq_t *temp,
		const int mismatches, const int min_len, const int max_len,
		const int ori);
int is_repetitive_q(const bwa_seq_t *query);
void pe_reverse_seqs(bwa_seq_t *seqs, const int n_seqs);
int is_paired(const bwa_seq_t *read, const int ori);
int get_mismatches_on_ol(const bwa_seq_t *query, const bwa_seq_t *seq,
		const int ol, const int max);
void free_read_seq(bwa_seq_t *p);
int has_rep_pattern(const bwa_seq_t *read);
int smith_waterman(const bwa_seq_t *seq_1, const bwa_seq_t *seq_2,
		const int score_mat, const int score_mis, const int score_gap,
		const int min_acceptable_score);

#ifdef __cplusplus
}
#endif

#endif /* PESEQ_H_ */
