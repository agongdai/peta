/*
 * tpl.h
 *
 *  Created on: 13-Oct-2011
 *      Author: carl
 */

#ifndef TPL_HPP_
#define TPL_HPP_
#include <stdio.h>
#include <glib.h>
#include "bwtaln.h"
#include "k_hash.h"

#define MAX_N_tpl_OUT 				7
#define MAX_N_tpl_IN 				7
#define INIT_N_GAPS 				3
#define INIT_N_READ_USED 			7
#define INIT_N_READ_PAIRED 			63
#define UNUSED_CONTIG_ID			-1
#define INVALID_CONTIG_ID			-2	// Means a reads has been tried to extend, but fail, never use later

using namespace std;

typedef GPtrArray tplarray;
typedef GPtrArray readarray;

typedef struct {
	bwa_seq_t *ctg;			// Sequence
	bwa_seq_t *r_tail;		// Right tail, virtual
	bwa_seq_t *l_tail;		// Left tail, virtual
	GPtrArray *b_juncs;		// Junctions with this template as branch
	GPtrArray *m_juncs;		// Junctions with this template as main
	GPtrArray *points;		// Branching points
	int id; 				// template id
	int32_t len;			// Length
	int8_t alive;			// Whether alive
	int8_t not_covered;		// Whether some region not covered by any reads.
	int8_t is_root;			// Whether it's a root node in the graph
	int8_t ori;				// Orientation
	float cov;				// Coverage
	uint64_t tid;			// Thread id
	bwa_seq_t *start_read;	// The starting kmer
	bwa_seq_t *last_read;	// The last read
	uint32_t kmer_freq;		// Sum of all kmer frequencies
	GPtrArray *reads;		// Reads on it. Used for paired validation
	GPtrArray *tried;		// Reads that loaded to pool once, then removed, marked as TRIED
	GPtrArray *vertexes;	// The template be broken into vertexes linearly
	testing_info *tinfo;	// Information for testing
} tpl;

typedef struct {
	int s_index;
	int size;
	int ori;
} eg_gap;

#ifdef __cplusplus
extern "C" {
#endif

	gint cmp_tpl_by_id(gpointer a, gpointer b);
	int has_nearby_pairs(hash_table *ht, GPtrArray *tpls, tpl *t, int n_pairs);
	GPtrArray *rm_dup_reads_on_tpl(GPtrArray *reads);
	void switch_tpl_fr(tpl *t);
	int rev_com_on_tpl(tpl *t, int start, bwa_seq_t *read);
	void p_tpl(tpl *t);
	eg_gap *init_gap(int s_index, int size, int ori);
	void p_tpl_reads(tpl *t);
	void free_eg_gap(eg_gap *gap);
	tpl *new_tpl();
	int is_high_cov(tpl *t);
	void reset_tid(tpl *t);
	void destroy_tpl(tpl *t, int status);
	void keep_ctg_only(tpl *t);
	void free_eg_seq(tpl *t);
	void mv_unpaired_to_tried(bwa_seq_t *seqs, tpl *t, const int n_tpls);
	bwa_seq_t *get_pure_tail(tpl *t, int len, int ori);
	bwa_seq_t *get_tail(tpl *t, int len, const int ori) ;
	bwa_seq_t *cut_tpl_tail(tpl *t, const int tail_len, int pos, const int ori);
	void set_tail(tpl *branch, tpl *parent_eg, const int shift, const int tail_len, const int ori);
	void save_tpls(tplarray *pfd_ctg_ids, FILE *ass_fa, const int ori,
			const int p_all, const int min_len);
	void find_reads_ahead(tpl *t, const int read_len, int ol_len, int *n_reads, const int ori);
	int find_reads_at_tail(tpl *t, int len, int min, int ori);
	void upd_reads_after_truncate(tpl *t, int trun_len);
	float calc_tpl_cov(tpl *t, int start, int end, int read_len);
	int vld_tpl_mates(tpl *t1, tpl *t2, int start_2, int end_2,
			const int min_n_pairs);
	void reverse_locus(tpl *t);
	void upd_locus_on_tpl(tpl *t, int pre_t_len, int pre_n_reads);
	void add2tpl(tpl *t, bwa_seq_t *r, const int locus);
	GPtrArray *reads_on_seq(bwa_seq_t *seq, hash_table *ht, const int n_mismatch);
	void refresh_reads_on_tail(hash_table *ht, tpl *t, int mismatches);
	void reset_boundary_reads(tpl *t, const int ori);
	int should_at_which_side(bwa_seq_t *seqs, tpl *t, bwa_seq_t *r);
	void clear_tpl_tails(tpl *t);
	void correct_tpl_base(bwa_seq_t *seqs, tpl *t, const int read_len, int start, int end);
	GPtrArray *check_branch_tail(hash_table *ht, tpl *t, bwa_seq_t *query,
			int shift, int mismatches, int8_t status, int ori);
	GPtrArray *align_tpl_tail(hash_table *ht, tpl *t, bwa_seq_t *tail, int max, int shift,
			int mismatches, int8_t status, int ori);
	int binary_exist(GPtrArray *reads, bwa_seq_t *q);
	int has_pairs_on_tpl(hash_table *ht, tpl *t, const int n_pairs);
	int find_pairs(GPtrArray *reads_1, GPtrArray *reads_2, int t1_id, int t2_id,
			int start_2, int end_2, const int min_n_pairs);
	int paired_by_reads(bwa_seq_t *seqs, tpl *t_1, tpl *t_2, int n_pairs);
	void mark_init_reads_used(hash_table *ht, tpl *t, bwa_seq_t *read, int mismatches);
	void unfrozen_tried(tpl *t);
	void reset_unpaired_reads(bwa_seq_t *seqs, tpl *t);
	int find_tried_tpl(tpl *t, const int tid);
	void reset_reads_status(GPtrArray *reads, int status);
	void add2tried(tpl *t, bwa_seq_t *r);
	void rm_from_tpl(tpl *t, int index);
	void truncate_tpl(tpl *t, int len, int ori);
	void rm_from_tried(tpl *t, const int rm_id);
	void unhold_reads_array(GPtrArray *reads);
	bwa_seq_t *get_tpl_ctg_wt(tpl *t, int *l_len, int *r_len, int *t_len);
	void refresh_tpl_reads(hash_table *ht, tpl *t, int mismatches);

#ifdef __cplusplus
}
#endif

#endif /* TPL_H_ */
