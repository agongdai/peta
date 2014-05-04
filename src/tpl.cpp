/*
 * tpl.c
 *
 *  Created on: 06-Mar-2012
 *      Author: carl
 */

#include <stdint.h>
#include <glib.h>
#include <cstdio>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>
#include "utils.h"
#include "tpl.hpp"
#include "peseq.h"
#include "bwtaln.h"
#include "k_hash.h"

using namespace std;

eg_gap *init_gap(int s_index, int size, int ori) {
	eg_gap *gap = (eg_gap*) malloc(sizeof(eg_gap));
	gap->s_index = s_index;
	gap->size = size;
	gap->ori = ori;
	return gap;
}

void p_anchor(char *header, anchor *a) {
	show_debug_msg(header, "HASH = %"ID64 "; Anchor: template [%d, %d] @ %d size %d; from %d \n", a->hash, a->t->id, a->t->len, a->locus, a->size, a->from);
}

void p_tpl(tpl *t) {
	if (!t) {
		show_debug_msg(__func__, "---- Template is NULL ----\n");
	}
	show_debug_msg(__func__, "---- Template %d alive: %d; visited %d ----\n",
			t->id, t->alive, t->visited);
	show_debug_msg(__func__, "\t Length: %d; Pairs: %.2f \n", t->len, t->pair_pc);
	if (t->reads)
		show_debug_msg(__func__, "\t Reads: %d\n", t->reads->len);
	show_debug_msg(__func__, "\t Coverage: %.2f\n", t->cov);
	if (t->tried)
		show_debug_msg(__func__, "\t Tried: %d\n", t->tried->len);
	if (t->m_juncs)
		show_debug_msg(__func__, "\t Junctions as main: %d \n", t->m_juncs->len);
	if (t->b_juncs)
		show_debug_msg(__func__, "\t Junctions as branch: %d \n",
				t->b_juncs->len);
	if (t->start_read)
		p_query("start", t->start_read);
	if (t->l_tail)
		p_query("ltail", t->l_tail);
	if (t->r_tail)
		p_query("rtail", t->r_tail);
	if (t->ctg)
		p_ctg_seq("contg", t->ctg);
	show_debug_msg(__func__, "---- Template %d ---- \n", t->id);
}

void p_tpl_reads(tpl *t) {
	GPtrArray *reads = t->reads;
	bwa_seq_t *r = NULL;
	int i = 0, j = 0;
	int read_len = 0;
	if (!t->reads || t->reads->len == 0) {
		show_debug_msg(__func__, "No reads on template [%d, %d] \n", t->id,
				t->len);
		return;
	}
	r = (bwa_seq_t*) g_ptr_array_index(t->reads, 0);
	read_len = r->len;
	g_ptr_array_sort(t->reads, (GCompareFunc) cmp_reads_by_contig_locus);
	printf("\n==== Here are the reads used on template [%d, %d] ==== \n",
			t->id, t->len);
	for (i = 0; i < reads->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(reads, i);
		//if (r->pos == IMPOSSIBLE_NEGATIVE)
		//	continue;
		if (i % 15 == 0) {
			for (j = 0; j < read_len; j++) {
				printf(" ");
			}
			for (j = 0; j < t->len; j++) {
				printf("%c", "ACGTN"[t->ctg->seq[j]]);
			}
			printf("\n");
		}
		for (j = 0; j < r->contig_locus + read_len; j++) {
			printf(" ");
		}
		for (j = 0; j < r->len; j++) {
			if (r->rev_com) {
				printf("%c", "ACGTN"[(int) r->rseq[j]]);
			} else {
				printf("%c", "ACGTN"[(int) r->seq[j]]);
			}
		}
		printf("\t [%s, %d] \t", r->name, r->len);
		printf(" [status: %d] ", r->status);
		if (r->rev_com)
			printf(" [    <<<<@%d]", r->pos);
		else
			printf(" [>>>>    @%d]", r->pos);
		printf(" [%d: %d]", r->contig_id, r->contig_locus);
		printf(" [cursor: %d]", r->cursor);
		printf("\n");
	}
	printf("==== End of printing reads used on template [%d, %d] ==== \n",
			t->id, t->len);
}

void free_eg_gap(eg_gap *gap) {
	if (gap)
		free(gap);
}

tpl *new_tpl() {
	tpl *t = (tpl*) malloc(sizeof(tpl));
	t->ctg = NULL;
	t->l_tail = NULL;
	t->r_tail = NULL;
	t->m_juncs = NULL;
	t->b_juncs = NULL;
	t->len = 0;
	t->id = 0;
	t->alive = 1;
	t->pair_pc = 0.0;
	t->visited = 0;
	t->ori = 0;
	t->start_read = NULL;
	t->last_read = NULL;
	t->tid = 0;
	t->kmer_freq = 0;
	t->cov = 0.0;
	t->vertexes = g_ptr_array_sized_new(0);
	t->reads = g_ptr_array_sized_new(0);
	t->tried = g_ptr_array_sized_new(0);
	return t;
}

/**
 * The t->tried stores the template ids which have been tried to merge with current template.
 * Used by merge_tpls_by_hash
 * @Desperate
 */
int find_tried_tpl(tpl *t, const int tid) {
	int i = 0, found = 0;
	bwa_seq_t *r = NULL;
	if (tid < 0 || !t->tried)
		return 0;
	for (i = 0; i < t->tried->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(t->tried, i);
		if (r->contig_id == tid)
			return 1;
	}
	return 0;
}

/**
 * Add to TRIED, pretend from using for this template
 */
void add2tried(tpl *t, bwa_seq_t *r) {
	r->contig_id = t->id;
	r->contig_locus = -1;
	r->status = TRIED;
	r->pos = IMPOSSIBLE_NEGATIVE;
	r->rev_com = 0;
	g_ptr_array_add(t->tried, r);
}

void reset_reads_status(GPtrArray *reads, int status) {
	index64 i = 0;
	bwa_seq_t *r = NULL;
	for (i = 0; i < reads->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(reads, i);
		r->status = status;
		reset_read(r);
	}
}

void reset_unpaired_reads(bwa_seq_t *seqs, tpl *t) {
	bwa_seq_t *mate = NULL, *read = NULL;
	int i = 0;
	for (i = 0; i < t->reads->len; i++) {
		read = (bwa_seq_t*) g_ptr_array_index(t->reads, i);
		mate = get_mate(read, seqs);
		if (mate->status != USED || mate->contig_id != t->id) {
			reset_to_fresh(mate);
		}
	}
}

/**
 * Reset those not USED reads to fresh.
 * And clear the 'tried' reads on the template
 */
void unfrozen_tried(tpl *t) {
	index64 i = 0;
	bwa_seq_t *r = NULL;
	if (t->tried) {
		for (i = 0; i < t->tried->len; i++) {
			r = (bwa_seq_t*) g_ptr_array_index(t->tried, i);
			if (r->status != USED) {
				reset_to_fresh(r);
			}
		}
		while (t->tried->len > 0)
			g_ptr_array_remove_index_fast(t->tried, 0);
	}
}

/**
 * For reads on the template, if its mate is used another template, add its mate to TRIED.
 * Used for merging. To save time.
 *
 * If there are multiple reads spanning two templates, add only one
 */
void mv_unpaired_to_tried(bwa_seq_t *seqs, tpl *t, const int n_tpls) {
	int i = 0;
	bwa_seq_t *r = NULL, *m = NULL;
	uint8_t *flag = NULL;
	if (!t || !t->reads || t->reads->len == 0)
		return;
	if (!t->tried)
		t->tried = g_ptr_array_sized_new(4);
	while (t->tried->len > 0)
		g_ptr_array_remove_index_fast(t->tried, 0);
	// flag[0] indicates template 0 has been tried with current template
	flag = (uint8_t*) calloc(n_tpls, sizeof(uint8_t));
	for (i = 0; i < t->reads->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(t->reads, i);
		m = get_mate(r, seqs);
		if (m->status == USED && m->contig_id != r->contig_id) {
			if (flag[m->contig_id] == 0) {
				g_ptr_array_add(t->tried, m);
				flag[m->contig_id] = 1;
			}
		}
	}
	free(flag);
}

/**
 * Check whether the template is with high coverage
 */
int is_high_cov(tpl *t) {
	float cov = 0.0;
	int read_len = 0;
	bwa_seq_t *r = NULL;
	if (!t || !t->reads || t->len <= 0 || t->reads->len <= 0)
		return 0;
	if (t->len < 200)
		return 0;
	r = (bwa_seq_t*) g_ptr_array_index(t->reads, 0);
	read_len = r->len;
	cov = ((float) read_len) * ((float) t->reads->len) / ((float) t->len);
	//show_debug_msg(__func__, "Coverage of %d: %d*%d/%d = %.2f \n", t->id, read_len, t->reads->len, t->len, cov);
	if (cov > HIHG_COV_THRE)
		return 1;
	return 0;
}

/**
 * Get pure head/tail of a template
 */
bwa_seq_t *get_pure_tail(tpl *t, int len, int ori) {
	bwa_seq_t *tail = NULL;
	if (t->len < len)
		return NULL;
	tail = ori ? new_seq(t->ctg, len, 0) : new_seq(t->ctg, len, t->len - len);
	return tail;
}

/**
 * Get the tail for extension.
 * The only chance that the template is shorter than len is after timmed by junction.
 */
bwa_seq_t *get_tail(tpl *t, int len, const int ori) {
	bwa_seq_t *tail = NULL;
	bwa_seq_t *real_seq = NULL;
	int l_len = 0, r_len = 0, t_len = 0;
	// Must be something wrong!
	if (!t || t->len < 0 || len < 0)
		return NULL;
	//p_tpl(t);
	real_seq = get_tpl_ctg_wt(t, &l_len, &r_len, &t_len);
	if (t->start_read->rev_com) switch_fr(t->start_read);
	//p_query(__func__, real_seq);
	if (real_seq->len < len) {
		show_debug_msg(__func__, "[WARNING] Tail of template [%d, %d] shows a wrong sequence.\n",
				t->id, t->len);
		tail = new_seq(real_seq, real_seq->len, 0);
		//p_tpl(t);
	} else {
		tail = ori ? new_seq(real_seq, len, 0) : new_seq(real_seq, len,
				real_seq->len - len);
	}
	if (t->start_read->rev_com) switch_fr(t->start_read);
	//p_query(__func__, tail);
	bwa_free_read_seq(1, real_seq);
	return tail;
}

gint cmp_tpl_by_id(gpointer a, gpointer b) {
	tpl *c_a = *((tpl**) a);
	tpl *c_b = *((tpl**) b);
	return (c_a->id - c_b->id);
}

gint cmp_tpl_by_rev_pair_pc(gpointer a, gpointer b) {
	tpl *c_a = *((tpl**) a);
	tpl *c_b = *((tpl**) b);
	gint larger = c_b->pair_pc - c_a->pair_pc >= 0 ? 1 : -1;
	return larger;
}

gint cmp_anchor_by_locus(gpointer a, gpointer b) {
	anchor *c_a = *((anchor**) a);
	anchor *c_b = *((anchor**) b);
	return (c_a->locus - c_b->locus);
}

gint cmp_anchor_by_rev_locus(gpointer a, gpointer b) {
	anchor *c_a = *((anchor**) a);
	anchor *c_b = *((anchor**) b);
	return (c_b->locus - c_a->locus);
}

gint cmp_anchor_by_from(gpointer a, gpointer b) {
	anchor *c_a = *((anchor**) a);
	anchor *c_b = *((anchor**) b);
	return (c_a->from - c_b->from);
}

gint cmp_anchor_by_hash(gpointer a, gpointer b) {
	anchor *c_a = *((anchor**) a);
	anchor *c_b = *((anchor**) b);
	gint rs = (c_a->hash >= c_b->hash) ? 1 : -1;
	return rs;
}

/**
 * Shrink an anchor with sized 'shifted'
 * Example: anchor 'AAAACCCCCCCC'
 * Shifted 4: anchor 'CCCCCCCC'
 */
void shrink_anchor_right(anchor *a, int shifted) {
	int s = min(shifted, a->size);
	a->locus += s; a->size -= s; a->from += s;
}

/**
 * Mark the anchor not usable anymore.
 */
void disable_anchor(anchor *a) {
	a->locus = -1; a->from = -1; a->size = 0; a->from = NULL;
}

/**
 * If two templates share regions longer than 11bp, there are multiple anchors nearby,
 * here merge them;
 * By 'cmp_anchor_by_hash', the anchors are ordered first by template id, then by locus, increasingly.
 */
void compact_anchors(GPtrArray *anchors) {
	int i = 0;
	anchor *a = NULL, *pre = NULL;
	g_ptr_array_sort(anchors, (GCompareFunc) cmp_anchor_by_hash);
	if (anchors && anchors->len > 0) {
		pre = (anchor*) g_ptr_array_index(anchors, anchors->len - 1);
		for (i = anchors->len - 2; i >= 0; i--) {
			a = (anchor*) g_ptr_array_index(anchors, i);
			if (a->t == pre->t && pre->locus - a->locus == 1 && pre->from - a->from == 1) {
				pre->locus--; pre->size++; pre->from--;
				free(a); g_ptr_array_remove_index_fast(anchors, i);
			} else pre = a;
		}
	}
	g_ptr_array_sort(anchors, (GCompareFunc) cmp_anchor_by_hash);
}

/**
 * Remove duplicates in an array
 */
GPtrArray *rm_dup_reads_on_tpl(GPtrArray *reads) {
	int i = 0;
	bwa_seq_t *pre = NULL, *post = NULL;
	GPtrArray *uni_reads = NULL;
	if (!reads || reads->len < 2) {
		return reads;
	}
	uni_reads = g_ptr_array_sized_new(reads->len);
	g_ptr_array_sort(reads, (GCompareFunc) cmp_reads_by_contig_id);
	pre = (bwa_seq_t*) g_ptr_array_index(reads, 0);
	g_ptr_array_add(uni_reads, pre);
	for (i = 1; i < reads->len; i++) {
		post = (bwa_seq_t*) g_ptr_array_index(reads, i);
		if (post->contig_id != pre->contig_id) {
			g_ptr_array_add(uni_reads, post);
			pre = post;
		}
	}
	g_ptr_array_free(reads, TRUE);
	return uni_reads;
}

/**
 * Test whether the read is reverse-complemented on the template.
 * Assumption: the read is for sure on the template, at least partially on it.
 */
int rev_com_on_tpl(tpl *t, int start, bwa_seq_t *read) {
	int i = 0, n_mis = 0;
	for (i = 0; i < read->len; i++) {
		if (start + i < 0 || start + i >= t->len)
			break;
		if (t->ctg->seq[start + i] != read->seq[i])
			n_mis++;
	}
	if (n_mis >= 10)
		return 1;
	else
		return 0;
}

/**
 * Get the nearest index of read with name as the query
 */
int binary_search_read(GPtrArray *reads, bwa_seq_t *q) {
	int start = 0, end = reads->len, middle = (start + end) / 2;
	bwa_seq_t *r = NULL, *m = NULL;
	//p_readarray(reads, 1);
	//p_query(__func__, q);
	while (start < end - 1) {
		middle = (start + end) / 2;
		//show_debug_msg(__func__, "[%d, %d, %d]\n", start, middle, end);
		m = (bwa_seq_t*) g_ptr_array_index(reads, middle);
		//p_query("MIDDLE", m);
		if (atoi(m->name) - atoi(q->name) == 0) {
			return middle;
		} else {
			if (atoi(m->name) - atoi(q->name) > 0) {
				end = middle;
			} else {
				start = middle + 1;
			}
		}
		//show_debug_msg(__func__, "[%d, %d, %d]\n", start, middle, end);
	}
	//show_debug_msg(__func__, "Middle: %d\n", middle);
	return start;
}

/**
 * Check whether a read exists in an array
 */
int binary_exist(GPtrArray *reads, bwa_seq_t *q) {
	bwa_seq_t *r = NULL;
	int index = binary_search_read(reads, q);
	if (index < 0 || index >= reads->len)
		return 0;
	r = (bwa_seq_t*) g_ptr_array_index(reads, index);
	if (r == q)
		return 1;
	return 0;
}

/**
 * Check whether there are at least n_pairs of paired reads on the template
 */
int has_pairs_on_tpl(hash_table *ht, tpl *t, const int n_pairs) {
	bwa_seq_t *r = NULL, *m = NULL;
	int i = 0;
	int pairs = 0;
	for (i = 0; i < t->reads->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(t->reads, i);
		m = get_mate(r, ht->seqs);
		if (m->status == USED && m->contig_id == t->id) {
			pairs++;
			//p_query("PAIRED READ", r);
			//p_query("PAIRED MATE", m);
		}
		// A pair would be counted twice
		if (pairs >= n_pairs * 2)
			return 1;
	}
	return 0;
}

/**
 * Find pairs between two lists
 * If t1_id is not 0, reads on reads_1 must has contig_id == t1_id
 */
int find_pairs(GPtrArray *reads_1, GPtrArray *reads_2, int t1_id, int t2_id,
		int start_2, int end_2, const int min_n_pairs) {
	bwa_seq_t *r1 = NULL, *r2 = NULL;
	uint32_t i = 0, j = 0, index_1 = 0, index_2 = 0;
	int n_pairs = 0;
	if (!reads_1 || !reads_2 || reads_1->len == 0 || reads_2->len == 0)
		return 0;
	r1 = (bwa_seq_t*) g_ptr_array_index(reads_2, 0);
	//show_debug_msg(__func__, "Searching %s in reads_1 \n", r1->name);
	index_1 = binary_search_read(reads_1, r1);
	//show_debug_msg(__func__, "Index_1: %d\n", index_1);

	r2 = (bwa_seq_t*) g_ptr_array_index(reads_1, 0);
	//show_debug_msg(__func__, "Searching %s in reads_2 \n", r2->name);
	index_2 = binary_search_read(reads_2, r2);
	//show_debug_msg(__func__, "Index_1: %d\n", index_2);

	//show_debug_msg(__func__, "Reads 1 \n");
	//p_readarray(reads_1, 1);
	//show_debug_msg(__func__, "Reads 2 \n");
	//p_readarray(reads_2, 1);

	for (i = index_1; i < reads_1->len; i++) {
		r1 = (bwa_seq_t*) g_ptr_array_index(reads_1, i);
		// Only check template id when the id is not 0
		if (t1_id && r1->contig_id != t1_id)
			continue;
		if (index_2 >= reads_2->len - 1)
			break;
		//p_query("READ_1", r1);
		for (j = index_2; j < reads_2->len; j++) {
			r2 = (bwa_seq_t*) g_ptr_array_index(reads_2, j);
			if (t2_id && (r2->contig_id != t2_id))
				continue;
			if (r2->contig_locus < start_2 || r2->contig_locus > end_2)
				continue;
			//p_query("READ_2", r2);
			if (is_mates(r1->name, r2->name)) {
				//p_query(__func__, r1);
				//p_query(__func__, r2);
				// Ignore if either one of them is like: "AAATAAAAAA"
				if (is_biased_q(r1) || is_biased_q(r2))
					continue;
				n_pairs++;
			}
			if (atoi(r2->name) - atoi(r1->name) >= 1) {
				index_2 = j;
				break;
			}
		}
		if (n_pairs >= min_n_pairs)
			return 1;
	}
	return 0;
}

/**
 * From the attributes on the reads, check whether two templates have pairs.
 */
int paired_by_reads(bwa_seq_t *seqs, tpl *t_1, tpl *t_2, int n_pairs) {
	int i = 0, found_pairs = 0;
	int loop_id = 0, pair_id = 0;
	tpl *t = NULL;
	bwa_seq_t *r = NULL, *m = NULL;
	if (t_1->reads->len > t_2->reads->len) {
		t = t_2;
		loop_id = t_2->id;
		pair_id = t_1->id;
	} else {
		t = t_1;
		loop_id = t_1->id;
		pair_id = t_2->id;
	}
	for (i = 0; i < t->reads->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(t->reads, i);
		m = get_mate(r, seqs);
		if (r->status == USED && m->status == USED && r->contig_id == loop_id
				&& m->contig_id == pair_id) {
			found_pairs++;
			//p_query(__func__, r);
			//p_query(__func__, m);
		}
		if (found_pairs >= n_pairs)
			return 1;
	}
	return 0;
}

int count_pairs_at_hit() {
	int n = 0;
	return n;
}

/**
 * Check whether there are reads before read_len
 */
void find_reads_ahead(tpl *t, const int read_len, int ol_len, int *n_reads,
		const int ori) {
	int i = 0, start = 0, end = 0, use_right_read = 0;
	int n = 0, left_ext_len = 0;
	bwa_seq_t *r = NULL;
	*n_reads = 0;
	// Refer to pool.cpp->forward to see the contig_locus during extension
	start = ori ? 0 : t->len - read_len - ol_len;
	end = ori ? read_len + ol_len : t->len;
	//show_debug_msg(__func__, "Start~End: [%d, %d] \n", start, end);
	for (i = t->reads->len - 1; i >= 0; i--) {
		r = (bwa_seq_t*) g_ptr_array_index(t->reads, i);
		if (r->status != USED || r->contig_id != t->id)
			continue;
		//show_debug_msg(__func__, "Start~End: [%d, %d]\n", start, end);
		if (r->contig_locus >= start && r->contig_locus <= end) {
			//p_query(__func__, r);
			n++;
		}
	}
	*n_reads = n;
}

/**
 *
 */
int find_reads_at_tail(tpl *t, int len, int min, int ori) {
	int i = 0, n = 0;
	bwa_seq_t *r = NULL;
	//p_tpl_reads(t);
	if (ori == 0) {
		for (i = t->reads->len - 1; i >= 0; i--) {
			r = (bwa_seq_t*) g_ptr_array_index(t->reads, i);
			if (r->contig_locus + r->len >= t->len - len) {
				//p_query(__func__, r);
				n++;
			}
			if (n > min)
				return 1;
		}
	} else {
		for (i = 0; i < t->reads->len; i++) {
			r = (bwa_seq_t*) g_ptr_array_index(t->reads, i);
			if (r->contig_locus < len) {
				//p_query(__func__, r);
				n++;
			}
			if (n > min)
				return 1;
		}
	}
	show_debug_msg(__func__, "[%d, %d]: Len %d; Ori %d; n_reads: %d\n", t->id,
			t->len, len, ori, n);
	return 0;
}

/**
 * Assumption: the reads on template are sorted by contig_locus
 */
float calc_tpl_cov(tpl *t, int start, int end, int read_len) {
	bwa_seq_t *r = NULL;
	int i = 0, j = 0, ol = 0;
	float sum = 0.0, cov = 0.0;
	//p_tpl_reads(t);
	//show_debug_msg(__func__, "Template [%d, %d] %d~%d\n", t->id, t->len, start,
	//		end);
	if (start == 0 && t->l_tail)
		start -= t->l_tail->len;
	if (!t || end < 0 || end > t->len || end <= start || read_len <= 0)
		return 0.0;
	for (i = 0; i < t->reads->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(t->reads, i);
		//p_query(__func__, r);
		if (r->contig_locus + r->len >= start) {
			for (j = i; j < t->reads->len; j++) {
				r = (bwa_seq_t*) g_ptr_array_index(t->reads, j);
				ol = r->len;
				ol = r->contig_locus < start ? ol - (start - r->contig_locus)
						: ol;
				ol = (r->contig_locus + read_len) >= end ? ol
						- ((r->contig_locus + read_len) - end) : ol;
				//show_debug_msg(__func__, "%.2f + %d \n", sum, ol);
				sum += ol;
			}
			break;
		}
	}
	cov = sum / (float) (end - start);
	//show_debug_msg(__func__, "Calculated template coverage [%d, %d]: %.2f \n",
	//		t->id, t->len, cov);
	return cov;
}

void upd_reads_after_truncate(tpl *t, int trun_len) {
	int i = 0;
	bwa_seq_t *r = NULL;
	if (!t || !t->alive || !t->reads || trun_len <= 0)
		return;
	for (i = 0; i < t->reads->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(t->reads, i);
		if (r->status == USED && r->contig_id == t->id) {
			r->contig_locus -= trun_len;
		}
	}
}

/**
 * Find whether the two templates contain at least min_n_pairs pairs of mates.
 * Assume the reads on them are already sorted by name
 */
int vld_tpl_mates(tpl *t1, tpl *t2, int start_2, int end_2,
		const int min_n_pairs) {
	int is_valid = 0;
	is_valid = find_pairs(t1->reads, t2->reads, t1->id, t2->id, start_2, end_2,
			min_n_pairs);
	return is_valid;
}

/**
 * Determine the tail of a branch template.
 *
 * This is used when a branch template is connected to the locus of '^'.
 * When another third-layer branch template is connected to this branch template,
 * 	its own length may be not long enough, then it cuts some length from the main template.
 * If the length of the branch template is long enough, the virtual tail is not used.
 *
 * Notation: '==': normal sequence; '--': tail
 * tpl: 	---==============================
 * pos: 	                   ^
 * ori: 	0 (to the right), using the left virtual tail
 * partial: ---================
 * If v_tail_len <= partial length, return:
 *                  ===========
 * Otherwise, try to return:
 *          ---================
 */
bwa_seq_t *cut_tpl_tail(tpl *t, int pos, const int tail_len, const int ori) {
	bwa_seq_t *tail = NULL, *partial = NULL, *main_tail = NULL;
	int v_tail_len = 0;
	if (t->len < pos)
		pos = t->len;
	//return new_seq(t->ctg, t->len, 0);
	//p_tpl(t);
	//show_debug_msg(__func__, "POS: %d; TailLen: %d\n", pos, tail_len);
	if (ori) {
		partial = new_seq(t->ctg, t->len - pos, pos);
		main_tail = t->r_tail;
	} else {
		partial = new_seq(t->ctg, pos, 0);
		main_tail = t->l_tail;
	}
	// If the tpl is long, cut the tail directly
	if (partial->len >= tail_len) {
		if (ori)
			tail = new_seq(partial, tail_len, 0);
		else
			tail = new_seq(partial, tail_len, partial->len - tail_len);
	} else {
		// If the tpl has a virtual tail, try to get it
		if (main_tail && main_tail->len > 0) {
			v_tail_len = main_tail->len + partial->len;
			v_tail_len = (v_tail_len > tail_len) ? tail_len : v_tail_len;
			tail = blank_seq(v_tail_len);
			if (ori) {
				memcpy(tail->seq, partial->seq, sizeof(ubyte_t) * partial->len);
				memcpy(tail->seq + partial->len, main_tail->seq, sizeof(ubyte_t)
						* (v_tail_len - partial->len));
			} else {
				memcpy(tail->seq, main_tail->seq + (main_tail->len
						+ partial->len - v_tail_len), sizeof(ubyte_t)
						* (v_tail_len - partial->len));
				memcpy(tail->seq + (v_tail_len - partial->len), partial->seq,
						sizeof(ubyte_t) * partial->len);
			}
			tail->len = v_tail_len;
			set_rev_com(tail);
		} else
			tail = new_seq(partial, partial->len, 0);
	}
	bwa_free_read_seq(1, partial);
	return tail;
}

/**
 * Set the left/right tail of the branch template
 * The tail starts position 'shift' of the main template with length 'tail_len',
 * The length of the tail is by default read_len - 1.
 */
void set_tail(tpl *branch, tpl *main_tpl, const int shift, const int tail_len,
		const int ori) {
	bwa_seq_t *tmp = NULL;
	// The right/left tail would be used in cut_tpl_tail function,
	//	but will be replaced later. So save to tmp first, free the old one later.
	if (ori) {
		tmp = branch->r_tail;
		branch->r_tail = cut_tpl_tail(main_tpl, shift, tail_len, ori);
	} else {
		tmp = branch->l_tail;
		branch->l_tail = cut_tpl_tail(main_tpl, shift, tail_len, ori);
	}
	bwa_free_read_seq(1, tmp);
}

/**
 * Set the contig locus of reads to be template_length - real_locus
 * When extending to left, the locus value is not determined until the extension stops
 * So it is hard to determine whether the distance of a pair is within range if:
 * 		1. Extend to left
 * 		2. Extend to right
 * 		3. Extend to left (!!!Problematic)
 */
void reverse_locus(tpl *t) {
	int i = 0;
	bwa_seq_t *r = NULL;
	for (i = 0; i < t->reads->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(t->reads, i);
		r->contig_locus = t->len - r->contig_locus;
	}
}

/**
 * Update the contig_locus value after left extension
 */
void upd_locus_on_tpl(tpl *t, int pre_t_len, int pre_n_reads) {
	index64 i = 0;
	bwa_seq_t *r = NULL;

	if (!t || !t->reads || pre_n_reads < 0 || pre_t_len < 0)
		return;
	for (i = 0; i < pre_n_reads; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(t->reads, i);
		r->contig_locus += t->len - pre_t_len;
	}
	for (i = pre_n_reads; i < t->reads->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(t->reads, i);
		r->contig_locus = t->len - r->contig_locus;
	}
}

/**
 * Add a read to the template
 */
void add2tpl(tpl *t, bwa_seq_t *r, const int locus) {
	r->contig_id = t->id;
	r->contig_locus = locus;
	r->pos = IMPOSSIBLE_NEGATIVE;
	r->status = USED;
	g_ptr_array_add(t->reads, r);
}

/**
 * Reset the reads on the template.
 * Align the template kmer by kmer to find reads, but do not remove existing reads on it.
 * To add more reads not obtained due to hashing limitation when extending
 */
void refresh_tpl_reads(hash_table *ht, tpl *t, int start, int end, int mismatches) {
	bwa_seq_t *r = NULL, *seq = NULL, *window = NULL;
	int left_len = 0, counted_len = 0, right_len = 0, n_mis = 0, rev_com = 0;
	int i = 0, j = 0, not_covered_len = 0, ol = 0;
	ubyte_t base_1 = 0, base_2 = 0;
	GPtrArray *refresh = NULL, *hits = NULL;
	if (!t || !t->alive || !t->reads || t->reads->len < 0 || t->len < 0)
		return;
	//unfrozen_tried(t);

	seq = get_tpl_ctg_wt(t, &left_len, &right_len, &counted_len);
	// If it happens, means something wrong
	if (seq->len < ht->o->read_len) {
		show_debug_msg("[WARNING]",
				"The sequence with tails shorter than read length: [%d, %d] \n",
				t->id, t->len);
		//p_tpl(t);
		bwa_free_read_seq(1, seq);
		return;
	}

	//p_tpl(t);
	//p_query(__func__, seq);
	for (i = max(0, start + right_len); i <= min(end + right_len, seq->len) - ht->o->read_len; i++) {
		window = new_seq(seq, ht->o->read_len, i);
		hits = g_ptr_array_sized_new(4);

		//if (i == 2911) {
		//	char *test_name = (char*) malloc(sizeof(char) * 8);
		//	sprintf(test_name, "SOME");
		//	window->name = test_name;
		//	p_query(__func__, window);
		//}

		hits = find_both_fr_full_reads(ht, window, hits, mismatches);
		if (hits->len == 0)
			not_covered_len++;
		else
			not_covered_len = 0;
		for (j = 0; j < hits->len; j++) {
			r = (bwa_seq_t*) g_ptr_array_index(hits, j);
			r->pos = IMPOSSIBLE_NEGATIVE;
			//p_query("ADDED", r);
			// For reads partially on left tail, the locus is negative
			if (r->status == FRESH || r->status == TRIED) {
				//if (strcmp(r->name, "15398") == 0)
				//p_query("ADDED", r);
				//p_query("TEMPL", window);
				ol = seq_ol(r, window, r->len, mismatches);
				//show_debug_msg(__func__, "Overlap: %d \n", ol);
				if (ol >= 0) {
					//p_query("NEW", r);
					add2tpl(t, r, i - left_len);
				} else {
					// For some read, its forward and reverse sequences are similar!
					r->rev_com = r->rev_com ? 0 : 1;
					if (seq_ol(r, window, r->len, mismatches) >= 0) {
						//p_query("NEW", r);
						add2tpl(t, r, i - left_len);
					} else
						r->rev_com = 0;
				}
			}
			if (r->contig_id == t->id)
				// If already on the template, just update the locus
				r->contig_locus = i - left_len;
		}
		g_ptr_array_free(hits, TRUE);
		bwa_free_read_seq(1, window);
		if (not_covered_len > ht->o->read_len * 2) {
			t->alive = 0;
			break;
		}
	}
	bwa_free_read_seq(1, seq);
	g_ptr_array_sort(t->reads, (GCompareFunc) cmp_reads_by_contig_locus);
	int n_pairs = count_pairs_on_tpl(t);
	t->pair_pc = (((float) n_pairs) * 2.0) / ((float) t->reads->len);
	t->cov = calc_tpl_cov(t, 0, t->len, ht->o->read_len);
}

/**
 * After connecting to existing, need to set the junction reads to USED,
 * 	to avoid multiple junctions
 */
void refresh_reads_on_tail(hash_table *ht, tpl *t, int mismatches) {
	bwa_seq_t *tail = NULL, *s = NULL, *window = NULL;
	bwa_seq_t *r = NULL;
	int i = 0, j = 0, len = 0, borrow_len = 0;
	GPtrArray *hits = g_ptr_array_sized_new(4);
	window = blank_seq(ht->o->read_len);
	if (t->l_tail) {
		borrow_len = (t->len >= (ht->o->read_len - 1)) ? (ht->o->read_len)
				: t->len;
		len = t->l_tail->len + borrow_len;
		s = blank_seq(len);
		memcpy(s->seq, t->l_tail->seq, t->l_tail->len * sizeof(ubyte_t));
		s->len += t->l_tail->len;
		memcpy(s->seq + s->len, t->ctg->seq, borrow_len * sizeof(ubyte_t));
		s->len += borrow_len;

		for (i = 0; i <= s->len - ht->o->read_len; i++) {
			copy_partial(s, window, i, ht->o->read_len);
			hits->len = 0;
			hits = find_both_fr_full_reads(ht, window, hits, mismatches);
			for (j = 0; j < hits->len; j++) {
				r = (bwa_seq_t*) g_ptr_array_index(hits, j);
				//if (strcmp(r->name, "15398") == 0)
				//			p_query(__func__, r);
				// For reads partially on left tail, the locus is negative
				if (r->status == FRESH)
					add2tpl(t, r, i - t->l_tail->len);
			}
		}
		bwa_free_read_seq(1, s);
	}

	if (t->r_tail) {
		borrow_len = (t->len >= (ht->o->read_len - 1)) ? (ht->o->read_len)
				: t->len;
		len = t->r_tail->len + borrow_len;
		s = blank_seq(len);
		memcpy(s->seq, t->ctg->seq + (t->len - borrow_len), borrow_len
				* sizeof(ubyte_t));
		s->len = borrow_len;
		memcpy(s->seq + s->len, t->r_tail->seq, t->r_tail->len
				* sizeof(ubyte_t));
		s->len += t->r_tail->len;

		for (i = 0; i <= s->len - ht->o->read_len; i++) {
			copy_partial(s, window, i, ht->o->read_len);
			hits->len = 0;
			hits = find_both_fr_full_reads(ht, window, hits, mismatches);
			for (j = 0; j < hits->len; j++) {
				r = (bwa_seq_t*) g_ptr_array_index(hits, j);
				// For reads partially on left tail, the locus is negative
				if (r->status == FRESH)
					add2tpl(t, r, t->len - borrow_len + i);
			}
		}
		bwa_free_read_seq(1, s);
	}
	bwa_free_read_seq(1, window);
	g_ptr_array_free(hits, TRUE);
}

/**
 * Assumption: the reads on the template are sorted by contig locus.
 * Correct template bases.
 */
void correct_tpl_base(bwa_seq_t *seqs, tpl *t, const int read_len, int start,
		int end) {
	int *cs = NULL, max = 0, weight = 0;
	GPtrArray *counters = NULL;
	int c = 0, max_c = 0;
	int i = 0, j = 0, locus = 0;
	bwa_seq_t *r = NULL, *m = NULL;
	if (!t->reads || t->reads->len == 0)
		return;
	start = max(0, start);
	end = min(t->len, end);
	counters = g_ptr_array_sized_new(t->len);
	//show_debug_msg(__func__,
	//		"Correcting template [%d, %d] at range [%d, %d] ...\n", t->id,
	//		t->len, start, end);
	//p_tpl_reads(t);
	//p_ctg_seq("BEFORE", t->ctg);

	for (i = 0; i < t->len; i++) {
		cs = (int*) calloc(5, sizeof(int));
		g_ptr_array_add(counters, cs);
	}

	// Count the occurrences of nucleotides at every position
	for (i = 0; i < t->reads->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(t->reads, i);
		m = get_mate(r, seqs);
		weight = 1;
		if (m->status == USED && m->contig_id == t->id)
			weight = MATE_MULTI;
		for (j = 0; j < r->len; j++) {
			locus = r->contig_locus + j;
			if (locus >= start && locus < end) {
				c = r->rev_com ? r->rseq[j] : r->seq[j];
				cs = (int*) g_ptr_array_index(counters, locus);
				cs[c] += weight;
			}
		}
	}

	for (i = 0; i < t->len; i++) {
		cs = (int*) g_ptr_array_index(counters, i);
		//show_debug_msg(__func__, "Locus %d: [%d,%d,%d,%d,%d]\n", i, cs[0], cs[1], cs[2], cs[3], cs[4]);
		max = -1;
		max_c = -1;
		for (j = 0; j < 4; j++) {
			if (cs[j] > 0 && cs[j] > max_c) {
				max = j;
				max_c = cs[j];
			}
		}
		if (max_c > 0) {
			if (t->ctg->seq[i] != max)
			show_debug_msg(__func__, "@%d Max: %c -> %c; Count: %d\n", i, "ACGTN"[t->ctg->seq[i]], "ACGTN"[max], max_c);
			t->ctg->seq[i] = max;
		}
	}
	set_rev_com(t->ctg);

	for (i = 0; i < t->len; i++) {
		cs = (int*) g_ptr_array_index(counters, i);
		free(cs);
	}
	g_ptr_array_free(counters, TRUE);
	//p_ctg_seq("AFTER", t->ctg);
}

void clear_tpl_tails(tpl *t) {
	if (t->r_tail) {
		bwa_free_read_seq(1, t->r_tail);
		t->r_tail = NULL;
	}
	if (t->l_tail) {
		bwa_free_read_seq(1, t->l_tail);
		t->l_tail = NULL;
	}
	if (t->b_juncs) {
		while (t->b_juncs->len > 0)
			g_ptr_array_remove_index_fast(t->b_juncs, 0);
	}
	if (t->m_juncs) {
		while (t->m_juncs->len > 0)
			g_ptr_array_remove_index_fast(t->m_juncs, 0);
	}
}

void reset_boundary_reads(tpl *t, const int ori) {
	int i = 0;
	bwa_seq_t *r = NULL;
	if (ori) {
		for (i = 0; i < t->reads->len; i++) {
			r = (bwa_seq_t*) g_ptr_array_index(t->reads, i);
			if (r->contig_locus > t->len - r->len)
				reset_to_fresh(r);
		}
	} else {
		for (i = 0; i < t->reads->len; i++) {
			r = (bwa_seq_t*) g_ptr_array_index(t->reads, i);
			if (r->contig_locus < 0)
				reset_to_fresh(r);
		}
	}
}

/**
 * Assumption: the input reads are in format of Forward-Forward
 * If the read should be at left side: -1
 * 						   right side: 1
 *                            unknown: 0
 * @r: the read whose mate is on the existing template
 */
int should_at_which_side(bwa_seq_t *seqs, tpl *t, bwa_seq_t *r) {
	bwa_seq_t *read = NULL, *mate = NULL, *m = NULL, *left = NULL, *right =
			NULL;
	int i = 0, side = UNKNOWN_SIDE;
	m = get_mate(r, seqs);
	if (m->contig_id != t->id)
		return UNKNOWN_SIDE;
	//p_tpl_reads(t);
	for (i = 0; i < t->reads->len; i++) {
		read = (bwa_seq_t*) g_ptr_array_index(t->reads, i);
		mate = get_mate(read, seqs);
		if (read->rev_com == mate->rev_com && mate->contig_id
				== read->contig_id && mate->contig_id == t->id) {
			left = read->contig_locus < mate->contig_locus ? read : mate;
			right = read->contig_locus < mate->contig_locus ? mate : read;

			//p_query("LEFT", left);
			//p_query("RIGHT", right);
			//p_query("JUMPED", r);
			//p_query("MATE", m);
			// Existing pair and the jumped pair are with same orientation
			if (is_left_mate(left->name)) {
				if (left->rev_com == m->rev_com) {
					if (is_left_mate(m->name))
						side = RIGHT_SIDE;
					else
						side = LEFT_SIDE;
				} else {
					if (is_left_mate(m->name))
						side = LEFT_SIDE;
					else
						side = RIGHT_SIDE;
				}
			} else {
				if (left->rev_com == m->rev_com) {
					if (is_left_mate(m->name))
						side = LEFT_SIDE;
					else
						side = RIGHT_SIDE;
				} else {
					if (is_left_mate(m->name))
						side = RIGHT_SIDE;
					else
						side = LEFT_SIDE;
				}
			}
			break;
		}
	}
	return side;
}

/**
 * Remove a read from the pool and reset the the read status
 */
void rm_from_tpl(tpl *t, int index) {
	if (index < 0 || !t || !t->reads || index >= t->reads->len)
		return;
	bwa_seq_t * r = (bwa_seq_t*) g_ptr_array_index(t->reads, index);
	reset_to_fresh(r);
	g_ptr_array_remove_index_fast(t->reads, index);
}

/**
 * From the tried, remove the reads with contig id as rm_id
 */
void rm_from_tried(tpl *t, const int rm_id) {
	int i = 0;
	bwa_seq_t *r = NULL;
	for (i = 0; i < t->tried->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(t->tried, i);
		if (r->contig_id == rm_id) {
			g_ptr_array_remove_index_fast(t->tried, i);
			break;
		}
	}
}

/**
 * If some read is already used by some template, do not reset
 */
void unhold_reads_array(GPtrArray *reads) {
	int i = 0;
	bwa_seq_t *r = NULL;
	for (i = 0; i < reads->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(reads, i);
		if (r->status == FRESH)
			reset_to_fresh(r);
	}
	g_ptr_array_free(reads, TRUE);
}

/**
 * Read extension starts from some read, mark this read and its close reads as used first
 */
void mark_init_reads_used(hash_table *ht, tpl *t, bwa_seq_t *read,
		int mismatches) {
	GPtrArray *hits = NULL;
	bwa_seq_t *r = NULL;
	index64 i = 0;

	if (read->rev_com)
		switch_fr(read);

	hits = find_both_fr_full_reads(ht, read, hits, mismatches);
	for (i = 0; i < hits->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(hits, i);
		//p_query(__func__, r);
		if (r->status == FRESH) {
			add2tpl(t, r, 0);
		}
	}
	if (read->rev_com)
		switch_fr(read);
	//p_readarray(t->reads, 1);
	g_ptr_array_free(hits, TRUE);
}

/**
 * Align tail-length query to find reads for branching
 *
 * At a position, align the 25bp portion to find FRESH reads overlapped with it
 * These reads should be partially overlapped with the template.
 *
 * For example, ori 0:
 * Template: ========----------===========
 * Align ---------- and find a read **----------**********
 * Then the read should overlap with the template for the portion **----------
 */
GPtrArray *check_branch_tail(hash_table *ht, tpl *t, bwa_seq_t *query,
		int shift, int mismatches, int8_t status, int ori) {
	bwa_seq_t *r = NULL, *tpl_seq = NULL;
	GPtrArray *hits = align_query(ht, query, status, mismatches);
	GPtrArray *picked = NULL;
	int ol_len = 0, n_mis = 0, start = 0;
	int i = 0, j = 0, x = 0;
	int is_picked = 0;
	ubyte_t read_c = 0, ref_c = 0;

	if (hits->len <= 0) {
		return hits;
	}
	tpl_seq = blank_seq(ht->o->read_len);
	picked = g_ptr_array_sized_new(hits->len);
	for (i = 0; i < hits->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(hits, i);
		if (has_n(r, 1)) {
			reset_to_fresh(r); continue;
		}

		ol_len = ori ? r->len - r->pos : query->len + r->pos;
		start = ori ? shift : shift - r->pos;
		//if (ol_len < query->len - 2)
		//	continue;
		//show_debug_msg(__func__, "Start: %d; ol: %d\n", start, ol_len);

		// Make sure that the mapped portion has acceptable mismatches
		is_picked = 0;
		if (start >= 0 && start + ol_len <= t->len) {
			copy_partial(t->ctg, tpl_seq, start, ol_len);
			//if (shift == 1076) {
			//p_query("TEMPLATE SEQ", tpl_seq);
			//p_query("HIT", r);
			//}
			if (ori) n_mis = seq_ol(r, tpl_seq, ol_len, mismatches);
			else n_mis = seq_ol(tpl_seq, r, ol_len, mismatches);
			//if (shift == 1076)
			//show_debug_msg(__func__, "n_mis: %d \n", n_mis);

			// Determine the cursor position, if the cursor is
			//	at head/tail (N_BAD_TAIL_SHIFT bp), not pick it
			if (n_mis >= 0) {
				if (ori) {
					for (j = r->pos - 1; j >= N_BAD_TAIL_SHIFT; j--) {
						read_c = r->rev_com ? r->rseq[j] : r->seq[j];
						if ((shift - (r->pos - j)) < 0
								|| t->ctg->seq[shift - (r->pos - j)] != read_c) {
							r->cursor = j;
							g_ptr_array_add(picked, r);
							is_picked = 1;
							break;
						}
					}
				} else {
					for (j = ol_len; j < r->len - N_BAD_TAIL_SHIFT; j++) {
						read_c = r->rev_com ? r->rseq[j] : r->seq[j];
						if ((start + j) >= t->len || t->ctg->seq[start + j] != read_c) {
							r->cursor = j;
							g_ptr_array_add(picked, r);
							is_picked = 1;
							break;
						}
					}
				}
			}
		}
		if (!is_picked) {
			reset_to_fresh(r);
		}
	}
	bwa_free_read_seq(1, tpl_seq);
	g_ptr_array_free(hits, TRUE);
	if (ori) g_ptr_array_sort(picked, (GCompareFunc) cmp_reads_by_cursor);
	else g_ptr_array_sort(picked, (GCompareFunc) cmp_reads_by_rev_cursor);
	return picked;
}

int has_nearby_pairs(hash_table *ht, GPtrArray *tpls, tpl *t, int n_pairs) {
	if (!t->reads || t->reads->len == 0 || n_pairs <= 0)
		return 0;
	int n = 0;
	bwa_seq_t *r = NULL, *m = NULL;
	int i = 0, j = 0;
	tpl *near = NULL;
	show_debug_msg(__func__, "Checking pairs for template [%d, %d] \n", t->id,
			t->len);
	for (i = 1; i < t->reads->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(t->reads, i);
		if (r->contig_locus < 0 || r->contig_locus > t->len - r->len) continue;
		m = get_mate(r, ht->seqs);
		if (m->status == USED) {
			for (j = 0; j < tpls->len; j++) {
				near = (tpl*) g_ptr_array_index(tpls, j);
				if (m->contig_id == near->id && m->contig_id != t->id ) {
					if (m->contig_locus < 0 || m->contig_locus > near->len - m->len) break;
					n++;
					p_query("READ", r);
					p_query("MATE", m);
					break;
				}
			}
		}
		if (n >= n_pairs)
			return 1;
	}
	return 0;
}

/**
 * Align the tail of template to find fresh reads.
 * Illustration:
 *
 * Contig read:    ==========		read_len: 10
 * Tail shift:         ^(shift=4)
 * Aligned tail:       ----			tail_len: 4
 *                     ||||
 * A hit read:        ==========	read_len: 10
 * Pos:                ^(pos=1)
 * Checked ol:        =======		read_len - (shift - pos) = 7
 */
GPtrArray *align_tpl_tail(hash_table *ht, tpl *t, bwa_seq_t *tail, int limit,
		int shift, int mismatches, int8_t status, int ori) {
	int64_t i = 0, cursor = 0, ol = 0, locus = 0;
	int n_mis = 0, added = 0;
	bwa_seq_t *r = NULL, *tpl_seq = NULL, *ol_seq = NULL;
	GPtrArray *hits = align_query(ht, tail, status, mismatches);
	GPtrArray *fresh_reads = g_ptr_array_sized_new(hits->len);

	tpl_seq = get_tail(t, ht->o->read_len, ori);
	if (t->len == -1)
		p_ctg_seq(__func__, tpl_seq);

	// These reads are not duplicated
	for (i = 0; i < hits->len; i++) {
		// If max <= 0, do not limit the reads to return
		//if (limit > 0 && fresh_reads->len >= limit)
		//	break;

		r = (bwa_seq_t*) g_ptr_array_index(hits, i);
		added = 0;
		// pos is the kmer position on the read
		cursor = ori ? (r->pos - 1) : (r->pos + tail->len);

		if (t->id == 108220) {
			printf("\n");
			p_query(__func__, tail);
			p_query(__func__, r);
			show_debug_msg(__func__, "CURSOR: %d\n", cursor);
		}

		if (ori) {
			cursor = r->pos - shift;
			ol = tpl_seq->len - cursor;
			locus = tpl_seq->len - ol;
			cursor--;
		} else {
			locus = shift - r->pos;
			ol = tpl_seq->len - locus;
			cursor = ol;
		}
		if (t->len == -1) {
			//	p_query(__func__, r);
			show_debug_msg(__func__, "LOCUS: %d; OL: %d; CURSOR: %d\n", locus,
					ol, cursor);
		}
		if (ol >= tail->len && locus >= 0 && locus < tpl_seq->len && (cursor
				>= 0 && cursor <= r->len - 1)) {
			//ol_seq = ori ? new_seq(tpl_seq, ol, 0)
			//		: new_seq(tpl_seq, ol, locus);
			if (ori) {
				n_mis = seq_ol(r, tpl_seq, ol, mismatches);
			} else {
				n_mis = seq_ol(tpl_seq, r, ol, mismatches);
			}
			if (t->len == -1 && t->id == -1) {
				p_query("CONTIG", tpl_seq);
				p_query("READ  ", r);
				p_ctg_seq("OVERLP", ol_seq);
				show_debug_msg(__func__, "CURSOR: %d\n", cursor);
				show_debug_msg(__func__, "Should have overlap: %d\n", ol);
				show_debug_msg(__func__, "N_MISMATCHES: %d\n --- \n", n_mis);
			}
			// n_mis >= 0 means similar with n_mis mismatches; -1 means not similar
			if (n_mis >= 0) {
				//show_debug_msg(__func__, "Cursor: %d\n", cursor);
				r->cursor = cursor;
				// In the pool, 'pos' stores how many mismatches between read and template
				r->pos = n_mis;
				//p_query("ADDED", r);
				g_ptr_array_add(fresh_reads, r);
				added = 1;
			} else {
				// For some read, its forward and reverse complement are highly similar!
				r->rev_com = r->rev_com ? 0 : 1;
				if (ori) {
					n_mis = seq_ol(r, tpl_seq, ol, mismatches);
				} else {
					n_mis = seq_ol(tpl_seq, r, ol, mismatches);
				}
				if (n_mis >= 0) {
					r->cursor = cursor;
					r->pos = n_mis;
					g_ptr_array_add(fresh_reads, r);
					added = 1;
				}
			}
			//bwa_free_read_seq(1, ol_seq);
		}
		// If not added, reset the pos to be -1
		if (!added) {
			r->pos = IMPOSSIBLE_NEGATIVE;
			r->rev_com = 0;
		}

		if (t->id == 108220)
			p_query(__func__, r);
	}
	//show_debug_msg(__func__, "Reads with the tail: %d\n", fresh_reads->len);
	g_ptr_array_free(hits, TRUE);
	bwa_free_read_seq(1, tpl_seq);
	//	if (t->len == 153 && t->id == 3) {
	//		p_readarray(fresh_reads, 1);
	//	}
	return fresh_reads;
}

/**
 * Reverse complement a template, need to update the read locus as well
 */
void switch_tpl_fr(tpl *t) {
	int i = 0;
	bwa_seq_t *r = NULL, *tmp = NULL;
	switch_fr(t->ctg);
	for (i = 0; i < t->reads->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(t->reads, i);
		r->contig_locus = t->len - r->contig_locus - r->len;
		r->rev_com = r->rev_com ? 0 : 1;
	}
	//	if (t->r_tail)
	//		tmp = t->r_tail;
	//	if (t->l_tail) {
	//		switch_fr(t->l_tail);
	//		t->r_tail = t->l_tail;
	//	}
	//	if (tmp) {
	//		switch_fr(tmp);
	//		t->l_tail = tmp;
	//	}
}

void save_tpls(tplarray *pfd_ctg_ids, FILE *ass_fa, const int ori,
		const int p_all, const int min_len) {
	int i = 0;
	tpl *t;
	char *h;
	bwa_seq_t *contig;
	if (!pfd_ctg_ids || pfd_ctg_ids->len <= 0)
		return;
	h = (char*) malloc(BUFSIZE);
	for (i = 0; i < pfd_ctg_ids->len; i++) {
		t = (tpl*) g_ptr_array_index(pfd_ctg_ids, i);
		//show_debug_msg(__func__, "Saving tpl %d length %d, alive %d \n", t->id, t->len, t->alive);
		if (p_all || (t && t->alive && t->ctg && t->len > min_len)) {
			contig = t->ctg;
			if (ori)
				seq_reverse(contig->len, contig->seq, 0);
			sprintf(h, ">%d\tlength: %d; start: %s; step: %s; pair_percentage: %.2f\n", t->id, contig->len,
						t->start_read->name, t->step, t->pair_pc);
			save_con(h, contig, ass_fa);
		}
	}
	free(h);
}

void destroy_tpl(tpl *t, int status) {
	index64 i = 0;
	bwa_seq_t *r = NULL;
	if (t) {
		show_debug_msg(__func__, "Freeing tpl [%d, %d] \n", t->id, t->len);
		if (t->ctg) bwa_free_read_seq(1, t->ctg);
		if (t->r_tail) bwa_free_read_seq(1, t->r_tail);
		if (t->l_tail) bwa_free_read_seq(1, t->l_tail);
		if (t->vertexes) g_ptr_array_free(t->vertexes, TRUE);
		if (t->reads) {
			for (i = 0; i < t->reads->len; i++) {
				r = (bwa_seq_t*) g_ptr_array_index(t->reads, i);
				r->contig_id = -1;
				r->contig_locus = -1;
				r->pos = IMPOSSIBLE_NEGATIVE;
				r->status = status;
				//reset_to_fresh(r);
			}
			g_ptr_array_free(t->reads, TRUE);
		}

		if (t->tried) {
			for (i = 0; i < t->tried->len; i++) {
				r = (bwa_seq_t*) g_ptr_array_index(t->tried, i);
				r->contig_id = -1;
				r->contig_locus = -1;
				r->pos = IMPOSSIBLE_NEGATIVE;
				r->status = status;
				//reset_to_fresh(r);
			}
			g_ptr_array_free(t->tried, TRUE);
		}

		if (t->m_juncs) g_ptr_array_free(t->m_juncs, TRUE);
		if (t->b_juncs) g_ptr_array_free(t->b_juncs, TRUE);
		if (t->step) free(t->step);
		free(t);
	}
}

/**
 * Keep only the template sequence
 */
void keep_ctg_only(tpl *t) {
	if (t) {
		if (t->r_tail) {
			bwa_free_read_seq(1, t->r_tail);
			t->r_tail = NULL;
		}
		if (t->l_tail) {
			bwa_free_read_seq(1, t->l_tail);
			t->l_tail = NULL;
		}
		if (t->vertexes) {
			g_ptr_array_free(t->vertexes, TRUE);
			t->vertexes = NULL;
		}
		if (t->reads) {
			g_ptr_array_free(t->reads, TRUE);
			t->reads = NULL;
		}
		if (t->tried) {
			g_ptr_array_free(t->tried, TRUE);
			t->tried = NULL;
		}
		if (t->m_juncs) {
			g_ptr_array_free(t->m_juncs, TRUE);
			t->m_juncs = NULL;
		}
		if (t->b_juncs) {
			g_ptr_array_free(t->b_juncs, TRUE);
			t->b_juncs = NULL;
		}
	}
}

/**
 * Get the sequence with left/right tails.
 */
bwa_seq_t *get_tpl_ctg_wt(tpl *t, int *l_len, int *r_len, int *t_len) {
	bwa_seq_t *s = NULL;
	int counted_len = 0;
	*l_len = (t->l_tail) ? (t->l_tail->len) : 0;
	*r_len = (t->r_tail) ? (t->r_tail->len) : 0;
	*t_len = *l_len + t->len + *r_len;
	s = blank_seq(*t_len);
	if (t->l_tail) {
		memcpy(s->seq, t->l_tail->seq, t->l_tail->len * sizeof(ubyte_t));
		s->len = *l_len;
	}
	memcpy(s->seq + s->len, t->ctg->seq, t->len);
	s->len += t->len;
	if (t->r_tail) {
		memcpy(s->seq + s->len, t->r_tail->seq, t->r_tail->len
				* sizeof(ubyte_t));
		s->len += t->r_tail->len;
	}
	set_rev_com(s);
	return s;
}

/**
 * Get the reads supporting the region on the template.
 * The caller need to free the memory of the return array
 */
GPtrArray *get_supporting_reads(tpl *t, int start, int end) {
	GPtrArray *reads = g_ptr_array_sized_new(4);
	int i = 0, rstart = 0, rend = 0;
	bwa_seq_t *r = NULL;
	for (i = 0; i < t->reads->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(t->reads, i);
		rend = min(t->len, r->contig_locus + r->len);
		rstart = max(0, r->contig_locus);
		if (start >= rstart && start < rend)
			g_ptr_array_add(reads, r);
		else if (end >= rend && end < rend)
			g_ptr_array_add(reads, r);
	}
	return reads;
}

int pairs_spanning_locus(bwa_seq_t *seqs, tpl *t, int locus) {
	bwa_seq_t *r = NULL, *m = NULL;
	int i = 0, n_pairs = 0;
	for (i = 0; i < t->reads->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(t->reads, i);
		m = get_mate(r, seqs);
		if (m->status == USED && m->contig_id == t->id) {
			if ((r->contig_locus >= locus && m->contig_locus) ||
				(r->contig_locus + r->len <= locus && m->contig_locus + m->len <= locus))
				continue;
			n_pairs++;
		}
	}
	return n_pairs;
}

/**
 * Count reads covering 'locus'; the reads should be within range between [s, e)
 */
int reads_covering_locus(tpl *t, int s, int locus, int e) {
	bwa_seq_t *r = NULL;
	int i = 0, n = 0;
	for (i = 0; i < t->reads->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(t->reads, i);
		if (r->contig_locus < s || r->contig_locus + r->len >= e) continue;
		if (r->contig_locus < locus && r->contig_locus + r->len > locus) n++;
	}
	return n;
}

int read_on_tpl(tpl *t, bwa_seq_t *r) {
	if (r->status == USED && r->contig_id == t->id) return 1;
	return 0;
}

int count_pairs_on_tpl(tpl *t) {
	int n = 0, i = 0;
	g_ptr_array_sort(t->reads, (GCompareFunc) cmp_reads_by_name);
	bwa_seq_t *r = NULL, *pre = NULL;
	if (t->reads->len < 2) return 0;
	pre = (bwa_seq_t*) g_ptr_array_index(t->reads, 0);
	for (i = 1; i < t->reads->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(t->reads, i);
		if (is_a_pair(pre, r)) n++;
		pre = r;
	}
	show_debug_msg(__func__, "Template [%d, %d]: %d / %d \n", t->id, t->len, n, t->reads->len);
	return n;
}

void destory_tpl_ht(hash_table *ht) {
	if (!ht) return;
	if (ht->k_mers_occ_acc) free(ht->k_mers_occ_acc);
	if (ht->n_kmers) free(ht->n_kmers);
	if (ht->o) free(ht->o);
	if (ht->pos) free(ht->pos);
	free(ht);
}

hash_table *hash_tpls(GPtrArray *tpls, int k, int interleaving) {
	hash_table *ht = (hash_table*) malloc(sizeof(hash_table));
	hash_key *k_mers_occ_acc, key = 0ULL;
	hash_value *pos = NULL;
	index64 n_pos = 0, pos_index = 0, n_seqs = tpls->len;
	hash_value value;
	int *n_occ = NULL;
	int tmp = 0, tmp_2 = 0;
	index64 n_k_mers = (1 << (k * 2)) + 1;
	tpl *t = NULL;
	int i = 0, j = 0, i_acc = 0;
	uint32_t *kmer_occ_on_reads = NULL;

	hash_opt *opt = init_hash_opt();
	opt->k = k;
	opt->interleaving = interleaving;
	opt->read_len = -1;
	opt->block_size = 1;

	show_msg(__func__, "Hashing %d templates ... \n", n_seqs);
	ht->seqs = NULL;
	ht->n_seqs = n_seqs;
	ht->o = opt;
	k_mers_occ_acc = (hash_value*) calloc(n_k_mers, sizeof(hash_value));

	show_debug_msg(__func__, "Round 1/2: Counting occurrences of k-mers %d templates ... \n", n_seqs);
	for (i = 0; i < n_seqs; i++) {
		t = (tpl*) g_ptr_array_index(tpls, i);
		if (t->len > 65536) {
			show_debug_msg(__func__, "[WARNING] Template [%d, %d] is too long to be hashed. \n", t->id, t->len);
			continue;
		}
		if (t->len < k) continue;
		for (j = 0; j <= t->ctg->len - k; j++) {
			key = get_hash_key(t->ctg->seq, j, interleaving, k);
			k_mers_occ_acc[key]++;
			n_pos++;
		}
	}

	tmp = 0;
	for (i = 1; i < n_k_mers - 1; i++) {
		tmp_2 = k_mers_occ_acc[i];
		k_mers_occ_acc[i] = k_mers_occ_acc[i - 1] + tmp;
		tmp = tmp_2;
	}
	k_mers_occ_acc[0] = 0;
	k_mers_occ_acc[n_k_mers - 1] = n_pos - 1;

	// Round 2: save the hashes.
	show_debug_msg(__func__, "Round 2/2: Placing k-mer pointers ... \n");
	pos = (hash_value*) calloc(n_pos, sizeof(hash_value));
	n_occ = (int*) calloc(n_k_mers, sizeof(int)); // Store how many occs through the iteration

	for (i = 0; i < n_seqs; i++) {
		t = (tpl*) g_ptr_array_index(tpls, i);
		if (t->len > 65536) {
			show_debug_msg(__func__, "[WARNING] Template [%d, %d] is too long to be hashed. \n", t->id, t->len);
			continue;
		}
		if (t->len < k) continue;
		for (j = 0; j <= t->ctg->len - k; j++) {
			key = get_hash_key(t->ctg->seq, j, interleaving, k);
			value = get_hash_value(t->id, j);
			pos_index = k_mers_occ_acc[key] + n_occ[key];
			pos[pos_index] = value;
			n_occ[key]++;
		}
	}
	opt->n_k_mers = n_k_mers;
	opt->n_pos = n_pos;
	ht->pos = pos;
	ht->k_mers_occ_acc = k_mers_occ_acc;
	ht->n_kmers = NULL;
	free(n_occ);
	show_debug_msg(__func__, "%d templates hashed. \n", n_seqs);
	return ht;
}

