#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <glib.h>
#include <math.h>
#include <pthread.h>
#include "rand.h"
#include "mate.h"
#include "clean.h"
#include "bwtaln.h"
#include "utils.h"
#include "peseq.h"
#include "pool.h"
#include "pechar.h"
#include "pealn.h"
#include "roadmap.h"
#include "edge.h"
#include "edgelist.h"
#include "readrm.h"
#include "pelib.h"
#include "pepath.h"
#include "scaffolding.h"
#include "pehash.h"
#include "correct.h"
#include "merge.h"
#include "hits.h"
#include "rnaseq.h"
#include "kmers.h"

int pair_ctg_id = 0;
int overlap_len = 0;
int insert_size = 0;
int sd_insert_size = 0;
char *out_root = NULL;
int n_threads = 1;
int stage = 0;
char *blat_exe = NULL;
int debug_mode = 0;
GMutex *id_mutex;
GMutex *sum_mutex;
struct timespec start_time, finish_time;

void pe_lib_help() {
	show_msg(__func__, "--------------------------------------------------");
	show_msg(__func__, "peta pelib <library>");
	show_msg(__func__, "--------------------------------------------------");
}

char *get_output_file(const char *file_name) {
	char *name = (char*) calloc(strlen(file_name) + strlen(out_root) + 4,
			sizeof(char));
	strcat(name, out_root);
	strcat(name, "/");
	strcat(name, file_name);
	return name;
}

void concat_doubles(double *base, int *n_total, double *partial, int n_part) {
	int i = 0;
	for (i = 0; i < n_part; i++) {
		base[i + *n_total] = partial[i];
		//		show_msg(__func__, "%d: %f \n", *n_total, partial[i]);
	}
	*n_total += n_part;
}

/**
 * Keep only those reads with its mates at the right direction and position
 *
 * ori = 0: For below, 'read' should be right, 'mate' should be left
 *   ---------------------->
 *          ---->        ---->
 *          mate         read
 */
void keep_mates_in_pool(edge *eg, pool *cur_pool, const hash_table *ht,
		const int ori, const int remove_single) {
	bwa_seq_t *read = NULL, *mate = NULL;
	int i = 0, to_remove = 1;
	for (i = 0; i < cur_pool->reads->len; i++) {
		read = g_ptr_array_index(cur_pool->reads, i);
		mate = get_mate(read, ht->seqs);
		//		// When extending to left, the initial mate pool and current pool overlap
		//		if (read->is_in_m_pool == eg->tid)
		//			continue;
		to_remove = 0;
		//p_query(__func__, mate);
		// The mate of 'read' should be used already
		//if (strcmp(read->name, "5896388") == 0) {
		//	show_debug_msg(__func__, "ORI: %d \n", ori);
		//	show_debug_msg(__func__, "Edge [%d, %d] \n", eg->id, eg->len);
		//	p_query(__func__, mate);
		//	p_query(__func__, read);
		//}
		if (is_paired(read, ori)) {
			// If the mate is not used, remove 'read'
			if ((read->rev_com != mate->rev_com) || (mate->status != TRIED)
					|| eg->tid != mate->tid) {
				//p_query(__func__, read);
				//p_query(__func__, mate);
				to_remove = 1;
			} else { // If the distance is not within range, remove 'read'
				if (get_abs(eg->len - mate->shift) > (insert_size
						+ sd_insert_size * SD_TIMES)) {
					to_remove = 1;
				}
			}
		} else {
			to_remove = remove_single;
		}

		// If the read is obtained by another thread, just remove it
		if (to_remove || read->is_in_c_pool != eg->tid) {
			read->status = TRIED;
			read->contig_id = eg->id;
			if (pool_rm_index(cur_pool, i)) {
				read->shift = eg->len - read->cursor + 1;
				i--;
			}
		}
	}
}

/**
 * Count the occurrences of next probable chars.
 * Result: set the array 'next' as the counters.
 */
void check_next_char(pool *cur_pool, edge *eg, int *next, const int ori) {
	int i = 0, pre_pos = 0, check_pre = 0;
	bwa_seq_t *s = NULL;
	if (cur_pool->n > 10)
		check_pre = 1;
	for (i = 0; i < cur_pool->n; i++) {
		s = g_ptr_array_index(cur_pool->reads, i);
		pre_pos = ori ? s->cursor + 1 : s->cursor - 1;
		// Only if read's last character is some as the contig's character, count it.
		if (s->rev_com) {
			if (check_pre) {
				if (pre_pos >= 0 && pre_pos < s->len && s->rseq[pre_pos]
						== eg->contig->seq[eg->contig->len - 1])
					check_c(next, s->rseq[s->cursor]);
			} else {
				check_c(next, s->rseq[s->cursor]);
			}
		} else {
			if (check_pre) {
				if (pre_pos >= 0 && pre_pos < s->len && s->seq[pre_pos]
						== eg->contig->seq[eg->contig->len - 1])
					check_c(next, s->seq[s->cursor]);
			} else {
				check_c(next, s->seq[s->cursor]);
			}
		}
	}
}

/**
 * From the alignment results, add/remove reads to/from the current pool.
 */
void maintain_pool(alignarray *aligns, const hash_table *ht, pool *cur_pool,
		edge *eg, bwa_seq_t *query, int *next, int ori) {
	alg *a = NULL;
	int i = 0, index = 0;
	bwa_seq_t *s = NULL, *mate = NULL, *seqs = ht->seqs;
	//show_debug_msg(__func__, "Iterating alignments... \n");
	// Add aligned reads to current pool
	for (i = 0; i < aligns->len; i++) {
		a = g_ptr_array_index(aligns, i);
		index = a->r_id;
		if (index >= ht->n_seqs)
			continue;
		s = &seqs[index];
		mate = get_mate(s, seqs);
		//p_query(__func__, s);
		//p_query(__func__, mate);
		//pre_cursor = s->cursor;

		if (s->status == USED || s->status == DEAD || (s->status == TRIED
				&& s->contig_id == eg->id) || s->tid != -1)
			continue;

		//          read              mate              read
		//        =========>        --------->        --------->
		// Edge: -------------------------------------------->
		// If the mate of read is TRIED by current edge,
		//    but from the mate, the read should be at the left, skip this read
		if (mate->status == TRIED && mate->contig_id == eg->id) {
			if (is_paired(mate, ori))
				continue;
		}
		s->rev_com = a->rev_comp;
		mate->rev_com = s->rev_com;
		if (s->rev_com)
			s->cursor = ori ? (s->len - query->len - 1 - a->pos) : (s->len
					- a->pos);
		else
			s->cursor = ori ? (a->pos - 1) : (a->pos + query->len);

		if (s->cursor >= s->len || s->cursor < 0) {
			s->cursor = 0;
			continue;
		}
		//if ((!ori && pre_cursor > s->cursor) || (ori && pre_cursor < s->cursor))
		//	continue;
		pool_add(cur_pool, s, eg->tid);
	}
	//show_debug_msg(__func__, "Removing partial... \n");
	// In current pool, if a read does not overlap with the tail properly, remove it
	//p_pool("BEFORE REMOVING PARTIAL", cur_pool, NULL);
	rm_partial(eg, cur_pool, ori, seqs, query, MISMATCHES);
	//p_pool("AFTER REMOVING PARTIAL", cur_pool, NULL);
	// Keep only the reads whose mate is used previously.
	if (eg->len >= (insert_size + sd_insert_size * SD_TIMES)) {
		//p_pool("BEFORE KEEPING MATES", cur_pool, NULL);
		keep_mates_in_pool(eg, cur_pool, ht, ori, 0);
		//p_pool("AFTER KEEPING MATES", cur_pool, NULL);
	}
	check_next_char(cur_pool, eg, next, ori);
}

/**
 * If return NULL, it indicates there are too many (>=4 by default) counter pairs.
 * In this case, the caller should not continue to do single extension.
 */
pool *get_start_pool(const hash_table *ht, bwa_seq_t *init_read, const int ori,
		const int tid) {
	alignarray *alns = NULL;
	int i = 0, to_free_query = 0;
	// In some cases, one mate is forward, but another is backward, which is not reasonable.
	// Here we introduce some toleration.
	int n_counter_pairs = 0;
	pool *init_pool = NULL;
	bwa_seq_t *seqs = ht->seqs, *s = NULL, *query = init_read, *mate = NULL;
	alg *a;
	if (is_repetitive_q(init_read) || has_rep_pattern(init_read)) {
		return NULL;
	}
	init_pool = new_pool();
	alns = g_ptr_array_sized_new(N_DEFAULT_ALIGNS);
	// If the initial read is reverse complement, make the query sequence accordingly
	if (init_read->rev_com) {
		query = new_mem_rev_seq(init_read, init_read->len, 0);
		to_free_query = 1;
	}
	//p_query(__func__, query);

	pe_aln_query(query, query->seq, ht, 0, overlap_len, 0, alns);
	pe_aln_query(query, query->rseq, ht, 0, overlap_len, 1, alns);
	//p_align(alns);
	for (i = 0; i < alns->len; i++) {
		a = (alg*) g_ptr_array_index(alns, i);
		s = &seqs[a->r_id];
		mate = get_mate(s, seqs);
		// If the read's mate is used, but the orientation is not the same, not add to pool
		if (mate->status == USED && (mate->rev_com != a->rev_comp)) {
			n_counter_pairs++;
			continue;
		}
		// If the pair of the initial read has wrong orientation, return NULL
		if (n_counter_pairs >= MAX_COUNTER_PAIRS) {
			break;
		}
		s->rev_com = a->rev_comp;
		if (s->rev_com)
			s->cursor = ori ? (0 - a->pos - 1) : (s->len - a->pos);
		else
			s->cursor = ori ? (a->pos - 1) : (a->pos + s->len);
		// p_query(__func__, s);
		if (s->contig_id == INVALID_CONTIG_ID || s->status != FRESH
				|| s->cursor < -1 || s->cursor > s->len || s->is_in_c_pool
				|| s->is_in_m_pool) {
			s->cursor = 0;
			continue;
		}
		pool_add(init_pool, s, tid);
	}
	free_alg(alns);
	if (to_free_query)
		bwa_free_read_seq(1, query);
	// Indicates that too many counter pairs found, not feasible to extend from init_read
	if (n_counter_pairs >= MAX_COUNTER_PAIRS) {
		free_pool(init_pool);
		return NULL;
	} else {
		//p_pool("INITIAL", init_pool, NULL);
		return init_pool;
	}
}

/**
 * Extend a template
 */
edge *pair_extension(edge *pre_eg, const hash_table *ht, bwa_seq_t *s,
		const int ori, const int tid) {
	bwa_seq_t *query = NULL;
	pool *cur_pool = NULL;
	int c = 0;
	int *next = NULL;
	alignarray *aligns = NULL;
	edge *eg = NULL;
	char *ctg_name;

	if (pre_eg) {
		eg = pre_eg;
	} else {
		eg = new_eg();
		eg->tid = tid;
		g_mutex_lock(id_mutex);
		eg->id = pair_ctg_id;
		pair_ctg_id++;
		show_msg(__func__, "---------- Processing read %d: %s ----------\n",
				pair_ctg_id, s->name);
		show_debug_msg(__func__,
				"---------- Processing read %d: %s ----------\n", pair_ctg_id,
				s->name);
		g_mutex_unlock(id_mutex);
		eg->name = strdup(s->name);
		eg->contig = new_seq(s, s->len, 0);
		ctg_name = (char*) calloc(16, sizeof(char));
		sprintf(ctg_name, "%d", eg->id);
		eg->contig->name = ctg_name;
		eg->len = s->len;
	}
	cur_pool = get_start_pool(ht, s, ori, eg->tid);
	if (!cur_pool)
		return eg;

	if (ori)
		query = new_seq(s, overlap_len, 0);
	else
		query = new_seq(s, overlap_len, s->len - overlap_len);
	aligns = g_ptr_array_sized_new(N_DEFAULT_ALIGNS);
	next = (int*) calloc(5, sizeof(int));
	if (ori)
		seq_reverse(eg->len, eg->contig->seq, 0);
	while (1) {
		p_query(__func__, query);
		reset_c(next, NULL); // Reset the counter
		if (is_repetitive_q(query) || has_rep_pattern(query)) {
			show_msg(__func__, "[%d, %d] Repetitive pattern, stop!\n", eg->id,
					eg->len);
			p_query("Repetitive pattern", query);
			// Try to add more reliable reads: the mate reads into the current pool
			if (eg->len > insert_size - SD_TIMES * sd_insert_size) {
				clear_pool(cur_pool);
				add_mates_by_ol(ht, eg, cur_pool, query->len - 4, 0, query,
						ori, insert_size, sd_insert_size);
				reset_c(next, NULL); // Reset the counter
				check_next_char(cur_pool, eg, next, ori);
			} else
				break;
		} else {
			pe_aln_query(query, query->seq, ht, MISMATCHES, query->len, 0,
					aligns);
			pe_aln_query(query, query->rseq, ht, MISMATCHES, query->len, 1,
					aligns);
			p_align(aligns);
			maintain_pool(aligns, ht, cur_pool, eg, query, next, ori);
			// If no read in current pool, try less stringent overlapped length and zero mismatches
			if (cur_pool->reads->len <= 0 && eg->len > insert_size - SD_TIMES
					* sd_insert_size) {
				//show_debug_msg(__func__, "Trying mate pool... \n");
				add_mates_by_ol(ht, eg, cur_pool, RELAX_MATE_OL_THRE,
						0, query, ori, insert_size, sd_insert_size);
				reset_c(next, NULL); // Reset the counter
				check_next_char(cur_pool, eg, next, ori);
				//p_pool("After adding mates", cur_pool, next);
			}
			reset_alg(aligns);
		}
		c = get_pure_most(next);
		show_debug_msg(__func__, "Ori %d, Edge %d, length %d \n", ori, eg->id,
				eg->len);
		p_ctg_seq("Contig", eg->contig);
		p_pool("Current Pool", cur_pool, next);
		//show_debug_msg(__func__, "Ori: %d, Next char: %d \n", ori, c);
		if (cur_pool->n <= 0) {
			show_msg(__func__, "[%d, %d] No hits, stop here. \n", eg->id,
					eg->len);
			show_debug_msg(__func__, "[%d, %d] No hits, stop here. \n", eg->id,
					eg->len);
			break;
		}
		forward(cur_pool, c, eg, ori);
		ext_con(eg->contig, c, 0);
		eg->len = eg->contig->len;
		ext_que(query, c, ori);
	}
	if (ori)
		seq_reverse(eg->len, eg->contig->seq, 0);
	free_alg(aligns);
	free(next);
	free_pool(cur_pool);
	bwa_free_read_seq(1, query);
	return eg;
}

/**
 * Extend a solid read.
 * First round:             ------------  =>
 * Second round:         <= --------------------------
 * Third round:     ---------------------------------- =>
 * Fourth round: <= --------------------------------------------
 *
 * That is because, after first round, there may not be enough mates to use in mate_pool because length is not enough.
 */
edge *pe_ext(const hash_table *ht, bwa_seq_t *query, const int tid) {
	int round_3_len = 0;
	edge *eg = NULL;
	pool *init_pool = NULL;
	bwa_seq_t *second_round_q = NULL;
	ubyte_t *rev = NULL;

	init_pool = get_start_pool(ht, query, 0, 0);
	if (!init_pool || init_pool->n == 0 || (stage <= 2 && bases_sup_branches(
			init_pool, 0, STRICT_BASES_SUP_THRE))) {
		show_msg(__func__, "Read %s may be in junction area, skip. \n",
				query->name);
		query->status = TRIED;
		query->tid = -1;
		free_pool(init_pool);
		return NULL;
	}
	free_pool(init_pool);
	//p_query(__func__, query);
	//show_debug_msg(__func__, "Extending to the right... \n");
	eg = pair_extension(NULL, ht, query, 0, tid);
	rev_reads_pos(eg);
	//show_debug_msg(__func__, "Edge after right extension: [%d, %d]\n", eg->id,
	//		eg->len);
	//show_debug_msg(__func__, "Got edge [%d, %d]. Extending to the left... \n",
	//		eg->id, eg->len);
	pair_extension(eg, ht, query, 1, tid);
	//show_debug_msg(__func__, "Edge after left extension: [%d, %d]\n", eg->id,
	//		eg->len);
	upd_reads(ht, eg, MISMATCHES);

	//show_debug_msg(__func__, "Second round: extending to the right... \n");
	second_round_q = new_seq(eg->contig, query->len, eg->len - query->len);
	pair_extension(eg, ht, second_round_q, 0, tid);
	round_3_len = eg->len;
	rev_reads_pos(eg);
	bwa_free_read_seq(1, second_round_q);

	//show_debug_msg(__func__, "Second round: extending to the left... \n");
	second_round_q = new_seq(eg->contig, query->len, 0);
	pair_extension(eg, ht, second_round_q, 1, tid);
	if (eg->len - round_3_len > 2) {
		upd_reads(ht, eg, MISMATCHES);
	} else {
		rev_reads_pos(eg);
	}

	rev = eg->contig->rseq;
	free(rev);
	rev = (ubyte_t*) malloc(eg->len + 1);
	memcpy(rev, eg->contig->seq, eg->len);
	seq_reverse(eg->len, rev, 1);
	eg->contig->rseq = rev;
	//log_edge(eg);
	//p_readarray(eg->reads, 1);
	return eg;
}

void est_insert_size(int n_max_pairs, char *lib_file, char *solid_file) {
	hash_table *ht = NULL;
	bwa_seq_t *query = NULL;
	FILE *solid = xopen(solid_file, "r");
	char line[80];
	int index = 0, n_pairs = 0, n_part_pairs = 0;
	int line_no = 0, start = 400;
	edge *eg = NULL;
	double *pairs = NULL, *partial_pairs = NULL, mean_ins_size = 0,
			sd_ins_size = 0;
	show_msg(__func__, "Library: %s \n", lib_file);
	show_msg(__func__, "Solid Reads: %s \n", solid_file);
	show_msg(__func__, "Maximum Pairs: %d \n", n_max_pairs);

	ht = pe_load_hash(lib_file);
	pairs = (double*) calloc(n_max_pairs + 1, sizeof(double));
	while (fgets(line, 80, solid) != NULL && n_pairs < n_max_pairs) {
		index = atoi(line);
		query = &ht->seqs[index];
		line_no++;
		if (query->status != FRESH || line_no < start)
			continue;
		show_msg(__func__,
				"---------- [%d] Processing read %d: %s ----------\n", line_no,
				pair_ctg_id, query->name);
		eg = pe_ext(ht, query, 0);
		if (!eg)
			continue;
		keep_pairs_only(eg, ht->seqs);
		g_ptr_array_sort(eg->pairs, (GCompareFunc) cmp_reads_by_name);
		partial_pairs = get_pair_dis_on_edge(eg, &n_part_pairs);
		if (n_part_pairs >= PAIRS_PER_EDGE)
			n_part_pairs = PAIRS_PER_EDGE;
		concat_doubles(pairs, &n_pairs, partial_pairs, n_part_pairs);
		free(partial_pairs);
		partial_pairs = NULL;
		p_ctg_seq("Contig now", eg->contig);
		show_msg(__func__, "Pairs of contig %d: [%d, %d] \n", eg->id, eg->len,
				eg->pairs->len);
		show_msg(__func__, "+%d, %d pairs now \n", n_part_pairs, n_pairs);
		n_part_pairs = 0;
		//p_readarray(eg->reads, 1);
		//log_edge(eg);
		destroy_eg(eg);
		eg = NULL;
	}
	mean_ins_size = mean(pairs, n_pairs);
	sd_ins_size = std_dev(pairs, n_pairs);
	show_msg(__func__, "Mean Insert Size: %f \n", mean_ins_size);
	show_msg(__func__, "Standard Deviation of Insert Size: %f \n", sd_ins_size);
}

/**
 * Keep only pairs on one edge
 */
int validate_edge(edgearray *all_edges, edge *eg, hash_table *ht,
		int *n_paired_reads, int *n_used_reads) {
	double pair_by_reads_perc = 0;
	if (stage == 1)
		pair_by_reads_perc = VALID_PAIR_PERC_STAGE_1;
	if (stage == 2)
		pair_by_reads_perc = VALID_PAIR_PERC_STAGE_2;
	if (eg) {
		//keep_pairs_only(eg, ht->seqs);
		if ((eg->len < 100 && eg->reads->len * ht->seqs->len < eg->len * 10)
				|| eg->pairs->len < eg->reads->len * pair_by_reads_perc
				|| eg->pairs->len > eg->reads->len || eg->reads->len == 0) { // || eg->pairs->len <= MIN_VALID_PAIRS) {
			show_msg(
					__func__,
					"ABANDONED [%d] %s: length %d, reads %d=>%d. Used reads %d/%d; Pair reads: %d/%d \n",
					eg->id, eg->name, eg->len, eg->reads->len, eg->pairs->len,
					*n_used_reads, ht->n_seqs, *n_paired_reads, ht->n_seqs);
			show_debug_msg(
					__func__,
					"ABANDONED [%d] %s: length %d, reads %d=>%d. Used reads %d/%d; Pair reads: %d/%d \n",
					eg->id, eg->name, eg->len, eg->reads->len, eg->pairs->len,
					*n_used_reads, ht->n_seqs, *n_paired_reads, ht->n_seqs);
			// upd_ctg_id(eg, -1, TRIED);
			eg->alive = 0;
			g_mutex_lock(sum_mutex);
			*n_used_reads += eg->reads->len;
			destroy_eg(eg);
			//g_ptr_array_add(all_edges, eg);
			g_mutex_unlock(sum_mutex);
			return 0;
		} else {
			g_mutex_lock(sum_mutex);
			*n_paired_reads += eg->pairs->len;
			*n_used_reads += eg->reads->len;
			g_ptr_array_add(all_edges, eg);
			g_mutex_unlock(sum_mutex);
			reset_tid(eg);
			show_msg(
					__func__,
					"[%d] %s: length %d, reads %d=>%d. Used reads %d/%d; Pair reads: %d/%d \n",
					eg->id, eg->name, eg->len, eg->reads->len, eg->pairs->len,
					*n_used_reads, ht->n_seqs, *n_paired_reads, ht->n_seqs);
			show_debug_msg(
					__func__,
					"[%d] %s: length %d, reads %d=>%d. Used reads %d/%d; Pair reads: %d/%d \n",
					eg->id, eg->name, eg->len, eg->reads->len, eg->pairs->len,
					*n_used_reads, ht->n_seqs, *n_paired_reads, ht->n_seqs);
			//log_edge(eg, ht->seqs);
			return 1;
		}
	}
	return 0;
}

typedef struct {
	edgearray *all_edges;
	hash_table *ht;
	int *n_paired_reads;
	int *n_used_reads;
	double stop_thre;
} thread_aux_t;

/**
 * Initiate a thread to assemble from some solid reads
 */
void *pe_lib_thread(gpointer solid_read, gpointer data) {
	int pre_n_single = 0, *n_paired_reads = NULL, *n_used_reads = NULL;
	bwa_seq_t *query = NULL, *seqs = NULL;
	edge *eg = NULL;
	thread_aux_t *d = (thread_aux_t*) data;
	hash_table *ht = NULL;
	ht = d->ht;
	seqs = ht->seqs;
	n_paired_reads = d->n_paired_reads;
	n_used_reads = d->n_used_reads;
	pre_n_single = *n_used_reads;
	uint32_t tid = 0;

	if (*n_used_reads > d->ht->n_seqs * d->stop_thre)
		return NULL;
	query = (bwa_seq_t*) solid_read;
	tid = atoi(query->name);
	if (query->status != FRESH) {
		return NULL;
	}

	// If this read is currently used by another thread
	if (query->is_in_c_pool > 0 && query->is_in_c_pool != tid)
		return NULL;
	if (query->is_in_m_pool > 0 && query->is_in_m_pool != tid)
		return NULL;
	if (has_n(query, 4) || is_biased_q(query) || has_rep_pattern(query)
			|| is_repetitive_q(query)) {
		query->status = TRIED;
		return NULL;
	}
	query->status = TRIED;
	query->tid = tid;
	eg = pe_ext(d->ht, query, tid);
	query->tid = -1;
	if (!eg) {
		g_mutex_lock(sum_mutex);
		*n_used_reads += 1;
		g_mutex_unlock(sum_mutex);
	} else {
		validate_edge(d->all_edges, eg, d->ht, d->n_paired_reads,
				d->n_used_reads);
	}
	//rescue_reads(seqs, ht->n_seqs);
	return NULL;
}

void run_threads(edgearray *all_edges, readarray *solid_reads, hash_table *ht,
		int *n_paired_reads, int *n_used_reads, double stop_thre) {
	thread_aux_t *data = NULL;
	bwa_seq_t *query = NULL;
	GThreadPool *thread_pool = NULL;
	int i = 0, j = 0, block_start = 0, nt = n_threads;
	data = (thread_aux_t*) calloc(1, sizeof(thread_aux_t));
	data->all_edges = all_edges;
	data->ht = ht;
	data->n_paired_reads = n_paired_reads;
	data->n_used_reads = n_used_reads;
	data->stop_thre = stop_thre;

	thread_pool
			= g_thread_pool_new((GFunc) pe_lib_thread, data, nt, TRUE, NULL);
	if (thread_pool == NULL) {
		err_fatal(__func__, "Failed to start the thread pool. \n");
	}
	while (block_start + JUMP_UNIT * nt < solid_reads->len) {
		for (i = 0; i < JUMP_UNIT; i++) {
			for (j = 0; j < nt; j++) {
				query = g_ptr_array_index(solid_reads,
						block_start + i + JUMP_UNIT * j);
				//if (should_start(query))
				g_thread_pool_push(thread_pool, (gpointer) query, NULL);
			}
		}
		block_start += JUMP_UNIT * nt;
	}
	for (i = block_start; i < solid_reads->len; i++) {
		query = g_ptr_array_index(solid_reads, i);
		//if (should_start(query))
		g_thread_pool_push(thread_pool, (gpointer) query, NULL);
	}
	g_thread_pool_free(thread_pool, 0, 1);
	free(data);
}

void post_process_edges(hash_table *ht, edgearray *all_edges, char *lib_file) {
	edgearray *final_paths = NULL;
	FILE *pair_contigs = NULL;
	FILE *merged_pair_contigs = NULL;
	char *name = NULL, *reads_name = NULL, *psl_name = NULL, blat_cmd[BUFSIZE];
	bwa_seq_t *seqs = NULL;
	GPtrArray *hits = NULL;
	hash_table *fresh_ht = NULL;

	seqs = ht->seqs;
	clean_edges(ht, all_edges);
	show_msg(__func__, "Total valid edges reported: %d \n", all_edges->len);
	show_debug_msg(__func__, "Total valid edges reported: %d \n",
			all_edges->len);

	name = get_output_file("paired.fa");
	pair_contigs = xopen(name, "w");
	//reset_read_ctg_id(ht->seqs, ht->n_seqs);
	reset_edge_ids(all_edges);
	save_edges(all_edges, pair_contigs, 0, 0, 0);
	fflush(pair_contigs);
	free(name);

	if (debug_mode) {
		name = get_output_file("roadmap.graph");
		reads_name = get_output_file("roadmap.reads");
		dump_rm(all_edges, name, reads_name);
		free(reads_name);
		free(name);
	}

	show_msg(__func__, "========================================== \n\n");
	show_msg(__func__, "Merging edges by overlapping... \n");
	// The merging assumes that the 'all_edges' are with contig ids 0,1,2,3...
	merge_ol_edges(all_edges, insert_size, sd_insert_size, ht, n_threads);
	//reset_read_ctg_id(ht->seqs, ht->n_seqs);
	reset_edge_ids(all_edges);
	name = get_output_file("merged.fa");
	merged_pair_contigs = xopen(name, "w");
	save_edges(all_edges, merged_pair_contigs, 0, 0, 0);
	fflush(merged_pair_contigs);
	fclose(merged_pair_contigs);
	free(name);

	if (debug_mode) {
		name = get_output_file("roadmap.1.graph");
		reads_name = get_output_file("roadmap.1.reads");
		dump_rm(all_edges, name, reads_name);
		free(reads_name);
		free(name);
	}

	clock_gettime(CLOCK_MONOTONIC, &finish_time);
	show_msg(__func__, "Template merging finished: %.2f sec\n",
			(float) (finish_time.tv_sec - start_time.tv_sec));

	show_msg(__func__, "========================================== \n\n ");
	show_msg(__func__, "Blatting merged contigs to find overlapping...\n");
	psl_name = get_output_file("merged.merged.psl");
	name = get_output_file("merged_pair_contigs.fa");
	show_msg(__func__, "%s %s %s %s ... \n", blat_exe, name, name, psl_name);
	sprintf(blat_cmd, "%s %s %s %s", blat_exe, name, name, psl_name);
	if (system(blat_cmd) != 0) {
		err_fatal(__func__, "Failed to call 'system' to execute %s \n",
				blat_cmd);
	}
	clock_gettime(CLOCK_MONOTONIC, &finish_time);
	show_msg(__func__, "Blat finished: %.2f sec\n", (float) (finish_time.tv_sec
			- start_time.tv_sec));
	free(name);

	show_msg(__func__, "========================================== \n\n ");
	show_msg(__func__, "Removing sub edges... \n");
	hits = read_blat_hits(psl_name);
	g_ptr_array_sort(hits, (GCompareFunc) cmp_hit_by_qname);
	mark_sub_edge(all_edges, hits);
	free_blat_hits(hits);
	// The kmer list and pos list are shrinked, not used anymore
	//destroy_ht(ht);

	show_msg(__func__, "========================================== \n\n ");
	show_msg(__func__, "Scaffolding %d edges... \n", all_edges->len);
	//fresh_ht = pe_load_hash(lib_file);
	// The merging assumes that the 'all_edges' are with contig ids 0,1,2,3...
	//realign_by_blat(all_edges, fresh_ht, n_threads);
	all_edges = scaffolding(all_edges, insert_size, sd_insert_size, fresh_ht,
			n_threads, psl_name);
	clock_gettime(CLOCK_MONOTONIC, &finish_time);
	show_msg(__func__, "Scaffolding finished: %.2f sec\n",
			(float) (finish_time.tv_sec - start_time.tv_sec));
	free(psl_name);

	if (debug_mode) {
		//realign_by_blat(all_edges, ht, n_threads);
		name = get_output_file("roadmap.2.graph");
		reads_name = get_output_file("roadmap.2.reads");
		dump_rm(all_edges, name, reads_name);
		free(reads_name);
		free(name);
	}

	show_msg(__func__, "========================================== \n\n ");
	if (debug_mode) {
		show_msg(__func__, "Drawing the roadmap... \n");
		name = get_output_file("roadmap.dot");
		graph_by_edges(all_edges, name);
		free(name);
	}
	show_msg(__func__, "Reporting combinatorial paths... \n");
	final_paths = report_paths(all_edges, seqs);
	name = get_output_file("peta.fa");
	save_paths(final_paths, name, 100);

	clock_gettime(CLOCK_MONOTONIC, &finish_time);
	show_msg(__func__, "Saving finished: %.2f sec\n",
			(float) (finish_time.tv_sec - start_time.tv_sec));

	free(name);
	fclose(pair_contigs);
	g_ptr_array_free(all_edges, TRUE);
	clock_gettime(CLOCK_MONOTONIC, &finish_time);
	show_msg(__func__, "Post processing finished: %.2f sec\n",
			(float) (finish_time.tv_sec - start_time.tv_sec));
}

void test_run(hash_table *ht, char *lib_file) {
	bwa_seq_t *seqs = ht->seqs, *query = NULL;
	edge *eg = NULL;
	edgearray *edges = g_ptr_array_sized_new(32);
	char *name = NULL;
	FILE *single_contig = NULL;

	pair_ctg_id = 1000;
	query = &seqs[5627228];
	p_query(__func__, query);
	eg = pe_ext(ht, query, 5627228);
	if (eg) {
		show_msg(__func__, "[%d] %s: length %d, reads %d=>%d. \n", eg->id,
				eg->name, eg->len, eg->reads->len, eg->pairs->len);
		g_ptr_array_add(edges, eg);
	} else {
		show_msg(__func__, "Hey, not extended. \n");
	}
	name = get_output_file("single.fa");
	single_contig = xopen(name, "w");
	save_edges(edges, single_contig, 0, 0, 0);
	exit(1);
}

void test_scaffolding(hash_table *ht) {
	GPtrArray *paths = NULL;
	FILE *merged_pair_contigs = NULL;
	GPtrArray *hits = read_blat_hits("../SRR027876_trim/merged.merged.psl");
	edgearray *all_edges = load_rm(ht, "../SRR027876_trim/roadmap.1.graph",
			"../SRR027876_trim/roadmap.1.reads",
			"../SRR027876_trim/merged_pair_contigs.fa");

	show_msg(__func__, "Removing sub edges... \n");
	g_ptr_array_sort(hits, (GCompareFunc) cmp_hit_by_qname);
	mark_sub_edge(all_edges, hits);
	free_blat_hits(hits);

	realign_by_blat(all_edges, ht, n_threads);
	merged_pair_contigs = xopen("../SRR027876_trim/validated.fa", "w");
	save_edges(all_edges, merged_pair_contigs, 0, 0, 0);

	show_msg(__func__, "Scaffolding...\n");
	all_edges = scaffolding(all_edges, insert_size, sd_insert_size, ht,
			n_threads, "../SRR027876_trim/merged.merged.psl");

	show_msg(__func__, "Reporting paths of %d edges... \n", all_edges->len);
	paths = report_paths(all_edges, ht->seqs);
	show_msg(__func__, "Saving %d paths... \n", paths->len);
	save_paths(paths, "../SRR027876_trim/peta.fa", 100);
	exit(1);
}

void test_merge(hash_table *ht) {
	FILE *merged_pair_contigs = NULL;
	//	GPtrArray *hits =
	//			read_blat_hits("../SRR097897_out/merged.merged.psl");
	edgearray *all_edges = load_rm(ht, "../SRR097897_out/roadmap.graph",
			"../SRR097897_out/roadmap.reads",
			"../SRR097897_out/pair_contigs.fa");
	reset_edge_ids(all_edges);
	//realign_by_blat(all_edges, ht, n_threads);
	merge_ol_edges(all_edges, insert_size, sd_insert_size, ht, n_threads);
	merged_pair_contigs = xopen("../SRR097897_out/merged_pair_contigs.fa", "w");
	save_edges(all_edges, merged_pair_contigs, 0, 0, 0);
	//	g_ptr_array_sort(hits, (GCompareFunc) cmp_hit_by_qname);
	//	mark_sub_edge(all_edges, hits);
	//	reset_edge_ids(all_edges);
	//	merged_pair_contigs = xopen("../SRR097897_out/validated.2.fa", "w");
	//	save_edges(all_edges, merged_pair_contigs, 0, 0, 0);
	//	realign_by_blat(all_edges, ht, n_threads);
	//	dump_rm(all_edges, "../SRR097897_out/roadmap.1.graph",
	//			"../SRR097897_out/roadmap.1.reads");
	exit(1);
}

readarray *load_solid_reads(const char *solid_fn, bwa_seq_t *seqs,
		const int n_seqs) {
	int i = 0;
	char line[80];
	readarray *solid_reads = NULL;
	bwa_seq_t *query = NULL;
	FILE *solid = xopen(solid_fn, "r");

	solid_reads = g_ptr_array_sized_new(n_seqs / 10);
	while (fgets(line, 80, solid) != NULL) {
		i = atoi(line);
		query = &seqs[i];
		g_ptr_array_add(solid_reads, query);
	}
	fclose(solid);
	return solid_reads;
}

void consume_solid_reads(hash_table *ht, const double stop_thre,
		edgearray *all_edges, readarray *solid_reads, int *n_used_reads,
		int *n_paired_reads) {
	double n_unit = N_SPLIT_UNIT, unit_perc = 0, max_perc = 0;
	int i = 0, n_used_reads_pre = 0;

	unit_perc = 1 / n_unit;
	for (i = 1; i <= n_unit; i++) {
		max_perc = i * unit_perc;
		if (max_perc >= stop_thre)
			max_perc = stop_thre;
		if (*n_used_reads > ht->n_seqs * max_perc)
			continue;
		n_used_reads_pre = *n_used_reads;
		run_threads(all_edges, solid_reads, ht, n_paired_reads, n_used_reads,
				max_perc);
		clock_gettime(CLOCK_MONOTONIC, &finish_time);
		show_msg(__func__, "Shrinking the hash table; %.2f sec... \n",
				(float) (finish_time.tv_sec - start_time.tv_sec));
		erase_reads_on_ht(ht);
		shrink_ht(ht);
		show_msg(__func__, "Stage %d, percentage %d/%.0f finished %.2f \n",
				stage, i, n_unit, i * unit_perc);
		show_msg(__func__,
				"Correcting counters: [single: %d], [paired: %d]... \n",
				*n_used_reads, *n_paired_reads);
		correct_used_reads_n(ht, n_used_reads, n_paired_reads);
		clock_gettime(CLOCK_MONOTONIC, &finish_time);
		show_msg(
				__func__,
				"Corrected counters: [single: %d], [paired: %d]; %.2f sec... \n",
				*n_used_reads, *n_paired_reads, (float) (finish_time.tv_sec
						- start_time.tv_sec));
		// Some extreme cases, the acutual used reads are not increased so much.
		i = *n_used_reads / (ht->n_seqs * max_perc);
		if (*n_used_reads >= stop_thre * ht->n_seqs || *n_used_reads
				== n_used_reads_pre)
			break;
	}
}

void pe_lib_core(int n_max_pairs, char *lib_file, char *solid_file) {
	hash_table *ht = NULL;
	bwa_seq_t *seqs = NULL, *s = NULL;
	int n_paired_reads = 0, n_used_reads = 0, n_rep = 0, i = 0;
	GPtrArray *all_edges = NULL;
	readarray *solid_reads = NULL;
	clean_opt *c_opt = NULL;

	show_msg(__func__, "Library: %s \n", lib_file);
	show_msg(__func__, "Solid Reads: %s \n", solid_file);

	all_edges = g_ptr_array_sized_new(BUFSIZ);
	ht = pe_load_hash(lib_file);
	test_run(ht, lib_file);
	//test_scaffolding(ht);
	//test_merge(ht);
	seqs = &ht->seqs[0];
	show_msg(__func__, "Removing repetitive reads... \n");
	n_rep = rm_repetitive_reads(seqs, ht->n_seqs);
	erase_reads_on_ht(ht);
	shrink_ht(ht);
	//n_used_reads = n_rep;
	//n_paired_reads = n_rep;
	show_msg(__func__, "# of repetitive reads are marked DEAD: %d\n", n_rep);
	//show_msg(__func__, "Correcting reads... \n");
	//correct_reads(ht, 6);
	//save_unpaired_seqs("../SRR027876_out/clean2.fa", ht->seqs, ht->n_seqs);

	solid_reads = load_solid_reads(solid_file, seqs, ht->n_seqs);
	clock_gettime(CLOCK_MONOTONIC, &finish_time);
	show_msg(__func__, "Solid reads loaded: %.2f sec\n",
			(float) (finish_time.tv_sec - start_time.tv_sec));

	stage = 1;
	show_msg(__func__, "========================================== \n");
	show_msg(__func__, "Stage 1/3: Trying to use up the paired reads... \n");
	consume_solid_reads(ht, STOP_THRE_STAGE_1, all_edges, solid_reads,
			&n_used_reads, &n_paired_reads);
	clock_gettime(CLOCK_MONOTONIC, &finish_time);
	show_msg(__func__, "Stage 1 finished: %.2f sec\n",
			(float) (finish_time.tv_sec - start_time.tv_sec));

	stage = 2;
	show_msg(__func__, "========================================== \n");
	show_msg(__func__,
			"Stage 2/3: Trying to assembly more unpaired reads... \n");
	rescue_reads(seqs, ht->n_seqs);
	correct_used_reads_n(ht, &n_used_reads, &n_paired_reads);
	c_opt = init_clean_opt();
	c_opt->kmer = 15;
	c_opt->stop_thre = 1;
	c_opt->n_threads = n_threads;
	g_ptr_array_free(solid_reads, TRUE);
	show_msg(__func__,
			"Stage 2/3: Calculating solid reads from unpaired reads... \n");
	solid_reads = calc_solid_reads(ht->seqs, ht->n_seqs, c_opt, (ht->n_seqs
			- n_paired_reads) * c_opt->stop_thre, 1, 1);
	consume_solid_reads(ht, STOP_THRE_STAGE_2, all_edges, solid_reads,
			&n_used_reads, &n_paired_reads);
	free(c_opt);
	clock_gettime(CLOCK_MONOTONIC, &finish_time);
	show_msg(__func__, "Stage 2 finished: %.2f sec\n",
			(float) (finish_time.tv_sec - start_time.tv_sec));

	stage = 3;
	show_msg(__func__, "========================================== \n");
	show_msg(__func__, "Stage 3/3: Trying to assembly unpaired reads... \n");
	rescue_reads(seqs, ht->n_seqs);
	correct_used_reads_n(ht, &n_used_reads, &n_paired_reads);

	while (solid_reads->len > 0) {
		g_ptr_array_remove_index_fast(solid_reads, 0);
	}
	for (i = 0; i < ht->n_seqs; i++) {
		s = &ht->seqs[i];
		if (s->status == FRESH)
			g_ptr_array_add(solid_reads, s);
	}
	consume_solid_reads(ht, STOP_THRE_STAGE_3, all_edges, solid_reads,
			&n_used_reads, &n_paired_reads);
	g_ptr_array_free(solid_reads, TRUE);
	clock_gettime(CLOCK_MONOTONIC, &finish_time);
	show_msg(__func__, "Stage 3 finished: %.2f sec\n",
			(float) (finish_time.tv_sec - start_time.tv_sec));

	// ht is destoryed in this function
	post_process_edges(ht, all_edges, lib_file);
}

int pe_lib_usage() {
	show_msg(__func__,
			"Command: ./peta pair -p MAX_PAIRS read_file starting_reads \n");
	return 1;
}

int pe_lib(int argc, char *argv[]) {
	int c = 0, n_max_pairs = 0;
	clock_gettime(CLOCK_MONOTONIC, &start_time);
	while ((c = getopt(argc, argv, "p:k:m:s:o:t:b:d:")) >= 0) {
		switch (c) {
		case 'p':
			n_max_pairs = atoi(optarg);
			break;
		case 'k':
			overlap_len = atoi(optarg);
			break;
		case 'o':
			out_root = optarg;
			break;
		case 'm':
			insert_size = atoi(optarg);
			break;
		case 's':
			sd_insert_size = atoi(optarg);
			break;
		case 't':
			n_threads = atoi(optarg);
			break;
		case 'b':
			blat_exe = optarg;
			break;
		case 'd':
			debug_mode = atoi(optarg);
			break;
		}
	}
	//test_smith_waterman(argv[optind]);
	if (optind + 2 > argc) {
		return pe_lib_usage();
	}
	if (!g_thread_supported())
		g_thread_init(NULL);

	id_mutex = g_mutex_new();
	sum_mutex = g_mutex_new();
	show_msg(__func__, "Maximum pairs: %d \n", n_max_pairs);
	show_msg(__func__, "Overlap length: %d \n", overlap_len);
	show_msg(__func__, "Output folder: %s \n", out_root);
	show_msg(__func__, "Insert size: %d \n", insert_size);
	show_msg(__func__, "Standard deviation: %d \n", sd_insert_size);
	if (n_max_pairs > 0) {
		est_insert_size(n_max_pairs, argv[optind], argv[optind + 1]);
	} else {
		pe_lib_core(n_max_pairs, argv[optind], argv[optind + 1]);
	}
	clock_gettime(CLOCK_MONOTONIC, &finish_time);
	show_msg(__func__, "Done: %.2f sec\n", (float) (finish_time.tv_sec
			- start_time.tv_sec));
	return 0;
}
