#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <glib.h>
#include <math.h>
#include "rand.h"
#include "ass.h"
#include "clean.h"
#include "bwase.h"
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
#include <pthread.h>
//#include <glib_global.h> // Always include as the last include.

int pair_ctg_id = 0;
int overlap_len = 0;
int insert_size = 0;
int sd_insert_size = 0;
char *out_root = NULL;
int n_threads = 1;
int stage = 1;
GMutex *id_mutex;
GMutex *sum_mutex;

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

gint cmp_reads_by_name(gpointer a, gpointer b) {
	bwa_seq_t *seq_a = *((bwa_seq_t**) a);
	bwa_seq_t *seq_b = *((bwa_seq_t**) b);
	return ((atoi(seq_a->name)) - atoi(seq_b->name));
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
					|| eg->id != mate->contig_id) {
				to_remove = 1;
			} else {// If the distance is not within range, remove 'read'
				if (abs(eg->len - mate->shift) > (insert_size + sd_insert_size
						* SD_TIMES)) {
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
			if (pool_rm_index(cur_pool, i))
				i--;
		}
	}
}

/**
 * From current edge, get all mates of the used reads.
 */
pool *get_mate_pool_from_edge(edge *eg, const hash_table *ht, const int ori) {
	int i = 0;
	bwa_seq_t *s = NULL, *mate = NULL, *seqs = NULL;
	pool *mate_pool = NULL;
	mate_pool = new_pool();
	seqs = ht->seqs;
	for (i = 0; i < eg->reads->len; i++) {
		s = g_ptr_array_index(eg->reads, i);
		mate = get_mate(s, seqs);
		//p_query("READ", s);
		//p_query("MATE", mate);
		// If the mate should have been used.
		if (is_paired(s, ori))
			continue;
		// If the insert size is not in the range
		if (abs(eg->len - s->shift) > (insert_size + sd_insert_size * SD_TIMES)) {
			continue;
		}
		// If the mate is already in use, either by current or another thread
		if (mate->is_in_c_pool || mate->is_in_m_pool || mate->status == USED)
			continue;
		// If the used read is used by another thread;
		//	or the mate has been used by this template before.
		if (!(s->status == TRIED && s->contig_id == eg->id) || (mate->status
				== TRIED && mate->contig_id == eg->id))
			continue;
		mate->rev_com = s->rev_com;
		mate_pool_add(mate_pool, mate, eg->tid);
	}
	return mate_pool;
}

/**
 * Add read to current pool from the mate pool.
 * Check whether a mate overlaps with the tail with length parameter 'nm'
 */
void add_mates_by_ol(const hash_table *ht, edge *eg, pool *cur_pool,
		const int ol, const int nm, bwa_seq_t *query, const int ori) {
	int i = 0, overlapped = 0;
	bwa_seq_t *mate = NULL, *tmp = NULL, *s = NULL;
	readarray *ol_mates = NULL;
	bwa_seq_t *template = NULL, *seqs = NULL;
	pool *mate_pool = NULL;
	reads_ht *rht = NULL;
	seqs = ht->seqs;
	// Copy read length of the end of the contig.
	template = new_seq(eg->contig, seqs->len, eg->len - seqs->len);
	if (ori) {
		seq_reverse(template->len, template->seq, 0);
	}
	mate_pool = get_mate_pool_from_edge(eg, ht, ori);
	//p_readarray(mate_pool->reads, 1);
	if (mate_pool->n >= N_BIG_MATE_POOL) {
		rht = build_reads_ht(ol, mate_pool->reads);
		ol_mates = find_reads_ol_template(rht, template, seqs, ori);
	} else {
		ol_mates = mate_pool->reads;
	}
	//p_readarray(ol_mates, 1);
	// Add the mate reads which overlap with the tail into the current pool
	for (i = 0; i < ol_mates->len; i++) {
		mate = g_ptr_array_index(ol_mates, i);
		s = get_mate(mate, seqs);
		// For the reads in the mate pool, these ones are not considered:
		//	1. Is already in c_pool (in this or another thread), or in the mate pool of other thread
		//	2. Is already used
		//	3. Its mate is not used
		//	4. Its mate is used, by by another edge
		//	5. The distance between the mates are out of range
		tmp = mate;
		if (mate->rev_com)
			tmp = new_mem_rev_seq(mate, mate->len, 0);
		overlapped = find_ol_within_k(tmp, template, nm, ol - 1, query->len - 1,
				ori);
		/*if (strcmp(mate->name, "2460877") == 0) {
			show_debug_msg("ORI", "ORI: %d \n", ori);
			p_ctg_seq("QUERY", query);
			p_query("MATE", tmp);
			p_query("ORIG", mate);
			p_query("USED", get_mate(mate, seqs));
			show_debug_msg(__func__, "OVERLAP 1: %d \n", overlapped);
		}*/
		if (overlapped >= ol) {
			// Only if this mate overlaps with some read in the cur_pool, add it.
			// It is important because sometimes it maybe added just for coincidence.
			mate->cursor = ori ? (mate->len - overlapped - 1) : overlapped;
			pool_add(cur_pool, mate, eg->tid);
			//show_debug_msg(__func__, "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");
			//show_debug_msg(__func__, "ORI%d \n", ori);
			//p_query("ADDED MATE", mate);
			//p_ctg_seq("TEMPLATE", template);
			//show_debug_msg(__func__, "Overlapped: %d \n", overlapped);
			//show_debug_msg(__func__, "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");
			rm_read_from_ht(rht, mate);
		}
		if (tmp != mate)
			bwa_free_read_seq(1, tmp);
	}
	if (mate_pool->n >= N_BIG_MATE_POOL)
		g_ptr_array_free(ol_mates, TRUE);
	bwa_free_read_seq(1, template);
	destroy_reads_ht(rht);
	free_mate_pool(mate_pool);
}

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
	int i = 0, index = 0, pre_cursor = 0;
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
		pre_cursor = s->cursor;

		if (s->is_in_c_pool || (s->is_in_m_pool && s->is_in_m_pool != eg->tid)
				|| s->status == USED || (s->status == TRIED && s->contig_id
				== eg->id))
			continue;

		// Stage 3 only use single reads, so TRIED reads are not used any more
		if (stage == 3) {
			if (s->status != FRESH)
				continue;
		} else {
			if (mate->status == TRIED && mate->contig_id == eg->id) {
				if (is_paired(mate, ori))
					continue;
			}
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
		if ((!ori && pre_cursor > s->cursor) || (ori && pre_cursor < s->cursor))
			continue;
		pool_add(cur_pool, s, eg->tid);
	}
	//show_debug_msg(__func__, "Removing partial... \n");
	// In current pool, if a read does not overlap with the tail properly, remove it
	rm_partial(eg, cur_pool, ori, seqs, query, 2);
	//p_pool("AFTER REMOVING PARTIAL", cur_pool, NULL);
	// Keep only the reads whose mate is used previously.
	if (stage != 3 && eg->len >= (insert_size + sd_insert_size * SD_TIMES)) {
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
		if (stage != 3 && mate->status == USED
				&& (mate->rev_com != a->rev_comp)) {
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
		if (s->contig_id == INVALID_CONTIG_ID || s->status == USED || s->cursor
				< -1 || s->cursor > s->len || s->is_in_c_pool
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
 * Correct the starting contig of an edge
 */
void correct_start(edge *eg, pool *cur_pool) {
	bwa_seq_t *read = NULL;
	int i = 0, j = 0, cursor = 0, has_more = 1;
	int *next = NULL;
	int correct_c = 0;
	if (!eg || !cur_pool)
		return;
	next = (int*) calloc(5, sizeof(int));
	for (j = 0; j < eg->contig->len; j++) {
		if (!has_more)
			break;
		has_more = 0;
		reset_c(next, NULL);
		for (i = 0; i < cur_pool->reads->len; i++) {
			read = g_ptr_array_index(cur_pool->reads, i);
			cursor = (read->cursor - read->len) + j;
			if (cursor >= 0 && cursor < read->len) {
				if (read->rev_com)
					check_c(next, read->rseq[cursor]);
				else
					check_c(next, read->seq[cursor]);
				has_more = 1;
			}
		}
		correct_c = get_pure_most(next);
		eg->contig->seq[j] = correct_c;
	}
	//p_ctg_seq("CORRECTED", eg->contig);
	free(next);
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
		//p_query(__func__, query);
		reset_c(next, NULL); // Reset the counter
		if (!cur_pool || is_repetitive_q(query)) {
			show_msg(__func__, "[%d, %d] Repetitive pattern, stop!\n", eg->id,
					eg->len);
			p_query("Repetitive pattern", query);
			break;
		}
		//p_query(__func__, query);
		pe_aln_query(query, query->seq, ht, MISMATCHES, query->len, 0, aligns);
		pe_aln_query(query, query->rseq, ht, MISMATCHES, query->len, 1, aligns);
		//p_align(aligns);
		maintain_pool(aligns, ht, cur_pool, eg, query, next, ori);
		// If no read in current pool, try less stringent overlapped length and zero mismatches
		if (stage != 3 && eg->len > insert_size - sd_insert_size
				&& cur_pool->reads->len == 0) {
			add_mates_by_ol(ht, eg, cur_pool, RELAX_MATE_OL_THRE,
					SHORT_MISMATCH, query, ori);
			reset_c(next, NULL); // Reset the counter
			check_next_char(cur_pool, eg, next, ori);
			//p_pool("After adding mates", cur_pool, next);
		}
		reset_alg(aligns);
		//show_debug_msg(__func__, "Edge %d, length %d \n", eg->id, eg->len);
		//p_ctg_seq("Contig", eg->contig);
		//p_pool("Current Pool", cur_pool, next);
		c = get_pure_most(next);
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
	clock_t t = clock();

	init_pool = get_start_pool(ht, query, 0, 0);
	if (!init_pool || init_pool->n == 0 || bases_sup_branches(init_pool, 0,
			STRICT_BASES_SUP_THRE)) {
		show_msg(__func__, "Read %s may be in junction area, skip. \n",
				query->name);
		query->status = TRIED;
		free_pool(init_pool);
		return NULL;
	}
	free_pool(init_pool);
	p_query(__func__, query);
	show_debug_msg(__func__, "Extending to the right... \n");
	eg = pair_extension(NULL, ht, query, 0, tid);
	rev_reads_pos(eg);
	show_debug_msg(__func__, "Edge after right extension: [%d, %d]\n", eg->id,
			eg->len);
	show_debug_msg(__func__, "Got edge [%d, %d]. Extending to the left... \n",
			eg->id, eg->len);
	pair_extension(eg, ht, query, 1, tid);
	show_debug_msg(__func__, "Edge after left extension: [%d, %d]\n", eg->id,
			eg->len);
	show_msg(__func__, "Before update: %d=>%d \n", eg->reads->len, eg->pairs->len);
	upd_reads(ht->seqs, eg, MISMATCHES);

	show_debug_msg(__func__, "Second round: extending to the right... \n");
	second_round_q = new_seq(eg->contig, query->len, eg->len - query->len);
	pair_extension(eg, ht, second_round_q, 0, tid);
	round_3_len = eg->len;
	rev_reads_pos(eg);
	bwa_free_read_seq(1, second_round_q);

	show_debug_msg(__func__, "Second round: extending to the left... \n");
	second_round_q = new_seq(eg->contig, query->len, 0);
	pair_extension(eg, ht, second_round_q, 1, tid);
	if (eg->len - round_3_len > 2) {
		upd_reads(ht->seqs, eg, MISMATCHES);
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
	show_msg(__func__, "Time for [%d, %d]: %.2f sec\n", eg->id, eg->len,
			(float) (clock() - t) / CLOCKS_PER_SEC);
	return eg;
}

/**
 * For the reads used on one edge, keep the pairs only.
 * For other reads, mark them as TRIED.
 * It means, these reads can be used again, but not used as a starting read.
 */
void keep_pairs_only(edge *eg, bwa_seq_t *seqs) {
	int i = 0, min_shift = 0, max_shift = 0, read_len = 0;
	bwa_seq_t *read = NULL, *mate = NULL;
	show_debug_msg(__func__, "Getting mate pairs...\n");
	//p_readarray(eg->reads, 1);
	for (i = 0; i < eg->reads->len; i++) {
		read = g_ptr_array_index(eg->reads, i);
		read_len = read->len;
		mate = get_mate(read, seqs);
		if (mate->status != FRESH && mate->contig_id == eg->id) {
			min_shift = read->shift < min_shift ? read->shift : min_shift;
			max_shift = read->shift > max_shift ? read->shift : max_shift;
			g_ptr_array_add(eg->pairs, read);
			g_ptr_array_add(eg->pairs, mate);
			mate->status = FRESH;
			read->status = FRESH;
		}
	}
	for (i = 0; i < eg->reads->len; i++) {
		read = g_ptr_array_index(eg->reads, i);
		read->status = TRIED;
	}
	for (i = 0; i < eg->pairs->len; i++) {
		read = g_ptr_array_index(eg->pairs, i);
		read->status = USED;
	}
}

void est_insert_size(int n_max_pairs, char *lib_file, char *solid_file) {
	hash_table *ht = NULL;
	bwa_seq_t *query = NULL, *seqs = NULL;
	FILE *solid = xopen(solid_file, "r");
	char line[80];
	int index = 0, ol = 0, n_pairs = 0, n_part_pairs = 0;
	int line_no = 0;
	edge *eg = NULL;
	double *pairs = NULL, *partial_pairs = NULL, mean_ins_size = 0,
			sd_ins_size = 0;
	show_msg(__func__, "Library: %s \n", lib_file);
	show_msg(__func__, "Solid Reads: %s \n", solid_file);
	show_msg(__func__, "Maximum Pairs: %d \n", n_max_pairs);

	ht = pe_load_hash(lib_file);
	seqs = &ht->seqs[0];
	ol = seqs->len / 2; // Read length
	pairs = (double*) calloc(n_max_pairs + 1, sizeof(double));
	while (fgets(line, 80, solid) != NULL && n_pairs < n_max_pairs) {
		index = atoi(line);
		query = &ht->seqs[index];

		if (query->status != FRESH)
			continue;
		show_msg(__func__,
				"---------- [%d] Processing read %d: %s ----------\n", line_no,
				pair_ctg_id, query->name);
		eg = pe_ext(ht, query, 0);
		if (!eg)
			continue;
		keep_pairs_only(eg, ht->seqs);
		if (eg->pairs->len * 2 < eg->reads->len) {
			destroy_eg(eg);
			continue;
		}
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
int validate_edge(edgearray *all_edges, edge *eg, const hash_table *ht,
		int *n_paired_reads, int *n_single_reads) {
	int n = 0;
	if (eg) {
		n = *n_paired_reads;
		if (stage == 3) {
			//p_ctg_seq(__func__, eg->contig);
			//p_readarray(eg->reads, 1);
			n = *n_single_reads;
		}
		// keep_pairs_only(eg, ht->seqs);
		if (eg->len < insert_size && eg->reads->len * ht->seqs->len < eg->len
				* 10) { // || eg->pairs->len <= MIN_VALID_PAIRS) {
			show_msg(
					__func__,
					"ABANDONED [%d] %s: length %d, reads %d=>%d. Total reads %d/%d \n",
					eg->id, eg->name, eg->len, eg->reads->len, eg->pairs->len,
					n, ht->n_seqs);
			show_debug_msg(
					__func__,
					"ABANDONED [%d] %s: length %d, reads %d=>%d. Total reads %d/%d \n",
					eg->id, eg->name, eg->len, eg->reads->len, eg->pairs->len,
					n, ht->n_seqs);
			upd_ctg_id(eg, -1, TRIED);
			eg->alive = 0;
			g_mutex_lock(sum_mutex);
			g_ptr_array_add(all_edges, eg);
			g_mutex_unlock(sum_mutex);
			return 0;
		} else {
			g_mutex_lock(sum_mutex);
			*n_paired_reads += eg->pairs->len;
			*n_single_reads += eg->reads->len;
			g_ptr_array_add(all_edges, eg);
			g_mutex_unlock(sum_mutex);
			show_msg(__func__,
					"[%d] %s: length %d, reads %d=>%d. Total reads %d/%d \n",
					eg->id, eg->name, eg->len, eg->reads->len, eg->pairs->len,
					n, ht->n_seqs);
			show_debug_msg(__func__,
					"[%d] %s: length %d, reads %d=>%d. Total reads %d/%d \n",
					eg->id, eg->name, eg->len, eg->reads->len, eg->pairs->len,
					n, ht->n_seqs);
			return 1;
		}
	}
	return 0;
}

void clean_edges(const hash_table *ht, edgearray *all_edges) {
	int i = 0;
	edge *eg = NULL;
	for (i = 0; i < all_edges->len; i++) {
		eg = g_ptr_array_index(all_edges, i);
		if (eg->alive == 0) {
			clear_used_reads(eg, 0);
			destroy_eg(eg);
			if (g_ptr_array_remove_index_fast(all_edges, i))
				i--;
		}
	}
}

typedef struct {
	edgearray *all_edges;
	readarray *solid_reads;
	hash_table *ht;
	int *n_paired_reads;
	int *n_single_reads;
	int start;
	int end;
	int tid;
	int stage;
	double stop_thre;
} thread_aux_t;

/**
 * Initiate a thread to assemble from some solid reads
 */
static void *pe_lib_thread(void *data) {
	int i = 0, *n_paired_reads = NULL, *n_single_reads = NULL;
	bwa_seq_t *query = NULL;
	edge *eg = NULL;
	thread_aux_t *d = (thread_aux_t*) data;
	hash_table *ht = NULL;
	show_debug_msg(__func__, "From %d to %d \n", d->start, d->end);
	ht = d->ht;
	for (i = d->start; i < d->end; i++) {
		query = g_ptr_array_index(d->solid_reads, i);
		n_paired_reads = d->n_paired_reads;
		n_single_reads = d->n_single_reads;
		//if (pair_ctg_id == 0)
		//	query = &ht->seqs[2592602];
		//if (pair_ctg_id == 1)
		//	query = &ht->seqs[3535214];
		//if (pair_ctg_id == 2)
		//	query = &ht->seqs[2738138];
		//		if (query->status != FRESH)
		//			continue;
		if (has_n(query) || is_biased_q(query) || has_rep_pattern(query)
				|| is_repetitive_q(query) || query->status == USED
				|| query->status == TRIED) {
			query->status = TRIED;
			continue;
		}
		query->status = TRIED;
		if (stage == 3) {
			if (*n_single_reads > d->ht->n_seqs * d->stop_thre)
				break;
		} else {
			if (*n_paired_reads > d->ht->n_seqs * d->stop_thre)
				break;
		}
		eg = pe_ext(d->ht, query, d->tid);

		validate_edge(d->all_edges, eg, d->ht, d->n_paired_reads,
				d->n_single_reads);
		eg = NULL;
		//if (pair_ctg_id >= 1)
		//	break;
	}
	return NULL;
}

void run_threads(edgearray *all_edges, readarray *solid_reads, hash_table *ht,
		int *n_paired_reads, int *n_single_reads, int start, int end,
		int n_per_threads, double stop_thre) {
	thread_aux_t *data;
	int i = 0;
	GThread *threads[n_threads];
	data = (thread_aux_t*) calloc(n_threads, sizeof(thread_aux_t));
	for (i = 0; i < n_threads; ++i) {
		data[i].all_edges = all_edges;
		data[i].solid_reads = solid_reads;
		data[i].ht = ht;
		data[i].n_paired_reads = n_paired_reads;
		data[i].n_single_reads = n_single_reads;
		data[i].start = start + n_per_threads * i;
		data[i].end = start + n_per_threads * (i + 1);
		if (data[i].end > end)
			data[i].end = end;
		data[i].tid = i + 1;
		data[i].stop_thre = stop_thre;
		threads[i]
				= g_thread_create((GThreadFunc) pe_lib_thread, data + i, TRUE, NULL);
	}

	/* wait for threads to finish */
	for (i = 0; i < n_threads; ++i) {
		g_thread_join(threads[i]);
	}
	free(data);
}

void reset_edge_ids(edgearray *all_edges) {
	int i = 0;
	edge *eg = NULL;
	char *ctg_name = malloc(16);
	for (i = 0; i < all_edges->len; i++) {
		eg = g_ptr_array_index(all_edges, i);
		if (eg->alive) {
			eg->id = i;
			// Set all 'contig_id' of reads to 'i', but not reset the status
			upd_ctg_id(eg, i, -1);
			set_rev_com(eg->contig);
			sprintf(ctg_name, "%d", eg->id);
			memcpy(eg->contig->name, ctg_name, 16);
		} else {
			if (g_ptr_array_remove_index_fast(all_edges, i))
				i--;
		}
	}
	free(ctg_name);
}

void post_process_edges(const hash_table *ht, edgearray *all_edges) {
	edgearray *final_paths = NULL;
	FILE *pair_contigs = NULL;
	FILE *merged_pair_contigs = NULL;
	char *name = NULL, *reads_name = NULL;
	clock_t t = clock();
	bwa_seq_t *seqs = NULL;
	edge *eg = NULL;
	int i = 0;

	seqs = ht->seqs;
	clean_edges(ht, all_edges);
	show_msg(__func__, "Total valid edges reported: %d \n", all_edges->len);
	show_debug_msg(__func__, "Total valid edges reported: %d \n",
			all_edges->len);

	name = get_output_file("pair_contigs.fa");
	pair_contigs = xopen(name, "w");
	save_edges(all_edges, pair_contigs, 0, 0, 100);
	fflush(pair_contigs);
	free(name);

	reset_edge_ids(all_edges);

	show_msg(__func__, "========================================== \n\n");
	show_msg(__func__, "Merging edges by overlapping... \n");
	merge_ol_edges(all_edges, insert_size, ht, n_threads);
	name = get_output_file("merged_pair_contigs.fa");
	merged_pair_contigs = xopen(name, "w");
	save_edges(all_edges, merged_pair_contigs, 0, 0, 100);
	fflush(merged_pair_contigs);
	free(name);
	show_msg(__func__, "Template merging finished: %.2f sec\n",
			(float) (clock() - t) / CLOCKS_PER_SEC);

	for (i = 0; i < all_edges->len; i++) {
		eg = g_ptr_array_index(all_edges, i);
		set_rev_com(eg->contig);
		g_ptr_array_sort(eg->reads, (GCompareFunc) cmp_read_by_name);
	}

	show_msg(__func__, "========================================== \n\n ");
	show_msg(__func__, "Scaffolding %d edges... \n", all_edges->len);
	name = get_output_file("roadmap.graph");
	reads_name = get_output_file("roadmap.reads");
	dump_rm(all_edges, name, reads_name);
	free(name);
	free(reads_name);
	scaffolding(all_edges, insert_size, ht, n_threads);
	show_msg(__func__, "Scaffolding finished: %.2f sec\n",
			(float) (clock() - t) / CLOCKS_PER_SEC);

	show_msg(__func__, "========================================== \n\n ");
	show_msg(__func__, "Saving the roadmap... \n");
	name = get_output_file("roadmap.dot");
	graph_by_edges(all_edges, name);
	free(name);
	show_msg(__func__, "Reporting combinatorial paths... \n");
	final_paths = report_paths(all_edges, seqs);
	name = get_output_file("peta.fa");
	save_paths(final_paths, name, 100);

	show_msg(__func__, "Saving finished: %.2f sec\n", (float) (clock() - t)
			/ CLOCKS_PER_SEC);

	free(name);
	fclose(pair_contigs);
	fclose(merged_pair_contigs);
	g_ptr_array_free(all_edges, TRUE);
	show_msg(__func__, "Post processing finished: %.2f sec\n", (float) (clock()
			- t) / CLOCKS_PER_SEC);
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

void pe_lib_single(const hash_table *ht, edgearray *all_edges,
		int *n_paired_reads, int *n_single_reads) {
	char *name = NULL, *lib_name = NULL, *solid_name = NULL;
	int n_unpaired = 0, n_per_threads = 0, n_single = 0, i = 0;
	hash_table *pht = NULL;
	hash_opt *h_opt = NULL;
	clean_opt *c_opt = NULL;
	readarray *solid_reads = NULL;
	edgearray *single_edges = NULL;
	edge *eg = NULL;
	clock_t t = clock();

	stage = 3;
	name = get_output_file("part.fa");
	n_unpaired = save_unpaired_seqs(name, ht->seqs, ht->n_seqs);
	show_msg(__func__, "Reads unpaired: %d, saved to %s \n", n_unpaired, name);

	// Hash the unpaired reads
	h_opt = init_hash_opt();
	h_opt->read_len = ht->seqs->len;
	h_opt->block_size = 9;
	h_opt->n_hash_block = (h_opt->read_len - h_opt->interleaving * h_opt->k)
			/ h_opt->block_size;
	pe_hash_core(name, h_opt);

	lib_name = get_output_file("part");
	c_opt = init_clean_opt();
	c_opt->kmer = 15;
	c_opt->lib_name = lib_name;
	pe_clean_core(name, c_opt);

	solid_name = get_output_file("part.solid");
	pht = pe_load_hash(name);
	solid_reads = load_solid_reads(solid_name, pht->seqs, pht->n_seqs);

	n_per_threads = solid_reads->len / n_threads;
	single_edges = g_ptr_array_sized_new(64);
	run_threads(single_edges, solid_reads, pht, n_paired_reads, &n_single, 0,
			solid_reads->len, n_per_threads, 0.9);

	show_msg(__func__, "Updating single edges... \n");
	show_msg(__func__, "Before: %.2f sec\n", (float) (clock() - t)
			/ CLOCKS_PER_SEC);
	for (i = 0; i < single_edges->len; i++) {
		eg = g_ptr_array_index(single_edges, i);
		if (eg->alive) {
			//show_msg(__func__, "Updating edge [%d, %d]: %d ", eg->id,
			//		eg->len, eg->reads->len);
			upd_reads_by_ht(ht, eg, MISMATCHES, stage);
			//show_msg(__func__, " to %d \n", eg->reads->len);
			g_ptr_array_add(all_edges, eg);
		} else {
			g_ptr_array_remove_index_fast(single_edges, i);
			i--;
		}
	}
	show_msg(__func__, "After: %.2f sec\n", (float) (clock() - t)
			/ CLOCKS_PER_SEC);

	free(solid_name);
	free(name);
	free(lib_name);
	free(h_opt);
	free(c_opt);
	destroy_ht(pht);
	g_ptr_array_free(solid_reads, TRUE);
	g_ptr_array_free(single_edges, TRUE);
}

void pe_lib_core(int n_max_pairs, char *lib_file, char *solid_file) {
	clock_t t = clock();
	hash_table *ht = NULL;
	bwa_seq_t *seqs = NULL;
	FILE *solid = NULL;
	int i = 0, n_paired_reads = 0, n_per_threads = 0, n_single_reads = 0;
	GPtrArray *all_edges = NULL;
	readarray *solid_reads = NULL;
	edge *eg = NULL;

	solid = xopen(solid_file, "r");
	show_msg(__func__, "Library: %s \n", lib_file);
	show_msg(__func__, "Solid Reads: %s \n", solid_file);

	all_edges = g_ptr_array_sized_new(BUFSIZ);
	ht = pe_load_hash(lib_file);
	seqs = &ht->seqs[0];
	solid_reads = load_solid_reads(solid_file, seqs, ht->n_seqs);
	show_msg(__func__, "Solid reads loaded: %.2f sec\n", (float) (clock() - t)
			/ CLOCKS_PER_SEC);

	n_per_threads = solid_reads->len / 2 / n_threads;
	show_msg(__func__, "========================================== \n");
	show_msg(__func__, "Stage 1/2: Trying to use up the paired reads... \n");
	run_threads(all_edges, solid_reads, ht, &n_paired_reads, &n_single_reads,
			0, solid_reads->len / 2, n_per_threads, 0.95);
	n_paired_reads = 0;
	n_single_reads = 0;
	for (i = 0; i < all_edges->len; i++) {
		eg = g_ptr_array_index(all_edges, i);
		if (eg->alive) {
			n_paired_reads += eg->pairs->len;
			n_single_reads += eg->reads->len;
		}
	}
	show_msg(__func__, "Total reads: [%d=>%d]/%d \n", n_single_reads,
			n_paired_reads, ht->n_seqs);
	show_msg(__func__, "Stage 1 finished: %.2f sec\n", (float) (clock() - t)
			/ CLOCKS_PER_SEC);
	reset_edge_ids(all_edges);

	show_msg(__func__, "========================================== \n");
	show_msg(__func__, "Stage 2/2: Trying to assembly unpaired reads... \n");
	pe_lib_single(ht, all_edges, &n_paired_reads, &n_single_reads);
	show_msg(__func__, "Stage 2 finished: %.2f sec\n", (float) (clock() - t)
			/ CLOCKS_PER_SEC);

	post_process_edges(ht, all_edges);
	g_ptr_array_free(solid_reads, TRUE);
	destroy_ht(ht);
}

int pe_lib_usage() {
	show_msg(__func__,
			"Command: ./peta pair -p MAX_PAIRS read_file starting_reads \n");
	return 1;
}

void test_smith_waterman(char *lib_file) {
	hash_table *ht = NULL;
	bwa_seq_t *seqs = NULL, *read_1 = NULL, *read_2 = NULL;
	int score = 0;

	ht = pe_load_hash(lib_file);
	seqs = ht->seqs;

	read_1 = &seqs[1646764];
	read_2 = &seqs[2299221];
	score = smith_waterman(read_1, read_2, 2, -1, -2, 0);
	p_query(__func__, read_1);
	p_query(__func__, read_2);
	show_debug_msg(__func__, "Similarity score: %d \n", score);
}

int pe_lib(int argc, char *argv[]) {
	clock_t t = clock();
	int c = 0, n_max_pairs = 0;
	while ((c = getopt(argc, argv, "p:k:m:d:o:t:")) >= 0) {
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
		case 'd':
			sd_insert_size = atoi(optarg);
			break;
		case 't':
			n_threads = atoi(optarg);
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
	show_msg(__func__, "Done: %.2f sec\n", (float) (clock() - t)
			/ CLOCKS_PER_SEC);
	return 0;
}
