#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <glib.h>
#include <math.h>
#include "rand.h"
#include "ass.h"
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

int pair_ctg_id = 1;
int overlap_len = 0;
int insert_size = 0;
int sd_insert_size = 0;

void pe_lib_help() {
	show_msg(__func__, "--------------------------------------------------");
	show_msg(__func__, "peta pelib <library>");
	show_msg(__func__, "--------------------------------------------------");
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
 * Return 1:
 * The mate of 'read' should be in the used area;
 * Return 0:
 * The mate of 'read' should be in the future
 */
int is_paired(const edge *eg, const bwa_seq_t *read, const int ori) {
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
 * Keep only those reads with its mates at the right direction and position
 *
 * ori = 0: For below, 'read' should be right, 'mate' should be left
 *   ---------------------->
 *          ---->        ---->
 *          mate         read
 */
void keep_mates_in_pool(edge *eg, pool *cur_pool, int *next,
		const hash_table *ht, const int ori, const int remove_single) {
	bwa_seq_t *read = NULL, *mate = NULL;
	int i = 0, to_remove = 1;
	for (i = 0; i < cur_pool->reads->len; i++) {
		read = g_ptr_array_index(cur_pool->reads, i);
		mate = get_mate(read, ht->seqs);
		to_remove = 0;
		// If the mate of 'read' should be used already
		if (is_paired(eg, read, ori)) {
			// If the mate is not used, remove 'read'
			if ((read->rev_com != mate->rev_com) || (mate->status != USED)
					|| (eg->id - mate->contig_id) > 1) {
				//p_query("MATE NOT USED", read);
				//p_query(__func__, mate);
				to_remove = 1;
			} else {// If the distance is not within range, remove 'read'
				if (abs(eg->len - mate->shift) > (insert_size + sd_insert_size
						* SD_TIMES)) {
					//p_query("DISTANCE OUT OF RANGE", read);
					//p_query(__func__, mate);
					to_remove = 1;
				}
			}
		} else {
			to_remove = remove_single;
		}

		if (to_remove) {
			if (mate_pool_rm_index(cur_pool, i))
				i--;
		}
	}
	reset_c(next, NULL); // Reset the counter
	for (i = 0; i < cur_pool->n; i++) {
		read = g_ptr_array_index(cur_pool->reads, i);
		if (read->rev_com)
			check_c(next, read->rseq[read->cursor]);
		else
			check_c(next, read->seq[read->cursor]);
	}
}

pool *get_mate_pool_from_edge(edge *eg, hash_table *ht) {
	int i = 0;
	bwa_seq_t *r = NULL, *mate = NULL, *seqs = NULL;
	pool *mate_pool = NULL;
	mate_pool = new_pool();
	seqs = ht->seqs;
	for (i = 0; i < eg->reads->len; i++) {
		r = g_ptr_array_index(eg->reads, i);
		mate = get_mate(r, seqs);
		if (mate->status == FRESH && mate->is_in_m_pool == 0) {
			mate->rev_com = r->rev_com;
			mate_pool_add(mate_pool, mate);
		}
	}
	return mate_pool;
}

void add_mates_by_ol(bwa_seq_t *seqs, edge *eg, pool *cur_pool,
		pool *mate_pool, const int ol, const int nm, bwa_seq_t *query,
		const int ori) {
	int i = 0, j = 0, overlapped = 0, read_mate_ol = 0, go_to_add = 0;
	bwa_seq_t *mate = NULL, *tmp = NULL, *s = NULL;
	// Add the mate reads which overlap with the tail into the current pool
	for (i = 0; i < mate_pool->n; i++) {
		mate = g_ptr_array_index(mate_pool->reads, i);
		s = get_mate(mate, seqs);
		if (mate->is_in_c_pool || mate->status == USED || (abs(eg->len
				- s->shift) > (insert_size + sd_insert_size * SD_TIMES)))
			continue;
		tmp = mate;
		if (mate->rev_com)
			tmp = new_mem_rev_seq(mate, mate->len, 0);

		overlapped = ori ? find_ol(tmp, query, nm) : find_ol(query, tmp, nm);
		//		if (strcmp(mate->name, "6334942") == 0) {
		//			show_debug_msg("ORI", "ORI: %d \n", ori);
		//			p_ctg_seq("QUERY", query);
		//			p_query("MATE", tmp);
		//			p_query("ORIG", mate);
		//			p_query("USED", get_mate(mate, seqs));
		//			show_debug_msg(__func__, "OVERLAP 1: %d \n", overlapped);
		//			show_debug_msg(__func__, "OVERLAP 2: %d \n", read_mate_ol);
		//		}
		if (overlapped >= ol) {
			// Only if this mate overlaps with some read in the cur_pool, add it.
			// It is important because sometimes it maybe added just for coincidence.
			go_to_add = 1; // If there is no read in current pool, go to add the mate
			// If some read in current pool, we require that at least some read overlaps with the mate
			for (j = 0; j < cur_pool->reads->len; j++) {
				go_to_add = 0;
				s = g_ptr_array_index(cur_pool->reads, j);
				read_mate_ol = ori ? find_ol(mate, s, nm)
						: find_ol(s, mate, nm);
				if (read_mate_ol >= ol) {
					go_to_add = 1;
					break;
				}
			}
			if (go_to_add) {
				mate->cursor = ori ? (mate->len - overlapped - 1) : overlapped;
				pool_add(cur_pool, mate);
				mate_pool_rm_index(mate_pool, i);
				i--;
			}
		}
		if (mate->rev_com)
			bwa_free_read_seq(1, tmp);
	}
}

void maintain_pool(alignarray *aligns, const hash_table *ht, pool *cur_pool,
		pool *mate_pool, edge *ass_eg, bwa_seq_t *query, int *next, int ori) {
	alg *a = NULL;
	int i = 0, index = 0;
	bwa_seq_t *s = NULL, *mate = NULL, *seqs = ht->seqs;
	for (i = 0; i < aligns->len; i++) {
		a = g_ptr_array_index(aligns, i);
		index = a->r_id;
		if (index >= ht->n_seqs)
			continue;
		s = &seqs[index];
		s->rev_com = a->rev_comp;
		mate = get_mate(s, seqs);
		mate->rev_com = s->rev_com;

		if (s->is_in_c_pool || (a->pos + query->len - 1) > s->len || s->status
				== USED || s->contig_id == ass_eg->id)
			continue;
		if (s->rev_com)
			s->cursor = ori ? (s->len - query->len - 1 - a->pos) : (s->len
					- a->pos);
		else
			s->cursor = ori ? (a->pos - 1) : (a->pos + query->len);

		if (s->cursor >= s->len || s->cursor < 0) {
			s->cursor = 0;
			continue;
		}
		pool_add(cur_pool, s);
		if (mate->status != USED && !mate->is_in_c_pool && !mate->is_in_m_pool) {
			mate_pool_add(mate_pool, mate);
		}
	}
	rm_partial(cur_pool, ori, query, 2);
	add_mates_by_ol(seqs, ass_eg, cur_pool, mate_pool, MATE_OVERLAP_THRE,
			MISMATCHES, query, ori);
	// Keep only the reads whose mate is used previously.
	if (ass_eg->len >= (insert_size + sd_insert_size * SD_TIMES))
		keep_mates_in_pool(ass_eg, cur_pool, next, ht, ori, 0);
	for (i = 0; i < cur_pool->n; i++) {
		s = g_ptr_array_index(cur_pool->reads, i);
		if (s->rev_com)
			check_c(next, s->rseq[s->cursor]);
		else
			check_c(next, s->seq[s->cursor]);
	}
}

/**
 * If return NULL, it indicates there are too many (>=4 by default) counter pairs.
 * In this case, the caller should not continue to do single extension.
 */
pool *get_start_pool(const hash_table *ht, const bwa_seq_t *init_read,
		const int ori) {
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
	p_query(__func__, query);

	pe_aln_query(query, query->seq, ht, 4, overlap_len, 0, alns);
	pe_aln_query(query, query->rseq, ht, 4, overlap_len, 1, alns);
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
		if (s->rev_com)
			s->cursor = ori ? (0 - a->pos - 1) : (s->len - a->pos);
		else
			s->cursor = ori ? (a->pos - 1) : (a->pos + s->len);
		if (s->contig_id == INVALID_CONTIG_ID || s->status == USED || s->cursor
				< -1 || s->cursor > s->len || s->is_in_c_pool) {
			s->cursor = 0;
			continue;
		}
		s->rev_com = a->rev_comp;
		mate->rev_com = a->rev_comp;
		pool_add(init_pool, s);
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

edge *pair_extension(edge *pre_eg, const hash_table *ht, const bwa_seq_t *s,
		const int ori) {
	bwa_seq_t *query = NULL;
	pool *cur_pool = NULL, *mate_pool = NULL;
	int *c = NULL;
	int *next = (int*) calloc(5, sizeof(int));
	alignarray *aligns = NULL;
	edge *eg = NULL;

	if (pre_eg) {
		eg = pre_eg;
	} else {
		eg = new_eg();
		eg->id = pair_ctg_id;
		eg->contig = new_seq(s, s->len, 0);
		eg->len = s->len;
	}

	if (ori)
		query = new_seq(s, overlap_len, 0);
	else
		query = new_seq(s, overlap_len, s->len - overlap_len);
	cur_pool = get_start_pool(ht, s, ori);
	mate_pool = pre_eg ? get_mate_pool_from_edge(pre_eg, ht)
			: get_init_mate_pool(cur_pool, ht->seqs);
	aligns = g_ptr_array_sized_new(N_DEFAULT_ALIGNS);
	if (ori)
		seq_reverse(eg->len, eg->contig->seq, 0);
	while (1) {
		//		p_query(__func__, query);
		reset_c(next, NULL); // Reset the counter
		if (!cur_pool || is_repetitive_q(query)) {
			show_msg(__func__, "[%d, %d] Repetitive pattern, stop!\n", eg->id,
					eg->len);
			p_query(__func__, query);
			break;
		}
		//p_query(__func__, query);
		pe_aln_query(query, query->seq, ht, MISMATCHES, query->len, 0, aligns);
		pe_aln_query(query, query->rseq, ht, MISMATCHES, query->len, 1, aligns);
		//p_align(aligns);
		maintain_pool(aligns, ht, cur_pool, mate_pool, eg, query, next, ori);
		// If no read in current pool, try less stringent overlapped length and zero mismatches
		if (cur_pool->reads->len == 0) {
			add_mates_by_ol(ht->seqs, eg, cur_pool, mate_pool,
					RELAX_MATE_OL_THRE, 0, query, ori);
		}
		reset_alg(aligns);
		//p_ctg_seq("Contig", eg->contig);
		//p_pool("Current Pool", cur_pool, next);
		//p_pool("Mate Pool", mate_pool, next);
		c = get_abs_most(next, STRICT_PERC);
		// show_debug_msg(__func__, "Next char: %d \n", c[0]);
		if (cur_pool->n <= 0) {
			show_msg(__func__, "[%d, %d] No hits, stop here. \n", eg->id,
					eg->len);
			show_debug_msg(__func__, "[%d, %d] No hits, stop here. \n", eg->id,
					eg->len);
			break;
		}
		if (c[0] == -1) {
			if (bases_sup_branches(cur_pool, ori, BASES_SUPPORT_THRE)) {
				// If there are multiple branches, keep those reads whose mate is used only.
				// Normally, this method is only called after the edge length is longer than
				// insert_size + sd_insert_size * 4. But here we keep those pairs whose distance
				// may be shorter than it.
				keep_mates_in_pool(eg, cur_pool, next, ht, ori, 1);
				//	p_pool("Current Pool After Removing", cur_pool, next);
				c = get_abs_most(next, STRICT_PERC);
				// If there are still multiple branches, stop!
				if (c[0] == -1) {
					show_debug_msg(__func__, "[%d]: length %d, reads %d \n",
							eg->id, eg->len, eg->reads->len);
					show_msg(__func__,
							"[%d, %d] Probable branching, stop here. \n",
							eg->id, eg->len);
					show_debug_msg(__func__,
							"[%d, %d] Probable branching, stop here. \n",
							eg->id, eg->len);
					break;
				}
			}
			c[0] = get_pure_most(next);
		}
		forward(cur_pool, c[0], eg, ori);
		ext_con(eg->contig, c[0], 0);
		//		show_debug_msg(__func__, "Edge %d, length %d \n", eg->id, eg->len);
		//		p_ctg_seq("Contig now", eg->contig);
		eg->len = eg->contig->len;
		ext_que(query, c[0], ori);
		free(c);
		c = NULL;
		if (eg->len % 50 == 0) {
			show_debug_msg(__func__, "Assembling... [%d, %d] \n", eg->id,
					eg->len);
			clean_mate_pool(mate_pool);
		}
	}
	if (ori)
		seq_reverse(eg->len, eg->contig->seq, 0);
	free_alg(aligns);
	free(next);
	free(c);
	free_pool(cur_pool);
	free_mate_pool(mate_pool);
	bwa_free_read_seq(1, query);
	return eg;
}

edge* pe_ext(hash_table *ht, bwa_seq_t *query) {
	edge *eg = NULL;
	pool *init_pool = NULL;
	bwa_seq_t *second_round_q = NULL;

	init_pool = get_start_pool(ht, query, 0);
	if (init_pool->n == 0 || bases_sup_branches(init_pool, 0,
			STRICT_BASES_SUP_THRE)) {
		show_msg(__func__, "Read %s may be in junction area, skip. \n",
				query->name);
		query->status = TRIED;
		free_pool(init_pool);
		return NULL;
	}
	free_pool(init_pool);

	show_debug_msg(__func__, "Extending to the right... \n");
	eg = pair_extension(NULL, ht, query, 0);
	show_debug_msg(__func__, "Edge after right extension: [%d, %d]\n", eg->id,
			eg->len);
	show_debug_msg(__func__, "Got edge [%d, %d]. Extending to the left... \n",
			eg->id, eg->len);
	pair_extension(eg, ht, query, 1);
	show_debug_msg(__func__, "Edge after left extension: [%d, %d]\n", eg->id,
			eg->len);

	show_debug_msg(__func__, "Second round: extending to the right... \n");
	second_round_q = new_seq(eg->contig, query->len, eg->len - query->len);
	pair_extension(eg, ht, second_round_q, 0);
	bwa_free_read_seq(1, second_round_q);

	show_debug_msg(__func__, "Second round: extending to the left... \n");
	second_round_q = new_seq(eg->contig, query->len, 0);
	pair_extension(eg, ht, second_round_q, 1);
	upd_reads(eg, MISMATCHES);
	//p_readarray(eg->reads, 1);
	return eg;
}

void keep_pairs_only(edge *eg, bwa_seq_t *seqs) {
	int i = 0, min_shift = 0, max_shift = eg->len, read_len = 0;
	bwa_seq_t *read = NULL, *mate = NULL;
	show_debug_msg(__func__, "Getting mate pairs...\n");
	p_readarray(eg->reads, 1);
	for (i = 0; i < eg->reads->len; i++) {
		read = g_ptr_array_index(eg->reads, i);
		read_len = read->len;
		mate = get_mate(read, seqs);
		if (mate->status == USED && mate->contig_id == eg->id) {
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
	if (max_shift < eg->len) {
		show_debug_msg(__func__, "Truncate edge [%d, %d] to [%d, %d] \n",
				eg->id, eg->len, min_shift, max_shift);
		max_shift += read_len;
		if (max_shift < eg->len) {
			eg->len = max_shift;
			eg->contig->seq[eg->len] = '\0';
			eg->contig->len = eg->len;
		}
		if (min_shift > 0) {
			trun_seq(eg->contig, min_shift);
			eg->len -= min_shift;
		}
	}
}

/**
 * order:
 * 	1: eg_1 is left
 *  -1: eg_1 is right
 */
int est_pair_gap(edge *eg_1, edge *eg_2, int order) {
	GPtrArray *paired_reads = NULL;
	int i = 0;
	bwa_seq_t *read = NULL, *mate = NULL;
	int len = 0, gap_sum = 0, gap_ave = -1, n_counted = 0;
	paired_reads = find_unconditional_paired_reads(eg_1->reads, eg_2->reads);
	if (paired_reads->len >= MIN_ORDER_PAIRS)
		for (i = 0; i < paired_reads->len; i += 2) {
			read = g_ptr_array_index(paired_reads, i);
			mate = g_ptr_array_index(paired_reads, i + 1);
			if (read->rev_com != mate->rev_com)
				continue;
			n_counted++;
			if (order == 1) {
				len = (eg_1->len - read->shift) + mate->shift;
			} else
				len = (eg_2->len - mate->shift) + read->shift;
			gap_sum += insert_size - len;
		}
	if (n_counted > 0) {
		gap_ave = gap_sum / n_counted;
		gap_ave = gap_ave > 0 ? gap_ave : 0;
	}
	g_ptr_array_free(paired_reads, TRUE);
	return gap_ave;
}

/**
 * Return:
 * 	0: They cannot be put together
 * 	1: eg_1 is left, eg_2 is right
 *  -1: eg_1 is right, eg_2 is left
 */
int order_two_edges(edge *eg_1, edge *eg_2) {
	GPtrArray *paired_reads = NULL;
	int order = 0;
	int i = 0;
	bwa_seq_t *read = NULL, *mate = NULL;
	paired_reads = find_unconditional_paired_reads(eg_1->reads, eg_2->reads);
	if (paired_reads->len < MIN_ORDER_PAIRS)
		order = 0;
	else {
		for (i = 0; i < paired_reads->len; i += 2) {
			read = g_ptr_array_index(paired_reads, i);
			mate = g_ptr_array_index(paired_reads, i + 1);
			if (read->rev_com != mate->rev_com)
				continue;
			if (is_left_mate(read->name)) {
				order = (read->rev_com) ? -1 : 1;
			} else {
				order = (read->rev_com) ? 1 : -1;
			}
			break;
		}
	}
	g_ptr_array_free(paired_reads, TRUE);
	return order;
}

GPtrArray *scaffolding(GPtrArray *single_edges) {
	edge *eg_i = NULL, *eg_j = NULL;
	int i = 0, j = 0;
	readarray *paired = NULL;
	GPtrArray *paths = NULL;
	rm_path *p = NULL;

	paths = g_ptr_array_sized_new(1023);
	for (i = 0; i < single_edges->len; i++) {
		eg_i = g_ptr_array_index(single_edges, i);
		p = get_single_edge_path(eg_i);
	}
}

void pe_lib_core(int n_max_pairs, char *lib_file, char *solid_file) {
	hash_table *ht = NULL;
	bwa_seq_t *query = NULL, *seqs = NULL, *r = NULL;
	FILE *solid = xopen(solid_file, "r");
	char line[80];
	int index = 0, ol = 0, n_pairs = 0, n_part_pairs = 0, n_total_reads = 0;
	int s_index = 0, e_index = 0, counter = -1, i = 0, line_no = 0;
	edge *eg = NULL;
	double *pairs = NULL, *partial_pairs = NULL, mean_ins_size = 0,
			sd_ins_size = 0;
	GPtrArray *all_edges = NULL;
	FILE *pair_contigs = xopen("pair_contigs.fa", "w");
	show_msg(__func__, "Library: %s \n", lib_file);
	show_msg(__func__, "Solid Reads: %s \n", solid_file);
	show_msg(__func__, "Maximum Pairs: %d \n", n_max_pairs);

	all_edges = g_ptr_array_sized_new(BUFSIZ);
	ht = pe_load_hash(lib_file);
	seqs = &ht->seqs[0];
	ol = seqs->len / 2; // Read length
	s_index = 1000;
	e_index = 1002;
	pairs = (double*) calloc(n_max_pairs + 1, sizeof(double));
	while (fgets(line, 80, solid) != NULL && n_total_reads < ht->n_seqs * 0.7) {
		line_no++;
		index = atoi(line);
		query = &ht->seqs[index];
		if (query->status != FRESH)
			continue;
		counter++;
		//		if (counter <= s_index)
		//			continue;
		//		if (counter > e_index)
		//			break;
		show_msg(__func__,
				"---------- [%d] Processing read %d: %s ----------\n", line_no,
				counter, query->name);
		show_debug_msg(__func__,
				"---------- [%d] Processing read %d: %s ----------\n", counter,
				line_no, query->name);
		eg = pe_ext(ht, query);
		//		g_ptr_array_sort(eg->reads, (GCompareFunc) cmp_reads_by_name);
		//		partial_pairs = get_pairs_on_edge(eg, &n_part_pairs);
		//		if (n_part_pairs >= PAIRS_PER_EDGE)
		//			n_part_pairs = PAIRS_PER_EDGE;
		//		if (n_part_pairs + n_pairs >= n_max_pairs) {
		//			break;
		//		}
		//		concat_doubles(pairs, &n_pairs, partial_pairs, n_part_pairs);
		//		free(partial_pairs);
		//		partial_pairs = NULL;
		//		p_ctg_seq("Contig now", eg->contig);
		//		show_msg(__func__, "reads of contig %d: [%d, %d] \n", eg->id, eg->len,
		//				eg->reads->len);
		//		show_msg(__func__, "+%d, %d pairs now \n", n_part_pairs, n_pairs);
		//		n_part_pairs = 0;
		//		//p_readarray(eg->reads, 1);
		//		log_edge(eg);
		//		destroy_eg(eg);
		//		eg = NULL;

		if (eg) {
			keep_pairs_only(eg, ht->seqs);
			//p_readarray(eg->pairs, 1);
			if (eg->len < 100 || eg->pairs->len <= 2) {
				show_msg(
						__func__,
						"ABANDONED [%d]: length %d, reads %d=>%d. Total reads %d/%d \n",
						eg->id, eg->len, eg->reads->len, eg->pairs->len,
						n_total_reads, ht->n_seqs);
				clear_used_reads(eg, 1);
				destroy_eg(eg);
				eg = NULL;
			} else {
				n_total_reads += eg->pairs->len;
				g_ptr_array_add(all_edges, eg);
				show_msg(__func__,
						"[%d]: length %d, reads %d=>%d. Total reads %d/%d \n",
						eg->id, eg->len, eg->reads->len, eg->pairs->len,
						n_total_reads, ht->n_seqs);
			}
			pair_ctg_id++;
		}
		//break;
	}
	save_edges(all_edges, pair_contigs, 0, 0, 100);
	free(partial_pairs);
	destroy_eg(eg);
	//	mean_ins_size = mean(pairs, n_pairs);
	//	sd_ins_size = std_dev(pairs, n_pairs);
	//	show_msg(__func__, "Mean Insert Size: %f \n", mean_ins_size);
	//	show_msg(__func__, "Standard Deviation of Insert Size: %f \n", sd_ins_size);
	free(pairs);
	fclose(solid);
	fclose(pair_contigs);
	destroy_ht(ht);
}

int pe_lib_usage() {
	show_msg(__func__,
			"Command: ./peta pair -p MAX_PAIRS read_file starting_reads \n");
	return 1;
}

int pe_lib(int argc, char *argv[]) {
	clock_t t = clock();
	int c = 0, n_max_pairs = 0;
	while ((c = getopt(argc, argv, "p:o:m:d:")) >= 0) {
		switch (c) {
		case 'p':
			n_max_pairs = atoi(optarg);
			break;
		case 'o':
			overlap_len = atoi(optarg);
			break;
		case 'm':
			insert_size = atoi(optarg);
			break;
		case 'd':
			sd_insert_size = atoi(optarg);
			break;
		}
	}
	if (optind + 2 > argc || n_max_pairs == 0) {
		return pe_lib_usage();
	}
	pe_lib_core(n_max_pairs, argv[optind], argv[optind + 1]);
	show_msg(__func__, "Done: %.2f sec\n", (float) (clock() - t)
			/ CLOCKS_PER_SEC);
	return 0;
}
