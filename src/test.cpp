#include <unordered_map>
#include <stdint.h>
#include <glib.h>
#include <cstdio>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>
#include "bwtaln.h"
#include "pechar.h"
#include "utils.h"
#include "peseq.h"
#include "rnaseq.h"
#include "tpl.hpp"
#include "junction.hpp"
#include "pool.hpp"
#include "hash.hpp"
#include "test.hpp"

using namespace std;

char *SPOMBE_TX =
		"/home/carl/Projects/peta_pair/scripts/spombe.broad.tx.fasta.rev";
int TESTING_CASES = 1;
char *ID_LIST =
		"/home/carl/Projects/peta_pair/SRR097897_half/idba.fa.full.only";
GPtrArray *MISSING_IDS = NULL;

void test_init() {
	if (!TESTING_CASES)
		return;
	read_ids(ID_LIST);
}

char *get_test_fa(int id, int suffix) {
	char *name = (char*) calloc(128, sizeof(char));
	sprintf(name, "test/%d.%d.fa", id, suffix);
	return name;
}

char *test_save_tpl(tpl *t, int suffix) {
	char *name = get_test_fa(t->id, suffix);
	FILE *fa = xopen(name, "w");
	char *h = (char*) malloc(BUFSIZE);
	sprintf(h, ">%d.%d\n", t->id, suffix);
	save_con(h, t->ctg, fa);
	fflush(fa);
	fclose(fa);
}

char *get_test_psl(int id, int suffix) {
	char *name = (char*) calloc(128, sizeof(char));
	sprintf(name, "test/%d.%d.psl", id, suffix);
	return name;
}

void blat_tpl_seq(tpl *t, int suffix) {
	if (!TESTING_CASES)
		return;
	char *fa_n = get_test_fa(t->id, suffix);
	test_save_tpl(t, suffix);
	char *psl_n = get_test_psl(t->id, suffix);
	char *cmd = (char*) malloc(sizeof(char) * 1024);
	sprintf(cmd, "blat %s %s %s > /dev/null", SPOMBE_TX, fa_n, psl_n);
	//show_msg(__func__, "%s\n", cmd);
	system(cmd);
	free(fa_n);
	free(cmd);
}

GPtrArray *read_ids(char *ids_fn) {
	char buf[1000];
	char *id = NULL;
	MISSING_IDS = g_ptr_array_sized_new(32);
	FILE *f = xopen(ids_fn, "r");
	while (fgets(buf, sizeof(buf), f)) {
		trim(buf);
		id = strdup(buf);
		//show_debug_msg("test_read_ids", "%s\n", id);
		g_ptr_array_add(MISSING_IDS, id);
	}
	fclose(f);
	return MISSING_IDS;
}

GPtrArray *read_blat_hits(char *blat_psl) {
	char buf[1000];
	int line_no = 0, i = 0;
	char *attr[32], *intstr[32];
	FILE *psl_fp = xopen(blat_psl, "r");
	GPtrArray *hits = g_ptr_array_sized_new(32);
	blat_hit *h = NULL;
	while (fgets(buf, sizeof(buf), psl_fp)) {
		line_no++;
		if (line_no <= 5)
			continue;
		i = 0;
		attr[0] = strtok(buf, "\t");
		while (attr[i] != NULL) { //ensure a pointer was found
			attr[++i] = strtok(NULL, "\t"); //continue to tokenize the string
		}
		if (i < 20)
			break;
		h = (blat_hit*) malloc(sizeof(blat_hit));
		h->n_match = atoi(attr[0]);
		h->n_mismatch = atoi(attr[1]);
		h->n_rep = atoi(attr[2]);
		h->n_n = atoi(attr[3]);
		h->n_query_gap = atoi(attr[4]);
		h->n_query_gap_base = atoi(attr[5]);
		h->n_ref_gap = atoi(attr[6]);
		h->n_ref_gap_base = atoi(attr[7]);
		h->strand = attr[8][0];
		h->query = strdup(attr[9]);
		h->q_len = atoi(attr[10]);
		h->q_start = atoi(attr[11]);
		h->q_end = atoi(attr[12]);
		h->ref = strdup(attr[13]);
		h->r_len = atoi(attr[14]);
		h->r_start = atoi(attr[15]);
		h->r_end = atoi(attr[16]);
		h->n_block = atoi(attr[17]);
		h->block_size = (int*) calloc(h->n_block, sizeof(int));
		h->block_size[0] = atoi(strtok(attr[18], ","));
		i = 0;
		while (i < h->n_block - 1) {
			h->block_size[++i] = atoi(strtok(NULL, ",")); //continue to tokenize the string
		}
		h->query_block_start = (int*) calloc(h->n_block, sizeof(int));
		h->query_block_start[0] = atoi(strtok(attr[19], ","));
		i = 0;
		while (i < h->n_block - 1) {
			h->query_block_start[++i] = atoi(strtok(NULL, ",")); //continue to tokenize the string
		}
		h->ref_block_start = (int*) calloc(h->n_block, sizeof(int));
		h->ref_block_start[0] = atoi(strtok(attr[20], ","));
		i = 0;
		while (i < h->n_block - 1) {
			h->ref_block_start[++i] = atoi(strtok(NULL, ",")); //continue to tokenize the string
		}
		g_ptr_array_add(hits, h);
	}
	fclose(psl_fp);
	return hits;
}

void test_print_blat_hit(blat_hit *h) {
	int i = 0;
	printf("[test_steps]\t");
	printf("%d\t", h->n_match);
	printf("%d\t", h->n_mismatch);
	printf("%d\t", h->n_rep);
	printf("%d\t", h->n_n);
	printf("%d\t", h->n_query_gap);
	printf("%d\t", h->n_query_gap_base);
	printf("%d\t", h->n_ref_gap);
	printf("%d\t", h->n_ref_gap_base);
	printf("%c\t", h->strand);
	printf("%s\t", h->query);
	printf("%d\t", h->q_len);
	printf("%d\t", h->q_start);
	printf("%d\t", h->q_end);
	printf("%s\t", h->ref);
	printf("%d\t", h->r_len);
	printf("%d\t", h->r_start);
	printf("%d\t", h->r_end);
	printf("%d\t", h->n_block);
	for (i = 0; i < h->n_block; i++) {
		printf("%d,", h->block_size[i]);
	}
	printf("\t");
	for (i = 0; i < h->n_block; i++) {
		printf("%d,", h->query_block_start[i]);
	}
	printf("\t");
	for (i = 0; i < h->n_block; i++) {
		printf("%d,", h->ref_block_start[i]);
	}
	printf("\t\n");
}

int count_bad_bases(blat_hit *h) {
	return h->n_mismatch + h->n_query_gap_base + h->n_ref_gap_base;
}

int test_is_full_length(blat_hit *h) {
	float similarity = float(h->n_match) / float(h->r_len);
	if (similarity >= 0.99 && count_bad_bases(h) <= 10)
		return 1;
	return 0;
}

void test_print_msg(const char *header, const char *fmt, ...) {
	if (!TESTING_CASES)
		return;
	va_list args;
	va_start(args, fmt);
	printf("[%s] ", header);
	vprintf(fmt, args);
	va_end(args);
}

void test_print_tag(char *msg) {
	if (!TESTING_CASES)
		return;
	show_debug_msg(__func__, "------ %s ------\n", msg);
}

int test_to_check(tpl *t, int suffix) {
	if (!TESTING_CASES)
		return 0;
	blat_tpl_seq(t, suffix);
	int to_check = 0;
	char *psl_n = get_test_psl(t->id, suffix);
	GPtrArray *hits = read_blat_hits(psl_n);
	blat_hit *h = NULL;
	char *id = NULL;
	int i = 0, j = 0;
	for (i = 0; i < hits->len; i++) {
		h = (blat_hit*) g_ptr_array_index(hits, i);
		for (j = 0; j < MISSING_IDS->len; j++) {
			id = (char*) g_ptr_array_index(MISSING_IDS, j);
			if (strcmp(id, h->ref) == 0) {
				to_check = 1;
				break;
			}
		}
	}
	return to_check;
}

int test_check_missing(tpl *t, int suffix) {
	if (!TESTING_CASES)
		return 0;
	blat_tpl_seq(t, suffix);
	char *psl_n = get_test_psl(t->id, suffix);
	int is_full_length = 0;

	GPtrArray *hits = read_blat_hits(psl_n);
	blat_hit *h = NULL;
	char *id = NULL;
	int i = 0, j = 0;
	if (hits->len == 0) {
		show_debug_msg(__func__, "No hits\n");
	}
	for (i = 0; i < hits->len; i++) {
		h = (blat_hit*) g_ptr_array_index(hits, i);
		for (j = 0; j < MISSING_IDS->len; j++) {
			id = (char*) g_ptr_array_index(MISSING_IDS, j);
			if (strcmp(id, h->ref) == 0) {
				test_print_blat_hit(h);
				show_debug_msg(__func__,
						"Template [%d, %d], Test %d: %s FULL-LENGTH\n", t->id,
						t->len, suffix, test_is_full_length(h) ? "IS " : "NOT");
				if (test_is_full_length(h))
					is_full_length = 1;
				break;
			}
		}
	}
	return is_full_length;
}

int test_align_tpl_seq(tpl *t, int suffix) {
	if (!TESTING_CASES)
		return 0;
	int is_full_length = 0;
	blat_tpl_seq(t, suffix);
	char *psl_n = get_test_psl(t->id, suffix);
	GPtrArray *hits = read_blat_hits(psl_n);
	blat_hit *h = NULL;
	int i = 0;
	for (i = 0; i < hits->len; i++) {
		h = (blat_hit*) g_ptr_array_index(hits, i);
		show_debug_msg(__func__,
				"Template [%d, %d], Test %d: %s FULL-LENGTH\n", t->id, t->len,
				suffix, test_is_full_length(h) ? "IS " : "NOT");
		if (test_is_full_length(h))
			is_full_length = 1;
		test_print_blat_hit(h);
	}
	if (hits->len == 0) {
		show_debug_msg(__func__, "Template [%d, %d], Test %d: NO HIT\n", t->id,
				t->len, suffix);
	}

	free(psl_n);
	return is_full_length;
}

void test_steps(tpl *t, int suffix, char *step) {
	if (!TESTING_CASES)
		return;
	blat_tpl_seq(t, suffix);
	char *psl_n = get_test_psl(t->id, suffix);
	GPtrArray *hits = read_blat_hits(psl_n);
	int i = 0, j = 0;
	blat_hit *h = NULL;
	char *id = NULL;
	for (i = 0; i < hits->len; i++) {
		h = (blat_hit*) g_ptr_array_index(hits, i);
		for (j = 0; j < MISSING_IDS->len; j++) {
			id = (char*) g_ptr_array_index(MISSING_IDS, j);
			if (strcmp(id, h->ref) == 0) {
				show_debug_msg(__func__,
						"[%s] %s:\t Template [%d, %d] %s \t %s \n", step, id,
						t->id, t->len, t->alive ? "ALIVE" : "DEAD",
						test_is_full_length(h) ? "FULL" : "PARTIAL");
				test_print_blat_hit(h);
			}
		}
	}
}

void validate_junctions(char *junc_fn, char *pair_fa, char *pair_psl,
		char *hash_fn) {
	int i = 0, j = 0, x = 0, has_hit = 0, has_valid_hit = 0;
	int left_len = 0, right_len = 0, half_len = 0, len = 0;
	tpl *t = NULL;
	junction *jun = NULL;
	GPtrArray *all_tpls = g_ptr_array_sized_new(32);
	GPtrArray *all_junctions = g_ptr_array_sized_new(32);
	GPtrArray *junc_seqs = g_ptr_array_sized_new(32);
	GPtrArray *paired_hits = NULL;
	blat_hit *h = NULL;
	bwa_seq_t *branch_part = NULL, *main_part = NULL, *junc_seq = NULL;

	paired_hits = read_blat_hits(pair_psl);

	read_juncs_from_file(junc_fn, pair_fa, all_tpls, all_junctions);
	hash_table *ht = load_k_hash(hash_fn);
	half_len = ht->o->read_len - JUNCTION_BOUNDARY_BASE;
	for (i = 0; i < all_tpls->len; i++) {
		t = (tpl*) g_ptr_array_index(all_tpls, i);
		t->reads = g_ptr_array_sized_new(64);
		refresh_tpl_reads(ht, t, N_MISMATCHES);
		g_ptr_array_sort(t->reads, (GCompareFunc) cmp_reads_by_name);
	}
	for (i = 0; i < all_tpls->len; i++) {
		t = (tpl*) g_ptr_array_index(all_tpls, i);
		printf("Inspecting template [%d, %d] ...\n", t->id, t->len);
		set_jun_reads(ht, t);
		if (t->m_juncs) {
			printf("\nMain template hits: \n");
			for (j = 0; j < paired_hits->len; j++) {
				h = (blat_hit*) g_ptr_array_index(paired_hits, j);
				if (atoi(h->query) == t->id) {
					test_print_blat_hit(h);
				}
			}
			printf("\nBranch template hits: \n");
			for (x = 0; x < t->m_juncs->len; x++) {
				jun = (junction*) g_ptr_array_index(t->m_juncs, x);
				p_junction(jun);
				has_hit = 0;
				has_valid_hit = 0;
				if (jun->ori) {
					len = min3(jun->branch_tpl->len, half_len, jun->locus);
					branch_part = new_seq(jun->branch_tpl->ctg, len,
							jun->branch_tpl->len - len);
					main_part = new_seq(t->ctg, len, jun->locus - len);
				} else {
					len = min3(jun->branch_tpl->len, half_len, t->len
							- jun->locus);
					branch_part = new_seq(jun->branch_tpl->ctg, len, 0);
					main_part = new_seq(t->ctg, len, jun->locus);
				}
				printf("Junction_pairs: %d \n", count_pairs(ht->seqs, jun));
				p_query("Branch_junc", branch_part);
				p_query("Main_junc", main_part);
				printf("Mismatches: %d\n", seq_ol(branch_part, main_part,
						main_part->len, main_part->len));
				//bwa_free_read_seq(1, branch_part);
				bwa_free_read_seq(1, main_part);

				if (jun->ori) {
					len = min3(jun->branch_tpl->len, half_len, t->len
							- jun->locus);
					main_part = new_seq(t->ctg, len, jun->locus);
					junc_seq = blank_seq(main_part->len + branch_part->len);
					memcpy(junc_seq->seq, branch_part->seq, branch_part->len);
					memcpy(junc_seq->seq + branch_part->len, main_part->seq,
							main_part->len);
				} else {
					len = min3(jun->branch_tpl->len, half_len, jun->locus);
					main_part = new_seq(t->ctg, len, jun->locus - len);
					junc_seq = blank_seq(main_part->len + branch_part->len);
					memcpy(junc_seq->seq, main_part->seq, main_part->len);
					memcpy(junc_seq->seq + main_part->len, branch_part->seq,
							branch_part->len);
				}
				junc_seq->len = branch_part->len + main_part->len;
				set_rev_com(junc_seq);
				junc_seq->name = (char*) malloc(sizeof(char) * 1024);
				sprintf(junc_seq->name, ">%d_%d_%d_%d\n", t->id,
						jun->branch_tpl->id, jun->locus, jun->ori);
				g_ptr_array_add(junc_seqs, junc_seq);

				for (j = 0; j < paired_hits->len; j++) {
					h = (blat_hit*) g_ptr_array_index(paired_hits, j);
					if (atoi(h->query) == jun->branch_tpl->id) {
						has_hit = 1;
						if (h->n_query_gap == 0 && h->n_ref_gap == 0
								&& h->n_mismatch < 5) {
							if (jun->ori == 0) {
								if (h->q_start <= 2) {
									has_valid_hit = 1;
								}
							} else {
								if (h->q_end >= h->q_len - 2) {
									has_valid_hit = 1;
								}
							}
						}
						test_print_blat_hit(h);
					}
				}
				if (!has_hit) {
					printf("VALID: No hit for branch template [%d, %d] \n",
							jun->branch_tpl->id, jun->branch_tpl->len);
				} else {
					if (has_valid_hit) {
						printf("VALID: The branch [%d, %d] is valid. \n",
								jun->branch_tpl->id, jun->branch_tpl->len);
					} else {
						printf("VALID: Wrong branch: [%d, %d] \n",
								jun->branch_tpl->id, jun->branch_tpl->len);
					}
				}
			}
			printf("\n+++\n");
		}
	}
	FILE *fp = xopen("../SRR097897_branch/junc.fa", "w");
	for (i = 0; i < junc_seqs->len; i++) {
		junc_seq = (bwa_seq_t*) g_ptr_array_index(junc_seqs, i);
		save_con(junc_seq->name, junc_seq, fp);
	}
	fclose(fp);
}

void read_juncs_from_file(char *junc_fn, char *pair_fa, GPtrArray *all_tpls,
		GPtrArray *all_junctions) {
	FILE *junc_fp = xopen(junc_fn, "r");
	bwa_seq_t *seqs = NULL, *ctg = NULL;
	uint64_t n_ctgs = 0, i = 0, id = 0;
	seqs = load_reads(pair_fa, &n_ctgs);
	tpl *t = NULL, *main_tpl = NULL, *branch = NULL;
	tpl_hash tpls;
	char buf[BUFSIZ];
	char *attr[18], *idstr[18];
	junction *jun = NULL;
	for (i = 0; i < n_ctgs; i++) {
		t = new_tpl();
		ctg = &seqs[i];
		t->id = atoi(ctg->name);
		t->ctg = new_seq(ctg, ctg->len, 0);
		id = t->id;

		tpls[t->id] = t;
		//show_debug_msg(__func__, "template %d \n", t->id);
		t->len = t->ctg->len;
		t->alive = 1;
		g_ptr_array_add(all_tpls, t);
	}

	int line = 0;
	while (fgets(buf, sizeof(buf), junc_fp)) {
		line++;
		if (line == 1)
			continue;
		i = 0;
		attr[0] = strtok(buf, "\t");
		while (attr[i] != NULL) { //ensure a pointer was found
			//			printf("fields[%d] = %s\n", i, fields[i]);
			attr[++i] = strtok(NULL, "\t"); //continue to tokenize the string
		}
		idstr[0] = strtok(attr[0], ", ");
		idstr[1] = strtok(NULL, ", ");
		i = 0;
		while (idstr[0][++i] != '\0') {
			idstr[0][i - 1] = idstr[0][i];
		}
		idstr[0][i - 1] = '\0';

		idstr[2] = strtok(attr[1], ", ");
		idstr[3] = strtok(NULL, ", ");
		i = 0;
		while (idstr[2][++i] != '\0') {
			idstr[2][i - 1] = idstr[2][i];
		}
		idstr[2][i - 1] = '\0';
		//printf("%d\n", atoi(idstr[0]));
		//printf("%d\n", atoi(idstr[2]));

		main_tpl = (tpl*) tpls[atoi(idstr[0])];
		branch = (tpl*) tpls[atoi(idstr[2])];
		jun
				= new_junction(main_tpl, branch, 0, atoi(attr[2]),
						atoi(attr[4]), 0);
		//atoi(attr[3]));
		//p_tpl(main_tpl);
		//p_tpl(branch);
		if (!main_tpl->m_juncs)
			main_tpl->m_juncs = g_ptr_array_sized_new(4);
		if (!branch->b_juncs)
			branch->b_juncs = g_ptr_array_sized_new(4);
		g_ptr_array_add(main_tpl->m_juncs, jun);
		g_ptr_array_add(branch->b_juncs, jun);
		g_ptr_array_add(all_junctions, jun);
	}
	bwa_free_read_seq(n_ctgs, seqs);
}

void blat_ref(char *joint_fa, char *joint_psl) {
	GPtrArray *paired_hits = NULL;
	int i = 0, j = 0, x = 0, has_hit = 0, has_valid_hit = 0;
	uint64_t n_joint = 0;
	int locus = 0, ori = 0;
	char *attr[18];
	bwa_seq_t *joints = load_reads(joint_fa, &n_joint), *s = NULL;
	blat_hit *h = NULL;
	char *good = NULL, *junc_name = NULL, *name_copy = NULL;
	paired_hits = read_blat_hits(joint_psl);
	for (i = 0; i < n_joint; i++) {
		s = &joints[i];
		p_query(__func__, s);
		has_valid_hit = 0;
		has_hit = 0;
		for (j = 0; j < paired_hits->len; j++) {
			h = (blat_hit*) g_ptr_array_index(paired_hits, j);
			if (strcmp(s->name, h->query) == 0) {
				has_hit = 1;
				x = 0;
				name_copy = strdup(s->name);
				attr[0] = strtok(name_copy, "_");
				while (attr[x] != NULL) {
					attr[++x] = strtok(NULL, "_");
				}
				locus = atoi(attr[2]);
				ori = atoi(attr[3]);
				printf("Locus: %d; ori: %d\n", locus, ori);
				if (h->q_start < locus && h->q_end > locus) {
					if (h->n_mismatch <= 4 && h->n_ref_gap == 0
							&& h->n_query_gap == 0) {
						has_valid_hit = 1;
						good = strdup(h->ref);
						junc_name = strdup(h->query);
					}
				}
				test_print_blat_hit(h);
			}
		}
		if (!has_hit) {
			printf("REPORT No hit for %s. \n", s->name);
		} else {
			if (has_valid_hit) {
				printf("REPORT valid hit for %s on %s \n", junc_name, good);
			} else {
				printf("REPORT invalid hit for %s \n", s->name);
			}
		}
	}

}
