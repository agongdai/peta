/*
 * peeva.c
 *
 *  Created on: Jan 1, 2012
 *      Author: Cai Shaojiang
 */

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <glib.h>
#include "peeva.h"
#include "utils.h"
#include "edge.h"
#include "roadmap.h"
#include "bwase.h"
#include "peseq.h"

int total_base = 0;
int n_gene = 0;
int n_exon = 0;
int n_stand_alone = 0;
int n_stand_alone_ass = 0;
int n_full_len = 0;
int n_target_egs = 0;
int stand_alone_base = 0;
float stand_alone_per = 0;
int ass_base = 0;
int n_ass_contig = 0;
int n_not_aligned;
int n_not_base;
int max_ctg_len;
int min_ctg_len;
int n50;
int max_ctg_not_aligned_len;
int second_ctg_not_aligned_len;
int third_ctg_not_aligned_len;
int opt_n50;
int best_n50;
float accuracy;
int f_draw_graph;
GPtrArray *all_len;
GPtrArray *not_aligned_ctgs;
GPtrArray *ori_tx;
exonarray *uni_ea;
txarray *ta;
txarray *tx_not_touched;
txarray *tx_stand_alone;

int eva_usage() {
	fprintf(stderr, "\n");
	fprintf(
			stderr,
			"Usage:   peta eva [options] <blast.sam> <result.txt> <exon files> <genes> \n\n");
	return 1;
}

void chomp(char *s) {
	if (s[strlen(s) - 1] == '\n')
		s[strlen(s) - 1] = 0;
}

eva_occ *new_occ() {
	eva_occ *o = (eva_occ*) malloc(sizeof(eva_occ));
	o->end = 0;
	o->start = 0;
	o->q_id = 0;
	o->ali_len = 0;
	o->evalue = 0;
	o->percentage = 0;
	o->r_len = 0;
	o->q_len = 0;
	o->r_id = 0;
	return o;
}

exon *new_exon() {
	exon *e = (exon*) malloc(sizeof(exon));
	e->id = 0;
	e->label = (char*) calloc(32, sizeof(char));
	e->len = 0;
	e->ta = g_ptr_array_sized_new(8);
	e->merged_to = e;
	e->drawed = 0;
	return e;
}

tx *new_tx() {
	tx *t = (tx*) malloc(sizeof(tx));
	t->ea = g_ptr_array_sized_new(32);
	t->len = 0;
	t->name = 0;
	t->touched = 0;
	return t;
}

void p_occ(eva_occ *o) {
	printf("[p_occ] %s:len=%d, %s:len=%d, %f, %d, %d, %d, %f \n", o->q_id,
			o->q_len, o->r_id, o->r_len, o->percentage, o->ali_len, o->start,
			o->end, o->evalue);
}

void p_exon(exon *e) {
	show_debug_msg(__func__, "%d: %d\n", e->id, e->len);
}

void p_tx(tx *t) {
	int i = 0;
	exon *ex = 0;
	show_debug_msg(__func__, "%s: %d \n", t->name, t->len);
	//	printf("\t");
	//	for (i = 0; i < t->ea->len; i++) {
	//		ex = g_ptr_array_index(t->ea, i);
	//		printf("%d[%d], ", ex->id, ex->len);
	//	}
	//	printf("\n");
}

gint cmp_occ(gpointer a, gpointer b) {
	eva_occ *occ_a = *((eva_occ**) a);
	eva_occ *occ_b = *((eva_occ**) b);
	return (occ_a->start - occ_b->start);
}

void get_tx_not_touched(gpointer key, gpointer value, gpointer user_data) {
	char *r_id = (char*) key;
	tx *t;
	int i = 0, t_len = 0;
	for (i = 0; i < ta->len; i++) {
		t = g_ptr_array_index(ta, i);
		if (strcmp(t->name, r_id) == 0) {
			t_len = t->len;
			//total_base += t_len;
			t->touched = 1;
		}
	}
}

void occ_p_iter(gpointer key, gpointer value, gpointer user_data) {
	char *r_id = (char*) key;
	occarray *occs = (occarray*) value;
	eva_occ *o = 0;
	tx *t;
	int i = 0, t_len = 0;
	for (i = 0; i < ta->len; i++) {
		t = g_ptr_array_index(ta, i);
		if (strcmp(t->name, r_id) == 0) {
			t_len = t->len;
			break;
		}
	}
	printf("%s [%d]: ", r_id, t_len);
	g_ptr_array_sort(occs, (GCompareFunc) cmp_occ);
	for (i = 0; i < occs->len; i++) {
		o = g_ptr_array_index(occs, i);
		printf("[%d, %d] ", o->start, o->end);
	}
	printf("\n");
}

void occ_c_iter(gpointer key, gpointer value, gpointer user_data) {
	occarray *occs = (occarray*) value;
	eva_occ *o = 0, *o_pre = 0;
	int i = 0;
	float tx_len = 0;
	g_ptr_array_sort(occs, (GCompareFunc) cmp_occ);
	tx *t;
	if (occs->len == 1) {
		o = g_ptr_array_index(occs, 0);
		for (i = 0; i < ta->len; i++) {
			t = g_ptr_array_index(ta, i);
			if (strcmp(t->name, (char*) key) == 0) {
				tx_len = t->len;
				break;
			}
		}
		if (abs(o->end - o->start) > tx_len * 0.9) {
			n_full_len++;
			for (i = 0; i < tx_stand_alone->len; i++) {
				t = g_ptr_array_index(tx_stand_alone, i);
				if (strcmp(t->name, key) == 0) {
					n_stand_alone_ass++;
					break;
				}
			}
		}
	}
	for (i = 0; i < occs->len; i++) {
		o = g_ptr_array_index(occs, i);
		// @TODO: here only consider forward strand.
		if (o->end < o->start) {
			o_pre = o;
			continue;
		}
		if (i == 0) {
			ass_base += o->end - o->start + 1;
		} else {
			// To avoid duplicate counting:
			// NM_001042362: [1, 771] [1, 770] [772, 1709] [772, 1709] [1709, 2014] [1709, 2014]
			if (o->end <= o_pre->end) {
				o_pre = o;
				continue;
			}
			if (o->start < o_pre->end) {
				ass_base += o->end - o_pre->end;
			} else {
				ass_base += o->end - o->start + 1;
			}
		}
		o_pre = o;
	}
}

gint cmp_ints(gpointer a, gpointer b) {
	int* x = (int*) a;
	int* y = (int*) b;
	return *y - *x;
}

gint cmp_ctgs(gpointer a, gpointer b) {
	edge* eg_a = *(edge**) a;
	edge* eg_b = *(edge**) b;
	return eg_b->len - eg_a->len;
}

gint cmp_exons(gpointer a, gpointer b) {
	exon* ex_a = *(exon**) a;
	exon* ex_b = *(exon**) b;
	return ex_a->id - ex_b->id;
}

gint cmp_exons_by_len(gpointer a, gpointer b) {
	exon* ex_a = *(exon**) a;
	exon* ex_b = *(exon**) b;
	return ex_b->len - ex_a->len;
}

int array_find(txarray *array, void *ptr) {
	int i = 0;
	if (!array || !ptr)
		return NOT_FOUND;
	for (i = 0; i < array->len; i++) {
		if (g_ptr_array_index(array, i) == ptr) {
			return i;
		}
	}
	return NOT_FOUND;
}

void get_tx_info(char *exon_fn, char *gene_fn) {
	char buf[BUFSIZ];
	exon *ex = 0, *ex_i;
	tx *t = 0, *tx_i = 0;
	FILE *exon_fp = xopen(exon_fn, "r");
	FILE *gene_fp = xopen(gene_fn, "r");
	char *attr[BUFSIZ];
	int exon_id = 0, i = 0, j = 0, tx_len = 0, is_stand_alone_gene = 0,
			max_exon_id = 0, is_new_ex = 1, is_new_tx = 1;
	show_debug_msg(__func__, "exon_fn: %s \n", exon_fn);
	show_debug_msg(__func__, "gene_fn: %s \n", gene_fn);
	while (fgets(buf, sizeof(buf), exon_fp)) {
		i = 0;
		attr[0] = strtok(buf, ":");
		exon_id = atoi(attr[0]);
		while (attr[i] != NULL) { //ensure a pointer was found
			attr[++i] = strtok(NULL, ":"); //continue to tokenize the string
		}
		ex = new_exon();
		ex->id = exon_id;
		ex->len = atoi(attr[1]);
		sprintf(ex->label, "%d:%d", exon_id, ex->len);
		for (j = 0; j < uni_ea->len; j++) {
			ex_i = g_ptr_array_index(uni_ea, j);
			if (ex_i->id == ex->id) {
				is_new_ex = 0;
				break;
			}
		}
		if (is_new_ex)
			g_ptr_array_add(uni_ea, ex);
		is_new_ex = 1;
		if (exon_id >= n_exon) {
			n_exon = exon_id;
		}
	}
	g_ptr_array_sort(uni_ea, (GCompareFunc) cmp_exons);
	//g_ptr_array_foreach(uni_ea, (GFunc) p_exon, NULL);
	fclose(exon_fp);
	n_exon++;
	while (fgets(buf, sizeof(buf), gene_fp)) {
		i = 0;
		is_stand_alone_gene = 1;
		// printf("%s \n", buf);
		if (strstr(buf, "#") != NULL) {
			attr[0] = strtok(buf, "=");
			while (attr[i] != NULL) { //ensure a pointer was found
				attr[++i] = strtok(NULL, "="); //continue to tokenize the string
			}
			tx_len = atoi(attr[1]);
			t = new_tx();
			t->len = tx_len;
			attr[0] = strtok(buf, "|");
			//printf("%s \n", attr[0]);
			t->name = strdup(attr[0]);
			memmove(t->name, &t->name[1], strlen(t->name) - 1);
			t->name[strlen(t->name) - 1] = '\0';
			is_new_tx = 1;
			// Some non-coding/micro RNAs have the same name, just treate them as one.
			for (i = 0; i < ta->len; i++) {
				tx_i = g_ptr_array_index(ta, i);
				if (strcmp(t->name, tx_i->name) == 0) {
					is_new_tx = 0;
				}
			}
			if (is_new_tx) {
				g_ptr_array_add(ta, t);
				total_base += tx_len;
				n_gene++;
			}
		} else {
			attr[0] = strtok(buf, ",");
			while (attr[i] != NULL) { //ensure a pointer was found
				attr[++i] = strtok(NULL, ","); //continue to tokenize the string
			}

			for (j = 0; j < i; j++) {
				// printf("i, j: %d, %d \n", i, j);
				exon_id = atoi(strtok(attr[j], "["));
				ex = g_ptr_array_index(uni_ea, exon_id);
				if (ex) {
					g_ptr_array_add(t->ea, ex);
					g_ptr_array_add(ex->ta, t);
				}
				if (atoi(attr[j]) == max_exon_id + 1) {
					max_exon_id++;
				} else {
					// show_debug_msg(__func__, "%s", meta);
					is_stand_alone_gene = 0;
				}
			}
			if (is_stand_alone_gene) {
				n_stand_alone++;
				stand_alone_base += tx_len;
				if (is_new_tx) {
					g_ptr_array_add(tx_stand_alone, t);
				}
			}
		}
	}
	//g_ptr_array_foreach(ta, (GFunc) p_tx, NULL);
	stand_alone_per = (float) stand_alone_base / total_base;
	fclose(gene_fp);
}

GHashTable *pe_eva_parse(char *sam_fn, char *res_fn) {
	char buf[BUFSIZ];
	char *id_str, *attr[10], *r_attr[5];
	int i = 0, q_len = 0, contig_id = 0;
	// Key is the transcript name, value is an array of pointers, each of which is hit
	GHashTable *o_arr = g_hash_table_new(g_str_hash, g_str_equal);
	occarray *o_arr_i = g_ptr_array_sized_new(16);
	FILE *sam_fp = xopen(sam_fn, "r");
	FILE *res_fp = xopen(res_fn, "w");
	eva_occ *o = 0;
	edge *eg = 0;
	while (fgets(buf, sizeof(buf), sam_fp)) {
		i = 0;
		// printf("%s", buf);
		if (strstr(buf, "#") != NULL) {
			if (strstr(buf, "Query") != NULL) {
				attr[0] = strtok(buf, "=");
				while (attr[i] != NULL) { //ensure a pointer was found
					attr[++i] = strtok(NULL, "="); //continue to tokenize the string
				}
				q_len = atoi(attr[1]);
				i = 0;

				attr[0] = strtok(buf, " ");
				while (attr[i] != NULL) { //ensure a pointer was found
					attr[++i] = strtok(NULL, " "); //continue to tokenize the string
				}
				contig_id = atoi(attr[2]);
				eg = new_eg();
				eg->id = contig_id;
				eg->len = q_len;
				g_ptr_array_add(all_len, eg);
			}
			if (strstr(buf, " 0 hits found") != NULL) {
				n_not_aligned++;
				n_not_base += q_len;
				g_ptr_array_add(not_aligned_ctgs, eg);
			}
			continue;
		}
		o = new_occ();
		o->q_len = q_len;
		chomp(buf);
		// printf("---------------------------------------\n %s \n", buf);
		attr[0] = strtok(buf, "\t");
		while (attr[i] != NULL) { //ensure a pointer was found
			attr[++i] = strtok(NULL, "\t"); //continue to tokenize the string
		}
		o->q_id = strdup(attr[0]);
		id_str = attr[1];
		i = 0;
		r_attr[0] = strtok(id_str, "|");
		while (r_attr[i] != NULL) { //ensure a pointer was found
			r_attr[++i] = strtok(NULL, "|"); //continue to tokenize the string
		}
		o->r_id = strdup(r_attr[1]);
		i = 0;
		r_attr[0] = strtok(r_attr[2], "=");
		while (r_attr[i] != NULL) { //ensure a pointer was found
			r_attr[++i] = strtok(NULL, "="); //continue to tokenize the string
		}
		o->r_len = atoi(r_attr[1]);
		o->percentage = atof(attr[2]) / 100;
		o->ali_len = atoi(attr[3]);
		o->start = atoi(attr[6]);
		o->end = atoi(attr[7]);
		o->evalue = atoi(attr[8]);
		// p_occ(o);
		if (o->ali_len <= 35 || o->end < o->start)
			continue;

		o_arr_i = g_hash_table_lookup(o_arr, o->r_id);
		if (o_arr_i == NULL) {
			o_arr_i = g_ptr_array_sized_new(16);
			g_ptr_array_add(o_arr_i, o);
			g_hash_table_insert(o_arr, o->r_id, o_arr_i);
		} else {
			g_ptr_array_add(o_arr_i, o);
		}
		//		if (counter++ > 10)
		//			break;
	}
	fclose(sam_fp);
	fclose(res_fp);
	return o_arr;
}

void p_ctgs(gpointer a, gpointer user_data) {
	edge *eg = (edge*) a;
	printf("[%d, %d]\n", eg->id, eg->len);
}

void p_summary(char *sam) {
	int i = 0;
	tx *t;
	fprintf(stderr,
			"\n------------------------------------------------------------------- \n");
	fprintf(stderr, "[p_summary] BLAST SAM file: \t\t %s \n", sam);
	fprintf(stderr, "[p_summary] Total base: \t\t %d \n", total_base);
	fprintf(stderr, "[p_summary] # of genes: \t\t %d \n", n_gene);
	fprintf(stderr, "[p_summary] # of exons: \t\t %d \n", n_exon);
	fprintf(stderr, "[p_summary] Optimal N50: \t\t %d \n", opt_n50);
	fprintf(stderr, "[p_summary] Best N50: \t\t\t %d \n", best_n50);
	fprintf(stderr, "[p_summary] # of stand alone gene: \t %d \n",
			n_stand_alone);
	fprintf(stderr, "[p_summary] # of stand alone base: \t %d \n",
			stand_alone_base);
	fprintf(stderr, "[p_summary] %c of stand alone base: \t %f \n", '%',
			stand_alone_per);
	fprintf(stderr, "[p_summary] # of target edges: \t %d \n", n_target_egs);
	fprintf(stderr, "[p_summary]  \n");

	fprintf(stderr, "[p_summary] Assembled base: \t\t %d \n", ass_base);
	fprintf(stderr, "[p_summary] # of contigs: \t\t %d \n", n_ass_contig);
	fprintf(stderr, "[p_summary] # of Full length: \t\t %d \n", n_full_len);
	fprintf(stderr, "[p_summary] Stand-alone assembled: \t\t %d \n",
			n_stand_alone_ass);
	fprintf(stderr, "[p_summary] Longest contig: \t\t %d \n", max_ctg_len);
	fprintf(stderr, "[p_summary] # of Ctgs not aligend: \t %d \n",
			n_not_aligned);
	fprintf(stderr, "[p_summary] Longest ctgs not aligend: \t %d, %d, %d \n",
			max_ctg_not_aligned_len, second_ctg_not_aligned_len,
			third_ctg_not_aligned_len);
	fprintf(stderr, "[p_summary] Bases not aligend: \t\t %d \n", n_not_base);
	fprintf(stderr, "[p_summary] N50 value:  \t\t %d \n", n50);
	fprintf(stderr, "[p_summary] Base coverage: \t\t %f \n", accuracy);
	printf("[p_summary] Transcripts not touched [tx id, tx length]: \n");
	for (i = 0; i < ta->len; i++) {
		t = g_ptr_array_index(ta, i);
		if (!t->touched)
			printf("[%s, %d]\n", t->name, t->len);
	}
	printf(
			"\n[p_summary] %d Contigs not aligned [contig id, contig length]: \n",
			not_aligned_ctgs->len);
	g_ptr_array_foreach(not_aligned_ctgs, (GFunc) p_ctgs, NULL);
	fprintf(stderr,
			"\n------------------------------------------------------------------- \n");
}

void set_ass_ctgs(char *ctgs_fn) {
	FILE *ctgs_fp = xopen(ctgs_fn, "r");
	char buf[BUFSIZ];
	char *attr[4];
	int i = 0;
	while (fgets(buf, sizeof(buf), ctgs_fp)) {
		chomp(buf);
		i = 0;
		if (strstr(buf, ">") != NULL) {
			attr[0] = strtok(buf, "=");
			while (attr[i] != NULL) { //ensure a pointer was found
				attr[++i] = strtok(NULL, "="); //continue to tokenize the string
			}
			n_ass_contig++;
			// ass_base += atoi(attr[1]);
			max_ctg_len = (max_ctg_len > atoi(attr[1])) ? max_ctg_len : atoi(
					attr[1]);
		}
	}
}

gint cmp_seqs(gpointer a, gpointer b) {
	bwa_seq_t* eg_a = *(bwa_seq_t**) a;
	bwa_seq_t* eg_b = *(bwa_seq_t**) b;
	return eg_b->len - eg_a->len;
}

void get_opt_info(char *tx_fn) {
	bwa_seqio_t *ks;
	int n_seqs = 0, i = 0, t_len = 0, len = 0;
	bwa_seq_t *seqs, *s;
	ks = bwa_open_reads(BWA_MODE, tx_fn);
	while ((seqs = bwa_read_seq(ks, 0xa00000, &n_seqs, BWA_MODE, 0)) != 0) {
		for (i = 0; i < n_seqs; i++) {
			s = &seqs[i];
			g_ptr_array_add(ori_tx, s);
		}
		break;
	}
	pe_reverse_seqs(seqs, n_seqs);
	g_ptr_array_sort(ori_tx, (GCompareFunc) cmp_seqs);
	for (i = 0; i < ori_tx->len; i++) {
		s = g_ptr_array_index(ori_tx, i);
		t_len += s->len;
	}
	for (i = 0; i < ori_tx->len; i++) {
		s = g_ptr_array_index(ori_tx, i);
		len += s->len;
		if (len > t_len * 0.5) {
			opt_n50 = s->len;
			break;
		}
	}
}

int exon_adj(exon *e1, exon *e2) {
	txarray *ta_1 = e1->ta;
	txarray *ta_2 = e2->ta;
	int i = 0, j = 0;
	tx *t;
	exon *e, *e_after;
	if (ta_1->len != ta_2->len)
		return 0;
	for (i = 0; i < ta_1->len; i++) {
		t = g_ptr_array_index(ta_1, i);
		for (j = 0; j < t->ea->len - 1; j++) {
			e = g_ptr_array_index(t->ea, j);
			e_after = g_ptr_array_index(t->ea, j + 1);
			if ((e == e1 && e_after != e2) || (e != e1 && e_after == e2)) {
				return 0;
			}
		}
	}
	return 1;
}

void draw_graph(char *s_read_fn, char *tx_fn) {
	FILE *s_read_fp = xopen(s_read_fn, "r");
	FILE *graph_fp = xopen("graph/target.dot", "w");
	FILE *tx_egs = xopen("read/edges.txt", "w");
	txarray *target_ta = g_ptr_array_sized_new(N_CTGS), *exon_ta = 0;
	GPtrArray *conns = g_ptr_array_sized_new(N_CTGS);
	GPtrArray *cpns = g_ptr_array_sized_new(N_CTGS);
	exon *e, *e_after;
	tx *t, *t_exon;
	bwa_seq_t *tx_seq = 0, *eg_seq = 0;
	char buf[BUFSIZ], shape_str[BUFSIZ], *shape_str_i;
	char *attr[BUFSIZ];
	int i = 0, j = 0, k = 0, merged = 0, sum_exon_len = 0, acc_exon_len = 0,
			acc_len = 0;
	// List all transcripts that should be assembled.
	while (fgets(buf, sizeof(buf), s_read_fp)) {
		chomp(buf);
		for (i = 0; i < ta->len; i++) {
			t = g_ptr_array_index(ta, i);
			if (strcmp(t->name, buf) == 0) {
				g_ptr_array_add(target_ta, t);
			}
		}
	}
	//g_ptr_array_foreach(target_ta, (GFunc) p_tx, NULL);
	// Get the "closure" of target transcripts
	while (1) {
		exon_ta = g_ptr_array_sized_new(256);
		// For each existing transcripts
		for (i = 0; i < target_ta->len; i++) {
			t = g_ptr_array_index(target_ta, i);
			// For each exons of the transcript
			for (j = 0; j < t->ea->len; j++) {
				e = g_ptr_array_index(t->ea, j);
				// For each transcript which contains current exon.
				for (k = 0; k < e->ta->len; k++) {
					t_exon = g_ptr_array_index(e->ta, k);
					// If not added, add it!
					if (array_find(target_ta, t_exon) == NOT_FOUND
							&& array_find(exon_ta, t_exon) == NOT_FOUND) {
						g_ptr_array_add(exon_ta, t_exon);
					}
				}
			}
		}
		if (exon_ta->len > 0)
			g_ptr_array_concat(target_ta, exon_ta);
		else
			break;
	}
	show_debug_msg(__func__, "-----------------------------------------\n");
	//g_ptr_array_foreach(target_ta, (GFunc) p_tx, NULL);

	fputs("digraph g { \n\trankdir = LR \n", graph_fp);
	for (i = 0; i < target_ta->len; i++) {
		t = g_ptr_array_index(target_ta, i);
		for (j = 0; j < t->ea->len - 1; j++) {
			e = g_ptr_array_index(t->ea, j);
			e_after = g_ptr_array_index(t->ea, j + 1);
			merged = exon_adj(e, e_after);
			if (merged) {
				// If the second exon is not merged yet, accumulate the length;
				// By default, e->merged_to = e;
				// Different transcripts may have the same exon, here avoids duplicate counting.
				if (e_after->merged_to == e_after) {
					e->merged_to->len += e_after->len;
					sprintf(e->merged_to->label, "%d-%d: %d", e->merged_to->id,
							e_after->id, e->merged_to->len);
					e_after->merged_to = e->merged_to;
				}
			} else {
				sprintf(shape_str, "\t%d -> %d \n", e->merged_to->id,
						e_after->id);
				for (k = 0; k < conns->len; k++) {
					shape_str_i = g_ptr_array_index(conns, k);
					if (strcmp(shape_str_i, shape_str) == 0)
						break;
				}
				if (k == conns->len) {
					fputs(shape_str, graph_fp);
					g_ptr_array_add(conns, strdup(shape_str));
				}
			}
		}
	}

	for (i = 0; i < target_ta->len; i++) {
		t = g_ptr_array_index(target_ta, i);
		for (j = 0; j < ori_tx->len; j++) {
			tx_seq = g_ptr_array_index(ori_tx, j);
			k = 0;
			attr[0] = strtok(tx_seq->name, "|");
			while (attr[k] != NULL) { //ensure a pointer was found
				attr[++k] = strtok(NULL, "|"); //continue to tokenize the string
			}
			if (attr[1] && strcmp(attr[1], t->name) == 0) {
				break;
			}
		}
		acc_len = 0;
		for (j = 0; j < t->ea->len; j++) {
			e = g_ptr_array_index(t->ea, j);
			// Only print the exons who are merged to.
			if (!e->merged_to->drawed) {
				sprintf(shape_str, "\t%d [shape=box, label=\"[%s] %s\"] \n",
						e->merged_to->id, t->name, e->merged_to->label);
				n_target_egs++;
				fputs(shape_str, graph_fp);
				e->merged_to->drawed = 1;
				g_ptr_array_add(cpns, e->merged_to);
				sum_exon_len += e->merged_to->len;
				eg_seq = new_seq(tx_seq, e->merged_to->len, acc_len);
				acc_len += e->merged_to->len;
				sprintf(shape_str, ">%s_%s\n", t->name, e->merged_to->label);
				save_con(shape_str, eg_seq, tx_egs);
			}
		}
	}
	fputs("}\n", graph_fp);

	g_ptr_array_sort(cpns, (GCompareFunc) cmp_exons_by_len);
	for (i = 0; i < cpns->len; i++) {
		e = g_ptr_array_index(cpns, i);
		acc_exon_len += e->len;
		show_debug_msg(__func__, "%s: %d => %d / %d\n", e->label, e->len,
				acc_exon_len, sum_exon_len);
		if (acc_exon_len * 2 >= sum_exon_len) {
			best_n50 = e->len;
			break;
		}
	}

	fclose(s_read_fp);
	fclose(graph_fp);
	fclose(tx_egs);
}

void pe_eva_core(char *sam_fn, char *res_fn, char *exon_fn, char *gene_fn,
		char *ctgs_fn, char *tx_fn, char *s_read_fn) {
	int i = 0, len = 0, t_len = 0;
	edge *eg = 0;
	all_len = g_ptr_array_sized_new(N_CTGS);
	ori_tx = g_ptr_array_sized_new(N_CTGS);
	uni_ea = g_ptr_array_sized_new(N_CTGS);
	ta = g_ptr_array_sized_new(N_CTGS);
	tx_not_touched = g_ptr_array_sized_new(N_CTGS);
	not_aligned_ctgs = g_ptr_array_sized_new(N_NOT_ALIGNED_CTGS);
	tx_stand_alone = g_ptr_array_sized_new(N_NOT_ALIGNED_CTGS);
	GHashTable *ht = pe_eva_parse(sam_fn, res_fn);

	get_opt_info(tx_fn);
	get_tx_info(exon_fn, gene_fn);
	if (f_draw_graph)
		draw_graph(s_read_fn, tx_fn);
	set_ass_ctgs(ctgs_fn);

	g_hash_table_foreach(ht, (GHFunc) occ_c_iter, NULL);
	show_msg(__func__,
			"-------------------------------------------------------------------\n");
	//total_base = 0; // tmp!!!
	g_hash_table_foreach(ht, (GHFunc) get_tx_not_touched, NULL);
	g_hash_table_foreach(ht, (GHFunc) occ_p_iter, NULL);
	show_msg(__func__,
			"-------------------------------------------------------------------\n");
	accuracy = (float) ass_base / total_base;
	g_ptr_array_sort(all_len, (GCompareFunc) cmp_ctgs);
	g_ptr_array_sort(not_aligned_ctgs, (GCompareFunc) cmp_ctgs);
	if (not_aligned_ctgs->len > 3) {
		eg = g_ptr_array_index(not_aligned_ctgs, 0);
		max_ctg_not_aligned_len = eg->len;
		eg = g_ptr_array_index(not_aligned_ctgs, 1);
		second_ctg_not_aligned_len = eg->len;
		eg = g_ptr_array_index(not_aligned_ctgs, 2);
		third_ctg_not_aligned_len = eg->len;
	}
	show_msg(__func__, "Assembled bases: %d \n", ass_base);
	for (i = 0; i < all_len->len; i++) {
		eg = g_ptr_array_index(all_len, i);
		t_len += eg->len;
	}
	for (i = 0; i < all_len->len; i++) {
		eg = g_ptr_array_index(all_len, i);
		if (i == 0) {
			max_ctg_len = eg->len;
		}
		len += eg->len;
		// printf("[%d, %d], \t", eg->id, eg->len);
		if (len > t_len * 0.5) {
			n50 = eg->len;
			break;
		}
	}
	// show_debug_msg(__func__, "\n");
	p_summary(sam_fn);
}

int pe_eva(int argc, char *argv[]) {
	clock_t t = clock();
	// int c;
	//	while ((c = getopt(argc, argv, "g:")) >= 0) {
	//		switch (c) {
	//		case 'g':
	//			f_draw_graph = atoi(optarg);
	//			break;
	//		}
	//	}
	f_draw_graph = 1;
	if (optind + 7 > argc) {
		return eva_usage();
	}
	pe_eva_core(argv[optind], argv[optind + 1], argv[optind + 2], argv[optind
			+ 3], argv[optind + 4], argv[optind + 5], argv[optind + 6]);
	fprintf(stderr, "[pe_ass] Done: %.2f sec\n", (float) (clock() - t)
			/ CLOCKS_PER_SEC);
	return 0;
}
