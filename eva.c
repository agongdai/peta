/*
 * eva.c
 *
 *  Created on: 16-Apr-2012
 *      Author: carl
 */

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <glib.h>
#include "utils.h"
#include "edge.h"
#include "roadmap.h"
#include "bwase.h"
#include "peseq.h"
#include "eva.h"
int compute_ori = 0;
char *sam_fn;
char *exons_fn;
char *edges_fn;
char *tx_fn;
char *contigs_fn;
char *start_reads_fn;
char *result_fn;
char *tx_sum_fn = "read/tx_sum.txt";
tx_info *ori_info;

int tmp_total_base = 0;
int tmp_got_base = 0;
int tmp_n_cufflinks = 0;

tx_info *init_info() {
	tx_info *i = (tx_info*) malloc(sizeof(tx_info));
	i->base_sd_total = 0.0;
	i->best_n50 = 0;
	i->graph_fn = "graph/target.dot";
	i->n_base = 0;
	i->n_uni_base = 0;
	i->n_exon = 0;
	i->n_graph_egs = 0;
	i->n_sd_base = 0;
	i->n_sd_tx = 0;
	i->n_tx = 0;
	i->opt_n50 = 0;
	i->sum_fn = tx_sum_fn;
	i->cpn_fn = "read/edges.txt";
	i->slots = (int*) calloc(MAX_LEN_FOR_PLOT / SLOT_SIZE + 1, sizeof(int));
	i->n_slot = MAX_LEN_FOR_PLOT / SLOT_SIZE + 1;
	i->txs = g_ptr_array_sized_new(N_CTGS);
	i->sd_txs = g_ptr_array_sized_new(N_CTGS);
	i->exons = g_ptr_array_sized_new(N_CTGS);
	i->tx_seqs = g_ptr_array_sized_new(N_CTGS);
	i->exon_seqs = g_ptr_array_sized_new(N_CTGS);
	return i;
}

rs_info *init_rs_info() {
	rs_info *i = (rs_info*) malloc(sizeof(rs_info));
	i->base_coverage = 0.0;
	i->ea_all = g_ptr_array_sized_new(N_CTGS);
	i->ea_full_len = g_ptr_array_sized_new(N_CTGS);
	i->ea_not_aligned = g_ptr_array_sized_new(N_CTGS);
	i->ea_not_touched = g_ptr_array_sized_new(N_CTGS);
	i->ea_sd = g_ptr_array_sized_new(N_CTGS);
	i->contigs = g_ptr_array_sized_new(N_CTGS);
	i->occ_all = g_ptr_array_sized_new(N_CTGS);
	i->max_len = 0;
	i->n50 = 0;
	i->n50_all = 0;
	i->n_base_shared = 0;
	i->n_base = 0;
	i->n_base_not_aligned = 0;
	i->n_ctgs = 0;
	i->n_full_len = 0;
	i->n_not_aligned = 0;
	i->n_sd = 0;
	i->slots = (int*) calloc(MAX_LEN_FOR_PLOT / SLOT_SIZE + 1, sizeof(int));
	i->n_slot = MAX_LEN_FOR_PLOT / SLOT_SIZE + 1;
	i->hits = g_hash_table_new(g_str_hash, g_str_equal);
	return i;
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

gint cmp_occ(gpointer a, gpointer b) {
	eva_occ *occ_a = *((eva_occ**) a);
	eva_occ *occ_b = *((eva_occ**) b);
	return (occ_a->start - occ_b->start);
}

void occ_p_iter(gpointer key, gpointer value, gpointer user_data) {
	char *r_id = (char*) key;
	occarray *occs = (occarray*) value;
	eva_occ *o = 0, *o_pre = 0;
	bwa_seq_t *t;
	int i = 0, t_len = 0;
	int n_base_got = 0;
	for (i = 0; i < ori_info->tx_seqs->len; i++) {
		t = g_ptr_array_index(ori_info->tx_seqs, i);
		if (strcmp(t->name, r_id) == 0) {
			t_len = t->len;
			break;
		}
	}
	tmp_n_cufflinks++;
	printf("%s [%d]: ", r_id, t_len);
	tmp_total_base += t_len;
	g_ptr_array_sort(occs, (GCompareFunc) cmp_occ);
	for (i = 0; i < occs->len; i++) {
		o = g_ptr_array_index(occs, i);
		printf("(%s, %d): [%d, %d], ", o->q_id, o->q_len, o->start, o->end);
		if (i == 0) {
			n_base_got += o->end - o->start;
		} else {
			if (o->start < o_pre->end)
				n_base_got += o->end - o_pre->end;
			else
				n_base_got += o->end - o->start;
		}
		o_pre = o;
	}
	tmp_got_base += n_base_got;
	printf("\t %d/%d:%f", n_base_got, t_len, ((float) (n_base_got)
			/ (float) t_len));
	printf("\t %d", tmp_n_cufflinks);
	printf("\n");
}

void occ_c_iter(gpointer key, gpointer value, gpointer user_data) {
	rs_info *info = (rs_info*) user_data;
	occarray *occs = (occarray*) value;
	eva_occ *o = 0, *o_pre = 0;
	int i = 0;
	float tx_len = 0;
	g_ptr_array_sort(occs, (GCompareFunc) cmp_occ);
	tx *t;
	bwa_seq_t *seq;
	if (occs->len == 1) {
		o = g_ptr_array_index(occs, 0);
		for (i = 0; i < ori_info->tx_seqs->len; i++) {
			seq = g_ptr_array_index(ori_info->tx_seqs, i);
			//show_debug_msg(__func__, "seq name: %s \n", seq->name);
			//show_debug_msg(__func__, "align name: %s \n", (char*) key);
			if (strcmp(seq->name, (char*) key) == 0) {
				tx_len = seq->len;
				break;
			}
		}
		if (abs(o->end - o->start) > tx_len * 0.9) {
			info->n_full_len++;
			for (i = 0; i < ori_info->n_sd_tx; i++) {
				t = g_ptr_array_index(ori_info->sd_txs, i);
				if (strcmp(t->name, key) == 0) {
					info->n_sd++;
					g_ptr_array_add(info->ea_sd, t);
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
			info->n_base_shared += o->end - o->start + 1;
		} else {
			// To avoid duplicate counting:
			// NM_001042362: [1, 771] [1, 770] [772, 1709] [772, 1709] [1709, 2014] [1709, 2014]
			if (o->end <= o_pre->end) {
				o_pre = o;
				continue;
			}
			if (o->start < o_pre->end) {
				info->n_base_shared += o->end - o_pre->end;
			} else {
				info->n_base_shared += o->end - o->start + 1;
			}
		}
		o_pre = o;
	}
}

char *fix_len(char *str, const int len) {
	int i = 0;
	char *new_str = (char*) malloc(sizeof(char) * len);
	if (strlen(str) > len)
		return str;
	memcpy(new_str, str, strlen(str));
	for (i = strlen(str); i < len; i++) {
		new_str[i] = ' ';
	}
	new_str[i] = '\0';
	return new_str;
}

void write_tx_sum(tx_info *info) {
	char item[BUFSIZ], *str;
	FILE *fp = xopen(info->sum_fn, "w");
	int i = 0;

	fputs(
			"-------------------------------------------------------------------\n",
			fp);
	str = fix_len("Original transcript file: ", ATTR_STR_LEN);
	sprintf(item, "%s\t%s\n", str, tx_fn);
	fputs(item, fp);

	str = fix_len("Exons file: ", ATTR_STR_LEN);
	sprintf(item, "%s\t%s\n", str, edges_fn);
	fputs(item, fp);

	str = fix_len("Total base: ", ATTR_STR_LEN);
	sprintf(item, "%s\t%d\n", str, info->n_base);
	fputs(item, fp);

	str = fix_len("Unique base: ", ATTR_STR_LEN);
	sprintf(item, "%s\t%d\n", str, info->n_uni_base);
	fputs(item, fp);

	str = fix_len("# of transcripts: ", ATTR_STR_LEN);
	sprintf(item, "%s\t%d\n", str, info->n_tx);
	fputs(item, fp);

	str = fix_len("# of exons: ", ATTR_STR_LEN);
	sprintf(item, "%s\t%d\n", str, info->n_exon);
	fputs(item, fp);

	str = fix_len("# of stand alone gene: ", ATTR_STR_LEN);
	sprintf(item, "%s\t%d\n", str, info->n_sd_tx);
	fputs(item, fp);

	str = fix_len("# of stand alone base: ", ATTR_STR_LEN);
	sprintf(item, "%s\t%d\n", str, info->n_sd_base);
	fputs(item, fp);

	str = fix_len("% of stand alone base: ", ATTR_STR_LEN);
	sprintf(item, "%s\t%f\n", str, info->base_sd_total);
	fputs(item, fp);

	str = fix_len("Optimal N50: ", ATTR_STR_LEN);
	sprintf(item, "%s\t%d\n", str, info->opt_n50);
	fputs(item, fp);

	fputs("\n", fp);

	str = fix_len("# of target edges: ", ATTR_STR_LEN);
	sprintf(item, "%s\t%d\n", str, info->n_graph_egs);
	fputs(item, fp);

	str = fix_len("Best N50: ", ATTR_STR_LEN);
	sprintf(item, "%s\t%d\n", str, info->best_n50);
	fputs(item, fp);

	fputs("Base range accumulation:\n", fp);
	for (i = 0; i < info->n_slot; i++) {
		sprintf(item, "%d\t%d\n", i * SLOT_SIZE, info->slots[i]);
		fputs(item, fp);
	}

	fputs(
			"-------------------------------------------------------------------\n\n\n",
			fp);
	fclose(fp);
}

gint cmp_ctgs(gpointer a, gpointer b) {
	edge* eg_a = *(edge**) a;
	edge* eg_b = *(edge**) b;
	return eg_b->len - eg_a->len;
}

void write_ass_rs(rs_info *info) {
	char buf[BUFSIZ], item[BUFSIZ], *str;
	FILE *ori_fp = xopen(tx_sum_fn, "r");
	FILE *rs_fp = xopen(result_fn, "w");
	int i = 0;
	edge *t = 0;

	while (fgets(buf, sizeof(buf), ori_fp)) {
		fputs(buf, rs_fp);
	}
	fputs("Assembled results by PETA: \n", rs_fp);
	fputs(
			"-------------------------------------------------------------------\n",
			rs_fp);

	str = fix_len("Assembled base: ", ATTR_STR_LEN);
	sprintf(item, "%s\t%d\n", str, info->n_base);
	fputs(item, rs_fp);

	str = fix_len("Assembled unique base: ", ATTR_STR_LEN);
	sprintf(item, "%s\t%d\n", str, info->n_base_shared);
	fputs(item, rs_fp);

	str = fix_len("# of contigs: ", ATTR_STR_LEN);
	sprintf(item, "%s\t%d\n", str, info->n_ctgs);
	fputs(item, rs_fp);

	str = fix_len("# of Full length: ", ATTR_STR_LEN);
	sprintf(item, "%s\t%d\n", str, info->n_full_len);
	fputs(item, rs_fp);

	str = fix_len("Stand-alone assembled: ", ATTR_STR_LEN);
	sprintf(item, "%s\t%d\n", str, info->n_sd);
	fputs(item, rs_fp);

	str = fix_len("Average length: ", ATTR_STR_LEN);
	sprintf(item, "%s\t%d\n", str, info->ave_len);
	fputs(item, rs_fp);

	str = fix_len("Longest contig: ", ATTR_STR_LEN);
	sprintf(item, "%s\t%d\n", str, info->max_len);
	fputs(item, rs_fp);

	str = fix_len("# of Ctgs not aligned: ", ATTR_STR_LEN);
	sprintf(item, "%s\t%d\n", str, info->n_not_aligned);
	fputs(item, rs_fp);

	str = fix_len("Bases not aligned: ", ATTR_STR_LEN);
	sprintf(item, "%s\t%d\n", str, info->n_base_not_aligned);
	fputs(item, rs_fp);

	str = fix_len("N50 value: ", ATTR_STR_LEN);
	sprintf(item, "%s\t%d\n", str, info->n50);
	fputs(item, rs_fp);

	str = fix_len("N50 value (Aligned based on all): ", ATTR_STR_LEN);
	sprintf(item, "%s\t%d\n", str, info->n50_aligned);
	fputs(item, rs_fp);

	str = fix_len("N50 value (based on all): ", ATTR_STR_LEN);
	sprintf(item, "%s\t%d\n", str, info->n50_all);
	fputs(item, rs_fp);

	str = fix_len("Base coverage: ", ATTR_STR_LEN);
	sprintf(item, "%s\t%f\n", str, info->base_coverage);
	fputs(item, rs_fp);

	str = fix_len("Contigs not aligned [id, length]: \n\t", ATTR_STR_LEN);
	fputs(str, rs_fp);
	g_ptr_array_sort(info->ea_not_aligned, (GCompareFunc) cmp_ctgs);
	for (i = 0; i < info->ea_not_aligned->len; i++) {
		t = g_ptr_array_index(info->ea_not_aligned, i);
		sprintf(item, "[%s:%d], ", t->name, t->len);
		fputs(item, rs_fp);
	}

	fputs("\nBase range accumulation:\n", rs_fp);
	for (i = 0; i < info->n_slot; i++) {
		sprintf(item, "%d\t%d\n", i * SLOT_SIZE, info->slots[i]);
		fputs(item, rs_fp);
	}

	fputs(
			"\n-------------------------------------------------------------------\n",
			rs_fp);

	fclose(ori_fp);
	fclose(rs_fp);
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

void p_exon(exon *e) {
	show_debug_msg(__func__, "%d: %d\n", e->id, e->len);
}

void p_tx(tx *t) {
	int i = 0;
	exon *ex = 0;
	show_debug_msg(__func__, "%s: %d \n", t->name, t->len);
	printf("\t");
	for (i = 0; i < t->ea->len; i++) {
		ex = g_ptr_array_index(t->ea, i);
		printf("%d[%d], ", ex->id, ex->len);
	}
	printf("\n");
}

void p_occ(eva_occ *o) {
	if (!o)
		return;
	show_debug_msg(__func__, "start: %d \n", o->start);
	show_debug_msg(__func__, "end: %d \n", o->end);
	show_debug_msg(__func__, "q_id: %s \n", o->q_id);
	show_debug_msg(__func__, "r_id: %s \n", o->r_id);
	show_debug_msg(__func__, "percentage: %f \n", o->percentage);
	show_debug_msg(__func__, "evalue: %l \n", o->evalue);
	show_debug_msg(__func__, "ali_len: %d \n", o->ali_len);
	show_debug_msg(__func__, "q_len: %d \n", o->q_len);
	show_debug_msg(__func__, "r_len: %d \n", o->r_len);
	show_debug_msg(__func__, "--------------------------------------------- \n");
}

void chomp(char *s) {
	if (s[strlen(s) - 1] == '\n')
		s[strlen(s) - 1] = 0;
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

tx *new_tx() {
	tx *t = (tx*) malloc(sizeof(tx));
	t->ea = g_ptr_array_sized_new(32);
	t->len = 0;
	t->name = 0;
	t->touched = 0;
	return t;
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

gint cmp_txs(gpointer a, gpointer b) {
	tx* eg_a = *(tx**) a;
	tx* eg_b = *(tx**) b;
	return eg_b->len - eg_a->len;
}

gint cmp_seqs_by_len(gpointer a, gpointer b) {
	bwa_seq_t* eg_a = *(bwa_seq_t**) a;
	bwa_seq_t* eg_b = *(bwa_seq_t**) b;
	return eg_b->len - eg_a->len;
}

gint cmp_occ_by_len(gpointer a, gpointer b) {
	eva_occ* eg_a = *(eva_occ**) a;
	eva_occ* eg_b = *(eva_occ**) b;
	return eg_b->ali_len - eg_a->ali_len;
}

int cal_n50(txarray *ta) {
	int i = 0, t_len = 0, acc_len = 0;
	tx *t = 0;
	g_ptr_array_sort(ta, (GCompareFunc) cmp_txs);
	for (i = 0; i < ta->len; i++) {
		t = g_ptr_array_index(ta, i);
		t_len += t->len;
	}
	for (i = 0; i < ta->len; i++) {
		t = g_ptr_array_index(ta, i);
		acc_len += t->len;
		if (acc_len * 2 >= t_len)
			return t->len;
	}
	return INVALID;
}

int cal_seq_n50(GPtrArray *seqs) {
	int i = 0, t_len = 0, acc_len = 0;
	bwa_seq_t *t = 0;
	g_ptr_array_sort(seqs, (GCompareFunc) cmp_seqs_by_len);
	for (i = 0; i < seqs->len; i++) {
		t = g_ptr_array_index(seqs, i);
		t_len += t->len;
	}
	for (i = 0; i < seqs->len; i++) {
		t = g_ptr_array_index(seqs, i);
		acc_len += t->len;
		if (acc_len * 2 >= t_len)
			return t->len;
	}
	return INVALID;
}

void read_tx_seqs(tx_info *info) {
	bwa_seqio_t *ks;
	int i = 0;
	bwa_seq_t *seqs, *s;
	ks = bwa_open_reads(BWA_MODE, tx_fn);
	while ((seqs = bwa_read_seq(ks, 0xa00000, &(info->n_tx), BWA_MODE, 0)) != 0) {
		for (i = 0; i < info->n_tx; i++) {
			s = &seqs[i];
			info->n_base += s->len;
			g_ptr_array_add(info->tx_seqs, s);
		}
		break;
	}
}

void read_exon_info(tx_info *info) {
	char buf[BUFSIZ];
	exon *ex = 0, *ex_i;
	FILE *exon_fp = xopen(exons_fn, "r");
	char *attr[BUFSIZ];
	int exon_id = 0, i = 0, j = 0, is_new_ex = 1;
	show_debug_msg(__func__, "exon_fn: %s \n", exons_fn);
	show_debug_msg(__func__, "gene_fn: %s \n", edges_fn);
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
		for (j = 0; j < info->exons->len; j++) {
			ex_i = g_ptr_array_index(info->exons, j);
			if (ex_i->id == ex->id) {
				is_new_ex = 0;
				break;
			}
		}
		if (is_new_ex)
			g_ptr_array_add(info->exons, ex);
		is_new_ex = 1;
		if (exon_id >= info->n_exon) {
			info->n_exon = exon_id;
		}
	}
	// Sort exons by id, such that it could be acting as kind of 'hash table'
	// index 0 points to the exon whose id is 0
	g_ptr_array_sort(info->exons, (GCompareFunc) cmp_exons);
	//g_ptr_array_foreach(info->exons, (GFunc) p_exon, NULL);
	fclose(exon_fp);
	info->n_exon++;
}

void read_tx_info(tx_info *info) {
	char buf[BUFSIZ];
	char *attr[BUFSIZ];
	FILE *gene_fp = xopen(edges_fn, "r");
	tx *t = 0, *tx_i;
	exon *ex;
	int i = 0, j = 0, is_sd_tx = 1, tx_len = 0, is_new_tx = 1, exon_id = 0,
			max_exon_id = 0;
	while (fgets(buf, sizeof(buf), gene_fp)) {
		i = 0;
		is_sd_tx = 1;
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
			// Some non-coding/micro RNAs have the same name, just treat them as one.
			for (i = 0; i < info->txs->len; i++) {
				tx_i = g_ptr_array_index(info->txs, i);
				if (strcmp(t->name, tx_i->name) == 0) {
					is_new_tx = 0;
				}
			}
			if (is_new_tx) {
				g_ptr_array_add(info->txs, t);
			}
		} else {
			attr[0] = strtok(buf, ",");
			while (attr[i] != NULL) { //ensure a pointer was found
				attr[++i] = strtok(NULL, ","); //continue to tokenize the string
			}

			for (j = 0; j < i; j++) {
				// printf("i, j: %d, %d \n", i, j);
				exon_id = atoi(strtok(attr[j], "["));
				ex = g_ptr_array_index(info->exons, exon_id);
				if (ex) {
					g_ptr_array_add(t->ea, ex);
					g_ptr_array_add(ex->ta, t);
				}
				if (atoi(attr[j]) == max_exon_id + 1) {
					max_exon_id++;
				} else {
					// show_debug_msg(__func__, "%s", meta);
					is_sd_tx = 0;
				}
			}
			if (is_sd_tx && is_new_tx) {
				info->n_sd_tx++;
				info->n_sd_base += tx_len;
				g_ptr_array_add(info->sd_txs, t);
			}
		}
	}
	//g_ptr_array_foreach(ta, (GFunc) p_tx, NULL);
	info->base_sd_total = (float) info->n_sd_base / info->n_base;
	fclose(gene_fp);
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

void draw_graph(tx_info *info) {
	FILE *s_read_fp = xopen(start_reads_fn, "r");
	FILE *graph_fp = xopen(info->graph_fn, "w");
	FILE *tx_egs = xopen(info->cpn_fn, "w");
	txarray *target_ta, *exon_ta = 0;
	GPtrArray *conns = g_ptr_array_sized_new(N_CTGS);
	GPtrArray *cpns = g_ptr_array_sized_new(N_CTGS);
	exon *e, *e_after;
	tx *t, *t_exon;
	bwa_seq_t *tx_seq = 0, *eg_seq = 0;
	char buf[BUFSIZ], shape_str[BUFSIZ], *shape_str_i;
	char *attr[BUFSIZ];
	int i = 0, j = 0, k = 0, merged = 0, sum_exon_len = 0, acc_exon_len = 0,
			acc_len = 0, slot_index = 0;
	// List all transcripts that should be assembled.
	// If it is to compute the original tx info, just set it to the full list
	if (compute_ori) {
		target_ta = info->txs;
	} else {
		target_ta = g_ptr_array_sized_new(N_CTGS);
		while (fgets(buf, sizeof(buf), s_read_fp)) {
			chomp(buf);
			for (i = 0; i < info->txs->len; i++) {
				t = g_ptr_array_index(info->txs, i);
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
		for (j = 0; j < info->tx_seqs->len; j++) {
			tx_seq = g_ptr_array_index(info->tx_seqs, j);
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
				info->n_graph_egs++;
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

	info->n_uni_base = sum_exon_len;
	g_ptr_array_sort(cpns, (GCompareFunc) cmp_exons_by_len);
	for (i = 0; i < cpns->len; i++) {
		e = g_ptr_array_index(cpns, i);
		acc_exon_len += e->len;
		//show_debug_msg(__func__, "%s: %d => %d / %d\n", e->label, e->len,
		//		acc_exon_len, sum_exon_len);
		// The n50 value only set once
		if (acc_exon_len * 2 >= sum_exon_len && !info->best_n50) {
			info->best_n50 = e->len;
		}
		slot_index = (e->len >= MAX_LEN_FOR_PLOT) ? (MAX_LEN_FOR_PLOT
				/ SLOT_SIZE) : e->len / SLOT_SIZE;
		info->slots[slot_index] += e->len;
	}

	fclose(s_read_fp);
	fclose(graph_fp);
	fclose(tx_egs);
}

void cal_opt_n50(tx_info *info) {
	info->opt_n50 = cal_seq_n50(info->tx_seqs);
}

void trim(char *str) {
	int i = 0;
	for (i = 0; i < strlen(str); i++) {
		if(str[i] == '\n') {
			str[i] = '\0';
			return;
		}
	}
}

void parse_sam(rs_info *info) {
	char buf[BUFSIZ];
	char *id_str, *attr[16], *r_attr[5];
	int i = 0, q_len = 0, contig_id = 0;
	occarray *o_arr_i = g_ptr_array_sized_new(16);
	FILE *sam_fp = xopen(sam_fn, "r");
	FILE *res_fp = xopen(result_fn, "w");
	eva_occ *o = 0;
	edge *eg = 0;
	bwa_seq_t *ctg = 0;
	while (fgets(buf, sizeof(buf), sam_fp)) {
		i = 0;
		trim(buf);
		// printf("%s", buf);
		if (strstr(buf, "#") != NULL) {
			if (strstr(buf, "Query") != NULL) {
				show_debug_msg(__func__, "%s\n", buf);
				i = 0;
				attr[0] = strtok(buf, " ");
				while (attr[i] != NULL) { //ensure a pointer was found
					attr[++i] = strtok(NULL, " "); //continue to tokenize the string
				}

				for (i = 0; i < info->contigs->len; i++) {
					ctg = g_ptr_array_index(info->contigs, i);
					if (strcmp(ctg->name, attr[2]) == 0) {
						q_len = ctg->len;
						break;
					}
				}

				//contig_id = atoi(attr[2]);
				eg = new_eg();
				//eg->id = contig_id;
				eg->name = strdup(attr[2]);
				eg->len = q_len;
				g_ptr_array_add(info->ea_all, eg);
			}
			if (strstr(buf, " 0 hits found") != NULL) {
				info->n_not_aligned++;
				info->n_base_not_aligned += q_len;
				g_ptr_array_add(info->ea_not_aligned, eg);
			}
			continue;
		}
		o = new_occ();
		o->q_len = q_len;
		chomp(buf);
		//printf("---------------------------------------\n%s \n", buf);
		attr[0] = strtok(buf, "\t");
		i = 0;
		while (attr[i] != NULL) { //ensure a pointer was found
			attr[++i] = strtok(NULL, "\t"); //continue to tokenize the string
		}
		o->q_id = strdup(attr[0]);
		id_str = attr[1];
		i = 0;
		//		r_attr[0] = strtok(id_str, "|");
		//		while (r_attr[i] != NULL) { //ensure a pointer was found
		//			r_attr[++i] = strtok(NULL, "|"); //continue to tokenize the string
		//		}
		o->r_id = strdup(id_str);
		//		i = 0;
		//		r_attr[0] = strtok(r_attr[2], "=");
		//		while (r_attr[i] != NULL) { //ensure a pointer was found
		//			r_attr[++i] = strtok(NULL, "="); //continue to tokenize the string
		//		}
		//		o->r_len = atoi(r_attr[1]);
		o->percentage = atof(attr[2]) / 100;
		o->ali_len = atoi(attr[3]);
		o->start = atoi(attr[6]);
		o->end = atoi(attr[7]);
		if (o->start > o->end) {
			o->start = atoi(attr[7]);
			o->end = atoi(attr[6]);
		}
		o->evalue = atoi(attr[8]);
		//p_occ(o);
		if (o->ali_len <= 60)
			continue;

		g_ptr_array_add(info->occ_all, o);

		o_arr_i = g_hash_table_lookup(info->hits, o->r_id);
		if (o_arr_i == NULL) {
			o_arr_i = g_ptr_array_sized_new(16);
			g_ptr_array_add(o_arr_i, o);
			g_hash_table_insert(info->hits, o->r_id, o_arr_i);
		} else {
			g_ptr_array_add(o_arr_i, o);
		}
		//		if (counter++ > 10)
		//			break;
	}
	fclose(sam_fp);
	fclose(res_fp);
}

void read_contigs(rs_info *info) {
	bwa_seqio_t *ks;
	int i = 0;
	bwa_seq_t *seqs, *s;
	ks = bwa_open_reads(BWA_MODE, contigs_fn);
	while ((seqs = bwa_read_seq(ks, 0xa00000, &(info->n_ctgs), BWA_MODE, 0))
			!= 0) {
		for (i = 0; i < info->n_ctgs; i++) {
			s = &seqs[i];
			info->max_len = (info->max_len > s->len) ? info->max_len : s->len;
			g_ptr_array_add(info->contigs, s);
			info->n_base += s->len;
		}
		break;
	}
}

void cal_ass_n50(rs_info *info) {
	int i = 0, acc_len = 0, slot_index = 0, acc_ali_len = 0;
	bwa_seq_t *s;
	eva_occ *occ = 0;
	g_ptr_array_sort(info->contigs, (GCompareFunc) cmp_seqs_by_len);
	for (i = 0; i < info->n_ctgs; i++) {
		s = g_ptr_array_index(info->contigs, i);
		acc_len += s->len;
		show_debug_msg(__func__, "[%s] %d: %d / %d \n", s->name, s->len,
				acc_len, ori_info->n_uni_base);
		if (acc_len * 2 >= info->n_base && !info->n50) {
			info->n50 = s->len;
		}
		if (acc_len * 2 >= ori_info->n_base && !info->n50_all) {
			info->n50_all = s->len;
		}
		slot_index = (s->len >= MAX_LEN_FOR_PLOT) ? (MAX_LEN_FOR_PLOT
				/ SLOT_SIZE) : s->len / SLOT_SIZE;
		info->slots[slot_index] += s->len;
	}
	info->ave_len = acc_len / info->n_ctgs;
	info->base_coverage = (float) info->n_base_shared / ori_info->n_base;

	g_ptr_array_sort(info->occ_all, (GCompareFunc) cmp_occ_by_len);
	for (i = 0; i < info->occ_all->len; i++) {
		occ = g_ptr_array_index(info->occ_all, i);
		acc_ali_len += occ->ali_len;
		if (acc_ali_len * 2 >= ori_info->n_base) {
			info->n50_aligned = occ->ali_len;
			break;
		}
	}
}

void get_ori_info() {
	tx_info *info = init_info();
	// Read the sequences of transcripts, get total base and # of txs
	read_tx_seqs(info);
	// Read the exon meta info, get # of exons and exon list
	//read_exon_info(info);
	// Read the exon list of txs
	//read_tx_info(info);
	// Calculate the optimal N50 value
	cal_opt_n50(info);
	// Draw the graph and calculate the best N50
	//draw_graph(info);
	write_tx_sum(info);
	ori_info = info;
}

void get_ass_info() {
	int i = 0, len = 0, j = 0;
	char name[BUFSIZ];
	bwa_seq_t *t;

	occarray *o_arr_i = 0;
	rs_info *info = init_rs_info();
	read_contigs(info);
	parse_sam(info);
	g_hash_table_foreach(info->hits, (GHFunc) occ_c_iter, info);
	//g_hash_table_foreach(info->hits, (GHFunc) occ_p_iter, info);
	for (i = 1; i < 3874; i++) {
		sprintf(name, "CUFF.%d.1", i);
		o_arr_i = g_hash_table_lookup(info->hits, name);
		if (o_arr_i == NULL) {
			for (j = 0; j < ori_info->tx_seqs->len; j++) {
				t = g_ptr_array_index(ori_info->tx_seqs, j);
				if (strcmp(t->name, name) == 0) {
					len = t->len;
					break;
				}
			}
			printf("%s: no touched! length %d\n", name, len);
		} else {
			occ_p_iter(name, o_arr_i, NULL);
		}
	}
	cal_ass_n50(info);
	write_ass_rs(info);
}

int eva_main(int argc, char *argv[]) {
	clock_t t = clock();
	int c;
	while ((c = getopt(argc, argv, "o:m:e:g:t:c:s:r:")) >= 0) {
		switch (c) {
		case 'o':
			compute_ori = atoi(optarg);
			break;
		case 'm':
			sam_fn = optarg;
			break;
		case 'e':
			exons_fn = optarg;
			break;
		case 'g':
			edges_fn = optarg;
			break;
		case 't':
			tx_fn = optarg;
			break;
		case 'c':
			contigs_fn = optarg;
			break;
		case 's':
			start_reads_fn = optarg;
			break;
		case 'r':
			result_fn = optarg;
			break;
		default:
			compute_ori = atoi(optarg);
			break;
		}
	}
	get_ori_info();
	get_ass_info();
	show_msg(__func__, "Done: %.2f sec\n", (float) (clock() - t)
			/ CLOCKS_PER_SEC);
	return 0;
}
