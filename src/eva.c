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
#include "rnaseq.h"

int compute_ori = 0;
char *sam_fn;
char *tx_fn;
char *contigs_fn;
char *result_dir;
tx_info *ori_info;

int tmp_total_base = 0;
int tmp_got_base = 0;
int tmp_n_cufflinks = 0;
GPtrArray *one_on_one_ids = NULL;
GPtrArray *full_length_ids = NULL;
GPtrArray *one_covered_ids = NULL;

char *get_result_file(const char *file_name) {
	char *name = (char*) calloc(512, sizeof(char));
	strcat(name, result_dir);
	strcat(name, "/");
	strcat(name, file_name);
	return name;
}

tx_info *init_info() {
	tx_info *i = (tx_info*) malloc(sizeof(tx_info));
	i->base_sd_total = 0.0;
	i->best_n50 = 0;
	i->n_base = 0;
	i->n_uni_base = 0;
	i->n_exon = 0;
	i->n_graph_egs = 0;
	i->n_sd_base = 0;
	i->n_sd_tx = 0;
	i->n_tx = 0;
	i->opt_n50 = 0;
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
	i->n_70_covered = 0;
	i->n_one_on_one = 0;
	i->n_one_covered = 0;
	i->n_not_aligned = 0;
	i->n_not_reached = 0;
	i->n_base_not_reached = 0;
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
	n_base_got++;
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
	int one_covered = 0, n_bases_on_contig = 0, is_full_length = 0,
			is_70_covered = 0, is_one_on_one = 0;
	g_ptr_array_sort(occs, (GCompareFunc) cmp_occ);
	bwa_seq_t *seq = NULL;

	for (i = 0; i < ori_info->tx_seqs->len; i++) {
		seq = g_ptr_array_index(ori_info->tx_seqs, i);
		//show_debug_msg(__func__, "seq name: %s \n", seq->name);
		//show_debug_msg(__func__, "align name: %s \n", (char*) key);
		if (strcmp(seq->name, (char*) key) == 0) {
			tx_len = seq->len;
			break;
		}
	}

	// Check whether is one-on-one transcript
	//	if (occs->len == 1) {
	for (i = 0; i < occs->len; i++) {
		o = g_ptr_array_index(occs, i);
		if (o->end >= o->start && o->ali_len > (tx_len * 0.9)) {
			if (o->ali_len > (o->q_len * 0.9)) {
				info->n_one_on_one++;
				is_one_on_one = 1;
				show_debug_msg(
						"==============================================",
						"ONE-ON-ONE\n");
				occ_p_iter(key, value, NULL);
				g_ptr_array_add(one_on_one_ids, (char *) key);
				break;
			}
		}
	}
	// Check whether is full-length
	if (!is_one_on_one) {
		for (i = 0; i < occs->len; i++) {
			o = g_ptr_array_index(occs, i);
			if (o->end >= o->start && o->ali_len > (tx_len * 0.9)) {
				info->n_full_len++;
				is_full_length = 1;
				show_debug_msg(
						"==============================================",
						"FULL-LENGTH\n");
				occ_p_iter(key, value, NULL);
				g_ptr_array_add(full_length_ids, (char *) key);
				break;
			}
		}
	}

	if (!is_one_on_one && !is_full_length) {
		for (i = 0; i < occs->len; i++) {
			o = g_ptr_array_index(occs, i);
			if (o->end >= o->start && o->ali_len > (tx_len * 0.7)) {
				if (o->ali_len > (o->q_len * 0.9)) {
					info->n_70_covered++;
					is_70_covered = 1;
					show_debug_msg(
							"==============================================",
							"70-COVERED\n");
					occ_p_iter(key, value, NULL);
					break;
				}
			}
		}
	}
	//	}

	one_covered = 1;
	for (i = 0; i < occs->len; i++) {
		o = g_ptr_array_index(occs, i);
		if (o->end < o->start) {
			o_pre = o;
			continue;
		}

		if (i == 0) {
			info->n_base_shared += o->end - o->start + 1;
			n_bases_on_contig += o->end - o->start + 1;
		} else {
			if (strcmp(o->q_id, o_pre->q_id) != 0 && one_covered)
				one_covered = 0;
			// To avoid duplicate counting:
			// NM_001042362: [1, 771] [1, 770] [772, 1709] [772, 1709] [1709, 2014] [1709, 2014]
			if (o->end <= o_pre->end) {
				o_pre = o;
				continue;
			}
			if (o->start < o_pre->end) {
				info->n_base_shared += o->end - o_pre->end;
				n_bases_on_contig += o->end - o_pre->end;
			} else {
				info->n_base_shared += o->end - o->start + 1;
				n_bases_on_contig += o->end - o->start + 1;
			}
		}
		o_pre = o;
	}
	if (one_covered && tx_len * 0.9 <= n_bases_on_contig) {
		if (occs->len > 1) {
			info->n_one_covered++;
			show_debug_msg("==============================================",
					"ONE-COVERED\n");
			occ_p_iter(key, value, user_data);
			if (!is_one_on_one && !is_full_length && !is_70_covered) {
				g_ptr_array_add(one_covered_ids, (char *) key);
			}
		}
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

gint cmp_ctgs(gpointer a, gpointer b) {
	edge* eg_a = *(edge**) a;
	edge* eg_b = *(edge**) b;
	return eg_b->len - eg_a->len;
}

void write_ass_rs(rs_info *info) {
	char item[BUFSIZ], *str;
	char *name = NULL;
	FILE *rs_fp = NULL;
	int i = 0;
	edge *t = 0;

	name = get_result_file("result.txt");
	rs_fp = xopen(name, "w");
	free(name);

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

	str = fix_len("# of 70% covered: ", ATTR_STR_LEN);
	sprintf(item, "%s\t%d\n", str, info->n_70_covered);
	fputs(item, rs_fp);

	str = fix_len("# of one-on-one transcripts: ", ATTR_STR_LEN);
	sprintf(item, "%s\t%d\n", str, info->n_one_on_one);
	fputs(item, rs_fp);

	str = fix_len("Transcripts covered by one contig: ", ATTR_STR_LEN);
	sprintf(item, "%s\t%d\n", str, info->n_one_covered);
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

	str = fix_len("Transcripts not reached: ", ATTR_STR_LEN);
	sprintf(item, "%s\t%d\n", str, info->n_not_reached);
	fputs(item, rs_fp);

	str = fix_len("Bases not reached: ", ATTR_STR_LEN);
	sprintf(item, "%s\t%d\n", str, info->n_base_not_reached);
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
	int i = 0;
	bwa_seq_t *seqs, *s;
	seqs = load_reads(tx_fn, &(info->n_tx));
	for (i = 0; i < info->n_tx; i++) {
		s = &seqs[i];
		info->n_base += s->len;
		g_ptr_array_add(info->tx_seqs, s);
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

void cal_opt_n50(tx_info *info) {
	info->opt_n50 = cal_seq_n50(info->tx_seqs);
}

void trim(char *str) {
	int i = 0;
	for (i = 0; i < strlen(str); i++) {
		if (str[i] == '\n') {
			str[i] = '\0';
			return;
		}
	}
}

void parse_sam(rs_info *info) {
	char buf[BUFSIZ];
	char *id_str, *attr[16];
	int i = 0, q_len = 0;
	occarray *o_arr_i = g_ptr_array_sized_new(16);
	FILE *sam_fp = xopen(sam_fn, "r");
	eva_occ *o = 0;
	edge *eg = 0;
	bwa_seq_t *ctg = 0;
	while (fgets(buf, sizeof(buf), sam_fp)) {
		i = 0;
		trim(buf);
		// printf("%s", buf);
		if (strstr(buf, "#") != NULL) {
			if (strstr(buf, "Query") != NULL) {
				//show_debug_msg(__func__, "%s\n", buf);
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
		if (o->ali_len <= 100)
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
}

void read_contigs(rs_info *info) {
	int i = 0;
	bwa_seq_t *seqs, *s;
	seqs = load_reads(contigs_fn, &(info->n_ctgs));
	pe_reverse_seqs(seqs, info->n_ctgs);
	for (i = 0; i < info->n_ctgs; i++) {
		s = &seqs[i];
		info->max_len = (info->max_len > s->len) ? info->max_len : s->len;
		g_ptr_array_add(info->contigs, s);
		info->n_base += s->len;
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

void write_not_aligned_ctgs(rs_info *info) {
	int i = 0, j = 0;
	bwa_seq_t *contig = 0;
	edge *eg = 0;
	FILE *not_aligned = NULL;
	char header[BUFSIZ], *name = NULL;

	name = get_result_file("not_aligned.fa");
	not_aligned = xopen(name, "w");
	for (i = 0; i < info->contigs->len; i++) {
		contig = g_ptr_array_index(info->contigs, i);
		for (j = 0; j < info->ea_not_aligned->len; j++) {
			eg = g_ptr_array_index(info->ea_not_aligned, j);
			if (strcmp(eg->name, contig->name) == 0) {
				sprintf(header, ">%s len=%d\n", contig->name, contig->len);
				save_con(header, contig, not_aligned);
				break;
			}
		}
	}
	fclose(not_aligned);
}

void get_ori_info() {
	tx_info *info = init_info();
	// Read the sequences of transcripts, get total base and # of txs
	read_tx_seqs(info);
	// Calculate the optimal N50 value
	cal_opt_n50(info);
	ori_info = info;
}

void write_hit_txs() {
	FILE *one_on_one_file = NULL;
	FILE *full_length_file = NULL;
	FILE *one_covered_file = NULL;
	char item[256], *hit_tx_id = NULL, *name = NULL;
	int i = 0;

	name = get_result_file("one_on_one.txt");
	one_on_one_file = xopen(name, "w");
	free(name);
	name = get_result_file("full_length.txt");
	full_length_file = xopen(name, "w");
	free(name);
	name = get_result_file("one_covered.txt");
	one_covered_file= xopen(name, "w");
	free(name);

	for (i = 0; i < one_on_one_ids->len; i++) {
		hit_tx_id = g_ptr_array_index(one_on_one_ids, i);
		sprintf(item, "%s\n", hit_tx_id);
		fputs(item, one_on_one_file);
	}
	for (i = 0; i < full_length_ids->len; i++) {
		hit_tx_id = g_ptr_array_index(full_length_ids, i);
		sprintf(item, "%s\n", hit_tx_id);
		fputs(item, full_length_file);
	}
	for (i = 0; i < one_covered_ids->len; i++) {
		hit_tx_id = g_ptr_array_index(one_covered_ids, i);
		sprintf(item, "%s\n", hit_tx_id);
		fputs(item, one_covered_file);
	}

	fclose(one_on_one_file);
	fclose(full_length_file);
	fclose(one_covered_file);
}

void get_ass_info() {
	int i = 0;
	char *name;
	bwa_seq_t *s;

	occarray *o_arr_i = 0;
	rs_info *info = init_rs_info();
	read_contigs(info);
	parse_sam(info);
	one_on_one_ids = g_ptr_array_sized_new(256);
	full_length_ids = g_ptr_array_sized_new(256);
	one_covered_ids = g_ptr_array_sized_new(256);
	g_hash_table_foreach(info->hits, (GHFunc) occ_c_iter, info);
	write_hit_txs();
	//g_hash_table_foreach(info->hits, (GHFunc) occ_p_iter, info);
	for (i = 1; i < ori_info->tx_seqs->len; i++) {
		s = g_ptr_array_index(ori_info->tx_seqs, i);
		name = s->name;
		o_arr_i = g_hash_table_lookup(info->hits, name);
		if (o_arr_i == NULL) {
			info->n_not_reached++;
			info->n_base_not_reached += s->len;
			printf("%s: no touched! length %d\n", name, s->len);
		} else {
			occ_p_iter(name, o_arr_i, NULL);
		}
	}
	cal_ass_n50(info);
	write_not_aligned_ctgs(info);
	write_ass_rs(info);
}

int eva_main(int argc, char *argv[]) {
	clock_t t = clock();
	int c;
	while ((c = getopt(argc, argv, "o:m:t:c:r:")) >= 0) {
		switch (c) {
		case 'o':
			compute_ori = atoi(optarg);
			break;
		case 'm':
			sam_fn = optarg;
			break;
		case 't':
			tx_fn = optarg;
			break;
		case 'c':
			contigs_fn = optarg;
			break;
		case 'r':
			result_dir = optarg;
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
