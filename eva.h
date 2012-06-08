/*
 * eva.h
 *
 *  Created on: 16-Apr-2012
 *      Author: carl
 */

#ifndef EVA_H_
#define EVA_H_

#define N_NOT_ALIGNED_CTGS	128
#define N_CTGS				4096
#define NOT_FOUND			-1
#define ATTR_STR_LEN		32
#define SLOT_SIZE			500
#define MAX_LEN_FOR_PLOT	10000

#include <glib.h>
#include "bwase.h"

typedef GPtrArray occarray;
typedef GPtrArray txarray;
typedef GPtrArray exonarray;
typedef struct exon exon;

typedef struct {
	int start;
	int end;
	char *q_id;
	char *r_id;
	float percentage;
	double evalue;
	int ali_len;
	int q_len;
	int r_len;
} eva_occ;

typedef struct {
	char *name;
	exonarray *ea;
	int len;
	int touched; // Whether any portion of this transcript is assembled.
} tx;

struct exon {
	int id;
	int len;
	char *label;
	txarray *ta;
	exon *merged_to;
	int drawed;
};

typedef struct {
	int n_base;
	int n_uni_base;
	int n_tx;
	int n_exon;
	int opt_n50;
	int best_n50;
	int n_sd_tx;
	int n_sd_base;
	int n_graph_egs;
	float base_sd_total;
	char *graph_fn;
	char *sum_fn;
	char *cpn_fn;
	int *slots;
	int n_slot;
	GPtrArray *txs;
	GPtrArray *sd_txs;
	GPtrArray *exons;
	GPtrArray *tx_seqs;
	GPtrArray *exon_seqs;
} tx_info;

typedef struct {
	int n_base_shared;
	int n_base;
	int n_ctgs;
	int n_full_len;
	int n_sd;
	int max_len;
	int n_not_aligned;
	int *slots;
	int n_slot;
	GPtrArray *contigs;
	GPtrArray *ea_all;
	GPtrArray *ea_full_len;
	GPtrArray *ea_sd;
	GPtrArray *ea_not_aligned;
	GPtrArray *ea_not_touched;
	GPtrArray *occ_all;
	int n_base_not_aligned;
	int n50;
	int n50_all;
	int n50_aligned;
	int ave_len;
	float base_coverage;
	// Key is the transcript name, value is an array of pointers, each is a hit
	GHashTable *hits;
} rs_info;

int eva_main(int argc, char *argv[]);

#endif /* EVA_H_ */
