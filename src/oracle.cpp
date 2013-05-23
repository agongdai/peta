/*
 * oracle.cpp
 *
 *  Created on: 21-May-2013
 *      Author: carl
 */

#include <vector>
#include <unordered_map>
#include <cassert>
#include <unistd.h>
#include <fstream>
#include <string>
#include <ctime>
#include <algorithm>
#include <stdint.h>
#include <glib.h>
#include <time.h>
#include <math.h>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>
#include <string.h>
#include <stdio.h>
#include "kmers.hpp"
#include "utils.h"
#include "peseq.h"
#include "rnaseq.h"
#include "bwtaln.h"
#include "pechar.h"
#include "hash.hpp"
#include "oracle.hpp"

float tx_is_expressed(bwa_seq_t *transcript, hash_map *hm) {
	uint64_t i = 0, n_blank = 0, n_ctu_blank = 0;
	bwa_seq_t *part = NULL;
	GPtrArray *hits = NULL;
	float rpkm = 0;
	float n_reads = 0;
	for (i = 0; i <= transcript->len - hm->o->read_len; i++) {
		part = new_seq(transcript, hm->o->read_len, i);
		//p_query(__func__, part);
		hits = align_full_seq(part, hm, 2);
		//p_readarray(hits, 1);
		if (hits->len == 0) {
			n_ctu_blank++;
			n_blank++;
		} else {
			n_ctu_blank = 0;
			n_reads += hits->len;
		}
		g_ptr_array_free(hits, TRUE);
		if (n_ctu_blank >= hm->o->read_len * 2 - 5) {
			show_debug_msg(__func__, "Transcript %s is not expressed! \n",
					transcript->name);
			return -1;
		}
	}
	show_debug_msg(__func__, "n_reads: %f; total reads: %d; tx len: %d\n",
			n_reads, hm->n_reads, transcript->len);
	rpkm = (float) (pow(10, 9) * n_reads) / ((float) (hm->n_reads
			* transcript->len));
	show_debug_msg(__func__, "Transcript %s is expressed! \n", transcript->name);
	return rpkm;
}

int oracle_set(int argc, char *argv[]) {
	int c;
	float rpkm = 0.0;
	uint32_t n_tx = 0, i = 0;
	mer_hash map;
	bwa_seq_t *transcripts = NULL, *tx = NULL;
	hash_map *hm = NULL;
	FILE *oracle_file = NULL;
	char fn[BUFSIZ], header[BUFSIZ];
	while ((c = getopt(argc, argv, "")) >= 0) {
	}
	if (optind + 1 > argc) {
		return 1;
	}
	transcripts = load_reads(argv[optind], &n_tx);
	hm = load_hash_map(argv[optind + 1], 0, map);
	hm->hash = &map;
	sprintf(fn, "%s.oracle", argv[optind]);
	oracle_file = xopen(fn, "w");
	for (i = 0; i < n_tx; i++) {
		tx = &transcripts[i];
		rpkm = tx_is_expressed(tx, hm);
		p_query(__func__, tx);
		if (rpkm > 0) {
			sprintf(header, ">%s rpkm:%.2f length:%d\n", tx->name, rpkm,
					tx->len);
			save_con(header, tx, oracle_file);
		}
	}
	fclose(oracle_file);
	return 0;
}

/**
 * Read in the exons in txt format: chr  start  end  weight
 */
GPtrArray *read_regions(char *fn) {
	FILE *exon_f = xopen(fn, "r");
	region *r = NULL;
	char buf[BUFSIZ];
	int i = 0;
	GPtrArray *exons = g_ptr_array_sized_new(BUFSIZ);
	char *attr[64];
	while (fgets(buf, sizeof(buf), exon_f)) {
		attr[0] = strtok(buf, "\t");
		while (attr[i] != NULL) {
			attr[++i] = strtok(NULL, "\t");
		}
		r = (region*) malloc(sizeof(region));
		r->chr = strdup(attr[0]);
		r->start = atoi(attr[1]);
		r->end = atoi(attr[2]);
		r->weight = atof(attr[3]);
		g_ptr_array_add(exons, r);
		i = 0;
	}
	fclose(exon_f);
	return exons;
}

void save_exons(GPtrArray *exons, char *fn) {
	FILE *fp = xopen(fn, "w");
	region *e = NULL;
	uint32_t i = 0;
	char line[BUFSIZ];
	for (i = 0; i < exons->len; i++) {
		e = (region*) g_ptr_array_index(exons, i);
		sprintf(line, "%s\t%d\t%d\t%.2f", e->chr, e->start, e->end, e->weight);
		fputs(line, fp);
	}
	fclose(fp);
}

void save_splicings(GPtrArray *splicings, char *fn) {
	FILE *fp = xopen(fn, "w");
	splicing *s = NULL;
	uint32_t i = 0;
	char line[BUFSIZ];
	for (i = 0; i < splicings->len; i++) {
		s = (splicing*) g_ptr_array_index(splicings, i);
		sprintf(line, "%s\t%d\t%d\t%.2f", s->chr, s->left->end, s->right->start, s->weight);
		fputs(line, fp);
	}
	fclose(fp);
}

int genome_splicings(int argc, char *argv[]) {
	int c;
	uint64_t kmer_int = 0;
	bwa_seq_t *chromosomes = NULL;
	uint32_t n_chr = 0, i = 0, j = 0;
	mer_hash map;
	hash_map *hm = NULL;
	GPtrArray *exons = NULL;
	region *r = NULL;
	while ((c = getopt(argc, argv, "")) >= 0) {
	}
	if (optind + 2 > argc) {
		return 1;
	}
	chromosomes = load_reads(argv[optind], &n_chr);
	exons = read_regions(argv[optind + 1]);
//	hm = load_hash_map(argv[optind + 2], 0, map);
	for (i = 0; i < exons->len; i++) {
		r = (region*) g_ptr_array_index(exons, i);
		show_debug_msg(__func__, "%s\t%d\t%d\t%.2f\n", r->chr, r->start, r->end, r->weight);
	}
	return 0;
}
