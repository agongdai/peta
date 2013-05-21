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

int tx_is_expressed(bwa_seq_t *transcript, hash_map *hm) {
	int i = 0, n_blank = 0, n_ctu_blank = 0;
	bwa_seq_t *part = NULL;
	GPtrArray *hits = NULL;
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
		}
		g_ptr_array_free(hits, TRUE);
		if (n_ctu_blank >= hm->o->read_len * 2 - 5) {
			show_debug_msg(__func__, "Transcript %s is not expressed! \n", transcript->name);
			return 0;
		}
	}
	show_debug_msg(__func__, "Transcript %s is expressed! \n", transcript->name);
	return 1;
}

int oracle_set(int argc, char *argv[]) {
	int c;
	uint32_t n_tx = 0, i = 0;
	mer_hash map;
	bwa_seq_t *transcripts = NULL, *tx = NULL;
	hash_map *hm = NULL;
	while ((c = getopt(argc, argv, "")) >= 0) {
	}
	if (optind + 1 > argc) {
		return 1;
	}
	transcripts = load_reads(argv[optind], &n_tx);
	hm = load_hash_map(argv[optind + 1], 0, map);
	hm->hash = &map;
	for (i = 0; i < n_tx; i++) {
		tx = &transcripts[i];
		p_query(__func__, tx);
		tx_is_expressed(tx, hm);
//		if (i > 10)
//			break;
	}
	return 0;
}
