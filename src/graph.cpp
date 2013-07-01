/*
 * graph.cpp
 *
 *  Created on: Jun 30, 2013
 *      Author: carl
 */

#include <unordered_map>
#include <stdint.h>
#include <glib.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>
#include "bwtaln.h"
#include "edge.h"
#include "kmers.hpp"
#include "ass.hpp"
#include "rnaseq.h"
#include "pechar.h"
#include "utils.h"
#include "peseq.h"
#include "junction.hpp"
#include "hash.hpp"

void assign_reads2tpl(edge *eg, hash_map *hm) {
	int j = 0, i = 0, read_len = hm->o->read_len;
	bwa_seq_t *seq = NULL, *read = NULL;
	GPtrArray *hits = NULL;

	for (i = 0; i <= eg->len - read_len; i++) {
		seq = new_seq(eg->ctg, read_len, i);
		hits = align_full_seq(seq, hm, 2);
		for (j = 0; j < hits->len; j++) {
			read = (bwa_seq_t*) g_ptr_array_index(hits, j);
			add_read2edge(eg, read);
		}
		bwa_free_read_seq(1, seq);
		g_ptr_array_free(hits, TRUE);
	}
}


