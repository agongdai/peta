/*
 * kmer.c
 *
 *  Created on: 13-Mar-2013
 *      Author: carl
 */

#include <glib.h>
#include "pehash.h"
#include "bwase.h"
#include "kmer.h"

GHashTable* count_frequency(hash_table *ht, const int k) {
	GHashTable* kmer_freq = g_hash_table_new(g_str_hash, g_int_equal);
	int i = 0;
	bwa_seq_t *read = NULL, *sub = NULL;
	mer *m = NULL;
}
