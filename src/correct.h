/*
 * correct.h
 *
 *  Created on: 06-Jan-2013
 *      Author: carl
 */

#ifndef CORRECT_H_
#define CORRECT_H_
#include <glib.h>
#include "pehash.h"

int rm_repetitive_reads(bwa_seq_t *seqs, const int n_seqs);
int correct_reads(hash_table *ht, const int n_threads);
GPtrArray *get_low_kmer_reads(hash_table *ht, const int k, const double thre);

#endif /* CORRECT_H_ */
