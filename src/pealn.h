/*
 * pealn.h
 *
 *  Created on: 26-Jun-2011
 *      Author: carl
 */

#ifndef PEALN_H_
#define PEALN_H_

#include <stdint.h>
#include "peseq.h"
#include "pehash.h"

#ifdef __cplusplus
extern "C" {
#endif

#define SHIFT 				6
#define N_DEFAULT_HITS		1024
#define N_DEFAULT_ALIGNS	128
#define MAX_SEQID_CHARS		31

typedef GPtrArray alignarray;

typedef struct {
	uint8_t nm; // Maximum edit distance
	uint8_t nsp; // Whether is not strand specific
	uint16_t ol; // Overlapping length
	uint16_t rl; // Read length
	uint16_t n_alg; // # of alternative alignments returned
	uint16_t mean; // Average span size of the RNA-seq lib
	uint16_t sd; // Standard deviation of the RNA-seq lib
	uint8_t pair; // Use pair-end info to extend
	char *solid_reads;
	char *out_root;
	int mode;
} ass_opt;

typedef struct {
	index64 q_id;
	index64 r_id;
	int pos;
	char strand;
	uint8_t diff;
	uint8_t rev_comp; // Whether is reverse complement
	uint32_t len;
	uint8_t flag;
} alg;

typedef struct {
	uint16_t n_algs; // # of total alignments
	alg *algs;
} alignment;

typedef struct {
	index64 seq_id;
	int shift; // Check SSAHA paper for meaning.
	int offset;
} pos_tuple;

void free_alg(alignarray *alns);
void reset_alg(alignarray *alns);
alg *aligned(const alignarray *alns, const index64 id);
void p_align(const alignarray *alns);
int pe_aln_test(int argc, char *argv[]);
void erase_reads_on_ht(hash_table *ht);
void pe_aln_query(const bwa_seq_t *query, const ubyte_t *q_seq, const hash_table *ht,
		const int mismatches, const int ol, const int rev, alignarray *aligns);

#ifdef __cplusplus
}
#endif

#endif /* PEALN_H_ */
