/*
 * kmer.h
 *
 *  Created on: 13-Mar-2013
 *      Author: carl
 */

#ifndef KMER_H_
#define KMER_H_

#include "bwase.h"

typedef struct {
	bwa_seq_t *s;
	int fre;
} mer;

#endif /* KMER_H_ */
