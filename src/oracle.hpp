/*
 * oracle.hpp
 *
 *  Created on: 21-May-2013
 *      Author: carl
 */

#ifndef ORACLE_HPP_
#define ORACLE_HPP_

#include <stdint.h>
#include <stdlib.h>
#include "bwtaln.h"

#define				EXON_OVERLAP	11

using namespace std;

// This is some region on the chromosome
typedef struct {
	char *chr;
	bwa_seq_t *ctg;
	int start;
	int end;
	float weight;
} region;

typedef struct {
	char *chr;
	region *left;
	region *right;
	float weight;
} splicing;

#ifdef __cplusplus
extern "C" {
#endif

	int oracle_set(int argc, char *argv[]);
	int genome_splicings(int argc, char *argv[]);

#ifdef __cplusplus
}
#endif

#endif /* ORACLE_HPP_ */
