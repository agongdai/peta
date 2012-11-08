/*
 * peeva.h
 *
 *  Created on: Jan 1, 2012
 *      Author: Cai Shaojiang
 */

#ifndef PEEVA_H_
#define PEEVA_H_

#define N_NOT_ALIGNED_CTGS	128
#define N_CTGS				4096
#define NOT_FOUND		-1

#include <glib.h>

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

int pe_eva(int argc, char *argv[]);

#endif /* PEEVA_H_ */
