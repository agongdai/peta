/*
 * psl.c
 *
 *  Created on: May 1, 2014
 *      Author: carl
 */

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <glib.h>
#include <inttypes.h>
#include "psl.h"
#include "utils.h"

hit *new_psl_hit() {
	hit *h = (hit*) malloc(sizeof(hit));
	h->matches = 0;
	h->misMatches = 0;
	h->repMatches = 0;
	h->nCount = 0;
	h->qNumInsert = 0;
	h->qBaseInsert = 0;
	h->tNumInsert = 0;
	h->tBaseInsert = 0;
	h->strand = '+';
	h->qName = NULL;
	h->qSize = 0;
	h->qStart = 0;
	h->qEnd = 0;
	h->tName = NULL;
	h->tSize = 0;
	h->tStart = 0;
	h->tEnd = 0;
	h->blockCount = 0;
	h->blockSizes = NULL;
	h->qStarts = NULL;
	h->tStarts = NULL;
	return h;
}

void destroy_hit(hit *h) {
	if (!h) return;
	if (h->qName) free(h->qName);
	if (h->tName) free(h->tName);
	if (h->blockSizes) free(h->blockSizes);
	if (h->qStarts) free(h->qStarts);
	if (h->tStarts) free(h->tStarts);
	free(h);
}

/**
 * fn: PSL file with 5-line header
 * hits: not-NULL array to hold hits
 * load_self: whether to load hits where qName == tName
 */
void read_hits(char *fn, GPtrArray *hits, int load_self) {
	FILE *psl = xopen(fn, "r");
	char buf[BUFSIZ];
	int line_no = 0, i = 0, n = 0;
	char *attr[32];
	hit *h = NULL;

	while (fgets(buf, sizeof(buf), psl)) {
		if (++line_no <= 5) continue;

		// Read attributes and look up hashtable
		attr[0] = strtok(buf, "\t");
		for (i = 1; i < 21; i++) attr[i] = strtok(NULL, "\t");
		if (!load_self && strcmp(attr[9], attr[13]) == 0) continue;

		h = new_psl_hit();
		h->matches = atoi(attr[0]);
		h->misMatches = atoi(attr[1]);
		h->repMatches = atoi(attr[2]);
		h->nCount = atoi(attr[3]);
		h->qNumInsert = atoi(attr[4]);
		h->qBaseInsert = atoi(attr[5]);
		h->tNumInsert = atoi(attr[6]);
		h->tBaseInsert = atoi(attr[7]);
		h->strand = attr[8];
		h->qName = strdup(attr[9]);
		h->qSize = atoi(attr[10]);
		h->qStart = atoi(attr[11]);
		h->qEnd = atoi(attr[12]);
		h->tName = strdup(attr[13]);
		h->tSize = atoi(attr[14]);
		h->tStart = atoi(attr[15]);
		h->tEnd = atoi(attr[16]);
		h->blockCount = atoi(attr[17]);
		h->blockSizes = (uint32_t*) calloc(h->blockCount, sizeof(uint32_t));
		h->qStarts = (uint32_t*) calloc(h->blockCount, sizeof(uint32_t));
		h->tStarts = (uint32_t*) calloc(h->blockCount, sizeof(uint32_t));

		h->blockSizes[0] = atoi(strtok(attr[18], ","));
		for (i = 1; i < atoi(attr[17]); i++) {
			h->blockSizes[i] = atoi(strtok(NULL, ","));
		}
		h->qStarts[0] = atoi(strtok(attr[19], ","));
		for (i = 1; i < atoi(attr[17]); i++) {
			h->qStarts[i] = atoi(strtok(NULL, ","));
		}
		h->tStarts[0] = atoi(strtok(attr[20], ","));
		for (i = 1; i < atoi(attr[17]); i++) {
			h->tStarts[i] = atoi(strtok(NULL, ","));
		}
		g_ptr_array_add(hits, h); n++;
	}
	show_msg(__func__, "Loaded %s hits ...\n", n);
	fclose(psl);
}
