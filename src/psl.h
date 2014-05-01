/*
 * psl.h
 *
 *  Created on: May 1, 2014
 *      Author: carl
 */

#ifndef PSL_H_
#define PSL_H_

#include <stdint.h>
#include <stdio.h>
#include <inttypes.h>

/**
 * PSL format: http://asia.ensembl.org/info/website/upload/psl.html
	1. matches - Number of matching bases that aren't repeats.
	2. misMatches - Number of bases that don't match.
	3. repMatches - Number of matching bases that are part of repeats.
	4. nCount - Number of 'N' bases.
	5. qNumInsert - Number of inserts in query.
	6. qBaseInsert - Number of bases inserted into query.
	7. tNumInsert - Number of inserts in target.
	8. tBaseInsert - Number of bases inserted into target.
	9. strand - defined as + (forward) or - (reverse) for query strand. In mouse, a second '+' or '-' indecates genomic strand.
	10. qName - Query sequence name.
	11. qSize - Query sequence size.
	12. qStart - Alignment start position in query.
	13. qEnd - Alignment end position in query.
	14. tName - Target sequence name.
	15. tSize - Target sequence size.
	16. tStart - Alignment start position in query.
	17. tEnd - Alignment end position in query.
	18. blockCount - Number of blocks in the alignment.
	19. blockSizes - Comma-separated list of sizes of each block.
	20. qStarts - Comma-separated list of start position of each block in query.
	21. tStarts - Comma-separated list of start position of each block in target.
**/
typedef struct {
	uint32_t matches;
	uint32_t misMatches;
	uint32_t repMatches;
	uint16_t nCount;
	uint8_t qNumInsert;
	uint16_t qBaseInsert;
	uint8_t tNumInsert;
	uint16_t tBaseInsert;
	char strand;
	char *qName;
	uint16_t qSize;
	uint32_t qStart;
	uint32_t qEnd;
	char *tName;
	uint32_t tSize;
	uint32_t tStart;
	uint32_t tEnd;
	uint8_t blockCount;
	uint32_t *blockSizes;
	uint32_t *qStarts;
	uint32_t *tStarts;
} hit;

void read_hits(char *fn, GPtrArray *hits, int load_self);

#endif /* PSL_H_ */
