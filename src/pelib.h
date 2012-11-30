#ifndef PELIB_H_
#define PELIB_H_

#define MISMATCHES			1
#define STRICT_PERC			0.8
#define PAIRS_PER_EDGE		100
#define MATE_OVERLAP_THRE	11
#define RELAX_MATE_OL_THRE	6
#define MAX_SINGLE_EDGES	100
#define SINGLE_EDGE_THRE	150

int pe_lib(int argc, char *argv[]);

#endif /* PELIB_H_ */
