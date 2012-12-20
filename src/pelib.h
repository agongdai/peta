#ifndef PELIB_H_
#define PELIB_H_

#define STRICT_PERC			0.8
#define PAIRS_PER_EDGE		100
#define MATE_OVERLAP_THRE	11
#define STRICT_MATCH_OL		8
#define RELAX_MATE_OL_THRE	8
#define MAX_SINGLE_EDGES	200
#define SINGLE_EDGE_THRE	150
#define N_BIG_MATE_POOL		1000
#define STOP_THRE_STAGE_1	0.8
#define STOP_THRE_STAGE_2	0.95
#define STOP_THRE_STAGE_3	0.99

int pe_lib(int argc, char *argv[]);

#endif /* PELIB_H_ */
