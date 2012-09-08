/*
 * ass.h
 *
 *  Created on: Jun 5, 2011
 *      Author: Carl
 */
#ifndef ASS_H_
#define ASS_H_
#include "bwase.h"
#define MAX_BRANCH_LEVEL		32
#define THRE_USE_PAIRED			5
#define INIT_ARR_SIZE			16
#define N_BACK_LEVEL			3
#define MIN_LEN_BF_CHECK		3		// Critical param, 5 results in very complicated graph
#define MAX_LEN_FOR_MATES		3
#define INIT_EDGE_SPACE			1024
#define R_POOL_INIT_SIZE		1024
#define N_MIN_VAL_EXT			1
#define TLR_LEN					16
#define INVALID					-1
#define STOP_THRE				0.95
#define GRACE_NM				2
#define MIN_VALID_PAIRS			3

enum RESULT_TYPE
{
    NO_ALIGN,
    NOT_EXTEND,
    REP_EXTEND,
    REP_QUE,
    MUL_PATH,
    QUE_PFD
};

// To store the single extension result
typedef struct {
	enum RESULT_TYPE type; // 0: no alignment; 1: cannot extend; 2: repetitive query; 3: multiple path; 4: query performed
	char *message;
	void *counter;
	bwa_seq_t *query;
	bwa_seq_t *used;
} ext_msg;

int pe_ass(int argc, char *argv[]);

#endif /* ASS_H_ */
