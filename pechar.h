/*
 * pechar.h
 *
 *  Created on: 20-Jun-2011
 *      Author: carl
 */

#ifndef PECHAR_H_
#define PECHAR_H_

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif
#define MULPATH 		0.25
#define LEAST_ALIGN 	0.5
#define ALIGH_THRE 		20
#define LEAST_ALIGN_N 	2
#define NOT_FOUND		-1
#define INVALID_CHAR	-1
#define N_NUCLITODE		4

typedef struct {
	unsigned int n_a;
	unsigned int n_c;
	unsigned int n_g;
	unsigned int n_t;
	unsigned int n_n;
} char_sta; // char statistics

char to_upper_lower(char c);
void reset_c(int *sta, int *c);
void check_c(int *sta, uint8_t c);
int sum_c(const int *sta);
int *get_most(const int *sta);
char compl(const char nt);
char mutate(const char c);
int *get_abs_most(const int *sta, const double threshold);

#ifdef __cplusplus
}
#endif
#endif /* PECHAR_H_ */
