/*
 * pechar.h
 *
 *  Created on: 20-Jun-2011
 *      Author: carl
 */

#ifndef PECHAR_H_
#define PECHAR_H_

#include <stdint.h>

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

#ifdef __cplusplus
extern "C" {
#endif

	char to_upper_lower(char c);
	void reset_c(int *sta, int *c);
	void check_c(int *sta, uint8_t c);
	int sum_c(const int *sta);
	int *get_most(const int *sta);
	char cpl(const char nt);
	char mutate(const char c);
	int *get_abs_most(const int *sta, const double threshold);
	int get_pure_most(const int *sta);
	uint64_t shift_bit_to_left(uint64_t kmer, const int new_c, const int n);
	uint64_t shift_bit_to_right(uint64_t kmer, const int new_c, const int n);
	uint64_t shift_bit(uint64_t kmer, const int new_c, const int n, const int ori);
	int get_max_index(const int *counters);
	int get_second_freq_char(const int *counters, const int max);

#ifdef __cplusplus
}
#endif
#endif /* PECHAR_H_ */
