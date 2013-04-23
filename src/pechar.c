/*
 * pechar.c
 *
 *  Created on: 20-Jun-2011
 *      Author: carl
 */
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include "pechar.h"
#include "rand.h"
char nucle_u[] = { 'A', 'C', 'G', 'T' };
char nucle_l[] = { 'a', 'c', 'g', 't' };

char to_upper_lower(char c) {
	return tolower(c);
}

char mutate(const char c) {
	int rand = -1;
	char mut;
	while (1) {
		rand = (int) (rand_f() * 3);
		if (isupper(c)) {
			if (c != nucle_u[rand]) {
				mut = nucle_u[rand];
				break;
			}
		} else {
			if (c != nucle_l[rand]) {
				mut = nucle_l[rand];
				break;
			}
		}
	}
	//	fprintf(stderr, "Mutated: %c -> %c \n", c, mut);
	return mut;
}

char compl(const char nt) {
	switch (nt) {
	case 'a':
		return 't';
	case 'c':
		return 'g';
	case 'g':
		return 'c';
	case 't':
		return 'a';
	case 'A':
		return 'T';
	case 'C':
		return 'G';
	case 'G':
		return 'C';
	case 'T':
		return 'A';
	case 'N':
		return 'N';
	case 'n':
		return 'n';
	default:
		return 0;
	}
}

void reset_c(int *sta, int *c) {
	if (sta)
		sta[0] = sta[1] = sta[2] = sta[3] = sta[4] = 0;
	if (c)
		c[0] = c[1] = c[2] = c[3] = c[4] = INVALID_CHAR;
}

void check_c(int *sta, uint8_t c) {
	sta[c]++;
}

int sum_c(const int *sta) {
	return sta[0] + sta[1] + sta[2] + sta[3] + sta[4];
}

/**
 * Only get the char with the largest count
 */
int *get_abs_most(const int *sta, const double threshold) {
	long thre = 0;
	int index = 0, most = 0;
	int *status = (int*) calloc(2, sizeof(int));
	int total = sta[0] + sta[1] + sta[2] + sta[3] + sta[4];
	thre = total * threshold;
	for (index = 0; index < 4; index++) {
		if (sta[index] > sta[most])
			most = index;
	}
	status[0] = most;
	//	if (sta[most] < thre)
	//		status[0] = -1;
	status[1] = most;
	return status;
}

int get_pure_most(const int *sta) {
	int index = 0, most = 0, all_0 = 1;
	for (index = 0; index < 4; index++) {
		if (sta[index] > 0)
			all_0 = 0;
		if (sta[index] > sta[most])
			most = index;
	}
	if (all_0)
		return -1;
	return most;
}

int *get_most(const int *sta) {
	int *next_c = (int*) calloc(6, sizeof(int));
	unsigned int most = 0, index = 0, total = 0;
	next_c[0] = next_c[1] = next_c[2] = next_c[3] = next_c[4] = next_c[5]
			= INVALID_CHAR;
	total = sta[0] + sta[1] + sta[2] + sta[3] + sta[4];
	float threshold = MULPATH * total;
	threshold = threshold > LEAST_ALIGN_N ? threshold : LEAST_ALIGN_N;
	if (!total)
		return next_c;
	if (sta[0] >= most) {
		most = sta[0];
		next_c[0] = 0;
	}
	if (sta[1] >= most) {
		most = sta[1];
		next_c[0] = 1;
	}
	if (sta[2] >= most) {
		most = sta[2];
		next_c[0] = 2;
	}
	if (sta[3] >= most) {
		most = sta[3];
		next_c[0] = 3;
	}
	if (sta[4] > most) {
		most = sta[4];
	}
	if (most < LEAST_ALIGN * total)
		return next_c;
	if ((sta[0] >= threshold) && next_c[0] != 0)
		next_c[++index] = 0;
	if ((sta[1] >= threshold) && next_c[0] != 1)
		next_c[++index] = 1;
	if ((sta[2] >= threshold) && next_c[0] != 2)
		next_c[++index] = 2;
	if ((sta[3] >= threshold) && next_c[0] != 3)
		next_c[++index] = 3;
	if ((sta[4] >= threshold) && next_c[0] != 4)
		next_c[++index] = 4;
	next_c[++index] = INVALID_CHAR;
	//	printf("[get_most] a:c:g:t:n => %d:%d:%d:%d:%d \n", sta->n_a, sta->n_c, sta->n_g, sta->n_t, sta->n_n);
	return next_c;
}

/**
 * Example:
 * ori =   101110
 * new_c =       11
 * n = 3
 *
 * Return: 00111011
 */
uint64_t shift_bit_to_left(uint64_t kmer, const int new_c, const int n) {
	uint64_t shifted = kmer;
	uint64_t all_ones = 1;
	all_ones <<= n * 2;
	all_ones--;
	shifted <<= 2;
	shifted |= (3 & new_c); // Set the new two lower bits
	shifted &= all_ones; // Set the old two upper bits
	return shifted;
}

uint64_t shift_bit_to_right(uint64_t kmer, const int new_c, const int n) {
	uint64_t shifted = kmer;
	uint64_t upper = new_c;
	shifted >>= 2;
	upper <<= (n - 1) * 2;
	shifted += upper;
	return shifted;
}

uint64_t shift_bit(uint64_t kmer, const int new_c, const int n, const int ori) {
	if (ori)
		return shift_bit_to_right(kmer, new_c, n);
	else
		return shift_bit_to_left(kmer, new_c, n);
}

/**
 In an array, if some value is maximum, return its index
 **/
int get_max_index(const int *counters) {
	int i = 0, max = 0, max_index = 0;
	for (i = 0; i < 4; i++) {
		if (counters[i] > max) {
			max = counters[i];
			max_index = i;
		}
	}
	// If no probable next char, return -1
	if (max == 0)
		return -1;
	return max_index;
}

int get_second_freq_char(const int *counters, const int max) {
	int *counters_without_max = (int*) calloc(4, sizeof(int));
	int i = 0, ret_val = 0;
	for (i = 0; i < 4; i++) {
		if (i == max)
			counters_without_max[i] = 0;
		else
			counters_without_max[i] = counters[i];
	}
	ret_val = get_max_index(counters_without_max);
	free(counters_without_max);
	return ret_val;
}
