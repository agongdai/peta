#ifndef LH3_UTILS_H
#define LH3_UTILS_H

#include <stdio.h>
#include <stdint.h>
#include <zlib.h>
#include <inttypes.h>
#include <glib.h>

#define err_fatal_simple(msg) err_fatal_simple_core(__func__, msg)
#define xopen(fn, mode) err_xopen_core(__func__, fn, mode)
#define xreopen(fn, mode, fp) err_xreopen_core(__func__, fn, mode, fp)
#define xzopen(fn, mode) err_xzopen_core(__func__, fn, mode)
#define xassert(cond, msg) if ((cond) == 0) err_fatal_simple_core(__func__, msg)

#define average(value1, value2)					( ((value1) & (value2)) + ((value1) ^ (value2)) / 2 )
#define min(value1, value2)						( ((value1) < (value2)) ? (value1) : (value2) )
#define max(value1, value2)						( ((value1) > (value2)) ? (value1) : (value2) )
#define med3(a, b, c)							( a<b ? (b<c ? b : a<c ? c : a) : (b>c ? b : a>c ? c : a))
#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

#define ID64 					PRId64
#define BUFSIZE 				1023
#define FNLEN 					1023
#define LINELEN					50
#define NO_REPEAT_BASES			4
#define NO_REPEAT_LEN 			15
#define	MATCH_SCORE				2
#define	MISMATCH_SCORE			-1
#define	INDEL_SCORE				-2
#define N_MISMATCHES			2
#define MIN_PAIRS				2
#define MIN_JUNCTION_READS		1
#define MIN_WEIGHT				2

typedef uint64_t index64;

#ifdef __cplusplus
extern "C" {
#endif

	void err_fatal(const char *header, const char *fmt, ...);
	void show_msg(const char *header, const char *fmt, ...);
	void show_debug_msg(const char *header, const char *fmt, ...);
	void err_fatal_simple_core(const char *func, const char *msg);
	FILE *err_xopen_core(const char *func, const char *fn, const char *mode);
	FILE *err_xreopen_core(const char *func, const char *fn, const char *mode, FILE *fp);
	gzFile err_xzopen_core(const char *func, const char *fn, const char *mode);
	double std_dev(double a[], double n);
	double mean(double a[], double n);
	float get_abs(float a);
	void trim(char *str);
	char *get_output_file(const char *file_name, const char *out_root);
	int find_in_array(GPtrArray *arr, gpointer value);

#ifdef __cplusplus
}
#endif

#endif
