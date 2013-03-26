#ifndef LH3_UTILS_H
#define LH3_UTILS_H

#include <stdio.h>
#include <stdint.h>
#include <zlib.h>
#include <inttypes.h>

#define err_fatal_simple(msg) err_fatal_simple_core(__func__, msg)
#define xopen(fn, mode) err_xopen_core(__func__, fn, mode)
#define xreopen(fn, mode, fp) err_xreopen_core(__func__, fn, mode, fp)
#define xzopen(fn, mode) err_xzopen_core(__func__, fn, mode)
#define xassert(cond, msg) if ((cond) == 0) err_fatal_simple_core(__func__, msg)

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
	int abs(int a);
	void trim(char *str);

#ifdef __cplusplus
}
#endif

#endif
