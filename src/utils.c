#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <math.h>
#include <glib.h>
#include "utils.h"

FILE *err_xopen_core(const char *func, const char *fn, const char *mode) {
	FILE *fp = 0;
	if (strcmp(fn, "-") == 0)
		return (strstr(mode, "r")) ? stdin : stdout;
	if ((fp = fopen(fn, mode)) == 0) {
		fprintf(stderr, "[%s] fail to open file '%s'. Abort!\n", func, fn);
		abort();
	}
	return fp;
}
FILE *err_xreopen_core(const char *func, const char *fn, const char *mode,
		FILE *fp) {
	if (freopen(fn, mode, fp) == 0) {
		fprintf(stderr, "[%s] fail to open file '%s': ", func, fn);
		perror(NULL);
		fprintf(stderr, "Abort!\n");
		abort();
	}
	return fp;
}
gzFile err_xzopen_core(const char *func, const char *fn, const char *mode) {
	gzFile fp;
	if (strcmp(fn, "-") == 0)
		return gzdopen(fileno((strstr(mode, "r")) ? stdin : stdout), mode);
	if ((fp = gzopen(fn, mode)) == 0) {
		fprintf(stderr, "[%s] fail to open file '%s'. Abort!\n", func, fn);
		abort();
	}
	return fp;
}
void err_fatal(const char *header, const char *fmt, ...) {
	va_list args;
	va_start(args, fmt);
	fprintf(stderr, "[%s] ", header);
	vfprintf(stderr, fmt, args);
	fprintf(stderr, " Abort!\n");
	va_end(args);
	abort();
}
void show_msg(const char *header, const char *fmt, ...) {
	va_list args;
	va_start(args, fmt);
	fprintf(stderr, "[%s] ", header);
	vfprintf(stderr, fmt, args);
	va_end(args);
}
void show_debug_msg(const char *header, const char *fmt, ...) {
	va_list args;
	va_start(args, fmt);
	printf("[%s] ", header);
	vprintf(fmt, args);
	va_end(args);
}

void err_fatal_simple_core(const char *func, const char *msg) {
	fprintf(stderr, "[%s] %s Abort!\n", func, msg);
	abort();
}

double mean(double a[], double n) {
	int i = 0;
	double sum = 0;
	for (i = 0; i < n; i++) {
		//printf("%f \n", a[i]);
		sum += a[i];
	}
	return sum / n;
}

double std_dev(double a[], double n) {
	if (n == 0)
		return 0.0;
	double sum = 0;
	double sq_sum = 0;
	int i = 0;
	for (i = 0; i < n; ++i) {
		sum += a[i];
		sq_sum += a[i] * a[i];
	}
	double mean = sum / n;
	double variance = sq_sum / n - mean * mean;
	return sqrt(variance);
}

float get_abs(float a) {
	return (a > 0.0) ? a : (0.0 - a);
}

void trim(char *str) {
	int i = 0;
	for (i = 0; i < strlen(str); i++) {
		if (str[i] == '\n') {
			str[i] = '\0';
			return;
		}
	}
}

char *get_output_file(const char *file_name, const char *out_root) {
	char *name = (char*) calloc(strlen(file_name) + strlen(out_root) + 4,
			sizeof(char));
	strcat(name, out_root);
	strcat(name, "/");
	strcat(name, file_name);
	return name;
}

int find_in_array(GPtrArray *arr, gpointer value) {
	uint32_t i = 0;
	gpointer v = NULL;
	for (i = 0; i < arr->len; i++) {
		v = g_ptr_array_index(arr, i);
		if (v == value)
			return i;
	}
	return -1;
}

int min3(int a, int b, int c) {
	if (a <= b && a <= c)
		return a;
	if (b <= a && b <= c)
		return b;
	if (c <= a && c <= b)
		return c;
}

char *str_dup(const char * s) {
  size_t len = 1+strlen(s);
  char *p = malloc(len);
  return p ? memcpy(p, s, len) : NULL;
}
