/* mystack.h -- Stack declaration and function prototypes:  */
#ifndef LIST_H_
#define LIST_H_

#include <glib.h>

typedef struct
{
	GPtrArray *values;
	GCompareFunc cmp_f; // Function for comparison
	int order; 			// 0: increasing; 1: reversed; -1: not sorted
} slist;

slist *new_slist(GCompareFunc cmp_f);
int slist_binary(slist *sl, gpointer query);
int slist_ins_pt(slist *sl, gpointer new_pt);
gint cmp_kmer_by_seq(gpointer a, gpointer b);

#endif
