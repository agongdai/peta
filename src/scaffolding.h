/*
 * scaffolding.h
 *
 *  Created on: 04-Feb-2013
 *      Author: carl
 */

#ifndef SCAFFOLDING_H_
#define SCAFFOLDING_H_

#include "edgelist.h"
#include <glib.h>

#ifdef __cplusplus
extern "C" {
#endif

void scaffolding(edgearray *all_edges, const int insert_size, const int sd_insert_size, hash_table *ht, const int n_threads);

#ifdef __cplusplus
}
#endif
#endif /* SCAFFOLDING_H_ */
