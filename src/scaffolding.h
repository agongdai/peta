/*
 * scaffolding.h
 *
 *  Created on: 18-Nov-2012
 *      Author: carl
 */

#ifndef SCAFFOLDING_H_
#define SCAFFOLDING_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include "roadmap.h"
#include "edgelist.h"
#include "bwase.h"
#include <glib.h>

GPtrArray *scaffolding(GPtrArray *single_edges, const int insert_size);

#ifdef __cplusplus
}
#endif

#endif /* SCAFFOLDING_H_ */
