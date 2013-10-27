/*
 * merge.hpp
 *
 *  Created on: Jul 29, 2013
 *      Author: carl
 */

#ifndef MERGE_HPP_
#define MERGE_HPP_

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <glib.h>
#include "utils.h"
#include "tpl.hpp"
#include "peseq.h"

void mv_reads_bt_tpls(tpl *from, tpl *to, int ol, int ori);
int merge_tpls(tpl *left, tpl *right, int ol, int rev_com);
int merged_jumped(hash_table *ht, tpl *t, tpl *jumped, int mis);
void merge_tpl_to_right(tpl *t, tpl *jumped, int ol, int rev_com);
void merge_tpl_to_left(tpl *t, tpl *jumped, int ol, int rev_com);

#endif /* MERGE_HPP_ */
