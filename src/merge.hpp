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

int merge_tpls(tpl *left, tpl *right, int ol, int rev_com);

#endif /* MERGE_HPP_ */
