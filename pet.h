/*
 * pet.h
 *
 *  Created on: 23-Jun-2011
 *      Author: carl
 */
#ifndef PET_H_
#define PET_H_

#include "bwtaln.h"
#include "pealn.h"

#ifdef __cplusplus
extern "C" {
#endif

pool *pet_cvg(const char *pet_fn, const ass_opt *opt);

#ifdef __cplusplus
}
#endif

#endif /* PET_H_ */
