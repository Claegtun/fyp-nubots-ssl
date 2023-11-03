/*
 * fir.h
 *
 *  Created on: Nov 3, 2023
 *      Author: clayton
 */

#ifndef INC_FIR_H_
#define INC_FIR_H_

#include "arm_math.h"

void fir_set_up();

void fir_preload();

void fir_filter(int16_t* input_data, uint16_t input_length, int16_t* output_data, uint16_t output_length);

#endif /* INC_FIR_H_ */
