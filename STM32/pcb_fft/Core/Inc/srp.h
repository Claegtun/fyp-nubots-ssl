/*
 * srp.h
 *
 *  Created on: Oct 25, 2023
 *      Author: clayton
 */

#ifndef INC_SRP_H_
#define INC_SRP_H_

#include "arm_math.h"
#include "main.h"

void compute_gcc_phat(
		arm_rfft_fast_instance_f32* fft,
		float32_t* R,
		float32_t* X_0,
		float32_t* X_1,
		float32_t* abs_X_0,
		float32_t* abs_X_1,
		float32_t* conj_X_1
);

float32_t compute_beam_energy(
		arm_rfft_fast_instance_f32* fft,
		uint16_t g,
		float32_t** R,
		const volatile int8_t** tau,
		uint8_t i_R[8][8]
);

#endif /* INC_SRP_H_ */
