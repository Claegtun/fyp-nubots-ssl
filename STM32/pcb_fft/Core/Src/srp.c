/*
 * srp.c
 *
 *  Created on: Oct 25, 2023
 *      Author: clayton
 */

#include "srp.h"
#include <stdlib.h>

//#define SET_TEST() LED_0_GPIO_Port->BSRR = LED_0_Pin;
//#define RESET_TEST() LED_0_GPIO_Port->BSRR = (uint32_t)LED_0_Pin << 16U;
#define SET_TEST() ;
#define RESET_TEST() ;

void compute_gcc_phat(
		arm_rfft_fast_instance_f32* fft,
		float32_t* R,
		float32_t* X_0,
		float32_t* X_1,
		float32_t* abs_X_0,
		float32_t* abs_X_1,
		float32_t* conj_X_1
) {
	float32_t* divisor = malloc(512*4);
	SET_TEST();
	arm_mult_f32(abs_X_0, abs_X_1, divisor, 512);
	RESET_TEST();

	float32_t* mult_X = malloc(512*4*2);
	SET_TEST();
	arm_cmplx_mult_cmplx_f32(X_0, conj_X_1, mult_X, 512);
	RESET_TEST();

	SET_TEST();
	float32_t* chi = malloc(512*4*2);
	for (int i = 0; i < 2*512; i++) {
		chi[i] = mult_X[i] / divisor[i/2];
	}
	RESET_TEST();

	SET_TEST();
	arm_rfft_fast_f32(fft, chi, R, 1);
	RESET_TEST();

	free(divisor);
	free(mult_X);
	free(chi);
}

float32_t compute_beam_energy(
		arm_rfft_fast_instance_f32* fft,
		uint16_t g,
		float32_t** R,
		const volatile int8_t** tau,
		uint8_t i_R[8][8]
) {
	float32_t E = 0;
	for (int i = 0; i < 8; i++) {
		for (int j = 0; j < 8; j++) {
			if (i == j) continue;
			if (i < j) {
				E += R[
				  i_R[i][j]
				][
					((uint16_t)tau[ i_R[i][j] ][g]) & 0x03FF
				];
			} else {
				E += R[
				  i_R[j][i]
				][
					(((uint16_t)tau[ i_R[j][i] ][g] ^ 0xFFFF) + 0x0001) & 0x03FF
				];
			}
		}
	}

	return E;
}
