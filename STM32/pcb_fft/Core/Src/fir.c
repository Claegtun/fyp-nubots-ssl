/*
 * fir.c
 *
 *  Created on: Nov 3, 2023
 *      Author: clayton
 */

#include "fir.h"
#include "fmac.h"
#include "arm_math.h"
#include "coefficients.h"

static volatile uint8_t fir_flags;

#define FIR_CONFIG_DONE 0x01
#define FIR_PRELOAD_DONE 0x02
#define FIR_FILTER_DONE 0x04

void fir_set_up() {

	fir_flags = 0x00;

	q15_t q_coefficients[64];
	arm_float_to_q15(coefficients, q_coefficients, 64);

	/*SET_BIT(
			FMAC->X1BUFCFG,
			64+1 << FMAC_X1BUFCFG_X1_BUF_SIZE_Pos |
			64 << FMAC_X1BUFCFG_X1_BASE_Pos
	);
	SET_BIT(
			FMAC->X2BUFCFG,
			64 << FMAC_X2BUFCFG_X2_BUF_SIZE_Pos |
			0 << FMAC_X2BUFCFG_X2_BASE_Pos
	);
	SET_BIT(
			FMAC->YBUFCFG,
			1 << FMAC_YBUFCFG_Y_BUF_SIZE_Pos |
			64 << FMAC_YBUFCFG_Y_BASE_Pos
	);
	SET_BIT(
			FMAC->CR,
			FMAC_CR_DMAWEN |
			FMAC_CR_DMAREN
	);

	SET_BIT(
			FMAC->PARAM,
			FMAC_PARAM_FUNC_2 |
			64 << FMAC_PARAM_P_Pos |
			0 << FMAC_PARAM_Q_Pos
	);

	HAL_DMA_Start(&hdma_fmac_rd, (uint32_t)q_coefficients, FMAC->WDATA, 64);*/

	FMAC_FilterConfigTypeDef fir_config = {
			.CoeffBSize = 64,
			.CoeffBaseAddress = 0,
			.CoeffASize = 0,
			.CoeffBufferSize = 64,
			.Filter = FMAC_FUNC_CONVO_FIR,
			.InputAccess = FMAC_BUFFER_ACCESS_DMA,
			.InputBaseAddress = 64,
			.InputBufferSize = 64+2,
			.InputThreshold = FMAC_THRESHOLD_1,
			.OutputAccess = FMAC_BUFFER_ACCESS_DMA,
			.OutputBaseAddress = 64,
			.OutputBufferSize = 64+2,
			.OutputThreshold = FMAC_THRESHOLD_1,
			.P = 64+1,
			.R = 0,
			.pCoeffA = NULL,
			.pCoeffB = q_coefficients
	};

	HAL_FMAC_FilterConfig_DMA(&hfmac, &fir_config);

	while(!(fir_flags & FIR_CONFIG_DONE));
}

void fir_preload() {
	q15_t state[64+2];
	for (int i = 0; i < 64+2; i++)
		state[i] = 0;

	HAL_FMAC_FilterPreload_DMA(&hfmac, state, 64+2, NULL, 0);

	while(!(fir_flags & FIR_PRELOAD_DONE));
}

void fir_filter(int16_t* input_data, uint16_t input_length, int16_t* output_data, uint16_t output_length) {
	HAL_FMAC_FilterStart(&hfmac, output_data, &output_length);
	HAL_FMAC_AppendFilterData(&hfmac, input_data, &input_length);

	while(!(fir_flags & FIR_FILTER_DONE));
}

void HAL_FMAC_FilterConfigCallback(FMAC_HandleTypeDef* _hfmac) {
	if (_hfmac == &hfmac) {
		fir_flags |= FIR_CONFIG_DONE;
	}
}

void HAL_FMAC_FilterPreloadCallback(FMAC_HandleTypeDef *_hfmac) {
	if (_hfmac == &hfmac) {
		fir_flags |= FIR_PRELOAD_DONE;
	}
}

void HAL_FMAC_OutputDataReadyCallback(FMAC_HandleTypeDef *_hfmac) {
	if (_hfmac == &hfmac) {
		fir_flags |= FIR_FILTER_DONE;
	}
}
