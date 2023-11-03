/* USER CODE BEGIN Header */
/**
  ******************************************************************************
  * @file    fmac.c
  * @brief   This file provides code for the configuration
  *          of the FMAC instances.
  ******************************************************************************
  * @attention
  *
  * Copyright (c) 2023 STMicroelectronics.
  * All rights reserved.
  *
  * This software is licensed under terms that can be found in the LICENSE file
  * in the root directory of this software component.
  * If no LICENSE file comes with this software, it is provided AS-IS.
  *
  ******************************************************************************
  */
/* USER CODE END Header */
/* Includes ------------------------------------------------------------------*/
#include "fmac.h"

/* USER CODE BEGIN 0 */
#include "dma.h"

/* USER CODE END 0 */

FMAC_HandleTypeDef hfmac;
DMA_HandleTypeDef hdma_fmac_rd;
DMA_HandleTypeDef hdma_fmac_wr;

/* FMAC init function */
void MX_FMAC_Init(void)
{

  /* USER CODE BEGIN FMAC_Init 0 */

  /* USER CODE END FMAC_Init 0 */

  /* USER CODE BEGIN FMAC_Init 1 */

  /* USER CODE END FMAC_Init 1 */
  hfmac.Instance = FMAC;
  if (HAL_FMAC_Init(&hfmac) != HAL_OK)
  {
    Error_Handler();
  }
  /* USER CODE BEGIN FMAC_Init 2 */

  /* USER CODE END FMAC_Init 2 */

}

void HAL_FMAC_MspInit(FMAC_HandleTypeDef* fmacHandle)
{

  if(fmacHandle->Instance==FMAC)
  {
  /* USER CODE BEGIN FMAC_MspInit 0 */

  /* USER CODE END FMAC_MspInit 0 */
    /* FMAC clock enable */
    __HAL_RCC_FMAC_CLK_ENABLE();

    /* FMAC DMA Init */
    /* FMAC_RD Init */
    hdma_fmac_rd.Instance = DMA1_Stream0;
    hdma_fmac_rd.Init.Request = DMA_REQUEST_FMAC_READ;
    hdma_fmac_rd.Init.Direction = DMA_PERIPH_TO_MEMORY;
    hdma_fmac_rd.Init.PeriphInc = DMA_PINC_DISABLE;
    hdma_fmac_rd.Init.MemInc = DMA_MINC_ENABLE;
    hdma_fmac_rd.Init.PeriphDataAlignment = DMA_PDATAALIGN_HALFWORD;
    hdma_fmac_rd.Init.MemDataAlignment = DMA_MDATAALIGN_HALFWORD;
    hdma_fmac_rd.Init.Mode = DMA_NORMAL;
    hdma_fmac_rd.Init.Priority = DMA_PRIORITY_LOW;
    hdma_fmac_rd.Init.FIFOMode = DMA_FIFOMODE_DISABLE;
    if (HAL_DMA_Init(&hdma_fmac_rd) != HAL_OK)
    {
      Error_Handler();
    }

    __HAL_LINKDMA(fmacHandle,hdmaIn,hdma_fmac_rd);

    /* FMAC_WR Init */
    hdma_fmac_wr.Instance = DMA1_Stream1;
    hdma_fmac_wr.Init.Request = DMA_REQUEST_FMAC_WRITE;
    hdma_fmac_wr.Init.Direction = DMA_MEMORY_TO_PERIPH;
    hdma_fmac_wr.Init.PeriphInc = DMA_PINC_DISABLE;
    hdma_fmac_wr.Init.MemInc = DMA_MINC_ENABLE;
    hdma_fmac_wr.Init.PeriphDataAlignment = DMA_PDATAALIGN_HALFWORD;
    hdma_fmac_wr.Init.MemDataAlignment = DMA_MDATAALIGN_HALFWORD;
    hdma_fmac_wr.Init.Mode = DMA_NORMAL;
    hdma_fmac_wr.Init.Priority = DMA_PRIORITY_LOW;
    hdma_fmac_wr.Init.FIFOMode = DMA_FIFOMODE_DISABLE;
    if (HAL_DMA_Init(&hdma_fmac_wr) != HAL_OK)
    {
      Error_Handler();
    }

    __HAL_LINKDMA(fmacHandle,hdmaOut,hdma_fmac_wr);

  /* USER CODE BEGIN FMAC_MspInit 1 */

	/* Connect the DMA channel to the FMAC handle */
  __HAL_LINKDMA(fmacHandle,hdmaIn,hdma_fmac_wr);
  __HAL_LINKDMA(fmacHandle,hdmaOut,hdma_fmac_rd);
	__HAL_LINKDMA(fmacHandle,hdmaPreload,hdma_memtomem_dma1_stream4);

  /* USER CODE END FMAC_MspInit 1 */
  }
}

void HAL_FMAC_MspDeInit(FMAC_HandleTypeDef* fmacHandle)
{

  if(fmacHandle->Instance==FMAC)
  {
  /* USER CODE BEGIN FMAC_MspDeInit 0 */

  /* USER CODE END FMAC_MspDeInit 0 */
    /* Peripheral clock disable */
    __HAL_RCC_FMAC_CLK_DISABLE();

    /* FMAC DMA DeInit */
    HAL_DMA_DeInit(fmacHandle->hdmaIn);
    HAL_DMA_DeInit(fmacHandle->hdmaOut);
  /* USER CODE BEGIN FMAC_MspDeInit 1 */

  /* USER CODE END FMAC_MspDeInit 1 */
  }
}

/* USER CODE BEGIN 1 */

/* USER CODE END 1 */
