/* USER CODE BEGIN Header */
/**
  ******************************************************************************
  * @file           : main.c
  * @brief          : Main program body
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
#include "main.h"
#include "dma.h"
#include "tim.h"
#include "usb_device.h"
#include "gpio.h"

/* Private includes ----------------------------------------------------------*/
/* USER CODE BEGIN Includes */
#include "arm_math.h"
#include "coefficients.h"
#include "usbd_cdc_if.h"
#include "string.h"
#include "srp.h"
#include "tables.h"

/* USER CODE END Includes */

/* Private typedef -----------------------------------------------------------*/
/* USER CODE BEGIN PTD */

/* USER CODE END PTD */

/* Private define ------------------------------------------------------------*/
/* USER CODE BEGIN PD */

/* USER CODE END PD */

/* Private macro -------------------------------------------------------------*/
/* USER CODE BEGIN PM */
//#define SET_TEST() LED_0_GPIO_Port->BSRR = LED_0_Pin
//#define RESET_TEST() LED_0_GPIO_Port->BSRR = (uint32_t)LED_0_Pin << 16U
#define SET_TEST() ;
#define RESET_TEST() ;

/* USER CODE END PM */

/* Private variables ---------------------------------------------------------*/

/* USER CODE BEGIN PV */
volatile static uint16_t data_buffer[4][1024];
volatile static uint32_t rx_buffer[1024*8];

volatile static uint8_t adc_parallel_dma_flag;

/* USER CODE END PV */

/* Private function prototypes -----------------------------------------------*/
void SystemClock_Config(void);
/* USER CODE BEGIN PFP */
void DMA_TIM8_callback(DMA_HandleTypeDef* hdma);

/* USER CODE END PFP */

/* Private user code ---------------------------------------------------------*/
/* USER CODE BEGIN 0 */

/* USER CODE END 0 */

/**
  * @brief  The application entry point.
  * @retval int
  */
int main(void)
{
  /* USER CODE BEGIN 1 */

  /* USER CODE END 1 */

  /* MCU Configuration--------------------------------------------------------*/

  /* Reset of all peripherals, Initializes the Flash interface and the Systick. */
  HAL_Init();

  /* USER CODE BEGIN Init */

  /* USER CODE END Init */

  /* Configure the system clock */
  SystemClock_Config();

  /* USER CODE BEGIN SysInit */

  /* USER CODE END SysInit */

  /* Initialize all configured peripherals */
  MX_GPIO_Init();
  MX_DMA_Init();
  MX_TIM3_Init();
  MX_TIM4_Init();
  MX_TIM5_Init();
  MX_TIM8_Init();
  MX_USB_DEVICE_Init();
  MX_TIM2_Init();
  /* USER CODE BEGIN 2 */

  // Allocate the memory for the DMA buffers.
	// + 256 KB
	volatile uint16_t* buffer_0 = malloc(8*8*1024*2);
	volatile uint16_t* buffer_1 = malloc(8*8*1024*2);

  // Make the RFFT instance.
  arm_rfft_fast_instance_f32 fft;
  arm_rfft_fast_init_f32(&fft, 1024);

  // Make the FIR instance.
	q15_t state[8320];
	for (int i = 0; i < 8320; i++)
		state[i] = 0;

	q15_t q_coefficients[64];
	arm_float_to_q15(coefficients, q_coefficients, 64);

	arm_fir_decimate_instance_q15 fir;
	arm_fir_decimate_init_q15(&fir, 64, 8, q_coefficients, state, 1024);

	// Allocate the memory for the deinterleaved data.
	// + 128 KB
	uint16_t* oversampled[8];
	for (int16_t i = 0; i < 8; i++)
		oversampled[i] = (uint16_t*)(D1_DTCMRAM_BASE + 8*1024*2*i);

	uint8_t i_R[8][8];
	uint8_t count = 0;
	for (int i = 0; i < 8; i++) {
		for (int j = i+1;  j < 8; j++) {
			i_R[i][j] = count++;
		}
	}

	/*// Allocate the memory for the TDOAs.
	// + 71 KB
	uint8_t* tau[28];
	for (int16_t i = 0; i < 28; i++)
		tau[i] = (uint8_t*)(D1_DTCMRAM_BASE + 0x00000000 + 2562*i);*/

	const volatile int8_t* tau[28];
	tau[0] = tau_table_0;
	tau[1] = tau_table_1;
	tau[2] = tau_table_2;
	tau[3] = tau_table_3;
	tau[4] = tau_table_4;
	tau[5] = tau_table_5;
	tau[6] = tau_table_6;
	tau[7] = tau_table_7;
	tau[8] = tau_table_8;
	tau[9] = tau_table_9;
	tau[10] = tau_table_10;
	tau[11] = tau_table_11;
	tau[12] = tau_table_12;
	tau[13] = tau_table_13;
	tau[14] = tau_table_14;
	tau[15] = tau_table_15;
	tau[16] = tau_table_16;
	tau[17] = tau_table_17;
	tau[18] = tau_table_18;
	tau[19] = tau_table_19;
	tau[20] = tau_table_20;
	tau[21] = tau_table_21;
	tau[22] = tau_table_22;
	tau[23] = tau_table_23;
	tau[24] = tau_table_24;
	tau[25] = tau_table_25;
	tau[26] = tau_table_26;
	tau[27] = tau_table_27;

	adc_parallel_dma_flag = 0;

  CLEAR_BIT(TIM8->CR1, TIM_CR1_ARPE);
	TIM8->ARR = 53; //27
	TIM8->CCR1 = 40; //14
	TIM8->RCR = 7;
	SET_BIT(TIM8->CR1, TIM_CR1_ARPE);

	//SET_BIT(TIM8->CCER, TIM_CCER_CC1P);

	// Reset the ADC.
	HAL_GPIO_WritePin(ADC_RESET_GPIO_Port, ADC_RESET_Pin, GPIO_PIN_SET);
	HAL_Delay(1); // > 3 us
	HAL_GPIO_WritePin(ADC_RESET_GPIO_Port, ADC_RESET_Pin, GPIO_PIN_RESET);
	HAL_Delay(1); // > 253 us

	RCC->AHB4ENR |= RCC_AHB4ENR_GPIODEN;
	WRITE_REG(GPIOD->MODER, 0x00000000UL);

	__HAL_TIM_ENABLE_DMA(&htim5, TIM_DMA_CC1);
	htim5.hdma[TIM_DMA_ID_CC1]->XferCpltCallback = DMA_TIM8_callback;
	HAL_DMA_Start_IT(htim5.hdma[TIM_DMA_ID_CC1], (uint32_t)&GPIOD->IDR + 0, (uint32_t)buffer_0, 1024*8*8-1);

	// Start the SCLK timer, the nCS timer, and the CONVST timer respectively.
	HAL_TIM_PWM_Start(&htim5, TIM_CHANNEL_1);
	HAL_TIM_PWM_Start(&htim8, TIM_CHANNEL_1);
	HAL_TIM_PWM_Start(&htim3, TIM_CHANNEL_1);
	HAL_TIM_PWM_Start(&htim4, TIM_CHANNEL_1);
	HAL_TIM_PWM_Start(&htim2, TIM_CHANNEL_1);

  /* USER CODE END 2 */

  /* Infinite loop */
  /* USER CODE BEGIN WHILE */
  while (1)
  {
  	//HAL_GPIO_TogglePin(LED_0_GPIO_Port, LED_0_Pin);
  	//HAL_Delay(500);

  	if (adc_parallel_dma_flag == 1) {
  		// Deinterleave the raw data.
  		SET_TEST();
			/*for (int32_t i = 0; i < 8*8*1024; i++) {
				oversampled[i%8][i/((int32_t)8)] = buffer_0[i];
			}*/
  		for (int32_t i = 0; i < 8*1024; i ++) {
  			oversampled[0][i] = buffer_0[i*8 + 0];
  		}
  		for (int32_t i = 0; i < 8*1024; i ++) {
				oversampled[1][i] = buffer_0[i*8 + 1];
			}
  		for (int32_t i = 0; i < 8*1024; i ++) {
				oversampled[2][i] = buffer_0[i*8 + 2];
			}
  		for (int32_t i = 0; i < 8*1024; i ++) {
				oversampled[3][i] = buffer_0[i*8 + 3];
			}
  		for (int32_t i = 0; i < 8*1024; i ++) {
				oversampled[4][i] = buffer_0[i*8 + 4];
			}
  		for (int32_t i = 0; i < 8*1024; i ++) {
				oversampled[5][i] = buffer_0[i*8 + 5];
			}
  		for (int32_t i = 0; i < 8*1024; i ++) {
				oversampled[6][i] = buffer_0[i*8 + 6];
			}
  		for (int32_t i = 0; i < 8*1024; i ++) {
				oversampled[7][i] = buffer_0[i*8 + 7];
			}
			RESET_TEST();

			// Allocate the memory for the decimated data.
			// + 16 KB
			q15_t* decimated[8];
			for (int16_t i = 0; i < 8; i++)
				decimated[i] = (q15_t*)(D3_SRAM_BASE + 1024*2*i);

			// Decimate the deinterleaved oversampled data.
			SET_TEST();
			for (int8_t i = 0; i < 8; i ++)
				arm_fir_decimate_q15(&fir, (q15_t*)oversampled[i], decimated[i], 8192);
			RESET_TEST();

			// Allocate the memory for the time-domain data.
			// + 32 KB
			float32_t* x[8];
			for (int16_t i = 0; i < 8; i++)
				x[i] = (float32_t*)(D1_DTCMRAM_BASE + 1024*4*i);

			// Convert the data from fixed-point to floating-point.
			SET_TEST();
			for (int8_t i = 0; i < 8; i ++)
				arm_q15_to_float(decimated[i], x[i], 1024);
			RESET_TEST();

			char header[3] = "\r\nH\r\n";
			char str[64] = "";
			/*for (int j = 0; j < 8; j++) {

			}*/

			// Allocate the memory for the frequency-domain data.
			// + 32 KB
			float32_t* X[8];
			for (int16_t i = 0; i < 8; i++)
				X[i] = (float32_t*)(D1_ITCMRAM_BASE + 0x00000000 + 512*4*2*i);

			// Calculate the FFT of each channel.
			SET_TEST();
			for (int i = 0; i < 8; i++)
				arm_rfft_fast_f32(&fft, x[i], X[i], 0);
			RESET_TEST();

			// Allocate the memory for the complex magnitudes.
			// + 16 KB
			float32_t* abs_X[8];
			for (int16_t i = 0; i < 8; i++)
				abs_X[i] = (float32_t*)(D1_DTCMRAM_BASE + 0x0001C000 + 512*4*i);

			// Compute the complex magnitude of the frequency-domain data.
			SET_TEST();
			for (int i = 0; i < 8; i++) {
				arm_cmplx_mag_f32(X[i], abs_X[i], 512);
			}
			RESET_TEST();

			CDC_Transmit_HS((uint8_t*)header, strlen(header));
			for (int i = 0; i < 512; i++) {
				sprintf(str, "%f,", abs_X[0][i]);
				CDC_Transmit_HS((uint8_t*)str, strlen(str));
				HAL_Delay(1);
			}

			// Allocate the memory for the complex conjugates.
			// + 32 KB
			float32_t* conj_X[8];
			for (int16_t i = 0; i < 8; i++)
				conj_X[i] = (float32_t*)(D1_ITCMRAM_BASE + 0x00008000 + 512*4*2*i);

			// Compute the complex conjugate of the frequency-domain data.
			SET_TEST();
			for (int i = 0; i < 8; i++) {
				arm_cmplx_conj_f32(X[i], conj_X[i], 512);
			}
			RESET_TEST();

			// Allocate the memory for the GCCs.
			// + 112 KB
			float32_t* R[28];
			for (int16_t i = 0; i < 28; i++)
				R[i] = (float32_t*)(D1_DTCMRAM_BASE + 0x00000000 + 1024*4*i);

			// Compute the GCCs for all pairs.
			for (int i = 0; i < 8; i++) {
				for (int j = i+1;  j < 8; j++) {
					SET_TEST();
					compute_gcc_phat(&fft, R[i_R[i][j]], X[i], X[j], abs_X[i], abs_X[j], conj_X[j]);
					RESET_TEST();
				}
			}

			CDC_Transmit_HS((uint8_t*)header, strlen(header));

			// Draw the energy-map.
			float32_t E, E_max = 0;
			uint16_t g_max = 0;
			SET_TEST();
			for (int g = 0; g < 2562; g++) {
				E = compute_beam_energy(&fft, g, R, tau, i_R);
				sprintf(str, "%f,", E);
				CDC_Transmit_HS((uint8_t*)str, strlen(str));
				HAL_Delay(1);
				if (E >= E_max) {
					E_max = E;
					g_max = g;
				}
			}
			RESET_TEST();

			/*float32_t theta, phi;
			theta = grid[0][g_max]*180.0/3.1415/2.0;
			phi = grid[1][g_max]*180.0/3.1415/2.0;

			float32_t robert_frost;
			robert_frost = theta;
			robert_frost = phi;*/

			adc_parallel_dma_flag = 1;
  	}

    /* USER CODE END WHILE */

    /* USER CODE BEGIN 3 */
  }
  /* USER CODE END 3 */
}

/**
  * @brief System Clock Configuration
  * @retval None
  */
void SystemClock_Config(void)
{
  RCC_OscInitTypeDef RCC_OscInitStruct = {0};
  RCC_ClkInitTypeDef RCC_ClkInitStruct = {0};
  RCC_CRSInitTypeDef RCC_CRSInitStruct = {0};

  /** Supply configuration update enable
  */
  HAL_PWREx_ConfigSupply(PWR_LDO_SUPPLY);

  /** Configure the main internal regulator output voltage
  */
  __HAL_PWR_VOLTAGESCALING_CONFIG(PWR_REGULATOR_VOLTAGE_SCALE0);

  while(!__HAL_PWR_GET_FLAG(PWR_FLAG_VOSRDY)) {}

  /** Initializes the RCC Oscillators according to the specified parameters
  * in the RCC_OscInitTypeDef structure.
  */
  RCC_OscInitStruct.OscillatorType = RCC_OSCILLATORTYPE_HSI48|RCC_OSCILLATORTYPE_HSI;
  RCC_OscInitStruct.HSIState = RCC_HSI_DIV1;
  RCC_OscInitStruct.HSICalibrationValue = 64;
  RCC_OscInitStruct.HSI48State = RCC_HSI48_ON;
  RCC_OscInitStruct.PLL.PLLState = RCC_PLL_ON;
  RCC_OscInitStruct.PLL.PLLSource = RCC_PLLSOURCE_HSI;
  RCC_OscInitStruct.PLL.PLLM = 4;
  RCC_OscInitStruct.PLL.PLLN = 34;
  RCC_OscInitStruct.PLL.PLLP = 1;
  RCC_OscInitStruct.PLL.PLLQ = 2;
  RCC_OscInitStruct.PLL.PLLR = 2;
  RCC_OscInitStruct.PLL.PLLRGE = RCC_PLL1VCIRANGE_3;
  RCC_OscInitStruct.PLL.PLLVCOSEL = RCC_PLL1VCOWIDE;
  RCC_OscInitStruct.PLL.PLLFRACN = 3072;
  if (HAL_RCC_OscConfig(&RCC_OscInitStruct) != HAL_OK)
  {
    Error_Handler();
  }

  /** Initializes the CPU, AHB and APB buses clocks
  */
  RCC_ClkInitStruct.ClockType = RCC_CLOCKTYPE_HCLK|RCC_CLOCKTYPE_SYSCLK
                              |RCC_CLOCKTYPE_PCLK1|RCC_CLOCKTYPE_PCLK2
                              |RCC_CLOCKTYPE_D3PCLK1|RCC_CLOCKTYPE_D1PCLK1;
  RCC_ClkInitStruct.SYSCLKSource = RCC_SYSCLKSOURCE_PLLCLK;
  RCC_ClkInitStruct.SYSCLKDivider = RCC_SYSCLK_DIV1;
  RCC_ClkInitStruct.AHBCLKDivider = RCC_HCLK_DIV2;
  RCC_ClkInitStruct.APB3CLKDivider = RCC_APB3_DIV2;
  RCC_ClkInitStruct.APB1CLKDivider = RCC_APB1_DIV2;
  RCC_ClkInitStruct.APB2CLKDivider = RCC_APB2_DIV2;
  RCC_ClkInitStruct.APB4CLKDivider = RCC_APB4_DIV2;

  if (HAL_RCC_ClockConfig(&RCC_ClkInitStruct, FLASH_LATENCY_3) != HAL_OK)
  {
    Error_Handler();
  }

  /** Enable the SYSCFG APB clock
  */
  __HAL_RCC_CRS_CLK_ENABLE();

  /** Configures CRS
  */
  RCC_CRSInitStruct.Prescaler = RCC_CRS_SYNC_DIV1;
  RCC_CRSInitStruct.Source = RCC_CRS_SYNC_SOURCE_USB2;
  RCC_CRSInitStruct.Polarity = RCC_CRS_SYNC_POLARITY_RISING;
  RCC_CRSInitStruct.ReloadValue = __HAL_RCC_CRS_RELOADVALUE_CALCULATE(48000000,1000);
  RCC_CRSInitStruct.ErrorLimitValue = 34;
  RCC_CRSInitStruct.HSI48CalibrationValue = 32;

  HAL_RCCEx_CRSConfig(&RCC_CRSInitStruct);
}

/* USER CODE BEGIN 4 */

void DMA_TIM8_callback(DMA_HandleTypeDef* hdma) {
	if (hdma == htim5.hdma[TIM_DMA_ID_CC1]) {
		adc_parallel_dma_flag = 1;
		//HAL_GPIO_TogglePin(SIG_0_GPIO_Port, SIG_0_Pin);
	}
}

/* USER CODE END 4 */

/**
  * @brief  This function is executed in case of error occurrence.
  * @retval None
  */
void Error_Handler(void)
{
  /* USER CODE BEGIN Error_Handler_Debug */
  /* User can add his own implementation to report the HAL error return state */
  __disable_irq();
  while (1)
  {
  }
  /* USER CODE END Error_Handler_Debug */
}

#ifdef  USE_FULL_ASSERT
/**
  * @brief  Reports the name of the source file and the source line number
  *         where the assert_param error has occurred.
  * @param  file: pointer to the source file name
  * @param  line: assert_param error line source number
  * @retval None
  */
void assert_failed(uint8_t *file, uint32_t line)
{
  /* USER CODE BEGIN 6 */
  /* User can add his own implementation to report the file name and line number,
     ex: printf("Wrong parameters value: file %s on line %d\r\n", file, line) */
  /* USER CODE END 6 */
}
#endif /* USE_FULL_ASSERT */
