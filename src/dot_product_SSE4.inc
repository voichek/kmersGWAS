///
///      @file  dot_product_SSE4.cpp
///     @brief  Very efficent dot product between bit-array (uint64_t) to floats array
//		using the specific archtecture of the proccessor (SSE4). 
//		Taken from: https://stackoverflow.com/questions/16051365/fast-dot-product-of-a-bit-vector-and-a-floating-point-vector
//		There is another version using the AVX2 archtecture which is 2 times faster.
//		But part of our cluster in MPI can't run AVX2 so we will stick to generality for now.
///
/// Detailed description starts here.
///
///    @author  Yoav Voichek (YV), yoav.voichek@tuebingen.mpg.de
///
///  @internal
///    Created  11/12/18
///   Revision  $Id: doxygen.cpp.templates,v 1.3 2010/07/06 09:20:12 mehner Exp $
///   Compiler  gcc/g++
///    Company  Max Planck Institute for Developmental Biology Dep 6
///  Copyright  Copyright (c) 2018, Yoav Voichek
///
/// This source code is released for free distribution under the terms of the
/// GNU General Public License as published by the Free Software Foundation.
///=====================================================================================
///


/* Includes */
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <smmintrin.h>  /* SSE 4.1 */
#include <time.h>

/* Defines */
#define ALIGNTO(n) __attribute__((aligned(n)))
#define USE_PINSRW 1



/*
 * Bit mask shuffle.
 * 
 * This version uses a loop to store eight u16 and reloads them as one __m128i.
 */

__m128 bitMaskShuffleStoreAndReload(__m128i mask){
	const __m128i perm    ALIGNTO(16) = _mm_set_epi8(15, 7, 14, 6, 13, 5, 12, 4,
			11, 3, 10, 2, 9,  1, 8,  0);
	int i;
	uint16_t interMask[8] ALIGNTO(16);

	/* Shuffle bitmask */
	/* Stage 1 */
	for(i=7;i>=0;i--){
		interMask[i] = _mm_movemask_epi8(mask);
		mask = _mm_slli_epi32(mask, 1);
	}

	/* Stage 2 */
	return _mm_castsi128_ps(
			_mm_shuffle_epi8(
				_mm_load_si128((const __m128i*)interMask),
				perm)
			);
}

/*
 * Bit mask shuffle.
 * 
 * This version uses the PINSTRW instruction.
 */

__m128 bitMaskShufflePINSRW(__m128i mask){
	const __m128i perm    ALIGNTO(16) = _mm_set_epi8(15, 7, 14, 6, 13, 5, 12, 4,
			11, 3, 10, 2, 9,  1, 8,  0);
	__m128i  imask ALIGNTO(16);

	/* Shuffle bitmask */
	/* Stage 1 */
	imask = _mm_insert_epi16(imask, _mm_movemask_epi8(mask), 7);
	mask = _mm_slli_epi16(mask, 1);
	imask = _mm_insert_epi16(imask, _mm_movemask_epi8(mask), 6);
	mask = _mm_slli_epi16(mask, 1);
	imask = _mm_insert_epi16(imask, _mm_movemask_epi8(mask), 5);
	mask = _mm_slli_epi16(mask, 1);
	imask = _mm_insert_epi16(imask, _mm_movemask_epi8(mask), 4);
	mask = _mm_slli_epi16(mask, 1);
	imask = _mm_insert_epi16(imask, _mm_movemask_epi8(mask), 3);
	mask = _mm_slli_epi16(mask, 1);
	imask = _mm_insert_epi16(imask, _mm_movemask_epi8(mask), 2);
	mask = _mm_slli_epi16(mask, 1);
	imask = _mm_insert_epi16(imask, _mm_movemask_epi8(mask), 1);
	mask = _mm_slli_epi16(mask, 1);
	imask = _mm_insert_epi16(imask, _mm_movemask_epi8(mask), 0);

	/* Stage 2 */
	return _mm_castsi128_ps(
			_mm_shuffle_epi8(
				imask,
				perm)
			);
}

/**
 *  * SSE 4.1 implementation.
 *   */

float dotSSE41(__m128 f[32], unsigned char maskArg[16]){
	int i, j, k;
	__m128i  mask      ALIGNTO(16) = _mm_load_si128((const __m128i*)maskArg);
	__m128   shufdMask ALIGNTO(16);
	__m128   zblended  ALIGNTO(16);
	__m128   sums      ALIGNTO(16) = _mm_setzero_ps();
	float    sumsf[4]  ALIGNTO(16);

	/* Shuffle bitmask */
#if USE_PINSRW
	shufdMask = bitMaskShufflePINSRW(mask);
#else
	shufdMask = bitMaskShuffleStoreAndReload(mask);
#endif

	/* Dot product */
	for(i=1;i>=0;i--){
		for(j=1;j>=0;j--){
			for(k=7;k>=0;k--){
				zblended  = _mm_setzero_ps();
				zblended  = _mm_blendv_ps(zblended, f[i*16+j+k*2], shufdMask);
				sums      = _mm_add_ps(sums, zblended);
				shufdMask = _mm_castsi128_ps(_mm_slli_epi32(_mm_castps_si128(shufdMask), 1));
			}
		}
	}

	/* Final Summation */
	_mm_store_ps(sumsf, sums);
	return sumsf[0] + sumsf[1] + sumsf[2] + sumsf[3];
}

/**
 *  * Reference C implementation
 *   */


float dotRefC(float f[128], unsigned char mask[16]){
	float sum = 0.0;
	int i;

	for(i=0;i<128;i++){
		sum += ((mask[i>>3]>>(i&7))&1) ? f[i] : 0.0;
	}

	return sum;
}
