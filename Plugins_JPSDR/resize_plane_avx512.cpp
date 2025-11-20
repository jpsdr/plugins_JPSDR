// Avisynth v2.5.  Copyright 2002 Ben Rudiak-Gould et al.
// http://www.avisynth.org

// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA, or visit
// http://www.gnu.org/copyleft/gpl.html .
//
// Linking Avisynth statically or dynamically with other modules is making a
// combined work based on Avisynth.  Thus, the terms and conditions of the GNU
// General Public License cover the whole combination.
//
// As a special exception, the copyright holders of Avisynth give you
// permission to link Avisynth with independent modules that communicate with
// Avisynth solely through the interfaces defined in avisynth.h, regardless of the license
// terms of these independent modules, and to copy and distribute the
// resulting combined work under terms of your choice, provided that
// every copy of the combined work is accompanied by a complete copy of
// the source code of Avisynth (the version of Avisynth used to produce the
// combined work), being distributed under the terms of the GNU General
// Public License plus this exception.  An independent module is a module
// which is not derived from or based on Avisynth, such as 3rd-party filters,
// import and export plugins, or graphical user interfaces.

// VS 2017 v15.3
#if _MSC_VER >= 1911

#include <immintrin.h>
#include "avs/minmax.h"
#include "JincResizeMT.h"

template <typename T, bool bFP16>
#if defined(CLANG)
__attribute__((__target__("avx512f")))
#endif
void resize_plane_avx512_1x(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *tab_coeff,
	const float Val_Min[], const float Val_Max[])
{
	const uint8_t idx = 0;

	const T *src = reinterpret_cast<const T*>(MT_DataGF->src[idx]);
	T *JincMT_RESTRICT dst = reinterpret_cast<T*>(MT_DataGF->dst[idx]);

	const ptrdiff_t src_pitch = (ptrdiff_t)MT_DataGF->src_pitch[idx] / sizeof(T);
	const ptrdiff_t dst_pitch = (ptrdiff_t)MT_DataGF->dst_pitch[idx] / sizeof(T);
	const ptrdiff_t src_pitchx2 = src_pitch << 1;

	const int Y_Min = MT_DataGF->dst_Y_h_min;
	const int Y_Max = MT_DataGF->dst_Y_h_max;
	const int dst_width = MT_DataGF->dst_Y_w;

	const T val_min = static_cast<T>(Val_Min[idx]);
	const T val_max = static_cast<T>(Val_Max[idx]);
	const __m512 min_val = _mm512_set1_ps(Val_Min[idx]);

	EWAPixelCoeffMeta *meta_y = tab_coeff->meta + (Y_Min*dst_width);

	const int filter_size = tab_coeff->filter_size, coeff_stride = tab_coeff->coeff_stride;
	const int coeff_stridex2 = coeff_stride << 1;

	const int filter_size_mod2 = (filter_size >> 1) << 1;
	const bool fs_notMod2 = filter_size_mod2 < filter_size;

    for (int y = Y_Min; y < Y_Max; y++)
    {
		EWAPixelCoeffMeta *meta = meta_y;

		for (int x = 0; x < dst_width; ++x)
        {
            const T *src_ptr = src + (meta->start_y * src_pitch + meta->start_x);
            const float *coeff_ptr = tab_coeff->factor + meta->coeff_meta;
            __m512 result = _mm512_setzero_ps();
			__m512 result2 = _mm512_setzero_ps();

            if JincMT_CONSTEXPR (std::is_same<T, uint8_t>::value)
            {
                for (int ly = 0; ly < filter_size_mod2; ly += 2)
                {
                    for (int lx = 0; lx < filter_size; lx += 16)
                    {
                        const __m512 src_ps = _mm512_cvtepi32_ps(_mm512_cvtepu8_epi32(_mm_loadu_si128(reinterpret_cast<const __m128i*>(src_ptr + lx))));
						const __m512 src_ps2 = _mm512_cvtepi32_ps(_mm512_cvtepu8_epi32(_mm_loadu_si128(reinterpret_cast<const __m128i*>(src_ptr + lx + src_pitch))));
						__m512 coeff, coeff2;

						if JincMT_CONSTEXPR (bFP16)
						{
							const int lx2 = lx >> 1;
							const __m256i coeff_fp16 = _mm256_load_si256((__m256i*)(coeff_ptr + lx2));
							const __m256i coeff2_fp16 = _mm256_load_si256((__m256i*)(coeff_ptr + lx2 + coeff_stride));

							coeff = _mm512_cvtph_ps(coeff_fp16); // separated grouped converts may be better for instructions grouping at some compilers
							coeff2 = _mm512_cvtph_ps(coeff2_fp16);
						}
						else
						{
							coeff = _mm512_load_ps(coeff_ptr + lx);
							coeff2 = _mm512_load_ps(coeff_ptr + lx + coeff_stride);
						}

                        result = _mm512_fmadd_ps(src_ps, coeff, result);
						result2 = _mm512_fmadd_ps(src_ps2, coeff2, result2);
					}

					coeff_ptr += coeff_stridex2;
                    src_ptr += src_pitchx2;
                }

				result = _mm512_add_ps(result, result2);

				// last row
				if (fs_notMod2)
				{
					for (int lx = 0; lx < filter_size; lx += 16)
					{
						const __m512 src_ps = _mm512_cvtepi32_ps(_mm512_cvtepu8_epi32(_mm_loadu_si128(reinterpret_cast<const __m128i*>(src_ptr + lx))));
						__m512 coeff;

						if JincMT_CONSTEXPR (bFP16)
							coeff = _mm512_cvtph_ps(_mm256_load_si256((__m256i*)(coeff_ptr + (lx >> 1))));
						else
							coeff = _mm512_load_ps(coeff_ptr + lx);

						result = _mm512_fmadd_ps(src_ps, coeff, result);
					}
				}

                const __m256 lo_hi_256 = _mm256_add_ps(_mm512_castps512_ps256(result), _mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(result), 1)));
                __m128 hsum = _mm_add_ps(_mm256_castps256_ps128(lo_hi_256), _mm256_extractf128_ps(lo_hi_256, 1));
                hsum = _mm_hadd_ps(_mm_hadd_ps(hsum, hsum), _mm_hadd_ps(hsum, hsum));
				const T final_res = _mm_cvtsi128_si32(_mm_packus_epi16(_mm_packus_epi32(_mm_cvtps_epi32(hsum), _mm_setzero_si128()), _mm_setzero_si128()));

				dst[x] = clamp(final_res,val_min,val_max);
            }
            else if JincMT_CONSTEXPR (std::is_same<T, uint16_t>::value)
            {
                for (int ly = 0; ly < filter_size_mod2; ly += 2)
                {
                    for (int lx = 0; lx < filter_size; lx += 16)
                    {
                        const __m512 src_ps = _mm512_cvtepi32_ps(_mm512_cvtepu16_epi32(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(src_ptr + lx))));
						const __m512 src_ps2 = _mm512_cvtepi32_ps(_mm512_cvtepu16_epi32(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(src_ptr + lx + src_pitch))));
						__m512 coeff, coeff2;

						if JincMT_CONSTEXPR (bFP16)
						{
							const int lx2 = lx >> 1;
							const __m256i coeff_fp16 = _mm256_load_si256((__m256i*)(coeff_ptr + lx2));
							const __m256i coeff2_fp16 = _mm256_load_si256((__m256i*)(coeff_ptr + lx2 + coeff_stride));

							coeff = _mm512_cvtph_ps(coeff_fp16); // separated grouped converts may be better for instructions grouping at some compilers
							coeff2 = _mm512_cvtph_ps(coeff2_fp16);
						}
						else
						{
							coeff = _mm512_load_ps(coeff_ptr + lx);
							coeff2 = _mm512_load_ps(coeff_ptr + lx + coeff_stride);
						}

						result = _mm512_fmadd_ps(src_ps, coeff, result);
						result2 = _mm512_fmadd_ps(src_ps2, coeff2, result2);
					}

                    coeff_ptr += coeff_stridex2;
                    src_ptr += src_pitchx2;
                }

				result = _mm512_add_ps(result, result2);

				// last row
				if (fs_notMod2)
				{
					for (int lx = 0; lx < filter_size; lx += 16)
					{
						const __m512 src_ps = _mm512_cvtepi32_ps(_mm512_cvtepu16_epi32(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(src_ptr + lx))));
						__m512 coeff;

						if JincMT_CONSTEXPR (bFP16)
							coeff = _mm512_cvtph_ps(_mm256_load_si256((__m256i*)(coeff_ptr + (lx >> 1))));
						else
							coeff = _mm512_load_ps(coeff_ptr + lx);

						result = _mm512_fmadd_ps(src_ps, coeff, result);
					}
				}

                const __m256 lo_hi_256 = _mm256_add_ps(_mm512_castps512_ps256(result), _mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(result), 1)));
                __m128 hsum = _mm_add_ps(_mm256_castps256_ps128(lo_hi_256), _mm256_extractf128_ps(lo_hi_256, 1));
                hsum = _mm_hadd_ps(_mm_hadd_ps(hsum, hsum), _mm_hadd_ps(hsum, hsum));
				const T final_res = _mm_cvtsi128_si32(_mm_packus_epi32(_mm_cvtps_epi32(hsum), _mm_setzero_si128()));

				dst[x] = clamp(final_res,val_min,val_max);
            }
            else
            {
                for (int ly = 0; ly < filter_size_mod2; ly += 2)
                {
                    for (int lx = 0; lx < filter_size; lx += 16)
                    {
                        const __m512 src_ps = _mm512_max_ps(_mm512_loadu_ps(src_ptr + lx), min_val);
						const __m512 src_ps2 = _mm512_max_ps(_mm512_loadu_ps(src_ptr + lx + src_pitch), min_val);
						__m512 coeff, coeff2;

						if JincMT_CONSTEXPR (bFP16)
						{
							const int lx2 = lx >> 1;
							const __m256i coeff_fp16 = _mm256_load_si256((__m256i*)(coeff_ptr + lx2));
							const __m256i coeff2_fp16 = _mm256_load_si256((__m256i*)(coeff_ptr + lx2 + coeff_stride));

							coeff = _mm512_cvtph_ps(coeff_fp16); // separated grouped converts may be better for instructions grouping at some compilers
							coeff2 = _mm512_cvtph_ps(coeff2_fp16);
						}
						else
						{
							coeff = _mm512_load_ps(coeff_ptr + lx);
							coeff2 = _mm512_load_ps(coeff_ptr + lx + coeff_stride);
						}

						result = _mm512_fmadd_ps(src_ps, coeff, result);
						result2 = _mm512_fmadd_ps(src_ps2, coeff2, result2);
                    }

                    coeff_ptr += coeff_stridex2;
                    src_ptr += src_pitchx2;
                }

				result = _mm512_add_ps(result, result2);

				// last row
				if (fs_notMod2)
				{
					for (int lx = 0; lx < filter_size; lx += 16)
					{
						const __m512 src_ps = _mm512_max_ps(_mm512_loadu_ps(src_ptr + lx), min_val);
						__m512 coeff;

						if JincMT_CONSTEXPR (bFP16)
							coeff = _mm512_cvtph_ps(_mm256_load_si256((__m256i*)(coeff_ptr + (lx >> 1))));
						else
							coeff = _mm512_load_ps(coeff_ptr + lx);

						result = _mm512_fmadd_ps(src_ps, coeff, result);
					}
				}

                const __m256 lo_hi_256 = _mm256_add_ps(_mm512_castps512_ps256(result), _mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(result), 1)));
                __m128 hsum = _mm_add_ps(_mm256_castps256_ps128(lo_hi_256), _mm256_extractf128_ps(lo_hi_256, 1));

				dst[x] = _mm_cvtss_f32(_mm_hadd_ps(_mm_hadd_ps(hsum, hsum), _mm_hadd_ps(hsum, hsum)));
            }
			meta++;
        } // for (x)
		meta_y += dst_width;
        dst += dst_pitch;
	} // for (y)
}


template <typename T, bool bFP16>
#if defined(CLANG)
__attribute__((__target__("avx512f")))
#endif
void resize_plane_avx512_2x(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *tab_coeff,
	const float Val_Min[], const float Val_Max[])
{
	const uint8_t idx1 = (PlaneYMode) ? 0 : 1;
	const uint8_t idx2 = (PlaneYMode) ? 3 : 2;

	const T *src1 = reinterpret_cast<const T*>(MT_DataGF->src[idx1]);
	const T *src2 = reinterpret_cast<const T*>(MT_DataGF->src[idx2]);
	T *JincMT_RESTRICT dst1 = reinterpret_cast<T*>(MT_DataGF->dst[idx1]);
	T *JincMT_RESTRICT dst2 = reinterpret_cast<T*>(MT_DataGF->dst[idx2]);

	const ptrdiff_t src_pitch1 = (ptrdiff_t)MT_DataGF->src_pitch[idx1] / sizeof(T);
	const ptrdiff_t src_pitch2 = (ptrdiff_t)MT_DataGF->src_pitch[idx2] / sizeof(T);
	const ptrdiff_t dst_pitch1 = (ptrdiff_t)MT_DataGF->dst_pitch[idx1] / sizeof(T);
	const ptrdiff_t dst_pitch2 = (ptrdiff_t)MT_DataGF->dst_pitch[idx2] / sizeof(T);
	const ptrdiff_t src_pitch1x2 = src_pitch1 << 1;
	const ptrdiff_t src_pitch2x2 = src_pitch2 << 1;

	const int Y_Min = (PlaneYMode) ? MT_DataGF->dst_Y_h_min : MT_DataGF->dst_UV_h_min;
	const int Y_Max = (PlaneYMode) ? MT_DataGF->dst_Y_h_max : MT_DataGF->dst_UV_h_max;
	const int dst_width = (PlaneYMode) ? MT_DataGF->dst_Y_w : MT_DataGF->dst_UV_w;

	const T val_min1 = static_cast<T>(Val_Min[idx1]);
	const T val_min2 = static_cast<T>(Val_Min[idx2]);
	const T val_max1 = static_cast<T>(Val_Max[idx1]);
	const T val_max2 = static_cast<T>(Val_Max[idx2]);
	const __m512 min_val1 = _mm512_set1_ps(Val_Min[idx1]);
	const __m512 min_val2 = _mm512_set1_ps(Val_Min[idx2]);

	EWAPixelCoeffMeta *meta_y = tab_coeff->meta + (Y_Min*dst_width);

	const int filter_size = tab_coeff->filter_size, coeff_stride = tab_coeff->coeff_stride;
	const int coeff_stridex2 = coeff_stride << 1;

	const int filter_size_mod2 = (filter_size >> 1) << 1;
	const bool fs_notMod2 = filter_size_mod2 < filter_size;

	for (int y = Y_Min; y < Y_Max; y++)
	{
		EWAPixelCoeffMeta *meta = meta_y;

		for (int x = 0; x < dst_width; ++x)
		{
			const T *src_ptr1 = src1 + (meta->start_y * src_pitch1 + meta->start_x);
			const T *src_ptr2 = src2 + (meta->start_y * src_pitch2 + meta->start_x);
			const float *coeff_ptr = tab_coeff->factor + meta->coeff_meta;
			__m512 result1 = _mm512_setzero_ps();
			__m512 result2 = _mm512_setzero_ps();
			__m512 result1_2 = _mm512_setzero_ps();
			__m512 result2_2 = _mm512_setzero_ps();

			if JincMT_CONSTEXPR (std::is_same<T, uint8_t>::value)
			{
				for (int ly = 0; ly < filter_size_mod2; ly += 2)
				{
					for (int lx = 0; lx < filter_size; lx += 16)
					{
						const __m512 src_ps1 = _mm512_cvtepi32_ps(_mm512_cvtepu8_epi32(_mm_loadu_si128(reinterpret_cast<const __m128i*>(src_ptr1 + lx))));
						const __m512 src_ps1_2 = _mm512_cvtepi32_ps(_mm512_cvtepu8_epi32(_mm_loadu_si128(reinterpret_cast<const __m128i*>(src_ptr1 + lx + src_pitch1))));
						const __m512 src_ps2 = _mm512_cvtepi32_ps(_mm512_cvtepu8_epi32(_mm_loadu_si128(reinterpret_cast<const __m128i*>(src_ptr2 + lx))));
						const __m512 src_ps2_2 = _mm512_cvtepi32_ps(_mm512_cvtepu8_epi32(_mm_loadu_si128(reinterpret_cast<const __m128i*>(src_ptr2 + lx + src_pitch2))));
						__m512 coeff, coeff2;

						if JincMT_CONSTEXPR (bFP16)
						{
							const int lx2 = lx >> 1;
							const __m256i coeff_fp16 = _mm256_load_si256((__m256i*)(coeff_ptr + lx2));
							const __m256i coeff2_fp16 = _mm256_load_si256((__m256i*)(coeff_ptr + lx2 + coeff_stride));

							coeff = _mm512_cvtph_ps(coeff_fp16); // separated grouped converts may be better for instructions grouping at some compilers
							coeff2 = _mm512_cvtph_ps(coeff2_fp16);
						}
						else
						{
							coeff = _mm512_load_ps(coeff_ptr + lx);
							coeff2 = _mm512_load_ps(coeff_ptr + lx + coeff_stride);
						}

						result1 = _mm512_fmadd_ps(src_ps1, coeff, result1);
						result1_2 = _mm512_fmadd_ps(src_ps1_2, coeff2, result1_2);
						result2 = _mm512_fmadd_ps(src_ps2, coeff, result2);
						result2_2 = _mm512_fmadd_ps(src_ps2_2, coeff2, result2_2);
					}

					coeff_ptr += coeff_stridex2;
					src_ptr1 += src_pitch1x2;
					src_ptr2 += src_pitch2x2;
				}

				result1 = _mm512_add_ps(result1, result1_2);
				result2 = _mm512_add_ps(result2, result2_2);

				// last row
				if (fs_notMod2)
				{
					for (int lx = 0; lx < filter_size; lx += 16)
					{
						const __m512 src_ps1 = _mm512_cvtepi32_ps(_mm512_cvtepu8_epi32(_mm_loadu_si128(reinterpret_cast<const __m128i*>(src_ptr1 + lx))));
						const __m512 src_ps2 = _mm512_cvtepi32_ps(_mm512_cvtepu8_epi32(_mm_loadu_si128(reinterpret_cast<const __m128i*>(src_ptr2 + lx))));
						__m512 coeff;

						if JincMT_CONSTEXPR (bFP16)
							coeff = _mm512_cvtph_ps(_mm256_load_si256((__m256i*)(coeff_ptr + (lx >> 1))));
						else
							coeff = _mm512_load_ps(coeff_ptr + lx);

						result1 = _mm512_fmadd_ps(src_ps1, coeff, result1);
						result2 = _mm512_fmadd_ps(src_ps2, coeff, result2);
					}
				}

				const __m256 lo_hi_256_1 = _mm256_add_ps(_mm512_castps512_ps256(result1), _mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(result1), 1)));
				const __m256 lo_hi_256_2 = _mm256_add_ps(_mm512_castps512_ps256(result2), _mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(result2), 1)));
				__m128 hsum1 = _mm_add_ps(_mm256_castps256_ps128(lo_hi_256_1), _mm256_extractf128_ps(lo_hi_256_1, 1));
				__m128 hsum2 = _mm_add_ps(_mm256_castps256_ps128(lo_hi_256_2), _mm256_extractf128_ps(lo_hi_256_2, 1));
				hsum1 = _mm_hadd_ps(_mm_hadd_ps(hsum1, hsum1), _mm_hadd_ps(hsum1, hsum1));
				hsum2 = _mm_hadd_ps(_mm_hadd_ps(hsum2, hsum2), _mm_hadd_ps(hsum2, hsum2));
				const T final_res1 = _mm_cvtsi128_si32(_mm_packus_epi16(_mm_packus_epi32(_mm_cvtps_epi32(hsum1), _mm_setzero_si128()), _mm_setzero_si128()));
				const T final_res2 = _mm_cvtsi128_si32(_mm_packus_epi16(_mm_packus_epi32(_mm_cvtps_epi32(hsum2), _mm_setzero_si128()), _mm_setzero_si128()));

				dst1[x] = clamp(final_res1, val_min1, val_max1);
				dst2[x] = clamp(final_res2, val_min2, val_max2);
			}
			else if JincMT_CONSTEXPR (std::is_same<T, uint16_t>::value)
			{
				for (int ly = 0; ly < filter_size_mod2; ly += 2)
				{
					for (int lx = 0; lx < filter_size; lx += 16)
					{
						const __m512 src_ps1 = _mm512_cvtepi32_ps(_mm512_cvtepu16_epi32(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(src_ptr1 + lx))));
						const __m512 src_ps1_2 = _mm512_cvtepi32_ps(_mm512_cvtepu16_epi32(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(src_ptr1 + lx + src_pitch1))));
						const __m512 src_ps2 = _mm512_cvtepi32_ps(_mm512_cvtepu16_epi32(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(src_ptr2 + lx))));
						const __m512 src_ps2_2 = _mm512_cvtepi32_ps(_mm512_cvtepu16_epi32(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(src_ptr2 + lx + src_pitch2))));
						__m512 coeff, coeff2;

						if JincMT_CONSTEXPR (bFP16)
						{
							const int lx2 = lx >> 1;
							const __m256i coeff_fp16 = _mm256_load_si256((__m256i*)(coeff_ptr + lx2));
							const __m256i coeff2_fp16 = _mm256_load_si256((__m256i*)(coeff_ptr + lx2 + coeff_stride));

							coeff = _mm512_cvtph_ps(coeff_fp16); // separated grouped converts may be better for instructions grouping at some compilers
							coeff2 = _mm512_cvtph_ps(coeff2_fp16);
						}
						else
						{
							coeff = _mm512_load_ps(coeff_ptr + lx);
							coeff2 = _mm512_load_ps(coeff_ptr + lx + coeff_stride);
						}

						result1 = _mm512_fmadd_ps(src_ps1, coeff, result1);
						result1_2 = _mm512_fmadd_ps(src_ps1_2, coeff2, result1_2);
						result2 = _mm512_fmadd_ps(src_ps2, coeff, result2);
						result2_2 = _mm512_fmadd_ps(src_ps2_2, coeff2, result2_2);
					}

					coeff_ptr += coeff_stridex2;
					src_ptr1 += src_pitch1x2;
					src_ptr2 += src_pitch2x2;
				}

				result1 = _mm512_add_ps(result1, result1_2);
				result2 = _mm512_add_ps(result2, result2_2);

				// last row
				if (fs_notMod2)
				{
					for (int lx = 0; lx < filter_size; lx += 16)
					{
						const __m512 src_ps1 = _mm512_cvtepi32_ps(_mm512_cvtepu16_epi32(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(src_ptr1 + lx))));
						const __m512 src_ps2 = _mm512_cvtepi32_ps(_mm512_cvtepu16_epi32(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(src_ptr2 + lx))));
						__m512 coeff;

						if JincMT_CONSTEXPR (bFP16)
							coeff = _mm512_cvtph_ps(_mm256_load_si256((__m256i*)(coeff_ptr + (lx >> 1))));
						else
							coeff = _mm512_load_ps(coeff_ptr + lx);

						result1 = _mm512_fmadd_ps(src_ps1, coeff, result1);
						result2 = _mm512_fmadd_ps(src_ps2, coeff, result2);
					}
				}

				const __m256 lo_hi_256_1 = _mm256_add_ps(_mm512_castps512_ps256(result1), _mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(result1), 1)));
				const __m256 lo_hi_256_2 = _mm256_add_ps(_mm512_castps512_ps256(result2), _mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(result2), 1)));
				__m128 hsum1 = _mm_add_ps(_mm256_castps256_ps128(lo_hi_256_1), _mm256_extractf128_ps(lo_hi_256_1, 1));
				__m128 hsum2 = _mm_add_ps(_mm256_castps256_ps128(lo_hi_256_2), _mm256_extractf128_ps(lo_hi_256_2, 1));
				hsum1 = _mm_hadd_ps(_mm_hadd_ps(hsum1, hsum1), _mm_hadd_ps(hsum1, hsum1));
				hsum2 = _mm_hadd_ps(_mm_hadd_ps(hsum2, hsum2), _mm_hadd_ps(hsum2, hsum2));
				const T final_res1 = _mm_cvtsi128_si32(_mm_packus_epi32(_mm_cvtps_epi32(hsum1), _mm_setzero_si128()));
				const T final_res2 = _mm_cvtsi128_si32(_mm_packus_epi32(_mm_cvtps_epi32(hsum2), _mm_setzero_si128()));

				dst1[x] = clamp(final_res1, val_min1, val_max1);
				dst2[x] = clamp(final_res2, val_min2, val_max2);
			}
			else
			{
				for (int ly = 0; ly < filter_size_mod2; ly += 2)
				{
					for (int lx = 0; lx < filter_size; lx += 16)
					{
						const __m512 src_ps1 = _mm512_max_ps(_mm512_loadu_ps(src_ptr1 + lx), min_val1);
						const __m512 src_ps1_2 = _mm512_max_ps(_mm512_loadu_ps(src_ptr1 + lx + src_pitch1), min_val1);
						const __m512 src_ps2 = _mm512_max_ps(_mm512_loadu_ps(src_ptr2 + lx), min_val2);
						const __m512 src_ps2_2 = _mm512_max_ps(_mm512_loadu_ps(src_ptr2 + lx + src_pitch2), min_val2);
						__m512 coeff, coeff2;

						if JincMT_CONSTEXPR (bFP16)
						{
							const int lx2 = lx >> 1;
							const __m256i coeff_fp16 = _mm256_load_si256((__m256i*)(coeff_ptr + lx2));
							const __m256i coeff2_fp16 = _mm256_load_si256((__m256i*)(coeff_ptr + lx2 + coeff_stride));

							coeff = _mm512_cvtph_ps(coeff_fp16); // separated grouped converts may be better for instructions grouping at some compilers
							coeff2 = _mm512_cvtph_ps(coeff2_fp16);
						}
						else
						{
							coeff = _mm512_load_ps(coeff_ptr + lx);
							coeff2 = _mm512_load_ps(coeff_ptr + lx + coeff_stride);
						}

						result1 = _mm512_fmadd_ps(src_ps1, coeff, result1);
						result1_2 = _mm512_fmadd_ps(src_ps1_2, coeff2, result1_2);
						result2 = _mm512_fmadd_ps(src_ps2, coeff, result2);
						result2_2 = _mm512_fmadd_ps(src_ps2_2, coeff2, result2_2);
					}

					coeff_ptr += coeff_stridex2;
					src_ptr1 += src_pitch1x2;
					src_ptr2 += src_pitch2x2;
				}

				result1 = _mm512_add_ps(result1, result1_2);
				result2 = _mm512_add_ps(result2, result2_2);

				// last row
				if (fs_notMod2)
				{
					for (int lx = 0; lx < filter_size; lx += 16)
					{
						const __m512 src_ps1 = _mm512_max_ps(_mm512_loadu_ps(src_ptr1 + lx), min_val1);
						const __m512 src_ps2 = _mm512_max_ps(_mm512_loadu_ps(src_ptr2 + lx), min_val2);
						__m512 coeff;

						if JincMT_CONSTEXPR (bFP16)
							coeff = _mm512_cvtph_ps(_mm256_load_si256((__m256i*)(coeff_ptr + (lx >> 1))));
						else
							coeff = _mm512_load_ps(coeff_ptr + lx);

						result1 = _mm512_fmadd_ps(src_ps1, coeff, result1);
						result2 = _mm512_fmadd_ps(src_ps2, coeff, result2);
					}
				}

				const __m256 lo_hi_256_1 = _mm256_add_ps(_mm512_castps512_ps256(result1), _mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(result1), 1)));
				const __m256 lo_hi_256_2 = _mm256_add_ps(_mm512_castps512_ps256(result2), _mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(result2), 1)));
				__m128 hsum1 = _mm_add_ps(_mm256_castps256_ps128(lo_hi_256_1), _mm256_extractf128_ps(lo_hi_256_1, 1));
				__m128 hsum2 = _mm_add_ps(_mm256_castps256_ps128(lo_hi_256_2), _mm256_extractf128_ps(lo_hi_256_2, 1));

				dst1[x] = _mm_cvtss_f32(_mm_hadd_ps(_mm_hadd_ps(hsum1, hsum1), _mm_hadd_ps(hsum1, hsum1)));
				dst2[x] = _mm_cvtss_f32(_mm_hadd_ps(_mm_hadd_ps(hsum2, hsum2), _mm_hadd_ps(hsum2, hsum2)));
			}
			meta++;
		} // for (x)
		meta_y += dst_width;
		dst1 += dst_pitch1;
		dst2 += dst_pitch2;
	} // for (y)
}


template <typename T, bool bFP16>
#if defined(CLANG)
__attribute__((__target__("avx512f")))
#endif
void resize_plane_avx512_3x(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *tab_coeff,
	const float Val_Min[], const float Val_Max[])
{
	const uint8_t idx1 = 0;
	const uint8_t idx2 = 1;
	const uint8_t idx3 = 2;

	const T *src1 = reinterpret_cast<const T*>(MT_DataGF->src[idx1]);
	const T *src2 = reinterpret_cast<const T*>(MT_DataGF->src[idx2]);
	const T *src3 = reinterpret_cast<const T*>(MT_DataGF->src[idx3]);
	T *JincMT_RESTRICT dst1 = reinterpret_cast<T*>(MT_DataGF->dst[idx1]);
	T *JincMT_RESTRICT dst2 = reinterpret_cast<T*>(MT_DataGF->dst[idx2]);
	T *JincMT_RESTRICT dst3 = reinterpret_cast<T*>(MT_DataGF->dst[idx3]);

	const ptrdiff_t src_pitch1 = (ptrdiff_t)MT_DataGF->src_pitch[idx1] / sizeof(T);
	const ptrdiff_t src_pitch2 = (ptrdiff_t)MT_DataGF->src_pitch[idx2] / sizeof(T);
	const ptrdiff_t src_pitch3 = (ptrdiff_t)MT_DataGF->src_pitch[idx3] / sizeof(T);
	const ptrdiff_t dst_pitch1 = (ptrdiff_t)MT_DataGF->dst_pitch[idx1] / sizeof(T);
	const ptrdiff_t dst_pitch2 = (ptrdiff_t)MT_DataGF->dst_pitch[idx2] / sizeof(T);
	const ptrdiff_t dst_pitch3 = (ptrdiff_t)MT_DataGF->dst_pitch[idx3] / sizeof(T);
	const ptrdiff_t src_pitch1x2 = src_pitch1 << 1;
	const ptrdiff_t src_pitch2x2 = src_pitch2 << 1;
	const ptrdiff_t src_pitch3x2 = src_pitch3 << 1;

	const int Y_Min = MT_DataGF->dst_Y_h_min;
	const int Y_Max = MT_DataGF->dst_Y_h_max;
	const int dst_width = MT_DataGF->dst_Y_w;

	const T val_min1 = static_cast<T>(Val_Min[idx1]);
	const T val_min2 = static_cast<T>(Val_Min[idx2]);
	const T val_min3 = static_cast<T>(Val_Min[idx3]);
	const T val_max1 = static_cast<T>(Val_Max[idx1]);
	const T val_max2 = static_cast<T>(Val_Max[idx2]);
	const T val_max3 = static_cast<T>(Val_Max[idx3]);
	const __m512 min_val1 = _mm512_set1_ps(Val_Min[idx1]);
	const __m512 min_val2 = _mm512_set1_ps(Val_Min[idx2]);
	const __m512 min_val3 = _mm512_set1_ps(Val_Min[idx3]);

	EWAPixelCoeffMeta *meta_y = tab_coeff->meta + (Y_Min*dst_width);

	const int filter_size = tab_coeff->filter_size, coeff_stride = tab_coeff->coeff_stride;
	const int coeff_stridex2 = coeff_stride << 1;

	const int filter_size_mod2 = (filter_size >> 1) << 1;
	const bool fs_notMod2 = filter_size_mod2 < filter_size;

	for (int y = Y_Min; y < Y_Max; y++)
	{
		EWAPixelCoeffMeta *meta = meta_y;

		for (int x = 0; x < dst_width; ++x)
		{
			const T *src_ptr1 = src1 + (meta->start_y * src_pitch1 + meta->start_x);
			const T *src_ptr2 = src2 + (meta->start_y * src_pitch2 + meta->start_x);
			const T *src_ptr3 = src3 + (meta->start_y * src_pitch3 + meta->start_x);
			const float *coeff_ptr = tab_coeff->factor + meta->coeff_meta;
			__m512 result1 = _mm512_setzero_ps();
			__m512 result2 = _mm512_setzero_ps();
			__m512 result3 = _mm512_setzero_ps();
			__m512 result1_2 = _mm512_setzero_ps();
			__m512 result2_2 = _mm512_setzero_ps();
			__m512 result3_2 = _mm512_setzero_ps();

			if JincMT_CONSTEXPR (std::is_same<T, uint8_t>::value)
			{
				for (int ly = 0; ly < filter_size_mod2; ly += 2)
				{
					for (int lx = 0; lx < filter_size; lx += 16)
					{
						const __m512 src_ps1 = _mm512_cvtepi32_ps(_mm512_cvtepu8_epi32(_mm_loadu_si128(reinterpret_cast<const __m128i*>(src_ptr1 + lx))));
						const __m512 src_ps1_2 = _mm512_cvtepi32_ps(_mm512_cvtepu8_epi32(_mm_loadu_si128(reinterpret_cast<const __m128i*>(src_ptr1 + lx + src_pitch1))));
						const __m512 src_ps2 = _mm512_cvtepi32_ps(_mm512_cvtepu8_epi32(_mm_loadu_si128(reinterpret_cast<const __m128i*>(src_ptr2 + lx))));
						const __m512 src_ps2_2 = _mm512_cvtepi32_ps(_mm512_cvtepu8_epi32(_mm_loadu_si128(reinterpret_cast<const __m128i*>(src_ptr2 + lx + src_pitch2))));
						const __m512 src_ps3 = _mm512_cvtepi32_ps(_mm512_cvtepu8_epi32(_mm_loadu_si128(reinterpret_cast<const __m128i*>(src_ptr3 + lx))));
						const __m512 src_ps3_2 = _mm512_cvtepi32_ps(_mm512_cvtepu8_epi32(_mm_loadu_si128(reinterpret_cast<const __m128i*>(src_ptr3 + lx + src_pitch3))));
						__m512 coeff, coeff2;

						if JincMT_CONSTEXPR (bFP16)
						{
							const int lx2 = lx >> 1;
							const __m256i coeff_fp16 = _mm256_load_si256((__m256i*)(coeff_ptr + lx2));
							const __m256i coeff2_fp16 = _mm256_load_si256((__m256i*)(coeff_ptr + lx2 + coeff_stride));

							coeff = _mm512_cvtph_ps(coeff_fp16); // separated grouped converts may be better for instructions grouping at some compilers
							coeff2 = _mm512_cvtph_ps(coeff2_fp16);
						}
						else
						{
							coeff = _mm512_load_ps(coeff_ptr + lx);
							coeff2 = _mm512_load_ps(coeff_ptr + lx + coeff_stride);
						}

						result1 = _mm512_fmadd_ps(src_ps1, coeff, result1);
						result1_2 = _mm512_fmadd_ps(src_ps1_2, coeff2, result1_2);
						result2 = _mm512_fmadd_ps(src_ps2, coeff, result2);
						result2_2 = _mm512_fmadd_ps(src_ps2_2, coeff2, result2_2);
						result3 = _mm512_fmadd_ps(src_ps3, coeff, result3);
						result3_2 = _mm512_fmadd_ps(src_ps3_2, coeff2, result3_2);
					}

					coeff_ptr += coeff_stridex2;
					src_ptr1 += src_pitch1x2;
					src_ptr2 += src_pitch2x2;
					src_ptr3 += src_pitch3x2;
				}

				result1 = _mm512_add_ps(result1, result1_2);
				result2 = _mm512_add_ps(result2, result2_2);
				result3 = _mm512_add_ps(result3, result3_2);

				// last row
				if (fs_notMod2)
				{
					for (int lx = 0; lx < filter_size; lx += 16)
					{
						const __m512 src_ps1 = _mm512_cvtepi32_ps(_mm512_cvtepu8_epi32(_mm_loadu_si128(reinterpret_cast<const __m128i*>(src_ptr1 + lx))));
						const __m512 src_ps2 = _mm512_cvtepi32_ps(_mm512_cvtepu8_epi32(_mm_loadu_si128(reinterpret_cast<const __m128i*>(src_ptr2 + lx))));
						const __m512 src_ps3 = _mm512_cvtepi32_ps(_mm512_cvtepu8_epi32(_mm_loadu_si128(reinterpret_cast<const __m128i*>(src_ptr3 + lx))));
						__m512 coeff;

						if JincMT_CONSTEXPR (bFP16)
							coeff = _mm512_cvtph_ps(_mm256_load_si256((__m256i*)(coeff_ptr + (lx >> 1))));
						else
							coeff = _mm512_load_ps(coeff_ptr + lx);

						result1 = _mm512_fmadd_ps(src_ps1, coeff, result1);
						result2 = _mm512_fmadd_ps(src_ps2, coeff, result2);
						result3 = _mm512_fmadd_ps(src_ps3, coeff, result3);
					}
				}

				const __m256 lo_hi_256_1 = _mm256_add_ps(_mm512_castps512_ps256(result1), _mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(result1), 1)));
				const __m256 lo_hi_256_2 = _mm256_add_ps(_mm512_castps512_ps256(result2), _mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(result2), 1)));
				const __m256 lo_hi_256_3 = _mm256_add_ps(_mm512_castps512_ps256(result3), _mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(result3), 1)));
				__m128 hsum1 = _mm_add_ps(_mm256_castps256_ps128(lo_hi_256_1), _mm256_extractf128_ps(lo_hi_256_1, 1));
				__m128 hsum2 = _mm_add_ps(_mm256_castps256_ps128(lo_hi_256_2), _mm256_extractf128_ps(lo_hi_256_2, 1));
				__m128 hsum3 = _mm_add_ps(_mm256_castps256_ps128(lo_hi_256_3), _mm256_extractf128_ps(lo_hi_256_3, 1));
				hsum1 = _mm_hadd_ps(_mm_hadd_ps(hsum1, hsum1), _mm_hadd_ps(hsum1, hsum1));
				hsum2 = _mm_hadd_ps(_mm_hadd_ps(hsum2, hsum2), _mm_hadd_ps(hsum2, hsum2));
				hsum3 = _mm_hadd_ps(_mm_hadd_ps(hsum3, hsum3), _mm_hadd_ps(hsum3, hsum3));
				const T final_res1 = _mm_cvtsi128_si32(_mm_packus_epi16(_mm_packus_epi32(_mm_cvtps_epi32(hsum1), _mm_setzero_si128()), _mm_setzero_si128()));
				const T final_res2 = _mm_cvtsi128_si32(_mm_packus_epi16(_mm_packus_epi32(_mm_cvtps_epi32(hsum2), _mm_setzero_si128()), _mm_setzero_si128()));
				const T final_res3 = _mm_cvtsi128_si32(_mm_packus_epi16(_mm_packus_epi32(_mm_cvtps_epi32(hsum3), _mm_setzero_si128()), _mm_setzero_si128()));

				dst1[x] = clamp(final_res1, val_min1, val_max1);
				dst2[x] = clamp(final_res2, val_min2, val_max2);
				dst3[x] = clamp(final_res3, val_min3, val_max3);
			}
			else if JincMT_CONSTEXPR (std::is_same<T, uint16_t>::value)
			{
				for (int ly = 0; ly < filter_size_mod2; ly += 2)
				{
					for (int lx = 0; lx < filter_size; lx += 16)
					{
						const __m512 src_ps1 = _mm512_cvtepi32_ps(_mm512_cvtepu16_epi32(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(src_ptr1 + lx))));
						const __m512 src_ps1_2 = _mm512_cvtepi32_ps(_mm512_cvtepu16_epi32(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(src_ptr1 + lx + src_pitch1))));
						const __m512 src_ps2 = _mm512_cvtepi32_ps(_mm512_cvtepu16_epi32(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(src_ptr2 + lx))));
						const __m512 src_ps2_2 = _mm512_cvtepi32_ps(_mm512_cvtepu16_epi32(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(src_ptr2 + lx + src_pitch2))));
						const __m512 src_ps3 = _mm512_cvtepi32_ps(_mm512_cvtepu16_epi32(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(src_ptr3 + lx))));
						const __m512 src_ps3_2 = _mm512_cvtepi32_ps(_mm512_cvtepu16_epi32(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(src_ptr3 + lx + src_pitch3))));
						__m512 coeff, coeff2;

						if JincMT_CONSTEXPR (bFP16)
						{
							const int lx2 = lx >> 1;
							const __m256i coeff_fp16 = _mm256_load_si256((__m256i*)(coeff_ptr + lx2));
							const __m256i coeff2_fp16 = _mm256_load_si256((__m256i*)(coeff_ptr + lx2 + coeff_stride));

							coeff = _mm512_cvtph_ps(coeff_fp16); // separated grouped converts may be better for instructions grouping at some compilers
							coeff2 = _mm512_cvtph_ps(coeff2_fp16);
						}
						else
						{
							coeff = _mm512_load_ps(coeff_ptr + lx);
							coeff2 = _mm512_load_ps(coeff_ptr + lx + coeff_stride);
						}

						result1 = _mm512_fmadd_ps(src_ps1, coeff, result1);
						result1_2 = _mm512_fmadd_ps(src_ps1_2, coeff2, result1_2);
						result2 = _mm512_fmadd_ps(src_ps2, coeff, result2);
						result2_2 = _mm512_fmadd_ps(src_ps2_2, coeff2, result2_2);
						result3 = _mm512_fmadd_ps(src_ps3, coeff, result3);
						result3_2 = _mm512_fmadd_ps(src_ps3_2, coeff2, result3_2);
					}

					coeff_ptr += coeff_stridex2;
					src_ptr1 += src_pitch1x2;
					src_ptr2 += src_pitch2x2;
					src_ptr3 += src_pitch3x2;
				}

				result1 = _mm512_add_ps(result1, result1_2);
				result2 = _mm512_add_ps(result2, result2_2);
				result3 = _mm512_add_ps(result3, result3_2);

				// last row
				if (fs_notMod2)
				{
					for (int lx = 0; lx < filter_size; lx += 16)
					{
						const __m512 src_ps1 = _mm512_cvtepi32_ps(_mm512_cvtepu16_epi32(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(src_ptr1 + lx))));
						const __m512 src_ps2 = _mm512_cvtepi32_ps(_mm512_cvtepu16_epi32(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(src_ptr2 + lx))));
						const __m512 src_ps3 = _mm512_cvtepi32_ps(_mm512_cvtepu16_epi32(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(src_ptr3 + lx))));
						__m512 coeff;

						if JincMT_CONSTEXPR (bFP16)
							coeff = _mm512_cvtph_ps(_mm256_load_si256((__m256i*)(coeff_ptr + (lx >> 1))));
						else
							coeff = _mm512_load_ps(coeff_ptr + lx);

						result1 = _mm512_fmadd_ps(src_ps1, coeff, result1);
						result2 = _mm512_fmadd_ps(src_ps2, coeff, result2);
						result3 = _mm512_fmadd_ps(src_ps3, coeff, result3);
					}
				}

				const __m256 lo_hi_256_1 = _mm256_add_ps(_mm512_castps512_ps256(result1), _mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(result1), 1)));
				const __m256 lo_hi_256_2 = _mm256_add_ps(_mm512_castps512_ps256(result2), _mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(result2), 1)));
				const __m256 lo_hi_256_3 = _mm256_add_ps(_mm512_castps512_ps256(result3), _mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(result3), 1)));
				__m128 hsum1 = _mm_add_ps(_mm256_castps256_ps128(lo_hi_256_1), _mm256_extractf128_ps(lo_hi_256_1, 1));
				__m128 hsum2 = _mm_add_ps(_mm256_castps256_ps128(lo_hi_256_2), _mm256_extractf128_ps(lo_hi_256_2, 1));
				__m128 hsum3 = _mm_add_ps(_mm256_castps256_ps128(lo_hi_256_3), _mm256_extractf128_ps(lo_hi_256_3, 1));
				hsum1 = _mm_hadd_ps(_mm_hadd_ps(hsum1, hsum1), _mm_hadd_ps(hsum1, hsum1));
				hsum2 = _mm_hadd_ps(_mm_hadd_ps(hsum2, hsum2), _mm_hadd_ps(hsum2, hsum2));
				hsum3 = _mm_hadd_ps(_mm_hadd_ps(hsum3, hsum3), _mm_hadd_ps(hsum3, hsum3));
				const T final_res1 = _mm_cvtsi128_si32(_mm_packus_epi32(_mm_cvtps_epi32(hsum1), _mm_setzero_si128()));
				const T final_res2 = _mm_cvtsi128_si32(_mm_packus_epi32(_mm_cvtps_epi32(hsum2), _mm_setzero_si128()));
				const T final_res3 = _mm_cvtsi128_si32(_mm_packus_epi32(_mm_cvtps_epi32(hsum3), _mm_setzero_si128()));

				dst1[x] = clamp(final_res1, val_min1, val_max1);
				dst2[x] = clamp(final_res2, val_min2, val_max2);
				dst3[x] = clamp(final_res3, val_min3, val_max3);
			}
			else
			{
				for (int ly = 0; ly < filter_size_mod2; ly += 2)
				{
					for (int lx = 0; lx < filter_size; lx += 16)
					{
						const __m512 src_ps1 = _mm512_max_ps(_mm512_loadu_ps(src_ptr1 + lx), min_val1);
						const __m512 src_ps1_2 = _mm512_max_ps(_mm512_loadu_ps(src_ptr1 + lx + src_pitch1), min_val1);
						const __m512 src_ps2 = _mm512_max_ps(_mm512_loadu_ps(src_ptr2 + lx), min_val2);
						const __m512 src_ps2_2 = _mm512_max_ps(_mm512_loadu_ps(src_ptr2 + lx + src_pitch2), min_val2);
						const __m512 src_ps3 = _mm512_max_ps(_mm512_loadu_ps(src_ptr3 + lx), min_val3);
						const __m512 src_ps3_2 = _mm512_max_ps(_mm512_loadu_ps(src_ptr3 + lx + src_pitch3), min_val3);
						__m512 coeff, coeff2;

						if JincMT_CONSTEXPR (bFP16)
						{
							const int lx2 = lx >> 1;
							const __m256i coeff_fp16 = _mm256_load_si256((__m256i*)(coeff_ptr + lx2));
							const __m256i coeff2_fp16 = _mm256_load_si256((__m256i*)(coeff_ptr + lx2 + coeff_stride));

							coeff = _mm512_cvtph_ps(coeff_fp16); // separated grouped converts may be better for instructions grouping at some compilers
							coeff2 = _mm512_cvtph_ps(coeff2_fp16);
						}
						else
						{
							coeff = _mm512_load_ps(coeff_ptr + lx);
							coeff2 = _mm512_load_ps(coeff_ptr + lx + coeff_stride);
						}

						result1 = _mm512_fmadd_ps(src_ps1, coeff, result1);
						result1_2 = _mm512_fmadd_ps(src_ps1_2, coeff2, result1_2);
						result2 = _mm512_fmadd_ps(src_ps2, coeff, result2);
						result2_2 = _mm512_fmadd_ps(src_ps2_2, coeff2, result2_2);
						result3 = _mm512_fmadd_ps(src_ps3, coeff, result3);
						result3_2 = _mm512_fmadd_ps(src_ps3_2, coeff2, result3_2);
					}

					coeff_ptr += coeff_stridex2;
					src_ptr1 += src_pitch1x2;
					src_ptr2 += src_pitch2x2;
					src_ptr3 += src_pitch3x2;
				}

				result1 = _mm512_add_ps(result1, result1_2);
				result2 = _mm512_add_ps(result2, result2_2);
				result3 = _mm512_add_ps(result3, result3_2);

				// last row
				if (fs_notMod2)
				{
					for (int lx = 0; lx < filter_size; lx += 16)
					{
						const __m512 src_ps1 = _mm512_max_ps(_mm512_loadu_ps(src_ptr1 + lx), min_val1);
						const __m512 src_ps2 = _mm512_max_ps(_mm512_loadu_ps(src_ptr2 + lx), min_val2);
						const __m512 src_ps3 = _mm512_max_ps(_mm512_loadu_ps(src_ptr3 + lx), min_val3);
						__m512 coeff;

						if JincMT_CONSTEXPR (bFP16)
							coeff = _mm512_cvtph_ps(_mm256_load_si256((__m256i*)(coeff_ptr + (lx >> 1))));
						else
							coeff = _mm512_load_ps(coeff_ptr + lx);

						result1 = _mm512_fmadd_ps(src_ps1, coeff, result1);
						result2 = _mm512_fmadd_ps(src_ps2, coeff, result2);
						result3 = _mm512_fmadd_ps(src_ps3, coeff, result3);
					}
				}

				const __m256 lo_hi_256_1 = _mm256_add_ps(_mm512_castps512_ps256(result1), _mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(result1), 1)));
				const __m256 lo_hi_256_2 = _mm256_add_ps(_mm512_castps512_ps256(result2), _mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(result2), 1)));
				const __m256 lo_hi_256_3 = _mm256_add_ps(_mm512_castps512_ps256(result3), _mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(result3), 1)));
				__m128 hsum1 = _mm_add_ps(_mm256_castps256_ps128(lo_hi_256_1), _mm256_extractf128_ps(lo_hi_256_1, 1));
				__m128 hsum2 = _mm_add_ps(_mm256_castps256_ps128(lo_hi_256_2), _mm256_extractf128_ps(lo_hi_256_2, 1));
				__m128 hsum3 = _mm_add_ps(_mm256_castps256_ps128(lo_hi_256_3), _mm256_extractf128_ps(lo_hi_256_3, 1));

				dst1[x] = _mm_cvtss_f32(_mm_hadd_ps(_mm_hadd_ps(hsum1, hsum1), _mm_hadd_ps(hsum1, hsum1)));
				dst2[x] = _mm_cvtss_f32(_mm_hadd_ps(_mm_hadd_ps(hsum2, hsum2), _mm_hadd_ps(hsum2, hsum2)));
				dst3[x] = _mm_cvtss_f32(_mm_hadd_ps(_mm_hadd_ps(hsum3, hsum3), _mm_hadd_ps(hsum3, hsum3)));
			}
			meta++;
		} // for (x)
		meta_y += dst_width;
		dst1 += dst_pitch1;
		dst2 += dst_pitch2;
		dst3 += dst_pitch3;
	} // for (y)
}


template <typename T, bool bFP16>
#if defined(CLANG)
__attribute__((__target__("avx512f")))
#endif
void resize_plane_avx512_4x(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *tab_coeff,
	const float Val_Min[], const float Val_Max[])
{
	const uint8_t idx1 = 0;
	const uint8_t idx2 = 1;
	const uint8_t idx3 = 2;
	const uint8_t idx4 = 3;

	const T *src1 = reinterpret_cast<const T*>(MT_DataGF->src[idx1]);
	const T *src2 = reinterpret_cast<const T*>(MT_DataGF->src[idx2]);
	const T *src3 = reinterpret_cast<const T*>(MT_DataGF->src[idx3]);
	const T *src4 = reinterpret_cast<const T*>(MT_DataGF->src[idx4]);
	T *JincMT_RESTRICT dst1 = reinterpret_cast<T*>(MT_DataGF->dst[idx1]);
	T *JincMT_RESTRICT dst2 = reinterpret_cast<T*>(MT_DataGF->dst[idx2]);
	T *JincMT_RESTRICT dst3 = reinterpret_cast<T*>(MT_DataGF->dst[idx3]);
	T *JincMT_RESTRICT dst4 = reinterpret_cast<T*>(MT_DataGF->dst[idx4]);

	const ptrdiff_t src_pitch1 = (ptrdiff_t)MT_DataGF->src_pitch[idx1] / sizeof(T);
	const ptrdiff_t src_pitch2 = (ptrdiff_t)MT_DataGF->src_pitch[idx2] / sizeof(T);
	const ptrdiff_t src_pitch3 = (ptrdiff_t)MT_DataGF->src_pitch[idx3] / sizeof(T);
	const ptrdiff_t src_pitch4 = (ptrdiff_t)MT_DataGF->src_pitch[idx4] / sizeof(T);
	const ptrdiff_t dst_pitch1 = (ptrdiff_t)MT_DataGF->dst_pitch[idx1] / sizeof(T);
	const ptrdiff_t dst_pitch2 = (ptrdiff_t)MT_DataGF->dst_pitch[idx2] / sizeof(T);
	const ptrdiff_t dst_pitch3 = (ptrdiff_t)MT_DataGF->dst_pitch[idx3] / sizeof(T);
	const ptrdiff_t dst_pitch4 = (ptrdiff_t)MT_DataGF->dst_pitch[idx4] / sizeof(T);
	const ptrdiff_t src_pitch1x2 = src_pitch1 << 1;
	const ptrdiff_t src_pitch2x2 = src_pitch2 << 1;
	const ptrdiff_t src_pitch3x2 = src_pitch3 << 1;
	const ptrdiff_t src_pitch4x2 = src_pitch4 << 1;

	const int Y_Min = MT_DataGF->dst_Y_h_min;
	const int Y_Max = MT_DataGF->dst_Y_h_max;
	const int dst_width = MT_DataGF->dst_Y_w;

	const T val_min1 = static_cast<T>(Val_Min[idx1]);
	const T val_min2 = static_cast<T>(Val_Min[idx2]);
	const T val_min3 = static_cast<T>(Val_Min[idx3]);
	const T val_min4 = static_cast<T>(Val_Min[idx4]);
	const T val_max1 = static_cast<T>(Val_Max[idx1]);
	const T val_max2 = static_cast<T>(Val_Max[idx2]);
	const T val_max3 = static_cast<T>(Val_Max[idx3]);
	const T val_max4 = static_cast<T>(Val_Max[idx4]);
	const __m512 min_val1 = _mm512_set1_ps(Val_Min[idx1]);
	const __m512 min_val2 = _mm512_set1_ps(Val_Min[idx2]);
	const __m512 min_val3 = _mm512_set1_ps(Val_Min[idx3]);
	const __m512 min_val4 = _mm512_set1_ps(Val_Min[idx4]);

	EWAPixelCoeffMeta *meta_y = tab_coeff->meta + (Y_Min*dst_width);

	const int filter_size = tab_coeff->filter_size, coeff_stride = tab_coeff->coeff_stride;
	const int coeff_stridex2 = coeff_stride << 1;

	const int filter_size_mod2 = (filter_size >> 1) << 1;
	const bool fs_notMod2 = filter_size_mod2 < filter_size;

	for (int y = Y_Min; y < Y_Max; y++)
	{
		EWAPixelCoeffMeta *meta = meta_y;

		for (int x = 0; x < dst_width; ++x)
		{
			const T *src_ptr1 = src1 + (meta->start_y * src_pitch1 + meta->start_x);
			const T *src_ptr2 = src2 + (meta->start_y * src_pitch2 + meta->start_x);
			const T *src_ptr3 = src3 + (meta->start_y * src_pitch3 + meta->start_x);
			const T *src_ptr4 = src4 + (meta->start_y * src_pitch4 + meta->start_x);
			const float *coeff_ptr = tab_coeff->factor + meta->coeff_meta;
			__m512 result1 = _mm512_setzero_ps();
			__m512 result2 = _mm512_setzero_ps();
			__m512 result3 = _mm512_setzero_ps();
			__m512 result4 = _mm512_setzero_ps();
			__m512 result1_2 = _mm512_setzero_ps();
			__m512 result2_2 = _mm512_setzero_ps();
			__m512 result3_2 = _mm512_setzero_ps();
			__m512 result4_2 = _mm512_setzero_ps();

			if JincMT_CONSTEXPR (std::is_same<T, uint8_t>::value)
			{
				for (int ly = 0; ly < filter_size_mod2; ly += 2)
				{
					for (int lx = 0; lx < filter_size; lx += 16)
					{
						const __m512 src_ps1 = _mm512_cvtepi32_ps(_mm512_cvtepu8_epi32(_mm_loadu_si128(reinterpret_cast<const __m128i*>(src_ptr1 + lx))));
						const __m512 src_ps1_2 = _mm512_cvtepi32_ps(_mm512_cvtepu8_epi32(_mm_loadu_si128(reinterpret_cast<const __m128i*>(src_ptr1 + lx + src_pitch1))));
						const __m512 src_ps2 = _mm512_cvtepi32_ps(_mm512_cvtepu8_epi32(_mm_loadu_si128(reinterpret_cast<const __m128i*>(src_ptr2 + lx))));
						const __m512 src_ps2_2 = _mm512_cvtepi32_ps(_mm512_cvtepu8_epi32(_mm_loadu_si128(reinterpret_cast<const __m128i*>(src_ptr2 + lx + src_pitch2))));
						const __m512 src_ps3 = _mm512_cvtepi32_ps(_mm512_cvtepu8_epi32(_mm_loadu_si128(reinterpret_cast<const __m128i*>(src_ptr3 + lx))));
						const __m512 src_ps3_2 = _mm512_cvtepi32_ps(_mm512_cvtepu8_epi32(_mm_loadu_si128(reinterpret_cast<const __m128i*>(src_ptr3 + lx + src_pitch3))));
						const __m512 src_ps4 = _mm512_cvtepi32_ps(_mm512_cvtepu8_epi32(_mm_loadu_si128(reinterpret_cast<const __m128i*>(src_ptr4 + lx))));
						const __m512 src_ps4_2 = _mm512_cvtepi32_ps(_mm512_cvtepu8_epi32(_mm_loadu_si128(reinterpret_cast<const __m128i*>(src_ptr4 + lx + src_pitch4))));
						__m512 coeff, coeff2;

						if JincMT_CONSTEXPR (bFP16)
						{
							const int lx2 = lx >> 1;
							const __m256i coeff_fp16 = _mm256_load_si256((__m256i*)(coeff_ptr + lx2));
							const __m256i coeff2_fp16 = _mm256_load_si256((__m256i*)(coeff_ptr + lx2 + coeff_stride));

							coeff = _mm512_cvtph_ps(coeff_fp16); // separated grouped converts may be better for instructions grouping at some compilers
							coeff2 = _mm512_cvtph_ps(coeff2_fp16);
						}
						else
						{
							coeff = _mm512_load_ps(coeff_ptr + lx);
							coeff2 = _mm512_load_ps(coeff_ptr + lx + coeff_stride);
						}

						result1 = _mm512_fmadd_ps(src_ps1, coeff, result1);
						result1_2 = _mm512_fmadd_ps(src_ps1_2, coeff2, result1_2);
						result2 = _mm512_fmadd_ps(src_ps2, coeff, result2);
						result2_2 = _mm512_fmadd_ps(src_ps2_2, coeff2, result2_2);
						result3 = _mm512_fmadd_ps(src_ps3, coeff, result3);
						result3_2 = _mm512_fmadd_ps(src_ps3_2, coeff2, result3_2);
						result4 = _mm512_fmadd_ps(src_ps4, coeff, result4);
						result4_2 = _mm512_fmadd_ps(src_ps4_2, coeff2, result4_2);
					}

					coeff_ptr += coeff_stridex2;
					src_ptr1 += src_pitch1x2;
					src_ptr2 += src_pitch2x2;
					src_ptr3 += src_pitch3x2;
					src_ptr4 += src_pitch4x2;
				}

				result1 = _mm512_add_ps(result1, result1_2);
				result2 = _mm512_add_ps(result2, result2_2);
				result3 = _mm512_add_ps(result3, result3_2);
				result4 = _mm512_add_ps(result4, result4_2);

				// last row
				if (fs_notMod2)
				{
					for (int lx = 0; lx < filter_size; lx += 16)
					{
						const __m512 src_ps1 = _mm512_cvtepi32_ps(_mm512_cvtepu8_epi32(_mm_loadu_si128(reinterpret_cast<const __m128i*>(src_ptr1 + lx))));
						const __m512 src_ps2 = _mm512_cvtepi32_ps(_mm512_cvtepu8_epi32(_mm_loadu_si128(reinterpret_cast<const __m128i*>(src_ptr2 + lx))));
						const __m512 src_ps3 = _mm512_cvtepi32_ps(_mm512_cvtepu8_epi32(_mm_loadu_si128(reinterpret_cast<const __m128i*>(src_ptr3 + lx))));
						const __m512 src_ps4 = _mm512_cvtepi32_ps(_mm512_cvtepu8_epi32(_mm_loadu_si128(reinterpret_cast<const __m128i*>(src_ptr4 + lx))));
						__m512 coeff;

						if JincMT_CONSTEXPR (bFP16)
							coeff = _mm512_cvtph_ps(_mm256_load_si256((__m256i*)(coeff_ptr + (lx >> 1))));
						else
							coeff = _mm512_load_ps(coeff_ptr + lx);

						result1 = _mm512_fmadd_ps(src_ps1, coeff, result1);
						result2 = _mm512_fmadd_ps(src_ps2, coeff, result2);
						result3 = _mm512_fmadd_ps(src_ps3, coeff, result3);
						result4 = _mm512_fmadd_ps(src_ps4, coeff, result4);
					}
				}

				const __m256 lo_hi_256_1 = _mm256_add_ps(_mm512_castps512_ps256(result1), _mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(result1), 1)));
				const __m256 lo_hi_256_2 = _mm256_add_ps(_mm512_castps512_ps256(result2), _mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(result2), 1)));
				const __m256 lo_hi_256_3 = _mm256_add_ps(_mm512_castps512_ps256(result3), _mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(result3), 1)));
				const __m256 lo_hi_256_4 = _mm256_add_ps(_mm512_castps512_ps256(result4), _mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(result4), 1)));
				__m128 hsum1 = _mm_add_ps(_mm256_castps256_ps128(lo_hi_256_1), _mm256_extractf128_ps(lo_hi_256_1, 1));
				__m128 hsum2 = _mm_add_ps(_mm256_castps256_ps128(lo_hi_256_2), _mm256_extractf128_ps(lo_hi_256_2, 1));
				__m128 hsum3 = _mm_add_ps(_mm256_castps256_ps128(lo_hi_256_3), _mm256_extractf128_ps(lo_hi_256_3, 1));
				__m128 hsum4 = _mm_add_ps(_mm256_castps256_ps128(lo_hi_256_4), _mm256_extractf128_ps(lo_hi_256_4, 1));
				hsum1 = _mm_hadd_ps(_mm_hadd_ps(hsum1, hsum1), _mm_hadd_ps(hsum1, hsum1));
				hsum2 = _mm_hadd_ps(_mm_hadd_ps(hsum2, hsum2), _mm_hadd_ps(hsum2, hsum2));
				hsum3 = _mm_hadd_ps(_mm_hadd_ps(hsum3, hsum3), _mm_hadd_ps(hsum3, hsum3));
				hsum4 = _mm_hadd_ps(_mm_hadd_ps(hsum4, hsum4), _mm_hadd_ps(hsum4, hsum4));
				const T final_res1 = _mm_cvtsi128_si32(_mm_packus_epi16(_mm_packus_epi32(_mm_cvtps_epi32(hsum1), _mm_setzero_si128()), _mm_setzero_si128()));
				const T final_res2 = _mm_cvtsi128_si32(_mm_packus_epi16(_mm_packus_epi32(_mm_cvtps_epi32(hsum2), _mm_setzero_si128()), _mm_setzero_si128()));
				const T final_res3 = _mm_cvtsi128_si32(_mm_packus_epi16(_mm_packus_epi32(_mm_cvtps_epi32(hsum3), _mm_setzero_si128()), _mm_setzero_si128()));
				const T final_res4 = _mm_cvtsi128_si32(_mm_packus_epi16(_mm_packus_epi32(_mm_cvtps_epi32(hsum4), _mm_setzero_si128()), _mm_setzero_si128()));

				dst1[x] = clamp(final_res1, val_min1, val_max1);
				dst2[x] = clamp(final_res2, val_min2, val_max2);
				dst3[x] = clamp(final_res3, val_min3, val_max3);
				dst4[x] = clamp(final_res4, val_min4, val_max4);
			}
			else if JincMT_CONSTEXPR (std::is_same<T, uint16_t>::value)
			{
				for (int ly = 0; ly < filter_size_mod2; ly += 2)
				{
					for (int lx = 0; lx < filter_size; lx += 16)
					{
						const __m512 src_ps1 = _mm512_cvtepi32_ps(_mm512_cvtepu16_epi32(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(src_ptr1 + lx))));
						const __m512 src_ps1_2 = _mm512_cvtepi32_ps(_mm512_cvtepu16_epi32(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(src_ptr1 + lx + src_pitch1))));
						const __m512 src_ps2 = _mm512_cvtepi32_ps(_mm512_cvtepu16_epi32(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(src_ptr2 + lx))));
						const __m512 src_ps2_2 = _mm512_cvtepi32_ps(_mm512_cvtepu16_epi32(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(src_ptr2 + lx + src_pitch2))));
						const __m512 src_ps3 = _mm512_cvtepi32_ps(_mm512_cvtepu16_epi32(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(src_ptr3 + lx))));
						const __m512 src_ps3_2 = _mm512_cvtepi32_ps(_mm512_cvtepu16_epi32(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(src_ptr3 + lx + src_pitch3))));
						const __m512 src_ps4 = _mm512_cvtepi32_ps(_mm512_cvtepu16_epi32(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(src_ptr4 + lx))));
						const __m512 src_ps4_2 = _mm512_cvtepi32_ps(_mm512_cvtepu16_epi32(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(src_ptr4 + lx + src_pitch4))));
						__m512 coeff, coeff2;

						if JincMT_CONSTEXPR (bFP16)
						{
							const int lx2 = lx >> 1;
							const __m256i coeff_fp16 = _mm256_load_si256((__m256i*)(coeff_ptr + lx2));
							const __m256i coeff2_fp16 = _mm256_load_si256((__m256i*)(coeff_ptr + lx2 + coeff_stride));

							coeff = _mm512_cvtph_ps(coeff_fp16); // separated grouped converts may be better for instructions grouping at some compilers
							coeff2 = _mm512_cvtph_ps(coeff2_fp16);
						}
						else
						{
							coeff = _mm512_load_ps(coeff_ptr + lx);
							coeff2 = _mm512_load_ps(coeff_ptr + lx + coeff_stride);
						}

						result1 = _mm512_fmadd_ps(src_ps1, coeff, result1);
						result1_2 = _mm512_fmadd_ps(src_ps1_2, coeff2, result1_2);
						result2 = _mm512_fmadd_ps(src_ps2, coeff, result2);
						result2_2 = _mm512_fmadd_ps(src_ps2_2, coeff2, result2_2);
						result3 = _mm512_fmadd_ps(src_ps3, coeff, result3);
						result3_2 = _mm512_fmadd_ps(src_ps3_2, coeff2, result3_2);
						result4 = _mm512_fmadd_ps(src_ps4, coeff, result4);
						result4_2 = _mm512_fmadd_ps(src_ps4_2, coeff2, result4_2);
					}

					coeff_ptr += coeff_stridex2;
					src_ptr1 += src_pitch1x2;
					src_ptr2 += src_pitch2x2;
					src_ptr3 += src_pitch3x2;
					src_ptr4 += src_pitch4x2;
				}

				result1 = _mm512_add_ps(result1, result1_2);
				result2 = _mm512_add_ps(result2, result2_2);
				result3 = _mm512_add_ps(result3, result3_2);
				result4 = _mm512_add_ps(result4, result4_2);

				// last row
				if (fs_notMod2)
				{
					for (int lx = 0; lx < filter_size; lx += 16)
					{
						const __m512 src_ps1 = _mm512_cvtepi32_ps(_mm512_cvtepu16_epi32(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(src_ptr1 + lx))));
						const __m512 src_ps2 = _mm512_cvtepi32_ps(_mm512_cvtepu16_epi32(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(src_ptr2 + lx))));
						const __m512 src_ps3 = _mm512_cvtepi32_ps(_mm512_cvtepu16_epi32(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(src_ptr3 + lx))));
						const __m512 src_ps4 = _mm512_cvtepi32_ps(_mm512_cvtepu16_epi32(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(src_ptr4 + lx))));
						__m512 coeff;

						if JincMT_CONSTEXPR (bFP16)
							coeff = _mm512_cvtph_ps(_mm256_load_si256((__m256i*)(coeff_ptr + (lx >> 1))));
						else
							coeff = _mm512_load_ps(coeff_ptr + lx);

						result1 = _mm512_fmadd_ps(src_ps1, coeff, result1);
						result2 = _mm512_fmadd_ps(src_ps2, coeff, result2);
						result3 = _mm512_fmadd_ps(src_ps3, coeff, result3);
						result4 = _mm512_fmadd_ps(src_ps4, coeff, result4);
					}
				}

				const __m256 lo_hi_256_1 = _mm256_add_ps(_mm512_castps512_ps256(result1), _mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(result1), 1)));
				const __m256 lo_hi_256_2 = _mm256_add_ps(_mm512_castps512_ps256(result2), _mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(result2), 1)));
				const __m256 lo_hi_256_3 = _mm256_add_ps(_mm512_castps512_ps256(result3), _mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(result3), 1)));
				const __m256 lo_hi_256_4 = _mm256_add_ps(_mm512_castps512_ps256(result4), _mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(result4), 1)));
				__m128 hsum1 = _mm_add_ps(_mm256_castps256_ps128(lo_hi_256_1), _mm256_extractf128_ps(lo_hi_256_1, 1));
				__m128 hsum2 = _mm_add_ps(_mm256_castps256_ps128(lo_hi_256_2), _mm256_extractf128_ps(lo_hi_256_2, 1));
				__m128 hsum3 = _mm_add_ps(_mm256_castps256_ps128(lo_hi_256_3), _mm256_extractf128_ps(lo_hi_256_3, 1));
				__m128 hsum4 = _mm_add_ps(_mm256_castps256_ps128(lo_hi_256_4), _mm256_extractf128_ps(lo_hi_256_4, 1));
				hsum1 = _mm_hadd_ps(_mm_hadd_ps(hsum1, hsum1), _mm_hadd_ps(hsum1, hsum1));
				hsum2 = _mm_hadd_ps(_mm_hadd_ps(hsum2, hsum2), _mm_hadd_ps(hsum2, hsum2));
				hsum3 = _mm_hadd_ps(_mm_hadd_ps(hsum3, hsum3), _mm_hadd_ps(hsum3, hsum3));
				hsum4 = _mm_hadd_ps(_mm_hadd_ps(hsum4, hsum4), _mm_hadd_ps(hsum4, hsum4));
				const T final_res1 = _mm_cvtsi128_si32(_mm_packus_epi32(_mm_cvtps_epi32(hsum1), _mm_setzero_si128()));
				const T final_res2 = _mm_cvtsi128_si32(_mm_packus_epi32(_mm_cvtps_epi32(hsum2), _mm_setzero_si128()));
				const T final_res3 = _mm_cvtsi128_si32(_mm_packus_epi32(_mm_cvtps_epi32(hsum3), _mm_setzero_si128()));
				const T final_res4 = _mm_cvtsi128_si32(_mm_packus_epi32(_mm_cvtps_epi32(hsum4), _mm_setzero_si128()));

				dst1[x] = clamp(final_res1, val_min1, val_max1);
				dst2[x] = clamp(final_res2, val_min2, val_max2);
				dst3[x] = clamp(final_res3, val_min3, val_max3);
				dst4[x] = clamp(final_res4, val_min4, val_max4);
			}
			else
			{
				for (int ly = 0; ly < filter_size_mod2; ly += 2)
				{
					for (int lx = 0; lx < filter_size; lx += 16)
					{
						const __m512 src_ps1 = _mm512_max_ps(_mm512_loadu_ps(src_ptr1 + lx), min_val1);
						const __m512 src_ps1_2 = _mm512_max_ps(_mm512_loadu_ps(src_ptr1 + lx + src_pitch1), min_val1);
						const __m512 src_ps2 = _mm512_max_ps(_mm512_loadu_ps(src_ptr2 + lx), min_val2);
						const __m512 src_ps2_2 = _mm512_max_ps(_mm512_loadu_ps(src_ptr2 + lx + src_pitch2), min_val2);
						const __m512 src_ps3 = _mm512_max_ps(_mm512_loadu_ps(src_ptr3 + lx), min_val3);
						const __m512 src_ps3_2 = _mm512_max_ps(_mm512_loadu_ps(src_ptr3 + lx + src_pitch3), min_val3);
						const __m512 src_ps4 = _mm512_max_ps(_mm512_loadu_ps(src_ptr4 + lx), min_val4);
						const __m512 src_ps4_2 = _mm512_max_ps(_mm512_loadu_ps(src_ptr4 + lx + src_pitch4), min_val4);
						__m512 coeff, coeff2;

						if JincMT_CONSTEXPR (bFP16)
						{
							const int lx2 = lx >> 1;
							const __m256i coeff_fp16 = _mm256_load_si256((__m256i*)(coeff_ptr + lx2));
							const __m256i coeff2_fp16 = _mm256_load_si256((__m256i*)(coeff_ptr + lx2 + coeff_stride));

							coeff = _mm512_cvtph_ps(coeff_fp16); // separated grouped converts may be better for instructions grouping at some compilers
							coeff2 = _mm512_cvtph_ps(coeff2_fp16);
						}
						else
						{
							coeff = _mm512_load_ps(coeff_ptr + lx);
							coeff2 = _mm512_load_ps(coeff_ptr + lx + coeff_stride);
						}

						result1 = _mm512_fmadd_ps(src_ps1, coeff, result1);
						result1_2 = _mm512_fmadd_ps(src_ps1_2, coeff2, result1_2);
						result2 = _mm512_fmadd_ps(src_ps2, coeff, result2);
						result2_2 = _mm512_fmadd_ps(src_ps2_2, coeff2, result2_2);
						result3 = _mm512_fmadd_ps(src_ps3, coeff, result3);
						result3_2 = _mm512_fmadd_ps(src_ps3_2, coeff2, result3_2);
						result4 = _mm512_fmadd_ps(src_ps4, coeff, result4);
						result4_2 = _mm512_fmadd_ps(src_ps4_2, coeff2, result4_2);
					}

					coeff_ptr += coeff_stridex2;
					src_ptr1 += src_pitch1x2;
					src_ptr2 += src_pitch2x2;
					src_ptr3 += src_pitch3x2;
					src_ptr4 += src_pitch4x2;
				}

				result1 = _mm512_add_ps(result1, result1_2);
				result2 = _mm512_add_ps(result2, result2_2);
				result3 = _mm512_add_ps(result3, result3_2);
				result4 = _mm512_add_ps(result4, result4_2);

				// last row
				if (fs_notMod2)
				{
					for (int lx = 0; lx < filter_size; lx += 16)
					{
						const __m512 src_ps1 = _mm512_max_ps(_mm512_loadu_ps(src_ptr1 + lx), min_val1);
						const __m512 src_ps2 = _mm512_max_ps(_mm512_loadu_ps(src_ptr2 + lx), min_val2);
						const __m512 src_ps3 = _mm512_max_ps(_mm512_loadu_ps(src_ptr3 + lx), min_val3);
						const __m512 src_ps4 = _mm512_max_ps(_mm512_loadu_ps(src_ptr4 + lx), min_val4);
						__m512 coeff;

						if JincMT_CONSTEXPR (bFP16)
							coeff = _mm512_cvtph_ps(_mm256_load_si256((__m256i*)(coeff_ptr + (lx >> 1))));
						else
							coeff = _mm512_load_ps(coeff_ptr + lx);

						result1 = _mm512_fmadd_ps(src_ps1, coeff, result1);
						result2 = _mm512_fmadd_ps(src_ps2, coeff, result2);
						result3 = _mm512_fmadd_ps(src_ps3, coeff, result3);
						result4 = _mm512_fmadd_ps(src_ps4, coeff, result4);
					}
				}

				const __m256 lo_hi_256_1 = _mm256_add_ps(_mm512_castps512_ps256(result1), _mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(result1), 1)));
				const __m256 lo_hi_256_2 = _mm256_add_ps(_mm512_castps512_ps256(result2), _mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(result2), 1)));
				const __m256 lo_hi_256_3 = _mm256_add_ps(_mm512_castps512_ps256(result3), _mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(result3), 1)));
				const __m256 lo_hi_256_4 = _mm256_add_ps(_mm512_castps512_ps256(result4), _mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(result4), 1)));
				__m128 hsum1 = _mm_add_ps(_mm256_castps256_ps128(lo_hi_256_1), _mm256_extractf128_ps(lo_hi_256_1, 1));
				__m128 hsum2 = _mm_add_ps(_mm256_castps256_ps128(lo_hi_256_2), _mm256_extractf128_ps(lo_hi_256_2, 1));
				__m128 hsum3 = _mm_add_ps(_mm256_castps256_ps128(lo_hi_256_3), _mm256_extractf128_ps(lo_hi_256_3, 1));
				__m128 hsum4 = _mm_add_ps(_mm256_castps256_ps128(lo_hi_256_4), _mm256_extractf128_ps(lo_hi_256_4, 1));

				dst1[x] = _mm_cvtss_f32(_mm_hadd_ps(_mm_hadd_ps(hsum1, hsum1), _mm_hadd_ps(hsum1, hsum1)));
				dst2[x] = _mm_cvtss_f32(_mm_hadd_ps(_mm_hadd_ps(hsum2, hsum2), _mm_hadd_ps(hsum2, hsum2)));
				dst3[x] = _mm_cvtss_f32(_mm_hadd_ps(_mm_hadd_ps(hsum3, hsum3), _mm_hadd_ps(hsum3, hsum3)));
				dst4[x] = _mm_cvtss_f32(_mm_hadd_ps(_mm_hadd_ps(hsum4, hsum4), _mm_hadd_ps(hsum4, hsum4)));
			}
			meta++;
		} // for (x)
		meta_y += dst_width;
		dst1 += dst_pitch1;
		dst2 += dst_pitch2;
		dst3 += dst_pitch3;
		dst4 += dst_pitch4;
	} // for (y)
}


template <typename T, bool bFP16>
#if defined(CLANG)
__attribute__((__target__("avx512f")))
#endif
void resize_plane_avx2_4x_v2(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *tab_coeff,
	const float Val_Min[], const float Val_Max[])
{
	const uint8_t idx1 = 0;
	const uint8_t idx2 = 1;
	const uint8_t idx3 = 2;
	const uint8_t idx4 = 3;

	const T *src1 = reinterpret_cast<const T*>(MT_DataGF->src[idx1]);
	const T *src2 = reinterpret_cast<const T*>(MT_DataGF->src[idx2]);
	const T *src3 = reinterpret_cast<const T*>(MT_DataGF->src[idx3]);
	const T *src4 = reinterpret_cast<const T*>(MT_DataGF->src[idx4]);
	T *JincMT_RESTRICT dst1 = reinterpret_cast<T*>(MT_DataGF->dst[idx1]);
	T *JincMT_RESTRICT dst2 = reinterpret_cast<T*>(MT_DataGF->dst[idx2]);
	T *JincMT_RESTRICT dst3 = reinterpret_cast<T*>(MT_DataGF->dst[idx3]);
	T *JincMT_RESTRICT dst4 = reinterpret_cast<T*>(MT_DataGF->dst[idx4]);

	const ptrdiff_t src_pitch1 = (ptrdiff_t)MT_DataGF->src_pitch[idx1] / sizeof(T);
	const ptrdiff_t src_pitch2 = (ptrdiff_t)MT_DataGF->src_pitch[idx2] / sizeof(T);
	const ptrdiff_t src_pitch3 = (ptrdiff_t)MT_DataGF->src_pitch[idx3] / sizeof(T);
	const ptrdiff_t src_pitch4 = (ptrdiff_t)MT_DataGF->src_pitch[idx4] / sizeof(T);
	const ptrdiff_t dst_pitch1 = (ptrdiff_t)MT_DataGF->dst_pitch[idx1] / sizeof(T);
	const ptrdiff_t dst_pitch2 = (ptrdiff_t)MT_DataGF->dst_pitch[idx2] / sizeof(T);
	const ptrdiff_t dst_pitch3 = (ptrdiff_t)MT_DataGF->dst_pitch[idx3] / sizeof(T);
	const ptrdiff_t dst_pitch4 = (ptrdiff_t)MT_DataGF->dst_pitch[idx4] / sizeof(T);
	const ptrdiff_t src_pitch1x2 = src_pitch1 << 1;
	const ptrdiff_t src_pitch2x2 = src_pitch2 << 1;
	const ptrdiff_t src_pitch3x2 = src_pitch3 << 1;
	const ptrdiff_t src_pitch4x2 = src_pitch4 << 1;

	const int Y_Min = MT_DataGF->dst_Y_h_min;
	const int Y_Max = MT_DataGF->dst_Y_h_max;
	const int dst_width = MT_DataGF->dst_Y_w;

	const T val_min1 = static_cast<T>(Val_Min[idx1]);
	const T val_min2 = static_cast<T>(Val_Min[idx2]);
	const T val_min3 = static_cast<T>(Val_Min[idx3]);
	const T val_min4 = static_cast<T>(Val_Min[idx4]);
	const T val_max1 = static_cast<T>(Val_Max[idx1]);
	const T val_max2 = static_cast<T>(Val_Max[idx2]);
	const T val_max3 = static_cast<T>(Val_Max[idx3]);
	const T val_max4 = static_cast<T>(Val_Max[idx4]);
	const __m256 min_val1 = _mm256_set1_ps(Val_Min[idx1]);
	const __m256 min_val2 = _mm256_set1_ps(Val_Min[idx2]);
	const __m256 min_val3 = _mm256_set1_ps(Val_Min[idx3]);
	const __m256 min_val4 = _mm256_set1_ps(Val_Min[idx4]);

	EWAPixelCoeffMeta *meta_y = tab_coeff->meta + (Y_Min*dst_width);

	const int filter_size = tab_coeff->filter_size, coeff_stride = tab_coeff->coeff_stride;
	const int coeff_stridex2 = coeff_stride << 1;

	const int filter_size_mod2 = (filter_size >> 1) << 1;
	const bool fs_notMod2 = filter_size_mod2 < filter_size;

	for (int y = Y_Min; y < Y_Max; y++)
	{
		EWAPixelCoeffMeta *meta = meta_y;

		for (int x = 0; x < dst_width; ++x)
		{
			const T *src_ptr1 = src1 + (meta->start_y * src_pitch1 + meta->start_x);
			const T *src_ptr2 = src2 + (meta->start_y * src_pitch2 + meta->start_x);
			const T *src_ptr3 = src3 + (meta->start_y * src_pitch3 + meta->start_x);
			const T *src_ptr4 = src4 + (meta->start_y * src_pitch4 + meta->start_x);
			const float *coeff_ptr = tab_coeff->factor + meta->coeff_meta;
			__m256 result1 = _mm256_setzero_ps();
			__m256 result2 = _mm256_setzero_ps();
			__m256 result3 = _mm256_setzero_ps();
			__m256 result4 = _mm256_setzero_ps();
			__m256 result1_2 = _mm256_setzero_ps();
			__m256 result2_2 = _mm256_setzero_ps();
			__m256 result3_2 = _mm256_setzero_ps();
			__m256 result4_2 = _mm256_setzero_ps();

			if JincMT_CONSTEXPR (std::is_same<T, uint8_t>::value)
			{
				for (int ly = 0; ly < filter_size_mod2; ly += 2)
				{
					for (int lx = 0; lx < filter_size; lx += 8)
					{
						const __m256 src_ps1 = _mm256_cvtepi32_ps(_mm256_cvtepu8_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr1 + lx)))));
						const __m256 src_ps1_2 = _mm256_cvtepi32_ps(_mm256_cvtepu8_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr1 + lx + src_pitch1)))));
						const __m256 src_ps2 = _mm256_cvtepi32_ps(_mm256_cvtepu8_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr2 + lx)))));
						const __m256 src_ps2_2 = _mm256_cvtepi32_ps(_mm256_cvtepu8_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr2 + lx + src_pitch2)))));
						const __m256 src_ps3 = _mm256_cvtepi32_ps(_mm256_cvtepu8_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr3 + lx)))));
						const __m256 src_ps3_2 = _mm256_cvtepi32_ps(_mm256_cvtepu8_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr3 + lx + src_pitch3)))));
						const __m256 src_ps4 = _mm256_cvtepi32_ps(_mm256_cvtepu8_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr4 + lx)))));
						const __m256 src_ps4_2 = _mm256_cvtepi32_ps(_mm256_cvtepu8_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr4 + lx + src_pitch4)))));
						__m256 coeff, coeff2;

						if JincMT_CONSTEXPR (bFP16)
						{
							const int lx2 = lx >> 1;
							__m128i coeff_fp16 = _mm_load_si128((__m128i*)(coeff_ptr + lx2));
							__m128i coeff2_fp16 = _mm_load_si128((__m128i*)(coeff_ptr + lx2 + coeff_stride));

							coeff = _mm256_cvtph_ps(coeff_fp16); // separated grouped converts may be better for instructions grouping at some compilers
							coeff2 = _mm256_cvtph_ps(coeff2_fp16);
						}
						else
						{
							coeff = _mm256_load_ps(coeff_ptr + lx);
							coeff2 = _mm256_load_ps(coeff_ptr + lx + coeff_stride);
						}

						result1 = _mm256_fmadd_ps(src_ps1, coeff, result1);
						result1_2 = _mm256_fmadd_ps(src_ps1_2, coeff2, result1_2);
						result2 = _mm256_fmadd_ps(src_ps2, coeff, result2);
						result2_2 = _mm256_fmadd_ps(src_ps2_2, coeff2, result2_2);
						result3 = _mm256_fmadd_ps(src_ps3, coeff, result3);
						result3_2 = _mm256_fmadd_ps(src_ps3_2, coeff2, result3_2);
						result4 = _mm256_fmadd_ps(src_ps4, coeff, result4);
						result4_2 = _mm256_fmadd_ps(src_ps4_2, coeff2, result4_2);
					}

					coeff_ptr += coeff_stridex2;
					src_ptr1 += src_pitch1x2;
					src_ptr2 += src_pitch2x2;
					src_ptr3 += src_pitch3x2;
					src_ptr4 += src_pitch4x2;
				}

				result1 = _mm256_add_ps(result1, result1_2);
				result2 = _mm256_add_ps(result2, result2_2);
				result3 = _mm256_add_ps(result3, result3_2);
				result4 = _mm256_add_ps(result4, result4_2);

				// last row
				if (fs_notMod2)
				{
					for (int lx = 0; lx < filter_size; lx += 8)
					{
						const __m256 src_ps1 = _mm256_cvtepi32_ps(_mm256_cvtepu8_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr1 + lx)))));
						const __m256 src_ps2 = _mm256_cvtepi32_ps(_mm256_cvtepu8_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr2 + lx)))));
						const __m256 src_ps3 = _mm256_cvtepi32_ps(_mm256_cvtepu8_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr3 + lx)))));
						const __m256 src_ps4 = _mm256_cvtepi32_ps(_mm256_cvtepu8_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr4 + lx)))));
						__m256 coeff;

						if JincMT_CONSTEXPR (bFP16)
							coeff = _mm256_cvtph_ps(_mm_load_si128((__m128i*)(coeff_ptr + (lx >> 1))));
						else
							coeff = _mm256_load_ps(coeff_ptr + lx);

						result1 = _mm256_fmadd_ps(src_ps1, coeff, result1);
						result2 = _mm256_fmadd_ps(src_ps2, coeff, result2);
						result3 = _mm256_fmadd_ps(src_ps3, coeff, result3);
						result4 = _mm256_fmadd_ps(src_ps4, coeff, result4);
					}
				}

				__m128 hsum1 = _mm_add_ps(_mm256_castps256_ps128(result1), _mm256_extractf128_ps(result1, 1));
				__m128 hsum2 = _mm_add_ps(_mm256_castps256_ps128(result2), _mm256_extractf128_ps(result2, 1));
				__m128 hsum3 = _mm_add_ps(_mm256_castps256_ps128(result3), _mm256_extractf128_ps(result3, 1));
				__m128 hsum4 = _mm_add_ps(_mm256_castps256_ps128(result4), _mm256_extractf128_ps(result4, 1));
				hsum1 = _mm_hadd_ps(_mm_hadd_ps(hsum1, hsum1), _mm_hadd_ps(hsum1, hsum1));
				hsum2 = _mm_hadd_ps(_mm_hadd_ps(hsum2, hsum2), _mm_hadd_ps(hsum2, hsum2));
				hsum3 = _mm_hadd_ps(_mm_hadd_ps(hsum3, hsum3), _mm_hadd_ps(hsum3, hsum3));
				hsum4 = _mm_hadd_ps(_mm_hadd_ps(hsum4, hsum4), _mm_hadd_ps(hsum4, hsum4));
				const T final_res1 = _mm_cvtsi128_si32(_mm_packus_epi16(_mm_packus_epi32(_mm_cvtps_epi32(hsum1), _mm_setzero_si128()), _mm_setzero_si128()));
				const T final_res2 = _mm_cvtsi128_si32(_mm_packus_epi16(_mm_packus_epi32(_mm_cvtps_epi32(hsum2), _mm_setzero_si128()), _mm_setzero_si128()));
				const T final_res3 = _mm_cvtsi128_si32(_mm_packus_epi16(_mm_packus_epi32(_mm_cvtps_epi32(hsum3), _mm_setzero_si128()), _mm_setzero_si128()));
				const T final_res4 = _mm_cvtsi128_si32(_mm_packus_epi16(_mm_packus_epi32(_mm_cvtps_epi32(hsum4), _mm_setzero_si128()), _mm_setzero_si128()));

				dst1[x] = clamp(final_res1, val_min1, val_max1);
				dst2[x] = clamp(final_res2, val_min2, val_max2);
				dst3[x] = clamp(final_res3, val_min3, val_max3);
				dst4[x] = clamp(final_res4, val_min4, val_max4);
			}
			else if JincMT_CONSTEXPR (std::is_same<T, uint16_t>::value)
			{
				for (int ly = 0; ly < filter_size_mod2; ly += 2)
				{
					for (int lx = 0; lx < filter_size; lx += 8)
					{
						const __m256 src_ps1 = _mm256_cvtepi32_ps(_mm256_cvtepu16_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr1 + lx)))));
						const __m256 src_ps1_2 = _mm256_cvtepi32_ps(_mm256_cvtepu16_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr1 + lx + src_pitch1)))));
						const __m256 src_ps2 = _mm256_cvtepi32_ps(_mm256_cvtepu16_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr2 + lx)))));
						const __m256 src_ps2_2 = _mm256_cvtepi32_ps(_mm256_cvtepu16_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr2 + lx + src_pitch2)))));
						const __m256 src_ps3 = _mm256_cvtepi32_ps(_mm256_cvtepu16_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr3 + lx)))));
						const __m256 src_ps3_2 = _mm256_cvtepi32_ps(_mm256_cvtepu16_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr3 + lx + src_pitch3)))));
						const __m256 src_ps4 = _mm256_cvtepi32_ps(_mm256_cvtepu16_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr4 + lx)))));
						const __m256 src_ps4_2 = _mm256_cvtepi32_ps(_mm256_cvtepu16_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr4 + lx + src_pitch4)))));
						__m256 coeff, coeff2;

						if JincMT_CONSTEXPR (bFP16)
						{
							const int lx2 = lx >> 1;
							__m128i coeff_fp16 = _mm_load_si128((__m128i*)(coeff_ptr + lx2));
							__m128i coeff2_fp16 = _mm_load_si128((__m128i*)(coeff_ptr + lx2 + coeff_stride));

							coeff = _mm256_cvtph_ps(coeff_fp16); // separated grouped converts may be better for instructions grouping at some compilers
							coeff2 = _mm256_cvtph_ps(coeff2_fp16);
						}
						else
						{
							coeff = _mm256_load_ps(coeff_ptr + lx);
							coeff2 = _mm256_load_ps(coeff_ptr + lx + coeff_stride);
						}

						result1 = _mm256_fmadd_ps(src_ps1, coeff, result1);
						result1_2 = _mm256_fmadd_ps(src_ps1_2, coeff2, result1_2);
						result2 = _mm256_fmadd_ps(src_ps2, coeff, result2);
						result2_2 = _mm256_fmadd_ps(src_ps2_2, coeff2, result2_2);
						result3 = _mm256_fmadd_ps(src_ps3, coeff, result3);
						result3_2 = _mm256_fmadd_ps(src_ps3_2, coeff2, result3_2);
						result4 = _mm256_fmadd_ps(src_ps4, coeff, result4);
						result4_2 = _mm256_fmadd_ps(src_ps4_2, coeff2, result4_2);
					}

					coeff_ptr += coeff_stridex2;
					src_ptr1 += src_pitch1x2;
					src_ptr2 += src_pitch2x2;
					src_ptr3 += src_pitch3x2;
					src_ptr4 += src_pitch4x2;
				}

				result1 = _mm256_add_ps(result1, result1_2);
				result2 = _mm256_add_ps(result2, result2_2);
				result3 = _mm256_add_ps(result3, result3_2);
				result4 = _mm256_add_ps(result4, result4_2);

				// last row
				if (fs_notMod2)
				{
					for (int lx = 0; lx < filter_size; lx += 8)
					{
						const __m256 src_ps1 = _mm256_cvtepi32_ps(_mm256_cvtepu16_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr1 + lx)))));
						const __m256 src_ps2 = _mm256_cvtepi32_ps(_mm256_cvtepu16_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr2 + lx)))));
						const __m256 src_ps3 = _mm256_cvtepi32_ps(_mm256_cvtepu16_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr3 + lx)))));
						const __m256 src_ps4 = _mm256_cvtepi32_ps(_mm256_cvtepu16_epi32(_mm_loadu_si128(const_cast<__m128i*>(reinterpret_cast<const __m128i*>(src_ptr4 + lx)))));
						__m256 coeff;

						if JincMT_CONSTEXPR (bFP16)
							coeff = _mm256_cvtph_ps(_mm_load_si128((__m128i*)(coeff_ptr + (lx >> 1))));
						else
							coeff = _mm256_load_ps(coeff_ptr + lx);

						result1 = _mm256_fmadd_ps(src_ps1, coeff, result1);
						result2 = _mm256_fmadd_ps(src_ps2, coeff, result2);
						result3 = _mm256_fmadd_ps(src_ps3, coeff, result3);
						result4 = _mm256_fmadd_ps(src_ps4, coeff, result4);
					}
				}

				__m128 hsum1 = _mm_add_ps(_mm256_castps256_ps128(result1), _mm256_extractf128_ps(result1, 1));
				__m128 hsum2 = _mm_add_ps(_mm256_castps256_ps128(result2), _mm256_extractf128_ps(result2, 1));
				__m128 hsum3 = _mm_add_ps(_mm256_castps256_ps128(result3), _mm256_extractf128_ps(result3, 1));
				__m128 hsum4 = _mm_add_ps(_mm256_castps256_ps128(result4), _mm256_extractf128_ps(result4, 1));
				hsum1 = _mm_hadd_ps(_mm_hadd_ps(hsum1, hsum1), _mm_hadd_ps(hsum1, hsum1));
				hsum2 = _mm_hadd_ps(_mm_hadd_ps(hsum2, hsum2), _mm_hadd_ps(hsum2, hsum2));
				hsum3 = _mm_hadd_ps(_mm_hadd_ps(hsum3, hsum3), _mm_hadd_ps(hsum3, hsum3));
				hsum4 = _mm_hadd_ps(_mm_hadd_ps(hsum4, hsum4), _mm_hadd_ps(hsum4, hsum4));
				const T final_res1 = _mm_cvtsi128_si32(_mm_packus_epi32(_mm_cvtps_epi32(hsum1), _mm_setzero_si128()));
				const T final_res2 = _mm_cvtsi128_si32(_mm_packus_epi32(_mm_cvtps_epi32(hsum2), _mm_setzero_si128()));
				const T final_res3 = _mm_cvtsi128_si32(_mm_packus_epi32(_mm_cvtps_epi32(hsum3), _mm_setzero_si128()));
				const T final_res4 = _mm_cvtsi128_si32(_mm_packus_epi32(_mm_cvtps_epi32(hsum4), _mm_setzero_si128()));

				dst1[x] = clamp(final_res1, val_min1, val_max1);
				dst2[x] = clamp(final_res2, val_min2, val_max2);
				dst3[x] = clamp(final_res3, val_min3, val_max3);
				dst4[x] = clamp(final_res4, val_min4, val_max4);
			}
			else
			{
				for (int ly = 0; ly < filter_size_mod2; ly += 2)
				{
					for (int lx = 0; lx < filter_size; lx += 8)
					{
						const __m256 src_ps1 = _mm256_max_ps(_mm256_loadu_ps(reinterpret_cast<const float*>(src_ptr1 + lx)), min_val1);
						const __m256 src_ps1_2 = _mm256_max_ps(_mm256_loadu_ps(reinterpret_cast<const float*>(src_ptr1 + lx + src_pitch1)), min_val1);
						const __m256 src_ps2 = _mm256_max_ps(_mm256_loadu_ps(reinterpret_cast<const float*>(src_ptr2 + lx)), min_val2);
						const __m256 src_ps2_2 = _mm256_max_ps(_mm256_loadu_ps(reinterpret_cast<const float*>(src_ptr2 + lx + src_pitch2)), min_val2);
						const __m256 src_ps3 = _mm256_max_ps(_mm256_loadu_ps(reinterpret_cast<const float*>(src_ptr3 + lx)), min_val3);
						const __m256 src_ps3_2 = _mm256_max_ps(_mm256_loadu_ps(reinterpret_cast<const float*>(src_ptr3 + lx + src_pitch3)), min_val3);
						const __m256 src_ps4 = _mm256_max_ps(_mm256_loadu_ps(reinterpret_cast<const float*>(src_ptr4 + lx)), min_val4);
						const __m256 src_ps4_2 = _mm256_max_ps(_mm256_loadu_ps(reinterpret_cast<const float*>(src_ptr4 + lx + src_pitch4)), min_val4);
						__m256 coeff, coeff2;

						if JincMT_CONSTEXPR (bFP16)
						{
							const int lx2 = lx >> 1;
							__m128i coeff_fp16 = _mm_load_si128((__m128i*)(coeff_ptr + lx2));
							__m128i coeff2_fp16 = _mm_load_si128((__m128i*)(coeff_ptr + lx2 + coeff_stride));

							coeff = _mm256_cvtph_ps(coeff_fp16); // separated grouped converts may be better for instructions grouping at some compilers
							coeff2 = _mm256_cvtph_ps(coeff2_fp16);
						}
						else
						{
							coeff = _mm256_load_ps(coeff_ptr + lx);
							coeff2 = _mm256_load_ps(coeff_ptr + lx + coeff_stride);
						}

						result1 = _mm256_fmadd_ps(src_ps1, coeff, result1);
						result1_2 = _mm256_fmadd_ps(src_ps1_2, coeff2, result1_2);
						result2 = _mm256_fmadd_ps(src_ps2, coeff, result2);
						result2_2 = _mm256_fmadd_ps(src_ps2_2, coeff2, result2_2);
						result3 = _mm256_fmadd_ps(src_ps3, coeff, result3);
						result3_2 = _mm256_fmadd_ps(src_ps3_2, coeff2, result3_2);
						result4 = _mm256_fmadd_ps(src_ps4, coeff, result4);
						result4_2 = _mm256_fmadd_ps(src_ps4_2, coeff2, result4_2);
					}

					coeff_ptr += coeff_stridex2;
					src_ptr1 += src_pitch1x2;
					src_ptr2 += src_pitch2x2;
					src_ptr3 += src_pitch3x2;
					src_ptr4 += src_pitch4x2;
				}

				result1 = _mm256_add_ps(result1, result1_2);
				result2 = _mm256_add_ps(result2, result2_2);
				result3 = _mm256_add_ps(result3, result3_2);
				result4 = _mm256_add_ps(result4, result4_2);

				// last row
				if (fs_notMod2)
				{
					for (int lx = 0; lx < filter_size; lx += 8)
					{
						const __m256 src_ps1 = _mm256_max_ps(_mm256_loadu_ps(reinterpret_cast<const float*>(src_ptr1 + lx)), min_val1);
						const __m256 src_ps2 = _mm256_max_ps(_mm256_loadu_ps(reinterpret_cast<const float*>(src_ptr2 + lx)), min_val2);
						const __m256 src_ps3 = _mm256_max_ps(_mm256_loadu_ps(reinterpret_cast<const float*>(src_ptr3 + lx)), min_val3);
						const __m256 src_ps4 = _mm256_max_ps(_mm256_loadu_ps(reinterpret_cast<const float*>(src_ptr4 + lx)), min_val4);
						__m256 coeff;

						if JincMT_CONSTEXPR (bFP16)
							coeff = _mm256_cvtph_ps(_mm_load_si128((__m128i*)(coeff_ptr + (lx >> 1))));
						else
							coeff = _mm256_load_ps(coeff_ptr + lx);

						result1 = _mm256_fmadd_ps(src_ps1, coeff, result1);
						result2 = _mm256_fmadd_ps(src_ps2, coeff, result2);
						result3 = _mm256_fmadd_ps(src_ps3, coeff, result3);
						result4 = _mm256_fmadd_ps(src_ps4, coeff, result4);
					}
				}

				__m128 hsum1 = _mm_add_ps(_mm256_castps256_ps128(result1), _mm256_extractf128_ps(result1, 1));
				__m128 hsum2 = _mm_add_ps(_mm256_castps256_ps128(result2), _mm256_extractf128_ps(result2, 1));
				__m128 hsum3 = _mm_add_ps(_mm256_castps256_ps128(result3), _mm256_extractf128_ps(result3, 1));
				__m128 hsum4 = _mm_add_ps(_mm256_castps256_ps128(result4), _mm256_extractf128_ps(result4, 1));

				dst1[x] = _mm_cvtss_f32(_mm_hadd_ps(_mm_hadd_ps(hsum1, hsum1), _mm_hadd_ps(hsum1, hsum1)));
				dst2[x] = _mm_cvtss_f32(_mm_hadd_ps(_mm_hadd_ps(hsum2, hsum2), _mm_hadd_ps(hsum2, hsum2)));
				dst3[x] = _mm_cvtss_f32(_mm_hadd_ps(_mm_hadd_ps(hsum3, hsum3), _mm_hadd_ps(hsum3, hsum3)));
				dst4[x] = _mm_cvtss_f32(_mm_hadd_ps(_mm_hadd_ps(hsum4, hsum4), _mm_hadd_ps(hsum4, hsum4)));
			}
			meta++;
		} // for (x)
		meta_y += dst_width;
		dst1 += dst_pitch1;
		dst2 += dst_pitch2;
		dst3 += dst_pitch3;
		dst4 += dst_pitch4;
	} // for (y)
}


template void resize_plane_avx512_1x<uint8_t, false>(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *tab_coeff,
	const float Val_Min[], const float Val_Max[]);
template void resize_plane_avx512_1x<uint16_t, false>(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *tab_coeff,
	const float Val_Min[], const float Val_Max[]);
template void resize_plane_avx512_1x<float, false>(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *tab_coeff,
	const float Val_Min[], const float Val_Max[]);

template void resize_plane_avx512_1x<uint8_t, true>(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *tab_coeff,
	const float Val_Min[], const float Val_Max[]);
template void resize_plane_avx512_1x<uint16_t, true>(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *tab_coeff,
	const float Val_Min[], const float Val_Max[]);
template void resize_plane_avx512_1x<float, true>(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *tab_coeff,
	const float Val_Min[], const float Val_Max[]);

template void resize_plane_avx512_2x<uint8_t, false>(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *tab_coeff,
	const float Val_Min[], const float Val_Max[]);
template void resize_plane_avx512_2x<uint16_t, false>(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *tab_coeff,
	const float Val_Min[], const float Val_Max[]);
template void resize_plane_avx512_2x<float, false>(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *tab_coeff,
	const float Val_Min[], const float Val_Max[]);

template void resize_plane_avx512_2x<uint8_t, true>(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *tab_coeff,
	const float Val_Min[], const float Val_Max[]);
template void resize_plane_avx512_2x<uint16_t, true>(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *tab_coeff,
	const float Val_Min[], const float Val_Max[]);
template void resize_plane_avx512_2x<float, true>(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *tab_coeff,
	const float Val_Min[], const float Val_Max[]);

template void resize_plane_avx512_3x<uint8_t, false>(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *tab_coeff,
	const float Val_Min[], const float Val_Max[]);
template void resize_plane_avx512_3x<uint16_t, false>(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *tab_coeff,
	const float Val_Min[], const float Val_Max[]);
template void resize_plane_avx512_3x<float, false>(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *tab_coeff,
	const float Val_Min[], const float Val_Max[]);

template void resize_plane_avx512_3x<uint8_t, true>(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *tab_coeff,
	const float Val_Min[], const float Val_Max[]);
template void resize_plane_avx512_3x<uint16_t, true>(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *tab_coeff,
	const float Val_Min[], const float Val_Max[]);
template void resize_plane_avx512_3x<float, true>(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *tab_coeff,
	const float Val_Min[], const float Val_Max[]);

template void resize_plane_avx512_4x<uint8_t, false>(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *tab_coeff,
	const float Val_Min[], const float Val_Max[]);
template void resize_plane_avx512_4x<uint16_t, false>(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *tab_coeff,
	const float Val_Min[], const float Val_Max[]);
template void resize_plane_avx512_4x<float, false>(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *tab_coeff,
	const float Val_Min[], const float Val_Max[]);

template void resize_plane_avx512_4x<uint8_t, true>(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *tab_coeff,
	const float Val_Min[], const float Val_Max[]);
template void resize_plane_avx512_4x<uint16_t, true>(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *tab_coeff,
	const float Val_Min[], const float Val_Max[]);
template void resize_plane_avx512_4x<float, true>(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *tab_coeff,
	const float Val_Min[], const float Val_Max[]);

template void resize_plane_avx2_4x_v2<uint8_t, false>(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *tab_coeff,
	const float Val_Min[], const float Val_Max[]);
template void resize_plane_avx2_4x_v2<uint16_t, false>(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *tab_coeff,
	const float Val_Min[], const float Val_Max[]);
template void resize_plane_avx2_4x_v2<float, false>(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *tab_coeff,
	const float Val_Min[], const float Val_Max[]);

template void resize_plane_avx2_4x_v2<uint8_t, true>(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *tab_coeff,
	const float Val_Min[], const float Val_Max[]);
template void resize_plane_avx2_4x_v2<uint16_t, true>(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *tab_coeff,
	const float Val_Min[], const float Val_Max[]);
template void resize_plane_avx2_4x_v2<float, true>(const MT_Data_Info_JincResizeMT *MT_DataGF, const bool PlaneYMode, const EWAPixelCoeff *tab_coeff,
	const float Val_Min[], const float Val_Max[]);

#endif