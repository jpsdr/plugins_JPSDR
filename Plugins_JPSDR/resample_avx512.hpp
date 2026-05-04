// AviSynth+.  Copyright 2026- AviSynth+ Project
// https://avs-plus.net
// http://avisynth.nl
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

/*

This is a common source cpp include file (not header) for multi-arch AVX512 functions.
Functions here are static, they will be compiled into each translation unit including this file.

*/

// Original function needed avx512vbmi feature flag, but we want to support also base AVX512 without VBMI.
// We use _mm512_permutex2var_epi8_SIMUL<UseVBMI> and _mm512_maskz_permutex2var_epi8_SIMUL<UseVBMI>
// Thus both Base AVX512 and ICL level arch is supported.
// We are using two separated source modules and include this hpp file templated
// with UseVBMI/UseVNNI


// helper function for simulating _mm512_permutex2var_epi8 when VBMI is not available
// The MSB bit (128) zeroing effect is _not_ considered here, the indices must be all positive and within 0-127 range.

#ifndef _M_X64
constexpr std::uint64_t make_low_mask64(int bit_count) noexcept
{
    return (bit_count <= 0)  ? 0ULL :
           (bit_count >= 64) ? ~0ULL :
                               ((1ULL << bit_count) - 1ULL);
}
#endif

template<bool UseVBMI>
#if defined(__clang__)
__attribute__((__target__("avx512f,avx512cd,avx512bw,avx512dq,avx512vl,avx512vnni,avx512vbmi,avx512vbmi2,avx512bitalg,avx512vpopcntdq")))
#endif
AVS_FORCEINLINE static __m512i _mm512_permutex2var_epi8_SIMUL(__m512i a, __m512i idx, __m512i b) {
  if constexpr (UseVBMI) {
    return _mm512_permutex2var_epi8(a, idx, b);
  }
  else {
    // Constants
    const __m512i v_one = _mm512_set1_epi16(1);

    // 1. Extract the byte indices for the first 32 and last 32 target pixels
    // We expand them to 16-bit so we can treat them as word-indices
    __m512i idx_lo = _mm512_cvtepu8_epi16(_mm512_castsi512_si256(idx));
    __m512i idx_hi = _mm512_cvtepu8_epi16(_mm512_extracti64x4_epi64(idx, 1));

    // Helper to process 32 bytes of the result at a time
    auto get_32_bytes = [&](__m512i target_idx) {
      // word_idx = byte_idx / 2
      __m512i word_idx = _mm512_srli_epi16(target_idx, 1);

      // VPERMT2W: Full 512-bit cross-lane word shuffle from 128-byte pool [a, b]
      __m512i words = _mm512_permutex2var_epi16(a, word_idx, b);

      // If the original byte index was odd, we need the High Byte of the word.
      // We shift those words right by 8 to put the High Byte into the Low Byte position.
      __mmask32 mask_odd = _mm512_test_epi16_mask(target_idx, v_one);
      words = _mm512_mask_srli_epi16(words, mask_odd, words, 8);

      // VPMOVWB: Truncates 32 words to 32 bytes LINEARLY (No lane scrambling)
      // Returns a __m256i
      return _mm512_cvtepi16_epi8(words);
      };

    // 2. Build the two 256-bit halves
    __m256i res_0_31 = get_32_bytes(idx_lo);
    __m256i res_32_63 = get_32_bytes(idx_hi);

    // 3. Combine into final __m512i
    return _mm512_inserti64x4(_mm512_castsi256_si512(res_0_31), res_32_63, 1);
  }
}

template<bool UseVBMI>
#if defined(__clang__)
__attribute__((__target__("avx512f,avx512cd,avx512bw,avx512dq,avx512vl,avx512vnni,avx512vbmi,avx512vbmi2,avx512bitalg,avx512vpopcntdq")))
#endif
AVS_FORCEINLINE static __m512i _mm512_maskz_permutex2var_epi8_SIMUL(
  __mmask64 k,
  __m512i a,
  __m512i idx,
  __m512i b)
{
  if constexpr (UseVBMI) {
    return _mm512_maskz_permutex2var_epi8(k, a, idx, b);
  }
  else {
    // 1. Run the base simulation to get the permuted bytes
    // Note: Using your existing logic to get the full 512-bit permuted result
    __m512i res = _mm512_permutex2var_epi8_SIMUL<false>(a, idx, b);

    // 2. Apply the zero-mask using AVX-512BW
    // _mm512_maskz_mov_epi8(k, src, src) returns 'src' where k[i]==1, and 0 where k[i]==0
    return _mm512_maskz_mov_epi8(k, res);
  }
}

// H-Float-Resampler: 16 pixels, filter size 4, transpose 4x (4x_m128) to 4x_m512
// Transposes a 4x4 matrix of 4-float vectors (16x16 float matrix effectively).
// Input/Output: Four 512-bit vectors (16 floats each) passed by reference.
#if defined(__clang__)
__attribute__((__target__("avx512f,avx512cd,avx512bw,avx512dq,avx512vl,avx512vnni,avx512vbmi,avx512vbmi2,avx512bitalg,avx512vpopcntdq")))
#endif
AVS_FORCEINLINE static void _MM_TRANSPOSE16_LANE4_PS(__m512& row0, __m512& row1, __m512& row2, __m512& row3) {
  // Stage 1: Interleave 32-bit (float) elements within 128-bit chunks (lanes)
  // t0 = (r0_lo, r1_lo) | t1 = (r0_hi, r1_hi)
  // t2 = (r2_lo, r3_lo) | t3 = (r2_hi, r3_hi)
  __m512 t0 = _mm512_unpacklo_ps(row0, row1);
  __m512 t1 = _mm512_unpackhi_ps(row0, row1);
  __m512 t2 = _mm512_unpacklo_ps(row2, row3);
  __m512 t3 = _mm512_unpackhi_ps(row2, row3);

  // Stage 2: Shuffle 128-bit chunks (lanes) to complete the transpose
  // We use _mm512_shuffle_ps which shuffles 64-bit blocks across the 512-bit register.
  // _MM_SHUFFLE(w, z, y, x) applies to the 64-bit pairs (4 floats) within each 128-bit lane.
  // Result: row0 = columns 0, 1, 2, 3
  row0 = _mm512_shuffle_ps(t0, t2, _MM_SHUFFLE(1, 0, 1, 0));
  // Result: row1 = columns 4, 5, 6, 7
  row1 = _mm512_shuffle_ps(t0, t2, _MM_SHUFFLE(3, 2, 3, 2));
  // Result: row2 = columns 8, 9, 10, 11
  row2 = _mm512_shuffle_ps(t1, t3, _MM_SHUFFLE(1, 0, 1, 0));
  // Result: row3 = columns 12, 13, 14, 15
  row3 = _mm512_shuffle_ps(t1, t3, _MM_SHUFFLE(3, 2, 3, 2));
}

// H-float-resampler: 16 pixels, filter size 8, transpose 8x (2x_m256) to 8x_m512
// Transposes an 8x8 matrix of 2-float vectors (16x16 float matrix).
// Input/Output: Eight 512-bit vectors (16 floats each) passed by reference.
#if defined(__clang__)
__attribute__((__target__("avx512f,avx512cd,avx512bw,avx512dq,avx512vl,avx512vnni,avx512vbmi,avx512vbmi2,avx512bitalg,avx512vpopcntdq")))
#endif
AVS_FORCEINLINE static void _MM_TRANSPOSE8x16_PS(
  __m512& r0, __m512& r1, __m512& r2, __m512& r3,
  __m512& r4, __m512& r5, __m512& r6, __m512& r7)
{
  // --- Stage 1: Unpack 32-bit (Pairs of rows) ---
  __m512 t0 = _mm512_unpacklo_ps(r0, r1);
  __m512 t1 = _mm512_unpackhi_ps(r0, r1);
  __m512 t2 = _mm512_unpacklo_ps(r2, r3);
  __m512 t3 = _mm512_unpackhi_ps(r2, r3);
  __m512 t4 = _mm512_unpacklo_ps(r4, r5);
  __m512 t5 = _mm512_unpackhi_ps(r4, r5);
  __m512 t6 = _mm512_unpacklo_ps(r6, r7);
  __m512 t7 = _mm512_unpackhi_ps(r6, r7);

  // --- Stage 2: Unpack 64-bit (Quads of rows) ---
  // Uses _mm512_unpacklo/hi_pd for 64-bit (double) to interleave pairs of __m512 floats
  __m512 u0 = _mm512_castpd_ps(_mm512_unpacklo_pd(_mm512_castps_pd(t0), _mm512_castps_pd(t2)));
  __m512 u1 = _mm512_castpd_ps(_mm512_unpackhi_pd(_mm512_castps_pd(t0), _mm512_castps_pd(t2)));
  __m512 u2 = _mm512_castpd_ps(_mm512_unpacklo_pd(_mm512_castps_pd(t1), _mm512_castps_pd(t3)));
  __m512 u3 = _mm512_castpd_ps(_mm512_unpackhi_pd(_mm512_castps_pd(t1), _mm512_castps_pd(t3)));
  __m512 u4 = _mm512_castpd_ps(_mm512_unpacklo_pd(_mm512_castps_pd(t4), _mm512_castps_pd(t6)));
  __m512 u5 = _mm512_castpd_ps(_mm512_unpackhi_pd(_mm512_castps_pd(t4), _mm512_castps_pd(t6)));
  __m512 u6 = _mm512_castpd_ps(_mm512_unpacklo_pd(_mm512_castps_pd(t5), _mm512_castps_pd(t7)));
  __m512 u7 = _mm512_castpd_ps(_mm512_unpackhi_pd(_mm512_castps_pd(t5), _mm512_castps_pd(t7)));

  // --- Stage 3: Shuffle 128-bit lanes (Octets of rows) ---
  // _mm512_shuffle_f32x4 shuffles the 128-bit (f32x4) sub-vectors within and between two __m512 vectors.
  // 0x88 = (10001000)_2: selects lane 0 from first input and lane 0 from second input for lo/hi 256 bits.
  // 0xDD = (11011101)_2: selects lane 3 from first input and lane 3 from second input for lo/hi 256 bits.
  __m512 v0 = _mm512_shuffle_f32x4(u0, u4, 0x88); // Col 0, 4 (interleaved)
  __m512 v1 = _mm512_shuffle_f32x4(u0, u4, 0xDD); // Col 1, 5 (interleaved)
  __m512 v2 = _mm512_shuffle_f32x4(u1, u5, 0x88); // Col 2, 6 (interleaved)
  __m512 v3 = _mm512_shuffle_f32x4(u1, u5, 0xDD); // Col 3, 7 (interleaved)
  __m512 v4 = _mm512_shuffle_f32x4(u2, u6, 0x88); // Col 8, 12 (interleaved)
  __m512 v5 = _mm512_shuffle_f32x4(u2, u6, 0xDD); // Col 9, 13 (interleaved)
  __m512 v6 = _mm512_shuffle_f32x4(u3, u7, 0x88); // Col 10, 14 (interleaved)
  __m512 v7 = _mm512_shuffle_f32x4(u3, u7, 0xDD); // Col 11, 15 (interleaved)

  // --- Stage 4: Permute to Linearize Indices ---
  // Corrects the order of the 128-bit lanes to linearize the columns.
  // The columns are currently: (0, 4, 8, 12, 1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11, 15)
  // The required order is: (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)
  __m512i idx = _mm512_setr_epi32(
    0, 4, 1, 5,    /* Lane 0: Rows 0, 1, 2, 3 */
    2, 6, 3, 7,    /* Lane 1: Rows 4, 5, 6, 7 */
    8, 12, 9, 13, /* Lane 2: Rows 8, 9, 10, 11 */
    10, 14, 11, 15 /* Lane 3: Rows 12, 13, 14, 15 */
  );

  // --- Final Assignment with Correct Mapping ---
  // Maps the permuted vector components back to the original row variables (now columns).
  r0 = _mm512_permutexvar_ps(idx, v0); // Col 0
  r1 = _mm512_permutexvar_ps(idx, v2); // Col 1
  r2 = _mm512_permutexvar_ps(idx, v4); // Col 2
  r3 = _mm512_permutexvar_ps(idx, v6); // Col 3
  r4 = _mm512_permutexvar_ps(idx, v1); // Col 4
  r5 = _mm512_permutexvar_ps(idx, v3); // Col 5
  r6 = _mm512_permutexvar_ps(idx, v5); // Col 6
  r7 = _mm512_permutexvar_ps(idx, v7); // Col 7
}

// Loads two 256-bit float vectors from registers (__m256) into a single 512-bit register.
// Equivalent to _mm512_insertf32x8(_mm512_castps256_ps512(lo), hi, 1)
#if defined(__clang__)
__attribute__((__target__("avx512f,avx512cd,avx512bw,avx512dq,avx512vl,avx512vnni,avx512vbmi,avx512vbmi2,avx512bitalg,avx512vpopcntdq")))
#endif
AVS_FORCEINLINE static __m512 _mm512_insert_2_m256(__m256 lo, __m256 hi) {
  return _mm512_insertf32x8(_mm512_castps256_ps512(lo), hi, 1);
}

// Loads four 128-bit float vectors (unaligned) into a single 512-bit register.
#if defined(__clang__)
__attribute__((__target__("avx512f,avx512cd,avx512bw,avx512dq,avx512vl,avx512vnni,avx512vbmi,avx512vbmi2,avx512bitalg,avx512vpopcntdq")))
#endif
AVS_FORCEINLINE static __m512 _mm512_loadu_4_m128(
  /* __m128 const* */ const float* addr1,
  /* __m128 const* */ const float* addr2,
  /* __m128 const* */ const float* addr3,
  /* __m128 const* */ const float* addr4)
{
  // The cast is needed for the first insertion to make the target a 512-bit register
  __m512 v = _mm512_castps128_ps512(_mm_loadu_ps(addr1));
  v = _mm512_insertf32x4(v, _mm_loadu_ps(addr2), 1);
  v = _mm512_insertf32x4(v, _mm_loadu_ps(addr3), 2);
  v = _mm512_insertf32x4(v, _mm_loadu_ps(addr4), 3);
  return v;
}

// Loads four 128-bit float vectors (aligned) into a single 512-bit register.
#if defined(__clang__)
__attribute__((__target__("avx512f,avx512cd,avx512bw,avx512dq,avx512vl,avx512vnni,avx512vbmi,avx512vbmi2,avx512bitalg,avx512vpopcntdq")))
#endif
AVS_FORCEINLINE static __m512 _mm512_load_4_m128(
  /* __m128 const* */ const float* addr1,
  /* __m128 const* */ const float* addr2,
  /* __m128 const* */ const float* addr3,
  /* __m128 const* */ const float* addr4)
{
  // The cast is needed for the first insertion to make the target a 512-bit register
  __m512 v = _mm512_castps128_ps512(_mm_load_ps(addr1));
  v = _mm512_insertf32x4(v, _mm_load_ps(addr2), 1);
  v = _mm512_insertf32x4(v, _mm_load_ps(addr3), 2);
  v = _mm512_insertf32x4(v, _mm_load_ps(addr4), 3);
  return v;
}

// Loads two 256 - bit unaligned integer vectors from registers(__m256i) into a single 512i register.
#if defined(__clang__)
__attribute__((__target__("avx512f,avx512cd,avx512bw,avx512dq,avx512vl,avx512vnni,avx512vbmi,avx512vbmi2,avx512bitalg,avx512vpopcntdq")))
#endif
AVS_FORCEINLINE static __m512i _mm512i_loadu_2_m256i(
  const __m256i* addr1,
  const __m256i* addr2)
{
  return _mm512_inserti64x4(_mm512_castsi256_si512(_mm256_loadu_si256(addr1)), _mm256_loadu_si256(addr2), 1);
}

// Loads two 256 - bit aligned integer vectors from registers(__m256) into a single 512 - bit register.
#if defined(__clang__)
__attribute__((__target__("avx512f,avx512cd,avx512bw,avx512dq,avx512vl,avx512vnni,avx512vbmi,avx512vbmi2,avx512bitalg,avx512vpopcntdq")))
#endif
AVS_FORCEINLINE static __m512i _mm512i_load_2_m256i(
  const __m256i* addr1,
  const __m256i* addr2)
{
  return _mm512_inserti64x4(_mm512_castsi256_si512(_mm256_load_si256(addr1)), _mm256_load_si256(addr2), 1);
}

// Integers
// Loads four 128-bit integer vectors (unaligned) into a single 512-bit integer register.
#if defined(__clang__)
__attribute__((__target__("avx512f,avx512cd,avx512bw,avx512dq,avx512vl,avx512vnni,avx512vbmi,avx512vbmi2,avx512bitalg,avx512vpopcntdq")))
#endif
AVS_FORCEINLINE static __m512i _mm512i_loadu_4_m128i(
  const __m128i* addr1,
  const __m128i* addr2,
  const __m128i* addr3,
  const __m128i* addr4)
{
  // The cast is needed for the first insertion to make the target a 512-bit register
  __m512i v = _mm512_castsi128_si512(_mm_loadu_si128(addr1));
  v = _mm512_inserti32x4(v, _mm_loadu_si128(addr2), 1);
  v = _mm512_inserti32x4(v, _mm_loadu_si128(addr3), 2);
  v = _mm512_inserti32x4(v, _mm_loadu_si128(addr4), 3);
  return v;
}

// filter size up to 4
// 64 target uint8_t pixels at a time
// 128-byte source loads (128 uint8_t pixels)
// maximum permute index is 128 for _mm512_permutex2var_epi8 (uint8_t)
template<bool UseVBMI>
#if defined(__clang__)
__attribute__((__target__("avx512f,avx512cd,avx512bw,avx512dq,avx512vl,avx512vnni,avx512vbmi,avx512vbmi2,avx512bitalg,avx512vpopcntdq")))
#endif
static void resize_h_planar_uint8_avx512_permutex_vstripe_ks4_internal(BYTE* dst8, const BYTE* src8, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height, int bits_per_pixel,const uint8_t range,const bool mode_YUY2)
{
  AVS_UNUSED(bits_per_pixel);

  const int filter_size = program->filter_size; // aligned, practically the coeff table stride

  constexpr int PIXELS_AT_A_TIME = 64;

  // 'source_overread_beyond_targetx' indicates if the filter kernel can read beyond the target width.
  // we load 2x64 source bytes at a time, so ensure safe overread if needed.
  // Our main loop processes calculates for 64 target pixels at a time.
  // Inside that, we load 128 source bytes (2x64) to be able to permutex from that.
  // This we have to check at each mod-PIXELS_AT_A_TIME boundary, the allowance of 128-byte source load.
  const int width_safe_mod = (program->safelimit_128_pixels_each64th_target.overread_possible ? program->safelimit_128_pixels_each64th_target.source_overread_beyond_targetx : width) / PIXELS_AT_A_TIME * PIXELS_AT_A_TIME;

  // Preconditions:

  // 'target_size_alignment' ensures we can safely access coefficients using offsets like
  // 'filter_size * 15' when processing 16 H pixels at a time
  // 'filter_size * 63' when processing 64 H pixels at a time

  // Ensure that coefficient loading beyond the valid target size is safe for 4x4 float loads.
  // We load 4x16bit coeffs at a time

  const int max_scanlines = program->max_scanlines;

  const int val_min = (range==1) ? 0 : 16;
  const int val_max = ((range==1) || (range==4)) ? 255 : (range==2) ? 235 : 240;

  __m256i clamp_limit_min = _mm256_set1_epi16((short)((val_min << 8)|val_min));
  __m256i clamp_limit_max = (mode_YUY2 && ((range>=2) && (range<=3))) ?
    _mm256_set1_epi16((short)(((int)240 << 8)|235)) : _mm256_set1_epi16((short)((val_max << 8)|val_max));	  

  __m512i rounder = _mm512_set1_epi32(1 << (FPScale8bits - 1));

  // Vertical stripe loop for L2 cache optimization
  for (int y_from = 0; y_from < height; y_from += max_scanlines)
  {
    int y_to = std::min(y_from + max_scanlines, height);

    // Reset current_coeff for the start of the stripe (points to start of row's coeffs)
    const short* AVS_RESTRICT current_coeff = program->pixel_coefficient;

    int x = 0;

    // Lambda to handle both safe (fast) and unsafe (masked/partial) loading paths
    auto do_h_integer_core = [&](auto partial_load) {

      // prepare coefs in transposed V-form
      // TODO: make storage in transposed form, 64 x uint16 transposition looks too slow
      // TO FIX: filter_size=4 resize uses filter_size=16 - lots or RAM/cache wasted, in the ready to use transposed form it will be much more optimized

      // 4coefs of 16bit is 64bits, can be loaded as set_epi64 ?
      __m512i coef_0_7 = _mm512_setr_epi64(
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 0),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 1),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 2),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 3),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 4),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 5),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 6),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 7)
      );

      __m512i coef_8_15 = _mm512_setr_epi64(
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 8),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 9),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 10),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 11),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 12),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 13),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 14),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 15)
      );

      __m512i coef_16_23 = _mm512_setr_epi64(
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 16),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 17),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 18),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 19),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 20),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 21),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 22),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 23)
      );

      __m512i coef_24_31 = _mm512_setr_epi64(
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 24),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 25),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 26),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 27),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 28),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 29),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 30),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 31)
      );

      __m512i coef_32_39 = _mm512_setr_epi64(
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 32),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 33),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 34),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 35),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 36),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 37),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 38),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 39)
      );

      __m512i coef_40_47 = _mm512_setr_epi64(
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 40),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 41),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 42),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 43),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 44),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 45),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 46),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 47)
      );

      __m512i coef_48_55 = _mm512_setr_epi64(
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 48),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 49),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 50),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 51),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 52),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 53),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 54),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 55)
      );

      __m512i coef_56_63 = _mm512_setr_epi64(
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 56),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 57),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 58),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 59),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 60),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 61),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 62),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 63)
      );

      // Transpose with permutex
      __m512i c_perm_0 = _mm512_set_epi16(28 + 32, 24 + 32, 20 + 32, 16 + 32, 12 + 32, 8 + 32, 4 + 32, 0 + 32, 28, 24, 20, 16, 12, 8, 4, 0,
        28 + 32, 24 + 32, 20 + 32, 16 + 32, 12 + 32, 8 + 32, 4 + 32, 0 + 32, 28, 24, 20, 16, 12, 8, 4, 0);
      __m512i one_epi16 = _mm512_set1_epi16(1);
      const __mmask32 k_high = 0xFFFF0000;
      // 0.0 .. 15.0 in low 256, 16.0 .. 31.0 in high 256
      __m512i coef_r0_0_31w = _mm512_mask_blend_epi16(k_high, _mm512_permutex2var_epi16(coef_0_7, c_perm_0, coef_8_15), _mm512_permutex2var_epi16(coef_16_23, c_perm_0, coef_24_31));
      // 32.0 .. 47.0  in low 256, 48.0 .. 63.0 in high 256
      __m512i coef_r0_32_63w = _mm512_mask_blend_epi16(k_high, _mm512_permutex2var_epi16(coef_32_39, c_perm_0, coef_40_47), _mm512_permutex2var_epi16(coef_48_55, c_perm_0, coef_56_63));
      c_perm_0 = _mm512_add_epi16(c_perm_0, one_epi16);
      __m512i coef_r1_0_31w = _mm512_mask_blend_epi16(k_high, _mm512_permutex2var_epi16(coef_0_7, c_perm_0, coef_8_15), _mm512_permutex2var_epi16(coef_16_23, c_perm_0, coef_24_31));
      __m512i coef_r1_32_63w = _mm512_mask_blend_epi16(k_high, _mm512_permutex2var_epi16(coef_32_39, c_perm_0, coef_40_47), _mm512_permutex2var_epi16(coef_48_55, c_perm_0, coef_56_63));
      c_perm_0 = _mm512_add_epi16(c_perm_0, one_epi16);
      __m512i coef_r2_0_31w = _mm512_mask_blend_epi16(k_high, _mm512_permutex2var_epi16(coef_0_7, c_perm_0, coef_8_15), _mm512_permutex2var_epi16(coef_16_23, c_perm_0, coef_24_31));
      __m512i coef_r2_32_63w = _mm512_mask_blend_epi16(k_high, _mm512_permutex2var_epi16(coef_32_39, c_perm_0, coef_40_47), _mm512_permutex2var_epi16(coef_48_55, c_perm_0, coef_56_63));
      c_perm_0 = _mm512_add_epi16(c_perm_0, one_epi16);
      __m512i coef_r3_0_31w = _mm512_mask_blend_epi16(k_high, _mm512_permutex2var_epi16(coef_0_7, c_perm_0, coef_8_15), _mm512_permutex2var_epi16(coef_16_23, c_perm_0, coef_24_31));
      __m512i coef_r3_32_63w = _mm512_mask_blend_epi16(k_high, _mm512_permutex2var_epi16(coef_32_39, c_perm_0, coef_40_47), _mm512_permutex2var_epi16(coef_48_55, c_perm_0, coef_56_63));

      // convert-transpose to H-pairs for madd ? better to do with single permutex in future
      // 8 to 8 512 registers - finally real working coeffs to store in the transposed resampling program for block of 64 target samples
      __m512i coef_r0r1_0_31lo = _mm512_unpacklo_epi16(coef_r0_0_31w, coef_r1_0_31w);
      __m512i coef_r0r1_0_31hi = _mm512_unpackhi_epi16(coef_r0_0_31w, coef_r1_0_31w);

      __m512i coef_r0r1_32_63lo = _mm512_unpacklo_epi16(coef_r0_32_63w, coef_r1_32_63w);
      __m512i coef_r0r1_32_63hi = _mm512_unpackhi_epi16(coef_r0_32_63w, coef_r1_32_63w);

      __m512i coef_r2r3_0_31lo = _mm512_unpacklo_epi16(coef_r2_0_31w, coef_r3_0_31w);
      __m512i coef_r2r3_0_31hi = _mm512_unpackhi_epi16(coef_r2_0_31w, coef_r3_0_31w);

      __m512i coef_r2r3_32_63lo = _mm512_unpacklo_epi16(coef_r2_32_63w, coef_r3_32_63w);
      __m512i coef_r2r3_32_63hi = _mm512_unpackhi_epi16(coef_r2_32_63w, coef_r3_32_63w);

      // TODO: store transposed resampling program coeffs to temp buffer for reusage at each line

      // convert resampling program in H-form into permuting indexes for src transposition in V-form
      __m512i perm_0_0_15 = _mm512_loadu_si512((__m512i*)(&program->pixel_offset[x])); // 16 offsets
      __m512i perm_0_16_31 = _mm512_loadu_si512((__m512i*)(&program->pixel_offset[x + 16])); //  16 offsets
      __m512i perm_0_32_47 = _mm512_loadu_si512((__m512i*)(&program->pixel_offset[x + 32])); //  16 offsets
      __m512i perm_0_48_63 = _mm512_loadu_si512((__m512i*)(&program->pixel_offset[x + 48])); //  16 offsets

      int iStart = program->pixel_offset[x];
      __m512i m512i_Start = _mm512_set1_epi32(iStart);

      perm_0_0_15 = _mm512_sub_epi32(perm_0_0_15, m512i_Start);
      perm_0_16_31 = _mm512_sub_epi32(perm_0_16_31, m512i_Start);
      perm_0_32_47 = _mm512_sub_epi32(perm_0_32_47, m512i_Start);
      perm_0_48_63 = _mm512_sub_epi32(perm_0_48_63, m512i_Start);

      __m256i m256i_perm_0_0_15 = _mm512_cvtepi32_epi16(perm_0_0_15);
      __m256i m256i_perm_0_16_31 = _mm512_cvtepi32_epi16(perm_0_16_31);
      __m256i m256i_perm_0_32_47 = _mm512_cvtepi32_epi16(perm_0_32_47);
      __m256i m256i_perm_0_48_63 = _mm512_cvtepi32_epi16(perm_0_48_63);

      __m128i mm128i_perm_0_0_15 = _mm256_cvtepi16_epi8(m256i_perm_0_0_15);
      __m128i mm128i_perm_0_16_31 = _mm256_cvtepi16_epi8(m256i_perm_0_16_31);
      __m128i mm128i_perm_0_32_47 = _mm256_cvtepi16_epi8(m256i_perm_0_32_47);
      __m128i mm128i_perm_0_48_63 = _mm256_cvtepi16_epi8(m256i_perm_0_48_63);

      // Insert each 128-bit register into the specific lane
      __m512i perm_0 = _mm512_inserti32x4(_mm512_setzero_si512(), mm128i_perm_0_0_15, 0); // Lane 0
      perm_0 = _mm512_inserti32x4(perm_0, mm128i_perm_0_16_31, 1); // Lane 1
      perm_0 = _mm512_inserti32x4(perm_0, mm128i_perm_0_32_47, 2); // Lane 2
      perm_0 = _mm512_inserti32x4(perm_0, mm128i_perm_0_48_63, 3); // Lane 3

      // Taps are contiguous (0, 1, 2, 3), so we increment perm indexes by 1.
      __m512i one_epi8 = _mm512_set1_epi8(1);
      __m512i perm_1 = _mm512_add_epi8(perm_0, one_epi8);
      __m512i perm_2 = _mm512_add_epi8(perm_1, one_epi8);
      __m512i perm_3 = _mm512_add_epi8(perm_2, one_epi8);

      uint8_t* AVS_RESTRICT dst_ptr = dst8 + x + y_from * dst_pitch;
      const uint8_t* src_ptr = src8 + iStart + y_from * src_pitch; // all permute offsets relative to this start offset

      // Calculate remaining pixels for bounds checking in partial_load mode. 1..128 remaining pixels possible.
      const int remaining = program->source_size - iStart;
	  // _bzhi_u64 creates a mask with the lower N bits set. If N >= 64, it returns all ones (~0ULL).
#ifndef _M_X64
      const __mmask64 k1 = make_low_mask64(remaining);
      const __mmask64 k2 = make_low_mask64(remaining - 64);
#else
      const __mmask64 k1 = _bzhi_u64(~0ULL, remaining);
      const __mmask64 k2 = _bzhi_u64(~0ULL, std::max(0, remaining - 64));
#endif

      for (int y = y_from; y < y_to; y++)
      {
        __m512i data_src, data_src2;

        if constexpr (partial_load) {
          // Safe masked loads for the image edge
          data_src = _mm512_maskz_loadu_epi8(k1, src_ptr);
          data_src2 = _mm512_maskz_loadu_epi8(k2, src_ptr + 64);
        }
        else {
          // Fast unaligned loads for the safe zone
          data_src = _mm512_loadu_si512(src_ptr);
          data_src2 = _mm512_loadu_si512(src_ptr + 64);
        }

        __m512i data_0 = _mm512_permutex2var_epi8_SIMUL<UseVBMI>(data_src, perm_0, data_src2);
        __m512i data_1 = _mm512_permutex2var_epi8_SIMUL<UseVBMI>(data_src, perm_1, data_src2);
        __m512i data_2 = _mm512_permutex2var_epi8_SIMUL<UseVBMI>(data_src, perm_2, data_src2);
        __m512i data_3 = _mm512_permutex2var_epi8_SIMUL<UseVBMI>(data_src, perm_3, data_src2);

        // unpack 8->16bit 4 to 8 512bit registers source
        __m512i src_r0_0_31 = _mm512_cvtepu8_epi16(_mm512_castsi512_si256(data_0));
        __m512i src_r0_32_63 = _mm512_cvtepu8_epi16(_mm512_extracti64x4_epi64(data_0, 1));

        __m512i src_r1_0_31 = _mm512_cvtepu8_epi16(_mm512_castsi512_si256(data_1));
        __m512i src_r1_32_63 = _mm512_cvtepu8_epi16(_mm512_extracti64x4_epi64(data_1, 1));

        __m512i src_r2_0_31 = _mm512_cvtepu8_epi16(_mm512_castsi512_si256(data_2));
        __m512i src_r2_32_63 = _mm512_cvtepu8_epi16(_mm512_extracti64x4_epi64(data_2, 1));

        __m512i src_r3_0_31 = _mm512_cvtepu8_epi16(_mm512_castsi512_si256(data_3));
        __m512i src_r3_32_63 = _mm512_cvtepu8_epi16(_mm512_extracti64x4_epi64(data_3, 1));

        // transposition to H-pairs 8 to 8 512bit registers (?)
        __m512i src_r0r1_0_31lo = _mm512_unpacklo_epi16(src_r0_0_31, src_r1_0_31);
        __m512i src_r0r1_0_31hi = _mm512_unpackhi_epi16(src_r0_0_31, src_r1_0_31);

        __m512i src_r0r1_32_63lo = _mm512_unpacklo_epi16(src_r0_32_63, src_r1_32_63);
        __m512i src_r0r1_32_63hi = _mm512_unpackhi_epi16(src_r0_32_63, src_r1_32_63);

        __m512i src_r2r3_0_31lo = _mm512_unpacklo_epi16(src_r2_0_31, src_r3_0_31);
        __m512i src_r2r3_0_31hi = _mm512_unpackhi_epi16(src_r2_0_31, src_r3_0_31);

        __m512i src_r2r3_32_63lo = _mm512_unpacklo_epi16(src_r2_32_63, src_r3_32_63);
        __m512i src_r2r3_32_63hi = _mm512_unpackhi_epi16(src_r2_32_63, src_r3_32_63);

        // making FMA in 32bits accs as in AVX256 V-resize
        __m512i result_0_31lo = _mm512_add_epi32(_mm512_madd_epi16(src_r0r1_0_31lo, coef_r0r1_0_31lo), _mm512_madd_epi16(src_r2r3_0_31lo, coef_r2r3_0_31lo));
        __m512i result_0_31hi = _mm512_add_epi32(_mm512_madd_epi16(src_r0r1_0_31hi, coef_r0r1_0_31hi), _mm512_madd_epi16(src_r2r3_0_31hi, coef_r2r3_0_31hi));

        __m512i result_32_63lo = _mm512_add_epi32(_mm512_madd_epi16(src_r0r1_32_63lo, coef_r0r1_32_63lo), _mm512_madd_epi16(src_r2r3_32_63lo, coef_r2r3_32_63lo));
        __m512i result_32_63hi = _mm512_add_epi32(_mm512_madd_epi16(src_r0r1_32_63hi, coef_r0r1_32_63hi), _mm512_madd_epi16(src_r2r3_32_63hi, coef_r2r3_32_63hi));

        // rounding
        result_0_31lo = _mm512_add_epi32(result_0_31lo, rounder);
        result_0_31hi = _mm512_add_epi32(result_0_31hi, rounder);
        result_32_63lo = _mm512_add_epi32(result_32_63lo, rounder);
        result_32_63hi = _mm512_add_epi32(result_32_63hi, rounder);
        // scaling down
        result_0_31lo = _mm512_srai_epi32(result_0_31lo, FPScale8bits);
        result_0_31hi = _mm512_srai_epi32(result_0_31hi, FPScale8bits);
        result_32_63lo = _mm512_srai_epi32(result_32_63lo, FPScale8bits);
        result_32_63hi = _mm512_srai_epi32(result_32_63hi, FPScale8bits);

        __m512i result_0_31_int16 = _mm512_packus_epi32(result_0_31lo, result_0_31hi);
        __m512i result_32_63_int16 = _mm512_packus_epi32(result_32_63lo, result_32_63hi);

        __m256i result_0_31_u8 = _mm512_cvtusepi16_epi8(result_0_31_int16);
        __m256i result_32_63_u8 = _mm512_cvtusepi16_epi8(result_32_63_int16);
		
		result_0_31_u8 = _mm256_min_epu8(result_0_31_u8, clamp_limit_max);
		result_0_31_u8 = _mm256_max_epu8(result_0_31_u8, clamp_limit_min);
		result_32_63_u8 = _mm256_min_epu8(result_32_63_u8, clamp_limit_max);
		result_32_63_u8 = _mm256_max_epu8(result_32_63_u8, clamp_limit_min);

        // cast is enough, no need to use zeroextend
        _mm512_stream_si512(reinterpret_cast<__m512i*>(dst_ptr),
          _mm512_inserti64x4(_mm512_castsi256_si512(result_0_31_u8), result_32_63_u8, 1));

        dst_ptr += dst_pitch;
        src_ptr += src_pitch;
      }

      current_coeff += filter_size * PIXELS_AT_A_TIME;
      };

    // Process the 'safe zone' where direct full unaligned loads are acceptable.
    for (; x < width_safe_mod; x += PIXELS_AT_A_TIME)
    {
      do_h_integer_core(std::false_type{});
    }

    // Process the potentially 'unsafe zone' near the image edge, using safe masked loading.
    for (; x < width; x += PIXELS_AT_A_TIME)
    {
      do_h_integer_core(std::true_type{});
    }
  }
}

// filter size up to 8
// 64 target uint8_t pixels at a time
// 128-byte source loads (128 uint8_t pixels)
// maximum permute index is 128 for _mm512_permutex2var_epi8(uint8_t)
template<bool UseVBMI>
#if defined(__clang__)
__attribute__((__target__("avx512f,avx512cd,avx512bw,avx512dq,avx512vl,avx512vnni,avx512vbmi,avx512vbmi2,avx512bitalg,avx512vpopcntdq")))
#endif
static void resize_h_planar_uint8_avx512_permutex_vstripe_ks8_internal(BYTE* dst8, const BYTE* src8, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height, int bits_per_pixel,const uint8_t range,const bool mode_YUY2)
{
  AVS_UNUSED(bits_per_pixel);

  const int filter_size = program->filter_size; // aligned, practically the coeff table stride

  constexpr int PIXELS_AT_A_TIME = 64;

  // 'source_overread_beyond_targetx' indicates if the filter kernel can read beyond the target width.
  // we load 2x64 source bytes at a time, so ensure safe overread if needed.
  // Our main loop processes calculates for 64 target pixels at a time.
  // Inside that, we load 128 source bytes (2x64) to be able to permutex from that.
  // This we have to check at each mod-PIXELS_AT_A_TIME boundary, the allowance of 128-byte source load.
  const int width_safe_mod = (program->safelimit_128_pixels_each64th_target.overread_possible ? program->safelimit_128_pixels_each64th_target.source_overread_beyond_targetx : width) / PIXELS_AT_A_TIME * PIXELS_AT_A_TIME;

  // Preconditions:

  // 'target_size_alignment' ensures we can safely access coefficients using offsets like
  // 'filter_size * 15' when processing 16 H pixels at a time
  // 'filter_size * 63' when processing 64 H pixels at a time


  // Ensure that coefficient loading beyond the valid target size is safe for 4x8 float loads.
  // We load 8x 'short' coeffs at a time
  // Loading is unaligned, but we fill __m128 registers before combining into __m512

  const int max_scanlines = program->max_scanlines;

  const int val_min = (range==1) ? 0 : 16;
  const int val_max = ((range==1) || (range==4)) ? 255 : (range==2) ? 235 : 240;

  __m256i clamp_limit_min = _mm256_set1_epi16((short)((val_min << 8)|val_min));
  __m256i clamp_limit_max = (mode_YUY2 && ((range>=2) && (range<=3))) ?
    _mm256_set1_epi16((short)(((int)240 << 8)|235)) : _mm256_set1_epi16((short)((val_max << 8)|val_max));	  

  __m512i rounder = _mm512_set1_epi32(1 << (FPScale8bits - 1));

  // Vertical stripe loop for L2 cache optimization
  for (int y_from = 0; y_from < height; y_from += max_scanlines)
  {
    int y_to = std::min(y_from + max_scanlines, height);

    // Reset current_coeff for the start of the stripe (points to start of row's coeffs)
    const short* AVS_RESTRICT current_coeff = program->pixel_coefficient;

    int x = 0;

    // Lambda to handle both safe (fast) and unsafe (masked/partial) loading paths
    auto do_h_integer_core = [&](auto partial_load) {

      // prepare coefs in transposed V-form
      // TODO: make storage in transposed form, 64 x uint16 transposition looks too slow

      // 8coefs of 16bit is 128bits 
      __m512i coef_0_3 = _mm512i_loadu_4_m128i(
        (__m128i*)(current_coeff + filter_size * 0), (__m128i*)(current_coeff + filter_size * 1), (__m128i*)(current_coeff + filter_size * 2), (__m128i*)(current_coeff + filter_size * 3));
      __m512i coef_4_7 = _mm512i_loadu_4_m128i(
        (__m128i*)(current_coeff + filter_size * 4), (__m128i*)(current_coeff + filter_size * 5), (__m128i*)(current_coeff + filter_size * 6), (__m128i*)(current_coeff + filter_size * 7));
      __m512i coef_8_11 = _mm512i_loadu_4_m128i(
        (__m128i*)(current_coeff + filter_size * 8), (__m128i*)(current_coeff + filter_size * 9), (__m128i*)(current_coeff + filter_size * 10), (__m128i*)(current_coeff + filter_size * 11));
      __m512i coef_12_15 = _mm512i_loadu_4_m128i(
        (__m128i*)(current_coeff + filter_size * 12), (__m128i*)(current_coeff + filter_size * 13), (__m128i*)(current_coeff + filter_size * 14), (__m128i*)(current_coeff + filter_size * 15));
      __m512i coef_16_19 = _mm512i_loadu_4_m128i(
        (__m128i*)(current_coeff + filter_size * 16), (__m128i*)(current_coeff + filter_size * 17), (__m128i*)(current_coeff + filter_size * 18), (__m128i*)(current_coeff + filter_size * 19));
      __m512i coef_20_23 = _mm512i_loadu_4_m128i(
        (__m128i*)(current_coeff + filter_size * 20), (__m128i*)(current_coeff + filter_size * 21), (__m128i*)(current_coeff + filter_size * 22), (__m128i*)(current_coeff + filter_size * 23));
      __m512i coef_24_27 = _mm512i_loadu_4_m128i(
        (__m128i*)(current_coeff + filter_size * 24), (__m128i*)(current_coeff + filter_size * 25), (__m128i*)(current_coeff + filter_size * 26), (__m128i*)(current_coeff + filter_size * 27));
      __m512i coef_28_31 = _mm512i_loadu_4_m128i(
        (__m128i*)(current_coeff + filter_size * 28), (__m128i*)(current_coeff + filter_size * 29), (__m128i*)(current_coeff + filter_size * 30), (__m128i*)(current_coeff + filter_size * 31));

      __m512i coef_32_35 = _mm512i_loadu_4_m128i(
        (__m128i*)(current_coeff + filter_size * 32), (__m128i*)(current_coeff + filter_size * 33), (__m128i*)(current_coeff + filter_size * 34), (__m128i*)(current_coeff + filter_size * 35));
      __m512i coef_36_39 = _mm512i_loadu_4_m128i(
        (__m128i*)(current_coeff + filter_size * 36), (__m128i*)(current_coeff + filter_size * 37), (__m128i*)(current_coeff + filter_size * 38), (__m128i*)(current_coeff + filter_size * 39));
      __m512i coef_40_43 = _mm512i_loadu_4_m128i(
        (__m128i*)(current_coeff + filter_size * 40), (__m128i*)(current_coeff + filter_size * 41), (__m128i*)(current_coeff + filter_size * 42), (__m128i*)(current_coeff + filter_size * 43));
      __m512i coef_44_47 = _mm512i_loadu_4_m128i(
        (__m128i*)(current_coeff + filter_size * 44), (__m128i*)(current_coeff + filter_size * 45), (__m128i*)(current_coeff + filter_size * 46), (__m128i*)(current_coeff + filter_size * 47));
      __m512i coef_48_51 = _mm512i_loadu_4_m128i(
        (__m128i*)(current_coeff + filter_size * 48), (__m128i*)(current_coeff + filter_size * 49), (__m128i*)(current_coeff + filter_size * 50), (__m128i*)(current_coeff + filter_size * 51));
      __m512i coef_52_55 = _mm512i_loadu_4_m128i(
        (__m128i*)(current_coeff + filter_size * 52), (__m128i*)(current_coeff + filter_size * 53), (__m128i*)(current_coeff + filter_size * 54), (__m128i*)(current_coeff + filter_size * 55));
      __m512i coef_56_59 = _mm512i_loadu_4_m128i(
        (__m128i*)(current_coeff + filter_size * 56), (__m128i*)(current_coeff + filter_size * 57), (__m128i*)(current_coeff + filter_size * 58), (__m128i*)(current_coeff + filter_size * 59));
      __m512i coef_60_63 = _mm512i_loadu_4_m128i(
        (__m128i*)(current_coeff + filter_size * 60), (__m128i*)(current_coeff + filter_size * 61), (__m128i*)(current_coeff + filter_size * 62), (__m128i*)(current_coeff + filter_size * 63));

      // Transpose with permutex
      __m512i c_perm_0_7 = _mm512_set_epi16(
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        24 + 32, 16 + 32, 8 + 32, 0 + 32, 24, 16, 8, 0);
      __m512i c_perm_8_15 = _mm512_set_epi16(
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        24 + 32, 16 + 32, 8 + 32, 0 + 32, 24, 16, 8, 0,
        0, 0, 0, 0, 0, 0, 0, 0);
      __m512i c_perm_16_23 = _mm512_set_epi16(
        0, 0, 0, 0, 0, 0, 0, 0,
        24 + 32, 16 + 32, 8 + 32, 0 + 32, 24, 16, 8, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0);
      __m512i c_perm_24_31 = _mm512_set_epi16(
        24 + 32, 16 + 32, 8 + 32, 0 + 32, 24, 16, 8, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0);

      __m512i one_epi16 = _mm512_set1_epi16(1);
      const __mmask32 k_8_15 = 0x0000FF00;
      const __mmask32 k_16_23 = 0x00FF0000;
      const __mmask32 k_24_31 = 0xFF000000;

      auto inc_perms = [&](
        __m512i& c0_7,
        __m512i& c8_15,
        __m512i& c16_23,
        __m512i& c24_31
        ) {
          c0_7 = _mm512_add_epi16(c0_7, one_epi16);
          c8_15 = _mm512_add_epi16(c8_15, one_epi16);
          c16_23 = _mm512_add_epi16(c16_23, one_epi16);
          c24_31 = _mm512_add_epi16(c24_31, one_epi16);
        };

      auto make_row_0_63 = [&](
        __m512i& row_0_31w, __m512i& row_32_63w,
        __m512i c0_7, __m512i c8_15, __m512i c16_23, __m512i c24_31
        ) {
          // 0..31
          row_0_31w = _mm512_mask_blend_epi16(
            k_8_15,
            _mm512_permutex2var_epi16(coef_0_3, c0_7, coef_4_7),
            _mm512_permutex2var_epi16(coef_8_11, c8_15, coef_12_15)
          );
          row_0_31w = _mm512_mask_blend_epi16(
            k_16_23,
            row_0_31w,
            _mm512_permutex2var_epi16(coef_16_19, c16_23, coef_20_23)
          );
          row_0_31w = _mm512_mask_blend_epi16(
            k_24_31,
            row_0_31w,
            _mm512_permutex2var_epi16(coef_24_27, c24_31, coef_28_31)
          );

          // 32..63
          row_32_63w = _mm512_mask_blend_epi16(
            k_8_15,
            _mm512_permutex2var_epi16(coef_32_35, c0_7, coef_36_39),
            _mm512_permutex2var_epi16(coef_40_43, c8_15, coef_44_47)
          );
          row_32_63w = _mm512_mask_blend_epi16(
            k_16_23,
            row_32_63w,
            _mm512_permutex2var_epi16(coef_48_51, c16_23, coef_52_55)
          );
          row_32_63w = _mm512_mask_blend_epi16(
            k_24_31,
            row_32_63w,
            _mm512_permutex2var_epi16(coef_56_59, c24_31, coef_60_63)
          );
        };

      __m512i coef_r0_0_31w, coef_r0_32_63w;
      __m512i coef_r1_0_31w, coef_r1_32_63w;
      __m512i coef_r2_0_31w, coef_r2_32_63w;
      __m512i coef_r3_0_31w, coef_r3_32_63w;
      __m512i coef_r4_0_31w, coef_r4_32_63w;
      __m512i coef_r5_0_31w, coef_r5_32_63w;
      __m512i coef_r6_0_31w, coef_r6_32_63w;
      __m512i coef_r7_0_31w, coef_r7_32_63w;

      // r0
      make_row_0_63(coef_r0_0_31w, coef_r0_32_63w, c_perm_0_7, c_perm_8_15, c_perm_16_23, c_perm_24_31);
      inc_perms(c_perm_0_7, c_perm_8_15, c_perm_16_23, c_perm_24_31);
      // r1
      make_row_0_63(coef_r1_0_31w, coef_r1_32_63w, c_perm_0_7, c_perm_8_15, c_perm_16_23, c_perm_24_31);
      inc_perms(c_perm_0_7, c_perm_8_15, c_perm_16_23, c_perm_24_31);
      // r2
      make_row_0_63(coef_r2_0_31w, coef_r2_32_63w, c_perm_0_7, c_perm_8_15, c_perm_16_23, c_perm_24_31);
      inc_perms(c_perm_0_7, c_perm_8_15, c_perm_16_23, c_perm_24_31);
      // r3
      make_row_0_63(coef_r3_0_31w, coef_r3_32_63w, c_perm_0_7, c_perm_8_15, c_perm_16_23, c_perm_24_31);
      inc_perms(c_perm_0_7, c_perm_8_15, c_perm_16_23, c_perm_24_31);
      // r4
      make_row_0_63(coef_r4_0_31w, coef_r4_32_63w, c_perm_0_7, c_perm_8_15, c_perm_16_23, c_perm_24_31);
      inc_perms(c_perm_0_7, c_perm_8_15, c_perm_16_23, c_perm_24_31);
      // r5
      make_row_0_63(coef_r5_0_31w, coef_r5_32_63w, c_perm_0_7, c_perm_8_15, c_perm_16_23, c_perm_24_31);
      inc_perms(c_perm_0_7, c_perm_8_15, c_perm_16_23, c_perm_24_31);
      // r6
      make_row_0_63(coef_r6_0_31w, coef_r6_32_63w, c_perm_0_7, c_perm_8_15, c_perm_16_23, c_perm_24_31);
      inc_perms(c_perm_0_7, c_perm_8_15, c_perm_16_23, c_perm_24_31);
      // r7
      make_row_0_63(coef_r7_0_31w, coef_r7_32_63w, c_perm_0_7, c_perm_8_15, c_perm_16_23, c_perm_24_31);
      /* // last one, not needed
      inc_perms(c_perm_0_7, c_perm_8_15, c_perm_16_23, c_perm_24_31);
      */

      // convert-transpose to H-pairs for madd ? better to do with single permutex in future
      // 16 to 16 512 registers - finally real working coeffs to store in the transposed resampling program for block of 64 target samples
      __m512i coef_r0r1_0_31lo = _mm512_unpacklo_epi16(coef_r0_0_31w, coef_r1_0_31w);
      __m512i coef_r0r1_0_31hi = _mm512_unpackhi_epi16(coef_r0_0_31w, coef_r1_0_31w);

      __m512i coef_r0r1_32_63lo = _mm512_unpacklo_epi16(coef_r0_32_63w, coef_r1_32_63w);
      __m512i coef_r0r1_32_63hi = _mm512_unpackhi_epi16(coef_r0_32_63w, coef_r1_32_63w);

      __m512i coef_r2r3_0_31lo = _mm512_unpacklo_epi16(coef_r2_0_31w, coef_r3_0_31w);
      __m512i coef_r2r3_0_31hi = _mm512_unpackhi_epi16(coef_r2_0_31w, coef_r3_0_31w);

      __m512i coef_r2r3_32_63lo = _mm512_unpacklo_epi16(coef_r2_32_63w, coef_r3_32_63w);
      __m512i coef_r2r3_32_63hi = _mm512_unpackhi_epi16(coef_r2_32_63w, coef_r3_32_63w);

      __m512i coef_r4r5_0_31lo = _mm512_unpacklo_epi16(coef_r4_0_31w, coef_r5_0_31w);
      __m512i coef_r4r5_0_31hi = _mm512_unpackhi_epi16(coef_r4_0_31w, coef_r5_0_31w);

      __m512i coef_r4r5_32_63lo = _mm512_unpacklo_epi16(coef_r4_32_63w, coef_r5_32_63w);
      __m512i coef_r4r5_32_63hi = _mm512_unpackhi_epi16(coef_r4_32_63w, coef_r5_32_63w);

      __m512i coef_r6r7_0_31lo = _mm512_unpacklo_epi16(coef_r6_0_31w, coef_r7_0_31w);
      __m512i coef_r6r7_0_31hi = _mm512_unpackhi_epi16(coef_r6_0_31w, coef_r7_0_31w);

      __m512i coef_r6r7_32_63lo = _mm512_unpacklo_epi16(coef_r6_32_63w, coef_r7_32_63w);
      __m512i coef_r6r7_32_63hi = _mm512_unpackhi_epi16(coef_r6_32_63w, coef_r7_32_63w);

      // TODO: store transposed resampling program coeffs to temp buffer for reusage at each line

      // convert resampling program in H-form into permuting indexes for src transposition in V-form
      __m512i perm_0_0_15 = _mm512_loadu_si512((__m512i*)(&program->pixel_offset[x])); // 16 offsets
      __m512i perm_0_16_31 = _mm512_loadu_si512((__m512i*)(&program->pixel_offset[x + 16])); //  16 offsets
      __m512i perm_0_32_47 = _mm512_loadu_si512((__m512i*)(&program->pixel_offset[x + 32])); //  16 offsets
      __m512i perm_0_48_63 = _mm512_loadu_si512((__m512i*)(&program->pixel_offset[x + 48])); //  16 offsets

      int iStart = program->pixel_offset[x];
      __m512i m512i_Start = _mm512_set1_epi32(iStart);

      perm_0_0_15 = _mm512_sub_epi32(perm_0_0_15, m512i_Start);
      perm_0_16_31 = _mm512_sub_epi32(perm_0_16_31, m512i_Start);
      perm_0_32_47 = _mm512_sub_epi32(perm_0_32_47, m512i_Start);
      perm_0_48_63 = _mm512_sub_epi32(perm_0_48_63, m512i_Start);

      __m256i m256i_perm_0_0_15 = _mm512_cvtepi32_epi16(perm_0_0_15);
      __m256i m256i_perm_0_16_31 = _mm512_cvtepi32_epi16(perm_0_16_31);
      __m256i m256i_perm_0_32_47 = _mm512_cvtepi32_epi16(perm_0_32_47);
      __m256i m256i_perm_0_48_63 = _mm512_cvtepi32_epi16(perm_0_48_63);

      __m128i mm128i_perm_0_0_15 = _mm256_cvtepi16_epi8(m256i_perm_0_0_15);
      __m128i mm128i_perm_0_16_31 = _mm256_cvtepi16_epi8(m256i_perm_0_16_31);
      __m128i mm128i_perm_0_32_47 = _mm256_cvtepi16_epi8(m256i_perm_0_32_47);
      __m128i mm128i_perm_0_48_63 = _mm256_cvtepi16_epi8(m256i_perm_0_48_63);

      // Insert each 128-bit register into the specific lane
      // __m512i perm_0 = _mm512_inserti32x4(_mm512_setzero_si512(), mm128i_perm_0_0_15, 0); // Lane 0
      __m512i perm_0 = _mm512_inserti32x4(_mm512_zextsi128_si512(mm128i_perm_0_0_15), mm128i_perm_0_16_31, 1); // Lane 0+1
      perm_0 = _mm512_inserti32x4(perm_0, mm128i_perm_0_32_47, 2); // Lane 2
      perm_0 = _mm512_inserti32x4(perm_0, mm128i_perm_0_48_63, 3); // Lane 3

      // Taps are contiguous (0, 1, 2, 3), so we increment perm indexes by 1 (in pairs by 2 each).
      const __m512i two_epi8 = _mm512_set1_epi8(2);

      uint8_t* AVS_RESTRICT dst_ptr = dst8 + x + y_from * dst_pitch;
      const uint8_t* src_ptr = src8 + iStart + y_from * src_pitch; // all permute offsets relative to this start offset

      // Calculate remaining pixels for bounds checking in partial_load mode. 1..128 remaining pixels possible.
      const int remaining = program->source_size - iStart;
	  // _bzhi_u64 creates a mask with the lower N bits set. If N >= 64, it returns all ones (~0ULL).
#ifndef _M_X64
      const __mmask64 k1 = make_low_mask64(remaining);
      const __mmask64 k2 = make_low_mask64(remaining - 64);
#else
      const __mmask64 k1 = _bzhi_u64(~0ULL, remaining);
      const __mmask64 k2 = _bzhi_u64(~0ULL, std::max(0, remaining - 64));
#endif

      for (int y = y_from; y < y_to; y++)
      {
        __m512i data_src, data_src2;

        // working permute indexes for advancing to save number of registers used
        __m512i perm_w_even = perm_0;
        __m512i perm_w_odd = _mm512_add_epi8(perm_0, _mm512_set1_epi8(1));

        if constexpr (partial_load) {
          // Safe masked loads for the image edge
          data_src = _mm512_maskz_loadu_epi8(k1, src_ptr);
          data_src2 = _mm512_maskz_loadu_epi8(k2, src_ptr + 64);
        }
        else {
          // Fast unaligned loads for the safe zone
          data_src = _mm512_loadu_si512(src_ptr);
          data_src2 = _mm512_loadu_si512(src_ptr + 64);
        }

        // rows 0..3
        __m512i data_0 = _mm512_permutex2var_epi8_SIMUL<UseVBMI>(data_src, perm_w_even, data_src2);
        __m512i data_1 = _mm512_permutex2var_epi8_SIMUL<UseVBMI>(data_src, perm_w_odd, data_src2);
        perm_w_even = _mm512_add_epi8(perm_w_even, two_epi8);
        perm_w_odd = _mm512_add_epi8(perm_w_odd, two_epi8);

        __m512i data_2 = _mm512_permutex2var_epi8_SIMUL<UseVBMI>(data_src, perm_w_even, data_src2);
        __m512i data_3 = _mm512_permutex2var_epi8_SIMUL<UseVBMI>(data_src, perm_w_odd, data_src2);
        perm_w_even = _mm512_add_epi8(perm_w_even, two_epi8);
        perm_w_odd = _mm512_add_epi8(perm_w_odd, two_epi8);

        // unpack 8->16bit 4 to 8 512bit registers source
        __m512i src_r0_0_31 = _mm512_cvtepu8_epi16(_mm512_castsi512_si256(data_0));
        __m512i src_r0_32_63 = _mm512_cvtepu8_epi16(_mm512_extracti64x4_epi64(data_0, 1));

        __m512i src_r1_0_31 = _mm512_cvtepu8_epi16(_mm512_castsi512_si256(data_1));
        __m512i src_r1_32_63 = _mm512_cvtepu8_epi16(_mm512_extracti64x4_epi64(data_1, 1));

        __m512i src_r2_0_31 = _mm512_cvtepu8_epi16(_mm512_castsi512_si256(data_2));
        __m512i src_r2_32_63 = _mm512_cvtepu8_epi16(_mm512_extracti64x4_epi64(data_2, 1));

        __m512i src_r3_0_31 = _mm512_cvtepu8_epi16(_mm512_castsi512_si256(data_3));
        __m512i src_r3_32_63 = _mm512_cvtepu8_epi16(_mm512_extracti64x4_epi64(data_3, 1));

        // transposition to H-pairs 8 to 8 512bit registers 
        __m512i src_r0r1_0_31lo = _mm512_unpacklo_epi16(src_r0_0_31, src_r1_0_31);
        __m512i src_r0r1_0_31hi = _mm512_unpackhi_epi16(src_r0_0_31, src_r1_0_31);

        __m512i src_r0r1_32_63lo = _mm512_unpacklo_epi16(src_r0_32_63, src_r1_32_63);
        __m512i src_r0r1_32_63hi = _mm512_unpackhi_epi16(src_r0_32_63, src_r1_32_63);

        __m512i src_r2r3_0_31lo = _mm512_unpacklo_epi16(src_r2_0_31, src_r3_0_31);
        __m512i src_r2r3_0_31hi = _mm512_unpackhi_epi16(src_r2_0_31, src_r3_0_31);

        __m512i src_r2r3_32_63lo = _mm512_unpacklo_epi16(src_r2_32_63, src_r3_32_63);
        __m512i src_r2r3_32_63hi = _mm512_unpackhi_epi16(src_r2_32_63, src_r3_32_63);

        // making FMA in 32bits accs as in AVX256 V-resize
        __m512i result_0_31lo = _mm512_add_epi32(_mm512_madd_epi16(src_r0r1_0_31lo, coef_r0r1_0_31lo), _mm512_madd_epi16(src_r2r3_0_31lo, coef_r2r3_0_31lo));
        __m512i result_0_31hi = _mm512_add_epi32(_mm512_madd_epi16(src_r0r1_0_31hi, coef_r0r1_0_31hi), _mm512_madd_epi16(src_r2r3_0_31hi, coef_r2r3_0_31hi));

        __m512i result_32_63lo = _mm512_add_epi32(_mm512_madd_epi16(src_r0r1_32_63lo, coef_r0r1_32_63lo), _mm512_madd_epi16(src_r2r3_32_63lo, coef_r2r3_32_63lo));
        __m512i result_32_63hi = _mm512_add_epi32(_mm512_madd_epi16(src_r0r1_32_63hi, coef_r0r1_32_63hi), _mm512_madd_epi16(src_r2r3_32_63hi, coef_r2r3_32_63hi));

        // rows 4..7
        __m512i data_4 = _mm512_permutex2var_epi8_SIMUL<UseVBMI>(data_src, perm_w_even, data_src2);
        __m512i data_5 = _mm512_permutex2var_epi8_SIMUL<UseVBMI>(data_src, perm_w_odd, data_src2);
        perm_w_even = _mm512_add_epi8(perm_w_even, two_epi8);
        perm_w_odd = _mm512_add_epi8(perm_w_odd, two_epi8);

        __m512i data_6 = _mm512_permutex2var_epi8_SIMUL<UseVBMI>(data_src, perm_w_even, data_src2);
        __m512i data_7 = _mm512_permutex2var_epi8_SIMUL<UseVBMI>(data_src, perm_w_odd, data_src2);

        // unpack 8->16bit 4 to 8 512bit registers source
        __m512i src_r4_0_31 = _mm512_cvtepu8_epi16(_mm512_castsi512_si256(data_4));
        __m512i src_r4_32_63 = _mm512_cvtepu8_epi16(_mm512_extracti64x4_epi64(data_4, 1));

        __m512i src_r5_0_31 = _mm512_cvtepu8_epi16(_mm512_castsi512_si256(data_5));
        __m512i src_r5_32_63 = _mm512_cvtepu8_epi16(_mm512_extracti64x4_epi64(data_5, 1));

        __m512i src_r6_0_31 = _mm512_cvtepu8_epi16(_mm512_castsi512_si256(data_6));
        __m512i src_r6_32_63 = _mm512_cvtepu8_epi16(_mm512_extracti64x4_epi64(data_6, 1));

        __m512i src_r7_0_31 = _mm512_cvtepu8_epi16(_mm512_castsi512_si256(data_7));
        __m512i src_r7_32_63 = _mm512_cvtepu8_epi16(_mm512_extracti64x4_epi64(data_7, 1));

        // transposition to H-pairs 8 to 8 512bit registers 
        __m512i src_r4r5_0_31lo = _mm512_unpacklo_epi16(src_r4_0_31, src_r5_0_31);
        __m512i src_r4r5_0_31hi = _mm512_unpackhi_epi16(src_r4_0_31, src_r5_0_31);

        __m512i src_r4r5_32_63lo = _mm512_unpacklo_epi16(src_r4_32_63, src_r5_32_63);
        __m512i src_r4r5_32_63hi = _mm512_unpackhi_epi16(src_r4_32_63, src_r5_32_63);

        __m512i src_r6r7_0_31lo = _mm512_unpacklo_epi16(src_r6_0_31, src_r7_0_31);
        __m512i src_r6r7_0_31hi = _mm512_unpackhi_epi16(src_r6_0_31, src_r7_0_31);

        __m512i src_r6r7_32_63lo = _mm512_unpacklo_epi16(src_r6_32_63, src_r7_32_63);
        __m512i src_r6r7_32_63hi = _mm512_unpackhi_epi16(src_r6_32_63, src_r7_32_63);

        result_0_31lo = _mm512_add_epi32(_mm512_madd_epi16(src_r4r5_0_31lo, coef_r4r5_0_31lo), result_0_31lo);
        result_0_31hi = _mm512_add_epi32(_mm512_madd_epi16(src_r4r5_0_31hi, coef_r4r5_0_31hi), result_0_31hi);

        result_32_63lo = _mm512_add_epi32(_mm512_madd_epi16(src_r4r5_32_63lo, coef_r4r5_32_63lo), result_32_63lo);
        result_32_63hi = _mm512_add_epi32(_mm512_madd_epi16(src_r4r5_32_63hi, coef_r4r5_32_63hi), result_32_63hi);

        result_0_31lo = _mm512_add_epi32(_mm512_madd_epi16(src_r6r7_0_31lo, coef_r6r7_0_31lo), result_0_31lo);
        result_0_31hi = _mm512_add_epi32(_mm512_madd_epi16(src_r6r7_0_31hi, coef_r6r7_0_31hi), result_0_31hi);

        result_32_63lo = _mm512_add_epi32(_mm512_madd_epi16(src_r6r7_32_63lo, coef_r6r7_32_63lo), result_32_63lo);
        result_32_63hi = _mm512_add_epi32(_mm512_madd_epi16(src_r6r7_32_63hi, coef_r6r7_32_63hi), result_32_63hi);

        // rounding
        result_0_31lo = _mm512_add_epi32(result_0_31lo, rounder);
        result_0_31hi = _mm512_add_epi32(result_0_31hi, rounder);
        result_32_63lo = _mm512_add_epi32(result_32_63lo, rounder);
        result_32_63hi = _mm512_add_epi32(result_32_63hi, rounder);
        // scaling down
        result_0_31lo = _mm512_srai_epi32(result_0_31lo, FPScale8bits);
        result_0_31hi = _mm512_srai_epi32(result_0_31hi, FPScale8bits);
        result_32_63lo = _mm512_srai_epi32(result_32_63lo, FPScale8bits);
        result_32_63hi = _mm512_srai_epi32(result_32_63hi, FPScale8bits);

        __m512i result_0_31_int16 = _mm512_packus_epi32(result_0_31lo, result_0_31hi);
        __m512i result_32_63_int16 = _mm512_packus_epi32(result_32_63lo, result_32_63hi);

        __m256i result_0_31_u8 = _mm512_cvtusepi16_epi8(result_0_31_int16);
        __m256i result_32_63_u8 = _mm512_cvtusepi16_epi8(result_32_63_int16);

		result_0_31_u8 = _mm256_min_epu8(result_0_31_u8, clamp_limit_max);
		result_0_31_u8 = _mm256_max_epu8(result_0_31_u8, clamp_limit_min);
		result_32_63_u8 = _mm256_min_epu8(result_32_63_u8, clamp_limit_max);
		result_32_63_u8 = _mm256_max_epu8(result_32_63_u8, clamp_limit_min);

        // cast is enough, no need to use zeroextend
        _mm512_stream_si512(reinterpret_cast<__m512i*>(dst_ptr),
          _mm512_inserti64x4(_mm512_castsi256_si512(result_0_31_u8), result_32_63_u8, 1));

        dst_ptr += dst_pitch;
        src_ptr += src_pitch;
      }

      current_coeff += filter_size * PIXELS_AT_A_TIME;
      };

    // Process the 'safe zone' where direct full unaligned loads are acceptable.
    for (; x < width_safe_mod; x += PIXELS_AT_A_TIME)
    {
      do_h_integer_core(std::false_type{});
    }

    // Process the potentially 'unsafe zone' near the image edge, using safe masked loading.
    for (; x < width; x += PIXELS_AT_A_TIME)
    {
      do_h_integer_core(std::true_type{});
    }
  }
}

// filter size up to 8
// 64 target uint8_t pixels at a time in 2 groups of 32 to support longer source loading to each group to support lower downsample ratios
// support /2 downsample ratios for resizers with no-resize kernel size of 4 (or support of 2 ?) (Bicubic, Bilinear, and others, also SinPowResize (?))
// 2 groups of 128-byte source loads (128 uint8_t pixels)
// maximum permute index is 128 for _mm512_permutex2var_epi8 (uint8_t)
template<bool UseVBMI>
#if defined(__clang__)
__attribute__((__target__("avx512f,avx512cd,avx512bw,avx512dq,avx512vl,avx512vnni,avx512vbmi,avx512vbmi2,avx512bitalg,avx512vpopcntdq")))
#endif
static void resize_h_planar_uint8_avx512_permutex_vstripe_2s32_ks8_internal(BYTE* dst8, const BYTE* src8, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height, int bits_per_pixel,const uint8_t range,const bool mode_YUY2)
{
  AVS_UNUSED(bits_per_pixel);

  const int filter_size = program->filter_size; // aligned, practically the coeff table stride

  constexpr int PIXELS_AT_A_TIME = 64;

  // 'source_overread_beyond_targetx' indicates if the filter kernel can read beyond the target width.
  // we load 2x64 source bytes at a time, so ensure safe overread if needed.
  // Our main loop processes calculates for 64 target pixels at a time.
  // Inside that, we load 128 source bytes (2x64) to be able to permutex from that.
  // This we have to check at each mod-PIXELS_AT_A_TIME boundary, the allowance of 128-byte source load.
  // Each group of 32 target samples loads 128 source with separate iStart offset (checker function arguments (32/*iSamplesInTheGroup*/, 128/*permutex_index_diff_limit*/, 8/*kernel_size*/)
  // safelimit of 128 each 64 expected is enough ?
  const int width_safe_mod = (program->safelimit_128_pixels_each64th_target.overread_possible ? program->safelimit_128_pixels_each64th_target.source_overread_beyond_targetx : width) / PIXELS_AT_A_TIME * PIXELS_AT_A_TIME;

  // Preconditions:

  // 'target_size_alignment' ensures we can safely access coefficients using offsets like
  // 'filter_size * 15' when processing 16 H pixels at a time
  // 'filter_size * 63' when processing 64 H pixels at a time

  // Ensure that coefficient loading beyond the valid target size is safe for 4x8 float loads.
  // We load 8x 'short' coeffs at a time

  const int max_scanlines = program->max_scanlines;

  const int val_min = (range==1) ? 0 : 16;
  const int val_max = ((range==1) || (range==4)) ? 255 : (range==2) ? 235 : 240;

  __m256i clamp_limit_min = _mm256_set1_epi16((short)((val_min << 8)|val_min));
  __m256i clamp_limit_max = (mode_YUY2 && ((range>=2) && (range<=3))) ?
    _mm256_set1_epi16((short)(((int)240 << 8)|235)) : _mm256_set1_epi16((short)((val_max << 8)|val_max));	  

  __m512i rounder = _mm512_set1_epi32(1 << (FPScale8bits - 1));

  // Vertical stripe loop for L2 cache optimization
  for (int y_from = 0; y_from < height; y_from += max_scanlines)
  {
    int y_to = std::min(y_from + max_scanlines, height);

    // Reset current_coeff for the start of the stripe (points to start of row's coeffs)
    const short* AVS_RESTRICT current_coeff = program->pixel_coefficient;

    int x = 0;

    // Lambda to handle both safe (fast) and unsafe (masked/partial) loading paths
    auto do_h_integer_core = [&](auto partial_load) {

      // prepare coefs in transposed V-form
      // TODO: make storage in transposed form, 64 x uint16 transposition looks too slow

      // 8coefs of 16bit is 128bits 
      __m512i coef_0_3 = _mm512i_loadu_4_m128i(
        (__m128i*)(current_coeff + filter_size * 0), (__m128i*)(current_coeff + filter_size * 1), (__m128i*)(current_coeff + filter_size * 2), (__m128i*)(current_coeff + filter_size * 3));
      __m512i coef_4_7 = _mm512i_loadu_4_m128i(
        (__m128i*)(current_coeff + filter_size * 4), (__m128i*)(current_coeff + filter_size * 5), (__m128i*)(current_coeff + filter_size * 6), (__m128i*)(current_coeff + filter_size * 7));
      __m512i coef_8_11 = _mm512i_loadu_4_m128i(
        (__m128i*)(current_coeff + filter_size * 8), (__m128i*)(current_coeff + filter_size * 9), (__m128i*)(current_coeff + filter_size * 10), (__m128i*)(current_coeff + filter_size * 11));
      __m512i coef_12_15 = _mm512i_loadu_4_m128i(
        (__m128i*)(current_coeff + filter_size * 12), (__m128i*)(current_coeff + filter_size * 13), (__m128i*)(current_coeff + filter_size * 14), (__m128i*)(current_coeff + filter_size * 15));
      __m512i coef_16_19 = _mm512i_loadu_4_m128i(
        (__m128i*)(current_coeff + filter_size * 16), (__m128i*)(current_coeff + filter_size * 17), (__m128i*)(current_coeff + filter_size * 18), (__m128i*)(current_coeff + filter_size * 19));
      __m512i coef_20_23 = _mm512i_loadu_4_m128i(
        (__m128i*)(current_coeff + filter_size * 20), (__m128i*)(current_coeff + filter_size * 21), (__m128i*)(current_coeff + filter_size * 22), (__m128i*)(current_coeff + filter_size * 23));
      __m512i coef_24_27 = _mm512i_loadu_4_m128i(
        (__m128i*)(current_coeff + filter_size * 24), (__m128i*)(current_coeff + filter_size * 25), (__m128i*)(current_coeff + filter_size * 26), (__m128i*)(current_coeff + filter_size * 27));
      __m512i coef_28_31 = _mm512i_loadu_4_m128i(
        (__m128i*)(current_coeff + filter_size * 28), (__m128i*)(current_coeff + filter_size * 29), (__m128i*)(current_coeff + filter_size * 30), (__m128i*)(current_coeff + filter_size * 31));

      __m512i coef_32_35 = _mm512i_loadu_4_m128i(
        (__m128i*)(current_coeff + filter_size * 32), (__m128i*)(current_coeff + filter_size * 33), (__m128i*)(current_coeff + filter_size * 34), (__m128i*)(current_coeff + filter_size * 35));
      __m512i coef_36_39 = _mm512i_loadu_4_m128i(
        (__m128i*)(current_coeff + filter_size * 36), (__m128i*)(current_coeff + filter_size * 37), (__m128i*)(current_coeff + filter_size * 38), (__m128i*)(current_coeff + filter_size * 39));
      __m512i coef_40_43 = _mm512i_loadu_4_m128i(
        (__m128i*)(current_coeff + filter_size * 40), (__m128i*)(current_coeff + filter_size * 41), (__m128i*)(current_coeff + filter_size * 42), (__m128i*)(current_coeff + filter_size * 43));
      __m512i coef_44_47 = _mm512i_loadu_4_m128i(
        (__m128i*)(current_coeff + filter_size * 44), (__m128i*)(current_coeff + filter_size * 45), (__m128i*)(current_coeff + filter_size * 46), (__m128i*)(current_coeff + filter_size * 47));
      __m512i coef_48_51 = _mm512i_loadu_4_m128i(
        (__m128i*)(current_coeff + filter_size * 48), (__m128i*)(current_coeff + filter_size * 49), (__m128i*)(current_coeff + filter_size * 50), (__m128i*)(current_coeff + filter_size * 51));
      __m512i coef_52_55 = _mm512i_loadu_4_m128i(
        (__m128i*)(current_coeff + filter_size * 52), (__m128i*)(current_coeff + filter_size * 53), (__m128i*)(current_coeff + filter_size * 54), (__m128i*)(current_coeff + filter_size * 55));
      __m512i coef_56_59 = _mm512i_loadu_4_m128i(
        (__m128i*)(current_coeff + filter_size * 56), (__m128i*)(current_coeff + filter_size * 57), (__m128i*)(current_coeff + filter_size * 58), (__m128i*)(current_coeff + filter_size * 59));
      __m512i coef_60_63 = _mm512i_loadu_4_m128i(
        (__m128i*)(current_coeff + filter_size * 60), (__m128i*)(current_coeff + filter_size * 61), (__m128i*)(current_coeff + filter_size * 62), (__m128i*)(current_coeff + filter_size * 63));

      // Transpose with permutex
      __m512i c_perm_0_7 = _mm512_set_epi16(
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        24 + 32, 16 + 32, 8 + 32, 0 + 32, 24, 16, 8, 0);
      __m512i c_perm_8_15 = _mm512_set_epi16(
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        24 + 32, 16 + 32, 8 + 32, 0 + 32, 24, 16, 8, 0,
        0, 0, 0, 0, 0, 0, 0, 0);
      __m512i c_perm_16_23 = _mm512_set_epi16(
        0, 0, 0, 0, 0, 0, 0, 0,
        24 + 32, 16 + 32, 8 + 32, 0 + 32, 24, 16, 8, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0);
      __m512i c_perm_24_31 = _mm512_set_epi16(
        24 + 32, 16 + 32, 8 + 32, 0 + 32, 24, 16, 8, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0);

      __m512i one_epi16 = _mm512_set1_epi16(1);
      const __mmask32 k_8_15 = 0x0000FF00;
      const __mmask32 k_16_23 = 0x00FF0000;
      const __mmask32 k_24_31 = 0xFF000000;

      auto inc_perms = [&](
        __m512i& c0_7,
        __m512i& c8_15,
        __m512i& c16_23,
        __m512i& c24_31
        ) {
          c0_7 = _mm512_add_epi16(c0_7, one_epi16);
          c8_15 = _mm512_add_epi16(c8_15, one_epi16);
          c16_23 = _mm512_add_epi16(c16_23, one_epi16);
          c24_31 = _mm512_add_epi16(c24_31, one_epi16);
        };

      auto make_row_0_63 = [&](
        __m512i& row_0_31w, __m512i& row_32_63w,
        __m512i c0_7, __m512i c8_15, __m512i c16_23, __m512i c24_31
        ) {
          // 0..31
          row_0_31w = _mm512_mask_blend_epi16(
            k_8_15,
            _mm512_permutex2var_epi16(coef_0_3, c0_7, coef_4_7),
            _mm512_permutex2var_epi16(coef_8_11, c8_15, coef_12_15)
          );
          row_0_31w = _mm512_mask_blend_epi16(
            k_16_23,
            row_0_31w,
            _mm512_permutex2var_epi16(coef_16_19, c16_23, coef_20_23)
          );
          row_0_31w = _mm512_mask_blend_epi16(
            k_24_31,
            row_0_31w,
            _mm512_permutex2var_epi16(coef_24_27, c24_31, coef_28_31)
          );

          // 32..63
          row_32_63w = _mm512_mask_blend_epi16(
            k_8_15,
            _mm512_permutex2var_epi16(coef_32_35, c0_7, coef_36_39),
            _mm512_permutex2var_epi16(coef_40_43, c8_15, coef_44_47)
          );
          row_32_63w = _mm512_mask_blend_epi16(
            k_16_23,
            row_32_63w,
            _mm512_permutex2var_epi16(coef_48_51, c16_23, coef_52_55)
          );
          row_32_63w = _mm512_mask_blend_epi16(
            k_24_31,
            row_32_63w,
            _mm512_permutex2var_epi16(coef_56_59, c24_31, coef_60_63)
          );
        };

      __m512i coef_r0_0_31w, coef_r0_32_63w;
      __m512i coef_r1_0_31w, coef_r1_32_63w;
      __m512i coef_r2_0_31w, coef_r2_32_63w;
      __m512i coef_r3_0_31w, coef_r3_32_63w;
      __m512i coef_r4_0_31w, coef_r4_32_63w;
      __m512i coef_r5_0_31w, coef_r5_32_63w;
      __m512i coef_r6_0_31w, coef_r6_32_63w;
      __m512i coef_r7_0_31w, coef_r7_32_63w;

      // r0
      make_row_0_63(coef_r0_0_31w, coef_r0_32_63w, c_perm_0_7, c_perm_8_15, c_perm_16_23, c_perm_24_31);
      inc_perms(c_perm_0_7, c_perm_8_15, c_perm_16_23, c_perm_24_31);
      // r1
      make_row_0_63(coef_r1_0_31w, coef_r1_32_63w, c_perm_0_7, c_perm_8_15, c_perm_16_23, c_perm_24_31);
      inc_perms(c_perm_0_7, c_perm_8_15, c_perm_16_23, c_perm_24_31);
      // r2
      make_row_0_63(coef_r2_0_31w, coef_r2_32_63w, c_perm_0_7, c_perm_8_15, c_perm_16_23, c_perm_24_31);
      inc_perms(c_perm_0_7, c_perm_8_15, c_perm_16_23, c_perm_24_31);
      // r3
      make_row_0_63(coef_r3_0_31w, coef_r3_32_63w, c_perm_0_7, c_perm_8_15, c_perm_16_23, c_perm_24_31);
      inc_perms(c_perm_0_7, c_perm_8_15, c_perm_16_23, c_perm_24_31);
      // r4
      make_row_0_63(coef_r4_0_31w, coef_r4_32_63w, c_perm_0_7, c_perm_8_15, c_perm_16_23, c_perm_24_31);
      inc_perms(c_perm_0_7, c_perm_8_15, c_perm_16_23, c_perm_24_31);
      // r5
      make_row_0_63(coef_r5_0_31w, coef_r5_32_63w, c_perm_0_7, c_perm_8_15, c_perm_16_23, c_perm_24_31);
      inc_perms(c_perm_0_7, c_perm_8_15, c_perm_16_23, c_perm_24_31);
      // r6
      make_row_0_63(coef_r6_0_31w, coef_r6_32_63w, c_perm_0_7, c_perm_8_15, c_perm_16_23, c_perm_24_31);
      inc_perms(c_perm_0_7, c_perm_8_15, c_perm_16_23, c_perm_24_31);
      // r7
      make_row_0_63(coef_r7_0_31w, coef_r7_32_63w, c_perm_0_7, c_perm_8_15, c_perm_16_23, c_perm_24_31);
      /* // last one, not needed
      inc_perms(c_perm_0_7, c_perm_8_15, c_perm_16_23, c_perm_24_31);
      */

      // convert-transpose to H-pairs for madd ? better to do with single permutex in future
      // 16 to 16 512 registers - finally real working coeffs to store in the transposed resampling program for block of 64 target samples
      __m512i coef_r0r1_0_31lo = _mm512_unpacklo_epi16(coef_r0_0_31w, coef_r1_0_31w);
      __m512i coef_r0r1_0_31hi = _mm512_unpackhi_epi16(coef_r0_0_31w, coef_r1_0_31w);

      __m512i coef_r0r1_32_63lo = _mm512_unpacklo_epi16(coef_r0_32_63w, coef_r1_32_63w);
      __m512i coef_r0r1_32_63hi = _mm512_unpackhi_epi16(coef_r0_32_63w, coef_r1_32_63w);

      __m512i coef_r2r3_0_31lo = _mm512_unpacklo_epi16(coef_r2_0_31w, coef_r3_0_31w);
      __m512i coef_r2r3_0_31hi = _mm512_unpackhi_epi16(coef_r2_0_31w, coef_r3_0_31w);

      __m512i coef_r2r3_32_63lo = _mm512_unpacklo_epi16(coef_r2_32_63w, coef_r3_32_63w);
      __m512i coef_r2r3_32_63hi = _mm512_unpackhi_epi16(coef_r2_32_63w, coef_r3_32_63w);

      __m512i coef_r4r5_0_31lo = _mm512_unpacklo_epi16(coef_r4_0_31w, coef_r5_0_31w);
      __m512i coef_r4r5_0_31hi = _mm512_unpackhi_epi16(coef_r4_0_31w, coef_r5_0_31w);

      __m512i coef_r4r5_32_63lo = _mm512_unpacklo_epi16(coef_r4_32_63w, coef_r5_32_63w);
      __m512i coef_r4r5_32_63hi = _mm512_unpackhi_epi16(coef_r4_32_63w, coef_r5_32_63w);

      __m512i coef_r6r7_0_31lo = _mm512_unpacklo_epi16(coef_r6_0_31w, coef_r7_0_31w);
      __m512i coef_r6r7_0_31hi = _mm512_unpackhi_epi16(coef_r6_0_31w, coef_r7_0_31w);

      __m512i coef_r6r7_32_63lo = _mm512_unpacklo_epi16(coef_r6_32_63w, coef_r7_32_63w);
      __m512i coef_r6r7_32_63hi = _mm512_unpackhi_epi16(coef_r6_32_63w, coef_r7_32_63w);

      // TODO: store transposed resampling program coeffs to temp buffer for reusage at each line

      // convert resampling program in H-form into permuting indexes for src transposition in V-form
      __m512i perm_0_0_15 = _mm512_loadu_si512((__m512i*)(&program->pixel_offset[x])); // 16 offsets
      __m512i perm_0_16_31 = _mm512_loadu_si512((__m512i*)(&program->pixel_offset[x + 16])); //  16 offsets
      __m512i perm_0_32_47 = _mm512_loadu_si512((__m512i*)(&program->pixel_offset[x + 32])); //  16 offsets
      __m512i perm_0_48_63 = _mm512_loadu_si512((__m512i*)(&program->pixel_offset[x + 48])); //  16 offsets

      int iStart = program->pixel_offset[x];
      __m512i m512i_Start = _mm512_set1_epi32(iStart);
      int iStart_2 = program->pixel_offset[x + 32]; // start offset for second group of 32 targets
      __m512i m512i_Start_2 = _mm512_set1_epi32(iStart_2);

      perm_0_0_15 = _mm512_sub_epi32(perm_0_0_15, m512i_Start);
      perm_0_16_31 = _mm512_sub_epi32(perm_0_16_31, m512i_Start);
      perm_0_32_47 = _mm512_sub_epi32(perm_0_32_47, m512i_Start_2);
      perm_0_48_63 = _mm512_sub_epi32(perm_0_48_63, m512i_Start_2);

      __m256i m256i_perm_0_0_15 = _mm512_cvtepi32_epi16(perm_0_0_15);
      __m256i m256i_perm_0_16_31 = _mm512_cvtepi32_epi16(perm_0_16_31);
      __m256i m256i_perm_0_32_47 = _mm512_cvtepi32_epi16(perm_0_32_47);
      __m256i m256i_perm_0_48_63 = _mm512_cvtepi32_epi16(perm_0_48_63);

      __m128i mm128i_perm_0_0_15 = _mm256_cvtepi16_epi8(m256i_perm_0_0_15);
      __m128i mm128i_perm_0_16_31 = _mm256_cvtepi16_epi8(m256i_perm_0_16_31);
      __m128i mm128i_perm_0_32_47 = _mm256_cvtepi16_epi8(m256i_perm_0_32_47);
      __m128i mm128i_perm_0_48_63 = _mm256_cvtepi16_epi8(m256i_perm_0_48_63);

      // Insert each 128-bit register into the specific lane
//      __m512i perm_0 = _mm512_inserti32x4(_mm512_setzero_si512(), mm128i_perm_0_0_15, 0); // Lane 0
      __m512i perm_0 = _mm512_inserti32x4(_mm512_zextsi128_si512(mm128i_perm_0_0_15), mm128i_perm_0_16_31, 1); // Lane 0+1
      perm_0 = _mm512_inserti32x4(perm_0, mm128i_perm_0_32_47, 2); // Lane 2
      perm_0 = _mm512_inserti32x4(perm_0, mm128i_perm_0_48_63, 3); // Lane 3

      // Taps are contiguous (0, 1, 2, 3), so we increment perm indexes by 1 (in pairs by 2 each).
      const __m512i two_epi8 = _mm512_set1_epi8(2);

      uint8_t* AVS_RESTRICT dst_ptr = dst8 + x + y_from * dst_pitch;
      const uint8_t* src_ptr = src8 + iStart + y_from * src_pitch; // all permute offsets relative to this start offset
      const uint8_t* src_ptr_2 = src8 + iStart_2 + y_from * src_pitch; // all permute offsets in second group of 32 relative to this start offset

      // Calculate remaining pixels for bounds checking in partial_load mode. 1..128 remaining pixels possible.
      const int remaining = program->source_size - iStart;
      const int remaining_2 = program->source_size - iStart_2;
	  // _bzhi_u64 creates a mask with the lower N bits set. If N >= 64, it returns all ones (~0ULL).
#ifndef _M_X64
      const __mmask64 k1 = make_low_mask64(remaining);
      const __mmask64 k2 = make_low_mask64(remaining - 64);

      const __mmask64 k1_2 = make_low_mask64(remaining_2);
      const __mmask64 k2_2 = make_low_mask64(remaining_2 - 64);
#else
      const __mmask64 k1 = _bzhi_u64(~0ULL, remaining);
      const __mmask64 k2 = _bzhi_u64(~0ULL, std::max(0, remaining - 64));

      const __mmask64 k1_2 = _bzhi_u64(~0ULL, remaining_2);
      const __mmask64 k2_2 = _bzhi_u64(~0ULL, std::max(0, remaining_2 - 64));
#endif

      for (int y = y_from; y < y_to; y++)
      {
        __m512i data_src, data_src2;
        __m512i data_src_2, data_src2_2;

        // working permute indexes for advancing to save number of registers used
        __m512i perm_w_even = perm_0;
        __m512i perm_w_odd = _mm512_add_epi8(perm_0, _mm512_set1_epi8(1));

        if constexpr (partial_load) {
          // Safe masked loads for the image edge
          data_src = _mm512_maskz_loadu_epi8(k1, src_ptr);
          data_src2 = _mm512_maskz_loadu_epi8(k2, src_ptr + 64);

          data_src_2 = _mm512_maskz_loadu_epi8(k1_2, src_ptr_2);
          data_src2_2 = _mm512_maskz_loadu_epi8(k2_2, src_ptr_2 + 64);
        }
        else {
          // Fast unaligned loads for the safe zone
          data_src = _mm512_loadu_si512(src_ptr);
          data_src2 = _mm512_loadu_si512(src_ptr + 64);

          data_src_2 = _mm512_loadu_si512(src_ptr_2);
          data_src2_2 = _mm512_loadu_si512(src_ptr_2 + 64);
        }

        // 64 target pixels into two groups: 0-31 and 32-63.
        const __mmask32 k_high = 0xFFFF0000;

        // rows 0..3
        __m512i data_0 = _mm512_mask_blend_epi16(k_high, _mm512_permutex2var_epi8_SIMUL<UseVBMI>(data_src, perm_w_even, data_src2), _mm512_permutex2var_epi8_SIMUL<UseVBMI>(data_src_2, perm_w_even, data_src2_2));
        __m512i data_1 = _mm512_mask_blend_epi16(k_high, _mm512_permutex2var_epi8_SIMUL<UseVBMI>(data_src, perm_w_odd, data_src2), _mm512_permutex2var_epi8_SIMUL<UseVBMI>(data_src_2, perm_w_odd, data_src2_2));
        perm_w_even = _mm512_add_epi8(perm_w_even, two_epi8);
        perm_w_odd = _mm512_add_epi8(perm_w_odd, two_epi8); // Use perm_w_odd here!
        __m512i data_2 = _mm512_mask_blend_epi16(k_high, _mm512_permutex2var_epi8_SIMUL<UseVBMI>(data_src, perm_w_even, data_src2), _mm512_permutex2var_epi8_SIMUL<UseVBMI>(data_src_2, perm_w_even, data_src2_2));
        __m512i data_3 = _mm512_mask_blend_epi16(k_high, _mm512_permutex2var_epi8_SIMUL<UseVBMI>(data_src, perm_w_odd, data_src2), _mm512_permutex2var_epi8_SIMUL<UseVBMI>(data_src_2, perm_w_odd, data_src2_2));

        // unpack 8->16bit 4 to 8 512bit registers source
        __m512i src_r0_0_31 = _mm512_cvtepu8_epi16(_mm512_castsi512_si256(data_0));
        __m512i src_r0_32_63 = _mm512_cvtepu8_epi16(_mm512_extracti64x4_epi64(data_0, 1));

        __m512i src_r1_0_31 = _mm512_cvtepu8_epi16(_mm512_castsi512_si256(data_1));
        __m512i src_r1_32_63 = _mm512_cvtepu8_epi16(_mm512_extracti64x4_epi64(data_1, 1));

        __m512i src_r2_0_31 = _mm512_cvtepu8_epi16(_mm512_castsi512_si256(data_2));
        __m512i src_r2_32_63 = _mm512_cvtepu8_epi16(_mm512_extracti64x4_epi64(data_2, 1));

        __m512i src_r3_0_31 = _mm512_cvtepu8_epi16(_mm512_castsi512_si256(data_3));
        __m512i src_r3_32_63 = _mm512_cvtepu8_epi16(_mm512_extracti64x4_epi64(data_3, 1));

        // transposition to H-pairs 8 to 8 512bit registers 
        __m512i src_r0r1_0_31lo = _mm512_unpacklo_epi16(src_r0_0_31, src_r1_0_31);
        __m512i src_r0r1_0_31hi = _mm512_unpackhi_epi16(src_r0_0_31, src_r1_0_31);

        __m512i src_r0r1_32_63lo = _mm512_unpacklo_epi16(src_r0_32_63, src_r1_32_63);
        __m512i src_r0r1_32_63hi = _mm512_unpackhi_epi16(src_r0_32_63, src_r1_32_63);

        __m512i src_r2r3_0_31lo = _mm512_unpacklo_epi16(src_r2_0_31, src_r3_0_31);
        __m512i src_r2r3_0_31hi = _mm512_unpackhi_epi16(src_r2_0_31, src_r3_0_31);

        __m512i src_r2r3_32_63lo = _mm512_unpacklo_epi16(src_r2_32_63, src_r3_32_63);
        __m512i src_r2r3_32_63hi = _mm512_unpackhi_epi16(src_r2_32_63, src_r3_32_63);

        // making FMA in 32bits accs as in AVX256 V-resize
        __m512i result_0_31lo = _mm512_add_epi32(_mm512_madd_epi16(src_r0r1_0_31lo, coef_r0r1_0_31lo), _mm512_madd_epi16(src_r2r3_0_31lo, coef_r2r3_0_31lo));
        __m512i result_0_31hi = _mm512_add_epi32(_mm512_madd_epi16(src_r0r1_0_31hi, coef_r0r1_0_31hi), _mm512_madd_epi16(src_r2r3_0_31hi, coef_r2r3_0_31hi));

        __m512i result_32_63lo = _mm512_add_epi32(_mm512_madd_epi16(src_r0r1_32_63lo, coef_r0r1_32_63lo), _mm512_madd_epi16(src_r2r3_32_63lo, coef_r2r3_32_63lo));
        __m512i result_32_63hi = _mm512_add_epi32(_mm512_madd_epi16(src_r0r1_32_63hi, coef_r0r1_32_63hi), _mm512_madd_epi16(src_r2r3_32_63hi, coef_r2r3_32_63hi));

        // rows 4..7
        perm_w_even = _mm512_add_epi8(perm_w_even, two_epi8);
        perm_w_odd = _mm512_add_epi8(perm_w_odd, two_epi8); // Use perm_w_odd here!
        __m512i data_4 = _mm512_mask_blend_epi16(k_high, _mm512_permutex2var_epi8_SIMUL<UseVBMI>(data_src, perm_w_even, data_src2), _mm512_permutex2var_epi8_SIMUL<UseVBMI>(data_src_2, perm_w_even, data_src2_2));
        __m512i data_5 = _mm512_mask_blend_epi16(k_high, _mm512_permutex2var_epi8_SIMUL<UseVBMI>(data_src, perm_w_odd, data_src2), _mm512_permutex2var_epi8_SIMUL<UseVBMI>(data_src_2, perm_w_odd, data_src2_2));
        perm_w_even = _mm512_add_epi8(perm_w_even, two_epi8);
        perm_w_odd = _mm512_add_epi8(perm_w_odd, two_epi8);
        __m512i data_6 = _mm512_mask_blend_epi16(k_high, _mm512_permutex2var_epi8_SIMUL<UseVBMI>(data_src, perm_w_even, data_src2), _mm512_permutex2var_epi8_SIMUL<UseVBMI>(data_src_2, perm_w_even, data_src2_2));
        __m512i data_7 = _mm512_mask_blend_epi16(k_high, _mm512_permutex2var_epi8_SIMUL<UseVBMI>(data_src, perm_w_odd, data_src2), _mm512_permutex2var_epi8_SIMUL<UseVBMI>(data_src_2, perm_w_odd, data_src2_2));

        // unpack 8->16bit 4 to 8 512bit registers source
        __m512i src_r4_0_31 = _mm512_cvtepu8_epi16(_mm512_castsi512_si256(data_4));
        __m512i src_r4_32_63 = _mm512_cvtepu8_epi16(_mm512_extracti64x4_epi64(data_4, 1));

        __m512i src_r5_0_31 = _mm512_cvtepu8_epi16(_mm512_castsi512_si256(data_5));
        __m512i src_r5_32_63 = _mm512_cvtepu8_epi16(_mm512_extracti64x4_epi64(data_5, 1));

        __m512i src_r6_0_31 = _mm512_cvtepu8_epi16(_mm512_castsi512_si256(data_6));
        __m512i src_r6_32_63 = _mm512_cvtepu8_epi16(_mm512_extracti64x4_epi64(data_6, 1));

        __m512i src_r7_0_31 = _mm512_cvtepu8_epi16(_mm512_castsi512_si256(data_7));
        __m512i src_r7_32_63 = _mm512_cvtepu8_epi16(_mm512_extracti64x4_epi64(data_7, 1));

        // transposition to H-pairs 8 to 8 512bit registers 
        __m512i src_r4r5_0_31lo = _mm512_unpacklo_epi16(src_r4_0_31, src_r5_0_31);
        __m512i src_r4r5_0_31hi = _mm512_unpackhi_epi16(src_r4_0_31, src_r5_0_31);

        __m512i src_r4r5_32_63lo = _mm512_unpacklo_epi16(src_r4_32_63, src_r5_32_63);
        __m512i src_r4r5_32_63hi = _mm512_unpackhi_epi16(src_r4_32_63, src_r5_32_63);

        __m512i src_r6r7_0_31lo = _mm512_unpacklo_epi16(src_r6_0_31, src_r7_0_31);
        __m512i src_r6r7_0_31hi = _mm512_unpackhi_epi16(src_r6_0_31, src_r7_0_31);

        __m512i src_r6r7_32_63lo = _mm512_unpacklo_epi16(src_r6_32_63, src_r7_32_63);
        __m512i src_r6r7_32_63hi = _mm512_unpackhi_epi16(src_r6_32_63, src_r7_32_63);
        result_0_31lo = _mm512_add_epi32(_mm512_madd_epi16(src_r4r5_0_31lo, coef_r4r5_0_31lo), result_0_31lo);
        result_0_31hi = _mm512_add_epi32(_mm512_madd_epi16(src_r4r5_0_31hi, coef_r4r5_0_31hi), result_0_31hi);

        result_32_63lo = _mm512_add_epi32(_mm512_madd_epi16(src_r4r5_32_63lo, coef_r4r5_32_63lo), result_32_63lo);
        result_32_63hi = _mm512_add_epi32(_mm512_madd_epi16(src_r4r5_32_63hi, coef_r4r5_32_63hi), result_32_63hi);

        result_0_31lo = _mm512_add_epi32(_mm512_madd_epi16(src_r6r7_0_31lo, coef_r6r7_0_31lo), result_0_31lo);
        result_0_31hi = _mm512_add_epi32(_mm512_madd_epi16(src_r6r7_0_31hi, coef_r6r7_0_31hi), result_0_31hi);

        result_32_63lo = _mm512_add_epi32(_mm512_madd_epi16(src_r6r7_32_63lo, coef_r6r7_32_63lo), result_32_63lo);
        result_32_63hi = _mm512_add_epi32(_mm512_madd_epi16(src_r6r7_32_63hi, coef_r6r7_32_63hi), result_32_63hi);

        // rounding
        result_0_31lo = _mm512_add_epi32(result_0_31lo, rounder);
        result_0_31hi = _mm512_add_epi32(result_0_31hi, rounder);
        result_32_63lo = _mm512_add_epi32(result_32_63lo, rounder);
        result_32_63hi = _mm512_add_epi32(result_32_63hi, rounder);
        // scaling down
        result_0_31lo = _mm512_srai_epi32(result_0_31lo, FPScale8bits);
        result_0_31hi = _mm512_srai_epi32(result_0_31hi, FPScale8bits);
        result_32_63lo = _mm512_srai_epi32(result_32_63lo, FPScale8bits);
        result_32_63hi = _mm512_srai_epi32(result_32_63hi, FPScale8bits);

        __m512i result_0_31_int16 = _mm512_packus_epi32(result_0_31lo, result_0_31hi);
        __m512i result_32_63_int16 = _mm512_packus_epi32(result_32_63lo, result_32_63hi);

        __m256i result_0_31_u8 = _mm512_cvtusepi16_epi8(result_0_31_int16);
        __m256i result_32_63_u8 = _mm512_cvtusepi16_epi8(result_32_63_int16);

		result_0_31_u8 = _mm256_min_epu8(result_0_31_u8, clamp_limit_max);
		result_0_31_u8 = _mm256_max_epu8(result_0_31_u8, clamp_limit_min);
		result_32_63_u8 = _mm256_min_epu8(result_32_63_u8, clamp_limit_max);
		result_32_63_u8 = _mm256_max_epu8(result_32_63_u8, clamp_limit_min);

        // cast is enough, no need to use zeroextend
        _mm512_stream_si512(reinterpret_cast<__m512i*>(dst_ptr),
          _mm512_inserti64x4(_mm512_castsi256_si512(result_0_31_u8), result_32_63_u8, 1));

        dst_ptr += dst_pitch;
        src_ptr += src_pitch;
        src_ptr_2 += src_pitch;
      }

      current_coeff += filter_size * PIXELS_AT_A_TIME;
      };

    // Process the 'safe zone' where direct full unaligned loads are acceptable.
    for (; x < width_safe_mod; x += PIXELS_AT_A_TIME)
    {
      do_h_integer_core(std::false_type{});
    }

    // Process the potentially 'unsafe zone' near the image edge, using safe masked loading.
    for (; x < width; x += PIXELS_AT_A_TIME)
    {
      do_h_integer_core(std::true_type{});
    }
  }
}

// filter size up to 16
// 32 target uint8_t pixels at a time
// 128-byte source loads (128 uint8_t pixels)
// maximum permute index is 128 for _mm512_permutex2var_epi8 (uint8_t)
// expect to support all upsampling ratios up to filter support of 8 (or 7..6 ?) and some downsampling ratios with filter support up to 3 (with downsample ratios from 0.5 or a bit lower)
template<bool UseVBMI>
#if defined(__clang__)
__attribute__((__target__("avx512f,avx512cd,avx512bw,avx512dq,avx512vl,avx512vnni,avx512vbmi,avx512vbmi2,avx512bitalg,avx512vpopcntdq")))
#endif
static void resize_h_planar_uint8_avx512_permutex_vstripe_ks16_internal(BYTE* dst8, const BYTE* src8, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height, int bits_per_pixel,const uint8_t range,const bool mode_YUY2)
{
  AVS_UNUSED(bits_per_pixel);

  const int filter_size = program->filter_size; // aligned, practically the coeff table stride

  constexpr int PIXELS_AT_A_TIME = 32;

  // 'source_overread_beyond_targetx' indicates if the filter kernel can read beyond the target width.
  // we load 2x64 source bytes at a time, so ensure safe overread if needed.
  // Our main loop processes calculates for 64 target pixels at a time.
  // Inside that, we load 128 source bytes (2x64) to be able to permutex from that.
  // This we have to check at each mod-PIXELS_AT_A_TIME boundary, the allowance of 128-byte source load.
  const int width_safe_mod = (program->safelimit_128_pixels_each64th_target.overread_possible ? program->safelimit_128_pixels_each64th_target.source_overread_beyond_targetx : width) / PIXELS_AT_A_TIME * PIXELS_AT_A_TIME;

  // Preconditions:

  // 'target_size_alignment' ensures we can safely access coefficients using offsets like
  // 'filter_size * 15' when processing 16 H pixels at a time
  // 'filter_size * 31' when processing 32 H pixels at a time

  // Ensure that coefficient loading beyond the valid target size is safe 
  // We load 16x 'short' coeffs at a time
  // We aligned_load from coeff positions, 32 bytes boundary needed, which is
  // guaranteed by filter_size_alignment >=16

  const int max_scanlines = program->max_scanlines;

  const int val_min = (range==1) ? 0 : 16;
  const int val_max = ((range==1) || (range==4)) ? 255 : (range==2) ? 235 : 240;

  __m256i clamp_limit_min = _mm256_set1_epi16((short)((val_min << 8)|val_min));
  __m256i clamp_limit_max = (mode_YUY2 && ((range>=2) && (range<=3))) ?
    _mm256_set1_epi16((short)(((int)240 << 8)|235)) : _mm256_set1_epi16((short)((val_max << 8)|val_max));	

  __m512i rounder = _mm512_set1_epi32(1 << (FPScale8bits - 1));

  // Vertical stripe loop for L2 cache optimization
  for (int y_from = 0; y_from < height; y_from += max_scanlines)
  {
    int y_to = std::min(y_from + max_scanlines, height);

    // Reset current_coeff for the start of the stripe (points to start of row's coeffs)
    const short* AVS_RESTRICT current_coeff = program->pixel_coefficient;

    int x = 0;

    // Lambda to handle both safe (fast) and unsafe (masked/partial) loading paths
    auto do_h_integer_core = [&](auto partial_load) {

      // prepare coefs in transposed V-form
      // TODO: make storage in transposed form, 32 x uint8 transposition looks too slow

      // 16coefs of 16bit is 256bits - load as pairs of _m256i
      // hope 16-coeffs blocks are 32-bytes aligned for m256i aligned loads ?
      __m512i coef_0_1 = _mm512i_load_2_m256i((__m256i*)(current_coeff + filter_size * 0), (__m256i*)(current_coeff + filter_size * 1));
      __m512i coef_2_3 = _mm512i_load_2_m256i((__m256i*)(current_coeff + filter_size * 2), (__m256i*)(current_coeff + filter_size * 3));
      __m512i coef_4_5 = _mm512i_load_2_m256i((__m256i*)(current_coeff + filter_size * 4), (__m256i*)(current_coeff + filter_size * 5));
      __m512i coef_6_7 = _mm512i_load_2_m256i((__m256i*)(current_coeff + filter_size * 6), (__m256i*)(current_coeff + filter_size * 7));
      __m512i coef_8_9 = _mm512i_load_2_m256i((__m256i*)(current_coeff + filter_size * 8), (__m256i*)(current_coeff + filter_size * 9));
      __m512i coef_10_11 = _mm512i_load_2_m256i((__m256i*)(current_coeff + filter_size * 10), (__m256i*)(current_coeff + filter_size * 11));
      __m512i coef_12_13 = _mm512i_load_2_m256i((__m256i*)(current_coeff + filter_size * 12), (__m256i*)(current_coeff + filter_size * 13));
      __m512i coef_14_15 = _mm512i_load_2_m256i((__m256i*)(current_coeff + filter_size * 14), (__m256i*)(current_coeff + filter_size * 15));
      __m512i coef_16_17 = _mm512i_load_2_m256i((__m256i*)(current_coeff + filter_size * 16), (__m256i*)(current_coeff + filter_size * 17));
      __m512i coef_18_19 = _mm512i_load_2_m256i((__m256i*)(current_coeff + filter_size * 18), (__m256i*)(current_coeff + filter_size * 19));
      __m512i coef_20_21 = _mm512i_load_2_m256i((__m256i*)(current_coeff + filter_size * 20), (__m256i*)(current_coeff + filter_size * 21));
      __m512i coef_22_23 = _mm512i_load_2_m256i((__m256i*)(current_coeff + filter_size * 22), (__m256i*)(current_coeff + filter_size * 23));
      __m512i coef_24_25 = _mm512i_load_2_m256i((__m256i*)(current_coeff + filter_size * 24), (__m256i*)(current_coeff + filter_size * 25));
      __m512i coef_26_27 = _mm512i_load_2_m256i((__m256i*)(current_coeff + filter_size * 26), (__m256i*)(current_coeff + filter_size * 27));
      __m512i coef_28_29 = _mm512i_load_2_m256i((__m256i*)(current_coeff + filter_size * 28), (__m256i*)(current_coeff + filter_size * 29));
      __m512i coef_30_31 = _mm512i_load_2_m256i((__m256i*)(current_coeff + filter_size * 30), (__m256i*)(current_coeff + filter_size * 31));

      // Transpose with permutex
      __m512i c_perm_0_3 = _mm512_set_epi16(
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        16 + 32, 0 + 32, 16, 0);

      __m512i c_perm_4_7 = _mm512_set_epi16(
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        16 + 32, 0 + 32, 16, 0,
        0, 0, 0, 0);

      __m512i c_perm_8_11 = _mm512_set_epi16(
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        16 + 32, 0 + 32, 16, 0,
        0, 0, 0, 0,
        0, 0, 0, 0);

      __m512i c_perm_12_15 = _mm512_set_epi16(
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        16 + 32, 0 + 32, 16, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0);

      __m512i c_perm_16_19 = _mm512_set_epi16(
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        16 + 32, 0 + 32, 16, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0);

      __m512i c_perm_20_23 = _mm512_set_epi16(
        0, 0, 0, 0,
        0, 0, 0, 0,
        16 + 32, 0 + 32, 16, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0);

      __m512i c_perm_24_27 = _mm512_set_epi16(
        0, 0, 0, 0,
        16 + 32, 0 + 32, 16, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0);

      __m512i c_perm_28_31 = _mm512_set_epi16(
        16 + 32, 0 + 32, 16, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0);


      __m512i one_epi16 = _mm512_set1_epi16(1);

      // Define masks for the 4-word (8-byte, 4x 16-bit word) segments within the 32-word vector.
      // Note: Each AVX-512 word (epi16) lane corresponds to one bit in the __mmask32.
      // The masks target chunks of 4 words (k_X_Y represents bits X through Y inclusive).
      const __mmask32 k_4_7 = 0x000000F0;
      const __mmask32 k_8_11 = 0x00000F00;
      const __mmask32 k_12_15 = 0x0000F000;
      const __mmask32 k_16_19 = 0x000F0000;
      const __mmask32 k_20_23 = 0x00F00000;
      const __mmask32 k_24_27 = 0x0F000000;
      const __mmask32 k_28_31 = 0xF0000000;

      // Helper lambda to increment all eight permutation vectors by one.
      auto inc_perms = [&](
        __m512i& p0_3, __m512i& p4_7, __m512i& p8_11, __m512i& p12_15,
        __m512i& p16_19, __m512i& p20_23, __m512i& p24_27, __m512i& p28_31
        ) {
          p0_3 = _mm512_add_epi16(p0_3, one_epi16);
          p4_7 = _mm512_add_epi16(p4_7, one_epi16);
          p8_11 = _mm512_add_epi16(p8_11, one_epi16);
          p12_15 = _mm512_add_epi16(p12_15, one_epi16);
          p16_19 = _mm512_add_epi16(p16_19, one_epi16);
          p20_23 = _mm512_add_epi16(p20_23, one_epi16);
          p24_27 = _mm512_add_epi16(p24_27, one_epi16);
          p28_31 = _mm512_add_epi16(p28_31, one_epi16);
        };

      // Helper lambda to construct one full 32-word (512-bit) coefficient row.
      // It uses mask blending to merge the results of _mm512_permutex2var_epi16
      // for different 4-word segments, using the current permutation vectors.
      auto make_coef_row = [&](
        __m512i& row_result,
        __m512i p0_3, __m512i p4_7, __m512i p8_11, __m512i p12_15,
        __m512i p16_19, __m512i p20_23, __m512i p24_27, __m512i p28_31
        ) {
          // Start with the first segment (words 0-3) and the fourth segment (words 4-7).
          row_result = _mm512_mask_blend_epi16(
            k_4_7,
            _mm512_permutex2var_epi16(coef_0_1, p0_3, coef_2_3),  // words 0-3 (unmasked)
            _mm512_permutex2var_epi16(coef_4_5, p4_7, coef_6_7)   // words 4-7 (masked)
          );

          // Merge segment 8-11
          row_result = _mm512_mask_blend_epi16(
            k_8_11,
            row_result,
            _mm512_permutex2var_epi16(coef_8_9, p8_11, coef_10_11)
          );

          // Merge segment 12-15
          row_result = _mm512_mask_blend_epi16(
            k_12_15,
            row_result,
            _mm512_permutex2var_epi16(coef_12_13, p12_15, coef_14_15)
          );

          // Merge segment 16-19
          row_result = _mm512_mask_blend_epi16(
            k_16_19,
            row_result,
            _mm512_permutex2var_epi16(coef_16_17, p16_19, coef_18_19)
          );

          // Merge segment 20-23
          row_result = _mm512_mask_blend_epi16(
            k_20_23,
            row_result,
            _mm512_permutex2var_epi16(coef_20_21, p20_23, coef_22_23)
          );

          // Merge segment 24-27
          row_result = _mm512_mask_blend_epi16(
            k_24_27,
            row_result,
            _mm512_permutex2var_epi16(coef_24_25, p24_27, coef_26_27)
          );

          // Merge segment 28-31
          row_result = _mm512_mask_blend_epi16(
            k_28_31,
            row_result,
            _mm512_permutex2var_epi16(coef_28_29, p28_31, coef_30_31)
          );
        };

      // Declare the coefficient row variables
      __m512i coef_r0_0_31w, coef_r1_0_31w, coef_r2_0_31w, coef_r3_0_31w;
      __m512i coef_r4_0_31w, coef_r5_0_31w, coef_r6_0_31w, coef_r7_0_31w;
      __m512i coef_r8_0_31w, coef_r9_0_31w, coef_r10_0_31w, coef_r11_0_31w;
      __m512i coef_r12_0_31w, coef_r13_0_31w, coef_r14_0_31w, coef_r15_0_31w;

      // Process Row 0
      make_coef_row(coef_r0_0_31w, c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);
      inc_perms(c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);

      // Process Row 1
      make_coef_row(coef_r1_0_31w, c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);
      inc_perms(c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);

      // Process Row 2
      make_coef_row(coef_r2_0_31w, c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);
      inc_perms(c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);

      // Process Row 3
      make_coef_row(coef_r3_0_31w, c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);
      inc_perms(c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);

      // Process Row 4
      make_coef_row(coef_r4_0_31w, c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);
      inc_perms(c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);

      // Process Row 5
      make_coef_row(coef_r5_0_31w, c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);
      inc_perms(c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);

      // Process Row 6
      make_coef_row(coef_r6_0_31w, c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);
      inc_perms(c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);

      // Process Row 7
      make_coef_row(coef_r7_0_31w, c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);
      inc_perms(c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);

      // Process Row 8
      make_coef_row(coef_r8_0_31w, c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);
      inc_perms(c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);

      // Process Row 9
      make_coef_row(coef_r9_0_31w, c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);
      inc_perms(c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);

      // Process Row 10
      make_coef_row(coef_r10_0_31w, c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);
      inc_perms(c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);

      // Process Row 11
      make_coef_row(coef_r11_0_31w, c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);
      inc_perms(c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);

      // Process Row 12
      make_coef_row(coef_r12_0_31w, c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);
      inc_perms(c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);

      // Process Row 13
      make_coef_row(coef_r13_0_31w, c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);
      inc_perms(c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);

      // Process Row 14
      make_coef_row(coef_r14_0_31w, c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);
      inc_perms(c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);

      // Process Row 15
      make_coef_row(coef_r15_0_31w, c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);
      // No inc_perms here

      // convert-transpose to H-pairs for madd ? better to do with single permutex in future
      // 16 to 16 512 registers - finally real working coeffs to store in the transposed resampling program for block of 64 target samples
      __m512i coef_r0r1_0_31lo = _mm512_unpacklo_epi16(coef_r0_0_31w, coef_r1_0_31w);
      __m512i coef_r0r1_0_31hi = _mm512_unpackhi_epi16(coef_r0_0_31w, coef_r1_0_31w);

      __m512i coef_r2r3_0_31lo = _mm512_unpacklo_epi16(coef_r2_0_31w, coef_r3_0_31w);
      __m512i coef_r2r3_0_31hi = _mm512_unpackhi_epi16(coef_r2_0_31w, coef_r3_0_31w);

      __m512i coef_r4r5_0_31lo = _mm512_unpacklo_epi16(coef_r4_0_31w, coef_r5_0_31w);
      __m512i coef_r4r5_0_31hi = _mm512_unpackhi_epi16(coef_r4_0_31w, coef_r5_0_31w);

      __m512i coef_r6r7_0_31lo = _mm512_unpacklo_epi16(coef_r6_0_31w, coef_r7_0_31w);
      __m512i coef_r6r7_0_31hi = _mm512_unpackhi_epi16(coef_r6_0_31w, coef_r7_0_31w);

      __m512i coef_r8r9_0_31lo = _mm512_unpacklo_epi16(coef_r8_0_31w, coef_r9_0_31w);
      __m512i coef_r8r9_0_31hi = _mm512_unpackhi_epi16(coef_r8_0_31w, coef_r9_0_31w);

      __m512i coef_r10r11_0_31lo = _mm512_unpacklo_epi16(coef_r10_0_31w, coef_r11_0_31w);
      __m512i coef_r10r11_0_31hi = _mm512_unpackhi_epi16(coef_r10_0_31w, coef_r11_0_31w);

      __m512i coef_r12r13_0_31lo = _mm512_unpacklo_epi16(coef_r12_0_31w, coef_r13_0_31w);
      __m512i coef_r12r13_0_31hi = _mm512_unpackhi_epi16(coef_r12_0_31w, coef_r13_0_31w);

      __m512i coef_r14r15_0_31lo = _mm512_unpacklo_epi16(coef_r14_0_31w, coef_r15_0_31w);
      __m512i coef_r14r15_0_31hi = _mm512_unpackhi_epi16(coef_r14_0_31w, coef_r15_0_31w);

      // TODO: store transposed resampling program coeffs to temp buffer for reusage at each line

      // convert resampling program in H-form into permuting indexes for src transposition in V-form
      __m512i perm_0_0_15 = _mm512_loadu_si512((__m512i*)(&program->pixel_offset[x])); // 16 offsets
      __m512i perm_0_16_31 = _mm512_loadu_si512((__m512i*)(&program->pixel_offset[x + 16])); //  16 offsets

      int iStart = program->pixel_offset[x];
      __m512i m512i_Start = _mm512_set1_epi32(iStart);

      perm_0_0_15 = _mm512_sub_epi32(perm_0_0_15, m512i_Start);
      perm_0_16_31 = _mm512_sub_epi32(perm_0_16_31, m512i_Start);

      __m256i m256i_perm_0_0_15 = _mm512_cvtepi32_epi16(perm_0_0_15);
      __m256i m256i_perm_0_16_31 = _mm512_cvtepi32_epi16(perm_0_16_31);

      __m128i mm128i_perm_0_0_15 = _mm256_cvtepi16_epi8(m256i_perm_0_0_15);
      __m128i mm128i_perm_0_16_31 = _mm256_cvtepi16_epi8(m256i_perm_0_16_31);

      __m512i perm_0 = _mm512_inserti32x4(_mm512_zextsi128_si512(mm128i_perm_0_0_15), mm128i_perm_0_16_31, 1); // Lane 0+1 - 32 offsets only

      // Taps are contiguous (0, 1, 2, 3), so we increment perm indexes by 1 (in pairs by 2 each).
      const __m512i two_epi8 = _mm512_set1_epi8(2);

      uint8_t* AVS_RESTRICT dst_ptr = dst8 + x + y_from * dst_pitch;
      const uint8_t* src_ptr = src8 + iStart + y_from * src_pitch; // all permute offsets relative to this start offset

      // Calculate remaining pixels for bounds checking in partial_load mode. 1..128 remaining pixels possible.
      const int remaining = program->source_size - iStart;
	  // _bzhi_u64 creates a mask with the lower N bits set. If N >= 64, it returns all ones (~0ULL).
#ifndef _M_X64
      const __mmask64 k1 = make_low_mask64(remaining);
      const __mmask64 k2 = make_low_mask64(remaining - 64);
#else
      const __mmask64 k1 = _bzhi_u64(~0ULL, remaining);
      const __mmask64 k2 = _bzhi_u64(~0ULL, std::max(0, remaining - 64));
#endif

      for (int y = y_from; y < y_to; y++)
      {
        __m512i data_src, data_src2;

        // working permute indexes for advancing to save number of registers used
        __m512i perm_w_even = perm_0;
        __m512i perm_w_odd = _mm512_add_epi8(perm_0, _mm512_set1_epi8(1));

        if constexpr (partial_load) {
          // Safe masked loads for the image edge
          data_src = _mm512_maskz_loadu_epi8(k1, src_ptr);
          data_src2 = _mm512_maskz_loadu_epi8(k2, src_ptr + 64);
        }
        else {
          // Fast unaligned loads for the safe zone
          data_src = _mm512_loadu_si512(src_ptr);
          data_src2 = _mm512_loadu_si512(src_ptr + 64);
        }

        // rows 0..3
        __m512i data_0 = _mm512_permutex2var_epi8_SIMUL<UseVBMI>(data_src, perm_w_even, data_src2);
        __m512i data_1 = _mm512_permutex2var_epi8_SIMUL<UseVBMI>(data_src, perm_w_odd, data_src2);
        perm_w_even = _mm512_add_epi8(perm_w_even, two_epi8);
        perm_w_odd = _mm512_add_epi8(perm_w_odd, two_epi8);

        __m512i data_2 = _mm512_permutex2var_epi8_SIMUL<UseVBMI>(data_src, perm_w_even, data_src2);
        __m512i data_3 = _mm512_permutex2var_epi8_SIMUL<UseVBMI>(data_src, perm_w_odd, data_src2);
        perm_w_even = _mm512_add_epi8(perm_w_even, two_epi8);
        perm_w_odd = _mm512_add_epi8(perm_w_odd, two_epi8);

        // unpack 8->16bit 4 to 4 512bit registers source
        __m512i src_r0_0_31 = _mm512_cvtepu8_epi16(_mm512_castsi512_si256(data_0));
        __m512i src_r1_0_31 = _mm512_cvtepu8_epi16(_mm512_castsi512_si256(data_1));
        __m512i src_r2_0_31 = _mm512_cvtepu8_epi16(_mm512_castsi512_si256(data_2));
        __m512i src_r3_0_31 = _mm512_cvtepu8_epi16(_mm512_castsi512_si256(data_3));

        // transposition to H-pairs 4 to 4 512bit registers 
        __m512i src_r0r1_0_31lo = _mm512_unpacklo_epi16(src_r0_0_31, src_r1_0_31);
        __m512i src_r0r1_0_31hi = _mm512_unpackhi_epi16(src_r0_0_31, src_r1_0_31);
        __m512i src_r2r3_0_31lo = _mm512_unpacklo_epi16(src_r2_0_31, src_r3_0_31);
        __m512i src_r2r3_0_31hi = _mm512_unpackhi_epi16(src_r2_0_31, src_r3_0_31);

        // making FMA in 32bits accs as in AVX256 V-resize
        __m512i result_0_31lo = _mm512_add_epi32(_mm512_madd_epi16(src_r0r1_0_31lo, coef_r0r1_0_31lo), _mm512_madd_epi16(src_r2r3_0_31lo, coef_r2r3_0_31lo));
        __m512i result_0_31hi = _mm512_add_epi32(_mm512_madd_epi16(src_r0r1_0_31hi, coef_r0r1_0_31hi), _mm512_madd_epi16(src_r2r3_0_31hi, coef_r2r3_0_31hi));

        // rows 4..7
        __m512i data_4 = _mm512_permutex2var_epi8_SIMUL<UseVBMI>(data_src, perm_w_even, data_src2);
        __m512i data_5 = _mm512_permutex2var_epi8_SIMUL<UseVBMI>(data_src, perm_w_odd, data_src2);
        perm_w_even = _mm512_add_epi8(perm_w_even, two_epi8);
        perm_w_odd = _mm512_add_epi8(perm_w_odd, two_epi8);

        __m512i data_6 = _mm512_permutex2var_epi8_SIMUL<UseVBMI>(data_src, perm_w_even, data_src2);
        __m512i data_7 = _mm512_permutex2var_epi8_SIMUL<UseVBMI>(data_src, perm_w_odd, data_src2);
        perm_w_even = _mm512_add_epi8(perm_w_even, two_epi8);
        perm_w_odd = _mm512_add_epi8(perm_w_odd, two_epi8);

        // unpack 8->16bit 4 to 8 512bit registers source
        __m512i src_r4_0_31 = _mm512_cvtepu8_epi16(_mm512_castsi512_si256(data_4));
        __m512i src_r5_0_31 = _mm512_cvtepu8_epi16(_mm512_castsi512_si256(data_5));
        __m512i src_r6_0_31 = _mm512_cvtepu8_epi16(_mm512_castsi512_si256(data_6));
        __m512i src_r7_0_31 = _mm512_cvtepu8_epi16(_mm512_castsi512_si256(data_7));

        // transposition to H-pairs 4 to 4 512bit registers 
        __m512i src_r4r5_0_31lo = _mm512_unpacklo_epi16(src_r4_0_31, src_r5_0_31);
        __m512i src_r4r5_0_31hi = _mm512_unpackhi_epi16(src_r4_0_31, src_r5_0_31);
        __m512i src_r6r7_0_31lo = _mm512_unpacklo_epi16(src_r6_0_31, src_r7_0_31);
        __m512i src_r6r7_0_31hi = _mm512_unpackhi_epi16(src_r6_0_31, src_r7_0_31);

        result_0_31lo = _mm512_add_epi32(_mm512_madd_epi16(src_r4r5_0_31lo, coef_r4r5_0_31lo), result_0_31lo);
        result_0_31hi = _mm512_add_epi32(_mm512_madd_epi16(src_r4r5_0_31hi, coef_r4r5_0_31hi), result_0_31hi);
        result_0_31lo = _mm512_add_epi32(_mm512_madd_epi16(src_r6r7_0_31lo, coef_r6r7_0_31lo), result_0_31lo);
        result_0_31hi = _mm512_add_epi32(_mm512_madd_epi16(src_r6r7_0_31hi, coef_r6r7_0_31hi), result_0_31hi);

        // rows 8..11
        __m512i data_8 = _mm512_permutex2var_epi8_SIMUL<UseVBMI>(data_src, perm_w_even, data_src2);
        __m512i data_9 = _mm512_permutex2var_epi8_SIMUL<UseVBMI>(data_src, perm_w_odd, data_src2);
        perm_w_even = _mm512_add_epi8(perm_w_even, two_epi8);
        perm_w_odd = _mm512_add_epi8(perm_w_odd, two_epi8);

        __m512i data_10 = _mm512_permutex2var_epi8_SIMUL<UseVBMI>(data_src, perm_w_even, data_src2);
        __m512i data_11 = _mm512_permutex2var_epi8_SIMUL<UseVBMI>(data_src, perm_w_odd, data_src2);
        perm_w_even = _mm512_add_epi8(perm_w_even, two_epi8);
        perm_w_odd = _mm512_add_epi8(perm_w_odd, two_epi8);

        // unpack 8->16bit 4 to 8 512bit registers source
        __m512i src_r8_0_31 = _mm512_cvtepu8_epi16(_mm512_castsi512_si256(data_8));
        __m512i src_r9_0_31 = _mm512_cvtepu8_epi16(_mm512_castsi512_si256(data_9));
        __m512i src_r10_0_31 = _mm512_cvtepu8_epi16(_mm512_castsi512_si256(data_10));
        __m512i src_r11_0_31 = _mm512_cvtepu8_epi16(_mm512_castsi512_si256(data_11));

        // transposition to H-pairs 4 to 4 512bit registers 
        __m512i src_r8r9_0_31lo = _mm512_unpacklo_epi16(src_r8_0_31, src_r9_0_31);
        __m512i src_r8r9_0_31hi = _mm512_unpackhi_epi16(src_r8_0_31, src_r9_0_31);
        __m512i src_r10r11_0_31lo = _mm512_unpacklo_epi16(src_r10_0_31, src_r11_0_31);
        __m512i src_r10r11_0_31hi = _mm512_unpackhi_epi16(src_r10_0_31, src_r11_0_31);

        result_0_31lo = _mm512_add_epi32(_mm512_madd_epi16(src_r8r9_0_31lo, coef_r8r9_0_31lo), result_0_31lo);
        result_0_31hi = _mm512_add_epi32(_mm512_madd_epi16(src_r8r9_0_31hi, coef_r8r9_0_31hi), result_0_31hi);
        result_0_31lo = _mm512_add_epi32(_mm512_madd_epi16(src_r10r11_0_31lo, coef_r10r11_0_31lo), result_0_31lo);
        result_0_31hi = _mm512_add_epi32(_mm512_madd_epi16(src_r10r11_0_31hi, coef_r10r11_0_31hi), result_0_31hi);

        // rows 12..15
        __m512i data_12 = _mm512_permutex2var_epi8_SIMUL<UseVBMI>(data_src, perm_w_even, data_src2);
        __m512i data_13 = _mm512_permutex2var_epi8_SIMUL<UseVBMI>(data_src, perm_w_odd, data_src2);
        perm_w_even = _mm512_add_epi8(perm_w_even, two_epi8);
        perm_w_odd = _mm512_add_epi8(perm_w_odd, two_epi8);

        __m512i data_14 = _mm512_permutex2var_epi8_SIMUL<UseVBMI>(data_src, perm_w_even, data_src2);
        __m512i data_15 = _mm512_permutex2var_epi8_SIMUL<UseVBMI>(data_src, perm_w_odd, data_src2);

        // unpack 8->16bit 4 to 8 512bit registers source
        __m512i src_r12_0_31 = _mm512_cvtepu8_epi16(_mm512_castsi512_si256(data_12));
        __m512i src_r13_0_31 = _mm512_cvtepu8_epi16(_mm512_castsi512_si256(data_13));
        __m512i src_r14_0_31 = _mm512_cvtepu8_epi16(_mm512_castsi512_si256(data_14));
        __m512i src_r15_0_31 = _mm512_cvtepu8_epi16(_mm512_castsi512_si256(data_15));

        // transposition to H-pairs 4 to 4 512bit registers 
        __m512i src_r12r13_0_31lo = _mm512_unpacklo_epi16(src_r12_0_31, src_r13_0_31);
        __m512i src_r12r13_0_31hi = _mm512_unpackhi_epi16(src_r12_0_31, src_r13_0_31);
        __m512i src_r14r15_0_31lo = _mm512_unpacklo_epi16(src_r14_0_31, src_r15_0_31);
        __m512i src_r14r15_0_31hi = _mm512_unpackhi_epi16(src_r14_0_31, src_r15_0_31);

        result_0_31lo = _mm512_add_epi32(_mm512_madd_epi16(src_r12r13_0_31lo, coef_r12r13_0_31lo), result_0_31lo);
        result_0_31hi = _mm512_add_epi32(_mm512_madd_epi16(src_r12r13_0_31hi, coef_r12r13_0_31hi), result_0_31hi);
        result_0_31lo = _mm512_add_epi32(_mm512_madd_epi16(src_r14r15_0_31lo, coef_r14r15_0_31lo), result_0_31lo);
        result_0_31hi = _mm512_add_epi32(_mm512_madd_epi16(src_r14r15_0_31hi, coef_r14r15_0_31hi), result_0_31hi);

        // rounding
        result_0_31lo = _mm512_add_epi32(result_0_31lo, rounder);
        result_0_31hi = _mm512_add_epi32(result_0_31hi, rounder);
        // scaling down
        result_0_31lo = _mm512_srai_epi32(result_0_31lo, FPScale8bits);
        result_0_31hi = _mm512_srai_epi32(result_0_31hi, FPScale8bits);

        __m512i result_0_31_int16 = _mm512_packus_epi32(result_0_31lo, result_0_31hi);

        __m256i result_0_31_u8 = _mm512_cvtusepi16_epi8(result_0_31_int16);

		result_0_31_u8 = _mm256_min_epu8(result_0_31_u8, clamp_limit_max);
		result_0_31_u8 = _mm256_max_epu8(result_0_31_u8, clamp_limit_min);

        _mm256_stream_si256(reinterpret_cast<__m256i*>(dst_ptr), result_0_31_u8);

        dst_ptr += dst_pitch;
        src_ptr += src_pitch;
      }

      current_coeff += filter_size * PIXELS_AT_A_TIME;
      };

    // Process the 'safe zone' where direct full unaligned loads are acceptable.
    for (; x < width_safe_mod; x += PIXELS_AT_A_TIME)
    {
      do_h_integer_core(std::false_type{});
    }

    // Process the potentially 'unsafe zone' near the image edge, using safe masked loading.
    for (; x < width; x += PIXELS_AT_A_TIME)
    {
      do_h_integer_core(std::true_type{});
    }
  }
}

// uint8_t "mp" versions

// filter size up to 4
// 64 target uint8_t pixels at a time
// 127-byte source loads (127 uint8_t pixels)
// maximum permute index is 127 for _mm512_maskz_permutex2var_epi8 (uint8_t) 
// more premutex version to create 8->16bit converted and low-hi unpacked registers in single permutex instruction
template<bool UseVNNI>
#if defined(__clang__)
__attribute__((__target__("avx512f,avx512cd,avx512bw,avx512dq,avx512vl,avx512vnni,avx512vbmi,avx512vbmi2,avx512bitalg,avx512vpopcntdq")))
#endif
void resize_h_planar_uint8_avx512_permutex_vstripe_mpz_ks4_internal(BYTE* dst8, const BYTE* src8, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height, int bits_per_pixel,const uint8_t range,const bool mode_YUY2)
{
  AVS_UNUSED(bits_per_pixel);

  const int filter_size = program->filter_size; // aligned, practically the coeff table stride

  constexpr int PIXELS_AT_A_TIME = 64;

  // 'source_overread_beyond_targetx' indicates if the filter kernel can read beyond the target width.
  // we load 2x64 source bytes at a time, so ensure safe overread if needed.
  // Our main loop processes calculates for 64 target pixels at a time.
  // Inside that, we load 128 source bytes (2x64) to be able to permutex from that.
  // This we have to check at each mod-PIXELS_AT_A_TIME boundary, the allowance of 128-byte source load.
  const int width_safe_mod = (program->safelimit_128_pixels_each64th_target.overread_possible ? program->safelimit_128_pixels_each64th_target.source_overread_beyond_targetx : width) / PIXELS_AT_A_TIME * PIXELS_AT_A_TIME;

  // Preconditions:

  // 'target_size_alignment' ensures we can safely access coefficients using offsets like
  // 'filter_size * 15' when processing 16 H pixels at a time
  // 'filter_size * 63' when processing 64 H pixels at a time

  // Ensure that coefficient loading beyond the valid target size is safe for 4x4 float loads.
  // We load 4x16bit coeffs at a time

  const int max_scanlines = program->max_scanlines;

  const int val_min = (range==1) ? 0 : 16;
  const int val_max = ((range==1) || (range==4)) ? 255 : (range==2) ? 235 : 240;

  __m256i clamp_limit_min = _mm256_set1_epi16((short)((val_min << 8)|val_min));
  __m256i clamp_limit_max = (mode_YUY2 && ((range>=2) && (range<=3))) ?
    _mm256_set1_epi16((short)(((int)240 << 8)|235)) : _mm256_set1_epi16((short)((val_max << 8)|val_max));	  

  __m512i rounder = _mm512_set1_epi32(1 << (FPScale8bits - 1));

  // Vertical stripe loop for L2 cache optimization
  for (int y_from = 0; y_from < height; y_from += max_scanlines)
  {
    int y_to = std::min(y_from + max_scanlines, height);

    // Reset current_coeff for the start of the stripe (points to start of row's coeffs)
    const short* AVS_RESTRICT current_coeff = program->pixel_coefficient;

    int x = 0;

    // Lambda to handle both safe (fast) and unsafe (masked/partial) loading paths
    auto do_h_integer_core = [&](auto partial_load) {

      // prepare coefs in transposed V-form
      // TODO: make storage in transposed form, 64 x uint16 transposition looks too slow
      // TO FIX: filter_size=4 resize uses filter_size=16 - lots or RAM/cache wasted, in the ready to use transposed form it will be much more optimized

      // 4coefs of 16bit is 64bits, can be loaded as set_epi64 ?
      __m512i coef_0_7 = _mm512_setr_epi64(
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 0),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 1),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 2),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 3),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 4),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 5),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 6),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 7)
      );

      __m512i coef_8_15 = _mm512_setr_epi64(
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 8),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 9),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 10),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 11),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 12),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 13),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 14),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 15)
      );

      __m512i coef_16_23 = _mm512_setr_epi64(
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 16),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 17),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 18),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 19),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 20),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 21),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 22),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 23)
      );

      __m512i coef_24_31 = _mm512_setr_epi64(
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 24),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 25),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 26),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 27),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 28),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 29),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 30),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 31)
      );

      __m512i coef_32_39 = _mm512_setr_epi64(
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 32),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 33),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 34),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 35),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 36),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 37),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 38),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 39)
      );

      __m512i coef_40_47 = _mm512_setr_epi64(
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 40),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 41),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 42),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 43),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 44),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 45),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 46),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 47)
      );

      __m512i coef_48_55 = _mm512_setr_epi64(
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 48),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 49),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 50),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 51),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 52),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 53),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 54),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 55)
      );

      __m512i coef_56_63 = _mm512_setr_epi64(
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 56),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 57),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 58),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 59),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 60),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 61),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 62),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 63)
      );

      // Transpose with permutex
      __m512i c_perm_0 = _mm512_set_epi16(28 + 32, 24 + 32, 20 + 32, 16 + 32, 12 + 32, 8 + 32, 4 + 32, 0 + 32, 28, 24, 20, 16, 12, 8, 4, 0,
        28 + 32, 24 + 32, 20 + 32, 16 + 32, 12 + 32, 8 + 32, 4 + 32, 0 + 32, 28, 24, 20, 16, 12, 8, 4, 0);
      __m512i one_epi16 = _mm512_set1_epi16(1);
      const __mmask32 k_high = 0xFFFF0000;
      // 0.0 .. 15.0 in low 256, 16.0 .. 31.0 in high 256
      __m512i coef_r0_0_31w = _mm512_mask_blend_epi16(k_high, _mm512_permutex2var_epi16(coef_0_7, c_perm_0, coef_8_15), _mm512_permutex2var_epi16(coef_16_23, c_perm_0, coef_24_31));
      // 32.0 .. 47.0  in low 256, 48.0 .. 63.0 in high 256
      __m512i coef_r0_32_63w = _mm512_mask_blend_epi16(k_high, _mm512_permutex2var_epi16(coef_32_39, c_perm_0, coef_40_47), _mm512_permutex2var_epi16(coef_48_55, c_perm_0, coef_56_63));
      c_perm_0 = _mm512_add_epi16(c_perm_0, one_epi16);
      __m512i coef_r1_0_31w = _mm512_mask_blend_epi16(k_high, _mm512_permutex2var_epi16(coef_0_7, c_perm_0, coef_8_15), _mm512_permutex2var_epi16(coef_16_23, c_perm_0, coef_24_31));
      __m512i coef_r1_32_63w = _mm512_mask_blend_epi16(k_high, _mm512_permutex2var_epi16(coef_32_39, c_perm_0, coef_40_47), _mm512_permutex2var_epi16(coef_48_55, c_perm_0, coef_56_63));
      c_perm_0 = _mm512_add_epi16(c_perm_0, one_epi16);
      __m512i coef_r2_0_31w = _mm512_mask_blend_epi16(k_high, _mm512_permutex2var_epi16(coef_0_7, c_perm_0, coef_8_15), _mm512_permutex2var_epi16(coef_16_23, c_perm_0, coef_24_31));
      __m512i coef_r2_32_63w = _mm512_mask_blend_epi16(k_high, _mm512_permutex2var_epi16(coef_32_39, c_perm_0, coef_40_47), _mm512_permutex2var_epi16(coef_48_55, c_perm_0, coef_56_63));
      c_perm_0 = _mm512_add_epi16(c_perm_0, one_epi16);
      __m512i coef_r3_0_31w = _mm512_mask_blend_epi16(k_high, _mm512_permutex2var_epi16(coef_0_7, c_perm_0, coef_8_15), _mm512_permutex2var_epi16(coef_16_23, c_perm_0, coef_24_31));
      __m512i coef_r3_32_63w = _mm512_mask_blend_epi16(k_high, _mm512_permutex2var_epi16(coef_32_39, c_perm_0, coef_40_47), _mm512_permutex2var_epi16(coef_48_55, c_perm_0, coef_56_63));

      // convert-transpose to H-pairs for madd ? better to do with single permutex in future
      // 8 to 8 512 registers - finally real working coeffs to store in the transposed resampling program for block of 64 target samples
      const __m512i coef_r0r1_0_31lo = _mm512_unpacklo_epi16(coef_r0_0_31w, coef_r1_0_31w);
      const __m512i coef_r0r1_0_31hi = _mm512_unpackhi_epi16(coef_r0_0_31w, coef_r1_0_31w);

      const __m512i coef_r0r1_32_63lo = _mm512_unpacklo_epi16(coef_r0_32_63w, coef_r1_32_63w);
      const __m512i coef_r0r1_32_63hi = _mm512_unpackhi_epi16(coef_r0_32_63w, coef_r1_32_63w);

      const __m512i coef_r2r3_0_31lo = _mm512_unpacklo_epi16(coef_r2_0_31w, coef_r3_0_31w);
      const __m512i coef_r2r3_0_31hi = _mm512_unpackhi_epi16(coef_r2_0_31w, coef_r3_0_31w);

      const __m512i coef_r2r3_32_63lo = _mm512_unpacklo_epi16(coef_r2_32_63w, coef_r3_32_63w);
      const __m512i coef_r2r3_32_63hi = _mm512_unpackhi_epi16(coef_r2_32_63w, coef_r3_32_63w);

      // TODO: store transposed resampling program coeffs to temp buffer for reusage at each line

      // convert resampling program in H-form into permuting indexes for src transposition in V-form
      __m512i perm_0_0_15 = _mm512_loadu_si512((__m512i*)(&program->pixel_offset[x])); // 16 offsets
      __m512i perm_0_16_31 = _mm512_loadu_si512((__m512i*)(&program->pixel_offset[x + 16])); //  16 offsets
      __m512i perm_0_32_47 = _mm512_loadu_si512((__m512i*)(&program->pixel_offset[x + 32])); //  16 offsets
      __m512i perm_0_48_63 = _mm512_loadu_si512((__m512i*)(&program->pixel_offset[x + 48])); //  16 offsets

      int iStart = program->pixel_offset[x];
      __m512i m512i_Start = _mm512_set1_epi32(iStart);

      perm_0_0_15 = _mm512_sub_epi32(perm_0_0_15, m512i_Start);
      perm_0_16_31 = _mm512_sub_epi32(perm_0_16_31, m512i_Start);
      perm_0_32_47 = _mm512_sub_epi32(perm_0_32_47, m512i_Start);
      perm_0_48_63 = _mm512_sub_epi32(perm_0_48_63, m512i_Start);

      __m256i m256i_perm_0_0_15 = _mm512_cvtepi32_epi16(perm_0_0_15);
      __m256i m256i_perm_0_16_31 = _mm512_cvtepi32_epi16(perm_0_16_31);
      __m256i m256i_perm_0_32_47 = _mm512_cvtepi32_epi16(perm_0_32_47);
      __m256i m256i_perm_0_48_63 = _mm512_cvtepi32_epi16(perm_0_48_63);

      __m512i perm_0_0_31 = _mm512_inserti64x4(_mm512_zextsi256_si512(m256i_perm_0_0_15), m256i_perm_0_16_31, 1);
      __m512i perm_0_32_63 = _mm512_inserti64x4(_mm512_zextsi256_si512(m256i_perm_0_32_47), m256i_perm_0_48_63, 1);

      // Taps are contiguous (0, 1, 2, 3), so we increment perm indexes by 1.
      __m512i perm_1_0_31 = _mm512_add_epi16(perm_0_0_31, one_epi16);
      __m512i perm_2_0_31 = _mm512_add_epi16(perm_1_0_31, one_epi16);
      __m512i perm_3_0_31 = _mm512_add_epi16(perm_2_0_31, one_epi16);

      __m512i perm_1_32_63 = _mm512_add_epi16(perm_0_32_63, one_epi16);
      __m512i perm_2_32_63 = _mm512_add_epi16(perm_1_32_63, one_epi16);
      __m512i perm_3_32_63 = _mm512_add_epi16(perm_2_32_63, one_epi16);

      const __m512i perm_r0r1_0_31lo = _mm512_unpacklo_epi16(perm_0_0_31, perm_1_0_31);
      const __m512i perm_r0r1_0_31hi = _mm512_unpackhi_epi16(perm_0_0_31, perm_1_0_31);

      const __m512i perm_r0r1_32_63lo = _mm512_unpacklo_epi16(perm_0_32_63, perm_1_32_63);
      const __m512i perm_r0r1_32_63hi = _mm512_unpackhi_epi16(perm_0_32_63, perm_1_32_63);

      const __m512i perm_r2r3_0_31lo = _mm512_unpacklo_epi16(perm_2_0_31, perm_3_0_31);
      const __m512i perm_r2r3_0_31hi = _mm512_unpackhi_epi16(perm_2_0_31, perm_3_0_31);

      const __m512i perm_r2r3_32_63lo = _mm512_unpacklo_epi16(perm_2_32_63, perm_3_32_63);
      const __m512i perm_r2r3_32_63hi = _mm512_unpackhi_epi16(perm_2_32_63, perm_3_32_63);

      uint8_t* AVS_RESTRICT dst_ptr = dst8 + x + y_from * dst_pitch;
      const uint8_t* src_ptr = src8 + iStart + y_from * src_pitch; // all permute offsets relative to this start offset

      // Calculate remaining pixels for bounds checking in partial_load mode. 1..128 remaining pixels possible.
      const int remaining = program->source_size - iStart;
	  // _bzhi_u64 creates a mask with the lower N bits set. If N >= 64, it returns all ones (~0ULL).
#ifndef _M_X64
      const __mmask64 k1 = make_low_mask64(remaining);
      const __mmask64 k2 = make_low_mask64(remaining - 64);
#else
      const __mmask64 k1 = _bzhi_u64(~0ULL, remaining);
      const __mmask64 k2 = _bzhi_u64(~0ULL, std::max(0, remaining - 64));
#endif

      // mask is used to zero out every odd byte, so that the result of the permute is a
      // vector of zero-extended 8-bit values in 16-bit lanes, preparing 8-bit data for 16-bit FMA operations in AVX512.
      const __mmask64 k_zh8 = 0x5555555555555555ULL;

      for (int y = y_from; y < y_to; y++)
      {
        // 8 coeffs + 8 permute_idx + 2 src + 4 temporal + 1 rounder ~= 23 regs (permute2var overwrite first source - really may be more needed)
        __m512i data_src, data_src2;

        if constexpr (partial_load) {
          // Safe masked loads for the image edge
          data_src = _mm512_maskz_loadu_epi8(k1, src_ptr);
          data_src2 = _mm512_maskz_loadu_epi8(k2, src_ptr + 64);
        }
        else {
          // Fast unaligned loads for the safe zone
          data_src = _mm512_loadu_si512(src_ptr);
          data_src2 = _mm512_loadu_si512(src_ptr + 64);
        }

        __m512i src_r0r1_0_31lo = _mm512_maskz_permutex2var_epi8_SIMUL<UseVNNI>(k_zh8, data_src, perm_r0r1_0_31lo, data_src2);
        __m512i src_r0r1_0_31hi = _mm512_maskz_permutex2var_epi8_SIMUL<UseVNNI>(k_zh8, data_src, perm_r0r1_0_31hi, data_src2);

        __m512i src_r0r1_32_63lo = _mm512_maskz_permutex2var_epi8_SIMUL<UseVNNI>(k_zh8, data_src, perm_r0r1_32_63lo, data_src2);
        __m512i src_r0r1_32_63hi = _mm512_maskz_permutex2var_epi8_SIMUL<UseVNNI>(k_zh8, data_src, perm_r0r1_32_63hi, data_src2);

        __m512i src_r2r3_0_31lo = _mm512_maskz_permutex2var_epi8_SIMUL<UseVNNI>(k_zh8, data_src, perm_r2r3_0_31lo, data_src2);
        __m512i src_r2r3_0_31hi = _mm512_maskz_permutex2var_epi8_SIMUL<UseVNNI>(k_zh8, data_src, perm_r2r3_0_31hi, data_src2);

        __m512i src_r2r3_32_63lo = _mm512_maskz_permutex2var_epi8_SIMUL<UseVNNI>(k_zh8, data_src, perm_r2r3_32_63lo, data_src2);
        __m512i src_r2r3_32_63hi = _mm512_maskz_permutex2var_epi8_SIMUL<UseVNNI>(k_zh8, data_src, perm_r2r3_32_63hi, data_src2);

        __m512i result_0_31lo, result_0_31hi;
        __m512i result_32_63lo, result_32_63hi;

        if constexpr (UseVNNI)
        {
          // Unlike _ks8, here in _ks4 LLVM did not use VNNI vpdpwssd, but decided to use madd and add. Perhaps it knows port pressure and latency?
          result_0_31lo = _mm512_dpwssd_epi32(rounder, src_r0r1_0_31lo, coef_r0r1_0_31lo);
          result_0_31lo = _mm512_dpwssd_epi32(result_0_31lo, src_r2r3_0_31lo, coef_r2r3_0_31lo);

          result_0_31hi = _mm512_dpwssd_epi32(rounder, src_r0r1_0_31hi, coef_r0r1_0_31hi);
          result_0_31hi = _mm512_dpwssd_epi32(result_0_31hi, src_r2r3_0_31hi, coef_r2r3_0_31hi);

          result_32_63lo = _mm512_dpwssd_epi32(rounder, src_r0r1_32_63lo, coef_r0r1_32_63lo);
          result_32_63lo = _mm512_dpwssd_epi32(result_32_63lo, src_r2r3_32_63lo, coef_r2r3_32_63lo);

          result_32_63hi = _mm512_dpwssd_epi32(rounder, src_r0r1_32_63hi, coef_r0r1_32_63hi);
          result_32_63hi = _mm512_dpwssd_epi32(result_32_63hi, src_r2r3_32_63hi, coef_r2r3_32_63hi);
        }
        else
        {
          // making FMA in 32bits accs as in AVX256 V-resize
          result_0_31lo = _mm512_add_epi32(_mm512_madd_epi16(src_r0r1_0_31lo, coef_r0r1_0_31lo), _mm512_madd_epi16(src_r2r3_0_31lo, coef_r2r3_0_31lo));
          result_0_31hi = _mm512_add_epi32(_mm512_madd_epi16(src_r0r1_0_31hi, coef_r0r1_0_31hi), _mm512_madd_epi16(src_r2r3_0_31hi, coef_r2r3_0_31hi));

          result_32_63lo = _mm512_add_epi32(_mm512_madd_epi16(src_r0r1_32_63lo, coef_r0r1_32_63lo), _mm512_madd_epi16(src_r2r3_32_63lo, coef_r2r3_32_63lo));
          result_32_63hi = _mm512_add_epi32(_mm512_madd_epi16(src_r0r1_32_63hi, coef_r0r1_32_63hi), _mm512_madd_epi16(src_r2r3_32_63hi, coef_r2r3_32_63hi));

          // rounding
          result_0_31lo = _mm512_add_epi32(result_0_31lo, rounder);
          result_0_31hi = _mm512_add_epi32(result_0_31hi, rounder);
          result_32_63lo = _mm512_add_epi32(result_32_63lo, rounder);
          result_32_63hi = _mm512_add_epi32(result_32_63hi, rounder);
        }

        // scaling down
        result_0_31lo = _mm512_srai_epi32(result_0_31lo, FPScale8bits);
        result_0_31hi = _mm512_srai_epi32(result_0_31hi, FPScale8bits);
        result_32_63lo = _mm512_srai_epi32(result_32_63lo, FPScale8bits);
        result_32_63hi = _mm512_srai_epi32(result_32_63hi, FPScale8bits);

        __m512i result_0_31_int16 = _mm512_packus_epi32(result_0_31lo, result_0_31hi);
        __m512i result_32_63_int16 = _mm512_packus_epi32(result_32_63lo, result_32_63hi);

        __m256i result_0_31_u8 = _mm512_cvtusepi16_epi8(result_0_31_int16);
        __m256i result_32_63_u8 = _mm512_cvtusepi16_epi8(result_32_63_int16);
		
		result_0_31_u8 = _mm256_min_epu8(result_0_31_u8, clamp_limit_max);
		result_0_31_u8 = _mm256_max_epu8(result_0_31_u8, clamp_limit_min);
		result_32_63_u8 = _mm256_min_epu8(result_32_63_u8, clamp_limit_max);
		result_32_63_u8 = _mm256_max_epu8(result_32_63_u8, clamp_limit_min);

        // cast is enough, no need to _mm512_zextsi256_si512
        _mm512_stream_si512(reinterpret_cast<__m512i*>(dst_ptr),
          _mm512_inserti64x4(_mm512_castsi256_si512(result_0_31_u8), result_32_63_u8, 1));

        dst_ptr += dst_pitch;
        src_ptr += src_pitch;
      }

      current_coeff += filter_size * PIXELS_AT_A_TIME;
      };

    // Process the 'safe zone' where direct full unaligned loads are acceptable.
    for (; x < width_safe_mod; x += PIXELS_AT_A_TIME)
    {
      do_h_integer_core(std::false_type{});
    }

    // Process the potentially 'unsafe zone' near the image edge, using safe masked loading.
    for (; x < width; x += PIXELS_AT_A_TIME)
    {
      do_h_integer_core(std::true_type{});
    }
  }
}

// filter size up to 8
// 64 target uint8_t pixels at a time
// 128-byte source loads (128 uint8_t pixels)
// maximum permute index is 128 for _mm512_maskz_permutex2var_epi8 (uint8_t)
// support VNNI and madd FMA
template<bool UseVNNI>
#if defined(__clang__)
__attribute__((__target__("avx512f,avx512cd,avx512bw,avx512dq,avx512vl,avx512vnni,avx512vbmi,avx512vbmi2,avx512bitalg,avx512vpopcntdq")))
#endif
void resize_h_planar_uint8_avx512_permutex_vstripe_mpz_ks8_internal(BYTE* dst8, const BYTE* src8, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height, int bits_per_pixel,const uint8_t range,const bool mode_YUY2)
{
  AVS_UNUSED(bits_per_pixel);

  const int filter_size = program->filter_size; // aligned, practically the coeff table stride

  constexpr int PIXELS_AT_A_TIME = 64;

  // 'source_overread_beyond_targetx' indicates if the filter kernel can read beyond the target width.
  // we load 2x64 source bytes at a time, so ensure safe overread if needed.
  // Our main loop processes calculates for 64 target pixels at a time.
  // Inside that, we load 128 source bytes (2x64) to be able to permutex from that.
  // This we have to check at each mod-PIXELS_AT_A_TIME boundary, the allowance of 128-byte source load.
  const int width_safe_mod = (program->safelimit_128_pixels_each64th_target.overread_possible ? program->safelimit_128_pixels_each64th_target.source_overread_beyond_targetx : width) / PIXELS_AT_A_TIME * PIXELS_AT_A_TIME;

  // Preconditions:

  // 'target_size_alignment' ensures we can safely access coefficients using offsets like
  // 'filter_size * 15' when processing 16 H pixels at a time
  // 'filter_size * 63' when processing 64 H pixels at a time


  // Ensure that coefficient loading beyond the valid target size is safe for 4x8 float loads.
  // We load 8x 'short' coeffs at a time
  // Loading is unaligned, but we fill __m128 registers before combining into __m512

  const int max_scanlines = program->max_scanlines;

  const int val_min = (range==1) ? 0 : 16;
  const int val_max = ((range==1) || (range==4)) ? 255 : (range==2) ? 235 : 240;

  __m256i clamp_limit_min = _mm256_set1_epi16((short)((val_min << 8)|val_min));
  __m256i clamp_limit_max = (mode_YUY2 && ((range>=2) && (range<=3))) ?
    _mm256_set1_epi16((short)(((int)240 << 8)|235)) : _mm256_set1_epi16((short)((val_max << 8)|val_max));	

  __m512i rounder = _mm512_set1_epi32(1 << (FPScale8bits - 1));

  // Vertical stripe loop for L2 cache optimization
  for (int y_from = 0; y_from < height; y_from += max_scanlines)
  {
    int y_to = std::min(y_from + max_scanlines, height);

    // Reset current_coeff for the start of the stripe (points to start of row's coeffs)
    const short* AVS_RESTRICT current_coeff = program->pixel_coefficient;

    int x = 0;

    // Lambda to handle both safe (fast) and unsafe (masked/partial) loading paths
    auto do_h_integer_core = [&](auto partial_load) {

      // prepare coefs in transposed V-form
      // TODO: make storage in transposed form, 64 x uint16 transposition looks too slow

      // 8coefs of 16bit is 128bits 
      __m512i coef_0_3 = _mm512i_loadu_4_m128i(
        (__m128i*)(current_coeff + filter_size * 0), (__m128i*)(current_coeff + filter_size * 1), (__m128i*)(current_coeff + filter_size * 2), (__m128i*)(current_coeff + filter_size * 3));
      __m512i coef_4_7 = _mm512i_loadu_4_m128i(
        (__m128i*)(current_coeff + filter_size * 4), (__m128i*)(current_coeff + filter_size * 5), (__m128i*)(current_coeff + filter_size * 6), (__m128i*)(current_coeff + filter_size * 7));
      __m512i coef_8_11 = _mm512i_loadu_4_m128i(
        (__m128i*)(current_coeff + filter_size * 8), (__m128i*)(current_coeff + filter_size * 9), (__m128i*)(current_coeff + filter_size * 10), (__m128i*)(current_coeff + filter_size * 11));
      __m512i coef_12_15 = _mm512i_loadu_4_m128i(
        (__m128i*)(current_coeff + filter_size * 12), (__m128i*)(current_coeff + filter_size * 13), (__m128i*)(current_coeff + filter_size * 14), (__m128i*)(current_coeff + filter_size * 15));
      __m512i coef_16_19 = _mm512i_loadu_4_m128i(
        (__m128i*)(current_coeff + filter_size * 16), (__m128i*)(current_coeff + filter_size * 17), (__m128i*)(current_coeff + filter_size * 18), (__m128i*)(current_coeff + filter_size * 19));
      __m512i coef_20_23 = _mm512i_loadu_4_m128i(
        (__m128i*)(current_coeff + filter_size * 20), (__m128i*)(current_coeff + filter_size * 21), (__m128i*)(current_coeff + filter_size * 22), (__m128i*)(current_coeff + filter_size * 23));
      __m512i coef_24_27 = _mm512i_loadu_4_m128i(
        (__m128i*)(current_coeff + filter_size * 24), (__m128i*)(current_coeff + filter_size * 25), (__m128i*)(current_coeff + filter_size * 26), (__m128i*)(current_coeff + filter_size * 27));
      __m512i coef_28_31 = _mm512i_loadu_4_m128i(
        (__m128i*)(current_coeff + filter_size * 28), (__m128i*)(current_coeff + filter_size * 29), (__m128i*)(current_coeff + filter_size * 30), (__m128i*)(current_coeff + filter_size * 31));

      __m512i coef_32_35 = _mm512i_loadu_4_m128i(
        (__m128i*)(current_coeff + filter_size * 32), (__m128i*)(current_coeff + filter_size * 33), (__m128i*)(current_coeff + filter_size * 34), (__m128i*)(current_coeff + filter_size * 35));
      __m512i coef_36_39 = _mm512i_loadu_4_m128i(
        (__m128i*)(current_coeff + filter_size * 36), (__m128i*)(current_coeff + filter_size * 37), (__m128i*)(current_coeff + filter_size * 38), (__m128i*)(current_coeff + filter_size * 39));
      __m512i coef_40_43 = _mm512i_loadu_4_m128i(
        (__m128i*)(current_coeff + filter_size * 40), (__m128i*)(current_coeff + filter_size * 41), (__m128i*)(current_coeff + filter_size * 42), (__m128i*)(current_coeff + filter_size * 43));
      __m512i coef_44_47 = _mm512i_loadu_4_m128i(
        (__m128i*)(current_coeff + filter_size * 44), (__m128i*)(current_coeff + filter_size * 45), (__m128i*)(current_coeff + filter_size * 46), (__m128i*)(current_coeff + filter_size * 47));
      __m512i coef_48_51 = _mm512i_loadu_4_m128i(
        (__m128i*)(current_coeff + filter_size * 48), (__m128i*)(current_coeff + filter_size * 49), (__m128i*)(current_coeff + filter_size * 50), (__m128i*)(current_coeff + filter_size * 51));
      __m512i coef_52_55 = _mm512i_loadu_4_m128i(
        (__m128i*)(current_coeff + filter_size * 52), (__m128i*)(current_coeff + filter_size * 53), (__m128i*)(current_coeff + filter_size * 54), (__m128i*)(current_coeff + filter_size * 55));
      __m512i coef_56_59 = _mm512i_loadu_4_m128i(
        (__m128i*)(current_coeff + filter_size * 56), (__m128i*)(current_coeff + filter_size * 57), (__m128i*)(current_coeff + filter_size * 58), (__m128i*)(current_coeff + filter_size * 59));
      __m512i coef_60_63 = _mm512i_loadu_4_m128i(
        (__m128i*)(current_coeff + filter_size * 60), (__m128i*)(current_coeff + filter_size * 61), (__m128i*)(current_coeff + filter_size * 62), (__m128i*)(current_coeff + filter_size * 63));

      // Transpose with permutex
      __m512i c_perm_0_7 = _mm512_set_epi16(
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        24 + 32, 16 + 32, 8 + 32, 0 + 32, 24, 16, 8, 0);
      __m512i c_perm_8_15 = _mm512_set_epi16(
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        24 + 32, 16 + 32, 8 + 32, 0 + 32, 24, 16, 8, 0,
        0, 0, 0, 0, 0, 0, 0, 0);
      __m512i c_perm_16_23 = _mm512_set_epi16(
        0, 0, 0, 0, 0, 0, 0, 0,
        24 + 32, 16 + 32, 8 + 32, 0 + 32, 24, 16, 8, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0);
      __m512i c_perm_24_31 = _mm512_set_epi16(
        24 + 32, 16 + 32, 8 + 32, 0 + 32, 24, 16, 8, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0);

      __m512i one_epi16 = _mm512_set1_epi16(1);
      const __mmask32 k_8_15 = 0x0000FF00;
      const __mmask32 k_16_23 = 0x00FF0000;
      const __mmask32 k_24_31 = 0xFF000000;

      auto inc_perms = [&](
        __m512i& c0_7,
        __m512i& c8_15,
        __m512i& c16_23,
        __m512i& c24_31
        ) {
          c0_7 = _mm512_add_epi16(c0_7, one_epi16);
          c8_15 = _mm512_add_epi16(c8_15, one_epi16);
          c16_23 = _mm512_add_epi16(c16_23, one_epi16);
          c24_31 = _mm512_add_epi16(c24_31, one_epi16);
        };

      auto make_row_0_63 = [&](
        __m512i& row_0_31w, __m512i& row_32_63w,
        __m512i c0_7, __m512i c8_15, __m512i c16_23, __m512i c24_31
        ) {
          // 0..31
          row_0_31w = _mm512_mask_blend_epi16(
            k_8_15,
            _mm512_permutex2var_epi16(coef_0_3, c0_7, coef_4_7),
            _mm512_permutex2var_epi16(coef_8_11, c8_15, coef_12_15)
          );
          row_0_31w = _mm512_mask_blend_epi16(
            k_16_23,
            row_0_31w,
            _mm512_permutex2var_epi16(coef_16_19, c16_23, coef_20_23)
          );
          row_0_31w = _mm512_mask_blend_epi16(
            k_24_31,
            row_0_31w,
            _mm512_permutex2var_epi16(coef_24_27, c24_31, coef_28_31)
          );

          // 32..63
          row_32_63w = _mm512_mask_blend_epi16(
            k_8_15,
            _mm512_permutex2var_epi16(coef_32_35, c0_7, coef_36_39),
            _mm512_permutex2var_epi16(coef_40_43, c8_15, coef_44_47)
          );
          row_32_63w = _mm512_mask_blend_epi16(
            k_16_23,
            row_32_63w,
            _mm512_permutex2var_epi16(coef_48_51, c16_23, coef_52_55)
          );
          row_32_63w = _mm512_mask_blend_epi16(
            k_24_31,
            row_32_63w,
            _mm512_permutex2var_epi16(coef_56_59, c24_31, coef_60_63)
          );
        };

      __m512i coef_r0_0_31w, coef_r0_32_63w;
      __m512i coef_r1_0_31w, coef_r1_32_63w;
      __m512i coef_r2_0_31w, coef_r2_32_63w;
      __m512i coef_r3_0_31w, coef_r3_32_63w;
      __m512i coef_r4_0_31w, coef_r4_32_63w;
      __m512i coef_r5_0_31w, coef_r5_32_63w;
      __m512i coef_r6_0_31w, coef_r6_32_63w;
      __m512i coef_r7_0_31w, coef_r7_32_63w;

      // r0
      make_row_0_63(coef_r0_0_31w, coef_r0_32_63w, c_perm_0_7, c_perm_8_15, c_perm_16_23, c_perm_24_31);
      inc_perms(c_perm_0_7, c_perm_8_15, c_perm_16_23, c_perm_24_31);
      // r1
      make_row_0_63(coef_r1_0_31w, coef_r1_32_63w, c_perm_0_7, c_perm_8_15, c_perm_16_23, c_perm_24_31);
      inc_perms(c_perm_0_7, c_perm_8_15, c_perm_16_23, c_perm_24_31);
      // r2
      make_row_0_63(coef_r2_0_31w, coef_r2_32_63w, c_perm_0_7, c_perm_8_15, c_perm_16_23, c_perm_24_31);
      inc_perms(c_perm_0_7, c_perm_8_15, c_perm_16_23, c_perm_24_31);
      // r3
      make_row_0_63(coef_r3_0_31w, coef_r3_32_63w, c_perm_0_7, c_perm_8_15, c_perm_16_23, c_perm_24_31);
      inc_perms(c_perm_0_7, c_perm_8_15, c_perm_16_23, c_perm_24_31);
      // r4
      make_row_0_63(coef_r4_0_31w, coef_r4_32_63w, c_perm_0_7, c_perm_8_15, c_perm_16_23, c_perm_24_31);
      inc_perms(c_perm_0_7, c_perm_8_15, c_perm_16_23, c_perm_24_31);
      // r5
      make_row_0_63(coef_r5_0_31w, coef_r5_32_63w, c_perm_0_7, c_perm_8_15, c_perm_16_23, c_perm_24_31);
      inc_perms(c_perm_0_7, c_perm_8_15, c_perm_16_23, c_perm_24_31);
      // r6
      make_row_0_63(coef_r6_0_31w, coef_r6_32_63w, c_perm_0_7, c_perm_8_15, c_perm_16_23, c_perm_24_31);
      inc_perms(c_perm_0_7, c_perm_8_15, c_perm_16_23, c_perm_24_31);
      // r7
      make_row_0_63(coef_r7_0_31w, coef_r7_32_63w, c_perm_0_7, c_perm_8_15, c_perm_16_23, c_perm_24_31);
      /* // last one, not needed
      inc_perms(c_perm_0_7, c_perm_8_15, c_perm_16_23, c_perm_24_31);
      */

      // convert-transpose to H-pairs for madd ? better to do with single permutex in future
      // 16 to 16 512 registers - finally real working coeffs to store in the transposed resampling program for block of 64 target samples
      __m512i coef_r0r1_0_31lo = _mm512_unpacklo_epi16(coef_r0_0_31w, coef_r1_0_31w);
      __m512i coef_r0r1_0_31hi = _mm512_unpackhi_epi16(coef_r0_0_31w, coef_r1_0_31w);

      __m512i coef_r0r1_32_63lo = _mm512_unpacklo_epi16(coef_r0_32_63w, coef_r1_32_63w);
      __m512i coef_r0r1_32_63hi = _mm512_unpackhi_epi16(coef_r0_32_63w, coef_r1_32_63w);

      __m512i coef_r2r3_0_31lo = _mm512_unpacklo_epi16(coef_r2_0_31w, coef_r3_0_31w);
      __m512i coef_r2r3_0_31hi = _mm512_unpackhi_epi16(coef_r2_0_31w, coef_r3_0_31w);

      __m512i coef_r2r3_32_63lo = _mm512_unpacklo_epi16(coef_r2_32_63w, coef_r3_32_63w);
      __m512i coef_r2r3_32_63hi = _mm512_unpackhi_epi16(coef_r2_32_63w, coef_r3_32_63w);

      __m512i coef_r4r5_0_31lo = _mm512_unpacklo_epi16(coef_r4_0_31w, coef_r5_0_31w);
      __m512i coef_r4r5_0_31hi = _mm512_unpackhi_epi16(coef_r4_0_31w, coef_r5_0_31w);

      __m512i coef_r4r5_32_63lo = _mm512_unpacklo_epi16(coef_r4_32_63w, coef_r5_32_63w);
      __m512i coef_r4r5_32_63hi = _mm512_unpackhi_epi16(coef_r4_32_63w, coef_r5_32_63w);

      __m512i coef_r6r7_0_31lo = _mm512_unpacklo_epi16(coef_r6_0_31w, coef_r7_0_31w);
      __m512i coef_r6r7_0_31hi = _mm512_unpackhi_epi16(coef_r6_0_31w, coef_r7_0_31w);

      __m512i coef_r6r7_32_63lo = _mm512_unpacklo_epi16(coef_r6_32_63w, coef_r7_32_63w);
      __m512i coef_r6r7_32_63hi = _mm512_unpackhi_epi16(coef_r6_32_63w, coef_r7_32_63w);

      // TODO: store transposed resampling program coeffs to temp buffer for reusage at each line

      // convert resampling program in H-form into permuting indexes for src transposition in V-form
      __m512i perm_0_0_15 = _mm512_loadu_si512((__m512i*)(&program->pixel_offset[x])); // 16 offsets
      __m512i perm_0_16_31 = _mm512_loadu_si512((__m512i*)(&program->pixel_offset[x + 16])); //  16 offsets
      __m512i perm_0_32_47 = _mm512_loadu_si512((__m512i*)(&program->pixel_offset[x + 32])); //  16 offsets
      __m512i perm_0_48_63 = _mm512_loadu_si512((__m512i*)(&program->pixel_offset[x + 48])); //  16 offsets

      int iStart = program->pixel_offset[x];
      __m512i m512i_Start = _mm512_set1_epi32(iStart);

      perm_0_0_15 = _mm512_sub_epi32(perm_0_0_15, m512i_Start);
      perm_0_16_31 = _mm512_sub_epi32(perm_0_16_31, m512i_Start);
      perm_0_32_47 = _mm512_sub_epi32(perm_0_32_47, m512i_Start);
      perm_0_48_63 = _mm512_sub_epi32(perm_0_48_63, m512i_Start);

      __m256i m256i_perm_0_0_15 = _mm512_cvtepi32_epi16(perm_0_0_15);
      __m256i m256i_perm_0_16_31 = _mm512_cvtepi32_epi16(perm_0_16_31);
      __m256i m256i_perm_0_32_47 = _mm512_cvtepi32_epi16(perm_0_32_47);
      __m256i m256i_perm_0_48_63 = _mm512_cvtepi32_epi16(perm_0_48_63);

      __m512i perm_0_0_31 = _mm512_inserti64x4(_mm512_zextsi256_si512(m256i_perm_0_0_15), m256i_perm_0_16_31, 1);
      __m512i perm_0_32_63 = _mm512_inserti64x4(_mm512_zextsi256_si512(m256i_perm_0_32_47), m256i_perm_0_48_63, 1);

      // no add_one to each perm group to save number of registers used in processing loop
      // need only to first pair (r0 and r1)
      __m512i perm_1_0_31 = _mm512_add_epi16(perm_0_0_31, one_epi16);
      __m512i perm_1_32_63 = _mm512_add_epi16(perm_0_32_63, one_epi16);

      const __m512i perm_r0r1_0_31lo = _mm512_unpacklo_epi16(perm_0_0_31, perm_1_0_31);
      const __m512i perm_r0r1_0_31hi = _mm512_unpackhi_epi16(perm_0_0_31, perm_1_0_31);

      const __m512i perm_r0r1_32_63lo = _mm512_unpacklo_epi16(perm_0_32_63, perm_1_32_63);
      const __m512i perm_r0r1_32_63hi = _mm512_unpackhi_epi16(perm_0_32_63, perm_1_32_63);

      // Taps are contiguous (0, 1, 2, 3), so we increment perm indexes by 1 (in pairs by 2 each).
      const __m512i two_epi16 = _mm512_set1_epi16(2);

      uint8_t* AVS_RESTRICT dst_ptr = dst8 + x + y_from * dst_pitch;
      const uint8_t* src_ptr = src8 + iStart + y_from * src_pitch; // all permute offsets relative to this start offset

      // Calculate remaining pixels for bounds checking in partial_load mode. 1..128 remaining pixels possible.
      const int remaining = program->source_size - iStart;
	  // _bzhi_u64 creates a mask with the lower N bits set. If N >= 64, it returns all ones (~0ULL).
#ifndef _M_X64
      const __mmask64 k1 = make_low_mask64(remaining);
      const __mmask64 k2 = make_low_mask64(remaining - 64);
#else
      const __mmask64 k1 = _bzhi_u64(~0ULL, remaining);
      const __mmask64 k2 = _bzhi_u64(~0ULL, std::max(0, remaining - 64));
#endif

      // mask is used to zero out every odd byte, so that the result of the permute is a
      // vector of zero-extended 8-bit values in 16-bit lanes, preparing 8-bit data for 16-bit FMA operations in AVX512.
      const __mmask64 k_zh8 = 0x5555555555555555ULL;

      for (int y = y_from; y < y_to; y++)
      {
        __m512i data_src, data_src2;

        // working permute indexes for advancing to save number of registers used
        __m512i perm_rNrNp1_0_31lo_w = perm_r0r1_0_31lo;
        __m512i perm_rNrNp1_0_31hi_w = perm_r0r1_0_31hi;

        __m512i perm_rNrNp1_32_63lo_w = perm_r0r1_32_63lo;
        __m512i perm_rNrNp1_32_63hi_w = perm_r0r1_32_63hi;

        if constexpr (partial_load) {
          // Safe masked loads for the image edge
          data_src = _mm512_maskz_loadu_epi8(k1, src_ptr);
          data_src2 = _mm512_maskz_loadu_epi8(k2, src_ptr + 64);
        }
        else {
          // Fast unaligned loads for the safe zone
          data_src = _mm512_loadu_si512(src_ptr);
          data_src2 = _mm512_loadu_si512(src_ptr + 64);
        }

        // rows 0..3
        __m512i src_r0r1_0_31lo = _mm512_maskz_permutex2var_epi8_SIMUL<UseVNNI>(k_zh8, data_src, perm_rNrNp1_0_31lo_w, data_src2);
        __m512i src_r0r1_0_31hi = _mm512_maskz_permutex2var_epi8_SIMUL<UseVNNI>(k_zh8, data_src, perm_rNrNp1_0_31hi_w, data_src2);

        __m512i src_r0r1_32_63lo = _mm512_maskz_permutex2var_epi8_SIMUL<UseVNNI>(k_zh8, data_src, perm_rNrNp1_32_63lo_w, data_src2);
        __m512i src_r0r1_32_63hi = _mm512_maskz_permutex2var_epi8_SIMUL<UseVNNI>(k_zh8, data_src, perm_rNrNp1_32_63hi_w, data_src2);

        // for r2r3
        perm_rNrNp1_0_31lo_w = _mm512_add_epi16(perm_rNrNp1_0_31lo_w, two_epi16);
        perm_rNrNp1_0_31hi_w = _mm512_add_epi16(perm_rNrNp1_0_31hi_w, two_epi16);
        perm_rNrNp1_32_63lo_w = _mm512_add_epi16(perm_rNrNp1_32_63lo_w, two_epi16);
        perm_rNrNp1_32_63hi_w = _mm512_add_epi16(perm_rNrNp1_32_63hi_w, two_epi16);

        __m512i src_r2r3_0_31lo = _mm512_maskz_permutex2var_epi8_SIMUL<UseVNNI>(k_zh8, data_src, perm_rNrNp1_0_31lo_w, data_src2);
        __m512i src_r2r3_0_31hi = _mm512_maskz_permutex2var_epi8_SIMUL<UseVNNI>(k_zh8, data_src, perm_rNrNp1_0_31hi_w, data_src2);

        __m512i src_r2r3_32_63lo = _mm512_maskz_permutex2var_epi8_SIMUL<UseVNNI>(k_zh8, data_src, perm_rNrNp1_32_63lo_w, data_src2);
        __m512i src_r2r3_32_63hi = _mm512_maskz_permutex2var_epi8_SIMUL<UseVNNI>(k_zh8, data_src, perm_rNrNp1_32_63hi_w, data_src2);

        // for r4r5
        perm_rNrNp1_0_31lo_w = _mm512_add_epi16(perm_rNrNp1_0_31lo_w, two_epi16);
        perm_rNrNp1_0_31hi_w = _mm512_add_epi16(perm_rNrNp1_0_31hi_w, two_epi16);
        perm_rNrNp1_32_63lo_w = _mm512_add_epi16(perm_rNrNp1_32_63lo_w, two_epi16);
        perm_rNrNp1_32_63hi_w = _mm512_add_epi16(perm_rNrNp1_32_63hi_w, two_epi16);

        __m512i result_0_31lo, result_0_31hi;
        __m512i result_32_63lo, result_32_63hi;

        if constexpr (UseVNNI)
        {
          result_0_31lo = _mm512_dpwssd_epi32(rounder, src_r0r1_0_31lo, coef_r0r1_0_31lo);
          result_0_31lo = _mm512_dpwssd_epi32(result_0_31lo, src_r2r3_0_31lo, coef_r2r3_0_31lo);

          result_0_31hi = _mm512_dpwssd_epi32(rounder, src_r0r1_0_31hi, coef_r0r1_0_31hi);
          result_0_31hi = _mm512_dpwssd_epi32(result_0_31hi, src_r2r3_0_31hi, coef_r2r3_0_31hi);

          result_32_63lo = _mm512_dpwssd_epi32(rounder, src_r0r1_32_63lo, coef_r0r1_32_63lo);
          result_32_63lo = _mm512_dpwssd_epi32(result_32_63lo, src_r2r3_32_63lo, coef_r2r3_32_63lo);

          result_32_63hi = _mm512_dpwssd_epi32(rounder, src_r0r1_32_63hi, coef_r0r1_32_63hi);
          result_32_63hi = _mm512_dpwssd_epi32(result_32_63hi, src_r2r3_32_63hi, coef_r2r3_32_63hi);
        }
        else
        {
          // making FMA in 32bits accs as in AVX256 V-resize
          result_0_31lo = _mm512_add_epi32(_mm512_madd_epi16(src_r0r1_0_31lo, coef_r0r1_0_31lo), _mm512_madd_epi16(src_r2r3_0_31lo, coef_r2r3_0_31lo));
          result_0_31hi = _mm512_add_epi32(_mm512_madd_epi16(src_r0r1_0_31hi, coef_r0r1_0_31hi), _mm512_madd_epi16(src_r2r3_0_31hi, coef_r2r3_0_31hi));

          result_32_63lo = _mm512_add_epi32(_mm512_madd_epi16(src_r0r1_32_63lo, coef_r0r1_32_63lo), _mm512_madd_epi16(src_r2r3_32_63lo, coef_r2r3_32_63lo));
          result_32_63hi = _mm512_add_epi32(_mm512_madd_epi16(src_r0r1_32_63hi, coef_r0r1_32_63hi), _mm512_madd_epi16(src_r2r3_32_63hi, coef_r2r3_32_63hi));
        }

        // rows 4..7
        __m512i src_r4r5_0_31lo = _mm512_maskz_permutex2var_epi8_SIMUL<UseVNNI>(k_zh8, data_src, perm_rNrNp1_0_31lo_w, data_src2);
        __m512i src_r4r5_0_31hi = _mm512_maskz_permutex2var_epi8_SIMUL<UseVNNI>(k_zh8, data_src, perm_rNrNp1_0_31hi_w, data_src2);

        __m512i src_r4r5_32_63lo = _mm512_maskz_permutex2var_epi8_SIMUL<UseVNNI>(k_zh8, data_src, perm_rNrNp1_32_63lo_w, data_src2);
        __m512i src_r4r5_32_63hi = _mm512_maskz_permutex2var_epi8_SIMUL<UseVNNI>(k_zh8, data_src, perm_rNrNp1_32_63hi_w, data_src2);

        // for r6r7
        perm_rNrNp1_0_31lo_w = _mm512_add_epi16(perm_rNrNp1_0_31lo_w, two_epi16);
        perm_rNrNp1_0_31hi_w = _mm512_add_epi16(perm_rNrNp1_0_31hi_w, two_epi16);
        perm_rNrNp1_32_63lo_w = _mm512_add_epi16(perm_rNrNp1_32_63lo_w, two_epi16);
        perm_rNrNp1_32_63hi_w = _mm512_add_epi16(perm_rNrNp1_32_63hi_w, two_epi16);

        __m512i src_r6r7_0_31lo = _mm512_maskz_permutex2var_epi8_SIMUL<UseVNNI>(k_zh8, data_src, perm_rNrNp1_0_31lo_w, data_src2);
        __m512i src_r6r7_0_31hi = _mm512_maskz_permutex2var_epi8_SIMUL<UseVNNI>(k_zh8, data_src, perm_rNrNp1_0_31hi_w, data_src2);

        __m512i src_r6r7_32_63lo = _mm512_maskz_permutex2var_epi8_SIMUL<UseVNNI>(k_zh8, data_src, perm_rNrNp1_32_63lo_w, data_src2);
        __m512i src_r6r7_32_63hi = _mm512_maskz_permutex2var_epi8_SIMUL<UseVNNI>(k_zh8, data_src, perm_rNrNp1_32_63hi_w, data_src2);

        if constexpr (UseVNNI)
        {
          result_0_31lo = _mm512_dpwssd_epi32(result_0_31lo, src_r4r5_0_31lo, coef_r4r5_0_31lo);
          result_0_31lo = _mm512_dpwssd_epi32(result_0_31lo, src_r6r7_0_31lo, coef_r6r7_0_31lo);

          result_0_31hi = _mm512_dpwssd_epi32(result_0_31hi, src_r4r5_0_31hi, coef_r4r5_0_31hi);
          result_0_31hi = _mm512_dpwssd_epi32(result_0_31hi, src_r6r7_0_31hi, coef_r6r7_0_31hi);

          result_32_63lo = _mm512_dpwssd_epi32(result_32_63lo, src_r4r5_32_63lo, coef_r4r5_32_63lo);
          result_32_63lo = _mm512_dpwssd_epi32(result_32_63lo, src_r6r7_32_63lo, coef_r6r7_32_63lo);

          result_32_63hi = _mm512_dpwssd_epi32(result_32_63hi, src_r4r5_32_63hi, coef_r4r5_32_63hi);
          result_32_63hi = _mm512_dpwssd_epi32(result_32_63hi, src_r6r7_32_63hi, coef_r6r7_32_63hi);

          // rounding VNNI in first FMA already summed
        }
        else
        {
          result_0_31lo = _mm512_add_epi32(_mm512_madd_epi16(src_r4r5_0_31lo, coef_r4r5_0_31lo), result_0_31lo);
          result_0_31hi = _mm512_add_epi32(_mm512_madd_epi16(src_r4r5_0_31hi, coef_r4r5_0_31hi), result_0_31hi);

          result_32_63lo = _mm512_add_epi32(_mm512_madd_epi16(src_r4r5_32_63lo, coef_r4r5_32_63lo), result_32_63lo);
          result_32_63hi = _mm512_add_epi32(_mm512_madd_epi16(src_r4r5_32_63hi, coef_r4r5_32_63hi), result_32_63hi);

          result_0_31lo = _mm512_add_epi32(_mm512_madd_epi16(src_r6r7_0_31lo, coef_r6r7_0_31lo), result_0_31lo);
          result_0_31hi = _mm512_add_epi32(_mm512_madd_epi16(src_r6r7_0_31hi, coef_r6r7_0_31hi), result_0_31hi);

          result_32_63lo = _mm512_add_epi32(_mm512_madd_epi16(src_r6r7_32_63lo, coef_r6r7_32_63lo), result_32_63lo);
          result_32_63hi = _mm512_add_epi32(_mm512_madd_epi16(src_r6r7_32_63hi, coef_r6r7_32_63hi), result_32_63hi);

          // rounding
          result_0_31lo = _mm512_add_epi32(result_0_31lo, rounder);
          result_0_31hi = _mm512_add_epi32(result_0_31hi, rounder);
          result_32_63lo = _mm512_add_epi32(result_32_63lo, rounder);
          result_32_63hi = _mm512_add_epi32(result_32_63hi, rounder);
        }

        // scaling down
        result_0_31lo = _mm512_srai_epi32(result_0_31lo, FPScale8bits);
        result_0_31hi = _mm512_srai_epi32(result_0_31hi, FPScale8bits);
        result_32_63lo = _mm512_srai_epi32(result_32_63lo, FPScale8bits);
        result_32_63hi = _mm512_srai_epi32(result_32_63hi, FPScale8bits);

        __m512i result_0_31_int16 = _mm512_packus_epi32(result_0_31lo, result_0_31hi);
        __m512i result_32_63_int16 = _mm512_packus_epi32(result_32_63lo, result_32_63hi);

        __m256i result_0_31_u8 = _mm512_cvtusepi16_epi8(result_0_31_int16);
        __m256i result_32_63_u8 = _mm512_cvtusepi16_epi8(result_32_63_int16);

		result_0_31_u8 = _mm256_min_epu8(result_0_31_u8, clamp_limit_max);
		result_0_31_u8 = _mm256_max_epu8(result_0_31_u8, clamp_limit_min);
		result_32_63_u8 = _mm256_min_epu8(result_32_63_u8, clamp_limit_max);
		result_32_63_u8 = _mm256_max_epu8(result_32_63_u8, clamp_limit_min);

        // cast is enough, no need to use zeroextend
        _mm512_stream_si512(reinterpret_cast<__m512i*>(dst_ptr),
          _mm512_inserti64x4(_mm512_castsi256_si512(result_0_31_u8), result_32_63_u8, 1));

        dst_ptr += dst_pitch;
        src_ptr += src_pitch;
      }

      current_coeff += filter_size * PIXELS_AT_A_TIME;
      };

    // Process the 'safe zone' where direct full unaligned loads are acceptable.
    for (; x < width_safe_mod; x += PIXELS_AT_A_TIME)
    {
      do_h_integer_core(std::false_type{});
    }

    // Process the potentially 'unsafe zone' near the image edge, using safe masked loading.
    for (; x < width; x += PIXELS_AT_A_TIME)
    {
      do_h_integer_core(std::true_type{});
    }
  }
}

// filter size up to 16
// 32 target uint8_t pixels at a time
// 128-byte source loads (128 uint8_t pixels)
// maximum permute index is 128 for _mm512_maskz_permutex2var_epi8 (uint8_t)
template<bool UseVNNI>
#if defined(__clang__)
__attribute__((__target__("avx512f,avx512cd,avx512bw,avx512dq,avx512vl,avx512vnni,avx512vbmi,avx512vbmi2,avx512bitalg,avx512vpopcntdq")))
#endif
void resize_h_planar_uint8_avx512_permutex_vstripe_mpz_ks16_internal(BYTE* dst8, const BYTE* src8, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height, int bits_per_pixel,const uint8_t range,const bool mode_YUY2)
{
  AVS_UNUSED(bits_per_pixel);

  const int filter_size = program->filter_size; // aligned, practically the coeff table stride

  constexpr int PIXELS_AT_A_TIME = 32;

  // 'source_overread_beyond_targetx' indicates if the filter kernel can read beyond the target width.
  // we load 2x64 source bytes at a time, so ensure safe overread if needed.
  // Our main loop processes calculates for 64 target pixels at a time.
  // Inside that, we load 128 source bytes (2x64) to be able to permutex from that.
  // This we have to check at each mod-PIXELS_AT_A_TIME boundary, the allowance of 128-byte source load.
  const int width_safe_mod = (program->safelimit_128_pixels_each64th_target.overread_possible ? program->safelimit_128_pixels_each64th_target.source_overread_beyond_targetx : width) / PIXELS_AT_A_TIME * PIXELS_AT_A_TIME;

  // Preconditions:

  // 'target_size_alignment' ensures we can safely access coefficients using offsets like
  // 'filter_size * 15' when processing 16 H pixels at a time
  // 'filter_size * 31' when processing 32 H pixels at a time

  // Ensure that coefficient loading beyond the valid target size is safe 
  // We load 16x 'short' coeffs at a time
  // We aligned_load from coeff positions, 32 bytes boundary needed, which is
  // guaranteed by filter_size_alignment >=16

  const int max_scanlines = program->max_scanlines;

  const int val_min = (range==1) ? 0 : 16;
  const int val_max = ((range==1) || (range==4)) ? 255 : (range==2) ? 235 : 240;

  __m256i clamp_limit_min = _mm256_set1_epi16((short)((val_min << 8)|val_min));
  __m256i clamp_limit_max = (mode_YUY2 && ((range>=2) && (range<=3))) ?
    _mm256_set1_epi16((short)(((int)240 << 8)|235)) : _mm256_set1_epi16((short)((val_max << 8)|val_max));	

  __m512i rounder = _mm512_set1_epi32(1 << (FPScale8bits - 1));

  // Vertical stripe loop for L2 cache optimization
  for (int y_from = 0; y_from < height; y_from += max_scanlines)
  {
    int y_to = std::min(y_from + max_scanlines, height);

    // Reset current_coeff for the start of the stripe (points to start of row's coeffs)
    const short* AVS_RESTRICT current_coeff = program->pixel_coefficient;

    int x = 0;

    // Lambda to handle both safe (fast) and unsafe (masked/partial) loading paths
    auto do_h_integer_core = [&](auto partial_load) {

      // prepare coefs in transposed V-form
      // TODO: make storage in transposed form, 32 x uint8 transposition looks too slow

      // 16coefs of 16bit is 256bits - load as pairs of _m256i
      // hope 16-coeffs blocks are 32-bytes aligned for m256i aligned loads ?
      __m512i coef_0_1 = _mm512i_load_2_m256i((__m256i*)(current_coeff + filter_size * 0), (__m256i*)(current_coeff + filter_size * 1));
      __m512i coef_2_3 = _mm512i_load_2_m256i((__m256i*)(current_coeff + filter_size * 2), (__m256i*)(current_coeff + filter_size * 3));
      __m512i coef_4_5 = _mm512i_load_2_m256i((__m256i*)(current_coeff + filter_size * 4), (__m256i*)(current_coeff + filter_size * 5));
      __m512i coef_6_7 = _mm512i_load_2_m256i((__m256i*)(current_coeff + filter_size * 6), (__m256i*)(current_coeff + filter_size * 7));
      __m512i coef_8_9 = _mm512i_load_2_m256i((__m256i*)(current_coeff + filter_size * 8), (__m256i*)(current_coeff + filter_size * 9));
      __m512i coef_10_11 = _mm512i_load_2_m256i((__m256i*)(current_coeff + filter_size * 10), (__m256i*)(current_coeff + filter_size * 11));
      __m512i coef_12_13 = _mm512i_load_2_m256i((__m256i*)(current_coeff + filter_size * 12), (__m256i*)(current_coeff + filter_size * 13));
      __m512i coef_14_15 = _mm512i_load_2_m256i((__m256i*)(current_coeff + filter_size * 14), (__m256i*)(current_coeff + filter_size * 15));
      __m512i coef_16_17 = _mm512i_load_2_m256i((__m256i*)(current_coeff + filter_size * 16), (__m256i*)(current_coeff + filter_size * 17));
      __m512i coef_18_19 = _mm512i_load_2_m256i((__m256i*)(current_coeff + filter_size * 18), (__m256i*)(current_coeff + filter_size * 19));
      __m512i coef_20_21 = _mm512i_load_2_m256i((__m256i*)(current_coeff + filter_size * 20), (__m256i*)(current_coeff + filter_size * 21));
      __m512i coef_22_23 = _mm512i_load_2_m256i((__m256i*)(current_coeff + filter_size * 22), (__m256i*)(current_coeff + filter_size * 23));
      __m512i coef_24_25 = _mm512i_load_2_m256i((__m256i*)(current_coeff + filter_size * 24), (__m256i*)(current_coeff + filter_size * 25));
      __m512i coef_26_27 = _mm512i_load_2_m256i((__m256i*)(current_coeff + filter_size * 26), (__m256i*)(current_coeff + filter_size * 27));
      __m512i coef_28_29 = _mm512i_load_2_m256i((__m256i*)(current_coeff + filter_size * 28), (__m256i*)(current_coeff + filter_size * 29));
      __m512i coef_30_31 = _mm512i_load_2_m256i((__m256i*)(current_coeff + filter_size * 30), (__m256i*)(current_coeff + filter_size * 31));

      // Transpose with permutex
      __m512i c_perm_0_3 = _mm512_set_epi16(
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        16 + 32, 0 + 32, 16, 0);

      __m512i c_perm_4_7 = _mm512_set_epi16(
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        16 + 32, 0 + 32, 16, 0,
        0, 0, 0, 0);

      __m512i c_perm_8_11 = _mm512_set_epi16(
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        16 + 32, 0 + 32, 16, 0,
        0, 0, 0, 0,
        0, 0, 0, 0);

      __m512i c_perm_12_15 = _mm512_set_epi16(
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        16 + 32, 0 + 32, 16, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0);

      __m512i c_perm_16_19 = _mm512_set_epi16(
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        16 + 32, 0 + 32, 16, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0);

      __m512i c_perm_20_23 = _mm512_set_epi16(
        0, 0, 0, 0,
        0, 0, 0, 0,
        16 + 32, 0 + 32, 16, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0);

      __m512i c_perm_24_27 = _mm512_set_epi16(
        0, 0, 0, 0,
        16 + 32, 0 + 32, 16, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0);

      __m512i c_perm_28_31 = _mm512_set_epi16(
        16 + 32, 0 + 32, 16, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0);


      __m512i one_epi16 = _mm512_set1_epi16(1);

      // Define masks for the 4-word (8-byte, 4x 16-bit word) segments within the 32-word vector.
      // Note: Each AVX-512 word (epi16) lane corresponds to one bit in the __mmask32.
      // The masks target chunks of 4 words (k_X_Y represents bits X through Y inclusive).
      const __mmask32 k_4_7 = 0x000000F0;
      const __mmask32 k_8_11 = 0x00000F00;
      const __mmask32 k_12_15 = 0x0000F000;
      const __mmask32 k_16_19 = 0x000F0000;
      const __mmask32 k_20_23 = 0x00F00000;
      const __mmask32 k_24_27 = 0x0F000000;
      const __mmask32 k_28_31 = 0xF0000000;

      // Helper lambda to increment all eight permutation vectors by one.
      auto inc_perms = [&](
        __m512i& p0_3, __m512i& p4_7, __m512i& p8_11, __m512i& p12_15,
        __m512i& p16_19, __m512i& p20_23, __m512i& p24_27, __m512i& p28_31
        ) {
          p0_3 = _mm512_add_epi16(p0_3, one_epi16);
          p4_7 = _mm512_add_epi16(p4_7, one_epi16);
          p8_11 = _mm512_add_epi16(p8_11, one_epi16);
          p12_15 = _mm512_add_epi16(p12_15, one_epi16);
          p16_19 = _mm512_add_epi16(p16_19, one_epi16);
          p20_23 = _mm512_add_epi16(p20_23, one_epi16);
          p24_27 = _mm512_add_epi16(p24_27, one_epi16);
          p28_31 = _mm512_add_epi16(p28_31, one_epi16);
        };

      // Helper lambda to construct one full 32-word (512-bit) coefficient row.
      // It uses mask blending to merge the results of _mm512_permutex2var_epi16
      // for different 4-word segments, using the current permutation vectors.
      auto make_coef_row = [&](
        __m512i& row_result,
        __m512i p0_3, __m512i p4_7, __m512i p8_11, __m512i p12_15,
        __m512i p16_19, __m512i p20_23, __m512i p24_27, __m512i p28_31
        ) {
          // Start with the first segment (words 0-3) and the fourth segment (words 4-7).
          row_result = _mm512_mask_blend_epi16(
            k_4_7,
            _mm512_permutex2var_epi16(coef_0_1, p0_3, coef_2_3),  // words 0-3 (unmasked)
            _mm512_permutex2var_epi16(coef_4_5, p4_7, coef_6_7)   // words 4-7 (masked)
          );

          // Merge segment 8-11
          row_result = _mm512_mask_blend_epi16(
            k_8_11,
            row_result,
            _mm512_permutex2var_epi16(coef_8_9, p8_11, coef_10_11)
          );

          // Merge segment 12-15
          row_result = _mm512_mask_blend_epi16(
            k_12_15,
            row_result,
            _mm512_permutex2var_epi16(coef_12_13, p12_15, coef_14_15)
          );

          // Merge segment 16-19
          row_result = _mm512_mask_blend_epi16(
            k_16_19,
            row_result,
            _mm512_permutex2var_epi16(coef_16_17, p16_19, coef_18_19)
          );

          // Merge segment 20-23
          row_result = _mm512_mask_blend_epi16(
            k_20_23,
            row_result,
            _mm512_permutex2var_epi16(coef_20_21, p20_23, coef_22_23)
          );

          // Merge segment 24-27
          row_result = _mm512_mask_blend_epi16(
            k_24_27,
            row_result,
            _mm512_permutex2var_epi16(coef_24_25, p24_27, coef_26_27)
          );

          // Merge segment 28-31
          row_result = _mm512_mask_blend_epi16(
            k_28_31,
            row_result,
            _mm512_permutex2var_epi16(coef_28_29, p28_31, coef_30_31)
          );
        };

      // Declare the coefficient row variables
      __m512i coef_r0_0_31w, coef_r1_0_31w, coef_r2_0_31w, coef_r3_0_31w;
      __m512i coef_r4_0_31w, coef_r5_0_31w, coef_r6_0_31w, coef_r7_0_31w;
      __m512i coef_r8_0_31w, coef_r9_0_31w, coef_r10_0_31w, coef_r11_0_31w;
      __m512i coef_r12_0_31w, coef_r13_0_31w, coef_r14_0_31w, coef_r15_0_31w;

      // Process Row 0
      make_coef_row(coef_r0_0_31w, c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);
      inc_perms(c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);

      // Process Row 1
      make_coef_row(coef_r1_0_31w, c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);
      inc_perms(c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);

      // Process Row 2
      make_coef_row(coef_r2_0_31w, c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);
      inc_perms(c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);

      // Process Row 3
      make_coef_row(coef_r3_0_31w, c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);
      inc_perms(c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);

      // Process Row 4
      make_coef_row(coef_r4_0_31w, c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);
      inc_perms(c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);

      // Process Row 5
      make_coef_row(coef_r5_0_31w, c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);
      inc_perms(c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);

      // Process Row 6
      make_coef_row(coef_r6_0_31w, c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);
      inc_perms(c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);

      // Process Row 7
      make_coef_row(coef_r7_0_31w, c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);
      inc_perms(c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);

      // Process Row 8
      make_coef_row(coef_r8_0_31w, c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);
      inc_perms(c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);

      // Process Row 9
      make_coef_row(coef_r9_0_31w, c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);
      inc_perms(c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);

      // Process Row 10
      make_coef_row(coef_r10_0_31w, c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);
      inc_perms(c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);

      // Process Row 11
      make_coef_row(coef_r11_0_31w, c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);
      inc_perms(c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);

      // Process Row 12
      make_coef_row(coef_r12_0_31w, c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);
      inc_perms(c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);

      // Process Row 13
      make_coef_row(coef_r13_0_31w, c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);
      inc_perms(c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);

      // Process Row 14
      make_coef_row(coef_r14_0_31w, c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);
      inc_perms(c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);

      // Process Row 15
      make_coef_row(coef_r15_0_31w, c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);
      // No inc_perms here

      // convert-transpose to H-pairs for madd ? better to do with single permutex in future
      // 16 to 16 512 registers - finally real working coeffs to store in the transposed resampling program for block of 64 target samples
      __m512i coef_r0r1_0_31lo = _mm512_unpacklo_epi16(coef_r0_0_31w, coef_r1_0_31w);
      __m512i coef_r0r1_0_31hi = _mm512_unpackhi_epi16(coef_r0_0_31w, coef_r1_0_31w);

      __m512i coef_r2r3_0_31lo = _mm512_unpacklo_epi16(coef_r2_0_31w, coef_r3_0_31w);
      __m512i coef_r2r3_0_31hi = _mm512_unpackhi_epi16(coef_r2_0_31w, coef_r3_0_31w);

      __m512i coef_r4r5_0_31lo = _mm512_unpacklo_epi16(coef_r4_0_31w, coef_r5_0_31w);
      __m512i coef_r4r5_0_31hi = _mm512_unpackhi_epi16(coef_r4_0_31w, coef_r5_0_31w);

      __m512i coef_r6r7_0_31lo = _mm512_unpacklo_epi16(coef_r6_0_31w, coef_r7_0_31w);
      __m512i coef_r6r7_0_31hi = _mm512_unpackhi_epi16(coef_r6_0_31w, coef_r7_0_31w);

      __m512i coef_r8r9_0_31lo = _mm512_unpacklo_epi16(coef_r8_0_31w, coef_r9_0_31w);
      __m512i coef_r8r9_0_31hi = _mm512_unpackhi_epi16(coef_r8_0_31w, coef_r9_0_31w);

      __m512i coef_r10r11_0_31lo = _mm512_unpacklo_epi16(coef_r10_0_31w, coef_r11_0_31w);
      __m512i coef_r10r11_0_31hi = _mm512_unpackhi_epi16(coef_r10_0_31w, coef_r11_0_31w);

      __m512i coef_r12r13_0_31lo = _mm512_unpacklo_epi16(coef_r12_0_31w, coef_r13_0_31w);
      __m512i coef_r12r13_0_31hi = _mm512_unpackhi_epi16(coef_r12_0_31w, coef_r13_0_31w);

      __m512i coef_r14r15_0_31lo = _mm512_unpacklo_epi16(coef_r14_0_31w, coef_r15_0_31w);
      __m512i coef_r14r15_0_31hi = _mm512_unpackhi_epi16(coef_r14_0_31w, coef_r15_0_31w);

      // TODO: store transposed resampling program coeffs to temp buffer for reusage at each line

      // convert resampling program in H-form into permuting indexes for src transposition in V-form
      __m512i perm_0_0_15 = _mm512_loadu_si512((__m512i*)(&program->pixel_offset[x])); // 16 offsets
      __m512i perm_0_16_31 = _mm512_loadu_si512((__m512i*)(&program->pixel_offset[x + 16])); //  16 offsets

      int iStart = program->pixel_offset[x];
      __m512i m512i_Start = _mm512_set1_epi32(iStart);

      perm_0_0_15 = _mm512_sub_epi32(perm_0_0_15, m512i_Start);
      perm_0_16_31 = _mm512_sub_epi32(perm_0_16_31, m512i_Start);

      __m256i m256i_perm_0_0_15 = _mm512_cvtepi32_epi16(perm_0_0_15);
      __m256i m256i_perm_0_16_31 = _mm512_cvtepi32_epi16(perm_0_16_31);

      __m512i perm_0_0_31 = _mm512_inserti64x4(_mm512_zextsi256_si512(m256i_perm_0_0_15), m256i_perm_0_16_31, 1);

      // no add_one to each perm group to save number of registers used in processing loop
      // need only to first pair (r0 and r1)
      __m512i perm_1_0_31 = _mm512_add_epi16(perm_0_0_31, one_epi16);

      const __m512i perm_r0r1_0_31lo = _mm512_unpacklo_epi16(perm_0_0_31, perm_1_0_31);
      const __m512i perm_r0r1_0_31hi = _mm512_unpackhi_epi16(perm_0_0_31, perm_1_0_31);

      // Taps are contiguous (0, 1, 2, 3), so we increment perm indexes by 1 (in pairs by 2 each).
      const __m512i two_epi16 = _mm512_set1_epi16(2);

      uint8_t* AVS_RESTRICT dst_ptr = dst8 + x + y_from * dst_pitch;
      const uint8_t* src_ptr = src8 + iStart + y_from * src_pitch; // all permute offsets relative to this start offset

      // Calculate remaining pixels for bounds checking in partial_load mode. 1..128 remaining pixels possible.
      const int remaining = program->source_size - iStart;
	  // _bzhi_u64 creates a mask with the lower N bits set. If N >= 64, it returns all ones (~0ULL).
#ifndef _M_X64
      const __mmask64 k1 = make_low_mask64(remaining);
      const __mmask64 k2 = make_low_mask64(remaining - 64);
#else
      const __mmask64 k1 = _bzhi_u64(~0ULL, remaining);
      const __mmask64 k2 = _bzhi_u64(~0ULL, std::max(0, remaining - 64));
#endif

      // mask is used to zero out every odd byte, so that the result of the permute is a
      // vector of zero-extended 8-bit values in 16-bit lanes, preparing 8-bit data for 16-bit FMA operations in AVX512.
      const __mmask64 k_zh8 = 0x5555555555555555ULL;

      for (int y = y_from; y < y_to; y++)
      {
        __m512i data_src, data_src2;

        // working permute indexes for advancing to save number of registers used
        __m512i perm_rNrNp1_0_31lo_w = perm_r0r1_0_31lo;
        __m512i perm_rNrNp1_0_31hi_w = perm_r0r1_0_31hi;

        if constexpr (partial_load) {
          // Safe masked loads for the image edge
          data_src = _mm512_maskz_loadu_epi8(k1, src_ptr);
          data_src2 = _mm512_maskz_loadu_epi8(k2, src_ptr + 64);
        }
        else {
          // Fast unaligned loads for the safe zone
          data_src = _mm512_loadu_si512(src_ptr);
          data_src2 = _mm512_loadu_si512(src_ptr + 64);
        }

        // rows 0..3
        __m512i src_r0r1_0_31lo = _mm512_maskz_permutex2var_epi8_SIMUL<UseVNNI>(k_zh8, data_src, perm_rNrNp1_0_31lo_w, data_src2); // even index filled from permute, odd index zeroed
        __m512i src_r0r1_0_31hi = _mm512_maskz_permutex2var_epi8_SIMUL<UseVNNI>(k_zh8, data_src, perm_rNrNp1_0_31hi_w, data_src2);

        // for r2r3
        perm_rNrNp1_0_31lo_w = _mm512_add_epi16(perm_rNrNp1_0_31lo_w, two_epi16);
        perm_rNrNp1_0_31hi_w = _mm512_add_epi16(perm_rNrNp1_0_31hi_w, two_epi16);

        __m512i src_r2r3_0_31lo = _mm512_maskz_permutex2var_epi8_SIMUL<UseVNNI>(k_zh8, data_src, perm_rNrNp1_0_31lo_w, data_src2);
        __m512i src_r2r3_0_31hi = _mm512_maskz_permutex2var_epi8_SIMUL<UseVNNI>(k_zh8, data_src, perm_rNrNp1_0_31hi_w, data_src2);

        // for r4r5
        perm_rNrNp1_0_31lo_w = _mm512_add_epi16(perm_rNrNp1_0_31lo_w, two_epi16);
        perm_rNrNp1_0_31hi_w = _mm512_add_epi16(perm_rNrNp1_0_31hi_w, two_epi16);

        __m512i result_0_31lo, result_0_31hi;

        if constexpr (UseVNNI)
        {
          result_0_31lo = _mm512_dpwssd_epi32(rounder, src_r0r1_0_31lo, coef_r0r1_0_31lo);
          result_0_31lo = _mm512_dpwssd_epi32(result_0_31lo, src_r2r3_0_31lo, coef_r2r3_0_31lo);
          result_0_31hi = _mm512_dpwssd_epi32(rounder, src_r0r1_0_31hi, coef_r0r1_0_31hi);
          result_0_31hi = _mm512_dpwssd_epi32(result_0_31hi, src_r2r3_0_31hi, coef_r2r3_0_31hi);
        }
        else
        {
          // making FMA in 32bits accs as in AVX256 V-resize
          result_0_31lo = _mm512_add_epi32(_mm512_madd_epi16(src_r0r1_0_31lo, coef_r0r1_0_31lo), _mm512_madd_epi16(src_r2r3_0_31lo, coef_r2r3_0_31lo));
          result_0_31hi = _mm512_add_epi32(_mm512_madd_epi16(src_r0r1_0_31hi, coef_r0r1_0_31hi), _mm512_madd_epi16(src_r2r3_0_31hi, coef_r2r3_0_31hi));
        }

        // rows 4..7
        __m512i src_r4r5_0_31lo = _mm512_maskz_permutex2var_epi8_SIMUL<UseVNNI>(k_zh8, data_src, perm_rNrNp1_0_31lo_w, data_src2);
        __m512i src_r4r5_0_31hi = _mm512_maskz_permutex2var_epi8_SIMUL<UseVNNI>(k_zh8, data_src, perm_rNrNp1_0_31hi_w, data_src2);

        // for r6r7
        perm_rNrNp1_0_31lo_w = _mm512_add_epi16(perm_rNrNp1_0_31lo_w, two_epi16);
        perm_rNrNp1_0_31hi_w = _mm512_add_epi16(perm_rNrNp1_0_31hi_w, two_epi16);

        __m512i src_r6r7_0_31lo = _mm512_maskz_permutex2var_epi8_SIMUL<UseVNNI>(k_zh8, data_src, perm_rNrNp1_0_31lo_w, data_src2);
        __m512i src_r6r7_0_31hi = _mm512_maskz_permutex2var_epi8_SIMUL<UseVNNI>(k_zh8, data_src, perm_rNrNp1_0_31hi_w, data_src2);

        // for r8r9
        perm_rNrNp1_0_31lo_w = _mm512_add_epi16(perm_rNrNp1_0_31lo_w, two_epi16);
        perm_rNrNp1_0_31hi_w = _mm512_add_epi16(perm_rNrNp1_0_31hi_w, two_epi16);

        if constexpr (UseVNNI)
        {
          result_0_31lo = _mm512_dpwssd_epi32(result_0_31lo, src_r4r5_0_31lo, coef_r4r5_0_31lo);
          result_0_31lo = _mm512_dpwssd_epi32(result_0_31lo, src_r6r7_0_31lo, coef_r6r7_0_31lo);
          result_0_31hi = _mm512_dpwssd_epi32(result_0_31hi, src_r4r5_0_31hi, coef_r4r5_0_31hi);
          result_0_31hi = _mm512_dpwssd_epi32(result_0_31hi, src_r6r7_0_31hi, coef_r6r7_0_31hi);
        }
        else
        {
          result_0_31lo = _mm512_add_epi32(_mm512_madd_epi16(src_r4r5_0_31lo, coef_r4r5_0_31lo), result_0_31lo);
          result_0_31hi = _mm512_add_epi32(_mm512_madd_epi16(src_r4r5_0_31hi, coef_r4r5_0_31hi), result_0_31hi);
          result_0_31lo = _mm512_add_epi32(_mm512_madd_epi16(src_r6r7_0_31lo, coef_r6r7_0_31lo), result_0_31lo);
          result_0_31hi = _mm512_add_epi32(_mm512_madd_epi16(src_r6r7_0_31hi, coef_r6r7_0_31hi), result_0_31hi);
        }

        // rows 8..11
        __m512i src_r8r9_0_31lo = _mm512_maskz_permutex2var_epi8_SIMUL<UseVNNI>(k_zh8, data_src, perm_rNrNp1_0_31lo_w, data_src2);
        __m512i src_r8r9_0_31hi = _mm512_maskz_permutex2var_epi8_SIMUL<UseVNNI>(k_zh8, data_src, perm_rNrNp1_0_31hi_w, data_src2);

        // for r10r11
        perm_rNrNp1_0_31lo_w = _mm512_add_epi16(perm_rNrNp1_0_31lo_w, two_epi16);
        perm_rNrNp1_0_31hi_w = _mm512_add_epi16(perm_rNrNp1_0_31hi_w, two_epi16);

        __m512i src_r10r11_0_31lo = _mm512_maskz_permutex2var_epi8_SIMUL<UseVNNI>(k_zh8, data_src, perm_rNrNp1_0_31lo_w, data_src2);
        __m512i src_r10r11_0_31hi = _mm512_maskz_permutex2var_epi8_SIMUL<UseVNNI>(k_zh8, data_src, perm_rNrNp1_0_31hi_w, data_src2);

        // for r12r13
        perm_rNrNp1_0_31lo_w = _mm512_add_epi16(perm_rNrNp1_0_31lo_w, two_epi16);
        perm_rNrNp1_0_31hi_w = _mm512_add_epi16(perm_rNrNp1_0_31hi_w, two_epi16);

        if constexpr (UseVNNI)
        {
          result_0_31lo = _mm512_dpwssd_epi32(result_0_31lo, src_r8r9_0_31lo, coef_r8r9_0_31lo);
          result_0_31lo = _mm512_dpwssd_epi32(result_0_31lo, src_r10r11_0_31lo, coef_r10r11_0_31lo);
          result_0_31hi = _mm512_dpwssd_epi32(result_0_31hi, src_r8r9_0_31hi, coef_r8r9_0_31hi);
          result_0_31hi = _mm512_dpwssd_epi32(result_0_31hi, src_r10r11_0_31hi, coef_r10r11_0_31hi);
        }
        else
        {
          result_0_31lo = _mm512_add_epi32(_mm512_madd_epi16(src_r8r9_0_31lo, coef_r8r9_0_31lo), result_0_31lo);
          result_0_31hi = _mm512_add_epi32(_mm512_madd_epi16(src_r8r9_0_31hi, coef_r8r9_0_31hi), result_0_31hi);
          result_0_31lo = _mm512_add_epi32(_mm512_madd_epi16(src_r10r11_0_31lo, coef_r10r11_0_31lo), result_0_31lo);
          result_0_31hi = _mm512_add_epi32(_mm512_madd_epi16(src_r10r11_0_31hi, coef_r10r11_0_31hi), result_0_31hi);
        }

        // rows 12..15
        __m512i src_r12r13_0_31lo = _mm512_maskz_permutex2var_epi8_SIMUL<UseVNNI>(k_zh8, data_src, perm_rNrNp1_0_31lo_w, data_src2);
        __m512i src_r12r13_0_31hi = _mm512_maskz_permutex2var_epi8_SIMUL<UseVNNI>(k_zh8, data_src, perm_rNrNp1_0_31hi_w, data_src2);

        // for r14r15
        perm_rNrNp1_0_31lo_w = _mm512_add_epi16(perm_rNrNp1_0_31lo_w, two_epi16);
        perm_rNrNp1_0_31hi_w = _mm512_add_epi16(perm_rNrNp1_0_31hi_w, two_epi16);

        __m512i src_r14r15_0_31lo = _mm512_maskz_permutex2var_epi8_SIMUL<UseVNNI>(k_zh8, data_src, perm_rNrNp1_0_31lo_w, data_src2);
        __m512i src_r14r15_0_31hi = _mm512_maskz_permutex2var_epi8_SIMUL<UseVNNI>(k_zh8, data_src, perm_rNrNp1_0_31hi_w, data_src2);

        if constexpr (UseVNNI)
        {
          result_0_31lo = _mm512_dpwssd_epi32(result_0_31lo, src_r12r13_0_31lo, coef_r12r13_0_31lo);
          result_0_31lo = _mm512_dpwssd_epi32(result_0_31lo, src_r14r15_0_31lo, coef_r14r15_0_31lo);
          result_0_31hi = _mm512_dpwssd_epi32(result_0_31hi, src_r12r13_0_31hi, coef_r12r13_0_31hi);
          result_0_31hi = _mm512_dpwssd_epi32(result_0_31hi, src_r14r15_0_31hi, coef_r14r15_0_31hi);

          // rounding in first FMA
        }
        else
        {
          result_0_31lo = _mm512_add_epi32(_mm512_madd_epi16(src_r12r13_0_31lo, coef_r12r13_0_31lo), result_0_31lo);
          result_0_31hi = _mm512_add_epi32(_mm512_madd_epi16(src_r12r13_0_31hi, coef_r12r13_0_31hi), result_0_31hi);
          result_0_31lo = _mm512_add_epi32(_mm512_madd_epi16(src_r14r15_0_31lo, coef_r14r15_0_31lo), result_0_31lo);
          result_0_31hi = _mm512_add_epi32(_mm512_madd_epi16(src_r14r15_0_31hi, coef_r14r15_0_31hi), result_0_31hi);

          // rounding
          result_0_31lo = _mm512_add_epi32(result_0_31lo, rounder);
          result_0_31hi = _mm512_add_epi32(result_0_31hi, rounder);
        }

        // scaling down
        result_0_31lo = _mm512_srai_epi32(result_0_31lo, FPScale8bits);
        result_0_31hi = _mm512_srai_epi32(result_0_31hi, FPScale8bits);

        __m512i result_0_31_int16 = _mm512_packus_epi32(result_0_31lo, result_0_31hi);

        __m256i result_0_31_u8 = _mm512_cvtusepi16_epi8(result_0_31_int16);

		result_0_31_u8 = _mm256_min_epu8(result_0_31_u8, clamp_limit_max);
		result_0_31_u8 = _mm256_max_epu8(result_0_31_u8, clamp_limit_min);

        _mm256_stream_si256(reinterpret_cast<__m256i*>(dst_ptr), result_0_31_u8);

        dst_ptr += dst_pitch;
        src_ptr += src_pitch;
      }

      current_coeff += filter_size * PIXELS_AT_A_TIME;
      };

    // Process the 'safe zone' where direct full unaligned loads are acceptable.
    for (; x < width_safe_mod; x += PIXELS_AT_A_TIME)
    {
      do_h_integer_core(std::false_type{});
    }

    // Process the potentially 'unsafe zone' near the image edge, using safe masked loading.
    for (; x < width; x += PIXELS_AT_A_TIME)
    {
      do_h_integer_core(std::true_type{});
    }
  }
}


// uint16_t

// filter size up to 4
// 32 target uint16_t pixels at a time
// 128-byte source loads (64 uint16_t pixels)
// maximum permute index is 64 for _mm512_permutex2var_epi16 (uint16_t)
// making lo-hi unpacking with single permutex2var operation and optional VNNI FMA
template<bool lessthan16bit, bool UseVNNI>
#if defined(__clang__)
__attribute__((__target__("avx512f,avx512cd,avx512bw,avx512dq,avx512vl,avx512vnni,avx512vbmi,avx512vbmi2,avx512bitalg,avx512vpopcntdq,bmi2")))
#endif
void resize_h_planar_uint16_avx512_permutex_vstripe_mp_ks4_internal(BYTE* dst8, const BYTE* src8, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height, int bits_per_pixel,const uint8_t range,const bool mode_YUY2)
{
  AVS_UNUSED(mode_YUY2);

  const int filter_size = program->filter_size; // aligned, practically the coeff table stride

  const uint16_t* src = (uint16_t*)src8;
  uint16_t* AVS_RESTRICT dst = (uint16_t * AVS_RESTRICT)dst8;
  dst_pitch = dst_pitch / sizeof(uint16_t);
  src_pitch = src_pitch / sizeof(uint16_t);

  constexpr int PIXELS_AT_A_TIME = 32;

  // 'source_overread_beyond_targetx' indicates if the filter kernel can read beyond the target width.
  // we load 2x32 source uint16_t pixels at a time, so ensure safe overread if needed.
  // Our main loop processes calculates for 32 target pixels at a time.
  // Inside that, we load 64 source uint16_t pixels (2x32) to be able to permutex from that.
  // This we have to check at each mod-PIXELS_AT_A_TIME boundary, the allowance of 64-element source load.
  const int width_safe_mod = (program->safelimit_64_pixels_each32th_target.overread_possible ? program->safelimit_64_pixels_each32th_target.source_overread_beyond_targetx : width) / PIXELS_AT_A_TIME * PIXELS_AT_A_TIME;

  // Preconditions:

  // 'target_size_alignment' ensures we can safely access coefficients using offsets like
  // 'filter_size * 31' when processing 32 source H pixels at a time

  // Ensure that coefficient loading beyond the valid target size is safe for 4 coeffs
  // We load 4x16bit coeffs at a time

  const int max_scanlines = program->max_scanlines;

  // for 16 bits only
  const __m512i shifttosigned = _mm512_set1_epi16(-32768);
  const __m512i shiftfromsigned = _mm512_set1_epi32(32768 << FPScale16bits);

  const uint16_t val_min = (range==1) ? 0 : (int)16 << (bits_per_pixel-8);
  const uint16_t val_max = ((range==1) || (range==4)) ? ((int)1 << bits_per_pixel)-1 : (range==2) ?
    ((int)235 << (bits_per_pixel-8)) : ((int)240 << (bits_per_pixel-8));
	
  __m512i clamp_limit_min = _mm512_set1_epi16(val_min);
  __m512i clamp_limit_max = _mm512_set1_epi16(val_max);

  __m512i rounder = _mm512_set1_epi32(1 << (FPScale16bits - 1));

  // Vertical stripe loop for L2 cache optimization
  for (int y_from = 0; y_from < height; y_from += max_scanlines)
  {
    int y_to = std::min(y_from + max_scanlines, height);

    // Reset current_coeff for the start of the stripe (points to start of row's coeffs)
    const short* AVS_RESTRICT current_coeff = program->pixel_coefficient;

    int x = 0;

    // Lambda to handle both safe (fast) and unsafe (masked/partial) loading paths
    auto do_h_integer_core = [&](auto partial_load) {

      // prepare coefs in transposed V-form
      // 32 source pixels, 32 coeff strides
      // TODO: make storage in transposed form, 64 x uint16 transposition looks too slow
      // TO FIX: filter_size=4 resize uses filter_size=16 - lots or RAM/cache wasted, in the ready to use transposed form it will be much more optimized

      // 4coefs of 16bit is 64bits, can be loaded as set_epi64 ?
      __m512i coef_0_7 = _mm512_setr_epi64(
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 0),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 1),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 2),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 3),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 4),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 5),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 6),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 7)
      );

      __m512i coef_8_15 = _mm512_setr_epi64(
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 8),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 9),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 10),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 11),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 12),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 13),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 14),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 15)
      );

      __m512i coef_16_23 = _mm512_setr_epi64(
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 16),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 17),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 18),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 19),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 20),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 21),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 22),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 23)
      );

      __m512i coef_24_31 = _mm512_setr_epi64(
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 24),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 25),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 26),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 27),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 28),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 29),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 30),
        *reinterpret_cast<const long long*>(current_coeff + filter_size * 31)
      );

      // Transpose with permutex
      __m512i c_perm_0 = _mm512_set_epi16(
        28 + 32, 24 + 32, 20 + 32, 16 + 32, 12 + 32, 8 + 32, 4 + 32, 0 + 32, 28, 24, 20, 16, 12, 8, 4, 0,
        28 + 32, 24 + 32, 20 + 32, 16 + 32, 12 + 32, 8 + 32, 4 + 32, 0 + 32, 28, 24, 20, 16, 12, 8, 4, 0);
      __m512i one_epi16 = _mm512_set1_epi16(1);

      const __mmask32 k_high = 0xFFFF0000;

      // 0.0 .. 15.0 in low 256, 16.0 .. 31.0 in high 256
      __m512i coef_r0_0_31w = _mm512_mask_blend_epi16(k_high, _mm512_permutex2var_epi16(coef_0_7, c_perm_0, coef_8_15), _mm512_permutex2var_epi16(coef_16_23, c_perm_0, coef_24_31));
      c_perm_0 = _mm512_add_epi16(c_perm_0, one_epi16);
      __m512i coef_r1_0_31w = _mm512_mask_blend_epi16(k_high, _mm512_permutex2var_epi16(coef_0_7, c_perm_0, coef_8_15), _mm512_permutex2var_epi16(coef_16_23, c_perm_0, coef_24_31));
      c_perm_0 = _mm512_add_epi16(c_perm_0, one_epi16);
      __m512i coef_r2_0_31w = _mm512_mask_blend_epi16(k_high, _mm512_permutex2var_epi16(coef_0_7, c_perm_0, coef_8_15), _mm512_permutex2var_epi16(coef_16_23, c_perm_0, coef_24_31));
      c_perm_0 = _mm512_add_epi16(c_perm_0, one_epi16);
      __m512i coef_r3_0_31w = _mm512_mask_blend_epi16(k_high, _mm512_permutex2var_epi16(coef_0_7, c_perm_0, coef_8_15), _mm512_permutex2var_epi16(coef_16_23, c_perm_0, coef_24_31));

      // convert-transpose to H-pairs for madd ? better to do with single permutex in future
      // 4 to 4 512 registers - finally real working coeffs to store in the transposed resampling program for block of 64 target samples
      __m512i coef_r0r1_0_31lo = _mm512_unpacklo_epi16(coef_r0_0_31w, coef_r1_0_31w);
      __m512i coef_r0r1_0_31hi = _mm512_unpackhi_epi16(coef_r0_0_31w, coef_r1_0_31w);

      __m512i coef_r2r3_0_31lo = _mm512_unpacklo_epi16(coef_r2_0_31w, coef_r3_0_31w);
      __m512i coef_r2r3_0_31hi = _mm512_unpackhi_epi16(coef_r2_0_31w, coef_r3_0_31w);

      // TODO: store transposed resampling program coeffs to temp buffer for reusage at each line

      // convert resampling program in H-form into permuting indexes for src transposition in V-form
      __m512i perm_0_0_15 = _mm512_loadu_si512((__m512i*)(&program->pixel_offset[x])); // 16 offsets
      __m512i perm_0_16_31 = _mm512_loadu_si512((__m512i*)(&program->pixel_offset[x + 16])); //  16 offsets

      int iStart = program->pixel_offset[x];
      __m512i m512i_Start = _mm512_set1_epi32(iStart);

      perm_0_0_15 = _mm512_sub_epi32(perm_0_0_15, m512i_Start);
      perm_0_16_31 = _mm512_sub_epi32(perm_0_16_31, m512i_Start);

      __m256i m256i_perm_0_0_15 = _mm512_cvtepi32_epi16(perm_0_0_15);
      __m256i m256i_perm_0_16_31 = _mm512_cvtepi32_epi16(perm_0_16_31);

      // Insert each 256-bit register into the specific lane
      __m512i perm_0 = _mm512_inserti64x4(_mm512_zextsi256_si512(m256i_perm_0_0_15), m256i_perm_0_16_31, 1);

      // Taps are contiguous (0, 1, 2, 3), so we increment perm indexes by 1.
      __m512i perm_1 = _mm512_add_epi16(perm_0, one_epi16);
      __m512i perm_2 = _mm512_add_epi16(perm_1, one_epi16);
      __m512i perm_3 = _mm512_add_epi16(perm_2, one_epi16);

      const __m512i perm_r0r1_0_31lo = _mm512_unpacklo_epi16(perm_0, perm_1);
      const __m512i perm_r0r1_0_31hi = _mm512_unpackhi_epi16(perm_0, perm_1);

      const __m512i perm_r2r3_0_31lo = _mm512_unpacklo_epi16(perm_2, perm_3);
      const __m512i perm_r2r3_0_31hi = _mm512_unpackhi_epi16(perm_2, perm_3);

      uint16_t* AVS_RESTRICT dst_ptr = dst + x + y_from * dst_pitch;
      const uint16_t* src_ptr = src + iStart + y_from * src_pitch; // all permute offsets relative to this start offset

      // Calculate remaining pixels for bounds checking in partial_load mode. 1..64 remaining uint16_t pixels possible.
      // only when partial_load
      const int remaining = program->source_size - iStart;
      // two masks for partial loads of 32 + 32 shorts
      const __mmask32 k1 = _bzhi_u32(~0UL, remaining);
      const __mmask32 k2 = _bzhi_u32(~0UL, remaining - 32);

      for (int y = y_from; y < y_to; y++)
      {
        __m512i data_src, data_src2;

        if constexpr (partial_load) {
          // Safe masked loads for the image edge 2x32 shorts
          data_src = _mm512_maskz_loadu_epi16(k1, src_ptr);
          data_src2 = _mm512_maskz_loadu_epi16(k2, src_ptr + 32);
        }
        else {
          // Fast unaligned loads for the safe zone
          data_src = _mm512_loadu_si512(src_ptr);
          data_src2 = _mm512_loadu_si512(src_ptr + 32);
        }

        __m512i src_r0r1_0_31lo = _mm512_permutex2var_epi16(data_src, perm_r0r1_0_31lo, data_src2);
        __m512i src_r0r1_0_31hi = _mm512_permutex2var_epi16(data_src, perm_r0r1_0_31hi, data_src2);

        __m512i src_r2r3_0_31lo = _mm512_permutex2var_epi16(data_src, perm_r2r3_0_31lo, data_src2);
        __m512i src_r2r3_0_31hi = _mm512_permutex2var_epi16(data_src, perm_r2r3_0_31hi, data_src2);

        if constexpr (!lessthan16bit) {
          // madd requires signed integers, so shift to signed range
          src_r0r1_0_31lo = _mm512_add_epi16(src_r0r1_0_31lo, shifttosigned);
          src_r0r1_0_31hi = _mm512_add_epi16(src_r0r1_0_31hi, shifttosigned);
          src_r2r3_0_31lo = _mm512_add_epi16(src_r2r3_0_31lo, shifttosigned);
          src_r2r3_0_31hi = _mm512_add_epi16(src_r2r3_0_31hi, shifttosigned);
        }

        __m512i result_0_31lo, result_0_31hi;

        if constexpr (UseVNNI)
        {
          result_0_31lo = _mm512_dpwssd_epi32(rounder, src_r0r1_0_31lo, coef_r0r1_0_31lo);
          result_0_31lo = _mm512_dpwssd_epi32(result_0_31lo, src_r2r3_0_31lo, coef_r2r3_0_31lo);

          result_0_31hi = _mm512_dpwssd_epi32(rounder, src_r0r1_0_31hi, coef_r0r1_0_31hi);
          result_0_31hi = _mm512_dpwssd_epi32(result_0_31hi, src_r2r3_0_31hi, coef_r2r3_0_31hi);
        }
        else
        {
          // making FMA in 32bits accs as in AVX256 V-resize
          result_0_31lo = _mm512_add_epi32(_mm512_madd_epi16(src_r0r1_0_31lo, coef_r0r1_0_31lo), _mm512_madd_epi16(src_r2r3_0_31lo, coef_r2r3_0_31lo));
          result_0_31hi = _mm512_add_epi32(_mm512_madd_epi16(src_r0r1_0_31hi, coef_r0r1_0_31hi), _mm512_madd_epi16(src_r2r3_0_31hi, coef_r2r3_0_31hi));

          // rounding
          result_0_31lo = _mm512_add_epi32(result_0_31lo, rounder);
          result_0_31hi = _mm512_add_epi32(result_0_31hi, rounder);
        }

        if constexpr (!lessthan16bit) {
          // return from signed range
          result_0_31lo = _mm512_add_epi32(result_0_31lo, shiftfromsigned);
          result_0_31hi = _mm512_add_epi32(result_0_31hi, shiftfromsigned);
        }

        // scale down
        result_0_31lo = _mm512_srai_epi32(result_0_31lo, FPScale16bits);
        result_0_31hi = _mm512_srai_epi32(result_0_31hi, FPScale16bits);

        // negative and over 16 bit values are clamped automatically
        __m512i result_0_31_int16 = _mm512_packus_epi32(result_0_31lo, result_0_31hi);

		result_0_31_int16 = _mm512_min_epu16(result_0_31_int16, clamp_limit_max);
		result_0_31_int16 = _mm512_max_epu16(result_0_31_int16, clamp_limit_min);

        _mm512_stream_si512(reinterpret_cast<__m512i*>(dst_ptr), result_0_31_int16);

        dst_ptr += dst_pitch;
        src_ptr += src_pitch;
      }

      current_coeff += filter_size * PIXELS_AT_A_TIME;
      };

    // Process the 'safe zone' where direct full unaligned loads are acceptable.
    for (; x < width_safe_mod; x += PIXELS_AT_A_TIME)
    {
      do_h_integer_core(std::false_type{});
    }

    // Process the potentially 'unsafe zone' near the image edge, using safe masked loading.
    for (; x < width; x += PIXELS_AT_A_TIME)
    {
      do_h_integer_core(std::true_type{});
    }
  }
}

// filter size up to 8
// 32 target uint16_t pixels at a time
// 128-byte source loads (64 uint16_t pixels)
// maximum permute index is 64 for _mm512_permutex2var_epi16 (uint16_t)
// support VNNI and madd FMA
template<bool lessthan16bit, bool UseVNNI>
#if defined(__clang__)
__attribute__((__target__("avx512f,avx512cd,avx512bw,avx512dq,avx512vl,avx512vnni,avx512vbmi,avx512vbmi2,avx512bitalg,avx512vpopcntdq,bmi2")))
#endif
void resize_h_planar_uint16_avx512_permutex_vstripe_mp_ks8_internal(BYTE* dst8, const BYTE* src8, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height, int bits_per_pixel,const uint8_t range,const bool mode_YUY2)
{
  AVS_UNUSED(mode_YUY2);

  const int filter_size = program->filter_size; // aligned, practically the coeff table stride

  const uint16_t* src = (uint16_t*)src8;
  uint16_t* AVS_RESTRICT dst = (uint16_t * AVS_RESTRICT)dst8;
  dst_pitch = dst_pitch / sizeof(uint16_t);
  src_pitch = src_pitch / sizeof(uint16_t);

  constexpr int PIXELS_AT_A_TIME = 32;

  // 'source_overread_beyond_targetx' indicates if the filter kernel can read beyond the target width.
  // we load 2x32 source bytes at a time, so ensure safe overread if needed.
  // Our main loop processes calculates for 64 target pixels at a time.
  // Inside that, we load 128 source bytes (2x64) to be able to permutex from that.
  // This we have to check at each mod-PIXELS_AT_A_TIME boundary, the allowance of 128-byte source load.
  const int width_safe_mod = (program->safelimit_64_pixels_each32th_target.overread_possible ? program->safelimit_64_pixels_each32th_target.source_overread_beyond_targetx : width) / PIXELS_AT_A_TIME * PIXELS_AT_A_TIME;

  // Preconditions:

  // 'target_size_alignment' ensures we can safely access coefficients using offsets like
  // 'filter_size * 15' when processing 16 H pixels at a time
  // 'filter_size * 63' when processing 64 H pixels at a time

  // Ensure that coefficient loading beyond the valid target size is safe for 4x8 float loads.
  // We load 8x 'short' coeffs at a time
  // Loading is unaligned, but we fill __m128 registers before combining into __m512

  const int max_scanlines = program->max_scanlines;

  // for 16 bits only
  const __m512i shifttosigned = _mm512_set1_epi16(-32768);
  const __m512i shiftfromsigned = _mm512_set1_epi32(32768 << FPScale16bits);


  const uint16_t val_min = (range==1) ? 0 : (int)16 << (bits_per_pixel-8);
  const uint16_t val_max = ((range==1) || (range==4)) ? ((int)1 << bits_per_pixel)-1 : (range==2) ?
    ((int)235 << (bits_per_pixel-8)) : ((int)240 << (bits_per_pixel-8));
	
  __m512i clamp_limit_min = _mm512_set1_epi16(val_min);
  __m512i clamp_limit_max = _mm512_set1_epi16(val_max);

  __m512i rounder = _mm512_set1_epi32(1 << (FPScale16bits - 1));

  // Vertical stripe loop for L2 cache optimization
  for (int y_from = 0; y_from < height; y_from += max_scanlines)
  {
    int y_to = std::min(y_from + max_scanlines, height);

    // Reset current_coeff for the start of the stripe (points to start of row's coeffs)
    const short* AVS_RESTRICT current_coeff = program->pixel_coefficient;

    int x = 0;

    // Lambda to handle both safe (fast) and unsafe (masked/partial) loading paths
    auto do_h_integer_core = [&](auto partial_load) {

      // prepare coefs in transposed V-form
      // TODO: make storage in transposed form, 64 x uint16 transposition looks too slow

      // 8coefs of 16bit is 128bits 
      __m512i coef_0_3 = _mm512i_loadu_4_m128i(
        (__m128i*)(current_coeff + filter_size * 0), (__m128i*)(current_coeff + filter_size * 1), (__m128i*)(current_coeff + filter_size * 2), (__m128i*)(current_coeff + filter_size * 3));
      __m512i coef_4_7 = _mm512i_loadu_4_m128i(
        (__m128i*)(current_coeff + filter_size * 4), (__m128i*)(current_coeff + filter_size * 5), (__m128i*)(current_coeff + filter_size * 6), (__m128i*)(current_coeff + filter_size * 7));
      __m512i coef_8_11 = _mm512i_loadu_4_m128i(
        (__m128i*)(current_coeff + filter_size * 8), (__m128i*)(current_coeff + filter_size * 9), (__m128i*)(current_coeff + filter_size * 10), (__m128i*)(current_coeff + filter_size * 11));
      __m512i coef_12_15 = _mm512i_loadu_4_m128i(
        (__m128i*)(current_coeff + filter_size * 12), (__m128i*)(current_coeff + filter_size * 13), (__m128i*)(current_coeff + filter_size * 14), (__m128i*)(current_coeff + filter_size * 15));
      __m512i coef_16_19 = _mm512i_loadu_4_m128i(
        (__m128i*)(current_coeff + filter_size * 16), (__m128i*)(current_coeff + filter_size * 17), (__m128i*)(current_coeff + filter_size * 18), (__m128i*)(current_coeff + filter_size * 19));
      __m512i coef_20_23 = _mm512i_loadu_4_m128i(
        (__m128i*)(current_coeff + filter_size * 20), (__m128i*)(current_coeff + filter_size * 21), (__m128i*)(current_coeff + filter_size * 22), (__m128i*)(current_coeff + filter_size * 23));
      __m512i coef_24_27 = _mm512i_loadu_4_m128i(
        (__m128i*)(current_coeff + filter_size * 24), (__m128i*)(current_coeff + filter_size * 25), (__m128i*)(current_coeff + filter_size * 26), (__m128i*)(current_coeff + filter_size * 27));
      __m512i coef_28_31 = _mm512i_loadu_4_m128i(
        (__m128i*)(current_coeff + filter_size * 28), (__m128i*)(current_coeff + filter_size * 29), (__m128i*)(current_coeff + filter_size * 30), (__m128i*)(current_coeff + filter_size * 31));

      // Transpose with permutex
      __m512i c_perm_0_7 = _mm512_set_epi16(
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        24 + 32, 16 + 32, 8 + 32, 0 + 32, 24, 16, 8, 0);
      __m512i c_perm_8_15 = _mm512_set_epi16(
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        24 + 32, 16 + 32, 8 + 32, 0 + 32, 24, 16, 8, 0,
        0, 0, 0, 0, 0, 0, 0, 0);
      __m512i c_perm_16_23 = _mm512_set_epi16(
        0, 0, 0, 0, 0, 0, 0, 0,
        24 + 32, 16 + 32, 8 + 32, 0 + 32, 24, 16, 8, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0);
      __m512i c_perm_24_31 = _mm512_set_epi16(
        24 + 32, 16 + 32, 8 + 32, 0 + 32, 24, 16, 8, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0);

      __m512i one_epi16 = _mm512_set1_epi16(1);
      const __mmask32 k_8_15 = 0x0000FF00;
      const __mmask32 k_16_23 = 0x00FF0000;
      const __mmask32 k_24_31 = 0xFF000000;

      auto inc_perms = [&](
        __m512i& c0_7,
        __m512i& c8_15,
        __m512i& c16_23,
        __m512i& c24_31
        ) {
          c0_7 = _mm512_add_epi16(c0_7, one_epi16);
          c8_15 = _mm512_add_epi16(c8_15, one_epi16);
          c16_23 = _mm512_add_epi16(c16_23, one_epi16);
          c24_31 = _mm512_add_epi16(c24_31, one_epi16);
        };

      auto make_row_0_31 = [&](
        __m512i& row_0_31w,
        __m512i c0_7, __m512i c8_15, __m512i c16_23, __m512i c24_31
        ) {
          // 0..31
          row_0_31w = _mm512_mask_blend_epi16(
            k_8_15,
            _mm512_permutex2var_epi16(coef_0_3, c0_7, coef_4_7),
            _mm512_permutex2var_epi16(coef_8_11, c8_15, coef_12_15)
          );
          row_0_31w = _mm512_mask_blend_epi16(
            k_16_23,
            row_0_31w,
            _mm512_permutex2var_epi16(coef_16_19, c16_23, coef_20_23)
          );
          row_0_31w = _mm512_mask_blend_epi16(
            k_24_31,
            row_0_31w,
            _mm512_permutex2var_epi16(coef_24_27, c24_31, coef_28_31)
          );
        };

      __m512i coef_r0_0_31w;
      __m512i coef_r1_0_31w;
      __m512i coef_r2_0_31w;
      __m512i coef_r3_0_31w;
      __m512i coef_r4_0_31w;
      __m512i coef_r5_0_31w;
      __m512i coef_r6_0_31w;
      __m512i coef_r7_0_31w;

      // r0
      make_row_0_31(coef_r0_0_31w, c_perm_0_7, c_perm_8_15, c_perm_16_23, c_perm_24_31);
      inc_perms(c_perm_0_7, c_perm_8_15, c_perm_16_23, c_perm_24_31);
      // r1
      make_row_0_31(coef_r1_0_31w, c_perm_0_7, c_perm_8_15, c_perm_16_23, c_perm_24_31);
      inc_perms(c_perm_0_7, c_perm_8_15, c_perm_16_23, c_perm_24_31);
      // r2
      make_row_0_31(coef_r2_0_31w, c_perm_0_7, c_perm_8_15, c_perm_16_23, c_perm_24_31);
      inc_perms(c_perm_0_7, c_perm_8_15, c_perm_16_23, c_perm_24_31);
      // r3
      make_row_0_31(coef_r3_0_31w, c_perm_0_7, c_perm_8_15, c_perm_16_23, c_perm_24_31);
      inc_perms(c_perm_0_7, c_perm_8_15, c_perm_16_23, c_perm_24_31);
      // r4
      make_row_0_31(coef_r4_0_31w, c_perm_0_7, c_perm_8_15, c_perm_16_23, c_perm_24_31);
      inc_perms(c_perm_0_7, c_perm_8_15, c_perm_16_23, c_perm_24_31);
      // r5
      make_row_0_31(coef_r5_0_31w, c_perm_0_7, c_perm_8_15, c_perm_16_23, c_perm_24_31);
      inc_perms(c_perm_0_7, c_perm_8_15, c_perm_16_23, c_perm_24_31);
      // r6
      make_row_0_31(coef_r6_0_31w, c_perm_0_7, c_perm_8_15, c_perm_16_23, c_perm_24_31);
      inc_perms(c_perm_0_7, c_perm_8_15, c_perm_16_23, c_perm_24_31);
      // r7
      make_row_0_31(coef_r7_0_31w, c_perm_0_7, c_perm_8_15, c_perm_16_23, c_perm_24_31);

      // convert-transpose to H-pairs for madd ? better to do with single permutex in future
      // 8 to 8 512 registers - finally real working coeffs to store in the transposed resampling program for block of 64 target samples
      __m512i coef_r0r1_0_31lo = _mm512_unpacklo_epi16(coef_r0_0_31w, coef_r1_0_31w);
      __m512i coef_r0r1_0_31hi = _mm512_unpackhi_epi16(coef_r0_0_31w, coef_r1_0_31w);

      __m512i coef_r2r3_0_31lo = _mm512_unpacklo_epi16(coef_r2_0_31w, coef_r3_0_31w);
      __m512i coef_r2r3_0_31hi = _mm512_unpackhi_epi16(coef_r2_0_31w, coef_r3_0_31w);

      __m512i coef_r4r5_0_31lo = _mm512_unpacklo_epi16(coef_r4_0_31w, coef_r5_0_31w);
      __m512i coef_r4r5_0_31hi = _mm512_unpackhi_epi16(coef_r4_0_31w, coef_r5_0_31w);

      __m512i coef_r6r7_0_31lo = _mm512_unpacklo_epi16(coef_r6_0_31w, coef_r7_0_31w);
      __m512i coef_r6r7_0_31hi = _mm512_unpackhi_epi16(coef_r6_0_31w, coef_r7_0_31w);

      // TODO: store transposed resampling program coeffs to temp buffer for reusage at each line

      // convert resampling program in H-form into permuting indexes for src transposition in V-form
      __m512i perm_0_0_15 = _mm512_loadu_si512((__m512i*)(&program->pixel_offset[x])); // 16 offsets
      __m512i perm_0_16_31 = _mm512_loadu_si512((__m512i*)(&program->pixel_offset[x + 16])); //  16 offsets

      int iStart = program->pixel_offset[x];
      __m512i m512i_Start = _mm512_set1_epi32(iStart);

      perm_0_0_15 = _mm512_sub_epi32(perm_0_0_15, m512i_Start);
      perm_0_16_31 = _mm512_sub_epi32(perm_0_16_31, m512i_Start);

      __m256i m256i_perm_0_0_15 = _mm512_cvtepi32_epi16(perm_0_0_15);
      __m256i m256i_perm_0_16_31 = _mm512_cvtepi32_epi16(perm_0_16_31);

      // Insert each 256-bit register into the specific lane
      __m512i perm_0_0_31 = _mm512_inserti64x4(_mm512_zextsi256_si512(m256i_perm_0_0_15), m256i_perm_0_16_31, 1);

      // no add_one to each perm group to save number of registers used in processing loop
      // need only to first pair (r0 and r1)
      __m512i perm_1_0_31 = _mm512_add_epi16(perm_0_0_31, one_epi16);

      const __m512i perm_r0r1_0_31lo = _mm512_unpacklo_epi16(perm_0_0_31, perm_1_0_31);
      const __m512i perm_r0r1_0_31hi = _mm512_unpackhi_epi16(perm_0_0_31, perm_1_0_31);

      // Taps are contiguous (0, 1, 2, 3), so we increment perm indexes by 1 (in pairs by 2 each).
      const __m512i two_epi16 = _mm512_set1_epi16(2);

      uint16_t* AVS_RESTRICT dst_ptr = dst + x + y_from * dst_pitch;
      const uint16_t* src_ptr = src + iStart + y_from * src_pitch; // all permute offsets relative to this start offset

      // Calculate remaining pixels for bounds checking in partial_load mode. 1..64 remaining uint16_t pixels possible.
      // only when partial_load
      const int remaining = program->source_size - iStart;
      // two masks for partial loads of 32 + 32 shorts
      const __mmask32 k1 = _bzhi_u32(~0UL, remaining);
      const __mmask32 k2 = _bzhi_u32(~0UL, remaining - 32);

      for (int y = y_from; y < y_to; y++)
      {
        __m512i data_src, data_src2;

        // working permute indexes for advancing to save number of registers used
        __m512i perm_rNrNp1_0_31lo_w = perm_r0r1_0_31lo;
        __m512i perm_rNrNp1_0_31hi_w = perm_r0r1_0_31hi;

        if constexpr (partial_load) {
          // Safe masked loads for the image edge
          data_src = _mm512_maskz_loadu_epi16(k1, src_ptr);
          data_src2 = _mm512_maskz_loadu_epi16(k2, src_ptr + 32);
        }
        else {
          // Fast unaligned loads for the safe zone
          data_src = _mm512_loadu_si512(src_ptr);
          data_src2 = _mm512_loadu_si512(src_ptr + 32);
        }

        // rows 0..3
        __m512i src_r0r1_0_31lo = _mm512_permutex2var_epi16(data_src, perm_rNrNp1_0_31lo_w, data_src2);
        __m512i src_r0r1_0_31hi = _mm512_permutex2var_epi16(data_src, perm_rNrNp1_0_31hi_w, data_src2);

        // for r2r3
        perm_rNrNp1_0_31lo_w = _mm512_add_epi16(perm_rNrNp1_0_31lo_w, two_epi16);
        perm_rNrNp1_0_31hi_w = _mm512_add_epi16(perm_rNrNp1_0_31hi_w, two_epi16);

        __m512i src_r2r3_0_31lo = _mm512_permutex2var_epi16(data_src, perm_rNrNp1_0_31lo_w, data_src2);
        __m512i src_r2r3_0_31hi = _mm512_permutex2var_epi16(data_src, perm_rNrNp1_0_31hi_w, data_src2);

        // for r4r5
        perm_rNrNp1_0_31lo_w = _mm512_add_epi16(perm_rNrNp1_0_31lo_w, two_epi16);
        perm_rNrNp1_0_31hi_w = _mm512_add_epi16(perm_rNrNp1_0_31hi_w, two_epi16);

        __m512i result_0_31lo, result_0_31hi;

        if constexpr (!lessthan16bit) {
          // madd requires signed integers, so shift to signed range
          src_r0r1_0_31lo = _mm512_add_epi16(src_r0r1_0_31lo, shifttosigned);
          src_r0r1_0_31hi = _mm512_add_epi16(src_r0r1_0_31hi, shifttosigned);
          src_r2r3_0_31lo = _mm512_add_epi16(src_r2r3_0_31lo, shifttosigned);
          src_r2r3_0_31hi = _mm512_add_epi16(src_r2r3_0_31hi, shifttosigned);
        }

        if constexpr (UseVNNI)
        {
          result_0_31lo = _mm512_dpwssd_epi32(rounder, src_r0r1_0_31lo, coef_r0r1_0_31lo);
          result_0_31lo = _mm512_dpwssd_epi32(result_0_31lo, src_r2r3_0_31lo, coef_r2r3_0_31lo);

          result_0_31hi = _mm512_dpwssd_epi32(rounder, src_r0r1_0_31hi, coef_r0r1_0_31hi);
          result_0_31hi = _mm512_dpwssd_epi32(result_0_31hi, src_r2r3_0_31hi, coef_r2r3_0_31hi);
        }
        else
        {
          // making FMA in 32bits accs as in AVX256 V-resize
          result_0_31lo = _mm512_add_epi32(_mm512_madd_epi16(src_r0r1_0_31lo, coef_r0r1_0_31lo), _mm512_madd_epi16(src_r2r3_0_31lo, coef_r2r3_0_31lo));
          result_0_31hi = _mm512_add_epi32(_mm512_madd_epi16(src_r0r1_0_31hi, coef_r0r1_0_31hi), _mm512_madd_epi16(src_r2r3_0_31hi, coef_r2r3_0_31hi));
        }

        // rows 4..7
        __m512i src_r4r5_0_31lo = _mm512_permutex2var_epi16(data_src, perm_rNrNp1_0_31lo_w, data_src2);
        __m512i src_r4r5_0_31hi = _mm512_permutex2var_epi16(data_src, perm_rNrNp1_0_31hi_w, data_src2);

        // for r6r7
        perm_rNrNp1_0_31lo_w = _mm512_add_epi16(perm_rNrNp1_0_31lo_w, two_epi16);
        perm_rNrNp1_0_31hi_w = _mm512_add_epi16(perm_rNrNp1_0_31hi_w, two_epi16);

        __m512i src_r6r7_0_31lo = _mm512_permutex2var_epi16(data_src, perm_rNrNp1_0_31lo_w, data_src2);
        __m512i src_r6r7_0_31hi = _mm512_permutex2var_epi16(data_src, perm_rNrNp1_0_31hi_w, data_src2);

        if constexpr (!lessthan16bit) {
          // madd requires signed integers, so shift to signed range
          src_r4r5_0_31lo = _mm512_add_epi16(src_r4r5_0_31lo, shifttosigned);
          src_r4r5_0_31hi = _mm512_add_epi16(src_r4r5_0_31hi, shifttosigned);
          src_r6r7_0_31lo = _mm512_add_epi16(src_r6r7_0_31lo, shifttosigned);
          src_r6r7_0_31hi = _mm512_add_epi16(src_r6r7_0_31hi, shifttosigned);
        }

        if constexpr (UseVNNI)
        {
          result_0_31lo = _mm512_dpwssd_epi32(result_0_31lo, src_r4r5_0_31lo, coef_r4r5_0_31lo);
          result_0_31lo = _mm512_dpwssd_epi32(result_0_31lo, src_r6r7_0_31lo, coef_r6r7_0_31lo);

          result_0_31hi = _mm512_dpwssd_epi32(result_0_31hi, src_r4r5_0_31hi, coef_r4r5_0_31hi);
          result_0_31hi = _mm512_dpwssd_epi32(result_0_31hi, src_r6r7_0_31hi, coef_r6r7_0_31hi);

          // rounding VNNI in first FMA already summed
        }
        else
        {
          result_0_31lo = _mm512_add_epi32(_mm512_madd_epi16(src_r4r5_0_31lo, coef_r4r5_0_31lo), result_0_31lo);
          result_0_31hi = _mm512_add_epi32(_mm512_madd_epi16(src_r4r5_0_31hi, coef_r4r5_0_31hi), result_0_31hi);

          result_0_31lo = _mm512_add_epi32(_mm512_madd_epi16(src_r6r7_0_31lo, coef_r6r7_0_31lo), result_0_31lo);
          result_0_31hi = _mm512_add_epi32(_mm512_madd_epi16(src_r6r7_0_31hi, coef_r6r7_0_31hi), result_0_31hi);

          // rounding
          result_0_31lo = _mm512_add_epi32(result_0_31lo, rounder);
          result_0_31hi = _mm512_add_epi32(result_0_31hi, rounder);
        }

        if constexpr (!lessthan16bit) {
          // return from signed range
          result_0_31lo = _mm512_add_epi32(result_0_31lo, shiftfromsigned);
          result_0_31hi = _mm512_add_epi32(result_0_31hi, shiftfromsigned);
        }

        // scaling down
        result_0_31lo = _mm512_srai_epi32(result_0_31lo, FPScale16bits);
        result_0_31hi = _mm512_srai_epi32(result_0_31hi, FPScale16bits);

        __m512i result_0_31_int16 = _mm512_packus_epi32(result_0_31lo, result_0_31hi);

		result_0_31_int16 = _mm512_min_epu16(result_0_31_int16, clamp_limit_max);
		result_0_31_int16 = _mm512_max_epu16(result_0_31_int16, clamp_limit_min);

        _mm512_stream_si512(reinterpret_cast<__m512i*>(dst_ptr), result_0_31_int16);

        dst_ptr += dst_pitch;
        src_ptr += src_pitch;
      }

      current_coeff += filter_size * PIXELS_AT_A_TIME;
      };

    // Process the 'safe zone' where direct full unaligned loads are acceptable.
    for (; x < width_safe_mod; x += PIXELS_AT_A_TIME)
    {
      do_h_integer_core(std::false_type{});
    }

    // Process the potentially 'unsafe zone' near the image edge, using safe masked loading.
    for (; x < width; x += PIXELS_AT_A_TIME)
    {
      do_h_integer_core(std::true_type{});
    }
  }
}

// filter size up to 16
// 32 target uint16_t pixels at a time
// 128-byte source loads (64 uint16_t pixels)
// maximum permute index is 64 for _mm512_permutex2var_epi16 (uint16_t)
// support VNNI and madd FMA
template<bool lessthan16bit, bool UseVNNI>
#if defined(__clang__)
__attribute__((__target__("avx512f,avx512cd,avx512bw,avx512dq,avx512vl,avx512vnni,avx512vbmi,avx512vbmi2,avx512bitalg,avx512vpopcntdq,bmi2")))
#endif
void resize_h_planar_uint16_avx512_permutex_vstripe_mp_ks16_internal(BYTE* dst8, const BYTE* src8, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height, int bits_per_pixel,const uint8_t range,const bool mode_YUY2)
{
  AVS_UNUSED(mode_YUY2);

  const int filter_size = program->filter_size; // aligned, practically the coeff table stride

  const uint16_t* src = (uint16_t*)src8;
  uint16_t* AVS_RESTRICT dst = (uint16_t * AVS_RESTRICT)dst8;
  dst_pitch = dst_pitch / sizeof(uint16_t);
  src_pitch = src_pitch / sizeof(uint16_t);

  constexpr int PIXELS_AT_A_TIME = 32;

  // 'source_overread_beyond_targetx' indicates if the filter kernel can read beyond the target width.
  // we load 2x32 source bytes at a time, so ensure safe overread if needed.
  // Our main loop processes calculates for 64 target pixels at a time.
  // Inside that, we load 128 source bytes (2x64) to be able to permutex from that.
  // This we have to check at each mod-PIXELS_AT_A_TIME boundary, the allowance of 128-byte source load.
  const int width_safe_mod = (program->safelimit_64_pixels_each32th_target.overread_possible ? program->safelimit_64_pixels_each32th_target.source_overread_beyond_targetx : width) / PIXELS_AT_A_TIME * PIXELS_AT_A_TIME;

  // Preconditions:

  // 'target_size_alignment' ensures we can safely access coefficients using offsets like
  // 'filter_size * 15' when processing 16 H pixels at a time
  // 'filter_size * 31' when processing 64 H pixels at a time


  // Ensure that coefficient loading beyond the valid target size is safe for 4x8 float loads.
  // We load 8x 'short' coeffs at a time
  // Loading is unaligned, but we fill __m128 registers before combining into __m512

  const int max_scanlines = program->max_scanlines;

  // for 16 bits only
  const __m512i shifttosigned = _mm512_set1_epi16(-32768);
  const __m512i shiftfromsigned = _mm512_set1_epi32(32768 << FPScale16bits);

  const uint16_t val_min = (range==1) ? 0 : (int)16 << (bits_per_pixel-8);
  const uint16_t val_max = ((range==1) || (range==4)) ? ((int)1 << bits_per_pixel)-1 : (range==2) ?
    ((int)235 << (bits_per_pixel-8)) : ((int)240 << (bits_per_pixel-8));
	
  __m512i clamp_limit_min = _mm512_set1_epi16(val_min);
  __m512i clamp_limit_max = _mm512_set1_epi16(val_max);

  __m512i rounder = _mm512_set1_epi32(1 << (FPScale16bits - 1));

  // Vertical stripe loop for L2 cache optimization
  for (int y_from = 0; y_from < height; y_from += max_scanlines)
  {
    int y_to = std::min(y_from + max_scanlines, height);

    // Reset current_coeff for the start of the stripe (points to start of row's coeffs)
    const short* AVS_RESTRICT current_coeff = program->pixel_coefficient;

    int x = 0;

    // Lambda to handle both safe (fast) and unsafe (masked/partial) loading paths
    auto do_h_integer_core = [&](auto partial_load) {

      // prepare coefs in transposed V-form
      // TODO: make storage in transposed form, 64 x uint16 transposition looks too slow

      // 16coefs of 16bit is 256bits - load as pairs of _m256i
      // hope 16-coeffs blocks are 32-bytes aligned for m256i aligned loads ?
      __m512i coef_0_1 = _mm512i_load_2_m256i((__m256i*)(current_coeff + filter_size * 0), (__m256i*)(current_coeff + filter_size * 1));
      __m512i coef_2_3 = _mm512i_load_2_m256i((__m256i*)(current_coeff + filter_size * 2), (__m256i*)(current_coeff + filter_size * 3));
      __m512i coef_4_5 = _mm512i_load_2_m256i((__m256i*)(current_coeff + filter_size * 4), (__m256i*)(current_coeff + filter_size * 5));
      __m512i coef_6_7 = _mm512i_load_2_m256i((__m256i*)(current_coeff + filter_size * 6), (__m256i*)(current_coeff + filter_size * 7));
      __m512i coef_8_9 = _mm512i_load_2_m256i((__m256i*)(current_coeff + filter_size * 8), (__m256i*)(current_coeff + filter_size * 9));
      __m512i coef_10_11 = _mm512i_load_2_m256i((__m256i*)(current_coeff + filter_size * 10), (__m256i*)(current_coeff + filter_size * 11));
      __m512i coef_12_13 = _mm512i_load_2_m256i((__m256i*)(current_coeff + filter_size * 12), (__m256i*)(current_coeff + filter_size * 13));
      __m512i coef_14_15 = _mm512i_load_2_m256i((__m256i*)(current_coeff + filter_size * 14), (__m256i*)(current_coeff + filter_size * 15));
      __m512i coef_16_17 = _mm512i_load_2_m256i((__m256i*)(current_coeff + filter_size * 16), (__m256i*)(current_coeff + filter_size * 17));
      __m512i coef_18_19 = _mm512i_load_2_m256i((__m256i*)(current_coeff + filter_size * 18), (__m256i*)(current_coeff + filter_size * 19));
      __m512i coef_20_21 = _mm512i_load_2_m256i((__m256i*)(current_coeff + filter_size * 20), (__m256i*)(current_coeff + filter_size * 21));
      __m512i coef_22_23 = _mm512i_load_2_m256i((__m256i*)(current_coeff + filter_size * 22), (__m256i*)(current_coeff + filter_size * 23));
      __m512i coef_24_25 = _mm512i_load_2_m256i((__m256i*)(current_coeff + filter_size * 24), (__m256i*)(current_coeff + filter_size * 25));
      __m512i coef_26_27 = _mm512i_load_2_m256i((__m256i*)(current_coeff + filter_size * 26), (__m256i*)(current_coeff + filter_size * 27));
      __m512i coef_28_29 = _mm512i_load_2_m256i((__m256i*)(current_coeff + filter_size * 28), (__m256i*)(current_coeff + filter_size * 29));
      __m512i coef_30_31 = _mm512i_load_2_m256i((__m256i*)(current_coeff + filter_size * 30), (__m256i*)(current_coeff + filter_size * 31));

      // Transpose with permutex
      __m512i c_perm_0_3 = _mm512_set_epi16(
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        16 + 32, 0 + 32, 16, 0);

      __m512i c_perm_4_7 = _mm512_set_epi16(
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        16 + 32, 0 + 32, 16, 0,
        0, 0, 0, 0);

      __m512i c_perm_8_11 = _mm512_set_epi16(
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        16 + 32, 0 + 32, 16, 0,
        0, 0, 0, 0,
        0, 0, 0, 0);

      __m512i c_perm_12_15 = _mm512_set_epi16(
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        16 + 32, 0 + 32, 16, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0);

      __m512i c_perm_16_19 = _mm512_set_epi16(
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        16 + 32, 0 + 32, 16, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0);

      __m512i c_perm_20_23 = _mm512_set_epi16(
        0, 0, 0, 0,
        0, 0, 0, 0,
        16 + 32, 0 + 32, 16, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0);

      __m512i c_perm_24_27 = _mm512_set_epi16(
        0, 0, 0, 0,
        16 + 32, 0 + 32, 16, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0);

      __m512i c_perm_28_31 = _mm512_set_epi16(
        16 + 32, 0 + 32, 16, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0);


      __m512i one_epi16 = _mm512_set1_epi16(1);

      // Define masks for the 4-word (8-byte, 4x 16-bit word) segments within the 32-word vector.
      // Note: Each AVX-512 word (epi16) lane corresponds to one bit in the __mmask32.
      // The masks target chunks of 4 words (k_X_Y represents bits X through Y inclusive).
      const __mmask32 k_4_7 = 0x000000F0;
      const __mmask32 k_8_11 = 0x00000F00;
      const __mmask32 k_12_15 = 0x0000F000;
      const __mmask32 k_16_19 = 0x000F0000;
      const __mmask32 k_20_23 = 0x00F00000;
      const __mmask32 k_24_27 = 0x0F000000;
      const __mmask32 k_28_31 = 0xF0000000;

      // Helper lambda to increment all eight permutation vectors by one.
      auto inc_perms = [&](
        __m512i& p0_3, __m512i& p4_7, __m512i& p8_11, __m512i& p12_15,
        __m512i& p16_19, __m512i& p20_23, __m512i& p24_27, __m512i& p28_31
        ) {
          p0_3 = _mm512_add_epi16(p0_3, one_epi16);
          p4_7 = _mm512_add_epi16(p4_7, one_epi16);
          p8_11 = _mm512_add_epi16(p8_11, one_epi16);
          p12_15 = _mm512_add_epi16(p12_15, one_epi16);
          p16_19 = _mm512_add_epi16(p16_19, one_epi16);
          p20_23 = _mm512_add_epi16(p20_23, one_epi16);
          p24_27 = _mm512_add_epi16(p24_27, one_epi16);
          p28_31 = _mm512_add_epi16(p28_31, one_epi16);
        };

      // Helper lambda to construct one full 32-word (512-bit) coefficient row.
      // It uses mask blending to merge the results of _mm512_permutex2var_epi16
      // for different 4-word segments, using the current permutation vectors.
      auto make_coef_row = [&](
        __m512i& row_result,
        __m512i p0_3, __m512i p4_7, __m512i p8_11, __m512i p12_15,
        __m512i p16_19, __m512i p20_23, __m512i p24_27, __m512i p28_31
        ) {
          // Start with the first segment (words 0-3) and the fourth segment (words 4-7).
          row_result = _mm512_mask_blend_epi16(
            k_4_7,
            _mm512_permutex2var_epi16(coef_0_1, p0_3, coef_2_3),  // words 0-3 (unmasked)
            _mm512_permutex2var_epi16(coef_4_5, p4_7, coef_6_7)   // words 4-7 (masked)
          );

          // Merge segment 8-11
          row_result = _mm512_mask_blend_epi16(
            k_8_11,
            row_result,
            _mm512_permutex2var_epi16(coef_8_9, p8_11, coef_10_11)
          );

          // Merge segment 12-15
          row_result = _mm512_mask_blend_epi16(
            k_12_15,
            row_result,
            _mm512_permutex2var_epi16(coef_12_13, p12_15, coef_14_15)
          );

          // Merge segment 16-19
          row_result = _mm512_mask_blend_epi16(
            k_16_19,
            row_result,
            _mm512_permutex2var_epi16(coef_16_17, p16_19, coef_18_19)
          );

          // Merge segment 20-23
          row_result = _mm512_mask_blend_epi16(
            k_20_23,
            row_result,
            _mm512_permutex2var_epi16(coef_20_21, p20_23, coef_22_23)
          );

          // Merge segment 24-27
          row_result = _mm512_mask_blend_epi16(
            k_24_27,
            row_result,
            _mm512_permutex2var_epi16(coef_24_25, p24_27, coef_26_27)
          );

          // Merge segment 28-31
          row_result = _mm512_mask_blend_epi16(
            k_28_31,
            row_result,
            _mm512_permutex2var_epi16(coef_28_29, p28_31, coef_30_31)
          );
        };

      // Declare the coefficient row variables
      __m512i coef_r0_0_31w, coef_r1_0_31w, coef_r2_0_31w, coef_r3_0_31w;
      __m512i coef_r4_0_31w, coef_r5_0_31w, coef_r6_0_31w, coef_r7_0_31w;
      __m512i coef_r8_0_31w, coef_r9_0_31w, coef_r10_0_31w, coef_r11_0_31w;
      __m512i coef_r12_0_31w, coef_r13_0_31w, coef_r14_0_31w, coef_r15_0_31w;

      // Process Row 0
      make_coef_row(coef_r0_0_31w, c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);
      inc_perms(c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);

      // Process Row 1
      make_coef_row(coef_r1_0_31w, c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);
      inc_perms(c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);

      // Process Row 2
      make_coef_row(coef_r2_0_31w, c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);
      inc_perms(c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);

      // Process Row 3
      make_coef_row(coef_r3_0_31w, c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);
      inc_perms(c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);

      // Process Row 4
      make_coef_row(coef_r4_0_31w, c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);
      inc_perms(c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);

      // Process Row 5
      make_coef_row(coef_r5_0_31w, c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);
      inc_perms(c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);

      // Process Row 6
      make_coef_row(coef_r6_0_31w, c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);
      inc_perms(c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);

      // Process Row 7
      make_coef_row(coef_r7_0_31w, c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);
      inc_perms(c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);

      // Process Row 8
      make_coef_row(coef_r8_0_31w, c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);
      inc_perms(c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);

      // Process Row 9
      make_coef_row(coef_r9_0_31w, c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);
      inc_perms(c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);

      // Process Row 10
      make_coef_row(coef_r10_0_31w, c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);
      inc_perms(c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);

      // Process Row 11
      make_coef_row(coef_r11_0_31w, c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);
      inc_perms(c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);

      // Process Row 12
      make_coef_row(coef_r12_0_31w, c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);
      inc_perms(c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);

      // Process Row 13
      make_coef_row(coef_r13_0_31w, c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);
      inc_perms(c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);

      // Process Row 14
      make_coef_row(coef_r14_0_31w, c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);
      inc_perms(c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);

      // Process Row 15
      make_coef_row(coef_r15_0_31w, c_perm_0_3, c_perm_4_7, c_perm_8_11, c_perm_12_15, c_perm_16_19, c_perm_20_23, c_perm_24_27, c_perm_28_31);
      // No inc_perms here

      // convert-transpose to H-pairs for madd ? better to do with single permutex in future
      // 16 to 16 512 registers - finally real working coeffs to store in the transposed resampling program for block of 64 target samples
      __m512i coef_r0r1_0_31lo = _mm512_unpacklo_epi16(coef_r0_0_31w, coef_r1_0_31w);
      __m512i coef_r0r1_0_31hi = _mm512_unpackhi_epi16(coef_r0_0_31w, coef_r1_0_31w);

      __m512i coef_r2r3_0_31lo = _mm512_unpacklo_epi16(coef_r2_0_31w, coef_r3_0_31w);
      __m512i coef_r2r3_0_31hi = _mm512_unpackhi_epi16(coef_r2_0_31w, coef_r3_0_31w);

      __m512i coef_r4r5_0_31lo = _mm512_unpacklo_epi16(coef_r4_0_31w, coef_r5_0_31w);
      __m512i coef_r4r5_0_31hi = _mm512_unpackhi_epi16(coef_r4_0_31w, coef_r5_0_31w);

      __m512i coef_r6r7_0_31lo = _mm512_unpacklo_epi16(coef_r6_0_31w, coef_r7_0_31w);
      __m512i coef_r6r7_0_31hi = _mm512_unpackhi_epi16(coef_r6_0_31w, coef_r7_0_31w);

      __m512i coef_r8r9_0_31lo = _mm512_unpacklo_epi16(coef_r8_0_31w, coef_r9_0_31w);
      __m512i coef_r8r9_0_31hi = _mm512_unpackhi_epi16(coef_r8_0_31w, coef_r9_0_31w);

      __m512i coef_r10r11_0_31lo = _mm512_unpacklo_epi16(coef_r10_0_31w, coef_r11_0_31w);
      __m512i coef_r10r11_0_31hi = _mm512_unpackhi_epi16(coef_r10_0_31w, coef_r11_0_31w);

      __m512i coef_r12r13_0_31lo = _mm512_unpacklo_epi16(coef_r12_0_31w, coef_r13_0_31w);
      __m512i coef_r12r13_0_31hi = _mm512_unpackhi_epi16(coef_r12_0_31w, coef_r13_0_31w);

      __m512i coef_r14r15_0_31lo = _mm512_unpacklo_epi16(coef_r14_0_31w, coef_r15_0_31w);
      __m512i coef_r14r15_0_31hi = _mm512_unpackhi_epi16(coef_r14_0_31w, coef_r15_0_31w);

      // TODO: store transposed resampling program coeffs to temp buffer for reusage at each line

      // convert resampling program in H-form into permuting indexes for src transposition in V-form
      __m512i perm_0_0_15 = _mm512_loadu_si512((__m512i*)(&program->pixel_offset[x])); // 16 offsets
      __m512i perm_0_16_31 = _mm512_loadu_si512((__m512i*)(&program->pixel_offset[x + 16])); //  16 offsets

      int iStart = program->pixel_offset[x];
      __m512i m512i_Start = _mm512_set1_epi32(iStart);

      perm_0_0_15 = _mm512_sub_epi32(perm_0_0_15, m512i_Start);
      perm_0_16_31 = _mm512_sub_epi32(perm_0_16_31, m512i_Start);

      __m256i m256i_perm_0_0_15 = _mm512_cvtepi32_epi16(perm_0_0_15);
      __m256i m256i_perm_0_16_31 = _mm512_cvtepi32_epi16(perm_0_16_31);

      // Insert each 256-bit register into the specific lane
      __m512i perm_0_0_31 = _mm512_inserti64x4(_mm512_zextsi256_si512(m256i_perm_0_0_15), m256i_perm_0_16_31, 1);

      // no add_one to each perm group to save number of registers used in processing loop
      // need only to first pair (r0 and r1)
      __m512i perm_1_0_31 = _mm512_add_epi16(perm_0_0_31, one_epi16);

      const __m512i perm_r0r1_0_31lo = _mm512_unpacklo_epi16(perm_0_0_31, perm_1_0_31);
      const __m512i perm_r0r1_0_31hi = _mm512_unpackhi_epi16(perm_0_0_31, perm_1_0_31);

      // Taps are contiguous (0, 1, 2, 3), so we increment perm indexes by 1 (in pairs by 2 each).
      const __m512i two_epi16 = _mm512_set1_epi16(2);

      uint16_t* AVS_RESTRICT dst_ptr = dst + x + y_from * dst_pitch;
      const uint16_t* src_ptr = src + iStart + y_from * src_pitch; // all permute offsets relative to this start offset

      // Calculate remaining pixels for bounds checking in partial_load mode. 1..64 remaining uint16_t pixels possible.
      // only when partial_load
      const int remaining = program->source_size - iStart;
      // two masks for partial loads of 32 + 32 shorts
      const __mmask32 k1 = _bzhi_u32(~0UL, remaining);
      const __mmask32 k2 = _bzhi_u32(~0UL, remaining - 32);

      for (int y = y_from; y < y_to; y++)
      {
        __m512i data_src, data_src2;

        // working permute indexes for advancing to save number of registers used
        __m512i perm_rNrNp1_0_31lo_w = perm_r0r1_0_31lo;
        __m512i perm_rNrNp1_0_31hi_w = perm_r0r1_0_31hi;

        if constexpr (partial_load) {
          // Safe masked loads for the image edge
          data_src = _mm512_maskz_loadu_epi16(k1, src_ptr);
          data_src2 = _mm512_maskz_loadu_epi16(k2, src_ptr + 32);
        }
        else {
          // Fast unaligned loads for the safe zone
          data_src = _mm512_loadu_si512(src_ptr);
          data_src2 = _mm512_loadu_si512(src_ptr + 32);
        }

        // rows 0..3
        __m512i src_r0r1_0_31lo = _mm512_permutex2var_epi16(data_src, perm_rNrNp1_0_31lo_w, data_src2);
        __m512i src_r0r1_0_31hi = _mm512_permutex2var_epi16(data_src, perm_rNrNp1_0_31hi_w, data_src2);

        // for r2r3
        perm_rNrNp1_0_31lo_w = _mm512_add_epi16(perm_rNrNp1_0_31lo_w, two_epi16);
        perm_rNrNp1_0_31hi_w = _mm512_add_epi16(perm_rNrNp1_0_31hi_w, two_epi16);

        __m512i src_r2r3_0_31lo = _mm512_permutex2var_epi16(data_src, perm_rNrNp1_0_31lo_w, data_src2);
        __m512i src_r2r3_0_31hi = _mm512_permutex2var_epi16(data_src, perm_rNrNp1_0_31hi_w, data_src2);

        // for r4r5
        perm_rNrNp1_0_31lo_w = _mm512_add_epi16(perm_rNrNp1_0_31lo_w, two_epi16);
        perm_rNrNp1_0_31hi_w = _mm512_add_epi16(perm_rNrNp1_0_31hi_w, two_epi16);

        __m512i result_0_31lo, result_0_31hi;

        if constexpr (!lessthan16bit) {
          // madd requires signed integers, so shift to signed range
          src_r0r1_0_31lo = _mm512_add_epi16(src_r0r1_0_31lo, shifttosigned);
          src_r0r1_0_31hi = _mm512_add_epi16(src_r0r1_0_31hi, shifttosigned);
          src_r2r3_0_31lo = _mm512_add_epi16(src_r2r3_0_31lo, shifttosigned);
          src_r2r3_0_31hi = _mm512_add_epi16(src_r2r3_0_31hi, shifttosigned);
        }

        if constexpr (UseVNNI)
        {
          result_0_31lo = _mm512_dpwssd_epi32(rounder, src_r0r1_0_31lo, coef_r0r1_0_31lo);
          result_0_31lo = _mm512_dpwssd_epi32(result_0_31lo, src_r2r3_0_31lo, coef_r2r3_0_31lo);

          result_0_31hi = _mm512_dpwssd_epi32(rounder, src_r0r1_0_31hi, coef_r0r1_0_31hi);
          result_0_31hi = _mm512_dpwssd_epi32(result_0_31hi, src_r2r3_0_31hi, coef_r2r3_0_31hi);
        }
        else
        {
          // making FMA in 32bits accs as in AVX256 V-resize
          result_0_31lo = _mm512_add_epi32(_mm512_madd_epi16(src_r0r1_0_31lo, coef_r0r1_0_31lo), _mm512_madd_epi16(src_r2r3_0_31lo, coef_r2r3_0_31lo));
          result_0_31hi = _mm512_add_epi32(_mm512_madd_epi16(src_r0r1_0_31hi, coef_r0r1_0_31hi), _mm512_madd_epi16(src_r2r3_0_31hi, coef_r2r3_0_31hi));
        }

        // rows 4..7
        __m512i src_r4r5_0_31lo = _mm512_permutex2var_epi16(data_src, perm_rNrNp1_0_31lo_w, data_src2);
        __m512i src_r4r5_0_31hi = _mm512_permutex2var_epi16(data_src, perm_rNrNp1_0_31hi_w, data_src2);

        // for r6r7
        perm_rNrNp1_0_31lo_w = _mm512_add_epi16(perm_rNrNp1_0_31lo_w, two_epi16);
        perm_rNrNp1_0_31hi_w = _mm512_add_epi16(perm_rNrNp1_0_31hi_w, two_epi16);

        __m512i src_r6r7_0_31lo = _mm512_permutex2var_epi16(data_src, perm_rNrNp1_0_31lo_w, data_src2);
        __m512i src_r6r7_0_31hi = _mm512_permutex2var_epi16(data_src, perm_rNrNp1_0_31hi_w, data_src2);

        // for r8r9
        perm_rNrNp1_0_31lo_w = _mm512_add_epi16(perm_rNrNp1_0_31lo_w, two_epi16);
        perm_rNrNp1_0_31hi_w = _mm512_add_epi16(perm_rNrNp1_0_31hi_w, two_epi16);

        if constexpr (!lessthan16bit) {
          // madd requires signed integers, so shift to signed range
          src_r4r5_0_31lo = _mm512_add_epi16(src_r4r5_0_31lo, shifttosigned);
          src_r4r5_0_31hi = _mm512_add_epi16(src_r4r5_0_31hi, shifttosigned);
          src_r6r7_0_31lo = _mm512_add_epi16(src_r6r7_0_31lo, shifttosigned);
          src_r6r7_0_31hi = _mm512_add_epi16(src_r6r7_0_31hi, shifttosigned);
        }

        if constexpr (UseVNNI)
        {
          result_0_31lo = _mm512_dpwssd_epi32(result_0_31lo, src_r4r5_0_31lo, coef_r4r5_0_31lo);
          result_0_31lo = _mm512_dpwssd_epi32(result_0_31lo, src_r6r7_0_31lo, coef_r6r7_0_31lo);

          result_0_31hi = _mm512_dpwssd_epi32(result_0_31hi, src_r4r5_0_31hi, coef_r4r5_0_31hi);
          result_0_31hi = _mm512_dpwssd_epi32(result_0_31hi, src_r6r7_0_31hi, coef_r6r7_0_31hi);
        }
        else
        {
          // making FMA in 32bits accs as in AVX256 V-resize
          result_0_31lo = _mm512_add_epi32(_mm512_madd_epi16(src_r4r5_0_31lo, coef_r4r5_0_31lo), result_0_31lo);
          result_0_31hi = _mm512_add_epi32(_mm512_madd_epi16(src_r4r5_0_31hi, coef_r4r5_0_31hi), result_0_31hi);

          result_0_31lo = _mm512_add_epi32(_mm512_madd_epi16(src_r6r7_0_31lo, coef_r6r7_0_31lo), result_0_31lo);
          result_0_31hi = _mm512_add_epi32(_mm512_madd_epi16(src_r6r7_0_31hi, coef_r6r7_0_31hi), result_0_31hi);
        }

        // rows 8..11
        __m512i src_r8r9_0_31lo = _mm512_permutex2var_epi16(data_src, perm_rNrNp1_0_31lo_w, data_src2);
        __m512i src_r8r9_0_31hi = _mm512_permutex2var_epi16(data_src, perm_rNrNp1_0_31hi_w, data_src2);

        // for r10r11
        perm_rNrNp1_0_31lo_w = _mm512_add_epi16(perm_rNrNp1_0_31lo_w, two_epi16);
        perm_rNrNp1_0_31hi_w = _mm512_add_epi16(perm_rNrNp1_0_31hi_w, two_epi16);

        __m512i src_r10r11_0_31lo = _mm512_permutex2var_epi16(data_src, perm_rNrNp1_0_31lo_w, data_src2);
        __m512i src_r10r11_0_31hi = _mm512_permutex2var_epi16(data_src, perm_rNrNp1_0_31hi_w, data_src2);

        // for r12r13
        perm_rNrNp1_0_31lo_w = _mm512_add_epi16(perm_rNrNp1_0_31lo_w, two_epi16);
        perm_rNrNp1_0_31hi_w = _mm512_add_epi16(perm_rNrNp1_0_31hi_w, two_epi16);

        if constexpr (!lessthan16bit) {
          // madd requires signed integers, so shift to signed range
          src_r8r9_0_31lo = _mm512_add_epi16(src_r8r9_0_31lo, shifttosigned);
          src_r8r9_0_31hi = _mm512_add_epi16(src_r8r9_0_31hi, shifttosigned);
          src_r10r11_0_31lo = _mm512_add_epi16(src_r10r11_0_31lo, shifttosigned);
          src_r10r11_0_31hi = _mm512_add_epi16(src_r10r11_0_31hi, shifttosigned);
        }

        if constexpr (UseVNNI)
        {
          result_0_31lo = _mm512_dpwssd_epi32(result_0_31lo, src_r8r9_0_31lo, coef_r8r9_0_31lo);
          result_0_31lo = _mm512_dpwssd_epi32(result_0_31lo, src_r10r11_0_31lo, coef_r10r11_0_31lo);

          result_0_31hi = _mm512_dpwssd_epi32(result_0_31hi, src_r8r9_0_31hi, coef_r8r9_0_31hi);
          result_0_31hi = _mm512_dpwssd_epi32(result_0_31hi, src_r10r11_0_31hi, coef_r10r11_0_31hi);
        }
        else
        {
          // making FMA in 32bits accs as in AVX256 V-resize
          result_0_31lo = _mm512_add_epi32(_mm512_madd_epi16(src_r8r9_0_31lo, coef_r8r9_0_31lo), result_0_31lo);
          result_0_31hi = _mm512_add_epi32(_mm512_madd_epi16(src_r8r9_0_31hi, coef_r8r9_0_31hi), result_0_31hi);

          result_0_31lo = _mm512_add_epi32(_mm512_madd_epi16(src_r10r11_0_31lo, coef_r10r11_0_31lo), result_0_31lo);
          result_0_31hi = _mm512_add_epi32(_mm512_madd_epi16(src_r10r11_0_31hi, coef_r10r11_0_31hi), result_0_31hi);
        }

        // rows 12..15
        __m512i src_r12r13_0_31lo = _mm512_permutex2var_epi16(data_src, perm_rNrNp1_0_31lo_w, data_src2);
        __m512i src_r12r13_0_31hi = _mm512_permutex2var_epi16(data_src, perm_rNrNp1_0_31hi_w, data_src2);

        // for r14r15
        perm_rNrNp1_0_31lo_w = _mm512_add_epi16(perm_rNrNp1_0_31lo_w, two_epi16);
        perm_rNrNp1_0_31hi_w = _mm512_add_epi16(perm_rNrNp1_0_31hi_w, two_epi16);

        __m512i src_r14r15_0_31lo = _mm512_permutex2var_epi16(data_src, perm_rNrNp1_0_31lo_w, data_src2);
        __m512i src_r14r15_0_31hi = _mm512_permutex2var_epi16(data_src, perm_rNrNp1_0_31hi_w, data_src2);

        if constexpr (!lessthan16bit) {
          // madd requires signed integers, so shift to signed range
          src_r12r13_0_31lo = _mm512_add_epi16(src_r12r13_0_31lo, shifttosigned);
          src_r12r13_0_31hi = _mm512_add_epi16(src_r12r13_0_31hi, shifttosigned);
          src_r14r15_0_31lo = _mm512_add_epi16(src_r14r15_0_31lo, shifttosigned);
          src_r14r15_0_31hi = _mm512_add_epi16(src_r14r15_0_31hi, shifttosigned);
        }

        if constexpr (UseVNNI)
        {
          result_0_31lo = _mm512_dpwssd_epi32(result_0_31lo, src_r12r13_0_31lo, coef_r12r13_0_31lo);
          result_0_31lo = _mm512_dpwssd_epi32(result_0_31lo, src_r14r15_0_31lo, coef_r14r15_0_31lo);

          result_0_31hi = _mm512_dpwssd_epi32(result_0_31hi, src_r12r13_0_31hi, coef_r12r13_0_31hi);
          result_0_31hi = _mm512_dpwssd_epi32(result_0_31hi, src_r14r15_0_31hi, coef_r14r15_0_31hi);

          // rounding VNNI in first FMA already summed
        }
        else
        {
          result_0_31lo = _mm512_add_epi32(_mm512_madd_epi16(src_r12r13_0_31lo, coef_r12r13_0_31lo), result_0_31lo);
          result_0_31hi = _mm512_add_epi32(_mm512_madd_epi16(src_r12r13_0_31hi, coef_r12r13_0_31hi), result_0_31hi);

          result_0_31lo = _mm512_add_epi32(_mm512_madd_epi16(src_r14r15_0_31lo, coef_r14r15_0_31lo), result_0_31lo);
          result_0_31hi = _mm512_add_epi32(_mm512_madd_epi16(src_r14r15_0_31hi, coef_r14r15_0_31hi), result_0_31hi);

          // rounding
          result_0_31lo = _mm512_add_epi32(result_0_31lo, rounder);
          result_0_31hi = _mm512_add_epi32(result_0_31hi, rounder);
        }

        if constexpr (!lessthan16bit) {
          // return from signed range
          result_0_31lo = _mm512_add_epi32(result_0_31lo, shiftfromsigned);
          result_0_31hi = _mm512_add_epi32(result_0_31hi, shiftfromsigned);
        }

        // scaling down
        result_0_31lo = _mm512_srai_epi32(result_0_31lo, FPScale16bits);
        result_0_31hi = _mm512_srai_epi32(result_0_31hi, FPScale16bits);

        __m512i result_0_31_int16 = _mm512_packus_epi32(result_0_31lo, result_0_31hi);

		result_0_31_int16 = _mm512_min_epu16(result_0_31_int16, clamp_limit_max);
		result_0_31_int16 = _mm512_max_epu16(result_0_31_int16, clamp_limit_min);

        _mm512_stream_si512(reinterpret_cast<__m512i*>(dst_ptr), result_0_31_int16);

        dst_ptr += dst_pitch;
        src_ptr += src_pitch;
      }

      current_coeff += filter_size * PIXELS_AT_A_TIME;
      };

    // Process the 'safe zone' where direct full unaligned loads are acceptable.
    for (; x < width_safe_mod; x += PIXELS_AT_A_TIME)
    {
      do_h_integer_core(std::false_type{});
    }

    // Process the potentially 'unsafe zone' near the image edge, using safe masked loading.
    for (; x < width; x += PIXELS_AT_A_TIME)
    {
      do_h_integer_core(std::true_type{});
    }
  }
}
