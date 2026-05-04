// Avisynth v2.5.  Copyright 2002 Ben Rudiak-Gould et al.
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

// VS 2017 v15.3
#if _MSC_VER >= 1911


#include "./avs/alignment.h"
#include "./avs/minmax.h"

#include "./resample_avx512.h"
#include <type_traits>

#include <immintrin.h> // Includes AVX-512 intrinsics

#if _MSC_VER < 1922 // Check for MSVC version less than 16.2 (VS 2019 16.2)
  // Define missing AVX-512BW mask intrinsics for older MSVC.
  // inline functions that perform the mask operations directly.
  // Since this is MSVC only, using specific __forceinline.
__forceinline __mmask64 _kand_mask64(__mmask64 a, __mmask64 b) { return a & b; }
__forceinline __mmask64 _kor_mask64(__mmask64 a, __mmask64 b) { return a | b; }
__forceinline __mmask32 _kand_mask32(__mmask32 a, __mmask32 b) { return a & b; }
__forceinline __mmask32 _kor_mask32(__mmask32 a, __mmask32 b) { return a | b; }
#endif

#include "./resample_avx512.hpp"

void resize_h_planar_uint8_avx512_permutex_vstripe_ks4_vbmi(BYTE* dst8, const BYTE* src8, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height, int bits_per_pixel,const uint8_t range,const bool mode_YUY2)
{
  // true template parameter: VBMI version
  resize_h_planar_uint8_avx512_permutex_vstripe_ks4_internal<true>(dst8, src8, dst_pitch, src_pitch, program, width, height, bits_per_pixel, range, mode_YUY2);
}

void resize_h_planar_uint8_avx512_permutex_vstripe_ks8_vbmi(BYTE* dst8, const BYTE* src8, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height, int bits_per_pixel,const uint8_t range,const bool mode_YUY2)
{
  // true template parameter: VBMI version
  resize_h_planar_uint8_avx512_permutex_vstripe_ks8_internal<true>(dst8, src8, dst_pitch, src_pitch, program, width, height, bits_per_pixel, range, mode_YUY2);
}

void resize_h_planar_uint8_avx512_permutex_vstripe_2s32_ks8_vbmi(BYTE* dst8, const BYTE* src8, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height, int bits_per_pixel,const uint8_t range,const bool mode_YUY2)
{
  // true template parameter: VBMI version
  resize_h_planar_uint8_avx512_permutex_vstripe_2s32_ks8_internal<true>(dst8, src8, dst_pitch, src_pitch, program, width, height, bits_per_pixel, range, mode_YUY2);
}

void resize_h_planar_uint8_avx512_permutex_vstripe_ks16_vbmi(BYTE* dst8, const BYTE* src8, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height, int bits_per_pixel,const uint8_t range,const bool mode_YUY2)
{
  // true template parameter: VBMI version
  resize_h_planar_uint8_avx512_permutex_vstripe_ks16_internal<true>(dst8, src8, dst_pitch, src_pitch, program, width, height, bits_per_pixel, range, mode_YUY2);
}

// uint8_t h "mpz" avx512base 4,8,16

void resize_h_planar_uint8_avx512_permutex_vstripe_mpz_ks4_vnni(BYTE* dst8, const BYTE* src8, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height, int bits_per_pixel,const uint8_t range,const bool mode_YUY2)
{
  // template parameter true: UseVNNI
  resize_h_planar_uint8_avx512_permutex_vstripe_mpz_ks4_internal<true>(dst8, src8, dst_pitch, src_pitch, program, width, height, bits_per_pixel, range, mode_YUY2);
}
void resize_h_planar_uint8_avx512_permutex_vstripe_mpz_ks8_vnni(BYTE* dst8, const BYTE* src8, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height, int bits_per_pixel,const uint8_t range,const bool mode_YUY2)
{
  // template parameter true: UseVNNI
  resize_h_planar_uint8_avx512_permutex_vstripe_mpz_ks8_internal<true>(dst8, src8, dst_pitch, src_pitch, program, width, height, bits_per_pixel, range, mode_YUY2);
}
void resize_h_planar_uint8_avx512_permutex_vstripe_mpz_ks16_vnni(BYTE* dst8, const BYTE* src8, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height, int bits_per_pixel,const uint8_t range,const bool mode_YUY2)
{
  // template parameter true: UseVNNI
  resize_h_planar_uint8_avx512_permutex_vstripe_mpz_ks16_internal<true>(dst8, src8, dst_pitch, src_pitch, program, width, height, bits_per_pixel, range, mode_YUY2);
}


// uint16_t h "mp" avx512fast vnni 4,8,16

template<bool lessthan16bit>
void resize_h_planar_uint16_avx512_permutex_vstripe_mp_ks4_vnni(BYTE* dst8, const BYTE* src8, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height, int bits_per_pixel,const uint8_t range,const bool mode_YUY2)
{
  // true template parameter: VNNI version
  resize_h_planar_uint16_avx512_permutex_vstripe_mp_ks4_internal<lessthan16bit, true>(dst8, src8, dst_pitch, src_pitch, program, width, height, bits_per_pixel, range, mode_YUY2);
}
template<bool lessthan16bit>
void resize_h_planar_uint16_avx512_permutex_vstripe_mp_ks8_vnni(BYTE* dst8, const BYTE* src8, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height, int bits_per_pixel,const uint8_t range,const bool mode_YUY2)
{
  // true template parameter: VNNI version
  resize_h_planar_uint16_avx512_permutex_vstripe_mp_ks8_internal<lessthan16bit, true>(dst8, src8, dst_pitch, src_pitch, program, width, height, bits_per_pixel, range, mode_YUY2);
}
template<bool lessthan16bit>
void resize_h_planar_uint16_avx512_permutex_vstripe_mp_ks16_vnni(BYTE* dst8, const BYTE* src8, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height, int bits_per_pixel,const uint8_t range,const bool mode_YUY2)
{
  // true template parameter: VNNI version
  resize_h_planar_uint16_avx512_permutex_vstripe_mp_ks16_internal<lessthan16bit, true>(dst8, src8, dst_pitch, src_pitch, program, width, height, bits_per_pixel, range, mode_YUY2);
}

// Explicit template instantiations
template void resize_h_planar_uint16_avx512_permutex_vstripe_mp_ks4_vnni<false>(BYTE* dst8, const BYTE* src8, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height, int bits_per_pixel,const uint8_t range,const bool mode_YUY2);
template void resize_h_planar_uint16_avx512_permutex_vstripe_mp_ks4_vnni<true>(BYTE* dst8, const BYTE* src8, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height, int bits_per_pixel,const uint8_t range,const bool mode_YUY2);
template void resize_h_planar_uint16_avx512_permutex_vstripe_mp_ks8_vnni<false>(BYTE* dst8, const BYTE* src8, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height, int bits_per_pixel,const uint8_t range,const bool mode_YUY2);
template void resize_h_planar_uint16_avx512_permutex_vstripe_mp_ks8_vnni<true>(BYTE* dst8, const BYTE* src8, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height, int bits_per_pixel,const uint8_t range,const bool mode_YUY2);
template void resize_h_planar_uint16_avx512_permutex_vstripe_mp_ks16_vnni<false>(BYTE* dst8, const BYTE* src8, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height, int bits_per_pixel,const uint8_t range,const bool mode_YUY2);
template void resize_h_planar_uint16_avx512_permutex_vstripe_mp_ks16_vnni<true>(BYTE* dst8, const BYTE* src8, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height, int bits_per_pixel,const uint8_t range,const bool mode_YUY2);

#endif
