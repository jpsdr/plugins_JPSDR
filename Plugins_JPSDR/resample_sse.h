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

#ifndef __Resample_SSE_H__
#define __Resample_SSE_H__

// Intrinsics for SSE4.1, SSSE3, SSE3, SSE2, ISSE and MMX
#include <emmintrin.h>
#include <smmintrin.h>

#include "avisynth.h"
#include "./avs/config.h"
#include "./resample.h"


#ifdef X86_32
void resize_v_mmx_planar(BYTE* dst8, const BYTE* src8, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int bits_per_pixel, int MinY, int MaxY, const int* pitch_table, const void* storage,const uint8_t range,const bool mode_YUY2);
#endif

void resize_v_sse2_planar(BYTE* dst8, const BYTE* src8, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int bits_per_pixel, int MinY, int MaxY, const int* pitch_table, const void* storage,const uint8_t range,const bool mode_YUY2);

template<bool lessthan16bit>
void resize_v_sse2_planar_uint16_t(BYTE* dst8, const BYTE* src8, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int bits_per_pixel, int MinY, int MaxY, const int* pitch_table, const void* storage,const uint8_t range,const bool mode_YUY2);

void resize_v_sse2_planar_float(BYTE* dst8, const BYTE* src8, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int bits_per_pixel, int MinY, int MaxY, const int* pitch_table, const void* storage,const uint8_t range,const bool mode_YUY2);

// =============================================================================

template<typename pixel_t, bool lessthan16bit>
#if defined(GCC) || defined(CLANG)
__attribute__((__target__("ssse3")))
#endif
void resizer_h_ssse3_generic_uint8_16(BYTE* dst8, const BYTE* src8, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height, int bits_per_pixel,const uint8_t range,const bool mode_YUY2);

#if defined(GCC) || defined(CLANG)
__attribute__((__target__("ssse3")))
#endif
void resizer_h_ssse3_generic_float(BYTE* dst8, const BYTE* src8, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height, int bits_per_pixel,const uint8_t range,const bool mode_YUY2);

#endif // __Resample_SSE_H__
