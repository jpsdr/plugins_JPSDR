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

#include <stdio.h>
#include "resample.h"
#include "./avs/config.h"
#include "./avs/alignment.h"

#define myfree(ptr) if (ptr!=NULL) { free(ptr); ptr=NULL;}
#define myCloseHandle(ptr) if (ptr!=NULL) { CloseHandle(ptr); ptr=NULL;}
#define mydelete(ptr) if (ptr!=NULL) { delete ptr; ptr=NULL;}
#define mydelete2(ptr) if (ptr!=NULL) { delete[] ptr; ptr=NULL;}

#include <type_traits>
// Intrinsics for SSE4.1, SSSE3, SSE3, SSE2, ISSE and MMX
#include <emmintrin.h>
#include <smmintrin.h>
#include <algorithm>

extern ThreadPoolInterface *poolInterface;

#if _MSC_VER >= 1900
  #define AVX_BUILD_POSSIBLE
#endif

/***************************************
 ********* Templated SSE Loader ********
 ***************************************/

typedef __m128i (SSELoader)(const __m128i*);
typedef __m128 (SSELoader_ps)(const float*);

__forceinline __m128i simd_load_aligned(const __m128i* adr)
{
  return _mm_load_si128(adr);
}

__forceinline __m128i simd_load_unaligned(const __m128i* adr)
{
  return _mm_loadu_si128(adr);
}

__forceinline __m128i simd_load_unaligned_sse3(const __m128i* adr)
{
  return _mm_lddqu_si128(adr);
}

__forceinline __m128i simd_load_streaming(const __m128i* adr)
{
  return _mm_stream_load_si128(const_cast<__m128i*>(adr));
}

// float loaders
__forceinline __m128 simd_loadps_aligned(const float * adr)
{
  return _mm_load_ps(adr);
}

__forceinline __m128 simd_loadps_unaligned(const float* adr)
{
  return _mm_loadu_ps(adr);
}


// fake _mm_packus_epi32 (orig is SSE4.1 only)
__forceinline __m128i _MM_PACKUS_EPI32( __m128i a, __m128i b )
{
  a = _mm_slli_epi32 (a, 16);
  a = _mm_srai_epi32 (a, 16);
  b = _mm_slli_epi32 (b, 16);
  b = _mm_srai_epi32 (b, 16);
  a = _mm_packs_epi32 (a, b);
  return a;
}

/***************************************
 ***** Vertical Resizer Assembly *******
 ***************************************/

template<typename pixel_t>
static void resize_v_planar_pointresize(BYTE* dst, const BYTE* src, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int bits_per_pixel, int MinY, int MaxY, const int* pitch_table, const void* storage,const uint8_t range,const bool mode_YUY2)
{
  pixel_t *src0 = (pixel_t *)src;
  pixel_t *dst0 = (pixel_t *)dst;
  dst_pitch/=sizeof(pixel_t);

  for (int y = MinY; y < MaxY; y++)
  {
	const pixel_t *src_ptr = src0 + pitch_table[program->pixel_offset[y]];
    
	memcpy(dst0,src_ptr,width*sizeof(pixel_t));

    dst0+=dst_pitch;
  }
}


static void resize_v_c_planar(BYTE* dst, const BYTE* src, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int bits_per_pixel, int MinY, int MaxY, const int* pitch_table, const void* storage,const uint8_t range,const bool mode_YUY2)
{
	const int filter_size = program->filter_size;
	const short *current_coeff = program->pixel_coefficient+filter_size*MinY;

	const int val_min = (range==1) ? 0 : 16;
	const int val_max = ((range==1) || (range==4)) ? 255 : (range==2) ? 235 : 240;

	if ((mode_YUY2) && ((range>=2) && (range<=3)))
	{
		const int TabMax[4] = {235,240,235,240};

		for (int y = MinY; y < MaxY; y++)
		{
			const BYTE *src_ptr = src + pitch_table[program->pixel_offset[y]];

			for (int x = 0; x < width; x++)
			{
				int result = 0;

				for (int i = 0; i < filter_size; i++)
					result += (src_ptr+pitch_table[i])[x] * current_coeff[i];

				result = (result+8192) >> 14;
				result = (result>TabMax[x & 0x03]) ? TabMax[x & 0x3] : (result<16) ? 16 : result;
				dst[x] = (BYTE) result;
			}

			dst += dst_pitch;
			current_coeff += filter_size;
		}
	}
	else
	{
		for (int y = MinY; y < MaxY; y++)
		{
			const BYTE *src_ptr = src + pitch_table[program->pixel_offset[y]];

			for (int x = 0; x < width; x++)
			{
				int result = 0;

				for (int i = 0; i < filter_size; i++)
					result += (src_ptr+pitch_table[i])[x] * current_coeff[i];

				result = (result+8192) >> 14;
				result = (result>val_max) ? val_max : (result<val_min) ? val_min : result;
				dst[x] = (BYTE) result;
			}

			dst += dst_pitch;
			current_coeff += filter_size;
		}
	}
}


static void resize_v_c_planar_f(BYTE* dst, const BYTE* src, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int bits_per_pixel, int MinY, int MaxY, const int* pitch_table, const void* storage,const uint8_t range,const bool mode_YUY2)
{
  const int filter_size = program->filter_size;
  const float *current_coeff = program->pixel_coefficient_float+filter_size*MinY;

  const float *src0 = (float *)src;
  float *dst0 = (float *)dst;

  dst_pitch>>=2;

  for (int y = MinY; y < MaxY; y++)
  {
	const float *src_ptr = src0 + pitch_table[program->pixel_offset[y]];

    for (int x = 0; x < width; x++)
	{
      float result = 0;

      for (int i = 0; i < filter_size; i++)
		result += (src_ptr+pitch_table[i])[x] * current_coeff[i];

      dst0[x] = result;
    }

    dst0 += dst_pitch;
    current_coeff += filter_size;
  }
}


static void resize_v_c_planar_s(BYTE* dst, const BYTE* src, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int bits_per_pixel, int MinY, int MaxY, const int* pitch_table, const void* storage,const uint8_t range,const bool mode_YUY2)
{
	const int filter_size = program->filter_size;
	const short *current_coeff = program->pixel_coefficient+filter_size*MinY;

	const uint16_t *src0 = (uint16_t *)src;
	uint16_t *dst0 = (uint16_t *)dst;
	const __int64 val_min = (range==1) ? 0 : (int)16 << (bits_per_pixel-8);
	const __int64 val235 = (int)235 << (bits_per_pixel-8);
	const __int64 val240 = (int)240 << (bits_per_pixel-8);
	const __int64 val_max = ((range==1) || (range==4)) ? ((int)1 << bits_per_pixel)-1 : (range==2) ? val235 : val240;

	dst_pitch>>=1;

	for (int y = MinY; y < MaxY; y++)
	{
		const uint16_t *src_ptr = src0 + pitch_table[program->pixel_offset[y]];

		for (int x = 0; x < width; x++)
		{
			__int64 result = 0;

			for (int i = 0; i < filter_size; i++)
				result += (src_ptr+pitch_table[i])[x] * current_coeff[i];

			result = (result+8192) >> 14;
			result = (result>val_max) ? val_max : (result<val_min) ? val_min : result;
			dst0[x] = (uint16_t) result;
		}
		dst0 += dst_pitch;
		current_coeff += filter_size;
	}
}


#ifdef X86_32
static void resize_v_mmx_planar(BYTE* dst, const BYTE* src, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int bits_per_pixel, int MinY, int MaxY, const int* pitch_table, const void* storage,const uint8_t range,const bool mode_YUY2)
{
  const int filter_size = program->filter_size;
  const short *current_coeff = program->pixel_coefficient + filter_size*MinY;

  const int wMod8 = (width >> 3) << 3;
  const int sizeMod2 = (filter_size >> 1) << 1;
  const bool notMod2 = sizeMod2 < filter_size;

  const __m64 zero = _mm_setzero_si64();

	const int val_min = (range==1) ? 0 : 16;
	const int val_max = ((range==1) || (range==4)) ? 255 : (range==2) ? 235 : 240;
	const int TabMax[4] = {235,240,235,240};

	const __m64 val_min_m64 = _mm_set1_pi16((short)((val_min << 8)|val_min));
	const __m64 val_max_m64 = (mode_YUY2 && ((range>=2) && (range<=3))) ? _mm_set1_pi16((short)(((int)240 << 8)|235)) : _mm_set1_pi16((short)((val_max << 8)|val_max));

  
  for (int y = MinY; y < MaxY; y++)
  {
    const BYTE *src_ptr = src + pitch_table[program->pixel_offset[y]];

    for (int x = 0; x < wMod8; x += 8)
	{
      __m64 result_1 = _mm_set1_pi32(8192); // Init. with rounder (16384/2 = 8192)
      __m64 result_2 = result_1;
      __m64 result_3 = result_1;
      __m64 result_4 = result_1;

      for (int i = 0; i < sizeMod2; i += 2)
	  {
        __m64 src_p1 = *(reinterpret_cast<const __m64*>(src_ptr+pitch_table[i]+x));   // For detailed explanation please see SSE2 version.
        __m64 src_p2 = *(reinterpret_cast<const __m64*>(src_ptr+pitch_table[i+1]+x));

        __m64 src_l = _mm_unpacklo_pi8(src_p1, src_p2);                                   
        __m64 src_h = _mm_unpackhi_pi8(src_p1, src_p2);                                   

        __m64 src_1 = _mm_unpacklo_pi8(src_l, zero);                                      
        __m64 src_2 = _mm_unpackhi_pi8(src_l, zero);                                      
        __m64 src_3 = _mm_unpacklo_pi8(src_h, zero);                                      
        __m64 src_4 = _mm_unpackhi_pi8(src_h, zero);                                      

        __m64 coeff = _mm_cvtsi32_si64(*reinterpret_cast<const int*>(current_coeff+i));   
        coeff = _mm_unpacklo_pi32(coeff, coeff);                                               

        __m64 dst_1 = _mm_madd_pi16(src_1, coeff);                                        
        __m64 dst_2 = _mm_madd_pi16(src_2, coeff);                                        
        __m64 dst_3 = _mm_madd_pi16(src_3, coeff);
        __m64 dst_4 = _mm_madd_pi16(src_4, coeff);

        result_1 = _mm_add_pi32(result_1, dst_1);
        result_2 = _mm_add_pi32(result_2, dst_2);
        result_3 = _mm_add_pi32(result_3, dst_3);
        result_4 = _mm_add_pi32(result_4, dst_4);
      }

      if (notMod2)
	  { // do last odd row
        __m64 src_p = *(reinterpret_cast<const __m64*>(src_ptr+pitch_table[sizeMod2]+x));

        __m64 src_l = _mm_unpacklo_pi8(src_p, zero);
        __m64 src_h = _mm_unpackhi_pi8(src_p, zero);

        __m64 coeff = _mm_set1_pi16(current_coeff[sizeMod2]);

        __m64 dst_ll = _mm_mullo_pi16(src_l, coeff);   // Multiply by coefficient
        __m64 dst_lh = _mm_mulhi_pi16(src_l, coeff);
        __m64 dst_hl = _mm_mullo_pi16(src_h, coeff);
        __m64 dst_hh = _mm_mulhi_pi16(src_h, coeff);

        __m64 dst_1 = _mm_unpacklo_pi16(dst_ll, dst_lh); // Unpack to 32-bit integer
        __m64 dst_2 = _mm_unpackhi_pi16(dst_ll, dst_lh);
        __m64 dst_3 = _mm_unpacklo_pi16(dst_hl, dst_hh);
        __m64 dst_4 = _mm_unpackhi_pi16(dst_hl, dst_hh);

        result_1 = _mm_add_pi32(result_1, dst_1);
        result_2 = _mm_add_pi32(result_2, dst_2);
        result_3 = _mm_add_pi32(result_3, dst_3);
        result_4 = _mm_add_pi32(result_4, dst_4);
      }

      // Divide by 16348 (FPRound)
      result_1  = _mm_srai_pi32(result_1, 14);
      result_2  = _mm_srai_pi32(result_2, 14);
      result_3  = _mm_srai_pi32(result_3, 14);
      result_4  = _mm_srai_pi32(result_4, 14);

      // Pack and store
      __m64 result_l = _mm_packs_pi32(result_1, result_2);
      __m64 result_h = _mm_packs_pi32(result_3, result_4);
      __m64 result   = _mm_packs_pu16(result_l, result_h);

		result = _mm_max_pu8(result,val_min_m64);
		result = _mm_min_pu8(result,val_max_m64);

      *(reinterpret_cast<__m64*>(dst+x)) = result;
    }

    // Leftover
	if ((mode_YUY2) && ((range>=2) && (range<=3)))
	{
		for (int x = wMod8; x < width; x++)
		{
			int result = 0;

		    for (int i = 0; i < filter_size; i++)
				result += (src_ptr+pitch_table[i])[x] * current_coeff[i];
			result = (result+8192) >> 14;
			result = (result>TabMax[x & 0x03]) ? TabMax[x & 0x03] : (result<val_min) ? val_min : result;
			dst[x] = (BYTE) result;
		}
	}
	else
	{
		for (int x = wMod8; x < width; x++)
		{
			int result = 0;

		    for (int i = 0; i < filter_size; i++)
				result += (src_ptr+pitch_table[i])[x] * current_coeff[i];
			result = (result+8192) >> 14;
			result = (result>val_max) ? val_max : (result<val_min) ? val_min : result;
			dst[x] = (BYTE) result;
		}
	}
	dst += dst_pitch;
	current_coeff += filter_size;

  }

  _mm_empty();
}
#endif

template<SSELoader load>
static void resize_v_sse2_planar(BYTE* dst, const BYTE* src, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int bits_per_pixel, int MinY, int MaxY, const int* pitch_table, const void* storage,const uint8_t range,const bool mode_YUY2)
{
  const int filter_size = program->filter_size;
  const short *current_coeff = program->pixel_coefficient + filter_size*MinY;
  
  const int wMod16 = (width >> 4) << 4;
  const int sizeMod2 = (filter_size >> 1) << 1;
  const bool notMod2 = sizeMod2 < filter_size;

  const __m128i zero = _mm_setzero_si128();


	const int val_min = (range==1) ? 0 : 16;
	const int val_max = ((range==1) || (range==4)) ? 255 : (range==2) ? 235 : 240;
	const int TabMax[4] = {235,240,235,240};

	const __m128i val_min_m128 = _mm_set1_epi16((short)((val_min << 8)|val_min));
	const __m128i val_max_m128 = (mode_YUY2 && ((range>=2) && (range<=3))) ? _mm_set1_epi16((short)(((int)240 << 8)|235)) : _mm_set1_epi16((short)((val_max << 8)|val_max));

  for (int y = MinY; y < MaxY; y++)
  {
    const BYTE *src_ptr = src + pitch_table[program->pixel_offset[y]];

    for (int x = 0; x < wMod16; x += 16)
	{
      __m128i result_1 = _mm_set1_epi32(8192); // Init. with rounder (16384/2 = 8192)
      __m128i result_2 = result_1;
      __m128i result_3 = result_1;
      __m128i result_4 = result_1;
      
      for (int i = 0; i < sizeMod2; i += 2)
	  {
        __m128i src_p1 = load(reinterpret_cast<const __m128i*>(src_ptr+pitch_table[i]+x));   // p|o|n|m|l|k|j|i|h|g|f|e|d|c|b|a
        __m128i src_p2 = load(reinterpret_cast<const __m128i*>(src_ptr+pitch_table[i+1]+x)); // P|O|N|M|L|K|J|I|H|G|F|E|D|C|B|A
         
        __m128i src_l = _mm_unpacklo_epi8(src_p1, src_p2);                                   // Hh|Gg|Ff|Ee|Dd|Cc|Bb|Aa
        __m128i src_h = _mm_unpackhi_epi8(src_p1, src_p2);                                   // Pp|Oo|Nn|Mm|Ll|Kk|Jj|Ii

        __m128i src_1 = _mm_unpacklo_epi8(src_l, zero);                                      // .D|.d|.C|.c|.B|.b|.A|.a
        __m128i src_2 = _mm_unpackhi_epi8(src_l, zero);                                      // .H|.h|.G|.g|.F|.f|.E|.e
        __m128i src_3 = _mm_unpacklo_epi8(src_h, zero);                                      // etc.
        __m128i src_4 = _mm_unpackhi_epi8(src_h, zero);                                      // etc.

        __m128i coeff = _mm_cvtsi32_si128(*reinterpret_cast<const int*>(current_coeff+i));   // XX|XX|XX|XX|XX|XX|CO|co
        coeff = _mm_shuffle_epi32(coeff, 0);                                                 // CO|co|CO|co|CO|co|CO|co
        
        __m128i dst_1 = _mm_madd_epi16(src_1, coeff);                                         // CO*D+co*d | CO*C+co*c | CO*B+co*b | CO*A+co*a
        __m128i dst_2 = _mm_madd_epi16(src_2, coeff);                                         // etc.
        __m128i dst_3 = _mm_madd_epi16(src_3, coeff);
        __m128i dst_4 = _mm_madd_epi16(src_4, coeff);

        result_1 = _mm_add_epi32(result_1, dst_1);
        result_2 = _mm_add_epi32(result_2, dst_2);
        result_3 = _mm_add_epi32(result_3, dst_3);
        result_4 = _mm_add_epi32(result_4, dst_4);
      }
      
      if (notMod2)
	  { // do last odd row
        __m128i src_p = load(reinterpret_cast<const __m128i*>(src_ptr+pitch_table[sizeMod2]+x));

        __m128i src_l = _mm_unpacklo_epi8(src_p, zero);
        __m128i src_h = _mm_unpackhi_epi8(src_p, zero);

        __m128i coeff = _mm_set1_epi16(current_coeff[sizeMod2]);

        __m128i dst_ll = _mm_mullo_epi16(src_l, coeff);   // Multiply by coefficient
        __m128i dst_lh = _mm_mulhi_epi16(src_l, coeff);
        __m128i dst_hl = _mm_mullo_epi16(src_h, coeff);
        __m128i dst_hh = _mm_mulhi_epi16(src_h, coeff);

        __m128i dst_1 = _mm_unpacklo_epi16(dst_ll, dst_lh); // Unpack to 32-bit integer
        __m128i dst_2 = _mm_unpackhi_epi16(dst_ll, dst_lh);
        __m128i dst_3 = _mm_unpacklo_epi16(dst_hl, dst_hh);
        __m128i dst_4 = _mm_unpackhi_epi16(dst_hl, dst_hh);

        result_1 = _mm_add_epi32(result_1, dst_1);
        result_2 = _mm_add_epi32(result_2, dst_2);
        result_3 = _mm_add_epi32(result_3, dst_3);
        result_4 = _mm_add_epi32(result_4, dst_4);
      }
      
      // Divide by 16348 (FPRound)
      result_1  = _mm_srai_epi32(result_1, 14);
      result_2  = _mm_srai_epi32(result_2, 14);
      result_3  = _mm_srai_epi32(result_3, 14);
      result_4  = _mm_srai_epi32(result_4, 14);

      // Pack and store
      __m128i result_l = _mm_packs_epi32(result_1, result_2);
      __m128i result_h = _mm_packs_epi32(result_3, result_4);
      __m128i result   = _mm_packus_epi16(result_l, result_h);

		result = _mm_max_epu8(result,val_min_m128);
		result = _mm_min_epu8(result,val_max_m128);

      _mm_store_si128(reinterpret_cast<__m128i*>(dst+x), result);
    }

    // Leftover
	if ((mode_YUY2) && ((range>=2) && (range<=3)))
	{
		for (int x = wMod16; x < width; x++)
		{
			int result = 0;

			for (int i = 0; i < filter_size; i++)
				result += (src_ptr+pitch_table[i])[x] * current_coeff[i];
			result = (result+8192) >> 14;
			result = (result>TabMax[x & 0x03]) ? TabMax[x & 0x03] : (result<val_min) ? val_min : result;
			dst[x] = (BYTE) result;
		}
	}
	else
	{
		for (int x = wMod16; x < width; x++)
		{
			int result = 0;

			for (int i = 0; i < filter_size; i++)
				result += (src_ptr+pitch_table[i])[x] * current_coeff[i];
			result = (result+8192) >> 14;
			result = (result>val_max) ? val_max : (result<val_min) ? val_min : result;
			dst[x] = (BYTE) result;
		}
	}

    dst += dst_pitch;
    current_coeff += filter_size;
  }
}


// for uint16_t and float. Both uses float arithmetic and coefficients
template<SSELoader_ps loadps>
static void resize_v_sseX_planar_32(BYTE* dst0, const BYTE* src0, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int bits_per_pixel, int MinY, int MaxY, const int* pitch_table, const void* storage,const uint8_t range,const bool mode_YUY2)
{
  const int filter_size = program->filter_size;
  const float *current_coeff_float = program->pixel_coefficient_float  + filter_size*MinY;

  const int wMod8 = (width >> 3) << 3; // uint16/float: 8 at a time (byte was 16 byte at a time)

  const float *src = (float *)src0;
  float *dst = (float *)dst0;
  dst_pitch >>= 2;
  src_pitch >>= 2;
  
  for (int y = MinY; y < MaxY; y++)
  {
    const float *src_ptr = src + pitch_table[program->pixel_offset[y]];

    for (int x = 0; x < wMod8; x+=8)
	{
      __m128 result_l_single = _mm_set1_ps(0.0f);
      __m128 result_h_single = result_l_single;

      const float *src2_ptr = src_ptr+x;

      for (int i = 0; i < filter_size; i++)
	  {
        __m128 src_l_single;
        __m128 src_h_single;
		
        // _mm_load_ps or _mm_loadu_ps template dependent
        src_l_single = loadps(reinterpret_cast<const float*>(src2_ptr)); // float  4*32=128 4 pixels at a time
        src_h_single = loadps(reinterpret_cast<const float*>(src2_ptr+4));
		  
        __m128 coeff = _mm_load1_ps(reinterpret_cast<const float*>(current_coeff_float+i)); // loads 1, fills all 4 floats
        __m128 dst_l = _mm_mul_ps(src_l_single, coeff); // Multiply by coefficient
        __m128 dst_h = _mm_mul_ps(src_h_single, coeff); // 4*(32bit*32bit=32bit)
        result_l_single = _mm_add_ps(result_l_single, dst_l); // accumulate result.
        result_h_single = _mm_add_ps(result_h_single, dst_h);

        src2_ptr += src_pitch;
      }

      _mm_store_ps(reinterpret_cast<float*>(dst+x), result_l_single);
      _mm_store_ps(reinterpret_cast<float*>(dst+x+4), result_h_single);
    }

    // Leftover
    for (int x = wMod8; x < width; x++)
	{
      float result = 0;
	  
      for (int i = 0; i < filter_size; i++)
        result += (src_ptr+pitch_table[i])[x] * current_coeff_float[i];	
      dst[x] =  result;
    }

    dst += dst_pitch;
    current_coeff_float += filter_size;
  }
}


// for uint16_t and float. Both uses float arithmetic and coefficients
template<SSELoader load, bool sse41>
static void resize_v_sseX_planar_16(BYTE* dst0, const BYTE* src0, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int bits_per_pixel, int MinY, int MaxY, const int* pitch_table, const void* storage,const uint8_t range,const bool mode_YUY2)
{
  const int filter_size = program->filter_size;
  const float *current_coeff_float = program->pixel_coefficient_float  + filter_size*MinY;

  const int wMod8 = (width >> 3) << 3; // uint16/float: 8 at a time (byte was 16 byte at a time)

  const __m128i zero = _mm_setzero_si128();

  const uint16_t *src = (uint16_t *)src0;
  uint16_t *dst = (uint16_t *)dst0;
  dst_pitch >>= 1;
  src_pitch >>= 1;

	const int val_min = (range==1) ? 0 : (int)16 << (bits_per_pixel-8);
	const int val_235 = (int)235 << (bits_per_pixel-8);
	const int val_240 = (int)240 << (bits_per_pixel-8);
	const int val_max = ((range==1) || (range==4)) ? ((int)1 << bits_per_pixel)-1 : (range==2) ? val_235 : val_240;
	const float val_minf = (float)val_min;
	const float val_maxf = (float)val_max;

  const __m128 val_min_128f = _mm_set1_ps(val_minf); // clamp limit
  const __m128 val_max_128f = _mm_set1_ps(val_maxf); // clamp limit
  const __m128i val_min_128 = _mm_set1_epi16((short)val_min);
  const __m128i val_max_128 = _mm_set1_epi16((short)val_max);
	  
  for (int y = MinY; y < MaxY; y++)
  {
    const uint16_t *src_ptr = src + pitch_table[program->pixel_offset[y]];

    for (int x = 0; x < wMod8; x+=8)
	{
      __m128 result_l_single = _mm_set1_ps(0.0f);
      __m128 result_h_single = result_l_single;

      const uint16_t *src2_ptr = src_ptr+x;

      for (int i = 0; i < filter_size; i++)
	  {
        __m128 src_l_single;
        __m128 src_h_single;
		
        __m128i src_p = load(reinterpret_cast<const __m128i*>(src2_ptr)); // uint16_t  8*16=128 8 pixels at a time
        __m128i src_l = _mm_unpacklo_epi16(src_p, zero); // spread lower  4*uint16_t pixel value -> 4*32 bit
        __m128i src_h = _mm_unpackhi_epi16(src_p, zero); // spread higher 4*uint16_t pixel value -> 4*32 bit
        src_l_single = _mm_cvtepi32_ps (src_l); // Converts the four signed 32-bit integer values of a to single-precision, floating-point values.
        src_h_single = _mm_cvtepi32_ps (src_h);
		  
        __m128 coeff = _mm_load1_ps(reinterpret_cast<const float*>(current_coeff_float+i)); // loads 1, fills all 4 floats
        __m128 dst_l = _mm_mul_ps(src_l_single, coeff); // Multiply by coefficient
        __m128 dst_h = _mm_mul_ps(src_h_single, coeff); // 4*(32bit*32bit=32bit)
        result_l_single = _mm_add_ps(result_l_single, dst_l); // accumulate result.
        result_h_single = _mm_add_ps(result_h_single, dst_h);

        src2_ptr += src_pitch;
      }

      // clamp!
	  if (!sse41)
	  {
		  result_l_single = _mm_min_ps(result_l_single, val_max_128f);  // mainly for 10-14 bit 
		  result_h_single = _mm_min_ps(result_h_single, val_max_128f); // mainly for 10-14 bit
		  result_l_single = _mm_max_ps(result_l_single, val_min_128f);  // mainly for 10-14 bit 
		  result_h_single = _mm_max_ps(result_h_single, val_min_128f); // mainly for 10-14 bit
		  // result = _mm_max_ps(result, zero); low limit through pack_us
          // Converts the four single-precision, floating-point values of a to signed 32-bit integer values.		  
		  __m128i result_l  = _mm_cvtps_epi32(result_l_single);
		  __m128i result_h  = _mm_cvtps_epi32(result_h_single);
		  //result_l = _mm_min_epi32(result_l, clamp_limit_i32); // mainly for 10-14 bit 
          //result_h = _mm_min_epi32(result_h, clamp_limit_i32); // mainly for 10-14 bit 
                                                                 // Pack and store
		  __m128i result = (_MM_PACKUS_EPI32(result_l, result_h)); // 4*32+4*32 = 8*16
		  _mm_stream_si128(reinterpret_cast<__m128i*>(dst+x), result);
		  
	  }
	  else
	  {
		  // result = _mm_max_ps(result, zero); low limit through pack_us
          // Converts the four single-precision, floating-point values of a to signed 32-bit integer values.		  
		  __m128i result_l  = _mm_cvtps_epi32(result_l_single);
		  __m128i result_h  = _mm_cvtps_epi32(result_h_single);
		  //result_l = _mm_min_epi32(result_l, clamp_limit_i32); // mainly for 10-14 bit 
          //result_h = _mm_min_epi32(result_h, clamp_limit_i32); // mainly for 10-14 bit 
                                                                 // Pack and store
		  __m128i result = _mm_packus_epi32(result_l, result_h);
		  result = _mm_min_epu16(result,val_max_128); // unsigned clamp here
		  result = _mm_max_epu16(result,val_min_128);
		  _mm_stream_si128(reinterpret_cast<__m128i*>(dst+x), result);
		  
	  }
	}

    // Leftover
	for (int x = wMod8; x < width; x++)
	{
		float result = 0;
	  
		for (int i = 0; i < filter_size; i++)
			result += (src_ptr+pitch_table[i])[x] * current_coeff_float[i];
		result = (result>val_maxf) ? val_maxf : (result<val_minf) ? val_minf : result;
		dst[x] = (uint16_t)result;
	}

    dst += dst_pitch;
    current_coeff_float += filter_size;
  }
}


template<SSELoader load>
static void resize_v_ssse3_planar(BYTE* dst, const BYTE* src, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int bits_per_pixel, int MinY, int MaxY, const int* pitch_table, const void* storage,const uint8_t range,const bool mode_YUY2)
{
  const int filter_size = program->filter_size;
  const short *current_coeff = program->pixel_coefficient + filter_size*MinY;
  
  const int wMod16 = (width >> 4) << 4;

	const int val_min = (range==1) ? 0 : 16;
	const int val_max = ((range==1) || (range==4)) ? 255 : (range==2) ? 235 : 240;
	const int TabMax[4] = {235,240,235,240};

	const __m128i val_min_m128 = _mm_set1_epi16((short)((val_min << 8)|val_min));
	const __m128i val_max_m128 = (mode_YUY2 && ((range>=2) && (range<=3))) ? _mm_set1_epi16((short)(((int)240 << 8)|235)) : _mm_set1_epi16((short)((val_max << 8)|val_max));

  const __m128i zero = _mm_setzero_si128();
  const __m128i coeff_unpacker = _mm_set_epi8(1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0);

  for (int y = MinY; y < MaxY; y++)
  {
    const BYTE* src_ptr = src + pitch_table[program->pixel_offset[y]];

    for (int x = 0; x < wMod16; x+=16)
	{
      __m128i result_l = _mm_set1_epi16(32); // Init. with rounder ((1 << 6)/2 = 32)
      __m128i result_h = result_l;

      const BYTE* src2_ptr = src_ptr+x;
      
      for (int i = 0; i < filter_size; i++)
	  {
        __m128i src_p = load(reinterpret_cast<const __m128i*>(src2_ptr));

        __m128i src_l = _mm_unpacklo_epi8(src_p, zero);
        __m128i src_h = _mm_unpackhi_epi8(src_p, zero);

        src_l = _mm_slli_epi16(src_l, 7);
        src_h = _mm_slli_epi16(src_h, 7);

        __m128i coeff = _mm_cvtsi32_si128(*reinterpret_cast<const int*>(current_coeff+i));
                coeff = _mm_shuffle_epi8(coeff, coeff_unpacker);

        __m128i dst_l = _mm_mulhrs_epi16(src_l, coeff);   // Multiply by coefficient (SSSE3)
        __m128i dst_h = _mm_mulhrs_epi16(src_h, coeff);

        result_l = _mm_add_epi16(result_l, dst_l);
        result_h = _mm_add_epi16(result_h, dst_h);

        src2_ptr += src_pitch;
      }

      // Divide by 64
      result_l  = _mm_srai_epi16(result_l, 6);
      result_h  = _mm_srai_epi16(result_h, 6);

      // Pack and store
      __m128i result   = _mm_packus_epi16(result_l, result_h);

		result = _mm_max_epu8(result,val_min_m128);
		result = _mm_min_epu8(result,val_max_m128);

      _mm_store_si128(reinterpret_cast<__m128i*>(dst+x), result);
    }

    // Leftover
	if ((mode_YUY2) && ((range>=2) && (range<=3)))
	{
		for (int x = wMod16; x < width; x++)
		{
			int result = 0;

			for (int i = 0; i < filter_size; i++)
				result += (src_ptr+pitch_table[i])[x] * current_coeff[i];
			result = (result+8192) >> 14;
			result = (result>TabMax[x & 0x03]) ? TabMax[x & 0x03] : (result<val_min) ? val_min : result;
			dst[x] = (BYTE) result;
		}
	}
	else
	{
		for (int x = wMod16; x < width; x++)
		{
			int result = 0;

			for (int i = 0; i < filter_size; i++)
				result += (src_ptr+pitch_table[i])[x] * current_coeff[i];
			result = (result+8192) >> 14;
			result = (result>val_max) ? val_max : (result<val_min) ? val_min : result;
			dst[x] = (BYTE) result;
		}
	}

    dst += dst_pitch;
    current_coeff += filter_size;
  }
}


#ifdef AVX_BUILD_POSSIBLE

// for uint16_t and float. Both uses float arithmetic and coefficients
// see the same in resample_avx2
void resize_v_avx_planar_32(BYTE* dst0, const BYTE* src0, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int bits_per_pixel, int MinY, int MaxY, const int* pitch_table, const void* storage,const uint8_t range,const bool mode_YUY2)
{
  _mm256_zeroupper();
  
  const int filter_size = program->filter_size;
  const float *current_coeff_float = program->pixel_coefficient_float + filter_size*MinY;

  const int wMod8 = (width >> 3) << 3; // uint16/float: 8 at a time (byte was 16 byte at a time)

  const float *src = (float *)src0;
  float *dst = (float *)dst0;
  dst_pitch>>=2;
  src_pitch>>=2;

  for (int y = MinY; y < MaxY; y++)
  {
    const float *src_ptr = src + pitch_table[program->pixel_offset[y]];

    for (int x = 0; x < wMod8; x+=8)
	{
      __m256 result_single = _mm256_set1_ps(0.0f);

      const float *src2_ptr = src_ptr+x;

      for (int i = 0; i < filter_size; i++)
	  {
        __m256 src_single;
		
        // float
        // avx solution is chosen is pointers are aligned
        __m128 src_l_single = _mm_load_ps(reinterpret_cast<const float*>(src2_ptr));   // float  4*32=128 4 pixels at a time
        __m128 src_h_single = _mm_load_ps(reinterpret_cast<const float*>(src2_ptr+4)); // float  4*32=128 4 pixels at a time
        src_single = _mm256_set_m128(src_h_single, src_l_single);
        // using one 256bit load instead of 2x128bit is slower on avx-only Ivy
        //src_single = _mm256_load_ps(reinterpret_cast<const float*>(src2_ptr)); // float  8*32=256 8 pixels at a time
		
        __m256 coeff = _mm256_broadcast_ss(reinterpret_cast<const float*>(current_coeff_float+i)); // loads 1, fills all 8 floats
        __m256 dst = _mm256_mul_ps(src_single, coeff); // Multiply by coefficient // 8*(32bit*32bit=32bit)
        result_single = _mm256_add_ps(result_single, dst); // accumulate result.

        src2_ptr += src_pitch;
      }

      // float
      _mm256_stream_ps(reinterpret_cast<float*>(dst+x), result_single);
    }

    // Leftover
    for (int x = wMod8; x < width; x++)
	{
      float result = 0;
      for (int i = 0; i < filter_size; i++)
        result += (src_ptr+pitch_table[i])[x] * current_coeff_float[i];
      dst[x] = (float) result;
    }

    dst += dst_pitch;
    current_coeff_float += filter_size;
  }
  _mm256_zeroupper();
}


// for uint16_t and float. Both uses float arithmetic and coefficients
// see the same in resample_avx2
void resize_v_avx_planar_16(BYTE* dst0, const BYTE* src0, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int bits_per_pixel, int MinY, int MaxY, const int* pitch_table, const void* storage,const uint8_t range,const bool mode_YUY2)
{
  _mm256_zeroupper();
  
  const int filter_size = program->filter_size;
  const float *current_coeff_float = program->pixel_coefficient_float + filter_size*MinY;

  const int wMod8 = (width >> 3) << 3; // uint16/float: 8 at a time (byte was 16 byte at a time)

  const __m128i zero = _mm_setzero_si128();

  const uint16_t *src = (uint16_t *)src0;
  uint16_t *dst = (uint16_t *)dst0;
  dst_pitch>>=1;
  src_pitch>>=1;

	const int val_min = (range==1) ? 0 : (int)16 << (bits_per_pixel-8);
	const int val_235 = (int)235 << (bits_per_pixel-8);
	const int val_240 = (int)240 << (bits_per_pixel-8);
	const int val_max = ((range==1) || (range==4)) ? ((int)1 << bits_per_pixel)-1 : (range==2) ? val_235 : val_240;
	const float val_minf = (float)val_min;
	const float val_maxf = (float)val_max;

  const __m128i val_min_128 = _mm_set1_epi16((short)val_min);
  const __m128i val_max_128 = _mm_set1_epi16((short)val_max);

  for (int y = MinY; y < MaxY; y++)
  {
    const uint16_t *src_ptr = src + pitch_table[program->pixel_offset[y]];

    for (int x = 0; x < wMod8; x+=8)
	{
      __m256 result_single = _mm256_set1_ps(0.0f);

      const uint16_t *src2_ptr = src_ptr+x;

      for (int i = 0; i < filter_size; i++)
	  {
        __m256 src_single;
		
        // word
        // avx solution is chosen is pointers are aligned
        __m128i src_p = _mm_load_si128(reinterpret_cast<const __m128i*>(src2_ptr)); // uint16_t  8*16=128 8 pixels at a time
        __m256i src256;
		
        // simple avx
        __m128i src_l = _mm_unpacklo_epi16(src_p, zero); // spread lower  4*uint16_t pixel value -> 4*32 bit
        __m128i src_h = _mm_unpackhi_epi16(src_p, zero); // spread higher 4*uint16_t pixel value -> 4*32 bit
        src256 = _mm256_set_m128i(src_h, src_l);
		  
        src_single = _mm256_cvtepi32_ps(src256); // Converts the eight signed 32-bit integer values of avx to single-precision, floating-point values.
		
        __m256 coeff = _mm256_broadcast_ss(reinterpret_cast<const float*>(current_coeff_float+i)); // loads 1, fills all 8 floats
        __m256 dst = _mm256_mul_ps(src_single, coeff); // Multiply by coefficient // 8*(32bit*32bit=32bit)
        result_single = _mm256_add_ps(result_single, dst); // accumulate result.

        src2_ptr += src_pitch;
      }

      // word
      // clamp! no! later at uint16 stage
      // result_single = _mm256_min_ps(result_single, clamp_limit_256); // mainly for 10-14 bit 
      // result = _mm_max_ps(result, zero); low limit through pack_us
      // Converts the 8 single-precision, floating-point values of a to signed 32-bit integer values.
      __m256i result256  = _mm256_cvtps_epi32(result_single);
      // Pack and store
      __m128i result = _mm_packus_epi32(_mm256_extractf128_si256(result256, 0), _mm256_extractf128_si256(result256, 1)); // 4*32+4*32 = 8*16
      result = _mm_min_epu16(result,val_max_128); // unsigned clamp here
	  result = _mm_max_epu16(result,val_min_128);
      _mm_stream_si128(reinterpret_cast<__m128i*>(dst+x), result);
    }

    // Leftover
	for (int x = wMod8; x < width; x++)
	{
		float result = 0;

		for (int i = 0; i < filter_size; i++)
			result += (src_ptr+pitch_table[i])[x] * current_coeff_float[i];
		result = (result>val_maxf) ? val_maxf : (result<val_minf) ? val_minf : result;
		dst[x] = (uint16_t) result;
	}

    dst += dst_pitch;
    current_coeff_float += filter_size;
  }
  _mm256_zeroupper();
}


// for uint16_t and float. Both uses float arithmetic and coefficients
// see the same in resample_avx
void resize_v_avx2_planar_32(BYTE* dst0, const BYTE* src0, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int bits_per_pixel, int MinY, int MaxY, const int* pitch_table, const void* storage,const uint8_t range,const bool mode_YUY2)
{
  _mm256_zeroupper();
  
  const int filter_size = program->filter_size;
  const float *current_coeff_float = program->pixel_coefficient_float + filter_size*MinY;

  const int wMod8 = (width >> 3) << 3; // uint16/float: 8 at a time (byte was 16 byte at a time)

  const float *src = (float *)src0;
  float *dst = (float *)dst0;
  dst_pitch>>=2;
  src_pitch>>=2;

  for (int y = MinY; y < MaxY; y++)
  {
    const float *src_ptr = src + pitch_table[program->pixel_offset[y]];

    for (int x = 0; x < wMod8; x+=8)
	{
      __m256 result_single = _mm256_set1_ps(0.0f);

      const float *src2_ptr = src_ptr+x;

      for (int i = 0; i < filter_size; i++)
	  {
        __m256 src_single;
		
        // float
        // avx solution is chosen is pointers are aligned
        __m128 src_l_single = _mm_load_ps(reinterpret_cast<const float*>(src2_ptr));   // float  4*32=128 4 pixels at a time
        __m128 src_h_single = _mm_load_ps(reinterpret_cast<const float*>(src2_ptr+4)); // float  4*32=128 4 pixels at a time
        src_single = _mm256_set_m128(src_h_single, src_l_single);
        // using one 256bit load instead of 2x128bit is slower on avx-only Ivy
        //src_single = _mm256_load_ps(reinterpret_cast<const float*>(src2_ptr)); // float  8*32=256 8 pixels at a time
		
        __m256 coeff = _mm256_broadcast_ss(reinterpret_cast<const float*>(current_coeff_float+i)); // loads 1, fills all 8 floats
        __m256 dst = _mm256_mul_ps(src_single, coeff); // Multiply by coefficient // 8*(32bit*32bit=32bit)
        result_single = _mm256_add_ps(result_single, dst); // accumulate result.

        src2_ptr += src_pitch;
      }

      // float
      _mm256_stream_ps(reinterpret_cast<float*>(dst+x), result_single);
    }

    // Leftover
    for (int x = wMod8; x < width; x++)
	{
      float result = 0;
      for (int i = 0; i < filter_size; i++)
        result += (src_ptr+pitch_table[i])[x] * current_coeff_float[i];
      dst[x] = (float) result;
    }

    dst += dst_pitch;
    current_coeff_float += filter_size;
  }
  _mm256_zeroupper();
}


// for uint16_t and float. Both uses float arithmetic and coefficients
// see the same in resample_avx
void resize_v_avx2_planar_16(BYTE* dst0, const BYTE* src0, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int bits_per_pixel, int MinY, int MaxY, const int* pitch_table, const void* storage,const uint8_t range,const bool mode_YUY2)
{
  _mm256_zeroupper();
  
  const int filter_size = program->filter_size;
  const float *current_coeff_float = program->pixel_coefficient_float + filter_size*MinY;

  const int wMod8 = (width >> 3) << 3; // uint16/float: 8 at a time (byte was 16 byte at a time)

  const __m128i zero = _mm_setzero_si128();

  const uint16_t *src = (uint16_t *)src0;
  uint16_t *dst = (uint16_t *)dst0;
  dst_pitch>>=1;
  src_pitch>>=1;

	const int val_min = (range==1) ? 0 : (int)16 << (bits_per_pixel-8);
	const int val_235 = (int)235 << (bits_per_pixel-8);
	const int val_240 = (int)240 << (bits_per_pixel-8);
	const int val_max = ((range==1) || (range==4)) ? ((int)1 << bits_per_pixel)-1 : (range==2) ? val_235 : val_240;
	const float val_minf = (float)val_min;
	const float val_maxf = (float)val_max;

  const __m128i val_min_128 = _mm_set1_epi16((short)val_min);
  const __m128i val_max_128 = _mm_set1_epi16((short)val_max);

  for (int y = MinY; y < MaxY; y++)
  {
    const uint16_t *src_ptr = src + pitch_table[program->pixel_offset[y]];

    for (int x = 0; x < wMod8; x+=8)
	{
      __m256 result_single = _mm256_set1_ps(0.0f);

      const uint16_t *src2_ptr = src_ptr+x;

      for (int i = 0; i < filter_size; i++)
	  {
        __m256 src_single;
		
        // word
        // avx solution is chosen is pointers are aligned
        __m128i src_p = _mm_load_si128(reinterpret_cast<const __m128i*>(src2_ptr)); // uint16_t  8*16=128 8 pixels at a time
        __m256i src256;
		
        // AVX2: 
        //_mm256_unpacklo is not good, because it works on the 2 * lower_64_bit of the two 128bit halves
        src256 = _mm256_cvtepu16_epi32(src_p);
		
        src_single = _mm256_cvtepi32_ps(src256); // Converts the eight signed 32-bit integer values of avx to single-precision, floating-point values.
		
        __m256 coeff = _mm256_broadcast_ss(reinterpret_cast<const float*>(current_coeff_float+i)); // loads 1, fills all 8 floats
        __m256 dst = _mm256_mul_ps(src_single, coeff); // Multiply by coefficient // 8*(32bit*32bit=32bit)
        result_single = _mm256_add_ps(result_single, dst); // accumulate result.

        src2_ptr += src_pitch;
      }

      // word
      // clamp! no! later at uint16 stage
      // result_single = _mm256_min_ps(result_single, clamp_limit_256); // mainly for 10-14 bit 
      // result = _mm_max_ps(result, zero); low limit through pack_us
      // Converts the 8 single-precision, floating-point values of a to signed 32-bit integer values.
      __m256i result256  = _mm256_cvtps_epi32(result_single);
      // Pack and store
      __m128i result = _mm_packus_epi32(_mm256_extractf128_si256(result256, 0), _mm256_extractf128_si256(result256, 1)); // 4*32+4*32 = 8*16
      result = _mm_min_epu16(result,val_max_128); // unsigned clamp here
	  result = _mm_max_epu16(result,val_min_128);
      _mm_stream_si128(reinterpret_cast<__m128i*>(dst+x), result);
    }

    // Leftover
	for (int x = wMod8; x < width; x++)
	{
		float result = 0;

		for (int i = 0; i < filter_size; i++)
			result += (src_ptr+pitch_table[i])[x] * current_coeff_float[i];
		result = (result>val_maxf) ? val_maxf : (result<val_minf) ? val_minf : result;
		dst[x] = (uint16_t) result;
	}

    dst += dst_pitch;
    current_coeff_float += filter_size;
  }
  _mm256_zeroupper();
}


#endif


__forceinline static void resize_v_create_pitch_table(int* table, int pitch, int height, uint8_t pixel_size)
{
  switch(pixel_size)
  {
	case 2 : pitch>>=1; break;
	case 4 : pitch>>=2; break;
	default : ;
  }
  table[0] = 0;
  for (int i = 1; i < height; i++)
    table[i] = table[i-1]+pitch;
}


/***************************************
 ********* Horizontal Resizer** ********
 ***************************************/

static void resize_h_pointresize(BYTE* dst, const BYTE* src, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height, int bits_per_pixel,const uint8_t range,const bool mode_YUY2)
{
  const int wMod4 = (width >> 2) << 2;

  for (int y = 0; y < height; y++)
  {
    for (int x = 0; x < wMod4; x+=4)
	{
#define pixel(a) src[program->pixel_offset[x+a]]
      unsigned int data = (pixel(3) << 24) + (pixel(2) << 16) + (pixel(1) << 8) + pixel(0);
#undef pixel
      *((unsigned int *)(dst+x)) = data;
    }

    for (int x = wMod4; x < width; x++)
      dst[x] = src[program->pixel_offset[x]];

    dst += dst_pitch;
    src += src_pitch;
  }
}

// make the resampling coefficient array mod8 friendly for simd, padding non-used coeffs with zeros
static void resize_h_prepare_coeff_8(ResamplingProgram* p,IScriptEnvironment* env)
{
  const int filter_size = AlignNumber(p->filter_size, 8);
  const int im0=p->target_size;
  short *new_coeff = (short*) _aligned_malloc(sizeof(short) * im0 * filter_size, 64);
  float *new_coeff_float = (float*) _aligned_malloc(sizeof(float) * im0 * filter_size, 64);
  if ((new_coeff==NULL) || (new_coeff_float==NULL))
  {
	myalignedfree(new_coeff_float);
    myalignedfree(new_coeff);
    env->ThrowError("ResizeHMT: Could not reserve memory in a resampler.");
  }

  memset(new_coeff, 0, sizeof(short) * im0 * filter_size);
  const int im=im0*filter_size;
  for (int i=0; i<im; i++)
	  new_coeff_float[i]=0.0f;
  
  // Copy coeff
  short *dst = new_coeff, *src = p->pixel_coefficient;
  float *dst_f = new_coeff_float, *src_f = p->pixel_coefficient_float;
  
  const int im1=p->filter_size;
  
  for (int i = 0; i < im0; i++)
  {
    for (int j = 0; j < im1; j++)
	{
      dst[j] = src[j];
      dst_f[j] = src_f[j];
    }

    dst += filter_size;
    src += im1;
	dst_f += filter_size;
	src_f += im1;
  }

  myalignedfree(p->pixel_coefficient_float);
  myalignedfree(p->pixel_coefficient);
  p->pixel_coefficient = new_coeff;
  p->pixel_coefficient_float = new_coeff_float;
}


static void resize_h_c_planar(BYTE* dst, const BYTE* src, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height, int bits_per_pixel,const uint8_t range,const bool mode_YUY2)
{
  const int filter_size = program->filter_size;
  int y_src_pitch=0,y_dst_pitch=0;
  
	const int val_min = (range==1) ? 0 : 16;
	const int val_max = ((range==1) || (range==4)) ? 255 : (range==2) ? 235 : 240;

  // external loop y is much faster

	if ((mode_YUY2) && ((range>=2) && (range<=3)))
	{
		const int TabMax[4] = {235,240,235,240};

		for (int y = 0; y < height; y++)
		{
			const short *current_coeff=program->pixel_coefficient;
	  
			for (int x = 0; x < width; x++)
			{
				const int begin = program->pixel_offset[x];
				int result = 0;
		  
				for (int i = 0; i < filter_size; i++)
	    			result+=(src+y_src_pitch)[(begin+i)]*current_coeff[i];
		
				result = (result + 8192) >> 14;
				result = (result>TabMax[x & 0x03]) ? TabMax[x & 0x03] : (result<16) ? 16 : result;
				(dst + y_dst_pitch)[x] = (BYTE)result;		  		  
				current_coeff+=filter_size;
			}
			y_dst_pitch+=dst_pitch;
			y_src_pitch+=src_pitch;	  
		}
	}
	else
	{
		for (int y = 0; y < height; y++)
		{
			const short *current_coeff=program->pixel_coefficient;
	  
			for (int x = 0; x < width; x++)
			{
				const int begin = program->pixel_offset[x];
				int result = 0;
		  
				for (int i = 0; i < filter_size; i++)
	    			result+=(src+y_src_pitch)[(begin+i)]*current_coeff[i];
		
				result = (result + 8192) >> 14;
				result = (result>val_max) ? val_max : (result<val_min) ? val_min : result;
				(dst + y_dst_pitch)[x] = (BYTE)result;		  		  
				current_coeff+=filter_size;
			}
			y_dst_pitch+=dst_pitch;
			y_src_pitch+=src_pitch;	  
		}
	}
 
}


static void resize_h_c_planar_s(BYTE* dst, const BYTE* src, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height, int bits_per_pixel,const uint8_t range,const bool mode_YUY2)
{
  const int filter_size = program->filter_size;
  int y_src_pitch=0,y_dst_pitch=0;
  const __int64 limit=(1 << bits_per_pixel) - 1;
  const uint16_t *src0 = (uint16_t *)src;
  uint16_t *dst0 = (uint16_t *)dst;
	const __int64 val_min = (range==1) ? 0 : (int)16 << (bits_per_pixel-8);
	const __int64 val235 = (int)235 << (bits_per_pixel-8);
	const __int64 val240 = (int)240 << (bits_per_pixel-8);
	const __int64 val_max = ((range==1) || (range==4)) ? ((int)1 << bits_per_pixel)-1 : (range==2) ? val235 : val240;

  src_pitch>>=1;
  dst_pitch>>=1;
  
	for (int y = 0; y < height; y++)
	{
		const short *current_coeff=program->pixel_coefficient;
  
		for (int x = 0; x < width; x++)
		{
			const int begin = program->pixel_offset[x];
			__int64 result = 0;
		  
			for (int i = 0; i < filter_size; i++)
				result+=(src0+y_src_pitch)[(begin+i)]*current_coeff[i];
		  
			result = (result + 8192) >> 14;
			result = (result>val_max) ? val_max : (result<val_min) ? val_min : result;
			(dst0 + y_dst_pitch)[x] = (uint16_t)result;
			current_coeff+=filter_size;
		}
		y_dst_pitch+=dst_pitch;
		y_src_pitch+=src_pitch;
	}
}


static void resize_h_c_planar_f(BYTE* dst, const BYTE* src, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height, int bits_per_pixel,const uint8_t range,const bool mode_YUY2)
{
  const int filter_size = program->filter_size;
  int y_src_pitch=0,y_dst_pitch=0;
  const float *src0=(float *)src;
  float *dst0=(float *)dst;

  src_pitch>>=2;
  dst_pitch>>=2;

  for (int y = 0; y < height; y++)
  {
	  const float *current_coeff=program->pixel_coefficient_float;
	  
	  for (int x = 0; x < width; x++)
	  {
		  const int begin = program->pixel_offset[x];
		  float result = 0;
		  
		  for (int i = 0; i < filter_size; i++)
			  result+=(src0+y_src_pitch)[(begin+i)]*current_coeff[i];
		  
		  (dst0 + y_dst_pitch)[x] = result;
		  current_coeff+=filter_size;
	  }
	  y_dst_pitch+=dst_pitch;
	  y_src_pitch+=src_pitch;
  }
}


#if 1
static void resizer_h_ssse3_generic_int16_float_32(BYTE* dst8, const BYTE* src8, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height, int bits_per_pixel,const uint8_t range,const bool mode_YUY2)
{
  const int filter_size = AlignNumber(program->filter_size, 8) >> 3;

  const float *src = reinterpret_cast<const float *>(src8);
  float *dst = reinterpret_cast<float *>(dst8);
  dst_pitch>>=2;
  src_pitch>>=2;

  for (int y = 0; y < height; y++)
  {
    const float *current_coeff = program->pixel_coefficient_float;
	
    for (int x = 0; x < width; x+=4)
	{
      __m128 result1 = _mm_set1_ps(0.0f);
      __m128 result2 = result1;
      __m128 result3 = result1;
      __m128 result4 = result1;

      const int begin1 = program->pixel_offset[x+0];
      const int begin2 = program->pixel_offset[x+1];
      const int begin3 = program->pixel_offset[x+2];
      const int begin4 = program->pixel_offset[x+3];

      // begin1, result1
      for (int i = 0; i < filter_size; i++)
	  {
        __m128 data_l_single, data_h_single;
		
        // float
        // unaligned
        data_l_single = _mm_loadu_ps(reinterpret_cast<const float*>(src+begin1+i*8)); // float  4*32=128 4 pixels at a time
        data_h_single = _mm_loadu_ps(reinterpret_cast<const float*>(src+begin1+i*8+4)); // float  4*32=128 4 pixels at a time
		
        __m128 coeff_l = /*loadps*/_mm_load_ps(reinterpret_cast<const float*>(current_coeff));    // always aligned
        __m128 coeff_h = /*loadps*/_mm_load_ps(reinterpret_cast<const float*>(current_coeff+4));  // always aligned
        __m128 dst_l = _mm_mul_ps(data_l_single, coeff_l); // Multiply by coefficient
        __m128 dst_h = _mm_mul_ps(data_h_single, coeff_h); // 4*(32bit*32bit=32bit)
        result1 = _mm_add_ps(result1, dst_l); // accumulate result.
        result1 = _mm_add_ps(result1, dst_h);

        current_coeff += 8;
      }

      // begin2, result2
      for (int i = 0; i < filter_size; i++)
	  {
        __m128 data_l_single, data_h_single;
		
        // float
        // unaligned
        data_l_single = _mm_loadu_ps(reinterpret_cast<const float*>(src+begin2+i*8)); // float  4*32=128 4 pixels at a time
        data_h_single = _mm_loadu_ps(reinterpret_cast<const float*>(src+begin2+i*8+4)); // float  4*32=128 4 pixels at a time
		
        __m128 coeff_l = /*loadps*/_mm_load_ps(reinterpret_cast<const float*>(current_coeff));    // always aligned
        __m128 coeff_h = /*loadps*/_mm_load_ps(reinterpret_cast<const float*>(current_coeff+4));  // always aligned
        __m128 dst_l = _mm_mul_ps(data_l_single, coeff_l); // Multiply by coefficient
        __m128 dst_h = _mm_mul_ps(data_h_single, coeff_h); // 4*(32bit*32bit=32bit)
        result2 = _mm_add_ps(result2, dst_l); // accumulate result.
        result2 = _mm_add_ps(result2, dst_h);

        current_coeff += 8;
      }

      // begin3, result3
      for (int i = 0; i < filter_size; i++)
	  {
        __m128 data_l_single, data_h_single;
		
        // float
        // unaligned
        data_l_single = _mm_loadu_ps(reinterpret_cast<const float*>(src+begin3+i*8)); // float  4*32=128 4 pixels at a time
        data_h_single = _mm_loadu_ps(reinterpret_cast<const float*>(src+begin3+i*8+4)); // float  4*32=128 4 pixels at a time
		
        __m128 coeff_l = /*loadps*/_mm_load_ps(reinterpret_cast<const float*>(current_coeff));    // always aligned
        __m128 coeff_h = /*loadps*/_mm_load_ps(reinterpret_cast<const float*>(current_coeff+4));  // always aligned
        __m128 dst_l = _mm_mul_ps(data_l_single, coeff_l); // Multiply by coefficient
        __m128 dst_h = _mm_mul_ps(data_h_single, coeff_h); // 4*(32bit*32bit=32bit)
        result3 = _mm_add_ps(result3, dst_l); // accumulate result.
        result3 = _mm_add_ps(result3, dst_h);

        current_coeff += 8;
      }

      // begin4, result4
      for (int i = 0; i < filter_size; i++)
	  {
        __m128 data_l_single, data_h_single;
		
        // float
        // unaligned
        data_l_single = _mm_loadu_ps(reinterpret_cast<const float*>(src+begin4+i*8)); // float  4*32=128 4 pixels at a time
        data_h_single = _mm_loadu_ps(reinterpret_cast<const float*>(src+begin4+i*8+4)); // float  4*32=128 4 pixels at a time
		
        __m128 coeff_l = /*loadps*/_mm_load_ps(reinterpret_cast<const float*>(current_coeff));    // always aligned
        __m128 coeff_h = /*loadps*/_mm_load_ps(reinterpret_cast<const float*>(current_coeff+4));  // always aligned
        __m128 dst_l = _mm_mul_ps(data_l_single, coeff_l); // Multiply by coefficient
        __m128 dst_h = _mm_mul_ps(data_h_single, coeff_h); // 4*(32bit*32bit=32bit)
        result4 = _mm_add_ps(result4, dst_l); // accumulate result.
        result4 = _mm_add_ps(result4, dst_h);

        current_coeff += 8;
      }
      
      __m128 result;

      // this part needs ssse3
      __m128 result12 = _mm_hadd_ps(result1, result2);
      __m128 result34 = _mm_hadd_ps(result3, result4);
      result = _mm_hadd_ps(result12, result34);

      // float
      // aligned
      _mm_store_ps(reinterpret_cast<float*>(dst+x), result); // 4 results at a time
    }

    dst += dst_pitch;
    src += src_pitch;
  }
}

template<bool hasSSE41>
static void resizer_h_ssse3_generic_int16_float_16(BYTE* dst8, const BYTE* src8, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height, int bits_per_pixel,const uint8_t range,const bool mode_YUY2)
{
  const int filter_size = AlignNumber(program->filter_size, 8) >> 3;
  const __m128i zero = _mm_setzero_si128();

  const uint16_t *src = reinterpret_cast<const uint16_t *>(src8);
  uint16_t *dst = reinterpret_cast<uint16_t *>(dst8);
  dst_pitch>>=1;
  src_pitch>>=1;

	const int val_min = (range==1) ? 0 : (int)16 << (bits_per_pixel-8);
	const int val_235 = (int)235 << (bits_per_pixel-8);
	const int val_240 = (int)240 << (bits_per_pixel-8);
	const int val_max = ((range==1) || (range==4)) ? ((int)1 << bits_per_pixel)-1 : (range==2) ? val_235 : val_240;
	const float val_minf = (float)val_min;
	const float val_maxf = (float)val_max;

  const __m128 val_min_128f = _mm_set1_ps(val_minf); // clamp limit
  const __m128 val_max_128f = _mm_set1_ps(val_maxf); // clamp limit

  for (int y = 0; y < height; y++)
  {
    const float *current_coeff = program->pixel_coefficient_float;
	
    for (int x = 0; x < width; x+=4)
	{
      __m128 result1 = _mm_set1_ps(0.0f);
      __m128 result2 = result1;
      __m128 result3 = result1;
      __m128 result4 = result1;

      const int begin1 = program->pixel_offset[x+0];
      const int begin2 = program->pixel_offset[x+1];
      const int begin3 = program->pixel_offset[x+2];
      const int begin4 = program->pixel_offset[x+3];

      // begin1, result1
      for (int i = 0; i < filter_size; i++)
	  {
        __m128 data_l_single, data_h_single;
		
        // unaligned
        __m128i src_p = _mm_loadu_si128(reinterpret_cast<const __m128i*>(src+begin1+(i<<3))); // uint16_t  8*16=128 8 pixels at a time
        __m128i src_l = _mm_unpacklo_epi16(src_p, zero); // spread lower  4*uint16_t pixel value -> 4*32 bit
        __m128i src_h = _mm_unpackhi_epi16(src_p, zero); // spread higher 4*uint16_t pixel value -> 4*32 bit
        data_l_single = _mm_cvtepi32_ps (src_l); // Converts the four signed 32-bit integer values of a to single-precision, floating-point values.
        data_h_single = _mm_cvtepi32_ps (src_h);
		
        __m128 coeff_l = /*loadps*/_mm_load_ps(reinterpret_cast<const float*>(current_coeff));    // always aligned
        __m128 coeff_h = /*loadps*/_mm_load_ps(reinterpret_cast<const float*>(current_coeff+4));  // always aligned
        __m128 dst_l = _mm_mul_ps(data_l_single, coeff_l); // Multiply by coefficient
        __m128 dst_h = _mm_mul_ps(data_h_single, coeff_h); // 4*(32bit*32bit=32bit)
        result1 = _mm_add_ps(result1, dst_l); // accumulate result.
        result1 = _mm_add_ps(result1, dst_h);

        current_coeff += 8;
      }

      // begin2, result2
      for (int i = 0; i < filter_size; i++)
	  {
        __m128 data_l_single, data_h_single;
		
        // unaligned
        __m128i src_p = _mm_loadu_si128(reinterpret_cast<const __m128i*>(src+begin2+(i<<3))); // uint16_t  8*16=128 8 pixels at a time
        __m128i src_l = _mm_unpacklo_epi16(src_p, zero); // spread lower  4*uint16_t pixel value -> 4*32 bit
        __m128i src_h = _mm_unpackhi_epi16(src_p, zero); // spread higher 4*uint16_t pixel value -> 4*32 bit
        data_l_single = _mm_cvtepi32_ps (src_l); // Converts the four signed 32-bit integer values of a to single-precision, floating-point values.
        data_h_single = _mm_cvtepi32_ps (src_h);
		
        __m128 coeff_l = /*loadps*/_mm_load_ps(reinterpret_cast<const float*>(current_coeff));    // always aligned
        __m128 coeff_h = /*loadps*/_mm_load_ps(reinterpret_cast<const float*>(current_coeff+4));  // always aligned
        __m128 dst_l = _mm_mul_ps(data_l_single, coeff_l); // Multiply by coefficient
        __m128 dst_h = _mm_mul_ps(data_h_single, coeff_h); // 4*(32bit*32bit=32bit)
        result2 = _mm_add_ps(result2, dst_l); // accumulate result.
        result2 = _mm_add_ps(result2, dst_h);

        current_coeff += 8;
      }

      // begin3, result3
      for (int i = 0; i < filter_size; i++)
	  {
        __m128 data_l_single, data_h_single;
		
        // unaligned
        __m128i src_p = _mm_loadu_si128(reinterpret_cast<const __m128i*>(src+begin3+(i<<3))); // uint16_t  8*16=128 8 pixels at a time
        __m128i src_l = _mm_unpacklo_epi16(src_p, zero); // spread lower  4*uint16_t pixel value -> 4*32 bit
        __m128i src_h = _mm_unpackhi_epi16(src_p, zero); // spread higher 4*uint16_t pixel value -> 4*32 bit
        data_l_single = _mm_cvtepi32_ps (src_l); // Converts the four signed 32-bit integer values of a to single-precision, floating-point values.
        data_h_single = _mm_cvtepi32_ps (src_h);
		
        __m128 coeff_l = /*loadps*/_mm_load_ps(reinterpret_cast<const float*>(current_coeff));    // always aligned
        __m128 coeff_h = /*loadps*/_mm_load_ps(reinterpret_cast<const float*>(current_coeff+4));  // always aligned
        __m128 dst_l = _mm_mul_ps(data_l_single, coeff_l); // Multiply by coefficient
        __m128 dst_h = _mm_mul_ps(data_h_single, coeff_h); // 4*(32bit*32bit=32bit)
        result3 = _mm_add_ps(result3, dst_l); // accumulate result.
        result3 = _mm_add_ps(result3, dst_h);

        current_coeff += 8;
      }

      // begin4, result4
      for (int i = 0; i < filter_size; i++)
	  {
        __m128 data_l_single, data_h_single;
		
        // unaligned
        __m128i src_p = _mm_loadu_si128(reinterpret_cast<const __m128i*>(src+begin4+(i<<3))); // uint16_t  8*16=128 8 pixels at a time
        __m128i src_l = _mm_unpacklo_epi16(src_p, zero); // spread lower  4*uint16_t pixel value -> 4*32 bit
        __m128i src_h = _mm_unpackhi_epi16(src_p, zero); // spread higher 4*uint16_t pixel value -> 4*32 bit
        data_l_single = _mm_cvtepi32_ps (src_l); // Converts the four signed 32-bit integer values of a to single-precision, floating-point values.
        data_h_single = _mm_cvtepi32_ps (src_h);
		
        __m128 coeff_l = /*loadps*/_mm_load_ps(reinterpret_cast<const float*>(current_coeff));    // always aligned
        __m128 coeff_h = /*loadps*/_mm_load_ps(reinterpret_cast<const float*>(current_coeff+4));  // always aligned
        __m128 dst_l = _mm_mul_ps(data_l_single, coeff_l); // Multiply by coefficient
        __m128 dst_h = _mm_mul_ps(data_h_single, coeff_h); // 4*(32bit*32bit=32bit)
        result4 = _mm_add_ps(result4, dst_l); // accumulate result.
        result4 = _mm_add_ps(result4, dst_h);

        current_coeff += 8;
      }
      
      __m128 result;

      // this part needs ssse3
      __m128 result12 = _mm_hadd_ps(result1, result2);
      __m128 result34 = _mm_hadd_ps(result3, result4);
      result = _mm_hadd_ps(result12, result34);

      result = _mm_min_ps(result,val_max_128f); // mainly for 10-14 bit
	  result = _mm_max_ps(result,val_min_128f);
      // result = _mm_max_ps(result, zero); low limit through pack_us

      // Converts the four single-precision, floating-point values of a to signed 32-bit integer values.
      __m128i result_4x_int32  = _mm_cvtps_epi32(result);  // 4 * 32 bit integers
      // SIMD Extensions 4 (SSE4) packus or simulation
      __m128i result_4x_uint16 = hasSSE41 ? _mm_packus_epi32(result_4x_int32, zero) : (_MM_PACKUS_EPI32(result_4x_int32, zero)) ; // 4*32+zeros = lower 4*16 OK
      _mm_storel_epi64(reinterpret_cast<__m128i *>(dst + x), result_4x_uint16);
    }

    dst += dst_pitch;
    src += src_pitch;
  }
}
#endif


static void resizer_h_ssse3_generic(BYTE* dst, const BYTE* src, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height, int bits_per_pixel,const uint8_t range,const bool mode_YUY2)
{
  const int filter_size = AlignNumber(program->filter_size, 8) >> 3;
  const __m128i zero = _mm_setzero_si128();

	const int val_min = (range==1) ? 0 : 16;
	const int val_max = ((range==1) || (range==4)) ? 255 : (range==2) ? 235 : 240;

	const __m128i val_min_m128 = _mm_set1_epi16((short)((val_min << 8)|val_min));
	const __m128i val_max_m128 = (mode_YUY2 && ((range>=2) && (range<=3))) ? _mm_set1_epi16((short)(((int)240 << 8)|235)) : _mm_set1_epi16((short)((val_max << 8)|val_max));

  for (int y = 0; y < height; y++)
  {
    const short *current_coeff = program->pixel_coefficient;
	
    for (int x = 0; x < width; x+=4)
	{
      __m128i result1 = _mm_setr_epi32(8192, 0, 0, 0);
      __m128i result2 = _mm_setr_epi32(8192, 0, 0, 0);
      __m128i result3 = _mm_setr_epi32(8192, 0, 0, 0);
      __m128i result4 = _mm_setr_epi32(8192, 0, 0, 0);

      const int begin1 = program->pixel_offset[x+0];
      const int begin2 = program->pixel_offset[x+1];
      const int begin3 = program->pixel_offset[x+2];
      const int begin4 = program->pixel_offset[x+3];

	  for (int i = 0; i < filter_size; i++)
	  {
	    __m128i data, coeff, current_result;
		data = _mm_loadl_epi64(reinterpret_cast<const __m128i*>(src+begin1+i*8)); // 8 * 8 bit pixels
        data = _mm_unpacklo_epi8(data, zero);
	    coeff = _mm_load_si128(reinterpret_cast<const __m128i*>(current_coeff)); // 8 coeffs 14 bit scaled -> ushort OK
		current_result = _mm_madd_epi16(data, coeff);
        result1 = _mm_add_epi32(result1, current_result);
			
	    current_coeff += 8;		
	  }

      for (int i = 0; i < filter_size; i++)
	  {
		__m128i data, coeff, current_result;
        data = _mm_loadl_epi64(reinterpret_cast<const __m128i*>(src+begin2+i*8));
	    data = _mm_unpacklo_epi8(data, zero);
		coeff = _mm_load_si128(reinterpret_cast<const __m128i*>(current_coeff));
        current_result = _mm_madd_epi16(data, coeff);
	    result2 = _mm_add_epi32(result2, current_result);

		current_coeff += 8;
      }

      for (int i = 0; i < filter_size; i++)
	  {
		__m128i data, coeff, current_result;
        data = _mm_loadl_epi64(reinterpret_cast<const __m128i*>(src+begin3+i*8));
	    data = _mm_unpacklo_epi8(data, zero);
		coeff = _mm_load_si128(reinterpret_cast<const __m128i*>(current_coeff));
        current_result = _mm_madd_epi16(data, coeff);
	    result3 = _mm_add_epi32(result3, current_result);

		current_coeff += 8;
	  }

      for (int i = 0; i < filter_size; i++)
	  {
		__m128i data, coeff, current_result;
        data = _mm_loadl_epi64(reinterpret_cast<const __m128i*>(src+begin4+i*8));
		data = _mm_unpacklo_epi8(data, zero);
	    coeff = _mm_load_si128(reinterpret_cast<const __m128i*>(current_coeff));
        current_result = _mm_madd_epi16(data, coeff);
	    result4 = _mm_add_epi32(result4, current_result);

        current_coeff += 8;
	  }

      __m128i result12 = _mm_hadd_epi32(result1, result2);
      __m128i result34 = _mm_hadd_epi32(result3, result4);
      __m128i result = _mm_hadd_epi32(result12, result34);

      result = _mm_srai_epi32(result, 14);

      result = _mm_packs_epi32(result, zero);
      result = _mm_packus_epi16(result, zero);

		result = _mm_max_epu8(result,val_min_m128);
		result = _mm_min_epu8(result,val_max_m128);

      *((int*)(dst+x)) = _mm_cvtsi128_si32(result);
    }

    dst += dst_pitch;
    src += src_pitch;
  }
}


static void resizer_h_ssse3_8(BYTE* dst, const BYTE* src, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height, int bits_per_pixel,const uint8_t range,const bool mode_YUY2)
{
  const int filter_size = AlignNumber(program->filter_size, 8) >> 3;

  const __m128i zero = _mm_setzero_si128();

	const int val_min = (range==1) ? 0 : 16;
	const int val_max = ((range==1) || (range==4)) ? 255 : (range==2) ? 235 : 240;
	const int TabMax[4] = {235,240,235,240};

	const __m128i val_min_m128 = _mm_set1_epi16((short)((val_min << 8)|val_min));
	const __m128i val_max_m128 = (mode_YUY2 && ((range>=2) && (range<=3))) ? _mm_set1_epi16((short)(((int)240 << 8)|235)) : _mm_set1_epi16((short)((val_max << 8)|val_max));

  for (int y = 0; y < height; y++)
  {
    short *current_coeff = program->pixel_coefficient;
	
    for (int x = 0; x < width; x+=4)
	{
      __m128i result1 = _mm_setr_epi32(8192, 0, 0, 0);
      __m128i result2 = _mm_setr_epi32(8192, 0, 0, 0);
      __m128i result3 = _mm_setr_epi32(8192, 0, 0, 0);
      __m128i result4 = _mm_setr_epi32(8192, 0, 0, 0);

      const int begin1 = program->pixel_offset[x+0];
      const int begin2 = program->pixel_offset[x+1];
      const int begin3 = program->pixel_offset[x+2];
      const int begin4 = program->pixel_offset[x+3];

      __m128i data, coeff, current_result;

      // Unroll 1
      data = _mm_loadl_epi64(reinterpret_cast<const __m128i*>(src+begin1));
      data = _mm_unpacklo_epi8(data, zero);
      coeff = _mm_load_si128(reinterpret_cast<const __m128i*>(current_coeff));
      current_result = _mm_madd_epi16(data, coeff);
      result1 = _mm_add_epi32(result1, current_result);

      current_coeff += 8;

      // Unroll 2
      data = _mm_loadl_epi64(reinterpret_cast<const __m128i*>(src+begin2));
      data = _mm_unpacklo_epi8(data, zero);
      coeff = _mm_load_si128(reinterpret_cast<const __m128i*>(current_coeff));
      current_result = _mm_madd_epi16(data, coeff);
      result2 = _mm_add_epi32(result2, current_result);

      current_coeff += 8;

      // Unroll 3
      data = _mm_loadl_epi64(reinterpret_cast<const __m128i*>(src+begin3));
      data = _mm_unpacklo_epi8(data, zero);
      coeff = _mm_load_si128(reinterpret_cast<const __m128i*>(current_coeff));
      current_result = _mm_madd_epi16(data, coeff);
      result3 = _mm_add_epi32(result3, current_result);

      current_coeff += 8;

      // Unroll 4
      data = _mm_loadl_epi64(reinterpret_cast<const __m128i*>(src+begin4));
      data = _mm_unpacklo_epi8(data, zero);
      coeff = _mm_load_si128(reinterpret_cast<const __m128i*>(current_coeff));
      current_result = _mm_madd_epi16(data, coeff);
      result4 = _mm_add_epi32(result4, current_result);

      current_coeff += 8;

      // Combine
      __m128i result12 = _mm_hadd_epi32(result1, result2);
      __m128i result34 = _mm_hadd_epi32(result3, result4);
      __m128i result = _mm_hadd_epi32(result12, result34);

      result = _mm_srai_epi32(result, 14);

      result = _mm_packs_epi32(result, zero);
      result = _mm_packus_epi16(result, zero);

		result = _mm_max_epu8(result,val_min_m128);
		result = _mm_min_epu8(result,val_max_m128);

      *((int*)(dst+x)) = _mm_cvtsi128_si32(result);
    }

    dst += dst_pitch;
    src += src_pitch;
  }
}


#ifdef AVX_BUILD_POSSIBLE

static void resizer_h_avx_generic_int16_float_32(BYTE* dst8, const BYTE* src8, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height, int bits_per_pixel,const uint8_t range,const bool mode_YUY2)
{
  _mm256_zeroupper();
  
  const int filter_size = AlignNumber(program->filter_size, 8) >> 3;

  const float *src = reinterpret_cast<const float *>(src8);
  float *dst = reinterpret_cast<float *>(dst8);
  dst_pitch>>=2;
  src_pitch>>=2;
  
  __m128 data_l_single, data_h_single;

  for (int y = 0; y < height; y++)
  {
    const float *current_coeff = program->pixel_coefficient_float;

    for (int x = 0; x < width; x+=4)
	{
      __m256 result1 = _mm256_setzero_ps();
      __m256 result2 = result1;
      __m256 result3 = result1;
      __m256 result4 = result1;

      const int begin1 = program->pixel_offset[x+0];
      const int begin2 = program->pixel_offset[x+1];
      const int begin3 = program->pixel_offset[x+2];
      const int begin4 = program->pixel_offset[x+3];

      // this part is repeated by x4
      // begin1, result1
      for (int i = 0; i < filter_size; i++)
	  {
        __m256 data_single;
		
        // float unaligned
        // using one 256bit load instead of 2x128bit is slower on avx-only Ivy
        data_l_single = _mm_loadu_ps(reinterpret_cast<const float*>(src + begin1 + i * 8)); // float  4*32=128 4 pixels at a time
        data_h_single = _mm_loadu_ps(reinterpret_cast<const float*>(src + begin1 + i * 8 + 4)); // float  4*32=128 4 pixels at a time
        data_single = _mm256_set_m128(data_h_single, data_l_single);
		
        __m256 coeff;
		
        __m128 coeff_l = _mm_load_ps(reinterpret_cast<const float*>(current_coeff));    // always aligned
        __m128 coeff_h = _mm_load_ps(reinterpret_cast<const float*>(current_coeff + 4));  // always aligned
        coeff = _mm256_set_m128(coeff_h, coeff_l);

        __m256 dst = _mm256_mul_ps(data_single, coeff); // Multiply by coefficient
        result1 = _mm256_add_ps(result1, dst);

        current_coeff += 8;
      }

      // begin2, result2
      for (int i = 0; i < filter_size; i++)
	  {
        __m256 data_single;
		
        // float
        // using one 256bit load instead of 2x128bit is slower on avx-only Ivy
        data_l_single = _mm_loadu_ps(reinterpret_cast<const float*>(src + begin2 + i * 8)); // float  4*32=128 4 pixels at a time
        data_h_single = _mm_loadu_ps(reinterpret_cast<const float*>(src + begin2 + i * 8 + 4)); // float  4*32=128 4 pixels at a time
        data_single = _mm256_set_m128(data_h_single, data_l_single);
		
        __m256 coeff;
        __m128 coeff_l = _mm_load_ps(reinterpret_cast<const float*>(current_coeff));    // always aligned
        __m128 coeff_h = _mm_load_ps(reinterpret_cast<const float*>(current_coeff + 4));  // always aligned
        coeff = _mm256_set_m128(coeff_h, coeff_l);
		
        __m256 dst = _mm256_mul_ps(data_single, coeff); // Multiply by coefficient
        result2 = _mm256_add_ps(result2, dst);

        current_coeff += 8;
      }

      // begin3, result3
      for (int i = 0; i < filter_size; i++)
	  {
        __m256 data_single;
		
        // float
        // using one 256bit load instead of 2x128bit is slower on avx-only Ivy
        data_l_single = _mm_loadu_ps(reinterpret_cast<const float*>(src + begin3 + i * 8)); // float  4*32=128 4 pixels at a time
        data_h_single = _mm_loadu_ps(reinterpret_cast<const float*>(src + begin3 + i * 8 + 4)); // float  4*32=128 4 pixels at a time
        data_single = _mm256_set_m128(data_h_single, data_l_single);
		
        __m256 coeff;
        __m128 coeff_l = _mm_load_ps(reinterpret_cast<const float*>(current_coeff));    // always aligned
        __m128 coeff_h = _mm_load_ps(reinterpret_cast<const float*>(current_coeff + 4));  // always aligned
        coeff = _mm256_set_m128(coeff_h, coeff_l);
		
        __m256 dst = _mm256_mul_ps(data_single, coeff); // Multiply by coefficient
        result3 = _mm256_add_ps(result3, dst);

        current_coeff += 8;
      }

      // begin4, result4
      for (int i = 0; i < filter_size; i++)
	  {
        __m256 data_single;
		
        // float
        // using one 256bit load instead of 2x128bit is slower on avx-only Ivy
        data_l_single = _mm_loadu_ps(reinterpret_cast<const float*>(src + begin4 + i * 8)); // float  4*32=128 4 pixels at a time
        data_h_single = _mm_loadu_ps(reinterpret_cast<const float*>(src + begin4 + i * 8 + 4)); // float  4*32=128 4 pixels at a time
        data_single = _mm256_set_m128(data_h_single, data_l_single);
        
        __m256 coeff;
        __m128 coeff_l = _mm_load_ps(reinterpret_cast<const float*>(current_coeff));    // always aligned
        __m128 coeff_h = _mm_load_ps(reinterpret_cast<const float*>(current_coeff + 4));  // always aligned
        coeff = _mm256_set_m128(coeff_h, coeff_l);
		
        __m256 dst = _mm256_mul_ps(data_single, coeff); // Multiply by coefficient
        result4 = _mm256_add_ps(result4, dst);

        current_coeff += 8;
      }
      
      __m128 result;

      const __m128 sumQuad1 = _mm_add_ps(_mm256_castps256_ps128(result1), _mm256_extractf128_ps(result1, 1));
      const __m128 sumQuad2 = _mm_add_ps(_mm256_castps256_ps128(result2), _mm256_extractf128_ps(result2, 1));
      __m128 result12 = _mm_hadd_ps(sumQuad1, sumQuad2);
      const __m128 sumQuad3 = _mm_add_ps(_mm256_castps256_ps128(result3), _mm256_extractf128_ps(result3, 1));
      const __m128 sumQuad4 = _mm_add_ps(_mm256_castps256_ps128(result4), _mm256_extractf128_ps(result4, 1));
      __m128 result34 = _mm_hadd_ps(sumQuad3, sumQuad4);
      result = _mm_hadd_ps(result12, result34);

      // float
      // aligned
      _mm_store_ps(reinterpret_cast<float*>(dst+x), result); // 4 results at a time
      
    }

    dst += dst_pitch;
    src += src_pitch;
  }
  _mm256_zeroupper();
}


static void resizer_h_avx_generic_int16_float_16(BYTE* dst8, const BYTE* src8, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height, int bits_per_pixel,const uint8_t range,const bool mode_YUY2)
{
  _mm256_zeroupper();
  
  const int filter_size = AlignNumber(program->filter_size, 8) >> 3;
  const __m128i zero128 = _mm_setzero_si128();

  const uint16_t *src = reinterpret_cast<const uint16_t *>(src8);
  uint16_t *dst = reinterpret_cast<uint16_t *>(dst8);
  dst_pitch>>=1;
  src_pitch>>=1;
  
	const int val_min = (range==1) ? 0 : (int)16 << (bits_per_pixel-8);
	const int val_235 = (int)235 << (bits_per_pixel-8);
	const int val_240 = (int)240 << (bits_per_pixel-8);
	const int val_max = ((range==1) || (range==4)) ? ((int)1 << bits_per_pixel)-1 : (range==2) ? val_235 : val_240;
	const float val_minf = (float)val_min;
	const float val_maxf = (float)val_max;

  const __m128 val_min_128f = _mm_set1_ps(val_minf); // clamp limit
  const __m128 val_max_128f = _mm_set1_ps(val_maxf); // clamp limit

  for (int y = 0; y < height; y++)
  {
    const float *current_coeff = program->pixel_coefficient_float;

    for (int x = 0; x < width; x+=4)
	{
      __m256 result1 = _mm256_setzero_ps();
      __m256 result2 = result1;
      __m256 result3 = result1;
      __m256 result4 = result1;

      const int begin1 = program->pixel_offset[x+0];
      const int begin2 = program->pixel_offset[x+1];
      const int begin3 = program->pixel_offset[x+2];
      const int begin4 = program->pixel_offset[x+3];

      // this part is repeated by x4
      // begin1, result1
      for (int i = 0; i < filter_size; i++)
	  {
        __m256 data_single;
		
        // word
        __m256i src256;
        __m128i src_p = _mm_loadu_si128(reinterpret_cast<const __m128i*>(src + begin1 + i * 8)); // uint16_t  8*16=128 8 pixels at a time
        __m128i src_l = _mm_unpacklo_epi16(src_p, zero128); // spread lower  4*uint16_t pixel value -> 4*32 bit
        __m128i src_h = _mm_unpackhi_epi16(src_p, zero128); // spread higher 4*uint16_t pixel value -> 4*32 bit
        src256 = _mm256_set_m128i(src_h, src_l);
		
        data_single = _mm256_cvtepi32_ps(src256); // Converts the 8x signed 32-bit integer values of a to single-precision, floating-point values.
		
        __m256 coeff;
        __m128 coeff_l = _mm_load_ps(reinterpret_cast<const float*>(current_coeff));    // always aligned
        __m128 coeff_h = _mm_load_ps(reinterpret_cast<const float*>(current_coeff + 4));  // always aligned
        coeff = _mm256_set_m128(coeff_h, coeff_l);

        __m256 dst = _mm256_mul_ps(data_single, coeff); // Multiply by coefficient
        result1 = _mm256_add_ps(result1, dst);

        current_coeff += 8;
      }

      // begin2, result2
      for (int i = 0; i < filter_size; i++)
	  {
        __m256 data_single;
		
        // word
        __m256i src256;
        __m128i src_p = _mm_loadu_si128(reinterpret_cast<const __m128i*>(src + begin2 + i * 8)); // uint16_t  8*16=128 8 pixels at a time
        __m128i src_l = _mm_unpacklo_epi16(src_p, zero128); // spread lower  4*uint16_t pixel value -> 4*32 bit
        __m128i src_h = _mm_unpackhi_epi16(src_p, zero128); // spread higher 4*uint16_t pixel value -> 4*32 bit
        src256 = _mm256_set_m128i(src_h, src_l);
		
        data_single = _mm256_cvtepi32_ps (src256); // Converts the 8x signed 32-bit integer values of a to single-precision, floating-point values.
		
        __m256 coeff;
        __m128 coeff_l = _mm_load_ps(reinterpret_cast<const float*>(current_coeff));    // always aligned
        __m128 coeff_h = _mm_load_ps(reinterpret_cast<const float*>(current_coeff + 4));  // always aligned
        coeff = _mm256_set_m128(coeff_h, coeff_l);
		
        __m256 dst = _mm256_mul_ps(data_single, coeff); // Multiply by coefficient
        result2 = _mm256_add_ps(result2, dst);

        current_coeff += 8;
      }

      // begin3, result3
      for (int i = 0; i < filter_size; i++)
	  {
        __m256 data_single;
		
        // word
        __m256i src256;
        __m128i src_p = _mm_loadu_si128(reinterpret_cast<const __m128i*>(src + begin3 + i * 8)); // uint16_t  8*16=128 8 pixels at a time
        __m128i src_l = _mm_unpacklo_epi16(src_p, zero128); // spread lower  4*uint16_t pixel value -> 4*32 bit
        __m128i src_h = _mm_unpackhi_epi16(src_p, zero128); // spread higher 4*uint16_t pixel value -> 4*32 bit
        src256 = _mm256_set_m128i(src_h, src_l);
		
        data_single = _mm256_cvtepi32_ps (src256); // Converts the 8x signed 32-bit integer values of a to single-precision, floating-point values.
		
        __m256 coeff;
        __m128 coeff_l = _mm_load_ps(reinterpret_cast<const float*>(current_coeff));    // always aligned
        __m128 coeff_h = _mm_load_ps(reinterpret_cast<const float*>(current_coeff + 4));  // always aligned
        coeff = _mm256_set_m128(coeff_h, coeff_l);
		
        __m256 dst = _mm256_mul_ps(data_single, coeff); // Multiply by coefficient
        result3 = _mm256_add_ps(result3, dst);

        current_coeff += 8;
      }

      // begin4, result4
      for (int i = 0; i < filter_size; i++)
	  {
        __m256 data_single;
		
        // word
        __m256i src256;
        __m128i src_p = _mm_loadu_si128(reinterpret_cast<const __m128i*>(src + begin4 + i * 8)); // uint16_t  8*16=128 8 pixels at a time
        __m128i src_l = _mm_unpacklo_epi16(src_p, zero128); // spread lower  4*uint16_t pixel value -> 4*32 bit
        __m128i src_h = _mm_unpackhi_epi16(src_p, zero128); // spread higher 4*uint16_t pixel value -> 4*32 bit
        src256 = _mm256_set_m128i(src_h, src_l);
		
        data_single = _mm256_cvtepi32_ps (src256); // Converts the 8x signed 32-bit integer values of a to single-precision, floating-point values.
		
        __m256 coeff;
        __m128 coeff_l = _mm_load_ps(reinterpret_cast<const float*>(current_coeff));    // always aligned
        __m128 coeff_h = _mm_load_ps(reinterpret_cast<const float*>(current_coeff + 4));  // always aligned
        coeff = _mm256_set_m128(coeff_h, coeff_l);
		
        __m256 dst = _mm256_mul_ps(data_single, coeff); // Multiply by coefficient
        result4 = _mm256_add_ps(result4, dst);

        current_coeff += 8;
      }
      
      __m128 result;

      const __m128 sumQuad1 = _mm_add_ps(_mm256_castps256_ps128(result1), _mm256_extractf128_ps(result1, 1));
      const __m128 sumQuad2 = _mm_add_ps(_mm256_castps256_ps128(result2), _mm256_extractf128_ps(result2, 1));
      __m128 result12 = _mm_hadd_ps(sumQuad1, sumQuad2);
      const __m128 sumQuad3 = _mm_add_ps(_mm256_castps256_ps128(result3), _mm256_extractf128_ps(result3, 1));
      const __m128 sumQuad4 = _mm_add_ps(_mm256_castps256_ps128(result4), _mm256_extractf128_ps(result4, 1));
      __m128 result34 = _mm_hadd_ps(sumQuad3, sumQuad4);
      result = _mm_hadd_ps(result12, result34);

      // clamp!
      result = _mm_min_ps(result,val_max_128f); // for 10-14 bit 
	  result = _mm_max_ps(result,val_min_128f);
      // low limit or 16 bit limit through pack_us
      // Converts the four single-precision, floating-point values of a to signed 32-bit integer values.
      __m128i result_4x_int32  = _mm_cvtps_epi32(result);  // 4 * 32 bit integers
      __m128i result_4x_uint16 = _mm_packus_epi32(result_4x_int32, zero128); // 4*32+zeros = lower 4*16 OK
      _mm_storel_epi64(reinterpret_cast<__m128i *>(dst + x), result_4x_uint16);
	  
	}

    dst += dst_pitch;
    src += src_pitch;
  }
  _mm256_zeroupper();
}


static void resizer_h_avx2_generic_int16_float_32(BYTE* dst8, const BYTE* src8, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height, int bits_per_pixel,const uint8_t range,const bool mode_YUY2)
{
  _mm256_zeroupper();
  
  const int filter_size = AlignNumber(program->filter_size, 8) >> 3;

  const float *src = reinterpret_cast<const float *>(src8);
  float *dst = reinterpret_cast<float *>(dst8);
  dst_pitch>>=2;
  src_pitch>>=2;
  
  for (int y = 0; y < height; y++)
  {
    const float *current_coeff = program->pixel_coefficient_float;

    for (int x = 0; x < width; x+=4)
	{
      __m256 result1 = _mm256_setzero_ps();
      __m256 result2 = result1;
      __m256 result3 = result1;
      __m256 result4 = result1;

      const int begin1 = program->pixel_offset[x+0];
      const int begin2 = program->pixel_offset[x+1];
      const int begin3 = program->pixel_offset[x+2];
      const int begin4 = program->pixel_offset[x+3];

      // this part is repeated by x4
      // begin1, result1
      for (int i = 0; i < filter_size; i++)
	  {
        __m256 data_single;
		
        // float unaligned
        data_single = _mm256_loadu_ps(reinterpret_cast<const float*>(src+begin1+i*8)); // float  8*32=256 8 pixels at a time
		
        __m256 coeff;
        // using one 256bit load instead of 2x128bit is slower on avx-only Ivy
        coeff = _mm256_load_ps(reinterpret_cast<const float*>(current_coeff));    // always aligned

        __m256 dst = _mm256_mul_ps(data_single, coeff); // Multiply by coefficient
        result1 = _mm256_add_ps(result1, dst);

        current_coeff += 8;
      }

      // begin2, result2
      for (int i = 0; i < filter_size; i++)
	  {
        __m256 data_single;
		
        // float
        data_single = _mm256_loadu_ps(reinterpret_cast<const float*>(src+begin2+i*8)); // float  8*32=256 8 pixels at a time
		
        __m256 coeff;
        // using one 256bit load instead of 2x128bit is slower on avx-only Ivy
        coeff = _mm256_load_ps(reinterpret_cast<const float*>(current_coeff));    // always aligned
		
        __m256 dst = _mm256_mul_ps(data_single, coeff); // Multiply by coefficient
        result2 = _mm256_add_ps(result2, dst);

        current_coeff += 8;
      }

      // begin3, result3
      for (int i = 0; i < filter_size; i++)
	  {
        __m256 data_single;
		
        // float
        data_single = _mm256_loadu_ps(reinterpret_cast<const float*>(src+begin3+i*8)); // float  8*32=256 8 pixels at a time
		
        __m256 coeff;
        // using one 256bit load instead of 2x128bit is slower on avx-only Ivy
        coeff = _mm256_load_ps(reinterpret_cast<const float*>(current_coeff));    // always aligned
		
        __m256 dst = _mm256_mul_ps(data_single, coeff); // Multiply by coefficient
        result3 = _mm256_add_ps(result3, dst);

        current_coeff += 8;
      }

      // begin4, result4
      for (int i = 0; i < filter_size; i++)
	  {
        __m256 data_single;
		
        // float
        data_single = _mm256_loadu_ps(reinterpret_cast<const float*>(src+begin4+i*8)); // float  8*32=256 8 pixels at a time
		
        __m256 coeff;
        // using one 256bit load instead of 2x128bit is slower on avx-only Ivy
        coeff = _mm256_load_ps(reinterpret_cast<const float*>(current_coeff));    // always aligned
		
        __m256 dst = _mm256_mul_ps(data_single, coeff); // Multiply by coefficient
        result4 = _mm256_add_ps(result4, dst);

        current_coeff += 8;
      }
      
      __m128 result;

      const __m128 sumQuad1 = _mm_add_ps(_mm256_castps256_ps128(result1), _mm256_extractf128_ps(result1, 1));
      const __m128 sumQuad2 = _mm_add_ps(_mm256_castps256_ps128(result2), _mm256_extractf128_ps(result2, 1));
      __m128 result12 = _mm_hadd_ps(sumQuad1, sumQuad2);
      const __m128 sumQuad3 = _mm_add_ps(_mm256_castps256_ps128(result3), _mm256_extractf128_ps(result3, 1));
      const __m128 sumQuad4 = _mm_add_ps(_mm256_castps256_ps128(result4), _mm256_extractf128_ps(result4, 1));
      __m128 result34 = _mm_hadd_ps(sumQuad3, sumQuad4);
      result = _mm_hadd_ps(result12, result34);

      // float
      // aligned
      _mm_store_ps(reinterpret_cast<float*>(dst+x), result); // 4 results at a time      
    }

    dst += dst_pitch;
    src += src_pitch;
  }
  _mm256_zeroupper();
}


static void resizer_h_avx2_generic_int16_float_16(BYTE* dst8, const BYTE* src8, int dst_pitch, int src_pitch, ResamplingProgram* program, int width, int height, int bits_per_pixel,const uint8_t range,const bool mode_YUY2)
{
  _mm256_zeroupper();
  
  const int filter_size = AlignNumber(program->filter_size, 8) >> 3;
  const __m128i zero128 = _mm_setzero_si128();

  const uint16_t *src = reinterpret_cast<const uint16_t *>(src8);
  uint16_t *dst = reinterpret_cast<uint16_t *>(dst8);
  dst_pitch>>=1;
  src_pitch>>=1;
  
	const int val_min = (range==1) ? 0 : (int)16 << (bits_per_pixel-8);
	const int val_235 = (int)235 << (bits_per_pixel-8);
	const int val_240 = (int)240 << (bits_per_pixel-8);
	const int val_max = ((range==1) || (range==4)) ? ((int)1 << bits_per_pixel)-1 : (range==2) ? val_235 : val_240;
	const float val_minf = (float)val_min;
	const float val_maxf = (float)val_max;

  const __m128 val_min_128f = _mm_set1_ps(val_minf); // clamp limit
  const __m128 val_max_128f = _mm_set1_ps(val_maxf); // clamp limit

  for (int y = 0; y < height; y++)
  {
    const float *current_coeff = program->pixel_coefficient_float;

    for (int x = 0; x < width; x+=4)
	{
      __m256 result1 = _mm256_setzero_ps();
      __m256 result2 = result1;
      __m256 result3 = result1;
      __m256 result4 = result1;

      const int begin1 = program->pixel_offset[x+0];
      const int begin2 = program->pixel_offset[x+1];
      const int begin3 = program->pixel_offset[x+2];
      const int begin4 = program->pixel_offset[x+3];

      // this part is repeated by x4
      // begin1, result1
      for (int i = 0; i < filter_size; i++)
	  {
        __m256 data_single;
		
        // word
        __m256i src256;
        src256 = _mm256_cvtepu16_epi32(_mm_loadu_si128(reinterpret_cast<const __m128i*>(src + begin1 + i * 8))); // 8*16->8*32 bits
		  
        data_single = _mm256_cvtepi32_ps(src256); // Converts the 8x signed 32-bit integer values of a to single-precision, floating-point values.
		
        __m256 coeff;
        // using one 256bit load instead of 2x128bit is slower on avx-only Ivy
        coeff = _mm256_load_ps(reinterpret_cast<const float*>(current_coeff));    // always aligned

        __m256 dst = _mm256_mul_ps(data_single, coeff); // Multiply by coefficient
        result1 = _mm256_add_ps(result1, dst);

        current_coeff += 8;
      }

      // begin2, result2
      for (int i = 0; i < filter_size; i++)
	  {
        __m256 data_single;
		
        // word
        __m256i src256;
        src256 = _mm256_cvtepu16_epi32(_mm_loadu_si128(reinterpret_cast<const __m128i*>(src + begin2 + i * 8))); // 8*16->8*32 bits
		
        data_single = _mm256_cvtepi32_ps (src256); // Converts the 8x signed 32-bit integer values of a to single-precision, floating-point values.
		
        __m256 coeff;
        // using one 256bit load instead of 2x128bit is slower on avx-only Ivy
        coeff = _mm256_load_ps(reinterpret_cast<const float*>(current_coeff));    // always aligned
		
        __m256 dst = _mm256_mul_ps(data_single, coeff); // Multiply by coefficient
        result2 = _mm256_add_ps(result2, dst);

        current_coeff += 8;
      }

      // begin3, result3
      for (int i = 0; i < filter_size; i++)
	  {
        __m256 data_single;
		
        // word
        __m256i src256;
        src256 = _mm256_cvtepu16_epi32(_mm_loadu_si128(reinterpret_cast<const __m128i*>(src + begin3 + i * 8))); // 8*16->8*32 bits
		
        data_single = _mm256_cvtepi32_ps (src256); // Converts the 8x signed 32-bit integer values of a to single-precision, floating-point values.
		
        __m256 coeff;
        // using one 256bit load instead of 2x128bit is slower on avx-only Ivy
        coeff = _mm256_load_ps(reinterpret_cast<const float*>(current_coeff));    // always aligned
		
        __m256 dst = _mm256_mul_ps(data_single, coeff); // Multiply by coefficient
        result3 = _mm256_add_ps(result3, dst);

        current_coeff += 8;
      }

      // begin4, result4
      for (int i = 0; i < filter_size; i++)
	  {
        __m256 data_single;
		
        // word
        __m256i src256;
        src256 = _mm256_cvtepu16_epi32(_mm_loadu_si128(reinterpret_cast<const __m128i*>(src + begin4 + i * 8))); // 8*16->8*32 bits
		
        data_single = _mm256_cvtepi32_ps (src256); // Converts the 8x signed 32-bit integer values of a to single-precision, floating-point values.
		
        __m256 coeff;
        // using one 256bit load instead of 2x128bit is slower on avx-only Ivy
        coeff = _mm256_load_ps(reinterpret_cast<const float*>(current_coeff));    // always aligned
		
        __m256 dst = _mm256_mul_ps(data_single, coeff); // Multiply by coefficient
        result4 = _mm256_add_ps(result4, dst);

        current_coeff += 8;
      }
      
      __m128 result;

      const __m128 sumQuad1 = _mm_add_ps(_mm256_castps256_ps128(result1), _mm256_extractf128_ps(result1, 1));
      const __m128 sumQuad2 = _mm_add_ps(_mm256_castps256_ps128(result2), _mm256_extractf128_ps(result2, 1));
      __m128 result12 = _mm_hadd_ps(sumQuad1, sumQuad2);
      const __m128 sumQuad3 = _mm_add_ps(_mm256_castps256_ps128(result3), _mm256_extractf128_ps(result3, 1));
      const __m128 sumQuad4 = _mm_add_ps(_mm256_castps256_ps128(result4), _mm256_extractf128_ps(result4, 1));
      __m128 result34 = _mm_hadd_ps(sumQuad3, sumQuad4);
      result = _mm_hadd_ps(result12, result34);

      // clamp!
      result = _mm_min_ps(result,val_max_128f); // for 10-14 bit 
	  result = _mm_max_ps(result,val_min_128f);
      // low limit or 16 bit limit through pack_us
      // Converts the four single-precision, floating-point values of a to signed 32-bit integer values.
      __m128i result_4x_int32  = _mm_cvtps_epi32(result);  // 4 * 32 bit integers
      __m128i result_4x_uint16 = _mm_packus_epi32(result_4x_int32, zero128); // 4*32+zeros = lower 4*16 OK
      _mm_storel_epi64(reinterpret_cast<__m128i *>(dst + x), result_4x_uint16);
    }

    dst += dst_pitch;
    src += src_pitch;
  }
  _mm256_zeroupper();
}

#endif



/********************************************************************
***** Declare index of new filters for Avisynth's filter engine *****
********************************************************************/



FilteredResizeH::FilteredResizeH( PClip _child, double subrange_left, double subrange_width,int target_width,
	uint8_t _threads,bool _sleep,int range_mode,bool desample,int accuracy,bool _avsp, ResamplingFunction* func, IScriptEnvironment* env )
  : GenericVideoFilter(_child),
  resampling_program_luma(NULL), resampling_program_chroma(NULL),
  filter_storage_luma(NULL), filter_storage_chroma(NULL),threads(_threads),sleep(_sleep),
  avsp(_avsp)
{
  src_width  = vi.width;
  src_height = vi.height;
  dst_width  = target_width;
  dst_height = vi.height;
  
  pixelsize = (uint8_t)vi.ComponentSize(); // AVS16
  grey = vi.IsY();
  bits_per_pixel = (uint8_t)vi.BitsPerComponent();
  isRGBPfamily = vi.IsPlanarRGB() || vi.IsPlanarRGBA();
  isAlphaChannel = vi.IsYUVA() || vi.IsPlanarRGBA();
  mode_YUY2 = vi.IsYUY2();

	if ((range_mode!=1) && (range_mode!=4))
	{
		if (vi.IsYUV())
		{
			plane_range[0]=2;
			plane_range[1]=3;
			plane_range[2]=3;
		}
		else
		{
			if (grey)
			{
				for (uint8_t i=0; i<3; i++)
					plane_range[i]=(range_mode==0) ? 2 : range_mode;
			}
			else
			{
				for (uint8_t i=0; i<3; i++)
					plane_range[i]=1;
			}
		}
	}
	else
	{
		if (vi.IsRGB()) range_mode=1;

		for (uint8_t i=0; i<3; i++)
			plane_range[i]=range_mode;
	}
	plane_range[3]=1;

	int16_t i;

	ResampleH_MT=StaticThreadpoolH;

	for (i=0; i<MAX_MT_THREADS; i++)
	{
		MT_Thread[i].pClass=this;
		MT_Thread[i].f_process=0;
		MT_Thread[i].thread_Id=(uint8_t)i;
		MT_Thread[i].pFunc=ResampleH_MT;
	}

	UserId=0;
	ghMutex=NULL;
	
	const int shift_w = (!grey && vi.IsPlanar() && !isRGBPfamily) ? vi.GetPlaneWidthSubsampling(PLANAR_U) : 0;
	const int shift_h = (!grey && vi.IsPlanar() && !isRGBPfamily) ? vi.GetPlaneHeightSubsampling(PLANAR_U) : 0;

	const int src_width = vi.IsPlanar() ? vi.width : vi.BytesFromPixels(vi.width)/pixelsize;
	const int dst_width = vi.IsPlanar() ? target_width : vi.BytesFromPixels(target_width)/pixelsize;

	if (vi.height<32) threads_number=1;
	else threads_number=threads;
	
	ghMutex=CreateMutex(NULL,FALSE,NULL);
	if (ghMutex==NULL) env->ThrowError("ResizeHMT: Unable to create Mutex!");

  // Main resampling program
  int SizeH;

  if (desample) resampling_program_luma = func->GetDesamplingProgram(target_width, subrange_left, subrange_width, vi.width, accuracy, 0, shift_w, SizeH, env);
  else
  {
	  resampling_program_luma = func->GetResamplingProgram(vi.width, subrange_left, subrange_width, target_width, env);
	  SizeH=dst_width;
  }

  threads_number=CreateMTData(threads_number,src_width,vi.height,SizeH,vi.height,shift_w,shift_h);

  if (vi.IsPlanar() && !grey && !isRGBPfamily)
  {
    const int div   = 1 << shift_w;

	if (desample)
	{
		int SizeOut;

	    resampling_program_chroma = func->GetDesamplingProgram(
		  target_width   >> shift_w,
	      subrange_left   / div,
		  subrange_width  / div,
	      vi.width   >> shift_w,
		  accuracy,SizeH,shift_w,SizeOut,
		  env);
	}
	else
	{
	    resampling_program_chroma = func->GetResamplingProgram(
		  vi.width       >> shift_w,
	      subrange_left   / div,
		  subrange_width  / div,
	      target_width   >> shift_w,
		  env);
	}
  }
  
  // Plannar + SSSE3 = use new horizontal resizer routines
  resampler_h_luma = GetResampler(env->GetCPUFlags(), true, pixelsize, bits_per_pixel, resampling_program_luma,env);

  if (vi.IsPlanar() && !grey && !isRGBPfamily) resampler_h_chroma = GetResampler(env->GetCPUFlags(), true, pixelsize, bits_per_pixel, resampling_program_chroma,env);

	if (threads_number>1)
	{
		if (!poolInterface->GetUserId(UserId))
		{
			poolInterface->DeAllocateAllThreads(true);
			FreeData();
			env->ThrowError("ResizeHMT: Error with the TheadPool while getting UserId!");
		}
	}
  
  // Change target video info size
  vi.width =SizeH;
}



int __stdcall FilteredResizeH::SetCacheHints(int cachehints,int frame_range)
{
  switch (cachehints)
  {
  case CACHE_GET_MTMODE :
    return MT_MULTI_INSTANCE;
  default :
    return 0;
  }
}


uint8_t FilteredResizeH::CreateMTData(uint8_t max_threads,int32_t src_size_x,int32_t src_size_y,int32_t dst_size_x,int32_t dst_size_y,int UV_w,int UV_h)
{
	int32_t _y_min,_dh;

	if ((max_threads<=1) || (max_threads>threads_number))
	{
		MT_Data[0].top=true;
		MT_Data[0].bottom=true;
		MT_Data[0].src_Y_h_min=0;
		MT_Data[0].dst_Y_h_min=0;
		MT_Data[0].src_Y_h_max=src_size_y;
		MT_Data[0].dst_Y_h_max=dst_size_y;
		MT_Data[0].src_UV_h_min=0;
		MT_Data[0].dst_UV_h_min=0;
		if (UV_h>0)
		{
			MT_Data[0].src_UV_h_max=src_size_y >> UV_h;
			MT_Data[0].dst_UV_h_max=dst_size_y >> UV_h;
		}
		else
		{
			MT_Data[0].src_UV_h_max=src_size_y;
			MT_Data[0].dst_UV_h_max=dst_size_y;
		}
		MT_Data[0].src_Y_w=src_size_x;
		MT_Data[0].dst_Y_w=dst_size_x;
		if (UV_w>0)
		{
			MT_Data[0].src_UV_w=src_size_x >> UV_w;
			MT_Data[0].dst_UV_w=dst_size_x >> UV_w;
		}
		else
		{
			MT_Data[0].src_UV_w=src_size_x;
			MT_Data[0].dst_UV_w=dst_size_x;
		}
		return(1);
	}

	int32_t src_dh_Y,src_dh_UV,dst_dh_Y,dst_dh_UV;
	int32_t h_y;
	uint8_t i,max=0;

	dst_dh_Y=(dst_size_y+(uint32_t)max_threads-1)/(uint32_t)max_threads;
	if (dst_dh_Y<16) dst_dh_Y=16;
	if ((dst_dh_Y & 3)!=0) dst_dh_Y=((dst_dh_Y+3) >> 2) << 2;

	if (src_size_y==dst_size_y) src_dh_Y=dst_dh_Y;
	else
	{
		src_dh_Y=(src_size_y+(uint32_t)max_threads-1)/(uint32_t)max_threads;
		if (src_dh_Y<16) src_dh_Y=16;
		if ((src_dh_Y & 3)!=0) src_dh_Y=((src_dh_Y+3) >> 2) << 2;
	}

	if (src_size_y<dst_size_y)
	{
		_y_min=src_size_y;
		_dh=src_dh_Y;
	}
	else
	{
		_y_min=dst_size_y;
		_dh=dst_dh_Y;
	}
	h_y=0;
	while (h_y<(_y_min-16))
	{
		max++;
		h_y+=_dh;
	}

	if (max==1)
	{
		MT_Data[0].top=true;
		MT_Data[0].bottom=true;
		MT_Data[0].src_Y_h_min=0;
		MT_Data[0].dst_Y_h_min=0;
		MT_Data[0].src_Y_h_max=src_size_y;
		MT_Data[0].dst_Y_h_max=dst_size_y;
		MT_Data[0].src_UV_h_min=0;
		MT_Data[0].dst_UV_h_min=0;
		if (UV_h>0)
		{
			MT_Data[0].src_UV_h_max=src_size_y >> UV_h;
			MT_Data[0].dst_UV_h_max=dst_size_y >> UV_h;
		}
		else
		{
			MT_Data[0].src_UV_h_max=src_size_y;
			MT_Data[0].dst_UV_h_max=dst_size_y;
		}
		MT_Data[0].src_Y_w=src_size_x;
		MT_Data[0].dst_Y_w=dst_size_x;
		if (UV_w>0)
		{
			MT_Data[0].src_UV_w=src_size_x >> UV_w;
			MT_Data[0].dst_UV_w=dst_size_x >> UV_w;
		}
		else
		{
			MT_Data[0].src_UV_w=src_size_x;
			MT_Data[0].dst_UV_w=dst_size_x;
		}
		return(1);
	}

	src_dh_UV= (UV_h>0) ? src_dh_Y>>UV_h : src_dh_Y;
	dst_dh_UV= (UV_h>0) ? dst_dh_Y>>UV_h : dst_dh_Y;

	MT_Data[0].top=true;
	MT_Data[0].bottom=false;
	MT_Data[0].src_Y_h_min=0;
	MT_Data[0].src_Y_h_max=src_dh_Y;
	MT_Data[0].dst_Y_h_min=0;
	MT_Data[0].dst_Y_h_max=dst_dh_Y;
	MT_Data[0].src_UV_h_min=0;
	MT_Data[0].src_UV_h_max=src_dh_UV;
	MT_Data[0].dst_UV_h_min=0;
	MT_Data[0].dst_UV_h_max=dst_dh_UV;

	i=1;
	while (i<max)
	{
		MT_Data[i].top=false;
		MT_Data[i].bottom=false;
		MT_Data[i].src_Y_h_min=MT_Data[i-1].src_Y_h_max;
		MT_Data[i].src_Y_h_max=MT_Data[i].src_Y_h_min+src_dh_Y;
		MT_Data[i].dst_Y_h_min=MT_Data[i-1].dst_Y_h_max;
		MT_Data[i].dst_Y_h_max=MT_Data[i].dst_Y_h_min+dst_dh_Y;
		MT_Data[i].src_UV_h_min=MT_Data[i-1].src_UV_h_max;
		MT_Data[i].src_UV_h_max=MT_Data[i].src_UV_h_min+src_dh_UV;
		MT_Data[i].dst_UV_h_min=MT_Data[i-1].dst_UV_h_max;
		MT_Data[i].dst_UV_h_max=MT_Data[i].dst_UV_h_min+dst_dh_UV;
		i++;
	}

	MT_Data[max-1].bottom=true;
	MT_Data[max-1].src_Y_h_max=src_size_y;
	MT_Data[max-1].dst_Y_h_max=dst_size_y;
	if (UV_h>0)
	{
		MT_Data[max-1].src_UV_h_max=src_size_y >> UV_h;
		MT_Data[max-1].dst_UV_h_max=dst_size_y >> UV_h;
	}
	else
	{
		MT_Data[max-1].src_UV_h_max=src_size_y;
		MT_Data[max-1].dst_UV_h_max=dst_size_y;
	}

	for (i=0; i<max; i++)
	{
		MT_Data[i].src_Y_w=src_size_x;
		MT_Data[i].dst_Y_w=dst_size_x;
		if (UV_w>0)
		{
			MT_Data[i].src_UV_w=src_size_x >> UV_w;
			MT_Data[i].dst_UV_w=dst_size_x >> UV_w;
		}
		else
		{
			MT_Data[i].src_UV_w=src_size_x;
			MT_Data[i].dst_UV_w=dst_size_x;
		}
	}

	return(max);
}


void FilteredResizeH::ResamplerLumaMT(uint8_t thread_num)
{
	const MT_Data_Info_ResampleMT mt_data_inf=MT_Data[thread_num];

	resampler_h_luma(mt_data_inf.dst1,mt_data_inf.src1,mt_data_inf.dst_pitch1,mt_data_inf.src_pitch1,
		mt_data_inf.resampling_program_luma,mt_data_inf.dst_Y_w,mt_data_inf.dst_Y_h_max-mt_data_inf.dst_Y_h_min,
		bits_per_pixel,plane_range[0],mode_YUY2);
}


void FilteredResizeH::ResamplerLumaMT2(uint8_t thread_num)
{
	const MT_Data_Info_ResampleMT mt_data_inf=MT_Data[thread_num];

	resampler_h_luma(mt_data_inf.dst2,mt_data_inf.src2,mt_data_inf.dst_pitch2,mt_data_inf.src_pitch2,
		mt_data_inf.resampling_program_luma,mt_data_inf.dst_Y_w,mt_data_inf.dst_Y_h_max-mt_data_inf.dst_Y_h_min,
		bits_per_pixel,plane_range[1],mode_YUY2);
}


void FilteredResizeH::ResamplerLumaMT3(uint8_t thread_num)
{
	const MT_Data_Info_ResampleMT mt_data_inf=MT_Data[thread_num];

	resampler_h_luma(mt_data_inf.dst3,mt_data_inf.src3,mt_data_inf.dst_pitch3,mt_data_inf.src_pitch3,
		mt_data_inf.resampling_program_luma,mt_data_inf.dst_Y_w,mt_data_inf.dst_Y_h_max-mt_data_inf.dst_Y_h_min,
		bits_per_pixel,plane_range[2],mode_YUY2);
}

void FilteredResizeH::ResamplerLumaMT4(uint8_t thread_num)
{
	const MT_Data_Info_ResampleMT mt_data_inf=MT_Data[thread_num];

	resampler_h_luma(mt_data_inf.dst4,mt_data_inf.src4,mt_data_inf.dst_pitch4,mt_data_inf.src_pitch4,
		mt_data_inf.resampling_program_luma,mt_data_inf.dst_Y_w,mt_data_inf.dst_Y_h_max-mt_data_inf.dst_Y_h_min,
		bits_per_pixel,plane_range[3],mode_YUY2);
}

void FilteredResizeH::ResamplerUChromaMT(uint8_t thread_num)
{
	const MT_Data_Info_ResampleMT mt_data_inf=MT_Data[thread_num];

	resampler_h_chroma(mt_data_inf.dst2,mt_data_inf.src2,mt_data_inf.dst_pitch2,mt_data_inf.src_pitch2,
		mt_data_inf.resampling_program_chroma,mt_data_inf.dst_UV_w,mt_data_inf.dst_UV_h_max-mt_data_inf.dst_UV_h_min,
		bits_per_pixel,plane_range[1],mode_YUY2);
}


void FilteredResizeH::ResamplerVChromaMT(uint8_t thread_num)
{
	const MT_Data_Info_ResampleMT mt_data_inf=MT_Data[thread_num];

	resampler_h_chroma(mt_data_inf.dst3,mt_data_inf.src3,mt_data_inf.dst_pitch3,mt_data_inf.src_pitch3,
		mt_data_inf.resampling_program_chroma,mt_data_inf.dst_UV_w,mt_data_inf.dst_UV_h_max-mt_data_inf.dst_UV_h_min,
		bits_per_pixel,plane_range[2],mode_YUY2);
}


void FilteredResizeH::StaticThreadpoolH(void *ptr)
{
	const Public_MT_Data_Thread *data=(const Public_MT_Data_Thread *)ptr;
	FilteredResizeH *ptrClass=(FilteredResizeH *)data->pClass;

	switch(data->f_process)
	{
		case 1 : ptrClass->ResamplerLumaMT(data->thread_Id);
			break;
		case 2 : ptrClass->ResamplerUChromaMT(data->thread_Id);
			break;
		case 3 : ptrClass->ResamplerVChromaMT(data->thread_Id);
			break;
		case 4 : ptrClass->ResamplerLumaMT2(data->thread_Id);
			break;
		case 5 : ptrClass->ResamplerLumaMT3(data->thread_Id);
			break;
		case 6 : ptrClass->ResamplerLumaMT4(data->thread_Id);
			break;		
		default : ;
	}
}


PVideoFrame __stdcall FilteredResizeH::GetFrame(int n, IScriptEnvironment* env)
{
  PVideoFrame src = child->GetFrame(n, env);
  PVideoFrame dst = env->NewVideoFrame(vi);
  
  const int src_pitch_1 = src->GetPitch();
  const int dst_pitch_1 = dst->GetPitch();
  const BYTE *srcp_1 = src->GetReadPtr();
        BYTE *dstp_1 = dst->GetWritePtr();

	const int src_pitch_2 = (!grey && vi.IsPlanar() && !isRGBPfamily) ? src->GetPitch(PLANAR_U) : (isRGBPfamily) ? src->GetPitch(PLANAR_B) : 0;
	const int dst_pitch_2 = (!grey && vi.IsPlanar() && !isRGBPfamily) ? dst->GetPitch(PLANAR_U) : (isRGBPfamily) ? dst->GetPitch(PLANAR_B) : 0;
	const BYTE *srcp_2 = (!grey && vi.IsPlanar() && !isRGBPfamily) ? src->GetReadPtr(PLANAR_U) : (isRGBPfamily) ? src->GetReadPtr(PLANAR_B) : NULL;
	BYTE *dstp_2 = (!grey && vi.IsPlanar() && !isRGBPfamily) ? dst->GetWritePtr(PLANAR_U) : (isRGBPfamily) ? dst->GetWritePtr(PLANAR_B) : NULL;

	const int src_pitch_3 = (!grey && vi.IsPlanar() && !isRGBPfamily) ? src->GetPitch(PLANAR_V) : (isRGBPfamily) ? src->GetPitch(PLANAR_R) : 0;
	const int dst_pitch_3 = (!grey && vi.IsPlanar() && !isRGBPfamily) ? dst->GetPitch(PLANAR_V) : (isRGBPfamily) ? dst->GetPitch(PLANAR_R) : 0;
	const BYTE *srcp_3 = (!grey && vi.IsPlanar() && !isRGBPfamily) ? src->GetReadPtr(PLANAR_V) : (isRGBPfamily) ? src->GetReadPtr(PLANAR_R) : NULL;
	BYTE *dstp_3 = (!grey && vi.IsPlanar() && !isRGBPfamily) ? dst->GetWritePtr(PLANAR_V) : (isRGBPfamily) ? dst->GetWritePtr(PLANAR_R) : NULL;
	
	const int src_pitch_4 = (isAlphaChannel) ? src->GetPitch(PLANAR_A) : 0;
	const int dst_pitch_4 = (isAlphaChannel) ? dst->GetPitch(PLANAR_A) : 0;
	const BYTE *srcp_4 = (isAlphaChannel) ? src->GetReadPtr(PLANAR_A) : NULL;
	BYTE *dstp_4 = (isAlphaChannel) ? dst->GetWritePtr(PLANAR_A) : NULL;
	
  WaitForSingleObject(ghMutex,INFINITE);
  
  if (threads_number>1)
  {
	if (!poolInterface->RequestThreadPool(UserId,threads_number,MT_Thread,-1,false))
	{
		ReleaseMutex(ghMutex);
		env->ThrowError("ResizeHMT: Error with the TheadPool while requesting threadpool!");
	}
  }
  
	for(uint8_t i=0; i<threads_number; i++)
	{
		MT_Data[i].src1=srcp_1+(MT_Data[i].src_Y_h_min*src_pitch_1);
		MT_Data[i].src2=srcp_2+(MT_Data[i].src_UV_h_min*src_pitch_2);
		MT_Data[i].src3=srcp_3+(MT_Data[i].src_UV_h_min*src_pitch_3);
		MT_Data[i].src4=srcp_4+(MT_Data[i].src_Y_h_min*src_pitch_4);
		MT_Data[i].src_pitch1=src_pitch_1;
		MT_Data[i].src_pitch2=src_pitch_2;
		MT_Data[i].src_pitch3=src_pitch_3;
		MT_Data[i].src_pitch4=src_pitch_4;
		MT_Data[i].dst1=dstp_1+(MT_Data[i].dst_Y_h_min*dst_pitch_1);
		MT_Data[i].dst2=dstp_2+(MT_Data[i].dst_UV_h_min*dst_pitch_2);
		MT_Data[i].dst3=dstp_3+(MT_Data[i].dst_UV_h_min*dst_pitch_3);
		MT_Data[i].dst4=dstp_4+(MT_Data[i].dst_Y_h_min*dst_pitch_4);
		MT_Data[i].dst_pitch1=dst_pitch_1;
		MT_Data[i].dst_pitch2=dst_pitch_2;
		MT_Data[i].dst_pitch3=dst_pitch_3;
		MT_Data[i].dst_pitch4=dst_pitch_4;
		MT_Data[i].filter_storage_luma=filter_storage_luma;
		MT_Data[i].resampling_program_luma=resampling_program_luma;
		MT_Data[i].resampling_program_chroma=resampling_program_chroma;
		MT_Data[i].filter_storage_chromaU=filter_storage_chroma;
		MT_Data[i].filter_storage_chromaV=filter_storage_chroma;
	}

	if (threads_number>1)
	{
		for(uint8_t i=0; i<threads_number; i++)
			MT_Thread[i].f_process=1;
		if (poolInterface->StartThreads(UserId)) poolInterface->WaitThreadsEnd(UserId);

		if (!grey && vi.IsPlanar() && !isRGBPfamily)
		{
			for(uint8_t i=0; i<threads_number; i++)
				MT_Thread[i].f_process=2;
			if (poolInterface->StartThreads(UserId)) poolInterface->WaitThreadsEnd(UserId);

			for(uint8_t i=0; i<threads_number; i++)
				MT_Thread[i].f_process=3;
			if (poolInterface->StartThreads(UserId)) poolInterface->WaitThreadsEnd(UserId);
		}
		else
		{
			if (isRGBPfamily)
			{
				for(uint8_t i=0; i<threads_number; i++)
					MT_Thread[i].f_process=4;
				if (poolInterface->StartThreads(UserId)) poolInterface->WaitThreadsEnd(UserId);

				for(uint8_t i=0; i<threads_number; i++)
					MT_Thread[i].f_process=5;
				if (poolInterface->StartThreads(UserId)) poolInterface->WaitThreadsEnd(UserId);								
			}
		}

		if (isAlphaChannel)
		{
			for(uint8_t i=0; i<threads_number; i++)
				MT_Thread[i].f_process=6;
			if (poolInterface->StartThreads(UserId)) poolInterface->WaitThreadsEnd(UserId);												
		}

		for(uint8_t i=0; i<threads_number; i++)
			MT_Thread[i].f_process=0;

		poolInterface->ReleaseThreadPool(UserId,sleep);
	}
	else
	{
		// Do resizing
		ResamplerLumaMT(0);
    
		if (!grey && vi.IsPlanar() && !isRGBPfamily)
		{
			// Plane U resizing   
			ResamplerUChromaMT(0);
			// Plane V resizing
			ResamplerVChromaMT(0);
		}
		else
		{
			if (isRGBPfamily)
			{
				// Plane B resizing
				ResamplerLumaMT2(0);
				// Plane R resizing
				ResamplerLumaMT3(0);
			}
		}
		// Plane A resizing
		if (isAlphaChannel) ResamplerLumaMT4(0);
	}

	ReleaseMutex(ghMutex);
	
  return dst;
}


ResamplerH FilteredResizeH::GetResampler(int CPU, bool aligned, int pixelsize, int bits_per_pixel, ResamplingProgram* program, IScriptEnvironment* env)
{
	if (pixelsize==1)
	{
		if ((CPU & CPUF_SSSE3)!=0)
		{
			resize_h_prepare_coeff_8(program,env);
			
			if (program->filter_size > 8) return resizer_h_ssse3_generic;
			else return resizer_h_ssse3_8;
		}
		else return resize_h_c_planar;
	}
	else if (pixelsize==2)
	{ 
		if ((CPU & CPUF_SSSE3)!=0)
		{
			resize_h_prepare_coeff_8(program,env);
#ifdef AVX_BUILD_POSSIBLE				
			if ((CPU & CPUF_AVX2)!=0) return resizer_h_avx2_generic_int16_float_16;
			if ((CPU & CPUF_AVX)!=0) return resizer_h_avx_generic_int16_float_16;
#endif			
			if ((CPU & CPUF_SSE4_1)!=0) return resizer_h_ssse3_generic_int16_float_16<true>;
			else return resizer_h_ssse3_generic_int16_float_16<false>;
		}
		else return resize_h_c_planar_s;
	}
	else
	{ //if (pixelsize == 4)
		if ((CPU & CPUF_SSSE3)!=0)
		{
			resize_h_prepare_coeff_8(program,env);			
#ifdef AVX_BUILD_POSSIBLE
			if ((CPU & CPUF_AVX2)!=0) return resizer_h_avx2_generic_int16_float_32;
			if ((CPU & CPUF_AVX)!=0) return resizer_h_avx_generic_int16_float_32;
#endif			
			return resizer_h_ssse3_generic_int16_float_32;
		}
		else return resize_h_c_planar_f;
	}
}

void FilteredResizeH::FreeData(void) 
{
	mydelete(resampling_program_luma);
	mydelete(resampling_program_chroma);
	
  myalignedfree(filter_storage_luma);
  myalignedfree(filter_storage_chroma);
  myCloseHandle(ghMutex);
}

FilteredResizeH::~FilteredResizeH(void)
{
	if (threads_number>1)
	{
		poolInterface->RemoveUserId(UserId);
		poolInterface->DeAllocateAllThreads(true);
	}
	FreeData();
}


/***************************************
 ***** Filtered Resize - Vertical ******
 ***************************************/

FilteredResizeV::FilteredResizeV( PClip _child, double subrange_top, double subrange_height, int target_height,
	uint8_t _threads,bool _sleep,int range_mode,bool desample,int accuracy,int ChromaS,uint8_t ShiftC,
	bool _avsp,ResamplingFunction* func, IScriptEnvironment* env )
  : GenericVideoFilter(_child),
    resampling_program_luma(NULL), resampling_program_chroma(NULL),
    src_pitch_table_luma(NULL), src_pitch_table_chromaU(NULL), src_pitch_table_chromaV(NULL),
    src_pitch_luma(-1), src_pitch_chromaU(-1), src_pitch_chromaV(-1),
    filter_storage_luma_aligned(NULL), filter_storage_luma_unaligned(NULL),
    filter_storage_chroma_aligned(NULL), filter_storage_chroma_unaligned(NULL),
	sleep(_sleep),threads(_threads),avsp(_avsp)
{
	int16_t i;

    pixelsize = (uint8_t)vi.ComponentSize(); // AVS16
	grey = vi.IsY();
	bits_per_pixel = (uint8_t)vi.BitsPerComponent();
	isRGBPfamily = vi.IsPlanarRGB() || vi.IsPlanarRGBA();
	isAlphaChannel = vi.IsYUVA() || vi.IsPlanarRGBA();
	mode_YUY2 = vi.IsYUY2();

	if ((range_mode!=1) && (range_mode!=4))
	{
		if (vi.IsYUV())
		{
			plane_range[0]=2;
			plane_range[1]=3;
			plane_range[2]=3;
		}
		else
		{
			if (grey)
			{
				for (uint8_t i=0; i<3; i++)
					plane_range[i]=(range_mode==0) ? 2 : range_mode;
			}
			else
			{
				for (uint8_t i=0; i<3; i++)
					plane_range[i]=1;
			}
		}
	}
	else
	{
		if (vi.IsRGB()) range_mode=1;

		for (uint8_t i=0; i<3; i++)
			plane_range[i]=range_mode;
	}
	plane_range[3]=1;
	
    ResampleV_MT=StaticThreadpoolV;
	
	for (i=0; i<MAX_MT_THREADS; i++)
	{
		MT_Thread[i].pClass=this;
		MT_Thread[i].f_process=0;
		MT_Thread[i].thread_Id=(uint8_t)i;
		MT_Thread[i].pFunc=ResampleV_MT;
	}
	UserId=0;
	ghMutex=NULL;
	
	const int shift_w = (!grey && vi.IsPlanar() && !isRGBPfamily) ? vi.GetPlaneWidthSubsampling(PLANAR_U) : 0;
	const int shift_h = (!grey && vi.IsPlanar() && !isRGBPfamily) ? vi.GetPlaneHeightSubsampling(PLANAR_U) : 0;

	const int work_width = vi.IsPlanar() ? vi.width : vi.BytesFromPixels(vi.width)/pixelsize;
	
	if (vi.height<32) threads_number=1;
	else threads_number=threads;

	ghMutex=CreateMutex(NULL,FALSE,NULL);
	if (ghMutex==NULL) env->ThrowError("ResizeVMT: Unable to create Mutex!");


  if (vi.IsRGB() && !isRGBPfamily)
  {
	  if (desample)
		subrange_top = target_height - subrange_top - subrange_height; // why?
	  else
		subrange_top = vi.height - subrange_top - subrange_height; // why?
  }

  // Create resampling program and pitch table
  int SizeV;

  if (ShiftC==0) ShiftC=shift_h;

  if (desample) resampling_program_luma  = func->GetDesamplingProgram(target_height, subrange_top, subrange_height, vi.height, accuracy, ChromaS, ShiftC, SizeV, env);
  else
  {
	  resampling_program_luma  = func->GetResamplingProgram(vi.height, subrange_top, subrange_height, target_height, env);
	  SizeV=target_height;
  }

  threads_number=CreateMTData(threads_number,work_width,vi.height,work_width,SizeV,shift_w,shift_h);

  src_pitch_table_luma = (int *)_aligned_malloc(sizeof(int) * vi.height, 64);
  if (src_pitch_table_luma==NULL)
  {
	  FreeData();
	  env->ThrowError("ResizeVMT: Could not reserve memory in a resampler.");
  }
  
  resampler_luma_aligned   = GetResampler(env->GetCPUFlags(), true ,pixelsize,bits_per_pixel, filter_storage_luma_aligned,   resampling_program_luma);
  resampler_luma_unaligned = GetResampler(env->GetCPUFlags(), false,pixelsize,bits_per_pixel, filter_storage_luma_unaligned, resampling_program_luma);

  if (vi.IsPlanar() && !grey && !isRGBPfamily)
  {
    const int div   = 1 << shift_h;

	if (desample)
	{
		int SizeOut;

	    resampling_program_chroma = func->GetDesamplingProgram(
		                              target_height  >> shift_h,
			                          subrange_top    / div,
				                      subrange_height / div,
					                  vi.height  >> shift_h,
									  accuracy,SizeV,shift_h,SizeOut,
						              env);
	}
	else
	{
	    resampling_program_chroma = func->GetResamplingProgram(
		                              vi.height      >> shift_h,
			                          subrange_top    / div,
				                      subrange_height / div,
					                  target_height  >> shift_h,
						              env);
	}
	src_pitch_table_chromaU = (int *)_aligned_malloc(sizeof(int) * (vi.height >> shift_h), 64);
	src_pitch_table_chromaV = (int *)_aligned_malloc(sizeof(int) * (vi.height >> shift_h), 64);
	if ((src_pitch_table_chromaU==NULL) || (src_pitch_table_chromaV==NULL))
	{
		FreeData();
		env->ThrowError("ResizeVMT: Could not reserve memory in a resampler.");
	}	
    resampler_chroma_aligned   = GetResampler(env->GetCPUFlags(), true ,pixelsize,bits_per_pixel, filter_storage_chroma_aligned,   resampling_program_chroma);
    resampler_chroma_unaligned = GetResampler(env->GetCPUFlags(), false,pixelsize,bits_per_pixel, filter_storage_chroma_unaligned, resampling_program_chroma);
  }

	if (threads_number>1)
	{
		if (!poolInterface->GetUserId(UserId))
		{
			poolInterface->DeAllocateAllThreads(true);
			FreeData();
			env->ThrowError("ResizeVMT: Error with the TheadPool while getting UserId!");
		}
	}

  // Change target video info size
  vi.height = SizeV;
}


int __stdcall FilteredResizeV::SetCacheHints(int cachehints,int frame_range)
{
  switch (cachehints)
  {
  case CACHE_GET_MTMODE :
    return MT_MULTI_INSTANCE;
  default :
    return 0;
  }
}


uint8_t FilteredResizeV::CreateMTData(uint8_t max_threads,int32_t src_size_x,int32_t src_size_y,int32_t dst_size_x,int32_t dst_size_y,int UV_w,int UV_h)
{
	int32_t _y_min,_dh;

	if ((max_threads<=1) || (max_threads>threads_number))
	{
		MT_Data[0].top=true;
		MT_Data[0].bottom=true;
		MT_Data[0].src_Y_h_min=0;
		MT_Data[0].dst_Y_h_min=0;
		MT_Data[0].src_Y_h_max=src_size_y;
		MT_Data[0].dst_Y_h_max=dst_size_y;
		MT_Data[0].src_UV_h_min=0;
		MT_Data[0].dst_UV_h_min=0;
		if (UV_h>0)
		{
			MT_Data[0].src_UV_h_max=src_size_y >> UV_h;
			MT_Data[0].dst_UV_h_max=dst_size_y >> UV_h;
		}
		else
		{
			MT_Data[0].src_UV_h_max=src_size_y;
			MT_Data[0].dst_UV_h_max=dst_size_y;
		}
		MT_Data[0].src_Y_w=src_size_x;
		MT_Data[0].dst_Y_w=dst_size_x;
		if (UV_w>0)
		{
			MT_Data[0].src_UV_w=src_size_x >> UV_w;
			MT_Data[0].dst_UV_w=dst_size_x >> UV_w;
		}
		else
		{
			MT_Data[0].src_UV_w=src_size_x;
			MT_Data[0].dst_UV_w=dst_size_x;
		}
		return(1);
	}

	int32_t src_dh_Y,src_dh_UV,dst_dh_Y,dst_dh_UV;
	int32_t h_y;
	uint8_t i,max=0;

	dst_dh_Y=(dst_size_y+(uint32_t)max_threads-1)/(uint32_t)max_threads;
	if (dst_dh_Y<16) dst_dh_Y=16;
	if ((dst_dh_Y & 3)!=0) dst_dh_Y=((dst_dh_Y+3) >> 2) << 2;

	if (src_size_y==dst_size_y) src_dh_Y=dst_dh_Y;
	else
	{
		src_dh_Y=(src_size_y+(uint32_t)max_threads-1)/(uint32_t)max_threads;
		if (src_dh_Y<16) src_dh_Y=16;
		if ((src_dh_Y & 3)!=0) src_dh_Y=((src_dh_Y+3) >> 2) << 2;
	}

	if (src_size_y<dst_size_y)
	{
		_y_min=src_size_y;
		_dh=src_dh_Y;
	}
	else
	{
		_y_min=dst_size_y;
		_dh=dst_dh_Y;
	}
	h_y=0;
	while (h_y<(_y_min-16))
	{
		max++;
		h_y+=_dh;
	}

	if (max==1)
	{
		MT_Data[0].top=true;
		MT_Data[0].bottom=true;
		MT_Data[0].src_Y_h_min=0;
		MT_Data[0].dst_Y_h_min=0;
		MT_Data[0].src_Y_h_max=src_size_y;
		MT_Data[0].dst_Y_h_max=dst_size_y;
		MT_Data[0].src_UV_h_min=0;
		MT_Data[0].dst_UV_h_min=0;
		if (UV_h>0)
		{
			MT_Data[0].src_UV_h_max=src_size_y >> UV_h;
			MT_Data[0].dst_UV_h_max=dst_size_y >> UV_h;
		}
		else
		{
			MT_Data[0].src_UV_h_max=src_size_y;
			MT_Data[0].dst_UV_h_max=dst_size_y;
		}
		MT_Data[0].src_Y_w=src_size_x;
		MT_Data[0].dst_Y_w=dst_size_x;
		if (UV_w>0)
		{
			MT_Data[0].src_UV_w=src_size_x >> UV_w;
			MT_Data[0].dst_UV_w=dst_size_x >> UV_w;
		}
		else
		{
			MT_Data[0].src_UV_w=src_size_x;
			MT_Data[0].dst_UV_w=dst_size_x;
		}
		return(1);
	}

	src_dh_UV= (UV_h>0) ? src_dh_Y>>UV_h : src_dh_Y;
	dst_dh_UV= (UV_h>0) ? dst_dh_Y>>UV_h : dst_dh_Y;

	MT_Data[0].top=true;
	MT_Data[0].bottom=false;
	MT_Data[0].src_Y_h_min=0;
	MT_Data[0].src_Y_h_max=src_dh_Y;
	MT_Data[0].dst_Y_h_min=0;
	MT_Data[0].dst_Y_h_max=dst_dh_Y;
	MT_Data[0].src_UV_h_min=0;
	MT_Data[0].src_UV_h_max=src_dh_UV;
	MT_Data[0].dst_UV_h_min=0;
	MT_Data[0].dst_UV_h_max=dst_dh_UV;

	i=1;
	while (i<max)
	{
		MT_Data[i].top=false;
		MT_Data[i].bottom=false;
		MT_Data[i].src_Y_h_min=MT_Data[i-1].src_Y_h_max;
		MT_Data[i].src_Y_h_max=MT_Data[i].src_Y_h_min+src_dh_Y;
		MT_Data[i].dst_Y_h_min=MT_Data[i-1].dst_Y_h_max;
		MT_Data[i].dst_Y_h_max=MT_Data[i].dst_Y_h_min+dst_dh_Y;
		MT_Data[i].src_UV_h_min=MT_Data[i-1].src_UV_h_max;
		MT_Data[i].src_UV_h_max=MT_Data[i].src_UV_h_min+src_dh_UV;
		MT_Data[i].dst_UV_h_min=MT_Data[i-1].dst_UV_h_max;
		MT_Data[i].dst_UV_h_max=MT_Data[i].dst_UV_h_min+dst_dh_UV;
		i++;
	}

	MT_Data[max-1].bottom=true;
	MT_Data[max-1].src_Y_h_max=src_size_y;
	MT_Data[max-1].dst_Y_h_max=dst_size_y;
	if (UV_h>0)
	{
		MT_Data[max-1].src_UV_h_max=src_size_y >> UV_h;
		MT_Data[max-1].dst_UV_h_max=dst_size_y >> UV_h;
	}
	else
	{
		MT_Data[max-1].src_UV_h_max=src_size_y;
		MT_Data[max-1].dst_UV_h_max=dst_size_y;
	}

	for (i=0; i<max; i++)
	{
		MT_Data[i].src_Y_w=src_size_x;
		MT_Data[i].dst_Y_w=dst_size_x;
		if (UV_w>0)
		{
			MT_Data[i].src_UV_w=src_size_x >> UV_w;
			MT_Data[i].dst_UV_w=dst_size_x >> UV_w;
		}
		else
		{
			MT_Data[i].src_UV_w=src_size_x;
			MT_Data[i].dst_UV_w=dst_size_x;
		}
	}

	return(max);
}

void FilteredResizeV::ResamplerLumaAlignedMT(uint8_t thread_num)
{
	const MT_Data_Info_ResampleMT mt_data_inf=MT_Data[thread_num];

	resampler_luma_aligned(mt_data_inf.dst1,mt_data_inf.src1,mt_data_inf.dst_pitch1,mt_data_inf.src_pitch1,
		mt_data_inf.resampling_program_luma,mt_data_inf.src_Y_w,bits_per_pixel,mt_data_inf.dst_Y_h_min,mt_data_inf.dst_Y_h_max,
		mt_data_inf.src_pitch_table_luma,mt_data_inf.filter_storage_luma,plane_range[0],mode_YUY2);
}


void FilteredResizeV::ResamplerLumaUnalignedMT(uint8_t thread_num)
{
	const MT_Data_Info_ResampleMT mt_data_inf=MT_Data[thread_num];

	resampler_luma_unaligned(mt_data_inf.dst1,mt_data_inf.src1,mt_data_inf.dst_pitch1,mt_data_inf.src_pitch1,
		mt_data_inf.resampling_program_luma,mt_data_inf.src_Y_w,bits_per_pixel,mt_data_inf.dst_Y_h_min,mt_data_inf.dst_Y_h_max,
		mt_data_inf.src_pitch_table_luma,mt_data_inf.filter_storage_luma,plane_range[0],mode_YUY2);
}

void FilteredResizeV::ResamplerLumaAlignedMT2(uint8_t thread_num)
{
	const MT_Data_Info_ResampleMT mt_data_inf=MT_Data[thread_num];

	resampler_luma_aligned(mt_data_inf.dst2,mt_data_inf.src2,mt_data_inf.dst_pitch2,mt_data_inf.src_pitch2,
		mt_data_inf.resampling_program_luma,mt_data_inf.src_Y_w,bits_per_pixel,mt_data_inf.dst_Y_h_min,mt_data_inf.dst_Y_h_max,
		mt_data_inf.src_pitch_table_luma,mt_data_inf.filter_storage_luma2,plane_range[1],mode_YUY2);
}


void FilteredResizeV::ResamplerLumaUnalignedMT2(uint8_t thread_num)
{
	const MT_Data_Info_ResampleMT mt_data_inf=MT_Data[thread_num];

	resampler_luma_unaligned(mt_data_inf.dst2,mt_data_inf.src2,mt_data_inf.dst_pitch2,mt_data_inf.src_pitch2,
		mt_data_inf.resampling_program_luma,mt_data_inf.src_Y_w,bits_per_pixel,mt_data_inf.dst_Y_h_min,mt_data_inf.dst_Y_h_max,
		mt_data_inf.src_pitch_table_luma,mt_data_inf.filter_storage_luma2,plane_range[1],mode_YUY2);
}


void FilteredResizeV::ResamplerLumaAlignedMT3(uint8_t thread_num)
{
	const MT_Data_Info_ResampleMT mt_data_inf=MT_Data[thread_num];

	resampler_luma_aligned(mt_data_inf.dst3,mt_data_inf.src3,mt_data_inf.dst_pitch3,mt_data_inf.src_pitch3,
		mt_data_inf.resampling_program_luma,mt_data_inf.src_Y_w,bits_per_pixel,mt_data_inf.dst_Y_h_min,mt_data_inf.dst_Y_h_max,
		mt_data_inf.src_pitch_table_luma,mt_data_inf.filter_storage_luma3,plane_range[2],mode_YUY2);
}


void FilteredResizeV::ResamplerLumaUnalignedMT3(uint8_t thread_num)
{
	const MT_Data_Info_ResampleMT mt_data_inf=MT_Data[thread_num];

	resampler_luma_unaligned(mt_data_inf.dst3,mt_data_inf.src3,mt_data_inf.dst_pitch3,mt_data_inf.src_pitch3,
		mt_data_inf.resampling_program_luma,mt_data_inf.src_Y_w,bits_per_pixel,mt_data_inf.dst_Y_h_min,mt_data_inf.dst_Y_h_max,
		mt_data_inf.src_pitch_table_luma,mt_data_inf.filter_storage_luma3,plane_range[2],mode_YUY2);
}


void FilteredResizeV::ResamplerLumaAlignedMT4(uint8_t thread_num)
{
	const MT_Data_Info_ResampleMT mt_data_inf=MT_Data[thread_num];

	resampler_luma_aligned(mt_data_inf.dst4,mt_data_inf.src4,mt_data_inf.dst_pitch4,mt_data_inf.src_pitch4,
		mt_data_inf.resampling_program_luma,mt_data_inf.src_Y_w,bits_per_pixel,mt_data_inf.dst_Y_h_min,mt_data_inf.dst_Y_h_max,
		mt_data_inf.src_pitch_table_luma,mt_data_inf.filter_storage_luma4,plane_range[3],mode_YUY2);
}


void FilteredResizeV::ResamplerLumaUnalignedMT4(uint8_t thread_num)
{
	const MT_Data_Info_ResampleMT mt_data_inf=MT_Data[thread_num];

	resampler_luma_unaligned(mt_data_inf.dst4,mt_data_inf.src4,mt_data_inf.dst_pitch4,mt_data_inf.src_pitch4,
		mt_data_inf.resampling_program_luma,mt_data_inf.src_Y_w,bits_per_pixel,mt_data_inf.dst_Y_h_min,mt_data_inf.dst_Y_h_max,
		mt_data_inf.src_pitch_table_luma,mt_data_inf.filter_storage_luma4,plane_range[3],mode_YUY2);
}


void FilteredResizeV::ResamplerUChromaAlignedMT(uint8_t thread_num)
{
	const MT_Data_Info_ResampleMT mt_data_inf=MT_Data[thread_num];

	resampler_chroma_aligned(mt_data_inf.dst2,mt_data_inf.src2,mt_data_inf.dst_pitch2,mt_data_inf.src_pitch2,
		mt_data_inf.resampling_program_chroma,mt_data_inf.src_UV_w,bits_per_pixel,mt_data_inf.dst_UV_h_min,mt_data_inf.dst_UV_h_max,
		mt_data_inf.src_pitch_table_chromaU,mt_data_inf.filter_storage_chromaU,plane_range[1],mode_YUY2);
}


void FilteredResizeV::ResamplerUChromaUnalignedMT(uint8_t thread_num)
{
	const MT_Data_Info_ResampleMT mt_data_inf=MT_Data[thread_num];

	resampler_chroma_unaligned(mt_data_inf.dst2,mt_data_inf.src2,mt_data_inf.dst_pitch2,mt_data_inf.src_pitch2,
		mt_data_inf.resampling_program_chroma,mt_data_inf.src_UV_w,bits_per_pixel,mt_data_inf.dst_UV_h_min,mt_data_inf.dst_UV_h_max,
		mt_data_inf.src_pitch_table_chromaU,mt_data_inf.filter_storage_chromaU,plane_range[1],mode_YUY2);
}


void FilteredResizeV::ResamplerVChromaAlignedMT(uint8_t thread_num)
{
	const MT_Data_Info_ResampleMT mt_data_inf=MT_Data[thread_num];

	resampler_chroma_aligned(mt_data_inf.dst3,mt_data_inf.src3,mt_data_inf.dst_pitch3,mt_data_inf.src_pitch3,
		mt_data_inf.resampling_program_chroma,mt_data_inf.src_UV_w,bits_per_pixel,mt_data_inf.dst_UV_h_min,mt_data_inf.dst_UV_h_max,
		mt_data_inf.src_pitch_table_chromaV,mt_data_inf.filter_storage_chromaV,plane_range[2],mode_YUY2);
}


void FilteredResizeV::ResamplerVChromaUnalignedMT(uint8_t thread_num)
{
	const MT_Data_Info_ResampleMT mt_data_inf=MT_Data[thread_num];

	resampler_chroma_unaligned(mt_data_inf.dst3,mt_data_inf.src3,mt_data_inf.dst_pitch3,mt_data_inf.src_pitch3,
		mt_data_inf.resampling_program_chroma,mt_data_inf.src_UV_w,bits_per_pixel,mt_data_inf.dst_UV_h_min,mt_data_inf.dst_UV_h_max,
		mt_data_inf.src_pitch_table_chromaV,mt_data_inf.filter_storage_chromaV,plane_range[2],mode_YUY2);
}


void FilteredResizeV::StaticThreadpoolV(void *ptr)
{
	const Public_MT_Data_Thread *data=(const Public_MT_Data_Thread *)ptr;
	FilteredResizeV *ptrClass=(FilteredResizeV *)data->pClass;
	
	switch(data->f_process)
	{
		case 1 : ptrClass->ResamplerLumaAlignedMT(data->thread_Id);
			break;
		case 2 : ptrClass->ResamplerLumaUnalignedMT(data->thread_Id);
			break;
		case 3 : ptrClass->ResamplerUChromaAlignedMT(data->thread_Id);
			break;
		case 4 : ptrClass->ResamplerUChromaUnalignedMT(data->thread_Id);
			break;
		case 5 : ptrClass->ResamplerVChromaAlignedMT(data->thread_Id);
			break;
		case 6 : ptrClass->ResamplerVChromaUnalignedMT(data->thread_Id);
			break;
		case 7 : ptrClass->ResamplerLumaAlignedMT2(data->thread_Id);
			break;
		case 8 : ptrClass->ResamplerLumaUnalignedMT2(data->thread_Id);
			break;			
		case 9 : ptrClass->ResamplerLumaAlignedMT3(data->thread_Id);
			break;
		case 10 : ptrClass->ResamplerLumaUnalignedMT3(data->thread_Id);
			break;			
		case 11 : ptrClass->ResamplerLumaAlignedMT4(data->thread_Id);
			break;
		case 12 : ptrClass->ResamplerLumaUnalignedMT4(data->thread_Id);
			break;			
		default : ;
	}
}


PVideoFrame __stdcall FilteredResizeV::GetFrame(int n, IScriptEnvironment* env)
{
  PVideoFrame src = child->GetFrame(n, env);
  PVideoFrame dst = env->NewVideoFrame(vi);
  
  const int src_pitch_1 = src->GetPitch();
  const int dst_pitch_1 = dst->GetPitch();
  const BYTE *srcp_1 = src->GetReadPtr();
        BYTE *dstp_1 = dst->GetWritePtr();

	const int src_pitch_2 = (!grey && vi.IsPlanar() && !isRGBPfamily) ? src->GetPitch(PLANAR_U) : (isRGBPfamily) ? src->GetPitch(PLANAR_B) : 0;
	const int dst_pitch_2 = (!grey && vi.IsPlanar() && !isRGBPfamily) ? dst->GetPitch(PLANAR_U) : (isRGBPfamily) ? dst->GetPitch(PLANAR_B) : 0;
	const BYTE *srcp_2 = (!grey && vi.IsPlanar() && !isRGBPfamily) ? src->GetReadPtr(PLANAR_U) : (isRGBPfamily) ? src->GetReadPtr(PLANAR_B) : NULL;
	BYTE *dstp_2 = (!grey && vi.IsPlanar() && !isRGBPfamily) ? dst->GetWritePtr(PLANAR_U) : (isRGBPfamily) ? dst->GetWritePtr(PLANAR_B) : NULL;

	const int src_pitch_3 = (!grey && vi.IsPlanar() && !isRGBPfamily) ? src->GetPitch(PLANAR_V) : (isRGBPfamily) ? src->GetPitch(PLANAR_R) : 0;
	const int dst_pitch_3 = (!grey && vi.IsPlanar() && !isRGBPfamily) ? dst->GetPitch(PLANAR_V) : (isRGBPfamily) ? dst->GetPitch(PLANAR_R) : 0;
	const BYTE *srcp_3 = (!grey && vi.IsPlanar() && !isRGBPfamily) ? src->GetReadPtr(PLANAR_V) : (isRGBPfamily) ? src->GetReadPtr(PLANAR_R) : NULL;
	BYTE *dstp_3 = (!grey && vi.IsPlanar() && !isRGBPfamily) ? dst->GetWritePtr(PLANAR_V) : (isRGBPfamily) ? dst->GetWritePtr(PLANAR_R) : NULL;
	
	const int src_pitch_4 = (isAlphaChannel) ? src->GetPitch(PLANAR_A) : 0;
	const int dst_pitch_4 = (isAlphaChannel) ? dst->GetPitch(PLANAR_A) : 0;
	const BYTE *srcp_4 = (isAlphaChannel) ? src->GetReadPtr(PLANAR_A) : NULL;
	BYTE *dstp_4 = (isAlphaChannel) ? dst->GetWritePtr(PLANAR_A) : NULL;  
 
  WaitForSingleObject(ghMutex,INFINITE);

  // Create pitch table
  if (src_pitch_luma != src->GetPitch())
  {
    src_pitch_luma = src->GetPitch();
    resize_v_create_pitch_table(src_pitch_table_luma, src_pitch_luma, src->GetHeight(),pixelsize);
  }

  if (!grey && vi.IsPlanar() && !isRGBPfamily)
  {
	if (src_pitch_chromaU != src->GetPitch(PLANAR_U))
	{
		src_pitch_chromaU = src->GetPitch(PLANAR_U);
		resize_v_create_pitch_table(src_pitch_table_chromaU, src_pitch_chromaU, src->GetHeight(PLANAR_U),pixelsize);
	}	  
	if (src_pitch_chromaV != src->GetPitch(PLANAR_V))
	{
		src_pitch_chromaV = src->GetPitch(PLANAR_V);
		resize_v_create_pitch_table(src_pitch_table_chromaV, src_pitch_chromaV, src->GetHeight(PLANAR_V),pixelsize);
	}	
  }

  if (threads_number>1)
  {
	if (!poolInterface->RequestThreadPool(UserId,threads_number,MT_Thread,-1,false))
	{
		ReleaseMutex(ghMutex);
		env->ThrowError("ResizeHMT: Error with the TheadPool while requesting threadpool!");
	}
  }

	for(uint8_t i=0; i<threads_number; i++)
	{		
		MT_Data[i].src1=srcp_1;
		MT_Data[i].src2=srcp_2;
		MT_Data[i].src3=srcp_3;
		MT_Data[i].src4=srcp_4;
		MT_Data[i].src_pitch1=src_pitch_1;
		MT_Data[i].src_pitch2=src_pitch_2;
		MT_Data[i].src_pitch3=src_pitch_3;
		MT_Data[i].src_pitch4=src_pitch_4;
		MT_Data[i].dst1=dstp_1+(MT_Data[i].dst_Y_h_min*dst_pitch_1);
		MT_Data[i].dst2=dstp_2+(MT_Data[i].dst_UV_h_min*dst_pitch_2);
		MT_Data[i].dst3=dstp_3+(MT_Data[i].dst_UV_h_min*dst_pitch_3);
		MT_Data[i].dst4=dstp_4+(MT_Data[i].dst_Y_h_min*dst_pitch_4);
		MT_Data[i].dst_pitch1=dst_pitch_1;
		MT_Data[i].dst_pitch2=dst_pitch_2;
		MT_Data[i].dst_pitch3=dst_pitch_3;
		MT_Data[i].dst_pitch4=dst_pitch_4;
		if (IsPtrAligned(srcp_1, 16) && ((src_pitch_1 & 15) == 0))
			MT_Data[i].filter_storage_luma=filter_storage_luma_aligned;
		else
			MT_Data[i].filter_storage_luma=filter_storage_luma_unaligned;
		if (IsPtrAligned(srcp_4, 16) && ((src_pitch_4 & 15) == 0))
			MT_Data[i].filter_storage_luma4=filter_storage_luma_aligned;
		else
			MT_Data[i].filter_storage_luma4=filter_storage_luma_unaligned;
		MT_Data[i].src_pitch_table_luma=src_pitch_table_luma;
		MT_Data[i].src_pitch_table_chromaU=src_pitch_table_chromaU;
		MT_Data[i].src_pitch_table_chromaV=src_pitch_table_chromaV;
		MT_Data[i].resampling_program_luma=resampling_program_luma;
		MT_Data[i].resampling_program_chroma=resampling_program_chroma;
		if (IsPtrAligned(srcp_2, 16) && ((src_pitch_2 & 15) == 0))
		{
			MT_Data[i].filter_storage_chromaU=filter_storage_chroma_aligned;
			MT_Data[i].filter_storage_luma2=filter_storage_luma_aligned;
		}
		else
		{
			MT_Data[i].filter_storage_chromaU=filter_storage_chroma_unaligned;
			MT_Data[i].filter_storage_luma2=filter_storage_luma_unaligned;
		}
		if (IsPtrAligned(srcp_3, 16) && ((src_pitch_3 & 15) == 0))
		{
			MT_Data[i].filter_storage_chromaV=filter_storage_chroma_aligned;
			MT_Data[i].filter_storage_luma3=filter_storage_luma_aligned;
		}
		else
		{
			MT_Data[i].filter_storage_chromaV=filter_storage_chroma_unaligned;
			MT_Data[i].filter_storage_luma3=filter_storage_luma_unaligned;
		}
	}

	if (threads_number>1)
	{
		uint8_t f_proc;

		if (IsPtrAligned(srcp_1, 16) && ((src_pitch_1 & 15) == 0)) f_proc=1;
		else f_proc=2;

		for(uint8_t i=0; i<threads_number; i++)
			MT_Thread[i].f_process=f_proc;
		if (poolInterface->StartThreads(UserId)) poolInterface->WaitThreadsEnd(UserId);

		if (!grey && vi.IsPlanar() && !isRGBPfamily)
		{
			if (IsPtrAligned(srcp_2, 16) && ((src_pitch_2 & 15) == 0)) f_proc=3;
			else f_proc=4;

			for(uint8_t i=0; i<threads_number; i++)
				MT_Thread[i].f_process=f_proc;
			if (poolInterface->StartThreads(UserId)) poolInterface->WaitThreadsEnd(UserId);

			if (IsPtrAligned(srcp_3, 16) && ((src_pitch_3 & 15) == 0)) f_proc=5;
			else f_proc=6;

			for(uint8_t i=0; i<threads_number; i++)
				MT_Thread[i].f_process=f_proc;
			if (poolInterface->StartThreads(UserId)) poolInterface->WaitThreadsEnd(UserId);
		}
		else
		{
			if (isRGBPfamily)
			{
				if (IsPtrAligned(srcp_2, 16) && ((src_pitch_2 & 15) == 0)) f_proc=7;
				else f_proc=8;

				for(uint8_t i=0; i<threads_number; i++)
					MT_Thread[i].f_process=f_proc;
				if (poolInterface->StartThreads(UserId)) poolInterface->WaitThreadsEnd(UserId);

				if (IsPtrAligned(srcp_3, 16) && ((src_pitch_3 & 15) == 0)) f_proc=9;
				else f_proc=10;

				for(uint8_t i=0; i<threads_number; i++)
					MT_Thread[i].f_process=f_proc;
				if (poolInterface->StartThreads(UserId)) poolInterface->WaitThreadsEnd(UserId);							
			}
		}
		
		if (isAlphaChannel)
		{
			if (IsPtrAligned(srcp_4, 16) && ((src_pitch_4 & 15) == 0)) f_proc=11;
			else f_proc=12;

			for(uint8_t i=0; i<threads_number; i++)
				MT_Thread[i].f_process=f_proc;
			if (poolInterface->StartThreads(UserId)) poolInterface->WaitThreadsEnd(UserId);			
		}

		for(uint8_t i=0; i<threads_number; i++)
			MT_Thread[i].f_process=0;

		poolInterface->ReleaseThreadPool(UserId,sleep);
	}
	else
	{
		// Do resizing
		if (IsPtrAligned(srcp_1, 16) && ((src_pitch_1 & 15) == 0))
			ResamplerLumaAlignedMT(0);
		else
			ResamplerLumaUnalignedMT(0);
    
		if (!grey && vi.IsPlanar() && !isRGBPfamily)
		{
			// Plane U resizing   
			if (IsPtrAligned(srcp_2, 16) && ((src_pitch_2 & 15) == 0))
				ResamplerUChromaAlignedMT(0);
			else
				ResamplerUChromaUnalignedMT(0);

			// Plane V resizing
			if (IsPtrAligned(srcp_3, 16) && ((src_pitch_3 & 15) == 0))
				ResamplerVChromaAlignedMT(0);
			else
				ResamplerVChromaUnalignedMT(0);
		}
		else
		{
			if (isRGBPfamily)
			{
				if (IsPtrAligned(srcp_2, 16) && ((src_pitch_2 & 15) == 0))
					ResamplerLumaAlignedMT2(0);
				else
					ResamplerLumaUnalignedMT2(0);		
				
				if (IsPtrAligned(srcp_3, 16) && ((src_pitch_3 & 15) == 0))
					ResamplerLumaAlignedMT3(0);
				else
					ResamplerLumaUnalignedMT3(0);								
			}			
		}
		
		if (isAlphaChannel)
		{
			if (IsPtrAligned(srcp_4, 16) && ((src_pitch_4 & 15) == 0))
				ResamplerLumaAlignedMT4(0);
			else
				ResamplerLumaUnalignedMT4(0);	
		}
	}

	ReleaseMutex(ghMutex);

  return dst;
}

ResamplerV FilteredResizeV::GetResampler(int CPU, bool aligned,int pixelsize, int bits_per_pixel, void*& storage, ResamplingProgram* program)
{
  if (program->filter_size == 1)
  {
    // Fast pointresize
    switch (pixelsize) // AVS16
    {
    case 1: return resize_v_planar_pointresize<uint8_t>;
    case 2: return resize_v_planar_pointresize<uint16_t>;
    default: // case 4:
      return resize_v_planar_pointresize<float>;
    }
  }
  else
  {
    // Other resizers
    if (pixelsize == 1)
    {
      if ((CPU & CPUF_SSSE3)!=0)
	  {
        if (aligned && ((CPU & CPUF_SSE4_1)!=0))
		{
          return resize_v_ssse3_planar<simd_load_streaming>;
        }
        else if (aligned)
		{ // SSSE3 aligned
          return resize_v_ssse3_planar<simd_load_aligned>;
        }
        else if ((CPU & CPUF_SSE3)!=0)
		{ // SSE3 lddqu
          return resize_v_ssse3_planar<simd_load_unaligned_sse3>;
        }
        else
		{ // SSSE3 unaligned
          return resize_v_ssse3_planar<simd_load_unaligned>;
        }
      }
      else if ((CPU & CPUF_SSE2)!=0)
	  {
        if (aligned && ((CPU & CPUF_SSE4_1)!=0))
		{ // SSE4.1 movntdqa constantly provide ~2% performance increase in my testing
          return resize_v_sse2_planar<simd_load_streaming>;
        }
        else if (aligned)
		{ // SSE2 aligned
          return resize_v_sse2_planar<simd_load_aligned>;
        }
        else if ((CPU & CPUF_SSE3)!=0)
		{ // SSE2 lddqu
          return resize_v_sse2_planar<simd_load_unaligned_sse3>;
        }
        else
		{ // SSE2 unaligned
          return resize_v_sse2_planar<simd_load_unaligned>;
        }
#ifdef X86_32
      }
      else if ((CPU & CPUF_MMX)!=0)
	  {
        return resize_v_mmx_planar;
#endif
      }
      else { // C version
        return resize_v_c_planar;
      }
    } 
    else if (pixelsize == 2)
	{
#ifdef AVX_BUILD_POSSIBLE		
		if (aligned && ((CPU & CPUF_AVX2)!=0)) return resize_v_avx2_planar_16;
		if (aligned && ((CPU & CPUF_AVX)!=0)) return resize_v_avx_planar_16;
		else
#endif			
		if ((CPU & CPUF_SSE4_1)!=0)
		{
			if (aligned)
			{
				return resize_v_sseX_planar_16<simd_load_streaming, true>;
			}
			else if ((CPU & CPUF_SSE3)!=0)
			{ // SSE3 lddqu
				return resize_v_sseX_planar_16<simd_load_unaligned_sse3, true>;
			}
			else
			{ // unaligned
				return resize_v_sseX_planar_16<simd_load_unaligned, true>;
			}
		}
		else if ((CPU & CPUF_SSE2)!=0)
		{
			if (aligned)
			{
				return resize_v_sseX_planar_16<simd_load_aligned, false>;
			}
			else
			{
				return resize_v_sseX_planar_16<simd_load_unaligned, false>;
			}
		}
		else
		{ // C version
			return resize_v_c_planar_s;
		}
    }
    else
	{ // if (pixelsize== 4) 
#ifdef AVX_BUILD_POSSIBLE			
		if (aligned && ((CPU & CPUF_AVX2)!=0)) return resize_v_avx2_planar_32;
		if (aligned && ((CPU & CPUF_AVX)!=0)) return resize_v_avx_planar_32;
		else
#endif			
		if ((CPU & CPUF_SSE2)!=0)
		{
			if (aligned)
			{
				return resize_v_sseX_planar_32<simd_loadps_aligned>;
			}
			else
			{
				return resize_v_sseX_planar_32<simd_loadps_unaligned>;
			}
		}
		else return resize_v_c_planar_f;
    }
  }
}


void FilteredResizeV::FreeData(void) 
{
	mydelete(resampling_program_luma);
	mydelete(resampling_program_chroma);
	myalignedfree(src_pitch_table_luma);
	myalignedfree(src_pitch_table_chromaU);
	myalignedfree(src_pitch_table_chromaV);

  myalignedfree(filter_storage_luma_aligned);
  myalignedfree(filter_storage_luma_unaligned);
  myalignedfree(filter_storage_chroma_aligned);
  myalignedfree(filter_storage_chroma_unaligned);
  myCloseHandle(ghMutex);
}


FilteredResizeV::~FilteredResizeV(void)
{
	if (threads_number>1)
	{
		poolInterface->RemoveUserId(UserId);
		poolInterface->DeAllocateAllThreads(true);
	}
	FreeData();
}


/**********************************************
 *******   Resampling Factory Methods   *******
 **********************************************/


PClip FilteredResizeMT::CreateResizeV(PClip clip, double subrange_top, double subrange_height, int target_height,
                    uint8_t _threads,bool _sleep,int range_mode,bool desample,int accuracy,int ChromaS,uint8_t ShiftC,
					bool _avsp,ResamplingFunction* func,IScriptEnvironment* env)
{
  return new FilteredResizeV(clip, subrange_top, subrange_height, target_height,_threads,_sleep,range_mode,desample,
	  accuracy,ChromaS,ShiftC,_avsp,func, env);
}


PClip FilteredResizeMT::CreateResizeH(PClip clip, double subrange_top, double subrange_height, int target_height,
                    uint8_t _threads,bool _sleep,int range_mode,bool desample,int accuracy,bool _avsp,ResamplingFunction* func,
					IScriptEnvironment* env)
{
  return new FilteredResizeH(clip, subrange_top, subrange_height, target_height,_threads,_sleep,range_mode,desample,
	  accuracy,_avsp,func, env);
}


PClip FilteredResizeMT::CreateResize(PClip clip, int target_width, int target_height,int _threads,
	bool _LogicalCores,bool _MaxPhysCores, bool _SetAffinity,bool _sleep,int prefetch,
	int range_mode,bool desample,int accuracy,int order,const AVSValue* args,ResamplingFunction* f,
	IScriptEnvironment* env)
{
  const VideoInfo& vi = clip->GetVideoInfo();
  const double subrange_left = args[0].AsFloat(0), subrange_top = args[1].AsFloat(0);
  
  if (target_height <= 0)
    env->ThrowError("ResizeMT: Height must be greater than 0.");

  if (target_width <= 0) {
    env->ThrowError("ResizeMT: Width must be greater than 0.");
  }

  if ((range_mode<0) || (range_mode>4)) env->ThrowError("ResizeMT: [range] must be between 0 and 4.");

  if ((prefetch<0) || (prefetch>MAX_THREAD_POOL)) env->ThrowError("ResizeMT: [prefetch] must be between 0 and %d.",MAX_THREAD_POOL);
  if (prefetch==0) prefetch=1;
  if ((_threads<0) || (_threads>MAX_MT_THREADS)) env->ThrowError("ResizeMT: [threads] must be between 0 and %d.",MAX_MT_THREADS);
  if ((accuracy<0) || (accuracy>2)) env->ThrowError("ResizeMT: [accuracy] must be between 0 and 2.");
  if ((order<0) || (order>2)) env->ThrowError("ResizeMT: [order] must be between 0 and 2.");

  const bool avsp=env->FunctionExists("ConvertBits");  
  const bool isRGBPfamily = vi.IsPlanarRGB() || vi.IsPlanarRGBA();
  const bool grey = vi.IsY();  
  const bool isAlphaChannel = vi.IsYUVA() || vi.IsPlanarRGBA();

  if (vi.IsPlanar() && !grey && !isRGBPfamily)
  {
    int  mask;
	
	mask = (1 << vi.GetPlaneHeightSubsampling(PLANAR_U)) - 1;
    if ((target_height & mask)!=0)
      env->ThrowError("ResizeMT: Planar destination height must be a multiple of %d.", mask+1);
  
    mask = (1 << vi.GetPlaneWidthSubsampling(PLANAR_U)) - 1;
    if ((target_width & mask)!=0)
      env->ThrowError("ResizeMT: Planar destination width must be a multiple of %d.", mask+1);
  
  }

  double subrange_width=(desample)?args[2].AsDblDef(target_width):args[2].AsDblDef(vi.width);
  double subrange_height=(desample)?args[3].AsDblDef(target_height):args[3].AsDblDef(vi.height);

  // Crop style syntax
  if (desample)
  {
	if (subrange_width  <= 0.0) subrange_width  = target_width  - subrange_left + subrange_width;
	if (subrange_height <= 0.0) subrange_height = target_height - subrange_top  + subrange_height;
  }
  else
  {
	if (subrange_width  <= 0.0) subrange_width  = vi.width  - subrange_left + subrange_width;
	if (subrange_height <= 0.0) subrange_height = vi.height - subrange_top  + subrange_height;
  }

  const int shift = (!grey && vi.IsPlanar() && !isRGBPfamily) ? vi.GetPlaneWidthSubsampling(PLANAR_U) : 0;
  int SizeH;

  if (desample) SizeH=f->GetDesamplingData(target_width, subrange_left, subrange_width, vi.width,shift,env);
  else SizeH=target_width;

  bool fast_resize=((env->GetCPUFlags() & CPUF_SSSE3) == CPUF_SSSE3 ) && vi.IsPlanar() && ((SizeH & 3) == 0);

	if (fast_resize && !grey && !isRGBPfamily)
	{
		const int dst_chroma_width = SizeH >> shift;

		if ((dst_chroma_width & 3) != 0) fast_resize = false;
	}  

  PClip result;
  // ensure that the intermediate area is maximal
  const double area_FirstH = (desample)?subrange_height*vi.width:subrange_height*target_width;
  const double area_FirstV = (desample)?subrange_width*target_height:subrange_width*vi.height;

  bool VFirst;

  if (desample)
  {
	  switch(order)
	  {
		case 0 : VFirst=(area_FirstH>=area_FirstV); break;
		case 1 : VFirst=true; break;
		case 2 : VFirst=false; break;
		default : VFirst=(area_FirstH>=area_FirstV); break;
	  }
  }
  else VFirst=(area_FirstH<area_FirstV);

  const bool FTurnL=(!avsp) && (env->FunctionExists("FTurnLeft") && ((env->GetCPUFlags() & CPUF_SSE2)!=0)) && (!vi.IsRGB());
  const bool FTurnR=(!avsp) && (env->FunctionExists("FTurnRight") && ((env->GetCPUFlags() & CPUF_SSE2)!=0)) && (!vi.IsRGB());

  auto turnRightFunction = (FTurnR) ? "FTurnRight" : "TurnRight";
  auto turnLeftFunction =  (FTurnL) ? "FTurnLeft" : "TurnLeft";

  uint8_t plane_range[4];

	if ((range_mode!=1) && (range_mode!=4))
	{
		if (vi.IsYUV())
		{
			plane_range[0]=2;
			plane_range[1]=3;
			plane_range[2]=3;
		}
		else
		{
			if (grey)
			{
				for (uint8_t i=0; i<3; i++)
					plane_range[i]=(range_mode==0) ? 2 : range_mode;
			}
			else
			{
				for (uint8_t i=0; i<3; i++)
					plane_range[i]=1;
			}
		}
	}
	else
	{
		if (vi.IsRGB()) range_mode=1;

		for (uint8_t i=0; i<3; i++)
			plane_range[i]=range_mode;
	}
	plane_range[3]=1;

	bool step1,step2;

  bool CropV,CropH;

  if (desample)
  {
	  CropV=((subrange_top==int(subrange_top)) && (subrange_height==vi.height)
		&& (subrange_top>=0) && ((subrange_top+subrange_height)<= target_height));
	  CropH=((subrange_left==int(subrange_left)) && (subrange_width==vi.width)
		  && (subrange_left>=0) && ((subrange_left+subrange_width)<=target_width));
  }
  else
  {
	  CropV=((subrange_top==int(subrange_top)) && (subrange_height==target_height)
		&& (subrange_top>=0) && ((subrange_top+subrange_height)<= vi.height));
	  CropH=((subrange_left==int(subrange_left)) && (subrange_width==target_width)
		  && (subrange_left>=0) && ((subrange_left+subrange_width)<=vi.width));
  }


  if (VFirst)
  {
	if ((subrange_top==0) && (subrange_height==target_height) && (subrange_height==vi.height)) step1=false;
	else
	{
		if (CropV)
		{
			const int mask = (vi.IsPlanar() && !grey && !isRGBPfamily) ? (1 << vi.GetPlaneHeightSubsampling(PLANAR_U)) - 1 : 0;

			if (((int(subrange_top) | int(subrange_height)) & mask) == 0) step1=false;
			else step1=true;
		}
		else step1=true;
	}
	if (!((subrange_left==0) && (subrange_width==target_width) && (subrange_width==vi.width)))
	{
		if (CropH)
		{
			const int mask = (vi.IsPlanar() && !grey && !isRGBPfamily) ? (1 << vi.GetPlaneWidthSubsampling(PLANAR_U)) - 1 : 0;

		    if (((int(subrange_left) | int(subrange_width)) & mask) == 0) step2=false;
			else step2=true;
		}
		else step2=true;
	}
	else step2=false;
  }
  else
  {
	if ((subrange_left==0) && (subrange_width==target_width) && (subrange_width==vi.width)) step1=false;
	else
	{
		if (CropH)
		{
			const int mask = (vi.IsPlanar() && !grey && !isRGBPfamily) ? (1 << vi.GetPlaneWidthSubsampling(PLANAR_U)) - 1 : 0;

		    if (((int(subrange_left) | int(subrange_width)) & mask) == 0) step1=false;
			else step1=true;
		}
		else step1=true;
	}
	if (!((subrange_top==0) && (subrange_height==target_height) && (subrange_height==vi.height)))
	{
		if (CropV)
		{
			const int mask = (vi.IsPlanar() && !grey && !isRGBPfamily) ? (1 << vi.GetPlaneHeightSubsampling(PLANAR_U)) - 1 : 0;
				
			if (((int(subrange_top) | int(subrange_height)) & mask) == 0) step2=false;
			else step2=true;
		}
		else step2=true;
	}
	else step2=false;
  }


	uint8_t threads_number=1;

	if ((_threads!=1) && (step1 || step2))
	{
		if (!poolInterface->CreatePool(prefetch)) env->ThrowError("ResizeMT: Unable to create ThreadPool!");

		threads_number=poolInterface->GetThreadNumber(_threads,_LogicalCores);

		if (threads_number==0) env->ThrowError("ResizeMT: Error with the TheadPool while getting CPU info!");

		if (threads_number>1)
		{
			if (!poolInterface->AllocateThreads(threads_number,0,0,_MaxPhysCores,_SetAffinity,true,-1))
				env->ThrowError("ResizeMT: Error with the TheadPool while allocating threadpool!");
		}
	}


  if (!fast_resize)
  {
	  if (VFirst)
	  {
		if ((subrange_top==0) && (subrange_height==target_height) && (subrange_height==vi.height)) result=clip;
		else
		{
			if (CropV)
			{
				const int mask = (vi.IsPlanar() && !grey && !isRGBPfamily) ? (1 << vi.GetPlaneHeightSubsampling(PLANAR_U)) - 1 : 0;

				if (((int(subrange_top) | int(subrange_height)) & mask) == 0)
				{
					if (desample) result=clip;
					else
					{
						AVSValue sargs[6] = {clip,0,int(subrange_top),vi.width,int(subrange_height),0};
						result=env->Invoke("Crop",AVSValue(sargs,6)).AsClip();
					}
				}
				else result = CreateResizeV(clip, subrange_top, subrange_height, target_height,threads_number,_sleep,(step2)?1:range_mode,desample,accuracy,0,0,avsp, f, env);
			}
			else result = CreateResizeV(clip, subrange_top, subrange_height, target_height,threads_number,_sleep,(step2)?1:range_mode,desample,accuracy,0,0,avsp, f, env);
		}
		if (!((subrange_left==0) && (subrange_width==target_width) && (subrange_width==vi.width)))
		{
			if (CropH)
			{
				const int mask = (vi.IsPlanar() && !grey && !isRGBPfamily) ? (1 << vi.GetPlaneWidthSubsampling(PLANAR_U)) - 1 : 0;

			    if (((int(subrange_left) | int(subrange_width)) & mask) == 0)
				{
					if (!desample)
					{
						AVSValue sargs[6] = {result,int(subrange_left),0,int(subrange_width),vi.height,0};
						result=env->Invoke("Crop",AVSValue(sargs,6)).AsClip();
					}
				}
				else
				{
					if (!vi.IsRGB() || isRGBPfamily)
					{
						if (vi.Is422() || vi.IsYUY2() || vi.IsYV411())
						{
							const int shift = vi.GetPlaneWidthSubsampling(PLANAR_U);
							const int div   = 1 << shift;

							AVSValue v,vv,vu,va;
							
							if (avsp)
							{
								AVSValue sargs[2] = {result,"Y"};
								
								v=env->Invoke("PlaneToY",AVSValue(sargs,2)).AsClip();
								sargs[1]="U";
								vu=env->Invoke("PlaneToY",AVSValue(sargs,2)).AsClip();
								sargs[1]="V";
								vv=env->Invoke("PlaneToY",AVSValue(sargs,2)).AsClip();
								if (isAlphaChannel)
								{
									sargs[1]="A";
									va=env->Invoke("PlaneToY",AVSValue(sargs,2)).AsClip();
								}
							}
							else
							{
								vu = env->Invoke("UtoY8",result).AsClip();
								vv = env->Invoke("VtoY8",result).AsClip();
								v = env->Invoke("ConvertToY8",result).AsClip();								
							}
							
							v = env->Invoke(turnRightFunction,v).AsClip();
							vu = env->Invoke(turnRightFunction,vu).AsClip();
							vv = env->Invoke(turnRightFunction,vv).AsClip();
							if (isAlphaChannel) va = env->Invoke(turnRightFunction,va).AsClip();
							v = CreateResizeV(v.AsClip(), subrange_left, subrange_width, target_width,threads_number,_sleep,plane_range[0],desample,accuracy,0,shift,avsp, f, env);

							VideoInfo vR = v.AsClip()->GetVideoInfo();
							int ChromaS=vR.height >> shift;

							vu = CreateResizeV(vu.AsClip(), subrange_left/div, subrange_width/div, target_width >> shift,threads_number,_sleep,plane_range[1],desample,accuracy,ChromaS,0,avsp, f, env);
							vv = CreateResizeV(vv.AsClip(), subrange_left/div, subrange_width/div, target_width >> shift,threads_number,_sleep,plane_range[2],desample,accuracy,ChromaS,0,avsp, f, env);
							if (isAlphaChannel) va = CreateResizeV(va.AsClip(), subrange_left, subrange_width, target_width,threads_number,_sleep,plane_range[3],desample,accuracy,0,shift,avsp, f, env);
							v = env->Invoke(turnLeftFunction,v).AsClip();
							vu = env->Invoke(turnLeftFunction,vu).AsClip();
							vv = env->Invoke(turnLeftFunction,vv).AsClip();
							if (isAlphaChannel) va = env->Invoke(turnLeftFunction,va).AsClip();

						    AVSValue ytouvargs[4] = {vu,vv,v,va};
						    if (isAlphaChannel) result=env->Invoke("YtoUV",AVSValue(ytouvargs,4)).AsClip();
							else result=env->Invoke("YtoUV",AVSValue(ytouvargs,3)).AsClip();
						    if (vi.IsYUY2()) result=env->Invoke("ConvertToYUY2",result).AsClip();
						}
						else
						{
							result=env->Invoke(turnRightFunction,result).AsClip();
							result=CreateResizeV(result, subrange_left, subrange_width, target_width,threads_number,_sleep,range_mode,desample,accuracy,0,0,avsp, f, env);
							result=env->Invoke(turnLeftFunction,result).AsClip();
						}
					}
					else
					{
						result=env->Invoke(turnLeftFunction,result).AsClip();
						result=CreateResizeV(result, subrange_left, subrange_width, target_width,threads_number,_sleep,range_mode,desample,accuracy,0,0,avsp, f, env);
						result=env->Invoke(turnRightFunction,result).AsClip();
					}
				}
			}
			else
			{
				if (!vi.IsRGB() || isRGBPfamily)
				{
				    if (vi.Is422() || vi.IsYUY2() || vi.IsYV411())
					{
						const int shift = vi.GetPlaneWidthSubsampling(PLANAR_U);
						const int div   = 1 << shift;

						AVSValue v,vv,vu,va;
						
						if (avsp)
						{
							AVSValue sargs[2] = {result,"Y"};
								
							v=env->Invoke("PlaneToY",AVSValue(sargs,2)).AsClip();
							sargs[1]="U";
							vu=env->Invoke("PlaneToY",AVSValue(sargs,2)).AsClip();
							sargs[1]="V";
							vv=env->Invoke("PlaneToY",AVSValue(sargs,2)).AsClip();
							if (isAlphaChannel)
							{
								sargs[1]="A";
								va=env->Invoke("PlaneToY",AVSValue(sargs,2)).AsClip();
							}
						}
						else
						{
							vu = env->Invoke("UtoY8",result).AsClip();
							vv = env->Invoke("VtoY8",result).AsClip();
							v = env->Invoke("ConvertToY8",result).AsClip();								
						}
			
						v = env->Invoke(turnRightFunction,v).AsClip();
						vu = env->Invoke(turnRightFunction,vu).AsClip();
						vv = env->Invoke(turnRightFunction,vv).AsClip();
						if (isAlphaChannel) va = env->Invoke(turnRightFunction,va).AsClip();
						v = CreateResizeV(v.AsClip(), subrange_left, subrange_width, target_width,threads_number,_sleep,plane_range[0],desample,accuracy,0,shift,avsp, f, env);

						VideoInfo vR = v.AsClip()->GetVideoInfo();
						int ChromaS=vR.height >> shift;

						vu = CreateResizeV(vu.AsClip(), subrange_left/div, subrange_width/div, target_width >> shift,threads_number,_sleep,plane_range[1],desample,accuracy,ChromaS,0,avsp, f, env);
						vv = CreateResizeV(vv.AsClip(), subrange_left/div, subrange_width/div, target_width >> shift,threads_number,_sleep,plane_range[2],desample,accuracy,ChromaS,0,avsp, f, env);
						if (isAlphaChannel) va = CreateResizeV(va.AsClip(), subrange_left, subrange_width, target_width,threads_number,_sleep,plane_range[3],desample,accuracy,0,shift,avsp, f, env);
						v = env->Invoke(turnLeftFunction,v).AsClip();
						vu = env->Invoke(turnLeftFunction,vu).AsClip();
						vv = env->Invoke(turnLeftFunction,vv).AsClip();
						if (isAlphaChannel) va = env->Invoke(turnLeftFunction,va).AsClip();
							
					    AVSValue ytouvargs[4] = {vu,vv,v,va};
					    if (isAlphaChannel) result=env->Invoke("YtoUV",AVSValue(ytouvargs,4)).AsClip();
						else result=env->Invoke("YtoUV",AVSValue(ytouvargs,3)).AsClip();
					    if (vi.IsYUY2()) result=env->Invoke("ConvertToYUY2",result).AsClip();
					}
					else
					{
						result=env->Invoke(turnRightFunction,result).AsClip();
						result=CreateResizeV(result, subrange_left, subrange_width, target_width,threads_number,_sleep,range_mode,desample,accuracy,0,0,avsp, f, env);
						result=env->Invoke(turnLeftFunction,result).AsClip();
					}
				}
				else
				{
					result=env->Invoke(turnLeftFunction,result).AsClip();
					result=CreateResizeV(result, subrange_left, subrange_width, target_width,threads_number,_sleep,range_mode,desample,accuracy,0,0,avsp, f, env);
					result=env->Invoke(turnRightFunction,result).AsClip();
				}
			}
		}
	  }
	  else
	  {
		if ((subrange_left==0) && (subrange_width==target_width) && (subrange_width==vi.width)) result=clip;
		else
		{
			if (CropH)
			{
				const int mask = (vi.IsPlanar() && !grey && !isRGBPfamily) ? (1 << vi.GetPlaneWidthSubsampling(PLANAR_U)) - 1 : 0;

			    if (((int(subrange_left) | int(subrange_width)) & mask) == 0)
				{
					if (desample) result=clip;
					else
					{
						AVSValue sargs[6] = {clip,int(subrange_left),0,int(subrange_width),vi.height,0};
						result=env->Invoke("Crop",AVSValue(sargs,6)).AsClip();
					}
				}
				else
				{
					if (!vi.IsRGB() || isRGBPfamily)
					{
					    if (vi.Is422() || vi.IsYUY2() || vi.IsYV411())
						{
							const int shift = vi.GetPlaneWidthSubsampling(PLANAR_U);
							const int div   = 1 << shift;

							AVSValue v,vv,vu,va;
							
							if (avsp)
							{
								AVSValue sargs[2] = {clip,"Y"};
								
								v=env->Invoke("PlaneToY",AVSValue(sargs,2)).AsClip();
								sargs[1]="U";
								vu=env->Invoke("PlaneToY",AVSValue(sargs,2)).AsClip();
								sargs[1]="V";
								vv=env->Invoke("PlaneToY",AVSValue(sargs,2)).AsClip();
								if (isAlphaChannel)
								{
									sargs[1]="A";
									va=env->Invoke("PlaneToY",AVSValue(sargs,2)).AsClip();
								}
							}
							else
							{
								vu = env->Invoke("UtoY8",clip).AsClip();
								vv = env->Invoke("VtoY8",clip).AsClip();
								v = env->Invoke("ConvertToY8",clip).AsClip();								
							}		

							v = env->Invoke(turnRightFunction,v).AsClip();
							vu = env->Invoke(turnRightFunction,vu).AsClip();
							vv = env->Invoke(turnRightFunction,vv).AsClip();
							if (isAlphaChannel) va = env->Invoke(turnRightFunction,va).AsClip();
							v = CreateResizeV(v.AsClip(), subrange_left, subrange_width, target_width,threads_number,_sleep,(step2)?1:plane_range[0],desample,accuracy,0,shift,avsp, f, env);

							VideoInfo vR = v.AsClip()->GetVideoInfo();
							int ChromaS=vR.height >> shift;

							vu = CreateResizeV(vu.AsClip(), subrange_left/div, subrange_width/div, target_width >> shift,threads_number,_sleep,(step2)?1:plane_range[1],desample,accuracy,ChromaS,0,avsp, f, env);
							vv = CreateResizeV(vv.AsClip(), subrange_left/div, subrange_width/div, target_width >> shift,threads_number,_sleep,(step2)?1:plane_range[2],desample,accuracy,ChromaS,0,avsp, f, env);
							if (isAlphaChannel) va = CreateResizeV(va.AsClip(), subrange_left, subrange_width, target_width,threads_number,_sleep,(step2)?1:plane_range[3],desample,accuracy,0,shift,avsp, f, env);
							v = env->Invoke(turnLeftFunction,v).AsClip();
							vu = env->Invoke(turnLeftFunction,vu).AsClip();
							vv = env->Invoke(turnLeftFunction,vv).AsClip();
							if (isAlphaChannel) va = env->Invoke(turnLeftFunction,va).AsClip();

						    AVSValue ytouvargs[4] = {vu,vv,v,va};
						    if (isAlphaChannel) result=env->Invoke("YtoUV",AVSValue(ytouvargs,4)).AsClip();
							else result=env->Invoke("YtoUV",AVSValue(ytouvargs,3)).AsClip();
						    if (vi.IsYUY2()) result=env->Invoke("ConvertToYUY2",result).AsClip();
						}
						else
						{
							result=env->Invoke(turnRightFunction,clip).AsClip();
							result=CreateResizeV(result, subrange_left, subrange_width, target_width,threads_number,_sleep,(step2)?1:range_mode,desample,accuracy,0,0,avsp, f, env);
							result=env->Invoke(turnLeftFunction,result).AsClip();
						}
					}
					else
					{
						result=env->Invoke(turnLeftFunction,clip).AsClip();
						result=CreateResizeV(result, subrange_left, subrange_width, target_width,threads_number,_sleep,(step2)?1:range_mode,desample,accuracy,0,0,avsp, f, env);
						result=env->Invoke(turnRightFunction,result).AsClip();
					}
				}
			}
			else
			{
				if (!vi.IsRGB() || isRGBPfamily)
				{
					if (vi.Is422() || vi.IsYUY2() || vi.IsYV411())
					{
						const int shift = vi.GetPlaneWidthSubsampling(PLANAR_U);
						const int div   = 1 << shift;

						AVSValue v,vv,vu,va;
						
						if (avsp)
						{
							AVSValue sargs[2] = {clip,"Y"};
								
							v=env->Invoke("PlaneToY",AVSValue(sargs,2)).AsClip();
							sargs[1]="U";
							vu=env->Invoke("PlaneToY",AVSValue(sargs,2)).AsClip();
							sargs[1]="V";
							vv=env->Invoke("PlaneToY",AVSValue(sargs,2)).AsClip();
							if (isAlphaChannel)
							{
								sargs[1]="A";
								va=env->Invoke("PlaneToY",AVSValue(sargs,2)).AsClip();
							}
						}
						else
						{
							vu = env->Invoke("UtoY8",clip).AsClip();
							vv = env->Invoke("VtoY8",clip).AsClip();
							v = env->Invoke("ConvertToY8",clip).AsClip();								
						}
						
						v = env->Invoke(turnRightFunction,v).AsClip();
						vu = env->Invoke(turnRightFunction,vu).AsClip();
						vv = env->Invoke(turnRightFunction,vv).AsClip();
						if (isAlphaChannel) va = env->Invoke(turnRightFunction,va).AsClip();
						v = CreateResizeV(v.AsClip(), subrange_left, subrange_width, target_width,threads_number,_sleep,(step2)?1:plane_range[0],desample,accuracy,0,shift,avsp, f, env);

						VideoInfo vR = v.AsClip()->GetVideoInfo();
						int ChromaS=vR.height >> shift;

						vu = CreateResizeV(vu.AsClip(), subrange_left/div, subrange_width/div, target_width >> shift,threads_number,_sleep,(step2)?1:plane_range[1],desample,accuracy,ChromaS,0,avsp, f, env);
						vv = CreateResizeV(vv.AsClip(), subrange_left/div, subrange_width/div, target_width >> shift,threads_number,_sleep,(step2)?1:plane_range[2],desample,accuracy,ChromaS,0,avsp, f, env);
						if (isAlphaChannel) va = CreateResizeV(va.AsClip(), subrange_left, subrange_width, target_width,threads_number,_sleep,(step2)?1:plane_range[3],desample,accuracy,0,shift,avsp, f, env);
						v = env->Invoke(turnLeftFunction,v).AsClip();
						vu = env->Invoke(turnLeftFunction,vu).AsClip();
						vv = env->Invoke(turnLeftFunction,vv).AsClip();
						if (isAlphaChannel) va = env->Invoke(turnLeftFunction,va).AsClip();

					    AVSValue ytouvargs[4] = {vu,vv,v,va};
					    if (isAlphaChannel) result=env->Invoke("YtoUV",AVSValue(ytouvargs,4)).AsClip();
						else result=env->Invoke("YtoUV",AVSValue(ytouvargs,3)).AsClip();
					    if (vi.IsYUY2()) result=env->Invoke("ConvertToYUY2",result).AsClip();
					}
					else
					{
						result=env->Invoke(turnRightFunction,clip).AsClip();
						result=CreateResizeV(result, subrange_left, subrange_width, target_width,threads_number,_sleep,(step2)?1:range_mode,desample,accuracy,0,0,avsp, f, env);
						result=env->Invoke(turnLeftFunction,result).AsClip();
					}
				}
				else
				{
					result=env->Invoke(turnLeftFunction,clip).AsClip();
					result=CreateResizeV(result, subrange_left, subrange_width, target_width,threads_number,_sleep,(step2)?1:range_mode,desample,accuracy,0,0,avsp, f, env);
					result=env->Invoke(turnRightFunction,result).AsClip();
				}
			}
		}
		if (!((subrange_top==0) && (subrange_height==target_height) && (subrange_height==vi.height)))
		{
			if (CropV)
			{
				const int mask = (vi.IsPlanar() && !grey && !isRGBPfamily) ? (1 << vi.GetPlaneHeightSubsampling(PLANAR_U)) - 1 : 0;
				
				if (((int(subrange_top) | int(subrange_height)) & mask) == 0)
				{
					if (!desample)
					{
						AVSValue sargs[6] = {result,0,int(subrange_top),vi.width,int(subrange_height),0};
						result=env->Invoke("Crop",AVSValue(sargs,6)).AsClip();
					}
				}
				else result = CreateResizeV(result, subrange_top, subrange_height, target_height,threads_number,_sleep,range_mode,desample,accuracy,0,0,avsp, f, env);
			}
			else result = CreateResizeV(result, subrange_top, subrange_height, target_height,threads_number,_sleep,range_mode,desample,accuracy,0,0,avsp, f, env);
		}
	  }
  }
  else
  {	  
	  if (VFirst)
	  {
		if ((subrange_top==0) && (subrange_height==target_height) && (subrange_height==vi.height)) result=clip;
		else
		{
			if (CropV)
			{
				const int mask = (vi.IsPlanar() && !grey && !isRGBPfamily) ? (1 << vi.GetPlaneHeightSubsampling(PLANAR_U)) - 1 : 0;

				if (((int(subrange_top) | int(subrange_height)) & mask) == 0)
				{
					if (desample) result=clip;
					else
					{
						AVSValue sargs[6] = {clip,0,int(subrange_top),vi.width,int(subrange_height),0};
						result=env->Invoke("Crop",AVSValue(sargs,6)).AsClip();
					}
				}
				else result = CreateResizeV(clip, subrange_top, subrange_height, target_height,threads_number,_sleep,(step2)?1:range_mode,desample,accuracy,0,0, avsp, f, env);
			}
			else result = CreateResizeV(clip, subrange_top, subrange_height, target_height,threads_number,_sleep,(step2)?1:range_mode,desample,accuracy,0,0,avsp, f, env);
		}
		if (!((subrange_left==0) && (subrange_width==target_width) && (subrange_width==vi.width)))
		{
			if (CropH)
			{
				const int mask = (vi.IsPlanar() && !grey && !isRGBPfamily) ? (1 << vi.GetPlaneWidthSubsampling(PLANAR_U)) - 1 : 0;

			    if (((int(subrange_left) | int(subrange_width)) & mask) == 0)
				{
					if (!desample)
					{
						AVSValue sargs[6] = {result,int(subrange_left),0,int(subrange_width),vi.height,0};
						result=env->Invoke("Crop",AVSValue(sargs,6)).AsClip();
					}
				}
				else result = CreateResizeH(result, subrange_left, subrange_width, target_width,threads_number,_sleep,range_mode,desample,accuracy,avsp, f, env);
			}
			else result = CreateResizeH(result, subrange_left, subrange_width, target_width,threads_number,_sleep,range_mode,desample,accuracy,avsp, f, env);
		}		
	  }
	  else
	  {
		if ((subrange_left==0) && (subrange_width==target_width) && (subrange_width==vi.width)) result=clip;
		else
		{
			if (CropH)
			{
				const int mask = (vi.IsPlanar() && !grey && !isRGBPfamily) ? (1 << vi.GetPlaneWidthSubsampling(PLANAR_U)) - 1 : 0;

			    if (((int(subrange_left) | int(subrange_width)) & mask) == 0)
				{
					if (desample) result=clip;
					else
					{
						AVSValue sargs[6] = {clip,int(subrange_left),0,int(subrange_width),vi.height,0};
						result=env->Invoke("Crop",AVSValue(sargs,6)).AsClip();
					}
				}
				else result = CreateResizeH(clip, subrange_left, subrange_width, target_width,threads_number,_sleep,(step2)?1:range_mode,desample,accuracy,avsp, f, env);
			}
			else result = CreateResizeH(clip, subrange_left, subrange_width, target_width,threads_number,_sleep,(step2)?1:range_mode,desample,accuracy,avsp, f, env);
		}
		if (!((subrange_top==0) && (subrange_height==target_height) && (subrange_height==vi.height)))
		{
			if (CropV)
			{
				const int mask = (vi.IsPlanar() && !grey && !isRGBPfamily) ? (1 << vi.GetPlaneHeightSubsampling(PLANAR_U)) - 1 : 0;
				
				if (((int(subrange_top) | int(subrange_height)) & mask) == 0)
				{
					if (!desample)
					{
						AVSValue sargs[6] = {result,0,int(subrange_top),vi.width,int(subrange_height),0};
						result=env->Invoke("Crop",AVSValue(sargs,6)).AsClip();
					}
				}
				else result = CreateResizeV(result, subrange_top, subrange_height, target_height,threads_number,_sleep,range_mode,desample,accuracy,0,0,avsp, f, env);
			}
			else result = CreateResizeV(result, subrange_top, subrange_height, target_height,threads_number,_sleep,range_mode,desample,accuracy,0,0,avsp, f, env);
		}
	  }
  }
  
  return result;
}


AVSValue __cdecl FilteredResizeMT::Create_PointResize(AVSValue args, void*, IScriptEnvironment* env)
{
  PointFilter f;
  return CreateResize( args[0].AsClip(), args[1].AsInt(), args[2].AsInt(),args[7].AsInt(0),
	  args[8].AsBool(true),args[9].AsBool(true),args[10].AsBool(false),args[11].AsBool(false),
	  args[12].AsInt(0),args[13].AsInt(1),false,0,0,&args[3],&f, env );
}


AVSValue __cdecl FilteredResizeMT::Create_BilinearResize(AVSValue args, void*, IScriptEnvironment* env)
{
  TriangleFilter f;
  return CreateResize( args[0].AsClip(), args[1].AsInt(), args[2].AsInt(),args[7].AsInt(0),
	  args[8].AsBool(true),args[9].AsBool(true),args[10].AsBool(false),args[11].AsBool(false),
	  args[12].AsInt(0),args[13].AsInt(1),false,0,0,&args[3],&f, env );
}


AVSValue __cdecl FilteredResizeMT::Create_BicubicResize(AVSValue args, void*, IScriptEnvironment* env)
{
  MitchellNetravaliFilter f(args[3].AsDblDef(1./3.), args[4].AsDblDef(1./3.));
  return CreateResize( args[0].AsClip(), args[1].AsInt(), args[2].AsInt(),args[9].AsInt(0),
	  args[10].AsBool(true),args[11].AsBool(true),args[12].AsBool(false),args[13].AsBool(false),
	  args[14].AsInt(0),args[15].AsInt(1),false,0,0,&args[5],&f, env );
}

AVSValue __cdecl FilteredResizeMT::Create_LanczosResize(AVSValue args, void*, IScriptEnvironment* env)
{
  LanczosFilter f(args[7].AsInt(3));
  return CreateResize( args[0].AsClip(), args[1].AsInt(), args[2].AsInt(),args[8].AsInt(0),
	  args[9].AsBool(true),args[10].AsBool(true),args[11].AsBool(false),args[12].AsBool(false),
	  args[13].AsInt(0),args[14].AsInt(1),false,0,0,&args[3],&f, env );
}

AVSValue __cdecl FilteredResizeMT::Create_Lanczos4Resize(AVSValue args, void*, IScriptEnvironment* env)
{
  LanczosFilter f(4);
  return CreateResize( args[0].AsClip(), args[1].AsInt(), args[2].AsInt(),args[7].AsInt(0),
	  args[8].AsBool(true),args[9].AsBool(true),args[10].AsBool(false),args[11].AsBool(false),
	  args[12].AsInt(0),args[13].AsInt(1),false,0,0,&args[3], &f, env );
}

AVSValue __cdecl FilteredResizeMT::Create_BlackmanResize(AVSValue args, void*, IScriptEnvironment* env)
{
  BlackmanFilter f(args[7].AsInt(4));
  return CreateResize( args[0].AsClip(), args[1].AsInt(), args[2].AsInt(),args[8].AsInt(0),
	  args[9].AsBool(true),args[10].AsBool(true),args[11].AsBool(false),args[12].AsBool(false),
	  args[13].AsInt(0),args[14].AsInt(1),false,0,0,&args[3],&f, env );
}

AVSValue __cdecl FilteredResizeMT::Create_Spline16Resize(AVSValue args, void*, IScriptEnvironment* env)
{
  Spline16Filter f;
  return CreateResize( args[0].AsClip(), args[1].AsInt(), args[2].AsInt(),args[7].AsInt(0),
	  args[8].AsBool(true),args[9].AsBool(true),args[10].AsBool(false),args[11].AsBool(false),
	  args[12].AsInt(0),args[13].AsInt(1),false,0,0,&args[3], &f, env );
}

AVSValue __cdecl FilteredResizeMT::Create_Spline36Resize(AVSValue args, void*, IScriptEnvironment* env)
{
  Spline36Filter f;
  return CreateResize( args[0].AsClip(), args[1].AsInt(), args[2].AsInt(),args[7].AsInt(0),
	  args[8].AsBool(true),args[9].AsBool(true),args[10].AsBool(false),args[11].AsBool(false),
	  args[12].AsInt(0),args[13].AsInt(1),false,0,0,&args[3], &f, env );
}

AVSValue __cdecl FilteredResizeMT::Create_Spline64Resize(AVSValue args, void*, IScriptEnvironment* env)
{
  Spline64Filter f;
  return CreateResize( args[0].AsClip(), args[1].AsInt(), args[2].AsInt(),args[7].AsInt(0),
	  args[8].AsBool(true),args[9].AsBool(true),args[10].AsBool(false),args[11].AsBool(false),
	  args[12].AsInt(0),args[13].AsInt(1),false,0,0,&args[3],&f, env );
}

AVSValue __cdecl FilteredResizeMT::Create_GaussianResize(AVSValue args, void*, IScriptEnvironment* env)
{
  GaussianFilter f(args[7].AsFloat(30.0f));
  return CreateResize( args[0].AsClip(), args[1].AsInt(), args[2].AsInt(),args[8].AsInt(0),
	  args[9].AsBool(true),args[10].AsBool(true),args[11].AsBool(false),args[12].AsBool(false),
	  args[13].AsInt(0),args[14].AsInt(1),false,0,0,&args[3],&f, env );
}

AVSValue __cdecl FilteredResizeMT::Create_SincResize(AVSValue args, void*, IScriptEnvironment* env)
{
  SincFilter f(args[7].AsInt(4));
  return CreateResize( args[0].AsClip(), args[1].AsInt(), args[2].AsInt(),args[8].AsInt(0),
	  args[9].AsBool(true),args[10].AsBool(true),args[11].AsBool(false),args[12].AsBool(false),
	  args[13].AsInt(0),args[14].AsInt(1),false,0,0,&args[3], &f, env );
}


AVSValue __cdecl FilteredResizeMT::Create_DeBilinearResize(AVSValue args, void*, IScriptEnvironment* env)
{
  TriangleFilter f;
  return CreateResize( args[0].AsClip(), args[1].AsInt(), args[2].AsInt(),args[7].AsInt(0),
	  args[8].AsBool(true),args[9].AsBool(true),args[10].AsBool(false),args[11].AsBool(false),
	  args[12].AsInt(0),args[13].AsInt(1),true,args[14].AsInt(0),args[15].AsInt(0),&args[3],&f, env );
}


AVSValue __cdecl FilteredResizeMT::Create_DeBicubicResize(AVSValue args, void*, IScriptEnvironment* env)
{
  MitchellNetravaliFilter f(args[3].AsDblDef(1./3.), args[4].AsDblDef(1./3.));
  return CreateResize( args[0].AsClip(), args[1].AsInt(), args[2].AsInt(),args[9].AsInt(0),
	  args[10].AsBool(true),args[11].AsBool(true),args[12].AsBool(false),args[13].AsBool(false),
	  args[14].AsInt(0),args[15].AsInt(1),true,args[16].AsInt(0),args[17].AsInt(0),&args[5],&f, env );
}

AVSValue __cdecl FilteredResizeMT::Create_DeLanczosResize(AVSValue args, void*, IScriptEnvironment* env)
{
  LanczosFilter f(args[7].AsInt(3));
  return CreateResize( args[0].AsClip(), args[1].AsInt(), args[2].AsInt(),args[8].AsInt(0),
	  args[9].AsBool(true),args[10].AsBool(true),args[11].AsBool(false),args[12].AsBool(false),
	  args[13].AsInt(0),args[14].AsInt(1),true,args[15].AsInt(0),args[16].AsInt(0),&args[3],&f, env );
}

AVSValue __cdecl FilteredResizeMT::Create_DeLanczos4Resize(AVSValue args, void*, IScriptEnvironment* env)
{
  LanczosFilter f(4);
  return CreateResize( args[0].AsClip(), args[1].AsInt(), args[2].AsInt(),args[7].AsInt(0),
	  args[8].AsBool(true),args[9].AsBool(true),args[10].AsBool(false),args[11].AsBool(false),
	  args[12].AsInt(0),args[13].AsInt(1),true,args[14].AsInt(0),args[15].AsInt(0),&args[3], &f, env );
}

AVSValue __cdecl FilteredResizeMT::Create_DeBlackmanResize(AVSValue args, void*, IScriptEnvironment* env)
{
  BlackmanFilter f(args[7].AsInt(4));
  return CreateResize( args[0].AsClip(), args[1].AsInt(), args[2].AsInt(),args[8].AsInt(0),
	  args[9].AsBool(true),args[10].AsBool(true),args[11].AsBool(false),args[12].AsBool(false),
	  args[13].AsInt(0),args[14].AsInt(1),true,args[15].AsInt(0),args[16].AsInt(0),&args[3],&f, env );
}

AVSValue __cdecl FilteredResizeMT::Create_DeSpline16Resize(AVSValue args, void*, IScriptEnvironment* env)
{
  Spline16Filter f;
  return CreateResize( args[0].AsClip(), args[1].AsInt(), args[2].AsInt(),args[7].AsInt(0),
	  args[8].AsBool(true),args[9].AsBool(true),args[10].AsBool(false),args[11].AsBool(false),
	  args[12].AsInt(0),args[13].AsInt(1),true,args[14].AsInt(0),args[15].AsInt(0),&args[3], &f, env );
}

AVSValue __cdecl FilteredResizeMT::Create_DeSpline36Resize(AVSValue args, void*, IScriptEnvironment* env)
{
  Spline36Filter f;
  return CreateResize( args[0].AsClip(), args[1].AsInt(), args[2].AsInt(),args[7].AsInt(0),
	  args[8].AsBool(true),args[9].AsBool(true),args[10].AsBool(false),args[11].AsBool(false),
	  args[12].AsInt(0),args[13].AsInt(1),true,args[14].AsInt(0),args[15].AsInt(0),&args[3], &f, env );
}

AVSValue __cdecl FilteredResizeMT::Create_DeSpline64Resize(AVSValue args, void*, IScriptEnvironment* env)
{
  Spline64Filter f;
  return CreateResize( args[0].AsClip(), args[1].AsInt(), args[2].AsInt(),args[7].AsInt(0),
	  args[8].AsBool(true),args[9].AsBool(true),args[10].AsBool(false),args[11].AsBool(false),
	  args[12].AsInt(0),args[13].AsInt(1),true,args[14].AsInt(0),args[15].AsInt(0),&args[3],&f, env );
}

AVSValue __cdecl FilteredResizeMT::Create_DeGaussianResize(AVSValue args, void*, IScriptEnvironment* env)
{
  GaussianFilter f(args[7].AsFloat(30.0f));
  return CreateResize( args[0].AsClip(), args[1].AsInt(), args[2].AsInt(),args[8].AsInt(0),
	  args[9].AsBool(true),args[10].AsBool(true),args[11].AsBool(false),args[12].AsBool(false),
	  args[13].AsInt(0),args[14].AsInt(1),true,args[15].AsInt(0),args[16].AsInt(0),&args[3],&f, env );
}

AVSValue __cdecl FilteredResizeMT::Create_DeSincResize(AVSValue args, void*, IScriptEnvironment* env)
{
  SincFilter f(args[7].AsInt(4));
  return CreateResize( args[0].AsClip(), args[1].AsInt(), args[2].AsInt(),args[8].AsInt(0),
	  args[9].AsBool(true),args[10].AsBool(true),args[11].AsBool(false),args[12].AsBool(false),
	  args[13].AsInt(0),args[14].AsInt(1),true,args[15].AsInt(0),args[16].AsInt(0),&args[3], &f, env );
}
